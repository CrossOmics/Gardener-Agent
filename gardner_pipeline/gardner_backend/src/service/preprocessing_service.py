import matplotlib.pyplot as plt
from typing import Optional, Dict, Any

from anndata import AnnData
from fastapi import Depends
from datetime import datetime

from dto.request.full_preprocessing_request import FullPreprocessingRequest
from dto.request.run_hvg_request import RunHVGRequest
from dto.response.full_preprocessing_response import FullPreprocessingResponse
from dto.response.hvg_result_dto import HVGResultDTO
from infrastructure.database.dao.analysis_snapshots_dao import AnalysisSnapshotsDAO
from infrastructure.filesystem.constants.filesystem_constants import QC_SUBSPACE, SNAPSHOTS_SUBSPACE, CACHE_SUBSPACE, \
    PREPROCESSING_SUBSPACE
import lotus
from lotus.preprocessing import calculate_qc_metrics, filter_cells, filter_genes, normalize_total, log1p, scale, pca, \
    neighbors
from lotus.visualization import violin, scatter, pca_variance_ratio
import lotus.visualization as visualization
from infrastructure.database.dao.project_meta_dao import ProjectMetaDAO

from infrastructure.database.dao.dataset_dao import DatasetDAO
from infrastructure.filesystem.storage import AssetStorage
from dto.response.qc_result_dto import QCResultDTO
from dto.request.filter_qc_request import FilterQCRequest
from dto.response.filter_qc_response import FilterQCResponse
from loguru import logger

from service.helper.service_constant import PREPROCESSING_BRANCH
from util.id_generate_utils import generate_business_id, generate_filename
from util.path_utils import get_dataset_relative


class PreprocessingService:
    def __init__(self, dataset_dao: DatasetDAO = Depends(),
                 snapshot_dao: AnalysisSnapshotsDAO = Depends(),
                 project_dao: ProjectMetaDAO = Depends(),
                 storage: AssetStorage = Depends()):
        self.dataset_dao = dataset_dao
        self.snapshot_dao = snapshot_dao
        self.project_dao = project_dao
        self.storage = storage

    def full_preprocessing(self, request: FullPreprocessingRequest) -> FullPreprocessingResponse:
        """
        Executes the complete preprocessing pipeline on a single AnnData object.
        """

        # Initialize Context & IDs
        logger.info(f"Starting Full Preprocessing for Dataset {request.dataset_id}...")

        snapshot_id = generate_business_id('s_pre')
        assets_base_path = f"{PREPROCESSING_SUBSPACE}/{snapshot_id}"

        # 1. Load Raw Data
        dataset = self.dataset_dao.get_dataset_by_business_id(request.dataset_id)
        if not dataset:
            raise ValueError(f"Dataset {request.dataset_id} not found.")

        try:
            adata = self.storage.load_anndata(dataset.dataset_path)
        except Exception as e:
            raise RuntimeError(f"Failed to load raw data: {e}")

        # Metrics Init
        # Capture raw cell count for later calculation
        n_cells_raw = adata.n_obs

        executed_steps = []
        skipped_steps = []
        thumbnail_map: Dict[str, Any] = {}

        start_time = datetime.utcnow()

        # Step 1: QC Calculation
        if not request.skip_qc_calculation:
            logger.info("Executing Step 1: QC Calculation")
            self._compute_qc_metrics_inplace(adata, request.organism, request.custom_prefixes)

            # Violin Plot
            # Naming Rule: Key must end with '_plot', filename matches key.
            plot_key = "qc_violin_plot"
            violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], jitter=0.4, show=False)
            fig_violin = plt.gcf()

            v_path = self.storage.save_file(
                fig_violin,
                get_dataset_relative(request.project_id, request.dataset_id, assets_base_path, f"{plot_key}.png")
            )
            thumbnail_map[plot_key] = v_path
            plt.close(fig_violin)

            # Scatter Plot (MT)
            plot_key = "qc_scatter_mt_plot"
            scatter(adata, x='total_counts', y='pct_counts_mt', show=False)
            fig_mt = plt.gcf()

            mt_path = self.storage.save_file(
                fig_mt,
                get_dataset_relative(request.project_id, request.dataset_id, assets_base_path, f"{plot_key}.png")
            )
            thumbnail_map[plot_key] = mt_path
            plt.close(fig_mt)

            executed_steps.append("QC Calculation")
        else:
            skipped_steps.append("QC Calculation")

        # Step 2: QC Filtering
        if not request.skip_qc_filter:
            logger.info("Executing Step 2: QC Filtering")

            if request.min_genes > 0: filter_cells(adata, min_genes=request.min_genes)
            if request.min_cells > 0: filter_genes(adata, min_cells=request.min_cells)

            if request.pct_mt_max is not None and 'pct_counts_mt' in adata.obs:
                adata = adata[adata.obs['pct_counts_mt'] <= request.pct_mt_max, :]

            if request.max_genes and request.max_genes > 0:
                adata = adata[adata.obs.n_genes_by_counts < request.max_genes, :]

            if request.cell_max_counts is not None:
                adata = adata[adata.obs.total_counts < request.cell_max_counts, :]

            logger.info(f"Filtered cells: {n_cells_raw} -> {adata.n_obs}")
            executed_steps.append("QC Filtering")
        else:
            skipped_steps.append("QC Filtering")

        # Step 3: HVG & Scaling
        # filter all 0 genes
        filter_genes(adata, min_cells=1)
        n_hvg_actual = 0
        if not request.skip_hvg:
            logger.info("Executing Step 3: HVG Selection")

            if 'counts' not in adata.layers: adata.layers['counts'] = adata.X.copy()
            normalize_total(adata, target_sum=request.target_sum)
            log1p(adata)

            lotus.preprocessing.highly_variable_genes(
                adata, n_top_genes=request.n_top_genes, flavor=request.flavor, subset=False
            )

            # Count actual HVGs
            n_hvg_actual = adata.var['highly_variable'].sum() if 'highly_variable' in adata.var else 0

            # HVG Plot
            plot_key = "hvg_dispersion_plot"
            lotus.visualization.highly_variable_genes(adata, show=False)
            fig_hvg = plt.gcf()

            hvg_path = self.storage.save_file(
                fig_hvg,
                get_dataset_relative(request.project_id, request.dataset_id, assets_base_path, f"{plot_key}.png")
            )
            thumbnail_map[plot_key] = hvg_path
            plt.close(fig_hvg)

            adata.raw = adata
            adata = adata[:, adata.var.highly_variable]
            scale(adata, max_value=10)

            executed_steps.append("HVG Selection")
        else:
            skipped_steps.append("HVG Selection")

        # Step 4: PCA
        if not request.skip_pca:
            logger.info("Executing Step 4: PCA")
            n_min = min(adata.n_obs, adata.n_vars)
            actual_comps = min(request.n_comps, n_min - 1)

            if actual_comps > 0:
                pca(adata, n_comps=actual_comps, svd_solver=request.svd_solver)

                # Elbow Plot
                plot_key = "pca_variance_plot"
                pca_variance_ratio(adata, log=True, show=False)
                fig_elbow = plt.gcf()

                elbow_path = self.storage.save_file(
                    fig_elbow,
                    get_dataset_relative(request.project_id, request.dataset_id, assets_base_path, f"{plot_key}.png")
                )
                thumbnail_map[plot_key] = elbow_path
                plt.close(fig_elbow)

                # PCA Scatter
                plot_key = "pca_scatter_plot"
                color_key = 'total_counts' if 'total_counts' in adata.obs else None
                visualization.pca(adata, color=color_key, show=False)
                fig_pca = plt.gcf()

                pca_path = self.storage.save_file(
                    fig_pca,
                    get_dataset_relative(request.project_id, request.dataset_id, assets_base_path, f"{plot_key}.png")
                )
                thumbnail_map[plot_key] = pca_path
                plt.close(fig_pca)

                executed_steps.append("PCA")
            else:
                skipped_steps.append("PCA (Data Small)")
        else:
            skipped_steps.append("PCA")

        # Step 5: Neighbors
        if not request.skip_neighbors:
            logger.info("Executing Step 5: Neighborhood Graph")
            if "X_pca" in adata.obsm:
                k = min(request.n_neighbors, adata.n_obs - 1)
                if k > 1:
                    neighbors(adata, n_neighbors=k, n_pcs=request.n_pcs)
                    executed_steps.append("Neighborhood Graph")
                else:
                    skipped_steps.append("Neighbors (Too few cells)")
            else:
                skipped_steps.append("Neighbors (Missing PCA)")
        else:
            skipped_steps.append("Neighborhood Graph")

        # Metrics Collection
        n_cells_final = int(adata.n_obs)
        retention = (n_cells_final / n_cells_raw * 100) if n_cells_raw > 0 else 0

        # Safe access for median metrics
        med_genes = int(adata.obs['n_genes_by_counts'].median()) if 'n_genes_by_counts' in adata.obs else 0
        med_counts = int(adata.obs['total_counts'].median()) if 'total_counts' in adata.obs else 0
        med_mt = round(float(adata.obs['pct_counts_mt'].median()), 2) if 'pct_counts_mt' in adata.obs else 0.0

        metrics_data = {
            "summary": {
                "n_cells_raw": n_cells_raw,
                "n_cells_final": n_cells_final,
                "retention_rate": f"{retention:.1f}%",
                "n_genes_final": int(adata.n_vars),
                "median_genes_per_cell": med_genes,
                "median_counts_per_cell": med_counts
            },
            "qc_details": {
                "median_pct_mt": med_mt,
                "filter_stats": {
                    "low_quality_cells_removed": n_cells_raw - n_cells_final
                }
            },
            "analysis_params": {
                "n_hvg_selected": int(n_hvg_actual),
                "pcs_used": request.n_pcs if not request.skip_pca else 0
            }
        }

        # Merge metrics into thumbnail_json
        thumbnail_map.update(metrics_data)

        # Persistence
        logger.info("Saving Final Snapshot...")

        file_name = f"{snapshot_id}.h5ad"
        relative_key = get_dataset_relative(request.project_id, request.dataset_id, SNAPSHOTS_SUBSPACE, file_name)

        try:
            saved_path = self.storage.save_anndata(adata, relative_key)
        except Exception as e:
            raise RuntimeError(f"Failed to save final snapshot: {e}")

        final_stage = executed_steps[-1] if executed_steps else "Raw Loaded"
        branch_name = PREPROCESSING_BRANCH
        current_count = AnalysisSnapshotsDAO.count_snapshots_by_branch(dataset_id=request.dataset_id,
                                                                       branch_name=branch_name)
        current_count += 1
        self.snapshot_dao.create_snapshot(
            dataset_id=request.dataset_id,
            snapshot_id=snapshot_id,
            branch_name=branch_name,
            create_time=start_time,
            snapshot_path=saved_path,
            snapshot_name=f"{branch_name} {current_count}",
            parent_snapshot_id=None,
            params_json=request.model_dump(),
            thumbnail_json=thumbnail_map,
            user_notes=f"Full pipeline execution. Steps: {', '.join(executed_steps)}"
        )

        return FullPreprocessingResponse(
            project_id=request.project_id,
            dataset_id=request.dataset_id,
            snapshot_id=snapshot_id,
            snapshot_path=saved_path,
            executed_steps=executed_steps,
            skipped_steps=skipped_steps,
            final_stage=final_stage,
            msg="Preprocessing pipeline completed successfully."
        )

    def qc_calculation(
            self,
            project_id: str,
            dataset_id: str,
            organism: str,
            custom_prefix: Optional[Dict] = None
    ) -> QCResultDTO:
        """
        Performs QC metric calculations and generates visualization assets.

        Steps:
        1. Load AnnData.
        2. Define gene prefixes based on organism (Human/Mouse).
        3. Calculate metrics via Scanpy.
        4. Save interactive JSON statistics for UI.
        5. Save static PNG plots for reports.
        """

        # Load Dataset
        dataset = self.dataset_dao.get_dataset_by_business_id(dataset_id)
        if not dataset:
            raise ValueError(f"Dataset {dataset_id} not found.")

        # Load raw data from storage
        adata = self.storage.load_anndata(dataset.dataset_path)

        cache_file_name = f"{dataset_id}_qc_metrics.h5ad"
        cache_key = get_dataset_relative(project_id, dataset_id, CACHE_SUBSPACE, cache_file_name)

        self._compute_qc_metrics_inplace(adata, organism, custom_prefix, cache_key)

        # Generate & Save UI Data (JSON)
        # Extract columns required for frontend Violin/Scatter plots
        obs_cols = ['n_genes_by_counts', 'total_counts', 'pct_counts_mt', 'pct_counts_ribo', 'pct_counts_hb']
        # Filter columns that actually exist (in case hb is missing)
        existing_cols = [col for col in obs_cols if col in adata.obs.columns]

        qc_df = adata.obs[existing_cols].copy()
        qc_df['barcode'] = qc_df.index.astype(str)

        # Convert to list of dicts for JSON serialization
        qc_json = qc_df.to_dict(orient='records')

        # Save JSON to workspace
        json_path = self.storage.save_file(
            content=qc_json,
            relative_key=get_dataset_relative(project_id, dataset_id, QC_SUBSPACE, 'qc_metrics.json'),
            project_id=project_id
        )

        # 5. Generate & Save Static Plots (PNG)
        # Set show=False to prevent display, capture fig via plt.gcf()

        # 5.1 Violin Plot
        violin(
            adata,
            ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
            jitter=0.4,
            multi_panel=True,
            show=False
        )
        fig_violin = plt.gcf()
        violin_filename = generate_filename("qc_violin_overview", 'png')
        violin_path = self.storage.save_file(
            fig_violin,
            get_dataset_relative(project_id, dataset_id, QC_SUBSPACE, violin_filename)
        )
        plt.close(fig_violin)

        # 5.2 Scatter: Counts vs MT
        scatter(
            adata,
            x='total_counts',
            y='pct_counts_mt',
            show=False
        )
        fig_mt = plt.gcf()
        qc_scatter_mt_name = generate_filename("qc_scatter_mt", "png")
        scatter_mt_path = self.storage.save_file(
            fig_mt,
            get_dataset_relative(project_id, dataset_id, QC_SUBSPACE, qc_scatter_mt_name)
        )
        plt.close(fig_mt)

        # 5.3 Scatter: Counts vs Genes
        scatter(
            adata,
            x='total_counts',
            y='n_genes_by_counts',
            show=False
        )
        fig_genes = plt.gcf()
        qc_scatter_genes_name = generate_filename("qc_scatter_genes", "png")
        scatter_genes_path = self.storage.save_file(
            fig_genes,
            get_dataset_relative(project_id, dataset_id, QC_SUBSPACE, qc_scatter_genes_name)
        )
        plt.close(fig_genes)

        # Skip saving intermediate AnnData state

        # Return DTO with all paths
        return QCResultDTO(
            dataset_id=dataset_id,
            metrics_json_path=json_path,
            violin_plot_path=violin_path,
            scatter_mt_path=scatter_mt_path,
            scatter_genes_path=scatter_genes_path
        )

    def apply_filter(self, request: FilterQCRequest) -> FilterQCResponse:
        """
        Applies QC filtering logic and persists the result.

        Strategy:
        1. Load Raw Data (Clean Slate).
        2. Attempt to merge QC metrics from the Lightweight Cache (Obs/Var only).
        3. If cache miss/invalid, fallback to re-calculating metrics.
        4. Apply filters.
        5. Save Snapshot.
        """

        # 1. Validation & Loading Raw Data
        dataset = self.dataset_dao.get_dataset_by_business_id(request.dataset_id)
        if not dataset:
            raise ValueError(f"Dataset {request.dataset_id} not found.")

        try:
            adata = self.storage.load_anndata(dataset.dataset_path)
        except FileNotFoundError:
            raise FileNotFoundError(f"Physical file missing for dataset: {dataset.dataset_path}")
        except Exception as e:
            raise RuntimeError(f"Failed to load source file: {str(e)}")

        # 2. [OPTIMIZED] Try Loading Lightweight Cache & Merge
        cache_file_name = f"{request.dataset_id}_qc_metrics.h5ad"
        cache_key = get_dataset_relative(request.project_id, request.dataset_id, CACHE_SUBSPACE, cache_file_name)
        # Use the generic merge method
        is_merged = self.storage.load_and_merge_anndata(adata, cache_key)
        # Verify if critical metrics exist (Business logic validation)
        metrics_ready = is_merged and ('total_counts' in adata.obs)

        # Fallback: Recalculate Metrics (If cache failed or is invalid)
        if not metrics_ready:
            logger.info("Fallback: Recalculating QC metrics...")
            # Retrieve Organism Info
            project = None
            try:
                project = self.project_dao.get_project_by_id(request.project_id)
            except:
                pass

            organism = project.organism if project else "Human"
            self._compute_qc_metrics_inplace(adata, organism, cache_key)

        # 4. Apply Filtering
        initial_cells = adata.n_obs
        logger.info(f"Applying filters. Initial cells: {initial_cells}")

        if request.min_genes > 0:
            filter_cells(adata, min_genes=request.min_genes)
        if request.min_cells > 0:
            filter_genes(adata, min_cells=request.min_cells)

        # Safe filtering now that metrics are guaranteed
        if request.max_counts is not None and 'total_counts' in adata.obs:
            adata = adata[adata.obs['total_counts'] <= request.max_counts, :]

        if request.pct_mt_max is not None and 'pct_counts_mt' in adata.obs:
            adata = adata[adata.obs['pct_counts_mt'] <= request.pct_mt_max, :]

        if request.pct_hb_max is not None and 'pct_counts_hb' in adata.obs:
            adata = adata[adata.obs['pct_counts_hb'] <= request.pct_hb_max, :]

        logger.info(f"Filtering complete. Remaining cells: {adata.n_obs}")

        # 5. Save Snapshot
        snapshot_id = generate_business_id('s_node_root_qc')

        file_name = f"{snapshot_id}.h5ad"
        relative_key = get_dataset_relative(request.project_id, request.dataset_id, SNAPSHOTS_SUBSPACE, file_name)

        try:
            # save the whole anndata after QC Filter
            saved_path = self.storage.save_anndata(adata, relative_key)
        except Exception as e:
            raise RuntimeError(f"Failed to save snapshot to storage: {str(e)}")

        # 6. Create DB Record
        params = request.model_dump()
        branch_name = "QC Filtered"
        current_count = AnalysisSnapshotsDAO.count_snapshots_by_branch(dataset_id=request.dataset_id,
                                                                       branch_name=branch_name)
        current_count += 1

        db_record = self.snapshot_dao.create_snapshot(
            dataset_id=dataset.dataset_id,
            snapshot_id=snapshot_id,
            branch_name=branch_name,
            snapshot_path=saved_path,
            snapshot_name=f"{branch_name} {current_count}",
            params_json=params,
            user_notes=f"Filtered from {initial_cells} to {adata.n_obs} cells."
        )

        if not db_record:
            raise RuntimeError("Database integrity error: Failed to create snapshot record.")

        return FilterQCResponse(
            snapshot_id=snapshot_id,
            snapshot_path=saved_path,
            n_obs_remaining=adata.n_obs,
            n_vars_remaining=adata.n_vars
        )

    # Internal Methods
    def _compute_qc_metrics_inplace(
            self,
            adata: AnnData,
            organism: str,
            custom_prefix: Optional[Dict] = None,
            cache_key: str = ''
    ) -> None:
        """
        Internal helper to determine prefixes based on organism and run scanpy.pp.calculate_qc_metrics.
        Modifies the AnnData object in-place.
        """
        # cache QC Calculation score
        # Record the column names before calculation
        obs_cols_before = set(adata.obs.columns)
        var_cols_before = set(adata.var.columns)

        # 1. Define Defaults
        mt_prefix = "MT-"
        ribo_prefix = ("RPS", "RPL")
        hb_pattern = "^HB[^(P)]"

        # 2. Organism Logic
        if organism and organism.lower() == "mouse":
            mt_prefix = "mt-"
            ribo_prefix = ("Rps", "Rpl")
            hb_pattern = "^Hb[^(p)]"

        # 3. Custom Overrides
        if custom_prefix:
            mt_prefix = custom_prefix.get("mt", mt_prefix)
            ribo_prefix = tuple(custom_prefix.get("ribo", ribo_prefix))
            hb_pattern = custom_prefix.get("hb", hb_pattern)

        logger.info(f"Computing QC metrics (Organism: {organism}). Prefixes: MT={mt_prefix}")

        # 4. Annotate Genes
        adata.var['mt'] = adata.var_names.str.startswith(mt_prefix)
        adata.var['ribo'] = adata.var_names.str.startswith(ribo_prefix)
        adata.var['hb'] = adata.var_names.str.contains(hb_pattern, regex=True)

        # 5. Run Scanpy Calculation
        qc_vars = ['mt', 'ribo', 'hb']
        calculate_qc_metrics(
            adata,
            qc_vars=qc_vars,
            percent_top=None,
            log1p=False,
            inplace=True
        )

        # Newly added columns = current columns - previous columns
        obs_new_cols = list(set(adata.obs.columns) - obs_cols_before)
        var_new_cols = list(set(adata.var.columns) - var_cols_before)

        #  Minimal storage Create an AnnData object containing only incremental data
        self.storage.save_incremental_anndata(
            adata_source=adata,
            obs_cols=obs_new_cols,
            var_cols=var_new_cols,
            relative_key=cache_key
        )

    def apply_hvg(self, request: RunHVGRequest) -> HVGResultDTO:
        """
        Executes Feature Selection & Scaling (Normalization -> HVG -> Scale).

        The flow follows the standard Scanpy "Golden Order":
        1. Normalize & Log1p.
        2. Identify Highly Variable Genes (HVG).
        3. Backup 'raw' state (normalized, full gene set) for DE analysis.
        4. Subset to HVGs and Scale (Z-score) for PCA/Clustering.
        """

        # Resolve Source Data Path
        # Logic: If snapshot_id is provided, use it.
        # Otherwise, find the latest "QC Filtered" snapshot from DB.

        source_path = None
        hint_message = "Please run QC calculation and QC filtering before proceeding to the Highly Variable Genes step."
        if request.snapshot_id:
            snapshot = self.snapshot_dao.get_snapshot_by_business_id(request.snapshot_id)
            if not snapshot:
                raise ValueError(
                    f"Source snapshot {request.snapshot_id} is missing. {hint_message}")
            source_path = snapshot.snapshot_path
        else:
            # Auto-find latest QC snapshot
            latest_qc = self.snapshot_dao.get_latest_snapshot(
                request.dataset_id, branch_name="QC Filtered"
            )
            if not latest_qc:
                raise ValueError(f"No QC snapshot found. Please provide snapshot_id. {hint_message}")
            source_path = latest_qc.snapshot_path
            # Update the source snapshot id
            request.snapshot_id = latest_qc.snapshot_id

        # Load Data
        try:
            adata = self.storage.load_anndata(source_path)
        except Exception as e:
            raise RuntimeError(f"Failed to load source file: {e}")

        # Normalization & Log1p
        logger.info("Step A: Normalizing and Log-transforming...")

        # Safety check: Ensure raw counts are preserved in layers if not already
        if 'counts' not in adata.layers:
            adata.layers['counts'] = adata.X.copy()

        # Normalize total counts per cell (default 1e4)
        normalize_total(adata, target_sum=request.target_sum)

        # Logarithmize the data (log(x+1))
        log1p(adata)

        # Identify Highly Variable Genes (HVG)
        logger.info(f"Step B: Identifying top {request.n_top_genes} HVGs...")

        lotus.preprocessing.highly_variable_genes(
            adata,
            n_top_genes=request.n_top_genes,
            flavor=request.flavor,
            subset=False
        )

        # Visualization (Dispersion Plot)
        # Generate the plot
        lotus.visualization.highly_variable_genes(adata, show=False)
        fig = plt.gcf()

        # Save plot to storage
        plot_filename = generate_filename("hvg_dispersion", 'png')
        hvg_plot_path = self.storage.save_file(
            fig,
            get_dataset_relative(request.project_id, request.dataset_id, QC_SUBSPACE, plot_filename)
        )
        plt.close(fig)

        # Backup Raw
        # We freeze the "Normalized + Full Gene" state into .raw
        # This is required for future Differential Expression (DE) analysis
        logger.info("Step D: Backing up full normalized data to .raw...")
        adata.raw = adata

        # Step E: Subset & Scale
        logger.info("Step E: Subsetting to HVGs and Scaling...")

        # 1. Physical subset: Keep only HVGs in X
        adata = adata[:, adata.var.highly_variable]

        # 2. Scale: Z-score normalization (make mean=0, std=1)
        # Note: This introduces negative numbers
        scale(adata, max_value=10)

        # Persistence (Snapshot)
        snapshot_id = generate_business_id('s_node_hvg')

        # Save the processed AnnData (PCA Ready)
        file_name = f"{snapshot_id}.h5ad"
        relative_key = get_dataset_relative(request.project_id, request.dataset_id, SNAPSHOTS_SUBSPACE, file_name)

        try:
            saved_path = self.storage.save_anndata(adata, relative_key)
        except Exception as e:
            raise RuntimeError(f"Failed to save snapshot: {e}")

        # Create DB Record
        n_genes = int(sum(adata.var['highly_variable'])) if 'highly_variable' in adata.var else adata.n_vars

        self.snapshot_dao.create_snapshot(
            dataset_id=request.dataset_id,
            snapshot_id=snapshot_id,
            branch_name="HVG Selected",
            snapshot_path=saved_path,
            parent_snapshot_id=request.snapshot_id,
            params_json=request.model_dump(),
            user_notes=f"Selected {n_genes} HVGs via {request.flavor}."
        )

        return HVGResultDTO(
            snapshot_id=snapshot_id,
            snapshot_path=saved_path,
            hvg_plot_path=hvg_plot_path,
            n_genes_found=adata.n_vars,  # Should be equal to n_top_genes
            msg="HVG selection and scaling complete."
        )
