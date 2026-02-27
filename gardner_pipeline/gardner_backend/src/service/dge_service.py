import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt

from fastapi import Depends
from loguru import logger
from dto.request.run_dge_request import RunDGERequest
from dto.response.dge_result_dto import DGEResultDTO

from infrastructure.database.dao.analysis_snapshots_dao import AnalysisSnapshotsDAO
from infrastructure.filesystem.storage import AssetStorage
from infrastructure.filesystem.constants.filesystem_constants import SNAPSHOTS_SUBSPACE, DGE_SUBSPACE

import lotus.annotation_analysis as aa
import lotus.visualization as visualization
from service.helper.project_id_resolver_service import ProjectIdResolverService
from service.helper.service_constant import CLUSTERING_BRANCH, DGE_BRANCH
from util.id_generate_utils import generate_filename, generate_business_id
from util.path_utils import get_dataset_relative
from util.snapshot_utils import resolve_source_snapshot


class DGEService:
    def __init__(
            self,
            snapshot_dao: AnalysisSnapshotsDAO = Depends(),
            storage: AssetStorage = Depends(),
            project_id_resolver: ProjectIdResolverService = Depends(),
    ):
        self.snapshot_dao = snapshot_dao
        self.storage = storage
        self.project_id_resolver = project_id_resolver

    def run_dge_analysis(self, request: RunDGERequest) -> DGEResultDTO:
        """
        Executes Differential Gene Expression (DGE) Analysis.
        """

        # Resolve Project ID
        request.project_id = self.project_id_resolver.resolve_project_id(
            request.project_id, request.dataset_id, request.snapshot_id
        )

        # Resolve Source Snapshot (Target: Clustering Branch)
        source_snapshot = resolve_source_snapshot(
            request.dataset_id,
            request.snapshot_id,
            target_branch=CLUSTERING_BRANCH,
            snapshot_dao=self.snapshot_dao
        )

        # Initialize Snapshot ID & Asset Paths
        snapshot_id = generate_business_id(f's_dge_{request.method}')
        assets_base_path = f"{DGE_SUBSPACE}/{snapshot_id}"

        # Load Data
        adata = self.storage.load_anndata(source_snapshot.snapshot_path)

        # Validation
        if request.groupby not in adata.obs:
            raise ValueError(
                f"Column '{request.groupby}' not found in observations. "
                f"Available columns: {list(adata.obs.columns)}"
            )

        # Clean categorical data
        if hasattr(adata.obs[request.groupby], "cat"):
            adata.obs[request.groupby] = adata.obs[request.groupby].cat.remove_unused_categories()

        n_groups = len(adata.obs[request.groupby].unique())
        if n_groups < 2:
            raise ValueError(
                f"Cannot perform DGE analysis: Grouping column '{request.groupby}' has only {n_groups} group(s). "
                "At least 2 groups are required."
            )

        # Core Calculation
        logger.info(f"Running DGE using {request.method} on group {request.groupby}...")

        aa.rank_genes_groups(
            adata,
            groupby=request.groupby,
            method=request.method,
            use_raw=request.use_raw,
            n_genes=request.n_top_genes
        )
        logger.debug("Finish DGE analysis.")
        # Metrics Collection: Top Markers
        # Extract top n_top_genes markers per group for quick preview in thumbnail/UI

        top_markers_summary = {}
        try:
            # Struct: [(gene, score), ...] per group or just list of genes
            # Simpler for JSON: Just list of gene names
            names = adata.uns['rank_genes_groups']['names']
            for group_name in names.dtype.names:
                genes = names[group_name][:request.n_top_genes].tolist()
                top_markers_summary[str(group_name)] = genes
        except Exception as e:
            logger.warning(f"Failed to extract top markers summary: {e}")

        # Export CSV (Full Results)
        csv_path = ""
        try:
            # None implies extracting all groups
            dge_df = aa.rank_genes_groups_df(adata, group=None)
            csv_filename = generate_filename('markers_all', 'csv')
            csv_path = self.storage.save_df_to_csv(
                dge_df,
                get_dataset_relative(request.project_id, request.dataset_id, assets_base_path, csv_filename)
            )
        except Exception as e:
            logger.error(f"Failed to export CSV: {e}")

        # Visualization
        thumbnail_map = {
            "n_groups": n_groups,
            "n_top_genes_calculated": request.n_top_genes,
            "top_markers": top_markers_summary,
            "results_file": csv_path
        }
        logger.debug("start saving the images")

        # Helper to save plot with strictly defined name
        def save_plot(fig, key_name):
            if not fig or not fig.get_axes():
                return None
            filename = f"{key_name}.png"
            path = self.storage.save_file(
                fig,
                get_dataset_relative(request.project_id, request.dataset_id, assets_base_path, filename),
                dpi=300
            )
            plt.close(fig)
            return path

        try:
            logger.debug("start ploting rank gene")
            # Rank Genes Plot
            plt.close('all')
            visualization.rank_genes_groups(adata, n_genes=20, sharey=False, show=False)
            thumbnail_map["dge_rank_plot"] = save_plot(plt.gcf(), "dge_rank_plot")

            logger.debug("start ploting Dot Plot")
            # Dot Plot
            plt.close('all')
            visualization.rank_genes_groups_dotplot(adata, n_genes=5, show=False)
            thumbnail_map["dge_dot_plot"] = save_plot(plt.gcf(), "dge_dot_plot")

            # Violin Plot
            logger.debug("start ploting Violin")
            plt.close('all')

            visualization.rank_genes_groups_violin(adata, n_genes=5, show=False)
            fig = plt.gcf()
            if fig.get_axes():
                thumbnail_map["dge_violin_plot"] = save_plot(fig, "dge_violin_plot")
            else:
                logger.warning("Violin plot generated an empty figure.")

            # 4. Heatmap
            # Key: dge_heatmap_plot
            # Only generate heatmap if not too many genes/groups to avoid clutter
            # plt.close('all')
            # visualization.rank_genes_groups_heatmap(adata, n_genes=5, show_gene_labels=True, show=False)
            # thumbnail_map["dge_heatmap_plot"] = save_plot(plt.gcf(), "dge_heatmap_plot")

        except Exception as e:
            logger.error(f"Error generating DGE plots: {e}")

        # Persist Snapshot
        file_name = f"{snapshot_id}.h5ad"
        relative_key = get_dataset_relative(request.project_id, request.dataset_id, SNAPSHOTS_SUBSPACE, file_name)

        logger.debug("start saving anndata")
        saved_path = self.storage.save_anndata(adata, relative_key)

        branch_name = DGE_BRANCH
        current_count = AnalysisSnapshotsDAO.count_snapshots_by_branch(dataset_id=request.dataset_id,
                                                                       branch_name=branch_name)
        current_count += 1

        # Update Database
        self.snapshot_dao.create_snapshot(
            dataset_id=request.dataset_id,
            snapshot_id=snapshot_id,
            branch_name=DGE_BRANCH,
            snapshot_path=saved_path,
            snapshot_name=f"{branch_name} {current_count}",
            parent_snapshot_id=source_snapshot.snapshot_id,
            params_json=request.model_dump(),
            thumbnail_json=thumbnail_map,
            user_notes=f"DGE ({request.method}) on '{request.groupby}'. Identified markers for {n_groups} groups."
        )
        logger.debug("Finish DGE Process")
        return DGEResultDTO(
            snapshot_id=snapshot_id,
            snapshot_path=saved_path,
            csv_path=csv_path,
            rank_genes_plot=thumbnail_map.get("dge_rank_plot", ""),
            dotplot=thumbnail_map.get("dge_dot_plot", ""),
            heatmap=thumbnail_map.get("dge_heatmap_plot", ""),
            violin=thumbnail_map.get("dge_violin_plot", ""),
            top_markers=top_markers_summary,
            msg="Differential gene expression analysis complete."
        )
