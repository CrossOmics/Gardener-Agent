import math
from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd
import scanpy as sc
from fastapi import Depends
from loguru import logger
from typing import Dict, Any, Optional
import matplotlib

import matplotlib.colors as mcolors
from lotus.annotation_analysis.gseapy_core import run_enrichr_analysis
from lotus.annotation_analysis.celltypist_core import run_celltypist_annotation

from infrastructure.database.dao.analysis_snapshots_dao import AnalysisSnapshotsDAO
from infrastructure.filesystem.storage import AssetStorage
from infrastructure.filesystem.constants.filesystem_constants import SNAPSHOTS_SUBSPACE, ANNOTATION_SUBSPACE
from util.id_generate_utils import generate_business_id, generate_filename
from util.path_utils import get_dataset_relative
from util.snapshot_utils import resolve_source_snapshot

from service.helper.project_id_resolver_service import ProjectIdResolverService
from service.helper.service_constant import DGE_BRANCH, CLUSTERING_BRANCH, ANNOTATION_BRANCH

from dto.request.annotation_request import RunAnnotationRequest, RunCellTypistRequest, RunFullAnnotationRequest, \
    UpdateAnnotationLabelRequest
from dto.response.annotation_result_dto import AnnotationResultDTO, UpdateAnnotationLabelResponse

matplotlib.use('Agg')


class AnnotationService:
    def __init__(
            self,
            snapshot_dao: AnalysisSnapshotsDAO = Depends(),
            storage: AssetStorage = Depends(),
            project_id_resolver: ProjectIdResolverService = Depends(),
    ):
        self.snapshot_dao = snapshot_dao
        self.storage = storage
        self.project_id_resolver = project_id_resolver

    def calculate_gseapy_confidence(self, adj_p_value: float) -> float:
        """
        Normalize GSEApy Adjusted P-value to a 0-1 confidence score.
        Uses Tanh(-log10(p)) to handle the wide dynamic range of P-values.
        """
        # 1. Safety: Prevent log(0) error. Enrichr sometimes returns 0 for p-value.
        # Use a small epsilon (1e-300) which is near the float limit.
        safe_p = max(float(adj_p_value), 1e-300)

        # 2. Calculate -log10 score
        # e.g., p=0.01 -> score=2; p=1e-10 -> score=10
        nlog_p = -np.log10(safe_p)

        # 3. Tanh Normalization
        # Scale=10 means p=1e-10 gives ~76% confidence.
        # p=1e-20 gives ~96% confidence.
        confidence = np.tanh(nlog_p / 10)

        return round(confidence, 4)

    def run_full_annotation(self, request: RunFullAnnotationRequest) -> AnnotationResultDTO:
        """
        [Refactored] Consolidated annotation pipeline with consistent color binding.
        It binds the new labels (T-cells) to the original Cluster ID colors (Cluster 0/Blue) immediately.
        """

        # 1. Validation & Setup
        run_gsea = request.categories is not None and len(request.categories) > 0
        run_ct = request.model_names is not None and len(request.model_names) > 0

        if not run_gsea and not run_ct:
            raise ValueError("No annotation tasks selected. Please provide 'categories' or 'model_names'.")

        # Resolve IDs
        request.project_id = self.project_id_resolver.resolve_project_id(
            request.project_id, request.dataset_id, request.snapshot_id
        )
        # Choose branch based on task
        target_branch = DGE_BRANCH if run_gsea else CLUSTERING_BRANCH
        source_snapshot = resolve_source_snapshot(
            request.dataset_id, request.snapshot_id,
            target_branch=target_branch, snapshot_dao=self.snapshot_dao
        )

        snapshot_id = generate_business_id('s_ann')
        assets_base_path = f"{ANNOTATION_SUBSPACE}/{snapshot_id}"

        # Load Data
        adata = self.storage.load_anndata(source_snapshot.snapshot_path)

        # Containers
        full_report_data = {}
        thumbnail_data = {}
        msg_parts = []

        # STEP 1: Identify Source of Truth for Clustering (ID & Colors)
        cluster_col = None

        # Try to find the clustering column used for DEG (Best Source)
        if 'rank_genes_groups' in adata.uns:
            cluster_col = adata.uns['rank_genes_groups']['params']['groupby']

        # Fallback guessing if not found in UNS
        if not cluster_col:
            if "leiden" in adata.obs:
                cluster_col = "leiden"
            elif "louvain" in adata.obs:
                cluster_col = "louvain"

        # Ensure Source Colors Exist (The "Binding" Preparation)
        if cluster_col and cluster_col in adata.obs:
            # Force categorical
            if not hasattr(adata.obs[cluster_col], 'cat'):
                adata.obs[cluster_col] = adata.obs[cluster_col].astype('category')

            # Generate colors if they don't exist yet (Critical for binding later)
            if f"{cluster_col}_colors" not in adata.uns:
                logger.info(f"Initializing fixed color palette for {cluster_col}")
                sc.pl.umap(adata, color=cluster_col, show=False)

        # Helper to extract colors for JSON response
        def get_category_colors(adata_obj, column_name):
            try:
                if f"{column_name}_colors" in adata_obj.uns and hasattr(adata_obj.obs[column_name], 'cat'):
                    categories = adata_obj.obs[column_name].cat.categories
                    colors = adata_obj.uns[f"{column_name}_colors"]
                    return dict(zip(categories, colors))
            except Exception:
                pass
            return {}

        # 2. GSEApy Step
        if run_gsea:
            if not cluster_col:
                raise ValueError("GSEApy requested but clustering data not found.")

            available_clusters = adata.obs[cluster_col].unique()  # Use obs to safely get list

            category_label_maps = {cat: {} for cat in request.categories}
            temp_gsea_summaries = {f"gseapy_{cat.replace(' ', '_').replace('-', '_')}": {} for cat in
                                   request.categories}

            # Initialize report keys
            for cat in request.categories:
                safe_cat = cat.replace(" ", "_").replace("-", "_")
                report_key = f"gseapy_{safe_cat}"
                full_report_data[report_key] = {}
                thumbnail_data[report_key] = {}

            # Iterate Clusters
            for cluster_id in available_clusters:
                cluster_str = str(cluster_id)

                # Extract markers
                df_markers = sc.get.rank_genes_groups_df(adata, group=cluster_id)
                gene_list = df_markers.head(request.top_n_genes)['names'].tolist()

                if not gene_list:
                    # Handle empty gene list
                    for cat in request.categories:
                        safe_cat = cat.replace(" ", "_").replace("-", "_")
                        report_key = f"gseapy_{safe_cat}"
                        category_label_maps[cat][cluster_str] = "Unknown"
                        temp_gsea_summaries[report_key][cluster_str] = {
                            "predicted_cell_type": "Unknown", "average_confidence": 0.0
                        }
                    continue

                # Run Enrichr
                enrich_result = run_enrichr_analysis(gene_list, request.categories)

                # Process results
                for cat in request.categories:
                    safe_cat = cat.replace(" ", "_").replace("-", "_")
                    report_key = f"gseapy_{safe_cat}"

                    if cat in enrich_result and enrich_result[cat].get("top_terms"):
                        top = enrich_result[cat]["top_terms"][0]
                        term_name = top.get("Term")
                        adj_p = top.get("Adjusted P-value", 1.0)
                        confidence = self.calculate_gseapy_confidence(adj_p)

                        # Save full report
                        full_report_data[report_key][cluster_str] = {
                            "predicted_cell_type": term_name,
                            "p_value": adj_p,
                            "combined_score": top.get("Combined Score", 0),
                            "average_confidence": confidence
                        }
                        # Save temp summary
                        temp_gsea_summaries[report_key][cluster_str] = {
                            "predicted_cell_type": term_name,
                            "average_confidence": confidence
                        }
                        category_label_maps[cat][cluster_str] = term_name
                    else:
                        full_report_data[report_key][cluster_str] = {"predicted_cell_type": "Unannotated"}
                        temp_gsea_summaries[report_key][cluster_str] = {
                            "predicted_cell_type": "Unannotated", "average_confidence": 0.0
                        }
                        category_label_maps[cat][cluster_str] = "Unannotated"

            # Apply Labels & BIND COLORS
            for cat in request.categories:
                safe_cat = cat.replace(" ", "_").replace("-", "_")
                report_key = f"gseapy_{safe_cat}"
                col_name = f"obs_{report_key}"

                # Map labels to new column
                adata.obs[col_name] = adata.obs[cluster_col].astype(str).map(category_label_maps[cat])
                adata.obs[col_name] = adata.obs[col_name].astype("category")

                # === BIND COLORS (Crucial Step) ===
                # This ensures "T-cells" gets the same color as "Cluster 0" (Blue)
                self._apply_cluster_colors_to_annotation(adata, cluster_col, col_name)

                # Plot UMAP (Now uses the synced colors)
                plot_key = report_key + "_plot"
                p_path = self._plot_annotated_umap(
                    request.project_id, request.dataset_id, assets_base_path,
                    adata, color_col=col_name, title=f"GSEApy: {cat}", suffix=report_key
                )
                if p_path:
                    thumbnail_data[plot_key] = p_path

                # Update Summary JSON with Colors
                label_color_map = get_category_colors(adata, col_name)
                current_summary = temp_gsea_summaries[report_key]
                for c_id, info in current_summary.items():
                    label = info.get("predicted_cell_type")
                    info["color"] = label_color_map.get(label, "#808080")

                thumbnail_data[report_key] = current_summary

            msg_parts.append(f"GSEApy({len(request.categories)})")

        # 3. CellTypist Step
        if run_ct:
            logger.info("Running CellTypist Annotation...")

            for model_name in request.model_names:
                try:
                    safe_model = model_name.replace('.pkl', '').replace(' ', '_').replace('-', '_').lower()
                    report_key = f"celltypist_{safe_model}"

                    res = run_celltypist_annotation(
                        adata=adata, model_name=model_name,
                        majority_voting=request.majority_voting,
                        target_cluster_col=request.target_cluster_col
                    )

                    if "error" not in res:
                        cluster_summary = res.get("cluster_summary", {})
                        full_report_data[report_key] = cluster_summary

                        # Determine label column
                        label_col = "majority_voting" if request.majority_voting else "predicted_labels"
                        col_name = f"obs_{report_key}"

                        # Save to Obs & Force Categorical
                        adata.obs[col_name] = adata.obs[label_col].astype("category")

                        # === BIND COLORS ===
                        # Sync colors from original clusters to CellTypist labels for consistency
                        if cluster_col:
                            self._apply_cluster_colors_to_annotation(adata, cluster_col, col_name)

                        # Plot
                        plot_key = report_key + "_plot"
                        p_path = self._plot_annotated_umap(
                            request.project_id, request.dataset_id, assets_base_path,
                            adata, color_col=col_name, title=f"CellTypist: {model_name}", suffix=report_key
                        )
                        if p_path:
                            thumbnail_data[plot_key] = p_path

                        # Extract Colors & Build Summary
                        label_color_map = get_category_colors(adata, col_name)
                        simple_summary = {}
                        for c_id, details in cluster_summary.items():
                            pred_label = details.get("predicted_cell_type", "Unknown")
                            simple_summary[c_id] = {
                                "predicted_cell_type": pred_label,
                                "average_confidence": details.get("average_confidence", 0.0),
                                "color": label_color_map.get(pred_label, "#808080")
                            }
                        thumbnail_data[report_key] = simple_summary
                    else:
                        full_report_data[report_key] = {"error": res["error"]}

                except Exception as e:
                    logger.error(f"CellTypist model {model_name} failed: {e}")
                    full_report_data[f"celltypist_{model_name}_error"] = {"error": str(e)}

            msg_parts.append(f"CellTypist({len(request.model_names)})")

        # 4. Persistence
        json_filename = "json_report.json"
        json_path = get_dataset_relative(request.project_id, request.dataset_id, assets_base_path, json_filename)
        self.storage.save_file(full_report_data, json_path)

        saved_path = self.storage.save_anndata(
            adata,
            get_dataset_relative(request.project_id, request.dataset_id, SNAPSHOTS_SUBSPACE, f"{snapshot_id}.h5ad")
        )

        thumbnail_data["report_json"] = json_path

        # IMPORTANT: Save 'target_cluster_col' in params
        # This tells 'update_cluster_labels' which column is the "Source of Truth" for colors.
        params_to_save = request.model_dump()
        if cluster_col:
            params_to_save['target_cluster_col'] = cluster_col

        branch_name = ANNOTATION_BRANCH
        current_count = AnalysisSnapshotsDAO.count_snapshots_by_branch(dataset_id=request.dataset_id,
                                                                       branch_name=branch_name)
        current_count += 1

        self.snapshot_dao.create_snapshot(
            dataset_id=request.dataset_id,
            snapshot_id=snapshot_id,
            branch_name=ANNOTATION_BRANCH,
            snapshot_path=saved_path,
            snapshot_name=f"{branch_name} {current_count}",
            parent_snapshot_id=source_snapshot.snapshot_id,
            thumbnail_json=thumbnail_data,
            params_json=params_to_save,
            user_notes=f"Combined Annotation: {', '.join(msg_parts)}"
        )

        api_results = {
            "report_path": json_path,
            "result_keys": list(full_report_data.keys()),
            "plots": [v for k, v in thumbnail_data.items() if k.endswith("_plot") and k != "report_json"]
        }

        return AnnotationResultDTO(
            snapshot_id=snapshot_id,
            snapshot_path=saved_path,
            cluster_id="ALL",
            enrichment_results=api_results,
            msg=f"Full annotation completed: {', '.join(msg_parts)}"
        )

    def update_cluster_labels(self, request: UpdateAnnotationLabelRequest) -> UpdateAnnotationLabelResponse:
        """
        Updates cluster labels in metadata AND regenerates the UMAP plot and .h5ad file.
        Ensures colors remain consistent with original clusters even after renaming.
        """
        snapshot_id = request.snapshot_id
        target_key = request.file_name  # e.g., "gseapy_CellMarker_2024"

        # 1. Retrieve Snapshot
        snapshot = self.snapshot_dao.get_snapshot_by_business_id(snapshot_id)
        if not snapshot:
            raise ValueError(f"Snapshot '{snapshot_id}' not found.")

        # 2. Update Database Metadata (JSON)
        current_data = snapshot.thumbnail_json or {}
        if target_key not in current_data:
            raise ValueError(f"Annotation category '{target_key}' not found in snapshot metadata.")

        annotation_group = current_data[target_key]
        changes_count = 0

        # Track updates for AnnData modification later
        updates_map = {}

        for cluster_id, new_label in request.updated_annotation.items():
            cid_str = str(cluster_id)
            if cid_str in annotation_group:
                old_label = annotation_group[cid_str].get("predicted_cell_type", "Unknown")

                # Update DB JSON
                annotation_group[cid_str]["predicted_cell_type"] = new_label

                # Record for AnnData update
                updates_map[cid_str] = new_label
                changes_count += 1

                logger.info(f"Updated Cluster {cid_str}: '{old_label}' -> '{new_label}'")

        # If nothing changed, return early
        if changes_count == 0:
            return UpdateAnnotationLabelResponse(
                msg="No matching clusters found to update.",
                snapshot_id=snapshot_id,
                updated_count=0,
                target_file=target_key
            )

        logger.info("Regenerating UMAP and updating AnnData file...")

        # 3. Load AnnData
        adata = self.storage.load_anndata(snapshot.snapshot_path)

        # 4. Identify Columns
        # Annotation col is "obs_{target_key}"
        annot_col = f"obs_{target_key}"

        # Try to find the original clustering column to map updates
        params = snapshot.params_json or {}
        cluster_col = params.get("target_cluster_col")

        if not cluster_col:
            # Fallback guessing
            if "leiden" in adata.obs:
                cluster_col = "leiden"
            elif "louvain" in adata.obs:
                cluster_col = "louvain"

        if not cluster_col or cluster_col not in adata.obs:
            logger.warning("Could not identify clustering column. UMAP update skipped, only DB updated.")
        else:
            # 5. Modify AnnData.obs
            # Convert clustering col to string for matching
            clusters_series = adata.obs[cluster_col].astype(str)

            # Get current annotations as series
            if annot_col in adata.obs:
                current_annots = adata.obs[annot_col].astype(object)
            else:
                current_annots = clusters_series.copy().astype(object)

            # Apply updates: For every cell in cluster X, update annotation to New Label
            for cid, new_lbl in updates_map.items():
                mask = (clusters_series == cid)
                current_annots[mask] = new_lbl

            # Assign back and force categorical (This usually resets colors)
            adata.obs[annot_col] = current_annots.astype("category")

            cluster_colors_key = f"{cluster_col}_colors"

            # map original cluster colors to the new labels
            if cluster_colors_key in adata.uns and hasattr(adata.obs[cluster_col], 'cat'):
                try:
                    # Map: { '0': '#1f77b4', '1': '#ff7f0e' }
                    src_cats = adata.obs[cluster_col].cat.categories
                    src_colors = adata.uns[cluster_colors_key]

                    if len(src_cats) == len(src_colors):
                        id_to_color = dict(zip(src_cats, src_colors))

                        # Build new color list matching the order of NEW categories
                        new_cats = adata.obs[annot_col].cat.categories
                        new_colors_list = []

                        for label in new_cats:
                            # Find a representative cell for this label
                            mask = (adata.obs[annot_col] == label)
                            if mask.any():
                                # Get the original numeric Cluster ID (e.g., '0')
                                rep_cluster_id = adata.obs[cluster_col][mask].iloc[0]
                                # Apply the original color (e.g., Blue)
                                color = id_to_color.get(str(rep_cluster_id), "#808080")
                                new_colors_list.append(color)
                            else:
                                new_colors_list.append("#808080")

                        # Assign manually constructed colors to .uns
                        adata.uns[f"{annot_col}_colors"] = new_colors_list

                except Exception as e:
                    logger.warning(f"Could not preserve colors: {e}")

            # 6. Regenerate UMAP Plot using helper
            plot_key = f"{target_key}_plot"
            assets_base_path = f"{ANNOTATION_SUBSPACE}/{snapshot_id}"

            project_id = self.project_id_resolver.resolve_project_id(snapshot_id=request.snapshot_id)
            dataset_id = self.project_id_resolver.resolve_dataset_id(snapshot_id=request.snapshot_id)

            # Generate Updated Plot
            new_plot_path = self._plot_annotated_umap(
                project_id, dataset_id, assets_base_path,
                adata, color_col=annot_col,
                title=f"{target_key} (Updated)",
                suffix=f"{target_key}"
            )

            # 7. Update image path in metadata
            if new_plot_path:
                current_data[plot_key] = new_plot_path

            # 8. Extract New Colors for JSON metadata
            def local_get_colors(adata_obj, col):
                try:
                    if f"{col}_colors" in adata_obj.uns and hasattr(adata_obj.obs[col], 'cat'):
                        cats = adata_obj.obs[col].cat.categories
                        colors = adata_obj.uns[f"{col}_colors"]
                        return dict(zip(cats, [mcolors.to_hex(c) for c in colors]))
                except Exception:
                    return {}
                return {}

            new_colors = local_get_colors(adata, annot_col)
            # Update 'color' field in JSON
            for cid, info in annotation_group.items():
                lbl = info.get("predicted_cell_type")
                if lbl in new_colors:
                    info["color"] = new_colors[lbl]

            # 9. Save AnnData (Overwrite existing .h5ad)
            self.storage.save_anndata(adata, snapshot.snapshot_path)

        # 10. Save DB changes
        snapshot.thumbnail_json = current_data
        snapshot.save()

        return UpdateAnnotationLabelResponse(
            msg="Annotation labels and UMAP updated successfully.",
            snapshot_id=snapshot_id,
            updated_count=changes_count,
            target_file=target_key
        )

    def run_single_gseapy_annotation(self, request: RunAnnotationRequest) -> AnnotationResultDTO:
        """
        Executes GSEApy annotation for a single cluster.
        """
        # Resolve Project ID
        request.project_id = self.project_id_resolver.resolve_project_id(
            request.project_id, request.dataset_id, request.snapshot_id
        )

        # Resolve Source Snapshot (Target: DGE Branch)
        source_snapshot = resolve_source_snapshot(
            request.dataset_id,
            request.snapshot_id,
            target_branch=DGE_BRANCH,
            snapshot_dao=self.snapshot_dao
        )

        # Initialize Snapshot ID & Asset Paths
        snapshot_id = generate_business_id(f's_gsea_c{request.cluster_id}')
        assets_base_path = f"{ANNOTATION_SUBSPACE}/{snapshot_id}"

        # Load Data
        adata = self.storage.load_anndata(source_snapshot.snapshot_path)

        if "rank_genes_groups" not in adata.uns:
            raise ValueError("Marker genes (rank_genes_groups) not found in the source snapshot.")

        try:
            df_markers = sc.get.rank_genes_groups_df(adata, group=request.cluster_id)
            gene_list = df_markers.head(request.top_n_genes)['names'].tolist()
        except KeyError:
            raise ValueError(f"Cluster ID '{request.cluster_id}' not found in marker data.")

        # 5. Core Analysis
        logger.info(f"Running Enrichr on {len(gene_list)} genes for Cluster {request.cluster_id}...")
        enrichment_results = run_enrichr_analysis(
            gene_list=gene_list,
            categories=request.categories
        )

        # Save JSON Report
        json_filename = f'{snapshot_id}.json'
        json_relative_path = get_dataset_relative(
            request.project_id, request.dataset_id, assets_base_path, json_filename
        )
        try:
            self.storage.save_file(enrichment_results, json_relative_path)
        except Exception as e:
            logger.error(f"Failed to save annotation JSON: {e}")
            raise e

        # Generate Plot
        plot_path = self._plot_enrichment_bar(
            request.project_id,
            request.dataset_id,
            assets_base_path,
            request.cluster_id,
            enrichment_results
        )

        # Persist Snapshot (H5AD)
        file_name = f"{snapshot_id}.h5ad"
        h5ad_relative_key = get_dataset_relative(
            request.project_id, request.dataset_id, SNAPSHOTS_SUBSPACE, file_name
        )
        saved_path = self.storage.save_anndata(adata, h5ad_relative_key)

        branch_name = ANNOTATION_BRANCH
        current_count = AnalysisSnapshotsDAO.count_snapshots_by_branch(dataset_id=request.dataset_id,
                                                                       branch_name=branch_name)
        current_count += 1
        # Create DB Record
        self.snapshot_dao.create_snapshot(
            dataset_id=request.dataset_id,
            snapshot_id=snapshot_id,
            branch_name=ANNOTATION_BRANCH,
            snapshot_path=saved_path,
            snapshot_name=f"{branch_name} {current_count}",
            parent_snapshot_id=source_snapshot.snapshot_id,
            thumbnail_json={
                "annotation_report_path": json_relative_path,
                "enrichment_plot_path": plot_path
            },
            params_json=request.model_dump(),
            user_notes=f"Annotated Cluster {request.cluster_id} using Enrichr."
        )

        return AnnotationResultDTO(
            snapshot_id=snapshot_id,
            snapshot_path=saved_path,
            cluster_id=request.cluster_id,
            enrichment_results=enrichment_results,
            msg="Annotation complete. Plot generated."
        )

    def run_all_gseapy_annotation(self, request: RunAnnotationRequest) -> AnnotationResultDTO:
        """
        [Optimized] Auto-annotates all clusters across MULTIPLE categories (databases).
        Saves a summary JSON where Key=ClusterID and Value=Metrics (Label, P-value, Log10 Score).
        """
        # 1. Resolve IDs
        request.project_id = self.project_id_resolver.resolve_project_id(
            request.project_id, request.dataset_id, request.snapshot_id
        )

        source_snapshot = resolve_source_snapshot(
            request.dataset_id, request.snapshot_id,
            target_branch=DGE_BRANCH, snapshot_dao=self.snapshot_dao
        )

        snapshot_id = generate_business_id('s_gsea')
        assets_base_path = f"{ANNOTATION_SUBSPACE}/{snapshot_id}"
        adata = self.storage.load_anndata(source_snapshot.snapshot_path)

        # 2. Determine Clusters
        try:
            groupby_key = adata.uns['rank_genes_groups']['params']['groupby']
            available_clusters = adata.uns["rank_genes_groups"]["names"].dtype.names
        except (KeyError, Exception):
            groupby_key = 'leiden' if 'leiden' in adata.obs else 'louvain'
            available_clusters = adata.obs[groupby_key].unique()

        # 3. Core Logic: Iterate Clusters -> Enrichr -> Calculate Scores
        summary_map: Dict[str, Any] = {}
        cluster_to_celltype: Dict[str, str] = {}
        target_categories = request.categories if request.categories else ["cellmarker"]

        for cluster_id in available_clusters:
            try:
                df_markers = sc.get.rank_genes_groups_df(adata, group=cluster_id)
                gene_list = df_markers.head(request.top_n_genes)['names'].tolist()

                if not gene_list:
                    cluster_to_celltype[str(cluster_id)] = "Unknown"
                    continue

                # Run Analysis
                enrich_result = run_enrichr_analysis(gene_list, target_categories)

                # Find best hit across databases
                best_overall_hit = {"cell_type": "Unknown", "p_value": 1.0, "source": "None", "neg_log10_p": 0.0,
                                    "confidence": 0.0}
                cluster_details = {}

                for cat in target_categories:
                    if cat in enrich_result and enrich_result[cat].get("top_terms"):
                        top_hit = enrich_result[cat]["top_terms"][0]
                        adj_p = top_hit.get("Adjusted P-value", 1.0)

                        # --- Logic Change: Calculate Scores ---
                        neg_log10 = -math.log10(adj_p) if adj_p > 0 else 300.0
                        confidence = round(max(0, 1.0 - adj_p), 4)

                        # Record details
                        cluster_details[cat] = {
                            "term": top_hit.get("Term"),
                            "adj_p": adj_p,
                            "neg_log10_p": round(neg_log10, 4),
                            "confidence": confidence
                        }

                        # Update Best Hit
                        if adj_p < best_overall_hit["p_value"]:
                            best_overall_hit = {
                                "cell_type": top_hit.get("Term"),
                                "p_value": adj_p,
                                "neg_log10_p": round(neg_log10, 4),
                                "confidence": confidence,
                                "source": cat
                            }

                # Store in Summary Map (This is what gets saved to JSON)
                summary_map[str(cluster_id)] = {
                    "best_match": best_overall_hit,
                    "all_categories": cluster_details
                }
                cluster_to_celltype[str(cluster_id)] = best_overall_hit["cell_type"]

            except Exception as e:
                logger.error(f"Failed annotating cluster {cluster_id}: {e}")
                cluster_to_celltype[str(cluster_id)] = "Error"

        # 4. Map labels to AnnData
        new_col_name = "gsea_auto_labels"
        adata.obs[new_col_name] = adata.obs[groupby_key].astype(str).map(cluster_to_celltype)

        # 5. Save JSON Asset (Key = Cluster ID, Value = { label, log10, conf ... })
        json_path = get_dataset_relative(request.project_id, request.dataset_id, assets_base_path,
                                         f"{snapshot_id}.json")
        self.storage.save_file(summary_map, json_path)

        # 6. Generate Plot
        plot_path = self._plot_annotated_umap(
            request.project_id, request.dataset_id, assets_base_path,
            adata, color_col=new_col_name, title="GSEApy Multi-Database Annotation", suffix="_summary"
        )

        # 7. Persist Snapshot
        saved_path = self.storage.save_anndata(adata, get_dataset_relative(request.project_id, request.dataset_id,
                                                                           SNAPSHOTS_SUBSPACE, f"{snapshot_id}.h5ad"))

        branch_name = ANNOTATION_BRANCH
        current_count = AnalysisSnapshotsDAO.count_snapshots_by_branch(dataset_id=request.dataset_id,
                                                                       branch_name=branch_name)
        current_count += 1
        self.snapshot_dao.create_snapshot(
            dataset_id=request.dataset_id, snapshot_id=snapshot_id, branch_name=ANNOTATION_BRANCH,
            snapshot_path=saved_path, parent_snapshot_id=source_snapshot.snapshot_id,
            thumbnail_json={"report": json_path, "umap": plot_path},
            snapshot_name=f"{branch_name} {current_count}",
            params_json=request.model_dump(), user_notes="Multi-database GSEApy annotation."
        )

        return AnnotationResultDTO(
            snapshot_id=snapshot_id, snapshot_path=saved_path,
            cluster_id="ALL", enrichment_results=summary_map, msg="Batch annotation complete."
        )

    def run_celltypist_annotations(self, request: RunCellTypistRequest) -> AnnotationResultDTO:
        """
        Executes CellTypist for MULTIPLE models.
        Generates a labeled UMAP for EACH model selected.
        """
        # Resolve IDs and load data
        request.project_id = self.project_id_resolver.resolve_project_id(
            request.project_id, request.dataset_id, request.snapshot_id
        )
        source_snapshot = resolve_source_snapshot(
            request.dataset_id, request.snapshot_id,
            target_branch=CLUSTERING_BRANCH, snapshot_dao=self.snapshot_dao
        )

        snapshot_id = generate_business_id('s_ct')
        assets_base_path = f"{ANNOTATION_SUBSPACE}/{snapshot_id}"
        adata = self.storage.load_anndata(source_snapshot.snapshot_path)

        overall_results = {}
        plot_paths = {}

        for model_name in request.model_names:
            logger.info(f"Applying CellTypist model: {model_name}")
            try:
                # Core annotation (In-place)
                res = run_celltypist_annotation(
                    adata=adata, model_name=model_name,
                    majority_voting=request.majority_voting,
                    target_cluster_col=request.target_cluster_col
                )

                if "error" in res:
                    overall_results[model_name] = {"error": res["error"]}
                    continue

                overall_results[model_name] = res.get("cluster_summary", {})

                # Resolve label column name used by CellTypist
                label_col = "majority_voting" if request.majority_voting else "predicted_labels"

                # Backup result to unique column to avoid overwrite in next iteration
                safe_name = model_name.replace('.pkl', '').replace(' ', '_').lower()
                unique_col = f"ct_{safe_name}"
                adata.obs[unique_col] = adata.obs[label_col]

                # Generate model-specific UMAP
                p_path = self._plot_annotated_umap(
                    request.project_id, request.dataset_id, assets_base_path,
                    adata, color_col=unique_col, title=f"CellTypist: {model_name}", suffix=f"_{safe_name}"
                )
                if p_path:
                    plot_paths[f"umap_{safe_name}"] = p_path

            except Exception as e:
                logger.error(f"Model {model_name} failed: {e}")
                overall_results[model_name] = {"error": str(e)}

        # Save assets and snapshot
        json_path = get_dataset_relative(request.project_id, request.dataset_id, assets_base_path,
                                         f"{snapshot_id}.json")
        self.storage.save_file(overall_results, json_path)

        saved_path = self.storage.save_anndata(adata, get_dataset_relative(request.project_id, request.dataset_id,
                                                                           SNAPSHOTS_SUBSPACE, f"{snapshot_id}.h5ad"))

        branch_name = ANNOTATION_BRANCH
        current_count = AnalysisSnapshotsDAO.count_snapshots_by_branch(dataset_id=request.dataset_id,
                                                                       branch_name=branch_name)
        current_count += 1
        self.snapshot_dao.create_snapshot(
            dataset_id=request.dataset_id, snapshot_id=snapshot_id, branch_name=ANNOTATION_BRANCH,
            snapshot_path=saved_path, parent_snapshot_id=source_snapshot.snapshot_id,
            thumbnail_json={"report": json_path, **plot_paths},
            snapshot_name=f"{branch_name} {current_count}",
            params_json=request.model_dump(), user_notes=f"Annotated with {len(request.model_names)} models."
        )

        return AnnotationResultDTO(
            snapshot_id=snapshot_id, snapshot_path=saved_path,
            cluster_id=list(overall_results.keys()), enrichment_results=overall_results,
            msg="Multi-model annotation complete."
        )

    def _plot_enrichment_bar(self, project_id: str, dataset_id: str, assets_base_path: str,
                             cluster_id: str, results: Dict[str, Any]) -> Optional[str]:

        try:
            category = "cellmarker" if "cellmarker" in results else list(results.keys())[0]
            if category not in results or not results[category].get("top_terms"):
                return None

            df = pd.DataFrame(results[category]["top_terms"][:10])
            df['log_p'] = df['Adjusted P-value'].apply(lambda x: -math.log10(x) if x > 0 else 50)

            fig, ax = plt.subplots(figsize=(10, 6))
            sns.barplot(data=df, x='log_p', y='Term', palette='viridis', ax=ax)
            ax.set_title(f"Top Enriched Terms - Cluster {cluster_id} ({category})")
            ax.set_xlabel("-log10(Adjusted P-value)")
            ax.set_ylabel("Term")
            fig.tight_layout()

            filename = generate_filename(f"enrich_plot_c{cluster_id}", "png")
            rel_path = get_dataset_relative(project_id, dataset_id, assets_base_path, filename)

            self.storage.save_file(fig, rel_path, dpi=600)
            plt.close(fig)
            return rel_path
        except Exception as e:
            logger.error(f"Failed to plot enrichment bar: {e}")
            return None

    def _plot_annotated_umap(self, project_id: str, dataset_id: str, assets_base_path: str,
                             adata: sc.AnnData, color_col: str, title: str, suffix: str = "") -> Optional[str]:
        """
        Generates an annotated UMAP.
        Filename is STRICTLY derived from 'suffix' -> {suffix}.png
        """
        try:
            if "X_umap" not in adata.obsm:
                logger.warning("X_umap not found. Skipping UMAP plot.")
                return None

            if color_col not in adata.obs:
                logger.warning(f"Column {color_col} not found for plotting.")
                return None

            fig, ax = plt.subplots(figsize=(10, 10))

            sc.pl.umap(
                adata,
                color=color_col,
                title=title,
                show=False,
                frameon=False,
                ax=ax,
                legend_loc="on data",
                legend_fontsize='small',
                legend_fontoutline=2
            )

            # Strict naming: No random hash, just the key + .png
            filename = f"{suffix}.png"

            rel_path = get_dataset_relative(project_id, dataset_id, assets_base_path, filename)

            self.storage.save_file(fig, rel_path, dpi=600)
            plt.close(fig)
            return rel_path
        except Exception as e:
            logger.error(f"Failed to plot UMAP for {suffix}: {e}")
            return None

    def _apply_cluster_colors_to_annotation(self, adata, source_cluster_col: str, target_annot_col: str):
        """
        Helper: Maps colors from the original numeric cluster ID to the new annotation label.
        Ensures consistency: If Cluster 0 is Blue, 'T-cells' (derived from 0) will be Blue.
        """
        try:
            # 1. Get Source Colors (e.g. leiden_colors)
            source_color_key = f"{source_cluster_col}_colors"
            if source_color_key not in adata.uns:
                # Force generate if missing
                sc.pl.umap(adata, color=source_cluster_col, show=False)

            if source_color_key in adata.uns and hasattr(adata.obs[source_cluster_col], 'cat'):
                source_cats = adata.obs[source_cluster_col].cat.categories
                source_colors = adata.uns[source_color_key]

                # Map: {'0': '#1f77b4', '1': '#ff7f0e', ...}
                id_to_color = dict(zip(source_cats, source_colors))

                # 2. Apply to Target (e.g. obs_gseapy_...)
                target_cats = adata.obs[target_annot_col].cat.categories
                new_colors = []

                for label in target_cats:
                    # Find the first cell with this label
                    mask = (adata.obs[target_annot_col] == label)
                    if mask.any():
                        # Find which original cluster ID this cell belonged to
                        original_id = adata.obs[source_cluster_col][mask].iloc[0]
                        # Inherit that cluster's color
                        new_colors.append(id_to_color.get(str(original_id), "#808080"))
                    else:
                        new_colors.append("#808080")  # Fallback grey

                # 3. Save to UNS so Scanpy uses it for plotting
                adata.uns[f"{target_annot_col}_colors"] = new_colors
                logger.info(f"Synced colors from {source_cluster_col} to {target_annot_col}")

        except Exception as e:
            logger.warning(f"Failed to sync colors: {e}")
