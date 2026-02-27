import matplotlib.pyplot as plt
from fastapi import Depends
from loguru import logger
import matplotlib.colors as mcolors

from dto.parameter.cluster_parameter import ClusterParameter
from dto.request.merge_clusters_request import MergeClustersRequest
from dto.request.run_clustering_request import RunClusteringRequest
from dto.request.sub_clustering_request import SubClusteringRequest
from dto.response.clustering_result_dto import ClusteringResultDTO

from infrastructure.database.dao.analysis_snapshots_dao import AnalysisSnapshotsDAO
from infrastructure.filesystem.storage import AssetStorage
from infrastructure.filesystem.constants.filesystem_constants import SNAPSHOTS_SUBSPACE, CLUSTERING_SUBSPACE
import lotus.clustering as clustering
import lotus.visualization as visualization
from lotus.core_analysis import corespect
from service.helper.project_id_resolver_service import ProjectIdResolverService
from service.helper.service_constant import PREPROCESSING_BRANCH, CLUSTERING_BRANCH
from util.id_generate_utils import generate_business_id
from util.path_utils import get_dataset_relative
from util.snapshot_utils import resolve_source_snapshot


class ClusteringService:
    def __init__(
            self,
            snapshot_dao: AnalysisSnapshotsDAO = Depends(),
            storage: AssetStorage = Depends(),
            project_id_resolver: ProjectIdResolverService = Depends(),
    ):
        self.snapshot_dao = snapshot_dao
        self.storage = storage
        self.project_id_resolver = project_id_resolver

    def _get_cluster_colors(self, adata, cluster_key: str) -> dict:
        """
        Helper function to extract cluster colors from AnnData.
        Returns a dictionary mapping cluster ID (str) to Hex Color (str).
        e.g., {"0": "#1f77b4", "1": "#ff7f0e"}
        """
        try:
            # 1. Check if color key exists in uns (Scanpy convention: {key}_colors)
            color_key = f"{cluster_key}_colors"
            if color_key not in adata.uns:
                logger.warning(f"Color key '{color_key}' not found in adata.uns. Returning empty map.")
                return {}

            # 2. Get categories (Cluster IDs)
            # Ensure categories are strings
            categories = adata.obs[cluster_key].cat.categories.astype(str)

            # 3. Get colors (Scanpy stores them as a sequence corresponding to categories)
            colors = adata.uns[color_key]

            # 4. Map them
            # Convert any non-hex colors to hex just in case
            color_map = {}
            for cat, color in zip(categories, colors):
                color_map[cat] = mcolors.to_hex(color)

            return color_map

        except Exception as e:
            logger.warning(f"Failed to extract cluster colors: {e}")
            return {}

    def _get_cluster_counts(self, adata, cluster_key: str) -> dict:
        """
        Helper function to extract cluster counts from AnnData.
        Returns a dictionary mapping cluster ID (str) to Count (int).
        e.g., {"0": 123, "1": 456}
        """
        try:
            counts = adata.obs[cluster_key].value_counts()
            # Convert index to string and values to int
            return {str(k): int(v) for k, v in counts.items()}
        except Exception as e:
            logger.warning(f"Failed to extract cluster counts: {e}")
            return {}

    def _calculate_stats(self, adata, cluster_key: str):
        """
        Helper to calculate common cluster statistics for thumbnail_json.
        """
        cluster_counts = adata.obs[cluster_key].value_counts()
        n_clusters = len(cluster_counts)
        min_size = int(cluster_counts.min())
        max_size = int(cluster_counts.max())
        median_size = int(cluster_counts.median())

        return n_clusters, min_size, max_size, median_size

    def run_clustering(self, request: RunClusteringRequest) -> ClusteringResultDTO:
        """
        Executes Clustering (Leiden/Louvain/CoreSpect) and UMAP for a single resolution.
        """
        # Resolve project id
        request.project_id = self.project_id_resolver.resolve_project_id(request.project_id, request.dataset_id,
                                                                         request.snapshot_id)

        # Validation
        if request.resolution is None:
            raise ValueError("Resolution must be provided for clustering.")

        snapshot_id = generate_business_id(f's_cluster_{request.method}')
        assets_base_path = f"{CLUSTERING_SUBSPACE}/{snapshot_id}"

        # Resolve Source Snapshot
        source_snapshot = resolve_source_snapshot(
            request.dataset_id,
            request.snapshot_id,
            target_branch=PREPROCESSING_BRANCH,
            snapshot_dao=self.snapshot_dao
        )

        # Load Data
        adata = self.storage.load_anndata(source_snapshot.snapshot_path)

        # Check for neighborhood graph
        if 'neighbors' not in adata.uns and request.method != 'cplearn':
            raise ValueError("Neighborhood graph not found. Please run 'Build Neighborhood Graph' first.")

        logger.info(f"Running {request.method} clustering with resolution={request.resolution}...")

        # CoreSpect placeholders
        actual_q = None
        actual_r = None

        # Directly use the method name as the key (e.g., 'leiden')
        key_added = request.method

        # Execute Clustering
        if request.method == 'leiden':
            clustering.leiden(
                adata,
                resolution=request.resolution,
                key_added=key_added
            )
        elif request.method == 'louvain':
            clustering.louvain(
                adata,
                resolution=request.resolution,
                key_added=key_added
            )
        elif request.method == 'cplearn':
            n_obs = adata.n_obs
            default_q = 20
            default_r = 10
            actual_q = default_q
            actual_r = default_r

            if n_obs <= default_q:
                actual_q = 2
            if n_obs <= default_r * 2:
                actual_r = 2

            corespect(
                adata,
                resolution=request.resolution,
                key_added=key_added,
                use_rep="X_pca",
                copy=False,
                q=actual_q,
                r=actual_r
            )
        else:
            raise ValueError(f"Unsupported clustering method: {request.method}")

        # Post-Processing
        if key_added not in adata.obs:
            raise RuntimeError(f"Clustering failed. Key {key_added} was not generated.")

        if hasattr(adata.obs[key_added], "cat"):
            adata.obs[key_added] = adata.obs[key_added].cat.remove_unused_categories()

        # UMAP Calculation
        if "X_umap" not in adata.obsm:
            logger.info("Computing UMAP embeddings...")
            clustering.umap(adata)

        # Visualization: Primary UMAP
        # This step implicitly generates the colors in adata.uns[{key}_colors]
        plot_key = "umap_plot"
        visualization.umap(
            adata,
            color=key_added,
            legend_loc='on data',
            title=f"{request.method} (res={request.resolution})",
            show=False
        )
        fig_primary = plt.gcf()

        primary_path = self.storage.save_file(
            fig_primary,
            get_dataset_relative(request.project_id, request.dataset_id, assets_base_path, f"{plot_key}.png"),
            dpi=600
        )
        plt.close(fig_primary)

        n_clusters, min_size, max_size, median_size = self._calculate_stats(adata, key_added)

        # New Feature: Extract color mapping
        cluster_colors = self._get_cluster_colors(adata, key_added)
        cluster_counts = self._get_cluster_counts(adata, key_added)

        thumbnail_map = {
            "summary": {
                "method": request.method,
                "resolution": request.resolution,
                "n_clusters": n_clusters,
                "n_cells_total": int(adata.n_obs)
            },
            "cluster_stats": {
                "min_cells_per_cluster": min_size,
                "max_cells_per_cluster": max_size,
                "median_cells_per_cluster": median_size
            },
            "cluster_ids": cluster_colors,  # Store ID -> Color mapping
            "cluster_counts": cluster_counts,
            plot_key: primary_path
        }

        # Prepare Parameter Object for DB
        cluster_parameter = ClusterParameter(
            method=request.method,
            resolution=request.resolution,
            run_hierarchical=request.run_hierarchical,
            qnn_size=actual_q,
            neighbor_radius=actual_r
        )

        # Hierarchical Clustering (Optional)
        dendrogram_path = None
        if request.run_hierarchical:
            if n_clusters < 2:
                logger.warning(
                    f"Hierarchical clustering skipped: Only {n_clusters} cluster(s) found."
                )
            else:
                logger.info("Computing hierarchical clustering dendrogram...")
                try:
                    clustering.dendrogram(adata, groupby=key_added)

                    # Plot Dendrogram
                    dend_plot_key = "dendrogram_plot"
                    visualization.dendrogram(adata, key_added, show=False)
                    fig_dend = plt.gcf()

                    dendrogram_path = self.storage.save_file(
                        fig_dend,
                        get_dataset_relative(request.project_id, request.dataset_id, assets_base_path,
                                             f"{dend_plot_key}.png"),
                        dpi=600
                    )
                    thumbnail_map[dend_plot_key] = dendrogram_path
                    plt.close(fig_dend)
                except Exception as e:
                    logger.warning(f"Dendrogram calculation failed: {e}")

        # Persistence
        file_name = f"{snapshot_id}.h5ad"
        relative_key = get_dataset_relative(request.project_id, request.dataset_id, SNAPSHOTS_SUBSPACE, file_name)

        saved_path = self.storage.save_anndata(adata, relative_key)

        branch_name = CLUSTERING_BRANCH
        current_count = AnalysisSnapshotsDAO.count_snapshots_by_branch(dataset_id=request.dataset_id,
                                                                       branch_name=branch_name)
        current_count += 1
        # Update Database
        self.snapshot_dao.create_snapshot(
            dataset_id=request.dataset_id,
            snapshot_id=snapshot_id,
            branch_name=branch_name,
            snapshot_path=saved_path,
            snapshot_name=f"{request.method}{request.resolution} {current_count}",
            parent_snapshot_id=source_snapshot.snapshot_id,
            params_json=cluster_parameter.model_dump(exclude_none=True),
            thumbnail_json=thumbnail_map,
            user_notes=f"Clustering ({request.method}) performed at res={request.resolution}. Found {n_clusters} clusters."
        )

        return ClusteringResultDTO(
            snapshot_id=snapshot_id,
            snapshot_path=saved_path,
            umap_plot_path=primary_path,
            dendrogram_path=dendrogram_path,
            clusters_summary={key_added: n_clusters},
            msg="Clustering and UMAP projection complete."
        )

    def merge_clusters(self, request: MergeClustersRequest) -> ClusteringResultDTO:
        """
        Merges specified clusters into a single new cluster, regenerates the UMAP,
        and re-calculates the dendrogram for the new grouping.

        Supports overwriting the current snapshot or creating a new branch.
        """
        # Resolve project id
        request.project_id = self.project_id_resolver.resolve_project_id(
            request.project_id, request.dataset_id, request.snapshot_id
        )

        # Load the specific source snapshot
        source_snapshot = self.snapshot_dao.get_snapshot_by_business_id(request.snapshot_id)
        if not source_snapshot:
            raise ValueError(f"Snapshot {request.snapshot_id} not found.")

        # Load AnnData
        adata = self.storage.load_anndata(source_snapshot.snapshot_path)

        # Validate Column exists
        if request.method not in adata.obs:
            # Check params if the key is missing (fallback logic)
            params = source_snapshot.params_json or {}
            cluster_key = params.get("method", "leiden")
            request.method = cluster_key
            if cluster_key not in adata.obs:
                raise ValueError(f"Clustering column '{cluster_key}' not found in data.")

        # Perform Merge (Modify adata.obs)
        clusters_to_merge_str = [str(c) for c in request.clusters_to_merge]

        # If a new label isn't provided, default to the smallest ID in the merge group
        if request.new_label is None or request.new_label.strip() == '' or request.new_label == 'None':
            try:
                # Convert to integers to find the minimum
                numeric_ids = [int(c) for c in clusters_to_merge_str]
                if not numeric_ids:
                    raise ValueError("Cannot merge an empty list of clusters.")
                # The new label is the smallest ID from the merged group
                request.new_label = str(min(numeric_ids))
            except (ValueError, TypeError):
                # Fallback for non-numeric labels, which shouldn't happen with cluster IDs
                logger.warning("Could not determine minimum numeric ID for merged cluster. Using concatenation as fallback.")
                request.new_label = '-'.join(sorted(clusters_to_merge_str))

        logger.info(f"Merging clusters {clusters_to_merge_str} into '{request.new_label}'...")

        obs_series = adata.obs[request.method].astype(str)
        existing_clusters = set(obs_series.unique())

        # Track valid and missing clusters
        valid_clusters_to_merge = []
        missing_clusters = []

        for c in clusters_to_merge_str:
            if c in existing_clusters:
                valid_clusters_to_merge.append(c)
            else:
                missing_clusters.append(c)
                logger.warning(f"Cluster '{c}' to merge not found in data. Skipping.")

        # If NO valid clusters were found, abort immediately.
        # This ensures the downstream service gets a clear error instead of a fake success.
        if not valid_clusters_to_merge:
            error_msg = (
                f"Merge aborted: None of the requested clusters {missing_clusters} exist in the current dataset. "
                f"Available clusters are: {sorted(list(existing_clusters))}. "
                "Please verify the cluster IDs from the previous step."
            )
            raise ValueError(error_msg)

        # Apply replacement only for valid clusters
        replace_map = {old: request.new_label for old in valid_clusters_to_merge}
        new_obs = obs_series.replace(replace_map)

        adata.obs[request.method] = new_obs.astype("category")

        # Recalculate stats
        n_clusters, min_size, max_size, median_size = self._calculate_stats(adata, request.method)

        # Visualization 1: Regenerate UMAP (This generates new colors)
        visualization.umap(
            adata,
            color=request.method,
            legend_loc='on data',
            title=f"{request.method} (Merged)",
            show=False
        )
        fig_primary = plt.gcf()

        # Determine Target Snapshot ID & Path
        if request.overwrite:
            target_snapshot_id = request.snapshot_id
            assets_base_path = f"{CLUSTERING_SUBSPACE}/{target_snapshot_id}"
            logger.info(f"Overwriting existing snapshot: {target_snapshot_id}")
        else:
            target_snapshot_id = generate_business_id(f's_merge')
            assets_base_path = f"{CLUSTERING_SUBSPACE}/{target_snapshot_id}"
            logger.info(f"Creating new snapshot: {target_snapshot_id}")

        primary_filename = "umap_merged_plot.png" if not request.overwrite else "umap_plot.png"
        plot_key_umap = "umap_plot"

        primary_path = self.storage.save_file(
            fig_primary,
            get_dataset_relative(request.project_id, request.dataset_id, assets_base_path, primary_filename),
            dpi=600
        )
        plt.close(fig_primary)

        cluster_colors = self._get_cluster_colors(adata, request.method)
        cluster_counts = self._get_cluster_counts(adata, request.method)

        # Visualization 2: Regenerate Dendrogram
        dendrogram_path = None
        plot_key_dend = "dendrogram_plot"

        if n_clusters >= 2:
            logger.info("Re-calculating hierarchical clustering dendrogram for merged groups...")
            try:
                clustering.dendrogram(adata, groupby=request.method)

                visualization.dendrogram(adata, request.method, show=False)
                fig_dend = plt.gcf()

                dendrogram_filename = "dendrogram_merged_plot.png" if not request.overwrite else "dendrogram_plot.png"

                dendrogram_path = self.storage.save_file(
                    fig_dend,
                    get_dataset_relative(request.project_id, request.dataset_id, assets_base_path,
                                         dendrogram_filename),
                    dpi=600
                )
                plt.close(fig_dend)
            except Exception as e:
                logger.warning(f"Dendrogram update failed: {e}")
        else:
            logger.info("Skipping dendrogram: Merged result has fewer than 2 clusters.")

        # Prepare FULL thumbnail map structure
        thumbnail_update = {
            "summary": {
                "method": request.method,
                "n_clusters": n_clusters,
                "n_cells_total": int(adata.n_obs),
                "resolution": source_snapshot.params_json.get("resolution"),
                "action": "merged"
            },
            "cluster_stats": {
                "min_cells_per_cluster": min_size,
                "max_cells_per_cluster": max_size,
                "median_cells_per_cluster": median_size
            },
            "cluster_ids": cluster_colors,
            "cluster_counts": cluster_counts,
            plot_key_umap: primary_path
        }
        if dendrogram_path:
            thumbnail_update[plot_key_dend] = dendrogram_path

        # Construct success message (include warnings about skipped clusters)
        success_msg = f"Clusters merged. New group '{request.new_label}' created."
        if missing_clusters:
            success_msg += f" Note: Clusters {missing_clusters} were skipped as they were not found."

        # Persistence
        if request.overwrite:
            file_name = f"{target_snapshot_id}.h5ad"
            relative_key = get_dataset_relative(request.project_id, request.dataset_id, SNAPSHOTS_SUBSPACE, file_name)
            saved_path = self.storage.save_anndata(adata, relative_key)

            self.snapshot_dao.update_snapshot(
                pk_id=source_snapshot.id,
                user_notes=f"Merged clusters {valid_clusters_to_merge} into '{request.new_label}'. Overwritten.",
                params_update={"merge_action": request.model_dump()},
                thumbnail_json=thumbnail_update
            )
        else:
            file_name = f"{target_snapshot_id}.h5ad"
            relative_key = get_dataset_relative(request.project_id, request.dataset_id, SNAPSHOTS_SUBSPACE, file_name)
            saved_path = self.storage.save_anndata(adata, relative_key)

            branch_name = CLUSTERING_BRANCH
            current_count = AnalysisSnapshotsDAO.count_snapshots_by_branch(dataset_id=request.dataset_id,
                                                                           branch_name=branch_name)
            current_count += 1
            self.snapshot_dao.create_snapshot(
                dataset_id=request.dataset_id,
                snapshot_id=target_snapshot_id,
                branch_name=CLUSTERING_BRANCH,
                snapshot_path=saved_path,
                snapshot_name=f"{request.method} Merged {current_count}",
                parent_snapshot_id=request.snapshot_id,
                params_json=request.model_dump(),
                thumbnail_json=thumbnail_update,
                user_notes=f"Merged clusters {valid_clusters_to_merge} into '{request.new_label}'."
            )

        return ClusteringResultDTO(
            snapshot_id=target_snapshot_id,
            snapshot_path=saved_path,
            umap_plot_path=primary_path,
            dendrogram_path=dendrogram_path,
            clusters_summary={request.method: n_clusters},
            msg=success_msg
        )

    def run_sub_clustering(self, request: SubClusteringRequest) -> ClusteringResultDTO:
        """
        Runs clustering only on the specified target clusters and merges the results
        back into the main dataset by assigning new, auto-incrementing integer IDs.
        """
        request.project_id = self.project_id_resolver.resolve_project_id(
            request.project_id, request.dataset_id, request.snapshot_id
        )

        source_snapshot = self.snapshot_dao.get_snapshot_by_business_id(request.snapshot_id)
        if not source_snapshot:
            raise ValueError(f"Snapshot {request.snapshot_id} not found.")

        adata = self.storage.load_anndata(source_snapshot.snapshot_path)

        # Validate the clustering key
        cluster_key = request.method
        if cluster_key not in adata.obs:
            params = source_snapshot.params_json or {}
            cluster_key = params.get("method", "leiden")
            request.method = cluster_key
            if cluster_key not in adata.obs:
                raise ValueError(f"Clustering column '{cluster_key}' not found in data.")

        # Validation: Check targets exist
        available_clusters = set(adata.obs[cluster_key].astype(str).unique())
        for t in request.target_clusters:
            if t not in available_clusters:
                raise ValueError(f"Target cluster '{t}' not found in dataset.")

        logger.info(
            f"Sub-clustering targets {request.target_clusters} using {request.method} (res={request.resolution})..."
        )

        # Create subset and re-cluster
        subset_mask = adata.obs[cluster_key].isin(request.target_clusters)
        adata_sub = adata[subset_mask].copy()

        sub_key = f"{cluster_key}_sub"
        if request.method == 'leiden':
            clustering.leiden(adata_sub, resolution=request.resolution, key_added=sub_key)
        elif request.method == 'louvain':
            clustering.louvain(adata_sub, resolution=request.resolution, key_added=sub_key)
        elif request.method == 'cplearn':
            n_obs = adata_sub.n_obs
            default_q = 20
            default_r = 10
            actual_q = default_q
            actual_r = default_r

            if n_obs <= default_q:
                actual_q = 2
            if n_obs <= default_r * 2:
                actual_r = 2

            corespect(
                adata_sub,
                resolution=request.resolution,
                key_added=sub_key,
                use_rep="X_pca",
                copy=False,
                q=actual_q,
                r=actual_r
            )
        else:
            raise ValueError(f"Method {request.method} not supported for sub-clustering.")

        # Merge results back by assigning new integer IDs
        adata.obs[cluster_key] = adata.obs[cluster_key].astype(str)

        # Find the maximum existing numeric cluster ID from the clusters NOT being sub-clustered.
        original_labels = adata.obs[~subset_mask][cluster_key].unique().tolist()
        numeric_labels = [int(l) for l in original_labels if l.isdigit()]
        max_id = max(numeric_labels) if numeric_labels else -1
        next_id = max_id + 1

        # Get the new sub-cluster labels generated in the subset (e.g., '0', '1', '2')
        # Ensure they are sorted to have a deterministic mapping to new global IDs
        sub_cluster_cats = sorted(adata_sub.obs[sub_key].cat.categories, key=lambda x: int(x) if x.isdigit() else x)

        # Create a mapping from the new sub-cluster labels to new global IDs (e.g., '0' -> '11', '1' -> '12')
        sub_id_map = {sub_id: str(next_id + i) for i, sub_id in enumerate(sub_cluster_cats)}
        logger.info(f"Mapping new sub-clusters to global IDs starting from {next_id}: {sub_id_map}")

        # Create a series with the new global labels for the subset of cells
        new_labels_for_subset = adata_sub.obs[sub_key].map(sub_id_map)

        # Update the original adata object for the sub-clustered cells
        adata.obs.loc[subset_mask, cluster_key] = new_labels_for_subset.values
        adata.obs[cluster_key] = adata.obs[cluster_key].astype("category")

        # Recalculate stats
        n_clusters, min_size, max_size, median_size = self._calculate_stats(adata, cluster_key)

        # Determine Target Snapshot ID & Path
        if request.overwrite:
            target_snapshot_id = request.snapshot_id
            assets_base_path = f"{CLUSTERING_SUBSPACE}/{target_snapshot_id}"
            logger.info(f"Overwriting existing snapshot: {target_snapshot_id}")
        else:
            target_snapshot_id = generate_business_id(f's_sub_{cluster_key}')
            assets_base_path = f"{CLUSTERING_SUBSPACE}/{target_snapshot_id}"
            logger.info(f"Creating new snapshot: {target_snapshot_id}")

        # Visualization 1: Regenerate UMAP (This generates new colors)
        plot_key_umap = "umap_plot"
        visualization.umap(
            adata,
            color=cluster_key,
            legend_loc='on data',
            title=f"{request.method} (Sub-clustered)",
            show=False
        )
        fig_primary = plt.gcf()

        primary_filename = "umap_subcluster_plot.png" if not request.overwrite else "umap_plot.png"

        primary_path = self.storage.save_file(
            fig_primary,
            get_dataset_relative(request.project_id, request.dataset_id, assets_base_path, primary_filename),
            dpi=600
        )
        plt.close(fig_primary)

        # New Feature: Extract color mapping
        cluster_colors = self._get_cluster_colors(adata, cluster_key)
        cluster_counts = self._get_cluster_counts(adata, cluster_key)

        # Visualization 2: Regenerate Dendrogram
        dendrogram_path = None
        plot_key_dend = "dendrogram_plot"

        if n_clusters >= 2:
            logger.info("Computing hierarchical clustering dendrogram for sub-clustered data...")
            try:
                clustering.dendrogram(adata, groupby=cluster_key)

                visualization.dendrogram(adata, cluster_key, show=False)
                fig_dend = plt.gcf()

                dend_filename = "dendrogram_subcluster_plot.png" if not request.overwrite else "dendrogram_plot.png"

                dendrogram_path = self.storage.save_file(
                    fig_dend,
                    get_dataset_relative(request.project_id, request.dataset_id, assets_base_path, dend_filename),
                    dpi=600
                )
                plt.close(fig_dend)
            except Exception as e:
                logger.warning(f"Dendrogram calculation failed: {e}")
        else:
            logger.info("Skipping dendrogram: Result has fewer than 2 clusters.")

        # Prepare FULL thumbnail map structure
        thumbnail_update = {
            "summary": {
                "method": request.method,
                "resolution": request.resolution,
                "n_clusters": n_clusters,
                "n_cells_total": int(adata.n_obs),
                "action": "sub-clustered"
            },
            "cluster_stats": {
                "min_cells_per_cluster": min_size,
                "max_cells_per_cluster": max_size,
                "median_cells_per_cluster": median_size
            },
            "cluster_ids": cluster_colors,
            "cluster_counts": cluster_counts,
            plot_key_umap: primary_path
        }
        if dendrogram_path:
            thumbnail_update[plot_key_dend] = dendrogram_path

        # Persistence
        if request.overwrite:
            file_name = f"{target_snapshot_id}.h5ad"
            relative_key = get_dataset_relative(request.project_id, request.dataset_id, SNAPSHOTS_SUBSPACE, file_name)
            saved_path = self.storage.save_anndata(adata, relative_key)

            self.snapshot_dao.update_snapshot(
                pk_id=source_snapshot.id,
                user_notes=f"Sub-clustered {request.target_clusters} (res={request.resolution}). Overwritten.",
                params_update={"sub_clustering": request.model_dump()},
                thumbnail_json=thumbnail_update
            )
        else:
            file_name = f"{target_snapshot_id}.h5ad"
            relative_key = get_dataset_relative(request.project_id, request.dataset_id, SNAPSHOTS_SUBSPACE, file_name)
            saved_path = self.storage.save_anndata(adata, relative_key)

            branch_name = CLUSTERING_BRANCH
            current_count = AnalysisSnapshotsDAO.count_snapshots_by_branch(dataset_id=request.dataset_id,
                                                                           branch_name=branch_name)
            current_count += 1
            self.snapshot_dao.create_snapshot(
                dataset_id=request.dataset_id,
                snapshot_id=target_snapshot_id,
                branch_name=CLUSTERING_BRANCH,
                snapshot_name=f"{request.method} Sub {current_count}",
                snapshot_path=saved_path,
                parent_snapshot_id=source_snapshot.snapshot_id,
                params_json=request.model_dump(),
                thumbnail_json=thumbnail_update,
                user_notes=f"Sub-clustered {request.target_clusters} from {source_snapshot.snapshot_id}."
            )

        return ClusteringResultDTO(
            snapshot_id=target_snapshot_id,
            snapshot_path=saved_path,
            umap_plot_path=primary_path,
            dendrogram_path=dendrogram_path,
            clusters_summary={cluster_key: n_clusters},
            msg="Sub-clustering complete."
        )
