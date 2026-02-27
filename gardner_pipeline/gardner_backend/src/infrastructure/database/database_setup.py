import json

from infrastructure.database.connection import get_default_db_manager
from infrastructure.database.model.agent_message_model import AgentMessage
from infrastructure.database.model.agent_session_model import AgentSession
from infrastructure.database.model.annotation_method import AnnotationMethod
from infrastructure.database.model.project_meta_model import ProjectMeta
from infrastructure.database.model.dataset_model import Dataset
from infrastructure.database.model.analysis_snapshots_model import AnalysisSnapshot
from infrastructure.database.model.user_preference import UserPreference


def db_setup(ensure_schema: bool = True):
    """
    Handles the full database initialization workflow:
    connecting, binding the proxy, and creating tables.
    """
    # 1. Acquire the database manager
    manager = get_default_db_manager()

    # 2. Bind the proxy to the actual database instance
    manager.initialize_proxy()

    # 3. Get the database connection and ensure it is open
    db = manager.get_connection()
    if db.is_closed():
        db.connect()

    # 4. Centralized table creation logic
    if ensure_schema:
        # Explicitly list all models that require table creation
        models_to_create = [
            ProjectMeta,
            Dataset,
            AnalysisSnapshot,
            AgentSession,
            AgentMessage,
            AnnotationMethod,
            UserPreference
        ]
        db.create_tables(models_to_create, safe=True)
        print(f"[DB Setup] Verified tables: {[m.__name__ for m in models_to_create]}")

    # 5. Initialize Default Data
    initialize_default_data()
    return db


def initialize_default_data():
    """
    Main entry point to seed database with initial data.
    """
    print("[DB Init] Checking default data integrity...")
    _init_default_project()
    _init_annotation_methods()
    _init_default_preferences()


def _init_default_project():
    """
    Ensures the system's default project ('p_default') and its initial
    chat session exist in the database upon startup.
    """
    try:
        # Ensure the default project exists
        # 'get_or_create' avoids IntegrityErrors if 'p_default' is already present
        project, p_created = ProjectMeta.get_or_create(
            project_id="p_default",
            defaults={
                "project_name": "Default Project",
                "project_path": "projects/p_default",
                "description": "System default project initialized automatically.",
                "is_deleted": False
            }
        )
        if p_created:
            print(f"[DB Init] Created default project: {project.project_name}")

        # Ensure a default chat session exists for this project
        # Linked via the 'project' Foreign Key (using the project object)
        session_id = "session_default_001"
        session, s_created = AgentSession.get_or_create(
            session_id=session_id,
            defaults={
                "project": project,  # Associates with the project instance
                "agent_name": "BioAgent-Pro",
                "session_name": "Default Analysis Session",
                "status": "INIT"
            }
        )

        if s_created:
            print(f"[DB Init] Created default session: {session.session_id} for {project.project_id}")

    except Exception as e:
        print(f"[DB Init] Error during default project/session initialization: {e}")


def _init_annotation_methods():
    """
    Load annotation methods from JSON files.
    Uses 'INSERT OR IGNORE' logic to skip duplicates and ensure successful startup.
    """
    # 1. Check if table is already populated (Optional optimization)
    # You can remove this check if you want to always try adding new records from JSON
    if AnnotationMethod.select().count() > 0:
        # print("[DB Init] Annotation methods already populated. Skipping.")
        return

    print("[DB Init] Seeding annotation methods from JSON files...")

    methods_to_insert = []

    from pathlib import Path
    current_file = Path(__file__).resolve().parent
    # Use relative paths resolved safely
    CELLTYPIST_JSON_PATH = current_file / "resources" / "celltypist_models.json"
    GSEAPY_JSON_PATH = current_file / "resources" / "gseapy_libraries.json"

    # 2. Process CellTypist Data
    import os
    if os.path.exists(CELLTYPIST_JSON_PATH):
        try:
            with open(CELLTYPIST_JSON_PATH, 'r', encoding='utf-8') as f:
                ct_data = json.load(f)

            for item in ct_data:
                methods_to_insert.append({
                    "method_name": item.get("model"),
                    "type": "celltypist",
                    "description": item.get("description"),
                    "species": "Human",
                    "organ": None
                })
        except Exception as e:
            print(f"[DB Init] Failed to load CellTypist JSON: {e}")
    else:
        print(f"[DB Init] Warning: CellTypist resource not found at {CELLTYPIST_JSON_PATH}")

    # 3. Process GSEApy Data
    if os.path.exists(GSEAPY_JSON_PATH):
        try:
            with open(GSEAPY_JSON_PATH, 'r', encoding='utf-8') as f:
                gp_data = json.load(f)
                records = gp_data.get('records', []) if isinstance(gp_data, dict) else gp_data

            for item in records:
                organs = item.get("organs", [])
                organ_str = ", ".join(organs) if isinstance(organs, list) else str(organs)
                if not organ_str:
                    organ_str = None

                methods_to_insert.append({
                    "method_name": item.get("name"),
                    "type": "gseapy",
                    "description": item.get("description"),
                    "species": item.get("species", "Human"),
                    "organ": organ_str
                })
        except Exception as e:
            print(f"[DB Init] Failed to load GSEApy JSON: {e}")
    else:
        print(f"[DB Init] Warning: GSEApy resource not found at {GSEAPY_JSON_PATH}")

    # 4. Safe Bulk Insert
    if methods_to_insert:
        try:
            with get_default_db_manager().get_connection().atomic():
                # Use .on_conflict_ignore() to skip rows that violate the UNIQUE constraint
                count = (AnnotationMethod
                         .insert_many(methods_to_insert)
                         .on_conflict_ignore()
                         .execute())

            print(f"[DB Init] Batch insert complete. {count} new records added (duplicates ignored).")

        except Exception as e:
            print(f"[DB Init] Unexpected error during annotation method seeding: {e}")


def _init_default_preferences():
    """
    Seeds the 'UserPreference' table with a standard default configuration.
    This ensures the frontend always has a baseline setting to load.
    """
    # 1. Define standard pipeline parameters
    standard_settings = {
        "preprocessing_params": {
            "skip_qc_calculation": False,
            "skip_qc_filter": False,
            "organism": "Human",
            "min_genes": 300,
            "min_cells": 5,
            "pct_mt_max": 15.0,
            "max_counts": 40000,

            "skip_hvg": False,
            "n_top_genes": 2500,
            "flavor": "seurat",
            "target_sum": 10000.0,

            "skip_pca": False,
            "n_comps": 50,
            "svd_solver": "arpack",

            "skip_neighbors": False,
            "n_neighbors": 20,
            "n_pcs": 30
        },

        "clustering_params": {
            "method": "leiden",
            "resolution": 0.5,
            "run_hierarchical": True
        },

        "deg_params": {
            "method": "wilcoxon",
            "n_top_genes": 10
        },

        "annotation_params": {
            "categories": [
                "CellMarker_2024",
                "BioPlanet_2019",
                "NCI-Nature_2016",
                "GO_Biological_Process_2023"
            ],
            "top_n_genes": 50,
            "model_names": ["Immune_All_Low.pkl"],
            "majority_voting": True,
            "target_cluster_col": "leiden"
        }
    }

    try:
        # 2. Insert if not exists (using 'Default Settings' as the key)
        pref, created = UserPreference.get_or_create(
            preference_name="Default Settings",
            defaults={
                "settings": standard_settings,
                "is_deleted": False
            }
        )

        if created:
            print(f"[DB Init] Created default user preference: {pref.preference_name}")
        else:
            print(f"[DB Init] Default preference '{pref.preference_name}' already exists.")

    except Exception as e:
        print(f"[DB Init] Failed to seed default preferences: {e}")
