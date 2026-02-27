from datetime import datetime
from typing import List, Optional
from peewee import DoesNotExist, IntegrityError

from ..model.project_meta_model import ProjectMeta


class ProjectMetaDAO:
    """
    Data Access Object (DAO) for managing ProjectMeta database operations.
    Encapsulates all database interactions related to projects to ensure
    consistent data handling and abstraction from the business logic.
    """

    @staticmethod
    def create_project(
            project_id: str,
            project_name: str,
            project_path: str,
            description: Optional[str] = None,
            ext_info: Optional[str] = None
    ) -> Optional[ProjectMeta]:
        """
        Creates a new ProjectMeta record in the database.

        Args:
            project_id (str): Unique identifier for the project.
            project_name (str): Human-readable name for the project.
            project_path (str): project file root path.
            description (str, optional): Project description.
            ext_info (str, optional): Additional metadata in JSON string format.

        Returns:
            ProjectMeta: The created project instance if successful.
        """
        try:
            now = datetime.utcnow()
            project = ProjectMeta.create(
                project_id=project_id,
                project_name=project_name,
                project_path=project_path,
                description=description,
                ext_info=ext_info,
                is_deleted=False,  # Active by default
                create_time=now,
                update_time=now
            )
            return project
        except IntegrityError as e:
            print(f"[Error] Failed to create project '{project_id}': {e}")
            return None
        except Exception as e:
            print(f"[Error] Unexpected error creating project: {e}")
            return None

    @staticmethod
    def get_project_by_id(pk_id: int) -> Optional[ProjectMeta]:
        """
        Retrieves a single ACTIVE project by its primary key ID.
        """
        try:
            return ProjectMeta.get((ProjectMeta.id == pk_id) & (ProjectMeta.is_deleted == False))
        except DoesNotExist:
            print(f"[Warning] Active Project with PK {pk_id} not found.")
            return None

    @staticmethod
    def get_project_file_path(project_id: str) -> Optional[str]:
        """
        Retrieves file path for an ACTIVE project.
        """
        try:
            project = ProjectMeta.get((ProjectMeta.project_id == project_id) & (ProjectMeta.is_deleted == False))
            return project.project_path
        except DoesNotExist:
            print(f"[Warning] Project with ID '{project_id}' not found.")
            return None

    @staticmethod
    def get_project_by_project_id(project_id: str) -> Optional[ProjectMeta]:
        """
        Retrieves a single ACTIVE project by its unique business string identifier.
        """
        try:
            return ProjectMeta.get((ProjectMeta.project_id == project_id) & (ProjectMeta.is_deleted == False))
        except DoesNotExist:
            print(f"[Warning] Project with ID '{project_id}' not found.")
            return None

    @staticmethod
    def get_all_projects(include_deleted: bool = False) -> List[ProjectMeta]:
        """
        Retrieves projects from the database, ordered by update time.

        Args:
            include_deleted (bool): If True, returns both active and deleted projects.
        """
        try:
            query = ProjectMeta.select().order_by(ProjectMeta.update_time.desc())
            if not include_deleted:
                query = query.where(ProjectMeta.is_deleted == False)
            return list(query)
        except Exception as e:
            print(f"[Error] Failed to retrieve all projects: {e}")
            return []

    @staticmethod
    def update_project(
            pk_id: int,
            project_name: Optional[str] = None,
            description: Optional[str] = None,
            project_path: Optional[str] = None,
            ext_info: Optional[str] = None
    ) -> bool:
        """
        Updates an existing project record. Automatically updates the 'update_time'.
        """
        try:
            project = ProjectMeta.get_by_id(pk_id)

            changed = False
            if project_name is not None:
                project.project_name = project_name
                changed = True

            if project_path is not None:
                project.project_path = project_path
                changed = True

            if description is not None:
                project.description = description
                changed = True

            if ext_info is not None:
                project.ext_info = ext_info
                changed = True

            if changed:
                project.update_time = datetime.utcnow()
                project.save()
                return True

            return True

        except DoesNotExist:
            print(f"[Error] Cannot update: Project ID {pk_id} does not exist.")
            return False
        except Exception as e:
            print(f"[Error] Failed to update project {pk_id}: {e}")
            return False

    @staticmethod
    def delete_project(pk_id: int, hard_delete: bool = False) -> bool:
        """
        Deletes a project record. Defaults to SOFT DELETE.

        Args:
            pk_id (int): Primary key of the project.
            hard_delete (bool): If True, physically removes the record.
        """
        try:
            project = ProjectMeta.get_by_id(pk_id)
            if hard_delete:
                project.delete_instance()
                print(f"[Info] Project ID {pk_id} HARD deleted successfully.")
            else:
                project.is_deleted = True
                project.update_time = datetime.utcnow()
                project.save()
                print(f"[Info] Project ID {pk_id} SOFT deleted successfully.")
            return True
        except DoesNotExist:
            print(f"[Warning] Cannot delete: Project ID {pk_id} does not exist.")
            return False
        except Exception as e:
            print(f"[Error] Failed to delete project {pk_id}: {e}")
            return False
