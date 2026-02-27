from peewee import (
    AutoField,
    CharField,
    DateTimeField,
    TextField,
    ForeignKeyField,
    BooleanField
)
from datetime import datetime
from .base_model import BaseModel
from .project_meta_model import ProjectMeta
from playhouse.sqlite_ext import JSONField


class Dataset(BaseModel):
    """
    Dataset table for managing raw single-cell data files.

    Links uploaded data files to projects and tracks their physical
    storage location, biological metadata, and statistics.

    Attributes:
        id: Auto-increment primary key
        dataset_id: Unique business identifier
        project_id: Foreign key to ProjectMeta
        dataset_path: Path relative to workspace root
        dataset_name: User-defined display name
        organism: Species information (e.g., Human, Mouse)
        tissue_type: Tissue origin (e.g., PBMC, Liver)
        details: Additional descriptive details about the dataset
        uploaded_time: ISO8601 timestamp of upload
        is_deleted: Soft deletion flag (0/False = Active, 1/True = Deleted)
        ext_info: Reserved JSON field for extensibility
    """

    # Primary Key
    id = AutoField(
        primary_key=True,
        help_text="Auto-increment primary key"
    )

    # dataset business id
    dataset_id = CharField(
        max_length=128,
        null=False,
        unique=True,
        help_text="dataset business id, naming rule: dataset_<timestamp>_<original_name>_<4_digits_random_suffix>"
    )

    # Foreign Key links to ProjectMeta's project id
    project_id = ForeignKeyField(
        ProjectMeta,
        field='project_id',
        backref='datasets',
        on_delete='CASCADE',
        column_name='project_id',
        help_text="Reference to the parent project's business ID"
    )

    dataset_path = TextField(
        null=False,
        help_text='Dataset relative file path in file workspace'
    )

    # Display Name
    dataset_name = CharField(
        max_length=255,
        null=False,
        help_text="User-defined dataset display name"
    )

    # Biological metadata (Moved from ProjectMeta)
    organism = CharField(
        max_length=50,
        null=True,
        index=True,
        help_text="Species (e.g., Human, Mouse)"
    )

    tissue_type = CharField(
        max_length=100,
        null=True,
        index=True,
        help_text="Tissue type (e.g., PBMC, Brain, Liver)"
    )

    # Detailed description
    details = TextField(
        null=True,
        help_text="Detailed description or notes about the dataset"
    )

    # Temporal Fields
    uploaded_time = DateTimeField(
        default=datetime.utcnow,
        null=False,
        help_text="Upload timestamp"
    )

    # Soft deletion flag
    is_deleted = BooleanField(
        default=False,
        null=False,
        index=True,
        help_text="Soft delete flag: 0 for active, 1 for deleted"
    )

    # Extensibility
    ext_info = JSONField(
        null=True,
        help_text="Reserved field for additional metadata (JSON format)"
    )

    class Meta:
        table_name = 'dataset'
        indexes = (
            # Index for quick lookup by project
            (('project_id',), False),
            # Composite index for biological queries
            (('organism', 'tissue_type'), False),
        )

    def __repr__(self) -> str:
        return (
            f"<Dataset id={self.id} "
            f"name='{self.dataset_name}' "
            f"project_id={self.project_id}>"
        )
