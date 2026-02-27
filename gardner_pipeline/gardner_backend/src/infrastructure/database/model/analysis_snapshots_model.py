from peewee import (
    AutoField,
    CharField,
    DateTimeField,
    TextField,
    ForeignKeyField
)
from datetime import datetime
from playhouse.sqlite_ext import JSONField
from .base_model import BaseModel
from .dataset_model import Dataset


class AnalysisSnapshot(BaseModel):
    """
    Analysis Snapshot table for tracking analysis history and branching.
    """

    # Primary Key
    id = AutoField(
        primary_key=True,
        help_text="Auto-increment primary key"
    )

    # Business Logic ID (System Generated)
    snapshot_id = CharField(
        max_length=128,
        unique=True,
        index=True,
        null=False,
        help_text="Unique system identifier for the snapshot"
    )

    snapshot_name = CharField(
        max_length=255,
        unique=False,
        null=False,
        default="New Snapshot",
        help_text="User-defined display name for the experiment"
    )

    # parent snapshot id
    parent_snapshot_id = CharField(
        max_length=128,
        unique=False,
        null=True,
        help_text="Parent process's snapshot id"
    )

    # Foreign Keys
    dataset_id = ForeignKeyField(
        Dataset,
        field='dataset_id',
        column_name='dataset_id',
        backref='snapshots',
        on_delete='CASCADE',
        help_text="Reference to the raw dataset business id"
    )

    # Metadata
    branch_name = CharField(
        max_length=255,
        null=False,
        help_text="Technical name of this analysis branch (e.g., 'Leiden Res 0.5')"
    )

    snapshot_path = TextField(
        null=False,
        help_text="Snapshot's relative path key within workspace"
    )

    # Scanpy Technical parameters (JSON Field)
    params_json = JSONField(
        null=True,
        help_text="Input parameters (QC thresholds, resolution, etc.) in JSON format"
    )

    thumbnail_json = JSONField(
        null=True,
        help_text="Map of visualization names and analysis results to relative file paths in JSON format"
    )

    # Temporal Fields
    create_time = DateTimeField(
        default=datetime.utcnow,
        null=False,
        help_text="Experiment start timestamp"
    )

    end_time = DateTimeField(
        null=True,
        help_text="Experiment completion timestamp"
    )

    user_notes = TextField(
        null=True,
        help_text="User's evaluation or notes for this result"
    )

    ext_info = JSONField(
        null=True,
        help_text="Reserved field for additional metadata"
    )

    class Meta:
        table_name = 'analysis_snapshots'
        indexes = (
            (('dataset_id', 'create_time'), False),
        )

    def __repr__(self) -> str:
        return (
            f"<Snapshot id={self.id} "
            f"id='{self.snapshot_id}' "
            f"display_name='{self.snapshot_name}'>"
        )