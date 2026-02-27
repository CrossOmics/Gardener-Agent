from peewee import (
    AutoField,
    CharField,
    DateTimeField,
    BooleanField
)
from datetime import datetime
from .base_model import BaseModel

from playhouse.sqlite_ext import JSONField


class UserPreference(BaseModel):
    """
    Model representing user global preferences.
    Stores custom named configurations (pipeline parameters) for quick reuse.
    """

    # Primary Key (Auto-incrementing)
    id = AutoField(primary_key=True)

    # Custom name for this preference setting (e.g., "Default High Sensitivity")
    preference_name = CharField(
        unique=True,
        null=False,
        index=True,
        help_text="Global custom name for the preference setting"
    )

    settings = JSONField(
        null=False,
        help_text="Pipeline parameters stored as a JSON object"
    )

    # Timestamp when this preference was created
    created_at = DateTimeField(
        default=datetime.utcnow,
        null=False,
        help_text="Creation timestamp"
    )

    # Timestamp when this preference was last updated
    updated_at = DateTimeField(
        default=datetime.utcnow,
        null=False,
        help_text="Last update timestamp"
    )

    # Soft deletion flag (0/False = Active, 1/True = Deleted)
    is_deleted = BooleanField(
        default=False,
        null=False,
        index=True,
        help_text="Soft delete flag"
    )

    class Meta:
        table_name = 'user_preference'

    def __str__(self):
        return f"[Preference] {self.preference_name} (Created: {self.created_at})"
