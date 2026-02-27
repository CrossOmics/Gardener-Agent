from datetime import datetime
from typing import Optional, Dict, Any
from peewee import IntegrityError

from infrastructure.database.model.user_preference import UserPreference


class UserPreferenceDAO:
    """
    Data Access Object for UserPreference.
    Handles CRUD operations for user global settings.
    """

    @staticmethod
    def get_by_name(name: str) -> Optional[UserPreference]:
        """
        Retrieve a UserPreference object by its unique preference_name.
        Excludes soft-deleted records.
        """
        return UserPreference.get_or_none(
            (UserPreference.preference_name == name) &
            (UserPreference.is_deleted == False)
        )

    @staticmethod
    def get_by_id(pk: int) -> Optional[UserPreference]:
        """
        Retrieve a preference by its primary key ID.
        """
        return UserPreference.get_or_none(UserPreference.id == pk)

    @staticmethod
    def update_timestamp(name: str) -> bool:
        """
        Directly update the 'updated_at' timestamp to now for a specific preference.
        Returns True if successful, False if not found.
        """
        now = datetime.utcnow()
        query = UserPreference.update(updated_at=now).where(
            (UserPreference.preference_name == name) &
            (UserPreference.is_deleted == False)
        )
        # execute() returns the number of rows modified
        return query.execute() > 0

    @staticmethod
    def create_preference(name: str, settings_json: Dict[str, Any]) -> Optional[UserPreference]:
        """
        Create a new preference.
        Sets 'created_at' and 'updated_at' to the exact same current timestamp.
        """
        now = datetime.utcnow()
        try:
            return UserPreference.create(
                preference_name=name,
                settings=settings_json,
                created_at=now,
                updated_at=now,
                is_deleted=False
            )
        except IntegrityError:
            # Handles Unique constraint violation (duplicate name)
            print(f"Error: Preference '{name}' already exists.")
            return None

    @staticmethod
    def delete_by_name(name: str) -> bool:
        """
        Delete a preference by name.
        Performs a 'Soft Delete' by setting is_deleted=True.
        """
        query = UserPreference.update(is_deleted=True).where(
            UserPreference.preference_name == name
        )
        return query.execute() > 0

    @staticmethod
    def get_latest_updated() -> Optional[UserPreference]:
        """
        Fetch the single most recently updated preference object.
        """
        return (UserPreference
                .select()
                .where(UserPreference.is_deleted == False)
                .order_by(UserPreference.updated_at.desc())
                .first())

    @staticmethod
    def count_all() -> int:
        """
        Get the total count of active (non-deleted) preferences.
        """
        return (UserPreference
                .select()
                .where(UserPreference.is_deleted == False)
                .count())

    @staticmethod
    def update_preference(name: str, settings_json: Dict[str, Any], new_name: Optional[str] = None) -> int:
        """
        Update settings and optionally rename.

        Logic:
        1. Find the target record by current name to get its ID.
        2. Perform the update using the ID (safer).
        3. Return the ID if successful.

        Returns:
            int: The ID of the updated row, or -1 if failure (not found or name conflict).
        """
        # Locate the record first to retrieve its ID
        target = UserPreference.get_or_none(
            (UserPreference.preference_name == name) &
            (UserPreference.is_deleted == False)
        )

        if not target:
            return -1

        target_id = target.id
        now = datetime.utcnow()

        # Prepare update data
        update_data = {
            UserPreference.settings: settings_json,
            UserPreference.updated_at: now
        }

        if new_name and new_name.strip():
            cleaned_new_name = new_name.strip()
            if cleaned_new_name != name:
                update_data[UserPreference.preference_name] = cleaned_new_name

        try:
            # Execute Update by ID
            query = UserPreference.update(update_data).where(
                UserPreference.id == target_id
            )

            rows_affected = query.execute()

            if rows_affected > 0:
                return target_id
            return -1

        except IntegrityError:
            print(f"[DAO Error] Rename failed. '{new_name}' might already exist.")
            return -1
