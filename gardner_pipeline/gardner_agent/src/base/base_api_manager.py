import keyring
import logging
from typing import Optional


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class CredentialManager:
    """
    A wrapper class to manage secure credentials using the system keyring.
    """

    def __init__(self, app_name: str):
        """
        Initialize with the application name to namespace the credentials.
        """
        self.app_name = app_name

    def save_secret(self, service_name: str, secret_value: str) -> None:
        """
        Securely saves a password or API key.
        """
        try:
            keyring.set_password(self.app_name, service_name, secret_value)
            logger.info(f"Successfully saved secret for: {service_name}")
        except Exception as e:
            logger.error(f"Failed to save secret for {service_name}: {e}")

    def get_secret(self, service_name: str) -> Optional[str]:
        """
        Retrieves a password or API key. Returns None if not found.
        """
        try:
            return keyring.get_password(self.app_name, service_name)
        except Exception as e:
            logger.error(f"Failed to retrieve secret for {service_name}: {e}")
            return None

    def delete_secret(self, service_name: str) -> None:
        """
        Removes a password or API key from the system.
        """
        try:
            keyring.delete_password(self.app_name, service_name)
            logger.info(f"Successfully deleted secret for: {service_name}")
        except keyring.errors.PasswordDeleteError:
            logger.warning(f"Secret for {service_name} not found, nothing to delete.")
        except Exception as e:
            logger.error(f"Error deleting secret for {service_name}: {e}")


if __name__ == "__main__":
    # 1. Initialize the manager
    manager = CredentialManager(app_name="gardener_agent")

    KEY_ID = "openai_api_key"
    FAKE_KEY = "sk-xxxxx"

    # 2. Save the key
    manager.save_secret(KEY_ID, FAKE_KEY)

    # 3. Retrieve and verify
    stored_key = manager.get_secret(KEY_ID)
    print(f"Retrieved Key: {stored_key}")

    # 4. Clean up
    manager.delete_secret(KEY_ID)