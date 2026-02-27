import random
import string
from datetime import datetime


def generate_business_id(prefix: str, random_length: int = 5) -> str:
    """
    Generate a unique business ID.
    Format: <prefix>_<unix_timestamp>_<random_digits>
    """
    timestamp = int(datetime.utcnow().timestamp())  # Unix seconds
    random_digits = ''.join(random.choices(string.digits, k=random_length))
    return f"{prefix}_{timestamp}_{random_digits}"


def generate_filename(prefix: str, suffix: str) -> str:
    """
    Generate a filename using business ID and file suffix.
    Format: <business_id>.<suffix>
    """
    business_id = generate_business_id(prefix)
    return f"{business_id}.{suffix}"


def generate_session_id_from_project(project_id: str, prefix: str) -> str:
    """
    Generate a session ID derived from project ID by replacing the first
    underscore-separated segment with the given prefix.

    Example:
        project_id: p_20260126150926_21762
        session_id: session_20260126150926_21762_48392
    """
    parts = project_id.split("_")

    # Replace the first segment (e.g., "p") with prefix
    parts[0] = prefix

    base = "_".join(parts)
    rand_suffix = ''.join(random.choices(string.digits, k=3))

    return f"{base}_{rand_suffix}"
