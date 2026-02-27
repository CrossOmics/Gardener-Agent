from base.base_call import call_backend


async def get_snapshot_info(snapshot_id: str):
    """
    Retrieve detailed information about a specific snapshot.

    Args:
        snapshot_id (str): The unique identifier for the snapshot.

    Returns:
        dict: A dictionary containing snapshot details including:
            - input: The parameters used for the experiment.
            - output: The biological metrics resulting from the experiment.
            - user_notes: Any notes added by the user.
    """
    # Request snapshot details from the backend
    # Endpoint: /snapshots/query/{snapshot_id}
    # This returns the 'input' (parameters) and 'output' (metrics)
    data = await call_backend("GET", f"/snapshots/query/{snapshot_id}")

    return data


async def get_snapshot_ancestors(snapshot_id: str):
    """
    Recursively traces the parent lineage of a given snapshot up to the root.

    Args:
        snapshot_id (str): The unique identifier for the snapshot.

    Returns:
        list: A list of ancestor nodes from the root down to the parent of the current snapshot.
    """
    # Request snapshot ancestors from the backend
    # Endpoint: /snapshots/lineage/ancestors/{snapshot_id}
    # This returns a list of SnapshotAncestorNode objects
    data = await call_backend("GET", f"/snapshots/lineage/ancestors/{snapshot_id}")

    return data


async def get_snapshot_ancestors_full(snapshot_id: str):
    """
    Recursively traces the parent lineage of a given snapshot up to the root, returning full details.

    Args:
        snapshot_id (str): The unique identifier for the snapshot.

    Returns:
        list: A list of full snapshot details from the root down to the parent of the current snapshot.
    """
    # Request full snapshot ancestors from the backend
    # Endpoint: /snapshots/lineage/ancestors/full/{snapshot_id}
    # This returns a list of SnapshotQueryResponse objects
    data = await call_backend("GET", f"/snapshots/lineage/ancestors/full/{snapshot_id}")

    return data


async def get_snapshot_peers(snapshot_id: str):
    """
    Retrieves all other snapshots that share the same dataset and branch as the given snapshot.

    Args:
        snapshot_id (str): The unique identifier for the reference snapshot.

    Returns:
        list: A list of full snapshot details for all peer snapshots in the same branch.
    """
    # Request peer snapshots from the backend
    # Endpoint: /snapshots/peers/{snapshot_id}
    # This returns a list of SnapshotQueryResponse objects
    data = await call_backend("GET", f"/snapshots/peers/{snapshot_id}")

    return data
