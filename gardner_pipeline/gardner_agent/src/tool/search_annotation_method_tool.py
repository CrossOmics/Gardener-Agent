from typing import Literal
from langchain_core.tools import tool
from base.base_call import call_backend


@tool(
    description="Search for available reference methods (GSEApy libraries or CellTypist models) for cell annotation.",
    parse_docstring=True
)
async def search_annotation_methods(
        keyword: str = "",
        method_type: Literal['gseapy', 'celltypist', 'all'] = "all"
) -> str:
    """
    Query the system's database for valid annotation references.

    # CORE RESPONSIBILITY
    Use this tool when you need to find valid parameter values for `run_full_annotation_pipeline`.
    - If the user asks "What CellTypist models do you have?", search with type='celltypist'.
    - If the user asks "Can we use the Human Lung Atlas?", search with keyword='Lung'.

    # STRATEGY
    1. **Don't Guess**: Exact spelling matters for annotation models. Use this tool to verify names like 'Immune_All_Low.pkl' or 'CellMarker_2024'.
    2. **Filter by Type**: Use 'gseapy' for pathway/enrichment libraries and 'celltypist' for automated classification models.

    Args:
        keyword: Search term for name, description, species, or organ (e.g., "Lung", "Immune").
        method_type: Filter by method type. Options: 'gseapy', 'celltypist', 'all' (default).
    """
    try:
        # The controller expects 'type' in the request, but 'type' is a reserved keyword in Python,
        # so we use 'method_type' in the function signature and map it.
        params = {
            "keyword": keyword,
            "type": method_type
        }

        # Call backend GET endpoint
        results = await call_backend("GET", "/methods/annotation/search", params)

        if not results:
            return f"No annotation methods found matching keyword='{keyword}' and type='{method_type}'."

        # Format the list of results into a readable string
        formatted_results = []
        for item in results:
            # item matches AnnotationMethodResponse DTO
            name = item.get("method_name")
            m_type = item.get("type")
            desc = item.get("description") or "No description"
            species = item.get("species")

            formatted_results.append(f"- [{m_type.upper()}] {name} ({species}): {desc}")

        return (
                f"Found {len(results)} matching annotation methods:\n" +
                "\n".join(formatted_results) +
                "\n\nObservation: You can use these exact names in the annotation pipeline."
        )

    except Exception as e:
        return f"Search for annotation methods failed: {str(e)}"
