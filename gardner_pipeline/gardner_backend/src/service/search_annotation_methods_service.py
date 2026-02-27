from typing import List
from fastapi import Depends

from infrastructure.database.dao.annotation_method_dao import AnnotationMethodDAO
from dto.request.annotation_method_request import SearchAnnotationMethodRequest
from dto.response.annotation_method_response import AnnotationMethodResponse


class AnnotationMethodService:
    """
    Service layer for managing Annotation Methods.
    Handles logic for routing search requests based on method type.
    """

    def __init__(self, dao: AnnotationMethodDAO = Depends()):
        """
        Dependency injection for the Data Access Object.
        """
        self.dao = dao

    def search_methods(self, request: SearchAnnotationMethodRequest) -> List[AnnotationMethodResponse]:
        """
        Orchestrates the search logic.

        Logic:
        1. If request.type is 'all' (case-insensitive), perform a broad keyword search.
        2. Otherwise, perform a search filtered by the specific type (e.g., 'gseapy').
        3. Convert the resulting database models into Pydantic Response DTOs.
        """

        # Normalize type to lowercase for consistent comparison
        search_type = request.type.lower() if request.type else "all"
        keyword = request.keyword.strip() if request.keyword else ""

        # Routing Logic
        if search_type == "all":
            # Scenario 1: Search everything regardless of type
            results = self.dao.search_by_keyword(keyword)
        else:
            # Scenario 2: Search within a specific category (e.g., 'gseapy' or 'celltypist')
            results = self.dao.search_by_keyword_and_type(keyword, search_type)

        # Convert Peewee/ORM models to Pydantic Response DTOs
        return [AnnotationMethodResponse.model_validate(item) for item in results]
