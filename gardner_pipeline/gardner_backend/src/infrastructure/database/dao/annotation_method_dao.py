from typing import List, Optional

from infrastructure.database.model.annotation_method import AnnotationMethod


class AnnotationMethodDAO:
    """
    Data Access Object for managing AnnotationMethod database operations.
    """

    @staticmethod
    def search_by_keyword(keyword: str) -> List[AnnotationMethod]:
        """
        Searches for annotation methods where the keyword appears in the name,
        description, species, or organ.

        Args:
            keyword (str): The search term.

        Returns:
            List[AnnotationMethod]: A list of matching annotation methods.
        """
        if not keyword:
            return list(AnnotationMethod.select().order_by(AnnotationMethod.method_name))

        search_pattern = f"%{keyword}%"

        # Using '**' for case-insensitive LIKE in Peewee (mapped to ILIKE or LIKE based on DB)
        query = (AnnotationMethod
                 .select()
                 .where(
            (AnnotationMethod.method_name ** search_pattern) |
            (AnnotationMethod.description ** search_pattern) |
            (AnnotationMethod.species ** search_pattern) |
            (AnnotationMethod.organ ** search_pattern)
        )
                 .order_by(AnnotationMethod.method_name))

        return list(query)

    @staticmethod
    def search_by_keyword_and_type(keyword: str, method_type: str) -> List[AnnotationMethod]:
        """
        Searches for annotation methods matching a keyword, restricted to a specific type
        (e.g., only 'gseapy' databases or only 'celltypist' models).

        Args:
            keyword (str): The search term.
            method_type (str): The type to filter by (e.g., "gseapy", "celltypist").

        Returns:
            List[AnnotationMethod]: A list of matching annotation methods of the specific type.
        """
        if not keyword:
            return list(AnnotationMethod.select().where(AnnotationMethod.type == method_type))

        search_pattern = f"%{keyword}%"

        query = (AnnotationMethod
                 .select()
                 .where(
            (AnnotationMethod.type == method_type) &
            (
                    (AnnotationMethod.method_name ** search_pattern) |
                    (AnnotationMethod.description ** search_pattern) |
                    (AnnotationMethod.species ** search_pattern) |
                    (AnnotationMethod.organ ** search_pattern)
            )
        )
                 .order_by(AnnotationMethod.method_name))

        return list(query)

    @staticmethod
    def get_all_by_type(method_type: str) -> List[AnnotationMethod]:
        """
        Retrieves all annotation methods for a specific type without keyword filtering.
        Useful for populating dropdown lists.

        Args:
            method_type (str): The type to filter by.

        Returns:
            List[AnnotationMethod]: All records matching the type.
        """
        query = (AnnotationMethod
                 .select()
                 .where(AnnotationMethod.type == method_type)
                 .order_by(AnnotationMethod.method_name))

        return list(query)

    @staticmethod
    def get_by_name(name: str) -> Optional[AnnotationMethod]:
        """
        Retrieves a single annotation method by its unique name.

        Args:
            name (str): The unique method_name.

        Returns:
            Optional[AnnotationMethod]: The object if found, else None.
        """
        return AnnotationMethod.get_or_none(AnnotationMethod.method_name == name)
