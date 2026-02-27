from peewee import AutoField, CharField, TextField
from .base_model import BaseModel


class AnnotationMethod(BaseModel):
    """
    Model representing reference information for annotation methods.
    Stores metadata for GSEApy databases (Enrichr libraries) and CellTypist models.
    """

    # Primary Key (Auto-incrementing)
    id = AutoField(primary_key=True)

    # The unique identifier/name for the GSEApy database or CellTypist model.
    # e.g., "CellMarker_2024" or "Immune_All_Low.pkl"
    method_name = CharField(unique=True, null=False, index=True, help_text="Unique name of the method")

    # The category of the method.
    # Enumeration: ["gseapy", "celltypist"]
    type = CharField(null=False, index=True, help_text="Type: 'gseapy' or 'celltypist'")

    # Detailed description of the annotation method/library.
    description = TextField(null=True, help_text="Specific description of the method")

    # The species associated with this reference (e.g., "Human", "Mouse").
    species = CharField(null=False, help_text="Target species")

    # The specific organ or tissue type associated with this model (e.g., "Lung", "Blood", "Pan-tissue").
    # Can be null if it is a general-purpose model.
    organ = CharField(null=True, help_text="Associated organ or tissue")

    class Meta:
        table_name = 'annotation_methods'

    def __str__(self):
        return f"[{self.type}] {self.method_name} ({self.species})"
