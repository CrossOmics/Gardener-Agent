ANNOTATION_PROMPT = """
You are an expert computational biologist specializing in single-cell RNA sequencing (scRNA-seq) analysis.
Your task is to resolve uncertain cell type annotations for specific clusters based on their top Differentially Expressed Genes (DGEs) and a set of candidate labels provided by other tools (e.g., CellTypist, GSEApy).

### INPUT DATA:
You will receive a JSON object where keys are Cluster IDs and values contain:
1. "dge": A list of top 20 marker genes for that cluster.
2. "uncertain_labels": A list of candidate cell types suggested by other methods (which had low confidence).

### YOUR GOAL:
For each cluster, analyze the "dge" list to determine the most likely biological cell type.
- Use the "uncertain_labels" as hints but do not be bound by them if the marker genes strongly suggest otherwise.
- The predicted cell type MUST be one of the names from the "CELL ONTOLOGY" list provided below.
- Assign a "predicted_cell_type" that is specific and biologically accurate (e.g., "CD8+ T Cell", "Goblet Cell", "Fibroblast").
- Assign an "average_confidence" score between 0.0 and 1.0 based on how well the marker genes match the cell type signature.

### CELL ONTOLOGY:
Cell type names must come from this list, one per line:
{valid_cell_types}

### OUTPUT FORMAT:
Return a STRICT JSON object. The keys must be the original Cluster IDs.
The values must be objects with:
- "predicted_cell_type": (string) The final annotated cell type.
- "average_confidence": (float) Your confidence score.

### EXAMPLES:

#### Input:
{{
  "9": {{
    "dge": ["Muc2", "Tff3", "Fcgbp", "Agr2", "Zg16", "Clca1"],
    "uncertain_labels": ["Epithelial", "Secretory"]
  }},
  "12": {{
    "dge": ["Cd3d", "Cd3e", "Cd8a", "Ccl5", "Gzmb", "Nkg7"],
    "uncertain_labels": ["T cell", "NK cell"]
  }}
}}

#### Output:
{{
  "9": {{
    "predicted_cell_type": "Goblet Cell",
    "average_confidence": 0.95
  }},
  "12": {{
    "predicted_cell_type": "CD8+ Cytotoxic T Cell",
    "average_confidence": 0.92
  }}
}}

### IMPORTANT CONSTRAINTS:
1. Return ONLY the JSON object. Do not include markdown formatting (like ```json ... ```).
2. Do not include explanations or extra text.
3. Ensure the JSON is valid and parsable.
4. The predicted cell type MUST be one of the names from the "CELL ONTOLOGY" list.

### CURRENT TASK:
Input Data:
{uncertain_data}
"""
