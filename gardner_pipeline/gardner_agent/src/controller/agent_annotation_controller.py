import json
from fastapi import APIRouter, HTTPException
from langchain_core.prompts import PromptTemplate
from langchain_core.output_parsers import JsonOutputParser
from loguru import logger

import base.base_model_manager as base_model_manager
from base.dtos.request.annotation_request import AnnotationRequest
from prompt.annotation_prompt import ANNOTATION_PROMPT

router = APIRouter(tags=["Agent Annotation"])


@router.post("/api/annotation")
async def annotate_clusters(request: AnnotationRequest):
    """
    Refines uncertain cell type annotations using LLM reasoning.
    Takes a dictionary of uncertain clusters (DGEs + candidate labels) and returns
    a JSON object with predicted cell types and confidence scores.
    """
    try:
        uncertain_data = request.uncertain_data
        valid_cell_types = request.valid_cell_types

        if not uncertain_data:
            return {"llm_annotation": {}}

        logger.info(f"Received annotation request for {len(uncertain_data)} clusters.")

        # 1. Prepare Prompt
        prompt_template = PromptTemplate(
            template=ANNOTATION_PROMPT,
            input_variables=["uncertain_data", "valid_cell_types"]
        )

        # 2. Convert input data to JSON string for the prompt
        input_json_str = json.dumps(uncertain_data, indent=2)

        # 3. Format the prompt and invoke LLM
        formatted_prompt = prompt_template.format_prompt(
            uncertain_data=input_json_str,
            valid_cell_types=valid_cell_types if valid_cell_types else "Not provided"
        )
        messages = formatted_prompt.to_messages()

        logger.info("Invoking LLM for annotation...")

        # Call the LLM directly with the messages
        llm_response = await base_model_manager.strong_llm.ainvoke(messages)

        # Extract token usage
        usage = llm_response.usage_metadata
        input_tokens = usage.get("input_tokens", 0)
        output_tokens = usage.get("output_tokens", 0)

        logger.debug(f"[LLM Token Usage] Input: {input_tokens} | Output: {output_tokens}")

        # Parse the output using JsonOutputParser
        parser = JsonOutputParser()
        result = parser.invoke(llm_response)

        logger.success("LLM annotation completed successfully.")

        # 4. Return Result
        # The backend expects a key "llm_annotation" containing the result dict.
        return {"llm_annotation": result}

    except Exception as e:
        logger.error(f"Error during LLM annotation: {e}")
        raise HTTPException(status_code=500, detail=f"LLM Annotation failed: {str(e)}")
