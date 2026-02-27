import os
import shelve
from typing import Tuple
from loguru import logger

# --- LangChain Imports ---
from langchain_openai import ChatOpenAI
from langchain_google_genai import ChatGoogleGenerativeAI
from langchain_anthropic import ChatAnthropic
from langchain_mistralai import ChatMistralAI
from langchain_huggingface import ChatHuggingFace, HuggingFaceEndpoint
from langchain_core.messages import AIMessage

# --- Internal Imports ---
from base.base_api_manager import CredentialManager
from base.eums.model_type import APIKey
from base.constants.local_config import CONFIG_DB, CUR_LLM


class MockLLM:
    """
    Fallback LLM used when no valid API key is configured.
    Returns a standard error message instructing the user to configure the agent.
    """

    def __init__(self):
        pass

    def invoke(self, *args, **kwargs):
        return AIMessage(
            content="**System Alert**: No valid API Key configured.\n\nPlease go to the **Model Config (Top Right)** in the application to set your API Key (e.g., OpenAI, Gemini, Anthropic) to enable the agent.")

    async def ainvoke(self, *args, **kwargs):
        return AIMessage(
            content="**System Alert**: No valid API Key configured.\n\nPlease go to the **Model Config (Top Right)** in the application to set your API Key (e.g., OpenAI, Gemini, Anthropic) to enable the agent.")

    def bind_tools(self, *args, **kwargs):
        # Return self to prevent attribute errors when the graph tries to bind tools
        return self


def _create_openai_llms(secret: str) -> Tuple[ChatOpenAI, ChatOpenAI]:
    base = ChatOpenAI(api_key=secret, model="gpt-4o-mini", temperature=0)
    strong = ChatOpenAI(api_key=secret, model="gpt-4o", temperature=0)
    return base, strong


def _create_google_llms(secret: str) -> Tuple[ChatGoogleGenerativeAI, ChatGoogleGenerativeAI]:
    base = ChatGoogleGenerativeAI(google_api_key=secret, model="gemini-3-flash-preview", temperature=0)
    strong = ChatGoogleGenerativeAI(google_api_key=secret, model="gemini-3-pro-preview", temperature=0)
    return base, strong


def _create_anthropic_llms(secret: str) -> Tuple[ChatAnthropic, ChatAnthropic]:
    base = ChatAnthropic(api_key=secret, model="claude-haiku-4-5", temperature=0)
    strong = ChatAnthropic(api_key=secret, model="claude-opus-4-6", temperature=0)
    return base, strong


def _create_mistral_llms(secret: str) -> Tuple[ChatMistralAI, ChatMistralAI]:
    base = ChatMistralAI(api_key=secret, model="mistral-small-latest", temperature=0)
    strong = ChatMistralAI(api_key=secret, model="mistral-large-latest", temperature=0)
    return base, strong


def _create_huggingface_llms(secret: str) -> Tuple[ChatHuggingFace, ChatHuggingFace]:
    # Base: Llama-3-8B
    llm_base = HuggingFaceEndpoint(
        repo_id="meta-llama/Meta-Llama-3-8B-Instruct",
        task="text-generation",
        max_new_tokens=512,
        do_sample=False,
        repetition_penalty=1.03,
        huggingfacehub_api_token=secret
    )
    base = ChatHuggingFace(llm=llm_base)

    # Strong: Llama-3-70B (if available/supported by endpoint) or same as base
    llm_strong = HuggingFaceEndpoint(
        repo_id="meta-llama/Meta-Llama-3-70B-Instruct",
        task="text-generation",
        max_new_tokens=1024,
        do_sample=False,
        repetition_penalty=1.03,
        huggingfacehub_api_token=secret
    )
    strong = ChatHuggingFace(llm=llm_strong)
    return base, strong


def _create_qwen_llms(secret: str) -> Tuple[ChatOpenAI, ChatOpenAI]:
    base = ChatOpenAI(
        api_key=secret,
        base_url="https://dashscope.aliyuncs.com/compatible-mode/v1",
        model="qwen-plus",
        temperature=0
    )
    strong = ChatOpenAI(
        api_key=secret,
        base_url="https://dashscope.aliyuncs.com/compatible-mode/v1",
        model="qwen-max",
        temperature=0
    )
    return base, strong


def _create_deepseek_llms(secret: str) -> Tuple[ChatOpenAI, ChatOpenAI]:
    base = ChatOpenAI(
        api_key=secret,
        base_url="https://api.deepseek.com",
        model="deepseek-chat",
        temperature=0
    )
    strong = ChatOpenAI(
        api_key=secret,
        base_url="https://api.deepseek.com",
        model="deepseek-reasoner",
        temperature=0
    )
    return base, strong


def init_llms() -> Tuple[object, object]:
    """
    Initializes the Base and Strong LLMs based on the following priority:
    1. Local Config (Shelve) + Keyring (CredentialManager)
    2. Environment Variables (Fallback)
    3. MockLLM (Safe Fallback)

    Returns:
        Tuple[base_llm, strong_llm]
    """
    logger.info("Initializing LLMs (Base & Strong)...")

    # 1. Try loading from Local Config (Shelve + Keyring)
    try:
        with shelve.open(CONFIG_DB) as db:
            key_type = db.get(CUR_LLM)

        if key_type:
            logger.info(f"Found configured LLM type: {key_type}")
            manager = CredentialManager(app_name="gardener_agent")

            if key_type in APIKey.__members__:
                secret = manager.get_secret(APIKey[key_type].value)

                if secret:
                    logger.info("Successfully retrieved API key from keyring.")

                    if key_type == "OPENAI_API_KEY":
                        return _create_openai_llms(secret)
                    elif key_type == "GOOGLE_API_KEY":
                        return _create_google_llms(secret)
                    elif key_type == "ANTHROPIC_API_KEY":
                        return _create_anthropic_llms(secret)
                    elif key_type == "MISTRAL_API_KEY":
                        return _create_mistral_llms(secret)
                    elif key_type == "HUGGINGFACE_TOKEN":
                        return _create_huggingface_llms(secret)
                    elif key_type == "QWEN_API_KEY":
                        return _create_qwen_llms(secret)
                    elif key_type == "DEEPSEEK_API_KEY":
                        return _create_deepseek_llms(secret)
                else:
                    logger.warning(f"Key type {key_type} is set, but no secret found in keyring.")
            else:
                logger.warning(f"Unknown key type in config: {key_type}")

    except Exception as e:
        logger.error(f"Failed to load LLM from local config: {e}")

    # 2. Fallback: Try Environment Variables
    logger.info("Attempting fallback to Environment Variables...")
    try:
        if os.getenv("GOOGLE_API_KEY"):
            return _create_google_llms(os.getenv("GOOGLE_API_KEY"))
        if os.getenv("OPENAI_API_KEY"):
            return _create_openai_llms(os.getenv("OPENAI_API_KEY"))
        if os.getenv("ANTHROPIC_API_KEY"):
            return _create_anthropic_llms(os.getenv("ANTHROPIC_API_KEY"))

    except Exception as e:
        logger.error(f"Failed to load LLM from environment: {e}")

    # 3. Final Fallback: MockLLM
    logger.warning("No valid LLM configuration found. Using MockLLM.")
    return MockLLM(), MockLLM()


def reload_llm():
    """
    Reloads the global base_llm and strong_llm instances.
    Call this function after updating the API key configuration to apply changes immediately.
    """
    global base_llm, strong_llm
    logger.info("Reloading LLM configurations...")
    base_llm, strong_llm = init_llms()
    return base_llm, strong_llm


# --- Global Instances ---
# These instances are imported by other modules (executor, planner, etc.)
base_llm, strong_llm = init_llms()
