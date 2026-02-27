from enum import Enum

class APIKey(str, Enum):
    """
    Enumeration of common API key identifiers.
    These string constants are used as keys when storing/retrieving secrets via CredentialManager.
    """
    
    # OpenAI API Key for accessing GPT models
    OPENAI_API_KEY = "OPENAI_API_KEY"

    # Anthropic API Key for accessing Claude models
    ANTHROPIC_API_KEY = "ANTHROPIC_API_KEY"

    # Google API Key for accessing Gemini or PaLM models
    GOOGLE_API_KEY = "GOOGLE_API_KEY"

    # Hugging Face Token for accessing models on the Hub
    HUGGINGFACE_TOKEN = "HUGGINGFACE_TOKEN"

    # Mistral AI API Key
    MISTRAL_API_KEY = "MISTRAL_API_KEY"

    # Qwen API Key (DashScope)
    QWEN_API_KEY = "QWEN_API_KEY"

    # DeepSeek API Key
    DEEPSEEK_API_KEY = "DEEPSEEK_API_KEY"
