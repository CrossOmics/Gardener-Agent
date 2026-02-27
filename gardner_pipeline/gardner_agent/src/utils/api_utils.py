import httpx

from base.eums.model_type import APIKey
from loguru import logger


async def validate_api_key_online(key_type: str, secret: str) -> bool:
    """
    Validates the API key against the provider's API.
    Returns True if valid, False otherwise.
    """
    try:
        async with httpx.AsyncClient(timeout=10.0) as client:
            if key_type == APIKey.OPENAI_API_KEY.value:
                resp = await client.get(
                    "https://api.openai.com/v1/models",
                    headers={"Authorization": f"Bearer {secret}"}
                )
                return resp.status_code == 200

            elif key_type == APIKey.ANTHROPIC_API_KEY.value:
                resp = await client.get(
                    "https://api.anthropic.com/v1/models",
                    headers={"x-api-key": secret, "anthropic-version": "2023-06-01"}
                )
                return resp.status_code == 200

            elif key_type == APIKey.GOOGLE_API_KEY.value:
                resp = await client.get(
                    f"https://generativelanguage.googleapis.com/v1beta/models?key={secret}"
                )
                return resp.status_code == 200

            elif key_type == APIKey.HUGGINGFACE_TOKEN.value:
                resp = await client.get(
                    "https://huggingface.co/api/whoami-v2",
                    headers={"Authorization": f"Bearer {secret}"}
                )
                return resp.status_code == 200

            elif key_type == APIKey.MISTRAL_API_KEY.value:
                resp = await client.get(
                    "https://api.mistral.ai/v1/models",
                    headers={"Authorization": f"Bearer {secret}"}
                )
                return resp.status_code == 200

            elif key_type == APIKey.QWEN_API_KEY.value:
                # Qwen (DashScope) validation
                resp = await client.get(
                    "https://dashscope.aliyuncs.com/compatible-mode/v1/models",
                    headers={"Authorization": f"Bearer {secret}"}
                )
                return resp.status_code == 200

            elif key_type == APIKey.DEEPSEEK_API_KEY.value:
                # DeepSeek validation
                resp = await client.get(
                    "https://api.deepseek.com/models",
                    headers={"Authorization": f"Bearer {secret}"}
                )
                return resp.status_code == 200


    except Exception as e:
        logger.warning(f"Online validation failed for {key_type}: {e}")
        return False

    return False
