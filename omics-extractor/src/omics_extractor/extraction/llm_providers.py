"""LLM provider backends for metadata extraction.

This module provides a pluggable interface for different LLM providers,
allowing easy switching between:
- Claude API (Anthropic)
- Gemini via Google Vertex AI
- Local models via VLLM (Llama, Mixtral, Qwen, etc.)
- Local models via Transformers
- Other APIs (OpenAI, etc.)

This enables benchmarking different models to find the best for your use case.
"""

from abc import ABC, abstractmethod
from typing import Dict, Any, Optional, List
import json
import os
from .llm_schemas import LLMExtractionResult
from ..enrich.providers.server_vllm import ServerVLLMProvider
import logging

logger = logging.getLogger(__name__)

# Unified default Claude model (override via CLAUDE_MODEL env)
DEFAULT_CLAUDE_MODEL = os.getenv("CLAUDE_MODEL", "claude-3-5-sonnet-20241022")


class LLMProvider(ABC):
    """Abstract base class for LLM providers."""

    @abstractmethod
    def extract(
        self,
        prompt: str,
        temperature: float = 0.0,
        max_tokens: int = 1000,
    ) -> LLMExtractionResult:
        """
        Extract metadata using the LLM.

        Args:
            prompt: Full extraction prompt
            temperature: Sampling temperature (0 = deterministic)
            max_tokens: Maximum tokens to generate

        Returns:
            LLMExtractionResult with extracted fields
        """
        pass

    @abstractmethod
    def get_model_name(self) -> str:
        """Get the model identifier."""
        pass

    @abstractmethod
    def supports_batching(self) -> bool:
        """Whether this provider supports batched inference."""
        pass

    def extract_batch(
        self,
        prompts: List[str],
        temperature: float = 0.0,
        max_tokens: int = 1000,
    ) -> List[LLMExtractionResult]:
        """
        Extract metadata for multiple prompts (batched).

        Default implementation calls extract() sequentially.
        Override for efficient batching.

        Args:
            prompts: List of extraction prompts
            temperature: Sampling temperature
            max_tokens: Maximum tokens per generation

        Returns:
            List of LLMExtractionResults
        """
        return [self.extract(p, temperature, max_tokens) for p in prompts]


class ClaudeProvider(LLMProvider):
    """Anthropic Claude API provider."""

    def __init__(
        self,
        api_key: Optional[str] = None,
        model: str = DEFAULT_CLAUDE_MODEL,
    ):
        """
        Initialize Claude provider.

        Args:
            api_key: Anthropic API key (or from ANTHROPIC_API_KEY env)
            model: Claude model to use (default: DEFAULT_CLAUDE_MODEL or CLAUDE_MODEL env)
        """
        self.api_key = api_key or os.environ.get("ANTHROPIC_API_KEY")
        if not self.api_key:
            raise ValueError("Anthropic API key required")

        self.model = model

        # Lazy import to avoid requiring anthropic for local models
        import anthropic
        self.client = anthropic.Anthropic(api_key=self.api_key)

    def extract(
        self,
        prompt: str,
        temperature: float = 0.0,
        max_tokens: int = 1000,
    ) -> LLMExtractionResult:
        """Extract using Claude API."""
        import time

        start = time.time()

        message = self.client.messages.create(
            model=self.model,
            max_tokens=max_tokens,
            temperature=temperature,
            messages=[{"role": "user", "content": prompt}],
        )

        latency_ms = (time.time() - start) * 1000

        # Parse JSON response
        response_text = message.content[0].text

        try:
            data = json.loads(response_text)
        except json.JSONDecodeError:
            # Try to extract JSON from markdown code blocks
            if "```json" in response_text:
                json_str = response_text.split("```json")[1].split("```")[0].strip()
                data = json.loads(json_str)
            else:
                raise ValueError(f"Failed to parse JSON from Claude response: {response_text}")

        # Build result
        result = LLMExtractionResult(
            organism=data.get("organism"),
            tissue=data.get("tissue"),
            cell_type=data.get("cell_type"),
            cell_line=data.get("cell_line"),
            treatment=data.get("treatment"),
            genotype=data.get("genotype"),
            strain=data.get("strain"),
            age=data.get("age"),
            sex=data.get("sex"),
            developmental_stage=data.get("developmental_stage"),
            confidence=data.get("confidence", {}),
            reasoning=data.get("reasoning"),
            model_name=self.model,
            tokens_used=message.usage.input_tokens + message.usage.output_tokens,
            latency_ms=latency_ms,
        )

        return result

    def get_model_name(self) -> str:
        return self.model

    def supports_batching(self) -> bool:
        return False  # Claude API doesn't support true batching


class GeminiProvider(LLMProvider):
    """Google Gemini via Vertex AI provider."""

    def __init__(
        self,
        project: Optional[str] = None,
        location: str = "us-central1",
        model: str = "gemini-2.5-flash",
    ):
        """
        Initialize Gemini provider via Vertex AI.

        Args:
            project: Google Cloud project ID (or from GOOGLE_CLOUD_PROJECT env)
            location: GCP region (default: us-central1)
            model: Gemini model name (default: gemini-2.5-flash)
        """
        self.project = project or os.environ.get("GOOGLE_CLOUD_PROJECT")
        if not self.project:
            raise ValueError(
                "Google Cloud project required. Set GOOGLE_CLOUD_PROJECT or pass --gemini-project"
            )
        self.location = location or os.environ.get("GOOGLE_CLOUD_LOCATION", "us-central1")
        self.model = model

        try:
            from google import genai
            from google.genai.types import HttpOptions
        except ImportError:
            raise ImportError(
                "google-genai not installed. Install with: pip install google-genai"
            )

        self.client = genai.Client(
            vertexai=True,
            project=self.project,
            location=self.location,
            http_options=HttpOptions(api_version="v1"),
        )

    def extract(
        self,
        prompt: str,
        temperature: float = 0.0,
        max_tokens: int = 1000,
    ) -> LLMExtractionResult:
        """Extract using Gemini via Vertex AI."""
        import time

        start = time.time()

        response = self.client.models.generate_content(
            model=self.model,
            contents=prompt,
        )
        logger.debug("Gemini raw response text length=%s", len(getattr(response, 'text', '') or ''))

        latency_ms = (time.time() - start) * 1000
        response_text = response.text

        try:
            data = json.loads(response_text)
        except json.JSONDecodeError:
            if "```json" in response_text:
                json_str = response_text.split("```json")[1].split("```")[0].strip()
                data = json.loads(json_str)
            else:
                raise ValueError(f"Failed to parse JSON from Gemini response: {response_text}")

        tokens_used = None
        try:
            tokens_used = response.usage_metadata.total_token_count
        except Exception:
            pass

        return LLMExtractionResult(
            organism=data.get("organism"),
            tissue=data.get("tissue"),
            cell_type=data.get("cell_type"),
            cell_line=data.get("cell_line"),
            treatment=data.get("treatment"),
            genotype=data.get("genotype"),
            strain=data.get("strain"),
            age=data.get("age"),
            sex=data.get("sex"),
            developmental_stage=data.get("developmental_stage"),
            confidence=data.get("confidence", {}),
            reasoning=data.get("reasoning"),
            model_name=self.model,
            tokens_used=tokens_used,
            latency_ms=latency_ms,
        )

    def get_model_name(self) -> str:
        return self.model

    def supports_batching(self) -> bool:
        return False


class TransformersProvider(LLMProvider):
    """Local Transformers provider with env-gated trust_remote_code."""

    def __init__(
        self,
        model_path: str = "meta-llama/Llama-3.1-70B-Instruct",
        tensor_parallel_size: int = 1,
        gpu_memory_utilization: float = 0.9,
        trust_remote_code: Optional[bool] = None,
    ):
        """
        Initialize VLLM provider.

        Args:
            model_path: HuggingFace model path or local path
            tensor_parallel_size: Number of GPUs for tensor parallelism (auto-adjusted if needed)
            gpu_memory_utilization: Fraction of GPU memory to use
        """
        self.model_path = model_path
        # Resolve device lazily; default to CUDA if available and requested
        requested_device = kwargs.get('device', 'cuda') if 'device' in kwargs else 'cuda'
        try:
            import torch
            self.device = 'cuda' if (requested_device == 'cuda' and torch.cuda.is_available()) else 'cpu'
        except Exception:
            self.device = 'cpu'

        try:
            from transformers import AutoTokenizer, AutoModelForCausalLM
            import torch
        except ImportError:
            raise ImportError("transformers not installed. Install with: pip install transformers torch")

        # Gate trust_remote_code via env or arg
        def _as_bool(x):
            if x is None:
                return None
            if isinstance(x, bool):
                return x
            s = str(x).strip().lower()
            return s in ('1','true','yes','y','on')
        trust = _as_bool(trust_remote_code)
        if trust is None:
            trust = _as_bool(os.getenv('VLLM_TRUST_REMOTE_CODE', '1'))

        logger.info("Loading %s with Transformers (trust_remote_code=%s)...", model_path, bool(trust))
        self.tokenizer = AutoTokenizer.from_pretrained(model_path, trust_remote_code=bool(trust))
        self.model = AutoModelForCausalLM.from_pretrained(
            model_path,
            trust_remote_code=bool(trust),
            device_map="auto" if self.device == "cuda" else None,
            torch_dtype=torch.float16 if self.device == "cuda" else torch.float32,
        )
        logger.info("Model loaded: %s", self.model_path)

    def extract(
        self,
        prompt: str,
        temperature: float = 0.0,
        max_tokens: int = 1000,
    ) -> LLMExtractionResult:
        """Extract using Transformers."""
        import time
        import torch

        start = time.time()

        # Tokenize
        inputs = self.tokenizer(prompt, return_tensors="pt").to(self.device)

        # Generate
        with torch.no_grad():
            outputs = self.model.generate(
                **inputs,
                max_new_tokens=max_tokens,
                temperature=temperature if temperature > 0 else None,
                do_sample=temperature > 0,
                pad_token_id=self.tokenizer.eos_token_id,
            )

        # Decode
        response_text = self.tokenizer.decode(outputs[0][inputs.input_ids.shape[1]:], skip_special_tokens=True)

        latency_ms = (time.time() - start) * 1000

        # Parse JSON
        try:
            data = json.loads(response_text)
        except json.JSONDecodeError:
            if "```json" in response_text:
                json_str = response_text.split("```json")[1].split("```")[0].strip()
                data = json.loads(json_str)
            else:
                data = {}

        result = LLMExtractionResult(
            organism=data.get("organism"),
            tissue=data.get("tissue"),
            cell_type=data.get("cell_type"),
            cell_line=data.get("cell_line"),
            treatment=data.get("treatment"),
            genotype=data.get("genotype"),
            strain=data.get("strain"),
            age=data.get("age"),
            sex=data.get("sex"),
            developmental_stage=data.get("developmental_stage"),
            confidence=data.get("confidence", {}),
            reasoning=data.get("reasoning"),
            model_name=self.model_path,
            tokens_used=outputs.shape[1],
            latency_ms=latency_ms,
        )

        return result

    def get_model_name(self) -> str:
        return self.model_path

    def supports_batching(self) -> bool:
        return False


# Factory function to create providers
def create_provider(
    provider_type: str = "claude",
    **kwargs
) -> LLMProvider:
    """
    Create an LLM provider.

    Args:
        provider_type: Type of provider ("claude", "gemini", "vllm", "transformers")
        **kwargs: Provider-specific arguments

    Returns:
        LLMProvider instance

    Examples:
        # Claude API
        provider = create_provider("claude", api_key="sk-ant-...")

        # Gemini via Vertex AI
        provider = create_provider("gemini", project="my-gcp-project", location="europe-west2")

        # VLLM on A100
        provider = create_provider("vllm", model_path="meta-llama/Llama-3.1-70B-Instruct")

        # Transformers (slower, but works everywhere)
        provider = create_provider("transformers", model_path="meta-llama/Llama-3.1-8B-Instruct")
    """
    if provider_type == "claude":
        return ClaudeProvider(**kwargs)
    elif provider_type == "gemini":
        return GeminiProvider(**kwargs)
    elif provider_type == "vllm":
        # For now, use Transformers as the local in-process default; server-vllm is available separately.
        return TransformersProvider(**kwargs)
    elif provider_type == "transformers":
        return TransformersProvider(**kwargs)
    elif provider_type == "server-vllm":
        return ServerVLLMProvider(**kwargs)
    else:
        raise ValueError(f"Unknown provider type: {provider_type}")
