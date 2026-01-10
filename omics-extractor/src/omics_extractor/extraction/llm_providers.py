"""LLM provider backends for metadata extraction.

This module provides a pluggable interface for different LLM providers,
allowing easy switching between:
- Claude API (Anthropic)
- Local models via VLLM (Llama, Mixtral, Qwen, etc.)
- Local models via Transformers
- Other APIs (OpenAI, etc.)

This enables benchmarking different models to find the best for your use case.
"""

from abc import ABC, abstractmethod
from typing import Dict, Any, Optional, List
from dataclasses import dataclass
import json
import os


@dataclass
class LLMExtractionResult:
    """Result from LLM extraction."""

    tissue: Optional[str] = None
    cell_type: Optional[str] = None
    cell_line: Optional[str] = None
    treatment: Optional[str] = None
    genotype: Optional[str] = None
    strain: Optional[str] = None
    age: Optional[str] = None
    sex: Optional[str] = None
    developmental_stage: Optional[str] = None

    # Per-field confidence scores (0-1)
    confidence: Dict[str, float] = None

    # Optional reasoning/explanation
    reasoning: Optional[str] = None

    # Metadata about the extraction
    model_name: Optional[str] = None
    tokens_used: Optional[int] = None
    latency_ms: Optional[float] = None

    def __post_init__(self):
        if self.confidence is None:
            self.confidence = {}


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
        model: str = "claude-sonnet-4-5-20250929",
    ):
        """
        Initialize Claude provider.

        Args:
            api_key: Anthropic API key (or from ANTHROPIC_API_KEY env)
            model: Claude model to use (default: claude-sonnet-4-5-20250929 - Claude 4.5 Sonnet)
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


class VLLMProvider(LLMProvider):
    """VLLM provider for local models on GPU (A100, etc.)."""

    def __init__(
        self,
        model_path: str = "meta-llama/Llama-3.1-70B-Instruct",
        tensor_parallel_size: int = 1,
        gpu_memory_utilization: float = 0.9,
    ):
        """
        Initialize VLLM provider.

        Args:
            model_path: HuggingFace model path or local path
            tensor_parallel_size: Number of GPUs for tensor parallelism
            gpu_memory_utilization: Fraction of GPU memory to use
        """
        self.model_path = model_path

        # Lazy import VLLM (only needed for local models)
        try:
            from vllm import LLM, SamplingParams
            self.LLM = LLM
            self.SamplingParams = SamplingParams
        except ImportError:
            raise ImportError(
                "VLLM not installed. Install with: pip install vllm\n"
                "For A100 deployment, VLLM provides the best performance."
            )

        # Initialize model (this caches after first load)
        print(f"Loading {model_path} with VLLM (this may take a few minutes)...")
        self.llm = self.LLM(
            model=model_path,
            tensor_parallel_size=tensor_parallel_size,
            gpu_memory_utilization=gpu_memory_utilization,
            trust_remote_code=True,  # For some models like Qwen
        )
        print(f"Model loaded successfully!")

    def extract(
        self,
        prompt: str,
        temperature: float = 0.0,
        max_tokens: int = 1000,
    ) -> LLMExtractionResult:
        """Extract using VLLM (single prompt)."""
        results = self.extract_batch([prompt], temperature, max_tokens)
        return results[0]

    def extract_batch(
        self,
        prompts: List[str],
        temperature: float = 0.0,
        max_tokens: int = 1000,
    ) -> List[LLMExtractionResult]:
        """Extract using VLLM (batched for efficiency)."""
        import time

        start = time.time()

        # VLLM sampling parameters
        sampling_params = self.SamplingParams(
            temperature=temperature,
            max_tokens=max_tokens,
            top_p=0.95 if temperature > 0 else 1.0,
        )

        # Generate (VLLM automatically batches)
        outputs = self.llm.generate(prompts, sampling_params)

        latency_ms = (time.time() - start) * 1000
        avg_latency = latency_ms / len(prompts)

        # Parse outputs
        results = []
        for output in outputs:
            response_text = output.outputs[0].text

            try:
                data = json.loads(response_text)
            except json.JSONDecodeError:
                # Try to extract JSON from markdown
                if "```json" in response_text:
                    json_str = response_text.split("```json")[1].split("```")[0].strip()
                    data = json.loads(json_str)
                else:
                    # Failed to parse - return empty result
                    data = {}

            result = LLMExtractionResult(
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
                tokens_used=len(output.outputs[0].token_ids),
                latency_ms=avg_latency,
            )
            results.append(result)

        return results

    def get_model_name(self) -> str:
        return self.model_path

    def supports_batching(self) -> bool:
        return True  # VLLM has excellent batching support


class TransformersProvider(LLMProvider):
    """HuggingFace Transformers provider (CPU/GPU, slower than VLLM)."""

    def __init__(
        self,
        model_path: str = "meta-llama/Llama-3.1-8B-Instruct",
        device: str = "cuda",
        load_in_8bit: bool = False,
    ):
        """
        Initialize Transformers provider.

        Args:
            model_path: HuggingFace model path
            device: Device to run on (cuda/cpu)
            load_in_8bit: Use 8-bit quantization for memory efficiency
        """
        self.model_path = model_path
        self.device = device

        try:
            from transformers import AutoTokenizer, AutoModelForCausalLM
            import torch
        except ImportError:
            raise ImportError("transformers not installed. Install with: pip install transformers torch")

        print(f"Loading {model_path} with Transformers...")
        self.tokenizer = AutoTokenizer.from_pretrained(model_path)
        self.model = AutoModelForCausalLM.from_pretrained(
            model_path,
            load_in_8bit=load_in_8bit,
            device_map="auto" if device == "cuda" else None,
            torch_dtype=torch.float16 if device == "cuda" else torch.float32,
        )
        print(f"Model loaded!")

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
        return False  # Could implement, but VLLM is better for batching


# Factory function to create providers
def create_provider(
    provider_type: str = "claude",
    **kwargs
) -> LLMProvider:
    """
    Create an LLM provider.

    Args:
        provider_type: Type of provider ("claude", "vllm", "transformers")
        **kwargs: Provider-specific arguments

    Returns:
        LLMProvider instance

    Examples:
        # Claude API
        provider = create_provider("claude", api_key="sk-ant-...")

        # VLLM on A100
        provider = create_provider("vllm", model_path="meta-llama/Llama-3.1-70B-Instruct")

        # Transformers (slower, but works everywhere)
        provider = create_provider("transformers", model_path="meta-llama/Llama-3.1-8B-Instruct")
    """
    if provider_type == "claude":
        return ClaudeProvider(**kwargs)
    elif provider_type == "vllm":
        return VLLMProvider(**kwargs)
    elif provider_type == "transformers":
        return TransformersProvider(**kwargs)
    else:
        raise ValueError(f"Unknown provider type: {provider_type}")
