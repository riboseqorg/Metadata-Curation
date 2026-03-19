from __future__ import annotations
from abc import ABC, abstractmethod
from typing import List
from ...extraction.llm_schemas import LLMExtractionResult

class LLMProvider(ABC):
    @abstractmethod
    def extract(self, prompt: str, temperature: float = 0.0, max_tokens: int = 1000) -> LLMExtractionResult:
        ...

    @abstractmethod
    def get_model_name(self) -> str:
        ...

    @abstractmethod
    def supports_batching(self) -> bool:
        ...

    def extract_batch(self, prompts: List[str], temperature: float = 0.0, max_tokens: int = 1000) -> List[LLMExtractionResult]:
        return [self.extract(p, temperature, max_tokens) for p in prompts]
