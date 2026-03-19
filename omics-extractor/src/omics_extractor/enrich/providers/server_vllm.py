from __future__ import annotations
from typing import List
import json
import os
import time
import urllib.request

from ...extraction.llm_schemas import LLMExtractionResult
from .base import LLMProvider

class ServerVLLMProvider(LLMProvider):
    def __init__(self, server_url: str, model: str, temperature: float = 0.0, max_tokens: int = 1000):
        self.server_url = server_url.rstrip("/")
        self.model = model
        self.temperature = temperature
        self.max_tokens = max_tokens

    def extract(self, prompt: str, temperature: float = 0.0, max_tokens: int = 1000) -> LLMExtractionResult:
        start = time.time()
        req = urllib.request.Request(
            url=f"{self.server_url}/generate",
            data=json.dumps({"model": self.model, "prompt": prompt, "temperature": temperature or self.temperature, "max_tokens": max_tokens or self.max_tokens}).encode("utf-8"),
            headers={"Content-Type": "application/json"},
            method="POST",
        )
        with urllib.request.urlopen(req, timeout=120) as resp:
            text = resp.read().decode("utf-8")
        latency_ms = (time.time() - start) * 1000
        try:
            data = json.loads(text)
        except json.JSONDecodeError:
            if "```json" in text:
                data = json.loads(text.split("```json")[1].split("```")[0].strip())
            else:
                data = {}
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
            tokens_used=None,
            latency_ms=latency_ms,
        )

    def get_model_name(self) -> str:
        return self.model

    def supports_batching(self) -> bool:
        return True
