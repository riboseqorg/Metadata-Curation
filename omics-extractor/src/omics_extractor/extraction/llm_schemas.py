"""Schemas for LLM extraction results."""

from typing import Optional, Dict, List
from pydantic import BaseModel, Field

class LLMExtractionResult(BaseModel):
    """Result from LLM extraction."""

    # Core biological metadata
    organism: Optional[str] = Field(None, description="Scientific name of the organism")
    tissue: Optional[str] = Field(None, description="Tissue type")
    cell_type: Optional[str] = Field(None, description="Cell type")
    cell_line: Optional[str] = Field(None, description="Cell line")
    developmental_stage: Optional[str] = Field(None, description="Developmental stage")
    strain: Optional[str] = Field(None, description="Strain/variety")
    genotype: Optional[str] = Field(None, description="Genotype/genetic background")
    sex: Optional[str] = Field(None, description="Biological sex")
    age: Optional[str] = Field(None, description="Age")

    # Experimental conditions
    condition: Optional[str] = Field(None, description="Experimental condition")
    treatment: Optional[str] = Field(None, description="Treatment/drug applied")
    timepoint: Optional[str] = Field(None, description="Time point")
    replicate: Optional[str] = Field(None, description="Biological replicate")
    batch: Optional[str] = Field(None, description="Batch number")

    # Disease/perturbation
    disease: Optional[str] = Field(None, description="Disease state")
    stress: Optional[str] = Field(None, description="Stress condition")
    temperature: Optional[str] = Field(None, description="Temperature")
    growth_condition: Optional[str] = Field(None, description="Growth conditions")

    # Confidence and reasoning
    confidence: Dict[str, float] = Field(default_factory=dict, description="Confidence per field")
    reasoning: Optional[str] = Field(None, description="LLM's reasoning")

    # Metadata about the extraction
    model_name: Optional[str] = Field(None, description="Name of the model used")
    tokens_used: Optional[int] = Field(None, description="Total tokens used (input + output)")
    latency_ms: Optional[float] = Field(None, description="Extraction latency in milliseconds")
