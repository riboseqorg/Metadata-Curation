"""LLM-based metadata extraction using Claude."""

from typing import Optional, Dict, Any, List
from pydantic import BaseModel, Field
import anthropic
import os
import json

from ..schemas.base import BaseProvenance, SampleMetadata


class LLMExtractionResult(BaseModel):
    """Result from LLM extraction."""

    tissue: Optional[str] = Field(None, description="Tissue type")
    cell_type: Optional[str] = Field(None, description="Cell type")
    cell_line: Optional[str] = Field(None, description="Cell line")
    treatment: Optional[str] = Field(None, description="Treatment/condition")
    genotype: Optional[str] = Field(None, description="Genotype")
    strain: Optional[str] = Field(None, description="Strain")
    age: Optional[str] = Field(None, description="Age")
    sex: Optional[str] = Field(None, description="Sex")
    developmental_stage: Optional[str] = Field(None, description="Developmental stage")

    # Confidence for each field
    confidence: Dict[str, float] = Field(default_factory=dict, description="Confidence per field")
    reasoning: Optional[str] = Field(None, description="LLM's reasoning")


def build_extraction_prompt(
    study_title: str,
    study_description: str,
    sample_title: Optional[str] = None,
    sample_description: Optional[str] = None,
    abstract: Optional[str] = None,
    existing_metadata: Optional[Dict[str, Any]] = None,
) -> str:
    """
    Build a prompt for Claude to extract sample metadata.

    Args:
        study_title: BioProject title
        study_description: BioProject description
        sample_title: BioSample title
        sample_description: BioSample description
        abstract: Publication abstract (optional)
        existing_metadata: Already extracted metadata to avoid redundancy

    Returns:
        Formatted prompt for Claude
    """
    prompt = f"""You are a bioinformatics metadata curator extracting sample-level metadata from study descriptions.

STUDY INFORMATION:
Title: {study_title}

Description:
{study_description}

"""

    if sample_title:
        prompt += f"""SAMPLE INFORMATION:
Title: {sample_title}
"""

    if sample_description:
        prompt += f"""Description: {sample_description}
"""

    if abstract:
        prompt += f"""
PUBLICATION ABSTRACT:
{abstract}
"""

    if existing_metadata:
        prompt += f"""
ALREADY EXTRACTED (do not duplicate):
{json.dumps(existing_metadata, indent=2)}
"""

    prompt += """
TASK:
Extract the following sample-level metadata fields. Only extract information you are confident about.
Return ONLY valid JSON with this structure:

{
  "tissue": "tissue type or null",
  "cell_type": "cell type or null",
  "cell_line": "cell line name or null",
  "treatment": "treatment/drug/condition or null",
  "genotype": "genotype/mutation or null",
  "strain": "strain name or null",
  "age": "age or null",
  "sex": "male/female or null",
  "developmental_stage": "developmental stage or null",
  "confidence": {
    "tissue": 0.0-1.0,
    "cell_type": 0.0-1.0
  },
  "reasoning": "brief explanation of extraction"
}

CONFIDENCE SCORING GUIDELINES:
- 1.0: Explicitly stated in structured field
- 0.8: Clearly implied from context
- 0.6: Inferred with reasonable certainty
- 0.4: Educated guess from limited information
- 0.0: No information available

IMPORTANT:
- Use lowercase, standardized terms
- If a field is not mentioned or cannot be inferred, use null
- Be conservative - only extract what you're confident about
- Provide confidence scores for each non-null field

Return only the JSON object, no other text.
"""

    return prompt


def extract_with_claude(
    study_title: str,
    study_description: str,
    sample_title: Optional[str] = None,
    sample_description: Optional[str] = None,
    abstract: Optional[str] = None,
    existing_metadata: Optional[Dict[str, Any]] = None,
    api_key: Optional[str] = None,
) -> LLMExtractionResult:
    """
    Extract metadata using Claude API.

    Args:
        study_title: BioProject title
        study_description: BioProject description
        sample_title: BioSample title
        sample_description: BioSample description
        abstract: Publication abstract
        existing_metadata: Already extracted fields
        api_key: Anthropic API key (or from ANTHROPIC_API_KEY env var)

    Returns:
        LLMExtractionResult with extracted metadata

    Raises:
        ValueError: If API key not provided
    """
    if api_key is None:
        api_key = os.environ.get("ANTHROPIC_API_KEY")

    if not api_key:
        raise ValueError("Anthropic API key required. Set ANTHROPIC_API_KEY environment variable or pass api_key parameter.")

    client = anthropic.Anthropic(api_key=api_key)

    prompt = build_extraction_prompt(
        study_title=study_title,
        study_description=study_description,
        sample_title=sample_title,
        sample_description=sample_description,
        abstract=abstract,
        existing_metadata=existing_metadata,
    )

    message = client.messages.create(
        model="claude-3-5-sonnet-20241022",
        max_tokens=2048,
        temperature=0,  # Deterministic for consistency
        messages=[
            {"role": "user", "content": prompt}
        ]
    )

    # Parse JSON response
    response_text = message.content[0].text

    # Extract JSON from markdown code blocks if present
    if "```json" in response_text:
        response_text = response_text.split("```json")[1].split("```")[0].strip()
    elif "```" in response_text:
        response_text = response_text.split("```")[1].split("```")[0].strip()

    result_dict = json.loads(response_text)

    return LLMExtractionResult(**result_dict)


def enrich_sample_metadata(
    sample: SampleMetadata,
    study_title: str,
    study_description: str,
    sample_title: Optional[str] = None,
    sample_description: Optional[str] = None,
    abstract: Optional[str] = None,
    api_key: Optional[str] = None,
    provider = None,  # LLMProvider instance
    source_id: str = "llm_extraction",
) -> SampleMetadata:
    """
    Enrich existing sample metadata with LLM extraction.

    Only extracts fields that are currently missing (None).
    Uses two-factor confidence: field (0.5 for LLM) × normalization (LLM's confidence).

    Args:
        sample: Existing SampleMetadata object
        study_title: Study title for context
        study_description: Study description for context
        sample_title: Sample title
        sample_description: Sample description
        abstract: Publication abstract (optional)
        api_key: Anthropic API key (for Claude, if provider not provided)
        provider: LLMProvider instance (if provided, uses this instead of api_key)
        source_id: Source identifier for provenance

    Returns:
        Enriched SampleMetadata with LLM-extracted fields
    """
    # Build dict of already-extracted fields
    existing = {}
    if sample.tissue and sample.tissue.value:
        existing["tissue"] = sample.tissue.value
    if sample.cell_type and sample.cell_type.value:
        existing["cell_type"] = sample.cell_type.value
    if sample.cell_line and sample.cell_line.value:
        existing["cell_line"] = sample.cell_line.value
    if sample.treatment and sample.treatment.value:
        existing["treatment"] = sample.treatment.value
    if sample.genotype and sample.genotype.value:
        existing["genotype"] = sample.genotype.value
    if sample.strain and sample.strain.value:
        existing["strain"] = sample.strain.value
    if sample.age and sample.age.value:
        existing["age"] = sample.age.value
    if sample.sex and sample.sex.value:
        existing["sex"] = sample.sex.value
    if sample.developmental_stage and sample.developmental_stage.value:
        existing["developmental_stage"] = sample.developmental_stage.value

    # Extract with LLM (provider or Claude)
    if provider is not None:
        # Use generic provider interface (returns LLMExtractionResult directly)
        prompt = build_extraction_prompt(
            study_title=study_title,
            study_description=study_description,
            sample_title=sample_title,
            sample_description=sample_description,
            abstract=abstract,
            existing_metadata=existing,
        )
        llm_result = provider.extract(prompt)
    else:
        # Use Claude directly (backward compatibility)
        llm_result = extract_with_claude(
            study_title=study_title,
            study_description=study_description,
            sample_title=sample_title,
            sample_description=sample_description,
            abstract=abstract,
            existing_metadata=existing,
            api_key=api_key,
        )

    # Field confidence for LLM extraction (source field is study description/abstract)
    # This is lower because we're inferring from high-level descriptions
    llm_field_confidence = 0.5

    # Add extracted fields that are missing
    if llm_result.tissue and not sample.tissue:
        norm_conf = llm_result.confidence.get("tissue", 0.5)
        sample.tissue = BaseProvenance(
            value=llm_result.tissue,
            source="llm",
            source_id=source_id,
            confidence=llm_field_confidence * norm_conf,
            field_confidence=llm_field_confidence,
            normalization_confidence=norm_conf,
            extraction_method="llm",
            extracted_text=f"Inferred from study context",
            notes=llm_result.reasoning,
        )

    if llm_result.cell_type and not sample.cell_type:
        norm_conf = llm_result.confidence.get("cell_type", 0.5)
        sample.cell_type = BaseProvenance(
            value=llm_result.cell_type,
            source="llm",
            source_id=source_id,
            confidence=llm_field_confidence * norm_conf,
            field_confidence=llm_field_confidence,
            normalization_confidence=norm_conf,
            extraction_method="llm",
            extracted_text=f"Inferred from study context",
            notes=llm_result.reasoning,
        )

    if llm_result.cell_line and not sample.cell_line:
        norm_conf = llm_result.confidence.get("cell_line", 0.5)
        sample.cell_line = BaseProvenance(
            value=llm_result.cell_line,
            source="llm",
            source_id=source_id,
            confidence=llm_field_confidence * norm_conf,
            field_confidence=llm_field_confidence,
            normalization_confidence=norm_conf,
            extraction_method="llm",
            extracted_text=f"Inferred from study context",
            notes=llm_result.reasoning,
        )

    if llm_result.treatment and not sample.treatment:
        norm_conf = llm_result.confidence.get("treatment", 0.5)
        sample.treatment = BaseProvenance(
            value=llm_result.treatment,
            source="llm",
            source_id=source_id,
            confidence=llm_field_confidence * norm_conf,
            field_confidence=llm_field_confidence,
            normalization_confidence=norm_conf,
            extraction_method="llm",
            extracted_text=f"Inferred from study context",
            notes=llm_result.reasoning,
        )

    if llm_result.genotype and not sample.genotype:
        norm_conf = llm_result.confidence.get("genotype", 0.5)
        sample.genotype = BaseProvenance(
            value=llm_result.genotype,
            source="llm",
            source_id=source_id,
            confidence=llm_field_confidence * norm_conf,
            field_confidence=llm_field_confidence,
            normalization_confidence=norm_conf,
            extraction_method="llm",
            extracted_text=f"Inferred from study context",
            notes=llm_result.reasoning,
        )

    if llm_result.strain and not sample.strain:
        norm_conf = llm_result.confidence.get("strain", 0.5)
        sample.strain = BaseProvenance(
            value=llm_result.strain,
            source="llm",
            source_id=source_id,
            confidence=llm_field_confidence * norm_conf,
            field_confidence=llm_field_confidence,
            normalization_confidence=norm_conf,
            extraction_method="llm",
            extracted_text=f"Inferred from study context",
            notes=llm_result.reasoning,
        )

    if llm_result.age and not sample.age:
        norm_conf = llm_result.confidence.get("age", 0.5)
        sample.age = BaseProvenance(
            value=llm_result.age,
            source="llm",
            source_id=source_id,
            confidence=llm_field_confidence * norm_conf,
            field_confidence=llm_field_confidence,
            normalization_confidence=norm_conf,
            extraction_method="llm",
            extracted_text=f"Inferred from study context",
            notes=llm_result.reasoning,
        )

    if llm_result.sex and not sample.sex:
        norm_conf = llm_result.confidence.get("sex", 0.5)
        sample.sex = BaseProvenance(
            value=llm_result.sex,
            source="llm",
            source_id=source_id,
            confidence=llm_field_confidence * norm_conf,
            field_confidence=llm_field_confidence,
            normalization_confidence=norm_conf,
            extraction_method="llm",
            extracted_text=f"Inferred from study context",
            notes=llm_result.reasoning,
        )

    if llm_result.developmental_stage and not sample.developmental_stage:
        norm_conf = llm_result.confidence.get("developmental_stage", 0.5)
        sample.developmental_stage = BaseProvenance(
            value=llm_result.developmental_stage,
            source="llm",
            source_id=source_id,
            confidence=llm_field_confidence * norm_conf,
            field_confidence=llm_field_confidence,
            normalization_confidence=norm_conf,
            extraction_method="llm",
            extracted_text=f"Inferred from study context",
            notes=llm_result.reasoning,
        )

    return sample
