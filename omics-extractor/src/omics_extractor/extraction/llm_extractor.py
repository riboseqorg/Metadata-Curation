"""LLM-based metadata extraction using Claude."""

from typing import Optional, Dict, Any, List
from pydantic import BaseModel, Field
import anthropic
import os
import json

from ..schemas.base import BaseProvenance, SampleMetadata
from ..ontologies.mapper import normalize_tissue, normalize_cell_type, normalize_organism


from pathlib import Path
import yaml
from .llm_schemas import LLMExtractionResult


def load_extraction_scheme(scheme_name: str) -> Dict[str, str]:
    """Load an extraction scheme from a YAML file."""
    # Schemes are located in the same directory as this file
    schemes_dir = Path(__file__).parent / "schemes"
    scheme_path = schemes_dir / f"{scheme_name}.yaml"
    
    if not scheme_path.exists():
        # Fallback to default if scheme doesn't exist
        scheme_path = schemes_dir / "default.yaml"
        
    if not scheme_path.exists():
        # Last resort - return a minimal hardcoded scheme
        return {
            "organism": "Scientific name of the organism",
            "tissue": "Tissue type",
            "cell_type": "Cell type"
        }
        
    with open(scheme_path, "r") as f:
        return yaml.safe_load(f)


def build_extraction_prompt(
    study_title: str,
    study_description: str,
    sample_title: Optional[str] = None,
    sample_description: Optional[str] = None,
    characteristics: Optional[Dict[str, str]] = None,
    run_metadata: Optional[Dict[str, Any]] = None,
    abstract: Optional[str] = None,
    journal: Optional[str] = None,
    authors: Optional[str] = None,
    publication_date: Optional[str] = None,
    existing_metadata: Optional[Dict[str, Any]] = None,
    fields: Optional[Dict[str, str]] = None,
) -> str:
    """
    Build a prompt for Claude to extract sample metadata.

    Args:
        study_title: BioProject title
        study_description: BioProject description
        sample_title: BioSample title
        sample_description: BioSample description
        characteristics: Raw sample characteristics/attributes
        run_metadata: SRA run metadata (library strategy, source, etc.)
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

    if journal or authors or publication_date:
        prompt += "PUBLICATION CONTEXT:\n"
        if journal:
            prompt += f"Journal: {journal}\n"
        if authors:
            prompt += f"Authors: {authors}\n"
        if publication_date:
            prompt += f"Date: {publication_date}\n"
        prompt += "\n"

    if sample_title:
        prompt += f"""SAMPLE INFORMATION:
Title: {sample_title}
"""

    if sample_description:
        prompt += f"""Description: {sample_description}
"""

    if characteristics or run_metadata:
        prompt += "\nSOURCE METADATA (raw fields from NCBI/GEO):\n"
        if characteristics:
            prompt += f"Biosample Attributes:\n{json.dumps(characteristics, indent=2)}\n"
        if run_metadata:
            prompt += f"Technical Context (SRA):\n{json.dumps(run_metadata, indent=2)}\n"

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

    # Use provided fields or default biosample attributes
    if not fields:
        fields = load_extraction_scheme("default")

    prompt += f"""
TASK:
Extract the following sample-level metadata fields. Only extract information you are confident about.
Return ONLY valid JSON with this structure:

{json.dumps(fields, indent=2)}

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
- Extract ALL relevant fields, not just the common ones

Return only the JSON object, no other text.
"""
    print(prompt)
    return prompt


def extract_with_claude(
    study_title: str,
    study_description: str,
    sample_title: Optional[str] = None,
    sample_description: Optional[str] = None,
    abstract: Optional[str] = None,
    journal: Optional[str] = None,
    authors: Optional[str] = None,
    publication_date: Optional[str] = None,
    characteristics: Optional[Dict[str, str]] = None,
    run_metadata: Optional[Dict[str, Any]] = None,
    existing_metadata: Optional[Dict[str, Any]] = None,
    api_key: Optional[str] = None,
    fields: Optional[Dict[str, str]] = None,
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
        characteristics=characteristics,
        run_metadata=run_metadata,
        abstract=abstract,
        journal=journal,
        authors=authors,
        publication_date=publication_date,
        existing_metadata=existing_metadata,
        fields=fields,
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
    journal: Optional[str] = None,
    authors: Optional[str] = None,
    publication_date: Optional[str] = None,
    provider: Optional[Any] = None,
    source_id: str = "llm",
    api_key: Optional[str] = None,
    scheme: str = "default",
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
    if sample.organism and sample.organism.value:
        existing["organism"] = sample.organism.value
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

    # Add raw characteristics if available (to provide full context to LLM)
    if sample.raw_characteristics:
        existing["raw_characteristics"] = sample.raw_characteristics

    # Load extraction scheme
    fields = load_extraction_scheme(scheme)

    # Extract with LLM (provider or Claude)
    if provider is not None:
        # Use generic provider interface (returns LLMExtractionResult directly)
        prompt = build_extraction_prompt(
            study_title=study_title,
            study_description=study_description,
            sample_title=sample_title,
            sample_description=sample_description,
            characteristics=sample.raw_characteristics,
            run_metadata=sample.technical_context,
            abstract=abstract,
            journal=journal,
            authors=authors,
            publication_date=publication_date,
            existing_metadata=existing,
            fields=fields,
        )
        llm_result = provider.extract(prompt)
    else:
        # Use Claude directly (backward compatibility)
        llm_result = extract_with_claude(
            study_title=study_title,
            study_description=study_description,
            sample_title=sample_title,
            sample_description=sample_description,
            characteristics=sample.raw_characteristics,
            run_metadata=sample.technical_context,
            abstract=abstract,
            journal=journal,
            authors=authors,
            publication_date=publication_date,
            existing_metadata=existing,
            api_key=api_key,
            fields=fields,
        )

    # Field confidence for LLM extraction (source field is study description/abstract)
    # This is lower because we're inferring from high-level descriptions
    llm_field_confidence = 0.5

    print(llm_result)

    # Apply ontology mapping to LLM-extracted values
    ontology_mappings = {}
    if llm_result.organism:
        normalized_org, org_id, onto_conf = normalize_organism(llm_result.organism)
        ontology_mappings["organism"] = (normalized_org, org_id, onto_conf)
        llm_result.organism = normalized_org

    if llm_result.tissue:
        normalized_tissue, tissue_id, onto_conf = normalize_tissue(llm_result.tissue)
        ontology_mappings["tissue"] = (normalized_tissue, tissue_id, onto_conf)
        llm_result.tissue = normalized_tissue  # Use normalized value

    if llm_result.cell_type:
        normalized_cell, cell_id, onto_conf = normalize_cell_type(llm_result.cell_type)
        ontology_mappings["cell_type"] = (normalized_cell, cell_id, onto_conf)
        llm_result.cell_type = normalized_cell  # Use normalized value

    # Add extracted fields that are missing
    # Fields to map from LLM result to sample
    mapping_fields = [
        "organism", "tissue", "cell_type", "cell_line", "strain", "genotype",
        "sex", "age", "developmental_stage", "condition", "treatment",
        "timepoint", "replicate", "batch", "disease", "stress",
        "temperature", "growth_condition", "library_strategy"
    ]

    for field in mapping_fields:
        llm_val = getattr(llm_result, field, None)
        existing_prov = getattr(sample, field, None)
        
        # Only enrich if sample doesn't have it or if it's empty
        if llm_val and (not existing_prov or not existing_prov.value):
            norm_conf = llm_result.confidence.get(field, 0.5)
            ontology_term = ontology_mappings.get(field, (None, None, None))[1]
            
            setattr(sample, field, BaseProvenance(
                value=llm_val,
                source="llm_enrichment",
                source_id=source_id,
                confidence=llm_field_confidence * norm_conf,
                field_confidence=llm_field_confidence,
                normalization_confidence=norm_conf,
                extraction_method="llm",
                extracted_text="Inferred from study context",
                notes=llm_result.reasoning,
                ontology_term=ontology_term,
            ))

    # Add any extra fields from LLM result to custom_fields
    if not sample.custom_fields:
        sample.custom_fields = {}
    
    # Get all fields defined in the scheme that were extracted but aren't in mapping_fields
    scheme_fields = set(fields.keys())
    standard_fields = set(mapping_fields)
    extra_fields = scheme_fields - standard_fields
    
    for field in extra_fields:
        llm_val = getattr(llm_result, field, None)
        if llm_val and field not in sample.custom_fields:
            sample.custom_fields[field] = llm_val

    return sample

