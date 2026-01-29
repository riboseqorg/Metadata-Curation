"""Enhanced sample extraction using field_mappings and optional LLM gap-filling.

This module provides two extraction modes:
1. Lightweight: Field mapping + ontologies (no LLM) - for production at scale
2. Enriched: Adds minimal LLM gap-filling - for samples with missing critical fields
"""

from typing import Dict, Optional
from ..schemas.base import BaseProvenance, SampleMetadata
from ..ontologies.mapper import normalize_tissue, normalize_cell_type, normalize_organism
from .field_mappings import extract_from_attributes, normalize_field_name


def extract_sample_with_field_mappings(
    biosample_meta,
    organism_from_sra: Optional[str] = None,
    sample_id: Optional[str] = None,
    bioproject_id: Optional[str] = None,
) -> SampleMetadata:
    """
    Extract sample metadata using field_mappings approach.

    This replaces hardcoded field lists with the flexible field_mappings dict.

    Args:
        biosample_meta: BioSampleMetadata object
        organism_from_sra: Optional organism from SRA (fallback)
        sample_id: Sample accession (if not biosample_id)
        bioproject_id: Parent BioProject accession

    Returns:
        SampleMetadata with extracted fields
    """
    # Start with organism (from BioSample or SRA)
    raw_organism = biosample_meta.organism or organism_from_sra or "unknown"
    organism_source = "biosample" if biosample_meta.organism else "sra"

    # Confidence scores for organism source
    organism_field_conf = 1.0 if biosample_meta.organism else (0.8 if organism_from_sra else 0.0)

    # Normalize organism
    normalized_organism, organism_ontology, organism_norm_conf = normalize_organism(raw_organism)
    organism_final_conf = organism_field_conf * organism_norm_conf

    # Create sample object
    sample = SampleMetadata(
        sample_id=sample_id or biosample_meta.biosample_id,
        bioproject_id=bioproject_id or "unknown",
        biosample_id=biosample_meta.biosample_id,
        organism=BaseProvenance(
            value=normalized_organism,
            source=organism_source,
            source_id=biosample_meta.biosample_id,
            confidence=organism_final_conf,
            field_confidence=organism_field_conf,
            normalization_confidence=organism_norm_conf,
            extraction_method="structured_field",
            ontology_term=organism_ontology,
        ),
        sample_title=BaseProvenance(
            value=biosample_meta.sample_name,
            source="biosample",
            source_id=biosample_meta.biosample_id,
            confidence=1.0,
            extraction_method="structured_field",
        ) if biosample_meta.sample_name else None,
        sample_description=BaseProvenance(
            value=biosample_meta.description,
            source="biosample",
            source_id=biosample_meta.biosample_id,
            confidence=1.0,
            extraction_method="structured_field",
        ) if biosample_meta.description else None,
        raw_characteristics=biosample_meta.attributes,
    )

    # Use field_mappings to extract standardized fields
    standardized = extract_from_attributes(biosample_meta.attributes)

    # Map each standardized field to SampleMetadata with confidence scoring
    _extract_with_confidence(sample, standardized, biosample_meta.biosample_id)

    return sample


def _extract_with_confidence(
    sample: SampleMetadata,
    standardized: Dict[str, str],
    biosample_id: str
):
    """
    Extract fields from standardized dict into SampleMetadata with confidence scores.

    Applies ontology mapping where appropriate and calculates two-factor confidence.
    """
    # TISSUE - with ontology mapping
    if "tissue" in standardized:
        raw_value = standardized["tissue"]
        normalized_value, ontology_id, norm_conf = normalize_tissue(raw_value)

        # Field confidence = 1.0 (came from structured BioSample field)
        field_conf = 1.0
        final_conf = field_conf * norm_conf

        sample.tissue = BaseProvenance(
            value=normalized_value,
            source="biosample",
            source_id=biosample_id,
            confidence=final_conf,
            field_confidence=field_conf,
            normalization_confidence=norm_conf,
            extracted_text=f"tissue: {raw_value}",
            ontology_term=ontology_id,
            ontology_label=normalized_value if ontology_id else None,
            extraction_method="structured_field",
        )

    # CELL_TYPE - with ontology mapping
    if "cell_type" in standardized:
        raw_value = standardized["cell_type"]
        normalized_value, ontology_id, norm_conf = normalize_cell_type(raw_value)

        field_conf = 1.0
        final_conf = field_conf * norm_conf

        sample.cell_type = BaseProvenance(
            value=normalized_value,
            source="biosample",
            source_id=biosample_id,
            confidence=final_conf,
            field_confidence=field_conf,
            normalization_confidence=norm_conf,
            extracted_text=f"cell_type: {raw_value}",
            ontology_term=ontology_id,
            ontology_label=normalized_value if ontology_id else None,
            extraction_method="structured_field",
        )

    # CELL_LINE - direct mapping
    if "cell_line" in standardized:
        sample.cell_line = BaseProvenance(
            value=standardized["cell_line"],
            source="biosample",
            source_id=biosample_id,
            confidence=1.0,
            extracted_text=f"cell_line: {standardized['cell_line']}",
            extraction_method="structured_field",
        )

    # STRAIN - direct mapping
    if "strain" in standardized:
        sample.strain = BaseProvenance(
            value=standardized["strain"],
            source="biosample",
            source_id=biosample_id,
            confidence=1.0,
            extracted_text=f"strain: {standardized['strain']}",
            extraction_method="structured_field",
        )

    # GENOTYPE - direct mapping
    if "genotype" in standardized:
        sample.genotype = BaseProvenance(
            value=standardized["genotype"],
            source="biosample",
            source_id=biosample_id,
            confidence=1.0,
            extracted_text=f"genotype: {standardized['genotype']}",
            extraction_method="structured_field",
        )

    # TREATMENT - direct mapping
    if "treatment" in standardized:
        sample.treatment = BaseProvenance(
            value=standardized["treatment"],
            source="biosample",
            source_id=biosample_id,
            confidence=1.0,
            extracted_text=f"treatment: {standardized['treatment']}",
            extraction_method="structured_field",
        )

    # AGE - direct mapping
    if "age" in standardized:
        sample.age = BaseProvenance(
            value=standardized["age"],
            source="biosample",
            source_id=biosample_id,
            confidence=1.0,
            extracted_text=f"age: {standardized['age']}",
            extraction_method="structured_field",
        )

    # SEX - direct mapping
    if "sex" in standardized:
        sample.sex = BaseProvenance(
            value=standardized["sex"],
            source="biosample",
            source_id=biosample_id,
            confidence=1.0,
            extracted_text=f"sex: {standardized['sex']}",
            extraction_method="structured_field",
        )

    # DEVELOPMENTAL_STAGE - direct mapping
    if "developmental_stage" in standardized:
        sample.developmental_stage = BaseProvenance(
            value=standardized["developmental_stage"],
            source="biosample",
            source_id=biosample_id,
            confidence=1.0,
            extracted_text=f"developmental_stage: {standardized['developmental_stage']}",
            extraction_method="structured_field",
        )


def needs_llm_enrichment(sample: SampleMetadata) -> bool:
    """
    Determine if a sample needs LLM gap-filling.

    Skip LLM if sample already has critical fields populated:
    - organism (always present)
    - tissue OR cell_line
    - strain (for model organisms)

    Args:
        sample: SampleMetadata object

    Returns:
        True if LLM enrichment would be beneficial, False otherwise
    """
    # Critical fields for RiboSeq metadata
    has_organism = bool(sample.organism and sample.organism.value)
    has_tissue_or_cell_line = bool(
        (sample.tissue and sample.tissue.value) or
        (sample.cell_line and sample.cell_line.value)
    )
    has_strain = bool(sample.strain and sample.strain.value)

    # Calculate completeness
    critical_fields_present = sum([
        1 if has_organism else 0,
        1 if has_tissue_or_cell_line else 0,
        1 if has_strain else 0
    ])

    # Need enrichment if less than 75% complete
    completeness = critical_fields_present / 3.0

    return completeness < 0.75


def build_minimal_llm_prompt(
    sample: SampleMetadata,
    study_title: str,
    study_description: str,
) -> str:
    """
    Build minimal LLM prompt for gap-filling.

    Only asks for fields that are missing. Keeps context small (<1000 tokens).

    Args:
        sample: SampleMetadata with baseline extraction
        study_title: Study title for context
        study_description: Study description for context

    Returns:
        Minimal prompt string
    """
    # Identify what we have
    existing = []
    if sample.organism:
        existing.append(f"organism: {sample.organism.value}")
    if sample.tissue:
        existing.append(f"tissue: {sample.tissue.value}")
    if sample.cell_line:
        existing.append(f"cell_line: {sample.cell_line.value}")
    if sample.cell_type:
        existing.append(f"cell_type: {sample.cell_type.value}")
    if sample.strain:
        existing.append(f"strain: {sample.strain.value}")
    if sample.treatment:
        existing.append(f"treatment: {sample.treatment.value}")

    # Identify what we need
    missing = []
    if not sample.tissue and not sample.cell_line:
        missing.append("tissue or cell_line")
    if not sample.cell_type:
        missing.append("cell_type")
    if not sample.strain:
        missing.append("strain")
    if not sample.treatment:
        missing.append("treatment")

    # Build minimal prompt
    prompt = f"""Study: {study_title}
Description: {study_description}

Sample: {sample.sample_title.value if sample.sample_title else 'N/A'}

Already extracted from structured fields:
{chr(10).join(existing) if existing else '(none)'}

Extract ONLY these missing fields: {', '.join(missing)}

Return JSON with just the missing fields. Only include fields you're confident about.
Use null for fields you can't determine.

{{
  "tissue": "value or null",
  "cell_type": "value or null",
  "strain": "value or null",
  "treatment": "value or null",
  "confidence": {{
    "tissue": 0.0-1.0,
    "cell_type": 0.0-1.0
  }}
}}
"""

    return prompt
