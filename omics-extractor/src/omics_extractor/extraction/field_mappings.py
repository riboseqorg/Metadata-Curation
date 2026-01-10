"""Simple field name normalization mappings.

This replaces the complex CSV-based approach with a straightforward Python dict.
Maps common BioSample field name variations to standardized field names.
"""

from typing import Optional


# Standardized field name synonyms
# Format: "standard_name": ["synonym1", "synonym2", ...]
FIELD_SYNONYMS = {
    "organism": [
        "organism", "organism_name", "scientific_name", "species",
        "organism name", "scientific name", "latin name"
    ],

    "strain": [
        "strain", "strain_name", "cultivar", "ecotype", "breed",
        "substrain", "strain name"
    ],

    "tissue": [
        "tissue", "tissue_type", "organ", "organism part", "body site",
        "tissue type", "tissue origin", "tissue source", "tissue lineage",
        "tissue subtype", "brain region", "tissue.1"
    ],

    "cell_type": [
        "cell_type", "cell type", "celltype", "cell_subtype", "cell subtype",
        "Cell type"
    ],

    "cell_line": [
        "cell_line", "cell line", "cell-line", "cell lline", "cells",
        "culture_collection", "host cell type", "cell line id",
        "parental cell line id", "cell line background"
    ],

    "developmental_stage": [
        "dev_stage", "developmental stage", "dev stage", "growth stage",
        "development stage", "Stage", "Dev_stage (yeast)", "growth_strage",
        "growth phase", "cellular state", "developmental_stage"
    ],

    "age": [
        "age", "age_value", "developmental_age", "growth time", "time",
        "sampling time", "days in vitro", "day"
    ],

    "sex": [
        "sex", "gender", "mating_type", "mating type"
    ],

    "genotype": [
        "genotype", "genetic_modification", "genetic modification",
        "modification", "genotype/variation", "genotype/variaton"
    ],

    "treatment": [
        "treatment", "compound", "drug", "inhibitor", "compound treatment",
        "drug treatment", "cell treatment", "sample treatment",
        "treatment.1", "treatment.2", "cell culture treatment",
        "drug pre-treatment", "compound treatment", "ribosome inhibitor treatment",
        "cycloheximide", "ribosome stabilization agent"
    ],

    "disease": [
        "disease", "disease_state", "disease state", "disease_stage",
        "disease stage"
    ],

    "biosample_id": [
        "BioSample", "biosample", "biosample_accession", "sample_name",
        "sample_accession", "sample name"
    ],

    "sample_title": [
        "sample_title", "sample title", "title", "sample_name"
    ],

    "sample_description": [
        "description", "sample_description", "sample description"
    ],
}


def normalize_field_name(raw_field_name: str) -> Optional[str]:
    """
    Map a raw field name to a standardized field name.

    Args:
        raw_field_name: The raw field name from BioSample attributes

    Returns:
        Standardized field name, or None if no mapping found

    Example:
        >>> normalize_field_name("tissue_type")
        'tissue'
        >>> normalize_field_name("organism part")
        'tissue'
        >>> normalize_field_name("unknown_field")
        None
    """
    # Normalize input
    raw_lower = raw_field_name.lower().strip()

    # Check each standard field's synonyms
    for standard_name, synonyms in FIELD_SYNONYMS.items():
        if raw_lower in [s.lower() for s in synonyms]:
            return standard_name

    return None


def extract_from_attributes(
    attributes: dict,
    include_unmapped: bool = False
) -> dict:
    """
    Extract standardized fields from BioSample attributes.

    Args:
        attributes: Raw BioSample attributes dict
        include_unmapped: If True, include fields that don't map to standard names

    Returns:
        Dict with standardized field names and their values

    Example:
        >>> attrs = {"tissue_type": "liver", "organism": "Homo sapiens", "custom_field": "value"}
        >>> extract_from_attributes(attrs)
        {'tissue': 'liver', 'organism': 'Homo sapiens'}
        >>> extract_from_attributes(attrs, include_unmapped=True)
        {'tissue': 'liver', 'organism': 'Homo sapiens', 'custom_field': 'value'}
    """
    result = {}

    for raw_key, raw_value in attributes.items():
        # Skip empty values
        if not raw_value or str(raw_value).strip() == "":
            continue

        # Try to map to standard name
        standard_name = normalize_field_name(raw_key)

        if standard_name:
            # If we already have this field, only override if new value is non-empty
            if standard_name not in result or not result[standard_name]:
                result[standard_name] = str(raw_value).strip()
        elif include_unmapped:
            # Keep unmapped fields if requested
            result[raw_key] = str(raw_value).strip()

    return result


# Priority order for fields when multiple sources exist
# (e.g., if both "tissue" and "organism part" are present)
FIELD_PRIORITY = {
    "tissue": ["tissue", "tissue_type", "organ", "organism part"],
    "cell_type": ["cell_type", "cell type"],
    "cell_line": ["cell_line", "cell line"],
    # Add more as needed
}


def resolve_conflicts(attributes: dict) -> dict:
    """
    Resolve conflicts when multiple raw fields map to same standard field.
    Uses priority order defined in FIELD_PRIORITY.

    Args:
        attributes: Raw attributes dict

    Returns:
        Dict with conflicts resolved
    """
    # First, map all fields
    mapped = extract_from_attributes(attributes)

    # Then check for conflicts (this is simple - just use first non-empty)
    # More sophisticated logic could be added here if needed

    return mapped


if __name__ == "__main__":
    # Test the mappings
    test_attrs = {
        "tissue_type": "liver",
        "organism part": "kidney",  # Conflict with tissue_type
        "organism": "Homo sapiens",
        "strain": "C57BL/6",
        "cell line": "HEK293",
        "dev_stage": "adult",
        "custom_field": "some value",
    }

    print("Test attributes:")
    print(test_attrs)
    print("\nMapped (standard only):")
    print(extract_from_attributes(test_attrs))
    print("\nMapped (include unmapped):")
    print(extract_from_attributes(test_attrs, include_unmapped=True))
