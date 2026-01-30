"""Traceable output format with all levels of detail in one structure.

This creates a hierarchical output where users can access:
1. Quick view (top level) - just the values
2. Standard view (one level down) - values + basic provenance
3. Detailed view (deep level) - complete extraction chain

All in the same output file, fully traceable from raw data to final value.
"""

from typing import Dict, Any, List, Optional
from datetime import datetime
from ..schemas.base import BaseProvenance, SampleMetadata, StudyMetadata, RunMetadata


def format_field_traceable(prov: BaseProvenance) -> Dict[str, Any]:
    """
    Format a single field with complete traceability.

    Returns hierarchical structure:
    {
        "value": "testis",  # <- Quick access
        "provenance": {     # <- Standard access
            "source": "biosample",
            "confidence": 0.95,
            "ontology_term": "UBERON:0000473"
        },
        "extraction_details": {  # <- Deep access for debugging
            "field_confidence": 1.0,
            "normalization_confidence": 0.95,
            "extraction_method": "structured_field",
            "timestamp": "2025-01-04T10:30:00",
            "original_text": "testis",
            "ontology_label": "testis",
            "notes": null
        },
        "raw_biosample_attributes": {...}  # <- Original source data
    }
    """
    output = {
        # Level 1: Quick access - just the value
        "value": prov.value,

        # Level 2: Standard provenance
        "provenance": {
            "source": prov.source,
            "source_id": prov.source_id,
            "confidence": round(prov.confidence, 3),
            "ontology_term": prov.ontology_term,
        },
    }

    # Level 3: Detailed extraction chain
    extraction_details = {}

    if prov.field_confidence is not None:
        extraction_details["field_confidence"] = round(prov.field_confidence, 3)
        extraction_details["field_confidence_explanation"] = _explain_field_confidence(prov)

    if prov.normalization_confidence is not None:
        extraction_details["normalization_confidence"] = round(prov.normalization_confidence, 3)
        extraction_details["normalization_explanation"] = _explain_normalization(prov)

    if prov.extraction_method:
        extraction_details["method"] = prov.extraction_method

    if prov.extraction_timestamp:
        extraction_details["timestamp"] = prov.extraction_timestamp.isoformat()

    if prov.extracted_text and prov.extracted_text != prov.value:
        extraction_details["original_text"] = prov.extracted_text
        extraction_details["transformation"] = f"{prov.extracted_text} → {prov.value}"

    if prov.ontology_label:
        extraction_details["ontology_label"] = prov.ontology_label

    if prov.notes:
        extraction_details["notes"] = prov.notes

    if extraction_details:
        output["extraction_details"] = extraction_details

    return output


def _explain_field_confidence(prov: BaseProvenance) -> str:
    """Generate human-readable explanation of field confidence."""
    if prov.field_confidence >= 0.95:
        return "Extracted from specific, high-quality field"
    elif prov.field_confidence >= 0.8:
        return "Extracted from reliable field"
    elif prov.field_confidence >= 0.5:
        return "Inferred from less specific field or LLM"
    else:
        return "Low confidence source - manual review recommended"


def _explain_normalization(prov: BaseProvenance) -> str:
    """Generate human-readable explanation of normalization."""
    if prov.normalization_confidence >= 0.95:
        return "Exact or near-exact match to standard term"
    elif prov.normalization_confidence >= 0.8:
        return "Good synonym match"
    elif prov.normalization_confidence >= 0.5:
        return "Partial match or derived value"
    else:
        return "Significant transformation required"


def format_sample_traceable(sample: SampleMetadata, include_raw: bool = True) -> Dict[str, Any]:
    """
    Format sample with complete traceability.

    Structure:
    {
        "quick_view": {...},           # Just values
        "identifiers": {...},
        "biological_metadata": {...},  # Full traceable format
        "experimental_metadata": {...},
        "raw_data": {...}              # Original fetched data (optional)
    }
    """
    output = {}

    # Level 1: Quick view - just values for fast access
    quick_view = {
        "sample_id": sample.sample_id,
        "organism": sample.organism.value if sample.organism else None,
    }

    for field in ["tissue", "cell_line", "cell_type", "strain", "genotype", "sex", "age",
                  "developmental_stage", "treatment", "condition", "disease", "timepoint",
                  "replicate", "batch", "stress", "temperature", "growth_condition", "library_strategy"]:
        attr = getattr(sample, field, None)
        if attr and isinstance(attr, BaseProvenance):
            quick_view[field] = attr.value

    output["quick_view"] = quick_view

    # Identifiers (always simple)
    output["identifiers"] = {
        "sample_id": sample.sample_id,
        "bioproject_id": sample.bioproject_id,
        "biosample_id": sample.biosample_id,
        "gsm_id": sample.gsm_id,
    }

    # Biological metadata - full traceable format
    output["biological_metadata"] = {
        "organism": format_field_traceable(sample.organism),
    }

    for field in ["tissue", "cell_line", "cell_type", "developmental_stage",
                  "strain", "genotype", "sex", "age"]:
        attr = getattr(sample, field, None)
        if attr and isinstance(attr, BaseProvenance):
            output["biological_metadata"][field] = format_field_traceable(attr)

    # Experimental metadata - full traceable format
    experimental = {}
    for field in ["condition", "treatment", "timepoint", "replicate", "batch"]:
        attr = getattr(sample, field, None)
        if attr and isinstance(attr, BaseProvenance):
            experimental[field] = format_field_traceable(attr)

    if experimental:
        output["experimental_metadata"] = experimental

    # Technical metadata
    if sample.library_strategy:
        output["library_strategy"] = format_field_traceable(sample.library_strategy)

    # Disease/perturbation metadata
    disease_perturbation = {}
    for field in ["disease", "stress", "temperature", "growth_condition"]:
        attr = getattr(sample, field, None)
        if attr and isinstance(attr, BaseProvenance):
            disease_perturbation[field] = format_field_traceable(attr)

    if disease_perturbation:
        output["disease_perturbation"] = disease_perturbation

    # Sample descriptions
    if sample.sample_title or sample.sample_description:
        output["descriptions"] = {}
        if sample.sample_title:
            output["descriptions"]["title"] = format_field_traceable(sample.sample_title)
        if sample.sample_description:
            output["descriptions"]["description"] = format_field_traceable(sample.sample_description)

    # Raw characteristics (original data from BioSample/GEO before processing)
    if include_raw and sample.raw_characteristics:
        output["raw_biosample_attributes"] = sample.raw_characteristics

    # Technical context for LLM
    if sample.technical_context:
        output["technical_context"] = sample.technical_context

    # Custom fields
    if sample.custom_fields:
        output["custom_fields"] = sample.custom_fields

    return output


def format_study_traceable(study: StudyMetadata) -> Dict[str, Any]:
    """Format study with complete traceability."""
    output = {
        # Quick view
        "quick_view": {
            "bioproject_id": study.bioproject_id,
            "title": study.title.value,
            "organism": study.organism.value,
        },

        # Identifiers
        "identifiers": {
            "bioproject_id": study.bioproject_id,
            "gse_id": study.gse_id,
            "sra_study_id": study.sra_study_id,
            "pmid": study.pmid,
            "pmc_id": study.pmc_id,
        },

        # Core metadata with full traceability
        "core_metadata": {
            "title": format_field_traceable(study.title),
            "description": format_field_traceable(study.description),
            "organism": format_field_traceable(study.organism),
        },
    }

    # Publication metadata
    if study.pmid:
        publication = {
            "pmid": study.pmid,
            "pmc_id": study.pmc_id,
            "publication_date": study.publication_date,
            "authors": study.authors,
        }

        if study.doi:
            publication["doi"] = format_field_traceable(study.doi)
        if study.publication_title:
            publication["title"] = format_field_traceable(study.publication_title)
        if study.journal:
            publication["journal"] = format_field_traceable(study.journal)
        if study.paper_abstract:
            publication["abstract"] = format_field_traceable(study.paper_abstract)

        output["publication"] = publication

    # Extraction metadata
    output["extraction_metadata"] = {
        "extraction_date": study.extraction_date.isoformat() if study.extraction_date else None,
        "extractor_version": study.extractor_version,
    }

    return output


def create_traceable_report(
    study: StudyMetadata,
    samples: Dict[str, SampleMetadata],
    runs: Optional[List[RunMetadata]] = None,
    include_raw: bool = True,
    include_statistics: bool = True,
) -> Dict[str, Any]:
    """
    Create a complete traceable report with all levels of detail.

    Args:
        study: Study metadata
        samples: Sample metadata dict
        runs: Optional run metadata
        include_raw: Include raw BioSample attributes
        include_statistics: Include extraction statistics

    Returns:
        Hierarchical report with:
        - Quick access views at top level
        - Standard provenance one level down
        - Detailed extraction chain deeper
        - Original raw data at bottom
    """
    report = {
        "metadata_version": "2.0.0",
        "format": "traceable",
        "extraction_timestamp": datetime.now().isoformat(),
        "description": "Hierarchical metadata with complete traceability from raw data to final values",

        # Study metadata
        "study": format_study_traceable(study),

        # Sample metadata
        "samples": {
            sample_id: format_sample_traceable(sample, include_raw=include_raw)
            for sample_id, sample in samples.items()
        },
    }

    # Run metadata (simpler - just full model dump)
    if runs:
        report["runs"] = [run.model_dump() for run in runs]

    # Extraction statistics
    if include_statistics:
        report["extraction_statistics"] = _compute_traceable_statistics(study, samples, runs)

    # Add navigation guide
    report["navigation_guide"] = {
        "description": "How to access different levels of detail",
        "quick_access": "samples[<id>].quick_view.<field> - Just the values",
        "standard_provenance": "samples[<id>].biological_metadata.<field>.provenance - Source + confidence",
        "detailed_debugging": "samples[<id>].biological_metadata.<field>.extraction_details - Full chain",
        "raw_data": "samples[<id>].raw_biosample_attributes - Original BioSample data",
    }

    return report


def _compute_traceable_statistics(
    study: StudyMetadata,
    samples: Dict[str, SampleMetadata],
    runs: Optional[List[RunMetadata]] = None
) -> Dict[str, Any]:
    """Compute extraction statistics."""
    stats = {
        "total_samples": len(samples),
        "total_runs": len(runs) if runs else 0,
    }

    # Field coverage
    field_coverage = {}
    critical_fields = [
        "organism", "tissue", "cell_line", "cell_type", "strain", "genotype", 
        "sex", "age", "developmental_stage", "condition", "treatment", 
        "timepoint", "replicate", "batch", "disease", "stress", 
        "temperature", "growth_condition", "library_strategy"
    ]

    for field in critical_fields:
        count = sum(
            1 for sample in samples.values()
            if getattr(sample, field, None) and getattr(sample, field).value
        )
        field_coverage[field] = {
            "count": count,
            "percentage": round(100 * count / len(samples), 1) if samples else 0,
        }

    stats["field_coverage"] = field_coverage

    # Confidence distribution
    confidence_bins = {"high": 0, "medium": 0, "low": 0}
    confidence_by_field = {}
    total_fields = 0

    for sample in samples.values():
        for field in critical_fields:
            attr = getattr(sample, field, None)
            if attr and isinstance(attr, BaseProvenance):
                total_fields += 1

                # Overall bins
                if attr.confidence >= 0.8:
                    confidence_bins["high"] += 1
                elif attr.confidence >= 0.5:
                    confidence_bins["medium"] += 1
                else:
                    confidence_bins["low"] += 1

                # Per-field tracking
                if field not in confidence_by_field:
                    confidence_by_field[field] = {"high": 0, "medium": 0, "low": 0, "total": 0}

                confidence_by_field[field]["total"] += 1
                if attr.confidence >= 0.8:
                    confidence_by_field[field]["high"] += 1
                elif attr.confidence >= 0.5:
                    confidence_by_field[field]["medium"] += 1
                else:
                    confidence_by_field[field]["low"] += 1

    stats["confidence_distribution"] = {
        "overall": {
            "high (≥0.8)": confidence_bins["high"],
            "medium (0.5-0.8)": confidence_bins["medium"],
            "low (<0.5)": confidence_bins["low"],
            "total_fields": total_fields,
        },
        "by_field": confidence_by_field,
    }

    # Source distribution
    source_counts = {}
    source_by_field = {}

    for sample in samples.values():
        for field in critical_fields:
            attr = getattr(sample, field, None)
            if attr and isinstance(attr, BaseProvenance):
                source = attr.source

                # Overall
                source_counts[source] = source_counts.get(source, 0) + 1

                # Per-field
                if field not in source_by_field:
                    source_by_field[field] = {}
                source_by_field[field][source] = source_by_field[field].get(source, 0) + 1

    stats["source_distribution"] = {
        "overall": source_counts,
        "by_field": source_by_field,
    }

    # Ontology mapping coverage
    ontology_mapped = 0
    ontology_total = 0
    ontology_by_field = {}

    for sample in samples.values():
        for field in ["organism", "tissue", "cell_type"]:
            attr = getattr(sample, field, None)
            if attr and isinstance(attr, BaseProvenance):
                ontology_total += 1

                # Per-field tracking
                if field not in ontology_by_field:
                    ontology_by_field[field] = {"mapped": 0, "total": 0}
                ontology_by_field[field]["total"] += 1

                if attr.ontology_term:
                    ontology_mapped += 1
                    ontology_by_field[field]["mapped"] += 1

    stats["ontology_mapping"] = {
        "overall": {
            "mapped": ontology_mapped,
            "total": ontology_total,
            "percentage": round(100 * ontology_mapped / ontology_total, 1) if ontology_total else 0,
        },
        "by_field": {
            field: {
                **data,
                "percentage": round(100 * data["mapped"] / data["total"], 1) if data["total"] else 0,
            }
            for field, data in ontology_by_field.items()
        },
    }

    # Completeness assessment
    complete_samples = 0
    for sample in samples.values():
        has_organism = bool(sample.organism and sample.organism.value)
        has_tissue_or_cell = bool(
            (sample.tissue and sample.tissue.value) or
            (sample.cell_line and sample.cell_line.value)
        )

        if has_organism and has_tissue_or_cell:
            complete_samples += 1

    stats["completeness"] = {
        "complete_samples": complete_samples,
        "percentage": round(100 * complete_samples / len(samples), 1) if samples else 0,
        "criteria": "organism + (tissue OR cell_line)",
    }

    # Extraction method breakdown
    extraction_methods = {}
    for sample in samples.values():
        for field in critical_fields:
            attr = getattr(sample, field, None)
            if attr and isinstance(attr, BaseProvenance) and attr.extraction_method:
                method = attr.extraction_method
                extraction_methods[method] = extraction_methods.get(method, 0) + 1

    stats["extraction_methods"] = extraction_methods

    return stats
