"""Output formatters for metadata with various detail levels.

Provides multiple output formats to help end users understand:
1. What was extracted
2. Where it came from (provenance)
3. How confident we are
4. What transformations were applied
"""

from typing import Dict, Any, List, Optional
from datetime import datetime
from ..schemas.base import BaseProvenance, SampleMetadata, StudyMetadata, RunMetadata


class MetadataFormatter:
    """Format metadata with various levels of detail."""

    @staticmethod
    def format_provenance_simple(prov: BaseProvenance) -> str:
        """Just the value - simplest output."""
        return prov.value

    @staticmethod
    def format_provenance_standard(prov: BaseProvenance) -> Dict[str, Any]:
        """Standard output with key provenance info."""
        return {
            "value": prov.value,
            "source": prov.source,
            "confidence": round(prov.confidence, 2),
            "ontology_term": prov.ontology_term,
        }

    @staticmethod
    def format_provenance_detailed(prov: BaseProvenance) -> Dict[str, Any]:
        """Detailed output with full provenance chain."""
        return {
            "value": prov.value,
            "source": prov.source,
            "source_id": prov.source_id,
            "confidence": {
                "total": round(prov.confidence, 3),
                "field": round(prov.field_confidence, 3) if prov.field_confidence else None,
                "normalization": round(prov.normalization_confidence, 3) if prov.normalization_confidence else None,
            },
            "ontology": {
                "term": prov.ontology_term,
                "label": prov.ontology_label,
            } if prov.ontology_term else None,
            "extraction": {
                "method": prov.extraction_method,
                "timestamp": prov.extraction_timestamp.isoformat() if prov.extraction_timestamp else None,
                "original_text": prov.extracted_text,
            },
            "notes": prov.notes,
        }

    @staticmethod
    def format_sample_simple(sample: SampleMetadata) -> Dict[str, str]:
        """Simple format - just values, no provenance."""
        output = {
            "sample_id": sample.sample_id,
            "organism": sample.organism.value if sample.organism else None,
        }

        # Add optional fields if present
        optional_fields = [
            "tissue", "cell_line", "cell_type", "developmental_stage",
            "strain", "genotype", "sex", "age", "condition", "treatment",
            "timepoint", "replicate", "disease", "sample_title"
        ]

        for field in optional_fields:
            attr = getattr(sample, field, None)
            if attr and isinstance(attr, BaseProvenance):
                output[field] = attr.value
            elif attr:
                output[field] = attr

        return output

    @staticmethod
    def format_sample_standard(sample: SampleMetadata) -> Dict[str, Any]:
        """Standard format - values with source and confidence."""
        output = {
            "sample_id": sample.sample_id,
            "bioproject_id": sample.bioproject_id,
            "biosample_id": sample.biosample_id,
            "organism": MetadataFormatter.format_provenance_standard(sample.organism),
        }

        # Add optional fields with provenance
        optional_fields = [
            "tissue", "cell_line", "cell_type", "developmental_stage",
            "strain", "genotype", "sex", "age", "condition", "treatment",
            "timepoint", "replicate", "disease", "sample_title", "sample_description"
        ]

        for field in optional_fields:
            attr = getattr(sample, field, None)
            if attr and isinstance(attr, BaseProvenance):
                output[field] = MetadataFormatter.format_provenance_standard(attr)

        return output

    @staticmethod
    def format_sample_detailed(sample: SampleMetadata) -> Dict[str, Any]:
        """Detailed format - complete provenance chain."""
        output = {
            "identifiers": {
                "sample_id": sample.sample_id,
                "bioproject_id": sample.bioproject_id,
                "biosample_id": sample.biosample_id,
                "gsm_id": sample.gsm_id,
            },
            "organism": MetadataFormatter.format_provenance_detailed(sample.organism),
        }

        # Biological context
        biological_fields = {
            "tissue": sample.tissue,
            "cell_line": sample.cell_line,
            "cell_type": sample.cell_type,
            "developmental_stage": sample.developmental_stage,
            "strain": sample.strain,
            "genotype": sample.genotype,
            "sex": sample.sex,
            "age": sample.age,
        }

        output["biological_context"] = {}
        for field, attr in biological_fields.items():
            if attr and isinstance(attr, BaseProvenance):
                output["biological_context"][field] = MetadataFormatter.format_provenance_detailed(attr)

        # Experimental conditions
        experimental_fields = {
            "condition": sample.condition,
            "treatment": sample.treatment,
            "timepoint": sample.timepoint,
            "replicate": sample.replicate,
            "batch": sample.batch,
        }

        output["experimental_conditions"] = {}
        for field, attr in experimental_fields.items():
            if attr and isinstance(attr, BaseProvenance):
                output["experimental_conditions"][field] = MetadataFormatter.format_provenance_detailed(attr)

        # Disease/perturbation
        if sample.disease or sample.stress or sample.temperature or sample.growth_condition:
            output["disease_perturbation"] = {}
            for field in ["disease", "stress", "temperature", "growth_condition"]:
                attr = getattr(sample, field, None)
                if attr and isinstance(attr, BaseProvenance):
                    output["disease_perturbation"][field] = MetadataFormatter.format_provenance_detailed(attr)

        # Sample descriptions
        if sample.sample_title or sample.sample_description:
            output["descriptions"] = {}
            if sample.sample_title:
                output["descriptions"]["title"] = MetadataFormatter.format_provenance_detailed(sample.sample_title)
            if sample.sample_description:
                output["descriptions"]["description"] = MetadataFormatter.format_provenance_detailed(sample.sample_description)

        # Custom fields
        if sample.custom_fields:
            output["custom_fields"] = sample.custom_fields

        return output

    @staticmethod
    def format_study_simple(study: StudyMetadata) -> Dict[str, str]:
        """Simple study format - just values."""
        return {
            "bioproject_id": study.bioproject_id,
            "title": study.title.value,
            "organism": study.organism.value,
            "pmid": study.pmid,
            "gse_id": study.gse_id,
        }

    @staticmethod
    def format_study_standard(study: StudyMetadata) -> Dict[str, Any]:
        """Standard study format with provenance."""
        output = {
            "bioproject_id": study.bioproject_id,
            "title": MetadataFormatter.format_provenance_standard(study.title),
            "description": MetadataFormatter.format_provenance_standard(study.description),
            "organism": MetadataFormatter.format_provenance_standard(study.organism),
        }

        # Optional IDs
        if study.gse_id:
            output["gse_id"] = study.gse_id
        if study.sra_study_id:
            output["sra_study_id"] = study.sra_study_id

        # Publication info
        if study.pmid:
            output["publication"] = {
                "pmid": study.pmid,
                "pmc_id": study.pmc_id,
                "doi": MetadataFormatter.format_provenance_standard(study.doi) if study.doi else None,
                "title": MetadataFormatter.format_provenance_standard(study.publication_title) if study.publication_title else None,
                "journal": MetadataFormatter.format_provenance_standard(study.journal) if study.journal else None,
                "publication_date": study.publication_date,
            }

        return output

    @staticmethod
    def format_study_detailed(study: StudyMetadata) -> Dict[str, Any]:
        """Detailed study format with complete provenance."""
        output = {
            "identifiers": {
                "bioproject_id": study.bioproject_id,
                "gse_id": study.gse_id,
                "sra_study_id": study.sra_study_id,
            },
            "core_metadata": {
                "title": MetadataFormatter.format_provenance_detailed(study.title),
                "description": MetadataFormatter.format_provenance_detailed(study.description),
                "organism": MetadataFormatter.format_provenance_detailed(study.organism),
            },
        }

        # Publication metadata
        if study.pmid:
            output["publication"] = {
                "pmid": study.pmid,
                "pmc_id": study.pmc_id,
                "doi": MetadataFormatter.format_provenance_detailed(study.doi) if study.doi else None,
                "title": MetadataFormatter.format_provenance_detailed(study.publication_title) if study.publication_title else None,
                "journal": MetadataFormatter.format_provenance_detailed(study.journal) if study.journal else None,
                "authors": study.authors,
                "publication_date": study.publication_date,
                "abstract": MetadataFormatter.format_provenance_detailed(study.paper_abstract) if study.paper_abstract else None,
            }

        # Extraction metadata
        output["extraction_metadata"] = {
            "extraction_date": study.extraction_date.isoformat() if study.extraction_date else None,
            "extractor_version": study.extractor_version,
        }

        # Custom fields
        if study.custom_fields:
            output["custom_fields"] = study.custom_fields

        return output

    @staticmethod
    def create_extraction_report(
        study: StudyMetadata,
        samples: Dict[str, SampleMetadata],
        runs: Optional[List[RunMetadata]] = None,
        detail_level: str = "standard"
    ) -> Dict[str, Any]:
        """
        Create a comprehensive extraction report with provenance tracking.

        Args:
            study: Study metadata
            samples: Sample metadata dict
            runs: Optional run metadata
            detail_level: "simple", "standard", or "detailed"

        Returns:
            Complete report with metadata and extraction statistics
        """
        # Select formatter based on detail level
        if detail_level == "simple":
            study_formatter = MetadataFormatter.format_study_simple
            sample_formatter = MetadataFormatter.format_sample_simple
        elif detail_level == "detailed":
            study_formatter = MetadataFormatter.format_study_detailed
            sample_formatter = MetadataFormatter.format_sample_detailed
        else:  # standard
            study_formatter = MetadataFormatter.format_study_standard
            sample_formatter = MetadataFormatter.format_sample_standard

        # Build report
        report = {
            "metadata_version": "1.0.0",
            "extraction_timestamp": datetime.now().isoformat(),
            "detail_level": detail_level,
            "study": study_formatter(study),
            "samples": {
                sample_id: sample_formatter(sample)
                for sample_id, sample in samples.items()
            },
        }

        # Add runs if provided
        if runs:
            report["runs"] = [run.model_dump() for run in runs]

        # Compute extraction statistics
        report["extraction_statistics"] = MetadataFormatter._compute_statistics(
            study, samples, runs
        )

        return report

    @staticmethod
    def _compute_statistics(
        study: StudyMetadata,
        samples: Dict[str, SampleMetadata],
        runs: Optional[List[RunMetadata]] = None
    ) -> Dict[str, Any]:
        """Compute statistics about extraction quality and coverage."""
        stats = {
            "total_samples": len(samples),
            "total_runs": len(runs) if runs else 0,
        }

        # Field coverage
        field_coverage = {}
        critical_fields = ["organism", "tissue", "cell_line", "cell_type", "strain"]

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
        total_fields = 0

        for sample in samples.values():
            for field in critical_fields:
                attr = getattr(sample, field, None)
                if attr and isinstance(attr, BaseProvenance):
                    total_fields += 1
                    if attr.confidence >= 0.8:
                        confidence_bins["high"] += 1
                    elif attr.confidence >= 0.5:
                        confidence_bins["medium"] += 1
                    else:
                        confidence_bins["low"] += 1

        stats["confidence_distribution"] = {
            "high (≥0.8)": confidence_bins["high"],
            "medium (0.5-0.8)": confidence_bins["medium"],
            "low (<0.5)": confidence_bins["low"],
            "total_fields": total_fields,
        }

        # Source distribution
        source_counts = {}
        for sample in samples.values():
            for field in critical_fields:
                attr = getattr(sample, field, None)
                if attr and isinstance(attr, BaseProvenance):
                    source = attr.source
                    source_counts[source] = source_counts.get(source, 0) + 1

        stats["source_distribution"] = source_counts

        # Ontology mapping coverage
        ontology_mapped = 0
        ontology_total = 0

        for sample in samples.values():
            for field in ["organism", "tissue", "cell_type"]:
                attr = getattr(sample, field, None)
                if attr and isinstance(attr, BaseProvenance):
                    ontology_total += 1
                    if attr.ontology_term:
                        ontology_mapped += 1

        stats["ontology_mapping"] = {
            "mapped": ontology_mapped,
            "total": ontology_total,
            "percentage": round(100 * ontology_mapped / ontology_total, 1) if ontology_total else 0,
        }

        # Completeness assessment
        complete_samples = 0
        for sample in samples.values():
            has_organism = bool(sample.organism and sample.organism.value)
            has_tissue_or_cell = bool(
                (sample.tissue and sample.tissue.value) or
                (sample.cell_line and sample.cell_line.value)
            )
            has_strain = bool(sample.strain and sample.strain.value)

            if has_organism and has_tissue_or_cell:
                complete_samples += 1

        stats["completeness"] = {
            "complete_samples": complete_samples,
            "percentage": round(100 * complete_samples / len(samples), 1) if samples else 0,
            "criteria": "organism + (tissue OR cell_line)",
        }

        return stats


def export_to_csv(samples: Dict[str, SampleMetadata], output_path: str, include_provenance: bool = False):
    """
    Export samples to CSV format.

    Args:
        samples: Sample metadata dict
        output_path: Output CSV file path
        include_provenance: If True, add source/confidence columns
    """
    import csv

    # Define columns
    base_columns = [
        "sample_id", "bioproject_id", "biosample_id",
        "organism", "tissue", "cell_line", "cell_type",
        "strain", "genotype", "sex", "age",
        "developmental_stage", "treatment", "condition"
    ]

    if include_provenance:
        # Add _source and _confidence columns for each field
        prov_columns = []
        for col in base_columns[3:]:  # Skip IDs
            prov_columns.extend([f"{col}_source", f"{col}_confidence", f"{col}_ontology"])
        columns = base_columns + prov_columns
    else:
        columns = base_columns

    # Write CSV
    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=columns)
        writer.writeheader()

        for sample in samples.values():
            row = {}

            # Base values
            for field in base_columns:
                if field in ["sample_id", "bioproject_id", "biosample_id"]:
                    row[field] = getattr(sample, field, "")
                else:
                    attr = getattr(sample, field, None)
                    row[field] = attr.value if attr and isinstance(attr, BaseProvenance) else ""

            # Provenance if requested
            if include_provenance:
                for field in base_columns[3:]:
                    attr = getattr(sample, field, None)
                    if attr and isinstance(attr, BaseProvenance):
                        row[f"{field}_source"] = attr.source
                        row[f"{field}_confidence"] = round(attr.confidence, 3)
                        row[f"{field}_ontology"] = attr.ontology_term or ""
                    else:
                        row[f"{field}_source"] = ""
                        row[f"{field}_confidence"] = ""
                        row[f"{field}_ontology"] = ""

            writer.writerow(row)


def create_provenance_summary(samples: Dict[str, SampleMetadata]) -> str:
    """
    Create a human-readable provenance summary.

    Args:
        samples: Sample metadata dict

    Returns:
        Markdown-formatted summary
    """
    lines = ["# Metadata Extraction Provenance Summary\n"]

    # Overall statistics
    stats = MetadataFormatter._compute_statistics(None, samples, None)
    lines.append(f"**Total Samples:** {stats['total_samples']}\n")

    # Field coverage
    lines.append("## Field Coverage\n")
    lines.append("| Field | Samples | Coverage |")
    lines.append("|-------|---------|----------|")
    for field, data in stats["field_coverage"].items():
        lines.append(f"| {field} | {data['count']} | {data['percentage']}% |")
    lines.append("")

    # Confidence distribution
    lines.append("## Confidence Distribution\n")
    conf_dist = stats["confidence_distribution"]
    total = conf_dist["total_fields"]
    lines.append(f"- **High confidence (≥0.8):** {conf_dist['high (≥0.8)']} / {total} ({100*conf_dist['high (≥0.8)']/total:.1f}%)")
    lines.append(f"- **Medium confidence (0.5-0.8):** {conf_dist['medium (0.5-0.8)']} / {total} ({100*conf_dist['medium (0.5-0.8)']/total:.1f}%)")
    lines.append(f"- **Low confidence (<0.5):** {conf_dist['low (<0.5)']} / {total} ({100*conf_dist['low (<0.5)']/total:.1f}%)\n")

    # Source distribution
    lines.append("## Data Sources\n")
    for source, count in sorted(stats["source_distribution"].items(), key=lambda x: x[1], reverse=True):
        lines.append(f"- **{source}:** {count} fields")
    lines.append("")

    # Ontology mapping
    ont_stats = stats["ontology_mapping"]
    lines.append("## Ontology Mapping\n")
    lines.append(f"**{ont_stats['mapped']} / {ont_stats['total']} fields mapped ({ont_stats['percentage']}%)**\n")

    # Completeness
    comp = stats["completeness"]
    lines.append("## Sample Completeness\n")
    lines.append(f"**{comp['complete_samples']} / {stats['total_samples']} samples complete ({comp['percentage']}%)**")
    lines.append(f"_Criteria: {comp['criteria']}_\n")

    return "\n".join(lines)
