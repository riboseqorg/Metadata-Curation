#!/usr/bin/env python3
"""
Demo script showing different output format levels.

This demonstrates how to use the MetadataFormatter to create
outputs at different detail levels for end users.
"""

import json
from pathlib import Path
from omics_extractor.extraction.builder import build_project_metadata
from omics_extractor.output.formatters import (
    MetadataFormatter,
    export_to_csv,
    create_provenance_summary,
)


def demo_simple_format(project_id: str = "PRJNA1170270"):
    """
    Simple format - just the values, no provenance.

    Use case: Quick data exploration, spreadsheet users who just want the metadata.
    """
    print("=" * 80)
    print("SIMPLE FORMAT - Values Only")
    print("=" * 80)
    print("\nUse case: Quick exploration, just show me the metadata\n")

    # Extract metadata
    metadata = build_project_metadata(project_id)

    # Format as simple
    report = MetadataFormatter.create_extraction_report(
        study=metadata["study"],
        samples=metadata["samples"],
        runs=metadata["runs"],
        detail_level="simple",
    )

    # Show sample output
    first_sample_id = list(report["samples"].keys())[0]
    first_sample = report["samples"][first_sample_id]

    print("Sample output structure:")
    print(json.dumps(first_sample, indent=2))

    print("\n✓ Simple format includes:")
    print("  - Just field values (strings)")
    print("  - No provenance information")
    print("  - Easy to read, compact")
    print("  - Perfect for quick exploration or CSV export")


def demo_standard_format(project_id: str = "PRJNA1170270"):
    """
    Standard format - values with source and confidence.

    Use case: Most common usage - understand what was extracted and how reliable it is.
    """
    print("\n" + "=" * 80)
    print("STANDARD FORMAT - Values + Provenance Summary")
    print("=" * 80)
    print("\nUse case: Understand where data came from and how confident we are\n")

    # Extract metadata
    metadata = build_project_metadata(project_id)

    # Format as standard
    report = MetadataFormatter.create_extraction_report(
        study=metadata["study"],
        samples=metadata["samples"],
        runs=metadata["runs"],
        detail_level="standard",
    )

    # Show sample output
    first_sample_id = list(report["samples"].keys())[0]
    first_sample = report["samples"][first_sample_id]

    print("Sample output structure:")
    print(json.dumps(first_sample, indent=2))

    print("\n✓ Standard format includes:")
    print("  - Field values")
    print("  - Data source (biosample, sra, llm, etc.)")
    print("  - Confidence score (0.0-1.0)")
    print("  - Ontology term IDs")
    print("  - Perfect for most users - balance of detail and readability")


def demo_detailed_format(project_id: str = "PRJNA1170270"):
    """
    Detailed format - complete extraction chain with all provenance.

    Use case: Debugging, auditing, understanding exactly how metadata was processed.
    """
    print("\n" + "=" * 80)
    print("DETAILED FORMAT - Complete Provenance Chain")
    print("=" * 80)
    print("\nUse case: Debugging, auditing, full traceability of extraction process\n")

    # Extract metadata
    metadata = build_project_metadata(project_id)

    # Format as detailed
    report = MetadataFormatter.create_extraction_report(
        study=metadata["study"],
        samples=metadata["samples"],
        runs=metadata["runs"],
        detail_level="detailed",
    )

    # Show sample output
    first_sample_id = list(report["samples"].keys())[0]
    first_sample = report["samples"][first_sample_id]

    # Show just one field in detail
    if "biological_context" in first_sample and "tissue" in first_sample["biological_context"]:
        print("Example: Tissue field with complete provenance:")
        print(json.dumps(first_sample["biological_context"]["tissue"], indent=2))

    print("\n✓ Detailed format includes:")
    print("  - Field values")
    print("  - Complete confidence breakdown:")
    print("    - Total confidence (field × normalization)")
    print("    - Field confidence (how specific the source field was)")
    print("    - Normalization confidence (how much transformation was needed)")
    print("  - Ontology mapping:")
    print("    - Ontology term ID")
    print("    - Canonical ontology label")
    print("  - Extraction details:")
    print("    - Extraction method (structured_field, llm, ner, etc.)")
    print("    - Timestamp")
    print("    - Original extracted text (before normalization)")
    print("  - Notes/warnings")


def demo_extraction_statistics(project_id: str = "PRJNA1170270"):
    """Show the extraction statistics included in all reports."""
    print("\n" + "=" * 80)
    print("EXTRACTION STATISTICS - Included in All Formats")
    print("=" * 80)
    print("\nProvides quality metrics for the extraction process\n")

    # Extract metadata
    metadata = build_project_metadata(project_id)

    # Format as standard (statistics are same across all formats)
    report = MetadataFormatter.create_extraction_report(
        study=metadata["study"],
        samples=metadata["samples"],
        runs=metadata["runs"],
        detail_level="standard",
    )

    stats = report["extraction_statistics"]

    print("Statistics included:")
    print(f"\n1. Field Coverage:")
    for field, data in stats["field_coverage"].items():
        print(f"   - {field}: {data['count']}/{stats['total_samples']} samples ({data['percentage']}%)")

    print(f"\n2. Confidence Distribution:")
    conf = stats["confidence_distribution"]
    total = conf["total_fields"]
    print(f"   - High confidence (≥0.8): {conf['high (≥0.8)']} / {total}")
    print(f"   - Medium confidence (0.5-0.8): {conf['medium (0.5-0.8)']} / {total}")
    print(f"   - Low confidence (<0.5): {conf['low (<0.5)']} / {total}")

    print(f"\n3. Data Source Distribution:")
    for source, count in stats["source_distribution"].items():
        print(f"   - {source}: {count} fields")

    print(f"\n4. Ontology Mapping:")
    ont = stats["ontology_mapping"]
    print(f"   - Mapped: {ont['mapped']} / {ont['total']} ({ont['percentage']}%)")

    print(f"\n5. Sample Completeness:")
    comp = stats["completeness"]
    print(f"   - Complete samples: {comp['complete_samples']} / {stats['total_samples']} ({comp['percentage']}%)")
    print(f"   - Criteria: {comp['criteria']}")


def demo_csv_export(project_id: str = "PRJNA1170270"):
    """Demonstrate CSV export with and without provenance."""
    print("\n" + "=" * 80)
    print("CSV EXPORT - Spreadsheet-Friendly Format")
    print("=" * 80)

    # Extract metadata
    metadata = build_project_metadata(project_id)

    # Export without provenance (simple)
    simple_csv = Path(f"{project_id}_simple.csv")
    export_to_csv(metadata["samples"], str(simple_csv), include_provenance=False)
    print(f"\n✓ Simple CSV exported to: {simple_csv}")
    print("  Columns: sample_id, organism, tissue, cell_line, strain, etc.")
    print("  Just values, no provenance")

    # Export with provenance (detailed)
    detailed_csv = Path(f"{project_id}_with_provenance.csv")
    export_to_csv(metadata["samples"], str(detailed_csv), include_provenance=True)
    print(f"\n✓ Detailed CSV exported to: {detailed_csv}")
    print("  Columns: same as simple, PLUS:")
    print("  - *_source columns (where each field came from)")
    print("  - *_confidence columns (confidence scores)")
    print("  - *_ontology columns (ontology term IDs)")


def demo_provenance_report(project_id: str = "PRJNA1170270"):
    """Demonstrate markdown provenance summary."""
    print("\n" + "=" * 80)
    print("PROVENANCE SUMMARY REPORT - Human-Readable Overview")
    print("=" * 80)

    # Extract metadata
    metadata = build_project_metadata(project_id)

    # Create summary
    summary = create_provenance_summary(metadata["samples"])

    # Save to file
    report_path = Path(f"{project_id}_provenance.md")
    with open(report_path, "w") as f:
        f.write(summary)

    print(f"\n✓ Provenance report exported to: {report_path}")
    print("\nReport includes:")
    print("  - Field coverage table")
    print("  - Confidence distribution")
    print("  - Data source distribution")
    print("  - Ontology mapping statistics")
    print("  - Sample completeness metrics")

    # Show preview
    print("\nPreview:")
    print("-" * 80)
    print(summary[:500] + "...")


def main():
    """Run all demos."""
    print("\n" + "=" * 80)
    print("METADATA OUTPUT FORMAT DEMONSTRATION")
    print("=" * 80)
    print("\nThis demo shows the different output formats available for metadata extraction.")
    print("Choose the format that best fits your use case:\n")

    try:
        # Run demos
        demo_simple_format()
        demo_standard_format()
        demo_detailed_format()
        demo_extraction_statistics()
        demo_csv_export()
        demo_provenance_report()

        print("\n" + "=" * 80)
        print("SUMMARY")
        print("=" * 80)

        print("\nWhen to use each format:")
        print("\n📊 SIMPLE FORMAT")
        print("  - Quick data exploration")
        print("  - Spreadsheet analysis")
        print("  - When you just need the values")

        print("\n📈 STANDARD FORMAT (Recommended)")
        print("  - Most common use case")
        print("  - Understand data quality")
        print("  - Balance of detail and readability")
        print("  - Good for downstream analysis with QC")

        print("\n🔍 DETAILED FORMAT")
        print("  - Debugging extraction issues")
        print("  - Auditing data provenance")
        print("  - Understanding transformation chain")
        print("  - Full traceability for compliance")

        print("\n📄 CSV EXPORT")
        print("  - Spreadsheet tools (Excel, Google Sheets)")
        print("  - R/pandas analysis")
        print("  - Quick filtering and sorting")

        print("\n📝 PROVENANCE REPORT")
        print("  - High-level quality assessment")
        print("  - Sharing extraction results with team")
        print("  - Documentation for publications")

        print("\n" + "=" * 80)

    except Exception as e:
        print(f"\n✗ Error: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()
