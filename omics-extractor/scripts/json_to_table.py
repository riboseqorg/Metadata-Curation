#!/usr/bin/env python3
"""
Convert enriched JSON files to tabular format (CSV/TSV).

For each metadata field, creates 3 columns:
  - field: original value
  - field_enriched: LLM-enriched value (if added)
  - field_ontology: ontology term (if available)

Usage:
    # Create TSV (default)
    python scripts/json_to_table.py enriched/*.json -o output.tsv

    # Create CSV
    python scripts/json_to_table.py enriched/*.json -o output.csv --format csv

    # Include all fields
    python scripts/json_to_table.py enriched/*.json -o output.tsv --all-fields
"""

import json
import csv
import sys
from pathlib import Path
from typing import List, Dict, Any, Optional
import argparse


def extract_field_value(field_obj: Any) -> tuple[Optional[str], Optional[str], Optional[str]]:
    """
    Extract value, source, and ontology from a field.

    Returns:
        (value, source, ontology_term)
    """
    if field_obj is None:
        return None, None, None

    if isinstance(field_obj, dict):
        value = field_obj.get("value")
        source = field_obj.get("source", "")
        ontology = field_obj.get("ontology_term")
        return value, source, ontology
    else:
        # Simple string value
        return str(field_obj), None, None


def process_sample(sample_data: dict, study_id: str, format_type: str) -> Dict[str, Any]:
    """
    Process a single sample and extract all fields.

    Returns:
        Dictionary with flattened fields including _enriched and _ontology variants
    """
    row = {}

    # Basic identifiers
    row["study_id"] = study_id

    # Detect format
    is_traceable = "biological_metadata" in sample_data

    if is_traceable:
        # Traceable format
        quick_view = sample_data.get("quick_view", {})
        bio_meta = sample_data.get("biological_metadata", {})

        row["sample_id"] = quick_view.get("sample_id")
        row["bioproject_id"] = quick_view.get("bioproject_id")

        # Extract all metadata fields (comprehensive list from SampleMetadata schema)
        fields_to_extract = [
            # Core biological metadata
            "organism", "tissue", "cell_type", "cell_line", "strain",
            "developmental_stage", "genotype", "age", "sex",
            # Experimental conditions
            "condition", "treatment", "timepoint", "replicate", "batch",
            # Disease/perturbation
            "disease", "stress", "temperature", "growth_condition"
        ]

        for field in fields_to_extract:
            if field in bio_meta:
                value, source, ontology = extract_field_value(bio_meta[field])

                # Base value
                row[field] = value

                # Mark if enriched
                row[f"{field}_enriched"] = "yes" if source == "llm_enrichment" else "no"

                # Ontology term
                row[f"{field}_ontology"] = ontology

        # Sample title/description
        sample_title = sample_data.get("sample_title", {})
        sample_desc = sample_data.get("sample_description", {})

        row["sample_title"] = sample_title.get("value") if isinstance(sample_title, dict) else sample_title
        row["sample_description"] = sample_desc.get("value") if isinstance(sample_desc, dict) else sample_desc

    else:
        # Flat format
        row["sample_id"] = sample_data.get("sample_id")
        row["bioproject_id"] = sample_data.get("bioproject_id")

        fields_to_extract = [
            # Core biological metadata
            "organism", "tissue", "cell_type", "cell_line", "strain",
            "developmental_stage", "genotype", "age", "sex",
            # Experimental conditions
            "condition", "treatment", "timepoint", "replicate", "batch",
            # Disease/perturbation
            "disease", "stress", "temperature", "growth_condition"
        ]

        for field in fields_to_extract:
            if field in sample_data:
                value, source, ontology = extract_field_value(sample_data.get(field))

                row[field] = value
                row[f"{field}_enriched"] = "yes" if source == "llm_enrichment" else "no"
                row[f"{field}_ontology"] = ontology

        # Sample title/description
        sample_title = sample_data.get("sample_title")
        sample_desc = sample_data.get("sample_description")

        if isinstance(sample_title, dict):
            row["sample_title"] = sample_title.get("value")
        else:
            row["sample_title"] = sample_title

        if isinstance(sample_desc, dict):
            row["sample_description"] = sample_desc.get("value")
        else:
            row["sample_description"] = sample_desc

    return row


def process_file(file_path: Path) -> List[Dict[str, Any]]:
    """Process a single JSON file and return list of row dicts."""

    with open(file_path) as f:
        data = json.load(f)

    study_id = file_path.stem.replace("_enriched", "")
    samples = data.get("samples", {})

    rows = []
    for sample_id, sample_data in samples.items():
        try:
            row = process_sample(sample_data, study_id, "traceable")
            rows.append(row)
        except Exception as e:
            print(f"Warning: Failed to process {study_id}/{sample_id}: {e}", file=sys.stderr)
            continue

    return rows


def write_table(rows: List[Dict[str, Any]], output_path: Path, delimiter: str = "\t"):
    """Write rows to CSV/TSV file."""

    if not rows:
        print("No rows to write!", file=sys.stderr)
        return

    # Determine all columns (union of all keys)
    all_columns = set()
    for row in rows:
        all_columns.update(row.keys())

    # Order columns logically
    ordered_columns = []

    # ID columns first
    for col in ["study_id", "bioproject_id", "sample_id"]:
        if col in all_columns:
            ordered_columns.append(col)
            all_columns.remove(col)

    # Then metadata fields (grouped with their _enriched and _ontology variants)
    base_fields = [
        # Core biological metadata
        "organism", "tissue", "cell_type", "cell_line", "strain",
        "developmental_stage", "genotype", "age", "sex",
        # Experimental conditions
        "condition", "treatment", "timepoint", "replicate", "batch",
        # Disease/perturbation
        "disease", "stress", "temperature", "growth_condition"
    ]

    for field in base_fields:
        if field in all_columns:
            ordered_columns.append(field)
            all_columns.remove(field)

        enriched_col = f"{field}_enriched"
        if enriched_col in all_columns:
            ordered_columns.append(enriched_col)
            all_columns.remove(enriched_col)

        ontology_col = f"{field}_ontology"
        if ontology_col in all_columns:
            ordered_columns.append(ontology_col)
            all_columns.remove(ontology_col)

    # Sample title/description
    for col in ["sample_title", "sample_description"]:
        if col in all_columns:
            ordered_columns.append(col)
            all_columns.remove(col)

    # Remaining columns
    ordered_columns.extend(sorted(all_columns))

    # Write to file
    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=ordered_columns, delimiter=delimiter, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

    print(f"✓ Wrote {len(rows)} rows to {output_path}")
    print(f"  Columns: {len(ordered_columns)}")
    print(f"  Studies: {len(set(row.get('study_id') for row in rows))}")


def main():
    parser = argparse.ArgumentParser(
        description="Convert enriched JSON files to tabular format",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Create TSV (default)
  %(prog)s enriched/*.json -o output.tsv

  # Create CSV
  %(prog)s enriched/*.json -o output.csv --format csv

  # Process specific files
  %(prog)s enriched/PRJNA*.json -o projects.tsv
        """
    )
    parser.add_argument("files", nargs="+", help="Enriched JSON files to process")
    parser.add_argument("--output", "-o", required=True, help="Output file (CSV or TSV)")
    parser.add_argument("--format", choices=["csv", "tsv"], default="tsv",
                       help="Output format (default: tsv)")

    args = parser.parse_args()

    # Determine delimiter
    delimiter = "," if args.format == "csv" else "\t"

    # Also check file extension
    output_path = Path(args.output)
    if output_path.suffix.lower() == ".csv":
        delimiter = ","

    # Process all files
    all_rows = []
    for file_path in args.files:
        path = Path(file_path)
        if not path.exists():
            print(f"Warning: {file_path} not found", file=sys.stderr)
            continue

        print(f"Processing {path.name}...", file=sys.stderr)
        rows = process_file(path)
        all_rows.extend(rows)

    if not all_rows:
        print("No data to write!", file=sys.stderr)
        sys.exit(1)

    # Write output
    write_table(all_rows, output_path, delimiter)


if __name__ == "__main__":
    main()
