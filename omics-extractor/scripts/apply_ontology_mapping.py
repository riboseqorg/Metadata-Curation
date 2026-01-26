#!/usr/bin/env python3
"""
Apply ontology normalization to enriched JSON files.

This script:
1. Reads enriched metadata JSON files
2. Maps tissue/cell_type/organism values to standardized ontology terms
3. Adds normalized values and ontology IDs
4. Reports mapping statistics and unmapped terms

Usage:
    # Process files and update them in-place
    python scripts/apply_ontology_mapping.py enriched/*.json --in-place

    # Create new normalized files
    python scripts/apply_ontology_mapping.py enriched/*.json -o normalized/

    # Dry run (show what would be mapped)
    python scripts/apply_ontology_mapping.py enriched/*.json --dry-run

    # Show unmapped terms only
    python scripts/apply_ontology_mapping.py enriched/*.json --show-unmapped
"""

import json
import sys
from pathlib import Path
from typing import Dict, List, Set, Tuple
from collections import defaultdict
import argparse

# Add parent directory to path to import omics_extractor
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from omics_extractor.ontologies.mapper import normalize_tissue, normalize_cell_type, normalize_organism


def apply_ontology_mapping(sample_data: dict, field: str, normalizer_func) -> Tuple[bool, float]:
    """
    Apply ontology normalization to a single field in a sample.

    Args:
        sample_data: Sample data dict
        field: Field name (e.g., 'tissue', 'cell_type', 'organism')
        normalizer_func: Function to normalize the field (e.g., normalize_tissue)

    Returns:
        (was_mapped, confidence_score)
    """
    # Determine format
    is_traceable = "biological_metadata" in sample_data

    if is_traceable:
        bio_meta = sample_data.get("biological_metadata", {})
        if field not in bio_meta:
            return False, 0.0

        field_obj = bio_meta[field]
    else:
        # Flat format
        if field not in sample_data:
            return False, 0.0
        field_obj = sample_data[field]

    # Extract value
    if isinstance(field_obj, dict):
        raw_value = field_obj.get("value")
    else:
        raw_value = field_obj

    if not raw_value or not isinstance(raw_value, str):
        return False, 0.0

    # Apply normalization
    normalized_value, ontology_id, confidence = normalizer_func(raw_value)

    # Only update if we got a mapping (ontology_id is not None)
    if ontology_id:
        # Update the field object
        if isinstance(field_obj, dict):
            # Already a dict - update it
            field_obj["value"] = normalized_value
            field_obj["ontology_term"] = ontology_id
            # Update confidence if this is higher than existing
            existing_conf = field_obj.get("confidence", 0.0)
            if confidence > existing_conf:
                field_obj["confidence"] = confidence
        else:
            # Convert to dict with ontology info
            new_field_obj = {
                "value": normalized_value,
                "source": "ontology_mapping",
                "ontology_term": ontology_id,
                "confidence": confidence,
            }

            if is_traceable:
                bio_meta[field] = new_field_obj
            else:
                sample_data[field] = new_field_obj

        return True, confidence

    return False, 0.0


def process_file(file_path: Path, dry_run: bool = False) -> Dict:
    """
    Process a single enriched JSON file and apply ontology mappings.

    Returns:
        Statistics dict
    """
    with open(file_path) as f:
        data = json.load(f)

    study_id = file_path.stem.replace("_enriched", "")
    samples = data.get("samples", {})

    stats = {
        "study_id": study_id,
        "total_samples": len(samples),
        "fields_mapped": defaultdict(int),
        "unmapped_values": defaultdict(set),
        "avg_confidence": defaultdict(list),
    }

    # Apply mappings to each sample
    for sample_id, sample_data in samples.items():
        # Tissue
        was_mapped, conf = apply_ontology_mapping(sample_data, "tissue", normalize_tissue)
        if was_mapped:
            stats["fields_mapped"]["tissue"] += 1
            stats["avg_confidence"]["tissue"].append(conf)
        else:
            # Track unmapped value
            is_traceable = "biological_metadata" in sample_data
            if is_traceable:
                bio_meta = sample_data.get("biological_metadata", {})
                if "tissue" in bio_meta:
                    val = bio_meta["tissue"].get("value") if isinstance(bio_meta["tissue"], dict) else bio_meta["tissue"]
                    if val:
                        stats["unmapped_values"]["tissue"].add(str(val))
            else:
                if "tissue" in sample_data:
                    val = sample_data["tissue"].get("value") if isinstance(sample_data["tissue"], dict) else sample_data["tissue"]
                    if val:
                        stats["unmapped_values"]["tissue"].add(str(val))

        # Cell type
        was_mapped, conf = apply_ontology_mapping(sample_data, "cell_type", normalize_cell_type)
        if was_mapped:
            stats["fields_mapped"]["cell_type"] += 1
            stats["avg_confidence"]["cell_type"].append(conf)
        else:
            is_traceable = "biological_metadata" in sample_data
            if is_traceable:
                bio_meta = sample_data.get("biological_metadata", {})
                if "cell_type" in bio_meta:
                    val = bio_meta["cell_type"].get("value") if isinstance(bio_meta["cell_type"], dict) else bio_meta["cell_type"]
                    if val:
                        stats["unmapped_values"]["cell_type"].add(str(val))
            else:
                if "cell_type" in sample_data:
                    val = sample_data["cell_type"].get("value") if isinstance(sample_data["cell_type"], dict) else sample_data["cell_type"]
                    if val:
                        stats["unmapped_values"]["cell_type"].add(str(val))

        # Organism
        was_mapped, conf = apply_ontology_mapping(sample_data, "organism", normalize_organism)
        if was_mapped:
            stats["fields_mapped"]["organism"] += 1
            stats["avg_confidence"]["organism"].append(conf)
        else:
            is_traceable = "biological_metadata" in sample_data
            if is_traceable:
                bio_meta = sample_data.get("biological_metadata", {})
                if "organism" in bio_meta:
                    val = bio_meta["organism"].get("value") if isinstance(bio_meta["organism"], dict) else bio_meta["organism"]
                    if val:
                        stats["unmapped_values"]["organism"].add(str(val))
            else:
                if "organism" in sample_data:
                    val = sample_data["organism"].get("value") if isinstance(sample_data["organism"], dict) else sample_data["organism"]
                    if val:
                        stats["unmapped_values"]["organism"].add(str(val))

    return data, stats


def print_summary(all_stats: List[Dict], show_unmapped: bool = False):
    """Print summary statistics."""

    print("\n" + "=" * 80)
    print("ONTOLOGY MAPPING SUMMARY")
    print("=" * 80)

    total_studies = len(all_stats)
    total_samples = sum(s["total_samples"] for s in all_stats)

    print(f"\nProcessed {total_studies} studies, {total_samples} samples")

    # Aggregate field stats
    print(f"\n{'Field':<15} {'Mapped':<10} {'Avg Confidence':<20} {'Unmapped Terms':<15}")
    print("-" * 70)

    for field in ["tissue", "cell_type", "organism"]:
        total_mapped = sum(s["fields_mapped"][field] for s in all_stats)

        all_confidences = []
        for s in all_stats:
            all_confidences.extend(s["avg_confidence"][field])

        avg_conf = sum(all_confidences) / len(all_confidences) if all_confidences else 0.0

        # Collect all unmapped values
        unmapped_set = set()
        for s in all_stats:
            unmapped_set.update(s["unmapped_values"][field])

        num_unmapped = len(unmapped_set)

        print(f"{field:<15} {total_mapped:<10} {avg_conf:.3f} ({len(all_confidences)})      {num_unmapped:<15}")

        if show_unmapped and unmapped_set:
            print(f"\n  Unmapped {field} terms:")
            for term in sorted(unmapped_set):
                print(f"    - {term}")
            print()

    print("\nMapping Quality:")
    print("  1.0   = Exact ontology match")
    print("  0.95  = Case-insensitive match")
    print("  0.9   = Synonym match")
    print("  0.7   = Fuzzy match")
    print("  0.5   = No match (original value kept)")


def main():
    parser = argparse.ArgumentParser(
        description="Apply ontology normalization to enriched metadata",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Update files in-place
  %(prog)s enriched/*.json --in-place

  # Create normalized copies
  %(prog)s enriched/*.json -o normalized/

  # Dry run (show what would be mapped)
  %(prog)s enriched/*.json --dry-run

  # Show unmapped terms
  %(prog)s enriched/*.json --show-unmapped
        """
    )
    parser.add_argument("files", nargs="+", help="Enriched JSON files to process")
    parser.add_argument("--in-place", action="store_true", help="Update files in-place")
    parser.add_argument("--output", "-o", help="Output directory for normalized files")
    parser.add_argument("--dry-run", action="store_true", help="Show stats without modifying files")
    parser.add_argument("--show-unmapped", action="store_true", help="Show unmapped terms")

    args = parser.parse_args()

    if args.in_place and args.output:
        print("Error: Cannot use both --in-place and --output", file=sys.stderr)
        sys.exit(1)

    if args.output:
        output_dir = Path(args.output)
        output_dir.mkdir(parents=True, exist_ok=True)

    all_stats = []

    for file_path in args.files:
        path = Path(file_path)
        if not path.exists():
            print(f"Warning: {file_path} not found", file=sys.stderr)
            continue

        print(f"Processing {path.name}...", file=sys.stderr)

        data, stats = process_file(path, dry_run=args.dry_run)
        all_stats.append(stats)

        # Write output if not dry run
        if not args.dry_run:
            if args.in_place:
                with open(path, "w") as f:
                    json.dump(data, f, indent=2)
                print(f"  Updated {path}", file=sys.stderr)
            elif args.output:
                output_path = output_dir / path.name
                with open(output_path, "w") as f:
                    json.dump(data, f, indent=2)
                print(f"  Wrote {output_path}", file=sys.stderr)

    if not all_stats:
        print("No files processed!", file=sys.stderr)
        sys.exit(1)

    # Print summary
    print_summary(all_stats, show_unmapped=args.show_unmapped)


if __name__ == "__main__":
    main()
