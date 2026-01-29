#!/usr/bin/env python3
"""
Assess enrichment performance from batch-enriched JSON files.

This script analyzes enriched metadata to provide:
1. Field coverage statistics (% of samples with each field)
2. Enrichment quality metrics (confidence scores)
3. Before/after comparison
4. Per-study breakdown

Usage:
    python scripts/assess_enrichment.py enriched/*.json
    python scripts/assess_enrichment.py enriched/*.json --verbose
    python scripts/assess_enrichment.py enriched/*.json --output report.json
"""

import json
import sys
from pathlib import Path
from typing import Dict, List, Optional
from collections import defaultdict
from dataclasses import dataclass, asdict


@dataclass
class FieldStats:
    """Statistics for a single field."""

    total_samples: int = 0
    filled_before: int = 0
    filled_after: int = 0
    enriched_count: int = 0  # Newly added by LLM
    avg_confidence: float = 0.0
    confidence_scores: List[float] = None

    def __post_init__(self):
        if self.confidence_scores is None:
            self.confidence_scores = []

    @property
    def coverage_before(self) -> float:
        """% coverage before enrichment."""
        return (self.filled_before / self.total_samples * 100) if self.total_samples > 0 else 0.0

    @property
    def coverage_after(self) -> float:
        """% coverage after enrichment."""
        return (self.filled_after / self.total_samples * 100) if self.total_samples > 0 else 0.0

    @property
    def improvement(self) -> float:
        """Percentage point improvement."""
        return self.coverage_after - self.coverage_before


@dataclass
class StudyStats:
    """Statistics for a single study."""

    study_id: str
    num_samples: int = 0
    fields_enriched: int = 0
    avg_confidence: float = 0.0
    field_stats: Dict[str, FieldStats] = None

    def __post_init__(self):
        if self.field_stats is None:
            self.field_stats = {}


def analyze_enriched_file(file_path: Path, verbose: bool = False) -> StudyStats:
    """Analyze a single enriched JSON file."""

    with open(file_path) as f:
        data = json.load(f)

    # Get study ID from filename or data
    study_id = file_path.stem.replace("_enriched", "")

    samples = data.get("samples", {})
    num_samples = len(samples)

    # Fields to track (comprehensive list from SampleMetadata schema)
    fields = [
        # Core biological metadata
        "organism", "tissue", "cell_type", "cell_line", "strain",
        "developmental_stage", "genotype", "age", "sex",
        # Experimental conditions
        "condition", "treatment", "timepoint", "replicate", "batch",
        # Disease/perturbation
        "disease", "stress", "temperature", "growth_condition"
    ]

    # Initialize field stats
    field_stats = {field: FieldStats(total_samples=num_samples) for field in fields}

    total_enriched = 0
    all_confidences = []

    # Analyze each sample
    for sample_id, sample_data in samples.items():
        # Determine format
        is_traceable = "biological_metadata" in sample_data

        if is_traceable:
            bio_meta = sample_data.get("biological_metadata", {})

            for field in fields:
                if field in bio_meta:
                    field_obj = bio_meta[field]
                    value = field_obj.get("value")
                    source = field_obj.get("source", "")
                    confidence = field_obj.get("confidence", 0.0)

                    if value:
                        field_stats[field].filled_after += 1

                        # Check if it was enriched (added by LLM)
                        if source == "llm_enrichment":
                            field_stats[field].enriched_count += 1
                            field_stats[field].confidence_scores.append(confidence)
                            all_confidences.append(confidence)
                            total_enriched += 1
                        else:
                            # Was already present
                            field_stats[field].filled_before += 1
        else:
            # Flat format
            for field in fields:
                if field in sample_data:
                    field_obj = sample_data.get(field)
                    if isinstance(field_obj, dict):
                        value = field_obj.get("value")
                        source = field_obj.get("source", "")
                        confidence = field_obj.get("confidence", 0.0)

                        if value:
                            field_stats[field].filled_after += 1

                            if source == "llm_enrichment":
                                field_stats[field].enriched_count += 1
                                field_stats[field].confidence_scores.append(confidence)
                                all_confidences.append(confidence)
                                total_enriched += 1
                            else:
                                field_stats[field].filled_before += 1
                    elif field_obj:  # Simple value
                        field_stats[field].filled_after += 1
                        field_stats[field].filled_before += 1

    # Calculate average confidences
    for field, stats in field_stats.items():
        if stats.confidence_scores:
            stats.avg_confidence = sum(stats.confidence_scores) / len(stats.confidence_scores)

    avg_confidence = sum(all_confidences) / len(all_confidences) if all_confidences else 0.0

    study_stats = StudyStats(
        study_id=study_id,
        num_samples=num_samples,
        fields_enriched=total_enriched,
        avg_confidence=avg_confidence,
        field_stats=field_stats,
    )

    if verbose:
        print(f"\n{study_id}:")
        print(f"  Samples: {num_samples}")
        print(f"  Fields enriched: {total_enriched}")
        print(f"  Avg confidence: {avg_confidence:.2f}")

    return study_stats


def print_summary_report(all_stats: List[StudyStats]):
    """Print a summary report of all studies."""

    print("\n" + "=" * 80)
    print("ENRICHMENT PERFORMANCE SUMMARY")
    print("=" * 80)

    # Overall stats
    total_studies = len(all_stats)
    total_samples = sum(s.num_samples for s in all_stats)
    total_enriched = sum(s.fields_enriched for s in all_stats)

    print(f"\nOverall:")
    print(f"  Studies processed: {total_studies}")
    print(f"  Total samples: {total_samples}")
    print(f"  Total fields enriched: {total_enriched}")
    print(f"  Avg fields per sample: {total_enriched / total_samples:.2f}")

    # Aggregate field stats
    print(f"\n{'Field':<20} {'Coverage Before':<18} {'Coverage After':<18} {'Enriched':<12} {'Avg Confidence':<15}")
    print("-" * 85)

    fields = [
        # Core biological metadata
        "organism", "tissue", "cell_type", "cell_line", "strain",
        "developmental_stage", "genotype", "age", "sex",
        # Experimental conditions
        "condition", "treatment", "timepoint", "replicate", "batch",
        # Disease/perturbation
        "disease", "stress", "temperature", "growth_condition"
    ]

    for field in fields:
        # Aggregate across all studies
        total_samples = sum(s.field_stats[field].total_samples for s in all_stats)
        filled_before = sum(s.field_stats[field].filled_before for s in all_stats)
        filled_after = sum(s.field_stats[field].filled_after for s in all_stats)
        enriched = sum(s.field_stats[field].enriched_count for s in all_stats)

        all_confidences = []
        for s in all_stats:
            all_confidences.extend(s.field_stats[field].confidence_scores)

        avg_conf = sum(all_confidences) / len(all_confidences) if all_confidences else 0.0

        cov_before = (filled_before / total_samples * 100) if total_samples > 0 else 0.0
        cov_after = (filled_after / total_samples * 100) if total_samples > 0 else 0.0
        improvement = cov_after - cov_before

        print(f"{field:<20} {cov_before:>6.1f}% ({filled_before:>3}) → {cov_after:>6.1f}% ({filled_after:>3})   +{enriched:<4} ({improvement:>+5.1f}%)   {avg_conf:.3f}")

    # Per-study breakdown
    print(f"\nPer-Study Breakdown:")
    print(f"{'Study ID':<30} {'Samples':<10} {'Enriched':<12} {'Avg Conf':<12}")
    print("-" * 65)

    for stats in sorted(all_stats, key=lambda s: s.fields_enriched, reverse=True):
        print(f"{stats.study_id:<30} {stats.num_samples:<10} {stats.fields_enriched:<12} {stats.avg_confidence:.3f}")


def main():
    import argparse

    parser = argparse.ArgumentParser(description="Assess enrichment performance")
    parser.add_argument("files", nargs="+", help="Enriched JSON files to analyze")
    parser.add_argument("--verbose", "-v", action="store_true", help="Verbose output")
    parser.add_argument("--output", "-o", help="Save report to JSON file")

    args = parser.parse_args()

    # Analyze all files
    all_stats = []
    for file_path in args.files:
        path = Path(file_path)
        if not path.exists():
            print(f"Warning: {file_path} not found", file=sys.stderr)
            continue

        stats = analyze_enriched_file(path, verbose=args.verbose)
        all_stats.append(stats)

    if not all_stats:
        print("No files analyzed!", file=sys.stderr)
        sys.exit(1)

    # Print summary
    print_summary_report(all_stats)

    # Save to JSON if requested
    if args.output:
        output_data = {
            "summary": {
                "total_studies": len(all_stats),
                "total_samples": sum(s.num_samples for s in all_stats),
                "total_enriched": sum(s.fields_enriched for s in all_stats),
            },
            "studies": [asdict(s) for s in all_stats],
        }

        with open(args.output, "w") as f:
            json.dump(output_data, f, indent=2)

        print(f"\nReport saved to {args.output}")


if __name__ == "__main__":
    main()
