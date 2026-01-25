#!/usr/bin/env python3
"""
Review enriched samples for quality assurance.

Shows enriched values alongside source text to verify accuracy.
Useful for manual spot-checking before scaling up.

Usage:
    # Review 10 random samples
    python scripts/review_enrichments.py enriched/*_enriched.json -n 10

    # Review 20 samples, compact view
    python scripts/review_enrichments.py enriched/*_enriched.json -n 20 --compact

    # Review only low-confidence extractions (<0.6)
    python scripts/review_enrichments.py enriched/*_enriched.json --max-confidence 0.6

    # Review specific study
    python scripts/review_enrichments.py enriched/PRJNA123456_enriched.json
"""

import json
import random
import sys
from pathlib import Path
from typing import List, Dict, Optional, Tuple
from dataclasses import dataclass
import argparse


@dataclass
class EnrichedField:
    """An enriched field with its source context."""

    field_name: str
    value: str
    confidence: float
    sample_id: str
    study_id: str

    # Source text
    study_title: str
    study_description: str
    study_abstract: str
    sample_title: str
    sample_description: str

    @property
    def all_source_text(self) -> str:
        """Combined source text."""
        parts = [
            self.study_title,
            self.study_description,
            self.study_abstract,
            self.sample_title,
            self.sample_description,
        ]
        return " ".join(p for p in parts if p)

    def appears_in_source(self) -> bool:
        """Check if value appears in source text."""
        return self.value.lower() in self.all_source_text.lower()


def extract_enriched_fields(file_path: Path, min_confidence: float = 0.0, max_confidence: float = 1.0) -> List[EnrichedField]:
    """Extract all enriched fields from a file."""

    with open(file_path) as f:
        data = json.load(f)

    study_id = file_path.stem.replace("_enriched", "")
    study = data.get("study", {})

    # Get study-level context
    study_title = study.get("title", {}).get("value", "")
    study_description = study.get("description", {}).get("value", "")
    study_abstract = study.get("paper_abstract", {}).get("value", "")

    enriched_fields = []

    samples = data.get("samples", {})
    for sample_id, sample_data in samples.items():
        # Get sample context
        sample_title = ""
        sample_description = ""

        # Check format
        is_traceable = "biological_metadata" in sample_data

        if is_traceable:
            sample_title = sample_data.get("sample_title", {}).get("value", "")
            sample_description = sample_data.get("sample_description", {}).get("value", "")
            bio_meta = sample_data.get("biological_metadata", {})
        else:
            # Flat format
            if "sample_title" in sample_data and isinstance(sample_data["sample_title"], dict):
                sample_title = sample_data["sample_title"].get("value", "")
            if "sample_description" in sample_data and isinstance(sample_data["sample_description"], dict):
                sample_description = sample_data["sample_description"].get("value", "")
            bio_meta = sample_data

        # Extract enriched fields
        fields_to_check = ["tissue", "cell_type", "cell_line", "strain", "treatment",
                          "disease", "developmental_stage", "age", "sex", "genotype"]

        for field in fields_to_check:
            if field in bio_meta and isinstance(bio_meta[field], dict):
                field_obj = bio_meta[field]

                # Only include LLM-enriched fields
                if field_obj.get("source") == "llm_enrichment":
                    value = field_obj.get("value")
                    confidence = field_obj.get("confidence", 0.0)

                    # Filter by confidence
                    if value and min_confidence <= confidence <= max_confidence:
                        enriched_fields.append(EnrichedField(
                            field_name=field,
                            value=value,
                            confidence=confidence,
                            sample_id=sample_id,
                            study_id=study_id,
                            study_title=study_title,
                            study_description=study_description,
                            study_abstract=study_abstract,
                            sample_title=sample_title,
                            sample_description=sample_description,
                        ))

    return enriched_fields


def highlight_value_in_text(text: str, value: str, max_length: int = 150) -> str:
    """Highlight where value appears in text, or show full text if not found."""

    if not text:
        return "[No text available]"

    value_lower = value.lower()
    text_lower = text.lower()

    if value_lower in text_lower:
        # Find position
        pos = text_lower.index(value_lower)

        # Extract context around the value
        start = max(0, pos - 50)
        end = min(len(text), pos + len(value) + 50)

        snippet = text[start:end]

        # Add ellipsis
        if start > 0:
            snippet = "..." + snippet
        if end < len(text):
            snippet = snippet + "..."

        # Highlight the value (case-insensitive)
        import re
        pattern = re.compile(re.escape(value), re.IGNORECASE)
        snippet = pattern.sub(f"**{value.upper()}**", snippet)

        return snippet
    else:
        # Value not found, show beginning of text
        if len(text) > max_length:
            return text[:max_length] + "..."
        return text


def print_text_review(fields: List[EnrichedField], compact: bool = False):
    """Print review in text format."""

    print("\n" + "=" * 80)
    print("ENRICHMENT QUALITY REVIEW")
    print("=" * 80)
    print(f"\nReviewing {len(fields)} enriched fields")

    # Summary stats
    found_count = sum(1 for f in fields if f.appears_in_source())
    avg_conf = sum(f.confidence for f in fields) / len(fields)

    print(f"Found in source: {found_count}/{len(fields)} ({found_count/len(fields)*100:.1f}%)")
    print(f"Avg confidence: {avg_conf:.3f}")
    print("\nLook for: hallucinations, contradictions, or unsupported inferences.\n")

    for i, field in enumerate(fields, 1):
        print("\n" + "─" * 80)

        # Header line
        found_marker = "✓" if field.appears_in_source() else "✗"
        conf_color = "HIGH" if field.confidence >= 0.8 else "MED" if field.confidence >= 0.6 else "LOW"

        print(f"[{i}/{len(fields)}] {field.field_name}: {field.value}")
        print(f"  Study: {field.study_id} | Sample: {field.sample_id}")
        print(f"  Confidence: {field.confidence:.3f} ({conf_color}) | In source: {found_marker}")

        if not compact:
            # Show source text
            if field.study_title:
                print(f"\n  Study: {field.study_title}")

            if field.study_description:
                highlighted = highlight_value_in_text(field.study_description, field.value, max_length=200)
                print(f"  Desc:  {highlighted}")

            if field.sample_title:
                print(f"\n  Sample: {field.sample_title}")

            if field.sample_description:
                highlighted = highlight_value_in_text(field.sample_description, field.value, max_length=200)
                print(f"  Desc:   {highlighted}")

            if field.study_abstract:
                highlighted = highlight_value_in_text(field.study_abstract, field.value, max_length=150)
                print(f"  Abstract: {highlighted}")

    print("\n" + "=" * 80)



def main():
    parser = argparse.ArgumentParser(
        description="Review enriched fields for quality assurance",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Review 10 random fields
  %(prog)s enriched/*.json -n 10

  # Review only low confidence (<0.6)
  %(prog)s enriched/*.json --max-confidence 0.6

  # Compact view (less text)
  %(prog)s enriched/*.json -n 20 --compact

  # Review specific study
  %(prog)s enriched/PRJNA123456_enriched.json
        """
    )
    parser.add_argument("files", nargs="+", help="Enriched JSON files to review")
    parser.add_argument("--num-samples", "-n", type=int, help="Number of random samples to review (default: all)")
    parser.add_argument("--min-confidence", type=float, default=0.0, help="Minimum confidence to include (default: 0.0)")
    parser.add_argument("--max-confidence", type=float, default=1.0, help="Maximum confidence to include (default: 1.0)")
    parser.add_argument("--compact", "-c", action="store_true", help="Compact view (header only, no source text)")
    parser.add_argument("--seed", type=int, help="Random seed for reproducible sampling")

    args = parser.parse_args()

    if args.seed:
        random.seed(args.seed)

    # Collect all enriched fields
    all_fields = []
    for file_path in args.files:
        path = Path(file_path)
        if not path.exists():
            print(f"Warning: {file_path} not found", file=sys.stderr)
            continue

        fields = extract_enriched_fields(path, args.min_confidence, args.max_confidence)
        all_fields.extend(fields)

    if not all_fields:
        print("No enriched fields found matching criteria!", file=sys.stderr)
        sys.exit(1)

    # Sample if requested
    if args.num_samples and args.num_samples < len(all_fields):
        all_fields = random.sample(all_fields, args.num_samples)

    # Print review
    print_text_review(all_fields, compact=args.compact)


if __name__ == "__main__":
    main()
