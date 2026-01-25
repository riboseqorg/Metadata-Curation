"""Command-line interface for omics-extractor."""

import sys
import json
import argparse
import os
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Optional

from .extraction.builder import build_project_metadata
from .extraction.llm_extractor import enrich_sample_metadata
from .extraction.batch_enricher import enrich_batch_from_files
from .extraction.llm_providers import create_provider
from .schemas.base import BaseProvenance, SampleMetadata
from .output.formatters import (
    MetadataFormatter,
    export_to_csv,
    create_provenance_summary,
)
from .output.traceable_format import create_traceable_report
import yaml


def serialize_metadata(obj):
    """Custom JSON serializer for Pydantic models and datetime."""
    if isinstance(obj, BaseProvenance):
        return obj.model_dump()
    elif isinstance(obj, datetime):
        return obj.isoformat()
    elif hasattr(obj, "model_dump"):
        return obj.model_dump()
    else:
        return str(obj)


def extract_command(args):
    """Extract metadata for a BioProject."""
    print(f"Extracting metadata for {args.bioproject}")
    print("-" * 60)

    try:
        # Build complete metadata
        metadata = build_project_metadata(args.bioproject)

        # Create traceable report (hierarchical format with all levels)
        output = create_traceable_report(
            study=metadata["study"],
            samples=metadata["samples"],
            runs=metadata["runs"],
            include_raw=True,
            include_statistics=True,
        )

        # Determine output path
        output_path = Path(args.output) if args.output else Path(f"{args.bioproject}_metadata.json")

        # Write JSON output
        with open(output_path, "w") as f:
            json.dump(output, f, indent=2, default=serialize_metadata)

        # Export to CSV if requested
        if getattr(args, 'csv', False):
            csv_path = output_path.with_suffix('.csv')
            export_to_csv(
                metadata["samples"],
                str(csv_path),
                include_provenance=getattr(args, 'csv_provenance', False)
            )
            print(f"  CSV export: {csv_path}")

        # Create provenance summary if requested
        if getattr(args, 'provenance_report', False):
            summary_path = output_path.with_suffix('.md')
            summary = create_provenance_summary(metadata["samples"])
            with open(summary_path, "w") as f:
                f.write(summary)
            print(f"  Provenance report: {summary_path}")

        print("-" * 60)
        print(f"✓ Metadata extracted successfully!")
        print(f"  Study: {metadata['study'].title.value}")
        print(f"  Samples: {len(metadata['samples'])}")
        print(f"  Runs: {len(metadata['runs'])}")
        print(f"  Format: traceable (hierarchical with all provenance levels)")
        print(f"  Output: {output_path}")

        # Show extraction statistics
        stats = output["extraction_statistics"]
        print(f"\n  Field coverage:")
        for field, data in stats["field_coverage"].items():
            if data["count"] > 0:
                print(f"    - {field}: {data['count']}/{len(metadata['samples'])} ({data['percentage']}%)")

        print(f"\n  Confidence distribution:")
        conf = stats["confidence_distribution"]["overall"]
        total = conf["total_fields"]
        if total > 0:
            print(f"    - High (≥0.8): {conf['high (≥0.8)']} ({100*conf['high (≥0.8)']/total:.0f}%)")
            print(f"    - Medium (0.5-0.8): {conf['medium (0.5-0.8)']} ({100*conf['medium (0.5-0.8)']/total:.0f}%)")
            print(f"    - Low (<0.5): {conf['low (<0.5)']} ({100*conf['low (<0.5)']/total:.0f}%)")

        print(f"\n  Ontology mapping: {stats['ontology_mapping']['overall']['mapped']}/{stats['ontology_mapping']['overall']['total']} ({stats['ontology_mapping']['overall']['percentage']}%)")

        # Show navigation guide
        print(f"\n  Access levels:")
        print(f"    - Quick: samples[<id>].quick_view.<field>")
        print(f"    - Standard: samples[<id>].biological_metadata.<field>.provenance")
        print(f"    - Detailed: samples[<id>].biological_metadata.<field>.extraction_details")

    except Exception as e:
        print(f"✗ Error: {e}", file=sys.stderr)
        if args.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)


def validate_command(args):
    """Validate extracted metadata."""
    print(f"Validating metadata from {args.input}")
    print("-" * 60)

    try:
        with open(args.input, "r") as f:
            data = json.load(f)

        # Simple validation checks
        issues = []
        warnings = []

        # Check study metadata
        study = data.get("study", {})
        if not study.get("title", {}).get("value"):
            issues.append("Study title is missing")

        if study.get("organism", {}).get("confidence", 0) < 0.5:
            warnings.append("Low confidence organism extraction")

        # Check samples
        samples = data.get("samples", {})
        for sample_id, sample in samples.items():
            # Check for tissue (priority field)
            if not sample.get("tissue"):
                warnings.append(f"Sample {sample_id}: No tissue information")

            # Check organism
            if not sample.get("organism", {}).get("value"):
                issues.append(f"Sample {sample_id}: Missing organism")

        # Check runs
        runs = data.get("runs", [])
        if not runs:
            issues.append("No runs found")

        # Report results
        print(f"Validated {len(samples)} samples and {len(runs)} runs")

        if issues:
            print(f"\n❌ Issues found ({len(issues)}):")
            for issue in issues:
                print(f"  - {issue}")

        if warnings:
            print(f"\n⚠️  Warnings ({len(warnings)}):")
            for warning in warnings[:10]:  # Show first 10
                print(f"  - {warning}")
            if len(warnings) > 10:
                print(f"  ... and {len(warnings) - 10} more")

        if not issues and not warnings:
            print("✓ All validation checks passed!")

        sys.exit(1 if issues else 0)

    except Exception as e:
        print(f"✗ Error: {e}", file=sys.stderr)
        if args.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)


def _batch_enrich_with_provider(
    input_files: List[Path],
    output_dir: Path,
    provider,
    resume: bool = False,
) -> Dict:
    """Batch enrich files using a pre-loaded provider (for local models)."""
    stats = {
        "total_files": len(input_files),
        "processed": 0,
        "failed": 0,
        "total_samples": 0,
        "total_enriched": 0,
    }

    for i, input_path in enumerate(input_files, 1):
        output_path = output_dir / f"{input_path.stem}_enriched.json"

        # Skip if resume and output exists
        if resume and output_path.exists():
            print(f"[{i}/{len(input_files)}] Skipping {input_path.name} (already exists)")
            stats["processed"] += 1
            continue

        print(f"[{i}/{len(input_files)}] Processing {input_path.name}...")

        try:
            with open(input_path) as f:
                data = json.load(f)

            # Get study context
            study = data.get("study", {})
            study_title = study.get("title", {}).get("value", "")
            study_description = study.get("description", {}).get("value", "")
            abstract = study.get("abstract", {}).get("value")

            samples = data.get("samples", {})
            stats["total_samples"] += len(samples)

            enriched_count = 0
            for sample_id, sample_data in samples.items():
                try:
                    # Extract from traceable format to flat SampleMetadata
                    quick_view = sample_data.get("quick_view", {})
                    bio_meta = sample_data.get("biological_metadata", {})

                    flat_data = {
                        "sample_id": quick_view.get("sample_id"),
                        "bioproject_id": quick_view.get("bioproject_id"),
                        "organism": quick_view.get("organism"),
                    }

                    # Extract biological fields
                    for field in ["tissue", "cell_type", "cell_line", "strain", "treatment", "disease", "developmental_stage", "age", "sex"]:
                        if field in bio_meta and bio_meta[field].get("value"):
                            flat_data[field] = bio_meta[field]["value"]

                    sample = SampleMetadata(**flat_data)
                    sample_title = sample_data.get("sample_title", {}).get("value")
                    sample_description = sample_data.get("sample_description", {}).get("value")

                    enriched_sample = enrich_sample_metadata(
                        sample=sample,
                        study_title=study_title,
                        study_description=study_description,
                        study_abstract=abstract,
                        sample_title=sample_title,
                        sample_description=sample_description,
                        provider=provider,
                    )

                    # Update traceable format with enriched values
                    for field in ["tissue", "cell_type", "cell_line", "strain", "treatment", "disease", "developmental_stage", "age", "sex"]:
                        new_val = getattr(enriched_sample, field, None)
                        old_val = bio_meta.get(field, {}).get("value")
                        if new_val and not old_val:
                            if field not in bio_meta:
                                bio_meta[field] = {}
                            bio_meta[field]["value"] = new_val
                            bio_meta[field]["source"] = "llm_enrichment"
                            bio_meta[field]["confidence"] = 0.8
                            quick_view[field] = new_val
                            enriched_count += 1

                except Exception as e:
                    print(f"  ⚠️  Failed to enrich sample {sample_id}: {e}")
                    continue

            # Save enriched data
            data["samples"] = samples
            data["enrichment_metadata"] = {
                "timestamp": datetime.now().isoformat(),
                "samples_enriched": enriched_count,
            }

            with open(output_path, "w") as f:
                json.dump(data, f, indent=2, default=serialize_metadata)

            stats["processed"] += 1
            stats["total_enriched"] += enriched_count
            print(f"  ✓ Enriched {enriched_count}/{len(samples)} samples")

        except Exception as e:
            print(f"  ✗ Failed: {e}")
            stats["failed"] += 1
            continue

    return stats


def load_provider_from_config(model_name: str, config_path: Path):
    """Load LLM provider from config file."""
    if not config_path.exists():
        print(f"✗ Error: Config file not found: {config_path}", file=sys.stderr)
        print(f"   Copy config/model_paths.example.yaml to config/model_paths.yaml", file=sys.stderr)
        sys.exit(1)

    with open(config_path) as f:
        config = yaml.safe_load(f)

    model_config = config.get("local_models", {}).get(model_name)
    if not model_config:
        print(f"✗ Error: Model '{model_name}' not found in config", file=sys.stderr)
        print(f"   Available models: {', '.join(config.get('local_models', {}).keys())}", file=sys.stderr)
        sys.exit(1)

    provider_type = model_config.get("provider", "vllm")
    print(f"Loading {model_name} with {provider_type}...")

    if provider_type == "vllm":
        return create_provider(
            "vllm",
            model_path=model_config["path"],
            tensor_parallel_size=model_config.get("tensor_parallel_size", 1),
            gpu_memory_utilization=model_config.get("gpu_memory_utilization", 0.9),
        )
    elif provider_type == "transformers":
        return create_provider(
            "transformers",
            model_path=model_config["path"],
            device=model_config.get("device", "cuda"),
        )
    else:
        print(f"✗ Error: Unknown provider type: {provider_type}", file=sys.stderr)
        sys.exit(1)


def enrich_command(args):
    """Enrich metadata using LLM (Phase 2 - GPU processing)."""
    print(f"Enriching metadata from {args.input}")
    print("-" * 60)

    try:
        # Load metadata
        input_path = Path(args.input)
        if not input_path.exists():
            print(f"✗ Error: File not found: {args.input}", file=sys.stderr)
            sys.exit(1)

        with open(input_path) as f:
            data = json.load(f)

        # Create provider based on mode
        if args.model:
            # Local model mode
            config_path = Path(args.config) if args.config else Path("config/model_paths.yaml")
            provider = load_provider_from_config(args.model, config_path)
        else:
            # Claude API mode (default)
            api_key = args.api_key or os.environ.get("ANTHROPIC_API_KEY")
            if not api_key:
                print(
                    "✗ Error: No API key provided. Set ANTHROPIC_API_KEY environment variable or use --api-key",
                    file=sys.stderr,
                )
                sys.exit(1)
            provider = create_provider("claude", api_key=api_key, model=args.claude_model or "claude-sonnet-4-5-20250929")

        # Extract study context
        study = data.get("study", {})
        study_title = study.get("title", {}).get("value", "")
        study_description = study.get("description", {}).get("value", "")
        abstract = study.get("abstract", {}).get("value")

        # Get samples
        samples = data.get("samples", {})
        total_samples = len(samples)

        if total_samples == 0:
            print("⚠️  No samples found in metadata file")
            sys.exit(0)

        # Determine which samples need enrichment
        samples_to_enrich = []
        for sample_id, sample_data in samples.items():
            # Check if any priority fields are missing
            missing_fields = []
            if not sample_data.get("tissue"):
                missing_fields.append("tissue")
            if not sample_data.get("cell_type"):
                missing_fields.append("cell_type")
            if not sample_data.get("cell_line"):
                missing_fields.append("cell_line")

            if missing_fields:
                samples_to_enrich.append((sample_id, sample_data, missing_fields))

        if args.only_if_missing and not samples_to_enrich:
            print("✓ All samples have complete metadata, no enrichment needed")
            sys.exit(0)

        if args.only_if_missing:
            print(f"Found {len(samples_to_enrich)}/{total_samples} samples with missing fields")
        else:
            samples_to_enrich = [(sid, sdata, []) for sid, sdata in samples.items()]
            print(f"Enriching all {total_samples} samples")

        # Track enrichment statistics
        enriched_count = 0
        fields_added = {"tissue": 0, "cell_type": 0, "cell_line": 0, "treatment": 0, "strain": 0}

        # Enrich each sample
        print("\nEnriching samples...")
        for i, (sample_id, sample_data, missing) in enumerate(samples_to_enrich, 1):
            if args.verbose:
                print(f"  [{i}/{len(samples_to_enrich)}] {sample_id} (missing: {', '.join(missing) or 'none'})")

            # Extract from traceable format to flat SampleMetadata
            quick_view = sample_data.get("quick_view", {})
            bio_meta = sample_data.get("biological_metadata", {})

            flat_data = {
                "sample_id": quick_view.get("sample_id"),
                "bioproject_id": quick_view.get("bioproject_id"),
                "organism": quick_view.get("organism"),
            }

            # Extract biological fields
            for field in ["tissue", "cell_type", "cell_line", "strain", "treatment", "disease", "developmental_stage", "age", "sex"]:
                if field in bio_meta and bio_meta[field].get("value"):
                    flat_data[field] = bio_meta[field]["value"]

            sample = SampleMetadata(**flat_data)

            # Get sample-specific context
            sample_title = sample_data.get("sample_title", {}).get("value")
            sample_description = sample_data.get("sample_description", {}).get("value")

            # Enrich with LLM
            try:
                enriched_sample = enrich_sample_metadata(
                    sample=sample,
                    study_title=study_title,
                    study_description=study_description,
                    sample_title=sample_title,
                    sample_description=sample_description,
                    abstract=abstract,
                    api_key=api_key,
                    source_id=f"llm_enrichment_{sample_id}",
                )

                # Update traceable format with enriched values
                for field in ["tissue", "cell_type", "cell_line", "strain", "treatment"]:
                    new_val = getattr(enriched_sample, field, None)
                    old_val = bio_meta.get(field, {}).get("value")
                    if new_val and not old_val:
                        if field not in bio_meta:
                            bio_meta[field] = {}
                        bio_meta[field]["value"] = new_val
                        bio_meta[field]["source"] = "llm_enrichment"
                        bio_meta[field]["confidence"] = 0.8
                        quick_view[field] = new_val
                        fields_added[field] += 1

                enriched_count += 1

            except Exception as e:
                if args.verbose:
                    print(f"    ⚠️  Failed to enrich {sample_id}: {e}")
                continue

        # Save enriched metadata
        output_path = Path(args.output) if args.output else input_path.with_stem(f"{input_path.stem}_enriched")

        data["samples"] = samples
        data["enrichment_metadata"] = {
            "timestamp": datetime.now().isoformat(),
            "samples_enriched": enriched_count,
            "fields_added": fields_added,
        }

        with open(output_path, "w") as f:
            json.dump(data, f, indent=2, default=serialize_metadata)

        # Report results
        print(f"\n✓ Enriched {enriched_count}/{len(samples_to_enrich)} samples")
        print(f"  Fields added:")
        for field, count in fields_added.items():
            if count > 0:
                print(f"    - {field}: {count}")
        print(f"\nSaved to: {output_path}")

    except Exception as e:
        print(f"✗ Error: {e}", file=sys.stderr)
        if args.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)


def batch_enrich_command(args):
    """Batch enrich multiple metadata files (optimized for GPU nodes)."""
    print(f"Batch enriching {len(args.inputs)} files")
    print("-" * 60)

    try:
        # Create provider based on mode
        if args.model:
            # Local model mode
            config_path = Path(args.config) if args.config else Path("config/model_paths.yaml")
            provider = load_provider_from_config(args.model, config_path)
        else:
            # Claude API mode (default)
            api_key = args.api_key or os.environ.get("ANTHROPIC_API_KEY")
            if not api_key:
                print(
                    "✗ Error: No API key provided. Set ANTHROPIC_API_KEY environment variable or use --api-key",
                    file=sys.stderr,
                )
                sys.exit(1)
            provider = create_provider("claude", api_key=api_key, model=args.claude_model or "claude-sonnet-4-5-20250929")

        # Parse input files
        input_files = []
        for pattern in args.inputs:
            path = Path(pattern)
            if path.is_file():
                input_files.append(path)
            elif "*" in pattern:
                # Glob pattern
                input_files.extend(Path().glob(pattern))
            else:
                print(f"⚠️  Skipping invalid input: {pattern}")

        if not input_files:
            print("✗ Error: No valid input files found", file=sys.stderr)
            sys.exit(1)

        # Create output directory
        output_dir = Path(args.output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        print(f"Processing {len(input_files)} files...")
        print(f"Output directory: {output_dir}")

        if args.model:
            # Local model mode - process sequentially with loaded model
            print(f"Using local model (sequential processing)")
            stats = _batch_enrich_with_provider(
                input_files=input_files,
                output_dir=output_dir,
                provider=provider,
                resume=args.resume,
            )
        else:
            # Claude API mode - use concurrent workers
            print(f"Workers: {args.workers}, Rate limit: {args.rate_limit} req/min")
            stats = enrich_batch_from_files(
                input_files=input_files,
                output_dir=output_dir,
                api_key=api_key,
                max_workers=args.workers,
                rate_limit_rpm=args.rate_limit,
                resume=args.resume,
            )

        # Report results
        print("\n" + "=" * 60)
        print("Batch enrichment complete!")
        print(f"  Files processed: {stats['processed']}/{stats['total_files']}")
        print(f"  Files failed: {stats['failed']}")
        print(f"  Total samples: {stats['total_samples']}")
        print(f"  Samples enriched: {stats['total_enriched']}")
        print(f"\nOutput files saved to: {output_dir}")

    except Exception as e:
        print(f"✗ Error: {e}", file=sys.stderr)
        if args.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)


def main():
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(
        description="Extract metadata from omics sequencing projects",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Phase 1: Extract metadata (CPU nodes)
  omics-extract extract PRJNA1170270

  # Phase 2: Enrich with local model (GPU nodes)
  omics-extract enrich PRJNA1170270_metadata.json --model qwen-2.5-7b

  # Enrich with Claude API
  omics-extract enrich PRJNA1170270_metadata.json

  # Batch enrich multiple projects with local model
  omics-extract batch-enrich *_metadata.json --output-dir enriched/ --model mistral-7b

  # Batch enrich with Claude API (GPU cluster)
  omics-extract batch-enrich *_metadata.json --output-dir enriched/ --workers 8

  # Resume interrupted batch job
  omics-extract batch-enrich *_metadata.json --output-dir enriched/ --model qwen-2.5-7b --resume

  # Full workflow with local model
  omics-extract extract PRJNA1170270 --output data.json
  omics-extract enrich data.json --model qwen-2.5-7b --only-if-missing
        """,
    )

    parser.add_argument("--verbose", "-v", action="store_true", help="Verbose output")

    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    # Extract command
    extract_parser = subparsers.add_parser("extract", help="Extract metadata for a BioProject (Phase 1 - CPU)")
    extract_parser.add_argument("bioproject", help="BioProject accession (PRJNA...)")
    extract_parser.add_argument(
        "--output", "-o", help="Output JSON file (default: <bioproject>_metadata.json)"
    )
    extract_parser.add_argument(
        "--detail",
        choices=["simple", "standard", "detailed"],
        default="standard",
        help="Detail level for output (default: standard). "
             "simple=values only, standard=values+provenance, detailed=complete extraction chain"
    )
    extract_parser.add_argument(
        "--csv",
        action="store_true",
        help="Also export to CSV format"
    )
    extract_parser.add_argument(
        "--csv-provenance",
        action="store_true",
        help="Include source/confidence columns in CSV (requires --csv)"
    )
    extract_parser.add_argument(
        "--provenance-report",
        action="store_true",
        help="Generate markdown provenance summary report"
    )
    extract_parser.set_defaults(func=extract_command)

    # Enrich command (Phase 2)
    enrich_parser = subparsers.add_parser("enrich", help="Enrich metadata with LLM (Phase 2 - GPU)")
    enrich_parser.add_argument("input", help="Input metadata JSON file from extract command")
    enrich_parser.add_argument(
        "--output", "-o", help="Output JSON file (default: <input>_enriched.json)"
    )
    enrich_parser.add_argument(
        "--model", help="Local model name from config (e.g., qwen-2.5-7b, mistral-7b)"
    )
    enrich_parser.add_argument(
        "--config", help="Path to model config YAML (default: config/model_paths.yaml)"
    )
    enrich_parser.add_argument(
        "--api-key", help="Anthropic API key (or set ANTHROPIC_API_KEY env var) - for Claude mode"
    )
    enrich_parser.add_argument(
        "--claude-model", help="Claude model name (default: claude-sonnet-4-5-20250929)"
    )
    enrich_parser.add_argument(
        "--only-if-missing",
        action="store_true",
        help="Only enrich samples with missing critical fields",
    )
    enrich_parser.set_defaults(func=enrich_command)

    # Batch enrich command (Phase 2 - optimized for GPU clusters)
    batch_parser = subparsers.add_parser("batch-enrich", help="Batch enrich multiple files (GPU-optimized)")
    batch_parser.add_argument("inputs", nargs="+", help="Input metadata files or glob patterns")
    batch_parser.add_argument(
        "--output-dir", "-o", required=True, help="Output directory for enriched files"
    )
    batch_parser.add_argument(
        "--model", help="Local model name from config (e.g., qwen-2.5-7b, mistral-7b)"
    )
    batch_parser.add_argument(
        "--config", help="Path to model config YAML (default: config/model_paths.yaml)"
    )
    batch_parser.add_argument(
        "--api-key", help="Anthropic API key (or set ANTHROPIC_API_KEY env var) - for Claude mode"
    )
    batch_parser.add_argument(
        "--claude-model", help="Claude model name (default: claude-sonnet-4-5-20250929)"
    )
    batch_parser.add_argument(
        "--workers", type=int, default=4, help="Concurrent API requests (default: 4, ignored for local models)"
    )
    batch_parser.add_argument(
        "--rate-limit", type=int, default=50, help="Max requests per minute (default: 50, ignored for local models)"
    )
    batch_parser.add_argument(
        "--resume", action="store_true", help="Resume from checkpoints"
    )
    batch_parser.set_defaults(func=batch_enrich_command)

    # Validate command
    validate_parser = subparsers.add_parser("validate", help="Validate extracted metadata")
    validate_parser.add_argument("input", help="Input JSON file to validate")
    validate_parser.set_defaults(func=validate_command)

    # Parse and execute
    args = parser.parse_args()

    if not args.command:
        parser.print_help()
        sys.exit(1)

    args.func(args)


if __name__ == "__main__":
    main()
