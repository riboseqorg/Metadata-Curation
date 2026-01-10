#!/usr/bin/env python3
"""
Batch extract metadata for 100 RiboSeq studies with traceable output.

This script demonstrates the traceable output format at scale:
- Processes 100 RiboSeq studies
- Outputs hierarchical format with all provenance levels
- Creates summary statistics across all studies
"""

import json
import sys
from pathlib import Path
from datetime import datetime
from typing import List, Dict, Any
import traceback

from omics_extractor.extraction.builder import build_project_metadata
from omics_extractor.output.traceable_format import create_traceable_report


# First 100 RiboSeq BioProjects from RiboSeq.org
RIBOSEQ_PROJECTS = [
    "PRJDB10544", "PRJDB10799", "PRJDB10882", "PRJDB11308", "PRJDB13770",
    "PRJDB15060", "PRJDB16373", "PRJDB2960", "PRJDB9278", "PRJDB9716",
    "PRJEB12126", "PRJEB14170", "PRJEB17636", "PRJEB18913", "PRJEB19102",
    "PRJEB19678", "PRJEB21099", "PRJEB21199", "PRJEB21224", "PRJEB23298",
    "PRJEB23398", "PRJEB23562", "PRJEB23773", "PRJEB25491", "PRJEB26279",
    "PRJEB26593", "PRJEB27418", "PRJEB27604", "PRJEB28203", "PRJEB28810",
    "PRJEB29208", "PRJEB30076", "PRJEB31666", "PRJEB31750", "PRJEB32121",
    "PRJEB32965", "PRJEB32969", "PRJEB33244", "PRJEB33323", "PRJEB35858",
    "PRJEB36134", "PRJEB36135", "PRJEB36148", "PRJEB36274", "PRJEB36467",
    "PRJEB36468", "PRJEB36473", "PRJEB36893", "PRJEB36902", "PRJEB37021",
    "PRJEB38391", "PRJEB39403", "PRJEB39904", "PRJEB39905", "PRJEB41027",
    "PRJEB42778", "PRJEB43647", "PRJEB43705", "PRJEB44238", "PRJEB45612",
    "PRJEB45876", "PRJEB46361", "PRJEB47140", "PRJEB4801", "PRJEB50305",
    "PRJEB5136", "PRJEB51486", "PRJEB5150", "PRJEB5263", "PRJEB54638",
    "PRJEB54714", "PRJEB56594", "PRJEB57800", "PRJEB5938", "PRJEB69686",
    "PRJEB7207", "PRJEB7261", "PRJEB7276", "PRJEB7282", "PRJEB7300",
    "PRJEB7301", "PRJEB7498", "PRJEB86645", "PRJEB9417", "PRJNA1002596",
    "PRJNA1002994", "PRJNA1003962", "PRJNA1004241", "PRJNA1004312", "PRJNA1007580",
    "PRJNA1008247", "PRJNA1008576", "PRJNA1009809", "PRJNA1013070", "PRJNA1016416",
    "PRJNA1017801", "PRJNA1019502", "PRJNA1019669", "PRJNA1021597", "PRJNA1022492",
]


def extract_single_project(project_id: str, output_dir: Path) -> Dict[str, Any]:
    """
    Extract metadata for a single project with traceable format.

    Returns:
        Summary dict with extraction results
    """
    print(f"\nProcessing {project_id}...")

    try:
        # Extract metadata
        metadata = build_project_metadata(project_id)

        # Create traceable report
        report = create_traceable_report(
            study=metadata["study"],
            samples=metadata["samples"],
            runs=metadata["runs"],
            include_raw=True,
            include_statistics=True,
        )

        # Save to file
        output_file = output_dir / f"{project_id}_traceable.json"
        with open(output_file, "w") as f:
            json.dump(report, f, indent=2, default=str)

        # Extract summary statistics
        stats = report["extraction_statistics"]
        summary = {
            "project_id": project_id,
            "status": "success",
            "samples": stats["total_samples"],
            "runs": stats["total_runs"],
            "completeness": stats["completeness"]["percentage"],
            "high_confidence_fields": stats["confidence_distribution"]["overall"]["high (≥0.8)"],
            "total_fields": stats["confidence_distribution"]["overall"]["total_fields"],
            "ontology_mapped": stats["ontology_mapping"]["overall"]["percentage"],
            "output_file": str(output_file),
        }

        print(f"  ✓ {stats['total_samples']} samples, {stats['total_runs']} runs")
        print(f"  ✓ {stats['completeness']['percentage']}% complete")
        print(f"  ✓ {stats['ontology_mapping']['overall']['percentage']}% ontology mapped")

        return summary

    except Exception as e:
        print(f"  ✗ Error: {e}")
        if "--verbose" in sys.argv:
            traceback.print_exc()

        return {
            "project_id": project_id,
            "status": "failed",
            "error": str(e),
        }


def compute_aggregate_statistics(summaries: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Compute aggregate statistics across all projects."""
    successful = [s for s in summaries if s["status"] == "success"]
    failed = [s for s in summaries if s["status"] == "failed"]

    if not successful:
        return {
            "total_projects": len(summaries),
            "successful": 0,
            "failed": len(failed),
        }

    total_samples = sum(s["samples"] for s in successful)
    total_runs = sum(s["runs"] for s in successful)
    total_fields = sum(s.get("total_fields", 0) for s in successful)
    high_conf_fields = sum(s.get("high_confidence_fields", 0) for s in successful)

    avg_completeness = sum(s["completeness"] for s in successful) / len(successful)
    avg_ontology_mapping = sum(s.get("ontology_mapped", 0) for s in successful) / len(successful)

    # Field coverage across all projects
    field_coverage = {}
    for s in successful:
        # Would need to read individual files for detailed breakdown
        pass

    return {
        "total_projects": len(summaries),
        "successful": len(successful),
        "failed": len(failed),
        "total_samples": total_samples,
        "total_runs": total_runs,
        "total_fields": total_fields,
        "high_confidence_fields": high_conf_fields,
        "high_confidence_percentage": round(100 * high_conf_fields / total_fields, 1) if total_fields else 0,
        "average_completeness": round(avg_completeness, 1),
        "average_ontology_mapping": round(avg_ontology_mapping, 1),
        "projects_by_completeness": {
            "complete (>90%)": len([s for s in successful if s["completeness"] > 90]),
            "good (70-90%)": len([s for s in successful if 70 <= s["completeness"] <= 90]),
            "partial (50-70%)": len([s for s in successful if 50 <= s["completeness"] < 70]),
            "low (<50%)": len([s for s in successful if s["completeness"] < 50]),
        },
    }


def main():
    """Process 100 RiboSeq studies with traceable output."""
    print("=" * 80)
    print("BATCH EXTRACTION: 100 RiboSeq Studies with Traceable Output")
    print("=" * 80)

    # Setup output directory
    output_dir = Path("riboseq_batch_output")
    output_dir.mkdir(exist_ok=True)
    print(f"\nOutput directory: {output_dir}")

    # Process all projects
    summaries = []
    start_time = datetime.now()

    for i, project_id in enumerate(RIBOSEQ_PROJECTS, 1):
        print(f"\n[{i}/{len(RIBOSEQ_PROJECTS)}] {project_id}")

        summary = extract_single_project(project_id, output_dir)
        summaries.append(summary)

    # Compute aggregate statistics
    print("\n" + "=" * 80)
    print("COMPUTING AGGREGATE STATISTICS")
    print("=" * 80)

    aggregate_stats = compute_aggregate_statistics(summaries)

    # Save summary
    summary_file = output_dir / "batch_summary.json"
    with open(summary_file, "w") as f:
        json.dump({
            "extraction_timestamp": start_time.isoformat(),
            "processing_time_seconds": (datetime.now() - start_time).total_seconds(),
            "aggregate_statistics": aggregate_stats,
            "project_summaries": summaries,
        }, f, indent=2, default=str)

    # Print results
    print(f"\n{'='*80}")
    print("BATCH EXTRACTION COMPLETE")
    print(f"{'='*80}\n")

    print(f"Total projects: {aggregate_stats['total_projects']}")
    print(f"  Successful: {aggregate_stats['successful']}")
    print(f"  Failed: {aggregate_stats['failed']}")

    print(f"\nTotal samples extracted: {aggregate_stats['total_samples']}")
    print(f"Total runs extracted: {aggregate_stats['total_runs']}")

    print(f"\nField extraction quality:")
    print(f"  Total fields extracted: {aggregate_stats['total_fields']}")
    print(f"  High confidence (≥0.8): {aggregate_stats['high_confidence_fields']} ({aggregate_stats['high_confidence_percentage']}%)")

    print(f"\nAverage completeness: {aggregate_stats['average_completeness']}%")
    print(f"Average ontology mapping: {aggregate_stats['average_ontology_mapping']}%")

    print(f"\nProjects by completeness:")
    for category, count in aggregate_stats["projects_by_completeness"].items():
        print(f"  {category}: {count}")

    print(f"\nOutput files:")
    print(f"  Individual reports: {output_dir}/*_traceable.json")
    print(f"  Batch summary: {summary_file}")

    print(f"\nProcessing time: {(datetime.now() - start_time).total_seconds():.1f} seconds")

    # Show example of traceable format
    if aggregate_stats['successful'] > 0:
        print(f"\n{'='*80}")
        print("EXAMPLE: Traceable Format Structure")
        print(f"{'='*80}\n")

        # Load first successful project
        first_success = next(s for s in summaries if s["status"] == "success")
        with open(first_success["output_file"]) as f:
            example = json.load(f)

        print("Top-level structure:")
        print(json.dumps({k: "..." for k in example.keys()}, indent=2))

        print("\nNavigation guide:")
        print(json.dumps(example["navigation_guide"], indent=2))

        if example["samples"]:
            first_sample_id = list(example["samples"].keys())[0]
            first_sample = example["samples"][first_sample_id]

            print(f"\nExample sample structure ({first_sample_id}):")
            print(json.dumps({k: "..." for k in first_sample.keys()}, indent=2))

            print(f"\nQuick view (just values):")
            print(json.dumps(first_sample["quick_view"], indent=2))

            if "tissue" in first_sample.get("biological_metadata", {}):
                print(f"\nExample field with full traceability (tissue):")
                print(json.dumps(first_sample["biological_metadata"]["tissue"], indent=2))

    print(f"\n{'='*80}\n")


if __name__ == "__main__":
    main()
