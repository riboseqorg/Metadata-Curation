#!/usr/bin/env python3
"""
Prepare batch enrichment jobs for SLURM.

Analyzes batch_summary.json and creates:
1. List of low-completeness projects needing enrichment
2. Optimal job array configuration
3. Estimated resource requirements
"""

import json
import sys
from pathlib import Path


def main():
    """Prepare enrichment batch configuration."""
    batch_summary = Path("riboseq_batch_output/batch_summary.json")

    if not batch_summary.exists():
        print(f"Error: {batch_summary} not found")
        print("Run batch extraction first:")
        print("  python scripts/batch_extract_riboseq_traceable.py")
        sys.exit(1)

    with open(batch_summary) as f:
        data = json.load(f)

    # Find projects needing enrichment (completeness < 50%)
    low_completeness = []
    medium_completeness = []
    high_completeness = []

    for project in data["project_summaries"]:
        if project["status"] != "success":
            continue

        completeness = project["completeness"]
        if completeness < 50:
            low_completeness.append(project)
        elif completeness < 90:
            medium_completeness.append(project)
        else:
            high_completeness.append(project)

    print("=" * 80)
    print("ENRICHMENT BATCH PREPARATION")
    print("=" * 80)

    print(f"\nProjects by completeness:")
    print(f"  Low (<50%):    {len(low_completeness)} projects - PRIORITY for enrichment")
    print(f"  Medium (50-90%): {len(medium_completeness)} projects - Optional enrichment")
    print(f"  High (>90%):   {len(high_completeness)} projects - Skip enrichment")

    # Calculate total samples needing enrichment
    total_samples = sum(p["samples"] for p in low_completeness)
    total_runs = sum(p["runs"] for p in low_completeness)

    print(f"\nLow-completeness projects contain:")
    print(f"  Total samples: {total_samples}")
    print(f"  Total runs: {total_runs}")

    # Estimate resources
    # Assume ~2 seconds per sample for LLM enrichment
    avg_samples_per_project = total_samples / len(low_completeness) if low_completeness else 0
    estimated_time_per_job = avg_samples_per_project * 2  # seconds
    estimated_time_minutes = estimated_time_per_job / 60

    print(f"\nEstimated resources per job:")
    print(f"  Avg samples/project: {avg_samples_per_project:.1f}")
    print(f"  Estimated time: {estimated_time_minutes:.1f} minutes")
    print(f"  Recommended time limit: 2 hours (safety margin)")

    # Save project lists
    output_dir = Path("riboseq_batch_output")

    # Low-completeness projects (priority)
    low_file = output_dir / "low_completeness_projects.txt"
    with open(low_file, "w") as f:
        for project in low_completeness:
            f.write(f"{project['project_id']}\n")

    print(f"\n✓ Created: {low_file}")
    print(f"  {len(low_completeness)} projects listed")

    # Medium-completeness projects (optional)
    medium_file = output_dir / "medium_completeness_projects.txt"
    with open(medium_file, "w") as f:
        for project in medium_completeness:
            f.write(f"{project['project_id']}\n")

    if medium_completeness:
        print(f"\n✓ Created: {medium_file}")
        print(f"  {len(medium_completeness)} projects listed (optional enrichment)")

    # Create enrichment config
    config = {
        "low_completeness": {
            "count": len(low_completeness),
            "total_samples": sum(p["samples"] for p in low_completeness),
            "projects": [p["project_id"] for p in low_completeness],
            "slurm_array": f"0-{len(low_completeness)-1}",
        },
        "medium_completeness": {
            "count": len(medium_completeness),
            "total_samples": sum(p["samples"] for p in medium_completeness),
            "projects": [p["project_id"] for p in medium_completeness],
            "slurm_array": f"0-{len(medium_completeness)-1}" if medium_completeness else None,
        },
        "resource_estimates": {
            "avg_samples_per_project": round(avg_samples_per_project, 1),
            "estimated_minutes_per_job": round(estimated_time_minutes, 1),
            "recommended_time_limit": "02:00:00",
            "recommended_memory": "32G",
            "recommended_gpus": 1,
        }
    }

    config_file = output_dir / "enrichment_config.json"
    with open(config_file, "w") as f:
        json.dump(config, f, indent=2)

    print(f"\n✓ Created: {config_file}")

    # Show SLURM command
    print("\n" + "=" * 80)
    print("READY TO SUBMIT")
    print("=" * 80)

    print("\nSubmit low-completeness projects (priority):")
    print(f"  sbatch --array=0-{len(low_completeness)-1} cluster/enrich_batch_slurm.sh")

    if medium_completeness:
        print("\nOptionally submit medium-completeness projects:")
        print(f"  # Edit cluster/enrich_batch_slurm.sh to use medium_completeness_projects.txt")
        print(f"  sbatch --array=0-{len(medium_completeness)-1} cluster/enrich_batch_slurm.sh")

    print("\nMonitor jobs:")
    print("  squeue -u $USER")
    print("  tail -f logs/enrich_*.out")

    print("\n" + "=" * 80)

    # Show example projects
    print("\nExample low-completeness projects:")
    for i, project in enumerate(low_completeness[:5], 1):
        print(f"  {i}. {project['project_id']}: {project['samples']} samples, {project['completeness']:.1f}% complete")

    if len(low_completeness) > 5:
        print(f"  ... and {len(low_completeness) - 5} more")


if __name__ == "__main__":
    main()
