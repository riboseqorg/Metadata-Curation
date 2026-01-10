#!/usr/bin/env python3
"""Demo: Project-level analysis with comprehensive LLM context.

This demonstrates giving the LLM full context about ALL samples in a project
to understand experimental design, replicates, and relationships.
"""

import sys
import json
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from omics_extractor.extraction.builder import build_project_metadata
from omics_extractor.extraction.llm_providers import create_provider
from omics_extractor.extraction.project_analyzer import analyze_project_design


def demo_project_analysis(project_id: str):
    """
    Demonstrate comprehensive project-level analysis.

    This gives the LLM:
    - ALL sample metadata (not just one sample)
    - ALL structured fields (organism, strain, etc.)
    - Study context
    - Run information (library strategy, etc.)

    The LLM then:
    - Identifies experimental variables
    - Groups replicates
    - Maps sample relationships
    - Characterizes experimental design
    """

    print(f"\n{'='*80}")
    print(f"PROJECT-LEVEL ANALYSIS: {project_id}")
    print(f"{'='*80}\n")

    # Phase 1: Get all metadata
    print("Phase 1: Extracting baseline metadata...")
    metadata = build_project_metadata(project_id)
    study = metadata["study"]
    samples = metadata["samples"]

    print(f"  ✓ Extracted {len(samples)} samples")
    print(f"  ✓ Study: {study.title.value if study.title else 'N/A'}")
    print()

    # Phase 2: LLM project analysis
    print("Phase 2: LLM project-level analysis...")
    print("  (Giving LLM ALL samples + ALL metadata fields)")
    print()

    import os
    if not os.environ.get("ANTHROPIC_API_KEY"):
        print("  ⚠️  ANTHROPIC_API_KEY not set - skipping LLM analysis")
        print("  Set it to run project analysis:")
        print("    export ANTHROPIC_API_KEY=your-key")
        return

    try:
        # Create provider
        provider = create_provider("claude", model="claude-sonnet-4-5-20250929")

        print(f"  Analyzing with {provider.get_model_name()}...")
        print(f"  (This may take 10-20 seconds for comprehensive analysis)")
        print()

        # Analyze project
        analysis = analyze_project_design(study, samples, provider)

        # Display results
        print("="*80)
        print("ANALYSIS RESULTS")
        print("="*80)
        print()

        print("Experimental Variables Identified:")
        if analysis.experimental_variables:
            for var in analysis.experimental_variables:
                print(f"  - {var}")
        else:
            print("  (none identified)")
        print()

        print("Replicate Groups:")
        if analysis.replicate_groups:
            for group in analysis.replicate_groups:
                print(f"  {group.group_id} ({group.replicate_type} replicates):")
                print(f"    Samples: {', '.join(group.sample_ids)}")
                print(f"    Conditions: {json.dumps(group.conditions, indent=6)}")
                print()
        else:
            print("  (no clear replicates identified)")
        print()

        print("Sample Relationships:")
        if analysis.relationships:
            for rel in analysis.relationships:
                print(f"  {rel.sample_id_1} ↔ {rel.sample_id_2}")
                print(f"    Type: {rel.relationship_type}")
                if rel.description:
                    print(f"    Description: {rel.description}")
                print()
        else:
            print("  (no relationships identified)")
        print()

        print("Experimental Design Summary:")
        print(f"  {analysis.design_summary}")
        print()

        print(f"Confidence: {analysis.confidence:.2f}")
        print()

    except Exception as e:
        print(f"  ✗ Analysis failed: {e}")
        import traceback
        traceback.print_exc()


def main():
    """Run project analysis demo."""

    print("\n" + "="*80)
    print("COMPREHENSIVE PROJECT ANALYSIS DEMO")
    print("="*80)
    print()
    print("This demo shows how to give the LLM comprehensive context:")
    print("  - ALL samples in the project (not just one)")
    print("  - ALL metadata fields (structured + unstructured)")
    print("  - Study description and context")
    print()
    print("The LLM then analyzes:")
    print("  - Experimental design and variables")
    print("  - Replicate groups (biological/technical)")
    print("  - Sample relationships (paired assays, time series, etc.)")
    print()

    # Test projects
    test_projects = [
        "PRJNA1170270",  # Ribo-seq mouse spermatocytes
    ]

    for project_id in test_projects:
        demo_project_analysis(project_id)

    print("="*80)
    print("KEY INSIGHTS")
    print("="*80)
    print()
    print("By giving the LLM ALL the context:")
    print("  ✓ Can identify experimental design patterns")
    print("  ✓ Can group biological/technical replicates")
    print("  ✓ Can find paired experiments (Ribo-Seq ↔ RNA-Seq)")
    print("  ✓ Can understand relationships between samples")
    print()
    print("This is MORE than just extracting tissue/cell_type per sample.")
    print("It's understanding the ENTIRE experimental design!")
    print()


if __name__ == "__main__":
    main()
