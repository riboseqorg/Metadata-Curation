#!/usr/bin/env python3
"""Demo: Two extraction modes - Lightweight vs Comprehensive.

This demonstrates:
1. Lightweight mode: Field mapping + ontologies (no LLM, for production scale)
2. Enriched mode: + minimal LLM gap-filling (for samples with gaps)
3. Comprehensive mode: Project-level analysis (run selectively on important studies)
"""

import sys
import os
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from omics_extractor.extraction.builder import build_project_metadata
from omics_extractor.extraction.enhanced_extractor import (
    needs_llm_enrichment,
    build_minimal_llm_prompt,
)
from omics_extractor.extraction.llm_providers import create_provider
from omics_extractor.extraction.project_analyzer import analyze_project_design


def demo_lightweight_mode(project_id: str):
    """
    Mode 1: Lightweight extraction (NO LLM).

    Uses:
    - Field mapping (field_mappings.py)
    - Ontology mapping
    - Direct structured field extraction

    Cost: $0
    Speed: Fast
    Use case: Production pipeline at scale
    """
    print(f"\n{'='*80}")
    print(f"MODE 1: LIGHTWEIGHT (No LLM)")
    print(f"{'='*80}\n")

    print(f"Extracting: {project_id}")
    metadata = build_project_metadata(project_id)
    study = metadata["study"]
    samples = metadata["samples"]

    print(f"  ✓ Extracted {len(samples)} samples")
    print(f"  ✓ Study: {study.title.value if study.title else 'N/A'}")
    print()

    # Show first sample
    if samples:
        sample_id = list(samples.keys())[0]
        sample = samples[sample_id]

        print(f"Sample: {sample_id}")
        print(f"  Baseline extraction (NO LLM):")

        fields = ["organism", "tissue", "cell_type", "cell_line", "strain", "treatment"]
        for field in fields:
            value = getattr(sample, field, None)
            if value and value.value:
                conf = value.confidence if hasattr(value, 'confidence') else 'N/A'
                ont = value.ontology_term if hasattr(value, 'ontology_term') else None
                print(f"    {field:15s}: {value.value:20s} (conf: {conf:.2f}, ont: {ont})")
            else:
                print(f"    {field:15s}: [MISSING]")

        # Check if LLM would help
        print()
        if needs_llm_enrichment(sample):
            print("  ⚠️  This sample has gaps - would benefit from LLM enrichment")
        else:
            print("  ✓ This sample is sufficiently complete - LLM not needed")


def demo_enriched_mode(project_id: str):
    """
    Mode 2: Enriched extraction (Minimal LLM gap-filling).

    Uses:
    - Everything from Mode 1
    - + Small LLM prompts for missing critical fields only

    Cost: ~$0.001-0.005 per sample (only for samples with gaps)
    Speed: Medium (LLM only when needed)
    Use case: Production with quality assurance
    """
    print(f"\n{'='*80}")
    print(f"MODE 2: ENRICHED (Minimal LLM for gaps)")
    print(f"{'='*80}\n")

    if not os.environ.get("ANTHROPIC_API_KEY"):
        print("  ⚠️  ANTHROPIC_API_KEY not set - skipping LLM enrichment")
        print("  This mode requires an API key to fill gaps with LLM")
        return

    print(f"Extracting: {project_id}")
    metadata = build_project_metadata(project_id)
    study = metadata["study"]
    samples = metadata["samples"]

    print(f"  ✓ Baseline extraction complete")
    print()

    # Check which samples need LLM
    samples_needing_llm = [
        (sid, s) for sid, s in samples.items()
        if needs_llm_enrichment(s)
    ]

    print(f"  Analysis:")
    print(f"    Total samples: {len(samples)}")
    print(f"    Complete without LLM: {len(samples) - len(samples_needing_llm)}")
    print(f"    Need LLM enrichment: {len(samples_needing_llm)}")
    print(f"    LLM usage rate: {len(samples_needing_llm) / len(samples) * 100:.1f}%")
    print()

    if samples_needing_llm:
        # Create provider
        provider = create_provider("claude")

        # Show first sample needing enrichment
        sample_id, sample = samples_needing_llm[0]

        print(f"Example: Enriching {sample_id}")
        print(f"  Building minimal prompt...")

        # Build minimal prompt
        prompt = build_minimal_llm_prompt(
            sample,
            study.title.value if study.title else "",
            study.description.value if study.description else "",
        )

        print(f"  Prompt size: ~{len(prompt.split())} words (~{len(prompt) // 4} tokens)")
        print()
        print(f"  Calling LLM...")

        try:
            result = provider.extract(prompt, max_tokens=200)

            print(f"  ✓ LLM enrichment complete")
            print(f"  Tokens used: ~{result.tokens_used if hasattr(result, 'tokens_used') else 'N/A'}")
            print(f"  Cost: ~${0.005:.4f}")

        except Exception as e:
            print(f"  ✗ LLM enrichment failed: {e}")

    else:
        print("  ✓ All samples complete - no LLM needed!")


def demo_comprehensive_mode(project_id: str):
    """
    Mode 3: Comprehensive project analysis (Large context).

    Uses:
    - Project-level LLM analysis
    - ALL samples + ALL metadata fields
    - Experimental design characterization
    - Replicate identification
    - Sample relationship mapping

    Cost: ~$0.10-0.20 per project (run once per project)
    Speed: Slower (comprehensive analysis)
    Use case: Important studies, quality datasets, experimental design analysis
    """
    print(f"\n{'='*80}")
    print(f"MODE 3: COMPREHENSIVE (Large context project analysis)")
    print(f"{'='*80}\n")

    if not os.environ.get("ANTHROPIC_API_KEY"):
        print("  ⚠️  ANTHROPIC_API_KEY not set - skipping comprehensive analysis")
        print("  This mode requires an API key for project-level analysis")
        return

    print(f"Analyzing: {project_id}")
    print(f"  This mode gives LLM ALL samples + ALL metadata for deep analysis")
    print()

    metadata = build_project_metadata(project_id)
    study = metadata["study"]
    samples = metadata["samples"]

    print(f"  ✓ Extracted {len(samples)} samples")
    print(f"  ✓ Study: {study.title.value if study.title else 'N/A'}")
    print()

    print(f"  Running comprehensive project analysis...")
    print(f"  (This uses large context - may take 15-30 seconds)")
    print()

    try:
        provider = create_provider("claude")

        analysis = analyze_project_design(study, samples, provider)

        print("="*80)
        print("COMPREHENSIVE ANALYSIS RESULTS")
        print("="*80)
        print()

        print("Experimental Variables:")
        if analysis.experimental_variables:
            for var in analysis.experimental_variables:
                print(f"  - {var}")
        else:
            print("  (none identified)")
        print()

        print("Replicate Groups:")
        if analysis.replicate_groups:
            for group in analysis.replicate_groups:
                print(f"  {group.group_id} ({group.replicate_type}):")
                print(f"    Samples: {', '.join(group.sample_ids)}")
        else:
            print("  (no replicates identified)")
        print()

        print("Sample Relationships:")
        if analysis.relationships:
            for rel in analysis.relationships:
                print(f"  {rel.sample_id_1} ↔ {rel.sample_id_2} ({rel.relationship_type})")
        else:
            print("  (no relationships identified)")
        print()

        print("Design Summary:")
        print(f"  {analysis.design_summary}")
        print()

        print(f"Confidence: {analysis.confidence:.2f}")
        print(f"Cost: ~$0.15 (estimated)")

    except Exception as e:
        print(f"  ✗ Analysis failed: {e}")
        import traceback
        traceback.print_exc()


def main():
    """Run all three modes."""

    print("\n" + "="*80)
    print("EXTRACTION MODES DEMONSTRATION")
    print("="*80)
    print()
    print("This demo shows three extraction approaches:")
    print()
    print("1. LIGHTWEIGHT: Field mapping + ontologies (NO LLM)")
    print("   → Fast, $0 cost, good for production at scale")
    print()
    print("2. ENRICHED: + Minimal LLM gap-filling")
    print("   → Medium speed, ~$0.001-0.005/sample, only when needed")
    print()
    print("3. COMPREHENSIVE: Project-level analysis (LARGE context)")
    print("   → Slower, ~$0.10-0.20/project, run selectively for important studies")
    print()

    test_project = "PRJNA1170270"

    # Mode 1: Lightweight (no LLM)
    demo_lightweight_mode(test_project)

    # Mode 2: Enriched (minimal LLM)
    demo_enriched_mode(test_project)

    # Mode 3: Comprehensive (large context)
    demo_comprehensive_mode(test_project)

    print("\n" + "="*80)
    print("SUMMARY")
    print("="*80)
    print()
    print("Production Pipeline Recommendation:")
    print("  1. Run Mode 1 (lightweight) on ALL samples")
    print("  2. Run Mode 2 (enriched) only on samples with gaps")
    print("  3. Run Mode 3 (comprehensive) only on important/curated studies")
    print()
    print("This gives you:")
    print("  ✓ Fast extraction for 70-80% of samples (no LLM)")
    print("  ✓ Quality assurance for 20-30% with gaps (minimal LLM)")
    print("  ✓ Deep analysis for selected studies (comprehensive LLM)")
    print("  ✓ Average cost: ~$0.001-0.002 per sample")
    print()


if __name__ == "__main__":
    main()
