#!/usr/bin/env python3
"""Quick benchmark: Compare baseline vs LLM enrichment on real RiboSeq data.

This script:
1. Fetches metadata from a few RiboSeq projects
2. Extracts with baseline (structured + ontology)
3. Enriches with LLM (Claude or local models)
4. Compares results and shows agreement
"""

import sys
import json
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from omics_extractor.extraction.builder import build_project_metadata
from omics_extractor.extraction.llm_providers import create_provider
from omics_extractor.extraction.llm_extractor import build_extraction_prompt

# Sample projects to test
TEST_PROJECTS = [
    "PRJNA1170270",  # Should have good metadata
]

def compare_extraction_methods(project_id: str):
    """Compare baseline vs LLM extraction."""
    
    print(f"\n{'='*80}")
    print(f"BENCHMARKING: {project_id}")
    print(f"{'='*80}\n")
    
    # Phase 1: Baseline extraction (structured + ontology)
    print("Phase 1: Baseline extraction (structured fields + ontology mapping)...")
    try:
        metadata = build_project_metadata(project_id)
        study = metadata["study"]
        samples = metadata["samples"]
        
        print(f"  ✓ Extracted {len(samples)} samples")
        
        # Show first sample baseline results
        if samples:
            first_sample_id = list(samples.keys())[0]
            first_sample = samples[first_sample_id]
            
            print(f"\n  Sample: {first_sample_id}")
            print(f"  Baseline extraction:")
            for field in ["tissue", "cell_type", "cell_line", "treatment", "strain"]:
                value = getattr(first_sample, field, None)
                if value:
                    conf = value.confidence if hasattr(value, 'confidence') else 'N/A'
                    ont = value.ontology_term if hasattr(value, 'ontology_term') else None
                    print(f"    {field:15s}: {value.value:20s} (conf: {conf:.2f}, ont: {ont})")
                else:
                    print(f"    {field:15s}: [MISSING]")
        
    except Exception as e:
        print(f"  ✗ Failed: {e}")
        return
    
    # Phase 2: LLM enrichment
    print(f"\nPhase 2: LLM enrichment (filling missing fields)...")
    
    # Check if we have API key
    import os
    if not os.environ.get("ANTHROPIC_API_KEY"):
        print("  ⚠️  ANTHROPIC_API_KEY not set - skipping LLM comparison")
        print("  Set it to compare with Claude:")
        print("    export ANTHROPIC_API_KEY=your-key")
        return
    
    try:
        # Create Claude provider (using Claude 4.5 Sonnet)
        provider = create_provider("claude", model="claude-sonnet-4-5-20250929")
        
        # Test on first sample
        if samples:
            first_sample_id = list(samples.keys())[0]
            first_sample = samples[first_sample_id]
            
            # Build prompt with baseline findings
            study_title = study.title.value if study.title else ""
            study_desc = study.description.value if study.description else ""
            sample_title = first_sample.sample_title.value if first_sample.sample_title else None
            sample_desc = first_sample.sample_description.value if first_sample.sample_description else None

            # Pass baseline findings to LLM so it knows what's already extracted
            existing = {}
            if first_sample.tissue:
                existing["tissue"] = first_sample.tissue.value
            if first_sample.cell_type:
                existing["cell_type"] = first_sample.cell_type.value
            if first_sample.cell_line:
                existing["cell_line"] = first_sample.cell_line.value
            if first_sample.treatment:
                existing["treatment"] = first_sample.treatment.value
            if first_sample.strain:
                existing["strain"] = first_sample.strain.value

            prompt = build_extraction_prompt(
                study_title=study_title,
                study_description=study_desc,
                sample_title=sample_title,
                sample_description=sample_desc,
                existing_metadata=existing,
            )
            
            print(f"  Extracting with {provider.get_model_name()}...")
            result = provider.extract(prompt)
            
            print(f"\n  LLM extraction:")
            for field in ["tissue", "cell_type", "cell_line", "treatment", "strain"]:
                value = getattr(result, field, None)
                conf = result.confidence.get(field, 0.0) if result.confidence else 0.0
                if value:
                    print(f"    {field:15s}: {value:20s} (conf: {conf:.2f})")
                else:
                    print(f"    {field:15s}: [MISSING]")
            
            # Compare
            print(f"\n  COMPARISON:")
            print(f"  {'Field':<15s} {'Baseline':<20s} {'LLM':<20s} {'Agreement'}")
            print(f"  {'-'*70}")
            
            for field in ["tissue", "cell_type", "cell_line", "treatment", "strain"]:
                baseline_val = getattr(first_sample, field, None)
                baseline_str = baseline_val.value if baseline_val else "[MISSING]"
                
                llm_val = getattr(result, field, None) 
                llm_str = llm_val if llm_val else "[MISSING]"
                
                # Check agreement
                if baseline_str == llm_str:
                    agreement = "✓ MATCH"
                elif baseline_str == "[MISSING]" and llm_str != "[MISSING]":
                    agreement = "+ LLM ADDED"
                elif baseline_str != "[MISSING]" and llm_str == "[MISSING]":
                    agreement = "! LLM MISSED"
                else:
                    agreement = "✗ DIFFER"
                
                print(f"  {field:<15s} {baseline_str:<20s} {llm_str:<20s} {agreement}")
            
            print(f"\n  Performance:")
            print(f"    Latency: {result.latency_ms:.0f}ms")
            print(f"    Tokens: {result.tokens_used}")
            
    except Exception as e:
        print(f"  ✗ LLM extraction failed: {e}")
        import traceback
        traceback.print_exc()

def main():
    """Run quick benchmark."""
    print("\n" + "="*80)
    print("QUICK BENCHMARK: Baseline vs LLM Metadata Extraction")
    print("="*80)
    print("\nThis compares:")
    print("  1. Baseline: Structured fields + ontology mapping")
    print("  2. LLM: Claude 3.5 Sonnet filling missing fields")
    print("\nModels to test (when running on A100):")
    print("  - Llama 3.3 70B (via VLLM)")
    print("  - Qwen 2.5 72B (via VLLM)")
    print("  - Mixtral 8x22B (via VLLM)")
    
    for project_id in TEST_PROJECTS:
        compare_extraction_methods(project_id)
    
    print("\n" + "="*80)
    print("NEXT STEPS:")
    print("="*80)
    print("\n1. On your Mac (development):")
    print("   - Run with Claude API to establish baseline")
    print("   - export ANTHROPIC_API_KEY=your-key")
    print("   - python scripts/quick_benchmark.py")
    print("\n2. On A100 (production):")
    print("   - Install VLLM: pip install vllm")
    print("   - Test local models:")
    print("     * meta-llama/Llama-3.3-70B-Instruct")
    print("     * Qwen/Qwen2.5-72B-Instruct")
    print("     * mistralai/Mixtral-8x22B-Instruct-v0.1")
    print("\n3. Compare results:")
    print("   - Which model matches Claude best?")
    print("   - Which is fastest?")
    print("   - Choose winner for production!")

if __name__ == "__main__":
    main()
