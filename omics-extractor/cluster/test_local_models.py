#!/usr/bin/env python3
"""Test local models on cluster with A100 GPUs.

This script:
1. Loads model configuration from config/model_paths.yaml
2. Tests each local model on real RiboSeq data
3. Compares results with Claude API baseline
4. Measures speed and quality

Usage:
    # Test all models
    python cluster/test_local_models.py

    # Test specific model
    python cluster/test_local_models.py --model llama-3.3-70b

    # Test with custom config
    python cluster/test_local_models.py --config /path/to/model_paths.yaml
"""

import argparse
import sys
import yaml
from pathlib import Path
from typing import Dict, Any
import time

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from omics_extractor.extraction.builder import build_project_metadata
from omics_extractor.extraction.llm_providers import create_provider
from omics_extractor.extraction.enhanced_extractor import needs_llm_enrichment, build_minimal_llm_prompt


def load_model_config(config_path: Path) -> Dict[str, Any]:
    """Load model configuration from YAML file."""
    if not config_path.exists():
        print(f"❌ Config file not found: {config_path}")
        print(f"   Copy config/model_paths.example.yaml to config/model_paths.yaml")
        print(f"   and update with your model paths")
        sys.exit(1)

    with open(config_path) as f:
        return yaml.safe_load(f)


def test_model(
    model_name: str,
    model_config: Dict[str, Any],
    project_id: str = "PRJNA1170270"
):
    """Test a single model on a RiboSeq project."""
    print(f"\n{'='*80}")
    print(f"Testing: {model_name}")
    print(f"{'='*80}\n")

    # Show config
    print(f"Configuration:")
    for key, value in model_config.items():
        print(f"  {key}: {value}")
    print()

    # Extract baseline metadata
    print(f"1. Extracting baseline metadata from {project_id}...")
    metadata = build_project_metadata(project_id)
    study = metadata["study"]
    samples = metadata["samples"]

    print(f"   ✓ Extracted {len(samples)} samples")
    print()

    # Check which samples need LLM
    samples_needing_llm = [
        (sid, s) for sid, s in samples.items()
        if needs_llm_enrichment(s)
    ]

    print(f"2. Gap analysis:")
    print(f"   Total samples: {len(samples)}")
    print(f"   Need LLM: {len(samples_needing_llm)}")
    print()

    if not samples_needing_llm:
        print("   ✓ All samples complete - no LLM needed!")
        return

    # Create provider
    print(f"3. Creating {model_config['provider']} provider...")
    try:
        if model_config['provider'] == 'vllm':
            provider = create_provider(
                "vllm",
                model_path=model_config['path'],
                tensor_parallel_size=model_config.get('tensor_parallel_size', 1),
                gpu_memory_utilization=model_config.get('gpu_memory_utilization', 0.9),
            )
        elif model_config['provider'] == 'transformers':
            provider = create_provider(
                "transformers",
                model_path=model_config['path'],
                device=model_config.get('device', 'cuda'),
            )
        else:
            print(f"   ❌ Unknown provider: {model_config['provider']}")
            return

        print(f"   ✓ Provider created: {provider.get_model_name()}")
    except Exception as e:
        print(f"   ❌ Failed to create provider: {e}")
        return

    # Test on first sample needing enrichment
    sample_id, sample = samples_needing_llm[0]

    print(f"\n4. Testing on sample: {sample_id}")

    # Build minimal prompt
    prompt = build_minimal_llm_prompt(
        sample,
        study.title.value if study.title else "",
        study.description.value if study.description else "",
    )

    print(f"   Prompt size: ~{len(prompt.split())} words")
    print()

    # Extract with LLM
    print(f"5. Running LLM extraction...")
    start_time = time.time()

    try:
        result = provider.extract(prompt, max_tokens=200)
        elapsed = time.time() - start_time

        print(f"   ✓ Extraction complete!")
        print(f"   Time: {elapsed:.2f}s")
        if hasattr(result, 'tokens_used'):
            print(f"   Tokens: {result.tokens_used}")
        print()

        # Show results
        print(f"6. Results:")
        print(f"   Before LLM:")
        for field in ["tissue", "cell_type", "strain", "treatment"]:
            value = getattr(sample, field, None)
            if value and value.value:
                print(f"     {field}: {value.value}")
            else:
                print(f"     {field}: [MISSING]")

        print(f"\n   After LLM:")
        for field in ["tissue", "cell_type", "strain", "treatment"]:
            value = getattr(result, field, None)
            if value:
                print(f"     {field}: {value}")
            else:
                print(f"     {field}: [MISSING]")

    except Exception as e:
        print(f"   ❌ Extraction failed: {e}")
        import traceback
        traceback.print_exc()


def main():
    parser = argparse.ArgumentParser(description="Test local models on cluster")
    parser.add_argument(
        "--config",
        default="config/model_paths.yaml",
        help="Path to model configuration YAML"
    )
    parser.add_argument(
        "--model",
        help="Test specific model (default: test all)"
    )
    parser.add_argument(
        "--project",
        default="PRJNA1170270",
        help="BioProject to test on"
    )

    args = parser.parse_args()

    # Load config
    config_path = Path(args.config)
    config = load_model_config(config_path)

    # Get models to test
    if args.model:
        if args.model not in config['local_models']:
            print(f"❌ Model '{args.model}' not found in config")
            print(f"   Available models: {', '.join(config['local_models'].keys())}")
            sys.exit(1)
        models_to_test = {args.model: config['local_models'][args.model]}
    else:
        models_to_test = config['local_models']

    print("\n" + "="*80)
    print("LOCAL MODEL TESTING ON CLUSTER")
    print("="*80)
    print(f"\nConfig: {config_path}")
    print(f"Project: {args.project}")
    print(f"Models to test: {len(models_to_test)}")
    print()

    # Test each model
    for model_name, model_config in models_to_test.items():
        test_model(model_name, model_config, args.project)

    print("\n" + "="*80)
    print("✅ Testing Complete!")
    print("="*80)


if __name__ == "__main__":
    main()
