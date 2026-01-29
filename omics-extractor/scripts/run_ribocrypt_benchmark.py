import json
from pathlib import Path
import sys
import os

# Add src to path
sys.path.append(str(Path(__file__).parent.parent / "src"))

from omics_extractor.evaluation.benchmark import load_gold_standard_dataset, evaluate_model, compare_models
from omics_extractor.extraction.llm_providers import create_provider

DATA_PATH = Path(__file__).parent.parent / "tests/data/ribocrypt_gold.json"
RESULTS_PATH = Path(__file__).parent.parent / "ribocrypt_benchmark_results.json"

def run_benchmark():
    if not DATA_PATH.exists():
        print(f"Error: Dataset not found at {DATA_PATH}. Run scripts/prepare_ribocrypt_benchmark.py first.")
        return

    print(f"Loading benchmark dataset from {DATA_PATH}...")
    test_dataset = load_gold_standard_dataset(DATA_PATH)
    print(f"Loaded {len(test_dataset)} samples.")

    # Initialize providers
    # We'll use Claude 3.5 Sonnet as the baseline
    providers = []
    try:
        # Using Haiku for the test run as suggested (faster and cheaper)
        haiku = create_provider("claude", model="claude-3-5-haiku-latest")
        providers.append(haiku)
    except Exception as e:
        print(f"Warning: Could not initialize Claude provider: {e}")

    if not providers:
        print("Error: No providers initialized. Check your API keys.")
        return

    # Run evaluation
    results = []
    for provider in providers:
        result = evaluate_model(
            provider=provider,
            test_dataset=test_dataset,
            fields_to_evaluate=["organism", "tissue", "cell_line"],
            verbose=True
        )
        results.append(result)

    # Compare and save
    compare_models(results, output_path=RESULTS_PATH)
    print(f"Full benchmark complete. Results saved to {RESULTS_PATH}")

if __name__ == "__main__":
    run_benchmark()
