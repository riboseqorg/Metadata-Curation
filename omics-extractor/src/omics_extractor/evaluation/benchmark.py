"""Benchmark and evaluate different LLM models for metadata extraction.

This module provides tools to:
1. Create gold-standard test datasets
2. Benchmark multiple models on the same data
3. Calculate accuracy metrics (precision, recall, F1)
4. Compare latency and cost
5. Generate comparison reports
"""

from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass, asdict
from pathlib import Path
import json
import time
from collections import defaultdict

from ..extraction.llm_providers import LLMProvider, LLMExtractionResult
from ..extraction.llm_extractor import build_extraction_prompt
from difflib import SequenceMatcher


def values_match(predicted: Optional[str], ground_truth: Optional[str], field: str = "") -> bool:
    """Check if two metadata values match with fuzzy logic and synonyms."""
    if not predicted or not ground_truth:
        return predicted == ground_truth

    p = str(predicted).lower().strip()
    g = str(ground_truth).lower().strip()

    if p == g:
        return True

    # 1. Organism synonyms
    if field == "strain" or field == "organism":
        organism_synonyms = {
            "yeast": ["saccharomyces cerevisiae", "s. cerevisiae", "schizosaccharomyces pombe", "s. pombe", "yeast"],
            "human": ["homo sapiens", "h. sapiens", "human"],
            "mouse": ["mus musculus", "m. musculus", "mouse"],
            "zebrafish": ["danio rerio", "d. rerio", "zebrafish"],
            "fruit fly": ["drosophila melanogaster", "d. melanogaster", "drosophila"],
            "arabidopsis": ["arabidopsis thaliana", "a. thaliana", "arabidopsis"],
        }
        for common, scientifics in organism_synonyms.items():
            if (p == common or p in scientifics) and (g == common or g in scientifics):
                return True
        
        # Substring match for organisms (e.g. "drosophila" matches "drosophila melanogaster")
        if p in g or g in p:
            if len(p) > 4 and len(g) > 4:
                return True

    # 2. Substring matching for tissues and cell types
    if field in ["tissue", "cell_type", "cell_line", "treatment", "strain"]:
        # Allow short matches for specific abbreviations
        common_abbreviations = ["wt", "tg", "ko", "kd", "oe", "mut"]
        if p in common_abbreviations or g in common_abbreviations:
            if p == g or p in g.split() or g in p.split():
                return True

        if p in g or g in p:
            # Check length to avoid trivial matches, but allow shorter for specific fields
            min_len = 3 if field != "strain" else 2
            if len(p) >= min_len and len(g) >= min_len:
                return True
        
        # Fuzzy matching using SequenceMatcher
        if SequenceMatcher(None, p, g).ratio() > 0.8:
            return True
        
        # Special case for CHX (cycloheximide)
        if ("chx" in p or "cycloheximide" in p) and ("chx" in g or "cycloheximide" in g):
            return True

    # 3. Fuzzy matching
    if SequenceMatcher(None, p, g).ratio() > 0.8:
        return True

    return False


@dataclass
class GoldStandardSample:
    """A sample with manually curated ground truth metadata."""

    sample_id: str
    study_title: str
    study_description: str
    sample_title: Optional[str] = None
    sample_description: Optional[str] = None
    characteristics: Optional[Dict[str, str]] = None
    run_metadata: Optional[Dict[str, Any]] = None
    abstract: Optional[str] = None

    # Ground truth (manually curated)
    ground_truth: Dict[str, Optional[str]] = None

    # Source information
    bioproject_id: Optional[str] = None
    biosample_id: Optional[str] = None

    def __post_init__(self):
        if self.ground_truth is None:
            self.ground_truth = {}


@dataclass
class ModelEvaluationResult:
    """Results from evaluating a model on test dataset."""

    model_name: str

    # Accuracy metrics per field
    field_metrics: Dict[str, Dict[str, float]] = None  # field -> {precision, recall, f1}

    # Overall metrics
    overall_precision: float = 0.0
    overall_recall: float = 0.0
    overall_f1: float = 0.0
    overall_accuracy: float = 0.0  # Exact match rate

    # Performance metrics
    avg_latency_ms: float = 0.0
    total_latency_ms: float = 0.0
    total_tokens: int = 0

    # Cost estimates (based on typical pricing)
    estimated_cost_usd: float = 0.0

    # Per-sample results
    sample_results: List[Dict] = None

    def __post_init__(self):
        if self.field_metrics is None:
            self.field_metrics = {}
        if self.sample_results is None:
            self.sample_results = []


def load_gold_standard_dataset(path: Path) -> List[GoldStandardSample]:
    """
    Load gold-standard test dataset from JSON.
    """
    with open(path) as f:
        data = json.load(f)

    samples = []
    for item in data["samples"]:
        sample = GoldStandardSample(
            sample_id=item["sample_id"],
            study_title=item["study_title"],
            study_description=item["study_description"],
            sample_title=item.get("sample_title"),
            sample_description=item.get("sample_description"),
            characteristics=item.get("characteristics"),
            run_metadata=item.get("run_metadata"),
            abstract=item.get("abstract"),
            ground_truth=item["ground_truth"],
            bioproject_id=item.get("bioproject_id"),
            biosample_id=item.get("biosample_id"),
        )
        samples.append(sample)

    return samples


def evaluate_model(
    provider: LLMProvider,
    test_dataset: List[GoldStandardSample],
    fields_to_evaluate: Optional[List[str]] = None,
    verbose: bool = False,
) -> ModelEvaluationResult:
    """
    Evaluate a model on gold-standard dataset.
    """
    if fields_to_evaluate is None:
        fields_to_evaluate = ["organism", "tissue", "cell_type", "cell_line", "treatment", "strain", "age", "sex"]

    if verbose:
        print(f"\nEvaluating {provider.get_model_name()} on {len(test_dataset)} samples...")
        print("=" * 60)

    # Track metrics
    field_counts = defaultdict(lambda: {"tp": 0, "fp": 0, "fn": 0, "tn": 0})
    total_exact_matches = 0
    total_latency = 0.0
    total_tokens = 0
    sample_results = []

    # Evaluate each sample
    for i, sample in enumerate(test_dataset, 1):
        if verbose:
            print(f"[{i}/{len(test_dataset)}] {sample.sample_id}...")

        # Build prompt
        prompt = build_extraction_prompt(
            study_title=sample.study_title,
            study_description=sample.study_description,
            sample_title=sample.sample_title,
            sample_description=sample.sample_description,
            characteristics=sample.characteristics,
            run_metadata=sample.run_metadata,
            abstract=sample.abstract,
        )

        # Extract
        try:
            result = provider.extract(prompt)
            total_latency += result.latency_ms or 0
            total_tokens += result.tokens_used or 0
        except Exception as e:
            if verbose:
                print(f"  ERROR: {e}")
            continue

        # Compare with ground truth
        sample_correct = True
        field_results = {}

        for field in fields_to_evaluate:
            predicted = getattr(result, field)
            ground_truth = sample.ground_truth.get(field)

            # Use fuzzy matching
            match = values_match(predicted, ground_truth, field)

            if predicted and ground_truth:
                if match:
                    field_counts[field]["tp"] += 1  # True positive
                    field_results[field] = f"MATCH ({predicted})"
                else:
                    field_counts[field]["fp"] += 1  # False positive
                    field_counts[field]["fn"] += 1  # Also missed the true value
                    field_results[field] = f"MISMATCH (pred: {predicted}, gt: {ground_truth})"
                    sample_correct = False
            elif predicted and not ground_truth:
                field_counts[field]["fp"] += 1  # False positive (hallucinated)
                field_results[field] = f"FP ({predicted})"
                sample_correct = False
            elif not predicted and ground_truth:
                field_counts[field]["fn"] += 1  # False negative (missed)
                field_results[field] = f"FN (missed {ground_truth})"
                sample_correct = False
            else:
                field_counts[field]["tn"] += 1  # True negative (correctly empty)
                field_results[field] = "(both empty)"

        if sample_correct:
            total_exact_matches += 1

        sample_results.append({
            "sample_id": sample.sample_id,
            "exact_match": sample_correct,
            "field_results": field_results,
            "latency_ms": result.latency_ms,
        })

        if verbose:
            status = "  EXACT MATCH" if sample_correct else "  PARTIAL"
            print(f"  {status}")

    # Calculate metrics per field
    field_metrics = {}
    for field, counts in field_counts.items():
        tp = counts["tp"]
        fp = counts["fp"]
        fn = counts["fn"]

        precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
        recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
        f1 = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0.0

        field_metrics[field] = {
            "precision": precision,
            "recall": recall,
            "f1": f1,
            "tp": tp,
            "fp": fp,
            "fn": fn,
        }

    # Overall metrics (micro-averaged)
    total_tp = sum(c["tp"] for c in field_counts.values())
    total_fp = sum(c["fp"] for c in field_counts.values())
    total_fn = sum(c["fn"] for c in field_counts.values())

    overall_precision = total_tp / (total_tp + total_fp) if (total_tp + total_fp) > 0 else 0.0
    overall_recall = total_tp / (total_tp + total_fn) if (total_tp + total_fn) > 0 else 0.0
    overall_f1 = 2 * (overall_precision * overall_recall) / (overall_precision + overall_recall) if (overall_precision + overall_recall) > 0 else 0.0
    overall_accuracy = total_exact_matches / len(test_dataset) if test_dataset else 0

    # Estimate cost
    estimated_cost = estimate_cost(
        model_name=provider.get_model_name(),
        total_tokens=total_tokens,
    )

    result = ModelEvaluationResult(
        model_name=provider.get_model_name(),
        field_metrics=field_metrics,
        overall_precision=overall_precision,
        overall_recall=overall_recall,
        overall_f1=overall_f1,
        overall_accuracy=overall_accuracy,
        avg_latency_ms=total_latency / len(test_dataset) if test_dataset else 0,
        total_latency_ms=total_latency,
        total_tokens=total_tokens,
        estimated_cost_usd=estimated_cost,
        sample_results=sample_results,
    )

    if verbose:
        print_evaluation_summary(result)

    return result


def estimate_cost(model_name: str, total_tokens: int) -> float:
    """
    Estimate cost based on model and token usage.
    """
    cost_per_1m_tokens = {
        "claude-3-5-sonnet": 9.0,
        "claude-3-sonnet": 6.0,
        "claude-3-5-haiku": 1.0,
        "gpt-4": 20.0,
        "gpt-3.5": 1.0,
    }

    # Check if local model
    if any(x in model_name.lower() for x in ["llama", "mixtral", "qwen", "mistral"]):
        return 0.0

    # Find matching pricing
    for key, price in cost_per_1m_tokens.items():
        if key in model_name.lower():
            return (total_tokens / 1_000_000) * price

    # Default estimate
    return (total_tokens / 1_000_000) * 10.0


def print_evaluation_summary(result: ModelEvaluationResult):
    """Print formatted evaluation summary."""
    print("\n" + "=" * 60)
    print(f"EVALUATION RESULTS: {result.model_name}")
    print("=" * 60)

    print(f"\nOVERALL METRICS:")
    print(f"  Exact Match Accuracy: {result.overall_accuracy:.1%}")
    print(f"  Precision: {result.overall_precision:.1%}")
    print(f"  Recall: {result.overall_recall:.1%}")
    print(f"  F1 Score: {result.overall_f1:.1%}")

    print(f"\nPER-FIELD METRICS:")
    for field, metrics in result.field_metrics.items():
        print(f"  {field:20s} - P: {metrics['precision']:.1%}, R: {metrics['recall']:.1%}, F1: {metrics['f1']:.1%}")

    print(f"\nPERFORMANCE:")
    print(f"  Avg Latency: {result.avg_latency_ms:.0f}ms")
    print(f"  Total Tokens: {result.total_tokens:,}")
    print(f"  Estimated Cost: ${result.estimated_cost_usd:.2f}")

    print("=" * 60 + "\n")


def compare_models(
    results: List[ModelEvaluationResult],
    output_path: Optional[Path] = None,
) -> None:
    """
    Compare multiple model evaluation results.
    """
    print("\n" + "=" * 80)
    print("MODEL COMPARISON")
    print("=" * 80)

    # Sort by F1 score
    results_sorted = sorted(results, key=lambda r: r.overall_f1, reverse=True)

    # Print table header
    print(f"\n{'Model':<40} {'Accuracy':<10} {'F1':<10} {'Latency':<12} {'Cost':<10}")
    print("-" * 80)

    for result in results_sorted:
        print(f"{result.model_name:<40} "
              f"{result.overall_accuracy:>8.1%}  "
              f"{result.overall_f1:>8.1%}  "
              f"{result.avg_latency_ms:>9.0f}ms  "
              f"${result.estimated_cost_usd:>8.2f}")

    print("\nRECOMMENDATION:")
    best = results_sorted[0]
    print(f"  Best Model: {best.model_name}")
    print(f"  - F1 Score: {best.overall_f1:.1%}")
    print(f"  - Accuracy: {best.overall_accuracy:.1%}")
    print(f"  - Avg Latency: {best.avg_latency_ms:.0f}ms")
    print(f"  - Cost: ${best.estimated_cost_usd:.2f}")

    # Save comparison
    if output_path:
        comparison_data = {
            "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
            "results": [asdict(r) for r in results_sorted],
        }
        with open(output_path, "w") as f:
            json.dump(comparison_data, f, indent=2)
        print(f"\nComparison saved to: {output_path}")

    print("=" * 80 + "\n")
