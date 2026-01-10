# LLM Model Comparison and Benchmarking

This guide explains how to evaluate and compare different LLM models for metadata extraction.

## Overview

The system supports **pluggable LLM backends**, allowing you to:
1. Test multiple models (Claude API, Llama, Mixtral, Qwen, etc.)
2. Benchmark them on gold-standard datasets
3. Compare accuracy, latency, and cost
4. Choose the best model for your use case

## Supported Model Providers

### 1. Claude API (Anthropic)
- **Models**: Claude 3.5 Sonnet, Claude 3 Sonnet, etc.
- **Pros**: Best quality, no GPU needed
- **Cons**: API costs, rate limits
- **Use case**: Baseline for comparison, small-scale production

### 2. VLLM (Recommended for A100)
- **Models**: Llama 3.1/3.3 70B, Mixtral 8x22B, Qwen 2.5 72B
- **Pros**: Fast batching, zero API cost, on-premise
- **Cons**: Requires GPU
- **Use case**: Large-scale production on A100 clusters

### 3. Transformers (Fallback)
- **Models**: Any HuggingFace model
- **Pros**: Works on CPU, flexible
- **Cons**: Slower than VLLM
- **Use case**: Development, small models

## Quick Start: Benchmark Models

### Step 1: Create Gold-Standard Test Dataset

```bash
# Copy the example template
cp config/gold_standard.example.json config/gold_standard.json

# Edit to add 20-50 manually curated samples
nano config/gold_standard.json
```

**What makes a good test dataset:**
- **Diverse**: Different organisms, tissues, treatments
- **Representative**: Matches your actual data
- **Carefully curated**: Read full publications to get ground truth
- **Challenging**: Include edge cases and ambiguous samples

See `config/gold_standard.example.json` for format.

### Step 2: Configure Models to Test

```bash
# Copy model config template
cp config/model_config.example.json config/model_config.json

# Edit to specify which models to benchmark
nano config/model_config.json
```

Example configuration:

```json
{
  "models_to_benchmark": [
    {
      "name": "Claude 3.5 Sonnet (Baseline)",
      "provider": "claude",
      "config": {
        "model": "claude-3-5-sonnet-20241022"
      }
    },
    {
      "name": "Llama 3.3 70B",
      "provider": "vllm",
      "config": {
        "model_path": "meta-llama/Llama-3.3-70B-Instruct",
        "tensor_parallel_size": 1
      }
    },
    {
      "name": "Qwen 2.5 72B",
      "provider": "vllm",
      "config": {
        "model_path": "Qwen/Qwen2.5-72B-Instruct"
      }
    }
  ]
}
```

### Step 3: Run Benchmark

```python
from omics_extractor.extraction.llm_providers import create_provider
from omics_extractor.evaluation.benchmark import (
    load_gold_standard_dataset,
    evaluate_model,
    compare_models,
)

# Load test dataset
test_data = load_gold_standard_dataset("config/gold_standard.json")

# Test Claude (baseline)
claude_provider = create_provider("claude", model="claude-3-5-sonnet-20241022")
claude_results = evaluate_model(claude_provider, test_data, verbose=True)

# Test Llama 3.3 70B (local on A100)
llama_provider = create_provider("vllm", model_path="meta-llama/Llama-3.3-70B-Instruct")
llama_results = evaluate_model(llama_provider, test_data, verbose=True)

# Compare
compare_models([claude_results, llama_results], output_path="benchmark_results.json")
```

## Understanding Metrics

### Accuracy Metrics

**Precision**: Of the fields the model extracted, what % were correct?
- High precision = Low false positives (rarely hallucinates)

**Recall**: Of the fields that exist, what % did the model extract?
- High recall = Low false negatives (rarely misses fields)

**F1 Score**: Harmonic mean of precision and recall
- Overall quality metric (best single metric)

**Exact Match Accuracy**: % of samples where ALL fields were correct
- Strictest metric (hardest to achieve)

### Performance Metrics

**Latency**: Time per sample (milliseconds)
- Lower = faster processing

**Tokens Used**: Total tokens consumed
- Affects API cost

**Estimated Cost**: Based on model pricing
- $0 for local models
- $X per 1M tokens for APIs

## Interpreting Results

Example output:

```
MODEL COMPARISON
================================================================================

Model                                    Accuracy   F1         Latency     Cost
--------------------------------------------------------------------------------
Llama-3.3-70B-Instruct                    82.5%     88.3%       450ms    $0.00
claude-3-5-sonnet-20241022                 85.0%     90.1%      1200ms    $0.45
Qwen-2.5-72B-Instruct                      80.0%     86.5%       420ms    $0.00

RECOMMENDATION:
  Best Model: claude-3-5-sonnet-20241022
  - F1 Score: 90.1%
  - Accuracy: 85.0%
  - Avg Latency: 1200ms
  - Cost: $0.45
```

**Decision Matrix:**

| Scenario | Recommended Model |
|----------|-------------------|
| Best accuracy at any cost | Claude 3.5 Sonnet |
| Best accuracy, zero cost | Llama 3.3 70B (if within 3% of Claude) |
| Fastest processing | Qwen 2.5 72B (usually faster than Llama) |
| Limited GPU memory | Llama 3.1 8B quantized |
| Hybrid (quality + cost) | Use Claude for low-confidence, Llama for high-confidence |

## Recommended Models for A100

Based on community benchmarks (your results may vary):

### Tier 1: Best Quality
1. **Claude 3.5 Sonnet** (API) - Baseline to beat
2. **Llama 3.3 70B** (VLLM) - Close to Claude, zero cost
3. **Qwen 2.5 72B** (VLLM) - Excellent for structured extraction

### Tier 2: Good Quality, Faster
4. **Mixtral 8x22B** (VLLM) - Fast, good quality
5. **Llama 3.1 70B** (VLLM) - Slightly older but proven

### Tier 3: Fast/Cheap
6. **Qwen 2.5 32B** (VLLM) - Great quality-to-size ratio
7. **Llama 3.1 8B** (quantized) - For limited resources

## Production Deployment

### Strategy 1: Pure Local (Recommended for A100)

```python
# Use best local model from benchmarks
provider = create_provider("vllm", model_path="meta-llama/Llama-3.3-70B-Instruct")

# Batch enrichment for efficiency
omics-extract batch-enrich projects/*.json \
  --output-dir enriched/ \
  --model llama-3.3-70b \
  --workers 8
```

**Pros:**
- Zero API costs
- No rate limits
- Data stays on-premise
- Batching = maximum GPU utilization

**Cons:**
- Need to download models (~140GB for 70B)
- Slightly lower quality than Claude (typically 2-5% F1)

### Strategy 2: Hybrid (Best of Both Worlds)

```python
# Use local model for bulk
local_provider = create_provider("vllm", model_path="meta-llama/Llama-3.3-70B-Instruct")

# Use Claude API for low-confidence cases
claude_provider = create_provider("claude")

# Extract with local model
result = local_provider.extract(prompt)

# If confidence < threshold, retry with Claude
if result.confidence.get("tissue", 1.0) < 0.7:
    result = claude_provider.extract(prompt)
```

**Pros:**
- 90% cost savings (most samples use local)
- Best accuracy on hard cases (Claude backup)
- Flexible scaling

## Creating Gold-Standard Datasets

### Best Practices

1. **Start small**: 20-30 samples for initial validation
2. **Diverse coverage**:
   - Multiple organisms (mouse, human, yeast, fly, etc.)
   - Various tissues (brain, liver, blood, cell lines, etc.)
   - Different assays (RNA-seq, Ribo-seq, ChIP-seq, etc.)
   - Edge cases (organoids, primary cells, treated samples)

3. **Careful curation**:
   - Read the full publication (not just title)
   - Check supplementary materials
   - Look at all metadata sources (BioSample + GEO + paper)
   - When unsure, mark as `null` (don't guess)

4. **Include challenging samples**:
   - Ambiguous descriptions
   - Missing information
   - Complex treatments
   - Multi-tissue studies

### Example Curation Process

For sample `SAMN12345678`:

1. **Read BioProject**: "Ribo-seq of mouse liver under fasting"
2. **Check BioSample attributes**:
   - `organism: Mus musculus`
   - `tissue: liver`
   - `strain: C57BL/6`
   - `treatment: fasted 24h`
3. **Read publication**: Confirms males, age 8-12 weeks
4. **Ground truth**:
   ```json
   {
     "tissue": "liver",
     "cell_type": null,
     "strain": "C57BL/6",
     "treatment": "fasting",
     "sex": "male",
     "age": "adult"
   }
   ```

## Troubleshooting

### VLLM Out of Memory

```
OutOfMemoryError: CUDA out of memory
```

**Solutions:**
- Reduce `gpu_memory_utilization` (try 0.8 instead of 0.9)
- Use smaller model (70B ’ 32B or 8B)
- Enable tensor parallelism across multiple GPUs
- Use quantization (8-bit or 4-bit)

### Low Model Accuracy

```
F1 Score: 45% (expected >80%)
```

**Possible causes:**
1. **Model too small**: Try 70B instead of 8B
2. **Prompt format mismatch**: Check if model expects specific format
3. **Temperature too high**: Use 0.0 for deterministic extraction
4. **Test dataset mismatch**: Ensure test data matches production data

### Slow Inference

```
Avg Latency: 5000ms (expected <1000ms)
```

**Solutions:**
- Use VLLM instead of Transformers (5-10x speedup)
- Enable batching (process multiple samples together)
- Check GPU utilization (should be >80%)
- Reduce `max_tokens` if generating too much

## Future Enhancements

1. **Automated test set expansion**: Use high-confidence predictions as pseudo-labels
2. **Active learning**: Focus curation effort on hardest samples
3. **Fine-tuning**: Train models specifically for metadata extraction
4. **Multi-task learning**: Train single model for extraction + normalization + validation
