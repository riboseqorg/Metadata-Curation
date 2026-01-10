# Running the Metadata Extraction Benchmark

This benchmark compares different approaches to metadata extraction:
1. **Baseline**: Structured fields + ontology mapping (no LLM)
2. **Claude API**: LLM enrichment filling missing fields
3. **Local Models** (on A100): Llama, Qwen, Mixtral

## Quick Start (2 steps)

### Step 1: Set your Anthropic API key

```bash
export ANTHROPIC_API_KEY=sk-ant-your-key-here
```

Get your key from: https://console.anthropic.com/settings/keys

### Step 2: Run the benchmark

```bash
./RUN_BENCHMARK.sh
```

That's it! The script will:
- Extract metadata from real RiboSeq projects using baseline
- Enrich with Claude 3.5 Sonnet
- Show side-by-side comparison
- Report agreement and differences

## What You'll See

```
BENCHMARKING: PRJNA1170270
Phase 1: Baseline extraction (structured fields + ontology mapping)...
  ✓ Extracted 1 samples

  Sample: SRS22569018
  Baseline extraction:
    tissue         : testis               (conf: 1.00, ont: UBERON:0000473)
    cell_type      : [MISSING]
    cell_line      : [MISSING]
    treatment      : [MISSING]
    strain         : C57BL/6              (conf: 1.00, ont: None)

Phase 2: LLM enrichment (filling missing fields)...
  Extracting with claude-3-5-sonnet-20241022...

  LLM extraction:
    tissue         : testis               (conf: 0.95)
    cell_type      : spermatogenic cell   (conf: 0.85)
    cell_line      : [MISSING]
    treatment      : [MISSING]
    strain         : C57BL/6              (conf: 0.90)

  COMPARISON:
  Field           Baseline             LLM                  Agreement
  ----------------------------------------------------------------------
  tissue          testis               testis               ✓ MATCH
  cell_type       [MISSING]            spermatogenic cell   + LLM ADDED
  cell_line       [MISSING]            [MISSING]            ✓ MATCH
  treatment       [MISSING]            [MISSING]            ✓ MATCH
  strain          C57BL/6              C57BL/6              ✓ MATCH
```

## Understanding Results

- **✓ MATCH**: Both methods agree
- **+ LLM ADDED**: LLM filled a missing field
- **! LLM MISSED**: LLM missed a field baseline found
- **✗ DIFFER**: Methods disagree (investigate!)

## Next: Testing Local Models on A100

Once you've established Claude as your baseline, test local models:

```python
from omics_extractor.extraction.llm_providers import create_provider

# Test Llama 3.3 70B
llama_provider = create_provider("vllm", 
    model_path="meta-llama/Llama-3.3-70B-Instruct")

# Test Qwen 2.5 72B  
qwen_provider = create_provider("vllm",
    model_path="Qwen/Qwen2.5-72B-Instruct")

# Compare all models...
```

See `docs/MODEL_COMPARISON.md` for detailed guide.

## Troubleshooting

### "ANTHROPIC_API_KEY not set"
- Make sure you exported the key in your current terminal
- Check: `echo $ANTHROPIC_API_KEY` (should show your key)

### "Rate limit exceeded"
- You're hitting API limits
- Wait a minute and try again
- Or use smaller test set

### Want to test more projects?
Edit `scripts/quick_benchmark.py` and add more project IDs to `TEST_PROJECTS` list.
