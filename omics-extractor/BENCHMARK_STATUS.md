# Benchmark Status & Model Options

## ✅ What's Ready to Run RIGHT NOW (on your Mac)

### Anthropic Claude API
- **Model**: `claude-sonnet-4-5-20250929` (Claude 4.5 Sonnet, Sept 2025)
- **Works on**: Mac (API calls, no GPU needed)
- **Setup**: Just add API key and run
  ```bash
  export ANTHROPIC_API_KEY=sk-ant-your-key-here
  ./RUN_BENCHMARK.sh
  ```
- **Cost**: Pricing varies by model (check Anthropic pricing page)

## 🔧 Local Model Options

### For Mac (CPU-only, development testing)

**Small models that CAN run on Mac CPU:**
- **Llama 3.1 8B** - via Transformers
  - Pros: Free, runs on CPU, decent quality
  - Cons: Slower (~10-30s per sample), lower accuracy than 70B
  - Setup: `uv add transformers torch`

- **Qwen 2.5 7B** - via Transformers
  - Similar to Llama 8B, sometimes better at structured extraction

**Reality check for Mac:**
- 8B models will be SLOW on CPU (~10-30 seconds per sample)
- Quality will be lower than Claude or 70B models
- Good for testing the infrastructure, not production use

### For A100 (production)

**Large models that REQUIRE A100 GPU:**
- **Llama 3.3 70B** - via VLLM (recommended)
  - Pros: Fast (400-600ms), high quality, free
  - Cons: Needs A100 GPU (40GB+ VRAM)

- **Qwen 2.5 72B** - via VLLM
  - Excellent at structured extraction tasks

- **Mixtral 8x22B** - via VLLM
  - Very fast, good quality

**These CANNOT run on Mac** - they need:
- 40-80GB GPU memory (A100 has 40/80GB)
- CUDA GPU (Mac has CPU/Metal)
- VLLM library (optimized for NVIDIA GPUs)

## 🎯 Current Benchmark Script

**What it does now:**
1. ✓ Fetches real RiboSeq metadata
2. ✓ Extracts with baseline (structured + ontology)
3. ✓ Enriches with Claude API (if key set)
4. ✓ Shows side-by-side comparison

**What it does NOT do yet:**
- ✗ Test local models (infrastructure exists, not integrated)
- ✗ Compare multiple models (need to add loop)
- ✗ Run on Mac CPU (could add 8B support, but slow)

## 📊 Recommended Workflow

### Phase 1: Mac Development (NOW)
```bash
# Test with Claude API to validate the system
export ANTHROPIC_API_KEY=sk-ant-your-key
./RUN_BENCHMARK.sh

# This will:
# - Show baseline extraction works
# - Show Claude enrichment works
# - Establish quality baseline
```

### Phase 2: A100 Production (LATER)
```bash
# On A100 cluster, test local models
python scripts/model_comparison.py \
  --models llama-3.3-70b qwen-2.5-72b mixtral-8x22b \
  --test-projects PRJNA1170270 PRJNA1001014 \
  --baseline claude

# This will:
# - Compare 3 local models vs Claude baseline
# - Show which is fastest
# - Show which matches Claude quality
# - Choose winner for production
```

## 🤔 Your Question: "Aren't you going to run local models on the Mac?"

**Short answer**: The 70B models you want (Llama 3.3, Qwen 2.5) **cannot** run on Mac - they require A100-class GPUs.

**Options**:

1. **Just use Claude for now** (recommended)
   - Test on Mac with API
   - Move to local 70B models when you have A100 access

2. **Test small 8B models on Mac CPU**
   - I can add Llama 8B support to the benchmark
   - Will be SLOW but proves the infrastructure works
   - Not representative of production quality

3. **Wait for A100 access**
   - Skip local testing on Mac entirely
   - Go straight to testing 70B models on A100

**Which do you prefer?**

## 📝 Model Name Updates

✅ **Now Using**: `claude-sonnet-4-5-20250929` (Claude 4.5 Sonnet, Sept 2025)

**Previous attempts that returned 404 errors:**
- `claude-3-5-sonnet-20241022` (Oct 2024) - ❌ Not found
- `claude-3-5-sonnet-20240620` (June 2024) - ❌ Not found

**Why the change:**
Anthropic has moved to Claude 4.x family with a new naming convention:
- Old format: `claude-3-5-sonnet-YYYYMMDD`
- New format: `claude-sonnet-4-5-YYYYMMDD`

The model identifier was updated based on recent Anthropic documentation showing Claude 4.5 Sonnet is the current model.
