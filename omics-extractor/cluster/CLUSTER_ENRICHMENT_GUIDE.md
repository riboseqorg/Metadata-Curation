# Cluster Enrichment Guide

## Overview

After extracting metadata from 100 RiboSeq projects, 61 have low completeness (<50%) and need LLM enrichment. This guide shows how to run enrichment on the GPU cluster using local 70B models.

## Two Approaches

### 1. Batch Job Array (Recommended)

**Best for:** Production runs, processing all 61 projects

**Advantages:**
- Parallel execution across multiple GPUs
- Fault tolerant - failed jobs don't affect others
- Better cluster scheduling
- Resume capability
- Automatic logging per project

**Run time:** ~2-4 hours (depending on GPU availability)

### 2. Interactive Session

**Best for:** Testing, debugging, manual processing

**Advantages:**
- Direct interaction for debugging
- See output in real-time
- Good for testing model loading
- Process 1-3 projects manually

**Run time:** Manual processing, slower

---

## Quick Start (Batch Approach)

### Step 1: Prepare Enrichment Batch

```bash
cd /path/to/Metadata-Curation/omics-extractor

# Analyze batch results and create project lists
python cluster/prepare_enrichment_batch.py
```

**Output:**
- `riboseq_batch_output/low_completeness_projects.txt` - 61 projects needing enrichment
- `riboseq_batch_output/enrichment_config.json` - Job configuration
- Resource estimates and SLURM commands

### Step 2: Configure Model Paths

```bash
# Copy example config
cp config/model_paths.example.yaml config/model_paths.yaml

# Edit with your model paths
nano config/model_paths.yaml
```

Example:
```yaml
local_models:
  llama-3.3-70b:
    path: "/data/models/Llama-3.3-70B-Instruct"
    provider: "vllm"
    tensor_parallel_size: 2
    gpu_memory_utilization: 0.9

default: "llama-3.3-70b"
```

### Step 3: Setup Log Directory

```bash
mkdir -p logs
```

### Step 4: Submit Job Array

```bash
# Edit cluster/enrich_batch_slurm.sh to set WORK_DIR
nano cluster/enrich_batch_slurm.sh

# Submit all 61 low-completeness projects
sbatch --array=0-60 cluster/enrich_batch_slurm.sh
```

### Step 5: Monitor Progress

```bash
# Check job status
squeue -u $USER

# Watch output from first job
tail -f logs/enrich_*_0.out

# Check completion
ls riboseq_batch_output/*_enriched.json | wc -l
# Should reach 61 when done
```

---

## Interactive Approach (For Testing)

### Step 1: Start Interactive Session

```bash
cd /path/to/Metadata-Curation/omics-extractor
bash cluster/interactive_enrich.sh
```

This requests an interactive GPU node.

### Step 2: Setup Environment

```bash
# You're now on a GPU node
cd /path/to/Metadata-Curation/omics-extractor
source .venv-vllm/bin/activate

# Verify GPU
nvidia-smi
```

### Step 3: Test Model Loading

```bash
# Test model loads correctly
python cluster/test_local_models.py --model llama-3.3-70b
```

### Step 4: Enrich Single Project

```bash
# Pick a low-completeness project
PROJECT_ID="PRJDB10544"

# Run enrichment
python -m omics_extractor.cli enrich \
    riboseq_batch_output/${PROJECT_ID}_traceable.json \
    --output riboseq_batch_output/${PROJECT_ID}_enriched.json \
    --only-if-missing \
    --verbose
```

### Step 5: Review Results

```bash
# Compare before/after
python -c "
import json
before = json.load(open('riboseq_batch_output/${PROJECT_ID}_traceable.json'))
after = json.load(open('riboseq_batch_output/${PROJECT_ID}_enriched.json'))

print(f\"Completeness: {before['extraction_statistics']['completeness']['percentage']}% → {after['extraction_statistics']['completeness']['percentage']}%\")
"
```

---

## Resource Requirements

### Per Job

- **GPU:** 1x A100 (40GB or 80GB)
- **CPU:** 4 cores
- **Memory:** 32GB
- **Time:** 2 hours (typical: 5-30 minutes)
- **Storage:** ~5-10MB per output file

### Total (61 projects)

- **If sequential:** ~2-10 hours on 1 GPU
- **If parallel (10 GPUs):** ~20-60 minutes
- **Storage:** ~500MB for all outputs

### Model Requirements

For 70B models with VLLM:
- **VRAM:** ~40GB (with quantization)
- **Tensor Parallel:** 2 GPUs recommended for full precision
- **Quantization:** Can use 1 GPU with 8-bit quantization

---

## Configuration

### SLURM Script (`cluster/enrich_batch_slurm.sh`)

Key parameters:
```bash
#SBATCH --partition=gpu        # GPU partition
#SBATCH --gres=gpu:1           # 1 GPU per job
#SBATCH --cpus-per-task=4      # 4 CPUs
#SBATCH --mem=32G              # 32GB RAM
#SBATCH --time=02:00:00        # 2 hour limit
#SBATCH --array=0-60           # 61 jobs (0-60)
```

### Model Configuration (`config/model_paths.yaml`)

For single GPU (quantized):
```yaml
local_models:
  llama-3.3-70b:
    path: "/data/models/Llama-3.3-70B-Instruct"
    provider: "vllm"
    tensor_parallel_size: 1
    gpu_memory_utilization: 0.9
    quantization: "awq"  # or "gptq"
```

For dual GPU (full precision):
```yaml
local_models:
  llama-3.3-70b:
    path: "/data/models/Llama-3.3-70B-Instruct"
    provider: "vllm"
    tensor_parallel_size: 2
    gpu_memory_utilization: 0.9
```

---

## Monitoring and Debugging

### Check Job Status

```bash
# All jobs
squeue -u $USER

# Specific job
squeue -j <job_id>

# Completed jobs
sacct -u $USER --format=JobID,JobName,State,ExitCode,Elapsed
```

### View Logs

```bash
# Real-time output
tail -f logs/enrich_*.out

# Check errors
grep -i error logs/enrich_*.err

# Find failed jobs
grep "exit code" logs/enrich_*.err
```

### Resume Failed Jobs

The script automatically skips completed projects:

```bash
# Resubmit - will only process failed/skipped projects
sbatch --array=0-60 cluster/enrich_batch_slurm.sh
```

### Check Completion Rate

```bash
# Count enriched files
ENRICHED=$(ls riboseq_batch_output/*_enriched.json 2>/dev/null | wc -l)
echo "Enriched: $ENRICHED / 61"

# Show which projects are done
comm -23 \
  <(cat riboseq_batch_output/low_completeness_projects.txt | sort) \
  <(ls riboseq_batch_output/*_enriched.json | xargs -n1 basename | sed 's/_enriched.json//' | sort)
```

---

## Troubleshooting

### Issue: Model doesn't load

```bash
# Check VLLM installation
python -c "import vllm; print(vllm.__version__)"

# Check CUDA
nvidia-smi

# Test model path
ls -lh /data/models/Llama-3.3-70B-Instruct/
```

### Issue: Out of memory

**Solution 1:** Use quantization
```yaml
quantization: "awq"  # or "gptq"
```

**Solution 2:** Reduce batch size
```yaml
max_num_seqs: 16  # default: 256
```

**Solution 3:** Request more memory
```bash
#SBATCH --mem=64G
```

### Issue: Jobs taking too long

Check:
```bash
# See which project is slow
tail logs/enrich_*_*.out | grep "Processing project"

# Check GPU utilization
nvidia-smi
```

If underutilized:
- Increase batch size in model config
- Check if CPU-bound (increase CPUs)

### Issue: Job killed/timeout

Increase time limit:
```bash
sbatch --time=04:00:00 --array=0-60 cluster/enrich_batch_slurm.sh
```

---

## After Enrichment

### Aggregate Results

```bash
# Count improvements
python -c "
import json
from pathlib import Path

improvements = []
for enriched_file in Path('riboseq_batch_output').glob('*_enriched.json'):
    project_id = enriched_file.stem.replace('_enriched', '')
    original_file = enriched_file.parent / f'{project_id}_traceable.json'

    if not original_file.exists():
        continue

    orig = json.load(open(original_file))
    enrich = json.load(open(enriched_file))

    before = orig['extraction_statistics']['completeness']['percentage']
    after = enrich['extraction_statistics']['completeness']['percentage']

    if after > before:
        improvements.append({
            'project': project_id,
            'before': before,
            'after': after,
            'improvement': after - before
        })

print(f'Projects improved: {len(improvements)}')
avg_improvement = sum(p['improvement'] for p in improvements) / len(improvements) if improvements else 0
print(f'Average improvement: {avg_improvement:.1f}%')

# Show top improvements
for p in sorted(improvements, key=lambda x: x['improvement'], reverse=True)[:5]:
    print(f\"  {p['project']}: {p['before']:.1f}% → {p['after']:.1f}% (+{p['improvement']:.1f}%)\")
"
```

### Create Summary Report

```bash
# Generate enrichment report
python scripts/analyze_enrichment_results.py \
    riboseq_batch_output/ \
    --output enrichment_report.md
```

---

## Cost Estimation

### Compute Cost (Example)

**Assumptions:**
- A100 GPU: $2/hour
- 61 projects @ 15 min avg = 15 hours GPU time
- Parallel on 10 GPUs = 1.5 hours wall time

**Cost:**
- Sequential: 15 hours × $2 = $30
- Parallel (10 GPUs): 1.5 hours × $2 × 10 = $30

**vs Claude API:**
- 1,916 samples × ~$0.02/sample = ~$38
- Similar cost, but local = reusable, private, customizable

---

## Best Practices

1. **Test First:** Run interactive session on 1-2 projects before batch
2. **Monitor Early:** Watch first few jobs for errors
3. **Resume Capability:** Script automatically skips completed projects
4. **Save Logs:** Keep logs for debugging and benchmarking
5. **Resource Right-Sizing:** Adjust based on first few jobs
6. **Parallel Limit:** Don't overwhelm cluster - use `--array=0-60%10` to limit to 10 concurrent jobs

---

## Quick Reference

```bash
# Prepare batch
python cluster/prepare_enrichment_batch.py

# Submit batch (limit to 10 concurrent)
sbatch --array=0-60%10 cluster/enrich_batch_slurm.sh

# Monitor
squeue -u $USER
tail -f logs/enrich_*.out

# Check completion
ls riboseq_batch_output/*_enriched.json | wc -l

# Resume if needed
sbatch --array=0-60%10 cluster/enrich_batch_slurm.sh
```

---

## Next Steps

After enrichment completes:
1. Analyze improvement statistics
2. Compare baseline vs LLM extraction
3. Identify which models work best
4. Benchmark against Claude API
5. Generate final metadata export for RiboSeq.org
