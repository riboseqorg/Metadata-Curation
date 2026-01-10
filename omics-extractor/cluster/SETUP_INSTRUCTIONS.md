# Cluster Setup Instructions

## Your Paths

- **Code directory**: `/hps/software/users/ensembl/genebuild/jackt/RiboSeq/Metadata-Curation/omics-extractor`
- **Work directory**: `/hps/nobackup/flicek/ensembl/genebuild/jackt/metadata-curation`

## Setup Steps

### 1. Create Work Directory Structure

```bash
cd /hps/nobackup/flicek/ensembl/genebuild/jackt/metadata-curation

# Create necessary directories
mkdir -p riboseq_batch_output
mkdir -p logs
mkdir -p config
mkdir -p models  # Optional: if downloading models locally
```

### 2. Copy Batch Results

You need the traceable JSON files from your local batch extraction:

```bash
# On your local machine, from Metadata-Curation/omics-extractor:
rsync -avz riboseq_batch_output/*.json \
    your-cluster:/hps/nobackup/flicek/ensembl/genebuild/jackt/metadata-curation/riboseq_batch_output/
```

### 3. Setup VLLM Environment

```bash
cd /hps/software/users/ensembl/genebuild/jackt/RiboSeq/Metadata-Curation/omics-extractor

# Run setup script
bash cluster/setup_vllm.sh

# This creates .venv-vllm/ and installs:
# - vllm
# - torch
# - transformers
# - your omics-extractor package
```

### 4. Download or Locate Model

Option A - Use existing cluster model:
```bash
# Find if Llama-3.3-70B is already available
ls /hps/software/models/  # or wherever models are stored
```

Option B - Download model to your work directory:
```bash
cd /hps/nobackup/flicek/ensembl/genebuild/jackt/metadata-curation/models

# Using huggingface-cli (if available)
huggingface-cli download meta-llama/Llama-3.3-70B-Instruct \
    --local-dir Llama-3.3-70B-Instruct

# Or using git-lfs
git clone https://huggingface.co/meta-llama/Llama-3.3-70B-Instruct
```

### 5. Configure Model Paths

```bash
cd /hps/nobackup/flicek/ensembl/genebuild/jackt/metadata-curation

# Copy example config
cp /hps/software/users/ensembl/genebuild/jackt/RiboSeq/Metadata-Curation/omics-extractor/config/model_paths.example.yaml \
   config/model_paths.yaml

# Edit with your model path
nano config/model_paths.yaml
```

Set your model path:
```yaml
local_models:
  llama-3.3-70b:
    path: "/hps/nobackup/flicek/ensembl/genebuild/jackt/metadata-curation/models/Llama-3.3-70B-Instruct"
    provider: "vllm"
    tensor_parallel_size: 1  # Use 2 for dual-GPU
    gpu_memory_utilization: 0.9
```

### 6. Test Interactive Session

```bash
# Request GPU session
srun --partition=gpu \
     --gres=gpu:1 \
     --cpus-per-task=4 \
     --mem=32G \
     --time=02:00:00 \
     --pty bash

# Once on GPU node:
cd /hps/software/users/ensembl/genebuild/jackt/RiboSeq/Metadata-Curation/omics-extractor
source .venv-vllm/bin/activate

# Set config path
export OMICS_EXTRACTOR_CONFIG=/hps/nobackup/flicek/ensembl/genebuild/jackt/metadata-curation/config/model_paths.yaml

# Verify GPU
nvidia-smi

# Test model loading
python cluster/test_local_models.py --model llama-3.3-70b
```

### 7. Run Test Enrichment

```bash
# Still in interactive session
cd /hps/nobackup/flicek/ensembl/genebuild/jackt/metadata-curation

# Get first project to test
PROJECT_ID=$(head -1 riboseq_batch_output/low_completeness_projects.txt)
echo "Testing with: $PROJECT_ID"

# Run enrichment
python -m omics_extractor.cli enrich \
    riboseq_batch_output/${PROJECT_ID}_traceable.json \
    --output riboseq_batch_output/${PROJECT_ID}_enriched.json \
    --only-if-missing \
    --verbose
```

### 8. Check Results

```bash
python -c "
import json
project_id = '${PROJECT_ID}'
before = json.load(open(f'riboseq_batch_output/{project_id}_traceable.json'))
after = json.load(open(f'riboseq_batch_output/{project_id}_enriched.json'))

print(f'Project: {project_id}')
print(f'Completeness: {before[\"extraction_statistics\"][\"completeness\"][\"percentage\"]}% → {after[\"extraction_statistics\"][\"completeness\"][\"percentage\"]}%')
print(f'High-confidence: {before[\"extraction_statistics\"][\"field_statistics\"][\"high_confidence\"]} → {after[\"extraction_statistics\"][\"field_statistics\"][\"high_confidence\"]}')
"
```

## Next: Submit Batch Jobs

Once the test works, prepare the batch submission:

### 1. Edit SLURM Script

```bash
nano /hps/software/users/ensembl/genebuild/jackt/RiboSeq/Metadata-Curation/omics-extractor/cluster/enrich_batch_slurm.sh
```

Update these lines:
```bash
WORK_DIR=/hps/nobackup/flicek/ensembl/genebuild/jackt/metadata-curation
CODE_DIR=/hps/software/users/ensembl/genebuild/jackt/RiboSeq/Metadata-Curation/omics-extractor
```

### 2. Submit Jobs

```bash
cd /hps/software/users/ensembl/genebuild/jackt/RiboSeq/Metadata-Curation/omics-extractor

# Submit with concurrency limit
sbatch --array=0-60%10 cluster/enrich_batch_slurm.sh
```

### 3. Monitor

```bash
# Check queue
squeue -u $USER

# Watch first job output
tail -f /hps/nobackup/flicek/ensembl/genebuild/jackt/metadata-curation/logs/enrich_*_0.out

# Count completed
ls /hps/nobackup/flicek/ensembl/genebuild/jackt/metadata-curation/riboseq_batch_output/*_enriched.json | wc -l
```

## Troubleshooting

### Model doesn't load
- Check GPU with `nvidia-smi`
- Verify model path exists
- Check VRAM requirements (70B needs ~40GB)

### Out of memory
- Use quantization in model_paths.yaml: `quantization: "awq"`
- Request more memory: `--mem=64G`
- Use 2 GPUs: `tensor_parallel_size: 2`

### Jobs fail
- Check logs: `grep -i error logs/enrich_*.err`
- Verify config path is correct
- Test single project in interactive mode first
