# LLM Integration for Metadata Enrichment

This document describes the LLM-based metadata enrichment system in omics-extractor.

## Overview

The LLM enrichment system uses a **two-phase architecture** designed to separate I/O-bound operations (data fetching) from compute-intensive operations (LLM inference):

- **Phase 1 (CPU nodes)**: Extract structured metadata from databases
- **Phase 2 (GPU nodes)**: Enrich missing fields using LLM inference

This separation ensures GPU nodes aren't wasted on network I/O and database queries.

## Architecture

### Phase 1: Data Extraction (CPU)

```bash
# Run on cheap CPU nodes
omics-extract extract PRJNA1170270 --output data.json
```

This phase:
- Fetches from BioProject, BioSample, SRA, GEO, PubMed APIs
- Applies ontology normalization (UBERON, CL, NCBITaxon)
- Calculates two-factor confidence scores
- Outputs JSON with structured metadata

**No GPU required** - pure I/O and data wrangling.

### Phase 2: LLM Enrichment (GPU)

```bash
# Run on GPU nodes (A100, etc.)
omics-extract enrich data.json --only-if-missing
```

This phase:
- Reads metadata JSON from Phase 1
- Identifies samples with missing critical fields
- Uses LLM to extract tissue, cell type, cell line, etc.
- Only fills gaps - never overwrites structured data

**GPU-optimized** for inference efficiency.

## Two-Factor Confidence Scoring

Every metadata field has two confidence components:

### 1. Field Confidence (WHERE did it come from?)

- **1.0**: Exact field match (e.g., "tissue" ’ tissue)
- **0.95**: Type variant (e.g., "tissue_type" ’ tissue)
- **0.7**: Generic source (e.g., "source_name" ’ tissue)
- **0.5**: Very generic (e.g., "isolation_source" ’ tissue)
- **0.5**: LLM inferred from study description

### 2. Normalization Confidence (HOW well does it match ontology?)

- **1.0**: Exact ontology match (e.g., "brain" ’ UBERON:0000955)
- **0.95**: Case-insensitive match (e.g., "Brain" ’ "brain")
- **0.9**: Synonym match (e.g., "breast" ’ "mammary gland")
- **0.7**: Fuzzy match
- **0.5**: No ontology match (keeps raw value)

### Final Confidence

```
final_confidence = field_confidence × normalization_confidence
```

**Example 1**: Exact tissue from exact field
```
"brain" from "tissue" field
’ Field: 1.0, Normalization: 1.0
’ Final: 1.0 (highest confidence)
’ UBERON:0000955
```

**Example 2**: LLM-inferred tissue with good normalization
```
"liver" inferred from study description via LLM
’ Field: 0.5 (LLM), Normalization: 1.0 (exact match to UBERON:0002107)
’ Final: 0.5 (medium confidence)
```

**Example 3**: Generic field, unmapped value
```
"weird_tissue_type" from "isolation_source"
’ Field: 0.5, Normalization: 0.5 (no match)
’ Final: 0.25 (low confidence - needs review)
```

## CLI Commands

### Single File Enrichment

```bash
# Enrich all samples
omics-extract enrich metadata.json

# Enrich only samples with missing critical fields
omics-extract enrich metadata.json --only-if-missing

# Specify output location
omics-extract enrich metadata.json --output enriched.json

# Provide API key directly
omics-extract enrich metadata.json --api-key sk-ant-...
```

### Batch Enrichment (GPU Cluster)

```bash
# Batch enrich multiple projects
omics-extract batch-enrich *_metadata.json --output-dir enriched/

# Increase concurrency for better GPU utilization
omics-extract batch-enrich projects/*.json \
  --output-dir enriched/ \
  --workers 8 \
  --rate-limit 100

# Resume interrupted batch job (uses checkpoints)
omics-extract batch-enrich projects/*.json \
  --output-dir enriched/ \
  --resume
```

### Batch Processing Features

1. **Concurrent API Requests**: Use `--workers` to parallelize (default: 4)
2. **Rate Limiting**: Automatic rate limiting to avoid API throttling
3. **Checkpointing**: Progress saved every 10 samples
4. **Resumption**: Use `--resume` to continue from checkpoint after interruption
5. **Statistics**: Detailed reporting of fields added per file

## LLM Provider Configuration

### Claude API (Default)

Set your API key:
```bash
export ANTHROPIC_API_KEY=sk-ant-...
```

Or pass directly:
```bash
omics-extract enrich data.json --api-key sk-ant-...
```

Uses **Claude 3.5 Sonnet** with:
- Temperature: 0 (deterministic)
- Max tokens: 1000
- Structured JSON output

### Local Models (Future)

For running on local GPUs without API calls, the system can be extended to support:
- **Llama 3** via transformers
- **Mixtral** via VLLM
- **Custom fine-tuned models**

To add local model support, extend `llm_extractor.py`:

```python
def extract_with_local_model(
    study_title: str,
    study_description: str,
    model_path: str = "meta-llama/Llama-3-70B",
) -> LLMExtractionResult:
    """Extract using local model."""
    # Load model with VLLM or transformers
    # Generate with same prompt template
    # Parse JSON response
    pass
```

## Performance Considerations

### macOS Development

On Mac (M1/M2/M3):
- Use Claude API for development/testing
- Small batches work fine on CPU
- Good for prototyping and validation

### Production on A100

On GPU clusters:
- **Batch enrichment** is critical for efficiency
- Use `--workers 8-16` for maximum throughput
- Set `--rate-limit` based on your API tier
- Monitor GPU utilization (should see high usage with batch mode)

### Cost Optimization

1. **Only enrich what's needed**: Use `--only-if-missing` to skip samples with complete metadata
2. **Batch multiple projects**: Use `batch-enrich` to process many projects in one GPU job
3. **Resume from checkpoints**: Don't re-process samples if job crashes
4. **Filter before enrichment**: Pre-filter samples that don't need LLM (e.g., GEO samples often have good metadata)

## Workflow Examples

### Simple Workflow

```bash
# Phase 1: Extract (CPU)
omics-extract extract PRJNA1170270

# Phase 2: Enrich (GPU)
omics-extract enrich PRJNA1170270_metadata.json --only-if-missing
```

### Batch Processing Workflow

```bash
# Phase 1: Extract multiple projects (CPU cluster)
for project in PRJNA1170270 PRJNA1176138 PRJNA1180945; do
  omics-extract extract $project --output projects/${project}.json
done

# Phase 2: Batch enrich (single GPU node)
omics-extract batch-enrich projects/*.json \
  --output-dir enriched/ \
  --workers 8 \
  --resume
```

### Large-Scale Pipeline

```bash
#!/bin/bash
# Extract 1000 projects on CPU cluster (parallel)
parallel -j 20 omics-extract extract {} --output data/{}.json ::: $(cat projects.txt)

# Batch enrich on GPU node (A100)
omics-extract batch-enrich data/*.json \
  --output-dir enriched/ \
  --workers 16 \
  --rate-limit 200 \
  --resume
```

## Output Format

Enriched metadata includes:

```json
{
  "bioproject_id": "PRJNA1170270",
  "samples": {
    "SAMN12345": {
      "tissue": {
        "value": "liver",
        "source": "llm",
        "confidence": 0.5,
        "field_confidence": 0.5,
        "normalization_confidence": 1.0,
        "ontology_term": "UBERON:0002107",
        "ontology_label": "liver",
        "extracted_text": "Inferred from study context",
        "extraction_method": "llm"
      }
    }
  },
  "enrichment_metadata": {
    "timestamp": "2025-01-03T12:00:00",
    "samples_enriched": 42,
    "fields_added": {
      "tissue": 15,
      "cell_type": 8,
      "cell_line": 3
    }
  }
}
```

## Monitoring and Debugging

### Verbose Mode

```bash
omics-extract enrich data.json --verbose
```

Shows:
- Per-sample enrichment progress
- API call latency
- Fields added per sample
- Errors and warnings

### Checkpoint Files

During batch enrichment, checkpoint files track progress:

```bash
enriched/PRJNA1170270_checkpoint.json
```

Contains:
- List of completed sample IDs
- Current statistics
- Timestamp

### Validation

After enrichment, validate metadata:

```bash
omics-extract validate enriched_metadata.json
```

Checks:
- Missing critical fields
- Confidence scores
- Ontology term validity

## Troubleshooting

### API Key Issues

```
Error: No API key provided
```

**Solution**: Set `ANTHROPIC_API_KEY` environment variable or use `--api-key`

### Rate Limiting

```
Error: Rate limit exceeded
```

**Solution**: Reduce `--workers` or increase `--rate-limit` delay

### Out of Memory (GPU)

```
CUDA out of memory
```

**Solution**: This shouldn't happen with Claude API. If using local models, reduce batch size.

### Checkpoint Corruption

```
Error: Invalid checkpoint file
```

**Solution**: Delete checkpoint file and use `--resume` to regenerate

## Future Enhancements

1. **Fine-tuned Models**: Train domain-specific models on curated metadata
2. **Active Learning**: Flag low-confidence extractions for human review
3. **Multi-modal**: Extract from figures and supplementary files
4. **Validation Loop**: LLM validates its own extractions
5. **Few-shot Examples**: Provide project-specific examples in prompt
