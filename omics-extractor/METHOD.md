# LLM-Enhanced Metadata Enrichment for Genomics Data

## Overview

This method combines structured database extraction with large language model (LLM) inference to comprehensively annotate biological samples from public sequencing repositories. The pipeline extracts metadata from NCBI databases (BioProject, BioSample, SRA), supplements with GEO and PubMed when available, then uses LLMs to infer missing biological context from study descriptions.

## Method

### 1. Structured Metadata Extraction

The pipeline first retrieves metadata from public databases using a waterfall strategy:

**Study-level metadata** is extracted from:
- **BioProject** (primary): Project title, description, organism, and publication links (100% coverage)
- **GEO Series** (supplemental): Enhanced project descriptions when available
- **PubMed**: Publication abstracts, authors, journal information, and DOI

**Sample-level metadata** is extracted from:
- **BioSample** (primary): Structured attributes including tissue, cell type, strain, genotype, treatment, age, and sex
- **SRA Run records**: Technical sequencing details (library strategy, platform, read counts) and links to parent BioSample records
- **GEO Samples** (supplemental): Legacy samples without BioSample records

Each extracted field is wrapped in a provenance structure tracking its source database, extraction method, and confidence score (0.0-1.0). Confidence scores reflect field specificity: exact structured fields receive 1.0, while inferred values receive lower scores.

### 2. LLM-Based Metadata Enrichment

After structured extraction, samples with incomplete metadata undergo LLM enrichment. The system identifies samples missing critical biological context (tissue, cell type, cell line, treatment, etc.) and constructs enrichment prompts from available study-level descriptions.

**Enrichment input** for each sample includes:
- BioProject title and description
- Sample title and description (when available)
- Publication abstract (when available)
- Already-extracted metadata to avoid duplication

**LLM extraction** uses a structured JSON prompt requesting 17 biological metadata fields:

*Core biological metadata:* organism, tissue, cell_type, cell_line, strain, developmental_stage, genotype, age, sex

*Experimental conditions:* condition, treatment, timepoint, replicate, batch

*Disease/perturbation:* disease, stress, temperature, growth_condition

The LLM is instructed to extract only information it can confidently infer from the provided context, returning null for unavailable fields. Each non-null field must include a confidence score (0.0-1.0) reflecting extraction certainty: 1.0 for explicitly stated facts, 0.8 for clear implications, 0.6 for reasonable inferences, and 0.4 for educated guesses.

**Final confidence calculation** combines field-level confidence (how specific the source text is, typically 0.5 for study-level descriptions) with normalization confidence (LLM's stated confidence) to produce conservative estimates: `final_confidence = field_confidence × llm_confidence`.

### 3. Ontology Normalization

Extracted tissue and cell type values undergo automatic ontology mapping to standardize terminology:

**Tissue normalization** maps to UBERON (Uber-anatomy Ontology) using a curated dictionary of ~40 common tissues. Case-insensitive exact matches receive confidence 1.0, synonyms receive 0.9 (e.g., "fat" → "adipose tissue", UBERON:0001013), and fuzzy matches above 85% similarity receive 0.7.

**Cell type normalization** maps to Cell Ontology (CL) using ~25 curated cell types with similar confidence scoring. Common variants are handled (e.g., "red blood cell" → "erythrocyte", CL:0000232).

Unmapped terms retain their original LLM-extracted values but receive no ontology identifier, flagging them for potential manual curation. All normalized values replace raw LLM outputs in the enriched metadata, ensuring consistent terminology across datasets.

### 4. Provenance Tracking

Every metadata field is stored with complete provenance:

```json
{
  "value": "brain",
  "source": "llm_enrichment",
  "source_id": "llm_enrichment_SRS12345",
  "confidence": 0.40,
  "field_confidence": 0.50,
  "normalization_confidence": 0.80,
  "extraction_method": "llm",
  "ontology_term": "UBERON:0000955",
  "notes": "Inferred from study description mentioning neural tissue analysis"
}
```

This enables downstream users to filter by confidence thresholds, trace metadata origins, and identify fields requiring manual review.

## Implementation

The pipeline is implemented as a command-line tool (`omics-extract`) with support for local or API-based LLM inference:

**Batch enrichment** processes multiple projects in parallel:
```bash
omics-extract batch-enrich extracted/*.json \
    --provider vllm \
    --model-name mixtralai/Mixtral-8x22B-Instruct-v0.1 \
    --base-url http://localhost:8000/v1
```

**Local LLM deployment** uses VLLM with automatic GPU detection and tensor parallelism for cost-effective processing. The system supports any OpenAI-compatible API, including Claude API, local VLLM, or Transformers.

**Quality assessment tools** analyze enrichment results:
- `assess_enrichment.py`: Statistical analysis of field coverage and confidence distributions
- `review_enrichments.py`: Manual spot-checking with source text highlighting
- `json_to_table.py`: Export to tabular format (CSV/TSV) with columns for value, enrichment status, and ontology ID

## Performance

In a test batch of 20 BioProject studies (259 samples), the system enriched 347 metadata fields in 40 minutes using Mixtral-8x22B (~9 seconds per sample). The multi-source extraction strategy achieved 100% metadata coverage: all samples received organism and project-level annotations from BioProject/BioSample, while LLM enrichment filled tissue, cell type, and treatment fields that were absent from structured databases.

## Advantages

**Comprehensive coverage**: Combines authoritative structured metadata with LLM-inferred context, maximizing field completeness without sacrificing accuracy.

**Provenance transparency**: Every field tracks its origin, extraction method, and confidence, enabling quality filtering and identifying fields needing manual review.

**Ontology standardization**: Automatic mapping to UBERON and Cell Ontology ensures consistent terminology while preserving original values for unmapped terms.

**Scalable inference**: Support for local LLMs enables cost-effective processing of thousands of projects without per-token API costs.

**Extensible architecture**: Generic provider interface supports any LLM backend, and the schema accommodates additional metadata fields without code changes.

## Limitations

**LLM hallucination risk**: Inferred metadata may contain errors when study descriptions are ambiguous or incomplete. Conservative confidence scoring and manual review tools mitigate this risk.

**Limited ontology coverage**: The curated ontology mappings cover common tissues and cell types but require expansion for specialized domains. Terms without mappings are flagged for manual curation.

**Context limitations**: LLM enrichment relies on study-level descriptions, which may not capture sample-specific variation within a project. Future work could incorporate publication full-text for finer-grained extraction.

**Inference cost**: While local LLMs reduce costs, batch processing still requires GPU resources. The system optimizes throughput with parallel processing and efficient prompt construction.
