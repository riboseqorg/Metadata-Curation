# Output Formats and Provenance Tracking

## Overview

The omics-extractor provides **multiple output formats** with **varying levels of detail** to help end users understand:

1. **What was extracted** - The metadata values themselves
2. **Where it came from** - Data source provenance (BioSample, SRA, LLM, etc.)
3. **How confident we are** - Two-factor confidence scoring
4. **What transformations were applied** - Ontology mapping, normalization, etc.

This allows users to choose the appropriate level of detail for their use case, from simple values-only CSV exports to complete extraction audit trails.

## The Three Detail Levels

### 1. Simple Format - Values Only

**Use case:** Quick data exploration, spreadsheet users who just need the metadata.

**What's included:**
- Just the metadata values (strings)
- No provenance information
- Compact and easy to read

**Example:**
```json
{
  "sample_id": "SRS12345",
  "organism": "Mus musculus",
  "tissue": "testis",
  "strain": "C57BL/6",
  "age": "8 weeks",
  "sex": "male"
}
```

**When to use:**
- Quick data exploration
- Generating simple reports
- Exporting to spreadsheets for non-technical users
- When provenance tracking is not needed

**CLI usage:**
```bash
omics-extract extract PRJNA1170270 --detail simple
```

---

### 2. Standard Format - Values + Provenance Summary (Recommended)

**Use case:** Most common usage - understand what was extracted and how reliable it is.

**What's included:**
- Field values
- Data source (biosample, sra, geo, llm, manual)
- Overall confidence score (0.0-1.0)
- Ontology term IDs (if applicable)

**Example:**
```json
{
  "sample_id": "SRS12345",
  "organism": {
    "value": "Mus musculus",
    "source": "biosample",
    "confidence": 0.95,
    "ontology_term": "NCBITaxon:10090"
  },
  "tissue": {
    "value": "testis",
    "source": "biosample",
    "confidence": 0.95,
    "ontology_term": "UBERON:0000473"
  },
  "strain": {
    "value": "C57BL/6",
    "source": "biosample",
    "confidence": 1.0,
    "ontology_term": null
  }
}
```

**When to use:**
- Production pipelines with quality control
- Downstream analysis requiring confidence scores
- When you need to know where data came from
- Most common use case - **recommended default**

**CLI usage:**
```bash
omics-extract extract PRJNA1170270 --detail standard  # Default
```

---

### 3. Detailed Format - Complete Provenance Chain

**Use case:** Debugging, auditing, understanding exactly how metadata was processed.

**What's included:**
- Everything from Standard format, PLUS:
- **Confidence breakdown:**
  - Total confidence (field × normalization)
  - Field confidence (how specific the source field was)
  - Normalization confidence (how much transformation was needed)
- **Ontology mapping details:**
  - Ontology term ID
  - Canonical ontology label
- **Extraction details:**
  - Extraction method (structured_field, llm, ner, regex)
  - Timestamp
  - Original extracted text (before normalization)
- **Notes/warnings** about extraction issues

**Example:**
```json
{
  "tissue": {
    "value": "testis",
    "source": "biosample",
    "source_id": "SAMN12345",
    "confidence": {
      "total": 0.950,
      "field": 1.0,
      "normalization": 0.95
    },
    "ontology": {
      "term": "UBERON:0000473",
      "label": "testis"
    },
    "extraction": {
      "method": "structured_field",
      "timestamp": "2025-01-04T10:30:00",
      "original_text": "testis"
    },
    "notes": null
  }
}
```

**When to use:**
- Debugging extraction issues
- Auditing data quality and provenance
- Understanding the full transformation chain
- Compliance requirements (full traceability)
- Research into extraction methods
- Quality assessment of LLM vs structured extraction

**CLI usage:**
```bash
omics-extract extract PRJNA1170270 --detail detailed
```

---

## Extraction Statistics

**All formats include extraction statistics** regardless of detail level. These provide quality metrics:

### 1. Field Coverage
How many samples have each field populated:
```json
"field_coverage": {
  "organism": {"count": 10, "percentage": 100.0},
  "tissue": {"count": 8, "percentage": 80.0},
  "cell_line": {"count": 2, "percentage": 20.0},
  "strain": {"count": 9, "percentage": 90.0}
}
```

### 2. Confidence Distribution
Distribution of confidence scores across all extracted fields:
```json
"confidence_distribution": {
  "high (≥0.8)": 45,
  "medium (0.5-0.8)": 12,
  "low (<0.5)": 3,
  "total_fields": 60
}
```

### 3. Source Distribution
Where the data came from:
```json
"source_distribution": {
  "biosample": 42,
  "sra": 8,
  "llm": 7,
  "geo": 3
}
```

### 4. Ontology Mapping Coverage
How many fields were successfully mapped to ontologies:
```json
"ontology_mapping": {
  "mapped": 28,
  "total": 30,
  "percentage": 93.3
}
```

### 5. Sample Completeness
How many samples have critical fields populated:
```json
"completeness": {
  "complete_samples": 9,
  "percentage": 90.0,
  "criteria": "organism + (tissue OR cell_line)"
}
```

---

## Additional Export Formats

### CSV Export

Export to CSV for spreadsheet tools:

**Simple CSV (values only):**
```bash
omics-extract extract PRJNA1170270 --csv
```

Produces:
```csv
sample_id,organism,tissue,cell_line,strain,age,sex
SRS001,Mus musculus,testis,,C57BL/6,8 weeks,male
SRS002,Mus musculus,liver,,C57BL/6,8 weeks,female
```

**CSV with Provenance:**
```bash
omics-extract extract PRJNA1170270 --csv --csv-provenance
```

Produces:
```csv
sample_id,organism,tissue,strain,organism_source,organism_confidence,organism_ontology,tissue_source,tissue_confidence,tissue_ontology,...
SRS001,Mus musculus,testis,C57BL/6,biosample,0.95,NCBITaxon:10090,biosample,0.95,UBERON:0000473,...
```

**Python API:**
```python
from omics_extractor.output.formatters import export_to_csv

export_to_csv(
    samples,
    "output.csv",
    include_provenance=True  # Add source/confidence columns
)
```

---

### Provenance Summary Report

Generate a human-readable markdown summary:

**CLI:**
```bash
omics-extract extract PRJNA1170270 --provenance-report
```

**Output:**
```markdown
# Metadata Extraction Provenance Summary

**Total Samples:** 10

## Field Coverage

| Field | Samples | Coverage |
|-------|---------|----------|
| organism | 10 | 100.0% |
| tissue | 8 | 80.0% |
| strain | 9 | 90.0% |

## Confidence Distribution

- **High confidence (≥0.8):** 45 / 60 (75.0%)
- **Medium confidence (0.5-0.8):** 12 / 60 (20.0%)
- **Low confidence (<0.5):** 3 / 60 (5.0%)

## Data Sources

- **biosample:** 42 fields
- **sra:** 8 fields
- **llm:** 7 fields

## Ontology Mapping

**28 / 30 fields mapped (93.3%)**

## Sample Completeness

**9 / 10 samples complete (90.0%)**
_Criteria: organism + (tissue OR cell_line)_
```

**Python API:**
```python
from omics_extractor.output.formatters import create_provenance_summary

summary_markdown = create_provenance_summary(samples)
with open("provenance_report.md", "w") as f:
    f.write(summary_markdown)
```

---

## Understanding Confidence Scores

### Two-Factor Confidence System

Each extracted field has a **two-factor confidence score**:

**1. Field Confidence** (0.0-1.0)
- How specific/reliable the source field is
- Examples:
  - `tissue_type` field → 1.0 (highly specific)
  - `organism` field → 1.0 (exact match)
  - Generic `description` field → 0.3 (might contain tissue info, but not specific)
  - LLM extraction from abstract → 0.5 (inferred)

**2. Normalization Confidence** (0.0-1.0)
- How much transformation was needed
- Examples:
  - Exact match to ontology: 1.0
  - Close synonym match: 0.95
  - Partial text match: 0.8
  - Inferred/derived value: 0.5

**Total Confidence = Field Confidence × Normalization Confidence**

Example:
```
Tissue extracted from specific "tissue_type" field: 1.0 × 0.95 = 0.95
Tissue inferred by LLM from abstract: 0.5 × 0.9 = 0.45
```

### Interpreting Confidence Scores

| Score | Interpretation | Action |
|-------|---------------|--------|
| ≥ 0.8 | **High confidence** | Trust and use directly |
| 0.5-0.8 | **Medium confidence** | Review for critical applications |
| < 0.5 | **Low confidence** | Manual review recommended |

---

## Python API Usage

### Creating Formatted Reports

```python
from omics_extractor.extraction.builder import build_project_metadata
from omics_extractor.output.formatters import MetadataFormatter

# Extract metadata
metadata = build_project_metadata("PRJNA1170270")

# Create formatted report
report = MetadataFormatter.create_extraction_report(
    study=metadata["study"],
    samples=metadata["samples"],
    runs=metadata["runs"],
    detail_level="standard",  # or "simple", "detailed"
)

# Save to JSON
import json
with open("output.json", "w") as f:
    json.dump(report, f, indent=2)
```

### Formatting Individual Fields

```python
from omics_extractor.output.formatters import MetadataFormatter

# Get a sample
sample = metadata["samples"]["SRS12345"]

# Format as simple (just values)
simple = MetadataFormatter.format_sample_simple(sample)
# Returns: {"sample_id": "SRS12345", "organism": "Mus musculus", ...}

# Format as standard (with provenance)
standard = MetadataFormatter.format_sample_standard(sample)
# Returns: {"sample_id": "SRS12345", "organism": {"value": "Mus musculus", "source": "biosample", ...}, ...}

# Format as detailed (complete provenance chain)
detailed = MetadataFormatter.format_sample_detailed(sample)
# Returns: Full nested structure with all provenance details
```

---

## Use Case Examples

### Use Case 1: Quick Data Exploration

**Goal:** Just see what metadata is available, don't care about provenance.

```bash
omics-extract extract PRJNA1170270 --detail simple --csv
```

Opens `PRJNA1170270.csv` in Excel, explores data.

---

### Use Case 2: Production Pipeline with QC

**Goal:** Extract metadata, filter by confidence, use high-quality data only.

```bash
omics-extract extract PRJNA1170270 --detail standard
```

```python
import json

with open("PRJNA1170270_metadata.json") as f:
    data = json.load(f)

# Filter high-confidence samples
high_quality_samples = {
    sample_id: sample
    for sample_id, sample in data["samples"].items()
    if sample["organism"]["confidence"] >= 0.8
    and (
        sample.get("tissue", {}).get("confidence", 0) >= 0.8 or
        sample.get("cell_line", {}).get("confidence", 0) >= 0.8
    )
}

print(f"High quality: {len(high_quality_samples)}/{len(data['samples'])} samples")
```

---

### Use Case 3: Debugging Extraction Issues

**Goal:** Understand why tissue wasn't extracted for certain samples.

```bash
omics-extract extract PRJNA1170270 --detail detailed
```

```python
import json

with open("PRJNA1170270_metadata.json") as f:
    data = json.load(f)

# Find samples missing tissue
for sample_id, sample in data["samples"].items():
    if "tissue" not in sample.get("biological_context", {}):
        print(f"\n{sample_id}: No tissue extracted")
        print(f"  Organism: {sample['organism']['value']}")
        print(f"  Extraction method: {sample['organism']['extraction']['method']}")
        print(f"  Original text: {sample['organism']['extraction']['original_text']}")
```

---

### Use Case 4: Publication-Ready Report

**Goal:** Generate summary statistics for a methods section.

```bash
omics-extract extract PRJNA1170270 --provenance-report
```

Produces `PRJNA1170270_provenance.md` with:
- Field coverage statistics
- Confidence distribution
- Data source breakdown
- Ontology mapping success rate

Include this in supplementary materials.

---

## CLI Command Reference

### Extract Command

```bash
omics-extract extract <BIOPROJECT> [OPTIONS]
```

**Options:**
- `--output PATH` - Output JSON file path
- `--detail {simple,standard,detailed}` - Detail level (default: standard)
- `--csv` - Also export to CSV
- `--csv-provenance` - Include source/confidence columns in CSV
- `--provenance-report` - Generate markdown provenance summary
- `--verbose` - Show detailed progress

**Examples:**

```bash
# Standard extraction (recommended)
omics-extract extract PRJNA1170270

# Simple format with CSV export
omics-extract extract PRJNA1170270 --detail simple --csv

# Detailed format for debugging
omics-extract extract PRJNA1170270 --detail detailed

# Full package: detailed JSON + CSV with provenance + report
omics-extract extract PRJNA1170270 \
  --detail detailed \
  --csv \
  --csv-provenance \
  --provenance-report
```

---

## Format Comparison Table

| Feature | Simple | Standard | Detailed |
|---------|--------|----------|----------|
| **Field values** | ✓ | ✓ | ✓ |
| **Data source** | ✗ | ✓ | ✓ |
| **Overall confidence** | ✗ | ✓ | ✓ |
| **Confidence breakdown** | ✗ | ✗ | ✓ |
| **Ontology term ID** | ✗ | ✓ | ✓ |
| **Ontology label** | ✗ | ✗ | ✓ |
| **Extraction method** | ✗ | ✗ | ✓ |
| **Timestamp** | ✗ | ✗ | ✓ |
| **Original text** | ✗ | ✗ | ✓ |
| **Notes/warnings** | ✗ | ✗ | ✓ |
| **Extraction statistics** | ✓ | ✓ | ✓ |
| **File size** | Small | Medium | Large |
| **Readability** | High | Medium | Low |
| **Best for** | Exploration | Production | Debugging |

---

## Best Practices

### 1. Choose the Right Format

- **Starting out?** Use `--detail simple --csv` to explore
- **Production pipeline?** Use `--detail standard` (default)
- **Debugging issues?** Use `--detail detailed`
- **Publication?** Use `--provenance-report`

### 2. Quality Control

```python
# Filter by confidence in production
high_confidence = {
    sid: s for sid, s in samples.items()
    if s["organism"]["confidence"] >= 0.8
}

# Flag samples needing review
review_needed = {
    sid: s for sid, s in samples.items()
    if s.get("tissue", {}).get("confidence", 0) < 0.5
}
```

### 3. Provenance Tracking

Always preserve the full provenance chain for:
- Published datasets
- Clinical/regulatory applications
- Long-term archives

Use `--detail detailed` for archival copies.

### 4. Export Multiple Formats

For maximum flexibility:
```bash
omics-extract extract PRJNA1170270 \
  --detail standard \
  --csv \
  --csv-provenance \
  --provenance-report
```

Produces:
- `PRJNA1170270_metadata.json` - Main output with provenance
- `PRJNA1170270_metadata.csv` - Spreadsheet-friendly version
- `PRJNA1170270_metadata.md` - Human-readable summary

---

## Next Steps

- See `scripts/demo_output_formats.py` for working examples
- See `PRACTICAL_LLM_STRATEGY.md` for extraction strategy
- See `CLUSTER_SETUP_COMPLETE.md` for local model usage
