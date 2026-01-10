# Traceable Output Format - Summary

## What We Built

A **single hierarchical output format** that contains ALL levels of detail in one file:

```
Quick View (top level)
   ↓
Standard Provenance (one level down)
   ↓
Detailed Extraction Chain (deep level)
   ↓
Raw BioSample Data (bottom level)
```

## Key Features

### 1. All Levels in One File

No need to choose between "simple", "standard", or "detailed" output formats. **Every file contains all levels**, and users access the level they need:

```python
# Level 1: Quick access
tissue = sample["quick_view"]["tissue"]  # Just the value: "testis"

# Level 2: Standard provenance
tissue = sample["biological_metadata"]["tissue"]["provenance"]
# Returns: {source: "biosample", confidence: 1.0, ontology_term: "UBERON:0000473"}

# Level 3: Detailed debugging
details = sample["biological_metadata"]["tissue"]["extraction_details"]
# Returns: {field_confidence: 1.0, normalization_confidence: 1.0, method: "structured_field", ...}

# Level 4: Raw data
raw = sample["raw_biosample_attributes"]["tissue"]  # Original: "tissue: testis"
```

### 2. Complete Traceability

Every field shows the full chain from raw data to final value:

```
Raw BioSample: "tissue: testis"
     ↓ (field_confidence: 1.0 - specific field)
Extracted: "testis"
     ↓ (normalization_confidence: 1.0 - exact match)
Ontology: UBERON:0000473
     ↓ (total_confidence: 1.0)
Final: "testis"
```

### 3. Two-Factor Confidence

**Total Confidence = Field Confidence × Normalization Confidence**

- **Field Confidence** - How reliable was the source?
  - 1.0 = Specific field like "tissue_type"
  - 0.5 = LLM inferred from abstract
  - 0.3 = Generic description field

- **Normalization Confidence** - How much transformation?
  - 1.0 = Exact match to ontology
  - 0.9 = Close synonym
  - 0.5 = Derived/transformed

### 4. Comprehensive Statistics

Every output includes quality metrics:

- **Field Coverage** - How many samples have each field
- **Confidence Distribution** - Overall and per-field
- **Source Distribution** - Where data came from (biosample, sra, llm, etc.)
- **Ontology Mapping** - Success rate and per-field breakdown
- **Completeness** - How many samples meet quality criteria
- **Extraction Methods** - Breakdown by method (structured_field, llm, etc.)

### 5. Navigation Guide

Every file includes a guide:

```json
{
  "navigation_guide": {
    "quick_access": "samples[<id>].quick_view.<field>",
    "standard_provenance": "samples[<id>].biological_metadata.<field>.provenance",
    "detailed_debugging": "samples[<id>].biological_metadata.<field>.extraction_details",
    "raw_data": "samples[<id>].raw_biosample_attributes"
  }
}
```

## File Structure

```json
{
  "metadata_version": "2.0.0",
  "format": "traceable",
  "extraction_timestamp": "...",
  "description": "Hierarchical metadata with complete traceability...",

  "study": {
    "quick_view": {...},      // Just values
    "identifiers": {...},
    "core_metadata": {...},   // Full traceable format
    "publication": {...},
    "extraction_metadata": {...}
  },

  "samples": {
    "<sample_id>": {
      "quick_view": {...},               // Level 1
      "identifiers": {...},
      "biological_metadata": {           // Levels 2 & 3
        "<field>": {
          "value": "...",
          "provenance": {...},           // Level 2
          "extraction_details": {...}    // Level 3
        }
      },
      "raw_biosample_attributes": {...}  // Level 4
    }
  },

  "runs": [...],
  "extraction_statistics": {...},
  "navigation_guide": {...}
}
```

## Usage Examples

### Quick Exploration

```python
import json

with open("PRJNA1170270_traceable.json") as f:
    data = json.load(f)

# Just show me the data
for sample_id, sample in data["samples"].items():
    print(f"{sample_id}: {sample['quick_view']}")
```

### Production Pipeline with QC

```python
# Filter by confidence for production use
high_quality = {
    sid: s for sid, s in data["samples"].items()
    if s.get("biological_metadata", {}).get("tissue", {})
      .get("provenance", {}).get("confidence", 0) >= 0.8
}
```

### Debugging Extraction Issues

```python
# Why wasn't tissue extracted for this sample?
sample = data["samples"]["SRS123"]
if "tissue" in sample.get("biological_metadata", {}):
    tissue = sample["biological_metadata"]["tissue"]
    print(f"Value: {tissue['value']}")
    print(f"Confidence: {tissue['provenance']['confidence']}")
    print(f"Method: {tissue['extraction_details']['method']}")
    print(f"Original: {tissue['extraction_details']['original_text']}")
else:
    print("No tissue extracted")
    print(f"Raw attributes: {sample['raw_biosample_attributes']}")
```

### Publication-Ready Statistics

```python
# Get extraction quality metrics for methods section
stats = data["extraction_statistics"]
print(f"Samples: {stats['total_samples']}")
print(f"Completeness: {stats['completeness']['percentage']}%")
print(f"High confidence fields: {stats['confidence_distribution']['overall']['high (≥0.8)']}")
print(f"Ontology mapping: {stats['ontology_mapping']['overall']['percentage']}%")
```

## Testing on 100 RiboSeq Studies

The system has been tested and works successfully on multiple RiboSeq studies. To run on 100 studies:

```bash
cd /Users/jackt/projects/Metadata-Curation/omics-extractor
~/.local/bin/uv run python scripts/batch_extract_riboseq_traceable.py
```

This will:
1. Process 100 RiboSeq BioProjects
2. Output traceable format for each
3. Compute aggregate statistics
4. Show example output structure

**Test run on 5 projects (completed successfully):**
- PRJNA1170270: 1 sample, 4 runs ✓
- PRJDB10544: 2 samples, 2 runs ✓
- PRJDB10799: 16 samples, 16 runs ✓
- PRJEB12126: 40 samples, 40 runs ✓
- PRJNA1002596: 12 samples, 12 runs ✓

## CLI Usage

The CLI now outputs traceable format by default:

```bash
# Extract with full traceability
omics-extract extract PRJNA1170270

# Also export to CSV
omics-extract extract PRJNA1170270 --csv

# Also create provenance summary report
omics-extract extract PRJNA1170270 --provenance-report
```

Output files:
- `PRJNA1170270_metadata.json` - Full traceable format
- `PRJNA1170270_metadata.csv` - Quick view as CSV (optional)
- `PRJNA1170270_metadata.md` - Human-readable summary (optional)

## Benefits

### For End Users

**Exploration:** Use `quick_view` for fast data browsing

**Production:** Use `provenance` for quality control and filtering

**Debugging:** Use `extraction_details` to understand extraction issues

**Compliance:** Use full traceable format for regulatory/audit requirements

### For Developers

**Single Format:** No need to maintain multiple output formats

**Backward Compatible:** Old code accessing `.value` still works

**Forward Compatible:** New provenance fields can be added without breaking existing code

**Self-Documenting:** Every file includes navigation guide

### For Research

**Reproducibility:** Complete provenance chain from raw to final

**Quality Assessment:** Detailed confidence scores and statistics

**Method Comparison:** Can see which extraction methods work best

**Publication:** Statistics suitable for methods sections

## Documentation

- **Complete Guide:** `docs/TRACEABLE_OUTPUT_FORMAT.md`
- **Implementation:** `src/omics_extractor/output/traceable_format.py`
- **Examples:** `scripts/test_traceable_output.py`
- **Batch Processing:** `scripts/batch_extract_riboseq_traceable.py`

## Next Steps

The traceable format is now the default output. For the 100 RiboSeq study batch extraction:

```bash
# Run full batch (may take a while)
~/.local/bin/uv run python scripts/batch_extract_riboseq_traceable.py

# Or run in smaller batches
~/.local/bin/uv run python scripts/batch_extract_riboseq_traceable.py --verbose
```

This will create:
- 100 individual traceable JSON files
- Aggregate statistics across all projects
- Quality assessment reports
- Example outputs showing the format
