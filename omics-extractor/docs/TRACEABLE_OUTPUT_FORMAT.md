# Traceable Output Format

## Overview

The traceable output format provides **complete provenance tracking** from raw data to final values in a **single hierarchical structure**. Users can access different levels of detail within the same file:

1. **Quick View** - Just values (top level)
2. **Standard Provenance** - Source + confidence (one level down)
3. **Detailed Extraction Chain** - Complete traceability (deep level)
4. **Raw Data** - Original BioSample attributes (bottom level)

All levels are present in every output file, allowing users to drill down from high-level summaries to low-level debugging details.

## File Structure

```json
{
  "metadata_version": "2.0.0",
  "format": "traceable",
  "extraction_timestamp": "2026-01-04T...",
  "description": "Hierarchical metadata with complete traceability...",

  "study": {...},           // Study-level metadata
  "samples": {...},         // Sample-level metadata
  "runs": [...],            // Run-level metadata
  "extraction_statistics": {...},  // Quality metrics
  "navigation_guide": {...}  // How to access different levels
}
```

## Access Patterns

### Level 1: Quick View (Just Values)

**Use case:** Quick exploration, spreadsheet-like view

```json
{
  "samples": {
    "SRS22569018": {
      "quick_view": {
        "sample_id": "SRS22569018",
        "organism": "Mus musculus",
        "tissue": "testis",
        "strain": "C57BL/6",
        "sex": "male",
        "age": "P14"
      }
    }
  }
}
```

**Access:** `samples[sample_id].quick_view.<field>`

**Contains:** Just the metadata values as strings

---

### Level 2: Standard Provenance

**Use case:** Understanding data quality and sources

```json
{
  "samples": {
    "SRS22569018": {
      "biological_metadata": {
        "tissue": {
          "value": "testis",
          "provenance": {
            "source": "biosample",
            "source_id": "SAMN42448935",
            "confidence": 1.0,
            "ontology_term": "UBERON:0000473"
          }
        }
      }
    }
  }
}
```

**Access:** `samples[sample_id].biological_metadata.<field>.provenance`

**Contains:**
- `value` - The extracted value
- `source` - Where it came from (biosample, sra, llm, etc.)
- `source_id` - Which record (SAMN..., SRR..., etc.)
- `confidence` - Overall confidence score (0.0-1.0)
- `ontology_term` - Ontology ID if mapped (UBERON:..., CL:..., etc.)

---

### Level 3: Detailed Extraction Chain

**Use case:** Debugging, auditing, full traceability

```json
{
  "samples": {
    "SRS22569018": {
      "biological_metadata": {
        "tissue": {
          "value": "testis",
          "provenance": {
            "source": "biosample",
            "source_id": "SAMN42448935",
            "confidence": 1.0,
            "ontology_term": "UBERON:0000473"
          },
          "extraction_details": {
            "field_confidence": 1.0,
            "field_confidence_explanation": "Extracted from specific, high-quality field",
            "normalization_confidence": 1.0,
            "normalization_explanation": "Exact or near-exact match to standard term",
            "method": "structured_field",
            "original_text": "tissue: testis",
            "transformation": "tissue: testis → testis",
            "ontology_label": "testis"
          }
        }
      }
    }
  }
}
```

**Access:** `samples[sample_id].biological_metadata.<field>.extraction_details`

**Contains:**
- `field_confidence` - How specific/reliable the source field was
- `field_confidence_explanation` - Human-readable explanation
- `normalization_confidence` - How much transformation was needed
- `normalization_explanation` - Human-readable explanation
- `method` - Extraction method (structured_field, llm, ner, etc.)
- `original_text` - Text before any normalization
- `transformation` - Shows the transformation applied
- `ontology_label` - Canonical ontology term label

**Total confidence = field_confidence × normalization_confidence**

---

### Level 4: Raw BioSample Attributes

**Use case:** Complete audit trail, regulatory compliance

```json
{
  "samples": {
    "SRS22569018": {
      "raw_biosample_attributes": {
        "organism": "Mus musculus",
        "tissue": "testis",
        "strain": "C57BL/6",
        "sex": "male",
        "age": "P14",
        "developmental_stage": "juvenile"
      }
    }
  }
}
```

**Access:** `samples[sample_id].raw_biosample_attributes`

**Contains:** Original key-value pairs from BioSample before any processing

---

## Sample Structure

Each sample has this hierarchical structure:

```json
{
  "quick_view": {...},          // Level 1: Just values
  "identifiers": {...},         // IDs (sample, bioproject, biosample, gsm)
  "biological_metadata": {      // Level 2 & 3: Full traceability
    "organism": {...},
    "tissue": {...},
    "cell_line": {...},
    "cell_type": {...},
    "developmental_stage": {...},
    "strain": {...},
    "genotype": {...},
    "sex": {...},
    "age": {...}
  },
  "experimental_metadata": {    // If present
    "condition": {...},
    "treatment": {...},
    "timepoint": {...},
    "replicate": {...},
    "batch": {...}
  },
  "disease_perturbation": {     // If present
    "disease": {...},
    "stress": {...},
    "temperature": {...},
    "growth_condition": {...}
  },
  "descriptions": {             // Sample title/description
    "title": {...},
    "description": {...}
  },
  "raw_biosample_attributes": {...}  // Level 4: Original data
}
```

---

## Study Structure

```json
{
  "quick_view": {
    "bioproject_id": "PRJNA1170270",
    "title": "Ribo-seq of mouse spermatocytes",
    "organism": "Mus musculus"
  },
  "identifiers": {
    "bioproject_id": "PRJNA1170270",
    "gse_id": null,
    "sra_study_id": null,
    "pmid": null,
    "pmc_id": null
  },
  "core_metadata": {
    "title": {...},        // Full traceable format
    "description": {...},  // Full traceable format
    "organism": {...}      // Full traceable format
  },
  "publication": {         // If PMID exists
    "pmid": "...",
    "doi": {...},
    "title": {...},
    "journal": {...},
    "authors": [...],
    "abstract": {...}
  },
  "extraction_metadata": {
    "extraction_date": "...",
    "extractor_version": "..."
  }
}
```

---

## Extraction Statistics

Every output includes comprehensive statistics:

```json
{
  "extraction_statistics": {
    "total_samples": 1,
    "total_runs": 4,

    "field_coverage": {
      "organism": {"count": 1, "percentage": 100.0},
      "tissue": {"count": 1, "percentage": 100.0},
      "strain": {"count": 1, "percentage": 100.0}
    },

    "confidence_distribution": {
      "overall": {
        "high (≥0.8)": 7,
        "medium (0.5-0.8)": 0,
        "low (<0.5)": 0,
        "total_fields": 7
      },
      "by_field": {
        "organism": {"high": 1, "medium": 0, "low": 0, "total": 1},
        "tissue": {"high": 1, "medium": 0, "low": 0, "total": 1}
      }
    },

    "source_distribution": {
      "overall": {"biosample": 7},
      "by_field": {
        "organism": {"biosample": 1},
        "tissue": {"biosample": 1}
      }
    },

    "ontology_mapping": {
      "overall": {
        "mapped": 2,
        "total": 2,
        "percentage": 100.0
      },
      "by_field": {
        "organism": {"mapped": 1, "total": 1, "percentage": 100.0},
        "tissue": {"mapped": 1, "total": 1, "percentage": 100.0}
      }
    },

    "completeness": {
      "complete_samples": 1,
      "percentage": 100.0,
      "criteria": "organism + (tissue OR cell_line)"
    },

    "extraction_methods": {
      "structured_field": 7
    }
  }
}
```

---

## Navigation Guide

Every output includes a navigation guide:

```json
{
  "navigation_guide": {
    "description": "How to access different levels of detail",
    "quick_access": "samples[<id>].quick_view.<field> - Just the values",
    "standard_provenance": "samples[<id>].biological_metadata.<field>.provenance - Source + confidence",
    "detailed_debugging": "samples[<id>].biological_metadata.<field>.extraction_details - Full chain",
    "raw_data": "samples[<id>].raw_biosample_attributes - Original BioSample data"
  }
}
```

---

## Usage Examples

### Example 1: Quick Data Exploration

```python
import json

with open("PRJNA1170270_traceable.json") as f:
    data = json.load(f)

# Quick view of all samples
for sample_id, sample in data["samples"].items():
    qv = sample["quick_view"]
    print(f"{sample_id}: {qv['organism']} - {qv.get('tissue', 'N/A')}")
```

### Example 2: Quality Control Filtering

```python
# Filter samples by confidence
high_quality = {}
for sample_id, sample in data["samples"].items():
    tissue = sample.get("biological_metadata", {}).get("tissue")
    if tissue and tissue["provenance"]["confidence"] >= 0.8:
        high_quality[sample_id] = sample

print(f"High quality: {len(high_quality)}/{len(data['samples'])} samples")
```

### Example 3: Provenance Audit

```python
# Show complete provenance chain for a field
sample = data["samples"]["SRS22569018"]
tissue = sample["biological_metadata"]["tissue"]

print(f"Value: {tissue['value']}")
print(f"\nProvenance:")
print(f"  Source: {tissue['provenance']['source']}")
print(f"  Confidence: {tissue['provenance']['confidence']}")
print(f"  Ontology: {tissue['provenance']['ontology_term']}")

print(f"\nExtraction Details:")
for key, value in tissue["extraction_details"].items():
    print(f"  {key}: {value}")
```

### Example 4: Compare Original vs Processed

```python
# Show transformation from raw to final
sample = data["samples"]["SRS22569018"]

raw_tissue = sample["raw_biosample_attributes"].get("tissue")
final_tissue = sample["biological_metadata"]["tissue"]["value"]
ontology = sample["biological_metadata"]["tissue"]["provenance"]["ontology_term"]

print(f"Original: {raw_tissue}")
print(f"Final: {final_tissue}")
print(f"Ontology: {ontology}")
```

---

## Confidence Score Interpretation

| Total Confidence | Meaning | Interpretation |
|-----------------|---------|----------------|
| ≥ 0.95 | Excellent | Exact match from specific field |
| 0.8 - 0.95 | High | Good synonym match or reliable source |
| 0.5 - 0.8 | Medium | Derived/inferred value, review for critical use |
| < 0.5 | Low | Uncertain extraction, manual review needed |

### Confidence Breakdown

**Field Confidence:**
- 1.0 = Specific field like "tissue_type"
- 0.8 = Generic structured field
- 0.5 = LLM inferred from text
- 0.3 = Extracted from generic description

**Normalization Confidence:**
- 1.0 = Exact ontology match
- 0.95 = Close synonym
- 0.8 = Partial match
- 0.5 = Derived/transformed

---

## Extraction Methods

| Method | Description |
|--------|-------------|
| `structured_field` | Direct extraction from BioSample/SRA fields |
| `llm` | Large language model inference |
| `ner` | Named entity recognition |
| `regex` | Pattern matching |
| `ontology_mapping` | Mapped via ontology lookup |
| `manual` | Manually curated |

---

## Best Practices

### 1. Start with Quick View

```python
# Quick exploration
quick_data = {
    sample_id: sample["quick_view"]
    for sample_id, sample in data["samples"].items()
}
```

### 2. Filter by Confidence for Production

```python
# Production-ready data
def is_high_quality(sample):
    tissue = sample.get("biological_metadata", {}).get("tissue")
    if not tissue:
        return False
    return tissue["provenance"]["confidence"] >= 0.8

prod_samples = {
    sid: s for sid, s in data["samples"].items()
    if is_high_quality(s)
}
```

### 3. Use Extraction Statistics

```python
# Check overall quality
stats = data["extraction_statistics"]
if stats["completeness"]["percentage"] < 70:
    print("Warning: Low completeness, consider LLM enrichment")

if stats["ontology_mapping"]["overall"]["percentage"] < 80:
    print("Warning: Low ontology mapping coverage")
```

### 4. Preserve Full Traceability for Archives

When archiving data or publishing datasets, preserve the complete traceable format for:
- Reproducibility
- Audit trails
- Regulatory compliance
- Quality assessment

---

## File Size Considerations

The traceable format is more verbose than simple value-only formats:

| Format | Size (per sample) | Use Case |
|--------|------------------|----------|
| Quick view only | ~200 bytes | Exploration |
| Standard provenance | ~1-2 KB | Production |
| Full traceable | ~3-5 KB | Archive |

For large-scale batch processing (1000+ samples):
- Individual files: Use full traceable format
- Aggregate summaries: Use extraction statistics
- Data exports: Use quick_view for CSV

---

## Next Steps

- See `scripts/test_traceable_output.py` for working examples
- See `scripts/batch_extract_riboseq_traceable.py` for batch processing
- See `docs/PRACTICAL_LLM_STRATEGY.md` for extraction strategy
