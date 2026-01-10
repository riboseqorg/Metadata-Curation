# Practical LLM Strategy for 70B Models

## The Problem

You're targeting **local 70B models** (Llama 3.3, Qwen 2.5) for production, not Claude 4.5. These models:
- Have smaller context windows (8K-32K vs Claude's 200K)
- Are less powerful at complex reasoning
- Can produce nonsense if overloaded with context
- Cost $0 but require careful prompt engineering

## The Solution: Minimal LLM Usage

**Philosophy**: Use structured extraction for 90%, LLM only for the hard 10%

### Three-Stage Pipeline

```
┌─────────────────────────────────────────────────────────────┐
│ Stage 1: Direct Structured Extraction (No LLM)             │
│ - Extract from BioSample attributes using Python dict      │
│ - Ontology mapping for known terms                         │
│ - Result: 80-90% of fields populated                       │
└─────────────────────────────────────────────────────────────┘
                            ↓
┌─────────────────────────────────────────────────────────────┐
│ Stage 2: Gap Analysis (No LLM)                             │
│ - Check which critical fields are still missing            │
│ - If completeness > 75%, SKIP LLM entirely                 │
│ - If completeness < 75%, proceed to Stage 3                │
└─────────────────────────────────────────────────────────────┘
                            ↓
┌─────────────────────────────────────────────────────────────┐
│ Stage 3: Minimal LLM Gap-Filling (Only when needed)        │
│ - Small prompt: just study title + sample title            │
│ - Ask ONLY for missing critical fields                     │
│ - Keep context under 1000 tokens                           │
└─────────────────────────────────────────────────────────────┘
```

## Implementation

### Stage 1: Direct Extraction (No LLM)

```python
from omics_extractor.extraction.field_mappings import extract_from_attributes
from omics_extractor.ontology.mapper import map_to_ontology

def extract_structured(biosample: dict) -> dict:
    """
    Extract from structured BioSample attributes.
    NO LLM NEEDED - just dict lookup + ontology mapping.
    """
    # 1. Map field names
    standardized = extract_from_attributes(biosample['attributes'])

    # 2. Ontology mapping
    if 'tissue' in standardized:
        ontology_term = map_to_ontology('tissue', standardized['tissue'])
        if ontology_term:
            standardized['tissue_ontology'] = ontology_term
            standardized['tissue_confidence'] = 1.0  # High confidence from structured

    if 'cell_type' in standardized:
        ontology_term = map_to_ontology('cell_type', standardized['cell_type'])
        if ontology_term:
            standardized['cell_type_ontology'] = ontology_term
            standardized['cell_type_confidence'] = 1.0

    return standardized
```

**Result**: Most samples will have `tissue`, `strain`, `organism` from structured fields!

### Stage 2: Gap Analysis (No LLM)

```python
def needs_llm_enrichment(extracted: dict) -> bool:
    """
    Decide if LLM is needed.

    Skip LLM if we have:
    - organism (always in BioSample)
    - tissue OR cell_line
    - strain (for model organisms)
    """
    critical_fields = ['organism', 'tissue', 'strain']

    # Count how many we have
    populated = sum(1 for field in critical_fields if extracted.get(field))
    completeness = populated / len(critical_fields)

    # Skip LLM if we're 75%+ complete
    return completeness < 0.75
```

**Result**: ~70-80% of samples skip LLM entirely!

### Stage 3: Minimal LLM (Only for gaps)

```python
def fill_gaps_with_llm(
    extracted: dict,
    study_title: str,
    sample_title: str
) -> dict:
    """
    Lightweight LLM call - only for missing critical fields.

    Token budget: ~500 input, ~100 output = ~600 total
    Cost: ~$0.003 per sample (vs $0.02 for comprehensive approach)
    """
    # Find what's missing
    missing = []
    if not extracted.get('tissue') and not extracted.get('cell_line'):
        missing.append('tissue or cell_line')
    if not extracted.get('cell_type'):
        missing.append('cell_type')
    if not extracted.get('strain'):
        missing.append('strain')

    if not missing:
        return extracted  # Nothing to do!

    # Minimal prompt
    prompt = f"""Study: {study_title}
Sample: {sample_title}

Already have: {', '.join(extracted.keys())}

Extract ONLY: {', '.join(missing)}

Return JSON with just these fields.
"""

    llm_result = llm.extract(prompt, max_tokens=150)  # Small output

    # Merge with existing (don't override structured fields!)
    for field, value in llm_result.items():
        if field not in extracted:  # Only fill gaps
            extracted[field] = value
            extracted[f'{field}_confidence'] = 0.5  # Lower than structured

    return extracted
```

**Token usage**: 500-700 tokens total (vs 5000-10000 for comprehensive)

## Example Scenarios

### Scenario 1: Good Structured Metadata (80% of cases)
```
BioSample attributes:
  organism: Mus musculus
  strain: C57BL/6
  tissue: testis

Stage 1: Extracts all three ✓
Stage 2: 100% complete → SKIP LLM
Stage 3: Not reached

Result: Perfect extraction, $0 LLM cost
```

### Scenario 2: Partial Structured Metadata (15% of cases)
```
BioSample attributes:
  organism: Mus musculus
  strain: C57BL/6
  sample_title: "spermatocytes from adult mouse testis"

Stage 1: Extracts organism, strain
Stage 2: Missing tissue → 66% complete → CALL LLM
Stage 3: LLM infers tissue=testis, cell_type=spermatocyte

Result: Complete extraction, ~$0.003 LLM cost
```

### Scenario 3: Poor Structured Metadata (5% of cases)
```
BioSample attributes:
  organism: Homo sapiens
  sample_title: "HepG2 cells treated with dexamethasone for 6 hours"

Stage 1: Extracts organism only
Stage 2: Missing tissue, cell_line → 33% complete → CALL LLM
Stage 3: LLM infers:
  - cell_line: HepG2
  - tissue: liver (from HepG2 knowledge)
  - cell_type: hepatocyte
  - treatment: dexamethasone

Result: Good extraction, ~$0.005 LLM cost
```

## Why This Works for 70B Models

**1. Small Context**
- 500-700 tokens vs 5000-10000 tokens
- Fits easily in any context window
- Less chance of "context confusion"

**2. Simple Task**
- "Extract tissue and cell_type" vs "Analyze entire experimental design"
- Clear, focused instructions
- Higher accuracy from smaller models

**3. Fallback Safety**
- Structured extraction is deterministic
- LLM only supplements, never overrides
- If LLM fails, you still have structured data

**4. Cost Effective**
- 70-80% of samples: $0 (no LLM)
- 15-20% of samples: ~$0.003 (minimal LLM)
- 5% of samples: ~$0.005 (more LLM help)
- **Average: ~$0.001 per sample** (50x cheaper than comprehensive)

## Comparison: Minimal vs Comprehensive

| Aspect | Comprehensive LLM | Minimal LLM (This) |
|--------|-------------------|-------------------|
| **Context size** | 5000-10000 tokens | 500-700 tokens |
| **LLM usage** | 100% of samples | 20-30% of samples |
| **Avg cost/sample** | $0.02 | $0.001 |
| **Works on 70B?** | Risky (context overload) | ✅ Yes (small prompts) |
| **Accuracy** | High (if model can handle) | High (structured + LLM) |
| **Speed** | Slow (always LLM) | Fast (mostly direct) |

## Migration Path

### Step 1: Implement Simple Field Mapping (Done!)
- ✅ `field_mappings.py` created
- Simple Python dict, no CSV parsing
- Easy to maintain

### Step 2: Update Baseline Extraction
- Use `extract_from_attributes()` instead of CSV approach
- Keep ontology mapping
- Test on current datasets

### Step 3: Add Minimal LLM Layer
- Implement gap analysis logic
- Add minimal prompt templates
- Test on samples with poor metadata

### Step 4: Benchmark 70B Models
- Test Llama 3.3 70B, Qwen 2.5 72B
- Compare with Claude baseline
- Optimize prompts for best 70B model

## Recommended Next Steps

1. **Test the simple field mapping**
   ```bash
   python -c "from omics_extractor.extraction.field_mappings import extract_from_attributes; \
              print(extract_from_attributes({'tissue_type': 'liver', 'strain': 'C57BL/6'}))"
   ```

2. **Integrate into baseline extraction**
   - Replace CSV-based column mapping
   - Use `field_mappings.py` instead

3. **Add gap-filling logic**
   - Implement `needs_llm_enrichment()`
   - Add minimal LLM prompts

4. **Benchmark on real data**
   - Test on RiboSeq projects
   - Measure: completeness without LLM vs with LLM
   - Compare costs

Would you like me to implement the baseline extraction update next?
