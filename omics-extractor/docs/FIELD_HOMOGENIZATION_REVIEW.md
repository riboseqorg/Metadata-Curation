# Field Homogenization Review: Current Approach vs LLM-Based Approach

## Current Manual CSV-Based Approach

### How It Works

The current system uses manually curated CSV files to map column names to standardized fields:

**1. Core.csv** - Main field mappings (11 categories)
```csv
Category,Main Name,All Names,Associated Columns
Cell Lines,CELL_LINE,"cell line, cell_line, cell-line, cells...",""cell type, cell_type..."
Tissue,TISSUE,"tissue, tissue_type, organ, organism part...",""cell type, cell_type..."
```

**2. Content.csv** - Value-based mappings (regex patterns)
```csv
Category,Column,Main Name,All Names
Library Names,LIBRARYTYPE,RNA,"rna-seq, RNA-Seq, ^rna_, _RNA$..."
Library Names,LIBRARYTYPE,RFP,"RFP, ribo-seq, Riboseq, ^Ribo_..."
```

**3. Column_value_mapping.csv** - Cross-field lookups
```csv
Column,Mapping
CELL_LINE,"CELL_LINE, TISSUE"
TISSUE,"CELL_LINE, TISSUE"
```

### Problems with This Approach

**1. Maintenance Burden**
- Manually adding every possible column name variant
- Regular expressions are brittle (e.g., `^rna_`, `_RNA$`, `.rna$`)
- New studies introduce new naming conventions not in the CSV

**2. Cross-Field Logic Is Hardcoded**
```r
# Example: Look in both CELL_LINE and TISSUE columns to populate TISSUE field
# This mapping is in Column_value_mapping.csv but the logic is hardcoded in R
```

**3. Context-Insensitive**
- "treatment" could mean:
  - Drug treatment (inhibitor)
  - Temperature treatment (experimental condition)
  - Time treatment (timepoint)
- CSV can't distinguish without context

**4. Limited to Exact Matches**
- If a study uses "cell culture type" instead of "cell line", it won't match
- Requires continuous manual updates to CSV

**5. No Reasoning**
- Can't infer "liver" from "hepatocyte"
- Can't understand "spermatocyte" is a cell type in testis tissue
- Can't resolve conflicts (e.g., if TISSUE says "liver" but CELL_LINE says "HepG2")

## LLM-Based Approach: What We're Building

### How It Works

**Stage 1: Structural Normalization (Keep Simple)**
```python
# Just basic field name standardization
# organism_name → organism
# tissue_type → tissue
# dev_stage → developmental_stage

FIELD_SYNONYMS = {
    "organism": ["organism_name", "species", "organism name"],
    "tissue": ["tissue_type", "organ", "organism part", "tissue origin"],
    "cell_type": ["cell_type", "celltype", "cell type"],
    # ... etc
}
```

**Stage 2: LLM Extraction with Context**
```python
# Give LLM ALL available metadata
prompt = f"""
Study: {study_title}
Description: {study_description}

Sample fields:
- organism: Mus musculus
- tissue_type: testis
- cell_line: [not specified]
- sample_title: spermatocytes from Mus musculus for ribo-seq
- treatment: [not specified]
- developmental_stage: adult

Extract standardized metadata for this sample.
"""

# LLM can reason:
# - "spermatocytes" in title → cell_type = spermatocyte
# - testis tissue + spermatocyte → biologically consistent
# - No cell line mentioned, and it's tissue sample → cell_line = null
```

**Stage 3: Ontology Mapping (Deterministic)**
```python
# After LLM extracts free text, map to ontologies
tissue = "testis" → UBERON:0000473
cell_type = "spermatocyte" → CL:0000017
```

### Advantages

**1. Context-Aware Reasoning**
```
Q: Is "HepG2" a tissue or cell line?
LLM: Cell line (human liver cancer cell line)

Q: If tissue=liver and cell_line=HepG2, what's the cell_type?
LLM: hepatocyte (liver parenchymal cell)
```

**2. Handles Variations Automatically**
```
"cell culture type" → Understands this means cell line
"organism part" → Understands this could be tissue
"subjected to" → Understands this could be treatment
```

**3. Biological Knowledge**
```
"spermatocyte" → Knows this is a cell type in testis
"hepatocyte" → Knows this is in liver
"T cell" → Knows this is from blood/spleen/lymph tissue
```

**4. Conflict Resolution**
```
If field X says "liver" and field Y says "kidney":
LLM can look at study description to determine which is correct
```

**5. Cross-Field Inference**
```
Given:
- tissue: not specified
- cell_line: HEK293

LLM infers:
- tissue: kidney (HEK293 = Human Embryonic Kidney)
- organism: Homo sapiens
- cell_type: epithelial
```

## Recommended Hybrid Approach

### Step 1: Simple Field Name Mapping (Deterministic)
Use a simple Python dict instead of CSV:

```python
# omics_extractor/extraction/field_synonyms.py
FIELD_SYNONYMS = {
    "organism": [
        "organism", "organism_name", "species", "organism name",
        "scientific_name", "scientific name", "latin name"
    ],
    "tissue": [
        "tissue", "tissue_type", "tissue type", "organ",
        "organism part", "body site", "tissue origin",
        "tissue source", "tissue lineage", "tissue subtype"
    ],
    "cell_type": [
        "cell_type", "cell type", "celltype", "cell_subtype",
        "cell subtype"
    ],
    "cell_line": [
        "cell_line", "cell line", "cell-line", "cell lline",
        "culture_collection", "cell line id", "cell line background"
    ],
    "strain": [
        "strain", "strain_name", "strain name", "cultivar",
        "ecotype", "breed", "substrain"
    ],
    "treatment": [
        "treatment", "compound", "drug", "inhibitor",
        "compound treatment", "drug treatment", "treatment.1",
        "cell treatment", "sample treatment"
    ],
    "developmental_stage": [
        "dev_stage", "developmental stage", "development stage",
        "dev stage", "growth stage", "stage"
    ],
    "age": [
        "age", "age_value", "developmental_age", "growth time"
    ],
    "sex": [
        "sex", "gender", "mating_type"
    ],
}

def normalize_field_name(raw_field_name: str) -> Optional[str]:
    """Map raw field name to standardized name."""
    raw_lower = raw_field_name.lower().strip()

    for standard_name, synonyms in FIELD_SYNONYMS.items():
        if raw_lower in [s.lower() for s in synonyms]:
            return standard_name

    return None
```

**Advantages:**
- ✅ Fast, deterministic
- ✅ Easy to maintain (Python dict vs CSV parsing)
- ✅ Version controlled
- ✅ Handles 90% of field name variations

### Step 2: LLM Populates Fields with Context

Give LLM ALL the raw fields (after normalization):

```python
def extract_with_llm_context(sample_data: dict, study_context: dict):
    """
    Extract metadata using LLM with full context.

    Args:
        sample_data: All raw fields from this sample
        study_context: Study title, description, etc.
    """

    # Build comprehensive prompt
    prompt = f"""
    Study: {study_context['title']}
    Description: {study_context['description']}

    Sample Raw Metadata:
    {json.dumps(sample_data, indent=2)}

    Task: Extract standardized metadata fields.
    - Look across ALL fields to find the best value for each target field
    - If multiple fields contain relevant info, use biological context to choose
    - Infer values when possible (e.g., cell_type from cell_line)

    Target fields to extract:
    - organism: species name (e.g., Homo sapiens, Mus musculus)
    - tissue: tissue/organ type (e.g., liver, brain, blood)
    - cell_type: specific cell type (e.g., hepatocyte, neuron)
    - cell_line: cell line name (e.g., HeLa, HEK293) or null if tissue sample
    - strain: strain/cultivar (e.g., C57BL/6, Columbia)
    - treatment: drug/compound treatment (e.g., dexamethasone, DMSO)
    - developmental_stage: life stage (e.g., adult, embryonic day 10)
    - age: age value with units (e.g., 8 weeks, postnatal day 3)
    - sex: male, female, or null

    Return JSON with extracted values and confidence scores.
    """

    return llm.extract(prompt)
```

**Advantages:**
- ✅ Sees ALL raw fields at once (not just pre-mapped ones)
- ✅ Can reason about which field to use if multiple contain same info
- ✅ Can infer missing values from context
- ✅ Handles novel field names automatically

### Step 3: Ontology Mapping (Deterministic)

After LLM extraction, map free text to ontology terms:

```python
from omics_extractor.ontology import OntologyMapper

mapper = OntologyMapper()

# LLM extracted: tissue = "testis"
ontology_term = mapper.map_tissue("testis")
# Returns: UBERON:0000473

# Confidence scoring
field_confidence = 0.5  # LLM source
normalization_confidence = 0.95  # High ontology match
total_confidence = field_confidence * normalization_confidence
```

## Comparison: Current vs LLM Approach

| Aspect | Current (CSV) | LLM-Based |
|--------|---------------|-----------|
| **Field name mapping** | Manual CSV maintenance | Simple Python dict |
| **Value extraction** | Regex patterns | Contextual reasoning |
| **Cross-field logic** | Hardcoded in R | Natural language instructions |
| **Novel field names** | Fails, needs CSV update | Handles automatically |
| **Biological inference** | None | Yes (cell_type from tissue, etc.) |
| **Conflict resolution** | First match wins | Contextual reasoning |
| **Maintenance** | High (continuous CSV updates) | Low (LLM generalizes) |
| **Explainability** | Limited (which regex matched?) | Good (LLM explains reasoning) |
| **Cost** | Free | ~$0.01-0.03 per sample |

## Migration Strategy

### Phase 1: Keep CSV for Field Name Mapping (Short Term)
- Convert CSV to Python dict for easier maintenance
- Use for initial field normalization only
- LLM handles value extraction

### Phase 2: LLM-Only with Few-Shot Examples (Long Term)
- Give LLM example mappings instead of exhaustive lists
- Let LLM generalize to new field names

```python
prompt = f"""
Examples of field name mappings:
- "tissue_type" → tissue
- "organism part" → tissue
- "cell line id" → cell_line
- "dev_stage" → developmental_stage

Raw fields from this sample:
{sample_fields}

Map each field to a standard name, then extract values.
"""
```

## Recommendation

**Start with Hybrid:**
1. Simple Python dict for field name normalization (replaces CSV)
2. LLM for value extraction with full context (replaces regex)
3. Deterministic ontology mapping for controlled terms

**This gives you:**
- ✅ Lower maintenance than CSV approach
- ✅ Better handling of novel field names
- ✅ Biological reasoning and inference
- ✅ Context-aware extraction
- ✅ Minimal cost (~$0.01-0.03 per sample)

**Eliminate:**
- ❌ Manual CSV maintenance
- ❌ Brittle regex patterns
- ❌ Hardcoded cross-field logic
- ❌ Context-insensitive matching
