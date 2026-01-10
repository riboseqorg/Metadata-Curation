# Two-Stage LLM Approach for Metadata Extraction

## Philosophy: Give the LLM Maximum Context

Rather than limiting what we send to the LLM to save tokens, we take the approach of **giving comprehensive context** so the LLM can do proper characterization and analysis.

## Two Stages

### Stage 1: Project-Level Analysis
**Goal**: Understand the experimental design as a whole

**Input to LLM**:
- Study title and description
- **ALL samples** (not just one)
- **ALL metadata fields** for each sample:
  - Structured fields (organism, strain, tissue, etc.)
  - Sample titles and descriptions
  - Run information (library strategy, platform, etc.)
  - Any baseline-extracted metadata

**LLM Tasks**:
1. Identify experimental variables (treatment, timepoint, tissue, etc.)
2. Group biological/technical replicates
3. Map relationships between samples:
   - Paired assays (Ribo-Seq ↔ RNA-Seq from same sample)
   - Time series / dose-response experiments
   - Multi-tissue/multi-condition studies
4. Characterize experimental design ("3x2 factorial with 4 replicates")

**Output**:
```python
ProjectAnalysisResult(
    experimental_variables=["tissue", "treatment"],
    replicate_groups=[
        ReplicateGroup(
            group_id="control_liver_replicates",
            sample_ids=["SRS123", "SRS124", "SRS125"],
            conditions={"tissue": "liver", "treatment": "control"},
            replicate_type="biological"
        )
    ],
    relationships=[
        SampleRelationship(
            sample_id_1="SRS123",
            sample_id_2="SRS456",
            relationship_type="paired_assay",
            description="Ribo-Seq and RNA-Seq from same liver sample"
        )
    ],
    design_summary="Paired Ribo-Seq/RNA-Seq from liver and heart tissue with and without drug treatment. 3 biological replicates per condition.",
    confidence=0.9
)
```

**Why This Matters**:
- You can't identify replicates by looking at ONE sample
- You can't map paired experiments without seeing ALL samples
- You can't characterize experimental design without the full picture

### Stage 2: Sample-Level Extraction
**Goal**: Extract specific metadata fields for each sample

**Input to LLM**:
- Study context (title, description)
- **This specific sample's** metadata:
  - Sample title and description
  - Structured fields already extracted (organism, strain, etc.)
  - Any baseline findings (from ontology mapping, etc.)
- **Results from Stage 1** (experimental design context)

**LLM Tasks**:
1. Extract missing metadata fields (tissue, cell_type, treatment, etc.)
2. Use experimental design context to make better inferences
3. Provide confidence scores
4. Don't duplicate what baseline extraction already found

**Output**:
```python
LLMExtractionResult(
    tissue="liver",
    cell_type="hepatocyte",
    treatment="dexamethasone",
    strain="C57BL/6",  # Already found by baseline, confirmed by LLM
    confidence={
        "tissue": 0.95,
        "cell_type": 0.85,
        "treatment": 0.90,
    },
    reasoning="Sample is from liver hepatocytes treated with dexamethasone based on sample title and experimental design analysis showing this is a liver drug treatment study."
)
```

## Why Two Stages?

### Why Not Just One Big Prompt?
- **Stage 1 is expensive but infrequent**: Run once per project (e.g., 1 project with 50 samples = 1 LLM call)
- **Stage 2 is cheaper and parallelizable**: Run per sample (50 samples = 50 parallel LLM calls)
- **Different tasks need different context**:
  - Finding replicates requires ALL samples
  - Extracting tissue for one sample doesn't need to see other samples

### Example Workflow

```python
# 1. Baseline extraction (structured fields + ontology)
metadata = build_project_metadata("PRJNA12345")
study = metadata["study"]
samples = metadata["samples"]

# 2. Stage 1: Project-level analysis (1 LLM call)
project_analysis = analyze_project_design(study, samples, llm_provider)
# Now we know: experimental variables, replicates, relationships

# 3. Stage 2: Sample-level enrichment (N parallel LLM calls)
for sample_id, sample in samples.items():
    enriched = enrich_sample_metadata(
        sample=sample,
        study=study,
        project_analysis=project_analysis,  # Context from Stage 1
        llm_provider=llm_provider
    )
```

## Benefits of Comprehensive Context

### 1. Better Inference
With full context, LLM knows:
- "This study has liver, heart, and kidney samples - so if sample title says 'cardiac tissue', that's definitely heart"
- "All samples have Ribo-Seq and RNA-Seq runs - so this SRR with RNA-Seq is likely paired with the Ribo-Seq from same sample"

### 2. Consistency
Project-level analysis ensures:
- Consistent terminology across samples ("treatment" vs "drug" vs "compound")
- Proper replicate grouping (biological vs technical)
- Unified experimental design understanding

### 3. Validation
Stage 1 can catch issues:
- "Sample 5 looks like an outlier - different organism?"
- "These 3 samples have identical metadata - likely technical replicates"

### 4. Rich Metadata
Beyond simple field extraction:
- Replicate relationships → Enable proper statistical analysis
- Paired assays → Enable integrated multi-omics analysis
- Experimental design → Enable proper experimental factor analysis

## Token Costs

### Is This Expensive?

**Stage 1 (Project-level)**:
- Input: ~5000-10000 tokens (study + all samples)
- Output: ~500-1000 tokens (structured analysis)
- Cost: ~$0.05-0.15 per project (Claude 4.5 Sonnet)
- Frequency: Once per project

**Stage 2 (Sample-level)**:
- Input: ~1000-2000 tokens per sample
- Output: ~200-500 tokens per sample
- Cost: ~$0.01-0.03 per sample
- Frequency: Once per sample

**Example Project (50 samples)**:
- Stage 1: $0.10 (once)
- Stage 2: $0.50-1.50 (50 × $0.01-0.03)
- **Total: $0.60-1.60 per project**

For comparison:
- Manual curation: Hours of human time per project
- Local 70B model: $0 but requires A100 GPU

**Verdict**: Token costs are negligible compared to value gained!

## Ontology Mapping

You asked about ontology mapping with LLM - two approaches:

### Approach A: Separate Ontology Tool
1. LLM extracts free text ("liver", "hepatocyte")
2. Separate ontology mapper finds controlled terms:
   - "liver" → UBERON:0002107
   - "hepatocyte" → CL:0000182
3. Confidence score combines LLM confidence × ontology match confidence

**Pros**: Deterministic ontology mapping, easier to debug
**Cons**: Two-step process

### Approach B: LLM Does Ontology Mapping
Give LLM the ontology:
```
"Extract tissue and provide ontology term.
Available tissue ontology terms:
- liver: UBERON:0002107
- kidney: UBERON:0002113
- heart: UBERON:0000948
..."
```

**Pros**: Single step, LLM can reason about best match
**Cons**: Need to provide ontology in prompt (token cost), LLM might hallucinate IDs

**Recommendation**: Start with Approach A (separate tools), consider B for future

## Demo

Run the project analysis demo:
```bash
export ANTHROPIC_API_KEY=your-key
python scripts/demo_project_analysis.py
```

This will show comprehensive project-level analysis with full context!
