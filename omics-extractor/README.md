# omics-extractor

Generic metadata extraction framework for omics data (RiboSeq, RNA-Seq, etc.)

Automatically extracts, validates, and structures metadata from SRA, GEO, and PubMed with full provenance tracking.

## Quick Start

### Installation

```bash
# Clone the repository
cd omics-extractor

# Install with uv (recommended - fast!)
curl -LsSf https://astral.sh/uv/install.sh | sh
source $HOME/.local/bin/env  # Add uv to PATH
uv sync --all-extras

# Or using pip
pip install -e ".[dev]"
```

### Extract Metadata for Any BioProject

```bash
# Activate virtual environment (if using uv)
source .venv/bin/activate

# Extract complete metadata
omics-extract extract PRJNA1170270

# Output: PRJNA1170270_metadata.json with:
# - Study metadata (title, description, publications)
# - Sample metadata (tissue, cell type, treatment)
# - Run metadata (library strategy, platform, stats)
# - Full provenance for every field

# Validate extracted metadata
omics-extract validate PRJNA1170270_metadata.json

# Extract to specific file
omics-extract extract PRJNA1170270 --output my_metadata.json
```

## Features

### 🎯 End-to-End Metadata Extraction

**One command extracts everything:**
- Study-level metadata from GEO and PubMed
- Sample biological characteristics (tissue, cell type, treatment)
- Run technical details (platform, library strategy, stats)
- Full provenance tracking (source, confidence, extraction method)

### 📦 Data Fetchers

- **BioProject Fetcher** - Study-level project metadata from NCBI BioProject
  - **Primary source for study metadata** - universal coverage (every project has BioProject record)
  - Project title, description, organism, data type, scope
  - **Direct publication links** to PubMed IDs
  - Authoritative project-level information
  - High-confidence extraction (0.95)

- **BioSample Fetcher** - Sample biological metadata from NCBI BioSample
  - **Primary source for sample metadata** - universal coverage (every SRA submission has BioSample)
  - Authoritative, structured, validated attributes
  - Extracts tissue, cell type, strain, genotype, treatment, age, sex
  - High-confidence extraction (0.9-1.0)
  - Full provenance tracking

- **SRA Fetcher** - Sequencing run metadata from NCBI SRA
  - All runs for a BioProject
  - Detailed run metadata (library strategy, platform, sample info)
  - **Extracts BioSample ID and organism** from SRA XML
  - Handles GSE to BioProject conversion

- **GEO Fetcher** - Project and sample metadata from GEO
  - **Supplemental source** for enhanced study descriptions
  - Series (GSE) metadata with project descriptions
  - Sample (GSM) metadata with characteristics
  - Partial coverage (~30-40% of projects)
  - Cross-references to BioProjects and BioSamples

- **PubMed Fetcher** - Publication metadata
  - Abstracts, authors, journal information
  - DOI extraction
  - **PMC ID detection** for full-text access
  - PMC URL generation for open access papers
  - MeSH keywords and affiliations

### 🔬 Schema System with Provenance

Every extracted value includes full provenance:

```json
{
  "tissue": {
    "value": "brain",
    "source": "geo",
    "source_id": "GSM3085447",
    "confidence": 0.9,
    "extraction_method": "structured_field",
    "ontology_term": null
  }
}
```

Models:
- `StudyMetadata` - Project-level information
- `SampleMetadata` - Biological characteristics
- `RunMetadata` / `RiboSeqRunMetadata` - Technical details
- `BaseProvenance` - Tracks every field's origin

## CLI Usage

### Extract Command

```bash
# Basic extraction
omics-extract extract PRJNA1170270

# Specify output file
omics-extract extract PRJNA1170270 --output results.json

# Verbose mode
omics-extract extract PRJNA1170270 --verbose
```

### Validate Command

```bash
# Validate extracted metadata
omics-extract validate PRJNA1170270_metadata.json

# Shows:
# - Missing required fields
# - Low confidence extractions
# - Tissue coverage warnings
```

## Python API Usage

### Complete Extraction

```python
from omics_extractor.extraction.builder import build_project_metadata

# Extract everything for a BioProject
metadata = build_project_metadata("PRJNA1170270")

# Access the results
study = metadata["study"]
samples = metadata["samples"]  # Dict[sample_id -> SampleMetadata]
runs = metadata["runs"]  # List[RunMetadata]

print(f"Study: {study.title.value}")
print(f"Samples: {len(samples)}")
print(f"Runs: {len(runs)}")

# Access sample tissue info
for sample_id, sample in samples.items():
    if sample.tissue:
        print(f"{sample_id}: {sample.tissue.value}")
```

### Individual Fetchers

```python
from omics_extractor.fetchers.bioproject import fetch_bioproject_metadata
from omics_extractor.fetchers.sra import fetch_project_runs, fetch_run_details
from omics_extractor.fetchers.geo import fetch_project_metadata, fetch_sample_metadata
from omics_extractor.fetchers.biosample import fetch_biosample_metadata
from omics_extractor.fetchers.pubmed import fetch_publication_metadata

# Fetch BioProject metadata (primary source for study info)
project = fetch_bioproject_metadata("PRJNA1176138")
print(f"Title: {project.title}")
print(f"Description: {project.description[:100]}...")
print(f"Organism: {project.organism}")
print(f"Publication PMID: {project.publication_id}")

# Fetch SRA runs
runs = fetch_project_runs("PRJNA1170270")
run_meta = fetch_run_details(runs[0])
print(f"BioSample: {run_meta.biosample_id}, Organism: {run_meta.organism}")

# Fetch BioSample metadata (primary source for sample info)
biosample = fetch_biosample_metadata("SAMN44381105")
print(f"Tissue: {biosample.attributes.get('tissue')}")
print(f"Strain: {biosample.attributes.get('strain')}")

# Fetch GEO metadata (supplemental source)
geo_project = fetch_project_metadata("GSE112882")
geo_sample = fetch_sample_metadata("GSM3085447")

# Fetch publication
pub = fetch_publication_metadata("29618526")
print(f"PMC URL: {pub.pmc_url}")  # Link to full text
```

## Testing

Comprehensive test suite with real RiboSeq data:

```bash
# Run all tests
pytest

# Run with coverage report
pytest --cov=src/omics_extractor

# Run specific test module
pytest tests/fetchers/test_sra.py -v

# Run integration tests only
pytest tests/test_integration.py -v
```

### Test Coverage

- **65 tests total** - All passing ✓
- **77% coverage for BioProject fetcher**
- **60% overall code coverage** (91% for fetchers, 98% for schemas)
- Tested with real BioProjects from curated RiboSeq dataset
- Integration tests verify cross-database linking (BioProject ↔ SRA ↔ GEO ↔ BioSample ↔ PubMed)

## Output Format

The CLI generates JSON with complete metadata hierarchy:

```json
{
  "bioproject_id": "PRJNA1170270",
  "extraction_timestamp": "2026-01-03T10:08:33.244444",
  "extractor_version": "0.1.0",
  "study": {
    "bioproject_id": "PRJNA1170270",
    "title": {...},
    "organism": {...},
    "pmid": "...",
    "pmc_id": "..."
  },
  "samples": {
    "SRS22569018": {
      "sample_id": "SRS22569018",
      "organism": {...},
      "tissue": {...},
      "cell_type": {...},
      "treatment": {...}
    }
  },
  "runs": [
    {
      "run_id": "SRR30910340",
      "library_strategy": {...},
      "platform": "ILLUMINA",
      "read_count": 10000000
    }
  ],
  "summary": {
    "total_samples": 1,
    "total_runs": 4,
    "riboseq_runs": 4
  }
}
```

## Key Findings from Real Data

### Multi-Source Metadata Strategy

The system uses a **waterfall approach** for maximum metadata coverage:

**Study-level metadata:**
1. **BioProject (primary)** - Universal coverage, authoritative project info (title, description, organism, publications)
2. **GEO (supplemental)** - Enhanced descriptions when available
3. **PubMed (publications)** - Full publication metadata via PMID links

**Sample-level metadata:**
1. **BioSample (primary)** - Universal coverage, authoritative structured attributes
2. **GEO (supplemental)** - Legacy samples without BioSample
3. **SRA (baseline)** - Organism and run-level technical details when other sources unavailable

**Result:** 100% metadata coverage across all tested projects
- BioProject provides study title, description, organism, and publication links for 100% of projects
- BioSample provides tissue, organism, genotype for all samples

### Library Strategy Inconsistencies

RiboSeq data in SRA has inconsistent labels:
- Correctly labeled: "Ribo-seq"
- Mislabeled as: "RNA-Seq", "OTHER"

**Solution:** Our system auto-detects RiboSeq runs and creates appropriate metadata models.

### Tissue Extraction

Tissue information is **prioritized** in extraction:
- Primary source: BioSample structured attributes
- Checks multiple field variations: tissue → tissue_type → source_name → isolation_source
- High confidence (0.9-1.0) for structured fields
- Falls back to GEO characteristics if BioSample unavailable
- Full provenance tracking from all sources

### PMC Availability

Many publications have PMC IDs for full-text access:
- PMC URLs automatically generated
- Ready for future full-text extraction
- Methods sections can be parsed later

## Project Structure

```
omics-extractor/
├── src/omics_extractor/
│   ├── fetchers/          # Data source fetchers
│   │   ├── bioproject.py  # BioProject fetcher (primary study metadata)
│   │   ├── biosample.py   # BioSample fetcher (primary sample metadata)
│   │   ├── sra.py         # SRA/NCBI fetcher (extracts BioSample ID)
│   │   ├── geo.py         # GEO fetcher (supplemental metadata)
│   │   └── pubmed.py      # PubMed fetcher (publications)
│   ├── schemas/           # Data models
│   │   └── base.py        # BaseProvenance, Study/Sample/Run models
│   ├── extraction/        # Aggregation layer
│   │   └── builder.py     # BioProject/BioSample → GEO → SRA waterfall strategy
│   └── cli.py             # Command-line interface
├── tests/                 # Comprehensive test suite
└── pyproject.toml         # Package configuration
```

## What's Next

Future enhancements:

1. **Ontology System** - Validate terms against NCBITaxon, UBERON, CL, EFO
2. **LLM Extraction** - Extract from abstracts and PMC full text
3. **Django Export** - Map to RiboSeq.Org portal schema
4. **Nextflow Pipeline** - Batch processing for thousands of projects
5. **Human Review Workflow** - Interface for curator validation

## Development

Built with:
- **uv** - Fast Python package manager
- **pytest** - Testing framework
- **pydantic** - Data validation and schemas
- **biopython** - NCBI Entrez API access
- **GEOparse** - GEO SOFT file parsing

### Contributing

```bash
# Install dev dependencies
uv sync --all-extras

# Run tests before committing
pytest

# Check code formatting
ruff check src/
```

## License

See LICENSE file.
