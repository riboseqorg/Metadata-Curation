"""Base schema definitions for metadata extraction."""

from typing import Optional, Dict, Any, List
from pydantic import BaseModel, Field, field_validator
from datetime import datetime


class BaseProvenance(BaseModel):
    """
    Track the origin and validation of a metadata field.

    This model wraps every extracted metadata value with:
    - Where it came from (source)
    - How confident we are (two-factor confidence: field × normalization)
    - Ontology mapping (if applicable)
    - Original raw text (for debugging)

    Confidence Scoring:
    - Final confidence = field_confidence × normalization_confidence
    - field_confidence: How specific the source field is (1.0 for exact match, lower for generic fields)
    - normalization_confidence: How much transformation was needed (1.0 for exact match, lower for inferred)
    """

    value: str = Field(description="The extracted/normalized value")
    source: str = Field(description="Data source: 'sra', 'geo', 'biosample', 'bioproject', 'pubmed', 'llm', 'manual'")
    source_id: str = Field(description="Source record ID (e.g., GSE123456, PMID:123456)")
    confidence: float = Field(ge=0.0, le=1.0, description="Final confidence score (field × normalization)")

    # Confidence components (optional, for transparency)
    field_confidence: Optional[float] = Field(None, ge=0.0, le=1.0, description="Field specificity confidence")
    normalization_confidence: Optional[float] = Field(None, ge=0.0, le=1.0, description="Normalization quality confidence")

    # Optional fields
    extracted_text: Optional[str] = Field(None, description="Original text before normalization")
    ontology_term: Optional[str] = Field(None, description="Ontology ID (e.g., 'UBERON:0000955')")
    ontology_label: Optional[str] = Field(None, description="Canonical ontology label")
    extraction_method: Optional[str] = Field(None, description="How extracted: 'structured_field', 'llm', 'ner', 'regex'")
    extraction_timestamp: Optional[datetime] = Field(None, description="When extracted")
    notes: Optional[str] = Field(None, description="Any additional notes or warnings")

    @field_validator('confidence')
    @classmethod
    def validate_confidence(cls, v):
        """Ensure confidence is between 0 and 1."""
        if not 0.0 <= v <= 1.0:
            raise ValueError(f"Confidence must be between 0.0 and 1.0, got {v}")
        return v

    def model_dump_simple(self) -> str:
        """Return just the value for simple display."""
        return self.value

    def has_ontology(self) -> bool:
        """Check if this field has an ontology mapping."""
        return self.ontology_term is not None

    def is_high_confidence(self, threshold: float = 0.8) -> bool:
        """Check if confidence exceeds threshold."""
        return self.confidence >= threshold


class StudyMetadata(BaseModel):
    """
    Study/Project level metadata.

    Maps to the portal's Study model.
    Represents a BioProject with overall project information.
    """

    # Required identifiers
    bioproject_id: str = Field(description="BioProject accession (PRJNA...)")

    # Core metadata with provenance
    title: BaseProvenance = Field(description="Study title")
    description: BaseProvenance = Field(description="Study description/summary")
    organism: BaseProvenance = Field(description="Primary organism studied (with NCBITaxon)")

    # Optional identifiers
    gse_id: Optional[str] = Field(None, description="GEO Series accession (GSE...)")
    sra_study_id: Optional[str] = Field(None, description="SRA Study ID")

    # Publication metadata
    pmid: Optional[str] = Field(None, description="PubMed ID")
    pmc_id: Optional[str] = Field(None, description="PubMed Central ID")
    doi: Optional[BaseProvenance] = Field(None, description="Publication DOI")
    authors: Optional[List[str]] = Field(None, description="Author list")
    publication_title: Optional[BaseProvenance] = Field(None, description="Paper title")
    journal: Optional[BaseProvenance] = Field(None, description="Journal name")
    publication_date: Optional[str] = Field(None, description="Publication date")
    paper_abstract: Optional[BaseProvenance] = Field(None, description="Paper abstract")

    # Metadata about the metadata
    extraction_date: Optional[datetime] = Field(None, description="When metadata was extracted")
    extractor_version: Optional[str] = Field(None, description="Version of extractor used")

    # Custom fields (for OpenColumns-like extensibility)
    custom_fields: Optional[Dict[str, Any]] = Field(None, description="Additional custom metadata")


class SampleMetadata(BaseModel):
    """
    Sample/Experiment level metadata.

    Represents biological sample characteristics.
    Multiple runs can come from the same sample.
    """

    # Required identifiers
    sample_id: str = Field(description="Sample accession (SRS/ERS/DRS or SAMN)")
    bioproject_id: str = Field(description="Parent BioProject")

    # Core biological metadata with provenance
    organism: BaseProvenance = Field(description="Organism (NCBITaxon)")

    # Biological context (all optional with provenance)
    tissue: Optional[BaseProvenance] = Field(None, description="Tissue type (UBERON)")
    cell_line: Optional[BaseProvenance] = Field(None, description="Cell line name (CLO)")
    cell_type: Optional[BaseProvenance] = Field(None, description="Cell type (CL)")
    developmental_stage: Optional[BaseProvenance] = Field(None, description="Dev stage (EFO)")
    strain: Optional[BaseProvenance] = Field(None, description="Strain/variety")
    genotype: Optional[BaseProvenance] = Field(None, description="Genetic background")
    sex: Optional[BaseProvenance] = Field(None, description="Biological sex")
    age: Optional[BaseProvenance] = Field(None, description="Age")

    # Experimental conditions (optional with provenance)
    condition: Optional[BaseProvenance] = Field(None, description="Experimental condition")
    treatment: Optional[BaseProvenance] = Field(None, description="Treatment applied (CHEBI/EFO)")
    timepoint: Optional[BaseProvenance] = Field(None, description="Time point")
    replicate: Optional[BaseProvenance] = Field(None, description="Biological replicate")
    batch: Optional[BaseProvenance] = Field(None, description="Batch number")

    # Disease/perturbation (optional with provenance)
    disease: Optional[BaseProvenance] = Field(None, description="Disease state (MONDO/DOID)")
    stress: Optional[BaseProvenance] = Field(None, description="Stress condition")
    temperature: Optional[BaseProvenance] = Field(None, description="Temperature")
    growth_condition: Optional[BaseProvenance] = Field(None, description="Growth conditions")

    # Optional identifiers
    gsm_id: Optional[str] = Field(None, description="GEO Sample accession (GSM...)")
    biosample_id: Optional[str] = Field(None, description="BioSample ID (SAMN...)")

    # Sample description
    sample_title: Optional[BaseProvenance] = Field(None, description="Sample title")
    sample_description: Optional[BaseProvenance] = Field(None, description="Sample description")

    # Raw characteristics from GEO (before parsing into structured fields)
    raw_characteristics: Optional[Dict[str, str]] = Field(None, description="Unparsed GEO characteristics")

    # Custom fields
    custom_fields: Optional[Dict[str, Any]] = Field(None, description="Additional custom metadata")


class RunMetadata(BaseModel):
    """
    Run level metadata.

    Represents a single sequencing run with technical details.
    This is the most granular level - actual data files.
    """

    # Required identifiers
    run_id: str = Field(description="Run accession (SRR/ERR/DRR)")
    sample_id: str = Field(description="Parent sample")
    bioproject_id: str = Field(description="Parent BioProject")
    experiment_id: str = Field(description="Experiment accession (SRX/ERX/DRX)")
    biosample_id: Optional[str] = Field(None, description="BioSample accession (SAMN/SAMEA/SAMD)")
    organism: Optional[str] = Field(None, description="Organism name")

    # Library metadata (with provenance where applicable)
    library_strategy: BaseProvenance = Field(description="Library strategy (RNA-Seq, RiboSeq, etc)")
    library_source: str = Field(description="Library source (TRANSCRIPTOMIC, etc)")
    library_selection: str = Field(description="Library selection method")
    library_layout: str = Field(description="SINGLE or PAIRED")
    library_name: Optional[str] = Field(None, description="Library name")

    # Sequencing platform
    platform: str = Field(description="Sequencing platform (ILLUMINA, etc)")
    instrument_model: Optional[str] = Field(None, description="Instrument model")

    # Run statistics
    read_count: Optional[int] = Field(None, description="Number of reads")
    base_count: Optional[int] = Field(None, description="Number of bases")
    avg_length: Optional[float] = Field(None, description="Average read length")

    # Dates
    run_date: Optional[str] = Field(None, description="Run publish date")

    # Custom fields for extensibility
    custom_fields: Optional[Dict[str, Any]] = Field(None, description="Additional custom metadata")


class RiboSeqRunMetadata(RunMetadata):
    """
    RiboSeq-specific run metadata.

    Extends RunMetadata with RiboSeq protocol details.
    """

    # Override library_strategy to enforce RiboSeq
    library_strategy: BaseProvenance = Field(
        description="Must be RiboSeq/Ribo-seq/ribosome profiling"
    )

    # RiboSeq protocol specifics (all optional with provenance)
    inhibitor: Optional[BaseProvenance] = Field(None, description="Translation inhibitor (CHEBI)")
    nuclease: Optional[BaseProvenance] = Field(None, description="Nuclease used (RNase I, MNase)")
    fraction: Optional[BaseProvenance] = Field(None, description="Fraction (monosome, polysome, etc)")

    # Footprint characteristics
    footprint_min_length: Optional[int] = Field(None, description="Min footprint length (nt)")
    footprint_max_length: Optional[int] = Field(None, description="Max footprint length (nt)")

    # Library prep details
    umi: Optional[BaseProvenance] = Field(None, description="UMI used")
    adapter: Optional[BaseProvenance] = Field(None, description="Adapter sequence")
    rrna_depletion: Optional[BaseProvenance] = Field(None, description="rRNA depletion method")
    kit: Optional[BaseProvenance] = Field(None, description="Library prep kit")

    # Processing hints
    monosome_purification: Optional[BaseProvenance] = Field(None, description="Monosome purification method")


# Field to ontology mapping configuration
FIELD_ONTOLOGY_MAPPING = {
    "organism": {
        "ontology": "NCBITaxon",
        "required": True,
        "description": "Organism taxonomy"
    },
    "tissue": {
        "ontology": "UBERON",
        "required": False,
        "description": "Anatomical structure"
    },
    "cell_type": {
        "ontology": "CL",
        "required": False,
        "description": "Cell type"
    },
    "cell_line": {
        "ontology": "CLO",
        "required": False,
        "description": "Cell line"
    },
    "developmental_stage": {
        "ontology": "EFO",
        "required": False,
        "description": "Developmental stage"
    },
    "treatment": {
        "ontology": ["CHEBI", "EFO"],
        "required": False,
        "description": "Chemical or treatment"
    },
    "inhibitor": {
        "ontology": "CHEBI",
        "required": False,
        "description": "Translation inhibitor"
    },
    "disease": {
        "ontology": ["MONDO", "DOID"],
        "required": False,
        "description": "Disease state"
    },
}
