"""Build complete metadata by aggregating data from multiple sources."""

from typing import Dict, List, Optional
from datetime import datetime

from ..fetchers.sra import fetch_project_runs, fetch_run_details
from ..fetchers.geo import fetch_project_metadata, fetch_sample_metadata
from ..fetchers.pubmed import fetch_publication_metadata, search_pubmed_for_project
from ..fetchers.biosample import fetch_biosample_metadata
from ..fetchers.bioproject import fetch_bioproject_metadata
from ..schemas.base import (
    BaseProvenance,
    StudyMetadata,
    SampleMetadata,
    RunMetadata,
    RiboSeqRunMetadata,
)
from ..ontologies.mapper import (
    normalize_tissue,
    normalize_cell_type,
    normalize_organism,
)
from .field_mappings import extract_from_attributes
from .enhanced_extractor import extract_sample_with_field_mappings


def build_study_metadata(bioproject_id: str, gse_id: Optional[str] = None) -> StudyMetadata:
    """
    Build study-level metadata from BioProject, GEO, and PubMed.

    Extraction priority: BioProject (primary, universal coverage) → GEO (supplemental) → PubMed (publications)

    Args:
        bioproject_id: BioProject accession (PRJNA...)
        gse_id: Optional GEO series ID (will search if not provided)

    Returns:
        StudyMetadata with provenance
    """
    # Fetch BioProject metadata first (primary source - universal coverage)
    bioproject_meta = None
    try:
        bioproject_meta = fetch_bioproject_metadata(bioproject_id)
    except ValueError:
        pass  # BioProject not available

    # Try to fetch GEO metadata if we have GSE ID (supplemental source)
    geo_meta = None
    if gse_id:
        try:
            geo_meta = fetch_project_metadata(gse_id)
        except ValueError:
            pass  # GEO not available

    # Build study metadata prioritizing BioProject
    if bioproject_meta:
        study = StudyMetadata(
            bioproject_id=bioproject_id,
            gse_id=gse_id,
            title=BaseProvenance(
                value=bioproject_meta.title or bioproject_id,
                source="bioproject",
                source_id=bioproject_id,
                confidence=0.95 if bioproject_meta.title else 0.5,
                extraction_method="structured_field" if bioproject_meta.title else "fallback",
                notes="Using BioProject ID as placeholder" if not bioproject_meta.title else None,
            ),
            description=BaseProvenance(
                value=bioproject_meta.description or "",
                source="bioproject",
                source_id=bioproject_id,
                confidence=0.9 if bioproject_meta.description else 0.0,
                extraction_method="structured_field" if bioproject_meta.description else None,
            ),
            organism=BaseProvenance(
                value=bioproject_meta.organism or "unknown",
                source="bioproject",
                source_id=bioproject_id,
                confidence=0.95 if bioproject_meta.organism else 0.0,
                extraction_method="structured_field" if bioproject_meta.organism else None,
            ),
        )

        # Enhance with GEO metadata if available (GEO descriptions sometimes better)
        if geo_meta:
            # Use GEO description if BioProject description is missing or GEO has more details
            if not bioproject_meta.description or (geo_meta.summary and len(geo_meta.summary) > len(bioproject_meta.description or "")):
                study.description = BaseProvenance(
                    value=geo_meta.summary,
                    source="geo",
                    source_id=geo_meta.gse_id,
                    confidence=1.0,
                    extraction_method="structured_field",
                )

    elif geo_meta:
        # Fall back to GEO if BioProject not available
        study = StudyMetadata(
            bioproject_id=bioproject_id,
            gse_id=geo_meta.gse_id,
            title=BaseProvenance(
                value=geo_meta.title,
                source="geo",
                source_id=geo_meta.gse_id,
                confidence=1.0,
                extraction_method="structured_field",
            ),
            description=BaseProvenance(
                value=geo_meta.summary,
                source="geo",
                source_id=geo_meta.gse_id,
                confidence=1.0,
                extraction_method="structured_field",
            ),
            organism=BaseProvenance(
                value=geo_meta.organism,
                source="geo",
                source_id=geo_meta.gse_id,
                confidence=1.0,
                extraction_method="structured_field",
            ),
        )
    else:
        # Minimal metadata if neither BioProject nor GEO available
        study = StudyMetadata(
            bioproject_id=bioproject_id,
            gse_id=gse_id,
            title=BaseProvenance(
                value=bioproject_id,
                source="bioproject",
                source_id=bioproject_id,
                confidence=0.5,
                notes="No metadata available, using BioProject ID as placeholder",
            ),
            description=BaseProvenance(
                value="",
                source="bioproject",
                source_id=bioproject_id,
                confidence=0.0,
            ),
            organism=BaseProvenance(
                value="unknown",
                source="bioproject",
                source_id=bioproject_id,
                confidence=0.0,
            ),
        )

    # Add publication metadata from BioProject or GEO
    pubmed_id = None
    if bioproject_meta and bioproject_meta.publication_id:
        pubmed_id = bioproject_meta.publication_id
    elif geo_meta and geo_meta.pubmed_id:
        pubmed_id = geo_meta.pubmed_id

    if pubmed_id:
        try:
            pub_meta = fetch_publication_metadata(pubmed_id)
            study.pmid = pub_meta.pmid
            study.pmc_id = pub_meta.pmc_id
            study.authors = pub_meta.authors
            study.publication_date = pub_meta.publication_date
            study.journal = BaseProvenance(
                value=pub_meta.journal,
                source="pubmed",
                source_id=pub_meta.pmid,
                confidence=1.0,
                extraction_method="structured_field",
            )
            study.publication_title = BaseProvenance(
                value=pub_meta.title,
                source="pubmed",
                source_id=pub_meta.pmid,
                confidence=1.0,
                extraction_method="structured_field",
            )
            study.paper_abstract = BaseProvenance(
                value=pub_meta.abstract,
                source="pubmed",
                source_id=pub_meta.pmid,
                confidence=1.0,
                extraction_method="structured_field",
            )
            if pub_meta.doi:
                study.doi = BaseProvenance(
                    value=pub_meta.doi,
                    source="pubmed",
                    source_id=pub_meta.pmid,
                    confidence=1.0,
                    extraction_method="structured_field",
                )
        except ValueError:
            pass  # Publication metadata not available

    study.extraction_date = datetime.now()
    study.extractor_version = "0.1.0"

    return study


def _build_from_biosample(
    sample_id: str,
    bioproject_id: str,
    biosample_meta,
    organism_from_sra: Optional[str] = None,
) -> SampleMetadata:
    """
    Build SampleMetadata from BioSample attributes using enhanced field mapping.

    Uses the new field_mappings approach instead of hardcoded field lists.

    Args:
        sample_id: Sample accession
        bioproject_id: Parent BioProject
        biosample_meta: BioSampleMetadata object
        organism_from_sra: Optional organism from SRA

    Returns:
        SampleMetadata with provenance from BioSample
    """
    # Use enhanced extractor with field_mappings
    return extract_sample_with_field_mappings(
        biosample_meta,
        organism_from_sra,
        sample_id=sample_id,
        bioproject_id=bioproject_id
    )


def build_sample_metadata(
    sample_id: str,
    bioproject_id: str,
    gsm_id: Optional[str] = None,
    biosample_id: Optional[str] = None,
    organism_from_sra: Optional[str] = None,
) -> SampleMetadata:
    """
    Build sample-level metadata from BioSample or GEO.

    Prioritizes BioSample (universal coverage) over GEO (partial coverage).
    Extraction priority: BioSample → GEO → SRA baseline.

    Args:
        sample_id: Sample accession (SRS/SAMN)
        bioproject_id: Parent BioProject
        gsm_id: Optional GEO sample ID
        biosample_id: Optional BioSample ID (preferred source)
        organism_from_sra: Optional organism from SRA metadata

    Returns:
        SampleMetadata with provenance
    """
    # Try to fetch BioSample metadata first (preferred - universal coverage)
    if biosample_id:
        try:
            biosample_meta = fetch_biosample_metadata(biosample_id)
            return _build_from_biosample(sample_id, bioproject_id, biosample_meta, organism_from_sra)
        except ValueError:
            pass  # BioSample not available, fall back to GEO

    # Fall back to GEO if BioSample unavailable
    geo_sample = None
    if gsm_id:
        try:
            geo_sample = fetch_sample_metadata(gsm_id)
        except ValueError:
            pass

    # If neither BioSample nor GEO available, create minimal metadata from SRA
    if not geo_sample:
        return SampleMetadata(
            sample_id=sample_id,
            bioproject_id=bioproject_id,
            organism=BaseProvenance(
                value=organism_from_sra or "unknown",
                source="sra",
                source_id=sample_id,
                confidence=0.8 if organism_from_sra else 0.0,
            ),
            gsm_id=gsm_id,
            biosample_id=biosample_id,
        )

    # Build from GEO sample
    sample = SampleMetadata(
        sample_id=sample_id,
        bioproject_id=bioproject_id,
        gsm_id=geo_sample.gsm_id,
        biosample_id=geo_sample.biosample_id,
        organism=BaseProvenance(
            value=geo_sample.organism or "unknown",
            source="geo",
            source_id=geo_sample.gsm_id,
            confidence=1.0 if geo_sample.organism else 0.0,
            extraction_method="structured_field",
        ),
        sample_title=BaseProvenance(
            value=geo_sample.title,
            source="geo",
            source_id=geo_sample.gsm_id,
            confidence=1.0,
            extraction_method="structured_field",
        ) if geo_sample.title else None,
        sample_description=BaseProvenance(
            value=geo_sample.description,
            source="geo",
            source_id=geo_sample.gsm_id,
            confidence=1.0,
            extraction_method="structured_field",
        ) if geo_sample.description else None,
        raw_characteristics=geo_sample.characteristics,
    )

    # Parse characteristics into structured fields
    # Prioritize tissue extraction as requested
    chars = geo_sample.characteristics

    # PRIORITY 1: Tissue (most important for RiboSeq)
    for key in ["tissue", "tissue type", "source_name", "source"]:
        if key in chars:
            sample.tissue = BaseProvenance(
                value=chars[key],
                source="geo",
                source_id=geo_sample.gsm_id,
                confidence=0.9,
                extracted_text=f"{key}: {chars[key]}",
                extraction_method="structured_field",
            )
            break

    # Cell line
    for key in ["cell line", "cell_line", "cell-line"]:
        if key in chars:
            sample.cell_line = BaseProvenance(
                value=chars[key],
                source="geo",
                source_id=geo_sample.gsm_id,
                confidence=0.9,
                extracted_text=f"{key}: {chars[key]}",
                extraction_method="structured_field",
            )
            break

    # Cell type
    for key in ["cell type", "cell_type", "cell-type"]:
        if key in chars:
            sample.cell_type = BaseProvenance(
                value=chars[key],
                source="geo",
                source_id=geo_sample.gsm_id,
                confidence=0.9,
                extracted_text=f"{key}: {chars[key]}",
                extraction_method="structured_field",
            )
            break

    # Treatment/condition
    for key in ["treatment", "condition", "compound"]:
        if key in chars:
            sample.treatment = BaseProvenance(
                value=chars[key],
                source="geo",
                source_id=geo_sample.gsm_id,
                confidence=0.9,
                extracted_text=f"{key}: {chars[key]}",
                extraction_method="structured_field",
            )
            break

    # Strain
    for key in ["strain", "genetic background"]:
        if key in chars:
            sample.strain = BaseProvenance(
                value=chars[key],
                source="geo",
                source_id=geo_sample.gsm_id,
                confidence=0.9,
                extracted_text=f"{key}: {chars[key]}",
                extraction_method="structured_field",
            )
            break

    # Genotype
    for key in ["genotype", "genetic modification"]:
        if key in chars:
            sample.genotype = BaseProvenance(
                value=chars[key],
                source="geo",
                source_id=geo_sample.gsm_id,
                confidence=0.9,
                extracted_text=f"{key}: {chars[key]}",
                extraction_method="structured_field",
            )
            break

    # Age
    for key in ["age", "developmental stage", "development stage"]:
        if key in chars:
            if "age" in key:
                sample.age = BaseProvenance(
                    value=chars[key],
                    source="geo",
                    source_id=geo_sample.gsm_id,
                    confidence=0.9,
                    extracted_text=f"{key}: {chars[key]}",
                    extraction_method="structured_field",
                )
            else:
                sample.developmental_stage = BaseProvenance(
                    value=chars[key],
                    source="geo",
                    source_id=geo_sample.gsm_id,
                    confidence=0.9,
                    extracted_text=f"{key}: {chars[key]}",
                    extraction_method="structured_field",
                )
            break

    # Sex
    for key in ["sex", "gender"]:
        if key in chars:
            sample.sex = BaseProvenance(
                value=chars[key],
                source="geo",
                source_id=geo_sample.gsm_id,
                confidence=0.9,
                extracted_text=f"{key}: {chars[key]}",
                extraction_method="structured_field",
            )
            break

    return sample


def build_run_metadata(run_id: str) -> RunMetadata:
    """
    Build run-level metadata from SRA.

    Args:
        run_id: Run accession (SRR/ERR/DRR)

    Returns:
        RunMetadata or RiboSeqRunMetadata with provenance
    """
    sra_meta = fetch_run_details(run_id)

    # Detect if this is RiboSeq
    is_riboseq = any(
        term in sra_meta.library_strategy.lower()
        for term in ["ribo", "ribosome profiling"]
    )

    # Base metadata fields
    base_fields = {
        "run_id": sra_meta.run_id,
        "sample_id": sra_meta.sample_id,
        "bioproject_id": sra_meta.bioproject_id,
        "experiment_id": sra_meta.experiment_id,
        "experiment_title": sra_meta.experiment_title,
        "biosample_id": sra_meta.biosample_id,
        "organism": sra_meta.organism,
        "library_strategy": BaseProvenance(
            value=sra_meta.library_strategy,
            source="sra",
            source_id=sra_meta.run_id,
            confidence=1.0,
            extraction_method="structured_field",
        ),
        "library_source": sra_meta.library_source,
        "library_selection": sra_meta.library_selection,
        "library_layout": sra_meta.library_layout,
        "library_name": sra_meta.library_name,
        "platform": sra_meta.platform,
        "instrument_model": sra_meta.instrument_model,
        "read_count": sra_meta.read_count,
        "base_count": sra_meta.base_count,
        "run_date": sra_meta.run_date,
    }

    # Calculate average length if we have counts
    if sra_meta.read_count and sra_meta.base_count and sra_meta.read_count > 0:
        base_fields["avg_length"] = sra_meta.base_count / sra_meta.read_count

    # Return appropriate metadata type
    if is_riboseq:
        return RiboSeqRunMetadata(**base_fields)
    else:
        return RunMetadata(**base_fields)


def build_project_metadata(bioproject_id: str) -> Dict:
    """
    Build complete metadata for a BioProject.

    Aggregates Study + Samples + Runs with full provenance.

    Args:
        bioproject_id: BioProject accession (PRJNA...)

    Returns:
        Dictionary with:
        - study: StudyMetadata
        - samples: Dict[sample_id -> SampleMetadata]
        - runs: List[RunMetadata]
    """
    print(f"Fetching runs for {bioproject_id}...")
    run_ids = fetch_project_runs(bioproject_id)
    print(f"Found {len(run_ids)} runs")

    # Build study metadata
    print("Building study metadata...")
    study = build_study_metadata(bioproject_id, gse_id=None)  # Could auto-detect GSE

    # Collect unique samples
    samples_dict = {}
    runs_list = []

    print("Building run and sample metadata...")
    for i, run_id in enumerate(run_ids, 1):
        print(f"  Processing run {i}/{len(run_ids)}: {run_id}")

        # Build run metadata
        run_meta = build_run_metadata(run_id)
        runs_list.append(run_meta)

        # Build sample metadata if we haven't seen this sample yet
        sample_id = run_meta.sample_id
        if sample_id not in samples_dict:
            sample_meta = build_sample_metadata(
                sample_id=sample_id,
                bioproject_id=bioproject_id,
                gsm_id=None,  # Could auto-detect from SRA
                biosample_id=run_meta.biosample_id,
                organism_from_sra=run_meta.organism,
            )
            # Add technical context for LLM extraction
            sample_meta.technical_context = {
                "experiment_title": run_meta.experiment_title,
                "library_name": run_meta.library_name,
                "library_selection": run_meta.library_selection,
                "library_source": run_meta.library_source,
                "library_layout": run_meta.library_layout,
                "platform": run_meta.platform,
                "instrument_model": run_meta.instrument_model,
                "read_count": run_meta.read_count,
                "run_date": run_meta.run_date,
            }
            samples_dict[sample_id] = sample_meta

    print(f"Complete! Study: 1, Samples: {len(samples_dict)}, Runs: {len(runs_list)}")

    return {
        "study": study,
        "samples": samples_dict,
        "runs": runs_list,
    }
