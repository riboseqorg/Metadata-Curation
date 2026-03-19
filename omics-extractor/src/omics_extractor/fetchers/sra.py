"""SRA metadata fetcher using NCBI Entrez API."""

from typing import List, Optional
from pydantic import BaseModel, Field
from Bio import Entrez
import xml.etree.ElementTree as ET
import time
from .entrez_config import configure as _configure_entrez, rate_limit

# Configure Entrez once per process
_configure_entrez()




class SRARunMetadata(BaseModel):
    """Metadata for a single SRA sequencing run."""

    run_id: str = Field(description="SRA run accession (SRR/ERR/DRR)")
    experiment_id: str = Field(description="Experiment accession (SRX/ERX/DRX)")
    sample_id: str = Field(description="Sample accession (SRS/ERS/DRS or SAMN)")
    bioproject_id: str = Field(description="BioProject accession (PRJNA/PRJEB/PRJDB)")
    biosample_id: Optional[str] = Field(None, description="BioSample accession (SAMN/SAMEA/SAMD)")
    organism: Optional[str] = Field(None, description="Organism name")

    # Library information
    library_strategy: str = Field(description="Library strategy (RNA-Seq, RiboSeq, etc)")
    library_source: str = Field(description="Library source (TRANSCRIPTOMIC, GENOMIC, etc)")
    library_selection: str = Field(description="Library selection method")
    library_name: Optional[str] = Field(None, description="Library name")
    library_layout: Optional[str] = Field(None, description="SINGLE or PAIRED")

    # Platform information
    platform: str = Field(description="Sequencing platform (ILLUMINA, etc)")
    instrument_model: Optional[str] = Field(None, description="Instrument model")

    # Run statistics
    read_count: Optional[int] = Field(None, description="Number of reads")
    base_count: Optional[int] = Field(None, description="Number of bases")

    # Additional metadata
    experiment_title: Optional[str] = Field(None, description="Experiment title")
    run_date: Optional[str] = Field(None, description="Run publish date")


def fetch_project_runs(project_id: str) -> List[str]:
    """
    Fetch all SRA run accessions for a given BioProject or GEO series.

    Args:
        project_id: BioProject ID (PRJNA...) or GEO series (GSE...)

    Returns:
        List of run accessions (SRR/ERR/DRR IDs)
    """
    # Convert GSE to BioProject if needed
    if project_id.startswith("GSE"):
        project_id = _convert_gse_to_bioproject(project_id)
        if not project_id:
            return []

    # Search SRA database for runs in this project
    try:
        # Use Entrez to search SRA database
        search_term = f"{project_id}[BioProject]"
        handle = Entrez.esearch(
            db="sra", term=search_term, retmax=10000, usehistory="y"
        )
        search_results = Entrez.read(handle)
        handle.close()

        rate_limit()

        if int(search_results["Count"]) == 0:
            return []

        # Fetch the run IDs
        webenv = search_results["WebEnv"]
        query_key = search_results["QueryKey"]

        handle = Entrez.efetch(
            db="sra",
            query_key=query_key,
            WebEnv=webenv,
            rettype="runinfo",
            retmode="text",
        )
        runinfo = handle.read()
        handle.close()

        rate_limit()

        # Decode bytes to string if needed
        if isinstance(runinfo, bytes):
            runinfo = runinfo.decode("utf-8")

        # Parse RunInfo CSV format
        run_ids = []
        for line in runinfo.strip().split("\n")[1:]:  # Skip header
            if line:
                parts = line.split(",")
                run_id = parts[0]  # First column is Run accession
                if run_id and (run_id.startswith("SRR") or run_id.startswith("ERR") or run_id.startswith("DRR")):
                    run_ids.append(run_id)

        return run_ids

    except Exception as e:
        # Handle API errors gracefully
        print(f"Error fetching runs for {project_id}: {e}")
        return []


def fetch_run_details(run_id: str) -> SRARunMetadata:
    """
    Fetch detailed metadata for a single SRA run.

    Args:
        run_id: SRA run accession (SRR/ERR/DRR)

    Returns:
        SRARunMetadata object with detailed information

    Raises:
        ValueError: If run ID not found
    """
    try:
        # Search for the run
        handle = Entrez.esearch(db="sra", term=run_id, retmax=1)
        search_results = Entrez.read(handle)
        handle.close()

        rate_limit()

        if int(search_results["Count"]) == 0:
            raise ValueError(f"Run {run_id} not found in SRA database")

        # Fetch detailed XML
        sra_id = search_results["IdList"][0]
        handle = Entrez.efetch(db="sra", id=sra_id, rettype="full", retmode="xml")
        xml_data = handle.read()
        handle.close()

        rate_limit()

        # Parse XML
        root = ET.fromstring(xml_data)

        # Extract metadata from XML structure
        metadata = _parse_sra_xml(root, run_id)

        # Rate limiting: NCBI allows 3 requests/second without API key
        rate_limit()

        return metadata

    except ValueError:
        raise
    except Exception as e:
        raise ValueError(f"Error fetching run {run_id}: {e}")


def _parse_sra_xml(root: ET.Element, run_id: str) -> SRARunMetadata:
    """Parse SRA XML to extract metadata."""

    # Navigate XML structure
    # SRA XML has structure: EXPERIMENT_PACKAGE_SET/EXPERIMENT_PACKAGE/...
    exp_pkg = root.find(".//EXPERIMENT_PACKAGE")
    if exp_pkg is None:
        raise ValueError("Invalid SRA XML structure")

    # Experiment info
    experiment = exp_pkg.find(".//EXPERIMENT")
    exp_accession = experiment.get("accession") if experiment is not None else "unknown"
    exp_title = experiment.findtext(".//TITLE", default=None)

    # Library info
    library_descriptor = experiment.find(".//LIBRARY_DESCRIPTOR")
    library_strategy = library_descriptor.findtext("LIBRARY_STRATEGY", default="unknown")
    library_source = library_descriptor.findtext("LIBRARY_SOURCE", default="unknown")
    library_selection = library_descriptor.findtext("LIBRARY_SELECTION", default="unknown")
    library_name = library_descriptor.findtext("LIBRARY_NAME", default=None)
    library_layout_elem = library_descriptor.find("LIBRARY_LAYOUT")
    library_layout = "SINGLE" if library_layout_elem.find("SINGLE") is not None else "PAIRED"

    # Platform info
    platform_elem = experiment.find(".//PLATFORM")
    platform = platform_elem.tag if platform_elem is not None else "unknown"
    for child in platform_elem:
        platform = child.tag  # Get the actual platform (ILLUMINA, etc)
        instrument_model = child.findtext("INSTRUMENT_MODEL", default=None)

    # Sample info
    sample = exp_pkg.find(".//SAMPLE")
    sample_accession = sample.get("accession") if sample is not None else "unknown"

    # Extract organism from sample
    organism = None
    if sample is not None:
        organism_elem = sample.find(".//SCIENTIFIC_NAME")
        if organism_elem is not None:
            organism = organism_elem.text

    # Extract BioSample ID from sample identifiers
    biosample_id = None
    if sample is not None:
        for identifier in sample.findall(".//EXTERNAL_ID"):
            if identifier.get("namespace") == "BioSample":
                biosample_id = identifier.text
                break

    # Run info
    run = exp_pkg.find(".//RUN")
    run_accession = run.get("accession") if run is not None else run_id
    run_date = run.get("published") if run is not None else None

    # Run statistics
    read_count = None
    base_count = None
    run_stats = run.find(".//Statistics") if run is not None else None
    if run_stats is not None:
        read_count_str = run_stats.get("nreads")
        base_count_str = run_stats.get("nbases")
        # Handle cases where count is "variable" or non-numeric
        try:
            read_count = int(read_count_str) if read_count_str else None
        except ValueError:
            read_count = None
        try:
            base_count = int(base_count_str) if base_count_str else None
        except ValueError:
            base_count = None

    # BioProject
    external_ids = exp_pkg.findall(".//EXTERNAL_ID")
    bioproject_id = None
    for ext_id in external_ids:
        if ext_id.get("namespace") == "BioProject":
            bioproject_id = ext_id.text
            break

    if not bioproject_id:
        bioproject_id = "unknown"

    return SRARunMetadata(
        run_id=run_accession,
        experiment_id=exp_accession,
        sample_id=sample_accession,
        bioproject_id=bioproject_id,
        biosample_id=biosample_id,
        organism=organism,
        library_strategy=library_strategy,
        library_source=library_source,
        library_selection=library_selection,
        library_name=library_name,
        library_layout=library_layout,
        platform=platform,
        instrument_model=instrument_model,
        read_count=read_count,
        base_count=base_count,
        experiment_title=exp_title,
        run_date=run_date,
    )


def _convert_gse_to_bioproject(gse_id: str) -> Optional[str]:
    """Convert GEO series ID to BioProject ID."""
    try:
        # Search GEO database
        handle = Entrez.esearch(db="gds", term=gse_id, retmax=1)
        search_results = Entrez.read(handle)
        handle.close()

        rate_limit()

        if int(search_results["Count"]) == 0:
            return None

        # Get GEO record
        geo_id = search_results["IdList"][0]
        handle = Entrez.esummary(db="gds", id=geo_id)
        summary = Entrez.read(handle)
        handle.close()

        rate_limit()

        # Extract BioProject from relations
        # This is a simplified approach - may need more robust parsing
        relations = summary[0].get("Relations", [])
        for relation in relations:
            if "BioProject" in relation.get("RelationType", ""):
                return relation.get("TargetObject", "").replace("bioproject/", "PRJNA")

        return None

    except Exception:
        return None
