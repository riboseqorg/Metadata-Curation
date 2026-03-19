"""GEO (Gene Expression Omnibus) metadata fetcher using GEOparse library."""

from typing import Optional, Dict
from pydantic import BaseModel, Field
import GEOparse
import re
from Bio import Entrez
from .entrez_config import configure as _configure_entrez, rate_limit, with_retries

_configure_entrez()


class GEOProjectMetadata(BaseModel):
    """Metadata for a GEO Series (project-level)."""

    gse_id: str = Field(description="GEO Series accession (GSE...)")
    title: str = Field(description="Project title")
    summary: str = Field(description="Project summary/description")
    overall_design: Optional[str] = Field(None, description="Series overall design")
    organism: str = Field(description="Organism(s) studied")

    # Optional fields
    pubmed_id: Optional[str] = Field(None, description="PubMed ID of associated publication")
    bioproject_id: Optional[str] = Field(None, description="Linked BioProject ID")
    contact_name: Optional[str] = Field(None, description="Contact person")
    contact_email: Optional[str] = Field(None, description="Contact email")
    submission_date: Optional[str] = Field(None, description="Submission date")
    last_update_date: Optional[str] = Field(None, description="Last update date")


class GEOSampleMetadata(BaseModel):
    """Metadata for a GEO Sample."""

    gsm_id: str = Field(description="GEO Sample accession (GSM...)")
    gse_id: str = Field(description="Parent GEO Series accession")
    title: str = Field(description="Sample title")
    description: str = Field(description="Sample description")
    characteristics: Dict[str, str] = Field(
        default_factory=dict, description="Sample characteristics (tissue, cell type, treatment, etc)"
    )

    # Optional fields
    source_name: Optional[str] = Field(None, description="Source name")
    organism: Optional[str] = Field(None, description="Organism")
    biosample_id: Optional[str] = Field(None, description="Linked BioSample ID")


def fetch_project_metadata(gse_id: str) -> GEOProjectMetadata:
    """
    Fetch metadata for a GEO Series (project).

    Args:
        gse_id: GEO Series accession (GSE...)

    Returns:
        GEOProjectMetadata object with project information

    Raises:
        ValueError: If GSE ID not found
    """
    try:
        # Use GEOparse to fetch and parse the series
        gse = GEOparse.get_GEO(geo=gse_id, silent=True)

        # Extract metadata from GEO object
        metadata = gse.metadata

        # Get basic info
        title = metadata.get("title", [""])[0]
        summary = metadata.get("summary", [""])[0]
        overall_design = metadata.get("overall_design", [None])[0]
        organism = metadata.get("organism", [""])[0]

        # Get optional fields
        pubmed_id = metadata.get("pubmed_id", [None])[0]
        submission_date = metadata.get("submission_date", [None])[0]
        last_update_date = metadata.get("last_update_date", [None])[0]
        contact_name = metadata.get("contact_name", [None])[0]
        contact_email = metadata.get("contact_email", [None])[0]

        # Extract BioProject from relations
        bioproject_id = None
        relations = metadata.get("relation", [])
        for relation in relations:
            if "BioProject" in relation:
                match = re.search(r"PRJ[NED][A-Z]\d+", relation)
                if match:
                    bioproject_id = match.group(0)
                    break

        return GEOProjectMetadata(
            gse_id=gse_id,
            title=title,
            summary=summary,
            overall_design=overall_design,
            organism=organism,
            pubmed_id=pubmed_id,
            bioproject_id=bioproject_id,
            contact_name=contact_name,
            contact_email=contact_email,
            submission_date=submission_date,
            last_update_date=last_update_date,
        )

    except Exception as e:
        error_msg = str(e).lower()
        if any(phrase in error_msg for phrase in ["does not exist", "not found", "no entries found", "download failed", "no such file"]):
            raise ValueError(f"GEO Series {gse_id} not found")
        raise ValueError(f"Error fetching GEO Series {gse_id}: {e}")


def fetch_sample_metadata(gsm_id: str) -> GEOSampleMetadata:
    """
    Fetch metadata for a GEO Sample.

    Args:
        gsm_id: GEO Sample accession (GSM...)

    Returns:
        GEOSampleMetadata object with sample information

    Raises:
        ValueError: If GSM ID not found
    """
    try:
        # Use GEOparse to fetch and parse the sample
        gsm = GEOparse.get_GEO(geo=gsm_id, silent=True)

        # Extract metadata from GSM object
        metadata = gsm.metadata

        # Get basic info
        title = metadata.get("title", [""])[0]
        description = metadata.get("description", [""])[0]
        source_name = metadata.get("source_name_ch1", [None])[0]
        organism = metadata.get("organism_ch1", [None])[0]

        # Get parent series
        series_id = metadata.get("series_id", ["unknown"])[0]

        # Extract characteristics
        characteristics = {}

        # Extract protocols and molecule
        protocols: Dict[str, str] = {}
        # Known protocol keys per channel
        for ch in ("ch1", "ch2"):
            for kind in ("growth_protocol", "extract_protocol", "treatment_protocol"):
                key = f"{kind}_{ch}"
                if key in metadata:
                    protocols[f"{kind.split('_')[0]}_{ch}"] = " ".join(metadata[key])
        molecule = None
        for ch in ("molecule_ch1", "molecule_ch2"):
            if ch in metadata and metadata[ch]:
                molecule = metadata[ch][0]
                break

        char_keys = [k for k in metadata.keys() if k.startswith("characteristics_ch")]
        for key in char_keys:
            for char in metadata[key]:
                # Parse "key: value" format
                if ":" in char:
                    k, v = char.split(":", 1)
                    characteristics[k.strip()] = v.strip()

        # Extract BioSample ID from relations
        biosample_id = None
        relations = metadata.get("relation", [])
        for relation in relations:
            if "BioSample" in relation:
                match = re.search(r"SAMN\d+", relation)
                if match:
                    biosample_id = match.group(0)
                    break

        return GEOSampleMetadata(
            gsm_id=gsm_id,
            gse_id=series_id,
            title=title,
            description=description,
            characteristics=characteristics,
            protocols=protocols,
            molecule=molecule,
            source_name=source_name,
            organism=organism,
            biosample_id=biosample_id,
        )

    except Exception as e:
        error_msg = str(e).lower()
        if any(phrase in error_msg for phrase in ["does not exist", "not found", "no entries found", "download failed", "no such file"]):
            raise ValueError(f"GEO Sample {gsm_id} not found")
        raise ValueError(f"Error fetching GEO Sample {gsm_id}: {e}")


def search_geo_series(organism: str, terms: Optional[list[str]] = None, boolean_query: Optional[str] = None, retmax: int = 2000) -> list[str]:
    """Search GEO series (GSE) IDs by organism + terms/boolean query via Entrez GDS.

    Returns a list of GSE accessions (e.g., GSE12345).
    """
    q_org = f"{organism}[Organism]"
    q_terms = ""
    if boolean_query:
        q_terms = f" AND ({boolean_query})"
    elif terms:
        ors = " OR ".join(terms)
        q_terms = f" AND ({ors})"
    # Limit to series
    q_full = f"{q_org}{q_terms} AND gse[ETYP]"
    try:
        res = with_retries(lambda: Entrez.esearch(db="gds", term=q_full, retmax=retmax), read=True)
        ids = res.get("IdList", [])
        if not ids:
            return []
        # ESummary to get accessions
        summ = with_retries(lambda: Entrez.esummary(db="gds", id=",".join(ids)), read=True)
        gses = []
        for doc in summ:
            acc = doc.get('Accession') or doc.get('GSE') or ''
            if acc and acc.startswith('GSE'):
                gses.append(acc)
        # unique
        seen = set(); out = []
        for g in gses:
            if g in seen: continue
            seen.add(g); out.append(g)
        return out
    except Exception:
        return []
