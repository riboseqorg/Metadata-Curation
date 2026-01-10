"""BioProject metadata fetcher using NCBI Entrez API."""

from typing import Optional
from pydantic import BaseModel, Field
from Bio import Entrez
import xml.etree.ElementTree as ET


# Set your email for NCBI Entrez (required by NCBI)
Entrez.email = "jackcurragh@gmail.com"


class BioProjectMetadata(BaseModel):
    """Metadata for a BioProject."""

    bioproject_id: str = Field(description="BioProject accession (PRJNA/PRJEB/PRJDB)")
    title: Optional[str] = Field(None, description="Project title")
    description: Optional[str] = Field(None, description="Project description")
    organism: Optional[str] = Field(None, description="Target organism(s)")
    publication_id: Optional[str] = Field(None, description="Associated PubMed ID")
    data_type: Optional[str] = Field(None, description="Data type (e.g., Transcriptome)")
    scope: Optional[str] = Field(None, description="Project scope (e.g., Multiisolate)")


def fetch_bioproject_metadata(bioproject_id: str) -> BioProjectMetadata:
    """
    Fetch metadata for a BioProject.

    BioProject is NCBI's organizing principle for large-scale genomic data.
    Each BioProject contains project-level metadata including title, description,
    organism, and publication information.

    Args:
        bioproject_id: BioProject accession (PRJNA/PRJEB/PRJDB)

    Returns:
        BioProjectMetadata object with project information

    Raises:
        ValueError: If BioProject ID not found
    """
    try:
        # Search for the BioProject
        handle = Entrez.esearch(db="bioproject", term=bioproject_id, retmax=1)
        search_results = Entrez.read(handle)
        handle.close()

        if int(search_results["Count"]) == 0:
            raise ValueError(f"BioProject {bioproject_id} not found")

        # Fetch detailed XML
        bioproject_uid = search_results["IdList"][0]
        handle = Entrez.efetch(db="bioproject", id=bioproject_uid, rettype="xml")
        xml_data = handle.read()
        handle.close()

        # Parse XML
        root = ET.fromstring(xml_data)

        # Extract metadata from XML structure
        metadata = _parse_bioproject_xml(root, bioproject_id)

        return metadata

    except ValueError:
        raise
    except Exception as e:
        raise ValueError(f"Error fetching BioProject {bioproject_id}: {e}")


def _parse_bioproject_xml(root: ET.Element, bioproject_id: str) -> BioProjectMetadata:
    """Parse BioProject XML to extract metadata."""

    # BioProject XML structure can be:
    # 1. RecordSet/DocumentSummary/Project (from efetch)
    # 2. Package/Project (from efetch with different params)
    # 3. DocumentSummarySet/DocumentSummary (from esummary)

    # Try to find Project element (most common from efetch)
    project = root.find(".//Project")
    if project is not None:
        return _parse_from_project(project, bioproject_id)

    # Fallback to DocumentSummary without Project (from esummary)
    doc_summary = root.find(".//DocumentSummary")
    if doc_summary is not None and doc_summary.find(".//Project") is None:
        return _parse_from_summary(doc_summary, bioproject_id)

    raise ValueError("Invalid BioProject XML structure")


def _parse_from_summary(doc_summary: ET.Element, bioproject_id: str) -> BioProjectMetadata:
    """Parse from DocumentSummary format."""

    title = doc_summary.findtext("Project_Title")
    description = doc_summary.findtext("Project_Description")

    # Get organism
    organism = None
    organism_elem = doc_summary.find(".//Organism")
    if organism_elem is not None:
        organism = organism_elem.findtext("OrganismName")

    # Get publication
    publication_id = None
    pub_elem = doc_summary.find(".//Publication")
    if pub_elem is not None:
        publication_id = pub_elem.get("id")

    return BioProjectMetadata(
        bioproject_id=bioproject_id,
        title=title,
        description=description,
        organism=organism,
        publication_id=publication_id,
    )


def _parse_from_project(project: ET.Element, bioproject_id: str) -> BioProjectMetadata:
    """Parse from Project element format."""

    # Get project descriptor
    descriptor = project.find(".//ProjectDescr")

    title = None
    description = None
    publication_id = None

    if descriptor is not None:
        # Title (can be in Title or Name element)
        title_elem = descriptor.find(".//Title")
        if title_elem is not None and title_elem.text:
            title = title_elem.text
        else:
            # Try Name element
            name_elem = descriptor.find(".//Name")
            if name_elem is not None and name_elem.text:
                title = name_elem.text

        # Description
        desc_elem = descriptor.find(".//Description")
        if desc_elem is not None and desc_elem.text:
            description = desc_elem.text

        # Publication - check for Publication element with DbType
        for pub_elem in descriptor.findall(".//Publication"):
            # Try to find PubMed ID
            db_elem = pub_elem.find(".//DbType[@db='pubmed']")
            if db_elem is not None and db_elem.text:
                publication_id = db_elem.text
                break
            # Also try id attribute directly on Publication
            if pub_elem.get("id"):
                publication_id = pub_elem.get("id")
                break

    # Get organism from ProjectType
    organism = None
    project_type = project.find(".//ProjectType")
    if project_type is not None:
        organism_elem = project_type.find(".//Organism")
        if organism_elem is not None:
            org_name = organism_elem.findtext(".//OrganismName")
            if org_name:
                organism = org_name

    # Get data type and scope
    data_type = None
    scope = None
    if project_type is not None:
        proj_type_submission = project_type.find(".//ProjectTypeSubmission")
        if proj_type_submission is not None:
            target = proj_type_submission.find(".//Target")
            if target is not None:
                data_type = target.get("sample_scope")

            method = proj_type_submission.find(".//Method")
            if method is not None:
                scope = method.get("method_type")

    return BioProjectMetadata(
        bioproject_id=bioproject_id,
        title=title,
        description=description,
        organism=organism,
        publication_id=publication_id,
        data_type=data_type,
        scope=scope,
    )
