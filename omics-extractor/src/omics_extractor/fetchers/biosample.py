"""BioSample metadata fetcher using NCBI Entrez API."""

from typing import Optional, Dict
from pydantic import BaseModel, Field
from Bio import Entrez
import xml.etree.ElementTree as ET
from .entrez_config import configure as _configure_entrez, rate_limit

# Configure Entrez once per process
_configure_entrez()




class BioSampleMetadata(BaseModel):
    """Metadata for a BioSample."""

    biosample_id: str = Field(description="BioSample accession (SAMN/SAMEA/SAMD)")
    sample_name: Optional[str] = Field(None, description="Sample name/title")
    organism: Optional[str] = Field(None, description="Organism name")
    attributes: Dict[str, str] = Field(default_factory=dict, description="Sample attributes")
    description: Optional[str] = Field(None, description="Sample description")
    package: Optional[str] = Field(None, description="BioSample package")
    bioproject_id: Optional[str] = Field(None, description="Parent BioProject")


def fetch_biosample_metadata(biosample_id: str) -> BioSampleMetadata:
    """
    Fetch metadata for a BioSample.

    BioSample is NCBI's database for biological sample metadata. Each BioSample
    record contains attributes like tissue, cell type, strain, treatment, etc.

    Args:
        biosample_id: BioSample accession (SAMN/SAMEA/SAMD)

    Returns:
        BioSampleMetadata object with sample attributes

    Raises:
        ValueError: If BioSample ID not found
    """
    try:
        # Search for the BioSample
        handle = Entrez.esearch(db="biosample", term=biosample_id, retmax=1)
        search_results = Entrez.read(handle)
        handle.close()

        rate_limit()

        if int(search_results["Count"]) == 0:
            raise ValueError(f"BioSample {biosample_id} not found")

        # Fetch detailed XML
        biosample_uid = search_results["IdList"][0]
        handle = Entrez.efetch(db="biosample", id=biosample_uid, rettype="full", retmode="xml")
        xml_data = handle.read()
        handle.close()

        rate_limit()

        # Parse XML
        root = ET.fromstring(xml_data)

        # Extract metadata from XML structure
        metadata = _parse_biosample_xml(root, biosample_id)

        return metadata

    except ValueError:
        raise
    except Exception as e:
        raise ValueError(f"Error fetching BioSample {biosample_id}: {e}")


def _parse_biosample_xml(root: ET.Element, biosample_id: str) -> BioSampleMetadata:
    """Parse BioSample XML to extract metadata."""

    # BioSample XML structure: BioSampleSet/BioSample
    biosample = root.find(".//BioSample")
    if biosample is None:
        raise ValueError("Invalid BioSample XML structure")

    # Extract description
    description_elem = biosample.find(".//Description")
    sample_name = None
    organism = None
    description = None

    if description_elem is not None:
        # Title
        title_elem = description_elem.find(".//Title")
        if title_elem is not None:
            sample_name = title_elem.text

        # Organism
        organism_elem = description_elem.find(".//Organism")
        if organism_elem is not None:
            organism_name = organism_elem.findtext("OrganismName")
            if organism_name:
                organism = organism_name

        # Comment/Description
        comment_elem = description_elem.find(".//Comment")
        if comment_elem is not None and comment_elem.find(".//Paragraph") is not None:
            description = comment_elem.findtext(".//Paragraph")

    # Extract package
    package = None
    pkg_elem = biosample.find(".//Package")
    if pkg_elem is not None and pkg_elem.text:
        package = pkg_elem.text

    # Extract attributes
    attributes = {}
    attributes_elem = biosample.find(".//Attributes")
    if attributes_elem is not None:
        for attr in attributes_elem.findall("Attribute"):
            attr_name = attr.get("attribute_name")
            attr_value = attr.text
            if attr_name and attr_value:
                # Normalize attribute names to lowercase with underscores
                normalized_name = attr_name.lower().replace(" ", "_")
                attributes[normalized_name] = attr_value

    # Extract BioProject links
    bioproject_id = None
    links = biosample.find(".//Links")
    if links is not None:
        for link in links.findall("Link"):
            if link.get("type") == "url" and link.get("label") == "BioProject":
                # Extract PRJNA from URL
                url = link.get("target", "")
                if "bioproject" in url.lower():
                    parts = url.split("/")
                    for part in parts:
                        if part.startswith("PRJ"):
                            bioproject_id = part
                            break

    return BioSampleMetadata(
        biosample_id=biosample_id,
        sample_name=sample_name,
        organism=organism,
        attributes=attributes,
        description=description,
        package=package,
        bioproject_id=bioproject_id,
    )
