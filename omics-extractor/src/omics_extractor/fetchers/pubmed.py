"""PubMed/PMC metadata fetcher using NCBI Entrez API."""

from typing import Optional, List
from pydantic import BaseModel, Field
from Bio import Entrez
import xml.etree.ElementTree as ET
from .entrez_config import configure as _configure_entrez, rate_limit

# Configure Entrez once per process
_configure_entrez()




class PubMedMetadata(BaseModel):
    """Metadata for a PubMed publication."""

    pmid: str = Field(description="PubMed ID")
    title: str = Field(description="Article title")
    abstract: str = Field(description="Article abstract")
    authors: List[str] = Field(description="List of authors")
    journal: str = Field(description="Journal name")
    publication_date: str = Field(description="Publication date")

    # Optional fields
    doi: Optional[str] = Field(None, description="Digital Object Identifier")
    pmc_id: Optional[str] = Field(None, description="PubMed Central ID (PMC...)")
    pmc_url: Optional[str] = Field(None, description="URL to full text in PMC")
    keywords: Optional[List[str]] = Field(None, description="MeSH keywords")
    affiliations: Optional[List[str]] = Field(None, description="Author affiliations")


def fetch_publication_metadata(pmid: str) -> PubMedMetadata:
    """
    Fetch metadata for a PubMed publication.

    Args:
        pmid: PubMed ID

    Returns:
        PubMedMetadata object with publication information

    Raises:
        ValueError: If PMID not found
    """
    try:
        # Fetch from PubMed database
        handle = Entrez.efetch(db="pubmed", id=pmid, rettype="xml", retmode="xml")
        xml_data = handle.read()
        handle.close()

        rate_limit()

        # Parse XML
        root = ET.fromstring(xml_data)

        # Check if article found
        article = root.find(".//PubmedArticle")
        if article is None:
            raise ValueError(f"PubMed ID {pmid} not found")

        # Extract basic info
        title = article.findtext(".//ArticleTitle", default="")
        abstract_text = article.findtext(".//Abstract/AbstractText", default="")

        # Extract authors
        authors = []
        author_list = article.findall(".//Author")
        for author in author_list:
            last_name = author.findtext("LastName", default="")
            initials = author.findtext("Initials", default="")
            if last_name and initials:
                authors.append(f"{last_name} {initials}")
            elif last_name:
                authors.append(last_name)

        # Extract journal info
        journal = article.findtext(".//Journal/Title", default="")
        pub_date_year = article.findtext(".//PubDate/Year", default="")
        pub_date_month = article.findtext(".//PubDate/Month", default="")
        pub_date_day = article.findtext(".//PubDate/Day", default="")

        # Format publication date
        pub_date_parts = [pub_date_year, pub_date_month, pub_date_day]
        publication_date = "-".join(p for p in pub_date_parts if p)

        # Extract DOI
        doi = None
        article_ids = article.findall(".//ArticleId")
        for article_id in article_ids:
            if article_id.get("IdType") == "doi":
                doi = article_id.text
                break

        # Extract PMC ID
        pmc_id = None
        pmc_url = None
        for article_id in article_ids:
            if article_id.get("IdType") == "pmc":
                pmc_id = article_id.text
                if pmc_id:
                    pmc_url = f"https://www.ncbi.nlm.nih.gov/pmc/articles/{pmc_id}/"
                break

        # Extract keywords/MeSH terms
        keywords = []
        mesh_headings = article.findall(".//MeshHeading/DescriptorName")
        for mesh in mesh_headings:
            keywords.append(mesh.text)

        # Extract affiliations
        affiliations = []
        affiliation_list = article.findall(".//Affiliation")
        for affil in affiliation_list:
            if affil.text:
                affiliations.append(affil.text)

        return PubMedMetadata(
            pmid=pmid,
            title=title,
            abstract=abstract_text,
            authors=authors,
            journal=journal,
            publication_date=publication_date,
            doi=doi,
            pmc_id=pmc_id,
            pmc_url=pmc_url,
            keywords=keywords if keywords else None,
            affiliations=affiliations if affiliations else None,
        )

    except ValueError:
        raise
    except Exception as e:
        raise ValueError(f"Error fetching PubMed ID {pmid}: {e}")


def search_pubmed_for_project(project_id: str) -> List[str]:
    """
    Search PubMed for publications linked to a BioProject or GEO series.

    Args:
        project_id: BioProject ID (PRJNA...) or GEO series (GSE...)

    Returns:
        List of PubMed IDs (may be empty)
    """
    try:
        # Search PubMed for the project ID
        search_term = f"{project_id}[All Fields]"
        handle = Entrez.esearch(db="pubmed", term=search_term, retmax=100)
        search_results = Entrez.read(handle)
        handle.close()

        rate_limit()

        pmids = search_results.get("IdList", [])
        return pmids

    except Exception:
        # If search fails, return empty list (not an error)
        return []


def fetch_full_text_pmc(pmc_id: str) -> Optional[str]:
    """
    Fetch full text XML from PubMed Central.

    This is a placeholder for future implementation where we extract
    methods/materials sections from PMC articles.

    Args:
        pmc_id: PMC ID (e.g., "PMC1234567")

    Returns:
        Full text XML string, or None if not available

    Note:
        This requires PMC Open Access subset access.
        Not all papers are available as full text.
    """
    try:
        # Remove "PMC" prefix if present
        if pmc_id.startswith("PMC"):
            pmc_id = pmc_id[3:]

        # Fetch from PMC
        handle = Entrez.efetch(db="pmc", id=pmc_id, rettype="xml", retmode="xml")
        xml_data = handle.read()
        handle.close()

        rate_limit()

        # Decode bytes to string if needed
        if isinstance(xml_data, bytes):
            xml_data = xml_data.decode("utf-8")

        return xml_data

    except Exception:
        # PMC full text may not be available
        return None
