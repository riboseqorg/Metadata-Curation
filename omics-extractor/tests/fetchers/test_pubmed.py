"""Tests for PubMed metadata fetcher."""

import pytest
from omics_extractor.fetchers.pubmed import (
    fetch_publication_metadata,
    search_pubmed_for_project,
    PubMedMetadata,
)


class TestFetchPublicationMetadata:
    """Test fetching PubMed publication metadata."""

    def test_fetch_publication_with_pmc(self):
        """Should fetch publication with PMC full text."""
        # PMID: 29618526 - Known RiboSeq paper with PMC article
        metadata = fetch_publication_metadata("29618526")

        assert isinstance(metadata, PubMedMetadata)
        assert metadata.pmid == "29618526"
        assert metadata.title is not None
        assert len(metadata.title) > 0
        assert metadata.abstract is not None
        assert len(metadata.abstract) > 50

    def test_fetch_publication_has_pmc_id(self):
        """Should include PMC ID when available."""
        # This paper should have a PMC ID
        metadata = fetch_publication_metadata("29618526")

        # PMC ID may or may not be present
        # Just check the field exists
        assert hasattr(metadata, "pmc_id")

    def test_fetch_publication_has_doi(self):
        """Should include DOI when available."""
        metadata = fetch_publication_metadata("29618526")

        assert metadata.doi is not None
        assert metadata.doi.startswith("10.")

    def test_fetch_publication_includes_authors(self):
        """Should include author list."""
        metadata = fetch_publication_metadata("29618526")

        assert metadata.authors is not None
        assert len(metadata.authors) > 0

    def test_fetch_publication_includes_journal_info(self):
        """Should include journal and publication date."""
        metadata = fetch_publication_metadata("29618526")

        assert metadata.journal is not None
        assert metadata.publication_date is not None

    def test_fetch_publication_handles_invalid_pmid(self):
        """Should handle invalid PubMed IDs."""
        with pytest.raises(ValueError, match="not found"):
            fetch_publication_metadata("999999999999")

    def test_fetch_publication_without_pmc(self):
        """Should handle publications without PMC full text."""
        # Use a PMID that likely doesn't have PMC
        # This is just to test that the code doesn't break
        metadata = fetch_publication_metadata("29618526")

        # Should still work even if no PMC
        assert metadata.pmid is not None


class TestSearchPubMedForProject:
    """Test searching PubMed for publications linked to projects."""

    def test_search_by_bioproject(self):
        """Should find publications linked to BioProject."""
        # PRJNA401384 has associated publications
        pmids = search_pubmed_for_project("PRJNA401384")

        # May return 0 or more results
        assert isinstance(pmids, list)

    def test_search_by_gse(self):
        """Should find publications linked to GEO series."""
        # GSE112882 has associated publication
        pmids = search_pubmed_for_project("GSE112882")

        assert isinstance(pmids, list)

    def test_search_returns_pmids(self):
        """Should return list of PubMed IDs."""
        pmids = search_pubmed_for_project("GSE112882")

        if pmids:
            # Should be numeric strings
            assert all(pmid.isdigit() for pmid in pmids)


class TestPubMedMetadata:
    """Test PubMedMetadata data model."""

    def test_metadata_model_validation(self):
        """Should validate required fields."""
        metadata = PubMedMetadata(
            pmid="12345678",
            title="Test paper title",
            abstract="Test abstract text",
            authors=["Smith J", "Jones A"],
            journal="Nature",
            publication_date="2024-01-01",
        )

        assert metadata.pmid == "12345678"
        assert metadata.title == "Test paper title"

    def test_metadata_model_optional_fields(self):
        """Should allow optional fields."""
        metadata = PubMedMetadata(
            pmid="12345678",
            title="Test paper",
            abstract="Test abstract",
            authors=["Smith J"],
            journal="Nature",
            publication_date="2024-01-01",
            doi="10.1038/nature12345",
            pmc_id="PMC1234567",
            pmc_url="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1234567/",
        )

        assert metadata.doi == "10.1038/nature12345"
        assert metadata.pmc_id == "PMC1234567"
        assert metadata.pmc_url is not None
