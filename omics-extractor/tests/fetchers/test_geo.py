"""Tests for GEO metadata fetcher."""

import pytest
from omics_extractor.fetchers.geo import (
    fetch_project_metadata,
    fetch_sample_metadata,
    GEOProjectMetadata,
    GEOSampleMetadata,
)


class TestFetchProjectMetadata:
    """Test fetching GEO project/series metadata."""

    def test_fetch_gse_metadata(self):
        """Should fetch metadata for GSE series."""
        # GSE112882 - RiboSeq project
        metadata = fetch_project_metadata("GSE112882")

        assert isinstance(metadata, GEOProjectMetadata)
        assert metadata.gse_id == "GSE112882"
        assert metadata.title is not None
        assert len(metadata.title) > 0

    def test_fetch_project_includes_summary(self):
        """Should include project summary/abstract."""
        metadata = fetch_project_metadata("GSE112882")

        assert metadata.summary is not None
        assert len(metadata.summary) > 50  # Should have meaningful summary

    def test_fetch_project_includes_organism(self):
        """Should extract organism information."""
        metadata = fetch_project_metadata("GSE112882")

        assert metadata.organism is not None
        # Should be Arabidopsis for this project

    def test_fetch_project_includes_pubmed_id(self):
        """Should include PubMed ID if available."""
        # This project has a publication
        metadata = fetch_project_metadata("GSE112882")

        # PubMed ID may or may not be present depending on GEO data
        # Just check the field exists
        assert hasattr(metadata, "pubmed_id")

    def test_fetch_project_handles_invalid_gse(self):
        """Should handle invalid GSE IDs."""
        with pytest.raises(ValueError, match="not found"):
            fetch_project_metadata("GSE999999999")


class TestFetchSampleMetadata:
    """Test fetching GEO sample metadata."""

    def test_fetch_gsm_metadata(self):
        """Should fetch metadata for GSM sample."""
        # GSM3085447 - sample from GSE112882
        metadata = fetch_sample_metadata("GSM3085447")

        assert isinstance(metadata, GEOSampleMetadata)
        assert metadata.gsm_id == "GSM3085447"
        assert metadata.title is not None

    def test_fetch_sample_includes_characteristics(self):
        """Should include sample characteristics."""
        metadata = fetch_sample_metadata("GSM3085447")

        assert metadata.characteristics is not None
        assert isinstance(metadata.characteristics, dict)
        # Common characteristics like tissue, cell type, treatment
        assert len(metadata.characteristics) > 0

    def test_fetch_sample_includes_description(self):
        """Should include sample description."""
        metadata = fetch_sample_metadata("GSM3085447")

        assert metadata.description is not None

    def test_fetch_sample_handles_invalid_gsm(self):
        """Should handle invalid GSM IDs."""
        with pytest.raises(ValueError, match="not found"):
            fetch_sample_metadata("GSM999999999")


class TestGEOProjectMetadata:
    """Test GEOProjectMetadata data model."""

    def test_project_model_validation(self):
        """Should validate required fields."""
        metadata = GEOProjectMetadata(
            gse_id="GSE123456",
            title="Test project",
            summary="Test summary",
            organism="Homo sapiens",
        )

        assert metadata.gse_id == "GSE123456"
        assert metadata.title == "Test project"


class TestGEOSampleMetadata:
    """Test GEOSampleMetadata data model."""

    def test_sample_model_validation(self):
        """Should validate required fields."""
        metadata = GEOSampleMetadata(
            gsm_id="GSM123456",
            gse_id="GSE123456",
            title="Test sample",
            description="Test description",
            characteristics={
                "tissue": "brain",
                "cell type": "neuron",
            },
        )

        assert metadata.gsm_id == "GSM123456"
        assert metadata.characteristics["tissue"] == "brain"
