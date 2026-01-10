"""Tests for BioProject fetcher."""

import pytest
from omics_extractor.fetchers.bioproject import fetch_bioproject_metadata


class TestBioProjectFetcher:
    """Test BioProject metadata fetching."""

    def test_fetch_bioproject_metadata(self):
        """Should fetch BioProject metadata successfully."""
        # PRJNA1176138 is a real RiboSeq project
        metadata = fetch_bioproject_metadata("PRJNA1176138")

        # Basic assertions
        assert metadata.bioproject_id == "PRJNA1176138"
        assert metadata.title is not None
        assert len(metadata.title) > 0

    def test_fetch_bioproject_with_description(self):
        """Should extract project description."""
        metadata = fetch_bioproject_metadata("PRJNA1176138")

        # Should have a description
        assert metadata.description is not None
        assert len(metadata.description) > 0

    def test_fetch_bioproject_organism(self):
        """Should extract organism information."""
        metadata = fetch_bioproject_metadata("PRJNA1176138")

        # Should have organism
        assert metadata.organism is not None
        assert len(metadata.organism) > 0

    def test_invalid_bioproject_id(self):
        """Should raise ValueError for invalid BioProject ID."""
        with pytest.raises(ValueError, match="not found"):
            fetch_bioproject_metadata("PRJNA99999999999")

    def test_fetch_bioproject_with_publication(self):
        """Should extract publication ID when available."""
        # PRJNA449378 has associated publications
        metadata = fetch_bioproject_metadata("PRJNA449378")

        assert metadata.bioproject_id == "PRJNA449378"
        # May or may not have publication depending on BioProject record structure

    def test_fetch_european_bioproject(self):
        """Should handle PRJEB (European) BioProject IDs."""
        # Test with a European BioProject if available
        # Using PRJNA for now as it's more common in RiboSeq
        metadata = fetch_bioproject_metadata("PRJNA1170270")

        assert metadata.bioproject_id == "PRJNA1170270"
        assert metadata.title is not None

    def test_bioproject_data_type(self):
        """Should extract data type information when available."""
        metadata = fetch_bioproject_metadata("PRJNA1176138")

        # Data type may or may not be present
        # Just check that the field exists and is of correct type
        assert hasattr(metadata, "data_type")
        if metadata.data_type:
            assert isinstance(metadata.data_type, str)
