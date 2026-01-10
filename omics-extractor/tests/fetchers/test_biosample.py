"""Tests for BioSample fetcher."""

import pytest
from omics_extractor.fetchers.biosample import fetch_biosample_metadata


class TestBioSampleFetcher:
    """Test BioSample metadata fetching."""

    def test_fetch_biosample_metadata(self):
        """Should fetch BioSample metadata successfully."""
        # SAMN44381105 is from PRJNA1176138 (RiboSeq project)
        metadata = fetch_biosample_metadata("SAMN44381105")

        # Basic assertions
        assert metadata.biosample_id == "SAMN44381105"
        assert metadata.organism is not None
        assert len(metadata.attributes) > 0

        # Should have sample attributes
        assert isinstance(metadata.attributes, dict)

    def test_fetch_biosample_with_tissue(self):
        """Should extract tissue information from BioSample attributes."""
        # SAMN44381105 is from PRJNA1176138, should have biological metadata
        metadata = fetch_biosample_metadata("SAMN44381105")

        # Check if we got tissue or related attributes
        attr_keys = [k.lower() for k in metadata.attributes.keys()]
        has_biological_attrs = any(
            keyword in " ".join(attr_keys)
            for keyword in ["tissue", "cell", "strain", "genotype", "treatment"]
        )
        assert has_biological_attrs, f"Expected biological attributes, got: {list(metadata.attributes.keys())}"

    def test_invalid_biosample_id(self):
        """Should raise ValueError for invalid BioSample ID."""
        with pytest.raises(ValueError, match="not found"):
            fetch_biosample_metadata("SAMN99999999999")

    def test_biosample_organism(self):
        """Should extract organism name."""
        metadata = fetch_biosample_metadata("SAMN44381105")
        assert metadata.organism is not None
        assert len(metadata.organism) > 0

    def test_biosample_attributes_dict(self):
        """Should return attributes as a dictionary."""
        metadata = fetch_biosample_metadata("SAMN44381105")
        assert isinstance(metadata.attributes, dict)
        # Should have some attributes
        assert len(metadata.attributes) > 0
        # All keys and values should be strings
        for key, value in metadata.attributes.items():
            assert isinstance(key, str)
            assert isinstance(value, str)
