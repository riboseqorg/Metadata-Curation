"""Tests for ontology mapping."""

import pytest
from pathlib import Path
from omics_extractor.ontologies.mapper import (
    OntologyMapper,
    normalize_tissue,
    normalize_cell_type,
    normalize_organism,
)


class TestOrganismMapping:
    """Test organism normalization (doesn't require ontology download)."""

    def test_exact_match_scientific_name(self):
        """Should match exact scientific names with highest confidence."""
        value, taxon_id, conf = normalize_organism("Mus musculus")
        assert value == "Mus musculus"
        assert taxon_id == "NCBITaxon:10090"
        assert conf == 1.0

    def test_common_name_match(self):
        """Should match common names with slightly lower confidence."""
        value, taxon_id, conf = normalize_organism("mouse")
        assert value == "Mus musculus"
        assert taxon_id == "NCBITaxon:10090"
        assert conf == 0.95

    def test_case_insensitive(self):
        """Should handle different capitalizations."""
        value, taxon_id, conf = normalize_organism("HOMO SAPIENS")
        assert value == "Homo sapiens"
        assert taxon_id == "NCBITaxon:9606"
        assert conf == 1.0

    def test_unknown_organism(self):
        """Should return original value for unknown organisms."""
        value, taxon_id, conf = normalize_organism("Unknown species")
        assert value == "Unknown species"
        assert taxon_id is None
        assert conf == 0.5

    def test_model_organisms(self):
        """Should recognize all common model organisms."""
        test_cases = [
            ("human", "Homo sapiens", "NCBITaxon:9606"),
            ("rat", "Rattus norvegicus", "NCBITaxon:10116"),
            ("zebrafish", "Danio rerio", "NCBITaxon:7955"),
            ("c. elegans", "Caenorhabditis elegans", "NCBITaxon:6239"),
            ("yeast", "Saccharomyces cerevisiae", "NCBITaxon:4932"),
        ]

        for raw, expected_name, expected_id in test_cases:
            value, taxon_id, conf = normalize_organism(raw)
            assert value == expected_name
            assert taxon_id == expected_id
            assert conf >= 0.9


@pytest.mark.integration
class TestTissueMapping:
    """
    Test tissue mapping with UBERON ontology.

    These tests require downloading the UBERON ontology (~50MB).
    Marked as integration tests to allow skipping in CI.
    """

    def test_exact_label_match(self):
        """Should match exact ontology labels."""
        # This will download UBERON on first run
        value, onto_id, conf = normalize_tissue("brain")

        # Should find UBERON term for brain
        assert onto_id is not None
        assert "UBERON" in onto_id
        assert conf >= 0.95  # High confidence for exact match

    def test_case_insensitive_match(self):
        """Should match regardless of case."""
        value1, id1, conf1 = normalize_tissue("Brain")
        value2, id2, conf2 = normalize_tissue("brain")
        value3, id3, conf3 = normalize_tissue("BRAIN")

        # All should map to same term
        assert id1 == id2 == id3
        # Case-insensitive gets slightly lower confidence
        assert conf1 >= 0.95

    def test_unknown_tissue(self):
        """Should return original value for unmapped terms."""
        value, onto_id, conf = normalize_tissue("fake_tissue_xyz")
        assert value == "fake_tissue_xyz"
        assert onto_id is None
        assert conf == 0.5

    def test_synonym_match(self):
        """Should match ontology synonyms."""
        # Try a tissue with known synonyms
        value, onto_id, conf = normalize_tissue("liver")

        assert onto_id is not None
        assert conf >= 0.9


@pytest.mark.integration
class TestCellTypeMapping:
    """
    Test cell type mapping with Cell Ontology.

    Requires downloading Cell Ontology.
    """

    def test_exact_cell_type_match(self):
        """Should match exact cell type labels."""
        value, onto_id, conf = normalize_cell_type("neuron")

        assert onto_id is not None
        assert "CL" in onto_id
        assert conf >= 0.95

    def test_unknown_cell_type(self):
        """Should return original for unknown cell types."""
        value, onto_id, conf = normalize_cell_type("fake_cell_type")
        assert value == "fake_cell_type"
        assert onto_id is None
        assert conf == 0.5


class TestOntologyMapper:
    """Test OntologyMapper class directly."""

    def test_mapper_singleton(self):
        """Should create singleton mapper."""
        from omics_extractor.ontologies.mapper import get_mapper

        mapper1 = get_mapper()
        mapper2 = get_mapper()

        assert mapper1 is mapper2

    def test_cache_directory_created(self):
        """Should create cache directory."""
        mapper = OntologyMapper()
        assert mapper.cache_dir.exists()
        assert mapper.cache_dir.is_dir()

    def test_custom_cache_dir(self):
        """Should use custom cache directory."""
        custom_dir = Path("/tmp/test_ontology_cache")
        mapper = OntologyMapper(cache_dir=custom_dir)

        assert mapper.cache_dir == custom_dir
        assert mapper.cache_dir.exists()
