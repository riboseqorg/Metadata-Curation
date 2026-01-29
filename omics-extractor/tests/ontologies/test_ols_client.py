
import pytest
import json
from unittest.mock import patch, MagicMock
from pathlib import Path
from omics_extractor.ontologies.mapper import OlsClient, OntologyMapper

class TestOlsClient:
    
    @pytest.fixture
    def mock_cache_dir(self, tmp_path):
        return tmp_path / "cache"

    def test_search_ols_api_success(self, mock_cache_dir):
        """Should query API and return result when not in cache."""
        client = OlsClient(cache_dir=mock_cache_dir)
        
        # Mock API response
        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_response.json.return_value = {
            "response": {
                "docs": [
                    {
                        "label": "liver",
                        "obo_id": "UBERON:0002107",
                        "ontology_name": "uberon"
                    }
                ]
            }
        }
        
        with patch("requests.get", return_value=mock_response) as mock_get:
            result = client.search("liver", "uberon")
            
            # Verify result
            assert result is not None
            assert result[0] == "liver"
            assert result[1] == "UBERON:0002107"
            assert result[2] == 1.0  # Exact match confidence
            
            # Verify API call
            mock_get.assert_called_once()
            args, kwargs = mock_get.call_args
            assert kwargs["params"]["q"] == "liver"
            assert kwargs["params"]["ontology"] == "uberon"

    def test_search_uses_cache(self, mock_cache_dir):
        """Should return cached result without hitting API."""
        client = OlsClient(cache_dir=mock_cache_dir)
        
        # Pre-populate cache
        cache_key = "uberon:brain"
        client.cache[cache_key] = {
            "value": ("brain", "UBERON:0000955", 1.0),
            "timestamp": 1234567890.0 # Old timestamp but we mocked time.time? No, let's use recent
        }
        
        # Use patch to mock time to ensure cache is valid
        with patch("time.time", return_value=1234567890.0 + 100):
            with patch("requests.get") as mock_get:
                result = client.search("brain", "uberon")
                
                assert result == ("brain", "UBERON:0000955", 1.0)
                mock_get.assert_not_called()

    def test_search_api_failure_returns_none(self, mock_cache_dir):
        """Should return None on API failure."""
        client = OlsClient(cache_dir=mock_cache_dir)
        
        with patch("requests.get", side_effect=Exception("API Error")):
            result = client.search("unknown", "uberon")
            assert result is None

class TestOntologyMapper:
    
    def test_map_tissue_uses_ols(self, tmp_path):
        """OntologyMapper should verify OLS is called."""
        mapper = OntologyMapper(cache_dir=tmp_path)
        
        # Mock OLS client search
        mapper.ols_client.search = MagicMock(return_value=("heart", "UBERON:0000948", 1.0))
        
        result = mapper.map_tissue("heart")
        assert result == ("heart", "UBERON:0000948", 1.0)
        mapper.ols_client.search.assert_called_with("heart", "uberon")

    @patch("omics_extractor.ontologies.mapper.pronto.Ontology")
    def test_map_tissue_uses_local_file(self, mock_ontology_cls, tmp_path):
        """Should fall back to local file if OLS returns None."""
        mapper = OntologyMapper(cache_dir=tmp_path)
        
        # Mock OLS return None (API fail)
        mapper.ols_client.search = MagicMock(return_value=None)
        
        # Mock local file existing
        (tmp_path / "uberon.obo").touch()
        
        # Mock pronto ontology loading
        mock_ontology = MagicMock()
        mock_term = MagicMock()
        mock_term.name = "Kidney"
        mock_term.id = "UBERON:0002113"
        mock_term.synonyms = []
        
        # Mock the terms() iterator
        mock_ontology.terms.return_value = [mock_term]
        mock_ontology_cls.return_value = mock_ontology
        
        # Search for "kidney"
        result = mapper.map_tissue("kidney")
        
        assert result == ("Kidney", "UBERON:0002113", 1.0)
        # Should have tried to load local ontology
        mock_ontology_cls.assert_called()

    def test_map_tissue_fallback(self, tmp_path):
        """Should fall back to legacy map if OLS returns None."""
        mapper = OntologyMapper(cache_dir=tmp_path)
        
        # Mock OLS return None
        mapper.ols_client.search = MagicMock(return_value=None)
        
        # Mock local search return None
        mapper._search_local = MagicMock(return_value=None)
        
        # "liver" is in legacy map
        result = mapper.map_tissue("liver")
        assert result[0] == "liver"
        assert result[1] == "UBERON:0002107"
