"""Ontology-based metadata normalization and mapping.

This module provides normalization confidence scoring based on how well
raw metadata values match standardized ontology terms.

Normalization Confidence Levels:
- 1.0: Exact match to ontology term ID or label
- 0.95: Case-insensitive match to label
- 0.9: Exact match to synonym
- 0.7: Fuzzy match to synonym (edit distance)
- 0.5: No ontology match found (keeps raw value)

Supported Ontologies:
- UBERON: Tissue/anatomy terms
- CL: Cell types
- NCBITaxon: Organism taxonomy
- EFO: Experimental Factor Ontology (treatment, disease, etc.)
"""

from typing import Optional, Tuple, Dict, Any, List
from pathlib import Path
import json
import time
import requests
import os
import pronto
from difflib import SequenceMatcher


class OlsClient:
    """Client for EBI Ontology Lookup Service (OLS)."""

    BASE_URL = "https://www.ebi.ac.uk/ols/api/search"

    def __init__(self, cache_dir: Optional[Path] = None):
        self.cache_dir = cache_dir or Path.home() / ".omics_extractor" / "ontologies"
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.cache_file = self.cache_dir / "ols_cache.json"
        self._load_cache()

    def _load_cache(self):
        """Load cache from disk."""
        if self.cache_file.exists():
            try:
                self.cache = json.loads(self.cache_file.read_text())
            except json.JSONDecodeError:
                self.cache = {}
        else:
            self.cache = {}

    def _save_cache(self):
        """Save cache to disk."""
        try:
            self.cache_file.write_text(json.dumps(self.cache, indent=2))
        except Exception:
            pass  # Ignore cache write errors

    def search(self, query: str, ontology: str) -> Optional[Tuple[str, str, float]]:
        """
        Search OLS for a term.
        
        Returns:
            Tuple(label, obo_id, score_confidence) or None
        """
        cache_key = f"{ontology}:{query.lower().strip()}"
        
        # Check cache (valid for 30 days)
        if cache_key in self.cache:
            cached_item = self.cache[cache_key]
            if time.time() - cached_item["timestamp"] < 30 * 24 * 3600:
                return cached_item["value"]
        
        # Query API
        try:
            params = {
                "q": query,
                "ontology": ontology.lower(),
                "rows": 1,
                "type": "class",
                "exact": "false"
            }
            # Timeout is short to quickly failover to local
            response = requests.get(self.BASE_URL, params=params, timeout=3)
            
            if response.status_code == 200:
                data = response.json()
                docs = data.get("response", {}).get("docs", [])
                
                if docs:
                    best_match = docs[0]
                    label = best_match.get("label")
                    obo_id = best_match.get("obo_id")
                    
                    # Calculate simple confidence
                    sim_score = SequenceMatcher(None, query.lower(), label.lower()).ratio()
                    
                    if sim_score == 1.0:
                        conf = 1.0
                    elif query.lower() == label.lower():
                        conf = 1.0
                    elif sim_score > 0.9:
                        conf = 0.95
                    elif sim_score > 0.8:
                        conf = 0.9
                    else:
                        conf = 0.7
                        
                    result = (label, obo_id, conf)
                    
                    # Cache result
                    self.cache[cache_key] = {
                        "value": result,
                        "timestamp": time.time()
                    }
                    self._save_cache()
                    
                    return result

        except Exception as e:
            # Silently fail to fallback
            pass
            
        return None


class OntologyMapper:
    """Maps raw metadata values to ontology terms with confidence scoring."""

    def __init__(self, cache_dir: Optional[Path] = None):
        """
        Initialize ontology mapper with API client and local fallback.
        """
        self.cache_dir = cache_dir or Path.home() / ".omics_extractor" / "ontologies"
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        
        self.ols_client = OlsClient(self.cache_dir)
        self._ontologies: Dict[str, pronto.Ontology] = {}
        
        # Define fallback OBO URLs
        self._obo_urls = {
            "uberon": "http://purl.obolibrary.org/obo/uberon.obo",
            "cl": "http://purl.obolibrary.org/obo/cl.obo",
            # ncbitaxon is too huge for local fallback usually, but we'll include logic
        }

        # Legacy/Fallback mappings (Hardcoded lists for offline safety)
        self._legacy_tissue_map = {
            "liver": ("liver", "UBERON:0002107", 1.0),
            "brain": ("brain", "UBERON:0000955", 1.0),
            "kidney": ("kidney", "UBERON:0002113", 1.0),
            "heart": ("heart", "UBERON:0000948", 1.0),
            "lung": ("lung", "UBERON:0002048", 1.0),
            "muscle": ("muscle tissue", "UBERON:0002385", 0.95),
            "blood": ("blood", "UBERON:0000178", 1.0),
            "spleen": ("spleen", "UBERON:0002106", 1.0),
            "skin": ("skin", "UBERON:0001003", 1.0),
        }

    def _get_local_ontology(self, name: str) -> Optional[pronto.Ontology]:
        """Lazy load local ontology if available."""
        if name in self._ontologies:
            return self._ontologies[name]
            
        # Check if file exists locally
        obo_path = self.cache_dir / f"{name}.obo"
        if obo_path.exists():
            try:
                print(f"Loading local ontology: {name} (this may take a moment)...")
                self._ontologies[name] = pronto.Ontology(str(obo_path))
                return self._ontologies[name]
            except Exception as e:
                print(f"Failed to load local ontology {name}: {e}")
                return None
        
        return None
        
    def _search_local(self, query: str, ontology_name: str) -> Optional[Tuple[str, str, float]]:
        """Search in locally loaded ontology."""
        ont = self._get_local_ontology(ontology_name)
        if not ont:
            return None
            
        query_lower = query.lower().strip()
        
        # 1. Exact Name/Synonym Match
        for term in ont.terms():
            if term.name and term.name.lower() == query_lower:
                return (term.name, term.id, 1.0)
            
            for synonym in term.synonyms:
                if synonym.description.lower() == query_lower:
                    return (term.name, term.id, 0.95)
                    
        return None

    def _fuzzy_match_score(self, text1: str, text2: str) -> float:
        return SequenceMatcher(None, text1.lower(), text2.lower()).ratio()

    def map_tissue(self, raw_value: str) -> Tuple[str, Optional[str], float]:
        """Map tissue using API -> Local -> Legacy."""
        # 1. API
        result = self.ols_client.search(raw_value, "uberon")
        if result:
            return result
            
        # 2. Local File
        result = self._search_local(raw_value, "uberon")
        if result:
            return result
            
        # 3. Legacy Map
        raw_lower = raw_value.lower().strip()
        if raw_lower in self._legacy_tissue_map:
             name, uberon_id, conf = self._legacy_tissue_map[raw_lower]
             return name, uberon_id, conf * 0.9
             
        return raw_value, None, 0.5

    def map_cell_type(self, raw_value: str) -> Tuple[str, Optional[str], float]:
        """Map cell type using API -> Local -> Legacy."""
        # 1. API
        result = self.ols_client.search(raw_value, "cl")
        if result:
            return result
            
        # 2. Local File
        result = self._search_local(raw_value, "cl")
        if result:
            return result
            
        return raw_value, None, 0.5

    def map_organism(self, raw_value: str) -> Tuple[str, Optional[str], float]:
        """Map organism using API -> (No local for taxon) -> Fallback."""
        result = self.ols_client.search(raw_value, "ncbitaxon")
        if result:
            return result
            
        # NCBITaxon is too big for pronto usually, so skip local loading
        return raw_value, None, 0.5


# Singleton mapper instance
_mapper: Optional[OntologyMapper] = None


def get_mapper() -> OntologyMapper:
    """Get or create singleton OntologyMapper instance."""
    global _mapper
    if _mapper is None:
        _mapper = OntologyMapper()
    return _mapper


def normalize_tissue(raw_value: str) -> Tuple[str, Optional[str], float]:
    """Normalize tissue term using UBERON ontology."""
    return get_mapper().map_tissue(raw_value)


def normalize_cell_type(raw_value: str) -> Tuple[str, Optional[str], float]:
    """Normalize cell type using Cell Ontology."""
    return get_mapper().map_cell_type(raw_value)


def normalize_organism(raw_value: str) -> Tuple[str, Optional[str], float]:
    """Normalize organism using NCBITaxon."""
    return get_mapper().map_organism(raw_value)

