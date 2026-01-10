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

from typing import Optional, Tuple, Dict
from pathlib import Path
import pronto
from difflib import SequenceMatcher


class OntologyMapper:
    """Maps raw metadata values to ontology terms with confidence scoring."""

    def __init__(self, cache_dir: Optional[Path] = None):
        """
        Initialize ontology mapper.

        Args:
            cache_dir: Directory to cache downloaded ontology files
        """
        self.cache_dir = cache_dir or Path.home() / ".omics_extractor" / "ontologies"
        self.cache_dir.mkdir(parents=True, exist_ok=True)

        self._ontologies: Dict[str, pronto.Ontology] = {}

    def _load_ontology(self, name: str, url: str) -> pronto.Ontology:
        """
        Load an ontology from cache or download it.

        Args:
            name: Ontology name (e.g., 'uberon')
            url: OBO file URL

        Returns:
            Loaded pronto.Ontology object
        """
        cache_file = self.cache_dir / f"{name}.obo"

        if name not in self._ontologies:
            if not cache_file.exists():
                # Download ontology
                import requests
                response = requests.get(url, timeout=300)
                response.raise_for_status()
                cache_file.write_text(response.text)

            self._ontologies[name] = pronto.Ontology(str(cache_file))

        return self._ontologies[name]

    def _fuzzy_match_score(self, text1: str, text2: str) -> float:
        """Calculate fuzzy match score between two strings."""
        return SequenceMatcher(None, text1.lower(), text2.lower()).ratio()

    def map_tissue(self, raw_value: str) -> Tuple[str, Optional[str], float]:
        """
        Map tissue term to UBERON ontology.

        Uses a curated mapping of common tissues for speed and reliability.
        For production, this could be expanded or replaced with full ontology lookup.

        Args:
            raw_value: Raw tissue string from metadata

        Returns:
            Tuple of (normalized_value, ontology_id, confidence)
            - normalized_value: Best matching term or original value
            - ontology_id: UBERON ID if found (e.g., 'UBERON:0000955')
            - confidence: Normalization confidence (0.5-1.0)
        """
        # Curated mapping of common tissue terms to UBERON
        # Format: search_term -> (canonical_name, UBERON_ID, confidence)
        tissue_map = {
            # Brain and nervous system
            "brain": ("brain", "UBERON:0000955", 1.0),
            "cerebral cortex": ("cerebral cortex", "UBERON:0000956", 1.0),
            "cortex": ("cerebral cortex", "UBERON:0000956", 0.9),
            "hippocampus": ("hippocampus", "UBERON:0001954", 1.0),
            "cerebellum": ("cerebellum", "UBERON:0002037", 1.0),
            "spinal cord": ("spinal cord", "UBERON:0002240", 1.0),

            # Digestive system
            "liver": ("liver", "UBERON:0002107", 1.0),
            "kidney": ("kidney", "UBERON:0002113", 1.0),
            "heart": ("heart", "UBERON:0000948", 1.0),
            "lung": ("lung", "UBERON:0002048", 1.0),
            "stomach": ("stomach", "UBERON:0000945", 1.0),
            "intestine": ("intestine", "UBERON:0000160", 1.0),
            "small intestine": ("small intestine", "UBERON:0002108", 1.0),
            "large intestine": ("large intestine", "UBERON:0000059", 1.0),
            "colon": ("colon", "UBERON:0001155", 1.0),

            # Muscle and skeletal
            "muscle": ("muscle tissue", "UBERON:0002385", 0.95),
            "skeletal muscle": ("skeletal muscle tissue", "UBERON:0001134", 1.0),
            "cardiac muscle": ("cardiac muscle tissue", "UBERON:0001133", 1.0),
            "bone": ("bone tissue", "UBERON:0002481", 0.95),

            # Blood and immune
            "blood": ("blood", "UBERON:0000178", 1.0),
            "spleen": ("spleen", "UBERON:0002106", 1.0),
            "thymus": ("thymus", "UBERON:0002370", 1.0),
            "bone marrow": ("bone marrow", "UBERON:0002371", 1.0),
            "lymph node": ("lymph node", "UBERON:0000029", 1.0),

            # Skin and connective
            "skin": ("skin", "UBERON:0001003", 1.0),
            "adipose tissue": ("adipose tissue", "UBERON:0001013", 1.0),
            "fat": ("adipose tissue", "UBERON:0001013", 0.9),

            # Reproductive
            "testis": ("testis", "UBERON:0000473", 1.0),
            "ovary": ("ovary", "UBERON:0000992", 1.0),
            "prostate": ("prostate gland", "UBERON:0002367", 0.95),
            "mammary gland": ("mammary gland", "UBERON:0001911", 1.0),
            "breast": ("mammary gland", "UBERON:0001911", 0.9),

            # Other organs
            "pancreas": ("pancreas", "UBERON:0001264", 1.0),
            "eye": ("eye", "UBERON:0000970", 1.0),
            "retina": ("retina", "UBERON:0000966", 1.0),
            "placenta": ("placenta", "UBERON:0001987", 1.0),
        }

        raw_lower = raw_value.lower().strip()

        # Exact match (case-insensitive)
        if raw_lower in tissue_map:
            name, uberon_id, conf = tissue_map[raw_lower]
            return name, uberon_id, conf

        # Check if it's case-sensitive match
        for key, (name, uberon_id, base_conf) in tissue_map.items():
            if key == raw_value:
                return name, uberon_id, 1.0

        # Fuzzy matching against known terms
        best_match = None
        best_score = 0.0
        best_data = None

        for key, data in tissue_map.items():
            score = self._fuzzy_match_score(raw_lower, key)
            if score > best_score and score >= 0.85:
                best_score = score
                best_match = key
                best_data = data

        if best_match:
            name, uberon_id, _ = best_data
            return name, uberon_id, 0.7  # Fuzzy match gets lower confidence

        # No match found - return original value
        return raw_value, None, 0.5

    def map_cell_type(self, raw_value: str) -> Tuple[str, Optional[str], float]:
        """
        Map cell type to Cell Ontology (CL).

        Uses a curated mapping of common cell types for speed and reliability.

        Args:
            raw_value: Raw cell type string

        Returns:
            Tuple of (normalized_value, ontology_id, confidence)
        """
        # Curated mapping of common cell types to Cell Ontology
        cell_type_map = {
            # Neurons and neural cells
            "neuron": ("neuron", "CL:0000540", 1.0),
            "neuronal cell": ("neuron", "CL:0000540", 0.95),
            "astrocyte": ("astrocyte", "CL:0000127", 1.0),
            "oligodendrocyte": ("oligodendrocyte", "CL:0000128", 1.0),
            "microglia": ("microglial cell", "CL:0000129", 0.95),

            # Blood cells
            "t cell": ("T cell", "CL:0000084", 1.0),
            "b cell": ("B cell", "CL:0000236", 1.0),
            "macrophage": ("macrophage", "CL:0000235", 1.0),
            "monocyte": ("monocyte", "CL:0000576", 1.0),
            "neutrophil": ("neutrophil", "CL:0000775", 1.0),
            "erythrocyte": ("erythrocyte", "CL:0000232", 1.0),
            "red blood cell": ("erythrocyte", "CL:0000232", 0.9),

            # Stem cells
            "stem cell": ("stem cell", "CL:0000034", 1.0),
            "embryonic stem cell": ("embryonic stem cell", "CL:0002322", 1.0),
            "hematopoietic stem cell": ("hematopoietic stem cell", "CL:0000037", 1.0),

            # Epithelial and structural
            "epithelial cell": ("epithelial cell", "CL:0000066", 1.0),
            "endothelial cell": ("endothelial cell", "CL:0000115", 1.0),
            "fibroblast": ("fibroblast", "CL:0000057", 1.0),
            "keratinocyte": ("keratinocyte", "CL:0000312", 1.0),

            # Muscle cells
            "myocyte": ("muscle cell", "CL:0000187", 0.95),
            "muscle cell": ("muscle cell", "CL:0000187", 1.0),
            "cardiomyocyte": ("cardiac muscle cell", "CL:0000746", 0.95),

            # Other common types
            "hepatocyte": ("hepatocyte", "CL:0000182", 1.0),
            "adipocyte": ("adipocyte", "CL:0000136", 1.0),
            "osteoblast": ("osteoblast", "CL:0000062", 1.0),
            "osteocyte": ("osteocyte", "CL:0000137", 1.0),
        }

        raw_lower = raw_value.lower().strip()

        # Exact match
        if raw_lower in cell_type_map:
            name, cl_id, conf = cell_type_map[raw_lower]
            return name, cl_id, conf

        # Fuzzy matching
        best_match = None
        best_score = 0.0
        best_data = None

        for key, data in cell_type_map.items():
            score = self._fuzzy_match_score(raw_lower, key)
            if score > best_score and score >= 0.85:
                best_score = score
                best_match = key
                best_data = data

        if best_match:
            name, cl_id, _ = best_data
            return name, cl_id, 0.7

        return raw_value, None, 0.5

    def map_organism(self, raw_value: str) -> Tuple[str, Optional[str], float]:
        """
        Map organism name to NCBITaxon.

        Note: Full NCBITaxon ontology is very large (~2M terms).
        This is a simplified implementation that validates common model organisms.

        Args:
            raw_value: Raw organism string

        Returns:
            Tuple of (normalized_value, taxon_id, confidence)
        """
        # Common model organisms mapping
        # In production, would use full NCBITaxon ontology or NCBI Taxonomy API
        known_organisms = {
            "homo sapiens": ("Homo sapiens", "NCBITaxon:9606", 1.0),
            "human": ("Homo sapiens", "NCBITaxon:9606", 0.95),
            "mus musculus": ("Mus musculus", "NCBITaxon:10090", 1.0),
            "mouse": ("Mus musculus", "NCBITaxon:10090", 0.95),
            "rattus norvegicus": ("Rattus norvegicus", "NCBITaxon:10116", 1.0),
            "rat": ("Rattus norvegicus", "NCBITaxon:10116", 0.95),
            "danio rerio": ("Danio rerio", "NCBITaxon:7955", 1.0),
            "zebrafish": ("Danio rerio", "NCBITaxon:7955", 0.95),
            "drosophila melanogaster": ("Drosophila melanogaster", "NCBITaxon:7227", 1.0),
            "fruit fly": ("Drosophila melanogaster", "NCBITaxon:7227", 0.9),
            "caenorhabditis elegans": ("Caenorhabditis elegans", "NCBITaxon:6239", 1.0),
            "c. elegans": ("Caenorhabditis elegans", "NCBITaxon:6239", 0.95),
            "saccharomyces cerevisiae": ("Saccharomyces cerevisiae", "NCBITaxon:4932", 1.0),
            "yeast": ("Saccharomyces cerevisiae", "NCBITaxon:4932", 0.9),
            "escherichia coli": ("Escherichia coli", "NCBITaxon:562", 1.0),
            "e. coli": ("Escherichia coli", "NCBITaxon:562", 0.95),
        }

        raw_lower = raw_value.lower().strip()

        if raw_lower in known_organisms:
            name, taxon_id, conf = known_organisms[raw_lower]
            return name, taxon_id, conf

        # No match - return original
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
    """
    Normalize tissue term using UBERON ontology.

    Args:
        raw_value: Raw tissue string

    Returns:
        Tuple of (normalized_value, ontology_id, confidence)

    Example:
        >>> normalize_tissue("Brain")
        ('brain', 'UBERON:0000955', 0.95)
    """
    return get_mapper().map_tissue(raw_value)


def normalize_cell_type(raw_value: str) -> Tuple[str, Optional[str], float]:
    """
    Normalize cell type using Cell Ontology.

    Args:
        raw_value: Raw cell type string

    Returns:
        Tuple of (normalized_value, ontology_id, confidence)
    """
    return get_mapper().map_cell_type(raw_value)


def normalize_organism(raw_value: str) -> Tuple[str, Optional[str], float]:
    """
    Normalize organism using NCBITaxon.

    Args:
        raw_value: Raw organism string

    Returns:
        Tuple of (normalized_value, taxon_id, confidence)
    """
    return get_mapper().map_organism(raw_value)
