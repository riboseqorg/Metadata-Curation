#!/usr/bin/env python3
"""
Look up ontology information for terms including synonyms and relationships.

This script provides a command-line interface to explore ontology mappings,
find synonyms, and discover related terms for tissues, cell types, and organisms.

Usage:
    # Look up a tissue term
    python scripts/ontology_lookup.py tissue "liver"

    # Look up a cell type
    python scripts/ontology_lookup.py cell_type "T cell"

    # Look up an organism
    python scripts/ontology_lookup.py organism "mouse"

    # Look up multiple terms
    python scripts/ontology_lookup.py tissue "brain" "liver" "heart"

    # Show all mapped terms for a category
    python scripts/ontology_lookup.py tissue --list-all
"""

import sys
from pathlib import Path
from typing import List, Tuple, Optional
import argparse

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from omics_extractor.ontologies.mapper import OntologyMapper


def lookup_term(mapper: OntologyMapper, category: str, term: str) -> None:
    """
    Look up a term and display its ontology information.

    Args:
        mapper: OntologyMapper instance
        category: 'tissue', 'cell_type', or 'organism'
        term: Term to look up
    """
    # Get the appropriate mapping function
    if category == "tissue":
        normalizer = mapper.map_tissue
        ontology_name = "UBERON"
    elif category == "cell_type":
        normalizer = mapper.map_cell_type
        ontology_name = "Cell Ontology (CL)"
    elif category == "organism":
        normalizer = mapper.map_organism
        ontology_name = "NCBITaxon"
    else:
        print(f"Unknown category: {category}", file=sys.stderr)
        return

    # Normalize the term
    normalized, ontology_id, confidence = normalizer(term)

    print(f"\n{'=' * 70}")
    print(f"Term: {term}")
    print(f"{'=' * 70}")

    if ontology_id:
        print(f"✓ Mapped to {ontology_name}")
        print(f"  Normalized: {normalized}")
        print(f"  Ontology ID: {ontology_id}")
        print(f"  Confidence: {confidence:.2f}")

        # Explain confidence level
        if confidence == 1.0:
            explanation = "Exact match to ontology term"
        elif confidence >= 0.95:
            explanation = "Case-insensitive match"
        elif confidence >= 0.9:
            explanation = "Match to synonym"
        elif confidence >= 0.7:
            explanation = "Fuzzy match"
        else:
            explanation = "Low confidence match"

        print(f"  Quality: {explanation}")
    else:
        print(f"✗ No mapping found in {ontology_name}")
        print(f"  Original value kept: {normalized}")
        print(f"  Confidence: {confidence:.2f}")

        # Suggest similar terms
        print(f"\n  Try looking up similar terms or check the ontology directly:")
        if category == "tissue":
            print(f"    UBERON: https://www.ebi.ac.uk/ols/ontologies/uberon")
        elif category == "cell_type":
            print(f"    Cell Ontology: https://www.ebi.ac.uk/ols/ontologies/cl")
        elif category == "organism":
            print(f"    NCBITaxon: https://www.ncbi.nlm.nih.gov/taxonomy")


def list_all_mappings(mapper: OntologyMapper, category: str) -> None:
    """
    List all available mappings for a category.

    Args:
        mapper: OntologyMapper instance
        category: 'tissue', 'cell_type', or 'organism'
    """
    # Access the internal mapping dictionaries
    # Note: This is a bit hacky but useful for exploration

    if category == "tissue":
        # Extract from map_tissue method source
        print("\n" + "=" * 70)
        print("AVAILABLE TISSUE MAPPINGS (UBERON)")
        print("=" * 70)
        print(f"\n{'Term':<25} {'Normalized':<25} {'Ontology ID':<20}")
        print("-" * 70)

        # Test common terms
        test_terms = [
            "brain", "cerebral cortex", "hippocampus", "cerebellum", "spinal cord",
            "liver", "kidney", "heart", "lung", "stomach", "intestine",
            "skeletal muscle", "cardiac muscle", "bone",
            "blood", "spleen", "thymus", "bone marrow",
            "skin", "adipose tissue", "fat",
            "testis", "ovary", "prostate", "breast",
            "pancreas", "eye", "retina", "placenta",
        ]

        for term in test_terms:
            normalized, ontology_id, _ = mapper.map_tissue(term)
            if ontology_id:
                print(f"{term:<25} {normalized:<25} {ontology_id:<20}")

    elif category == "cell_type":
        print("\n" + "=" * 70)
        print("AVAILABLE CELL TYPE MAPPINGS (Cell Ontology)")
        print("=" * 70)
        print(f"\n{'Term':<25} {'Normalized':<25} {'Ontology ID':<20}")
        print("-" * 70)

        test_terms = [
            "neuron", "astrocyte", "oligodendrocyte", "microglia",
            "T cell", "B cell", "macrophage", "monocyte", "neutrophil",
            "stem cell", "embryonic stem cell",
            "epithelial cell", "endothelial cell", "fibroblast",
            "muscle cell", "cardiomyocyte",
            "hepatocyte", "adipocyte", "osteoblast",
        ]

        for term in test_terms:
            normalized, ontology_id, _ = mapper.map_cell_type(term)
            if ontology_id:
                print(f"{term:<25} {normalized:<25} {ontology_id:<20}")

    elif category == "organism":
        print("\n" + "=" * 70)
        print("AVAILABLE ORGANISM MAPPINGS (NCBITaxon)")
        print("=" * 70)
        print(f"\n{'Term':<25} {'Normalized':<25} {'Taxonomy ID':<20}")
        print("-" * 70)

        test_terms = [
            "human", "homo sapiens",
            "mouse", "mus musculus",
            "rat", "rattus norvegicus",
            "zebrafish", "danio rerio",
            "fruit fly", "drosophila melanogaster",
            "c. elegans", "caenorhabditis elegans",
            "yeast", "saccharomyces cerevisiae",
            "e. coli", "escherichia coli",
        ]

        for term in test_terms:
            normalized, ontology_id, _ = mapper.map_organism(term)
            if ontology_id:
                print(f"{term:<25} {normalized:<25} {ontology_id:<20}")

    print(f"\nTotal mappings shown: {len([t for t in test_terms if mapper.map_tissue(t)[1] if category == 'tissue' else mapper.map_cell_type(t)[1] if category == 'cell_type' else mapper.map_organism(t)[1]])}")
    print("\nNote: These are curated common terms. For comprehensive coverage,")
    print("      consider downloading full ontologies (UBERON, CL, NCBITaxon).")


def main():
    parser = argparse.ArgumentParser(
        description="Look up ontology mappings and synonyms",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Look up specific terms
  %(prog)s tissue "liver" "brain" "heart"
  %(prog)s cell_type "T cell" "macrophage"
  %(prog)s organism "mouse" "human"

  # List all available mappings
  %(prog)s tissue --list-all
  %(prog)s cell_type --list-all
  %(prog)s organism --list-all
        """
    )
    parser.add_argument("category", choices=["tissue", "cell_type", "organism"],
                       help="Category to look up")
    parser.add_argument("terms", nargs="*", help="Terms to look up")
    parser.add_argument("--list-all", action="store_true",
                       help="List all available mappings for this category")

    args = parser.parse_args()

    # Initialize mapper
    mapper = OntologyMapper()

    if args.list_all:
        list_all_mappings(mapper, args.category)
    elif args.terms:
        for term in args.terms:
            lookup_term(mapper, args.category, term)
    else:
        print("Error: Provide terms to look up or use --list-all", file=sys.stderr)
        parser.print_help()
        sys.exit(1)


if __name__ == "__main__":
    main()
