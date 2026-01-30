"""Test dynamic extraction schemes."""
import pytest
from omics_extractor.extraction.llm_extractor import load_extraction_scheme, build_extraction_prompt


def test_load_default_scheme():
    """Test loading the default extraction scheme."""
    scheme = load_extraction_scheme("default")
    
    assert "organism" in scheme
    assert "tissue" in scheme
    assert "cell_type" in scheme
    assert "library_strategy" in scheme


def test_load_ribo_seq_scheme():
    """Test loading the Ribo-Seq extraction scheme."""
    scheme = load_extraction_scheme("ribo_seq")
    
    # Standard fields
    assert "organism" in scheme
    assert "tissue" in scheme
    
    # Ribo-Seq specific fields
    assert "inhibitor" in scheme
    assert "nuclease" in scheme
    assert "digestion_temperature" in scheme
    assert "digestion_time" in scheme


def test_load_nanopore_scheme():
    """Test loading the Nanopore RNA extraction scheme."""
    scheme = load_extraction_scheme("nanopore_rna")
    
    # Standard fields
    assert "organism" in scheme
    assert "tissue" in scheme
    
    # Nanopore specific fields
    assert "flow_cell_type" in scheme
    assert "kit_type" in scheme
    assert "basecalling_model" in scheme


def test_build_prompt_with_custom_scheme():
    """Test that build_extraction_prompt uses custom schemes."""
    ribo_scheme = load_extraction_scheme("ribo_seq")
    
    prompt = build_extraction_prompt(
        study_title="Test Ribo-Seq Study",
        study_description="A study using cycloheximide and RNase I",
        fields=ribo_scheme
    )
    
    # Verify that custom fields appear in the prompt
    assert "inhibitor" in prompt
    assert "nuclease" in prompt
    assert "digestion_temperature" in prompt


def test_fallback_to_default():
    """Test that non-existent schemes fall back to default."""
    scheme = load_extraction_scheme("nonexistent_scheme")
    
    # Should still return a valid scheme (default)
    assert "organism" in scheme
    assert "tissue" in scheme
