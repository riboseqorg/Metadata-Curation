"""Tests for LLM extraction."""

import pytest
import os
from unittest.mock import patch, MagicMock

from omics_extractor.extraction.llm_extractor import (
    build_extraction_prompt,
    extract_with_claude,
    enrich_sample_metadata,
    LLMExtractionResult,
)
from omics_extractor.schemas.base import SampleMetadata, BaseProvenance


class TestPromptBuilder:
    """Test prompt generation."""

    def test_build_basic_prompt(self):
        """Should build a basic extraction prompt."""
        prompt = build_extraction_prompt(
            study_title="Test Study",
            study_description="This is a test study about brain tissue.",
        )

        assert "Test Study" in prompt
        assert "brain tissue" in prompt
        assert "TASK:" in prompt
        assert "tissue" in prompt

    def test_build_prompt_with_sample_info(self):
        """Should include sample information in prompt."""
        prompt = build_extraction_prompt(
            study_title="Test Study",
            study_description="Test description",
            sample_title="Sample 1",
            sample_description="Brain sample from adult mouse",
        )

        assert "Sample 1" in prompt
        assert "Brain sample from adult mouse" in prompt

    def test_build_prompt_with_abstract(self):
        """Should include abstract in prompt."""
        prompt = build_extraction_prompt(
            study_title="Test Study",
            study_description="Test description",
            abstract="This study examines hippocampal neurons in C57BL/6 mice."
        )

        assert "hippocampal neurons" in prompt
        assert "C57BL/6" in prompt

    def test_build_prompt_with_existing_metadata(self):
        """Should show existing metadata to avoid duplication."""
        existing = {"tissue": "brain", "strain": "C57BL/6"}
        prompt = build_extraction_prompt(
            study_title="Test Study",
            study_description="Test description",
            existing_metadata=existing,
        )

        assert "ALREADY EXTRACTED" in prompt
        assert "brain" in prompt
        assert "C57BL/6" in prompt


class TestLLMExtraction:
    """Test LLM extraction (mocked)."""

    def test_extract_with_claude_no_api_key(self):
        """Should raise error if no API key provided."""
        with pytest.raises(ValueError, match="API key required"):
            extract_with_claude(
                study_title="Test",
                study_description="Test",
                api_key=None,
            )

    @patch.dict(os.environ, {"ANTHROPIC_API_KEY": "test-key"})
    @patch("omics_extractor.extraction.llm_extractor.anthropic.Anthropic")
    def test_extract_with_claude_success(self, mock_anthropic):
        """Should extract metadata using Claude API."""
        # Mock Claude response
        mock_client = MagicMock()
        mock_anthropic.return_value = mock_client

        mock_message = MagicMock()
        mock_message.content = [MagicMock()]
        mock_message.content[0].text = '''```json
{
  "tissue": "brain",
  "cell_type": "neuron",
  "cell_line": null,
  "treatment": null,
  "genotype": null,
  "strain": "C57BL/6",
  "age": "adult",
  "sex": null,
  "developmental_stage": null,
  "confidence": {
    "tissue": 0.9,
    "cell_type": 0.8,
    "strain": 0.85,
    "age": 0.7
  },
  "reasoning": "Tissue and strain explicitly mentioned in study description"
}
```'''
        mock_client.messages.create.return_value = mock_message

        result = extract_with_claude(
            study_title="Brain Study",
            study_description="Examining hippocampal neurons in adult C57BL/6 mice",
        )

        assert result.tissue == "brain"
        assert result.cell_type == "neuron"
        assert result.strain == "C57BL/6"
        assert result.age == "adult"
        assert result.confidence["tissue"] == 0.9
        assert result.confidence["cell_type"] == 0.8

    @patch.dict(os.environ, {"ANTHROPIC_API_KEY": "test-key"})
    @patch("omics_extractor.extraction.llm_extractor.anthropic.Anthropic")
    def test_enrich_sample_metadata(self, mock_anthropic):
        """Should enrich sample metadata without overwriting existing fields."""
        # Create sample with tissue already extracted
        sample = SampleMetadata(
            sample_id="test",
            bioproject_id="test",
            organism=BaseProvenance(
                value="Mus musculus",
                source="biosample",
                source_id="test",
                confidence=1.0,
            ),
            tissue=BaseProvenance(
                value="brain",
                source="biosample",
                source_id="test",
                confidence=1.0,
            ),
        )

        # Mock Claude to extract strain (which is missing)
        mock_client = MagicMock()
        mock_anthropic.return_value = mock_client

        mock_message = MagicMock()
        mock_message.content = [MagicMock()]
        mock_message.content[0].text = '''
{
  "tissue": null,
  "cell_type": null,
  "cell_line": null,
  "treatment": null,
  "genotype": null,
  "strain": "C57BL/6",
  "age": "adult",
  "sex": "male",
  "developmental_stage": null,
  "confidence": {
    "strain": 0.85,
    "age": 0.7,
    "sex": 0.6
  },
  "reasoning": "Strain and age inferred from study context"
}
'''
        mock_client.messages.create.return_value = mock_message

        # Enrich
        enriched = enrich_sample_metadata(
            sample=sample,
            study_title="Mouse Brain Study",
            study_description="Adult male C57BL/6 mice were used",
        )

        # Original tissue should remain unchanged
        assert enriched.tissue.value == "brain"
        assert enriched.tissue.source == "biosample"
        assert enriched.tissue.confidence == 1.0

        # New strain should be added from LLM
        assert enriched.strain is not None
        assert enriched.strain.value == "C57BL/6"
        assert enriched.strain.source == "llm"
        # Field confidence (0.5) × normalization confidence (0.85) = 0.425
        assert enriched.strain.confidence == pytest.approx(0.5 * 0.85)
        assert enriched.strain.field_confidence == 0.5
        assert enriched.strain.normalization_confidence == 0.85

        # Age should be added
        assert enriched.age is not None
        assert enriched.age.value == "adult"
        assert enriched.age.confidence == pytest.approx(0.5 * 0.7)

        # Sex should be added
        assert enriched.sex is not None
        assert enriched.sex.value == "male"
        assert enriched.sex.confidence == pytest.approx(0.5 * 0.6)
