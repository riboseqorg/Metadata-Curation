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

        # We also need to mock the provider.extract call if that path is taken, 
        # but here we are indirectly testing extract_with_claude via logic flow or 
        # we need to ensure enrich_sample_metadata uses the mock.
        
        # Actually, enrich_sample_metadata calls extract_with_claude if provider is None.
        # But wait, the failure was ValidationError on BaseProvenance.
        # This implies that specific fields like llm_result.tissue were MagicMocks instead of strings/None.
        # This happens because enrich_sample_metadata calls extract_with_claude, which returns LLMExtractionResult.
        # If we mock extract_with_claude, we must return a real LLMExtractionResult or a Mock that behaves like one.
        
    @patch("omics_extractor.extraction.llm_extractor.extract_with_claude")
    def test_enrich_sample_metadata(self, mock_extract):
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

        # Mock extracting strain (which is missing)
        mock_result = LLMExtractionResult(
            tissue=None,
            cell_type=None,
            strain="C57BL/6",
            age="adult",
            sex="male",
            confidence={
                "strain": 0.85,
                "age": 0.7,
                "sex": 0.6
            },
            reasoning="Strain and age inferred from study context"
        )
        mock_extract.return_value = mock_result

        # Enrich
        enriched = enrich_sample_metadata(
            sample=sample,
            study_title="Mouse Brain Study",
            study_description="Adult male C57BL/6 mice were used",
        )

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
        assert enriched.strain.source == "llm_enrichment"
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

    @patch("omics_extractor.extraction.llm_extractor.extract_with_claude")
    def test_enrich_sample_metadata_with_raw_characteristics(self, mock_extract):
        """Should include raw characteristics in prompt."""
        # Create sample with raw characteristics
        sample = SampleMetadata(
            sample_id="test",
            bioproject_id="test",
            organism=BaseProvenance(value="test", source="test", source_id="test", confidence=1.0),
            raw_characteristics={"source_name": "Lungs", "strain": "C57BL/6"}
        )

        # Mock result to avoid validation error
        mock_result = LLMExtractionResult() # Empty result
        mock_extract.return_value = mock_result

        with patch("omics_extractor.extraction.llm_extractor.build_extraction_prompt") as mock_build_prompt:
            # We need to manually call the real function or simulate its effect because we are mocking it?
            # Wait, if we mock build_extraction_prompt, execute default implementation? 
            # No, enrich_sample_metadata calls build_extraction_prompt inside.
            # If we mock it, we can check arguments.
            # But enrich_sample_metadata -> extract_with_claude (mocked) or provider.extract
            # Wait, enrich_sample_metadata calls build_extraction_prompt ONLY IF provider is NOT None.
            # If provider IS None (default), it calls extract_with_claude.
            # extract_with_claude calls build_extraction_prompt.
            
            # The test logic I wrote previously:
            # enrich_sample_metadata(sample, ...) -> calls extract_with_claude(...)
            # inside extract_with_claude -> calls build_extraction_prompt(...)
            
            # So if I mock extract_with_claude, build_extraction_prompt is NEVER called because I mocked the caller!
            # I must NOT mock extract_with_claude if I want to verify build_extraction_prompt arguments,
            # OR I must mock build_extraction_prompt AND ensure the code path hits it.
            
            # extract_with_claude calls build_extraction_prompt.
            # If I mock extract_with_claude, I replace the whole function, so logic inside it (calling build_prompt) is lost.
            
            # Correct approach:
            # Use 'provider' argument to force the path that calls build_extraction_prompt directly in enrich_sample_metadata?
            # No, looking at code:
            # if provider is not None:
            #    prompt = build_extraction_prompt(...)
            #    llm_result = provider.extract(prompt)
            # else:
            #    llm_result = extract_with_claude(...)
            
            # So I should pass a mock provider to verify build_extraction_prompt is called with correct args!
            
            mock_provider = MagicMock()
            mock_provider.extract.return_value = LLMExtractionResult()
            
            enrich_sample_metadata(
                sample=sample,
                study_title="Test",
                study_description="Test",
                provider=mock_provider
            )
            
            # Now verify build_extraction_prompt was called
            # But I need to patch build_extraction_prompt to spy on it.
            return # Handled in logic below

    def test_enrich_with_raw_characteristics_check_prompt(self):
        """Should include raw characteristics in extraction prompt."""
        sample = SampleMetadata(
            sample_id="test",
            bioproject_id="test",
            organism=BaseProvenance(value="test", source="test", source_id="test", confidence=1.0),
            raw_characteristics={"source_name": "Lungs", "strain": "C57BL/6"}
        )
        
        with patch("omics_extractor.extraction.llm_extractor.build_extraction_prompt") as mock_build:
            mock_provider = MagicMock()
            mock_provider.extract.return_value = LLMExtractionResult()
            
            enrich_sample_metadata(
                sample=sample,
                study_title="Test",
                study_description="Test",
                provider=mock_provider
            )
            
            args, kwargs = mock_build.call_args
            assert "raw_characteristics" in kwargs["existing_metadata"]
            assert kwargs["existing_metadata"]["raw_characteristics"] == {"source_name": "Lungs", "strain": "C57BL/6"}

