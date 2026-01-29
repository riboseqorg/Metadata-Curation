import json
import sys
from pathlib import Path
from datetime import datetime

# Add src to path
sys.path.insert(0, str(Path.cwd() / "src"))

from omics_extractor.schemas.base import SampleMetadata, BaseProvenance
from omics_extractor.extraction.llm_schemas import LLMExtractionResult
from omics_extractor.extraction.llm_extractor import enrich_sample_metadata
from omics_extractor.output.traceable_format import create_traceable_report, format_study_traceable

def verify_synchronization():
    print("Testing Metadata Synchronization...")
    
    # 1. Create a dummy SampleMetadata with only organism
    sample = SampleMetadata(
        sample_id="SAMN12345",
        bioproject_id="PRJNA67890",
        organism=BaseProvenance(
            value="Mus musculus",
            source="biosample",
            source_id="SAMN12345",
            confidence=1.0,
            extraction_method="structured_field"
        )
    )
    
    # 2. Mock LLM result with ALL fields
    llm_result = LLMExtractionResult(
        organism="Mus musculus",
        tissue="brain",
        cell_line="Neuro-2a",
        cell_type="neuron",
        strain="C57BL/6",
        genotype="Wild-type",
        sex="female",
        age="8 weeks",
        developmental_stage="adult",
        condition="heat stress",
        treatment="caffeine",
        timepoint="4h",
        replicate="1",
        batch="B1",
        disease="none",
        stress="heat",
        temperature="42C",
        growth_condition="standard",
        confidence={f: 0.9 for f in ["tissue", "cell_line", "cell_type", "strain", "condition"]},
        reasoning="Extracted all fields for test."
    )
    
    print("  ✓ Created mock LLM result with 15+ fields")
    
    # 3. Enrich sample (this tests llm_extractor.py logic)
    # We mock the provider by just passing the result to a helper if needed, 
    # but enrich_sample_metadata takes the sample and llm_result directly in some versions?
    # Actually, enrich_sample_metadata calls provider.extract. 
    # Let's look at it again.
    
    from unittest.mock import MagicMock
    mock_provider = MagicMock()
    mock_provider.extract.return_value = llm_result
    
    enriched_sample = enrich_sample_metadata(
        sample=sample,
        study_title="Test Study",
        study_description="Test Description",
        provider=mock_provider
    )
    
    # 4. Verify enriched_sample has all fields
    fields_to_check = [
        "tissue", "cell_line", "cell_type", "strain", "genotype", 
        "sex", "age", "developmental_stage", "condition", "treatment", 
        "timepoint", "replicate", "batch", "disease", "stress", 
        "temperature", "growth_condition"
    ]
    
    missing_in_sample = []
    for field in fields_to_check:
        val = getattr(enriched_sample, field, None)
        if not val or not val.value:
            missing_in_sample.append(field)
    
    if missing_in_sample:
        print(f"  ✗ FAILED: Fields missing in enriched SampleMetadata: {missing_in_sample}")
        return False
    else:
        print("  ✓ All fields correctly merged into SampleMetadata")

    # 5. Verify Traceable Report (this tests traceable_format.py)
    from omics_extractor.schemas.base import StudyMetadata
    study = StudyMetadata(
        bioproject_id="PRJNA67890",
        title=BaseProvenance(value="Test Study", source="bioproject", source_id="PRJNA67890", confidence=1.0),
        description=BaseProvenance(value="Test Description", source="bioproject", source_id="PRJNA67890", confidence=1.0),
        organism=BaseProvenance(value="Mus musculus", source="bioproject", source_id="PRJNA67890", confidence=1.0)
    )
    
    report = create_traceable_report(
        study=study,
        samples={"SAMN12345": enriched_sample},
        include_statistics=True
    )
    
    # Check quick_view
    quick_view = report["samples"]["SAMN12345"]["quick_view"]
    missing_in_quick = []
    for field in fields_to_check:
        if field not in quick_view:
            missing_in_quick.append(field)
            
    if missing_in_quick:
        print(f"  ✗ FAILED: Fields missing in report quick_view: {missing_in_quick}")
        return False
    else:
        print("  ✓ All fields present in report quick_view")
        
    # Check statistics coverage
    coverage = report["extraction_statistics"]["field_coverage"]
    missing_in_stats = []
    for field in fields_to_check:
        if field not in coverage:
            missing_in_stats.append(field)
            
    if missing_in_stats:
        print(f"  ✗ FAILED: Fields missing in report coverage stats: {missing_in_stats}")
        return False
    else:
        print("  ✓ All fields present in report coverage statistics")

    print("\n✓ SUCCESS: Metadata synchronization verified across all layers!")
    return True

if __name__ == "__main__":
    if verify_synchronization():
        sys.exit(0)
    else:
        sys.exit(1)
