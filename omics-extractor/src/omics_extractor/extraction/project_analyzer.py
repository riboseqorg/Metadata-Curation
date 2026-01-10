"""Project-level LLM analysis for experimental design characterization."""

from typing import Dict, List, Optional, Any
from pydantic import BaseModel, Field
import json

from ..schemas.base import SampleMetadata, StudyMetadata


class ReplicateGroup(BaseModel):
    """A group of replicate samples."""
    group_id: str
    sample_ids: List[str]
    conditions: Dict[str, Any]
    replicate_type: str  # "biological", "technical", "mixed"


class SampleRelationship(BaseModel):
    """Relationship between two samples."""
    sample_id_1: str
    sample_id_2: str
    relationship_type: str  # "paired_assay", "time_series", "dose_response", etc.
    description: Optional[str] = None


class ProjectAnalysisResult(BaseModel):
    """Result from project-level LLM analysis."""

    # Experimental variables identified
    experimental_variables: List[str] = Field(
        default_factory=list,
        description="Key variables that distinguish samples (e.g., 'treatment', 'timepoint', 'tissue')"
    )

    # Replicate groups
    replicate_groups: List[ReplicateGroup] = Field(
        default_factory=list,
        description="Groups of samples that are replicates"
    )

    # Sample relationships
    relationships: List[SampleRelationship] = Field(
        default_factory=list,
        description="Relationships between samples (e.g., paired Ribo-Seq/RNA-Seq)"
    )

    # Experimental design summary
    design_summary: Optional[str] = Field(
        None,
        description="Natural language summary of the experimental design"
    )

    # Confidence
    confidence: float = Field(
        0.0,
        description="Overall confidence in the analysis (0.0-1.0)"
    )


def build_project_analysis_prompt(
    study: StudyMetadata,
    samples: Dict[str, SampleMetadata],
) -> str:
    """
    Build a prompt for comprehensive project-level analysis.

    This prompt is inspired by experimental design analysis and asks the LLM to:
    1. Identify experimental variables and conditions
    2. Group replicates
    3. Map relationships between samples
    4. Characterize the experimental design

    Args:
        study: Study-level metadata
        samples: Dictionary of sample_id -> SampleMetadata

    Returns:
        Formatted prompt for LLM
    """

    # Build comprehensive metadata representation
    metadata_lines = []

    # Study-level info
    metadata_lines.append("=== STUDY METADATA ===")
    metadata_lines.append(f"BioProject ID: {study.bioproject_id if study.bioproject_id else 'N/A'}")
    metadata_lines.append(f"Title: {study.title.value if study.title else 'N/A'}")
    metadata_lines.append(f"Description: {study.description.value if study.description else 'N/A'}")
    metadata_lines.append("")

    # Sample-level info
    metadata_lines.append("=== SAMPLE METADATA ===")
    metadata_lines.append(f"Total Samples: {len(samples)}")
    metadata_lines.append("")

    for sample_id, sample in samples.items():
        metadata_lines.append(f"Sample: {sample_id}")

        # Add all available fields
        fields_to_show = [
            ("BioSample ID", sample.biosample_id if sample.biosample_id else None),
            ("Sample Title", sample.sample_title.value if sample.sample_title else None),
            ("Sample Description", sample.sample_description.value if sample.sample_description else None),
            ("Organism", sample.organism.value if sample.organism else None),
            ("Tissue", sample.tissue.value if sample.tissue else None),
            ("Cell Type", sample.cell_type.value if sample.cell_type else None),
            ("Cell Line", sample.cell_line.value if sample.cell_line else None),
            ("Treatment", sample.treatment.value if sample.treatment else None),
            ("Strain", sample.strain.value if sample.strain else None),
            ("Genotype", sample.genotype.value if sample.genotype else None),
            ("Age", sample.age.value if sample.age else None),
            ("Sex", sample.sex.value if sample.sex else None),
            ("Developmental Stage", sample.developmental_stage.value if sample.developmental_stage else None),
        ]

        for field_name, field_value in fields_to_show:
            if field_value:
                metadata_lines.append(f"  {field_name}: {field_value}")

        # Add run information
        if hasattr(sample, 'runs') and sample.runs:
            metadata_lines.append(f"  Number of Runs: {len(sample.runs)}")
            for run in sample.runs[:3]:  # Show first 3 runs
                if run.library_strategy:
                    metadata_lines.append(f"    Run {run.run_id.value}: {run.library_strategy.value}")

        metadata_lines.append("")

    metadata_str = "\n".join(metadata_lines)

    prompt = f"""You will be analyzing experimental metadata from a scientific study (likely genomics or sequencing) to characterize the experimental design and understand relationships between samples.

Here is the metadata to analyze:

<metadata>
{metadata_str}
</metadata>

Your goal is to produce a structured, machine-readable output that:
1. Identifies experimental variables and conditions that distinguish samples
2. Groups samples that are replicates (identical or nearly identical conditions)
3. Maps relationships between related samples (e.g., which Ribo-Seq sample corresponds to which RNA-Seq sample)
4. Characterizes the overall experimental design

## Analysis Process

Before producing your final output, work through your analysis systematically inside <analysis> tags:

1. **List all samples and key fields**: Identify all sample IDs and the metadata fields that are present.

2. **Identify experimental variables**: Determine which metadata fields represent experimental conditions that distinguish samples (e.g., treatment, tissue, timepoint) vs. technical metadata (e.g., batch IDs, run IDs).

3. **Extract condition values**: For each experimental variable, list the unique values across all samples.

4. **Group replicates**: Identify groups of samples with identical (or nearly identical) experimental conditions. These are likely biological or technical replicates.

5. **Map relationships**: Look for samples that are related but not identical:
   - Paired experiments (same sample, different assay: Ribo-Seq vs RNA-Seq)
   - Time series or dose-response experiments
   - Different tissues/cell types from same organism/condition

6. **Characterize design**: Describe the overall experimental design (e.g., "2x3 factorial design with 3 replicates" or "time course with 5 timepoints").

## Output Format

After your analysis, provide your results as valid JSON inside <json_output> tags:

{{
  "experimental_variables": ["list", "of", "key", "variables"],
  "replicate_groups": [
    {{
      "group_id": "rep_group_1",
      "sample_ids": ["sample1", "sample2"],
      "conditions": {{"variable1": "value1", "variable2": "value2"}},
      "replicate_type": "biological"
    }}
  ],
  "relationships": [
    {{
      "sample_id_1": "sample1",
      "sample_id_2": "sample2",
      "relationship_type": "paired_assay",
      "description": "Ribo-Seq and RNA-Seq from same sample"
    }}
  ],
  "design_summary": "Natural language summary of experimental design",
  "confidence": 0.85
}}

IMPORTANT:
- Be conservative - only identify relationships and replicates you're confident about
- If you can't identify clear replicate groups, return empty list
- Confidence should reflect how clear the experimental design is from the metadata
- Focus on biological relationships, not just technical groupings

Return only valid JSON that matches the schema above.
"""

    return prompt


def analyze_project_design(
    study: StudyMetadata,
    samples: Dict[str, SampleMetadata],
    provider,  # LLMProvider instance
) -> ProjectAnalysisResult:
    """
    Analyze project-level experimental design using LLM.

    Args:
        study: Study metadata
        samples: Dictionary of sample metadata
        provider: LLM provider instance (Claude, VLLM, etc.)

    Returns:
        ProjectAnalysisResult with experimental design analysis
    """

    prompt = build_project_analysis_prompt(study, samples)

    # Extract using provider
    llm_result = provider.extract(prompt)

    # Parse the response text to extract JSON
    response_text = llm_result.raw_response if hasattr(llm_result, 'raw_response') else str(llm_result)

    # Extract JSON from tags
    if "<json_output>" in response_text:
        json_text = response_text.split("<json_output>")[1].split("</json_output>")[0].strip()
    elif "```json" in response_text:
        json_text = response_text.split("```json")[1].split("```")[0].strip()
    elif "```" in response_text:
        json_text = response_text.split("```")[1].split("```")[0].strip()
    else:
        json_text = response_text

    result_dict = json.loads(json_text)

    return ProjectAnalysisResult(**result_dict)
