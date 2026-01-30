"""Unified reporting and assessment logic for omics metadata.

Integrates table generation, enrichment assessment, and quality review.
"""

import json
import csv
import sys
import random
from pathlib import Path
from typing import List, Dict, Any, Optional, Tuple
from dataclasses import dataclass, asdict
from datetime import datetime

from ..schemas.base import SampleMetadata, BaseProvenance


# =============================================================================
# 1. TABLE GENERATION LOGIC (from json_to_table.py)
# =============================================================================

def extract_field_value(field_obj: Any) -> tuple[Optional[str], Optional[str], Optional[float], Optional[str]]:
    """Extract value, source, confidence, and ontology from a field."""
    if field_obj is None:
        return None, None, None, None

    if isinstance(field_obj, dict):
        # Recursive extraction in case of nested values or Pydantic-to-dict artifacts
        value = field_obj.get("value")
        if isinstance(value, dict) and "value" in value:
            # Handle rare double-nesting or complex provenance
            value = value.get("value")
            
        # Handle traceable format (nested provenance) or flat format
        prov = field_obj.get("provenance", {}) if "provenance" in field_obj else field_obj
        source = prov.get("source", "")
        confidence = prov.get("confidence")
        ontology = prov.get("ontology_term")
        return str(value) if value is not None else None, source, confidence, ontology
    else:
        # Simple string value
        return str(field_obj) if field_obj is not None else None, None, None, None


def process_sample_to_row(sample_data: dict, study_id: str, library_strategies: Optional[Dict[str, str]] = None) -> Dict[str, Any]:
    """Process a single sample and extract all fields into a flat row."""
    row = {}
    row["study_id"] = study_id
    
    sample_id = None

    # Detect format
    is_traceable = "biological_metadata" in sample_data

    # Canonical fields to extract
    fields_to_extract = [
        "organism", "tissue", "cell_type", "cell_line", "strain",
        "developmental_stage", "genotype", "age", "sex",
        "condition", "treatment", "timepoint", "replicate", "batch",
        "disease", "stress", "temperature", "growth_condition"
    ]

    if is_traceable:
        quick_view = sample_data.get("quick_view", {})
        bio_meta = sample_data.get("biological_metadata", {})
        experimental = sample_data.get("experimental_metadata", {})
        disease_info = sample_data.get("disease_perturbation", {})

        sample_id = quick_view.get("sample_id")
        row["sample_id"] = sample_id
        row["bioproject_id"] = quick_view.get("bioproject_id")
        # Include library strategy if provided
        sra_strategy = library_strategies.get(sample_id) if library_strategies else None
        inferred_strategy = quick_view.get("library_strategy")
        
        # Prioritize inferred if SRA is unhelpful or missing
        if inferred_strategy and (not sra_strategy or str(sra_strategy).lower() in ["other", "unknown", "rna-seq"]):
             row["library_strategy"] = inferred_strategy
        elif sra_strategy:
            row["library_strategy"] = sra_strategy
        else:
            row["library_strategy"] = "unknown"

        for field in fields_to_extract:
            # Look in all traceable categories
            field_obj = bio_meta.get(field) or experimental.get(field) or disease_info.get(field)
            
            value, source, conf, ontology = extract_field_value(field_obj)
            row[field] = value
            row[f"{field}_enriched"] = "yes" if source == "llm_enrichment" else "no"
            row[f"{field}_ontology"] = ontology

        # Descriptions
        desc = sample_data.get("descriptions", {})
        row["sample_title"] = desc.get("title", {}).get("value")
        row["sample_description"] = desc.get("description", {}).get("value")

    else:
        # Flat format
        sample_id = sample_data.get("sample_id")
        row["sample_id"] = sample_id
        row["bioproject_id"] = sample_data.get("bioproject_id")
        
        # Include library strategy if provided
        sra_strategy = library_strategies.get(sample_id) if library_strategies else None
        
        # In flat format, it might be a direct field or nested
        ls_obj = sample_data.get("library_strategy")
        inferred_strategy, _, _, _ = extract_field_value(ls_obj)
        
        # Prioritize inferred if SRA is unhelpful or missing
        if inferred_strategy and (not sra_strategy or str(sra_strategy).lower() in ["other", "unknown", "rna-seq"]):
             row["library_strategy"] = inferred_strategy
        elif sra_strategy:
            row["library_strategy"] = sra_strategy
        else:
            row["library_strategy"] = "unknown"

        for field in fields_to_extract:
            if field in sample_data:
                value, source, conf, ontology = extract_field_value(sample_data.get(field))
                row[field] = value
                row[f"{field}_enriched"] = "yes" if source == "llm_enrichment" else "no"
                row[f"{field}_ontology"] = ontology

        # Descriptions
        st = sample_data.get("sample_title")
        sd = sample_data.get("sample_description")
        row["sample_title"] = st.get("value") if isinstance(st, dict) else st
        row["sample_description"] = sd.get("value") if isinstance(sd, dict) else sd

    return row


def generate_tabular_report(input_files: List[Path], output_path: Path, format_type: str = "tsv"):
    """Generate a tabular CSV/TSV report from multiple JSON files."""
    all_rows = []
    
    for file_path in input_files:
        if not file_path.exists():
            continue
            
        with open(file_path) as f:
            data = json.load(f)
            
        study_id = file_path.stem.replace("_enriched", "").replace("_metadata", "")
        samples = data.get("samples", {})
        runs = data.get("runs", [])
        
        # Build map of sample_id -> library_strategy
        library_map = {}
        for run in runs:
            s_id = run.get("sample_id")
            strategy = run.get("library_strategy")
            if s_id and strategy:
                # If strategy is a dict (BaseProvenance), extract value
                if isinstance(strategy, dict):
                    strategy = strategy.get("value", "unknown")
                
                # If we have multiple runs, prioritize Ribo-seq then RNA-Seq
                current = library_map.get(s_id, "")
                if "ribo" in str(strategy).lower():
                    library_map[s_id] = str(strategy)
                elif not current or "rna" in str(strategy).lower():
                    library_map[s_id] = str(strategy)

        for sample_id, sample_data in samples.items():
            all_rows.append(process_sample_to_row(sample_data, study_id, library_map))

    if not all_rows:
        return False

    # Determine unique columns and order them
    all_columns = set()
    for row in all_rows:
        all_columns.update(row.keys())

    ordered_columns = []
    for col in ["study_id", "bioproject_id", "sample_id", "library_strategy"]:
        if col in all_columns:
            ordered_columns.append(col)
            all_columns.remove(col)

    base_fields = [
        "organism", "tissue", "cell_type", "cell_line", "strain",
        "developmental_stage", "genotype", "age", "sex",
        "condition", "treatment", "timepoint", "replicate", "batch",
        "disease", "stress", "temperature", "growth_condition"
    ]

    for field in base_fields:
        for suffix in ["", "_enriched", "_ontology"]:
            col = f"{field}{suffix}"
            if col in all_columns:
                ordered_columns.append(col)
                all_columns.remove(col)

    for col in ["sample_title", "sample_description"]:
        if col in all_columns:
            ordered_columns.append(col)
            all_columns.remove(col)

    ordered_columns.extend(sorted(all_columns))

    delimiter = "," if format_type == "csv" else "\t"
    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=ordered_columns, delimiter=delimiter, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(all_rows)
        
    return True


# =============================================================================
# 2. ASSESSMENT LOGIC (from assess_enrichment.py)
# =============================================================================

@dataclass
class FieldStats:
    total_samples: int = 0
    filled_before: int = 0
    filled_after: int = 0
    enriched_count: int = 0
    avg_confidence: float = 0.0
    confidence_scores: List[float] = None

    def __post_init__(self):
        if self.confidence_scores is None:
            self.confidence_scores = []


def generate_summary_report(input_files: List[Path]):
    """Generate high-level performance metrics across multiple studies."""
    fields = [
        "organism", "tissue", "cell_type", "cell_line", "strain",
        "developmental_stage", "genotype", "age", "sex",
        "condition", "treatment", "timepoint", "replicate", "batch",
        "disease", "stress", "temperature", "growth_condition"
    ]

    total_samples_global = 0
    global_field_stats = {f: {"before": 0, "after": 0, "enriched": 0, "confidences": []} for f in fields}
    
    study_summaries = []

    for file_path in input_files:
        with open(file_path) as f:
            data = json.load(f)
            
        samples = data.get("samples", {})
        num_samples = len(samples)
        total_samples_global += num_samples
        
        study_enriched = 0
        study_confidences = []

        for sample_id, sample_data in samples.items():
            is_traceable = "biological_metadata" in sample_data
            
            for field in fields:
                if is_traceable:
                    bio = sample_data.get("biological_metadata", {})
                    exp = sample_data.get("experimental_metadata", {})
                    dis = sample_data.get("disease_perturbation", {})
                    field_obj = bio.get(field) or exp.get(field) or dis.get(field)
                else:
                    field_obj = sample_data.get(field)

                if field_obj:
                    val, source, conf, ont = extract_field_value(field_obj)
                    if val:
                        global_field_stats[field]["after"] += 1
                        if source == "llm_enrichment":
                            global_field_stats[field]["enriched"] += 1
                            global_field_stats[field]["confidences"].append(conf or 0.0)
                            study_enriched += 1
                            study_confidences.append(conf or 0.0)
                        else:
                            global_field_stats[field]["before"] += 1
        
        study_summaries.append({
            "id": file_path.stem.replace("_enriched", ""),
            "samples": num_samples,
            "enriched": study_enriched,
            "avg_conf": sum(study_confidences)/len(study_confidences) if study_confidences else 0.0
        })

    # Print output
    print("\n" + "=" * 80)
    print(f"{'ENRICHMENT PERFORMANCE SUMMARY':^80}")
    print("=" * 80)
    print(f"\nOverall: {len(input_files)} studies, {total_samples_global} samples")
    print(f"{'Field':<20} {'Coverage Before':<18} {'Coverage After':<18} {'Enriched':<12} {'Avg Conf':<15}")
    print("-" * 85)

    for field in fields:
        stats = global_field_stats[field]
        cov_before = (stats["before"] / total_samples_global * 100) if total_samples_global > 0 else 0
        cov_after = (stats["after"] / total_samples_global * 100) if total_samples_global > 0 else 0
        avg_conf = sum(stats["confidences"]) / len(stats["confidences"]) if stats["confidences"] else 0
        
        print(f"{field:<20} {cov_before:>6.1f}% ({stats['before']:>3}) → {cov_after:>6.1f}% ({stats['after']:>3})   +{stats['enriched']:<4} ({cov_after-cov_before:>+5.1f}%)   {avg_conf:.3f}")

    print("\nPer-Study Breakdown:")
    print(f"{'Study ID':<30} {'Samples':<10} {'Enriched':<12} {'Avg Conf':<12}")
    print("-" * 65)
    for s in sorted(study_summaries, key=lambda x: x["enriched"], reverse=True):
        print(f"{s['id']:<30} {s['samples']:<10} {s['enriched']:<12} {s['avg_conf']:.3f}")


# =============================================================================
# 3. REVIEW LOGIC (from review_enrichments.py)
# =============================================================================

@dataclass
class EnrichedField:
    field_name: str
    value: str
    confidence: float
    sample_id: str
    study_id: str
    context: Dict[str, str]

    def appears_in_source(self) -> bool:
        all_text = " ".join(self.context.values()).lower()
        return str(self.value).lower() in all_text if self.value else False


def perform_quality_review(input_files: List[Path], num_to_review: int = 10):
    """Interactively review a sample of enriched fields against source text."""
    all_enriched = []
    
    for file_path in input_files:
        with open(file_path) as f:
            data = json.load(f)
            
        study_id = file_path.stem.replace("_enriched", "")
        study = data.get("study", {})
        context = {
            "st": study.get("title", {}).get("value", ""),
            "sd": study.get("description", {}).get("value", ""),
            "sa": study.get("publication", {}).get("abstract", {}).get("value", "")
        }
        
        fields_to_check = [
            "organism", "tissue", "cell_type", "cell_line", "strain",
            "developmental_stage", "genotype", "age", "sex",
            "condition", "treatment", "timepoint", "replicate", "batch",
            "disease", "stress", "temperature", "growth_condition"
        ]

        for sample_id, sample_data in data.get("samples", {}).items():
            is_traceable = "biological_metadata" in sample_data
            
            # Add sample-specific context
            s_ctx = context.copy()
            if is_traceable:
                s_ctx["title"] = sample_data.get("descriptions", {}).get("title", {}).get("value", "")
                s_ctx["desc"] = sample_data.get("descriptions", {}).get("description", {}).get("value", "")
                bio = sample_data.get("biological_metadata", {})
                exp = sample_data.get("experimental_metadata", {})
                dis = sample_data.get("disease_perturbation", {})
                
                for f in fields_to_check:
                    field_obj = bio.get(f) or exp.get(f) or dis.get(f)
                    if field_obj:
                        val, source, conf, ont = extract_field_value(field_obj)
                        if source == "llm_enrichment":
                            all_enriched.append(EnrichedField(f, val, conf or 0.0, sample_id, study_id, s_ctx))
            else:
                s_ctx["title"] = sample_data.get("sample_title", {}).get("value") if isinstance(sample_data.get("sample_title"), dict) else sample_data.get("sample_title", "")
                s_ctx["desc"] = sample_data.get("sample_description", {}).get("value") if isinstance(sample_data.get("sample_description"), dict) else sample_data.get("sample_description", "")
                
                for f in fields_to_check:
                    field_obj = sample_data.get(f)
                    if field_obj:
                        val, source, conf, ont = extract_field_value(field_obj)
                        if source == "llm_enrichment":
                            all_enriched.append(EnrichedField(f, val, conf or 0.0, sample_id, study_id, s_ctx))

    if not all_enriched:
        print("No enriched fields found to review.")
        return

    to_review = random.sample(all_enriched, min(num_to_review, len(all_enriched)))
    
    print("\n" + "=" * 80)
    print(f"{'ENRICHMENT QUALITY REVIEW (Sample of ' + str(len(to_review)) + ')':^80}")
    print("=" * 80)

    for i, f in enumerate(to_review, 1):
        print(f"\n[{i}/{len(to_review)}] {f.field_name.upper()}: {f.value}")
        print(f"  Study: {f.study_id} | Sample: {f.sample_id} | Confidence: {f.confidence:.2f}")
        
        found = "YES" if f.appears_in_source() else "NO"
        print(f"  Appears in source text: {found}")
        
        # Show a snippet of where it might be
        all_text = " ".join(f.context.values())
        if f.value and str(f.value).lower() in all_text.lower():
            idx = all_text.lower().index(str(f.value).lower())
            start = max(0, idx - 60)
            end = min(len(all_text), idx + len(str(f.value)) + 60)
            snippet = all_text[start:end].replace("\n", " ")
            print(f"  Context: ...{snippet}...")
        else:
            print(f"  Context: [Value not found in direct text - inferred from broader study description]")
        print("─" * 40)
