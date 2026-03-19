from __future__ import annotations
from typing import Dict, Any, List, Tuple
from pathlib import Path
import json

from ..schemas.base import StudyMetadata, SampleMetadata, RunMetadata, BaseProvenance

# Utilities to bridge between on-disk JSON (traceable or flat) and Pydantic models.

PRIORITY_FIELDS = [
    "tissue", "cell_type", "cell_line", "strain", "treatment", "disease",
    "developmental_stage", "age", "sex", "genotype", "condition",
    "timepoint", "replicate", "batch", "stress", "temperature",
    "growth_condition", "sample_title", "sample_description"
]


def _ensure_base_provenance(field_obj: Any, default_source: str = "structured_field") -> Any:
    if not isinstance(field_obj, dict):
        return field_obj
    if "value" in field_obj:
        field_obj.setdefault("source", default_source)
        field_obj.setdefault("source_id", "unknown")
        field_obj.setdefault("confidence", 0.9)
    return field_obj


def parse_traceable_sample(sample_data: Dict[str, Any], study: Dict[str, Any]) -> Tuple[str, SampleMetadata]:
    quick_view = sample_data.get("quick_view", {})
    bio_meta = sample_data.get("biological_metadata", {})

    # bioproject_id fallbacks
    bioproject_id = quick_view.get("bioproject_id") or \
                    sample_data.get("identifiers", {}).get("bioproject_id") or \
                    study.get("identifiers", {}).get("bioproject_id") or \
                    study.get("quick_view", {}).get("bioproject_id") or \
                    study.get("bioproject_id") or \
                    study.get("bioproject", {}).get("value")

    flat_data: Dict[str, Any] = {
        "sample_id": quick_view.get("sample_id") or sample_data.get("sample_id"),
        "bioproject_id": bioproject_id,
        "organism": _ensure_base_provenance(bio_meta.get("organism", {})),
        "raw_characteristics": sample_data.get("raw_biosample_attributes") or sample_data.get("raw_characteristics"),
        "technical_context": sample_data.get("technical_context"),
    }

    for field in PRIORITY_FIELDS:
        if field in bio_meta:
            flat_data[field] = _ensure_base_provenance(bio_meta[field])
        elif field in sample_data and isinstance(sample_data.get(field), dict):
            flat_data[field] = _ensure_base_provenance(sample_data[field])

    sample = SampleMetadata(**flat_data)
    return sample.sample_id, sample


def parse_flat_sample(sample_data: Dict[str, Any], study: Dict[str, Any]) -> Tuple[str, SampleMetadata]:
    # ensure bioproject_id if missing
    if sample_data.get("bioproject_id") is None:
        bioproject_id = study.get("bioproject_id") or study.get("bioproject", {}).get("value")
        if bioproject_id:
            sample_data["bioproject_id"] = bioproject_id
    sample = SampleMetadata(**sample_data)
    return sample.sample_id, sample


def load_bundle_from_file(path: Path) -> Dict[str, Any]:
    """Load a study/samples/runs bundle from a JSON file, accepting traceable or flat formats."""
    with open(path) as f:
        data = json.load(f)

    study_raw = data.get("study", {})
    # Fast path: if already serialized via model_dump, try to reconstruct simply
    if study_raw and isinstance(study_raw, dict) and (study_raw.get("title") or study_raw.get("core_metadata")):
        # Build StudyMetadata with safe defaults
        def to_prov(obj: Dict[str, Any]) -> BaseProvenance:
            return BaseProvenance(
                value=obj.get("value") or obj.get("model_dump", {}).get("value", ""),
                source=obj.get("source", "unknown"),
                source_id=obj.get("source_id", "unknown"),
                confidence=float(obj.get("confidence", 0.0)),
                field_confidence=obj.get("field_confidence"),
                normalization_confidence=obj.get("normalization_confidence"),
                extracted_text=obj.get("extracted_text"),
                ontology_term=obj.get("ontology_term"),
                ontology_label=obj.get("ontology_label"),
                extraction_method=obj.get("extraction_method"),
            )
        title = study_raw.get("quick_view", {}).get("title") or study_raw.get("core_metadata", {}).get("title", {}).get("value") or study_raw.get("title", {}).get("value") or ""
        organism_obj = study_raw.get("organism") or study_raw.get("core_metadata", {}).get("organism")
        organism = to_prov(organism_obj) if isinstance(organism_obj, dict) else BaseProvenance(value=str(organism_obj or "unknown"), source="unknown", source_id="unknown", confidence=0.0)
        description_obj = study_raw.get("description") or study_raw.get("core_metadata", {}).get("description") or {"value": ""}
        study = StudyMetadata(
            bioproject_id=study_raw.get("bioproject_id") or study_raw.get("quick_view", {}).get("bioproject_id") or data.get("bioproject_id", ""),
            title=to_prov(study_raw.get("title", {}) or study_raw.get("core_metadata", {}).get("title", {}) or {"value": title, "confidence": 0.0, "source": "unknown", "source_id": "unknown"}),
            description=to_prov(description_obj),
            organism=organism,
            gse_id=study_raw.get("gse_id"),
            sra_study_id=study_raw.get("sra_study_id"),
            pmid=study_raw.get("pmid"),
            pmc_id=study_raw.get("pmc_id"),
            publication_date=study_raw.get("publication_date"),
        )
    else:
        # Minimal fallback
        study = StudyMetadata(
            bioproject_id=str(data.get("bioproject_id") or ""),
            title=BaseProvenance(value=str(data.get("title") or ""), source="unknown", source_id="unknown", confidence=0.0),
            description=BaseProvenance(value=str(data.get("description") or ""), source="unknown", source_id="unknown", confidence=0.0),
            organism=BaseProvenance(value=str(data.get("organism") or "unknown"), source="unknown", source_id="unknown", confidence=0.0),
        )

    # Samples
    samples_out: Dict[str, SampleMetadata] = {}
    for sid, sdata in (data.get("samples") or {}).items():
        if "quick_view" in sdata or "biological_metadata" in sdata:
            sid2, sample = parse_traceable_sample(sdata, study_raw)
        else:
            sid2, sample = parse_flat_sample(sdata, study_raw)
        samples_out[sid2] = sample

    # Runs: accept dicts already in simple dict form
    runs_out: List[RunMetadata] = []
    for r in data.get("runs", []):
        try:
            # library_strategy can be dict; normalize
            lib = r.get("library_strategy")
            if isinstance(lib, dict):
                r = dict(r)
                r["library_strategy"] = BaseProvenance(
                    value=lib.get("value", "unknown"),
                    source=lib.get("source", "unknown"),
                    source_id=lib.get("source_id", "unknown"),
                    confidence=float(lib.get("confidence", 0.0)),
                )
            runs_out.append(RunMetadata(**r))
        except Exception:
            # Skip invalid entries rather than failing the whole load
            continue

    return {"study": study, "samples": samples_out, "runs": runs_out}

