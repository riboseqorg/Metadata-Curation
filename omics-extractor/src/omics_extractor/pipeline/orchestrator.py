from __future__ import annotations
from typing import Any, Dict
from pathlib import Path
import json

from ..extraction.builder import build_project_metadata
from ..output.traceable_format import create_traceable_report
# Compatible import whether orchestrator is imported as a package or module
try:
    from ..normalize.service import re_normalize_file  # when package root is omics_extractor
except ImportError:
    from omics_extractor.normalize.service import re_normalize_file  # fallback absolute
from ..io.format_adapter import load_bundle_from_file
from ..output.reports import generate_tabular_report


def run_extract(bioproject_id: str, out_path: Path) -> Dict[str, Any]:
    bundle = build_project_metadata(bioproject_id)
    report = create_traceable_report(
        study=bundle["study"],
        samples=bundle["samples"],
        runs=bundle["runs"],
        include_raw=True,
        include_statistics=True,
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(report, f, indent=2, default=lambda o: getattr(o, "model_dump", lambda: str(o))())
    return {"study": 1, "samples": len(bundle["samples"]), "runs": len(bundle["runs"]) }


def run_normalize(input_path: Path, out_path: Path, schema_version: str | None, strict: bool, ontology_cache: Path | None) -> Dict[str, Any]:
    stats = re_normalize_file(input_path, out_path, schema_version=schema_version, strict=strict, ontology_cache=ontology_cache)
    return stats


def run_export(inputs: list[Path], out_path: Path, fmt: str) -> bool:
    return generate_tabular_report(inputs, out_path, format_type=fmt)
