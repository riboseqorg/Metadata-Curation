from __future__ import annotations
from typing import Dict, Any, Optional
from pathlib import Path
import json

from ..io.format_adapter import load_bundle_from_file
from ..ontologies.mapper import normalize_tissue, normalize_cell_type, normalize_organism
from ..schemas.base import BaseProvenance

FIELDS = (
    ("organism", normalize_organism),
    ("tissue", normalize_tissue),
    ("cell_type", normalize_cell_type),
)


def _apply_norm(value: Optional[str], mapper) -> tuple[Optional[str], Optional[str], float]:
    if not value:
        return None, None, 0.0
    v, term, conf = mapper(value)
    return v, term, conf


def re_normalize_file(input_path: Path, out_path: Path, schema_version: Optional[str] = None, strict: bool = False, ontology_cache: Optional[Path] = None) -> Dict[str, Any]:
    bundle = load_bundle_from_file(input_path)

    total = {"samples": len(bundle["samples"]), "updated": 0}

    for sid, sample in bundle["samples"].items():
        updated = False
        for field, mapper in FIELDS:
            prov: Optional[BaseProvenance] = getattr(sample, field, None)
            if prov and prov.value:
                new_val, term, nconf = _apply_norm(prov.value, mapper)
                if new_val:
                    prov.value = new_val
                    prov.ontology_term = term
                    # Keep existing field_confidence if present, recompute combined confidence when normalization_confidence available
                    prov.normalization_confidence = nconf
                    if prov.field_confidence is not None:
                        prov.confidence = float(prov.field_confidence) * float(nconf)
                    elif prov.confidence and nconf:
                        # best effort
                        prov.confidence = float(min(1.0, prov.confidence * nconf))
                    updated = True
                elif strict:
                    setattr(sample, field, None)
                    updated = True
        if updated:
            total["updated"] += 1

    # Write new traceable bundle (study unchanged)
    from ..output.traceable_format import create_traceable_report
    out = create_traceable_report(bundle["study"], bundle["samples"], bundle["runs"], include_raw=True, include_statistics=True)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2, default=lambda o: getattr(o, "model_dump", lambda: str(o))())

    return total

