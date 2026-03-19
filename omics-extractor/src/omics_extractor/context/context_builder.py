from __future__ import annotations
from typing import Dict, Any, Optional
from dataclasses import dataclass
import hashlib

from ..schemas.base import StudyMetadata, SampleMetadata

@dataclass
class ContextOpts:
    include_geo_overall: bool = True
    include_geo_protocols: bool = True
    include_geo_source: bool = True
    merge_channels: bool = True
    include_bioproject_type: bool = True
    include_bioproject_scope: bool = True
    include_pubmed_mesh: bool = True
    include_biosample_package: bool = False
    include_biosample_description: bool = True
    truncate_chars: int = 0  # 0 = no truncation


def _truncate(txt: Optional[str], n: int) -> Optional[str]:
    if not txt or n <= 0:
        return txt
    return txt if len(txt) <= n else txt[:n]


def _merge_protocol_block(protocols: Dict[str, str], kind: str, opts: ContextOpts) -> Optional[str]:
    ch1 = protocols.get(f"{kind}_ch1")
    ch2 = protocols.get(f"{kind}_ch2")
    if not ch1 and not ch2:
        return None
    if not opts.merge_channels:
        parts = []
        if ch1:
            parts.append(f"CH1: {ch1}")
        if ch2:
            parts.append(f"CH2: {ch2}")
        return " \n".join(parts)
    parts = []
    if ch1:
        parts.append(f"CH1: {ch1}")
    if ch2:
        parts.append(f"CH2: {ch2}")
    return " ".join(parts)


def compose_llm_context(
    study: StudyMetadata,
    sample: SampleMetadata,
    profile: str = "extended",
    opts: Optional[ContextOpts] = None,
) -> Dict[str, Any]:
    opts = opts or ContextOpts()

    ctx: Dict[str, Any] = {
        "study_title": getattr(study.title, 'value', None),
        "study_description": getattr(study.description, 'value', None),
        "abstract": getattr(study, 'paper_abstract', None).value if getattr(study, 'paper_abstract', None) else None,
        "journal": getattr(study, 'journal', None).value if getattr(study, 'journal', None) else None,
        "authors": getattr(study, 'authors', None),
        "publication_date": getattr(study, 'publication_date', None),
        "sample_title": getattr(sample, 'sample_title', None).value if getattr(sample, 'sample_title', None) else None,
        "sample_description": getattr(sample, 'sample_description', None).value if getattr(sample, 'sample_description', None) else None,
        "raw_characteristics": sample.raw_characteristics,
        "technical_context": sample.technical_context,
    }

    if profile in ("extended", "full"):
        if opts.include_geo_overall:
            overall = None
            cf = getattr(study, 'custom_fields', None)
            if isinstance(cf, dict):
                overall = cf.get('geo_overall_design')
            ctx["geo_overall_design"] = _truncate(overall, opts.truncate_chars) if overall else overall
        if opts.include_geo_source:
            src_name = None
            if sample.raw_characteristics and isinstance(sample.raw_characteristics, dict):
                src_name = sample.raw_characteristics.get('source_name') or sample.raw_characteristics.get('source')
            ctx["geo_source_name"] = _truncate(src_name, opts.truncate_chars) if src_name else src_name
        if opts.include_geo_protocols:
            # Prefer explicit attribute if present; otherwise fall back to custom_fields["protocols"].
            protocols = getattr(sample, 'protocols', None) or (
                (getattr(sample, 'custom_fields', None) or {}).get('protocols') if getattr(sample, 'custom_fields', None) else {}
            ) or {}
            block: Dict[str, str] = {}
            for kind in ("growth", "extract", "treatment"):
                merged = _merge_protocol_block(protocols, kind, opts)
                if merged:
                    block[kind] = _truncate(merged, opts.truncate_chars)
            if block:
                ctx["geo_protocols"] = block
        # BioProject type/scope
        if opts.include_bioproject_type or opts.include_bioproject_scope:
            cf = getattr(study, "custom_fields", None) or {}
            bt = cf.get("bioproject_data_type") if opts.include_bioproject_type else None
            bs = cf.get("bioproject_scope") if opts.include_bioproject_scope else None
            if bt or bs:
                if bt and bs:
                    ctx["bioproject_summary"] = f"Type: {bt}; Scope: {bs}"
                elif bt:
                    ctx["bioproject_summary"] = f"Type: {bt}"
                else:
                    ctx["bioproject_summary"] = f"Scope: {bs}"
        # PubMed MeSH
        if opts.include_pubmed_mesh:
            cf = getattr(study, "custom_fields", None) or {}
            mesh = cf.get("pubmed_mesh")
            if mesh:
                ctx["pubmed_mesh"] = ", ".join(mesh) if isinstance(mesh, list) else str(mesh)
        # BioSample package/description
        if opts.include_biosample_package:
            cf_s = getattr(sample, "custom_fields", None) or {}
            pkg = cf_s.get("biosample_package")
            if pkg:
                ctx["biosample_package"] = pkg
        if opts.include_biosample_description:
            try:
                desc = getattr(sample, "sample_description", None)
                if desc and getattr(desc, "source", None) == "biosample":
                    ctx["biosample_description"] = desc.value
            except Exception:
                pass


    snapshot_blocks: Dict[str, Any] = {}
    for k in ("geo_overall_design", "geo_source_name"):
        if ctx.get(k):
            txt = ctx[k]
            snapshot_blocks[k] = {
                "length": len(txt),
                "sha256": hashlib.sha256(txt.encode('utf-8')).hexdigest(),
            }
    if ctx.get("geo_protocols"):
        for pk, pv in ctx["geo_protocols"].items():
            snapshot_blocks[f"geo_protocols.{pk}"] = {
                "length": len(pv),
                "sha256": hashlib.sha256(pv.encode('utf-8')).hexdigest(),
            }

    ctx["context_snapshot"] = {
        "profile": profile,
        "opts": {
            "include_geo_overall": opts.include_geo_overall,
            "include_geo_protocols": opts.include_geo_protocols,
            "include_geo_source": opts.include_geo_source,
            "merge_channels": opts.merge_channels,
            "truncate_chars": opts.truncate_chars,
        },
        "blocks": snapshot_blocks,
    }

    return ctx
