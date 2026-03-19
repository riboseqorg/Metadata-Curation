from __future__ import annotations
from typing import Dict, Any, List, Optional, Tuple

from ...schemas.base import SampleMetadata
from ...ontologies.mapper import get_mapper
from ...extraction.llm_extractor import load_extraction_scheme, build_extraction_prompt

FEW_SHOT: Dict[str, List[Tuple[str,str]]] = {
    "ribo_seq": [
        (
            "Study: ribosome profiling in mouse liver with cycloheximide pretreatment; Sample: liver tissue;",
            '{"tissue": "liver", "inhibitor": "cycloheximide", "library_strategy": "ribo-seq", "confidence": {"tissue": 0.9, "inhibitor": 0.8, "library_strategy": 0.9}}'
        ),
        (
            "Study: Ribo-seq in human HEK293 cells using RNase I digestion; Sample: HEK293 cell line;",
            '{"cell_line": "hek293", "nuclease": "rnase i", "library_strategy": "ribo-seq", "confidence": {"cell_line": 0.9, "nuclease": 0.8, "library_strategy": 0.9}}'
        ),
        (
            "Study: Yeast ribosome profiling with micrococcal nuclease; Sample: whole cells;",
            '{"organism": "saccharomyces cerevisiae", "nuclease": "micrococcal nuclease", "confidence": {"organism": 0.9, "nuclease": 0.8}}'
        ),
        (
            "Study: Ribo-seq in mouse hepatocytes (25°C, 45 min digestion); Sample: liver tissue;",
            '{"tissue": "liver", "digestion_temperature": "25 c", "digestion_time": "45 min", "confidence": {"tissue": 0.9, "digestion_temperature": 0.7, "digestion_time": 0.7}}'
        ),
    ],
    "nanopore_rna": [
        ("Direct RNA-seq on human HEK293 cell line;", '{"cell_line": "hek293", "cell_type": null, "confidence": {"cell_line": 0.9}}')
    ],
}

SHORTLIST_FIELDS = ("tissue", "cell_type", "treatment")

def infer_strategy(sample: SampleMetadata) -> Optional[str]:
    tc = sample.technical_context or {}
    ls = (tc.get("library_strategy") or "").lower()
    if any(k in ls for k in ("ribo", "ribosome")):
        return "ribo_seq"
    if "nanopore" in (tc.get("platform") or "").lower():
        return "nanopore_rna"
    return None


def build_prompt(
    study_title: str,
    study_description: str,
    sample: SampleMetadata,
    abstract: Optional[str],
    journal: Optional[str],
    authors: Optional[str],
    publication_date: Optional[str],
    scheme: str = "default",
    shortlist_k: int = 0,
    few_shot: str = "off",
    extra_context: Optional[Dict[str, Any]] = None,
) -> str:
    fields = load_extraction_scheme(scheme)

    base = build_extraction_prompt(
        study_title=study_title,
        study_description=study_description,
        sample_title=sample.sample_title.value if sample.sample_title else None,
        sample_description=sample.sample_description.value if sample.sample_description else None,
        characteristics=sample.raw_characteristics,
        run_metadata=sample.technical_context,
        abstract=abstract,
        journal=journal,
        authors=authors,
        publication_date=publication_date,
        existing_metadata={k: getattr(getattr(sample, k), 'value', None) for k in fields.keys() if getattr(sample, k, None)},
        fields=fields,
        extra_context=extra_context,
    )

    extra: List[str] = []

    if shortlist_k and shortlist_k > 0:
        mapper = get_mapper()
        sl_texts: List[str] = []
        for fname in SHORTLIST_FIELDS:
            cand: List[tuple] = []
            raw = (sample.raw_characteristics or {}).get(fname) or (sample.raw_characteristics or {}).get(fname.replace('_',' '))
            if raw and isinstance(raw, str):
                onto = 'uberon' if fname == 'tissue' else ('cl' if fname == 'cell_type' else 'efo')
                try:
                    # Fetch up to shortlist_k candidates when available
                    res = mapper.ols_client.search(raw, onto, rows=shortlist_k)
                    if isinstance(res, list):
                        cand.extend(res[:shortlist_k])
                    elif res:
                        cand.append(res)
                except Exception:
                    pass
            seen = set()
            picked: List[tuple] = []
            for label, obo_id, _conf in cand:
                key = (label.lower().strip(), obo_id)
                if key in seen:
                    continue
                seen.add(key)
                picked.append((label, obo_id))
                if len(picked) >= shortlist_k:
                    break
            if picked:
                sl_texts.append(f"- {fname}: " + ", ".join([f"{lbl} ({oid})" for lbl, oid in picked]))
        if sl_texts:
            extra.append("\nALLOWED VALUES (ontology shortlists; choose one or null):\n" + "\n".join(sl_texts) + "\n")

    if few_shot and few_shot != "off":
        strat = infer_strategy(sample) or (scheme if scheme != "auto" else None)
        shots = FEW_SHOT.get(strat or "", [])
        if shots:
            ex_lines: List[str] = []
            for prompt, js in shots:
                ex_lines.append("Example:\n" + prompt + "\nOutput JSON:\n" + js)
            extra.append("\nFEW-SHOT EXAMPLES:\n" + "\n\n".join(ex_lines) + "\n")

    return base + ("\n".join(extra) if extra else "")
