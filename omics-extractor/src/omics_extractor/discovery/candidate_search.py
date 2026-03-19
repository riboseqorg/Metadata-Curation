from __future__ import annotations

"""Cross-source candidate discovery for Ribo‑seq studies.

Strategy:
- Broad organism filter across BioProject and GEO
- Heuristic text matching on titles/abstracts/overall design for Ribo‑seq terms
- Collect candidate BioProject IDs
- For each candidate, fetch runs via SRA and keep Run→study_accession rows
- Provide an evidence score and debugging columns when requested
"""

from typing import Dict, List, Tuple, Iterable, Optional
import re
from dataclasses import dataclass, asdict

from ..fetchers.bioproject import fetch_bioproject_metadata
from ..fetchers.geo import fetch_project_metadata as fetch_geo_project, search_geo_series
from ..fetchers.sra import fetch_project_runs
from ..fetchers.search import search_sra_runinfo
from ..fetchers.pubmed import search_pubmed_for_project, fetch_publication_metadata
from ..fetchers.entrez_config import configure as _configure_entrez, rate_limit, with_retries
from Bio import Entrez
from ..fetchers.ena import search_ena_runs


RIBO_TERMS = [
    r"ribo[-\s]?seq",
    r"ribosome\s+profil(ing|e)",
    r"ribosome\s+footprint(ing)?",
    r"RPF\b",
]


def _score_text(*texts: Iterable[str]) -> Tuple[int, List[str]]:
    """Return (score, hits) for Ribo‑seq terms in provided texts."""
    t = "\n".join([x for x in texts if x])
    hits = []
    score = 0
    for pat in RIBO_TERMS:
        if re.search(pat, t, flags=re.IGNORECASE):
            hits.append(pat)
            score += 1
    return score, hits


@dataclass
class Candidate:
    source: str  # bioproject|geo
    id: str      # PRJNA…
    organism: str | None
    title: str | None
    description: str | None
    pubmed: str | None
    score: int
    hits: List[str]

    def to_row(self) -> Dict[str, str]:
        d = asdict(self)
        d["hits"] = ";".join(self.hits)
        return d


def _bioproject_ids_from_query(organism: str, extra_terms: Optional[List[str]] = None, retmax: int = 10000) -> List[str]:
    """Search BioProject DB for organism + ribo terms and return PRJ accessions.

    Uses esearch to get UIDs, then efetch XML per UID and extracts the accession
    robustly (from ArchiveID or by regex fallback).
    """
    _configure_entrez()
    terms = extra_terms or ["ribo-seq", "ribosome profiling", "ribosome footprinting", "RPF"]
    term_q = " OR ".join([f'"{t}"' if ' ' in t else t for t in terms])
    query = f"{organism}[All Fields] AND ({term_q})"
    res = with_retries(lambda: Entrez.esearch(db="bioproject", term=query, retmax=retmax), read=True)
    uids = res.get("IdList", [])
    if not uids:
        return []
    prjs: List[str] = []
    import xml.etree.ElementTree as ET
    import re as _re
    for uid in uids[:retmax]:
        try:
            xml_txt = with_retries(lambda: Entrez.efetch(db="bioproject", id=uid, rettype="xml"), read=False)
            if isinstance(xml_txt, bytes):
                xml_txt = xml_txt.decode("utf-8", errors="replace")
            acc = None
            try:
                root = ET.fromstring(xml_txt)
                # Try ArchiveID node attribute 'accession'
                for elem in root.iter():
                    if elem.tag.endswith("ArchiveID"):
                        acc = elem.attrib.get("accession") or elem.text
                        if acc:
                            break
                # Fallback: regex for PRJ accessions in the XML text
                if not acc:
                    m = _re.search(r"(PRJ[A-Z]{1,2}\d+)", xml_txt)
                    if m:
                        acc = m.group(1)
            except Exception:
                # If XML parse fails, fallback to regex search
                m = _re.search(r"(PRJ[A-Z]{1,2}\d+)", xml_txt)
                acc = m.group(1) if m else None
            if acc and acc.startswith("PRJ"):
                prjs.append(acc)
        except Exception:
            continue
    # unique & stable order
    seen = set(); out = []
    for p in prjs:
        if p in seen: 
            continue
        seen.add(p); out.append(p)
    return out


def discover_candidates(
    organism: str,
    geo_ids: List[str] | None = None,
    include_bioproject_search: bool = True,
    include_pubmed_evidence: bool = True,
    bioproject_retmax: int = 500,
    extra_terms: Optional[List[str]] = None,
) -> List[Candidate]:
    """Find candidate projects from BioProject and GEO for an organism.

    For now we rely on user-provided GEO IDs if any, otherwise BioProject only.
    This can be extended with Entrez organism searches for GEO in future.
    """
    cands: List[Candidate] = []

    # BioProject pass: simple heuristic — search a coarse list via Entrez utilities
    # We do not have a direct Entrez search here; callers can enumerate PRJ ids.
    # In this repo, builder fetches by explicit PRJ id. So we provide a helper here
    # that can be fed with known PRJs (e.g., from an upstream list) later.
    # For minimal viable discovery, we simply skip automatic enumeration and rely
    # on GEO-derived BioProject IDs + optional manual PRJs.

    # GEO: find series and convert to BioProject if possible; or accept provided seeds
    if geo_ids is None:
        # try a light GEO search for the organism and terms
        try:
            found = search_geo_series(organism, terms=extra_terms)
        except Exception:
            found = []
        geo_ids = found
    if geo_ids:
        for gse in geo_ids:
            try:
                gp = fetch_geo_project(gse)
            except Exception:
                continue
            if gp.organism and organism.lower() not in gp.organism.lower():
                continue
            score, hits = _score_text(gp.title, gp.summary, gp.overall_design or "")
            prj = gp.bioproject_id
            if not prj:
                # keep GEO-only candidate for inspection but it won’t yield runs until mapped
                continue
            cands.append(Candidate(
                source="geo",
                id=prj,
                organism=gp.organism,
                title=gp.title,
                description=gp.summary,
                pubmed=gp.pubmed_id,
                score=score,
                hits=hits,
            ))

    # BioProject search pass
    if include_bioproject_search:
        for prj in _bioproject_ids_from_query(organism, extra_terms=extra_terms, retmax=bioproject_retmax):
            try:
                bp = fetch_bioproject_metadata(prj)
            except Exception:
                continue
            if bp.organism and organism.lower() not in (bp.organism or "").lower():
                continue
            score, hits = _score_text(bp.title or "", bp.description or "")
            cands.append(Candidate(
                source="bioproject",
                id=prj,
                organism=bp.organism,
                title=bp.title,
                description=bp.description,
                pubmed=bp.publication_id,
                score=score,
                hits=hits,
            ))

    # De-duplicate by id
    uniq: Dict[str, Candidate] = {}
    for c in cands:
        if (prev := uniq.get(c.id)) is None or c.score > prev.score:
            uniq[c.id] = c
    out = list(uniq.values())

    # Optional PubMed evidence augmentation
    if include_pubmed_evidence:
        for c in out:
            try:
                pmids = search_pubmed_for_project(c.id)
            except Exception:
                pmids = []
            if not pmids:
                continue
            # Take first PMID for light scoring
            try:
                pm = fetch_publication_metadata(pmids[0])
                add, hits = _score_text(pm.title, pm.abstract)
                c.score += add
                c.hits.extend(hits)
                if not c.pubmed:
                    c.pubmed = pm.pmid
            except Exception:
                pass

    return out


def runs_from_candidates(cands: List[Candidate], include_ena: bool = False, terms: Optional[List[str]] = None, boolean_query: Optional[str] = None):
    """Fetch runs per BioProject candidate and return (rows_by_prj, prj_sources).

    rows_by_prj: Dict[PRJ -> List[ {Run,study_accession} ]]
    prj_sources: Dict[PRJ -> Set[str]] where sources include 'sra' and optionally 'ena'.
    """
    rows_by_prj: Dict[str, List[Dict[str, str]]] = {}
    prj_sources: Dict[str, set] = {}
    for c in cands:
        if not c.id.startswith("PRJ"):
            continue
        # Pull full runinfo to capture SRAStudy if present
        _, run_rows = search_sra_runinfo(f"{c.id}[BioProject]")
        if not run_rows:
            # fallback to minimal list
            s_rows = []
            for rid in fetch_project_runs(c.id):
                s_rows.append({"Run": rid, "study_accession": c.id})
            if s_rows:
                rows_by_prj.setdefault(c.id, []).extend(s_rows)
                prj_sources.setdefault(c.id, set()).add('sra')
            continue
        s_rows = []
        for r in run_rows:
            run = (r.get("Run") or r.get("run") or "").strip()
            study = (r.get("SRAStudy") or r.get("Study") or r.get("study_accession") or "").strip()
            study = study or c.id
            if run:
                s_rows.append({"Run": run, "study_accession": study})
        if s_rows:
            rows_by_prj.setdefault(c.id, []).extend(s_rows)
            prj_sources.setdefault(c.id, set()).add('sra')
        if include_ena:
            # Pull ENA as well and merge
            ena_rows = search_ena_runs(c.organism or "", terms=terms, boolean_query=boolean_query)
            e_rows = []
            for er in ena_rows:
                rid = (er.get('run_accession') or '').strip()
                # keep only this PRJ
                prj = (er.get('bioproject_accession') or '').strip()
                if prj != c.id:
                    continue
                st = (er.get('study_accession') or prj or '').strip()
                if rid and st:
                    e_rows.append({"Run": rid, "study_accession": st})
            if e_rows:
                rows_by_prj.setdefault(c.id, []).extend(e_rows)
                prj_sources.setdefault(c.id, set()).add('ena')
    # de-dup runs per PRJ
    for prj, lst in list(rows_by_prj.items()):
        seen = set(); uniq = []
        for r in lst:
            if r['Run'] in seen: continue
            seen.add(r['Run']); uniq.append(r)
        rows_by_prj[prj] = uniq
    return rows_by_prj, prj_sources
