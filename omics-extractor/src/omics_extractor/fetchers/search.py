"""Lightweight SRA search helpers (run discovery).

Builds SRA E-utilities queries and returns runinfo rows, reusing the
central Entrez configuration and polite rate limiting.
"""

from __future__ import annotations

from typing import List, Dict, Tuple, Optional
from Bio import Entrez
from .entrez_config import configure as _configure_entrez, rate_limit, with_retries
import io
import csv

_configure_entrez()


DEFAULT_RIBO_TERMS = [
    '"ribosome profiling"',
    '"Ribo-Seq"',
    'ribo-seq',
    'RiboSeq',
    'RPF',
    '"ribosome footprinting"',
]


def build_query(organism: str, extra_terms: List[str] | None = None) -> str:
    """Build an SRA boolean query for an organism and (optional) terms.

    If extra_terms is provided, they are OR-joined and ANDed with the organism.
    """
    org_q = f"({organism})[Organism]"
    if extra_terms:
        term_q = " OR ".join(extra_terms)
        return f"{org_q} AND ({term_q})"
    return org_q


def search_sra_runinfo(query: str, retmax: int = 200000) -> Tuple[str, List[Dict[str, str]]]:
    """Execute an SRA search and return raw runinfo CSV plus parsed rows.

    Returns a tuple of (raw_csv_text, list_of_row_dicts).
    """
    results = with_retries(lambda: Entrez.esearch(db="sra", term=query, retmax=retmax, usehistory="y"), read=True)

    count = int(results.get("Count", 0))
    if count == 0:
        return ("", [])

    webenv = results["WebEnv"]
    query_key = results["QueryKey"]

    text = with_retries(lambda: Entrez.efetch(
        db="sra",
        query_key=query_key,
        WebEnv=webenv,
        rettype="runinfo",
        retmode="text",
    ), read=False)

    if isinstance(text, bytes):
        text = text.decode("utf-8", errors="replace")

    reader = csv.DictReader(io.StringIO(text))
    rows = [row for row in reader]
    return (text, rows)


def search_riboseq_runs_for_organism(organism: str) -> Tuple[str, List[Dict[str, str]]]:
    """Convenience: Ribo‑seq term bundle for a given organism."""
    q = build_query(organism, DEFAULT_RIBO_TERMS)
    return search_sra_runinfo(q)


def build_terms_for_strategy(strategy: str) -> List[str]:
    """Heuristic term list for a given library strategy.

    - For Ribo-Seq, reuse the curated defaults.
    - Otherwise, include the quoted strategy and a lowercase hyphenated variant.
    """
    s = (strategy or '').strip()
    if not s:
        return []
    low = s.lower()
    if low.startswith('ribo') or 'ribosome' in low:
        return DEFAULT_RIBO_TERMS
    terms = [f'"{s}"']
    terms.append(low.replace(' ', '-').replace('_', '-'))
    return list(dict.fromkeys(terms))


def search_runs_for_organism(organism: str, terms: List[str] | None, retmax: int = 200000) -> Tuple[str, List[Dict[str, str]]]:
    """Generic SRA search for organism with optional terms."""
    q = build_query(organism, terms)
    return search_sra_runinfo(q, retmax=retmax)


def search_runs_with_query(organism: str, boolean_query: str, retmax: int = 200000) -> Tuple[str, List[Dict[str, str]]]:
    """Search SRA runs with a boolean query string (one call).

    boolean_query should be an SRA-compatible boolean expression.
    """
    q = f"({organism})[Organism] AND ({boolean_query})"
    return search_sra_runinfo(q, retmax=retmax)


def fetch_runs_for_bioprojects(prjs: List[str], chunk_size: int = 50, retmax: int = 200000) -> List[Dict[str, str]]:
    """Fetch SRA runinfo rows for a list of PRJ accessions using chunked OR queries."""
    rows: List[Dict[str, str]] = []
    if not prjs:
        return rows
    # Build OR-chunks of PRJ[BioProject]
    for i in range(0, len(prjs), chunk_size):
        chunk = prjs[i:i+chunk_size]
        term = " OR ".join([f"{p}[BioProject]" for p in chunk])
        _, r = search_sra_runinfo(term, retmax=retmax)
        rows.extend(r)
    # de-dup by Run
    seen = set(); uniq = []
    for r in rows:
        rid = (r.get("Run") or r.get("run") or "").strip()
        if not rid or rid in seen:
            continue
        seen.add(rid); uniq.append(r)
    return uniq

def bioproject_ids_via_esummary(
    organism: str,
    terms: Optional[list[str]] | None = None,
    boolean_query: str | None = None,
    retmax: int = 1000,
) -> list[str]:
    """Enumerate BioProject accessions via esearch+esummary (chunked, fast)."""
    _configure_entrez()
    from Bio import Entrez
    # Build query
    q_org = f"{organism}[All Fields]"
    if boolean_query:
        q = f"{q_org} AND ({boolean_query})"
    elif terms:
        ors = " OR ".join([f'"{t}"' if ' ' in t else t for t in terms])
        q = f"{q_org} AND ({ors})"
    else:
        q = q_org
    res = Entrez.read(Entrez.esearch(db='bioproject', term=q, retmax=retmax))
    ids = res.get('IdList', [])
    if not ids:
        return []
    prjs: list[str] = []
    for i in range(0, len(ids), 400):
        chunk = ids[i:i+400]
        summ = Entrez.read(Entrez.esummary(db='bioproject', id=','.join(chunk)))
        # Biopython may return a dict with DocumentSummarySet, or a list
        docs = []
        if isinstance(summ, dict) and 'DocumentSummarySet' in summ:
            docs = summ['DocumentSummarySet'].get('DocumentSummary', [])
        elif isinstance(summ, list):
            docs = summ
        for doc in docs:
            if isinstance(doc, dict):
                acc = (doc.get('Project_Acc') or doc.get('Project_Accession') or doc.get('Accession') or doc.get('ProjectAcc'))
                if acc and str(acc).startswith('PRJ'):
                    prjs.append(str(acc))
    # unique stable
    seen=set(); out=[]
    for x in prjs:
        if x in seen:
            continue
        seen.add(x); out.append(x)
    return out
