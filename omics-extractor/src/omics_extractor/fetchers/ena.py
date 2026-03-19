"""ENA Portal API helpers for run discovery.

Docs: https://www.ebi.ac.uk/ena/portal/api/
"""

from __future__ import annotations

from typing import List, Dict, Optional
import requests


def _build_ena_query(organism: str, terms: Optional[List[str]] = None, boolean_query: Optional[str] = None) -> str:
    # Constrain by organism scientific name
    base = f'scientific_name="{organism}"'
    # Prefer boolean query if provided; map to library_strategy terms where possible by simple heuristic
    if boolean_query:
        # Try to map typical strategy phrases to library_strategy field
        # Fall back to raw text match via library_strategy OR study_title OR experiment_title
        mapped = (
            boolean_query
            .replace('"', '\\"')
            .replace('library_strategy', 'library_strategy')
        )
        return f"{base} AND ({mapped})"
    if terms:
        ors = []
        for t in terms:
            t = t.replace('"', '\\"')
            ors.append(f'library_strategy="{t}"')
        return f"{base} AND (" + " OR ".join(ors) + ")"
    return base


def search_ena_runs(organism: str, terms: Optional[List[str]] = None, boolean_query: Optional[str] = None, limit: int = 20000) -> List[Dict[str, str]]:
    q = _build_ena_query(organism, terms=terms, boolean_query=boolean_query)
    url = "https://www.ebi.ac.uk/ena/portal/api/search"
    fields = [
        "run_accession",
        "study_accession",
        "bioproject_accession",
        "library_strategy",
        "scientific_name",
    ]
    params = {
        "result": "read_run",
        "query": q,
        "format": "tsv",
        "fields": ",".join(fields),
        "limit": str(limit),
    }
    try:
        r = requests.get(url, params=params, timeout=60)
        r.raise_for_status()
        text = r.text
        lines = [l for l in text.splitlines() if l.strip()]
        if not lines:
            return []
        headers = lines[0].split('\t')
        out: List[Dict[str, str]] = []
        for line in lines[1:]:
            parts = line.split('\t')
            row = {headers[i]: (parts[i] if i < len(parts) else "") for i in range(len(headers))}
            out.append(row)
        return out
    except Exception:
        return []

