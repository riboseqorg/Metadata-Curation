from __future__ import annotations
from pathlib import Path
from typing import Dict, Any, List
import json

from fastapi import FastAPI, HTTPException, Request
from fastapi.responses import HTMLResponse, JSONResponse


def _load_records(p: Path) -> List[Dict[str, Any]]:
    if not p.exists():
        raise HTTPException(status_code=404, detail="input not found")
    if p.suffix.lower() == ".jsonl":
        return [json.loads(line) for line in p.read_text().splitlines() if line.strip()]
    # Minimal TSV/CSV support: convert to JSON records with llm_context limited
    import csv
    delim = "\t" if p.suffix.lower() == ".tsv" else ","
    rows = list(csv.DictReader(p.open(), delimiter=delim))
    recs = []
    for r in rows:
        recs.append({
            "ids": {"bioproject_id": r.get("bioproject_id"), "sample_id": r.get("sample_id")},
            "llm_context": {
                "study_title": r.get("study_title"),
                "study_description": r.get("study_description"),
                "raw_characteristics": r.get("raw_characteristics_json"),
                "technical_context": r.get("technical_context_json"),
            },
            "targets": {},
            "current_values": json.loads(r.get("current_values_json")) if r.get("current_values_json") else {},
        })
    return recs


def build_app(input_path: Path, out_path: Path, read_only: bool = False) -> FastAPI:
    app = FastAPI()
    db: List[Dict[str, Any]] = _load_records(input_path)
    # Optional curated preload
    curated: Dict[tuple, Dict[str, Any]] = {}
    if out_path.exists():
        try:
            for line in out_path.read_text().splitlines():
                if not line.strip():
                    continue
                rec = json.loads(line)
                key = (rec.get('ids',{}).get('bioproject_id'), rec.get('ids',{}).get('sample_id'))
                curated[key] = rec.get('curated', {})
        except Exception:
            curated = {}

    TEMPLATE = """
<!doctype html>
<meta charset="utf-8" />
<title>Gold Curation</title>
<link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/@picocss/pico@2/css/pico.min.css">
<script src="https://unpkg.com/htmx.org@1.9.10"></script>
<main class="container">
  <hgroup>
    <h1>Gold Curation</h1>
    <p id="status">Samples: {{N}}</p>
  </hgroup>
  <div class="grid">
    <div>
      <details open>
        <summary>Targets</summary>
        <div id="targets" hx-get="/targets" hx-trigger="load"></div>
      </details>
      <button id="prev" hx-get="/prev" hx-target="#page" hx-swap="outerHTML">Prev (j)</button>
      <button id="next" hx-get="/next" hx-target="#page" hx-swap="outerHTML">Next (k)</button>
      <small>Use j/k to navigate</small>
    </div>
    <article id="page" hx-get="/sample/0" hx-trigger="load"></article>
  </div>
</main>
<script>
  document.addEventListener('keydown', e => { if (e.key==='j') document.getElementById('prev').click(); if (e.key==='k') document.getElementById('next').click(); });
</script>
"""

    state = {"idx": 0}

    def _sample(idx: int) -> Dict[str, Any]:
        if idx < 0: idx = 0
        if idx >= len(db): idx = len(db)-1
        state["idx"] = idx
        return db[idx]

    def _targets_html(r: Dict[str, Any]) -> str:
        # prefer curated values if present; else current_values
        key = (r['ids']['bioproject_id'], r['ids']['sample_id'])
        cur_saved = curated.get(key, {})
        cur = r.get("current_values", {})
        t_fields = [
            "organism", "tissue", "cell_type", "cell_line", "treatment", "disease",
            "genotype", "strain", "age", "sex", "developmental_stage"
        ]
        inputs = []
        for f in t_fields:
            val = cur_saved.get(f) or cur.get(f, "")
            disabled = "disabled" if read_only else ""
            inputs.append(f"<label>{f}<input name='{f}' value='{(val or '')}' {disabled}></label>")
        save_btn = "" if read_only else "<button hx-post='/save' hx-include='#targets' hx-target='#saveStatus'>Save (Enter)</button>"
        return f"<form>{''.join(inputs)}<div id='saveStatus'></div>{save_btn}</form>"

    def _render(idx: int) -> str:
        r = _sample(idx)
        ctx = r.get("llm_context", {})
        # Group important blocks explicitly
        blocks = []
        def block(name, value):
            if value:
                blocks.append(f"<details open><summary>{name}</summary><pre style='white-space:pre-wrap'>{json.dumps(value, indent=2) if isinstance(value, (dict,list)) else (value or '')}</pre></details>")
        block('Study Title', ctx.get('study_title'))
        block('Study Description', ctx.get('study_description'))
        block('Raw Characteristics', ctx.get('raw_characteristics'))
        block('Technical (SRA)', ctx.get('technical_context'))
        block('BioProject Summary', ctx.get('bioproject_summary'))
        block('GEO Overall Design', ctx.get('geo_overall_design'))
        block('GEO Protocols', ctx.get('geo_protocols'))
        block('PubMed MeSH', ctx.get('pubmed_mesh'))
        block('BioSample Description', ctx.get('biosample_description'))
        block('SRA Runs', ctx.get('sra_runs'))
        # Any remaining blocks
        for k,v in ctx.items():
            if k in {'study_title','study_description','raw_characteristics','technical_context','bioproject_summary','geo_overall_design','geo_protocols','pubmed_mesh','biosample_description','sra_runs'}:
                continue
            block(k, v)
        ctx_html = ''.join(blocks)
        # Include an out-of-band swap to refresh the left targets form
        targets_oob = f"<div id='targets' hx-swap-oob='true'>{_targets_html(r)}</div>"
        return f"""
        {targets_oob}
        <article>
          <header><strong>{r['ids']['bioproject_id']} · {r['ids']['sample_id']}</strong></header>
          {ctx_html}
        </article>
        """

    @app.get("/", response_class=HTMLResponse)
    def index():
        return HTMLResponse(TEMPLATE.replace("{{N}}", str(len(db))))

    @app.get("/sample/{idx}", response_class=HTMLResponse)
    def sample(idx: int):
        return HTMLResponse(_render(idx))

    @app.get("/next", response_class=HTMLResponse)
    def nxt():
        return HTMLResponse(_render(state["idx"]+1))

    @app.get("/prev", response_class=HTMLResponse)
    def prv():
        return HTMLResponse(_render(state["idx"]-1))

    @app.get("/targets", response_class=HTMLResponse)
    def targets():
        return HTMLResponse(_targets_html(_sample(state["idx"])))

    @app.post("/save")
    async def save(request: Request):
        data = dict((await request.form()).items()) if hasattr(request, 'form') else {}
        if read_only:
            return HTMLResponse("<small>Read-only</small>")
        # Persist as JSONL (append or update line by line by sample_id)
        r = _sample(state["idx"])
        rec = {"ids": r["ids"], "curated": data}
        out_path.parent.mkdir(parents=True, exist_ok=True)
        with out_path.open("a") as f:
            f.write(json.dumps(rec) + "\n")
        # Update in-memory curated values for immediate reflection
        curated[(r['ids']['bioproject_id'], r['ids']['sample_id'])] = data
        return HTMLResponse("<small>Saved</small>")

    return app
