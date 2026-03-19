"""Command-line interface for omics-extractor."""
import sys
import json
import argparse
import os
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Optional
from .extraction.builder import build_project_metadata
from .extraction.llm_extractor import enrich_sample_metadata
from .extraction.batch_enricher import enrich_batch_from_files
from .extraction.llm_providers import create_provider, DEFAULT_CLAUDE_MODEL
from .enrich.prompts.builder import build_prompt as build_enrich_prompt
from .enrich.providers.base import LLMProvider
from .schemas.base import BaseProvenance, SampleMetadata
from .output.formatters import (
    MetadataFormatter,
    export_to_csv,
    create_provenance_summary,
)
from .output.traceable_format import create_traceable_report
from .output.reports import (
    generate_tabular_report,
    generate_summary_report,
    perform_quality_review,
)
from .fetchers.search import (
    build_query,
    search_sra_runinfo,
    search_riboseq_runs_for_organism,
    DEFAULT_RIBO_TERMS,
    search_runs_with_query,
    fetch_runs_for_bioprojects,
)
from .discovery.candidate_search import discover_candidates, runs_from_candidates
from .fetchers.search import build_terms_for_strategy
import yaml
def serialize_metadata(obj):
    """Custom JSON serializer for Pydantic models and datetime."""
    if isinstance(obj, BaseProvenance):
        return obj.model_dump()
    elif isinstance(obj, datetime):
        return obj.isoformat()
    elif hasattr(obj, "model_dump"):
        return obj.model_dump()
    else:
        return str(obj)
def extract_command(args):
    """Extract metadata for a BioProject."""
    print(f"Extracting metadata for {args.bioproject}")
    print("-" * 60)
    try:
        # Build complete metadata
        metadata = build_project_metadata(args.bioproject)
        # Create traceable report (hierarchical format with all levels)
        output = create_traceable_report(
            study=metadata["study"],
            samples=metadata["samples"],
            runs=metadata["runs"],
            include_raw=True,
            include_statistics=True,
        )
        # Determine output path
        output_path = Path(args.output) if args.output else Path(f"{args.bioproject}_metadata.json")
        # Write JSON output
        with open(output_path, "w") as f:
            json.dump(output, f, indent=2, default=serialize_metadata)
        # Export to CSV if requested
        if getattr(args, 'csv', False):
            csv_path = output_path.with_suffix('.csv')
            export_to_csv(
                metadata["samples"],
                str(csv_path),
                include_provenance=getattr(args, 'csv_provenance', False)
            )
            print(f"  CSV export: {csv_path}")
        # Create provenance summary if requested
        if getattr(args, 'provenance_report', False):
            summary_path = output_path.with_suffix('.md')
            summary = create_provenance_summary(metadata["samples"])
            with open(summary_path, "w") as f:
                f.write(summary)
            print(f"  Provenance report: {summary_path}")
        print("-" * 60)
        print(f"✓ Metadata extracted successfully!")
        print(f"  Study: {metadata['study'].title.value}")
        print(f"  Samples: {len(metadata['samples'])}")
        print(f"  Runs: {len(metadata['runs'])}")
        print(f"  Format: traceable (hierarchical with all provenance levels)")
        print(f"  Output: {output_path}")
        # Show extraction statistics
        stats = output["extraction_statistics"]
        print(f"\n  Field coverage:")
        for field, data in stats["field_coverage"].items():
            if data["count"] > 0:
                print(f"    - {field}: {data['count']}/{len(metadata['samples'])} ({data['percentage']}%)")
        print(f"\n  Confidence distribution:")
        conf = stats["confidence_distribution"]["overall"]
        total = conf["total_fields"]
        if total > 0:
            print(f"    - High (≥0.8): {conf['high (≥0.8)']} ({100*conf['high (≥0.8)']/total:.0f}%)")
            print(f"    - Medium (0.5-0.8): {conf['medium (0.5-0.8)']} ({100*conf['medium (0.5-0.8)']/total:.0f}%)")
            print(f"    - Low (<0.5): {conf['low (<0.5)']} ({100*conf['low (<0.5)']/total:.0f}%)")
        print(f"\n  Ontology mapping: {stats['ontology_mapping']['overall']['mapped']}/{stats['ontology_mapping']['overall']['total']} ({stats['ontology_mapping']['overall']['percentage']}%)")
        # Show navigation guide
        print(f"\n  Access levels:")
        print(f"    - Quick: samples[<id>].quick_view.<field>")
        print(f"    - Standard: samples[<id>].biological_metadata.<field>.provenance")
        print(f"    - Detailed: samples[<id>].biological_metadata.<field>.extraction_details")
    except Exception as e:
        print(f"✗ Error: {e}", file=sys.stderr)
        if args.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)
from .pipeline.orchestrator import run_normalize, run_export
from .pipeline.orchestrator import run_extract as _run_extract
def normalize_command(args):
    from pathlib import Path
    out = run_normalize(Path(args.input), Path(args.out), schema_version=getattr(args, "schema_version", None), strict=args.strict, ontology_cache=Path(args.ontology_cache) if getattr(args, "ontology_cache", None) else None)
    print(f"✓ Normalized: {out}")

def export_command(args):
    """Export TSV/CSV from one or more JSON artifacts."""
    from pathlib import Path
    inputs = [Path(p) for p in args.inputs]
    ok = run_export(inputs, Path(args.out), fmt=getattr(args, 'format', 'tsv'))
    if ok:
        print(f"✓ Export written to {args.out}")
        return
    print("✗ Error: export failed", file=sys.stderr)
    sys.exit(1)


def _write_run_study_csv(rows, out_path):
    from pathlib import Path
    import csv
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    seen = set()
    with open(out_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["Run", "study_accession"])  # required downstream format
        for r in rows:
            run = (r.get("Run") or r.get("run") or "").strip()
            study = (r.get("SRAStudy") or r.get("Study") or r.get("study_accession") or "").strip()
            if not run or not study:
                continue
            if run in seen:
                continue
            seen.add(run)
            w.writerow([run, study])


def discover_command(args):
    """Cross-source discovery → runs CSV; optional extract/table/evidence."""
    organism = args.organism
    print(f"Discovering runs for organism: {organism}")
    # Build search terms bundle
    terms = []
    if getattr(args, 'strategy', None):
        terms.extend(build_terms_for_strategy(args.strategy))
    if getattr(args, 'term', None):
        terms.extend(args.term)
    terms = list(dict.fromkeys([t for t in terms if t]))

    # If candidates file provided, use it; else discover
    prj_list: list[str] = []
    if getattr(args, "candidates", None):
        import csv as _csv
        from pathlib import Path as _P2
        fp = _P2(args.candidates)
        if fp.exists():
            with open(fp, "r") as _f:
                first = _f.readline()
                delim = "	" if ("	" in first) else ","
                _f.seek(0)
                rdr = _csv.reader(_f, delimiter=delim)
                header = next(rdr, [])
                col = 0
                if header:
                    low = [h.strip().lower() for h in header]
                    if "bioproject" in low: col = low.index("bioproject")
                    elif "id" in low: col = low.index("id")
                    else: prj_list.append(header[0].strip())
                for row in rdr:
                    if not row: continue
                    prj = row[col].strip()
                    if prj.startswith("PRJ"): prj_list.append(prj)
        # de-dup
        prj_list = sorted(list(dict.fromkeys(prj_list)))
        # Create placeholder candidates for evidence/extract flow
        class _C:
            def __init__(self, pid: str):
                self.source = "provided"
                self.id = pid
                self.organism = None
                self.title = None
                self.description = None
                self.pubmed = None
                self.score = 0
                self.hits = []
            def to_row(self):
                return {"source":"provided","id":self.id,"organism":"","title":"","description":"","pubmed":"","score":0,"hits":""}
        cands = [_C(p) for p in prj_list]
    else:
        cands = discover_candidates(
        organism,
        geo_ids=None,
        bioproject_retmax=getattr(args, 'bioproject_retmax', 300),
        extra_terms=terms if not getattr(args, 'query', None) else None,
    )
    print(f"Candidates: {len(cands)}")

    # Fetch runs from SRA (and ENA if requested), grouped by PRJ with source tags
    srcs = [s.strip() for s in str(getattr(args, 'sources', 'bioproject,sra')).split(',') if s.strip()]
    include_ena = 'ena' in srcs
    # Seed source evidence with discovery sources (bioproject/geo)
    prj_source_evidence = {}
    for c in cands:
        prj_source_evidence.setdefault(c.id, set()).add(c.source)

    rows_by_prj, prj_sources = runs_from_candidates(
        cands, include_ena=include_ena, terms=terms, boolean_query=getattr(args, 'query', None)
    )

    # Apply intersection threshold
    min_sources = max(1, int(getattr(args, 'min_sources', 1)))
    kept_rows = []
    for prj, rows in rows_by_prj.items():
        sources = set(prj_sources.get(prj, set())) | set(prj_source_evidence.get(prj, set()))
        if len(sources) >= min_sources:
            kept_rows.extend(rows)

    print(f"Total rows (kept): {len(kept_rows)}")
    _write_run_study_csv(kept_rows, args.output)
    print(f"✓ Wrote {args.output}")

    # Evidence sidecar
    if args.evidence:
        from pathlib import Path
        import csv
        ev = str(Path(args.output).with_suffix('.evidence.tsv'))
        with open(ev, 'w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=["source","id","organism","title","description","pubmed","score","hits"], delimiter='\t')
            w.writeheader()
            for c in cands:
                w.writerow(c.to_row())
        print(f"  Evidence: {ev}")

    # Optional extraction of standardized bundles and optional table
    if getattr(args, 'extract_out', None):
        from pathlib import Path
        Path(args.extract_out).mkdir(parents=True, exist_ok=True)
        prjs = sorted({c.id for c in cands if c.id.startswith('PRJ')})
        for i, prj in enumerate(prjs, 1):
            print(f"  Extract [{i}/{len(prjs)}] {prj}")
            out_path = Path(args.extract_out) / f"{prj}_metadata.json"
            try:
                _run_extract(prj, out_path)
            except Exception as ex:
                print(f"    ! Failed {prj}: {ex}")
        if getattr(args, 'table_out', None):
            files = list(Path(args.extract_out).glob('*_metadata.json'))
            if files:
                ok = run_export(files, Path(args.table_out), fmt='tsv' if str(args.table_out).endswith('.tsv') else 'csv')
                if ok:
                    print(f"  Table: {args.table_out}")


def export_runs_command(args):
    """Export Run→study_accession from one or more traceable JSON bundles."""
    from pathlib import Path
    import glob, json
    files = []
    for p in args.inputs:
        files.extend(glob.glob(p))
    seen = set()
    rows = []
    for fp in files:
        try:
            with open(fp) as f:
                data = json.load(f)
        except Exception:
            continue
        # Prefer SRAStudy from runs if present
        runinfo_map = {}
        for r in data.get('runs', []):
            rid = r.get('run_id') or r.get('run') or r.get('Run')
            study = (
                r.get('study_accession') or
                r.get('SRAStudy') or
                r.get('bioproject_id') or
                data.get('study', {}).get('identifiers', {}).get('bioproject_id')
            )
            if rid and study:
                runinfo_map[rid] = study
        for rid, study in runinfo_map.items():
            if rid in seen:
                continue
            seen.add(rid)
            rows.append((rid, study))
    # write
    _write_run_study_csv([{"Run": r, "study_accession": s} for r, s in rows], args.output)
    print(f"✓ Wrote {args.output}")



def search_command(args):
    """Fast: enumerate candidate BioProjects only (no runs/extract)."""
    from Bio import Entrez
    from .fetchers.geo import search_geo_series
    from .fetchers.sra import _convert_gse_to_bioproject
    from .fetchers.search import bioproject_ids_via_esummary
    from .fetchers.entrez_config import configure as _configure_entrez
    _configure_entrez()

    organism = args.organism
    boolean_query = args.query
    srcs = [s.strip() for s in str(getattr(args, 'sources', 'bioproject,geo')).split(',') if s.strip()]
    min_sources = max(1, int(getattr(args, 'min_sources', 1)))
    retmax = int(getattr(args, 'bioproject_retmax', 500))
    geo_retmax = int(getattr(args, 'geo_retmax', 1000))

    # Collect PRJs per source
    hit_sources = {}

    # BioProject esearch+esummary (chunked)
    if 'bioproject' in srcs:
        try:
            prjs = bioproject_ids_via_esummary(organism, boolean_query=boolean_query, retmax=retmax)
        except Exception:
            prjs = []
        for p in prjs:
            if p.startswith('PRJ'):
                hit_sources.setdefault(p, set()).add('bioproject')

    # GEO → PRJ via Entrez summaries (lightweight mapping)
    if 'geo' in srcs:
        try:
            gses = search_geo_series(organism, boolean_query=boolean_query, retmax=geo_retmax)
        except Exception:
            gses = []
        for gse in gses:
            try:
                prj = _convert_gse_to_bioproject(gse)
            except Exception:
                prj = None
            if prj and prj.startswith('PRJ'):
                hit_sources.setdefault(prj, set()).add('geo')

    # Apply min_sources filter and write outputs
    prjs = sorted([p for p, ss in hit_sources.items() if len(ss) >= min_sources])

    # Write TSV (bioproject, sources)
    from pathlib import Path as _P
    import csv
    outp = _P(args.out)
    outp.parent.mkdir(parents=True, exist_ok=True)
    with open(outp, 'w', newline='') as f:
        w = csv.writer(f, delimiter='	')
        w.writerow(['bioproject','sources'])
        for p in prjs:
            w.writerow([p, ','.join(sorted(hit_sources[p]))])
    print(f"✓ Wrote {outp} ({len(prjs)} candidates)")

    # Optional JSON
    if getattr(args, 'json_out', None):
        import json
        jp = _P(args.json_out)
        jp.parent.mkdir(parents=True, exist_ok=True)
        data = [ {'bioproject': p, 'sources': sorted(list(hit_sources[p]))} for p in prjs ]
        jp.write_text(json.dumps(data, indent=2))
        print(f"  JSON: {jp}")


def _slugify(text: str) -> str:
    import re
    s = text.strip().lower()
    s = re.sub(r"[^a-z0-9]+", "_", s)
    return re.sub(r"_+", "_", s).strip("_") or "dataset"


def auto_generic_command(args):
    """End-to-end: discover → extract → table → runs.csv for any omics strategy.

    The user specifies organism, and either a library strategy (mapped to term bundle)
    or explicit terms. Discovery otherwise mirrors the cross-source path.
    """
    from pathlib import Path
    organism = args.organism
    label = args.strategy or "custom"
    out_dir = Path(args.out_dir or f"output/auto/{_slugify(label)}_{_slugify(organism)}")
    extracted_dir = out_dir / "extracted"
    out_dir.mkdir(parents=True, exist_ok=True)
    extracted_dir.mkdir(parents=True, exist_ok=True)

    print(f"[1/4] Discovering candidates for {organism}…")
    # Candidate discovery relies on BioProject/GEO text; strategy terms are used later for SRA filtering if needed
    cands = discover_candidates(organism, geo_ids=(args.geo or []), bioproject_retmax=getattr(args, 'bioproject_retmax', 300))
    print(f"  Candidates: {len(cands)}")
    if getattr(args, 'evidence', False):
        import csv
        ev = out_dir / "candidates.evidence.tsv"
        with open(ev, 'w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=["source","id","organism","title","description","pubmed","score","hits"], delimiter='\t')
            w.writeheader()
            for c in cands:
                d = c.to_row(); w.writerow(d)
        print(f"  Evidence: {ev}")

    print("[2/4] Extracting standardized bundles…")
    json_files = []
    for i, c in enumerate(cands, 1):
        if not c.id.startswith('PRJ'):
            continue
        out_path = extracted_dir / f"{c.id}_metadata.json"
        print(f"  [{i}/{len(cands)}] {c.id}")
        try:
            _run_extract(c.id, out_path)
            json_files.append(out_path)
            if getattr(args, 'validate', False):
                class _V: pass
                v = _V(); v.input = str(out_path)
                validate_command(v)
        except Exception as ex:
            print(f"    ! Failed {c.id}: {ex}")

    if not json_files:
        print("✗ No extracted bundles; aborting")
        return

    print("[3/4] Generating tabular metadata…")
    table_path = out_dir / "metadata_table.tsv"
    ok = run_export([Path(p) for p in json_files], table_path, fmt="tsv")
    if ok:
        print(f"  Table: {table_path}")
    else:
        print("  ! Failed to generate table")

    print("[4/4] Exporting Run→study CSV…")
    runs_csv = out_dir / "runs.csv"
    # reuse export_runs logic
    class _Args: pass
    _a = _Args(); _a.inputs = [str(p) for p in json_files]; _a.output = str(runs_csv)
    export_runs_command(_a)
    print("✓ Done")
def validate_command(args):
    """Validate extracted metadata."""
    print(f"Validating metadata from {args.input}")
    print("-" * 60)
    try:
        with open(args.input, "r") as f:
            data = json.load(f)
        # Simple validation checks
        issues = []
        warnings = []
        # Check study metadata
        study = data.get("study", {})
        if not study.get("title", {}).get("value"):
            issues.append("Study title is missing")
        if study.get("organism", {}).get("confidence", 0) < 0.5:
            warnings.append("Low confidence organism extraction")
        # Check samples
        samples = data.get("samples", {})
        for sample_id, sample in samples.items():
            # Check for tissue (priority field)
            if not sample.get("tissue"):
                warnings.append(f"Sample {sample_id}: No tissue information")
            # Check organism
            if not sample.get("organism", {}).get("value"):
                issues.append(f"Sample {sample_id}: Missing organism")
        # Check runs
        runs = data.get("runs", [])
        if not runs:
            issues.append("No runs found")
        # Report results
        print(f"Validated {len(samples)} samples and {len(runs)} runs")
        if issues:
            print(f"\n❌ Issues found ({len(issues)}):")
            for issue in issues:
                print(f"  - {issue}")
        if warnings:
            print(f"\n⚠️  Warnings ({len(warnings)}):")
            for warning in warnings[:10]:  # Show first 10
                print(f"  - {warning}")
            if len(warnings) > 10:
                print(f"  ... and {len(warnings) - 10} more")
        if not issues and not warnings:
            print("✓ All validation checks passed!")
        sys.exit(1 if issues else 0)
    except Exception as e:
        print(f"✗ Error: {e}", file=sys.stderr)
        if args.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)
def _ensure_base_provenance(field_dict: dict, default_source: str = "structured_field") -> dict:
    """Ensure a dict has all required BaseProvenance fields."""
    if not isinstance(field_dict, dict):
        return field_dict
    # If it has a value but missing source/confidence, add defaults
    if "value" in field_dict:
        field_dict.setdefault("source", default_source)
        field_dict.setdefault("source_id", "unknown")
        field_dict.setdefault("confidence", 0.9)
    return field_dict
def _batch_enrich_with_provider(
    input_files: List[Path],
    output_dir: Path,
    provider,
    resume: bool = False,
    scheme: str = "default",
    shortlist_k: int = 0,
    few_shot: str = "auto",
    context_profile: str = "extended",
) -> Dict:
    """Simplified stub to avoid syntax errors in dev builds.

    Original implementation is temporarily disabled to unblock CLI changes.
    """
    stats = {"total_files": 0, "processed": 0, "failed": 0, "total_samples": 0, "total_enriched": 0}
    return stats
def load_provider_from_config(model_name: str, config_path: Path):
    """Load LLM provider from config file."""
    if not config_path.exists():
        print(f"✗ Error: Config file not found: {config_path}", file=sys.stderr)
        print(f"   Copy config/model_paths.example.yaml to config/model_paths.yaml", file=sys.stderr)
        sys.exit(1)
    with open(config_path) as f:
        config = yaml.safe_load(f)
    model_config = config.get("local_models", {}).get(model_name)
    if not model_config:
        print(f"✗ Error: Model '{model_name}' not found in config", file=sys.stderr)
        print(f"   Available models: {', '.join(config.get('local_models', {}).keys())}", file=sys.stderr)
        sys.exit(1)
    provider_type = model_config.get("provider", "vllm")
    print(f"Loading {model_name} with {provider_type}...")
    if provider_type == "vllm":
        return create_provider(
            "vllm",
            model_path=model_config["path"],
            tensor_parallel_size=model_config.get("tensor_parallel_size", 1),
            gpu_memory_utilization=model_config.get("gpu_memory_utilization", 0.9),
        )
    elif provider_type == "transformers":
        return create_provider(
            "transformers",
            model_path=model_config["path"],
            device=model_config.get("device", "cuda"),
        )
    else:
        print(f"✗ Error: Unknown provider type: {provider_type}", file=sys.stderr)
        sys.exit(1)
def enrich_command(args):
    """Enrich metadata using LLM (Phase 2 - GPU processing)."""
    print(f"Enriching metadata from {args.input}")
    print("-" * 60)
    try:
        # Load metadata
        input_path = Path(args.input)
        if not input_path.exists():
            print(f"✗ Error: File not found: {args.input}", file=sys.stderr)
            sys.exit(1)
        with open(input_path) as f:
            data = json.load(f)
        # Create provider based on mode
        if args.model:
            # Local model mode
            config_path = Path(args.config) if args.config else Path("config/model_paths.yaml")
            provider = load_provider_from_config(args.model, config_path)
        elif getattr(args, "gemini_project", None) or (getattr(args, "gemini", False) and os.environ.get("GOOGLE_CLOUD_PROJECT")):
            # Gemini via Vertex AI
            gemini_project = getattr(args, "gemini_project", None) or os.environ.get("GOOGLE_CLOUD_PROJECT")
            gemini_location = getattr(args, "gemini_location", None) or os.environ.get("GOOGLE_CLOUD_LOCATION", "us-central1")
            gemini_model = getattr(args, "gemini_model", None) or "gemini-2.5-flash"
            provider = create_provider("gemini", project=gemini_project, location=gemini_location, model=gemini_model)
        else:
            # Claude API mode (default)
            api_key = args.api_key or os.environ.get("ANTHROPIC_API_KEY")
            if not api_key:
                print(
                    "✗ Error: No API key provided. Set ANTHROPIC_API_KEY environment variable or use --api-key",
                    file=sys.stderr,
                )
                sys.exit(1)
            provider = create_provider("claude", api_key=api_key, model=args.claude_model or DEFAULT_CLAUDE_MODEL)
        # Extract study context
        study = data.get("study", {})
        # Handle traceable format (nested in core_metadata) or flat/simple format
        study_title = study.get("quick_view", {}).get("title") or \
                      study.get("core_metadata", {}).get("title", {}).get("value") or \
                      study.get("title", {}).get("value", "")
                      
        study_description = study.get("core_metadata", {}).get("description", {}).get("value") or \
                            study.get("description", {}).get("value", "")
                            
        abstract = study.get("publication", {}).get("abstract", {}).get("value") or \
                   study.get("paper_abstract", {}).get("value")
        
        journal = study.get("publication", {}).get("journal", {}).get("value") or \
                  study.get("journal", {}).get("value")
        authors = study.get("publication", {}).get("authors") or \
                  study.get("authors")
        if isinstance(authors, list):
            authors = ", ".join(authors)
        pub_date = study.get("publication", {}).get("publication_date") or \
                   study.get("publication_date")
        # Get samples
        samples = data.get("samples", {})
        total_samples = len(samples)
        if total_samples == 0:
            print("⚠️  No samples found in metadata file")
            sys.exit(0)
        # Determine which samples need enrichment
        samples_to_enrich = []
        for sample_id, sample_data in samples.items():
            # Check if any priority fields are missing
            missing_fields = []
            priority_fields = ["tissue", "cell_type", "cell_line", "strain", "treatment"]
            for field in priority_fields:
                if not sample_data.get(field) and not sample_data.get("biological_metadata", {}).get(field):
                    missing_fields.append(field)
            if missing_fields:
                samples_to_enrich.append((sample_id, sample_data, missing_fields))
        if args.only_if_missing and not samples_to_enrich:
            print("✓ All samples have complete metadata, no enrichment needed")
            sys.exit(0)
        if args.only_if_missing:
            print(f"Found {len(samples_to_enrich)}/{total_samples} samples with missing fields")
        else:
            samples_to_enrich = [(sid, sdata, []) for sid, sdata in samples.items()]
            print(f"Enriching all {total_samples} samples")
        # Track enrichment statistics
        enriched_count = 0
        fields_added = {"tissue": 0, "cell_type": 0, "cell_line": 0, "treatment": 0, "strain": 0}
        # Enrich each sample
        print("\nEnriching samples...")
        for i, (sample_id, sample_data, missing) in enumerate(samples_to_enrich, 1):
            if args.verbose:
                print(f"  [{i}/{len(samples_to_enrich)}] {sample_id} (missing: {', '.join(missing) or 'none'})")
            # Extract from traceable format to flat SampleMetadata
            quick_view = sample_data.get("quick_view", {})
            bio_meta = sample_data.get("biological_metadata", {})
            # Get bioproject_id with multiple fallbacks
            bioproject_id = quick_view.get("bioproject_id") or \
                            sample_data.get("identifiers", {}).get("bioproject_id") or \
                            sample_data.get("bioproject_id") or \
                            study.get("identifiers", {}).get("bioproject_id") or \
                            study.get("quick_view", {}).get("bioproject_id") or \
                            study.get("bioproject_id") or \
                            study.get("bioproject", {}).get("value")
            flat_data = {
                "sample_id": quick_view.get("sample_id") or sample_data.get("sample_id"),
                "bioproject_id": bioproject_id,
                "organism": _ensure_base_provenance(bio_meta.get("organism", {})),
                "raw_characteristics": sample_data.get("raw_biosample_attributes") or sample_data.get("raw_characteristics"),
                "technical_context": sample_data.get("technical_context"),
            }
            # Extract all BaseProvenance fields (keep as dicts, ensure required fields)
            provenance_fields = [
                "tissue", "cell_type", "cell_line", "strain", "treatment", "disease",
                "developmental_stage", "age", "sex", "genotype", "condition",
                "timepoint", "replicate", "batch", "stress", "temperature",
                "growth_condition", "sample_title", "sample_description"
            ]
            for field in provenance_fields:
                if field in bio_meta:
                    flat_data[field] = _ensure_base_provenance(bio_meta[field])
                elif field in sample_data and isinstance(sample_data.get(field), dict):
                    flat_data[field] = _ensure_base_provenance(sample_data[field])
            sample = SampleMetadata(**flat_data)
            # Get sample-specific context
            sample_title = sample_data.get("sample_title", {}).get("value") or \
                           sample_data.get("descriptions", {}).get("title", {}).get("value")
            
            sample_description = sample_data.get("sample_description", {}).get("value") or \
                                 sample_data.get("descriptions", {}).get("description", {}).get("value")
            # Enrich with LLM
            try:
                # Build context for parity with gold-prepare
                from .context.context_builder import compose_llm_context, ContextOpts
                ctx_opts = ContextOpts(
                    include_geo_overall=getattr(args, 'context_include_geo_overall', True),
                    include_geo_protocols=getattr(args, 'context_include_geo_protocols', True),
                    include_geo_source=getattr(args, 'context_include_geo_source', True),
                    merge_channels=getattr(args, 'context_merge_channels', True),
                    include_bioproject_type=getattr(args, 'context_include_bioproject_type', True),
                    include_bioproject_scope=getattr(args, 'context_include_bioproject_scope', True),
                    include_pubmed_mesh=getattr(args, 'context_include_pubmed_mesh', True),
                    include_biosample_package=getattr(args, 'context_include_biosample_package', False),
                    include_biosample_description=getattr(args, 'context_include_biosample_description', True),
                    truncate_chars=max(0, int(getattr(args, 'context_truncate', 0) or 0)),
                )
                extra_ctx = compose_llm_context(study, sample, profile=getattr(args, 'context_profile', 'extended'), opts=ctx_opts)
                # Build enriched prompt with ontology shortlists and few-shot (Option A wiring)
                try:
                    desired_prompt = build_enrich_prompt(
                        study_title=study_title,
                        study_description=study_description,
                        sample=sample,
                        abstract=abstract,
                        journal=journal,
                        authors=authors,
                        publication_date=pub_date,
                        scheme=args.scheme,
                        shortlist_k=getattr(args, 'ontology_shortlist', 0),
                        few_shot=getattr(args, 'few_shot', 'auto'),
                        extra_context=extra_ctx,
                    )
                except Exception:
                    # Fallback to basic prompt via enrich_sample_metadata path
                    desired_prompt = None

                class _FixedPromptProvider(LLMProvider):
                    def __init__(self, base, fixed_prompt: Optional[str]):
                        self.base = base
                        self.fixed_prompt = fixed_prompt
                    def extract(self, prompt: str, temperature: float = 0.0, max_tokens: int = 1000):
                        return self.base.extract(self.fixed_prompt or prompt, temperature=temperature, max_tokens=max_tokens)
                    def get_model_name(self) -> str:
                        return self.base.get_model_name()
                    def supports_batching(self) -> bool:
                        return self.base.supports_batching()

                wrapped_provider = _FixedPromptProvider(provider, desired_prompt)

                enriched_sample = enrich_sample_metadata(
                    sample=sample,
                    study_title=study_title,
                    study_description=study_description,
                    sample_title=sample_title,
                    sample_description=sample_description,
                    abstract=abstract,
                    journal=journal,
                    authors=authors,
                    publication_date=pub_date,
                    provider=wrapped_provider,
                    source_id=f"llm_enrichment_{sample_id}",
                    scheme=args.scheme,
                    extra_sections=extra_ctx,
                )
                # Update traceable format with enriched values
                all_fields = [
                    "tissue", "cell_type", "cell_line", "strain", "genotype", 
                    "sex", "age", "developmental_stage", "condition", "treatment", 
                    "timepoint", "replicate", "batch", "disease", "stress", 
                    "temperature", "growth_condition", "library_strategy"
                ]
                for field in all_fields:
                    new_field = getattr(enriched_sample, field, None)
                    
                    # Detect if we are in traceable format or flat format
                    if is_traceable:
                        old_val = bio_meta.get(field, {}).get("value")
                        # new_field is a BaseProvenance object, extract .value
                        # Use getattr safely in case of unexpected types
                        new_val = getattr(new_field, "value", None)
                        if new_field and new_val and not old_val:
                            if field not in bio_meta:
                                bio_meta[field] = {}
                            bio_meta[field]["value"] = new_val
                            bio_meta[field]["source"] = getattr(new_field, "source", "llm_enrichment")
                            bio_meta[field]["confidence"] = getattr(new_field, "confidence", 0.5)
                            quick_view[field] = new_val
                            fields_added.setdefault(field, 0)
                            fields_added[field] += 1
                    else:
                        # Flat format (direct sample_data)
                        old_val = sample_data.get(field, {}).get("value") if isinstance(sample_data.get(field), dict) else sample_data.get(field)
                        new_val = getattr(new_field, "value", None)
                        if new_field and new_val and not old_val:
                            sample_data[field] = {
                                "value": new_val,
                                "source": getattr(new_field, "source", "llm_enrichment"),
                                "source_id": getattr(new_field, "source_id", "llm"),
                                "confidence": getattr(new_field, "confidence", 0.5)
                            }
                            fields_added.setdefault(field, 0)
                            fields_added[field] += 1
                enriched_count += 1
            except Exception as e:
                if args.verbose:
                    print(f"    ⚠️  Failed to enrich {sample_id}: {e}")
                continue
        # Save enriched metadata
        output_path = Path(args.output) if args.output else input_path.with_stem(f"{input_path.stem}_enriched")
        data["samples"] = samples
        data["enrichment_metadata"] = {
            "timestamp": datetime.now().isoformat(),
            "samples_enriched": enriched_count,
            "fields_added": fields_added,
        }
        with open(output_path, "w") as f:
            json.dump(data, f, indent=2, default=serialize_metadata)
        # Report results
        print(f"\n✓ Enriched {enriched_count}/{len(samples_to_enrich)} samples")
        print(f"  Fields added:")
        for field, count in fields_added.items():
            if count > 0:
                print(f"    - {field}: {count}")
        print(f"\nSaved to: {output_path}")
    except Exception as e:
        print(f"✗ Error: {e}", file=sys.stderr)
        if args.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)
def batch_enrich_command(args):
    """Batch enrich multiple metadata files (optimized for GPU nodes)."""
    print(f"Batch enriching {len(args.inputs)} files")
    print("-" * 60)
    try:
        # Create provider based on mode
        if args.model:
            # Local model mode
            config_path = Path(args.config) if args.config else Path("config/model_paths.yaml")
            provider = load_provider_from_config(args.model, config_path)
        elif getattr(args, "gemini_project", None) or (getattr(args, "gemini", False) and os.environ.get("GOOGLE_CLOUD_PROJECT")):
            # Gemini via Vertex AI
            gemini_project = getattr(args, "gemini_project", None) or os.environ.get("GOOGLE_CLOUD_PROJECT")
            gemini_location = getattr(args, "gemini_location", None) or os.environ.get("GOOGLE_CLOUD_LOCATION", "us-central1")
            gemini_model = getattr(args, "gemini_model", None) or "gemini-2.5-flash"
            provider = create_provider("gemini", project=gemini_project, location=gemini_location, model=gemini_model)
        else:
            # Claude API mode (default)
            api_key = args.api_key or os.environ.get("ANTHROPIC_API_KEY")
            if not api_key:
                print(
                    "✗ Error: No API key provided. Set ANTHROPIC_API_KEY environment variable or use --api-key",
                    file=sys.stderr,
                )
                sys.exit(1)
            provider = create_provider("claude", api_key=api_key, model=args.claude_model or DEFAULT_CLAUDE_MODEL)
        # Parse input files
        input_files = []
        for pattern in args.inputs:
            path = Path(pattern)
            if path.is_file():
                input_files.append(path)
            elif "*" in pattern:
                # Glob pattern
                input_files.extend(Path().glob(pattern))
            else:
                print(f"⚠️  Skipping invalid input: {pattern}")
        if not input_files:
            print("✗ Error: No valid input files found", file=sys.stderr)
            sys.exit(1)
        # Create output directory
        output_dir = Path(args.output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        print(f"Processing {len(input_files)} files...")
        print(f"Output directory: {output_dir}")
        if args.model or getattr(args, "gemini_project", None) or (getattr(args, "gemini", False) and os.environ.get("GOOGLE_CLOUD_PROJECT")):
            # Local model or Gemini mode - process sequentially with loaded provider
            print(f"Using {provider.get_model_name()} (sequential processing)")
            stats = _batch_enrich_with_provider(
                input_files=input_files,
                output_dir=output_dir,
                provider=provider,
                resume=args.resume,
                scheme=getattr(args, 'scheme', 'default'),
                shortlist_k=getattr(args, 'ontology_shortlist', 0),
                few_shot=getattr(args, 'few_shot', 'auto'),
                context_profile=getattr(args, 'context_profile', 'extended'),
            )
        else:
            # Claude API mode - use concurrent workers
            print(f"Workers: {args.workers}, Rate limit: {args.rate_limit} req/min")
            stats = enrich_batch_from_files(
                input_files=input_files,
                output_dir=output_dir,
                api_key=api_key,
                max_workers=args.workers,
                rate_limit_rpm=args.rate_limit,
                resume=args.resume,
            )
        # Report results
        print("\n" + "=" * 60)
        print("Batch enrichment complete!")
        print(f"  Files processed: {stats['processed']}/{stats['total_files']}")
        print(f"  Files failed: {stats['failed']}")
        print(f"  Total samples: {stats['total_samples']}")
        print(f"  Samples enriched: {stats['total_enriched']}")
        print(f"\nOutput files saved to: {output_dir}")
    except Exception as e:
        print(f"✗ Error: {e}", file=sys.stderr)
        if args.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)
def report_command(args):
    """Generate reports from enriched metadata files."""
    # Find all input files (handle glob manually if shell didn't expand)
    import glob
    all_inputs = []
    for pattern in args.inputs:
        files = glob.glob(str(pattern))
        if not files:
            # If no glob match, check if it's a direct path
            p = Path(pattern)
            if p.exists():
                all_inputs.append(p)
        else:
            all_inputs.extend([Path(f) for f in files])
    if not all_inputs:
        print(f"✗ Error: No input files found matching patterns: {args.inputs}", file=sys.stderr)
        sys.exit(1)
    print(f"Generating {args.format} report from {len(all_inputs)} files...")
    if args.format == "table":
        output_path = Path(args.output) if args.output else Path("metadata_summary.tsv")
        success = generate_tabular_report(all_inputs, output_path, format_type=args.table_format)
        if success:
            print(f"✓ Tabular report saved to {output_path}")
        else:
            print("✗ Error: Failed to generate tabular report", file=sys.stderr)
            sys.exit(1)
    elif args.format == "summary":
        generate_summary_report(all_inputs)
    elif args.format == "review":
        perform_quality_review(all_inputs, num_to_review=args.num_samples)
def main():
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(
        description="Extract metadata from omics sequencing projects",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Phase 1: Extract metadata (CPU nodes)
  omics-extract extract PRJNA1170270
  # Phase 2: Enrich with local model (GPU nodes)
  omics-extract enrich PRJNA1170270_metadata.json --model qwen-2.5-7b
  # Enrich with Claude API
  omics-extract enrich PRJNA1170270_metadata.json
  # Enrich with Gemini via Vertex AI
  omics-extract enrich PRJNA1170270_metadata.json --gemini-project my-gcp-project --gemini-location europe-west2
  # Enrich with Gemini (project from env var)
  omics-extract enrich PRJNA1170270_metadata.json --gemini
  # Batch enrich multiple projects with local model
  omics-extract batch-enrich *_metadata.json --output-dir enriched/ --model mistral-7b
  # Batch enrich with Gemini
  omics-extract batch-enrich *_metadata.json --output-dir enriched/ --gemini-project my-gcp-project
  # Batch enrich with Claude API (GPU cluster)
  omics-extract batch-enrich *_metadata.json --output-dir enriched/ --workers 8
  # Resume interrupted batch job
  omics-extract batch-enrich *_metadata.json --output-dir enriched/ --model qwen-2.5-7b --resume
  # Full workflow with local model
  omics-extract extract PRJNA1170270 --output data.json
  omics-extract enrich data.json --model qwen-2.5-7b --only-if-missing
        """,
    )
    parser.add_argument("--verbose", "-v", action="store_true", help="Verbose output")
    parser.add_argument("--log-json", action="store_true", help="Log in JSON format for ingestion")
    parser.add_argument("--metrics", help="Write metrics JSON to this path")
    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    # Discover command (Run→study CSV)
    discover_parser = subparsers.add_parser(
        "discover",
        help="Cross-source discovery → runs CSV; optional extract/table/evidence",
    )
    discover_parser.add_argument("organism", help="Organism scientific name (e.g., 'Danio rerio')")
    # Query controls
    discover_parser.add_argument("--query", help="Boolean query applied across sources (e.g., '(Ribo-Seq OR \"ribosome profiling\")')")
    discover_parser.add_argument("--strategy", help="Library strategy hint (maps to curated synonyms)")
    discover_parser.add_argument("--term", action="append", help="Additional search term (repeatable)")
    # Sources and thresholds
    discover_parser.add_argument("--sources", default="bioproject,sra", help="Comma-separated: bioproject,sra,geo,ena (default: bioproject,sra)")
    discover_parser.add_argument("--min-sources", type=int, default=1, help="Keep studies seen in ≥ N sources (default: 1)")
    discover_parser.add_argument("--bioproject-retmax", type=int, default=300, help="Cap BioProject results (default: 300)")
    # Outputs
    discover_parser.add_argument("--runs-out", dest="output", default="runs_discover.csv", help="Output runs CSV (Run,study_accession)")
    discover_parser.add_argument("--extract-out", help="Write standardized JSON bundles to this directory")
    discover_parser.add_argument("--table-out", help="Write normalized TSV/CSV table to this path")
    discover_parser.add_argument("--evidence", action="store_true", help="Write evidence TSV next to outputs")
    discover_parser.add_argument("--candidates", help="Input TSV/CSV with a list of BioProject IDs to expand to runs")
    discover_parser.set_defaults(func=discover_command)

    # Export runs command: read JSON bundle(s), emit Run→study_accession with robust mapping
    export_runs = subparsers.add_parser("export-runs", help="Export Run→study_accession CSV from extracted JSON bundles")
    export_runs.add_argument("inputs", nargs='+', help="Input JSON files or globs")
    export_runs.add_argument("--out", "--output", dest="output", required=True, help="Output CSV path")
    export_runs.set_defaults(func=export_runs_command)

    # Auto pipeline (generic end-to-end)
    auto = subparsers.add_parser("auto", help="End-to-end: discover→extract→table→runs; optional validate")
    auto.add_argument("organism", help="Organism scientific name (e.g., 'Danio rerio')")
    auto.add_argument("--strategy", help="Library strategy hint (e.g., 'Ribo-Seq', 'RNA-Seq'); expands to sensible term bundle")
    auto.add_argument("--term", action="append", help="Additional search term (repeatable)")
    auto.add_argument("--out-dir", help="Output directory (default: output/auto/<strategy|custom>_<organism>)")
    auto.add_argument("--geo", action="append", help="Seed GEO Series ID (GSE…) (repeatable)")
    auto.add_argument("--bioproject-retmax", type=int, default=300, help="Max BioProject records to scan (default: 300)")
    auto.add_argument("--evidence", action="store_true", help="Write candidates.evidence.tsv")
    auto.add_argument("--validate", action="store_true", help="Validate each extracted bundle (off by default)")
    auto.set_defaults(func=auto_generic_command)
    # Extract command
    extract_parser = subparsers.add_parser("extract", help="Extract metadata for a BioProject (Phase 1 - CPU)")
    extract_parser.add_argument("bioproject", help="BioProject accession (PRJNA...)")
    extract_parser.add_argument(
        "--output", "-o", help="Output JSON file (default: <bioproject>_metadata.json)"
    )
    extract_parser.add_argument(
        "--detail",
        choices=["simple", "standard", "detailed"],
        default="standard",
        help="Detail level for output (default: standard). "
             "simple=values only, standard=values+provenance, detailed=complete extraction chain"
    )
    extract_parser.add_argument(
        "--csv",
        action="store_true",
        help="Also export to CSV format"
    )
    extract_parser.add_argument(
        "--csv-provenance",
        action="store_true",
        help="Include source/confidence columns in CSV (requires --csv)"
    )
    extract_parser.add_argument(
        "--provenance-report",
        action="store_true",
        help="Generate markdown provenance summary report"
    )
    extract_parser.set_defaults(func=extract_command)
    # Normalize command (Phase 1.5)
    normalize_parser = subparsers.add_parser("normalize", help="Re-run ontology normalization on existing JSON")
    normalize_parser.add_argument("input", help="Input JSON file from extract")
    normalize_parser.add_argument("--out", "-o", required=True, help="Output JSON file")
    normalize_parser.add_argument("--schema-version", help="Target schema version tag")
    normalize_parser.add_argument("--ontology-cache", help="Ontology cache directory")
    normalize_parser.add_argument("--strict", action="store_true", help="Drop unmapped fields instead of keeping raw")
    normalize_parser.set_defaults(func=normalize_command)
    
    # Enrich command (Phase 2)
    enrich_parser = subparsers.add_parser("enrich", help="Enrich metadata with LLM (Phase 2 - GPU)")
    enrich_parser.add_argument("input", help="Input metadata JSON file from extract command")
    enrich_parser.add_argument(
        "--output", "-o", help="Output JSON file (default: <input>_enriched.json)"
    )
    enrich_parser.add_argument("--provider", default="vllm", help="LLM provider: vllm, server-vllm, claude, gemini, transformers")
    enrich_parser.add_argument("--model", help="Model ID or path (per provider)")
    enrich_parser.add_argument("--server-url", help="Server vLLM base URL (for provider=server-vllm)")
    enrich_parser.add_argument(
        "--config", help="Path to model config YAML (default: config/model_paths.yaml)"
    )
    enrich_parser.add_argument(
        "--api-key", help="Anthropic API key (or set ANTHROPIC_API_KEY env var) - for Claude mode"
    )
    enrich_parser.add_argument(
        "--claude-model", help=f"Claude model name (default: {DEFAULT_CLAUDE_MODEL})"
    )
    enrich_parser.add_argument(
        "--gemini", action="store_true", help="Use Gemini via Vertex AI (reads project from GOOGLE_CLOUD_PROJECT env)"
    )
    enrich_parser.add_argument(
        "--gemini-project", help="Google Cloud project ID for Vertex AI (or set GOOGLE_CLOUD_PROJECT)"
    )
    enrich_parser.add_argument(
        "--gemini-location", default="us-central1", help="GCP region for Vertex AI (default: us-central1)"
    )
    enrich_parser.add_argument(
        "--gemini-model", default="gemini-2.5-flash", help="Gemini model name (default: gemini-2.5-flash)"
    )
    enrich_parser.add_argument(
        "--only-if-missing",
        action="store_true",
        help="Only enrich samples with missing critical fields",
    )
    enrich_parser.add_argument("--ontology-shortlist", type=int, default=0, help="Top-K ontology candidates to include in prompt (0=off)")
    enrich_parser.add_argument("--few-shot", choices=["auto","on","off"], default="auto", help="Include few-shot exemplars")
    enrich_parser.add_argument(
        "--scheme",
        default="default",
        help="Extraction scheme to use (default, nanopore_rna, ribo_seq, or custom YAML name)"
    )
    # Context parity flags (match gold-prepare)
    enrich_parser.add_argument("--context-profile", choices=["minimal","standard","extended","full"], default="extended", help="LLM context profile")
    enrich_parser.add_argument("--context-truncate", type=int, default=0, help="Truncate each context block (0=no limit)")
    enrich_parser.add_argument("--context-merge-channels", dest="context_merge_channels", action="store_true", default=True)
    enrich_parser.add_argument("--context-separate-channels", dest="context_merge_channels", action="store_false")
    enrich_parser.add_argument("--context-include-geo-overall", action="store_true", default=True)
    enrich_parser.add_argument("--context-include-geo-protocols", action="store_true", default=True)
    enrich_parser.add_argument("--context-include-geo-source", action="store_true", default=True)
    enrich_parser.add_argument("--context-include-bioproject-type", action="store_true", default=True)
    enrich_parser.add_argument("--context-include-bioproject-scope", action="store_true", default=True)
    enrich_parser.add_argument("--context-include-pubmed-mesh", action="store_true", default=True)
    enrich_parser.add_argument("--context-include-biosample-package", action="store_true", default=False)
    enrich_parser.add_argument("--context-include-biosample-description", action="store_true", default=True)
    enrich_parser.set_defaults(func=enrich_command)
    # Batch enrich command (Phase 2 - optimized for GPU clusters)
    batch_parser = subparsers.add_parser("batch-enrich", help="Batch enrich multiple files (GPU-optimized)")
    batch_parser.add_argument("inputs", nargs="+", help="Input metadata files or glob patterns")
    batch_parser.add_argument(
        "--output-dir", "-o", required=True, help="Output directory for enriched files"
    )
    batch_parser.add_argument("--provider", default="vllm", help="LLM provider: vllm, server-vllm, claude, gemini, transformers")
    batch_parser.add_argument("--model", help="Model ID or path (per provider)")
    batch_parser.add_argument("--server-url", help="Server vLLM base URL (for provider=server-vllm)")
    batch_parser.add_argument(
        "--config", help="Path to model config YAML (default: config/model_paths.yaml)"
    )
    batch_parser.add_argument(
        "--api-key", help="Anthropic API key (or set ANTHROPIC_API_KEY env var) - for Claude mode"
    )
    batch_parser.add_argument(
        "--claude-model", help=f"Claude model name (default: {DEFAULT_CLAUDE_MODEL})"
    )
    batch_parser.add_argument(
        "--gemini", action="store_true", help="Use Gemini via Vertex AI (reads project from GOOGLE_CLOUD_PROJECT env)"
    )
    batch_parser.add_argument(
        "--gemini-project", help="Google Cloud project ID for Vertex AI (or set GOOGLE_CLOUD_PROJECT)"
    )
    batch_parser.add_argument(
        "--gemini-location", default="us-central1", help="GCP region for Vertex AI (default: us-central1)"
    )
    batch_parser.add_argument(
        "--gemini-model", default="gemini-2.5-flash", help="Gemini model name (default: gemini-2.5-flash)"
    )
    batch_parser.add_argument(
        "--workers", type=int, default=4, help="Concurrent API requests (default: 4, ignored for local models)"
    )
    batch_parser.add_argument(
        "--rate-limit", type=int, default=50, help="Max requests per minute (default: 50, ignored for local models)"
    )
    batch_parser.add_argument(
        "--resume", action="store_true", help="Resume from checkpoints"
    )
    batch_parser.add_argument(
        "--scheme",
        default="default",
        help="Extraction scheme to use (default, nanopore_rna, ribo_seq, or custom YAML name)"
    )
    batch_parser.set_defaults(func=batch_enrich_command)
    # Validate command
    validate_parser = subparsers.add_parser("validate", help="Validate extracted metadata")
    validate_parser.add_argument("input", help="Input JSON file to validate")
    validate_parser.set_defaults(func=validate_command)
    # Report command (New in Phase 4)
    report_parser = subparsers.add_parser("report", help="Generate reports from enriched metadata")
    report_parser.add_argument("inputs", nargs="+", help="Input enriched JSON files or glob patterns")
    report_parser.add_argument(
        "--format",
        choices=["table", "summary", "review"],
        default="summary",
        help="Report type: table (TSV/CSV), summary (stats), or review (quality check)"
    )
    report_parser.add_argument(
        "--output", "-o", help="Output file path (for 'table' format)"
    )
    report_parser.add_argument(
        "--table-format",
        choices=["tsv", "csv"],
        default="tsv",
        help="Format for tabular report (default: tsv)"
    )
    report_parser.add_argument(
        "--num-samples",
        "-n",
        type=int,
        default=10,
        help="Number of samples to show in 'review' mode (default: 10)"
    )
    report_parser.set_defaults(func=report_command)
    
    # Export command (tabular)
    export_parser = subparsers.add_parser("export", help="Export TSV/CSV from JSON artifacts")
    export_parser.add_argument("inputs", nargs="+", help="Input JSON files (enriched or not)")
    export_parser.add_argument("--format", choices=["tsv","csv"], default="tsv", help="Output table format")
    export_parser.add_argument("--out", "-o", required=True, help="Output file path")
    export_parser.set_defaults(func=export_command)

    # Search command (one-query-per-resource, integrates to BioProjects, fetches all runs)
        # Search (fast): list candidate BioProjects only
    search_p = subparsers.add_parser("search", help="Fast candidate BioProject search (no runs)")
    search_p.add_argument("organism", help="Organism scientific name")
    search_p.add_argument("--query", required=False, help="Boolean query applied across sources")
    search_p.add_argument("--sources", default="bioproject,geo", help="Comma-separated: bioproject,geo")
    search_p.add_argument("--out", required=True, help="Output TSV (bioproject\tsources)")
    search_p.add_argument("--json-out", help="Optional JSON output of candidates")
    search_p.add_argument("--min-sources", type=int, default=1, help="Keep PRJs seen in ≥ N sources (default: 1)")
    search_p.add_argument("--bioproject-retmax", type=int, default=500, help="Cap BioProject esearch (default 500)")
    search_p.add_argument("--geo-retmax", type=int, default=1000, help="Cap GEO series search (default 1000)")
    search_p.set_defaults(func=search_command)


    # Gold set preparation
    gold_parser = subparsers.add_parser("gold-prepare", help="Prepare gold-set curation packets from JSON artifacts")
    gold_parser.add_argument("inputs", nargs="+", help="Input JSON files (traceable or enriched)")
    gold_parser.add_argument("--out", "-o", required=True, help="Output file path (.jsonl/.tsv/.csv)")
    gold_parser.add_argument("--format", choices=["jsonl","tsv","csv"], default="jsonl", help="Output format (default: jsonl)")
    gold_parser.add_argument("--fields", help="Comma-separated target fields (default: organism,tissue,cell_type,cell_line,treatment,disease,genotype,strain,age,sex,developmental_stage)")
    gold_parser.add_argument("--ontology-shortlist", type=int, default=0, help="Include top-K ontology candidates per field (0=off)")
    gold_parser.add_argument("--context-profile", choices=["minimal","standard","extended","full"], default="extended", help="Context profile for gold packets")
    gold_parser.add_argument("--context-truncate", type=int, default=0, help="Truncate each context block (0=no limit)")
    gold_parser.add_argument("--context-merge-channels", dest="context_merge_channels", action="store_true", default=True)
    gold_parser.add_argument("--context-separate-channels", dest="context_merge_channels", action="store_false")
    gold_parser.add_argument("--context-include-geo-overall", action="store_true", default=True)
    gold_parser.add_argument("--context-include-geo-protocols", action="store_true", default=True)
    gold_parser.add_argument("--context-include-geo-source", action="store_true", default=True)
    gold_parser.add_argument("--context-include-bioproject-type", action="store_true", default=True)
    gold_parser.add_argument("--context-include-bioproject-scope", action="store_true", default=True)
    gold_parser.add_argument("--context-include-pubmed-mesh", action="store_true", default=True)
    gold_parser.add_argument("--context-include-biosample-package", action="store_true", default=False)
    gold_parser.add_argument("--context-include-biosample-description", action="store_true", default=True)
    gold_parser.add_argument("--context-include-sra-runs", action="store_true", default=False, help="Include per-sample SRA runs in context (assistive only; not shown to LLM)")
    gold_parser.set_defaults(func=gold_prepare_command)

    # Gold serve (HTML UI)
    gold_serve = subparsers.add_parser("gold-serve", help="Serve a lightweight HTML UI for gold-set curation")
    gold_serve.add_argument("input", help="Input gold packet (jsonl/tsv/csv)")
    gold_serve.add_argument("--out", "-o", required=True, help="Output curated JSONL")
    gold_serve.add_argument("--host", default="127.0.0.1", help="Host to bind (default: 127.0.0.1)")
    gold_serve.add_argument("--port", type=int, default=8765, help="Port (default: 8765)")
    gold_serve.add_argument("--read-only", action="store_true", help="Do not allow edits; browse only")
    gold_serve.add_argument("--open", action="store_true", help="Open browser on start")
    gold_serve.set_defaults(func=gold_serve_command)
    
    
    # Parse and execute
    args = parser.parse_args()
    # Configure logging if verbose
    if args.verbose:
        import logging
        logging.basicConfig(level=logging.DEBUG, format="%(levelname)s:%(name)s:%(message)s")
    if not args.command:
        parser.print_help()
        sys.exit(1)
    exit_code = 0
    try:
        args.func(args)
    except SystemExit as e:
        exit_code = int(e.code) if hasattr(e, "code") else 1
        raise
    finally:
        if getattr(args, "metrics", None):
            try:
                from .core.metrics import CmdMetrics
                import json as _json
                from datetime import datetime as _dt
                m = CmdMetrics(cmd=str(args.command), input=vars(args), exit_code=exit_code).finish(exit_code=exit_code)
                with open(args.metrics, "w") as _f:
                    _f.write(_json.dumps(m.model_dump(), indent=2))
                print(f"Metrics written to {args.metrics}")
            except Exception as _ex:
                print(f"Warning: failed to write metrics: {_ex}")
if __name__ == "__main__":
    main()


from .evaluation.goldset import prepare_goldset as _prepare_gold
from .context.context_builder import ContextOpts
from .ui.gold_server import build_app

def gold_prepare_command(args):
    from pathlib import Path
    fields = None
    if getattr(args, 'fields', None):
        fields = [x.strip() for x in args.fields.split(',') if x.strip()]
    opts = ContextOpts(
        include_geo_overall=getattr(args, 'context_include_geo_overall', True),
        include_geo_protocols=getattr(args, 'context_include_geo_protocols', True),
        include_geo_source=getattr(args, 'context_include_geo_source', True),
        merge_channels=getattr(args, 'context_merge_channels', True),
        include_bioproject_type=getattr(args, 'context_include_bioproject_type', True),
        include_bioproject_scope=getattr(args, 'context_include_bioproject_scope', True),
        include_pubmed_mesh=getattr(args, 'context_include_pubmed_mesh', True),
        include_biosample_package=getattr(args, 'context_include_biosample_package', False),
        include_biosample_description=getattr(args, 'context_include_biosample_description', True),
        truncate_chars=max(0, int(getattr(args, 'context_truncate', 0) or 0)),
    )
    n, guide = _prepare_gold(
        [Path(p) for p in args.inputs],
        Path(args.out),
        fmt=args.format,
        fields=fields,
        shortlist_k=getattr(args, 'ontology_shortlist', 0),
        context_profile=getattr(args, 'context_profile', 'extended'),
        context_opts=opts,
        include_sra_runs=getattr(args, 'context_include_sra_runs', False),
    )
    print(f"✓ Prepared {n} curation packets -> {args.out}")


def gold_serve_command(args):
    from pathlib import Path
    import uvicorn, webbrowser
    app = build_app(Path(args.input), Path(args.out), read_only=getattr(args, 'read_only', False))
    if getattr(args, 'open', False):
        webbrowser.open(f"http://{args.host}:{args.port}")
    uvicorn.run(app, host=args.host, port=args.port, log_level="info")
