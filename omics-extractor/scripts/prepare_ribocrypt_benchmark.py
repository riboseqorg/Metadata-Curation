import csv
import json
from pathlib import Path
from typing import List, Dict, Any
import sys
import os

# Add src to path to import omics_extractor
sys.path.append(str(Path(__file__).parent.parent / "src"))

from omics_extractor.extraction.builder import build_study_metadata, build_sample_metadata, build_run_metadata
from omics_extractor.evaluation.benchmark import GoldStandardSample

CSV_PATH = "/Users/jackt/projects/Metadata-Curation/resources/RiboCrypt_Metadata_13_09_24.csv"
OUTPUT_PATH = Path(__file__).parent.parent / "tests/data/ribocrypt_gold.json"

def prepare_benchmark(limit_projects=10):
    print(f"Reading RiboCrypt metadata from {CSV_PATH}...")
    
    projects = {} # bioproject -> [samples]
    
    with open(CSV_PATH, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            bp = row.get("BioProject") or row.get("study_accession")
            if not bp:
                continue
            
            if bp not in projects:
                projects[bp] = []
            
            # Simplified ground truth mapping
            # RiboCrypt often uses NA or empty
            def clean(val):
                if not val or val.upper() == "NA":
                    return None
                return val.lower()

            # RiboCrypt specific mapping
            # INHIBITOR usually contains ribo-seq specific like CHX
            # CONDITION often contains the more interesting experimental variable
            hib = clean(row.get("INHIBITOR"))
            cond = clean(row.get("CONDITION"))
            
            # RiboCrypt specific mapping
            hib = clean(row.get("INHIBITOR"))
            cond = clean(row.get("CONDITION"))
            
            # Combine condition and inhibitor for a broader treatment field
            treatment_list = []
            if cond and cond.lower() not in ["none", "wt", "na", "-"]:
                treatment_list.append(cond)
            if hib and hib.lower() not in ["none", "wt", "na", "-", "none (chx in lysis and gradient)"]:
                treatment_list.append(hib)
            
            treatment = " ".join(treatment_list) if treatment_list else None
            if not treatment and cond and cond.lower() != "wt":
                treatment = cond

            strain = clean(row.get("GENE"))
            if strain and strain.lower() in ["none", "wt", "na", "-", "wt/nt", "normal"]:
                strain = None

            tissue = clean(row.get("TISSUE"))
            if tissue and tissue.lower() in ["none", "na", "-"]:
                tissue = None
            
            cell_line = clean(row.get("CELL_LINE"))
            if cell_line and cell_line.lower() in ["none", "na", "-"]:
                cell_line = None

            gt = {
                "tissue": tissue,
                "cell_type": None,  # RiboCrypt CSV doesn't have a distinct cell_type column
                "cell_line": cell_line,
                "condition": cond,
                "treatment": treatment,
                "timepoint": clean(row.get("TIMEPOINT")),
                "replicate": clean(row.get("REPLICATE")),
                "organism": clean(row.get("ScientificName")).lower() if row.get("ScientificName") else None,
                "strain": strain,
            }
            
            # Only keep if has some ground truth
            if any(gt.values()):
                projects[bp].append({
                    "sample_id": row.get("Run"),
                    "biosample_id": row.get("BioSample"),
                    "ground_truth": gt,
                    "row": dict(row) # Capture the current row state
                })

    print(f"Found {len(projects)} projects with some metadata.")
    
    # Select a diverse set of samples
    selected_samples: List[GoldStandardSample] = []
    
    # Randomly select projects
    import random
    bp_list = list(projects.keys())
    random.shuffle(bp_list)
    
    count = 0
    for bp in bp_list:
        samples = projects[bp]
        if not samples:
            continue
        
        print(f"Processing project {bp}...")
        try:
            # Fetch study info
            study = build_study_metadata(bp)
            study_title = getattr(study.title, 'value', bp)
            study_desc = getattr(study.description, 'value', "")
            abstract = getattr(study.paper_abstract, 'value', None)
            
            # Pick up to 2 samples per project
            for s_data in samples[:2]:
                run_id = s_data["sample_id"]
                bs_id = s_data["biosample_id"]
                
                print(f"  Fetching sample {run_id} ({bs_id})...")
                # Fetch baseline sample info to get titles/descriptions
                sample_meta = build_sample_metadata(run_id, bp, biosample_id=bs_id)
                
                s_title = getattr(sample_meta.sample_title, 'value', None)
                s_desc = getattr(sample_meta.sample_description, 'value', None)
                
                # Fetch run metadata for additional context (library strategy, platform, etc.)
                try:
                    rmeta = build_run_metadata(run_id)
                    run_info = rmeta.model_dump()
                    # Convert BaseProvenance objects to values for readability in prompt
                    for k, v in run_info.items():
                        if isinstance(v, dict) and 'value' in v:
                            run_info[k] = v['value']
                except Exception as e:
                    print(f"    Warning: Could not fetch run metadata for {run_id}: {e}")
                    run_info = None

                # Use only the raw characteristics from the source (BioSample/GEO)
                chars = sample_meta.raw_characteristics or {}

                # Backfill missing ground truth from NCBI metadata if possible
                # This makes the benchmark fairer when RiboCrypt is missing info that is clearly in the raw context
                gt = s_data["ground_truth"].copy()
                
                # Helper to backfill from SampleMetadata fields
                def backfill(field_name, sm_attr):
                    if not gt.get(field_name) and sm_attr:
                        gt[field_name] = sm_attr.value.lower()

                backfill("organism", sample_meta.organism)
                backfill("tissue", sample_meta.tissue)
                backfill("strain", sample_meta.strain)
                backfill("cell_line", sample_meta.cell_line)
                backfill("cell_type", sample_meta.cell_type)
                backfill("age", sample_meta.age)
                backfill("sex", sample_meta.sex)
                
                # Also try to get organism from run_metadata if still missing
                if not gt.get("organism") and rmeta and getattr(rmeta, 'organism', None):
                    gt["organism"] = rmeta.organism.lower()

                gold = GoldStandardSample(
                    sample_id=run_id,
                    study_title=study_title,
                    study_description=study_desc,
                    sample_title=s_title,
                    sample_description=s_desc,
                    characteristics=chars,
                    run_metadata=run_info,
                    abstract=abstract,
                    ground_truth=gt,
                    bioproject_id=bp,
                    biosample_id=bs_id
                )
                selected_samples.append(gold)
            
            count += 1
            if count >= limit_projects:
                break
                
        except Exception as e:
            print(f"  Error processing {bp}: {e}")
            continue

    # Save to JSON
    output_data = {
        "samples": [as_dict(s) for s in selected_samples]
    }
    
    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    with open(OUTPUT_PATH, "w") as f:
        json.dump(output_data, f, indent=2)
    
    print(f"Created benchmark dataset with {len(selected_samples)} samples at {OUTPUT_PATH}")

def as_dict(obj):
    # Simple recursive dict conversion for dataclasses
    from dataclasses import is_dataclass, asdict
    if is_dataclass(obj):
        return asdict(obj)
    return obj

if __name__ == "__main__":
    prepare_benchmark()
