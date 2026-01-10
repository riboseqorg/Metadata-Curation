#!/usr/bin/env python3
"""
Quick test of traceable output format on 5 RiboSeq studies.
"""

import json
import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from omics_extractor.extraction.builder import build_project_metadata
from omics_extractor.output.traceable_format import create_traceable_report


def main():
    """Test traceable format on a few projects."""
    test_projects = [
        "PRJNA1170270",  # Known working project
        "PRJDB10544",
        "PRJDB10799",
        "PRJEB12126",
        "PRJNA1002596",
    ]

    output_dir = Path("test_traceable_output")
    output_dir.mkdir(exist_ok=True)

    for i, project_id in enumerate(test_projects, 1):
        print(f"\n[{i}/{len(test_projects)}] Processing {project_id}...")

        try:
            # Extract metadata
            metadata = build_project_metadata(project_id)

            # Create traceable report
            report = create_traceable_report(
                study=metadata["study"],
                samples=metadata["samples"],
                runs=metadata["runs"],
                include_raw=True,
                include_statistics=True,
            )

            # Save
            output_file = output_dir / f"{project_id}_traceable.json"
            with open(output_file, "w") as f:
                json.dump(report, f, indent=2, default=str)

            print(f"  ✓ {len(metadata['samples'])} samples, {len(metadata['runs'])} runs")
            print(f"  ✓ Saved to: {output_file}")

            # Show example structure from first project
            if i == 1:
                print(f"\n  Example structure:")
                print(f"    Top level keys: {list(report.keys())}")

                if report["samples"]:
                    first_sample_id = list(report["samples"].keys())[0]
                    first_sample = report["samples"][first_sample_id]
                    print(f"\n    Sample structure: {list(first_sample.keys())}")
                    print(f"\n    Quick view: {first_sample['quick_view']}")

                    if "tissue" in first_sample.get("biological_metadata", {}):
                        tissue_field = first_sample["biological_metadata"]["tissue"]
                        print(f"\n    Tissue field structure:")
                        print(f"      - value: {tissue_field['value']}")
                        print(f"      - provenance keys: {list(tissue_field['provenance'].keys())}")
                        if "extraction_details" in tissue_field:
                            print(f"      - extraction_details keys: {list(tissue_field['extraction_details'].keys())}")

        except Exception as e:
            print(f"  ✗ Error: {e}")
            import traceback
            traceback.print_exc()

    print(f"\n{'='*60}")
    print(f"Test complete! Check {output_dir}/ for outputs")


if __name__ == "__main__":
    main()
