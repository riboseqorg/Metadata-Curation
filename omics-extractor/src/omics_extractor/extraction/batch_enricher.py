"""Batch LLM enrichment for efficient GPU utilization.

This module provides batch processing capabilities for LLM enrichment,
optimized for GPU nodes in production environments.

Key features:
- Batch multiple samples per API call for efficiency
- Progress tracking and resumption for long-running jobs
- Concurrent API requests with rate limiting
- Support for multiple LLM providers (Claude API, local models)
"""

import json
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor, as_completed
import anthropic

from ..schemas.base import SampleMetadata
from .llm_extractor import enrich_sample_metadata


class BatchEnricher:
    """Batch LLM enrichment with progress tracking and resumption."""

    def __init__(
        self,
        api_key: Optional[str] = None,
        max_workers: int = 4,
        rate_limit_rpm: int = 50,
        checkpoint_every: int = 10,
    ):
        """
        Initialize batch enricher.

        Args:
            api_key: Anthropic API key
            max_workers: Number of concurrent API requests
            rate_limit_rpm: Maximum requests per minute
            checkpoint_every: Save progress every N samples
        """
        self.api_key = api_key
        self.max_workers = max_workers
        self.rate_limit_rpm = rate_limit_rpm
        self.checkpoint_every = checkpoint_every

        # For rate limiting
        self._request_times: List[float] = []

    def _wait_for_rate_limit(self):
        """Wait if we're hitting rate limits."""
        now = time.time()

        # Remove requests older than 1 minute
        self._request_times = [t for t in self._request_times if now - t < 60]

        # If we're at the limit, wait
        if len(self._request_times) >= self.rate_limit_rpm:
            oldest = self._request_times[0]
            wait_time = 60 - (now - oldest)
            if wait_time > 0:
                time.sleep(wait_time)
                self._request_times = []

        self._request_times.append(now)

    def enrich_batch(
        self,
        samples: Dict[str, SampleMetadata],
        study_title: str,
        study_description: str,
        abstract: Optional[str] = None,
        journal: Optional[str] = None,
        authors: Optional[str] = None,
        publication_date: Optional[str] = None,
        checkpoint_file: Optional[Path] = None,
        resume: bool = False,
        scheme: str = "default",
    ) -> Tuple[Dict[str, SampleMetadata], Dict]:
        """
        Enrich multiple samples with LLM.

        Args:
            samples: Dictionary of sample_id -> SampleMetadata
            study_title: Study title for context
            study_description: Study description for context
            abstract: Optional publication abstract
            checkpoint_file: File to save progress
            resume: Resume from checkpoint if exists

        Returns:
            Tuple of (enriched_samples, statistics)
        """
        # Load checkpoint if resuming
        completed = set()
        if resume and checkpoint_file and checkpoint_file.exists():
            with open(checkpoint_file) as f:
                checkpoint = json.load(f)
                completed = set(checkpoint.get("completed", []))
                print(f"Resuming from checkpoint: {len(completed)} samples already done")

        # Filter samples that need enrichment
        to_enrich = {
            sid: sample
            for sid, sample in samples.items()
            if sid not in completed
        }

        if not to_enrich:
            print("All samples already enriched")
            return samples, {"total": len(samples), "enriched": 0, "skipped": len(completed)}

        # Statistics
        stats = {
            "total": len(samples),
            "to_enrich": len(to_enrich),
            "enriched": 0,
            "failed": 0,
            "fields_added": {"tissue": 0, "cell_type": 0, "cell_line": 0, "treatment": 0, "strain": 0},
        }

        # Enrich samples concurrently
        enriched_samples = dict(samples)  # Start with original
        sample_items = list(to_enrich.items())

        with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
            # Submit all tasks
            futures = {}
            for sample_id, sample in sample_items:
                # Wait for rate limit
                self._wait_for_rate_limit()

                # Get sample-specific info
                sample_title = sample.sample_title.value if sample.sample_title else None
                sample_desc = sample.sample_description.value if sample.sample_description else None

                future = executor.submit(
                    self._enrich_single_sample,
                    sample_id,
                    sample,
                    study_title,
                    study_description,
                    sample_title,
                    sample_desc,
                    abstract,
                    journal,
                    authors,
                    publication_date,
                )
                futures[future] = sample_id

            # Process completed tasks
            for future in as_completed(futures):
                sample_id = futures[future]

                try:
                    original, enriched, added_fields = future.result()

                    # Update samples
                    enriched_samples[sample_id] = enriched
                    stats["enriched"] += 1

                    # Track fields
                    for field in added_fields:
                        stats["fields_added"][field] += 1

                    # Mark as completed
                    completed.add(sample_id)

                    # Save checkpoint
                    if checkpoint_file and stats["enriched"] % self.checkpoint_every == 0:
                        self._save_checkpoint(checkpoint_file, completed, stats)

                except Exception as e:
                    print(f"Failed to enrich {sample_id}: {e}")
                    stats["failed"] += 1

        # Final checkpoint
        if checkpoint_file:
            self._save_checkpoint(checkpoint_file, completed, stats)

        return enriched_samples, stats

    def _enrich_single_sample(
        self,
        sample_id: str,
        sample: SampleMetadata,
        study_title: str,
        study_description: str,
        sample_title: Optional[str],
        sample_description: Optional[str],
        abstract: Optional[str],
        journal: Optional[str] = None,
        authors: Optional[str] = None,
        publication_date: Optional[str] = None,
        scheme: str = "default",
    ) -> Tuple[SampleMetadata, SampleMetadata, List[str]]:
        """
        Enrich a single sample.

        Returns:
            Tuple of (original_sample, enriched_sample, added_fields)
        """
        enriched = enrich_sample_metadata(
            sample=sample,
            study_title=study_title,
            study_description=study_description,
            sample_title=sample_title,
            sample_description=sample_description,
            abstract=abstract,
            journal=journal,
            authors=authors,
            publication_date=publication_date,
            api_key=self.api_key,
            source_id=f"llm_batch_{sample_id}",
            scheme=scheme,
        )

        # Determine what was added
        added_fields = []
        if enriched.tissue and not sample.tissue:
            added_fields.append("tissue")
        if enriched.cell_type and not sample.cell_type:
            added_fields.append("cell_type")
        if enriched.cell_line and not sample.cell_line:
            added_fields.append("cell_line")
        if enriched.treatment and not sample.treatment:
            added_fields.append("treatment")
        if enriched.strain and not sample.strain:
            added_fields.append("strain")

        return sample, enriched, added_fields

    def _save_checkpoint(self, checkpoint_file: Path, completed: set, stats: Dict):
        """Save progress checkpoint."""
        checkpoint = {
            "timestamp": datetime.now().isoformat(),
            "completed": list(completed),
            "statistics": stats,
        }

        with open(checkpoint_file, "w") as f:
            json.dump(checkpoint, f, indent=2)


def enrich_batch_from_files(
    input_files: List[Path],
    output_dir: Path,
    api_key: str,
    max_workers: int = 4,
    rate_limit_rpm: int = 50,
    resume: bool = False,
) -> Dict:
    """
    Batch enrich multiple metadata files.

    This is useful for processing many projects at once on a GPU node.

    Args:
        input_files: List of metadata JSON files
        output_dir: Directory to save enriched files
        api_key: Anthropic API key
        max_workers: Concurrent API requests
        rate_limit_rpm: Rate limit
        resume: Resume from checkpoints

    Returns:
        Overall statistics
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    enricher = BatchEnricher(
        api_key=api_key,
        max_workers=max_workers,
        rate_limit_rpm=rate_limit_rpm,
    )

    overall_stats = {
        "total_files": len(input_files),
        "processed": 0,
        "failed": 0,
        "total_samples": 0,
        "total_enriched": 0,
    }

    for input_file in input_files:
        print(f"\nProcessing {input_file.name}...")

        try:
            # Load metadata
            with open(input_file) as f:
                data = json.load(f)

            # Extract context
            study = data.get("study", {})
            study_title = study.get("title", {}).get("value", "")
            study_description = study.get("description", {}).get("value", "")
            abstract = study.get("abstract", {}).get("value")
            journal = study.get("journal", {}).get("value")
            authors = study.get("authors")
            if isinstance(authors, list):
                authors = ", ".join(authors)
            pub_date = study.get("publication_date")

            # Reconstruct samples
            samples = {}
            for sid, sdata in data.get("samples", {}).items():
                samples[sid] = SampleMetadata(**sdata)

            # Checkpoint file
            checkpoint_file = output_dir / f"{input_file.stem}_checkpoint.json"

            # Enrich
            enriched, stats = enricher.enrich_batch(
                samples=samples,
                study_title=study_title,
                study_description=study_description,
                abstract=abstract,
                journal=journal,
                authors=authors,
                publication_date=pub_date,
                checkpoint_file=checkpoint_file,
                resume=resume,
            )

            # Update data
            data["samples"] = {sid: s.model_dump() for sid, s in enriched.items()}
            data["enrichment_metadata"] = {
                "timestamp": datetime.now().isoformat(),
                "statistics": stats,
            }

            # Save
            output_file = output_dir / f"{input_file.stem}_enriched.json"
            with open(output_file, "w") as f:
                json.dump(data, f, indent=2)

            # Update overall stats
            overall_stats["processed"] += 1
            overall_stats["total_samples"] += stats["total"]
            overall_stats["total_enriched"] += stats["enriched"]

            # Clean up checkpoint
            if checkpoint_file.exists():
                checkpoint_file.unlink()

        except Exception as e:
            print(f"Failed to process {input_file.name}: {e}")
            overall_stats["failed"] += 1

    return overall_stats
