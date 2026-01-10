"""Tests for SRA metadata fetcher."""

import pytest
from unittest.mock import Mock, patch
from omics_extractor.fetchers.sra import (
    fetch_project_runs,
    fetch_run_details,
    SRARunMetadata,
)


class TestFetchProjectRuns:
    """Test fetching all runs for a BioProject."""

    def test_fetch_project_runs_returns_list(self):
        """Should return list of run IDs for valid project."""
        # Use a small, stable test project
        # PRJNA401384 - small RiboSeq project with known runs
        runs = fetch_project_runs("PRJNA401384")

        assert isinstance(runs, list)
        assert len(runs) > 0
        assert all(run.startswith("SRR") or run.startswith("ERR") or run.startswith("DRR") for run in runs)

    def test_fetch_project_runs_handles_invalid_project(self):
        """Should handle invalid project IDs gracefully."""
        runs = fetch_project_runs("PRJNA999999999")
        assert runs == []

    def test_fetch_project_runs_accepts_gse_id(self):
        """Should accept GSE IDs and convert to BioProject."""
        # GSE112882 corresponds to PRJNA401384
        runs = fetch_project_runs("GSE112882")
        assert isinstance(runs, list)


class TestFetchRunDetails:
    """Test fetching detailed metadata for a single run."""

    def test_fetch_run_details_returns_metadata(self):
        """Should return structured metadata for valid run."""
        # SRR6426586 - known RiboSeq run from PRJNA401384
        metadata = fetch_run_details("SRR6426586")

        assert isinstance(metadata, SRARunMetadata)
        assert metadata.run_id == "SRR6426586"
        assert metadata.library_strategy is not None
        assert metadata.library_source is not None
        assert metadata.platform is not None

    def test_fetch_run_details_includes_sample_info(self):
        """Should include sample and experiment metadata."""
        metadata = fetch_run_details("SRR6426586")

        assert metadata.sample_id is not None
        assert metadata.experiment_id is not None
        assert metadata.bioproject_id is not None

    def test_fetch_run_details_handles_invalid_run(self):
        """Should raise error for invalid run IDs."""
        with pytest.raises(ValueError, match="Run .* not found"):
            fetch_run_details("SRR999999999999")


class TestSRARunMetadata:
    """Test the SRARunMetadata data model."""

    def test_metadata_model_validation(self):
        """Should validate required fields."""
        metadata = SRARunMetadata(
            run_id="SRR123456",
            experiment_id="SRX123456",
            sample_id="SRS123456",
            bioproject_id="PRJNA123456",
            library_strategy="RNA-Seq",
            library_source="TRANSCRIPTOMIC",
            library_selection="cDNA",
            platform="ILLUMINA",
        )

        assert metadata.run_id == "SRR123456"
        assert metadata.library_strategy == "RNA-Seq"

    def test_metadata_model_optional_fields(self):
        """Should allow optional fields."""
        metadata = SRARunMetadata(
            run_id="SRR123456",
            experiment_id="SRX123456",
            sample_id="SRS123456",
            bioproject_id="PRJNA123456",
            library_strategy="RNA-Seq",
            library_source="TRANSCRIPTOMIC",
            library_selection="cDNA",
            platform="ILLUMINA",
            instrument_model="Illumina HiSeq 2500",
            read_count=10000000,
            base_count=1000000000,
        )

        assert metadata.instrument_model == "Illumina HiSeq 2500"
        assert metadata.read_count == 10000000
