"""Integration tests using real RiboSeq BioProjects from the curated dataset."""

import pytest
from omics_extractor.fetchers.sra import fetch_project_runs, fetch_run_details
from omics_extractor.fetchers.geo import fetch_project_metadata, fetch_sample_metadata


# Real RiboSeq projects from ../bioprojects.csv
REAL_PROJECTS = {
    "PRJNA1170270": {
        "expected_runs_min": 1,
        "description": "Ribo-seq of mouse spermatocytes",
    },
    "PRJNA1176138": {
        "expected_runs_min": 20,
        "description": "Cytosolic Ribosomal Protein Haploinsufficiency",
    },
}


class TestRealDataIntegration:
    """Test fetchers with real curated RiboSeq projects."""

    @pytest.mark.parametrize("project_id,expected", REAL_PROJECTS.items())
    def test_fetch_real_project_runs(self, project_id, expected):
        """Should fetch runs from real RiboSeq projects."""
        runs = fetch_project_runs(project_id)

        assert len(runs) >= expected["expected_runs_min"]
        assert all(
            run.startswith("SRR") or run.startswith("ERR") or run.startswith("DRR")
            for run in runs
        )

    def test_fetch_real_run_details(self):
        """Should fetch detailed metadata for a known RiboSeq run."""
        # PRJNA1170270 - Ribo-seq of mouse spermatocytes
        # SRR30910340 - first run from this project
        metadata = fetch_run_details("SRR30910340")

        assert metadata.run_id == "SRR30910340"
        assert metadata.bioproject_id == "PRJNA1170270"
        assert metadata.library_strategy in ["Ribo-seq", "RiboSeq", "OTHER"]
        assert metadata.platform == "ILLUMINA"
        assert metadata.sample_id is not None

    def test_library_strategy_detection(self):
        """Should detect various library strategies in curated dataset."""
        # Test projects with known library strategies
        test_cases = [
            ("PRJNA1170270", ["Ribo-seq", "RiboSeq"]),  # Correctly labeled
            ("PRJNA1176138", ["OTHER", "Ribo-seq"]),  # Mislabeled as OTHER
            ("PRJNA1177760", ["RNA-Seq", "Ribo-seq"]),  # Mislabeled as RNA-Seq
        ]

        for project_id, acceptable_strategies in test_cases:
            runs = fetch_project_runs(project_id)
            if runs:
                metadata = fetch_run_details(runs[0])
                # Check if library strategy is one of the acceptable values
                # (some projects are mislabeled in SRA)
                assert metadata.library_strategy in acceptable_strategies, (
                    f"{project_id} has unexpected library strategy: "
                    f"{metadata.library_strategy}"
                )


class TestGEOIntegration:
    """Test GEO fetcher with real data."""

    def test_fetch_geo_series_with_bioproject(self):
        """Should fetch GEO series and link to BioProject."""
        # GSE112882 corresponds to PRJNA401384 (from our earlier manual test)
        metadata = fetch_project_metadata("GSE112882")

        assert metadata.gse_id == "GSE112882"
        assert metadata.title is not None
        assert len(metadata.title) > 0
        assert metadata.organism is not None


class TestCrossReference:
    """Test cross-referencing between SRA and GEO."""

    def test_bioproject_to_geo_series(self):
        """Should be able to cross-reference BioProject to GEO."""
        # GSE112882 has BioProject PRJNA449378 (not PRJNA401384)
        # Fetch GEO series
        geo_metadata = fetch_project_metadata("GSE112882")
        assert geo_metadata.bioproject_id == "PRJNA449378"

        # Verify we can fetch runs from that BioProject
        runs = fetch_project_runs("PRJNA449378")
        assert len(runs) > 0

        # Get run details to verify linkage
        run_metadata = fetch_run_details(runs[0])
        assert run_metadata.bioproject_id == "PRJNA449378"

    def test_sample_metadata_linkage(self):
        """Should link SRA sample to GEO sample."""
        # Get a run from a project
        runs = fetch_project_runs("PRJNA401384")
        run_metadata = fetch_run_details(runs[0])

        # Sample ID should exist
        assert run_metadata.sample_id is not None
        assert run_metadata.sample_id.startswith("SRS")


class TestDataQuality:
    """Test data quality and completeness."""

    def test_run_has_essential_metadata(self):
        """All runs should have essential metadata fields."""
        runs = fetch_project_runs("PRJNA1170270")
        metadata = fetch_run_details(runs[0])

        # Essential fields should be present
        assert metadata.run_id is not None
        assert metadata.experiment_id is not None
        assert metadata.sample_id is not None
        assert metadata.bioproject_id is not None
        assert metadata.library_strategy is not None
        assert metadata.library_source is not None
        assert metadata.platform is not None

    def test_project_metadata_completeness(self):
        """Project metadata should have key fields."""
        metadata = fetch_project_metadata("GSE112882")

        # Key fields
        assert metadata.gse_id is not None
        assert metadata.title is not None
        assert metadata.summary is not None
        assert metadata.organism is not None

        # Optional but commonly present
        # Note: not all projects have these
        if metadata.pubmed_id:
            assert len(metadata.pubmed_id) > 0
