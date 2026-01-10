"""Tests for base schema definitions."""

import pytest
from datetime import datetime
from omics_extractor.schemas.base import (
    BaseProvenance,
    StudyMetadata,
    SampleMetadata,
    RunMetadata,
    RiboSeqRunMetadata,
)


class TestBaseProvenance:
    """Test BaseProvenance model."""

    def test_create_basic_provenance(self):
        """Should create provenance with required fields."""
        prov = BaseProvenance(
            value="Mus musculus",
            source="geo",
            source_id="GSE112882",
            confidence=0.95,
        )

        assert prov.value == "Mus musculus"
        assert prov.source == "geo"
        assert prov.confidence == 0.95

    def test_provenance_with_ontology(self):
        """Should store ontology mapping."""
        prov = BaseProvenance(
            value="brain",
            source="geo",
            source_id="GSE112882",
            confidence=0.90,
            ontology_term="UBERON:0000955",
            ontology_label="brain",
        )

        assert prov.has_ontology()
        assert prov.ontology_term == "UBERON:0000955"

    def test_confidence_validation(self):
        """Should validate confidence is between 0 and 1."""
        with pytest.raises(ValueError):
            BaseProvenance(
                value="test",
                source="geo",
                source_id="GSE123",
                confidence=1.5,  # Invalid
            )

    def test_high_confidence_check(self):
        """Should check if confidence is high."""
        high_conf = BaseProvenance(
            value="test",
            source="sra",
            source_id="SRR123",
            confidence=0.95,
        )
        low_conf = BaseProvenance(
            value="test",
            source="llm",
            source_id="GSE123",
            confidence=0.5,
        )

        assert high_conf.is_high_confidence()
        assert not low_conf.is_high_confidence()

    def test_extraction_metadata(self):
        """Should store extraction metadata."""
        now = datetime.now()
        prov = BaseProvenance(
            value="test",
            source="llm",
            source_id="GSE123",
            confidence=0.8,
            extracted_text="original text",
            extraction_method="gpt-4",
            extraction_timestamp=now,
            notes="Extracted from abstract",
        )

        assert prov.extracted_text == "original text"
        assert prov.extraction_method == "gpt-4"
        assert prov.extraction_timestamp == now


class TestStudyMetadata:
    """Test StudyMetadata model."""

    def test_create_minimal_study(self):
        """Should create study with required fields."""
        study = StudyMetadata(
            bioproject_id="PRJNA123456",
            title=BaseProvenance(
                value="Test Study",
                source="geo",
                source_id="GSE112882",
                confidence=1.0,
            ),
            description=BaseProvenance(
                value="Study description",
                source="geo",
                source_id="GSE112882",
                confidence=1.0,
            ),
            organism=BaseProvenance(
                value="Homo sapiens",
                source="geo",
                source_id="GSE112882",
                confidence=1.0,
                ontology_term="NCBITaxon:9606",
            ),
        )

        assert study.bioproject_id == "PRJNA123456"
        assert study.organism.ontology_term == "NCBITaxon:9606"

    def test_study_with_publication(self):
        """Should include publication metadata."""
        study = StudyMetadata(
            bioproject_id="PRJNA123456",
            title=BaseProvenance(
                value="Test",
                source="geo",
                source_id="GSE123",
                confidence=1.0,
            ),
            description=BaseProvenance(
                value="Desc",
                source="geo",
                source_id="GSE123",
                confidence=1.0,
            ),
            organism=BaseProvenance(
                value="Mouse",
                source="geo",
                source_id="GSE123",
                confidence=1.0,
            ),
            pmid="29618526",
            doi=BaseProvenance(
                value="10.1038/nature12345",
                source="pubmed",
                source_id="29618526",
                confidence=1.0,
            ),
            authors=["Smith J", "Jones A"],
        )

        assert study.pmid == "29618526"
        assert len(study.authors) == 2


class TestSampleMetadata:
    """Test SampleMetadata model."""

    def test_create_minimal_sample(self):
        """Should create sample with required fields."""
        sample = SampleMetadata(
            sample_id="SRS123456",
            bioproject_id="PRJNA123456",
            organism=BaseProvenance(
                value="Mus musculus",
                source="sra",
                source_id="SRS123456",
                confidence=1.0,
                ontology_term="NCBITaxon:10090",
            ),
        )

        assert sample.sample_id == "SRS123456"
        assert sample.organism.has_ontology()

    def test_sample_with_biological_context(self):
        """Should include biological metadata."""
        sample = SampleMetadata(
            sample_id="SRS123456",
            bioproject_id="PRJNA123456",
            organism=BaseProvenance(
                value="Mouse",
                source="geo",
                source_id="GSM123",
                confidence=1.0,
            ),
            tissue=BaseProvenance(
                value="brain",
                source="geo",
                source_id="GSM123",
                confidence=0.9,
                ontology_term="UBERON:0000955",
            ),
            cell_type=BaseProvenance(
                value="neuron",
                source="geo",
                source_id="GSM123",
                confidence=0.85,
                ontology_term="CL:0000540",
            ),
            age=BaseProvenance(
                value="8 weeks",
                source="geo",
                source_id="GSM123",
                confidence=1.0,
            ),
        )

        assert sample.tissue.value == "brain"
        assert sample.tissue.has_ontology()
        assert sample.age.value == "8 weeks"

    def test_sample_with_experimental_conditions(self):
        """Should include experimental metadata."""
        sample = SampleMetadata(
            sample_id="SRS123456",
            bioproject_id="PRJNA123456",
            organism=BaseProvenance(
                value="Mouse",
                source="geo",
                source_id="GSM123",
                confidence=1.0,
            ),
            treatment=BaseProvenance(
                value="doxycycline",
                source="geo",
                source_id="GSM123",
                confidence=0.9,
            ),
            timepoint=BaseProvenance(
                value="24h",
                source="geo",
                source_id="GSM123",
                confidence=1.0,
            ),
            replicate=BaseProvenance(
                value="1",
                source="geo",
                source_id="GSM123",
                confidence=1.0,
            ),
        )

        assert sample.treatment.value == "doxycycline"
        assert sample.timepoint.value == "24h"


class TestRunMetadata:
    """Test RunMetadata model."""

    def test_create_basic_run(self):
        """Should create run with required fields."""
        run = RunMetadata(
            run_id="SRR123456",
            sample_id="SRS123456",
            bioproject_id="PRJNA123456",
            experiment_id="SRX123456",
            library_strategy=BaseProvenance(
                value="RNA-Seq",
                source="sra",
                source_id="SRR123456",
                confidence=1.0,
            ),
            library_source="TRANSCRIPTOMIC",
            library_selection="cDNA",
            library_layout="PAIRED",
            platform="ILLUMINA",
        )

        assert run.run_id == "SRR123456"
        assert run.library_strategy.value == "RNA-Seq"
        assert run.platform == "ILLUMINA"

    def test_run_with_statistics(self):
        """Should include sequencing statistics."""
        run = RunMetadata(
            run_id="SRR123456",
            sample_id="SRS123456",
            bioproject_id="PRJNA123456",
            experiment_id="SRX123456",
            library_strategy=BaseProvenance(
                value="RNA-Seq",
                source="sra",
                source_id="SRR123456",
                confidence=1.0,
            ),
            library_source="TRANSCRIPTOMIC",
            library_selection="cDNA",
            library_layout="SINGLE",
            platform="ILLUMINA",
            read_count=10000000,
            base_count=1500000000,
            avg_length=150.0,
        )

        assert run.read_count == 10000000
        assert run.avg_length == 150.0


class TestRiboSeqRunMetadata:
    """Test RiboSeqRunMetadata model."""

    def test_create_riboseq_run(self):
        """Should create RiboSeq run with protocol details."""
        run = RiboSeqRunMetadata(
            run_id="SRR123456",
            sample_id="SRS123456",
            bioproject_id="PRJNA123456",
            experiment_id="SRX123456",
            library_strategy=BaseProvenance(
                value="RiboSeq",
                source="sra",
                source_id="SRR123456",
                confidence=1.0,
            ),
            library_source="TRANSCRIPTOMIC",
            library_selection="other",
            library_layout="SINGLE",
            platform="ILLUMINA",
            inhibitor=BaseProvenance(
                value="cycloheximide",
                source="geo",
                source_id="GSM123",
                confidence=0.9,
                ontology_term="CHEBI:27641",
            ),
            nuclease=BaseProvenance(
                value="RNase I",
                source="geo",
                source_id="GSM123",
                confidence=0.85,
            ),
            footprint_min_length=25,
            footprint_max_length=35,
        )

        assert run.inhibitor.value == "cycloheximide"
        assert run.inhibitor.has_ontology()
        assert run.nuclease.value == "RNase I"
        assert run.footprint_min_length == 25
