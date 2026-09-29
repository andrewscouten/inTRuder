"""Reading reference FASTA, including the header and case variation it arrives with."""

from __future__ import annotations

import gzip

import pytest

from intruder.trcore.io import fetch_region, open_indexed, read_fasta

pytest.importorskip("pyfastx")

HEADER = ">NZ_CP012345.1 Burkholderia cepacia chromosome 1, complete sequence"
BASES = "ATGCGGCATCAGCAGCAGCAGCAGTTACGG"


def write_fasta(path, text):
    path.write_text(text)
    return path


@pytest.fixture
def genome(tmp_path):
    return write_fasta(
        tmp_path / "genome.fna",
        f"{HEADER}\n{BASES}\nATCCGGGTTT\n>contig2 a second one\nAAAACCCCGGGGTTTT\n",
    )


def test_reads_every_contig(genome):
    records = list(read_fasta(genome))
    assert [r.name for r in records] == ["NZ_CP012345.1", "contig2"]


def test_wrapped_lines_are_joined(genome):
    first = next(iter(read_fasta(genome)))
    assert first.sequence == BASES + "ATCCGGGTTT"


def test_header_splits_into_name_and_description(genome):
    """The organism string lives in the description and is carried through unparsed."""
    first = next(iter(read_fasta(genome)))
    assert first.name == "NZ_CP012345.1"
    assert first.description == "Burkholderia cepacia chromosome 1, complete sequence"


def test_header_with_no_description(tmp_path):
    path = write_fasta(tmp_path / "bare.fna", ">contig1\nACGT\n")
    assert next(iter(read_fasta(path))).description == ""


def test_soft_masked_bases_are_upper_cased(tmp_path):
    """Motif comparison is case-sensitive, so lower case would make acg != ACG."""
    path = write_fasta(tmp_path / "masked.fna", ">c1\nacgtACGTacgt\n")
    assert next(iter(read_fasta(path))).sequence == "ACGTACGTACGT"


def test_case_is_preserved_when_asked(tmp_path):
    path = write_fasta(tmp_path / "masked.fna", ">c1\nacgtACGT\n")
    record = next(iter(read_fasta(path, uppercase=False)))
    assert record.sequence == "acgtACGT"


def test_reads_gzipped_fasta(tmp_path, genome):
    path = tmp_path / "genome.fna.gz"
    path.write_bytes(gzip.compress(genome.read_bytes()))
    assert [r.name for r in read_fasta(path)] == ["NZ_CP012345.1", "contig2"]


def test_origin_is_genomic_and_names_its_contig(genome):
    first = next(iter(read_fasta(genome)))
    assert first.origin.in_reference is True
    assert first.origin.locate(2, 6) == ("NZ_CP012345.1", 2, 6)


def test_streaming_writes_no_index_beside_the_input(genome):
    """The genomes this scans are checksummed against a manifest lock."""
    list(read_fasta(genome))
    assert not (genome.parent / "genome.fna.fxi").exists()


def test_fetch_region_offsets_are_absolute(genome):
    fasta = open_indexed(genome)
    window = fetch_region(fasta, "NZ_CP012345.1", 4, 12)
    assert window.sequence == BASES[4:12]
    assert window.origin.locate(0, 3) == ("NZ_CP012345.1", 4, 7)
