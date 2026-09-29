"""The two output schemas, and why a call cannot be written in either one."""

from __future__ import annotations

import csv

import pytest

from intruder.pipeline.trf.find import RepeatCall
from intruder.pipeline.trf.writers import (
    BED_COLUMNS,
    GENOME_COLUMNS,
    INSERTION_COLUMNS,
    genome_row,
    write_bed,
    write_tsv,
)
from intruder.trcore.io import Genomic, Inserted


def call(origin=None, **overrides):
    base = {"name": "contig_1", "origin": origin or Genomic("contig_1"),
            "start": 100, "end": 190, "motif": "CAG", "purity": 0.912}
    return RepeatCall(**{**base, **overrides})


def test_genome_row_carries_the_reference_interval():
    row = genome_row(call(), genome="GCF_000009045")
    assert (row["contig"], row["start"], row["end"]) == ("contig_1", 100, 190)
    assert row["genome"] == "GCF_000009045"
    assert set(row) == set(GENOME_COLUMNS)


def test_genome_row_refuses_an_insertion_call():
    """An insertion has no interval, so a genome row would index nothing."""
    inserted = call(origin=Inserted(contig="chr1", position=99, length=90))
    with pytest.raises(ValueError, match="no reference interval"):
        genome_row(inserted, genome="x")


def test_genome_row_reports_window_offsets_on_the_contig():
    windowed = call(origin=Genomic("chr1", start=1_000))
    row = genome_row(windowed, genome="x")
    assert (row["start"], row["end"]) == (1_100, 1_190)


def test_the_insertion_schema_is_unchanged():
    """novelty, strchive and the isolation forest all read stage 01's columns."""
    assert INSERTION_COLUMNS == (
        "chrom", "ins_coord", "SVID", "depth", "insert_size", "sample", "allele",
        "rep_start", "rep_end", "motif", "purity", "motif_length", "rep_length",
        "rep_units",
    )


def test_the_two_schemas_do_not_share_a_coordinate_column():
    """`start` must not mean a contig position in one table and an offset in another."""
    assert {"start", "end"} & set(INSERTION_COLUMNS) == set()
    assert {"rep_start", "rep_end"} & set(GENOME_COLUMNS) == set()


def test_purity_is_rounded_for_output():
    assert genome_row(call(purity=0.9123456), genome="x")["purity"] == 0.912


def test_write_tsv_round_trips(tmp_path):
    path = tmp_path / "calls.tsv"
    rows = [genome_row(call(), genome="g")]
    assert write_tsv(path, GENOME_COLUMNS, rows) == 1
    with open(path, newline="") as handle:
        read_back = list(csv.DictReader(handle, delimiter="\t"))
    assert read_back[0]["motif"] == "CAG"
    assert list(read_back[0]) == list(GENOME_COLUMNS)


def test_write_tsv_appends_extra_columns(tmp_path):
    path = tmp_path / "calls.tsv"
    rows = [{**genome_row(call(), genome="g"), "filter_status": "PASS"}]
    write_tsv(path, GENOME_COLUMNS, rows, extras=["filter_status"])
    header = path.read_text().splitlines()[0].split("\t")
    assert header[-1] == "filter_status"


def test_write_bed_is_headerless_bed4(tmp_path):
    path = tmp_path / "loci.bed"
    assert write_bed(path, [genome_row(call(), genome="g")]) == 1
    line = path.read_text().splitlines()
    assert len(line) == 1
    assert line[0].split("\t") == ["contig_1", "100", "190", "CAG"]
    assert BED_COLUMNS == ("contig", "start", "end", "motif")


def test_a_vcf_without_read_support_tags_is_an_ordinary_input():
    """`DR`/`DV` are Sniffles' tags, not the VCF spec's.

    cyvcf2 raises a bare `KeyError: b'DV'` from inside its header lookup, which
    came out as a traceback from the middle of a whole-VCF scan.
    """
    from intruder.pipeline.trf.writers import _depth, _format_field

    class NoTags:
        def format(self, tag):
            raise KeyError(tag.encode())

    assert _format_field(NoTags(), "DV") is None
    # Empty, not "0,0": zero is a measurement, and this is the absence of one.
    assert _depth(None, None, 0) == ""
    assert _depth([[15]], [[20]], 0) == "15,20"


# htslib's missing-integer sentinel, which is what cyvcf2 hands back for a
# sample whose DR:DV is `.` in a VCF that declares both tags.
BCF_INT32_MISSING = -2_147_483_648


def test_a_sample_with_no_depth_recorded_is_empty_not_a_sentinel():
    """`DR:DV` of `.` is the absence of a measurement, like a VCF without the tags.

    Written through, the column read `-2147483648,-2147483648`, which
    `as_total` sums to a large negative depth: not an unreadable cell a
    threshold declines to judge, but a measurement failing every `--min-depth`.
    """
    from intruder.pipeline.trf.writers import _depth

    missing = [[BCF_INT32_MISSING]]
    assert _depth(missing, missing, 0) == ""
    assert _depth(missing, [[20]], 0) == ""        # half a pair is not an answer
    assert _depth([[15]], missing, 0) == ""
    assert _depth([[0]], [[0]], 0) == "0,0"        # a measured zero still is one


def test_a_sentinel_depth_would_not_survive_the_reader():
    """Why the empty cell matters: the pair the sentinel writes parses to a depth."""
    from intruder.trcore.utils.parse import as_total

    assert as_total(f"{BCF_INT32_MISSING},{BCF_INT32_MISSING}") < 0
    assert as_total("") is None


# --------------------------------------------------------------------------- #
# genotypes
# --------------------------------------------------------------------------- #

def test_alleles_drops_the_phase_flag_whatever_the_ploidy():
    """cyvcf2 appends a phased bool, and `True == 1` in Python.

    A haploid `1` arrives as `[1, True]`; reading the first two entries counted
    that bool as a second copy of allele 1, so a haploid carrier emitted two
    rows and wrote `True` into the numeric `allele` column. Every chrX and chrY
    call in a male sample is haploid.
    """
    from intruder.pipeline.trf.writers import _alleles

    assert _alleles([1, True]) == [1]              # haploid  1
    assert _alleles([0, 1, False]) == [0, 1]       # diploid  0/1
    assert _alleles([1, 0, True]) == [1, 0]        # phased   1|0
    assert _alleles([-1, -1, False]) == [-1, -1]   # no call  ./.
    assert _alleles([1, 1, 1, False]) == [1, 1, 1]  # triploid
