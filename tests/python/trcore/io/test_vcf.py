"""Reading insertion sequence out of an SV VCF."""

from __future__ import annotations

import pytest

from intruder.trcore.io import read_insertions

pytest.importorskip("cyvcf2")

HEADER = """\
##fileformat=VCFv4.2
##contig=<ID=chr1,length=100000>
##ALT=<ID=INS,Description="Insertion">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of structural variant">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
"""

INSERT = "CAG" * 20          # 60 bp, over the default 50 bp floor


def write_vcf(tmp_path, *rows):
    path = tmp_path / "sv.vcf"
    path.write_text(HEADER + "".join(row + "\n" for row in rows))
    return path


def row(pos=100, vid="INS0", ref="A", alt=None, svtype="INS", svlen=60):
    alt = ref + INSERT if alt is None else alt
    info = f"SVTYPE={svtype}" + (f";SVLEN={svlen}" if svlen is not None else "")
    return f"chr1\t{pos}\t{vid}\t{ref}\t{alt}\t.\tPASS\t{info}"


def test_padding_base_is_stripped(tmp_path):
    """ALT repeats the REF base at POS; the inserted sequence is what follows it."""
    record = next(iter(read_insertions(write_vcf(tmp_path, row()))))
    assert record.sequence == INSERT


def test_alt_without_padding_is_left_alone(tmp_path):
    path = write_vcf(tmp_path, row(ref="A", alt="G" + INSERT[1:]))
    assert next(iter(read_insertions(path))).sequence == "G" + INSERT[1:]


def test_position_is_zero_based(tmp_path):
    """VCF POS 100 is the 1-based anchor, so the 0-based base it follows is 99."""
    record = next(iter(read_insertions(write_vcf(tmp_path, row(pos=100)))))
    assert record.origin.anchor() == ("chr1", 99)


def test_origin_carries_a_length_and_claims_no_reference_interval(tmp_path):
    record = next(iter(read_insertions(write_vcf(tmp_path, row()))))
    assert record.origin.length == 60
    assert record.origin.in_reference is False
    assert record.origin.locate(0, 10) is None


def test_origin_length_measures_the_alt_not_the_declared_svlen(tmp_path):
    """A merged callset's SVLEN is a consensus and need not describe its ALT.

    Every offset in this package indexes ALT, so a length that describes
    something else is not a bound those offsets fall inside -- and it is the
    denominator `novelty.insertions` and `pipeline.trf.filters` divide repeat
    coverage by. SURVIVOR's merged VCF in `data/sv_output` disagrees with its
    own ALT on 297 of 481 insertions.
    """
    record = next(iter(read_insertions(write_vcf(tmp_path, row(svlen=300)))))
    assert record.origin.length == len(record.sequence) == 60


def test_the_length_gate_still_reads_svlen(tmp_path):
    """The gate is about the variant the caller declared; the origin is about ALT."""
    path = write_vcf(tmp_path, row(svlen=300))
    assert list(read_insertions(path, max_length=200)) == []
    assert len(list(read_insertions(path, max_length=400))) == 1


def test_symbolic_alt_is_skipped(tmp_path):
    """<INS> declares an insertion without supplying sequence to scan."""
    path = write_vcf(tmp_path, row(alt="<INS>"))
    assert list(read_insertions(path)) == []


def test_non_insertions_are_skipped(tmp_path):
    path = write_vcf(tmp_path, row(svtype="DEL"))
    assert list(read_insertions(path)) == []


def test_missing_svlen_falls_back_to_the_length_of_the_alt(tmp_path):
    """SVLEN is optional in VCF 4.x, and the sequence is right there.

    Skipping the record lost the only thing this function returns while its
    length was in hand -- and a whole callset written without SVLEN came back
    empty, which reads as a file holding no insertions.
    """
    record = next(iter(read_insertions(write_vcf(tmp_path, row(svlen=None)))))
    assert len(record.sequence) == 60
    assert record.origin.length == 60


def test_the_fallback_length_is_still_gated(tmp_path):
    """It is a length like any other, not a way past ``min_length``."""
    short = "A" + "CAG" * 5      # 15 bp of insertion
    path = write_vcf(tmp_path, row(alt=short, svlen=None))
    assert list(read_insertions(path)) == []
    assert len(list(read_insertions(path, min_length=10))) == 1


def test_a_declared_svlen_still_wins_over_the_alt(tmp_path):
    """The fallback fires only when there is nothing to fall back from.

    A merged callset's SVLEN describes the consensus rather than the ALT it
    kept, and the gate is deliberately about the variant -- so a record that
    declares 300 is gated on 300 even though its ALT is 60 bp.
    """
    path = write_vcf(tmp_path, row(svlen=300))
    assert list(read_insertions(path, max_length=100)) == []
    assert len(list(read_insertions(path, max_length=500))) == 1


def test_length_bounds(tmp_path):
    short = "A" + "CAG" * 5      # 15 bp
    path = write_vcf(tmp_path, row(alt=short, svlen=15))
    assert list(read_insertions(path)) == []
    assert len(list(read_insertions(path, min_length=10))) == 1


def test_no_upper_bound_when_max_length_is_none(tmp_path):
    big = "A" + "CAG" * 5_000    # 15 kb, over the 10 kb default ceiling
    path = write_vcf(tmp_path, row(alt=big, svlen=15_000))
    assert list(read_insertions(path)) == []
    assert len(list(read_insertions(path, max_length=None))) == 1


def test_one_record_per_variant_with_the_variant_available(tmp_path):
    path = write_vcf(tmp_path, row(vid="INS0"), row(pos=500, vid="INS1"))
    pairs = list(read_insertions(path, with_variant=True))
    assert [r.name for r, _ in pairs] == ["INS0", "INS1"]
    assert [v.POS for _, v in pairs] == [100, 500]


def test_name_falls_back_to_position_when_id_is_missing(tmp_path):
    record = next(iter(read_insertions(write_vcf(tmp_path, row(vid=".")))))
    assert record.name == "chr1:100"
