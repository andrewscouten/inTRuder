"""Running the finder over records, and what a call's offsets mean."""

from __future__ import annotations

import pytest

from intruder.pipeline.trf.find import RepeatCall, find_in_record, find_repeats
from intruder.pipeline.trf.params import GENOME, FinderParams
from intruder.trcore.io import Genomic, Inserted, SeqRecord
from intruder.trcore.motifs import primitive_unit

pytest.importorskip("pytrf")

FLANK = "GTCCATGAGCTTACGGATCCTAGGCATTCAGGACTTAGCCA"      # non-repetitive


def genomic(sequence, contig="contig_1", start=0):
    return SeqRecord(name=contig, sequence=sequence,
                     origin=Genomic(contig=contig, start=start))


def test_a_planted_repeat_is_found_with_its_motif():
    record = genomic(FLANK + "CAG" * 30 + FLANK)
    calls = find_in_record(record, GENOME)
    assert calls, "a 90 bp perfect CAG tract should be found"
    assert any(call.motif in ("CAG", "AGC", "GCA") for call in calls)


def test_offsets_are_zero_based_half_open_into_the_record():
    """pytrf reports 1-based inclusive; the conversion happens once, in find.

    Pinned with ``max_extension=0`` so the call stops at the tract's real edge:
    with any extension budget the finder deliberately runs past it through
    mismatches, and a boundary that has moved cannot show an off-by-one.
    """
    record = genomic(FLANK + "CAG" * 30 + FLANK)
    call = max(find_in_record(record, GENOME.replace(max_extension=0)),
               key=lambda c: c.rep_length)
    assert call.start == len(FLANK)
    assert call.end == len(FLANK) + 90
    assert record.sequence[call.start:call.end] == "CAG" * 30
    assert call.rep_length == call.end - call.start


def test_extension_runs_past_the_tract_into_the_flank():
    """The trade `max_extension` makes, asserted rather than assumed."""
    record = genomic(FLANK + "CAG" * 30 + FLANK)
    exact = max(find_in_record(record, GENOME.replace(max_extension=0)),
                key=lambda c: c.rep_length)
    extended = max(find_in_record(record, GENOME.replace(max_extension=20)),
                   key=lambda c: c.rep_length)
    assert extended.start < exact.start and extended.end > exact.end
    assert extended.purity < exact.purity


def test_purity_is_a_fraction_not_a_percentage():
    record = genomic(FLANK + "CAG" * 30 + FLANK)
    for call in find_in_record(record, GENOME):
        assert 0.0 <= call.purity <= 1.0


def test_identity_floor_is_honoured():
    """With the scale conversion in place, the floor actually excludes calls."""
    record = genomic(FLANK + "CAG" * 30 + FLANK)
    assert find_in_record(record, GENOME.replace(min_identity=0.5))
    assert not find_in_record(record, GENOME.replace(min_identity=1.01))


def test_a_genomic_call_locates_to_a_reference_interval():
    record = genomic(FLANK + "CAG" * 30 + FLANK, contig="chrX")
    call = max(find_in_record(record, GENOME), key=lambda c: c.rep_length)
    contig, start, end = call.locate()
    assert (contig, start, end) == ("chrX", call.start, call.end)


def test_a_windowed_record_reports_positions_on_the_full_contig():
    """A window starting at 1000 must not report its repeat at offset 40."""
    record = genomic(FLANK + "CAG" * 30 + FLANK, contig="chrX", start=1_000)
    call = max(find_in_record(record, GENOME), key=lambda c: c.rep_length)
    _, start, end = call.locate()
    assert start == call.start + 1_000
    assert end == call.end + 1_000


def test_an_inserted_call_has_no_reference_interval():
    """The whole point of an insertion: the reference has nothing there."""
    record = SeqRecord(name="INS0", sequence=FLANK + "CAG" * 30 + FLANK,
                       origin=Inserted(contig="chr1", position=99, length=90))
    call = max(find_in_record(record, GENOME), key=lambda c: c.rep_length)
    assert call.locate() is None


def test_calls_do_not_carry_the_record_sequence():
    """A contig scan emits many calls; pinning the contig to each would not scale."""
    record = genomic(FLANK + "CAG" * 30 + FLANK)
    call = find_in_record(record, GENOME)[0]
    assert not any(isinstance(value, str) and len(value) > 100
                   for value in vars(call).values())


def test_find_repeats_streams_across_records():
    records = [genomic(FLANK + "CAG" * 30 + FLANK, contig="a"),
               genomic(FLANK + "AT" * 45 + FLANK, contig="b")]
    contigs = {call.origin.contig for call in find_repeats(records, GENOME)}
    assert contigs == {"a", "b"}


# --------------------------------------------------------------------------- #
# the motif, reduced to the unit the rest of the project keys on
# --------------------------------------------------------------------------- #

def test_a_homopolymer_is_reported_as_its_primitive_unit():
    """ATRFinder reports the unit it seeded on, which can be a power of the real one.

    A poly-G came back as ``GGGG``: `motif_length` 4 rather than 1, `rep_units`
    a quarter of the copies, and `--drop-motif-sizes 1` no longer dropping a
    homopolymer.
    """
    record = genomic(FLANK + "G" * 40 + FLANK)
    call = max(find_in_record(record, GENOME), key=lambda c: c.rep_length)
    assert call.motif == "G"
    assert call.motif_length == 1


def test_a_dinucleotide_tract_is_not_reported_as_its_square():
    record = genomic(FLANK + "GT" * 30 + FLANK)
    call = max(find_in_record(record, GENOME), key=lambda c: c.rep_length)
    assert call.motif in ("GT", "TG")


def test_the_reported_motif_is_never_a_repetition_of_a_shorter_unit():
    """The property, over every call of a scan rather than a chosen one."""
    record = genomic(FLANK + "G" * 40 + FLANK + "GT" * 30 + FLANK + "CAG" * 30 + FLANK)
    for call in find_in_record(record, GENOME):
        assert primitive_unit(call.motif) == call.motif


def test_rep_units_counts_copies_of_the_reduced_motif():
    """Reducing the motif without re-deriving the count would halve the copies."""
    record = genomic(FLANK + "G" * 40 + FLANK)
    call = max(find_in_record(record, GENOME), key=lambda c: c.rep_length)
    assert call.rep_units == call.rep_length // len(call.motif)
    assert call.rep_units >= 40


def test_a_genuine_motif_is_left_alone():
    """`primitive_unit` fires only on an exact power, so CAG stays CAG."""
    record = genomic(FLANK + "CAG" * 30 + FLANK)
    call = max(find_in_record(record, GENOME), key=lambda c: c.rep_length)
    assert len(call.motif) == 3


def test_rep_units_floors_whole_copies():
    call = RepeatCall(name="x", origin=Genomic("c"), start=0, end=10,
                      motif="CAG", purity=1.0)
    assert call.rep_units == 3          # 10 bp spans 3 whole copies of a 3-mer


def test_empty_sequence_yields_nothing():
    assert find_in_record(genomic(""), FinderParams()) == []
