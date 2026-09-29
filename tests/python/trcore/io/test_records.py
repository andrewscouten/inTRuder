"""The coordinate origins, which are the whole reason this record type exists."""

from __future__ import annotations

from intruder.trcore.io import Genomic, Inserted, SeqRecord


def test_genomic_offsets_are_reference_positions():
    origin = Genomic(contig="NZ_CP012345.1")
    assert origin.locate(10, 25) == ("NZ_CP012345.1", 10, 25)
    assert origin.in_reference is True


def test_genomic_window_reports_position_on_the_full_contig():
    """A repeat at offset 40 of a window starting at 1000 sits at 1040, not 40."""
    origin = Genomic(contig="chr1", start=1_000)
    assert origin.locate(40, 55) == ("chr1", 1_040, 1_055)


def test_inserted_has_no_reference_interval():
    """The reference has nothing there, so there is no interval to hand back."""
    origin = Inserted(contig="chr1", position=99, length=138)
    assert origin.locate(10, 25) is None
    assert origin.in_reference is False


def test_inserted_anchor_is_the_base_the_insertion_follows():
    origin = Inserted(contig="chr1", position=99, length=138)
    assert origin.anchor() == ("chr1", 99)


def test_record_length_is_its_sequence_length():
    record = SeqRecord(name="c1", sequence="ACGTACGT", origin=Genomic(contig="c1"))
    assert len(record) == 8
    assert record.description == ""


def test_origins_are_frozen():
    """Hashable, so a scan can key results by origin without copying defensively."""
    assert Genomic("chr1", 0) == Genomic("chr1", 0)
    assert len({Genomic("chr1"), Genomic("chr1"), Genomic("chr2")}) == 2
