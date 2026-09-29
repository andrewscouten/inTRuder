"""A named sequence to scan, and what its offsets are measured against.

Every repeat finder in this repository hands a string to ``pytrf`` and gets back
offsets into that string. What those offsets *mean* depends entirely on where
the string came from, and the two sources disagree in a way that does not
announce itself:

``Genomic``
    The string is a contig of the reference. An offset is a position on that
    contig, and a repeat call is a reference interval you can write to BED,
    intersect with another catalogue, or screen against CHM13.

``Inserted``
    The string is sequence a sample carries and the reference does not, lifted
    out of a VCF ``ALT``. An offset is a position *inside the insertion*. There
    is no reference interval to write down: the whole point of the record is
    that the reference has nothing there. All the reference can supply is the
    anchor base the insertion follows.

Conflating them produces a BED file of coordinates that index nothing, which
reads as valid and lands wherever the insertion happened to sit. So the origin
travels with the sequence, and :meth:`locate` returns ``None`` rather than a
plausible number when no reference interval exists.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class Genomic:
    """Offsets into reference sequence, at ``start`` along ``contig``.

    ``start`` is non-zero only when the record holds part of a contig -- a
    window fetched by region rather than a whole scan -- so that a repeat found
    in the window still reports where it sits on the full contig.
    """

    contig: str
    start: int = 0

    #: A call here is a reference interval.
    in_reference = True

    def locate(self, start: int, end: int) -> tuple[str, int, int]:
        """Offsets within the record -> a 0-based half-open reference interval."""
        return (self.contig, self.start + start, self.start + end)


@dataclass(frozen=True)
class Inserted:
    """Offsets into sequence inserted after ``position`` on ``contig``.

    ``position`` is 0-based: the reference base the insertion follows, which is
    :func:`trcore.coords.to_internal` applied to a VCF ``POS``. ``length`` is
    the length of the inserted sequence itself, kept here because it belongs to
    the insertion rather than to any one repeat call inside it -- it is the
    denominator ``novelty.insertions`` divides repeat coverage by, and the bound
    every offset measured in this origin falls inside.

    It is deliberately not the ``SVLEN`` the VCF declared. A merged callset's
    ``SVLEN`` is a consensus over the calls it merged and need not describe the
    one ``ALT`` it kept, and dividing a union of offsets into ``ALT`` by a
    number describing something else is not a coverage fraction at all. See
    :mod:`trcore.io.vcf` for the measurement.
    """

    contig: str
    position: int
    length: int | None = None

    #: A call here sits in sequence the reference does not have.
    in_reference = False

    def locate(self, start: int, end: int) -> None:
        """Always ``None``. See the module docstring for why this is not an oversight.

        The signature matches :meth:`Genomic.locate` so a caller can hold either
        origin, but there is no reference interval to return: the reference has
        no sequence here at all.
        """

    def anchor(self) -> tuple[str, int]:
        """The reference base the insertion follows."""
        return (self.contig, self.position)


#: Either kind of origin. Both carry ``contig``, ``in_reference`` and ``locate``.
Origin = Genomic | Inserted


@dataclass(frozen=True)
class SeqRecord:
    """One named sequence, and where its offsets are measured from.

    ``description`` is whatever the source had left over after the name: the
    rest of a FASTA header, empty for a VCF insertion. It is carried rather
    than parsed because no two references agree on what goes in there --
    NCBI, Ensembl, prokka and the HMP reference genomes all write it
    differently, and the organism string this project needs is in it.
    """

    name: str
    sequence: str
    origin: Origin
    description: str = ""

    def __len__(self) -> int:
        return len(self.sequence)
