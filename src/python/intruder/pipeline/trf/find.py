"""Running the repeat finder over records, whatever produced them.

This is the whole of the calling step, and it is deliberately small. It takes
:class:`trcore.io.SeqRecord` objects and yields :class:`RepeatCall` objects; it
opens no files, knows no formats, and never asks where a record came from. A
VCF insertion and a bacterial contig reach it as the same thing, which is the
point -- ``pytrf.ATRFinder`` only ever wanted a name and a string.

What the finder returns is in *its* conventions, and they are converted here,
once: ``repeat.start`` is 1-based inclusive, and ``repeat.identity`` is a
percentage. Everything downstream of this module sees 0-based half-open offsets
and a purity in ``0..1``. The coordinate half of that goes through
:func:`trcore.coords.to_internal`, the same call ``trcore.io.vcf`` makes on a
VCF ``POS`` -- a second 1-based source converted by a second hand-written ``-
1`` is exactly the drift that module exists to prevent.

The motif is reduced the same way, through
:func:`trcore.motifs.primitive_unit`. ``ATRFinder`` reports the unit it happened
to seed on, which is sometimes a whole-number power of the real one: a 32 bp
poly-G comes back as ``GGGG``, and a ``GT`` dinucleotide tract as ``GTGT``. Rare
-- 25 of the 18166 calls in a scan of this repository's insertion sequences --
and wrong in three places at once when it happens. ``motif_length`` is a
multiple of the period, so ``--drop-motif-sizes 1`` does not drop a homopolymer
reported as ``GGGG``; ``rep_units`` is divided down by that same multiple; and
``novelty`` and ``strchive`` both key the call through
:func:`trcore.motifs.canonical_motif`, which reduces to the primitive unit --
so the table said ``GGGG`` about the very call they had already agreed was
``G``. Reducing here is safe in a way a general rewrite would not be: it fires
only when the reported motif is exactly ``unit * k``, and such a motif is a
repeat of ``unit`` by definition.

A call's offsets are measured inside the record it was found in. Turning them
into a reference interval is :meth:`trcore.io.Origin.locate`'s job and only
succeeds for a genomic record -- see :mod:`trcore.io.records` for why an
insertion has no interval to give.
"""

from __future__ import annotations

import os
import sys
from collections.abc import Iterable, Iterator
from contextlib import contextmanager
from dataclasses import dataclass

from intruder.trcore.coords import to_internal
from intruder.trcore.io import Origin, SeqRecord
from intruder.trcore.motifs import primitive_unit

from .params import FinderParams

#: pytrf reports ``repeat.start`` 1-based inclusive, as a VCF reports POS.
_PYTRF_COORD_BASE = 1


@dataclass(frozen=True)
class RepeatCall:
    """One tandem repeat, with offsets measured inside the record holding it.

    The record's sequence is deliberately not carried. A genome scan holds one
    contig at a time but can emit tens of thousands of calls from it, and
    attaching the string to each would pin megabytes per call; the origin is
    enough to say where the call sits, and anything needing the bases again
    (flank features, for instance) re-reads them from the source.
    """

    name: str
    origin: Origin
    start: int          # 0-based half-open, within the record
    end: int
    motif: str
    purity: float       # 0..1

    @property
    def motif_length(self) -> int:
        return len(self.motif)

    @property
    def rep_length(self) -> int:
        return self.end - self.start

    @property
    def rep_units(self) -> int:
        """Whole copies of the motif spanned, floored, as stage 01 has always reported it."""
        return self.rep_length // self.motif_length if self.motif else 0

    def locate(self) -> tuple[str, int, int] | None:
        """The reference interval this call sits on, or ``None`` inside an insertion."""
        return self.origin.locate(self.start, self.end)


def find_in_record(record: SeqRecord, params: FinderParams) -> list[RepeatCall]:
    """Every repeat in one record, as a list.

    Returned eagerly rather than streamed because a caller that needs the
    record's own context -- the VCF writer needs the variant's genotypes to fan
    a call out across samples -- has to hold the calls and the source together
    anyway, and scanning once per record rather than once per sample is exactly
    the saving that motivates one record per insertion.
    """
    finder = _import_pytrf()
    with _quiet(params.quiet):
        found = list(finder.ATRFinder(record.name, record.sequence, **params.as_kwargs()))
    return [
        RepeatCall(
            name=record.name,
            origin=record.origin,
            start=to_internal(repeat.start, _PYTRF_COORD_BASE),
            end=repeat.end,                     # 1-based inclusive end == 0-based exclusive
            motif=primitive_unit(repeat.motif),
            purity=float(repeat.identity) / 100,
        )
        for repeat in found
    ]


def find_repeats(records: Iterable[SeqRecord],
                 params: FinderParams) -> Iterator[RepeatCall]:
    """Every repeat in every record, streamed."""
    for record in records:
        yield from find_in_record(record, params)


@contextmanager
def _quiet(enabled: bool):
    """Redirect fds 1 and 2 to devnull while the finder runs.

    pytrf writes directly to the C-level descriptors, so a Python-level
    redirect does not reach it. The cost is that *everything* written while it
    is in force disappears, including a genuine failure -- tolerable around a
    short insertion, much less so around a multi-megabyte contig scan, which is
    why ``--finder-errors`` turns it off and why the suppression is scoped to
    one record rather than wrapped around the whole run.
    """
    if not enabled:
        yield
        return

    sys.stdout.flush()
    sys.stderr.flush()
    saved = (os.dup(1), os.dup(2))
    devnull = os.open(os.devnull, os.O_WRONLY)
    try:
        os.dup2(devnull, 1)
        os.dup2(devnull, 2)
        os.close(devnull)
        yield
    finally:
        sys.stdout.flush()
        sys.stderr.flush()
        os.dup2(saved[0], 1)
        os.dup2(saved[1], 2)
        os.close(saved[0])
        os.close(saved[1])


def _import_pytrf():
    """Import pytrf, or explain which extra supplies it."""
    try:
        import pytrf
    except ImportError as exc:                              # pragma: no cover - import guard
        raise ImportError(
            "finding repeats needs pytrf, which is not installed by default.\n"
            "  in this repo:  pixi run -e pipeline <your command>\n"
            "  to add it:     pixi add --feature trf pytrf\n"
            "  installed:     pip install 'intruder[trf]'"
        ) from exc
    return pytrf
