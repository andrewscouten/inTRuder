"""Turning calls into rows, in the schema their origin allows.

The finder's output is the same either way; what can be *written down* about a
call is not, and the split is the one :mod:`trcore.io.records` draws. A call
inside an insertion has no reference interval -- ``Origin.locate`` returns
``None`` -- but it does have the variant it came from, so its row carries the
SV's identity and the sample's read support. A call on a contig has the
opposite: a real interval to write to BED, and no sample at all.

Trying to serve both from one column set is what produces a table with a
``start`` column that indexes nothing for half its rows, so there are two.
:data:`INSERTION_COLUMNS` is unchanged from ``sv_trfcaller.py``, because
``novelty``, ``strchive`` and the isolation forest all read tables in it.
"""

from __future__ import annotations

import csv
from collections.abc import Iterable, Iterator, Sequence
from pathlib import Path

from .find import RepeatCall

#: Stage 01's long-standing schema for calls inside SV insertions.
INSERTION_COLUMNS = (
    "chrom", "ins_coord", "SVID", "depth", "insert_size", "sample", "allele",
    "rep_start", "rep_end", "motif", "purity", "motif_length", "rep_length",
    "rep_units",
)

#: Calls on reference sequence. No sample, depth or insertion columns.
GENOME_COLUMNS = (
    "genome", "contig", "start", "end", "motif", "purity", "motif_length",
    "rep_length", "rep_units",
)

#: BED4, the schema ``novelty.platforms`` reads a catalogue in.
BED_COLUMNS = ("contig", "start", "end", "motif")

#: The ALT index the scanned sequence came from. ``trcore.io.read_insertions``
#: lifts ``ALT[0]`` and nothing else, so allele 1 is the only one a call found
#: in that sequence is evidence about -- see :func:`insertion_rows`.
SCANNED_ALLELE = 1


def insertion_rows(call: RepeatCall, variant, samples: Sequence[str]) -> Iterator[dict]:
    """One row per carried allele per sample, for a call inside ``variant``.

    Fans a single call out across the samples genotyped for the variant. The
    scan happens once per insertion upstream of this, not once per carrier --
    every sample on a VCF row is genotyped against the identical ALT sequence,
    so scanning per sample only ever re-derived the same calls.

    Homozygous-alt genotypes yield two identical rows, one per allele copy.
    That is stage 01's existing behaviour and the duplicate-row issue the
    reporting scripts compensate for by counting distinct ``(sample, SVID)``
    pairs; it is preserved here so tables stay comparable across the rewrite.

    Only allele :data:`SCANNED_ALLELE` is written. A multi-allelic row carries
    one ALT per alternate and the reader scanned the first, so a ``1/2``
    genotype used to emit a second row claiming this call sits in an allele
    whose sequence was never looked at -- a fabricated observation, and one that
    counts again in every distinct-``(sample, SVID)`` tally downstream. Sniffles
    and SURVIVOR both write biallelic INS records, so no table in this
    repository changes; what changes is what happens when one does not.

    A VCF with no sample columns yields one *unattributed* row per call rather
    than none. Fanning out across samples is how a call reaches the table, so a
    sites-only callset used to scan normally, find its repeats and then drop
    every one of them -- an empty table that reads exactly like a scan which
    found nothing. ``sample``, ``allele`` and ``depth`` are left empty, which is
    the "no value to judge" every reader of this table already understands: it
    is what :meth:`filters.Check.fails` passes rather than fails, and what the
    distinct-``(sample, SVID)`` tallies count once. A zero or a placeholder name
    would instead be a sample fact the VCF never stated.
    """
    shared = {
        "chrom": variant.CHROM,
        "ins_coord": variant.POS,
        "SVID": variant.ID,
        "insert_size": getattr(call.origin, "length", None),
        "rep_start": call.start,
        "rep_end": call.end,
        "motif": call.motif,
        "purity": round(call.purity, 3),
        "motif_length": call.motif_length,
        "rep_length": call.rep_length,
        "rep_units": call.rep_units,
    }

    if not samples:
        yield {**shared, "depth": "", "sample": "", "allele": ""}
        return

    depth_variant = _format_field(variant, "DV")
    depth_reference = _format_field(variant, "DR")
    genotypes = _genotypes(variant)

    for index, sample in enumerate(samples):
        for allele in _alleles(genotypes[index]):
            if allele != SCANNED_ALLELE:
                continue
            yield {
                **shared,
                "depth": _depth(depth_variant, depth_reference, index),
                "sample": sample,
                "allele": allele,
            }


def _genotypes(variant):
    """``variant.genotypes``, or one line saying which FORMAT field is absent.

    A VCF may carry sample columns and no ``GT`` -- a caller reporting read
    support per sample and leaving the genotype to a later step writes one --
    and cyvcf2 answers that with a bare ``Exception("error parsing
    genotypes")``. Bare, so neither of :func:`cli.main`'s handlers catches it:
    an ordinary input came out as a traceback through this module naming
    neither the file nor the field. ``RuntimeError`` is what the parsers raise
    about an input they cannot read, which is what makes this one of them.

    The ``except`` is as wide as the thing it translates: cyvcf2 raises the
    base class here, and nothing else happens inside the ``try``.
    """
    try:
        return variant.genotypes
    except Exception as exc:                # cyvcf2 raises the base class here
        raise RuntimeError(
            f"{variant.CHROM}:{variant.POS} has sample columns but no readable GT; "
            f"the insertion schema attributes a call to the samples carrying it, "
            f"which needs genotypes ({exc})"
        ) from exc


def _alleles(genotype: Sequence) -> Sequence:
    """The allele indices of a cyvcf2 genotype, without its trailing phase flag.

    ``variant.genotypes`` yields one allele per copy followed by a bool saying
    whether the call is phased -- ``[0, 1, False]`` for ``0/1``, and ``[1, True]``
    for the bare ``1`` a haploid call writes. Taking ``genotype[:2]`` reads that
    bool as a second allele, and ``True == 1`` in Python, so a haploid carrier of
    :data:`SCANNED_ALLELE` emitted two rows instead of one -- a second copy of an
    allele the sample does not have, with ``True`` written into the numeric
    ``allele`` column. Every chrX and chrY call in a male sample is haploid, and
    the duplicate counts again in the distinct-``(sample, SVID)`` tallies the
    reporting scripts run.

    Dropping the last element rather than keeping the first two is what makes
    that independent of ploidy, which is the thing actually varying here.
    """
    return genotype[:-1]


def _format_field(variant, tag: str):
    """A FORMAT column, or ``None`` when this VCF does not declare it.

    ``DR``/``DV`` are Sniffles' read-support tags, not part of the VCF spec, and
    a caller that does not emit them is an ordinary input rather than a broken
    one. cyvcf2 raises a bare ``KeyError: b'DV'`` from inside its header lookup
    for a tag it cannot find, which surfaces as a traceback out of a whole-VCF
    scan and says nothing about which file or which field.
    """
    try:
        return variant.format(tag)
    except KeyError:
        return None


def _depth(variant_support, reference_support, index: int) -> str:
    """``"DV,DR"`` for one sample, or empty unless both halves were measured.

    Both, because :func:`trcore.utils.parse.as_total` reads the cell as a sum and
    declines half a pair -- writing one measured half beside a missing one would
    produce a cell that parses to nothing anyway, spelled as though it held data.

    Empty rather than ``0,0``: zero is a measurement and this is the absence of
    one. :func:`trcore.utils.parse.as_total` reads the pair back as a sum, and an
    empty cell is what :meth:`filters.Check.fails` passes as "nothing to judge"
    -- so a depth threshold against a VCF without these tags filters nothing,
    instead of silently rejecting every call for having no support recorded.

    A sample whose ``DR:DV`` is ``.`` reaches here the same way, one level down:
    the tags are declared, so cyvcf2 returns an array, and the cell it holds is
    htslib's missing-integer sentinel. Written through, the column read
    ``-2147483648,-2147483648``, which sums to a large negative depth -- not an
    unreadable cell a threshold declines to judge, but a measurement that fails
    every ``--min-depth`` and trains the isolation forest on the sentinel. Read
    support cannot be negative, so any negative value is that sentinel.
    """
    if variant_support is None or reference_support is None:
        return ""
    variant_reads, reference_reads = variant_support[index][0], reference_support[index][0]
    if variant_reads < 0 or reference_reads < 0:
        return ""
    return f"{variant_reads},{reference_reads}"


def genome_row(call: RepeatCall, genome: str) -> dict:
    """One row for a call on reference sequence.

    ``genome`` names the assembly the call came from and is supplied by the
    caller, from the file being scanned -- not read out of the FASTA header,
    which every reference source lays out differently (see
    :class:`trcore.io.SeqRecord`).
    """
    located = call.locate()
    if located is None:                     # pragma: no cover - guarded by the CLI
        raise ValueError(
            f"{call.name!r} has no reference interval: a genome row needs a Genomic "
            "origin, and this call came from an insertion"
        )
    contig, start, end = located
    return {
        "genome": genome,
        "contig": contig,
        "start": start,
        "end": end,
        "motif": call.motif,
        "purity": round(call.purity, 3),
        "motif_length": call.motif_length,
        "rep_length": call.rep_length,
        "rep_units": call.rep_units,
    }


def write_tsv(path: str | Path, columns: Sequence[str], rows: Iterable[dict],
              *, extras: Sequence[str] = ()) -> int:
    """Write ``rows`` as a TSV with a header; returns how many were written.

    ``extras`` are columns appended after ``columns`` -- the filter step's
    ``filter_status``, for instance -- so a schema can be extended downstream
    without redefining it here. An extra already present in ``columns`` is not
    appended twice: filtering an already-tagged table is an ordinary thing to
    do, and ``csv.DictWriter`` would otherwise emit the header twice over.
    """
    fields = list(dict.fromkeys([*columns, *extras]))
    written = 0
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t",
                                extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
            written += 1
    return written


def write_bed(path: str | Path, rows: Iterable[dict]) -> int:
    """Write genome rows as headerless BED4, for the read filter and the CHM13 screen."""
    written = 0
    with open(path, "w", newline="") as handle:
        for row in rows:
            handle.write("\t".join(str(row[column]) for column in BED_COLUMNS) + "\n")
            written += 1
    return written
