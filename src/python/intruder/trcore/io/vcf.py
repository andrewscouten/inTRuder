"""Reading insertion sequence out of an SV VCF, through cyvcf2.

The sequence a long-read SV caller reports for an insertion is the ``ALT``
field, usually with the reference base at ``POS`` repeated as its first
character. Stripping that padding base is the one piece of real parsing here,
and it was previously inlined in ``pipeline/trf/sv_trfcaller.py``.

One record per insertion, not per sample. A VCF row carries a single ``ALT``,
so every sample genotyped for that row is scanned against identical sequence;
emitting one record per sample made the caller run ``pytrf`` over the same
string once per carrier. Sample-level facts -- genotype, ``DR``/``DV`` depth --
belong to the caller that needs them and stay on the ``cyvcf2`` variant, which
:func:`read_insertions` can hand back alongside each record.

Symbolic alternates (``<INS>``) are skipped. They declare an insertion without
supplying its sequence, so there is nothing to scan, and treating the literal
string ``<INS>`` as bases would find a repeat in punctuation.

``Inserted.length`` measures the sequence that was read, not ``SVLEN``. The two
are the same thing for a single-sample caller and routinely are not for a merged
one: in ``data/sv_output/survivor_multi_sample_vcf``, SURVIVOR's ``SVLEN``
disagrees with its own ``ALT`` on 297 of 481 insertions, by as much as 117 vs 53
bp, because the merged record keeps one representative ``ALT`` and a consensus
length. Every offset this package produces indexes ``ALT``, so ``SVLEN`` is the
wrong denominator for repeat coverage and the wrong bound to compare an offset
against -- with it, 1222 of the 29858 calls that file yields reported a
``rep_end`` past their own ``insert_size``. ``SVLEN`` is still what the length
gate below reads, because that gate is about the variant rather than about the
string -- falling back to the length of ``ALT`` only when the record declares
no ``SVLEN`` at all.
"""

from __future__ import annotations

from collections.abc import Iterator
from pathlib import Path

from ..coords import to_internal
from .records import Inserted, SeqRecord

#: VCF POS is 1-based; every offset in this package is 0-based.
_VCF_COORD_BASE = 1


def _insert_sequence(ref: str, alt: str) -> str:
    """``ALT`` minus the padding copy of ``REF`` at its head, when there is one."""
    if ref and alt[:len(ref)] == ref:
        return alt[len(ref):]
    return alt


def _svlen(variant) -> int | None:
    """``SVLEN`` as a positive int. Some callers write it as a list, some signed."""
    raw = variant.INFO.get("SVLEN")
    if raw is None:
        return None
    if isinstance(raw, (list, tuple)):
        if not raw:
            return None
        raw = raw[0]
    try:
        return abs(int(raw))
    except (TypeError, ValueError):
        return None


def sample_names(path: str | Path) -> list[str]:
    """The samples genotyped in ``path``, in column order.

    A caller fanning one call out across carriers needs these alongside the
    records, and asking for them here keeps cyvcf2 an implementation detail of
    this module: a step that imports it directly gets a bare ``ImportError``
    instead of :func:`_import_cyvcf2`'s message naming the extra to install.
    """
    VCF = _import_cyvcf2()
    return list(VCF(str(path)).samples)


def read_insertions(path: str | Path, *, min_length: int = 50,
                    max_length: int | None = 10_000,
                    with_variant: bool = False) -> Iterator:
    """Stream every ``SVTYPE=INS`` record of ``path`` as a :class:`SeqRecord`.

    Yields a record whose :class:`Inserted` origin carries the anchor position
    and ``SVLEN``. With ``with_variant`` the yield is ``(record, variant)``, so
    a caller that needs genotypes or read support has the ``cyvcf2`` variant
    without reopening the file.

    ``min_length``/``max_length`` gate on ``SVLEN`` and default to the bounds
    ``sv_trfcaller.py`` used; ``max_length=None`` removes the upper bound. The
    gate is the one place ``SVLEN`` is read: it is the caller's declared variant
    size, and :attr:`Inserted.length` is the length of ``ALT`` -- see below.

    A record that declares no ``SVLEN`` is gated on the length of the sequence
    it does carry, not skipped. ``SVLEN`` is optional in VCF 4.x and several
    callers leave it off a record whose ``ALT`` spells the insertion out in
    full; dropping those lost the one thing this function exists to return
    while the sequence -- and therefore its length -- was in hand. The fallback
    cannot change what a record with ``SVLEN`` does, so no callset in this
    repository moves: it only stops a callset without one coming back empty.
    """
    VCF = _import_cyvcf2()
    for variant in VCF(str(path)):
        if not variant.is_sv or variant.INFO.get("SVTYPE") != "INS":
            continue
        if not variant.ALT:
            continue
        alt = variant.ALT[0]
        if alt.startswith("<"):
            continue
        sequence = _insert_sequence(variant.REF, alt)
        if not sequence:
            continue

        length = _svlen(variant)
        if length is None:
            length = len(sequence)
        if length < min_length or (max_length is not None and length > max_length):
            continue

        record = SeqRecord(
            name=variant.ID or f"{variant.CHROM}:{variant.POS}",
            sequence=sequence.upper(),
            origin=Inserted(
                contig=variant.CHROM,
                position=to_internal(variant.POS, _VCF_COORD_BASE),
                length=len(sequence),
            ),
        )
        yield (record, variant) if with_variant else record


def _import_cyvcf2():
    """Import cyvcf2's ``VCF``, or explain which extra supplies it."""
    try:
        from cyvcf2 import VCF
    except ImportError as exc:                              # pragma: no cover - import guard
        raise ImportError(
            "reading VCF needs cyvcf2, which is not installed by default.\n"
            "  in this repo:  pixi run -e pipeline <your command>\n"
            "  to add it:     pixi add --feature sv cyvcf2\n"
            "  installed:     pip install 'intruder[sv]'"
        ) from exc
    return VCF
