"""Reading the table columns that do not hold a plain number.

Stage 01's table is passed between steps as a TSV, and two of its columns carry
something other than a bare integer:

``insert_size``
    ``SVLEN``, which a caller may write signed, and which ``sv_trfcaller.py``
    wrote straight through from a VCF list field -- so the column reads
    ``[415]`` as often as ``415``.

``depth``
    ``"DV,DR"``: variant and reference read support in one cell, which is a
    total only once the pair is summed.

Every step that filters or features a stage 01 table has to read those, and
each one that re-derived the parse drifted from the others: the scalar reader
in ``pipeline.trf.filters`` returned ``-415`` where the vectorised
``novelty.insertions.parse_sizes`` returned ``415``, and returned ``None``
where that one recovered ``415`` from ``[415,200]``. A size one step can read
and another cannot is not a rounding difference -- it is the denominator of
repeat coverage, so the two steps disagree about which rows pass a
``--min-coverage`` threshold while both look correct in isolation.

So the scalar definition lives here, the same arrangement
:func:`trcore.coords.union_length` already has: pure standard library at this
level, and a vectorised wrapper beside the step that needs catalogue-scale
throughput. The wrappers are held to these semantics by the agreement tests in
``tests/python/trcore/utils/test_parse.py``.

A value none of these can read comes back as ``None``, never as ``0``. Zero is
a measurement and fails every threshold; ``None`` is the absence of one, which
is what a filter reads as "nothing to judge" and passes. Collapsing the two
drops exactly the rows the pipeline knows least about.
"""

from __future__ import annotations

import re

#: The first signed integer in a field, for a value that will not parse whole.
_FIRST_INTEGER = re.compile(r"-?\d+")


def as_number(raw) -> float | int | None:
    """``raw`` as a number, or ``None`` when it is not one.

    Integral text stays an ``int`` rather than becoming a float, so a motif
    length read back out of a TSV still compares equal to the ``int`` in a
    ``--keep-motif-sizes`` set.
    """
    if isinstance(raw, (int, float)):
        return raw
    try:
        text = str(raw).strip()
        return int(text) if text.lstrip("-").isdigit() else float(text)
    except (TypeError, ValueError):
        return None


def as_size(raw) -> int | None:
    """An insertion length, however the stage that wrote it spelled ``SVLEN``.

    Tolerates the ``[415]`` a VCF list field writes, recovers the first integer
    from anything else that will not parse whole (``[415,200]``, ``415bp``),
    and takes the absolute value -- a deletion-signed ``SVLEN`` is a length, and
    :func:`trcore.io.vcf._svlen` already normalises it the same way on the way
    in. Without the ``abs`` a negative size makes repeat coverage negative, so
    the row fails every ``--min-coverage`` rather than being judged on what it
    actually covers.

    This is the scalar form of ``novelty.insertions.parse_sizes``; the two are
    pinned to each other by test_parse.py, because a size only one of them can
    read means the filter step and the screen step disagree about the same row.
    """
    value = as_number(raw)
    if value is None:
        match = _FIRST_INTEGER.search(str(raw))
        value = int(match.group()) if match else None
    if value is None:
        return None
    try:
        return abs(int(value))
    except (OverflowError, ValueError):     # nan and inf are not lengths
        return None


def as_total(raw) -> float | int | None:
    """``"15,20"`` -> ``35``: total read support from a ``DV,DR`` pair.

    A bare number is already a total and passes through, which is what the
    column holds once an earlier filter pass has collapsed it.

    Half a pair is not two thirds of an answer: if any part will not parse the
    whole field is unreadable, and the caller is told nothing rather than a
    number built from the parts that happened to survive.
    """
    if isinstance(raw, (int, float)):
        return raw
    parts = [part for part in str(raw).strip().split(",") if part.strip()]
    values = [as_number(part) for part in parts]
    if not values or any(value is None for value in values):
        return None
    return sum(values)
