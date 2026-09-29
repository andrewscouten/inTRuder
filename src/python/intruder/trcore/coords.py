"""Coordinate conventions, in one place because getting them wrong is silent.

Every step here works in **0-based half-open** coordinates internally and
converts once at the edges. The inputs do not agree: VCF ``POS`` is 1-based, a
BED interval is 0-based half-open, and UCSC, TRExplorer and STRchive are all
BED-style. Mixing the two shifts every interval by exactly one base, which
produces output that is plausible, self-consistent, and wrong.
"""

from __future__ import annotations


def normalize_chrom(chrom: str) -> str:
    """Map a contig name onto UCSC style (``1`` -> ``chr1``, ``MT`` -> ``chrM``).

    Every catalogue this project reads -- UCSC, TRExplorer, STRchive -- names
    contigs this way, so normalising on the way in means a query never misses
    only because its caller wrote ``1`` instead of ``chr1``.
    """
    name = str(chrom).strip()
    if not name.startswith("chr"):
        name = "chr" + name
    if name in ("chrMT", "chrmt"):
        name = "chrM"
    return name


def to_internal(pos, coord_base: int):
    """Input coordinate -> 0-based. A 1-based VCF POS is the base before the insert."""
    return pos - 1 if coord_base == 1 else pos


def to_external(start, end, coord_base: int):
    """0-based half-open interval -> the caller's convention."""
    return (start + 1, end) if coord_base == 1 else (start, end)


def union_length(intervals) -> int:
    """Bases covered by the union of 0-based half-open ``(start, end)`` intervals.

    A repeat finder reports overlapping calls over the same stretch of sequence
    -- one insertion in the sample data carries 64 of them -- so adding up each
    call's own length double-counts bases. As the numerator of "what fraction of
    this insertion is tandem repeat at all" that is not merely imprecise: it
    puts the fraction above 1, and a ``>= 0.8`` threshold on it then passes rows
    it was written to reject.

    Empty and reversed intervals contribute nothing. This is the scalar form of
    ``novelty.insertions.union_length``, which does the same sweep grouped over
    a whole DataFrame; the definition lives here so the two cannot drift.
    """
    ordered = sorted((int(start), int(end)) for start, end in intervals if end > start)
    total = 0
    open_start = open_end = None
    for start, end in ordered:
        if open_end is None:
            open_start, open_end = start, end
        elif start > open_end:                  # gap: bank the run and start a new one
            total += open_end - open_start
            open_start, open_end = start, end
        elif end > open_end:                    # overlapping or nested: extend
            open_end = end
    if open_end is not None:
        total += open_end - open_start
    return total


def interval_distance(start: int, end: int, other_start: int, other_end: int) -> int:
    """Distance in bp between two 0-based half-open intervals; ``0`` when they overlap.

    Directly adjacent intervals are 1 bp apart, so a query landing on the base
    immediately after a repeat is reported at distance 1 rather than 0 -- the
    caller decides with a window whether that is close enough.
    """
    if end > other_start and start < other_end:
        return 0
    if start >= other_end:                  # query lies to the right
        return start - other_end + 1
    return other_start - end + 1            # query lies to the left


def coverage_fraction(covered, size, *, digits: int = 3) -> float | None:
    """``covered / size`` as a fraction rounded to ``digits``, or ``None``.

    The numerator is a union of repeat calls and the denominator the length of
    the sequence they were found in -- what fraction of an insertion is tandem
    repeat at all. ``trf`` writes it as ``repeat_coverage`` and ``novelty`` as
    ``insertion_purity``, from the same two numbers, and the two must agree on
    every row: they are compared against the same threshold at different stages,
    so a row either step rounds the other way is a row the pipeline includes or
    excludes depending on which column was asked.

    Rounding is the part that drifted. ``round(308 / 320, 3)`` is ``0.963`` and
    ``numpy.round`` of the same ratio is ``0.962`` -- Python rounds the decimal
    value of a double that sits just above the tie, while numpy scales and
    rounds half to even. Both are defensible and neither is reproducible from
    the other, so the ratio is never formed as a float here. ``(2 * covered *
    scale + size) // (2 * size)`` is ``floor(covered / size * scale + 1/2)`` in
    exact integer arithmetic: ties round up, and the result is a property of the
    two integers rather than of a binary approximation to their quotient. A
    vectorised implementation over int64 reproduces it exactly, which is what
    lets the agreement tests pin one to the other.

    Capped at 1: a call can run past the reported length of the sequence it was
    found in. ``None`` when there is no denominator to divide by -- that is the
    absence of a measurement rather than a coverage of zero, and the two must
    stay distinguishable, since zero fails every threshold and a missing value
    is judged by none.
    """
    if covered is None or size is None:
        return None
    covered, size = int(covered), int(size)
    if size <= 0 or covered < 0:
        return None
    scale = 10 ** digits
    scaled = (2 * covered * scale + size) // (2 * size)
    return min(scaled / scale, 1.0)
