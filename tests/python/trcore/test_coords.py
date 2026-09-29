"""Coordinate conventions -- the shared definition of 'where' and 'how far'.

Both steps convert VCF coordinates and measure distance to an interval. When
those two disagree by a base, the tables they produce look comparable and are
not, so the rules live in one place and are tested once.
"""

from __future__ import annotations

import pytest

from intruder.trcore.coords import (
    coverage_fraction,
    interval_distance,
    normalize_chrom,
    to_external,
    to_internal,
    union_length,
)


@pytest.mark.parametrize("raw,expected", [
    ("chr1", "chr1"), ("1", "chr1"), (" 1 ", "chr1"), ("X", "chrX"),
    ("chrX", "chrX"), ("MT", "chrM"), ("chrMT", "chrM"), ("M", "chrM"),
    (1, "chr1"),
])
def test_normalize_chrom(raw, expected):
    assert normalize_chrom(raw) == expected


def test_to_internal_shifts_only_one_based_input():
    assert to_internal(10001, coord_base=1) == 10000
    assert to_internal(10000, coord_base=0) == 10000


def test_to_external_is_the_inverse():
    assert to_external(10000, 10468, coord_base=1) == (10001, 10468)
    assert to_external(10000, 10468, coord_base=0) == (10000, 10468)
    start = to_internal(10001, coord_base=1)
    assert to_external(start, start + 1, coord_base=1)[0] == 10001


# The interval under test is [100, 200): first base 100, last base 199.
@pytest.mark.parametrize("start,end,expected", [
    (100, 101, 0),      # first base
    (199, 200, 0),      # last base
    (100, 200, 0),      # exactly the interval
    (50, 150, 0),       # straddles the left edge
    (150, 250, 0),      # straddles the right edge
    (0, 1000, 0),       # contains it
    (99, 100, 1),       # one base to the left
    (200, 201, 1),      # one base to the right (end is exclusive)
    (0, 100, 1),        # abuts the left edge
    (90, 91, 10),
    (209, 210, 10),
])
def test_interval_distance_is_symmetric_and_zero_on_overlap(start, end, expected):
    assert interval_distance(start, end, 100, 200) == expected


def test_interval_distance_is_commutative():
    for a, b in [((0, 10), (100, 200)), ((150, 160), (100, 200)), ((300, 400), (100, 200))]:
        assert interval_distance(*a, *b) == interval_distance(*b, *a)


def test_union_length_counts_overlapping_intervals_once():
    """The correction this exists for: a repeat finder emits overlapping calls."""
    assert union_length([(0, 40), (10, 40), (20, 40)]) == 40


def test_union_length_adds_disjoint_intervals():
    assert union_length([(0, 20), (60, 90)]) == 50


def test_union_length_merges_adjacent_intervals():
    assert union_length([(0, 20), (20, 40)]) == 40


def test_union_length_handles_nesting():
    assert union_length([(0, 100), (10, 20), (30, 40)]) == 100


def test_union_length_ignores_empty_and_reversed_intervals():
    assert union_length([]) == 0
    assert union_length([(10, 10), (30, 20)]) == 0


def test_union_length_does_not_depend_on_input_order():
    intervals = [(60, 90), (0, 20), (10, 40)]
    assert union_length(intervals) == union_length(sorted(intervals, reverse=True))


# --------------------------------------------------------------------------- #
# coverage_fraction
# --------------------------------------------------------------------------- #

def test_coverage_fraction_is_the_ratio_rounded_to_three_places():
    assert coverage_fraction(80, 100) == 0.8
    assert coverage_fraction(1, 3) == 0.333


def test_coverage_fraction_is_capped_at_one():
    """A call can run past the reported length of the sequence holding it."""
    assert coverage_fraction(150, 100) == 1.0


def test_zero_coverage_is_still_a_measurement():
    assert coverage_fraction(0, 100) == 0.0


def test_no_denominator_is_not_a_coverage_of_zero():
    """``None``, which a filter reads as nothing to judge and passes.

    Zero fails every ``--min-coverage``; the absence of a measurement must not.
    """
    assert coverage_fraction(50, 0) is None
    assert coverage_fraction(50, None) is None
    assert coverage_fraction(None, 100) is None
    assert coverage_fraction(50, -100) is None


def test_a_tie_rounds_up_rather_than_to_even():
    """The row `trf` and `novelty` used to disagree about.

    308/320 is 0.9625 exactly. ``round`` takes the double to 0.963 and
    ``numpy.round`` takes it to 0.962, so the two steps wrote different numbers
    for one insertion of every callset in ``data/sv_output``.
    """
    assert coverage_fraction(308, 320) == 0.963


def _exact_ties(limit: int = 400):
    """``(covered, size)`` pairs whose ratio is exactly a third-decimal tie.

    ``covered / size * 1000`` is an integer plus a half exactly when
    ``2000 * covered`` leaves a remainder of ``size`` against ``2 * size``.
    """
    for size in range(1, limit + 1):
        for covered in range(size + 1):
            if (2000 * covered) % (2 * size) == size:
                yield covered, size


def test_every_tie_rounds_up_whatever_the_denominator():
    """The property a float expression cannot promise.

    Whether ``covered / size`` lands above or below the tie depends on the
    binary approximation to a number neither operand is, so the direction
    varies with the denominator. In exact integer arithmetic it does not.
    """
    ties = list(_exact_ties())
    assert len(ties) > 100, "expected the search to actually find ties"
    for covered, size in ties:
        expected = (covered * 1000 * 2 // size + 1) // 2 / 1000
        assert coverage_fraction(covered, size) == min(expected, 1.0), (
            f"{covered}/{size}")


def test_a_tie_is_where_the_float_expression_actually_differs():
    """Not a hypothetical: `round` and the integer rule disagree on real pairs."""
    disagreed = [(c, s) for c, s in _exact_ties()
                 if round(c / s, 3) != coverage_fraction(c, s)]
    assert disagreed, "no tie distinguishes the two rules; the test proves nothing"


# --------------------------------------------------------------------------- #
# agreement with the vectorised twins
#
# `novelty.insertions` computes both of these over a DataFrame, at catalogue
# scale, and writes the result as `insertion_purity` beside the `repeat_coverage`
# `trf.filters` writes from the scalars here. Both sides were correct on their
# own tests and disagreed with each other at the third decimal -- 308 bases of a
# 320 bp insertion reached one as 0.963 and the other as 0.962 -- which is
# invisible until the same row is scored by two stages.
# --------------------------------------------------------------------------- #

COVERAGE_PAIRS = [(80, 100), (308, 320), (1, 3), (0, 100), (150, 100), (1, 8),
                  (5, 8), (2, 16), (7, 80), (999, 1000), (1, 2000), (13, 160)]


@pytest.mark.parametrize("covered, size", COVERAGE_PAIRS)
def test_coverage_fraction_agrees_with_the_vectorised_twin(covered, size):
    import pandas as pd

    from intruder.pipeline.novelty.insertions import coverage_fractions

    scalar = coverage_fraction(covered, size)
    vectorised = coverage_fractions([covered], [size]).iloc[0]

    if scalar is None:
        assert vectorised is pd.NA
    else:
        assert vectorised is not pd.NA
        assert scalar == vectorised, f"{covered}/{size}: {scalar} != {vectorised}"


def test_the_two_coverage_implementations_agree_on_every_tie():
    """The ties are the only place they can differ, so check all of them."""
    from intruder.pipeline.novelty.insertions import coverage_fractions

    ties = list(_exact_ties(200))
    covered = [c for c, _ in ties]
    sizes = [s for _, s in ties]
    vectorised = coverage_fractions(covered, sizes)

    for i, (c, s) in enumerate(ties):
        assert coverage_fraction(c, s) == vectorised.iloc[i], f"{c}/{s}"


def test_a_missing_denominator_is_missing_on_both_sides():
    import pandas as pd

    from intruder.pipeline.novelty.insertions import coverage_fractions

    assert coverage_fraction(50, 0) is None
    assert coverage_fractions([50], [0]).iloc[0] is pd.NA
    assert coverage_fractions([50], [pd.NA]).iloc[0] is pd.NA


UNION_CASES = [
    [(0, 40), (10, 40), (20, 40)],
    [(0, 20), (60, 90)],
    [(0, 20), (20, 40)],
    [(0, 100), (10, 20), (30, 40)],
    [(10, 10), (30, 20)],
    [(5, 7)],
    [(0, 1), (2, 3), (4, 5)],
]


@pytest.mark.parametrize("intervals", UNION_CASES)
def test_union_length_agrees_with_the_vectorised_twin(intervals):
    """Nothing asserted this before, though both docstrings promise it."""
    import pandas as pd

    from intruder.pipeline.novelty.insertions import union_length as vectorised

    frame = pd.DataFrame({"svid": ["A"] * len(intervals),
                          "rep_start": [s for s, _ in intervals],
                          "rep_end": [e for _, e in intervals]})
    got = int(vectorised(frame, ["svid"], "rep_start", "rep_end").iloc[0])
    assert got == union_length(intervals), intervals
