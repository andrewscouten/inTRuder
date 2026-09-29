"""Reading a table cell that is not a plain number.

The agreement tests at the bottom are the point of the module. Both parsers
were correct in isolation and disagreed with each other, which is invisible
until two steps score the same row differently.
"""

from __future__ import annotations

import pytest

from intruder.trcore.utils.parse import as_number, as_size, as_total

# --------------------------------------------------------------------------- #
# as_number
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("raw, expected", [
    ("40", 40),
    ("-40", -40),
    ("0.75", 0.75),
    (" 40 ", 40),
    (40, 40),
    (0.75, 0.75),
    ("", None),
    ("unknown", None),
    (None, None),
])
def test_as_number(raw, expected):
    assert as_number(raw) == expected


def test_integral_text_stays_an_int():
    """A motif length read back out of a TSV has to match the int in a keep set."""
    assert as_number("4") in frozenset({4})


# --------------------------------------------------------------------------- #
# as_size
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("raw, expected", [
    ("415", 415),
    ("[415]", 415),             # a VCF list-valued SVLEN written straight through
    ("415.0", 415),
    (415, 415),
    ("-415", 415),              # SVLEN may be written signed; a length is not
    ("[-415]", 415),
    ("415bp", 415),
    ("[415,200]", 415),         # the first entry is this sample's insertion
    ("", None),
    ("[]", None),
    ("unknown", None),
    ("nan", None),              # not a length, and must not become one
])
def test_as_size(raw, expected):
    assert as_size(raw) == expected


def test_a_signed_size_does_not_make_coverage_negative():
    """The reason for the ``abs``: a negative denominator fails every threshold."""
    assert as_size("-415") > 0


# --------------------------------------------------------------------------- #
# as_total
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("raw, expected", [
    ("15,20", 35),              # DV,DR as the trf writer emits it
    ("50", 50),                 # already collapsed by an earlier filter pass
    (50, 50),
    (" 15 , 20 ", 35),
    ("", None),
    ("unknown", None),
    ("15,unknown", None),       # half a pair is not two thirds of an answer
])
def test_as_total(raw, expected):
    assert as_total(raw) == expected


# --------------------------------------------------------------------------- #
# agreement with the vectorised twin
#
# `novelty.insertions.parse_sizes` reads the same column over a DataFrame. The
# two drifted: the scalar side returned -415 where this one returned 415, and
# None where this one recovered 415 from "[415,200]" and "415bp". Both looked
# right on their own tests; what nothing asserted was that they match, so the
# filter step and the novelty screen quietly disagreed about which rows had a
# judgeable insertion length at all.
# --------------------------------------------------------------------------- #

SIZE_FORMS = ["415", "[415]", "415.0", "-415", "[-415]", "415bp", "[415,200]",
              "", "[]", "unknown", "nan"]


@pytest.mark.parametrize("raw", SIZE_FORMS)
def test_as_size_agrees_with_the_vectorised_parse_sizes(raw):
    import pandas as pd

    from intruder.pipeline.novelty.insertions import parse_sizes

    scalar = as_size(raw)
    vectorised = parse_sizes([raw]).iloc[0]

    if scalar is None:
        assert vectorised is pd.NA, f"{raw!r}: scalar cannot read it, vectorised can"
    else:
        assert vectorised is not pd.NA, f"{raw!r}: vectorised cannot read it, scalar can"
        assert scalar == vectorised, f"{raw!r}: {scalar} != {vectorised}"
