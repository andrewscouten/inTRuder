"""Named thresholds, the presets grouping them, and the coverage they correct."""

from __future__ import annotations

import pytest

from intruder.pipeline.trf import filters
from intruder.pipeline.trf.filters import PASS, Check
from intruder.trcore.utils.parse import as_size, as_total


def row(**overrides):
    base = {"chrom": "chr1", "ins_coord": 100, "SVID": "INS0", "sample": "S1",
            "purity": 0.9, "motif_length": 4, "rep_units": 10, "rep_length": 40,
            "rep_start": 0, "rep_end": 40, "insert_size": 100, "depth": 50}
    return {**base, **overrides}


# --------------------------------------------------------------------------- #
# checks
# --------------------------------------------------------------------------- #

def test_an_unbounded_check_is_inactive():
    """"Not filtering on depth" and "filtering at zero" stay distinguishable."""
    assert not Check("depth", "depth").active
    assert not Check("depth", "depth").fails(row(depth=0))
    assert Check("depth", "depth", minimum=0).active


def test_a_missing_value_is_not_evidence_of_failure():
    check = Check("depth", "depth", minimum=30)
    assert not check.fails(row(depth=""))
    assert not check.fails({k: v for k, v in row().items() if k != "depth"})


def test_a_non_numeric_value_passes_rather_than_crashing():
    assert not Check("depth", "depth", minimum=30).fails(row(depth="unknown"))


def test_numeric_bounds():
    assert Check("purity", "purity", minimum=0.7).fails(row(purity=0.5))
    assert not Check("purity", "purity", minimum=0.7).fails(row(purity=0.7))
    assert Check("insert-size", "insert_size", maximum=100).fails(row(insert_size=101))


def test_membership_keeps_and_drops():
    drop = Check("motif-size", "motif_length", drop=frozenset({1, 2, 3}))
    assert drop.fails(row(motif_length=2))
    assert not drop.fails(row(motif_length=4))

    keep = Check("motif-size", "motif_length", keep=frozenset({4, 5}))
    assert not keep.fails(row(motif_length=4))
    assert keep.fails(row(motif_length=6))


def test_values_are_compared_as_numbers_not_strings():
    """A TSV read with csv.DictReader hands every column back as text."""
    assert not Check("motif-size", "motif_length", keep=frozenset({4})).fails(
        row(motif_length="4"))


# --------------------------------------------------------------------------- #
# presets
# --------------------------------------------------------------------------- #

def test_build_returns_only_active_checks():
    assert all(check.active for check in filters.build("ins"))


def test_preset_thresholds_are_overridable():
    check = next(c for c in filters.build("ins", {"purity": {"minimum": 0.95}})
                 if c.name == "purity")
    assert check.minimum == 0.95


def test_a_preset_threshold_can_be_turned_off():
    """``--min-purity none`` has to be expressible, not just adjustable."""
    checks = filters.build("ins", {"purity": {"minimum": None}})
    assert not any(c.name == "purity" for c in checks)


def test_unknown_names_are_rejected():
    with pytest.raises(KeyError):
        filters.build("nope")
    with pytest.raises(KeyError):
        filters.build("ins", {"not-a-filter": {"minimum": 1}})


def test_insertion_filters_are_named_inapplicable_to_genome_calls():
    """A filter on a column the table cannot have is a mistake, not a no-op."""
    unusable = filters.inapplicable(filters.build("ins"), "genome")
    assert set(unusable) == {"coverage", "insert-size", "depth"}
    assert filters.inapplicable(filters.build("genome"), "genome") == []


# --------------------------------------------------------------------------- #
# coverage
# --------------------------------------------------------------------------- #

def test_coverage_counts_the_union_not_the_sum():
    """The bug this replaces: overlapping calls double-counted shared bases.

    Three calls over the same 40 bp of a 100 bp insertion cover 40%, not 120%.
    """
    rows = [row(rep_start=0, rep_end=40), row(rep_start=10, rep_end=40),
            row(rep_start=20, rep_end=40)]
    filters.add_repeat_coverage(rows)
    assert all(r["repeat_coverage"] == 0.4 for r in rows)


def test_coverage_never_exceeds_one():
    rows = [row(rep_start=0, rep_end=400, insert_size=100)]
    filters.add_repeat_coverage(rows)
    assert rows[0]["repeat_coverage"] == 1.0


def test_coverage_is_computed_per_insertion():
    rows = [row(SVID="A", rep_start=0, rep_end=50),
            row(SVID="B", rep_start=0, rep_end=10)]
    filters.add_repeat_coverage(rows)
    assert rows[0]["repeat_coverage"] == 0.5
    assert rows[1]["repeat_coverage"] == 0.1


def test_disjoint_calls_add_up():
    rows = [row(rep_start=0, rep_end=20), row(rep_start=60, rep_end=90)]
    filters.add_repeat_coverage(rows)
    assert rows[0]["repeat_coverage"] == 0.5


def test_two_insertions_sharing_an_svid_are_not_one_group():
    """A merged callset can put two insertions under the same chrom/pos/SVID.

    SURVIVOR writes a 441 bp and a 2809 bp insertion at chr1:3502017 under one
    ID. Grouped together, offsets measured in the long one land in the short
    one's union and are divided by its length -- a fraction built from two
    different strings, which the cap at 1.0 then makes look like a clean 100%.
    """
    rows = [row(insert_size=100, rep_start=0, rep_end=20),
            row(insert_size=1_000, rep_start=900, rep_end=1_000)]
    filters.add_repeat_coverage(rows)

    assert rows[0]["repeat_coverage"] == 0.2      # 20 / 100, not (20 + 100) / 100
    assert rows[1]["repeat_coverage"] == 0.1      # 100 / 1000


def test_the_grouping_key_includes_the_denominator():
    """`novelty.cli`'s --insertion-key default has to be this same list.

    The two steps compute the same coverage fraction over the same table, so a
    grouping only one of them applies is a row they disagree about while each
    looks right in isolation.
    """
    assert filters.INSERTION_KEYS == ("chrom", "ins_coord", "SVID", "sample",
                                      "insert_size")

    from intruder.pipeline.novelty.cli import build_parser

    args = build_parser().parse_args(["annotate", "calls.tsv", "out.tsv"])
    assert tuple(args.insertion_key) == filters.INSERTION_KEYS


# --------------------------------------------------------------------------- #
# applying
# --------------------------------------------------------------------------- #

def test_reasons_accumulate_rather_than_short_circuiting():
    checks = filters.build("ins", {"purity": {"minimum": 0.9}, "depth": {"minimum": 60}})
    tagged, _ = filters.apply([row(purity=0.1, depth=1, repeat_coverage=1.0)], checks)
    assert set(tagged[0]["filter_status"].split(",")) >= {"purity", "depth"}


def test_a_passing_row_is_tagged_pass():
    tagged, _ = filters.apply([row(repeat_coverage=1.0)], filters.build("ins"))
    assert tagged[0]["filter_status"] == PASS


def test_every_row_comes_back_tagged():
    """The caller decides what to drop; apply only labels."""
    rows = [row(purity=0.1), row(purity=0.99)]
    tagged, _ = filters.apply(rows, filters.build("genome"))
    assert len(tagged) == 2


def test_funnel_counts_do_not_depend_on_check_order():
    rows = [row(purity=0.1, depth=1, repeat_coverage=1.0) for _ in range(3)]
    checks = filters.build("ins")
    _, forward = filters.apply([dict(r) for r in rows], checks)
    _, reverse = filters.apply([dict(r) for r in rows], list(reversed(checks)))
    assert {s["step"]: s["failed"] for s in forward if s["step"] not in ("input", "kept")} \
        == {s["step"]: s["failed"] for s in reverse if s["step"] not in ("input", "kept")}


def test_funnel_reports_rows_kept():
    rows = [row(purity=0.1), row(purity=0.99)]
    _, funnel = filters.apply(rows, filters.build("genome"))
    kept = next(step for step in funnel if step["step"] == "kept")
    assert kept["rows"] == 1 and kept["failed"] == 1


# --------------------------------------------------------------------------- #
# reading a column that is not a plain number
#
# The readers live in `trcore.utils.parse` and are tested there, against the
# vectorised twin they have to agree with. What belongs here is the wiring:
# which column each one is reached through.
# --------------------------------------------------------------------------- #

def test_the_columns_that_are_not_plain_numbers_route_to_the_shared_readers():
    assert filters.READERS == {"depth": as_total, "insert_size": as_size}


def test_every_other_column_is_read_as_a_plain_number():
    columns = {check.column for check in filters.CATALOG.values()}
    assert columns - set(filters.READERS) == {
        "purity", "motif_length", "rep_units", "rep_length", "repeat_coverage"}


def test_depth_is_read_as_the_total_of_its_parts():
    """`--min-depth` compares the sum of the DV,DR pair, not either half."""
    assert Check("depth", "depth", minimum=30).fails(row(depth="15,10"))       # 25
    assert not Check("depth", "depth", minimum=30).fails(row(depth="15,20"))   # 35


def test_insert_size_is_read_through_the_tolerant_size_parser():
    """`[415]` is what a VCF list-valued SVLEN looks like once written through."""
    check = Check("insert-size", "insert_size", maximum=400)
    assert check.fails(row(insert_size="[415]"))
    assert not check.fails(row(insert_size="[380]"))


def test_a_thin_depth_pair_fails_a_threshold_it_used_to_slip_past():
    assert Check("depth", "depth", minimum=30).fails(row(depth="1,2"))
    assert not Check("depth", "depth", minimum=30).fails(row(depth="15,20"))


def test_coverage_is_empty_rather_than_zero_when_the_size_will_not_parse():
    """Empty passes every check; zero fails every one. Only one of those is a judgement."""
    rows = [row(insert_size="", rep_start=0, rep_end=40),
            row(SVID="INS1", insert_size="[100]", rep_start=0, rep_end=40)]
    filters.add_repeat_coverage(rows)

    assert rows[0]["repeat_coverage"] == ""
    assert rows[1]["repeat_coverage"] == 0.4
    assert not Check("coverage", "repeat_coverage", minimum=0.8).fails(rows[0])


def test_coverage_is_empty_rather_than_zero_when_the_offsets_will_not_parse():
    """The same rule on the other side of the division.

    A table with no ``rep_start``/``rep_end`` at all used to yield a union of
    zero bases over a perfectly readable size, so every row scored 0.0 and
    ``--min-coverage`` silently dropped the lot. An empty group is the absence
    of a numerator, not a measurement that nothing was covered.
    """
    rows = [{"chrom": "chr1", "ins_coord": 10, "SVID": "INS0", "sample": "S1",
             "insert_size": "100", "purity": "0.9"}]
    filters.add_repeat_coverage(rows)

    assert rows[0]["repeat_coverage"] == ""
    assert not Check("coverage", "repeat_coverage", minimum=0.8).fails(rows[0])


def test_one_unreadable_row_does_not_erase_its_insertions_coverage():
    """A group keeps the offsets that *did* parse; only an empty group is unmeasured."""
    rows = [row(rep_start=0, rep_end=40), row(rep_start="", rep_end="")]
    filters.add_repeat_coverage(rows)
    assert all(r["repeat_coverage"] == 0.4 for r in rows)


def test_coverage_rounds_the_way_novelty_does():
    """308 bases of a 320 bp insertion: the row the two steps used to split on.

    `novelty.insertions.add_insertion_purity` writes the same two integers to
    `insertion_purity`, and wrote 0.962 for this row while this one wrote 0.963.
    """
    rows = [row(rep_start=0, rep_end=308, insert_size=320)]
    filters.add_repeat_coverage(rows)
    assert rows[0]["repeat_coverage"] == 0.963
