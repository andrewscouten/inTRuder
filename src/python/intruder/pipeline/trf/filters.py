"""Named thresholds on a table of repeat calls, and the presets that group them.

Every filter has a name. That name is what the CLI flag is called, what the
``filter_status`` column says when a row fails, and what the funnel counts --
so a row that was dropped can always be traced back to the one rule that
dropped it, and a threshold can be moved without reading any code.

Rows are tagged rather than deleted as they go. The difference matters: a row
failing three checks accumulates all three reasons, whereas dropping in
sequence attributes it only to whichever check happened to run first, and the
funnel then changes shape when the checks are reordered. Tagging makes the
funnel a summary of the reason column instead of an artefact of evaluation
order, which is what lets a threshold sweep reuse one scan of the input.

This module is pure standard library and works on dicts, because that is what
its input is: a TSV read with :class:`csv.DictReader`. The vectorised
equivalent over a DataFrame is ``novelty.insertions``, at a later stage and a
different scale. The two must not disagree about any of the arithmetic they
share, so none of it is defined here: interval union comes from
:func:`trcore.coords.union_length`, and reading back a cell that is not a
plain number -- ``insert_size``, ``depth`` -- from :mod:`trcore.utils.parse`.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from dataclasses import dataclass, replace

from intruder.trcore.coords import coverage_fraction, union_length
from intruder.trcore.utils.parse import as_number, as_size, as_total

#: ``filter_status`` for a row that failed nothing.
PASS = "PASS"

#: Grouping columns identifying one insertion in one sample.
#:
#: ``insert_size`` is in the key because it is the denominator the union this
#: groups becomes the numerator of, and a merged callset can put two different
#: insertions under the same first four columns: SURVIVOR writes 257 duplicate
#: ``SVID``s into ``data/sv_output/survivor_multi_sample_vcf``, one pair of them
#: sharing a ``CHROM``/``POS`` as well -- a 441 bp insertion and a 2809 bp one at
#: chr1:3502017. Without it their calls union together, so offsets measured in
#: the longer string are divided by the shorter one's length. Grouping by the
#: denominator is what keeps the fraction about a single sequence.
#:
#: ``novelty.cli``'s ``--insertion-key`` default is the same list, and has to
#: stay that way: the two compute the same quantity over the same table.
INSERTION_KEYS = ("chrom", "ins_coord", "SVID", "sample", "insert_size")


@dataclass(frozen=True)
class Check:
    """One named rule about one column.

    Either a numeric bound (``minimum``/``maximum``) or a set membership test
    (``keep``/``drop``); a check with none of them set is inactive and is
    skipped entirely, so "not filtering on depth" and "filtering on depth with
    a threshold of zero" stay distinguishable.

    A missing or unreadable value is not evidence of failure and passes, the
    same rule ``novelty.insertions.filter_reasons`` follows -- a filter should
    remove rows it can judge, not rows it cannot. That holds for a membership
    test as much as for a bound: a ``motif_length`` of ``"?"`` is not a motif
    length outside ``--keep-motif-sizes``, it is an unread cell, and failing it
    drops exactly the rows the filter knows least about.
    """

    name: str
    column: str
    minimum: float | None = None
    maximum: float | None = None
    keep: frozenset | None = None
    drop: frozenset | None = None
    help: str = ""

    @property
    def active(self) -> bool:
        return any(bound is not None
                   for bound in (self.minimum, self.maximum, self.keep, self.drop))

    def fails(self, row: dict) -> bool:
        """Whether ``row`` violates this check."""
        if not self.active:
            return False
        raw = row.get(self.column)
        if raw is None or raw == "":
            return False

        value = READERS.get(self.column, as_number)(raw)
        if value is None:                   # nothing this check can judge
            return False

        if self.keep is not None or self.drop is not None:
            if self.drop is not None and value in self.drop:
                return True
            return self.keep is not None and value not in self.keep

        if self.minimum is not None and value < self.minimum:
            return True
        return self.maximum is not None and value > self.maximum


#: Columns whose cell is not a plain number, and how to read each one. The
#: readers themselves live in :mod:`trcore.utils.parse`, because `novelty` and
#: the isolation forest read the same two columns out of the same tables and a
#: size only one of them can parse is a row they silently disagree about.
READERS = {"depth": as_total, "insert_size": as_size}


# --------------------------------------------------------------------------- #
# the registry
# --------------------------------------------------------------------------- #

# Declared inactive. A preset (below) or a CLI flag supplies the bound; what
# lives here is the name, the column it reads and what it means.
CATALOG: dict[str, Check] = {
    check.name: check
    for check in (
        Check("purity", "purity",
              help="minimum purity of the repeat itself, 0-1 (the finder's own identity)"),
        Check("motif-size", "motif_length",
              help="motif lengths to keep or drop, e.g. --drop-motif-size 1,2,3"),
        Check("rep-units", "rep_units",
              help="minimum whole copies of the motif"),
        Check("rep-length", "rep_length",
              help="minimum span of the repeat in bp"),
        Check("coverage", "repeat_coverage",
              help="minimum fraction of the insertion covered by repeat (union of calls)"),
        Check("insert-size", "insert_size",
              help="maximum insertion size in bp"),
        Check("depth", "depth",
              help="minimum supporting reads"),
    )
}

#: Which named checks apply to calls from each kind of record.
APPLICABLE = {
    "ins": ("purity", "motif-size", "rep-units", "rep-length",
            "coverage", "insert-size", "depth"),
    "genome": ("purity", "motif-size", "rep-units", "rep-length"),
}

# The bounds each preset actually sets. `ins` reproduces what filter_ins_trf.py
# defaulted to, so a run through the new CLI is comparable with the old tables.
PRESETS: dict[str, dict[str, dict]] = {
    "ins": {
        "purity": {"minimum": 0.7},
        "coverage": {"minimum": 0.8},
        "insert-size": {"maximum": 10_000},
        "depth": {"minimum": 30},
    },
    "genome": {
        "purity": {"minimum": 0.7},
    },
    "none": {},
}


def build(preset: str = "none", overrides: dict[str, dict] | None = None) -> list[Check]:
    """The active checks for ``preset``, with per-check ``overrides`` applied.

    An override of ``{"purity": {"minimum": None}}`` turns a preset's threshold
    back off, so ``--filter ins --min-purity none`` is expressible.
    """
    if preset not in PRESETS:
        raise KeyError(f"unknown filter preset {preset!r}; choose from {sorted(PRESETS)}")

    bounds: dict[str, dict] = {name: dict(v) for name, v in PRESETS[preset].items()}
    for name, override in (overrides or {}).items():
        if name not in CATALOG:
            raise KeyError(f"unknown filter {name!r}; choose from {sorted(CATALOG)}")
        bounds.setdefault(name, {}).update(override)

    checks = [replace(CATALOG[name], **bound) for name, bound in bounds.items()]
    return [check for check in checks if check.active]


def inapplicable(checks: Sequence[Check], source: str) -> list[str]:
    """Names of ``checks`` that mean nothing for calls from ``source``.

    Mirrors ``novelty.catalog.RepeatFilter.inapplicable``: a filter on a column
    the table cannot have is a mistake worth naming once and up front, rather
    than one that silently passes every row.
    """
    allowed = APPLICABLE.get(source, ())
    return [check.name for check in checks if check.name not in allowed]


# --------------------------------------------------------------------------- #
# derived columns
# --------------------------------------------------------------------------- #

def add_repeat_coverage(rows: Sequence[dict], keys: Sequence[str] = INSERTION_KEYS,
                        *, start_col: str = "rep_start", end_col: str = "rep_end",
                        size_col: str = "insert_size") -> Sequence[dict]:
    """Add ``repeat_coverage``: repeat bases over insertion length, per insertion.

    The numerator is the *union* of the calls belonging to one insertion. The
    finder emits overlapping calls over the same stretch, so the previous
    per-call ``rep_length / insert_size`` counted shared bases once per call and
    could exceed 1 -- which quietly defeated the ``>= 0.8`` threshold it fed.

    The fraction itself is :func:`trcore.coords.coverage_fraction`, not an
    expression here. ``novelty.insertions.add_insertion_purity`` writes the same
    two integers to ``insertion_purity``, the two columns are compared against
    the same thresholds at different stages, and they used to round the tie at
    the third decimal in opposite directions -- 308 bases of a 320 bp insertion
    reached ``trf`` as 0.963 and ``novelty`` as 0.962, on one row of every
    callset in ``data/sv_output``. Sharing the definition is what makes that
    unrepeatable; the cap at 1 and the rounding rule live there, and
    ``tests/python/pipeline/trf/test_real_data.py`` compares the two columns row
    for row over both callsets rather than taking this paragraph's word for it.

    *Unmeasurable* coverage is empty rather than zero, on either side of the
    division. Zero is a judgement -- it fails every ``--min-coverage`` -- and
    this function has not made one; empty is what :meth:`Check.fails` reads as
    "no value to judge" and passes, which is the rule every other check here
    obeys. That covers an insertion length that will not parse, and equally an
    insertion none of whose offsets will: a table without ``rep_start`` and
    ``rep_end`` at all used to yield a union of zero bases over a perfectly
    readable size, so every row scored 0.0 and ``--min-coverage`` silently
    dropped the lot. An empty group is the absence of a numerator, not a
    measurement that nothing was covered.
    """
    grouped: dict[tuple, list[tuple[int, int]]] = {}
    for row in rows:
        key = tuple(row.get(k) for k in keys)
        start, end = as_number(row.get(start_col)), as_number(row.get(end_col))
        if start is not None and end is not None:
            grouped.setdefault(key, []).append((int(start), int(end)))

    for row in rows:
        key = tuple(row.get(k) for k in keys)
        measured = grouped.get(key)
        covered = union_length(measured) if measured is not None else None
        fraction = coverage_fraction(covered, as_size(row.get(size_col)))
        row["repeat_coverage"] = "" if fraction is None else fraction
    return rows


# --------------------------------------------------------------------------- #
# applying them
# --------------------------------------------------------------------------- #

def reasons(row: dict, checks: Sequence[Check]) -> str:
    """``PASS``, or the names of every check ``row`` failed, comma separated."""
    failed = [check.name for check in checks if check.fails(row)]
    return ",".join(failed) if failed else PASS


def apply(rows: Iterable[dict], checks: Sequence[Check],
          *, column: str = "filter_status") -> tuple[list[dict], list[dict]]:
    """Tag every row with its verdict, and summarise what each check removed.

    Returns ``(tagged, funnel)``. Every input row comes back, carrying
    ``column``; the caller decides whether to write the failures out. The funnel
    counts rows failing each check *independently*, so the numbers do not depend
    on the order the checks are listed in and can overlap -- a row failing both
    ``purity`` and ``depth`` is counted by both, and ``removed`` is the count of
    rows failing at least one.
    """
    tagged = []
    per_check = {check.name: 0 for check in checks}
    removed = 0

    for row in rows:
        failed = [check.name for check in checks if check.fails(row)]
        for name in failed:
            per_check[name] += 1
        if failed:
            removed += 1
        row[column] = ",".join(failed) if failed else PASS
        tagged.append(row)

    funnel = [{"step": "input", "rows": len(tagged), "failed": 0}]
    funnel += [{"step": name, "rows": len(tagged), "failed": count}
               for name, count in per_check.items()]
    funnel.append({"step": "kept", "rows": len(tagged) - removed, "failed": removed})
    return tagged, funnel
