"""The find/filter pair over this repository's own SV callsets.

Every bug these assertions cover passed the synthetic tests either side of it.
A VCF written for a test carries one caller's conventions, and the ones that
matter here are the conventions two callers *disagree* about: Sniffles declares
an ``SVLEN`` that matches its ``ALT`` and SURVIVOR, merging across 69 samples,
does not; Sniffles emits ``DR``/``DV`` and SURVIVOR emits neither. Both files
are in ``data/sv_output``, so the disagreement is checked rather than imagined.

What is asserted is the table's internal contract -- offsets inside the sequence
they index, columns that agree with each other, a motif the rest of the project
would key the same way -- not particular call counts, which are a property of
the finder's untuned hyperparameters (see :mod:`pipeline.trf.params`) and would
make this a change-detector.
"""

from __future__ import annotations

import csv

import pytest

from intruder.pipeline.trf import filters
from intruder.pipeline.trf.cli import main
from intruder.trcore.motifs import primitive_unit
from intruder.trcore.paths import repo_root
from intruder.trcore.utils.parse import as_size, as_total

pytest.importorskip("pytrf")
pytest.importorskip("cyvcf2")

ROOT = repo_root(__file__)

#: One insertion callset per caller, because the two disagree about the fields
#: this step reads: SURVIVOR's SVLEN does not describe its ALT, and it writes no
#: read-support tags at all.
CALLSETS = {
    "sniffles": "data/sv_output/sniffles/first_500_INS.vcf",
    "survivor": "data/sv_output/survivor_multi_sample_vcf/first_500_INS.vcf",
}

pytestmark = pytest.mark.skipif(
    ROOT is None or not all((ROOT / path).exists() for path in CALLSETS.values()),
    reason="the bundled SV callsets are not available outside a checkout",
)


def insertions(vcf):
    """``(locus, svlen, alt_length)`` per scannable INS record, read independently.

    Deliberately a second, blunter reading of the file rather than a call into
    ``trcore.io``: an assertion about what the reader produced is worth little
    if it is checked against the reader.
    """
    from cyvcf2 import VCF

    for variant in VCF(str(vcf)):
        if variant.INFO.get("SVTYPE") != "INS" or not variant.ALT:
            continue
        alt = variant.ALT[0]
        if alt.startswith("<"):
            continue
        declared = variant.INFO.get("SVLEN")
        declared = declared[0] if isinstance(declared, (list, tuple)) else declared
        yield ((variant.CHROM, str(variant.POS), variant.ID),
               abs(int(declared)),
               len(alt.removeprefix(variant.REF)))


@pytest.fixture(scope="module", params=sorted(CALLSETS))
def calls(request, tmp_path_factory):
    """``(vcf_path, rows)`` for one callset, scanned once for the whole module."""
    vcf = ROOT / CALLSETS[request.param]
    out = tmp_path_factory.mktemp(request.param) / "calls.tsv"
    assert main(["find", str(vcf), str(out)]) == 0
    with open(out, newline="") as handle:
        return vcf, list(csv.DictReader(handle, delimiter="\t"))


def test_the_scan_produces_calls(calls):
    _, rows = calls
    assert rows, "both callsets carry insertions the sv profile should call"


# --------------------------------------------------------------------------- #
# offsets index the sequence that was scanned
# --------------------------------------------------------------------------- #

def test_every_offset_lies_inside_the_insertion_it_is_measured_in(calls):
    """With SVLEN as `insert_size`, 1222 of SURVIVOR's 29858 rows failed this.

    `rep_start`/`rep_end` are offsets into ALT. A merged record keeps one
    representative ALT and a consensus SVLEN, so the two describe different
    strings -- and an offset past its own reported size is the visible end of
    that, with a silently wrong coverage denominator as the useful one.
    """
    _, rows = calls
    for row in rows:
        size = as_size(row["insert_size"])
        assert size is not None
        assert 0 <= int(row["rep_start"]) < int(row["rep_end"]) <= size


def test_insert_size_is_the_length_of_a_scanned_alt(calls):
    """Read back off the VCF, so this is pinned to the file rather than to itself.

    A *set* per locus rather than one value: SURVIVOR writes two insertions of
    441 and 2809 bp under the same CHROM, POS and ID, so an ID does not name one
    insertion in a merged callset.
    """
    vcf, rows = calls
    scanned: dict[tuple, set[int]] = {}
    for locus, _svlen, alt_length in insertions(vcf):
        scanned.setdefault(locus, set()).add(alt_length)

    for row in rows:
        locus = (row["chrom"], row["ins_coord"], row["SVID"])
        assert as_size(row["insert_size"]) in scanned[locus]


def test_the_two_callsets_really_do_disagree_about_svlen(calls):
    """The premise of the test above, asserted rather than assumed.

    If SURVIVOR's SVLEN ever matched its ALT, these tables would be green for a
    reason that has nothing to do with the bug.
    """
    vcf, _ = calls
    mismatched = sum(svlen != alt_length for _, svlen, alt_length in insertions(vcf))

    if "survivor" in str(vcf):
        assert mismatched > 100, "the merged callset is the one that disagrees"
    else:
        assert mismatched == 0, "a single-sample caller's SVLEN describes its ALT"


# --------------------------------------------------------------------------- #
# the columns agree with each other
# --------------------------------------------------------------------------- #

def test_a_merged_callset_reuses_svids_across_insertions(calls):
    """Why `insert_size` is part of `filters.INSERTION_KEYS`.

    Coverage groups a union of offsets and divides by a length. If two
    insertions share a grouping key, offsets measured in one string are divided
    by the other's length. SURVIVOR reuses 257 SVIDs, and one pair collides on
    CHROM and POS too, so the first four columns do not identify an insertion.
    """
    vcf, _ = calls
    seen: dict[tuple, set[int]] = {}
    for locus, _svlen, alt_length in insertions(vcf):
        seen.setdefault(locus, set()).add(alt_length)

    collisions = {locus: sizes for locus, sizes in seen.items() if len(sizes) > 1}
    if "survivor" in str(vcf):
        assert collisions, "the merged callset is the one that collides"
        assert "insert_size" in filters.INSERTION_KEYS
    else:
        assert not collisions


def test_the_derived_columns_follow_from_the_motif_and_the_span(calls):
    _, rows = calls
    for row in rows:
        motif, span = row["motif"], int(row["rep_end"]) - int(row["rep_start"])
        assert int(row["motif_length"]) == len(motif)
        assert int(row["rep_length"]) == span
        assert int(row["rep_units"]) == span // len(motif)


def test_no_motif_is_a_repetition_of_a_shorter_one(calls):
    """`novelty` and `strchive` key every call through `canonical_motif`.

    That reduces to the primitive unit, so a table reporting `GGGG` disagreed
    with both of them about a call all three had in hand -- and carried a
    `motif_length` of 4 for a homopolymer, which `--drop-motif-sizes 1` then
    declined to drop.
    """
    _, rows = calls
    assert [row["motif"] for row in rows
            if primitive_unit(row["motif"]) != row["motif"]] == []


def test_purity_is_a_fraction_above_the_profiles_floor(calls):
    from intruder.pipeline.trf.params import SV

    _, rows = calls
    assert all(SV.min_identity <= float(row["purity"]) <= 1.0 for row in rows)


# --------------------------------------------------------------------------- #
# what each caller supplies, and what it does not
# --------------------------------------------------------------------------- #

def test_depth_is_either_a_readable_total_or_absent(calls):
    """Sniffles writes DR/DV, SURVIVOR writes neither. Both are ordinary inputs.

    What must not appear is a cell that parses to a number nobody measured:
    htslib's missing-integer sentinel written through reads as a depth of
    -4294967296, which fails every `--min-depth` instead of declining to be
    judged by one.
    """
    vcf, rows = calls
    depths = {row["depth"] for row in rows}

    if "survivor" in str(vcf):
        assert depths == {""}, "no DR/DV tags, so no read support to record"
    else:
        assert "" not in depths
    for depth in depths:
        total = as_total(depth) if depth else None
        assert total is None or total >= 0


def test_only_the_scanned_allele_is_written(calls):
    """The reader lifts ALT[0], so allele 1 is the only one a call is evidence about."""
    _, rows = calls
    assert {row["allele"] for row in rows} == {"1"}


def test_rows_are_attributed_to_samples_the_vcf_genotypes(calls):
    from cyvcf2 import VCF

    vcf, rows = calls
    assert {row["sample"] for row in rows} <= set(VCF(str(vcf)).samples)


# --------------------------------------------------------------------------- #
# and the filter step on top
# --------------------------------------------------------------------------- #

def test_the_ins_preset_filters_the_scan_and_its_funnel_adds_up(calls, tmp_path):
    vcf, _ = calls
    scanned = tmp_path / "calls.tsv"
    filtered = tmp_path / "filtered.tsv"
    assert main(["find", str(vcf), str(scanned)]) == 0
    assert main(["filter", str(scanned), str(filtered), "--filter", "ins"]) == 0

    with open(tmp_path / "filtered.tsv.stats.tsv", newline="") as handle:
        funnel = {step["step"]: step for step in csv.DictReader(handle, delimiter="\t")}
    with open(filtered, newline="") as handle:
        kept = list(csv.DictReader(handle, delimiter="\t"))

    assert int(funnel["kept"]["rows"]) == len(kept)
    assert int(funnel["kept"]["rows"]) + int(funnel["kept"]["failed"]) \
        == int(funnel["input"]["rows"])
    assert all(row["filter_status"] == "PASS" for row in kept)


def test_coverage_is_a_fraction_of_the_insertion_it_was_measured_in(calls, tmp_path):
    """The cap at 1.0 used to be load-bearing; now nothing should reach it by error."""
    vcf, _ = calls
    scanned, filtered = tmp_path / "calls.tsv", tmp_path / "filtered.tsv"
    main(["find", str(vcf), str(scanned)])
    main(["filter", str(scanned), str(filtered), "--filter", "none",
          "--min-coverage", "0.8", "--keep-failed"])

    with open(filtered, newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    coverages = [float(row["repeat_coverage"]) for row in rows if row["repeat_coverage"]]
    assert coverages
    assert all(0.0 <= coverage <= 1.0 for coverage in coverages)


def test_repeat_coverage_agrees_with_novelty_row_for_row(calls, tmp_path):
    """The two columns the pipeline compares against the same thresholds.

    `trf filter` writes `repeat_coverage` through :mod:`csv` and
    `novelty annotate` writes `insertion_purity` through pandas, from the same
    three columns of the same table. They are the same quantity under two names,
    so every row has to carry the same number -- and one row of each callset did
    not: 308 bases of a 320 bp insertion is 0.9625, which `round` takes to 0.963
    and `numpy.round` takes to 0.962. Both sides passed their own tests.

    Synthetic rows cannot catch this. A tie at the third decimal needs a
    numerator and denominator that a hand-written fixture has no reason to pick;
    these two callsets supply one each in ~35000 rows. So the agreement is
    checked here, over both files, rather than asserted in a docstring.
    """
    pd = pytest.importorskip("pandas")

    from intruder.pipeline.novelty.insertions import add_insertion_purity

    vcf, _ = calls
    scanned, filtered = tmp_path / "calls.tsv", tmp_path / "filtered.tsv"
    assert main(["find", str(vcf), str(scanned)]) == 0
    assert main(["filter", str(scanned), str(filtered), "--filter", "none",
                 "--min-coverage", "0.8", "--keep-failed"]) == 0

    with open(filtered, newline="") as handle:
        scalar = [row["repeat_coverage"] for row in csv.DictReader(handle, delimiter="\t")]

    vectorised = add_insertion_purity(
        pd.read_csv(scanned, sep="\t"), keys=list(filters.INSERTION_KEYS),
    )["insertion_purity"]

    assert len(scalar) == len(vectorised)
    disagreed = [
        (i, a, b) for i, (a, b) in enumerate(zip(scalar, vectorised))
        if not ((a == "" and pd.isna(b)) or (a != "" and float(a) == float(b)))
    ]
    assert not disagreed, f"{len(disagreed)} row(s) differ, first: {disagreed[:3]}"


def test_the_callsets_really_do_contain_a_third_decimal_tie(calls, tmp_path):
    """Otherwise the agreement test above would pass without exercising anything.

    A tie is where the two rounding rules can differ at all: ``covered / size``
    is exactly a third-decimal half-step when ``2000 * covered`` leaves a
    remainder of ``size`` against ``2 * size``.
    """
    from intruder.trcore.coords import union_length

    vcf, _ = calls
    scanned = tmp_path / "calls.tsv"
    assert main(["find", str(vcf), str(scanned)]) == 0
    with open(scanned, newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    grouped: dict[tuple, list[tuple[int, int]]] = {}
    for row in rows:
        key = tuple(row[column] for column in filters.INSERTION_KEYS)
        grouped.setdefault(key, []).append((int(row["rep_start"]), int(row["rep_end"])))

    ties = 0
    for key, intervals in grouped.items():
        size = as_size(dict(zip(filters.INSERTION_KEYS, key))["insert_size"])
        if size and (2000 * union_length(intervals)) % (2 * size) == size:
            ties += 1
    assert ties, "no insertion lands on a tie; this callset cannot detect the drift"


# --------------------------------------------------------------------------- #
# the FASTA reader, over the same real sequences
#
# Everything above reads a VCF. The other half of this step reads a FASTA, and
# the claim the package is built on is that the finder cannot tell the
# difference -- `trcore.io` normalises both to one `SeqRecord` and `find` scans
# whichever it is handed. Asserting that on a FASTA written for a test proves
# little: the file would carry this package's own conventions. So the sequences
# here are the real insertions, taken from the real callsets and written out as
# FASTA, which makes the two readers comparable on identical bases.
#
# The bundled `first_500_INS.fa` is deliberately not the input: its first line is
# a stray TSV header sitting before any `>`, so pyfastx refuses the whole file
# (see the last test here).
# --------------------------------------------------------------------------- #

def insertion_fasta(vcf, path):
    """The callset's insertions as a FASTA; returns ``[(contig, SVID, length)]``.

    Contig names are made unique, which a merged callset forces: SURVIVOR reuses
    an ``SVID`` across insertions of different lengths, so naming records after
    the ``SVID`` alone writes two different sequences under one contig name --
    and a reader then has no way to say which of them a call sits in.
    """
    from intruder.trcore.io import read_insertions

    records = []
    seen: dict[str, int] = {}
    with open(path, "w") as handle:
        for record in read_insertions(vcf):
            seen[record.name] = seen.get(record.name, 0) + 1
            contig = (record.name if seen[record.name] == 1
                      else f"{record.name}.{seen[record.name]}")
            records.append((contig, record.name, len(record.sequence)))
            handle.write(f">{contig}\n{record.sequence}\n")
    return records


def test_the_two_readers_find_the_same_repeats_in_the_same_sequence(calls, tmp_path):
    """One finder, two sources -- checked on real bases rather than asserted.

    The same insertions, once through cyvcf2 and once through pyfastx, under one
    profile. Every call has to match on offset, motif and purity: the reader
    decides what a record *is*, and must not decide what is in it.

    Grouped by ``SVID`` rather than per record, because the VCF schema has no
    column that distinguishes two insertions sharing one -- which is the same
    ambiguity `insertion_fasta` had to name its way around.
    """
    vcf, _ = calls
    fasta = tmp_path / "insertions.fa"
    records = insertion_fasta(vcf, fasta)
    assert records, "no insertions to compare"
    svid_of = {contig: svid for contig, svid, _ in records}

    from_vcf, from_fasta = tmp_path / "vcf.tsv", tmp_path / "fasta.tsv"
    assert main(["find", "--format", "vcf", "--profile", "sv",
                 str(vcf), str(from_vcf)]) == 0
    assert main(["find", "--format", "fasta", "--profile", "sv",
                 str(fasta), str(from_fasta)]) == 0

    def grouped(path, name_column, start_column, end_column, rename=None):
        found: dict[str, set] = {}
        with open(path, newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                name = row[name_column]
                key = rename[name] if rename else name
                found.setdefault(key, set()).add(
                    (row[start_column], row[end_column], row["motif"], row["purity"]))
        return found

    by_vcf = grouped(from_vcf, "SVID", "rep_start", "rep_end")
    by_fasta = grouped(from_fasta, "contig", "start", "end", rename=svid_of)

    assert set(by_fasta) == set(by_vcf)
    for svid in sorted(by_vcf):
        assert by_fasta[svid] == by_vcf[svid], svid


def test_a_genome_scan_of_real_sequence_keeps_its_rows_and_bed_in_step(calls, tmp_path):
    """The genome schema and its BED export, over real bases.

    `--bed` is re-derived from the written table rather than held in memory, so
    the two can fall out of step without either looking wrong on its own.
    """
    vcf, _ = calls
    fasta = tmp_path / "insertions.fa"
    lengths = {contig: length for contig, _, length in insertion_fasta(vcf, fasta)}

    table, bed = tmp_path / "calls.tsv", tmp_path / "calls.bed"
    assert main(["find", "--format", "fasta", str(fasta), str(table),
                 "--bed", str(bed)]) == 0

    with open(table, newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    with open(bed) as handle:
        intervals = [line.rstrip("\n").split("\t") for line in handle if line.strip()]

    assert rows
    assert len(intervals) == len(rows), "BED and table disagree on how many calls there are"

    for row, interval in zip(rows, intervals):
        assert interval == [row["contig"], row["start"], row["end"], row["motif"]]
        start, end = int(row["start"]), int(row["end"])
        assert 0 <= start < end <= lengths[row["contig"]], row["contig"]
        assert end - start == int(row["rep_length"])
        assert primitive_unit(row["motif"]) == row["motif"]
        assert 0.0 <= float(row["purity"]) <= 1.0


def test_the_bundled_insertion_fasta_scans(tmp_path):
    """The repository's own FASTA, read directly rather than rebuilt from the VCF.

    It carried a stray TSV header on line 1 for a while, before any ``>``, which
    pyfastx refuses -- so the one real FASTA here could not be handed to this
    step at all, while `find.py`'s docstring cited a scan of it. The count below
    is that docstring's number.

    535 of its 6148 records hold no sequence. That is legal, the finder finds
    nothing in them, and the scan reports the tally so it is not mistaken for
    having found fewer repeats.
    """
    fasta = ROOT / "data/sv_output/survivor_multi_sample_vcf/first_500_INS.fa"
    if not fasta.exists():                      # pragma: no cover - checkout only
        pytest.skip("the bundled insertion FASTA is not available")

    table = tmp_path / "calls.tsv"
    assert main(["find", "--format", "fasta", str(fasta), str(table)]) == 0
    with open(table, newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    assert len(rows) == 18166
    assert all(row["genome"] == "first_500_INS" for row in rows)


def test_a_scan_reports_records_that_hold_no_sequence(tmp_path, capsys):
    """Otherwise a mostly-empty FASTA is indistinguishable from a quiet one."""
    fasta = tmp_path / "some_empty.fa"
    fasta.write_text(">a\n\n>b\n" + "CAG" * 20 + "\n>c\n\n")

    assert main(["find", "--format", "fasta", str(fasta), str(tmp_path / "o.tsv")]) == 0
    assert "3 record(s) scanned (2 empty)" in capsys.readouterr().err


def test_a_scan_with_nothing_empty_says_nothing_about_it(tmp_path, capsys):
    fasta = tmp_path / "none_empty.fa"
    fasta.write_text(">a\n" + "CAG" * 20 + "\n")

    assert main(["find", "--format", "fasta", str(fasta), str(tmp_path / "o.tsv")]) == 0

    # The summary line alone: the line above it echoes the output path, and
    # pytest builds tmp_path out of the test's own name.
    summary = [line for line in capsys.readouterr().err.splitlines()
               if "record(s) scanned" in line]
    assert summary == ["[trf] 1 record(s) scanned, 1 call row(s) written"]


def test_an_unreadable_input_is_reported_rather_than_traced(tmp_path, capsys):
    """Every ordinary input failure, as one line naming the file.

    pyfastx refuses anything before the first ``>`` except a blank line -- a
    ``;`` comment from the original Pearson format included -- and cyvcf2
    refuses a file that is not VCF or BCF. Each announced it from inside its own
    internals, so a mistyped path came out as four lines of parser stack saying
    neither which file was wrong nor what it was read as.
    """
    preamble = tmp_path / "with_preamble.fa"
    preamble.write_text("sample\tSVID\trep_start\n>a\n" + "CAG" * 20 + "\n")
    missing_vcf = tmp_path / "absent.vcf"
    not_a_vcf = tmp_path / "prose.vcf"
    not_a_vcf.write_text("this is not a vcf\n")

    for argv, expected in (
        (["find", "--format", "fasta", str(preamble), str(tmp_path / "a.tsv")], "fasta"),
        (["find", "--format", "vcf", str(missing_vcf), str(tmp_path / "b.tsv")], "vcf"),
        (["find", "--format", "vcf", str(not_a_vcf), str(tmp_path / "c.tsv")], "vcf"),
    ):
        assert main(argv) == 1
        err = capsys.readouterr().err
        assert "cannot read" in err and f"as {expected}" in err
        assert not (tmp_path / argv[-1].rsplit("/", 1)[-1]).exists(), \
            "a failed scan left an empty table behind"
