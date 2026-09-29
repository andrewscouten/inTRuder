"""The command line: choosing a reader, and the two commands end to end."""

from __future__ import annotations

import csv

import pytest

from intruder.pipeline.trf.cli import (
    OFF,
    build_parser,
    int_set,
    main,
    optional,
    sniff_format,
)

pytest.importorskip("pytrf")
pytest.importorskip("pyfastx")

FLANK = "GTCCATGAGCTTACGGATCCTAGGCATTCAGGACTTAGCCA"


def write_fasta(tmp_path, **contigs):
    path = tmp_path / "genome.fna"
    with open(path, "w") as handle:
        handle.writelines(f">{name} some organism, complete genome\n{sequence}\n"
                          for name, sequence in contigs.items())
    return path


VCF_HEADER = """\
##fileformat=VCFv4.2
##contig=<ID=chr1,length=100000>
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="type">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="len">
##FORMAT=<ID=GT,Number=1,Type=String,Description="gt">
##FORMAT=<ID=DR,Number=1,Type=Integer,Description="ref depth">
##FORMAT=<ID=DV,Number=1,Type=Integer,Description="alt depth">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3
"""


def write_vcf(tmp_path, insert="CAG" * 40):
    """One insertion carried heterozygously, homozygously, and not at all."""
    path = tmp_path / "sv.vcf"
    row = (f"chr1\t500\tINS0\tA\tA{insert}\t.\tPASS\t"
           f"SVTYPE=INS;SVLEN={len(insert)}\tGT:DR:DV\t"
           f"0/1:20:15\t1/1:2:40\t0/0:30:0\n")
    path.write_text(VCF_HEADER + row)
    return path


def read_tsv(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


# --------------------------------------------------------------------------- #
# format sniffing
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("name, expected", [
    ("x.vcf", "vcf"), ("x.vcf.gz", "vcf"), ("x.bcf", "vcf"),
    ("x.fa", "fasta"), ("x.fna.gz", "fasta"), ("x.fasta", "fasta"),
    ("GCF_000009045.2.fna", "fasta"),
])
def test_sniff_format(name, expected):
    assert sniff_format(name) == expected


def test_sniff_refuses_to_guess():
    """A wrong guess surfaces as a parser error deep inside a reader."""
    with pytest.raises(SystemExit, match="--format"):
        sniff_format("calls.tsv")


def test_optional_accepts_none():
    """"none" is OFF, not None -- None is what argparse leaves for an absent flag."""
    assert optional(int)("none") is OFF
    assert optional(int)("50") == 50


def test_int_set():
    assert int_set("1,2,3") == frozenset({1, 2, 3})
    assert int_set("none") is OFF
    assert int_set("") is OFF


# --------------------------------------------------------------------------- #
# find
# --------------------------------------------------------------------------- #

def test_find_on_a_fasta_writes_genome_rows(tmp_path):
    fasta = write_fasta(tmp_path, contig_1=FLANK + "CAG" * 30 + FLANK)
    out = tmp_path / "calls.tsv"
    assert main(["find", str(fasta), str(out), "--genome", "TOY"]) == 0

    rows = read_tsv(out)
    assert rows, "the planted tract should be called"
    assert rows[0]["genome"] == "TOY"
    assert rows[0]["contig"] == "contig_1"
    assert float(rows[0]["purity"]) >= 0.8


def test_find_names_the_genome_from_the_file_when_not_given(tmp_path):
    fasta = write_fasta(tmp_path, c=FLANK + "CAG" * 30 + FLANK)
    out = tmp_path / "calls.tsv"
    main(["find", str(fasta), str(out)])
    assert read_tsv(out)[0]["genome"] == "genome"


def test_find_scans_every_contig(tmp_path):
    fasta = write_fasta(tmp_path,
                        c1=FLANK + "CAG" * 30 + FLANK,
                        c2=FLANK + "AT" * 45 + FLANK)
    out = tmp_path / "calls.tsv"
    main(["find", str(fasta), str(out)])
    assert {row["contig"] for row in read_tsv(out)} == {"c1", "c2"}


def test_find_writes_bed_alongside(tmp_path):
    fasta = write_fasta(tmp_path, c1=FLANK + "CAG" * 30 + FLANK)
    out, bed = tmp_path / "calls.tsv", tmp_path / "loci.bed"
    main(["find", str(fasta), str(out), "--bed", str(bed)])
    assert len(bed.read_text().splitlines()) == len(read_tsv(out))


def test_bed_is_rejected_for_insertions(tmp_path):
    """An insertion has no interval, so there is nothing to write to BED."""
    fasta = write_fasta(tmp_path, c1=FLANK)
    with pytest.raises(SystemExit, match="--bed"):
        main(["find", str(fasta), str(tmp_path / "o.tsv"),
              "--format", "vcf", "--bed", str(tmp_path / "x.bed")])


def test_finder_flags_override_the_profile(tmp_path):
    """A seed longer than the contig must actually empty the table."""
    fasta = write_fasta(tmp_path, c1=FLANK + "CAG" * 30 + FLANK)
    out = tmp_path / "calls.tsv"
    main(["find", str(fasta), str(out), "--min-seed-length", "500"])
    assert read_tsv(out) == []


def test_min_identity_refuses_a_percent(tmp_path, capsys):
    """`--min-identity 80` reaches the finder as 8000% and quietly finds nothing.

    The 0-1/percent conversion is deliberately invisible from the command line,
    which is exactly what makes 80 the obvious thing to type. Refusing it is
    what keeps the hidden conversion safe.
    """
    fasta = write_fasta(tmp_path, c1=FLANK + "CAG" * 30 + FLANK)
    with pytest.raises(SystemExit):
        main(["find", str(fasta), str(tmp_path / "o.tsv"), "--min-identity", "80"])
    assert "0.8" in capsys.readouterr().err


def test_min_identity_still_accepts_the_bounds_themselves(tmp_path):
    fasta = write_fasta(tmp_path, c1=FLANK + "CAG" * 30 + FLANK)
    for value in ("0", "0.8", "1.0"):
        assert main(["find", str(fasta), str(tmp_path / "o.tsv"),
                     "--min-identity", value]) == 0


# --------------------------------------------------------------------------- #
# find, from a VCF
# --------------------------------------------------------------------------- #

def test_find_on_a_vcf_writes_the_insertion_schema(tmp_path):
    pytest.importorskip("cyvcf2")
    out = tmp_path / "calls.tsv"
    assert main(["find", str(write_vcf(tmp_path)), str(out)]) == 0

    rows = read_tsv(out)
    assert rows
    assert rows[0]["SVID"] == "INS0"
    assert rows[0]["insert_size"] == "120"
    assert rows[0]["motif"] == "CAG"
    assert float(rows[0]["purity"]) == 1.0


def test_calls_are_fanned_out_to_carriers_only(tmp_path):
    """0/1 gives one row, 1/1 two, 0/0 none -- stage 01's existing behaviour."""
    pytest.importorskip("cyvcf2")
    out = tmp_path / "calls.tsv"
    main(["find", str(write_vcf(tmp_path)), str(out)])

    samples = [row["sample"] for row in read_tsv(out)]
    assert samples.count("S1") == 1
    assert samples.count("S2") == 2
    assert "S3" not in samples


def test_depth_is_written_as_variant_then_reference(tmp_path):
    pytest.importorskip("cyvcf2")
    out = tmp_path / "calls.tsv"
    main(["find", str(write_vcf(tmp_path)), str(out)])
    row = next(r for r in read_tsv(out) if r["sample"] == "S1")
    assert row["depth"] == "15,20"          # DV,DR


def test_insertion_offsets_are_inside_the_insertion(tmp_path):
    """Not reference coordinates: ins_coord is the anchor, rep_start is an offset."""
    pytest.importorskip("cyvcf2")
    out = tmp_path / "calls.tsv"
    main(["find", str(write_vcf(tmp_path)), str(out)])
    row = read_tsv(out)[0]
    assert row["ins_coord"] == "500"
    assert int(row["rep_start"]) < int(row["insert_size"])


def test_sv_length_gates_are_applied(tmp_path):
    pytest.importorskip("cyvcf2")
    out = tmp_path / "calls.tsv"
    main(["find", str(write_vcf(tmp_path)), str(out), "--min-sv-length", "500"])
    assert read_tsv(out) == []


def test_max_sv_length_none_lifts_the_ceiling(tmp_path):
    """The documented way to scan every insertion, whatever its size.

    `--max-sv-length none` parses to the OFF sentinel, which is a CLI-layer
    idea; handing it to `read_insertions` as if it were a bound died comparing
    an int to it, so the flag crashed rather than lifting anything.
    """
    pytest.importorskip("cyvcf2")
    out = tmp_path / "calls.tsv"
    vcf = write_vcf(tmp_path, insert="CAG" * 4_000)         # 12 kb, over the 10 kb default

    assert main(["find", str(vcf), str(out), "--max-sv-length", "none"]) == 0
    assert read_tsv(out), "an insertion above the default ceiling should be scanned"

    assert main(["find", str(vcf), str(out)]) == 0
    assert read_tsv(out) == [], "and excluded by the default ceiling"


def test_a_genotype_citing_an_unscanned_alternate_is_not_written(tmp_path):
    """Only ALT[0] is scanned, so only allele 1 is evidence about anything.

    A `1/2` genotype used to emit a row for allele 2 as well, attributing calls
    found in ALT[0] to sequence that was never read -- and counting the sample
    twice in the distinct-(sample, SVID) tallies the reporting scripts run.
    """
    pytest.importorskip("cyvcf2")
    path = tmp_path / "multi.vcf"
    insert = "CAG" * 40
    path.write_text(VCF_HEADER + (
        f"chr1\t500\tINS0\tA\tA{insert},A{'GT' * 60}\t.\tPASS\t"
        f"SVTYPE=INS;SVLEN={len(insert)}\tGT:DR:DV\t1/2:20:15\t1/1:2:40\t0/0:30:0\n"))

    out = tmp_path / "calls.tsv"
    main(["find", str(path), str(out)])
    rows = read_tsv(out)

    assert rows, "the first alternate is still scanned"
    assert {row["allele"] for row in rows} == {"1"}
    assert [row["sample"] for row in rows].count("S1") == 1      # 1/2 carries one copy


def test_a_haploid_carrier_gets_one_row_per_call(tmp_path):
    """A haploid genotype has one allele copy, so it is evidence once.

    cyvcf2 returns `[1, True]` for a bare `1` -- the allele, then the phased
    flag -- and `True == 1`, so reading the first two entries emitted a second
    row for a copy the sample does not carry, with `True` in the numeric
    `allele` column. chrX and chrY calls in a male sample are all haploid.
    """
    pytest.importorskip("cyvcf2")
    path = tmp_path / "haploid.vcf"
    insert = "CAG" * 40
    path.write_text(VCF_HEADER + (
        f"chr1\t500\tINS0\tA\tA{insert}\t.\tPASS\t"
        f"SVTYPE=INS;SVLEN={len(insert)}\tGT:DR:DV\t1:20:15\t1/1:2:40\t0:30:0\n"))

    out = tmp_path / "calls.tsv"
    main(["find", str(path), str(out)])
    rows = read_tsv(out)

    samples = [row["sample"] for row in rows]
    assert samples.count("S1") == 1          # haploid carrier: one copy, one row
    assert samples.count("S2") == 2          # 1/1 still yields two
    assert "S3" not in samples               # haploid reference carries nothing
    assert {row["allele"] for row in rows} == {"1"}, "no allele written as 'True'"


def test_a_sample_without_read_support_has_an_empty_depth(tmp_path):
    """`DR:DV` of `.` must not reach the table as htslib's missing sentinel."""
    pytest.importorskip("cyvcf2")
    path = tmp_path / "missing_depth.vcf"
    insert = "CAG" * 40
    path.write_text(VCF_HEADER + (
        f"chr1\t500\tINS0\tA\tA{insert}\t.\tPASS\t"
        f"SVTYPE=INS;SVLEN={len(insert)}\tGT:DR:DV\t0/1:.:.\t0/1:20:15\t0/0:30:0\n"))

    out = tmp_path / "calls.tsv"
    main(["find", str(path), str(out)])
    depths = {row["sample"]: row["depth"] for row in read_tsv(out)}
    assert depths["S1"] == ""
    assert depths["S2"] == "15,20"

    # And the empty cell is what a depth threshold declines to judge, rather
    # than a negative depth that fails every one of them.
    filtered = tmp_path / "filtered.tsv"
    main(["filter", str(out), str(filtered), "--filter", "ins", "--keep-failed"])
    status = {row["sample"]: row["filter_status"] for row in read_tsv(filtered)}
    assert "depth" not in status["S1"]


# --------------------------------------------------------------------------- #
# insert_size: the length of the sequence the offsets index
# --------------------------------------------------------------------------- #

def test_insert_size_measures_the_alt_not_the_declared_svlen(tmp_path):
    """A merged callset's SVLEN need not describe the one ALT it kept.

    `rep_start`/`rep_end` index ALT, so SVLEN is the wrong denominator for
    repeat coverage and the wrong bound to compare an offset against. SURVIVOR
    disagrees with its own ALT on 297 of the 481 insertions in this
    repository's merged VCF.
    """
    pytest.importorskip("cyvcf2")
    path = tmp_path / "merged.vcf"
    insert = "CAG" * 40                                  # 120 bp of actual sequence
    path.write_text(VCF_HEADER + (
        f"chr1\t500\tINS0\tA\tA{insert}\t.\tPASS\t"
        f"SVTYPE=INS;SVLEN=300\tGT:DR:DV\t0/1:20:15\t0/0:2:0\t0/0:30:0\n"))

    out = tmp_path / "calls.tsv"
    main(["find", str(path), str(out)])
    rows = read_tsv(out)

    assert rows
    assert {row["insert_size"] for row in rows} == {"120"}
    assert all(int(row["rep_end"]) <= int(row["insert_size"]) for row in rows)


def test_coverage_of_a_wholly_repetitive_insertion_is_one(tmp_path):
    """The consequence downstream: an all-repeat insertion must not read as 40%.

    With SVLEN as the denominator this insertion covered 120/300, failing the
    `ins` preset's `--min-coverage 0.8` although every base of the sequence
    scanned is inside a call.
    """
    pytest.importorskip("cyvcf2")
    path = tmp_path / "merged.vcf"
    insert = "CAG" * 40
    path.write_text(VCF_HEADER + (
        f"chr1\t500\tINS0\tA\tA{insert}\t.\tPASS\t"
        f"SVTYPE=INS;SVLEN=300\tGT:DR:DV\t0/1:20:40\t0/0:2:0\t0/0:30:0\n"))

    calls, filtered = tmp_path / "calls.tsv", tmp_path / "filtered.tsv"
    main(["find", str(path), str(calls)])
    main(["filter", str(calls), str(filtered), "--filter", "ins", "--keep-failed"])

    rows = read_tsv(filtered)
    assert rows
    assert {row["repeat_coverage"] for row in rows} == {"1.0"}
    assert all("coverage" not in row["filter_status"] for row in rows)


def test_the_svlen_gate_still_reads_svlen(tmp_path):
    """`--max-sv-length` bounds the variant the caller declared, not the string.

    The two are separate questions, and only one of them moved.
    """
    pytest.importorskip("cyvcf2")
    path = tmp_path / "merged.vcf"
    insert = "CAG" * 40                                  # 120 bp
    path.write_text(VCF_HEADER + (
        f"chr1\t500\tINS0\tA\tA{insert}\t.\tPASS\t"
        f"SVTYPE=INS;SVLEN=300\tGT:DR:DV\t0/1:20:15\t0/0:2:0\t0/0:30:0\n"))

    out = tmp_path / "calls.tsv"
    main(["find", str(path), str(out), "--max-sv-length", "200"])
    assert read_tsv(out) == [], "SVLEN 300 is over the ceiling, whatever ALT measures"


# --------------------------------------------------------------------------- #
# filter
# --------------------------------------------------------------------------- #

def test_filter_drops_failing_rows_and_writes_a_funnel(tmp_path):
    fasta = write_fasta(tmp_path, c1=FLANK + "CAG" * 30 + FLANK)
    calls, filtered = tmp_path / "calls.tsv", tmp_path / "filtered.tsv"
    main(["find", str(fasta), str(calls)])
    assert main(["filter", str(calls), str(filtered),
                 "--filter", "genome", "--min-purity", "1.0"]) == 0

    assert read_tsv(filtered) == []
    funnel = read_tsv(tmp_path / "filtered.tsv.stats.tsv")
    assert {step["step"] for step in funnel} >= {"input", "purity", "kept"}


def test_filter_keep_failed_tags_without_dropping(tmp_path):
    fasta = write_fasta(tmp_path, c1=FLANK + "CAG" * 30 + FLANK)
    calls, filtered = tmp_path / "calls.tsv", tmp_path / "filtered.tsv"
    main(["find", str(fasta), str(calls)])
    main(["filter", str(calls), str(filtered), "--filter", "genome",
          "--min-purity", "1.0", "--keep-failed"])

    rows = read_tsv(filtered)
    assert rows and all(row["filter_status"] == "purity" for row in rows)


def test_filter_rejects_checks_that_cannot_apply(tmp_path):
    """`--filter ins` against a genome table names the columns it needs."""
    fasta = write_fasta(tmp_path, c1=FLANK + "CAG" * 30 + FLANK)
    calls = tmp_path / "calls.tsv"
    main(["find", str(fasta), str(calls)])
    with pytest.raises(SystemExit, match="depth"):
        main(["filter", str(calls), str(tmp_path / "f.tsv"), "--filter", "ins"])


def test_a_genome_scan_streams_rather_than_collecting_every_row(tmp_path):
    """One contig and one row at a time, not a list of the whole scan.

    The ``--bed`` export re-reads the written table instead of keeping the rows
    it just wrote, on the grounds that a genome scan emits more than fits
    comfortably in memory. That was only worth doing if the rows were never
    collected in the first place.
    """
    from intruder.pipeline.trf.cli import _find_in_fasta, _Tally
    from intruder.pipeline.trf.params import GENOME

    fasta = write_fasta(tmp_path, **{f"c{n}": FLANK + "CAG" * 30 + FLANK
                                     for n in range(6)})
    args = build_parser().parse_args(["find", str(fasta), str(tmp_path / "o.tsv")])

    tally = _Tally()
    rows = _find_in_fasta(args, GENOME, tally)
    assert next(iter(rows))["contig"] == "c0"
    assert tally.records < 6, "a row was available before every contig was scanned"


def test_a_missing_parser_leaves_no_empty_table_behind(tmp_path, monkeypatch):
    """A failed run must not look like a scan that found nothing.

    Streaming the rows straight into the writer would open the output and write
    its header before the reader's ImportError surfaced -- see `_primed`.
    """
    from intruder.pipeline.trf import cli as cli_module

    def missing(*_args, **_kwargs):
        raise ImportError("reading VCF needs cyvcf2")

    monkeypatch.setattr(cli_module, "sample_names", missing)
    out = tmp_path / "calls.tsv"
    assert main(["find", "--format", "vcf", str(write_fasta(tmp_path, c1=FLANK)),
                 str(out)]) == 1
    assert not out.exists()


def test_a_missing_parser_exits_cleanly(tmp_path, monkeypatch, capsys):
    """Asking for a format this environment cannot read is an outcome, not a crash.

    The guard's message names the pixi command to install it; a traceback on
    top of that only buries the one line worth reading.
    """
    from intruder.pipeline.trf import cli as cli_module

    def missing(*_args, **_kwargs):
        raise ImportError("reading VCF needs cyvcf2, which is not installed "
                          "by default.\n  in this repo:  pixi run -e pipeline")

    monkeypatch.setattr(cli_module, "sample_names", missing)
    fasta = write_fasta(tmp_path, c1=FLANK)
    assert main(["find", "--format", "vcf", str(fasta), str(tmp_path / "o.tsv")]) == 1

    printed = capsys.readouterr().err
    assert "pixi run -e pipeline" in printed
    assert "Traceback" not in printed


def test_filters_listing_runs(capsys):
    assert main(["filters"]) == 0
    printed = capsys.readouterr().out
    assert "purity" in printed and "ins" in printed


# --------------------------------------------------------------------------- #
# the seam between the parsers and the filters
#
# Both regressions below passed every unit test either side of them:
# `optional("none") is None` and `filters.build(..., {"purity": {"minimum":
# None}})` were each tested in isolation and each correct, and an `ins` table's
# depth column was only ever filtered in tests that wrote it as a bare int. What
# was untested was the handover, so these go through `main`.
# --------------------------------------------------------------------------- #

INSERTION_HEADER = (
    "chrom\tins_coord\tSVID\tdepth\tinsert_size\tsample\tallele\trep_start\t"
    "rep_end\tmotif\tpurity\tmotif_length\trep_length\trep_units\n"
)


def write_insertion_calls(tmp_path, *rows, name="calls.tsv"):
    """An `ins`-schema table, written the way `writers.insertion_rows` writes one."""
    path = tmp_path / name
    with open(path, "w") as handle:
        handle.write(INSERTION_HEADER)
        handle.writelines("\t".join(str(field) for field in row) + "\n" for row in rows)
    return path


def test_depth_is_filtered_on_the_total_of_the_dv_dr_pair(tmp_path):
    """The depth column holds "DV,DR", and the threshold is on the sum.

    `filter_ins_trf.py` summed the pair; a plain numeric parse cannot read it at
    all and returns None, which `Check.fails` passes as "cannot judge" -- so the
    `ins` preset's `depth >= 30` silently kept every row however thin.
    """
    calls = write_insertion_calls(
        tmp_path,
        ("chr1", 100, "INS0", "1,2", 200, "S1", 1, 0, 190, "CAG", 0.95, 3, 190, 63),
        ("chr1", 300, "INS1", "50,40", 200, "S1", 1, 0, 190, "CAG", 0.95, 3, 190, 63),
    )
    out = tmp_path / "filtered.tsv"
    main(["filter", str(calls), str(out), "--filter", "ins", "--keep-failed"])

    status = {row["SVID"]: row["filter_status"] for row in read_tsv(out)}
    assert "depth" in status["INS0"]            # 1 + 2 = 3, under 30
    assert status["INS1"] == "PASS"             # 50 + 40 = 90


def test_none_turns_a_presets_threshold_back_off(tmp_path):
    """`--min-purity none` has to reach `build` as a bound, not vanish en route."""
    calls = write_insertion_calls(
        tmp_path,
        ("chr1", 100, "INS0", "50,40", 200, "S1", 1, 0, 190, "CAG", 0.40, 3, 190, 63),
    )
    kept = tmp_path / "kept.tsv"
    main(["filter", str(calls), str(kept), "--filter", "ins", "--min-purity", "none"])
    assert [row["SVID"] for row in read_tsv(kept)] == ["INS0"]

    # ... and without it the preset's own 0.7 still applies.
    dropped = tmp_path / "dropped.tsv"
    main(["filter", str(calls), str(dropped), "--filter", "ins"])
    assert read_tsv(dropped) == []


def test_an_unreadable_insert_size_is_not_scored_as_zero_coverage(tmp_path):
    """A row nothing can judge passes; only a computed coverage can fail one.

    Scoring an unparseable size as 0.0 coverage inverts the rule every other
    check here follows, and drops exactly the rows the filter knows least about.
    `[200]` is what a VCF list-valued SVLEN looks like once written through, and
    `novelty.insertions.parse_sizes` reads it -- so this side must too.
    """
    calls = write_insertion_calls(
        tmp_path,
        ("chr1", 100, "INS0", "50,40", "[200]", "S1", 1, 0, 190, "CAG", 0.95, 3, 190, 63),
        ("chr1", 300, "INS1", "50,40", "", "S1", 1, 0, 190, "CAG", 0.95, 3, 190, 63),
    )
    out = tmp_path / "filtered.tsv"
    main(["filter", str(calls), str(out), "--filter", "none",
          "--min-coverage", "0.8", "--keep-failed"])

    rows = {row["SVID"]: row for row in read_tsv(out)}
    assert rows["INS0"]["repeat_coverage"] == "0.95"    # 190 / 200
    assert rows["INS1"]["repeat_coverage"] == ""        # no size, no judgement
    assert all(row["filter_status"] == "PASS" for row in rows.values())


def test_refiltering_a_tagged_table_writes_one_status_column(tmp_path):
    calls = write_insertion_calls(
        tmp_path,
        ("chr1", 100, "INS0", "50,40", 200, "S1", 1, 0, 190, "CAG", 0.95, 3, 190, 63),
    )
    once, twice = tmp_path / "once.tsv", tmp_path / "twice.tsv"
    main(["filter", str(calls), str(once), "--filter", "ins", "--keep-failed"])
    main(["filter", str(once), str(twice), "--filter", "ins", "--keep-failed"])

    with open(twice) as handle:
        header = handle.readline().rstrip("\n").split("\t")
    assert header.count("filter_status") == 1
    assert len(header) == len(set(header))


def test_a_table_without_offsets_is_not_silently_emptied(tmp_path):
    """The mirror of the unparseable-size rule, one stage further out.

    A table carrying `insert_size` but no `rep_start`/`rep_end` used to give
    every row a union of zero bases over a readable size, so `--min-coverage`
    scored the lot 0.0 and dropped all of them without a word. Coverage nothing
    could measure has to reach `Check.fails` as "no value to judge".
    """
    calls = tmp_path / "no_offsets.tsv"
    with open(calls, "w", newline="") as handle:
        handle.write("chrom\tins_coord\tSVID\tsample\tinsert_size\tpurity\n")
        handle.write("chr1\t100\tINS0\tS1\t200\t0.95\n")

    out = tmp_path / "filtered.tsv"
    main(["filter", str(calls), str(out), "--filter", "none",
          "--min-coverage", "0.8"])

    rows = read_tsv(out)
    assert len(rows) == 1, "a row nothing could judge was dropped"
    assert rows[0]["repeat_coverage"] == ""
    assert rows[0]["filter_status"] == "PASS"


# --------------------------------------------------------------------------- #
# VCFs that do not carry a genotyped sample
# --------------------------------------------------------------------------- #

#: A sites-only header: the same INFO tags, and no sample columns at all.
SITES_HEADER = """\
##fileformat=VCFv4.2
##contig=<ID=chr1,length=100000>
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="type">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="len">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
"""


def test_a_sites_only_vcf_keeps_its_calls_unattributed(tmp_path):
    """A call is evidence about the insertion whether or not anyone is genotyped.

    Fanning out across samples is how a call reaches the table, so a VCF with
    no sample columns scanned normally, found its repeats and then dropped
    every one of them -- an empty table and exit 0, indistinguishable from a
    file holding nothing repetitive.
    """
    pytest.importorskip("cyvcf2")
    path = tmp_path / "sites_only.vcf"
    insert = "CAG" * 40
    path.write_text(SITES_HEADER + (
        f"chr1\t500\tINS0\tA\tA{insert}\t.\tPASS\tSVTYPE=INS;SVLEN={len(insert)}\n"))

    out = tmp_path / "calls.tsv"
    assert main(["find", str(path), str(out)]) == 0
    rows = read_tsv(out)

    assert rows, "the calls are real; there is simply nobody to attribute them to"
    assert {row["SVID"] for row in rows} == {"INS0"}
    assert {row["insert_size"] for row in rows} == {"120"}
    # Empty, not zero and not a placeholder name: the VCF states no sample fact.
    assert all(row["sample"] == "" and row["allele"] == "" and row["depth"] == ""
               for row in rows)


def test_an_unattributed_row_is_judged_on_what_it_does_carry(tmp_path):
    """The empty cells have to read as "nothing to judge", like a missing depth.

    Otherwise the rows arrive only to be dropped by the next step for lacking
    a sample -- which is the original bug one stage later.
    """
    pytest.importorskip("cyvcf2")
    path = tmp_path / "sites_only.vcf"
    insert = "CAG" * 40
    path.write_text(SITES_HEADER + (
        f"chr1\t500\tINS0\tA\tA{insert}\t.\tPASS\tSVTYPE=INS;SVLEN={len(insert)}\n"))

    calls, filtered = tmp_path / "calls.tsv", tmp_path / "filtered.tsv"
    assert main(["find", str(path), str(calls)]) == 0
    assert main(["filter", str(calls), str(filtered), "--filter", "ins"]) == 0

    rows = read_tsv(filtered)
    assert rows, "a depth threshold must not reject a row for having no sample"
    assert all(row["filter_status"] == "PASS" for row in rows)
    # Coverage is still measurable: it needs the offsets and the size, not a sample.
    assert all(float(row["repeat_coverage"]) >= 0.8 for row in rows)


def test_a_vcf_with_samples_but_no_genotypes_is_reported_not_traced(tmp_path, capsys):
    """cyvcf2 raises a bare `Exception`, which neither handler in `main` caught.

    A caller may report read support per sample and leave the genotype to a
    later step. That is an ordinary input, and it came out as a traceback
    through `writers.insertion_rows` naming neither the file nor the field.
    """
    pytest.importorskip("cyvcf2")
    path = tmp_path / "no_gt.vcf"
    insert = "CAG" * 40
    path.write_text(
        '##fileformat=VCFv4.2\n'
        '##contig=<ID=chr1,length=100000>\n'
        '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="type">\n'
        '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="len">\n'
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="dp">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
        f"chr1\t500\tINS0\tA\tA{insert}\t.\tPASS\t"
        f"SVTYPE=INS;SVLEN={len(insert)}\tDP\t30\n")

    out = tmp_path / "calls.tsv"
    assert main(["find", str(path), str(out)]) == 1
    err = capsys.readouterr().err
    assert "GT" in err and str(path) in err
    assert not out.exists(), "a failed scan left an empty table behind"


# --------------------------------------------------------------------------- #
# which file the error is about
# --------------------------------------------------------------------------- #

def test_an_unwritable_output_names_the_output(tmp_path, capsys):
    """It used to name the input, which is the one file that was fine.

    `exc.filename` is what separates them: Python fills it in on an OSError it
    raised itself, and pyfastx and cyvcf2 leave it empty when they refuse an
    input.
    """
    fasta = write_fasta(tmp_path, c1=FLANK + "CAG" * 30 + FLANK)
    missing = tmp_path / "no_such_dir" / "calls.tsv"

    assert main(["find", str(fasta), str(missing)]) == 1
    err = capsys.readouterr().err
    assert "cannot write" in err and str(missing) in err
    assert "cannot read" not in err


def test_an_unreadable_input_names_the_format_it_was_read_as(tmp_path, capsys):
    """Sniffed or not, the message says what the run actually did.

    `args.format` is still "auto" whenever the format came from the file name,
    so the one case worth explaining was reported as "cannot read x.vcf as
    auto" -- the resolved format is the useful half of that sentence.
    """
    not_a_vcf = tmp_path / "prose.vcf"
    not_a_vcf.write_text("this is not a vcf\n")

    assert main(["find", str(not_a_vcf), str(tmp_path / "o.tsv")]) == 1
    err = capsys.readouterr().err
    assert "as vcf" in err and "as auto" not in err


def test_a_headerless_table_is_not_given_an_invented_schema(tmp_path):
    """An empty file has no columns, and the filter step makes none up.

    It used to write a table whose only column was `filter_status` and report a
    successful run of 0/0 rows -- a header for a schema no scan produces.
    """
    empty = tmp_path / "empty.tsv"
    empty.write_text("")

    out = tmp_path / "filtered.tsv"
    with pytest.raises(SystemExit, match="no header"):
        main(["filter", str(empty), str(out)])
    assert not out.exists()


def test_a_broken_pipe_out_of_the_listing_is_not_an_attribute_error(capsys):
    """`trf filters | head` has no `args.input` for the handler to name.

    The handler written to keep an ordinary OSError tidy raised AttributeError
    from inside itself on the one subcommand that takes no input.
    """
    import builtins

    real = builtins.print

    def closed_stdout(*args, **kwargs):
        if kwargs.get("file") is None:              # the listing writes to stdout
            raise BrokenPipeError(32, "Broken pipe")
        return real(*args, **kwargs)

    builtins.print = closed_stdout
    try:
        assert main(["filters"]) == 1
    finally:
        builtins.print = real
    assert "Broken pipe" in capsys.readouterr().err
