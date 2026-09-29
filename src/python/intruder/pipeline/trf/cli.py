"""Command line for finding tandem repeats and filtering the calls.

``pixi install`` installs this as the ``trf`` command; ``python -m
intruder.pipeline.trf`` is the same program without installing anything.

    # what can be filtered, and what each preset sets
    pixi run trf filters

    # SV insertions out of a VCF (the old `svpytrf` behaviour)
    pixi run trf find --format vcf HPRC.ins.vcf calls.tsv

    # every contig of a bacterial reference, plus a BED for the read filter
    pixi run trf find --format fasta --genome GCF_000009045 \\
        genome.fna.gz calls.tsv --bed loci.bed

    # tag calls against a preset, overriding one threshold
    pixi run trf filter --filter ins --min-purity 0.8 calls.tsv filtered.tsv

``--format`` chooses which :mod:`trcore.io` reader produces the records and
``--profile`` chooses the finder settings sized for them; both default to
matching the input, so naming one is only necessary to override the other.
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections.abc import Iterator
from dataclasses import dataclass
from itertools import chain
from pathlib import Path

from intruder.trcore.io import read_fasta, read_insertions, sample_names

from . import filters
from .find import find_in_record, find_repeats
from .params import DEFAULT_PROFILE, PROFILES
from .writers import (
    GENOME_COLUMNS,
    INSERTION_COLUMNS,
    genome_row,
    insertion_rows,
    write_bed,
    write_tsv,
)

FORMATS = ("auto", "vcf", "fasta")

_VCF_SUFFIXES = (".vcf", ".bcf")
_FASTA_SUFFIXES = (".fa", ".fna", ".fasta", ".fas", ".ffn")


def sniff_format(path: str | Path) -> str:
    """``vcf`` or ``fasta`` from the file name, ignoring a ``.gz``.

    Named after ``novelty.platforms.sniff_format``, which does the same job for
    catalogue layouts. Raises rather than guessing: a misread format produces a
    confusing parser error deep inside a reader instead of a clear one here.
    """
    suffixes = [s.lower() for s in Path(path).suffixes]
    if suffixes and suffixes[-1] in (".gz", ".bgz"):
        suffixes.pop()
    tail = suffixes[-1] if suffixes else ""
    if tail in _VCF_SUFFIXES:
        return "vcf"
    if tail in _FASTA_SUFFIXES:
        return "fasta"
    raise SystemExit(
        f"cannot tell what {Path(path).name!r} is from its name; "
        f"pass --format vcf or --format fasta"
    )


class _Off:
    """``--min-purity none``: the flag *was* given, and it turns the bound off.

    Distinct from ``None``, which is what argparse leaves behind for a flag
    nobody passed. Collapsing the two is what made ``--min-purity none``
    unexpressible: :func:`run_filter` drops the flags that were not given, so a
    "none" parsed to ``None`` was dropped along with them and the preset's own
    threshold survived the very flag written to remove it.
    """

    def __repr__(self) -> str:                  # pragma: no cover - debugging aid
        return "none"


#: The sentinel :func:`optional` and :func:`int_set` return for an explicit "none".
OFF = _Off()

_OFF_WORDS = ("none", "off")


def optional(kind):
    """An argparse type accepting ``none`` to mean "no bound", as :data:`OFF`."""
    def parse(text: str):
        return OFF if text.strip().lower() in _OFF_WORDS else kind(text)
    parse.__name__ = getattr(kind, "__name__", "value")
    return parse


def int_set(text: str) -> frozenset[int] | _Off:
    """``1,2,3`` -> ``{1, 2, 3}``; empty or ``none`` -> :data:`OFF`."""
    stripped = text.strip()
    if not stripped or stripped.lower() in _OFF_WORDS:
        return OFF
    return frozenset(int(part) for part in stripped.split(",") if part.strip())


def unset(value):
    """:data:`OFF` -> ``None``, for handing a bound to something outside this layer.

    :data:`OFF` only means anything to :func:`run_filter`, which has to tell
    "flag absent" from "flag set to none". Everywhere else it is an ordinary
    "no bound" and has to arrive as ``None``: ``--max-sv-length none`` reached
    ``trcore.io.read_insertions`` as the sentinel itself and died comparing an
    int to it, so the one documented way to lift the SVLEN ceiling was the one
    that could not run.
    """
    return None if value is OFF else value


def fraction(text: str) -> float:
    """An argparse type for a purity in ``0..1``, rejecting a percent.

    The scale conversion :meth:`params.FinderParams.as_kwargs` performs is
    invisible from the command line, so ``--min-identity 80`` looks like the
    obvious way to ask for 80%. It reaches the finder as 8000%, which no repeat
    can satisfy, and the run reports a clean, empty table. Refusing the value
    is what makes the conversion safe to keep hidden.
    """
    value = float(text)
    if not 0.0 <= value <= 1.0:
        raise argparse.ArgumentTypeError(
            f"{text} is not a purity in 0-1. This is converted to the finder's "
            f"percent scale for you, so 80 percent is --min-identity 0.8")
    return value


# --------------------------------------------------------------------------- #
# find
# --------------------------------------------------------------------------- #

@dataclass
class _Tally:
    """Records a streaming scan has reached, counted as its rows are consumed.

    The scan is written down as a generator so a genome run holds one contig
    and one row at a time, which means the record count is not known until the
    writer has drained it. Passing a tally in is what lets the summary line
    still report it.
    """

    records: int = 0

    #: Records carrying no sequence at all. Legal, and the finder correctly finds
    #: nothing in them -- but a scan whose input is 9% empty records looks
    #: identical to one that simply found fewer repeats, and this repository's own
    #: `first_500_INS.fa` is exactly that: 535 of its 6148 records are empty.
    #: Counting them is what makes the difference visible without dividing one
    #: number by another. Only a FASTA can produce them; `read_insertions` drops
    #: an insertion with no sequence before it becomes a record.
    empty: int = 0


def _primed(rows: Iterator[dict]) -> Iterator[dict]:
    """Pull the first row now, before the caller opens an output file.

    The readers import their parsers on first use, so an environment without
    cyvcf2 or pyfastx only discovers it when a record is actually read. Handing
    a lazy generator straight to :func:`write_tsv` would move that discovery to
    after the output has been created and its header written, turning "this
    needs cyvcf2" into the same message plus an empty table that looks like a
    scan which found nothing.
    """
    rows = iter(rows)
    try:
        first = next(rows)
    except StopIteration:
        return iter(())
    return chain([first], rows)


def run_find(args) -> int:
    source = args.format if args.format != "auto" else sniff_format(args.input)
    profile = args.profile if args.profile != "auto" else DEFAULT_PROFILE[source]

    # Written back so that the error handler in `main` names the format the run
    # actually used. It read `args.format`, which is still "auto" whenever the
    # format was sniffed -- so the one case where the reader's complaint is
    # worth explaining was the one reported as "cannot read x.vcf as auto".
    args.format = source

    # Checked before anything is read: a genome scan is minutes of work, and a
    # flag combination that cannot succeed should not cost them.
    if args.bed and source != "fasta":
        raise SystemExit("--bed needs reference intervals; it is only valid with "
                         "--format fasta (an insertion has no interval to write)")

    params = PROFILES[profile].replace(
        min_motif=args.min_motif,
        max_motif=args.max_motif,
        min_identity=args.min_identity,
        min_seed_units=args.min_seed_units,
        min_seed_length=args.min_seed_length,
        max_extension=args.max_extension,
        quiet=False if args.finder_errors else None,
    )

    print(f"[trf] {source} -> {args.output} (profile {profile})", file=sys.stderr)
    tally = _Tally()
    if source == "vcf":
        rows, columns = _find_in_vcf(args, params, tally), INSERTION_COLUMNS
    else:
        rows, columns = _find_in_fasta(args, params, tally), GENOME_COLUMNS

    written = write_tsv(args.output, columns, _primed(rows))
    counted = f"{tally.records} record(s) scanned"
    if tally.empty:
        counted += f" ({tally.empty} empty)"
    print(f"[trf] {counted}, {written} call row(s) written", file=sys.stderr)

    if args.bed:
        # Re-derived rather than held: a genome scan can emit far more rows than
        # fit comfortably in memory, so the table on disk is the source.
        with open(args.output, newline="") as handle:
            n_bed = write_bed(args.bed, csv.DictReader(handle, delimiter="\t"))
        print(f"[trf] {n_bed} BED interval(s) -> {args.bed}", file=sys.stderr)
    return 0


def _find_in_vcf(args, params, tally: _Tally) -> Iterator[dict]:
    """Calls from SV insertions, fanned out across the samples carrying them.

    ``sample_names`` runs before the generator is built, not inside it, so a
    missing cyvcf2 is raised here rather than on the first row -- see
    :func:`_primed`.
    """
    samples = sample_names(args.input)
    records = read_insertions(args.input, min_length=args.min_sv_length,
                              max_length=unset(args.max_sv_length), with_variant=True)

    def rows():
        for record, variant in records:
            tally.records += 1
            for call in find_in_record(record, params):
                yield from insertion_rows(call, variant, samples)
    return rows()


def _genome_name(path: str | Path) -> str:
    """An assembly name from a FASTA's file name, keeping the accession whole.

    Strips a compression suffix and then the FASTA suffix, and nothing else:
    splitting on the first ``.`` instead turns
    ``GCF_000009045.2_ASM904v1_genomic.fna.gz`` into ``GCF_000009045`` and drops
    the assembly version, which is the part that distinguishes two builds of
    the same accession in the manifest.
    """
    name = Path(path).name
    for suffix in (".gz", ".bgz"):
        if name.lower().endswith(suffix):
            name = name[: -len(suffix)]
            break
    stem = Path(name).stem
    return stem or name


def _find_in_fasta(args, params, tally: _Tally) -> Iterator[dict]:
    """Calls on reference contigs, one row each.

    Streamed rather than collected. A whole-genome scan emits far more rows
    than a list of them wants to be -- the ``--bed`` export below already
    re-reads the written table for exactly that reason, which was worth doing
    only if the table was never in memory to begin with.
    """
    genome = args.genome or _genome_name(args.input)
    records = read_fasta(args.input)

    def scan():
        for record in records:
            tally.records += 1
            if not len(record):
                tally.empty += 1
            yield record

    return (genome_row(call, genome) for call in find_repeats(scan(), params))


# --------------------------------------------------------------------------- #
# filter
# --------------------------------------------------------------------------- #

def run_filter(args) -> int:
    overrides = {
        "purity": {"minimum": args.min_purity},
        "rep-units": {"minimum": args.min_rep_units},
        "rep-length": {"minimum": args.min_rep_length},
        "coverage": {"minimum": args.min_coverage},
        "insert-size": {"maximum": args.max_insert_size},
        "depth": {"minimum": args.min_depth},
        "motif-size": {"keep": args.keep_motif_sizes, "drop": args.drop_motif_sizes},
    }
    # Two passes, and the order matters: drop the flags argparse defaulted to
    # None because nobody passed them, *then* turn the ones that were passed as
    # "none" back into the None that `filters.build` reads as "no bound".
    overrides = {name: {k: unset(v) for k, v in bound.items() if v is not None}
                 for name, bound in overrides.items()}
    checks = filters.build(args.filter, {k: v for k, v in overrides.items() if v})

    with open(args.input, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
        columns = list(reader.fieldnames or ())

    if not columns:
        # An empty file has no schema, and the filter step invents none. It used
        # to write a table whose only column was `filter_status` and exit 0 --
        # a header for a schema no scan produces, reported as a successful run.
        raise SystemExit(f"{args.input}: no header, so there is no table to filter")

    source = "ins" if "insert_size" in columns else "genome"
    unusable = filters.inapplicable(checks, source)
    if unusable:
        raise SystemExit(
            f"these filters do not apply to {source} calls: {', '.join(unusable)}"
        )

    if any(check.name == "coverage" for check in checks):
        filters.add_repeat_coverage(rows)
        if "repeat_coverage" not in columns:
            columns.append("repeat_coverage")

    tagged, funnel = filters.apply(rows, checks)
    kept = tagged if args.keep_failed else [r for r in tagged
                                            if r["filter_status"] == filters.PASS]
    write_tsv(args.output, columns, kept, extras=["filter_status"])

    stats = args.stats_output or f"{args.output}.stats.tsv"
    write_tsv(stats, ("step", "rows", "failed"), funnel)
    print(f"[trf] {len(kept)}/{len(tagged)} row(s) kept -> {args.output}; "
          f"funnel -> {stats}", file=sys.stderr)
    return 0


# --------------------------------------------------------------------------- #
# filters
# --------------------------------------------------------------------------- #

def run_filters(args) -> int:
    print("checks:")
    for name, check in filters.CATALOG.items():
        print(f"  {name:<12} {check.column:<16} {check.help}")
    print("\npresets (--filter):")
    for preset, bounds in filters.PRESETS.items():
        if not bounds:
            print(f"  {preset:<12} (nothing set)")
            continue
        described = ", ".join(
            f"{name} {'>=' if 'minimum' in bound else '<='}"
            f"{bound.get('minimum', bound.get('maximum'))}"
            for name, bound in bounds.items()
        )
        print(f"  {preset:<12} {described}")
    print("\napplicable by source:")
    for source, names in filters.APPLICABLE.items():
        print(f"  {source:<12} {', '.join(names)}")
    return 0


# --------------------------------------------------------------------------- #
# wiring
# --------------------------------------------------------------------------- #

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="trf",
        description="Find tandem repeats in sequence from any source, and filter the calls.",
    )
    subcommands = parser.add_subparsers(dest="command", required=True)

    find = subcommands.add_parser(
        "find", help="run the repeat finder over a VCF's insertions or a FASTA's contigs",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    find.add_argument("input", help="SV VCF or reference FASTA (optionally gzipped)")
    find.add_argument("output", help="TSV of repeat calls")
    find.add_argument("--format", choices=FORMATS, default="auto",
                      help="which reader produces the records")
    find.add_argument("--profile", choices=("auto", *PROFILES), default="auto",
                      help="finder settings sized for the input")
    find.add_argument("--genome", default=None,
                      help="assembly name for the genome column (default: from the file name)")
    find.add_argument("--bed", default=None,
                      help="also write BED4 of the loci (FASTA only)")
    find.add_argument("--finder-errors", action="store_true",
                      help="let the finder write to stderr instead of silencing it")

    sv = find.add_argument_group("VCF options")
    sv.add_argument("--min-sv-length", type=int, default=50, help="minimum SVLEN")
    sv.add_argument("--max-sv-length", type=optional(int), default=10_000,
                    help="maximum SVLEN, or none")

    finder = find.add_argument_group(
        "finder options (override the profile; all are untuned hyperparameters)")
    finder.add_argument("--min-motif", type=int, default=None,
                        help="shortest motif to look for")
    finder.add_argument("--max-motif", type=int, default=None,
                        help="longest motif to look for")
    finder.add_argument("--min-identity", type=fraction, default=None,
                        help="minimum purity of a call, 0-1 (converted to the finder's %%)")
    finder.add_argument("--min-seed-units", type=int, default=None,
                        help="exact copies needed to seed a repeat")
    finder.add_argument("--min-seed-length", type=int, default=None,
                        help="bp needed to seed a repeat")
    finder.add_argument("--max-extension", type=int, default=None,
                        help="bp the finder may extend a seed through mismatches; "
                             "trades boundary accuracy against catching a decayed tract")
    find.set_defaults(handler=run_find)

    filter_cmd = subcommands.add_parser(
        "filter", help="tag and drop calls against named thresholds",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    filter_cmd.add_argument("input", help="TSV of repeat calls")
    filter_cmd.add_argument("output", help="filtered TSV")
    filter_cmd.add_argument("--filter", choices=tuple(filters.PRESETS), default="none",
                            help="preset of named thresholds to start from")
    filter_cmd.add_argument("-s", "--stats-output", default=None,
                            help="funnel TSV (default: <output>.stats.tsv)")
    filter_cmd.add_argument("--keep-failed", action="store_true",
                            help="write every row with its filter_status, dropping none")
    filter_cmd.add_argument("--min-purity", type=optional(float), default=None)
    filter_cmd.add_argument("--min-rep-units", type=optional(int), default=None)
    filter_cmd.add_argument("--min-rep-length", type=optional(int), default=None)
    filter_cmd.add_argument("--min-coverage", type=optional(float), default=None)
    filter_cmd.add_argument("--max-insert-size", type=optional(int), default=None)
    filter_cmd.add_argument("--min-depth", type=optional(int), default=None)
    filter_cmd.add_argument("--keep-motif-sizes", type=int_set, default=None,
                            help="motif lengths to keep exclusively, e.g. 4,5,6")
    filter_cmd.add_argument("--drop-motif-sizes", type=int_set, default=None,
                            help="motif lengths to drop, e.g. 1,2,3")
    filter_cmd.set_defaults(handler=run_filter)

    listing = subcommands.add_parser("filters", help="list the named checks and presets")
    listing.set_defaults(handler=run_filters)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        return args.handler(args)
    except ImportError as exc:
        # The parsers and the finder are optional, so a missing one is an
        # ordinary outcome of asking for a format this environment cannot read,
        # not a crash. The guards in `trcore.io` and `find` already say which
        # command installs it; a traceback on top only buries that.
        print(f"error: {exc}", file=sys.stderr)
        return 1
    except (OSError, RuntimeError) as exc:
        # An unreadable input is the other ordinary failure, and it used to be
        # the only one that came out as a traceback. A file that is missing, or
        # not the format it was said to be, is a thing the person running this
        # can fix -- but each parser announces it in its own vocabulary, from
        # inside its own internals: pyfastx raises `RuntimeError: ... is not
        # plain or gzip compressed fasta formatted file`, and `FileExistsError`
        # (an OSError, despite the name) when the file is absent, while cyvcf2
        # raises `OSError: ... is not valid bcf or vcf (format: 2 mode: r)`.
        # Four lines of pyfastx or cyvcf2 stack say less about which file was
        # wrong than one line naming it does.
        #
        # Deliberately narrow: these two are what the readers raise about their
        # input. A TypeError or KeyError from this package is a bug and should
        # still come out as a traceback rather than a tidy message.
        #
        # Which file went wrong is `exc.filename`, which Python fills in on an
        # OSError it raised itself and the parsers leave empty. So an output
        # this process could not create names the output, while pyfastx or
        # cyvcf2 refusing an input names the input -- a mistyped output
        # directory used to be reported as "cannot read <input>", which sends
        # the reader to the one file that was fine. `args.input` is read
        # defensively for the same reason: the `filters` subcommand takes no
        # input at all, and a broken pipe out of it raised AttributeError from
        # inside the handler written to keep it tidy.
        source = getattr(args, "input", None)
        target = getattr(exc, "filename", None)
        if target is not None and str(target) != str(source):
            print(f"error: cannot write {target}: {exc}", file=sys.stderr)
        elif source is None:
            print(f"error: {exc}", file=sys.stderr)
        else:
            print(f"error: cannot read {source} as {getattr(args, 'format', 'input')}: "
                  f"{exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
