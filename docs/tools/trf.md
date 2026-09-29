# Repeat finding (`trf`)

**Finds tandem repeats in sequence, whatever file the sequence arrived in, and
filters the calls against named thresholds.**

This is stage 01, and everything downstream reads the table it writes.

```bash
pixi run trf find --format vcf  HPRC.ins.vcf  calls.tsv     # SV insertions
pixi run trf find --format fasta genome.fna.gz calls.tsv     # reference contigs
pixi run trf filter --filter ins calls.tsv filtered.tsv      # tag and drop
pixi run trf filters                                         # what can be filtered
```

## One finder, two sources

`pytrf.ATRFinder` takes a name and a string. It has no idea where the string
came from, and that is the design the rest of this step is built around:
[`trcore.io`](../../src/python/intruder/trcore/io/) normalises an SV VCF's
insertions and a reference FASTA's contigs into the same `SeqRecord`, and
`find.py` scans whichever it is handed.

What differs is on either side of the finder:

| | `--format vcf` | `--format fasta` |
|---|---|---|
| reader | `read_insertions` (cyvcf2) | `read_fasta` (pyfastx) |
| origin | `Inserted` | `Genomic` |
| offsets mean | position inside the insertion | position on the contig |
| profile | `sv` | `genome` |
| output schema | `chrom, ins_coord, SVID, depth, insert_size, sample, allele, rep_start, rep_end, …` | `genome, contig, start, end, …` |
| BED export | no — there is no interval | yes, `--bed` |

The schemas are separate on purpose. A call inside an insertion has no
reference interval: `Origin.locate()` returns `None`, because the whole point
of the record is that the reference has nothing there. Writing both kinds into
one table would produce a `start` column that indexes nothing for half its rows.

`insert_size` measures the `ALT` that was scanned, **not** the declared `SVLEN`.
For a single-sample caller the two agree; for a merged one they routinely do
not, because the merged record keeps one representative `ALT` alongside a
consensus length — SURVIVOR disagrees with its own `ALT` on 297 of the 481
insertions in `data/sv_output`, by as much as 117 bp against 53. Since
`rep_start`/`rep_end` index `ALT`, `SVLEN` is the wrong bound to compare an
offset against and the wrong denominator for `repeat_coverage`. `--min-sv-length`
and `--max-sv-length` still gate on `SVLEN`: that is a question about the
variant rather than about the string. A record that declares no `SVLEN` at all
is gated on the length of the sequence it does carry — `SVLEN` is optional in
VCF 4.x, and skipping such a record dropped the insertion while its length was
in hand.

The `motif` column is reduced to its **primitive unit** — the finder reports the
unit it seeded on, which is occasionally a whole-number power of the real one
(`GGGG` for a poly-G, `GTGT` for a `GT` tract). Left alone it makes
`motif_length` a multiple of the period, so `--drop-motif-sizes 1` does not drop
a homopolymer; divides `rep_units` down by that same multiple; and disagrees
with `novelty` and `strchive`, which key every call through
`trcore.motifs.canonical_motif` and had already reduced it.

`--format` defaults to sniffing the file name and refuses to guess when it
cannot tell. `--profile` defaults to matching the source, so naming one is only
needed to override the other.

## Parameters

All six finder settings are untuned hyperparameters, exposed on the CLI and
grouped into profiles in [`params.py`](../../src/python/intruder/pipeline/trf/params.py).

| flag | what it does |
|---|---|
| `--min-motif` / `--max-motif` | shortest / longest motif to look for |
| `--min-identity` | minimum purity of a call, **0–1** |
| `--min-seed-units` / `--min-seed-length` | what it takes to seed a repeat |
| `--max-extension` | bp the finder may extend a seed *through mismatches* |

Two of these have been renamed and one has been corrected. Both changes were
load-bearing:

**`--min-identity` is a fraction, and it is now applied.** ATRFinder measures
identity in percent. `sv_trfcaller.py` passed its `0.8` straight through, so
the finder was told to accept anything above 0.8 **percent** — a gate that
admits everything — while the flag was documented as 80%. Calls at 55% identity
came back from runs nominally filtered at 80%, and all purity control in the
pipeline actually came from the post-hoc `--min-purity` filter. The conversion
now happens in `FinderParams.as_kwargs`, the one place the two vocabularies meet.

Because that conversion is invisible from the command line, `--min-identity 80`
is the obvious thing to type — and it used to reach the finder as 8000%, which
no repeat can satisfy, so the run reported a clean, empty table. Values outside
0–1 are now refused by the argument parser.

**`--max-extension` was `--max_rep_length`, which it never was.** It caps
nothing: with it set to `0`, a 12 bp seed still returns a 90 bp call. It is the
budget for extending a seed through mismatches while identity holds up. The old
name read as a bound on the call and collided with the `rep_length` output
column, which is a different quantity. `--min-seed-length` and
`--min-seed-units` were renamed for the same reason.

These two interact, and the interaction is why they move together. With the
identity gate working, a long extension always decays below the floor before it
stops — `max_extension=10000` plus `min_identity=0.8` returns **nothing at
all**. Smaller values trade tolerance of a decayed tract against boundary
accuracy:

| `--max-extension` | a planted 90 bp CAG tract at 400–490 |
|---|---|
| 0 | `400–490` at 100% — exact, but stops at the first substitution |
| 20 | `381–510` at 84% — catches decay, bleeds ~20 bp into the flank |
| 50 | motif itself destabilises (`CAG` collapses into a 15-mer) |
| 10000 | no calls survive the identity gate |

Both profiles currently use 20. That is provisional, and worth stating
precisely: on synthetic sequence it recovers *that* a tract is there, not the
tract. Three planted 120 bp tracts come back as 160/155/159 bp — inflated by the
full extension at each end — with purity dragged from 1.00 to ~0.85 by the
absorbed background. The only firm result is the upper bound: from about 100 up,
a scan finds nothing at all, because an extension that long always decays
through the identity gate before it stops. Real loci have degenerate edges that
no synthetic tract has, so tuning belongs against real genomes.

> **Stage 01 tables produced before this change are not reproducible with the
> `sv` profile.** They contain calls down to ~55% identity that the corrected
> gate rejects.

## Filtering

Every threshold has a name. That name is the CLI flag, the value in the
`filter_status` column when a row fails, and the row in the funnel — so a
dropped row can always be traced to the rule that dropped it.

```bash
pixi run trf filter --filter ins --min-purity 0.8 --drop-motif-sizes 1,2,3 \
    calls.tsv filtered.tsv
```

| check | column | applies to |
|---|---|---|
| `purity` | `purity` | both |
| `motif-size` | `motif_length` | both |
| `rep-units` / `rep-length` | `rep_units` / `rep_length` | both |
| `coverage` | `repeat_coverage` | insertions |
| `insert-size` | `insert_size` | insertions |
| `depth` | `depth` | insertions |

A VCF with **no sample columns** yields one unattributed row per call rather
than none: `sample`, `allele` and `depth` are left empty. Fanning a call out
across its carriers is how it reaches the table, so a sites-only callset used
to scan normally, find its repeats and drop every one of them — an empty table
that reads as a file holding nothing repetitive. A VCF that has sample columns
but no readable `GT` is a different thing and is refused, naming the record:
the schema attributes a call to the samples carrying it, and there is no way to
say which those are.

`depth` holds `DV,DR` and is filtered on the **total** of the pair, the way
`filter_ins_trf.py` read it. A VCF that does not declare those tags — they are
Sniffles', not the VCF spec's — leaves the column empty rather than failing the
scan, and an empty column is not judged. A sample whose `DR:DV` is `.` in a VCF
that *does* declare them leaves it empty for the same reason: written through,
htslib's missing-integer sentinel reads as a depth of `-4294967296`, which is a
measurement failing every `--min-depth` rather than an absent one.

`--filter` picks a preset (`ins`, `genome`, `none`) and individual flags
override it; `--min-purity none` turns a preset's threshold back off. Asking for
a check the table has no column for fails up front and names it, rather than
silently passing every row.

**A value the filter cannot read is not a failure.** A missing or unparseable
cell passes every check — the rule is that a filter removes rows it can judge,
not rows it cannot. `repeat_coverage` follows it too: an insertion whose size
will not parse gets an empty coverage rather than `0.0`, since zero is a
judgement and would fail every `--min-coverage`. Sizes written as `[415]` by a
list-valued VCF `SVLEN` do parse, as do a signed `-415` and a trailing-unit
`415bp`. That parse is
[`trcore.utils.parse`](../../src/python/intruder/trcore/utils/parse.py), shared
with `novelty.insertions.parse_sizes` so the two cannot disagree — they did,
on all three of those forms, which meant the filter step and the novelty screen
silently differed about which rows had a judgeable insertion length at all.

Rows are **tagged, not dropped as they go**. A row failing three checks
accumulates all three reasons, so the funnel is a summary of the reason column
rather than an artefact of evaluation order, and reordering the checks cannot
change the counts. `--keep-failed` writes every row with its `filter_status`
instead of dropping the failures.

**`coverage` counts the union of the calls over an insertion, not their sum.**
The finder emits overlapping calls over the same stretch — one insertion in the
sample data carries 64 of them — so the previous per-call
`rep_length / insert_size` counted shared bases once per call and could exceed
1, which quietly defeated the `>= 0.8` threshold it fed. The interval sweep is
[`trcore.coords.union_length`](../../src/python/intruder/trcore/coords.py),
shared with `novelty.insertions` so the two cannot disagree.

The union is grouped by `chrom, ins_coord, SVID, sample, insert_size`. The size
is in the key because it is the denominator: a merged callset can put two
insertions under one ID — SURVIVOR writes a 441 bp and a 2809 bp insertion at
chr1:3502017 under `Sniffles2.INS.EDS0` — and grouping those together divides
offsets measured in one string by the other's length. `novelty`'s
`--insertion-key` default is the same list, for the same reason.

## Installing only what you use

Nothing that reads or writes sequence is installed by default. The finder and
both parsers are pixi features, imported when they are *used* rather than when
the package is imported, so `pixi install` builds nothing that needs htslib:

| feature / extra | supplies | needed for |
|---|---|---|
| `trf` | pytrf | finding repeats at all |
| `sv` | cyvcf2 | `--format vcf` |
| `fasta` | pyfastx | `--format fasta` |

The `pipeline` environment carries all three, and is where the full test suite
runs:

```bash
pixi run -e pipeline trf find --format vcf calls.vcf out.tsv
pixi add --feature sv cyvcf2        # or add one to an environment you own
```

Asking `default` for a format it cannot read exits 1 with the command to fix it,
rather than a traceback:

```
$ pixi run -e default trf find --format vcf sv.vcf out.tsv
error: reading VCF needs cyvcf2, which is not installed by default.
  in this repo:  pixi run -e pipeline <your command>
  to add it:     pixi add --feature sv cyvcf2
  installed:     pip install 'intruder[sv]'
```

One caveat: `fasta` is not an independent gate. bioconda's pytrf depends on
`pyfastx >=2.1` (`pixi tree --invert pyfastx`), so any environment with the
finder has pyfastx whether or not `fasta` was asked for. What the feature buys
is a FASTA reader *without* the finder — `trcore.io` on its own.

## Layout

| file | what it does |
|---|---|
| `params.py` | finder settings, and the `sv` / `genome` profiles |
| `find.py` | the finder; records in, `RepeatCall` out |
| `filters.py` | named checks, presets, the coverage correction |
| `writers.py` | the two schemas, plus BED export |
| `cli.py` | `find`, `filter`, `filters` |

`find.py` takes any iterable of `SeqRecord` and opens no files, so it is
testable without touching disk. The package imports no other pipeline step —
only `trcore`, for the record types, coordinate conventions and motif
primitives every step must agree on.
