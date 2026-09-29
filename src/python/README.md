# Python source (`src/python`)

Python code for the project. Environments are managed with [pixi](https://pixi.sh);
`[tool.pixi]` in `../../pyproject.toml` declares all of them.

## Layout

Everything lives under one installed package, `intruder`:

```
src/python/intruder/
├── trcore/       shared primitives
│   ├── coords/motifs/flanks   the domain: what this project means by a repeat
│   ├── io/          read sequence off disk — FASTA, VCF insertions, one record type
│   └── utils/       infrastructure that knows nothing about repeats — parse, fetch
├── pipeline/     the pipeline steps
│   ├── trf/         call repeats inside SV insertions, then filter the calls
│   ├── novelty/     is this repeat absent from the reference and the catalogues?
│   ├── strchive/    is it a known disease locus?
│   ├── compression/ how compressible is the insertion? a cheap repetitiveness proxy
│   └── annotation/  genic / clinical context
└── analysis/     post-hoc analysis of pipeline output — never a pipeline step
```

Two rules hold this together:

- **Steps do not import each other.** Each reads files and writes files, so a
  Nextflow process, a shell loop and a notebook can all drive it unchanged.
  Anything genuinely common goes in `trcore`, which is the only cross-cutting
  import allowed.
- **One owned top-level name.** `intruder` is what gets installed into
  site-packages. Bare `pipeline`, `analysis` or `filter` packages would collide
  with real distributions on PyPI.

Shell scripts live in `../../scripts/`, not here. R lives in `../R/`.

## Setup

```bash
pixi install                  # build the default environment from the lockfile
```

## Environments

The default environment is the pipeline plus ruff and pytest, and nothing else.
Everything heavy is a separate environment, so a Nextflow worker that only writes
a TSV never installs a plotting stack — and installing one environment can never
disturb another:

| Environment | Install | What it covers |
|---|---|---|
| `default` | `pixi install` | every pipeline step — the core (pytrf, pandas, numpy, tqdm) plus the `sv`, `compression`, `search`, `annotation` and `modeling` features — plus ruff and pytest |
| `modeling` | `pixi install -e modeling` | shap, scikit-learn, scipy, matplotlib — for `intruder.modeling`; solved outside the pipeline group, and not on osx-64 (see pyproject.toml) |
| `analysis` | `pixi install -e analysis` | matplotlib, seaborn, scikit-learn, umap-learn — for `intruder.analysis` |
| `notebooks` | `pixi install -e notebooks` | JupyterLab and ipykernel on top of `analysis` |
| `dx` | `pixi install -e dx` | dxpy, for `scripts/dnanexus/` — see [DNAnexus docs](../../docs/scripts/DNANexus.md) |
| `backend` | `pixi install -e backend` | the web service: FastAPI, LangGraph, DuckDB, torch — shares nothing with the pipeline |

Put a dependency in `[project.dependencies]` only if *every* pipeline step
imports it. If one step needs something heavy, give it an extra under
`[project.optional-dependencies]` and a matching feature under `[tool.pixi]`, so
an image can install that step alone.

The compiled packages come from conda rather than PyPI, so no install builds
htslib from source.

## Common tasks

```bash
pixi add pandas                       # a runtime dependency, from conda
pixi add --pypi some-pypi-only-pkg    # one that only exists on PyPI
pixi add --feature analysis seaborn   # into another environment instead
pixi run lint-pipeline                # lint
pixi run test-pipeline                # run tests

# identify repeats using pyTRF from a multisample SV file
pixi run trf find --format vcf multisample.vcf trf_output.tsv

# annotate TRF output with novelty verdicts
pixi run novelty -i trf_output.tsv -o trf_novelty.tsv

# filter novelty output by motif purity and repeat coverage, etc.
pixi run filter -i trf_novelty.tsv -o trf_novelty_filtered.tsv

# annotate a VCF with per-ALT compressibility (SVCOMP) -- see
# ../../docs/scripts/annotate_compression.md
pixi run compression -i multisample.vcf -o multisample_comp.vcf
```

The command names above are unchanged by the move to `intruder/` — only the
module paths behind them shifted. To run a step without the console script, use
its module path: `pixi run python -m intruder.pipeline.novelty --help`.

## Tests

Tests live in `../../tests/python/`, mirroring the package layout — never beside
the code. `pixi run test-pipeline` from anywhere in the repo runs them all; CI runs
the same task on every push and pull request.

The Python version and every dependency are pinned in `../../pixi.lock`. Commit
it along with `pyproject.toml`.
