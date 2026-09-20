# Getting Started

Everything needed to set up inTRuder, run the web interface, and run the pipeline — both through
Nextflow and as standalone commands.

## Prerequisites

- [pixi](https://pixi.sh) — manages every environment in the repo, Python and Node included
- Nothing else: Node comes from the `frontend` environment, and every command below is a pixi
  task. Run `pixi task list` to see them all.
- [Docker](https://docs.docker.com/get-docker/) — only needed for the containerized setup or for
  running the Nextflow pipeline

## 1. Clone and set up

```bash
git clone https://github.com/collaborativebioinformatics/inTRuder.git
cd inTRuder
pixi run setup     # every environment, plus npm install for the web app and the
                   # synthetic demo dataset the interface opens on
```

`pixi run setup` also copies `backend/.env.example` to `backend/.env` — add a model credential there
to enable chat in the web interface. Skip this step entirely if you only want the containers.

## 2. Run the web interface

```bash
pixi run dev       # backend on :8000, frontend on :3000
```

Or in containers, with no toolchain to install:

```bash
docker compose up --build      # same two ports; your data/ is bind-mounted, not baked in
```

Either way, open [http://localhost:3000](http://localhost:3000). `pixi run backend` and
`pixi run frontend` run the two halves on their own. Backend internals — model providers, the agent's
tools, the SQL sandbox — are covered in [`backend/README.md`](../backend/README.md).

## 3. Run the pipeline

The [Nextflow](https://www.nextflow.io/) entrypoint is `workflows/main.nf`, and its processes run
inside the image built from [`docker/pipeline.Dockerfile`](../docker/pipeline.Dockerfile) — so build
that once first:

```bash
docker build -f docker/pipeline.Dockerfile -t novel-tr-pipeline:latest .
nextflow run workflows/main.nf -profile docker
```

With no `--input_vcf`, that runs stage 01 (TR detection) over the 500 committed insertions in
`data/sv_output/sniffles/first_500_INS.vcf`. Point it at your own SV calls and turn the optional
stages on with flags:

```bash
nextflow run workflows/main.nf -profile docker \
    --input_vcf path/to/insertions.vcf \
    --run_novelty \
    --run_annotation \
    --run_validation --tr_catalogue_bed path/to/tandem_repeats.bed
```

Results are published under `results/`, with `results/final/final_output.tsv` as the one
predictable endpoint and run reports in `reports/`. Which processes are real and which are still
placeholders is tracked in the
[Methods outline](Methods_overview.md#5-pipeline-orchestration).

Each step is also a standalone CLI, file in and file out, so nothing forces you through Nextflow:

```bash
pixi run svpytrf -i multisample.vcf -o trf.tsv                    # 01  TRs inside inserted alleles
pixi run novelty --platform ucsc,trexplorer annotate trf.tsv trf.novelty.tsv   # 02  known or novel?
```

See the [Novelty screen](tools/NOVELTY_SCREEN.md) and [STRchive comparison](tools/STRCHIVE_COMPARE.md)
docs for what each stage actually does, and the [Methods outline](Methods_overview.md) for the
full pipeline write-up.

## Other environments

- [Python source](../src/python/README.md) — pixi environments, adding dependencies, running
  scripts, linting and tests
- [R source](../src/R/README.md) — renv-managed environment, `renv::restore()`, snapshotting new
  packages
- [Run programs on DNAnexus](scripts/DNANexus.md) — `scripts/dnanexus/dx-*.sh`: start a machine,
  get a terminal or run one program on it, then stop the machine
