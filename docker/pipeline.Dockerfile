# Minimal environment for the novel TR pipeline's Python-based processes.
# Extend this as new steps (samtools, minimap2, etc.) get added -- a tool with a
# bioconda package is one line in pyproject.toml rather than a second package
# manager in here.

# Pinned for reproducible builds. The image ships pixi and nothing else of note.
FROM ghcr.io/prefix-dev/pixi:0.78.0-bookworm-slim

# procps provides `ps`, which Nextflow needs from inside the container to
# collect task resource metrics. Everything else the pipeline used to compile
# against -- htslib for cyvcf2, the autotools chain for parasail -- now arrives
# prebuilt from conda, so no build toolchain is installed here at all.
RUN apt-get update && apt-get install -y --no-install-recommends \
    procps \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# pixi.lock pins every package for linux-64 exactly as it is pinned on a
# contributor's machine; --locked refuses to resolve anything else. pyproject
# and src/python come along because `intruder` is installed from this directory,
# not downloaded -- the build is building this project too, not just its
# dependencies. README.md is there because pyproject.toml declares it.
COPY pyproject.toml pixi.lock README.md ./
COPY src/python ./src/python

RUN pixi install --locked --environment default \
    && pixi clean cache --yes

# The environment on PATH instead of an activation hook: Nextflow runs each
# process body with its own shell and never sees an ENTRYPOINT, so a container
# that needs activating is a container whose commands are not found.
ENV PATH=/app/.pixi/envs/default/bin:$PATH

# Bundle sv_trfcaller.py at a stable, simple path for FIND_TRS to call
# (in addition to it already being installed as part of the novelty
# package above).
COPY src/python/intruder/pipeline/trf/sv_trfcaller.py /opt/scripts/sv_trfcaller.py

# Bundle normalize_svtype.py (from the annotation team's sv_preprocess
# pipeline - not yet merged to their main branch, so this is a local
# copy for now) for the PREPROCESS process to call directly. It's pure
# standard-library Python (argparse + gzip only), so no extra
# dependencies are needed for it.
COPY pipelines/sv_preprocess/scripts/normalize_svtype.py /opt/scripts/normalize_svtype.py

# Bundle filter_ins_trf.py (from the teammate who built the 02B
# filtering stage) for the FILTER_BY_COVERAGE process to call directly.
# It's pure standard-library Python (no external dependencies), so no
# extra packages are needed for it.
COPY src/python/intruder/pipeline/trf/filter_ins_trf.py /opt/scripts/filter_ins_trf.py

# --- Pre-bake novelty's reference catalogs (UCSC simpleRepeat ~30MB,
# TRExplorer ~45MB) at build time, so runs never need network access
# and are instant.
#
# IMPORTANT: novelty's default cache dir ("data/reference") is a
# RELATIVE path. Nextflow runs each task in its own work directory
# (not /app), so a relative path baked in at /app/data/reference would
# NOT be found at runtime - the container would silently look in the
# wrong place and re-download every run, defeating the whole point.
# Setting NOVELTY_CACHE as an absolute path fixes the location
# regardless of what directory a command is actually run from.
ENV NOVELTY_CACHE=/opt/novelty_cache

# Dummy query forces both catalogs to download and cache now, at build
# time, rather than on first real use.
RUN mkdir -p /opt/novelty_cache && \
    novelty --platform ucsc,trexplorer query --chrom chr1 --pos 1000000 --motif AT
# More dependencies go in pyproject.toml (or pixi.toml, for a conda one) as
# later pipeline stages need them - then re-run `pixi lock` locally, commit the
# updated pixi.lock, and rebuild.
