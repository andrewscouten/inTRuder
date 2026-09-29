"""Primitives shared by every tandem-repeat pipeline step.

Steps in this repository talk to each other through files, not imports, and
that stays true: nothing here knows about a catalogue or a caller. What lives
at this level is the small set of definitions that were provably identical
across steps and that are wrong in the same way when they are wrong -- how a
motif is reduced to a comparable key, how far off two motifs may be and still
be the same repeat, how a coordinate is converted and measured, how a table
column that is not a plain number is read back, and where a downloaded
catalogue is cached.

Duplicating those is not a harmless copy. Two steps that disagree by one base
on what "overlapping" means, or that fold strands differently, produce tables
that look comparable and are not.

Reading sequence off disk lives in ``trcore.io``: FASTA through pyfastx, VCF
insertions through cyvcf2, both yielding one record type that carries the
sequence and the coordinate origin its offsets are measured against. Two steps
need that now, the insertion caller and the genome caller, and the origin is
what separates them: a genome offset is a position on a contig, while an
insertion offset is a position inside sequence the reference does not have. A
second parser written against the same formats would diverge on exactly that,
the way a second coordinate convention would.

pyfastx and cyvcf2 are optional extras, imported inside ``trcore.io`` and
nowhere above it, so ``import trcore`` still needs nothing outside the standard
library and a step that only compares motifs does not acquire htslib to do it.

Deliberately *not* here: anything catalogue-shaped. ``novelty.catalog`` is a
numpy columnar index over millions of reference repeats with an on-disk cache;
``strchive.catalog`` is 82 JSON records carrying disease semantics. They share
three field names and no behaviour, and merging them would force pure-stdlib
code to import pandas to look up a disease locus.

Everything outside ``trcore.io`` is pure standard library. Vectorised wrappers
for catalogue-scale data live with the step that needs them (see
``novelty.platforms``).
"""

from .coords import (
    coverage_fraction,
    interval_distance,
    normalize_chrom,
    to_external,
    to_internal,
    union_length,
)
from .flanks import Flank, Flanks, flank_identity, flanks_of
from .motifs import (
    DEFAULT_EQUIVALENCE,
    DEFAULT_TOLERANCE,
    MATCH_EXACT,
    MATCH_FUZZY,
    MATCH_KINDS,
    MATCH_NONE,
    MATCH_SUBREPEAT,
    MATCH_VNTR,
    MAX_FUZZY_MOTIF,
    STR_MAX_MOTIF,
    MotifEquivalence,
    MotifMatch,
    MotifTolerance,
    canonical_motif,
    edit_budget,
    gc_fraction,
    least_rotation,
    motif_distance,
    primitive_unit,
    reverse_complement,
    shannon_entropy,
    tiling_distance,
)
from .utils.fetch import cache_root, download_bytes, download_file
from .utils.parse import as_number, as_size, as_total

__all__ = [
    "DEFAULT_EQUIVALENCE",
    "DEFAULT_TOLERANCE",
    "MATCH_EXACT",
    "MATCH_FUZZY",
    "MATCH_KINDS",
    "MATCH_NONE",
    "MATCH_SUBREPEAT",
    "MATCH_VNTR",
    "MAX_FUZZY_MOTIF",
    "STR_MAX_MOTIF",
    "Flank",
    "Flanks",
    "MotifEquivalence",
    "MotifMatch",
    "MotifTolerance",
    "as_number",
    "as_size",
    "as_total",
    "cache_root",
    "canonical_motif",
    "coverage_fraction",
    "download_bytes",
    "download_file",
    "edit_budget",
    "flank_identity",
    "flanks_of",
    "gc_fraction",
    "interval_distance",
    "least_rotation",
    "motif_distance",
    "normalize_chrom",
    "primitive_unit",
    "reverse_complement",
    "shannon_entropy",
    "tiling_distance",
    "to_external",
    "to_internal",
    "union_length",
]
