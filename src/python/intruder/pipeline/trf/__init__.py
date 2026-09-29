"""Finding tandem repeats in sequence, and filtering the calls.

One step of a file-in/file-out pipeline, and the first one: everything
downstream reads the table this writes. It imports no other pipeline step --
only :mod:`trcore`, for the record types, coordinate conventions and motif
primitives every step must agree on.

The finder does not know where its sequence came from. ``trcore.io`` normalises
an SV VCF's insertions and a reference FASTA's contigs to the same
:class:`trcore.io.SeqRecord`, and :func:`find.find_repeats` scans whichever it
is handed. What differs is on either side of that: which reader runs, which
parameter profile suits the input, and -- because an insertion has no reference
interval and a contig has no sample -- which schema the calls are written in.

    params      finder settings, and the sv/genome profiles
    find        the finder itself; records in, RepeatCall out
    filters     named thresholds on a call table, and the presets grouping them
    writers     the two output schemas, plus BED export
    cli         the `python -m intruder.pipeline.trf` command line
"""

from .filters import CATALOG, PASS, PRESETS, Check, add_repeat_coverage
from .find import RepeatCall, find_in_record, find_repeats
from .params import GENOME, PROFILES, SV, FinderParams
from .writers import (
    BED_COLUMNS,
    GENOME_COLUMNS,
    INSERTION_COLUMNS,
    genome_row,
    insertion_rows,
    write_bed,
    write_tsv,
)

__all__ = [
    "BED_COLUMNS",
    "CATALOG",
    "GENOME",
    "GENOME_COLUMNS",
    "INSERTION_COLUMNS",
    "PASS",
    "PRESETS",
    "PROFILES",
    "SV",
    "Check",
    "FinderParams",
    "RepeatCall",
    "add_repeat_coverage",
    "find_in_record",
    "find_repeats",
    "genome_row",
    "insertion_rows",
    "write_bed",
    "write_tsv",
]
