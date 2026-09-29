"""Reading sequence off disk, normalised to one record type.

Two sources supply sequence to scan for repeats: a reference FASTA (pyfastx)
and the insertion calls in an SV VCF (cyvcf2). Both land as a
:class:`SeqRecord` carrying the sequence and an :class:`Origin` that says what
its offsets mean, so a repeat scan works the same either way and cannot quietly
treat an offset inside an insertion as a position on a contig.

The parsers are optional extras, imported when a reader is *called* rather than
when this package is imported. ``import intruder.trcore.io`` therefore costs
nothing but the standard library, and a step that only reads FASTA never needs
htslib on the machine. Calling a reader whose parser is missing raises an
``ImportError`` naming the extra that supplies it.

    from intruder.trcore.io import read_fasta

    for record in read_fasta("genome.fna.gz"):
        for repeat in pytrf.ATRFinder(record.name, record.sequence):
            contig, start, end = record.origin.locate(repeat.start - 1, repeat.end)
"""

from .fasta import fetch_region, open_indexed, read_fasta
from .records import Genomic, Inserted, Origin, SeqRecord
from .vcf import read_insertions, sample_names

__all__ = [
    "Genomic",
    "Inserted",
    "Origin",
    "SeqRecord",
    "fetch_region",
    "open_indexed",
    "read_fasta",
    "read_insertions",
    "sample_names",
]
