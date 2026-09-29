"""Reading reference sequence out of a FASTA, through pyfastx.

FASTA has no specification. Line width, soft-masking case, ambiguity codes and
line endings all vary between sources, and the header after ``>`` is free text
that NCBI, Ensembl, prokka and the HMP reference genomes each lay out their own
way. pyfastx absorbs the file-level variation; the header is carried through
unparsed on :attr:`SeqRecord.description` because guessing at its structure is
how a caller ends up with an organism name that is right for one reference and
silently wrong for the next.

Sequence is upper-cased on the way out. Soft-masked lower-case bases are common
in reference assemblies, and every motif comparison in :mod:`trcore.motifs` is
case-sensitive, so leaving them alone makes ``acg`` and ``ACG`` different
repeats.

Nothing is written next to the input. pyfastx can build a ``.fxi`` index beside
a FASTA for random access, which is worth having for region fetches but drops a
new file into whatever directory the genome lives in -- and in this project
those directories are checksummed against a manifest lock. So the default
streams, and indexing is opt-in through :func:`open_indexed` for callers that
want ``fetch`` and can afford the sidecar.
"""

from __future__ import annotations

from collections.abc import Iterator
from pathlib import Path

from .records import Genomic, SeqRecord


def _split_header(header: str) -> tuple[str, str]:
    """``>`` line -> ``(name, description)``, split on the first run of whitespace.

    Whitespace, not a space: some references separate the contig name from its
    description with a tab, and splitting on ``" "`` alone folds the whole line
    into the name -- which then becomes the ``contig`` value in every row and
    every BED interval the scan writes.
    """
    parts = header.split(None, 1)
    if not parts:
        return "", ""
    name, description = parts[0], parts[1] if len(parts) > 1 else ""
    return name, description.strip()


def read_fasta(path: str | Path, *, uppercase: bool = True) -> Iterator[SeqRecord]:
    """Stream every contig of ``path`` as a :class:`SeqRecord` with a genomic origin.

    Handles plain and gzip-compressed FASTA. Reads one record at a time and
    builds no index, so a 5 Mb bacterial chromosome costs one contig of memory
    and a 198-genome catalogue run touches nothing on disk but the inputs.
    """
    pyfastx = _import_pyfastx()
    for header, sequence in pyfastx.Fasta(str(path), build_index=False, full_name=True):
        name, description = _split_header(header)
        yield SeqRecord(
            name=name,
            sequence=sequence.upper() if uppercase else sequence,
            origin=Genomic(contig=name),
            description=description,
        )


def open_indexed(path: str | Path, *, uppercase: bool = True):
    """A pyfastx ``Fasta`` with a ``.fxi`` index, for random access by contig or region.

    Writes ``<path>.fxi`` next to the FASTA on first call. Use it when the
    access pattern is "fetch these regions", and prefer :func:`read_fasta` for a
    whole-genome scan, which needs no index at all.
    """
    pyfastx = _import_pyfastx()
    return pyfastx.Fasta(str(path), uppercase=uppercase)


def fetch_region(fasta, contig: str, start: int, end: int) -> SeqRecord:
    """A 0-based half-open window of ``fasta``, as a record that knows where it came from.

    ``fasta`` is what :func:`open_indexed` returned. The window's
    :class:`Genomic` origin carries ``start``, so a repeat found at offset 40 of
    a window beginning at 1_000 reports position 1_040 rather than 40.
    """
    sequence = fasta.fetch(contig, (start + 1, end))   # pyfastx regions are 1-based inclusive
    return SeqRecord(
        name=f"{contig}:{start}-{end}",
        sequence=sequence,
        origin=Genomic(contig=contig, start=start),
    )


def _import_pyfastx():
    """Import pyfastx, or explain which extra supplies it."""
    try:
        import pyfastx
    except ImportError as exc:                              # pragma: no cover - import guard
        raise ImportError(
            "reading FASTA needs pyfastx, which is not installed by default.\n"
            "  in this repo:  pixi run -e pipeline <your command>\n"
            "  to add it:     pixi add --feature fasta pyfastx\n"
            "  installed:     pip install 'intruder[fasta]'"
        ) from exc
    return pyfastx
