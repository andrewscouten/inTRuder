"""Infrastructure the steps share, none of which knows what a tandem repeat is.

``trcore`` holds two different kinds of thing, and they read very differently.
:mod:`trcore.coords`, :mod:`trcore.motifs` and :mod:`trcore.flanks` are the
*domain*: a coordinate, a motif, the sequence either side of a call. Change one
and you have changed what this project means by a repeat.

What lives here is the other kind -- plumbing that would look the same in a
project about anything else:

    parse     a table cell that is not a plain number -> the number it holds
    fetch     move bytes over the network, and pick a directory to cache them in

The admission rule, so this does not become the drawer everything ends up in:
a module belongs here only if it is shared by more than one step *and* its
docstring could be written without the words "repeat", "motif" or "locus". A
module that fails the second half is domain logic and belongs a level up, next
to ``coords.py``, however utility-shaped it looks. ``paths.py`` stays up there
for a different reason: it is the marker-walk that finds the checkout, and
:mod:`trcore.utils.fetch` depends on it rather than the other way round.

Everything here is pure standard library, like the rest of ``trcore`` outside
``trcore.io``.
"""

from .fetch import cache_root, download_bytes, download_file
from .parse import as_number, as_size, as_total

__all__ = [
    "as_number",
    "as_size",
    "as_total",
    "cache_root",
    "download_bytes",
    "download_file",
]
