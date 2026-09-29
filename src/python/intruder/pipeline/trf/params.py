"""Finder settings, and the profiles that make them fit their input.

``pytrf.ATRFinder`` takes six numbers and has no opinion about what it is
scanning. That is the right design for a finder and the wrong one for a
pipeline: the same six numbers mean different things against a 10 kb insertion
and against a 5 Mb bacterial chromosome, and the failure is silent -- you get
calls either way, just far too many or far too few.

So the numbers are named once here and selected as a *profile*, rather than
being six defaults on an argument parser that every caller silently inherits.

    sv        VCF insertions, corrected -- see the note above :data:`SV`
    genome    provisional, for whole-contig FASTA scans

None of these are settled values. They are hyperparameters with no ground
truth to optimise against, which is the same position ``novelty`` is in with
its screen thresholds: the useful move there was to enumerate a grid and see
which conclusions survive the choice, rather than to pick an optimum. Doing
that here needs real genomes rather than the synthetic checks these defaults
came from, so the values stand until then and the CLI exposes every one of them.
"""

from __future__ import annotations

from dataclasses import dataclass, replace


@dataclass(frozen=True)
class FinderParams:
    """One setting per ``ATRFinder`` knob, under this repository's names.

    The finder's keyword names are abbreviations of its seed-and-extend
    algorithm (``min_seedrep``, ``min_seedlen``, ``max_extend``); the names here
    spell the same concepts out. :meth:`as_kwargs` is the one place the two
    vocabularies meet, and the only place a unit is converted.
    """

    min_motif: int = 1
    max_motif: int = 500

    #: Minimum purity, 0..1 -- the same scale as the ``purity`` output column.
    #: ``ATRFinder`` wants this as a percentage; :meth:`as_kwargs` converts.
    min_identity: float = 0.8

    # The finder seeds on an exact repeat, then extends through mismatches
    # while identity holds up. These three describe *that*, not the repeat that
    # comes out: `sv_trfcaller.py` called them min_rep_units / min_rep_length /
    # max_rep_length, which reads as bounds on the call and collides with the
    # `rep_units` and `rep_length` output columns, which are different
    # quantities. `max_extension` in particular caps nothing: with it set to 0
    # a 12 bp seed still returns a 90 bp call.
    min_seed_units: int = 3
    min_seed_length: int = 10
    max_extension: int = 10_000

    #: Send the finder's own stderr chatter to devnull while it scans.
    quiet: bool = True

    def as_kwargs(self) -> dict[str, object]:
        """The keyword arguments ``pytrf.ATRFinder`` actually takes.

        ``min_identity`` crosses a scale boundary here. The finder measures
        identity as a percentage (``repeat.identity`` comes back as ``69.0``),
        while every purity in this repository is a fraction. ``sv_trfcaller.py``
        passed its ``0.8`` straight through, so the finder was told to accept
        anything above 0.8 *percent* and the gate never fired -- calls at 55%
        identity came back from a run nominally filtered at 80%. Converting in
        the one place the two vocabularies meet is what stops that recurring.
        """
        return {
            "min_motif": self.min_motif,
            "max_motif": self.max_motif,
            "min_identity": self.min_identity * 100,
            "min_seedrep": self.min_seed_units,
            "min_seedlen": self.min_seed_length,
            "max_extend": self.max_extension,
        }

    def replace(self, **overrides) -> FinderParams:
        """A copy with ``overrides`` applied, ignoring the ones that are ``None``.

        Lets the CLI layer per-flag overrides onto a profile without having to
        know which flags the user actually passed.
        """
        return replace(self, **{k: v for k, v in overrides.items() if v is not None})


# Settings for VCF insertion sequence. NOT the defaults `sv_trfcaller.py`
# shipped: those passed `min_identity=0.8` to a finder measuring identity in
# percent, so the gate was set to 0.8% and never fired, and every call it made
# reached the table however impure. Correcting the scale alone returns *zero*
# calls at `max_extension=10000` -- an unbounded extension always decays below
# 80% before it stops -- so the two move together or not at all.
#
# Consequence, deliberately accepted: stage 01 tables produced before this are
# not reproducible with this profile. They contain calls down to ~55% identity.
SV = FinderParams(
    min_identity=0.8,
    max_extension=20,
)

# Provisional, and deliberately left that way. `max_motif=500` is sized for a
# long-read insertion of at most a few kb; against a 5 Mb contig it invites the
# finder to extend a seed across a large fraction of a gene, so the ceiling
# comes down and the seed floor comes up.
#
# `max_extension` is the knob that matters, and the inherited 10000 is not
# merely untuned but unusable: with the identity gate working, an extension
# that large always decays below 80% before it stops, and a scan of a genome
# with three planted tracts returns nothing at all. So it takes the same value
# as SV and stays flagged.
#
# What that value is evidenced for is narrow, and worth stating exactly. On a
# synthetic scan it recovers *that* a tract is there; it does not recover the
# tract. Three planted 120 bp tracts come back as 160/155/159 bp -- inflated by
# the full extension at each end -- with purity dragged from 1.00 to ~0.85 by
# the absorbed background. The value that is exactly right on that test is 0,
# and the only firm result is the upper bound: anything from about 100 up finds
# nothing. It trades boundary accuracy against tolerance of a decayed tract, and
# which way to trade is a property of real bacterial repeats, not of a synthetic
# test where every planted tract is perfect by construction: at 0 a decayed CAG
# tract is truncated at its first substitution, at 50 the reported motif itself
# destabilises (CAG collapses into a 15-mer). Tuning it belongs against the
# synced contaminant and HMP reference genomes -- see the parameter-profile item
# in docs/todo.md in the parent project -- so it stays small, stays provisional,
# and the CLI exposes it.
GENOME = FinderParams(
    max_motif=100,
    min_seed_length=12,
    max_extension=20,      # provisional; see above
)

PROFILES = {"sv": SV, "genome": GENOME}

#: The profile to use when the reader was picked by sniffing rather than named.
DEFAULT_PROFILE = {"vcf": "sv", "fasta": "genome"}
