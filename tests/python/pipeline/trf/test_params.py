"""Finder settings, and the unit conversion that was silently wrong."""

from __future__ import annotations

import pytest

from intruder.pipeline.trf.params import DEFAULT_PROFILE, GENOME, PROFILES, SV, FinderParams

pytest.importorskip("pytrf")


def test_identity_is_converted_to_the_finders_percent_scale():
    """The regression this module exists for.

    ``ATRFinder`` measures identity in percent. ``sv_trfcaller.py`` passed its
    ``0.8`` straight through, so the finder was told "at least 0.8%" -- a gate
    that admits everything -- while the flag was documented as 80%.
    """
    assert FinderParams(min_identity=0.8).as_kwargs()["min_identity"] == 80.0
    assert FinderParams(min_identity=1.0).as_kwargs()["min_identity"] == 100.0


def test_seed_names_map_onto_the_finders_abbreviations():
    kwargs = FinderParams(min_seed_units=4, min_seed_length=14,
                          max_extension=25).as_kwargs()
    assert kwargs["min_seedrep"] == 4
    assert kwargs["min_seedlen"] == 14
    assert kwargs["max_extend"] == 25


def test_as_kwargs_passes_only_what_the_finder_accepts():
    """``quiet`` is ours, not the finder's; sending it through would raise."""
    import pytrf

    pytrf.ATRFinder("x", "ACGT" * 30, **FinderParams().as_kwargs())


def test_replace_ignores_unset_overrides():
    """The CLI layers every flag on, most of them ``None`` for "not given"."""
    params = SV.replace(min_motif=None, max_motif=7)
    assert params.min_motif == SV.min_motif
    assert params.max_motif == 7


def test_both_profiles_gate_on_identity():
    """Neither profile may ship the disabled gate the old defaults amounted to."""
    for name, profile in PROFILES.items():
        assert profile.min_identity >= 0.5, f"{name} would admit near-random sequence"


def test_every_source_has_a_default_profile():
    assert set(DEFAULT_PROFILE.values()) <= set(PROFILES)


def test_effective_min_units_is_whichever_floor_asks_for_more():
    """The two seed floors are one constraint, and the unit one is often inert.

    A seed must clear ``min_seed_units`` and ``min_seed_length`` together, so a
    trinucleotide under the GENOME profile needs four units rather than the three
    the profile names. Getting this wrong reads a grid point as a change when it
    is not one.
    """
    params = FinderParams(min_seed_units=3, min_seed_length=12)
    assert params.effective_min_units(3) == 4      # 3 units would be a 9 bp seed
    assert params.effective_min_units(4) == 3      # 12 bp exactly; the unit floor binds
    assert params.effective_min_units(1) == 12     # mononucleotide: all length
    assert params.effective_min_units(100) == 3    # long motif: all units

    assert params.seed_units_bind_from == 4
    assert FinderParams(min_seed_units=3, min_seed_length=9).seed_units_bind_from == 3
    assert FinderParams(min_seed_units=3, min_seed_length=6).seed_units_bind_from == 2


def test_effective_min_units_rejects_a_motif_shorter_than_a_base():
    with pytest.raises(ValueError, match="at least 1"):
        FinderParams().effective_min_units(0)


def test_effective_min_units_agrees_with_the_finder():
    """The arithmetic above, checked against ``ATRFinder`` rather than asserted.

    A CAG tract of exactly three units is a 9 bp seed: returned while the length
    floor is 9 or below, dropped at 10 and 12, with ``min_seed_units`` fixed at 3
    throughout. This is the step that makes a seed-length sweep move in jumps.
    """
    import pytrf

    flank = "GTCATGCCTAGGACTTGACCATGCAGTTACGGATCCTAAGCTTGGCCAATTGCA"
    for units in (3, 4):
        seq = flank + "CAG" * units + flank
        for seed_length in (6, 8, 9, 10, 12):
            params = FinderParams(min_seed_units=3, min_seed_length=seed_length,
                                  max_extension=0)
            found = bool(list(pytrf.ATRFinder("c", seq, **params.as_kwargs())))
            assert found == (units >= params.effective_min_units(3)), (
                f"{units} units of CAG at min_seed_length={seed_length}"
            )


def test_genome_profiles_unit_floor_is_inert_for_short_motifs():
    """Documenting the profile's own trap, so a change to it is deliberate.

    GENOME asks for 3 seed units but sets a 12 bp floor, which no mono-, di- or
    trinucleotide repeat can satisfy at 3 units. Those are the repeats a bacterial
    genome is mostly made of, which is why the profile stays flagged provisional.
    """
    assert GENOME.seed_units_bind_from == 4
    for motif_length in (1, 2, 3):
        assert GENOME.effective_min_units(motif_length) > GENOME.min_seed_units


def test_genome_profile_extension_is_survivable():
    """A large extension plus a working gate finds nothing at all.

    The inherited 10000 is not merely untuned: an extension that long always
    decays below the identity floor before it stops, so the scan returns empty.
    Whatever this is tuned to later, it has to stay small enough to survive the
    gate -- which is the coupling that makes these two a single choice.
    """
    assert GENOME.max_extension <= 100
