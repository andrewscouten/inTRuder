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


def test_genome_profile_extension_is_survivable():
    """A large extension plus a working gate finds nothing at all.

    The inherited 10000 is not merely untuned: an extension that long always
    decays below the identity floor before it stops, so the scan returns empty.
    Whatever this is tuned to later, it has to stay small enough to survive the
    gate -- which is the coupling that makes these two a single choice.
    """
    assert GENOME.max_extension <= 100
