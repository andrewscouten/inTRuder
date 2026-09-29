"""The package imports without its parsers, which is what makes them optional.

``trcore.io`` is imported by steps that read one format, or neither. If a
parser were imported at module level, adding this package to ``trcore`` would
quietly make htslib a requirement of comparing two motifs.
"""

from __future__ import annotations

import builtins
import importlib
import sys

import pytest

PARSERS = ("pyfastx", "cyvcf2")


@pytest.fixture
def without_parsers(monkeypatch):
    """Make importing pyfastx and cyvcf2 fail, for this test only."""
    real_import = builtins.__import__

    def blocked(name, *args, **kwargs):
        if name.split(".")[0] in PARSERS:
            raise ImportError(f"blocked for test: {name}")
        return real_import(name, *args, **kwargs)

    for module in [m for m in sys.modules if m.split(".")[0] in PARSERS]:
        monkeypatch.delitem(sys.modules, module)
    for module in [m for m in sys.modules if m.startswith("intruder.trcore.io")]:
        monkeypatch.delitem(sys.modules, module)
    monkeypatch.setattr(builtins, "__import__", blocked)


def test_package_imports_with_no_parser_installed(without_parsers):
    io = importlib.import_module("intruder.trcore.io")
    assert io.Genomic("chr1").locate(1, 4) == ("chr1", 1, 4)


@pytest.mark.parametrize("parser, feature, call", [
    ("pyfastx", "fasta", lambda io: next(iter(io.read_fasta("genome.fna")))),
    ("cyvcf2", "sv", lambda io: next(iter(io.read_insertions("calls.vcf")))),
    ("cyvcf2", "sv", lambda io: io.sample_names("calls.vcf")),
])
def test_a_missing_parser_says_how_to_install_it(without_parsers, parser, feature, call):
    """Naming the package is not enough -- the error has to be actionable.

    pixi is this repository's environment manager, so the first thing the
    message offers is a pixi command. Saying only "pyfastx is missing" leaves
    the reader to work out which of four environments to ask for.
    """
    io = importlib.import_module("intruder.trcore.io")
    with pytest.raises(ImportError) as raised:
        call(io)

    message = str(raised.value)
    assert parser in message
    assert "pixi run -e pipeline" in message
    assert f"pixi add --feature {feature} {parser}" in message
