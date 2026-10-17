"""Tests for command-line argument validation."""

import sys

import pytest

from pdb2pqr import main


@pytest.mark.parametrize("entrypoint", ["main", "run_pdb2pqr"])
def test_retired_pkaani_fails_before_processing(
    entrypoint, tmp_path, monkeypatch, capsys
):
    """Retired backend requests must not touch inputs, outputs, or PROPKA."""

    def unexpected_processing(*args, **kwargs):
        pytest.fail("A retired backend request reached structure processing")

    monkeypatch.setattr(main, "main_driver", unexpected_processing)
    monkeypatch.setattr(main.io, "setup_logger", unexpected_processing)
    monkeypatch.setattr(main, "run_propka", unexpected_processing)
    args = [
        "--titration-state-method=pkaani",
        str(tmp_path / "missing.pdb"),
        str(tmp_path / "output.pqr"),
    ]
    with pytest.raises(SystemExit) as exc:
        if entrypoint == "main":
            monkeypatch.setattr(sys, "argv", ["pdb2pqr", *args])
            main.main()
        else:
            main.run_pdb2pqr(args)

    assert exc.value.code == 2
    error = capsys.readouterr().err
    assert "pKa-ANI support has been retired" in error
    assert "--titration-state-method=propka" in error
    assert "https://pdb2pqr.readthedocs.io/en/latest/releases.html" in error
    assert list(tmp_path.iterdir()) == []


def test_propka_selection_is_unchanged():
    """PROPKA remains opt-in with the requested pH."""
    parser = main.build_main_parser()
    args = parser.parse_args(
        [
            "--titration-state-method=propka",
            "--with-ph=6.5",
            "input.pdb",
            "output.pqr",
        ]
    )
    assert args.pka_method == "propka"
    assert args.ph == 6.5
    assert parser.parse_args(["input.pdb", "output.pqr"]).pka_method is None
