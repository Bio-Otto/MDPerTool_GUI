import argparse
from pathlib import Path
import pytest

from mdpertool.no_gui.cli_window import _resolve_output_directory, add_arguments_tu_subparsers


def test_resolve_output_directory_returns_cwd_when_none(monkeypatch, tmp_path):
    monkeypatch.chdir(tmp_path)

    resolved = _resolve_output_directory(None)

    assert resolved == str(tmp_path)


def test_resolve_output_directory_creates_target_directory(tmp_path):
    target = tmp_path / "nested" / "run_output"

    resolved = _resolve_output_directory(str(target))

    assert resolved == str(target.resolve())
    assert Path(resolved).exists()
    assert Path(resolved).is_dir()


def test_cli_parser_default_output_and_flags():
    parser = argparse.ArgumentParser()
    add_arguments_tu_subparsers(parser)

    parsed = parser.parse_args([
        "-p", "protein.pdb",
        "-pert_res", "ASN17A",
        "-speed_factor", "4",
    ])

    assert parsed.output is None
    assert parsed.write_group_summary is True
    assert parsed.write_run_manifest is True
    assert parsed.use_switching_distance is True


def test_cli_parser_rejects_invalid_boolean_value():
    parser = argparse.ArgumentParser()
    add_arguments_tu_subparsers(parser)

    with pytest.raises(SystemExit):
        parser.parse_args([
            "-p", "protein.pdb",
            "-pert_res", "ASN17A",
            "-speed_factor", "4",
            "-swd-use", "not_bool",
        ])


def test_cli_parser_requires_perturbation_arguments():
    parser = argparse.ArgumentParser()
    add_arguments_tu_subparsers(parser)

    with pytest.raises(SystemExit):
        parser.parse_args([
            "-p", "protein.pdb",
        ])
