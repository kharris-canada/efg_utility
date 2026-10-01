from __future__ import annotations

import os
from pathlib import Path
import subprocess
import sys

from efg_utility.cli import main


def test_console_entry_handler_prints_result(tmp_path: Path, capsys: object) -> None:
    input_file = tmp_path / "charge.txt"
    input_file.write_text("X 0 0 1 +1\n", encoding="utf-8")

    exit_code = main(["calculate", str(input_file)])
    captured = capsys.readouterr()  # type: ignore[attr-defined]

    assert exit_code == 0
    assert captured.out.startswith("VZZ = +")
    assert " V/m^2\neta = 0\n" in captured.out
    assert captured.err == ""


def test_module_invocation(tmp_path: Path) -> None:
    input_file = tmp_path / "charge.txt"
    input_file.write_text("X 0 0 1 +1\n", encoding="utf-8")
    environment = os.environ.copy()
    source_root = str(Path(__file__).parents[1] / "src")
    environment["PYTHONPATH"] = os.pathsep.join(
        filter(None, (source_root, environment.get("PYTHONPATH")))
    )

    completed = subprocess.run(
        [sys.executable, "-m", "efg_utility", "calculate", str(input_file)],
        check=False,
        capture_output=True,
        text=True,
        env=environment,
    )

    assert completed.returncode == 0
    assert completed.stdout.startswith("VZZ = +")
    assert "eta = 0" in completed.stdout
    assert completed.stderr == ""


def test_input_error_has_no_traceback(tmp_path: Path, capsys: object) -> None:
    input_file = tmp_path / "bad.txt"
    input_file.write_text("X 0 0 1 fractional\n", encoding="utf-8")

    exit_code = main(["calculate", str(input_file)])
    captured = capsys.readouterr()  # type: ignore[attr-defined]

    assert exit_code == 2
    assert captured.out == ""
    assert "charge must be an integer" in captured.err
    assert "Traceback" not in captured.err


def test_empty_file_reports_undefined_eta(tmp_path: Path, capsys: object) -> None:
    input_file = tmp_path / "empty.txt"
    input_file.write_text("# no charges\n", encoding="utf-8")

    exit_code = main(["calculate", str(input_file)])
    captured = capsys.readouterr()  # type: ignore[attr-defined]

    assert exit_code == 0
    assert "VZZ = +0.000000000000e+00 V/m^2" in captured.out
    assert "eta = undefined" in captured.out
