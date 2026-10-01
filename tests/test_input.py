from __future__ import annotations

from pathlib import Path

import pytest

from efg_utility.input import ANGSTROM_IN_METRES, InputFileError, read_point_charges


def test_reads_comments_scientific_notation_and_signed_integer(tmp_path: Path) -> None:
    input_file = tmp_path / "charges.txt"
    input_file.write_text(
        "# name x y z charge\n\nNa 1e0 -2.5 3.0 +1 # nearby ion\n",
        encoding="utf-8",
    )

    charges = read_point_charges(input_file)

    assert len(charges) == 1
    assert charges[0].position_m == pytest.approx(
        (1.0 * ANGSTROM_IN_METRES, -2.5 * ANGSTROM_IN_METRES, 3.0 * ANGSTROM_IN_METRES)
    )
    assert charges[0].charge_e == 1


@pytest.mark.parametrize(
    ("row", "message"),
    [
        ("Na 1 2 +1", "expected 5 fields"),
        ("Na nope 2 3 +1", "invalid x coordinate"),
        ("Na nan 2 3 +1", "x coordinate must be finite"),
        ("Na 1 2 3 0.5", "charge must be an integer"),
        ("Na 1 2 3 3s", "charge must be an integer"),
    ],
)
def test_reports_filename_and_line_for_invalid_rows(
    tmp_path: Path, row: str, message: str
) -> None:
    input_file = tmp_path / "bad.txt"
    input_file.write_text(f"# comment\n{row}\n", encoding="utf-8")

    with pytest.raises(InputFileError) as error:
        read_point_charges(input_file)

    assert f"{input_file}:2:" in str(error.value)
    assert message in str(error.value)


def test_reports_unreadable_file_without_traceback_details(tmp_path: Path) -> None:
    missing = tmp_path / "missing.txt"

    with pytest.raises(InputFileError, match="missing.txt"):
        read_point_charges(missing)
