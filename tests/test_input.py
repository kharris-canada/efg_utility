from __future__ import annotations

from pathlib import Path

import pytest

from efg_utility.constants import BOHR_RADIUS_M
from efg_utility.input import (
    ANGSTROM_IN_METRES,
    InputFileError,
    read_hydrogenic_orbitals,
    read_point_charges,
)


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


def test_reads_default_explicit_and_zeff_orbital_scales(tmp_path: Path) -> None:
    input_file = tmp_path / "orbitals.txt"
    input_file.write_text(
        "H1 0 0 1 1s - - 1\n"
        "H2 1 0 0 1s 0.4 - 0.5\n"
        "H3 0 1 0 1s - 2.0 2\n",
        encoding="utf-8",
    )

    orbitals = read_hydrogenic_orbitals(input_file)

    assert len(orbitals) == 3
    assert orbitals[0].radius_m == pytest.approx(BOHR_RADIUS_M)
    assert orbitals[1].radius_m == pytest.approx(0.4 * ANGSTROM_IN_METRES)
    assert orbitals[1].occupation == 0.5
    assert orbitals[2].radius_m == pytest.approx(BOHR_RADIUS_M / 2.0)
    assert orbitals[2].occupation == 2.0


def test_reads_2s_and_real_2p_orbitals(tmp_path: Path) -> None:
    input_file = tmp_path / "orbitals.txt"
    input_file.write_text(
        "A 0 0 1 2s - - 1\n"
        "B 0 0 1 2px - - 1\n"
        "C 0 0 1 2p_y - - 1\n"
        "D 0 0 1 2PZ - - 1\n",
        encoding="utf-8",
    )

    orbitals = read_hydrogenic_orbitals(input_file)

    assert [orbital.orbital_type for orbital in orbitals] == ["2s", "2p", "2p", "2p"]
    assert [orbital.axis for orbital in orbitals[1:]] == [
        (1.0, 0.0, 0.0),
        (0.0, 1.0, 0.0),
        (0.0, 0.0, 1.0),
    ]


def test_reads_and_normalizes_arbitrary_2p_axis(tmp_path: Path) -> None:
    input_file = tmp_path / "orbital.txt"
    input_file.write_text("P 0 0 1 2p(1,1,0) - - 1\n", encoding="utf-8")

    orbital = read_hydrogenic_orbitals(input_file)[0]

    assert orbital.orbital_type == "2p"
    assert orbital.axis == pytest.approx((2.0**-0.5, 2.0**-0.5, 0.0))


@pytest.mark.parametrize(
    ("row", "message"),
    [
        ("H 0 0 1 1s - -", "expected 8 fields"),
        ("H 0 0 1 2p - - 1", "requires an axis"),
        ("H 0 0 1 2p(0,0,0) - - 1", "axis must be nonzero"),
        ("H 0 0 1 1s 0.5 2 1", "mutually exclusive"),
        ("H 0 0 1 1s 0 - 1", "radius must be positive"),
        ("H 0 0 1 1s - 0 1", "effective nuclear charge must"),
        ("H 0 0 1 1s - - -1", "occupation must be nonnegative"),
    ],
)
def test_reports_invalid_orbital_rows(
    tmp_path: Path, row: str, message: str
) -> None:
    input_file = tmp_path / "bad_orbital.txt"
    input_file.write_text(f"# comment\n{row}\n", encoding="utf-8")

    with pytest.raises(InputFileError) as error:
        read_hydrogenic_orbitals(input_file)

    assert f"{input_file}:2:" in str(error.value)
    assert message in str(error.value)
