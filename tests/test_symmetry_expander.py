from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

import expand_symmetry
from efg_utility import calculate_efg, calculate_orbital_field
from efg_utility.input import read_hydrogenic_orbitals, read_point_charges


@pytest.mark.parametrize(
    ("name", "order"),
    [
        ("C1", 1),
        ("Ci", 2),
        ("Cs", 2),
        ("C5", 5),
        ("C4v", 8),
        ("C3h", 6),
        ("S6", 6),
        ("D5", 10),
        ("D4h", 16),
        ("D3d", 12),
        ("T", 12),
        ("Th", 24),
        ("Td", 24),
        ("O", 24),
        ("O_h", 48),
        ("I", 60),
        ("Ih", 120),
    ],
)
def test_point_group_orders_and_matrix_properties(name: str, order: int) -> None:
    operations = expand_symmetry.point_group_operations(name)
    keys = {tuple(np.round(item, 8).ravel()) for item in operations}

    assert len(operations) == order
    assert len(keys) == order
    for operation in operations:
        assert operation.T @ operation == pytest.approx(np.identity(3), abs=2.0e-12)
        assert abs(float(np.linalg.det(operation))) == pytest.approx(1.0)
        assert tuple(np.round(operation.T, 8).ravel()) in keys
    for first in operations:
        for second in operations:
            assert tuple(np.round(first @ second, 8).ravel()) in keys


def _record(
    position: tuple[float, float, float],
    *,
    payload: tuple[str, ...] = ("-2",),
    axis: tuple[float, float, float] | None = None,
) -> expand_symmetry.SeedRecord:
    return expand_symmetry.SeedRecord(
        "X",
        np.asarray(position, dtype=np.float64),
        payload,
        np.asarray(axis, dtype=np.float64) if axis is not None else None,
    )


def test_oh_general_special_and_origin_orbit_sizes() -> None:
    operations = expand_symmetry.point_group_operations("Oh")

    assert len(
        expand_symmetry.expand_record(
            _record((0.321, 0.547, 0.918)),
            operations,
            position_tolerance=1.0e-9,
        )
    ) == 48
    assert len(
        expand_symmetry.expand_record(
            _record((0.0, 0.0, 1.0)),
            operations,
            position_tolerance=1.0e-9,
        )
    ) == 6
    assert len(
        expand_symmetry.expand_record(
            _record((0.0, 0.0, 0.0)),
            operations,
            position_tolerance=1.0e-9,
        )
    ) == 1


def test_oh_p_orbital_at_origin_generates_three_axes() -> None:
    operations = expand_symmetry.point_group_operations("Oh")
    record = _record(
        (0.0, 0.0, 0.0),
        payload=("2pz", "-", "-", "1"),
        axis=(0.0, 0.0, 1.0),
    )
    expanded = expand_symmetry.expand_record(
        record, operations, position_tolerance=1.0e-9
    )

    assert len(expanded) == 3
    assert {
        expand_symmetry.format_p_descriptor(axis)
        for _position, axis in expanded
        if axis is not None
    } == {"2px", "2py", "2pz"}


def test_td_p_axis_orbits_at_tetrahedral_vertex() -> None:
    operations = expand_symmetry.point_group_operations("Td")
    position = np.array((1.0, 1.0, 1.0))

    laboratory_z = expand_symmetry.SeedRecord(
        "X",
        position,
        ("2pz", "-", "6.0", "0.333"),
        np.array((0.0, 0.0, 1.0)),
    )
    transverse_orbit = expand_symmetry.expand_record(
        laboratory_z, operations, position_tolerance=1.0e-9
    )
    assert len(transverse_orbit) == 12
    by_position: dict[tuple[int, int, int], set[str]] = {}
    for center, axis in transverse_orbit:
        assert axis is not None
        key = tuple(round(value) for value in center)
        by_position.setdefault(key, set()).add(
            expand_symmetry.format_p_descriptor(axis)
        )
    assert len(by_position) == 4
    assert all(axes == {"2px", "2py", "2pz"} for axes in by_position.values())
    formatted = expand_symmetry.format_expanded_records(
        [laboratory_z], operations, position_tolerance=1.0e-9
    )
    output_positions = [tuple(line.split()[1:4]) for line in formatted.splitlines()]
    assert all(
        len(set(output_positions[index : index + 3])) == 1
        for index in (0, 3, 6, 9)
    )

    radial = expand_symmetry.SeedRecord(
        "X",
        position,
        ("2p(1,1,1)", "-", "6.0", "0.333"),
        position.copy(),
    )
    radial_orbit = expand_symmetry.expand_record(
        radial, operations, position_tolerance=1.0e-9
    )
    assert len(radial_orbit) == 4
    for center, axis in radial_orbit:
        assert axis is not None
        assert abs(float(np.dot(center, axis))) == pytest.approx(
            float(np.linalg.norm(center)), abs=1.0e-12
        )


def test_c3_generates_general_p_axis_descriptors() -> None:
    operations = expand_symmetry.point_group_operations("C3")
    record = _record(
        (0.0, 0.0, 1.0),
        payload=("2px", "-", "-", "1"),
        axis=(1.0, 0.0, 0.0),
    )

    text = expand_symmetry.format_expanded_records(
        [record], operations, position_tolerance=1.0e-9
    )

    assert len(text.splitlines()) == 3
    assert "2px" in text
    assert "2p(" in text


def test_generated_oh_charge_file_is_valid_and_has_zero_efg(tmp_path: Path) -> None:
    seed_file = tmp_path / "seed.txt"
    output_file = tmp_path / "expanded.txt"
    seed_file.write_text("X 0 0 1 -2\n", encoding="utf-8")
    records = expand_symmetry.read_seed_records(seed_file, "charge")
    text = expand_symmetry.format_expanded_records(
        records,
        expand_symmetry.point_group_operations("Oh"),
        position_tolerance=1.0e-9,
    )
    output_file.write_text(text, encoding="utf-8")

    charges = read_point_charges(output_file)
    result = calculate_efg(charges)

    assert len(charges) == 6
    assert {charge.position_m for charge in charges} == {
        (1.0e-10, 0.0, 0.0),
        (-1.0e-10, 0.0, 0.0),
        (0.0, 1.0e-10, 0.0),
        (0.0, -1.0e-10, 0.0),
        (0.0, 0.0, 1.0e-10),
        (0.0, 0.0, -1.0e-10),
    }
    assert result.vzz == 0.0
    assert result.eta is None


def test_generated_general_orbitals_round_trip_to_efg_parser(tmp_path: Path) -> None:
    seed_file = tmp_path / "seed.txt"
    output_file = tmp_path / "expanded.txt"
    seed_file.write_text("P 0 0 1 2p(1,1,0) - 2 0.5\n", encoding="utf-8")
    records = expand_symmetry.read_seed_records(seed_file, "orbital")
    output_file.write_text(
        expand_symmetry.format_expanded_records(
            records,
            expand_symmetry.point_group_operations("C3"),
            position_tolerance=1.0e-9,
        ),
        encoding="utf-8",
    )

    orbitals = read_hydrogenic_orbitals(output_file)

    assert len(orbitals) == 3
    assert all(orbital.axis is not None for orbital in orbitals)
    assert all(orbital.occupation == 0.5 for orbital in orbitals)


def test_oh_expansion_of_complete_2p_shell_has_negligible_vzz(
    tmp_path: Path,
) -> None:
    seed_file = tmp_path / "p_shell_seed.txt"
    expanded_file = tmp_path / "p_shell_oh.txt"
    seed_file.write_text(
        "Px 0 0 2 2px - - 1\n"
        "Py 0 0 2 2py - - 1\n"
        "Pz 0 0 2 2pz - - 1\n",
        encoding="utf-8",
    )
    seeds = expand_symmetry.read_seed_records(seed_file, "orbital")
    expanded_file.write_text(
        expand_symmetry.format_expanded_records(
            seeds,
            expand_symmetry.point_group_operations("Oh"),
            position_tolerance=1.0e-9,
        ),
        encoding="utf-8",
    )
    orbitals = read_hydrogenic_orbitals(expanded_file)
    result = calculate_orbital_field(orbitals)
    contribution_scale = sum(
        float(
            np.max(
                np.abs(
                    calculate_orbital_field([orbital]).efg_tensor_v_per_m2
                )
            )
        )
        for orbital in orbitals
    )

    # The radial 2pz seed generates 6 partners. Each transverse seed generates
    # a 12-member orbit; separate seed rows intentionally remain independent.
    assert len(orbitals) == 30
    assert abs(result.parameters.vzz) <= 1.0e-12 * contribution_scale
    assert float(np.max(np.abs(result.efg_tensor_v_per_m2))) <= (
        1.0e-12 * contribution_scale
    )


def test_tolerance_merges_nearly_coincident_positions() -> None:
    operations = expand_symmetry.point_group_operations("C2")
    record = _record((1.0e-10, 0.0, 1.0))

    assert len(
        expand_symmetry.expand_record(
            record, operations, position_tolerance=1.0e-9
        )
    ) == 1
    assert len(
        expand_symmetry.expand_record(
            record, operations, position_tolerance=1.0e-12
        )
    ) == 2


def test_separate_seed_rows_remain_independent(tmp_path: Path) -> None:
    seed_file = tmp_path / "seed.txt"
    seed_file.write_text("A 0 0 1 1\nB 0 0 1 1\n", encoding="utf-8")
    records = expand_symmetry.read_seed_records(seed_file, "charge")
    text = expand_symmetry.format_expanded_records(
        records,
        expand_symmetry.point_group_operations("Ci"),
        position_tolerance=1.0e-9,
    )

    assert len(text.splitlines()) == 4


def test_cli_stdout_and_output_file(tmp_path: Path, capsys: object) -> None:
    seed_file = tmp_path / "seed.txt"
    output_file = tmp_path / "expanded.txt"
    seed_file.write_text("X 0 0 1 -2\n", encoding="utf-8")

    assert expand_symmetry.main(
        ["--group", "Oh", "--mode", "charge", str(seed_file)]
    ) == 0
    captured = capsys.readouterr()  # type: ignore[attr-defined]
    assert len(captured.out.splitlines()) == 6
    assert captured.err == ""

    assert expand_symmetry.main(
        [
            "--group",
            "Oh",
            "--mode",
            "charge",
            str(seed_file),
            "-o",
            str(output_file),
        ]
    ) == 0
    assert len(output_file.read_text(encoding="utf-8").splitlines()) == 6


def test_script_invocation_and_user_error_have_no_traceback(tmp_path: Path) -> None:
    seed_file = tmp_path / "seed.txt"
    seed_file.write_text("X 0 0 1 -2\n", encoding="utf-8")
    script = Path(__file__).parents[1] / "expand_symmetry.py"

    completed = subprocess.run(
        [
            sys.executable,
            str(script),
            "--group",
            "Cinfv",
            "--mode",
            "charge",
            str(seed_file),
        ],
        check=False,
        capture_output=True,
        text=True,
    )

    assert completed.returncode == 2
    assert completed.stdout == ""
    assert "infinite" in completed.stderr
    assert "Traceback" not in completed.stderr


@pytest.mark.parametrize("name", ["S3", "D1", "banana", "C361"])
def test_invalid_group_names_are_rejected(name: str) -> None:
    with pytest.raises(expand_symmetry.SymmetryExpansionError):
        expand_symmetry.point_group_operations(name)
