from __future__ import annotations

import math

import numpy as np
import pytest

from efg_utility import PointCharge, calculate_efg
from efg_utility.point_charge import (
    COULOMB_CONSTANT,
    ELEMENTARY_CHARGE_C,
)

ANGSTROM = 1.0e-10


@pytest.mark.parametrize(
    "position",
    [
        (ANGSTROM, 0.0, 0.0),
        (0.0, ANGSTROM, 0.0),
        (0.0, 0.0, ANGSTROM),
    ],
)
def test_single_positive_charge_on_axis(position: tuple[float, float, float]) -> None:
    result = calculate_efg([PointCharge(position, +1)])
    expected_vzz = 2.0 * COULOMB_CONSTANT * ELEMENTARY_CHARGE_C / ANGSTROM**3

    assert result.vzz == pytest.approx(expected_vzz)
    assert result.eta == 0.0


def test_reversing_charge_reverses_vzz_but_not_eta() -> None:
    positive = calculate_efg([PointCharge((0.0, 0.0, ANGSTROM), +1)])
    negative = calculate_efg([PointCharge((0.0, 0.0, ANGSTROM), -1)])

    assert negative.vzz == pytest.approx(-positive.vzz)
    assert negative.eta == positive.eta == 0.0


def test_off_axis_charge_matches_equation_7() -> None:
    position = np.array((1.0, 2.0, 3.0)) * ANGSTROM
    radius_squared = float(position @ position)
    expected_tensor = (
        COULOMB_CONSTANT
        * ELEMENTARY_CHARGE_C
        * (3.0 * np.outer(position, position) - radius_squared * np.identity(3))
        / radius_squared**2.5
    )
    eigenvalues = np.linalg.eigvalsh(expected_tensor)
    expected_vzz = float(eigenvalues[np.argmax(np.abs(eigenvalues))])

    result = calculate_efg([PointCharge(tuple(position), +1)])

    assert result.vzz == pytest.approx(expected_vzz)
    assert result.eta == pytest.approx(0.0, abs=1.0e-12)


def test_rectangular_configuration_has_expected_principal_parameters() -> None:
    a = 2.0 * ANGSTROM
    b = 1.0 * ANGSTROM
    charges = [
        PointCharge((x, y, 0.0), +1)
        for x in (-a / 2.0, a / 2.0)
        for y in (-b / 2.0, b / 2.0)
    ]
    result = calculate_efg(charges)

    radius_squared = (a / 2.0) ** 2 + (b / 2.0) ** 2
    factor = 4.0 * COULOMB_CONSTANT * ELEMENTARY_CHARGE_C / radius_squared**2.5
    diagonal = np.array(
        (
            factor * (2.0 * (a / 2.0) ** 2 - (b / 2.0) ** 2),
            factor * (-(a / 2.0) ** 2 + 2.0 * (b / 2.0) ** 2),
            factor * (-radius_squared),
        )
    )
    z_index = int(np.argmax(np.abs(diagonal)))
    expected_vzz = float(diagonal[z_index])
    remaining = np.delete(diagonal, z_index)
    y_index = int(np.argmax(np.abs(remaining)))
    expected_vyy = float(remaining[y_index])
    expected_vxx = float(np.delete(remaining, y_index)[0])

    assert result.vzz == pytest.approx(expected_vzz)
    assert result.eta == pytest.approx((expected_vxx - expected_vyy) / expected_vzz)


@pytest.mark.parametrize(
    "charges",
    [
        [
            PointCharge(
                tuple(
                    sign * ANGSTROM if coordinate == axis else 0.0
                    for coordinate in range(3)
                ),
                +1,
            )
            for axis in range(3)
            for sign in (-1.0, 1.0)
        ],
        [
            PointCharge((sx * ANGSTROM, sy * ANGSTROM, sz * ANGSTROM), +1)
            for sx, sy, sz in (
                (+1.0, +1.0, +1.0),
                (-1.0, -1.0, +1.0),
                (+1.0, -1.0, -1.0),
                (-1.0, +1.0, -1.0),
            )
        ],
        [
            PointCharge((sx * ANGSTROM, sy * ANGSTROM, sz * ANGSTROM), +1)
            for sx in (-1.0, 1.0)
            for sy in (-1.0, 1.0)
            for sz in (-1.0, 1.0)
        ],
    ],
)
def test_high_symmetry_configurations_are_zero(charges: list[PointCharge]) -> None:
    result = calculate_efg(charges)

    assert result.vzz == 0.0
    assert result.eta is None


def test_principal_parameters_are_rotation_invariant() -> None:
    charges = [
        PointCharge((0.7 * ANGSTROM, -1.2 * ANGSTROM, 2.0 * ANGSTROM), +2),
        PointCharge((-1.5 * ANGSTROM, 0.2 * ANGSTROM, 0.4 * ANGSTROM), -1),
    ]
    angle = math.radians(37.0)
    rotation = np.array(
        (
            (math.cos(angle), -math.sin(angle), 0.0),
            (math.sin(angle), math.cos(angle), 0.0),
            (0.0, 0.0, 1.0),
        )
    )
    rotated = [
        PointCharge(tuple(rotation @ np.array(charge.position_m)), charge.charge_e)
        for charge in charges
    ]

    original_result = calculate_efg(charges)
    rotated_result = calculate_efg(rotated)

    assert rotated_result.vzz == pytest.approx(original_result.vzz)
    assert rotated_result.eta == pytest.approx(original_result.eta)


def test_principal_parameters_are_axis_permutation_invariant() -> None:
    charges = [
        PointCharge((0.4 * ANGSTROM, -1.3 * ANGSTROM, 2.1 * ANGSTROM), +1),
        PointCharge((-0.8 * ANGSTROM, 0.5 * ANGSTROM, 1.7 * ANGSTROM), -2),
    ]
    permuted = [
        PointCharge(
            (charge.position_m[2], charge.position_m[0], charge.position_m[1]),
            charge.charge_e,
        )
        for charge in charges
    ]

    original_result = calculate_efg(charges)
    permuted_result = calculate_efg(permuted)

    assert permuted_result.vzz == pytest.approx(original_result.vzz)
    assert permuted_result.eta == pytest.approx(original_result.eta)


def test_zero_charge_at_origin_is_permitted() -> None:
    result = calculate_efg([PointCharge((0.0, 0.0, 0.0), 0)])

    assert result.vzz == 0.0
    assert result.eta is None


def test_nonzero_charge_at_origin_is_rejected() -> None:
    with pytest.raises(ValueError, match="singular"):
        calculate_efg([PointCharge((0.0, 0.0, 0.0), 1)])
