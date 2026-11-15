from __future__ import annotations

import math

import numpy as np
import pytest

from efg_utility import (
    BOHR_RADIUS_M,
    HydrogenicOrbital,
    PointCharge,
    calculate_orbital_field,
    calculate_point_charge_tensor,
    hydrogenic_1s_hessian,
    hydrogenic_1s_potential,
    hydrogenic_orbital_hessian,
    hydrogenic_orbital_potential,
    radius_from_zeff,
)
from efg_utility.constants import COULOMB_CONSTANT, ELEMENTARY_CHARGE_C


def _orbital(
    center_m: tuple[float, float, float],
    *,
    radius_m: float = BOHR_RADIUS_M,
    occupation: float = 1.0,
) -> HydrogenicOrbital:
    return HydrogenicOrbital(center_m, "1s", radius_m, occupation)


def _typed_orbital(
    orbital_type: str,
    center_m: tuple[float, float, float],
    *,
    radius_m: float = BOHR_RADIUS_M,
    occupation: float = 1.0,
) -> HydrogenicOrbital:
    return HydrogenicOrbital(center_m, orbital_type, radius_m, occupation)


def _numerical_density_normalization() -> float:
    nodes, weights = np.polynomial.legendre.leggauss(240)
    radial = (nodes + 1.0) * 12.5
    radial_weights = weights * 12.5
    return float(np.sum(radial_weights * 4.0 * radial**2 * np.exp(-2.0 * radial)))


def _numerical_axial_efg_zz(
    distance_over_radius: float,
    radius_m: float,
    occupation: float = 1.0,
) -> float:
    """Test-only spherical quadrature of the traceless EFG kernel.

    Coordinates are centered at the EFG observation point. Gauss-Legendre
    angular integration performs the cancellation of the singular-looking
    rank-2 kernel, and radial integration is split at the orbital center.
    """

    radial_nodes, radial_weights = np.polynomial.legendre.leggauss(500)
    mu, mu_weights = np.polynomial.legendre.leggauss(400)
    ratio = distance_over_radius

    def integrate_interval(lower: float, upper: float) -> float:
        radial = lower + (radial_nodes + 1.0) * (upper - lower) / 2.0
        weights = radial_weights * (upper - lower) / 2.0
        orbital_radius = np.sqrt(
            radial[:, None] ** 2
            + ratio**2
            - 2.0 * radial[:, None] * ratio * mu[None, :]
        )
        angular = (
            np.exp(-2.0 * orbital_radius)
            * (3.0 * mu[None, :] ** 2 - 1.0)
        ) @ mu_weights
        return float(np.sum(weights * angular / radial))

    cutoff = max(25.0, ratio + 25.0)
    dimensionless = 2.0 * (
        integrate_interval(0.0, ratio)
        + integrate_interval(ratio, cutoff)
    )
    electron_charge = -occupation * ELEMENTARY_CHARGE_C
    return COULOMB_CONSTANT * electron_charge / radius_m**3 * dimensionless


def _numerical_n2_axial_efg_zz(
    orbital_type: str,
    distance_over_radius: float,
    radius_m: float,
) -> float:
    """Independent spherical quadrature for 2s or axial 2pz density."""

    radial_nodes, radial_weights = np.polynomial.legendre.leggauss(500)
    mu, mu_weights = np.polynomial.legendre.leggauss(400)
    ratio = distance_over_radius

    def integrate_interval(lower: float, upper: float) -> float:
        radial = lower + (radial_nodes + 1.0) * (upper - lower) / 2.0
        weights = radial_weights * (upper - lower) / 2.0
        relative_radius = np.sqrt(
            radial[:, None] ** 2
            + ratio**2
            - 2.0 * radial[:, None] * ratio * mu[None, :]
        )
        kernel = 3.0 * mu[None, :] ** 2 - 1.0
        if orbital_type == "2s":
            density_shape = (2.0 - relative_radius) ** 2 * np.exp(
                -relative_radius
            )
        elif orbital_type == "2px":
            density_shape = (
                0.5
                * radial[:, None] ** 2
                * (1.0 - mu[None, :] ** 2)
                * np.exp(-relative_radius)
            )
        else:
            relative_z = radial[:, None] * mu[None, :] - ratio
            density_shape = relative_z**2 * np.exp(-relative_radius)
        angular = (density_shape * kernel) @ mu_weights
        return float(np.sum(weights * angular / radial))

    cutoff = max(40.0, ratio + 40.0)
    dimensionless = (
        integrate_interval(0.0, ratio)
        + integrate_interval(ratio, cutoff)
    ) / 16.0
    common = -COULOMB_CONSTANT * ELEMENTARY_CHARGE_C / radius_m**3
    return common * dimensionless


def test_charge_density_is_numerically_normalized() -> None:
    assert _numerical_density_normalization() == pytest.approx(1.0, abs=2.0e-12)


def test_origin_has_finite_potential_isotropic_hessian_and_zero_efg() -> None:
    orbital = _orbital((0.0, 0.0, 0.0))
    result = calculate_orbital_field([orbital])
    electron_charge = -ELEMENTARY_CHARGE_C
    expected_potential = COULOMB_CONSTANT * electron_charge / BOHR_RADIUS_M
    expected_trace = (
        -4.0 * COULOMB_CONSTANT * electron_charge / BOHR_RADIUS_M**3
    )

    assert hydrogenic_1s_potential(orbital) == pytest.approx(expected_potential)
    assert result.potential_v == pytest.approx(expected_potential)
    assert hydrogenic_1s_hessian(orbital) == pytest.approx(
        np.identity(3) * expected_trace / 3.0
    )
    assert result.hessian_v_per_m2 == pytest.approx(
        np.identity(3) * expected_trace / 3.0
    )
    assert result.efg_tensor_v_per_m2 == pytest.approx(np.zeros((3, 3)))
    assert result.parameters.vzz == 0.0
    assert result.parameters.eta is None


@pytest.mark.parametrize("ratio", [0.5, 1.0, 2.0, 5.0])
def test_full_hessian_trace_obeys_poisson_equation(ratio: float) -> None:
    radius = BOHR_RADIUS_M
    orbital = _orbital((0.0, 0.0, ratio * radius))
    hessian = hydrogenic_1s_hessian(orbital)
    electron_charge = -ELEMENTARY_CHARGE_C
    expected_trace = (
        -4.0
        * COULOMB_CONSTANT
        * electron_charge
        / radius**3
        * math.exp(-2.0 * ratio)
    )

    assert float(np.trace(hessian)) == pytest.approx(expected_trace)


def test_small_displacement_uses_stable_series() -> None:
    ratio = 1.0e-6
    radius = BOHR_RADIUS_M
    result = calculate_orbital_field([_orbital((0.0, 0.0, ratio * radius))])
    common = -COULOMB_CONSTANT * ELEMENTARY_CHARGE_C / radius**3
    leading_efg_zz = common * (4.0 / 3.0) * ratio

    assert result.efg_tensor_v_per_m2[2, 2] == pytest.approx(
        leading_efg_zz,
        rel=2.0e-6,
    )
    assert np.all(np.isfinite(result.hessian_v_per_m2))


def test_large_distance_matches_point_electron_tensor() -> None:
    radius = BOHR_RADIUS_M
    center = (0.0, 0.0, 50.0 * radius)
    orbital_result = calculate_orbital_field([_orbital(center)])
    point_tensor = calculate_point_charge_tensor([PointCharge(center, -1)])

    assert orbital_result.efg_tensor_v_per_m2 == pytest.approx(
        point_tensor,
        rel=2.0e-12,
    )


def test_axial_symmetry() -> None:
    result = calculate_orbital_field(
        [_orbital((0.0, 0.0, 2.0 * BOHR_RADIUS_M))]
    )
    tensor = result.efg_tensor_v_per_m2

    assert tensor[0, 0] == pytest.approx(tensor[1, 1])
    assert tensor - np.diag(np.diag(tensor)) == pytest.approx(np.zeros((3, 3)))
    assert result.parameters.eta == 0.0


def test_rotational_covariance() -> None:
    radius = BOHR_RADIUS_M
    axial = calculate_orbital_field([_orbital((0.0, 0.0, 2.0 * radius))])
    direction = np.array((1.0, 2.0, 3.0))
    direction /= np.linalg.norm(direction)
    first_axis = np.cross((0.0, 0.0, 1.0), direction)
    first_axis /= np.linalg.norm(first_axis)
    second_axis = np.cross(direction, first_axis)
    rotation = np.column_stack((first_axis, second_axis, direction))
    rotated = calculate_orbital_field(
        [_orbital(tuple(2.0 * radius * direction))]
    )

    assert rotated.hessian_v_per_m2 == pytest.approx(
        rotation @ axial.hessian_v_per_m2 @ rotation.T,
        rel=2.0e-14,
    )
    assert rotated.efg_tensor_v_per_m2 == pytest.approx(
        rotation @ axial.efg_tensor_v_per_m2 @ rotation.T,
        rel=2.0e-14,
    )


@pytest.mark.parametrize("ratio", [0.5, 1.0, 2.0, 5.0, 10.0])
def test_analytic_efg_matches_independent_spherical_quadrature(ratio: float) -> None:
    radius = BOHR_RADIUS_M
    result = calculate_orbital_field(
        [_orbital((0.0, 0.0, ratio * radius))]
    )
    numerical_zz = _numerical_axial_efg_zz(ratio, radius)

    assert result.efg_tensor_v_per_m2[2, 2] == pytest.approx(
        numerical_zz,
        rel=5.0e-6,
    )


def test_multiple_orbitals_sum_before_tensor_analysis() -> None:
    first = _orbital((0.0, 0.0, BOHR_RADIUS_M), occupation=0.5)
    second = _orbital((BOHR_RADIUS_M, 0.0, 0.0), occupation=1.5)
    combined = calculate_orbital_field([first, second])
    first_result = calculate_orbital_field([first])
    second_result = calculate_orbital_field([second])

    assert combined.potential_v == pytest.approx(
        first_result.potential_v + second_result.potential_v
    )
    assert combined.hessian_v_per_m2 == pytest.approx(
        first_result.hessian_v_per_m2 + second_result.hessian_v_per_m2
    )


def test_radius_from_zeff() -> None:
    assert radius_from_zeff(2.0) == pytest.approx(BOHR_RADIUS_M / 2.0)


def test_n2_densities_are_numerically_normalized() -> None:
    nodes, weights = np.polynomial.legendre.leggauss(240)
    radial = (nodes + 1.0) * 25.0
    radial_weights = weights * 25.0
    two_s = np.sum(
        radial_weights
        * radial**2
        * (2.0 - radial) ** 2
        * np.exp(-radial)
        / 8.0
    )
    two_p = np.sum(radial_weights * radial**4 * np.exp(-radial) / 24.0)

    assert float(two_s) == pytest.approx(1.0, abs=2.0e-12)
    assert float(two_p) == pytest.approx(1.0, abs=2.0e-12)


@pytest.mark.parametrize("orbital_type", ["2s", "2px", "2py", "2pz"])
def test_n2_large_distance_recovers_point_charge(orbital_type: str) -> None:
    radius = BOHR_RADIUS_M
    center = (0.0, 0.0, 1.0e5 * radius)
    result = calculate_orbital_field([_typed_orbital(orbital_type, center)])
    point_tensor = calculate_point_charge_tensor([PointCharge(center, -1)])

    assert result.efg_tensor_v_per_m2 == pytest.approx(
        point_tensor,
        rel=2.0e-8,
        abs=1.0e-4,
    )


def test_2s_origin_limit_and_poisson_trace() -> None:
    radius = BOHR_RADIUS_M
    orbital = _typed_orbital("2s", (0.0, 0.0, 0.0))
    result = calculate_orbital_field([orbital])
    common = -COULOMB_CONSTANT * ELEMENTARY_CHARGE_C / radius**3

    assert hydrogenic_orbital_potential(orbital) == pytest.approx(
        -COULOMB_CONSTANT * ELEMENTARY_CHARGE_C / (4.0 * radius)
    )
    assert hydrogenic_orbital_hessian(orbital) == pytest.approx(
        np.identity(3) * (-common / 2.0) / 3.0
    )
    assert result.efg_tensor_v_per_m2 == pytest.approx(np.zeros((3, 3)))
    assert result.parameters.eta is None


@pytest.mark.parametrize("ratio", [0.5, 1.0, 2.0, 5.0, 10.0])
def test_2s_analytic_efg_matches_independent_quadrature(ratio: float) -> None:
    radius = BOHR_RADIUS_M
    result = calculate_orbital_field(
        [_typed_orbital("2s", (0.0, 0.0, ratio * radius))]
    )

    assert result.efg_tensor_v_per_m2[2, 2] == pytest.approx(
        _numerical_n2_axial_efg_zz("2s", ratio, radius),
        rel=8.0e-6,
    )


@pytest.mark.parametrize("ratio", [0.5, 1.0, 2.0, 5.0, 10.0])
def test_2pz_analytic_efg_matches_independent_quadrature(ratio: float) -> None:
    radius = BOHR_RADIUS_M
    result = calculate_orbital_field(
        [_typed_orbital("2pz", (0.0, 0.0, ratio * radius))]
    )

    assert result.efg_tensor_v_per_m2[2, 2] == pytest.approx(
        _numerical_n2_axial_efg_zz("2pz", ratio, radius),
        rel=8.0e-6,
    )


@pytest.mark.parametrize("ratio", [0.5, 1.0, 2.0, 5.0, 10.0])
def test_perpendicular_2px_efg_matches_independent_quadrature(ratio: float) -> None:
    radius = BOHR_RADIUS_M
    result = calculate_orbital_field(
        [_typed_orbital("2px", (0.0, 0.0, ratio * radius))]
    )

    assert result.efg_tensor_v_per_m2[2, 2] == pytest.approx(
        _numerical_n2_axial_efg_zz("2px", ratio, radius),
        rel=8.0e-6,
    )


def test_2p_origin_has_intrinsic_quadrupolar_efg() -> None:
    radius = BOHR_RADIUS_M
    common = -COULOMB_CONSTANT * ELEMENTARY_CHARGE_C / radius**3
    result = calculate_orbital_field(
        [_typed_orbital("2pz", (0.0, 0.0, 0.0))]
    )
    expected = common * np.diag((-1.0, -1.0, 2.0)) / 60.0

    assert result.potential_v == pytest.approx(
        -COULOMB_CONSTANT * ELEMENTARY_CHARGE_C / (4.0 * radius)
    )
    assert result.hessian_v_per_m2 == pytest.approx(expected)
    assert result.efg_tensor_v_per_m2 == pytest.approx(expected)
    assert float(np.trace(result.hessian_v_per_m2)) == pytest.approx(0.0, abs=1.0)
    assert result.parameters.eta == 0.0


def test_2p_axis_and_center_rotate_together() -> None:
    radius = BOHR_RADIUS_M
    z_result = calculate_orbital_field(
        [_typed_orbital("2pz", (0.0, 0.0, 2.0 * radius))]
    )
    x_result = calculate_orbital_field(
        [_typed_orbital("2px", (2.0 * radius, 0.0, 0.0))]
    )
    rotation = np.array(((0.0, 0.0, 1.0), (0.0, 1.0, 0.0), (1.0, 0.0, 0.0)))

    assert x_result.hessian_v_per_m2 == pytest.approx(
        rotation @ z_result.hessian_v_per_m2 @ rotation.T
    )
    assert x_result.efg_tensor_v_per_m2 == pytest.approx(
        rotation @ z_result.efg_tensor_v_per_m2 @ rotation.T
    )


def test_arbitrary_2p_axis_is_rotationally_covariant() -> None:
    radius = BOHR_RADIUS_M
    center = np.array((0.3, -0.7, 1.1)) * radius
    axis = np.array((1.0, 2.0, -0.5))
    axis /= np.linalg.norm(axis)
    angle = 0.713
    rotation = np.array(
        (
            (math.cos(angle), -math.sin(angle), 0.0),
            (math.sin(angle), math.cos(angle), 0.0),
            (0.0, 0.0, 1.0),
        )
    )
    original = calculate_orbital_field(
        [HydrogenicOrbital(tuple(center), "2p", radius, 1.0, tuple(axis))]
    )
    rotated = calculate_orbital_field(
        [
            HydrogenicOrbital(
                tuple(rotation @ center),
                "2p",
                radius,
                1.0,
                tuple(rotation @ axis),
            )
        ]
    )

    assert rotated.potential_v == pytest.approx(original.potential_v)
    assert rotated.hessian_v_per_m2 == pytest.approx(
        rotation @ original.hessian_v_per_m2 @ rotation.T,
        rel=3.0e-14,
    )


def test_arbitrary_axis_descriptor_matches_programmatic_axis() -> None:
    radius = BOHR_RADIUS_M
    center = (0.2 * radius, 0.4 * radius, radius)
    descriptor = calculate_orbital_field(
        [HydrogenicOrbital(center, "2p(1,2,3)", radius)]
    )
    programmatic = calculate_orbital_field(
        [HydrogenicOrbital(center, "2p", radius, 1.0, (1.0, 2.0, 3.0))]
    )

    assert descriptor.hessian_v_per_m2 == pytest.approx(
        programmatic.hessian_v_per_m2
    )


@pytest.mark.parametrize("orbital_type", ["2s", "2px", "2py", "2pz"])
def test_n2_occupation_scales_fields_and_zero_is_empty(orbital_type: str) -> None:
    center = (0.4 * BOHR_RADIUS_M, 0.2 * BOHR_RADIUS_M, BOHR_RADIUS_M)
    full = calculate_orbital_field([_typed_orbital(orbital_type, center)])
    quarter = calculate_orbital_field(
        [_typed_orbital(orbital_type, center, occupation=0.25)]
    )
    empty = calculate_orbital_field(
        [_typed_orbital(orbital_type, center, occupation=0.0)]
    )

    assert quarter.potential_v == pytest.approx(full.potential_v * 0.25)
    assert quarter.hessian_v_per_m2 == pytest.approx(
        full.hessian_v_per_m2 * 0.25
    )
    assert empty.potential_v == 0.0
    assert empty.hessian_v_per_m2 == pytest.approx(np.zeros((3, 3)))
    assert empty.parameters.eta is None


@pytest.mark.parametrize(
    ("orbital_type", "center"),
    [
        ("2s", (0.0, 0.0, 1.3)),
        ("2px", (0.7, 0.2, 1.1)),
        ("2py", (0.7, 0.2, 1.1)),
        ("2pz", (0.7, 0.2, 1.1)),
    ],
)
def test_n2_hessian_trace_obeys_poisson_equation(
    orbital_type: str,
    center: tuple[float, float, float],
) -> None:
    radius = BOHR_RADIUS_M
    center_m = tuple(component * radius for component in center)
    orbital = _typed_orbital(orbital_type, center_m)
    hessian = hydrogenic_orbital_hessian(orbital)
    common = -COULOMB_CONSTANT * ELEMENTARY_CHARGE_C / radius**3
    x = float(np.linalg.norm(center))
    if orbital_type == "2s":
        expected_trace = -common * (2.0 - x) ** 2 * math.exp(-x) / 8.0
    else:
        axis_index = {"2px": 0, "2py": 1, "2pz": 2}[orbital_type]
        expected_trace = (
            -common * center[axis_index] ** 2 * math.exp(-x) / 8.0
        )

    assert float(np.trace(hessian)) == pytest.approx(expected_trace, rel=2.0e-14)
    assert float(np.trace(calculate_orbital_field([orbital]).efg_tensor_v_per_m2)) == pytest.approx(
        0.0,
        abs=2.0e5,
    )


@pytest.mark.parametrize(
    "orbital",
    [
        HydrogenicOrbital((0.0, 0.0, 1.0), "3s", BOHR_RADIUS_M, 1.0),
        HydrogenicOrbital((0.0, 0.0, 1.0), "1s", 0.0, 1.0),
        HydrogenicOrbital((0.0, 0.0, 1.0), "1s", BOHR_RADIUS_M, -1.0),
    ],
)
def test_invalid_orbital_parameters_are_rejected(
    orbital: HydrogenicOrbital,
) -> None:
    with pytest.raises(ValueError):
        calculate_orbital_field([orbital])
