"""Analytic fields from displaced hydrogenic orbitals."""

from __future__ import annotations

import math
import re
from collections.abc import Iterable
from functools import lru_cache

import numpy as np

from .constants import BOHR_RADIUS_M, COULOMB_CONSTANT, ELEMENTARY_CHARGE_C
from .errors import CalculationError
from .models import (
    HydrogenicOrbital,
    OrbitalCombination,
    OrbitalFieldResult,
)
from .tensor import efg_parameters_from_tensor, traceless_tensor

_SMALL_X_THRESHOLD = 1.0e-3
_SERIES_X_THRESHOLD = 0.1
_IDENTITY = np.identity(3)
_P_ORBITAL_AXES = {
    "2px": np.array((1.0, 0.0, 0.0)),
    "2py": np.array((0.0, 1.0, 0.0)),
    "2pz": np.array((0.0, 0.0, 1.0)),
}
_ORBITAL_ALIASES = {"2p_x": "2px", "2p_y": "2py", "2p_z": "2pz"}
_P_DESCRIPTOR_PATTERN = re.compile(
    r"2p\(([^,()]+),([^,()]+),([^,()]+)\)",
    re.IGNORECASE,
)
SUPPORTED_ORBITAL_TYPES = ("1s", "2s", "2px", "2py", "2pz", "2p(ux,uy,uz)")


def normalize_orbital_type(orbital_type: str) -> str:
    """Return the canonical, case-insensitive orbital type name."""

    normalized = orbital_type.casefold()
    return _ORBITAL_ALIASES.get(normalized, normalized)


def _normalized_axis(axis: object) -> np.ndarray:
    try:
        values = np.asarray(axis, dtype=np.float64)
    except (TypeError, ValueError):
        raise CalculationError("a 2p axis must contain 3 finite numbers") from None
    if values.shape != (3,) or not np.all(np.isfinite(values)):
        raise CalculationError("a 2p axis must contain 3 finite numbers")
    magnitude = float(np.linalg.norm(values))
    if not math.isfinite(magnitude) or magnitude == 0.0:
        raise CalculationError("a 2p axis must be nonzero")
    return values / magnitude


def parse_orbital_descriptor(
    descriptor: str,
) -> tuple[str, tuple[float, float, float] | None]:
    """Parse an orbital token into its family and optional normalized p axis."""

    normalized = normalize_orbital_type(descriptor)
    if normalized in ("1s", "2s"):
        return normalized, None
    if normalized in _P_ORBITAL_AXES:
        axis = _P_ORBITAL_AXES[normalized]
        return "2p", tuple(float(value) for value in axis)
    if normalized == "2p":
        return "2p", None

    match = _P_DESCRIPTOR_PATTERN.fullmatch(normalized)
    if match is not None:
        try:
            raw_axis = tuple(float(value) for value in match.groups())
        except ValueError:
            raise CalculationError(
                f"invalid 2p axis in orbital descriptor {descriptor!r}"
            ) from None
        axis = _normalized_axis(raw_axis)
        return "2p", tuple(float(value) for value in axis)

    expected = ", ".join(repr(item) for item in SUPPORTED_ORBITAL_TYPES)
    raise CalculationError(
        f"unsupported orbital type {descriptor!r}; expected one of {expected}"
    )


def radius_from_zeff(zeff: float) -> float:
    """Return the hydrogenic length scale ``a0 / Zeff`` in metres."""

    if not math.isfinite(zeff) or zeff <= 0.0:
        raise CalculationError("effective nuclear charge must be finite and positive")
    radius = BOHR_RADIUS_M / zeff
    if radius == 0.0:
        raise CalculationError("effective nuclear charge produces an unsupported radius")
    return radius


def _validated_orbital(
    orbital: HydrogenicOrbital,
) -> tuple[np.ndarray, str, float, float, np.ndarray | None]:
    orbital_type, descriptor_axis = parse_orbital_descriptor(orbital.orbital_type)
    supplied_axis = (
        _normalized_axis(orbital.axis) if orbital.axis is not None else None
    )
    if orbital_type == "2p":
        if descriptor_axis is not None and supplied_axis is not None:
            raise CalculationError(
                "specify a 2p axis in either orbital_type or axis, not both"
            )
        raw_axis = supplied_axis if supplied_axis is not None else descriptor_axis
        if raw_axis is None:
            raise CalculationError("orbital type '2p' requires an axis")
        axis = np.asarray(raw_axis, dtype=np.float64)
    else:
        if supplied_axis is not None:
            raise CalculationError(f"orbital type {orbital_type!r} does not use an axis")
        axis = None

    center = np.asarray(orbital.center_m, dtype=np.float64)
    if center.shape != (3,) or not np.all(np.isfinite(center)):
        raise CalculationError("an orbital center must contain 3 finite values")
    if not math.isfinite(orbital.radius_m) or orbital.radius_m <= 0.0:
        raise CalculationError("orbital radius must be finite and positive")
    if not math.isfinite(orbital.occupation) or orbital.occupation < 0.0:
        raise CalculationError("orbital occupation must be finite and nonnegative")
    return (
        center,
        orbital_type,
        float(orbital.radius_m),
        float(orbital.occupation),
        axis,
    )


def _potential_shape_small_x(x: float) -> float:
    x2 = x * x
    return 1.0 + x2 * (
        -2.0 / 3.0
        + x
        * (
            2.0 / 3.0
            + x
            * (
                -2.0 / 5.0
                + x
                * (
                    8.0 / 45.0
                    + x * (-4.0 / 63.0 + x * (2.0 / 105.0 - 2.0 * x / 405.0))
                )
            )
        )
    )


def _anisotropy_shape_small_x(x: float) -> float:
    return x * (
        2.0
        + x
        * (
            -16.0 / 5.0
            + x
            * (
                8.0 / 3.0
                + x
                * (-32.0 / 21.0 + x * (2.0 / 3.0 - 32.0 * x / 135.0))
            )
        )
    )


def _potential_shape(x: float) -> float:
    if x == 0.0:
        return 1.0
    if x <= _SMALL_X_THRESHOLD:
        return _potential_shape_small_x(x)
    exponential = math.exp(-2.0 * x)
    return (1.0 - (1.0 + x) * exponential) / x


def _anisotropy_shape(x: float) -> float:
    if x == 0.0:
        return 0.0
    if x <= _SMALL_X_THRESHOLD:
        return _anisotropy_shape_small_x(x)
    inverse_x = 1.0 / x
    inverse_x2 = inverse_x * inverse_x
    inverse_x3 = inverse_x2 * inverse_x
    exponential = math.exp(-2.0 * x)
    return 3.0 * inverse_x3 - exponential * (
        3.0 * inverse_x3 + 6.0 * inverse_x2 + 6.0 * inverse_x + 4.0
    )


def _integer_lower_gamma(order: int, x: float) -> float:
    """Lower incomplete gamma for a positive integer order."""

    if x <= _SERIES_X_THRESHOLD:
        term = 1.0 / order
        total = term
        for index in range(1, 80):
            term *= -x / index * (order + index - 1) / (order + index)
            total += term
            if abs(term) <= abs(total) * 2.0e-16:
                break
        return x**order * total
    return math.factorial(order - 1) - _integer_upper_gamma(order, x)


def _integer_upper_gamma(order: int, x: float) -> float:
    """Upper incomplete gamma for a positive integer order."""

    term = 1.0
    total = term
    for index in range(1, order):
        term *= x / index
        total += term
    return math.factorial(order - 1) * math.exp(-x) * total


def _radial_direction(center: np.ndarray, radius: float) -> tuple[float, np.ndarray]:
    distance = float(np.linalg.norm(center))
    if not math.isfinite(distance):
        raise CalculationError("orbital center is outside the supported numeric range")
    if distance == 0.0:
        return 0.0, np.zeros(3)
    return distance / radius, center / distance


def _hydrogenic_1s_fields(
    orbital: HydrogenicOrbital,
) -> tuple[float, np.ndarray, np.ndarray]:
    center, orbital_type, radius, occupation, _axis = _validated_orbital(orbital)
    if orbital_type != "1s":
        raise CalculationError("the 1s analytic function requires orbital type '1s'")
    electron_charge = -occupation * ELEMENTARY_CHARGE_C
    x, direction = _radial_direction(center, radius)
    try:
        common_scale = COULOMB_CONSTANT * electron_charge / radius**3
        potential = (
            COULOMB_CONSTANT * electron_charge / radius * _potential_shape(x)
        )
    except (OverflowError, ZeroDivisionError):
        raise CalculationError(
            "orbital field is outside the supported numeric range"
        ) from None
    trace = -4.0 * common_scale * math.exp(-2.0 * x)
    anisotropy = common_scale * _anisotropy_shape(x)

    if x == 0.0:
        efg = np.zeros((3, 3), dtype=np.float64)
    else:
        efg = anisotropy * (np.outer(direction, direction) - _IDENTITY / 3.0)
        efg[2, 2] -= float(np.trace(efg))

    hessian = efg + _IDENTITY * (trace / 3.0)
    _ensure_finite_fields(potential, hessian, efg)
    return potential, hessian, efg


def _two_s_shell_coefficient(power: int) -> float:
    """Coefficient of x**power in x^2(2-x)^2 exp(-x) / 8."""

    coefficient = 0.0
    for offset, multiplier in ((2, 4.0), (3, -4.0), (4, 1.0)):
        index = power - offset
        if index >= 0:
            coefficient += multiplier * (-1.0) ** index / math.factorial(index)
    return coefficient / 8.0


def _two_s_shapes(x: float) -> tuple[float, float]:
    """Return dimensionless 2s potential and radial EFG coefficients."""

    if x == 0.0:
        return 0.25, 0.0

    exponential = math.exp(-x)
    outer_integral = exponential * (2.0 + 2.0 * x - x * x + x**3) / 8.0
    if x <= _SERIES_X_THRESHOLD:
        enclosed = 0.0
        anisotropy = 0.0
        for power in range(2, 24):
            coefficient = _two_s_shell_coefficient(power)
            enclosed += coefficient * x ** (power + 1) / (power + 1)
            anisotropy += (
                coefficient
                * (2.0 - power)
                / (power + 1)
                * x ** (power - 2)
            )
    else:
        enclosed = 1.0 - exponential * (
            1.0 + x + x * x / 2.0 + x**4 / 8.0
        )
        shell_density = x * x * (2.0 - x) ** 2 * exponential / 8.0
        anisotropy = 3.0 * enclosed / x**3 - shell_density / x**2
    return enclosed / x + outer_integral, anisotropy


def _hydrogenic_2s_fields(
    orbital: HydrogenicOrbital,
) -> tuple[float, np.ndarray, np.ndarray]:
    center, orbital_type, radius, occupation, _axis = _validated_orbital(orbital)
    if orbital_type != "2s":
        raise CalculationError("the 2s analytic function requires orbital type '2s'")
    electron_charge = -occupation * ELEMENTARY_CHARGE_C
    x, direction = _radial_direction(center, radius)
    try:
        common_scale = COULOMB_CONSTANT * electron_charge / radius**3
        potential_shape, anisotropy_shape = _two_s_shapes(x)
        potential = COULOMB_CONSTANT * electron_charge / radius * potential_shape
    except (OverflowError, ZeroDivisionError):
        raise CalculationError(
            "orbital field is outside the supported numeric range"
        ) from None

    trace = (
        -common_scale * (2.0 - x) ** 2 * math.exp(-x) / 8.0
    )
    if x == 0.0:
        efg = np.zeros((3, 3), dtype=np.float64)
    else:
        efg = common_scale * anisotropy_shape * (
            np.outer(direction, direction) - _IDENTITY / 3.0
        )
        efg[2, 2] -= float(np.trace(efg))
    hessian = efg + _IDENTITY * (trace / 3.0)
    _ensure_finite_fields(potential, hessian, efg)
    return potential, hessian, efg


def _two_p_radial_shapes(
    x: float,
) -> tuple[float, float, float, float, float, float]:
    """Return f, f', f'', g, g', g'' for the 2p potential expansion."""

    if x <= _SERIES_X_THRESHOLD:
        f = 0.25
        f_prime = 0.0
        f_second = 0.0
        g = x * x / 60.0
        g_prime = x / 30.0
        g_second = 1.0 / 30.0
        for index in range(20):
            power = index + 4
            f_coefficient = -(-1.0) ** index / (
                24.0
                * math.factorial(index)
                * (index + 4)
                * (index + 5)
            )
            g_coefficient = -(-1.0) ** index / (
                12.0
                * math.factorial(index)
                * (index + 2)
                * (index + 7)
            )
            f += f_coefficient * x**power
            f_prime += power * f_coefficient * x ** (power - 1)
            f_second += power * (power - 1) * f_coefficient * x ** (power - 2)
            g += g_coefficient * x**power
            g_prime += power * g_coefficient * x ** (power - 1)
            g_second += power * (power - 1) * g_coefficient * x ** (power - 2)
        return f, f_prime, f_second, g, g_prime, g_second

    lower_five = _integer_lower_gamma(5, x)
    lower_seven = _integer_lower_gamma(7, x)
    upper_two = _integer_upper_gamma(2, x)
    exponential = math.exp(-x)
    f = lower_five / (24.0 * x) + _integer_upper_gamma(4, x) / 24.0
    f_prime = -lower_five / (24.0 * x * x)
    f_second = lower_five / (12.0 * x**3) - x * x * exponential / 24.0
    g = lower_seven / (60.0 * x**3) + x * x * upper_two / 60.0
    g_prime = (2.0 * x * upper_two - 3.0 * lower_seven / x**4) / 60.0
    g_second = (
        2.0 * upper_two
        - 5.0 * x * x * exponential
        + 12.0 * lower_seven / x**5
    ) / 60.0
    return f, f_prime, f_second, g, g_prime, g_second


def _two_p_dimensionless_hessian(
    x: float,
    direction: np.ndarray,
    axis: np.ndarray,
    f_prime: float,
    f_second: float,
    g: float,
    g_prime: float,
    g_second: float,
) -> np.ndarray:
    if x == 0.0:
        return (3.0 * np.outer(axis, axis) - _IDENTITY) / 60.0

    projection = float(axis @ direction)
    axial_coordinate = x * projection
    radial_outer = np.outer(direction, direction)
    tangent = _IDENTITY - radial_outer

    if x <= _SERIES_X_THRESHOLD:
        alpha = 1.0 / 40.0
        alpha_prime = 0.0
        alpha_second = 0.0
        for index in range(20):
            power = index + 4
            coefficient = -(-1.0) ** index / (
                12.0
                * math.factorial(index)
                * (index + 2)
                * (index + 7)
            )
            alpha_coefficient = 1.5 * coefficient
            alpha += alpha_coefficient * x ** (power - 2)
            alpha_prime += (
                (power - 2) * alpha_coefficient * x ** (power - 3)
            )
            alpha_second += (
                (power - 2)
                * (power - 3)
                * alpha_coefficient
                * x ** (power - 4)
            )
    else:
        alpha = 1.5 * g / x**2
        alpha_prime = 1.5 * (g_prime / x**2 - 2.0 * g / x**3)
        alpha_second = 1.5 * (
            g_second / x**2 - 4.0 * g_prime / x**3 + 6.0 * g / x**4
        )

    beta_prime = -g_prime / 2.0
    beta_second = -g_second / 2.0
    radial_hessian = (f_second + beta_second) * radial_outer + (
        (f_prime + beta_prime) / x
    ) * tangent
    angular_hessian = (
        alpha_second * axial_coordinate**2 * radial_outer
        + alpha_prime / x * axial_coordinate**2 * tangent
        + 2.0
        * alpha_prime
        * axial_coordinate
        * (np.outer(direction, axis) + np.outer(axis, direction))
        + 2.0 * alpha * np.outer(axis, axis)
    )
    return radial_hessian + angular_hessian


def _hydrogenic_2p_fields(
    orbital: HydrogenicOrbital,
) -> tuple[float, np.ndarray, np.ndarray]:
    center, orbital_type, radius, occupation, axis = _validated_orbital(orbital)
    if orbital_type != "2p" or axis is None:
        raise CalculationError("the 2p analytic function requires a 2p orbital")
    electron_charge = -occupation * ELEMENTARY_CHARGE_C
    x, direction = _radial_direction(center, radius)
    try:
        common_scale = COULOMB_CONSTANT * electron_charge / radius**3
        f, f_prime, f_second, g, g_prime, g_second = _two_p_radial_shapes(x)
        projection = float(axis @ direction) if x != 0.0 else 0.0
        legendre_two = (3.0 * projection**2 - 1.0) / 2.0
        potential = (
            COULOMB_CONSTANT * electron_charge / radius * (f + g * legendre_two)
        )
        dimensionless_hessian = _two_p_dimensionless_hessian(
            x,
            direction,
            axis,
            f_prime,
            f_second,
            g,
            g_prime,
            g_second,
        )
    except (OverflowError, ZeroDivisionError):
        raise CalculationError(
            "orbital field is outside the supported numeric range"
        ) from None

    efg = common_scale * traceless_tensor(dimensionless_hessian)
    trace = (
        -common_scale * x * x * projection**2 * math.exp(-x) / 8.0
        if x != 0.0
        else 0.0
    )
    hessian = efg + _IDENTITY * (trace / 3.0)
    _ensure_finite_fields(potential, hessian, efg)
    return potential, hessian, efg


def _ensure_finite_fields(
    potential: float,
    hessian: np.ndarray,
    efg: np.ndarray,
) -> None:
    if not (
        math.isfinite(potential)
        and np.all(np.isfinite(hessian))
        and np.all(np.isfinite(efg))
    ):
        raise CalculationError("orbital field is outside the supported numeric range")


def _orbital_fields(
    orbital: HydrogenicOrbital,
) -> tuple[float, np.ndarray, np.ndarray]:
    orbital_type, _axis = parse_orbital_descriptor(orbital.orbital_type)
    if orbital_type == "1s":
        return _hydrogenic_1s_fields(orbital)
    if orbital_type == "2s":
        return _hydrogenic_2s_fields(orbital)
    if orbital_type == "2p":
        return _hydrogenic_2p_fields(orbital)
    _validated_orbital(orbital)
    raise AssertionError("validated orbital type has no field backend")


def hydrogenic_orbital_potential(orbital: HydrogenicOrbital) -> float:
    """Return an orbital's analytic electrostatic potential at the origin."""

    potential, _hessian, _efg = _orbital_fields(orbital)
    return potential


def hydrogenic_orbital_hessian(orbital: HydrogenicOrbital) -> np.ndarray:
    """Return an orbital's full potential Hessian in V/m^2."""

    _potential, hessian, _efg = _orbital_fields(orbital)
    return hessian


def hydrogenic_1s_potential(orbital: HydrogenicOrbital) -> float:
    """Return a 1s orbital's analytic electrostatic potential in volts."""

    potential, _hessian, _efg = _hydrogenic_1s_fields(orbital)
    return potential


def hydrogenic_1s_hessian(orbital: HydrogenicOrbital) -> np.ndarray:
    """Return a 1s orbital's full potential Hessian in V/m^2."""

    _potential, hessian, _efg = _hydrogenic_1s_fields(orbital)
    return hessian


def calculate_orbital_field(
    orbitals: Iterable[HydrogenicOrbital],
) -> OrbitalFieldResult:
    """Sum orbital fields and derive the traceless EFG PAS parameters."""

    potential = 0.0
    hessian = np.zeros((3, 3), dtype=np.float64)
    contribution_scale = 0.0

    for orbital in orbitals:
        orbital_potential, orbital_hessian, orbital_efg = _orbital_fields(orbital)
        potential += orbital_potential
        hessian += orbital_hessian
        contribution_scale += float(np.max(np.abs(orbital_efg)))

    if not math.isfinite(potential) or not np.all(np.isfinite(hessian)):
        raise CalculationError("summed orbital field is outside the supported numeric range")

    efg = traceless_tensor(hessian)
    parameters = efg_parameters_from_tensor(
        efg,
        cancellation_scale=contribution_scale,
    )
    return OrbitalFieldResult(
        potential_v=potential,
        hessian_v_per_m2=hessian,
        efg_tensor_v_per_m2=efg,
        parameters=parameters,
    )


def _orbital_wavefunction(orbital: HydrogenicOrbital, points_m: np.ndarray) -> np.ndarray:
    """Evaluate a real normalized 1s, 2s, or 2p AO at Cartesian points."""
    center, orbital_type, radius, _occupation, axis = _validated_orbital(orbital)
    displacement = points_m - center
    radial_distance = np.linalg.norm(displacement, axis=-1)
    scaled = radial_distance / radius
    if orbital_type == "1s":
        return np.exp(-scaled) / (math.sqrt(math.pi) * radius**1.5)
    exponential = np.exp(-0.5 * scaled)
    if orbital_type == "2s":
        return ((2.0 - scaled) * exponential
                / (4.0 * math.sqrt(2.0 * math.pi) * radius**1.5))
    assert axis is not None
    directional = displacement @ axis
    return (directional * exponential
            / (4.0 * math.sqrt(2.0 * math.pi) * radius**2.5))


def _orbital_identity(orbital: HydrogenicOrbital) -> tuple[object, ...]:
    center, orbital_type, radius, _occupation, axis = _validated_orbital(orbital)
    axis_key = None if axis is None else tuple(float(value) for value in axis)
    return (tuple(float(value) for value in center), orbital_type, radius, axis_key)


@lru_cache(maxsize=32768)
def _cross_density_unit_fields(
    identity_a: tuple[object, ...], identity_b: tuple[object, ...]
) -> tuple[float, np.ndarray, np.ndarray, float]:
    """Integrate fields from ``-e * psi_a * psi_b`` in spherical coordinates.

    The quadrupolar kernel is integrated with angular cancellation performed
    at each radius. Gauss-Legendre radial and polar quadrature plus a periodic
    azimuthal rule keep this dependency-free (NumPy only).
    """
    def make_orbital(identity: tuple[object, ...]) -> HydrogenicOrbital:
        center, orbital_type, radius, axis = identity
        return HydrogenicOrbital(
            center_m=center,  # type: ignore[arg-type]
            orbital_type=str(orbital_type),
            radius_m=float(radius),
            axis=axis,  # type: ignore[arg-type]
        )

    orbital_a = make_orbital(identity_a)
    orbital_b = make_orbital(identity_b)
    center_a = np.asarray(orbital_a.center_m, dtype=np.float64)
    center_b = np.asarray(orbital_b.center_m, dtype=np.float64)
    radii = (orbital_a.radius_m, orbital_b.radius_m)
    outer_radius = max(
        float(np.linalg.norm(center_a)), float(np.linalg.norm(center_b))
    ) + 38.0 * max(radii)

    radial_nodes, radial_weights = np.polynomial.legendre.leggauss(8)
    radial_step = 2.0 * min(radii)
    boundaries = [0.0, outer_radius]
    boundaries.extend(
        distance for distance in (
            float(np.linalg.norm(center_a)), float(np.linalg.norm(center_b))
        ) if 0.0 < distance < outer_radius
    )
    boundaries.extend(np.arange(radial_step, outer_radius, radial_step).tolist())
    boundaries = sorted(set(boundaries))
    radial_parts: list[np.ndarray] = []
    weight_parts: list[np.ndarray] = []
    for lower, upper in zip(boundaries[:-1], boundaries[1:], strict=True):
        radial_parts.append(0.5 * (upper - lower) * radial_nodes + 0.5 * (upper + lower))
        weight_parts.append(0.5 * (upper - lower) * radial_weights)
    radial = np.concatenate(radial_parts)
    radial_weights = np.concatenate(weight_parts)
    polar_nodes, polar_weights = np.polynomial.legendre.leggauss(80)
    phi = 2.0 * math.pi * np.arange(160, dtype=np.float64) / 160.0
    mu_grid, phi_grid = np.meshgrid(polar_nodes, phi, indexing="ij")
    angular_weights = np.broadcast_to(
        polar_weights[:, None] * (2.0 * math.pi / phi.size), mu_grid.shape
    ).ravel()
    sin_theta = np.sqrt(np.maximum(0.0, 1.0 - mu_grid**2))
    local_directions = np.column_stack((
        (sin_theta * np.cos(phi_grid)).ravel(),
        (sin_theta * np.sin(phi_grid)).ravel(),
        mu_grid.ravel(),
    ))
    alignment_axis = center_a + center_b
    if float(np.linalg.norm(alignment_axis)) <= 1.0e-14 * max(radii):
        alignment_axis = center_a if np.linalg.norm(center_a) >= np.linalg.norm(center_b) else center_b
    if float(np.linalg.norm(alignment_axis)) == 0.0:
        alignment_axis = np.array((0.0, 0.0, 1.0))
    axis_z = alignment_axis / np.linalg.norm(alignment_axis)
    reference = np.array((0.0, 0.0, 1.0))
    if abs(float(np.dot(axis_z, reference))) > 0.9:
        reference = np.array((1.0, 0.0, 0.0))
    axis_x = np.cross(reference, axis_z)
    axis_x /= np.linalg.norm(axis_x)
    axis_y = np.cross(axis_z, axis_x)
    rotation = np.column_stack((axis_x, axis_y, axis_z))
    directions = local_directions @ rotation.T
    dyad_kernel = 3.0 * directions[:, :, None] * directions[:, None, :]
    dyad_kernel -= np.identity(3)[None, :, :]

    potential_integral = 0.0
    tensor_integral = np.zeros((3, 3), dtype=np.float64)
    overlap_integral = 0.0
    origin_points = np.zeros((1, 3), dtype=np.float64)
    psi_a_origin = float(_orbital_wavefunction(orbital_a, origin_points)[0])
    psi_b_origin = float(_orbital_wavefunction(orbital_b, origin_points)[0])
    for radius, weight in zip(radial, radial_weights, strict=True):
        points = radius * directions
        product = _orbital_wavefunction(orbital_a, points) * _orbital_wavefunction(
            orbital_b, points
        )
        weighted_product = angular_weights * product
        overlap_integral += weight * radius * radius * float(np.sum(weighted_product))
        potential_integral += weight * radius * float(np.sum(weighted_product))
        angular_tensor = np.einsum(
            "a,aij->ij", weighted_product, dyad_kernel, optimize=True
        )
        tensor_integral += (weight / radius) * angular_tensor

    charge_scale = -ELEMENTARY_CHARGE_C * COULOMB_CONSTANT
    potential = charge_scale * potential_integral
    efg = traceless_tensor(charge_scale * tensor_integral)
    trace = 4.0 * math.pi * COULOMB_CONSTANT * ELEMENTARY_CHARGE_C * (
        psi_a_origin * psi_b_origin
    )
    hessian = efg + _IDENTITY * (trace / 3.0)
    _ensure_finite_fields(potential, hessian, efg)
    return potential, hessian, efg, overlap_integral


def calculate_orbital_combinations(
    combinations: Iterable[OrbitalCombination],
) -> OrbitalFieldResult:
    """Calculate the total field including coherent AO interference terms.

    Every combination contributes ``occupation * |sum(c_i psi_i)|**2``.
    Distinct combinations contribute additively; cross terms are included only
    between AO terms sharing the same combination.
    """
    potential = 0.0
    hessian = np.zeros((3, 3), dtype=np.float64)
    contribution_scale = 0.0
    for combination in combinations:
        if not math.isfinite(combination.occupation) or combination.occupation < 0.0:
            raise CalculationError("combination occupation must be finite and nonnegative")
        terms = combination.terms
        if not terms:
            raise CalculationError(f"{combination.name}: combination has no AO terms")
        validated: list[tuple[HydrogenicOrbital, float, tuple[object, ...]]] = []
        for term in terms:
            if not math.isfinite(term.coefficient):
                raise CalculationError("AO coefficient must be finite")
            # Validate even zero-amplitude terms so malformed input is not hidden.
            _validated_orbital(term.orbital)
            validated.append((term.orbital, term.coefficient, _orbital_identity(term.orbital)))

        # Normalize the SALC in the actual nonorthogonal AO basis. The tables'
        # coefficients are normalized in the conventional orthogonal-AO basis;
        # displaced hydrogenic AOs generally have nonzero overlap.
        norm_squared = 0.0
        for index, (_orbital, coefficient, identity) in enumerate(validated):
            norm_squared += coefficient * coefficient
            for _other_orbital, other_coefficient, other_identity in validated[index + 1:]:
                if identity == other_identity:
                    overlap = 1.0
                else:
                    first, second = sorted((identity, other_identity), key=repr)
                    _p, _h, _v, overlap = _cross_density_unit_fields(first, second)
                norm_squared += 2.0 * coefficient * other_coefficient * overlap
        if not math.isfinite(norm_squared) or norm_squared <= 1.0e-14:
            raise CalculationError(
                f"{combination.name}: AO combination has nonpositive or negligible norm"
            )
        density_scale = combination.occupation / norm_squared

        for index, (orbital, coefficient, identity) in enumerate(validated):
            scale = density_scale * coefficient * coefficient
            if scale == 0.0:
                continue
            unit_orbital = HydrogenicOrbital(
                center_m=orbital.center_m,
                orbital_type=orbital.orbital_type,
                radius_m=orbital.radius_m,
                occupation=1.0,
                axis=orbital.axis,
            )
            term_potential, term_hessian, term_efg = _orbital_fields(unit_orbital)
            potential += scale * term_potential
            hessian += scale * term_hessian
            contribution_scale += scale * float(np.max(np.abs(term_efg)))

            for other_orbital, other_coefficient, other_identity in validated[index + 1:]:
                pair_scale = 2.0 * density_scale * coefficient * other_coefficient
                if pair_scale == 0.0:
                    continue
                if identity == other_identity:
                    pair_potential, pair_hessian, pair_efg = term_potential, term_hessian, term_efg
                else:
                    first, second = sorted((identity, other_identity), key=repr)
                    pair_potential, pair_hessian, pair_efg, _overlap = _cross_density_unit_fields(first, second)
                potential += pair_scale * pair_potential
                hessian += pair_scale * pair_hessian
                contribution_scale += abs(pair_scale) * float(np.max(np.abs(pair_efg)))

    if not math.isfinite(potential) or not np.all(np.isfinite(hessian)):
        raise CalculationError("summed orbital field is outside the supported numeric range")
    hessian = 0.5 * (hessian + hessian.T)
    efg = traceless_tensor(hessian)
    parameters = efg_parameters_from_tensor(efg, cancellation_scale=contribution_scale)
    return OrbitalFieldResult(potential, hessian, efg, parameters)
