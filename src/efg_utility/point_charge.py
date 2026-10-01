"""Numerical point-charge EFG calculation."""

from __future__ import annotations

from collections.abc import Iterable
import math

import numpy as np

from .models import EFGResult, PointCharge


# Elementary charge is exact in the SI. Vacuum permittivity is the 2022 CODATA
# recommended value in F/m.
ELEMENTARY_CHARGE_C = 1.602176634e-19
VACUUM_PERMITTIVITY_F_PER_M = 8.8541878188e-12
COULOMB_CONSTANT = 1.0 / (4.0 * math.pi * VACUUM_PERMITTIVITY_F_PER_M)

_CANCELLATION_RTOL = 1.0e-12
_ORDERING_RTOL = 1.0e-12
_ETA_ENDPOINT_ATOL = 1.0e-12


class CalculationError(ValueError):
    """Raised when a requested EFG calculation is singular or non-finite."""


def _magnitude_tied(value: float, reference: float) -> bool:
    return math.isclose(
        abs(value),
        abs(reference),
        rel_tol=_ORDERING_RTOL,
        abs_tol=0.0,
    )


def _largest_magnitude_index(values: np.ndarray, indices: list[int]) -> int:
    """Return the largest-magnitude index, preferring positive values on ties."""

    largest = max(indices, key=lambda index: abs(float(values[index])))
    tied = [
        index
        for index in indices
        if _magnitude_tied(float(values[index]), float(values[largest]))
    ]
    return max(tied, key=lambda index: float(values[index]))


def _principal_components(eigenvalues: np.ndarray) -> tuple[float, float, float]:
    indices = [0, 1, 2]
    z_index = _largest_magnitude_index(eigenvalues, indices)
    indices.remove(z_index)
    y_index = _largest_magnitude_index(eigenvalues, indices)
    indices.remove(y_index)
    x_index = indices[0]
    return (
        float(eigenvalues[x_index]),
        float(eigenvalues[y_index]),
        float(eigenvalues[z_index]),
    )


def _normalize_eta(eta: float) -> float:
    if math.isclose(eta, 0.0, rel_tol=0.0, abs_tol=_ETA_ENDPOINT_ATOL):
        return 0.0
    if math.isclose(eta, 1.0, rel_tol=0.0, abs_tol=_ETA_ENDPOINT_ATOL):
        return 1.0
    return eta


def calculate_efg(charges: Iterable[PointCharge]) -> EFGResult:
    """Calculate signed ``VZZ`` and ``eta`` at the origin.

    Positions must be supplied in metres and charges in integer multiples of
    the elementary charge. A nonzero charge at the origin is singular.
    """

    tensor = np.zeros((3, 3), dtype=np.float64)
    contribution_scale = 0.0

    for charge in charges:
        if charge.charge_e == 0:
            continue

        position = np.asarray(charge.position_m, dtype=np.float64)
        if position.shape != (3,) or not np.all(np.isfinite(position)):
            raise CalculationError("each point-charge position must contain 3 finite values")

        radius_squared = float(position @ position)
        if radius_squared == 0.0:
            raise CalculationError("a nonzero point charge at the origin is singular")
        if not math.isfinite(radius_squared):
            raise CalculationError("point-charge position is outside the supported numeric range")

        try:
            with np.errstate(over="raise", divide="raise", invalid="raise"):
                numerator = 3.0 * np.outer(position, position)
                numerator -= radius_squared * np.identity(3)
                scale = (
                    COULOMB_CONSTANT
                    * (charge.charge_e * ELEMENTARY_CHARGE_C)
                    / (radius_squared ** 2.5)
                )
                contribution = scale * numerator
        except (FloatingPointError, OverflowError):
            raise CalculationError(
                "point-charge contribution is outside the supported numeric range"
            ) from None

        if not np.all(np.isfinite(contribution)):
            raise CalculationError(
                "point-charge contribution is outside the supported numeric range"
            )

        try:
            with np.errstate(over="raise", invalid="raise"):
                tensor += contribution
        except FloatingPointError:
            raise CalculationError(
                "summed EFG tensor is outside the supported numeric range"
            ) from None

        contribution_scale += float(np.max(np.abs(contribution)))
        if not math.isfinite(contribution_scale):
            raise CalculationError(
                "summed EFG tensor is outside the supported numeric range"
            )

    tensor = 0.5 * (tensor + tensor.T)
    tensor_magnitude = float(np.max(np.abs(tensor)))

    if contribution_scale == 0.0 or (
        tensor_magnitude <= _CANCELLATION_RTOL * contribution_scale
    ):
        return EFGResult(vzz=0.0, eta=None)

    try:
        eigenvalues = np.linalg.eigvalsh(tensor)
    except np.linalg.LinAlgError as error:
        raise CalculationError(f"could not diagonalize the EFG tensor: {error}") from None

    vxx, vyy, vzz = _principal_components(eigenvalues)
    eta = _normalize_eta((vxx - vyy) / vzz)
    return EFGResult(vzz=vzz, eta=eta)
