"""Numerical point-charge EFG calculation."""

from __future__ import annotations

import math
from collections.abc import Iterable

import numpy as np

from .constants import (
    COULOMB_CONSTANT,
    ELEMENTARY_CHARGE_C,
)
from .constants import (
    VACUUM_PERMITTIVITY_F_PER_M as _VACUUM_PERMITTIVITY_F_PER_M,
)
from .errors import CalculationError
from .models import EFGResult, PointCharge
from .tensor import efg_parameters_from_tensor

# Preserve the constants that were available from this module in version 0.1.
VACUUM_PERMITTIVITY_F_PER_M = _VACUUM_PERMITTIVITY_F_PER_M


def _point_charge_tensor_and_scale(
    charges: Iterable[PointCharge],
) -> tuple[np.ndarray, float]:
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

    return 0.5 * (tensor + tensor.T), contribution_scale


def calculate_point_charge_tensor(charges: Iterable[PointCharge]) -> np.ndarray:
    """Calculate the traceless point-charge EFG tensor at the origin."""

    tensor, _contribution_scale = _point_charge_tensor_and_scale(charges)
    return tensor


def calculate_efg(charges: Iterable[PointCharge]) -> EFGResult:
    """Calculate signed ``VZZ`` and ``eta`` at the origin.

    Positions must be supplied in metres and charges in integer multiples of
    the elementary charge. A nonzero charge at the origin is singular.
    """

    tensor, contribution_scale = _point_charge_tensor_and_scale(charges)
    return efg_parameters_from_tensor(
        tensor,
        cancellation_scale=contribution_scale,
    )
