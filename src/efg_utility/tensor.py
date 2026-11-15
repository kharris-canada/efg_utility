"""Shared tensor transformations and principal-axis conventions."""

from __future__ import annotations

import math

import numpy as np

from .errors import CalculationError
from .models import EFGResult

_CANCELLATION_RTOL = 1.0e-12
_ORDERING_RTOL = 1.0e-12
_ETA_ENDPOINT_ATOL = 1.0e-12


def _validated_symmetric_tensor(tensor: np.ndarray) -> np.ndarray:
    array = np.asarray(tensor, dtype=np.float64)
    if array.shape != (3, 3) or not np.all(np.isfinite(array)):
        raise CalculationError("an EFG tensor must be a finite 3 by 3 matrix")
    return 0.5 * (array + array.T)


def traceless_tensor(tensor: np.ndarray) -> np.ndarray:
    """Return the symmetric traceless part of a 3 by 3 tensor."""

    symmetric = _validated_symmetric_tensor(tensor)
    result = symmetric - np.identity(3) * (float(np.trace(symmetric)) / 3.0)

    # Removing the residual from one diagonal entry makes the numerical trace
    # exactly zero without changing the physical precision of the result.
    result[2, 2] -= float(np.trace(result))
    return result


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


def efg_parameters_from_tensor(
    tensor: np.ndarray,
    *,
    cancellation_scale: float | None = None,
) -> EFGResult:
    """Return signed ``VZZ`` and ``eta`` for a symmetric EFG tensor.

    ``cancellation_scale`` may be the sum of individual contribution
    magnitudes. It lets a field that cancels to floating-point noise be treated
    as zero without treating every physically small tensor as zero.
    """

    symmetric = _validated_symmetric_tensor(tensor)
    tensor_magnitude = float(np.max(np.abs(symmetric)))
    scale = tensor_magnitude if cancellation_scale is None else cancellation_scale
    if not math.isfinite(scale) or scale < 0.0:
        raise CalculationError("the EFG cancellation scale must be finite and nonnegative")

    if scale == 0.0 or tensor_magnitude <= _CANCELLATION_RTOL * scale:
        return EFGResult(vzz=0.0, eta=None)

    try:
        eigenvalues = np.linalg.eigvalsh(symmetric)
    except np.linalg.LinAlgError as error:
        raise CalculationError(f"could not diagonalize the EFG tensor: {error}") from None

    vxx, vyy, vzz = _principal_components(eigenvalues)
    return EFGResult(vzz=vzz, eta=_normalize_eta((vxx - vyy) / vzz))
