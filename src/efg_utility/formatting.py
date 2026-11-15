"""Human-readable formatting for EFG results."""

from __future__ import annotations

import numpy as np

from .models import EFGResult, OrbitalFieldResult


def format_result(result: EFGResult) -> str:
    """Format an EFG result using the stable version-1 terminal contract."""

    eta = "undefined" if result.eta is None else f"{result.eta:.12g}"
    return f"VZZ = {result.vzz:+.12e} V/m^2\neta = {eta}"


def _format_matrix(tensor: np.ndarray) -> str:
    return "\n".join(
        "[" + " ".join(f"{float(value):+.12e}" for value in row) + "]"
        for row in tensor
    )


def format_orbital_result(result: OrbitalFieldResult) -> str:
    """Format the full and traceless orbital-field tensors."""

    eta = (
        "undefined"
        if result.parameters.eta is None
        else f"{result.parameters.eta:.12g}"
    )
    trace = float(np.trace(result.hessian_v_per_m2))
    return "\n".join(
        (
            f"potential = {result.potential_v:+.12e} V",
            "full Hessian (V/m^2):",
            _format_matrix(result.hessian_v_per_m2),
            f"trace = {trace:+.12e} V/m^2",
            "traceless EFG (V/m^2):",
            _format_matrix(result.efg_tensor_v_per_m2),
            f"VZZ = {result.parameters.vzz:+.12e} V/m^2",
            f"eta = {eta}",
        )
    )
