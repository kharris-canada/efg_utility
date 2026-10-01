"""Human-readable formatting for EFG results."""

from __future__ import annotations

from .models import EFGResult


def format_result(result: EFGResult) -> str:
    """Format an EFG result using the stable version-1 terminal contract."""

    eta = "undefined" if result.eta is None else f"{result.eta:.12g}"
    return f"VZZ = {result.vzz:+.12e} V/m^2\neta = {eta}"
