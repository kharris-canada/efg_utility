"""Public data types used by EFG calculations."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class PointCharge:
    """A point charge whose position is expressed in metres.

    ``charge_e`` is an integer multiple of the elementary charge.
    """

    position_m: tuple[float, float, float]
    charge_e: int


@dataclass(frozen=True, slots=True)
class EFGResult:
    """The orientation-independent principal-axis EFG parameters."""

    vzz: float
    eta: float | None
