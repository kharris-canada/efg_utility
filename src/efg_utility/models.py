"""Public data types used by EFG calculations."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


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


@dataclass(frozen=True, slots=True)
class HydrogenicOrbital:
    """A normalized hydrogenic orbital centered away from the origin.

    Positions and ``radius_m`` use SI metres. ``occupation`` multiplies the
    negative electron charge density and may be fractional. A general real
    2p orbital uses ``orbital_type="2p"`` and a nonzero Cartesian ``axis``;
    the axis is normalized by the calculation backend.
    """

    center_m: tuple[float, float, float]
    orbital_type: str
    radius_m: float
    occupation: float = 1.0
    axis: tuple[float, float, float] | None = None


@dataclass(frozen=True, slots=True)
class OrbitalCombinationTerm:
    """One signed AO amplitude in a coherent orbital combination."""

    orbital: HydrogenicOrbital
    coefficient: float


@dataclass(frozen=True, slots=True)
class OrbitalCombination:
    """A coherent linear combination with a total electron occupation."""

    terms: tuple[OrbitalCombinationTerm, ...]
    occupation: float = 2.0
    name: str = "combination"


@dataclass(frozen=True, slots=True)
class OrbitalFieldResult:
    """Potential, full Hessian, and traceless EFG for one or more orbitals."""

    potential_v: float
    hessian_v_per_m2: np.ndarray
    efg_tensor_v_per_m2: np.ndarray
    parameters: EFGResult
