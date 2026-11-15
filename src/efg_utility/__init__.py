"""Point-charge and hydrogenic-orbital electric field gradients."""

from .constants import BOHR_RADIUS_M
from .errors import CalculationError
from .models import (
    EFGResult,
    HydrogenicOrbital,
    OrbitalCombination,
    OrbitalCombinationTerm,
    OrbitalFieldResult,
    PointCharge,
)
from .orbital import (
    SUPPORTED_ORBITAL_TYPES,
    calculate_orbital_field,
    calculate_orbital_combinations,
    hydrogenic_1s_hessian,
    hydrogenic_1s_potential,
    hydrogenic_orbital_hessian,
    hydrogenic_orbital_potential,
    normalize_orbital_type,
    parse_orbital_descriptor,
    radius_from_zeff,
)
from .point_charge import calculate_efg, calculate_point_charge_tensor
from .tensor import efg_parameters_from_tensor, traceless_tensor

__all__ = [
    "BOHR_RADIUS_M",
    "SUPPORTED_ORBITAL_TYPES",
    "CalculationError",
    "EFGResult",
    "HydrogenicOrbital",
    "OrbitalCombination",
    "OrbitalCombinationTerm",
    "OrbitalFieldResult",
    "PointCharge",
    "calculate_efg",
    "calculate_orbital_field",
    "calculate_orbital_combinations",
    "calculate_point_charge_tensor",
    "efg_parameters_from_tensor",
    "hydrogenic_1s_hessian",
    "hydrogenic_1s_potential",
    "hydrogenic_orbital_hessian",
    "hydrogenic_orbital_potential",
    "normalize_orbital_type",
    "parse_orbital_descriptor",
    "radius_from_zeff",
    "traceless_tensor",
]
__version__ = "0.5.0"
