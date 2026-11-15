"""Physical constants and unit conversions used by :mod:`efg_utility`."""

from __future__ import annotations

import math

# Exact by definition in the SI.
ELEMENTARY_CHARGE_C = 1.602176634e-19

# 2022 CODATA recommended values.
VACUUM_PERMITTIVITY_F_PER_M = 8.8541878188e-12
BOHR_RADIUS_M = 5.29177210544e-11

COULOMB_CONSTANT = 1.0 / (4.0 * math.pi * VACUUM_PERMITTIVITY_F_PER_M)
ANGSTROM_IN_METRES = 1.0e-10
