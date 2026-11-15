from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from efg_utility.errors import CalculationError
from efg_utility.input import InputFileError, read_orbital_input
from efg_utility.models import (
    HydrogenicOrbital,
    OrbitalCombination,
    OrbitalCombinationTerm,
)
from efg_utility.constants import BOHR_RADIUS_M
from efg_utility.orbital import (
    _cross_density_unit_fields,
    _orbital_identity,
    calculate_orbital_combinations,
)


def _ao(center: tuple[float, float, float]) -> HydrogenicOrbital:
    return HydrogenicOrbital(center, "1s", 0.8e-10)


def test_coherent_identical_terms_normalize_to_two_electrons() -> None:
    orbital = _ao((0.0, 0.0, 1.0e-10))
    combination = OrbitalCombination(
        terms=(
            OrbitalCombinationTerm(orbital, 0.5),
            OrbitalCombinationTerm(orbital, 0.5),
        ),
        occupation=2.0,
        name="same-AO",
    )
    expected = calculate_orbital_combinations((
        OrbitalCombination((OrbitalCombinationTerm(orbital, 1.0),), 2.0),
    ))
    actual = calculate_orbital_combinations((combination,))
    np.testing.assert_allclose(actual.potential_v, expected.potential_v, rtol=2e-13)
    np.testing.assert_allclose(
        actual.hessian_v_per_m2, expected.hessian_v_per_m2, rtol=2e-13, atol=1e-5
    )


def test_relative_ao_phase_changes_coherent_density_and_efg() -> None:
    first = _ao((0.8e-10, 0.0, 0.0))
    second = _ao((0.0, 0.8e-10, 0.0))
    constructive = OrbitalCombination(
        (OrbitalCombinationTerm(first, 1.0), OrbitalCombinationTerm(second, 1.0)),
        occupation=2.0,
        name="constructive",
    )
    destructive = OrbitalCombination(
        (OrbitalCombinationTerm(first, 1.0), OrbitalCombinationTerm(second, -1.0)),
        occupation=2.0,
        name="destructive",
    )
    plus = calculate_orbital_combinations((constructive,))
    minus = calculate_orbital_combinations((destructive,))
    assert not np.allclose(plus.efg_tensor_v_per_m2, minus.efg_tensor_v_per_m2)
    assert plus.parameters.vzz != pytest.approx(minus.parameters.vzz, rel=1e-5)


def test_overlap_quadrature_matches_analytic_1s_overlap() -> None:
    first = HydrogenicOrbital((0.0, 0.0, 0.0), "1s", BOHR_RADIUS_M)
    second = HydrogenicOrbital((BOHR_RADIUS_M, 0.0, 0.0), "1s", BOHR_RADIUS_M)
    identities = sorted((_orbital_identity(first), _orbital_identity(second)), key=repr)
    overlap = _cross_density_unit_fields(*identities)[3]
    expected = np.exp(-1.0) * (1.0 + 1.0 + 1.0 / 3.0)
    assert overlap == pytest.approx(expected, rel=1e-7)


def test_grouped_ten_column_input_preserves_signed_coefficients(tmp_path: Path) -> None:
    input_file = tmp_path / "salc.txt"
    input_file.write_text(
        "A 0 0 1 1s - 4 +0.6 salc-A-1 2\n"
        "B 1 0 0 2pz - 4 -0.8 salc-A-1 2\n",
        encoding="utf-8",
    )
    independent, combinations = read_orbital_input(input_file)
    assert independent == []
    assert len(combinations) == 1
    assert combinations[0].name == "salc-A-1"
    assert combinations[0].occupation == 2.0
    assert [term.coefficient for term in combinations[0].terms] == [0.6, -0.8]
    result = calculate_orbital_combinations(combinations)
    assert np.all(np.isfinite(result.efg_tensor_v_per_m2))


def test_inconsistent_occupation_in_group_is_rejected(tmp_path: Path) -> None:
    input_file = tmp_path / "salc.txt"
    input_file.write_text(
        "A 0 0 1 1s - 4 0.6 salc-A-1 2\n"
        "B 1 0 0 2pz - 4 0.8 salc-A-1 1\n",
        encoding="utf-8",
    )
    with pytest.raises(InputFileError, match="occupation.*does not match"):
        read_orbital_input(input_file)


def test_zero_norm_coherent_combination_is_rejected() -> None:
    orbital = _ao((0.0, 0.0, 1.0e-10))
    combination = OrbitalCombination(
        (OrbitalCombinationTerm(orbital, 1.0), OrbitalCombinationTerm(orbital, -1.0)),
        occupation=2.0,
        name="cancelled",
    )
    with pytest.raises(CalculationError, match="nonpositive or negligible norm"):
        calculate_orbital_combinations((combination,))
