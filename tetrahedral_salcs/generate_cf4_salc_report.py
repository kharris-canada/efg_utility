#!/usr/bin/env python3
"""Generate the SALC-resolved, closed-shell CF4 teaching report."""

from __future__ import annotations

import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))

from efg_utility import (  # noqa: E402
    HydrogenicOrbital,
    OrbitalCombination,
    OrbitalCombinationTerm,
    calculate_orbital_combinations,
    efg_parameters_from_tensor,
    radius_from_zeff,
)
from efg_utility.constants import ANGSTROM_IN_METRES  # noqa: E402
from salc_to_efg import (  # noqa: E402
    _evaluate_expression,
    _parse_salc_id,
    read_table,
    select_salcs,
)


DISTANCE_ANGSTROM = 1.32
ZEFF = 6.0
ELECTRONS_PER_MO = 2
VALENCE_ELECTRONS = 32
ROUND_OFF_RELATIVE_TOLERANCE = 1.0e-12
TABLE = ROOT / "tetrahedral_salcs" / "salcs_cartesian.tsv"
REPORT = ROOT / "tetrahedral_salcs" / "cf4_salc_efg_report.md"


def cf_bond_frame() -> np.ndarray:
    """Columns are new x,y,z axes; the (+,+,+) C-F bond lies along new z."""
    x_axis = np.array((1.0, -1.0, 0.0)) / np.sqrt(2.0)
    z_axis = np.array((1.0, 1.0, 1.0)) / np.sqrt(3.0)
    y_axis = np.cross(z_axis, x_axis)
    return np.column_stack((x_axis, y_axis, z_axis))


def no_d_mo_model() -> dict[str, tuple[np.ndarray, np.ndarray, int]]:
    """Return eigenvalues, eigenvectors, and occupied count for each irrep."""
    a1 = np.array(
        ((-5.0, 0.0, 1.2), (0.0, -3.7, 1.5), (1.2, 1.5, -2.0)),
        dtype=np.float64,
    )
    t2 = np.diag((-5.0, -3.7, -2.0, 0.5)).astype(np.float64)
    t2[:3, 3] = (1.0, 0.8, 0.7)
    t2[3, :3] = (1.0, 0.8, 0.7)
    a1_values, a1_vectors = np.linalg.eigh(a1)
    t2_values, t2_vectors = np.linalg.eigh(t2)
    return {
        "a1": (a1_values, a1_vectors, 2),
        "t2": (t2_values, t2_vectors, 3),
    }


def salc_projection(irrep: str, rank: str, model: dict[str, tuple[np.ndarray, np.ndarray, int]]) -> tuple[float, float]:
    """Return (occupied SALC weight, dominant occupied MO energy)."""
    if irrep == "a1":
        row = {"A": 0, "B": 1}[rank]
        values, vectors, occupied = model["a1"]
    elif irrep == "t2":
        row = {"A": 0, "B": 1, "C": 2}[rank]
        values, vectors, occupied = model["t2"]
    elif irrep in {"e", "t1"}:
        return 1.0, -2.0
    else:
        raise ValueError(f"unexpected SALC irrep/rank {irrep}-{rank}")
    weights = vectors[row, :occupied] ** 2
    return float(np.sum(weights)), float(values[int(np.argmax(weights))])


def _format_matrix_cell(tensor: np.ndarray, threshold: float) -> str:
    """Render a matrix as plain text with Markdown line breaks inside a cell."""
    rows = []
    for row in tensor:
        values = [0.0 if abs(value) < threshold else value for value in row]
        rows.append("(" + ", ".join(f"{value:+.3e}" for value in values) + ")")
    return "[" + "<br>".join(rows) + "]"


def _format_matrix_text(tensor: np.ndarray, threshold: float) -> str:
    rows = []
    for row in tensor:
        values = [0.0 if abs(value) < threshold else value for value in row]
        rows.append("[" + "  ".join(f"{value:+.3e}" for value in values) + "]")
    return "\n".join(rows)


def _mo_label(irrep: str, energy: float, model: dict[str, tuple[np.ndarray, np.ndarray, int]]) -> str:
    if irrep == "a1":
        index = int(np.argmin(np.abs(model["a1"][0] - energy))) + 1
    elif irrep == "t2":
        index = int(np.argmin(np.abs(model["t2"][0] - energy))) + 1
    else:
        return "1e" if irrep == "e" else "1t₁"
    return f"{index}{'a₁' if irrep == 'a1' else 't₂'}"


def main() -> None:
    model = no_d_mo_model()
    table = read_table(TABLE)
    all_salcs = select_salcs(table, "all")

    levels = []
    for irrep, (energies, _vectors, _occupied_count) in model.items():
        degeneracy = 1 if irrep == "a1" else 3
        levels.extend((float(energy), irrep, degeneracy) for energy in energies)
    levels.extend(((-2.0, "e", 2), (-2.0, "t1", 3)))
    levels.sort(key=lambda item: (item[0], item[1]))
    orbitals_to_fill = VALENCE_ELECTRONS // ELECTRONS_PER_MO
    occupied = 0
    occupied_levels: set[tuple[str, float]] = set()
    for energy, irrep, degeneracy in levels:
        if occupied + degeneracy <= orbitals_to_fill:
            occupied += degeneracy
            occupied_levels.add((irrep, energy))
    if occupied != orbitals_to_fill:
        raise RuntimeError("toy MO scheme did not yield an integer closed-shell filling")

    radius_m = radius_from_zeff(ZEFF)
    calculated: dict[str, dict[str, object]] = {}
    group_keys: dict[tuple[str, str], tuple[str, float]] = {}
    for salc_id, rows in all_salcs.items():
        irrep, rank, _partner = _parse_salc_id(
            salc_id, path=TABLE, line_number=rows[0].line_number
        )
        weight, dominant_energy = salc_projection(irrep, rank, model)
        if (irrep, dominant_energy) not in occupied_levels:
            raise RuntimeError(f"{salc_id} was assigned to an unoccupied MO")
        terms = []
        for row in rows:
            position_angstrom = tuple(
                _evaluate_expression(expression, DISTANCE_ANGSTROM)
                for expression in row.coordinate_expressions
            )
            orbital = HydrogenicOrbital(
                center_m=tuple(value * ANGSTROM_IN_METRES for value in position_angstrom),
                orbital_type=row.ao,
                radius_m=radius_m,
            )
            terms.append(OrbitalCombinationTerm(orbital, row.coefficient))
        occupied_electrons = ELECTRONS_PER_MO * weight
        field = calculate_orbital_combinations((OrbitalCombination(
            terms=tuple(terms),
            occupation=occupied_electrons,
            name=salc_id,
        ),))
        calculated[salc_id] = {
            "irrep": irrep,
            "rank": rank,
            "weight": weight,
            "occupation": occupied_electrons,
            "energy": dominant_energy,
            "tensor": field.efg_tensor_v_per_m2,
        }
        group_keys[(irrep, rank)] = (irrep, dominant_energy)

    grouped: dict[tuple[str, str], list[str]] = defaultdict(list)
    for salc_id, item in calculated.items():
        grouped[(str(item["irrep"]), str(item["rank"]))].append(salc_id)
    for ids in grouped.values():
        ids.sort(key=lambda value: _parse_salc_id(
            value, path=TABLE, line_number=1
        )[2])

    total_tensor = np.sum(
        [np.asarray(item["tensor"]) for item in calculated.values()], axis=0
    )
    scale = sum(
        float(np.max(np.abs(np.asarray(item["tensor"]))))
        for item in calculated.values()
    )
    total_parameters = efg_parameters_from_tensor(total_tensor, cancellation_scale=scale)
    frame = cf_bond_frame()
    display_threshold = scale * ROUND_OFF_RELATIVE_TOLERANCE

    # Rotate tensors into a right-handed coordinate frame whose z axis is the
    # C-F bond to ligand (+,+,+). This is a coordinate change, not a PAS claim.
    for item in calculated.values():
        tensor = np.asarray(item["tensor"])
        item["tensor"] = frame.T @ tensor @ frame
    total_tensor = frame.T @ total_tensor @ frame

    energy_rows = []
    for energy, irrep, degeneracy in levels:
        label = _mo_label(irrep, energy, model)
        status = "occupied" if (irrep, energy) in occupied_levels else "virtual"
        energy_rows.append(
            f"| {label} ({irrep.upper()}) | {energy:+.3f} | {degeneracy} | {status} |"
        )

    sections = []
    order = sorted(
        grouped,
        key=lambda key: (group_keys[key][1], key[0], key[1]),
    )
    for irrep, rank in order:
        ids = grouped[(irrep, rank)]
        energy = float(calculated[ids[0]]["energy"])
        degeneracy = len(ids)
        mo_label = _mo_label(irrep, energy, model)
        tensors = [np.asarray(calculated[salc_id]["tensor"]) for salc_id in ids]
        row_sum = np.sum(tensors, axis=0)
        paired_tables = []
        for start in range(0, degeneracy, 2):
            pair_ids = ids[start:start + 2]
            pair_tensors = tensors[start:start + 2]
            headers = []
            occupations = []
            matrix_cells = []
            for salc_id, tensor in zip(pair_ids, pair_tensors, strict=True):
                headers.append(salc_id)
                occupations.append(
                    f"Projected occupation: {calculated[salc_id]['occupation']:.4f} e"
                )
                matrix_cells.append(_format_matrix_cell(tensor, display_threshold))
            table_header = "| " + " | ".join(headers) + " |"
            table_rule = "| " + " | ".join(":--" for _ in headers) + " |"
            occupation_row = "| " + " | ".join(occupations) + " |"
            matrix_row = "| " + " | ".join(matrix_cells) + " |"
            paired_tables.append("\n".join((
                table_header,
                table_rule,
                occupation_row,
                matrix_row,
            )))

        sum_table = "\n".join((
            "| Quantity | EFG tensor (V/m²) |",
            "| :-- | :-- |",
            "| Sum over " + str(degeneracy) + " partner(s) | "
            + _format_matrix_cell(row_sum, display_threshold) + " |",
        ))
        markdown_table = "\n\n".join((*paired_tables, sum_table))
        sections.append(
            f"### {mo_label} ({irrep.upper()}), energy {energy:+.3f} toy units; "
            f"{degeneracy}-fold\n\n{markdown_table}"
        )

    total_electrons_in_salc_projection = sum(
        float(item["occupation"]) for item in calculated.values()
    )
    eta = "undefined" if total_parameters.eta is None else f"{total_parameters.eta:.6g}"
    energy_table = "\n".join(energy_rows)
    section_text = "\n\n".join(sections)
    report = f"""# CF₄: SALC-resolved electric-field-gradient contributions

## Total from the occupied SALC projections

The summed traceless EFG tensor at the central carbon is

```text
{_format_matrix_text(total_tensor, display_threshold)}
```

`VZZ = {total_parameters.vzz:+.6e} V/m²`; `eta = {eta}`.

## Electron count and illustrative MO filling

Neutral CF₄ has 42 electrons. Removing the ten 1s core electrons leaves 32
valence electrons. The model contains the 16 ligand 2s/2p SALCs and central
carbon 2s/2p functions (20 spatial MOs total); no d functions are included.
Filling the lowest 16 spatial MOs with two electrons each gives the following
toy ordering. Only the central s/p basis functions needed to estimate MO
mixing and occupation projections enter this table; **their EFG fields are
not calculated or included**.

| MO level | Toy energy | Degeneracy | Filling |
|:--|--:|--:|:--|
{energy_table}

For A₁ and T₂, each SALC's projected occupation is twice the sum of its squared
coefficients in the occupied MOs of that symmetry. These are fractional
projections, not integer electrons assigned to isolated SALCs: the toy MOs
mix ligand SALCs with central-carbon s/p functions, so an occupied two-electron
MO can place only a fraction of its population in a particular ligand SALC.
Its row is placed at the
lowest-energy occupied MO where it has its largest SALC weight. E and T₁ are
nonbonding in this no-d model and therefore retain two electrons per partner.
This bookkeeping assigns {total_electrons_in_salc_projection:.3f} electrons to
the ligand SALC projections; the remainder of the 32-electron population is
on central-atom components and is intentionally omitted from the EFG sum.

The energy matrix is the explicitly illustrative model from the tetrahedral
bonding diagram with its d orbitals removed: ligand s, radial pσ, and
tangential pπ start at −5.0, −3.7, and −2.0 arbitrary units; carbon 2s and 2p
start at −2.0 and +0.5. The A₁ couplings are 1.2 and 1.5; the T₂ couplings
from ligand s/radial-p/tangential-p to carbon 2p are 1.0/0.8/0.7. These are
teaching values, not a fitted or spectroscopic CF₄ calculation. The geometry
uses `d = {DISTANCE_ANGSTROM:.2f} Å` (representative of the reported
C–F distance near 1.32 Å) and `Zeff = {ZEFF:.1f}` for all ligand hydrogenic
2s/2p functions. The latter is a modeling choice, not an experimentally
determined effective charge. [A CF₄ structure reference](https://atom.uni-frankfurt.de/publications/files/Heck2017Master.pdf)
reports 132.3 pm.

## SALC tensors, ordered by their dominant occupied MO

Each SALC field is calculated from its **complete signed AO combination**,
including AO interference and finite-basis overlap normalization, at the
projected occupation shown. The central-atom AO fields and cross terms between
different SALCs in the same MO are omitted by design. Degenerate SALC partners
are shown no more than two per row, and the sum of their full Cartesian tensors
appears in a separate table below them.
Tensor entries are in V/m². The reported coordinate frame is right-handed with
the C–F bond to the `(+, +, +)` ligand along z. Because the total Td EFG is
zero, this z axis is a convenient symmetry axis, not a uniquely defined PAS.
Displayed entries smaller than {ROUND_OFF_RELATIVE_TOLERANCE:.0e} times the
largest individual SALC tensor scale are set to zero; this suppresses numerical
roundoff, not physical contributions.

{section_text}

## Interpretation

The rows are SALC contributions arranged using a simple MO filling model, not
the EFG of the actual molecular orbitals. The per-SALC approximation omits
inter-SALC and central-atom interference, so it is not an exact decomposition
of a quantum-chemical density. Still, equal occupation of every partner in
each degenerate irrep makes the summed density tetrahedrally symmetric. Td
has no nonzero invariant, symmetric, traceless rank-2 tensor, so the EFG must
vanish in every coordinate system—not only in a PAS. Individual matrices
need not be diagonal in the bond-aligned frame, and the tensor assigned to an
individual member of a degenerate set depends on the chosen partner basis;
the sum over the full degenerate set is the symmetry-respecting quantity.
Before the display roundoff threshold is applied, the total's largest residual entry is
{np.max(np.abs(total_tensor)):.3e} V/m², roughly {np.max(np.abs(total_tensor)) / scale:.1e}
of the sum-of-contributions scale. This is numerical residue, not a tensor
that needs diagonalization into a PAS.
"""
    REPORT.write_text(report, encoding="utf-8")
    print(f"Wrote {REPORT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
