#!/usr/bin/env python3
"""Validate and draw a 16-function ligand s/p SALC basis for Td symmetry.

Run from any directory with ``python tetrahedral_salcs/generate_figures.py``.
Drawings are schematic signed-orbital diagrams, not density isosurfaces.
"""

from __future__ import annotations

import argparse
import itertools
import math
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
FIGURES = ROOT / "figures"
DERIVATION = ROOT / "derivation.md"
BLUE = "#175f91"
ORANGE = "#d46b36"
VERTICES = np.array(((1, 1, 1), (1, -1, -1),
                     (-1, 1, -1), (-1, -1, 1)), dtype=float)
LABELS = ("A", "B", "C", "D")
CHARACTERS = {
    # E, 8C3, 3C2, 6S4, 6 sigma_d
    "A1": (1, 1, 1, 1, 1),
    "A2": (1, 1, 1, -1, -1),
    "E": (2, -1, 2, 0, 0),
    "T1": (3, 0, -1, 1, -1),
    "T2": (3, 0, -1, -1, 1),
}


@dataclass(frozen=True)
class SALC:
    label: str
    description: str
    irrep: str
    family: str
    factor: str
    integer_coefficients: np.ndarray  # Four rows: s scalar or (px, py, pz)
    coefficients: np.ndarray  # Shape (4, 4): local (s, px, py, pz)


def _record(label: str, description: str, irrep: str, family: str,
            factor: str, integers: np.ndarray) -> SALC:
    scales = {"1/2": 1 / 2, "1/(2√3)": 1 / (2 * math.sqrt(3)),
              "1/(2√2)": 1 / (2 * math.sqrt(2)),
              "1/(2√6)": 1 / (2 * math.sqrt(6))}
    vectors = np.zeros((4, 4))
    vectors[:, 0 if family == "s" else slice(1, 4)] = integers * scales[factor]
    return SALC(label, description, irrep, family, factor, integers, vectors)


def all_salcs() -> list[SALC]:
    x, y, z = VERTICES.T.astype(int)
    ones = np.ones(4, dtype=int)

    def s(label: str, description: str, irrep: str, numbers: np.ndarray) -> SALC:
        return _record(label, description, irrep, "s", "1/2", numbers)

    def p(label: str, description: str, irrep: str, family: str,
          factor: str, *columns: np.ndarray) -> SALC:
        return _record(label, description, irrep, family, factor,
                       np.column_stack(columns))

    return [
        s("Ψ(s;0)", "uniform ligand s", "A1", ones),
        s("Ψ(s;x)", "x-weighted ligand s", "T2", x),
        s("Ψ(s;y)", "y-weighted ligand s", "T2", y),
        s("Ψ(s;z)", "z-weighted ligand s", "T2", z),
        p("Ψ(σ;0)", "uniform radial pσ", "A1", "radial pσ",
          "1/(2√3)", x, y, z),
        p("Ψ(σ;x)", "x-weighted radial pσ", "T2", "radial pσ",
          "1/(2√3)", ones, z, y),
        p("Ψ(σ;y)", "y-weighted radial pσ", "T2", "radial pσ",
          "1/(2√3)", z, ones, x),
        p("Ψ(σ;z)", "z-weighted radial pσ", "T2", "radial pσ",
          "1/(2√3)", y, x, ones),
        p("Ψ(E;θ)", "tangential E, θ partner", "E", "tangential E",
          "1/(2√6)", -x, -y, 2 * z),
        p("Ψ(E;ε)", "tangential E, ε partner", "E", "tangential E",
          "1/(2√2)", x, -y, np.zeros(4, dtype=int)),
        p("Ψ(T₁;x)", "tangential T₁, x partner", "T1", "tangential T₁",
          "1/(2√2)", np.zeros(4, dtype=int), -z, y),
        p("Ψ(T₁;y)", "tangential T₁, y partner", "T1", "tangential T₁",
          "1/(2√2)", z, np.zeros(4, dtype=int), -x),
        p("Ψ(T₁;z)", "tangential T₁, z partner", "T1", "tangential T₁",
          "1/(2√2)", -y, x, np.zeros(4, dtype=int)),
        p("Ψ(πT₂;x)", "tangential T₂, x partner", "T2", "tangential T₂",
          "1/(2√6)", 2 * ones, -z, -y),
        p("Ψ(πT₂;y)", "tangential T₂, y partner", "T2", "tangential T₂",
          "1/(2√6)", -z, 2 * ones, -x),
        p("Ψ(πT₂;z)", "tangential T₂, z partner", "T2", "tangential T₂",
          "1/(2√6)", -y, -x, 2 * ones),
    ]


def group_matrices() -> list[np.ndarray]:
    """All signed permutation matrices preserving xyz = +1."""
    matrices = []
    for permutation in itertools.permutations(range(3)):
        for signs in itertools.product((-1, 1), repeat=3):
            if math.prod(signs) != 1:
                continue
            operation = np.zeros((3, 3), dtype=int)
            operation[np.arange(3), permutation] = signs
            matrices.append(operation)
    return matrices


def representation(operation: np.ndarray) -> np.ndarray:
    """Act on four site s orbitals and four sets of Cartesian p orbitals."""
    matrix = np.zeros((16, 16))
    for j, vertex in enumerate(VERTICES):
        target = next(k for k, other in enumerate(VERTICES)
                      if np.array_equal(operation @ vertex, other))
        matrix[4 * target, 4 * j] = 1
        matrix[4 * target + 1:4 * target + 4,
               4 * j + 1:4 * j + 4] = operation
    return matrix


def class_representatives() -> tuple[np.ndarray, ...]:
    identity = np.eye(3, dtype=int)
    c3 = np.array(((0, 1, 0), (0, 0, 1), (1, 0, 0)))
    c2 = np.diag((1, -1, -1))
    s4 = np.array(((0, -1, 0), (1, 0, 0), (0, 0, -1)))
    mirror = np.array(((0, 1, 0), (1, 0, 0), (0, 0, 1)))
    return identity, c3, c2, s4, mirror


def class_name(operation: np.ndarray) -> str:
    if np.array_equal(operation, np.eye(3)):
        return "E"
    determinant = round(np.linalg.det(operation))
    if determinant == 1:
        return "C3" if np.array_equal(np.linalg.matrix_power(operation, 3),
                                       np.eye(3)) else "C2"
    return "S4" if np.array_equal(np.linalg.matrix_power(operation, 4),
                                   np.eye(3)) and not np.array_equal(
                                       np.linalg.matrix_power(operation, 2),
                                       np.eye(3)) else "mirror"


def validate() -> None:
    """Check the group, all invariant irrep spaces, and published coefficients."""
    matrices = group_matrices()
    assert len(matrices) == 24
    assert len({tuple(g.flat) for g in matrices}) == 24
    assert all(np.array_equal(g.T @ g, np.eye(3)) for g in matrices)
    assert all(abs(round(np.linalg.det(g))) == 1 for g in matrices)
    assert all(all(any(np.array_equal(g @ v, w) for w in VERTICES)
                   for v in VERTICES) for g in matrices)
    assert all(any(np.array_equal(a @ b, c) for c in matrices)
               for a in matrices for b in matrices)
    assert all(any(np.array_equal(a @ b, np.eye(3)) for b in matrices)
               for a in matrices)
    assert {name: sum(class_name(g) == name for g in matrices)
            for name in ("E", "C3", "C2", "S4", "mirror")} == {
                "E": 1, "C3": 8, "C2": 3, "S4": 6, "mirror": 6}
    assert all(np.allclose(representation(a @ b),
                           representation(a) @ representation(b))
               for a in matrices for b in matrices)

    salcs = all_salcs()
    assert len(salcs) == 16
    basis = np.vstack([salc.coefficients.reshape(16) for salc in salcs])
    assert np.allclose(basis @ basis.T, np.eye(16))
    assert all(np.isfinite(salc.coefficients).all() for salc in salcs)
    for salc in salcs[8:]:
        assert np.allclose(np.sum(salc.coefficients[:, 1:] * VERTICES,
                                  axis=1), 0.0)

    blocks = ((0,), (1, 2, 3), (4,), (5, 6, 7), (8, 9),
              (10, 11, 12), (13, 14, 15))
    for block in blocks:
        selected = basis[list(block)]
        others = basis[[i for i in range(16) if i not in block]]
        for operation in matrices:
            action = representation(operation)
            assert np.allclose(selected @ action @ others.T, 0.0, atol=1e-12)
        for class_index, representative in enumerate(class_representatives()):
            character = np.trace(selected @ representation(representative)
                                 @ selected.T)
            assert math.isclose(character, CHARACTERS[salcs[block[0]].irrep][
                class_index], abs_tol=1e-12)

    validate_table(salcs)


def validate_table(salcs: list[SALC]) -> None:
    """Compare every printed coefficient with the actual plotted basis."""
    lines = DERIVATION.read_text(encoding="utf-8").splitlines()
    rows = [line for line in lines if line.startswith("| ")
            and line[2:4].isdigit()]
    assert len(rows) == len(salcs) == 16
    for index, (line, salc) in enumerate(zip(rows, salcs), 1):
        cells = [cell.strip() for cell in line.strip("|").split("|")]
        assert len(cells) == 9
        assert int(cells[0]) == index
        assert cells[2] == salc.label
        assert cells[3] == salc.irrep.translate(str.maketrans("12", "₁₂"))
        assert cells[4] == salc.factor
        parsed = []
        for item in cells[5:]:
            if salc.family == "s":
                parsed.append(int(item.replace("−", "-")))
            else:
                assert item.startswith("(") and item.endswith(")")
                parsed.append(tuple(int(x.replace("−", "-"))
                                    for x in item[1:-1].split(",")))
        assert np.array_equal(np.asarray(parsed), salc.integer_coefficients)


def _sphere(ax, center: np.ndarray, radius: float, color: str) -> None:
    polar = np.linspace(0, math.pi, 12)
    azimuth = np.linspace(0, 2 * math.pi, 18)
    x = center[0] + radius * np.outer(np.sin(polar), np.cos(azimuth))
    y = center[1] + radius * np.outer(np.sin(polar), np.sin(azimuth))
    z = center[2] + radius * np.outer(np.cos(polar), np.ones_like(azimuth))
    ax.plot_surface(x, y, z, color=color, linewidth=0, shade=True,
                    antialiased=False, alpha=0.96)


def _lobe(ax, center: np.ndarray, axis: np.ndarray, size: float,
          color: str) -> None:
    reference = np.array((0.0, 0.0, 1.0))
    if abs(axis @ reference) > 0.93:
        reference = np.array((0.0, 1.0, 0.0))
    side1 = np.cross(axis, reference)
    side1 /= np.linalg.norm(side1)
    side2 = np.cross(axis, side1)
    polar = np.linspace(0, math.pi, 10)
    azimuth = np.linspace(0, 2 * math.pi, 16)
    along = 0.23 * size * np.outer(np.cos(polar), np.ones_like(azimuth))
    across1 = 0.13 * size * np.outer(np.sin(polar), np.cos(azimuth))
    across2 = 0.13 * size * np.outer(np.sin(polar), np.sin(azimuth))
    surface = (center[:, None, None] + axis[:, None, None] * along
               + side1[:, None, None] * across1
               + side2[:, None, None] * across2)
    ax.plot_surface(*surface, color=color, linewidth=0, shade=True,
                    antialiased=False, alpha=0.94)


def draw(salc: SALC, index: int) -> Path:
    fig = plt.figure(figsize=(7.5, 7.0), dpi=145)
    grid = fig.add_gridspec(2, 1, height_ratios=(4.6, 1.5), hspace=0.05)
    ax = fig.add_subplot(grid[0], projection="3d")
    lower = fig.add_subplot(grid[1])

    for i, j in itertools.combinations(range(4), 2):
        ax.plot(*VERTICES[[i, j]].T, color="#bac3ca", linewidth=0.9)
    ax.scatter(*VERTICES.T, color="#87929c", s=12, depthshade=False)

    peak = max(np.linalg.norm(salc.coefficients[:, 1:], axis=1))
    if salc.family == "s":
        peak = max(abs(salc.coefficients[:, 0]))
    for j, vertex in enumerate(VERTICES):
        if salc.family == "s":
            value = salc.coefficients[j, 0]
            if abs(value) > 1e-12:
                _sphere(ax, vertex, 0.25 * math.sqrt(abs(value) / peak),
                        BLUE if value > 0 else ORANGE)
        else:
            vector = salc.coefficients[j, 1:]
            amplitude = np.linalg.norm(vector)
            if amplitude > 1e-12:
                direction = vector / amplitude
                size = math.sqrt(amplitude / peak)
                for phase in (+1, -1):
                    center = vertex + phase * 0.31 * direction
                    _lobe(ax, center, direction, size,
                          BLUE if phase > 0 else ORANGE)
        ax.text(*(vertex + np.array((0.0, 0.0, 0.25))), LABELS[j],
                fontsize=10, color="#25313b", ha="center")
    ax.set_xlim(-1.65, 1.65)
    ax.set_ylim(-1.65, 1.65)
    ax.set_zlim(-1.65, 1.65)
    ax.set_box_aspect((1, 1, 1), zoom=1.55)
    ax.view_init(elev=19, azim=-56)
    ax.set_axis_off()

    if salc.family == "s":
        vals = salc.coefficients[:, 0]
        lower.bar(np.arange(4), vals,
                  color=[BLUE if v >= 0 else ORANGE for v in vals],
                  width=0.65)
        lower.axhline(0, color="#59656e", linewidth=0.8)
        lower.set_ylim(-0.6, 0.6)
        lower.set_xticks(np.arange(4), LABELS)
        lower.set_ylabel("s coefficient", fontsize=9)
        lower.spines[["top", "right", "left"]].set_visible(False)
        lower.tick_params(axis="both", labelsize=9, length=0)
        lower.set_yticks((-0.5, 0, 0.5))
        lower.grid(axis="y", color="#e3e8eb", linewidth=0.5)
    else:
        values = salc.coefficients[:, 1:].T
        lower.imshow(values, cmap="RdBu", vmin=-0.62, vmax=0.62,
                     aspect="auto", origin="upper")
        lower.set_xticks(np.arange(4), LABELS)
        lower.set_yticks(np.arange(3), ("p(x)", "p(y)", "p(z)"))
        lower.tick_params(axis="both", labelsize=9, length=0)
        for row in range(3):
            for column in range(4):
                value = values[row, column]
                lower.text(column, row, f"{value:+.3f}" if value else "0",
                           ha="center", va="center", fontsize=8,
                           color="#19242c")
        for spine in lower.spines.values():
            spine.set_visible(False)
    lower.set_title("Cartesian AO coefficients at each ligand",
                    fontsize=9, pad=6)

    irrep = salc.irrep.translate(str.maketrans("12", "₁₂"))
    fig.suptitle(f"{index:02d}  {salc.description}  ·  {irrep}",
                 fontsize=13, y=0.96)
    fig.text(0.5, 0.02, "Blue: positive phase    Orange: negative phase"
             "   ·   Schematic lobes", ha="center", fontsize=8,
             color="#4a5660")

    FIGURES.mkdir(parents=True, exist_ok=True)
    safe = (salc.label.replace("Ψ(", "").replace(")", "")
            .replace(";", "_").replace("₁", "1").replace("₂", "2")
            .replace("σ", "sigma").replace("π", "pi")
            .replace("θ", "theta").replace("ε", "epsilon"))
    path = FIGURES / f"{index:02d}_{safe}.png"
    fig.savefig(path, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check-only", action="store_true",
                        help="validate the SALCs without drawing PNGs")
    args = parser.parse_args()
    validate()
    if args.check_only:
        print("Validated 16 Td SALCs, coefficients, and group characters.")
        return
    for index, salc in enumerate(all_salcs(), 1):
        draw(salc, index)
    print(f"Validated and generated 16 images in {FIGURES}")


if __name__ == "__main__":
    main()
