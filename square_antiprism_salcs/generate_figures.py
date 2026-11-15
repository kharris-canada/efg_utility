#!/usr/bin/env python3
"""Validate and draw the 32 s/p ligand SALCs of a D4d square antiprism.

Run from any directory with ``python square_antiprism_salcs/generate_figures.py``.
The PNGs are schematic signed-orbital drawings, not computed isodensity surfaces.
"""

from __future__ import annotations

import argparse
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
RADIUS = 1.0
HALF_HEIGHT = 2.0 ** (-0.75)  # equal top-square and lateral edge lengths
POSITIVE = "#175f91"
NEGATIVE = "#d46b36"
NEUTRAL = "#87929c"


@dataclass(frozen=True)
class Family:
    key: str
    title: str
    s_sign: int
    m_sign: int
    axis: str | None


@dataclass(frozen=True)
class SALC:
    family: Family
    mode: str
    m: int
    component: str | None
    irrep: str
    coefficients: np.ndarray


FAMILIES = (
    Family("s", "ligand s", +1, +1, None),
    Family("sigma", "bond-directed pσ", +1, +1, "sigma"),
    Family("meridional", "meridional pπ", -1, +1, "meridional"),
    Family("tangential", "tangential pπ", +1, -1, "tangential"),
)
MODES = ((0, None), (4, None), (1, "c"), (1, "s"),
         (2, "c"), (2, "s"), (3, "c"), (3, "s"))
CHARACTERS = {
    # E, S8, C4, S8^3, C2, C2', sigma_d
    "A1": (1, 1, 1, 1, 1, 1, 1),
    "A2": (1, 1, 1, 1, 1, -1, -1),
    "B1": (1, -1, 1, -1, 1, 1, -1),
    "B2": (1, -1, 1, -1, 1, -1, 1),
    "E1": (2, math.sqrt(2), 0, -math.sqrt(2), -2, 0, 0),
    "E2": (2, 0, -2, 0, 2, 0, 0),
    "E3": (2, -math.sqrt(2), 0, math.sqrt(2), -2, 0, 0),
}


def geometry() -> tuple[np.ndarray, dict[str, np.ndarray]]:
    angles = np.arange(8) * math.pi / 4.0
    radial = np.column_stack((np.cos(angles), np.sin(angles), np.zeros(8)))
    tangent = np.column_stack((-np.sin(angles), np.cos(angles), np.zeros(8)))
    z = np.where(np.arange(8) % 2 == 0, HALF_HEIGHT, -HALF_HEIGHT)
    ez = np.array((0.0, 0.0, 1.0))
    length = math.hypot(RADIUS, HALF_HEIGHT)
    sigma = (RADIUS * radial + z[:, None] * ez) / length
    meridional = (z[:, None] * radial - RADIUS * ez) / length
    points = RADIUS * radial + z[:, None] * ez
    return points, {
        "sigma": sigma,
        "meridional": meridional,
        "tangential": tangent,
    }


def coefficients(m: int, component: str | None) -> np.ndarray:
    j = np.arange(8)
    if m == 0:
        return np.ones(8) / math.sqrt(8)
    if m == 4:
        return (-1.0) ** j / math.sqrt(8)
    angles = m * j * math.pi / 4.0
    return np.cos(angles) / 2.0 if component == "c" else np.sin(angles) / 2.0


def irrep_for(family: Family, m: int) -> str:
    if m in (0, 4):
        s_eigenvalue = family.s_sign * (1 if m == 0 else -1)
        return {
            (+1, +1): "A1",
            (+1, -1): "A2",
            (-1, -1): "B1",
            (-1, +1): "B2",
        }[(s_eigenvalue, family.m_sign)]
    s_character = 2.0 * family.s_sign * math.cos(m * math.pi / 4.0)
    if s_character > 0.5:
        return "E1"
    if s_character < -0.5:
        return "E3"
    return "E2"


def all_salcs() -> list[SALC]:
    return [
        SALC(family, f"m{m}{component or ''}", m, component,
             irrep_for(family, m), coefficients(m, component))
        for family in FAMILIES
        for m, component in MODES
    ]


def operation_matrices() -> tuple[np.ndarray, np.ndarray]:
    c = 1.0 / math.sqrt(2.0)
    s8 = np.array(((c, -c, 0.0), (c, c, 0.0), (0.0, 0.0, -1.0)))
    mirror = np.diag((1.0, -1.0, 1.0))
    return s8, mirror


def ao_action(family: Family, generator: str) -> np.ndarray:
    matrix = np.zeros((8, 8))
    for j in range(8):
        if generator == "S":
            matrix[(j + 1) % 8, j] = family.s_sign
        else:
            matrix[(-j) % 8, j] = family.m_sign
    return matrix


def validate() -> None:
    """Check geometry, local p signs, orthonormality, and every irrep character."""
    points, axes = geometry()
    s8, mirror = operation_matrices()
    eye3 = np.eye(3)
    eye8 = np.eye(8)
    assert np.allclose(s8.T @ s8, eye3)
    assert np.allclose(mirror.T @ mirror, eye3)
    assert np.allclose(np.linalg.matrix_power(s8, 8), eye3)
    assert np.allclose(mirror @ s8 @ mirror, s8.T)
    group = [np.linalg.matrix_power(s8, k) for k in range(8)]
    group += [mirror @ np.linalg.matrix_power(s8, k) for k in range(8)]
    assert all(not np.allclose(a, b) for i, a in enumerate(group)
               for b in group[i + 1:])
    assert all(any(np.allclose(a @ b, c) for c in group)
               for a in group for b in group)

    for j in range(8):
        nxt = (j + 1) % 8
        reflected = (-j) % 8
        assert np.allclose(s8 @ points[j], points[nxt])
        assert np.allclose(mirror @ points[j], points[reflected])
        for family in FAMILIES[1:]:
            local = axes[family.axis]
            assert np.allclose(s8 @ local[j], family.s_sign * local[nxt])
            assert np.allclose(mirror @ local[j],
                               family.m_sign * local[reflected])

    salcs = all_salcs()
    assert len(salcs) == 32
    class_reps = (
        lambda s, m: eye8,
        lambda s, m: s,
        lambda s, m: s @ s,
        lambda s, m: np.linalg.matrix_power(s, 3),
        lambda s, m: np.linalg.matrix_power(s, 4),
        lambda s, m: m @ s,
        lambda s, m: m,
    )
    for family in FAMILIES:
        records = [record for record in salcs if record.family == family]
        basis = np.vstack([record.coefficients for record in records])
        assert np.allclose(basis @ basis.T, eye8)
        s = ao_action(family, "S")
        m = ao_action(family, "M")
        assert np.allclose(np.linalg.matrix_power(s, 8), eye8)
        assert np.allclose(m @ s @ m, s.T)
        blocks = ((0,), (1,), (2, 3), (4, 5), (6, 7))
        for block in blocks:
            record = records[block[0]]
            selected = basis[list(block)]
            for class_index, representative in enumerate(class_reps):
                transformed = representative(s, m)
                projected = selected @ transformed @ selected.T
                actual = np.trace(projected)
                expected = CHARACTERS[record.irrep][class_index]
                assert math.isclose(actual, expected, abs_tol=1e-12)
                for other in range(8):
                    if other not in block:
                        assert np.allclose(selected @ transformed
                                           @ basis[other], 0.0, atol=1e-12)

    counts: dict[str, int] = {}
    for record in salcs:
        counts[record.irrep] = counts.get(record.irrep, 0) + 1
    assert counts == {"A1": 3, "B2": 3, "E1": 8, "E2": 8,
                      "E3": 8, "A2": 1, "B1": 1}
    validate_coefficient_table(salcs)


def validate_coefficient_table(salcs: list[SALC]) -> None:
    """Ensure the published 32-row table matches the plotted coefficients."""
    lines = DERIVATION.read_text(encoding="utf-8").splitlines()
    rows = [line for line in lines if line.startswith("| ")
            and line[2:4].isdigit()]
    assert len(rows) == len(salcs) == 32
    family_names = ("s", "pσ", "pπ meridional", "pπ tangential")
    for index, (line, record) in enumerate(zip(rows, salcs), 1):
        cells = [cell.strip() for cell in line.strip("|").split("|")]
        assert len(cells) == 13
        assert int(cells[0]) == index
        assert cells[1] == family_names[(index - 1) // 8]
        assert cells[3] == record.irrep.translate(
            str.maketrans("123", "₁₂₃"))
        factor = 1.0 / math.sqrt(8) if cells[4] == "1/√8" else 0.5
        assert cells[4] in ("1/√8", "1/2")
        values = {"1": 1.0, "−1": -1.0, "0": 0.0,
                  "t": 1.0 / math.sqrt(2), "−t": -1.0 / math.sqrt(2)}
        actual = np.array([values[cell] for cell in cells[5:]]) * factor
        assert np.allclose(actual, record.coefficients, atol=1e-12)


def _sphere(ax, center: np.ndarray, radius: float, color: str) -> None:
    polar = np.linspace(0.0, math.pi, 12)
    azimuth = np.linspace(0.0, 2.0 * math.pi, 18)
    x = center[0] + radius * np.outer(np.sin(polar), np.cos(azimuth))
    y = center[1] + radius * np.outer(np.sin(polar), np.sin(azimuth))
    z = center[2] + radius * np.outer(np.cos(polar), np.ones_like(azimuth))
    ax.plot_surface(x, y, z, color=color, linewidth=0, shade=True,
                    antialiased=False, alpha=0.95)


def _lobe(ax, center: np.ndarray, axis: np.ndarray, scale: float,
          color: str) -> None:
    reference = np.array((0.0, 0.0, 1.0))
    if abs(axis @ reference) > 0.93:
        reference = np.array((0.0, 1.0, 0.0))
    transverse1 = np.cross(axis, reference)
    transverse1 /= np.linalg.norm(transverse1)
    transverse2 = np.cross(axis, transverse1)
    polar = np.linspace(0.0, math.pi, 10)
    azimuth = np.linspace(0.0, 2.0 * math.pi, 16)
    along = 0.20 * scale * np.outer(np.cos(polar), np.ones_like(azimuth))
    across1 = 0.11 * scale * np.outer(np.sin(polar), np.cos(azimuth))
    across2 = 0.11 * scale * np.outer(np.sin(polar), np.sin(azimuth))
    surface = (center[:, None, None] + axis[:, None, None] * along
               + transverse1[:, None, None] * across1
               + transverse2[:, None, None] * across2)
    ax.plot_surface(*surface, color=color, linewidth=0, shade=True,
                    antialiased=False, alpha=0.93)


def draw(record: SALC, index: int) -> Path:
    points, axes = geometry()
    fig = plt.figure(figsize=(8.0, 7.0), dpi=145)
    grid = fig.add_gridspec(2, 1, height_ratios=(5.1, 1.2), hspace=0.08)
    ax = fig.add_subplot(grid[0], projection="3d")
    bars = fig.add_subplot(grid[1])

    top = (0, 2, 4, 6, 0)
    bottom = (1, 3, 5, 7, 1)
    for ring in (top, bottom):
        ax.plot(*points[list(ring)].T, color="#b6bec5", linewidth=1.0)
    for j in (0, 2, 4, 6):
        for k in ((j - 1) % 8, (j + 1) % 8):
            ax.plot(*points[[j, k]].T, color="#cbd1d6", linewidth=0.8)
    ax.scatter(*points.T, s=11, color=NEUTRAL, depthshade=False)

    peak = max(abs(record.coefficients))
    for j, coefficient in enumerate(record.coefficients):
        if abs(coefficient) < 1e-10:
            continue
        amplitude = abs(coefficient) / peak
        if record.family.axis is None:
            color = POSITIVE if coefficient > 0 else NEGATIVE
            _sphere(ax, points[j], 0.17 * math.sqrt(amplitude), color)
        else:
            direction = axes[record.family.axis][j]
            for sign in (+1, -1):
                lobe_center = points[j] + sign * 0.22 * direction
                color = POSITIVE if coefficient * sign > 0 else NEGATIVE
                _lobe(ax, lobe_center, direction, math.sqrt(amplitude), color)

    for j, point in enumerate(points):
        label = f"{'T' if j % 2 == 0 else 'B'}{j // 2}"
        ax.text(*(point + np.array((0.0, 0.0, 0.18))), label,
                fontsize=8, color="#25313b", ha="center")
    ax.set_xlim(-1.55, 1.55)
    ax.set_ylim(-1.55, 1.55)
    ax.set_zlim(-1.25, 1.25)
    ax.set_box_aspect((1, 1, 0.8), zoom=1.62)
    ax.view_init(elev=22, azim=-57)
    ax.set_axis_off()

    labels = [f"{'T' if j % 2 == 0 else 'B'}{j // 2}" for j in range(8)]
    colors = [POSITIVE if value >= 0 else NEGATIVE
              for value in record.coefficients]
    bars.bar(np.arange(8), record.coefficients, color=colors, width=0.65)
    bars.axhline(0, color="#59656e", linewidth=0.8)
    bars.set_ylim(-0.59, 0.59)
    bars.set_xticks(np.arange(8), labels)
    bars.set_ylabel("coefficient", fontsize=9)
    bars.tick_params(axis="both", labelsize=8, length=0)
    bars.spines[["top", "right", "left"]].set_visible(False)
    bars.spines["bottom"].set_color("#aeb8bf")
    bars.set_yticks((-0.5, 0, 0.5))
    bars.grid(axis="y", color="#e3e8eb", linewidth=0.5)
    bars.set_axisbelow(True)
    irrep_label = record.irrep.translate(str.maketrans("123", "₁₂₃"))
    mode_label = (f"m = {record.m}" if record.component is None else
                  f"m = {record.m}, {'cosine' if record.component == 'c' else 'sine'}")
    fig.suptitle(f"{index:02d}  {record.family.title}  ·  "
                 f"{irrep_label}  ·  {mode_label}", fontsize=13, y=0.96)
    fig.text(0.5, 0.02, "Blue: positive phase    Orange: negative phase  "
             "·  Schematic lobes; coefficient scale shown below",
             ha="center", fontsize=8, color="#4a5660")

    FIGURES.mkdir(parents=True, exist_ok=True)
    path = FIGURES / f"{index:02d}_{record.family.key}_{record.irrep}_{record.mode}.png"
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
        print("Validated 32 SALCs, their D4d characters, and the geometry.")
        return
    for index, record in enumerate(all_salcs(), 1):
        draw(record, index)
    print(f"Validated and generated 32 images in {FIGURES}")


if __name__ == "__main__":
    main()
