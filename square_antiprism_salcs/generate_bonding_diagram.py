#!/usr/bin/env python3
"""Draw an illustrative D4d ligand-SALC / central-s,p,d MO diagram.

Energies and couplings are teaching choices, not a calculation for a compound.
Each plotted MO bar represents an entire symmetry-degenerate multiplet.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
FIGURE_PNG = ROOT / "bonding_diagram.png"
FIGURE_SVG = ROOT / "bonding_diagram.svg"

LIGAND_COLOR = "#087b85"
MO_COLOR = "#273f78"
CENTRAL_COLOR = "#ad6535"
LINK_COLOR = "#728393"
ENERGY_LIMITS = (-6.15, 2.03)


@dataclass(frozen=True)
class Sector:
    name: str
    label: str
    dimension: int
    ligand_labels: tuple[str, ...]
    ligand_energies: tuple[float, ...]
    central_labels: tuple[str, ...]
    central_energies: tuple[float, ...]
    couplings: tuple[tuple[float, ...], ...]

    def hamiltonian(self) -> np.ndarray:
        left = np.asarray(self.ligand_energies)
        right = np.asarray(self.central_energies)
        coupling = np.asarray(self.couplings, dtype=float).reshape(
            (len(left), len(right)))
        matrix = np.diag(np.concatenate((left, right)))
        matrix[:len(left), len(left):] = coupling
        matrix[len(left):, :len(left)] = coupling.T
        return matrix

    def solutions(self) -> tuple[np.ndarray, np.ndarray]:
        return np.linalg.eigh(self.hamiltonian())


SECTORS = (
    Sector(
        name="A1", label="A₁ · singlets", dimension=1,
        ligand_labels=("s SALC · A₁", "radial pσ · A₁", "meridional pπ · A₁"),
        ligand_energies=(-5.0, -3.7, -2.6),
        central_labels=("s · A₁", "d(z²) · A₁"),
        central_energies=(-2.0, -1.0),
        couplings=((1.1, 0.25), (1.3, 0.4), (0.5, 0.65)),
    ),
    Sector(
        name="B2", label="B₂ · singlets", dimension=1,
        ligand_labels=("s SALC · B₂", "radial pσ · B₂", "meridional pπ · B₂"),
        ligand_energies=(-5.0, -3.7, -2.6),
        central_labels=("p(z) · B₂",), central_energies=(0.5,),
        couplings=((0.65,), (0.9,), (0.6,)),
    ),
    Sector(
        name="A2", label="A₂ · singlet", dimension=1,
        ligand_labels=("tangential pπ · A₂",), ligand_energies=(-2.0,),
        central_labels=(), central_energies=(), couplings=((),),
    ),
    Sector(
        name="B1", label="B₁ · singlet", dimension=1,
        ligand_labels=("tangential pπ · B₁",), ligand_energies=(-2.0,),
        central_labels=(), central_energies=(), couplings=((),),
    ),
    Sector(
        name="E1", label="E₁ · doublets", dimension=2,
        ligand_labels=("s SALC · E₁", "radial pσ · E₁", "meridional pπ · E₁",
                       "tangential pπ · E₁"),
        ligand_energies=(-5.0, -3.7, -2.6, -2.0),
        central_labels=("p(x,y) · E₁",), central_energies=(0.5,),
        couplings=((0.75,), (0.95,), (0.55,), (0.4,)),
    ),
    Sector(
        name="E2", label="E₂ · doublets", dimension=2,
        ligand_labels=("s SALC · E₂", "radial pσ · E₂", "meridional pπ · E₂",
                       "tangential pπ · E₂"),
        ligand_energies=(-5.0, -3.7, -2.6, -2.0),
        central_labels=("d(x²−y²,xy) · E₂",), central_energies=(-1.0,),
        couplings=((0.25,), (0.6,), (0.8,), (0.7,)),
    ),
    Sector(
        name="E3", label="E₃ · doublets", dimension=2,
        ligand_labels=("s SALC · E₃", "radial pσ · E₃", "meridional pπ · E₃",
                       "tangential pπ · E₃"),
        ligand_energies=(-5.0, -3.7, -2.6, -2.0),
        central_labels=("d(xz,yz) · E₃",), central_energies=(-1.0,),
        couplings=((0.25,), (0.65,), (0.85,), (0.7,)),
    ),
)


def connector_endpoints(source_index: int, ligand_count: int,
                        source_energy: float, mo_energy: float,
                        left_x: float, mo_x: float, right_x: float,
                        half_width: float) -> tuple[tuple[float, float],
                                                    tuple[float, float]]:
    """Return the (x, energy) endpoints in left-to-right plot order."""
    if source_index < ligand_count:
        return ((left_x + half_width, source_energy),
                (mo_x - half_width, mo_energy))
    return ((mo_x + half_width, mo_energy),
            (right_x - half_width, source_energy))


def validate() -> None:
    """Check symmetry counts, levels, eigensystems, and connector direction."""
    sectors = {sector.name: sector for sector in SECTORS}
    assert set(sectors) == {"A1", "A2", "B1", "B2", "E1", "E2", "E3"}
    ligand_orbitals = sum(s.dimension * len(s.ligand_labels) for s in SECTORS)
    central_orbitals = sum(s.dimension * len(s.central_labels) for s in SECTORS)
    mo_orbitals = sum(s.dimension * len(s.hamiltonian()) for s in SECTORS)
    assert (ligand_orbitals, central_orbitals, mo_orbitals) == (32, 9, 41)
    assert {name: len(s.hamiltonian()) for name, s in sectors.items()} == {
        "A1": 5, "A2": 1, "B1": 1, "B2": 4,
        "E1": 5, "E2": 5, "E3": 5}
    for sector in SECTORS:
        hamiltonian = sector.hamiltonian()
        values, vectors = sector.solutions()
        assert np.allclose(hamiltonian, hamiltonian.T)
        assert np.isfinite(hamiltonian).all()
        assert np.allclose(vectors.T @ vectors, np.eye(len(values)))
        assert np.allclose(hamiltonian @ vectors, vectors * values)
        assert np.isclose(sum(values), np.trace(hamiltonian))
        assert ENERGY_LIMITS[0] < min(values) <= max(values) < ENERGY_LIMITS[1]
    assert np.isclose(sectors["A2"].solutions()[0][0], -2.0)
    assert np.isclose(sectors["B1"].solutions()[0][0], -2.0)
    ligand_link = connector_endpoints(0, 3, -5.0, -5.5,
                                      0.43, 1.53, 2.63, 0.215)
    central_link = connector_endpoints(3, 3, -2.0, -0.5,
                                       0.43, 1.53, 2.63, 0.215)
    assert ligand_link[0][1] == -5.0 and ligand_link[1][1] == -5.5
    assert central_link[0][1] == -0.5 and central_link[1][1] == -2.0


def _draw_level(ax, center_x: float, energy: float, label: str,
                color: str, *, width: float = 0.43) -> None:
    ax.hlines(energy, center_x - width / 2, center_x + width / 2,
              linewidth=3.0, color=color, zorder=3)
    ax.text(center_x, energy + 0.13, label, ha="center", va="bottom",
            fontsize=8.4, color="#263442", zorder=4)


def draw_sector(ax, sector: Sector, *, show_ylabel: bool) -> None:
    eigenvalues, vectors = sector.solutions()
    ligand_count = len(sector.ligand_labels)
    inputs = np.asarray(sector.ligand_energies + sector.central_energies)
    left_x, mo_x, right_x = 0.43, 1.53, 2.63
    half_width = 0.215

    ax.set_xlim(0.0, 3.05)
    ax.set_ylim(*ENERGY_LIMITS)
    ax.set_xticks(())
    ax.set_yticks((-6, -4, -2, 0, 2))
    ax.grid(axis="y", color="#e5e9ec", linewidth=0.65, zorder=0)
    ax.set_axisbelow(True)
    ax.set_title(sector.label, fontsize=13, loc="left", pad=8,
                 fontweight="medium", color="#1e2d38")
    ax.text(left_x, 1.76, "LIGAND SALCs", ha="center", fontsize=8.6,
            color=LIGAND_COLOR, fontweight="medium")
    ax.text(mo_x, 1.76, "MO LEVELS", ha="center", fontsize=8.6,
            color=MO_COLOR, fontweight="medium")
    ax.text(right_x, 1.76, "CENTRAL ATOM", ha="center", fontsize=8.6,
            color=CENTRAL_COLOR, fontweight="medium")
    if show_ylabel:
        ax.set_ylabel("Illustrative energy (arbitrary units)", fontsize=9)
    else:
        ax.tick_params(axis="y", labelleft=False)
    for spine in ("top", "right", "bottom"):
        ax.spines[spine].set_visible(False)
    ax.spines["left"].set_color("#b3bcc4")
    ax.tick_params(axis="y", labelsize=8, length=0)

    # Lines represent components of at least 8% weight in the eigenvector.
    for i, source_energy in enumerate(inputs):
        for j, target_energy in enumerate(eigenvalues):
            weight = vectors[i, j] ** 2
            if weight < 0.08:
                continue
            (x0, y0), (x1, y1) = connector_endpoints(
                i, ligand_count, source_energy, target_energy,
                left_x, mo_x, right_x, half_width)
            ax.plot((x0, x1), (y0, y1),
                    color=LINK_COLOR, alpha=0.22 + 0.43 * weight,
                    linewidth=0.65 + 1.65 * weight, zorder=1)

    for label, energy in zip(sector.ligand_labels, sector.ligand_energies):
        _draw_level(ax, left_x, energy, label, LIGAND_COLOR)
    for label, energy in zip(sector.central_labels, sector.central_energies):
        _draw_level(ax, right_x, energy, label, CENTRAL_COLOR)
    irrep = sector.name.translate(str.maketrans("123", "₁₂₃")).lower()
    for j, energy in enumerate(eigenvalues, 1):
        _draw_level(ax, mo_x, energy,
                    f"{j}{irrep} (×{sector.dimension})", MO_COLOR)

    if not sector.central_labels:
        ax.text(right_x, -2.02,
                f"No central s/p/d\norbital of {sector.name.translate(str.maketrans('12', '₁₂'))} symmetry",
                ha="center", va="center", fontsize=9, color="#66737e")


def draw() -> None:
    fig, axes = plt.subplots(4, 2, figsize=(14.4, 18.6), dpi=180,
                             sharey=True, layout="constrained")
    left_names = ("A1", "B2", "A2", "B1")
    right_names = ("E1", "E2", "E3")
    sectors = {sector.name: sector for sector in SECTORS}
    for row, name in enumerate(left_names):
        draw_sector(axes[row, 0], sectors[name], show_ylabel=True)
    for row, name in enumerate(right_names):
        draw_sector(axes[row, 1], sectors[name], show_ylabel=False)
    axes[3, 1].axis("off")
    fig.suptitle("Square-antiprismatic bonding: ligand s/p SALCs + central s/p/d",
                 fontsize=17, color="#1e2d38", fontweight="medium")
    fig.supxlabel("Every MO bar represents ×1 or ×2 symmetry partners; connectors show ≥8% component weight.",
                  fontsize=9, color="#56636e")
    fig.savefig(FIGURE_PNG, facecolor="white", bbox_inches="tight")
    fig.savefig(FIGURE_SVG, facecolor="white", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check-only", action="store_true",
                        help="validate the model without plotting")
    args = parser.parse_args()
    validate()
    if args.check_only:
        print("Validated 32 ligand, 9 central, and 41 MO spatial orbitals.")
    else:
        draw()
        print(f"Generated {FIGURE_PNG} and {FIGURE_SVG}")


if __name__ == "__main__":
    main()
