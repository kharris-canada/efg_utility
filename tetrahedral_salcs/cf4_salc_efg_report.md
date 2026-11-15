# CF₄: SALC-resolved electric-field-gradient contributions

## Total from the occupied SALC projections

The summed traceless EFG tensor at the central carbon is

```text
[+0.000e+00  +0.000e+00  +0.000e+00]
[+0.000e+00  +0.000e+00  +0.000e+00]
[+0.000e+00  +0.000e+00  +0.000e+00]
```

`VZZ = +0.000000e+00 V/m²`; `eta = undefined`.

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
| 1a₁ (A1) | -5.597 | 1 | occupied |
| 1t₂ (T2) | -5.196 | 3 | occupied |
| 2a₁ (A1) | -4.243 | 1 | occupied |
| 2t₂ (T2) | -3.830 | 3 | occupied |
| 3t₂ (T2) | -2.144 | 3 | occupied |
| 1e (E) | -2.000 | 2 | occupied |
| 1t₁ (T1) | -2.000 | 3 | occupied |
| 3a₁ (A1) | -0.860 | 1 | virtual |
| 4t₂ (T2) | +0.970 | 3 | virtual |

For A₁ and T₂, each SALC's projected occupation is twice the sum of its squared
coefficients in the occupied MOs of that symmetry. These are fractional
projections, not integer electrons assigned to isolated SALCs: the toy MOs
mix ligand SALCs with central-carbon s/p functions, so an occupied two-electron
MO can place only a fraction of its population in a particular ligand SALC.
Its row is placed at the
lowest-energy occupied MO where it has its largest SALC weight. E and T₁ are
nonbonding in this no-d model and therefore retain two electrons per partner.
This bookkeeping assigns 30.858 electrons to
the ligand SALC projections; the remainder of the 32-electron population is
on central-atom components and is intentionally omitted from the EFG sum.

The energy matrix is the explicitly illustrative model from the tetrahedral
bonding diagram with its d orbitals removed: ligand s, radial pσ, and
tangential pπ start at −5.0, −3.7, and −2.0 arbitrary units; carbon 2s and 2p
start at −2.0 and +0.5. The A₁ couplings are 1.2 and 1.5; the T₂ couplings
from ligand s/radial-p/tangential-p to carbon 2p are 1.0/0.8/0.7. These are
teaching values, not a fitted or spectroscopic CF₄ calculation. The geometry
uses `d = 1.32 Å` (representative of the reported
C–F distance near 1.32 Å) and `Zeff = 6.0` for all ligand hydrogenic
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
Displayed entries smaller than 1e-12 times the
largest individual SALC tensor scale are set to zero; this suppresses numerical
roundoff, not physical contributions.

### 1a₁ (A1), energy -5.597 toy units; 1-fold

| a1-A-1 |
| :-- |
| Projected occupation: 1.8767 e |
| [(+0.000e+00, +0.000e+00, +0.000e+00)<br>(+0.000e+00, +0.000e+00, +0.000e+00)<br>(+0.000e+00, +0.000e+00, +0.000e+00)] |

| Quantity | EFG tensor (V/m²) |
| :-- | :-- |
| Sum over 1 partner(s) | [(+0.000e+00, +0.000e+00, +0.000e+00)<br>(+0.000e+00, +0.000e+00, +0.000e+00)<br>(+0.000e+00, +0.000e+00, +0.000e+00)] |

### 1t₂ (T2), energy -5.196 toy units; 3-fold

| t2-A-1 | t2-A-2 |
| :-- | :-- |
| Projected occupation: 1.9496 e | Projected occupation: 1.9496 e |
| [(-1.868e+19, -3.235e+19, -4.576e+19)<br>(-3.235e+19, +1.868e+19, -2.642e+19)<br>(-4.576e+19, -2.642e+19, +0.000e+00)] | [(-1.868e+19, +3.235e+19, +4.576e+19)<br>(+3.235e+19, +1.868e+19, -2.642e+19)<br>(+4.576e+19, -2.642e+19, +0.000e+00)] |

| t2-A-3 |
| :-- |
| Projected occupation: 1.9496 e |
| [(+3.736e+19, +0.000e+00, +0.000e+00)<br>(+0.000e+00, -3.736e+19, +5.283e+19)<br>(+0.000e+00, +5.283e+19, +0.000e+00)] |

| Quantity | EFG tensor (V/m²) |
| :-- | :-- |
| Sum over 3 partner(s) | [(+0.000e+00, +0.000e+00, +0.000e+00)<br>(+0.000e+00, +0.000e+00, +0.000e+00)<br>(+0.000e+00, +0.000e+00, +0.000e+00)] |

### 2a₁ (A1), energy -4.243 toy units; 1-fold

| a1-B-1 |
| :-- |
| Projected occupation: 1.5907 e |
| [(+0.000e+00, +0.000e+00, +0.000e+00)<br>(+0.000e+00, +0.000e+00, +0.000e+00)<br>(+0.000e+00, +0.000e+00, +0.000e+00)] |

| Quantity | EFG tensor (V/m²) |
| :-- | :-- |
| Sum over 1 partner(s) | [(+0.000e+00, +0.000e+00, +0.000e+00)<br>(+0.000e+00, +0.000e+00, +0.000e+00)<br>(+0.000e+00, +0.000e+00, +0.000e+00)] |

### 2t₂ (T2), energy -3.830 toy units; 3-fold

| t2-B-1 | t2-B-2 |
| :-- | :-- |
| Projected occupation: 1.9473 e | Projected occupation: 1.9473 e |
| [(-1.943e+19, -3.365e+19, -4.759e+19)<br>(-3.365e+19, +1.943e+19, -2.748e+19)<br>(-4.759e+19, -2.748e+19, +0.000e+00)] | [(-1.943e+19, +3.365e+19, +4.759e+19)<br>(+3.365e+19, +1.943e+19, -2.748e+19)<br>(+4.759e+19, -2.748e+19, +0.000e+00)] |

| t2-B-3 |
| :-- |
| Projected occupation: 1.9473 e |
| [(+3.886e+19, +0.000e+00, +0.000e+00)<br>(+0.000e+00, -3.886e+19, +5.495e+19)<br>(+0.000e+00, +5.495e+19, +0.000e+00)] |

| Quantity | EFG tensor (V/m²) |
| :-- | :-- |
| Sum over 3 partner(s) | [(+0.000e+00, +0.000e+00, +0.000e+00)<br>(+0.000e+00, +0.000e+00, +0.000e+00)<br>(+0.000e+00, +0.000e+00, +0.000e+00)] |

### 3t₂ (T2), energy -2.144 toy units; 3-fold

| t2-C-1 | t2-C-2 |
| :-- | :-- |
| Projected occupation: 1.9001 e | Projected occupation: 1.9001 e |
| [(-2.038e+19, -3.531e+19, -4.993e+19)<br>(-3.531e+19, +2.038e+19, -2.883e+19)<br>(-4.993e+19, -2.883e+19, +0.000e+00)] | [(-2.038e+19, +3.531e+19, +4.993e+19)<br>(+3.531e+19, +2.038e+19, -2.883e+19)<br>(+4.993e+19, -2.883e+19, +0.000e+00)] |

| t2-C-3 |
| :-- |
| Projected occupation: 1.9001 e |
| [(+4.077e+19, +0.000e+00, +0.000e+00)<br>(+0.000e+00, -4.077e+19, +5.765e+19)<br>(+0.000e+00, +5.765e+19, +0.000e+00)] |

| Quantity | EFG tensor (V/m²) |
| :-- | :-- |
| Sum over 3 partner(s) | [(+0.000e+00, +0.000e+00, +0.000e+00)<br>(+0.000e+00, +0.000e+00, +0.000e+00)<br>(+0.000e+00, +0.000e+00, +0.000e+00)] |

### 1e (E), energy -2.000 toy units; 2-fold

| e-C-1 | e-C-2 |
| :-- | :-- |
| Projected occupation: 2.0000 e | Projected occupation: 2.0000 e |
| [(+3.791e+19, +0.000e+00, +0.000e+00)<br>(+0.000e+00, -3.791e+19, +5.361e+19)<br>(+0.000e+00, +5.361e+19, +0.000e+00)] | [(-3.791e+19, +0.000e+00, +0.000e+00)<br>(+0.000e+00, +3.791e+19, -5.361e+19)<br>(+0.000e+00, -5.361e+19, +0.000e+00)] |

| Quantity | EFG tensor (V/m²) |
| :-- | :-- |
| Sum over 2 partner(s) | [(+0.000e+00, +0.000e+00, +0.000e+00)<br>(+0.000e+00, +0.000e+00, +0.000e+00)<br>(+0.000e+00, +0.000e+00, +0.000e+00)] |

### 1t₁ (T1), energy -2.000 toy units; 3-fold

| t₁-C-1 | t₁-C-2 |
| :-- | :-- |
| Projected occupation: 2.0000 e | Projected occupation: 2.0000 e |
| [(+1.512e+19, +2.619e+19, +3.705e+19)<br>(+2.619e+19, -1.512e+19, +2.139e+19)<br>(+3.705e+19, +2.139e+19, +0.000e+00)] | [(+1.512e+19, -2.619e+19, -3.705e+19)<br>(-2.619e+19, -1.512e+19, +2.139e+19)<br>(-3.705e+19, +2.139e+19, +0.000e+00)] |

| t₁-C-3 |
| :-- |
| Projected occupation: 2.0000 e |
| [(-3.025e+19, +0.000e+00, +0.000e+00)<br>(+0.000e+00, +3.025e+19, -4.278e+19)<br>(+0.000e+00, -4.278e+19, +0.000e+00)] |

| Quantity | EFG tensor (V/m²) |
| :-- | :-- |
| Sum over 3 partner(s) | [(+0.000e+00, +0.000e+00, +0.000e+00)<br>(+0.000e+00, +0.000e+00, +0.000e+00)<br>(+0.000e+00, +0.000e+00, +0.000e+00)] |

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
7.564e+05 V/m², roughly 7.3e-16
of the sum-of-contributions scale. This is numerical residue, not a tensor
that needs diagonalization into a PAS.
