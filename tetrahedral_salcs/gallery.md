# Illustrated T<sub>d</sub> tetrahedral ligand SALCs

These are the **16 ligand s/p SALCs** derived in [the mathematical analysis](derivation.md). Four sites (A–D) each carry one s and three real p atomic orbitals. The central atom's orbitals are not part of this set; their possible mixing is shown in the separate [s/p/d bonding diagram](bonding_diagram.md). The 3D drawings are schematic; they show wavefunction **phase**, not electron density or charge. Blue is positive phase, orange is negative phase. A p orbital has one lobe of each phase.

At each ligand, a linear combination of pₓ, pᵧ, and p<sub>z</sub> with a shared radial shape is drawn as **one oriented p orbital**, not as three separate lobed functions. The lower panel shows s coefficients as bars or p coefficients as a signed Cartesian heatmap. The heatmap numbers are rounded; the [complete table](derivation.md#complete-coefficient-table) gives exact values. A zero heatmap entry means that Cartesian AO is absent. The displayed tetrahedron uses the site-sign vectors A = (+,+,+), B = (+,−,−), C = (−,+,−), D = (−,−,+).

To reproduce all figures from the project root, install the isolated figure dependencies with `python -m pip install -r tetrahedral_salcs/requirements.txt`, then run `python tetrahedral_salcs/generate_figures.py`. The generator also validates the symmetry and published coefficient table. Run it with `--check-only` to validate without redrawing.

## Four ligand-s SALCs: A₁ + T₂

### 01 · A₁ · Ψ(s;0)

All four s orbitals enter with coefficient +1/2. This is the totally symmetric ligand-s combination.

![Uniform A1 ligand-s SALC, all four phases positive](figures/01_s_0.png)

### 02 · T₂ · Ψ(s;x)

Sites A and B, whose x coordinates are positive, contribute positively; C and D contribute negatively. This is the x partner of the ligand-s T₂ triplet.

![T2 x partner from ligand s orbitals, A and B positive](figures/02_s_x.png)

### 03 · T₂ · Ψ(s;y)

Sites A and C have positive phase; B and D have negative phase. A Td operation can mix this partner with the x and z partners.

![T2 y partner from ligand s orbitals, A and C positive](figures/03_s_y.png)

### 04 · T₂ · Ψ(s;z)

Sites A and D have positive phase; B and C have negative phase. Together, figures 02–04 form one three-dimensional irrep.

![T2 z partner from ligand s orbitals, A and D positive](figures/04_s_z.png)

## Four radial pσ SALCs: A₁ + T₂

Each local positive pσ direction points **outward along that ligand's bond**. Its Cartesian AO coefficients therefore change from site to site even when the local radial coefficient is the same.

### 05 · A₁ · Ψ(σ;0)

All four outward radial pσ orbitals enter with the same positive local coefficient +1/2. The 3D blue lobe is outward at every vertex; the heatmap resolves each into global pₓ, pᵧ, p_z components.

![Uniform radial p sigma A1 SALC with outward positive lobes](figures/05_sigma_0.png)

### 06 · T₂ · Ψ(σ;x)

The local radial coefficient follows the x-coordinate signs: positive at A and B, negative at C and D. The p-axis itself rotates with each bond.

![Radial p sigma T2 x partner](figures/06_sigma_x.png)

### 07 · T₂ · Ψ(σ;y)

The radial sign pattern is positive at A and C, negative at B and D. This is the y partner of the radial T₂ triplet.

![Radial p sigma T2 y partner](figures/07_sigma_y.png)

### 08 · T₂ · Ψ(σ;z)

The radial sign pattern is positive at A and D, negative at B and C. This completes the four-dimensional radial pσ subspace.

![Radial p sigma T2 z partner](figures/08_sigma_z.png)

## Two tangential E SALCs

The local p direction in these figures is **perpendicular to the bond at every site**. The two E partners may be orthogonally recombined, but one partner alone need not be invariant under a Td rotation.

### 09 · E · Ψ(E;θ)

The local Cartesian coefficient vector is proportional to (−x<sub>j</sub>, −y<sub>j</sub>, 2z<sub>j</sub>). The z component has twice the magnitude of either other Cartesian component. This is the θ member of our chosen E basis; the singled-out z direction is a **basis convention**, not a physical lowering of tetrahedral symmetry.

![Tangential E theta SALC with coefficients minus x minus y plus twice z](figures/09_E_theta.png)

### 10 · E · Ψ(E;ε)

The local vector is proportional to (x<sub>j</sub>, −y<sub>j</sub>, 0). There is no p_z contribution. Together with figure 09, this spans the E irrep.

![Tangential E epsilon SALC with no p z coefficients](figures/10_E_epsilon.png)

## Three tangential T₁ SALCs

These local p directions arise from a cross product of a Cartesian partner axis with the bond vector. As a set, they transform like an **axial vector**: reflections and improper rotations distinguish T₁ from the polar-vector T₂ sets.

### 11 · T₁ · Ψ(T₁;x)

The coefficient vector at site j is proportional to (0, −z<sub>j</sub>, y<sub>j</sub>). No pₓ AO contributes anywhere, though the resulting local p direction changes with the site.

![Tangential T1 x partner with zero p x coefficients](figures/11_T1_x.png)

### 12 · T₁ · Ψ(T₁;y)

Here the local vector is proportional to (z<sub>j</sub>, 0, −x<sub>j</sub>). The pᵧ row of the coefficient heatmap is zero.

![Tangential T1 y partner with zero p y coefficients](figures/12_T1_y.png)

### 13 · T₁ · Ψ(T₁;z)

The local vector is proportional to (−y<sub>j</sub>, x<sub>j</sub>, 0). There is no p_z contribution. Figures 11–13 make one T₁ triplet.

![Tangential T1 z partner with zero p z coefficients](figures/13_T1_z.png)

## Three tangential T₂ SALCs

These are the T₂ combinations **orthogonal to the radial T₂ triplet**. Their local vectors are also perpendicular to each bond, but they transform like a polar rather than an axial vector.

### 14 · T₂ · Ψ(πT₂;x)

The local coefficient vector is proportional to (2, −z<sub>j</sub>, −y<sub>j</sub>). Its dot product with the bond-sign vector is zero at every ligand because y<sub>j</sub>z<sub>j</sub> = x<sub>j</sub>.

![Tangential T2 x partner orthogonal to each bond](figures/14_piT2_x.png)

### 15 · T₂ · Ψ(πT₂;y)

The cyclic partner has local vector (−z<sub>j</sub>, 2, −x<sub>j</sub>). Its y component is twice the magnitude of either cross-component.

![Tangential T2 y partner orthogonal to each bond](figures/15_piT2_y.png)

### 16 · T₂ · Ψ(πT₂;z)

The last partner has local vector (−y<sub>j</sub>, −x<sub>j</sub>, 2). These three figures complete the eight tangential p SALCs and the full 16-function ligand basis.

![Tangential T2 z partner orthogonal to each bond](figures/16_piT2_z.png)
