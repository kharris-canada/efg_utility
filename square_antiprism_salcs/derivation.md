# Ligand s- and p-orbital SALCs of a square antiprism

This is a complete, explicit basis for **eight equivalent ligand sites**, each carrying one s orbital and three p orbitals: 32 ligand atomic orbitals (AOs) and 32 symmetry-adapted linear combinations (SALCs). The central atom is at the origin but its own orbitals are **not** included in the 32. The geometry has exact D₄d symmetry for any positive radius and half-height. The companion [figure gallery](gallery.md) draws every SALC.

The coefficient signs describe **orbital phase**, not electron charge or occupation. We treat AOs on different sites as orthonormal when normalizing the tabulated combinations. Real overlapping AOs require an overlap-aware normalization, but the symmetry labels and relative phase patterns remain the same.

## 1. Label the eight sites

Let ρ > 0 be the radius of either square and h > 0 its distance above or below the centre. Number the sites in the order visited by an eightfold improper rotation:

$$
\theta_j=\frac{j\pi}{4},\qquad
z_j=(-1)^j h,\qquad
\mathbf R_j=(\rho\cos\theta_j,\,\rho\sin\theta_j,\,z_j),
\quad j=0,\ldots,7.
$$

Even j are on the **top** square and odd j on the **bottom** square. The bottom is rotated 45° relative to the top. We name the sites T0 = j0, B0 = j1, T1 = j2, B1 = j3, T2 = j4, B2 = j5, T3 = j6, and B3 = j7. Here `B0` as a *site name* must not be confused with the B₁/B₂ *irrep names*.

The dimensions ρ and h affect the shapes and overlaps of real orbitals, but not these symmetry coefficients. For the drawings we choose ρ = 1 and h/ρ = 2<sup>−3/4</sup> ≈ 0.5946, which makes the top-square and side edges equal; the derivation does not require that particular ratio.

## 2. Generate D₄d with two matrices

Take S to be a 45° rotation about z followed by reflection across the xy plane, and M to be reflection across the xz plane:

$$
S=\begin{pmatrix}
1/\sqrt2&-1/\sqrt2&0\\
1/\sqrt2& 1/\sqrt2&0\\
0&0&-1
\end{pmatrix},
\qquad
M=\begin{pmatrix}1&0&0\\0&-1&0\\0&0&1\end{pmatrix}.
$$

These operations send site j to j + 1 and −j, respectively, with indices taken modulo 8. They satisfy S⁸ = M² = I and MSM = S⁻¹. Hence the eight powers of S and the eight operations MSᵏ make all 16 members of D₄d. In the usual class notation, S is S₈, S² is C₄, S⁴ is C₂, the four MS<sup>even</sup> are mirror planes σ<sub>d</sub>, and the four MS<sup>odd</sup> are perpendicular twofold rotations C₂′. Our xz mirror plane is a rotated Cartesian setting of the standard group; the irrep labels and characters are unchanged.

## 3. Choose a local p frame at each ligand

At site j, let **e**<sub>r,j</sub> be the horizontal outward radial direction, **e**<sub>φ,j</sub> the horizontal tangential direction, and **e**<sub>z</sub> the vertical unit vector. Define L = √(ρ² + h²) and three mutually perpendicular p axes:

$$
\begin{aligned}
\mathbf e_{r,j}&=(\cos\theta_j,\sin\theta_j,0),
&\mathbf e_{\phi,j}&=(-\sin\theta_j,\cos\theta_j,0),\\
\mathbf u_j&=\frac{\rho\mathbf e_{r,j}+z_j\mathbf e_z}{L}
&&\text{(bond-directed, outward)},\\
\mathbf v_j&=\frac{z_j\mathbf e_{r,j}-\rho\mathbf e_z}{L}
&&\text{(meridional, perpendicular to the bond)},\\
\mathbf w_j&=\mathbf e_{\phi,j}
&&\text{(tangential, perpendicular to the bond)}.
\end{aligned}
$$

The four AO families are s<sub>j</sub>, σ<sub>j</sub> (p along **u**), μ<sub>j</sub> (p along **v**), and τ<sub>j</sub> (p along **w**). Each p function is **positive on its specified +axis** and negative on the opposite lobe. This phase convention is needed to state unambiguous coefficients. Choosing all bond-directed p axes inward instead would multiply every σ SALC by −1 without changing its symmetry.

Applying the two matrices to both the **centre** and the **p direction** gives the crucial sign rules:

$$
\begin{array}{c|cc|cc}
\text{AO family }F & S F_j & M F_j & a_F & b_F\\\hline
s      & s_{j+1}      & s_{-j}      &+1&+1\\
\sigma & \sigma_{j+1} & \sigma_{-j} &+1&+1\\
\mu    &-\mu_{j+1}    & \mu_{-j}    &-1&+1\\
\tau   & \tau_{j+1}   &-\tau_{-j}   &+1&-1
\end{array}
$$

The minus signs are *not* optional: under S, the meridional p axis reverses relative to the chosen local axis at the new site; under M, the tangential p axis reverses. In shorthand, S F<sub>j</sub> = a<sub>F</sub>F<sub>j+1</sub> and M F<sub>j</sub> = b<sub>F</sub>F<sub>−j</sub>. Notice that a symmetry operation may move a p orbital **and** change its phase.

For instance, S changes z<sub>j</sub> to −z<sub>j</sub> and rotates the horizontal radial vector to the next site. Therefore

$$
S\mathbf v_j
=\frac{z_j\mathbf e_{r,j+1}+\rho\mathbf e_z}{L}
=-\mathbf v_{j+1}.
$$

This is why the meridional row carries a minus sign even though its new centre is simply site j + 1.

## 4. Count the required irreps

Here is the D₄d character table in the class order used below ([Quanty point-group reference](https://www.quanty.org/physics_chemistry/point_groups/d4d)):

| Irrep | E (1) | S₈ (2) | C₄ (2) | S₈³ (2) | C₂ (1) | C₂′ (4) | σ<sub>d</sub> (4) |
|:--|--:|--:|--:|--:|--:|--:|--:|
| A₁ | 1 | 1 | 1 | 1 | 1 | 1 | 1 |
| A₂ | 1 | 1 | 1 | 1 | 1 | −1 | −1 |
| B₁ | 1 | −1 | 1 | −1 | 1 | 1 | −1 |
| B₂ | 1 | −1 | 1 | −1 | 1 | −1 | 1 |
| E₁ | 2 | √2 | 0 | −√2 | −2 | 0 | 0 |
| E₂ | 2 | 0 | −2 | 0 | 2 | 0 | 0 |
| E₃ | 2 | −√2 | 0 | √2 | −2 | 0 | 0 |

For any one of our eight-AO families, the identity fixes all eight sites, so its representation character is 8. Every *other* non-mirror operation moves all sites, giving character 0. A mirror fixes exactly two sites. On those two AOs, s, σ, and μ are unchanged, while τ reverses sign. Thus, in the same class order,

$$
\chi_{s}=\chi_{\sigma}=\chi_{\mu}=(8,0,0,0,0,0,2),
\qquad
\chi_{\tau}=(8,0,0,0,0,0,-2).
$$

The number of copies of irrep Γ follows from the character inner product (the numbers in parentheses atop the character table are class sizes):

$$
n_{\Gamma}=\frac1{16}\sum_{\text{classes }C}
|C|\,\chi_F(C)\,\chi_{\Gamma}(C)^{*}.
$$

For example, the s-family has one A₁ SALC because (8 × 1 + 4 × 2 × 1)/16 = 1, but no A₂ SALC because (8 × 1 + 4 × 2 × −1)/16 = 0. Repeating this for the seven irreps yields

$$
\begin{aligned}
\Gamma_s=\Gamma_{\sigma}=\Gamma_{\mu}
 &= A_1\oplus B_2\oplus E_1\oplus E_2\oplus E_3,\\
\Gamma_{\tau}
 &= A_2\oplus B_1\oplus E_1\oplus E_2\oplus E_3,\\
\Gamma_{\rm all}
 &= 3A_1\oplus A_2\oplus B_1\oplus3B_2
    \oplus4E_1\oplus4E_2\oplus4E_3.
\end{aligned}
$$

The final dimensions add to 3 + 1 + 1 + 3 + 8 + 8 + 8 = **32**, so no AO is lost. Counting characters tells us *which* irreps occur, but not their explicit coefficients. We construct those next.

## 5. Construct all 32 SALCs by a discrete Fourier basis

Because S walks through the sites in order, the natural basis consists of waves around the eight-site orbit. Let F stand for one AO family (s, σ, μ, or τ), and let F<sub>j</sub> mean its **local AO at site j**. To avoid confusing a site number with a Fourier mode, call the resulting SALCs Ψ<sup>(F)</sup> with a mode subscript:

$$
\begin{aligned}
\Psi^{(F)}_0&=\frac1{\sqrt8}\sum_{j=0}^{7}F_j,\\
\Psi^{(F)}_4&=\frac1{\sqrt8}\sum_{j=0}^{7}(-1)^jF_j,\\
\Psi^{(F)}_{m,c}&=\frac12\sum_{j=0}^{7}
\cos\!\left(\frac{mj\pi}{4}\right)F_j,\\
\Psi^{(F)}_{m,s}&=\frac12\sum_{j=0}^{7}
\sin\!\left(\frac{mj\pi}{4}\right)F_j,
\qquad m=1,2,3.
\end{aligned}
$$

The factors normalize the functions in the orthonormal-AO approximation: the constant and alternating coefficient vectors each have squared length 8, while every sine and cosine vector shown has squared length 4. Different Fourier vectors are orthogonal. Therefore each family contributes **two single functions and three two-partner sets**, for 1 + 1 + 2 + 2 + 2 = 8 SALCs.

The sign rules above identify their irreps without guessing from pictures. The single functions obey

$$
S\Psi^{(F)}_0=a_F\Psi^{(F)}_0,\quad
S\Psi^{(F)}_4=-a_F\Psi^{(F)}_4,\quad
M\Psi^{(F)}_0=b_F\Psi^{(F)}_0,\quad
M\Psi^{(F)}_4=b_F\Psi^{(F)}_4.
$$

For a pair with α = mπ/4,

$$
\begin{pmatrix}S\Psi^{(F)}_{m,c}\\S\Psi^{(F)}_{m,s}\end{pmatrix}
=a_F\begin{pmatrix}\cos\alpha&\sin\alpha\\-\sin\alpha&\cos\alpha\end{pmatrix}
\begin{pmatrix}\Psi^{(F)}_{m,c}\\\Psi^{(F)}_{m,s}\end{pmatrix},
\qquad
\begin{pmatrix}M\Psi^{(F)}_{m,c}\\M\Psi^{(F)}_{m,s}\end{pmatrix}
=\begin{pmatrix}b_F&0\\0&-b_F\end{pmatrix}
\begin{pmatrix}\Psi^{(F)}_{m,c}\\\Psi^{(F)}_{m,s}\end{pmatrix}.
$$

In particular, the trace of S in that pair is 2a<sub>F</sub>cos(mπ/4): +√2 means E₁, 0 means E₂, and −√2 means E₃. The complete assignment is:

| AO family | a<sub>F</sub> | b<sub>F</sub> | m = 0 | m = 4 | m = 1 pair | m = 2 pair | m = 3 pair |
|:--|--:|--:|:--|:--|:--|:--|:--|
| ligand s | +1 | +1 | A₁ | B₂ | E₁ | E₂ | E₃ |
| bond pσ | +1 | +1 | A₁ | B₂ | E₁ | E₂ | E₃ |
| meridional pπ | −1 | +1 | B₂ | A₁ | E₃ | E₂ | E₁ |
| tangential pπ | +1 | −1 | A₂ | B₁ | E₁ | E₂ | E₃ |

The meridional family's reversal of E₁ and E₃ is a useful test: rotating **positions alone** would miss the negative p-axis sign under S and label these states incorrectly. Under D₄d, a central s AO has A₁ symmetry, central p<sub>z</sub> has B₂ symmetry, and central (p<sub>x</sub>, p<sub>y</sub>) form E₁. Thus only ligand SALCs of A₁, B₂, or E₁ symmetry can mix with central s/p AOs; symmetry alone does not say how strongly they mix or which MOs are occupied.

This is **one symmetry-adapted basis**, not a unique set of molecular orbitals. For example, the three A₁ SALCs from different AO families can be mixed with one another and still have A₁ symmetry. The same is true of repeated E₁, E₂, E₃, and B₂ copies, provided entire partner spaces are transformed consistently. An energy calculation, not the character table alone, chooses the molecular-orbital combinations.

## 6. How to read the coefficient table

Every row below is **one complete normalized SALC**, not just a schematic sign pattern. The AO-family column tells you which local orbital at each vertex is used. Multiply each of the eight entries by the row's prefactor, then sum those eight AOs. For example, row 03 is

$$
\Psi^{(s)}_{1,c}=\tfrac12(s_0+t s_1-t s_3-s_4-t s_5+t s_7),
\qquad t=1/\sqrt2.
$$

Zero entries mean no contribution at that site. The c/s members of each E irrep are **partners**: they transform into linear combinations of one another, not individually into themselves. Their particular phase and basis orientation are conventions; any orthogonal rotation within a partner pair gives an equally valid SALC pair. The coefficient table fixes one reproducible choice.

The plotting code additionally checks the site transformations, all 16 group matrices, orthogonality, and the characters of every one- and two-dimensional block. The group character labels follow the [Quanty D₄d reference](https://www.quanty.org/physics_chemistry/point_groups/d4d).

## Complete coefficient table

Let **t = 1/√2**. The eight site columns are in S₈-orbit order. An entry `t` in a row with prefactor 1/2 means an actual coefficient of 1/(2√2); `−1` means minus the prefactor. In the SALC column, `Ψ(F;m)` or `Ψ(F;m,c/s)` identifies AO family F and **Fourier mode m**, not a site. The figure number links to the corresponding gallery image.

| Fig. | AO family | SALC | Irrep | Prefactor | T0 (j0) | B0 (j1) | T1 (j2) | B1 (j3) | T2 (j4) | B2 (j5) | T3 (j6) | B3 (j7) |
|--:|:--|:--|:--|:--|--:|--:|--:|--:|--:|--:|--:|--:|
| 01 | s | Ψ(s;0) | A₁ | 1/√8 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 |
| 02 | s | Ψ(s;4) | B₂ | 1/√8 | 1 | −1 | 1 | −1 | 1 | −1 | 1 | −1 |
| 03 | s | Ψ(s;1,c) | E₁ | 1/2 | 1 | t | 0 | −t | −1 | −t | 0 | t |
| 04 | s | Ψ(s;1,s) | E₁ | 1/2 | 0 | t | 1 | t | 0 | −t | −1 | −t |
| 05 | s | Ψ(s;2,c) | E₂ | 1/2 | 1 | 0 | −1 | 0 | 1 | 0 | −1 | 0 |
| 06 | s | Ψ(s;2,s) | E₂ | 1/2 | 0 | 1 | 0 | −1 | 0 | 1 | 0 | −1 |
| 07 | s | Ψ(s;3,c) | E₃ | 1/2 | 1 | −t | 0 | t | −1 | t | 0 | −t |
| 08 | s | Ψ(s;3,s) | E₃ | 1/2 | 0 | t | −1 | t | 0 | −t | 1 | −t |
| 09 | pσ | Ψ(σ;0) | A₁ | 1/√8 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 |
| 10 | pσ | Ψ(σ;4) | B₂ | 1/√8 | 1 | −1 | 1 | −1 | 1 | −1 | 1 | −1 |
| 11 | pσ | Ψ(σ;1,c) | E₁ | 1/2 | 1 | t | 0 | −t | −1 | −t | 0 | t |
| 12 | pσ | Ψ(σ;1,s) | E₁ | 1/2 | 0 | t | 1 | t | 0 | −t | −1 | −t |
| 13 | pσ | Ψ(σ;2,c) | E₂ | 1/2 | 1 | 0 | −1 | 0 | 1 | 0 | −1 | 0 |
| 14 | pσ | Ψ(σ;2,s) | E₂ | 1/2 | 0 | 1 | 0 | −1 | 0 | 1 | 0 | −1 |
| 15 | pσ | Ψ(σ;3,c) | E₃ | 1/2 | 1 | −t | 0 | t | −1 | t | 0 | −t |
| 16 | pσ | Ψ(σ;3,s) | E₃ | 1/2 | 0 | t | −1 | t | 0 | −t | 1 | −t |
| 17 | pπ meridional | Ψ(μ;0) | B₂ | 1/√8 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 |
| 18 | pπ meridional | Ψ(μ;4) | A₁ | 1/√8 | 1 | −1 | 1 | −1 | 1 | −1 | 1 | −1 |
| 19 | pπ meridional | Ψ(μ;1,c) | E₃ | 1/2 | 1 | t | 0 | −t | −1 | −t | 0 | t |
| 20 | pπ meridional | Ψ(μ;1,s) | E₃ | 1/2 | 0 | t | 1 | t | 0 | −t | −1 | −t |
| 21 | pπ meridional | Ψ(μ;2,c) | E₂ | 1/2 | 1 | 0 | −1 | 0 | 1 | 0 | −1 | 0 |
| 22 | pπ meridional | Ψ(μ;2,s) | E₂ | 1/2 | 0 | 1 | 0 | −1 | 0 | 1 | 0 | −1 |
| 23 | pπ meridional | Ψ(μ;3,c) | E₁ | 1/2 | 1 | −t | 0 | t | −1 | t | 0 | −t |
| 24 | pπ meridional | Ψ(μ;3,s) | E₁ | 1/2 | 0 | t | −1 | t | 0 | −t | 1 | −t |
| 25 | pπ tangential | Ψ(τ;0) | A₂ | 1/√8 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 |
| 26 | pπ tangential | Ψ(τ;4) | B₁ | 1/√8 | 1 | −1 | 1 | −1 | 1 | −1 | 1 | −1 |
| 27 | pπ tangential | Ψ(τ;1,c) | E₁ | 1/2 | 1 | t | 0 | −t | −1 | −t | 0 | t |
| 28 | pπ tangential | Ψ(τ;1,s) | E₁ | 1/2 | 0 | t | 1 | t | 0 | −t | −1 | −t |
| 29 | pπ tangential | Ψ(τ;2,c) | E₂ | 1/2 | 1 | 0 | −1 | 0 | 1 | 0 | −1 | 0 |
| 30 | pπ tangential | Ψ(τ;2,s) | E₂ | 1/2 | 0 | 1 | 0 | −1 | 0 | 1 | 0 | −1 |
| 31 | pπ tangential | Ψ(τ;3,c) | E₃ | 1/2 | 1 | −t | 0 | t | −1 | t | 0 | −t |
| 32 | pπ tangential | Ψ(τ;3,s) | E₃ | 1/2 | 0 | t | −1 | t | 0 | −t | 1 | −t |
