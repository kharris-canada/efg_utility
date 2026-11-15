# Ligand s- and p-orbital SALCs for a tetrahedron (T<sub>d</sub>)

This is one **complete basis of 16 ligand SALCs**: four equivalent ligand sites, each with one s orbital and three real p orbitals. The central atom's orbitals are not included. The companion [illustrated gallery](gallery.md) draws all 16 combinations, and the [bonding diagram](bonding_diagram.md) shows how they can mix with a central atom's s/p/d orbitals. A SALC is a combination of *wavefunctions*; a positive or negative coefficient denotes orbital phase, not charge or occupation.

We assume that the 16 starting atomic orbitals (AOs) are orthonormal when specifying normalized coefficients. Real AOs on different ligands overlap, so a quantitative molecular-orbital calculation must account for the overlap matrix. The symmetry labels and relative coefficient patterns below do not depend on the tetrahedron's size.

## 1. Put the four ligands in a definite Cartesian orientation

Write the position of ligand j as the vector of signs below multiplied by R/√3. Thus every ligand is a distance R from the central atom. We use the same letters A–D for the sites and the labels in the figures.

| Site | x<sub>j</sub> | y<sub>j</sub> | z<sub>j</sub> | Position |
|:--|--:|--:|--:|:--|
| A | +1 | +1 | +1 | (R/√3)(+1,+1,+1) |
| B | +1 | −1 | −1 | (R/√3)(+1,−1,−1) |
| C | −1 | +1 | −1 | (R/√3)(−1,+1,−1) |
| D | −1 | −1 | +1 | (R/√3)(−1,−1,+1) |

The sign vector at site j is **t**<sub>j</sub> = (x<sub>j</sub>, y<sub>j</sub>, z<sub>j</sub>). Every sign vector obeys x<sub>j</sub>y<sub>j</sub>z<sub>j</sub> = +1. This gives the useful identities x<sub>j</sub>y<sub>j</sub> = z<sub>j</sub>, y<sub>j</sub>z<sub>j</sub> = x<sub>j</sub>, and z<sub>j</sub>x<sub>j</sub> = y<sub>j</sub>.

At each site, use s<sub>j</sub> and three p functions p<sub>xj</sub>, p<sub>yj</sub>, p<sub>zj</sub>. Assume the three p functions on one ligand share the same radial shape, so they transform as Cartesian components of a vector. The positive lobe of p<sub>xj</sub> points along the **global +x axis**, and similarly for y and z. These Cartesian p functions are a convenient *basis*, not three symmetry-preserved individual orbitals: a tetrahedral operation generally moves a site **and** mixes its p directions.

## 2. See the site and p-orbital transformations explicitly

Here are four representative T<sub>d</sub> operations in this orientation. The first cycles x → z → y → x and is a threefold rotation about the A–origin axis; the second is a twofold rotation about x; the third is a fourfold improper rotation about z; the last is a mirror exchanging x and y.

$$
C_3=\begin{pmatrix}0&1&0\\0&0&1\\1&0&0\end{pmatrix},\qquad
C_2=\begin{pmatrix}1&0&0\\0&-1&0\\0&0&-1\end{pmatrix},\qquad
S_4=\begin{pmatrix}0&-1&0\\1&0&0\\0&0&-1\end{pmatrix},\qquad
M=\begin{pmatrix}0&1&0\\1&0&0\\0&0&1\end{pmatrix}.
$$

For example, M fixes A and D but exchanges B and C. It also exchanges the p<sub>x</sub> and p<sub>y</sub> directions. Consequently **the AO p<sub>xB</sub> becomes p<sub>yC</sub>**, not p<sub>xC</sub>. Similarly, the C₃ matrix fixes site A, cycles B → D → C → B, and sends a p<sub>x</sub> direction to p<sub>z</sub>. This is why simply copying the site coefficients of an s SALC cannot generate all p SALCs.

More generally, all 24 T<sub>d</sub> matrices in this orientation are **signed permutation matrices with an even number of sign reversals**. Permuting coordinates leaves xyz unchanged, and an even number of sign reversals also leaves xyz unchanged; these matrices therefore permute the four vertices above. The 24 distinct matrices exhaust the tetrahedron's 24 vertex permutations. The identity, together with the four matrices above, represents the five classes E, 8C₃, 3C₂, 6S₄, and 6σ<sub>d</sub>. The [Quanty T<sub>d</sub> reference](https://www.quanty.org/physics_chemistry/point_groups/td) gives the same class order and character table:

| Irrep | E (1) | C₃ (8) | C₂ (3) | S₄ (6) | σ<sub>d</sub> (6) |
|:--|--:|--:|--:|--:|--:|
| A₁ | 1 | 1 | 1 | 1 | 1 |
| A₂ | 1 | 1 | 1 | −1 | −1 |
| E | 2 | −1 | 2 | 0 | 0 |
| T₁ | 3 | 0 | −1 | 1 | −1 |
| T₂ | 3 | 0 | −1 | −1 | 1 |

The parenthesized numbers in the table heading are the numbers of operations in each class. In this orientation, an ordinary Cartesian vector (x,y,z) transforms as T₂, whereas an axial vector (such as a rotation axis) transforms as T₁. These differ under improper operations.

## 3. Count the SALCs before constructing them

For the four ligand s orbitals, the character of an operation is the number of sites it fixes. The identity fixes four; a C₃ fixes one; a C₂ or S₄ fixes none; a mirror fixes two. Thus, in the class order above,

$$
\chi_s=(4,1,0,0,2).
$$

For all 12 p AOs, a fixed site contributes the trace of that operation's **3 × 3 Cartesian p-direction matrix**. The identity contributes 4 × 3 = 12; C₃ fixes one site but its matrix has trace zero; the mirror fixes two sites and has matrix trace one. The other classes fix no sites. Therefore

$$
\chi_p=(12,0,0,0,2).
$$

Use the character inner product to count each irrep Γ:

$$
n_\Gamma=\frac1{24}\sum_{\text{classes }K}|K|\,\chi(K)\,\chi_\Gamma(K)^*.
$$

For example, the ligand-s representation contains one T₂ because [4(3) + 8(1)(0) + 6(2)(1)]/24 = 1. Carrying out the five inner products gives

$$
\begin{aligned}
\Gamma_s&=A_1\oplus T_2,\\
\Gamma_p&=A_1\oplus E\oplus T_1\oplus 2T_2,\\
\Gamma_{s+p}&=2A_1\oplus E\oplus T_1\oplus 3T_2.
\end{aligned}
$$

The dimensions are 2(1) + 2 + 3 + 3(3) = **16**, as required. A₂ does not occur in this s/p ligand basis. The character calculation predicts the symmetry content, but not the actual coefficients.

For comparison, a central-atom s orbital has A₁ symmetry and its three p orbitals form T₂. Therefore ligand A₁ and T₂ SALCs can mix with central s/p orbitals, while ligand E and T₁ SALCs have no central s/p partner of matching symmetry. Symmetry permits mixing; it does not require it.

## 4. Construct the four s SALCs

The four site-sign columns 1, x<sub>j</sub>, y<sub>j</sub>, z<sub>j</sub> are mutually orthogonal. Each has squared length four, so division by two normalizes it:

$$
\Psi(s;0)=\frac12\sum_j s_j,\qquad
\Psi(s;\alpha)=\frac12\sum_j t_{j\alpha}s_j,
\quad \alpha=x,y,z.
$$

The uniform combination is A₁. The x, y, z partners transform exactly like the Cartesian components of a polar vector, so they form one T₂ triplet. For example, Ψ(s;x) has the signs (+,+,−,−) at A,B,C,D. Unlike a one-dimensional SALC, one T₂ partner need not map to itself: a symmetry operation can turn Ψ(s;x) into Ψ(s;y) or another linear combination within the triplet.

## 5. Organize the 12 p AOs before writing their coefficients

Let **p**<sub>j</sub> = (p<sub>xj</sub>, p<sub>yj</sub>, p<sub>zj</sub>) denote the three p functions at site j. Define a uniform Cartesian triplet and, for any real 3 × 3 matrix K, a matrix-weighted combination:

$$
P_\alpha=\frac12\sum_j p_{\alpha j},\qquad
\Phi_K=\frac12\sum_j (K\mathbf t_j)\mathbin{\cdot}\mathbf p_j.
$$

The four sign vectors satisfy

$$
\sum_j\mathbf t_j=\mathbf0,\qquad
\sum_j\mathbf t_j\mathbf t_j^{\mathsf T}=4I.
$$

Consequently, the three P functions are orthonormal and perpendicular to all nine independent Φ functions. More precisely, assuming orthonormal AOs,

$$
\langle P_\alpha|P_\beta\rangle=\delta_{\alpha\beta},\quad
\langle P_\alpha|\Phi_K\rangle=0,\quad
\langle\Phi_K|\Phi_L\rangle=\operatorname{Tr}(K^{\mathsf T}L).
$$

That last identity makes the problem a familiar one: find an orthonormal basis of 3 × 3 matrices. Under a point-group matrix g, both the site signs and the p directions transform. Their simultaneous transformation gives

$$
gP_\alpha=\sum_\beta g_{\beta\alpha}P_\beta,
\qquad g\Phi_K=\Phi_{gKg^{\mathsf T}}.
$$

The nine matrix functions separate into an isotropic trace (A₁), two diagonal traceless matrices (E), three antisymmetric matrices (T₁), and three symmetric off-diagonal matrices (T₂). The uniform P triplet is another T₂. This reproduces the p-space decomposition calculated from characters. The T₁ assignment for antisymmetric matrices follows explicitly by writing such a matrix as K(**a**) **v** = **a** × **v**: under an improper operation, its axis acquires the extra determinant sign,

$$
gK(\mathbf a)g^{\mathsf T}=K\!\left((\det g)g\mathbf a\right).
$$

Thus the antisymmetric triplet is axial, whereas P and the symmetric off-diagonal triplet are polar.

To make the coefficients explicit, define the symmetric off-diagonal matrices by their nonzero entries:

$$
Q_x=\frac{E_{yz}+E_{zy}}{\sqrt2},\qquad
Q_y=\frac{E_{zx}+E_{xz}}{\sqrt2},\qquad
Q_z=\frac{E_{xy}+E_{yx}}{\sqrt2}.
$$

Here E<sub>ab</sub> is the matrix with a single 1 in row a, column b. The Q triplet and P triplet both transform as T₂, so they may be mixed *within corresponding x, y, z partners*.

### The radial A₁ and T₂ combinations

At site j the outward radial p orbital is

$$
\sigma_j=\frac{\mathbf t_j\mathbin{\cdot}\mathbf p_j}{\sqrt3}.
$$

Its four-site representation is identical to that of the s orbitals: one radial A₁ plus one radial T₂. Their normalized combinations are

$$
\Psi(\sigma;0)=\frac12\sum_j\sigma_j=\Phi_{I/\sqrt3},
\qquad
\Psi(\sigma;\alpha)=\frac12\sum_j t_{j\alpha}\sigma_j.
$$

For x, the tetrahedral identity x<sub>j</sub>y<sub>j</sub> = z<sub>j</sub> and its companions give

$$
\Psi(\sigma;x)
=\frac1{2\sqrt3}\sum_j
\bigl(p_{xj}+z_jp_{yj}+y_jp_{zj}\bigr)
=\frac{P_x+\sqrt2\Phi_{Q_x}}{\sqrt3}.
$$

Cyclic permutation of x, y, z gives the other two partners. This equation reveals why the radial T₂ is *not* simply the uniform Cartesian p triplet: it mixes P and Q. The orthogonal T₂ triplet is

$$
\Psi(\pi T_2;\alpha)
=\frac{\sqrt2P_\alpha-\Phi_{Q_\alpha}}{\sqrt3}.
$$

For x, this becomes

$$
\Psi(\pi T_2;x)=\frac1{2\sqrt6}\sum_j
\bigl(2p_{xj}-z_jp_{yj}-y_jp_{zj}\bigr).
$$

It is tangential **at every site**: its local coefficient vector has zero dot product with **t**<sub>j</sub>, since 2x<sub>j</sub> − 2y<sub>j</sub>z<sub>j</sub> = 0. The y and z partners work the same way.

### The tangential E and T₁ combinations

Two orthonormal traceless diagonal matrices and three orthonormal antisymmetric matrices give the remaining five p SALCs:

$$
\begin{aligned}
\Psi(E;\theta)&=\frac1{2\sqrt6}\sum_j
(-x_jp_{xj}-y_jp_{yj}+2z_jp_{zj}),\\
\Psi(E;\epsilon)&=\frac1{2\sqrt2}\sum_j
(x_jp_{xj}-y_jp_{yj}),\\
\Psi(T_1;x)&=\frac1{2\sqrt2}\sum_j
(-z_jp_{yj}+y_jp_{zj}),\\
\Psi(T_1;y)&=\frac1{2\sqrt2}\sum_j
(z_jp_{xj}-x_jp_{zj}),\\
\Psi(T_1;z)&=\frac1{2\sqrt2}\sum_j
(-y_jp_{xj}+x_jp_{yj}).
\end{aligned}
$$

Each E or T₁ coefficient vector is perpendicular to its site's bond direction. The E pair is diagonal and traceless; the T₁ matrices are antisymmetric, so they behave as an **axial** rather than a polar triplet under reflections and improper rotations. All five are tangential p combinations. Along with the tangential T₂ triplet above, they account for the eight-dimensional p space perpendicular to the four radial pσ orbitals.

These are **symmetry-adapted basis functions, not uniquely determined molecular orbitals**. In particular, the three T₂ copies (from s, radial pσ, and tangential p) can mix with one another while retaining T₂ symmetry, provided their partner triplets are mixed consistently. The two A₁ copies can also mix. A Hamiltonian and overlap model, not the character table by itself, determine any actual molecular-orbital coefficients and occupations.

## 6. Read the complete coefficient table

Every row below is **one complete normalized SALC**. Multiply each site entry by that row's prefactor, then sum the indicated local AOs. In an s row, a site entry is the coefficient of s<sub>j</sub>. In a p row, a triple is the coefficient vector of (p<sub>xj</sub>, p<sub>yj</sub>, p<sub>zj</sub>) at that site, in exactly that order. For example, the A entry in row 09, with prefactor 1/(2√6), means −p<sub>xA</sub>/(2√6) − p<sub>yA</sub>/(2√6) + 2p<sub>zA</sub>/(2√6). A zero component is an absent AO contribution.

The x/y/z entries in a T₁ or T₂ set are **partners of one multidimensional irrep**, not three individually invariant functions. Our phase and partner-axis choices are conventions; changing the sign of a complete SALC or orthogonally recombining partners gives an equally valid basis. The [gallery](gallery.md) uses exactly the coefficients below. Its Python generator checks all 24 group matrices, every invariant subspace, the character of each irrep block, orthonormality, local tangency where claimed, and every entry of this table.

## Complete coefficient table

The four site columns follow A, B, C, D as defined above. For each p row, **each triple is (pₓ,pᵧ,p_z)** before multiplication by the prefactor. Minus signs are orbital phases, not electronic charges.

| Fig. | AO sector | SALC | Irrep | Prefactor | A (+,+,+) | B (+,−,−) | C (−,+,−) | D (−,−,+) |
|--:|:--|:--|:--|:--|:--|:--|:--|:--|
| 01 | s | Ψ(s;0) | A₁ | 1/2 | 1 | 1 | 1 | 1 |
| 02 | s | Ψ(s;x) | T₂ | 1/2 | 1 | 1 | −1 | −1 |
| 03 | s | Ψ(s;y) | T₂ | 1/2 | 1 | −1 | 1 | −1 |
| 04 | s | Ψ(s;z) | T₂ | 1/2 | 1 | −1 | −1 | 1 |
| 05 | radial pσ | Ψ(σ;0) | A₁ | 1/(2√3) | (1,1,1) | (1,−1,−1) | (−1,1,−1) | (−1,−1,1) |
| 06 | radial pσ | Ψ(σ;x) | T₂ | 1/(2√3) | (1,1,1) | (1,−1,−1) | (1,−1,1) | (1,1,−1) |
| 07 | radial pσ | Ψ(σ;y) | T₂ | 1/(2√3) | (1,1,1) | (−1,1,1) | (−1,1,−1) | (1,1,−1) |
| 08 | radial pσ | Ψ(σ;z) | T₂ | 1/(2√3) | (1,1,1) | (−1,1,1) | (1,−1,1) | (−1,−1,1) |
| 09 | tangential E | Ψ(E;θ) | E | 1/(2√6) | (−1,−1,2) | (−1,1,−2) | (1,−1,−2) | (1,1,2) |
| 10 | tangential E | Ψ(E;ε) | E | 1/(2√2) | (1,−1,0) | (1,1,0) | (−1,−1,0) | (−1,1,0) |
| 11 | tangential T₁ | Ψ(T₁;x) | T₁ | 1/(2√2) | (0,−1,1) | (0,1,−1) | (0,1,1) | (0,−1,−1) |
| 12 | tangential T₁ | Ψ(T₁;y) | T₁ | 1/(2√2) | (1,0,−1) | (−1,0,−1) | (−1,0,1) | (1,0,1) |
| 13 | tangential T₁ | Ψ(T₁;z) | T₁ | 1/(2√2) | (−1,1,0) | (1,1,0) | (−1,−1,0) | (1,−1,0) |
| 14 | tangential T₂ | Ψ(πT₂;x) | T₂ | 1/(2√6) | (2,−1,−1) | (2,1,1) | (2,1,−1) | (2,−1,1) |
| 15 | tangential T₂ | Ψ(πT₂;y) | T₂ | 1/(2√6) | (−1,2,−1) | (1,2,−1) | (1,2,1) | (−1,2,1) |
| 16 | tangential T₂ | Ψ(πT₂;z) | T₂ | 1/(2√6) | (−1,−1,2) | (1,−1,2) | (−1,1,2) | (1,1,2) |
