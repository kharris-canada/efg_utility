# Why a tetrahedral site has zero electric-field gradient

For every symmetry operation that fixes the nucleus at the origin, the potential Hessian must satisfy

$$
\boxed{H=RHR^{T}.}
$$

Here R is the 3 × 3 matrix representing the operation. H is the Hessian of the electrostatic potential Φ at the nucleus: its entry H<sub>ij</sub> is the second derivative with respect to coordinates i and j. We will use this one condition to prove that the **traceless electric-field-gradient (EFG) tensor** is zero at a site with exact tetrahedral (T<sub>d</sub>) symmetry. No point-charge approximation is needed.

Why should the condition hold? A genuine symmetry operation takes the *entire source-charge distribution*—ligand nuclei and electron density—onto an indistinguishable distribution. In symbols, ρ(Rr) = ρ(r). Because Coulomb's law respects rotations and reflections, that distribution generates the same potential after the operation: Φ(Rr) = Φ(r). The Hessian records the local curvature of the *potential* at the nucleus, not the curvature of the *electron density* there: the potential at one point depends on charges throughout space. Nonetheless, that derived curvature must inherit the distribution's symmetry.

The boxed invariance equation, **H = RHRᵀ**, follows from a **second-order Taylor expansion** of the potential about the origin, which every operation fixes:

$$
\Phi(\mathbf r)=\Phi(0)+\mathbf g^{T}\mathbf r
+\tfrac12\mathbf r^{T}H\mathbf r+\cdots.
$$

Here r = (x, y, z)ᵀ is a column vector describing a small displacement from the origin; its transpose rᵀ = (x, y, z) is a row vector. Thus rᵀHr is a single number; half of it is the quadratic curvature term in the expansion. The vector g = ∇Φ(0) is the gradient of the potential at the origin—the negative of the electric-field vector. The three displayed terms give the potential's value, linear change, and quadratic curvature, respectively. H is the Hessian defined above.

Replacing r by Rr changes the quadratic term to

$$
\tfrac12\mathbf r^{T}R^{T}HR\mathbf r.
$$

Because R is a linear transformation, it does not mix terms of different orders. Since the potential is unchanged for every nearby r, equality of the quadratic terms gives H = RᵀHR. **Every point-group operation is orthogonal**: proper rotations, mirror reflections, improper rotations (rotation followed by reflection), and inversion all satisfy R⁻¹ = Rᵀ. Multiplying H = RᵀHR on the left by R and on the right by Rᵀ therefore gives the boxed equation H = RHRᵀ. Inversion alone has R = −I, making that equation true for *any* H; it imposes no constraint on a rank-two Hessian. This argument assumes the potential due to the *other* charges is twice differentiable at the central nucleus; its singular self-potential is excluded.

## Choose coordinates and write a generic Hessian

Put four equivalent ligand sites at alternating corners of a cube of half-edge length a > 0:

$$
\begin{aligned}
A&=a( 1, 1, 1), & B&=a( 1,-1,-1),\\
C&=a(-1, 1,-1), & D&=a(-1,-1, 1).
\end{aligned}
$$

This tetrahedron has three C<sub>2</sub> (180°) axes along the Cartesian x, y, and z directions. Its four C<sub>3</sub> (120°) axes run through the centre and a vertex, along body diagonals such as (1, 1, 1)—**not** along Cartesian z in this orientation. These are the axes listed for the xyz setting in [Quanty's tetrahedral operation table](https://www.quanty.org/physics_chemistry/point_groups/td/orientation_xyz). We need only two of the C<sub>2</sub> operations and one C<sub>3</sub> operation: if they already force H to be proportional to the identity, every remaining operation automatically preserves it.

Because mixed second derivatives agree, H is symmetric. Before imposing symmetry, it has six independent entries:

$$
H=\begin{pmatrix}
H_{xx}&H_{xy}&H_{xz}\\
H_{xy}&H_{yy}&H_{yz}\\
H_{xz}&H_{yz}&H_{zz}
\end{pmatrix}.
$$

## First half-turn: eliminate two off-diagonal entries

A 180° rotation about x leaves x unchanged and reverses y and z:

$$
(x,y,z)\longmapsto(x,-y,-z),\qquad
R_x=\begin{pmatrix}1&0&0\\0&-1&0\\0&0&-1\end{pmatrix}.
$$

This is a symmetry of our tetrahedron: for example, it exchanges A ↔ B and C ↔ D. Multiplying the generic Hessian on both sides gives

$$
R_xHR_x^{T}=\begin{pmatrix}
H_{xx}&-H_{xy}&-H_{xz}\\
-H_{xy}&H_{yy}&H_{yz}\\
-H_{xz}&H_{yz}&H_{zz}
\end{pmatrix}.
$$

The original and transformed matrices must be equal **entry by entry**. For any diagonal sign-change matrix R = diag(s<sub>x</sub>, s<sub>y</sub>, s<sub>z</sub>), the transformed ij entry is s<sub>i</sub>s<sub>j</sub>H<sub>ij</sub>: an entry changes sign when exactly one of its indices is reversed. Here H<sub>xy</sub> = −H<sub>xy</sub> and H<sub>xz</sub> = −H<sub>xz</sub>, so both entries are zero. The yz entry survives because *both* its indices are reversed. We now have

$$
H=\begin{pmatrix}H_{xx}&0&0\\0&H_{yy}&H_{yz}\\0&H_{yz}&H_{zz}\end{pmatrix}.
$$

## Second half-turn: eliminate the remaining off-diagonal entry

For a 180° rotation about y, the coordinate map and matrix are

$$
(x,y,z)\longmapsto(-x,y,-z),\qquad
R_y=\begin{pmatrix}-1&0&0\\0&1&0\\0&0&-1\end{pmatrix}.
$$

Applying this operation to the *already reduced* Hessian gives

$$
R_yHR_y^{T}=\begin{pmatrix}
H_{xx}&0&0\\
0&H_{yy}&-H_{yz}\\
0&-H_{yz}&H_{zz}
\end{pmatrix}.
$$

Equality now requires H<sub>yz</sub> = −H<sub>yz</sub>, hence H<sub>yz</sub> = 0. The two perpendicular half-turns have proved that H is diagonal. They have **not** yet proved that its diagonal entries are equal:

$$
H=\operatorname{diag}(H_{xx},H_{yy},H_{zz}).
$$

## A threefold rotation: equate the diagonal entries

Take the C<sub>3</sub> axis through A = a(1, 1, 1). One of its 120° rotations cyclically sends the coordinate basis vectors as follows:

$$
\mathbf e_x\longmapsto\mathbf e_z,\qquad
\mathbf e_y\longmapsto\mathbf e_x,\qquad
\mathbf e_z\longmapsto\mathbf e_y.
$$

The images of the basis vectors are the **columns** of the matrix. Thus the rotation takes (x, y, z) to (y, z, x) and has matrix

$$
P=\begin{pmatrix}
0&1&0\\
0&0&1\\
1&0&0
\end{pmatrix}.
$$

We can check that this really is the desired tetrahedral operation: P(1, 1, 1)ᵀ = (1, 1, 1)ᵀ, while it cycles the other sites B → D → C → B. Also P³ = I and det P = +1, so it is a proper rotation of order three rather than a reflection. The inverse cycle would work equally well.

Applying P to the diagonal Hessian simply cycles its entries:

$$
P\begin{pmatrix}H_{xx}&0&0\\0&H_{yy}&0\\0&0&H_{zz}\end{pmatrix}P^{T}
=\begin{pmatrix}H_{yy}&0&0\\0&H_{zz}&0\\0&0&H_{xx}\end{pmatrix}.
$$

Entry-by-entry equality therefore gives H<sub>xx</sub> = H<sub>yy</sub> = H<sub>zz</sub> ≡ λ. Every allowed full Hessian has the form

$$
\boxed{H=\lambda I.}
$$

## From the full Hessian to the EFG

The quadrupolar EFG is the **traceless** part of the Hessian, not necessarily the entire Hessian. With the usual convention,

$$
V^{\mathrm{EFG}}=H-\frac{\operatorname{Tr}H}{3}I.
$$

Since H = λI has trace 3λ,

$$
\boxed{V^{\mathrm{EFG}}=\lambda I-\lambda I=0.}
$$

Thus every EFG component vanishes, including the largest principal-axis component V<sub>ZZ</sub>. The asymmetry parameter η = (V<sub>XX</sub> − V<sub>YY</sub>)/V<sub>ZZ</sub> is **undefined**, not zero, because its denominator is zero. The [EFG definition used in quadrupolar spectroscopy](https://pubs.acs.org/doi/10.1021/acs.chemmater.4c00883) likewise subtracts the isotropic trace.

Does H = λI also imply H = 0? Only if the observation point is charge-free. Poisson's equation gives

$$
\operatorname{Tr}H=\nabla^2\Phi(0)=-\frac{\rho(0)}{\varepsilon_0},
\qquad
\lambda=-\frac{\rho(0)}{3\varepsilon_0}.
$$

If ρ(0) = 0, as in an ideal point-charge model whose charges are all away from the origin, the **full** Hessian is zero. If diffuse electron density reaches the central nucleus, the full Hessian may instead be a nonzero isotropic matrix. Either way, its traceless EFG is exactly zero. This distinction matters for calculations with spatial atomic orbitals.

## What the result does—and does not—assume

The proof applies to any sufficiently smooth potential whose **total source distribution** is invariant under T<sub>d</sub>. Four equivalent radial p-orbital *densities* satisfy this condition: a tetrahedral rotation carries each bond-directed cloud onto the corresponding cloud at another ligand. Their sum need not be spherical. It only needs tetrahedral symmetry, which is enough to forbid a rank-two traceless curvature at the centre.

As a direct check, let û<sub>i</sub> point from the centre to each of the four sites above. The tetrahedral directions obey

$$
\sum_{i=1}^{4}\hat{\mathbf u}_i\hat{\mathbf u}_i^{T}=\frac43 I.
$$

Identical radial p clouds have the same axially symmetric traceless contribution, with a coefficient β that depends on their size and charge:

$$
V_i=\beta\bigl(\hat{\mathbf u}_i\hat{\mathbf u}_i^{T}-I/3\bigr).
$$

Consequently, their sum is zero even though **each** contribution is generally nonzero. This direct cancellation agrees with, but is less general than, the symmetry proof above.

Finally, tetrahedrally arranged *nuclei* alone are not the entire premise. If unequal orbital occupations, a structural distortion, an external field, or another effect makes the **total electronic density** less symmetric than T<sub>d</sub>, the conclusion need not hold. Exact T<sub>d</sub> symmetry of the total charge distribution is the sufficient condition; the point-charge picture is only one way to realize it.
