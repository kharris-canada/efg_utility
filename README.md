# efg_utility

`efg_utility` calculates electrostatic potential derivatives at the origin
from point charges or displaced hydrogenic orbitals. Point-charge calculations
report the signed principal-axis component `VZZ` and asymmetry parameter `eta`.
Orbital calculations additionally distinguish the full potential Hessian from
its conventional traceless EFG part.

## Installation

```bash
python -m pip install .
```

For development and testing:

```bash
python -m pip install -e '.[dev]'
pytest
```

## Point-charge input

The input is a whitespace-delimited text file with five fields per data row:

```text
atom_name  x_angstrom  y_angstrom  z_angstrom  integer_charge
```

Coordinates are relative to the target nucleus at the origin and are measured
in angstroms. The charge is an integer multiple of the elementary charge. The
atom name is accepted for readability and discarded. Blank lines and `#`
comments are allowed, including comments after a data row.

For example:

```text
# Six +1 charges in an octahedral arrangement
X  1.0  0.0  0.0  +1
X -1.0  0.0  0.0  +1
X  0.0  1.0  0.0  +1
X  0.0 -1.0  0.0  +1
X  0.0  0.0  1.0  +1
X  0.0  0.0 -1.0  +1
```

## Usage

After installation:

```bash
efg_utility calculate charges.txt
```

The package can also be run directly as a module:

```bash
python -m efg_utility calculate charges.txt
```

Output is expressed in SI units:

```text
VZZ = +2.879929...e+21 V/m^2
eta = 0
```

`VZZ` is the signed eigenvalue with the greatest absolute magnitude. The other
two principal components are ordered so that `|VZZ| >= |VYY| >= |VXX|`, and
`eta = (VXX - VYY) / VZZ`. If the total EFG is zero within floating-point
cancellation tolerance, `VZZ` is reported as zero and `eta` as `undefined`.

A nonzero point charge at the origin is singular and is rejected. Zero-charge
rows are permitted and make no contribution.

The command-line interface uses subcommands so later input-generation and
additional orbital types can be added independently.

## Hydrogenic orbitals

Use the `orbital` subcommand with eight-column independent-orbital rows:

```text
name  x_A  y_A  z_A  orbital  radius_A  zeff  occupation
```

The center coordinates and explicit radius are in angstroms. Supported orbital
types are `1s`, `2s`, `2px`, `2py`, and `2pz`. The real p-orbital suffix gives
the lobe axis in the input coordinate system. `2p_x`, `2p_y`, and `2p_z` are
also accepted as aliases. An arbitrary real p axis can be written as a single
whitespace-free token such as `2p(1,1,0)`; its three finite, nonzero components
are normalized automatically. A bare `2p` is rejected because its orientation
would be ambiguous.

Use `-` for an unused radius or effective nuclear charge field:

- `- -` uses the Bohr radius.
- `0.40 -` uses an explicit radius of 0.40 angstrom.
- `- 2.0` uses `a = a0 / Zeff` with `Zeff = 2.0`.

Radius and `Zeff` cannot both be specified. Occupation is any finite,
nonnegative number and scales the negative electron density. Atom names are
accepted for readability and discarded.

For a single occupied 1s orbital centered at `(0, 0, 1)` angstrom:

```text
H  0.0  0.0  1.0  1s  -  -  1.0
```

For coherent AO combinations such as SALCs, use ten columns:

```text
name  x_A  y_A  z_A  orbital  radius_A  zeff  coefficient  group  occupation
```

Signed coefficients on rows with the same `group` label are amplitudes in one
linear combination. The backend evaluates the complete density
`occupation * |sum(c_i * chi_i)|^2`, including cross terms between different
AOs in that group. Distinct groups add as separate occupied orbitals; no
cross terms are added between them. The backend normalizes each combination
with its calculated AO overlap matrix. Eight-column rows continue to mean
independent orbital densities and do not interfere with one another.

Run:

```bash
efg_utility orbital examples/orbital_1s.txt
```

The output contains the electrostatic potential, the full Hessian and its
trace, the traceless EFG tensor, and PAS values derived from the traceless
tensor. The full Hessian generally has a nonzero trace because the displaced
electron density may be nonzero at the observation point. The traceless EFG
removes the isotropic Poisson contribution:

```text
potential = ... V
full Hessian (V/m^2):
[...]
[...]
[...]
trace = ... V/m^2
traceless EFG (V/m^2):
[...]
[...]
[...]
VZZ = ... V/m^2
eta = ...
```

### Recovering the point-charge limit

The file [`examples/orbital_1s_zeff4.txt`](examples/orbital_1s_zeff4.txt)
places one occupied 1s orbital at `(0, 0, 1)` angstrom with `Zeff = 4`:

```text
H  0.0  0.0  1.0  1s  -  4.0  1.0
```

Run it with:

```bash
efg_utility orbital examples/orbital_1s_zeff4.txt
```

This is an interesting validation case because `a = a0 / 4` is approximately
0.1323 angstrom, so the observation point is at `R/a = 7.56`. The electron
density is therefore compact compared with its distance from the origin and
should look almost like a point electron. The calculation gives:

```text
trace = +6.764058212208e+17 V/m^2
VZZ = -2.879376050505e+21 V/m^2
eta = 0
```

A point electron at the same position gives
`VZZ = -2.879929093734e+21 V/m^2`. The orbital result differs in magnitude by
only about 0.0192%, demonstrating that the analytic density calculation
correctly recovers the existing point-charge limit as `R/a` becomes large.

For a normalized 1s electron density with `x = R/a` and total charge
`q = -occupation * e`, the analytic potential is

```text
Phi(R) = k*q/a * [1 - (1 + x) exp(-2x)] / x.
```

At large `R/a`, the traceless tensor approaches the result for a point charge
`q` at the orbital center. At `R = 0`, the full Hessian is finite and isotropic,
while the traceless EFG is zero and `eta` is undefined.

## 2s and 2p orbitals

The n=2 implementation uses the standard normalized hydrogenic functions with
the same scale `a = a0 / Zeff`:

```text
|psi_2s|^2   = (2 - r/a)^2 exp(-r/a) / (32*pi*a^3)
|psi_2p_u|^2 = (u . r)^2 exp(-r/a) / (32*pi*a^5)
```

Here `u` is the x, y, or z unit vector selected by `2px`, `2py`, or `2pz`.
The 2s density is spherical, so its field is evaluated from its enclosed radial
charge. A real 2p density contains only monopole and quadrupole angular terms;
the code evaluates both terms analytically and differentiates their potential
to obtain the full Cartesian Hessian. This avoids spatial grids and remains
well behaved at the orbital center.

The example files place occupied 2s and 2pz orbitals one angstrom along z:

```bash
efg_utility orbital examples/orbital_2s.txt
efg_utility orbital examples/orbital_2pz.txt
```

These cases are useful because the 2s result must remain axially symmetric
about the displacement, whereas the 2pz result combines that displacement
with an intrinsically directional density. At an orbital center, 2s has a zero
traceless EFG, but 2p retains a finite quadrupolar EFG aligned with its lobe
axis. Far from either density, both recover the point-electron tensor.

The analytic n=2 tensors are tested against independent spherical
Gauss-Legendre integration at `R/a = 0.5, 1, 2, 5, 10`, as well as density
normalization, Poisson traces, axis rotations, center limits, and the
point-charge limit. Future orbital families can use the same input and result
interfaces through separate orbital-type backends.

## Point-group symmetry expansion

The standalone [`expand_symmetry.py`](expand_symmetry.py) script turns one or
more symmetry-unique charge or orbital rows into a complete EFG input file.
It applies every operation of the requested finite molecular point group and
removes duplicate members within each seed orbit:

```bash
python expand_symmetry.py \
    --group Oh \
    --mode charge \
    examples/symmetry_oh_seed.txt \
    -o octahedral_charges.txt

python -m efg_utility calculate octahedral_charges.txt
```

Without `-o`, generated rows are printed to standard output, so piping and
shell redirection also work:

```bash
python expand_symmetry.py --group Oh --mode charge seed.txt > charges.txt
```

For the seed

```text
X  0  0  1  -2
```

`Oh` produces the six positions `(+/-1,0,0)`, `(0,+/-1,0)`, and
`(0,0,+/-1)`. No special-position declaration is needed. Although `Oh` has 48
operations, operations in the seed's stabilizer produce duplicates that are
removed automatically.

The accepted finite Schoenflies groups are:

- `C1`, `Ci`, `Cs`
- `Cn`, `Cnv`, and `Cnh`
- even-order `Sn`
- `Dn`, `Dnh`, and `Dnd`
- `T`, `Th`, `Td`, `O`, `Oh`, `I`, and `Ih`

Names are case-insensitive and underscores are optional. Infinite groups such
as `Cinfv` cannot be enumerated and are rejected. The standard orientation has
the principal axis along z, a vertical mirror plane in xz, and a perpendicular
twofold axis along x. Tetrahedral and octahedral groups use their conventional
Cartesian orientation. Icosahedral groups have a fivefold axis along z and a
twofold axis in the xz plane.

For charge and s-orbital rows, only the position is transformed. For a 2p row,
both its center and its axis are transformed. General rotated axes are emitted
as `2p(ux,uy,uz)`, while axes equivalent to x, y, or z use `2px`, `2py`, or
`2pz`. Axes `u` and `-u` describe the same density and are deduplicated. This
also means that symmetry may generate several differently oriented p orbitals
at the same special-position center.

The default position deduplication tolerance is `1e-9` angstrom and can be
changed with `--tolerance`. Separate input rows remain independent, allowing
distinct orbitals at the same center.

For a four-vertex tetrahedron of real p orbitals pointing outward from the
origin, use a p axis along the seed vertex's body diagonal:

```text
X  1  1  1  2p(1,1,1)  -  6.0  0.333
```

Expanding this row with `--group Td --mode orbital` gives four centers and one
radial p orbital at each. A `2pz` seed at `(1,1,1)` has a different meaning:
its axis is the laboratory z direction. The threefold rotation about that
vertex leaves its center fixed but rotates the p axis among x, y, and z.
Its symmetry orbit therefore has 12 distinct orbital rows: three orientations
at each of four centers. The apparent repeated positions are separate orbital
densities, each retaining the seed occupation.

## Converting a SALC table to EFG orbital input

[`salc_to_efg.py`](salc_to_efg.py) selects one SALC from either Cartesian SALC
table and prints rows accepted by `efg_utility orbital`:

```bash
python salc_to_efg.py tetrahedral_salcs/salcs_cartesian.tsv t2-A-1 \
    --distance 1.5 --zeff 4
python salc_to_efg.py tetrahedral_salcs/salcs_cartesian.tsv t2-A-1 \
    --distance 1.5 --zeff 4 > t2_A_1_orbitals.txt
python -m efg_utility orbital t2_A_1_orbitals.txt
```

The distance `d` and effective nuclear charge are required. An optional
`--occupation` (default `2`, appropriate for a closed-shell SALC) sets the
total occupation of each selected SALC
partner. Use a full id such as `t2-A-1` for one partner, omit the final number
(`t2-A`) for every partner at that irrep and energy rank, or use `all` to emit
the entire table. Thus, for the closed-shell sum over every listed SALC:

```bash
python salc_to_efg.py tetrahedral_salcs/salcs_cartesian.tsv all \
    --distance 1.5 --zeff 4 > tetrahedral_salcs.txt
python -m efg_utility orbital tetrahedral_salcs.txt
```

The helper retains each signed coefficient and emits ten-column rows grouped
by SALC id. Thus constructive/destructive interference within each selected
SALC is included, and the backend renormalizes coefficients for overlap
between the finite-size hydrogenic AOs. To calculate a total from multiple
SALCs, append their rows to one input file (using their group identifiers) and
run the orbital command once. Each SALC partner is assigned occupation 2 by
default; `--occupation` can override this when needed. For example:

```bash
python salc_to_efg.py square_antiprism_salcs/salcs_cartesian.tsv e1-A \
    --distance 2.0 --zeff 6 > e1_A_orbitals.txt
python -m efg_utility orbital e1_A_orbitals.txt
```
