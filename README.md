# efg_utility

`efg_utility` calculates the electric field gradient (EFG) at the origin from
a collection of point charges. It reports the signed principal-axis component
`VZZ` and the asymmetry parameter `eta`.

## Installation

```bash
python -m pip install .
```

For development and testing:

```bash
python -m pip install -e '.[dev]'
pytest
```

## Input format

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
orbital-based calculations can be added independently. Those modes are not
implemented in this version.
