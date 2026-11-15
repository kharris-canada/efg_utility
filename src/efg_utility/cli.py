"""Command-line interface for :mod:`efg_utility`."""

from __future__ import annotations

import argparse
import sys
from collections.abc import Sequence
from pathlib import Path

from .formatting import format_orbital_result, format_result
from .input import InputFileError, read_orbital_input, read_point_charges
from .models import HydrogenicOrbital, OrbitalCombination, OrbitalCombinationTerm
from .orbital import calculate_orbital_combinations
from .point_charge import CalculationError, calculate_efg


def _run_calculate(arguments: argparse.Namespace) -> int:
    try:
        charges = read_point_charges(arguments.input_file)
        result = calculate_efg(charges)
    except (InputFileError, CalculationError) as error:
        print(f"efg_utility: error: {error}", file=sys.stderr)
        return 2

    print(format_result(result))
    return 0


def _run_orbital(arguments: argparse.Namespace) -> int:
    try:
        orbitals, combinations = read_orbital_input(arguments.input_file)
        combinations.extend(
            OrbitalCombination(
                terms=(OrbitalCombinationTerm(
                    orbital=HydrogenicOrbital(
                        center_m=orbital.center_m,
                        orbital_type=orbital.orbital_type,
                        radius_m=orbital.radius_m,
                        axis=orbital.axis,
                    ),
                    coefficient=1.0,
                ),),
                occupation=orbital.occupation,
                name=f"independent-{index}",
            )
            for index, orbital in enumerate(orbitals, start=1)
        )
        result = calculate_orbital_combinations(combinations)
    except (InputFileError, CalculationError) as error:
        print(f"efg_utility: error: {error}", file=sys.stderr)
        return 2

    print(format_orbital_result(result))
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="efg_utility",
        description="Calculate electric field gradient parameters.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    calculate_parser = subparsers.add_parser(
        "calculate",
        help="calculate the EFG at the origin from point charges",
    )
    calculate_parser.add_argument(
        "input_file",
        type=Path,
        metavar="INPUT_FILE",
        help="five-column point-charge input file",
    )
    calculate_parser.set_defaults(handler=_run_calculate)

    orbital_parser = subparsers.add_parser(
        "orbital",
        help="calculate the field at the origin from hydrogenic orbitals/SALCs",
    )
    orbital_parser.add_argument(
        "input_file",
        type=Path,
        metavar="INPUT_FILE",
        help="eight-column orbital rows and/or ten-column grouped SALC AO terms",
    )
    orbital_parser.set_defaults(handler=_run_orbital)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    arguments = parser.parse_args(argv)
    return int(arguments.handler(arguments))
