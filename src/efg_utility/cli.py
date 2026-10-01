"""Command-line interface for :mod:`efg_utility`."""

from __future__ import annotations

import argparse
from collections.abc import Sequence
from pathlib import Path
import sys

from .formatting import format_result
from .input import InputFileError, read_point_charges
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
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    arguments = parser.parse_args(argv)
    return int(arguments.handler(arguments))
