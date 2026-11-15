#!/usr/bin/env python3
"""Select a Cartesian SALC and print orbital input for ``efg_utility``.

Examples::

    python salc_to_efg.py tetrahedral_salcs/salcs_cartesian.tsv t2-A-1 \
        --distance 1.5 --zeff 4
    python salc_to_efg.py square_antiprism_salcs/salcs_cartesian.tsv e1-A \
        --distance 2.0 --zeff 6 --occupation 2

A complete identifier selects one partner; omitting its final partner number
selects every partner at that symmetry and energy rank. Signed AO amplitudes
are retained, and rows sharing a SALC identifier are treated coherently by the
EFG backend.
"""

from __future__ import annotations

import argparse
import ast
import math
import operator
import re
import sys
from collections import OrderedDict
from dataclasses import dataclass
from pathlib import Path
from typing import Callable


_NUMBER = r"[+-]?(?:(?:\d+(?:\.\d*)?)|(?:\.\d+))(?:[eE][+-]?\d+)?"
_SALC_ID = re.compile(r"(?P<irrep>[a-z][a-z0-9]*)-(?P<energy>[A-Z])-(?P<partner>[1-9][0-9]*)", re.IGNORECASE)
_SELECTION = re.compile(r"(?P<irrep>[a-z][a-z0-9]*)-(?P<energy>[A-Z])(?:-(?P<partner>[1-9][0-9]*))?", re.IGNORECASE)
_AO_CONTRIBUTION = re.compile(rf"(?P<coefficient>{_NUMBER})\*(?P<ao>2s|2px|2py|2pz)", re.IGNORECASE)
_SUBSCRIPT_DIGITS = str.maketrans("₀₁₂₃₄₅₆₇₈₉", "0123456789")


class SalcInputError(ValueError):
    """Raised for malformed SALC tables, selections, or numeric arguments."""


@dataclass(frozen=True, slots=True)
class Contribution:
    salc_id: str
    coefficient: float
    ao: str
    coordinate_expressions: tuple[str, str, str]
    line_number: int


def _parse_salc_id(text: str, *, path: Path, line_number: int) -> tuple[str, str, int]:
    normalized_text = text.translate(_SUBSCRIPT_DIGITS)
    match = _SALC_ID.fullmatch(normalized_text)
    if match is None:
        raise SalcInputError(
            f"{path}:{line_number}: invalid SALC id {text!r}; expected irrep-energy-partner"
        )
    return (match.group("irrep").lower(), match.group("energy").upper(),
            int(match.group("partner")))


def read_table(path: str | Path) -> list[Contribution]:
    """Read the five-column Cartesian SALC TSV format."""
    table_path = Path(path)
    contributions: list[Contribution] = []
    try:
        with table_path.open("r", encoding="utf-8") as stream:
            for line_number, raw_line in enumerate(stream, start=1):
                line = raw_line.strip()
                if not line or line.startswith("#"):
                    continue
                fields = line.split("\t")
                if len(fields) != 5:
                    raise SalcInputError(
                        f"{table_path}:{line_number}: expected 5 tab-separated columns, found {len(fields)}"
                    )
                salc_id, ao_contribution, *coordinate_expressions = fields
                _parse_salc_id(salc_id, path=table_path, line_number=line_number)
                match = _AO_CONTRIBUTION.fullmatch(ao_contribution)
                if match is None:
                    raise SalcInputError(
                        f"{table_path}:{line_number}: invalid AO contribution {ao_contribution!r}"
                    )
                coefficient = float(match.group("coefficient"))
                if not math.isfinite(coefficient):
                    raise SalcInputError(
                        f"{table_path}:{line_number}: AO coefficient must be finite"
                    )
                contributions.append(Contribution(
                    salc_id=salc_id,
                    coefficient=coefficient,
                    ao=match.group("ao").lower(),
                    coordinate_expressions=(
                        coordinate_expressions[0],
                        coordinate_expressions[1],
                        coordinate_expressions[2],
                    ),
                    line_number=line_number,
                ))
    except SalcInputError:
        raise
    except OSError as error:
        detail = error.strerror or str(error)
        raise SalcInputError(f"{table_path}: {detail}") from None
    if not contributions:
        raise SalcInputError(f"{table_path}: no SALC contribution rows found")
    return contributions


def _evaluate_expression(expression: str, distance: float) -> float:
    """Evaluate a coordinate expression with only d, sqrt, and arithmetic."""
    try:
        tree = ast.parse(expression.replace("^", "**"), mode="eval")
    except SyntaxError:
        raise SalcInputError(f"invalid coordinate expression {expression!r}") from None

    binary_operators: dict[type[ast.operator], Callable[[float, float], float]] = {
        ast.Add: operator.add,
        ast.Sub: operator.sub,
        ast.Mult: operator.mul,
        ast.Div: operator.truediv,
        ast.Pow: operator.pow,
    }
    unary_operators: dict[type[ast.unaryop], Callable[[float], float]] = {
        ast.UAdd: operator.pos,
        ast.USub: operator.neg,
    }

    def evaluate(node: ast.AST) -> float:
        if isinstance(node, ast.Expression):
            return evaluate(node.body)
        if isinstance(node, ast.Constant) and isinstance(node.value, (int, float)):
            return float(node.value)
        if isinstance(node, ast.Name) and node.id == "d":
            return distance
        if isinstance(node, ast.BinOp) and type(node.op) in binary_operators:
            left = evaluate(node.left)
            right = evaluate(node.right)
            if isinstance(node.op, ast.Pow) and abs(right) > 8.0:
                raise SalcInputError("coordinate expression exponent is out of range")
            return binary_operators[type(node.op)](left, right)
        if isinstance(node, ast.UnaryOp) and type(node.op) in unary_operators:
            return unary_operators[type(node.op)](evaluate(node.operand))
        if (isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
                and node.func.id == "sqrt" and len(node.args) == 1
                and not node.keywords):
            return math.sqrt(evaluate(node.args[0]))
        raise SalcInputError(f"unsupported coordinate expression {expression!r}")

    try:
        result = evaluate(tree)
    except (ArithmeticError, OverflowError, ValueError):
        raise SalcInputError(f"could not evaluate coordinate expression {expression!r}") from None
    if not math.isfinite(result):
        raise SalcInputError(f"coordinate expression {expression!r} is not finite")
    return 0.0 if abs(result) < 5.0e-14 else result


def select_salcs(contributions: list[Contribution], selection: str) -> OrderedDict[str, list[Contribution]]:
    """Select one SALC id or all partners matching an irrep and energy rank."""
    by_id: OrderedDict[str, list[Contribution]] = OrderedDict()
    for contribution in contributions:
        by_id.setdefault(contribution.salc_id, []).append(contribution)
    if selection.casefold() in {"all", "*"}:
        return by_id

    match = _SELECTION.fullmatch(selection.translate(_SUBSCRIPT_DIGITS))
    if match is None:
        raise SalcInputError(
            f"invalid selection {selection!r}; use irrep-energy or irrep-energy-partner, e.g. t2-A-1"
        )
    wanted_irrep = match.group("irrep").lower()
    wanted_energy = match.group("energy").upper()
    wanted_partner = match.group("partner")
    selected: OrderedDict[str, list[Contribution]] = OrderedDict()
    for salc_id, rows in by_id.items():
        irrep, energy, partner = _parse_salc_id(
            salc_id, path=Path("SALC table"), line_number=rows[0].line_number
        )
        if irrep != wanted_irrep or energy != wanted_energy:
            continue
        if wanted_partner is not None and partner != int(wanted_partner):
            continue
        selected[salc_id] = rows
    if not selected:
        raise SalcInputError(f"no SALCs found for selection {selection!r}")
    return selected


def _format_number(value: float) -> str:
    if abs(value) < 5.0e-14:
        value = 0.0
    return f"{value:.12f}".rstrip("0").rstrip(".") or "0"


def _combination_rows(
    salc_id: str,
    rows: list[Contribution],
    *,
    distance: float,
    zeff: float,
    total_occupation: float,
) -> list[str]:
    output: list[str] = []
    for term_number, row in enumerate(rows, start=1):
        coordinate = tuple(
            _evaluate_expression(expr, distance)
            for expr in row.coordinate_expressions
        )
        fields = (
            f"{salc_id}_term{term_number:03d}",
            *(_format_number(value) for value in coordinate),
            row.ao,
            "-",
            f"{zeff:.15g}",
            f"{row.coefficient:.15g}",
            salc_id,
            f"{total_occupation:.15g}",
        )
        output.append(" ".join(fields))
    if not output:
        raise SalcInputError(f"{salc_id}: no AO contributions")
    return output


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Select one SALC (or all partners) and emit orbital input rows."
    )
    parser.add_argument("table", type=Path, help="Cartesian SALC TSV file")
    parser.add_argument(
        "selection", help="SALC id, irrep-energy (all partners), or 'all' for the full table"
    )
    parser.add_argument(
        "--distance", type=float, required=True, metavar="D_ANGSTROM",
        help="center-to-ligand distance in angstroms, substituted for d",
    )
    parser.add_argument(
        "--zeff", type=float, required=True,
        help="effective nuclear charge used for every emitted orbital",
    )
    parser.add_argument(
        "--occupation", type=float, default=2.0, metavar="N",
        help="total occupation assigned to each selected SALC partner (default: 2)",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    arguments = parser.parse_args(argv)
    if not math.isfinite(arguments.distance) or arguments.distance <= 0.0:
        parser.error("--distance must be a positive finite number")
    if not math.isfinite(arguments.zeff) or arguments.zeff <= 0.0:
        parser.error("--zeff must be a positive finite number")
    if not math.isfinite(arguments.occupation) or arguments.occupation < 0.0:
        parser.error("--occupation must be a nonnegative finite number")

    try:
        table = read_table(arguments.table)
        selected = select_salcs(table, arguments.selection)
        output = [
            "# Signed AO amplitudes are retained and interfere within each SALC group.",
            "# name x_A y_A z_A orbital radius_A zeff coefficient group occupation",
            "# Each selected SALC partner is assigned the requested total occupation.",
        ]
        for salc_id, rows in selected.items():
            output.extend(_combination_rows(
                salc_id,
                rows,
                distance=arguments.distance,
                zeff=arguments.zeff,
                total_occupation=arguments.occupation,
            ))
    except SalcInputError as error:
        print(f"salc_to_efg: error: {error}", file=sys.stderr)
        return 2

    print("\n".join(output))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
