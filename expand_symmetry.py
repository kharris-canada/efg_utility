#!/usr/bin/env python3
"""Expand EFG input records using a finite molecular point group."""

from __future__ import annotations

import argparse
import math
import re
import sys
from collections import deque
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

import numpy as np

DEFAULT_POSITION_TOLERANCE_ANGSTROM = 1.0e-9
_AXIS_TOLERANCE = 1.0e-10
_MATRIX_TOLERANCE = 2.0e-10
_MAX_N = 360
_INTEGER_PATTERN = re.compile(r"[+-]?[0-9]+")
_P_DESCRIPTOR_PATTERN = re.compile(
    r"2p\(([^,()]+),([^,()]+),([^,()]+)\)",
    re.IGNORECASE,
)
_P_SHORTHANDS = {
    "2px": np.array((1.0, 0.0, 0.0)),
    "2py": np.array((0.0, 1.0, 0.0)),
    "2pz": np.array((0.0, 0.0, 1.0)),
    "2p_x": np.array((1.0, 0.0, 0.0)),
    "2p_y": np.array((0.0, 1.0, 0.0)),
    "2p_z": np.array((0.0, 0.0, 1.0)),
}


class SymmetryExpansionError(ValueError):
    """Raised for invalid point groups, source rows, or output requests."""


@dataclass(frozen=True, slots=True)
class SeedRecord:
    """One validated charge or orbital source row."""

    name: str
    position: np.ndarray
    payload: tuple[str, ...]
    orbital_axis: np.ndarray | None = None


def _rotation(axis: Sequence[float], angle: float) -> np.ndarray:
    vector = np.asarray(axis, dtype=np.float64)
    vector /= np.linalg.norm(vector)
    x, y, z = vector
    cross = np.array(((0.0, -z, y), (z, 0.0, -x), (-y, x, 0.0)))
    outer = np.outer(vector, vector)
    return (
        math.cos(angle) * np.identity(3)
        + math.sin(angle) * cross
        + (1.0 - math.cos(angle)) * outer
    )


def _reflection(normal: Sequence[float]) -> np.ndarray:
    vector = np.asarray(normal, dtype=np.float64)
    vector /= np.linalg.norm(vector)
    return np.identity(3) - 2.0 * np.outer(vector, vector)


def _equivalent_matrix(first: np.ndarray, second: np.ndarray) -> bool:
    return bool(np.max(np.abs(first - second)) <= _MATRIX_TOLERANCE)


def _group_closure(
    generators: Sequence[np.ndarray],
    *,
    expected_order: int,
) -> tuple[np.ndarray, ...]:
    identity = np.identity(3)
    operations = [identity]
    pending: deque[np.ndarray] = deque((identity,))

    while pending:
        current = pending.popleft()
        for generator in generators:
            candidate = current @ generator
            candidate[np.abs(candidate) < 1.0e-14] = 0.0
            if any(_equivalent_matrix(candidate, item) for item in operations):
                continue
            operations.append(candidate)
            pending.append(candidate)
            if len(operations) > expected_order:
                raise SymmetryExpansionError(
                    "internal point-group construction exceeded its expected order"
                )

    if len(operations) != expected_order:
        raise SymmetryExpansionError(
            "internal point-group construction produced "
            f"{len(operations)} operations; expected {expected_order}"
        )
    return tuple(operations)


def _polyhedral_group(name: str) -> tuple[np.ndarray, ...]:
    inversion = -np.identity(3)
    tetrahedral_generators = (
        _rotation((1.0, 1.0, 1.0), 2.0 * math.pi / 3.0),
        _rotation((1.0, 0.0, 0.0), math.pi),
    )
    octahedral_generators = (
        _rotation((0.0, 0.0, 1.0), math.pi / 2.0),
        _rotation((1.0, 1.0, 1.0), 2.0 * math.pi / 3.0),
    )

    if name == "t":
        return _group_closure(tetrahedral_generators, expected_order=12)
    if name == "th":
        return _group_closure(
            (*tetrahedral_generators, inversion), expected_order=24
        )
    if name == "td":
        return _group_closure(
            (*tetrahedral_generators, _reflection((1.0, -1.0, 0.0))),
            expected_order=24,
        )
    if name == "o":
        return _group_closure(octahedral_generators, expected_order=24)
    if name == "oh":
        return _group_closure(
            (*octahedral_generators, inversion), expected_order=48
        )

    if name in ("i", "ih"):
        fivefold = _rotation((0.0, 0.0, 1.0), 2.0 * math.pi / 5.0)
        cos_beta = 1.0 / (2.0 * math.sin(math.pi / 5.0))
        twofold_axis = (math.sqrt(1.0 - cos_beta * cos_beta), 0.0, cos_beta)
        generators: tuple[np.ndarray, ...] = (
            fivefold,
            _rotation(twofold_axis, math.pi),
        )
        expected_order = 60
        if name == "ih":
            generators = (*generators, inversion)
            expected_order = 120
        return _group_closure(generators, expected_order=expected_order)

    raise AssertionError(f"unknown polyhedral group {name}")


def point_group_operations(group_name: str) -> tuple[np.ndarray, ...]:
    """Return all Cartesian operations for a finite Schoenflies point group."""

    name = group_name.replace("_", "").casefold()
    if name == "c1":
        return (np.identity(3),)
    if name == "ci":
        return (np.identity(3), -np.identity(3))
    if name == "cs":
        return (np.identity(3), _reflection((0.0, 1.0, 0.0)))
    if name in {"t", "th", "td", "o", "oh", "i", "ih"}:
        return _polyhedral_group(name)
    if "inf" in name or "∞" in name:
        raise SymmetryExpansionError(
            f"point group {group_name!r} is infinite and cannot be enumerated"
        )

    cyclic_match = re.fullmatch(r"c([1-9][0-9]*)([vh]?)", name)
    if cyclic_match is not None:
        order = _validated_family_order(cyclic_match.group(1), group_name)
        suffix = cyclic_match.group(2)
        generators = [_rotation((0.0, 0.0, 1.0), 2.0 * math.pi / order)]
        expected_order = order
        if suffix == "v":
            generators.append(_reflection((0.0, 1.0, 0.0)))
            expected_order *= 2
        elif suffix == "h":
            generators.append(_reflection((0.0, 0.0, 1.0)))
            expected_order *= 2
        return _group_closure(generators, expected_order=expected_order)

    improper_match = re.fullmatch(r"s([1-9][0-9]*)", name)
    if improper_match is not None:
        order = _validated_family_order(improper_match.group(1), group_name)
        if order % 2 != 0:
            raise SymmetryExpansionError(
                "S_n point groups use an even n (for example S4 or S6)"
            )
        improper_rotation = (
            _reflection((0.0, 0.0, 1.0))
            @ _rotation((0.0, 0.0, 1.0), 2.0 * math.pi / order)
        )
        return _group_closure((improper_rotation,), expected_order=order)

    dihedral_match = re.fullmatch(r"d([1-9][0-9]*)([hd]?)", name)
    if dihedral_match is not None:
        order = _validated_family_order(dihedral_match.group(1), group_name)
        if order < 2:
            raise SymmetryExpansionError("D_n point groups require n >= 2")
        suffix = dihedral_match.group(2)
        generators = [
            _rotation((0.0, 0.0, 1.0), 2.0 * math.pi / order),
            _rotation((1.0, 0.0, 0.0), math.pi),
        ]
        expected_order = 2 * order
        if suffix == "h":
            generators.append(_reflection((0.0, 0.0, 1.0)))
            expected_order *= 2
        elif suffix == "d":
            diagonal_angle = math.pi / (2.0 * order)
            generators.append(
                _reflection(
                    (-math.sin(diagonal_angle), math.cos(diagonal_angle), 0.0)
                )
            )
            expected_order *= 2
        return _group_closure(generators, expected_order=expected_order)

    raise SymmetryExpansionError(
        f"unsupported point group {group_name!r}; use a finite Schoenflies name"
    )


def _validated_family_order(text: str, group_name: str) -> int:
    order = int(text)
    if order > _MAX_N:
        raise SymmetryExpansionError(
            f"point group {group_name!r} exceeds the supported maximum n={_MAX_N}"
        )
    return order


def _parse_finite(text: str, *, field: str, path: Path, line: int) -> float:
    try:
        value = float(text)
    except ValueError:
        raise SymmetryExpansionError(
            f"{path}:{line}: {field} must be a number, got {text!r}"
        ) from None
    if not math.isfinite(value):
        raise SymmetryExpansionError(f"{path}:{line}: {field} must be finite")
    return value


def _parse_p_axis(
    descriptor: str,
    *,
    path: Path,
    line: int,
) -> np.ndarray | None:
    normalized = descriptor.casefold()
    if normalized in ("1s", "2s"):
        return None
    if normalized in _P_SHORTHANDS:
        return _P_SHORTHANDS[normalized].copy()
    if normalized == "2p":
        raise SymmetryExpansionError(
            f"{path}:{line}: orbital type '2p' requires an axis"
        )
    match = _P_DESCRIPTOR_PATTERN.fullmatch(normalized)
    if match is None:
        raise SymmetryExpansionError(
            f"{path}:{line}: unsupported orbital type {descriptor!r}"
        )
    values = np.array(
        [
            _parse_finite(item, field="2p axis component", path=path, line=line)
            for item in match.groups()
        ]
    )
    magnitude = float(np.linalg.norm(values))
    if not math.isfinite(magnitude) or magnitude == 0.0:
        raise SymmetryExpansionError(f"{path}:{line}: a 2p axis must be nonzero")
    return values / magnitude


def read_seed_records(path: str | Path, mode: str) -> list[SeedRecord]:
    """Read and validate charge or orbital seed rows."""

    input_path = Path(path)
    expected_fields = 5 if mode == "charge" else 8
    records: list[SeedRecord] = []
    try:
        with input_path.open("r", encoding="utf-8") as stream:
            for line_number, raw_line in enumerate(stream, start=1):
                content = raw_line.partition("#")[0].strip()
                if not content:
                    continue
                fields = content.split()
                if len(fields) != expected_fields:
                    raise SymmetryExpansionError(
                        f"{input_path}:{line_number}: expected {expected_fields} "
                        f"fields, found {len(fields)}"
                    )
                coordinates = np.array(
                    [
                        _parse_finite(
                            value,
                            field=f"{'xyz'[index]} coordinate",
                            path=input_path,
                            line=line_number,
                        )
                        for index, value in enumerate(fields[1:4])
                    ]
                )
                payload = tuple(fields[4:])
                orbital_axis = None
                if mode == "charge":
                    if _INTEGER_PATTERN.fullmatch(payload[0]) is None:
                        raise SymmetryExpansionError(
                            f"{input_path}:{line_number}: charge must be an integer, "
                            f"got {payload[0]!r}"
                        )
                else:
                    orbital_axis = _parse_p_axis(
                        payload[0], path=input_path, line=line_number
                    )
                    _validate_orbital_payload(payload, input_path, line_number)
                records.append(
                    SeedRecord(fields[0], coordinates, payload, orbital_axis)
                )
    except SymmetryExpansionError:
        raise
    except OSError as error:
        raise SymmetryExpansionError(
            f"{input_path}: {error.strerror or str(error)}"
        ) from None
    return records


def _validate_orbital_payload(
    payload: tuple[str, ...], path: Path, line: int
) -> None:
    _orbital_type, radius_text, zeff_text, occupation_text = payload
    if radius_text != "-" and zeff_text != "-":
        raise SymmetryExpansionError(
            f"{path}:{line}: radius and Zeff are mutually exclusive; use '-' for one"
        )
    if radius_text != "-":
        radius = _parse_finite(radius_text, field="orbital radius", path=path, line=line)
        if radius <= 0.0:
            raise SymmetryExpansionError(
                f"{path}:{line}: orbital radius must be positive"
            )
    if zeff_text != "-":
        zeff = _parse_finite(zeff_text, field="Zeff", path=path, line=line)
        if zeff <= 0.0:
            raise SymmetryExpansionError(f"{path}:{line}: Zeff must be positive")
    occupation = _parse_finite(
        occupation_text, field="occupation", path=path, line=line
    )
    if occupation < 0.0:
        raise SymmetryExpansionError(
            f"{path}:{line}: occupation must be nonnegative"
        )


def _canonical_axis(axis: np.ndarray) -> np.ndarray:
    result = axis / np.linalg.norm(axis)
    for value in result:
        if abs(value) > _AXIS_TOLERANCE:
            if value < 0.0:
                result = -result
            break
    result[np.abs(result) < _AXIS_TOLERANCE] = 0.0
    return result


def _equivalent_source(
    first_position: np.ndarray,
    first_axis: np.ndarray | None,
    second_position: np.ndarray,
    second_axis: np.ndarray | None,
    position_tolerance: float,
) -> bool:
    if np.max(np.abs(first_position - second_position)) > position_tolerance:
        return False
    if first_axis is None or second_axis is None:
        return first_axis is None and second_axis is None
    return bool(
        min(
            np.max(np.abs(first_axis - second_axis)),
            np.max(np.abs(first_axis + second_axis)),
        )
        <= _AXIS_TOLERANCE
    )


def expand_record(
    record: SeedRecord,
    operations: Sequence[np.ndarray],
    *,
    position_tolerance: float,
) -> list[tuple[np.ndarray, np.ndarray | None]]:
    """Expand one record into its unique symmetry orbit."""

    expanded: list[tuple[np.ndarray, np.ndarray | None]] = []
    for operation in operations:
        position = operation @ record.position
        axis = (
            _canonical_axis(operation @ record.orbital_axis)
            if record.orbital_axis is not None
            else None
        )
        if any(
            _equivalent_source(
                position,
                axis,
                known_position,
                known_axis,
                position_tolerance,
            )
            for known_position, known_axis in expanded
        ):
            continue
        position[np.abs(position) < position_tolerance] = 0.0
        expanded.append((position, axis))

    def sort_key(item: tuple[np.ndarray, np.ndarray | None]) -> tuple[int | float, ...]:
        position, axis = item
        axis_values = axis if axis is not None else np.zeros(3)
        # Group positions that differ only by matrix roundoff before ordering
        # their distinct p axes. Raw floats can interleave visually identical
        # centers when several orbital orientations share a special position.
        position_key = tuple(round(float(value) / position_tolerance) for value in position)
        axis_key = tuple(float(value) for value in axis_values)
        return (*position_key, *axis_key)

    return sorted(expanded, key=sort_key)


def _format_number(value: float) -> str:
    if value == 0.0:
        return "0"
    return f"{value:.12g}"


def format_p_descriptor(axis: np.ndarray) -> str:
    """Format an unoriented normalized p axis as a valid orbital token."""

    canonical = _canonical_axis(axis.copy())
    for name, candidate in (
        ("2px", np.array((1.0, 0.0, 0.0))),
        ("2py", np.array((0.0, 1.0, 0.0))),
        ("2pz", np.array((0.0, 0.0, 1.0))),
    ):
        if np.max(np.abs(canonical - candidate)) <= _AXIS_TOLERANCE:
            return name
    components = ",".join(_format_number(float(value)) for value in canonical)
    return f"2p({components})"


def format_expanded_records(
    records: Sequence[SeedRecord],
    operations: Sequence[np.ndarray],
    *,
    position_tolerance: float,
) -> str:
    """Return complete EFG input text for all seed records."""

    lines: list[str] = []
    for record in records:
        for position, axis in expand_record(
            record, operations, position_tolerance=position_tolerance
        ):
            payload = list(record.payload)
            if axis is not None:
                payload[0] = format_p_descriptor(axis)
            fields = [
                record.name,
                *(_format_number(float(value)) for value in position),
                *payload,
            ]
            lines.append("  ".join(fields))
    return "\n".join(lines) + ("\n" if lines else "")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Expand EFG source rows using a finite molecular point group."
    )
    parser.add_argument("input_file", type=Path, metavar="INPUT_FILE")
    parser.add_argument("--group", required=True, help="finite Schoenflies point group")
    parser.add_argument(
        "--mode", required=True, choices=("charge", "orbital"), help="input row grammar"
    )
    parser.add_argument("-o", "--output", type=Path, help="output file (default: stdout)")
    parser.add_argument(
        "--tolerance",
        type=float,
        default=DEFAULT_POSITION_TOLERANCE_ANGSTROM,
        metavar="ANGSTROM",
        help="absolute position deduplication tolerance (default: 1e-9)",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    arguments = build_parser().parse_args(argv)
    try:
        if not math.isfinite(arguments.tolerance) or arguments.tolerance <= 0.0:
            raise SymmetryExpansionError("tolerance must be finite and positive")
        operations = point_group_operations(arguments.group)
        records = read_seed_records(arguments.input_file, arguments.mode)
        output = format_expanded_records(
            records,
            operations,
            position_tolerance=arguments.tolerance,
        )
        if arguments.output is None:
            sys.stdout.write(output)
        else:
            if arguments.output.resolve() == arguments.input_file.resolve():
                raise SymmetryExpansionError("input and output files must be different")
            arguments.output.write_text(output, encoding="utf-8")
    except (SymmetryExpansionError, OSError) as error:
        detail = error.strerror if isinstance(error, OSError) else str(error)
        print(f"expand_symmetry.py: error: {detail or error}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
