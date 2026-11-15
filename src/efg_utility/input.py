"""Shared text input handling for spatial source records."""

from __future__ import annotations

import math
import re
from dataclasses import dataclass
from pathlib import Path

from .constants import ANGSTROM_IN_METRES, BOHR_RADIUS_M
from .models import (
    HydrogenicOrbital,
    OrbitalCombination,
    OrbitalCombinationTerm,
    PointCharge,
)
from .orbital import parse_orbital_descriptor, radius_from_zeff

_INTEGER_PATTERN = re.compile(r"[+-]?[0-9]+")


class InputFileError(ValueError):
    """Raised when an input file cannot be read or contains an invalid row."""


@dataclass(frozen=True, slots=True)
class SpatialRecord:
    """A position plus an uninterpreted source descriptor.

    The leading atom-name field is intentionally not retained. Keeping the
    final field lexical here allows future calculation modes to interpret
    orbital descriptors without changing coordinate parsing.
    """

    position_angstrom: tuple[float, float, float]
    source_spec: str
    line_number: int


@dataclass(frozen=True, slots=True)
class _SpatialRow:
    position_angstrom: tuple[float, float, float]
    payload: tuple[str, ...]
    line_number: int


def _row_error(path: Path, line_number: int, message: str) -> InputFileError:
    return InputFileError(f"{path}:{line_number}: {message}")


def _read_spatial_rows(
    path: str | Path,
    *,
    expected_fields: int | tuple[int, ...],
) -> list[_SpatialRow]:
    input_path = Path(path)
    rows: list[_SpatialRow] = []

    try:
        with input_path.open("r", encoding="utf-8") as stream:
            for line_number, raw_line in enumerate(stream, start=1):
                content = raw_line.partition("#")[0].strip()
                if not content:
                    continue

                fields = content.split()
                allowed_fields = (
                    (expected_fields,) if isinstance(expected_fields, int)
                    else expected_fields
                )
                if len(fields) not in allowed_fields:
                    expected_description = (
                        str(expected_fields) if isinstance(expected_fields, int)
                        else " or ".join(str(value) for value in expected_fields)
                    )
                    raise _row_error(
                        input_path,
                        line_number,
                        f"expected {expected_description} fields, found {len(fields)}",
                    )

                _atom_name = fields[0]
                coordinate_text = fields[1:4]
                coordinates: list[float] = []
                for axis, text in zip("xyz", coordinate_text, strict=True):
                    try:
                        value = float(text)
                    except ValueError:
                        raise _row_error(
                            input_path,
                            line_number,
                            f"invalid {axis} coordinate {text!r}",
                        ) from None
                    if not math.isfinite(value):
                        raise _row_error(
                            input_path,
                            line_number,
                            f"{axis} coordinate must be finite",
                        )
                    coordinates.append(value)

                rows.append(
                    _SpatialRow(
                        position_angstrom=(
                            coordinates[0],
                            coordinates[1],
                            coordinates[2],
                        ),
                        payload=tuple(fields[4:]),
                        line_number=line_number,
                    )
                )
    except InputFileError:
        raise
    except OSError as error:
        detail = error.strerror or str(error)
        raise InputFileError(f"{input_path}: {detail}") from None

    return rows


def read_spatial_records(path: str | Path) -> list[SpatialRecord]:
    """Read shared five-column spatial records from ``path``."""

    input_path = Path(path)
    records: list[SpatialRecord] = []
    for row in _read_spatial_rows(input_path, expected_fields=5):
        records.append(
            SpatialRecord(
                position_angstrom=row.position_angstrom,
                source_spec=row.payload[0],
                line_number=row.line_number,
            )
        )
    return records


def read_point_charges(path: str | Path) -> list[PointCharge]:
    """Read point charges, validating the fifth field as an integer."""

    input_path = Path(path)
    charges: list[PointCharge] = []
    for record in read_spatial_records(input_path):
        if _INTEGER_PATTERN.fullmatch(record.source_spec) is None:
            raise _row_error(
                input_path,
                record.line_number,
                f"charge must be an integer, got {record.source_spec!r}",
            )

        charges.append(
            PointCharge(
                position_m=tuple(
                    coordinate * ANGSTROM_IN_METRES
                    for coordinate in record.position_angstrom
                ),
                charge_e=int(record.source_spec, 10),
            )
        )

    return charges


def _parse_finite_float(
    path: Path,
    line_number: int,
    text: str,
    field_name: str,
) -> float:
    try:
        value = float(text)
    except ValueError:
        raise _row_error(
            path,
            line_number,
            f"{field_name} must be a number, got {text!r}",
        ) from None
    if not math.isfinite(value):
        raise _row_error(path, line_number, f"{field_name} must be finite")
    return value


def read_hydrogenic_orbitals(path: str | Path) -> list[HydrogenicOrbital]:
    """Read exact eight-column hydrogenic-orbital records."""

    input_path = Path(path)
    orbitals: list[HydrogenicOrbital] = []
    for row in _read_spatial_rows(input_path, expected_fields=8):
        orbital_type, radius_text, zeff_text, occupation_text = row.payload
        try:
            canonical_type, axis = parse_orbital_descriptor(orbital_type)
        except ValueError as error:
            raise _row_error(
                input_path,
                row.line_number,
                str(error),
            ) from None
        if canonical_type == "2p" and axis is None:
            raise _row_error(
                input_path,
                row.line_number,
                "orbital type '2p' requires an axis; use 2px, 2py, 2pz, "
                "or 2p(ux,uy,uz)",
            )
        if radius_text != "-" and zeff_text != "-":
            raise _row_error(
                input_path,
                row.line_number,
                "radius and Zeff are mutually exclusive; use '-' for one",
            )

        if radius_text == "-" and zeff_text == "-":
            radius_m = BOHR_RADIUS_M
        elif radius_text != "-":
            radius_angstrom = _parse_finite_float(
                input_path,
                row.line_number,
                radius_text,
                "orbital radius",
            )
            if radius_angstrom <= 0.0:
                raise _row_error(
                    input_path,
                    row.line_number,
                    "orbital radius must be positive",
                )
            radius_m = radius_angstrom * ANGSTROM_IN_METRES
        else:
            zeff = _parse_finite_float(
                input_path,
                row.line_number,
                zeff_text,
                "Zeff",
            )
            try:
                radius_m = radius_from_zeff(zeff)
            except ValueError as error:
                raise _row_error(input_path, row.line_number, str(error)) from None

        occupation = _parse_finite_float(
            input_path,
            row.line_number,
            occupation_text,
            "occupation",
        )
        if occupation < 0.0:
            raise _row_error(
                input_path,
                row.line_number,
                "occupation must be nonnegative",
            )

        orbitals.append(
            HydrogenicOrbital(
                center_m=tuple(
                    coordinate * ANGSTROM_IN_METRES
                    for coordinate in row.position_angstrom
                ),
                orbital_type=canonical_type,
                radius_m=radius_m,
                occupation=occupation,
                axis=axis,
            )
        )

    return orbitals


def read_orbital_input(
    path: str | Path,
) -> tuple[list[HydrogenicOrbital], list[OrbitalCombination]]:
    """Read independent 8-column orbitals and/or grouped 10-column SALCs.

    The grouped form is ``name x y z orbital radius zeff coefficient group
    occupation``. Rows sharing ``group`` form one coherent spatial orbital;
    their signed amplitudes interfere. Occupation is specified per group.
    """
    input_path = Path(path)
    rows = _read_spatial_rows(input_path, expected_fields=(8, 10))
    independent: list[HydrogenicOrbital] = []
    groups: dict[str, list[OrbitalCombinationTerm]] = {}
    occupations: dict[str, float] = {}

    for row in rows:
        if len(row.payload) == 4:
            orbital_type, radius_text, zeff_text, occupation_text = row.payload
            coefficient = None
            group_name = None
        else:
            (orbital_type, radius_text, zeff_text, coefficient_text,
             group_name, occupation_text) = row.payload
            coefficient = _parse_finite_float(
                input_path, row.line_number, coefficient_text, "coefficient"
            )
            if not group_name:
                raise _row_error(input_path, row.line_number, "SALC group name cannot be empty")

        try:
            canonical_type, axis = parse_orbital_descriptor(orbital_type)
        except ValueError as error:
            raise _row_error(input_path, row.line_number, str(error)) from None
        if canonical_type == "2p" and axis is None:
            raise _row_error(
                input_path, row.line_number,
                "orbital type '2p' requires an axis; use 2px, 2py, 2pz, or 2p(ux,uy,uz)",
            )
        if radius_text != "-" and zeff_text != "-":
            raise _row_error(
                input_path, row.line_number,
                "radius and Zeff are mutually exclusive; use '-' for one",
            )
        if radius_text == "-" and zeff_text == "-":
            radius_m = BOHR_RADIUS_M
        elif radius_text != "-":
            radius_angstrom = _parse_finite_float(
                input_path, row.line_number, radius_text, "orbital radius"
            )
            if radius_angstrom <= 0.0:
                raise _row_error(input_path, row.line_number, "orbital radius must be positive")
            radius_m = radius_angstrom * ANGSTROM_IN_METRES
        else:
            zeff = _parse_finite_float(input_path, row.line_number, zeff_text, "Zeff")
            if zeff <= 0.0:
                raise _row_error(input_path, row.line_number, "Zeff must be positive")
            try:
                radius_m = radius_from_zeff(zeff)
            except ValueError as error:
                raise _row_error(input_path, row.line_number, str(error)) from None

        occupation = _parse_finite_float(
            input_path, row.line_number, occupation_text, "occupation"
        )
        if occupation < 0.0:
            raise _row_error(input_path, row.line_number, "occupation must be nonnegative")
        orbital = HydrogenicOrbital(
            center_m=tuple(value * ANGSTROM_IN_METRES for value in row.position_angstrom),
            orbital_type=canonical_type,
            radius_m=radius_m,
            occupation=occupation if coefficient is None else 1.0,
            axis=axis,
        )
        if coefficient is None:
            independent.append(orbital)
            continue

        assert group_name is not None
        if group_name in occupations and occupations[group_name] != occupation:
            raise _row_error(
                input_path, row.line_number,
                f"occupation for SALC group {group_name!r} does not match its other rows",
            )
        occupations[group_name] = occupation
        groups.setdefault(group_name, []).append(
            OrbitalCombinationTerm(orbital=orbital, coefficient=coefficient)
        )

    combinations = [
        OrbitalCombination(
            terms=tuple(terms), occupation=occupations[name], name=name
        )
        for name, terms in groups.items()
    ]
    return independent, combinations
