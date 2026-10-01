"""Shared text input handling for spatial source records."""

from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path
import re

from .models import PointCharge


ANGSTROM_IN_METRES = 1.0e-10
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


def _row_error(path: Path, line_number: int, message: str) -> InputFileError:
    return InputFileError(f"{path}:{line_number}: {message}")


def read_spatial_records(path: str | Path) -> list[SpatialRecord]:
    """Read shared five-column spatial records from ``path``."""

    input_path = Path(path)
    records: list[SpatialRecord] = []

    try:
        with input_path.open("r", encoding="utf-8") as stream:
            for line_number, raw_line in enumerate(stream, start=1):
                content = raw_line.partition("#")[0].strip()
                if not content:
                    continue

                fields = content.split()
                if len(fields) != 5:
                    raise _row_error(
                        input_path,
                        line_number,
                        f"expected 5 fields, found {len(fields)}",
                    )

                _atom_name, *coordinate_text, source_spec = fields
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

                records.append(
                    SpatialRecord(
                        position_angstrom=(
                            coordinates[0],
                            coordinates[1],
                            coordinates[2],
                        ),
                        source_spec=source_spec,
                        line_number=line_number,
                    )
                )
    except InputFileError:
        raise
    except OSError as error:
        detail = error.strerror or str(error)
        raise InputFileError(f"{input_path}: {detail}") from None

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
