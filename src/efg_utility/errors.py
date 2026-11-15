"""Exceptions shared by the numerical calculation backends."""


class CalculationError(ValueError):
    """Raised when a requested field calculation is invalid or non-finite."""
