"""Tests for associated legendre functions."""

import pytest
from sympy import Rational, Symbol, simplify

from polyharmonics import associated_legendre

x = Symbol("x")


@pytest.mark.parametrize(
    ("n", "m", "expected"),
    [
        (0, 0, 1),
        (10, 100, 0),
        (1, -1, -2 * (1 - x**2) ** Rational(1, 2)),
        (
            [3, 3, 4],
            [3, 2, 4],
            [
                -15 * x**2 * (1 - x**2) ** Rational(1, 2)
                + 15 * (1 - x**2) ** Rational(1, 2),
                -15 * x**3 + 15 * x,
                105 * x**4 - 210 * x**2 + 105,
            ],
        ),
    ],
)
def test_associated_legendre(n, m, expected):
    """Test the associated Legendre function with parametrization."""
    if isinstance(n, int) and isinstance(m, int):
        assert simplify(associated_legendre(n, m) - expected) == 0
    else:
        pols = associated_legendre(n, m)
        assert all(simplify(pol - exp) == 0 for pol, exp in zip(pols, expected))
