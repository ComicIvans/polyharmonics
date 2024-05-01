"""Tests for legendre polynomials."""

import pytest
from sympy import Symbol, simplify

from polyharmonics import legendre

x = Symbol("x")


@pytest.mark.parametrize(
    ("n", "expected"),
    [
        (0, 1),
        (1, x),
        (2, 3 / 2 * x**2 - 1 / 2),
        ([3, 2], [5 / 2 * x**3 - 3 / 2 * x, 3 / 2 * x**2 - 1 / 2]),
    ],
)
def test_legendre(n, expected):
    """Test the Legendre function with parametrization."""
    if isinstance(n, int):
        assert simplify(legendre(n) - expected) == 0
    else:
        pols = legendre(n)
        assert all(simplify(pol - exp) == 0 for pol, exp in zip(pols, expected))
