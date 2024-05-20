"""Tests for legendre polynomials."""

import pytest
from sympy import Symbol, simplify

from polyharmonics import legendre
from polyharmonics.legendre_polynomials import (
    legendre_def,
    legendre_exp,
    legendre_rec,
    legendre_store,
)

x = Symbol("x")


@pytest.mark.parametrize(
    ("n", "expected"),
    [
        (0, 1),
        (1, x),
        (2, 3 / 2 * x**2 - 1 / 2),
        (3, 5 / 2 * x**3 - 3 / 2 * x),
    ],
)
def test_legendre(n, expected):
    """Test the Legendre function with parametrization."""
    assert simplify(legendre(n) - expected) == 0
    for i in range(5, 10):
        legendre_store.reset()
        assert simplify(legendre_rec(i, store=True) - legendre_def(i, store=True)) == 0
        assert simplify(legendre_rec(i, store=False) - legendre_def(i, store=False)) == 0
        assert simplify(legendre_rec(i, store=False) - legendre_def(i, store=True)) == 0
        assert simplify(legendre_rec(i, store=True) - legendre_exp(i)) == 0
