"""Tests for legendre polynomials."""

from random import random

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
    for i in range(5, 10, 2):
        val = random() * 2 - 1
        legendre_store.reset()
        assert simplify(legendre_rec(i, store=True) - legendre_def(i, store=True)) == 0
        assert legendre_rec(i, eval=val, store=True) - legendre_def(
            i, eval=val, store=False
        ) == pytest.approx(0)
        assert simplify(legendre_rec(i, store=False) - legendre_def(i, store=False)) == 0
        assert legendre_def(i, eval=val, store=False) - legendre_def(
            i, eval=val, store=False
        ) == pytest.approx(0)
        assert simplify(legendre_rec(i, store=False) - legendre_def(i, store=True)) == 0
        assert legendre_exp(i, eval=val) - legendre_rec(i, eval=val) == pytest.approx(0)
        assert simplify(legendre_rec(i, store=True) - legendre_exp(i)) == 0
        assert legendre_rec(i, eval=val) - legendre_exp(i).subs(x, val) == pytest.approx(0)
