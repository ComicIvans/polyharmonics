"""Tests for associated legendre functions."""

from typing import List

import pytest
from sympy import Expr, Rational, Symbol, acos, cos, simplify, sin

from polyharmonics import associated_legendre
from polyharmonics.associated_legendre_functions import (
    ass_legendre_store,
    associated_legendre_def,
    associated_legendre_rec,
    associated_legendre_rec_alt,
)

x = Symbol("x")
th = Symbol("θ")


@pytest.mark.parametrize(
    ("n", "m", "polar", "expected"),
    [
        (0, 0, False, 1),
        (2, 2, True, 3 * sin(th) ** 2),
        (3, 2, True, 15 * cos(th) * sin(th) ** 2),
        (10, 100, False, 0),
        (1, -1, False, -2 * (1 - x**2) ** Rational(1, 2)),
        (4, 4, False, 105 * x**4 - 210 * x**2 + 105),
        (
            3,
            3,
            False,
            -15 * x**2 * (1 - x**2) ** Rational(1, 2) + 15 * (1 - x**2) ** Rational(1, 2),
        ),
    ],
)
def test_associated_legendre(n, m, polar, expected):
    """Test the associated Legendre function with parametrization."""
    assert simplify(associated_legendre(n, m, polar=polar) - expected) == 0
    for i in range(5, 10):
        for j in range(-2, 3):
            ass_legendre_store.reset()
            fun_pairs: List[Expr] = [
                associated_legendre_def(i, j, store=True),
                associated_legendre_rec(i, j, store=False),
                associated_legendre_def(i, j, store=False),
                associated_legendre_rec(i, j, store=True),
                associated_legendre_rec_alt(i, j),
                associated_legendre_def(i, j),
            ]
            x_vals = [-0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75]
            for x_val in x_vals:
                for i in range(0, len(fun_pairs), 2):
                    assert (
                        fun_pairs[i].subs(x, x_val) - fun_pairs[i + 1].subs(x, x_val)
                    ) == 0
            polar_fun_pairs: List[Expr] = [
                associated_legendre_def(i, j, store=True, polar=True),
                associated_legendre_rec_alt(i, j, polar=True),
                associated_legendre_rec(i, j, store=True, polar=True),
                associated_legendre_rec_alt(i, j, polar=True),
            ]
            for th_val in [acos(x_val) for x_val in x_vals]:
                for i in range(0, len(polar_fun_pairs), 2):
                    assert (
                        fun_pairs[i].subs(th, th_val) - fun_pairs[i + 1].subs(th, th_val)
                    ) == 0
