"""Tests for spherical harmonics."""

import pytest
from sympy import I, Symbol, cos, exp, pi, simplify, sin, sqrt

from polyharmonics import spherical_harmonic

th = Symbol("θ")
phi = Symbol("φ")


@pytest.mark.parametrize(
    ("n", "m", "evaluate", "expected"),
    [
        (0, 0, None, 1 / sqrt(4 * pi)),
        (2, 2, None, sqrt(15 / (32 * pi)) * sin(th) ** 2 * exp(2 * I * phi)),
        (
            3,
            2,
            None,
            sqrt(105 / (32 * pi)) * cos(th) * sin(th) ** 2 * exp(2 * I * phi),
        ),
        (0, 0, (0, 1), 1 / sqrt(4 * pi)),
        (
            2,
            2,
            (pi / 2, pi),
            sqrt(15 / (32 * pi)),
        ),
        (3, 2, (-pi / 4, -pi), sqrt(105 / (32 * pi)) / (2 * sqrt(2))),
    ],
)
def test_spherical_harmonic(n, m, evaluate, expected):
    """Test the spherical harmonic with parametrization."""
    if evaluate is None:
        assert simplify(spherical_harmonic(n, m, eval=evaluate) - expected) == 0
    else:
        res = complex(spherical_harmonic(n, m, eval=evaluate))
        expected = complex(expected)
        assert simplify(res.real - expected.real) == pytest.approx(0)
        assert simplify(res.imag - expected.imag) == pytest.approx(0)
