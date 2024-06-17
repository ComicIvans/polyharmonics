"""Tests for spherical harmonics."""

from math import pi as pi_num

import pytest
from sympy import I, Symbol, cos, exp, pi, simplify, sin, sqrt

from polyharmonics import spherical_harmonic

th = Symbol("θ")
phi = Symbol("φ")


@pytest.mark.parametrize(
    ("n", "m", "th", "phi", "expected"),
    [
        (0, 0, None, None, 1 / sqrt(4 * pi)),
        (2, 2, None, None, sqrt(15 / (32 * pi)) * sin(th) ** 2 * exp(2 * I * phi)),
        (
            3,
            2,
            None,
            None,
            sqrt(105 / (32 * pi)) * cos(th) * sin(th) ** 2 * exp(2 * I * phi),
        ),
        (0, 0, 0, 1, 1 / sqrt(4 * pi_num)),
        (
            2,
            2,
            pi_num / 2,
            pi_num,
            sqrt(15 / (32 * pi_num)),
        ),
        (3, 2, -pi_num / 4, -pi_num, sqrt(105 / (32 * pi_num)) / (2 * sqrt(2))),
    ],
)
def test_spherical_harmonic(n, m, th, phi, expected):
    """Test the spherical harmonic with parametrization."""
    if th is None and phi is None:
        assert simplify(spherical_harmonic(n, m) - expected) == 0
    else:
        res = complex(spherical_harmonic(n, m, th=th, phi=phi))
        expected = complex(expected)
        assert simplify(res.real - expected.real) == pytest.approx(0)
        assert simplify(res.imag - expected.imag) == pytest.approx(0)
