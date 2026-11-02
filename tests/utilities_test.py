"""Tests for 3D geometry helpers in pdb2pqr.utilities."""

import math

import numpy as np
import pytest

from pdb2pqr import utilities as util


@pytest.mark.parametrize("wrap", [list, tuple, np.array])
def test_vector_helpers_accept_sequence_types(wrap):
    """3D helpers accept lists, tuples, and ndarrays and match NumPy."""
    a, b = wrap([1.0, -2.0, 3.5]), wrap([-0.5, 4.0, 1.0])
    assert util.distance(a, b) == pytest.approx(
        np.linalg.norm(np.array(a) - np.array(b))
    )
    assert util.subtract(a, b) == [1.5, -6.0, 2.5]
    assert util.cross(a, b) == [-16.0, -2.75, 3.0]
    assert util.dot(a, b) == pytest.approx(-5.0)
    assert np.allclose(util.normalize(a), np.array(a) / np.linalg.norm(a))


def test_dihedral_signed_and_degenerate_inputs():
    """Dihedral preserves chirality, planar limits, and degenerate NaNs."""
    p1, p2, p3 = [1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0]
    right_deg = 0.5 * math.pi * util.RADIANS_TO_DEGREES
    assert util.dihedral(p1, p2, p3, [0.0, 1.0, -1.0]) == pytest.approx(
        right_deg
    )
    assert util.dihedral(p1, p2, p3, [0.0, 1.0, 1.0]) == pytest.approx(
        -right_deg
    )
    assert util.dihedral(p1, p2, p3, [-1.0, 1.0, 0.0]) == 180.0
    assert util.dihedral(p1, p2, p3, [1.0, 1.0, 0.0]) == 0.0
    assert all(math.isnan(v) for v in util.normalize([0.0, 0.0, 0.0]))
    assert math.isnan(
        util.dihedral(
            [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [2.0, 1.0, 0.0]
        )
    )
