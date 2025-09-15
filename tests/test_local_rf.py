# ================================== LICENSE ===================================
# Magnopy - Python package for magnons.
# Copyright (C) 2023-2025 Magnopy Team
#
# e-mail: anry@uv.es, web: magnopy.org
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
# ================================ END LICENSE =================================


import numpy as np
import pytest
from hypothesis import given
from hypothesis import strategies as st
from hypothesis.extra.numpy import arrays as harrays

from magnopy import span_local_rf


def test_span_local_rf_along_z():
    direction = [0, 0, 1]

    x_a, y_a, z_a = span_local_rf(direction)

    assert np.allclose(x_a, [1, 0, 0])
    assert np.allclose(y_a, [0, 1, 0])
    assert np.allclose(z_a, [0, 0, 1])


def test_span_local_rf_opposite_z():
    direction_vector = [0, 0, -1]

    x_a, y_a, z_a = span_local_rf(direction_vector)

    assert np.allclose(x_a, [0, -1, 0])
    assert np.allclose(y_a, [-1, 0, 0])
    assert np.allclose(z_a, [0, 0, -1])


@given(
    harrays(
        np.float64,
        (3,),
        elements=st.floats(min_value=-1e8, max_value=1e8, allow_subnormal=False),
    )
)
def test_span_local_rf(direction_vector):
    if np.allclose(direction_vector, np.zeros(3)):
        with pytest.raises(ValueError):
            x_a, y_a, z_a = span_local_rf(direction_vector)

    else:
        x_a, y_a, z_a = span_local_rf(direction_vector)

        assert np.linalg.det([x_a, y_a, z_a]) > 0.0

        assert np.allclose(
            [np.linalg.norm(x_a), np.linalg.norm(y_a), np.linalg.norm(z_a)], np.ones(3)
        )

        assert np.allclose([x_a @ y_a, x_a @ z_a, y_a @ z_a], np.zeros(3))
