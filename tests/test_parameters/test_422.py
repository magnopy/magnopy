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
from hypothesis import given
from hypothesis import strategies as st
from hypothesis.extra.numpy import arrays as harrays

from magnopy import converter422

MAX_MODULUS = 1e4
ARRAY_3X3X3X3 = harrays(
    np.float64,
    (3, 3, 3, 3),
    elements=st.floats(min_value=-MAX_MODULUS, max_value=MAX_MODULUS),
)

INDICES = [
    [0, 0, 0, 0],
    [0, 1, 0, 1],
    [0, 2, 0, 2],
    [1, 0, 1, 0],
    [1, 1, 1, 1],
    [1, 2, 1, 2],
    [2, 0, 2, 0],
    [2, 1, 2, 1],
    [2, 2, 2, 2],
]


@given(st.floats(min_value=-MAX_MODULUS, max_value=MAX_MODULUS))
def test_from_biquadratic(B):
    parameter = converter422.from_biquadratic(B=B)

    for i in range(3):
        for j in range(3):
            for u in range(3):
                for v in range(3):
                    if [i, j, u, v] not in INDICES:
                        assert parameter[i, j, u, v] == 0
                    else:
                        assert parameter[i, j, u, v] == B


@given(ARRAY_3X3X3X3)
def test_to_biquadratic(parameter):
    B = converter422.to_biquadratic(parameter=parameter)

    real_B = np.sum([parameter[i, j, u, v] for (i, j, u, v) in INDICES]) / 9

    assert abs(B - real_B) < 1e-8
