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

from magnopy import converter22

MAX_MODULUS = 1e8
ARRAY_3X3 = harrays(
    np.float64,
    (3, 3),
    elements=st.floats(min_value=-MAX_MODULUS, max_value=MAX_MODULUS),
)
ARRAY_3 = harrays(
    np.float64,
    (3,),
    elements=st.floats(min_value=-MAX_MODULUS, max_value=MAX_MODULUS),
)


@given(st.floats(min_value=-MAX_MODULUS, max_value=MAX_MODULUS))
def test_from_iso(iso):
    parameter = converter22.from_iso(iso=iso)
    assert np.allclose(parameter, iso * np.eye(3))


@given(ARRAY_3)
def test_from_dmi(dmi):
    parameter = converter22.from_dmi(dmi=dmi)
    assert np.allclose(
        parameter,
        [
            [0, dmi[2], -dmi[1]],
            [-dmi[2], 0, dmi[0]],
            [dmi[1], -dmi[0], 0],
        ],
    )


@given(ARRAY_3X3)
def test_to_iso(parameter):
    iso = converter22.to_iso(parameter=parameter)
    assert np.allclose(iso, np.trace(parameter) / 3.0)

    iso_parameter = converter22.to_iso(parameter=parameter, matrix_form=True)
    assert np.allclose(iso_parameter, np.trace(parameter) / 3.0 * np.eye(3))


@given(ARRAY_3X3)
def test_to_symm_anisotropy(parameter):
    aniso = converter22.to_symm_anisotropy(parameter=parameter)
    assert np.allclose(
        aniso, (parameter + parameter.T) / 2.0 - np.trace(parameter) / 3.0 * np.eye(3)
    )


@given(ARRAY_3X3)
def test_to_dmi(parameter):
    asymm = (parameter - parameter.T) / 2.0

    dmi = converter22.to_dmi(parameter=parameter)
    assert np.allclose(dmi, [asymm[1][2], asymm[2][0], asymm[0][1]])

    dmi_parameter = converter22.to_dmi(parameter=parameter, matrix_form=True)
    assert np.allclose(dmi_parameter, asymm)
