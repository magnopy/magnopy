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

from magnopy import Convention, SpinHamiltonian, make_supercell

MAX_MODULUS = 1e8
ARRAY = harrays(
    np.float64,
    (3, 3, 3),
    elements=st.floats(min_value=-MAX_MODULUS, max_value=MAX_MODULUS),
)

CONVENTION = Convention(
    spin_normalized=False,
    multiple_counting=True,
    c1=1,
    c21=1,
    c22=1,
    c31=1,
    c32=1,
    c33=1,
    c41=1,
    c421=1,
    c422=1,
    c43=1,
    c44=1,
)


@given(st.integers(), ARRAY)
def test_add_31(alpha, parameter):
    atoms = {"names": ["Cr" for _ in range(9)], "spins": [1 for _ in range(9)]}

    spinham = SpinHamiltonian(cell=np.eye(3), atoms=atoms, convention=CONVENTION)

    if 0 <= alpha < len(spinham.atoms.names):
        spinham.add_31(alpha, parameter)
    else:
        with pytest.raises(ValueError):
            spinham.add_31(alpha, parameter)


@given(
    st.integers(min_value=0, max_value=8),
    st.integers(min_value=0, max_value=8),
    st.integers(min_value=0, max_value=8),
    st.integers(min_value=0, max_value=8),
    ARRAY,
)
def test_add_31_sorting(alpha1, alpha2, alpha3, alpha4, parameter):
    atoms = {"names": ["Cr" for _ in range(9)], "spins": [1 for _ in range(9)]}

    spinham = SpinHamiltonian(cell=np.eye(3), atoms=atoms, convention=CONVENTION)

    spinham.add_31(alpha1, parameter)

    if alpha2 == alpha1:
        with pytest.raises(ValueError):
            spinham.add_31(alpha2, parameter)
    else:
        spinham.add_31(alpha2, parameter)

    spinham.add_31(alpha3, parameter, replace=True)
    spinham.add_31(alpha4, parameter, replace=True)

    for i in range(len(spinham._31) - 1):
        assert spinham._31[i][:-1] <= spinham._31[i + 1][:-1]


@given(st.integers())
def test_remove_31(r_alpha):
    atoms = {"names": ["Cr" for _ in range(9)], "spins": [1 for _ in range(9)]}

    spinham = SpinHamiltonian(cell=np.eye(3), atoms=atoms, convention=CONVENTION)

    for alpha in range(len(spinham.atoms.names)):
        spinham.add_31(alpha, np.eye(3))

    bond = [r_alpha]
    if 0 <= r_alpha < len(spinham.atoms.names):
        original_bonds = [tmp[:-1] for tmp in spinham._31]
        original_length = len(spinham._31)

        spinham.remove_31(*bond)

        if bond in original_bonds:
            updated_bonds = [tmp[:-1] for tmp in spinham._31]

            assert len(spinham._31) == original_length - 1
            assert bond not in updated_bonds
        else:
            assert len(spinham._31) == original_length
    else:
        with pytest.raises(ValueError):
            spinham.remove_31(*bond)


@given(ARRAY, st.floats(min_value=0.1, max_value=1e4))
def test_mul(parameter, number):
    atoms = {"names": ["Cr" for _ in range(9)], "spins": [1 for _ in range(9)]}

    spinham = SpinHamiltonian(cell=np.eye(3), atoms=atoms, convention=CONVENTION)

    spinham.add_31(alpha=0, parameter=parameter)
    spinham.add_31(alpha=4, parameter=parameter * 1.32)
    spinham.add_31(alpha=7, parameter=parameter / 3)

    m_spinham = spinham * number

    assert spinham.M == m_spinham.M

    assert len(spinham.p31) == len(m_spinham.p31)

    params = list(spinham.p31)
    m_params = list(m_spinham.p31)
    for i in range(len(params)):
        assert params[i][0] == m_params[i][0]

        assert np.allclose(number * params[i][1], m_params[i][1])


@given(ARRAY, st.floats(min_value=0.1, max_value=1e4))
def test_rmul(parameter, number):
    atoms = {"names": ["Cr" for _ in range(9)], "spins": [1 for _ in range(9)]}

    spinham = SpinHamiltonian(cell=np.eye(3), atoms=atoms, convention=CONVENTION)

    spinham.add_31(alpha=0, parameter=parameter)
    spinham.add_31(alpha=4, parameter=parameter * 1.32)
    spinham.add_31(alpha=7, parameter=parameter / 3)

    m_spinham = number * spinham

    assert spinham.M == m_spinham.M

    assert len(spinham.p31) == len(m_spinham.p31)

    params = list(spinham.p31)
    m_params = list(m_spinham.p31)
    for i in range(len(params)):
        assert params[i][0] == m_params[i][0]

        assert np.allclose(number * params[i][1], m_params[i][1])


@given(ARRAY, ARRAY)
def test_add(parameter1, parameter2):
    atoms = dict(
        names=["Cr" for _ in range(9)],
        spins=[1 for _ in range(9)],
        positions=[[0.1 * i, 0, 0] for i in range(9)],
        g_factors=[2 for _ in range(9)],
    )

    spinham1 = SpinHamiltonian(cell=np.eye(3), atoms=atoms, convention=CONVENTION)
    spinham2 = SpinHamiltonian(cell=np.eye(3), atoms=atoms, convention=CONVENTION)

    spinham1.add_31(alpha=0, parameter=parameter1)
    spinham1.add_31(alpha=4, parameter=parameter1 * 1.32)
    spinham1.add_31(alpha=7, parameter=parameter1 / 3)

    spinham2.add_31(alpha=0, parameter=parameter2)
    spinham2.add_31(alpha=4, parameter=parameter2 * 1.32)
    spinham2.add_31(alpha=8, parameter=parameter2 / 3)

    m_spinham = spinham1 + spinham2

    assert m_spinham.M == 4

    assert len(m_spinham.p31) == 4

    for i in range(2):
        assert np.allclose(m_spinham._31[i][1], spinham1._31[i][1] + spinham2._31[i][1])

    assert np.allclose(m_spinham._31[2][1], spinham1._31[2][1])
    assert np.allclose(m_spinham._31[3][1], spinham2._31[2][1])


@given(
    ARRAY,
    st.integers(min_value=1, max_value=5),
    st.integers(min_value=1, max_value=5),
    st.integers(min_value=1, max_value=5),
)
def test_make_supercell(parameter1, i, j, k):
    atoms = dict(
        names=["Cr1", "Cr2"],
        spins=[1, 2],
        positions=[[0, 0, 0], [0.5, 0.5, 0.5]],
        g_factors=[2, 2],
    )

    spinham = SpinHamiltonian(cell=np.eye(3), atoms=atoms, convention=CONVENTION)

    spinham.add_31(alpha=0, parameter=parameter1)
    spinham.add_31(alpha=1, parameter=parameter1 * 1.42)

    new_spinham = make_supercell(spinham=spinham, supercell=(i, j, k))

    assert len(new_spinham.p31) == i * j * k * 2
