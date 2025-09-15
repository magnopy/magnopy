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
from magnopy._spinham._c43 import _get_primary_p43

MAX_MODULUS = 1e8
ARRAY = harrays(
    np.float64,
    (3, 3, 3, 3),
    elements=st.floats(min_value=-MAX_MODULUS, max_value=MAX_MODULUS),
)
RANDOM_UC = harrays(int, (4, 3), elements=st.integers(min_value=-1000, max_value=1000))

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


@given(
    st.integers(),
    st.integers(),
    st.integers(),
    st.tuples(st.integers(), st.integers(), st.integers()),
    st.tuples(st.integers(), st.integers(), st.integers()),
    ARRAY,
)
def test_add_43(alpha, beta, gamma, nu, _lambda, parameter):
    atoms = {"names": ["Cr" for _ in range(9)], "spins": [1 for _ in range(9)]}

    spinham = SpinHamiltonian(cell=np.eye(3), atoms=atoms, convention=CONVENTION)

    if (
        0 <= alpha < len(spinham.atoms.names)
        and 0 <= beta < len(spinham.atoms.names)
        and 0 <= gamma < len(spinham.atoms.names)
    ):
        spinham.add_43(alpha, beta, gamma, nu, _lambda, parameter)
    else:
        with pytest.raises(ValueError):
            spinham.add_43(alpha, beta, gamma, nu, _lambda, parameter)


@given(
    st.integers(min_value=0, max_value=8),
    st.integers(min_value=0, max_value=8),
    st.integers(min_value=0, max_value=8),
    st.tuples(st.integers(), st.integers(), st.integers()),
    st.tuples(st.integers(), st.integers(), st.integers()),
    st.integers(min_value=0, max_value=8),
    st.integers(min_value=0, max_value=8),
    st.integers(min_value=0, max_value=8),
    st.tuples(st.integers(), st.integers(), st.integers()),
    st.tuples(st.integers(), st.integers(), st.integers()),
    st.integers(min_value=0, max_value=8),
    st.integers(min_value=0, max_value=8),
    st.integers(min_value=0, max_value=8),
    st.tuples(st.integers(), st.integers(), st.integers()),
    st.tuples(st.integers(), st.integers(), st.integers()),
    st.integers(min_value=0, max_value=8),
    st.integers(min_value=0, max_value=8),
    st.integers(min_value=0, max_value=8),
    st.tuples(st.integers(), st.integers(), st.integers()),
    st.tuples(st.integers(), st.integers(), st.integers()),
    ARRAY,
)
def test_add_43_sorting(
    alpha1,
    beta1,
    gamma1,
    nu1,
    _lambda1,
    alpha2,
    beta2,
    gamma2,
    nu2,
    _lambda2,
    alpha3,
    beta3,
    gamma3,
    nu3,
    _lambda3,
    alpha4,
    beta4,
    gamma4,
    nu4,
    _lambda4,
    parameter,
):
    atoms = {"names": ["Cr" for _ in range(9)], "spins": [1 for _ in range(9)]}

    spinham = SpinHamiltonian(cell=np.eye(3), atoms=atoms, convention=CONVENTION)

    spinham.add_43(alpha1, beta1, gamma1, nu1, _lambda1, parameter)

    if [alpha1, beta1, gamma1, nu1, _lambda1] == [alpha2, beta2, gamma2, nu2, _lambda2]:
        with pytest.raises(ValueError):
            spinham.add_43(alpha2, beta2, gamma2, nu2, _lambda2, parameter)
    else:
        spinham.add_43(alpha2, beta2, gamma2, nu2, _lambda2, parameter)

    spinham.add_43(alpha3, beta3, gamma3, nu3, _lambda3, parameter, replace=True)
    spinham.add_43(alpha4, beta4, gamma4, nu4, _lambda4, parameter, replace=True)

    for i in range(len(spinham._43) - 1):
        assert spinham._43[i][:-1] <= spinham._43[i + 1][:-1]


@given(
    st.integers(min_value=0, max_value=2),
    st.integers(min_value=0, max_value=2),
    st.integers(min_value=0, max_value=2),
    st.tuples(st.integers(), st.integers(), st.integers()),
    st.tuples(st.integers(), st.integers(), st.integers()),
    RANDOM_UC,
    RANDOM_UC,
)
def test_remove_43(r_alpha, r_beta, r_gamma, r_nu, r_lambda, nus, lambdas):
    atoms = {"names": ["Cr" for _ in range(4)], "spins": [1 for _ in range(4)]}

    spinham = SpinHamiltonian(cell=np.eye(3), atoms=atoms, convention=CONVENTION)

    for alpha in range(len(spinham.atoms.names)):
        for beta in range(alpha, len(spinham.atoms.names)):
            for gamma in range(beta, len(spinham.atoms.names)):
                for nu in nus:
                    for _lambda in lambdas:
                        nu = (int(nu[0]), int(nu[1]), int(nu[2]))
                        _lambda = (int(_lambda[0]), int(_lambda[1]), int(_lambda[2]))
                        spinham.add_43(
                            alpha,
                            beta,
                            gamma,
                            nu,
                            _lambda,
                            np.ones((3, 3, 3, 3)),
                            replace=True,
                        )

    bond = [r_alpha, r_beta, r_gamma, r_nu, r_lambda]
    if (
        0 <= r_alpha < len(spinham.atoms.names)
        and 0 <= r_beta < len(spinham.atoms.names)
        and 0 <= r_gamma < len(spinham.atoms.names)
    ):
        original_bonds = [tmp[:-1] for tmp in spinham._43]
        original_length = len(spinham._43)

        spinham.remove_43(*bond)

        primary_bond = list(_get_primary_p43(*bond))
        if primary_bond in original_bonds:
            updated_bonds = [tmp[:-1] for tmp in spinham._43]

            assert len(spinham._43) == original_length - 1
            assert bond not in updated_bonds
            assert primary_bond not in updated_bonds
        else:
            assert len(spinham._43) == original_length

    else:
        with pytest.raises(ValueError):
            spinham.remove_43(*bond)


@given(ARRAY)
def test_iterator_index_bug(parameter):
    atoms = {"names": ["Cr" for _ in range(9)], "spins": [1 for _ in range(9)]}

    spinham = SpinHamiltonian(
        cell=np.eye(3),
        atoms=atoms,
        convention=Convention(multiple_counting=True, spin_normalized=False),
    )

    spinham.add_43(
        alpha=7, beta=8, gamma=4, nu=(0, 3, -4), _lambda=(1, 3, 6), parameter=parameter
    )

    list(spinham.p43)


@given(ARRAY, st.floats(min_value=0.1, max_value=1e4))
def test_mul(parameter, number):
    atoms = {"names": ["Cr" for _ in range(9)], "spins": [1 for _ in range(9)]}

    spinham = SpinHamiltonian(
        cell=np.eye(3),
        atoms=atoms,
        convention=Convention(multiple_counting=True, spin_normalized=False),
    )

    spinham.add_43(
        alpha=0, beta=1, gamma=4, nu=(0, 3, -4), _lambda=(1, 3, 6), parameter=parameter
    )
    spinham.add_43(
        alpha=4,
        beta=2,
        gamma=3,
        nu=(1, 0, 0),
        _lambda=(0, 0, 5),
        parameter=parameter * 1.43,
    )
    spinham.add_43(
        alpha=7,
        beta=5,
        gamma=2,
        nu=(0, 0, 0),
        _lambda=(-1, 2, 3),
        parameter=parameter / 3,
    )

    m_spinham = spinham * number

    assert spinham.M == m_spinham.M

    assert len(spinham.p43) == len(m_spinham.p43)

    params = list(spinham.p43)
    m_params = list(m_spinham.p43)
    for i in range(len(params)):
        assert params[i][0] == m_params[i][0]
        assert params[i][1] == m_params[i][1]
        assert params[i][2] == m_params[i][2]
        assert params[i][3] == m_params[i][3]
        assert params[i][4] == m_params[i][4]

        assert np.allclose(number * params[i][5], m_params[i][5])


@given(ARRAY, st.floats(min_value=0.1, max_value=1e4))
def test_rmul(parameter, number):
    atoms = {"names": ["Cr" for _ in range(9)], "spins": [1 for _ in range(9)]}

    spinham = SpinHamiltonian(
        cell=np.eye(3),
        atoms=atoms,
        convention=Convention(multiple_counting=True, spin_normalized=False),
    )

    spinham.add_43(
        alpha=0, beta=1, gamma=4, nu=(0, 3, -4), _lambda=(1, 3, 6), parameter=parameter
    )
    spinham.add_43(
        alpha=4,
        beta=2,
        gamma=3,
        nu=(1, 0, 0),
        _lambda=(0, 0, 5),
        parameter=parameter * 1.43,
    )
    spinham.add_43(
        alpha=7,
        beta=5,
        gamma=2,
        nu=(0, 0, 0),
        _lambda=(-1, 2, 3),
        parameter=parameter / 3,
    )

    m_spinham = number * spinham

    assert spinham.M == m_spinham.M

    assert len(spinham.p43) == len(m_spinham.p43)

    params = list(spinham.p43)
    m_params = list(m_spinham.p43)
    for i in range(len(params)):
        assert params[i][0] == m_params[i][0]
        assert params[i][1] == m_params[i][1]
        assert params[i][2] == m_params[i][2]
        assert params[i][3] == m_params[i][3]
        assert params[i][4] == m_params[i][4]

        assert np.allclose(number * params[i][5], m_params[i][5])


@given(ARRAY, ARRAY)
def test_add(parameter1, parameter2):
    atoms = dict(
        names=["Cr" for _ in range(9)],
        spins=[1 for _ in range(9)],
        positions=[[0.1 * i, 0, 0] for i in range(9)],
        g_factors=[2 for _ in range(9)],
    )

    spinham1 = SpinHamiltonian(
        cell=np.eye(3), atoms=atoms, convention=Convention(multiple_counting=True)
    )
    spinham2 = SpinHamiltonian(
        cell=np.eye(3), atoms=atoms, convention=Convention(multiple_counting=True)
    )

    spinham1.add_43(
        alpha=0, beta=1, gamma=4, nu=(0, 3, -4), _lambda=(1, 3, 6), parameter=parameter1
    )
    spinham1.add_43(
        alpha=4,
        beta=2,
        gamma=3,
        nu=(1, 0, 0),
        _lambda=(0, 0, 5),
        parameter=parameter1 * 1.43,
    )
    spinham1.add_43(
        alpha=7,
        beta=5,
        gamma=2,
        nu=(0, 0, 0),
        _lambda=(1, 2, 3),
        parameter=parameter1 / 3,
    )

    spinham2.add_43(
        alpha=0, beta=1, gamma=4, nu=(0, 3, -4), _lambda=(1, 3, 6), parameter=parameter2
    )
    spinham2.add_43(
        alpha=4,
        beta=2,
        gamma=3,
        nu=(1, 0, 0),
        _lambda=(0, 0, 5),
        parameter=parameter2 * 1.43,
    )
    spinham2.add_43(
        alpha=8,
        beta=5,
        gamma=2,
        nu=(0, 0, 0),
        _lambda=(1, 2, 3),
        parameter=parameter2 / 3,
    )

    m_spinham = spinham1 + spinham2

    assert m_spinham.M == 8

    assert len(m_spinham.p43) == 24

    for i in range(2):
        assert np.allclose(
            m_spinham._43[i][-1], spinham1._43[i][-1] + spinham2._43[i][-1]
        )

    assert np.allclose(m_spinham._43[2][-1], spinham1._43[2][-1])
    assert np.allclose(m_spinham._43[3][-1], spinham2._43[2][-1])


@given(
    ARRAY,
    st.integers(min_value=1, max_value=5),
    st.integers(min_value=1, max_value=5),
    st.integers(min_value=1, max_value=5),
)
def test_make_supercell(parameter1, i, j, k):
    atoms = dict(
        names=["Cr1", "Cr2", "Cr3"],
        spins=[1, 2, 3],
        positions=[[0, 0, 0], [0.33, 0.33, 0.33], [0.66, 0.66, 0.66]],
        g_factors=[2, 2, 2],
    )

    spinham = SpinHamiltonian(cell=np.eye(3), atoms=atoms, convention=CONVENTION)

    spinham.add_43(
        alpha=0, beta=1, gamma=2, nu=(0, 0, 0), _lambda=(0, 0, 0), parameter=parameter1
    )
    spinham.add_43(
        alpha=0,
        beta=1,
        gamma=2,
        nu=(1, 0, 0),
        _lambda=(0, 1, 0),
        parameter=parameter1 * 1.42,
    )

    new_spinham = make_supercell(spinham=spinham, supercell=(i, j, k))

    assert len(new_spinham.p43) == i * j * k * len(spinham.p43)
