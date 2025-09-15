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

from magnopy import Convention, Energy, SpinHamiltonian

MAX_MODULUS = 1e4
ARRAY = harrays(
    np.float64,
    (3, 3),
    elements=st.floats(min_value=-MAX_MODULUS, max_value=MAX_MODULUS),
)


@given(ARRAY)
def test_energy_p22(parameter):
    cell = np.eye(3, dtype=float)

    atoms = dict(
        names=["A1", "A2"],
        positions=[[0, 0, 0], [0.5, 0.5, 0.5]],
        spins=[0.5, 2],
        g_factors=[2, 2],
    )
    convention = Convention(multiple_counting=True, spin_normalized=False, c22=1)

    spinham = SpinHamiltonian(cell=cell, atoms=atoms, convention=convention)

    spinham.add_22(alpha=0, beta=1, nu=(0, 0, 0), parameter=parameter)

    energy = Energy(spinham=spinham)

    for i1, sd1 in enumerate([[1, 0, 0], [0, 1, 0], [0, 0, 1]]):
        for i2, sd2 in enumerate([[1, 0, 0], [0, 1, 0], [0, 0, 1]]):
            assert abs(energy.E_0([sd1, sd2]) - 2 * parameter[i1][i2]) < 1e-8


@given(ARRAY)
def test_gradient_p22(parameter):
    cell = np.eye(3, dtype=float)

    atoms = dict(
        names=["A1", "A2"],
        positions=[[0, 0, 0], [0.5, 0.5, 0.5]],
        spins=[0.5, 2],
        g_factors=[2, 2],
    )
    convention = Convention(multiple_counting=True, spin_normalized=False, c22=1)

    spinham = SpinHamiltonian(cell=cell, atoms=atoms, convention=convention)

    spinham.add_22(alpha=0, beta=1, nu=(0, 0, 0), parameter=parameter)

    energy = Energy(spinham=spinham)
    h = 1e-3

    for sd1 in [[1, 0, 0], [0, 1, 0], [0, 0, 1], np.array([1, 1, 1]) / np.sqrt(3)]:
        for sd2 in [[1, 0, 0], [0, 1, 0], [0, 0, 1], np.array([1, 1, 1]) / np.sqrt(3)]:
            numerical_derivative = [[0, 0, 0], [0, 0, 0]]
            derivative = energy.gradient(spin_directions=[sd1, sd2])
            for i in range(3):
                sd1[i] += h
                energy_plus = energy.E_0(spin_directions=[sd1, sd2], _normalize=False)
                sd1[i] -= 2 * h
                energy_minus = energy.E_0(spin_directions=[sd1, sd2], _normalize=False)
                sd1[i] += h

                numerical_derivative[0][i] = (energy_plus - energy_minus) / 2 / h

                sd2[i] += h
                energy_plus = energy.E_0(spin_directions=[sd1, sd2], _normalize=False)
                sd2[i] -= 2 * h
                energy_minus = energy.E_0(spin_directions=[sd1, sd2], _normalize=False)
                sd2[i] += h

                numerical_derivative[1][i] = (energy_plus - energy_minus) / 2 / h

            assert np.allclose(derivative, numerical_derivative)
