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

from magnopy import Convention, SpinHamiltonian, make_supercell

MAX_MODULUS = 1e8
RANDOM_CELL = harrays(
    np.float64,
    (3, 3),
    elements=st.floats(min_value=-MAX_MODULUS, max_value=MAX_MODULUS),
)


@given(
    RANDOM_CELL,
    st.integers(min_value=1, max_value=5),
    st.integers(min_value=1, max_value=5),
    st.integers(min_value=1, max_value=5),
)
def test_cell(cell, i, j, k):
    atoms = dict(
        names=["Cr"],
        spins=[1],
        positions=[[0, 0, 0]],
        g_factors=[2],
    )

    spinham = SpinHamiltonian(cell=cell, atoms=atoms, convention=Convention(c1=1))

    new_spinham = make_supercell(spinham=spinham, supercell=(i, j, k))

    assert np.allclose(new_spinham.cell, [i * cell[0], j * cell[1], k * cell[2]])


@given(
    st.integers(min_value=1, max_value=5),
    st.integers(min_value=1, max_value=5),
    st.integers(min_value=1, max_value=5),
    st.integers(min_value=1, max_value=5),
)
def test_atoms(i, j, k, N):
    atoms = dict(
        names=[f"Cr{alpha + 1}" for alpha in range(N)],
        spins=[(1 + alpha) / 2 for alpha in range(N)],
        positions=[[alpha / N, 0, 0] for alpha in range(N)],
        g_factors=[alpha + 1 for alpha in range(N)],
    )

    spinham = SpinHamiltonian(cell=np.eye(3), atoms=atoms, convention=Convention())

    new_spinham = make_supercell(spinham=spinham, supercell=(i, j, k))

    assert len(new_spinham.atoms.names) == i * j * k * N
    assert len(new_spinham.atoms.positions) == i * j * k * N
    assert len(new_spinham.atoms.g_factors) == i * j * k * N
    assert len(new_spinham.atoms.spins) == i * j * k * N

    for kk in range(k):
        for jj in range(j):
            for ii in range(i):
                for alpha in range(N):
                    assert (
                        new_spinham.atoms.names[alpha + N * (ii + jj * i + kk * j * i)]
                        == f"Cr{alpha + 1}_{ii}_{jj}_{kk}"
                    )

                    assert (
                        new_spinham.atoms.spins[alpha + N * (ii + jj * i + kk * j * i)]
                        == (alpha + 1) / 2
                    )

                    assert (
                        new_spinham.atoms.g_factors[
                            alpha + N * (ii + jj * i + kk * j * i)
                        ]
                        == alpha + 1
                    )

                    new_position = (
                        spinham.atoms.positions[alpha]
                        + ii * spinham.cell[0]
                        + jj * spinham.cell[1]
                        + kk * spinham.cell[2]
                    )

                    new_position[0] /= i
                    new_position[1] /= j
                    new_position[2] /= k
                    assert np.allclose(
                        new_spinham.atoms.positions[
                            alpha + N * (ii + jj * i + kk * j * i)
                        ],
                        new_position,
                    )
