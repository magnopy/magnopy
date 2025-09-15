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
from hypothesis import given, settings
from hypothesis import strategies as st
from wulfric.cell import sc_get_example_cell

from magnopy import Convention, SpinHamiltonian

LATTICE_VARIATIONS = [
    "cub",
    "fcc",
    "bcc",
    "tet",
    "bct1",
    "bct2",
    "orc",
    "orcf1",
    "orcf2",
    "orcf3",
    "orci",
    "orcc",
    "hex",
    "rhl1",
    "rhl2",
    "mcl",
    "mclc1",
    "mclc2",
    "mclc3",
    "mclc4",
    "mclc5",
    # Temporary - FIXME later
    # "tri1a",
    # "tri2a",
    # "tri1b",
    # "tri2b",
]


def test_raises():
    cell = np.eye(3)
    atoms = dict(
        names=["Cr1", "Cr2"],
        spins=[3 / 2, 3 / 2],
        positions=[[0, 0, 0], [0.5, 0.5, 0.5]],
        g_factors=[2, 2],
    )
    convention = Convention(multiple_counting=True, spin_normalized=False, c22=1)

    spinham = SpinHamiltonian(cell=cell, atoms=atoms, convention=convention)

    with pytest.raises(ValueError):
        spinham.add_dipole_dipole()

    with pytest.raises(ValueError):
        spinham.add_dipole_dipole(R_cut=-10)

    with pytest.raises(ValueError):
        spinham.add_dipole_dipole(E_cut=-10)

    with pytest.raises(ValueError):
        spinham.add_dipole_dipole(E_cut=0)


@pytest.mark.parametrize("lattice_variation", LATTICE_VARIATIONS)
@given(
    st.floats(min_value=0, max_value=7),
)
@settings(deadline=1000)
def test_R_cut(lattice_variation, R_cut):
    cell = sc_get_example_cell(lattice_variation=lattice_variation)
    atoms = dict(
        names=["Cr1", "Cr2"],
        spins=[3 / 2, 3 / 2],
        positions=[[0, 0, 0], [0.5, 0.5, 0.5]],
        g_factors=[2, 2],
    )
    convention = Convention(multiple_counting=True, spin_normalized=False, c22=1)

    spinham = SpinHamiltonian(cell=cell, atoms=atoms, convention=convention)

    spinham_double = spinham.copy()

    spinham.add_dipole_dipole(
        R_cut=R_cut, alphas=[i for i in range(len(spinham.atoms.names))]
    )

    spinham_double.add_dipole_dipole(
        R_cut=2 * R_cut, alphas=[i for i in range(len(spinham_double.atoms.names))]
    )

    if len(spinham.p22) > 0:
        assert len(spinham.p22) < len(spinham_double.p22)


@pytest.mark.parametrize("lattice_variation", LATTICE_VARIATIONS)
@given(
    st.floats(min_value=0.01, max_value=7),
)
@settings(deadline=1000)
def test_E_cut(lattice_variation, E_cut):
    cell = sc_get_example_cell(lattice_variation=lattice_variation)
    atoms = dict(
        names=["Cr1", "Cr2"],
        spins=[3 / 2, 3 / 2],
        positions=[[0, 0, 0], [0.5, 0.5, 0.5]],
        g_factors=[2, 2],
    )
    convention = Convention(multiple_counting=True, spin_normalized=False, c22=1)

    spinham = SpinHamiltonian(cell=cell, atoms=atoms, convention=convention)

    spinham.add_dipole_dipole(
        E_cut=E_cut, alphas=[i for i in range(len(spinham.atoms.names))]
    )
