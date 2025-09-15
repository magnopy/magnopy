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

from magnopy import Convention, Energy, SpinHamiltonian


def test_energy_p1():
    cell = np.eye(3, dtype=float)

    atoms = dict(
        names=["A1", "A2"],
        positions=[[0, 0, 0], [0.5, 0.5, 0.5]],
        spins=[1, 1],
        g_factors=[2, 2],
    )
    convention = Convention(multiple_counting=True, spin_normalized=False, c1=1, c22=1)

    spinham = SpinHamiltonian(cell=cell, atoms=atoms, convention=convention)

    BOHR_MAGNETON = 0.057883818060  # meV / Tesla

    spinham.add_magnetic_field(h=[0, 0, 1], alphas=[0])

    spinham.add_22(alpha=0, beta=1, nu=(0, 0, 0), parameter=np.zeros((3, 3)))

    assert len(spinham.p1) == 1

    energy = Energy(spinham=spinham)

    assert abs(energy.E_0([[0, 0, 1]]) - 2 * BOHR_MAGNETON) < 1e-8
