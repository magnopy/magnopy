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
from hypothesis import strategies as st
from hypothesis.extra.numpy import arrays as harrays

from magnopy import Convention, SpinHamiltonian, converter22

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
RANDOM_UC = harrays(int, (4, 3), elements=st.integers(min_value=-1000, max_value=1000))


def test_multiple_counting():
    spinham = SpinHamiltonian(
        cell=np.eye(3),
        atoms={"names": ["Cr1", "Cr2"]},
        convention=Convention(multiple_counting=True),
    )

    spinham.add_21(0, np.eye(3))
    spinham.add_22(0, 1, (0, 0, 0), 1 * np.eye(3))
    spinham.add_22(0, 1, (1, 0, 0), 2 * np.eye(3))
    spinham.add_22(1, 0, (0, -1, 0), 2 * np.eye(3))

    assert len(spinham._22) == 3
    assert len(spinham.p22) == 6

    assert len(spinham._21) == 1
    assert len(spinham.p21) == 1

    params = list(spinham.p22)

    assert (params[0][3] == params[3][3]).all()
    assert (params[1][3] == params[4][3]).all()
    assert (params[2][3] == params[5][3]).all()

    spinham.convention = Convention(multiple_counting=False)

    assert len(spinham._22) == 3
    assert len(spinham.p22) == 3

    assert len(spinham._21) == 1
    assert len(spinham.p21) == 1


def test_magnetic_atoms():
    spinham = SpinHamiltonian(
        cell=np.eye(3),
        atoms={"names": ["Cr1", "Cr2"]},
        convention=Convention(multiple_counting=True),
    )

    spinham.add_21(0, np.eye(3))
    assert len(spinham.magnetic_atoms.names) == 1
    assert spinham.M == 1
    assert spinham.magnetic_atoms.names[0] == "Cr1"

    spinham.add_22(0, 1, (0, 0, 0), 1 * np.eye(3))
    assert len(spinham.magnetic_atoms.names) == 2
    assert spinham.M == 2
    assert spinham.magnetic_atoms.names[0] == "Cr1"
    assert spinham.magnetic_atoms.names[1] == "Cr2"

    spinham.add_22(0, 1, (1, 0, 0), 2 * np.eye(3))
    assert len(spinham.magnetic_atoms.names) == 2
    assert spinham.M == 2
    assert spinham.magnetic_atoms.names[0] == "Cr1"
    assert spinham.magnetic_atoms.names[1] == "Cr2"

    spinham.add_22(1, 0, (0, -1, 0), 2 * np.eye(3))
    assert len(spinham.magnetic_atoms.names) == 2
    assert spinham.M == 2
    assert spinham.magnetic_atoms.names[0] == "Cr1"
    assert spinham.magnetic_atoms.names[1] == "Cr2"


def test_convention_manipulation():
    atoms = dict(names=["Cr"], postions=[[0, 0, 0]], spins=[3 / 2])

    spinham = SpinHamiltonian(
        cell=np.eye(3), atoms=atoms, convention=Convention(True, False, c21=1, c22=0.5)
    )
    spinham.add_21(0, np.eye(3))
    spinham.add_22(0, 0, (1, 0, 0), np.eye(3))

    assert len(spinham.p22) == 2
    assert converter22.to_iso(parameter=list(spinham.p22)[0][3]) == 1
    assert converter22.to_iso(parameter=spinham.p21[0][1]) == 1

    assert spinham.convention.multiple_counting

    spinham.convention = Convention(False, False, c21=1, c22=0.5)
    assert len(spinham.p22) == 1
    assert converter22.to_iso(parameter=list(spinham.p22)[0][3]) == 2
    assert converter22.to_iso(parameter=spinham.p21[0][1]) == 1

    spinham.convention = Convention(True, False, c21=1, c22=0.5)
    assert len(spinham.p22) == 2
    assert converter22.to_iso(parameter=list(spinham.p22)[0][3]) == 1
    assert converter22.to_iso(parameter=spinham.p21[0][1]) == 1

    assert not spinham.convention.spin_normalized

    spinham.convention = Convention(True, True, c21=1, c22=0.5)
    assert converter22.to_iso(parameter=list(spinham.p22)[0][3]) == 9 / 4
    assert converter22.to_iso(parameter=spinham.p21[0][1]) == 9 / 4

    spinham.convention = Convention(True, False, c21=1, c22=0.5)
    assert converter22.to_iso(parameter=list(spinham.p22)[0][3]) == 1
    assert converter22.to_iso(parameter=spinham.p21[0][1]) == 1

    assert spinham.convention.c22 == 0.5

    spinham.convention = Convention(True, False, c21=1, c22=-1 / 2)
    assert converter22.to_iso(parameter=list(spinham.p22)[0][3]) == -1
    assert converter22.to_iso(parameter=spinham.p21[0][1]) == 1

    spinham.convention = Convention(True, False, c21=1, c22=-1)
    assert converter22.to_iso(parameter=list(spinham.p22)[0][3]) == -0.5
    assert converter22.to_iso(parameter=spinham.p21[0][1]) == 1

    spinham.convention = Convention(True, False, c21=1, c22=-2)
    assert converter22.to_iso(parameter=list(spinham.p22)[0][3]) == -0.25
    assert converter22.to_iso(parameter=spinham.p21[0][1]) == 1

    spinham.convention = Convention(True, False, c21=1, c22=-1 / 2)
    assert converter22.to_iso(parameter=list(spinham.p22)[0][3]) == -1
    assert converter22.to_iso(parameter=spinham.p21[0][1]) == 1

    spinham.convention = Convention(True, False, c21=1, c22=-2)
    assert converter22.to_iso(parameter=list(spinham.p22)[0][3]) == -0.25
    assert converter22.to_iso(parameter=spinham.p21[0][1]) == 1

    spinham.convention = Convention(True, False, c21=1, c22=-1)
    assert converter22.to_iso(parameter=list(spinham.p22)[0][3]) == -0.5
    assert converter22.to_iso(parameter=spinham.p21[0][1]) == 1

    assert spinham.convention.c22 == -1

    spinham.convention = Convention(True, False, c21=1, c22=1)
    assert converter22.to_iso(parameter=list(spinham.p22)[0][3]) == 0.5
    assert converter22.to_iso(parameter=spinham.p21[0][1]) == 1

    spinham.convention = Convention(True, False, c21=1, c22=-1)
    assert converter22.to_iso(parameter=list(spinham.p22)[0][3]) == -0.5
    assert converter22.to_iso(parameter=spinham.p21[0][1]) == 1

    spinham.convention = Convention(True, False, c21=-1, c22=-1)
    assert converter22.to_iso(parameter=list(spinham.p22)[0][3]) == -0.5
    assert converter22.to_iso(parameter=spinham.p21[0][1]) == -1

    spinham.convention = Convention(True, False, c21=0.5, c22=-1)
    assert converter22.to_iso(parameter=list(spinham.p22)[0][3]) == -0.5
    assert converter22.to_iso(parameter=spinham.p21[0][1]) == 2

    spinham.convention = Convention(True, False, c21=1, c22=-1)

    spinham.convention = Convention.get_predefined("spinw")
    assert converter22.to_iso(parameter=list(spinham.p22)[0][3]) == 0.5
    assert converter22.to_iso(parameter=spinham.p21[0][1]) == 1

    spinham.convention = Convention.get_predefined("tb2j")
    assert converter22.to_iso(parameter=list(spinham.p22)[0][3]) == -9 / 8
    assert converter22.to_iso(parameter=spinham.p21[0][1]) == -9 / 4

    spinham.convention = Convention.get_predefined("vampire")
    assert converter22.to_iso(parameter=list(spinham.p22)[0][3]) == -9 / 4
    assert converter22.to_iso(parameter=spinham.p21[0][1]) == -9 / 4


def get_empty():
    cell = 4.2 * np.eye(3)
    atoms = dict(
        names=[f"Cr{i + 1}" for i in range(9)],
        g_factors=[2 for _ in range(9)],
        positions=[[i, i / 2, 0] for i in range(9)],
        spins=[3 / 2 for _ in range(9)],
    )
    convention = Convention(
        multiple_counting=True, spin_normalized=False, c1=1, c22=0.5, c21=1
    )

    spinham = SpinHamiltonian(cell=cell, atoms=atoms, convention=convention)

    spinham.add_21(alpha=3, parameter=np.eye(3))
    spinham.add_22(alpha=3, beta=5, nu=[1, 0, 0], parameter=np.eye(3))
    spinham.add_22(alpha=0, beta=1, nu=[0, 0, 0], parameter=np.eye(3))

    assert spinham.M == 4
    assert len(spinham.p21) == 1
    assert len(spinham.p22) == 4

    empty_spinham = spinham.get_empty()

    assert empty_spinham.M == 0
    assert len(empty_spinham.p1) == 0
    assert len(empty_spinham.p21) == 0
    assert len(empty_spinham.p22) == 0
    assert len(empty_spinham.p31) == 0
    assert len(empty_spinham.p32) == 0
    assert len(empty_spinham.p33) == 0
    assert len(empty_spinham.p41) == 0
    assert len(empty_spinham.p421) == 0
    assert len(empty_spinham.p422) == 0
    assert len(empty_spinham.p43) == 0
    assert len(empty_spinham.p44) == 0
