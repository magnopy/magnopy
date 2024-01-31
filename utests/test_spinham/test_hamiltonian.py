# MAGNOPY - Python package for magnons.
# Copyright (C) 2023-2024 Magnopy Team
#
# e-mail: anry@uv.es, web: magnopy.com
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

from math import sqrt

import numpy as np
import pytest
from wulfric import ABS_TOL, REL_TOL, Atom, compare_numerically

from magnopy.exceptions import NotationError
from magnopy.spinham.hamiltonian import SpinHamiltonian
from magnopy.spinham.parameter import MatrixParameter


################################################################################
#                                 Legacy Tests                                 #
################################################################################
def test_iteration():
    model = SpinHamiltonian()
    Cr1 = Atom("Cr1", (0.25, 0.25, 0))
    Cr2 = Atom("Cr2", (0.75, 0.75, 0))
    model.add_atom(Cr1)
    model.add_atom(Cr2)
    bonds = [
        (12, Cr1, Cr2, (0, 0, 0)),
        (2, Cr2, Cr1, (0, 0, 0)),
        (6, Cr1, Cr1, (1, 0, 0)),
        (3.425, Cr1, Cr1, (-1, 0, 0)),
        (5.3, Cr2, Cr2, (1, 0, 0)),
        (7.34, Cr2, Cr2, (-1, 0, 0)),
        (12.4, Cr1, Cr1, (0, 2, 0)),
        (34, Cr1, Cr1, (0, -2, 0)),
        (1.098, Cr2, Cr2, (0, 2, 0)),
        (0.0054, Cr2, Cr2, (0, -2, 0)),
        (0.35, Cr2, Cr1, (2, 2, 0)),
        (-2.35, Cr1, Cr2, (-2, -2, 0)),
    ]
    for iso, atom1, atom2, R in bonds:
        model.add_bond(atom1, atom2, R, iso=iso)
    for i, (atom1, atom2, R, J) in enumerate(model):
        assert isinstance(atom1, Atom)
        assert isinstance(atom2, Atom)
        assert isinstance(R, tuple)
        assert J == MatrixParameter(iso=bonds[i][0])


def test_contains():
    model = SpinHamiltonian()
    Cr1 = Atom("Cr1", (0.25, 0.25, 0))
    Cr2 = Atom("Cr2", (0.75, 0.75, 0))
    model.add_atom(Cr1)
    model.add_atom(Cr2)
    bonds = [
        (12, Cr1, Cr2, (0, 0, 0)),
        (12, Cr2, Cr1, (0, 0, 0)),
        (12, Cr1, Cr1, (1, 0, 0)),
        (12, Cr1, Cr1, (-1, 0, 0)),
        (12, Cr2, Cr2, (1, 0, 0)),
        (12, Cr2, Cr2, (-1, 0, 0)),
        (12, Cr1, Cr1, (0, 2, 0)),
        (12, Cr1, Cr1, (0, -2, 0)),
        (12, Cr2, Cr2, (0, 2, 0)),
        (12, Cr2, Cr2, (0, -2, 0)),
        (12, Cr2, Cr1, (2, 2, 0)),
        (12, Cr1, Cr2, (-2, -2, 0)),
    ]
    for iso, atom1, atom2, R in bonds:
        model.add_bond(atom1, atom2, R, iso=iso)
    assert (Cr1, Cr2, (0, 0, 0)) in model
    assert (Cr2, Cr2, (1, 0, 0)) in model
    assert (Cr1, Cr2, (4, 0, 0)) not in model
    assert (Cr1, Cr2, (0, 8, 0)) not in model


def test_getitem():
    model = SpinHamiltonian()
    Cr1 = Atom("Cr1", (0.25, 0.25, 0))
    Cr2 = Atom("Cr2", (0.75, 0.75, 0))
    model.add_atom(Cr1)
    model.add_atom(Cr2)
    bonds = [
        (12, Cr1, Cr2, (0, 0, 0)),
        (12, Cr2, Cr1, (0, 0, 0)),
        (12, Cr1, Cr1, (1, 0, 0)),
        (12, Cr1, Cr1, (-1, 0, 0)),
        (12, Cr2, Cr2, (1, 0, 0)),
        (12, Cr2, Cr2, (-1, 0, 0)),
        (12, Cr1, Cr1, (0, 2, 0)),
        (12, Cr1, Cr1, (0, -2, 0)),
        (12, Cr2, Cr2, (0, 2, 0)),
        (12, Cr2, Cr2, (0, -2, 0)),
        (12, Cr2, Cr1, (2, 2, 0)),
        (12, Cr1, Cr2, (-2, -2, 0)),
    ]
    for iso, atom1, atom2, R in bonds:
        model.add_bond(atom1, atom2, R, iso=iso)
    assert compare_numerically(model[(Cr1, Cr2, (0, 0, 0))].iso, "==", 12)
    assert compare_numerically(model[(Cr2, Cr2, (1, 0, 0))].iso, "==", 12)
    assert compare_numerically(model[(Cr1, Cr2, (-2, -2, 0))].iso, "==", 12)
    with pytest.raises(KeyError):
        assert compare_numerically(model[(Cr1, Cr2, (0, 8, 0))].iso, "==", 12)


def test_cell_error():
    model = SpinHamiltonian()
    with pytest.raises(ValueError):
        model.cell = [3, 4, 5]


def test_I():
    Cr1 = Atom("Cr1", (1, 6, 2))
    Cr2 = Atom("Cr2", (1, 3, 5))
    Cr3 = Atom("Cr3", (1, 3, 3))
    model = SpinHamiltonian()
    model.add_atom(Cr1)
    assert model.I == 0
    model.add_bond("Cr1", "Cr1", (1, 0, 0), iso=1)
    assert model.I == 1
    model.add_atom(Cr2)
    assert model.I == 1
    model.add_atom(Cr3)
    assert model.I == 1
    model.add_bond("Cr2", "Cr3", (0, 0, 0), iso=1)
    assert model.I == 3
    model.remove_atom(Cr1)
    assert model.I == 2
    model.remove_atom(Cr2)
    assert model.I == 0
    model.remove_atom(Cr3)
    assert model.I == 0


def test_add_bond():
    model = SpinHamiltonian()
    Cr1 = Atom("Cr1", (2, 5, 1))
    Cr2 = Atom("Cr2", (4, 2, 1))
    Cr3 = Atom("Cr3", (5, 1, 8))
    bond12 = MatrixParameter(iso=12)
    bond13 = MatrixParameter(iso=13)
    bond23 = MatrixParameter(iso=23)
    bond31 = MatrixParameter(iso=31)
    model.add_bond(Cr1, Cr2, (0, 0, 0), J=bond12)
    model.add_bond(Cr2, Cr3, (0, 0, 0), J=bond23)
    model.add_bond(Cr3, Cr1, (0, 0, 0), J=bond31)
    assert len(model) == 3
    assert (Cr1, Cr2, (0, 0, 0)) in model
    assert (Cr2, Cr3, (0, 0, 0)) in model
    assert (Cr3, Cr1, (0, 0, 0)) in model
    model.add_bond(Cr1, Cr3, (0, 0, 0), J=bond13)
    assert len(model) == 4
    assert (Cr1, Cr2, (0, 0, 0)) in model
    assert (Cr2, Cr3, (0, 0, 0)) in model
    assert (Cr3, Cr1, (0, 0, 0)) in model
    assert (Cr1, Cr3, (0, 0, 0)) in model


def test_remove_bond():
    model = SpinHamiltonian()
    Cr1 = Atom("Cr1", (2, 5, 1))
    Cr2 = Atom("Cr2", (4, 2, 1))
    Cr3 = Atom("Cr3", (5, 1, 8))
    bond12 = MatrixParameter(iso=12)
    bond13 = MatrixParameter(iso=13)
    bond23 = MatrixParameter(iso=23)
    bond31 = MatrixParameter(iso=31)
    model.add_bond(Cr1, Cr2, (0, 0, 0), J=bond12)
    model.add_bond(Cr2, Cr3, (0, 0, 0), J=bond23)
    model.add_bond(Cr3, Cr1, (0, 0, 0), J=bond31)
    model.add_bond(Cr1, Cr3, (0, 0, 0), J=bond13)
    assert len(model) == 4
    assert (Cr1, Cr2, (0, 0, 0)) in model
    assert (Cr2, Cr3, (0, 0, 0)) in model
    assert (Cr3, Cr1, (0, 0, 0)) in model
    assert (Cr1, Cr3, (0, 0, 0)) in model
    model.remove_bond(Cr1, Cr2, (0, 0, 0))
    assert len(model) == 3
    assert (Cr1, Cr2, (0, 0, 0)) not in model
    assert (Cr2, Cr3, (0, 0, 0)) in model
    assert (Cr3, Cr1, (0, 0, 0)) in model
    assert (Cr1, Cr3, (0, 0, 0)) in model
    model.remove_bond(Cr3, Cr1, (0, 0, 0))
    assert len(model) == 2
    assert (Cr1, Cr2, (0, 0, 0)) not in model
    assert (Cr2, Cr3, (0, 0, 0)) in model
    assert (Cr3, Cr1, (0, 0, 0)) not in model
    assert (Cr1, Cr3, (0, 0, 0)) in model


def test_add_atom():
    model = SpinHamiltonian()
    assert len(model.magnetic_atoms) == 0
    Cr1 = Atom("Cr1", (2, 5, 1))
    model.add_atom(Cr1)
    assert len(model.magnetic_atoms) == 0
    assert len(model.atoms) == 1
    Cr2 = Atom("Cr2", (4, 2, 1))
    model.add_atom(Cr2)
    assert len(model.magnetic_atoms) == 0
    assert len(model.atoms) == 2
    assert compare_numerically(model.get_atom("Cr1").position[0], "==", 2)
    assert compare_numerically(model.get_atom("Cr1").position[1], "==", 5)
    assert compare_numerically(model.get_atom("Cr1").position[2], "==", 1)
    assert compare_numerically(model.get_atom("Cr2").position[0], "==", 4)
    assert compare_numerically(model.get_atom("Cr2").position[1], "==", 2)
    assert compare_numerically(model.get_atom("Cr2").position[2], "==", 1)
    model.remove_atom(Cr2)
    Cr2 = Atom("Cr2", (4, 3, 1))
    model.add_atom(Cr2)
    assert compare_numerically(model.get_atom("Cr2").position[0], "==", 4)
    assert compare_numerically(model.get_atom("Cr2").position[1], "==", 3)
    assert compare_numerically(model.get_atom("Cr2").position[2], "==", 1)


def test_remove_atom():
    model = SpinHamiltonian()
    Cr1 = Atom("Cr1", (2, 5, 1))
    Cr2 = Atom("Cr2", (4, 2, 1))
    Cr3 = Atom("Cr3", (5, 1, 8))
    bond12 = MatrixParameter(iso=12)
    bond13 = MatrixParameter(iso=13)
    bond23 = MatrixParameter(iso=23)
    bond31 = MatrixParameter(iso=31)
    model.add_bond(Cr1, Cr2, (0, 0, 0), J=bond12)
    model.add_bond(Cr2, Cr3, (0, 0, 0), J=bond23)
    model.add_bond(Cr3, Cr1, (0, 0, 0), J=bond31)
    model.add_bond(Cr1, Cr3, (0, 0, 0), J=bond13)
    assert len(model) == 4
    assert (Cr1, Cr2, (0, 0, 0)) in model
    assert (Cr2, Cr3, (0, 0, 0)) in model
    assert (Cr3, Cr1, (0, 0, 0)) in model
    assert (Cr1, Cr3, (0, 0, 0)) in model
    assert len(model.magnetic_atoms) == 3
    model.remove_atom(Cr1)
    assert len(model.magnetic_atoms) == 2
    assert len(model) == 1
    assert (Cr1, Cr2, (0, 0, 0)) not in model
    assert (Cr2, Cr3, (0, 0, 0)) in model
    assert (Cr3, Cr1, (0, 0, 0)) not in model
    assert (Cr1, Cr3, (0, 0, 0)) not in model


def test_get_atom_coordinates():
    model = SpinHamiltonian()
    model.cell = [[10, 0, 0], [0, 10, 0], [0, 0, 10]]
    Cr1 = Atom("Cr1", (0.1, 0.6, 0.2))
    model.add_atom(Cr1)
    x, y, z = model.get_atom_coordinates(Cr1, relative=False)
    assert compare_numerically(x, "==", 1)
    assert compare_numerically(y, "==", 6)
    assert compare_numerically(z, "==", 2)
    x, y, z = model.get_atom_coordinates(Cr1, R=[1, 0, 0], relative=False)
    assert compare_numerically(x, "==", 11)
    assert compare_numerically(y, "==", 6)
    assert compare_numerically(z, "==", 2)
    x, y, z = model.get_atom_coordinates(Cr1, R=[0, -1, 0], relative=False)
    assert compare_numerically(x, "==", 1)
    assert compare_numerically(y, "==", -4)
    assert compare_numerically(z, "==", 2)
    x, y, z = model.get_atom_coordinates(Cr1, R=[0, 0, 2], relative=False)
    assert compare_numerically(x, "==", 1)
    assert compare_numerically(y, "==", 6)
    assert compare_numerically(z, "==", 22)
    x, y, z = model.get_atom_coordinates(Cr1, R=[0, -3, 2], relative=False)
    assert compare_numerically(x, "==", 1)
    assert compare_numerically(y, "==", -24)
    assert compare_numerically(z, "==", 22)
    x, y, z = model.get_atom_coordinates(Cr1, R=[3, -3, 2], relative=False)
    assert compare_numerically(x, "==", 31)
    assert compare_numerically(y, "==", -24)
    assert compare_numerically(z, "==", 22)
    model.cell = [[10, 10, 0], [0, 10, 10], [0, 0, 10]]
    x, y, z = model.get_atom_coordinates(Cr1, relative=False)
    assert compare_numerically(x, "==", 1)
    assert compare_numerically(y, "==", 7)
    assert compare_numerically(z, "==", 8)
    x, y, z = model.get_atom_coordinates(Cr1, R=[1, 0, 0], relative=False)
    assert compare_numerically(x, "==", 11)
    assert compare_numerically(y, "==", 17)
    assert compare_numerically(z, "==", 8)
    x, y, z = model.get_atom_coordinates(Cr1, R=[0, -1, 0], relative=False)
    assert compare_numerically(x, "==", 1)
    assert compare_numerically(y, "==", -3)
    assert compare_numerically(z, "==", -2)
    x, y, z = model.get_atom_coordinates(Cr1, R=[0, 0, 2], relative=False)
    assert compare_numerically(x, "==", 1)
    assert compare_numerically(y, "==", 7)
    assert compare_numerically(z, "==", 28)
    x, y, z = model.get_atom_coordinates(Cr1, R=[0, -3, 2], relative=False)
    assert compare_numerically(x, "==", 1)
    assert compare_numerically(y, "==", -23)
    assert compare_numerically(z, "==", -2)
    x, y, z = model.get_atom_coordinates(Cr1, R=[3, -3, 2], relative=False)
    assert compare_numerically(x, "==", 31)
    assert compare_numerically(y, "==", 7)
    assert compare_numerically(z, "==", -2)


def test_get_bond_vector():
    model = SpinHamiltonian()
    model.cell = [[10, 0, 0], [0, 10, 0], [0, 0, 10]]
    Cr1 = Atom("Cr1", (0.1, 0.6, 0.2))
    Cr2 = Atom("Cr2", (0.1, 0.3, 0.4))
    model.add_atom(Cr1)
    model.add_atom(Cr2)
    vector = model.get_vector(Cr1, Cr2)
    assert np.allclose(vector, np.array([0, -3, 2]))
    vector = model.get_vector(Cr1, Cr2, R=[3, 2, -5])
    assert np.allclose(vector, np.array([30, 17, -48]))


def test_get_distance():
    model = SpinHamiltonian()
    model.cell = [[10, 0, 0], [0, 10, 0], [0, 0, 10]]
    Cr1 = Atom("Cr1", (0.1, 0.6, 0.2))
    Cr2 = Atom("Cr2", (0.1, 0.3, 0.4))
    model.add_atom(Cr1)
    model.add_atom(Cr2)
    d = model.get_distance(Cr1, Cr2)
    assert compare_numerically(d, "==", sqrt(13))
    d = model.get_distance(Cr1, Cr2, R=[3, 2, -5])
    assert compare_numerically(d, "==", sqrt(900 + 17**2 + 48**2))


def test_filter():
    model = SpinHamiltonian()
    Cr1 = Atom("Cr1", (0.25, 0.25, 0))
    Cr2 = Atom("Cr2", (0.75, 0.75, 0))
    bonds = [
        (12, Cr1, Cr2, (0, 0, 0)),
        (12, Cr2, Cr1, (0, 0, 0)),
        (12, Cr1, Cr1, (1, 0, 0)),
        (12, Cr1, Cr1, (-1, 0, 0)),
        (12, Cr2, Cr2, (1, 0, 0)),
        (12, Cr2, Cr2, (-1, 0, 0)),
        (12, Cr1, Cr1, (0, 2, 0)),
        (12, Cr1, Cr1, (0, -2, 0)),
        (12, Cr2, Cr2, (0, 2, 0)),
        (12, Cr2, Cr2, (0, -2, 0)),
        (12, Cr2, Cr1, (2, 2, 0)),
        (12, Cr1, Cr2, (-2, -2, 0)),
    ]
    for iso, atom1, atom2, R in bonds:
        model.add_bond(atom1, atom2, R, iso=iso)
    assert len(model) == 12
    assert (Cr1, Cr2, (0, 0, 0)) in model
    assert (Cr2, Cr1, (0, 0, 0)) in model
    assert (Cr1, Cr1, (1, 0, 0)) in model
    assert (Cr1, Cr1, (-1, 0, 0)) in model
    assert (Cr2, Cr2, (1, 0, 0)) in model
    assert (Cr2, Cr2, (-1, 0, 0)) in model
    assert (Cr1, Cr1, (0, 2, 0)) in model
    assert (Cr1, Cr1, (0, -2, 0)) in model
    assert (Cr2, Cr2, (0, 2, 0)) in model
    assert (Cr2, Cr2, (0, -2, 0)) in model
    assert (Cr2, Cr1, (2, 2, 0)) in model
    assert (Cr1, Cr2, (-2, -2, 0)) in model
    filtered_model = model.filtered(max_distance=1)
    assert len(filtered_model) == 6
    assert (Cr1, Cr2, (0, 0, 0)) in model
    assert (Cr2, Cr1, (0, 0, 0)) in model
    assert (Cr1, Cr1, (1, 0, 0)) in model
    assert (Cr1, Cr1, (-1, 0, 0)) in model
    assert (Cr2, Cr2, (1, 0, 0)) in model
    assert (Cr2, Cr2, (-1, 0, 0)) in model
    filtered_model = model.filtered(min_distance=1)
    assert len(filtered_model) == 10
    assert (Cr1, Cr1, (1, 0, 0)) in model
    assert (Cr1, Cr1, (-1, 0, 0)) in model
    assert (Cr2, Cr2, (1, 0, 0)) in model
    assert (Cr2, Cr2, (-1, 0, 0)) in model
    assert (Cr1, Cr1, (0, 2, 0)) in model
    assert (Cr1, Cr1, (0, -2, 0)) in model
    assert (Cr2, Cr2, (0, 2, 0)) in model
    assert (Cr2, Cr2, (0, -2, 0)) in model
    assert (Cr2, Cr1, (2, 2, 0)) in model
    assert (Cr1, Cr2, (-2, -2, 0)) in model
    filtered_model = model.filtered(min_distance=1, max_distance=2)
    assert len(filtered_model) == 8
    assert (Cr1, Cr1, (1, 0, 0)) in model
    assert (Cr1, Cr1, (-1, 0, 0)) in model
    assert (Cr2, Cr2, (1, 0, 0)) in model
    assert (Cr2, Cr2, (-1, 0, 0)) in model
    assert (Cr1, Cr1, (0, 2, 0)) in model
    assert (Cr1, Cr1, (0, -2, 0)) in model
    assert (Cr2, Cr2, (0, 2, 0)) in model
    assert (Cr2, Cr2, (0, -2, 0)) in model


def test_notation_manipulation():
    model = SpinHamiltonian()
    Cr = Atom("Cr", (0, 0, 0), spin=3 / 2)
    model[Cr, Cr, (1, 0, 0)] = MatrixParameter(iso=1)
    assert compare_numerically(model[Cr, Cr, (1, 0, 0)].iso, "==", 1)

    with pytest.raises(NotationError):
        model.double_counting
    with pytest.raises(NotationError):
        model.spin_normalized
    with pytest.raises(NotationError):
        model.exchange_factor
    with pytest.raises(NotationError):
        model.on_site_factor

    model.notation = "magnopy"
    assert compare_numerically(model[Cr, Cr, (1, 0, 0)].iso, "==", 1)
    assert len(model) == 2
    model.notation = "magnopy"
    assert compare_numerically(model[Cr, Cr, (1, 0, 0)].iso, "==", 1)
    assert len(model) == 2

    assert model.double_counting
    assert len(model) == 2
    model.double_counting = False
    assert compare_numerically(model[Cr, Cr, (1, 0, 0)].iso, "==", 2)
    assert len(model) == 1
    model.double_counting = True
    assert compare_numerically(model[Cr, Cr, (1, 0, 0)].iso, "==", 1)
    assert len(model) == 2

    assert not model.spin_normalized
    model.spin_normalized = True
    assert compare_numerically(model[Cr, Cr, (1, 0, 0)].iso, "==", 9 / 4)
    model.spin_normalized = False
    assert compare_numerically(model[Cr, Cr, (1, 0, 0)].iso, "==", 1)

    assert model.exchange_factor == 0.5
    model.exchange_factor = -1 / 2
    assert compare_numerically(model[Cr, Cr, (1, 0, 0)].iso, "==", -1)
    model.exchange_factor = -1
    assert compare_numerically(model[Cr, Cr, (1, 0, 0)].iso, "==", -0.5)
    model.exchange_factor = -2
    assert compare_numerically(model[Cr, Cr, (1, 0, 0)].iso, "==", -0.25)
    model.exchange_factor = -1 / 2
    assert compare_numerically(model[Cr, Cr, (1, 0, 0)].iso, "==", -1)
    model.exchange_factor = -2
    assert compare_numerically(model[Cr, Cr, (1, 0, 0)].iso, "==", -0.25)
    model.exchange_factor = -1
    assert compare_numerically(model[Cr, Cr, (1, 0, 0)].iso, "==", -0.5)

    assert model.exchange_factor == -1
    model.exchange_factor = 1
    assert compare_numerically(model[Cr, Cr, (1, 0, 0)].iso, "==", 0.5)
    model.exchange_factor = -1
    assert compare_numerically(model[Cr, Cr, (1, 0, 0)].iso, "==", -0.5)

    model.notation = "SpinW"
    assert compare_numerically(model[Cr, Cr, (1, 0, 0)].iso, "==", 0.5)
    model.notation = "TB2J"
    assert compare_numerically(model[Cr, Cr, (1, 0, 0)].iso, "==", -9 / 8)
    model.notation = "Vampire"
    assert compare_numerically(model[Cr, Cr, (1, 0, 0)].iso, "==", -9 / 4)


def test_predefined_notations():
    model = SpinHamiltonian()
    Cr = Atom("Cr", (0, 0, 0), spin=3 / 2)
    model[Cr, Cr, (1, 0, 0)] = MatrixParameter(iso=1)
    assert compare_numerically(model[Cr, Cr, (1, 0, 0)].iso, "==", 1)

    with pytest.raises(NotationError):
        model.double_counting
    with pytest.raises(NotationError):
        model.spin_normalized
    with pytest.raises(NotationError):
        model.exchange_factor
    with pytest.raises(NotationError):
        model.on_site_factor

    model.notation = "magnopy"
    assert compare_numerically(model[Cr, Cr, (1, 0, 0)].iso, "==", 1)
    assert (
        model.double_counting
        and not model.spin_normalized
        and model.exchange_factor == 0.5
        and model.on_site_factor == 1
    )

    model.notation = "TB2J"
    assert compare_numerically(model[Cr, Cr, (1, 0, 0)].iso, "==", -9 / 8)
    assert (
        model.double_counting
        and model.spin_normalized
        and model.exchange_factor == -1
        and model.on_site_factor == -1
    )

    model.notation = "SpinW"
    assert compare_numerically(model[Cr, Cr, (1, 0, 0)].iso, "==", 0.5)
    assert (
        model.double_counting
        and not model.spin_normalized
        and model.exchange_factor == 1
        and model.on_site_factor == 1
    )

    model.notation = "Vampire"
    assert compare_numerically(model[Cr, Cr, (1, 0, 0)].iso, "==", -9 / 4)
    assert (
        model.double_counting
        and model.spin_normalized
        and model.exchange_factor == -0.5
        and model.on_site_factor == -1
    )


def test_add_remove_bond_with_notation():
    model = SpinHamiltonian()
    Cr1 = Atom("Cr1", (2, 5, 1))
    Cr2 = Atom("Cr2", (4, 2, 1))
    bond = MatrixParameter(iso=1)
    model.double_counting = False
    model.add_bond(Cr1, Cr2, (0, 0, 0), J=bond)
    assert (Cr1, Cr2, (0, 0, 0)) in model
    assert (Cr2, Cr1, (0, 0, 0)) not in model
    model.remove_bond(Cr1, Cr2, (0, 0, 0))
    assert (Cr1, Cr2, (0, 0, 0)) not in model
    assert (Cr2, Cr1, (0, 0, 0)) not in model

    model.double_counting = True
    model.add_bond(Cr1, Cr2, (0, 0, 0), J=bond)
    assert (Cr1, Cr2, (0, 0, 0)) in model
    assert (Cr2, Cr1, (0, 0, 0)) in model
    model.remove_bond(Cr2, Cr1, (0, 0, 0))
    assert (Cr1, Cr2, (0, 0, 0)) not in model
    assert (Cr2, Cr1, (0, 0, 0)) not in model
