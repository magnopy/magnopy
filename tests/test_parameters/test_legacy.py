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

from magnopy import converter22


def test_legacy_test_1():
    full_parameter = converter22.from_iso(iso=23)
    assert converter22.to_iso(parameter=full_parameter) == 23
    assert (
        converter22.to_symm_anisotropy(parameter=full_parameter)
        == np.zeros((3, 3), dtype=float)
    ).all()
    assert (
        converter22.to_dmi(parameter=full_parameter) == np.zeros(3, dtype=float)
    ).all()


def test_legacy_test_2():
    full_parameter = converter22.from_iso(iso=23) + converter22.from_dmi(dmi=(1, 1, 1))
    assert converter22.to_iso(parameter=full_parameter) == 23
    assert (
        converter22.to_symm_anisotropy(parameter=full_parameter)
        == np.zeros((3, 3), dtype=float)
    ).all()
    assert (
        converter22.to_dmi(parameter=full_parameter) == np.ones(3, dtype=float)
    ).all()


def test_legacy_test_3():
    full_parameter = (
        converter22.from_iso(iso=23)
        + converter22.from_dmi(dmi=(1, 1, 1))
        + [[1, 1, 0], [1, -0.5, 0], [0, 0, -0.5]]
    )

    assert converter22.to_iso(parameter=full_parameter) == 23
    assert (
        converter22.to_symm_anisotropy(parameter=full_parameter)
        == np.array([[1, 1, 0], [1, -0.5, 0], [0, 0, -0.5]])
    ).all()
    assert (
        converter22.to_dmi(parameter=full_parameter) == np.ones(3, dtype=float)
    ).all()


def test_legacy_test_4():
    parameter1 = np.array([[1, 5, 2], [5, 8, 4], [2, 6, 3]])
    parameter2 = np.array([[1, 2, 0], [1, 1, 0], [0, 0, 1]])

    assert (
        converter22.to_dmi(parameter=parameter1, matrix_form=True)
        == np.array([[0, 0, 0], [0, 0, -1], [0, 1, 0]])
    ).all()

    assert (
        converter22.to_dmi(parameter=parameter2, matrix_form=True)
        == np.array([[0, 0.5, 0], [-0.5, 0, 0], [0, 0, 0]])
    ).all()


def test_legacy_test_5():
    parameter = np.array([[1, 5, 2], [5, 8, 4], [2, 6, 3]])
    assert converter22.to_iso(parameter=parameter) == 4

    assert (
        converter22.to_iso(parameter=parameter, matrix_form=True)
        == np.array([[4, 0, 0], [0, 4, 0], [0, 0, 4]])
    ).all()

    assert (
        converter22.to_symm_anisotropy(parameter=parameter)
        == np.array([[-3, 5, 2], [5, 4, 5], [2, 5, -1]])
    ).all()

    assert (converter22.to_dmi(parameter=parameter) == np.array([-1, 0, 0])).all()

    assert (
        converter22.to_dmi(parameter=parameter, matrix_form=True)
        == np.array([[0, 0, 0], [0, 0, -1], [0, 1, 0]])
    ).all()
