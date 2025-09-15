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


import pytest
from hypothesis import given
from hypothesis import strategies as st

from magnopy import Convention, ConventionError


def tes_raises():
    convention = Convention()

    with pytest.raises(ConventionError):
        convention.multiple_counting

    with pytest.raises(ConventionError):
        convention.spin_normalized

    with pytest.raises(ConventionError):
        convention.c1

    with pytest.raises(ConventionError):
        convention.c21

    with pytest.raises(ConventionError):
        convention.c22

    with pytest.raises(ConventionError):
        convention.c31

    with pytest.raises(ConventionError):
        convention.c32

    with pytest.raises(ConventionError):
        convention.c33

    with pytest.raises(ConventionError):
        convention.c41

    with pytest.raises(ConventionError):
        convention.c421

    with pytest.raises(ConventionError):
        convention.c422

    with pytest.raises(ConventionError):
        convention.c43

    with pytest.raises(ConventionError):
        convention.c44


@given(st.text())
def test_name(name):
    convention = Convention()

    assert convention.name == "custom"

    convention.name = name

    assert convention.name == name.lower()


def test_get_predefined():
    _ = Convention.get_predefined(name="tb2j")
    _ = Convention.get_predefined(name="vampire")
    _ = Convention.get_predefined(name="spinW")
