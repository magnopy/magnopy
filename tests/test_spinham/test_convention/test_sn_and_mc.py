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

from magnopy import Convention


################################################################################
#                                Spin normalized                               #
################################################################################
@given(st.booleans(), st.booleans())
def test_spin_normalized(value, new_value):
    convention = Convention(spin_normalized=value)

    assert convention.spin_normalized == value

    with pytest.raises(AttributeError):
        convention.spin_normalized = new_value


@given(st.booleans())
def test_eq_spin_normalized(sn):
    convention1 = Convention(spin_normalized=sn)
    convention2 = Convention(spin_normalized=sn)
    convention3 = Convention(spin_normalized=not sn)
    convention4 = Convention()

    assert convention1 == convention2
    assert convention1 != convention3
    assert convention1 != convention4
    assert convention2 != convention3
    assert convention2 != convention4
    assert convention3 != convention4


def test_get_modified_spin_normalized():
    convention = Convention(
        spin_normalized=True,
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

    mod_convention = convention.get_modified(spin_normalized=False)

    assert not mod_convention.spin_normalized
    assert mod_convention.multiple_counting
    assert mod_convention.c1 == 1
    assert mod_convention.c21 == 1
    assert mod_convention.c22 == 1
    assert mod_convention.c31 == 1
    assert mod_convention.c32 == 1
    assert mod_convention.c33 == 1
    assert mod_convention.c41 == 1
    assert mod_convention.c421 == 1
    assert mod_convention.c422 == 1
    assert mod_convention.c43 == 1
    assert mod_convention.c44 == 1


################################################################################
#                               Multiple counting                              #
################################################################################
@given(st.booleans(), st.booleans())
def test_multiple_counting(value, new_value):
    convention = Convention(multiple_counting=value)

    assert convention.multiple_counting == value

    with pytest.raises(AttributeError):
        convention.multiple_counting = new_value


@given(st.booleans())
def test_eq_multiple_counting(mc):
    convention1 = Convention(multiple_counting=mc)
    convention2 = Convention(multiple_counting=mc)
    convention3 = Convention(multiple_counting=not mc)
    convention4 = Convention()

    assert convention1 == convention2
    assert convention1 != convention3
    assert convention1 != convention4
    assert convention2 != convention3
    assert convention2 != convention4
    assert convention3 != convention4


def test_get_modified_multiple_counting():
    convention = Convention(
        spin_normalized=True,
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

    mod_convention = convention.get_modified(multiple_counting=False)

    assert mod_convention.spin_normalized
    assert not mod_convention.multiple_counting
    assert mod_convention.c1 == 1
    assert mod_convention.c21 == 1
    assert mod_convention.c22 == 1
    assert mod_convention.c31 == 1
    assert mod_convention.c32 == 1
    assert mod_convention.c33 == 1
    assert mod_convention.c41 == 1
    assert mod_convention.c421 == 1
    assert mod_convention.c422 == 1
    assert mod_convention.c43 == 1
    assert mod_convention.c44 == 1
