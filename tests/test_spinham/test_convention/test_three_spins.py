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


@given(st.floats(allow_nan=False), st.floats(allow_nan=False))
def test_c31(value, new_value):
    convention = Convention(c31=value)

    assert convention.c31 == value

    with pytest.raises(AttributeError):
        convention.c31 = new_value


@given(st.floats(allow_nan=False))
def test_eq_c31(c31):
    convention1 = Convention(c31=c31)
    convention2 = Convention(c31=c31)
    convention3 = Convention(c31=3 if c31 != 3 else 4)
    convention4 = Convention()

    assert convention1 == convention2
    assert convention1 != convention3
    assert convention1 != convention4
    assert convention2 != convention3
    assert convention2 != convention4
    assert convention3 != convention4


def test_get_modified_c31():
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

    mod_convention = convention.get_modified(c31=2)

    assert mod_convention.spin_normalized
    assert mod_convention.multiple_counting
    assert mod_convention.c1 == 1
    assert mod_convention.c21 == 1
    assert mod_convention.c22 == 1
    assert mod_convention.c31 == 2
    assert mod_convention.c32 == 1
    assert mod_convention.c33 == 1
    assert mod_convention.c41 == 1
    assert mod_convention.c421 == 1
    assert mod_convention.c422 == 1
    assert mod_convention.c43 == 1
    assert mod_convention.c44 == 1


@given(st.floats(allow_nan=False), st.floats(allow_nan=False))
def test_c32(value, new_value):
    convention = Convention(c32=value)

    assert convention.c32 == value

    with pytest.raises(AttributeError):
        convention.c32 = new_value


@given(st.floats(allow_nan=False))
def test_eq_c32(c32):
    convention1 = Convention(c32=c32)
    convention2 = Convention(c32=c32)
    convention3 = Convention(c32=3 if c32 != 3 else 4)
    convention4 = Convention()

    assert convention1 == convention2
    assert convention1 != convention3
    assert convention1 != convention4
    assert convention2 != convention3
    assert convention2 != convention4
    assert convention3 != convention4


def test_get_modified_c32():
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

    mod_convention = convention.get_modified(c32=2)

    assert mod_convention.spin_normalized
    assert mod_convention.multiple_counting
    assert mod_convention.c1 == 1
    assert mod_convention.c21 == 1
    assert mod_convention.c22 == 1
    assert mod_convention.c31 == 1
    assert mod_convention.c32 == 2
    assert mod_convention.c33 == 1
    assert mod_convention.c41 == 1
    assert mod_convention.c421 == 1
    assert mod_convention.c422 == 1
    assert mod_convention.c43 == 1
    assert mod_convention.c44 == 1


@given(st.floats(allow_nan=False), st.floats(allow_nan=False))
def test_c33(value, new_value):
    convention = Convention(c33=value)

    assert convention.c33 == value

    with pytest.raises(AttributeError):
        convention.c33 = new_value


@given(st.floats(allow_nan=False))
def test_eq_c33(c33):
    convention1 = Convention(c33=c33)
    convention2 = Convention(c33=c33)
    convention3 = Convention(c33=3 if c33 != 3 else 4)
    convention4 = Convention()

    assert convention1 == convention2
    assert convention1 != convention3
    assert convention1 != convention4
    assert convention2 != convention3
    assert convention2 != convention4
    assert convention3 != convention4


def test_get_modified_c33():
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

    mod_convention = convention.get_modified(c33=2)

    assert mod_convention.spin_normalized
    assert mod_convention.multiple_counting
    assert mod_convention.c1 == 1
    assert mod_convention.c21 == 1
    assert mod_convention.c22 == 1
    assert mod_convention.c31 == 1
    assert mod_convention.c32 == 1
    assert mod_convention.c33 == 2
    assert mod_convention.c41 == 1
    assert mod_convention.c421 == 1
    assert mod_convention.c422 == 1
    assert mod_convention.c43 == 1
    assert mod_convention.c44 == 1
