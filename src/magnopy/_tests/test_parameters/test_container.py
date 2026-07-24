# ================================== LICENSE ===================================
# Magnopy - Python package for magnons.
#
# Copyright (C) 2023 Magnopy Team
#
# e-mail: anry@uv.es, web: magnopy.org
#
# This program is free software: you  can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the  Free Software
# Foundation,  either  version 3  of the License,  or (at your option) any later
# version.
#
# This program is distributed in the  hope  that it will be useful,  but WITHOUT
# ANY WARRANTY;  without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the  GNU General Public License  along with
# this program.  If not, see <https://www.gnu.org/licenses/>.
# ================================ END LICENSE =================================
import pytest

from magnopy._parameters._container import _npn

################################################################################
#                                     _npn                                     #
################################################################################
_NPN_REFERENCE = (
    ((), (0,), 1, 1),  # partition (1,)
    (((0, 0, 0),), (0, 0), 2, 1),  # partition (2, 0)
    (((0, 0, 0),), (0, 1), 2, 2),  # partition (1, 1)
    (((1, 0, 0),), (0, 0), 2, 2),  # partition (1, 1)
    (((0, 0, 0), (0, 0, 0)), (0, 0, 0), 3, 1),  # partition (3, 0, 0)
    (((0, 0, 0), (0, 0, 0)), (0, 1, 1), 3, 2),  # partition (2, 1, 0)
    (((1, 0, 0), (1, 0, 0)), (0, 0, 0), 3, 2),  # partition (2, 1, 0)
    (((0, 0, 0), (0, 0, 0)), (0, 0, 1), 3, 2),  # partition (2, 1, 0)
    (((0, 0, 0), (1, 0, 0)), (0, 0, 0), 3, 2),  # partition (2, 1, 0)
    (((0, 0, 0), (0, 0, 0)), (1, 0, 1), 3, 2),  # partition (2, 1, 0)
    (((-1, 0, 0), (0, 0, 0)), (0, 0, 0), 3, 2),  # partition (2, 1, 0)
    (((0, 0, 0), (0, 0, 0)), (0, 1, 2), 3, 3),  # partition (1, 1, 1)
    (((1, 0, 0), (0, 1, 0)), (0, 0, 0), 3, 3),  # partition (1, 1, 1)
    (((0, 0, 0), (0, 0, 0), (0, 0, 0)), (0, 0, 0, 0), 4, 1),  # partition (4, 0, 0, 0)
    (((0, 0, 0), (0, 0, 0), (0, 0, 0)), (0, 1, 1, 1), 4, 2),  # partition (3, 1, 0, 0)
    (((1, 0, 0), (1, 0, 0), (1, 0, 0)), (0, 0, 0, 0), 4, 2),  # partition (3, 1, 0, 0)
    (((0, 0, 0), (0, 0, 0), (0, 0, 0)), (0, 0, 1, 1), 4, 3),  # partition (2, 2, 0, 0)
    (((0, 0, 0), (1, 0, 0), (1, 0, 0)), (0, 0, 0, 0), 4, 3),  # partition (2, 2, 0, 0)
    (((0, 0, 0), (0, 0, 0), (0, 0, 0)), (1, 0, 1, 1), 4, 2),  # partition (3, 1, 0, 0)
    (((-1, 0, 0), (0, 0, 0), (0, 0, 0)), (0, 0, 0, 0), 4, 2),  # partition (3, 1, 0, 0)
    (((0, 0, 0), (0, 0, 0), (0, 0, 0)), (0, 1, 2, 2), 4, 4),  # partition (2, 1, 1, 0)
    (((1, 0, 0), (0, 1, 0), (0, 1, 0)), (0, 0, 0, 0), 4, 4),  # partition (2, 1, 1, 0)
    (((0, 0, 0), (0, 0, 0), (0, 0, 0)), (0, 0, 0, 1), 4, 2),  # partition (3, 1, 0, 0)
    (((0, 0, 0), (0, 0, 0), (1, 0, 0)), (0, 0, 0, 0), 4, 2),  # partition (3, 1, 0, 0)
    (((0, 0, 0), (0, 0, 0), (0, 0, 0)), (1, 0, 0, 1), 4, 3),  # partition (2, 2, 0, 0)
    (((-1, 0, 0), (-1, 0, 0), (0, 0, 0)), (0, 0, 0, 0), 4, 3),  # partition (2, 2, 0, 0)
    (((0, 0, 0), (0, 0, 0), (0, 0, 0)), (0, 1, 1, 2), 4, 4),  # partition (2, 1, 1, 0)
    (((1, 0, 0), (1, 0, 0), (0, 1, 0)), (0, 0, 0, 0), 4, 4),  # partition (2, 1, 1, 0)
    (((0, 0, 0), (0, 0, 0), (0, 0, 0)), (0, 1, 0, 1), 4, 3),  # partition (2, 2, 0, 0)
    (((1, 0, 0), (0, 0, 0), (1, 0, 0)), (0, 0, 0, 0), 4, 3),  # partition (2, 2, 0, 0)
    (((0, 0, 0), (0, 0, 0), (0, 0, 0)), (1, 1, 0, 1), 4, 2),  # partition (3, 1, 0, 0)
    (((0, 0, 0), (-1, 0, 0), (0, 0, 0)), (0, 0, 0, 0), 4, 2),  # partition (3, 1, 0, 0)
    (((0, 0, 0), (0, 0, 0), (0, 0, 0)), (0, 2, 1, 2), 4, 4),  # partition (2, 1, 1, 0)
    (((0, 1, 0), (1, 0, 0), (0, 1, 0)), (0, 0, 0, 0), 4, 4),  # partition (2, 1, 1, 0)
    (((0, 0, 0), (0, 0, 0), (0, 0, 0)), (0, 0, 1, 2), 4, 4),  # partition (2, 1, 1, 0)
    (((0, 0, 0), (1, 0, 0), (0, 1, 0)), (0, 0, 0, 0), 4, 4),  # partition (2, 1, 1, 0)
    (((0, 0, 0), (0, 0, 0), (0, 0, 0)), (1, 0, 1, 2), 4, 4),  # partition (2, 1, 1, 0)
    (((-1, 0, 0), (0, 0, 0), (-1, 1, 0)), (0, 0, 0, 0), 4, 4),  # partition (2, 1, 1, 0)
    (((0, 0, 0), (0, 0, 0), (0, 0, 0)), (2, 0, 1, 2), 4, 4),  # partition (2, 1, 1, 0)
    (((0, -1, 0), (1, -1, 0), (0, 0, 0)), (0, 0, 0, 0), 4, 4),  # partition (2, 1, 1, 0)
    (((0, 0, 0), (0, 0, 0), (0, 0, 0)), (0, 1, 2, 3), 4, 5),  # partition (1, 1, 1, 1)
    (((1, 0, 0), (0, 1, 0), (0, 0, 1)), (0, 0, 0, 0), 4, 5),  # partition (1, 1, 1, 1)
)


@pytest.mark.parametrize("nus, alphas, n, p_n", _NPN_REFERENCE)
def test_npn_by_reference(nus, alphas, n, p_n):
    assert _npn(nus, alphas) == (n, p_n)


@pytest.mark.parametrize(
    "nus, alphas",
    (
        ((), ()),
        (((0, 0, 0),) * 4, tuple(range(5))),
        (((0, 0, 0),) * 10, tuple(range(10))),
    ),
)
def test_npn_rejects_unsupported_rank(nus, alphas):
    with pytest.raises(ValueError, match="got n ="):
        _npn(nus, alphas)
