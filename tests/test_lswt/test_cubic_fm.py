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


from math import cos

import numpy as np
import pytest
from hypothesis import given
from hypothesis import strategies as st
from hypothesis.extra.numpy import arrays as harrays

from magnopy import LSWT, ColpaFailed, Convention, SpinHamiltonian


@pytest.mark.parametrize(
    "spin_direction",
    ([0, 0, 1],),  # , [0, 0, -1], [0, 1, 0], [1, 0, 0], [1, 1, 1]),
)
@given(
    harrays(
        np.float64,
        (3,),
        elements=st.floats(min_value=-1e8, max_value=1e8, allow_subnormal=False),
    )
)
def test_ferromagnet_one_spin_cubic(spin_direction, k):
    a = 1
    J = -1
    S = 0.5

    atoms = {"names": ["Fe"], "spins": [S], "positions": [[0, 0, 0]]}

    omega_anal = (
        2 * S * J * 6 * (1 / 3 * (cos(k[0] * a) + cos(k[1] * a) + cos(k[2] * a)) - 1)
    )

    spinham = SpinHamiltonian(
        cell=[[a, 0, 0], [0, a, 0], [0, 0, a]],
        atoms=atoms,
        convention=Convention(spin_normalized=False, multiple_counting=True, c22=1),
    )

    spinham.add_22(alpha=0, beta=0, nu=(1, 0, 0), parameter=J * np.eye(3, dtype=float))
    spinham.add_22(alpha=0, beta=0, nu=(0, 1, 0), parameter=J * np.eye(3, dtype=float))
    spinham.add_22(alpha=0, beta=0, nu=(0, 0, 1), parameter=J * np.eye(3, dtype=float))

    disp = LSWT(spinham=spinham, spin_directions=[spin_direction])

    assert disp.M == 1

    assert len(disp._J2) == 6
    assert len(disp._J1) == 1

    assert len(disp.A2) == 6
    assert len(disp.B2) == 6

    assert (
        abs(disp.A(k=k) - (3 - (cos(k[0] * a) + cos(k[1] * a) + cos(k[2] * a)))) < 1e-8
    )

    assert abs(disp.B(k=k)) < 1e-8

    assert np.allclose(
        disp.GDM(k=k),
        [
            [3 - (cos(k[0] * a) + cos(k[1] * a) + cos(k[2] * a)), 0],
            [0, 3 - (cos(k[0] * a) + cos(k[1] * a) + cos(k[2] * a))],
        ],
    )

    if cos(k[0] * a) - 1 == 0 and cos(k[1] * a) - 1 == 0 and cos(k[2] * a) - 1 == 0:
        try:
            omega = disp.omega(k=k)
        except ColpaFailed:
            return
    else:
        omega = disp.omega(k=k)

    assert len(omega) == 1

    omega = omega[0]

    assert omega.imag == 0

    assert (omega.real - omega_anal) < 1e-8
