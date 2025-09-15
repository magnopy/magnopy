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

from magnopy._spinham._validators import _validate_atom_index


@property
def _p1(spinham) -> list:
    r"""
    Parameters of (one spin & one site) term of the Hamiltonian.

    .. math::

        \boldsymbol{J}_1(\boldsymbol{r}_{\alpha})

    of the term

    .. math::

        C_1
        \sum_{\mu, \alpha, i}
        J_1^i(\boldsymbol{r}_{\alpha})
        S_{\mu,\alpha}^i

    Returns
    -------
    parameters : list
        List of parameters. The list has a form of

        .. code-block:: python

            [[alpha, J], ...]

        ``0 <= len(parameters) <= len(spinham.atoms.names)``.

        where ``alpha`` is an index of the atom to which the parameter is assigned and
        ``J``  is a (3, ) :numpy:`ndarray`. The parameters are sorted by the index of an
        atom ``alpha`` in the ascending order.

    See Also
    --------
    add_1
    remove_1
    """

    return spinham._1


def _add_1(spinham, alpha: int, parameter, replace=False) -> None:
    r"""
    Adds a (one spin & one site) parameter to the Hamiltonian.

    Parameters
    ----------
    alpha : int
        Index of an atom, with which the parameter is associated.

        ``0 <= alpha < len(spinham.atoms.names)``.
    parameter : (3, ) |array-like|_
        Value of the parameter (:math:`3\times1` vector).
    replace : bool, default False
        Whether to replace the value of the parameter if an atom already has a
        parameter associated with it.

    Raises
    ------
    ValueError
        If an atom already has a parameter associated with it.

    See Also
    --------
    p1
    remove_1
    """

    _validate_atom_index(index=alpha, atoms=spinham.atoms)
    spinham._reset_internals()

    parameter = np.array(parameter)

    # TD-BINARY_SEARCH
    # Try to find the place for the new one inside the list
    index = 0
    while index < len(spinham._1):
        # If already present in the model
        if spinham._1[index][0] == alpha:
            # Either replace
            if replace:
                spinham._1[index] = [alpha, parameter]
                return
            # Or raise an error
            raise ValueError(
                f"(One spin & one site) parameter is already set "
                f"for atom {alpha} ('{spinham.atoms.names[alpha]}')"
            )

        # If it should be inserted before current element
        if spinham._1[index][0] > alpha:
            spinham._1.insert(index, [alpha, parameter])
            return

        index += 1

    # If it should be inserted at the end or at the beginning of the list
    spinham._1.append([alpha, parameter])


def _remove_1(spinham, alpha: int) -> None:
    r"""
    Removes a (one spin & one site) parameter from the Hamiltonian.

    Parameters
    ----------
    alpha : int
        Index of an atom, with which the parameter is associated.

        ``0 <= alpha < len(spinham.atoms.names)``.

    See Also
    --------
    p1
    add_1
    """

    _validate_atom_index(index=alpha, atoms=spinham.atoms)

    for i in range(len(spinham._1)):
        # As the list is sorted, there is no point in resuming the search
        # when a larger element is found
        if spinham._1[i][0] > alpha:
            return

        if spinham._1[i][0] == alpha:
            del spinham._1[i]
            spinham._reset_internals()
            return
