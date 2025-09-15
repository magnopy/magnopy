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

from magnopy._spinham._validators import (
    _spins_ordered,
    _validate_atom_index,
    _validate_unit_cell_index,
)


def _get_primary_p422(alpha, beta, nu, parameter=None):
    r"""
    Return the primary version of the parameter.

    For the definition of the primary version see
    :ref:`user-guide_theory-behind_multiple-counting`.

    Parameters
    ----------
    alpha : int
        Index of the first atom.
    beta : int
        Index of the second atom.
    nu : tuple of 3 int
        Unit cell for the second atom.
    parameter : (3, 3, 3, 3) :numpy:`ndarray`, optional
        Full matrix of the parameter.

    Returns
    -------
    alpha : int
        Index of the first atom.
    beta : int
        Index of the second atom.
    nu : tuple of 3 int
        Unit cell for the second atom.
    parameter : (3, 3, 3, 3) :numpy:`ndarray`
        Full matrix of the parameter. It is returned only if ``parameter is not None``.
    """

    if _spins_ordered(mu1=(0, 0, 0), alpha1=alpha, mu2=nu, alpha2=beta):
        pass
    else:
        i, j, k = nu
        alpha, beta, nu = beta, alpha, (-i, -j, -k)
        if parameter is not None:
            parameter = np.transpose(parameter, (2, 3, 0, 1))

    if parameter is None:
        return alpha, beta, nu
    return alpha, beta, nu, parameter


class _P422_iterator:
    R"""
    Iterator over the (four spins & two sites (2+2)) parameters of the spin Hamiltonian.
    """

    def __init__(self, spinham) -> None:
        self.container = spinham._422
        self.mc = spinham.convention.multiple_counting
        self.length = len(self.container)
        self.index = 0

    def __next__(self):
        if self.index < self.length:
            self.index += 1
            return self.container[self.index - 1]

        elif self.mc and self.index < 2 * self.length:
            self.index += 1
            alpha, beta, (i, j, k), parameter = self.container[
                self.index - 1 - self.length
            ]
            return [beta, alpha, (-i, -j, -k), np.transpose(parameter, (2, 3, 0, 1))]

        raise StopIteration

    def __len__(self):
        return self.length * (1 + int(self.mc))

    def __iter__(self):
        return self


@property
def _p422(spinham):
    r"""
    Parameters of (four spins & two sites (2+2)) term of the Hamiltonian.

    .. math::

        \boldsymbol{J}_{4,2,2}(\boldsymbol{r}_{\nu,\alpha\beta})

    of the term

    .. math::

        C_{4,2,2}
        \sum_{\substack{\mu, \nu, \alpha, \beta,\\ i, j, u, v}}
        J^{ijuv}_{4,2,2}(\boldsymbol{r}_{\nu,\alpha\beta})
        S_{\mu,\alpha}^i
        S_{\mu,\alpha}^j
        S_{\mu+\nu,\beta}^u
        S_{\mu+\nu,\beta}^v

    Returns
    -------
    parameters : iterator
        List of parameters. The list has a form of

        .. code-block:: python

            [[alpha, beta, nu, J], ...]

        where

        ``alpha`` is an index of the atom located in the (0,0,0) unit cell.

        ``beta`` is an index of the atom located in the  nu unit cell.

        ``nu`` defines the unit cell of the second atom (beta). It is a tuple of 3
        integers.

        ``J`` is a (3, 3, 3, 3) :numpy:`ndarray`.

    See Also
    --------
    add_422
    remove_422
    """

    return _P422_iterator(spinham)


def _add_422(
    spinham, alpha: int, beta: int, nu: tuple, parameter, replace=False
) -> None:
    r"""
    Adds a (four spins & two sites (2+2)) parameter to the Hamiltonian.

    Doubles of the bonds are managed automatically (independently of the convention of the
    Hamiltonian).


    Raises
    ------
    ValueError
        If an atom already has a parameter associated with it.

    Parameters
    ----------
    alpha : int
        Index of an atom from the (0, 0, 0) unit cell.

        ``0 <= alpha < len(spinham.atoms.names)``.
    beta : int
        Index of an atom from the nu unit cell.

        ``0 <= beta < len(spinham.atoms.names)``.
    nu : tuple of 3 int
        Three relative coordinates with respect to the three lattice vectors, that
        specify the unit cell for the second atom.

        .. math::

            \nu
            =
            (x_{\boldsymbol{a}_1}, x_{\boldsymbol{a}_2}, x_{\boldsymbol{a}_3})

    parameter : (3, 3, 3, 3) |array-like|_
        Value of the parameter (:math:`3\times3\times3\times3` matrix).
    replace : bool, default False
        Whether to replace the value of the parameter if the pair of atoms
        ``alpha, beta, nu`` or its double already have a parameter associated
        with it.

    See Also
    --------
    p422
    remove_422

    Notes
    -----
    If ``spinham.convention.multiple_counting`` is ``True``, then this function adds both
    the bond and its double to the Hamiltonian. It will cause an ``ValueError`` to
    add the double of the bond after the bond is added.

    If ``spinham.convention.multiple_counting`` is ``False``, then only the primary
    version of the bond is added to the Hamiltonian.

    For the definition of the primary version see
    :ref:`user-guide_theory-behind_multiple-counting`.
    """

    _validate_atom_index(index=alpha, atoms=spinham.atoms)
    _validate_atom_index(index=beta, atoms=spinham.atoms)
    _validate_unit_cell_index(ijk=nu)
    spinham._reset_internals()

    parameter = np.array(parameter)

    alpha, beta, nu, parameter = _get_primary_p422(
        alpha=alpha, beta=beta, nu=nu, parameter=parameter
    )

    # TD-BINARY_SEARCH

    # Try to find the place for the new one inside the list
    index = 0
    while index < len(spinham._422):
        # If already present in the model
        if spinham._422[index][:3] == [alpha, beta, nu]:
            # Either replace
            if replace:
                spinham._422[index] = [alpha, beta, nu, parameter]
                return
            # Or raise an error
            raise ValueError(
                f"Exchange like parameter is already set for the pair of atoms "
                f"{alpha} and {beta} ({nu}). Or for their double bond."
            )

        # If it should be inserted before current element
        if spinham._422[index][:3] > [alpha, beta, nu]:
            spinham._422.insert(index, [alpha, beta, nu, parameter])
            return

        index += 1

    # If it should be inserted at the end or at the beginning of the list
    spinham._422.append([alpha, beta, nu, parameter])


def _remove_422(spinham, alpha: int, beta: int, nu: tuple) -> None:
    r"""
    Removes a (four spins & two sites (2+2)) parameter from the Hamiltonian.

    Doubles of the bonds are managed automatically (independently of the convention of the
    Hamiltonian).

    Parameters
    ----------
    alpha : int
        Index of an atom from the (0, 0, 0) unit cell.

        ``0 <= alpha < len(spinham.atoms.names)``.
    beta : int
        Index of an atom from the nu unit cell.

        ``0 <= beta < len(spinham.atoms.names)``.
    nu : tuple of 3 int
        Three relative coordinates with respect to the three lattice vectors, that
        specify the unit cell for the second atom.

        .. math::

            \nu
            =
            (x_{\boldsymbol{a}_1}, x_{\boldsymbol{a}_2}, x_{\boldsymbol{a}_3})

    See Also
    --------
    p422
    add_422

    Notes
    -----
    If ``spinham.convention.multiple_counting`` is ``True``, then this function removes
    all versions of the bond from the Hamiltonian.

    If ``spinham.convention.multiple_counting`` is ``False``, then this function removes
    the primary version of the given bond.

    For the definition of the primary version see
    :ref:`user-guide_theory-behind_multiple-counting`.

    For instance, if ``(1, 0, (0, 0, 0))`` is given, then this function attempts to
    remove either both ``(1, 0, (0, 0, 0))`` and ``(0, 1, (0, 0, 0))`` if
    ``spinham.convention.multiple_counting == True`` or the primary version
    ``(0, 1, (0, 0, 0))`` if ``spinham.convention.multiple_counting == False``.
    """

    _validate_atom_index(index=alpha, atoms=spinham.atoms)
    _validate_atom_index(index=beta, atoms=spinham.atoms)
    _validate_unit_cell_index(ijk=nu)

    alpha, beta, nu = _get_primary_p422(alpha=alpha, beta=beta, nu=nu)

    # TD-BINARY_SEARCH

    for index in range(len(spinham._422)):
        # As the list is sorted, there is no point in resuming the search
        # when a larger element is found
        if spinham._422[index][:3] > [alpha, beta, nu]:
            return

        if spinham._422[index][:3] == [alpha, beta, nu]:
            del spinham._422[index]
            spinham._reset_internals()
            return
