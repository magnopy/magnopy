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


def _get_primary_p33(alpha, beta, gamma, nu, _lambda, parameter=None):
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
    gamma : int
        Index of the third atom.
    nu : tuple of 3 int
        Unit cell for the second atom.
    _lambda : tuple of 3 int
        Unit cell for the third atom.
    parameter : (3, 3, 3) :numpy:`ndarray`, optional
        Full matrix of the parameter.

    Returns
    -------
    alpha : int
        Index of the first atom.
    beta : int
        Index of the second atom.
    gamma : int
        Index of the third atom.
    nu : tuple of 3 int
        Unit cell for the second atom.
    _lambda : tuple of 3 int
        Unit cell for the third atom.
    parameter : (3, 3, 3) :numpy:`ndarray`, optional
        Full matrix of the parameter. It is returned only if ``parameter is not None``.
    """

    def _ordered(mu1, alpha1, mu2, alpha2, mu3, alpha3):
        return _spins_ordered(
            mu1=mu1, alpha1=alpha1, mu2=mu2, alpha2=alpha2
        ) and _spins_ordered(mu1=mu2, alpha1=alpha2, mu2=mu3, alpha2=alpha3)

    # Case 1
    if _ordered(
        mu1=(0, 0, 0), alpha1=alpha, mu2=nu, alpha2=beta, mu3=_lambda, alpha3=gamma
    ):
        pass
    # Case 2
    elif _ordered(
        mu1=(0, 0, 0), alpha1=alpha, mu2=_lambda, alpha2=gamma, mu3=nu, alpha3=beta
    ):
        alpha, beta, gamma = alpha, gamma, beta
        nu, _lambda = _lambda, nu
        if parameter is not None:
            parameter = np.transpose(parameter, (0, 2, 1))
    # Case 3
    elif _ordered(
        mu1=nu, alpha1=beta, mu2=(0, 0, 0), alpha2=alpha, mu3=_lambda, alpha3=gamma
    ):
        alpha, beta, gamma = beta, alpha, gamma
        nu1, nu2, nu3 = nu
        lambda1, lambda2, lambda3 = _lambda
        nu = (-nu1, -nu2, -nu3)
        _lambda = (lambda1 - nu1, lambda2 - nu2, lambda3 - nu3)
        if parameter is not None:
            parameter = np.transpose(parameter, (1, 0, 2))
    # Case 4
    elif _ordered(
        mu1=nu, alpha1=beta, mu2=_lambda, alpha2=gamma, mu3=(0, 0, 0), alpha3=alpha
    ):
        alpha, beta, gamma = beta, gamma, alpha
        nu1, nu2, nu3 = nu
        lambda1, lambda2, lambda3 = _lambda
        nu = (lambda1 - nu1, lambda2 - nu2, lambda3 - nu3)
        _lambda = (-nu1, -nu2, -nu3)
        if parameter is not None:
            parameter = np.transpose(parameter, (2, 0, 1))
    # Case 5
    elif _ordered(
        mu1=_lambda, alpha1=gamma, mu2=(0, 0, 0), alpha2=alpha, mu3=nu, alpha3=beta
    ):
        alpha, beta, gamma = gamma, alpha, beta
        nu1, nu2, nu3 = nu
        lambda1, lambda2, lambda3 = _lambda
        nu = (-lambda1, -lambda2, -lambda3)
        _lambda = (nu1 - lambda1, nu2 - lambda2, nu3 - lambda3)
        if parameter is not None:
            parameter = np.transpose(parameter, (1, 2, 0))
    # Case 6
    elif _ordered(
        mu1=_lambda, alpha1=gamma, mu2=nu, alpha2=beta, mu3=(0, 0, 0), alpha3=alpha
    ):
        alpha, beta, gamma = gamma, beta, alpha
        nu1, nu2, nu3 = nu
        lambda1, lambda2, lambda3 = _lambda
        nu = (nu1 - lambda1, nu2 - lambda2, nu3 - lambda3)
        _lambda = (-lambda1, -lambda2, -lambda3)
        if parameter is not None:
            parameter = np.transpose(parameter, (2, 1, 0))

    if parameter is None:
        return alpha, beta, gamma, nu, _lambda

    return alpha, beta, gamma, nu, _lambda, parameter


class _P33_iterator:
    R"""
    Iterator over the (three spins & three sites) parameters of the spin Hamiltonian.
    """

    def __init__(self, spinham) -> None:
        self.container = spinham._33
        self.mc = spinham.convention.multiple_counting
        self.length = len(self.container)
        self.index = 0

    def __next__(self):
        # Case 1
        if self.index < self.length:
            self.index += 1
            return self.container[self.index - 1]
        # Case 2
        elif self.mc and self.index < 2 * self.length:
            self.index += 1
            alpha, beta, gamma, nu, _lambda, parameter = self.container[
                self.index - 1 - self.length
            ]
            return [
                alpha,
                gamma,
                beta,
                _lambda,
                nu,
                np.transpose(parameter, (0, 2, 1)),
            ]
        # Case 3
        elif self.mc and self.index < 3 * self.length:
            self.index += 1
            (
                alpha,
                beta,
                gamma,
                (nu1, nu2, nu3),
                (lambda1, lambda2, lambda3),
                parameter,
            ) = self.container[self.index - 1 - 2 * self.length]
            return [
                beta,
                alpha,
                gamma,
                (-nu1, -nu2, -nu3),
                (lambda1 - nu1, lambda2 - nu2, lambda3 - nu3),
                np.transpose(parameter, (1, 0, 2)),
            ]
        # Case 4
        elif self.mc and self.index < 4 * self.length:
            self.index += 1
            (
                alpha,
                beta,
                gamma,
                (nu1, nu2, nu3),
                (lambda1, lambda2, lambda3),
                parameter,
            ) = self.container[self.index - 1 - 3 * self.length]
            return [
                beta,
                gamma,
                alpha,
                (lambda1 - nu1, lambda2 - nu2, lambda3 - nu3),
                (-nu1, -nu2, -nu3),
                np.transpose(parameter, (2, 0, 1)),
            ]
        # Case 5
        elif self.mc and self.index < 5 * self.length:
            self.index += 1
            (
                alpha,
                beta,
                gamma,
                (nu1, nu2, nu3),
                (lambda1, lambda2, lambda3),
                parameter,
            ) = self.container[self.index - 1 - 4 * self.length]
            return [
                gamma,
                alpha,
                beta,
                (-lambda1, -lambda2, -lambda3),
                (nu1 - lambda1, nu2 - lambda2, nu3 - lambda3),
                np.transpose(parameter, (1, 2, 0)),
            ]
        # Case 6
        elif self.mc and self.index < 6 * self.length:
            self.index += 1
            (
                alpha,
                beta,
                gamma,
                (nu1, nu2, nu3),
                (lambda1, lambda2, lambda3),
                parameter,
            ) = self.container[self.index - 1 - 5 * self.length]
            return [
                gamma,
                beta,
                alpha,
                (nu1 - lambda1, nu2 - lambda2, nu3 - lambda3),
                (-lambda1, -lambda2, -lambda3),
                np.transpose(parameter, (2, 1, 0)),
            ]

        raise StopIteration

    def __len__(self):
        return self.length * (1 + 5 * int(self.mc))

    def __iter__(self):
        return self


@property
def _p33(spinham):
    r"""
    Parameters of (three spins & three sites) term of the Hamiltonian.

    .. math::

        \boldsymbol{J}_{3,3}(\boldsymbol{r}_{\nu,\alpha\beta}, \boldsymbol{r}_{\lambda,\alpha\gamma})

    of the term

    .. math::

        C_{3,3}
        \sum_{\substack{\mu, \nu, \alpha, \beta,\\ i, j, u}}
        J^{iju}_{3,3}(\boldsymbol{r}_{\nu,\alpha\beta}, \boldsymbol{r}_{\lambda,\alpha\gamma})
        S_{\mu,\alpha}^i
        S_{\mu+\nu,\beta}^j
        S_{\mu+\lambda, \gamma}^u

    Returns
    -------
    parameters : iterator
        List of parameters. The list has a form of

        .. code-block:: python

            [[alpha, beta, gamma, nu, lambda, J], ...]

        where

        ``alpha`` is an index of the atom located in the (0,0,0) unit cell.

        ``beta`` is an index of the atom located in the  nu unit cell.

        ``gamma`` is an index of the atom located in the  lambda unit cell.

        ``nu`` defines the unit cell of the second atom (beta). It is a tuple of 3
        integers.

        ``lambda`` defines the unit cell of the third atom (gamma). It is a tuple of 3
        integers.

        ``J`` is a (3, 3, 3) :numpy:`ndarray`.

    See Also
    --------
    add_33
    remove_33
    """

    return _P33_iterator(spinham)


def _add_33(
    spinham,
    alpha: int,
    beta: int,
    gamma: int,
    nu: tuple,
    _lambda: tuple,
    parameter,
    replace=False,
) -> None:
    r"""
    Adds a (three spins & three sites) parameter to the Hamiltonian.

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
    gamma : int
        Index of an atom from the _lambda unit cell.

        ``0 <= gamma < len(spinham.atoms.names)``.
    nu : tuple of 3 int
        Three relative coordinates with respect to the three lattice vectors, that
        specify the unit cell for the second atom.

        .. math::

            \nu
            =
            (x_{\boldsymbol{a}_1}, x_{\boldsymbol{a}_2}, x_{\boldsymbol{a}_3})
    _lambda : tuple of 3 int
        Three relative coordinates with respect to the three lattice vectors, that
        specify the unit cell for the third atom.

        .. math::

            \lambda
            =
            (x_{\boldsymbol{a}_1}, x_{\boldsymbol{a}_2}, x_{\boldsymbol{a}_3})

    parameter : (3, 3, 3) |array-like|_
        Value of the parameter (:math:`3\times3\times3` matrix).
    replace : bool, default False
        Whether to replace the value of the parameter if the triplet of atoms
        ``alpha, beta, gamma, nu, lambda`` or one of its duplicates already have a
        parameter associated with it.

    See Also
    --------
    p33
    remove_33

    Notes
    -----
    If ``spinham.convention.multiple_counting`` is ``True``, then this function adds
    the bond and all its duplicates to the Hamiltonian. It will cause an ``ValueError``
    to add the duplicate of the bond after the bond is added.

    If ``spinham.convention.multiple_counting`` is ``False``, then only the primary
    version of the bond is added to the Hamiltonian.

    For the definition of the primary version see
    :ref:`user-guide_theory-behind_multiple-counting`.
    """

    _validate_atom_index(index=alpha, atoms=spinham.atoms)
    _validate_atom_index(index=beta, atoms=spinham.atoms)
    _validate_atom_index(index=gamma, atoms=spinham.atoms)
    _validate_unit_cell_index(ijk=nu)
    _validate_unit_cell_index(ijk=_lambda)
    spinham._reset_internals()

    parameter = np.array(parameter)

    alpha, beta, gamma, nu, _lambda, parameter = _get_primary_p33(
        alpha=alpha, beta=beta, gamma=gamma, nu=nu, _lambda=_lambda, parameter=parameter
    )

    # TD-BINARY_SEARCH

    # Try to find the place for the new one inside the list
    index = 0
    while index < len(spinham._33):
        # If already present in the model
        if spinham._33[index][:5] == [alpha, beta, gamma, nu, _lambda]:
            # Either replace
            if replace:
                spinham._33[index] = [alpha, beta, gamma, nu, _lambda, parameter]
                return
            # Or raise an error
            raise ValueError(
                f"Parameter is already set for the triple of atoms "
                f"{alpha}, {beta} {nu}, {gamma} {_lambda}. Or for their duplicate."
            )

        # If it should be inserted before current element
        if spinham._33[index][:5] > [alpha, beta, gamma, nu, _lambda]:
            spinham._33.insert(index, [alpha, beta, gamma, nu, _lambda, parameter])
            return

        index += 1

    # If it should be inserted at the end or at the beginning of the list
    spinham._33.append([alpha, beta, gamma, nu, _lambda, parameter])


def _remove_33(
    spinham, alpha: int, beta: int, gamma: int, nu: tuple, _lambda: tuple
) -> None:
    r"""
    Removes a (three spins & three sites) parameter from the Hamiltonian.

    Duplicates of the bonds are managed automatically (independently of the convention of
    the Hamiltonian).

    Parameters
    ----------
    alpha : int
        Index of an atom from the (0, 0, 0) unit cell.

        ``0 <= alpha < len(spinham.atoms.names)``.
    beta : int
        Index of an atom from the nu unit cell.

        ``0 <= beta < len(spinham.atoms.names)``.
    gamma : int
        Index of an atom from the _lambda unit cell.

        ``0 <= gamma < len(spinham.atoms.names)``.
    nu : tuple of 3 int
        Three relative coordinates with respect to the three lattice vectors, that
        specify the unit cell for the second atom.

        .. math::

            \nu
            =
            (x_{\boldsymbol{a}_1}, x_{\boldsymbol{a}_2}, x_{\boldsymbol{a}_3})
    _lambda : tuple of 3 int
        Three relative coordinates with respect to the three lattice vectors, that
        specify the unit cell for the third atom.

        .. math::

            \lambda
            =
            (x_{\boldsymbol{a}_1}, x_{\boldsymbol{a}_2}, x_{\boldsymbol{a}_3})

    See Also
    --------
    p33
    add_33

    Notes
    -----
    If ``spinham.convention.multiple_counting`` is ``True``, then this function removes
    all versions of the bond from the Hamiltonian.

    If ``spinham.convention.multiple_counting`` is ``False``, then this function removes
    the primary version of the given bond.

    For the definition of the primary version see
    :ref:`user-guide_theory-behind_multiple-counting`.
    """

    _validate_atom_index(index=alpha, atoms=spinham.atoms)
    _validate_atom_index(index=beta, atoms=spinham.atoms)
    _validate_atom_index(index=gamma, atoms=spinham.atoms)
    _validate_unit_cell_index(ijk=nu)
    _validate_unit_cell_index(ijk=_lambda)

    alpha, beta, gamma, nu, _lambda = _get_primary_p33(
        alpha=alpha, beta=beta, gamma=gamma, nu=nu, _lambda=_lambda
    )

    # TD-BINARY_SEARCH

    for index in range(len(spinham._33)):
        # As the list is sorted, there is no point in resuming the search
        # when a larger element is found
        if spinham._33[index][:5] > [alpha, beta, gamma, nu, _lambda]:
            return

        if spinham._33[index][:5] == [alpha, beta, gamma, nu, _lambda]:
            del spinham._33[index]
            spinham._reset_internals()
            return
