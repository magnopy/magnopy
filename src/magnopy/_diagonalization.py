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


from functools import cmp_to_key

import numpy as np
from numpy.linalg import LinAlgError

from magnopy._exceptions import ColpaFailed

# Save local scope at this moment
old_dir = set(dir())
old_dir.add("old_dir")


def _check_grand_dynamical_matrix(D):
    r"""
    Check that grand dynamical matrix is (2N, 2N) matrix

    Parameters
    ----------
    D : |array-like|_
        Candidate for the grand dynamical matrix

    Returns
    -------
    D : (2N, 2N) :numpy:`ndarray`
        Grand dynamical matrix.
    N : int

    Raises
    ------
    ValueError
        If the check is not passed
    """

    D = np.array(D)

    if len(D.shape) != 2:
        raise ValueError(f"Grand dynamical matrix is not 2-dimensional, got {D.shape}.")

    if D.shape[0] != D.shape[1]:
        raise (f"Grand dynamical matrix is not square, got {D.shape}.")

    if D.shape[0] % 2 != 0:
        raise (f"Size of the grand dynamical matrix is not even, got {D.shape}.")

    return D, D.shape[0] // 2


def _inverse_by_colpa(matrix):
    # Compute G from G^-1 (or vise versa) following Colpa, see equation (3.7) for details

    N = matrix.shape[0] // 2
    matrix = np.conjugate(matrix).T
    matrix[:N, N:] *= -1
    matrix[N:, :N] *= -1

    return matrix


def solve_via_colpa(D, sort_by_first_N=True):
    r"""
    Diagonalize grand-dynamical matrix following the method of Colpa (section 3, remark
    1 of [1]_).

    Solves the Bogoliubov Hamiltonian of the form

    .. math::

        \hat{H} = \sum_{r^{\prime}, r = 1}^m
        a_{r^{\prime}}^{\dagger}\boldsymbol{\Delta}_1^{r^{\prime}r}a_r +
        a_{r^{\prime}}^{\dagger}\boldsymbol{\Delta}_2^{r^{\prime}r}a_{m+r}^{\dagger} +
        a_{m+r^{\prime}}^{\dagger}\boldsymbol{\Delta}_3^{r^{\prime}r}a_r +
        a_{m+r^{\prime}}^{\dagger}\boldsymbol{\Delta}_4^{r^{\prime}r}a_{m+r}^{\dagger}

    ensuring the bosonic commutation relations.

    In a matrix form the Hamiltonian can be written as

    .. math::

        \hat{H} = \boldsymbol{\cal A}^{\dagger} \boldsymbol{D} \boldsymbol{\cal A}

    where

    .. math::

        \boldsymbol{\cal A} =
        \begin{pmatrix}
            a_1 \\
            \cdots \\
            a_m \\
            a_{m+1} \\
            \cdots \\
            a_{2m} \\
        \end{pmatrix}

    After diagonalization the Hamiltonian has the form

    .. math::

        \hat{H}
        =
        \boldsymbol{\cal B}^{\dagger} \boldsymbol{\mathcal{E}} \boldsymbol{\cal B}

    where

    .. math::

        \boldsymbol{\cal B} =
        \begin{pmatrix}
            b_1 \\
            \cdots \\
            b_m \\
            b_{m+1} \\
            \cdots \\
            b_{2m} \\
        \end{pmatrix}

    Parameters
    ----------
    D : (2N, 2N) |array-like|_
        Grand dynamical matrix. If it is Hermitian and positive-defined, then obtained
        eigenvalues are positive and real.

        .. math::

            \boldsymbol{\mathcal{D}} = \begin{pmatrix}
                \boldsymbol{\Delta_1} & \boldsymbol{\Delta_2} \\
                \boldsymbol{\Delta_3} & \boldsymbol{\Delta_4}
            \end{pmatrix}

    Returns
    -------
    E : (2N,) :numpy:`ndarray`
        The eigenvalues. It is an array of the diagonal elements of the diagonal matrix
        :math:`\boldsymbol{\mathcal{E}}`. First N elements correspond to the
        :math:`b^{\dagger}_1b_1, \dots, b^{\dagger}_mb_m` and last N elements - to
        the :math:`b^{\dagger}_{m+1}b_{m+1}, \dots, b^{\dagger}_{2m}b_{2m}`.

        Eigenvalues are sorted individually for the first N and the last N elements,
        based on the transformation matrix and not on the values of E itself.

        .. math::

            \boldsymbol{\mathcal{E}}
            =
            (\boldsymbol{G}^{\dagger})^{-1} \boldsymbol{D} \boldsymbol{G}^{-1}

    G : (2N, 2N) :numpy:`ndarray`
        Transformation matrix, that changes the basis from the original set of bosonic
        operators :math:`\boldsymbol{a}` to the set of new bosonic operators
        :math:`\boldsymbol{b}` which diagonalize the Hamiltonian:

        .. math::
            \boldsymbol{\cal B} = \boldsymbol{G} \boldsymbol{\cal A}

        The rows are ordered in the same way as the eigenvalues.

    Raises
    ------
    ColpaFailed
        If the algorithm fails. Typically it means that the grand dynamical matrix
        :math:`\boldsymbol{D}` is not positive-defined.
    ValueError
        If the grand dynamical matrix is not square or its shape is not even.

    References
    ----------
    .. [1] Colpa, J.H.P., 1978.
        Diagonalization of the quadratic boson hamiltonian.
        Physica A: Statistical Mechanics and its Applications,
        93(3-4), pp.327-353.

    Examples
    --------

    For already diagonal matrix this function does not do much

    .. doctest::

        >>> import magnopy
        >>> D = [[1, 0], [0, 2]]
        >>> E, G = magnopy.solve_via_colpa(D)
        >>> E
        array([1., 2.])
        >>> G
        array([[ 1., -0.],
               [-0.,  1.]])

    .. doctest::

        >>> import magnopy
        >>> D = [[1, 1j], [-1j, 2]]
        >>> E, G = magnopy.solve_via_colpa(D)
        >>> E
        array([0.61803399+0.j, 1.61803399+0.j])
        >>> G
        array([[ 1.08204454-0.j        ,  0.        +0.41330424j],
               [-0.        -0.41330424j,  1.08204454-0.j        ]])
        >>> E, G = magnopy.solve_via_colpa(D)
        >>> E
        array([0.61803399+0.j, 1.61803399+0.j])
        >>> G # doctest: +SKIP
        array([[1.08204454+0.j        , 0.        -0.41330424j],
               [0.        +0.41330424j, 1.08204454+0.j        ]])
    """

    D, N = _check_grand_dynamical_matrix(D)

    g = np.diag(np.concatenate((np.ones(N), -np.ones(N))))

    try:
        # In Colpa article decomposition is K^{\dag}K, while numpy gives KK^{\dag}
        K = np.conjugate(np.linalg.cholesky(D)).T
    except LinAlgError:
        raise ColpaFailed

    L, U = np.linalg.eig(K @ g @ np.conjugate(K).T)

    # Sort with respect to L, in descending order
    U = np.concatenate((L[:, None], U.T), axis=1).T
    # U = np.concatenate((L[np.newaxis, :], U), axis=0)
    U = U[:, np.argsort(U[0])]
    L = U[0, ::-1]
    U = U[1:, ::-1]

    E = g @ L

    G_inv = np.linalg.inv(K) @ U @ np.sqrt(np.diag(E))

    G = _inverse_by_colpa(G_inv)

    # Sort first N and second N individually based on the transformation matrix
    tmp = np.concatenate((E[:, np.newaxis], G), axis=1)

    def compare(array1, array2):
        if sort_by_first_N:
            difference = np.round(
                array1[1 : N + 1].real - array2[1 : N + 1].real, decimals=15
            )
        else:
            difference = np.round(
                array1[N + 1 :].real - array2[N + 1 :].real, decimals=15
            )

        if np.allclose(difference, np.zeros(difference.shape)):
            return 0.0

        return difference[np.nonzero(difference)[0][0]]

    upper_part = np.array(sorted(tmp[:N], key=cmp_to_key(compare)))
    lower_part = np.array(sorted(tmp[N:], key=cmp_to_key(compare)))

    E = np.concatenate((upper_part[:, 0], lower_part[:, 0]))
    G = np.concatenate((upper_part[:, 1:], lower_part[:, 1:]), axis=0)

    return E, G


# Populate __all__ with objects defined in this file
__all__ = list(set(dir()) - old_dir)
# Remove all semi-private objects
__all__ = [i for i in __all__ if not i.startswith("_")]
del old_dir
