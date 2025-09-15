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

# Save local scope at this moment
old_dir = set(dir())
old_dir.add("old_dir")


def to_iso(parameter, matrix_form=False):
    r"""
    Extracts isotropic part of the full matrix parameter.

    .. math::

        C_{2,2}
        \boldsymbol{S}_{\mu}
        \boldsymbol{J}
        \boldsymbol{S}_{\nu}
        =
        C_{2,2}
        J_{iso}
        \boldsymbol{S}_{\mu}
        \cdot
        \boldsymbol{S}_{\nu}
        +
        C_{2,2}
        \boldsymbol{S}_{\mu}
        \boldsymbol{J}_S
        \boldsymbol{S}_{\nu}
        +
        C_{2,2}
        \boldsymbol{D}
        \cdot
        \left(
        \boldsymbol{S}_{\mu}
        \times
        \boldsymbol{S}_{\nu}
        \right)

    where :math:`J_{iso}` is defined as

    .. math::
        J_{iso} = \dfrac{tr(\boldsymbol{J})}{3}

    This term can be written in the matrix form as

    .. math::

        C_{2,2}
        J_{iso}
        \boldsymbol{S}_{\mu}
        \cdot
        \boldsymbol{S}_{\nu}
        &=
        C_{2,2}
        \boldsymbol{S}_{\mu}
        \boldsymbol{J}_{iso}
        \boldsymbol{S}_{\nu}
        \\
        \boldsymbol{J}_{iso}
        &=
        \begin{pmatrix}
            J_{iso} & 0 & 0 \\
            0 & J_{iso} & 0 \\
            0 & 0 & J_{iso} \\
        \end{pmatrix}

    Parameters
    ----------
    parameter : (3, 3) |array-like|_
        Full matrix of the exchange parameter (:math:`\boldsymbol{J}`).
    matrix_form : bool, default False
        Whether to return isotropic part of the matrix instead of isotropic parameter.

    Returns
    -------
    iso : float or (3, 3) :numpy:`ndarray`
        Isotropic parameter.

        * If ``matrix_form == False``, then returns a number :math:`J_{iso}`.
        * If ``matrix_form == True``, then returns a matrix :math:`\boldsymbol{J}_{iso}`.

    See Also
    --------
    from_iso
    to_dmi
    to_symm_anisotropy

    Examples
    --------

    .. doctest::

        >>> from magnopy import converter22
        >>> matrix = [[1, 3, 4], [-1, -2, 3], [4, 0, 10]]
        >>> converter22.to_iso(matrix)
        3.0
        >>> converter22.to_iso(matrix, matrix_form=True)
        array([[3., 0., 0.],
               [0., 3., 0.],
               [0., 0., 3.]])
    """

    iso = (parameter[0][0] + parameter[1][1] + parameter[2][2]) / 3

    if matrix_form:
        return iso * np.eye(3, dtype=float)

    return float(iso)


def to_symm_anisotropy(parameter):
    r"""
    Extracts traceless, symmetric anisotropic part of the full matrix parameter.

    .. math::

        C_{2,2}
        \boldsymbol{S}_{\mu}
        \boldsymbol{J}
        \boldsymbol{S}_{\nu}
        =
        C_{2,2}
        J_{iso}
        \boldsymbol{S}_{\mu}
        \cdot
        \boldsymbol{S}_{\nu}
        +
        C_{2,2}
        \boldsymbol{S}_{\mu}
        \boldsymbol{J}_S
        \boldsymbol{S}_{\nu}
        +
        C_{2,2}
        \boldsymbol{D}
        \cdot
        \left(
        \boldsymbol{S}_{\mu}
        \times
        \boldsymbol{S}_{\nu}
        \right)

    where matrix :math:`\boldsymbol{J}_S` is defined as

    .. math::

        \boldsymbol{J}_S
        =
        \dfrac{\boldsymbol{J}
        +
        \boldsymbol{J}^T}{2}

    Parameters
    ----------
    parameter : (3, 3) |array-like|_
        Full matrix of the exchange parameter (:math:`\boldsymbol{J}`).

    Returns
    -------
    aniso : float or (3, 3) :numpy:`ndarray`
        Matrix of a traceless, symmetric anisotropy.

    See Also
    --------
    to_iso
    to_dmi

    Examples
    --------

    .. doctest::

        >>> import magnopy
        >>> matrix = [[1, 3, 4], [-1, -2, 0], [4, 0, 10]]
        >>> magnopy.converter22.to_symm_anisotropy(matrix)
        array([[-2.,  1.,  4.],
               [ 1., -5.,  0.],
               [ 4.,  0.,  7.]])
    """

    parameter = np.array(parameter)

    return (parameter + parameter.T) / 2 - to_iso(parameter=parameter, matrix_form=True)


def to_dmi(parameter, matrix_form=False):
    r"""
    Extracts antisymmetric part (DMI) of the full matrix parameter.


    .. math::

        C_{2,2}
        \boldsymbol{S}_{\mu}
        \boldsymbol{J}
        \boldsymbol{S}_{\nu}
        =
        C_{2,2}
        J_{iso}
        \boldsymbol{S}_{\mu}
        \cdot
        \boldsymbol{S}_{\nu}
        +
        C_{2,2}
        \boldsymbol{S}_{\mu}
        \boldsymbol{J}_S
        \boldsymbol{S}_{\nu}
        +
        C_{2,2}
        \boldsymbol{D}
        \cdot
        \left(
        \boldsymbol{S}_{\mu}
        \times
        \boldsymbol{S}_{\nu}
        \right)

    where vector :math:`\boldsymbol{D}` is defined as

    .. math::

        \begin{pmatrix}
            0 & D^z & -D^y \\
            -D^z & 0 & D^x \\
            D^y & -D^x & 0 \\
        \end{pmatrix}
        =
        \dfrac{\boldsymbol{J}
        -
        \boldsymbol{J}^T}{2}
        =
        \boldsymbol{J}_A


    Parameters
    ----------
    parameter : (3, 3) |array-like|_
        Full matrix of the exchange parameter (:math:`\boldsymbol{J}`).
    matrix_form : bool, default False
        Whether to return dmi in a matrix form instead of a vector.

    Returns
    -------
    dmi : (3,) or (3, 3) :numpy:`ndarray`
        Antisymmetric exchange (DMI).

        * If ``matrix_form == False``, then returns a vector :math:`\boldsymbol{D}`.
        * If ``matrix_form == True``, then returns a matrix :math:`\boldsymbol{J}_A`.

    See Also
    --------
    from_dmi
    to_iso
    to_symm_anisotropy

    Examples
    --------

    .. doctest::

        >>> from magnopy import converter22
        >>> parameter = [[1, 3, 0], [-1, -2, 3], [0, 3, 9]]
        >>> converter22.to_dmi(parameter)
        array([0., 0., 2.])
        >>> converter22.to_dmi(parameter, matrix_form = True)
        array([[ 0.,  2.,  0.],
               [-2.,  0.,  0.],
               [ 0.,  0.,  0.]])
    """

    parameter = np.array(parameter)

    asymm_matrix = (parameter - parameter.T) / 2

    if matrix_form:
        return asymm_matrix

    return np.array(
        [asymm_matrix[1][2], asymm_matrix[2][0], asymm_matrix[0][1]],
        dtype=float,
    )


def from_iso(iso):
    r"""
    Computes matrix form of the isotropic exchange parameter.

    .. math::

        C_{2,2}
        J_{iso}
        \boldsymbol{S}_{\mu}
        \cdot
        \boldsymbol{S}_{\nu}
        =
        C_{2,2}
        \boldsymbol{S}_{\mu}
        \boldsymbol{J}_{iso}
        \boldsymbol{S}_{\nu}

    where matrix :math:`\boldsymbol{J}_{iso}` is defined as

    .. math::

        \boldsymbol{J}_{iso}
        =
        \begin{pmatrix}
            J_{iso} & 0 & 0 \\
            0 & J_{iso} & 0 \\
            0 & 0 & J_{iso} \\
        \end{pmatrix}

    Parameters
    ----------
    iso : int or float
        Isotropic exchange parameter.

    Returns
    -------
    parameter : (3, 3) :numpy:`ndarray`
        Matrix form of the isotropic exchange parameter (:math:`\boldsymbol{J}_{iso}`).

    See Also
    --------
    to_iso
    from_dmi

    Examples
    --------

    .. doctest::

        >>> from magnopy import converter22
        >>> converter22.from_iso(iso=1)
        array([[1., 0., 0.],
               [0., 1., 0.],
               [0., 0., 1.]])
    """

    return iso * np.eye(3, dtype=float)


def from_dmi(dmi):
    r"""
    Computes matrix form of the Dzyaloshinsky-Moriya interaction parameter.

    .. math::


        C_{2,2}
        \boldsymbol{D}
        \cdot
        \left(
        \boldsymbol{S}_{\mu}
        \times
        \boldsymbol{S}_{\nu}
        \right)
        =
        C_{2,2}
        \boldsymbol{S}_{\mu}
        \boldsymbol{J}_A
        \boldsymbol{S}_{\nu}


    where matrix :math:`\boldsymbol{J}_A` is defined as

    .. math::

        \boldsymbol{J}_A
        =
        \begin{pmatrix}
            0 & D^z & -D^y \\
            -D^z & 0 & D^x \\
            D^y & -D^x & 0 \\
        \end{pmatrix}

    Parameters
    ----------
    dmi : (3,) |array-like|_
        Vector of Dzyaloshinsky-Moriya interaction parameter :math:`(D_x, D_y, D_z)`.

    Returns
    -------
    parameter : (3, 3) :numpy:`ndarray`
        Matrix form of the Dzyaloshinsky-Moriya interaction parameter
        (:math:`\boldsymbol{J}_A`).

    See Also
    --------
    to_dmi
    from_iso


    Examples
    --------

    .. doctest::

        >>> from magnopy import converter22
        >>> converter22.from_dmi(dmi = (1, 2, 0))
        array([[ 0.,  0., -2.],
               [ 0.,  0.,  1.],
               [ 2., -1.,  0.]])
    """

    return np.array(
        [
            [0, dmi[2], -dmi[1]],
            [-dmi[2], 0, dmi[0]],
            [dmi[1], -dmi[0], 0],
        ],
        dtype=float,
    )


# Populate __all__ with objects defined in this file
__all__ = list(set(dir()) - old_dir)
# Remove all semi-private objects
__all__ = [i for i in __all__ if not i.startswith("_")]
del old_dir
