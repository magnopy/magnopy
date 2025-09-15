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

INDICES = [
    (0, 0, 0, 0),
    (0, 1, 0, 1),
    (0, 2, 0, 2),
    (1, 0, 1, 0),
    (1, 1, 1, 1),
    (1, 2, 1, 2),
    (2, 0, 2, 0),
    (2, 1, 2, 1),
    (2, 2, 2, 2),
]

# Save local scope at this moment
old_dir = set(dir())
old_dir.add("old_dir")


def to_biquadratic(parameter, tensor_form=False):
    r"""
    Computes biquadratic exchange parameter from full matrix parameter.


    .. math::

        C_{4,2,2}
        \sum_{i,j,u,v}
        J^{ijuv}
        S_{\mu}^i
        S_{\mu}^j
        S_{\nu}^u
        S_{\nu}^v
        =
        C_{4,2,2}B
        \left(
        \boldsymbol{S}_{\mu}
        \cdot
        \boldsymbol{S}_{\nu}
        \right)^2
        +
        \dots

    where :math:`B` is defined as

    .. math::

        B = \dfrac{J^{xxxx} + J^{xyxy} + J^{xzxz} +
        J^{yxyx} + J^{yyyy} + J^{yzyz} +
        J^{zxzx} + J^{zyzy} + J^{zzzz}}{9}


    Parameters
    ----------
    parameter : (3, 3, 3, 3) |array-like|_
        Full tensor parameter (:math:`\boldsymbol{J}`).
    tensor_form : bool, default False
        Whether to return tensor form of biquadratic exchange parameter instead of the scalar.

    Returns
    -------
    B : float or (3, 3, 3, 3) :numpy:`ndarray`
        Biquadratic exchange parameter.

        * If ``tensor_form == False``, then returns a number :math:`B`.
        * If ``tensor_form == True``, then returns an array :math:`\boldsymbol{J}_B`.

    See Also
    --------
    from_biquadratic

    Examples
    --------

    .. doctest::

        >>> import numpy as np
        >>> from magnopy import converter422
        >>> parameter = np.ones((3, 3, 3, 3))
        >>> B = converter422.to_biquadratic(parameter=parameter)
        >>> B
        1.0
        >>> J_B = converter422.to_biquadratic(parameter=parameter, tensor_form=True)
        >>> J_B
        array([[[[1., 0., 0.],
                 [0., 0., 0.],
                 [0., 0., 0.]],
        <BLANKLINE>
                [[0., 1., 0.],
                 [0., 0., 0.],
                 [0., 0., 0.]],
        <BLANKLINE>
                [[0., 0., 1.],
                 [0., 0., 0.],
                 [0., 0., 0.]]],
        <BLANKLINE>
        <BLANKLINE>
               [[[0., 0., 0.],
                 [1., 0., 0.],
                 [0., 0., 0.]],
        <BLANKLINE>
                [[0., 0., 0.],
                 [0., 1., 0.],
                 [0., 0., 0.]],
        <BLANKLINE>
                [[0., 0., 0.],
                 [0., 0., 1.],
                 [0., 0., 0.]]],
        <BLANKLINE>
        <BLANKLINE>
               [[[0., 0., 0.],
                 [0., 0., 0.],
                 [1., 0., 0.]],
        <BLANKLINE>
                [[0., 0., 0.],
                 [0., 0., 0.],
                 [0., 1., 0.]],
        <BLANKLINE>
                [[0., 0., 0.],
                 [0., 0., 0.],
                 [0., 0., 1.]]]])

    """

    parameter = np.array(parameter)

    if parameter.shape != (3, 3, 3, 3):
        raise ValueError(f"Wrong shape of the parameter, got {parameter.shape}.")

    B = 0.0

    for i, j, u, v in INDICES:
        B += parameter[i, j, u, v]

    B /= 9

    if tensor_form:
        return from_biquadratic(B=B)

    return float(B)


def from_biquadratic(B):
    r"""
    Computes tensor form of the biquadratic exchange parameter.


    .. math::

        C_{4,2,2}
        B
        \left(
        \boldsymbol{S}_{\mu}
        \cdot
        \boldsymbol{S}_{\nu}
        \right)^2
        =
        C_{4,2,2}
        \sum_{i,j,u,v}
        J_B^{ijuv}
        S_{\mu}^i
        S_{\mu}^j
        S_{\nu}^u
        S_{\nu}^v


    where tensor :math:`\boldsymbol{J}_B` is defined as

    *   :math:`J_B^{ijuv} = B` if :math:`(ijuv)` is one of

        .. math::
            \begin{matrix}
                (xxxx), & (xyxy), & (xzxz), \\
                (yxyx), & (yyyy), & (yzyz), \\
                (zxzx), & (zyzy), & (yyyy)
            \end{matrix}

    *   :math:`J_B^{ijuv} = 0` otherwise.


    Parameters
    ----------
    B : float
        Biquadratic exchange parameter.

    Returns
    -------
    parameter : (3, 3, 3, 3) :numpy:`ndarray`
        Tensor form of the biquadratic exchange parameter.

    See Also
    --------
    to_biquadratic

    Examples
    --------

    .. doctest::

        >>> from magnopy import converter422
        >>> parameter = converter422.from_biquadratic(B=1)
        >>> parameter.shape
        (3, 3, 3, 3)
        >>> parameter
        array([[[[1., 0., 0.],
                 [0., 0., 0.],
                 [0., 0., 0.]],
        <BLANKLINE>
                [[0., 1., 0.],
                 [0., 0., 0.],
                 [0., 0., 0.]],
        <BLANKLINE>
                [[0., 0., 1.],
                 [0., 0., 0.],
                 [0., 0., 0.]]],
        <BLANKLINE>
        <BLANKLINE>
               [[[0., 0., 0.],
                 [1., 0., 0.],
                 [0., 0., 0.]],
        <BLANKLINE>
                [[0., 0., 0.],
                 [0., 1., 0.],
                 [0., 0., 0.]],
        <BLANKLINE>
                [[0., 0., 0.],
                 [0., 0., 1.],
                 [0., 0., 0.]]],
        <BLANKLINE>
        <BLANKLINE>
               [[[0., 0., 0.],
                 [0., 0., 0.],
                 [1., 0., 0.]],
        <BLANKLINE>
                [[0., 0., 0.],
                 [0., 0., 0.],
                 [0., 1., 0.]],
        <BLANKLINE>
                [[0., 0., 0.],
                 [0., 0., 0.],
                 [0., 0., 1.]]]])

    """

    parameter = np.zeros((3, 3, 3, 3), dtype=float)

    for i, j, u, v in INDICES:
        parameter[i, j, u, v] = B

    return parameter


# Populate __all__ with objects defined in this file
__all__ = list(set(dir()) - old_dir)
# Remove all semi-private objects
__all__ = [i for i in __all__ if not i.startswith("_")]
del old_dir
