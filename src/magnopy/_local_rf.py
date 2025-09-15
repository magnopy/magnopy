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


def span_local_rf(direction_vector, hybridize=False):
    r"""
    Span local right-handed reference frame based on the direction vector.

    Parameters
    ----------
    direction_vector : (3, ) |array-like|_
        Direction of the z axis of the local reference frame.
    hybridize : bool, default False

        * If ``hybridize == True``, then returns ``p_alpha, z_alpha``.
        * If ``hybridize == False``, then returns ``x_alpha, y_alpha, z_alpha``.

    Returns
    -------
    x_alpha : (3, ) :numpy:`ndarray`
    y_alpha : (3, ) :numpy:`ndarray`
    p_alpha : (3, ) :numpy:`ndarray`
        ``p_alpha = x_alpha + 1j * y_alpha``.
    z_alpha : (3, ) :numpy:`ndarray`

    See Also
    --------
    span_local_rfs

    Examples
    --------

    Two special cases are handled as

    .. doctest::

        >>> import magnopy
        >>> x, y, z = magnopy.span_local_rf([0, 0, 1])
        >>> x
        array([1., 0., 0.])
        >>> y
        array([0., 1., 0.])
        >>> z
        array([0., 0., 1.])
        >>> p, z = magnopy.span_local_rf([0, 0, 1], hybridize=True)
        >>> p
        array([1.+0.j, 0.+1.j, 0.+0.j])
        >>> z
        array([0., 0., 1.])

    .. doctest::

        >>> import magnopy
        >>> x, y, z = magnopy.span_local_rf([0, 0, -1])
        >>> x
        array([ 0., -1.,  0.])
        >>> y
        array([-1.,  0.,  0.])
        >>> z
        array([ 0.,  0., -1.])
        >>> p, z = magnopy.span_local_rf([0, 0, -1], hybridize=True)
        >>> p
        array([ 0.-1.j, -1.+0.j,  0.+0.j])
        >>> z
        array([ 0.,  0., -1.])

    For the arbitrary direction the global reference frame is rotated as a whole

    .. doctest::

        >>> import magnopy
        >>> x, y, z = magnopy.span_local_rf([1, 1, 1])
        >>> x
        array([ 0.78867513, -0.21132487, -0.57735027])
        >>> y
        array([-0.21132487,  0.78867513, -0.57735027])
        >>> z
        array([0.57735027, 0.57735027, 0.57735027])
        >>> p, z = magnopy.span_local_rf([1, 1, 1], hybridize=True)
        >>> p
        array([ 0.78867513-0.21132487j, -0.21132487+0.78867513j,
               -0.57735027-0.57735027j])
        >>> z
        array([0.57735027, 0.57735027, 0.57735027])

    """

    direction_vector = np.array(direction_vector, dtype=float)

    if np.allclose(direction_vector, np.zeros(3)):
        raise ValueError("Zero vector.")

    direction_vector /= np.linalg.norm(direction_vector)

    if np.allclose(direction_vector, [0, 0, 1]):
        x_alpha = np.array([1, 0, 0], dtype=float)
        y_alpha = np.array([0, 1, 0], dtype=float)
    elif np.allclose(direction_vector, [0, 0, -1]):
        x_alpha = np.array([0, -1, 0], dtype=float)
        y_alpha = np.array([-1, 0, 0], dtype=float)
    else:
        z_dir = [0, 0, 1]

        sin_rot_angle = np.linalg.norm(np.cross(z_dir, direction_vector))
        cos_rot_angle = np.dot(z_dir, direction_vector)
        # direction_vector and z_dir are unit vectors
        ux, uy, uz = np.cross(z_dir, direction_vector) / sin_rot_angle

        x_alpha = np.array(
            [
                ux**2 * (1 - cos_rot_angle) + cos_rot_angle,
                ux * uy * (1 - cos_rot_angle) + uz * sin_rot_angle,
                ux * uz * (1 - cos_rot_angle) - uy * sin_rot_angle,
            ]
        )

        y_alpha = np.array(
            [
                ux * uy * (1 - cos_rot_angle) - uz * sin_rot_angle,
                uy**2 * (1 - cos_rot_angle) + cos_rot_angle,
                uy * uz * (1 - cos_rot_angle) + ux * sin_rot_angle,
            ]
        )

    if hybridize:
        return x_alpha + 1j * y_alpha, direction_vector
    return x_alpha, y_alpha, direction_vector


def span_local_rfs(directional_vectors, hybridize=False):
    r"""
    Span a series of local right-handed reference frames based on a series of the
    direction vectors.

    Parameters
    ----------
    direction_vectors : (M, 3) |array-like|_
        Direction of the z axis of the local reference frames.
    hybridize : bool, default False

        * If ``hybridize == True``, then returns ``p_alphas, z_alphas``.
        * If ``hybridize == False``, then returns ``x_alphas, y_alphas, z_alphas``.

    Returns
    -------
    x_alphas : (M, 3) :numpy:`ndarray`
    y_alphas : (M, 3) :numpy:`ndarray`
    p_alphas : (M, 3) :numpy:`ndarray`
        ``p_alpha = x_alpha + 1j * y_alpha``.
    z_alphas : (M, 3) :numpy:`ndarray`

    See Also
    --------
    span_local_rf

    Examples
    --------


    .. doctest::

        >>> import magnopy
        >>> x, y, z = magnopy.span_local_rfs([[0, 0, 1], [0, 0, -1], [1, 1, 1]])
        >>> x
        array([[ 1.        ,  0.        ,  0.        ],
               [ 0.        , -1.        ,  0.        ],
               [ 0.78867513, -0.21132487, -0.57735027]])
        >>> y
        array([[ 0.        ,  1.        ,  0.        ],
               [-1.        ,  0.        ,  0.        ],
               [-0.21132487,  0.78867513, -0.57735027]])
        >>> z
        array([[ 0.        ,  0.        ,  1.        ],
               [ 0.        ,  0.        , -1.        ],
               [ 0.57735027,  0.57735027,  0.57735027]])
    """

    results = []

    for directional_vector in directional_vectors:
        results.append(
            span_local_rf(direction_vector=directional_vector, hybridize=hybridize)
        )

    results = np.array(results)

    if hybridize:
        return results[:, 0], results[:, 1]

    return results[:, 0], results[:, 1], results[:, 2]


# Populate __all__ with objects defined in this file
__all__ = list(set(dir()) - old_dir)
# Remove all semi-private objects
__all__ = [i for i in __all__ if not i.startswith("_")]
del old_dir
