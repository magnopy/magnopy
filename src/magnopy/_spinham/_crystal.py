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
from dataclasses import dataclass
import numpy as np


def _as_array(value, name, shape, dtype) -> np.ndarray:
    """
    Coerce ``value`` to a plain ``numpy.ndarray`` of ``dtype`` and ``shape``.

    Coercion and shape check only -- no value level validation.

    Any failure is re-raised as ``ValueError``.

    Parameters
    ----------
    value : |array-like|_
    name : str
        Field name, used in error messages.
    shape : tuple of int
        Required shape. Checked exactly.
    dtype : type
        Target type passed directly to ``numpy.asarray``.

    Returns
    -------
    array : (``shape``) :numpy:`ndarray` of type ``dtype``
        May alias ``value`` if it was already a matching ndarray.

    Raises
    ------
    ValueError
        If ``value`` cannot be interpreted as ``dtype``, or its shape differs from
        ``shape``.
    """

    try:
        array = np.asarray(value, dtype=dtype)
    except (ValueError, TypeError) as error:
        raise ValueError(
            f"{name} could not be interpreted as an array of {np.dtype(dtype)} of "
            f"shape {shape}. Got {type(value).__name__}. (numpy: {error})"
        ) from error

    if array.shape != shape:
        raise ValueError(f"{name} must have shape {shape}, got {array.shape}.")

    return array


# Numerical tolerance for checking whether cell is singular
# Scale-free
_SINGULAR_TOL = 1e-10


# eq=False, generated __eq__ compares fields with ==,
# which returns an array for ndarrays, not a bool.
# __eq__ is supplied
@dataclass(frozen=True, eq=False)
class _Crystal:
    r"""
    Immutable crystal.

    Parameters
    ----------
    cell : (3, 3) |array-like|_ of float
        Lattice vectors, stored as rows. Checked to be finite and non-singular.
    names : iterable of str
        Labels for atoms. Fixes the number of atoms ``N``. Other fields are checked
        against ``len(names)``. Names enter ``__eq__`` and ``__hash__``;
        ``is_close`` compares them exactly.

        |wulfric|_ derives symmetry-equivalent sites from names when ``spglib_types``
        are not supplied as user-data on :py:class:`.SpinHamiltonian`.
    positions : (N, 3) |array-like|_ of float
        Relative (i.e. fractional) coordinates of atoms. Absolute coordinates in real
        space would be ``positions @ cell``. Not wrapped into [0, 1), an atom sitting
        outside of the (0, 0, 0) unit cell is allowed.
        Checked to be finite.
    spins : (N,) |array-like|_ of float
        Spin magnitudes, checked to be finite and non-negative.
    g_factors : (N,) |array-like|_ of float
        Lande g-factors, checked to be finite and non-zero.
    magnetic : (N,) |array-like|_ of bool or 0/1
        Which sites enter a :py:class:`.SpinHamiltonian`. Declared by the user/loader
        explicitly, not inferred from interaction parameters. Non-magnetic sites are
        kept, since they define the crystal symmetry, even though they do not enter
        Magnopy's calculations (e.g. :py:class:`.LSWT` or :py:class:`.Energy`).

    Raises
    ------
    ValueError
        - If ``cell`` is singular.
        - If any of ``cell``, ``positions``, ``spins`` or ``g_factors`` contains a
          non-finite value.
        - If any spin is negative.
        - If any g-factor is zero.
        - If no atom is magnetic.
        - If the crystal is empty (``N==0``).
        - If length of ``positions``, ``spins``, ``g_factors`` or ``magnetic`` is not
          equal to ``len(names)`` (N).
        - If any of ``cell``, ``positions``, ``spins``, ``g_factors`` or ``magnetic``
          has an unexpected shape.
        - If any element of ``magnetic`` is not a bool or not 0/1.

    Notes
    -----
    Constructor arguments are stored as attributes of the same name.

    ``names`` is stored as tuple of str.

    All array fields are copied at construction, so mutation of the array you passed
    in has no effect on the crystal. The stored arrays are read-only:
    ``crystal.positions[0] = [0, 0, 0]`` raises ``ValueError``.

    ``__eq__`` and ``__hash__`` are exact (bitwise). ``is_close`` is a tolerant
    comparison. ``__eq__`` answers *is this the same value?*, ``is_close`` answers
    *is this the same material?*.

    Two crystals describing the same material, built through different routes, may
    differ in their float-valued attributes due to rounding errors (with no physical
    significance), therefore ``__eq__`` compares them unequal. ``is_close`` is
    introduced to handle that case.

    ``is_close`` cannot serve as ``__eq__``. Dictionaries and sets require equality to
    be an equivalence relation, and a tolerant comparison is not one: three crystals
    might satisfy ``a.is_close(b)`` and ``b.is_close(c)`` while ``a`` and ``c`` are not
    close.
    """

    cell: np.ndarray
    names: tuple[str, ...]
    positions: np.ndarray  # (N, 3) float, relative
    spins: np.ndarray  # (N,)   float
    g_factors: np.ndarray  # (N,)   float
    magnetic: np.ndarray  # (N,)   bool

    def __post_init__(self):

        ################################## cell ################################
        cell = _as_array(value=self.cell, name="cell", shape=(3, 3), dtype=float)
        if not np.isfinite(cell).all():
            raise ValueError("Cell contains non-finite elements.")

        scale = float(np.abs(cell).max())
        det = float(np.linalg.det(cell))
        if scale == 0.0 or abs(det) < _SINGULAR_TOL * scale**3:
            raise ValueError(
                f"cell must be non-singular, got det = {det:.3e} "
                f"|det|/max|cell|^3 = {abs(det) / scale**3 if scale else 0.0:.3e}; "
                f"the three lattice vectors must be linearly independent."
            )

        # Own it before freezing it
        cell = np.array(cell, copy=True)
        cell += 0.0  # -0.0 -> +0.0
        cell.flags.writeable = False
        object.__setattr__(self, "cell", cell)
