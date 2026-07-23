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
from functools import cached_property
import numpy as np


# Dimensionless numerical tolerance for checking whether the cell is singular.
# |det(cell)| is a sizable fraction of the |max(cell)|^3 (its bounding box).
# if the ration of det/box is too small (< _SINGULAR_TOL), then the cell is
# considered to be singular. The tolerance is relative by design.
_SINGULAR_TOL = 1e-10


def _as_float_array(value, name, shape) -> np.ndarray:
    """
    Coerce ``value`` to a plain ``numpy.ndarray`` of float and ``shape``.

    Coercion and shape check only -- no value level validation.

    Any failure is re-raised as ``ValueError``.

    Parameters
    ----------
    value : |array-like|_
    name : str
        Field name, used in error messages.
    shape : tuple of int
        Required shape. Checked exactly.

    Returns
    -------
    array : (``shape``) :numpy:`ndarray` of type float
        May alias ``value`` if it was already a matching ndarray.

    Raises
    ------
    ValueError
        If ``value`` cannot be interpreted as array of floats, or its shape differs
        from ``shape``, or it contains non-finite elements.
    """

    try:
        array = np.asarray(value, dtype=float)
    except (ValueError, TypeError) as error:
        raise ValueError(
            f"{name} could not be interpreted as an array of float. "
            f"Got {type(value).__name__}. (numpy: {error})"
        ) from error

    if array.shape != shape:
        raise ValueError(f"{name} must have shape {shape}, got {array.shape}.")

    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} contains non-finite elements.")

    return array


def _as_bool_array(value, name, shape) -> np.ndarray:
    """
    Coerce ``value`` to a plain ``numpy.ndarray`` of bool and ``shape``.

    Coercion and shape check only -- no value level validation.

    Any failure is re-raised as ``ValueError``.

    Parameters
    ----------
    value : |array-like|_
    name : str
        Field name, used in error messages.
    shape : tuple of int
        Required shape. Checked exactly.

    Returns
    -------
    array : (``shape``) :numpy:`ndarray` of type bool
        May alias ``value`` if it was already a matching ndarray.

    Raises
    ------
    ValueError
        If ``value`` cannot be interpreted as array of bools, or its shape differs
        from ``shape``, or if its elements are not in ``[True, False, 1, 0]``.
    """

    ALLOWED_ELEMENTS = [True, False, 1, 0]

    try:
        array = np.asarray(value)
    except (ValueError, TypeError) as error:
        raise ValueError(
            f"{name} could not be interpreted as an array (numpy: {error})."
        ) from error

    if array.shape != shape:
        raise ValueError(f"{name} must have shape {shape}, got {array.shape}.")

    if not np.all(np.isin(array, ALLOWED_ELEMENTS)):
        bad = np.unique(array[~np.isin(array, ALLOWED_ELEMENTS)])
        raise ValueError(f"{name} must be 0, 1 or bool, got {bad.tolist()}.")

    return array.astype(bool)


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

    """

    cell: np.ndarray
    names: tuple[str, ...]
    positions: np.ndarray  # (N, 3) float, relative
    spins: np.ndarray  # (N,)   float
    g_factors: np.ndarray  # (N,)   float
    magnetic: np.ndarray  # (N,)   bool

    def __post_init__(self):

        ############################# types & shape ############################
        names = tuple(self.names)
        N = len(names)
        if N == 0:
            raise ValueError(
                "The crystal is empty (N == 0), must have at least one atom."
            )

        cell = _as_float_array(value=self.cell, name="cell", shape=(3, 3))
        positions = _as_float_array(
            value=self.positions, name="positions", shape=(N, 3)
        )
        spins = _as_float_array(value=self.spins, name="spins", shape=(N,))
        g_factors = _as_float_array(value=self.g_factors, name="g-factors", shape=(N,))

        magnetic = _as_bool_array(self.magnetic, name="magnetic", shape=(N,))

        ############################### validate ###############################
        scale = float(np.abs(cell).max())
        det = float(np.linalg.det(cell))
        if scale == 0.0 or abs(det) < _SINGULAR_TOL * scale**3:
            raise ValueError(
                f"cell must be non-singular, got det = {det:.3e}; "
                f"|det|/max|cell|^3 = {abs(det) / scale**3 if scale else 0.0:.3e}; "
                f"the three lattice vectors must be linearly independent."
            )

        for index, name in enumerate(names):
            if not isinstance(name, str):
                raise ValueError(
                    f"names[{index}] is not a str, got {type(name).__name__}"
                )

        if np.any(spins < 0):
            raise ValueError("spins contain negative elements")

        if np.any(g_factors == 0):
            raise ValueError("g_factors contain zero-valued elements")

        if not np.any(magnetic):
            raise ValueError("At least one atom must be magnetic, got none.")

        ############################# copy & freeze ############################
        object.__setattr__(self, "names", names)
        for name, array in [
            ("cell", cell),
            ("positions", positions),
            ("spins", spins),
            ("g_factors", g_factors),
            ("magnetic", magnetic),
        ]:
            # Own it before freezing it
            if name == "magnetic":
                array = np.array(array, dtype=bool, copy=True)
            else:
                array = np.array(array, copy=True, dtype=float)
                # -0.0 -> +0.0
                array += 0.0
            array.flags["WRITEABLE"] = False
            object.__setattr__(self, name, array)

    ############################################################################
    #                                 __repr__                                 #
    ############################################################################

    def __repr__(self):
        a, b, c = np.linalg.norm(self.cell, axis=1)
        return (
            f"_Crystal(a={a:.4e}, b={b:.4e}, c={c:.4e}; "
            f"{self.M_prime} atoms in total, {self.M} magnetic atoms)"
        )

    ############################################################################
    #                                identities                                #
    ############################################################################

    @cached_property
    def _hash_key(self):
        return hash(
            (
                self.names,
                self.cell.tobytes(),
                self.positions.tobytes(),
                self.spins.tobytes(),
                self.g_factors.tobytes(),
                self.magnetic.tobytes(),
            )
        )

    def __hash__(self):
        return self._hash_key

    def __eq__(self, other):
        r"""
        Exact (bitwise) comparison, consistent with :py:meth:`.__hash__`.

        Answers *is this the same value?*. For *is this the same material?* use
        :py:meth:`.is_close`.
        """

        if not isinstance(other, _Crystal):
            return NotImplemented

        return (
            self._hash_key == other._hash_key  # cheap fail
            and self.names == other.names
            and self.cell.tobytes() == other.cell.tobytes()
            and self.positions.tobytes() == other.positions.tobytes()
            and self.spins.tobytes() == other.spins.tobytes()
            and self.g_factors.tobytes() == other.g_factors.tobytes()
            and self.magnetic.tobytes() == other.magnetic.tobytes()
        )

    def is_close(self, other, cell_rtol=1e-10, atol=1e-8) -> bool:
        r"""
        Tolerant comparison of two crystals.

        Answers the question *Is it the same crystal?*.

        Parameters
        ----------
        other : :py:class:`._Crystal`
            Other crystal that will be compared.
        cell_rtol : float, default 1e-10
            Relative tolerance for cell.
        atol : float, default 1e-8
            Absolute tolerance for positions, spins and g-factors

        Returns
        -------
        result : bool
            ``True`` if two crystals appear to be the same within tolerance.
            ``False`` if not.

        Notes
        -----
        See :py:meth:`.diff` for the details on how exactly the comparison is made.

        ``__eq__`` and ``__hash__`` are exact (bitwise). ``is_close`` is a tolerant
        comparison. ``__eq__`` answers *is this the same value?*, ``is_close`` answers
        *is this the same material?*.

        Two crystals describing the same material, built through different routes, may
        differ in their float-valued attributes due to rounding errors (with no physical
        significance), therefore ``__eq__`` compares them unequal. ``is_close`` is
        introduced to handle that case.

        ``is_close`` cannot serve as ``__eq__``. Dictionaries and sets require equality
        to be an equivalence relation, and a tolerant comparison is not one: three
        crystals might satisfy ``a.is_close(b)`` and ``b.is_close(c)`` while ``a`` and
        ``c`` are not close.

        See Also
        --------
        diff
        """

        return self.diff(other=other, cell_rtol=cell_rtol, atol=atol) == ()

    def diff(self, other, cell_rtol=1e-10, atol=1e-8) -> tuple:
        r"""
        Describes every field in which two crystals differ physically.

        * ``cell`` -- relative only (``atol=0``). It is dimensioned and magnopy
          does not know the units, so only relative agreement is meaningful.
        * ``positions`` -- absolute, after wrapping the difference into
          [-0.5, 0.5): fractional coordinates are periodic, so 0.999999... and
          -1e-15 are the same site.
        * ``spins``, ``g_factors`` -- absolute. Bounded quantities where zero is
          a legal value, so a relative tolerance would be undefined there.
        * ``names``, ``magnetic`` -- exact. ``is_close`` is tolerant about
          floats, not about what the atoms are.

        Parameters
        ----------
        other : :py:class:`._Crystal`
            Other crystal that will be compared.
        cell_rtol : float, default 1e-10
            Relative tolerance for cell.
        atol : float, default 1e-8
            Absolute tolerance for positions, spins and g-factors

        Returns
        -------
        report : tuple of str
            A tuple with the messages about crystal differences.

        See Also
        --------
        is_close
        """

        if not isinstance(other, _Crystal):
            raise TypeError(
                f"Cannot compare a crystal to non-crystal, got {type(other)}"
            )

        report = []

        if len(self.names) != len(other.names):
            return (
                "Two crystals have different amount of atoms: "
                f"{len(self.names)} and {len(other.names)}",
            )

        if not np.allclose(self.cell, other.cell, atol=0.0, rtol=cell_rtol):
            worst = float(np.abs(self.cell - other.cell).max())
            report.append(f"cell: differs by up to {worst:.3e} (absolute)")

        if self.names != other.names:
            for index, (mine, theirs) in enumerate(zip(self.names, other.names)):
                if mine != theirs:
                    report.append(f"names: atom {index} is {mine!r} vs {theirs!r}")
                    break

        if not np.array_equal(self.magnetic, other.magnetic):
            mismatch = np.flatnonzero(self.magnetic != other.magnetic).tolist()
            report.append(f"magnetic: differs at atoms {mismatch}")

        for name in ("spins", "g_factors"):
            mine, theirs = getattr(self, name), getattr(other, name)
            if not np.allclose(mine, theirs, rtol=0.0, atol=atol):
                index = int(np.abs(mine - theirs).argmax())
                report.append(
                    f"{name}: atom {index} is {mine[index]} vs {theirs[index]}"
                )

        delta = self.positions - other.positions
        # Wrap the difference into [-0.5, 0.5)
        delta -= np.round(delta)
        if np.abs(delta).max() > atol:
            index = int(np.abs(delta).max(axis=1).argmax())
            report.append(
                f"positions: atom {index} differs by "
                f"{np.abs(delta[index]).max():.3e} (fractional)"
            )

        return tuple(report)

    ############################################################################
    #                             Derived properties                           #
    ############################################################################
    @cached_property
    def M_prime(self) -> int:
        r"""
        Total amount of atoms in the crystal.

        Returns
        -------
        M_prime : int
        """

        return len(self.names)

    def __len__(self):
        r"""Total amount of atoms, see M_prime."""
        return self.M_prime

    @cached_property
    def M(self) -> int:
        r"""
        Amount of *magnetic* atoms in the crystal.

        Returns
        -------
        M : int

        Notes
        -----
        Atom ``i`` is magnetic if ``crystal.magnetic[i]`` is ``True``.
        """

        return int(np.count_nonzero(self.magnetic))

    @cached_property
    def map_to_all(self) -> np.ndarray:
        r"""
        Which atom is the i-th magnetic atom?

        Returns
        -------
        map_to_all : :numpy:`ndarray` of int
            Array of indices, such that ``map_to_all[i]`` is an index in atom of the
            i-th magnetic atom.
        """

        result = np.flatnonzero(self.magnetic)
        result.flags["WRITEABLE"] = False

        return result

    @cached_property
    def map_to_magnetic(self) -> np.ndarray:
        r"""
        Which magnetic atom is the i-th atom?

        Returns
        -------
        map_to_magnetic : :numpy:`ndarray` of int
            Array of indices, such that ``map_to_magnetic[i]`` is an index of i-th
            atom among the magnetic atoms. -1 if the atom is not magnetic.
        """

        result = np.full(self.M_prime, -1, dtype=int)
        result[self.map_to_all] = np.arange(self.M)
        result.flags["WRITEABLE"] = False

        return result
