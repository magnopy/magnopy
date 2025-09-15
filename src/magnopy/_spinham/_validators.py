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


def _validate_atom_index(index, atoms) -> None:
    r"""
    Validate that the atom index is in agreement with the amount of atoms

    Parameter
    ---------
    index
        Potential index of an atom.
    atoms : dict
        Dictionary with the atoms. This function relies on ``atoms["names"]``.

    Raises
    ------
    TypeError
        If ``index`` is not an integer.
    ValueError
        If ``index`` is out of range.
    """

    if not isinstance(index, int):
        raise TypeError(
            f"Only integers are supported as atom indices, "
            f"got '{type(index)}' from '{index}'"
        )

    if not 0 <= index < len(atoms["names"]):
        raise ValueError(
            "Index should be greater or equal to 0 and less than "
            f"{len(atoms['names'])}', got {index}."
        )


def _validate_unit_cell_index(ijk) -> None:
    r"""
    Validate that ijk can specify unit cell.

    Parameters
    ----------
    ijk
        Potential index of the unit cell.

    Raises
    ------
    TypeError
        If ``ijk`` is not a ``tuple``.
    TypeError
        If either ``i``, ``j`` or ``k`` is not an ``int``.
    """

    if not len(ijk) == 3:
        raise TypeError(
            f"Unit cell index has to be a tuple or a list of the length 3, "
            f"got '{len(ijk)}'"
        )

    if not isinstance(ijk[0], int):
        raise TypeError(
            f"First element of the unit cell index is not an 'int', got "
            f"{type(ijk[0])} from '{ijk[0]}'"
        )

    if not isinstance(ijk[1], int):
        raise TypeError(
            f"Second element of the unit cell index is not an 'int', got "
            f"{type(ijk[1])} from '{ijk[1]}'"
        )

    if not isinstance(ijk[2], int):
        raise TypeError(
            f"Third element of the unit cell index is not an 'int', got "
            f"{type(ijk[2])} from '{ijk[2]}'"
        )


def _spins_ordered(mu1, alpha1, mu2, alpha2) -> bool:
    r"""
    Compare two spins based on their positions.

    For the definition of comparison see
    :ref:`user-guide_theory-behind_multiple-counting`.

    Parameters
    ----------
    mu1 : tuple of 3 int
    alpha1 : int
    mu2 : tuple of 3 int
    alpha2 : int

    Returns
    -------
    result : bool
    """

    i1, j1, k1 = mu1
    i2, j2, k2 = mu2

    i = i2 - i1
    j = j2 - j1
    k = k2 - k1

    if (
        i > 0
        or (i == 0 and j > 0)
        or (i == 0 and j == 0 and k > 0)
        or (i == 0 and j == 0 and k == 0 and alpha1 < alpha2)
    ):
        return True

    return False
