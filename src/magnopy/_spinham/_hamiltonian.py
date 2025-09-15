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


from copy import deepcopy
from math import ceil

import numpy as np
from wulfric import add_sugar

from magnopy._spinham._c1 import _add_1, _p1, _remove_1
from magnopy._spinham._c21 import _add_21, _p21, _remove_21
from magnopy._spinham._c22 import _add_22, _get_primary_p22, _p22, _remove_22
from magnopy._spinham._c31 import _add_31, _p31, _remove_31
from magnopy._spinham._c32 import _add_32, _p32, _remove_32
from magnopy._spinham._c33 import _add_33, _p33, _remove_33
from magnopy._spinham._c41 import _add_41, _p41, _remove_41
from magnopy._spinham._c43 import _add_43, _p43, _remove_43
from magnopy._spinham._c44 import _add_44, _p44, _remove_44
from magnopy._spinham._c421 import _add_421, _p421, _remove_421
from magnopy._spinham._c422 import _add_422, _p422, _remove_422
from magnopy._spinham._convention import Convention

# Save local scope at this moment
old_dir = set(dir())
old_dir.add("old_dir")


def _merge(list1: list, list2: list) -> list:
    r"""
    Merge two sorted parameter lists for any term.

    Lists of parameters have the form

    .. code-block:: python

        list = [[specs, parameter], ...]

    Comparison is based on specs.

    Parameter
    ---------
    list1 : list
        First list of parameters.
    list2 : list
        Second list of parameters.

    Returns
    -------
    merged_list : list
        Merged list of parameters.
    """

    list1 = deepcopy(list1)
    list2 = deepcopy(list2)

    merged_list = []

    i1 = 0
    i2 = 0

    while i1 < len(list1) or i2 < len(list2):
        if i1 >= len(list1):
            merged_list.append(list2[i2])
            i2 += 1
        elif i2 >= len(list2):
            merged_list.append(list1[i1])
            i1 += 1
        elif list1[i1][:-1] == list2[i2][:-1]:
            merged_list.append(list1[i1])
            merged_list[-1][-1] = merged_list[-1][-1] + list2[i2][-1]
            i1 += 1
            i2 += 1
        elif list1[i1][:-1] < list2[i2][:-1]:
            merged_list.append(list1[i1])
            i1 += 1
        else:
            merged_list.append(list2[i2])
            i2 += 1

    return merged_list


class SpinHamiltonian:
    r"""
    Spin Hamiltonian.

    Parameters
    ----------
    convention : :py:class:`.Convention` or str
        A convention of the spin Hamiltonian.
    cell : (3, 3) |array-like|_
        Matrix of a cell, rows are interpreted as vectors.
    atoms : dict
        Dictionary with atoms.

    Examples
    --------

    For example of usage see the page in the user guide -
    :ref:`user-guide_usage_spin-hamiltonian`.

    """

    def __init__(self, cell, atoms, convention) -> None:
        self._cell = np.array(cell)

        self._atoms = add_sugar(atoms)

        # Only the magnetic sites
        self._magnetic_atoms = None
        self._map_to_magnetic = None
        self._map_to_all = None

        self._convention = convention

        # [[alpha, parameter], ...]
        self._1 = []

        # [[alpha, parameter], ...]
        self._21 = []

        # [[alpha, beta, nu, parameter], ...]
        self._22 = []

        # [[alpha, parameter], ...]
        self._31 = []

        # [[alpha, beta, nu, parameter], ...]
        self._32 = []

        # [[alpha, beta, gamma, nu, lambda, parameter], ...]
        self._33 = []

        # [[alpha, parameter], ...]
        self._41 = []

        # [[alpha, beta, nu, parameter], ...]
        self._421 = []

        # [[alpha, beta, nu, parameter], ...]
        self._422 = []

        # [[alpha, beta, gamma, nu, lambda, parameter], ...]
        self._43 = []

        # [[alpha, beta, gamma, epsilon, nu, lambda, rho, parameter], ...]
        self._44 = []

    ############################################################################
    #                              Cell and Atoms                              #
    ############################################################################

    @property
    def cell(self):
        r"""
        Cell of the crystal on which the Hamiltonian is build.

        Returns
        -------
        cell : (3, 3) :numpy:`ndarray`
            Matrix of the cell, rows are lattice vectors.

        Notes
        -----
        If is not recommended to change the ``cell`` property after the creation of
        :py:class:`.SpinHamiltonian`. In fact an attempt to do so will raise an
        ``AttributeError``:

        .. doctest::

            >>> import numpy as np
            >>> import magnopy
            >>> convention = magnopy.Convention()
            >>> spinham = magnopy.SpinHamiltonian(
            ...     cell=np.eye(3), atoms={}, convention=convention
            ... )
            >>> spinham.cell = 2 * np.eye(3)
            Traceback (most recent call last):
            ...
            AttributeError: Change of the cell attribute is not supported after the creation of SpinHamiltonian instance. If you need to modify cell, then use pre-defined methods of SpinHamiltonian or create a new one.

        Use pre-defined methods of the :py:class:`.SpinHamiltonian` class to safely
        modify the cell.

        If you need to change the cell attribute, then use

        .. doctest::

            >>> import numpy as np
            >>> import magnopy
            >>> convention = magnopy.Convention()
            >>> spinham = magnopy.SpinHamiltonian(
            ...     cell=np.eye(3), atoms={}, convention=convention
            ... )
            >>> spinham.cell
            array([[1., 0., 0.],
                   [0., 1., 0.],
                   [0., 0., 1.]])
            >>> spinham._cell = 2 * np.eye(3)
            >>> spinham.cell
            array([[2., 0., 0.],
                   [0., 2., 0.],
                   [0., 0., 2.]])


        In the latter case correct behavior of magnopy **is not** guaranteed. Use only
        if you have a deep understanding of the magnopy source code.
        """

        return self._cell

    @cell.setter
    def cell(self, new_value):
        raise AttributeError(
            "Change of the cell attribute is not supported after "
            "the creation of SpinHamiltonian instance. If you need to modify cell, "
            "then use pre-defined methods of SpinHamiltonian or create a new one."
        )

    @property
    def atoms(self) -> dict:
        r"""
        Atoms of the crystal on which the Hamiltonian is build.

        Returns
        -------
        atoms : dict (with added sugar)
            Dictionary with the atoms.

        Notes
        -----
        If is not recommended to change the atoms property after the creation of
        :py:class:`.SpinHamiltonian`. In fact an attempt to do so will raise an
        ``AttributeError``:

        .. doctest::

            >>> import numpy as np
            >>> import magnopy
            >>> convention = magnopy.Convention()
            >>> spinham = magnopy.SpinHamiltonian(
            ...     cell=np.eye(3), atoms={}, convention=convention
            ... )
            >>> spinham.atoms = {"names": ["Cr"]}
            Traceback (most recent call last):
            ...
            AttributeError: Change of the atoms dictionary is not supported after the creation of SpinHamiltonian instance. If you need to modify atoms, then use pre-defined methods of SpinHamiltonian or create a new one.

        Use pre-defined methods of the :py:class:`.SpinHamiltonian` class to safely
        modify atoms.

        If you need to change the whole dictionary at once, then use

        .. doctest::

            >>> import numpy as np
            >>> import magnopy
            >>> convention = magnopy.Convention()
            >>> spinham = magnopy.SpinHamiltonian(
            ...     cell=np.eye(3), atoms={}, convention=convention
            ... )
            >>> spinham.atoms
            {}
            >>> spinham._atoms = {"names": ["Cr"]}
            >>> spinham.atoms
            {'names': ['Cr']}

        In the latter case correct behavior of magnopy **is not** guaranteed. Use only
        if you have a deep understanding of the magnopy source code.
        """

        return self._atoms

    @atoms.setter
    def atoms(self, new_value):
        raise AttributeError(
            "Change of the atoms dictionary is not supported after "
            "the creation of SpinHamiltonian instance. If you need to modify atoms, "
            "then use pre-defined methods of SpinHamiltonian or create a new one."
        )

    def _reset_internals(self):
        self._map_to_magnetic = None
        self._map_to_all = None
        self._magnetic_atoms = None

    def _update_internals(self):
        # Identify magnetic sites
        indices = set()

        for alpha, _ in self._1:
            indices.add(alpha)

        for alpha, _ in self._21:
            indices.add(alpha)

        for alpha, beta, _, _ in self._22:
            indices.add(alpha)
            indices.add(beta)

        for alpha, _ in self._31:
            indices.add(alpha)

        for alpha, beta, _, _ in self._32:
            indices.add(alpha)
            indices.add(beta)

        for alpha, beta, gamma, _, _, _ in self._33:
            indices.add(alpha)
            indices.add(beta)
            indices.add(gamma)

        for alpha, _ in self._41:
            indices.add(alpha)

        for alpha, beta, _, _ in self._421:
            indices.add(alpha)
            indices.add(beta)

        for alpha, beta, _, _ in self._422:
            indices.add(alpha)
            indices.add(beta)

        for alpha, beta, gamma, _, _, _ in self._43:
            indices.add(alpha)
            indices.add(beta)
            indices.add(gamma)

        for alpha, beta, gamma, epsilon, _, _, _, _ in self._44:
            indices.add(alpha)
            indices.add(beta)
            indices.add(gamma)
            indices.add(epsilon)

        indices = sorted(list(indices))

        # Create index map from all to magnetic
        self._map_to_magnetic = [None for _ in range(len(self.atoms.names))]
        for i in range(len(indices)):
            self._map_to_magnetic[indices[i]] = i

        # Create index map from magnetic to all
        self._map_to_all = indices

        # Create magnetic atoms dictionary
        self._magnetic_atoms = add_sugar({})
        for key in self.atoms:
            self._magnetic_atoms[key] = []

            for full_index in indices:
                self._magnetic_atoms[key].append(self.atoms[key][full_index])

    @property
    def map_to_magnetic(self):
        r"""
        Index map from all atoms to the magnetic ones.

        Returns
        -------
        map_to_magnetic (L, ) list of int
            Index map. Integers. ``L = len(spinham.atoms.names)``
        """

        if self._map_to_magnetic is None:
            self._update_internals()

        return self._map_to_magnetic

    @property
    def map_to_all(self):
        r"""
        Index map from magnetic atoms to all atoms.

        Returns
        -------
        map_to_all (M, ) list of int
            Index map. Integers. ``M = len(spinham.magnetic_atoms.names)``
        """

        if self._map_to_all is None:
            self._update_internals()

        return self._map_to_all

    @property
    def magnetic_atoms(self):
        r"""
        Magnetic atoms of the spin Hamiltonian.

        Magnetic atom is defined as an atom with at least one parameter associated with
        it.

        This property is dynamically computed at every call.

        Returns
        -------
        magnetic_atoms : list of int
            Indices of magnetic atoms in the ``spinham.atoms``. Sorted.

        See Also
        --------
        M
        """

        if self._magnetic_atoms is None:
            self._update_internals()

        return self._magnetic_atoms

    @property
    def M(self):
        r"""
        Number of spins (magnetic atoms) in the unit cell.

        Returns
        -------
        M : int
            Number of spins (magnetic atoms) in the unit cell.

        See Also
        --------
        magnetic_atoms
        """

        return len(self.magnetic_atoms.names)

    ############################################################################
    #                          Convention properties                           #
    ############################################################################
    @property
    def convention(self) -> Convention:
        r"""
        Convention of the spin Hamiltonian.

        See Also
        --------
        Convention
        """

        return self._convention

    @convention.setter
    def convention(self, new_convention: Convention):
        self._set_multiple_counting(new_convention._multiple_counting)

        self._set_spin_normalization(new_convention._spin_normalized)

        self._set_c1(new_convention._c1)

        self._set_c21(new_convention._c21)
        self._set_c22(new_convention._c22)

        self._set_c31(new_convention._c31)
        self._set_c32(new_convention._c32)
        self._set_c33(new_convention._c33)

        self._set_c41(new_convention._c41)
        self._set_c421(new_convention._c421)
        self._set_c422(new_convention._c422)
        self._set_c43(new_convention._c43)
        self._set_c44(new_convention._c44)

        self._convention = new_convention

    def _set_multiple_counting(self, multiple_counting: bool) -> None:
        if multiple_counting is None or self.convention._multiple_counting is None:
            return

        multiple_counting = bool(multiple_counting)

        if self.convention.multiple_counting == multiple_counting:
            return

        # For (two spins & two sites)

        # It was absent before
        if multiple_counting:
            factor = 0.5
        # It was present before
        else:
            factor = 2.0

        for index in range(len(self._22)):
            self._22[index][3] = self._22[index][3] * factor

        # For (three spins & two sites)

        for index in range(len(self._32)):
            self._32[index][3] = self._32[index][3] * factor

        # For (four spins & two sites (3+1))

        for index in range(len(self._421)):
            self._421[index][3] = self._421[index][3] * factor

        # For (four spins & two sites (2+2))

        for index in range(len(self._422)):
            self._422[index][3] = self._422[index][3] * factor

        # For (three spins & three sites)

        # It was absent before
        if multiple_counting:
            factor = 1 / 6
        # It was present before
        else:
            factor = 6

        for index in range(len(self._33)):
            self._33[index][5] = self._33[index][5] * factor

        # For (four spins & three sites)

        for index in range(len(self._43)):
            self._43[index][5] = self._43[index][5] * factor

        # For (four spins & four sites)

        # It was absent before
        if multiple_counting:
            factor = 1 / 24
        # It was present before
        else:
            factor = 24

        for index in range(len(self._44)):
            self._44[index][7] = self._44[index][7] * factor

    def _set_spin_normalization(self, spin_normalized: bool) -> None:
        if spin_normalized is None or self.convention._spin_normalized is None:
            return

        spin_normalized = bool(spin_normalized)

        if self.convention.spin_normalized == spin_normalized:
            return

        # Before it was not normalized
        if spin_normalized:
            # For (one spin & one site)
            for index in range(len(self._1)):
                alpha = self._1[index][0]
                self._1[index][1] = self._1[index][1] * self.atoms.spins[alpha]
            # For (two spins & one site)
            for index in range(len(self._21)):
                alpha = self._21[index][0]
                self._21[index][1] = self._21[index][1] * self.atoms.spins[alpha] ** 2
            # For (two spins & two sites)
            for index in range(len(self._22)):
                alpha = self._22[index][0]
                beta = self._22[index][1]
                self._22[index][3] = self._22[index][3] * (
                    self.atoms.spins[alpha] * self.atoms.spins[beta]
                )
            # For (three spins & one site)
            for index in range(len(self._31)):
                alpha = self._31[index][0]
                self._31[index][1] = self._31[index][1] * self.atoms.spins[alpha] ** 3
            # For (three spins & two sites)
            for index in range(len(self._32)):
                alpha = self._32[index][0]
                beta = self._32[index][1]
                self._32[index][3] = self._32[index][3] * (
                    self.atoms.spins[alpha] ** 2 * self.atoms.spins[beta]
                )
            # For (three spins & three sites)
            for index in range(len(self._33)):
                alpha = self._33[index][0]
                beta = self._33[index][1]
                gamma = self._33[index][2]
                self._33[index][5] = self._33[index][5] * (
                    self.atoms.spins[alpha]
                    * self.atoms.spins[beta]
                    * self.atoms.spins[gamma]
                )
            # For (four spins & one site)
            for index in range(len(self._41)):
                alpha = self._41[index][0]
                self._41[index][1] = self._41[index][1] * self.atoms.spins[alpha] ** 4
            # For (four spins & two sites (3+1))
            for index in range(len(self._421)):
                alpha = self._421[index][0]
                beta = self._421[index][1]
                self._421[index][3] = self._421[index][3] * (
                    self.atoms.spins[alpha] ** 3 * self.atoms.spins[beta]
                )
            # For (four spins & two sites (2+2))
            for index in range(len(self._422)):
                alpha = self._422[index][0]
                beta = self._422[index][1]
                self._422[index][3] = self._422[index][3] * (
                    self.atoms.spins[alpha] ** 2 * self.atoms.spins[beta] ** 2
                )
            # For (four spins & three sites)
            for index in range(len(self._43)):
                alpha = self._43[index][0]
                beta = self._43[index][1]
                gamma = self._43[index][2]
                self._43[index][5] = self._43[index][5] * (
                    self.atoms.spins[alpha] ** 2
                    * self.atoms.spins[beta]
                    * self.atoms.spins[gamma]
                )
            # For (four spins & four sites)
            for index in range(len(self._44)):
                alpha = self._44[index][0]
                beta = self._44[index][1]
                gamma = self._44[index][2]
                epsilon = self._44[index][3]
                self._44[index][7] = self._44[index][7] * (
                    self.atoms.spins[alpha]
                    * self.atoms.spins[beta]
                    * self.atoms.spins[gamma]
                    * self.atoms.spins[epsilon]
                )
        # Before it was normalized
        else:
            # For (one spin & one site)
            for index in range(len(self._1)):
                alpha = self._1[index][0]
                self._1[index][1] = self._1[index][1] / self.atoms.spins[alpha]
            # For (two spins & one site)
            for index in range(len(self._21)):
                alpha = self._21[index][0]
                self._21[index][1] = self._21[index][1] / self.atoms.spins[alpha] ** 2
            # For (two spins & two sites)
            for index in range(len(self._22)):
                alpha = self._22[index][0]
                beta = self._22[index][1]
                self._22[index][3] = self._22[index][3] / (
                    self.atoms.spins[alpha] * self.atoms.spins[beta]
                )
            # For (three spins & one site)
            for index in range(len(self._31)):
                alpha = self._31[index][0]
                self._31[index][1] = self._31[index][1] / self.atoms.spins[alpha] ** 3
            # For (three spins & two sites)
            for index in range(len(self._32)):
                alpha = self._32[index][0]
                beta = self._32[index][1]
                self._32[index][3] = self._32[index][3] / (
                    self.atoms.spins[alpha] ** 2 * self.atoms.spins[beta]
                )
            # For (three spins & three sites)
            for index in range(len(self._33)):
                alpha = self._33[index][0]
                beta = self._33[index][1]
                gamma = self._33[index][2]
                self._33[index][5] = self._33[index][5] / (
                    self.atoms.spins[alpha]
                    * self.atoms.spins[beta]
                    * self.atoms.spins[gamma]
                )
            # For (four spins & one site)
            for index in range(len(self._41)):
                alpha = self._41[index][0]
                self._41[index][1] = self._41[index][1] / self.atoms.spins[alpha] ** 4
            # For (four spins & two sites (3+1))
            for index in range(len(self._421)):
                alpha = self._421[index][0]
                beta = self._421[index][1]
                self._421[index][3] = self._421[index][3] / (
                    self.atoms.spins[alpha] ** 3 * self.atoms.spins[beta]
                )
            # For (four spins & two sites (2+2))
            for index in range(len(self._422)):
                alpha = self._422[index][0]
                beta = self._422[index][1]
                self._422[index][3] = self._422[index][3] / (
                    self.atoms.spins[alpha] ** 2 * self.atoms.spins[beta] ** 2
                )
            # For (four spins & three sites)
            for index in range(len(self._43)):
                alpha = self._43[index][0]
                beta = self._43[index][1]
                gamma = self._43[index][2]
                self._43[index][5] = self._43[index][5] / (
                    self.atoms.spins[alpha] ** 2
                    * self.atoms.spins[beta]
                    * self.atoms.spins[gamma]
                )
            # For (four spins & four sites)
            for index in range(len(self._44)):
                alpha = self._44[index][0]
                beta = self._44[index][1]
                gamma = self._44[index][2]
                epsilon = self._44[index][3]
                self._44[index][7] = self._44[index][7] / (
                    self.atoms.spins[alpha]
                    * self.atoms.spins[beta]
                    * self.atoms.spins[gamma]
                    * self.atoms.spins[epsilon]
                )

    def _set_c1(self, new_c1: float) -> None:
        if new_c1 is None or self.convention._c1 is None:
            return

        new_c1 = float(new_c1)

        if self.convention.c1 == new_c1:
            return

        # If factor is changing one has to scale parameters.
        for index in range(len(self._1)):
            self._1[index][1] = self._1[index][1] * self.convention.c1 / new_c1

    def _set_c21(self, new_c21: float) -> None:
        if new_c21 is None or self.convention._c21 is None:
            return

        new_c21 = float(new_c21)

        if self.convention.c21 == new_c21:
            return

        # If factor is changing one has to scale parameters.
        for index in range(len(self._21)):
            self._21[index][1] = self._21[index][1] * self.convention.c21 / new_c21

    def _set_c22(self, new_c22: float) -> None:
        if new_c22 is None or self.convention._c22 is None:
            return

        new_c22 = float(new_c22)

        if self.convention.c22 == new_c22:
            return

        # If factor is changing one has to scale parameters.
        for index in range(len(self._22)):
            self._22[index][3] = self._22[index][3] * self.convention.c22 / new_c22

    def _set_c31(self, new_c31: float) -> None:
        if new_c31 is None or self.convention._c31 is None:
            return

        new_c31 = float(new_c31)

        if self.convention.c31 == new_c31:
            return

        # If factor is changing one has to scale parameters.
        for index in range(len(self._31)):
            self._31[index][1] = self._31[index][1] * self.convention.c31 / new_c31

    def _set_c32(self, new_c32: float) -> None:
        if new_c32 is None or self.convention._c32 is None:
            return

        new_c32 = float(new_c32)

        if self.convention.c32 == new_c32:
            return

        # If factor is changing one has to scale parameters.
        for index in range(len(self._32)):
            self._32[index][3] = self._32[index][3] * self.convention.c32 / new_c32

    def _set_c33(self, new_c33: float) -> None:
        if new_c33 is None or self.convention._c33 is None:
            return

        new_c33 = float(new_c33)

        if self.convention.c33 == new_c33:
            return

        # If factor is changing one has to scale parameters.
        for index in range(len(self._33)):
            self._33[index][5] = self._33[index][5] * self.convention.c33 / new_c33

    def _set_c41(self, new_c41: float) -> None:
        if new_c41 is None or self.convention._c41 is None:
            return

        new_c41 = float(new_c41)

        if self.convention.c41 == new_c41:
            return

        # If factor is changing one has to scale parameters.
        for index in range(len(self._41)):
            self._41[index][1] = self._41[index][1] * self.convention.c41 / new_c41

    def _set_c421(self, new_c421: float) -> None:
        if new_c421 is None or self.convention._c421 is None:
            return

        new_c421 = float(new_c421)

        if self.convention.c421 == new_c421:
            return

        # If factor is changing one has to scale parameters.
        for index in range(len(self._421)):
            self._421[index][3] = self._421[index][3] * self.convention.c421 / new_c421

    def _set_c422(self, new_c422: float) -> None:
        if new_c422 is None or self.convention._c422 is None:
            return

        new_c422 = float(new_c422)

        if self.convention.c422 == new_c422:
            return

        # If factor is changing one has to scale parameters.
        for index in range(len(self._422)):
            self._422[index][3] = self._422[index][3] * self.convention.c422 / new_c422

    def _set_c43(self, new_c43: float) -> None:
        if new_c43 is None or self.convention._c43 is None:
            return

        new_c43 = float(new_c43)

        if self.convention.c43 == new_c43:
            return

        # If factor is changing one has to scale parameters.
        for index in range(len(self._43)):
            self._43[index][5] = self._43[index][5] * self.convention.c43 / new_c43

    def _set_c44(self, new_c44: float) -> None:
        if new_c44 is None or self.convention._c44 is None:
            return

        new_c44 = float(new_c44)

        if self.convention.c44 == new_c44:
            return

        # If factor is changing one has to scale parameters.
        for index in range(len(self._44)):
            self._44[index][7] = self._44[index][7] * self.convention.c44 / new_c44

    ############################################################################
    #                          External magnetic field                         #
    ############################################################################
    def add_magnetic_field(self, h, alphas=None) -> None:
        r"""
        Adds external magnetic field to the Hamiltonian in the form of one spin
        parameters.

        .. math::

            \mu_B  g_{\alpha} \boldsymbol{h}\cdot\boldsymbol{S}_{\mu,\alpha}
            =
            C_1
            \boldsymbol{S}_{\mu,\alpha}
            \cdot
            \boldsymbol{J}_{Zeeman}(\boldsymbol{r}_{\alpha})

        where :math:`\boldsymbol{J}_{Zeeman}(\boldsymbol{r}_{\alpha})` is defined as

        .. math::

            \boldsymbol{J}_{Zeeman}(\boldsymbol{r}_{\alpha})
            =
            \dfrac{\mu_B g_{\alpha}}{C_1}\boldsymbol{h}

        Parameters
        ----------
        h : (3, ) |array-like|_
            Vector of magnetic field given in the units of Tesla.
        alphas : list of int, optional
            Indices of atoms, to which the magnetic field effect should be added.

        Notes
        -----
        To minimize the energy the magnetic moment will be aligned with the
        direction of the external field. But spin vector will be directed opposite to the
        direction of the magnetic field.

        * If ``alphas is None``, then parameters of the magnetic field added
          only to the magnetic atoms. In other words only to atoms that already have
          at least one other parameter (any) associated with it.
        * If ``alpha is not None``, then parameters of magnetic field are added
          to the atoms with the provided indices (based on the order in
          :py:attr:`.SpinHamiltonian.atoms`)
        """

        if self.convention._c1 is None:
            self.convention._c1 = 1.0

        h = np.array(h, dtype=float)

        BOHR_MAGNETON = 0.057883818060  # meV / Tesla

        if alphas is None:
            alphas = self.map_to_all

        zeeman_parameters = [
            BOHR_MAGNETON * self.atoms.g_factors[alpha] * h / self.convention.c1
            for alpha in alphas
        ]

        i = 0
        j = 0
        new_p1 = []
        while i < len(alphas) or j < len(self._1):
            if i >= len(alphas) or (j < len(self._1) and alphas[i] > self._1[j][0]):
                new_p1.append(self._1[j])
                j += 1
            elif j >= len(self._1) or (i < len(alphas) and alphas[i] < self._1[j][0]):
                new_p1.append([alphas[i], zeeman_parameters[i]])
                i += 1
            elif alphas[i] == self._1[j][0]:
                new_p1.append(
                    [
                        alphas[i],
                        zeeman_parameters[i] + self._1[j][1],
                    ]
                )
                i += 1
                j += 1

        self._1 = new_p1
        self._reset_internals()

    ############################################################################
    #                    Magnetic dipole-dipole interaction                    #
    ############################################################################

    def add_dipole_dipole(self, R_cut=None, E_cut=None, alphas=None):
        r"""
        Add magnetic dipole dipole interaction to the Hamiltonian.

        Magnetic dipole dipole interaction is added in the form of two-spin & two-sites
        parameter

        .. math::

            C_{2,2}
            \sum_{\mu,\nu,\alpha,\beta,i,j}
            J_{dd}(\boldsymbol{r}_{\nu,\alpha\beta})^{ij}
            S_{\mu,\alpha}^i
            S_{\mu+\nu,\beta}^j

        where the parameter is defined as

        .. math::

            J_{dd}(\boldsymbol{r}_{\nu,\alpha\beta})^{ij}
            =
            \dfrac{\mu_0\mu_B^2}{8\pi C_{2,2}}
            \dfrac{g_{\alpha}g_{\beta}}{\vert\boldsymbol{r}_{\nu,\alpha\beta}\vert^3}
            (\delta_{k,l} - 3\hat{r}_{\nu,\alpha\beta}^i\hat{r}_{\nu,\alpha\beta}^j)

        if :py:attr:`.SpinHamiltonian.convention.multiple_counting` is ``True`` and as

        .. math::

            J_{dd}(\boldsymbol{r}_{\nu,\alpha\beta})^{ij}
            =
            \dfrac{\mu_0\mu_B^2}{4\pi C_{2,2}}
            \dfrac{g_{\alpha}g_{\beta}}{\vert\boldsymbol{r}_{\nu,\alpha\beta}\vert^3}
            (\delta_{k,l} - 3\hat{r}_{\nu,\alpha\beta}^i\hat{r}_{\nu,\alpha\beta}^j)

        if :py:attr:`.SpinHamiltonian.convention.multiple_counting` is ``False``.

        where :math:`g_{\alpha}` is a g-factor, :math:`\boldsymbol{\hat{r}}_{\nu,\alpha\beta}`
        is a unit vector.

        Parameters
        ----------
        R_cut : float, optional
            Cut off radius for the distance between a pair of atoms.
            :math:`R_{cut} \ge 0`.
        E_cut : float, optional
            Cut off value for the maximum value of the parameter.
            :math:`E_{cut} > 0`.
        alphas : list of int, optional
            Indices of atoms, to which the magnetic field effect should be added.

        Raises
        ------
        ValueError
            * If none of the  ``R_cut`` or ``E_cut`` are provided.
            * If ``R_cut < 0``
            * If ``E_cut <= 0``

        Notes
        -----

        *   If only ``R_cut`` is given, then the dipole dipole term between the pair of
            spins :math:`S_{\mu,\alpha}^i` and :math:`S_{\mu+\nu,\beta}^j` is added if
            :math:`\vert\boldsymbol{r}_{\nu,\alpha\beta}\vert <= R_{cut}`.

        *   If only ``E_cut`` is given, then the ``R_cut`` is estimated as

            .. math::

                R_{cut}
                =
                \left(
                3\sqrt{2}
                \dfrac{\mu_0\mu_B^2g_{max}^2}{8\pi C_{2,2}E_{cut}}
                \right)^{\dfrac{1}{3}}

            if :py:attr:`.SpinHamiltonian.convention.multiple_counting` is ``True`` and as

            .. math::

                R_{cut}
                =
                \left(
                3\sqrt{2}
                \dfrac{\mu_0\mu_B^2g_{max}^2}{4\pi C_{2,2}E_{cut}}
                \right)^{\dfrac{1}{3}}

            if :py:attr:`.SpinHamiltonian.convention.multiple_counting` is ``False``.

            The dipole dipole term between the pair of spins :math:`S_{\mu,\alpha}^i` and
            :math:`S_{\mu+\nu,\beta}^j` is added if
            :math:`\vert\boldsymbol{r}_{\nu,\alpha\beta}\vert \le R_{cut}` and
            :math:`\vert J_{dd}(\boldsymbol{r}_{\nu,\alpha\beta})^{ij}\vert\ge E_{cut}`
            for some :math:`i, j`.

        *   If both ``R_cut`` and ``E_cut`` are provided, then the dipole dipole term
            between the pair of spins :math:`S_{\mu,\alpha}^i` and
            :math:`S_{\mu+\nu,\beta}^j` is added if
            :math:`\vert\boldsymbol{r}_{\nu,\alpha\beta}\vert \le R_{cut}` and
            :math:`\vert J_{dd}(\boldsymbol{r}_{\nu,\alpha\beta})^{ij}\vert \ge E_{cut}`
            for some :math:`i, j`.

        Magnetic dipole-dipole interaction is added either to magnetic atoms or
        to the list of the atoms provided by user.

        * If ``alphas is None``, then parameters of the magnetic field added
          only to the magnetic atoms. In other words only to atoms that already have
          at least one other parameter (any) associated with it.
        * If ``alpha is not None``, then parameters of magnetic field are added
          to the atoms with the provided indices (based on the order in
          :py:attr:`.SpinHamiltonian.atoms`)
        """
        # Constants
        MU_0_MU_B = (
            1.256637061 * 9.2740100657**2 * 6.241509074 / 1000
        )  # meV * Angstrom^3

        if E_cut is None and R_cut is None:
            raise ValueError("Expected either E_cut or R_cut, got neither.")

        if E_cut is not None:
            if E_cut <= 0:
                raise ValueError(f"Expected positive cut-off energy, got {E_cut}.")

            R_cut = (
                3
                * np.sqrt(2)
                * MU_0_MU_B
                * max(self.atoms.g_factors) ** 2
                / 4
                / np.pi
                / self.convention.c22
                / E_cut
            )

            if self.convention.multiple_counting:
                R_cut = R_cut / 2

            R_cut = R_cut ** (1 / 3)
        else:
            R_cut = float(R_cut)

        if R_cut < 0:
            raise ValueError(f"Expected positive cut-off radius, got {R_cut}.")

        # Get indices for unit cells of interest
        a1, a2, a3 = self.cell
        a_3_perp = abs(np.cross(a1, a2) @ a3 / np.linalg.norm(np.cross(a1, a2)))
        m_3_max = ceil(R_cut / a_3_perp)
        a_2_perp = np.cross(a2, a3) @ a1 / np.linalg.norm(np.cross(a2, a3))
        m_2_max = ceil(R_cut / a_2_perp)

        m_1_max = ceil(R_cut / np.linalg.norm(a1))

        # Run over all pairs of atoms between (0, 0, 0) and all unit cells of
        # interest
        tmp_parameters = []

        if alphas is None:
            alphas = self.map_to_all

        for alpha in alphas:
            for k in range(-m_3_max, m_3_max + 1):
                for j in range(-m_2_max, m_2_max + 1):
                    for i in range(-m_1_max, m_1_max + 1):
                        for beta in alphas:
                            if k == 0 and j == 0 and i == 0 and alpha == beta:
                                continue

                            if _get_primary_p22(
                                alpha=alpha, beta=beta, nu=(i, j, k)
                            ) != (alpha, beta, (i, j, k)):
                                continue

                            vector = (
                                np.array([i, j, k])
                                + self.atoms.positions[beta]
                                - self.atoms.positions[alpha]
                            ) @ self.cell
                            distance = np.linalg.norm(vector)

                            if distance <= R_cut:
                                parameter = (
                                    MU_0_MU_B
                                    / 4
                                    / np.pi
                                    / self.convention.c22
                                    * self.atoms.g_factors[alpha]
                                    * self.atoms.g_factors[beta]
                                    / distance**3
                                ) * (
                                    np.eye(3, dtype=float)
                                    - 3 * np.outer(vector, vector) / distance**2
                                )
                                if self.convention.multiple_counting:
                                    parameter = parameter / 2

                                if E_cut is None or (np.abs(parameter) >= E_cut).any():
                                    tmp_parameters.append(
                                        [alpha, beta, (i, j, k), parameter]
                                    )

        if len(tmp_parameters) > 0:
            tmp_parameters.sort(key=lambda x: x[:-1])

            self._22 = _merge(list1=self._22, list2=tmp_parameters)

    ############################################################################
    #                                Copy getter                               #
    ############################################################################
    def copy(self):
        R"""
        Returns a new, independent copy of the same Hamiltonian.

        Returns
        -------
        spinham : :py:class:`.SpinHamiltonian`
            A new instance of the same Hamiltonian.
        """

        return deepcopy(self)

    def get_empty(self):
        r"""
        Returns the Hamiltonian with the same cell, atoms and convention, but with no
        parameters present.

        Returns
        -------
        spinham : py:class:`.SpinHamiltonian`
            New instance of the spin Hamiltonian.

        Notes
        -----
        Note that in the new Hamiltonian ``spinham.M == 0`` - as there is no parameters
        present, then no atoms are considered to be magnetic.
        """

        return SpinHamiltonian(
            cell=self.cell, atoms=self.atoms, convention=self.convention
        )

    ############################################################################
    #                           Arithmetic operations                          #
    ############################################################################
    def __mul__(self, number):
        if not isinstance(number, int) and not isinstance(number, float):
            raise TypeError(
                f"unsupported operand type(s) for *: '{type(number)}' and 'SpinHamiltonian'"
            )

        spinham = self.copy()

        # One spin
        for i in range(len(spinham._1)):
            spinham._1[i][1] *= number

        # Two spins
        for i in range(len(spinham._21)):
            spinham._21[i][1] *= number

        for i in range(len(spinham._22)):
            spinham._22[i][3] *= number

        # Three spins
        for i in range(len(spinham._31)):
            spinham._31[i][1] *= number

        for i in range(len(spinham._32)):
            spinham._32[i][3] *= number

        for i in range(len(spinham._33)):
            spinham._33[i][5] *= number

        # Four spins
        for i in range(len(spinham._41)):
            spinham._41[i][1] *= number

        for i in range(len(spinham._421)):
            spinham._421[i][3] *= number

        for i in range(len(spinham._422)):
            spinham._422[i][3] *= number

        for i in range(len(spinham._43)):
            spinham._43[i][5] *= number

        for i in range(len(spinham._44)):
            spinham._44[i][7] *= number

        return spinham

    def __rmul__(self, number):
        return self.__mul__(number=number)

    def __add__(self, other):
        if not isinstance(other, SpinHamiltonian):
            raise NotImplementedError

        # Check that unit cells are the same
        if not np.allclose(self.cell, other.cell):
            raise ValueError(
                "Unit cells of two Hamiltonians are different, "
                "summation is not supported"
            )

        # Check that atoms are the same
        same_atoms = True
        if len(self.atoms.names) != len(other.atoms.names):
            same_atoms = False
        else:
            for i in range(len(self.atoms.names)):
                if (
                    self.atoms.names[i] != other.atoms.names[i]
                    or not np.allclose(
                        self.atoms.positions[i], other.atoms.positions[i]
                    )
                    or abs(self.atoms.spins[i] - other.atoms.spins[i]) > 1e-8
                    or abs(self.atoms.g_factors[i] - other.atoms.g_factors[i]) > 1e-8
                ):
                    same_atoms = False

        if not same_atoms:
            raise ValueError(
                "Atoms of two spin Hamiltonians are different, "
                "summation is not supported."
            )

        # Make sure that conventions are the same
        other_convention = other.convention
        other.convention = self.convention

        result = self.get_empty()

        # One spin terms
        result._1 = _merge(list1=self._1, list2=other._1)

        # Two spin terms
        result._21 = _merge(list1=self._21, list2=other._21)
        result._22 = _merge(list1=self._22, list2=other._22)

        # Three spin terms
        result._31 = _merge(list1=self._31, list2=other._31)
        result._32 = _merge(list1=self._32, list2=other._32)
        result._33 = _merge(list1=self._33, list2=other._33)

        # Four spin terms
        result._41 = _merge(list1=self._41, list2=other._41)
        result._421 = _merge(list1=self._421, list2=other._421)
        result._422 = _merge(list1=self._422, list2=other._422)
        result._43 = _merge(list1=self._43, list2=other._43)
        result._44 = _merge(list1=self._44, list2=other._44)

        # Restore convention of other Hamiltonian
        other.convention = other_convention

        return result

    def __sub__(self, other):
        return self + (-1) * other

    ############################################################################
    #                            One spin & one site                           #
    ############################################################################
    p1 = _p1
    add_1 = _add_1
    remove_1 = _remove_1

    ############################################################################
    #                           Two spins & one site                           #
    ############################################################################
    p21 = _p21
    add_21 = _add_21
    remove_21 = _remove_21

    ############################################################################
    #                           Two spins & two sites                          #
    ############################################################################
    p22 = _p22
    add_22 = _add_22
    remove_22 = _remove_22

    ############################################################################
    #                          Three spins & one site                          #
    ############################################################################
    p31 = _p31
    add_31 = _add_31
    remove_31 = _remove_31

    ############################################################################
    #                          Three spins & two sites                         #
    ############################################################################
    p32 = _p32
    add_32 = _add_32
    remove_32 = _remove_32

    ############################################################################
    #                         Three spins & three sites                        #
    ############################################################################
    p33 = _p33
    add_33 = _add_33
    remove_33 = _remove_33

    ############################################################################
    #                           Four spins & one site                          #
    ############################################################################
    p41 = _p41
    add_41 = _add_41
    remove_41 = _remove_41

    ############################################################################
    #                          Four spins & two sites (3+1)                    #
    ############################################################################
    p421 = _p421
    add_421 = _add_421
    remove_421 = _remove_421

    ############################################################################
    #                          Four spins & two sites (2+2)                    #
    ############################################################################
    p422 = _p422
    add_422 = _add_422
    remove_422 = _remove_422

    ############################################################################
    #                         Four spins & three sites                         #
    ############################################################################
    p43 = _p43
    add_43 = _add_43
    remove_43 = _remove_43

    ############################################################################
    #                          Four spins & four sites                         #
    ############################################################################
    p44 = _p44
    add_44 = _add_44
    remove_44 = _remove_44


# Populate __all__ with objects defined in this file
__all__ = list(set(dir()) - old_dir)
# Remove all semi-private objects
__all__ = [i for i in __all__ if not i.startswith("_")]
del old_dir
