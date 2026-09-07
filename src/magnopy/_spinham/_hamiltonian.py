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


from copy import deepcopy
from math import ceil

import numpy as np
from wulfric import add_sugar

from magnopy._spinham._legacy_code import (
    _add_1,
    _add_21,
    _add_22,
    _add_31,
    _add_41,
    _remove_1,
    _remove_21,
    _remove_22,
    _remove_31,
    _remove_41,
)


from magnopy._spinham._convention import Convention
from magnopy._data_validation import (
    _validate_atom_index,
    _validate_unit_cell_index,
    _validated_units,
)
from magnopy._parameters._interaction_parameters import (
    _InteractionParameters,
    _InteractionParametersIterator,
    _get_specs,
)
from magnopy._parameters._equivalent_sets import (
    _get_equivalent,
    _get_missing_parameters,
    _set_distribution,
)


from magnopy._constants._units import _PARAMETER_UNITS, _PARAMETER_UNITS_MAKEUP
from magnopy._constants._si import BOHR_MAGNETON, ANGSTROM, VACUUM_MAGNETIC_PERMEABILITY

# Save local scope at this moment
old_dir = set(dir())
old_dir.add("old_dir")


class SpinHamiltonian:
    r"""
    Spin Hamiltonian.

    Implements the :ref:`user-guide_theory-behind_spin-hamiltonian`.

    Parameters
    ----------

    convention : :py:class:`.Convention` or str
        A convention of the spin Hamiltonian. See :ref:`user-guide_usage_convention` for
        more details.

    cell : (3, 3) |array-like|_
        Matrix of a cell, rows are interpreted as vectors. See
        :ref:`user-guide_usage_cell` for more details.

    atoms : dict
        Dictionary with atoms. See :ref:`user-guide_usage_atoms` for more details.

    units : str, default "meV"
        .. versionadded:: 0.3.0

        Units of the Hamiltonian's parameters. See :py:attr:`.SpinHamiltonian.units`
        for more details. Case-insensitive.

    Examples
    --------

    See :ref:`user-guide_usage_spin-hamiltonian` for examples.

    """

    def __init__(self, cell, atoms, convention, units="meV") -> None:
        self._cell = np.array(cell)

        self._atoms = add_sugar(atoms)

        # Only the magnetic sites
        self._magnetic_atoms = None
        self._map_to_magnetic = None
        self._map_to_all = None

        # To keep track of the external magnetic field
        self._magnetic_field = np.zeros(3, dtype=float)
        self._zeeman_parameters = _InteractionParameters()

        self._convention = convention

        if units.lower() not in _PARAMETER_UNITS:
            raise ValueError(
                f'Given units ("{units}") are not supported. Please use one of\n  * '
                + "\n  * ".join(list(_PARAMETER_UNITS))
            )

        self._units = units.lower()

        self._parameters = _InteractionParameters()

    ############################################################################
    #                              Cell and Atoms                              #
    ############################################################################

    @property
    def cell(self):
        r"""
        Cell of the crystal on which the Hamiltonian is built.

        See :ref:`user-guide_usage_cell` for more details.

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
            AttributeError: Change of the cell attribute is not allowed after the creation of SpinHamiltonian instance. SpinHamiltonian.cell is immutable.
        """

        return self._cell

    @cell.setter
    def cell(self, new_value):
        raise AttributeError(
            "Change of the cell attribute is not allowed after the creation of SpinHamiltonian instance. SpinHamiltonian.cell is immutable."
        )

    @property
    def atoms(self) -> dict:
        r"""
        Atoms of the crystal on which the Hamiltonian is build.

        See :ref:`user-guide_usage_atoms` for more details.

        Returns
        -------

        atoms : dict (with added sugar)
            Dictionary with the atoms.

        See Also
        --------
        M_prime

        Notes
        -----

        If is not recommended to change the atoms property after the creation of
        :py:class:`.SpinHamiltonian`. In fact an attempt to do so raises an
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
            AttributeError: Change of the atoms dictionary is not allowed after the creation of SpinHamiltonian instance. SpinHamiltonian.atoms is immutable.

        """

        return self._atoms

    @atoms.setter
    def atoms(self, new_value):
        raise AttributeError(
            "Change of the atoms dictionary is not allowed after the creation of SpinHamiltonian instance. SpinHamiltonian.atoms is immutable."
        )

    @property
    def M_prime(self):
        r"""
        Amount of atoms in the unit cell.

        .. versionadded:: 0.5.2

        Both magnetic and non-magnetic atoms are counted.

        Returns
        -------

        M_prime : int
            Amount of atoms (magnetic and non-magnetic) in the unit cell.

        See Also
        --------
        atoms

        """

        return len(self.atoms.names)

    def _reset_internals(self):
        self._map_to_magnetic = None
        self._map_to_all = None
        self._magnetic_atoms = None

    def _update_internals(self):
        # Identify magnetic sites
        indices = set()

        for specs, _ in self._parameters._container:
            for alpha in specs[3]:
                indices.add(alpha)

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

        map_to_magnetic (M_prime, ) list of int
            Index map. Integers. ``M_prime = len(spinham.atoms.names)``

        See also
        --------
        map_to_all

        Examples
        --------

        See :ref:`user-guide_usage_spin-hamiltonian_magnetic-vs-non-magnetic` for more
        details on how to use this property.
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

        See also
        --------
        map_to_magnetic

        Examples
        --------

        See :ref:`user-guide_usage_spin-hamiltonian_magnetic-vs-non-magnetic` for more
        details on how to use this property.
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

        magnetic_atoms : dict
            Dictionary of magnetic atoms. Have the same structure as
            :py:attr:`.SpinHamiltonian.atoms`.

        See Also
        --------

        M

        Examples
        --------

        See :ref:`user-guide_usage_spin-hamiltonian_magnetic-vs-non-magnetic` for more
        details on how to use this property.
        """

        if self._magnetic_atoms is None:
            self._update_internals()

        return self._magnetic_atoms

    @property
    def M(self):
        r"""
        Amount of magnetic atoms in the unit cell.

        Returns
        -------

        M : int
            Amount of magnetic atoms in the unit cell.

        See Also
        --------
        magnetic_atoms
        """

        return len(self.magnetic_atoms.names)

    ############################################################################
    #                        Distribution in equal sets                        #
    ############################################################################

    def restore_missing_parameters(self, strategy="zeros") -> None:
        r"""
        Checks that all interactions from the equivalent sets are present in the
        Hamiltonian and adds the missing ones (if any).

        .. versionadded:: 0.5.2

        See :ref:`user-guide_theory-behind_equivalent-parameters` for more details on sets
        of equivalent parameters.

        Parameters
        ----------

        strategy : str, default "zeros"

            When some parameters from the equivalent set are missing, this argument
            defines the value of the added parameters. Case-insensitive. Supported options
            are

            * "zeros" (default)
              Values of missing parameters are set to zeros. This option does not change
              the physics of the Hamiltonian.

            * "mean"
              Values of the missing parameters are set to the mean value of the parameters
              from the equivalent set. This value might change the physics of the
              Hamiltonian (as it adds new non-zero parameters).

        Raises
        ------
        ValueError
            If ``strategy`` is not one of the supported options.

        Notes
        -----
        If ``spinham.convention.multiple_counting`` is False, then this function does
        nothing as only one parameter per set is stored and no parameters can be missing.
        """

        if not self.convention._multiple_counting:
            return

        missing_parameters = _get_missing_parameters(
            parameters=self._parameters, strategy=strategy.lower()
        )

        if len(missing_parameters) > 0:
            self._parameters = self._parameters + missing_parameters

    # DEPRECATED 0.6.0
    # REMOVE in November of 2026
    def symmetrize(self) -> None:
        r"""
        Symmetrize interaction parameters as specified in the SI note 3 of |paper-2026|_.

        .. deprecated:: 0.6.0 Use :py:meth:`.SpinHamiltonian.set_distribution` with ``strategy="symmetrize"`` instead. This method will be removed in November of 2026.


        Please use :py:meth:`.SpinHamiltonian.set_distribution` instead.
        """

        import warnings

        warnings.warn(
            "The method SpinHamiltonian.symmetrize is deprecated since version 0.6.0 and will be removed in November of 2026. Please use SpinHamiltonian.set_distribution with strategy='symmetrize' instead.",
            DeprecationWarning,
            stacklevel=2,
        )

        self.set_distribution(strategy="symmetrized")

    def set_distribution(self, strategy="symmetrize") -> None:
        """
        Enforces one of the supported distributions of parameters within the sets of
        equivalent parameters.

        .. versionadded:: 0.5.2

        See :ref:`user-guide_theory-behind_equivalent-parameters` for more details.

        Parameters
        ----------
        strategy : str, default "symmetrize"

            Strategy for distributing overall value between parameters of the equivalent
            set. Case-insensitive. Supported options are

            * "symmetrize" (default)
              All parameters from the set are equal to each other.

            * "one-for-all"
              The sum is assigned to the single representative parameter. All other
              parameters from the set are set to zeros.

        Raises
        ------
        ValueError
            If ``strategy`` is not one of the supported options.
        ValueError
            If ``spinham.convention.multiple_counting`` is False. See Notes below.

        See Also
        --------
        restore_missing_parameters

        Notes
        -----

        If ``spinham.convention.multiple_counting`` is False, then this function raises
        a ``ValueError`` as the distribution of parameters is fixed and can not be
        changed.

        On contrary to :py:meth:`.SpinHamiltonian.restore_missing_parameters`, this
        method never changes the physics of the Hamiltonian.

        Examples
        --------

        See :ref:`user-guide_usage_spin-hamiltonian_symmetrization` for more details on
        how to use this function.
        """

        # Note: self._zeeman_parameters never change, as there is always one
        # parameter per set

        if not self.convention._multiple_counting:
            raise ValueError(
                "When spinham.convention.multiple_counting is False, the distribution of parameters is fixed and cannot be changed."
            )

        new_parameters = _set_distribution(
            parameters=self._parameters, strategy=strategy.lower()
        )

        self._parameters = new_parameters

    ############################################################################
    #                                Convention                                #
    ############################################################################

    @property
    def convention(self) -> Convention:
        r"""
        Convention of the spin Hamiltonian.

        Returns
        -------

        convention : :py:class:`.Convention`

        See Also
        --------

        Convention

        Examples
        --------

        See :ref:`user-guide_usage_spin-hamiltonian_convention` for more details on how to
        use this property.

        """

        return self._convention

    @convention.setter
    def convention(self, new_convention: Convention):
        self._set_multiple_counting(new_convention._multiple_counting)

        self._set_spin_normalization(new_convention._spin_normalized)

        self._set_convention_constant(
            old_value=self.convention._c1, new_value=new_convention._c1, n=1, p_n=1
        )
        self._set_convention_constant(
            old_value=self.convention._c21, new_value=new_convention._c21, n=2, p_n=1
        )
        self._set_convention_constant(
            old_value=self.convention._c22, new_value=new_convention._c22, n=2, p_n=2
        )
        self._set_convention_constant(
            old_value=self.convention._c31, new_value=new_convention._c31, n=3, p_n=1
        )
        self._set_convention_constant(
            old_value=self.convention._c32, new_value=new_convention._c32, n=3, p_n=2
        )
        self._set_convention_constant(
            old_value=self.convention._c33, new_value=new_convention._c33, n=3, p_n=3
        )
        self._set_convention_constant(
            old_value=self.convention._c41, new_value=new_convention._c41, n=4, p_n=1
        )
        self._set_convention_constant(
            old_value=self.convention._c42, new_value=new_convention._c42, n=4, p_n=2
        )
        self._set_convention_constant(
            old_value=self.convention._c43, new_value=new_convention._c43, n=4, p_n=3
        )
        self._set_convention_constant(
            old_value=self.convention._c44, new_value=new_convention._c44, n=4, p_n=4
        )
        self._set_convention_constant(
            old_value=self.convention._c45, new_value=new_convention._c45, n=4, p_n=5
        )

        self._convention = new_convention

    def _set_multiple_counting(self, multiple_counting: bool) -> None:
        if multiple_counting is None or self.convention._multiple_counting is None:
            return

        multiple_counting = bool(multiple_counting)

        if self.convention.multiple_counting == multiple_counting:
            return

        new_parameters = _InteractionParameters()
        for (n, p_n, nus, alphas), parameter in self._parameters._container:
            equivalent_parameters = _get_equivalent(
                n=n, p_n=p_n, nus=nus, alphas=alphas, parameter=parameter
            )

            degeneracy = len(equivalent_parameters)

            # It was absent before
            if multiple_counting:
                for nus, alphas, parameter in equivalent_parameters:
                    new_parameters.add(
                        specs=(n, p_n, nus, alphas), parameter=parameter / degeneracy
                    )
            # It was present before
            else:
                nus, alphas, parameter = equivalent_parameters[0]
                new_parameters.add(
                    specs=(n, p_n, nus, alphas), parameter=parameter, when_present="sum"
                )

        self._parameters = new_parameters

    def _set_spin_normalization(self, spin_normalized: bool) -> None:
        if spin_normalized is None or self.convention._spin_normalized is None:
            return

        spin_normalized = bool(spin_normalized)

        if self.convention.spin_normalized == spin_normalized:
            return

        for i, (specs, _) in enumerate(self._parameters._container):
            factor = 1.0
            for alpha in specs[3]:
                factor = factor * self.atoms.spins[alpha]
            # Before it was not normalized
            if spin_normalized:
                self._parameters._container[i][1] *= factor
            # Before it was normalized
            else:
                self._parameters._container[i][1] /= factor

    def _set_convention_constant(self, old_value, new_value, n, p_n) -> None:
        if new_value is None or old_value is None:
            return

        new_value = float(new_value)

        if old_value == new_value:
            return

        # If factor is changing one has to scale parameters.
        for i, (specs, _) in enumerate(self._parameters._container):
            if specs[:2] == (n, p_n):
                self._parameters._container[i][1] *= old_value / new_value
            elif specs[:2] > (n, p_n):
                break

    ############################################################################
    #                                   Units                                  #
    ############################################################################
    @property
    def units(self) -> str:
        r"""
        Units of the Hamiltonian's interaction parameters.

        .. versionadded:: 0.3.0

        The parameters of the Hamiltonian are stored in some units of energy (or
        energy-like).

        When user adds a parameters to the Hamiltonian (:py:meth:`.SpinHamiltonian.add`)
        the ``parameter`` argument is understood to be given in the units of
        :py:attr:`.SpinHamiltonian.units`.

        By default the Hamiltonian stores and expects parameters in "meV", but the user
        can choose out of the supported units. See
        :ref:`user-guide_usage_units_parameters` for the full list of supported units.

        When Hamiltonian already has some parameters in it, then the change of
        :py:attr:`.SpinHamiltonian.units` converts all parameter to the new units.
        The parameters that the user tries to add afterwards are expected to be in the new
        units.

        Returns
        -------

        units : str

        Notes
        -----

        List of supported units can be fount in :ref:`user-guide_usage_units` page.

        Examples
        --------

        See :ref:`user-guide_usage_spin-hamiltonian_units` for more details on how to use
        this property.
        """

        return _PARAMETER_UNITS_MAKEUP[self._units]

    @units.setter
    def units(self, new_units: str):
        new_units = _validated_units(units=new_units, supported_units=_PARAMETER_UNITS)

        conversion_factor = _PARAMETER_UNITS[self._units] / _PARAMETER_UNITS[new_units]

        for i in range(len(self._parameters._container)):
            self._parameters._container[i][1] *= conversion_factor

        for i in range(len(self._zeeman_parameters._container)):
            self._zeeman_parameters._container[i][1] *= conversion_factor

        self._units = new_units.lower()

    ############################################################################
    #                          External magnetic field                         #
    ############################################################################

    @property
    def magnetic_field(self):
        r"""
        External magnetic field applied to the Hamiltonian.

        Returns
        -------
        B : (3, ) :numpy:`ndarray`
            Vector of magnetic field (magnetic flux density, B) in the units of
            Tesla.

        See Also
        --------
        set_magnetic_field

        Notes
        -----

        See :py:meth:`.SpinHamiltonian.set_magnetic_field` for more details.

        Examples
        --------

        See :ref:`user-guide_usage_spin-hamiltonian_magnetic-field` for more details on
        how to use this property.
        """

        if len(self._zeeman_parameters) == 0:
            self._magnetic_field = np.zeros(3, dtype=float)

        return self._magnetic_field

    @magnetic_field.setter
    def magnetic_field(self, new_value):
        self.set_magnetic_field(B=new_value)

    def set_magnetic_field(self, B, alphas=None) -> None:
        r"""
        Sets a uniform external magnetic field applied to the Hamiltonian.

        .. versionadded:: 0.5.2

        Parameters
        ----------

        B : (3, ) |array-like|_
            Vector of magnetic field (magnetic flux density, B) given in the units of
            Tesla.

        alphas : list of int, optional
            Indices of atoms, to which the magnetic field effect should be added. If
            none give, then magnetic field is added only to the
            :py:attr:`.SpinHamiltonian.magnetic_atoms`. See
            :ref:`user-guide_usage_spin-hamiltonian_magnetic-vs-non-magnetic` for more
            details on the difference between magnetic and non-magnetic atoms.

        See Also
        --------

        magnetic_field

        Notes
        -----

        The Zeeman term is defined as

        .. math::
            \sum_{\mu, \alpha}
            \mu_B  g_{\alpha} \boldsymbol{B}\cdot\boldsymbol{S}_{\mu,\alpha}


        We map this term to the :ref:`ug_tb_sh_1-1`, which are written as

        .. math::

            C_1
            \sum_{\mu, \alpha}
            \boldsymbol{S}_{\mu,\alpha}
            \cdot
            \boldsymbol{J}^{Zeeman}_{\alpha}

        where :math:`\boldsymbol{J}^{Zeeman}_{\alpha}` is defined as

        .. math::

            \boldsymbol{J}^{Zeeman}_{\alpha}
            =
            \dfrac{\mu_B g_{\alpha}}{C_1}\boldsymbol{B}

        Zeeman energy is minimal when the magnetic moment is aligned with the direction of
        the external field. Therefore, the spin vector shall be directed opposite to the
        direction of the magnetic field. In other words, in the ground state
        :math:`\ket{0}` of the Zeeman Hamiltonian, the eigenvalue of the spin operator
        :math:`\boldsymbol{S}_{\mu,\alpha}^{||\boldsymbol{B}}` shall be equal to
        :math:`-S_{\alpha}`. Therefore, the g-factors (``spinham.atoms.g_factors``) are
        expected to be positive when the Zeeman term is written as we do in Magnopy.


        Examples
        --------

        See :ref:`user-guide_usage_spin-hamiltonian_magnetic-field` for more details on
        how to use this method.
        """

        B = np.array(B, dtype=float)

        if B.shape != (3,):
            raise ValueError(
                f"Expected magnetic field to be a vector of shape (3,), got {B.shape}."
            )

        if self.convention._c1 is None:
            self.convention._c1 = 1.0

        mu_B = BOHR_MAGNETON / _PARAMETER_UNITS[self._units]  # spinham.units / Tesla

        if alphas is None:
            alphas = self.map_to_all
        else:
            alphas = sorted(alphas)

        self._magnetic_field = B

        zeeman_parameters = [
            mu_B * self.atoms.g_factors[alpha] * B / self.convention.c1
            for alpha in alphas
        ]

        zeeman_term = _InteractionParameters()

        for alpha, parameter in zip(alphas, zeeman_parameters):
            zeeman_term.add(specs=(1, 1, (), (alpha,)), parameter=parameter)

        if len(self._zeeman_parameters) != 0:
            self._parameters = self._parameters - self._zeeman_parameters

        self._zeeman_parameters = zeeman_term

        self._parameters = self._parameters + zeeman_term

        self._reset_internals()

    @property
    def zeeman_parameters(self):
        r"""
        Returns an iterator over Zeeman parameters of the Hamiltonian.

        .. versionadded:: 0.5.2

        Returns
        -------
        parameters : iterator
            Iterator over parameters of the Hamiltonian that originate from the Zeeman
            interaction. See :py:attr:`.SpinHamiltonian.parameters` for more details on
            the returned iterator.

        See Also
        --------
        set_magnetic_field
        magnetic_field
        parameters
        """

        return _InteractionParametersIterator(
            parameters=self._zeeman_parameters, n=1, p_n=1
        )

    def add_magnetic_field(self, B=None, alphas=None) -> None:
        r"""
        Adds external magnetic field to the Hamiltonian.

        Parameters
        ----------

        B : (3, ) |array-like|_
            See :py:attr:`.SpinHamiltonian.set_magnetic_field` for details.

        alphas : list of int, optional
            See :py:attr:`.SpinHamiltonian.set_magnetic_field` for details.

        Notes
        -----

        The call ``spinham.add_magnetic_field(B = B, alphas = alphas)`` is equivalent to
        the call
        ``spinham.set_magnetic_field(B = B + spinham.magnetic_field, alphas = alphas)``.

        In other words, this method "adds" the magnetic field to the Hamiltonian whether
        there was some magnetic field before or not.

        On contrary, :py:attr:`.SpinHamiltonian.set_magnetic_field` "sets" the magnetic
        field, i.e. replaces the previous magnetic field (if any).
        """

        self.set_magnetic_field(B=B + self.magnetic_field, alphas=alphas)

    ############################################################################
    #                    Magnetic dipole-dipole interaction                    #
    ############################################################################

    def add_dipole_dipole(self, R_cut=None, E_cut=None, alphas=None):
        r"""
        Adds magnetic dipole dipole interaction to the Hamiltonian.

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
            :math:`E_{cut} > 0`. Expected in the same units as
            :py:attr:`.SpinHamiltonian.units`.

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
            VACUUM_MAGNETIC_PERMEABILITY
            * BOHR_MAGNETON**2
            / ANGSTROM**3
            / _PARAMETER_UNITS[self._units]
        )  # spinham.units * Angstrom^3

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
        dd_parameters = _InteractionParameters()

        if alphas is None:
            alphas = self.map_to_all

        for alpha in alphas:
            for k in range(-m_3_max, m_3_max + 1):
                for j in range(-m_2_max, m_2_max + 1):
                    for i in range(-m_1_max, m_1_max + 1):
                        for beta in alphas:
                            if ((i, j, k), beta) == ((0, 0, 0), alpha):
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
                                    if not self.convention.multiple_counting:
                                        equivalent_parameters = _get_equivalent(
                                            n=2,
                                            p_n=2,
                                            nus=((i, j, k),),
                                            alphas=(alpha, beta),
                                            parameter=parameter,
                                        )

                                        nu, alphas, parameter = equivalent_parameters[0]
                                    else:
                                        nu = (i, j, k)
                                        alphas = (alpha, beta)

                                    dd_parameters.add(
                                        specs=(2, 2, (nu,), alphas),
                                        parameter=parameter,
                                        when_present="skip",
                                    )

        self._parameters = self._parameters + dd_parameters

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
        Returns the Hamiltonian with the same cell, atoms, units and convention, but with no
        parameters present.

        Returns
        -------

        spinham : :py:class:`.SpinHamiltonian`
            New instance of the spin Hamiltonian.

        Notes
        -----
        Note that in the new Hamiltonian ``spinham.M == 0`` - as there is no parameters
        present, meaning that no atoms are considered to be magnetic.
        """

        return SpinHamiltonian(
            cell=self.cell,
            atoms=self.atoms,
            convention=self.convention,
            units=self.units,
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

        spinham._parameters = spinham._parameters * number

        return spinham

    def __rmul__(self, number):
        return self.__mul__(number=number)

    def __add__(self, other):
        if not isinstance(other, SpinHamiltonian):
            raise NotImplementedError

        # Check that unit cells are the same
        if not np.allclose(self.cell, other.cell):
            raise ValueError(
                "Unit cells of two Hamiltonians are different, summation is not supported"
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
                "Atoms of two spin Hamiltonians are different, summation is not supported."
            )

        result = other.copy()

        # Make sure that units are the same
        result.units = self.units

        # Make sure that conventions are the same
        result.convention = self.convention

        result._parameters = self._parameters + result._parameters

        # Make sure that Zeeman parameters and magnetic field are correctly tracked.
        result._magnetic_field = self._magnetic_field + result._magnetic_field
        result._zeeman_parameters = self._zeeman_parameters + result._zeeman_parameters

        return result

    def __sub__(self, other):
        return self + (-1) * other

    ############################################################################
    #                          Interaction parameters                          #
    ############################################################################
    def add(
        self,
        nus,
        alphas,
        parameter,
        units=None,
        populate_equivalent=False,
        when_present="raise error",
    ):
        r"""
        Adds any parameter with at most four components of spin operator to the
        Hamiltonian.

        .. versionadded:: 0.5.0

        See :ref:`user-guide_theory-behind_spin-hamiltonian` for the definition of the
        Hamiltonian.

        Parameters
        ----------

        nus : (n, 3) or (n-1, 3) list/tuple of tuple of int
            List of unit cell indices associated with the parameter. Each unit cell index
            is a tuple of three integers (t_1, t_2, t_3) corresponding to the translation
            by :math:`t_1 \boldsymbol{a}_1 + t_2 \boldsymbol{a}_2 + t_3 \boldsymbol{a}_3`.

        alphas : (n,) list/tuple of int
            List of indices of atoms associated with the parameter. Based on the order in
            :py:attr:`.SpinHamiltonian.atoms`.

        units : str, optional
            Units in which the  ``parameter``  is given.  Parameters have the the units of
            energy. By default assumes :py:attr:`.SpinHamiltonian.units`.  For the list of
            the  supported  units  see  :ref:`user-guide_usage_units_parameters`.  If
            given  ``units``  are different from  :py:attr:`.SpinHamiltonian.units`,  then
            the  parameter's  value  will be  converted  automatically from  ``units``  to
            :py:attr:`.SpinHamiltonian.units`.

        populate_equivalent : bool, default False
            Whether to automatically populate all equivalent parameters related by the
            symmetrization procedure. Ignored if ``convention.multiple_counting`` is
            ``False``. See :ref:`user-guide_usage_spin-hamiltonian_equivalent-parameters`
            for more details.

        when_present : str, default "raise error"
            Action to take if such parameter is already present in the Hamiltonian.
            Case-insensitive. Supported values are:

            - ``"raise error"`` (default): raises an error.
            - ``"replace"``: replace existing value of the parameter with the new one.
            - ``"sum"``: add the value of the parameter to the existing one.
            - ``"mean"``: replace the value of the parameter with the arithmetic mean of
              existing and new parameters.
            - ``"skip"``: Leave existing parameter unchanged and continue without raising
              an error.

        Notes
        -----

        If ``len(nus) == len(alphas) - 1``, then

        * ``alphas[0]`` always located in the unit cell (0, 0, 0).
        * ``alphas[i]`` is located in the unit cell ``nus[i-1]`` for ``i >= 1``.

        If ``len(nus) == len(alphas)``, then

        * ``alphas[i]`` is located in the unit cell ``nus[i]`` for all ``i``.

        Note that the translational symmetry is always enforced, so in the second case the
        ``nus`` are updated as ``nus[i] = nus[i] - nus[0]`` for all ``i`` before the
        parameter is added to the Hamiltonian.

        Examples
        --------

        See :ref:`user-guide_usage_spin-hamiltonian_adding-parameters` for more details
        on how to use this method.

        """

        self._reset_internals()

        for alpha in alphas:
            _validate_atom_index(index=alpha, atoms=self.atoms)

        for nu in nus:
            _validate_unit_cell_index(ijk=nu)

        if not (len(alphas) - 1 == len(nus) or len(alphas) == len(nus)):
            raise ValueError(
                f"Expected number of unit cell indices to be equal to or one less than the number of atom indices, got {len(alphas)} atom indices and {len(nus)} unit cell indices."
            )

        if len(nus) == len(alphas) - 1:
            nus = tuple([(0, 0, 0)] + list(nus))

        parameter = np.array(parameter)

        if not len(parameter.shape) == len(alphas) or any(
            _ != 3 for _ in parameter.shape
        ):
            raise ValueError(
                f"Expected parameter to be a tensor with {len(alphas)} components of size 3, got {parameter.shape}."
            )

        if units is not None:
            units = _validated_units(units=units, supported_units=_PARAMETER_UNITS)
            parameter = (
                parameter * _PARAMETER_UNITS[units] / _PARAMETER_UNITS[self._units]
            )

        specs = _get_specs(nus=nus, alphas=alphas)

        if self.convention.multiple_counting:
            if populate_equivalent:
                equivalent_parameters = _get_equivalent(
                    n=specs[0],
                    p_n=specs[1],
                    nus=specs[2],
                    alphas=specs[3],
                    parameter=parameter,
                )

                for nus, alphas, parameter in equivalent_parameters:
                    self._parameters.add(
                        specs=(specs[0], specs[1], nus, alphas),
                        parameter=parameter,
                        when_present=when_present,
                    )
            else:
                self._parameters.add(
                    specs=specs, parameter=parameter, when_present=when_present
                )
        else:
            equivalent_parameters = _get_equivalent(
                n=specs[0],
                p_n=specs[1],
                nus=specs[2],
                alphas=specs[3],
                parameter=parameter,
            )
            nus, alphas, parameter = equivalent_parameters[0]

            self._parameters.add(
                specs=(specs[0], specs[1], nus, alphas),
                parameter=parameter,
                when_present=when_present,
            )

    def remove(self, nus, alphas, remove_equivalent=False):
        r"""
        Removes any parameter with at most four components of spin operator from the
        Hamiltonian.

        .. versionadded:: 0.5.0

        See :ref:`user-guide_theory-behind_spin-hamiltonian` for the definition of the
        Hamiltonian.

        Parameters
        ----------

        nus : (n, 3) or (n-1, 3) list/tuple of tuple of int
            List of unit cell indices associated with the parameter. Each unit cell index
            is a tuple of three integers (t_1, t_2, t_3) corresponding to the translation
            by :math:`t_1 \boldsymbol{a}_1 + t_2 \boldsymbol{a}_2 + t_3 \boldsymbol{a}_3`.

        alphas : (n,) list/tuple of int
            List of indices of atoms associated with the parameter. Based on the order in
            :py:attr:`.SpinHamiltonian.atoms`.

        remove_equivalent : bool, default False
            Whether to automatically remove all equivalent parameters related by the
            symmetrization procedure. Ignored if ``convention.multiple_counting`` is
            ``False``. See :ref:`user-guide_usage_spin-hamiltonian_equivalent-parameters`
            for more details.


        Notes
        -----

        See notes of :py:meth:`.SpinHamiltonian.add` for the details on ``nus`` and
        ``alphas``.

        Be careful when removing :py:attr:`.SpinHamiltonian.p1` parameters when
        :py:attr:`.SpinHamiltonian.magnetic_field` is not zero, as the Zeeman term is
        stored as a :py:attr:`.SpinHamiltonian.p1` parameters.

        Examples
        --------

        See :ref:`user-guide_usage_spin-hamiltonian_removing-parameters` for more details
        on how to use this method.

        """

        for alpha in alphas:
            _validate_atom_index(index=alpha, atoms=self.atoms)

        for nu in nus:
            _validate_unit_cell_index(ijk=nu)

        if not (len(alphas) - 1 == len(nus) or len(alphas) == len(nus)):
            raise ValueError(
                f"Expected number of unit cell indices to be equal to or one less than the number of atom indices, got {len(alphas)} atom indices and {len(nus)} unit cell indices."
            )

        if len(nus) == len(alphas) - 1:
            nus = tuple([(0, 0, 0)] + list(nus))

        specs = _get_specs(nus=nus, alphas=alphas)

        if (
            specs[0] == 1
            and specs[1] == 1
            and not np.allclose(self.magnetic_field, np.zeros(3))
        ):
            import warnings

            warnings.warn(
                "Attempt to remove p1 parameter when magnetic field is not zero.",
                UserWarning,
                stacklevel=2,
            )

        # Need to remove from the Zeeman term as well, if it is present.
        # Multiple counting does not affect the Zeeman term.
        if specs[0] == 1 and specs[1] == 1 and len(self._zeeman_parameters) != 0:
            self._zeeman_parameters.remove(specs=specs)

        # Remove from the main container of parameters.
        if self.convention.multiple_counting:
            if remove_equivalent:
                equivalent_parameters = _get_equivalent(
                    n=specs[0],
                    p_n=specs[1],
                    nus=specs[2],
                    alphas=specs[3],
                    parameter=None,
                )

                for nus, alphas, _ in equivalent_parameters:
                    self._parameters.remove(specs=(specs[0], specs[1], nus, alphas))
            else:
                self._parameters.remove(specs=specs)
        else:
            equivalent_parameters = _get_equivalent(
                n=specs[0], p_n=specs[1], nus=specs[2], alphas=specs[3]
            )
            nus, alphas, _ = equivalent_parameters[0]

            self._parameters.remove(specs=(specs[0], specs[1], nus, alphas))

    def parameters(self, n=None, p_n=None):
        r"""
        Returns an iterator over parameters of the Hamiltonian.

        .. versionadded:: 0.5.0

        See :ref:`user-guide_theory-behind_spin-hamiltonian` for the definition of the
        Hamiltonian and the meaning of ``n`` and ``p_n``.

        Parameters
        ----------
        n : int, optional
            Number of spins in the terms of the spin Hamiltonian. Expected to be between
            1 and 4. If not given, then all parameters are returned and ``p_n`` is
            ignored.
        p_n : int, optional
            Index of integer partition of ``n``. Expected to be between 1 and maximal
            number of integer partitions for the given ``n``. If not given, then
            parameters of all respective terms of the spin Hamiltonian with ``n`` spins
            are returned. Ignored if ``n`` is not given.

        Returns
        -------
        parameters : iterator
            Iterator over parameters of the Hamiltonian. Each element of the iterator is
            a tuple ``(nus, alphas, parameter)`` where ``nus`` is a
            tuple of unit cell indices, ``alphas`` is a tuple of atom indices, and
            ``parameter`` is a vector/matrix/tensor of the interaction parameter.

            - ``len(alphas) = n``
            - ``len(nus) == len(alphas) - 1``
            - ``len(parameter.shape) == len(alphas)``
            - ``parameter.shape[i] == 3`` for all ``0 <= i < len(alphas)``
            - ``alphas[0]`` is always located in the unit cell (0, 0, 0)
            - ``alphas[i]`` is located in the unit cell ``nus[i-1]`` for
              ``1 <= i < len(alphas)``

        See Also
        --------
        add
        remove
        p1
        p21
        p22
        p31
        p32
        p33
        p41
        p42
        p43
        p44
        p45
        """

        return _InteractionParametersIterator(parameters=self._parameters, n=n, p_n=p_n)

    def purge(self, tolerance=None):
        r"""
        Removes parameters of the Hamiltonian that are smaller than the given tolerance.

        .. versionadded:: 0.5.2

        Parameters
        ----------
        tolerance : float, default 1e-8 meV
            Parameters with absolute value of all components smaller than the tolerance
            are removed. Expected to be non-negative. Expected to be
            in the same units as :py:attr:`.SpinHamiltonian.units`. Default
            value is :math:`10^{-8}` meV (converted to the units of the Hamiltonian).

        Examples
        --------

        First lets create a Hamiltonian

        .. doctest::

            >>> import numpy as np
            >>> import magnopy
            >>> cell = np.eye(3)
            >>> atoms = dict(
            ...     names=["Fe"], positions=[[0, 0, 0]], spins=[5 / 2], g_factors=[2]
            ... )
            >>> convention = magnopy.Convention(
            ...     multiple_counting=True, spin_normalized=False, c22=1
            ... )
            >>> spinham = magnopy.SpinHamiltonian(
            ...     cell=cell, atoms=atoms, convention=convention
            ... )

        and add some parameters to it

        .. doctest::

            >>> spinham.add(
            ...     nus=[(1, 0, 0)],
            ...     alphas=[0, 0],
            ...     parameter=np.eye(3) * 1e-9,
            ...     populate_equivalent=True,
            ... )
            >>> spinham.add(
            ...     nus=[(0, 1, 0)],
            ...     alphas=[0, 0],
            ...     parameter=np.eye(3) * 1e-7,
            ...     populate_equivalent=True,
            ... )

        There are four interaction parameters in the Hamiltonian

        .. doctest::

            >>> len(spinham.parameters())
            4
            >>> for nus, alphas, parameter in spinham.parameters():
            ...     print(nus, alphas)
            ((-1, 0, 0),) (0, 0)
            ((0, -1, 0),) (0, 0)
            ((0, 1, 0),) (0, 0)
            ((1, 0, 0),) (0, 0)


        Now if we purge the Hamiltonian with the default tolerance of 1e-8 meV, then
        only two parameters remain

        .. doctest::

            >>> spinham.purge()
            >>> len(spinham.parameters())
            2
            >>> for nus, alphas, parameter in spinham.parameters():
            ...     print(nus, alphas)
            ((0, -1, 0),) (0, 0)
            ((0, 1, 0),) (0, 0)
        """

        if tolerance is None:
            tolerance = 1e-8 * _PARAMETER_UNITS["mev"] / _PARAMETER_UNITS[self._units]
        tolerance = abs(float(tolerance))

        # Remove from main container
        new_parameters = _InteractionParameters()
        for specs, parameter in self._parameters._container:
            if not np.all(np.abs(parameter) < tolerance):
                new_parameters._container.append([specs, parameter])
        self._parameters = new_parameters

        # Remove from Zeeman parameters
        new_zeeman_parameters = _InteractionParameters()
        for specs, parameter in self._zeeman_parameters._container:
            if not np.all(np.abs(parameter) < tolerance):
                new_zeeman_parameters._container.append([specs, parameter])
        self._zeeman_parameters = new_zeeman_parameters

    ############################################################################
    #                          One spin & one site (1)                         #
    ############################################################################
    @property
    def p1(self):
        r"""
        Parameters of one spin & one site term of the Hamiltonian (:ref:`ug_tb_sh_1-1`).

        See :ref:`user-guide_theory-behind_spin-hamiltonian` for more details.

        Returns
        -------
        parameters : iterator
            Iterator over parameters of the Hamiltonian. See notes of
            :py:meth:`.SpinHamiltonian.parameters` for more details.

        See Also
        --------
        add
        remove
        parameters
        """

        return self.parameters(n=1, p_n=1)

    ############################################################################
    #                       Two spins & one site (2 + 0)                       #
    ############################################################################

    @property
    def p21(self):
        r"""
        Parameters of two spins & one site term of the Hamiltonian (:ref:`ug_tb_sh_2-1`).

        See :ref:`user-guide_theory-behind_spin-hamiltonian` for more details.

        Returns
        -------
        parameters : iterator
            Iterator over parameters of the Hamiltonian. See notes of
            :py:meth:`.SpinHamiltonian.parameters` for more details.

        See Also
        --------
        add
        remove
        parameters
        """

        return self.parameters(n=2, p_n=1)

    ############################################################################
    #                        Two spins & two sites (1+1)                       #
    ############################################################################

    @property
    def p22(self):
        r"""
        Parameters of two spins & two sites term of the Hamiltonian (:ref:`ug_tb_sh_2-2`).

        See :ref:`user-guide_theory-behind_spin-hamiltonian` for more details.

        Returns
        -------
        parameters : iterator
            Iterator over parameters of the Hamiltonian. See notes of
            :py:meth:`.SpinHamiltonian.parameters` for more details.

        See Also
        --------
        add
        remove
        parameters
        """

        return self.parameters(n=2, p_n=2)

    ############################################################################
    #                      Three spins & one site (3+0+0)                      #
    ############################################################################

    @property
    def p31(self):
        r"""
        Parameters of three spins & one site term of the Hamiltonian
        (:ref:`ug_tb_sh_3-1`).

        See :ref:`user-guide_theory-behind_spin-hamiltonian` for more details.

        Returns
        -------
        parameters : iterator
            Iterator over parameters of the Hamiltonian. See notes of
            :py:meth:`.SpinHamiltonian.parameters` for more details.

        See Also
        --------
        add
        remove
        parameters
        """

        return self.parameters(n=3, p_n=1)

    ############################################################################
    #                      Three spins & two sites (2+1+0)                     #
    ############################################################################
    @property
    def p32(self):
        r"""
        Parameters of three spins & two sites term of the Hamiltonian
        (:ref:`ug_tb_sh_3-2`).

        See :ref:`user-guide_theory-behind_spin-hamiltonian` for more details.

        Returns
        -------
        parameters : iterator
            Iterator over parameters of the Hamiltonian. See notes of
            :py:meth:`.SpinHamiltonian.parameters` for more details.

        See Also
        --------
        add
        remove
        parameters
        """

        return self.parameters(n=3, p_n=2)

    ############################################################################
    #                     Three spins & three sites (1+1+1)                    #
    ############################################################################
    @property
    def p33(self):
        r"""
        Parameters of three spins & three sites term of the Hamiltonian
        (:ref:`ug_tb_sh_3-3`).

        See :ref:`user-guide_theory-behind_spin-hamiltonian` for more details.

        Returns
        -------
        parameters : iterator
            Iterator over parameters of the Hamiltonian. See notes of
            :py:meth:`.SpinHamiltonian.parameters` for more details.

        See Also
        --------
        add
        remove
        parameters
        """

        return self.parameters(n=3, p_n=3)

    ############################################################################
    #                     Four spins & one site (4+0+0+0)                      #
    ############################################################################
    @property
    def p41(self):
        r"""
        Parameters of four spins & one site term of the Hamiltonian (:ref:`ug_tb_sh_4-1`).

        See :ref:`user-guide_theory-behind_spin-hamiltonian` for more details.

        Returns
        -------
        parameters : iterator
            Iterator over parameters of the Hamiltonian. See notes of
            :py:meth:`.SpinHamiltonian.parameters` for more details.

        See Also
        --------
        add
        remove
        parameters
        """

        return self.parameters(n=4, p_n=1)

    ############################################################################
    #                     Four spins & two sites (3+1+0+0)                     #
    ############################################################################
    @property
    def p42(self):
        r"""
        Parameters of four spins & two sites term of the Hamiltonian
        (:ref:`ug_tb_sh_4-2`).

        See :ref:`user-guide_theory-behind_spin-hamiltonian` for more details.

        Returns
        -------
        parameters : iterator
            Iterator over parameters of the Hamiltonian. See notes of
            :py:meth:`.SpinHamiltonian.parameters` for more details.

        See Also
        --------
        add
        remove
        parameters
        """

        return self.parameters(n=4, p_n=2)

    ############################################################################
    #                     Four spins & two sites (2+2+0+0)                     #
    ############################################################################
    @property
    def p43(self):
        r"""
        Parameters of four spins & two sites term of the Hamiltonian
        (:ref:`ug_tb_sh_4-3`).

        See :ref:`user-guide_theory-behind_spin-hamiltonian` for more details.

        Returns
        -------
        parameters : iterator
            Iterator over parameters of the Hamiltonian. See notes of
            :py:meth:`.SpinHamiltonian.parameters` for more details.

        See Also
        --------
        add
        remove
        parameters
        """

        return self.parameters(n=4, p_n=3)

    ############################################################################
    #                    Four spins & three sites (2+1+1+0)                    #
    ############################################################################
    @property
    def p44(self):
        r"""
        Parameters of four spins & three sites term of the Hamiltonian
        (:ref:`ug_tb_sh_4-4`).

        See :ref:`user-guide_theory-behind_spin-hamiltonian` for more details.

        Returns
        -------
        parameters : iterator
            Iterator over parameters of the Hamiltonian. See notes of
            :py:meth:`.SpinHamiltonian.parameters` for more details.

        See Also
        --------
        add
        remove
        parameters
        """

        return self.parameters(n=4, p_n=4)

    ############################################################################
    #                    Four spins & four sites (1+1+1+1)                     #
    ############################################################################
    @property
    def p45(self):
        r"""
        Parameters of four spins & four sites term of the Hamiltonian
        (:ref:`ug_tb_sh_4-5`).

        See :ref:`user-guide_theory-behind_spin-hamiltonian` for more details.

        Returns
        -------
        parameters : iterator
            Iterator over parameters of the Hamiltonian. See notes of
            :py:meth:`.SpinHamiltonian.parameters` for more details.

        See Also
        --------
        add
        remove
        parameters
        """

        return self.parameters(n=4, p_n=5)

    ############################################################################
    #                              Legacy methods                              #
    ############################################################################

    add_1 = _add_1
    remove_1 = _remove_1

    add_21 = _add_21
    remove_21 = _remove_21

    add_22 = _add_22
    remove_22 = _remove_22

    add_31 = _add_31
    remove_31 = _remove_31

    add_41 = _add_41
    remove_41 = _remove_41


# Populate __all__ with objects defined in this file
__all__ = list(set(dir()) - old_dir)
# Remove all semi-private objects
__all__ = [i for i in __all__ if not i.startswith("_")]
del old_dir
