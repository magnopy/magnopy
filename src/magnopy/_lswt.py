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


import numpy as np

from magnopy._diagonalization import solve_via_colpa
from magnopy._exceptions import ColpaFailed
from magnopy._local_rf import span_local_rfs

from magnopy._data_validation import _validated_units
from magnopy._constants._units import _ENERGY_UNITS, _MAGNON_ENERGY_UNITS
from magnopy._parameters._interaction_parameters import _InteractionParametersIterator
from magnopy._parameters._renormalization import _renormalized_parameters


# Save local scope at this moment
old_dir = set(dir())
old_dir.add("old_dir")


class LSWT:
    r"""
    Linear Spin Wave theory.

    It is created from some spin Hamiltonian and set of direction vectors, that defines
    the ground state.

    Parameters
    ----------

    spinham : :py:class:`.SpinHamiltonian`
        Spin Hamiltonian.

    spin_directions : (M, 3) |array-like|_
        Directions of spin vectors. Only directions of vectors are used, modulus is
        ignored. If spin Hamiltonian contains non-magnetic atom, then only the spin
        directions for the magnetic atoms are expected. The order of spin directions is
        the same as the order of magnetic atoms in
        :py:attr:`SpinHamiltonian.magnetic_atoms`. See Notes for more details.

    Attributes
    ----------

    z : (M, 3) :numpy:`ndarray`
        Spin directions (directions of local quantization axes).

    p : (M, 3) :numpy:`ndarray`
        Hybridized x and y components of the local coordinate system
        :math:`\mathbf{p} = \mathbf{x} + i \mathbf{y}`.

    M : int
        Number of spins in the unit cell

    cell : (3, 3) :numpy:`ndarray`
        Unit cell. Rows are vectors, columns are cartesian components.

    spins : (M, ) :numpy:`ndarray`
        Spin values of the magnetic centers.

    Notes
    -----

    Let the spin Hamiltonian contain three atoms Cr1, Br, Cr3 in that order. Assume that
    two atoms are magnetic (Cr1 and Cr3), one atom is not (Br). Then ``spin_directions``
    is a (2, 3) array with ``spin_directions[0]`` being the direction for spin of Cr1 and
    ``spin_directions[1]`` being the direction of spin for Cr3.

    Examples
    --------

    .. doctest::

        >>> import magnopy
        >>> spinham = magnopy.examples.cubic_ferro_nn()
        >>> lswt = magnopy.LSWT(spinham=spinham, spin_directions=[[0, 0, 1]])
    """

    def __init__(self, spinham, spin_directions):
        spin_directions = np.array(spin_directions, dtype=float)
        spin_directions /= np.linalg.norm(spin_directions, axis=1)[:, np.newaxis]

        x, y, self.z = span_local_rfs(
            directional_vectors=spin_directions, hybridize=False
        )
        self.p = x + 1j * y

        self.spins = np.array(spinham.magnetic_atoms.spins, dtype=float)
        self.M = spinham.M
        self.cell = spinham.cell

        initial_units = spinham.units
        initial_convention = spinham.convention

        self._convention = initial_convention.get_modified(
            spin_normalized=False,
            multiple_counting=True,
            c1=1,
            c21=1,
            c22=1,
            c31=1,
            c32=1,
            c33=1,
            c41=1,
            c42=1,
            c43=1,
            c44=1,
            c45=1,
            keep_undefined=True,
        )

        spinham.units = "mev"
        spinham.convention = self._convention
        self._parameters = spinham._parameters.copy()
        spinham.units = initial_units
        spinham.convention = initial_convention

        for i, ((n, p_n, nus, alphas), _) in enumerate(self._parameters._container):
            alphas = tuple(spinham.map_to_magnetic[alpha] for alpha in alphas)
            self._parameters._container[i][0] = (n, p_n, nus, alphas)

        self._parameters = _renormalized_parameters(
            parameters=self._parameters,
            convention=self._convention,
            spin_directions=spin_directions,
            spin_values=self.spins,
        )

        self.A1 = np.zeros((self.M), dtype=float)

        for _, alphas, parameter in _InteractionParametersIterator(
            self._parameters, n=1, p_n=1
        ):
            self.A1[alphas[0]] = 0.5 * parameter @ self.z[alphas[0]]

        self.A2 = {}
        self.B2 = {}
        for nus, alphas, parameter in _InteractionParametersIterator(
            self._parameters, n=2
        ):
            nu = nus[0]
            if nu not in self.A2:
                self.A2[nu] = np.zeros((self.M, self.M), dtype=complex)
                self.B2[nu] = np.zeros((self.M, self.M), dtype=complex)

            self.A2[nu][alphas[0], alphas[1]] = (
                0.5
                * np.sqrt(self.spins[alphas[0]] * self.spins[alphas[1]])
                * (self.p[alphas[0]] @ parameter @ np.conjugate(self.p[alphas[1]]))
            )
            self.B2[nu][alphas[0], alphas[1]] = (
                0.5
                * np.sqrt(self.spins[alphas[0]] * self.spins[alphas[1]])
                * (self.p[alphas[0]] @ parameter @ self.p[alphas[1]])
            )

    def E_2(self, units="meV") -> float:
        r"""
        Computes the correction to the ground state energy that arises from the LSWT.

        Parameters
        ----------

        units : str, default "meV"
            .. versionadded:: 0.3.0

            Units of energy. See :ref:`user-guide_usage_units_energy` for the full
            list of supported units.

        Returns
        -------

        E_2 : float

        Examples
        --------

        .. doctest::

            >>> import magnopy
            >>> spinham = magnopy.examples.cubic_ferro_nn()
            >>> lswt = magnopy.LSWT(spinham=spinham, spin_directions=[[0, 0, 1]])
            >>> # Default units are meV
            >>> lswt.E_2()
            -1.5
        """

        result = float(np.sum(self.A1))

        # Convert units if necessary
        if units != "meV":
            units = _validated_units(units=units, supported_units=_ENERGY_UNITS)
            result = result * _ENERGY_UNITS["mev"] / _ENERGY_UNITS[units]

        return result

    def O(self, units="meV"):  # noqa E743
        r"""
        Computes coefficient of the one-operator terms.

        Parameters
        ----------

        units : str, default "meV"
            .. versionadded:: 0.3.0

            Units of energy. See :ref:`user-guide_usage_units_energy` for the full
            list of supported units.


        Returns
        -------

        O : (M, ) :numpy:`ndarray`
            Elements are complex numbers.

        Notes
        -----

        Before the diagonalization, the magnon Hamiltonian has the form

        .. math::

            \mathcal{H}
            =
            \dots
            +
            \sqrt{N}
            \sum_{\alpha}
            \Bigl(
            O_{\alpha}
            a_{\alpha}(\boldsymbol{0})
            +
            \overline{O_{\alpha}}
            a^{\dagger}_{\alpha}(\boldsymbol{0})
            \Bigr)
            +
            \dots

        where overline denotes complex conjugation. This function computes the
        coefficients :math:`O_{\alpha}`.

        Examples
        --------

        .. doctest::

            >>> import magnopy
            >>> spinham = magnopy.examples.cubic_ferro_nn()
            >>> lswt = magnopy.LSWT(spinham=spinham, spin_directions=[[0, 0, 1]])
            >>> lswt.O()
            array([0.+0.j])
        """

        result = np.zeros((self.M), dtype=complex)
        for _, alphas, parameter in _InteractionParametersIterator(
            self._parameters, n=1, p_n=1
        ):
            result[alphas[0]] = (
                np.sqrt(self.spins[alphas[0]] / 2)
                * np.conjugate(self.p[alphas[0]])
                @ parameter
            )

        # Convert units if necessary
        if units != "meV":
            units = _validated_units(units=units, supported_units=_ENERGY_UNITS)
            result = result * _ENERGY_UNITS["mev"] / _ENERGY_UNITS[units]

        return result

    def A(self, k, relative=False, units="meV"):
        r"""
        Computes part of the grand dynamical matrix.

        Parameters
        ----------

        k : (3,) |array-like|_
            Reciprocal vector

        relative : bool, default False
            If ``relative=True``, then ``k`` is interpreted as given relative to the
            reciprocal unit cell. Otherwise it is interpreted as given in absolute
            coordinates.

        units : str, default "meV"
            .. versionadded:: 0.3.0

            Units of energy. See :ref:`user-guide_usage_units_energy` for the full
            list of supported units.


        Returns
        -------

        A : (M, M) :numpy:`ndarray`
            :math:`A_{\alpha\beta}(\boldsymbol{k})`.

        Notes
        -----

        Before the diagonalization, the magnon Hamiltonian has the form

        .. math::

            \mathcal{H}
            =
            \dots
            +
            \sum_{\boldsymbol{k}, \alpha}
            \boldsymbol{\mathcal{A}}(\boldsymbol{k})^{\dagger}
            \begin{pmatrix}
            \boldsymbol{A}(\boldsymbol{k}) & \boldsymbol{B}^{\dagger}(\boldsymbol{k}) \\
            \boldsymbol{B}(\boldsymbol{k}) & \overline{\boldsymbol{A}(-\boldsymbol{k})}
            \end{pmatrix}
            \boldsymbol{\mathcal{A}}(\boldsymbol{k})

        where

        .. math::
            \boldsymbol{\mathcal{A}}(\boldsymbol{k})
            =
            \begin{pmatrix}
            a_1(\boldsymbol{k}),
            \dots,
            a_M(\boldsymbol{k}),
            a^{\dagger}_1(-\boldsymbol{k}),
            \dots,
            a^{\dagger}_M(-\boldsymbol{k}),
            \end{pmatrix}

        This function computes the matrix :math:`\boldsymbol{A}(\boldsymbol{k})`.

        See Also
        --------

        LSWT.B
        LSWT.GDM

        Examples
        --------

        .. doctest::

            >>> import magnopy
            >>> spinham = magnopy.examples.cubic_ferro_nn()
            >>> lswt = magnopy.LSWT(spinham=spinham, spin_directions=[[0, 0, 1]])
            >>> lswt.A(k=[0, 0, 0.5], relative=True)
            array([[1.+0.j]])
        """

        k = np.array(k)

        result = np.zeros((self.M, self.M), dtype=complex)

        for nu in self.A2:
            if relative:
                phase = 2 * np.pi * (k @ nu)
            else:
                phase = k @ (nu @ self.cell)
            result = result + self.A2[nu] * np.exp(1j * phase)

        result = result - np.diag(self.A1)

        # Convert units if necessary
        if units != "meV":
            units = _validated_units(units=units, supported_units=_ENERGY_UNITS)
            result = result * _ENERGY_UNITS["mev"] / _ENERGY_UNITS[units]

        return result

    def B(self, k, relative=False, units="meV"):
        r"""
        Computes part of the grand dynamical matrix.

        Parameters
        ----------

        k : (3,) |array-like|_
            Reciprocal vector.

        relative : bool, default False
            If ``relative=True``, then ``k`` is interpreted as given relative to the
            reciprocal unit cell. Otherwise it is interpreted as given in absolute
            coordinates.

        units : str, default "meV"
            .. versionadded:: 0.3.0

            Units of energy. See :ref:`user-guide_usage_units_energy` for the full
            list of supported units.


        Returns
        -------

        B : (M, M) :numpy:`ndarray`
            :math:`B_{\alpha\beta}(\boldsymbol{k})`.

        Notes
        -----

        Before the diagonalization, the magnon Hamiltonian has the form

        .. math::

            \mathcal{H}
            =
            \dots
            +
            \sum_{\boldsymbol{k}, \alpha}
            \boldsymbol{\mathcal{A}}(\boldsymbol{k})^{\dagger}
            \begin{pmatrix}
            \boldsymbol{A}(\boldsymbol{k}) & \boldsymbol{B}^{\dagger}(\boldsymbol{k}) \\
            \boldsymbol{B}(\boldsymbol{k}) & \overline{\boldsymbol{A}(-\boldsymbol{k})}
            \end{pmatrix}
            \boldsymbol{\mathcal{A}}(\boldsymbol{k})

        where

        .. math::
            \boldsymbol{\mathcal{A}}(\boldsymbol{k})
            =
            \begin{pmatrix}
            a_1(\boldsymbol{k}),
            \dots,
            a_M(\boldsymbol{k}),
            a^{\dagger}_1(-\boldsymbol{k}),
            \dots,
            a^{\dagger}_M(-\boldsymbol{k}),
            \end{pmatrix}

        This function computes the matrix :math:`\boldsymbol{B}(\boldsymbol{k})`.

        See Also
        --------

        LSWT.A
        LSWT.GDM

        Examples
        --------

        .. doctest::

            >>> import magnopy
            >>> spinham = magnopy.examples.cubic_ferro_nn()
            >>> lswt = magnopy.LSWT(spinham=spinham, spin_directions=[[0, 0, 1]])
            >>> lswt.B(k=[0, 0, 0.5], relative=True)
            array([[0.+0.j]])
        """

        k = np.array(k)

        result = np.zeros((self.M, self.M), dtype=complex)

        for nu in self.B2:
            if relative:
                phase = 2 * np.pi * (k @ nu)
            else:
                phase = k @ (nu @ self.cell)
            result = result + self.B2[nu] * np.exp(1j * phase)

        # Convert units if necessary
        if units != "meV":
            units = _validated_units(units=units, supported_units=_ENERGY_UNITS)
            result = result * _ENERGY_UNITS["mev"] / _ENERGY_UNITS[units]

        return result

    def GDM(self, k, relative=False, units="meV"):
        r"""
        Computes grand dynamical matrix.

        Parameters
        ----------

        k : (3,) |array-like|_
            Reciprocal vector.

        relative : bool, default False
            If ``relative=True``, then ``k`` is interpreted as given relative to the
            reciprocal unit cell. Otherwise it is interpreted as given in absolute
            coordinates.

        units : str, default "meV"
            .. versionadded:: 0.3.0

            Units of energy. See :ref:`user-guide_usage_units_energy` for the full
            list of supported units.


        Returns
        -------

        gdm : (2M, 2M) :numpy:`ndarray`
            Gran dynamical matrix.

        Notes
        -----

        Before the diagonalization, the magnon Hamiltonian has the form

        .. math::

            \mathcal{H}
            =
            \dots
            +
            \sum_{\boldsymbol{k}, \alpha}
            \boldsymbol{\mathcal{A}}(\boldsymbol{k})^{\dagger}
            \begin{pmatrix}
            \boldsymbol{A}(\boldsymbol{k}) & \boldsymbol{B}^{\dagger}(\boldsymbol{k}) \\
            \boldsymbol{B}(\boldsymbol{k}) & \overline{\boldsymbol{A}(-\boldsymbol{k})}
            \end{pmatrix}
            \boldsymbol{\mathcal{A}}(\boldsymbol{k})

        where

        .. math::
            \boldsymbol{\mathcal{A}}(\boldsymbol{k})
            =
            \begin{pmatrix}
            a_1(\boldsymbol{k}),
            \dots,
            a_M(\boldsymbol{k}),
            a^{\dagger}_1(-\boldsymbol{k}),
            \dots,
            a^{\dagger}_M(-\boldsymbol{k}),
            \end{pmatrix}

        This function computes the grand dynamical matrix
        :math:`\boldsymbol{D}(\boldsymbol{k})`

        .. math::

            \boldsymbol{D}(\boldsymbol{k})
            =
            \begin{pmatrix}
            \boldsymbol{A}(\boldsymbol{k}) & \boldsymbol{B}^{\dagger}(\boldsymbol{k}) \\
            \boldsymbol{B}(\boldsymbol{k}) & \overline{\boldsymbol{A}(-\boldsymbol{k})}
            \end{pmatrix}

        See Also
        --------

        LSWT.A
        LSWT.B

        Examples
        --------

        .. doctest::

            >>> import magnopy
            >>> spinham = magnopy.examples.cubic_ferro_nn()
            >>> lswt = magnopy.LSWT(spinham=spinham, spin_directions=[[0, 0, 1]])
            >>> lswt.GDM(k=[0,0,0.5], relative=True) # doctest: +SKIP
            array([[1.+0.j, 0.-0.j],
                   [0.+0.j, 1.-0.j]])
        """

        k = np.array(k, dtype=float)

        A = self.A(k=k, relative=relative, units=units)
        A_m = self.A(k=-k, relative=relative, units=units)

        B = self.B(k=k, relative=relative, units=units)

        left = np.concatenate((A, np.conjugate(B).T), axis=0)
        right = np.concatenate((B, np.conjugate(A_m)), axis=0)
        gdm = np.concatenate((left, right), axis=1)

        return gdm

    def diagonalize(self, k, relative=False, units="meV"):
        r"""
        Diagonalizes the Hamiltonian for the given ``k`` point.

        Parameters
        ----------

        k : (3,) |array-like|_
            Reciprocal vector.

        relative : bool, default False
            If ``relative=True``, then ``k`` is interpreted as given relative to the
            reciprocal unit cell. Otherwise it is interpreted as given in absolute
            coordinates.

        units : str, default "meV"
            .. versionadded:: 0.3.0

            Units of energy. See :ref:`user-guide_usage_units_magnon-energy` for the
            full list of supported units.

        Returns
        -------

        omegas : (M, ) :numpy:`ndarray`
            Array of omegas. Note, that data type is complex. If the ground state is
            correct, then the complex part should be zero.

        delta : float
            Constant energy term that results from diagonalization. Note, that data type
            is complex. If the ground state is correct, then the complex part should be
            zero.

        G : (M, 2M) :numpy:`ndarray`
            Transformation matrix from the original boson operators.

            .. math::

                \begin{pmatrix}
                    b_1(\boldsymbol{k}) \\
                    \dots \\
                    b_M(\boldsymbol{k}) \\
                \end{pmatrix}
                =
                \mathcal{G}
                \begin{pmatrix}
                    a_1(\boldsymbol{k}) \\
                    \dots \\
                    a_M(\boldsymbol{k}) \\
                    a^{\dagger}_1(-\boldsymbol{k}) \\
                    \dots \\
                    a^{\dagger}_M(-\boldsymbol{k}) \\
                \end{pmatrix}

        See Also
        --------

        LSWT.omega
        LSWT.delta
        LSWT.G

        Examples
        --------

        .. doctest::

            >>> import magnopy
            >>> spinham = magnopy.examples.cubic_ferro_nn()
            >>> lswt = magnopy.LSWT(spinham=spinham, spin_directions=[[0, 0, 1]])
            >>> lswt.diagonalize(k=[0, 0, 0.5], relative=True) # doctest: +SKIP
            (array([2.+0.j]), 0j, array([[1.+0.j, 0.+0.j]]))
        """

        GDM = self.GDM(np.array(k, dtype=float), relative=relative)

        # Diagonalize via Colpa's method
        try:
            E, G = solve_via_colpa(GDM)
        except ColpaFailed:
            # Try to diagonalize with suspected Goldstone mode
            try:
                E, G = solve_via_colpa(
                    GDM + (1e-10) * np.eye(GDM.shape[0], dtype=float),
                )
            # Return NaNs if it still fails
            except ColpaFailed:
                # Try to diagonalize for the negative GDMs
                # Note: solve_via_colpa will return positive eigenvalues,
                # so we need to negate them back
                try:
                    E, G = solve_via_colpa(-GDM)
                    E = -E
                except ColpaFailed:
                    return (
                        [np.nan for _ in range(self.M)],
                        np.nan,
                        [[np.nan for _ in range(2 * self.M)] for _ in range(self.M)],
                    )

        # Convert units if necessary
        if units != "meV":
            units = _validated_units(units=units, supported_units=_MAGNON_ENERGY_UNITS)
            tmp_factor = _MAGNON_ENERGY_UNITS["mev"] / _MAGNON_ENERGY_UNITS[units]
            E = E * tmp_factor

        # Factor of two explained in the paper (TODO: add doi after publication)
        energies = E[: self.M] * 2
        transformation_matrices = G[: self.M]

        return (
            energies,  # energies (M)
            complex(0.5 * (np.sum(E[self.M :]) - np.sum(E[: self.M]))),  # delta term
            transformation_matrices,  # transformation matrix (M x 2M)
        )

    def omega(self, k, relative=False, units="meV"):
        r"""
        Computes magnon's eigenenergies at the given ``k`` point.

        Parameters
        ----------

        k : (3,) |array-like|_
            Reciprocal vector.

        relative : bool, default False
            If ``relative=True``, then ``k`` is interpreted as given relative to the
            reciprocal unit cell. Otherwise it is interpreted as given in absolute
            coordinates.

        units : str, default "meV"
            .. versionadded:: 0.3.0

            Units of energy. See :ref:`user-guide_usage_units_magnon-energy` for the
            full list of supported units.


        Returns
        -------
        omegas : (M, ) :numpy:`ndarray`
            Array of omegas. Note, that data type is complex. If the ground state is correct,
            then the complex part should be zero.

        See Also
        --------

        LSWT.diagonalize
        LSWT.delta

        Examples
        --------

        .. doctest::

            >>> import magnopy
            >>> spinham = magnopy.examples.cubic_ferro_nn()
            >>> lswt = magnopy.LSWT(spinham=spinham, spin_directions=[[0, 0, 1]])
            >>> lswt.omega(k=[0, 0, 0.5], relative=True)
            array([2.])
        """

        return self.diagonalize(k=k, relative=relative, units=units)[0]

    def delta(self, k, relative=False, units="meV"):
        r"""
        Computes constant delta term of the diagonalized Hamiltonian.

        .. math::

            \sum_{\boldsymbol{k}}\Delta(\boldsymbol{k})

        Parameters
        ----------

        k : (3,) |array-like|_
            Reciprocal vector.

        relative : bool, default False
            If ``relative=True``, then ``k`` is interpreted as given relative to the
            reciprocal unit cell. Otherwise it is interpreted as given in absolute
            coordinates.

        units : str, default "meV"
            .. versionadded:: 0.3.0

            Units of energy. See :ref:`user-guide_usage_units_magnon-energy` for the
            full list of supported units.


        Returns
        -------

        delta : float
            Constant energy term that results from diagonalization. Note, that data type is complex. If the ground state is correct,
            then the complex part should be zero.

        See Also
        --------

        LSWT.diagonalize
        LSWT.omega

        Examples
        --------

        .. doctest::

            >>> import magnopy
            >>> spinham = magnopy.examples.cubic_ferro_nn()
            >>> lswt = magnopy.LSWT(spinham=spinham, spin_directions=[[0, 0, 1]])
            >>> lswt.delta(k=[0, 0, 0.5], relative=True)
            0j
        """
        return self.diagonalize(k=k, relative=relative, units=units)[1]

    def G(self, k, relative=False):
        r"""
        Computes transformation matrix to the new bosonic operators.

        .. math::

            b_{\alpha}(\boldsymbol{k})
            =
            \sum_{\beta}
            (\mathcal{G})_{\alpha, \beta}(\boldsymbol{k})
            \mathcal{A}_{\beta}(\boldsymbol{k})

        Parameters
        ----------

        k : (3,) |array-like|_
            Reciprocal vector

        relative : bool, default False
            If ``relative=True``, then ``k`` is interpreted as given relative to the
            reciprocal unit cell. Otherwise it is interpreted as given in absolute
            coordinates.

        Returns
        -------

        G : (M, 2M) :numpy:`ndarray`
            Transformation matrix from the original boson operators.

        See Also
        --------

        LSWT.diagonalize
        LSWT.omega
        LSWT.delta

        Examples
        --------

        .. doctest::

            >>> import magnopy
            >>> spinham = magnopy.examples.cubic_ferro_nn()
            >>> lswt = magnopy.LSWT(spinham=spinham, spin_directions=[[0, 0, 1]])
            >>> lswt.G(k=[0, 0, 0.5], relative=True)  # doctest: +SKIP
            array([[1.+0.j, 0.+0.j]])
        """
        return self.diagonalize(k=k, relative=relative)[2]


# Populate __all__ with objects defined in this file
__all__ = list(set(dir()) - old_dir)
# Remove all semi-private objects
__all__ = [i for i in __all__ if not i.startswith("_")]
del old_dir
