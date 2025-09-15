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

from magnopy._diagonalization import solve_via_colpa
from magnopy._exceptions import ColpaFailed
from magnopy._local_rf import span_local_rfs

# Save local scope at this moment
old_dir = set(dir())
old_dir.add("old_dir")


class LSWT:
    r"""
    Linear Spin Wave theory.

    It is created from some spin Hamitlonian and set of direction vectors,
    that define the orientation of spins.

    Parameters
    ----------
    spinham : :py:class:`.SpinHamiltonian`
        Spin Hamiltonian.
    spin_directions : (M, 3) |array-like|_
        Directions of spin vectors. Only directions of vectors are used, modulus is
        ignored.
        If spin Hamiltonian contains non-magnetic atom, then only the spin directions
        for the magnetic atoms are expected. The order of spin directions is the same as
        the order of magnetic atoms in ``spinham.atoms.spins``.

    Attributes
    ----------
    z : (M, 3) :numpy:`ndarray`
        Spin directions.
    p : (M, 3) :numpy:`ndarray`
        Hybridized x and y components of the local coordinate system.
    M : int
        Number of spins in the unit cell
    cell : (3, 3) :numpy:`ndarray`
        Unit cell. Rows are vectors, columns are cartesian components.
    spins : (M, ) :numpy:`ndarray`
        Spin values.

    Notes
    -----
    If spin Hamiltonian contains three atoms Cr1, Br, Cr3 in that order. Assume that two
    atoms are magnetic (Cr1 and Cr2), one atom is not (Br). Then ``spin_directions`` is
    a (2, 3) array with ``spin_directions[0]`` being the direction for spin of Cr1 and
    ``spin_directions[1]`` being the direction of spin for Cr2.

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

        initial_convention = spinham.convention

        spinham.convention = initial_convention.get_modified(
            spin_normalized=False, multiple_counting=True
        )

        self.M = spinham.M
        self.cell = spinham.cell

        ########################################################################
        #                    Renormalized one-spin parameter                   #
        ########################################################################
        self._J1 = np.zeros((self.M, 3), dtype=float)

        # One spin
        for alpha, parameter in spinham.p1:
            alpha = spinham.map_to_magnetic[alpha]
            self._J1[alpha] = self._J1[alpha] + spinham.convention.c1 * parameter

        # Two spins & one site
        for alpha, parameter in spinham.p21:
            alpha = spinham.map_to_magnetic[alpha]
            self._J1[alpha] = self._J1[alpha] + (
                2
                * spinham.convention.c21
                * np.einsum("ij,j->i", parameter, self.z[alpha])
                * self.spins[alpha]
            )

        # Two spins & two sites
        for alpha, beta, _, parameter in spinham.p22:
            alpha = spinham.map_to_magnetic[alpha]
            beta = spinham.map_to_magnetic[beta]
            self._J1[alpha] = self._J1[alpha] + (
                2
                * spinham.convention.c22
                * (parameter @ self.z[beta])
                * self.spins[beta]
            )

        # Three spins & one site
        for alpha, parameter in spinham.p31:
            alpha = spinham.map_to_magnetic[alpha]
            self._J1[alpha] = self._J1[alpha] + (
                3
                * spinham.convention.c31
                * np.einsum("iju,j,u->i", parameter, self.z[alpha], self.z[alpha])
                * self.spins[alpha] ** 2
            )

        # Three spins & two sites
        for alpha, beta, _, parameter in spinham.p32:
            alpha = spinham.map_to_magnetic[alpha]
            beta = spinham.map_to_magnetic[beta]
            self._J1[alpha] = self._J1[alpha] + (
                3
                * spinham.convention.c32
                * np.einsum("iju,j,u->i", parameter, self.z[alpha], self.z[beta])
                * self.spins[alpha]
                * self.spins[beta]
            )

        # Three spins & three sites
        for alpha, beta, gamma, _, _, parameter in spinham.p33:
            alpha = spinham.map_to_magnetic[alpha]
            beta = spinham.map_to_magnetic[beta]
            gamma = spinham.map_to_magnetic[gamma]
            self._J1[alpha] = self._J1[alpha] + (
                3
                * spinham.convention.c33
                * np.einsum("iju,j,u->i", parameter, self.z[beta], self.z[gamma])
                * self.spins[beta]
                * self.spins[gamma]
            )

        # Four spins & one site
        for alpha, parameter in spinham.p41:
            alpha = spinham.map_to_magnetic[alpha]
            self._J1[alpha] = self._J1[alpha] + (
                4
                * spinham.convention.c41
                * np.einsum(
                    "ijuv,j,u,v->i",
                    parameter,
                    self.z[alpha],
                    self.z[alpha],
                    self.z[alpha],
                )
                * self.spins[alpha] ** 3
            )

        # Four spins & two sites (1+3)
        for alpha, beta, _, parameter in spinham.p421:
            alpha = spinham.map_to_magnetic[alpha]
            beta = spinham.map_to_magnetic[beta]
            self._J1[alpha] = self._J1[alpha] + (
                4
                * spinham.convention.c421
                * np.einsum(
                    "ijuv,j,u,v->i",
                    parameter,
                    self.z[alpha],
                    self.z[alpha],
                    self.z[beta],
                )
                * self.spins[alpha] ** 2
                * self.spins[beta]
            )

        # Four spins & two sites (2+2)
        for alpha, beta, _, parameter in spinham.p422:
            alpha = spinham.map_to_magnetic[alpha]
            beta = spinham.map_to_magnetic[beta]
            self._J1[alpha] = self._J1[alpha] + (
                4
                * spinham.convention.c422
                * np.einsum(
                    "ijuv,j,u,v->i",
                    parameter,
                    self.z[alpha],
                    self.z[beta],
                    self.z[beta],
                )
                * self.spins[alpha]
                * self.spins[beta] ** 2
            )

        # Four spins & three sites
        for alpha, beta, gamma, _, _, parameter in spinham.p43:
            alpha = spinham.map_to_magnetic[alpha]
            beta = spinham.map_to_magnetic[beta]
            gamma = spinham.map_to_magnetic[gamma]
            self._J1[alpha] = self._J1[alpha] + (
                4
                * spinham.convention.c43
                * np.einsum(
                    "ijuv,j,u,v->i",
                    parameter,
                    self.z[alpha],
                    self.z[beta],
                    self.z[gamma],
                )
                * self.spins[alpha]
                * self.spins[beta]
                * self.spins[gamma]
            )

        # Four spins & four sites
        for alpha, beta, gamma, epsilon, _, _, _, parameter in spinham.p44:
            alpha = spinham.map_to_magnetic[alpha]
            beta = spinham.map_to_magnetic[beta]
            gamma = spinham.map_to_magnetic[gamma]
            epsilon = spinham.map_to_magnetic[epsilon]
            self._J1[alpha] = self._J1[alpha] + (
                4
                * spinham.convention.c44
                * np.einsum(
                    "ijuv,j,u,v->i",
                    parameter,
                    self.z[beta],
                    self.z[gamma],
                    self.z[epsilon],
                )
                * self.spins[beta]
                * self.spins[gamma]
                * self.spins[epsilon]
            )

        ########################################################################
        #                   Renormalized two-spins parameter                   #
        ########################################################################
        self._J2 = {}

        # First - terms with delta in from of them

        # Two spins & one site
        for alpha, parameter in spinham.p21:
            alpha = spinham.map_to_magnetic[alpha]
            if (0, 0, 0) not in self._J2:
                self._J2[(0, 0, 0)] = np.zeros((self.M, self.M, 3, 3), dtype=float)

            self._J2[(0, 0, 0)][alpha, alpha] += 2 * spinham.convention.c21 * parameter

        # Three spins & one site
        for alpha, parameter in spinham.p31:
            alpha = spinham.map_to_magnetic[alpha]
            if (0, 0, 0) not in self._J2:
                self._J2[(0, 0, 0)] = np.zeros((self.M, self.M, 3, 3), dtype=float)

            self._J2[(0, 0, 0)][alpha, alpha] += (
                3
                * spinham.convention.c31
                * np.einsum("iju,u->ij", parameter, self.z[alpha])
                * self.spins[alpha]
            )

        # Four spins & one site
        for alpha, parameter in spinham.p41:
            alpha = spinham.map_to_magnetic[alpha]
            if (0, 0, 0) not in self._J2:
                self._J2[(0, 0, 0)] = np.zeros((self.M, self.M, 3, 3), dtype=float)

            self._J2[(0, 0, 0)][alpha, alpha] += (
                6
                * spinham.convention.c41
                * np.einsum("ijuv,u,v->ij", parameter, self.z[alpha], self.z[alpha])
                * self.spins[alpha] ** 2
            )

        # Then all other parameters

        # Two spins & two sites
        for alpha, beta, nu, parameter in spinham.p22:
            alpha = spinham.map_to_magnetic[alpha]
            beta = spinham.map_to_magnetic[beta]
            if nu not in self._J2:
                self._J2[nu] = np.zeros((self.M, self.M, 3, 3), dtype=float)

            self._J2[nu][alpha, beta] += spinham.convention.c22 * parameter

        # Three spins & two sites
        for alpha, beta, nu, parameter in spinham.p32:
            alpha = spinham.map_to_magnetic[alpha]
            beta = spinham.map_to_magnetic[beta]
            if nu not in self._J2:
                self._J2[nu] = np.zeros((self.M, self.M, 3, 3), dtype=float)

            self._J2[nu][alpha, beta] += (
                3
                * spinham.convention.c32
                * np.einsum("iuj,u->ij", parameter, self.z[alpha])
                * self.spins[alpha]
            )

        # Three spins & three sites
        for alpha, beta, gamma, nu, _, parameter in spinham.p33:
            alpha = spinham.map_to_magnetic[alpha]
            beta = spinham.map_to_magnetic[beta]
            gamma = spinham.map_to_magnetic[gamma]
            if nu not in self._J2:
                self._J2[nu] = np.zeros((self.M, self.M, 3, 3), dtype=float)

            self._J2[nu][alpha, beta] += (
                3
                * spinham.convention.c33
                * np.einsum("iju,u->ij", parameter, self.z[gamma])
                * self.spins[gamma]
            )

        # Four spins & two sites (1+3)
        for alpha, beta, nu, parameter in spinham.p421:
            alpha = spinham.map_to_magnetic[alpha]
            beta = spinham.map_to_magnetic[beta]
            if nu not in self._J2:
                self._J2[nu] = np.zeros((self.M, self.M, 3, 3), dtype=float)

            self._J2[nu][alpha, beta] += (
                6
                * spinham.convention.c421
                * np.einsum("iuvj,u,v->ij", parameter, self.z[alpha], self.z[alpha])
                * self.spins[alpha] ** 2
            )

        # Four spins & two sites (2+2)
        for alpha, beta, nu, parameter in spinham.p422:
            alpha = spinham.map_to_magnetic[alpha]
            beta = spinham.map_to_magnetic[beta]
            if nu not in self._J2:
                self._J2[nu] = np.zeros((self.M, self.M, 3, 3), dtype=float)

            self._J2[nu][alpha, beta] += (
                6
                * spinham.convention.c422
                * np.einsum("iujv,u,v->ij", parameter, self.z[alpha], self.z[beta])
                * self.spins[alpha]
                * self.spins[beta]
            )

        # Four spins & three sites
        for alpha, beta, gamma, nu, _, parameter in spinham.p43:
            alpha = spinham.map_to_magnetic[alpha]
            beta = spinham.map_to_magnetic[beta]
            gamma = spinham.map_to_magnetic[gamma]
            if nu not in self._J2:
                self._J2[nu] = np.zeros((self.M, self.M, 3, 3), dtype=float)

            self._J2[nu][alpha, beta] += (
                6
                * spinham.convention.c43
                * np.einsum("iujv,u->ij", parameter, self.z[alpha], self.z[gamma])
                * self.spins[alpha]
                * self.spins[gamma]
            )

        # Four spins & four sites
        for alpha, beta, gamma, epsilon, nu, _, _, parameter in spinham.p44:
            alpha = spinham.map_to_magnetic[alpha]
            beta = spinham.map_to_magnetic[beta]
            gamma = spinham.map_to_magnetic[gamma]
            epsilon = spinham.map_to_magnetic[epsilon]
            if nu not in self._J2:
                self._J2[nu] = np.zeros((self.M, self.M, 3, 3), dtype=float)

            self._J2[nu][alpha, beta] += (
                6
                * spinham.convention.c44
                * np.einsum("ijuv,u->ij", parameter, self.z[gamma], self.z[epsilon])
                * self.spins[gamma]
                * self.spins[epsilon]
            )

        spinham.convention = initial_convention

        self.A1 = 0.5 * np.sum(self._J1 * self.z, axis=1)

        self.A2 = {}
        self.B2 = {}

        for nu in self._J2:
            self.A2[nu] = 0.5 * np.einsum(
                "abij,a,b,ai,bj->ab",
                self._J2[nu],
                np.sqrt(self.spins),
                np.sqrt(self.spins),
                self.p,
                np.conjugate(self.p),
            )
            self.B2[nu] = 0.5 * np.einsum(
                "abij,a,b,ai,bj->ab",
                self._J2[nu],
                np.sqrt(self.spins),
                np.sqrt(self.spins),
                np.conjugate(self.p),
                np.conjugate(self.p),
            )

    @property
    def E_2(self) -> float:
        r"""
        Correction to the ground state energy that arises from the LSWT.

        Returns
        -------
        E_2 : float

        Examples
        --------

        .. doctest::

            >>> import magnopy
            >>> spinham = magnopy.examples.cubic_ferro_nn()
            >>> lswt = magnopy.LSWT(spinham=spinham, spin_directions=[[0, 0, 1]])
            >>> lswt.E_2
            -1.5
        """

        return float(0.5 * np.sum(self._J1 * self.z))

    @property
    def O(self):  # noqa E743
        r"""
        Coefficient before the one-operator terms.

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
            >>> lswt.O
            array([0.+0.j])
        """

        return np.einsum(
            "a,ai,ai->a",
            np.sqrt(self.spins) / np.sqrt(2),
            np.conjugate(self.p),
            self._J1,
        )

    def A(self, k, relative=False):
        r"""
        Part of the Grand dynamical matrix.

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

        return result

    def B(self, k, relative=False):
        r"""
        Part of the Grand dynamical matrix.

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

        result = result

        return result

    def GDM(self, k, relative=False):
        r"""
        Grand dynamical matrix.

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

        A = self.A(k=k, relative=relative)
        A_m = self.A(k=-k, relative=relative)

        B = self.B(k=k, relative=relative)

        left = np.concatenate((A, np.conjugate(B).T), axis=0)
        right = np.concatenate((B, np.conjugate(A_m)), axis=0)
        gdm = np.concatenate((left, right), axis=1)

        return gdm

    def diagonalize(self, k, relative=False):
        r"""
        Diagonalize the Hamiltonian for the given ``k`` point and return all possible
        quantities at once.

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
        omegas : (M, ) :numpy:`ndarray`
            Array of omegas. Note, that data type is complex. If the ground state is correct,
            then the complex part should be zero.
        delta : float
            Constant energy term that results from diagonalization. Note, that data type is complex. If the ground state is correct,
            then the complex part should be zero.
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

        Examples
        --------

        .. doctest::

            >>> import magnopy
            >>> spinham = magnopy.examples.cubic_ferro_nn()
            >>> lswt = magnopy.LSWT(spinham=spinham, spin_directions=[[0, 0, 1]])
            >>> lswt.diagonalize(k=[0, 0, 0.5], relative=True) # doctest: +SKIP
            (array([2.+0.j]), 0j, array([[1.+0.j, 0.+0.j]]))
        """

        k_plus = np.array(k)
        k_minus = -k_plus

        GDM_plus = self.GDM(k_plus, relative=relative)
        GDM_minus = self.GDM(k_minus, relative=relative)

        try:
            E_plus, G_plus = solve_via_colpa(GDM_plus, sort_by_first_N=True)
            E_minus, G_minus = solve_via_colpa(GDM_minus, sort_by_first_N=False)
        except ColpaFailed:
            try:
                E_plus, G_plus = solve_via_colpa(-GDM_plus, sort_by_first_N=True)
                E_minus, G_minus = solve_via_colpa(-GDM_minus, sort_by_first_N=False)
            except ColpaFailed:
                try:
                    E_plus, G_plus = solve_via_colpa(
                        GDM_plus + (1e-8) * np.ones(GDM_plus.shape, dtype=float),
                        sort_by_first_N=True,
                    )
                    E_minus, G_minus = solve_via_colpa(
                        GDM_minus + (1e-8) * np.ones(GDM_minus.shape, dtype=float),
                        sort_by_first_N=False,
                    )
                except ColpaFailed:
                    return (
                        [np.nan for _ in range(self.M)],
                        np.nan,
                        [[np.nan for _ in range(2 * self.M)] for _ in range(self.M)],
                    )

        # Sort by energy values

        energies = E_plus[: self.M] + E_minus[self.M :]
        transformation_matrices = G_plus[: self.M]

        sorting_indices = np.argsort(energies)

        return (
            energies[sorting_indices],
            complex(0.5 * (np.sum(E_plus[self.M :]) - np.sum(E_plus[: self.M]))),
            transformation_matrices[sorting_indices],
        )

    def omega(self, k, relative=False):
        r"""
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
            array([2.+0.j])
        """

        return self.diagonalize(k=k, relative=relative)[0]

    def delta(self, k, relative=False):
        r"""
        Constant energy term of the diagonalized Hamiltonian.

        .. math::

            \sum_{\boldsymbol{k}}\Delta(\boldsymbol{k})

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
        return self.diagonalize(k=k, relative=relative)[1]

    def G(self, k, relative=False):
        r"""
        Transformation matrix to the new bosonic operators.

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

    # Removed in v0.2.0. Warning will be removed in March of 2026
    def G_inv(self, *args, **kwargs):
        raise DeprecationWarning("This method was removed in v0.2.0 in favor of LSWT.G")


# Populate __all__ with objects defined in this file
__all__ = list(set(dir()) - old_dir)
# Remove all semi-private objects
__all__ = [i for i in __all__ if not i.startswith("_")]
del old_dir


if __name__ == "__main__":
    import magnopy

    spinham = magnopy.io.load_grogu(
        "/Users/rybakov.ad/Codes/magnopy/tmp-09.06.25/Ni-U4/converted_NiPS3.txt"
    )
    sd = np.loadtxt(
        "/Users/rybakov.ad/Codes/magnopy/tmp-09.06.25/Ni-U4/optimize/SPIN_DIRECTIONS.txt"
    )
    lswt = LSWT(spinham, sd)
    k = np.array([0, 0, 0], dtype=float)

    D = lswt.GDM(k)

    E, G = magnopy.solve_via_colpa(D)

    print(np.linalg.inv(np.conjugate(G).T) @ D @ np.linalg.inv(G) == E)
