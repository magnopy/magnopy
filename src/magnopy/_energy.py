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


from math import log10

import numpy as np

# Save local scope at this moment
old_dir = set(dir())
old_dir.add("old_dir")

_C1 = 1e-4
_C2 = 0.9


def _cubic_interpolation(alpha_l, alpha_h, phi_l, phi_h, der_l, der_h):
    r"""
    Computes the minimum of a cubic interpolation for the function f(alpha) with
    values f_l, f_h and derivatives fp_l, fp_h at the points alpha_l, alpha_h.

    Parameters
    ----------
    alpha_l : float
        Lower bound of the interval.
    alpha_h : float
        Upper bound of the interval.
    phi_l : float
        Value of the function at alpha_l.
    phi_h : float
        Value of the function at alpha_h.
    der_l : float
        Derivative of the function at alpha_l.
    der_h : float
        Derivative of the function at alpha_h.

    Returns
    -------
    alpha_min : float
        Position of the minimum of the cubic interpolation for the function f(alpha).
    """

    if abs(alpha_l - alpha_h) < np.finfo(float).eps:
        return alpha_l

    d_1 = der_l + der_h - 3 * (phi_l - phi_h) / (alpha_l - alpha_h)

    if d_1**2 - der_l * der_h < 0:
        if phi_l <= phi_h:
            return alpha_l
        else:
            return alpha_h

    d_2 = np.sign(alpha_h - alpha_l) * np.sqrt(d_1**2 - der_l * der_h)

    return alpha_h - (alpha_h - alpha_l) * (der_h + d_2 - d_1) / (
        der_h - der_l + 2 * d_2
    )


def _rotate_sd(reference_sd, rotation):
    r"""

    Parameters
    ----------
    reference_sd : (M, 3) :numpy:`ndarray`
        Reference direction of the spin vectors.
    rotation : (M*3,) :numpy:`ndarray`
        Rotation of the spin vectors parameterized with the skew-symmetric matrix.

    Returns
    -------
    directions : (I, 3) :numpy:`ndarray`
        Rotated set of direction vectors.
    """

    directions = reference_sd.copy()

    ax = rotation[::3]
    ay = rotation[1::3]
    az = rotation[2::3]

    thetas = np.sqrt(ax**2 + ay**2 + az**2)

    r = []
    for alpha in range(len(thetas)):
        theta = thetas[alpha]

        if theta < np.finfo(float).eps:
            continue

        r = np.array([ax[alpha], ay[alpha], az[alpha]]) / theta

        directions[alpha] = (
            np.cos(theta) * directions[alpha]
            + np.sin(theta) * np.cross(r, directions[alpha])
            + (1 - np.cos(theta)) * r * (r @ directions[alpha])
        )

    return directions


class Energy:
    r"""
    Ground state energy of the given spin Hamiltonian.

    This class is optimized for the computation of the energy for any spin
    directions for the given Hamiltonian.

    If the spin Hamiltonian is modified, then a new instance of the energy class
    should be created from it.

    Parameters
    ----------
    spinham : :py:class:`.SpinHamiltonian`
        Spin Hamiltonian for the calculation of energy.

    Examples
    --------

    First, one need to create some spin Hamiltonian

    .. doctest::

        >>> import numpy as np
        >>> import magnopy
        >>> cell = np.eye(3)
        >>> atoms = dict(
        ...     names=["Fe"], spins=[1.5], g_factors=[2], positions=[[0, 0, 0]]
        ... )
        >>> convention = magnopy.Convention(
        ...     multiple_counting=True, spin_normalized=False, c21=1, c22=-1
        ... )
        >>> spinham = magnopy.SpinHamiltonian(
        ...     cell=cell, atoms=atoms, convention=convention
        ... )

    Then, add some parameters to the Hamiltonian

    .. doctest::

        >>> spinham.add_21(alpha=0, parameter=np.diag([0, 0, -1]))
        >>> spinham.add_22(
        ...     alpha=0,
        ...     beta=0,
        ...     nu=(1, 0, 0),
        ...     parameter=magnopy.converter22.from_iso(iso=1),
        ... )

    Now everything is ready to create an instance of the Energy class

    .. doctest::

        >>> energy = magnopy.Energy(spinham)

    Finally, ``energy`` an be used to compute classical energy of the Hamiltonian,
    its gradient, torque or search for the local minima.

    .. doctest::

        >>> sd1 = [[1, 0, 0]]
        >>> sd2 = [[0, 1, 0]]
        >>> sd3 = [[0, 0, 1]]
        >>> energy(sd1), energy(sd2), energy(sd3)
        (-4.5, -4.5, -6.75)
    """

    def __init__(self, spinham):
        initial_convention = spinham.convention

        magnopy_convention = initial_convention.get_modified(
            spin_normalized=False, multiple_counting=True
        )

        spinham.convention = magnopy_convention

        self.spins = np.array(spinham.magnetic_atoms.spins, dtype=float)
        self.M = spinham.M

        ########################################################################
        #                               One spin                               #
        ########################################################################

        self.J_1 = np.zeros((spinham.M, 3), dtype=float)

        for atom, parameter in spinham.p1:
            alpha = spinham.map_to_magnetic[atom]

            self.J_1[alpha] += spinham.convention.c1 * parameter

        ########################################################################
        #                               Two spins                              #
        ########################################################################

        self.J_21 = np.zeros((spinham.M, 3, 3), dtype=float)

        for atom, parameter in spinham.p21:
            alpha = spinham.map_to_magnetic[atom]

            self.J_21[alpha] += spinham.convention.c21 * parameter

        self.J_22 = {}

        for atom1, atom2, _, parameter in spinham.p22:
            alpha = spinham.map_to_magnetic[atom1]
            beta = spinham.map_to_magnetic[atom2]

            if (alpha, beta) not in self.J_22:
                self.J_22[(alpha, beta)] = np.zeros((3, 3), dtype=float)

            self.J_22[(alpha, beta)] += spinham.convention.c22 * parameter

        ########################################################################
        #                              Three spins                             #
        ########################################################################

        self.J_31 = np.zeros((spinham.M, 3, 3, 3), dtype=float)

        for atom, parameter in spinham.p31:
            alpha = spinham.map_to_magnetic[atom]

            self.J_31[alpha] += spinham.convention.c31 * parameter

        self.J_32 = {}

        for atom1, atom2, _, parameter in spinham.p32:
            alpha = spinham.map_to_magnetic[atom1]
            beta = spinham.map_to_magnetic[atom2]

            if (alpha, beta) not in self.J_32:
                self.J_32[(alpha, beta)] = np.zeros((3, 3, 3), dtype=float)

            self.J_32[(alpha, beta)] += spinham.convention.c32 * parameter

        self.J_33 = {}

        for atom1, atom2, atom3, _, _, parameter in spinham.p33:
            alpha = spinham.map_to_magnetic[atom1]
            beta = spinham.map_to_magnetic[atom2]
            gamma = spinham.map_to_magnetic[atom3]

            if (alpha, beta, gamma) not in self.J_33:
                self.J_33[(alpha, beta, gamma)] = np.zeros((3, 3, 3), dtype=float)

            self.J_33[(alpha, beta, gamma)] += spinham.convention.c33 * parameter

        ########################################################################
        #                              Four spins                              #
        ########################################################################

        self.J_41 = np.zeros((spinham.M, 3, 3, 3, 3), dtype=float)

        for atom, parameter in spinham.p41:
            alpha = spinham.map_to_magnetic[atom]

            self.J_41[alpha] += spinham.convention.c41 * parameter

        self.J_421 = {}

        for atom1, atom2, _, parameter in spinham.p421:
            alpha = spinham.map_to_magnetic[atom1]
            beta = spinham.map_to_magnetic[atom2]

            if (alpha, beta) not in self.J_421:
                self.J_421[(alpha, beta)] = np.zeros((3, 3, 3, 3), dtype=float)

            self.J_421[(alpha, beta)] += spinham.convention.c421 * parameter

        self.J_422 = {}

        for atom1, atom2, _, parameter in spinham.p422:
            alpha = spinham.map_to_magnetic[atom1]
            beta = spinham.map_to_magnetic[atom2]

            if (alpha, beta) not in self.J_422:
                self.J_422[(alpha, beta)] = np.zeros((3, 3, 3, 3), dtype=float)

            self.J_422[(alpha, beta)] += spinham.convention.c422 * parameter

        self.J_43 = {}

        for atom1, atom2, atom3, _, _, parameter in spinham.p43:
            alpha = spinham.map_to_magnetic[atom1]
            beta = spinham.map_to_magnetic[atom2]
            gamma = spinham.map_to_magnetic[atom3]

            if (alpha, beta, gamma) not in self.J_43:
                self.J_43[(alpha, beta, gamma)] = np.zeros((3, 3, 3, 3), dtype=float)

            self.J_43[(alpha, beta, gamma)] += spinham.convention.c43 * parameter

        self.J_44 = {}

        for atom1, atom2, atom3, atom4, _, _, _, parameter in spinham.p44:
            alpha = spinham.map_to_magnetic[atom1]
            beta = spinham.map_to_magnetic[atom2]
            gamma = spinham.map_to_magnetic[atom3]
            epsilon = spinham.map_to_magnetic[atom4]

            if (alpha, beta, gamma, epsilon) not in self.J_44:
                self.J_44[(alpha, beta, gamma, epsilon)] = np.zeros(
                    (3, 3, 3, 3), dtype=float
                )

            self.J_44[(alpha, beta, gamma, epsilon)] += (
                spinham.convention.c44 * parameter
            )

        spinham.convention = initial_convention

    def __call__(self, spin_directions, _normalize=True) -> float:
        return self.E_0(spin_directions=spin_directions, _normalize=_normalize)

    def E_0(self, spin_directions, _normalize=True) -> float:
        r"""
        Computes classical energy of the spin Hamiltonian.

        Parameters
        ----------
        spin_directions : (M, 3) |array-like|_
            Directions of spin vectors. Only directions of vectors are used,
            modulus is ignored. ``M`` is the amount of magnetic atoms in the
            Hamiltonian. The order of spin directions is the same as the order
            of magnetic atoms in ``spinham.magnetic_atoms.spins``.
        _normalize : bool, default True
            Whether to normalize the spin_directions or use the provided vectors as is.
            This parameter is technical and we do not recommend to use it at all.

        Returns
        -------
        E_0 : float
            Classic energy of state with ``spin_directions``.


        Examples
        --------

        First, one need to create some spin Hamiltonian

        .. doctest::

            >>> import numpy as np
            >>> import magnopy
            >>> cell = np.eye(3)
            >>> atoms = dict(
            ...     names=["Fe"], spins=[1.5], g_factors=[2], positions=[[0, 0, 0]]
            ... )
            >>> convention = magnopy.Convention(
            ...     multiple_counting=True, spin_normalized=False, c21=1, c22=-1
            ... )
            >>> spinham = magnopy.SpinHamiltonian(
            ...     cell=cell, atoms=atoms, convention=convention
            ... )

        Then, add some parameters to the Hamiltonian

        .. doctest::

            >>> spinham.add_21(alpha=0, parameter=np.diag([0, 0, -1]))
            >>> spinham.add_22(
            ...     alpha=0,
            ...     beta=0,
            ...     nu=(1, 0, 0),
            ...     parameter=magnopy.converter22.from_iso(iso=1),
            ... )

        Now everything is ready to create an instance of the Energy class

        .. doctest::

            >>> energy = magnopy.Energy(spinham)

        Finally, ``energy`` an be used to compute classical energy of the Hamiltonian
        for arbitrary spin configuration.

        .. doctest::

            >>> sd1 = [[1, 0, 0]]
            >>> sd2 = [[0, 1, 0]]
            >>> sd3 = [[0, 0, 1]]
            >>> energy.E_0(sd1), energy.E_0(sd2), energy.E_0(sd3)
            (-4.5, -4.5, -6.75)
            >>> # The command above is equivalent to
            >>> energy(sd1), energy(sd2), energy(sd3)
            (-4.5, -4.5, -6.75)
        """

        spin_directions = np.array(spin_directions, dtype=float)

        if _normalize:
            spin_directions = (
                spin_directions / np.linalg.norm(spin_directions, axis=1)[:, np.newaxis]
            )
        spins = spin_directions * self.spins[:, np.newaxis]

        energy = 0

        energy += np.diag(self.J_1 @ spins.T).sum()

        energy += np.einsum("mij,mi,mj->m", self.J_21, spins, spins).sum()

        energy += np.einsum("miju,mi,mj,mu->m", self.J_31, spins, spins, spins).sum()

        energy += np.einsum(
            "mijuv,mi,mj,mu,mv->m", self.J_41, spins, spins, spins, spins
        ).sum()

        for alpha, beta in self.J_22:
            energy += spins[alpha] @ self.J_22[(alpha, beta)] @ spins[beta]

        for alpha, beta in self.J_32:
            energy += np.einsum(
                "iju,i,j,u",
                self.J_32[(alpha, beta)],
                spins[alpha],
                spins[alpha],
                spins[beta],
            )

        for alpha, beta in self.J_421:
            energy += np.einsum(
                "ijuv,i,j,u,v",
                self.J_421[(alpha, beta)],
                spins[alpha],
                spins[alpha],
                spins[alpha],
                spins[beta],
            )

        for alpha, beta in self.J_422:
            energy += np.einsum(
                "ijuv,i,j,u,v",
                self.J_422[(alpha, beta)],
                spins[alpha],
                spins[alpha],
                spins[beta],
                spins[beta],
            )

        for alpha, beta, gamma in self.J_33:
            energy += np.einsum(
                "iju,i,j,u",
                self.J_33[(alpha, beta, gamma)],
                spins[alpha],
                spins[beta],
                spins[gamma],
            )

        for alpha, beta, gamma in self.J_43:
            energy += np.einsum(
                "ijuv,i,j,u,v",
                self.J_43[(alpha, beta, gamma)],
                spins[alpha],
                spins[alpha],
                spins[beta],
                spins[gamma],
            )

        for alpha, beta, gamma, epsilon in self.J_44:
            energy += np.einsum(
                "ijuv,i,j,u,v",
                self.J_44[(alpha, beta, gamma, epsilon)],
                spins[alpha],
                spins[beta],
                spins[gamma],
                spins[epsilon],
            )

        return float(energy)

    def gradient(self, spin_directions, _normalize=True):
        r"""
        Computes first derivatives of energy (:math:`E^{(0)}`) with respect to the
        components of the spin directional vectors.

        Parameters
        ----------
        spin_directions : (M, 3) |array-like|_
            Directions of spin vectors. Only directions of vectors are used,
            modulus is ignored. ``M`` is the amount of magnetic atoms in the
            Hamiltonian. The order of spin directions is the same as the order
            of magnetic atoms in ``spinham.magnetic_atoms.spins``.
        _normalize : bool, default True
            Whether to normalize the spin_directions or use the provided vectors as is.
            This parameter is technical and we do not recommend to use it at all.

        Returns
        -------
        gradient : (M, 3) :numpy:`ndarray`
            Gradient of energy.

            .. code-block:: python

                [
                    [dE / dz1x, dE / dz1y, dE / dz1z],
                    [dE / dz2x, dE / dz2y, dE / dz2z],
                    ...[dE / dzMx, dE / dzMy, dE / dzMz],
                ]
        """

        spin_directions = np.array(spin_directions, dtype=float)

        if _normalize:
            spin_directions = (
                spin_directions / np.linalg.norm(spin_directions, axis=1)[:, np.newaxis]
            )

        gradient = np.zeros((self.M, 3), dtype=float)

        gradient += self.J_1 * self.spins[:, np.newaxis]

        gradient += 2 * np.einsum(
            "mtj,mj,m->mt", self.J_21, spin_directions, self.spins**2
        )

        gradient += 3 * np.einsum(
            "mtju,mj,mu,m->mt",
            self.J_31,
            spin_directions,
            spin_directions,
            self.spins**3,
        )

        gradient += 4 * np.einsum(
            "mtjuv,mj,mu,mv,m->mt",
            self.J_41,
            spin_directions,
            spin_directions,
            spin_directions,
            self.spins**4,
        )

        for alpha, beta in self.J_22:
            gradient[alpha] += 2 * (
                self.J_22[(alpha, beta)]
                @ spin_directions[beta]
                * self.spins[alpha]
                * self.spins[beta]
            )

        for alpha, beta in self.J_32:
            gradient[alpha] += 3 * (
                np.einsum(
                    "tju,j,u->t",
                    self.J_32[(alpha, beta)],
                    spin_directions[alpha],
                    spin_directions[beta],
                )
                * self.spins[alpha] ** 2
                * self.spins[beta]
            )

        for alpha, beta in self.J_421:
            gradient[alpha] += 4 * (
                np.einsum(
                    "tjuv,j,u,v->t",
                    self.J_421[(alpha, beta)],
                    spin_directions[alpha],
                    spin_directions[alpha],
                    spin_directions[beta],
                )
                * self.spins[alpha] ** 3
                * self.spins[beta]
            )

        for alpha, beta in self.J_422:
            gradient[alpha] += 4 * (
                np.einsum(
                    "tjuv,j,u,v->t",
                    self.J_422[(alpha, beta)],
                    spin_directions[alpha],
                    spin_directions[beta],
                    spin_directions[beta],
                )
                * self.spins[alpha] ** 2
                * self.spins[beta] ** 2
            )

        for alpha, beta, gamma in self.J_33:
            gradient[alpha] += 3 * (
                np.einsum(
                    "tju,j,u->t",
                    self.J_33[(alpha, beta, gamma)],
                    spin_directions[beta],
                    spin_directions[gamma],
                )
                * self.spins[alpha]
                * self.spins[beta]
                * self.spins[gamma]
            )

        for alpha, beta, gamma in self.J_43:
            gradient[alpha] += 4 * (
                np.einsum(
                    "tjuv,j,u,v->t",
                    self.J_43[(alpha, beta, gamma)],
                    spin_directions[alpha],
                    spin_directions[beta],
                    spin_directions[gamma],
                )
                * self.spins[alpha] ** 2
                * self.spins[beta]
                * self.spins[gamma]
            )

        for alpha, beta, gamma, epsilon in self.J_44:
            gradient[alpha] += 4 * (
                np.einsum(
                    "tjuv,j,u,v->t",
                    self.J_44[(alpha, beta, gamma, epsilon)],
                    spin_directions[beta],
                    spin_directions[gamma],
                    spin_directions[epsilon],
                )
                * self.spins[alpha]
                * self.spins[beta]
                * self.spins[gamma]
                * self.spins[epsilon]
            )

        return gradient

    def torque(self, spin_directions, _normalize=True):
        r"""
        Computes torque on each spin.

        Parameters
        ----------
        spin_directions : (M, 3) |array-like|_
            Directions of spin vectors. Only directions of vectors are used,
            modulus is ignored. ``M`` is the amount of magnetic atoms in the
            Hamiltonian. The order of spin directions is the same as the order
            of magnetic atoms in ``spinham.magnetic_atoms.spins``.
        _normalize : bool, default True
            Whether to normalize the spin_directions or use the provided vectors as is.
            This parameter is technical and we do not recommend to use it at all.

        Returns
        -------
        torque : (M, 3) :numpy:`ndarray`

            .. code-block:: python

                [[t1x, t1y, t1z], [t2x, t2y, t2z], ...[tMx, tMy, tMz]]
        """
        return np.cross(
            spin_directions,
            self.gradient(spin_directions=spin_directions, _normalize=_normalize),
        )

    def _zoom(
        self,
        reference_sd,
        search_direction,
        phi_0,
        der_0,
        alpha_lo,
        alpha_hi,
        c1=_C1,
        c2=_C2,
    ):
        sd_lo = _rotate_sd(
            reference_sd=reference_sd, rotation=alpha_lo * search_direction
        )
        sd_hi = _rotate_sd(
            reference_sd=reference_sd, rotation=alpha_hi * search_direction
        )

        phi_lo = self.E_0(spin_directions=sd_lo)
        phi_hi = self.E_0(spin_directions=sd_hi)

        der_lo = self.torque(spin_directions=sd_lo).flatten() @ search_direction
        der_hi = self.torque(spin_directions=sd_hi).flatten() @ search_direction

        trial_steps = 0
        phi_min = None
        while True:
            # Interpolate
            alpha_j = _cubic_interpolation(
                alpha_l=alpha_lo,
                alpha_h=alpha_hi,
                phi_l=phi_lo,
                phi_h=phi_hi,
                der_l=der_lo,
                der_h=der_hi,
            )
            # Evaluate \phi(\alpha_i)
            sd_j = _rotate_sd(
                reference_sd=reference_sd, rotation=alpha_j * search_direction
            )
            phi_j = self.E_0(spin_directions=sd_j)

            # Safeguard
            if phi_min is None:
                phi_min = phi_j
            else:
                if phi_j < phi_min:
                    phi_min = phi_j
                    trial_steps = 0
            trial_steps += 1
            if trial_steps > 10:
                return alpha_j

            # Evaluate \phi^{\prime}(\alpha_i)
            der_j = self.torque(spin_directions=sd_j).flatten() @ search_direction

            if phi_j > phi_0 + c1 * alpha_j * der_0 or phi_j >= phi_lo:
                alpha_hi = alpha_j
                phi_hi = phi_j
                der_hi = der_j
            else:
                if abs(der_j) <= -c2 * der_0:
                    return alpha_j
                if der_j * (alpha_hi - alpha_lo) >= 0:
                    alpha_hi = alpha_lo
                    phi_hi = phi_lo
                    der_hi = der_lo
                alpha_lo = alpha_j
                phi_lo = phi_j
                der_lo = der_j

    def _line_search(
        self,
        reference_sd,
        search_direction,
        phi_0,
        der_0,
        c1=_C1,
        c2=_C2,
        alpha_max=2.0,
        max_iterations=10000,
    ):
        # First check if step alpha=1 is good to go:
        sd_1 = _rotate_sd(reference_sd=reference_sd, rotation=search_direction)
        phi_1 = self.E_0(spin_directions=sd_1)
        der_1 = self.torque(spin_directions=sd_1).flatten() @ search_direction

        if phi_1 <= phi_0 + c1 * der_0 and abs(der_1) <= c2 * abs(der_0):
            return 1.0

        # If not, then proceed to the algorithm

        alpha_prev = 0
        phi_prev = phi_0
        der_prev = der_0

        sd_max = _rotate_sd(
            reference_sd=reference_sd, rotation=alpha_max * search_direction
        )
        phi_max = self.E_0(spin_directions=sd_max)
        der_max = self.torque(spin_directions=sd_max).flatten() @ search_direction

        alpha_i = _cubic_interpolation(
            alpha_l=alpha_prev,
            alpha_h=alpha_max,
            phi_l=phi_prev,
            phi_h=phi_max,
            der_l=der_prev,
            der_h=der_max,
        )

        for i in range(1, max_iterations):
            # Evaluate \phi(\alpha_i)
            sd_i = _rotate_sd(
                reference_sd=reference_sd, rotation=alpha_i * search_direction
            )
            phi_i = self.E_0(spin_directions=sd_i)

            if phi_i > phi_0 + c1 * alpha_i * der_0 or (i > 1 and phi_i >= phi_prev):
                return self._zoom(
                    reference_sd=reference_sd,
                    search_direction=search_direction,
                    phi_0=phi_0,
                    der_0=der_0,
                    alpha_lo=alpha_prev,
                    alpha_hi=alpha_i,
                )

            # Evaluate \phi^{\prime}(\alpha_i)
            der_i = self.torque(spin_directions=sd_i).flatten() @ search_direction

            if abs(der_i) <= -c2 * der_0:
                return alpha_i

            if der_i >= 0:
                return self._zoom(
                    reference_sd=reference_sd,
                    search_direction=search_direction,
                    phi_0=phi_0,
                    der_0=der_0,
                    alpha_lo=alpha_i,
                    alpha_hi=alpha_prev,
                )

            # Choose alpha_{i+1}
            alpha_next = _cubic_interpolation(
                alpha_l=alpha_i,
                alpha_h=alpha_max,
                phi_l=phi_i,
                phi_h=phi_max,
                der_l=der_i,
                der_h=der_max,
            )

            # Set i <- i+1
            alpha_prev = alpha_i
            phi_prev = phi_i
            der_prev = der_i

            alpha_i = alpha_next

        raise ValueError(
            f"Line search did not converge in {max_iterations} iterations."
        )

    def optimize(
        self,
        initial_guess=None,
        energy_tolerance=1e-5,
        torque_tolerance=1e-5,
        quiet=False,
    ):
        r"""
        Optimize classical energy by varying the directions of spins in the unit cell.

        Parameters
        ----------
        initial_guess : (M, 3) or (3,) |array-like|_, optional
            Initial guess for the direction of the spin vectors.
        energy_tolerance : float, default 1e-5
            Energy tolerance for the two consecutive steps of the optimization.
        torque_tolerance : float, default 1e-5
            Torque tolerance for the two consecutive steps of the optimization.
        quiet : bool, default False
            Whether to suppress the output of the progress.

        Returns
        -------
        optimized_directions : (M, 3) :numpy:`ndarray`
            Optimized direction of the spin vectors.

        See Also
        --------
        optimize_generator
        """

        if initial_guess is None:
            initial_guess = np.random.uniform(low=-1, high=1, size=(self.M, 3))

        sd_k = initial_guess / np.linalg.norm(initial_guess, axis=1)[:, np.newaxis]

        if not quiet:
            n_energy = max(-(int(log10(energy_tolerance)) - 2), 6)
            n_torque = max(-(int(log10(torque_tolerance)) - 2), 6)
            print(
                "─" * 5
                + "┬"
                + "─" * 13
                + "┬"
                + "─" * (6 + n_energy)
                + "┬"
                + "─" * (6 + n_torque)
            )
            print(
                f"{'step':^4} │ "
                f"{'E_0':^11} │ "
                f"{'delta E_0':^{n_energy + 4}} │ "
                f"{'max torque':^{n_torque + 4}}"
            )
            print(
                "─" * 5
                + "┴"
                + "─" * 13
                + "┴"
                + "─" * (6 + n_energy)
                + "┴"
                + "─" * (6 + n_torque)
            )

        tolerance = np.array([energy_tolerance, torque_tolerance], dtype=float)

        delta = 2 * tolerance

        hessinv_k = np.eye(3 * self.M, dtype=float)

        energy_k = self.E_0(spin_directions=sd_k)
        gradient_k = self.torque(spin_directions=sd_k).flatten()

        first_iteration = True
        step_counter = 1
        while (delta >= tolerance).any():
            search_direction = -hessinv_k @ gradient_k
            # print(f"search = {search_direction}")

            alpha_k = self._line_search(
                reference_sd=sd_k,
                search_direction=search_direction,
                phi_0=energy_k,
                der_0=gradient_k @ search_direction,
            )

            # alpha_k = max(alpha_k, 1e-3)
            # print(f"alpha_k = {alpha_k}")

            s_k = alpha_k * search_direction
            # print(f"s_k = {s_k}")

            sd_next = _rotate_sd(reference_sd=sd_k, rotation=s_k)
            # print(f"sd_next = {sd_next}")
            energy_next = self.E_0(spin_directions=sd_next)
            gradient_next = self.torque(spin_directions=sd_next).flatten()

            delta = np.array(
                [
                    abs(energy_next - energy_k),
                    # Pay attention to the np.reshape keywords
                    np.linalg.norm(
                        np.reshape(gradient_next, (self.M, 3)), axis=1
                    ).max(),
                ]
            )
            # print(f"deltas: {delta[0]:11.7f} {delta[1]:11.7f}")
            if not quiet:
                print(
                    f"{step_counter:<4}   "
                    f"{energy_next:>11.7f}   "
                    f"{delta[0]:>{n_energy + 4}.{n_energy}f}   "
                    f"{delta[1]:>{n_torque + 4}.{n_torque}f}"
                )

            if (delta < tolerance).all():
                break

            y_k = gradient_next - gradient_k
            # print(f"y_k = {y_k}")

            rho_k = 1 / (y_k @ s_k)
            # print(f"rho = {rho_k}")

            EYE = np.eye(3 * self.M)
            OUTER = np.outer(y_k, s_k)

            if first_iteration:
                first_iteration = False
                hessinv_k = (y_k @ s_k) / (y_k @ y_k) * hessinv_k

            hessinv_k = (EYE - rho_k * OUTER.T) @ hessinv_k @ (
                EYE - rho_k * OUTER
            ) + rho_k * s_k @ s_k

            sd_k = sd_next
            energy_k = energy_next
            gradient_k = gradient_next

            step_counter += 1
            # print("=" * 40)

        if not quiet:
            print("─" * (33 + n_energy + n_torque))
        return sd_next

    def optimize_generator(
        self,
        initial_guess=None,
        energy_tolerance=1e-5,
        torque_tolerance=1e-5,
    ):
        r"""
        Optimize classical energy by varying the directions of spins in the unit cell.

        .. versionadded:: 0.2.0

        Parameters
        ----------
        initial_guess : (M, 3) or (3,) |array-like|_, optional
            Initial guess for the direction of the spin vectors.
        energy_tolerance : float, default 1e-5
            Energy tolerance for the two consecutive steps of the optimization.
        torque_tolerance : float, default 1e-5
            Torque tolerance for the two consecutive steps of the optimization.

        Yields
        ------
        energy : float
            Classical energy of the iteration step
        gradient : (M, 3) :numpy:`ndarray`
            Gradient vectors for each spin of the iteration step.
        spin_directions : (M, 3) :numpy:`ndarray`
            Directions of the spin vectors of the iteration step.

        See Also
        --------
        optimize
        """

        if initial_guess is None:
            initial_guess = np.random.uniform(low=-1, high=1, size=(self.M, 3))

        sd_k = initial_guess / np.linalg.norm(initial_guess, axis=1)[:, np.newaxis]

        tolerance = np.array([energy_tolerance, torque_tolerance], dtype=float)

        delta = 2 * tolerance

        hessinv_k = np.eye(3 * self.M, dtype=float)

        energy_k = self.E_0(spin_directions=sd_k)
        gradient_k = self.torque(spin_directions=sd_k).flatten()

        first_iteration = True
        step_counter = 1

        yield (energy_k, gradient_k, sd_k)

        while (delta >= tolerance).any():
            search_direction = -hessinv_k @ gradient_k

            alpha_k = self._line_search(
                reference_sd=sd_k,
                search_direction=search_direction,
                phi_0=energy_k,
                der_0=gradient_k @ search_direction,
            )

            s_k = alpha_k * search_direction

            sd_next = _rotate_sd(reference_sd=sd_k, rotation=s_k)

            energy_next = self.E_0(spin_directions=sd_next)
            gradient_next = self.torque(spin_directions=sd_next).flatten()

            yield (energy_next, gradient_next, sd_next)

            delta = np.array(
                [
                    abs(energy_next - energy_k),
                    # Pay attention to the np.reshape keywords
                    np.linalg.norm(
                        np.reshape(gradient_next, (self.M, 3)), axis=1
                    ).max(),
                ]
            )

            if (delta < tolerance).all():
                break

            y_k = gradient_next - gradient_k

            rho_k = 1 / (y_k @ s_k)

            EYE = np.eye(3 * self.M)
            OUTER = np.outer(y_k, s_k)

            if first_iteration:
                first_iteration = False
                hessinv_k = (y_k @ s_k) / (y_k @ y_k) * hessinv_k

            hessinv_k = (EYE - rho_k * OUTER.T) @ hessinv_k @ (
                EYE - rho_k * OUTER
            ) + rho_k * s_k @ s_k

            sd_k = sd_next
            energy_k = energy_next
            gradient_k = gradient_next

            step_counter += 1

        return sd_next


# Populate __all__ with objects defined in this file
__all__ = list(set(dir()) - old_dir)
# Remove all semi-private objects
__all__ = [i for i in __all__ if not i.startswith("_")]
del old_dir
