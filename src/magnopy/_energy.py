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


from math import log10
import warnings

import numpy as np

from magnopy._data_validation import _validated_units
from magnopy._constants._units import _ENERGY_UNITS
from magnopy._parameters._interaction_parameters import (
    _InteractionParameters,
    _InteractionParametersIterator,
)
from magnopy._parameters._renormalization import _renormalized_parameters
from magnopy._local_rf import span_local_rfs
from magnopy._exceptions import ConventionError


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
    Classical energy of the spin Hamiltonian.

    This class is optimized for calculation of the energy for any spin directions for the
    given Hamiltonian.

    .. important::
        If the spin Hamiltonian is modified, then a new instance of the energy class
        should be created.

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

    Finally, ``energy`` can be used to compute classical energy of the Hamiltonian,
    its gradient, torque or search for the local minimum.

    .. doctest::

        >>> sd1 = [[1, 0, 0]]
        >>> sd2 = [[0, 1, 0]]
        >>> sd3 = [[0, 0, 1]]
        >>> # Default units are meV
        >>> energy(sd1), energy(sd2), energy(sd3)
        (-4.5, -4.5, -6.75)
        >>> # You can request other units
        >>> print(f"{energy(sd1, units='Joule'):.4e}")
        -7.2098e-22
    """

    def __init__(self, spinham):

        if len(spinham.p1) > 0 and spinham.convention._c1 is None:
            raise ConventionError(convention=spinham.convention, property="c1")

        if len(spinham.p21) > 0 and spinham.convention._c21 is None:
            raise ConventionError(convention=spinham.convention, property="c21")

        if len(spinham.p22) > 0 and spinham.convention._c22 is None:
            raise ConventionError(convention=spinham.convention, property="c22")

        if len(spinham.p31) > 0 and spinham.convention._c31 is None:
            raise ConventionError(convention=spinham.convention, property="c31")

        if len(spinham.p32) > 0 and spinham.convention._c32 is None:
            raise ConventionError(convention=spinham.convention, property="c32")

        if len(spinham.p33) > 0 and spinham.convention._c33 is None:
            raise ConventionError(convention=spinham.convention, property="c33")

        if len(spinham.p41) > 0 and spinham.convention._c41 is None:
            raise ConventionError(convention=spinham.convention, property="c41")

        if len(spinham.p42) > 0 and spinham.convention._c42 is None:
            raise ConventionError(convention=spinham.convention, property="c42")

        if len(spinham.p43) > 0 and spinham.convention._c43 is None:
            raise ConventionError(convention=spinham.convention, property="c43")

        if len(spinham.p44) > 0 and spinham.convention._c44 is None:
            raise ConventionError(convention=spinham.convention, property="c44")

        if len(spinham.p45) > 0 and spinham.convention._c45 is None:
            raise ConventionError(convention=spinham.convention, property="c45")

        if (
            len(spinham.parameters()) > 0
            and spinham.convention._spin_normalized is None
        ):
            raise ConventionError(
                convention=spinham.convention, property="spin_normalized"
            )

        if (
            len(spinham.parameters()) > 0
            and spinham.convention._multiple_counting is None
        ):
            raise ConventionError(
                convention=spinham.convention, property="multiple_counting"
            )

        initial_units = spinham.units
        initial_convention = spinham.convention

        self.convention = initial_convention.get_modified(
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

        spinham.units = "meV"
        spinham.convention = self.convention

        self.spins = np.array(spinham.magnetic_atoms.spins, dtype=float)
        self.M = spinham.M

        self._parameters = _InteractionParameters()

        for (n, p_n, _, alphas), parameter in spinham._parameters._container:
            n = len(alphas)
            alphas = tuple([spinham.map_to_magnetic[alpha] for alpha in alphas])
            self._parameters.add(
                specs=(n, p_n, tuple([(0, 0, 0)] * (n - 1)), alphas),
                parameter=parameter,
                when_present="sum",
            )

        spinham.units = initial_units
        spinham.convention = initial_convention

    def __call__(
        self, spin_directions, units="meV", _normalize=True, quantum_correction=False
    ) -> float:
        result = self.E_0(
            spin_directions=spin_directions, units=units, _normalize=_normalize
        )

        if quantum_correction:
            result += self.E_corr(
                spin_directions=spin_directions, units=units, _normalize=_normalize
            )
        return result

    def E_0(self, spin_directions, units="meV", _normalize=True) -> float:
        r"""
        Computes classical energy of the spin Hamiltonian.

        Parameters
        ----------

        spin_directions : (M, 3) |array-like|_
            Directions of spin vectors. Only directions of vectors are used, modulus is
            ignored. ``M`` is the amount of magnetic atoms in the Hamiltonian. The order
            of spin directions is the same as the order of magnetic atoms in
            ``spinham.magnetic_atoms.spins``.

        units : str, default "meV"
            .. versionadded:: 0.3.0

            Units of energy. See :ref:`user-guide_usage_units_energy` for the full
            list of supported units.


        _normalize : bool, default True
            Whether to normalize the spin_directions or use the provided vectors as is.
            This parameter is technical and we do not recommend using it at all.

        Returns
        -------

        E_0 : float
            Classic energy of state with ``spin_directions``. Returned in the units of
            ``units``.


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
            >>> # Default units are meV
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

        for specs, parameter in self._parameters._container:
            n, _, _, alphas = specs

            if n == 1:
                energy += parameter @ spins[alphas[0]]
            elif n == 2:
                energy += spins[alphas[0]] @ parameter @ spins[alphas[1]]
            elif n == 3:
                energy += np.einsum(
                    "iju,i,j,u",
                    parameter,
                    spins[alphas[0]],
                    spins[alphas[1]],
                    spins[alphas[2]],
                )
            elif n == 4:
                energy += np.einsum(
                    "ijuv,i,j,u,v",
                    parameter,
                    spins[alphas[0]],
                    spins[alphas[1]],
                    spins[alphas[2]],
                    spins[alphas[3]],
                )
            else:
                raise ValueError(f"Unsupported n={n} in energy calculation.")

        # Convert units if necessary
        if units != "meV":
            units = _validated_units(units=units, supported_units=_ENERGY_UNITS)
            energy = energy * _ENERGY_UNITS["mev"] / _ENERGY_UNITS[units]

        return float(energy)

    def E_corr(self, spin_directions, units="mev", _normalize=True) -> float:
        r"""
        Computes quantum correction to the classical energy of the spin Hamiltonian.

        See eq. S.26 in supplementary material of |paper-2026|_.

        Parameters
        ----------

        spin_directions : (M, 3) |array-like|_
            Directions of spin vectors. Only directions of vectors are used,
            modulus is ignored. ``M`` is the amount of magnetic atoms in the
            Hamiltonian. The order of spin directions is the same as the order
            of magnetic atoms in ``spinham.magnetic_atoms.spins``.

        units : str, default "meV"
            .. versionadded:: 0.3.0

            Units of energy. See :ref:`user-guide_usage_units_energy` for the full
            list of supported units.
        """
        x, y, z = span_local_rfs(
            directional_vectors=spin_directions, hybridize=False, _normalize=_normalize
        )
        p = x + 1j * y
        result = 0

        renormalized_parameters = _renormalized_parameters(
            parameters=self._parameters,
            convention=self.convention,
            spin_directions=z,
            spin_values=self.spins,
        )

        # S.36b
        for _, alphas, parameter in _InteractionParametersIterator(
            renormalized_parameters, n=2, p_n=1
        ):
            alpha_1 = alphas[0]
            result += (
                0.5
                * self.spins[alpha_1]
                * np.conjugate(p[alpha_1])
                @ parameter
                @ p[alpha_1]
            )

        # S.36c
        for _, alphas, parameter in _InteractionParametersIterator(
            renormalized_parameters, n=3, p_n=1
        ):
            alpha_1 = alphas[0]
            result -= (
                0.5
                * self.spins[alpha_1]
                * np.einsum(
                    "ijk,i,j,k",
                    parameter,
                    np.conjugate(p[alpha_1]),
                    z[alphas[1]],
                    p[alphas[2]],
                )
            )

        # S.36d (first term)
        for _, alphas, parameter in _InteractionParametersIterator(
            renormalized_parameters, n=4, p_n=1
        ):
            alpha_1 = alphas[0]

            middle_matrix = np.zeros((3, 3), dtype=complex)
            middle_matrix += np.einsum("i,j->ij", z[alpha_1], z[alpha_1])
            if round(2 * self.spins[alpha_1]) > 1:
                middle_matrix += (self.spins[alpha_1] - 0.5) * np.einsum(
                    "i,j->ij", np.conjugate(p[alpha_1]), p[alpha_1]
                )
            middle_matrix += (
                0.5
                * self.spins[alpha_1]
                * np.einsum("i,j->ij", p[alpha_1], np.conjugate(p[alpha_1]))
            )

            result += (
                0.5
                * self.spins[alpha_1]
                * np.einsum(
                    "ijkl,i,jk,l",
                    parameter,
                    np.conjugate(p[alpha_1]),
                    middle_matrix,
                    p[alpha_1],
                )
            )

        # S.36d (second term)
        for nus, alphas, parameter in _InteractionParametersIterator(
            renormalized_parameters, n=4, p_n=3
        ):
            # Symmetrization is used in the analytical formula,
            # so need to skip equivalent terms
            if nus[0] != (0, 0, 0):
                continue

            alpha_1 = alphas[0]
            alpha_2 = alphas[2]

            result += (
                0.75
                * self.spins[alpha_1]
                * self.spins[alpha_2]
                * np.einsum(
                    "ijkl,i,j,k,l",
                    parameter,
                    np.conjugate(p[alpha_1]),
                    p[alpha_1],
                    np.conjugate(p[alpha_2]),
                    p[alpha_2],
                )
            )

        # Convert units if necessary
        if units != "meV":
            units = _validated_units(units=units, supported_units=_ENERGY_UNITS)
            result = result * _ENERGY_UNITS["mev"] / _ENERGY_UNITS[units]

        return float(result.real)

    def gradient(
        self, spin_directions, units="meV", _normalize=True, quantum_correction=False
    ):
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

        units : str, default "meV"
            .. versionadded:: 0.3.0

            Units of energy. See :ref:`user-guide_usage_units_energy` for the full
            list of supported units.


        _normalize : bool, default True
            Whether to normalize the spin_directions or use the provided vectors as is.
            This parameter is technical and we do not recommend using it at all.

        quantum_correction : bool, default False
            Whether to include quantum correction to the energy in the optimization. If
            ``True``, then it optimizes :py:func:`.Energy.E_0` + :py:func:`.Energy.E_corr`.
            If ``False``, then it optimizes :py:func:`.Energy.E_0` only.

            .. warning::
                This option is experimental. Will be improved and tested in future
                releases. Use with caution.

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

        for _, alphas, parameter in _InteractionParametersIterator(
            _renormalized_parameters(
                parameters=self._parameters,
                convention=self.convention,
                spin_directions=spin_directions,
                spin_values=self.spins,
            ),
            n=1,
            p_n=1,
        ):
            gradient[alphas[0]] = parameter * self.spins[alphas[0]]

        if quantum_correction:
            sd = spin_directions.copy()
            h = 1e-6
            for alpha in range(self.M):
                for i in range(3):
                    sd[alpha][i] += h
                    energy_plus = self.E_corr(spin_directions=sd, _normalize=False)
                    sd[alpha][i] -= 2 * h
                    energy_minus = self.E_corr(spin_directions=sd, _normalize=False)
                    sd[alpha][i] += h

                    gradient[alpha][i] += (energy_plus - energy_minus) / 2 / h

        # Convert units if necessary
        if units != "meV":
            units = _validated_units(units=units, supported_units=_ENERGY_UNITS)
            gradient = gradient * _ENERGY_UNITS["mev"] / _ENERGY_UNITS[units]

        return gradient

    def torque(
        self, spin_directions, units="meV", _normalize=True, quantum_correction=False
    ):
        r"""
        Computes torque on each spin.

        Parameters
        ----------

        spin_directions : (M, 3) |array-like|_
            Directions of spin vectors. Only directions of vectors are used,
            modulus is ignored. ``M`` is the amount of magnetic atoms in the
            Hamiltonian. The order of spin directions is the same as the order
            of magnetic atoms in ``spinham.magnetic_atoms.spins``.

        units : str, default "meV"
            .. versionadded:: 0.3.0

            Units of energy. See :ref:`user-guide_usage_units_energy` for the full
            list of supported units.

        _normalize : bool, default True
            Whether to normalize the spin_directions or use the provided vectors as is.
            This parameter is technical and we do not recommend using it at all.

        quantum_correction : bool, default False
            Whether to include quantum correction to the energy in the optimization. If
            ``True``, then it optimizes :py:func:`.Energy.E_0` + :py:func:`.Energy.E_corr`.
            If ``False``, then it optimizes :py:func:`.Energy.E_0` only.

            .. warning::
                This option is experimental. Will be improved and tested in future
                releases. Use with caution.

        Returns
        -------

        torque : (M, 3) :numpy:`ndarray`

            .. code-block:: python

                [[t1x, t1y, t1z], [t2x, t2y, t2z], ...[tMx, tMy, tMz]]
        """
        return np.cross(
            spin_directions,
            self.gradient(
                spin_directions=spin_directions,
                units=units,
                _normalize=_normalize,
                quantum_correction=quantum_correction,
            ),
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
        quantum_correction=False,
    ):
        sd_lo = _rotate_sd(
            reference_sd=reference_sd, rotation=alpha_lo * search_direction
        )
        sd_hi = _rotate_sd(
            reference_sd=reference_sd, rotation=alpha_hi * search_direction
        )

        phi_lo = self(spin_directions=sd_lo, quantum_correction=quantum_correction)
        phi_hi = self(spin_directions=sd_hi, quantum_correction=quantum_correction)

        der_lo = (
            self.torque(
                spin_directions=sd_lo, quantum_correction=quantum_correction
            ).flatten()
            @ search_direction
        )
        der_hi = (
            self.torque(
                spin_directions=sd_hi, quantum_correction=quantum_correction
            ).flatten()
            @ search_direction
        )

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
            phi_j = self(spin_directions=sd_j, quantum_correction=quantum_correction)

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
            der_j = (
                self.torque(
                    spin_directions=sd_j, quantum_correction=quantum_correction
                ).flatten()
                @ search_direction
            )

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
        quantum_correction=False,
    ):
        # First check if step alpha=1 is good to go:
        sd_1 = _rotate_sd(reference_sd=reference_sd, rotation=search_direction)
        phi_1 = self(spin_directions=sd_1, quantum_correction=quantum_correction)
        der_1 = (
            self.torque(
                spin_directions=sd_1, quantum_correction=quantum_correction
            ).flatten()
            @ search_direction
        )

        if phi_1 <= phi_0 + c1 * der_0 and abs(der_1) <= c2 * abs(der_0):
            return 1.0

        # If not, then proceed to the algorithm

        alpha_prev = 0
        phi_prev = phi_0
        der_prev = der_0

        sd_max = _rotate_sd(
            reference_sd=reference_sd, rotation=alpha_max * search_direction
        )
        phi_max = self(spin_directions=sd_max, quantum_correction=quantum_correction)
        der_max = (
            self.torque(
                spin_directions=sd_max, quantum_correction=quantum_correction
            ).flatten()
            @ search_direction
        )

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
            phi_i = self(spin_directions=sd_i, quantum_correction=quantum_correction)

            if phi_i > phi_0 + c1 * alpha_i * der_0 or (i > 1 and phi_i >= phi_prev):
                return self._zoom(
                    reference_sd=reference_sd,
                    search_direction=search_direction,
                    phi_0=phi_0,
                    der_0=der_0,
                    alpha_lo=alpha_prev,
                    alpha_hi=alpha_i,
                    quantum_correction=quantum_correction,
                )

            # Evaluate \phi^{\prime}(\alpha_i)
            der_i = (
                self.torque(
                    spin_directions=sd_i, quantum_correction=quantum_correction
                ).flatten()
                @ search_direction
            )

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
                    quantum_correction=quantum_correction,
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
        quantum_correction=False,
    ):
        r"""
        Optimizes classical energy by varying the directions of spins in the unit cell.

        Parameters
        ----------

        initial_guess : (M, 3) or (3,) |array-like|_, optional
            Initial guess for the direction of the spin vectors.

        energy_tolerance : float, default 1e-5
            Energy tolerance for the two consecutive steps of the optimization. In the units of meV.

        torque_tolerance : float, default 1e-5
            Torque tolerance for the two consecutive steps of the optimization. In the units of meV.

        quiet : bool, default False
            Whether to suppress the output of the progress.

        quantum_correction : bool, default False
            Whether to include quantum correction to the energy in the optimization. If
            ``True``, then it optimizes :py:func:`.Energy.E_0` + :py:func:`.Energy.E_corr`.
            If ``False``, then it optimizes :py:func:`.Energy.E_0` only.

            .. warning::
                This option is experimental. Will be improved and tested in future
                releases. Use with caution.

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
            n_energy = max(-(int(log10(energy_tolerance)) - 2), 12)
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
                f"{'Energy':^11} │ "
                f"{'delta Energy':^{n_energy + 4}} │ "
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

        energy_k = self(spin_directions=sd_k, quantum_correction=quantum_correction)
        gradient_k = self.torque(
            spin_directions=sd_k, quantum_correction=quantum_correction
        ).flatten()

        first_iteration = True
        step_counter = 1
        # Curvature failure guard to avoid silent looping in degenerate manifolds
        curv_fail_run = 0
        max_curv_fails = 10

        while (delta >= tolerance).any():
            search_direction = -hessinv_k @ gradient_k

            alpha_k = self._line_search(
                reference_sd=sd_k,
                search_direction=search_direction,
                phi_0=energy_k,
                der_0=gradient_k @ search_direction,
                quantum_correction=quantum_correction,
            )

            s_k = alpha_k * search_direction

            sd_next = _rotate_sd(reference_sd=sd_k, rotation=s_k)

            energy_next = self(
                spin_directions=sd_next, quantum_correction=quantum_correction
            )
            gradient_next = self.torque(
                spin_directions=sd_next, quantum_correction=quantum_correction
            ).flatten()

            delta = np.array(
                [
                    abs(energy_next - energy_k),
                    # Pay attention to the np.reshape keywords
                    np.linalg.norm(
                        np.reshape(gradient_next, (self.M, 3)), axis=1
                    ).max(),
                ]
            )

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
            # Curvature safeguard: avoid divide-by-zero / non-finite ys
            ys = float(y_k @ s_k)
            EYE = np.eye(hessinv_k.shape[0], dtype=float)

            if (not np.isfinite(ys)) or (abs(ys) < 1e-12):
                # Degenerate case: skip the inverse-BFGS update (or reset H if preferred)
                curv_fail_run += 1
                if curv_fail_run > max_curv_fails:
                    raise RuntimeError(
                        f"BFGS curvature failure repeated {curv_fail_run} times: "
                        f"s^T y={ys}, ||s||={np.linalg.norm(s_k)}, ||y||={np.linalg.norm(y_k)}. "
                        "Try a different initial guess (i.e. re-run magnopy-optimize-sd or energy.optimize())."
                    )
                # Reset hessian:
                hessinv_k = EYE
                warnings.warn(
                    f"BFGS curvature failure repeated {curv_fail_run} times: s^T y={ys}, ||s||={np.linalg.norm(s_k)}, ||y||={np.linalg.norm(y_k)}. Hessian was reset.",
                    RuntimeWarning,
                    stacklevel=2,
                )
            else:
                curv_fail_run = 0
                rho_k = 1.0 / ys
                OUTER = np.outer(y_k, s_k)

                # Safe initial scaling of H^{-1}
                if first_iteration:
                    first_iteration = False
                    denom = float(y_k @ y_k)
                    if np.isfinite(denom) and denom > 0.0:
                        hessinv_k = (ys / denom) * hessinv_k

                # Stable inverse-BFGS update
                hessinv_k = (EYE - rho_k * OUTER.T) @ hessinv_k @ (
                    EYE - rho_k * OUTER
                ) + rho_k * np.outer(s_k, s_k)

            sd_k = sd_next
            energy_k = energy_next
            gradient_k = gradient_next

            step_counter += 1
            # print("=" * 40)

        if not quiet:
            print("─" * (33 + n_energy + n_torque))
        return sd_next


# Populate __all__ with objects defined in this file
__all__ = list(set(dir()) - old_dir)
# Remove all semi-private objects
__all__ = [i for i in __all__ if not i.startswith("_")]
del old_dir
