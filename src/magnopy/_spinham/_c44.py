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

from magnopy._spinham._validators import (
    _spins_ordered,
    _validate_atom_index,
    _validate_unit_cell_index,
)


def _get_primary_p44(alpha, beta, gamma, epsilon, nu, _lambda, rho, parameter=None):
    r"""
    Return the primary version of the parameter.

    For the definition of the primary version see
    :ref:`user-guide_theory-behind_multiple-counting`.

    Parameters
    ----------
    alpha : int
        Index of the first atom.
    beta : int
        Index of the second atom.
    gamma : int
        Index of the third atom.
    epsilon : int
        Index of the fourth atom.
    nu : tuple of 3 int
        Unit cell for the second atom.
    _lambda : tuple of 3 int
        Unit cell for the third atom.
    rho : tuple of 3 int
        Unit cell for the fourth atom.
    parameter : (3, 3, 3, 3) :numpy:`ndarray`, optional
        Full matrix of the parameter.

    Returns
    -------
    alpha : int
        Index of the first atom.
    beta : int
        Index of the second atom.
    gamma : int
        Index of the third atom.
    epsilon : int
        Index of the fourth atom.
    nu : tuple of 3 int
        Unit cell for the second atom.
    _lambda : tuple of 3 int
        Unit cell for the third atom.
    rho : tuple of 3 int
        Unit cell for the fourth atom.
    parameter : (3, 3, 3, 3) :numpy:`ndarray`, optional
        Full matrix of the parameter. It is returned only if ``parameter is not None``.
    """

    def _ordered(mu1, alpha1, mu2, alpha2, mu3, alpha3, mu4, alpha4):
        return (
            _spins_ordered(mu1=mu1, alpha1=alpha1, mu2=mu2, alpha2=alpha2)
            and _spins_ordered(mu1=mu2, alpha1=alpha2, mu2=mu3, alpha2=alpha3)
            and _spins_ordered(mu1=mu3, alpha1=alpha3, mu2=mu4, alpha2=alpha4)
        )

    # Case 1
    if _ordered(
        mu1=(0, 0, 0),
        alpha1=alpha,
        mu2=nu,
        alpha2=beta,
        mu3=_lambda,
        alpha3=gamma,
        mu4=rho,
        alpha4=epsilon,
    ):
        pass
    # Case 2
    elif _ordered(
        mu1=(0, 0, 0),
        alpha1=alpha,
        mu2=nu,
        alpha2=beta,
        mu3=rho,
        alpha3=epsilon,
        mu4=_lambda,
        alpha4=gamma,
    ):
        alpha, beta, gamma, epsilon = alpha, beta, epsilon, gamma
        nu, _lambda, rho = nu, rho, _lambda
        if parameter is not None:
            parameter = np.transpose(parameter, (0, 1, 3, 2))
    # Case 3
    elif _ordered(
        mu1=(0, 0, 0),
        alpha1=alpha,
        mu2=_lambda,
        alpha2=gamma,
        mu3=nu,
        alpha3=beta,
        mu4=rho,
        alpha4=epsilon,
    ):
        alpha, beta, gamma, epsilon = alpha, gamma, beta, epsilon
        nu, _lambda, rho = _lambda, nu, rho
        if parameter is not None:
            parameter = np.transpose(parameter, (0, 2, 1, 3))
    # Case 4
    elif _ordered(
        mu1=(0, 0, 0),
        alpha1=alpha,
        mu2=_lambda,
        alpha2=gamma,
        mu3=rho,
        alpha3=epsilon,
        mu4=nu,
        alpha4=beta,
    ):
        alpha, beta, gamma, epsilon = alpha, gamma, epsilon, beta
        nu, _lambda, rho = _lambda, rho, nu
        if parameter is not None:
            parameter = np.transpose(parameter, (0, 3, 1, 2))
    # Case 5
    elif _ordered(
        mu1=(0, 0, 0),
        alpha1=alpha,
        mu2=rho,
        alpha2=epsilon,
        mu3=nu,
        alpha3=beta,
        mu4=_lambda,
        alpha4=gamma,
    ):
        alpha, beta, gamma, epsilon = alpha, epsilon, beta, gamma
        nu, _lambda, rho = rho, nu, _lambda
        if parameter is not None:
            parameter = np.transpose(parameter, (0, 2, 3, 1))
    # Case 6
    elif _ordered(
        mu1=(0, 0, 0),
        alpha1=alpha,
        mu2=rho,
        alpha2=epsilon,
        mu3=_lambda,
        alpha3=gamma,
        mu4=nu,
        alpha4=beta,
    ):
        alpha, beta, gamma, epsilon = alpha, epsilon, gamma, beta
        nu, _lambda, rho = rho, _lambda, nu
        if parameter is not None:
            parameter = np.transpose(parameter, (0, 3, 2, 1))
    # Case 7
    elif _ordered(
        mu1=nu,
        alpha1=beta,
        mu2=(0, 0, 0),
        alpha2=alpha,
        mu3=_lambda,
        alpha3=gamma,
        mu4=rho,
        alpha4=epsilon,
    ):
        alpha, beta, gamma, epsilon = beta, alpha, gamma, epsilon
        nu1, nu2, nu3 = nu
        lambda1, lambda2, lambda3 = _lambda
        rho1, rho2, rho3 = rho
        nu = (-nu1, -nu2, -nu3)
        _lambda = (lambda1 - nu1, lambda2 - nu2, lambda3 - nu3)
        rho = (rho1 - nu1, rho2 - nu2, rho3 - nu3)
        if parameter is not None:
            parameter = np.transpose(parameter, (1, 0, 2, 3))
    # Case 8
    elif _ordered(
        mu1=nu,
        alpha1=beta,
        mu2=(0, 0, 0),
        alpha2=alpha,
        mu3=rho,
        alpha3=epsilon,
        mu4=_lambda,
        alpha4=gamma,
    ):
        alpha, beta, gamma, epsilon = beta, alpha, epsilon, gamma
        nu1, nu2, nu3 = nu
        lambda1, lambda2, lambda3 = _lambda
        rho1, rho2, rho3 = rho
        nu = (-nu1, -nu2, -nu3)
        _lambda = (rho1 - nu1, rho2 - nu2, rho3 - nu3)
        rho = (lambda1 - nu1, lambda2 - nu2, lambda3 - nu3)
        if parameter is not None:
            parameter = np.transpose(parameter, (1, 0, 3, 2))
    # Case 9
    elif _ordered(
        mu1=nu,
        alpha1=beta,
        mu2=_lambda,
        alpha2=gamma,
        mu3=(0, 0, 0),
        alpha3=alpha,
        mu4=rho,
        alpha4=epsilon,
    ):
        alpha, beta, gamma, epsilon = beta, gamma, alpha, epsilon
        nu1, nu2, nu3 = nu
        lambda1, lambda2, lambda3 = _lambda
        rho1, rho2, rho3 = rho
        nu = (lambda1 - nu1, lambda2 - nu2, lambda3 - nu3)
        _lambda = (-nu1, -nu2, -nu3)
        rho = (rho1 - nu1, rho2 - nu2, rho3 - nu3)
        if parameter is not None:
            parameter = np.transpose(parameter, (2, 0, 1, 3))
    # Case 10
    elif _ordered(
        mu1=nu,
        alpha1=beta,
        mu2=_lambda,
        alpha2=gamma,
        mu3=rho,
        alpha3=epsilon,
        mu4=(0, 0, 0),
        alpha4=alpha,
    ):
        alpha, beta, gamma, epsilon = beta, gamma, epsilon, alpha
        nu1, nu2, nu3 = nu
        lambda1, lambda2, lambda3 = _lambda
        rho1, rho2, rho3 = rho
        nu = (lambda1 - nu1, lambda2 - nu2, lambda3 - nu3)
        _lambda = (rho1 - nu1, rho2 - nu2, rho3 - nu3)
        rho = (-nu1, -nu2, -nu3)
        if parameter is not None:
            parameter = np.transpose(parameter, (3, 0, 1, 2))
    # Case 11
    elif _ordered(
        mu1=nu,
        alpha1=beta,
        mu2=rho,
        alpha2=epsilon,
        mu3=(0, 0, 0),
        alpha3=alpha,
        mu4=_lambda,
        alpha4=gamma,
    ):
        alpha, beta, gamma, epsilon = beta, epsilon, alpha, gamma
        nu1, nu2, nu3 = nu
        lambda1, lambda2, lambda3 = _lambda
        rho1, rho2, rho3 = rho
        nu = (rho1 - nu1, rho2 - nu2, rho3 - nu3)
        _lambda = (-nu1, -nu2, -nu3)
        rho = (lambda1 - nu1, lambda2 - nu2, lambda3 - nu3)
        if parameter is not None:
            parameter = np.transpose(parameter, (2, 0, 3, 1))
    # Case 12
    elif _ordered(
        mu1=nu,
        alpha1=beta,
        mu2=rho,
        alpha2=epsilon,
        mu3=_lambda,
        alpha3=gamma,
        mu4=(0, 0, 0),
        alpha4=alpha,
    ):
        alpha, beta, gamma, epsilon = beta, epsilon, gamma, alpha
        nu1, nu2, nu3 = nu
        lambda1, lambda2, lambda3 = _lambda
        rho1, rho2, rho3 = rho
        nu = (rho1 - nu1, rho2 - nu2, rho3 - nu3)
        _lambda = (lambda1 - nu1, lambda2 - nu2, lambda3 - nu3)
        rho = (-nu1, -nu2, -nu3)
        if parameter is not None:
            parameter = np.transpose(parameter, (3, 0, 2, 1))
    # Case 13
    elif _ordered(
        mu1=_lambda,
        alpha1=gamma,
        mu2=(0, 0, 0),
        alpha2=alpha,
        mu3=nu,
        alpha3=beta,
        mu4=rho,
        alpha4=epsilon,
    ):
        alpha, beta, gamma, epsilon = gamma, alpha, beta, epsilon
        nu1, nu2, nu3 = nu
        lambda1, lambda2, lambda3 = _lambda
        rho1, rho2, rho3 = rho
        nu = (-lambda1, -lambda2, -lambda3)
        _lambda = (nu1 - lambda1, nu2 - lambda2, nu3 - lambda3)
        rho = (rho1 - lambda1, rho2 - lambda2, rho3 - lambda3)
        if parameter is not None:
            parameter = np.transpose(parameter, (1, 2, 0, 3))
    # Case 14
    elif _ordered(
        mu1=_lambda,
        alpha1=gamma,
        mu2=(0, 0, 0),
        alpha2=alpha,
        mu3=rho,
        alpha3=epsilon,
        mu4=nu,
        alpha4=beta,
    ):
        alpha, beta, gamma, epsilon = gamma, alpha, epsilon, beta
        nu1, nu2, nu3 = nu
        lambda1, lambda2, lambda3 = _lambda
        rho1, rho2, rho3 = rho
        nu = (-lambda1, -lambda2, -lambda3)
        _lambda = (rho1 - lambda1, rho2 - lambda2, rho3 - lambda3)
        rho = (nu1 - lambda1, nu2 - lambda2, nu3 - lambda3)
        if parameter is not None:
            parameter = np.transpose(parameter, (1, 3, 0, 2))
    # Case 15
    elif _ordered(
        mu1=_lambda,
        alpha1=gamma,
        mu2=nu,
        alpha2=beta,
        mu3=(0, 0, 0),
        alpha3=alpha,
        mu4=rho,
        alpha4=epsilon,
    ):
        alpha, beta, gamma, epsilon = gamma, beta, alpha, epsilon
        nu1, nu2, nu3 = nu
        lambda1, lambda2, lambda3 = _lambda
        rho1, rho2, rho3 = rho
        nu = (nu1 - lambda1, nu2 - lambda2, nu3 - lambda3)
        _lambda = (-lambda1, -lambda2, -lambda3)
        rho = (rho1 - lambda1, rho2 - lambda2, rho3 - lambda3)
        if parameter is not None:
            parameter = np.transpose(parameter, (2, 1, 0, 3))
    # Case 16
    elif _ordered(
        mu1=_lambda,
        alpha1=gamma,
        mu2=nu,
        alpha2=beta,
        mu3=rho,
        alpha3=epsilon,
        mu4=(0, 0, 0),
        alpha4=alpha,
    ):
        alpha, beta, gamma, epsilon = gamma, beta, epsilon, alpha
        nu1, nu2, nu3 = nu
        lambda1, lambda2, lambda3 = _lambda
        rho1, rho2, rho3 = rho
        nu = (nu1 - lambda1, nu2 - lambda2, nu3 - lambda3)
        _lambda = (rho1 - lambda1, rho2 - lambda2, rho3 - lambda3)
        rho = (-lambda1, -lambda2, -lambda3)
        if parameter is not None:
            parameter = np.transpose(parameter, (3, 1, 0, 2))
    # Case 17
    elif _ordered(
        mu1=_lambda,
        alpha1=gamma,
        mu2=rho,
        alpha2=epsilon,
        mu3=(0, 0, 0),
        alpha3=alpha,
        mu4=nu,
        alpha4=beta,
    ):
        alpha, beta, gamma, epsilon = gamma, epsilon, alpha, beta
        nu1, nu2, nu3 = nu
        lambda1, lambda2, lambda3 = _lambda
        rho1, rho2, rho3 = rho
        nu = (rho1 - lambda1, rho2 - lambda2, rho3 - lambda3)
        _lambda = (-lambda1, -lambda2, -lambda3)
        rho = (nu1 - lambda1, nu2 - lambda2, nu3 - lambda3)
        if parameter is not None:
            parameter = np.transpose(parameter, (2, 3, 0, 1))
    # Case 18
    elif _ordered(
        mu1=_lambda,
        alpha1=gamma,
        mu2=rho,
        alpha2=epsilon,
        mu3=nu,
        alpha3=beta,
        mu4=(0, 0, 0),
        alpha4=alpha,
    ):
        alpha, beta, gamma, epsilon = gamma, epsilon, beta, alpha
        nu1, nu2, nu3 = nu
        lambda1, lambda2, lambda3 = _lambda
        rho1, rho2, rho3 = rho
        nu = (rho1 - lambda1, rho2 - lambda2, rho3 - lambda3)
        _lambda = (nu1 - lambda1, nu2 - lambda2, nu3 - lambda3)
        rho = (-lambda1, -lambda2, -lambda3)
        if parameter is not None:
            parameter = np.transpose(parameter, (3, 2, 0, 1))
    # Case 19
    elif _ordered(
        mu1=rho,
        alpha1=epsilon,
        mu2=(0, 0, 0),
        alpha2=alpha,
        mu3=nu,
        alpha3=beta,
        mu4=_lambda,
        alpha4=gamma,
    ):
        alpha, beta, gamma, epsilon = epsilon, alpha, beta, gamma
        nu1, nu2, nu3 = nu
        lambda1, lambda2, lambda3 = _lambda
        rho1, rho2, rho3 = rho
        nu = (-rho1, -rho2, -rho3)
        _lambda = (nu1 - rho1, nu2 - rho2, nu3 - rho3)
        rho = (lambda1 - rho1, lambda2 - rho2, lambda3 - rho3)
        if parameter is not None:
            parameter = np.transpose(parameter, (1, 2, 3, 0))
    # Case 20
    elif _ordered(
        mu1=rho,
        alpha1=epsilon,
        mu2=(0, 0, 0),
        alpha2=alpha,
        mu3=_lambda,
        alpha3=gamma,
        mu4=nu,
        alpha4=beta,
    ):
        alpha, beta, gamma, epsilon = epsilon, alpha, gamma, beta
        nu1, nu2, nu3 = nu
        lambda1, lambda2, lambda3 = _lambda
        rho1, rho2, rho3 = rho
        nu = (-rho1, -rho2, -rho3)
        _lambda = (lambda1 - rho1, lambda2 - rho2, lambda3 - rho3)
        rho = (nu1 - rho1, nu2 - rho2, nu3 - rho3)
        if parameter is not None:
            parameter = np.transpose(parameter, (1, 3, 2, 0))
    # Case 21
    elif _ordered(
        mu1=rho,
        alpha1=epsilon,
        mu2=nu,
        alpha2=beta,
        mu3=(0, 0, 0),
        alpha3=alpha,
        mu4=_lambda,
        alpha4=gamma,
    ):
        alpha, beta, gamma, epsilon = epsilon, beta, alpha, gamma
        nu1, nu2, nu3 = nu
        lambda1, lambda2, lambda3 = _lambda
        rho1, rho2, rho3 = rho
        nu = (nu1 - rho1, nu2 - rho2, nu3 - rho3)
        _lambda = (-rho1, -rho2, -rho3)
        rho = (lambda1 - rho1, lambda2 - rho2, lambda3 - rho3)
        if parameter is not None:
            parameter = np.transpose(parameter, (2, 1, 3, 0))
    # Case 22
    elif _ordered(
        mu1=rho,
        alpha1=epsilon,
        mu2=nu,
        alpha2=beta,
        mu3=_lambda,
        alpha3=gamma,
        mu4=(0, 0, 0),
        alpha4=alpha,
    ):
        alpha, beta, gamma, epsilon = epsilon, beta, gamma, alpha
        nu1, nu2, nu3 = nu
        lambda1, lambda2, lambda3 = _lambda
        rho1, rho2, rho3 = rho
        nu = (nu1 - rho1, nu2 - rho2, nu3 - rho3)
        _lambda = (lambda1 - rho1, lambda2 - rho2, lambda3 - rho3)
        rho = (-rho1, -rho2, -rho3)
        if parameter is not None:
            parameter = np.transpose(parameter, (3, 1, 2, 0))
    # Case 23
    elif _ordered(
        mu1=rho,
        alpha1=epsilon,
        mu2=_lambda,
        alpha2=gamma,
        mu3=(0, 0, 0),
        alpha3=alpha,
        mu4=nu,
        alpha4=beta,
    ):
        alpha, beta, gamma, epsilon = epsilon, beta, gamma, alpha
        nu1, nu2, nu3 = nu
        lambda1, lambda2, lambda3 = _lambda
        rho1, rho2, rho3 = rho
        nu = (lambda1 - rho1, lambda2 - rho2, lambda3 - rho3)
        _lambda = (-rho1, -rho2, -rho3)
        rho = (nu1 - rho1, nu2 - rho2, nu3 - rho3)
        if parameter is not None:
            parameter = np.transpose(parameter, (2, 3, 1, 0))
    # Case 24
    elif _ordered(
        mu1=rho,
        alpha1=epsilon,
        mu2=_lambda,
        alpha2=gamma,
        mu3=nu,
        alpha3=beta,
        mu4=(0, 0, 0),
        alpha4=alpha,
    ):
        alpha, beta, gamma, epsilon = epsilon, beta, gamma, alpha
        nu1, nu2, nu3 = nu
        lambda1, lambda2, lambda3 = _lambda
        rho1, rho2, rho3 = rho
        nu = (lambda1 - rho1, lambda2 - rho2, lambda3 - rho3)
        _lambda = (nu1 - rho1, nu2 - rho2, nu3 - rho3)
        rho = (-rho1, -rho2, -rho3)
        if parameter is not None:
            parameter = np.transpose(parameter, (3, 2, 1, 0))

    if parameter is None:
        return alpha, beta, gamma, epsilon, nu, _lambda, rho

    return alpha, beta, gamma, epsilon, nu, _lambda, rho, parameter


class _P44_iterator:
    R"""
    Iterator over the (four spins & four sites) parameters of the spin Hamiltonian.
    """

    def __init__(self, spinham) -> None:
        self.container = spinham._44
        self.mc = spinham.convention.multiple_counting
        self.length = len(self.container)
        self.index = 0

    def __next__(self):
        # Case 1
        if self.index < self.length:
            self.index += 1
            return self.container[self.index - 1]
        # Case 2
        elif self.mc and self.index < 2 * self.length:
            self.index += 1
            alpha, beta, gamma, epsilon, nu, _lambda, rho, parameter = self.container[
                self.index - 1 - self.length
            ]
            return [
                alpha,
                beta,
                epsilon,
                gamma,
                nu,
                rho,
                _lambda,
                np.transpose(parameter, (0, 1, 3, 2)),
            ]
        # Case 3
        elif self.mc and self.index < 3 * self.length:
            self.index += 1
            alpha, beta, gamma, epsilon, nu, _lambda, rho, parameter = self.container[
                self.index - 1 - 2 * self.length
            ]
            return [
                alpha,
                gamma,
                beta,
                epsilon,
                _lambda,
                nu,
                rho,
                np.transpose(parameter, (0, 2, 1, 3)),
            ]
        # Case 4
        elif self.mc and self.index < 4 * self.length:
            self.index += 1
            alpha, beta, gamma, epsilon, nu, _lambda, rho, parameter = self.container[
                self.index - 1 - 3 * self.length
            ]
            return [
                alpha,
                gamma,
                epsilon,
                beta,
                _lambda,
                rho,
                nu,
                np.transpose(parameter, (0, 3, 1, 2)),
            ]
        # Case 5
        elif self.mc and self.index < 5 * self.length:
            self.index += 1
            alpha, beta, gamma, epsilon, nu, _lambda, rho, parameter = self.container[
                self.index - 1 - 4 * self.length
            ]
            return [
                alpha,
                epsilon,
                beta,
                gamma,
                rho,
                nu,
                _lambda,
                np.transpose(parameter, (0, 2, 3, 1)),
            ]
        # Case 6
        elif self.mc and self.index < 6 * self.length:
            self.index += 1
            alpha, beta, gamma, epsilon, nu, _lambda, rho, parameter = self.container[
                self.index - 1 - 5 * self.length
            ]
            return [
                alpha,
                epsilon,
                gamma,
                beta,
                rho,
                _lambda,
                nu,
                np.transpose(parameter, (0, 3, 2, 1)),
            ]
        # Case 7
        elif self.mc and self.index < 7 * self.length:
            self.index += 1
            (
                alpha,
                beta,
                gamma,
                epsilon,
                (nu1, nu2, nu3),
                (lambda1, lambda2, lambda3),
                (rho1, rho2, rho3),
                parameter,
            ) = self.container[self.index - 1 - 6 * self.length]
            return [
                beta,
                alpha,
                gamma,
                epsilon,
                (-nu1, -nu2, -nu3),
                (lambda1 - nu1, lambda2 - nu2, lambda3 - nu3),
                (rho1 - nu1, rho2 - nu2, rho3 - nu3),
                np.transpose(parameter, (1, 0, 2, 3)),
            ]
        # Case 8
        elif self.mc and self.index < 8 * self.length:
            self.index += 1
            (
                alpha,
                beta,
                gamma,
                epsilon,
                (nu1, nu2, nu3),
                (lambda1, lambda2, lambda3),
                (rho1, rho2, rho3),
                parameter,
            ) = self.container[self.index - 1 - 7 * self.length]
            return [
                beta,
                alpha,
                epsilon,
                gamma,
                (-nu1, -nu2, -nu3),
                (rho1 - nu1, rho2 - nu2, rho3 - nu3),
                (lambda1 - nu1, lambda2 - nu2, lambda3 - nu3),
                np.transpose(parameter, (1, 0, 3, 2)),
            ]
        # Case 9
        elif self.mc and self.index < 9 * self.length:
            self.index += 1
            (
                alpha,
                beta,
                gamma,
                epsilon,
                (nu1, nu2, nu3),
                (lambda1, lambda2, lambda3),
                (rho1, rho2, rho3),
                parameter,
            ) = self.container[self.index - 1 - 8 * self.length]
            return [
                beta,
                gamma,
                alpha,
                epsilon,
                (lambda1 - nu1, lambda2 - nu2, lambda3 - nu3),
                (-nu1, -nu2, -nu3),
                (rho1 - nu1, rho2 - nu2, rho3 - nu3),
                np.transpose(parameter, (2, 0, 1, 3)),
            ]
        # Case 10
        elif self.mc and self.index < 10 * self.length:
            self.index += 1
            (
                alpha,
                beta,
                gamma,
                epsilon,
                (nu1, nu2, nu3),
                (lambda1, lambda2, lambda3),
                (rho1, rho2, rho3),
                parameter,
            ) = self.container[self.index - 1 - 9 * self.length]
            return [
                beta,
                gamma,
                epsilon,
                alpha,
                (lambda1 - nu1, lambda2 - nu2, lambda3 - nu3),
                (rho1 - nu1, rho2 - nu2, rho3 - nu3),
                (-nu1, -nu2, -nu3),
                np.transpose(parameter, (3, 0, 1, 2)),
            ]
        # Case 11
        elif self.mc and self.index < 11 * self.length:
            self.index += 1
            (
                alpha,
                beta,
                gamma,
                epsilon,
                (nu1, nu2, nu3),
                (lambda1, lambda2, lambda3),
                (rho1, rho2, rho3),
                parameter,
            ) = self.container[self.index - 1 - 10 * self.length]
            return [
                beta,
                epsilon,
                alpha,
                gamma,
                (rho1 - nu1, rho2 - nu2, rho3 - nu3),
                (-nu1, -nu2, -nu3),
                (lambda1 - nu1, lambda2 - nu2, lambda3 - nu3),
                np.transpose(parameter, (2, 0, 3, 1)),
            ]
        # Case 12
        elif self.mc and self.index < 12 * self.length:
            self.index += 1
            (
                alpha,
                beta,
                gamma,
                epsilon,
                (nu1, nu2, nu3),
                (lambda1, lambda2, lambda3),
                (rho1, rho2, rho3),
                parameter,
            ) = self.container[self.index - 1 - 11 * self.length]
            return [
                beta,
                epsilon,
                gamma,
                alpha,
                (rho1 - nu1, rho2 - nu2, rho3 - nu3),
                (lambda1 - nu1, lambda2 - nu2, lambda3 - nu3),
                (-nu1, -nu2, -nu3),
                np.transpose(parameter, (3, 0, 2, 1)),
            ]
        # Case 13
        elif self.mc and self.index < 13 * self.length:
            self.index += 1
            (
                alpha,
                beta,
                gamma,
                epsilon,
                (nu1, nu2, nu3),
                (lambda1, lambda2, lambda3),
                (rho1, rho2, rho3),
                parameter,
            ) = self.container[self.index - 1 - 12 * self.length]
            return [
                gamma,
                alpha,
                beta,
                epsilon,
                (-lambda1, -lambda2, -lambda3),
                (nu1 - lambda1, nu2 - lambda2, nu3 - lambda3),
                (rho1 - lambda1, rho2 - lambda2, rho3 - lambda3),
                np.transpose(parameter, (1, 2, 0, 3)),
            ]
        # Case 14
        elif self.mc and self.index < 14 * self.length:
            self.index += 1
            (
                alpha,
                beta,
                gamma,
                epsilon,
                (nu1, nu2, nu3),
                (lambda1, lambda2, lambda3),
                (rho1, rho2, rho3),
                parameter,
            ) = self.container[self.index - 1 - 13 * self.length]
            return [
                gamma,
                alpha,
                epsilon,
                beta,
                (-lambda1, -lambda2, -lambda3),
                (rho1 - lambda1, rho2 - lambda2, rho3 - lambda3),
                (nu1 - lambda1, nu2 - lambda2, nu3 - lambda3),
                np.transpose(parameter, (1, 3, 0, 2)),
            ]
        # Case 15
        elif self.mc and self.index < 15 * self.length:
            self.index += 1
            (
                alpha,
                beta,
                gamma,
                epsilon,
                (nu1, nu2, nu3),
                (lambda1, lambda2, lambda3),
                (rho1, rho2, rho3),
                parameter,
            ) = self.container[self.index - 1 - 14 * self.length]
            return [
                gamma,
                beta,
                alpha,
                epsilon,
                (nu1 - lambda1, nu2 - lambda2, nu3 - lambda3),
                (-lambda1, -lambda2, -lambda3),
                (rho1 - lambda1, rho2 - lambda2, rho3 - lambda3),
                np.transpose(parameter, (2, 1, 0, 3)),
            ]
        # Case 16
        elif self.mc and self.index < 16 * self.length:
            self.index += 1
            (
                alpha,
                beta,
                gamma,
                epsilon,
                (nu1, nu2, nu3),
                (lambda1, lambda2, lambda3),
                (rho1, rho2, rho3),
                parameter,
            ) = self.container[self.index - 1 - 15 * self.length]
            return [
                gamma,
                beta,
                epsilon,
                alpha,
                (nu1 - lambda1, nu2 - lambda2, nu3 - lambda3),
                (rho1 - lambda1, rho2 - lambda2, rho3 - lambda3),
                (-lambda1, -lambda2, -lambda3),
                np.transpose(parameter, (3, 1, 0, 2)),
            ]
        # Case 17
        elif self.mc and self.index < 17 * self.length:
            self.index += 1
            (
                alpha,
                beta,
                gamma,
                epsilon,
                (nu1, nu2, nu3),
                (lambda1, lambda2, lambda3),
                (rho1, rho2, rho3),
                parameter,
            ) = self.container[self.index - 1 - 16 * self.length]
            return [
                gamma,
                epsilon,
                alpha,
                beta,
                (rho1 - lambda1, rho2 - lambda2, rho3 - lambda3),
                (-lambda1, -lambda2, -lambda3),
                (nu1 - lambda1, nu2 - lambda2, nu3 - lambda3),
                np.transpose(parameter, (2, 3, 0, 1)),
            ]
        # Case 18
        elif self.mc and self.index < 18 * self.length:
            self.index += 1
            (
                alpha,
                beta,
                gamma,
                epsilon,
                (nu1, nu2, nu3),
                (lambda1, lambda2, lambda3),
                (rho1, rho2, rho3),
                parameter,
            ) = self.container[self.index - 1 - 17 * self.length]
            return [
                gamma,
                epsilon,
                beta,
                alpha,
                (rho1 - lambda1, rho2 - lambda2, rho3 - lambda3),
                (nu1 - lambda1, nu2 - lambda2, nu3 - lambda3),
                (-lambda1, -lambda2, -lambda3),
                np.transpose(parameter, (3, 2, 0, 1)),
            ]
        # Case 19
        elif self.mc and self.index < 19 * self.length:
            self.index += 1
            (
                alpha,
                beta,
                gamma,
                epsilon,
                (nu1, nu2, nu3),
                (lambda1, lambda2, lambda3),
                (rho1, rho2, rho3),
                parameter,
            ) = self.container[self.index - 1 - 18 * self.length]
            return [
                epsilon,
                alpha,
                beta,
                gamma,
                (-rho1, -rho2, -rho3),
                (nu1 - rho1, nu2 - rho2, nu3 - rho3),
                (lambda1 - rho1, lambda2 - rho2, lambda3 - rho3),
                np.transpose(parameter, (1, 2, 3, 0)),
            ]
        # Case 20
        elif self.mc and self.index < 20 * self.length:
            self.index += 1
            (
                alpha,
                beta,
                gamma,
                epsilon,
                (nu1, nu2, nu3),
                (lambda1, lambda2, lambda3),
                (rho1, rho2, rho3),
                parameter,
            ) = self.container[self.index - 1 - 19 * self.length]
            return [
                epsilon,
                alpha,
                gamma,
                beta,
                (-rho1, -rho2, -rho3),
                (lambda1 - rho1, lambda2 - rho2, lambda3 - rho3),
                (nu1 - rho1, nu2 - rho2, nu3 - rho3),
                np.transpose(parameter, (1, 3, 2, 0)),
            ]
        # Case 21
        elif self.mc and self.index < 21 * self.length:
            self.index += 1
            (
                alpha,
                beta,
                gamma,
                epsilon,
                (nu1, nu2, nu3),
                (lambda1, lambda2, lambda3),
                (rho1, rho2, rho3),
                parameter,
            ) = self.container[self.index - 1 - 20 * self.length]
            return [
                epsilon,
                beta,
                alpha,
                gamma,
                (nu1 - rho1, nu2 - rho2, nu3 - rho3),
                (-rho1, -rho2, -rho3),
                (lambda1 - rho1, lambda2 - rho2, lambda3 - rho3),
                np.transpose(parameter, (2, 1, 3, 0)),
            ]
        # Case 22
        elif self.mc and self.index < 22 * self.length:
            self.index += 1
            (
                alpha,
                beta,
                gamma,
                epsilon,
                (nu1, nu2, nu3),
                (lambda1, lambda2, lambda3),
                (rho1, rho2, rho3),
                parameter,
            ) = self.container[self.index - 1 - 21 * self.length]
            return [
                epsilon,
                beta,
                gamma,
                alpha,
                (nu1 - rho1, nu2 - rho2, nu3 - rho3),
                (lambda1 - rho1, lambda2 - rho2, lambda3 - rho3),
                (-rho1, -rho2, -rho3),
                np.transpose(parameter, (3, 1, 2, 0)),
            ]
        # Case 23
        elif self.mc and self.index < 23 * self.length:
            self.index += 1
            (
                alpha,
                beta,
                gamma,
                epsilon,
                (nu1, nu2, nu3),
                (lambda1, lambda2, lambda3),
                (rho1, rho2, rho3),
                parameter,
            ) = self.container[self.index - 1 - 22 * self.length]
            return [
                epsilon,
                beta,
                gamma,
                alpha,
                (lambda1 - rho1, lambda2 - rho2, lambda3 - rho3),
                (-rho1, -rho2, -rho3),
                (nu1 - rho1, nu2 - rho2, nu3 - rho3),
                np.transpose(parameter, (2, 3, 1, 0)),
            ]
        # Case 24
        elif self.mc and self.index < 24 * self.length:
            self.index += 1
            (
                alpha,
                beta,
                gamma,
                epsilon,
                (nu1, nu2, nu3),
                (lambda1, lambda2, lambda3),
                (rho1, rho2, rho3),
                parameter,
            ) = self.container[self.index - 1 - 23 * self.length]
            return [
                epsilon,
                beta,
                gamma,
                alpha,
                (lambda1 - rho1, lambda2 - rho2, lambda3 - rho3),
                (nu1 - rho1, nu2 - rho2, nu3 - rho3),
                (-rho1, -rho2, -rho3),
                np.transpose(parameter, (3, 2, 1, 0)),
            ]

        raise StopIteration

    def __len__(self):
        return self.length * (1 + 23 * int(self.mc))

    def __iter__(self):
        return self


@property
def _p44(spinham):
    r"""
    Parameters of (four spins & four sites) term of the Hamiltonian.

    .. math::

        \boldsymbol{J}_{4,4}(\boldsymbol{r}_{\nu,\alpha\beta}, \boldsymbol{r}_{\lambda,\alpha\gamma}, \boldsymbol{r}_{\rho,\alpha\varepsilon})

    of the term

    .. math::

        C_{4,4}
        \sum_{\substack{\mu, \nu, \alpha, \beta,\\ i, j, u, v}}
        J^{ijuv}_{4,4}(\boldsymbol{r}_{\nu,\alpha\beta}, \boldsymbol{r}_{\lambda,\alpha\gamma}, \boldsymbol{r}_{\rho,\alpha\varepsilon})
        S_{\mu,\alpha}^i
        S_{\mu+\nu,\beta}^j
        S_{\mu+\lambda, \gamma}^u
        S_{\mu+\rho, \varepsilon}^v

    Returns
    -------
    parameters : iterator
        List of parameters. The list has a form of

        .. code-block:: python

            [[alpha, beta, gamma, epsilon, nu, lambda, rho, J], ...]

        where

        ``alpha`` is an index of the atom located in the (0,0,0) unit cell.

        ``beta`` is an index of the atom located in the nu unit cell.

        ``gamma`` is an index of the atom located in the lambda unit cell.

        ``epsilon`` is an index of the atom located in the rho unit cell.

        ``nu`` defines the unit cell of the second atom (beta). It is a tuple of 3
        integers.

        ``lambda`` defines the unit cell of the third atom (gamma). It is a tuple of 3
        integers.

        ``rho`` defines the unit cell of the fourth atom (gamma). It is a tuple of 3
        integers.

        ``J`` is a (3, 3, 3, 3) :numpy:`ndarray`.

    See Also
    --------
    add_44
    remove_44
    """

    return _P44_iterator(spinham)


def _add_44(
    spinham,
    alpha: int,
    beta: int,
    gamma: int,
    epsilon: int,
    nu: tuple,
    _lambda: tuple,
    rho: tuple,
    parameter,
    replace=False,
) -> None:
    r"""
    Adds a (four spins & four sites) parameter to the Hamiltonian.

    Doubles of the bonds are managed automatically (independently of the convention of the
    Hamiltonian).


    Raises
    ------
    ValueError
        If an atom already has a parameter associated with it.

    Parameters
    ----------
    alpha : int
        Index of an atom from the (0, 0, 0) unit cell.

        ``0 <= alpha < len(spinham.atoms.names)``.
    beta : int
        Index of an atom from the nu unit cell.

        ``0 <= beta < len(spinham.atoms.names)``.
    gamma : int
        Index of an atom from the _lambda unit cell.

        ``0 <= gamma < len(spinham.atoms.names)``.
    epsilon : int
        Index of an atom from the rho unit cell.

        ``0 <= epsilon < len(spinham.atoms.names)``.
    nu : tuple of 3 int
        Three relative coordinates with respect to the three lattice vectors, that
        specify the unit cell for the second atom.

        .. math::

            \nu
            =
            (x_{\boldsymbol{a}_1}, x_{\boldsymbol{a}_2}, x_{\boldsymbol{a}_3})
    _lambda : tuple of 3 int
        Three relative coordinates with respect to the three lattice vectors, that
        specify the unit cell for the third atom.

        .. math::

            \lambda
            =
            (x_{\boldsymbol{a}_1}, x_{\boldsymbol{a}_2}, x_{\boldsymbol{a}_3})
    rho : tuple of 3 int
        Three relative coordinates with respect to the three lattice vectors, that
        specify the unit cell for the fourth atom.

        .. math::

            \rho
            =
            (x_{\boldsymbol{a}_1}, x_{\boldsymbol{a}_2}, x_{\boldsymbol{a}_3})

    parameter : (3, 3, 3, 3) |array-like|_
        Value of the parameter (:math:`3\times3\times3\times3` matrix).
    replace : bool, default False
        Whether to replace the value of the parameter if the triplet of atoms
        ``alpha, beta, gamma, nu, lambda`` or one of its duplicates already have a
        parameter associated with it.

    See Also
    --------
    p44
    remove_44

    Notes
    -----
    If ``spinham.convention.multiple_counting`` is ``True``, then this function adds
    the bond and all its duplicates to the Hamiltonian. It will cause an ``ValueError``
    to add the duplicate of the bond after the bond is added.

    If ``spinham.convention.multiple_counting`` is ``False``, then only the primary
    version of the bond is added to the Hamiltonian.

    For the definition of the primary version see
    :ref:`user-guide_theory-behind_multiple-counting`.
    """

    _validate_atom_index(index=alpha, atoms=spinham.atoms)
    _validate_atom_index(index=beta, atoms=spinham.atoms)
    _validate_atom_index(index=gamma, atoms=spinham.atoms)
    _validate_atom_index(index=epsilon, atoms=spinham.atoms)
    _validate_unit_cell_index(ijk=nu)
    _validate_unit_cell_index(ijk=_lambda)
    _validate_unit_cell_index(ijk=rho)
    spinham._reset_internals()

    parameter = np.array(parameter)

    alpha, beta, gamma, epsilon, nu, _lambda, rho, parameter = _get_primary_p44(
        alpha=alpha,
        beta=beta,
        gamma=gamma,
        epsilon=epsilon,
        nu=nu,
        _lambda=_lambda,
        rho=rho,
        parameter=parameter,
    )

    # TD-BINARY_SEARCH

    # Try to find the place for the new one inside the list
    index = 0
    while index < len(spinham._44):
        # If already present in the model
        if spinham._44[index][:7] == [alpha, beta, gamma, epsilon, nu, _lambda, rho]:
            # Either replace
            if replace:
                spinham._44[index] = [
                    alpha,
                    beta,
                    gamma,
                    epsilon,
                    nu,
                    _lambda,
                    rho,
                    parameter,
                ]
                return
            # Or raise an error
            raise ValueError(
                f"Parameter is already set for the quartet of atoms "
                f"{alpha}, {beta} {nu}, {gamma} {_lambda}, {epsilon} {rho}. Or for their duplicate."
            )

        # If it should be inserted before current element
        if spinham._44[index][:7] > [alpha, beta, gamma, epsilon, nu, _lambda, rho]:
            spinham._44.insert(
                index, [alpha, beta, gamma, epsilon, nu, _lambda, rho, parameter]
            )
            return

        index += 1

    # If it should be inserted at the end or at the beginning of the list
    spinham._44.append([alpha, beta, gamma, epsilon, nu, _lambda, rho, parameter])


def _remove_44(
    spinham,
    alpha: int,
    beta: int,
    gamma: int,
    epsilon: int,
    nu: tuple,
    _lambda: tuple,
    rho: tuple,
) -> None:
    r"""
    Removes a (four spins & four sites) parameter from the Hamiltonian.

    Duplicates of the bonds are managed automatically (independently of the convention of
    the Hamiltonian).

    Parameters
    ----------
    alpha : int
        Index of an atom from the (0, 0, 0) unit cell.

        ``0 <= alpha < len(spinham.atoms.names)``.
    beta : int
        Index of an atom from the nu unit cell.

        ``0 <= beta < len(spinham.atoms.names)``.
    gamma : int
        Index of an atom from the _lambda unit cell.

        ``0 <= gamma < len(spinham.atoms.names)``.
    epsilon : int
        Index of an atom from the rho unit cell.

        ``0 <= epsilon < len(spinham.atoms.names)``.
    nu : tuple of 3 int
        Three relative coordinates with respect to the three lattice vectors, that
        specify the unit cell for the second atom.

        .. math::

            \nu
            =
            (x_{\boldsymbol{a}_1}, x_{\boldsymbol{a}_2}, x_{\boldsymbol{a}_3})
    _lambda : tuple of 3 int
        Three relative coordinates with respect to the three lattice vectors, that
        specify the unit cell for the third atom.

        .. math::

            \lambda
            =
            (x_{\boldsymbol{a}_1}, x_{\boldsymbol{a}_2}, x_{\boldsymbol{a}_3})
    rho : tuple of 3 int
        Three relative coordinates with respect to the three lattice vectors, that
        specify the unit cell for the fourth atom.

        .. math::

            \rho
            =
            (x_{\boldsymbol{a}_1}, x_{\boldsymbol{a}_2}, x_{\boldsymbol{a}_3})

    See Also
    --------
    p44
    add_44

    Notes
    -----
    If ``spinham.convention.multiple_counting`` is ``True``, then this function removes
    all versions of the bond from the Hamiltonian.

    If ``spinham.convention.multiple_counting`` is ``False``, then this function removes
    the primary version of the given bond.

    For the definition of the primary version see
    :ref:`user-guide_theory-behind_multiple-counting`.
    """

    _validate_atom_index(index=alpha, atoms=spinham.atoms)
    _validate_atom_index(index=beta, atoms=spinham.atoms)
    _validate_atom_index(index=gamma, atoms=spinham.atoms)
    _validate_atom_index(index=epsilon, atoms=spinham.atoms)
    _validate_unit_cell_index(ijk=nu)
    _validate_unit_cell_index(ijk=_lambda)
    _validate_unit_cell_index(ijk=rho)

    alpha, beta, gamma, epsilon, nu, _lambda, rho = _get_primary_p44(
        alpha=alpha,
        beta=beta,
        gamma=gamma,
        epsilon=epsilon,
        nu=nu,
        _lambda=_lambda,
        rho=rho,
    )

    # TD-BINARY_SEARCH

    for index in range(len(spinham._44)):
        # As the list is sorted, there is no point in resuming the search
        # when a larger element is found
        if spinham._44[index][:7] > [alpha, beta, gamma, epsilon, nu, _lambda, rho]:
            return

        if spinham._44[index][:7] == [alpha, beta, gamma, epsilon, nu, _lambda, rho]:
            del spinham._44[index]
            spinham._reset_internals()
            return
