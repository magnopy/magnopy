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


R"""
Convention of spin Hamiltonian
"""

from magnopy._spinham._hamiltonian import SpinHamiltonian

# Save local scope at this moment
old_dir = set(dir())
old_dir.add("old_dir")


def make_supercell(spinham: SpinHamiltonian, supercell):
    r"""
    Creates a spin Hamiltonian on the supercell of the given one.

    Parameters
    ----------
    spinham : :py:class:`.SpinHamiltonian`
        Original spin Hamiltonian. ``spinham.cell`` is interpreted as the original
        unit cell.
    supercell : (3, ) tuple or list of int
        Repetitions of the unit cell (``spinham.cell``) along each lattice
        vector that define the unit cell. If :math:`(i, j, k)` is given, then the
        supercell is formally defined as
        :math:`(i\cdot\boldsymbol{a}_1, j\cdot\boldsymbol{a}_2, k\cdot\boldsymbol{a}_3)`,
        where :math:`(\boldsymbol{a}_1, \boldsymbol{a}_2, \boldsymbol{a}_3)` is the
        original cell (``spinham.cell``).

    Returns
    -------
    new_spinham : :py:class:`.SpinHamiltonian`
        Spin Hamiltonian that is defined on a ``supercell`` and has the same parameters
        as the given ``spinham`` propagated over the whole supercell.

    Examples
    --------

    First a simple example with only on-site parameters present in the Hamiltonian (i.e
    ``p21``) and two atoms per unit cell.

    .. doctest::

        >>> import numpy as np
        >>> import magnopy
        >>> # First create the original spin Hamiltonian
        >>> cell = np.eye(3)
        >>> atoms = dict(
        ...     names=["Fe1", "Fe2"],
        ...     positions=[[0, 0, 0], [0.5, 0.5, 0.5]],
        ...     spins=[1, 1],
        ...     g_factors=[2, 2],
        ... )
        >>> convention = magnopy.Convention(
        ...     spin_normalized=False, multiple_counting=True, c21=1
        ... )
        >>> spinham = magnopy.SpinHamiltonian(
        ...     cell=cell, atoms=atoms, convention=convention
        ... )
        >>> # Add an on-site parameter for both atoms
        >>> spinham.add_21(alpha=0, parameter=np.eye(3))
        >>> spinham.add_21(alpha=1, parameter=np.eye(3))
        >>> # Now create a spin Hamiltonian on the (2, 2, 2) supercell
        >>> new_spinham = magnopy.make_supercell(spinham=spinham, supercell=(2, 2, 2))
        >>> # Note that the unit cell of the new Hamiltonian is the supercell of the original one
        >>> new_spinham.cell
        array([[2., 0., 0.],
               [0., 2., 0.],
               [0., 0., 2.]])
        >>> # All the atoms of the supercell are now contained in the unit cell of the
        >>> # new Hamiltonian
        >>> len(spinham.atoms.names)
        2
        >>> len(new_spinham.atoms.names)
        16
        >>> for i in range(2):
        ...     print(
        ...         spinham.atoms.names[i],
        ...         spinham.atoms.positions[i],
        ...         spinham.atoms.spins[i],
        ...     )
        Fe1 [0, 0, 0] 1
        Fe2 [0.5, 0.5, 0.5] 1
        >>> for i in range(16):
        ...     print(
        ...         new_spinham.atoms.names[i],
        ...         new_spinham.atoms.positions[i],
        ...         new_spinham.atoms.spins[i],
        ...     )
        Fe1_0_0_0 [0.0, 0.0, 0.0] 1
        Fe2_0_0_0 [0.25, 0.25, 0.25] 1
        Fe1_1_0_0 [0.5, 0.0, 0.0] 1
        Fe2_1_0_0 [0.75, 0.25, 0.25] 1
        Fe1_0_1_0 [0.0, 0.5, 0.0] 1
        Fe2_0_1_0 [0.25, 0.75, 0.25] 1
        Fe1_1_1_0 [0.5, 0.5, 0.0] 1
        Fe2_1_1_0 [0.75, 0.75, 0.25] 1
        Fe1_0_0_1 [0.0, 0.0, 0.5] 1
        Fe2_0_0_1 [0.25, 0.25, 0.75] 1
        Fe1_1_0_1 [0.5, 0.0, 0.5] 1
        Fe2_1_0_1 [0.75, 0.25, 0.75] 1
        Fe1_0_1_1 [0.0, 0.5, 0.5] 1
        Fe2_0_1_1 [0.25, 0.75, 0.75] 1
        Fe1_1_1_1 [0.5, 0.5, 0.5] 1
        Fe2_1_1_1 [0.75, 0.75, 0.75] 1
        >>> # The parameters are updated as well
        >>> len(spinham.p21)
        2
        >>> for alpha, _ in spinham.p21:
        ...     print(alpha)
        0
        1
        >>> len(new_spinham.p21)
        16
        >>> for alpha, _ in new_spinham.p21:
        ...     print(alpha)
        0
        1
        2
        3
        4
        5
        6
        7
        8
        9
        10
        11
        12
        13
        14
        15

    Second example with one atom per unit cell, but with the bilinear interaction that
    connects two different atoms (``p22``).

    .. doctest::

        >>> import numpy as np
        >>> import magnopy
        >>> # First create the original spin Hamiltonian
        >>> cell = np.eye(3)
        >>> atoms = dict(names=["Fe"], positions=[[0, 0, 0]], spins=[1], g_factors=[2])
        >>> convention = magnopy.Convention(
        ...     spin_normalized=False, multiple_counting=True, c22=1
        ... )
        >>> spinham = magnopy.SpinHamiltonian(
        ...     cell=cell, atoms=atoms, convention=convention
        ... )
        >>> # Add an on-site parameter for both atoms
        >>> spinham.add_22(alpha=0, beta=0, nu=(1, 0, 0), parameter=np.eye(3))
        >>> # Now create a spin Hamiltonian on the (2, 2, 2) supercell
        >>> new_spinham = magnopy.make_supercell(spinham=spinham, supercell=(2, 2, 2))
        >>> len(new_spinham.atoms.names)
        8
        >>> for i in range(8):
        ...     print(
        ...         new_spinham.atoms.names[i],
        ...         new_spinham.atoms.positions[i],
        ...         new_spinham.atoms.spins[i],
        ...     )
        Fe_0_0_0 [0.0, 0.0, 0.0] 1
        Fe_1_0_0 [0.5, 0.0, 0.0] 1
        Fe_0_1_0 [0.0, 0.5, 0.0] 1
        Fe_1_1_0 [0.5, 0.5, 0.0] 1
        Fe_0_0_1 [0.0, 0.0, 0.5] 1
        Fe_1_0_1 [0.5, 0.0, 0.5] 1
        Fe_0_1_1 [0.0, 0.5, 0.5] 1
        Fe_1_1_1 [0.5, 0.5, 0.5] 1
        >>> # The bonds were recalculated automatically
        >>> for alpha, beta, nu, _ in spinham.p22:
        ...     print(alpha, beta, nu)
        0 0 (1, 0, 0)
        0 0 (-1, 0, 0)
        >>> for alpha, beta, nu, _ in new_spinham.p22:
        ...     print(alpha, beta, nu)
        0 1 (0, 0, 0)
        1 0 (1, 0, 0)
        2 3 (0, 0, 0)
        3 2 (1, 0, 0)
        4 5 (0, 0, 0)
        5 4 (1, 0, 0)
        6 7 (0, 0, 0)
        7 6 (1, 0, 0)
        1 0 (0, 0, 0)
        0 1 (-1, 0, 0)
        3 2 (0, 0, 0)
        2 3 (-1, 0, 0)
        5 4 (0, 0, 0)
        4 5 (-1, 0, 0)
        7 6 (0, 0, 0)
        6 7 (-1, 0, 0)


    """

    if supercell[0] < 1 or supercell[1] < 1 or supercell[2] < 1:
        raise ValueError(
            f"Supercell repetitions should be larger or equal to 1, got {supercell}"
        )

    new_cell = [supercell[i] * spinham.cell[i] for i in range(3)]

    new_atoms = {}

    for key in spinham.atoms:
        new_atoms[key] = []

    for k in range(supercell[2]):
        for j in range(supercell[1]):
            for i in range(supercell[0]):
                for atom_index in range(len(spinham.atoms.names)):
                    for key in spinham.atoms:
                        if key == "positions":
                            position = spinham.atoms.positions[atom_index]
                            new_position = [
                                (position[0] + i) / supercell[0],
                                (position[1] + j) / supercell[1],
                                (position[2] + k) / supercell[2],
                            ]
                            new_atoms["positions"].append(new_position)
                        elif key == "names":
                            new_atoms["names"].append(
                                f"{spinham.atoms.names[atom_index]}_{i}_{j}_{k}"
                            )
                        else:
                            new_atoms[key].append(spinham.atoms[key][atom_index])

    new_spinham = SpinHamiltonian(
        cell=new_cell, atoms=new_atoms, convention=spinham.convention
    )

    def get_new_indices(alpha, nu, ijk):
        nu = [nu[index] + ijk[index] for index in range(3)]

        i, j, k = [nu[index] % supercell[index] for index in range(3)]

        nu = [nu[index] // supercell[index] for index in range(3)]

        alpha = alpha + (i + j * supercell[0] + k * supercell[1] * supercell[0]) * len(
            spinham.atoms.names
        )

        return alpha, tuple(nu)

    for k in range(supercell[2]):
        for j in range(supercell[1]):
            for i in range(supercell[0]):
                # One spin
                for alpha, parameter in spinham._1:
                    alpha, _ = get_new_indices(alpha=alpha, nu=(0, 0, 0), ijk=(i, j, k))

                    new_spinham.add_1(alpha=alpha, parameter=parameter)

                # Two spins
                for alpha, parameter in spinham._21:
                    alpha, _ = get_new_indices(alpha=alpha, nu=(0, 0, 0), ijk=(i, j, k))

                    new_spinham.add_21(alpha=alpha, parameter=parameter)

                for alpha, beta, nu, parameter in spinham._22:
                    alpha, _ = get_new_indices(alpha=alpha, nu=(0, 0, 0), ijk=(i, j, k))
                    beta, nu = get_new_indices(alpha=beta, nu=nu, ijk=(i, j, k))

                    new_spinham.add_22(
                        alpha=alpha, beta=beta, nu=nu, parameter=parameter, replace=True
                    )

                # Three spins
                for alpha, parameter in spinham._31:
                    alpha, _ = get_new_indices(alpha=alpha, nu=(0, 0, 0), ijk=(i, j, k))

                    new_spinham.add_31(alpha=alpha, parameter=parameter)

                for alpha, beta, nu, parameter in spinham._32:
                    alpha, _ = get_new_indices(alpha=alpha, nu=(0, 0, 0), ijk=(i, j, k))
                    beta, nu = get_new_indices(alpha=beta, nu=nu, ijk=(i, j, k))

                    new_spinham.add_32(
                        alpha=alpha, beta=beta, nu=nu, parameter=parameter, replace=True
                    )

                for alpha, beta, gamma, nu, _lambda, parameter in spinham._33:
                    alpha, _ = get_new_indices(alpha=alpha, nu=(0, 0, 0), ijk=(i, j, k))
                    beta, nu = get_new_indices(alpha=beta, nu=nu, ijk=(i, j, k))
                    gamma, _lambda = get_new_indices(
                        alpha=gamma, nu=_lambda, ijk=(i, j, k)
                    )

                    new_spinham.add_33(
                        alpha=alpha,
                        beta=beta,
                        gamma=gamma,
                        nu=nu,
                        _lambda=_lambda,
                        parameter=parameter,
                        replace=True,
                    )

                # Four spins
                for alpha, parameter in spinham._41:
                    alpha, _ = get_new_indices(alpha=alpha, nu=(0, 0, 0), ijk=(i, j, k))

                    new_spinham.add_41(alpha=alpha, parameter=parameter)

                for alpha, beta, nu, parameter in spinham._421:
                    alpha, _ = get_new_indices(alpha=alpha, nu=(0, 0, 0), ijk=(i, j, k))
                    beta, nu = get_new_indices(alpha=beta, nu=nu, ijk=(i, j, k))

                    new_spinham.add_421(
                        alpha=alpha, beta=beta, nu=nu, parameter=parameter, replace=True
                    )

                for alpha, beta, nu, parameter in spinham._422:
                    alpha, _ = get_new_indices(alpha=alpha, nu=(0, 0, 0), ijk=(i, j, k))
                    beta, nu = get_new_indices(alpha=beta, nu=nu, ijk=(i, j, k))

                    new_spinham.add_422(
                        alpha=alpha, beta=beta, nu=nu, parameter=parameter, replace=True
                    )

                for alpha, beta, gamma, nu, _lambda, parameter in spinham._43:
                    alpha, _ = get_new_indices(alpha=alpha, nu=(0, 0, 0), ijk=(i, j, k))
                    beta, nu = get_new_indices(alpha=beta, nu=nu, ijk=(i, j, k))
                    gamma, _lambda = get_new_indices(
                        alpha=gamma, nu=_lambda, ijk=(i, j, k)
                    )

                    new_spinham.add_43(
                        alpha=alpha,
                        beta=beta,
                        gamma=gamma,
                        nu=nu,
                        _lambda=_lambda,
                        parameter=parameter,
                        replace=True,
                    )

                for (
                    alpha,
                    beta,
                    gamma,
                    epsilon,
                    nu,
                    _lambda,
                    rho,
                    parameter,
                ) in spinham._44:
                    alpha, _ = get_new_indices(alpha=alpha, nu=(0, 0, 0), ijk=(i, j, k))
                    beta, nu = get_new_indices(alpha=beta, nu=nu, ijk=(i, j, k))
                    gamma, _lambda = get_new_indices(
                        alpha=gamma, nu=_lambda, ijk=(i, j, k)
                    )
                    epsilon, rho = get_new_indices(alpha=epsilon, nu=rho, ijk=(i, j, k))
                    new_spinham.add_44(
                        alpha=alpha,
                        beta=beta,
                        gamma=gamma,
                        epsilon=epsilon,
                        nu=nu,
                        _lambda=_lambda,
                        rho=rho,
                        parameter=parameter,
                        replace=True,
                    )

    return new_spinham


# Populate __all__ with objects defined in this file
__all__ = list(set(dir()) - old_dir)
# Remove all semi-private objects
__all__ = [i for i in __all__ if not i.startswith("_")]
del old_dir
