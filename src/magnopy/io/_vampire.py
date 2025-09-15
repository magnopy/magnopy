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


import os

import numpy as np
from wulfric.cell import get_params
from wulfric.crystal import get_atom_species

from magnopy._constants._internal_units import ENERGY
from magnopy._package_info import logo
from magnopy._parameters._p22 import to_dmi, to_symm_anisotropy
from magnopy._spinham._convention import Convention
from magnopy._spinham._hamiltonian import SpinHamiltonian

# Save local scope at this moment
old_dir = set(dir())
old_dir.add("old_dir")


def dump_vampire(
    spinham: SpinHamiltonian,
    seedname="vampire",
    anisotropic=True,
    dmi=True,
    custom_mask=None,
    decimals=5,
    materials=None,
    no_logo=False,
):
    """
    Save the Hamiltonian in the format suitable for |Vampire|_.

    Parameters
    ----------
    spinham : :py:class:`.SpinHamiltonian`
        Spin Hamiltonian to be saved.
    seedname : str, default "vampire"
        Seedname for the .UCF and .mat files. Extensions are added automatically.
        Input file always have the name "input".
    anisotropic : bool, default True
        Whether to output anisotropic exchange.
    dmi : bool, default True
        Whether to output DMI exchange.
    custom_mask : func, optional
        Custom mask for the exchange parameter. Function which take (3,3) numpy:`ndarray`
        as an input and returns (3,3) numpy:`ndarray` as an output. If given, then
        ``anisotropic`` and ``dmi`` parameters are ignored.
    decimals : int, default 4
        Number of decimals to be printed (only for the exchange values).
    materials : list of int, optional
        List of materials for the atoms. Has to have the same length as the number of
        magnetic atoms in the ``spinham``. Order is the same as in :py:attr:`.SpinHamiltonian.magnetic_atoms`.
        If not given, each atom will be considered as a separate material. Starting from 0.
    no_logo : bool, default False
        Whether to print the logo in the output files.

    Returns
    -------
    content : str
        Content of the .UCF file if ``filename`` is not given.
    """

    head, _ = os.path.split(seedname)

    if head != "":
        os.makedirs(head, exist_ok=True)

    dump_vampire_ucf(
        spinham,
        filename=seedname + ".UCF",
        anisotropic=anisotropic,
        dmi=dmi,
        custom_mask=custom_mask,
        decimals=decimals,
        materials=materials,
        no_logo=no_logo,
    )
    dump_vampire_mat(
        spinham,
        filename=seedname + ".mat",
        materials=materials,
        no_logo=no_logo,
    )
    with open(os.path.join(head, "input-template"), "w", encoding="utf-8") as file:
        if not no_logo:
            file.write(logo(comment=True, date_time=True) + "\n")

        file.write(
            "\n".join(
                [
                    "#------------------------------------------",
                    f"material:file={seedname}.mat",
                    f"material:unit-cell-file = {seedname}.UCF",
                    "#------------------------------------------",
                    "# TODO: your simulation parameters",
                ]
            )
        )


def dump_vampire_mat(
    spinham: SpinHamiltonian, filename=None, materials=None, no_logo=False
):
    """
    Write .mat file for |Vampire|_.

    Parameters
    ----------
    spinham : :py:class:`.SpinHamiltonian`
        Spin Hamiltonian to be saved.
    filename : str, optional
        Name for the .mat file. No extensions is added automatically.
        If not given, the output is returned as a string.
    materials : list of int, optional
        List of materials for the atoms. Has to have the same length as the number of
        magnetic atoms in the ``spinham``. If not given, each atom will be considered
        as a separate material. Material index starts from 0 and should contain all
        consecutive integers between 0 and number of materials. Number of materials
        cannot be higher than number of magnetic atoms.
    no_logo : bool, default False
        Whether to print the logo in the output files.

    Notes
    -----
    Examples of the correct ``materials`` list for 5 magnetic atoms

    .. code-block:: python

        [0, 0, 0, 0, 0]
        [1, 3, 2, 1, 0]
        [0, 1, 2, 3, 4]
    """

    if materials is None:
        materials = [i for i in range(len(spinham.magnetic_atoms.names))]
    else:
        if len(materials) != len(spinham.magnetic_atoms.names):
            raise ValueError(
                f"Expected {len(spinham.magnetic_atoms.names)} materials, got "
                f"{len(materials)}."
            )
        materials_pool = set(materials)
        higher_material = max(materials_pool)
        for i in range(0, higher_material + 1):
            if i not in materials_pool:
                raise ValueError(
                    f"Materials indices should be consecutive integers between 0 and "
                    f"{higher_material}. Missing {i}."
                )

    if no_logo:
        text = []
    else:
        text = [logo(comment=True, date_time=True)]

    text.append(f"material:num-materials = {max(materials) + 1}")

    for i in range(spinham.M):
        if materials[i] not in materials[:i]:
            m_i = materials[i] + 1
            text.append("#---------------------------------------------------")
            text.append(f"# Material {m_i}")
            text.append("#---------------------------------------------------")
            text.append(
                f"material[{m_i}]:material-name = {spinham.magnetic_atoms.names[i]}"
            )
            text.append(
                f"material[{m_i}]:material-element = {get_atom_species(spinham.magnetic_atoms.names[i])}"
            )
            text.append(
                f"material[{m_i}]:atomic-spin-moment={spinham.magnetic_atoms.spins[i] * spinham.magnetic_atoms.g_factors[i]} ! muB"
            )
            text.append(f"material[{m_i}]:initial-spin-direction = random")
            text.append(f"material[{m_i}]:damping-constant = 0.1")
            text.append(f"material[{m_i}]:uniaxial-anisotropy-constant = 0.0")

    text.append("#---------------------------------------------------")

    text = "\n".join(text)

    if filename is None:
        return "".join(text)

    with open(filename, "w", encoding="utf-8") as file:
        file.write("".join(text))


def dump_vampire_ucf(
    spinham: SpinHamiltonian,
    filename=None,
    anisotropic=True,
    dmi=True,
    custom_mask=None,
    decimals=5,
    materials=None,
    no_logo=False,
):
    """
    Write .UCF file for |Vampire|_.

    Parameters
    ----------
    spinham : :py:class:`.SpinHamiltonian`
        Spin Hamiltonian to be saved.
    filename : str, optional
        Name for the .UCF file. No extension is added automatically.
        If not given, the output is returned as a string.
    anisotropic : bool, default True
        Whether to output anisotropic exchange.
    dmi : bool, default True
        Whether to output DMI exchange.
    custom_mask : func, optional
        Custom mask for the exchange parameter. Function which take (3,3) numpy:`ndarray`
        as an input and returns (3,3) numpy:`ndarray` as an output. If given, then
        ``anisotropic`` and ``dmi`` parameters are ignored.
    decimals : int, default 4
        Number of decimals to be printed (only for the exchange values).
    materials : list of int, optional
        List of materials for the atoms. Has to have the same length as the number of
        magnetic atoms in the ``spinham``. Order is the same as in :py:attr:`.SpinHamiltonian.magnetic_atoms`.
        If not given, each atom will be considered as a separate material.
    no_logo : bool, default False
        Whether to print the logo in the output files.

    Returns
    -------
    content : str
        Content of the .UCF file if ``filename`` is not given.
    """
    if materials is None:
        materials = [i for i in range(len(spinham.magnetic_atoms.names))]

    original_convention = spinham.convention
    spinham.convention = Convention.get_predefined(name="Vampire")

    if no_logo:
        text = []
    else:
        text = [logo(comment=True, date_time=True)]

    a, b, c, _, _, _ = get_params(spinham.cell)
    text.append("# Unit cell size:")
    text.append(f"{a:.8f} {b:.8f} {c:.8f}")
    text.append("# Unit cell lattice vectors:")
    text.append(
        f"{spinham.cell[0][0]:15.8f} {spinham.cell[0][1]:15.8f} {spinham.cell[0][2]:15.8f}"
    )
    text.append(
        f"{spinham.cell[1][0]:15.8f} {spinham.cell[1][1]:15.8f} {spinham.cell[1][2]:15.8f}"
    )
    text.append(
        f"{spinham.cell[2][0]:15.8f} {spinham.cell[2][1]:15.8f} {spinham.cell[2][2]:15.8f}"
    )
    text.append("# Atoms")
    text.append(f"{len(spinham.magnetic_atoms.names)} {len(np.unique(materials))}")

    for alpha in range(spinham.M):
        position = spinham.magnetic_atoms.positions[alpha]
        text.append(
            f"{alpha:<5} {position[0]:15.8f} {position[1]:15.8f} {position[2]:15.8f} {materials[alpha]:>5}"
        )

    text.append("# Interactions")
    text.append(f"{len(spinham.p22)} tensorial")

    IID = 0
    fmt = f"{7 + decimals}.{decimals}e"

    # Write (two spins & one site)
    for alpha, J in spinham.p21:
        alpha = spinham.map_to_magnetic[alpha]
        if custom_mask is not None:
            J = custom_mask(J)
        else:
            if not dmi:
                J -= to_dmi(J, matrix_form=True)
            if not anisotropic:
                J -= to_symm_anisotropy(J)
        J = J * ENERGY
        text.append(
            f"{IID:<5} {alpha:>3} {alpha:>3}  {0:>2} {0:>2} {0:>2}  "
            f"{J[0][0]:{fmt}} {J[0][1]:{fmt}} {J[0][2]:{fmt}} "
            f"{J[1][0]:{fmt}} {J[1][1]:{fmt}} {J[1][2]:{fmt}} "
            f"{J[2][0]:{fmt}} {J[2][1]:{fmt}} {J[2][2]:{fmt}}"
        )
        IID += 1

    # Write (two spins & two sites)
    bonds = []
    for alpha, beta, nu, J in spinham.p22:
        alpha = spinham.map_to_magnetic[alpha]
        beta = spinham.map_to_magnetic[beta]
        if custom_mask is not None:
            J = custom_mask(J)
        else:
            if not dmi:
                J -= to_dmi(J, matrix_form=True)
            if not anisotropic:
                J -= to_symm_anisotropy(J)
        # print(alpha, beta, nu)
        # print(J, end="\n\n")
        J = J * ENERGY
        r_alpha = np.array(spinham.magnetic_atoms.positions[alpha])
        r_beta = np.array(spinham.magnetic_atoms.positions[beta])

        distance = np.linalg.norm((r_beta - r_alpha + nu) @ spinham.cell)
        bonds.append([alpha, beta, nu, J, distance])

    bonds = sorted(bonds, key=lambda x: x[4])
    for alpha, beta, (i, j, k), J, _ in bonds:
        text.append(
            f"{IID:<5} {alpha:>3} {beta:>3}  {i:>2} {j:>2} {k:>2}  "
            f"{J[0][0]:{fmt}} {J[0][1]:{fmt}} {J[0][2]:{fmt}} "
            f"{J[1][0]:{fmt}} {J[1][1]:{fmt}} {J[1][2]:{fmt}} "
            f"{J[2][0]:{fmt}} {J[2][1]:{fmt}} {J[2][2]:{fmt}}"
        )
        IID += 1

    spinham.convention = original_convention

    text = "\n".join(text)

    if filename is None:
        return text

    with open(filename, "w", encoding="utf-8") as file:
        file.write(text)


# Populate __all__ with objects defined in this file
__all__ = list(set(dir()) - old_dir)
# Remove all semi-private objects
__all__ = [i for i in __all__ if not i.startswith("_")]
del old_dir
