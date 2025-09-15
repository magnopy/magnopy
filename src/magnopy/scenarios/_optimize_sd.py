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

from magnopy._energy import Energy
from magnopy._package_info import logo
from magnopy._spinham._supercell import make_supercell
from magnopy._plotly_engine import PlotlyEngine


# Save local scope at this moment
old_dir = set(dir())
old_dir.add("old_dir")


def optimize_sd(
    spinham,
    supercell=(1, 1, 1),
    magnetic_field=None,
    energy_tolerance=1e-5,
    torque_tolerance=1e-5,
    output_folder="magnopy-results",
    comment=None,
    no_html=False,
    hide_personal_data=False,
) -> None:
    r"""
    Optimizes classical energy of spin Hamiltonian and finds a set of spin directions
    that describe local minima on the energy landscape.

    Parameters
    ----------
    spinham : :py:class:`.SpinHamiltonian`
        Spin Hamiltonian.
    supercell : (3, ) tuple of int
        If different from ``(1, 1, 1)``, then a supercell Hamiltonian is constructed and
        spins are varied within the supercell and not within a unit cell.

        .. versionadded:: 0.2.0
    magnetic_field : (3, ) |array-like|_
        Vector of external magnetic field, given in Tesla.
    energy_tolerance : float, default 1e-5
        Tolerance parameter. Difference between classical energies of two consecutive
        optimization steps.
    torque_tolerance : float, default 1e-5
        Tolerance parameter. Maximum torque among all spins.
    output_folder : str, default "magnopy-results"
        Name for the folder where to save the output files. If the folder does not exist
        then it will be created.
    comment : str, optional
        Any comment to output right after the logo.
    no_html : bool, default False
        Whether to produce .html files with interactive representation of the data.
        If ``no_html=False``, then requires |plotly|_ to be installed.

        .. versionadded:: 0.2.0
    hide_personal_data : bool, default False
        Whether to use ``os.path.abspath()`` when printing the paths to the output and
        input files.

        .. versionadded:: 0.2.0

    Raises
    ------
    ValueError
        If ``len(supercell) != 3``.
    ValueError
        If ``supercell[0] < 1`` or ``supercell[1] < 1`` or ``supercell[2] < 1``.
    """

    def envelope_path(pathname):
        if hide_personal_data:
            return pathname
        else:
            return os.path.abspath(pathname)

    # Create the output directory if it does not exist
    os.makedirs(output_folder, exist_ok=True)

    # Validate supercell
    supercell = tuple(supercell)
    if len(supercell) != 3:
        raise ValueError(
            f"Expected a tuple of three int, got {len(supercell)} elements."
        )
    if supercell[0] < 1 or supercell[1] < 1 or supercell[2] < 1:
        raise ValueError(f"Supercell repetitions must be >=1, got {supercell}.")

    # Print logo and a comment
    print(logo(date_time=True))
    if comment is not None:
        print(f"\n{' Comment ':=^90}\n")
        print(comment)

    # Output optimization parameters
    print(f"\n{' Start optimization ':=^90}\n")
    print(f"Energy tolerance : {energy_tolerance:.5e}")
    print(f"Torque tolerance : {torque_tolerance:.5e}")

    # Add magnetic field if any
    if magnetic_field is not None:
        spinham.add_magnetic_field(h=magnetic_field)

    # Make a supercell if needed
    original_spinham = spinham
    if supercell != (1, 1, 1):
        spinham = make_supercell(spinham=spinham, supercell=supercell)
        print(
            f"Minimizing on the supercell of {supercell[0]} x {supercell[1]} x {supercell[2]} unit cells."
        )
    else:
        print("Minimizing on the original unit cell of the Hamiltonian.")

    # Make an initial guess
    initial_guess = np.random.uniform(low=-1, high=1, size=(spinham.M, 3))
    initial_guess = initial_guess / np.linalg.norm(initial_guess, axis=1)[:, np.newaxis]
    # Save an initial guess to the .txt file
    filename = os.path.join(output_folder, "INITIAL_GUESS.TXT")
    np.savetxt(
        filename,
        initial_guess,
        fmt="%12.8f %12.8f %12.8f",
    )
    print(
        f"\nSpin directions of the initial guess are saved in file\n  {envelope_path(filename)}"
    )

    # Optimize spin directions
    energy = Energy(spinham=spinham)
    spin_directions = energy.optimize(
        initial_guess=initial_guess,
        energy_tolerance=energy_tolerance,
        torque_tolerance=torque_tolerance,
        quiet=False,
    )
    print("Optimization is done.")

    # Output classical energy
    E_0 = energy.E_0(spin_directions=spin_directions)
    print(f"\nClassic ground state energy (E_0) is {E_0:>15.6f} meV")

    # Save spin directions to a .txt file
    filename = os.path.join(output_folder, "SPIN_DIRECTIONS.txt")
    np.savetxt(
        filename,
        spin_directions,
        fmt="%12.8f %12.8f %12.8f",
    )
    print(f"\nOptimized spin directions are saved in file\n  {envelope_path(filename)}")

    # Compute spin's positions
    positions = np.array(spinham.magnetic_atoms.positions) @ spinham.cell
    # Save spin positions to a .txt file
    filename = os.path.join(output_folder, "SPIN_POSITIONS.txt")
    np.savetxt(
        filename,
        positions,
        fmt="%12.8f %12.8f %12.8f",
    )
    print(f"\nSpin positions are saved in file\n  {envelope_path(filename)}")

    if not no_html:
        filename = os.path.join(output_folder, "SPIN_DIRECTIONS.html")

        pe = PlotlyEngine()

        pe.plot_cell(
            original_spinham.cell,
            color="Black",
            legend_label="Unit cell",
        )

        original_uc_spins_indices = [i for i in range(original_spinham.M)]

        pe.plot_spin_directions(
            positions=positions[original_uc_spins_indices],
            spin_directions=spin_directions[original_uc_spins_indices],
            colors="#A47864",
            legend_label="Spins of the unit cell",
        )

        if supercell != (1, 1, 1):
            other_uc_spins_indices = [i for i in range(original_spinham.M, spinham.M)]
            pe.plot_spin_directions(
                positions=positions[other_uc_spins_indices],
                spin_directions=spin_directions[other_uc_spins_indices],
                colors="#535FCF",
                legend_label="Spins of other unit cells",
            )

        pe.plot_spin_directions(
            positions=positions,
            spin_directions=initial_guess,
            colors="#0DB00D",
            legend_label="Initial guess",
        )

        pe.save(
            output_name=filename,
            axes_visible=True,
            legend_position="top",
            kwargs_write_html=dict(include_plotlyjs=True, full_html=True),
        )

        print(
            f"\nImage of spin directions is saved in file\n  {envelope_path(filename)}"
        )

    print(f"\n{' Finished ':=^90}")


# Populate __all__ with objects defined in this file
__all__ = list(set(dir()) - old_dir)
# Remove all semi-private objects
__all__ = [i for i in __all__ if not i.startswith("_")]
del old_dir
