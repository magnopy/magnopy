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


import os

import numpy as np

from magnopy._energy import Energy
from magnopy._package_info import logo
from magnopy._spinham._supercell import make_supercell
from magnopy._plotly_engine import PlotlyEngine
from magnopy._constants._icons import ICON_OUT_FILE


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
    quantum_correction=False,
) -> None:
    r"""
    Optimizes classical energy of spin Hamiltonian and finds a set of spin directions
    that describes local minimum of the energy landscape.

    Progress of calculation is shown in the standard output (``print()``). A bunch of the
    output files is created and saved on the disk inside the ``output_folder``.

    Parameters
    ----------

    spinham : :py:class:`.SpinHamiltonian`
        Spin Hamiltonian object.

    supercell : (3, ) tuple of int, default (1, 1, 1)
        .. versionadded:: 0.2.0

        Tuple with number of repetitions of the unit cell along each lattice vector.
        Spin directions are varied within the supercell during the optimization. By
        default, optimization is done within the unit cell.

    magnetic_field : (3, ) |array-like|_
        Vector of external magnetic field (magnetic flux density, B), given in Tesla.

    energy_tolerance : float, default 1e-5
        Tolerance parameter. Difference between classical energies of two consecutive
        optimization steps.

    torque_tolerance : float, default 1e-5
        Tolerance parameter. Maximum torque among all spins at the current optimization
        step.

    output_folder : str, default "magnopy-results"
        Name for the folder where to save the output files. The folder is created if it
        does not exist.

    comment : str, optional
        Any comment, that will be shown in the standard output right after the Magnopy's
        logo.

    no_html : bool, default False
        .. versionadded:: 0.2.0

        Whether to produce .html files. If ``no_html = False``, then |plotly|_ is expected
        to be available.

    hide_personal_data : bool, default False
        .. versionadded:: 0.2.0

        If ``False``, then ``os.path.abspath(pathname)`` is used to show full paths to the output
        and input files. If ``True``, then only ``pathname`` is used.

    quantum_correction : bool, default False
        Whether to include quantum correction to the energy in the optimization. If
        ``True``, then it optimizes :py:func:`.Energy.E_0` + :py:func:`.Energy.E_corr`.
        If ``False``, then it optimizes :py:func:`.Energy.E_0` only.

        .. warning::
            This option is experimental. Will be improved and tested in future
            releases. Use with caution.

    Raises
    ------

    ValueError
        If ``len(supercell) != 3``.

    ValueError
        If ``supercell[0] < 1`` or ``supercell[1] < 1`` or ``supercell[2] < 1``.
    """
    ################################################################################
    ##                                  Filenames                                 ##
    ################################################################################
    # Create the output directory if it does not exist
    os.makedirs(output_folder, exist_ok=True)

    def envelope_path(pathname):
        if hide_personal_data:
            return pathname
        else:
            return os.path.abspath(pathname)

    INITIAL_GUESS_TXT = envelope_path(os.path.join(output_folder, "INITIAL_GUESS.txt"))
    SPIN_DIRECTIONS_TXT = envelope_path(
        os.path.join(output_folder, "SPIN_DIRECTIONS.txt")
    )
    SPIN_POSITIONS_TXT = envelope_path(
        os.path.join(output_folder, "SPIN_POSITIONS.txt")
    )
    SPIN_DIRECTIONS_HTML = envelope_path(
        os.path.join(output_folder, "SPIN_DIRECTIONS.html")
    )
    E_0_TXT = envelope_path(os.path.join(output_folder, "E_0.txt"))
    E_CORR_TXT = envelope_path(os.path.join(output_folder, "E_corr.txt"))

    ################################################################################
    ##                           Input data verification                          ##
    ################################################################################

    # Supercell
    supercell = tuple(supercell)
    if len(supercell) != 3:
        raise ValueError(
            f"Expected a tuple of three int for supercell, got {len(supercell)} elements."
        )
    if supercell[0] < 1 or supercell[1] < 1 or supercell[2] < 1:
        raise ValueError(f"Supercell repetitions must be >=1, got {supercell}.")

    ################################################################################
    ##                              Logo and comment                              ##
    ################################################################################
    print(logo(date_time=True, logo_width=80))
    if comment is not None:
        print(f"\n{' Comment ':=^80}\n")
        print(comment)

    ################################################################################
    ##                          Optimization parameters                           ##
    ################################################################################
    print(f"\n{' Optimization parameters ':=^80}\n")

    # Tolerance parameters
    print(f"Energy tolerance       : {energy_tolerance:.5e} meV")
    print(f"Torque tolerance       : {torque_tolerance:.5e}")

    # Target energy function
    print(
        f"Target energy function : {'E^{(0)}' if not quantum_correction else 'E^{(0)} + E^{corr}'}"
    )

    # Magnetic field
    if magnetic_field is not None:
        print(
            f"Magnetic flux density : "
            f"|{magnetic_field[0]:.5f}, {magnetic_field[1]:.5f}, {magnetic_field[2]:.5f}|"
            f" = {np.linalg.norm(magnetic_field):.5f} Tesla"
        )
        spinham.add_magnetic_field(B=magnetic_field)
    else:
        print("Magnetic flux density : None")

    # Supercell
    print(
        f"Supercell             : {supercell[0]} x {supercell[1]} x {supercell[2]}",
        end="",
    )
    original_spinham = spinham
    if supercell != (1, 1, 1):
        spinham = make_supercell(spinham=spinham, supercell=supercell)
        print(" (constructed supercell of the Hamiltonian)")
    else:
        print(" (original unit cell of the Hamiltonian)")

    # Initial guess
    initial_guess = np.random.uniform(low=-1, high=1, size=(spinham.M, 3))
    initial_guess = initial_guess / np.linalg.norm(initial_guess, axis=1)[:, np.newaxis]
    np.savetxt(INITIAL_GUESS_TXT, initial_guess, fmt="%12.8f %12.8f %12.8f")
    print(
        f"\nSpin directions of the initial guess are saved in file\n{ICON_OUT_FILE} {INITIAL_GUESS_TXT}"
    )

    ################################################################################
    ##                                Optimization                                ##
    ################################################################################
    print(f"\n{' Start optimization ':=^80}\n")

    energy = Energy(spinham=spinham)
    spin_directions = energy.optimize(
        initial_guess=initial_guess,
        energy_tolerance=energy_tolerance,
        torque_tolerance=torque_tolerance,
        quiet=False,
        quantum_correction=quantum_correction,
    )
    print("Optimization is done.")

    ################################################################################
    ##                                 Text output                                ##
    ################################################################################
    print(f"\n{' Output ':=^80}\n")

    # Classical energy
    E_0 = energy.E_0(spin_directions=spin_directions)
    with open(E_0_TXT, "w", encoding="utf-8") as f:
        f.write(f"{E_0:.8f} meV\n")
    print(
        f"Classical energy of optimized state (E_0 = {E_0:.3f} meV) is saved in file\n{ICON_OUT_FILE} {E_0_TXT}"
    )

    # Correction to classical energy
    E_corr = energy.E_corr(spin_directions=spin_directions)
    with open(E_CORR_TXT, "w", encoding="utf-8") as f:
        f.write(f"{E_corr:.8f} meV\n")
    print(
        f"Correction to the classical energy of optimized state (E_corr = {E_corr:.3f} meV) is saved in file\n{ICON_OUT_FILE} {E_CORR_TXT}"
    )

    # Optimized spin directions
    np.savetxt(SPIN_DIRECTIONS_TXT, spin_directions, fmt="%12.8f %12.8f %12.8f")
    print(
        f"\nOptimized spin directions are saved in file\n{ICON_OUT_FILE} {SPIN_DIRECTIONS_TXT}"
    )

    # Real-space absolute positions of magnetic centers
    positions = np.array(spinham.magnetic_atoms["positions"]) @ spinham.cell
    np.savetxt(SPIN_POSITIONS_TXT, positions, fmt="%12.8f %12.8f %12.8f")
    print(f"\nSpin positions are saved in file\n{ICON_OUT_FILE} {SPIN_POSITIONS_TXT}")

    ################################################################################
    ##                                 HTML output                                ##
    ################################################################################
    if not no_html:
        try:
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
                other_uc_spins_indices = [
                    i for i in range(original_spinham.M, spinham.M)
                ]
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
                output_name=SPIN_DIRECTIONS_HTML,
                axes_visible=True,
                legend_position="top",
                kwargs_write_html=dict(include_plotlyjs=True, full_html=True),
            )

            print(
                f"\nImage of spin directions is saved in file\n{ICON_OUT_FILE} {SPIN_DIRECTIONS_HTML}"
            )
        except ImportError:
            print(
                "\nCannot produce .html files because plotly is not available.\n"
                "You can install plotly with 'pip install plotly'"
            )
    else:
        print("\nHTML output is disabled by user (no_html=True).\n")

    print(f"\n{' Finished ':=^80}")


# Populate __all__ with objects defined in this file
__all__ = list(set(dir()) - old_dir)
# Remove all semi-private objects
__all__ = [i for i in __all__ if not i.startswith("_")]
del old_dir
