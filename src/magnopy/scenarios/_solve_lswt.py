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
import wulfric

from magnopy._energy import Energy
from magnopy._lswt import LSWT
from magnopy._package_info import logo
from magnopy._parallelization import multiprocess_over_k
from magnopy.io._k_resolved import plot_k_resolved
from magnopy._plotly_engine import PlotlyEngine

try:
    import scipy  # noqa F401

    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False


# Save local scope at this moment
old_dir = set(dir())
old_dir.add("old_dir")


def solve_lswt(
    spinham,
    spin_directions=None,
    k_path=None,
    kpoints=None,
    relative=False,
    magnetic_field=None,
    output_folder="magnopy-results",
    number_processors=None,
    comment=None,
    no_html=False,
    hide_personal_data=False,
    spglib_symprec=1e-5,
) -> None:
    r"""
    Solves the spin Hamiltonian at the level of Linear Spin Wave theory.
    Outputs progress in the standard output (``print()``) and saves some data to
    the files on the disk.

    Parameters
    ----------
    spinham : :py:class:`.SpinHamiltonian`
        Spin Hamiltonian.
    spin_directions : (M, 3) |array-like|_, optional.
        Directions of the local quantization axis for each spin. Magnitude of the vector
        is ignored, only the direction is considered. If ``None``, then magnopy attempts
        to optimize classical energy of spin Hamiltonian to determine spin directions.
    k_path : str, optional
        Specification of the k-path. The format is "G-X-Y|G-Z" For more details
        on the format see documentation of |wulfric|_. If nothing given, then the
        k-path is computed by |wulfric|_ automatically based on the lattice type.
        Ignored if ``kpoints`` are given.
    kpoints : (N, 3) |array-like|_, optional
        Explicit list of k-points to be used instead of automatically generated.
    relative : bool, default False
        If ``relative == True``, then ``kpoints`` are interpreted as given relative to
        the reciprocal unit cell. Otherwise it is interpreted as given in absolute
        coordinates.
    magnetic_field : (3, ) |array-like|_
        Vector of external magnetic field, given in Tesla.
    output_folder : str, default "magnopy-results"
        Name for the folder where to save the output files. If the folder does not exist
        then it will be created.
    number_processors : int, optional
        Number of processors to be used in computation. By default magnopy uses all
        available processes. Use ``number_processors=1`` to run in serial mode.
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
    spglib_symprec : float, default 1e-5
        Tolerance parameter for the space group symmetry search by |spglib|_. Reduce it
        if the space group is not the one you expected.

        .. versionadded:: 0.2.0

    Notes
    -----

    When using this function of magnopy in your Python scripts make sure to safeguard
    your script with the

    .. code-block:: python

        import magnopy

        # Import more stuff
        # or
        # Define your functions, classes

        if __name__ == "__main__":

            # Write your executable code here

    For more information refer to the  "Safe importing of main module" section in
    |multiprocessing|_ docs.

    """

    ################################################################################
    ##                   Data verification and envelope function                  ##
    ################################################################################
    def envelope_path(pathname):
        if hide_personal_data:
            return pathname
        else:
            return os.path.abspath(pathname)

    # Create the output directory if it does not exist
    os.makedirs(output_folder, exist_ok=True)

    all_good = True

    ################################################################################
    ##                              Logo and comment                              ##
    ################################################################################
    # Print logo and a comment
    print(logo(date_time=True))
    if comment is not None:
        print(f"\n{' Comment ':=^90}\n")
        print(comment)

    ################################################################################
    ##                                Ground state                                ##
    ################################################################################
    # Print header
    print(f"\n{' Ground state ':=^90}\n")

    # Add magnetic field if any
    if magnetic_field is not None:
        spinham.add_magnetic_field(h=magnetic_field)

    # Get energy class
    energy = Energy(spinham=spinham)

    # Optimize spin directions
    if spin_directions is None:
        print("Spin directions are not given, start to optimize ...")

        spin_directions = energy.optimize(
            energy_tolerance=1e-5, torque_tolerance=1e-5, quiet=False
        )
        print("Optimization is done.")
    # Or normalize them
    else:
        print("Spin directions of the ground state are provided by the user.")

        spin_directions = np.array(spin_directions, dtype=float)
        spin_directions = (
            spin_directions / np.linalg.norm(spin_directions, axis=1)[:, np.newaxis]
        )

    # Save spin directions and spin values to the .txt file
    filename = os.path.join(output_folder, "SPIN_VECTORS.txt")
    np.savetxt(
        filename,
        np.concatenate(
            (spin_directions, np.array(spinham.magnetic_atoms["spins"])[:, np.newaxis]),
            axis=1,
        ),
        fmt="%12.8f %12.8f %12.8f   %12.8f",
    )
    print(
        f"\nDirections of spin vectors of the ground state and spin values are saved in file\n  {envelope_path(filename)}"
    )
    # Save spin directions as a .html file
    if not no_html:
        filename = os.path.join(output_folder, "SPIN_DIRECTIONS.html")

        pe = PlotlyEngine()

        pe.plot_cell(
            spinham.cell,
            color="Black",
            legend_label="Unit cell",
        )

        pe.plot_spin_directions(
            positions=np.array(spinham.magnetic_atoms.positions) @ spinham.cell,
            spin_directions=spin_directions,
            colors="#A47864",
            legend_label="Spins of the unit cell",
        )

        pe.save(
            output_name=filename,
            axes_visible=True,
            legend_position="top",
            kwargs_write_html=dict(include_plotlyjs=True, full_html=True),
        )

        print(
            f"\nImage of spin directions is saved in file\n  {envelope_path(filename)}.html\n"
        )

    # Output order of atoms and their spglib types as understood by wulfric
    spglib_types = wulfric.get_spglib_types(atoms=spinham.magnetic_atoms)
    name_n = max([4] + [len(name) for name in spinham.magnetic_atoms["names"]])
    print("Order of the atoms is")
    print(f"{'Name':<{name_n}} spglib_type")
    for n_i, name in enumerate(spinham.magnetic_atoms["names"]):
        print(f"{name:{name_n}} {spglib_types[n_i]:>11}")

    # Output classical energy
    E_0 = energy.E_0(spin_directions=spin_directions)
    print(f"\n{'Classic ground state energy (E_0)':<51} is {E_0:>15.6f} meV\n")

    ################################################################################
    ##                            K-points and k-path                             ##
    ################################################################################
    # Treat kpoints
    print(f"\n{' K-points and k-path ':=^90}\n")

    if kpoints is not None:
        if relative:
            kpoints_relative = np.array(kpoints, dtype=float)
            kpoints_absolute = kpoints_relative @ wulfric.cell.get_reciprocal(
                cell=spinham.cell
            )
        else:
            kpoints_absolute = np.array(kpoints, dtype=float)
            kpoints_relative = kpoints_absolute @ np.linalg.inv(
                wulfric.cell.get_reciprocal(cell=spinham.cell)
            )

        flat_indices = np.concatenate(
            (
                [0.0],
                np.linalg.norm(kpoints_absolute[1:] - kpoints_absolute[:-1], axis=1),
            )
        )

        print("K-points are provided by the user.")

    else:
        spglib_data = wulfric.get_spglib_data(
            cell=spinham.cell, atoms=spinham.atoms, spglib_symprec=spglib_symprec
        )
        print(
            "Deducing k-points based on the crystal symmetry.",
            f"spglib_symprec is  {spglib_symprec:.5e}.",
            f"Space group is     {spglib_data.space_group_number}",
            f"Bravais lattice is {spglib_data.crystal_family + spglib_data.centring_type}",
            "Using convention of HPKOT paper. See docs of wulfric for details (wulfric.org).",
            sep="\n",
        )
        kp = wulfric.Kpoints.from_crystal(
            cell=spinham.cell,
            atoms=spinham.atoms,
            convention="HPKOT",
            spglib_data=spglib_data,
        )

        # Save pre-defined high-symmetry points in a .txt file
        filename = os.path.join(output_folder, "HIGH-SYMMETRY_POINTS.txt")
        label_n = max([5] + [len(label) for label in kp.hs_names])
        with open(filename, "w") as f:
            f.write(
                f"{'Label':{label_n}} {'k_x':>12} {'k_y':>12} {'k_z':>12}    {'r_b1':>12} {'r_b2':>12} {'r_b3':>12}\n"
            )
            for label in kp.hs_names:
                k_rel = kp.hs_coordinates[label]
                k_abs = k_rel @ kp.rcell
                f.write(
                    f"{label:<{label_n}} {k_abs[0]:12.8f} {k_abs[1]:12.8f} {k_abs[2]:12.8f}    {k_rel[0]:12.8f} {k_rel[1]:12.8f} {k_rel[2]:12.8f}\n"
                )
        print(
            f"\nFull list of pre-defined high-symmetry points is saved in file\n  {envelope_path(filename)}\n"
        )

        # Try to set custom k path
        if k_path is not None:
            print("K-path is provided by the user.")
            try:
                kp.path = k_path
            except ValueError:
                all_good = False
                print(f"\n{'  WARNING  ':!^90}")
                print(
                    "User-provided k-path contains labels of high-symmetry point that are not defined.",
                    "See file",
                    f"  {envelope_path(filename)}",
                    "for the list of pre-defined high-symmetry points.",
                    "Using recommended k-path instead:",
                    f"  {kp.path_string}",
                    sep="\n",
                )
                print(f"{'  END OF WARNING  ':!^90}\n")
        else:
            print("Using recommended k-path:", f"  {kp.path_string}")
        kpoints_relative = kp.points(relative=True)
        kpoints_absolute = kpoints_relative @ kp.rcell
        flat_indices = kp.flat_points(relative=False)

        # Produce .html file with the hs points, k-path and brillouin zones
        if not no_html:
            if SCIPY_AVAILABLE:
                filename = os.path.join(output_folder, "K-POINTS.html")
                pe = wulfric.PlotlyEngine()

                prim_cell, _ = wulfric.crystal.get_primitive(
                    cell=spinham.cell,
                    atoms=spinham.atoms,
                    convention="SC",
                    spglib_data=spglib_data,
                )
                pe.plot_brillouin_zone(
                    cell=prim_cell,
                    color="red",
                    legend_label="Brillouin zone of the primitive cell",
                )
                pe.plot_brillouin_zone(
                    cell=spinham.cell,
                    color="chocolate",
                    legend_label="Brillouin zone of the spinham.cell",
                )
                pe.plot_kpath(kp=kp)
                pe.plot_kpoints(kp=kp, only_from_kpath=True)

                pe.save(output_name=filename)
                print(
                    f"\nHigh-symmetry points and chosen k-path are plotted in\n  {envelope_path(filename)}"
                )
            else:
                print(
                    "\nCan not plot Brillouin zone without scipy. Please install it with\n  pip install scipy"
                )

    # Save k-points info to the .txt file
    filename = os.path.join(output_folder, "K-POINTS.txt")
    np.savetxt(
        filename,
        np.concatenate(
            (
                kpoints_absolute,
                kpoints_relative,
                flat_indices[:, np.newaxis],
            ),
            axis=1,
        ),
        fmt="%12.8f %12.8f %12.8f   %12.8f %12.8f %12.8f   %12.8f",
        header=f"{'k_x':>12} {'k_y':>12} {'k_z':>12}   {'r_b1':>12} {'r_b2':>12} {'r_b3':>12}   {'flat index':>12}",
        comments="",
    )
    print(f"\nExplicit list of k-points is saved in file\n  {envelope_path(filename)}")

    ################################################################################
    ##                                    LSWT                                    ##
    ################################################################################
    print(f"\n{' Start LSWT ':=^90}\n")
    lswt = LSWT(spinham=spinham, spin_directions=spin_directions)

    # Output correction energy
    print(
        f"{'Correction to the classic ground state energy (E_2)':<50} is {lswt.E_2:>15.6f} meV\n"
    )

    # Output one-operator coefficients
    print(
        "Coefficients before one-operator terms (shall be zero if the ground state is correct)"
    )
    print("  " + "\n  ".join([f"{o:12.8f}" for o in lswt.O]))
    if not np.allclose(lswt.O, np.zeros(lswt.O.shape)):
        all_good = False
        print(f"\n{'  WARNING  ':!^90}")
        print(
            "Coefficients before the one-operator terms are not zero. It might indicate  that",
            "the ground state (spin directions) is not a ground state of the considered spin",
            "Hamiltonian. The results might not be meaningful. If coefficients are << 1, that might",
            "be an artifact of the finite point arithmetic and the results might be just fine.",
            sep="\n",
        )
        print(f"{'  END OF WARNING  ':!^90}\n")

    # Compute data for each k-point
    print("\nStart calculations over k-points ... ", end="")
    results = multiprocess_over_k(
        kpoints=kpoints_absolute,
        function=lswt.diagonalize,
        relative=False,
        number_processors=number_processors,
    )
    omegas = np.array([i[0] for i in results])
    deltas = np.array([i[1] for i in results])
    n_modes = len(omegas[0])
    print("Done")

    # Save omegas to the .txt file
    filename = os.path.join(output_folder, "OMEGAS.txt")
    np.savetxt(
        filename,
        omegas.real,
        fmt=("%15.6e " * n_modes)[:-1],
        header=" ".join([f"{f'mode {i + 1}':>15}" for i in range(n_modes)]),
        comments="",
    )
    print(f"\nOmegas are saved in file\n  {envelope_path(filename)}")

    # Plot omegas as a .png
    filename = filename[:-4] + ".png"
    # TODO: REFACTOR
    plot_k_resolved(
        data=omegas.real,
        kp=kp,
        output_filename=filename,
        ylabel=R"$\omega_{\alpha}(\boldsymbol{k})$, meV",
    )
    print(f"Plot is saved in file\n  {envelope_path(filename)}")

    # Check for the imaginary part
    if not np.allclose(omegas.imag, np.zeros(omegas.imag.shape)):
        all_good = False
        print(f"\n{'  WARNING  ':!^90}")
        print(
            "Eigenfrequiencies has non-zero imaginary component for some k vectors. It might\n"
            "indicate that the ground state (spin directions) is not a ground state of the\n"
            "considered spin Hamiltonian. The results might not be meaningful.\n"
        )
        filename = os.path.join(output_folder, "OMEGAS-IMAG.txt")
        np.savetxt(
            filename,
            omegas.imag,
            fmt=("%15.6e " * n_modes)[:-1],
            header=" ".join([f"{f'mode {i + 1}':>15}" for i in range(n_modes)]),
            comments="",
        )
        print(f"Imaginary part of omegas is saved in file\n  {envelope_path(filename)}")

        filename = filename[:-4] + ".png"
        # TODO: REFACTOR
        plot_k_resolved(
            data=omegas.imag,
            kp=kp,
            output_filename=filename,
            ylabel=R"$\mathcal{Im}(\omega_{\alpha}(\boldsymbol{k}))$, meV",
        )
        print(f"Plot of imaginary part is saved in file\n  {envelope_path(filename)}")
        print(f"{'  END OF WARNING  ':!^90}\n")

    # Save deltas to the .txt file
    filename = os.path.join(output_folder, "DELTAS.txt")
    np.savetxt(filename, deltas.real, fmt="%10.6e", header="Delta", comments="")
    print(f"Deltas are saved in file\n  {envelope_path(filename)}")

    # Plot deltas as a .png
    filename = filename[:-4] + ".png"
    # TODO: REFACTOR
    plot_k_resolved(
        data=deltas.real,
        kp=kp,
        output_filename=filename,
        ylabel=R"$\Delta(\boldsymbol{k})$, meV",
    )
    print(f"Plot is saved in file\n  {envelope_path(filename)}")

    if all_good:
        print(f"\n{' Finished OK ':=^90}")
    else:
        print(f"\n{' Finished with WARNINGS ':=^90}")


# Populate __all__ with objects defined in this file
__all__ = list(set(dir()) - old_dir)
# Remove all semi-private objects
__all__ = [i for i in __all__ if not i.startswith("_")]
del old_dir
