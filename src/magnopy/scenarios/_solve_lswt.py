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
import wulfric

from magnopy._energy import Energy
from magnopy._lswt import LSWT
from magnopy._package_info import logo
from magnopy._parallelization import multiprocess_over_k
from magnopy.io import plot_dispersion
from magnopy._plotly_engine import PlotlyEngine
from magnopy._constants._icons import ICON_OUT_FILE
from magnopy.io._warnings import _envelope_warning


try:
    import scipy  # noqa F401

    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False

try:
    import plotly  # noqa F401

    PLOTLY_AVAILABLE = True
except ImportError:
    PLOTLY_AVAILABLE = False

try:
    import matplotlib.pyplot as plt  # noqa F401

    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False


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
    Computes magnon Hamiltonian at the level of Linear Spin Wave theory.

    Progress of calculation is shown in the standard output (``print()``). A bunch of the
    output files is created and saved on the disk inside the ``output_folder``.

    Parameters
    ----------

    spinham : :py:class:`.SpinHamiltonian`
        Spin Hamiltonian object.

    spin_directions : (M, 3) |array-like|_, optional.
        Directions of the local quantization axis for each spin. Magnitude of the vector
        is ignored, only the direction is considered. If ``None``, then Magnopy attempts
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
        Vector of external magnetic field (magnetic flux density, B), given in Tesla.

    output_folder : str, default "magnopy-results"
        Name for the folder where to save the output files. The folder is created if it
        does not exist.

    number_processors : int, optional
        Number of processors to be used in computation. By default Magnopy uses all
        available processes. Use ``number_processors=1`` to run in serial mode.

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

    spglib_symprec : float, default 1e-5
        .. versionadded:: 0.2.0

        Tolerance parameter for the space group symmetry search by |spglib|_. Reduce it
        if the space group is not the one you expected.


    Notes
    -----

    When using this function of Magnopy in your Python scripts make sure to safeguard
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
    ##                                  Filenames                                 ##
    ################################################################################
    # Create the output directory if it does not exist
    os.makedirs(output_folder, exist_ok=True)

    def envelope_path(pathname):
        if hide_personal_data:
            return pathname
        else:
            return os.path.abspath(pathname)

    SPIN_VECTORS_TXT = envelope_path(os.path.join(output_folder, "SPIN_VECTORS.txt"))
    SPIN_DIRECTIONS_HTML = envelope_path(
        os.path.join(output_folder, "SPIN_DIRECTIONS.html")
    )
    HIGH_SYMMETRY_POINTS_TXT = envelope_path(
        os.path.join(output_folder, "HIGH-SYMMETRY_POINTS.txt")
    )
    K_POINTS_HTML = envelope_path(os.path.join(output_folder, "K-POINTS.html"))
    K_POINTS_TXT = envelope_path(os.path.join(output_folder, "K-POINTS.txt"))
    OMEGAS_TXT = envelope_path(os.path.join(output_folder, "OMEGAS.txt"))
    OMEGAS_PNG = envelope_path(os.path.join(output_folder, "OMEGAS.png"))
    OMEGAS_IMAG_TXT = envelope_path(os.path.join(output_folder, "OMEGAS-IMAG.txt"))
    OMEGAS_IMAG_PNG = envelope_path(os.path.join(output_folder, "OMEGAS-IMAG.png"))
    DELTAS_TXT = envelope_path(os.path.join(output_folder, "DELTAS.txt"))
    DELTAS_PNG = envelope_path(os.path.join(output_folder, "DELTAS.png"))
    E_0_TXT = envelope_path(os.path.join(output_folder, "E_0.txt"))
    E_2_TXT = envelope_path(os.path.join(output_folder, "E_2.txt"))
    E_CORR_TXT = envelope_path(os.path.join(output_folder, "E_corr.txt"))
    ONE_OPERATOR_TERMS_TXT = envelope_path(
        os.path.join(output_folder, "ONE_OPERATOR_TERMS.txt")
    )

    all_good = True

    ################################################################################
    ##                       Logo, comment and plotly check                       ##
    ################################################################################
    # Print logo and a comment
    print(logo(date_time=True, logo_width=80))
    if comment is not None:
        print(f"\n{' Comment ':=^80}\n")
        print(comment)

    if no_html:
        print("\nHTML output is disabled by user (no_html=True).")

    if not PLOTLY_AVAILABLE and not no_html:
        print(
            _envelope_warning(
                "Cannot produce files\n  - "
                + "\n  - ".join(
                    [os.path.basename(_) for _ in [K_POINTS_HTML, SPIN_DIRECTIONS_HTML]]
                )
                + "\nbecause plotly is not available.\n"
                "You can install plotly with 'pip install plotly'"
            )
        )

    if not SCIPY_AVAILABLE and not no_html:
        print(
            _envelope_warning(
                "Cannot produce files\n  - "
                + "\n  - ".join([os.path.basename(_) for _ in [K_POINTS_HTML]])
                + "\nbecause scipy is not available."
                "\nYou can install scipy with 'pip install scipy'"
            )
        )

    if not MATPLOTLIB_AVAILABLE:
        print(
            _envelope_warning(
                "Cannot produce files\n  - "
                + "\n  - ".join(
                    [
                        os.path.basename(_)
                        for _ in [OMEGAS_PNG, OMEGAS_IMAG_PNG, DELTAS_PNG]
                    ]
                )
                + "\nbecause matplotlib is not available.\n"
                "You can install matplotlib with 'pip install matplotlib'"
            )
        )

    ################################################################################
    ##                              External effects                              ##
    ################################################################################
    print(f"\n{' External effects ':=^80}\n")

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
    ################################################################################
    ##                       Optimization of spin directions                      ##
    ################################################################################

    # Energy class will be needed later on as well,
    # thus it is outside of the if-block
    energy = Energy(spinham=spinham)

    if spin_directions is None:
        print(f"\n{' Optimization of spin directions ':=^80}\n")

        print(
            "Spin directions for the ground state are not given, attempt to optimize.\n"
        )

        # Tolerance parameters
        print(f"Energy tolerance       : {1e-5:.5e} meV")
        print(f"Torque tolerance       : {1e-5:.5e}")
        print(R"Target energy function : E^{(0)}")

        print(
            "Supercell             : 1 x 1 x 1 (original unit cell of the Hamiltonian)"
        )
        print(
            "\nNote: we recommend obtaining ground state outside of the magnopy-lswt program\n"
            "and provide --spin-directions argument to it. See magnopy-optimize-sd, for\n"
            "dedicated spin optimization of magnopy."
        )

        spin_directions = energy.optimize(
            energy_tolerance=1e-5, torque_tolerance=1e-5, quiet=False
        )
    else:
        # Normalize spin directions to unity
        spin_directions = np.array(spin_directions, dtype=float)
        spin_directions = (
            spin_directions / np.linalg.norm(spin_directions, axis=1)[:, np.newaxis]
        )

    ################################################################################
    ##                                Ground state                                ##
    ################################################################################
    print(f"\n{' Ground state ':=^80}\n")

    # Even if they are present in spinham.atoms
    # this function will get the correct ones
    spglib_types = wulfric.get_spglib_types(atoms=spinham.magnetic_atoms)

    name_n = max([4] + [len(name) for name in spinham.magnetic_atoms["names"]])
    print("Order of the atoms and their spglib types are")
    print(f"{'Name':<{name_n}} spglib_type")
    for n_i, name in enumerate(spinham.magnetic_atoms["names"]):
        print(f"{name:{name_n}} {spglib_types[n_i]:^11}")

    # Save spin directions and spin values to the .txt file
    np.savetxt(
        SPIN_VECTORS_TXT,
        np.concatenate(
            (spin_directions, np.array(spinham.magnetic_atoms["spins"])[:, np.newaxis]),
            axis=1,
        ),
        fmt="%12.8f %12.8f %12.8f   %12.8f",
    )
    print(
        f"\nDirections of spin vectors of the ground state and spin values are saved in "
        f"file\n{ICON_OUT_FILE} {SPIN_VECTORS_TXT}"
    )
    # Save spin directions as a .html file
    if not no_html and PLOTLY_AVAILABLE:
        pe = PlotlyEngine()
        pe.plot_cell(
            spinham.cell,
            color="Black",
            legend_label="Unit cell",
        )
        pe.plot_spin_directions(
            positions=np.array(spinham.magnetic_atoms["positions"]) @ spinham.cell,
            spin_directions=spin_directions,
            colors="#A47864",
            legend_label="Spins of the unit cell",
        )
        pe.save(
            output_name=SPIN_DIRECTIONS_HTML,
            axes_visible=True,
            legend_position="top",
            kwargs_write_html=dict(include_plotlyjs=True, full_html=True),
        )
        print(
            f"\nImage of the spin directions is saved in file\n{ICON_OUT_FILE} {SPIN_DIRECTIONS_HTML}"
        )

    # Classical energy
    E_0 = energy.E_0(spin_directions=spin_directions)
    with open(E_0_TXT, "w", encoding="utf-8") as f:
        f.write(f"{E_0:.8f} meV\n")
    print(
        f"\nClassical energy of optimized state (E_0 = {E_0:.3f} meV) is saved in file\n{ICON_OUT_FILE} {E_0_TXT}"
    )

    # Correction to classical energy
    E_corr = energy.E_corr(spin_directions=spin_directions)
    with open(E_CORR_TXT, "w", encoding="utf-8") as f:
        f.write(f"{E_corr:.8f} meV\n")
    print(
        f"Correction to the classical energy of optimized state (E_corr = {E_corr:.3f} meV) is saved in file\n{ICON_OUT_FILE} {E_CORR_TXT}"
    )

    ################################################################################
    ##                            K-points and k-path                             ##
    ################################################################################
    print(f"\n{' K-points and k-path ':=^80}\n")

    ticks = None
    labels = None
    kpoints_relative = np.array([[]])
    kpoints_absolute = np.array([[]])
    x_data = np.array([])
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

        x_data = np.cumsum(
            np.concatenate(
                (
                    [0.0],
                    np.linalg.norm(
                        kpoints_absolute[1:] - kpoints_absolute[:-1], axis=1
                    ),
                )
            )
        )

        print("K-points are provided by user.")

    else:
        print("Deducing k-points based on the crystal symmetry.")
        print("See wulfric.org for more details on procedure and conventions.")
        spglib_data = wulfric.get_spglib_data(
            cell=spinham.cell, atoms=spinham.atoms, spglib_symprec=spglib_symprec
        )
        print(
            f"\nspglib_symprec  : {spglib_symprec:.5e}.",
            f"Space group     : {spglib_data.space_group_number}",
            f"Bravais lattice : {spglib_data.crystal_family + spglib_data.centring_type}",
            "Convention      : HPKOT",
            sep="\n",
        )
        kp = wulfric.Kpoints.from_crystal(
            cell=spinham.cell,
            atoms=spinham.atoms,
            convention="HPKOT",
            spglib_data=spglib_data,
        )

        # Try to set custom k path
        if k_path is not None:
            try:
                kp.path = k_path
            except ValueError:
                all_good = False
                print(
                    _envelope_warning(
                        "User-provided k-path contains undefined labels of high-symmetry points."
                        "\nPre-defined points are\n  - "
                        + "\n  - ".join(kp.hs_names)
                        + "\nUsing recommended k-path instead.",
                    )
                )

        print(f"K-path          : {kp.path_string}")

        # Save pre-defined high-symmetry points in a .txt file
        label_n = max([5] + [len(label) for label in kp.hs_names])
        with open(HIGH_SYMMETRY_POINTS_TXT, "w", encoding="utf-8") as f:
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
            f"\nFull list of pre-defined high-symmetry points is saved in file\n{ICON_OUT_FILE} {HIGH_SYMMETRY_POINTS_TXT}"
        )

        kpoints_relative = kp.points(relative=True)
        kpoints_absolute = kpoints_relative @ kp.rcell
        x_data = kp.flat_points(relative=False)
        ticks = kp.ticks(relative=False)
        labels = kp.labels

        # Produce .html file with the hs points, k-path and brillouin zones
        if not no_html and PLOTLY_AVAILABLE and SCIPY_AVAILABLE:
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

            pe.save(output_name=K_POINTS_HTML)
            print(
                f"\nHigh-symmetry points and chosen k-path are plotted in\n{ICON_OUT_FILE} {K_POINTS_HTML}"
            )

    # Save k-points info to the .txt file
    np.savetxt(
        K_POINTS_TXT,
        np.concatenate(
            (kpoints_absolute, kpoints_relative, x_data[:, np.newaxis]), axis=1
        ),
        fmt="%12.8f %12.8f %12.8f   %12.8f %12.8f %12.8f   %12.8f",
        header=f"{'k_x':>12} {'k_y':>12} {'k_z':>12}   {'r_b1':>12} {'r_b2':>12} {'r_b3':>12}   {'flat index':>12}",
        comments="",
    )
    print(
        f"\nExplicit list of k-points is saved in file\n{ICON_OUT_FILE} {K_POINTS_TXT}"
    )

    ################################################################################
    ##                                    LSWT                                    ##
    ################################################################################
    print(f"\n{' LSWT ':=^80}\n")
    lswt = LSWT(spinham=spinham, spin_directions=spin_directions)

    # Correction to the classical energy
    E_2 = lswt.E_2()
    with open(E_2_TXT, "w", encoding="utf-8") as f:
        f.write(f"{E_2:.8f} meV\n")
    print(
        f"\nCorrection to the classical ground state energy (E_2 = {E_2:.3f} meV) is saved in file\n{ICON_OUT_FILE} {E_2_TXT}"
    )

    # One-operator coefficients
    one_operator_coefficients = lswt.O()
    np.savetxt(
        ONE_OPERATOR_TERMS_TXT,
        np.concatenate(
            (
                one_operator_coefficients.real[:, np.newaxis],
                one_operator_coefficients.imag[:, np.newaxis],
            ),
            axis=1,
        ),
        fmt="%12.8f %12.8f",
        header=f"{'Re(O_alpha)':>12} {'Im(O_alpha)':>12}",
        comments="",
    )
    print(
        f"\nCoefficients before one-operator term are saved in file\n{ICON_OUT_FILE} {ONE_OPERATOR_TERMS_TXT}"
    )
    print("(shall be zero if the ground state is correct)")

    if not np.allclose(
        one_operator_coefficients, np.zeros(one_operator_coefficients.shape)
    ):
        all_good = False
        print(
            _envelope_warning(
                "Coefficients before the one-operator terms are not zero. It might "
                "indicate that the ground state (spin directions) is not a ground state "
                "of the considered spin Hamiltonian and the results might not be "
                "meaningful.  If O_alpha << 1 then the problem can also be numerical "
                "(due to the finite point arithmetic) and the results are just fine in "
                "that case. Contact developers if you are in doubts: magnopy.org."
            )
        )

    # Compute data for each k-point
    print("\nStart calculations over k-points ... ", end="")
    results = multiprocess_over_k(
        kpoints=kpoints_absolute,
        function=lswt.diagonalize,
        relative=False,
        number_processors=number_processors,
    )
    omegas = np.array([i[0] for i in results]).T
    deltas = np.array([i[1] for i in results])
    n_modes = len(omegas)
    has_imaginary = not np.allclose(omegas.imag, np.zeros(omegas.imag.shape))
    has_nans = np.any(np.isnan(omegas))
    print("Done")

    if has_nans:
        all_good = False
        print(
            _envelope_warning(
                "Some eigenfrequiencies could not be computed (NaN values). It might "
                "indicate that the ground state (spin directions) is not a ground state "
                "of the considered spin Hamiltonian and the results might not be "
                "meaningful. The problem can also be numerical (due to the finite point "
                "arithmetic) and the results are just fine in that case. Contact "
                "developers if you are in doubts: magnopy.org."
            )
        )

    if has_imaginary:
        all_good = False
        print(
            _envelope_warning(
                "Eigenfrequiencies has non-zero imaginary component for some k vectors. "
                "It might indicate that the ground state (spin directions) is not a "
                "ground state of the considered spin Hamiltonian and the results might "
                "not be meaningful. If Im(omega) << 1 then the problem can also be "
                "numerical (due to the finite point arithmetic) and the results are just "
                "fine in that case. Contact developers if you are in doubts: magnopy.org."
            )
        )

    if np.any(omegas.real < -1e-8):
        all_good = False
        print(
            _envelope_warning(
                "Some eigenfrequiencies are negative. This might indicate that the "
                "ground state is not a ground state of the considered spin Hamiltonian "
                "and the results might not be meaningful. Minimum of the spectrum can "
                "indicate a better ground state. Stable ground state should always "
                "have non-negative eigenfrequiencies of the excitations."
            )
        )

    ################################################################################
    ##                                 Text output                                ##
    ################################################################################
    print(f"\n{' Output ':=^80}\n")

    # Omegas
    np.savetxt(
        OMEGAS_TXT,
        omegas.real.T,
        fmt=("%15.6e " * n_modes)[:-1],
        header=" ".join([f"{f'mode {i + 1}':>15}" for i in range(n_modes)]),
        comments="",
    )
    print(f"\nOmegas are saved in file\n{ICON_OUT_FILE} {OMEGAS_TXT}")

    # Deltas
    np.savetxt(DELTAS_TXT, deltas.real, fmt="%10.6e", header="Delta", comments="")
    print(f"Deltas are saved in file\n{ICON_OUT_FILE} {DELTAS_TXT}")

    # Imaginary omegas
    if has_imaginary or has_nans:
        np.savetxt(
            OMEGAS_IMAG_TXT,
            omegas.imag.T,
            fmt=("%15.6e " * n_modes)[:-1],
            header=" ".join([f"{f'mode {i + 1}':>15}" for i in range(n_modes)]),
            comments="",
        )
        print(
            f"Imaginary part of omegas is saved in file\n{ICON_OUT_FILE} {OMEGAS_IMAG_TXT}"
        )

    ################################################################################
    ##                                 png output                                 ##
    ################################################################################

    if MATPLOTLIB_AVAILABLE:
        # Omegas
        plot_dispersion(
            modes=omegas.real,
            x_data=x_data,
            ticks=ticks,
            labels=labels,
            output_filename=OMEGAS_PNG,
            ylabel=R"$\omega_{\alpha}(\boldsymbol{k})$",
        )
        print(f"Plot is saved in file\n{ICON_OUT_FILE} {OMEGAS_PNG}")

        # Deltas
        plot_dispersion(
            modes=deltas.real,
            x_data=x_data,
            ticks=ticks,
            labels=labels,
            output_filename=DELTAS_PNG,
            ylabel=R"$\Delta(\boldsymbol{k})$",
        )
        print(f"Plot is saved in file\n{ICON_OUT_FILE} {DELTAS_PNG}")

        # Imaginary omegas
        if has_imaginary or has_nans:
            plot_dispersion(
                modes=omegas.imag,
                x_data=x_data,
                ticks=ticks,
                labels=labels,
                output_filename=OMEGAS_IMAG_PNG,
                ylabel=R"$\mathcal{Im}(\omega_{\alpha}(\boldsymbol{k}))$",
            )
            print(
                f"Plot of imaginary part is saved in file\n{ICON_OUT_FILE} {OMEGAS_IMAG_PNG}"
            )

    if all_good:
        print(f"\n{' Finished OK ':=^80}")
    else:
        print(f"\n{' Finished with WARNINGS ':=^80}")


# Populate __all__ with objects defined in this file
__all__ = list(set(dir()) - old_dir)
# Remove all semi-private objects
__all__ = [i for i in __all__ if not i.startswith("_")]
del old_dir
