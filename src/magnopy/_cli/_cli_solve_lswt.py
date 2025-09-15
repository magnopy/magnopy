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
import warnings
import sys
from argparse import ArgumentParser, RawDescriptionHelpFormatter

import numpy as np

from magnopy._package_info import logo
from magnopy.io._grogu import load_grogu
from magnopy.io._spin_directions import read_spin_directions
from magnopy.io._tb2j import load_tb2j
from magnopy.scenarios._solve_lswt import solve_lswt


def manager():
    parser = get_parser()

    args = parser.parse_args()

    # Handle deprecations
    if args.make_sd_image is not None:
        warnings.warn(
            "This argument was deprecated in the release v0.2.0. The spin direction image is now plotted by default, please use --no-html if you want to disable it. This argument will be removed from magnopy in March of 2026"
        )

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    # Load spin directions
    if args.spin_directions is None:
        pass
    elif len(args.spin_directions) == 1:
        args.spin_directions = read_spin_directions(filename=args.spin_directions[0])
    else:
        args.spin_directions = np.array(args.spin_directions)
        args.spin_directions = args.spin_directions.reshape(
            (len(args.spin_directions) // 3, 3)
        )

    # Process spin values
    if args.spin_values is not None:
        args.spin_values = [float(tmp) for tmp in args.spin_values]

    # Load kpoints
    kpoints = []
    if args.kpoints is not None:
        with open(args.kpoints, "r") as f:
            for i, line in enumerate(f):
                # Remove comment lines
                if line.startswith("#"):
                    continue
                # Remove inline comments and leading/trailing whitespaces
                line = line.split("#")[0].strip()
                # Check for empty lines empty lines
                if line:
                    line = line.split()
                    if len(line) != 3:
                        raise ValueError(
                            f"Expected three numbers per line (in line{i}),"
                            f"got: {len(line)}."
                        )

                    kpoints.append(list(map(float, line)))

        args.kpoints = kpoints

    # Load spin Hamiltonian
    if args.spinham_source.lower() == "tb2j":
        spinham = load_tb2j(
            filename=args.spinham_filename, spin_values=args.spin_values
        )
    elif args.spinham_source.lower() == "grogu":
        spinham = load_grogu(filename=args.spinham_filename)
    else:
        raise ValueError(
            'Supported sources of spin Hamiltonian are "GROGU" and "TB2J", '
            f'got "{args.spinham_source}".'
        )
    if args.hide_personal_data:
        spinham_filename = args.spinham_filename
    else:
        spinham_filename = os.path.abspath(args.spinham_filename)
    comment = (
        f'Source of the parameters is "{args.spinham_source}".\n'
        f"Loaded parameters of the spin Hamiltonian from the file\n  "
        f"{spinham_filename}."
    )

    solve_lswt(
        spinham=spinham,
        spin_directions=args.spin_directions,
        k_path=args.k_path,
        kpoints=args.kpoints,
        relative=args.relative,
        magnetic_field=args.magnetic_field,
        output_folder=args.output_folder,
        number_processors=args.number_processors,
        comment=comment,
        no_html=args.no_html,
        hide_personal_data=args.hide_personal_data,
        spglib_symprec=args.spglib_symprec,
    )


def get_parser():
    parser = ArgumentParser(
        description=logo()
        + "\n\nThis script solves the spin Hamiltonian at the level of "
        "Linear Spin Wave Theory (LSWT) and outputs (almost) every possible quantity.",
        formatter_class=RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "-sf",
        "--spinham-filename",
        type=str,
        metavar="filename",
        default=None,
        required=True,
        help="Path to the spin Hamiltonian file, from where the parameters would be read.",
    )
    parser.add_argument(
        "-ss",
        "--spinham-source",
        type=str,
        metavar="name",
        default=None,
        required=True,
        choices=["GROGU", "TB2J"],
        help='Source of the spin Hamiltonian. Either "GROGU" or "TB2J"',
    )
    parser.add_argument(
        "-sd",
        "--spin-directions",
        nargs="*",
        type=str,
        default=None,
        metavar="S1_x S2_y S3_z ...",
        help="To fully define the system for the calculations of magnons one need the "
        "information about the ground state in addition to the parameters of the "
        "Hamiltonian. There are two ways to give this information to magnopy:\n"
        " * Give a path to the file. In the file there should be M lines with three "
        "numbers in each. The order of the lines would match the order of magnetic "
        "atoms in the spin Hamiltonian."
        " * Give a sequence of 3*M numbers directly to this parameter.\n"
        "If none provided, then magnopy attempts to optimize the spin directions prior "
        "to the LSWT calculations.",
    )
    parser.add_argument(
        "-sv",
        "--spin-values",
        nargs="*",
        type=str,
        metavar="S1 S2 S3 ...",
        help="In the case when the parameters of spin Hamiltonian comes from TB2J, one "
        "might want to change the values of spins to be closer to half-integers. This "
        "option allows that. Order of the M numbers should match the order of magnetic "
        "atoms in the spin Hamiltonian. Note that those numbers are always positive. To "
        "specify AFM order use opposite spin directions and not spin values of the "
        "opposite sign.",
    )
    parser.add_argument(
        "-kp",
        "--k-path",
        default=None,
        metavar="GAMMA-X-S|GAMMA-Y",
        type=str,
        help="Path of high symmetry k-points for the plots of dispersion and other "
        "quantities.",
    )
    parser.add_argument(
        "-kps",
        "--kpoints",
        type=str,
        default=None,
        help="Alternatively one could provide an explicit list of k-points for calculation. "
        "In that case provide a path to the file, in which each k-point is given in a "
        "separate line with three numbers per line.",
    )
    parser.add_argument(
        "-r",
        "--relative",
        default=False,
        action="store_true",
        help="When an explicit list of k-points is given, this option specify whether "
        "to consider them as relative or absolute coordinates. Absolute by default.",
    )
    parser.add_argument(
        "-mf",
        "--magnetic-field",
        default=None,
        nargs=3,
        type=float,
        help="Vector of external magnetic field, given in the units of Tesla.",
    )
    parser.add_argument(
        "-of",
        "--output-folder",
        type=str,
        default="magnopy-results",
        help="Folder where all output files of magnopy wil be saved.",
    )
    parser.add_argument(
        "-np",
        "--number-processors",
        type=int,
        default=None,
        help="Number of processes for multithreading. Uses all available processors by "
        "default. Pass 1 to run in serial.",
    )
    parser.add_argument(
        "-no-html",
        "--no-html",
        action="store_true",
        default=False,
        help="html files are generally heavy (~> 5 Mb). This option allows to disable "
        "their production to save disk space.",
    )
    parser.add_argument(
        "-hpd",
        "--hide-personal-data",
        action="store_true",
        default=False,
        help="Whether to strip the parts of the paths as to hide the file structure of "
        "you personal computer.",
    )
    parser.add_argument(
        "-spg-s",
        "--spglib-symprec",
        type=float,
        default=1e-5,
        help="Tolerance parameter for the space group symmetry search by spglib.",
    )

    # Deprecated arguments
    parser.add_argument(
        "-msdi",
        "--make-sd-image",
        nargs=3,
        type=int,
        default=None,
        help="make_sd_image is deprecated, use --no-html instead. This arguments will be removed from magnopy in March of 2026",
    )

    return parser
