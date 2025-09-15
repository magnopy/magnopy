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
import sys
from argparse import ArgumentParser, RawDescriptionHelpFormatter
import warnings

from magnopy._package_info import logo
from magnopy.io._grogu import load_grogu
from magnopy.io._tb2j import load_tb2j
from magnopy.scenarios._optimize_sd import optimize_sd


def manager():
    parser = get_parser()

    args = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    # Handle deprecations
    if args.make_sd_image is not None:
        warnings.warn(
            "This argument was deprecated in the release v0.2.0. The spin direction image is now plotted by default, please use --no-html if you want to disable it. This argument will be removed from magnopy in March of 2026"
        )

    # Process spin values
    if args.spin_values is not None:
        args.spin_values = [float(tmp) for tmp in args.spin_values]

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

    optimize_sd(
        spinham=spinham,
        supercell=args.supercell,
        magnetic_field=args.magnetic_field,
        energy_tolerance=args.energy_tolerance,
        torque_tolerance=args.torque_tolerance,
        output_folder=args.output_folder,
        comment=comment,
        no_html=args.no_html,
        hide_personal_data=args.hide_personal_data,
    )


def get_parser():
    parser = ArgumentParser(
        description=logo()
        + "\n\nThis script optimizes classical energy of the spin Hamiltonian and "
        "finds the spin directions that describe a local minima of the energy landscape.",
        formatter_class=RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "-sf",
        "--spinham-filename",
        type=str,
        metavar="FILENAME",
        default=None,
        required=True,
        help="Path to the spin Hamiltonian file, from where the parameters would be read.",
    )
    parser.add_argument(
        "-ss",
        "--spinham-source",
        type=str,
        metavar="KEYWORD",
        default=None,
        required=True,
        choices=["GROGU", "TB2J"],
        help='Source of the spin Hamiltonian. Either "GROGU" or "TB2J"',
    )
    parser.add_argument(
        "-sv",
        "--spin-values",
        nargs="*",
        type=str,
        metavar=("S1", "S2"),
        help="In the case when the parameters of spin Hamiltonian comes from TB2J, one "
        "might want to change the values of spins to be closer to half-integers. This "
        "option allows that. Order of the M numbers should match the order of magnetic "
        "atoms in the spin Hamiltonian. Note that those numbers are always positive. To "
        "specify AFM order use opposite spin directions and not spin values of the "
        "opposite sign.",
    )
    parser.add_argument(
        "-s",
        "--supercell",
        nargs=3,
        type=int,
        default=(1, 1, 1),
        metavar=("xa_1", "xa_2", "xa_3"),
        help="Specification of the supercell for the spin optimization. Expects three "
        "integers as an input. Pass 1 1 1 to optimize within the original unit cell.",
    )
    parser.add_argument(
        "-et",
        "--energy-tolerance",
        default=1e-5,
        type=float,
        help="Tolerance parameter. Difference between classical energies of two "
        "consecutive optimization steps.",
    )
    parser.add_argument(
        "-tt",
        "--torque-tolerance",
        default=1e-5,
        type=float,
        help="Maximum torque among all spins.",
    )
    parser.add_argument(
        "-mf",
        "--magnetic-field",
        default=None,
        nargs=3,
        metavar=("h_x", "h_y", "h_z"),
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
        "-no-html",
        "--no-html",
        action="store_true",
        default=False,
        help="Disable plotting of the spin direction image in the .html format. html "
        "files are generally heavy (~> 5 Mb). This option allows to disable their "
        "production to save disk space.",
    )
    parser.add_argument(
        "-hpd",
        "--hide-personal-data",
        action="store_true",
        default=False,
        help="Whether to strip the parts of the paths as to hide the file structure of "
        "you personal computer.",
    )

    # Deprecated in the version v0.2.0
    # Will be removed in March 2026
    parser.add_argument(
        "-msdi",
        "--make-sd-image",
        nargs=3,
        type=int,
        default=None,
        help="make_sd_image is deprecated, image is made by default, use --no-html to suppress. This arguments will be removed from magnopy in March of 2026",
    )

    return parser
