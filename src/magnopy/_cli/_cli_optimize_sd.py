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
import sys
from argparse import ArgumentParser, RawDescriptionHelpFormatter

from magnopy._package_info import logo
from magnopy.io._grogu import load_grogu
from magnopy.io._tb2j import load_tb2j
from magnopy.scenarios._optimize_sd import optimize_sd
from magnopy._cli._arguments_library import (
    _get_command_info,
    _add_spinham_input,
    _add_spin_values,
    _add_magnetic_field,
    _add_output_folder,
    _add_no_html,
    _add_hide_personal_data,
    _add_supercell,
    _add_energy_tolerance,
    _add_torque_tolerance,
    _add_quantum_correction,
)
from magnopy._constants._icons import ICON_IN_FILE


def manager():
    # Configure a parser
    parser = ArgumentParser(
        description=logo(logo_width=80)
        + "\n\nThis script optimizes classical energy of the spin Hamiltonian and "
        "finds spin directions that define a local minimum of the energy landscape.",
        formatter_class=RawDescriptionHelpFormatter,
    )

    _add_spinham_input(parser=parser)
    _add_energy_tolerance(parser=parser)
    _add_torque_tolerance(parser=parser)
    _add_supercell(parser=parser)
    _add_magnetic_field(parser=parser)
    _add_output_folder(parser=parser)
    _add_spin_values(parser=parser)
    _add_no_html(parser=parser)
    _add_hide_personal_data(parser=parser)
    _add_quantum_correction(parser=parser)

    # Parse arguments
    args = parser.parse_args()

    # Save executed command and arguments into a comment
    comment = _get_command_info(args=args, hide_personal_data=args.hide_personal_data)

    # Handle execution with no arguments
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    # Process spin values
    if args.spin_values is not None:
        args.spin_values = [float(tmp) for tmp in args.spin_values]

    # Load spin Hamiltonian
    if args.spinham_source.lower() == "tb2j":
        spinham = load_tb2j(
            filename=args.spinham_filename, spin_values=args.spin_values
        )
    elif args.spinham_source.lower() == "grogu":
        spinham = load_grogu(
            filename=args.spinham_filename, spin_values=args.spin_values
        )
    else:
        raise ValueError(
            'Supported sources of spin Hamiltonian are "GROGU" and "TB2J", '
            f'got "{args.spinham_source}".'
        )

    if args.hide_personal_data:
        spinham_filename = args.spinham_filename
    else:
        spinham_filename = os.path.abspath(args.spinham_filename)

    comment += f'\n\nParameters are loaded from "{args.spinham_source.upper()}", source file\n{ICON_IN_FILE} {spinham_filename}'

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
        quantum_correction=args.quantum_correction,
    )
