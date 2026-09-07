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

import numpy as np

from magnopy._package_info import logo
from magnopy.io._grogu import load_grogu
from magnopy.io._spin_directions import read_spin_directions
from magnopy.io._tb2j import load_tb2j
from magnopy.scenarios._solve_lswt import solve_lswt
from magnopy._cli._arguments_library import (
    _get_command_info,
    _add_spinham_input,
    _add_spin_values,
    _add_magnetic_field,
    _add_output_folder,
    _add_no_html,
    _add_hide_personal_data,
    _add_spin_directions,
    _add_k_path,
    _add_k_points,
    _add_relative,
    _add_spglib_symprec,
    _add_number_processors,
    _add_spglib_types,
)
from magnopy._constants._icons import ICON_IN_FILE


def manager():
    # Configure a parser
    parser = ArgumentParser(
        description=logo(logo_width=80)
        + "\n\nThis script solves the spin Hamiltonian at the level of "
        "Linear Spin Wave Theory (LSWT).",
        formatter_class=RawDescriptionHelpFormatter,
    )

    # Add common arguments
    _add_spinham_input(parser=parser)
    _add_spin_directions(parser=parser)
    _add_magnetic_field(parser=parser)
    _add_output_folder(parser=parser)
    _add_k_path(parser=parser)
    _add_k_points(parser=parser)
    _add_relative(parser=parser)
    _add_spglib_symprec(parser=parser)
    _add_number_processors(parser=parser)
    _add_spin_values(parser=parser)
    _add_no_html(parser=parser)
    _add_hide_personal_data(parser=parser)
    _add_spglib_types(parser=parser)

    args = parser.parse_args()

    # Save executed command and arguments into a comment
    comment = _get_command_info(args=args, hide_personal_data=args.hide_personal_data)

    # Handle execution with no arguments
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
        with open(args.kpoints, "r", encoding="utf-8") as f:
            for i, line in enumerate(f):
                # Remove comment lines
                if line.startswith("#"):
                    continue
                # Remove inline comments and leading/trailing whitespaces
                line = line.split("#")[0].strip()
                # Check for empty lines
                if line:
                    line = line.split()
                    if len(line) != 3:
                        raise ValueError(
                            f"Expected three numbers per line (in line {i}),"
                            f"got: {len(line)}."
                        )

                    kpoints.append(list(map(float, line)))

        args.kpoints = kpoints

    # Load spin Hamiltonian
    if args.spinham_source.lower() == "tb2j":
        spinham = load_tb2j(
            filename=args.spinham_filename,
            spin_values=args.spin_values,
            spglib_types=args.spglib_types,
        )
    elif args.spinham_source.lower() == "grogu":
        spinham = load_grogu(
            filename=args.spinham_filename,
            spin_values=args.spin_values,
            spglib_types=args.spglib_types,
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
