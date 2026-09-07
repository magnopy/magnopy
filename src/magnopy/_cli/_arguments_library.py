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
import sys
from magnopy._constants._icons import ICON_IN_ARG


def _get_command_info(args, hide_personal_data=False):
    if hide_personal_data:
        comment = f"Executed the command\n\n  magnopy-optimize-sd {' '.join(sys.argv[1:])}\n\n"
    else:
        comment = f"Executed the command\n\n  {' '.join(sys.argv)}\n\n"

    # Save arguments to the comment
    comment += "which resulted in argument values\n\n" + "\n".join(
        [f"{ICON_IN_ARG} {key:20} : {value}" for key, value in vars(args).items()]
    )

    return comment


def _add_spinham_input(parser):
    parser.add_argument(
        "-ss",
        "--spinham-source",
        type=str,
        # metavar="KEYWORD",
        required=True,
        choices=["GROGU", "TB2J", "tb2j", "grogu"],
        help="Magnopy supports two sources of spin Hamiltonian: GROGU or TB2J. Source of "
        "the spin Hamiltonian can not be deduced reliably from the Hamiltonian file. Thus, "
        "Magnopy asks you to specify what was the source of the file that you pass to "
        '"--spinham-filename".',
    )
    parser.add_argument(
        "-sf",
        "--spinham-filename",
        type=str,
        metavar="FILENAME",
        required=True,
        help="Path to the file with the parameters of spin Hamiltonian. If you use "
        '"--spinham-source TB2J", then this is an "exchange.out" file, produced by TB2J. '
        'If you use "--spinham-source GROGU", then this is a .txt file with the '
        "Hamiltonian parameters produced by GROGU.",
    )


def _add_spin_values(parser):
    parser.add_argument(
        "-sv",
        "--spin-values",
        nargs="*",
        type=str,
        metavar=("S1",),
        default=None,
        help='Only works with "--spinham-source TB2J". This option changes spin values of '
        "the magnetic atoms. It expects M numbers in the same order as the order of "
        "magnetic atoms in the input TB2J file. Magnetic atoms are defined as the atoms "
        "that have at least one interaction associated with them. Spin values are always "
        "positive. If you want to specify spin-down state or AFM order use spin directions "
        "and not the spin values.",
    )


def _add_magnetic_field(parser):
    parser.add_argument(
        "-mf",
        "--magnetic-field",
        default=None,
        nargs=3,
        metavar=("B_X", "B_Y", "B_Z"),
        type=float,
        help="Vector of the external magnetic field (magnetic flux density, B), given in "
        "the units of Tesla.",
    )


def _add_output_folder(parser):
    parser.add_argument(
        "-of",
        "--output-folder",
        type=str,
        default="magnopy-results",
        help="Path to the folder where all output files of Magnopy wil be saved. Defaults "
        'to "magnopy-results". If the folder does not exist, then Magnopy tries to '
        "create it.",
    )


def _add_no_html(parser):
    parser.add_argument(
        "-no-html",
        "--no-html",
        action="store_true",
        default=False,
        help=".html files are generally heavy (~ 5 Mb). This option allows to disable "
        "their generation to save disk space.",
    )


def _add_hide_personal_data(parser):
    parser.add_argument(
        "-hpd",
        "--hide-personal-data",
        action="store_true",
        default=False,
        help="Whether not to use os.path.abspath. Hides the file structure of your "
        "personal computer. Note, that the relative paths can still expose the "
        "structure of you filesystem.",
    )


def _add_supercell(parser):
    parser.add_argument(
        "-s",
        "--supercell",
        nargs=3,
        type=int,
        default=(1, 1, 1),
        metavar=("xA_1", "xA_2", "xA_3"),
        help="Definition of the supercell for the spin optimization. Expects three "
        "integers as an input, each integer is a multiple of the original unit cell "
        "along one of the three lattice vectors. By default optimizes on the unit cell.",
    )


def _add_energy_tolerance(parser):
    parser.add_argument(
        "-et",
        "--energy-tolerance",
        default=1e-5,
        type=float,
        help="Tolerance parameter. Difference between classical energies of two "
        "consecutive optimization steps.",
    )


def _add_torque_tolerance(parser):
    parser.add_argument(
        "-tt",
        "--torque-tolerance",
        default=1e-5,
        type=float,
        help="Tolerance parameter. Maximum value of torque among all spins at current "
        "optimization step.",
    )


def _add_spin_directions(parser):
    parser.add_argument(
        "-sd",
        "--spin-directions",
        nargs="*",
        type=str,
        default=None,
        metavar=("FILENAME",),
        help="Set of spin directions for every magnetic site in the unit cell defines "
        "the ground state of the Hamiltonian. Magnopy can try to identify the ground "
        "state on the fly via energy minimization (see also magnopy-optimize-sd). "
        "However, we recommend providing the set of spin directions explicitly. This "
        "argument expects a filename for the file that contains three numbers per line "
        "and has M lines in total. M is the number of magnetic sites. Alternatively, you "
        "can give a set of 3*M numbers directly to this argument. Then, first three "
        "numbers define spin direction of the first magnetic site, next three numbers - "
        "of the second and so on. Magnetic sites are defined as atoms that have at least "
        "one parameter associated with them.",
    )


def _add_k_path(parser):
    parser.add_argument(
        "-kp",
        "--k-path",
        default=None,
        metavar="GAMMA-X-S|GAMMA-Y",
        type=str,
        help="Path in reciprocal space between high symmetry points. Used for the plots "
        'of magnon dispersion and other quantities. Ignored if "--kpoints FILENAME" is '
        "used. List of available high symmetry points is calculated automatically by "
        "wulfric package using convention of HPKOT paper based on the space group of the "
        "crystal structure (see wulfric.org for more details). Numerical tolerance for "
        "the space group search via spglib package can be controlled with "
        '"--spglib-symprec NUMBER" argument. Note: high symmetry points are defined by '
        "wulfric based on the reciprocal lattice of the primitive cell, but their "
        "relative coordinates are given in the basis of reciprocal cell of the input "
        "cell (i. e. the same unit cell as in the TB2J or GROGU file). Primitive cell "
        "might be different from the input cell. If you want to control what atoms are "
        "considered to be equivalent in the space group search - manually change "
        "their names in the input file. To learn how a name is translated into the "
        "spglib_type see documentation of the wulfric package (wulfric.org).",
    )


def _add_k_points(parser):
    parser.add_argument(
        "-kps",
        "--kpoints",
        type=str,
        default=None,
        help="Explicit list of k-points in reciprocal space. A filename is expected. The "
        "file shall contain three numbers per line. Each line specifies one k-point. By "
        "default three numbers are understood as absolute coordinates in reciprocal "
        'space. If "--relative" is used, than three numbers are understood as relative '
        "coordinates in the basis of the reciprocal cell of the input unit cell (i. e. "
        "of the same unit cell as in the TB2J or GROGU file, not of the primitive cell).",
    )


def _add_relative(parser):
    parser.add_argument(
        "-r",
        "--relative",
        default=False,
        action="store_true",
        help='Relevant only if "--kpoints FILENAME" is used. See description of '
        '"--kpoints".',
    )


def _add_number_processors(parser):
    parser.add_argument(
        "-np",
        "--number-processors",
        type=int,
        default=None,
        help="Magnopy is parallelized over the k-points. By default it uses all "
        "available processors. Pass 1 to run in serial.",
    )


def _add_spglib_symprec(parser):
    parser.add_argument(
        "-spg-s",
        "--spglib-symprec",
        metavar="NUMBER",
        type=float,
        default=1e-5,
        help="Tolerance parameter for the space group symmetry search by spglib.",
    )


def _add_spglib_types(parser):
    parser.add_argument(
        "-spg-t",
        "--spglib-types",
        metavar=("TYPE1", "TYPE2"),
        nargs="*",
        type=int,
        default=None,
        help="Set of spglib types for the space group search by spglib. Expects N non-zero "
        "integers in total. N is the amount of all sites/atoms in the input files (not"
        "just magnetic, but all).  If not given, then guesses automatically by wulfric "
        "(see wulfric.get_spglib_types() at wulfric.org types for more details).",
    )


def _add_quantum_correction(parser):
    parser.add_argument(
        "-qc",
        "--quantum-correction",
        default=False,
        action="store_true",
        help="Whether to include quantum correction to the classical energy for the "
        "optimization of spin directions. This parameter is experimental. It will be "
        "improved and tested in future versions. Use with caution.",
    )
