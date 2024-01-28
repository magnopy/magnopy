# MAGNOPY - Python package for magnons.
# Copyright (C) 2023-2024 Magnopy Team
#
# e-mail: anry@uv.es, web: magnopy.com
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

import logging

import numpy as np
import wulfric
from termcolor import colored
from wulfric import Atom, print_2d_array
from wulfric.constants import TORADIANS

from magnopy.io.verify import verify_model_file

_logger = logging.getLogger(__name__)

from magnopy._pinfo import logo
from magnopy.spinham.hamiltonian import SpinHamiltonian
from magnopy.spinham.parameter import ExchangeParameter
from magnopy.units import (
    ANGSTROM,
    BOHR,
    ELECTRON_VOLT,
    JOULE,
    KELVIN,
    MILLI_ELECTRON_VOLT,
    RYDBERG,
    TRUE_KEYWORDS,
)

__all__ = ["load_model", "dump_model"]


SEPARATOR = "=" * 80
SUBSEPARATOR = "-" * 80


def _write_bond(name1: str, name2: str, ijk: tuple, J: ExchangeParameter):
    r"""
    Write a bond to the output file.

    Parameters
    ----------
    name1 : str
        Name of the first atom in the bond.
    name2 : str
        Name of the second atom in the bond.
    ijk : (3,) tuple
        Index of the bond.
    J : :py:class:`.ExchangeParameter`
        Exchange parameter for the bond.
    """

    text = [SUBSEPARATOR]
    text.append(
        f"{name1:<6} "
        + f"{name2:<6} "
        + f"{ijk[0]:>3} {ijk[1]:>3} {ijk[2]:>3} "
        + f"{J.iso:>8.4}"
    )
    text.append("Matrix")
    text.append(print_2d_array(J.matrix, fmt="8.4f", print_result=False, borders=False))

    return "\n".join(text)


def dump_model(spinham: SpinHamiltonian, filename):
    """
    Dump a SpinHamiltonian object to a .txt file.

    Parameters
    ----------
    spinham : :py:class`.SpinHamiltonian`
        SpinHamiltonian object to dump.
    filename : str
        Filename to dump SpinHamiltonian object to.
    """

    # Write logo
    text = [logo(comment=True)]
    text.append(SEPARATOR)

    # Write unit cell
    text.append("Cell: Angstrom")
    text.append(f"# {'x':>10} {'y':>12} {'z':>12}")
    text.append(
        wulfric.print_2d_array(
            spinham.cell,
            fmt="12.8f",
            print_result=False,
            borders=False,
            footer_column=[
                "# <- a (first lattice vector)",
                "# <- b (second lattice vector)",
                "# <- c (third lattice vector)",
            ],
        )
    )
    text.append(SEPARATOR)

    # Write atom's position and spins
    text.append("Atoms: relative")
    text.append(f"# Name {'x':>10} {'y':>12} {'z':>12}")
    for atom in spinham.atoms:
        text.append(
            f"{atom.name:4} "
            + f"{atom.position[0]:12.8f} "
            + f"{atom.position[1]:12.8f} "
            + f"{atom.position[2]:12.8f} "
        )
        try:
            text[-1] += (
                f"{atom.spin_vector[0]:8.4f} "
                + f"{atom.spin_vector[1]:8.4f} "
                + f"{atom.spin_vector[2]:8.4f}"
            )
        except ValueError:
            pass
    text.append(SEPARATOR)

    # Write notation
    text.append("Notation: ")
    try:
        text.append(f"Spin normalized: {spinham.spin_normalized}")
        text.append(f"Double counting: {spinham.double_counting}")
        text.append(f"Factor: {spinham.factor}")
        text.append("# Latex-styled:")
        text.append("# " + "\n# ".join(spinham.notation_string.split("\n")))
    except ValueError:
        text.append("Undefined")
    text.append(SEPARATOR)

    # Define which atom's names are unique
    names = [atom.name for atom in spinham.atoms]
    unique_names = [
        atom.name if names.count(atom.name) == 1 else atom.fullname
        for atom in spinham.atoms
    ]
    unique_names = dict(zip(spinham.atoms, unique_names))

    # Write exchange parameters
    text.append("Parameters: meV")
    text.append(f"{'# Atom1 Atom2':<13} {'i':>3} {'j':>3} {'k':>3} {'Jiso':>8}")
    for atom1, atom2, (i, j, k), J in spinham:
        text.append(_write_bond(unique_names[atom1], unique_names[atom2], (i, j, k), J))
    text.append(SEPARATOR)

    with open(filename, "w", encoding="utf-8") as f:
        f.write("\n".join(text))


def _read_cell(lines, spinham: SpinHamiltonian):
    R"""
    Read information from the cell section as described in the documentation.

    Input lines have to follow the format::

        Cell <Units> <Scale>
        a_x a_y a_z
        b_x b_y b_z
        c_x c_y c_z

    where optional keywords are:

    * <Units> - starts either from "b" or "a".
    * <Scale> - either one float or three floats separated by at least one space.

    Parameters
    ==========
    lines : (4,) list of str
        Cell section of the input file.
    spinham : :py:class:`.SpinHamiltonian`
        Spin Hamiltonian, where the data are saved.

    Returns
    =======
    spinham : :py:class:`.SpinHamiltonian`
        Spin Hamiltonian with updated cell.
    """

    # Cell <Units> <Scale>
    # or
    # Cell <Scale> <Units>
    line = lines[0].lower().split()

    # Search for the <Units> keyword
    if len(line) == 1:
        _logger.info("No <Units> nor <Scale> keywords are detected.")
        units = "angstroms"
    if len(line) >= 2:
        # <Units> keyword can be first in the keyword list. Only the <Units> keyword is alphabetic.
        if str.isalpha(line[1]):
            units = line.pop(1).lower()
            _logger.info(f'"{units}" keyword is detected.')
        # <Units> keyword can be last in the keyword list. Only the <Units> keyword is alphabetic.
        elif str.isalpha(line[-1]):
            units = line.pop(-1).lower()
            _logger.info(f'"{units}" keyword is detected.')
        else:
            _logger.info("No <Units> keywords is detected.")
            units = "angstroms"

    # Process <Units> keyword
    # Only those two cases are possible, since the input file is pre-verified
    if units.startswith("b"):
        units_conversion = BOHR
    elif units.startswith("a"):
        units_conversion = ANGSTROM

    # Process <Scale> keyword
    # Cell <s>
    # or
    # Cell <s_x s_y s_z>
    if len(line) == 2:
        scale_factors = [float(line[1]) for _ in range(3)]
    elif len(line) == 4:
        scale_factors = [float(line[i]) for i in range(1, 4)]
    else:
        scale_factors = [1.0, 1.0, 1.0]

    _logger.info(
        f"Units conversion factor: {units_conversion}; scale factors: {' '.join([str(f) for f in scale_factors])}"
    )

    cell = np.zeros((3, 3), dtype=float)
    for i in range(1, 4):
        line = lines[i]
        cell[i - 1] = [float(x) for x in line.split()[:3]]

    cell *= units_conversion
    for i in range(3):
        cell[i] *= scale_factors[i]

    spinham.cell = cell

    return spinham


def _read_atoms(lines, spinham: SpinHamiltonian):
    R"""
    Read information from the atoms section as described in the documentation.

    Input lines have to follow the format::

        Atoms <Units>
        A1 i j k <spin1>
        ...

    where optional keywords are:

    * <Units> - starts either from "r" or "a" or "b".
    * <spin> - either one float or three floats or four
    floats separated by at least one space. Alternatively it can have the format
    "p<float> t<float> float"

    Parameters
    ==========
    lines : (N,) list of str
        Atoms section of the input file.
    spinham : :py:class:`.SpinHamiltonian`
        Spin Hamiltonian, where the data are saved.

    Returns
    =======
    spinham : :py:class:`.SpinHamiltonian`
        Spin Hamiltonian with added atoms.
    """

    line = lines[0].lower().split()

    # Search for the <Units>
    if len(line) == 2:
        units = line[1]
        _logger.info(f'"{units}" keyword is detected.')
    else:
        units = "relative"
        _logger.info(
            f"No <Units> keyword is detected. Default value (relative) is used."
        )

    # Decide the case based on <Units>
    # Only three cases are possible, since the input lines are verified.
    if units.startswith("r"):
        relative = True
        units_conversion = ANGSTROM
    elif units.startswith("b"):
        relative = False
        units_conversion = BOHR
    elif units.startswith("a"):
        relative = False
        units_conversion = ANGSTROM

    _logger.info(
        f"Units conversion factor: {units_conversion}; coordinates are relative: {relative}"
    )

    # Read atoms's information
    # A i j k <spin>
    for line in lines[1:]:
        line = line.split()
        label = line[0]
        coordinates = np.array([float(x) for x in line[1:4]])
        spin_vector = None

        # Only spin value is given
        if len(line) == 5:
            spin_vector = [0, 0, float(line[4])]

        # Either spin vector is given or theta, phi and spin value
        if len(line) == 7:
            spin_data = " ".join(line[4:7])
            # Atom i j k p<angle> t<angle> S
            # or
            # Atom i j k p<angle> S t<angle>
            # or
            # Atom i j k t<angle> p<angle> S
            # or
            # Atom i j k t<angle> S p<angle>
            # or
            # Atom i j k S p<angle> t<angle>
            # or
            # Atom i j k S t<angle> p<angle>
            if "p" in spin_data:
                spin_data = spin_data.split()
                for i in range(3):
                    if spin_data[i].lower().startswith("p"):
                        phi = float(spin_data[i][1:]) * TORADIANS
                    elif spin_data[i].lower().startswith("t"):
                        theta = float(spin_data[i][1:]) * TORADIANS
                    else:
                        S = float(spin_data[i])
                    spin_vector = [
                        S * np.cos(phi) * np.sin(theta),
                        S * np.sin(phi) * np.sin(theta),
                        S * np.cos(theta),
                    ]
            # Atom i j k Sx Sy Sz
            else:
                spin_vector = np.array([float(x) for x in line[4:7]])
        # Spin direction and spin value is given
        # A i j k sdirx sdiry sdirz S
        if len(line) == 8:
            spin_vector = np.array([float(x) for x in line[4:7]])
            spin_vector = spin_vector / np.linalg.norm(spin_vector) * float(line[7])

        # Add atom to the Hamiltonian
        if spin_vector is None:
            spinham.add_atom(
                new_atom=label,
                position=coordinates * units_conversion,
                relative=relative,
            )
        else:
            spinham.add_atom(
                new_atom=label,
                position=coordinates * units_conversion,
                spin=spin_vector,
                relative=relative,
            )

    return spinham


def _read_parameters(lines, spinham):
    R"""
    Read information from the parameter section as described in the documentation.

    Input lines have to follow the format::

        Parameters <Units>
        ----------
        bond
        ----------
        ...

    where optional keywords are:

    * <Units> - starts either from "m" or "e" or "k" or "j" or "r".

    Parameters
    ==========
    lines : (N,) list of str
        Parameters section of the input file.
    spinham : :py:class:`.SpinHamiltonian`
        Spin Hamiltonian, where the data are saved.

    Returns
    =======
    spinham : :py:class:`.SpinHamiltonian`
        Spin Hamiltonian with added bonds.
    """

    # Search for the <Units>
    # Parameters <Units>
    line = lines[0].lower().split()
    if len(line) >= 2:
        units = line[1]
        _logger.info(f'"{units}" keyword is detected.')
    else:
        units = "meV"
        _logger.info(f"No <Units> keyword is detected. Default value (meV) is used.")

    # Decide the case based on <Units>
    # Only those cases are possible, since the input lines are verified.
    if units.startswith("r"):
        units_conversion = RYDBERG
        _logger.info(
            "Parameters are provided in Rydberg energy units. Will be converted to meV."
        )
    elif units.startswith("j"):
        units_conversion = JOULE
        _logger.info(f"Parameters are provided in Joule. Will be converted to meV.")
    elif units.startswith("k"):
        units_conversion = KELVIN
        _logger.info(f"Parameters are provided in Kelvin. Will be converted to meV.")
    elif units.startswith("e"):
        units_conversion = ELECTRON_VOLT
        _logger.info(
            f"Parameters are provided in electron-Volts. Will be converted to meV."
        )
    elif units.startswith("m"):
        units_conversion = MILLI_ELECTRON_VOLT
        _logger.info(f"Parameters are provided in meV.")

    _logger.info(f"Units conversion factor: {units_conversion}")

    # Skip first line with the section header
    i = 1
    while i < len(lines):
        # Skip subsection separators
        while i < len(lines) and lines[i].startswith("-" * 10):
            i += 1

        # Check if we reached the end of the file
        if i >= len(lines):
            break

        # Detect the beginning and end of the bond data
        bond_start = i
        while i < len(lines) and not lines[i].startswith("-" * 10):
            i += 1
        bond_end = i

        # Read bond and add it to the Hamiltonian
        bond = _read_bond(
            lines[bond_start:bond_end], spinham, units_conversion=units_conversion
        )

    return spinham


def _read_bond(lines, spinham: SpinHamiltonian, units_conversion=1):
    R"""
    Read information from the bond subsection as described in the documentation.

    Input lines have to follow the format::

        A1 A2 i j k <Isotropic parameter>
        <Matrix
        Jxx Jxy Jxz
        Jyx Jyy Jyz
        Jzx Jzy Jzz>
        <Symmetric anisotropy
        Jxx Jxy Jxz
        Jyx Jyy Jyz
        Jzx Jzy Jzz>
        <DMI Dx Dy Dz>

    Parameters
    ==========
    lines : (N,) list of str
        Parameters section of the input file.
    spinham : :py:class:`.SpinHamiltonian`
        Spin Hamiltonian, where the data are saved.

    Returns
    =======
    spinham : :py:class:`.SpinHamiltonian`
        Spin Hamiltonian with added bond.
    """
    line = lines[0].split()
    label1, label2 = line[:2]
    atom1 = spinham.get_atom(label1)
    atom2 = spinham.get_atom(label2)
    R = tuple([int(x) for x in line[2:5]])
    if len(line) == 6:
        iso = float(line[5])
    else:
        iso = None

    matrix = None
    symm = None
    dmi = None
    i = 1
    while i < len(lines):
        if lines[i].lower().startswith("m"):
            matrix = np.zeros((3, 3), dtype=float)
            for j in range(3):
                i += 1
                matrix[j] = [float(x) for x in lines[i].split()]
        if lines[i].lower().startswith("d"):
            print(lines[i])
            dmi = [float(x) for x in lines[i].split()[1:]]
        parameter = ExchangeParameter()
        if matrix is not None:
            parameter.matrix = matrix
        if iso is not None:
            parameter.iso = iso
        if dmi is not None:
            parameter.dmi = dmi
        if symm is not None:
            parameter.aniso = symm
        spinham.add_bond(atom1, atom2, R, J=parameter)
        i += 1
    return spinham


def _read_notation(lines, spinham):
    # Skip first line with the section header
    i = 1
    while i < len(lines):
        line = lines[i]
        # Whether spins are normalized
        if line.startswith("s"):
            spinham.spin_normalized = (
                line.split("=")[1].strip().lower() in TRUE_KEYWORDS
            )
        # Whether double counting is present
        elif line.startswith("d"):
            spinham.double_counting = (
                line.split("=")[1].strip().lower() in TRUE_KEYWORDS
            )
        # Exchange factor
        elif line.startswith("e"):
            spinham.exchange_factor = float(line.split("=")[1])
        # On-site factor
        elif line.startswith("o"):
            spinham.on_site_factor = float(line.split("=")[1])
        i += 1
    return spinham


def filter_model_file(filename, save_filtered=False):
    R"""
    Filter out all comments and blank lines from the model input file.

    Parameters
    ==========
    filename : str
        Path to the file.
    save_filtered : bool, default False
        Whether to save filtered copy as a separate file.
        A name is the same as of the original file with ".filtered" added to the end.

    Returns
    =======
    filtered_lines : (N,) list of str
        Content of the file without comments and blank lines.
    lines_indices : (N,) list
        Indices of filtered lines in the original file.
    """
    # Read the content of the file
    with open(filename, "r", encoding="utf-8") as f:
        lines = f.readlines()

    # Filter comments and blank lines
    filtered_lines = []
    line_indices = []
    for l_i, line in enumerate(lines):
        if line.startswith("#"):
            continue
        line = line.strip().split("#")[0]
        if line:
            filtered_lines.append(line)
            line_indices.append(l_i + 1)

    if save_filtered:
        with open(filename + ".filtered", "w", encoding="utf-8") as f:
            f.write("\n".join(filtered_lines))
        _logger.debug(
            f"Filtered input file is saved in {os.path.abspath(filename + '.filtered')}"
        )

    return filtered_lines, line_indices


def load_model(filename, save_filtered=False, verbose=False) -> SpinHamiltonian:
    r"""
    Load a SpinHamiltonian object from a .txt file.

    Parameters
    ----------
    filename : str
        Filename to load SpinHamiltonian object from.
    save_filtered : bool
        Whether to save the pre-processed file.
    verbose : bool, default False
        Whether to output verbose comments on the progress.

    Returns
    -------
    spinham :py:class:`.SpinHamiltonian`
        SpinHamiltonian object loaded from file.
    """

    lines, indices = filter_model_file(filename=filename, save_filtered=save_filtered)

    # Verify input
    sections = verify_model_file(
        lines, indices, raise_on_fail=True, return_sections=True
    )

    # Construct spin Hamiltonian:
    spinham = SpinHamiltonian()

    _read_cell(lines[slice(*sections["cell"])], spinham)

    _read_atoms(lines[slice(*sections["atoms"])], spinham)

    _read_parameters(lines[slice(*sections["parameters"])], spinham)

    _read_notation(lines[slice(*sections["notation"])], spinham)

    return spinham


if __name__ == "__main__":
    # logging.basicConfig(filename="log.log", level=logging.INFO)
    logging.info("Started")
    model = load_model("../../../tmp/test.txt", verbose=True)
    logging.info("Finished")
