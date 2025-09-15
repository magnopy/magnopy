usage: magnopy-lswt [-h] -sf filename -ss name [-sd [S1_x S2_y S3_z ... ...]]
                    [-sv [S1 S2 S3 ... ...]] [-kp GAMMA-X-S|GAMMA-Y]
                    [-kps KPOINTS] [-r]
                    [-mf MAGNETIC_FIELD MAGNETIC_FIELD MAGNETIC_FIELD]
                    [-of OUTPUT_FOLDER] [-np NUMBER_PROCESSORS] [-no-html]
                    [-hpd] [-spg-s SPGLIB_SYMPREC]
                    [-msdi MAKE_SD_IMAGE MAKE_SD_IMAGE MAKE_SD_IMAGE]

███╗   ███╗  █████╗   ██████╗  ███╗   ██╗  ██████╗  ██████╗  ██╗   ██╗
████╗ ████║ ██╔══██╗ ██╔════╝  ████╗  ██║ ██╔═══██╗ ██╔══██╗ ╚██╗ ██╔╝
██╔████╔██║ ███████║ ██║  ███╗ ██╔██╗ ██║ ██║   ██║ ██████╔╝  ╚████╔╝
██║╚██╔╝██║ ██╔══██║ ██║  ╚██║ ██║╚██╗██║ ██║   ██║ ██╔═══╝    ╚██╔╝
██║ ╚═╝ ██║ ██║  ██║ ╚██████╔╝ ██║ ╚████║ ╚██████╔╝ ██║         ██║
╚═╝     ╚═╝ ╚═╝  ╚═╝  ╚═════╝  ╚═╝  ╚═══╝  ╚═════╝  ╚═╝         ╚═╝
                                                           ▄   ▄
                      Version: 0.2.0                       █▀█▀█
                Documentation: magnopy.org                 █▄█▄█
             Release date: 10 September 2025                ███   ▄▄
                    License: GNU GPLv3                      ████ █  █
          Copyright (C) 2023-2025  Magnopy Team             ████    █
                                                            ▀▀▀▀▀▀▀▀

This script solves the spin Hamiltonian at the level of Linear Spin Wave Theory (LSWT) and outputs (almost) every possible quantity.

options:
  -h, --help            show this help message and exit
  -sf, --spinham-filename filename
                        Path to the spin Hamiltonian file, from where the
                        parameters would be read.
  -ss, --spinham-source name
                        Source of the spin Hamiltonian. Either "GROGU" or
                        "TB2J"
  -sd, --spin-directions [S1_x S2_y S3_z ... ...]
                        To fully define the system for the calculations of
                        magnons one need the information about the ground
                        state in addition to the parameters of the
                        Hamiltonian. There are two ways to give this
                        information to magnopy: * Give a path to the file. In
                        the file there should be M lines with three numbers in
                        each. The order of the lines would match the order of
                        magnetic atoms in the spin Hamiltonian. * Give a
                        sequence of 3*M numbers directly to this parameter. If
                        none provided, then magnopy attempts to optimize the
                        spin directions prior to the LSWT calculations.
  -sv, --spin-values [S1 S2 S3 ... ...]
                        In the case when the parameters of spin Hamiltonian
                        comes from TB2J, one might want to change the values
                        of spins to be closer to half-integers. This option
                        allows that. Order of the M numbers should match the
                        order of magnetic atoms in the spin Hamiltonian. Note
                        that those numbers are always positive. To specify AFM
                        order use opposite spin directions and not spin values
                        of the opposite sign.
  -kp, --k-path GAMMA-X-S|GAMMA-Y
                        Path of high symmetry k-points for the plots of
                        dispersion and other quantities.
  -kps, --kpoints KPOINTS
                        Alternatively one could provide an explicit list of
                        k-points for calculation. In that case provide a path
                        to the file, in which each k-point is given in a
                        separate line with three numbers per line.
  -r, --relative        When an explicit list of k-points is given, this
                        option specify whether to consider them as relative or
                        absolute coordinates. Absolute by default.
  -mf, --magnetic-field MAGNETIC_FIELD MAGNETIC_FIELD MAGNETIC_FIELD
                        Vector of external magnetic field, given in the units
                        of Tesla.
  -of, --output-folder OUTPUT_FOLDER
                        Folder where all output files of magnopy wil be saved.
  -np, --number-processors NUMBER_PROCESSORS
                        Number of processes for multithreading. Uses all
                        available processors by default. Pass 1 to run in
                        serial.
  -no-html, --no-html   html files are generally heavy (~> 5 Mb). This option
                        allows to disable their production to save disk space.
  -hpd, --hide-personal-data
                        Whether to strip the parts of the paths as to hide the
                        file structure of you personal computer.
  -spg-s, --spglib-symprec SPGLIB_SYMPREC
                        Tolerance parameter for the space group symmetry
                        search by spglib.
  -msdi, --make-sd-image MAKE_SD_IMAGE MAKE_SD_IMAGE MAKE_SD_IMAGE
                        make_sd_image is deprecated, use --no-html instead.
                        This arguments will be removed from magnopy in March
                        of 2026
