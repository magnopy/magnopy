usage: magnopy-optimize-sd [-h] -sf FILENAME -ss KEYWORD [-sv [S1 [S2 ...]]]
                           [-s xa_1 xa_2 xa_3] [-et ENERGY_TOLERANCE]
                           [-tt TORQUE_TOLERANCE] [-mf h_x h_y h_z]
                           [-of OUTPUT_FOLDER] [-no-html] [-hpd]
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

This script optimizes classical energy of the spin Hamiltonian and finds the spin directions that describe a local minima of the energy landscape.

options:
  -h, --help            show this help message and exit
  -sf, --spinham-filename FILENAME
                        Path to the spin Hamiltonian file, from where the
                        parameters would be read.
  -ss, --spinham-source KEYWORD
                        Source of the spin Hamiltonian. Either "GROGU" or
                        "TB2J"
  -sv, --spin-values [S1 [S2 ...]]
                        In the case when the parameters of spin Hamiltonian
                        comes from TB2J, one might want to change the values
                        of spins to be closer to half-integers. This option
                        allows that. Order of the M numbers should match the
                        order of magnetic atoms in the spin Hamiltonian. Note
                        that those numbers are always positive. To specify AFM
                        order use opposite spin directions and not spin values
                        of the opposite sign.
  -s, --supercell xa_1 xa_2 xa_3
                        Specification of the supercell for the spin
                        optimization. Expects three integers as an input. Pass
                        1 1 1 to optimize within the original unit cell.
  -et, --energy-tolerance ENERGY_TOLERANCE
                        Tolerance parameter. Difference between classical
                        energies of two consecutive optimization steps.
  -tt, --torque-tolerance TORQUE_TOLERANCE
                        Maximum torque among all spins.
  -mf, --magnetic-field h_x h_y h_z
                        Vector of external magnetic field, given in the units
                        of Tesla.
  -of, --output-folder OUTPUT_FOLDER
                        Folder where all output files of magnopy wil be saved.
  -no-html, --no-html   Disable plotting of the spin direction image in the
                        .html format. html files are generally heavy (~> 5
                        Mb). This option allows to disable their production to
                        save disk space.
  -hpd, --hide-personal-data
                        Whether to strip the parts of the paths as to hide the
                        file structure of you personal computer.
  -msdi, --make-sd-image MAKE_SD_IMAGE MAKE_SD_IMAGE MAKE_SD_IMAGE
                        make_sd_image is deprecated, image is made by default,
                        use --no-html to suppress. This arguments will be
                        removed from magnopy in March of 2026
