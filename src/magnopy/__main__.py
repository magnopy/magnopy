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


from argparse import ArgumentParser, RawDescriptionHelpFormatter

from magnopy import __version__
from magnopy._package_info import _warranty, logo


def main():
    parser = ArgumentParser(
        description=logo() + "\n\nAvailable scripts are:\n\n"
        "* magnopy-optimize-sd\n\n"
        "* magnopy-lswt\n\n"
        "To call for help for each script type <script name> --help\n"
        "Information below is relevant only to 'magnopy' command.",
        formatter_class=RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "commands",
        default=None,
        help="command/commands on what to do. Use to display some information about package.",
        choices=["logo", "warranty"],
        nargs="*",
    )
    parser.add_argument(
        "-v",
        "--version",
        action="store_true",
        help="print version of magnopy",
    )

    args = parser.parse_args()

    if args.version:
        print(f"magnopy v{__version__}")

    if len(args.commands) == 0:
        parser.print_help()
        return

    for command in args.commands:
        if command == "logo":
            print(logo())
        elif command == "warranty":
            print("\n" + _warranty() + "\n")
        else:
            raise ValueError(f"Sub-command {args.command} is not recognized.")


if __name__ == "__main__":
    main()
