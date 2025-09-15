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
import re
from argparse import ArgumentParser
from calendar import month_name
from datetime import datetime


class VersionError(Exception):
    pass


class MetadataError(Exception):
    pass


class ReleaseNotesError(Exception):
    pass


def envelope(message: str):
    """
    Decorator for printing a message before and "All good" after a function call.

    Parameters
    ----------
    message : str
        Message to print before the function.
    """

    def wrapper(func):
        def inner(*args, **kwargs):
            print(f"{message} ... ", end="")
            func(*args, **kwargs)
            print("All good")

        return inner

    return wrapper


@envelope(message="Checking __init__.py")
def check_init(version, root_dir: str):
    """
    Check the __release_date__ and __version__ variables.

    Parameters
    ----------
    version : str
        Target version for the release.
    """

    cd = datetime.now()

    variables = ["__release_date__", "__version__", "__doclink__"]
    expected_values = [
        f"{cd.day} {month_name[cd.month]} {cd.year}",
        version,
        "magnopy.org",
    ]
    found = [False, False, False]
    good = [True, True, True]

    # Read __init__.py
    with open(os.path.join(root_dir, "src", "magnopy", "__init__.py"), "r") as f:
        init_file_content = f.readlines()

    # Run trough the lines
    for line in init_file_content:
        for v_i, variable in enumerate(variables):
            # Check presence
            if line.startswith(variable):
                found[v_i] = True
                # Check value
                if line != f'{variable} = "{expected_values[v_i]}"\n':
                    good[v_i] = False

    # Check if all variables were present
    if not all(found):
        raise MetadataError(
            "".join(
                [
                    "\nFailed to find some variables in '__init__.py':\n",
                    "red",
                    *[
                        f"    {variables[i]:20} : {'found' if found[i] else 'not found'}\n"
                        for i in range(len(variables))
                    ],
                ]
            )
        )

    # Check that all values are correct
    if not all(good):
        raise MetadataError(
            "".join(
                [
                    "\nSome variables are not right in '__init__.py':\n",
                    *[
                        f"    {variables[i]:20} : {'good' if good[i] else f'Expected {expected_values[i]}'}\n"
                        for i in range(len(variables))
                    ],
                ]
            )
        )


@envelope(message="Checking release notes")
def check_release_notes(version: str, root_dir: str):
    """
    Check if the release notes are up to date.

    Parameters
    ----------
    version : str
        Target version for the release.
    """
    path = os.path.join(root_dir, "docs", "source", "release-notes")

    # (major, minor, micro)
    major, minor, micro = tuple(map(int, version.split(".")[:3]))

    for _, _, filenames in os.walk(path):
        break

    # (major, minor)
    files = []
    for filename in sorted(filenames):
        if re.fullmatch(R"[0-9]*\.[0-9]*\.rst", filename):
            files.append(tuple(map(int, filename.split(".")[:2])))

    # Check the minor version file
    if (major, minor) not in files:
        raise ReleaseNotesError(
            "\n".join(
                [
                    "",
                    f"Release notes for the minor version {major}.{minor} are not available",
                    "Please add the file",
                    "",
                    f"    {major}.{minor}.rst",
                    "",
                    "to the directory",
                    "",
                    "    docs/source/release-notes/",
                    "",
                ]
            )
        )

    # Check the toctree in the index file
    # Read index.rst
    found_pages = []
    with open(os.path.join(path, "index.rst"), "r") as f:
        for line in f:
            line = line.strip().split(".")

            if len(line) == 2:
                try:
                    found_pages.append(tuple(map(int, line)))
                except ValueError:
                    pass

    if (major, minor) not in found_pages:
        raise ReleaseNotesError(
            "\n".join(
                [
                    f'Did not find reference to "{major}.{minor}" in the toctree of "docs/source/release-notes/index.rst"',
                    "Found references to",
                    *[f"    {_[0]}.{_[1]}" for _ in found_pages],
                ]
            )
        )

    # Check the version note in the minor version file
    # For the x.x.0 versions there are no note required, but 'Whats new?' section
    with open(os.path.join(path, f"{major}.{minor}.rst"), "r") as f:
        found_note = False
        date_line = None
        for line in f:
            line = line.strip()
            if line == version:
                found_note = True

            if found_note and line.startswith("**Date**:"):
                date_line = line
                break

    if not found_note:
        raise ReleaseNotesError(
            "\n".join(
                [
                    "",
                    f"Release notes for the {version} version are not found. Please add the note for the {version} version to the file",
                    "",
                    f"    docs/source/release-notes/{major}.{minor}.rst",
                    "",
                    "Follow the style:",
                    "",
                    f"    {version}",
                    f"    {'':=^{len(version)}}",
                    "    # Note text",
                    "",
                ]
            )
        )

    if date_line is None:
        raise ReleaseNotesError(
            f'Did not find a date of the release in the release notes for the version {version} in the file "docs/source/release-notes/{major}.{minor}.rst"'
        )

    cd = datetime.now()
    if date_line != f"**Date**: {cd.day} {month_name[cd.month]} {cd.year}":
        raise ReleaseNotesError(
            f'Date is not correct or format is not correct. Expected "**Date**: {cd.day} {month_name[cd.month]} {cd.year}", got "{date_line}"'
        )


def main(version: str, root_dir: str):
    if version == "None":
        raise VersionError(
            "\n".join(
                [
                    "",
                    "Version is undefined",
                    "For the tags use the syntax:",
                    "",
                    "    v<major>.<minor>.<micro>",
                    "",
                    "For the make command use the syntax:",
                    "",
                    "    make prepare-release VERSION=v<major>.<minor>.<micro>",
                ]
            )
        )
    elif version[0] != "v" or len(version.split(".")) != 3:
        raise VersionError(
            f'Version shall has the form "v<major>.<minor>.<micro>", got {version}'
        )

    print(f"Checking for {version} release")

    version = version[1:]

    check_init(version=version, root_dir=root_dir)

    check_release_notes(version=version, root_dir=root_dir)


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument(
        "-v",
        "--version",
        metavar="vx.x.x",
        type=str,
        required=True,
        help="Version to release",
    )
    parser.add_argument(
        "-rd",
        "--root-dir",
        metavar="PATH",
        type=str,
        required=True,
        help="Path to the root directory of the project",
    )

    args = parser.parse_args()
    main(**vars(args))
