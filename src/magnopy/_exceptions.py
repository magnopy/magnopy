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


R"""Exceptions"""

# Save local scope at this moment
old_dir = set(dir())
old_dir.add("old_dir")


class ColpaFailed(Exception):
    r"""
    Raised when diagonalization via Colpa fails.
    """

    def __init__(self):
        message = "Diagonalization via Colpa failed."
        super().__init__(message)


class ConventionError(Exception):
    r"""
    Raised if convention or part of the convention of spin Hamiltonian is not defined.
    """

    def __init__(self, convention, property):
        message = (
            f"Convention of spin Hamiltonian has an undefined property '{property}':\n"
            + convention.summary(return_as_string=True)
        )
        super().__init__(message)


# Populate __all__ with objects defined in this file
__all__ = list(set(dir()) - old_dir)
# Remove all semi-private objects
__all__ = [i for i in __all__ if not i.startswith("_")]
del old_dir
