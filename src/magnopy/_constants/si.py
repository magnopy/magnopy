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


R"""
Constants and custom units in the International system of units.
"""

# Save local scope at this moment
old_dir = set(dir())
old_dir.add("old_dir")

################################################################################
#                                    Units                                     #
################################################################################
# Length, in Metre
ANGSTROM = 1e-10
# Energy, in Joule
ELECTRON_VOLT = 1.602176634e-19
# Time, in Second
SECOND = 1
# Temperature, in Kelvin
KELVIN = 1
# Magnetic field, in Tesla
TESLA = 1


################################################################################
#                                  Constants                                   #
################################################################################
K_BOLTZMANN = 1.380649e-23  # Joule / Kelvin
BOHR_RADIUS = 5.29177210903e-11  # Metre
BOHR_MAGNETON = 9.2740100783e-24  # Joule / Tesla
PLANK_CONSTANT = 6.62607015e-34  # Joule * Second
SPEED_OF_LIGHT = 299792458  # Metre / Second
RYDBERG_CONSTANT = 10973731.568160  # 1 / Metre
RYDBERG_ENERGY = PLANK_CONSTANT * SPEED_OF_LIGHT * RYDBERG_CONSTANT  # Joule
ELEMENTARY_CHARGE = 1.602176634e-19  # Coulomb

# Populate __all__ with objects defined in this file
__all__ = list(set(dir()) - old_dir)
# Remove all semi-private objects
__all__ = [i for i in __all__ if not i.startswith("_")]
del old_dir
