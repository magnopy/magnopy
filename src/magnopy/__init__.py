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


__version__ = "0.2.1"
__doclink__ = "magnopy.org"
__release_date__ = "10 September 2025"

from . import _constants, examples, io, scenarios, experimental
from ._diagonalization import *
from ._energy import *
from ._exceptions import *
from ._local_rf import *
from ._lswt import *
from ._package_info import *
from ._parallelization import *
from ._parameters import *
from ._spinham import *
from ._plotly_engine import *
