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
from wulfric import PlotlyEngine as W_PlotlyEngine
import numpy as np
from random import choices
from string import ascii_lowercase as ASCII_LOWERCASE

# Save local scope at this moment
old_dir = set(dir())
old_dir.add("old_dir")


class PlotlyEngine(W_PlotlyEngine):
    r"""
    Plotting engine based on |plotly|_.

    Parameters
    ----------
    fig : plotly graph object, optional
        Figure to plot on. If not provided, then a new one is created as
        ``fig = go.Figure()``.
    _sphinx_gallery_fix : bool, default  False
        Fixes display issues when building documentation using sphinx gallery.
        Please, always ignore this argument

    Attributes
    ----------
    fig : plotly graph object
        Figure to plot on.

    Notes
    -----
    This class is a part of ``magnopy[visual]``
    """

    def plot_spin_directions(
        self,
        positions,
        spin_directions,
        colors="#000000",
        legend_label=None,
        legend_group=None,
        row=1,
        col=1,
    ):
        r"""
        Plots a set of spin directions.

        Parameters
        ----------
        positions : (N, 3) |array-like|_
            Positions of magnetic centers.
        spin_directions : (N, 3) |array-like|_
            Direction of spin vectors for each magnetic center. Only direction is used,
            the modulus is ignored. Spin directions are all normalized to the same length
            but not necessary to 1.
        colors : str or list of str, default "#000000"
            Color of the line. Any value that is supported by |plotly|_.
        legend_label : str, optional
            Label of the line that is displayed in the figure.
        legend_group : str, optional
            Legend's group. If ``None``, then defaults to the random string of 10
            characters.
        row : int, default 1
            Row of the subplot.
        col : int, default 1
            Column of the subplot.
        """

        if legend_group is None:
            legend_group = "".join(choices(ASCII_LOWERCASE, k=10))

        pos = np.array(positions, dtype=float)
        sd = np.array(spin_directions, dtype=float)

        # Get normalization length
        tmp = pos[:, np.newaxis] - pos[np.newaxis, :]
        tmp = np.linalg.norm(tmp, axis=2)
        norm_length = (tmp + np.eye(tmp.shape[0]) * tmp.max()).min()

        if norm_length == 0:
            norm_length = 1

        if isinstance(colors, str):
            colors = [colors for _ in range(len(pos))]

        sd = sd / np.linalg.norm(sd, axis=1)[:, np.newaxis] * norm_length

        for i in range(len(pos)):
            self.plot_vector(
                start_point=pos[i],
                end_point=pos[i] + sd[i],
                color=colors[i],
                vector_label=None,
                legend_label=legend_label,
                legend_group=legend_group,
                row=row,
                col=col,
            )
            if legend_label is not None:
                legend_label = None


# Populate __all__ with objects defined in this file
__all__ = list(set(dir()) - old_dir)
# Remove all semi-private objects
__all__ = [i for i in __all__ if not i.startswith("_")]
del old_dir
