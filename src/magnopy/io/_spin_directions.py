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


import numpy as np

try:
    import plotly.graph_objects as go

    PLOTLY_AVAILABLE = True
    PLOTLY_ERROR_MESSAGE = "If you see this message, please contact developers of the code (see magnopy.org)."
except ImportError:
    PLOTLY_AVAILABLE = False
    PLOTLY_ERROR_MESSAGE = "\n".join(
        [
            "Installation of plotly is not found, can not produce .html pictures.",
            "Either install plotly with",
            "",
            "    pip install plotly",
            "",
            "Or disable html output with",
            "",
            "    no_html=True (when using magnopy as Python library)",
            "    --no-html (when using magnopy through command line interface)",
            "",
        ]
    )


# Save local scope at this moment
old_dir = set(dir())
old_dir.add("old_dir")


def read_spin_directions(filename: str):
    r"""
    Read directions of the spins from the file.

    Parameters
    ----------
    filename : str or (3*M, ) |array-like|_
        File with the spin directions. See notes for the specification of the file
        format.

    Returns
    -------
    spin_directions : (M, ) :numpy:`ndarray`
        If ``spin_directions`` is an |array-like|_, then first three elements are

    Notes
    -----
    The file is expected to contain three numbers per line, here is an example for two
    spins

    .. code-block:: text

        S1_x S1_y S1_z
        S2_x S2_y S2_z

    Only the direction of the spin vector is recognized, the modulus is ignored.
    Comments are allowed at any place of the file and preceded by the symbol "#".
    If the symbol "#" is found, the part of the line after it is ignored. Here
    are examples of valid use of the comments

    .. code-block:: text

        # Spin vectors for the material XX
        S1_x S1_y S1_z # Atom X1
        # This comments is here by some reason
        S2_x S2_y S2_z # Atom X2

    """

    spin_directions = []
    with open(filename, "r") as f:
        for i, line in enumerate(f):
            # Remove comment lines
            if line.startswith("#"):
                continue
            # Remove inline comments and leading/trailing whitespaces
            line = line.split("#")[0].strip()
            # Check for empty lines empty lines
            if line:
                line = line.split()
                if len(line) != 3:
                    raise ValueError(
                        f"Expected three numbers per line (in line{i}),"
                        f"got: {len(line)}."
                    )
                for tmp in line:
                    spin_directions.append(float(tmp))

    if len(spin_directions) % 3 != 0:
        raise ValueError(
            f"Length of the spin list should be dividable by three, got: {len(spin_directions)}."
        )

    spin_directions = np.array(spin_directions, dtype=float)
    # Pay attention to the np.reshape keywords
    spin_directions = np.reshape(spin_directions, (len(spin_directions) // 3, 3))
    spin_directions = (
        spin_directions / np.linalg.norm(spin_directions, axis=1)[:, np.newaxis]
    )
    return spin_directions


def _plot_cones(fig, positions, spin_directions, color, name=None):
    scale = 0.5

    if not PLOTLY_AVAILABLE:
        raise ImportError(PLOTLY_ERROR_MESSAGE)

    # Prepare data
    x, y, z = np.transpose(positions, axes=(1, 0))
    u, v, w = np.transpose(spin_directions, axes=(1, 0))

    fig.add_traces(
        data=go.Cone(
            x=x + u * scale,
            y=y + v * scale,
            z=z + w * scale,
            u=u * (1 - scale),
            v=v * (1 - scale),
            w=w * (1 - scale),
            sizemode="raw",
            anchor="tail",
            legendgroup=name,
            name=name,
            showlegend=name is not None,
            showscale=False,
            colorscale=[color, color],
            hoverinfo="none",
        )
    )

    fig.add_traces(
        data=go.Scatter3d(
            mode="markers",
            x=x,
            y=y,
            z=z,
            marker=dict(size=10, color=color),
            hoverinfo="none",
            showlegend=False,
            legendgroup=name,
        )
    )

    for i in range(0, len(x)):
        fig.add_traces(
            dict(
                x=[x[i], x[i] + u[i] * scale],
                y=[y[i], y[i] + v[i] * scale],
                z=[z[i], z[i] + w[i] * scale],
                mode="lines",
                type="scatter3d",
                hoverinfo="none",
                line={"color": color, "width": 10},
                legendgroup=name,
                showlegend=False,
            )
        )

    return fig


# Populate __all__ with objects defined in this file
__all__ = list(set(dir()) - old_dir)
# Remove all semi-private objects
__all__ = [i for i in __all__ if not i.startswith("_")]
del old_dir
