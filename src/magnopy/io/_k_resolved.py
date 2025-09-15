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


from typing import Iterable
import warnings

try:
    import matplotlib.pyplot as plt

    MATPLOTLIB_AVAILABLE = True
    MATPLOTLIB_ERROR_MESSAGE = "If you see this message, please contact developers of the code (see magnopy.org)."
except ImportError:
    MATPLOTLIB_AVAILABLE = False
    MATPLOTLIB_ERROR_MESSAGE = "\n".join(
        [
            "Installation of matplotlib is not found, can not produce .png files",
            "Please install matplotlib with",
            "",
            "    pip install matplotlib",
            "",
        ]
    )
import numpy as np

# Save local scope at this moment
old_dir = set(dir())
old_dir.add("old_dir")


def output_k_resolved(
    data,
    data_headers=None,
    output_filename=None,
    kpoints=None,
    relative=False,
    rcell=None,
    flat_indices=None,
    digits=4,
    scientific_notation=True,
):
    r"""

    .. deprecated:: 0.2.0
        Always meant to be a temporary solution. This function will be removed in March of 2026

    Outputs any k-resolved data.

    Parameters
    ----------
    data : (N, M) |array-like|_
        Some k-resolved data. N (:math:`\ge 1`) is the amount of kpoints. M is the
        number of data modes/entries.
    data_headers : (N, ) list of str, optional
        Header for the data columns. ``[f"data {i+1} for i in range(len(data))]`` by
        default.
    output_filename : str, optional
        Name of the file for saving the data. If ``None``, then return a list of lines.
    kpoints : (N, 3) |array-like|_, optional
        List of the kpoints. ``kpoints[i]`` correspond to the ``data[i]``.
    relative : bool, default False
        If ``relative == True``, then given ``kpoints`` are interpreted as relative to
        the reciprocal unit cell. Otherwise they are interpreted as absolute coordinates.
    rcell : (3, 3) |array-like|_, optional
        Reciprocal unit cell. Rows are interpreted as vectors, columns as Cartesian
        coordinates.
    flat_indices : (N, 3) |array-like|_, optional
        kpoints converted to the list of floats according to some rule. Typically used
        for plotting the band structure. If ``None`` provided, then computed automatically,
        based on the distance of the given ``kpoints`` (whether absolute or relative).
    digits : int, default 4
        Number of digits after the comma for the data elements.
    scientific_notation : bool, default True
        Whether to use scientific notation for the data elements.

    Returns
    -------
    lines : str
        Only returned if ``output_filename is None``. Use ``print("\n".join(lines))``
        to output the results to the standard output stream.
    """

    warnings.warn(
        "This function was deprecated in the release v0.2.0. This function will be removed from magnopy in March of 2026"
    )

    # Prepare format for the data elements
    chars = digits + 1 + 1 + 3
    if scientific_notation:
        chars += 1 + 1 + 3
    chars = max(chars, 8)

    # Format string for data elements
    if scientific_notation:
        fmt = f">{chars}.{digits}e"
    else:
        fmt = f">{chars}.{digits}f"

    # Format string for kpoints
    kp_fmt = ">12.8f"

    # Prepare the header
    header = ["#"]
    if flat_indices is not None or kpoints is not None:
        header.append(f"{'flat index':>12}")
    header = header + [f"{tmp:>{chars}}" for tmp in data_headers]

    if kpoints is not None:
        if rcell is not None:
            header.append(
                " ".join([f"{f'k_{comp}':>12}" for comp in ["b1", "b2", "b3"]])
            )
            header.append(" ".join([f"{f'k_{comp}':>12}" for comp in "xyz"]))
        elif relative:
            header.append(
                " ".join([f"{f'k_{comp}':>12}" for comp in ["b1", "b2", "b3"]])
            )
        else:
            header.append(" ".join([f"{f'k_{comp}':>12}" for comp in "xyz"]))

    header = " ".join(header)

    # Output header
    if output_filename is not None:
        f = open(output_filename, "w")
        f.write(header + "\n")
    else:
        lines = [header]

    # Output data for each kpoint
    for i in range(len(data)):
        line = [" "]

        if flat_indices is not None:
            line.append(f"{flat_indices[i]:{kp_fmt}}")
        elif kpoints is not None:
            if i == 0:
                line.append(f"{0.0:{kp_fmt}}")
            else:
                line.append(
                    f"{np.linalg.norm(np.array(kpoints[i], dtype=float) - kpoints[i - 1]):{kp_fmt}}"
                )

        if isinstance(data[i], Iterable):
            for entry in data[i]:
                line.append(f"{entry:{fmt}}")
        else:
            line.append(f"{data[i]:{fmt}}")

        if kpoints is not None:
            if rcell is not None:
                if relative:
                    k_vec = kpoints[i] @ np.array(rcell)
                    line.append(
                        " ".join([f"{kpoints[i][comp]:{kp_fmt}}" for comp in range(3)])
                    )
                    line.append(
                        " ".join([f"{k_vec[comp]:{kp_fmt}}" for comp in range(3)])
                    )
                else:
                    k_vec_rel = kpoints[i] @ np.linalg.inv(rcell)
                    line.append(
                        " ".join([f"{k_vec_rel[comp]:{kp_fmt}}" for comp in range(3)])
                    )
                    line.append(
                        " ".join([f"{kpoints[i][comp]:{kp_fmt}}" for comp in range(3)])
                    )
            else:
                line.append(
                    " ".join([f"{kpoints[i][comp]:{kp_fmt}}" for comp in range(3)])
                )

        line = " ".join(line)
        if output_filename is not None:
            f.write(line + "\n")
        else:
            lines.append(line)

    if output_filename is not None:
        f.close()
    else:
        return lines


def plot_k_resolved(data, kp=None, output_filename=None, ylabel=None):
    r"""
    Plot some k-resolved data.

    If only the ``data`` are given, then an index of the omegas is used for abscissa (x
    axis).

    Parameters
    ----------
    data : (N, M) |array-like|_
        Some k-resolved data. N (:math:`\ge 1`) is the amount of kpoints. M is the
        number of data modes/entries.
    kp : :py:class:`wulfric.Kpoints`, optional.
        Instance of the :py:class:`wulfric.Kpoints` class. It should be the same
        instance that were used in the preparation of the ``data``.
    output_filename : str, optional
        Name of the file for saving the image. If ``None``, then the graph would be
        opened in the interactive matplotlib window.
    ylabel : str, optional
        Label for the ordinate (y axis).
    """

    if not MATPLOTLIB_AVAILABLE:
        import warnings

        warnings.warn(MATPLOTLIB_ERROR_MESSAGE, RuntimeWarning)

        return

    data = np.array(data).T

    fig, ax = plt.subplots()

    if len(data.shape) == 2:
        for entry in data:
            if kp is not None:
                ax.plot(kp.flat_points(), entry, lw=1, color="#A47864")
            else:
                ax.plot(entry, lw=1, color="#A47864")
                ax.set_xlim(0, len(entry))
    else:
        if kp is not None:
            ax.plot(kp.flat_points(), data, lw=1, color="#A47864")
        else:
            ax.plot(data, lw=1, color="#A47864")
            ax.set_xlim(0, len(data))

    if kp is not None:
        ax.set_xticks(kp.ticks(), kp.labels, fontsize=13)
        ax.set_xlim(kp.ticks()[0], kp.ticks()[-1])
        ax.vlines(
            kp.ticks(),
            0,
            1,
            lw=0.5,
            color="grey",
            ls="dashed",
            zorder=0,
            transform=ax.get_xaxis_transform(),
        )

    if ylabel is not None:
        ax.set_ylabel(ylabel, fontsize=15)

    ax.hlines(
        0,
        0,
        1,
        lw=0.5,
        color="grey",
        linestyle="dashed",
        transform=ax.get_yaxis_transform(),
    )

    if output_filename is not None:
        plt.savefig(output_filename, dpi=400, bbox_inches="tight")
    else:
        plt.show()

    plt.close()


# Populate __all__ with objects defined in this file
__all__ = list(set(dir()) - old_dir)
# Remove all semi-private objects
__all__ = [i for i in __all__ if not i.startswith("_")]
del old_dir
