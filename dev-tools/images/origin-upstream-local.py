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
from argparse import ArgumentParser
from math import atan, pi

import matplotlib.pyplot as plt
import numpy as np


def arrow_with_text(
    ax, x, y, u, v, text=None, n_shift=2.1, quiver_kwargs=None, text_kwargs=None
):
    if quiver_kwargs is None:
        quiver_kwargs = {}
    if text_kwargs is None:
        text_kwargs = {}

    center = np.array((x + u / 2, y + v / 2), dtype=float)

    if u != 0:
        normal = (-v / u, 1)
        angle = atan(v / u) / pi * 180
    else:
        normal = (1, -u / v)
        if v > 0:
            angle = 90
        else:
            angle = 270

    normal = np.array(normal, dtype=float)
    normal /= np.linalg.norm(normal)

    if text is not None:
        ax.quiver(x, y, u, v, **quiver_kwargs, zorder=6)

    ax.text(
        *(center + n_shift * normal),
        text,
        **text_kwargs,
        rotation=angle,
        rotation_mode="anchor",
        zorder=5,
    )


def main(root_directory):
    fig, ax = plt.subplots(figsize=(5, 5))
    bbox = dict(
        boxstyle="round",
        ec=(35 / 255, 106 / 255, 211 / 255, 1),
        fc=(35 / 255, 106 / 255, 211 / 255, 0.5),
    )
    bbox_arrow = dict(
        boxstyle="round",
        ec=(1, 1, 1, 1),
        fc=(1, 1, 1, 1),
    )
    text_style = dict(bbox=bbox, size=15, va="center", ha="center")
    text_arrow = dict(bbox=bbox_arrow, size=13, va="center", ha="center")
    arrow_style = dict(
        angles="xy", scale_units="xy", scale=1, headlength=3, headaxislength=2.7
    )

    ax.text(25, 10, "Local repository", **text_style)

    ax.text(13, 45, "Origin", **text_style)
    ax.text(37, 45, "Upstream", **text_style)

    arrow_with_text(
        ax,
        8,
        15,
        0,
        10,
        n_shift=-2.1,
        text="git push",
        quiver_kwargs=arrow_style,
        text_kwargs=text_arrow,
    )
    arrow_with_text(
        ax,
        18,
        25,
        0,
        -10,
        text="git pull",
        quiver_kwargs=arrow_style,
        text_kwargs=text_arrow,
    )
    arrow_with_text(
        ax,
        40,
        30,
        0,
        -20,
        text="git pull upstream",
        quiver_kwargs=arrow_style,
        text_kwargs=text_arrow,
    )
    arrow_with_text(
        ax,
        20,
        35,
        10,
        0,
        text="pull request",
        quiver_kwargs=arrow_style,
        text_kwargs=text_arrow,
    )
    arrow_with_text(
        ax,
        30,
        30,
        -10,
        0,
        n_shift=-2.1,
        text="sync fork",
        quiver_kwargs=arrow_style,
        text_kwargs=text_arrow,
    )

    ax.text(25, 5, "git commit", **text_arrow)

    ax.set_xlim(0, 50)
    ax.set_ylim(0, 50)

    ax.hlines([1, 20, 49], 1, 49, color="black", lw=1)
    ax.vlines(25, 20, 49, color="black", lw=1)
    ax.vlines([1, 49], 1, 49, color="black", lw=1)

    # ax.vlines([i for i in range(0, 51)], 0, 50, color="grey", lw=0.1)
    # ax.hlines([i for i in range(0, 51)], 0, 50, color="grey", lw=0.1)
    # ax.vlines([10 * i for i in range(0, 5)], 0, 50, ls="--", color="black", lw=0.5)
    # ax.hlines([10 * i for i in range(0, 5)], 0, 50, ls="--", color="black", lw=0.5)

    ax.axis("off")
    filename = os.path.join(
        root_directory,
        "docs",
        "images",
        "origin-upstream-local.png",
    )
    plt.savefig(filename, dpi=600, bbox_inches="tight")
    print(f"File is saved in {os.path.abspath(filename)}")


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument(
        "-rd",
        "--root_directory",
        type=str,
        help="Root directory of magnopy",
        required=True,
    )

    main(**vars(parser.parse_args()))
