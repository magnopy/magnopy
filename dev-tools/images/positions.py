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

import matplotlib.pyplot as plt
import numpy as np


def get_ax_size(ax, fig):
    bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    width, height = bbox.width, bbox.height
    width *= fig.dpi
    height *= fig.dpi
    return width, height


def arrow_with_text(
    ax, origin, target, text, shift=(0, 0), color="black", arrow_style=None
):
    shift = np.array(shift)
    origin = np.array(origin)
    target = np.array(target)
    if arrow_style is None:
        arrow_style = dict(
            angles="xy",
            scale_units="xy",
            scale=1,
            width=0.003,
            headlength=8,
            headwidth=3,
            headaxislength=7,
        )
    ax.quiver(*origin, *(target - origin), **arrow_style, color=color)
    ax.text(
        *((origin + target) / 2 + shift),
        text,
        ha="center",
        va="center",
        fontsize=20,
        color=color,
    )


def plot_origin_point(ax, x, y):
    ax.scatter(x, y, color="black", s=40, zorder=10)
    ax.text(
        x,
        y,
        "\nReference point",
        color="black",
        zorder=1,
        ha="center",
        va="top",
        fontsize=20,
    )


def prepare_fig_ax(sites=False):
    fig, ax = plt.subplots(figsize=(18, 5))

    # Cell's borders

    ax_width, _ = get_ax_size(ax=ax, fig=fig)
    linewidth = 2
    lw_in_pixels = linewidth * fig.dpi / 72.0
    lw_relative = lw_in_pixels / ax_width

    xlim = (-0.1, 18.1)
    offset = lw_relative * (xlim[1] - xlim[0])
    ax.plot(
        [offset, 6 - offset, 6 - offset, offset, offset],
        [0, 0, 5, 5, 0],
        lw=linewidth,
        color="tab:red",
        zorder=0,
    )
    ax.plot(
        [6 + offset, 12 - offset, 12 - offset, 6 + offset, 6 + offset],
        [0, 0, 5, 5, 0],
        lw=linewidth,
        color="tab:green",
        zorder=0,
    )
    ax.plot(
        [12 + offset, 18 - offset, 18 - offset, 12 + offset, 12 + offset],
        [0, 0, 5, 5, 0],
        lw=linewidth,
        color="tab:blue",
        zorder=0,
    )

    ax.text(
        3,
        5.5,
        R"Cell index: $\mu-\nu$",
        ha="center",
        va="bottom",
        fontsize=20,
        color="tab:red",
    )
    ax.text(
        9,
        5.5,
        R"Cell index: $\mu$",
        ha="center",
        va="bottom",
        fontsize=20,
        color="tab:green",
    )
    ax.text(
        15,
        5.5,
        R"Cell index: $\mu+\nu$",
        ha="center",
        va="bottom",
        fontsize=20,
        color="tab:blue",
    )

    # Reference frame
    origin = np.array((6, -2))
    plot_origin_point(ax, *origin)

    if sites:
        # Magnetic sites
        site_1 = [1.5, 4]
        site_2 = [4.5, 2]
        sites = []

        for i in range(3):
            for j in range(1):
                sites.append(np.array((site_1[0] + i * 6, site_1[1] + j * 5)))

                ax.scatter(
                    sites[-1][0],
                    sites[-1][1],
                    s=400,
                    color="tab:orange",
                    zorder=4,
                )

                sites.append(np.array((site_2[0] + i * 6, site_2[1] + j * 5)))

                ax.scatter(
                    sites[-1][0],
                    sites[-1][1],
                    s=400,
                    color="tab:olive",
                    zorder=4,
                )

    return origin, sites, offset, fig, ax


def plot_unit_cells(root_directory):
    origin, _, offset, fig, ax = prepare_fig_ax()

    arrow_style = dict(
        angles="xy",
        scale_units="xy",
        scale=1,
        width=0.004,
        headlength=5,
        headaxislength=4,
        headwidth=4,
    )

    arrow_with_text(
        ax,
        origin,
        (offset, 0),
        R"$\boldsymbol{r}_{\mu - \nu}$",
        shift=(-0.35, -0.35),
        color="tab:red",
        arrow_style=arrow_style,
    )
    arrow_with_text(
        ax,
        origin,
        (6 + offset, 0),
        R"$\boldsymbol{r}_{\mu}$",
        shift=(-0.45, 0),
        color="tab:green",
        arrow_style=arrow_style,
    )
    arrow_with_text(
        ax,
        origin,
        (12 + offset, 0),
        R"$\boldsymbol{r}_{\mu+\nu}$",
        shift=(-0.45, 0.35),
        color="tab:blue",
        arrow_style=arrow_style,
    )

    # Between unit cells
    arrow_with_text(
        ax,
        (6 + offset, 0.5),
        (12 - offset, 0.5),
        R"$\boldsymbol{r}_{\nu}$",
        shift=(0, 0.5),
        color="black",
        arrow_style=arrow_style,
    )
    arrow_with_text(
        ax,
        (6 - offset, 1),
        (offset, 1),
        R"$\boldsymbol{r}_{-\nu}$",
        shift=(0, 0.5),
        color="black",
        arrow_style=arrow_style,
    )
    # arrow_with_text(
    #     ax,
    #     (6, 2.5),
    #     (18, 2.5),
    #     R"$\boldsymbol{r}_{2\nu}$",
    #     shift=(-0.5, 0.3),
    #     color="black",
    # )

    ax.axis("off")
    ax.set_aspect(1)
    filename = os.path.join(
        root_directory,
        "docs",
        "images",
        "positions_cells.png",
    )
    fig.savefig(filename, dpi=400, bbox_inches="tight")
    print(f"File is saved in {os.path.abspath(filename)}")
    plt.close()


def plot_sites(root_directory):
    origin, sites, offset, fig, ax = prepare_fig_ax(sites=True)

    arrow_style = dict(
        angles="xy",
        scale_units="xy",
        scale=1,
        width=0.004,
        headlength=5,
        headaxislength=4,
        headwidth=4,
        zorder=5,
    )

    scale = 0.045

    # Sites from unit cells
    arrow_with_text(
        ax,
        (offset, 0),
        (sites[0][0] * (1 - scale), sites[0][1] * (1 - scale)),
        R"$\boldsymbol{r}_{\alpha}$",
        shift=(-0.2, 0.5),
        color="tab:orange",
        arrow_style=arrow_style,
    )
    arrow_with_text(
        ax,
        (offset, 0),
        (sites[1][0] * (1 - scale), sites[1][1] * (1 - scale)),
        R"$\boldsymbol{r}_{\beta}$",
        shift=(-0.2, 0.5),
        color="tab:olive",
        arrow_style=arrow_style,
    )
    ax.text(
        sites[0][0] - 0.7,
        sites[0][1] + 0.3,
        R"Atom index: $\alpha$",
        color="tab:orange",
        va="bottom",
        ha="left",
        fontsize=20,
    )
    ax.text(
        sites[1][0] - 3,
        sites[1][1] + 0.3,
        R"Atom index: $\beta$",
        color="tab:olive",
        va="bottom",
        ha="left",
        fontsize=20,
    )

    # More vectors
    arrow_with_text(
        ax,
        origin,
        (
            sites[2][0] - 0.04 * (sites[2][0] - origin[0]),
            sites[2][1] - 0.04 * (sites[2][1] - origin[1]),
        ),
        R"$\boldsymbol{r}_{\mu,\alpha}$",
        shift=(-1.3, -1.8),
        color="black",
        arrow_style=arrow_style,
    )
    arrow_with_text(
        ax,
        origin,
        (
            sites[5][0] - 0.03 * (sites[5][0] - origin[0]),
            sites[5][1] - 0.03 * (sites[5][1] - origin[1]),
        ),
        R"$\boldsymbol{r}_{\mu+\nu,\beta}$",
        shift=(-1.5, -1.5),
        color="black",
        arrow_style=arrow_style,
    )
    arrow_with_text(
        ax,
        (
            sites[2][0] + 0.03 * (sites[5][0] - sites[2][0]),
            sites[2][1] + 0.03 * (sites[5][1] - sites[2][1]),
        ),
        (
            sites[5][0] - 0.03 * (sites[5][0] - sites[2][0]),
            sites[5][1] - 0.03 * (sites[5][1] - sites[2][1]),
        ),
        R"$\boldsymbol{r}_{\nu,\alpha\beta}$",
        shift=(-1, 0.8),
        color="black",
        arrow_style=arrow_style,
    )

    ax.axis("off")
    ax.set_aspect(1)
    filename = os.path.join(
        root_directory,
        "docs",
        "images",
        "positions_sites.png",
    )
    fig.savefig(filename, dpi=400, bbox_inches="tight")
    print(f"File is saved in {os.path.abspath(filename)}")
    plt.close()


def main(root_directory):
    plot_unit_cells(root_directory)
    plot_sites(root_directory)


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
