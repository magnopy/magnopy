# ================================== LICENSE ===================================
# Magnopy - Python package for magnons.
#
# Copyright (C) 2023 Magnopy Team
#
# e-mail: anry@uv.es, web: magnopy.org
#
# This program is free software: you  can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the  Free Software
# Foundation,  either  version 3  of the License,  or (at your option) any later
# version.
#
# This program is distributed in the  hope  that it will be useful,  but WITHOUT
# ANY WARRANTY;  without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the  GNU General Public License  along with
# this program.  If not, see <https://www.gnu.org/licenses/>.
# ================================ END LICENSE =================================

import matplotlib.pyplot as plt

# Style parameters
COLOR_BORDER = "#017acb"  # "#265DC3"
COLOR_FILL = "#e6f5ff"  # "#E6F4FD"  # "#E9EEF7"
FONTSIZE = 10


# Display the cat of magnopy's logo
def plot_logo(ax, x0, y0, size):
    sx = size
    sy = size / 5 * 5.5

    style = dict(color="#ffffff", lw=0, zorder=1)
    ax.fill_between(
        [x0 - 9.5 * sx, x0 - 5.5 * sx],
        [y0 + 8.5 * sy, y0 + 8.5 * sy],
        [y0 + 11.5 * sy, y0 + 11.5 * sy],
        **style,
    )

    style["color"] = COLOR_BORDER
    ax.fill_between(
        [
            x0 - 10 * sx,
            x0 - 9 * sx,
            x0 - 9 * sx,
            x0 - 8 * sx,
            x0 - 8 * sx,
            x0 - 7 * sx,
            x0 - 7 * sx,
            x0 - 6 * sx,
            x0 - 6 * sx,
            x0 - 5 * sx,
            x0 - 5 * sx,
            x0 - sx,
            x0 - sx,
            x0,
        ],
        [
            y0 + 8 * sy,
            y0 + 8 * sy,
            y0,
            y0,
            y0,
            y0,
            y0,
            y0,
            y0,
            y0,
            y0,
            y0,
            y0 + sy,
            y0 + sy,
        ],
        [
            y0 + 13 * sy,
            y0 + 13 * sy,
            y0 + 9 * sy,
            y0 + 9 * sy,
            y0 + 11 * sy,
            y0 + 11 * sy,
            y0 + 9 * sy,
            y0 + 9 * sy,
            y0 + 6 * sy,
            y0 + 6 * sy,
            y0 + sy,
            y0 + sy,
            y0 + 6 * sy,
            y0 + 6 * sy,
        ],
        **style,
    )

    ax.fill_between(
        [x0 - 9 * sx, x0 - 6 * sx, x0 - 6 * sx, x0 - 5 * sx],
        [y0 + 11 * sy, y0 + 11 * sy, y0 + 8 * sy, y0 + 8 * sy],
        [y0 + 12 * sy, y0 + 12 * sy, y0 + 13 * sy, y0 + 13 * sy],
        **style,
    )

    ax.fill_between(
        [x0 - 4 * sx, x0 - 3 * sx, x0 - 3 * sx, x0 - sx],
        [y0 + 4 * sy, y0 + 4 * sy, y0 + 6 * sy, y0 + 6 * sy],
        [y0 + 6 * sy, y0 + 6 * sy, y0 + 7 * sy, y0 + 7 * sy],
        **style,
    )


def main():
    ##################
    # Prepare canvas #
    ##################
    cm = 1 / 2.54
    fig = plt.figure(figsize=(5.5 * cm, 5 * cm))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    ax.axis("off")

    ####################
    # Draw a rectangle #
    ####################
    left = 5
    right = 95
    bottom = 15
    top = 85
    gap_left = 47
    gap_right = 53
    # Fill
    ax.fill_between(
        [left, right], [bottom, bottom], [top, top], color=COLOR_FILL, zorder=0, lw=0
    )
    # Border
    border = dict(color=COLOR_BORDER, lw=3, solid_capstyle="round", zorder=1)

    ax.plot([gap_left, left, left, gap_left], [top, top, bottom, bottom], **border)
    ax.plot([gap_right, right, right, gap_right], [bottom, bottom, top, top], **border)

    #################
    # Draw the logo #
    #################
    plot_logo(ax, 92, 17, 1)

    #################
    # Define styles #
    #################
    # Style of text
    text_out = dict(
        fontsize=FONTSIZE, color="black", ha="center", va="center", zorder=1
    )
    text_in = dict(
        color=COLOR_BORDER, va="center", ha="center", fontsize=FONTSIZE, zorder=1
    )
    text_in_small = dict(
        color=COLOR_BORDER, va="center", fontsize=0.75 * FONTSIZE, zorder=1
    )

    # Style of the arrows
    def arrow(x0, y0, x1, y1, color):
        ax.quiver(
            x0,
            y0,
            x1 - x0,
            y1 - y0,
            angles="xy",
            scale_units="xy",
            width=0.006,
            scale=1,
            headlength=3,
            headaxislength=2.7,
            color=color,
            zorder=1,
        )

    # Style of the lines
    def line(x0, y0, x1, y1):
        ax.plot(
            [x0, x1],
            [y0, y1],
            color=COLOR_BORDER,
            lw=1,
            zorder=1,
            solid_capstyle="round",
        )

    ########################
    # Text outside the box #
    ########################
    ax.text(50, (100 + top) / 2, R"parameters & configuration", **text_out)
    ax.text(50, bottom / 2, R"output (visual & numerical)", **text_out)

    #######################
    # Text inside the box #
    #######################
    # Spin Hamiltonian
    ax.text(50, 72.5, "Spin Hamiltonian", **text_in)
    ax.text(84, 72.5, R"$\hat{\mathcal{H}}$", **text_in_small, ha="left")
    # Vacuum state
    ax.text(50, 56, "Vacuum state", **text_in)
    ax.text(80, 56, R"$\boldsymbol{z}_{\alpha}$", **text_in_small, ha="left")
    # LSWT
    ax.text(64, 41, "LSWT", **text_in)
    ax.text(79, 44, R"$E^{(2)}$", **text_in_small, ha="left")
    ax.text(79, 38, R"$\omega_{\beta}(\boldsymbol{k})$", **text_in_small, ha="left")
    # Energy
    ax.text(35, 41, "Energy", **text_in)
    ax.text(8, 44, R"$E^{(0)}$", **text_in_small, ha="left")
    ax.text(8, 38, R"$E^{corr}$", **text_in_small, ha="left")
    # Post-processing
    ax.text(50, 25, "Post-processing", **text_in)

    #################
    # Content lines #
    #################
    # Spin Hamiltonian
    line(79, 72.5, 83, 72.5)
    # Vacuum state
    line(74, 56, 78, 56)
    # LSWT
    line(23, 41, 19, 44)
    line(23, 41, 19, 38)
    # Energy
    line(74, 41, 78, 44)
    line(74, 41, 78, 38)

    ##########
    # Arrows #
    ##########
    # In
    arrow(50, top + 4, 50, top - 4, "black")
    # Spin Hamiltonian -> Vacuum state
    arrow(50, 68.25, 50, 60.25, COLOR_BORDER)
    # Vacuum state -> LSWT
    arrow(46, 52.5, 35, 45.5, COLOR_BORDER)
    # Vacuum state -> Energy
    arrow(54, 52.5, 65, 45.5, COLOR_BORDER)
    # LSWT -> Post-processing
    arrow(35, 37, 46, 29, COLOR_BORDER)
    # Energy -> Post-processing
    arrow(65, 37, 54, 29, COLOR_BORDER)
    # Out
    arrow(50, bottom + 4, 50, bottom - 4, "black")

    ###############
    # Save figure #
    ###############
    plt.savefig("toc.eps", dpi=600)
    plt.savefig("toc.png", dpi=600)
    plt.savefig("toc.pdf", dpi=600)
    plt.close()


if __name__ == "__main__":
    main()
