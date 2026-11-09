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


# Single-column
WIDTH = 3 + 3 / 8  # inches
FONTSIZE = 10


def main():
    # Prepare canvas
    fig = plt.figure(figsize=(WIDTH, WIDTH * 0.75))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    # Define styles
    text_box = dict(
        ha="center",
        va="center",
        fontsize=FONTSIZE,
        bbox=dict(boxstyle="round,pad=0.5", fc="white", ec="black", lw=1.5),
    )
    text_arrow = dict(ha="center", va="center", fontsize=FONTSIZE)
    arrow_style = dict(
        angles="xy", scale_units="xy", scale=1, headlength=3, headaxislength=2.7
    )

    # Add text in boxes
    ax.text(0.5, 0.8, "Experimental data", **text_box)
    ax.text(0.5, 0.5, "Spin Hamiltonian", **text_box)
    ax.text(0.5, 0.2, "ab initio methods", **text_box)

    # Add text for arrows
    ax.text(0.65, 0.37, "(2) Calculate", **text_arrow)
    ax.text(0.8, 0.67, "(3) Simulate", **text_arrow)
    ax.text(0.25, 0.63, "(1) Fit", **text_arrow)

    # Draw arrows
    ax.quiver(0.35, 0.73, 0, -0.16, **arrow_style)
    ax.quiver(0.65, 0.57, 0, +0.16, **arrow_style)
    ax.quiver(0.5, 0.27, 0, +0.16, **arrow_style)

    # Save figure
    plt.savefig("figure-1.png", dpi=400)
    plt.savefig("figure-1.eps", dpi=400)
    plt.savefig("figure-1.pdf", dpi=400)
    plt.close()


if __name__ == "__main__":
    main()
