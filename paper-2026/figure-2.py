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

import numpy as np
import matplotlib.pyplot as plt
# import matplotlib.patheffects as pe

import magnopy
import wulfric

####################
# Style parameters #
####################
# Colors
COLOR_NUMERICAL = "#006dcc"
COLOR_ANALYTICAL = "#ffcc33"
COLOR_RATIO = "#00800b"
COLOR_PARAMETERS = "#13ec95"
COLOR_PARAMETERS_BACK = "#000000"
RED = "#ec2213"
BLUE = "#1395ec"
GREY = "#bbbbbb"
LIGHT_GREY = "#cccccc"
# Fontsize baseline
FONTSIZE = 10
# Two-column
WIDTH = 2 * (3 + 3 / 8)  # inches
# Arrow style
ARROW_STYLE_1 = dict(
    angles="xy",
    scale_units="xy",
    width=0.01,
    scale=1,
    headlength=3,
    headaxislength=2.7,
)
ARROW_STYLE_2 = dict(
    angles="xy",
    scale_units="xy",
    width=0.01,
    scale=1,
    headlength=3,
    headaxislength=2.7,
)

RATIO_LINESTYLE = dict(
    color=COLOR_NUMERICAL,
    zorder=1,
    lw=1.5,
)

ANALYTICAL_LINESTYLE = dict(
    linestyle="-",
    color=COLOR_ANALYTICAL,
    lw=3,
    zorder=1,
)
NUMERICAL_LINESTYLE = dict(
    linestyle="--",
    color=COLOR_NUMERICAL,
    lw=1.5,
    zorder=2,
    dash_capstyle="round",
    dashes=[3.5, 2.5],
)


###############################
# Spin Hamiltonian parameters #
###############################
# Isotropic exchange parameters
J_FERRO = -2  # meV
J_ANTIFERRO = -J_FERRO  # meV
# MAgnetic flux density
B = 5  # T
# Easy-axis anisotropy
A = 0.5  # meV
# Antysymetric Dzyaloshinskii-Moriya interaction
D = 0.5  # meV
# Spin value
S = 3 / 2
# g-factor and Bohr magneton
G_FACTOR = 2
MU_B = 0.057883817982  # meV/T

# Reciprocal-space path
K_PATH = "GAMMA-X-M-GAMMA"


# Plots a hord of a circle
# dmi = [x,y,scale]
def hord(
    ax, offset, p1, p2, dmi=None, cap_angle=0.1, color_parameters=COLOR_PARAMETERS
):
    p1 = np.array(p1)
    p2 = np.array(p2)

    x, y = p1 - p2

    if y == 0:
        if x > 0:
            n = [0, -1]
        else:
            n = [0, 1]
    else:
        if y > 0:
            nx = -np.sqrt(y**2 / (x**2 + y**2))
        else:
            nx = np.sqrt(y**2 / (x**2 + y**2))
        n = [-nx, x / y * nx]

    n = np.array(n) * offset

    shift = n + (p1 + p2) / 2

    p1 -= shift
    p2 -= shift

    cosangle = p1 @ p2 / (np.linalg.norm(p1) * np.linalg.norm(p2))
    max_angle = np.arccos(cosangle)

    points = []
    rot_matrix = np.array(
        [
            [np.cos(max_angle), -np.sin(max_angle)],
            [np.sin(max_angle), np.cos(max_angle)],
        ]
    )
    if np.allclose(rot_matrix @ p2, p1):
        start = p2
    else:
        start = p1
    cap_angle = cap_angle * max_angle
    for angle in np.linspace(cap_angle, max_angle - cap_angle, 100):
        rot_matrix = np.array(
            [
                [np.cos(angle), -np.sin(angle)],
                [np.sin(angle), np.cos(angle)],
            ]
        )
        p = rot_matrix @ start
        points.append(p + shift)

    points = np.array(points).T

    ax.plot(
        points[0],
        points[1],
        color=color_parameters,
        lw=1.5,
        zorder=1,
        solid_capstyle="round",
    )

    if dmi is not None:
        dx = dmi[0] * dmi[2]
        dy = dmi[1] * dmi[2]
        dmi_pos = max_angle / 2.5
        rot_matrix = np.array(
            [
                [np.cos(dmi_pos), -np.sin(dmi_pos)],
                [np.sin(dmi_pos), np.cos(dmi_pos)],
            ]
        )
        x0, y0 = rot_matrix @ start + shift
        ax.quiver(
            x0 - dx * 0.4,
            y0 - dy * 0.4,
            dx,
            dy,
            angles="xy",
            scale_units="xy",
            width=0.012,
            scale=1,
            headlength=3,
            headaxislength=2.7,
            color=color_parameters,
            zorder=1,
        )


# Ferromagnet
def schema_1(ax):
    # Configure axes
    ax.set_xlim(-1.2, 2.2)
    ax.set_ylim(-1.2, 2.2)
    ax.set_xticks([])
    ax.set_yticks([])

    # Draw lattice
    lattice_style = dict(color=LIGHT_GREY, lw=1, zorder=0)
    ax.vlines([-1, 0, 1, 2], 0, 1, transform=ax.get_xaxis_transform(), **lattice_style)
    ax.hlines([-1, 0, 1, 2], 0, 1, transform=ax.get_yaxis_transform(), **lattice_style)

    # Draw direction vectors
    points = []
    for i in range(-1, 2):
        for j in range(-1, 2):
            points.append([0.5 + i, 0.5 + j])
    points = np.array(points)
    shift = 0.02
    for x, y in points:
        ax.quiver(
            x - 0.1 + shift,
            y + 0.1 - shift,
            +0.2,
            -0.2,
            color="black",
            zorder=2,
            **ARROW_STYLE_1,
        )

    # Draw interaction parameters
    offset = 0.7
    dmi_scale = 0.4
    cap_angle = 0.17
    hord(
        ax=ax,
        offset=offset,
        p1=[0.5, 0.5],
        p2=[1.5, 0.5],
        dmi=[0, 1, dmi_scale],
        cap_angle=cap_angle,
    )
    hord(
        ax=ax,
        offset=offset,
        p1=[0.5, 0.5],
        p2=[0.5, 1.5],
        dmi=[-1, 0, dmi_scale],
        cap_angle=cap_angle,
    )
    hord(
        ax=ax,
        offset=offset,
        p1=[0.5, 0.5],
        p2=[-0.5, 0.5],
        dmi=[0, -1, dmi_scale],
        cap_angle=cap_angle,
    )
    hord(
        ax=ax,
        offset=offset,
        p1=[0.5, 0.5],
        p2=[0.5, -0.5],
        dmi=[1, 0, dmi_scale],
        cap_angle=cap_angle,
    )


# Antiferromagnet
def schema_2(ax):
    # Configure axes
    ax.set_xlim(-1.7, 2.7)
    ax.set_ylim(-1.7, 2.7)
    ax.set_xticks([])
    ax.set_yticks([])

    # Draw crystal lattice
    lattice_style = dict(color=LIGHT_GREY, lw=1, zorder=0)
    ax.vlines([-1, 0, 1, 2], 0, 1, transform=ax.get_xaxis_transform(), **lattice_style)
    ax.hlines([-1, 0, 1, 2], 0, 1, transform=ax.get_yaxis_transform(), **lattice_style)

    # Draw magnetic lattice
    lattice_style = dict(
        color=LIGHT_GREY,
        lw=1,
        linestyle="dashed",
        dash_capstyle="round",
        dashes=[4, 3],
        zorder=0,
    )
    ax.plot([-2.5, 0.5], [0.5, -2.5], **lattice_style)
    ax.plot([-2.5, 2.5], [2.5, -2.5], **lattice_style)
    ax.plot([-1.5, 3.5], [3.5, -1.5], **lattice_style)
    ax.plot([0.5, 3.5], [3.5, 0.5], **lattice_style)

    ax.plot([-2.5, 0.5], [0.5, 3.5], **lattice_style)
    ax.plot([-2.5, 2.5], [-1.5, 3.5], **lattice_style)
    ax.plot([-1.5, 3.5], [-2.5, 2.5], **lattice_style)
    ax.plot([0.5, 3.5], [-2.5, 0.5], **lattice_style)

    # Draw direction vectors
    points_up = []
    points_down = []
    for i in range(-2, 3):
        for j in range(-2, 3):
            if (i + j) % 2 == 0:
                points_down.append([i + 0.5, j + 0.5])
            else:
                points_up.append([i + 0.5, j + 0.5])
    points_up = np.array(points_up)
    points_down = np.array(points_down)
    shift = 0.0
    for x, y in points_up:
        ax.quiver(
            x + 0.1 - shift,
            y - 0.1 + shift,
            -0.2,
            0.2,
            color=BLUE,
            zorder=2,
            **ARROW_STYLE_2,
        )
    for x, y in points_down:
        ax.quiver(
            x - 0.1 + shift,
            y + 0.1 - shift,
            0.2,
            -0.2,
            color=RED,
            zorder=2,
            **ARROW_STYLE_2,
        )

    # Draw interaction parameters
    offset = 0.7
    cap_angle = 0.17
    # First set
    hord(ax=ax, offset=offset, p1=[0.5, 0.5], p2=[1.5, 0.5], cap_angle=cap_angle)
    hord(ax=ax, offset=offset, p1=[0.5, 0.5], p2=[0.5, 1.5], cap_angle=cap_angle)
    hord(ax=ax, offset=offset, p1=[0.5, 0.5], p2=[-0.5, 0.5], cap_angle=cap_angle)
    hord(ax=ax, offset=offset, p1=[0.5, 0.5], p2=[0.5, -0.5], cap_angle=cap_angle)
    # Second set
    hord(ax=ax, offset=offset, p1=[0.5, -0.5], p2=[1.5, -0.5], cap_angle=cap_angle)
    hord(ax=ax, offset=offset, p1=[0.5, -0.5], p2=[0.5, 0.5], cap_angle=cap_angle)
    hord(ax=ax, offset=offset, p1=[0.5, -0.5], p2=[-0.5, -0.5], cap_angle=cap_angle)
    hord(ax=ax, offset=offset, p1=[0.5, -0.5], p2=[0.5, -1.5], cap_angle=cap_angle)


def configure_axes(ax, ax_ratio, kp):
    # Main axis
    ax.set_ylim(0, None)
    ax.tick_params(axis="y", which="major", labelsize=FONTSIZE)
    ax.set_ylabel("Energy (meV)", fontsize=FONTSIZE)
    ax.set_xticks(kp.ticks(), kp.labels, fontsize=FONTSIZE)
    ax.set_xlim(kp.ticks()[0], kp.ticks()[-1])
    ax.vlines(
        kp.ticks(),
        0,
        1,
        lw=0.5,
        transform=ax.get_xaxis_transform(),
        color=GREY,
        zorder=0,
    )

    # Ratio axis
    ax_ratio.set_xticks([])
    ax_ratio.set_yticks([0, 1, 2])
    ax_ratio.yaxis.tick_right()
    ax_ratio.yaxis.set_label_position("right")
    ax_ratio.set_ylim(0, 2)
    ax_ratio.tick_params(axis="y", which="major", labelsize=FONTSIZE)
    ax_ratio.set_ylabel("Ratio", fontsize=FONTSIZE)
    ax_ratio.set_xlim(kp.ticks()[0], kp.ticks()[-1])
    ax_ratio.vlines(
        kp.ticks(),
        0,
        1,
        lw=0.5,
        transform=ax_ratio.get_xaxis_transform(),
        color=GREY,
        zorder=0,
    )


# Ferromagnet
def plot_model_1(ax, ax_ratio):
    ################
    # Get the data #
    ################

    # Prepare a spin Hamiltonian
    cell = np.eye(3)
    atoms = dict(
        names=["Fe"], positions=[[0.5, 0.5, 0.0]], spins=[S], g_factors=[G_FACTOR]
    )
    convention = magnopy.Convention(
        multiple_counting=True, spin_normalized=False, c22=1, c1=1, c21=1
    )
    spinham = magnopy.SpinHamiltonian(cell=cell, atoms=atoms, convention=convention)

    # Add parameters to the Hamiltonian
    spinham.add_22(
        alpha=0,
        beta=0,
        nu=(1, 0, 0),
        parameter=magnopy.converter22.from_dmi(dmi=[0, D, 0])
        + magnopy.converter22.from_iso(J_FERRO),
    )
    spinham.add_22(
        alpha=0,
        beta=0,
        nu=(0, 1, 0),
        parameter=magnopy.converter22.from_dmi(dmi=[-D, 0, 0])
        + magnopy.converter22.from_iso(J_FERRO),
    )
    spinham.add_21(
        alpha=0, parameter=-A / 2 * np.array([[1, -1, 0], [-1, 1, 0], [0, 0, 0]])
    )

    # Add the Zeeman term
    spinham.add_magnetic_field(B=B / np.sqrt(2) * np.array([-1, 1, 0]))

    # Define the vacuum state
    spin_directions = np.array([[1, -1, 0]]) / np.sqrt(2)

    # Prepare LSWT
    lswt = magnopy.LSWT(spinham, spin_directions=spin_directions)

    # Configure the reciprocal-space
    kp = wulfric.Kpoints.from_crystal(cell=cell, atoms=atoms)
    kp.path = K_PATH

    # Compute dispersions
    omega_analytical = []
    omega_numerical = []
    for k in kp.points(relative=False):
        omega_analytical.append(
            G_FACTOR * MU_B * B
            + 2 * S * A
            - 4 * S * J_FERRO * (2 - np.cos(k @ cell[0]) - np.cos(k @ cell[1]))
            - 2 * np.sqrt(2) * D * S * (np.sin(k @ cell[0]) + np.sin(k @ cell[1]))
        )
        omega_numerical.append(lswt.omega(k=k, relative=False)[0].real)
    omega_analytical = np.array(omega_analytical)
    omega_numerical = np.array(omega_numerical)

    #################
    # Plot the data #
    #################

    # Analytical
    ax.plot(
        kp.flat_points(), omega_analytical, label="Analytical", **ANALYTICAL_LINESTYLE
    )
    # Numerical
    ax.plot(kp.flat_points(), omega_numerical, label="Magnopy", **NUMERICAL_LINESTYLE)
    # Ratio
    ax_ratio.plot(
        kp.flat_points(),
        omega_numerical / omega_analytical,
        **RATIO_LINESTYLE,
    )

    # Display the legend
    ax.legend(framealpha=1, fontsize=FONTSIZE, edgecolor="inherit", loc="lower center")

    # Apply the axes configuration
    configure_axes(ax, ax_ratio, kp)


# Antiferromagnet
def plot_model_2(ax, ax_ratio):
    ################
    # Get the data #
    ################

    # Prepare a spin Hamiltonian
    cell = np.array([[1, 1, 0], [-1, 1, 0], [0, 0, 1]], dtype=float)
    atoms = dict(
        names=["Fe1", "Fe2"],
        positions=[[0, 0, 0], [0.5, 0.5, 0.0]],
        spins=[S, S],
        g_factors=[G_FACTOR, G_FACTOR],
    )
    convention = magnopy.Convention(
        multiple_counting=True, spin_normalized=False, c22=1, c1=1, c21=1
    )
    spinham = magnopy.SpinHamiltonian(cell=cell, atoms=atoms, convention=convention)

    # Add parameters to the Hamiltonian
    spinham.add_22(
        alpha=1,
        beta=0,
        nu=(0, 0, 0),
        parameter=magnopy.converter22.from_iso(J_ANTIFERRO),
    )
    spinham.add_22(
        alpha=1,
        beta=0,
        nu=(1, 0, 0),
        parameter=magnopy.converter22.from_iso(J_ANTIFERRO),
    )
    spinham.add_22(
        alpha=1,
        beta=0,
        nu=(0, 1, 0),
        parameter=magnopy.converter22.from_iso(J_ANTIFERRO),
    )
    spinham.add_22(
        alpha=1,
        beta=0,
        nu=(1, 1, 0),
        parameter=magnopy.converter22.from_iso(J_ANTIFERRO),
    )

    spinham.add_21(
        alpha=0, parameter=-A / 2 * np.array([[1, -1, 0], [-1, 1, 0], [0, 0, 0]])
    )
    spinham.add_21(
        alpha=1, parameter=-A / 2 * np.array([[1, -1, 0], [-1, 1, 0], [0, 0, 0]])
    )

    # Add the Zeeman term
    spinham.add_magnetic_field(B=B / np.sqrt(2) * np.array([-1, 1, 0]))

    # Define the vacuum state
    spin_directions = np.array([[-1, 1, 0], [1, -1, 0]]) / np.sqrt(2)

    # Prepare LSWT
    lswt = magnopy.LSWT(spinham, spin_directions=spin_directions)

    # Configure the reciprocal-space
    kp = wulfric.Kpoints.from_crystal(cell=cell, atoms=atoms)
    kp.path = K_PATH

    # Compute dispersions
    omega_analytical = []
    omega_numerical = []
    for k in kp.points(relative=False):
        omega_analytical.append(
            [
                G_FACTOR * MU_B * B
                + 8
                * J_ANTIFERRO
                * S
                * np.sqrt(
                    (1 + A / 4 / J_ANTIFERRO) ** 2
                    - np.cos(k @ cell[0] / 2) ** 2 * np.cos(k @ cell[1] / 2) ** 2
                ),
                -G_FACTOR * MU_B * B
                + 8
                * J_ANTIFERRO
                * S
                * np.sqrt(
                    (1 + A / 4 / J_ANTIFERRO) ** 2
                    - np.cos(k @ cell[0] / 2) ** 2 * np.cos(k @ cell[1] / 2) ** 2
                ),
            ]
        )
        omega_numerical.append(lswt.omega(k=k, relative=False))
    omega_analytical = np.array(omega_analytical).T
    omega_numerical = np.array(omega_numerical).real.T

    #################
    # Plot the data #
    #################

    for mode_index in range(2):
        # Analytical
        ax.plot(
            kp.flat_points(),
            omega_analytical[mode_index],
            label="Analytical" if mode_index == 0 else None,
            **ANALYTICAL_LINESTYLE,
        )
        # Numerical
        ax.plot(
            kp.flat_points(),
            omega_numerical[mode_index],
            label="Magnopy" if mode_index == 0 else None,
            **NUMERICAL_LINESTYLE,
        )
        # Ratio
        ax_ratio.plot(
            kp.flat_points(),
            omega_numerical[mode_index] / omega_analytical[mode_index],
            **RATIO_LINESTYLE,
        )

    # Display the legend
    ax.legend(framealpha=1, fontsize=FONTSIZE, edgecolor="inherit", loc="lower left")

    # Apply the axes configuration
    configure_axes(ax, ax_ratio, kp)


def main():
    ##################
    # Prepare canvas #
    ##################

    # Layout parameters
    sqz = 0.85
    dx = 0.08
    dy = (sqz - 1 + 4 * dx) / 3
    zoom = 1

    ratio = 1 / 5

    # Figure
    fig = plt.figure(figsize=(WIDTH, WIDTH * sqz))

    # Global reference frame
    frameax_size = 0.08
    frameax_x = 0.5
    frameax_y = 0.6
    frameax = fig.add_axes(
        [
            frameax_x - frameax_size / 2,
            frameax_y - frameax_size / 2 / sqz,
            frameax_size,
            frameax_size / sqz,
        ]
    )
    frame_arrow_style = dict(
        angles="xy",
        scale_units="xy",
        width=0.05,
        scale=1,
        headlength=3,
        headaxislength=2.7,
        capstyle="round",
    )
    frameax.quiver(0.1, 0.1, 0.8, 0, **frame_arrow_style, color="tab:red")
    frameax.text(
        0.9, 0.2, "x", fontsize=FONTSIZE, ha="center", va="bottom", color="tab:red"
    )
    frameax.quiver(0.1, 0.1, 0, 0.8, **frame_arrow_style, color="tab:green")
    frameax.text(
        0.2, 0.9, "y", fontsize=FONTSIZE, ha="left", va="center", color="tab:green"
    )

    frameax.set_xlim(0, 1)
    frameax.set_ylim(0, 1)
    frameax.axis("off")

    # Top left axis
    ax_a = fig.add_axes(
        [
            dx - (zoom - 1) * (1 - 4 * dx) / 4,
            (1 + dy) / 2 - (zoom - 1) * (1 - 3 * dy) / 4,
            zoom * (1 - 4 * dx) / 2,
            zoom * (1 - 3 * dy) / 2,
        ]
    )
    # Top right axis
    ax_b = fig.add_axes(
        [
            (1 + 2 * dx) / 2 - (zoom - 1) * (1 - 4 * dx) / 4,
            (1 + dy) / 2 - (zoom - 1) * (1 - 3 * dy) / 4,
            zoom * (1 - 4 * dx) / 2,
            zoom * (1 - 3 * dy) / 2,
        ]
    )
    # Bottom left axis
    ax_c = fig.add_axes(
        [
            dx,
            dy,
            (1 - 4 * dx) / 2,
            (1 - ratio) * (1 - 3 * dy) / 2,
        ]
    )
    # Bottom right axis
    ax_d = fig.add_axes(
        [
            (1 + 2 * dx) / 2,
            dy,
            (1 - 4 * dx) / 2,
            (1 - ratio) * (1 - 3 * dy) / 2,
        ]
    )
    # Ratio for bottom left axis
    ax_c_r = fig.add_axes(
        [
            dx,
            dy + (1 - ratio) * (1 - 3 * dy) / 2,
            (1 - 4 * dx) / 2,
            ratio * (1 - 3 * dy) / 2,
        ]
    )
    # Ratio for bottom right axis
    ax_d_r = fig.add_axes(
        [
            (1 + 2 * dx) / 2,
            dy + (1 - ratio) * (1 - 3 * dy) / 2,
            (1 - 4 * dx) / 2,
            ratio * (1 - 3 * dy) / 2,
        ]
    )

    # subfigure labels
    style = dict(fontsize=1.2 * FONTSIZE, ha="right", va="top")
    fig.text(dx / 2, 1 - dy + (zoom - 1) * (1 - 4 * dx) / 4, "(a)", **style)
    fig.text((1 + dx) / 2, 1 - dy + (zoom - 1) * (1 - 4 * dx) / 4, "(b)", **style)
    fig.text(dx / 2, (1 - dy) / 2, "(c)", **style)
    fig.text((1 + dx) / 2, (1 - dy) / 2, "(d)", **style)

    ####################################
    # Plot schematic of the two models #
    ####################################

    # Ferromagnet
    schema_1(ax_a)

    # Antiferromagnet
    schema_2(ax_b)
    ###################
    # Plot comparison #
    ###################

    # Ferromagnet
    plot_model_1(ax=ax_c, ax_ratio=ax_c_r)

    # Antiferromagnet
    plot_model_2(ax=ax_d, ax_ratio=ax_d_r)

    ###############
    # Save figure #
    ###############
    plt.savefig("figure-2.eps", dpi=600)
    plt.savefig("figure-2.png", dpi=600)
    plt.savefig("figure-2.pdf", dpi=600)
    plt.close()


if __name__ == "__main__":
    main()
