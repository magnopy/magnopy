# MAGNOPY - Python package for magnons.
# Copyright (C) 2023-2025 Magnopy Team
#
# e-mail: anry@uv.es, web: magnopy.com
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

from _energy import Energy
from tqdm import tqdm
import numpy as np

class Energy_modified(Energy):

    def __init__(self, spinham, beta):
        super().__init__(spinham)
        self.beta = beta

    def generate_possible_flippings(self, number_flippings=0, thetas_flag=False, phis_flag=False):
        r"""
        Generate random sites, and random versors (pointing to specific thetas and phis in polar coordinates),
            and random angles of rotation with respect to the random versors.

        Parameters
        ----------
        number_flippings : int, default 0
            Number of successive flippings considered in the Monte Carlo simulation
        thetas : bool, default False
            The polar angles theta are not considered
        phis : bool, default False
            The polar angles phi are not considered

        Returns
        -------
        sites : (M, 1) |array_like|_
        thetas : (M, 1) |array_like|_
        phis : (M, 1) |array_like|_
        alphas : (M, 1) |array_lie|_
        """
        sites = np.random.randint(0,self.M, number_flippings)
        
        if phis_flag is False: 
            thetas = np.random.choice([0,np.pi], size=number_flippings)
            phis = np.zeros((number_flippings), dtype=float)
        elif thetas_flag is False:
            phis = np.random.uniform(0, 2*np.pi, number_flippings)
            thetas = np.zeros((number_flippings), dtype=float)
        else:
            raise ValueError(f"No angular degrees of freedom.")
        
        alphas = np.random.uniform(0, 2*np.pi, number_flippings)

        return sites, phis, thetas, alphas

    def metropolis_acceptance(self,spin_directions,site,theta,phi,alpha,old_E_0=None):
        r"""
        Metropolis acceptance move: evaluation of the new spin configuration and comparison of the relative energy with the one of the old spin configuration.

        Parameters
        ----------
        spin_directions: (M,3) |array_like|_
        site: int 
            Position of the rotated spin
        theta: float
            Polar coordinate of the axis of rotation
        phi: float
            Polar coordinate of the axis of rotation
        alpha: float
            Angle of rotation of the spin with respect to the the axis of rotation
        old_E_0: float
            Energy of a precedent spin configuration

        Returns
        -------
        new_spin_directions : (M, 3) |array_like|_
            New spin configuration
        new_E_0: float
            Energy related to the new spin configuration
        """
        u=np.array([np.sin(theta)*np.cos(phi),np.sin(theta)*np.sin(phi),np.cos(theta)])
        
        ux=np.zeros((3,3),dtype=float)
        ux[0,1]=-u[2]
        ux[0,2]=u[1]
        ux[1,2]=-u[0]
        ux=-ux.T
        
        uu=np.zeros((3,3),dtype=float)
        for i in range(3):
            for j in range(3):
                uu[i,j]=u[i]*u[j]

        R=np.cos(alpha)*np.eye(3)+np.sin(alpha)*ux+(1-np.cos(alpha))*uu

        new_spin_directions=spin_directions.copy()
        new_spin_directions[site]=R@new_spin_directions[site]

        if old_E_0 is None:
            old_E_0=self.E_0(spin_directions)
        
        new_E_0=self.E_0(new_spin_directions)

        if new_E_0<old_E_0:
            return new_spin_directions, new_E_0
        else:
            if np.random.rand()<np.exp(-self.beta*(new_E_0-old_E_0)):
                return new_spin_directions, new_E_0
        return spin_directions, old_E_0

    def warming_up(self, initial_guess=None, number_flippings=0, thetas_flag=False, phis_flag=False):
        r"""
        Warming up of the spin lattice, following Metropolis acceptance rules.

        Parameters
        ----------
        spin_directions: (M,3) |array_like|_
        site: int 
            Position of the rotated spin
        theta: float
            Polar coordinate of the axis of rotation
        phi: float
            Polar coordinate of the axis of rotation
        alpha: float
            Angle of rotation of the spin with respect to the the axis of rotation
        old_E_0: float
            Energy of a precedent spin configuration

        Returns
        -------
        new_spin_directions : (M, 3) |array_like|_
            New spin configuration
        """
        sites, phis, thetas, alphas=self.generate_possible_flippings(number_flippings, thetas_flag, phis_flag)
        
        if initial_guess is None:
            initial_guess = np.random.uniform(low=-1, high=1, size=(self.M, 3))

        sd = initial_guess / np.linalg.norm(initial_guess, axis=1)[:, np.newaxis]
    
        old_E_0=self.E_0(sd)

        for n in tqdm(range(number_flippings),desc="warming up..."):
            sd, old_E_0=self.metropolis_acceptance(sd,sites[n],thetas[n],phis[n],alphas[n],old_E_0)

        return initial_guess

    def measuring(self, initial_guess, number_flippings, thetas_flag=False, phis_flag=False):
        data = {
            'E': np.zeros(number_flippings,dtype=float),
            'E2' : np.zeros(number_flippings,dtype=float)}
        
        sites, phis, thetas, alphas=self.generate_possible_flippings(number_flippings, thetas_flag, phis_flag)
        
        if initial_guess is None:
            initial_guess = np.random.uniform(low=-1, high=1, size=(self.M, 3))

        sd = initial_guess / np.linalg.norm(initial_guess, axis=1)[:, np.newaxis]
    
        old_E_0=self.E_0(sd)

        for n in tqdm(range(number_flippings),desc="warming up..."):
            sd, old_E_0=self.metropolis_acceptance(sd,sites[n],thetas[n],phis[n],alphas[n],old_E_0)
            data['E']=old_E_0
            data['E2']=old_E_0**2

        return sd, data
    

if __name__ == "__main__":
    import magnopy.io as mio
    from magnopy.examples import cubic_ferro_nn, ivuzjo

    # spinham = cubic_ferro_nn(
    #     a=1,
    #     J_iso=1,
    #     J_21=(-1, 0, 0),
    #     S=0.5,
    #     dimensions=3,
    # )
    # spinham.add_magnetic_field(h=[0, 1, 0])

    spinham = ivuzjo(N=20)

    beta=1000
    number_flippings=100000

    energy = Energy_modified(spinham=spinham,beta=beta)

    warmed_sd=energy.warming_up(None,number_flippings)
    optimized_sd,data=energy.measuring(warmed_sd,number_flippings)

    print(optimized_sd)

    from mpl_toolkits.axes_grid1 import make_axes_locatable

    import matplotlib
    matplotlib.use("Agg")   # non-interactive backend (no GUI needed)
    import matplotlib.pyplot as plt

    
    fig, axs = plt.subplots(2, 2, figsize=(9, 9))
    fig.subplots_adjust(hspace=0.25, wspace=0.5)

    axs = axs.flatten()

    positions = np.array(spinham.atoms.positions)
    for i in range(3):
        im = axs[i].scatter(
            positions[:, 0],
            positions[:, 1],
            c=optimized_sd[:, i],
            vmin=-1,
            vmax=1,
            cmap="bwr",
        )
        divider = make_axes_locatable(axs[i])
        cax = divider.append_axes("right", size="5%", pad=0.05)
        axs[i].set_aspect(1)
        axs[i].set_title(f"$S_{'xyz'[i]}$")
        plt.colorbar(im, cax=cax)

    im = axs[3].quiver(
        positions[:, 0],
        positions[:, 1],
        0.5 * optimized_sd[:, 0] / np.linalg.norm(optimized_sd[:, :2], axis=1),
        0.5 * optimized_sd[:, 1] / np.linalg.norm(optimized_sd[:, :2], axis=1),
        optimized_sd[:, 2],
        angles="xy",
        scale_units="xy",
        scale=1,
        cmap="bwr",
        headlength=8,
        headaxislength=7,
        headwidth=5
    )

    # Now set normalization for colors explicitly:
    im.set_clim(vmin=-1, vmax=1)

    divider = make_axes_locatable(axs[3])
    cax = divider.append_axes("right", size="5%", pad=0.05)
    axs[3].set_aspect(1)
    axs[3].set_title(f"Vectors")
    plt.colorbar(im, cax=cax)

    plt.savefig("test.png", dpi=400, bbox_inches="tight")



