.. _user-guide_usage_energy:

******
Energy
******

The :py:class:`.Energy` class is used as an interface to everything about the classical energy
of the Hamiltonian. It is created from a spin Hamiltonian

.. doctest::

    >>> import numpy as np
    >>> import magnopy
    >>> # Hamiltonian for the nearest-neighbor ferromagnet on a cubic lattice
    >>> spinham = magnopy.examples.cubic_ferro_nn(S=5/2, J_iso=1, J_21=np.diag([2, -1, -1]))
    >>> energy = magnopy.Energy(spinham=spinham)

And then takes a set of spin directions to compute :math:`E^{(0)}`
(:py:meth:`.Energy.E_0`).

Classical energy
================

The classical energy :math:`E^{(0)}` is computed for a given set of spin directions

.. doctest::

    >>> # In meV by default
    >>> energy.E_0(spin_directions = [[0, 0, 1]])
    -25.0
    >>> energy.E_0(spin_directions = [[0, 1, 0]])
    -25.0
    >>> energy.E_0(spin_directions = [[1, 0, 0]])
    -6.25

Optimal spin directions
=======================

The optimal spin directions (at least a local minimum) can be found by numerical
minimization

.. doctest::

    >>> optimal_directions = energy.optimize(initial_guess = [[0.5, 0.3, 0.1]]) # doctest: +SKIP
