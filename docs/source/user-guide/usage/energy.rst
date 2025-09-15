.. _user-guide_usage_energy:

******
Energy
******

:py:class:`.Energy` class is used to compute the classical ground state energy of
the Hamiltonian. It is created based on some spin Hamiltonian

.. doctest::

    >>> import numpy as np
    >>> import magnopy
    >>> cell = np.eye(3)
    >>> atoms = {
    ...     "names" : ["Fe1"],
    ...     "species" : ["Fe"],
    ...     "positions" : [[0.0, 0.0, 0.0]],
    ...     "spins" : [5/2],
    ...     "g_factors" : [2]
    ... }
    >>> convention = magnopy.Convention(
    ...     multiple_counting=True, spin_normalized=False, c1=1, c21=1, c22=-1 / 2, c31=1, c41=1
    ... )
    >>> spinham = magnopy.SpinHamiltonian(cell=cell, atoms=atoms, convention=convention)
    >>> spinham.add_21(alpha=0, parameter=np.diag([2, -1, -1]))
    >>> spinham.add_22(alpha = 0, beta = 0, nu = (1, 0, 0), parameter = np.eye(3))
    >>> spinham.add_22(alpha = 0, beta = 0, nu = (0, 1, 0), parameter = np.eye(3))
    >>> spinham.add_22(alpha = 0, beta = 0, nu = (0, 0, 1), parameter = np.eye(3))
    >>> energy = magnopy.Energy(spinham=spinham)

And then takes a set of spin directions to compute :math:`E^{(0)}` (:py:meth:`.Energy.E_0`)

.. doctest::

    >>> energy.E_0(spin_directions = [[0, 0, 1]])
    -25.0
    >>> energy.E_0(spin_directions = [[0, 1, 0]])
    -25.0
    >>> energy.E_0(spin_directions = [[1, 0, 0]])
    -6.25
