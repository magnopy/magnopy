.. _user-guide_usage_spin-directions:

***************
Spin directions
***************

Spin direction :math:`\boldsymbol{z}_1, ..., \boldsymbol{z}_M` is a rather simple data
structure: a list of 3D unit vectors (|array-like|_). It is used to define a state of the
Hamiltonian, for example its ground state.

The number of unit vectors in the list must be equal to the number of the
magnetic atoms in the unit cell of the spin Hamiltonian (i.e.
:py:attr:`.SpinHamiltonian.M`). The order of the unit vectors in the list must
correspond to the order of the atoms in
:py:attr:`.SpinHamiltonian.magnetic_atoms`.

* ``len(spin_directions) == spinham.M``
* ``spin_directions[i]`` is the spin direction (or directional vector) of the ``i``-th
  magnetic atom in the unit cell (i. e. of the ``spinham.magnetic_atoms.names[i]`` atom).


See also :ref:`user-guide_usage_spin-hamiltonian_magnetic-vs-non-magnetic`.

Here is an example of how to create a list of spin directions for a system with 3 atoms in
total.

First, we create a Hamiltonian with three atoms

.. doctest::

    >>> import numpy as np
    >>> import magnopy
    >>> # Create a Hamiltonian with three atoms
    >>> atoms = dict(
    ...     names=["Fe1", "Fe2", "Fe3"],
    ...     spins = [5/2, 5/2, 5/2],
    ...     g_factors=[2, 2, 2],
    ...     positions=[[0, 0, 0],[0.5, 0, 0],[0, 0.5, 0]]
    ... )
    >>> convention = magnopy.Convention(
    ...     multiple_counting=True,
    ...     spin_normalized=False,
    ...     c21=1
    ... )
    >>> spinham = magnopy.SpinHamiltonian(
    ...     cell=np.eye(3),
    ...     atoms=atoms,
    ...     convention=convention
    ... )

Then add interactions to the Hamiltonian that involve only "Fe1" and "Fe3" atoms

.. doctest::

    >>> # Add an on-site quadratic anisotropy to the first and third atom
    >>> spinham.add(nus = [(0, 0, 0)], alphas=[0, 0], parameter = np.diag([1, 2, 3]))
    >>> spinham.add(nus = [(0, 0, 0)], alphas=[2, 2], parameter = np.diag([1, 2, 3]))

Now there are two magnetic atoms in the Hamiltonian

.. doctest::

    >>> spinham.magnetic_atoms.names
    ['Fe1', 'Fe3']

Therefore, the list of spin directions must contain two unit vectors, for example

.. doctest::

    >>> spin_directions = [
    ...     [1, 0, 0],  # spin direction for Fe1
    ...     [0, 0, 1]   # spin direction for Fe3
    ... ]

If we add one more parameter that involves "Fe2" atom as well

.. doctest::

    >>> # Isotropic exchange interaction between the first and second atoms
    >>> spinham.add(nus=[(0, 0, 0)], alphas=[0, 1], parameter=np.eye(3))

Now there are three magnetic atoms in the Hamiltonian

.. doctest::

    >>> spinham.magnetic_atoms.names
    ['Fe1', 'Fe2', 'Fe3']

and the list of spin directions must contain three unit vectors

.. doctest::

    >>> spin_directions = [
    ...     [1, 0, 0],  # spin direction for Fe1
    ...     [0, 1, 0],  # spin direction for Fe2
    ...     [0, 0, 1]   # spin direction for Fe3
    ... ]
