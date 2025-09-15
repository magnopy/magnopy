.. _user-guide_usage_spin-hamiltonian:

****************
Spin Hamiltonian
****************

For the theoretical background on the spin Hamiltonian see
:ref:`user-guide_theory-behind_spin-hamiltonian`.

:py:class:`.SpinHamiltonian` class is at the heart of magnopy. Every calculation starts
with the definition of some spin Hamiltonian. This class stores the crystal structure,
convention and all parameters in it.

Note, that it stores the values of spins, but not the spin direction. The motivation
behind is that for any spin Hamiltonian the direction of spin vectors would be defined
by its parameters (it can be done by finding the ground state of the Hamiltonian).
Nevertheless, a calculation is possible for an arbitrary set of spin directions, that
are not necessary describe the ground state of the Hamiltonian, therefore, we decided
not to include them in this class to avoid confusion and unnecessary bookkeeping.


To create spin Hamiltonian one need three objects :ref:`user-guide_usage_cell`,
:ref:`user-guide_usage_atoms` and :ref:`user-guide_usage_convention`.

.. doctest::

    >>> import numpy as np
    >>> import magnopy
    >>> cell = np.eye(3)
    >>> atoms = {
    ...     "names" : ["Fe1", "Fe2"],
    ...     "species" : ["Fe", "Fe"],
    ...     "positions" : [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]],
    ...     "spins" : [5/2, 5/2],
    ...     "g_factors" : [2, 2]
    ... }
    >>> convention = magnopy.Convention(
    ...     multiple_counting=True, spin_normalized=False, c1=1, c21=1, c22=1 / 2, c31=1, c41=1
    ... )
    >>> spinham = magnopy.SpinHamiltonian(cell=cell, atoms=atoms, convention=convention)


Adding and removing parameters
==============================

Spin Hamiltonian stores parameters of the
:ref:`user-guide_theory-behind_spin-hamiltonian_expanded-form`. There are eleven types
of parameters. For each type there is a property that starts with ``p`` (i.e.
:py:attr:`.SpinHamiltonian.p1`), that provides access to the parameters of the
Hamiltonian and two functions that start with ``add_`` (i.e.
:py:meth:`.SpinHamiltonian.add_1`) or ``remove_`` (i.e.
:py:meth:`.SpinHamiltonian.remove_1`), that add or remove a parameter from the
Hamiltonian.

.. doctest::

    >>> import numpy as np
    >>> # Add on-site anisotropy (two spins & one site)
    >>> # Atoms are identified by their index in the spinham.atoms: 0 for Fe1
    >>> spinham.add_21(alpha=0, parameter=np.diag([2, -1, -1]))
    >>> # Add nearest-neighbor bilinear exchange (two spins & two sites)
    >>> spinham.add_22(alpha = 0, beta = 0, nu = (1, 0, 0), parameter = np.eye(3))
    >>> spinham.add_22(alpha = 0, beta = 0, nu = (0, 1, 0), parameter = np.eye(3))
    >>> spinham.add_22(alpha = 0, beta = 0, nu = (0, 0, 1), parameter = np.eye(3))

Each parameter property behaves as a list with parameters (technically it is either a
list or an iterator)

.. doctest::

    >>> for alpha, parameter in spinham.p21:
    ...     print(spinham.atoms.names[alpha], parameter, sep="\n")
    ...
    Fe1
    [[ 2  0  0]
     [ 0 -1  0]
     [ 0  0 -1]]

Note that there are 6 parameters in the ``p22``, as ``multiple_counting`` is ``True``

.. doctest::

    >>> for alpha, beta, nu, parameter in spinham.p22:
    ...     print(spinham.atoms.names[alpha], spinham.atoms.names[beta], nu)
    ...     print(parameter)
    ...
    Fe1 Fe1 (0, 0, 1)
    [[1. 0. 0.]
     [0. 1. 0.]
     [0. 0. 1.]]
    Fe1 Fe1 (0, 1, 0)
    [[1. 0. 0.]
     [0. 1. 0.]
     [0. 0. 1.]]
    Fe1 Fe1 (1, 0, 0)
    [[1. 0. 0.]
     [0. 1. 0.]
     [0. 0. 1.]]
    Fe1 Fe1 (0, 0, -1)
    [[1. 0. 0.]
     [0. 1. 0.]
     [0. 0. 1.]]
    Fe1 Fe1 (0, -1, 0)
    [[1. 0. 0.]
     [0. 1. 0.]
     [0. 0. 1.]]
    Fe1 Fe1 (-1, 0, 0)
    [[1. 0. 0.]
     [0. 1. 0.]
     [0. 0. 1.]]

Cell and atoms
==============

Spin Hamiltonian class stores cell :py:attr:`.SpinHamiltonian.cell` and atoms
:py:attr:`.SpinHamiltonian.atoms` as attributes

.. doctest::

    >>> spinham.cell
    array([[1., 0., 0.],
           [0., 1., 0.],
           [0., 0., 1.]])
    >>> spinham.atoms
    {'names': ['Fe1', 'Fe2'], 'species': ['Fe', 'Fe'], 'positions': [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]], 'spins': [2.5, 2.5], 'g_factors': [2, 2]}
    >>> # Magnopy adds syntactic sugar to the atoms dictionary inside the SpinHamiltonian class:
    >>> # a command
    >>> spinham.atoms.names
    ['Fe1', 'Fe2']
    >>> # is equivalent to
    >>> spinham.atoms["names"]
    ['Fe1', 'Fe2']
    >>> # It works with any key of atoms dictionary
    >>> spinham.atoms.spins
    [2.5, 2.5]

Cell and atoms are not meant to be changed once the Hamiltonian is created

.. doctest::

    >>> spinham.cell = 2 * np.eye(3)
    Traceback (most recent call last):
    ...
    AttributeError: Change of the cell attribute is not supported after the creation of SpinHamiltonian instance. If you need to modify cell, then use pre-defined methods of SpinHamiltonian or create a new one.
    >>> spinham.atoms = {}
    Traceback (most recent call last):
    ...
    AttributeError: Change of the atoms dictionary is not supported after the creation of SpinHamiltonian instance. If you need to modify atoms, then use pre-defined methods of SpinHamiltonian or create a new one.


Convention
==========

Convention of the Hamiltonian is stored as its attribute (:py:attr:`.SpinHamiltonian.convention`).

.. doctest::

    >>> spinham.convention.summary()
    custom convention where
      * Bonds are counted multiple times in the sum;
      * Spin vectors are not normalized;
      * c1 = 1.0;
      * c21 = 1.0;
      * c22 = 0.5;
      * c31 = 1.0;
      * Undefined c32 factor;
      * Undefined c33 factor;
      * c41 = 1.0;
      * Undefined c421 factor;
      * Undefined c422 factor;
      * Undefined c43 factor;
      * Undefined c44 factor.

The convention of the Hamiltonian can be changed. If the convention is being changed, then
the parameters will be adjusted accordingly. For example, if we change the numerical
factor before the two spins & one site term or remove multiple counting

.. doctest::

    >>> new_convention = spinham.convention.get_modified(multiple_counting=False)
    >>> spinham.convention = new_convention
    >>> for alpha, parameter in spinham.p21:
    ...     print(spinham.atoms.names[alpha], parameter, sep="\n")
    ...
    Fe1
    [[ 2  0  0]
     [ 0 -1  0]
     [ 0  0 -1]]
    >>> for alpha, beta, nu, parameter in spinham.p22:
    ...     print(spinham.atoms.names[alpha], spinham.atoms.names[beta], nu)
    ...     print(parameter)
    ...
    Fe1 Fe1 (0, 0, 1)
    [[2. 0. 0.]
     [0. 2. 0.]
     [0. 0. 2.]]
    Fe1 Fe1 (0, 1, 0)
    [[2. 0. 0.]
     [0. 2. 0.]
     [0. 0. 2.]]
    Fe1 Fe1 (1, 0, 0)
    [[2. 0. 0.]
     [0. 2. 0.]
     [0. 0. 2.]]

The main principle of changing convention can be formulated as "Energy of the Hamiltonian
should not change with its convention".

Magnetic vs non-magnetic atoms
==============================

Magnopy defines magnetic atom as an atom that has at least one parameter of the spin
Hamiltonian associated with it. Each spin Hamiltonian contains :math:`M` magnetic atoms
(:py:attr:`.SpinHamiltonian.M`). However, the crystal (cell & atoms) that are used for
the definition of the spin Hamiltonian can contain :math:`M^{\prime} \ne M` atoms.

Attribute :py:attr:`.SpinHamiltonian.atoms` returns a dictionary with all atoms of the
crystal, while :py:attr:`.SpinHamiltonian.magnetic_atoms` returns dictionary with only
magnetic atoms. The order of atoms is the same in both.

The indices in the specification of parameters correspond to the
:py:attr:`.SpinHamiltonian.atoms`. If you need to convert an index of
:py:attr:`.SpinHamiltonian.atoms` to an index of
:py:attr:`.SpinHamiltonian.magnetic_atoms` use the property :py:attr:`.SpinHamiltonian.map_to_magnetic`

.. doctest::

    >>> index_in_atoms = 0
    >>> index_in_magnetic_atoms = spinham.map_to_magnetic[index_in_atoms]

To convert from an index of :py:attr:`.SpinHamiltonian.magnetic_atoms` to the index of
:py:attr:`.SpinHamiltonian.atoms` use

.. doctest::

    >>> index_in_magnetic_atoms = 0
    >>> index_in_atoms = spinham.map_to_all[index_in_magnetic_atoms]


Now look at the example

.. doctest::

    >>> # Create a Hamiltonian with three atoms
    >>> atoms = dict(
    ... names=["Cr1", "Cr2", "Cr3"],
    ... spins = [3/2, 3/2, 3/2],
    ... g_factors=[2, 2, 2],
    ... positions=[[0, 0, 0],[0.5, 0, 0],[0, 0.5, 0]])
    >>> conv = magnopy.Convention(
    ... multiple_counting=True,
    ... spin_normalized=False,
    ... c21=1)
    >>> spinham = magnopy.SpinHamiltonian(
    ... cell=np.eye(3),
    ... atoms=atoms,
    ... convention=conv)

At this moment there is no magnetic atoms in the Hamiltonian (in the magnopy's context),
even though all atoms of the crystal have non-zero spin value.

.. doctest::

    >>> spinham.M
    0
    >>> spinham.magnetic_atoms
    {'names': [], 'spins': [], 'g_factors': [], 'positions': []}
    >>> spinham.map_to_magnetic
    [None, None, None]
    >>> spinham.map_to_all
    []

Lets add an on-site quadratic anisotropy to the second atom

.. doctest::

    >>> spinham.add_21(alpha=1, parameter = np.diag([1, 2, 3]))
    >>> spinham.M
    1
    >>> spinham.magnetic_atoms
    {'names': ['Cr2'], 'spins': [1.5], 'g_factors': [2], 'positions': [[0.5, 0, 0]]}
    >>> spinham.map_to_magnetic
    [None, 0, None]
    >>> spinham.map_to_all
    [1]
    >>> # Note that in the specification of the parameter the index
    >>> # of spinham.atoms is used
    >>> print(spinham.p21[0][0])
    1

Now second atom has a parameter associated with it, hence it is considered magnetic.
Two mapping lists can be used to safely convert from one to another

.. doctest::

    >>> # From index of magnetic atom to the index of non-magnetic atom
    >>> for i in range(spinham.M):
    ...     print(spinham.magnetic_atoms.names[i],
    ...     spinham.atoms.names[spinham.map_to_all[i]])
    ...
    Cr2 Cr2
    >>> # From index of non-magnetic atom to the index of magnetic atom
    >>> for i in range(len(spinham.atoms.names)):
    ...     if spinham.map_to_magnetic[i] is not None:
    ...         print(spinham.magnetic_atoms.names[spinham.map_to_magnetic[i]],
    ...               spinham.atoms.names[i])
    ...
    Cr2 Cr2
