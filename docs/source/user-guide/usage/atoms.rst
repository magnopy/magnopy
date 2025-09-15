.. _user-guide_usage_atoms:

*****
Atoms
*****

Cell and atoms are the two data structures that are implemented in a number of codes (if
not in each one of them). Magnopy does not re-implement them, but rather use the same
data structure as in another python package called |wulfric|_.

Here we briefly describe how the unit cell is stored, if you want to read the original
description of |wulfric|_ we refer you to |wulfric-key-concepts|_.

Atoms are stored as a plain python dictionary. Keys of the ``atoms`` are
properties of atoms. Values are the lists of :math:`M^{\prime}` elements each,
where :math:`M^{\prime}` is an amount of atoms.

.. doctest::

    >>> atoms = {
    ...     "names": ["Cr1", "Br1", "S1", "Cr2", "Br2", "S2"],
    ...     "species": ["Cr", "Br", "S", "Cr", "Br", "S"],
    ...     "positions": [
    ...         [0.500000, 0.000000, 0.882382],
    ...         [0.000000, 0.000000, 0.677322],
    ...         [0.500000, 0.500000, 0.935321],
    ...         [0.000000, 0.500000, 0.117618],
    ...         [0.500000, 0.500000, 0.322678],
    ...         [0.000000, 0.000000, 0.064679],
    ...     ],
    ... }


Keys recognized by magnopy:

*   "names" :
    ``list`` of ``str``. Inherited from |wulfric|_. This property is an informally a
    a main property and, as a rule of thumb, should be always present in the atoms
    dictionary.
*   "species" :
    ``list`` of ``str``. Inherited from |wulfric|_.
*   "positions" :
    ``list`` of **relative** coordinates of atoms. Inherited from |wulfric|_. Each
    element is an |array-like|_ of length :math:`3`.
*   "spins" :
    ``list`` of ``float``. Spin values for each atom.
*   "g_factors" :
    ``list`` of ``float``.

Magnopy extends the definition of ``atoms`` from |wulfric|_, therefore,
any function of |wulfric|_ may be used on it.

In magnopy atoms are used in the definition of
:ref:`user-guide_theory-behind_spin-hamiltonian`
and is one of three objects that are required for the creation of the
:py:class:`.SpinHamiltonian` class. It is stored as an its attribute
:py:attr:`.SpinHamiltonian.atoms`.

Magnetic vs non-magnetic atoms
==============================

Magnopy defines magnetic atom as an atom that has at least one parameter of the spin
Hamiltonian associated with it. Each spin Hamiltonian contains :math:`M` magnetic atoms
(:py:attr:`.SpinHamiltonian.M`). However, the crystal (cell & atoms) that are used for
the definition of the spin Hamiltonian can contain :math:`M^{\prime} \ne M` atoms.

Attribute :py:attr:`.SpinHamiltonian.atoms` returns a dictionary with all atoms of the
crystal, while :py:attr:`.SpinHamiltonian.magnetic_atoms` returns dictionary with only
magnetic atoms. The order of atoms is the same in both.
