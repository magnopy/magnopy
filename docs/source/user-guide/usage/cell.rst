.. _user-guide_usage_cell:

****
Cell
****

Cell and atoms are the two data structures that are implemented in a number of codes (if
not in each one of them). Magnopy does not re-implement them, but rather use the same
data structure as in another python package called |wulfric|_.

Here we briefly describe how the unit cell is stored, if you want to read the original
description of |wulfric|_ we refer you to |wulfric-key-concepts|_.

``cell`` is a two-dimensional :math:`3\times3` matrix (|array-like|_), that defines
three lattice vectors. The rows of the ``cell`` are vectors, while the columns are
cartesian coordinates. Here is an example of a simple orthorhombic cell

.. doctest::

    >>> cell = [
    ...     [3.553350, 0.000000, 0.000000],
    ...     [0.000000, 4.744935, 0.000000],
    ...     [0.000000, 0.000000, 8.760497],
    ... ]


with the three lattice vectors being

.. doctest::

    >>> a1 = cell[0]
    >>> a2 = cell[1]
    >>> a3 = cell[2]

In magnopy cell is used in the definition of
:ref:`user-guide_theory-behind_spin-hamiltonian`
and is one of three objects that are required for the creation of the
:py:class:`.SpinHamiltonian` class. It is stored as an its attribute
:py:attr:`.SpinHamiltonian.cell`.
