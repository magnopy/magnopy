.. _user-guide_usage_cell:

****
Cell
****

A cell is a set of three vectors that defines a periodic lattice in real space.

Magnopy stores it in the same way as many Python codes do (|spglib|_, |wulfric|_, ...).

``cell`` is a two-dimensional :math:`3\times3` matrix (|array-like|_) that groups
three lattice vectors as its rows. Here is an example of an orthorhombic cell

.. doctest::

    >>> cell = [
    ...     [3.5534, 0.0000, 0.0000],
    ...     [0.0000, 4.7449, 0.0000],
    ...     [0.0000, 0.0000, 8.7605],
    ... ]

with the three lattice vectors being

.. doctest::

    >>> cell[0] # a_1
    [3.5534, 0.0, 0.0]

.. doctest::

    >>> cell[1] # a_2
    [0.0, 4.7449, 0.0]

.. doctest::

    >>> cell[2] # a_3
    [0.0, 0.0, 8.7605]

``cell`` is one of the three objects that are required for creation of an instance of the
:py:class:`.SpinHamiltonian` class. It is then stored as an immutable attribute
:py:attr:`.SpinHamiltonian.cell`.

.. hint::
    Magnopy does not define any helper functions to manipulate the cell, as it is out
    of the scope of this package. Instead we depend on |wulfric|_ for all manipulations
    with the cell. For example, to compute the matrix of the cell from six lattice
    parameters one can use :py:func:`wulfric.cell.from_params` function.
