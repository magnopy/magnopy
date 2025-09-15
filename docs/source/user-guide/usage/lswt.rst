.. _user-guide_usage_lswt:

****
LSWT
****

:py:class:`.LSWT` class is the main point for the computation at the level of
the linear spin wave theory. It is created from some spi Hamiltonian. The parameters of
the spin Hamiltonian are optimized for the calculations of the LSWT upon creation of the
:py:class:`.LSWT` instance and the spin Hamiltonian is not stored within the :py:class:`.LSWT`
class.

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
    >>> spinham.add_21(alpha=0, parameter=np.diag([2, -1, -2]))
    >>> spinham.add_22(alpha = 0, beta = 0, nu = (1, 0, 0), parameter = np.eye(3))
    >>> spinham.add_22(alpha = 0, beta = 0, nu = (0, 1, 0), parameter = np.eye(3))
    >>> spinham.add_22(alpha = 0, beta = 0, nu = (0, 0, 1), parameter = np.eye(3))
    >>> lswt = magnopy.LSWT(spinham=spinham, spin_directions = [[0, 0, 1]])

Once created, it can be used to access the properties of the LSWT Hamiltonian

.. doctest::

    >>> lswt.M
    1
    >>> lswt.z
    array([[0., 0., 1.]])
    >>> lswt.p
    array([[1.+0.j, 0.+1.j, 0.+0.j]])
    >>> lswt.spins
    array([2.5])
    >>> lswt.cell
    array([[1., 0., 0.],
           [0., 1., 0.],
           [0., 0., 1.]])

K-independent properties
========================

Correction to the energy (:py:attr:`.LSWT.E_2`)

.. doctest::

    >>> lswt.E_2
    -12.5

Coefficients before the one-operator terms (:py:attr:`.LSWT.O`)

.. doctest::

    >>> lswt.O
    array([0.+0.j])

K-dependant properties
======================

Part of the spin Hamiltonian that depends on the wave vector :math:`\boldsymbol{k}`.
They can be computed one by one

.. doctest::

    >>> omega = lswt.omega(k = [0.5, 0, 0])
    >>> delta = lswt.delta(k = [0.5, 0, 0])

or all at once

.. doctest::

    >>> omega, delta, G = lswt.diagonalize(k = [0.5, 0, 0])

Note that call of :py:meth:`.LSWT.omega`, :py:meth:`.LSWT.delta` or :py:meth:`.LSWT.G`
invokes the call of :py:meth:`.LSWT.diagonalize`. Therefore, we recommend to use
:py:meth:`.LSWT.diagonalize` to avoid duplicate calculations.
