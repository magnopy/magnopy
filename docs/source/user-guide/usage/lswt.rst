.. _user-guide_usage_lswt:

****
LSWT
****

All computations related to the linear spin wave theory (LSWT) are performed
by the :py:class:`.LSWT` class. It is created from a spin Hamiltonian. The parameters of
the spin Hamiltonian are optimized for the calculations of the LSWT upon creation of the
:py:class:`.LSWT` instance and the spin Hamiltonian is not stored within the
:py:class:`.LSWT` class.

.. doctest::

    >>> import numpy as np
    >>> import magnopy
    >>> # Hamiltonian for the nearest-neighbor ferromagnet on a cubic lattice
    >>> spinham = magnopy.examples.cubic_ferro_nn(S=5/2, J_iso=1, J_21=np.diag([2, -1, -2]))
    >>> lswt = magnopy.LSWT(spinham=spinham, spin_directions = [[0, 0, 1]])

Once created, it can be used to access the properties of the LSWT Hamiltonian.

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

Correction to the energy (:py:meth:`.LSWT.E_2`)

.. doctest::

    >>> lswt.E_2()
    -12.5

Coefficients before the one-operator terms (:py:meth:`.LSWT.O`)

.. doctest::

    >>> lswt.O()
    array([0.+0.j])

K-dependent properties
======================

Part of the spin Hamiltonian that depends on the wave vector :math:`\boldsymbol{k}`.
It can be computed one by one

.. doctest::

    >>> omega = lswt.omega(k = [0.5, 0, 0])
    >>> delta = lswt.delta(k = [0.5, 0, 0])

or all at once

.. doctest::

    >>> omega, delta, G = lswt.diagonalize(k = [0.5, 0, 0])

Note that call of :py:meth:`.LSWT.omega`, :py:meth:`.LSWT.delta` or :py:meth:`.LSWT.G`
invokes the call of :py:meth:`.LSWT.diagonalize`. Therefore, we recommend using
:py:meth:`.LSWT.diagonalize` to avoid duplicate calculations.
