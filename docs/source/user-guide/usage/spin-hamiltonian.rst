.. _user-guide_usage_spin-hamiltonian:

****************
Spin Hamiltonian
****************

.. hint::

    To visualize the one- and two-spin interaction of the Hamiltonian you can use the
    experimental function of Magnopy :py:func:`.experimental.plot_spinham`.



For the theoretical background on the spin Hamiltonian see
:ref:`user-guide_theory-behind_spin-hamiltonian`.

Spin Hamiltonian in Magnopy is an instance of the :py:class:`.SpinHamiltonian` class.

Three objects are required to create it: :ref:`user-guide_usage_cell`,
:ref:`user-guide_usage_atoms`, and :ref:`user-guide_usage_convention`.

.. doctest::

    >>> import numpy as np
    >>> import magnopy
    >>> cell = np.eye(3)
    >>> atoms = {
    ...     "names" : ["Fe1", "Fe2"],
    ...     "positions" : [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]],
    ...     "spins" : [5/2, 5/2],
    ...     "g_factors" : [2, 2],
    ...     "spglib_types" : [1, 1],
    ... }
    >>> convention = magnopy.Convention(
    ...     multiple_counting=True, spin_normalized=False, c1=1, c21=1, c22=1, c33=1, c45=1,
    ... )

Once those three objects are defined, the spin Hamiltonian can be created as

.. doctest::

    >>> spinham = magnopy.SpinHamiltonian(cell=cell, atoms=atoms, convention=convention)

The cell and atoms of the Hamiltonian can be viewed at any time as

.. doctest::

    >>> spinham.cell
    array([[1., 0., 0.],
           [0., 1., 0.],
           [0., 0., 1.]])
    >>> spinham.atoms
    {'names': ['Fe1', 'Fe2'], 'positions': [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]], 'spins': [2.5, 2.5], 'g_factors': [2, 2], 'spglib_types': [1, 1]}

However, they cannot be changed, to avoid inconsistencies later on

.. doctest::

    >>> spinham.cell = 2 * np.eye(3)
    Traceback (most recent call last):
    ...
    AttributeError: Change of the cell attribute is not allowed after the creation of SpinHamiltonian instance. SpinHamiltonian.cell is immutable.
    >>> spinham.atoms = {}
    Traceback (most recent call last):
    ...
    AttributeError: Change of the atoms dictionary is not allowed after the creation of SpinHamiltonian instance. SpinHamiltonian.atoms is immutable.

.. hint::

    .. doctest::

        >>> spinham.atoms.names
        ['Fe1', 'Fe2']

    is equivalent to

    .. doctest::

        >>> spinham.atoms["names"]
        ['Fe1', 'Fe2']

Accessing the parameters
========================

First of all, one needs to be able to access the parameters of the Hamiltonian. The main
method for that task is :py:meth:`.SpinHamiltonian.parameters`. It returns an iterator
over the parameters of the Hamiltonian:

.. doctest::

    >>> # For the explanation of nus, alphas, and parameter see next sections
    >>> for nus, alphas, parameter in spinham.parameters():
    ...     print(nus, alphas, parameter)

Each element of the iterator is a tuple ``(nus, alphas, parameter)`` with the following
properties

* ``len(alphas) == n``
* ``len(nus) == n - 1``
* atom ``alphas[0]`` is located in the ``(0, 0, 0)`` unit cell.
* atom ``alphas[i]`` is located in the ``nus[i-1]`` unit cell for all ``1 <= i < n``

``nus`` is a tuple of indices :math:`(\nu_2, ..., \nu_n)`, ``alphas`` is a tuple
of indices :math:`(\alpha_1, ..., \alpha_n)`, ``parameter`` is the
vector/matrix/tensor of the interaction parameter :math:`J^{i_1, ...,
i_n}_{\nu_2, ..., \nu_n; \alpha_1, ..., \alpha_n}`, where ``n`` is the number of
spin operators involved in the term. More on the meaning of ``nus``, ``alphas``,
and ``parameter`` can be found in the next sections.

:py:meth:`.SpinHamiltonian.parameters` can take two optional arguments ``n`` and
``p_n`` that filter the parameters by eleven types as described in
:ref:`user-guide_theory-behind_spin-hamiltonian` page.

*   ``n`` selects the number of spin operators in the term of the Hamiltonian.
*   ``p_n`` for each ``n`` selects one of the sub-types of parameters. For example, for
    ``n=2`` there are two sub-types of parameters: ``p_n = 1`` and ``p_n = 2``. The first
    one filters interaction where both spin operators reside at the exact same spot in the
    real space. The second one filters interactions where spin operators are located at
    different positions in the real space.

.. hint::
    There are eleven properties of :py:class:`.SpinHamiltonian` that are equivalent to
    the call of :py:meth:`.SpinHamiltonian.parameters` with the right values of ``n`` and
    ``p_n`` as summarized in the table below.

    =============================== ==========================================
    property                        same as
    =============================== ==========================================
    :py:attr:`.SpinHamiltonian.p1`  ``SpinHamiltonian.parameters(n=1, p_n=1)``
    :py:attr:`.SpinHamiltonian.p21` ``SpinHamiltonian.parameters(n=2, p_n=1)``
    :py:attr:`.SpinHamiltonian.p22` ``SpinHamiltonian.parameters(n=2, p_n=2)``
    :py:attr:`.SpinHamiltonian.p31` ``SpinHamiltonian.parameters(n=3, p_n=1)``
    :py:attr:`.SpinHamiltonian.p32` ``SpinHamiltonian.parameters(n=3, p_n=2)``
    :py:attr:`.SpinHamiltonian.p33` ``SpinHamiltonian.parameters(n=3, p_n=3)``
    :py:attr:`.SpinHamiltonian.p41` ``SpinHamiltonian.parameters(n=4, p_n=1)``
    :py:attr:`.SpinHamiltonian.p42` ``SpinHamiltonian.parameters(n=4, p_n=2)``
    :py:attr:`.SpinHamiltonian.p43` ``SpinHamiltonian.parameters(n=4, p_n=3)``
    :py:attr:`.SpinHamiltonian.p44` ``SpinHamiltonian.parameters(n=4, p_n=4)``
    :py:attr:`.SpinHamiltonian.p45` ``SpinHamiltonian.parameters(n=4, p_n=5)``
    =============================== ==========================================

So far the spin Hamiltonian that we created does not have any parameters in it

.. doctest::

    >>> len(spinham.parameters())
    0

.. _user-guide_usage_spin-hamiltonian_adding-parameters:

Adding parameters
=================

The method :py:meth:`.SpinHamiltonian.add` can be used to add any interaction parameter to
the spin Hamiltonian.

This method follows the notation of the algebraic form of the Hamiltonian described in
the :ref:`user-guide_theory-behind_spin-hamiltonian` page and expects two bits of
information

* Positions of all spin operators involved in a term:

  - ``nus``, that are expected to be either :math:`(\mu, \mu+\nu_2, ..., \mu+\nu_n)` or
    :math:`(\nu_2, ..., \nu_n)`. More on those two options in the
    :ref:`user-guide_usage_spin-hamiltonian_translational-symmetry` section.

  - ``alphas``, that are expected to be :math:`(\alpha_1, ..., \alpha_n)`.

* Vector/matrix/tensor of the interaction parameter: ``parameter``.



Let us give examples for every value of ``1 <= n <= 4``.

:math:`n = 1`
-------------

The algebraic form of a single term (omitting the convention constant) is

.. math::

    J^{i_1}_{\alpha_1}
    S_{\mu, \alpha_1}^{i_1}

Let us pick some values for the indices and parameters

.. doctest::

    >>> alpha_1 = 0
    >>> mu = (0, 0, 0)
    >>> parameter = np.array([1, 2, 3])

Then the term can be added to the Hamiltonian in two equivalent ways

.. doctest::

    >>> spinham.add(nus = [mu], alphas = [alpha_1], parameter = parameter)

.. doctest::
    :hide:

    >>> spinham.remove(nus = [mu], alphas = [alpha_1])

or

.. doctest::

    >>> spinham.add(nus = [], alphas = [alpha_1], parameter = parameter)

:math:`n = 2`
-------------

The algebraic form of a single term (omitting the convention constant) is

.. math::

    J^{i_1, i_2}_{\nu_2; \alpha_1, \alpha_2}
    S_{\mu, \alpha_1}^{i_1}
    S_{\mu + \nu_2, \alpha_2}^{i_2}

Let us pick some values for the indices and parameters

.. doctest::

    >>> alpha_1 = 0
    >>> alpha_2 = 1
    >>> mu = (0, 0, 0)
    >>> nu_2 = (-1, -1, -1)
    >>> parameter = np.eye(3)

Then the term can be added to the Hamiltonian in two equivalent ways

.. doctest::

    >>> mu_plus_nu_2 = tuple([mu[i] + nu_2[i] for i in range(3)])
    >>> spinham.add(
    ...     nus = [mu, mu_plus_nu_2],
    ...     alphas = [alpha_1, alpha_2],
    ...     parameter = parameter,
    ... )

.. doctest::
    :hide:

    >>> spinham.remove(nus = [mu, mu_plus_nu_2], alphas = [alpha_1, alpha_2])

or

.. doctest::

    >>> spinham.add(
    ...     nus = [nu_2],
    ...     alphas = [alpha_1, alpha_2],
    ...     parameter = parameter,
    ... )

:math:`n = 3`
-------------

The algebraic form of a single term (omitting the convention constant) is

.. math::

    J^{i_1, i_2, i_3}_{\nu_2, \nu_3; \alpha_1, \alpha_2, \alpha_3}
    S_{\mu, \alpha_1}^{i_1}
    S_{\mu + \nu_2, \alpha_1}^{i_2}
    S_{\mu + \nu_3, \alpha_3}^{i_3}

Let us pick some values for the indices and parameters

.. doctest::

    >>> alpha_1 = 0
    >>> alpha_2 = 1
    >>> alpha_3 = 0
    >>> mu = (0, 0, 0)
    >>> nu_2 = (-1, -1, -1)
    >>> nu_3 = (1, 0, 0)
    >>> parameter = np.ones((3, 3, 3))

Then the term can be added to the Hamiltonian in two equivalent ways

.. doctest::

    >>> mu_plus_nu_2 = tuple([mu[i] + nu_2[i] for i in range(3)])
    >>> mu_plus_nu_3 = tuple([mu[i] + nu_3[i] for i in range(3)])
    >>> spinham.add(
    ...     nus = [mu, mu_plus_nu_2, mu_plus_nu_3],
    ...     alphas = [alpha_1, alpha_2, alpha_3],
    ...     parameter = parameter
    ... )

.. doctest::
    :hide:

    >>> spinham.remove(
    ...     nus = [mu, mu_plus_nu_2, mu_plus_nu_3],
    ...     alphas = [alpha_1, alpha_2, alpha_3]
    ... )

or

.. doctest::

    >>> spinham.add(
    ...     nus = [nu_2, nu_3],
    ...     alphas = [alpha_1, alpha_2, alpha_3],
    ...     parameter = parameter
    ... )

:math:`n = 4`
-------------

The algebraic form of a single term (omitting the convention constant) is

.. math::

    J^{i_1, i_2, i_3, i_4}_{\nu_2, \nu_3, \nu_4; \alpha_1, \alpha_2, \alpha_3, \alpha_4}
    S_{\mu, \alpha_1}^{i_1}
    S_{\mu + \nu_2, \alpha_2}^{i_2}
    S_{\mu + \nu_3, \alpha_3}^{i_3}
    S_{\mu + \nu_4, \alpha_4}^{i_4}

Let us pick some values for the indices and parameters

.. doctest::

    >>> alpha_1 = 0
    >>> alpha_2 = 1
    >>> alpha_3 = 0
    >>> alpha_4 = 1
    >>> mu = (0, 0, 0)
    >>> nu_2 = (-1, -1, -1)
    >>> nu_3 = (1, 0, 0)
    >>> nu_4 = (0, 1, 0)
    >>> parameter = np.ones((3, 3, 3, 3))

Then the term can be added to the Hamiltonian in two equivalent ways

.. doctest::

    >>> mu_plus_nu_2 = tuple([mu[i] + nu_2[i] for i in range(3)])
    >>> mu_plus_nu_3 = tuple([mu[i] + nu_3[i] for i in range(3)])
    >>> mu_plus_nu_4 = tuple([mu[i] + nu_4[i] for i in range(3)])
    >>> spinham.add(
    ...     nus = [mu, mu_plus_nu_2, mu_plus_nu_3, mu_plus_nu_4],
    ...     alphas = [alpha_1, alpha_2, alpha_3, alpha_4],
    ...     parameter = parameter
    ... )

.. doctest::
    :hide:

    >>> spinham.remove(
    ...     nus = [mu, mu_plus_nu_2, mu_plus_nu_3, mu_plus_nu_4],
    ...     alphas = [alpha_1, alpha_2, alpha_3, alpha_4]
    ... )

or

.. doctest::

    >>> spinham.add(
    ...     nus = [nu_2, nu_3, nu_4],
    ...     alphas = [alpha_1, alpha_2, alpha_3, alpha_4],
    ...     parameter = parameter
    ... )




.. _user-guide_usage_spin-hamiltonian_translational-symmetry:

Translational symmetry
======================

As you can see in the :ref:`user-guide_usage_spin-hamiltonian_adding-parameters` section,
there are two equivalent methods to add a parameter to the Hamiltonian due to the
translational symmetry of the underlying lattice.

Note that the interaction parameters do not depend on the unit cell index :math:`\mu`.
Therefore, it is enough to store only the parameters with a single fixed
value of the unit cell index :math:`\mu`. We choose that value to be :math:`(0, 0, 0)`.

Therefore, when the user provides ``nus`` with the length of ``n - 1``, Magnopy interprets
``nus`` as :math:`(\nu_2, ..., \nu_n)` and :math:`\mu = (0, 0, 0)` is implied.

However, when the user provides ``nus`` with the length of ``n``, Magnopy interprets
``nus`` as :math:`(\mu, \mu+\nu_2, ..., \mu+\nu_n)`. In this case the user is free
to use any value for :math:`\mu`. Magnopy, will automatically shift all elements of
``nus`` to enforce :math:`\mu = (0, 0, 0)` (i.e. ``nus[0] == (0, 0, 0)``).

Let us demonstrate the latter with an example. First, we create two copies of the spin
Hamiltonian with no parameters in it

.. doctest::

    >>> spinham_v1 = spinham.get_empty()
    >>> spinham_v2 = spinham.get_empty()

Then we add single interaction to the first Hamiltonian

.. doctest::

    >>> spinham_v1.add(nus=[(0,0,0), (1,0,0)], alphas = [0,1], parameter = np.eye(3))

and add transitionally equivalent interaction to the second Hamiltonian

.. doctest::

    >>> spinham_v2.add(nus=[(-3, 34, 21), (-2, 34, 21)], alphas = [0,1], parameter = np.eye(3))

Now both Hamiltonians have the same parameter stored, with :math:`\mu` forced to be
:math:`(0, 0, 0)`.

.. doctest::

    >>> for nus, alphas, parameter in spinham_v1.p22:
    ...     print(nus, alphas, parameter, sep="\n")
    ((1, 0, 0),)
    (0, 1)
    [[1. 0. 0.]
     [0. 1. 0.]
     [0. 0. 1.]]
    >>> for nus, alphas, parameter in spinham_v2.p22:
    ...     print(nus, alphas, parameter, sep="\n")
    ((1, 0, 0),)
    (0, 1)
    [[1. 0. 0.]
     [0. 1. 0.]
     [0. 0. 1.]]


.. _user-guide_usage_spin-hamiltonian_equivalent-parameters:

Equivalent parameters
=====================

See :ref:`user-guide_theory-behind_equivalent-parameters` for the introduction to the
concept of equivalent parameters.

We use the term with two spin operators that are located at different positions in real
space as an example in this section. Similar logic applies to other terms of the
Hamiltonian.

First, we create an empty copy of the Hamiltonian

.. doctest::

    >>> spinham = spinham.get_empty()

Equivalent set contains two parameters in the case of the term with two spin operators
located at different positions. Consider the interaction from the "Fe1" atom in the unit
cell :math:`(0, 0, 0)` to the "Fe2" atom in the unit cell :math:`(-1, -1, -1)` with the
matrix of the interaction parameter

.. doctest::

    >>> parameter1 = [
    ...     [1, 0.5, 0],
    ...     [-0.5, 1, 0],
    ...     [0, 0, 1],
    ... ]

Its equivalent interaction is the one from the "Fe2" atom in the unit cell
:math:`(0, 0, 0)` to the "Fe1" atom in the unit cell :math:`(1, 1, 1)` with the matrix of
the interaction parameter

.. doctest::

    >>> # Note that parameter2 == parameter1.T
    >>> parameter2 = [
    ...     [1, -0.5, 0],
    ...     [0.5, 1, 0],
    ...     [0, 0, 1],
    ... ]

The behavior of Magnopy depends on the value of ``multiple_counting`` property of the
Hamiltonian's convention.

multiple_counting = True
------------------------

.. doctest::

    >>> spinham.convention.multiple_counting
    True

In this case the user is expected to add both parameters to the Hamiltonian by hand

.. doctest::

    >>> spinham.add(
    ...     nus = [(-1, -1, -1)],
    ...     alphas = [0, 1],
    ...     parameter = parameter1
    ... )
    >>> spinham.add(
    ...     nus = [(1, 1, 1)],
    ...     alphas = [1, 0],
    ...     parameter = parameter2
    ... )

The amount of equivalent parameters in the set grows with the increasing amount of
involved spin operators. For example, sets of the :ref:`ug_tb_sh_4-5` contain 24
interactions. To avoid the necessity of adding all of them by hand you can pass the
keyword argument ``populate_equivalent=True`` to the :py:meth:`.SpinHamiltonian.add`
method.

.. doctest::

    >>> spinham_pe = spinham.get_empty()
    >>> spinham_pe.add(
    ...     nus = [(-1, -1, -1)],
    ...     alphas = [0, 1],
    ...     parameter = parameter1,
    ...     populate_equivalent=True
    ... )

The results are equivalent

.. doctest::

    >>> for nus, alphas, parameter in spinham.p22:
    ...     print(spinham.atoms.names[alphas[0]], spinham.atoms.names[alphas[1]], nus[0])
    ...     print(nus, alphas, parameter, sep="\n")
    Fe1 Fe2 (-1, -1, -1)
    ((-1, -1, -1),)
    (0, 1)
    [[ 1.   0.5  0. ]
     [-0.5  1.   0. ]
     [ 0.   0.   1. ]]
    Fe2 Fe1 (1, 1, 1)
    ((1, 1, 1),)
    (1, 0)
    [[ 1.  -0.5  0. ]
     [ 0.5  1.   0. ]
     [ 0.   0.   1. ]]

.. doctest::

    >>> for nus, alphas, parameter in spinham_pe.p22:
    ...     print(spinham.atoms.names[alphas[0]], spinham.atoms.names[alphas[1]], nus[0])
    ...     print(nus, alphas, parameter, sep="\n")
    Fe1 Fe2 (-1, -1, -1)
    ((-1, -1, -1),)
    (0, 1)
    [[ 1.   0.5  0. ]
     [-0.5  1.   0. ]
     [ 0.   0.   1. ]]
    Fe2 Fe1 (1, 1, 1)
    ((1, 1, 1),)
    (1, 0)
    [[ 1.  -0.5  0. ]
     [ 0.5  1.   0. ]
     [ 0.   0.   1. ]]

multiple_counting = False
-------------------------

.. doctest::

    >>> spinham = spinham.get_empty()
    >>> spinham.convention = spinham.convention.get_modified(multiple_counting=False)
    >>> spinham.convention.multiple_counting
    False

In the user is expected to add only one parameter from the set to the Hamiltonian (any
parameter from the set will work).

.. doctest::

    >>> spinham.add(
    ...     nus = [(-1, -1, -1)],
    ...     alphas = [0, 1],
    ...     parameter = np.array(parameter1) + np.array(parameter2).T
    ... )

Then, if you try to add the second parameter from the set, Magnopy will raise an error
as the representative parameter of the set is already added to the Hamiltonian

.. doctest::

    >>> spinham.add(
    ...     nus = [(1, 1, 1)],
    ...     alphas = [1, 0],
    ...     parameter = np.array(parameter2) + np.array(parameter1).T
    ... )
    Traceback (most recent call last):
    ...
    ValueError: Parameter with such specs is already present.

Note that if you change the convention back to ``multiple_counting=True``, the parameters
are the same as before

.. doctest::

    >>> spinham.convention = spinham.convention.get_modified(multiple_counting=True)
    >>> for nus, alphas, parameter in spinham.p22:
    ...     print(spinham.atoms.names[alphas[0]], spinham.atoms.names[alphas[1]], nus[0])
    ...     print(nus, alphas, parameter, sep="\n")
    Fe1 Fe2 (-1, -1, -1)
    ((-1, -1, -1),)
    (0, 1)
    [[ 1.   0.5  0. ]
     [-0.5  1.   0. ]
     [ 0.   0.   1. ]]
    Fe2 Fe1 (1, 1, 1)
    ((1, 1, 1),)
    (1, 0)
    [[ 1.  -0.5  0. ]
     [ 0.5  1.   0. ]
     [ 0.   0.   1. ]]


.. _user-guide_usage_spin-hamiltonian_symmetrization:

Symmetrization
==============

See :ref:`user-guide_theory-behind_equivalent-parameters` for the introduction to the
concept of equivalent parameters and their symmetrization.

This section is only relevant when ``multiple_counting=True``. In that case, if you
do not use ``populate_equivalent=True`` while adding parameters to the Hamiltonian with
:py:meth:`.SpinHamiltonian.add` method, then you a free to add non-symmetrized
matrices/tensors of parameters to the Hamiltonian. Magnopy silently symmetrizes those
parameters when necessary. In other words, while the Hamiltonian can contain
non-symmetrized parameters, we **do not guarantee** that they stay non-symmetrized if some
operations are performed with the Hamiltonian.

Let us illustrate that with an example. First, we create an empty copy of the Hamiltonian

.. doctest::

    >>> spinham = spinham.get_empty()
    >>> spinham.convention = spinham.convention.get_modified(multiple_counting=True)

Then we consider the same pair of equivalent parameters as in the previous section:

An interaction from the "Fe1" atom in the unit cell :math:`(0, 0, 0)` to the "Fe2" atom in
the unit cell :math:`(-1, -1, -1)` and an interaction from the "Fe2" atom in the unit cell
:math:`(0, 0, 0)` to the "Fe1" atom in the unit cell :math:`(1, 1, 1)`. However, this time
we choose the matrices of interaction parameters to be non-symmetrized (while keeping the
same physics as in the previous section)

.. doctest::

    >>> parameter1 = [
    ...     [0, 0.3, 0],
    ...     [-0.5, 1, 0],
    ...     [0, 0, 1.68],
    ... ]
    >>> parameter2 = [
    ...     [2, -0.5, 0],
    ...     [0.7, 1, 0],
    ...     [0, 0, 0.32],
    ... ]

Now we add those parameters to the Hamiltonian

.. doctest::

    >>> spinham.add(
    ...     nus = [(-1, -1, -1)],
    ...     alphas = [0, 1],
    ...     parameter = parameter1
    ... )
    >>> spinham.add(
    ...     nus = [(1, 1, 1)],
    ...     alphas = [1, 0],
    ...     parameter = parameter2
    ... )

And check that the parameters are indeed non-symmetrized

.. doctest::

    >>> for nus, alphas, parameter in spinham.p22:
    ...     print(spinham.atoms.names[alphas[0]], spinham.atoms.names[alphas[1]], nus[0])
    ...     print(nus, alphas, parameter, sep="\n")
    Fe1 Fe2 (-1, -1, -1)
    ((-1, -1, -1),)
    (0, 1)
    [[ 0.    0.3   0.  ]
     [-0.5   1.    0.  ]
     [ 0.    0.    1.68]]
    Fe2 Fe1 (1, 1, 1)
    ((1, 1, 1),)
    (1, 0)
    [[ 2.   -0.5   0.  ]
     [ 0.7   1.    0.  ]
     [ 0.    0.    0.32]]

Now one can symmetrize those parameters by calling the method
:py:meth:`.SpinHamiltonian.set_distribution`.

.. doctest::

    >>> spinham.set_distribution(strategy="symmetrize")

And now the parameters are symmetrized and are the same as in the previous section

.. doctest::

    >>> for nus, alphas, parameter in spinham.p22:
    ...     print(spinham.atoms.names[alphas[0]], spinham.atoms.names[alphas[1]], nus[0])
    ...     print(nus, alphas, parameter, sep="\n")
    Fe1 Fe2 (-1, -1, -1)
    ((-1, -1, -1),)
    (0, 1)
    [[ 1.   0.5  0. ]
     [-0.5  1.   0. ]
     [ 0.   0.   1. ]]
    Fe2 Fe1 (1, 1, 1)
    ((1, 1, 1),)
    (1, 0)
    [[ 1.  -0.5  0. ]
     [ 0.5  1.   0. ]
     [ 0.   0.   1. ]]

.. note::
    :py:meth:`.SpinHamiltonian.set_distribution` is different from using
    ``populate_equivalent=True`` in the :py:meth:`.SpinHamiltonian.add`. The former
    takes the parameters that are already present in the Hamiltonian and redistribute
    them, while the latter adds extra parameters to the Hamiltonian.

    Here is an example

    .. doctest::

        >>> spinham_v1 = spinham.get_empty()
        >>> spinham_v2 = spinham.get_empty()
        >>> parameter = [
        ...     [0, 0.3, 0],
        ...     [-0.5, 1, 0],
        ...     [0, 0, 2],
        ... ]

    In the first case we add a parameter and then symmetrize the Hamiltonian

    .. doctest::

        >>> spinham_v1.add(
        ...     nus=[(-1, -1, -1)],
        ...     alphas = [0, 1],
        ...     parameter = parameter,
        ... )
        >>> spinham_v1.set_distribution(strategy="symmetrize")

    In the second case we request the code to populate equivalent parameters

    .. doctest::

        >>> spinham_v2.add(
        ...     nus=[(-1, -1, -1)],
        ...     alphas = [0, 1],
        ...     parameter = parameter,
        ...     populate_equivalent=True,
        ... )

    In both cases there are two interaction parameters

    .. doctest::

        >>> len(spinham_v1.p22)
        2
        >>> len(spinham_v2.p22)
        2

    However, the values of the parameters are different (meaning that two Hamiltonians
    describe different physics)

    .. doctest::

        >>> for nus, alphas, parameter in spinham_v1.p22:
        ...     print(spinham.atoms.names[alphas[0]], spinham.atoms.names[alphas[1]], nus[0])
        ...     print(nus, alphas, parameter, sep="\n")
        Fe1 Fe2 (-1, -1, -1)
        ((-1, -1, -1),)
        (0, 1)
        [[ 0.    0.15  0.  ]
         [-0.25  0.5   0.  ]
         [ 0.    0.    1.  ]]
        Fe2 Fe1 (1, 1, 1)
        ((1, 1, 1),)
        (1, 0)
        [[ 0.   -0.25  0.  ]
         [ 0.15  0.5   0.  ]
         [ 0.    0.    1.  ]]

    .. doctest::

        >>> for nus, alphas, parameter in spinham_v2.p22:
        ...     print(spinham.atoms.names[alphas[0]], spinham.atoms.names[alphas[1]], nus[0])
        ...     print(nus, alphas, parameter, sep="\n")
        Fe1 Fe2 (-1, -1, -1)
        ((-1, -1, -1),)
        (0, 1)
        [[ 0.   0.3  0. ]
         [-0.5  1.   0. ]
         [ 0.   0.   2. ]]
        Fe2 Fe1 (1, 1, 1)
        ((1, 1, 1),)
        (1, 0)
        [[ 0.  -0.5  0. ]
         [ 0.3  1.   0. ]
         [ 0.   0.   2. ]]


.. _user-guide_usage_spin-hamiltonian_removing-parameters:

Removing a parameter
====================

To remove a parameter from the Hamiltonian use :py:meth:`.SpinHamiltonian.remove` method
that expects ``nus`` and ``alphas`` as arguments.

First, we add a single parameter to the Hamiltonian

.. doctest::

    >>> spinham = spinham.get_empty()
    >>> spinham.add(
    ...     nus = [],
    ...     alphas = [0],
    ...     parameter = [1,2,3]
    ... )

There is one parameter in the Hamiltonian now

.. doctest::

    >>> len(spinham.parameters())
    1

Now we can remove it

.. doctest::

    >>> spinham.remove(
    ...     nus = [],
    ...     alphas = [0],
    ... )

And there is no parameters in the Hamiltonian anymore

.. doctest::

    >>> len(spinham.parameters())
    0


.. _user-guide_usage_spin-hamiltonian_convention:

Convention
==========

Convention of the Hamiltonian is stored as its attribute
(:py:attr:`.SpinHamiltonian.convention`).

.. doctest::

    >>> print(spinham.convention)
    "custom" convention where
      * Bonds are counted multiple times in the sum;
      * Spin vectors are not normalized;
      * c1 = 1.0;
      * c21 = 1.0;
      * c22 = 1.0;
      * Undefined c31 factor;
      * Undefined c32 factor;
      * c33 = 1.0;
      * Undefined c41 factor;
      * Undefined c42 factor;
      * Undefined c43 factor;
      * Undefined c44 factor;
      * c45 = 1.0.

The convention of the Hamiltonian can be changed. When the convention is being changed,
the parameters adjust accordingly.

First, we add some parameters to the Hamiltonian

.. doctest::

    >>> spinham = spinham.get_empty()
    >>> # On-site anisotropy
    >>> spinham.add(nus=[(0,0,0)], alphas=[0,0], parameter = np.diag([2, 1, 1]))
    >>> spinham.add(nus=[(0,0,0)], alphas=[1,1], parameter = np.diag([2, 1, 1]))
    >>> # Some of the nearest-neighbors exchange bonds
    >>> spinham.add(nus=[(0,0,0)], alphas=[0,1], parameter = np.eye(3), populate_equivalent=True)


Now you can change the constant before the sum

.. doctest::

    >>> for nus, alphas, parameter in spinham.p21:
    ...     print(spinham.atoms.names[alphas[0]], parameter, sep="\n")
    Fe1
    [[2. 0. 0.]
     [0. 1. 0.]
     [0. 0. 1.]]
    Fe2
    [[2. 0. 0.]
     [0. 1. 0.]
     [0. 0. 1.]]
    >>> spinham.convention = spinham.convention.get_modified(c21=2)
    >>> for nus, alphas, parameter in spinham.p21:
    ...     print(spinham.atoms.names[alphas[0]], parameter, sep="\n")
    Fe1
    [[1.  0.  0. ]
     [0.  0.5 0. ]
     [0.  0.  0.5]]
    Fe2
    [[1.  0.  0. ]
     [0.  0.5 0. ]
     [0.  0.  0.5]]

Or the multiple counting of the parameters

.. doctest::

    >>> len(spinham.p22)
    2
    >>> for nus, alphas, parameter in spinham.p22:
    ...     print(spinham.atoms.names[alphas[0]], spinham.atoms.names[alphas[1]], nus[0])
    ...     print(nus, alphas, parameter, sep="\n")
    Fe1 Fe2 (0, 0, 0)
    ((0, 0, 0),)
    (0, 1)
    [[1. 0. 0.]
     [0. 1. 0.]
     [0. 0. 1.]]
    Fe2 Fe1 (0, 0, 0)
    ((0, 0, 0),)
    (1, 0)
    [[1. 0. 0.]
     [0. 1. 0.]
     [0. 0. 1.]]
    >>> spinham.convention = spinham.convention.get_modified(multiple_counting=False)
    >>> len(spinham.p22)
    1
    >>> for nus, alphas, parameter in spinham.p22:
    ...     print(spinham.atoms.names[alphas[0]], spinham.atoms.names[alphas[1]], nus[0])
    ...     print(nus, alphas, parameter, sep="\n")
    Fe2 Fe1 (0, 0, 0)
    ((0, 0, 0),)
    (1, 0)
    [[2. 0. 0.]
     [0. 2. 0.]
     [0. 0. 2.]]

Or any other property of the convention.

.. hint::

    The main principle of changing the convention can be formulated as "Physical
    properties of the Hamiltonian do not depend on its convention".


.. _user-guide_usage_spin-hamiltonian_units:

Units
=====

:py:class:`.SpinHamiltonian` supports a number of units for its parameters. See
:ref:`user-guide_usage_units_parameters` for the full list.

First, we prepare an empty Hamiltonian for the examples in this section

.. doctest::

    >>> spinham = spinham.get_empty()
    >>> spinham.convention = spinham.convention.get_modified(multiple_counting=True)

When you use :py:meth:`.SpinHamiltonian.add` method the vector/matrix/tensor that you pass
as the ``parameter`` argument is interpreted in the units of
:py:attr:`.SpinHamiltonian.units`. By default, those are milli-electronvolts (``meV``).

.. doctest::

    >>> spinham.units
    'meV'

There is a number of ways to control the units of the parameters. Let us add a single
exchange interaction to the Hamiltonian to illustrate those.

.. doctest::

    >>> spinham.add(nus=[(0,0,0)], alphas=[0,1], parameter = np.eye(3))

Now the interaction parameter is an isotropic exchange with the value of 1 meV.

.. doctest::

    >>> for _, _, parameter in spinham.p22:
    ...     print(parameter)
    [[1. 0. 0.]
     [0. 1. 0.]
     [0. 0. 1.]]

You can change the units at the hamiltonian level as

.. doctest::

    >>> spinham.units = 'Joule'
    >>> spinham.units
    'Joule'

All parameters of the Hamiltonian are automatically converted to the new units

.. doctest::

    >>> for _, _, parameter in spinham.p22:
    ...     print(parameter)
    [[1.60217663e-22 0.00000000e+00 0.00000000e+00]
     [0.00000000e+00 1.60217663e-22 0.00000000e+00]
     [0.00000000e+00 0.00000000e+00 1.60217663e-22]]

and all parameters that will be added later on are expected to be given in the units
of Joule.

Alternatively, you can specify the units of the parameter while adding it to the
Hamiltonian as

.. doctest::

    >>> spinham.add(
    ...     nus=[(0,0,0)],
    ...     alphas=[1,0],
    ...     parameter=np.eye(3),
    ...     units='meV'
    ... )

Then the parameter is automatically converted to the units of the Hamiltonian

.. doctest::

    >>> for _, _, parameter in spinham.p22:
    ...     print(parameter)
    [[1.60217663e-22 0.00000000e+00 0.00000000e+00]
     [0.00000000e+00 1.60217663e-22 0.00000000e+00]
     [0.00000000e+00 0.00000000e+00 1.60217663e-22]]
    [[1.60217663e-22 0.00000000e+00 0.00000000e+00]
     [0.00000000e+00 1.60217663e-22 0.00000000e+00]
     [0.00000000e+00 0.00000000e+00 1.60217663e-22]]



.. _user-guide_usage_spin-hamiltonian_magnetic-vs-non-magnetic:

Magnetic vs non-magnetic atoms
==============================

The Hamiltonian can contain any amount of atoms (:py:attr:`.SpinHamiltonian.atoms`).
However, not all of them necessarily have parameters associated with them.

In magnopy we classify atoms into two types

*   Magnetic atoms

    Those are the atoms that have at least one parameter of the spin Hamiltonian associated
    with them. Even if all elements of the corresponding vector/matrix/tensor of the
    parameter are zeros.

*   Non-magnetic atoms

    All other atoms.

.. note::

    Concept of magnetic and non-magnetic atoms **is not related** to the physical
    properties of the atom (such as spin, g-factor, etc.). It only appears in the context
    of the spin Hamiltonian.

Now, magnetic atoms are the ones that contribute to the physics of the spin Hamiltonian.
In fact in the |paper-2026|_ only magnetic atoms are considered.

However, non-magnetic atoms together with the magnetic ones are typically used
to define the crystal structure and relevant high-symmetry points and paths in the
reciprocal space. Thus, Magnopy keeps track of all atoms and of the magnetic atoms.

*   :py:attr:`.SpinHamiltonian.atoms` contains all atoms of the crystal: both magnetic and
    non-magnetic.

*   :py:attr:`.SpinHamiltonian.magnetic_atoms` contains only magnetic atoms.

From the perspective of the user the following is important:

*   Atom's indices in the context of interaction parameters of
    :py:class:`.SpinHamiltonian` are always the indices of
    :py:attr:`.SpinHamiltonian.atoms`.

*   ``spin_directions`` are always given in the order of
    :py:attr:`.SpinHamiltonian.magnetic_atoms`.

*   Named interactions (i. e. Zeeman interaction or magnetic dipole-dipole interaction)
    are always added only to magnetic atoms by default.

*   The order of atoms in :py:attr:`.SpinHamiltonian.magnetic_atoms` is the same as the
    order of atoms in :py:attr:`.SpinHamiltonian.atoms`.

Let us give an example to illustrate this concept.

.. doctest::

    >>> spinham = spinham.get_empty()
    >>> spinham.atoms.names
    ['Fe1', 'Fe2']


There are no parameters in the Hamiltonian

.. doctest::

    >>> len(spinham.parameters())
    0

Therefore, there are no magnetic atoms in the Hamiltonian

.. doctest::

    >>> spinham.magnetic_atoms.names
    []

Now we add an exchange interaction between "Fe2" from the unit cell :math:`(0, 0, 0)`
and "Fe2" from the unit cell :math:`(1, 0, 0)`.

.. doctest::

    >>> spinham.add(nus=[(0,0,0)], alphas=[1,1], parameter = np.eye(3))

Now there is a single interaction in the Hamiltonian that involves only "Fe2" atom.

.. doctest::

    >>> len(spinham.parameters())
    1
    >>> for _, alphas, _ in spinham.parameters():
    ...     print(alphas)
    (1, 1)

Therefore, "Fe2" is the only magnetic atom in the Hamiltonian now, while "Fe1" is still
non-magnetic.

.. doctest::

    >>> spinham.magnetic_atoms.names
    ['Fe2']


.. hint::
    You can always convert the index of an atom in :py:attr:`.SpinHamiltonian.atoms` to
    the index of the same atom in :py:attr:`.SpinHamiltonian.magnetic_atoms` and vice
    versa using the properties :py:attr:`.SpinHamiltonian.map_to_magnetic` and
    :py:attr:`.SpinHamiltonian.map_to_all`.

    Recall that there are two atoms in the Hamiltonian: "Fe1" and "Fe2".
    "Fe1" is non-magnetic and "Fe2" is magnetic.

    .. doctest::

        >>> spinham.atoms.names
        ['Fe1', 'Fe2']
        >>> spinham.magnetic_atoms.names
        ['Fe2']

    The index of the "Fe2" atom in :py:attr:`.SpinHamiltonian.atoms` is 1, while its index
    in :py:attr:`.SpinHamiltonian.magnetic_atoms` is 0, in other words

    .. doctest::

        >>> index_in_atoms = 1
        >>> index_in_magnetic_atoms = 0

    To convert from one to another use

    .. doctest::

        >>> spinham.map_to_magnetic[index_in_atoms]
        0

    and

    .. doctest::

        >>> spinham.map_to_all[index_in_magnetic_atoms]
        1

    :py:attr:`.SpinHamiltonian.map_to_magnetic` and :py:attr:`.SpinHamiltonian.map_to_all`
    are simply lists of indices

    .. doctest::

        >>> spinham.map_to_magnetic
        [None, 0]
        >>> spinham.map_to_all
        [1]

    Note that you cannot convert the index of the "Fe1" atom in
    :py:attr:`.SpinHamiltonian.atoms` to the index in
    :py:attr:`.SpinHamiltonian.magnetic_atoms` as atom "Fe1" is non-magnetic.


.. _user-guide_usage_spin-hamiltonian_magnetic-field:

Magnetic field
==============

Zeeman interaction can be added by hand using the :py:meth:`.SpinHamiltonian.add`, as
it has the form of :ref:`ug_tb_sh_1-1`. However, we recommend using pre-defined method
that does it automatically: :py:meth:`.SpinHamiltonian.set_magnetic_field`.

Zeeman term in Magnopy is stored as part of the :py:attr:`.SpinHamiltonian.p1`
parameters. The value of the magnetic flux density :math:`\boldsymbol{B}` can be accessed
via :py:attr:`.SpinHamiltonian.magnetic_field` property.

First, we create an empty Hamiltonian

.. doctest::

    >>> spinham = spinham.get_empty()
    >>> spinham.units = "meV"

Which atoms?
------------

By default magnetic field is being added only to the
:ref:`magnetic atoms <user-guide_usage_spin-hamiltonian_magnetic-vs-non-magnetic>`.
Since there is no magnetic atoms in the Hamiltonian

.. doctest::

    >>> spinham.magnetic_atoms.names
    []

addition of the magnetic field has no effect

.. doctest::

    >>> spinham.set_magnetic_field([1, 0, 0])
    >>> len(spinham.p1)
    0
    >>> spinham.magnetic_field
    array([0., 0., 0.])

Typically, there are some interaction parameters already added to the Hamiltonian and the
magnetic atoms are clearly identified. For example, exchange interaction

.. doctest::

    >>> spinham.add(nus=[(0, 0, 0)], alphas=[0, 1], parameter = np.eye(3))
    >>> spinham.magnetic_atoms.names
    ['Fe1', 'Fe2']

If we add the magnetic field now

.. doctest::

    >>> spinham.set_magnetic_field([1, 0, 0])
    >>> len(spinham.p1)
    2
    >>> for nus, alphas, parameter in spinham.p1:
    ...     print(spinham.atoms.names[alphas[0]], parameter, sep="   ")
    Fe1   [0.11576764 0.         0.        ]
    Fe2   [0.11576764 0.         0.        ]

We see that the Zeeman term is being added to both atoms, as they are both magnetic.


Alternatively, you can control explicitly to which atoms the magnetic field is
being added by supplying the atom's indices of :py:attr:`.SpinHamiltonian.atoms`

.. doctest::

    >>> # First, reset the parameters of the Hamiltonian
    >>> spinham = spinham.get_empty()
    >>> # There is no magnetic atoms in the Hamiltonian
    >>> spinham.magnetic_atoms.names
    []

.. doctest::

    >>> spinham.set_magnetic_field([1, 0, 0], alphas=[0, 1])


We can check that the magnetic field have been added to both atoms, even though they were
non-magnetic

.. doctest::

    >>> len(spinham.p1)
    2
    >>> for nus, alphas, parameter in spinham.p1:
    ...     print(spinham.atoms.names[alphas[0]], parameter, sep="   ")
    Fe1   [0.11576764 0.         0.        ]
    Fe2   [0.11576764 0.         0.        ]

Note that now both atoms are magnetic afterwards, as they both have Zeeman interaction
associated with them.

.. doctest::

    >>> spinham.magnetic_atoms.names
    ['Fe1', 'Fe2']

Shortcut
--------

The value of the magnetic field can be checked at any time as

.. doctest::

    >>> spinham.magnetic_field
    array([1., 0., 0.])

Moreover, you can use the same property to change its value

.. doctest::

    >>> spinham.magnetic_field = [0, 2, 0]
    >>> spinham.magnetic_field
    array([0., 2., 0.])
    >>> for nus, alphas, parameter in spinham.p1:
    ...     print(spinham.atoms.names[alphas[0]], parameter, sep="   ")
    Fe1   [0.         0.23153527 0.        ]
    Fe2   [0.         0.23153527 0.        ]

The latter is equivalent to

.. doctest::

    >>> spinham.set_magnetic_field(B=[0, 2, 0])
    >>> spinham.magnetic_field
    array([0., 2., 0.])
    >>> for nus, alphas, parameter in spinham.p1:
    ...     print(spinham.atoms.names[alphas[0]], parameter, sep="   ")
    Fe1   [0.         0.23153527 0.        ]
    Fe2   [0.         0.23153527 0.        ]

Note that the method :py:meth:`.SpinHamiltonian.set_magnetic_field` is more powerful as it
allows controlling with which atoms the magnetic field is interacting, while the property
:py:attr:`.SpinHamiltonian.magnetic_field` always adds Zeeman term only for the magnetic
atoms.

Zeeman parameters
-----------------

Zeeman interaction is store as part of the :ref:`ug_tb_sh_1-1` of the Hamiltonian.
The combined effect of the Zeeman term and any other :ref:`ug_tb_sh_1-1` can always be
checked by looking at :py:attr:`.SpinHamiltonian.p1` parameters.

However, if you would like to check the Zeeman parameters by themselves, you can use the
property :py:attr:`.SpinHamiltonian.zeeman_parameters`.

Let us illustrate with an example. First, create an empty Hamiltonian

.. doctest::

    >>> spinham = spinham.get_empty()

Then add some non-Zeeman :ref:`ug_tb_sh_1-1` to the Hamiltonian

.. doctest::

    >>> spinham.add(nus=[], alphas=[0], parameter = [-1, 0, 2])
    >>> spinham.add(nus=[], alphas=[1], parameter = [0, 0, 7])

Now add the Zeeman term

.. doctest::

    >>> spinham.magnetic_field = [1, 0, 0]

If we check the :py:attr:`.SpinHamiltonian.p1` parameters now, we see that the
Zeeman term is being added to the existing :ref:`ug_tb_sh_1-1`

.. doctest::

    >>> for nus, alphas, parameter in spinham.p1:
    ...     print(spinham.atoms.names[alphas[0]], parameter, sep="   ")
    Fe1   [-0.88423236  0.          2.        ]
    Fe2   [0.11576764 0.         7.        ]

However, one can access the Zeeman parameters separately as

.. doctest::

    >>> for nus, alphas, parameter in spinham.zeeman_parameters:
    ...     print(spinham.atoms.names[alphas[0]], parameter, sep="   ")
    Fe1   [0.11576764 0.         0.        ]
    Fe2   [0.11576764 0.         0.        ]

Incremental change
------------------

The method :py:meth:`.SpinHamiltonian.set_magnetic_field` and the property
:py:attr:`.SpinHamiltonian.magnetic_field` set the value of the magnetic field. However,
if you'd like to change the value of the magnetic field incrementally, you can use the
method :py:meth:`.SpinHamiltonian.add_magnetic_field`.

First, get an empty Hamiltonian

.. doctest::

    >>> spinham = spinham.get_empty()

Then, increase the magnetic field from 0 to 1 Tesla in steps of 0.1 Tesla

.. doctest::

    >>> for i in range(10):
    ...     if i == 0: print("-"*40)
    ...
    ...     spinham.add_magnetic_field(B = [0.1, 0, 0], alphas=[0, 1])
    ...
    ...     print(f"B = {np.round(spinham.magnetic_field, decimals=1)} Tesla")
    ...     for nus, alphas, parameter in spinham.p1:
    ...         print(spinham.atoms.names[alphas[0]], parameter, sep="   ")
    ...     print("-"*40)
    ----------------------------------------
    B = [0.1 0.  0. ] Tesla
    Fe1   [0.01157676 0.         0.        ]
    Fe2   [0.01157676 0.         0.        ]
    ----------------------------------------
    B = [0.2 0.  0. ] Tesla
    Fe1   [0.02315353 0.         0.        ]
    Fe2   [0.02315353 0.         0.        ]
    ----------------------------------------
    B = [0.3 0.  0. ] Tesla
    Fe1   [0.03473029 0.         0.        ]
    Fe2   [0.03473029 0.         0.        ]
    ----------------------------------------
    B = [0.4 0.  0. ] Tesla
    Fe1   [0.04630705 0.         0.        ]
    Fe2   [0.04630705 0.         0.        ]
    ----------------------------------------
    B = [0.5 0.  0. ] Tesla
    Fe1   [0.05788382 0.         0.        ]
    Fe2   [0.05788382 0.         0.        ]
    ----------------------------------------
    B = [0.6 0.  0. ] Tesla
    Fe1   [0.06946058 0.         0.        ]
    Fe2   [0.06946058 0.         0.        ]
    ----------------------------------------
    B = [0.7 0.  0. ] Tesla
    Fe1   [0.08103735 0.         0.        ]
    Fe2   [0.08103735 0.         0.        ]
    ----------------------------------------
    B = [0.8 0.  0. ] Tesla
    Fe1   [0.09261411 0.         0.        ]
    Fe2   [0.09261411 0.         0.        ]
    ----------------------------------------
    B = [0.9 0.  0. ] Tesla
    Fe1   [0.10419087 0.         0.        ]
    Fe2   [0.10419087 0.         0.        ]
    ----------------------------------------
    B = [1. 0. 0.] Tesla
    Fe1   [0.11576764 0.         0.        ]
    Fe2   [0.11576764 0.         0.        ]
    ----------------------------------------

The same can be achieved with :py:attr:`.SpinHamiltonian.magnetic_field` or
:py:meth:`.SpinHamiltonian.set_magnetic_field`, see the dropdown below.

.. dropdown:: Alternative styles

    .. doctest::

        >>> spinham = spinham.get_empty()

    .. doctest::

        >>> fields = np.linspace([0.1, 0, 0], [1, 0, 0], 10)
        >>> for field in fields:
        ...
        ...     spinham.set_magnetic_field(B = field, alphas=[0, 1])
        ...
        ...     print(f"B = {np.round(spinham.magnetic_field, decimals=1)} Tesla")
        ...     for nus, alphas, parameter in spinham.p1:
        ...         print(spinham.atoms.names[alphas[0]], parameter, sep="   ")
        ...     print("-"*40)
        B = [0.1 0.  0. ] Tesla
        Fe1   [0.01157676 0.         0.        ]
        Fe2   [0.01157676 0.         0.        ]
        ----------------------------------------
        B = [0.2 0.  0. ] Tesla
        Fe1   [0.02315353 0.         0.        ]
        Fe2   [0.02315353 0.         0.        ]
        ----------------------------------------
        B = [0.3 0.  0. ] Tesla
        Fe1   [0.03473029 0.         0.        ]
        Fe2   [0.03473029 0.         0.        ]
        ----------------------------------------
        B = [0.4 0.  0. ] Tesla
        Fe1   [0.04630705 0.         0.        ]
        Fe2   [0.04630705 0.         0.        ]
        ----------------------------------------
        B = [0.5 0.  0. ] Tesla
        Fe1   [0.05788382 0.         0.        ]
        Fe2   [0.05788382 0.         0.        ]
        ----------------------------------------
        B = [0.6 0.  0. ] Tesla
        Fe1   [0.06946058 0.         0.        ]
        Fe2   [0.06946058 0.         0.        ]
        ----------------------------------------
        B = [0.7 0.  0. ] Tesla
        Fe1   [0.08103735 0.         0.        ]
        Fe2   [0.08103735 0.         0.        ]
        ----------------------------------------
        B = [0.8 0.  0. ] Tesla
        Fe1   [0.09261411 0.         0.        ]
        Fe2   [0.09261411 0.         0.        ]
        ----------------------------------------
        B = [0.9 0.  0. ] Tesla
        Fe1   [0.10419087 0.         0.        ]
        Fe2   [0.10419087 0.         0.        ]
        ----------------------------------------
        B = [1. 0. 0.] Tesla
        Fe1   [0.11576764 0.         0.        ]
        Fe2   [0.11576764 0.         0.        ]
        ----------------------------------------

    .. note::

        The following two examples only work because both atoms are already
        magnetic. If they are not, then one can promote them to be magnetic
        with, for example,

        .. doctest::

            >>> spinham.set_magnetic_field(B = [0, 0, 0], alphas=[0, 1])


    .. doctest::

        >>> for i in range(10):
        ...     if i == 0: print("-"*40)
        ...
        ...     spinham.magnetic_field += [0.1, 0, 0]
        ...
        ...     print(f"B = {np.round(spinham.magnetic_field, decimals=1)} Tesla")
        ...     for nus, alphas, parameter in spinham.p1:
        ...         print(spinham.atoms.names[alphas[0]], parameter, sep="   ")
        ...     print("-"*40)
        ----------------------------------------
        B = [0.1 0.  0. ] Tesla
        Fe1   [0.01157676 0.         0.        ]
        Fe2   [0.01157676 0.         0.        ]
        ----------------------------------------
        B = [0.2 0.  0. ] Tesla
        Fe1   [0.02315353 0.         0.        ]
        Fe2   [0.02315353 0.         0.        ]
        ----------------------------------------
        B = [0.3 0.  0. ] Tesla
        Fe1   [0.03473029 0.         0.        ]
        Fe2   [0.03473029 0.         0.        ]
        ----------------------------------------
        B = [0.4 0.  0. ] Tesla
        Fe1   [0.04630705 0.         0.        ]
        Fe2   [0.04630705 0.         0.        ]
        ----------------------------------------
        B = [0.5 0.  0. ] Tesla
        Fe1   [0.05788382 0.         0.        ]
        Fe2   [0.05788382 0.         0.        ]
        ----------------------------------------
        B = [0.6 0.  0. ] Tesla
        Fe1   [0.06946058 0.         0.        ]
        Fe2   [0.06946058 0.         0.        ]
        ----------------------------------------
        B = [0.7 0.  0. ] Tesla
        Fe1   [0.08103735 0.         0.        ]
        Fe2   [0.08103735 0.         0.        ]
        ----------------------------------------
        B = [0.8 0.  0. ] Tesla
        Fe1   [0.09261411 0.         0.        ]
        Fe2   [0.09261411 0.         0.        ]
        ----------------------------------------
        B = [0.9 0.  0. ] Tesla
        Fe1   [0.10419087 0.         0.        ]
        Fe2   [0.10419087 0.         0.        ]
        ----------------------------------------
        B = [1. 0. 0.] Tesla
        Fe1   [0.11576764 0.         0.        ]
        Fe2   [0.11576764 0.         0.        ]
        ----------------------------------------



    .. doctest::

        >>> fields = np.linspace([0.1, 0, 0], [1, 0, 0], 10)
        >>> for field in fields:
        ...
        ...     spinham.magnetic_field = field
        ...
        ...     print(f"B = {np.round(spinham.magnetic_field, decimals=1)} Tesla")
        ...     for nus, alphas, parameter in spinham.p1:
        ...         print(spinham.atoms.names[alphas[0]], parameter, sep="   ")
        ...     print("-"*40)
        B = [0.1 0.  0. ] Tesla
        Fe1   [0.01157676 0.         0.        ]
        Fe2   [0.01157676 0.         0.        ]
        ----------------------------------------
        B = [0.2 0.  0. ] Tesla
        Fe1   [0.02315353 0.         0.        ]
        Fe2   [0.02315353 0.         0.        ]
        ----------------------------------------
        B = [0.3 0.  0. ] Tesla
        Fe1   [0.03473029 0.         0.        ]
        Fe2   [0.03473029 0.         0.        ]
        ----------------------------------------
        B = [0.4 0.  0. ] Tesla
        Fe1   [0.04630705 0.         0.        ]
        Fe2   [0.04630705 0.         0.        ]
        ----------------------------------------
        B = [0.5 0.  0. ] Tesla
        Fe1   [0.05788382 0.         0.        ]
        Fe2   [0.05788382 0.         0.        ]
        ----------------------------------------
        B = [0.6 0.  0. ] Tesla
        Fe1   [0.06946058 0.         0.        ]
        Fe2   [0.06946058 0.         0.        ]
        ----------------------------------------
        B = [0.7 0.  0. ] Tesla
        Fe1   [0.08103735 0.         0.        ]
        Fe2   [0.08103735 0.         0.        ]
        ----------------------------------------
        B = [0.8 0.  0. ] Tesla
        Fe1   [0.09261411 0.         0.        ]
        Fe2   [0.09261411 0.         0.        ]
        ----------------------------------------
        B = [0.9 0.  0. ] Tesla
        Fe1   [0.10419087 0.         0.        ]
        Fe2   [0.10419087 0.         0.        ]
        ----------------------------------------
        B = [1. 0. 0.] Tesla
        Fe1   [0.11576764 0.         0.        ]
        Fe2   [0.11576764 0.         0.        ]
        ----------------------------------------
