.. _user-guide_usage_convention:

**********
Convention
**********

For the theoretical background on the convention problem see
:ref:`user-guide_theory-behind_convention`.

Magnopy implements a compact, hashable :py:class:`.Convention` class to combine all
parameters of the spin Hamiltonian's convention.

Here is an example of the custom convention where constants :math:`C_1`, :math:`C_{2,1}`,
:math:`C_{2, 2}`, multiple counting and spin normalization are defined, while the other
constants are left undefined.

.. doctest::

    >>> from magnopy import Convention
    >>> convention = Convention(
    ...    multiple_counting=True,
    ...    spin_normalized=False,
    ...    c1=1,
    ...    c21=1,
    ...    c22=1,
    ...    name="Convention example",
    ... )

The properties of the convention can be checked at any time

.. doctest::

    >>> convention.multiple_counting
    True
    >>> convention.spin_normalized
    False
    >>> convention.c1
    1.0
    >>> convention.c21
    1.0
    >>> convention.c22
    1.0

As well as the full summary about the convention

.. doctest::

    >>> print(convention)
    "convention example" convention where
      * Bonds are counted multiple times in the sum;
      * Spin vectors are not normalized;
      * c1 = 1.0;
      * c21 = 1.0;
      * c22 = 1.0;
      * Undefined c31 factor;
      * Undefined c32 factor;
      * Undefined c33 factor;
      * Undefined c41 factor;
      * Undefined c42 factor;
      * Undefined c43 factor;
      * Undefined c44 factor;
      * Undefined c45 factor.


:py:class:`.Convention` is one of the three objects that are required for the creation of
an instance of the :py:class:`.SpinHamiltonian` class. It is stored as an attribute
:py:attr:`.SpinHamiltonian.convention` and can be accessed or modified through it.

Modifying convention
====================

Every instance of the convention class is meant to be static, therefore, the properties of
the existing instance cannot be changed for example

.. doctest::

    >>> convention.multiple_counting = False
    Traceback (most recent call last):
    ...
    AttributeError: It is intentionally forbidden to set properties of convention. Use convention.get_modified(...) and/or spinham.convention = spinham.convention.get_modified(...)

Thus Magnopy forces you to create a new instance of the convention class, even if you
need to change only one property of the convention.

For example, if you would like to change the value of :math:`C_{2,1}` from 1.0 to -1.0,
you can do the following

.. doctest::

    >>> mod_1 = Convention(
    ...     multiple_counting=convention.multiple_counting,
    ...     spin_normalized=convention.spin_normalized,
    ...     c1=convention.c1,
    ...     c21=-1.0,
    ...     c22=convention.c22,
    ...     name="Modified convention (verbose method)",
    ... )

However, this approach is unnecessarily verbose. Therefore, we provide a shortcut
(:py:meth:`.Convention.get_modified`) to create a modified convention based on the
existing one

.. doctest::

    >>> mod_2 = convention.get_modified(c21=-1.0, name="Modified convention (shortcut method)")

The resulting conventions are the same

.. doctest::

    >>> mod_1 == mod_2
    True

Comparing conventions
=====================

As you noticed, two instances of the :py:class:`.Convention` class can be compared for
equality

.. doctest::

    >>> mod_1 == mod_2
    True
    >>> mod_2 == convention
    False

We note that only the properties of the convention are compared, while the name of the
convention is not taken into account

.. doctest::

    >>> mod_1.name
    'modified convention (verbose method)'
    >>> mod_2.name
    'modified convention (shortcut method)'
    >>> # Therefore
    >>> mod_1.name == mod_2.name
    False
    >>> # But
    >>> mod_1 == mod_2
    True


Pre-defined conventions
=======================

Magnopy gives access to the predefined conventions of the spin Hamiltonian from other
codes

.. doctest::

    >>> tb2j_convention = Convention.get_predefined("TB2J")
    >>> print(tb2j_convention)
    "tb2j" convention where
      * Bonds are counted multiple times in the sum;
      * Spin vectors are normalized to 1;
      * Undefined c1 factor;
      * c21 = -1.0;
      * c22 = -1.0;
      * Undefined c31 factor;
      * Undefined c32 factor;
      * Undefined c33 factor;
      * Undefined c41 factor;
      * Undefined c42 factor;
      * Undefined c43 factor;
      * Undefined c44 factor;
      * Undefined c45 factor.

.. doctest::

    >>> grogu_convention = Convention.get_predefined("GROGU")
    >>> print(grogu_convention)
    "grogu" convention where
      * Bonds are counted multiple times in the sum;
      * Spin vectors are normalized to 1;
      * Undefined c1 factor;
      * c21 = 1.0;
      * c22 = 0.5;
      * Undefined c31 factor;
      * Undefined c32 factor;
      * Undefined c33 factor;
      * Undefined c41 factor;
      * Undefined c42 factor;
      * Undefined c43 factor;
      * Undefined c44 factor;
      * Undefined c45 factor.

To see all supported codes see :py:meth:`.Convention.get_predefined`.
