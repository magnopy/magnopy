.. _user-guide_usage_convention:

**********
Convention
**********

For the theoretical background of the convention problem and why it is important see
:ref:`user-guide_theory-behind_convention-problem`.

Magnopy implements a compact :py:class:`.Convention` class to store the convention of the
spin Hamiltonian. It stores all parameters that define convention in one data structure.

.. doctest::

    >>> from magnopy import Convention
    >>> convention = Convention(multiple_counting = True, spin_normalized = False, c1 = 1, c21 = 1)
    >>> convention.name
    'custom'
    >>> convention.multiple_counting
    True
    >>> convention.c21
    1.0
    >>> convention.summary()
    custom convention where
      * Bonds are counted multiple times in the sum;
      * Spin vectors are not normalized;
      * c1 = 1.0;
      * c21 = 1.0;
      * Undefined c22 factor;
      * Undefined c31 factor;
      * Undefined c32 factor;
      * Undefined c33 factor;
      * Undefined c41 factor;
      * Undefined c421 factor;
      * Undefined c422 factor;
      * Undefined c43 factor;
      * Undefined c44 factor.

There are two boolean properties (:py:attr:`.Convention.multiple_counting` and
:py:attr:`.Convention.spin_normalized`) and eleven constants that fully define the
convention of the spin Hamiltonian.

Convention is one of three objects that are required for the creation of the
:py:class:`.SpinHamiltonian` class. It is stored as an its attribute
:py:attr:`.SpinHamiltonian.convention`.

Modified convention
===================

Every instance of the convention class is meant to be static, therefore the properties of
the convention can not be changed

.. doctest::

    >>> convention.multiple_counting = False
    Traceback (most recent call last):
    ...
    AttributeError: It is intentionally forbidden to set properties of convention. Use correct methods of SpinHamiltonian class to change convention.

If you need to have a new convention, then create a new instance of the
:py:class:`.Convention` class.

For that purpose convention class has a method :py:meth:`.Convention.get_modified` that
allows you to change several properties of convention while keeping the other ones intact

.. doctest::

    >>> new_convention = convention.get_modified(c1=-1.0, c33 = 1.0)
    >>> new_convention.summary()
    custom convention where
      * Bonds are counted multiple times in the sum;
      * Spin vectors are not normalized;
      * c1 = -1.0;
      * c21 = 1.0;
      * Undefined c22 factor;
      * Undefined c31 factor;
      * Undefined c32 factor;
      * c33 = 1.0;
      * Undefined c41 factor;
      * Undefined c421 factor;
      * Undefined c422 factor;
      * Undefined c43 factor;
      * Undefined c44 factor.


Pre-defined conventions
=======================

Magnopy gives access to the predefined conventions of the spin Hamiltonian from other
popular codes

.. doctest::

    >>> tb2j_convention = Convention.get_predefined("tb2j")
    >>> vampire_convention = Convention.get_predefined("vampire")
    >>> tb2j_convention.summary()
    tb2j convention where
      * Bonds are counted multiple times in the sum;
      * Spin vectors are normalized to 1;
      * Undefined c1 factor;
      * c21 = -1.0;
      * c22 = -1.0;
      * Undefined c31 factor;
      * Undefined c32 factor;
      * Undefined c33 factor;
      * Undefined c41 factor;
      * Undefined c421 factor;
      * Undefined c422 factor;
      * Undefined c43 factor;
      * Undefined c44 factor.
    >>> vampire_convention.summary()
    vampire convention where
      * Bonds are counted multiple times in the sum;
      * Spin vectors are normalized to 1;
      * Undefined c1 factor;
      * c21 = -1.0;
      * c22 = -0.5;
      * Undefined c31 factor;
      * Undefined c32 factor;
      * Undefined c33 factor;
      * Undefined c41 factor;
      * Undefined c421 factor;
      * Undefined c422 factor;
      * Undefined c43 factor;
      * Undefined c44 factor.

To see all supported codes see :py:meth:`.Convention.get_predefined`.
