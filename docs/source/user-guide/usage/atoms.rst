.. _user-guide_usage_atoms:

***********
Atoms/Sites
***********

.. note::

    Historically, we called the collection of positions in the real space with a number
    of properties defined for each position as "atoms".

    In the |paper-2026|_ we do not use the term "atoms", but rather talk about "magnetic
    sites". However, not every atom of the crystal would be a magnetic site in general.
    Moreover, a cluster of several atoms can form a single magnetic site in some models.

    Altogether, the situation is a bit confusing.

    We continue to call a position in the real space with a number of properties (such as
    name, spin value, g-factor, etc.) an "atom" in the documentation and code of Magnopy.
    Still we want to stress that this is just a historical name. It is for the user to
    decide how to interpret an "atom" in their particular model.


Atoms are stored as a plain Python dictionary. Keys of the ``atoms`` are plural and
describe what kind of information is stored in the corresponding value. Values are lists
of the same length :math:`M^{\prime}`.


Here is an example of the ``atoms`` dictionary with six atoms in it

.. doctest::

    >>> atoms = {
    ...     "names": ["Cr1", "Br1", "S1", "Cr2", "Br2", "S2"],
    ...     "positions": [
    ...         [0.5, 0.0, 0.882382],
    ...         [0.0, 0.0, 0.677322],
    ...         [0.5, 0.5, 0.935321],
    ...         [0.0, 0.5, 0.117618],
    ...         [0.5, 0.5, 0.322678],
    ...         [0.0, 0.0, 0.064679],
    ...     ],
    ...     "spins" : [3/2, None, None, 3/2, None, None],
    ...     "g_factors" : [2.0, None, None, 2.0, None, None],
    ...     "spglib_types": [1, 2, 3, 1, 2, 3],
    ... }

.. doctest::

    >>> f"{len(atoms['names'])} atoms"
    '6 atoms'

.. doctest::

    >>> for name, position, spin, g_factor, spglib_type in zip(*atoms.values()):
    ...     print(f'Atom "{name}"')
    ...     print(f"  position is {position}")
    ...     print(f"  S = {spin}")
    ...     print(f"  g = {g_factor}")
    ...     print(f"  spglib type is {spglib_type}")
    Atom "Cr1"
      position is [0.5, 0.0, 0.882382]
      S = 1.5
      g = 2.0
      spglib type is 1
    Atom "Br1"
      position is [0.0, 0.0, 0.677322]
      S = None
      g = None
      spglib type is 2
    Atom "S1"
      position is [0.5, 0.5, 0.935321]
      S = None
      g = None
      spglib type is 3
    Atom "Cr2"
      position is [0.0, 0.5, 0.117618]
      S = 1.5
      g = 2.0
      spglib type is 1
    Atom "Br2"
      position is [0.5, 0.5, 0.322678]
      S = None
      g = None
      spglib type is 2
    Atom "S2"
      position is [0.0, 0.0, 0.064679]
      S = None
      g = None
      spglib type is 3

``atoms`` is one of the three objects that are required for the creation of an instance of
the :py:class:`.SpinHamiltonian` class. It is then stored as an immutable attribute
:py:attr:`.SpinHamiltonian.atoms`.

We specify keys of the ``atoms`` dictionary that are expected by Magnopy in the docstrings
of the relevant functions and classes.

.. hint::
    The ``atoms`` dictionary is defined in the same way as in |wulfric|_. Magnopy does not
    define any helper functions to manipulate the atoms, as it is out of the scope of
    this package. Instead we depend on |wulfric|_ for all manipulations with atoms.


Expected keys
=============


*   "names" :
    Inherited from |wulfric|_. Arbitrary labels of atoms. Each element is a string.

*   "positions" :
    **Relative** coordinates of atoms. Inherited from |wulfric|_. Each element is an
    |array-like|_ of length :math:`3`.

*   "spins" :
    Spin values for each atom. Each element is a ``float`` or ``int``.

*   "g_factors" :
    g-factors for each atom. Used in the definition of the Zeeman term. Each element is a
    ``float`` or ``int``.

*   "spglib_types" :
    Inherited from |wulfric|_. Atomic types as defined in |spglib|_. Used in the symmetry
    analysis (via |wulfric|_ and |spglib|_). See :py:func:`wulfric.get_spglib_data` and
    :py:func:`wulfric.get_spglib_types` for details. Each element is an ``int``
    (:math:`\ge 1`).

.. hint::

    The dictionary structure of ``atoms`` allows adding any other property simply by
    adding a new key and a corresponding list of values to the ``atoms`` dictionary. This
    addition would not break anything in Magnopy.

Magnetic vs non-magnetic atoms
==============================

By itself the dictionary ``atoms`` does not differentiate between magnetic and
non-magnetic atoms. The concept of magnetic and non-magnetic appears only when ``atoms``
are used to create a :py:class:`.SpinHamiltonian` object. The magnetic atoms are defined
as *atoms that have at least one interaction parameter associated with them (even if the
value of the interaction parameter is zero)*.

More details are given in the
:ref:`user-guide_usage_spin-hamiltonian_magnetic-vs-non-magnetic` page.

.. note::
    In the example above one might be tempted to say that "Cr1" and "Cr2" are magnetic
    atoms, while "Br1", "Br2", "S1", and "S2" are non-magnetic based on the values of the
    "spins" key. However, this is **not** how the magnetic atoms are defined in Magnopy.
