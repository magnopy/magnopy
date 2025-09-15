.. _user-guide_input_model-file_atoms:

Atoms
=====

.. literalinclude:: full-hdf5.inc
  :linenos:
  :lines: 1,2,7-24,50
  :language: text

.. dropdown:: Description

  All names of datasets and groups are lowercase.

  Group's name: "atoms".

  * Groups

    At least one group is required. Every group describes an atom.
    Name of the first group is "1", the second "2", etc. Name of the last
    group is "n", where "n" is the number of atoms.

    Inside each group following datasets are recognized:

    * "name" - a string
    * "position" - three numbers (relative position)
    * "spin" - three numbers
    * "g-factor" - one number
    * "charge" - one number
    * "orbital-moment" - three numbers
    * "total-moment" - three numbers
