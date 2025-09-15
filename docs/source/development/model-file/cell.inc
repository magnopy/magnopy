.. _user-guide_input_model-file_cell:

Cell
====

.. literalinclude:: full-hdf5.inc
  :linenos:
  :lines: 1-6,50
  :language: text

.. dropdown:: Description

  All names of datasets, groups and attributes are lowercase.

  Group's name: "cell".

  * Attributes

    * "units" - either "Angstrom" or "Bohr". Case-insensitive. Full word is
      checked.

  * Datasets

    * "a1" - three numbers.
    * "a2" - three numbers.
    * "a3" - three numbers.
