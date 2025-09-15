.. _user-guide_input_model-file_on-site:

On-site anisotropy parameters
=============================

.. literalinclude:: full-hdf5.inc
  :linenos:
  :lines: 1,2,42-50
  :language: text

.. dropdown:: Description

  All names of datasets and groups are lowercase.

  Group's name: "on-site".

  * Attributes

    * "units" - one of the following: "milli-electron-Volt", "electron-Volt", "Joule", "Kelvin", "Rydberg". Case-insensitive. Full word is
      checked.

  * Groups

    * <number from 1 to n>

      Every group describes one on-site parameter.
      There may be no groups.
      Inside each group two datasets are required:

      * "atom" - a string.
      * "matrix" - :math:`3\times3` matrix of numbers.
