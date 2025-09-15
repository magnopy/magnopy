.. _user-guide_input_model-file_exchange:

Exchange parameters (bilinear)
==============================

.. literalinclude:: full-hdf5.inc
  :linenos:
  :lines: 1,2,30-41,50
  :language: text

.. dropdown:: Description

  All names of datasets and groups are lowercase.

  Group's name: "exchange".

  * Attributes

    * "units" - one of the following: "milli-electron-Volt", "electron-Volt", "Joule", "Kelvin", "Rydberg". Case-insensitive. Full word is
      checked.

  * Groups

    * <number from 1 to n>

      Every group describes one exchange parameter.
      There may be no groups.
      Inside each group four datasets are required:

      * "atom-1" - a string.
      * "atom-2" - a string.
      * "ijk" - three integers.
      * "matrix" - :math:`3\times3` matrix of numbers.
