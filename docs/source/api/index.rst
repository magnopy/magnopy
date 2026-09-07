.. module:: magnopy

.. _api:

*************
API reference
*************

:Release: |version|

The main interface to the package should be imported as

.. doctest::

  >>> import magnopy

Sub-modules
===========
.. toctree::
  :maxdepth: 1

  io
  scenarios
  examples
  si
  experimental


Classes
=======

.. autosummary::
    :toctree: generated/

    Convention
    SpinHamiltonian
    Energy
    LSWT
    PlotlyEngine

Functions
=========

.. autosummary::
  :toctree: generated/

  solve_via_colpa
  span_local_rf
  span_local_rfs
  logo
  multiprocess_over_k
  make_supercell
  is_eigenstate

Interaction parameters
======================

.. autosummary::
  :toctree: generated/

  converter22
  converter43
  get_equivalent_parameters

Unit tests
==========

.. autosummary::
  :toctree: generated/

  test

Exceptions
==========

.. autosummary::
  :toctree: generated/

  ConventionError
  ColpaFailed
