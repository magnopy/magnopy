.. toctree::
    :maxdepth: 1
    :hidden:

    User Guide <user-guide/index>
    user-support/index
    api/index
    release-notes/index
    cite
    development/index

:Release: |version|
:Date: |release_date|

**Useful links**:
|issue-tracker|_ |
:ref:`Cite <cite>` |
:ref:`support`

.. hint::

  Install magnopy with optional dependencies (|plotly|_ and |matplotlib|_) to get magnopy
  to produce graphics

  .. code-block:: bash

    pip install "magnopy[visual]"

What is magnopy?
================

Magnopy is a python code that, given :ref:`user-guide_theory-behind_spin-hamiltonian` in **any**
:ref:`convention <user-guide_theory-behind_convention-problem>`, computes bosonic Hamiltonian
of the form

.. include:: core-formulas/bosonic-hamiltonian.inc

where

* :math:`E^{(0)}` is a classical energy of the ground state;

and the next three terms are derived with the Linear Spin Wave Theory (LSWT):

* :math:`E^{(2)}` is a quantum correction to the ground state energy;
* :math:`\omega_{\alpha}(\boldsymbol{k})` is magnon dispersion at the level of LSWT;
* :math:`\sum_{\boldsymbol{k}}\Delta(\boldsymbol{k})` is an energy measure of
  :math:`\pm \boldsymbol{k}` asymmetry.


Planned to be implemented
=========================

(Please indicate to the developers if you would like to use those features or if you
would like to participate in their implementation)

* Renormalized Spin Wave theory (RSWT)
* Magnon decay rate
* ...
