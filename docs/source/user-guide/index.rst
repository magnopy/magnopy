.. _user-guide:

********************
Magnopy's user guide
********************


.. toctree::
    :caption: Installation
    :maxdepth: 1
    :hidden:

    installation

.. toctree::
    :caption: Theory behind
    :maxdepth: 1
    :hidden:

    theory-behind/spin-hamiltonian
    theory-behind/convention
    theory-behind/equivalent-parameters
    theory-behind/energy-minimization

.. toctree::
    :caption: Magnopy as a Python library
    :maxdepth: 1
    :hidden:

    usage/overview
    usage/cell
    usage/atoms
    usage/convention
    usage/units
    usage/spin-hamiltonian
    usage/spin-directions
    usage/energy
    usage/lswt

.. toctree::
    :caption: Magnopy as a black box
    :maxdepth: 1
    :hidden:

    cli/how-to-execute
    cli/magnopy/index
    cli/magnopy-optimize-sd/index
    cli/magnopy-lswt/index

.. toctree::
    :caption: How-to guides
    :maxdepth: 1
    :hidden:

    how-to/index


Installation
============

Magnopy is a Python library that is distributed via |PYPI|_ and can be installed
as any other Python library.

:doc:`installation`
    Details on how to install Magnopy with ``pip`` or from source code.

Theory behind
=============

As with any tool, it is important to understand how Magnopy operates and what it can do.
In addition to the |paper-2026|_ that describes the method in detail, we prepared a few
pages with brief summaries of the key concepts needed to understand the code.

:doc:`theory-behind/spin-hamiltonian`
    What is the spin Hamiltonian and how is it defined in Magnopy.

:doc:`theory-behind/convention`
    What is a convention of the spin Hamiltonian and how is it treated in Magnopy.

:doc:`theory-behind/energy-minimization`
    Technical details of the energy minimization procedure implemented in Magnopy.

.. _user-guide_magnopy-as-a-python-library:

Magnopy as a Python library
===========================

The first way to use Magnopy is to use it within Python scripts.

For deeper understanding of Magnopy you can read the overview that graphically
summarizes data structures of Magnopy. Then you can read materials below where
we explain each concept and give code examples.

:doc:`usage/overview`
    An overview of the core objects in Magnopy. A good place to start.

:doc:`usage/cell`
    Cell that defines a lattice. Same as
    :math:`(\boldsymbol{a}_1, \boldsymbol{a}_2, \boldsymbol{a}_3)` in the |paper-2026|_.

:doc:`usage/atoms`
    Magnetic sites that are placed in the cell. They are called "atoms" due to the
    historical reasons, however in the |paper-2026|_ we refer to them as "magnetic sites".

:doc:`usage/convention`
    Convention of the spin Hamiltonian. Also discussed in
    :ref:`user-guide_theory-behind_convention`.

:doc:`usage/units`
    Physical units of the Hamiltonian's interaction parameters and other quantities.

:doc:`usage/spin-hamiltonian`
    Spin Hamiltonian that is defined based on the cell, atoms and convention.

:doc:`usage/spin-directions`
    Spin directions that define the state of the system (in particular, the vacuum state).
    Same as :math:`\boldsymbol{z}_{\alpha}` in the |paper-2026|_.

:doc:`usage/energy`
    Energy of the Hamiltonian in a particular state defined by the spin directions.

:doc:`usage/lswt`
    Main entry point for the linear spin wave theory calculations, defined on some
    spin Hamiltonian and some spin directions.



Magnopy as a black box
======================

The second way to use Magnopy is via a command line interface.

There are a number of scripts that take some files as input and produce some other
files as an output. You can use those scripts without any knowledge of python, but you
will need to run a command in a terminal.

:doc:`cli/how-to-execute`
    A guide on how to execute a script from the command line.

Below is the list of the scripts that are available in Magnopy today.

:doc:`cli/magnopy/index`
    Interface for stats about the package and tests of the package.

:doc:`cli/magnopy-optimize-sd/index`
    A script that performs the minimization of classical energy with respect to spin directions.

:doc:`cli/magnopy-lswt/index`
    A script that computes all terms of the magnon Hamiltonian at the level of linear
    spin wave theory.

Please reach out to the developer of Magnopy if there is a common calculation that can be
performed with it, but there is no dedicated script.


How-to guides
=============

How to guides are a collection of small code examples that demonstrate how to achieve
a small task with Magnopy. They are less verbose than the tutorials (|magnopy-tutorials|_)
and assume that you can read the materials from the
:ref:`user-guide_magnopy-as-a-python-library` and :ref:`api` for details.

:doc:`how-to/index`
    A list of how-to guides that are available in Magnopy today.

If you do not see the how-to guide that answers your question, please do not hesitate to
ask it via our :ref:`support` channels.
