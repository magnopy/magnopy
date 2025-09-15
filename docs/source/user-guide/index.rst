.. _user-guide:

******************
Magnopy user guide
******************


.. toctree::
    :caption: Installation
    :maxdepth: 1

    installation

Magnopy is a python library and to get the most of it it one need to use it as one.
Nevertheless, we do not want to limit it accessability only to the ones who is familiar
with python. For that purpose we introduce a number of the scripts that perform
calculation with magnopy and can be executed as any other command line tool.


.. toctree::
    :caption: Magnopy as a black box
    :maxdepth: 1

    cli/how-to-execute
    cli/magnopy-optimize-sd/index
    cli/magnopy-lswt/index

The rest of the user guide describe how to use magnopy within python.

.. toctree::
    :caption: How-to guides
    :maxdepth: 2

    how-to/index

For the more in-deep understanding on how to use magnopy you can read the material below
where we explain how different elements of magnopy python library can be used

.. toctree::
    :caption: Fundamentals and usage
    :maxdepth: 1

    usage/data-structures
    usage/cell
    usage/atoms
    usage/convention
    usage/spin-hamiltonian
    usage/energy
    usage/lswt


With any tool it is important to understand what it does exactly. Here we prepared some
explanatory material about theory behind magnopy.

First of all, at the heart of magnopy there is a spin Hamiltonian. In the pages below we
explain how the Hamiltonian is defined and stored within magnopy.

.. toctree::
    :caption: Theory behind
    :maxdepth: 1

    theory-behind/spin-hamiltonian
    theory-behind/convention-problem
    theory-behind/multiple-counting
    theory-behind/energy-minimization
