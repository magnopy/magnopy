.. _user-guide_cli_optimize-sd:

*******************
magnopy-optimize-sd
*******************

This scenario optimizes classical energy of the spin Hamiltonian and finds a set of spin
directions that describes a local minimum of the energy landscape. It implements an
algorithm described in [1]_.

Visit |tutorial-optimize-sd|_ for examples of input and output files.

.. _user-guide_cli_optimize-sd_help:

Getting help
============

We recommend getting the accurate and full list of the script's parameters, that
reflects the installed version of Magnopy with the command

.. code-block::

    magnopy-optimize-sd --help

which outputs to the standard output channel (console or terminal) Magnopy's metadata and
*full* list of script's arguments. Here is an example of this output

.. hint::

    Go :ref:`here <user-guide_cli_common-notes_read-help>` to learn how to read this help
    message.

.. literalinclude:: help.inc
    :language: text

.. _user-guide_cli_optimize-sd_output-files:


Output files
============

Human-readable text with the progress of calculations and compact output data is printed
directly to the console. This output is meant to explain itself, thus we do not document
it here.

In addition, a number of .txt and/or .html files is produced.

Visit |tutorial-optimize-sd|_ for examples of the output text and files.

E_0.txt
-------

.. versionadded:: 0.4.0

A file with the value of classical energy of the optimized state. It contains a number and
unit string separated by a space symbol.

INITIAL_GUESS.txt
-----------------

.. versionadded:: 0.2.0

A file with the spin directions of the initial guess.

Magnopy makes a new random guess for the spin directions prior to each optimization (in
other words, for each execution of the script).

There are M lines in the file. M is a number of magnetic atoms in the spin Hamiltonian.
Each line has three numbers on it, separated by at least one space symbol.

First number is the x component, second the y, third the z of the spin direction vector.

SPIN_DIRECTIONS.txt
-------------------

A file with the optimized spin directions.

There are M lines in the file. M is a number of magnetic atoms in the spin Hamiltonian.
Each line has three numbers on it, separated by at least one space symbol.

First number is an x component, second - y, third - z of the spin direction vector.

SPIN_DIRECTIONS.html
--------------------

.. versionadded:: 0.2.0

**Requires** : Installation of |plotly|_ or ``magnopy[visual]``.

**Options** : Use ``--no-html`` to disable an output of this file.

An interactive .html file with 3D image of minimized spin directions and initial guess.

Data can be found in "SPIN_DIRECTIONS.txt", "INITIAL_GUESS.txt" and "SPIN_POSITIONS.txt".

SPIN_POSITIONS.txt
------------------

.. versionadded:: 0.2.0

A file with the *absolute* positions of the magnetic sites.

There are M lines in the file. M is a number of magnetic atoms in the spin Hamiltonian.
Each line has three numbers on it, separated by at least one space symbol.

First number is an x coordinate, second - y, third - z.

References
==========

.. [1] Ivanov, A.V., Uzdin, V.M. and Jónsson, H., 2021.
    Fast and robust algorithm for energy minimization of spin systems applied
    in an analysis of high temperature spin configurations in terms of skyrmion
    density.
    Computer Physics Communications, 260, p.107749.
