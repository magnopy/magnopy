.. _user-guide_cli_lswt:

************
magnopy-lswt
************

This scenario runs a calculation for the given spin Hamiltonian at the level of the
linear spin wave theory and outputs the majority of the results that Magnopy can compute.

Visit |tutorial-lswt|_ for examples of input and output files.

.. _user-guide_cli_lswt_help:

Getting help
============

We recommend getting the accurate and full list of the script's parameters that
reflect the installed version of Magnopy with the command

.. code-block::

    magnopy-lswt --help

which outputs to the standard output channel (console or terminal) Magnopy's metadata and
*full* list of script's arguments. Here is an example of this output

.. hint::

    Go :ref:`here <user-guide_cli_common-notes_read-help>` to learn how to read this help
    message.

.. literalinclude:: help.inc
    :language: text


Output files
============

Human-readable text with the progress of calculations and compact output data is printed
directly to the console. This output is meant to explain itself, thus we do not document
it here.

In addition, a number of .txt and/or .html files is produced.

Visit |tutorial-lswt|_ for examples of the output text and files.

DELTAS.png
----------

**Requires** : Installation of |matplotlib|_ or ``magnopy[visual]``.

Static image with the delta term of the magnon Hamiltonian.

Data can be found in "DELTAS.txt" and "K-POINTS.txt".

DELTAS.txt
----------

A file with the values of the delta term of the magnon Hamiltonian.

There are :math:`L + 1` lines in the file. First line is a header that indicates the
meaning of each column. Then, there are :math:`L` lines with values of magnon energies for
each of :math:`L` k-points.

Each line has one number on it. The number is a delta term of the magnon Hamiltonian.


E_0.txt
-------

.. versionadded:: 0.4.0

A file with the value of classical energy of the ground state.

There is one line in the file.

It contains a number and unit string separated by a space symbol.


E_2.txt
-------

.. versionadded:: 0.4.0

A file with the value of quantum correction to the classical energy of the ground state,
that results from linear spin wave theory.

There is one line in the file.

It contains a number and unit string separated by a space symbol.


HIGH-SYMMETRY_POINTS.txt
------------------------

.. versionadded:: 0.2.0

**Options** : Not produced if ``--kpoints`` is used.

There are :math:`N + 1` lines in the file. First line is a header, that indicates the
meaning of each column. Then, there are :math:`N` lines for :math:`N` high-symmetry
points.

Each line has one string followed by six numbers on it, separated by at least one space
symbol.

The string is a label of the high-symmetry point that can be used in the specification of
the k-path.

First three numbers are the *absolute* coordinates of the high-symmetry point in the
reciprocal space.

Last three numbers are the relative coordinates of the high-symmetry point in the basis of
the reciprocal cell of the given  unit cell (i.e. same unit cell as in the input file).

K-POINTS.html
-------------

.. versionadded:: 0.2.0

**Requires** : Installation of |plotly|_ and |scipy|_ or ``magnopy[visual]``.

**Options** : Not produced if ``--kpoints`` is used. Use ``--no-html`` to disable an
output of this file.

An interactive .html file with 3D image of the chosen k-path, high-symmetry points and
first Brillouin zones of the given unit cell (i.e. same unit cell as in the input file)
and of the primitive cell.

Part of the data can be found in "HIGH-SYMMETRY_POINTS.txt".

K-POINTS.txt
------------

.. versionadded:: 0.2.0

A file with the full list of the k-points that were used in the calculations.

There are :math:`L + 1` lines in the file. First line is a header, that indicates the
meaning of each column. Then, there are :math:`L` lines with :math:`L` k-points.

Each line has seven numbers on it, separated by at least one space symbol.

First three numbers are the *absolute* coordinates of the high-symmetry point in the
reciprocal space.

Next three numbers are the relative coordinates of the high-symmetry point in the basis of
the reciprocal cell of the given  unit cell (i. e. same unit cell as in the input file).

Last number is a single index for the k-point, that can be used for the plots (for example
band plots).

OMEGAS.png
----------

**Requires** : Installation of |matplotlib|_ or ``magnopy[visual]``.

Static image with the magnon dispersion.

Data can be found in "OMEGAS.txt" and "K-POINTS.txt".

OMEGAS.txt
----------

A file with the values of magnon energies.

There are :math:`L + 1` lines in the file. First line is a header, that indicates the
meaning of each column. Then, there are :math:`L` lines with values of magnon energies for
each of :math:`L` k-points.

Each line has :math:`M` numbers on it, separated by at least one space symbol.

Each number is a magnon energy of :math:`i`-th magnon mode.

OMEGAS-IMAG.png
---------------

**Warning**: If this file appears in the output, then something might be wrong
with the set-up of the calculations (wrong ground state, ...)

**Requires** : Installation of |matplotlib|_ or ``magnopy[visual]``.

Static image with the imaginary part of the magnon dispersion.

Data can be found in "OMEGAS-IMAG.txt" and "K-POINTS.txt".

OMEGAS-IMAG.txt
---------------

**Warning** If this file appeared, then something might be wrong with the set-up of the
calculations (wrong ground state, ...)

A file with the imaginary part of the values of magnon energies.

There are :math:`L + 1` lines in the file. First line is a header, that indicates the
meaning of each column. Then, there are :math:`L` lines with imaginary part of the values
of magnon energies for each of :math:`L` k-points.

Each line has :math:`M` numbers on it, separated by at least one space symbol.

Each number is an imaginary part of the magnon energy of :math:`i`-th magnon mode.

ONE_OPERATOR_TERMS.txt
-----------------------------

Coefficients before the one-operator terms of the magnon Hamiltonian.

There are :math:`M + 1` lines in the file. M is a number of magnetic atoms in the spin
Hamiltonian.

First line is a header, that indicates the meaning of each column. Then, there are :math:`M`
lines with values of the coefficients.

Each line has two numbers on it, separated by at least one space symbol.

First number is the real part of the coefficient, second - imaginary part.

SPIN_DIRECTIONS.html
--------------------

.. versionadded:: 0.2.0

**Requires** : Installation of |plotly|_ or ``magnopy[visual]``.

**Options** : Use ``--no-html`` to disable an output of this file.

An interactive .html file with 3D image of the spin directions that were used as the
ground state.

Part of the data can be found in "SPIN_VECTORS.txt".

SPIN_VECTORS.txt
----------------

.. versionadded:: 0.2.0

A file with the spin vectors that were used as the ground state.

There are M lines in the file. M is a number of magnetic atoms in the spin Hamiltonian.

Each line has four numbers on it, separated by at least one space symbol.

First number is the x component, second the y, third the z of the spin direction vector. Fourth
number is the spin value.
