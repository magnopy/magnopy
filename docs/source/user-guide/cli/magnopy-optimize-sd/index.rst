.. _user-guide_cli_optimize-sd:

*******************
magnopy-optimize-sd
*******************

This scenario optimizes classical energy of the spin Hamiltonian and finds the spin
directions that describe a local minima of the energy landscape.

First we give an example of what this script can compute and then go through its usage
and parameters in details.

.. _user-guide_cli_optimize-sd_help:

Help message
============
This page of documentation is written by hand and might become outdated due to the
human error. Moreover, we do not intend to cover all possible parameters of the script
in this page of the documentation. To get the automatically generated description of **all**
input parameters, that is produced by the actual version of magnopy that is installed
in you environment use

.. code-block::

    magnopy-optimize-sd --help

That should output something similar to

.. literalinclude:: help.inc
    :language: text

At the very beginning there is a syntax for the usage of the script, where required
arguments are given without brackets and optional arguments are written within brackets.
Then there is a logo of magnopy, followed by the list of supported arguments with their
short and long names and explanation of what they represent.

.. hint::

    The short (i.e. ``-sf``) and long (i.e. ``--spinham-filename``) are absolutely equivalent.
    Feel free to use either of them. The long name is usually self-explanatory and
    the short one is added purely for the convenience of the user.

.. _user-guide_cli_optimize-sd_spinham:

Spin Hamiltonian and its source
===============================

This script works with the spin Hamiltonian that is coming from some third-party software.
At the moment magnopy supports |TB2J|_ and |GROGU|_.

.. hint::
    There is number of ways to use this script with the hand-made Hamiltonian:

    * Prepare the file that mimics the format of |TB2J|_.
    * Prepare the file that mimics the |GROGU-FF|_.
    * Prepare the spin Hamiltonian programmatically and use the scenario of this
      command-line script from within your python scripts: :py:func:`.scenarios.optimize_sd`.

To tell the script what spin Hamiltonian to use provide

* Source of the spin Hamiltonian (``-ss`` or ``--spinham-source``);
* Path to the file with the spin Hamiltonian (``-sf`` or ``--spinham-filename``)

For example, if the file with the spin Hamiltonian is located in the
"data/hamiltonians/trial1/TB2J/exchange.out" and the source of the file is |TB2J|_,
then pass to the script two parameters either in the short form

.. code-block:: bash

    magnopy-optimize-sd  -ss TB2J -sf data/hamiltonians/trial1/TB2J/exchange.out ...

or in the long form

.. code-block:: bash

    magnopy-optimize-sd  -spinham-source TB2J -spinham-filename data/hamiltonians/trial1/TB2J/exchange.out ...

.. note::
    The dots ``...`` are not a part of the syntax. They are used only to highlight the
    parameters that are described in the particular chapter of the documentation and
    hide all other parameters that might or might not be passed to the script.


.. _user-guide_cli_optimize-sd_supercell:

Minimization domain
===================

By default magnopy vary only the spins within the original unit cell of the Hamiltonian.
In that way it can miss the true ground state if it spans over several unit cells that
can not be transformed into one another by a simple translation. To address this issue,
we offer an option of minimization on the supercell. The supercell is produced by a number of
translations of the original unit cell (``-s`` or ``--superell``). For example, to ask
for a minimization of the :math:`3\times7\times2` supercell one can use the command, in
the short form

.. code-block:: bash

    magnopy-optimize-sd ... -s 3 7 2 ...

or in the long form

.. code-block:: bash

    magnopy-optimize-sd ... --supercell 3 7 2 ...

In that case every spin in the :math:`3\times7\times2` supercell is treated as an
individual variable. Note, that the computational cost of minimization will grow with
the size of the supercell.


.. note::
    The dots ``...`` are not a part of the syntax. They are used only to highlight the
    parameters that are described in the particular chapter of the documentation and
    hide all other parameters that might or might not be passed to the script.


.. _user-guide_cli_optimize-sd_tolerance:

Accuracy or tolerance conditions
================================

In theory numerical optimization can continues indefinitely, improving accuracy of
some target value with each step. In reality an algorithm should reach some values of the
tolerance after the finite amount of steps.

The minimization algorithm implemented in magnopy [1]_ traces two values:

* Absolute value of the difference in total energy between two consecutive steps
  of the minimization (``-et`` or ``--energy-tolerance``).
* Maximum (among all spins of the unit cell or supercell) value of the torque vector (``-tt`` or
  ``--torque-tolerance``).

An algorithm stops and output the obtained spin directions when both tolerance parameters
are reached. The default values, that magnopy uses should lead to the reasonable results
in most of the cases.

However, if you want to increase accuracy of one of the parameters or both, then try to
pass the corresponding parameters to the script. In the short form

.. code-block:: bash

    magnopy-optimize-sd ... -et 0.000001 -tt 0.001 ...

or in the long form

.. code-block:: bash

    magnopy-optimize-sd ... --energy-tolerance 0.000001 --torque-tolerance 0.001 ...


.. note::
    The dots ``...`` are not a part of the syntax. They are used only to highlight the
    parameters that are described in the particular chapter of the documentation and
    hide all other parameters that might or might not be passed to the script.

.. _user-guide_cli_optimize-sd_magnetic-field:

External magnetic field
=======================

The file with the :ref:`spin Hamiltonian <user-guide_cli_optimize-sd_spinham>`
specifies the interaction parameters that are intrinsic to the material.

In order to add additional effects, for instance an external magnetic field one
can use the ``-mf`` or ``--magnetic-field`` parameter.

This parameter expects three numbers, that specify three Cartesian components of the
external magnetic field. The value of the provided vector is interpreted in Tesla.

For example to add magnetic field of 2.42 Tesla along the direction :math:`(1, 1, 0)`
(i.e. in the :math:`xy` plane, right in between the :math:`x` and :math:`y` axis) pass
to the script the parameter, in the short form

.. code-block:: bash

    magnopy-optimize-sd ... -mf 1.7112 1.7112 0 ...

or in the long form

.. code-block:: bash

    magnopy-optimize-sd ... --magnetic-field 1.7112 1.7112 0 ...


.. note::
    The dots ``...`` are not a part of the syntax. They are used only to highlight the
    parameters that are described in the particular chapter of the documentation and
    hide all other parameters that might or might not be passed to the script.

.. _user-guide_cli_optimize-sd_output:

Output of the script
====================
The script have two types of the output:

*   Text output to the console

    Magnopy outputs the progress of the calculation in the standard output stream, that is
    typically printed directly to the terminal. If you would like to save this text in a
    file, we recommend to use stream redirect ``>`` operator as

    .. code-block:: bash

        magnopy-optimize-sd ... > output.txt

    In that way there will be no output to the console, but all the information will be
    saved in the file "output.txt".

*   Output that is saved in the separated files.

    A number of the files will be saved in the folder that is named "magnopy-results"
    by default. If you would like to change its name, for instance to "magnopy-SO-trial-1",
    then you can use the parameter ``-of`` or ``--output-folder``. In the short form

    .. code-block:: bash

        magnopy-optimize-sd ... -of magnopy-SO-trial-1 ...

    or in the long form

    .. code-block:: bash

        magnopy-optimize-sd ... --output-folder magnopy-SO-trial-1 ...

    .. note::

        The visual capabilities of magnopy require a third-party plotting library
        |plotly|_. It is not included as a default dependency of magnopy and therefore,
        have to be installed manually. It can be installed with ``pip``, in the same
        way as magnopy:

        .. code-block:: bash

            pip install plotly

        or

        .. code-block:: bash

            pip3 install plotly


.. note::
    The dots ``...`` are not a part of the syntax. They are used only to highlight the
    parameters that are described in the particular chapter of the documentation and
    hide all other parameters that might or might not be passed to the script.


References
==========

.. [1] Ivanov, A.V., Uzdin, V.M. and Jónsson, H., 2021.
    Fast and robust algorithm for energy minimization of spin systems applied
    in an analysis of high temperature spin configurations in terms of skyrmion
    density.
    Computer Physics Communications, 260, p.107749.
