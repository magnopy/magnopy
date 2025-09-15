.. _user-guide_cli_magnopy-lswt:

************
magnopy-lswt
************

This scenario runs a calculation for the given spin Hamiltonian at the level of the
linear spin wave theory and outputs majority of the results that magnopy can compute.



Help message
============
This page of documentation is written by hand and might become outdated due to the
human error. Moreover, we do not intend to cover all possible parameters of the script
in this page of the documentation. To get the automatically generated description of **all**
input parameters, that is produced by the actual version of magnopy that is installed
in you environment use

.. code-block::

    magnopy-lswt --help

That should output something similar to

.. literalinclude:: help.inc
    :language: text

At the very beginning there is a syntax for the usage of the script, where required
arguments are given without brackets and optional arguments are written within brackets.
Then there is a logo of magnopy, followed by the list of supported arguments with their
short and long names and explanation of what they represent.

.. hint::

    The short (i.e. ``-sf``) and long (i.e. ``--spinham-filename``) are absolutely equivalent.
    Feel free to use either of them. The long name usually is self-explanatory and
    the short one is added purely for the convenience of the user.


.. _user-guide_cli_lswt_spinham:

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

    magnopy-lswt  -ss TB2J -sf data/hamiltonians/trial1/TB2J/exchange.out ...

or in the long form

.. code-block:: bash

    magnopy-lswt  -spinham-source TB2J -spinham-filename data/hamiltonians/trial1/TB2J/exchange.out ...

.. note::
    The dots ``...`` are not a part of the syntax. They are used only to highlight the
    parameters that are described in the particular chapter of the documentation and
    hide all other parameters that might or might not be passed to the script.

.. _user-guide_cli_lswt_ground-state:

Ground state
============

For the calculation of exited states (what magnons are) one need to knows the
ground state - spin directions for every spin in the Hamiltonian. There are two ways for
magnopy to know the spin directions

*   Input from the user

    User can provide a file with spin directions using either the short form of the argument

    .. code-block:: bash

        magnopy-lswt ... -sd SPIN-DIRECTIONS.txt

    or long form

    .. code-block:: bash

        magnopy-lswt ... --spin-directions SPIN-DIRECTIONS.txt

*   Internal optimization

    If user do not provide any input, then magnopy tries to optimize spin directions within
    unit cell.

.. note::
    The dots ``...`` are not a part of the syntax. They are used only to highlight the
    parameters that are described in the particular chapter of the documentation and
    hide all other parameters that might or might not be passed to the script.



.. _user-guide_cli_lswt_reciprocal-space:

K-path and k-points
===================

Magnopy solves magnon problem for a set of points in reciprocal space. Therefore, it needs
to know a set of k-points to perform the calculations. User is left with two options

*   Provide explicit list of k-points, using short form of the arguments

    .. code-block:: bash

        magnopy-lswt ... -kps  K-POINTS.txt

    or long form

    .. code-block:: bash

        magnopy-lswt ... --kpoints  K-POINTS.txt

*   Let magnopy deduce the set of high-symmetry points based on the space group of the
    crystal and use recommended k-path (see documentation of |wulfric|_ for more details,
    magnopy uses ``convention="HPKOT"``). In that second case the user can still control
    the k-path, but limited to the list of the predefined high-symmetry points.

    .. code-block:: bash

        magnopy-lswt ... -kp  GAMMA-X-S|GAMMA-Y

    or in the long form

    .. code-block:: bash

        magnopy-lswt ... --k-path  GAMMA-X-S|GAMMA-Y

.. note::
    The dots ``...`` are not a part of the syntax. They are used only to highlight the
    parameters that are described in the particular chapter of the documentation and
    hide all other parameters that might or might not be passed to the script.

.. _user-guide_cli_lswt_magnetic-field:

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

    magnopy-lswt ... -mf 1.7112 1.7112 0 ...

or in the long form

.. code-block:: bash

    magnopy-lswt ... --magnetic-field 1.7112 1.7112 0 ...


.. note::
    The dots ``...`` are not a part of the syntax. They are used only to highlight the
    parameters that are described in the particular chapter of the documentation and
    hide all other parameters that might or might not be passed to the script.

.. _user-guide_cli_lswt_output:

Output of the script
====================
The script have two types of the output:

*   Text output to the console

    Magnopy outputs the progress of the calculation in the standard output stream, that is
    typically printed directly to the terminal. If you would like to save this text in a
    file, we recommend to use stream redirect ``>`` operator as

    .. code-block:: bash

        magnopy-lswt ... > output.txt

    In that way there will be no output to the console, but all the information will be
    saved in the file "output.txt".

*   Output that is saved in the separated files.

    A number of the files will be saved in the folder that is named "magnopy-results"
    by default. If you would like to change its name, for instance to "magnopy-LSWT-trial-1",
    then you can use the parameter ``-of`` or ``--output-folder``. In the short form

    .. code-block:: bash

        magnopy-lswt ... -of magnopy-LSWT-trial-1 ...

    or in the long form

    .. code-block:: bash

        magnopy-lswt ... --output-folder magnopy-LSWT-trial-1 ...

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
