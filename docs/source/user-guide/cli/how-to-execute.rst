.. _user-guide_cli_common-notes:

***********************
How to execute scripts?
***********************

On this page we explain how the command-line interface of Magnopy works.

A number of command-line scripts are defined in Magnopy. The name of every
single one of them starts with ``magnopy-``. Examples on this page use
``magnopy-scenario`` as a placeholder for the script's name.

.. _user-guide_cli_common-notes_help:

Getting help
============

To display the help message and check what parameters are available for each script use

.. code-block:: bash

    magnopy-scenario --help

See :ref:`user-guide_cli_common-notes_read-help` to read more about help messages.

.. _user-guide_cli_common-notes_arguments:

Arguments
=========

Every script expects some "arguments" (or "parameters" or "options") as an input. There
are three types of arguments.


.. _user-guide_cli_common-notes_arguments_positional:

Positional arguments
--------------------

Only the *value* of an argument is expected. Script recognizes its meaning based on the
position of that value within the command.

For example, assume that a script expects two positional arguments:

- first argument with the name of the input file;
- second argument with the name of the output file.

If you execute

.. code-block:: bash

    magnopy-scenario input.txt output.txt

a script will use the file "input.txt" as an input file and write output to the
"output.txt" file. However, if you execute

.. code-block:: bash

    magnopy-scenario output.txt input.txt

then "output.txt" will be used as an input file and output will be written to the
"input.txt" file.

.. important::

    Order of positional arguments matters.

.. _user-guide_cli_common-notes_arguments_keyword_with_value:

Keyword arguments with value
----------------------------

Both the *keyword* and the *value* (one or several) of an argument are expected. The
keyword is fixed and indicates the meaning of the argument. The value is the data that are
passed to the script. Keyword arguments are given *after* positional ones. Order of
keyword arguments does not matter.

For example, assume that a script expects an argument with the keyword ``--input-file``,
that requires a single value. If you execute

.. code-block:: bash

    magnopy-scenario --input-file input.txt

then "--input-file" is a keyword and "input.txt" is the value.

Next, assume that a script expects an argument with the keyword ``--magnetic-field``, that
requires three values. If you execute

.. code-block:: bash

    magnopy-scenario --magnetic-field 0 0 1

then "--magnetic-field" is a keyword and "0", "0" and "1" are the values. In this example
the values describe components of the magnetic flux density :math:`B = (0, 0, 1)`.

.. important::

    The keywords always start with "-" or "--".

.. _user-guide_cli_common-notes_arguments_keyword_without_value:

Keyword arguments without value
-------------------------------

Only the *keyword* is expected. The keyword is fixed and indicates both the meaning of the
argument and its value. This type of arguments is typically used when the argument has
only two possible values: ``True`` or ``False``. When such an argument is given to
the script it switches the default value to its opposite.

For example, assume that a script expects an argument with the keyword ``--relative`` and
default value ``False``. If you execute

.. code-block:: bash

    magnopy-scenario

a script will use ``False`` as a value for the argument with the keyword ``--relative``.
However, if you execute

.. code-block:: bash

    magnopy-scenario --relative

then a script will use ``True`` as a value for the argument with the keyword
``--relative``.

.. important::

    The keywords always start with "-" or "--".


.. _user-guide_cli_common-notes_long_vs_short:

Long vs short keywords
======================

The majority of arguments in Magnopy's scripts have two equivalent keywords: a
long one and a short one. You are free to use either of them. The long version
of the keyword starts with ``--`` and the short version of the keyword starts
with a single ``-``.

The purpose of having both long and short keywords is to provide descriptive keywords
(i.e. "long" ones) and to allow experienced users an option of using the short ones.

For example, assume that a set of arguments is defined for a script

==================== ============= ==============
Long keyword         Short keyword Value
==================== ============= ==============
``--input-file``     ``-if``       One string
``--output-file``    ``-of``       One string
``--magnetic-field`` ``-mf``       Three numbers
``--relative``       ``-r``        No value
==================== ============= ==============

Then, two following commands are equivalent

.. code-block:: bash

    magnopy-scenario --input-file input.txt --output-file output.txt --magnetic-field 0 0 1 --relative

.. code-block:: bash

    magnopy-scenario -if input.txt -of output.txt -mf 0 0 1 -r

The first one is descriptive and the second one is compact.

.. hint::

    You can use long keywords for some of the arguments and short ones for the other. The
    arguments are independent in this context.


.. _user-guide_cli_common-notes_read-help:

How to read help message
========================

Assume that you executed a command

.. code-block:: bash

    magnopy-scenario --help

and the script printed in the terminal the following

.. code-block:: text
    :linenos:

    usage: magnopy-scenario [-h] -if FILENAME -of FILENAME [-mf h_x h_y h_z] [-r]
                            [-sv [S1 ...]]

                  #   #
                  #####
                  # # #
                  # # #
                  #####
                   ###
                   ###   #   #  ###   ###  #   #  ###  ####  #   #
                   ####  ## ## #   # #     ##  # #   # #   # ## ##
                   ####  # # # ##### #  ## # # # #   # ####   ###
                   ####  #   # #   # #   # #  ## #   # #       #
                   ####  #   # #   #  ###  #   #  ###  #       #
                   #############################################
                   ##                                         ##
                   ##        Version: major.minor.micro       ##
                   ##       Release date: DD Month YYYY       ##
                   ##                                         ##
                   ##           License: GNU GPLv3            ##
                   ##       Documentation: magnopy.org        ##
                   ##     Copyright (C) 2023 Magnopy Team     ##
                   ##                                         ##
                   #############################################

    This script is doing a thing, when provided a thing.

    options:
      -h, --help            Show this help message and exit
      -if, --input-file     FILENAME
                            Input file for the script.
      -of, --output-file     FILENAME
                            Output file for the script.
      -mf, --magnetic-field H_X H_Y H_Z
                            Vector of external magnetic field, given in the units
                            of Tesla.
      -r, --relative
                            Whether to consider a thing to be a relative thing.
      -sv, --spin-values [S1 ...]
                            Spin values for the input thing.
      -om, --optimization-mode {memory,speed}
                            What kind of optimization shall be used.

A number of things can be deduced from this message.

Lines 1-2
---------

A draft of the command for using the script. Arguments that are enclosed in "[]" are
optional, other arguments are mandatory. In this example ``--input-file`` and
``--output-file`` are mandatory and all other arguments are optional.

Lines 4-24
----------

Magnopy's logo and metadata.

*   Version of Magnopy that is installed. For example "0.3.0", which would mean "0" major
    version, "3" minor version and "0" micro version;
*   Release date;
*   License;
*   Link to the website with documentation;
*   Copyright message;

Line 26
-------

Short description of what this script can do.

Lines 29-42
-----------

Full list of all supported arguments and their description.

=========== ==============================================================================
=========== ==============================================================================
lines 30-31 *   "-if" is a short keyword of the argument.
            *   "--input-file" is a long keyword of the argument.
            *   "FILENAME" is a placeholder for its value. Substitute "FILENAME" by an
                actual value.
            *   "Input file for the script." is a description of what this argument means.
lines 34-36 *   "-mf" is a short keyword of the argument.
            *   "--magnetic-field" is a long keyword of the argument.
            *   "H_X H_Y H_Z" are the placeholders for its values. Three placeholders
                indicate that this argument expects three values.
            *   "Vector of external magnetic field, given in the units of Tesla." is a
                description of what this argument means.
lines 37-38 *   "-r" is a short keyword of the argument.
            *   "--relative" is a long keyword of the argument.
            *   There are no placeholder for the value, which means that this is a
                :ref:`user-guide_cli_common-notes_arguments_keyword_without_value`.
            *   "Whether to consider a thing to be a relative thing." is a description of
                what this argument means.
lines 39-40 *   "-sv" is a short keyword of the argument.
            *   "--spin-values" is a long keyword of the argument.
            *   "[S1 ...]" is a placeholder for the values. Brackets and "..." indicate
                that this argument expects several values. For example, substitute
                "[S1 ...]" by "1 0.5 1.5" to pass three values to this argument.
            *   "Spin values for the input thing." is a description of what this argument
                means.
lines 41-42 *   "-om" is a short keyword of the argument.
            *   "--optimization-mode" is a long keyword of the argument.
            *   "{memory,speed}" is a placeholder for the values. Curly braces
                indicate that one of the pre-defined values is expected. Use either
                ``--optimization-mode memory`` or ``--optimization-mode speed``.
            *   "What kind of optimization shall be used." is a description of what this
                argument means.
=========== ==============================================================================
