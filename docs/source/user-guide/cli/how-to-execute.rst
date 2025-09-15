.. _user-guide_cli_common-notes:

**************************
How to execute the script?
**************************

On this page we explain how the command-line interface of magnopy works.

There is a number of calculation scenarios defined within magnopy. Each scenario
correspond to one individual script. For example :ref:`user-guide_cli_magnopy-lswt`
performs calculation at the level of linear spin-wave theory and outputs all possible
results from it.

The examples on this page use ``magnopy-scenario`` as a placeholder for the script
name.

To display the help message and check what parameters are available for each script use
one of the commands

.. code-block:: bash

    magnopy-scenario

.. code-block:: bash

    magnopy-scenario -h

.. code-block:: bash

    magnopy-scenario --help

Options
=======

Every script accepts one or more "arguments" (or "parameters" or "options") as an input.

There are three types of arguments

*   Positional arguments

    For positional arguments only the value of an argument is provided and the script
    recognizes its meaning based on the position of that value. For example assume that
    the script expects two positional arguments: first one for the input filename and
    second one for the output filename. Then if user runs the command

    .. code-block:: bash

        magnopy-scenario input_file.txt output_file.txt

    the script will use the file "input_file.txt" as an input source and "output_file.txt"
    as an output. However, if user runs the command

    .. code-block:: bash

        magnopy-scenario output_file.txt input_file.txt

    then the situation will be the opposite: "input_file.txt" will be used as an output
    and "output_file.txt" will be used as input source.

*   Keyword arguments with value

    For this type of argument user should provide a keyword and the a value of an
    argument. The keyword is used by the script to understand how to interpret the
    value. The order of the keyword arguments does not matter, but they should be given
    after the positional arguments. For example, if user runs the command

    .. code-block:: bash

        magnopy-scenario --input_file input_file.txt

    then "--input_file" is a keyword and "input_filename.txt" is the value.

    The keywords are always preceded either by "-" or "--". In magnopy every (or almost
    every) keyword argument has two available keywords: a short one (preceded with "-")
    and a long one (preceded with "--"). You can use either of them, the commands


    .. code-block:: bash

        magnopy-scenario --input_file input_file.txt

    and

    .. code-block:: bash

        magnopy-scenario -if input_file.txt

    would be identical with "--input_file" being the long name and "-if" - the short one.

*   Keyword arguments without value

    This type of argument is very similar to the keyword argument with value. It is
    typically used for the True/False values. If the argument is not given, then the
    default value is assumed, if it is given, then the opposite is understood. For
    example, if the script has an argument defined with the keyword "--relative" and
    default value "False", then when user runs

    .. code-block:: bash

        magnopy-scenario

    the script will use ``relative = False``. However, if user runs the command

    .. code-block:: bash

        magnopy-scenario --relative

    the script will use ``relative = True``.
