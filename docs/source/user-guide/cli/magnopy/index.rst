.. _user-guide_cli_magnopy:

*******
magnopy
*******

The CLI command ``magnopy`` can display information about the package or run its test suite.

Getting help
============

We recommend getting the accurate and full list of the script's parameters that
reflect the installed version of Magnopy with the command

.. code-block::

    magnopy --help

which outputs to the standard output channel (console or terminal) Magnopy's metadata and
*full* list of the script's arguments. Here is an example of this output

.. hint::

    Go :ref:`here <user-guide_cli_common-notes_read-help>` to learn how to read this help
    message.

.. literalinclude:: help.inc
    :language: text

.. _user-guide_cli_magnopy_version:

Version
=======

To get the version of the Magnopy's installation, use

.. code-block:: bash

   magnopy --version

which shall output something similar to

.. literalinclude:: version.inc
    :language: text

.. _user-guide_cli_magnopy_test:

magnopy test
============

To run the test suite of Magnopy, use

.. code-block:: bash

   magnopy test

The runtime of the test suite depends on the hardware and can take several minutes.
Progress is printed to the console.

.. _user-guide_cli_magnopy_logo:

magnopy logo
============

To display Magnopy's logo along with the information about the package, use

.. code-block:: bash

   magnopy logo

which shall output something similar to

.. literalinclude:: logo.inc
    :language: text
