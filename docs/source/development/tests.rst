.. _development_tests:

*******
Testing
*******

In magnopy we rely on |pytest|_, |hypothesis|_ and |doctest|_ for testing. Ideally, all
functions in magnopy should be covered by unit tests and every code snippet in the
documentations should be a doctest.

Unit tests
==========

All unit tests are located in the "tests/" directory. To run the tests, you can use the
following command (provided that the |GNU-make|_ command is available)

.. code-block:: bash

  make test

Alternatively, you can run

.. code-block:: bash

  pytest -s

Structure of the "test/" directory loosely follows the structure of the "src/magnopy/"
directory.

For example if you've added a new function named "rotate()" to the
"src/magnopy/magnons/_dispersion.py" file, then you should add a new test function
named "test_rotate()" to the file "tests/test_magnons/test_dispersion.py".

Doctests
========

Across the documentation there are many examples of how to use magnopy with code
snippets. These code snippets are tested using |doctest|_, which ensure that the
documentation correctly reflects an actual behavior of the code. In order to run
doctests you need to build the :ref:`documentation <development_documentation>` and then
run doctests with the command (provided that the |GNU-make|_ command is available)

.. code-block:: bash

  make doctest

Alternatively, you can use the command

.. code-block:: bash

  sphinx-build -b doctest "docs/source" "docs/_build"
