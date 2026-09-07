.. _development_tests:

*******
Testing
*******

In Magnopy we rely on |pytest|_, |hypothesis|_ and |doctest|_ for testing. Ideally, all
functions in Magnopy should be covered by unit tests and every code snippet in the
documentation should be a doctest.

Unit tests
==========

All unit tests are located in the "srs/magnopy/_tests/" directory. To run the
tests, you can use the following command (provided that the |GNU-make|_ command
is available)

.. code-block:: bash

  make test

Alternatively, you can run

.. code-block:: bash

  magnopy test

The structure of the "_tests/" directory loosely follows the structure of the
"src/magnopy/" directory.

For example if you've added a new function named "rotate()" to the
"src/magnopy/magnons/_dispersion.py" file, then you should add a new test function
named "test_rotate()" to the file "_tests/test_magnons/test_dispersion.py".

Doctests
========

Across the documentation there are many examples of how to use Magnopy with code
snippets. These code snippets are tested using |doctest|_, which ensures that the
documentation correctly reflects the actual behavior of the code. In order to run
doctests you need to build the :ref:`documentation <development_documentation>` and then
run doctests with the command (provided that the |GNU-make|_ command is available)

.. code-block:: bash

  make doctest

Alternatively, you can use the command

.. code-block:: bash

  sphinx-build -b doctest "docs/source" "docs/_build"
