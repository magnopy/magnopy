.. _development_documentation:

*************
Documentation
*************

The documentation of magnopy is build with |sphinx|_.

The best way to get a feeling about how the documentation of magnopy is structured is
to read the source code in the "docs/source/" directory and compare it's content and
structure with this webpage. If you have any doubts we encourage you to
:ref:`contact us <support>`.

Building the documentation
==========================

To build documentation simply run (provided that |GNU-make|_ command is available)

.. code-block:: bash

  make html

Alternatively, you can use the command

.. code-block:: bash

  sphinx-build -M html "docs/source" "docs/_build"

Documentation structure
=======================

Documentation of magnopy has two main parts

* User guide ("docs/source/user-guide/" directory)

  Hand-written |reStructuredText|_ files with the examples and explanation of the
  magnopy's capabilities and theory behind.

* API ("docs/source/api/" directory)

  Semi-automatically generated documentation of the source code, it is mostly build
  based on the docstrings of the source code using |sphinx-autodoc|_ and
  |sphinx-autosummary|_.

The rest of the documentation is located in the "docs/source/" directory and includes,
among other things

* "docs/source/conf.py" file

  The configuration file for |sphinx|_.

* "docs/source/index.rst" file

  The main page of the documentation. It includes the table of contents and the
  introduction to the magnopy.

* "docs/source/user-support/" directory

  The page with the information about how to get support for the users of magnopy.

* "docs/source/release-notes/" directory

  The release notes for each version of magnopy.

* "docs/source/development/" directory

  Root folder for the documentation for the development of magnopy.


Docstrings
==========

All public classes and functions have to have a docstring.
The docstring has to be written in the |numpydoc|_ style guide.

To get a feeling about the style you can read examples in the source code of magnopy.
