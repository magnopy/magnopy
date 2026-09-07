.. _user-guide_installation:

**********************
How to install Magnopy
**********************

Requirements for Magnopy installation are

* |Python|_ (3.10 or higher)

And several third-party libraries

.. literalinclude:: ../../../requirements.txt


Magnopy can be installed with :ref:`pip <user-guide_installation_pip>` or from
:ref:`source <user-guide_installation_source>`.

Do you have Python?
===================

Most likely Python is already installed on your machine (if not, check
|Python-installation|_). One of the ways to check if you have Python installed
is to run the command in your terminal

.. code-block:: bash

    python

If you see something like

.. code-block:: bash

    Python 3.13.5 (main, Jun 11 2025, 15:36:57) [Clang 15.0.0 (clang-1500.1.0.2.5)] on darwin
    Type "help", "copyright", "credits" or "license" for more information.
    >>>

then you have it.

In most cases ``python`` command launches Python 3, however if it launches Python 2, then
you may need to use ``python3`` instead (and ``pip3`` instead of ``pip``).

.. hint::
    Use ``exit()`` or press ``Ctrl+D`` to close Python console.

.. _user-guide_installation_pip:

.. hint::

    On linux and macOS one can create a virtual environment for the Magnopy's
    installation with (use your version of Python instead of ``python3.13`` if
    needed)

    .. code-block:: bash

        python3.13 -m venv .venv

    And activate it with

    .. code-block:: bash

        source .venv/bin/activate

Installation with pip
=====================

Magnopy is published and distributed via |PYPI|_. To install it use the command
(you may need to use ``pip3``)

.. code-block:: bash

    pip install magnopy

.. hint::
    If you execute

    .. code-block::

        magnopy

    and see something like

    .. literalinclude:: memo.inc

    then the installation worked and Magnopy is ready to use.

If there are any bugs with the installation, please drop a message to the developers
through one of our :ref:`support <support>` channels and we would be happy to help.

Optionally, if you want Magnopy to be able to produce .png and .html files you can
install |plotly|_ and |matplotlib|_ manually or install them as Magnopy's dependencies

.. code-block:: bash

    pip install magnopy[visual]

.. note::
    You may need to escape the ``[`` and ``]`` characters, because they are special
    characters in most shells. For example, in bash you can use backslash

    .. code-block:: bash

        pip install magnopy\[visual\]

    Or enclose full name in quotes

    .. code-block:: bash

        pip install "magnopy[visual]"

.. hint::
    If you are using |jupyter|_, then Magnopy can be installed as

    .. code-block::

        %pip install magnopy

    or

    .. code-block::

        %pip install "magnopy[visual]"

    within it. You may need to restart you kernel as well.

.. _user-guide_installation_source:

Installation from source
========================

Source code of Magnopy is publicly available (see |repository|_),
therefore Magnopy can be installed from source.

*   Clone the project (in other words, download the source code)

    .. code-block:: bash

        git clone git@github.com:magnopy/magnopy.git

*   Change the directory

    .. code-block:: bash

        cd magnopy

*   Install the requirements (you may need to use ``pip3``)

    .. code-block:: bash

        pip install -r requirements.txt

*   To install Magnopy, run (you may need to use ``pip3``)

    .. code-block:: bash

        pip install .

    or

    .. code-block:: bash

        pip install ".[visual]"

Update
======

New versions of Magnopy are often published. We recommend updating from time to time to
get the new features and bug fixes.

To update Magnopy to the latest available version (|version|) use the command (you may
need to use ``pip3``)

.. code-block:: bash

    pip install magnopy --upgrade

or

.. code-block:: bash

    pip install "magnopy[visual]" --upgrade
