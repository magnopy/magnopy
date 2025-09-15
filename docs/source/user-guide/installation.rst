.. _user-guide_installation:

**********************
How to install magnopy
**********************

Requirement for magnopy installation are:

* |Python|_ (>=3.9)

And several libraries:

.. literalinclude:: ../../../requirements.txt


Magnopy can be installed with :ref:`pip <user-guide_installation_pip>` or from
:ref:`source <user-guide_installation_source>`.

Do you have Python?
===================

Most likely Python is already installed on your machine (if not check these links:
|Python-installation|_). One of the ways to check if you have python installed is to
run the command in your terminal

.. code-block:: bash

    python

If you see something like

.. code-block:: bash

    Python 3.10.9 (main, Dec 15 2022, 18:25:35) [Clang 14.0.0 (clang-1400.0.29.202)] on darwin
    Type "help", "copyright", "credits" or "license" for more information.
    >>>

then you have it.

In most cases ``python`` command launches python3, however if it launches python2, then
you may need to use ``python3`` instead (and ``pip3`` instead of ``pip``).

.. hint::
    Use ``exit()`` or press ``ctrl+D`` to close python console.

.. _user-guide_installation_pip:

.. hint::
    On linux and OSX systems one can create a virtual environment for the magnopy's installation with
    (use your version of python instead of ``python3.13`` if needed)

    .. code-block:: bash

        python3.13 -m venv .venv

    And activate it with

    .. code-block:: bash

        source .venv/bin/activate

Installation with pip
=====================

Magnopy is published and distributed from the |PYPI|_. To install it use the command
(you may need to use ``pip3``)

.. code-block:: bash

    pip install magnopy

.. hint::
    You can use the command

    .. code-block::

        magnopy

    to check whether the installation worked.

If there are any bugs with the installation, please drop a message to the developers
through one of our :ref:`support` channels and we would be happy to help.

Optionally, if you want magnopy to produce .png and .html files with graphics you can
install |plotly|_ and |matplotlib|_ manually or install them with magnopy as

.. code-block:: bash

    pip install magnopy[visual]

.. note::
    You may need to escape the ``[`` and ``]`` characters, because they are special
    characters in most shells. For example, in bash you can use backslash to escape them

    .. code-block:: bash

        pip install magnopy\[visual\]

    Or enclose full name in quotes

    .. code-block:: bash

        pip install "magnopy[visual]"

.. hint::
    If you are using |jupyter|_, then magnopy can be installed with

    .. code-block::

        %pip install magnopy

    or

    .. code-block::

        %pip install "magnopy[visual]"

    within it.

.. _user-guide_installation_source:

Installation from source
========================

Alternatively, source code of magnopy is publicly available under the GPL-3.0 license,
therefore it is possible to install it from source.

*   Clone the project to your local computer (in other words: download the source code)

    .. code-block:: bash

        git clone git@github.com:magnopy/magnopy.git

*   Change the directory

    .. code-block:: bash

        cd magnopy

*   Install the requirements (you may need to use ``pip3``)

    .. code-block:: bash

        pip install -r requirements.txt

*   To install magnopy, run (you may need to use ``pip3``)

    .. code-block:: bash

        pip install .

Update
======

If you want to update the package to the latest available version (|version|),
then use the command (you may need to use ``pip3``)

.. code-block:: bash

    pip install magnopy --upgrade
