.. _development:

*****************
Developer's guide
*****************

If you're interested in seeing who has already contributed to this project, please visit
our :ref:`Contributors page <development_contributors>`. We appreciate all contributions
and look forward to see your name on that list.

It is not necessary to be a programmer to contribute. You can help us with the
:ref:`new features <user-support_feature>`, :ref:`finding bugs <user-support_bug>` or
:ref:`documentation <development_documentation>`.

Please visit the topical guides listed below.


Topical guides
==============

.. toctree::
    :hidden:

    contributors

.. toctree::
    :maxdepth: 2

    deprecation-policy
    documentation
    tests
    origin-upstream
    branches
    theory-notes




Development process
===================

Here we summarize main steps for contributing code to magnopy. We assume that you have an
account on `<https://github.com>`_, and familiar with |Git|_.

Fork and clone
--------------

*   Go to the |repository|_ (:ref:`upstream <development_origin-upstream_upstream>`) and
    click on the "Fork" button. Now you have your own copy of the magnopy repository
    under your GitHub account (:ref:`origin <development_origin-upstream_origin>`).
*   Clone your fork to your local machine:

    -   If you are using ssh-key

        .. code-block:: bash

            git clone git@github.com:your-username/magnopy.git

    -    If you are not using ssh-key

        .. code-block:: bash

            git clone https://github.com/your-username/magnopy.git

*   Change the directory

    .. code-block:: bash

        cd magnopy

    Now you are in the root folder of you local repository
    (:ref:`local <development_origin-upstream_local>`).

*   Add the :ref:`upstream <development_origin-upstream>` repository to your
    :ref:`local <development_origin-upstream_local>` repository

    .. code-block:: bash

        git remote add upstream https://github.com/magnopy/magnopy.git

*   Pull the latest changes from the magnopy repository if necessary

    .. code-block:: bash

        git pull upstream main

    or

    .. code-block:: bash

        git fetch upstream
        git merge upstream/main

Set up the environment
----------------------

We recommend to use virtual environment (with |venv|_, for example). Once the virtual
environment is created and activated, you can install the requirements

*   Package dependencies

    .. code-block:: bash

        pip install -r requirements.txt

*   For the package development

    .. code-block:: bash

        pip install -r requirements-dev.txt

*   For the documentation

    .. code-block:: bash

        pip install -r docs/requirements.txt

*   For the tests

    .. code-block:: bash

        pip install -r tests/requirements.txt

.. note::
    For the linux and OSX systems there is a scenario defined. It installs all
    requirements. Note: it does not create a virtual environment for you

    .. code-block:: bash

        make requirements

Enable pre-commit
-----------------

We use |pre-commit|_ to enforce some rules on the code style.
To enable it, run the following command

.. code-block:: bash

    pre-commit install

Now, every time you commit something, pre-commit will check it for you. If some of the
checks fail, pre-commit will automatically fix them and abort the commit. You need to
add the fixed files to the staging area and commit again.

.. hint::
    If you want to run pre-commit manually, you can use

    .. code-block:: bash

        pre-commit run --all-files

Develop your contribution
-------------------------

*   Create a :ref:`dedicated branch <development_branches>` for your feature, that you are
    going to develop

    .. code-block:: bash

        git checkout -b feature-name

*   Develop your contribution. Commit your progress locally
    (|git-add|_ and |git-commit|_). Use |good-commit-messages|_. Write
    :ref:`tests <development_tests>`. Write
    :ref:`documentation <development_documentation>`.

Submit your contribution
------------------------

*   Push the changes to your forked repository

    .. code-block:: bash

      git push origin feature-name

*   Go to your forked repository on GitHub and click on the green "Compare & pull request"
    button. Describe your contribution and submit the pull request. Please mention the
    issue number if your contribution is related to any.

Review and merge
----------------

*   Once the pull request is submitted, the code will be reviewed. If there are any
    comments, please fix them. You can make the changes locally, commit them and push to
    the same branch of :ref:`origin <development_origin-upstream_origin>` repository and
    they will be added to the pull request automatically.
*   Once the pull request is approved, it will be merged to the
    main branch of the |repository|_.
