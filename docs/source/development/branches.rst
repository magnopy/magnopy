.. _development_branches:

************
Branch rules
************

Branches of the upstream repository
===================================

In the upstream repository the releases are made from the ``main`` branch. The
``main`` branch is the default branch in the upstream repository.

For each release we create a tag named ``vMajor.Minor.Micro`` in the upstream
repository.

Once new ``Minor`` release is made, the ``Major.Minor`` branch is created for the latest
``vMajor.Minor.Micro`` tag. For example, if the release history looks like

* 0.1.0
* 0.1.1
* 0.1.2
* 0.1.3
* 0.1.4
* 0.1.5
* 0.1.6
* 0.1.7
* 0.1.8
* 0.1.9

then  at the time of the ``v0.2.0`` release the ``0.1`` branch is created from the
``v0.1.9`` tag.

Branches of your local repository
=================================

When you clone the upstream repository, you will have the ``main`` branch checked out by
default.

When you want to contribute to the project, you should create a new branch from ``main``
and name it according to the feature you are working on. If you are fixing a bug, it is
wise to name the branch ``bugfix-issue_number`` when possible.

Branches of origin repository
=============================

At any moment of your work on the feature or bugfix, you can push changes to your forked
repository. With the first push, you should create a new branch in your forked
repository. For example, if you named your local branch as ``new-cool-feature``, type

.. code-block:: bash

  git push -u origin new-cool-feature

It will create a new branch in your forked repository and link your local branch to it.
For all next pushes you can simply type ``git push``.
