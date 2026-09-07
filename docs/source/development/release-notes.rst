.. _development_release-notes:

**************************
Template for release notes
**************************

Here you can find the guide for writing the release notes.

Overall, always include links to the issues that are relevant to the release.


.. _development_release-notes_minor:

Minor release
=============

New minor release requires creation of the new file named "0.<minor>.rst" in the
"docs/source/release-notes" folder. The file should be based on the template below. Do
not forget to add this new file in the toctree inside the
"docs/source/release-notes/index.rst" file, at the top of the list.

.. code-block:: text

    .. _release-notes_0.<minor>:

    ***********
    Version 0.<minor>
    ***********


    <Follow template for micro release to describe 0.<minor>.0 release>


    Breaking changes wrt 0.<minor-1>
    ========================

    Deprecated features
    -------------------

    * First deprecated feature. Identify alternative if possible.
    * Second deprecated feature. Identify alternative if possible.
    * ...

    Removed features
    ----------------

    * First removed feature. Identify alternative if possible.
    * Second removed feature. Identify alternative if possible.
    * ...

    Change of behavior
    ------------------

    * First changed behavior.
    * Second changed behavior.
    * ...

    Transition guide
    ================

    <Describe how to transition from 0.<minor-1> to 0.<minor>. Include code snippets
    if possible.>



.. _development_release-notes_micro:

Micro release
=============

New micro release requires addition of the new section in the existing file
"docs/source/release-notes/0.<minor>.rst". The new section shall be added at the top of
the file. The section should be based on the template below.

.. code-block:: text

    0.<minor>.<micro>
    =====

    **Date**: <DD Month YYYY>

    New features
    ------------

    * First new feature.
    * Second new feature.
    * ...

    Bugfix
    ------

    * First fixed bug.
    * Second fixed bug.
    * ...

    Improvements
    ------------

    * First improvement.
    * Second improvement.
    * ...
