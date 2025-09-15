.. _development_deprecation_policy:

***********************************
Deprecation policy (for beta stage)
***********************************

**Created**: 04 September 2025

**Status**: Active

**Author**: Andrey Rybakov (`adrybakov <https://github.com/adrybakov>`_)

All versions of magnopy code that fit ``0.*`` are considered as beta stage of magnopy's
lifecycle. This stage allows for more dynamic evolution of the code with frequent
changes and supply of new features.

Once the release ``v1.0.0`` is published the magnopy will be considered stable and new,
more restrictive deprecation police will be published.


Deprecation policy
==================
Removal of public features shall be done in two stages. At the beta stage


Stage 1: Deprecation (desired, but optional)
--------------------------------------------

This first stage is optional for the beta stage of magnopy's development.

*   Shall be skipped **only** if its realization is incompatible with the new code.
*   Shall be mentioned in the release notes of the version when deprecation was
    introduced.
*   Shall issue a ``DeprecationWarning`` with

    * Information on what user should use instead
    * Version in which it has been deprecated

*   Shall not be introduced in the micro release. Can be introduced in major or minor
    release.
*   Shall be mentioned in the documentation (``.. deprecated::``).

Stage 2: Removal
----------------

*   Shall be mentioned in the release notes of the version where the removal was done.
*   Shall be done after six months since the release of the version in which deprecation
    was introduced.
*   Shall be done only in the major or minor release. Shall not be done in the micro
    release.


New features
============

*   Any new public feature both in python library and in command line interface have to
    have ``.. versionadded:: v<major>.<minor>.<micro>`` in its documentation.
