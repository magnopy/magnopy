# ================================== LICENSE ===================================
# Magnopy - Python package for magnons.
#
# Copyright (C) 2023 Magnopy Team
#
# e-mail: anry@uv.es, web: magnopy.org
#
# This program is free software: you  can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the  Free Software
# Foundation,  either  version 3  of the License,  or (at your option) any later
# version.
#
# This program is distributed in the  hope  that it will be useful,  but WITHOUT
# ANY WARRANTY;  without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the  GNU General Public License  along with
# this program.  If not, see <https://www.gnu.org/licenses/>.
# ================================ END LICENSE =================================


def _npn(nus, alphas) -> tuple:
    r"""
    Deduce ``n`` and ``p_n`` from canonical ``nus`` and ``alphas``.

    Does not perform validation of the input types or shapes.

    Parameters
    ----------
    nus : (n - 1,) tuple of (3,) tuple of int
        Unit cells of each site, excluding the first one that sits in ``(0, 0, 0)``.
    alphas : (n,) tuple of int
        Index of each site.

    Returns
    -------
    n : int
        Number of spin operators associated with the interaction defined by ``nus`` and
        ``alphas``.
    p_n : int
        Partition class of the interaction defined by ``nus`` and ``alphas``.
        ``1 <= p_n <= p(n)`` where ``p(n)`` is the number of integer partitions of
        ``n``.

    Notes
    -----
    A site is defined by the ``(nu, alpha)`` pair.

    ``p_n`` indexes integer partitions of ``n`` in descending orger (when sorted
    lexicographically).

    Full definition of ``p_n`` is writted in the docs, see
    :ref:`user-guide_theory-behind_spin-hamiltonian`.

    See Also
    --------
    get_npn
    """


class _InteractionParameters:
    r"""
    Sorted container for interaction parameters of a spin Hamiltonian.

    Knows nothing about units, convention or crystal. Stores whatever numbers it is
    given.

    Each parameter is identified by its ``nus`` and ``alphas``.
    From those the container derives ``specs = (n, p_n, nus, alphas)``,
    sorts by ``specs`` and stores it.

    Iterating over the container yields ``nus, alphas, parameter``.
    ``parameters.items()`` yields full form of the stored data: ``specs, parameter``.

    Notes
    -----

    Internally the container stores the data as ``[(specs, parameter), ...]``.

    **Promises.** Must hold for any container at any moment:

    * ``specs``: ``n, p_n`` exactly match their ``nus, alphas``.
    * ``len(nus) == len(alphas) - 1``, i.e. ``alphas[0]`` sits in ``(0, 0, 0)`` unit
      cell, whose nu is not stored. ``alphas[i]`` sits in ``nus[i-1]`` unit cell
      for ``i > 0``.
    * The container is sorted by ``specs``. Thus, each group of parameters, defined by
      the pair ``n, p_n`` is contiguous and can be found by bisection and iterated over.
    * ``specs`` are unique within the container. Two parameters describing the same
      interaction cannot coexist.
    * Shape of parameters: ``parameter.shape == (3,)*n`` at all times.

    **Bookkeeping.** ``n, p_n`` are stored (for speed), but computed by the container
    from ``nus, alphas``. These can not be supplied to the container from outside, nor
    edited. Boundaries of ``(n, p_n)`` groups are derived by bisection, not stored.

    **Ownership.** Parameters are copied on input. Iteration yields read-only views.

    **Duplicates.** Duplicate ``specs`` are excluded by design. Every operation that
    could produce one takes an explicit rule (for example ``on_duplicate``) for
    resolving it, so no repair step is needed.
    """
