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

    Raises
    ------
    ValueError
       If the given ``nus, alphas`` describe a term with ``n >= 5`` or ``alphas`` is
       empty.

    Notes
    -----
    A site is defined by the ``(nu, alpha)`` pair.

    ``p_n`` indexes integer partitions of ``n`` in descending order (when sorted
    lexicographically).

    Full definition of ``p_n`` is written in the docs, see
    :ref:`user-guide_theory-behind_spin-hamiltonian`.

    See Also
    --------
    get_npn
    """
    n = len(alphas)

    # optimization: n=2 are more frequent
    if n == 2:
        return 2, (1 if (nus[0] == (0, 0, 0) and alphas[0] == alphas[1]) else 2)

    if n == 1:
        return 1, 1

    if n == 3:
        return 3, len(
            {
                ((0, 0, 0), alphas[0]),
                (nus[0], alphas[1]),
                (nus[1], alphas[2]),
            }
        )

    if n == 4:
        sites = (
            ((0, 0, 0), alphas[0]),
            (nus[0], alphas[1]),
            (nus[1], alphas[2]),
            (nus[2], alphas[3]),
        )
        distinct = len({*sites})

        if distinct == 1:
            return 4, 1
        if distinct == 3:
            return 4, 4
        if distinct == 4:
            return 4, 5
        return 4, (2 if sites.count(sites[0]) in (1, 3) else 3)

    raise ValueError(f"Expected 1 <= n <= 4 spin operators, got n = {n}.")


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

    def __init__(self) -> None:
        self._container = []

    ############################################################################
    #                                Bisecting                                 #
    ############################################################################

    def _index(self, nus, alphas) -> tuple:
        r"""
        Finds an index of the parameter in container.

        Does not perform validation of the input types or shapes.

        Derives ``specs`` via :py:func:`._npn`, then bisects.

        Parameters
        ----------
        nus : (n - 1,) tuple of (3,) tuple of int
            Unit cells of each site.
        alphas : (n,) tuple of int
            Index of each site.

        Returns
        -------
        index : int
            Position of the parameter in the container. If ``found`` is ``True``,
            the parameter is stored as ``self._container[index]``.
            If ``found`` is ``False``, ``index`` is the position at which it would
            be inserted to keep the container sorted by ``specs`` (i.e.
            ``self._container.insert(index, ...)`` preserves the order).
            ``0 <= index <= len(self)``.
        found : bool
            ``True`` if a parameter with these canonical ``specs`` is present.
            ``False`` otherwise.

        """

    def get_index(self, nus, alphas):
        r"""
        Finds an index of the parameter in container.

        Parameters
        ----------
        nus : (n,) or (n-1,) iterable of (3,) iterable of int
            Unit cells for each site.
        alphas : (n,) iterable of int
            Indices of sites within each unit cell. See notes of :py:meth:`.add` for
            details.

        Returns
        -------
        index : int
            Position of the parameter in the container.

        Raises
        ------
        ValueError
            If there is no such item.
        """
        raise NotImplementedError

    ############################################################################
    #                         Single-parameter methods                         #
    ############################################################################

    def add(self, nus, alphas, parameter) -> None:
        r"""
        Adds a single parameter to the container.

        Parameters
        ----------
        nus : (n,) or (n-1,) iterable of (3,) iterable of int
            Unit cells for each site.
        alphas : (n,) iterable of int
            Indices of sites within unit cell. If ``len(nus) == len(alphas)``,
            then no condition is enforced on ``nus[0]``, all ``nus`` are shifted so
            ``nus[0] == (0,0,0)`` for storage. If ``len(nus) == len(alphas) - 1``, then
            ``alphas[0]`` sits in ``(0,0,0)``, ``alphas[i]`` sits in ``nus[i-1]`` for
            ``i >= 1``.
        parameter : (3,)*n |array-like|_
            Tensor of the interaction parameter. User input is copied, thus mutation
            of the original array does not affect stored parameter.

        Raises
        ------
        ValueError
            If canonical ``specs`` generated from ``nus, alphas`` are already present
            in the container.

        See Also
        --------
        extend
            For bulk additions.
        set
            For rewriting parameters.
        """
        raise NotImplementedError

    def set(self, nus, alphas, parameter) -> None:
        r"""
        Sets a single parameter.

        Adds if absent, rewrites if present.

        Parameters
        ----------
        nus : (n,) or (n-1,) iterable of (3,) iterable of int
            Unit cells for each site.
        alphas : (n,) iterable of int
            Indices of sites within each unit cell. See notes of :py:meth:`.add` for
            details.
        parameter : (3,)*n |array-like|_
            Tensor of the interaction parameter. User's input is copied, thus mutation
            of the original array does not affect stored parameter.

        See Also
        --------
        add
        extend
            For bulk additions.
        """
        raise NotImplementedError

    def remove(self, nus, alphas) -> bool:
        r"""
        Removes a single parameter from the container.

        Parameters
        ----------
        nus : (n,) or (n-1,) iterable of (3,) iterable of int
            Unit cells for each site.
        alphas : (n,) iterable of int
            Indices of sites within each unit cell. See notes of :py:meth:`.add` for
            details.

        Returns
        -------
        status : bool
            ``True`` if parameter with such spec was present and is now removed.
            ``False`` if the parameter with such specs was not present.
        """
        raise NotImplementedError

    ############################################################################
    #                         Many-parameters methods                          #
    ############################################################################

    def extend(self, items, on_duplicate="raise", on_existing="raise") -> None:
        r"""
        Extend the container with many parameters at once.

        More efficient than :py:meth:`.add` on many parameters.

        Parameters
        ----------
        items : iterable
            Iterable of the interaction parameters of the form
            ``[[nus, alphas, parameters], ...]`` (list is used for illustration).
            See :py:meth:`.add` for description of the elements.
        on_duplicate : str, default "raise"
            What to do if ``items`` contains duplicate specs after canonilisation.
            Case insensitive. Supported values are:

            * "raise" (default). Raises ``ValueError``.
            * "sum". Sums tensors for the parameters with the same specs.
            * "mean". Computes arithmetic mean between tensors with the same specs.

        on_existing : str, default "raise"
            What to do if specs from ``items`` already present in the container.
            Applied after ``on_duplicate`` is resolved. Case insensitive. Supported
            values are:

            * "raise" (default). Raises ``ValueError``.
            * "sum". Sums the incoming tensor from ``items`` with the existing tensor
              from the container.
            * "replace". Replaces existing tensor in the container with the incoming
              tensor from ``items``.

        Raises
        ------
        ValueError
            See ``on_duplicate`` and ``on_existing``.

        See Also
        --------
        add
        set
        """
        raise NotImplementedError

    ############################################################################
    #                                  Other                                   #
    ############################################################################

    def __len__(self) -> int:
        r"""
        Total amount of parameters in the container.
        """
        return len(self._container)

    def __contains__(self, item) -> bool:
        r"""
        Checks if such parameter is in the container.

        Parameters
        ----------
        item : tuple
            ``item = (nus, alphas)``. See notes of :py:meth:`.add` for details.

        Returns
        -------
        answer: bool
            ``True`` if such parameter is present. ``False`` otherwise.
        """
        raise NotImplementedError

    ############################################################################
    #                          Arithmetic operations                           #
    ############################################################################

    def __add__(self, other):
        r"""
        Merge two containers, sums parameters with the same specs.
        """
        raise NotImplementedError

    def __mul__(self, number):
        r"""
        Multiply all parameters by a number.

        Parameters
        ----------
        number : int | float
        """
        raise NotImplementedError

    def __rmul__(self, number):
        r"""
        Multiply all parameters by a number.

        Parameters
        ----------
        number : int | float
        """
        return self.__mul__(number=number)

    def __sub__(self, other):
        r"""
        Merge two containers, subtract parameters of ``other`` from parameters
        of ``self`` with the same specs.
        """
        return self + (-1) * other

    ############################################################################
    #                                 Copying                                  #
    ############################################################################

    def __deepcopy__(self, memo):
        r"""
        Deep copy of the container.

        memo unused, since the structure is acyclic. Left to satisfy the protocol.
        """

        return self.copy()

    def copy(self):
        r"""
        Copy of the object (deep).

        Returns
        -------
        copied_container : _InteractionParameters
            Deep copy of the container.
        """

        raise NotImplementedError
