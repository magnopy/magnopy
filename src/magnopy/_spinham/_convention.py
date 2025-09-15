# ================================== LICENSE ===================================
# Magnopy - Python package for magnons.
# Copyright (C) 2023-2025 Magnopy Team
#
# e-mail: anry@uv.es, web: magnopy.org
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
# ================================ END LICENSE =================================


R"""
Convention of spin Hamiltonian
"""

from magnopy._constants._spinham_conventions import _CONVENTIONS
from magnopy._exceptions import ConventionError

# Save local scope at this moment
old_dir = set(dir())
old_dir.add("old_dir")


class Convention:
    R"""
    Convention of the spin Hamiltonian.

    For the detailed description of the convention problem see :ref:`user-guide_theory-behind_convention-problem`.

    Parameters
    ----------
    multiple_counting : bool, optional
        Whether the pairs of spins are counted multiple times in the Hamiltonian's sums.
    spin_normalized : bool, optional
        Whether spin vectors/operators are normalized to 1. If ``True``, then spin
        vectors/operators are normalized.
    c1 : float, optional
        Numerical factor before the (one spin & one site) term of the Hamiltonian.
    c21 : float, optional
        Numerical factor before the (two spins & one site) term of the Hamiltonian.
    c22 : float, optional
        Numerical factor before the (two spins & two sites) term of the Hamiltonian.
    c31 : float, optional
        Numerical factor before the (three spins & one site) term of the Hamiltonian.
    c32 : float, optional
        Numerical factor before the (three spins & two sites) term of the Hamiltonian.
    c33 : float, optional
        Numerical factor before the (three spins & three sites) term of the Hamiltonian.
    c41 : float, optional
        Numerical factor before the (four spins & one site) term of the Hamiltonian.
    c421 : float, optional
        Numerical factor before the (four spins & two sites & 1+3) term of the Hamiltonian.
    c422 : float, optional
        Numerical factor before the (four spins & two sites & 2+2) term of the Hamiltonian.
    c43 : float, optional
        Numerical factor before the (four spins & three sites) term of the Hamiltonian.
    c44 : float, optional
        Numerical factor before the (four spins & four sites) term of the Hamiltonian.
    name : str, default "custom"
        A label for the convention. Any string, case-insensitive.

    Examples
    --------

    .. doctest::

        >>> from magnopy import Convention
        >>> n1 = Convention(True, True, c21=1, c22=-0.5)
        >>> n2 = Convention(False, True, c21=1, c22=-0.5)
        >>> n3 = Convention(False, True, c22=-0.5)
        >>> n1.multiple_counting
        True
        >>> n1 == n2
        False
        >>> n3.c21
        Traceback (most recent call last):
        ...
        magnopy._exceptions.ConventionError: Convention of spin Hamiltonian has an undefined property 'c21':
        custom convention where
          * Bonds are counted once in the sum;
          * Spin vectors are normalized to 1;
          * Undefined c1 factor;
          * Undefined c21 factor;
          * c22 = -0.5;
          * Undefined c31 factor;
          * Undefined c32 factor;
          * Undefined c33 factor;
          * Undefined c41 factor;
          * Undefined c421 factor;
          * Undefined c422 factor;
          * Undefined c43 factor;
          * Undefined c44 factor.
        >>> n3.name
        'custom'

    """

    __slots__ = (
        "_multiple_counting",
        "_spin_normalized",
        "_c1",
        "_c21",
        "_c22",
        "_c31",
        "_c32",
        "_c33",
        "_c41",
        "_c421",
        "_c422",
        "_c43",
        "_c44",
        "_name",
    )

    def __init__(
        self,
        multiple_counting: bool = None,
        spin_normalized: bool = None,
        c1: float = None,
        c21: float = None,
        c22: float = None,
        c31: float = None,
        c32: float = None,
        c33: float = None,
        c41: float = None,
        c421: float = None,
        c422: float = None,
        c43: float = None,
        c44: float = None,
        name: str = "custom",
    ) -> None:
        if multiple_counting is not None:
            self._multiple_counting = bool(multiple_counting)
        else:
            self._multiple_counting = None

        if spin_normalized is not None:
            self._spin_normalized = bool(spin_normalized)
        else:
            self._spin_normalized = None

        if c1 is not None:
            self._c1 = float(c1)
        else:
            self._c1 = None

        if c21 is not None:
            self._c21 = float(c21)
        else:
            self._c21 = None

        if c22 is not None:
            self._c22 = float(c22)
        else:
            self._c22 = None

        if c31 is not None:
            self._c31 = float(c31)
        else:
            self._c31 = None

        if c32 is not None:
            self._c32 = float(c32)
        else:
            self._c32 = None

        if c33 is not None:
            self._c33 = float(c33)
        else:
            self._c33 = None

        if c41 is not None:
            self._c41 = float(c41)
        else:
            self._c41 = None

        if c421 is not None:
            self._c421 = float(c421)
        else:
            self._c421 = None

        if c422 is not None:
            self._c422 = float(c422)
        else:
            self._c422 = None

        if c43 is not None:
            self._c43 = float(c43)
        else:
            self._c43 = None

        if c44 is not None:
            self._c44 = float(c44)
        else:
            self._c44 = None

        self._name = str(name).lower()

    def summary(self, return_as_string=False):
        r"""
        Gives human-readable summary of the convention.

        Parameters
        ----------
        return_as_string : bool, default False
            Whether to print or return a ``str``. If ``True``, then return an ``str``.
            If ``False``, then print it.

        Examples
        --------

        .. doctest::

            >>> from magnopy import Convention
            >>> n1 = Convention(True, True, c21=1, c22=-0.5)
            >>> n1.summary()
            custom convention where
              * Bonds are counted multiple times in the sum;
              * Spin vectors are normalized to 1;
              * Undefined c1 factor;
              * c21 = 1.0;
              * c22 = -0.5;
              * Undefined c31 factor;
              * Undefined c32 factor;
              * Undefined c33 factor;
              * Undefined c41 factor;
              * Undefined c421 factor;
              * Undefined c422 factor;
              * Undefined c43 factor;
              * Undefined c44 factor.
        """

        summary = [f"{self.name} convention where"]

        if self._multiple_counting is None:
            summary.append("  * Undefined multiple counting;")
        elif self._multiple_counting:
            summary.append("  * Bonds are counted multiple times in the sum;")
        else:
            summary.append("  * Bonds are counted once in the sum;")

        if self._spin_normalized is None:
            summary.append("  * Undefined spin normalization;")
        elif self._spin_normalized:
            summary.append("  * Spin vectors are normalized to 1;")
        else:
            summary.append("  * Spin vectors are not normalized;")

        # One spin
        if self._c1 is None:
            summary.append("  * Undefined c1 factor;")
        else:
            summary.append(f"  * c1 = {self._c1};")

        # Two spins
        if self._c21 is None:
            summary.append("  * Undefined c21 factor;")
        else:
            summary.append(f"  * c21 = {self._c21};")

        if self._c22 is None:
            summary.append("  * Undefined c22 factor;")
        else:
            summary.append(f"  * c22 = {self._c22};")

        # Three spins
        if self._c31 is None:
            summary.append("  * Undefined c31 factor;")
        else:
            summary.append(f"  * c31 = {self._c31};")

        if self._c32 is None:
            summary.append("  * Undefined c32 factor;")
        else:
            summary.append(f"  * c32 = {self._c32};")

        if self._c33 is None:
            summary.append("  * Undefined c33 factor;")
        else:
            summary.append(f"  * c33 = {self._c33};")

        # Four spins
        if self._c41 is None:
            summary.append("  * Undefined c41 factor;")
        else:
            summary.append(f"  * c41 = {self._c41};")

        if self._c421 is None:
            summary.append("  * Undefined c421 factor;")
        else:
            summary.append(f"  * c421 = {self._c421};")

        if self._c422 is None:
            summary.append("  * Undefined c422 factor;")
        else:
            summary.append(f"  * c422 = {self._c422};")

        if self._c43 is None:
            summary.append("  * Undefined c43 factor;")
        else:
            summary.append(f"  * c43 = {self._c43};")

        if self._c44 is None:
            summary.append("  * Undefined c44 factor.")
        else:
            summary.append(f"  * c44 = {self._c44}.")

        summary = ("\n").join(summary)

        if return_as_string:
            return summary

        print(summary)

    @property
    def name(self) -> str:
        r"""
        A label for the convention. Any string, case-insensitive.
        """

        return self._name

    @name.setter
    def name(self, new_value: str):
        self._name = str(new_value).lower()

    ################################################################################
    #                               Multiple counting                              #
    ################################################################################
    @property
    def multiple_counting(self) -> bool:
        r"""
        Whether the pairs of spins are counted multiple times in the Hamiltonian's sums.

        If ``True``, then pairs are counted multiple times.
        """
        if self._multiple_counting is None:
            raise ConventionError(convention=self, property="multiple_counting")
        return self._multiple_counting

    @multiple_counting.setter
    def multiple_counting(self, new_value: bool):
        raise AttributeError(
            "It is intentionally forbidden to set properties of convention. "
            "Use correct methods of SpinHamiltonian class to change convention."
        )

    ################################################################################
    #                            Normalization of spins                            #
    ################################################################################
    @property
    def spin_normalized(self) -> bool:
        r"""
        Whether spin vectors/operators are normalized to 1.

        If ``True``, then spin vectors/operators are normalized.
        """
        if self._spin_normalized is None:
            raise ConventionError(convention=self, property="spin_normalized")
        return self._spin_normalized

    @spin_normalized.setter
    def spin_normalized(self, new_value: bool):
        raise AttributeError(
            "It is intentionally forbidden to set properties of convention. "
            "Use correct methods of SpinHamiltonian class to change convention."
        )

    ################################################################################
    #                                   One spin                                   #
    ################################################################################
    @property
    def c1(self) -> float:
        r"""
        Numerical factor before the (one spin & one site) sum of the Hamiltonian.
        """
        if self._c1 is None:
            raise ConventionError(convention=self, property="c1")
        return self._c1

    @c1.setter
    def c1(self, new_value: float):
        raise AttributeError(
            "It is intentionally forbidden to set properties of convention. "
            "Use correct methods of SpinHamiltonian class to change convention."
        )

    ################################################################################
    #                                   Two spins                                  #
    ################################################################################
    @property
    def c21(self) -> float:
        r"""
        Numerical factor before the (two spins & one site) sum of the Hamiltonian.
        """
        if self._c21 is None:
            raise ConventionError(convention=self, property="c21")
        return self._c21

    @c21.setter
    def c21(self, new_value: float):
        raise AttributeError(
            "It is intentionally forbidden to set properties of convention. "
            "Use correct methods of SpinHamiltonian class to change convention."
        )

    @property
    def c22(self) -> float:
        r"""
        Numerical factor before the (two spins & two sites) sum of the Hamiltonian.
        """
        if self._c22 is None:
            raise ConventionError(convention=self, property="c22")
        return self._c22

    @c22.setter
    def c22(self, new_value: float):
        raise AttributeError(
            "It is intentionally forbidden to set properties of convention. "
            "Use correct methods of SpinHamiltonian class to change convention."
        )

    ################################################################################
    #                                  Three spins                                 #
    ################################################################################
    @property
    def c31(self) -> float:
        r"""
        Numerical factor before the (three spins & one site) sum of the Hamiltonian.
        """
        if self._c31 is None:
            raise ConventionError(convention=self, property="c31")
        return self._c31

    @c31.setter
    def c31(self, new_value: float):
        raise AttributeError(
            "It is intentionally forbidden to set properties of convention. "
            "Use correct methods of SpinHamiltonian class to change convention."
        )

    @property
    def c32(self) -> float:
        r"""
        Numerical factor before the (three spins & two sites) sum of the Hamiltonian.
        """
        if self._c32 is None:
            raise ConventionError(convention=self, property="c32")
        return self._c32

    @c32.setter
    def c32(self, new_value: float):
        raise AttributeError(
            "It is intentionally forbidden to set properties of convention. "
            "Use correct methods of SpinHamiltonian class to change convention."
        )

    @property
    def c33(self) -> float:
        r"""
        Numerical factor before the (three spins & three sites) sum of the Hamiltonian.
        """
        if self._c33 is None:
            raise ConventionError(convention=self, property="c33")
        return self._c33

    @c33.setter
    def c33(self, new_value: float):
        raise AttributeError(
            "It is intentionally forbidden to set properties of convention. "
            "Use correct methods of SpinHamiltonian class to change convention."
        )

    ################################################################################
    #                                  Four spins                                  #
    ################################################################################
    @property
    def c41(self) -> float:
        r"""
        Numerical factor before the (four spins & one site) sum of the Hamiltonian.
        """
        if self._c41 is None:
            raise ConventionError(convention=self, property="c41")
        return self._c41

    @c41.setter
    def c41(self, new_value: float):
        raise AttributeError(
            "It is intentionally forbidden to set properties of convention. "
            "Use correct methods of SpinHamiltonian class to change convention."
        )

    @property
    def c421(self) -> float:
        r"""
        Numerical factor before the (four spins & two sites (1+3)) sum of the Hamiltonian.
        """
        if self._c421 is None:
            raise ConventionError(convention=self, property="c421")
        return self._c421

    @c421.setter
    def c421(self, new_value: float):
        raise AttributeError(
            "It is intentionally forbidden to set properties of convention. "
            "Use correct methods of SpinHamiltonian class to change convention."
        )

    @property
    def c422(self) -> float:
        r"""
        Numerical factor before the (four spins & two sites (2+2)) sum of the Hamiltonian.
        """
        if self._c422 is None:
            raise ConventionError(convention=self, property="c422")
        return self._c422

    @c422.setter
    def c422(self, new_value: float):
        raise AttributeError(
            "It is intentionally forbidden to set properties of convention. "
            "Use correct methods of SpinHamiltonian class to change convention."
        )

    @property
    def c43(self) -> float:
        r"""
        Numerical factor before the (four spins & three sites) sum of the Hamiltonian.
        """
        if self._c43 is None:
            raise ConventionError(convention=self, property="c43")
        return self._c43

    @c43.setter
    def c43(self, new_value: float):
        raise AttributeError(
            "It is intentionally forbidden to set properties of convention. "
            "Use correct methods of SpinHamiltonian class to change convention."
        )

    @property
    def c44(self) -> float:
        r"""
        Numerical factor before the (four spins & four sites) sum of the Hamiltonian.
        """
        if self._c44 is None:
            raise ConventionError(convention=self, property="c44")
        return self._c44

    @c44.setter
    def c44(self, new_value: float):
        raise AttributeError(
            "It is intentionally forbidden to set properties of convention. "
            "Use correct methods of SpinHamiltonian class to change convention."
        )

    def __eq__(self, other):
        # Note semi-private attributes are compared intentionally, as
        # public ones will raise an error if not defined
        # If attributes are not defined in both conventions,
        # then that attribute is considered equal.
        return (
            self._multiple_counting == other._multiple_counting
            and self._spin_normalized == other._spin_normalized
            and self._c1 == other._c1
            and self._c21 == other._c21
            and self._c22 == other._c22
            and self._c31 == other._c31
            and self._c32 == other._c32
            and self._c33 == other._c33
            and self._c41 == other._c41
            and self._c421 == other._c421
            and self._c422 == other._c422
            and self._c43 == other._c43
            and self._c44 == other._c44
        )

    @staticmethod
    def get_predefined(name: str):
        r"""
        Returns one of the pre-defined conventions.

        Parameters
        ----------
        name : str
            Name of the desired pre-defined convention. Supported are

            * "tb2j"
            * "grogu"
            * vampire"
            * "spinw"

            Case-insensitive.

        Returns
        -------
        convention : :py:class:`.Convention`

        Examples
        --------

        .. doctest::

            >>> import magnopy
            >>> tb2j = magnopy.Convention.get_predefined("TB2J")
            >>> tb2j.summary()
            tb2j convention where
              * Bonds are counted multiple times in the sum;
              * Spin vectors are normalized to 1;
              * Undefined c1 factor;
              * c21 = -1.0;
              * c22 = -1.0;
              * Undefined c31 factor;
              * Undefined c32 factor;
              * Undefined c33 factor;
              * Undefined c41 factor;
              * Undefined c421 factor;
              * Undefined c422 factor;
              * Undefined c43 factor;
              * Undefined c44 factor.
            >>> grogu = magnopy.Convention.get_predefined("GROGU")
            >>> grogu.summary()
            grogu convention where
              * Bonds are counted multiple times in the sum;
              * Spin vectors are normalized to 1;
              * Undefined c1 factor;
              * c21 = 1.0;
              * c22 = 0.5;
              * Undefined c31 factor;
              * Undefined c32 factor;
              * Undefined c33 factor;
              * Undefined c41 factor;
              * Undefined c421 factor;
              * Undefined c422 factor;
              * Undefined c43 factor;
              * Undefined c44 factor.
            >>> vampire = magnopy.Convention.get_predefined("Vampire")
            >>> vampire.summary()
            vampire convention where
              * Bonds are counted multiple times in the sum;
              * Spin vectors are normalized to 1;
              * Undefined c1 factor;
              * c21 = -1.0;
              * c22 = -0.5;
              * Undefined c31 factor;
              * Undefined c32 factor;
              * Undefined c33 factor;
              * Undefined c41 factor;
              * Undefined c421 factor;
              * Undefined c422 factor;
              * Undefined c43 factor;
              * Undefined c44 factor.
            >>> spinW = magnopy.Convention.get_predefined("spinW")
            >>> spinW.summary()
            spinw convention where
              * Bonds are counted multiple times in the sum;
              * Spin vectors are not normalized;
              * Undefined c1 factor;
              * c21 = 1.0;
              * c22 = 1.0;
              * Undefined c31 factor;
              * Undefined c32 factor;
              * Undefined c33 factor;
              * Undefined c41 factor;
              * Undefined c421 factor;
              * Undefined c422 factor;
              * Undefined c43 factor;
              * Undefined c44 factor.
        """

        name = name.lower()

        if name not in _CONVENTIONS:
            raise ValueError(
                f'"{name}" convention is undefined. Supported are\n - '
                + "\n - ".join([f'"{key}"' for key in _CONVENTIONS])
            )

        return Convention(name=name, **_CONVENTIONS[name])

    def get_modified(
        self,
        multiple_counting: bool = None,
        spin_normalized: bool = None,
        c1: float = None,
        c21: float = None,
        c22: float = None,
        c31: float = None,
        c32: float = None,
        c33: float = None,
        c41: float = None,
        c421: float = None,
        c422: float = None,
        c43: float = None,
        c44: float = None,
        name: str = None,
    ):
        r"""
        Returns the new instance of the :py:class:`.Convention` class based on the called
        one with changed given properties.

        Parameters
        ----------
        multiple_counting : bool, optional
            Whether the pairs of spins are counted multiple times in the Hamiltonian's sums.
            Modified to the given value, if None, then kept the same as in the original convention.
        spin_normalized : bool, optional
            Whether spin vectors/operators are normalized to 1. If ``True``, then spin
            vectors/operators are normalized.
            Modified to the given value, if None, then kept the same as in the original convention.
        c1 : float, optional
            Numerical factor before the (one spin & one site) term of the Hamiltonian.
            Modified to the given value, if None, then kept the same as in the original convention.
        c21 : float, optional
            Numerical factor before the (two spins & one site) term of the Hamiltonian.
            Modified to the given value, if None, then kept the same as in the original convention.
        c22 : float, optional
            Numerical factor before the (two spins & two sites) term of the Hamiltonian.
            Modified to the given value, if None, then kept the same as in the original convention.
        c31 : float, optional
            Numerical factor before the (three spins & one site) term of the Hamiltonian.
            Modified to the given value, if None, then kept the same as in the original convention.
        c32 : float, optional
            Numerical factor before the (three spins & two sites) term of the Hamiltonian.
            Modified to the given value, if None, then kept the same as in the original convention.
        c33 : float, optional
            Numerical factor before the (three spins & three sites) term of the Hamiltonian.
            Modified to the given value, if None, then kept the same as in the original convention.
        c41 : float, optional
            Numerical factor before the (four spins & one site) term of the Hamiltonian.
            Modified to the given value, if None, then kept the same as in the original convention.
        c421 : float, optional
            Numerical factor before the (four spins & two sites & 1+3) term of the Hamiltonian.
            Modified to the given value, if None, then kept the same as in the original convention.
        c422 : float, optional
            Numerical factor before the (four spins & two sites & 2+2) term of the Hamiltonian.
            Modified to the given value, if None, then kept the same as in the original convention.
        c43 : float, optional
            Numerical factor before the (four spins & three sites) term of the Hamiltonian.
            Modified to the given value, if None, then kept the same as in the original convention.
        c44 : float, optional
            Numerical factor before the (four spins & four sites) term of the Hamiltonian.
            Modified to the given value, if None, then kept the same as in the original convention.
        name : str, optional
            A label for the convention. Any string, case-insensitive.
            Modified to the given value, if None, then kept the same as in the original convention.

        Examples
        --------

        .. doctest::

            >>> import magnopy
            >>> conv = magnopy.Convention(
            ...     name="original",
            ...     multiple_counting=True,
            ...     spin_normalized=False,
            ...     c1=1,
            ...     c21=1,
            ...     c22=-1,
            ... )
            >>> conv.summary()
            original convention where
              * Bonds are counted multiple times in the sum;
              * Spin vectors are not normalized;
              * c1 = 1.0;
              * c21 = 1.0;
              * c22 = -1.0;
              * Undefined c31 factor;
              * Undefined c32 factor;
              * Undefined c33 factor;
              * Undefined c41 factor;
              * Undefined c421 factor;
              * Undefined c422 factor;
              * Undefined c43 factor;
              * Undefined c44 factor.
            >>> mod_conv = conv.get_modified(name="modified", c22=1, c33=-3)
            >>> mod_conv.summary()
            modified convention where
              * Bonds are counted multiple times in the sum;
              * Spin vectors are not normalized;
              * c1 = 1.0;
              * c21 = 1.0;
              * c22 = 1.0;
              * Undefined c31 factor;
              * Undefined c32 factor;
              * c33 = -3.0;
              * Undefined c41 factor;
              * Undefined c421 factor;
              * Undefined c422 factor;
              * Undefined c43 factor;
              * Undefined c44 factor.
        """

        if multiple_counting is None:
            multiple_counting = self._multiple_counting

        if spin_normalized is None:
            spin_normalized = self._spin_normalized

        if c1 is None:
            c1 = self._c1

        if c21 is None:
            c21 = self._c21

        if c22 is None:
            c22 = self._c22

        if c31 is None:
            c31 = self._c31

        if c32 is None:
            c32 = self._c32

        if c33 is None:
            c33 = self._c33

        if c41 is None:
            c41 = self._c41

        if c421 is None:
            c421 = self._c421

        if c422 is None:
            c422 = self._c422

        if c43 is None:
            c43 = self._c43

        if c44 is None:
            c44 = self._c44

        if name is None:
            name = self.name

        return Convention(
            spin_normalized=spin_normalized,
            multiple_counting=multiple_counting,
            c1=c1,
            c21=c21,
            c22=c22,
            c31=c31,
            c32=c32,
            c33=c33,
            c41=c41,
            c421=c421,
            c422=c422,
            c43=c43,
            c44=c44,
            name=name,
        )


# Populate __all__ with objects defined in this file
__all__ = list(set(dir()) - old_dir)
# Remove all semi-private objects
__all__ = [i for i in __all__ if not i.startswith("_")]
del old_dir
