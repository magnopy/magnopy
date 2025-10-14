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
import copy

from magnopy._constants._conventions import _SPINHAM_CONVENTIONS
from magnopy._exceptions import ConventionError

# Save local scope at this moment
old_dir = set(dir())
old_dir.add("old_dir")


class Convention:
    """
    Convention for the Hamiltonian.

    Parameters
    ----------
    name : str, default "custom"
        Name of the convention.
    spin_normalized : bool, optional
        Whether spins are normalized to unit vector or their actual spin value.
    multiple_counting : bool, optional
        Whether bonds are counted multiple times in the Hamiltonian.
    c1 : float, optional
        Coefficient for the one-spin term
    c21 : float, optional
        Coefficient for the (two spins & one site) term.
    c22 : float, optional
        Coefficient for the (two spins & two sites) term.
    c31 : float, optional
        Coefficient for the (three spins & one site) term.
    c32 : float, optional
        Coefficient for the (three spins & two sites) term.
    c33 : float, optional
        Coefficient for the (three spins & three sites) term.
    c41 : float, optional
        Coefficient for the (four spins & one site) term.
    c421 : float, optional
        Coefficient for the (four spins & two sites (1+3)) term.
    c422 : float, optional
        Coefficient for the (four spins & two sites (2+2)) term.
    c43 : float, optional
        Coefficient for the (four spins & three sites) term.
    c44 : float, optional
        Coefficient for the (four spins & four sites) term.

    Attributes
    ----------
    name : str
    spin_normalized : bool
    multiple_counting : bool
    c1 : float
    c21 : float
    c22 : float
    c31 : float
    c32 : float
    c33 : float
    c41 : float
    c421 : float
    c422 : float
    c43 : float
    c44 : float

    Examples
    --------

    To create a convention for some arbitrary Hamiltonian use

    .. doctest::

        >>> import magnopy
        >>> convention = magnopy.Convention(
        ...     spin_normalized=False,
        ...     multiple_counting=True,
        ...     c1=1,
        ...     c21=1,
        ...     c22=0.5,
        ...     c31=1,
        ...     c32=1,
        ...     c33=1,
        ...     c41=1,
        ...     c421=1,
        ...     c422=1,
        ...     c43=1,
        ...     c44=1,
        ... )
        >>> convention
        Convention(name='custom', spin_normalized=False, multiple_counting=True, c1=1, c21=1, c22=0.5, c31=1, c32=1, c33=1, c41=1, c421=1, c422=1, c43=1, c44=1)

    One can print the summary of the convention with

    .. doctest::

        >>> print(convention.summary(return_as_string=True))
        Convention: custom
          spin_normalized: False
          multiple_counting: True
          c1: 1
          c21: 1
          c22: 0.5
          c31: 1
          c32: 1
          c33: 1
          c41: 1
          c421: 1
          c422: 1
          c43: 1
          c44: 1
        >>> # Now, the following will also work:
        >>> print(convention)
        Convention: custom
          spin_normalized: False
          multiple_counting: True
          c1: 1
          c21: 1
          c22: 0.5
          c31: 1
          c32: 1
          c33: 1
          c41: 1
          c421: 1
          c422: 1
          c43: 1
          c44: 1

    You can define parameters partially

    .. doctest::

        >>> convention = magnopy.Convention(c22=-1)
        >>> print(convention.c22)
        -1
        >>> # The properties that are not defined can not be accessed
        >>> # convention.spin_normalized # doctest: +SKIP
        # Traceback (most recent call last):
        # ...
        # magnopy._exceptions.ConventionError: Convention of spin Hamiltonian has an undefined property 'spin_normalized':
        # Convention: custom
        #   spin_normalized: undefined
        #   multiple_counting: undefined
        #   c1: undefined
        #   c21: undefined
        #   c22: -1
        #   c31: undefined
        #   c32: undefined
        #   c33: undefined
        #   c41: undefined
        #   c421: undefined
        #   c422: undefined
        #   c43: undefined
        #   c44: undefined

    Magnopy also supports a few predefined conventions

    .. doctest::

        >>> convention = magnopy.Convention.get_predefined("tb2j")
        >>> print(convention.summary(return_as_string=True))
        Convention: tb2j
          spin_normalized: True
          multiple_counting: True
          c1: 1
          c21: -1
          c22: -1
          c31: 1
          c32: 1
          c33: 1
          c41: 1
          c421: 1
          c422: 1
          c43: 1
          c44: 1
    """

    def __init__(
        self,
        name="custom",
        spin_normalized=None,
        multiple_counting=None,
        c1=None,
        c21=None,
        c22=None,
        c31=None,
        c32=None,
        c33=None,
        c41=None,
        c421=None,
        c422=None,
        c43=None,
        c44=None,
    ):
        self._name = name.lower()
        self._spin_normalized = spin_normalized
        self._multiple_counting = multiple_counting
        self._c1 = c1
        self._c21 = c21
        self._c22 = c22
        self._c31 = c31
        self._c32 = c32
        self._c33 = c33
        self._c41 = c41
        self._c421 = c421
        self._c422 = c422
        self._c43 = c43
        self._c44 = c44

    def __eq__(self, other):
        if not isinstance(other, type(self)):
            return NotImplemented

        if self.name != other.name:
            return False

        params = [
            "spin_normalized",
            "multiple_counting",
            "c1",
            "c21",
            "c22",
            "c31",
            "c32",
            "c33",
            "c41",
            "c421",
            "c422",
            "c43",
            "c44",
        ]
        for p in params:
            if getattr(self, f"_{p}") != getattr(other, f"_{p}"):
                return False
        return True

    def __repr__(self):
        repr_str = f"Convention(name='{self._name}'"
        attrs = [
            ("spin_normalized", self._spin_normalized),
            ("multiple_counting", self._multiple_counting),
            ("c1", self._c1),
            ("c21", self._c21),
            ("c22", self._c22),
            ("c31", self._c31),
            ("c32", self._c32),
            ("c33", self._c33),
            ("c41", self._c41),
            ("c421", self._c421),
            ("c422", self._c422),
            ("c43", self._c43),
            ("c44", self._c44),
        ]
        for name, value in attrs:
            if value is not None:
                repr_str += f", {name}={value}"
        repr_str += ")"
        return repr_str

    def __str__(self):
        """
        Returns a human-readable string representation of the convention.
        """
        return self.summary(return_as_string=True)

    def summary(self, return_as_string=False):
        """
        Returns a summary of the convention.

        Parameters
        ----------
        return_as_string : bool, default False
            If ``True``, then a string is returned. Otherwise, it is printed to
            the console.

        Returns
        -------
        summary : str, optional
        """
        summary_lines = [f"Convention: {self.name}"]
        for attr in [
            "spin_normalized",
            "multiple_counting",
            "c1",
            "c21",
            "c22",
            "c31",
            "c32",
            "c33",
            "c41",
            "c421",
            "c422",
            "c43",
            "c44",
        ]:
            value = getattr(self, f"_{attr}")
            summary_lines.append(
                f"  {attr}: {value if value is not None else 'undefined'}"
            )
        summary_str = "\n".join(summary_lines)

        if return_as_string:
            return summary_str
        print(summary_str)

    @classmethod
    def get_predefined(cls, name: str):
        """
        Returns one of the predefined conventions.

        Parameters
        ----------
        name : str
            Name of the convention. Possible values are:

            * ``tb2j``
            * ``vampire``
            * ``grogu``
            * ``spinw``

        Returns
        -------
        convention : :py:class:`~.Convention`
            Instance of the class :py:class:`~.Convention` with the predefined values.

        Raises
        ------
        ValueError
            If given name is not supported.
        """

        name = name.lower()

        if name not in _SPINHAM_CONVENTIONS:
            raise ValueError(
                f'Given name "{name}" is not in the list of supported conventions: '
                + ", ".join(list(_SPINHAM_CONVENTIONS))
            )

        kwargs = dict(name=name, c1=1, c31=1, c32=1, c33=1, c41=1, c421=1, c422=1, c43=1, c44=1)
        kwargs.update(_SPINHAM_CONVENTIONS[name])
        return cls(**kwargs)

    def get_modified(self, **kwargs):
        """
        Returns a modified version of the current convention.

        Parameters
        ----------
        name : str, optional
            Name of the convention.
        spin_normalized : bool, optional
            Whether spins are normalized to unit vector or their actual spin value.
        multiple_counting : bool, optional
            Whether bonds are counted multiple times in the Hamiltonian.
        c1 : float, optional
            Coefficient for the one-spin term
        c21 : float, optional
            Coefficient for the (two spins & one site) term.
        c22 : float, optional
            Coefficient for the (two spins & two sites) term.
        c31 : float, optional
            Coefficient for the (three spins & one site) term.
        c32 : float, optional
            Coefficient for the (three spins & two sites) term.
        c33 : float, optional
            Coefficient for the (three spins & three sites) term.
        c41 : float, optional
            Coefficient for the (four spins & one site) term.
        c421 : float, optional
            Coefficient for the (four spins & two sites (1+3)) term.
        c422 : float, optional
            Coefficient for the (four spins & two sites (2+2)) term.
        c43 : float, optional
            Coefficient for the (four spins & three sites) term.
        c44 : float, optional
            Coefficient for the (four spins & four sites) term.

        Returns
        -------
        convention : :py:class:`~.Convention`
            Instance of the class :py:class:`~.Convention` with the requested modifications.

        Examples
        --------

        To create a modified version of some existing convention use

        .. doctest::

            >>> import magnopy
            >>> convention = magnopy.Convention.get_predefined("tb2j")
            >>> convention_modified = convention.get_modified(c21=0.5)
            >>> print(convention_modified.summary(return_as_string=True))
            Convention: tb2j
              spin_normalized: True
              multiple_counting: True
              c1: 1
              c21: 0.5
              c22: -1
              c31: 1
              c32: 1
              c33: 1
              c41: 1
              c421: 1
              c422: 1
              c43: 1
              c44: 1
        """
        new_convention = copy.copy(self)
        for key, value in kwargs.items():
            setattr(new_convention, f"_{key}", value)
        return new_convention

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, name):
        self._name = name.lower()

    @property
    def spin_normalized(self):
        if self._spin_normalized is None:
            raise ConventionError(self, "spin_normalized")
        return self._spin_normalized

    @spin_normalized.setter
    def spin_normalized(self, value):
        raise AttributeError(f"{self.__class__.__name__}'s spin_normalized is immutable")

    @property
    def multiple_counting(self):
        if self._multiple_counting is None:
            raise ConventionError(self, "multiple_counting")
        return self._multiple_counting

    @multiple_counting.setter
    def multiple_counting(self, value):
        raise AttributeError(f"{self.__class__.__name__}'s multiple_counting is immutable")

    @property
    def c1(self):
        if self._c1 is None:
            raise ConventionError(self, "c1")
        return self._c1

    @c1.setter
    def c1(self, value):
        raise AttributeError(f"{self.__class__.__name__}'s c1 is immutable")

    @property
    def c21(self):
        if self._c21 is None:
            raise ConventionError(self, "c21")
        return self._c21

    @c21.setter
    def c21(self, value):
        raise AttributeError(f"{self.__class__.__name__}'s c21 is immutable")

    @property
    def c22(self):
        if self._c22 is None:
            raise ConventionError(self, "c22")
        return self._c22

    @c22.setter
    def c22(self, value):
        raise AttributeError(f"{self.__class__.__name__}'s c22 is immutable")

    @property
    def c31(self):
        if self._c31 is None:
            raise ConventionError(self, "c31")
        return self._c31

    @c31.setter
    def c31(self, value):
        raise AttributeError(f"{self.__class__.__name__}'s c31 is immutable")

    @property
    def c32(self):
        if self._c32 is None:
            raise ConventionError(self, "c32")
        return self._c32

    @c32.setter
    def c32(self, value):
        raise AttributeError(f"{self.__class__.__name__}'s c32 is immutable")

    @property
    def c33(self):
        if self._c33 is None:
            raise ConventionError(self, "c33")
        return self._c33

    @c33.setter
    def c33(self, value):
        raise AttributeError(f"{self.__class__.__name__}'s c33 is immutable")

    @property
    def c41(self):
        if self._c41 is None:
            raise ConventionError(self, "c41")
        return self._c41

    @c41.setter
    def c41(self, value):
        raise AttributeError(f"{self.__class__.__name__}'s c41 is immutable")

    @property
    def c421(self):
        if self._c421 is None:
            raise ConventionError(self, "c421")
        return self._c421

    @c421.setter
    def c421(self, value):
        raise AttributeError(f"{self.__class__.__name__}'s c421 is immutable")

    @property
    def c422(self):
        if self._c422 is None:
            raise ConventionError(self, "c422")
        return self._c422

    @c422.setter
    def c422(self, value):
        raise AttributeError(f"{self.__class__.__name__}'s c422 is immutable")

    @property
    def c43(self):
        if self._c43 is None:
            raise ConventionError(self, "c43")
        return self._c43

    @c43.setter
    def c43(self, value):
        raise AttributeError(f"{self.__class__.__name__}'s c43 is immutable")

    @property
    def c44(self):
        if self._c44 is None:
            raise ConventionError(self, "c44")
        return self._c44

    @c44.setter
    def c44(self, value):
        raise AttributeError(f"{self.__class__.__name__}'s c44 is immutable")


# Populate __all__ with objects defined in this file
__all__ = list(set(dir()) - old_dir)
# Remove all semi-private objects
__all__ = [i for i in __all__ if not i.startswith("_")]
del old_dir