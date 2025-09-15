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


from multiprocessing import Pool

# Save local scope at this moment
old_dir = set(dir())
old_dir.add("old_dir")


def multiprocess_over_k(kpoints, function, relative=False, number_processors=None):
    r"""
    Parallelize calculation over the kpoints using |multiprocessing|_ module.

    Parameters
    ----------
    kpoints : (N, 3) |array-like|_
        List of the kpoints.
    function : callable
        Function that process one kpoint and will be called as

        .. code-block:: python

            result = function(kpoints[i], relative)

    relative : bool, default False
        If ``relative=True``, then ``k`` is interpreted as given relative to the
        reciprocal unit cell. Otherwise it is interpreted as given in absolute
        coordinates.
    number_processors : int, optional
        By default magnopy uses all available processes. Pass ``number_processors=1`` to
        run in serial.

    Returns
    -------
    results : (N, ) list
        List of objects that are returned by the ``function``.

    Notes
    -----

    When using this function of magnopy in your Python scripts make sure to safeguard
    your script with the

    .. code-block:: python

        import magnopy

        # Import more stuff
        # or
        # Define your functions, classes

        if __name__ == "__main__":

            # Write your executable code here

    For more information refer to the  "Safe importing of main module" section in
    |multiprocessing|_ docs.
    """

    relative = [relative for _ in kpoints]

    if number_processors == 1:
        results = list(map(function, kpoints, relative))
    else:
        with Pool(number_processors) as p:
            results = p.starmap(function, zip(kpoints, relative))

    return results


# Populate __all__ with objects defined in this file
__all__ = list(set(dir()) - old_dir)
# Remove all semi-private objects
__all__ = [i for i in __all__ if not i.startswith("_")]
del old_dir
