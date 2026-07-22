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

from dataclasses import fields
import numpy as np
import pytest

from magnopy._spinham._crystal import _Crystal

################################################################################
#                      Configuration and helper functions                      #
################################################################################


def make(**overrides):
    r"""
    Two-atom crystal. Both atoms are magnetic.
    """

    kwargs = dict(
        cell=np.eye(3),
        names=["Fe", "Fe"],
        positions=[[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]],
        spins=[2.5, 2.5],
        g_factors=[2.0, 2.0],
        magnetic=[True, True],
    )
    kwargs.update(overrides)
    return _Crystal(**kwargs)


FLOAT_ARRAY_FIELDS = ("cell", "positions", "spins", "g_factors")
BOOL_ARRAY_FIELDS = ("magnetic",)
ARRAY_FIELDS = FLOAT_ARRAY_FIELDS + BOOL_ARRAY_FIELDS
PER_ATOM_FIELDS = ("positions", "spins", "g_factors", "magnetic")
ALL_FIELDS = ARRAY_FIELDS + ("names",)


################################################################################
#                                _as_array tests                               #
################################################################################
def test_non_finite_rejected_after_coercion():
    # isfinite must run on the coerced array; on an object-type array it raises
    # TypeError instead of validation
    with pytest.raises(ValueError, match="non-finite"):
        make(spins=np.array([2.5, np.inf], dtype=object))


################################################################################
#                                  Field tests                                 #
################################################################################


def test_field_constants_match_dataclass():
    assert set(ALL_FIELDS) == {_.name for _ in fields(_Crystal)}


################################################################################
#                               Validation tests                               #
#                                                                              #
# One test per clause of raises. These tests validate that the incorrect       #
# user's input is handled correctly. They do not test whether the              #
# constructor processes them correctly into store attributes.                  #
################################################################################


def test_singular_cell_rejected():
    with pytest.raises(ValueError, match="singular"):
        make(cell=[[1, 0, 0], [2, 0, 0], [0, 0, 1]])


def test_near_singular_cell_rejected():
    with pytest.raises(ValueError, match="singular"):
        make(cell=[[1, 0, 0], [1, 1e-14, 0], [0, 0, 1]])


def test_singularity_check_is_scale_free():
    for a in [1e-3, 1e-1, 1.0, 1e2, 1e4]:
        make(cell=a * np.eye(3))  # must not raise


def test_non_str_names_rejected():
    with pytest.raises(ValueError, match="str"):
        make(names=(1, 2))


@pytest.mark.parametrize("field", FLOAT_ARRAY_FIELDS)
@pytest.mark.parametrize("bad_value", (np.nan, np.inf, -np.inf))
def test_non_finite_rejected(field, bad_value):
    good_crystal = make()
    arr = np.array(getattr(good_crystal, field), dtype=float)
    arr.flat[0] = bad_value
    with pytest.raises(ValueError, match="non-finite"):
        make(**{field: arr})


def test_negative_spin_rejected():
    with pytest.raises(ValueError, match="spin"):
        make(spins=[2.5, -1.0])


def test_zero_spin_allowed():
    make(spins=[2.5, 0.0])


def test_g_factor_rejected():
    with pytest.raises(ValueError, match="g.factor"):
        make(g_factors=[2.0, 0.0])


def test_no_magnetic_atom_rejected():
    with pytest.raises(ValueError, match="magnetic"):
        make(magnetic=[False, False])


def test_empty_crystal_rejected():
    with pytest.raises(ValueError, match="empty"):
        _Crystal(
            cell=np.eye(3),
            names=(),
            positions=np.zeros((0, 3)),
            spins=np.zeros(0),
            g_factors=np.zeros(0),
            magnetic=np.zeros(0, dtype=bool),
        )


@pytest.mark.parametrize("field", PER_ATOM_FIELDS)
def test_length_mismatch_rejected(field):
    good_crystal = make()
    truncated_value = np.asarray(getattr(good_crystal, field))[:-1]
    with pytest.raises(ValueError):
        make(**{field: truncated_value})


@pytest.mark.parametrize(
    "field, bad_value",
    (
        ("cell", np.eye(2)),
        ("cell", np.ones((3, 3, 3))),
        ("positions", [[0.0, 0.0], [0.5, 0.5]]),  # (N, 2)
        ("spins", [[2.5], [2.5]]),  # (N, 1)
        ("g_factors", [[2.0, 2.0]]),  # (1, N)
        ("magnetic", [[True], [True]]),  # (N, 1)
    ),
)
def test_wrong_shape_rejected(field, bad_value):
    with pytest.raises(ValueError):
        make(**{field: bad_value})


def test_magnetic_accepts_zero_one():
    crystal = make(magnetic=[1, 0])
    assert crystal.magnetic.dtype == bool
    assert crystal.magnetic.tolist() == [True, False]


def test_magnetic_other_integers_rejected():
    with pytest.raises(ValueError):
        make(magnetic=[1, 2])


def test_magnetic_near_one_rejected():
    # dtype=int would truncate 1.000042 -> 1 before validation could see it
    with pytest.raises(ValueError):
        make(magnetic=[1.0, 1.000042])


def test_magnetic_accepts_negative_zero():
    assert make(magnetic=[-0.0, 1.0]).magnetic.tolist() == [False, True]


def test_cell_not_array_like_raises_value_error():
    with pytest.raises(ValueError):
        make(cell="not a cell")


def test_ragged_positions_rejected():
    with pytest.raises(ValueError):
        make(positions=[[0, 0, 0], [1, 1]])


@pytest.mark.parametrize("missing", ALL_FIELDS)
def test_missing_field_is_type_error(missing):
    # No defaults in _Crystal.
    kwargs = dict(
        cell=np.eye(3),
        names=("Fe", "Fe"),
        positions=[[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]],
        spins=[2.5, 2.5],
        g_factors=[2.0, 2.0],
        magnetic=[True, True],
    )
    del kwargs[missing]
    with pytest.raises(TypeError):
        _Crystal(**kwargs)


################################################################################
#                                Storage types                                 #
################################################################################


def test_names_stored_as_tuple_of_str():
    crystal = make(names=["Fe", "Fe"])
    assert isinstance(crystal.names, tuple)
    assert all(isinstance(_, str) for _ in crystal.names)


def test_names_accepts_any_iterable():
    # An iterator is one-shot: an implementation that measures `names` before
    # consuming it, or walks it twice, silently yields an empty crystal.
    # zip/map/filter and generator comprehensions are all iterators.
    assert make(names=iter(["Fe", "Fe"])).names == ("Fe", "Fe")


@pytest.mark.parametrize("field", FLOAT_ARRAY_FIELDS)
def test_float_fields_are_float64_ndarray(field):
    crystal = make(spins=[2, 3], g_factors=[2, 2])  # ints on input
    arr = getattr(crystal, field)
    assert isinstance(arr, np.ndarray) and arr.dtype == np.float64


@pytest.mark.parametrize("field", BOOL_ARRAY_FIELDS)
def test_bool_fields_are_bool_ndarrays(field):
    crystal = make()
    arr = getattr(crystal, field)
    assert isinstance(arr, np.ndarray)
    assert arr.dtype == bool


@pytest.mark.parametrize("data", ([1, 0], [True, False], [1.0, 0.0]))
def test_magnetic_stored_as_bool_ndarrays(data):
    crystal = make(magnetic=data)
    arr = getattr(crystal, "magnetic")
    assert isinstance(arr, np.ndarray)
    assert arr.dtype == bool


@pytest.mark.parametrize("seq", (list, tuple))
def test_any_sequence_type_is_accepted(seq):
    _Crystal(
        names=seq(["Fe", "Fe"]),
        cell=seq([seq([1, 0, 0]), seq([0, 1, 0]), seq([0, 0, 1])]),
        positions=seq([seq([0.0, 0.0, 0.0]), seq([0.5, 0.5, 0.5])]),
        spins=seq([2.5, 2.5]),
        g_factors=seq([2.0, 2.0]),
        magnetic=seq([True, True]),
    )


def test_mixed_containers_are_accepted():
    _Crystal(
        names=iter(["Fe", "Fe"]),
        cell=np.eye(3),
        positions=[(0.0, 0.0, 0.0), np.array([0.5, 0.5, 0.5])],
        spins=(2.5, 2.5),
        g_factors=np.array([2.0, 2.0]),
        magnetic=[True, False],
    )


################################################################################
#                                 Immutability                                 #
################################################################################


def test_stored_names_are_read_only():
    crystal = make()
    with pytest.raises(TypeError):
        crystal.names[0] = "New name"


@pytest.mark.parametrize("field", ARRAY_FIELDS)
def test_stored_arrays_are_read_only(field):
    crystal = make()
    with pytest.raises(ValueError):
        getattr(crystal, field)[0] = 0


@pytest.mark.parametrize("field", ALL_FIELDS)
def test_fields_cannot_be_rebound(field):
    crystal = make()
    with pytest.raises(AttributeError):
        setattr(crystal, field, "new value")


def test_caller_array_stays_writable():
    positions = np.array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]])
    make(positions=positions)
    positions[0, 0] = 0.42  # must not raise


def test_constructor_copies_its_arguments():
    # copying is what makes freeze safe to apply
    positions = np.array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]])
    crystal = _Crystal(
        cell=np.eye(3),
        names=("Fe", "Fe"),
        positions=positions,
        spins=[2.5, 2.5],
        g_factors=[2.0, 2.0],
        magnetic=[True, True],
    )
    positions[0, 0] = 0.42  # user's array still writeable
    assert crystal.positions[0, 0] == 0.0  # crystal is unaffected


################################################################################
#                                  Identities                                  #
################################################################################


def test_crystal_equals_itself():
    crystal = make()
    assert crystal == crystal
    assert hash(crystal) == hash(crystal)


def test_equal_crystals_are_equal_and_hash_alike():
    c1, c2 = make(), make()
    assert c1 == c2
    assert hash(c1) == hash(c2)


def test_crystal_is_usable_as_dict_key():
    d = {make(): "value"}
    assert d[make()] == "value"


@pytest.mark.parametrize(
    "field, value",
    (
        ("cell", 2 * np.eye(3)),
        ("names", ("Fe", "Ni")),
        ("positions", [[0.0, 0.0, 0.0], [0.25, 0.25, 0.25]]),
        ("spins", [2.5, 1.5]),
        ("g_factors", [2.0, 2.1]),
        ("magnetic", [True, False]),
    ),
)
def test_any_core_field_change_breaks_equality(field, value):
    assert make() != make(**{field: value})


def test_names_are_part_of_the_crystal():
    assert make(names=("Fe", "Fe")) != make(names=("Fe", "Ni"))


def test_equality_is_exact_not_tolerant():
    assert make() != make(
        positions=[[0.0, 0.0, 0.0], [0.5 + np.spacing(0.5), 0.5, 0.5]]
    )


def test_negative_zero_hashes_and_compares_as_positive_zero():
    c1 = make(positions=[[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]])
    c2 = make(positions=[[-0.0, 0.0, 0.0], [0.5, 0.5, 0.5]])
    assert c1 == c2
    assert hash(c1) == hash(c2)
    assert c1.positions.tobytes() == c2.positions.tobytes()


def test_eq_with_non_crystal_returns_not_implemented():
    assert make().__eq__(object()) is NotImplemented
    assert make() != 42


################################################################################
#                                   is_close                                   #
################################################################################


def test_is_close_accepts_last_bit_differences():
    x, i, L = 0.3, 1, 3

    route_1 = (x + i) / L
    route_2 = x / L + i / L
    assert route_1 != route_2  # Guard the test's assumption

    c1 = make(positions=[[0.0, 0.0, 0.0], [route_1, 0.5, 0.5]])
    c2 = make(positions=[[0.0, 0.0, 0.0], [route_2, 0.5, 0.5]])
    assert c1 != c2  # exact comparison
    assert c1.is_close(c2)  # tolerant comparison


def test_is_close_wraps_fractional_positions():
    c1 = make(positions=[[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]])
    c2 = make(positions=[[1.0, 0.0, 0.0], [0.5, 0.5, 0.5]])
    assert c1.is_close(c2)


def test_is_close_rejects_different_names():
    # is_close is tolerant about floats, not about what the atoms are
    assert not make(names=("Fe", "Fe")).is_close(make(names=("Fe", "Ni")))


def test_is_close_rejects_different_magnetic_sites():
    assert not make(magnetic=[True, True]).is_close(make(magnetic=[True, False]))


def test_is_close_rejects_physically_different_positions():
    assert not make(positions=[[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]).is_close(
        make(positions=[[0.0, 0.0, 0.0], [0.25, 0.25, 0.25]])
    )


def test_is_close_rejects_different_g_factors():
    assert not make().is_close(make(g_factors=[2.0, 2.5]))


def test_is_close_rejects_different_spins():
    assert not make().is_close(make(spins=[2.5, 1.5]))


def test_is_close_rejects_different_length():
    # Needs to be handled before any allclose, which would raise on shape mismatch
    c_3 = make(
        names=("Fe", "Fe", "Fe"),
        positions=[[0.0, 0.0, 0.0], [0.5, 0.5, 0.5], [0.25, 0.25, 0.25]],
        spins=[2.5, 2.5, 2.5],
        g_factors=[2.0, 2.0, 2.0],
        magnetic=[True, True, True],
    )
    assert not make().is_close(c_3)


def test_is_close_is_reflexive_and_symmetric():
    c1 = make(positions=[[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]])
    c2 = make(positions=[[0.0, 0.0, 0.0], [0.5 + np.spacing(0.5), 0.5, 0.5]])
    assert c1.is_close(c1)
    assert c1.is_close(c2) == c2.is_close(c1)


def test_is_close_is_not_transitive():
    # If this test ever fails, someone has made is_close exact -- and the fix is
    # NOT to delete this test

    atol = 1e-12
    d = 0.9 * atol

    c1 = make(positions=[[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]])
    c2 = make(positions=[[d, 0.0, 0.0], [0.5, 0.5, 0.5]])
    c3 = make(positions=[[2 * d, 0.0, 0.0], [0.5, 0.5, 0.5]])

    assert c1.is_close(c2, atol=atol)
    assert c2.is_close(c3, atol=atol)
    assert not c1.is_close(c3, atol=atol)
