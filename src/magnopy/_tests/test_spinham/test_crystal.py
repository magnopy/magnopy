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
from inspect import signature
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


def make_with_ligand(**overrides):
    # The ligand carries non-zero spin on purpose
    kwargs = dict(
        cell=np.eye(3),
        names=["Fe", "Fe", "O"],
        positions=[[0.0, 0.0, 0.0], [0.5, 0.5, 0.0], [0.5, 0.0, 0.5]],
        spins=[2.5, 2.5, 0.1],
        g_factors=[2.0, 2.0, 2.0],
        magnetic=[True, True, False],
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


def test_is_close_uses_atol_boundary():
    # Just under and just over: pins that atol is actually applied to positions
    atol = 1e-9
    under = make(positions=[[0.9 * atol, 0.0, 0.0], [0.5, 0.5, 0.5]])
    over = make(positions=[[1.1 * atol, 0.0, 0.0], [0.5, 0.5, 0.5]])
    assert make().is_close(under, atol=atol)
    assert not make().is_close(over, atol=atol)


def test_is_close_atol_applies_to_spins_and_g_factors():
    atol = 1e-6
    assert make().is_close(make(spins=[2.5, 2.5 + 0.9 * atol]), atol=atol)
    assert not make().is_close(make(spins=[2.5, 2.5 + 1.1 * atol]), atol=atol)
    assert make().is_close(make(g_factors=[2.0, 2.0 + 0.9 * atol]), atol=atol)
    assert not make().is_close(make(g_factors=[2.0, 2.0 + 1.1 * atol]), atol=atol)


def test_is_close_cell_uses_relative_tolerance():
    # cell has physical dimensions and magnopy does not make strict assumption about
    # its units. Thus, the same relative discrepancy must be accepted (or rejected)
    # whether the cell is given in angstrom, nanometer or else. An absolute tolerance
    # can not achieve that.
    rtol = 1e-11
    for a in (1.0, 1e-3, 1e3, 1e-10):
        crystal1 = make(cell=a * np.eye(3))
        crystal2 = make(cell=a * (1.0 + rtol) * np.eye(3))
        assert crystal1.is_close(crystal2, cell_rtol=10 * rtol), (
            f"Failed with scale {a} (close)"
        )
        assert not crystal1.is_close(crystal2, cell_rtol=0.1 * rtol), (
            f"Failed with scale {a} (not close)"
        )


def test_is_close_cell_ignores_atol():
    a = 1e6
    c1 = make(cell=a * np.eye(3))
    c2 = make(cell=(a + 1e-3) * np.eye(3))
    assert c1.is_close(c2, atol=1e-8, cell_rtol=1e-8)


def test_is_close_wraps_only_integer_offsets():
    # 0.4, 0.6 are not the same site
    # guard against over-eager wrappers
    c1 = make(positions=[[0.0, 0.0, 0.0], [0.4, 0.5, 0.5]])
    c2 = make(positions=[[0.0, 0.0, 0.0], [0.6, 0.5, 0.5]])
    assert not c1.is_close(c2)


def test_is_close_wraps_boundary_half_cell():
    # Whatever the implementation of wrapper does it must not call these two close
    c1 = make(positions=[[0.0, 0.0, 0.0], [0.0, 0.5, 0.5]])
    c2 = make(positions=[[0.5, 0.0, 0.0], [0.0, 0.5, 0.5]])
    assert not c1.is_close(c2)


def test_is_close_with_non_crystal():
    with pytest.raises(TypeError):
        make().is_close(42)


################################################################################
#                                     diff                                     #
################################################################################


# iff == "if and only if"
def test_diff_empty_iff_is_close():
    assert make().diff(make()) == ()
    assert make().diff(make(spins=[2.5, 1.5])) != ()


@pytest.mark.parametrize(
    "overrides",
    (
        {},
        {"g_factors": [2.0, 2.1]},
        {"names": ("Fe", "Ni")},
    ),
)
def test_diff_and_is_close_agree(overrides):
    other = make(**overrides)
    assert make().is_close(other) == (make().diff(other) == ())


@pytest.mark.parametrize(
    "field, value, mentioned",
    (
        ("cell", 2 * np.eye(3), "cell"),
        ("names", ("Fe", "Ni"), "name"),
        ("positions", [[0.0, 0.0, 0.0], [0.25, 0.25, 0.25]], "position"),
        ("spins", [2.5, 1.5], "spin"),
        ("g_factors", [2.0, 2.1], "factor"),
        ("magnetic", [True, False], "magnetic"),
    ),
)
def test_diff_names_the_mismatched_field(field, value, mentioned):
    report = " ".join(make().diff(make(**{field: value})))
    assert mentioned in report.lower()


def test_diff_reports_every_different_field():
    report = make().diff(make(spins=[2.5, 1.5], g_factors=[2.0, 2.5]))
    joined = " ".join(report).lower()
    assert "spin" in joined and "factor" in joined


def test_diff_on_different_length():
    c_3 = make(
        names=("Fe", "Fe", "Fe"),
        positions=[[0.0, 0.0, 0.0], [0.5, 0.5, 0.5], [0.25, 0.25, 0.25]],
        spins=[2.5, 2.5, 2.5],
        g_factors=[2.0, 2.0, 2.0],
        magnetic=[True, True, True],
    )
    report = make().diff(c_3)
    assert report != ()


def test_diff_with_non_crystal_raises():
    with pytest.raises(TypeError):
        make().diff(42)


def test_is_close_and_diff_share_a_signature():
    # is_close forwards to diff, so their signature shall not drift apart
    assert (
        signature(_Crystal.is_close).parameters == signature(_Crystal.diff).parameters
    )


def test_diff_returns_a_tuple_of_strings():
    c_3 = make(
        names=("Fe", "Fe", "Fe"),
        positions=[[0.0, 0.0, 0.0], [0.5, 0.5, 0.5], [0.25, 0.25, 0.25]],
        spins=[2.5, 2.5, 2.5],
        g_factors=[2.0, 2.0, 2.0],
        magnetic=[True, True, True],
    )

    # Gurads the lenght ismatch branch
    report = make().diff(c_3)
    assert isinstance(report, tuple)
    assert all(isinstance(line, str) for line in report)


################################################################################
#                                   select()                                   #
################################################################################


def test_select_with_boolean_mask():
    crystal = make_with_ligand()
    subset = crystal.select(crystal.magnetic)
    assert subset.names == ("Fe", "Fe")
    np.testing.assert_equal(subset.spins, [2.5, 2.5])


def tes_select_with_integer_array():
    crystal = make_with_ligand()
    subset = crystal.select(np.array([2, 0]))
    assert subset.names == ("O", "Fe")
    np.testing.assert_equa(subset.positions, crystal.positions[[2, 0]])


def test_select_with_list_of_integers():
    assert make_with_ligand().select([1, 2]).names == ("Fe", "O")


def test_select_with_slice():
    assert make_with_ligand().select(slice(0, 2)).names == ("Fe", "Fe")


def test_select_reorders():
    # index array is a permutation, not just a filter
    crystal = make_with_ligand()
    subset = crystal.select([2, 1, 0])
    assert subset.names == ("O", "Fe", "Fe")
    np.testing.assert_equal(subset.spins, [0.1, 2.5, 2.5])


def test_select_can_repeat_indices():
    crystal = make_with_ligand()
    tiled = crystal.select([0, 2, 0, 2, 0, 2])
    assert len(tiled) == 6
    assert tiled.names == ("Fe", "O") * 3


def test_select_identity_reproduces_the_crystal():
    crystal = make()
    assert crystal.select(np.arange(len(crystal))) == crystal


def test_select_returns_a_crystal():
    assert isinstance(make_with_ligand().select([0, 1]), _Crystal)


def test_select_keeps_the_cell():
    # Subset of atoms still leaves in the same cell.
    # If the caller needs the supercell - its their job to construct one explicitly
    crystal = make_with_ligand(cell=2 * np.eye(3))
    np.testing.assert_equal(crystal.cell, crystal.select([0, 1]).cell)


def test_select_carries_every_field():
    crystal = make_with_ligand()
    subset = crystal.select([2, 0])
    for name in ("spins", "g_factors", "magnetic", "positions"):
        subset_field = getattr(subset, name)
        original_field = getattr(crystal, name)[[2, 0]]
        np.testing.assert_equal(subset_field, original_field)


def test_select_result_is_frozen():
    subset = make_with_ligand().select([0, 1])
    with pytest.raises(ValueError):
        subset.positions[0] = 0


def test_select_revalidates():
    # Every construction shall go through __post__init__
    # Selecting only the ligand leaves no magnetic atoms
    with pytest.raises(ValueError, match="magnetic"):
        make_with_ligand().select([2])


def test_select_empty_rejected():
    with pytest.raises(ValueError, match="empty"):
        make_with_ligand().select([])


def select_scalar_index_rejected():
    # A scalar would invite for `for i in range(len(c)): c.select(i)`, which is
    # not intended use of select.
    with pytest.raises(TypeError):
        make_with_ligand().select(0)


def select_out_of_range_rejected():
    with pytest.raises(IndexError):
        make_with_ligand().select([0, 99])


def test_select_boolean_mask_of_wrong_length_rejected():
    # Its ok for indices, but not for a mask
    with pytest.raises(ValueError):
        make_with_ligand().select(np.array([True, False]))


def test_crystal_is_not_a_sequence():
    # __getitem__ was replaced by select precisely so that __len__ alone does
    # not make Python treat a crystal as a sequence of atoms. Without this,
    # iter(), list(), `in` and np.asarray() all fall back to crystal[0],
    # crystal[1], ... and fail with an incomprehensible message.
    crystal = make_with_ligand()
    with pytest.raises(TypeError):
        list(crystal)
    with pytest.raises(TypeError):
        crystal[0]


################################################################################
#                              Derived properties                              #
################################################################################


def test_M_prie_counts_all_atoms():
    assert make().M_prime == 2
    assert make_with_ligand().M_prime == 3


@pytest.mark.parametrize("func", (make, make_with_ligand))
def test_len_is_M_prime(func):
    crystal = func()
    assert len(crystal) == crystal.M_prime


def test_M_counts_magnetic_atoms():
    assert make().M == 2
    assert make_with_ligand().M == 2


def test_M_ignores_spin_magnitude():
    # A magnetic site with S=0 still counts;
    # A non-magnetic site with S!=0 stil does not.
    assert make(spins=[2.5, 0.0]).M == 2
    assert make_with_ligand(magnetic=[True, False, False]).M == 1


@pytest.mark.parametrize("name", ("map_to_all", "map_to_magnetic"))
def test_maps_are_integer_array(name):
    arr = getattr(make_with_ligand(), name)
    assert isinstance(arr, np.ndarray)
    assert arr.dtype.kind == "i"


def test_map_to_all_lists_magnetic_atom_indices():
    np.testing.assert_equal(make_with_ligand().map_to_all, [0, 1])
    np.testing.assert_equal(make().map_to_all, [0, 1])


def test_map_to_all_is_sorted():
    crystal = make_with_ligand(magnetic=[False, True, True])
    np.testing.assert_equal(crystal.map_to_all, [1, 2])


@pytest.mark.parametrize("func", (make, make_with_ligand))
def test_maps_has_correct_lengths(func):
    crystal = func()
    assert len(crystal.map_to_magnetic) == crystal.M_prime
    assert len(crystal.map_to_all) == crystal.M


def test_map_to_magnetic_gives_positions_among_magnetic():
    np.testing.assert_equal(make_with_ligand().map_to_magnetic, [0, 1, -1])


def test_map_to_magnetic_is_minus_one_for_non_magnetic():
    np.testing.assert_equal(
        make_with_ligand(magnetic=[False, True, False]).map_to_magnetic,
        [-1, 0, -1],
    )


@pytest.mark.parametrize(
    "func, overrides",
    (
        (make, {}),
        (make_with_ligand, {}),
        (make_with_ligand, {"magnetic": [False, True, True]}),
    ),
)
def test_maps_are_mutual_inverses(func, overrides):
    crystal = func(**overrides)
    np.testing.assert_equal(
        crystal.map_to_magnetic[crystal.map_to_all], np.arange(crystal.M)
    )


@pytest.mark.parametrize("name", ("map_to_all", "map_to_magnetic"))
def test_maps_are_read_only(name):
    arr = getattr(make_with_ligand(), name)
    with pytest.raises(ValueError):
        arr[0] = 0


@pytest.mark.parametrize("name", ("map_to_all", "map_to_magnetic"))
def test_maps_are_cached(name):
    crystal = make_with_ligand()
    assert getattr(crystal, name) is getattr(crystal, name)


@pytest.mark.parametrize("name", ("map_to_all", "map_to_magnetic"))
def test_all_magnetic_crystal_has_identity_maps(name):
    # degenerate case, when all atoms are magnetic
    crystal = make()
    np.testing.assert_equal(getattr(crystal, name), np.arange(len(crystal)))


def test_single_magnetic_atom():
    crystal = make_with_ligand(magnetic=[False, False, True])
    assert crystal.M == 1 and crystal.M_prime == 3
    np.testing.assert_equal(crystal.map_to_all, [2])
    np.testing.assert_equal(crystal.map_to_magnetic, [-1, -1, 0])


################################################################################
#                                   __repr__                                   #
################################################################################


def test_repr_is_a_single_line():
    assert "\n" not in repr(make())


def test_repr_does_not_dump_arrays():
    crystal = make(
        names=["Fe"] * 20,
        positions=[[i / 20, 0.0, 0.0] for i in range(20)],
        spins=[1.0] * 20,
        g_factors=[2.0] * 20,
        magnetic=[True] * 20,
    )

    assert len(repr(crystal)) < 100


def test_repr_reports_atom_counts():
    N = 9
    M = 4
    crystal = make(
        names=["Fe"] * N,
        positions=[[i / N, 0.0, 0.0] for i in range(N)],
        spins=[1.0] * N,
        g_factors=[2.0] * N,
        magnetic=[True] * M + [False] * (N - M),
    )
    representation = repr(crystal)
    assert "9 atoms" in representation
    assert "4 magnetic" in representation


def test_repr_reports_lattice_parameters():
    # lengths only
    crystal = make(cell=np.diag([3.0, 4.0, 12.0]))
    representation = repr(crystal)
    assert "3.0" in representation
    assert "4.0" in representation
    assert "1.20" in representation


def test_repr_lattice_parameters_are_vector_lengths():
    # and not a diagonal entries
    crystal = make(cell=[[3.0, 4.0, 0.0], [0.0, 7.0, 0.0], [0.0, 0.0, 8.0]])
    assert "5.00" in repr(crystal)


def test_repr_distinguishes_different_crystals():
    assert repr(make()) != repr(make(cell=2 * np.eye(3)))
    assert repr(make()) != repr(make_with_ligand())


def test_repr_names_the_class():
    assert repr(make()).startswith("_Crystal")


def test_magnetic_atoms_is_a_crystal():
    assert isinstance(make_with_ligand().magnetic_atoms, _Crystal)


def tets_magnetic_atoms_keeps_only_magnetic_ones():
    subset = make_with_ligand().magnetic_atoms
    assert subset.names == ("Fe", "Fe")
    assert len(subset) == subset.M == 2
    assert np.all(subset.magnetic)


def test_magnetic_atoms_carries_the_right_values():
    crystal = make_with_ligand()
    subset = crystal.magnetic_atoms
    np.testing.assert_equal(subset.spins, crystal.spins[crystal.map_to_all])
    np.testing.assert_equal(subset.positions, crystal.positions[crystal.map_to_all])


def test_magnetic_atoms_of_all_magnetic_crystal_is_itself():
    crystal = make()
    assert crystal == crystal.magnetic_atoms


def test_magnetic_atoms_is_idempotent():
    crystal = make_with_ligand()
    assert crystal.magnetic_atoms.magnetic_atoms == crystal.magnetic_atoms


def test_magnetic_atoms_are_all_magnetic():
    subset = make_with_ligand().magnetic_atoms
    assert subset.magnetic.all()


def test_magnetic_atoms_is_cached():
    crystal = make_with_ligand()
    assert crystal.magnetic_atoms is crystal.magnetic_atoms


def test_magnetic_atoms_agree_with_select():
    crystal = make_with_ligand()
    assert crystal.magnetic_atoms == crystal.select(crystal.map_to_all)
