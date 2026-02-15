"""Tests for python-flint's `fmpz_vec` type."""

from __future__ import annotations

from flint import fmpz, fmpz_vec
from flint.test.helpers import raises


def test_fmpz_vec_construct_int_and_iter() -> None:
    v = fmpz_vec(3)
    assert len(v) == 3
    assert str(v) == "['0', '0', '0']"
    assert repr(v) == "fmpz_vec(['0', '0', '0'], 3)"

    v[0] = 1
    v[1] = fmpz(2)
    v[2] = -3
    assert int(v[0]) == 1
    assert int(v[1]) == 2
    assert int(v[2]) == -3
    assert [int(x) for x in v] == [1, 2, -3]
    assert v.str() == "['1', '2', '-3']"

    w = fmpz_vec([4, 5, 6])
    assert len(w) == 3
    assert [int(x) for x in w] == [4, 5, 6]
    assert str(w) == "['4', '5', '6']"
    assert repr(w) == "fmpz_vec(['4', '5', '6'], 3)"


def test_fmpz_vec_double_indirect_and_errors() -> None:
    v = fmpz_vec(2, double_indirect=True)
    v[0] = 10
    v[1] = 11
    assert [int(x) for x in v] == [10, 11]

    assert raises(lambda: v[2], IndexError)
    assert raises(lambda: v[-1], IndexError)

    def get_bad_index() -> None:
        _ = v[1.25]  # type: ignore[index]

    def set_bad_index() -> None:
        v[1.25] = 7  # type: ignore[index]

    def set_bad_value() -> None:
        v[0] = object()  # type: ignore[assignment]

    def set_oob() -> None:
        v[2] = 1

    assert raises(get_bad_index, TypeError)
    assert raises(set_bad_index, TypeError)
    assert raises(set_bad_value, TypeError)
    assert raises(set_oob, IndexError)


def test_fmpz_vec_constructor_invalid_inputs() -> None:
    assert raises(lambda: fmpz_vec(-1), ValueError)
    assert raises(lambda: fmpz_vec(-1, double_indirect=True), ValueError)
    assert raises(lambda: fmpz_vec(1 << 62), OverflowError)
    assert raises(lambda: fmpz_vec(object()), TypeError)  # type: ignore[arg-type]
