"""Tests for python-flint's `fmpq_vec` type."""

from __future__ import annotations

from flint import fmpq, fmpq_vec
from flint.test.helpers import raises


def test_fmpq_vec_construct_and_iter() -> None:
    v = fmpq_vec(3)
    assert len(v) == 3
    assert str(v) == "['0', '0', '0']"
    assert repr(v) == "fmpq_vec(['0', '0', '0'], 3)"

    v[0] = fmpq(1, 2)
    v[1] = 2
    v[2] = fmpq(-3, 4)
    assert str(v[0]) == "1/2"
    assert str(v[1]) == "2"
    assert str(v[2]) == "-3/4"
    assert [str(x) for x in v] == ["1/2", "2", "-3/4"]
    assert v.str() == "['1/2', '2', '-3/4']"

    w = fmpq_vec([fmpq(5, 6), 7, fmpq(8, 9)])
    assert len(w) == 3
    assert [str(x) for x in w] == ["5/6", "7", "8/9"]
    assert str(w) == "['5/6', '7', '8/9']"
    assert repr(w) == "fmpq_vec(['5/6', '7', '8/9'], 3)"


def test_fmpq_vec_double_indirect_and_errors() -> None:
    v = fmpq_vec(2, double_indirect=True)
    v[0] = fmpq(10, 3)
    v[1] = fmpq(11, 5)
    assert [str(x) for x in v] == ["10/3", "11/5"]

    assert raises(lambda: v[2], IndexError)
    assert raises(lambda: v[-1], IndexError)

    def get_bad_index() -> None:
        _ = v[1.25]  # type: ignore[index]

    def set_bad_index() -> None:
        v[1.25] = fmpq(7, 3)  # type: ignore[index]

    def set_bad_value() -> None:
        v[0] = object()  # type: ignore[assignment]

    def set_oob() -> None:
        v[2] = 1

    assert raises(get_bad_index, TypeError)
    assert raises(set_bad_index, TypeError)
    assert raises(set_bad_value, TypeError)
    assert raises(set_oob, IndexError)
