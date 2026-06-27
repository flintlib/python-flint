from __future__ import annotations

from flint import acb, acb_mat
from flint.test.helpers import is_close_acb_mat as is_close, raises
from flint.types.acb_theta import acb_theta, acb_theta_jets


def test_acb_theta_basic() -> None:
    if acb_theta is None:
        return
    theta = acb_theta

    z = acb(1 + 1j)
    tau = acb(1.25 + 3j)
    zmat = acb_mat([[z]])
    taumat = acb_mat([[tau]])

    direct = theta(zmat, taumat)
    via_method = taumat.theta(zmat)
    assert is_close(direct, via_method, tol=1e-12, rel_tol=1e-12, max_width=1e-12)
    assert direct.nrows() == 1
    assert direct.ncols() == 4

    squared = theta(zmat, taumat, square=True)
    assert squared.nrows() == 1
    assert squared.ncols() == 4


def test_acb_theta_shape_assertions() -> None:
    if acb_theta is None:
        return

    z = acb_mat([[0]])
    tau_bad = acb_mat([[1j, 0]])
    assert raises(lambda: acb_theta(z, tau_bad), AssertionError)

    tau = acb_mat([[1j]])
    z_bad_rows = acb_mat([[0], [0]])
    assert raises(lambda: acb_theta(z_bad_rows, tau), AssertionError)

    z_bad_cols = acb_mat([[0, 0]])
    assert raises(lambda: acb_theta(z_bad_cols, tau), AssertionError)

    assert raises(lambda: acb_theta(object(), tau), TypeError)  # type: ignore[arg-type]
    assert raises(lambda: acb_theta(z, object()), TypeError)  # type: ignore[arg-type]


def test_acb_theta_jets_basic() -> None:
    if acb_theta_jets is None:
        return

    z = acb(1 + 1j)
    tau = acb(1.25 + 3j)
    zmat = acb_mat([[z]])
    taumat = acb_mat([[tau]])
    ord = 2

    direct = acb_theta_jets(zmat, taumat, ord)
    via_method = taumat.theta_jets(zmat, ord)
    assert is_close(direct, via_method, tol=1e-12, rel_tol=1e-12, max_width=1e-12)
    assert direct.nrows() == 4
    assert direct.ncols() == 3
