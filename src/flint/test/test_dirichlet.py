from __future__ import annotations

from flint import dirichlet_char, dirichlet_group, fmpz
from flint.test.helpers import is_close_acb, raises


def test_dirichlet_group_basics() -> None:
    g = dirichlet_group(5)
    assert g.size() == 4
    assert g.q == 5
    assert g.exponent() == 4
    assert repr(g) == "Dirichlet group mod q = 5"
    assert str(g) == "Dirichlet group mod q = 5"

    assert raises(lambda: dirichlet_group(0), AssertionError)


def test_dirichlet_char_properties() -> None:
    chi = dirichlet_char(7, 3)

    assert isinstance(chi.group(), dirichlet_group)
    assert chi.index() == 1
    assert chi.modulus() == 7
    assert chi.number() == 3
    assert chi.order() == 6
    assert chi.is_real() is False
    assert chi.is_primitive() is True
    assert chi.conductor() == 7
    assert chi.is_principal() is False
    assert chi.parity() == 1

    assert repr(chi) == "dirichlet_char(7, 3)"
    assert str(chi) == "dirichlet_char(7, 3)"

    assert raises(lambda: dirichlet_char(8, 2), AssertionError)


def test_dirichlet_char_values_and_mul() -> None:
    chi = dirichlet_char(7, 3)
    principal = dirichlet_char(7, 1)
    chi5 = dirichlet_char(7, 5)

    assert is_close_acb(chi(7), 0)
    assert is_close_acb(chi(2), -0.5 + 0.866025403784439j, tol=1e-12, rel_tol=1e-12)
    assert is_close_acb(chi(fmpz(2)), -0.5 + 0.866025403784439j, tol=1e-12, rel_tol=1e-12)

    assert chi.chi_exponent(7) is None
    assert chi.chi_exponent(2) == 2
    assert principal.chi_exponent(2) == 0

    prod = chi * chi5
    assert prod.number() == 1
    for n in [1, 2, 3, 4, 5, 6, 7, 8]:
        assert is_close_acb(prod(n), chi(n) * chi5(n), tol=1e-12, rel_tol=1e-12)

    assert raises(lambda: chi * dirichlet_char(5, 1), AssertionError)
    assert raises(lambda: chi(object()), TypeError)  # type: ignore[arg-type]
    assert raises(lambda: chi.chi_exponent(object()), TypeError)  # type: ignore[arg-type]


def test_dirichlet_char_special_functions() -> None:
    chi = dirichlet_char(7, 3)
    assert is_close_acb(chi.l_function(2), 0.902247025301257 + 0.232548981277895j, tol=1e-12, rel_tol=1e-12)
    assert is_close_acb(chi.l(2), chi.l_function(2), tol=1e-12, rel_tol=1e-12)

    zeta_chi = dirichlet_char(1, 1)
    assert is_close_acb(zeta_chi.hardy_z(1), -0.736305462867318, tol=1e-12, rel_tol=1e-12)
