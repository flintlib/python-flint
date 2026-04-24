from __future__ import annotations

from flint.types.acb_mat import acb_mat


def acb_theta(z: acb_mat, tau: acb_mat, square: bool | int = False) -> acb_mat: ...

def acb_theta_jets(z: acb_mat, ord: int, square: bool = False) -> list[acb_mat]: ...
