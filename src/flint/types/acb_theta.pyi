from __future__ import annotations

from typing import Protocol

from flint.types.acb_mat import acb_mat


class _AcbThetaFunc(Protocol):
    def __call__(self, z: acb_mat, tau: acb_mat, square: bool | int = False) -> acb_mat: ...


class _AcbThetaJetsFunc(Protocol):
    def __call__(self, z: acb_mat, tau: acb_mat, ord: int) -> acb_mat: ...


acb_theta: _AcbThetaFunc | None

acb_theta_jets: _AcbThetaJetsFunc | None
