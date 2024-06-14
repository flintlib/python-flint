from flint.flint_base.flint_context cimport getprec
from flint.types.acb cimport acb
from flint.types.acb_mat cimport acb_mat
from flint.flintlib.acb cimport *
from flint.flintlib.acb_mat cimport *
from flint.flintlib.acb_theta cimport *

def acb_mat_theta(acb_mat z, acb_mat tau, ulong square=False):
    r"""
    Computes the vector valued Riemann theta function `(\theta_{a,b}{z, tau) : a, b \in \{0,1\}^{g}\)` or its squares.
    This is a wrapper for the function `acb_theta_all` and it follows the same conventions for the ordering of the theta characteristics.

    This should be used via method `acb_mat.theta` with the order of `z` and `tau` swapped,

        >>> from flint import acb, acb_mat, showgood
        >>> z = acb(1+1j); tau = acb(1.25+3j)
        >>> t0, t1, t2, t3 = acb_mat([[tau]]).theta(acb_mat([[z]]))
        >>> sum([abs(x) for x in acb_mat([z.modular_theta(tau)]) - acb_mat([[-t3,t2,t0,t1]])])
        [+/- 3.82e-14]
        >>> for i in range(4):showgood(lambda: acb_mat([[tau]]).theta(acb_mat([[z]]))[i], dps=25)
        ...
        0.9694430387796704100046143 - 0.03055696120816803328582847j
        1.030556961196006476576271 + 0.03055696120816803328582847j
        -1.220790267576967690128359 - 1.827055516791154669091679j
        -1.820235910124989594900076 + 1.216251950154477951760042j
        >>> acb_mat([[1j,0],[0,2*1j]]).theta(acb_mat([[0],[0]]))
        ([1.09049252082308 +/- 3.59e-15] + [+/- 2.43e-16]j, [1.08237710165638 +/- 4.15e-15] + [+/- 2.43e-16]j, [0.916991251621117 +/- 6.30e-16] + [+/- 2.43e-16]j, [0.910167024735558 +/- 7.93e-16] + [+/- 2.43e-16]j, [0.451696791791346 +/- 5.46e-16] + [+/- 2.43e-16]j, [+/- 2.43e-16] + [+/- 2.43e-16]j, [0.379830212998946 +/- 4.47e-16] + [+/- 2.43e-16]j, [+/- 2.43e-16] + [+/- 2.43e-16]j, [0.916991251621117 +/- 6.30e-16] + [+/- 2.43e-16]j, [0.910167024735558 +/- 7.93e-16] + [+/- 2.43e-16]j, [+/- 2.43e-16] + [+/- 2.43e-16]j, [+/- 2.43e-16] + [+/- 2.43e-16]j, [0.379830212998946 +/- 4.47e-16] + [+/- 2.43e-16]j, [+/- 2.43e-16] + [+/- 2.43e-16]j, [+/- 2.43e-16] + [+/- 2.43e-16]j, [+/- 2.43e-16] + [+/- 2.43e-16]j)

    """
    g = tau.nrows()
    assert tau.ncols() == g
    assert z.nrows() == g
    assert z.ncols() == 1

    # convert input
    cdef acb_ptr zvec
    zvec = _acb_vec_init(g)
    cdef long i
    for i in range(g):
        acb_set(zvec + i, acb_mat_entry(z.val, i, 0))

    # initialize the output
    cdef slong nb = 1 << (2 * g)
    cdef acb_ptr theta = _acb_vec_init(nb)

    acb_theta_all(theta, zvec, tau.val, square, getprec())
    _acb_vec_clear(zvec, g)
    # copy the output
    res = tuple()
    cdef acb r
    for i in range(nb):
        r = acb.__new__(acb)
        acb_set(r.val, theta + i)
        res += (r,)
    _acb_vec_clear(theta, nb)
    return res
