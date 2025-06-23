from flint.flint_base.flint_context cimport getprec
from flint.types.acb cimport acb
from flint.types.acb_mat cimport acb_mat
from flint.flintlib.functions.acb cimport *
from flint.flintlib.types.acb cimport (
    acb_mat_entry,
)
from flint.flintlib.functions.acb_mat cimport *
from flint.flintlib.functions.acb_theta cimport *


def acb_theta(acb_mat z, acb_mat tau, ulong square=False):
    r"""
    Computes the vector valued Riemann theta function
    `(\theta_{a,b}(z, \tau) : a, b \in \{0,1\}^{g})` or its squares.

    This is a wrapper for the C-function
    `acb_theta_all <https://flintlib.org/doc/acb_theta.html#c.acb_theta_all>`_
    and it follows the same conventions for the ordering of the theta characteristics.

    This should be used via the method :meth:`.acb_mat.theta`, explicitly ``tau.theta(z)``.

        >>> from flint import acb, acb_mat, showgood, ctx
        >>> z = acb(1+1j); tau = acb(1.25+3j)
        >>> t0, t1, t2, t3 = acb_mat([[tau]]).theta(acb_mat([[z]])).entries()
        >>> all([abs(x) < 1e-13 for x in acb_mat([z.modular_theta(tau)]) - acb_mat([[-t3,t2,t0,t1]])])
        True
        >>> showgood(lambda: acb_mat([[tau]]).theta(acb_mat([[z]])).transpose(), dps=25)
        [0.9694430387796704100046143 - 0.03055696120816803328582847j]
        [ 1.030556961196006476576271 + 0.03055696120816803328582847j]
        [  -1.220790267576967690128359 - 1.827055516791154669091679j]
        [  -1.820235910124989594900076 + 1.216251950154477951760042j]
        >>> acb_mat([[1j,0],[0,2*1j]]).theta(acb_mat([[0],[0]])).transpose()
        [ [1.09049252082308 +/- 5.07e-15] + [+/- 1.73e-15]j]
        [ [1.08237710165638 +/- 5.64e-15] + [+/- 1.74e-15]j]
        [ [0.91699125162112 +/- 4.17e-15] + [+/- 1.33e-15]j]
        [ [0.91016702473556 +/- 3.01e-15] + [+/- 1.34e-15]j]
        [[0.451696791791346 +/- 4.89e-16] + [+/- 1.86e-16]j]
        [                  [+/- 5.42e-16] + [+/- 5.42e-16]j]
        [[0.379830212998946 +/- 3.47e-16] + [+/- 1.43e-16]j]
        [                  [+/- 4.62e-16] + [+/- 4.62e-16]j]
        [ [0.91699125162112 +/- 6.83e-15] + [+/- 4.00e-15]j]
        [ [0.91016702473556 +/- 5.70e-15] + [+/- 4.02e-15]j]
        [                  [+/- 2.44e-16] + [+/- 2.44e-16]j]
        [                  [+/- 2.45e-16] + [+/- 2.45e-16]j]
        [[0.379830212998946 +/- 6.39e-16] + [+/- 4.35e-16]j]
        [                  [+/- 5.86e-16] + [+/- 5.86e-16]j]
        [                  [+/- 2.60e-17] + [+/- 2.60e-17]j]
        [                  [+/- 1.16e-16] + [+/- 1.16e-16]j]
        >>> ctx.prec = 10000
        >>> print(acb_mat([[1j, 0],[0,1j]]).theta(acb_mat([[0],[0]])).transpose().str(25))
        [ [1.180340599016096226045338 +/- 5.95e-26] + [+/- 2.35e-3010]j]
        [[0.9925441784910574194770081 +/- 3.15e-26] + [+/- 2.79e-3010]j]
        [[0.9925441784910574194770081 +/- 3.15e-26] + [+/- 1.88e-3010]j]
        [[0.8346268416740731862814297 +/- 3.29e-26] + [+/- 2.23e-3010]j]
        [[0.9925441784910574194770081 +/- 3.15e-26] + [+/- 5.80e-3011]j]
        [                          [+/- 1.35e-3009] + [+/- 1.35e-3009]j]
        [[0.8346268416740731862814297 +/- 3.29e-26] + [+/- 4.64e-3011]j]
        [                          [+/- 1.14e-3009] + [+/- 1.14e-3009]j]
        [[0.9925441784910574194770081 +/- 3.15e-26] + [+/- 3.16e-3009]j]
        [[0.8346268416740731862814297 +/- 3.29e-26] + [+/- 3.75e-3009]j]
        [                          [+/- 4.05e-3011] + [+/- 4.05e-3011]j]
        [                          [+/- 4.81e-3011] + [+/- 4.81e-3011]j]
        [[0.8346268416740731862814297 +/- 3.29e-26] + [+/- 7.80e-3010]j]
        [                          [+/- 5.19e-3009] + [+/- 5.19e-3009]j]
        [                          [+/- 1.00e-3011] + [+/- 1.00e-3011]j]
        [                          [+/- 2.81e-3010] + [+/- 2.81e-3010]j]
        >>> ctx.default()

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
    res = []
    cdef acb r
    for i in range(nb):
        r = acb.__new__(acb)
        acb_set(r.val, theta + i)
        res.append(r)
    _acb_vec_clear(theta, nb)
    return acb_mat([res])
