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
        >>> acb_mat([[1j,0],[0,2*1j]]).theta(acb_mat([[0],[0]])).transpose() # doctest: +SKIP
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
        >>> print(acb_mat([[1j, 0],[0,1j]]).theta(acb_mat([[0],[0]])).transpose().str(25)) # doctest: +SKIP
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


def acb_theta_jets(acb_mat z, acb_mat tau, slong ord):
    r"""
    Computes the coefficients of the Taylor expansion of the vector valued Riemann
    theta function `(\theta_{a,b}(z, \tau) : a, b \in \{0,1\}^{g})` or its squares.

    This is a wrapper for the C-function
    `acb_theta_jet <https://flintlib.org/doc/acb_theta.html#c.acb_theta_jet>`_
    and it follows the same conventions for the ordering of the theta characteristics.

    This should be used via the method :meth:`.acb_mat.theta_jets`, explicitly ``tau.theta_jets(z, ord)``.

        >>> from flint import acb, acb_mat, showgood, ctx
        >>> z = acb(1+1j); tau = acb(1.25+3j)
        >>> acb_mat([[tau]]).theta_jets(acb_mat([[z]]),2)  # doctest: +SKIP
        [[0.969443038779670 +/- 5.67e-16] + [-0.0305569612081680 +/- 5.13e-17]j,
        [-0.191993710594950 +/- 4.89e-16] + [0.191993710747776 +/- 7.42e-16]j,
        [0.60317023860834 +/- 2.93e-15] + [0.60317023764810 +/- 4.86e-15]j]
        [[1.03055696119601 +/- 3.89e-15] + [0.0305569612081680 +/- 5.13e-17]j,
        [0.191993710594950 +/- 4.89e-16] + [-0.191993710442123 +/- 4.62e-16]j,
        [-0.60317023668787 +/- 5.90e-15] + [-0.60317023764810 +/- 4.86e-15]j]
        [[-1.22079026757697 +/- 4.36e-15] + [-1.82705551679115 +/- 5.17e-15]j,
        [-5.71849316258739 +/- 7.02e-15] + [3.82088827346268 +/- 5.75e-15]j,
        [6.0241074288587 +/- 3.88e-14] + [9.0163253443780 +/- 2.05e-14]j]
        [[-1.82023591012499 +/- 2.67e-15] + [1.21625195015448 +/- 4.14e-15]j,
        [3.8353056542516 +/- 4.99e-14] + [5.73981078971270 +/- 6.74e-15]j,
        [8.9823364151977 +/- 2.44e-14] + [-6.0022138700195 +/- 3.72e-14]j]
    """
    g = tau.nrows()
    if g == 0:
        return acb_mat(0, 0)

    # Calculate the length of the jet for one characteristic
    # This is the number of multi-indices (alpha) such that |alpha| < ord
    cdef slong nj = acb_theta_jet_nb(g, ord)

    # Total number of characteristics
    cdef slong nb = 1 << (2 * g)

    # Total number of acb elements to allocate
    # FLINT stores nj coefficients for each of the nb characteristics
    cdef slong total_size = nb * nj

    cdef acb_ptr zvec = _acb_vec_init(g)
    cdef slong i, j
    for i in range(g):
        acb_set(zvec + i, acb_mat_entry(z.val, i, 0))

    # Parameters for characteristics
    cdef slong nb_in = 1  # Number of input z vectors
    cdef ulong ab = 0     # Base characteristic
    cdef ulong all = True  # Compute all 2^2g characteristics
    cdef ulong square = False  # Don't compute the squares of the thetas.

    # Initialize the output buffer
    cdef acb_ptr theta = _acb_vec_init(total_size)

    # Call the FLINT C function
    # Note: Computes all partial derivatives up to total order 'ord'
    acb_theta_jet(theta, zvec, nb_in, tau.val, ord, ab, all, square, getprec())

    # Copy the output into a structured format
    res_mat = acb_mat(nb, nj)
    for i in range(nb):
        for j in range(nj):
            acb_set(acb_mat_entry(res_mat.val, i, j), theta + (i * nj + j))

    # Cleanup
    _acb_vec_clear(zvec, g)
    _acb_vec_clear(theta, total_size)

    return res_mat

