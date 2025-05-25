from flint.flintlib.functions.fmpz cimport fmpz_abs, fmpz_root, fmpz_set
from flint.types.fmpz cimport fmpz, any_as_fmpz
from flint.utils.typecheck cimport typecheck

cdef class qfb:
    """
    The qfb type represents definite binary quadratic forms
    over Z, with composition, inverse and power operations
    compatible with the class group of a given discriminant.

    Some operations require the form to be primitive.
    """
    def __cinit__(self):
        qfb_init(self.val)
        self.D = fmpz(0)

    def __dealloc__(self):
        qfb_clear(self.val)

    def __init__(self, a, b, c):
        a_fmpz = any_as_fmpz(a)
        b_fmpz = any_as_fmpz(b)
        c_fmpz = any_as_fmpz(c)
        if a_fmpz is NotImplemented:
            raise TypeError(f"Incorrect type {type(a)} for qfb coefficient")
        if b_fmpz is NotImplemented:
            raise TypeError(f"Incorrect type {type(b)} for qfb coefficient")
        if c_fmpz is NotImplemented:
            raise TypeError(f"Incorrect type {type(c)} for qfb coefficient")
        fmpz_set(self.val[0].a, (<fmpz>a_fmpz).val)
        fmpz_set(self.val[0].b, (<fmpz>b_fmpz).val)
        fmpz_set(self.val[0].c, (<fmpz>c_fmpz).val)
        D = fmpz()
        qfb_discriminant(D.val, self.val)
        self.D = D

    def __repr__(self):
        a, b, c = self.coefficients()
        return f"qfb({a}, {b}, {c})"

    def __eq__(self, other):
        if self is other:
            return True

        if typecheck(other, qfb):
            return bool(qfb_equal(self.val, (<qfb>other).val))

        return False

    def __mul__(q1, q2):
        "Returns a reduced form equivalent to the composition of q1 and q2"
        if not q1.is_primitive():
            raise ValueError(f"{q1} is not primitive")

        cdef qfb res = qfb.__new__(qfb)
        cdef fmpz_t L
        fmpz_abs(L, q1.D.val)
        fmpz_root(L, L, 4)
        qfb_nucomp(res.val, q1.val, (<qfb>q2).val, q1.D.val, L)
        qfb_reduce(res.val, res.val, q1.D.val)
        res.D = q1.D
        return res

    def __pow__(q, e, mod):
        "Returns a reduced form equivalent to the e-th iterated composition of q"
        if mod is not None:
            raise NotImplementedError("modular exponentiation")

        if not q.is_primitive():
            raise ValueError(f"{q} is not primitive")

        e_fmpz = any_as_fmpz(e)
        if e_fmpz is NotImplemented:
            raise TypeError(f"exponent cannot be cast to an fmpz type: {e}")

        # qfb_pow does not support negative exponents and will loop forever
        # if a negative integer is provided.
        e_abs = abs(e_fmpz)

        cdef qfb res = qfb.__new__(qfb)
        qfb_pow(res.val, q.val, q.D.val, (<fmpz>e_abs).val)
        if e_fmpz < 0:
            qfb_inverse(res.val, res.val)
        res.D = q.D
        return res

    def coefficients(self):
        """
        Returns coefficients (a, b, c) of the form as a polynomial q(x,y)=ax²+bxy+cy²
        """
        a = fmpz()
        fmpz_set(a.val, self.val[0].a)
        b = fmpz()
        fmpz_set(b.val, self.val[0].b)
        c = fmpz()
        fmpz_set(c.val, self.val[0].c)
        return a, b, c

    def discriminant(self):
        return self.D

    def is_reduced(self):
        return bool(qfb_is_reduced(self.val))

    def is_primitive(self):
        return bool(qfb_is_primitive(self.val))

    def inverse(self):
        cdef qfb res = qfb.__new__(qfb)
        qfb_inverse(res.val, self.val)
        res.D = self.D
        return res

    def reduce(self):
        cdef qfb res = qfb.__new__(qfb)
        qfb_reduce(res.val, self.val, self.D.val)
        res.D = self.D
        return res

    @classmethod
    def prime_form(cls, D, p):
        """
        Returns the unique reduced form with 0 < b ≤ p. Requires that p is prime.
        """

        d_fmpz = any_as_fmpz(D)
        p_fmpz = any_as_fmpz(p)

        cdef qfb res = qfb.__new__(qfb)
        qfb_prime_form(res.val, (<fmpz>d_fmpz).val, (<fmpz>p_fmpz).val)
        res.D = d_fmpz
        return res
