cdef dict _dirichlet_group_cache = {}

cdef class dirichlet_group(object):
    """
    Represents the group of Dirichlet characters modulo a given q.

        >>> G = dirichlet_group(5)
        >>> G
        Dirichlet group mod q = 5
    """

    cdef dirichlet_group_t val
    cdef int _init

    def __cinit__(self):
        self._init = 0

    def __dealloc__(self):
        if self._init != 0:
            dirichlet_group_clear(self.val)

    def __init__(self, ulong q):
        assert q >= 1
        dirichlet_group_init(self.val, q)
        self._init = 1

    cpdef long size(self):
        return dirichlet_group_size(self.val)

    @property
    def q(self):
        return fmpz(self.val.q)

    def exponent(self):
        return fmpz(self.val.expo)

    def __repr__(self):
        return "Dirichlet group mod q = %s" % self.q

    def __str__(self):
        return repr(self)


cdef class dirichlet_char(object):
    """
    Represents a Dirichlet character.

    Calling the character evaluates the character at the
    given integer, returning an acb.

        >>> chi = dirichlet_char(7, 1)
        >>> for n in range(7):
        ...     print(chi(n))
        ... 
        0
        1.00000000000000
        1.00000000000000
        1.00000000000000
        1.00000000000000
        1.00000000000000
        1.00000000000000
        >>> chi = dirichlet_char(7, 3)
        >>> for n in range(7):
        ...     print(chi(n))
        ... 
        0
        1.00000000000000
        -0.500000000000000 + [0.866025403784439 +/- 5.15e-16]j
        0.500000000000000 + [0.866025403784439 +/- 5.15e-16]j
        -0.500000000000000 + [-0.866025403784439 +/- 5.15e-16]j
        0.500000000000000 + [-0.866025403784439 +/- 5.15e-16]j
        -1.00000000000000
    """

    cdef dirichlet_char_t val
    cdef dirichlet_group G

    def __cinit__(self):
        pass

    def __dealloc__(self):
        dirichlet_char_clear(self.val)

    def __init__(self, ulong q, ulong l=1):
        """
        Creates the Dirichlet character with Conrey number (q, l).

            >>> chi = dirichlet_char(5, 1)

        """
        if q in _dirichlet_group_cache:
            G = _dirichlet_group_cache[q]
        else:
            G = dirichlet_group(q)
            _dirichlet_group_cache[q] = G
        self.G = G
        dirichlet_char_init(self.val, self.G.val)
        assert 1 <= l <= max(q,2)-1 and n_gcd(q, l) == 1
        dirichlet_char_log(self.val, self.G.val, l)

    def group(self):
        return self.G

    def index(self):
        return dirichlet_index_char(self.G.val, self.val)

    def modulus(self):
        return self.G.q

    def number(self):
        return dirichlet_char_exp(self.G.val, self.val)

    def order(self):
        return dirichlet_order_char(self.G.val, self.val)

    def is_real(self):
        return bool(dirichlet_char_is_real(self.G.val, self.val))

    def is_primitive(self):
        return bool(dirichlet_char_is_primitive(self.G.val, self.val))

    def conductor(self):
        return dirichlet_conductor_char(self.G.val, self.val)

    def is_principal(self):
        return bool(dirichlet_char_is_principal(self.G.val, self.val))

    def parity(self):
        return dirichlet_parity_char(self.G.val, self.val)

    def __repr__(self):
        return "dirichlet_char(%i, %i)" % (self.modulus(), self.number())

    def __str__(self):
        return repr(self)

    def __mul__(dirichlet_char self, dirichlet_char other):
        cdef dirichlet_char res
        assert self.G.q == other.G.q
        res = dirichlet_char(self.G.q, 1)
        dirichlet_char_mul(res.val, self.G.val, self.val, other.val)
        return res

    def __call__(self, n):
        cdef acb v
        cdef fmpz m
        m = fmpz(n) % self.G.q
        v = acb.__new__(acb)
        acb_dirichlet_chi((<acb>v).val, self.G.val, self.val, fmpz_get_ui(m.val), getprec())
        return v

    def chi_exponent(self, n):
        cdef fmpz m
        cdef ulong v
        m = fmpz(n) % self.G.q
        expo = self.G.exponent()
        v = dirichlet_chi(self.G.val, self.val, fmpz_get_ui(m.val))
        if v == DIRICHLET_CHI_NULL:
            return None
        else:
            return fmpz(v)

    def l(self, s):
        """
        Evaluates the Dirichlet L-function of this character at the given
        complex number s.

            >>> chi = dirichlet_char(1, 1)
            >>> showgood(lambda: chi.l(2), dps=25)
            1.644934066848226436472415
            >>> chi = dirichlet_char(7, 3)
            >>> showgood(lambda: chi.l(2+3j), dps=25)
            1.273313649440490751755284 - 0.07432329442559421607102118j

        """
        s = any_as_acb(s)
        cdef acb v
        v = acb.__new__(acb)
        acb_dirichlet_l((<acb>v).val, (<acb>s).val, self.G.val, self.val, getprec())
        return v

    def hardy_z(self, s):
        """
        Evaluates the Hardy Z-function of this character at the given
        complex number s.

            >>> chi = dirichlet_char(1, 1)
            >>> showgood(lambda: chi.hardy_z(1), dps=25)
            -0.7363054628673177346778998

        """
        s = any_as_acb(s)
        cdef acb v
        v = acb.__new__(acb)
        acb_dirichlet_hardy_z((<acb>v).val, (<acb>s).val, self.G.val, self.val, 1, getprec())
        return v

