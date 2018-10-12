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

    def __repr__(self):
        return "Dirichlet group mod q = %s" % self.q

    def __str__(self):
        return repr(self)


cdef class dirichlet_char(object):
    """
    Represents a Dirichlet character.

    Calling the character evaluates the character at the
    given integer, returning an acb.

        >>> chi = dirichlet_char(dirichlet_group(7), 0)
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
        >>> chi = dirichlet_char(dirichlet_group(7), 1)
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

    def __init__(self, dirichlet_group G, ulong index=0):
        """
        Creates the Dirichlet character of given index belonging
        to the group G.

            >>> G = dirichlet_group(5)
            >>> chi = dirichlet_char(G, 0)

        """
        self.G = G
        dirichlet_char_init(self.val, self.G.val)
        assert index < dirichlet_group_size(self.G.val)
        dirichlet_char_index(self.val, self.G.val, index)

    @property
    def index(self):
        return dirichlet_index_char(self.G.val, self.val)

    def __repr__(self):
        return "Dirichlet character index %s in group mod q = %s" % (self.index, self.G.q)

    def __str__(self):
        return repr(self)

    def __call__(self, n):
        cdef acb v
        cdef fmpz m
        m = fmpz(n) % self.G.q
        v = acb.__new__(acb)
        acb_dirichlet_chi((<acb>v).val, self.G.val, self.val, fmpz_get_ui(m.val), getprec())
        return v

    def l(self, s):
        """
        Evaluates the Dirichlet L-function of this character at the given
        complex number s.

            >>> chi = dirichlet_char(dirichlet_group(1), 0)
            >>> showgood(lambda: chi.l(2), dps=25)
            1.644934066848226436472415
            >>> chi = dirichlet_char(dirichlet_group(7), 1)
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

            >>> chi = dirichlet_char(dirichlet_group(1), 0)
            >>> showgood(lambda: chi.hardy_z(1), dps=25)
            -0.7363054628673177346778998

        """
        s = any_as_acb(s)
        cdef acb v
        v = acb.__new__(acb)
        acb_dirichlet_hardy_z((<acb>v).val, (<acb>s).val, self.G.val, self.val, 1, getprec())
        return v
