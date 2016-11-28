cdef class dirichlet_group(object):

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

    @property
    def q(self):
        return fmpz(self.val.q)

    def __repr__(self):
        return "Dirichlet group mod q = %s" % self.q

    def __str__(self):
        return repr(self)


cdef class dirichlet_char(object):

    cdef dirichlet_char_t val
    cdef dirichlet_group G

    def __cinit__(self):
        pass

    def __dealloc__(self):
        dirichlet_char_clear(self.val)

    def __init__(self, dirichlet_group G, ulong index=0):
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
        s = any_as_acb(s)
        cdef acb v
        v = acb.__new__(acb)
        acb_dirichlet_l((<acb>v).val, (<acb>s).val, self.G.val, self.val, getprec())
        return v

    def hardy_z(self, s):
        s = any_as_acb(s)
        cdef acb v
        v = acb.__new__(acb)
        acb_dirichlet_hardy_z((<acb>v).val, (<acb>s).val, self.G.val, self.val, 1, getprec())
        return v


