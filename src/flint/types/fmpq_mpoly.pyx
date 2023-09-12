

cdef dict _fmpq_mpoly_ctx_cache = {}

@cython.auto_pickle(False)
cdef class fmpq_mpoly_ctx(flint_mpoly_context):
    """
    A class for storing the polynomial context

    :param nvars: The number of variables in the ring
    :param ordering:  The term order for the ring
    :param names:  A tuple containing the names of the variables of the ring.

    Do not construct one of these directly, use `get_fmpz_mpoly_context`.
    """
#    cdef fmpz_mpoly_ctx_t val

    def __init__(self, slong nvars, ordering, names):
        if ordering == "lex":
            fmpq_mpoly_ctx_init(self.val, nvars, ordering_t.ORD_LEX)
        elif ordering == "deglex":
            fmpq_mpoly_ctx_init(self.val, nvars, ordering_t.ORD_DEGLEX)
        elif ordering == "degrevlex":
            fmpq_mpoly_ctx_init(self.val, nvars, ordering_t.ORD_DEGREVLEX)
        else:
            raise ValueError("Unimplemented term order %s" % ordering)

        super().__init__(nvars, names)

    cpdef slong nvars(self):
        """
        Return the number of variables in the context

            >>> ctx = get_fmpz_mpoly_context(4, "lex", 'x')
            >>> ctx.nvars()
            4
        """
        return self.val.zctx.minfo.nvars

    cpdef ordering(self):
        """
        Return the term order of the context object.

            >>> ctx = get_fmpz_mpoly_context(4, "deglex", 'w')
            >>> ctx.ordering()
            'deglex'
        """
        if self.val.zctx.minfo.ord == ordering_t.ORD_LEX:
            return "lex"
        if self.val.zctx.minfo.ord == ordering_t.ORD_DEGLEX:
            return "deglex"
        if self.val.zctx.minfo.ord == ordering_t.ORD_DEGREVLEX:
            return "degrevlex"

    def gen(self, slong i):
        """
        Return the `i`th generator of the polynomial ring

            >>> ctx = get_fmpz_mpoly_context(3, 'degrevlex', 'z')
            >>> ctx.gen(1)
            z1
        """
        cdef fmpq_mpoly res
        assert i >= 0 and i < self.val.zctx.minfo.nvars
        res = fmpq_mpoly.__new__(fmpz_mpoly)
        res.ctx = self
        fmpq_mpoly_init(res.val, res.ctx.val)
        res._init = True
        fmpq_mpoly_gen(res.val, i, res.ctx.val)
        return res

    def constant(self, z):
        """
        Create a constant polynomial in this context
        """
        cdef fmpq_mpoly res
        z = any_as_fmpq(z)
        if z is NotImplemented:
            raise ValueError("A constant fmpq_mpoly is a fmpq")
        res = fmpz_mpoly.__new__(fmpz_mpoly)
        res.ctx = self
        fmpz_mpoly_init(res.val, res.ctx.val)
        res._init = True
        fmpz_mpoly_set_fmpz(res.val, (<fmpz>z).val, res.ctx.val)
        return res

    def fmpz_mpoly_from_dict(self, d):
         """
         Create a fmpz_mpoly from a dictionary.

         The dictionary's keys are tuples of ints (or anything that implicitly converts
         to fmpz) representing exponents, and corresponding values of fmpz.

             >>> ctx = get_fmpz_mpoly_context(2,'lex','x,y')
             >>> ctx.fmpz_mpoly_from_dict({(1,0):2, (1,1):3, (0,1):1})
             3*x*y + 2*x + y
         """
         cdef long n
         cdef fmpz_t coefficient
         cdef fmpz_struct *exponents
         cdef int xtype
         cdef int nvars = self.nvars()
         cdef int i,j
         cdef int count
         cdef fmpz_mpoly res

         if not PyDict_Check(d):
             raise ValueError("expected a dictionary")
         n = PyDict_Size(d)
         fmpz_init(coefficient)
         exponents = <fmpz_struct *> libc.stdlib.calloc(nvars, sizeof(fmpz_struct))
         if exponents == NULL:
             raise MemoryError()
         for i in range(nvars):
             fmpz_init(exponents + i)
         fmpz_init(coefficient)
         res = fmpz_mpoly.__new__(fmpz_mpoly)
         res.ctx = self
         fmpz_mpoly_init(res.val, res.ctx.val)
         res._init = True
         count = 0
         for k,v in d.items():
             xtype = fmpz_set_any_ref(coefficient, v)
             if xtype == FMPZ_UNKNOWN:
                 libc.stdlib.free(exponents)
                 raise TypeError("invalid coefficient type %s" % type(v))
             if not PyTuple_Check(k):
                 libc.stdlib.free(exponents)
                 raise TypeError("Expected tuple of ints as key not %s" % type(k))
             if PyTuple_GET_SIZE(k) != nvars:
                 libc.stdlib.free(exponents)
                 raise TypeError("Expected exponent tuple of length %d" % nvars)
             for i,tup in enumerate(k):
                 xtype = fmpz_set_any_ref(exponents + i, tup)
                 if xtype == FMPZ_UNKNOWN:
                     libc.stdlib.free(exponents)
                     raise TypeError("Invalid exponent type %s" % type(tup))
             #Todo lobby for fmpz_mpoly_push_term_fmpz_ffmpz
             if not fmpz_is_zero(coefficient):
                 _fmpz_mpoly_push_exp_ffmpz(res.val, exponents, self.val)
                 fmpz_mpoly_set_term_coeff_fmpz(res.val, count, coefficient, self.val)
                 count += 1
         for i in range(nvars):
             fmpz_clear(exponents + i)
         fmpz_clear(coefficient)
         fmpz_mpoly_sort_terms(res.val, self.val)
         return res


def get_fmpz_mpoly_context(slong nvars=1, ordering="lex", names='x'):
    if nvars <= 0:
        nvars = 1
    nametup = tuple(name.strip() for name in names.split(','))
    if len(nametup) != nvars:
        if len(nametup) != 1:
            raise ValueError("Number of variables does not equal number of names")
        nametup = tuple(nametup[0] + str(i) for i in range(nvars))
    key = (nvars, ordering, nametup)
    ctx = _fmpz_mpoly_ctx_cache.get(key)
    if ctx is None:
        ctx = fmpz_mpoly_ctx(nvars, ordering, nametup)
        _fmpz_mpoly_ctx_cache[key] = ctx
    return ctx
