from flint._flint cimport (
    ARF_RND_DOWN,
    arf_rnd_t,
    flint_cleanup,
    flint_get_num_threads,
    flint_set_num_threads
)
from flint.utils.conversion cimport prec_to_dps, dps_to_prec

cdef class FlintContext:
    cdef public bint pretty
    cdef public long _prec
    cdef public long _dps
    cdef arf_rnd_t rnd
    cdef public bint unicode
    cdef public long _cap

    def __init__(self):
        self.default()

    def default(self):
        self.pretty = True
        self.rnd = ARF_RND_DOWN
        self.prec = 53
        self.unicode = False
        self.threads = 1
        self.cap = 10

    property prec:

        def __set__(self, prec):
            cdef long cprec = prec
            if cprec < 2:
                raise ValueError("prec must be >= 2")
            self._prec = cprec
            self._dps = prec_to_dps(cprec)

        def __get__(self):
            return self._prec

    property dps:

        def __set__(self, prec):
            self.prec = dps_to_prec(prec)

        def __get__(self):
            return self._dps

    property cap:

        def __set__(self, long cap):
            if cap < 0:
                raise ValueError("cap must be >= 0")
            self._cap = cap

        def __get__(self):
            return self._cap

    property threads:

        def __set__(self, long num):
            assert num >= 1 and num <= 64
            flint_set_num_threads(num)

        def __get__(self):
            return flint_get_num_threads()

    def __repr__(self):
        return "pretty = %-8s  # pretty-print repr() output\n" \
               "unicode = %-8s # use unicode characters in output\n" \
               "prec = %-8s    # real/complex precision (in bits)\n"   \
               "dps = %-8s     # real/complex precision (in digits)\n"    \
               "cap = %-8s     # power series precision\n"    \
               "threads = %-8s # max number of threads used internally\n" % \
            (self.pretty, self.unicode, self.prec, self.dps, self.cap, self.threads)

    def cleanup(self):
        flint_cleanup()
