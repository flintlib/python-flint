from flint.flintlib.arf cimport ARF_RND_DOWN
from flint.flintlib.flint cimport (
    flint_cleanup,
    flint_get_num_threads,
    flint_set_num_threads
)
from flint.utils.conversion cimport prec_to_dps, dps_to_prec

cdef class FlintContext:
    def __init__(self):
        self.default()

    def default(self):
        self.pretty = True
        self.rnd = ARF_RND_DOWN
        self.prec = 53
        self.unicode = False
        self.threads = 1
        self.cap = 10

    @property 
    def prec(self):
        return self._prec

    @prec.setter
    def prec(self, prec):
        cdef long cprec = prec
        if cprec < 2:
            raise ValueError("prec must be >= 2")
        self._prec = cprec
        self._dps = prec_to_dps(cprec)

    @property
    def dps(self):
        return self._dps

    @dps.setter
    def dps(self, prec):
        self.prec = dps_to_prec(prec)

    @property
    def cap(self):
        return self._cap

    @cap.setter
    def cap(self, long cap):
        if cap < 0:
            raise ValueError("cap must be >= 0")
        self._cap = cap

    @property
    def threads(self):
        return flint_get_num_threads()
    
    @threads.setter
    def threads(self, long num):
        assert num >= 1 and num <= 64
        flint_set_num_threads(num)

    def __repr__(self):
        return "pretty = %-8s  # pretty-print repr() output\n" \
               "unicode = %-8s # use unicode characters in output\n" \
               "prec = %-8s    # real/complex precision (in bits)\n"   \
               "dps = %-8s     # real/complex precision (in digits)\n"    \
               "cap = %-8s     # power series precision\n"    \
               "threads = %-8s # max number of threads used internally" % \
            (self.pretty, self.unicode, self.prec, self.dps, self.cap, self.threads)

    def cleanup(self):
        flint_cleanup()

cdef FlintContext thectx = FlintContext()
