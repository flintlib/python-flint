from flint.flintlib.types.arf cimport arf_rnd_t
from flint.flintlib.types.flint cimport (
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
        self.rnd = arf_rnd_t.ARF_RND_DOWN
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

    def extraprec(self, n):
        return self.workprec(n + self.prec)

    def extradps(self, n):
        return self.workdps(n + self.dps)

    def workprec(self, n):
        return PrecisionManager(self, eprec=n)

    def workdps(self, n):
        return PrecisionManager(self, edps=n)

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


class PrecisionManager:
    def __init__(self, ctx, eprec=None, edps=None):
        if eprec is not None and edps is not None:
            raise ValueError("two different precisions requested")

        self.ctx = ctx

        self.eprec = eprec
        self.edps = edps

    def __enter__(self):
        self._oldprec = self.ctx.prec

        if self.eprec is not None:
            self.ctx.prec = self.eprec

        if self.edps is not None:
            self.ctx.dps = self.edps

    def __exit__(self, type, value, traceback):
        self.ctx.prec = self._oldprec


cdef FlintContext thectx = FlintContext()
