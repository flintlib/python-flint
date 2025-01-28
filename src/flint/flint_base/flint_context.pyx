from flint.flintlib.types.arf cimport arf_rnd_t
from flint.flintlib.types.flint cimport (
    flint_cleanup,
    flint_get_num_threads,
    flint_set_num_threads
)
from flint.utils.conversion cimport prec_to_dps, dps_to_prec

from functools import wraps

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
        """
        Adds n bits of precision to the current flint context.

            >>> from flint import arb, ctx
            >>> with ctx.extraprec(5): x = arb(2).sqrt().str()
            >>> x
            '[1.414213562373095 +/- 5.53e-17]'

        This function also works as a wrapper:

            >>> from flint import arb, ctx
            >>> @ctx.extraprec(10)
            ... def f(x):
            ...     return x.sqrt().str()
            >>> f(arb(2))
            '[1.41421356237309505 +/- 1.46e-18]'
        """
        return self.workprec(n + self.prec)

    def extradps(self, n):
        """
        Adds n digits of precision to the current flint context.

            >>> from flint import arb, ctx
            >>> with ctx.extradps(5): x = arb(2).sqrt().str()
            >>> x
            '[1.4142135623730950488 +/- 2.76e-21]'

        This function also works as a wrapper:

            >>> from flint import arb, ctx
            >>> @ctx.extradps(10)
            ... def f(x):
            ...     return x.sqrt().str()
            >>> f(arb(2))
            '[1.414213562373095048801689 +/- 3.13e-25]'
        """
        return self.workdps(n + self.dps)

    def workprec(self, n):
        """
        Sets the working precision for the current flint context,
        using a python context manager.

            >>> from flint import arb, ctx
            >>> with ctx.workprec(5): x = arb(2).sqrt().str()
            >>> x
            '[1e+0 +/- 0.438]'

        This function also works as a wrapper:

            >>> from flint import arb, ctx
            >>> @ctx.workprec(24)
            ... def f(x):
            ...     return x.sqrt().str()
            >>> f(arb(2))
            '[1.41421 +/- 3.66e-6]'
        """
        return PrecisionManager(self, eprec=n)

    def workdps(self, n):
        """
        Sets the working precision in digits for the current
        flint context, using a python context manager.

            >>> from flint import arb, ctx
            >>> with ctx.workdps(5): x = arb(2).sqrt().str()
            >>> x
            '[1.4142 +/- 1.51e-5]'

        This function also works as a wrapper:

            >>> from flint import arb, ctx
            >>> @ctx.workdps(10)
            ... def f(x):
            ...     return x.sqrt().str()
            >>> f(arb(2))
            '[1.414213562 +/- 3.85e-10]'
        """
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


cdef class PrecisionManager:
    cdef FlintContext ctx
    cdef int eprec
    cdef int edps
    cdef int _oldprec

    def __init__(self, ctx, eprec=-1, edps=-1):
        if eprec != -1 and edps != -1:
            raise ValueError("two different precisions requested")

        self.ctx = ctx

        self.eprec = eprec
        self.edps = edps

    def __call__(self, func):
        @wraps(func)
        def wrapped(*args, **kwargs):
            _oldprec = self.ctx.prec

            try:
                if self.eprec != -1:
                    self.ctx.prec = self.eprec

                if self.edps != -1:
                    self.ctx.dps = self.edps

                return func(*args, **kwargs)
            finally:
                self.ctx.prec = _oldprec

        return wrapped

    def __enter__(self):
        self._oldprec = self.ctx.prec

        if self.eprec != -1:
            self.ctx.prec = self.eprec

        if self.edps != -1:
            self.ctx.dps = self.edps

    def __exit__(self, type, value, traceback):
        self.ctx.prec = self._oldprec


cdef FlintContext thectx = FlintContext()
