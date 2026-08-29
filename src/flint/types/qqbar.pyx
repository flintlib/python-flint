from flint.types.fmpz cimport fmpz
from flint.types.fmpq cimport fmpq

cdef class qqbar(flint_scalar):
    """
    The qqbar class represents algebraic numbers.
    """

    def __cinit__(self):
        qqbar_init(self.val)

    def __dealloc__(self):
        qqbar_clear(self.val)

    def __init__(self, val=None):
        if val is not None:
            if isinstance(val, qqbar):
                qqbar_set(self.val, (<qqbar>val).val)
            elif isinstance(val, int):
                qqbar_set_si(self.val, val)
            elif isinstance(val, fmpz):
                qqbar_set_fmpz(self.val, (<fmpz>val).val)
            elif isinstance(val, fmpq):
                qqbar_set_fmpq(self.val, (<fmpq>val).val)
            else:
                raise TypeError("Cannot initialize qqbar type from input")

    def str(self, Py_ssize_t digits=15):
        """
        Converts the algebraic number to a decimal string representation.
        """
        import sys
        import io

        # Temporarily intercept stdout using a clean, cross-platform string buffer
        old_stdout = sys.stdout
        sys.stdout = io.StringIO()

        try:
            qqbar_printn(self.val, digits)
            sys.stdout.flush()
            output = sys.stdout.getvalue()
        finally:
            sys.stdout = old_stdout

        return output.strip()

    def __str__(self):
        return self.str()

    def __repr__(self):
        return f"qqbar({self.str()})"

    def __add__(x, y):
        cdef qqbar res = qqbar()
        cdef qqbar c_x = x if isinstance(x, qqbar) else qqbar(x)
        cdef qqbar c_y = y if isinstance(y, qqbar) else qqbar(y)
        qqbar_add(res.val, c_x.val, c_y.val)
        return res

    def __sub__(x, y):
        cdef qqbar res = qqbar()
        cdef qqbar c_x = x if isinstance(x, qqbar) else qqbar(x)
        cdef qqbar c_y = y if isinstance(y, qqbar) else qqbar(y)
        qqbar_sub(res.val, c_x.val, c_y.val)
        return res

    def __mul__(x, y):
        cdef qqbar res = qqbar()
        cdef qqbar c_x = x if isinstance(x, qqbar) else qqbar(x)
        cdef qqbar c_y = y if isinstance(y, qqbar) else qqbar(y)
        qqbar_mul(res.val, c_x.val, c_y.val)
        return res

    def __truediv__(x, y):
        cdef qqbar res = qqbar()
        cdef qqbar c_x = x if isinstance(x, qqbar) else qqbar(x)
        cdef qqbar c_y = y if isinstance(y, qqbar) else qqbar(y)
        qqbar_div(res.val, c_x.val, c_y.val)
        return res
