# src/flint/types/qqbar.pyx
from flint.flint_base.flint_base cimport flint_scalar
from flint.types.fmpz cimport fmpz
from flint.types.fmpq cimport fmpq
from flint.types.qqbar cimport (
    qqbar_init, qqbar_clear, qqbar_set, qqbar_set_si,
    qqbar_add, qqbar_sub, qqbar_mul, qqbar_div, qqbar
)

cdef class qqbar(flint_scalar):
    """
    The qqbar class represents algebraic numbers inside python-flint.
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

    # Explicit Python Operator Overloads
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
    def str(self, Py_ssize_t digits=15):
        """
        Converts the algebraic number to a decimal string representation by capturing C-level stdout.
        """
        import os
        import sys

        # Save current stdout file descriptor
        cdef int stdout_fd = 1
        cdef int saved_stdout = os.dup(stdout_fd)
        
        # Create a pipe to catch stream output
        pipe_read, pipe_write = os.pipe()
        
        # Redirect stdout to our pipe
        os.dup2(pipe_write, stdout_fd)
        
        try:
            # Execute raw C-level printing function
            qqbar_printn(self.val, digits)
            sys.stdout.flush()
        finally:
            # Revert stdout back to original descriptor
            os.close(pipe_write)
            os.dup2(saved_stdout, stdout_fd)
            os.close(saved_stdout)
            
        # Read the captured string from the pipe
        cdef bytes captured = os.read(pipe_read, 4096)
        os.close(pipe_read)
        
        return captured.decode('utf-8').strip()

    def __str__(self):
        return self.str()

    def __repr__(self):
        # Displays a clean numerical representation
        return f"qqbar({self.str()})"