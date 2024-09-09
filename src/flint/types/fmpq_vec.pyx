from flint.flintlib.types.fmpq cimport fmpq_struct
from flint.flintlib.functions.fmpq cimport fmpq_set, fmpq_init
from flint.flintlib.functions.fmpq_vec cimport _fmpq_vec_init, _fmpq_vec_clear

from flint.types.fmpq cimport fmpq, any_as_fmpq

cimport libc.stdlib

cdef class fmpq_vec:
    def __cinit__(self, iterable_or_len, bint double_indirect=False):
        if isinstance(iterable_or_len, int):
            self.length = iterable_or_len
        else:
            self.length = len(iterable_or_len)

        self.val = _fmpq_vec_init(self.length)

        if double_indirect:
            self.double_indirect = <fmpq_struct **> libc.stdlib.malloc(self.length * sizeof(fmpq_struct *))
            if self.double_indirect is NULL:
                raise MemoryError("malloc returned a null pointer")

            for i in range(self.length):
                self.double_indirect[i] = &self.val[i]
        else:
            self.double_indirect = NULL

    def __init__(self, iterable_or_len, double_indirect: bool = False):
        if not isinstance(iterable_or_len, int):
            for i, x in enumerate(iterable_or_len):
                self[i] = x

    def __dealloc__(self):
        libc.stdlib.free(self.double_indirect)
        if self.val is not NULL:
            _fmpq_vec_clear(self.val, self.length)
            self.val = NULL

    def __getitem__(self, x):
        if not isinstance(x, int):
            raise TypeError("index is not integer")
        elif not 0 <= x < self.length:
            raise IndexError("index out of range")

        cdef fmpq z = fmpq.__new__(fmpq)
        fmpq_init(z.val)
        fmpq_set(z.val, &self.val[x])
        return z

    def __setitem__(self, x, y):
        if not isinstance(x, int):
            raise TypeError("index is not integer")
        elif not 0 <= x < self.length:
            raise IndexError("index out of range")

        y = any_as_fmpq(y)
        if y is NotImplemented:
            raise TypeError("argument is not coercible to fmpq")

        fmpq_set(&self.val[x], (<fmpq>y).val)

    def __len__(self):
        return self.length

    def __str__(self):
        return self.str()

    def __repr__(self):
        return self.repr()

    def __iter__(self):
        cdef fmpq z
        for i in range(self.length):
            z = fmpq.__new__(fmpq)
            fmpq_set(z.val, &self.val[i])
            yield z

    def str(self, *args):
        s = [None] * self.length
        for i in range(self.length):
            x = <fmpq>fmpq.__new__(fmpq)
            fmpq_init(x.val)
            fmpq_set(x.val, &self.val[i])
            s[i] = x.str(*args)
        return str(s)

    def repr(self, *args):
        return f"fmpq_vec({self.str(*args)}, {self.length})"
