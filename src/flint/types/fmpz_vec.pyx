from flint.flintlib.fmpz cimport fmpz_struct, fmpz_set, fmpz_init_set
from flint.flintlib.flint cimport slong
from flint.flintlib.fmpz_vec cimport _fmpz_vec_init, _fmpz_vec_clear

from flint.types.fmpz cimport fmpz, any_as_fmpz

cimport libc.stdlib

cdef class fmpz_vec:
    def __cinit__(self, slong length, bint double_indirect=False):
        self.val = _fmpz_vec_init(length)
        self.length = length
        if double_indirect:
            self.double_indirect = <fmpz_struct **> libc.stdlib.malloc(length * sizeof(fmpz_struct *))
            if self.double_indirect is NULL:
                raise MemoryError("malloc returned a null pointer")

            for i in range(length):
                self.double_indirect[i] = &self.val[i]
        else:
            self.double_indirect = NULL

    def __dealloc__(self):
        libc.stdlib.free(self.double_indirect)
        _fmpz_vec_clear(self.val, self.length)

    def __getitem__(self, x):
        if not isinstance(x, int):
            raise TypeError("index is not integer")
        elif not 0 <= x < self.length:
            raise IndexError("index out of range")

        cdef fmpz z = fmpz.__new__(fmpz)
        fmpz_init_set(z.val, &self.val[x])
        return z

    def __setitem__(self, x, y):
        if not isinstance(x, int):
            raise TypeError("index is not integer")
        elif not 0 <= x < self.length:
            raise IndexError("index out of range")

        y = any_as_fmpz(y)
        if y is NotImplemented:
            raise TypeError("argument is not coercible to fmpz")

        fmpz_set(&self.val[x], (<fmpz>y).val)

    def __len__(self):
        return self.length

    def str(self, *args):
        s = [None] * self.length
        for i in range(self.length):
            x = <fmpz>fmpz.__new__(fmpz)
            fmpz_init_set(x.val, &self.val[i])
            s[i] = x.str(*args)
        return str(s)

    def repr(self, *args):
        return f"fmpz_vec({self.str(*args)}, {self.length})"

    @staticmethod
    def from_iterable(iterable, double_indirect: bool = False):
        length = len(iterable)

        vec = fmpz_vec(length, double_indirect)
        for i, x in enumerate(iterable):
            x = any_as_fmpz(x)
            if x is NotImplemented:
                raise TypeError("argument is not coercible to fmpz")

            fmpz_set(&vec.val[i], (<fmpz>x).val)
        return vec

    def to_tuple(self):
        t = tuple(fmpz.__new__(fmpz) for _ in range(self.length))
        for i in range(self.length):
            fmpz_init_set((<fmpz>t[i]).val, &self.val[i])
        return t
