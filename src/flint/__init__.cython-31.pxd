from cpython.ref cimport PyObject, Py_DECREF

cdef extern from "python_flint/fmpz.h":
    int PyFlint_Import() except -1
    PyObject *PyFlint_FMPZ_Add(PyObject *, PyObject *) except NULL
    const void *PyFlint_FMPZ_GetValue(PyObject *) except NULL

cdef extern class flint.types.fmpz.fmpz [object PyFlint_FMPZ_Object, check_size ignore]:
    pass

cdef inline fmpz fmpz_add(fmpz a, fmpz b):
    cdef PyObject *result_ptr = PyFlint_FMPZ_Add(<PyObject *>a, <PyObject *>b)
    cdef fmpz result = <fmpz>result_ptr
    # Assigning the pointer to an extension-type variable adds a reference.
    # Consume the new reference returned by the C API before returning it.
    Py_DECREF(<object>result_ptr)
    return result

cdef inline const void *fmpz_get_value(fmpz value) except NULL:
    return PyFlint_FMPZ_GetValue(<PyObject *>value)
