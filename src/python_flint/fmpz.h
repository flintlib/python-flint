#ifndef PYTHON_FLINT_FMPZ_H
#define PYTHON_FLINT_FMPZ_H

#include <Python.h>
#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#define PYFLINT_FMPZ_ABI_VERSION 1
#define PYFLINT_FMPZ_CAPSULE_NAME "flint.types.fmpz._C_API"

typedef struct PyFlint_FMPZ_API_v1 {
    uint32_t abi_version;
    size_t struct_size;
    PyObject *(*fmpz_add)(PyObject *, PyObject *);
    const void *(*fmpz_get_value)(PyObject *);
} PyFlint_FMPZ_API_v1;

#ifndef PYFLINT_FMPZ_MODULE
static PyFlint_FMPZ_API_v1 *PyFlint_FMPZ_API = NULL;

static int
PyFlint_Import(void)
{
    PyFlint_FMPZ_API_v1 *api;

    api = (PyFlint_FMPZ_API_v1 *)PyCapsule_Import(
        PYFLINT_FMPZ_CAPSULE_NAME, 0);
    if (api == NULL)
        return -1;
    if (api->abi_version != PYFLINT_FMPZ_ABI_VERSION) {
        PyErr_Format(PyExc_ImportError,
                     "unsupported python-flint C API version %u",
                     (unsigned int)api->abi_version);
        return -1;
    }
    if (api->struct_size < sizeof(PyFlint_FMPZ_API_v1)) {
        PyErr_SetString(PyExc_ImportError,
                        "python-flint C API table is too small");
        return -1;
    }
    PyFlint_FMPZ_API = api;
    return 0;
}

static PyObject *
PyFlint_FMPZ_Add(PyObject *a, PyObject *b)
{
    if (PyFlint_FMPZ_API == NULL && PyFlint_Import() < 0)
        return NULL;
    return PyFlint_FMPZ_API->fmpz_add(a, b);
}

static const void *
PyFlint_FMPZ_GetValue(PyObject *value)
{
    if (PyFlint_FMPZ_API == NULL && PyFlint_Import() < 0)
        return NULL;
    return PyFlint_FMPZ_API->fmpz_get_value(value);
}
#endif

#ifdef __cplusplus
}
#endif
#endif
