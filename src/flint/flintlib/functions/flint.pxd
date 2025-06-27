from flint.flintlib.types.flint cimport flint_rand_t, slong

# unknown type ...
# unknown type FILE
# unknown type flint_err_t
# unknown type int ( cmp) ( void
# unknown type size_t
# unknown type va_list
# unknown type va_list)
# unknown type void (

# .. macro:: __FLINT_VERSION
# .. macro:: __FLINT_RELEASE
# .. macro:: FLINT_VERSION
# .. macro:: FLINT_BITS
# .. macro:: FLINT_D_BITS
# .. macro:: FLINT_ABS(x)
# .. macro:: FLINT_UABS(x)
# .. macro:: FLINT_MIN(x, y)
# .. macro:: FLINT_SWAP(T, x, y)
# .. macro:: FLINT_SGN(x)
# .. macro:: UWORD_MIN

cdef extern from "flint/flint.h":
    # void * flint_malloc(size_t size)
    # void * flint_realloc(void * ptr, size_t size)
    # void * flint_calloc(size_t num, size_t size)
    void flint_free(void * ptr)
    void flint_rand_init(flint_rand_t state)
    void flint_rand_clear(flint_rand_t state)
    void flint_set_num_threads(int num_threads)
    # int flint_get_num_threads(void)
    int flint_set_num_workers(int num_workers)
    void flint_reset_num_workers(int num_workers)
    # int flint_printf(const char * format, ...)
    # int flint_fprintf(FILE * fs, const char * format, ...)
    # int flint_vprintf(const char * format, va_list vlist)
    # int flint_vfprintf(FILE * fs, const char * format, va_list vlist)
    # int flint_sprintf(char * s, const char * str, ...)
    # int flint_scanf(const char * str, ...)
    # int flint_fscanf(FILE * f, const char * str, ...)
    # int flint_sscanf(const char * s, const char * str, ...)
    # void flint_abort(void)
    # void flint_throw(flint_err_t exc, const char * msg, ...)
    # void flint_set_abort(void (* func)(void))
    # void flint_set_throw(void (* func)(flint_err_t, const char *, va_list))
    # void flint_merge_sort(void * buf, slong len, slong size, int (* cmp) (const void *, const void *, void *), void * data)
    # void flint_sort(void * buf, slong len, slong size, int (* cmp) (const void *, const void *, void *), void * data)
