from cpython.version cimport PY_MAJOR_VERSION

cdef inline long prec_to_dps(n):
    return max(1, int(round(int(n)/3.3219280948873626)-1))

cdef inline long dps_to_prec(n):
    return max(1, int(round((int(n)+1)*3.3219280948873626)))

cdef inline chars_from_str(s):
    if PY_MAJOR_VERSION < 3:
        return s
    else:
        return s.encode('ascii')

cdef inline str_from_chars(s):
    if PY_MAJOR_VERSION < 3:
        return str(s)
    else:
        return bytes(s).decode('ascii')

cdef inline matrix_to_str(tab):
    if len(tab) == 0 or len(tab[0]) == 0:
        return "[]"
    tab = [[str(c) for c in row] for row in tab]
    widths = []
    for i in xrange(len(tab[0])):
        w = max([len(row[i]) for row in tab])
        widths.append(w)
    for i in xrange(len(tab)):
        tab[i] = [s.rjust(widths[j]) for j, s in enumerate(tab[i])]
        tab[i] = "[" + (", ".join(tab[i])) + "]"
    return "\n".join(tab)

cdef inline _str_trunc(s, trunc=0):
    if trunc > 0 and len(s) > 3 * trunc:
        left = right = trunc
        omitted = len(s) - left - right
        return s[:left] + ("{...%s digits...}" % omitted) + s[-right:]
    return s

