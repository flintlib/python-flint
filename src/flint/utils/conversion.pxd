cdef inline long prec_to_dps(n):
    return max(1, int(round(int(n)/3.3219280948873626)-1))

cdef inline long dps_to_prec(n):
    return max(1, int(round((int(n)+1)*3.3219280948873626)))

cdef inline chars_from_str(s):
    return s.encode('ascii')

cdef inline str_from_chars(s):
    return bytes(s).decode('ascii')

cdef inline _str_trunc(s, trunc=0):
    if trunc > 0 and len(s) > 3 * trunc:
        left = right = trunc
        omitted = len(s) - left - right
        return s[:left] + ("{...%s digits...}" % omitted) + s[-right:]
    return s
