cdef inline long prec_to_dps(n):
    return max(1, int(round(int(n)/3.3219280948873626)-1))

cdef inline long dps_to_prec(n):
    return max(1, int(round((int(n)+1)*3.3219280948873626)))