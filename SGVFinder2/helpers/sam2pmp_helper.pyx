cimport cython
from libc.math cimport pow, log10
from cpython.mem cimport PyMem_Malloc, PyMem_Free
cdef int am1 = ord('a') - 1
cdef double *ps = <double *> PyMem_Malloc(42 * cython.sizeof(double))
cdef double *qs = <double *> PyMem_Malloc(42 * cython.sizeof(double))
if ps is NULL or qs is NULL:
    raise MemoryError()
cdef int x = 1
ps[0] = -100
qs[0] = 0
while x < 42:
    ps[x] = log10(1-pow(10, -0.1 * x))
    qs[x] = -0.1 * x
    x += 1

@cython.boundscheck(False)
cpdef double calc_quality(alp, quals, ref):
    cdef int alpl = len(alp)
    cdef int ql = len(quals)
    cdef int rl = len(ref)
    cdef int *_ref = <int *> PyMem_Malloc(rl * cython.sizeof(int))
    cdef int *_alp0 = <int *> PyMem_Malloc(alpl * cython.sizeof(int))
    cdef int *_alp1 = <int *> PyMem_Malloc(alpl * cython.sizeof(int))
    cdef int *_quals = <int *> PyMem_Malloc(ql * cython.sizeof(int))
    if _alp0 is NULL or _alp1 is NULL or _quals is NULL or _ref is NULL:
        raise MemoryError()
    cdef int i=0
    while i < alpl:
        if alp[i][0] is None:
            _alp0[i] = -100
        else:
            _alp0[i] = alp[i][0]
        if alp[i][1] is None:
            _alp1[i] = -100
        else:
            _alp1[i] = alp[i][1]
        i += 1
    i = 0
    while i < ql:
        _quals[i] = quals[i]
        i += 1
    i = 0
    while i < rl:
        _ref[i] = ord(ref[i])
        i += 1
    return _calc_quality(_alp0, _alp1, _quals, _ref, alpl)

@cython.boundscheck(False)
cdef double _calc_quality(int *alp0, int *alp1, int *quals, int *ref, int alpl):
    cdef int cntr_ref = 0
    cdef int cntr_que = 0
    cdef double res = 0
    cdef int i = 0
    while i < alpl:
        if alp0[i] == -100:
            cntr_ref += 1
        elif alp1[i] == -100: 
            res += qs[quals[cntr_que]]
            cntr_que += 1
        else:
            if ref[cntr_ref] > am1:
                res += qs[quals[cntr_que]]
            else:
                res += ps[quals[cntr_que]]
            cntr_ref += 1
            cntr_que += 1
        i += 1
    PyMem_Free(alp0)
    PyMem_Free(alp1)
    PyMem_Free(quals)
    PyMem_Free(ref)
    return res


