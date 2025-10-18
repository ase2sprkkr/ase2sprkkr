import cython
import numpy as np
cimport numpy as np

cdef extern from "aaa.h":
    cdef void aaa_(double complex* cdata,
                   double* data,
                   int* a,
                   int* b
                   );


def call_a(double complex[:,::1] cdata, double[:,::1] data, int a, int b):
    aaa_(&cdata[0,0], &data[0,0], &a, &b)

def call_b(cdata, data, int a, int b):
    cdef double complex[::1,:] _cdata = np.asfortranarray(cdata)
    cdef double[::1,:] _data = np.asfortranarray(data)
    aaa_(&_cdata[0,0], &_data[0,0], &a, &b)

