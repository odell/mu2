import os
import sys
from ctypes import cdll, c_int, c_double, POINTER

import numpy as np

LIB = np.ctypeslib.load_library('libkcd.so', os.environ['HOME'])

kcotdelta = LIB.kcotdelta
kcotdelta.argtypes = (
    c_double,
    POINTER(c_double),
    POINTER(c_double),
    POINTER(c_double),
    c_int,
    c_double,
    c_int,
    c_double
)
kcotdelta.restype = c_double

def kcotdelta_py(k, V, q, wq, qmax, l, mass):
    Vp = V.copy()
    qp = q.copy()
    wqp = wq.copy()
    return kcotdelta(k,
            Vp.ctypes.data_as(POINTER(c_double)), 
            qp.ctypes.data_as(POINTER(c_double)),
            wqp.ctypes.data_as(POINTER(c_double)),
            qp.size, qmax, l, mass)