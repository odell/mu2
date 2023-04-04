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

kcotdelta_pert1 = LIB.kcotdelta_pert1
kcotdelta_pert1.argtypes = (
    c_double,
    POINTER(c_double),
    POINTER(c_double),
    POINTER(c_double),
    POINTER(c_double),
    c_int,
    c_double,
    c_int,
    c_double
)
kcotdelta_pert1.restype = c_double

def kcotdelta_pert1_py(k, v0, v1, q, wq, qmax, l, mass):
    v0p = v0.copy()
    v1p = v1.copy()
    qp = q.copy()
    wqp = wq.copy()
    return kcotdelta_pert1(k,
            v0p.ctypes.data_as(POINTER(c_double)), 
            v1p.ctypes.data_as(POINTER(c_double)), 
            qp.ctypes.data_as(POINTER(c_double)),
            wqp.ctypes.data_as(POINTER(c_double)),
            qp.size, qmax, l, mass)

tps = LIB.t_pert1_sum
tps.argtypes = (
    POINTER(c_double),
    POINTER(c_double),
    c_double,
    POINTER(c_double),
    POINTER(c_double),
    POINTER(c_double),
    POINTER(c_double),
    c_int,
    c_double,
    c_double
)
def t_pert1_sum_py(k, v0, v1, q, wq, qmax, mass):
    tr = np.zeros(1)
    ti = np.zeros(1)
    tps(tr.ctypes.data_as(POINTER(c_double)),
        ti.ctypes.data_as(POINTER(c_double)),
        k,
        v0.ctypes.data_as(POINTER(c_double)),
        v1.ctypes.data_as(POINTER(c_double)),
        q.ctypes.data_as(POINTER(c_double)),
        wq.ctypes.data_as(POINTER(c_double)),
        q.size, qmax, mass)
    return tr[0], ti[0]