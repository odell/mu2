'''
Defines common constants and functions.
'''

import numpy as np
from scipy.special import roots_legendre, spherical_jn

def kd(i, j):
    '''
    Kronecker delta
    '''
    return 1.0 if i == j else 0.0

def mesh(a, b, n):
    x, w = roots_legendre(n)
    x *= (b - a)/2
    x += (b + a)/2
    w *= (b - a)/2
    return x, w


def log_mesh(a, b, n):
    '''
    Returns a logarithmicaly spaced mesh: nodes (x) and weights (w), from a to
    b (noninclusive).
    '''
    x, w = roots_legendre(n)
    x *= (np.log(b+1)-np.log(a+1))/2
    x += (np.log(b+1)+np.log(a+1))/2
    w *= (np.log(b+1)-np.log(a+1))/2
    w = np.exp(x)*w
    x = np.exp(x)-1
    return x, w


class Mesh:
    def __init__(self, a, b, n, log=True):
        self.lower_bound = a
        self.upper_bound = b
        self.size = n
        if log:
            self.nodes, self.weights = log_mesh(self.lower_bound, self.upper_bound, self.size)
        else:
            self.nodes, self.weights = mesh(self.lower_bound, self.upper_bound,
                    self.n)


def ft_matrix_gen(pot, l, lp, p, r, wr):
    '''
    Given a coordinate-space potential (pot), returns the partial-wave
    projected, momentum-space matrix elements.
    '''
    nq = np.size(p)
    nr = np.size(r)
    A = np.zeros((nq, nr))
    B = np.zeros((nr, nq))
    for (i, pi) in enumerate(p):
        for ((j, rj), wj) in zip(enumerate(r), wr):
            A[i, j] = 2/np.pi * wj * rj**2 * spherical_jn(l, pi*rj)
            B[j, i] = pot(rj) * spherical_jn(lp, pi*rj)
    return np.matmul(A, B)
