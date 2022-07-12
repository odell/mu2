import numpy as np
from .utility import ft_matrix_gen, Mesh


def local_lo_counterterm(r, R):
    return np.exp(-(r/R)**4)


def local_nlo_counterterm(r, R):
    derivative_factor = (16*r**6 - 12*r**2*R**4) / R**8
    return derivative_factor * local_lo_counterterm(r, R)


def nonlocal_regulator(q, L, n, l):
    return (q/L)**l*np.exp(-(q/L)**n)


def nonlocal_lo_term(p, k, L, n, l):
    return nonlocal_regulator(p, L, n, l) * nonlocal_regulator(k, L, n, l)


def nonlocal_nlo_term(p, k, L, n, l):
    return ((p/L)**2 + (k/L)**2)/2 * nonlocal_lo_term(p, k, L, n, l)


class Counterterm:
    def __init__(self, qmesh, R, ell, nonloc=True):
        self.qmesh = qmesh
        self.R = R
        self.Lambda = 2/self.R
        self.ell = ell
        self.nonloc = nonloc

    def gen(self, glo, gnlo):
        '''
        To be implemented in subclasses.
        '''
        raise ValueError

        
class LocalCounterterm(Counterterm):
    def __init__(self, lo_term, nlo_term, rmesh: Mesh, qmesh: Mesh, R, ell):
        super().__init__(qmesh, R, ell, nonloc=False)
        self.rmesh = rmesh
        self.x_tilde_lo = ft_matrix_gen(lambda r: lo_term(r, R),
                ell, ell, self.qmesh.nodes, self.rmesh.nodes, self.rmesh.weights)
        self.x_tilde_nlo = ft_matrix_gen(lambda r: nlo_term(r, R),
                ell, ell, self.qmesh.nodes, self.rmesh.nodes, self.rmesh.weights)
    

    def gen(self, glo, gnlo):
        return glo*self.x_tilde_lo + gnlo*self.x_tilde_nlo


class NonlocalCounterterm(Counterterm):
    def __init__(self, lo_term, nlo_term, nonloc_reg, qmesh, R, ell):
        super().__init__(qmesh, R, ell, nonloc=True)
        self.lo_term = lo_term # V(p, k, Lambda)
        self.nlo_term = nlo_term # V(p, k, Lambda)
        self.nonloc_reg = nonloc_reg # V(p, k, Lambda)

        self.x_reg = np.array(
            [[self.nonloc_reg(p, k, self.Lambda) for p in self.qmesh.nodes] for
                k in self.qmesh.nodes]
        )

        self.x_tilde_lo = np.array(
            [[self.lo_term(p, k, self.Lambda) for p in self.qmesh.nodes] for k in
                self.qmesh.nodes]
        )
        self.x_tilde_nlo = np.array(
            [[self.nlo_term(p, k, self.Lambda) for p in
                self.qmesh.nodes] for k in self.qmesh.nodes]
        )

    
    def gen(self, glo, gnlo):
        return glo * self.x_tilde_lo + gnlo * self.x_tilde_nlo
