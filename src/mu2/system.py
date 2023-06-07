import numpy as np
from scipy import optimize

from . import scatter
from .bind import spectrum
from .interaction import Interaction

from .cc import cscatter

class System:
    '''
    The system consists of two identical particles.
    '''
    def __init__(self,
                 interaction: Interaction,
                 mu: float,
                 ell: int,
    ):
        self.interaction = interaction
        self.r_c = interaction.R
        self.mu = mu
        self.ell = ell
        self.q = self.interaction.counterterm.qmesh.nodes
        self.wq = self.interaction.counterterm.qmesh.weights
        self.qmax = self.interaction.counterterm.qmesh.upper_bound
        self.v_tilde = self.interaction.matrix_elements(self.interaction.counterterm.qmesh, self.ell, self.ell)
    
    
    def phase_shifts(self, glo, gnlo, ks):
        '''
        Calculates phase shifts (in radians).
        '''
        xterm = self.interaction.counterterm.gen(glo, gnlo)
        v = self.v_tilde + xterm
        return np.array(
            [scatter.phase_shift(ki, v, self.q, self.wq, self.qmax,
                2*self.mu, degrees=False) for ki in ks]
        )

        
    def cross_sections(self, glo, gnlo, ks):
        '''
        Calculates the P-wave term in the cross section.
        '''
        deltas = self.phase_shifts(glo, gnlo, ks)
        return np.array(
            [4*np.pi/k**2 * 3 * np.sin(delta)**2 for (k, delta) in zip(ks, deltas)]
        )


    def effective_range_expansion(self, glo, gnlo, ks, use_c=False):
        '''
        Calculates k^(2l+1) cot(delta).
        '''
        xterm = self.interaction.counterterm.gen(glo, gnlo)
        v = self.v_tilde + xterm
        if use_c:
            return self.kcotd_gen_fast(ks, glo, gnlo)
        else:
            return np.array(
                [ki**(2*self.ell)*scatter.kcotdelta(ki, v, self.q, self.wq,
                    self.qmax, 2*self.mu) for ki in ks]
            )
    
    
    def kcotd(self, p, glo, gnlo):    
        '''
        Calculates k^(2l+1) cot(delta) as at a specified momentum, p.
        '''
        xterm = self.interaction.counterterm.gen(glo, gnlo)
        v = self.v_tilde + xterm
        scatter_temp = scatter.kcotdelta(p, v, self.q, self.wq, self.qmax, 2*self.mu)
        return p**(2*self.ell)*scatter_temp
    
    
    def kcotd_gen(self, ks, glo, gnlo):    
        '''
        Calculates k^(2l+1) cot(delta) as function of ks (array) so that we can
        curve_fit the phaseshifts.
        '''
        xterm = self.interaction.counterterm.gen(glo, gnlo)
        v = self.v_tilde + xterm
        return np.array(
            [ki**(2*self.ell)*scatter.kcotdelta(ki, v, self.q, self.wq,
                self.qmax, 2*self.mu) for ki in ks]
        )
    

    def kcotd_gen_pert1(self, ks, glo, gnlo):
        v0 = self.v_tilde + self.interaction.counterterm.gen(glo, 0)
        v1 = self.interaction.counterterm.gen(0, gnlo)
        return np.array([scatter.kcotdelta_pert1(k, v0, v1, self.q, self.wq, self.qmax, 2*self.mu) for k in ks])
    
    
    def kcotd_gen_pert1a(self, ks, glo, gnlo):
        '''
        Long and Yang implementation of NLO, perturbative scattering
        '''
        v0 = self.v_tilde + self.interaction.counterterm.gen(glo, 0)
        v1 = self.interaction.counterterm.gen(0, gnlo)
        return np.array([scatter.kcotdelta_pert1a(k, v0, v1, self.q, self.wq, self.qmax, 2*self.mu) for k in ks])
    
    
    def T_gen_fast(self, ks, glo, gnlo):    
        '''
        Calculates T (see Thompson&Nunes Table  3.1) as function of ks (array) so that we can
        curve_fit the phaseshifts. Uses the C code inside the cc dir. Needs to
        be compiled (try `make libs` inside cc) to a shared library first.
        '''
        xterm = self.interaction.counterterm.gen(glo, gnlo)
        v = self.v_tilde + xterm
        return np.array(
            [-self.mu*ki*np.pi*cscatter.t_on_shell_py(ki, v, self.q, self.wq, self.qmax, 2*self.mu)
                for ki in ks],
            dtype=np.complex128
        )


    def kcotd_gen_fast(self, ks, glo, gnlo):    
        '''
        Calculates k^(2l+1) cot(delta) as function of ks (array) so that we can
        curve_fit the phaseshifts. Uses the C code inside the cc dir. Needs to
        be compiled (try `make libs` inside cc) to a shared library first.
        '''
        xterm = self.interaction.counterterm.gen(glo, gnlo)
        v = self.v_tilde + xterm
        return np.array(
            [cscatter.kcotdelta_py(ki, v, self.q, self.wq, self.qmax, self.ell, 2*self.mu) for ki in ks]
        )
    

    def kcotd_gen_pert1_fast(self, ks, glo, gnlo):    
        '''
        Calculates k^(2l+1) cot(delta) as function of ks (array) so that we can
        curve_fit the phaseshifts. Uses the C code inside the cc dir. Needs to
        be compiled (try `make libs` inside cc) to a shared library first.
        '''
        v0 = self.v_tilde + self.interaction.counterterm.gen(glo, 0)
        v1 = self.interaction.counterterm.gen(0, gnlo)
        return np.array(
            [cscatter.kcotdelta_pert1_py(ki, v0, v1, self.q, self.wq, self.qmax, self.ell, 2*self.mu) for ki in ks]
        )
    
    
    def t_gen_pert1_fast(self, ks, glo, gnlo):
        '''
        Computes the real and imaginary components of t = t^(0) + t^(1).
        '''
        v0 = self.v_tilde + self.interaction.counterterm.gen(glo, 0)
        v1 = self.interaction.counterterm.gen(0, gnlo)
        return np.array(
            [cscatter.t_pert1_sum_py(ki, v0, v1, self.q, self.wq, self.qmax, 2*self.mu) for ki in ks]
        )
    
    
    def a0_and_r0_pert1_fast(self, glo, gnlo, ks, p0=None, return_sigma=False):
        '''
        Returns a_0 and r_0 after fitting the effective range parameters.
        '''
        assert self.ell == 0, 'This is not an S-wave (l = 0) system.'
        kcds = self.kcotd_gen_pert1_fast(ks, glo, gnlo)
        result = optimize.curve_fit(
            lambda x, c0, c2, c3: c0 + c2*x**2 + c3*x**3,
            ks, kcds, p0=p0, maxfev=20000
        )
        pars, cov = result
        sig = np.sqrt(np.diag(cov))
        if return_sigma:
            return -1/pars[0], 2*pars[1], sig
        else:
            return -1/pars[0], 2*pars[1]
    
    
    def effective_range_parameters(self, glo, gnlo, ks, expansion, p0=None, use_c=False):
        '''
        Returns the effective range parameters for partial wave, ell, at given
        LO (glo) and NLO (gnlo) strengths. The expected analytical form of the
        effective range function must be provided (see expansion).
        '''
        kcds = self.effective_range_expansion(glo, gnlo, ks, use_c=use_c)
        result = optimize.curve_fit(expansion, ks, kcds, p0=p0, maxfev=20000)
        pars = result[0]
        cov = result[1]
        sig = np.sqrt(np.diag(cov))
        return pars, sig
    
    
    def a0_and_r0(self, glo, gnlo, ks, p0=None, use_c=False):
        '''
        Returns a_0 and r_0 after fitting the effective range parameters.
        '''
        assert self.ell == 0, 'This is not an S-wave (l = 0) system.'
        pars, _ = self.effective_range_parameters(glo, gnlo, ks,
            lambda x, c0, c2, c3: c0 + c2*x**2 + c3*x**3,
            p0=p0, use_c=use_c
        )
        return -1/pars[0], 2*pars[1]
    
    
    def a1_and_r1(self, glo, gnlo, ks, p0=None, use_c=False):
        '''
        Returns a_1 and r_1 after fitting the effective range parameters.
        '''
        assert self.ell == 1, 'This is not a P-wave (l = 1) system.'
        pars, _ = self.effective_range_parameters(glo, gnlo, ks,
            lambda x, c0, c1, c2, c3: c0 + c1*x + c2*x**2 + c3*x**3,
            p0=p0, use_c=use_c
        )
        return -1/pars[0], 2*pars[2]
    
    
    def bound_state_spectrum(self, glo, gnlo):
        xterm = self.interaction.counterterm.gen(glo, gnlo)
        v = self.v_tilde + xterm
        return spectrum(v, 0, 0, self.q, self.wq, self.mu)
