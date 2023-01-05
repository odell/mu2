'''
Functions for computing scattering observables.
'''

import numpy as np
from scipy.interpolate import interp2d

def generate_new_v_matrix(v, q, q0):
    nq = len(q)
    vp = np.zeros((nq+1, nq+1))
    f = interp2d(q, q, v, kind='cubic')
    for i in range(nq):
        for j in range(nq):
            vp[i, j] = v[i, j]
        vp[-1, i] = f(q[i], q0)[0]
        vp[i, -1] = f(q0, q[i])[0]
    vp[nq, nq] = f(q0, q0)[0]
    return vp


def kernel_matrix_gen(q0, vp, q_nodes, q_weights, qmax, mass):
    nq = np.size(q_nodes)
    kernel = np.zeros((nq+1, nq+1), dtype=complex)
    logdiff_numerical = np.dot(q_weights, list(map(lambda x: x**2/(q0**2-x**2),
        q_nodes)))
    logdiff_analytical = -qmax-q0/2*np.log(1-2*q0/(q0+qmax)) - 1j*np.pi*q0/2
    logdiff = logdiff_analytical - logdiff_numerical
    for i in range(nq+1):
        for j in range(nq):
            kernel[i, j] = float(i == j) - mass*q_weights[j]*q_nodes[j]**2*vp[i, j] / (q0**2-q_nodes[j]**2)
        kernel[i, nq] = float(i == nq) - mass*vp[i, nq]*logdiff
    return kernel


def t_on_shell(q0, v_matrix, q_nodes, q_weights, qmax, mass):
    nq = np.size(q_nodes)
    vp = generate_new_v_matrix(v_matrix, q_nodes, q0)
    kernel = kernel_matrix_gen(q0, vp, q_nodes, q_weights, qmax, mass)
    t = np.linalg.solve(kernel, vp[:, nq])
    return t[-1]


def kcotdelta(q0, v_matrix, q_nodes, q_weights, qmax, mass):
    t_onsh = t_on_shell(q0, v_matrix, q_nodes, q_weights, qmax, mass)
    tr = np.real(t_onsh)
    ti = np.imag(t_onsh)
    return -2*tr/(mass*np.pi*(tr**2+ti**2))


def phase_shift(q0, v_matrix, q_nodes, q_weights, qmax, mass, degrees=True):
    kcd = kcotdelta(q0, v_matrix, q_nodes, q_weights, qmax, mass)
    delta = np.arctan(q0/kcd)
    if degrees:
        delta *= 180/np.pi
    return delta


# First-order perturbation functions

def driving_term_pert1(q0, v1, t0_half_shell, q, wq, qmax, mass):
    nq = np.size(q)
    kernel = np.zeros((nq+1, nq+1), dtype=complex)
    logdiff_numerical = np.dot(wq, list(map(lambda x: x**2/(q0**2-x**2), q)))
    logdiff_analytical = -qmax - q0/2*np.log(1-2*q0/(q0+qmax)) - 1j*np.pi*q0/2
    logdiff = logdiff_analytical - logdiff_numerical
    for i in range(nq+1):
        for j in range(nq):
            kernel[i, j] = mass*wq[j]*q[j]**2*v1[i, j] / (q0**2-q[j]**2)
        kernel[i, nq] = mass*v1[i, nq]*logdiff
    return v1[:, -1] + kernel @ t0_half_shell
    

def t_on_shell_pert1(q0, v0, v1, q_nodes, q_weights, qmax, mass):
    v0p = generate_new_v_matrix(v0, q_nodes, q0)
    K = kernel_matrix_gen(q0, v0p, q_nodes, q_weights, qmax, mass)
    t0_half_shell = np.linalg.solve(K, v0p[:, -1])

    v1p = generate_new_v_matrix(v1, q_nodes, q0)
    d1 = driving_term_pert1(q0, v1p, t0_half_shell, q_nodes, q_weights, qmax, mass)

    return t0_half_shell[-1]  + np.linalg.solve(K, d1)[-1]


def kcotdelta_pert1(q0, v0, v1, q_nodes, q_weights, qmax, mass):
    t_onsh = t_on_shell_pert1(q0, v0, v1, q_nodes, q_weights, qmax, mass)
    tr = np.real(t_onsh)
    ti = np.imag(t_onsh)
    return -2*tr/(mass*np.pi*(tr**2+ti**2))
