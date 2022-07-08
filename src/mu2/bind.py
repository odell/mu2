'''
Functions for computing binding energies.
'''
import numpy as np

def kernel_matrix_gen(en, v_matrix, q, wq, mu):
    nq = np.size(q)
    kernel = np.zeros((nq, nq))
    for i in range(nq):
        qi = q[i]
        G0 = 1/(2*mu*en-qi**2)
        for j in range(nq):
            wj = wq[j]
            qj = q[j]
            mwq2 = 2*mu*wj*qj**2
            kernel[i, j] = float(i == j) - mwq2*v_matrix[i, j]*G0
    return kernel


def det_kernel(en, v_matrix, x_matrix_bare, gi, q, wq, mu):
    v = v_matrix + gi*x_matrix_bare
    k = kernel_matrix_gen(en, v, q, wq, mu)
    return np.linalg.det(k)


def diagonalize(v_bare_matrix, x_bare_matrix, g, p, wp, mu):
    nq = np.size(p)
    v_matrix = v_bare_matrix+g*x_bare_matrix
    ham = np.zeros((nq, nq))
    for (i, pi) in enumerate(p):
        for (j, pj) in enumerate(p):
            ham[i, j] = (i == j)*pi*pj/(2*mu) + wp[j]*pj**2*v_matrix[i, j]
    
    return np.linalg.eig(ham)


def bound_states(v_bare_matrix, x_bare_matrix, g, p, wp, mu):
    evals, evecs = diagonalize(v_bare_matrix, x_bare_matrix, g, p, wp, mu)
    ii = np.sort(np.intersect1d(
        np.where(np.real(evals) < 0)[0],
        np.where(np.imag(evals) == 0)[0]
    ))
    return evals[ii], evecs[:, ii]


def spectrum(v_bare_matrix, x_bare_matrix, g, p, wp, mu):
    evals, _ = diagonalize(v_bare_matrix, x_bare_matrix, g, p, wp, mu)
    bound_state_energies = np.real(np.array(sorted(list(filter(lambda j: np.imag(j) == 0,
        filter(lambda i: i < 0, evals)))))[::-1])
    return bound_state_energies
