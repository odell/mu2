// Copyright 2022 Daniel Odell

#include "scatter_pert.h"

/*
First-order perturbation functions
*/

void driving_term_pert1(
    gsl_vector_complex* d1, // NLO driving term
    double k,
    double* v1, // NLO interaction (exclusive)
    gsl_vector_complex* t_on_shell, // LO, half-shell t-matrix elements, t^{(0)}
    double* q, // momentum nodes
    double* wq, // momentum weights
    int nq, // number of momentum nodes
    double qmax,
    double mass
) {
    int nqp = nq+1;
    gsl_matrix_complex* kernel = gsl_matrix_complex_alloc(nqp, nqp);
    double logdiff_numerical = 0.0;
    for (int i = 0; i < nq; i++)
      logdiff_numerical += wq[i] * q[i]*q[i]/(k*k - q[i]*q[i]);
    gsl_complex logdiff_analytical = gsl_complex_rect(-qmax - k/2*log(1-2*k/(k+qmax)), -M_PI*k/2);
    gsl_complex logdiff = gsl_complex_sub(logdiff_analytical, gsl_complex_rect(logdiff_numerical, 0));
    for (int i = 0; i < nqp; i++) {
      for (int j = 0; j < nq; j++) {
        gsl_matrix_complex_set(
          kernel, i, j,
          gsl_complex_rect(mass*wq[j]*q[j]*q[j]*v1[i*nqp+j] / (k*k - q[j]*q[j]), 0)
        );
      }
      gsl_matrix_complex_set(kernel, i, nq, gsl_complex_mul_real(logdiff, mass*v1[i*nqp+nq]));
    }
    // construct a complex vector of the last row of v1
    for (int i = 0; i < nqp; i++)
      gsl_vector_complex_set(d1, i, gsl_complex_rect(v1[i*nqp+nq], 0));
    
    // multiply kernel @ t0_half_shell
    gsl_vector_complex* Kt = gsl_vector_complex_alloc(nq+1);
    gsl_blas_zgemv(CblasNoTrans, gsl_complex_rect(1.0, 0.0), kernel, t_on_shell, gsl_complex_rect(0.0, 0.0), Kt);

    // add kernel @ t0_half_shell to the driving term
    gsl_vector_complex_add(
      d1,
      Kt
    );

    gsl_matrix_complex_free(kernel);
    gsl_vector_complex_free(Kt);
}

gsl_complex t_on_shell_pert1(
    double k,
    double* v0, // NLO interaction (exclusive)
    double* v1, // NLO interaction (exclusive)
    double* q, // momentum nodes
    double* wq, // momentum weights
    int nq, // number of momentum nodes
    double qmax,
    double mass

) {
    // Interpolate v0 at half-shell elements corresponding to k
    int nqp = nq+1;
    double v0p[nqp*nqp]; 
    augment_potential_matrix(v0p, k, v0, q, nq);

    // Generate the kernel (only depends on v0)
    gsl_matrix_complex* kernel = gsl_matrix_complex_alloc(nqp, nqp);
    kernel_matrix_gen(kernel, k, v0p, q, wq, nq, qmax, mass);

    // Solve for t0
    gsl_vector_complex* t0 = gsl_vector_complex_alloc(nqp);
    gsl_permutation* perm = gsl_permutation_alloc(nqp);

    gsl_vector_complex* v0p_col = gsl_vector_complex_alloc(nqp);
    for (int i = 0; i < nqp; i++)
      gsl_vector_complex_set(v0p_col, i,
        gsl_complex_rect(v0p[i*nqp+nq], 0.0)
      );

    int s;
    gsl_linalg_complex_LU_decomp(kernel, perm, &s);
    // gsl_complex k00 = gsl_matrix_complex_get(kernel, 0, 0);
    // printf("%.8e %.8e\n", GSL_REAL(k00), GSL_IMAG(k00));
    gsl_linalg_complex_LU_solve(kernel, perm, v0p_col, t0);
    // printf("%.8e %.8e\n", GSL_REAL(k00), GSL_IMAG(k00));
    // printf("\n");

    double v1p[nqp*nqp];
    augment_potential_matrix(v1p, k, v1, q, nq);
    gsl_vector_complex* d1 = gsl_vector_complex_alloc(nqp);
    driving_term_pert1(d1, k, v1p, t0, q, wq, nq, qmax, mass);
    gsl_vector_complex* t1 = gsl_vector_complex_alloc(nqp);
    gsl_linalg_complex_LU_solve(kernel, perm, d1, t1);

    // printf("%.8e %.8e\n", GSL_REAL(gsl_vector_complex_get(t0, nq)), GSL_IMAG(gsl_vector_complex_get(t0, nq)));
    // printf("%.8e %.8e\n", GSL_REAL(gsl_vector_complex_get(t1, nq)), GSL_IMAG(gsl_vector_complex_get(t1, nq)));

    gsl_vector_complex_add(t0, t1);
    gsl_complex result = gsl_vector_complex_get(t0, nq);

    // printf("%.8e %.8e\n", GSL_REAL(gsl_vector_complex_get(t0, nq)), GSL_IMAG(gsl_vector_complex_get(t0, nq)));
    // printf("%.8e %.8e\n", GSL_REAL(gsl_vector_complex_get(t1, nq)), GSL_IMAG(gsl_vector_complex_get(t1, nq)));
    // printf("\n");

    gsl_vector_complex_free(t0);
    gsl_vector_complex_free(v0p_col);
    gsl_vector_complex_free(t1);
    gsl_vector_complex_free(d1);
    gsl_permutation_free(perm);
    gsl_matrix_complex_free(kernel);

    return result;
}

double kcotdelta_pert1(
    double k,
    double* v0, // LO interaction
    double* v1, // NLO interaction (exclusive)
    double* q, // momentum nodes
    double* wq, // momentum weights
    int nq, // number of momentum nodes
    double qmax,
    int l,
    double mass
) {
  gsl_complex t = t_on_shell_pert1(k, v0, v1, q, wq, nq, qmax, mass);
  double tr = GSL_REAL(t), ti = GSL_IMAG(t);
  double kcd = -2*tr/(mass*M_PI*(tr*tr + ti*ti));
  double k2l = pow(k, 2.0*l);
  return k2l*kcd;
}
