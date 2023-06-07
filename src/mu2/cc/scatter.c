// Copyright 2022 Daniel Odell

#include "scatter.h"

void augment_potential_matrix(
    double* vp_matrix,
    double q0,
    double* v_matrix,
    double* p,
    int np
) {
  int nq=np, nqp=nq+1;
  double x[nq];
  double y[nq];
  gsl_interp2d* w = gsl_interp2d_alloc(gsl_interp2d_bicubic, nq, nq);
  gsl_interp_accel* xacc = gsl_interp_accel_alloc();
  gsl_interp_accel* yacc = gsl_interp_accel_alloc();

  // copy pmesh points to double array
  for (int i = 0; i < nq; i++) {
    x[i] = p[i];
    y[i] = p[i];
  }
  gsl_interp2d_init(w, x, y, v_matrix, nq, nq);

  for (int i = 0; i < nq; i++) {
    for (int j = 0; j < nq; j++) {
      // copy v_matrix to the interior of vp_matrix
      vp_matrix[i*nqp+j] = v_matrix[i*nq+j];
    }
    // interpolate (extrapolate) the last column
    vp_matrix[i*nqp+nq] = gsl_interp2d_eval_extrap(w, x, y, v_matrix,
        p[i], q0, xacc, yacc);
    // interpolate (extrapolate) the last row
    vp_matrix[nq*nqp+i] = gsl_interp2d_eval_extrap(w, x, y, v_matrix, q0,
      p[i], xacc, yacc);
  }
  vp_matrix[nq*nqp+nq] = gsl_interp2d_eval_extrap(w, x, y, v_matrix, q0, q0,
    xacc, yacc);

  gsl_interp2d_free(w);
  gsl_interp_accel_free(xacc);
  gsl_interp_accel_free(yacc);
}

void kernel_matrix_gen(
    gsl_matrix_complex* kernel,
    double q0,
    double* vp_matrix,
    double* p,
    double* wp,
    int np,
    double qmax,
    double mass
) {
  int nq = np, nqp = nq+1;
  double qi, wi, wj, qj;
  double logdiff_numerical = 0.0;
  gsl_complex logdiff_analytical, logdiff;
  for (int i = 0; i < nq; i++) {
    qi = p[i];
    wi = wp[i];
    logdiff_numerical += wi*qi*qi/(q0*q0 - qi*qi);
  }
  logdiff_analytical = gsl_complex_sub(
      gsl_complex_rect(-qmax-q0/2*log(1-2*q0/(q0+qmax)), 0),
      gsl_complex_rect(0, M_PI*q0/2.0)
    );
  logdiff = gsl_complex_sub(
      logdiff_analytical,
      gsl_complex_rect(logdiff_numerical, 0)
  );
  for (int i = 0; i < nqp; i++) {
    for (int j = 0; j < nqp-1; j++) {
      qj = p[j];
      wj = wp[j];
      gsl_matrix_complex_set(kernel, i, j, 
        gsl_complex_rect((double)(i == j) -
        mass*wj*qj*qj*vp_matrix[i*nqp+j]/(q0*q0-qj*qj), 0)
      );
    }
    gsl_matrix_complex_set(kernel, i, nqp-1,
      gsl_complex_sub(
        gsl_complex_rect((double)(i == nq), 0),
        gsl_complex_mul(
          gsl_complex_rect(mass*vp_matrix[i*nqp+(nqp-1)], 0),
          logdiff
        )
      )
    );
  }
}

gsl_complex t_on_shell(
    double q0,
    double* v_matrix,
    double* p,
    double* wp,
    int np,
    double qmax,
    double mass
) {
  int nq=np, nqp=nq+1, s;
  // double* vp_matrix = (double*)(malloc(sizeof(double)*nqp*nqp));
  double vp_matrix[nqp*nqp];
  gsl_matrix_complex* kernel = gsl_matrix_complex_alloc(nqp, nqp);
  gsl_vector_complex* t      = gsl_vector_complex_alloc(nqp);
  gsl_vector_complex* vp     = gsl_vector_complex_alloc(nqp);
  gsl_permutation* perm = gsl_permutation_alloc(nqp);

  augment_potential_matrix(vp_matrix, q0, v_matrix, p, np);
  for (int i = 0; i < nqp; i++)
    gsl_vector_complex_set(vp, i, 
      gsl_complex_rect(vp_matrix[i*nqp+(nqp-1)], 0)
    );

  kernel_matrix_gen(kernel, q0, vp_matrix, p, wp, np, qmax, mass);

  gsl_linalg_complex_LU_decomp(kernel, perm, &s);
  gsl_linalg_complex_LU_solve(kernel, perm, vp, t);

  gsl_complex result = gsl_vector_complex_get(t, nqp-1);

  // free(vp_matrix);
  gsl_matrix_complex_free(kernel);
  gsl_vector_complex_free(t);
  gsl_vector_complex_free(vp);
  gsl_permutation_free(perm);

  return result;
}

void t_on_shell_ref(
    double* tr,
    double* ti,
    double q0,
    double* v_matrix,
    double* p,
    double* wp,
    int np,
    double qmax,
    double mass
) {
  gsl_complex t = t_on_shell(q0, v_matrix, p, wp, np, qmax, mass);
  *tr = GSL_REAL(t);
  *ti = GSL_IMAG(t);
}

double kcotdelta(
    double k,
    double* v_matrix,
    double* p, // momentum nodes
    double* wp, // momentum weights
    int np, // number of momentum nodes
    double qmax,
    int l,
    double mass
) {
  gsl_complex t = t_on_shell(k, v_matrix, p, wp, np, qmax, mass);
  double tr = GSL_REAL(t), ti = GSL_IMAG(t);
  double kcd = -2*tr/(mass*M_PI*(tr*tr+ti*ti));
  double k2l = pow(k, 2.0*l);
  return k2l*kcd;
}
