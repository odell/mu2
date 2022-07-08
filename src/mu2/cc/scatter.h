// Copyright 2022 Daniel Odell

#define _USE_MATH_DEFINES
#include <math.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <stdio.h>

void augment_potential_matrix(
    double* vp_matrix,
    double q0,
    double* v_matrix,
    double* p,
    int np
);

void kernel_matrix_gen(
    gsl_matrix_complex* kernel,
    double q0,
    double* vp_matrix,
    double* p,
    double* wp,
    int np,
    double qmax,
    double mass
);

gsl_complex t_on_shell(
    double q0,
    double* v_matrix,
    double* p,
    double* wp,
    int np,
    double qmax,
    double mass
);

double kcotdelta(
    double k,
    double* v_matrix,
    double* p, // momentum nodes
    double* wp, // momentum weights
    int np, // number of momentum nodes
    double qmax,
    int l,
    double mass
);