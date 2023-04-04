// Copyright 2022 Daniel Odell

#define _USE_MATH_DEFINES
#include <math.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <stdio.h>
#include "scatter.h"


/*
First-order perturbation functions
*/
void driving_term_pert1(
    gsl_vector_complex* d1, // NLO driving term
    double k,
    double* v1, // NLO interaction (exclusive)
    gsl_vector_complex* t_on_shell, // LO, half-shell t-matrix elements, t^{(0)}
    double* p, // momentum nodes
    double* wp, // momentum weights
    int np, // number of momentum nodes
    double qmax,
    double mass
);

gsl_complex t_on_shell_pert1(
    double k,
    double* v0, // NLO interaction (exclusive)
    double* v1, // NLO interaction (exclusive)
    double* p, // momentum nodes
    double* wp, // momentum weights
    int np, // number of momentum nodes
    double qmax,
    double mass
);

void t_pert1_sum(
    double* t_real, // Re(t_0 + t_1)
    double* t_imag, // Im(t_0 + t_1)
    double k, // scattering momentum
    double* v0, // NLO interaction (exclusive)
    double* v1, // NLO interaction (exclusive)
    double* p, // momentum nodes
    double* wp, // momentum weights
    int np, // number of momentum nodes
    double qmax,
    double mass
);

double kcotdelta_pert1(
    double k,
    double* v0, // LO interaction
    double* v1, // NLO interaction (exclusive)
    double* p, // momentum nodes
    double* wp, // momentum weights
    int np, // number of momentum nodes
    double qmax,
    int l,
    double mass
);
