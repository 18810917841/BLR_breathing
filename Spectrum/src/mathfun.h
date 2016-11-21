#ifndef MATHFUN_H
#define MATHFUN_H

#include <gsl/gsl_histogram.h>

double *array_malloc(int n);
double ** matrix_malloc(int n1, int n2);
void matrix_free(double **p,int n);

void m_t(double *A,double *B,int n);
void m_i(double * a,int n);

void m_m(double *A,double *B,double *C,int n);
void m_v(double *A,double *B,double *C,int n);
double v_v(double *A,double *B,int n);
double v_m_v(double *a,double *M,double *b,int n);
double det_mat(double *a, int n);
int Chol_decomp_L(double *a, int n);
void display_mat(double *a, int m, int n);

double *get_cov_matrix(double *theta, int nstep, int ntheta);
void get_histogram_val_error(gsl_histogram *hist, double *val, double *err_lo, double *err_hi);

double fmax(double a ,double b);
double fmin(double a ,double b);
double mod(double y, double x);
void wrap(double *x, double min, double max);
int mod_int(int y, int x);
double logdiffexp(double x1, double x2);
double logsumexp(double *x, int n);
#endif
