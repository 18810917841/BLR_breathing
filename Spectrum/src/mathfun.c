#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_cblas.h>
#include <lapacke.h>

#include "mathfun.h"

//linear algebra
double *array_malloc(int n)
{
    double *array;
    if(!(array=malloc(n*sizeof(double))))
    {
        printf("uable to allocate the matrix! \n");
        exit(-1);
    }
    return array;
}

double ** matrix_malloc(int n1, int n2)
{
    double ** mat;
    int i;
    if(!(mat = malloc(n1*sizeof(double*))))
    {
        fprintf(stderr, "Unable to allocate the matrix!\n");
        exit(-1);
    }
    for(i=0; i<n1; i++)
    {
        if(!(mat[i] = malloc(n2*sizeof(double))))
        {
            fprintf(stderr, "Unable to allocate the matrix!\n");
            exit(-1);
        }
    }
    return mat;
}

void matrix_free(double **p,int n)
{
  int i=0;
  for(i=0;i<n;i++)
    free(p[i]);
  free(p);
}

void m_t(double *A,double *B,int n) //matrix transpose
{
    int i,j;
    for(i=0;i<n;i++)
    	for(j=0;j<n;j++)
    		B[i*n+j]=A[j*n+i];
}

void m_i(double * a,int n) //matrix inverse
{
    int * ipiv;
    ipiv=malloc(n*sizeof(int));
    LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n, n, a, n, ipiv);
    LAPACKE_dgetri(LAPACK_ROW_MAJOR, n, a, n, ipiv);
    free(ipiv);
}

void m_m(double *A,double *B,double *C,int n) //matrix multiply matrix
{
    const  enum CBLAS_ORDER Order=CblasRowMajor;
    const  enum CBLAS_TRANSPOSE TransA=CblasNoTrans;
    const  enum CBLAS_TRANSPOSE TransB=CblasNoTrans;
    const int M=n;//A的行数，C的行数
    const int N=n;//B的列数，C的列数
    const int K=n;//A的列数，B的行数
    const float alpha=1;
    const float beta=0;
    const int lda=K;//A的列
    const int ldb=N;//B的列
    const int ldc=N;//C的列

    cblas_dgemm(Order, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
}

void m_v(double *A,double *B,double *C,int n)
{
    const  enum CBLAS_ORDER Order=CblasRowMajor;
    const  enum CBLAS_TRANSPOSE TransA=CblasNoTrans;
    const  enum CBLAS_TRANSPOSE TransB=CblasNoTrans;
    const int M=n;//A的行数，C的行数
    const int N=1;//B的列数，C的列数
    const int K=n;//A的列数，B的行数
    const float alpha=1;
    const float beta=0;
    const int lda=K;//A的列
    const int ldb=N;//B的列
    const int ldc=N;//C的列

    cblas_dgemm(Order, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
}

double v_v(double *A,double *B,int n)
{
    int k;
    double l;
    for(k=0,l=0;k<n;k++)
    	l=l+A[k]*B[k];
    return l;
}

double v_m_v(double *a,double *M,double *b,int n)
{
    double *c;
    double re;
    c=array_malloc(n);
    m_v(M,b,c,n);
    re=v_v(a,c,n);
    free(c);
    return re;
}

double det_mat(double *a, int n)
{
    int * ipiv;
    int * info;
    int i;
    double det;
    ipiv=malloc(n*sizeof(int));
    info=malloc(8*sizeof(int));
    *info=LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n, n, a, n, ipiv);
    if(*info!=0)
    {
        printf("Wrong!\n");
        exit(-1);
    }

    det = 1.0;
    for(i=0; i<n; i++)
    {
        det *= a[i*n+i];
        if (ipiv[i] != i)
        {
            ipiv[ipiv[i]] = ipiv[i];
            det = -det;
        }
    }
    det=fabs(det);
    free(ipiv);
    free(info);
    return det;
}

int Chol_decomp_L(double *a, int n)
{
    int i,j;
    char uplo = 'L';
    int *info;
    info=malloc(8*sizeof(int));
    *info=LAPACKE_dpotrf(LAPACK_ROW_MAJOR, uplo, n, a, n);
    if(*info<0)
    {
        fprintf(stderr, "The %d-th argument had an illegal value!\n", *info);
        return *info;
    }
    else if (*info>0)
    {
        fprintf(stderr, "The leading minor of order %d is not positive definite, and the factorization could not be completed.\n", *info);
        return *info;
    }
    for(i=0;i<n;i++)
    	for(j=0;j<n;j++)
    		if(!(a[i*n+j]<1e4 && a[i*n+j]>-1e4))
    			return -1;

    for(i=0;i<n;i++)
    	for(j=n-1;j>i;j--)
    		a[i*n+j] = 0.0;
    return 0;
}

void display_mat(double *a, int m, int n)
{
    printf("\n");
    int i, j;
    for(i=0;i<m;i++)
    {
        for(j=0;j<n;j++)
        {
            printf("%e\t", a[i*n+j]);
        }
        printf("\n");
    }
    printf("\n");
}

// statistics
double *get_cov_matrix(double *theta, int nstep, int ntheta)
{
    int i, j;
    double cov;
    double *cov_matrix;
    cov_matrix=array_malloc(ntheta*ntheta);
    for(i=0; i<ntheta; i++)
    {
        cov = gsl_stats_covariance(&theta[i*nstep], 1, &theta[i*nstep], 1, nstep);
        cov = cov>1.0e-6 ? cov : 1.0e-6;
        cov_matrix[i*ntheta + i] = cov;
        for(j=0; j<i; j++)
        {
            cov = gsl_stats_covariance(&theta[i*nstep], 1, &theta[j*nstep], 1, nstep);
            cov = cov>0.0001?cov:0.0;
            cov_matrix[i*ntheta + j] = cov_matrix[j*ntheta + i] = cov;
        }
    }
    return cov_matrix;
}


void get_histogram_val_error(gsl_histogram *hist, double *val, double *err_lo, double *err_hi)
{
    int idx, i, nh;
    double max;

    nh = gsl_histogram_bins(hist);
    idx = gsl_histogram_max_bin(hist);
    max = gsl_histogram_max_val(hist);
    *val = 0.5*(hist->range[idx] + hist->range[idx+1]);

    max *= exp(-0.5);

    for(i=idx-1; i>0; i--)
    {
        if(hist->bin[i]- max >= 0 && hist->bin[i-1] - max <=0 )
        {
            *err_lo = 0.25 * ( hist->range[i-1] + hist->range[i+1] + 2.0*hist->range[i] );
            break;
        }
    }

    if(i==0)
    *err_lo = hist->range[0];

    for(i=idx; i<nh-1; i++)
    {
        if(hist->bin[i]-max >= 0 && hist->bin[i+1] - max <=0 )
        {
            *err_hi = 0.25 * ( hist->range[i] + hist->range[i+2] + 2.0*hist->range[i+1] );
            break;
        }
    }
    if(i==nh-1)
    *err_hi = hist->range[nh];
}

// useful functions
double fmax(double a ,double b)
{
    return a>b?a:b;
}

double fmin(double a ,double b)
{
    return a<b?a:b;
}

double mod(double y, double x)
{
  if(x <= 0)
  printf("Warning in mod(double, double) (Utils.cpp)");
  return (y/x - floor(y/x))*x;
}

void wrap(double *x, double min, double max)
{
  *x = mod(*x - min, max - min) + min;
}

int mod_int(int y, int x)
{
  if(y >= 0)
    return y - (y/x)*x;
  else
    return (x-1) - mod_int(-y-1, x);
}

double logdiffexp(double x1, double x2) // x1 is larger
{
  double biggest = x1;
  double xx1 = x1 - biggest, xx2 = x2 - biggest;
  return log(exp(xx1) - exp(xx2)) + biggest;
}

double logsumexp(double *x, int n)
{
  int j;
  double sum, max;
  max = x[0];
  for(j = 0; j < n; j++)
    max = fmax(max, x[j]);
  sum = 0.0;
  for(j=0; j< n; j++)
    sum += exp( x[j] - max);
  return log(sum) + max;
}
