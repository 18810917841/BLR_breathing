#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <time.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_statistics.h>

#include "main.h"
#include "mathfun.h"

#define M_SUN 1.989e33
#define G 6.67259e-8
#define C 29979245800
#define SECOND_PER_DAY 86400
#define LIGHT_DAY (SECOND_PER_DAY*C)
#define SECOND_PER_YEAR (365.25*86400)
#define LIGHT_YEAR (SECOND_PER_YEAR*C)
#define PC (3.26*LIGHT_YEAR)

const gsl_rng_type *gsl_T;
gsl_rng *gsl_r;
static char BLR_breathing_file[100] = "./data/BLR_breathing";
static int thistask ,totaltask, namelen;

static double r_u = 20.0*LIGHT_DAY;
static double m_u = 1.0e8*M_SUN;
static double v_u;
static double t_u;

static double r_mean_0 = 1.0;
static double vphi_mean_0 = 1.0;
static double vr_mean_0 = 0.0;

static double s = 1.0;
static double alpha = 22.845/10;
static double Gamma_0 = 0.10;
static double sigma = 0.04;
static double tau = 12;
static double nu = 0.0;
static double beta = 0.003;
static double theta_open = M_PI/2;
static double theta_inc = 0.0;
static double r_min;
static double r_max;

static int n_cloud;
static const int n_cloud_pertask = 500;
static double t_min = 0.0;
static double t_max = 100.0;
static int t_step = 1000;
static double *PSmat = NULL;
static double *random_array = NULL;
static double *flux_array = NULL;

double Gamma_t(double t)
{
    int i = (t-t_min)*(t_step-1)/(t_max-t_min);
    return flux_array[i];
}

double weight(double r)
{
  double w = pow(r,2*s/3-2);
  if(r<r_min)
    return 0;
  else if(r>r_max)
    return 0;
  else
    return w;
}

double accer_r(double r,double vr,double t)
{
  double Gamma = Gamma_t(t);
  double a_r = 0.0;
  a_r += alpha * Gamma/pow(r,2-2*s/3); //Radiation Pressure
  a_r += -1/pow(r,2); //Gravitation
  a_r += -beta * pow(r,nu+s/3)*vr; //Drag Force
  return a_r;
}

double accer_phi(double r, double vphi)
{
    double a_phi = 0;
    a_phi += -beta * pow(r,nu+s/3)*vphi; //Drag Force
    return a_phi;
}


int func (double t, const double w[], double f[], void *params)
{
  (void)(t); /* avoid unused parameter warning */
  double x=w[0];
  double vx=w[1];
  double y=w[2];
  double vy=w[3];
  double r = sqrt(x*x+y*y);
  double vr = (vx*x + vy*y)/r;
  double vphi = (vy*x - vx*y)/r;
  double a_r = accer_r(r,vr,t);
  double a_phi = accer_phi(r, vphi);
  f[0] = vx;
  f[1] = (a_r*x - a_phi*y)/r;
  f[2] = vy;
  f[3] = (a_r*y + a_phi*x)/r;
  return GSL_SUCCESS;
}

int main(int argc, char **argv)
{
  char proc_name[MPI_MAX_PROCESSOR_NAME];
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &thistask);
  MPI_Comm_size(MPI_COMM_WORLD, &totaltask);
  MPI_Get_processor_name(proc_name, &namelen);
  MPI_Barrier(MPI_COMM_WORLD);

  gsl_rng_env_setup();
  gsl_T=gsl_rng_default;
  gsl_r = gsl_rng_alloc(gsl_T);
  gsl_rng_set(gsl_r, time(NULL)+thistask);

  int i=0,j=0;
  int flag1,flag2;
  n_cloud = n_cloud_pertask * totaltask;
  v_u = sqrt(G*m_u/r_u);
  t_u = r_u/v_u;
  r_min = 2*G*m_u/C/C/r_u;
  r_max = 100*r_mean_0;
  if(thistask==0)
  {
    printf("r_u = %e light day\n",r_u/LIGHT_DAY);
    printf("t_u = %e year\n",t_u/SECOND_PER_YEAR);
    printf("v_u = %e km/s\n",v_u/1e5);
    printf("m_u = %e M_SUN\n",m_u/M_SUN);
  }

  FILE *fp,*fp_orbit[n_cloud_pertask];
  if(thistask==0)
  {
    fp = fopen(BLR_breathing_file,"w");
    for(j=0;j<n_cloud_pertask;j++)
    {
      char file_name_1[100];
      char file_name_2[100] = "./data/orbit_";
      sprintf(file_name_1,"%d",j);
      strcat(file_name_2,file_name_1);
      fp_orbit[j] = fopen(file_name_2,"w");
    }
  }

  double dt = (t_max-t_min)/(t_step-1);
  int *ncloud_array = malloc(t_step*sizeof(int));
  double *time_array = array_malloc(t_step);
  double *r_mean_array = array_malloc(t_step);
  double *r_median_array = array_malloc(t_step);
  double *lag_array = array_malloc(t_step);
  random_array = array_malloc(t_step);
  PSmat = array_malloc(t_step * t_step);
  flux_array = array_malloc(t_step);
  double **radius_array_2d = matrix_malloc(t_step,n_cloud_pertask);
  double **weight_array_2d = matrix_malloc(t_step,n_cloud_pertask);
  double **lag_array_2d = matrix_malloc(t_step,n_cloud_pertask);

  for(i=0;i<t_step;i++)
  {
    time_array[i] = t_min + dt*i;
    r_mean_array[i] = 0.0;
    r_median_array[i] = 0.0;
    lag_array[i] = 0.0;
    ncloud_array[i] = 0;
    random_array[i] = gsl_ran_gaussian(gsl_r, 1.0);
    for(j=0;j<n_cloud_pertask;j++)
    {
      radius_array_2d[i][j] = 0.0;
      lag_array_2d[i][j] = 0.0;
      weight_array_2d[i][j] = 0.0;
    }
  }

  for(i=0; i<t_step; i++)
    for(j=0; j<=i; j++)
    {
      PSmat[i*t_step+j] = sigma*sigma* exp (- pow (fabs(time_array[i] -time_array[j]) / tau, 1.0));
      PSmat[j*t_step+i] = PSmat[i*t_step+j];
    }
  Chol_decomp_L(PSmat, t_step);
  m_v(PSmat,random_array,flux_array,t_step);
  for(i=0; i<t_step; i++)
    flux_array[i] += Gamma_0;


  gsl_odeiv2_system sys = {func, NULL, 4, NULL};
  gsl_odeiv2_driver *d = NULL;
  for(j=0;j<n_cloud_pertask;j++)
  {
    if(thistask==0) printf("ncloud = %d\n",j);
    double phi_0 = 2.0*M_PI*gsl_rng_uniform(gsl_r);
    double l_phi = 2.0*M_PI*gsl_rng_uniform(gsl_r);
    double l_theta = acos(cos(theta_open)+(1-cos(theta_open))*gsl_rng_uniform(gsl_r));
    double r_0 = (1+gsl_ran_gaussian(gsl_r, 0.1))*r_mean_0;
    double vphi_0 = (1+gsl_ran_gaussian(gsl_r, 0.1))*vphi_mean_0;
    double vr_0 = gsl_ran_gaussian(gsl_r, 0.1);

    double x_0 = r_0*cos(phi_0);
    double y_0 = r_0*sin(phi_0);
    double vx_0 = -vphi_0 * sin(phi_0) + vr_0 * cos(phi_0);
    double vy_0 = vphi_0 * cos(phi_0) + vr_0 * sin(phi_0);

    d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd, 1e-6, 1e-6, 0.0);
    double w[4] = {x_0, vx_0, y_0, vy_0};
    double t = t_min;

    flag2 = 1;
    for (i = 0; i < t_step; i++)
    {
      flag1 = 1;
      double ti = time_array[i];
      int status = gsl_odeiv2_driver_apply (d, &t, ti, w);
      if (status != GSL_SUCCESS)
      {
        printf ("error, return value=%d\n", status);
        break;
      }
      double x,y,z,xb,yb,zb,r;

      x=w[0];y=w[2];z=0.0;
      r = sqrt(x*x+y*y);
      if(r>r_max || r<r_min)
      {
        flag1 = 0;
        if(flag2 && thistask)
        {
          printf("The cloud %d is out of BLR at r=%lf light_day,t = %lf year\n", \
          j, r*r_u/LIGHT_DAY,time_array[i]*t_u/SECOND_PER_YEAR);
        }
        flag2 = 0;
      }
      xb = cos(l_theta)*cos(l_phi) * x - sin(l_phi) * y + sin(l_theta)*cos(l_phi) * z;
      yb = cos(l_theta)*sin(l_phi) * x + cos(l_phi) * y + sin(l_theta)*sin(l_phi) * z;
      zb =-sin(l_theta) * x + cos(l_theta)* z;
      x=xb*cos(theta_inc)-zb*sin(theta_inc);
      y=yb;
      z=xb*sin(theta_inc)+zb*cos(theta_inc);

      //double vx,vy,vz,vxb,vyb,vzb;
      //vx=w[1];vy=w[3];vz=0.0;
      // vxb = cos(l_theta)*cos(l_phi) * vx - sin(l_phi) * vy + sin(l_theta)*cos(l_phi) * vz;
      // vyb = cos(l_theta)*sin(l_phi) * vx + cos(l_phi) * vy + sin(l_theta)*sin(l_phi) * vz;
      // vzb =-sin(l_theta) * vx + cos(l_theta)* vz;
      // vx=vxb*cos(theta_inc)-vzb*sin(theta_inc);
      // vy=vyb;
      // vz=vxb*sin(theta_inc)+vzb*cos(theta_inc);
      radius_array_2d[i][j] = r;
      lag_array_2d[i][j] = r -z;
      weight_array_2d[i][j] = weight(r);
      if(thistask==0)
        fprintf(fp_orbit[j],"%lf %lf %lf %lf %lf\n",ti,x,y,z,r);
      if(flag1)
        ncloud_array[i]++;
    }
    gsl_odeiv2_driver_free(d);
  }

  for(i=0;i<t_step;i++)
  {
    if(ncloud_array[i]==0)
      return 0;
    r_mean_array[i] = gsl_stats_wmean(weight_array_2d[i],1,radius_array_2d[i],1,n_cloud_pertask);
    lag_array[i] = gsl_stats_wmean(weight_array_2d[i],1,lag_array_2d[i],1,n_cloud_pertask);
    gsl_sort (radius_array_2d[i], 1, n_cloud_pertask);
    r_median_array[i] = gsl_stats_median_from_sorted_data(radius_array_2d[i],1,n_cloud_pertask);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  double *buf1 = array_malloc(t_step);
  double *buf2 = array_malloc(t_step);
  double *buf3 = array_malloc(t_step);
  MPI_Reduce(r_mean_array, buf1, t_step, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(r_median_array, buf2, t_step, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(lag_array, buf3, t_step, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if(thistask==0)
  {
    memcpy(r_mean_array,buf1,t_step*sizeof(double));
    memcpy(r_median_array,buf2,t_step*sizeof(double));
    memcpy(lag_array,buf3,t_step*sizeof(double));
    for(i=0;i<t_step;i++)
    {
        if(i % 2 == 0)
        {
            fprintf(fp, "%.5e %.5e %.5e %.5e %.5e\n", \
            time_array[i]*t_u/SECOND_PER_YEAR, \
            r_mean_array[i]*r_u/LIGHT_DAY/totaltask, \
            r_median_array[i]*r_u/LIGHT_DAY/totaltask, \
            lag_array[i]*r_u/LIGHT_DAY/totaltask,
            Gamma_t(time_array[i]));
        }
    }
  }
  MPI_Finalize();
  return 0;
}
