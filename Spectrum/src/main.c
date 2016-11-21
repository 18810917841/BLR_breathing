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
static char flie_line_profile[100] = "./data/line_profile_0.1";
static int thistask ,totaltask, namelen;

static double r_u = 20*LIGHT_DAY;
static double m_u = 1.0e8;
static double v_u;
static double t_u;

static double v_0 = 0.5;
static double vr_0 = 0.0;
static double r_23 = 10.0;
static double s = 1.2;
static double Gamma_0 = 0.1;
static double alpha_0 = 0.5;
static double theta_open = M_PI/2;
static double theta_inc = 0.0;
static double c1 = 1.14e-11;
static double c2 = 8.8e-13;

static int n_cloud_pertask = 100;
static double t_min = 0.0;
static double t_max = 200.0;
static double v_min = -3.0;
static double v_max = 3.0;
static int t_step = 1000;
static int v_step = 100;

double alpha_r(double r)
{
  return alpha_0;
}

double Gamma_t(double t)
{
  return Gamma_0;
}

double accer(double r,double t)
{
  double alpha = alpha_r(r);
  double Gamma = Gamma_t(t);
  double beta = Gamma*c1/c2*alpha/pow(r/r_23,-2*s/3);
  return (beta-1)*1/r/r;
}

int func (double t, const double w[], double f[], void *params)
{
  (void)(t); /* avoid unused parameter warning */
  double x=w[0];
  double vx=w[1];
  double y=w[2];
  double vy=w[3];
  double r = sqrt(x*x+y*y);
  double a = accer(r,t);
  f[0] = vx;
  f[1] = a*x/r;
  f[2] = vy;
  f[3] = a*y/r;
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
  gsl_odeiv2_system sys = {func, NULL, 4, NULL};
  gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd, 1e-6, 1e-6, 0.0);

  int i=0,j=0;
  v_u = sqrt(G*m_u*M_SUN/r_u);
  t_u = r_u/v_u;
  printf("r_u = %e light day\n",r_u/LIGHT_DAY);
  printf("t_u = %e year\n",t_u/SECOND_PER_YEAR);
  printf("v_u = %e km/s\n",v_u/1e5);
  printf("m_u = %e M_SUN\n",m_u);
  // theta_open = theta_open / 180.0 * M_PI;
  // theta_inc = theta_inc / 180.0 * M_PI;
  FILE *fp;
  if(thistask==0)
  {
    fp = fopen(flie_line_profile,"w");
  }
  double dv = (v_max-v_min)/(v_step-1);
  double dt = (t_max-t_min)/(t_step-1);
  double *line_profile_velocity = array_malloc(v_step);
  double *line_profile_flux = array_malloc(v_step);
  for(i=0;i<v_step;i++)
  {
    line_profile_velocity[i] = v_min + dv*i;
    line_profile_flux[i] = 0.0;
  }
  for(i=0;i<n_cloud_pertask;i++)
  {
    if(thistask==0) printf("n_cloud = %d\n",i);
    double phi_0 = 2.0*M_PI*gsl_rng_uniform(gsl_r);
    double l_phi = 2.0*M_PI*gsl_rng_uniform(gsl_r);
    double l_theta = acos(cos(theta_open)+(1-cos(theta_open))*gsl_rng_uniform(gsl_r));
    //double l_theta = acos(gsl_rng_uniform(gsl_r)*cos(theta_open));
    double x_0 = cos(phi_0);
    double y_0 = sin(phi_0);
    double vx_0 = -v_0 * sin(phi_0) + vr_0 * cos(phi_0);
    double vy_0 = v_0 * cos(phi_0) + vr_0 * sin(phi_0);
    double w[4] = {x_0, vx_0, y_0, vy_0};
    for (j = 0; j <= t_step; j++)
    {
      double tj = t_min + dt*j;
      int status = gsl_odeiv2_driver_apply (d, &t_min, tj, w);
      if (status != GSL_SUCCESS)
      {
        printf ("error, return value=%d\n", status);
        break;
      }
      double x,y,z,vx,vy,vz,r;
      double xb,yb,zb,vxb,vyb,vzb;
      x=w[0];y=w[2];z=0.0;
      vx=w[1];vy=w[3];vz=0.0;
      r = sqrt(x*x+y*y);

      xb = cos(l_theta)*cos(l_phi) * x - sin(l_phi) * y + sin(l_theta)*cos(l_phi) * z;
      yb = cos(l_theta)*sin(l_phi) * x + cos(l_phi) * y + sin(l_theta)*sin(l_phi) * z;
      zb =-sin(l_theta) * x + cos(l_theta)* z;
      x=xb*cos(theta_inc)-zb*sin(theta_inc);
      y=yb;
      z=xb*sin(theta_inc)+zb*cos(theta_inc);

      vxb = cos(l_theta)*cos(l_phi) * vx - sin(l_phi) * vy + sin(l_theta)*cos(l_phi) * vz;
      vyb = cos(l_theta)*sin(l_phi) * vx + cos(l_phi) * vy + sin(l_theta)*sin(l_phi) * vz;
      vzb =-sin(l_theta) * vx + cos(l_theta)* vz;
      vx=vxb*cos(theta_inc)-vzb*sin(theta_inc);
      vy=vyb;
      vz=vxb*sin(theta_inc)+vzb*cos(theta_inc);

      double V = -vz;
      if(V<v_min || V>=v_max+dv)
        continue;
      int idv= (V - v_min)/dv;
      line_profile_flux[idv] += pow(r,2*s/3-2.0);
    }
  }
  gsl_odeiv2_driver_free (d);
  MPI_Barrier(MPI_COMM_WORLD);
  double *buf = array_malloc(v_step);
  MPI_Reduce(line_profile_flux, buf, v_step, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if(thistask==0)
  {
    memcpy(line_profile_flux,buf,v_step*sizeof(double));
    double scale=0.0;
    for(i=0;i<v_step;i++)
    {
      if(scale<line_profile_flux[i])
        scale = line_profile_flux[i];
    }
    scale = 1.0/scale;
    for(i=0;i<v_step;i++)
    {
      line_profile_flux[i] *= scale;
      fprintf(fp, "%.5e %.5e \n",line_profile_velocity[i],line_profile_flux[i]);
    }
  }
  return 0;
}
