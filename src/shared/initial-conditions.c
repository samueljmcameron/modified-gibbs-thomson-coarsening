#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


#define MCHN_ERR 1e-15

void gaussian_distribution(double *R, double R_avg, double sigma_Ravg,int N,gsl_rng *RNG)
{

  void delta_distribution(double *R,double R_avg,int N);

  if (fabs(sigma_Ravg)<MCHN_ERR) {

    delta_distribution(R,R_avg,N);

  } else {
  
    for (int i = 0; i < N; i++) {
      R[i] = gsl_ran_gaussian(RNG,sigma_Ravg)+R_avg;
      while (R[i] <= 0) R[i] = gsl_ran_gaussian(RNG,sigma_Ravg)+R_avg;
    }
    
  }

  return;
}

void delta_distribution(double *R,double R_avg,int N)
{
  int i;

  for (i = 0; i < N; i++) R[i] = R_avg;

  return;
}


void nucleatedDistribution(double *R,double R_avg,int N,gsl_rng *RNG)
{

  double maximum_ofG(double x0,double xavg,int Npoints);
  double G(double x);


  int i;
  int Found;
  double z, u;
  double xavg = 1.0,x0 = 1.5;
  double z0 = x0/xavg;
  double c = 2*z0*xavg*maximum_ofG(x0,xavg,1000);
  double u_dist = 1./1.5;

  for (i = 0; i < N; i++) {
    Found = 0;
    while (Found != 1) {
      z = 1.5*gsl_rng_uniform(RNG);
      u = gsl_rng_uniform(RNG);
      if (u <= xavg*G(xavg*z)/(c*u_dist)) {
	R[i] = z*R_avg;
	Found = 1;
      }
    }
  }
  return;
}


double maximum_ofG(double x0,double xavg,int Npoints)
{

  double G(double x);

  int i;
  double lower = 0, upper = x0-1e-5;
  double step = (upper-lower)/Npoints;
  double max1 = xavg*G(0);
  double max2;

  for (i = 1; i <= Npoints; i++) {
    max2 = xavg*G(lower + i*step);
    if (max2 > max1) max1 = max2;
  }
  return max1;
}


double G(double x)
{

  if (x > 1.5) return 0;

  double c = 3*(3*3*3)*2.718281828/3.174802104;

  return c*x*x*exp(1.5/(x-1.5))/(pow(1.5-x,11.0/3.0)*pow(x+3,7.0/3.0));
}
