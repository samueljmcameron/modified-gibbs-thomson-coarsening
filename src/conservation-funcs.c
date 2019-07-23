#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "headerfile.h"

double compute_alpha(double beta, double init_volFrac,double chi)
{
  return init_volFrac*(1-(1+chi)/beta)+(1+chi)/beta;
}

double change_in_alpha(double drop_volume,double chi, double alpha,
		       double beta,double L_sys)
{
  double V = L_sys*L_sys*L_sys;

  double new_alpha;

  new_alpha = (1+chi)/beta*(1-drop_volume/V)+drop_volume/V;

  return alpha-new_alpha;
}


double set_L_sys(double init_volFrac, double *R, int n)
{
  int i;
  double sum = 0;

  for (i = 0; i < n; i++) sum += R[i]*R[i]*R[i];

  return cbrt(4*M_PI/3.0*sum/init_volFrac);
}



double totalVolume(double *R, int n)
{
  int i;
  double sum = 0;

  for (i = 0; i < n; i++) sum += R[i]*R[i]*R[i];

  return sum*4.0*M_PI/3.0;

}

double volumeFraction_in_drops(double *R,double L_sys,int n)
{
  int i;
  double sum = 0;

  for (i = 0; i < n; i++) sum += R[i]*R[i]*R[i];

  return sum*4.0*M_PI/3.0/(L_sys*L_sys*L_sys);

}

