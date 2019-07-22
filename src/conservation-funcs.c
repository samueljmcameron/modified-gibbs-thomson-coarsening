#include <stdio.h>
#include <stdlib.h>
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
