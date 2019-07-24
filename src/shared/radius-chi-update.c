#include <stdio.h>
#include <math.h>


void rmv_zero(double *R, double *R_last, int k, int n)
{
  int i;

  double Rtemp = R[k];
  double R_lasttemp = R_last[k];
  
  // Shift all Rs down so that R[k+1]->R[k], ..., //
  for (i = k; i < n-1; i++) {
    R[i] = R[i+1];
    R_last[i] = R_last[i+1];
  }
  // Save the kth radius in the n-1 spot, R[k]->R[n-1]. //
  R[n-1] = Rtemp;
  R_last[n-1] = R_lasttemp;

  return;
}


double dtChoose(double *R, double *B,double dt_max, int *k,int n)
// Choose the timestep dt such that only a single droplet shrinks to zero! //
// This is done to conserve volume fractions (by avoiding negative radii). //
// Only applicable for 3d case. Also, k is pointing to the radius value    //
// which is to be shrunk in the next time step choosen by the function     //
// (assuming the timestep is < dt_max).
{

  double dt_function(double R_i, double B_i);

  int i,j;
  double m1;
  double m2;

  // ignore B[i] >=0, since that means dt < 0 required to shrink R[i] //
  i = 0;
  do {
    i+= 1;
    if (i == n+1) {
      *k = 2*n;
      printf("k = %d\n n = %d\n",*k,n);
      return dt_max;
    }
  }
  while (B[i-1] >= 0); 

  m1 = dt_function(R[i-1],B[i-1]);

  *k = i-1;

  for (j = i; j < n; j++) {
    if (B[j] < 0) {
      m2 = dt_function(R[j],B[j]);
      if (m2 < m1) {
	m1 = m2;
	*k = j;
      }
    }
  }

  if (m1 > dt_max) {
    m1 = dt_max;
    *k = 2*n;;
  }
  return m1;
}

double dt_function(double R_i, double B_i)
// Function to minimize in dtChoose. Only applicable for //
// 3d case and B_i<0. //
{
  return -1*R_i*R_i*R_i/(3*B_i);
}



void updateRadius_and_chi(double *R, double *R_last, double *B,
			  double dt,double *chi,double alpha,
			  double beta,double L_sys,int k,int n)
{
  double chi_calculation(double alpha,double beta,double dropvolume,
			 double V);

  int i;
  double sum1 = 0;

  for (i = 0; i < n; i++) {
    R_last[i] = R[i];
    if (i==k) R[i] = 0;
    else R[i] = cbrt(R[i]*R[i]*R[i] + dt*3*B[i]);
    sum1 += R[i]*R[i]*R[i];
  }

  *chi = chi_calculation(alpha,beta,4*M_PI/3.0*sum1,
			 L_sys*L_sys*L_sys);

  return;
}

double chi_calculation(double alpha,double beta,double dropvolume,
		       double V)
{
  double phi = dropvolume/V;

  return (alpha-phi)*beta/(1-phi)-1;
}


