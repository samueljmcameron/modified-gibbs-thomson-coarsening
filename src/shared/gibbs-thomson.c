#include <stdio.h>


void set_GibbsThomson(double *c, double *R, double chi,
		      double Req, int n)
{

  double GT_eq(double R_i,double Req);

  int i;

  for (i = 0; i < n; i++) {
    c[i] = GT_eq(R[i],Req)-chi;
  }
  
  return;
}

double GT_eq(double R_i,double Req)
// GibbsThomson effect with an equilibrium radius, Req. //
// The form of the free energy is approximated via a    //
// Taylor series expansion near Req. gppReq:=g''(Req),  //
// and gReq := g(Req), where g(R) is the free energy    //
// per unit volume. //
{
  return (R_i-Req)*(5.0/3.0*R_i-Req);
}

