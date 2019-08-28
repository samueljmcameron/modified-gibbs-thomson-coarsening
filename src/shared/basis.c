#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

void generateCubicBasis(double *basis_matrix,double *R,double L_sys,int N)
{

  bool checkOverlap(double *basis_matrix,int i1, int i2,double *R,
		    double L_sys,double buffer);
  
  double s0 = 0.0;
  double sf = 0.9;

  int Ncbrt = round(cbrt(N));

  if (Ncbrt*Ncbrt*Ncbrt != N) {
    printf("need a droplet number that has integer cube root!");
    exit(1);
  }
  
  double sh = (sf-s0)/(Ncbrt-1);

  for (int i = 0, num = 0; i < Ncbrt; i ++) {

    for (int j = 0; j < Ncbrt; j ++ ) {

      for (int k = 0; k < Ncbrt; k ++ ) {

	basis_matrix[num*3] = i*sh;
	basis_matrix[num*3+1] = j*sh;
	basis_matrix[num*3+2] = k*sh;

	for (int n = 0; n < num; n++) {
	  if (checkOverlap(basis_matrix,num,n,R,L_sys,1.0)) {
	    printf("overlapping at R[%d] = %e, R[%d] = %e!\n",num,R[num],
		   n,R[n]);
	  }
	}
	
	num += 1;

      }
    }
  }
  return;
}

/*
void generateHexagonalBasis(double *basis_matrix,double *R,double L_sys, double Ravg, int N)
{


  double ae0_x = 1.0;
  double a0_y = 0.0;
  double 
  double a1_x = cos(2*M_PI/3);
  double a1_y = sin(2*M_PI/3);

  int Ncbrt = round(cbrt(N));

  if (Ncbrt*Ncbrt*Ncbrt != N) {
    printf("need a droplet number that has integer cube root!");
    exit(1);
  }

  double a1x = 1;
  double a1y = 0;
  double a1z = 0;

  double a2x = 0.5;
  double a2y = sqrt(3)/2.0;
  double a2z = 0;

  double zup = sqrt(2.0/3.0);

  double bx = 0.5;
  double by = 0.5*sqrt(1.0/3.0);
  
  for (int n3 = 0, double cx, double cy,int i = 0; n3 < Ncbrt; n3 ++) {

    if (n3 % 2 != 0) {

      cx = bx;
      cy = by;

    } else {
      cx = 0;
      cy = 0;
    }

    for (int n2 = 0,int noff; n2 < Ncbrt; n2 ++) {

      noff = (n2-1)/2;

      for (int n1 = -noff; n1 < Ncbrt - noff; n1 ++) {

	basis_matrix[i*3] = n1*a1x+n2*a2x+n3*zup+cx;
	basis_matrix[i*3+1] = n1*a1y+n2*a2y+n3*zup+cy;
	basis_matrix[i*3+2] = gsl_rng_uniform(RNG);

	i ++ 1;
    
    for (int j = 0; j < i; j++) {
      if (checkOverlap(basis_matrix,i,j,R,L_sys,buffer)) {
	i -= 1;
	break;
      }
    }
  }
  return;
}
*/

void generateBasis(double *basis_matrix,double *R,double L_sys, double buffer, int N,
		   gsl_rng *RNG)
{


  bool checkOverlap(double *basis_matrix,int i1, int i2,double *R,
		    double L_sys,double buffer); 

  for (int i = 0; i < N; i++) {
    basis_matrix[i*3] = gsl_rng_uniform(RNG);
    basis_matrix[i*3+1] = gsl_rng_uniform(RNG);
    basis_matrix[i*3+2] = gsl_rng_uniform(RNG);
    
    for (int j = 0; j < i; j++) {
      if (checkOverlap(basis_matrix,i,j,R,L_sys,buffer)) {
	i -= 1;
	break;
      }
    }
  }
  return;
}

bool overlapping(double *basis_matrix,double *R, double L_sys, int n)
{

  bool checkOverlap(double *basis_matrix,int i1, int i2,double *R,
		    double L_sys,double buffer) ;

  bool overlap = false;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < i; j++) {
      if (checkOverlap(basis_matrix,i,j,R,L_sys,1)) {
	printf("overlapping particles at (%d,%d)!!\n",i,j);
	overlap = true;
      }
    }
  }
  return overlap; 
}


bool checkOverlap(double *basis_matrix,int i1, int i2,double *R,
		  double L_sys,double buffer) 
{

  double distanceBasis(double *basis_matrix,int i1, int i2, int u,
		       int v, int w);

  if (distanceBasis(basis_matrix,i1,i2,0,0,0)*L_sys/buffer<=R[i1]+R[i2]) return true;
  else if (distanceBasis(basis_matrix,i1,i2,1,0,0)*L_sys/buffer<=R[i1]+R[i2]) return true;
  else if (distanceBasis(basis_matrix,i1,i2,1,1,0)*L_sys/buffer<=R[i1]+R[i2]) return true;
  else if (distanceBasis(basis_matrix,i1,i2,0,1,0)*L_sys/buffer<=R[i1]+R[i2]) return true;
  else if (distanceBasis(basis_matrix,i1,i2,0,0,1)*L_sys/buffer<=R[i1]+R[i2]) return true;
  else if (distanceBasis(basis_matrix,i1,i2,1,0,1)*L_sys/buffer<=R[i1]+R[i2]) return true;
  else if (distanceBasis(basis_matrix,i1,i2,0,1,1)*L_sys/buffer<=R[i1]+R[i2]) return true;
  else if (distanceBasis(basis_matrix,i1,i2,1,1,1)*L_sys/buffer<=R[i1]+R[i2]) return true;
  else if (distanceBasis(basis_matrix,i2,i1,1,0,0)*L_sys/buffer<=R[i1]+R[i2]) return true;
  else if (distanceBasis(basis_matrix,i2,i1,1,1,0)*L_sys/buffer<=R[i1]+R[i2]) return true;
  else if (distanceBasis(basis_matrix,i2,i1,0,1,0)*L_sys/buffer<=R[i1]+R[i2]) return true;
  else if (distanceBasis(basis_matrix,i2,i1,0,0,1)*L_sys/buffer<=R[i1]+R[i2]) return true;
  else if (distanceBasis(basis_matrix,i2,i1,1,0,1)*L_sys/buffer<=R[i1]+R[i2]) return true;
  else if (distanceBasis(basis_matrix,i2,i1,0,1,1)*L_sys/buffer<=R[i1]+R[i2]) return true;
  else if (distanceBasis(basis_matrix,i2,i1,1,1,1)*L_sys/buffer<=R[i1]+R[i2]) return true;
  else return false;
}

double distanceBasis(double *basis_matrix,int i1, int i2, int u,
		     int v, int w)
{
  double delx = basis_matrix[i1*3]+u-basis_matrix[i2*3];
  double dely = basis_matrix[i1*3+1]+v-basis_matrix[i2*3+1];
  double delz = basis_matrix[i1*3+2]+w-basis_matrix[i2*3+2];

  return sqrt(delx*delx+dely*dely+delz*delz);
}
