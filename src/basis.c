#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


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
