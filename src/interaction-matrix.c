#include <stdio.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_permutation.h>

void update_matrix(double *copy, double *coeff_matrix,double *R,
		   int n,int N)
{
  int i,j;
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      if (i == j) {
	copy[i*(N)+j] = -(1.0/R[i]);
      }
      coeff_matrix[i*(N)+j] = copy[i*(N)+j];
    }
  }
  return;
}

void set_LU(gsl_matrix_view *m, gsl_vector_view *c_vec, 
	    gsl_vector_view *B_vec, gsl_permutation **p,
	    double *c, double *B, double *coeff_matrix,int s, int n,int N)
{

  gsl_permutation *gsl_permutation_realloc (gsl_permutation *p,const size_t n);

  // Put B and c vectors into appropriate form for gsl library. //
  *B_vec = gsl_vector_view_array(B, (n));                         // B is TBD by solver
  *c_vec = gsl_vector_view_array(c,(n));

  // Set up and perform LU decomposition. //
  *m = gsl_matrix_view_array_with_tda(coeff_matrix, (n), (n),N);
  *p = gsl_permutation_realloc(*p,(n));
  gsl_linalg_LU_decomp(&(*m).matrix, *p, &s);

  return;
}

void rearrange_matrix(double *coeff_matrix, double *copy, double *basis_matrix,
		      int k, int lastn,int N)
{

  double offDiagonal(double L_sys,double *basis_matrix,int i, int j,int N);

  int i, j;

  for (i = k; i < lastn-1; i++) {
      for (j = 0; j < k; j++) {
	copy[i*(N)+j] = copy[(i+1)*(N)+j];
	//	coeff_matrix[i*(N)+j] = copy[i*(N)+j];
	copy[j*(N)+i] = copy[j*(N)+(i+1)];
	//coeff_matrix[j*(N)+i] = copy[j*(N)+i];
      }
      for (j = k; j < lastn-1; j++) {
	copy[i*(N)+j] = copy[(i+1)*(N)+(j+1)];
	//coeff_matrix[i*(N)+j] = copy[i*(N)+j];
      }
      basis_matrix[i*3] = basis_matrix[(i+1)*3];
      basis_matrix[i*3+1] = basis_matrix[(i+1)*3+1];
      basis_matrix[i*3+2] = basis_matrix[(i+1)*3+2];
  }

  for (i = 0; i < lastn-1; i++) {
    for (j = 0; j < lastn-1; j++) {
      coeff_matrix[i*(N)+j] = copy[i*(N)+j];
    }
  }

  return;
}

void set_matrices(double *coeff_matrix,double *copy,double *basis_matrix,
		  double *R,double L_sys, int N)
{

  double offDiagonal(double L_sys,double *basis_matrix,int i, int j,int N);

  int i,j;

  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++) {
      if (i == j) {
	coeff_matrix[i*(N)+j] = -1*(1./R[i]);
	copy[i*(N)+j] = coeff_matrix[i*(N)+j];
      }
      else if (i > j) {
	coeff_matrix[i*(N)+j] = coeff_matrix[j*(N)+i]; // symmetry of the matrix allows this
	copy[i*(N)+j] = coeff_matrix[i*(N)+j];
      }
      else {
	coeff_matrix[i*(N)+j] = offDiagonal(L_sys,basis_matrix,i,j,N);
	copy[i*(N)+j] = coeff_matrix[i*(N)+j];
      }
    }
  }
  return;
}


double offDiagonal(double L_sys,double *basis_matrix,int i, int j,int N)
{

  double distanceBasis(double *basis_matrix,int i1, int i2, int u,
		       int v, int w);

  return -1.0/(L_sys*distanceBasis(basis_matrix,i,j,0,0,0));
}



gsl_permutation *gsl_permutation_realloc (gsl_permutation *p,const size_t n)
{

  if (n == 0)
    {
      //      GSL_ERROR_VAL ("permutation length n must be positive integer",
      //	     GSL_EDOM, 0);
    }                                                         

  p->data = (size_t *) realloc (p->data,n * sizeof (size_t));

  if (p->data == 0)
    {
      free (p);         /* exception in constructor, avoid memory leak */
      //GSL_ERROR_VAL ("failed to reallocate space for permutation data",
      //	     GSL_ENOMEM, 0);
    }

  p->size = n;

  return p;
}
