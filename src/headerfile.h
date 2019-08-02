#ifndef _HEADERFILE_H_
#define _HEADERFILE_H_

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include <stdbool.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_permutation.h>

struct params{

  double R_avg0;
  double sigma_R0;
  double R_eq;
  double volFrac_0;
  double beta;
  double chi_0;
  int N;
  double dt_max;
  double t_final;
  int t_intervals;
};

unsigned long int generate_seed(void);

void set_params(char **args,struct params *p);

double compute_alpha(double beta, double init_volFrac,double chi);


void allocate_vectors(int N,int t_intervals,double **R,double **R_last,double **c,
		      double **B,double **workingVector,double **coeff_matrix,
		      double **copy,double **basis_matrix,double **tvals);

void build_tvals(double *tvals, double t_final,int t_intervals);

void gaussian_distribution(double *R, double R_avg, double sigma_R,int N,gsl_rng *RNG);


void generateBasis(double *basis_matrix,double *R,double L_sys, double buffer, int N,
		   gsl_rng *RNG);

double set_L_sys(double init_volFrac, double *R, int n);

void set_matrices(double *coeff_matrix,double *copy,double *basis_matrix,
		  double *R,double L_sys, int N);


void initialize_file(FILE **output,char *path,char *fname,struct params p,
		     bool overlap);

void set_GibbsThomson(double *c, double *R, double chi,
		      double Req, int n);

void update_matrix(double *copy, double *coeff_matrix,double *R,
		   int n,int N);

void set_LU(gsl_matrix_view *m, gsl_vector_view *c_vec, 
	    gsl_vector_view *B_vec, gsl_permutation **p,
	    double *c, double *B, double *coeff_matrix,int s, int n,int N);

double dtChoose(double *R, double *B,double dt_max, int *k,int n);

void updateRadius_and_chi(double *R, double *R_last, double *B,
			  double dt,double *chi,double alpha,
			  double beta,double L_sys,int k,int n);

void save_basis(char *path,int jset,double *basis_matrix,double L_sys,
		int n,bool overlap,struct params p);

void save_Rdistribution(char *path,int jset,double *R,double *R_last,
			double dt,int n,bool overlap,struct params p);

void rmv_zero(double *R, double *R_last, int k, int n);

void rearrange_matrix(double *coeff_matrix, double *copy, double *basis_matrix,
		      int k, int lastn,int N);

double change_in_alpha(double drop_volume,double chi, double alpha,
		       double beta,double L_sys);

double totalVolume(double *R, int n);

bool overlapping(double *basis_matrix,double *R, double L_sys, int n);

void free_vectors(double **R,double **R_last,double **c,double **B,
		  double **workingVector,double **coeff_matrix,double **copy,
		  double **basis_matrix,double **tvals);




#endif
