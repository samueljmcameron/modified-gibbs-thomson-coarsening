#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include <stdbool.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_permutation.h>
#include <string.h>
#include "headerfile.h"

int main(int argc, char **argv)
{

  unsigned long int generate_seed(void);

  void set_params(char **args,struct params *p);

  double compute_alpha(double beta, double init_volFrac,double chi);


  void allocate_vectors(int N,int t_intervals,double **R,double **R_last,double **c,
			double **B,double **workingVector,double **coeff_matrix,
			double **copy,double **basis_matrix,double **tvals);

  void build_tvals(double *tvals, double t_final,int t_intervals);

  void gaussian_distribution(double *R, double R_avg, double sigma_R,int N,gsl_rng *RNG);

  void generateCubicBasis(double *basis_matrix,double *R,double L_sys,int N);
  
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


  // Set up random number generator //

  gsl_rng *RNG;

  const gsl_rng_type * T;

  gsl_rng_env_setup();       // e.g. "GSL_RNG_SEED=124" before ./ in cmd line changes seed from 0 to 124
  T = gsl_rng_default;
  RNG = gsl_rng_alloc (T);
  gsl_rng_set(RNG,generate_seed());

  // define all starting values. //

  struct params p;

  set_params(argv,&p);

  int n=p.N,lastn=p.N;                                     // current # of drops, n-1
  double chi = p.chi_0;                                    // supersaturation at time t

  double alpha=compute_alpha(p.beta,p.volFrac_0,p.chi_0); // dimensionless solute number




  double *R,*R_last;            // radius stuff
  double *B;                    // sources and sinks 
  double *tvals;                // list of tvals to save

  // stuff to invert the matrix //

  double *workingVector;
  double *coeff_matrix,*copy;  // coeff_matrix used to solve for Bs
  double *basis_matrix;        // positions of droplets (flattened matrix)
  double *c;                   // vector with Gibbs Thomson effect, and c[n] = 1

  // malloc stuff //

  allocate_vectors(p.N,p.t_intervals,&R,&R_last,&c,&B,&workingVector,
		   &coeff_matrix,&copy,&basis_matrix,&tvals);

  build_tvals(tvals,p.t_final,p.t_intervals); // list of times to save R distribution


  // structures for linear algebra stuff.

  gsl_matrix_view m, refine;
  gsl_vector_view work,c_vec,B_vec;
  int s;
  gsl_permutation *perm;

  perm = gsl_permutation_alloc(p.N);

  // initial conditions //
  // Create gaussian distribution of drops around Ravg //

  gaussian_distribution(R,p.R_avg0,p.sigma_R0,p.N,RNG);

  double L_sys;                                 // side length of system (V = L_sys^3)
  L_sys = set_L_sys(p.volFrac_0,R,p.N);
  printf("L_sys = %lf\n",L_sys);


  generateCubicBasis(basis_matrix,R,L_sys,p.N);
  
  // Calculate coeff_matrix, and a copy of coeff_matrix. //
  set_matrices(coeff_matrix,copy,basis_matrix,R,L_sys,p.N);


  FILE *chi_file;                 // supersaturation data


  int jset = 0;                     // number of R_files which have already been saved

  double t=0,dt=0;                  // time, dt, max allowed dt, time, final time
  double save_time=tvals[0];        // last time at which the R_files were saved.
  int overlap = false;              // whether droplets are overlapping or not



  initialize_file(&chi_file,argv[1],"chi_vs_t",p,overlap);

  fprintf(chi_file,"%.12e\t%.12e\n",t,chi);


  // Rs, matrix, and chi are all initialized. Time to run simulation! //

  clock_t begin2 = clock();

  printf("\n\n");

  double sum = 0;

  int drop_removed_index;


  while (t < p.t_final && n>1) {

    printf("t = %lf\n",t);
    // Calculate the c vector. //
    set_GibbsThomson(c,R,chi,p.R_eq,n);
    update_matrix(copy,coeff_matrix,R,n,p.N);
    
    // solve for the B coefficients (B[n] is the supersaturation) //
    set_LU(&m,&c_vec,&B_vec,&perm,c,B,coeff_matrix,s,n,p.N);
    gsl_linalg_LU_solve (&m.matrix, perm, &c_vec.vector, &B_vec.vector);
    refine = gsl_matrix_view_array_with_tda(copy,n,n,p.N);
    work = gsl_vector_view_array(workingVector,n);
    gsl_linalg_LU_refine(&refine.matrix,&m.matrix,perm,&c_vec.vector,&B_vec.vector,&work.vector);


    if (isnan(dt)) printf("dt is NAN\n");

    dt = dtChoose(R,B,p.dt_max,&drop_removed_index,n);

    updateRadius_and_chi(R,R_last,B,dt,&chi,alpha,p.beta,L_sys,drop_removed_index,n);
    fprintf(chi_file,"%.12e\t%.12e\n",t,chi);

    // save initial R distribution
    if (t == 0) {

      save_Rdistribution(argv[1],jset,R,R_last,dt,n,overlap,p);

      save_basis(argv[1],jset,basis_matrix,L_sys,n,overlap,p);

      jset += 1;
      save_time = tvals[jset];
    }
    t += dt;

    // next IF statement only executed if dt is such that a radius has shrunk. //
    if (drop_removed_index != 2*n) {
      rmv_zero(R,R_last,drop_removed_index,n);
      lastn = n;
      n = n-1;
      rearrange_matrix(coeff_matrix,copy,basis_matrix,drop_removed_index,lastn,p.N);
      sum = 0;
      for (int i = 0; i < n; i++) {
	sum += B[i];
      }
    }

    if (t > save_time) {



      //if (change_in_alpha(totalVolume(R,n),chi,alpha,p.beta,L_sys) != 0) {
      //printf("not zero!\n");
      //}

      // save when t reaches save_time


      save_Rdistribution(argv[1],jset,R,R_last,dt,n,overlap,p);

      save_basis(argv[1],jset,basis_matrix,L_sys,n,overlap,p);

      jset += 1;

      if (jset == p.t_intervals) save_time = 1e30;
      else save_time = tvals[jset];
    }

  }
  
  // save the final configuration //

  jset = p.t_intervals;

  save_Rdistribution(argv[1],jset,R,R_last,dt,n,overlap,p);

  save_basis(argv[1],jset,basis_matrix,L_sys,n,overlap,p);

  clock_t end2 = clock();
  
  printf("(alpha(t)-alpha(0))/alpha(0)=%6.6e\n",
	 change_in_alpha(totalVolume(R,n),chi,alpha,
			 p.beta,L_sys)/alpha);
  printf("time taken = %lf s\n",(double)(end2-begin2)/CLOCKS_PER_SEC);


  free_vectors(&R,&R_last,&c,&B,&workingVector,&coeff_matrix,&copy,&basis_matrix,&tvals);

  fclose(chi_file);
  
  gsl_permutation_free(perm);
  gsl_rng_free(RNG);
  return 0;
}
