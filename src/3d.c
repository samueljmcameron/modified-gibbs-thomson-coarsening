#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_integration.h>
#include <string.h>

gsl_rng *RNG;

/*#define COEFF_MATRIX(i,j) (coeff_matrix[(N+1)*(i)+(j)])
#define COPY(i,j) (copy[(N+1)*(i)+(j)])
#define BASIS_MATRIX(i,j) (basis_matrix[(2)*(i)+(j)])
*/

#define T_INTERVALS 19              // how many different time values (t) are saved to files.
#define NAME_LENGTH 200           // space for input file names


// saving stuff functions //
void directory_build(int argnum,double *Ravg,double *Rsigma,
		     double *Req,double *dt_max,double *init_volFrac,
		     double *beta, double *chi,int *N,
		     char **prefix,char **path,char *argv2,char *argv3);
void save_BASIS(char **fname,char **path,char **prefix,
		FILE *basis_files_JSET,int runNumber,int JSET,
		double *basis_matrix,double L_sys, int n,
		double current_t_val,int OVERLAP);
void save_Rdistribution(char **fname,char **path,char **prefix,
			FILE *R_files_JSET,int runNumber, int JSET,
			double *R,double *R_last,double dt,int n,
			double current_t_val,int OVERLAP);
void make_chi_file(char **fname,char **path, char **prefix,
		   FILE **chi_file,double Ravg,double Rsigma,
		   double Req, double dt_max, double init_volFrac,
		   double alpha,double beta, double chi,double L_sys,
		   int N,int runNumber,unsigned long seed);
void TVALS_2_save(double *tvals, double last_tval);

// initialization functions //
void onevalue_distribution(double *R,double Ravg,int N);
void gaussianDistribution(double *R, double Ravg, double Rsigma,int N);
void nucleatedDistribution(double *R,double Ravg,int N);
void generateBasis(double *basis_matrix,double *R,double L_sys, double buffer, int N);
double set_L_sys(double init_volFrac, double *R, int n);
double G(double x);
double maximum_ofG(double x0,double xavg,int Npoints);
double distanceBasis(double *basis_matrix,int i1, int i2, int u,
		     int v, int w);
int checkOverlap(double *basis_matrix,int i1, int i2,double *R,double L_sys,double buffer);

// calculating volume, volume fraction //
double totalVolume(double *R, int n);
double volumeFraction_in_drops(double *R,double L_sys,int n);

// calculating initial matrix, used to solve for B[i]s //
void set_matrices(double *coeff_matrix,double *copy,double *basis_matrix,double *R,
		  double L_sys,int N);
double Delement(double L_sys);
double offDiagonal(double L_sys,double *basis_matrix,int i, int j,int N);
double exp_func1(double x, void *params);
double exp_func2(double x, void *params);
double Dterm1(int h, int k,int l);
double Dterm2(int h, int k,int l);
double complex offDiagterm1(int h, int k, int l,double a_i0,double a_i1,
			    double a_i2, double a_j0, double a_j1, double a_j2);
double offDiagterm2(int h, int k, int l, double a_i0, double a_i1,double a_i2,
		    double a_j0,double a_j1, double a_j2);

// setting up to solve for the growth and shrink coefficients B[i] //
void rearrange_matrix(double *coeff_matrix,double *copy, double *basis_matrix,
		      int k, int lastn,int N);
void update_matrix(double *copy, double *coeff_matrix,double *R,
		   int n,int N);
void set_GibbsThomson(double *c, double *R, double chi,
		      double Req, int n);
void set_LU(gsl_matrix_view *m, gsl_vector_view *c_vec, gsl_vector_view *B_vec,
	    gsl_permutation **p, double *c, double *B, double *coeff_matrix,
	    int s, int n,int N);
double GT_eq(double R_i,double Req);
gsl_permutation *gsl_permutation_realloc(gsl_permutation *p,const size_t n);


// updating the R[i]s //
double dtChoose(double *R, double *B, double dt_max,int *k,int n);
void updateRadius_and_chi(double *R, double *R_last, double *B,
			  double dt,double *chi,double alpha,
			  double beta,double L_sys,int k,int n);
double dt_function(double R_i, double B_i);

// change in the supersaturation and conservation stuff //
double compute_alpha(double beta, double init_volFrac,double chi);
double change_in_alpha(double drop_volume,double chi, double alpha,
		       double beta,double L_sys);
double chi_calculation(double alpha,double beta,double dropvolume,
		       double V);
// remove all shrunk R[i]s //
void rmv_zero(double *R,double *R_last, int k, int n);

// conditions for loop to keep going //
int overlapping(double *basis_matrix,double *R, double L_sys, int n);



int main(int argc, char **argv)
{
  // Set up random number generator //
  const gsl_rng_type * T;

  gsl_rng_env_setup();       // e.g. "GSL_RNG_SEED=124" before ./ in cmd line changes seed from 0 to 124
  T = gsl_rng_default;
  RNG = gsl_rng_alloc (T);

  unsigned long int randval;
  FILE *seedFile;
  seedFile = fopen("/dev/urandom","r");
  fread(&randval,sizeof(randval),1,seedFile);
  fclose(seedFile);
  randval = 13185379641497912115;
  gsl_rng_set(RNG,randval);
  printf("randval = %lu\n",randval);
  printf("first value = %lu\n",gsl_rng_get(RNG));


  // define all stuff in the system //
  
  int i,k;
  double Ravg,Rsigma,Req,*R,*R_last;            // radius stuff
  double *B;                                    // sources and sinks 
  int runNumber = atof(argv[1]);
  int N,n,lastn;                                // init # of drops, # of drops, n-1
  double t=0,dt,dt_max,current_t_val,last_tval; // time, dt, max allowed dt, time, final time
  double init_volFrac,chi;                      // droplet volfrac, chi 
  double alpha,beta;                            // dimensionless quantities from conservation laws
  double L_sys;                                 // side length of system (V = L_sys^3)
  int OVERLAP = 0;                              // whether droplets are overlapping or not



  // stuff to invert the matrix //

  double *workingVector;
  double *coeff_matrix,*copy;  // coeff_matrix used to solve for Bs
  double *basis_matrix;        // positions of droplets (flattened matrix)
  double *c;                   // vector with Gibbs Thomson effect, and c[n] = 1
  gsl_matrix_view m, refine;
  gsl_vector_view work,c_vec,B_vec;
  int s;
  gsl_permutation *p;

  // file saving stuff //

  FILE *R_files[T_INTERVALS];     // array of files to save R at certain times
  FILE *R_final_file;             // final file to save R at end time
  FILE *basis_files[T_INTERVALS]; // array of files to save 
  FILE *basis_final_file;         // final file to save.
  FILE *chi_file;                 // supersaturation data
  double *tvals;                  // list of tvals to save
  char **path,**fname,**prefix;   // path name, file names TBD below
  int JSET = 0;                   // number of R_files which have already been saved

  path = malloc(sizeof(char*));
  prefix = malloc(sizeof(char*));
  path[0] = malloc(sizeof(char)*NAME_LENGTH);
  prefix[0] = malloc(sizeof(char)*NAME_LENGTH);
  fname = malloc(sizeof(char*)*2);
  fname[0] = malloc(sizeof(char)*NAME_LENGTH);
  fname[1] = malloc(sizeof(char)*NAME_LENGTH);
  

  // define all starting values, and input user chooses which one to change from default //

  dt = 0.0;
  Ravg = 6.0;
  Rsigma = 0.1;
  Req = 10.0;
  dt_max = 0.1;
  init_volFrac = 0.01;
  beta = 8.0;
  chi = 0.01;
  N = 1000;

  directory_build(argc,&Ravg,&Rsigma,&Req,
		  &dt_max,&init_volFrac,&beta,&chi,
		  &N,prefix,path,argv[2],argv[3]);

  n = N;
  lastn = N;
  last_tval = 4000.0;//11.72;// Ravg*Ravg*Ravg;
  alpha = compute_alpha(beta,init_volFrac,chi);

  printf("alpha = %lf\n",alpha);
  printf("beta = %lf\n",beta);
  printf("prefix = %s, path = %s\n",prefix[0],path[0]);

  // malloc stuff //

  p = gsl_permutation_alloc((N));
  R = malloc(N*sizeof(double));
  R_last = malloc(N*sizeof(double));
  c = malloc((N)*sizeof(double));
  B = malloc((N)*sizeof(double));
  workingVector = malloc((N)*sizeof(double));
  coeff_matrix = malloc(sizeof(double)*(N)*(N));
  copy = malloc(sizeof(double)*(N)*(N));
  basis_matrix = malloc(sizeof(double)*N*3);
  tvals = malloc(T_INTERVALS*sizeof(double));

  TVALS_2_save(tvals,last_tval); // list of times to save R distribution


  // initial conditions //

  clock_t begin1 = clock(); 
  
  // Create nucleated distribution of Rvalues with mean Ravg, std deviation Rsigma //
  //  nucleatedDistribution(R,Ravg,N);

  // Create gaussian distribution of drops around Ravg //
  gaussianDistribution(R,Ravg,Rsigma,N);

  L_sys = set_L_sys(init_volFrac,R,n);
  printf("L_sys = %lf\n",L_sys);


  generateBasis(basis_matrix,R,L_sys,2,N); // this takes a while for large vol frac
  init_volFrac = volumeFraction_in_drops(R,L_sys,n);
  printf("basis generated, phi_0 = %lf\n",init_volFrac);
  
  // Calculate coeff_matrix, and a copy of coeff_matrix. //
  set_matrices(coeff_matrix,copy,basis_matrix,R,L_sys,N);
  printf("matrices computed\n");
  current_t_val = tvals[0];

  make_chi_file(fname,path,prefix,&chi_file,Ravg,Rsigma,Req,
		dt_max,init_volFrac,alpha,beta,chi,L_sys,N,
		runNumber,randval);

  fprintf(chi_file,"%.12e\t%.12e\n",t,chi);

  clock_t end1 = clock();

  printf("time taken = %lf s\n",(double)(end1-begin1)/CLOCKS_PER_SEC);


  // Rs, matrix, and chi are all initialized. Time to run simulation! //

  clock_t begin2 = clock();

  printf("\n\n");

  double sum = 0;

  while (t < last_tval && n>1) {
    // Calculate the c vector. //
    set_GibbsThomson(c,R,chi,Req,n);
    update_matrix(copy,coeff_matrix,R,n,N);
    
    // solve for the B coefficients (B[n] is the supersaturation) //
    set_LU(&m,&c_vec,&B_vec,&p,c,B,coeff_matrix,s,n,N);
    gsl_linalg_LU_solve (&m.matrix, p, &c_vec.vector, &B_vec.vector);
    refine = gsl_matrix_view_array_with_tda(copy, (n), (n),N);
    work = gsl_vector_view_array(workingVector, (n));
    gsl_linalg_LU_refine(&refine.matrix,&m.matrix,p,&c_vec.vector,&B_vec.vector,&work.vector);


    if (isnan(dt)) printf("dt is NAN\n");

    dt = dtChoose(R,B,dt_max,&k,n);

    updateRadius_and_chi(R,R_last,B,dt,&chi,alpha,beta,L_sys,k,n);
    fprintf(chi_file,"%.12e\t%.12e\n",t,chi);

    // save initial R distribution
    if (t == 0) {
      save_Rdistribution(fname,path,prefix,R_files[JSET],runNumber,
			 JSET,R,R_last,dt,n,current_t_val,0);
      save_BASIS(fname,path,prefix,basis_files[JSET],runNumber,JSET,
		 basis_matrix,L_sys,n,current_t_val,0);
      JSET += 1;
      current_t_val = tvals[JSET];
    }
    t += dt;

    // next IF statement only executed if dt is such that a radius has shrunk. //
    if (k != 2*n) {
      rmv_zero(R,R_last,k,n);
      lastn = n;
      n = n-1;
      rearrange_matrix(coeff_matrix,copy,basis_matrix,k,lastn,N);
      sum = 0;
      for (i = 0; i < n; i++) {
	sum += B[i];
      }
    }

    if (t > current_t_val) {

      printf("t = %e, t%d = %e\n",t,JSET,tvals[JSET]);
      printf("(alpha(t)-alpha(0))/alpha(0)=%6.6e\n",change_in_alpha(totalVolume(R,n),chi,alpha,
							 beta,L_sys)/alpha);

      if (change_in_alpha(totalVolume(R,n),chi,alpha,beta,L_sys) != 0) {
	printf("not zero!\n");
      }

      if (!OVERLAP) {
	printf("checking for overlap...\n");
	OVERLAP = overlapping(basis_matrix,R,L_sys,n);
      }
      // save when t reaches current_t_val
      save_Rdistribution(fname,path,prefix,R_files[JSET],runNumber,
			 JSET,R,R_last,dt,n,current_t_val,OVERLAP);
      save_BASIS(fname,path,prefix,basis_files[JSET],runNumber,JSET,
		 basis_matrix,L_sys,n,current_t_val,OVERLAP);

      JSET += 1;
      if (JSET == T_INTERVALS) current_t_val = 1e30;
      else current_t_val = tvals[JSET];
    }

  }
  
  // save the final configuration //
  save_Rdistribution(fname,path,prefix,R_final_file,runNumber,
		     JSET,R,R_last,dt,n,t,OVERLAP);
  save_BASIS(fname,path,prefix,basis_final_file,runNumber,JSET,
	     basis_matrix,L_sys,n,t,OVERLAP);

  clock_t end2 = clock();
  
  printf("(alpha(t)-alpha(0))/alpha(0)=%6.6e\n",change_in_alpha(totalVolume(R,n),chi,alpha,
								beta,L_sys)/alpha);
  printf("time taken = %lf s\n",(double)(end2-begin2)/CLOCKS_PER_SEC);


  free(R);
  free(R_last);
  free(B);
  free(c);
  free(coeff_matrix);
  free(copy);
  free(basis_matrix);
  free(workingVector);

  free(prefix[0]);
  free(prefix);
  free(path[0]);
  free(path);
  free(fname[0]);
  free(fname[1]);
  free(fname);
  fclose(chi_file);
  
  gsl_permutation_free (p);
  gsl_rng_free(RNG);
  return 0;
}


void directory_build(int argnum,double *Ravg,double *Rsigma,
		     double *Req,double *dt_max,double *init_volFrac,
		     double *beta, double *chi,
		     int *N,char **prefix,char **path,
		     char *argv2,char *argv3)
{
  
  char *variable_is = "double";
  
  if (argnum <4) {
    printf("forgot final argument, exiting to system\n");
    exit(1);
  }
  else if (!strcmp(argv3,"test")) {
    variable_is = "blank";
  }
  else if (!strcmp(argv3,"varyRavg")) {
    *Ravg = atof(argv2);
    printf("%s %lf\n",argv3,*Ravg);
  }
  else if (!strcmp(argv3,"varyRsigma")) {
    *Rsigma = atof(argv2);
    printf("%s %lf\n",argv3,*Rsigma);
  }
  else if (!strcmp(argv3,"varyReq")) {
    *Req = atof(argv2);
    printf("%s %lf\n",argv3,*Req);
  }
  else if (!strcmp(argv3,"varydt_max")) {
    *dt_max = atof(argv2);
    printf("%s %lf\n",argv3,*dt_max);
  }
  else if (!strcmp(argv3,"varyvolFrac")) {
    *init_volFrac = atof(argv2);
    printf("%s %lf\n",argv3,*init_volFrac);
  }
  else if (!strcmp(argv3,"varybeta")) {
    *beta = atof(argv2);
    printf("%s %lf\n",argv3,*beta);
  }
  else if (!strcmp(argv3,"varychi")) {
    *chi = atof(argv2);
    printf("%s %lf\n",argv3,*chi);
  }
  else if (!strcmp(argv3,"varyN")) {
    *N = atoi(argv2);
    printf("%s %d\n",argv3,*N);
    variable_is = "int";
  }
  else {
    printf("bad final argument! Exiting to system\n");
    exit(1);
  }

  if (!strcmp(variable_is,"double")) {
    snprintf(prefix[0],NAME_LENGTH,"%.3e",atof(argv2));
    // path name for saved files //
  }
  else if (!strcmp(variable_is,"int")) {
    snprintf(prefix[0],NAME_LENGTH,"%d",atoi(argv2));
    // path name for saved files //
  }

  snprintf(path[0],NAME_LENGTH,"../%s/",argv3);
  return;
}

void save_BASIS(char **fname,char **path,char **prefix,
		FILE *basis_files_JSET,int runNumber,int JSET,
		double *basis_matrix,double L_sys, int n,
		double current_t_val,int OVERLAP)
{

  int i;
  
  int precision = 6;
  char xs[10] = "x";
  char ys[10] = "y";
  char zs[10] = "z";
  // open file for BASIS positions//    
  if (OVERLAP) {
    snprintf(fname[0],NAME_LENGTH,"%s%srun_%d_t%d_%.3e_BASIS_FILE_OVERLAP.txt",
	     path[0],prefix[0],runNumber,JSET,current_t_val);
  }
  else {
    snprintf(fname[0],NAME_LENGTH,"%s%srun_%d_t%d_%.3e_BASIS_FILE.txt",
	     path[0],prefix[0],runNumber,JSET,current_t_val);
  }
  basis_files_JSET = fopen(fname[0],"w");
  fprintf(basis_files_JSET,"# system size length scale L_sys = %6.6e\n",L_sys);
  fprintf(basis_files_JSET,"#%*s\t%*s\t%*s\n",precision+6,xs,precision+6,ys,
	  precision+6,zs);
  // save BASIS distribution //                                                      
  for (i = 0; i < n; i++) {
    fprintf(basis_files_JSET,"%.*e\t%.*e\t%.*e\n",precision,
	    L_sys*basis_matrix[i*3],precision,L_sys*basis_matrix[i*3+1],
	    precision,L_sys*basis_matrix[i*3+2]);
  }
  fclose(basis_files_JSET);

  return;
}

void save_Rdistribution(char **fname,char **path,char **prefix,
			FILE *R_files_JSET,int runNumber, int JSET,
			double *R,double *R_last,double dt,int n,
			double current_t_val,int OVERLAP)
{

  int i;

  if (OVERLAP) {
    printf("overlapping at t = %lf!\n",current_t_val);
    snprintf(fname[0],NAME_LENGTH,"%s%srun_%d_t%d_%.3eOVERLAP.txt",
	     path[0],prefix[0],runNumber,JSET,current_t_val);
  }
  else {
    snprintf(fname[0],NAME_LENGTH,"%s%srun_%d_t%d_%.3e.txt",
	     path[0],prefix[0],runNumber,JSET,current_t_val);
  }
  R_files_JSET = fopen(fname[0],"w");
  // save initial distribution at t = 0 //
  fprintf(R_files_JSET,"#dt=%.12e\n",dt);
  fprintf(R_files_JSET,"# tau_crslnk = inf\n");
  fprintf(R_files_JSET,"# R(t-dt)\tR(t)\n");
  for (i = 0; i < n; i++) fprintf(R_files_JSET,
				  "%.12e\t%.12e\n",
				  R_last[i],R[i]);
  fclose(R_files_JSET);
  return;
}

void make_chi_file(char **fname,char **path, char **prefix,
		   FILE **chi_file,double Ravg,double Rsigma,
		   double Req, double dt_max, double init_volFrac,
		   double alpha,double beta, double chi,double L_sys,
		   int N,int runNumber,unsigned long seed)
{
  snprintf(fname[1],NAME_LENGTH,"%s%srun_%d_chivst.txt",
	   path[0],prefix[0],runNumber);
  *chi_file = fopen(fname[1],"w");
  fprintf(*chi_file,"#RNG seed = %lu\n",
	  seed);
  fprintf(*chi_file,"#Ravg = %6.6e\n",Ravg);
  fprintf(*chi_file,"#Rsigma = %6.6e\n",Rsigma);
  fprintf(*chi_file,"#Req = %6.6e\n",Req);
  fprintf(*chi_file,"#dt_max = %6.6e\n",dt_max);
  fprintf(*chi_file,"#init_volFrac = %6.6e\n",init_volFrac);
  fprintf(*chi_file,"#alpha = %6.6e\n",alpha);
  fprintf(*chi_file,"#beta = %6.6e\n",beta);
  fprintf(*chi_file,"#chi(0) = %6.6e\n",chi);
  fprintf(*chi_file,"#L_sys = %6.6e\n",L_sys);
  fprintf(*chi_file,"#N(0) = %d\n",N);

  return;
}

void TVALS_2_save(double *tvals, double last_tval)
{
  int i;
  double spacing = log10(last_tval)/(T_INTERVALS-1);

  tvals[0] = 0;
  for (i = 1; i < T_INTERVALS; i++) {
    tvals[i] = pow(10,(i-1)*spacing);
  }
  return;
}

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

double set_L_sys(double init_volFrac, double *R, int n)
{
  int i;
  double sum = 0;

  for (i = 0; i < n; i++) sum += R[i]*R[i]*R[i];

  return cbrt(4*M_PI/3.0*sum/init_volFrac);
}

void set_LU(gsl_matrix_view *m, gsl_vector_view *c_vec, 
	    gsl_vector_view *B_vec, gsl_permutation **p,
	    double *c, double *B, double *coeff_matrix,int s, int n,int N)
{

  // Put B and c vectors into appropriate form for gsl library. //
  *B_vec = gsl_vector_view_array(B, (n));                         // B is TBD by solver
  *c_vec = gsl_vector_view_array(c,(n));

  // Set up and perform LU decomposition. //
  *m = gsl_matrix_view_array_with_tda(coeff_matrix, (n), (n),N);
  *p = gsl_permutation_realloc(*p,(n));
  gsl_linalg_LU_decomp(&(*m).matrix, *p, &s);

  return;
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


void rearrange_matrix(double *coeff_matrix, double *copy, double *basis_matrix,
		      int k, int lastn,int N)
{
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

void set_matrices(double *coeff_matrix,double *copy,double *basis_matrix,
		  double *R,double L_sys, int N)
{
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

void set_GibbsThomson(double *c, double *R, double chi,
		      double Req, int n)
{
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


void onevalue_distribution(double *R,double Ravg,int N)
{
  int i;

  for (i = 0; i < N; i++) R[i] = Ravg;

  return;
}

void gaussianDistribution(double *R, double Ravg, double Rsigma,int N)
{
  int i;

  for (i = 0; i < N; i++) {
    R[i] = gsl_ran_gaussian(RNG,Rsigma)+Ravg;
    while (R[i] <= 0) R[i] = gsl_ran_gaussian(RNG,Rsigma)+Ravg;
  }

  return;
}

void updateRadius_and_chi(double *R, double *R_last, double *B,
			  double dt,double *chi,double alpha,
			  double beta,double L_sys,int k,int n)
{
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

double offDiagonal(double L_sys,double *basis_matrix,int i, int j,int N)
{
  return -1.0/(L_sys*distanceBasis(basis_matrix,i,j,0,0,0));
}

double G(double x)
{

  if (x > 1.5) return 0;

  double c = 3*(3*3*3)*2.718281828/3.174802104;

  return c*x*x*exp(1.5/(x-1.5))/(pow(1.5-x,11.0/3.0)*pow(x+3,7.0/3.0));
}

double maximum_ofG(double x0,double xavg,int Npoints)
{
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

void nucleatedDistribution(double *R,double Ravg,int N)
{
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
	R[i] = z*Ravg;
	Found = 1;
      }
    }
  }
  return;
}

double distanceBasis(double *basis_matrix,int i1, int i2, int u,
		     int v, int w)
{
  double delx = basis_matrix[i1*3]+u-basis_matrix[i2*3];
  double dely = basis_matrix[i1*3+1]+v-basis_matrix[i2*3+1];
  double delz = basis_matrix[i1*3+2]+w-basis_matrix[i2*3+2];

  return sqrt(delx*delx+dely*dely+delz*delz);
}

int checkOverlap(double *basis_matrix,int i1, int i2,double *R,double L_sys,double buffer) 
{
  if (distanceBasis(basis_matrix,i1,i2,0,0,0)*L_sys/buffer<=R[i1]+R[i2]) return 1;
  else if (distanceBasis(basis_matrix,i1,i2,1,0,0)*L_sys/buffer<=R[i1]+R[i2]) return 1;
  else if (distanceBasis(basis_matrix,i1,i2,1,1,0)*L_sys/buffer<=R[i1]+R[i2]) return 1;
  else if (distanceBasis(basis_matrix,i1,i2,0,1,0)*L_sys/buffer<=R[i1]+R[i2]) return 1;
  else if (distanceBasis(basis_matrix,i1,i2,0,0,1)*L_sys/buffer<=R[i1]+R[i2]) return 1;
  else if (distanceBasis(basis_matrix,i1,i2,1,0,1)*L_sys/buffer<=R[i1]+R[i2]) return 1;
  else if (distanceBasis(basis_matrix,i1,i2,0,1,1)*L_sys/buffer<=R[i1]+R[i2]) return 1;
  else if (distanceBasis(basis_matrix,i1,i2,1,1,1)*L_sys/buffer<=R[i1]+R[i2]) return 1;
  else if (distanceBasis(basis_matrix,i2,i1,1,0,0)*L_sys/buffer<=R[i1]+R[i2]) return 1;
  else if (distanceBasis(basis_matrix,i2,i1,1,1,0)*L_sys/buffer<=R[i1]+R[i2]) return 1;
  else if (distanceBasis(basis_matrix,i2,i1,0,1,0)*L_sys/buffer<=R[i1]+R[i2]) return 1;
  else if (distanceBasis(basis_matrix,i2,i1,0,0,1)*L_sys/buffer<=R[i1]+R[i2]) return 1;
  else if (distanceBasis(basis_matrix,i2,i1,1,0,1)*L_sys/buffer<=R[i1]+R[i2]) return 1;
  else if (distanceBasis(basis_matrix,i2,i1,0,1,1)*L_sys/buffer<=R[i1]+R[i2]) return 1;
  else if (distanceBasis(basis_matrix,i2,i1,1,1,1)*L_sys/buffer<=R[i1]+R[i2]) return 1;
  else return 0;
}



void generateBasis(double *basis_matrix,double *R,double L_sys, double buffer, int N)
{
  int i,j;

  for (i = 0; i < N; i++) {
    basis_matrix[i*3] = gsl_rng_uniform(RNG);
    basis_matrix[i*3+1] = gsl_rng_uniform(RNG);
    basis_matrix[i*3+2] = gsl_rng_uniform(RNG);
    
    for (j = 0; j < i; j++) {
      if (checkOverlap(basis_matrix,i,j,R,L_sys,buffer)) {
	i -= 1;
	break;
      }
    }
  }
  return;
}


int overlapping(double *basis_matrix,double *R, double L_sys, int n)
{
  int i,j;
  int confirm_overlap = 0;
  for (i = 0; i < n; i++) {
    for (j = 0; j < i; j++) {
      if (checkOverlap(basis_matrix,i,j,R,L_sys,1)) {
	printf("overlapping particles at (%d,%d)!!\n",i,j);
	confirm_overlap = 1;
      }
    }
  }
  if (confirm_overlap) return 1;
  return 0;
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
