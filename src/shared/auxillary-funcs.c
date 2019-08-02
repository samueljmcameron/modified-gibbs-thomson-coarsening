#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../headerfile.h"



void save_basis(char *path,int jset,
		double *basis_matrix,double L_sys, int n,bool overlap,
		struct params p)
{

  void initialize_file(FILE **output,char *path,char *fname,struct params p,
		       bool overlap);

  FILE *basis_file;
  
  int precision = 6;
  char xs[10] = "x";
  char ys[10] = "y";
  char zs[10] = "z";

  char fname[100];


  // open file for BASIS positions//    

  snprintf(fname,sizeof(fname),"basis_%d",jset);
  initialize_file(&basis_file,path,fname,p,overlap);
  fprintf(basis_file,"# system size length scale L_sys = %6.6e\n",L_sys);
  fprintf(basis_file,"#%*s\t%*s\t%*s\n",precision+6,xs,precision+6,ys,
	  precision+6,zs);
  // save BASIS distribution //                                                      
  for (int i = 0; i < n; i++) {
    fprintf(basis_file,"%.*e\t%.*e\t%.*e\n",precision,
	    L_sys*basis_matrix[i*3],precision,L_sys*basis_matrix[i*3+1],
	    precision,L_sys*basis_matrix[i*3+2]);
  }
  fclose(basis_file);

  return;
}

void save_Rdistribution(char *path,int jset,double *R,double *R_last,
			double dt,int n,bool overlap,struct params p)
{

  void initialize_file(FILE **output,char *path,char *fname,struct params p,
		       bool overlap);

  FILE *R_file;

  int i;

  char fname[100];

  snprintf(fname,sizeof(fname),"radius_%d",jset);

  initialize_file(&R_file,path,fname,p,overlap);

  fprintf(R_file,"#dt=%.12e\n",dt);
  fprintf(R_file,"# tau_crslnk = inf\n");
  fprintf(R_file,"# R(t-dt)\tR(t)\n");
  for (i = 0; i < n; i++) fprintf(R_file,
				  "%.12e\t%.12e\n",
				  R_last[i],R[i]);
  fclose(R_file);
  return;
}


void build_tvals(double *tvals, double last_tval,int t_intervals)
{

  double spacing = log10(last_tval)/(t_intervals-1);

  tvals[0] = 0;
  for (int i = 1; i < t_intervals; i++) {
    tvals[i] = pow(10,(i-1)*spacing);
  }
  return;
}

unsigned long int generate_seed(void)
{

  unsigned long int randval;
  FILE *seedFile;
  seedFile = fopen("/dev/urandom","r");
  
  int CHECK = fread(&randval,sizeof(randval),1,seedFile);
  
  fclose(seedFile);


  return randval;

}

void free_vectors(double **R,double **R_last,double **c,double **B,
		  double **workingVector,double **coeff_matrix,double **copy,
		  double **basis_matrix,double **tvals)
{

  free(*R);
  free(*R_last);
  free(*B);
  free(*c);
  free(*coeff_matrix);
  free(*copy);
  free(*basis_matrix);
  free(*workingVector);

  return;
}

void allocate_vectors(int N,int t_intervals,double **R,double **R_last,double **c,
		      double **B,double **workingVector,double **coeff_matrix,
		      double **copy,double **basis_matrix,double **tvals)
{

  *R = malloc(N*sizeof(double));
  *R_last = malloc(N*sizeof(double));
  *c = malloc((N)*sizeof(double));
  *B = malloc((N)*sizeof(double));
  *workingVector = malloc((N)*sizeof(double));
  *coeff_matrix = malloc(sizeof(double)*(N)*(N));
  *copy = malloc(sizeof(double)*(N)*(N));
  *basis_matrix = malloc(sizeof(double)*N*3);
  *tvals = malloc(t_intervals*sizeof(double));

  return;

}

void set_params(char **args,struct params *p)
{

  void print_params(struct params *p);

  sscanf(args[2],"%lf",&p->R_avg0);
  sscanf(args[3],"%lf",&p->sigma_R0);
  sscanf(args[4],"%lf",&p->R_eq);
  sscanf(args[5],"%lf",&p->volFrac_0);
  sscanf(args[6],"%lf",&p->beta);
  sscanf(args[7],"%lf",&p->chi_0);
  sscanf(args[8],"%d",&p->N);
  sscanf(args[9],"%lf",&p->dt_max);
  sscanf(args[10],"%lf",&p->t_final);
  sscanf(args[11],"%d",&p->t_intervals);

  print_params(p);

  return;
}


void print_params(struct params *p)
{
  
  printf("R_avg0=%e\n",p->R_avg0);
  printf("sigma_R0=%e\n",p->sigma_R0);
  printf("R_eq=%e\n",p->R_eq);
  printf("volFrac_0=%e\n",p->volFrac_0);
  printf("beta=%e\n",p->beta);
  printf("chi_0=%e\n",p->chi_0);
  printf("N=%d\n",p->N);
  printf("dt_max=%e\n",p->dt_max);

  return;
}



void initialize_file(FILE **output,char *path,char *fname,struct params p,
		     bool overlap)
{

  int num_chars = 400;
  char suffix[num_chars];
  char f[num_chars];

  if (overlap) {

    snprintf(suffix,sizeof(suffix),"%1.4e_%1.4e_%1.4e_%1.4e_%1.4e_"
	     "%1.4e-overlap.txt",p.R_avg0,p.sigma_R0,p.R_eq,
	     p.volFrac_0,p.beta,p.chi_0);
  } else {

    snprintf(suffix,sizeof(suffix),"%1.4e_%1.4e_%1.4e_%1.4e_%1.4e_"
	     "%1.4e.txt",p.R_avg0,p.sigma_R0,p.R_eq,
	     p.volFrac_0,p.beta,p.chi_0);

  }

  snprintf(f,sizeof(f),"%s%s_%s",path,fname,suffix);

  *output = fopen(f,"w");

  return;

}
