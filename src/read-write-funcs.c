#include <stdio.h>
#include <stdlib.h>
#include "headerfile.h"

void set_params(char **args,struct params *p)
{

  void print_params(struct params *p);

  sscanf(args[2],"%lf",&p->R_avg0);
  sscanf(args[3],"%lf",&p->sigma_Ravg0);
  sscanf(args[4],"%lf",&p->R_eq);
  sscanf(args[5],"%lf",&p->volFrac_0);
  sscanf(args[6],"%lf",&p->beta);
  sscanf(args[7],"%lf",&p->chi_0);
  sscanf(args[8],"%d",&p->N);
  sscanf(args[9],"%lf",&p->dt_max);
  sscanf(args[10],"%lf",&p->t_final);

  print_params(p);

  return;
}


void print_params(struct params *p)
{
  
  printf("R_avg0=%e\n",p->R_avg0);
  printf("sigma_Ravg0=%e\n",p->sigma_Ravg0);
  printf("R_eq=%e\n",p->R_eq);
  printf("volFrac_0=%e\n",p->volFrac_0);
  printf("beta=%e\n",p->beta);
  printf("chi_0=%e\n",p->chi_0);
  printf("N=%e\n",p->N);
  printf("dt_max=%e\n",p->dt_max);

  return;
}



void initialize_file(FILE **output,char *path,char *fname,struct params p)
{

  int num_chars = 400;
  char suffix[num_chars];
  char f[num_chars];

  snprintf(suffix,sizeof(suffix),"%1.4e_%1.4e_%1.4e_%1.4e_%1.4e"
	   "%1.4e.txt",p.R_avg0,p.sigma_Ravg0,p.R_eq,
	   p.volFrac_0,p.beta,p.chi_0);
  snprintf(f,sizeof(f),"%s%s_%s",path,fname,suffix);

  *output = fopen(f,"w");

  return;

}
