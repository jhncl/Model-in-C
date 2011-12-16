#ifndef _functions_h
#define _functions_h

/* public common function prototypes */
double gauss_sample(gsl_rng *RNG, struct_data *D ,int start, int N,double x[],double tau,double mu_0,double tau_0);
double MCMC_base(gsl_rng *RNG, struct_data *D,struct_para *D_para,struct_priors *D_priors,double *accept,double *h,double para,double (*foo)(struct struct_data *D,struct struct_para *D_para,struct struct_priors *D_priors,double,int,int,int),int c,int l, int m);
double MCMC_alpha_c(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int c,int l, int m);
int gibbsandMHloop(int iter,int thin,gsl_rng *RNG,struct_data *D,struct_para *D_para,struct_priors *D_priors ,struct_MH *D_MH,int print);
#endif
