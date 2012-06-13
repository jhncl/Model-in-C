#ifndef _functions_h
#define _functions_h

/* public common function prototypes */
double logistic_function(double t,double K, double r, double P);
double logistic_function_E(double t,double K, double r, double P);
double gauss_sample(gsl_rng *RNG, struct_data *D ,int start, int N,double x[],double tau,double mu_0,double tau_0);
double gamma_sample(gsl_rng *RNG, struct_data *D ,int start, int N, double x[], double mu,double alpha,double beta);
double MCMC_base(gsl_rng *RNG, struct_data *D,struct_para *D_para,struct_priors *D_priors,double *accept,double *h,double para,double (*foo)(struct struct_data *D,struct struct_para *D_para,struct struct_priors *D_priors,double,int,int),int l, int m);
double MCMC_P(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int l, int m);
double MCMC_nu_l(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int l, int m);
double MCMC_K_lm(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int l, int m);
double MCMC_r_lm(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int l, int m);
int gibbsandMHloop(int iter,int thin,gsl_rng *RNG,struct_data *D,struct_para *D_para,struct_priors *D_priors ,struct_MH *D_MH,int CAPL,int print,double *OUT,char **HEADER);/*********************/

double MCMC_sigma_K_o(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int l, int m);
double MCMC_sigma_r_o(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int l, int m);
double MCMC_tau_K_l(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int l, int m);
double MCMC_tau_r_l(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int l, int m);
double MCMC_sigma_nu(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int l, int m);
#endif
