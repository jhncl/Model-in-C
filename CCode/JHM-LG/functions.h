#ifndef _functions_h
#define _functions_h

/* public common function prototypes */
double logistic_function(double t,double K, double r, double P);
double logistic_function_E(double t,double K, double r, double P);

double gauss_sample(gsl_rng *RNG, struct_data *D ,int start, int N,double x[],double tau,double mu_0,double tau_0);

double MCMC_base(gsl_rng *RNG, struct_data *D,struct_para *D_para,struct_priors *D_priors,double *accept,double *h,double para,double (*foo)(struct struct_data *D,struct struct_para *D_para,struct struct_priors *D_priors,double,int,int,int),int c,int l, int m);

double MCMC_P(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int c,int l, int m);
double MCMC_nu_l(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int c,int l, int m);
double MCMC_K_clm(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int c,int l, int m);
double MCMC_r_clm(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int c,int l, int m);
double MCMC_K_o_l(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int c,int l, int m);
double MCMC_r_o_l(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int c,int l, int m);
double MCMC_alpha_c(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int c,int l, int m);
double MCMC_beta_c(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int c,int l, int m);
double aug_delta_l(gsl_rng *RNG,struct_data *D,struct_para *D_para,struct_priors *D_priors,int l);
double MCMC_gamma_cl(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int c,int l, int m);
double MCMC_omega_cl(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int c,int l, int m);
double MCMC_sigma_gamma(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int c,int l, int m);
double MCMC_sigma_omega(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int c,int l, int m);
double MCMC_upsilon_c(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int c,int l, int m);
double MCMC_sigma_upsilon(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int c,int l, int m);
double MCMC_sigma_K_o(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int c,int l, int m);
double MCMC_sigma_r_o(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int c,int l, int m);
double MCMC_tau_K_cl(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int c,int l, int m);
double MCMC_tau_r_cl(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int c,int l, int m);

int gibbsandMHloop(int iter,int thin,gsl_rng *RNG,struct_data *D,struct_para *D_para,struct_priors *D_priors,struct_MH *D_MH, int print);

void update_lik_a(double para_mean[],double para_prec[],double scale[],double interaction[],struct_data *D,struct_para *D_para,double *mu,double *sig);
void update_lik_b(double pass_mean[],double para_var,double pass_var[],struct_data *D,struct_para *D_para,double *mu,double *sig);
double sim_para(gsl_rng *RNG,double pass_mean,double pass_var,double para_mean,double para_var);

#endif
