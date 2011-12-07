#include <stdio.h>
#include <math.h>
#include <gsl/gsl_randist.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h> 
#include <gsl/gsl_sf_log.h>      
#include <unistd.h>
#include <time.h>

typedef struct struct_data {
  double *y, *x;
  int L, M,N,*NoORF,*NoSUM,*NoTIME,SHIFTlmn,SHIFTmn, maxy,maxTIMEa,maxTIMEb;
} struct_data;


typedef struct struct_para {
  double
    *alpha,
    *beta,
    *delta,
    *gamma,	        sigma_gamma,
    *omega,            sigma_omega,
    *upsilon_c,       sigma_upsilon,
    upsilon_p,
 
    *K_clm,            *tau_K_cl,
    *r_clm,            *tau_r_cl,

    *K_o_l,            sigma_K_o,
    *r_o_l,            sigma_r_o,
    *nu_l,             sigma_nu,

    K_p,
    r_p,
    nu_p,
    P;
} struct_para;


typedef struct struct_priors {
  double
    alpha_mu,               eta_alpha,
    beta_mu,                eta_beta,
    p,    
    eta_gamma,              psi_gamma,
    eta_omega,              psi_omega,
    eta_upsilon,	    phi_upsilon,	    
    upsilon_mu,            eta_upsilon_mu,

    sigma_K,                phi_K,
    sigma_r,                phi_r,

    eta_K_p,                psi_K_o,
    eta_r_p,                psi_r_o,
    eta_nu_p,               psi_nu,

    K_mu,                   eta_K_mu,
    r_mu,                   eta_r_mu,
    nu_mu,                  eta_nu_mu,
    P_mu,                   eta_P_mu;
} struct_priors;

typedef struct struct_MH {
  double hP,hnu,hK,hr,accept_P,accept_nu,accept_K,accept_r;
} struct_MH;

