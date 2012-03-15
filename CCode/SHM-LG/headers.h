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
  int L, M,N,*NoORF,*NoSUM,*NoTIME,SHIFTlmn,maxy,maxNoTIME;
} struct_data;

typedef struct struct_MH {
  double hP,hnu,hK,hr,halpha,accept_P,accept_nu,accept_K,accept_r;
} struct_MH;

typedef struct struct_para {
  double
    *K_lm,            *tau_K_l,
    *r_lm,            *tau_r_l,

    *K_o_l,            sigma_K_o,
    *r_o_l,            sigma_r_o,
    *nu_l,             sigma_nu,

    K_p,
    r_p,
    nu_p,
    P,
	*star_K_lm, *star_tau_K_l,
	star_K_l, star_sigma_K_o,
	*star_r_lm, *star_tau_r_l,
	star_r_l, star_sigma_r_o;
} struct_para;

typedef struct struct_priors {
  double
    sigma_K,                phi_K,
    sigma_r,                phi_r,

    eta_K_o,                psi_K_o,
    eta_r_o,                psi_r_o,
    eta_nu,                 psi_nu,

    K_mu,                   eta_K_p,
    r_mu,                   eta_r_p,
    nu_mu,                  eta_nu_p,
    P_mu,                   eta_P;
} struct_priors;


