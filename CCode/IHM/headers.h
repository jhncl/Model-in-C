#include <stdio.h>
#include <math.h>
#include <gsl/gsl_randist.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h> 
#include <gsl/gsl_sf_log.h>      
#include <unistd.h>
#include <time.h>
#include <gsl/gsl_matrix.h>
#include <string.h>

typedef struct struct_data {
  double *y;
  int L, M,N,*NoORF,*NoSUM,*NoTIME,SHIFTmn,MAXmn,CAPiter;
} struct_data;

typedef struct struct_MH {
  double hZ,hnu,hup,accept_Z,accept_nu,accept_up;
} struct_MH;

typedef struct struct_para {
  double
    *alpha_c,
    *delta_l,
    *gamma_cl,	        sigma_gamma,
    *upsilon_c,       sigma_upsilon,
    upsilon_p,

    *Z_l,              sigma_Z,         
    *nu_l,             sigma_nu,
    Z_p,
    nu_p;

} struct_para;

typedef struct struct_priors {
  double
	Z_mu,			eta_Z_p,
	eta_Z,			psi_Z,

	eta_nu,			psi_nu,
	nu_mu,			eta_nu_p,
	alpha_mu,		eta_alpha,
	p,   
	eta_gamma,		psi_gamma,
	eta_upsilon,		psi_upsilon,	    
	upsilon_mu,		eta_upsilon_p;
} struct_priors;


