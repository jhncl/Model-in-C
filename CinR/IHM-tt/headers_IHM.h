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
#include <R.h>
#include <Rmath.h>

typedef struct struct_data_IHM{
  double *y;
  int L, M,N,*NoORF,*NoSUM,*NoTIME,SHIFTmn,MAXmn;
} struct_data_IHM;

typedef struct struct_MH_IHM{
  double 
halpha_c,
hsigma_gamma,
hupsilon_c,
hsigma_upsilon,
hsigma_nu,
hsigma_Z,
hnu_p,
hgamma_cl,
hZ_l,
hnu_l,
accept_Z,accept_nu,accept_up;
} struct_MH_IHM;

typedef struct struct_para_IHM{
  double
    *alpha_c,
    *delta_l,
    *gamma_cl,	        sigma_gamma,
    *upsilon_c,       sigma_upsilon,

    *Z_l,              sigma_Z,         
    *nu_l,             sigma_nu,
    Z_p,
    nu_p;

} struct_para_IHM;

typedef struct struct_priors_IHM{
  double
	Z_mu,			eta_Z_p,
	eta_Z,			psi_Z,

	eta_nu,			psi_nu,
	nu_mu,			eta_nu_p,
	alpha_mu,		eta_alpha,
	p,   
	eta_gamma,		psi_gamma,
	eta_upsilon,		psi_upsilon,	    
    upsilon_mu,
    df;
} struct_priors_IHM;


