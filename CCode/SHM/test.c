#include "test.h"

     int
     main (void)
     {
/*struct_data data;*/
int y[386460];
int x[386460];
double K_lm[187000];
double tau_K_l[4294]; 
double r_lm[187000];
double tau_r_l[4294];
double K_o_l[4294];
double r_o_l[4294];		
double nu_l[4294];		
	const gsl_rng_type * T;
	
double test[100000];
	gsl_rng * RNG;
	gsl_rng_env_setup ();
	T = gsl_rng_default;
	RNG = gsl_rng_alloc (T);
     
 	gsl_rng_set (RNG, 1); /*seed*/

test[1]=1;
y[286460]=1;
x[386460]=1;
K_lm[187000]=1;
tau_K_l[4294]=1; 
r_lm[187000]=1;
tau_r_l[4294]=1;
K_o_l[4294]=1;
r_o_l[4294]=1;		
nu_l[4294]=1;
       	gsl_rng_free(RNG);
	return 0;
}
