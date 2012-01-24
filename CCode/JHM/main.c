#include "headers.h"
#include "datain.h"
#include "functions.h"
#include "print.h"
int
main (int argc,char *argv[])
{
	struct_data *data= malloc(sizeof(struct_data));
	struct_para *para= malloc(sizeof(struct_para));
	struct_priors *priors= malloc(sizeof(struct_priors));
	struct_MH *MH = malloc(sizeof(struct_MH));

	int burn,iters,thin;
	long seed;
	const gsl_rng_type * T;
	gsl_rng * RNG;

	testargc(argc);

	gsl_rng_env_setup ();
	T = gsl_rng_default;
	RNG = gsl_rng_alloc (T);
	seed = time (NULL) * getpid();    
  	gsl_rng_set (RNG, seed); /*seed*/

	burn=atoi(argv[1]);   /*Burn in*/
	iters=atoi(argv[2]);    /*iterations*/
	thin=atoi(argv[3]);        /*thining*/

        inzstruct_data(data);
      	inzstruct_para(para,data);
	inzstruct_priors(priors);
	inzstruct_MH(MH);
/*************************/
para->K_o_l[0]=-1.28022656;
para->r_o_l[0]=-0.42613390;
para->K_clm[0]=-3.58148222;

para->P=-8.59382058;
para->upsilon_c[1]=-1;
para->nu_l[0]=13.29172915;
para->tau_K_cl[0]=-2.17165073;
para->tau_r_cl[0]=3.97912256;
para->alpha_c[1]=-1.84746155;    
para->beta_c[1]=0.42673077;
para->r_clm[0]=0.52425651;
para->r_clm[1]=0.63300187;
para->r_clm[2]=0.70384475;
para->r_clm[3]=0.64524218 ;
para->r_clm[4]=0.69444077 ;
para->r_clm[5]=0.85463257;   

para->K_clm[0]=-3.58148222 ;
para->K_clm[1]=-2.90553615;
para->K_clm[2]=-2.65059866;
para->K_clm[3]=-2.92556808 ;
para->K_clm[4]=-2.92762388  ;
para->K_clm[5]=-2.23997120; 
       
para->r_clm[24]=1.00955579 ;
para->r_clm[25]=1.02423878;
para->r_clm[26]= 1.00212593;
para->r_clm[27]=0.94895349   ;
para->r_clm[28]=0.94584712  ;
para->r_clm[29]=1.05514850 ;   
para->K_clm[24]=-1.98346830 ;
para->K_clm[25]= -1.85942304 ;
para->K_clm[26]=   -2.06285725;
para->K_clm[27]=-2.02502261  ;
para->K_clm[28]= -2.19092948  ;
para->K_clm[29]=-1.96633827;  


/*************************/
	gibbsandMHloop(burn,1,RNG,data,para,priors,MH,0);
	gibbsandMHloop(iters,thin,RNG,data,para,priors,MH,1);

       	gsl_rng_free(RNG);
	return 0;
}
