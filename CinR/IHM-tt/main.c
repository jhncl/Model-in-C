#include "headers.h"
#include "datain.h"
#include "functions.h"
#include "print.h"
     int
     main (int *arga,int *argb,int *argc,double *OUT,char **HEADER,int *QFAIA,double *QFADyA,int *QFADNoORFA,int *QFAIB,double *QFADyB,int *QFADNoORFB,double *PRIORS)/****/
     {
	struct_data *data= malloc(sizeof(struct_data));
	struct_para *para= malloc(sizeof(struct_para));
	struct_priors *priors= malloc(sizeof(struct_priors));
	struct_MH *MH = malloc(sizeof(struct_MH));


	int burn,iters,thin;
	long seed;
	const gsl_rng_type * T;
	gsl_rng * RNG;

	/*testargc(argc);*/

	gsl_rng_env_setup ();
	T = gsl_rng_default;
	RNG = gsl_rng_alloc (T);
	seed = time (NULL) * getpid();    
  	gsl_rng_set (RNG, seed); /*seed*/
	burn=*arga;   /*Burn in*/
	iters=*argb;    /*iterations*/
	thin=*argc;        /*thining*/
       inzstruct_data(data,QFAIA,QFADyA,QFADNoORFA,QFAIB,QFADyB,QFADNoORFB);
	inzstruct_priors(priors,PRIORS);
	inzstruct_para(para,data,priors);
	inzstruct_MH(MH);

	gibbsandMHloop(burn,1,RNG,data,para,priors,MH,0,OUT,HEADER);
	gibbsandMHloop(iters,thin,RNG,data,para,priors,MH,1,OUT,HEADER);

       	gsl_rng_free(RNG);
	return 0;
}
