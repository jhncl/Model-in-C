#include "headers.h"
#include "datain.h"
#include "functions.h"
#include "print.h"
     int
     main (int *arga,int *argb,int *argc,int *argd,double *OUT, char **HEADER,int *QFAI,double *QFADy,double *QFADx,int *QFADNoORF,int *QFADNoTIME,double *PRIORS)
     {
	struct_data *data= malloc(sizeof(struct_data));
	struct_para *para= malloc(sizeof(struct_para));
	struct_priors *priors= malloc(sizeof(struct_priors));
	struct_MH *MH = malloc(sizeof(struct_MH));

	int burn,iters,thin, CAPL;

	burn=*arga;    /*Burn in*/
	iters=*argb;    /*iterations*/
	thin=*argc;         /*thining*/

	CAPL=*argd;              /*CAP D->L*/

        inzstruct_data(data,QFAI,QFADy,QFADx,QFADNoORF,QFADNoTIME);
	inzstruct_priors(priors,PRIORS);
	inzstruct_para(para,data,priors);
	inzstruct_MH(MH);
	gibbsandMHloop(burn,1,data,para,priors,MH,CAPL,0,OUT,HEADER);
	gibbsandMHloop(iters,thin,data,para,priors,MH,CAPL,1,OUT,HEADER);

	return 0;
}
