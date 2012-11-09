#include "headers_JHM.h"
#include "datain_JHM.h"
#include "functions_JHM.h"
#include "print_JHM.h"
int
     main_JHM (int *arga,int *argb,int *argc,double *OUT, char **HEADER,int *QFAIA,double *QFADyA,double *QFADxA,int *QFADNoORFA,int *QFANoTIMEA,int *QFAIB,double *QFADyB,double *QFADxB,int *QFADNoORFB,int *QFANoTIMEB,double *PRIORS)/****/
{
  GetRNGstate();
	struct_data_JHM *data= malloc(sizeof(struct_data_JHM));
	struct_para_JHM *para= malloc(sizeof(struct_para_JHM));
	struct_priors_JHM *priors= malloc(sizeof(struct_priors_JHM));
	struct_MH_JHM *MH = malloc(sizeof(struct_MH_JHM));

	int burn,iters,thin;

	burn=*arga;   /*Burn in*/
	iters=*argb;    /*iterations*/
	thin=*argc;        /*thining*/

        inzstruct_data_JHM(data,QFAIA,QFADyA,QFADxA,QFADNoORFA,QFANoTIMEA,QFAIB,QFADyB,QFADxB,QFADNoORFB,QFANoTIMEB);
	inzstruct_priors_JHM(priors,PRIORS);
      	inzstruct_para_JHM(para,data,priors);
	inzstruct_MH_JHM(MH);

	gibbsandMHloop_JHM(burn,1,data,para,priors,MH,0,OUT,HEADER);
	gibbsandMHloop_JHM(iters,thin,data,para,priors,MH,1,OUT,HEADER);
  PutRNGstate();
	return 0;
}
