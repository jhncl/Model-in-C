#include "headers_IHM.h"
#include "datain_IHM.h"
#include "functions_IHM.h"
#include "print_IHM.h"
     int
    main_IHM(int *arga,int *argb,int *argc,double *OUT,char **HEADER,int *QFAIA,double *QFADyA,int *QFADNoORFA,int *QFAIB,double *QFADyB,int *QFADNoORFB,double *PRIORS)/****/
     {
	struct_data_IHM *data= malloc(sizeof(struct_data_IHM));
	struct_para_IHM *para= malloc(sizeof(struct_para_IHM));
	struct_priors_IHM *priors= malloc(sizeof(struct_priors_IHM));
	struct_MH_IHM *MH = malloc(sizeof(struct_MH_IHM));

	int burn,iters,thin;
	burn=*arga;   /*Burn in*/
	iters=*argb;    /*iterations*/
	thin=*argc;        /*thining*/
       inzstruct_data_IHM(data,QFAIA,QFADyA,QFADNoORFA,QFAIB,QFADyB,QFADNoORFB);
	inzstruct_priors_IHM(priors,PRIORS);
	inzstruct_para_IHM(para,data,priors);
	inzstruct_MH_IHM(MH);

	gibbsandMHloop_IHM(burn,1,data,para,priors,MH,0,OUT,HEADER);
	gibbsandMHloop_IHM(iters,thin,data,para,priors,MH,1,OUT,HEADER);


	return 0;
}
