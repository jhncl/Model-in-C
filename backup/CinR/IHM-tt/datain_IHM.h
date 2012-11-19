#ifndef _datain_h
#define _datain_h

int testargc_IHM(int argc);

int datadouble_IHM(char filename[],char filename2[],struct_data_IHM *D,double *QFADyA,double *QFADyB);
int dataint_IHM(char filename[],char filename2[], int datavec[] ,int length);
int dataLMN_IHM(char filename[],char filename2[], int *datavecL);

int inzstruct_MH_IHM(struct_MH_IHM *MH);
int inzstruct_priors_IHM(struct_priors_IHM *D_priors,double *PRIORS);
int inzstruct_data_IHM(struct_data_IHM *data,int *QFAIA,double *QFADyA,int *QFADNoORFA,int *QFAIB,double *QFADyB,int *QFADNoORFB);
int inzstruct_para_IHM(struct_para_IHM *para,struct_data_IHM *data,struct_priors_IHM *priors);

int fillMH_IHM(struct_MH_IHM *MH);
int filldata_IHM(struct_data_IHM *D);
int fillpara_IHM(struct_para_IHM *D_para,struct_data_IHM *D,struct_priors_IHM *D_priors);
int fillpriors_IHM(struct_priors_IHM *D_priors);

#endif
