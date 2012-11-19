#ifndef _datain_h
#define _datain_h

int testargc(int argc);
int testsame(int a,int b);

int datadouble(char filename[], double datavec[], int length);
int dataint(char filename[], int datavec[], int length);
int dataLMN(char filename[], int *datavecL, int *datavecM,int *datavecN,int *datavecmaxy,int *datavecmaxTIME);

int inzstruct_MH(struct_MH *MH);
int inzstruct_priors(struct_priors *D_priors,double *PRIORS);
int inzstruct_data(struct_data *data,int *QFAI,double *QFADy,double *QFADx,int *QFADNoORF,int *QFADNoTIME);
int inzstruct_para(struct_para *para,struct_data *data,struct_priors *D_priors);

int fillMH(struct_MH *MH);
int filldata(struct_data *D);
int fillpara(struct_para *D_para,struct_data *D,struct_priors *D_priors);
int fillpriors(struct_priors *D_priors);

#endif
