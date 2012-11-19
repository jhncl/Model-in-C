#ifndef _datain_h
#define _datain_h

int testargc(int argc);

int datadouble(char filename[],char filename2[],struct_data *D);
int dataint(char filename[],char filename2[], int datavec[] ,int length);
int dataLMN(char filename[],char filename2[], int *datavecL);

int inzstruct_MH(struct_MH *MH);
int inzstruct_priors(struct_priors *priors);
int inzstruct_data(struct_data *data);
int inzstruct_para(struct_para *para,struct_data *data,struct_priors *priors);

int fillMH(struct_MH *MH);
int filldata(struct_data *D);
int fillpara(struct_para *D_para,struct_data *D,struct_priors *D_priors);
int fillpriors(struct_priors *D_priors);

#endif
