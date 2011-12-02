#ifndef _datain_h
#define _datain_h

int testargc(int argc);

int datadouble(char filename[], double datavec[] );
int dataint(char filename[], int datavec[] );
int dataLMN(char filename[], int *datavecL, int *datavecM,int *datavecN );

int inzstruct_MH(struct_MH *MH);
int inzstruct_priors(struct_priors *priors);
int inzstruct_data(struct_data *data);
int inzstruct_para(struct_para *para,struct_data *data);

int fillMH(struct_MH *MH);
int filldata(struct_data *D);
int fillpara(struct_para *D_para,struct_data *D);
int fillpriors(struct_priors *D_priors);

#endif
