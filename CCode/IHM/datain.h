#ifndef _datain_h
#define _datain_h

int testargc(int argc);

int datadouble(char filename[], char filename2[], double datavec[],int length,struct_data *D);
int dataint(char filename[],char filename2[], int datavec[] ,int length);
int dataLMN(char filename[],char filename2[], int *datavecL);

int inzstruct_MH(struct_MH *MH);
int inzstruct_priors(struct_priors *priors);
int inzstruct_data(struct_data *data);
int inzstruct_para(struct_para *para,struct_data *data);

int fillMH(struct_MH *MH);
int filldata(struct_data *D);
int fillpara(struct_para *D_para,struct_data *D);
int fillpriors(struct_priors *D_priors);
double meanParticleMatrix(gsl_matrix*,int col,int start,int end);
gsl_matrix *getParticleMatrix(char *filename);
#endif
