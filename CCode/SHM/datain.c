#include "headers.h"
#include "datain.h"


/*TEST*/

int testargc(int argc)
{
 	if (argc!=5) {
    		perror("argc failed");
    		exit(EXIT_FAILURE);
  	}
return 0;
}

int testsame(int a,int b)
{
 	if (a!=b) {
    		perror("data int failed");
    		exit(EXIT_FAILURE);
  	}
return 0;
}



/*READ IN*/


int datadouble(char filename[], double datavec[],int length )
{
	int i;
	char number[20];
	double data;
	FILE *file;
	file=fopen(filename, "r");
	i=0;
	if ( file != NULL ){
		fscanf(file, "%s %lf",number,&data);
		while ( fscanf(file, "%s %lf",number,&data)!=-1){
			datavec[i]=data;   
			i++;
		}
	}
	else{perror(filename);}
	testsame(length,i);
	fclose(file);
return 0;
}

int dataint(char filename[], int datavec[], int length)
{
	int i;
	char number[20];
	double data;
	FILE *file;
	file = fopen(filename, "r");
	i=0;
	if ( file != NULL ){
		fscanf(file, "%s %lf",number,&data);
		while (fscanf(file, "%s %lf",number,&data)!=-1){
			datavec[i]=data;
			i++;
        	}
    	}
	else{perror(filename);}
	testsame(length,i);
	fclose(file);
return 0;
}

int dataLMN(char filename[], int *datavecL,int *datavecM,int *datavecN,int *datavecmaxy,int *datavecmaxTIME)
{
	char number[20];
	double data;
	FILE *file = fopen(filename, "r");
	if ( file != NULL ){
		fscanf(file, "%s %lf",number,&data);
		fscanf(file, "%s %lf",number,&data);
		*datavecL=data;     
		fscanf(file, "%s %lf",number,&data);
		*datavecM=data;
		fscanf(file, "%s %lf",number,&data);
		*datavecN=data;
		fscanf(file, "%s %lf",number,&data);
		*datavecmaxy=data;
		fscanf(file, "%s %lf",number,&data);
		*datavecmaxTIME=data;
	}
	else{perror(filename);}
	fclose(file);
return 0;
}

/*INZ*/

int inzstruct_MH(struct_MH *MH)
{
	fillMH(MH);
return 0;
}

int inzstruct_priors(struct_priors *priors)
{
	fillpriors(priors);
return 0;
}

int inzstruct_data(struct_data *data)
{
	long size;
	dataLMN("LMNmaxdata.txt",&data->L,&data->M,&data->N,&data->maxy,&data->maxNoTIME);  
	testsame(data->L*data->M*data->N,data->maxy);

	size=data->L*data->M*data->N; /*input from file*/ 
  	data->y=malloc(size*sizeof(double));        /*Cycle with SHIFTlmn*/
        data->x=malloc(size*sizeof(double));        /*Cycle with SHIFTlmn*/
	size=data->L;
	data->NoORF=malloc(size*sizeof(double));    /*Cycle with data->L*/
	data->NoSUM=malloc(size*sizeof(double));    /*Cycle with data->L*/

	/*if (data->y==NULL||data->x==NULL||data->NoORF==NULL||data->NoSUM==NULL||data->NoTIME==NULL) {
		perror("malloc failed");
    		exit(EXIT_FAILURE);
  	}*/

	datadouble("ydata.txt",data->y,data->L*data->M*data->N);
        datadouble("xdata.txt",data->x,data->L*data->M*data->N);
        dataint("NoORFdata.txt",data->NoORF,data->L);

	filldata(data);
	testsame(data->maxNoTIME,data->SHIFTlmn);

	size=data->SHIFTlmn;/*inputfromfile*/
	data->NoTIME=malloc(size*sizeof(double));   /*Cycle with SHIFTlm*/
        dataint("NoTIMEdata.txt",data->NoTIME,data->SHIFTlmn);


return 0;
}

int inzstruct_para(struct_para *para,struct_data *data)
{
	long size;

	size=data->SHIFTlmn;
	para->K_lm=malloc(size*sizeof(double));
	para->r_lm=malloc(size*sizeof(double));
	size=data->L;
	para->tau_K_l=malloc(size*sizeof(double));
	para->tau_r_l=malloc(size*sizeof(double));
	para->K_o_l=malloc(size*sizeof(double));
	para->r_o_l=malloc(size*sizeof(double));
	para->nu_l=malloc(size*sizeof(double));

	fillpara(para,data);
return 0;
}

/*FILL*/

int fillMH(struct_MH *MH)
{
	MH->hK=0.1;	MH->accept_K=0;
	MH->hr=0.1;	MH->accept_r=0;
	MH->hnu=0.1;	MH->accept_nu=0;
	MH->hP=0.2;	MH->accept_P=0;  /*h sd; accept=0*/
return 0;
}

int filldata(struct_data *D)
{
	int l;

	D->NoSUM[0]=0;
	for (l=1;l<(D->L);l++){
		D->NoSUM[l]=D->NoSUM[l-1]+D->NoORF[l-1];
	}
	D->SHIFTlmn=D->NoSUM[D->L-1]+D->NoORF[D->L-1];/*create mnSHIFT*/
return 0;
}

int fillpara(struct_para *D_para, struct_data *D)
{
int l,m,mm;
	/*initials*/
	/*K*/
	for (l=0;l<D->L;l++){
		for (m=0;m<D->NoORF[l];m++){
			mm=D->NoSUM[l]+m;
			D_para->K_lm[mm]=D->y[l*D->M*D->N + m*D->N + D->NoTIME[mm]-1];
			if(D_para->K_lm[mm]>0){D_para->K_lm[mm]=gsl_sf_log(D_para->K_lm[mm]);}
		}
	}

	for (l=0;l<D->L;l++)          {D_para->tau_K_l[l]=7;}                  /*Precision*/

	for (l=0;l<D->L;l++)          {D_para->K_o_l[l]=gsl_sf_log(0.25);}        /*LMean*/
	D_para->sigma_K_o=6;               /*Precision*/
	D_para->K_p=gsl_sf_log(0.1);       /*LMean*/

	/*r*/
	for (l=0;l<D->L;l++){
		for (m=0;m<D->NoORF[l];m++){
			mm=D->NoSUM[l]+m;
			D_para->r_lm[mm]=gsl_sf_log(2.5);
		}
	}                          /*LMean*/

	for (l=0;l<D->L;l++)          {D_para->tau_r_l[l]=-1;}                  /*Precision*/

	for (l=0;l<D->L;l++)          {D_para->r_o_l[l]=gsl_sf_log(2.5);}        /*LMean*/
	D_para->sigma_r_o=0;               /*Precision*/

	D_para->r_p=gsl_sf_log(2.5);       /*LMean*/

	/*nu*/
	for (l=0;l<D->L;l++)          {D_para->nu_l[l]=18;}                      /*LMean*/
	D_para->sigma_nu=-2.5;   /*Precision for lMean*/

	D_para->nu_p=18;   /*LMean*/
	/*P*/
	D_para->P=gsl_sf_log(0.0001);      /*LMean*/
return 0;
}

int fillpriors(struct_priors *D_priors)
{
	/*Priors*/
	/*K*/
	D_priors->sigma_K=7;               D_priors->phi_K=1.3;               /*Gamma  Shape; Scale */
	D_priors->eta_K_o=8;               D_priors->psi_K_o=1;             /*Gamma  Shape; Scale */
	/*r*/
	D_priors->sigma_r=-1;               D_priors->phi_r=1.2;               /*Gamma  Shape; Scale */
	D_priors->eta_r_o=1;               D_priors->psi_r_o=1;             /*Gamma  Shape; Scale */
	/*nu*/
	D_priors->eta_nu=-1;              D_priors->psi_nu=1;              /*Gamma  Shape; Scale */

	/*K*//*r*//*nu*//*P*/
	D_priors->K_mu=gsl_sf_log(0.2192928);      D_priors->eta_K_p=1;      /*Normal  LMean; Precisions */
	D_priors->r_mu=gsl_sf_log(2.5);            D_priors->eta_r_p=1;      /*Normal  LMean; Precisions */
	D_priors->nu_mu=gsl_sf_log(31);            D_priors->eta_nu_p=1;     /*Normal  LMean; Precisions */
	D_priors->P_mu=gsl_sf_log(0.0002);         D_priors->eta_P=1/0.01;   /*Normal  LMean; Precisions */
return 0;
}

/* eof  */
