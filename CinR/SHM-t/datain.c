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

  return 0;
}

int dataint(char filename[], int datavec[], int length)
{

  return 0;
}

int dataLMN(char filename[], int *datavecL,int *datavecM,int *datavecN,int *datavecmaxy,int *datavecmaxTIME)
{
 
  return 0;
}

/*INZ*/

int inzstruct_MH(struct_MH *MH)
{
  fillMH(MH);
  return 0;
}

int inzstruct_priors(struct_priors *D_priors,double *PRIORS)
{
    /*K*/
    D_priors->sigma_K=D_priors->tau_K_mu=PRIORS[0];
   
    D_priors->phi_K=D_priors->eta_tau_K_p=PRIORS[1];          

    D_priors->eta_K_o=PRIORS[2]; 

    D_priors->psi_K_o=PRIORS[3];            
    /*r*/

    D_priors->sigma_r=D_priors->tau_r_mu=PRIORS[4];              

    D_priors->phi_r=D_priors->eta_tau_r_p=PRIORS[5];          

    D_priors->eta_r_o=PRIORS[6];              

    D_priors->psi_r_o=PRIORS[7];       
    /*nu*/

    D_priors->eta_nu=PRIORS[8];         

    D_priors->psi_nu=PRIORS[9];      
    /*K*//*r*//*nu*//*P*/

    D_priors->K_mu=PRIORS[10];     

    D_priors->eta_K_p=PRIORS[11];      /*Normal  LMean; Precisions */

    D_priors->r_mu=PRIORS[12];   

    D_priors->eta_r_p=PRIORS[13];      /*Normal  LMean; Precisions */

    D_priors->nu_mu=PRIORS[14];      

    D_priors->eta_nu_p=PRIORS[15];     /*Normal  LMean; Precisions */

    D_priors->P_mu=PRIORS[16];   

    D_priors->eta_P=PRIORS[17];   /*Normal  LMean; Precisions */
    D_priors->df=3;
    D_priors->df2=5;
    D_priors->eta_tau_K=D_priors->eta_tau_K_p;  D_priors->psi_tau_K=D_priors ->eta_tau_K_p;
    D_priors->eta_tau_r=D_priors->eta_tau_r_p;  D_priors->psi_tau_r=D_priors->eta_tau_r_p;
  /*fillpriors(priors);*/
  return 0;
}

int inzstruct_data(struct_data *data,int *QFAI,double *QFADy,double *QFADx,int *QFADNoORF,int *QFADNoTIME)
{
	int i;
  long size;
	data->L=QFAI[0];
	data->M=QFAI[1];
	data->N=QFAI[2];
	data->maxy=QFAI[3];
	data->maxNoTIME=QFAI[4];
  /*dataLMN("LMNmaxdata.txt",&data->L,&data->M,&data->N,&data->maxy,&data->maxNoTIME);  */
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

 for (i=0;i<(data->L*data->M*data->N);i++){
data->y[i]=QFADy[i];
data->x[i]=QFADx[i];
}
 for (i=0;i<(data->L);i++){
data->NoORF[i]=QFADNoORF[i];
}
  /*datadouble("ydata.txt",data->y,data->L*data->M*data->N);
  datadouble("xdata.txt",data->x,data->L*data->M*data->N);
  dataint("NoORFdata.txt",data->NoORF,data->L);*/

  filldata(data);
  testsame(data->maxNoTIME,data->SHIFTlmn);

  size=data->SHIFTlmn;/*inputfromfile*/
  data->NoTIME=malloc(size*sizeof(double));   /*Cycle with SHIFTlm*/

 for (i=0;i<(data->SHIFTlmn);i++){
data->NoTIME[i]=QFADNoTIME[i];
}
  /*dataint("NoTIMEdata.txt",data->NoTIME,data->SHIFTlmn);*/


  return 0;
}

int inzstruct_para(struct_para *para,struct_data *data,struct_priors *priors)
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

  fillpara(para,data,priors);
  return 0;
}

/*FILL*/

int fillMH(struct_MH *MH)
{
  MH->hK=0.1;MH->accept_K=0;
  MH->hr=0.1;MH->accept_r=0;
  MH->hnu=0.1;MH->accept_nu=0;
  MH->hP=0.2;MH->accept_P=0;  /*h sd; accept=0*/
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

int fillpara(struct_para *D_para, struct_data *D,struct_priors *D_priors)
{
  int l,m,mm;
  double SUM=0,SUMa=0;
  /*initials*/
  /*K*/
  D_para->A=0;
  D_para->B=0;

 for (l=0;l<D->L;l++){
    for (m=0;m<D->NoORF[l];m++){
      mm=D->NoSUM[l]+m;
      if(D->y[l*D->M*D->N + m*D->N + D->NoTIME[mm]-1]<=0){D_para->K_lm[mm]=D_priors->P_mu;SUM+=D_para->K_lm[mm];}
	else{     
	  D_para->K_lm[mm]=gsl_sf_log(D->y[l*D->M*D->N + m*D->N + D->NoTIME[mm]-1]);SUM+=D_para->K_lm[mm];
	}
    }
    D_para->K_o_l[l]=SUM/D->NoORF[l];
    SUM=0;
    SUMa+=D_para->K_o_l[l];
  }
 D_para->K_p=SUMa/D->L;       /*LMean*/

  for (l=0;l<D->L;l++)          {D_para->tau_K_l[l]=D_priors->sigma_K;}                  /*Precision*/
  
  D_para->sigma_K_o=D_para->sigma_K_o_b=D_priors->eta_K_o;               /*Precision*/
  /*r*/
  for (l=0;l<D->L;l++){
    for (m=0;m<D->NoORF[l];m++){
      mm=D->NoSUM[l]+m;
      D_para->r_lm[mm]=D_priors->r_mu;
    }
  }                          /*LMean*/

  for (l=0;l<D->L;l++)          {D_para->tau_r_l[l]=D_priors->sigma_r;}                  /*Precision*/

  for (l=0;l<D->L;l++)          {D_para->r_o_l[l]=D_priors->r_mu;}        /*LMean*/
  D_para->sigma_r_o=D_priors->eta_r_o;               /*Precision*/

  D_para->r_p=D_priors->r_mu;       /*LMean*/
  /*nu*/
  for (l=0;l<D->L;l++)          {D_para->nu_l[l]=D_priors->nu_mu;}                      /*LMean*/
  D_para->sigma_nu=D_priors->eta_nu;   /*Precision for lMean*/

  D_para->nu_p=D_priors->nu_mu;   /*LMean*/
  /*P*/
  D_para->P=D_priors->P_mu;      /*LMean*/


  D_para->tau_K_p=D_priors->tau_K_mu;
  D_para->sigma_tau_K=D_priors->eta_tau_K;
  D_para->tau_r_p=D_priors->tau_r_mu;
  D_para->sigma_tau_r=D_priors->eta_tau_r;

  return 0;
}

int fillpriors(struct_priors *D_priors)
{

  return 0;
}

/* eof  */
