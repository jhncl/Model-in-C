#include "headers_JHM.h"
#include "datain_JHM.h"


/*TEST*/

int testargc_JHM(int argc)
{
 	if (argc!=4) {
    		perror("argc failed");
    		exit(EXIT_FAILURE);
  	}
return 0;
}

int testsame_JHM(int a,int b)
{
 	if (a!=b) {
    		perror("data int failed");
    		exit(EXIT_FAILURE);
  	}
return 0;
}

/*READ IN*/


int datadouble_JHM(char filename[], char filename2[], double datavec[],int length)
{

return 0;
}

int dataint_JHM(char filename[],char filename2[], int datavec[] ,int lengtha,int lengthb )
{


return 0;
}

int dataLMN_JHM(char filename[],char filename2[], int *datavecL,int *datavecM,int *datavecN,int *datavecmaxy,int *datavecmaxTIMEa,int *datavecmaxTIMEb)
{
	
return 0;
}

/*INZ*/

int inzstruct_MH_JHM(struct_MH_JHM *MH)
{
	fillMH_JHM(MH);
return 0;
}

int inzstruct_priors_JHM(struct_priors_JHM *D_priors,double *PRIORS)
{
	/*K*/
	D_priors->tau_K_mu=PRIORS[0];        

       D_priors->eta_tau_K_p=PRIORS[1];             

	D_priors->eta_K_o=PRIORS[2];   

            D_priors->psi_K_o=PRIORS[3];        
	/*r*/

	D_priors->tau_r_mu=PRIORS[4];           

    D_priors->eta_tau_r_p=PRIORS[5];              

	D_priors->eta_r_o=PRIORS[6];            

   D_priors->psi_r_o=PRIORS[7];         
	/*nu*/

	D_priors->eta_nu=PRIORS[8];    

          D_priors->psi_nu=PRIORS[9];           

	/*K*//*r*//*nu*//*P*/

 D_priors->K_mu=PRIORS[10];   

   D_priors->eta_K_p=PRIORS[11];     

 D_priors->r_mu=PRIORS[12];       

     D_priors->eta_r_p=PRIORS[13];      

	D_priors->nu_mu=PRIORS[14];           

 D_priors->eta_nu_p=PRIORS[15];    

	D_priors->P_mu=PRIORS[16];     

   D_priors->eta_P=PRIORS[17];   

	/*data2.c*/       

	D_priors->alpha_mu=PRIORS[19];     

     D_priors->eta_alpha=PRIORS[20];

	D_priors->beta_mu=PRIORS[21];        

   D_priors->eta_beta=PRIORS[22];

	D_priors->p=PRIORS[23];    

	D_priors->eta_gamma=PRIORS[24];

	D_priors->psi_gamma=PRIORS[25];

	D_priors->eta_omega=PRIORS[26];	

D_priors->psi_omega=PRIORS[27];

	D_priors->eta_upsilon=PRIORS[28];	

 	   D_priors->psi_upsilon=PRIORS[29];	    

	D_priors->upsilon_mu=PRIORS[30];	
	D_priors->df=3;
	D_priors->df2=5;
	D_priors->eta_tau_K=D_priors->eta_tau_K_p;  D_priors->psi_tau_K=D_priors->eta_tau_K_p;
	D_priors->eta_tau_r=D_priors->eta_tau_r_p;  D_priors->psi_tau_r=D_priors->eta_tau_r_p;
	/*fillpriors(priors);*/
return 0;
}

int inzstruct_data_JHM(struct_data_JHM *data,int *QFAIA,double *QFADyA,double *QFADxA,int *QFADNoORFA,int *QFADNoTIMEA,int *QFAIB,double *QFADyB,double *QFADxB,int *QFADNoORFB,int *QFADNoTIMEB)
{
	int i;
	long size;

		data->L=QFAIA[0];     

		data->M=QFAIA[1];
		data->N=QFAIA[2];
		data->maxy=QFAIA[3];
		data->maxTIMEa=QFAIA[4];
		data->L=QFAIB[0];
		data->M=QFAIB[1];
		data->N=QFAIB[2];
		data->maxy=QFAIB[3];
		data->maxTIMEb=QFAIB[4];

		/*dataLMN("LMNmaxdataA1.txt","LMNmaxdataB1.txt",
		&data->L,&data->M,&data->N,&data->maxy,&data->maxTIMEa,&data->maxTIMEb);*/     
		
	data->SHIFTlmn=data->maxy;
	size=data->maxy*2;
  	data->y=malloc(size*sizeof(double));  
        data->x=malloc(size*sizeof(double));  
	size=data->L*2;
	data->NoORF=malloc(size*sizeof(double));  
	data->NoSUM=malloc(size*sizeof(double));
	size=data->maxTIMEa+data->maxTIMEb;
	data->NoTIME=malloc(size*sizeof(double));  
	
	if (data->y==NULL||data->x==NULL||data->NoORF==NULL||data->NoSUM==NULL||data->NoTIME==NULL) {
		perror("malloc failed");
    		exit(EXIT_FAILURE);
  	}
 for (i=0;i<(data->maxy);i++){
data->y[i]=QFADyA[i];
data->x[i]=QFADxA[i];
data->y[i+data->maxy]=QFADyB[i];
data->x[i+data->maxy]=QFADxB[i];
}

 for (i=0;i<(data->L);i++){
data->NoORF[i]=QFADNoORFA[i];
data->NoORF[i+data->L]=QFADNoORFB[i];
}
 for (i=0;i<(data->maxTIMEa);i++){
data->NoTIME[i]=QFADNoTIMEA[i];
}
 for (i=0;i<(data->maxTIMEb);i++){
data->NoTIME[i+data->maxTIMEa]=QFADNoTIMEA[i];
}
	/*datadouble("ydataA1.txt","ydataB1.txt",data->y,data->maxy);
        datadouble("xdataA1.txt","xdataB1.txt",data->x,data->maxy);
        dataint("NoORFdataA1.txt","NoORFdataB1.txt",data->NoORF,data->L,data->L);
        dataint("NoTIMEdataA1.txt","NoTIMEdataB1.txt",data->NoTIME,data->maxTIMEa,data->maxTIMEb);*/

filldata_JHM(data);
return 0;
}

int inzstruct_para_JHM(struct_para_JHM *para,struct_data_JHM *data,struct_priors_JHM *priors)
{
	long size;
	size=data->L*2;
	para->tau_K_cl=malloc(size*sizeof(double));
	para->tau_r_cl=malloc(size*sizeof(double));
	size=data->maxTIMEa+data->maxTIMEb;/*inputfromfile*/
	para->K_clm=malloc(size*sizeof(double));
	para->r_clm=malloc(size*sizeof(double));
	size=data->L;
	para->delta_l=malloc(size*sizeof(double));
	para->gamma_cl=malloc(size*sizeof(double));
	para->omega_cl=malloc(size*sizeof(double));
	para->K_o_l=malloc(size*sizeof(double));
	para->r_o_l=malloc(size*sizeof(double));
	para->nu_l=malloc(size*sizeof(double));
	size=2;
	para->alpha_c=malloc(size*sizeof(double));
	para->beta_c=malloc(size*sizeof(double));
	para->upsilon_c=malloc(size*sizeof(double));

	para->tau_K_p=malloc(size*sizeof(double));
	  para->tau_r_p=malloc(size*sizeof(double));
	  para->sigma_tau_K=malloc(size*sizeof(double));
	  para->sigma_tau_r=malloc(size*sizeof(double));

	fillpara_JHM(para,data,priors);
return 0;
}

/*FILL*/

int fillMH_JHM(struct_MH_JHM *MH)
{
	MH->hK=0.1;	MH->accept_K=0;
	MH->hr=0.1;	MH->accept_r=0;
	MH->hnu=0.4;	MH->accept_nu=0;
	MH->hP=0.2;	MH->accept_P=0;  /*h sd; accept=0*/
	MH->halpha=0.01;
return 0;
}

int filldata_JHM(struct_data_JHM *D)
{
	int l;
	D->NoSUM[0]=0;
	for (l=1;l<(D->L);l++){
		D->NoSUM[l]=D->NoSUM[l-1]+D->NoORF[l-1];
	}

	D->NoSUM[D->L]=D->SHIFTmn=D->NoSUM[D->L-1]+D->NoORF[D->L-1];/*create mnSHIFT*/

	for (l=(1+D->L);l<(2*D->L);l++){
		D->NoSUM[l]=D->NoSUM[l-1]+D->NoORF[l-1];
	}


return 0;
}

int fillpara_JHM(struct_para_JHM *D_para, struct_data_JHM *D,struct_priors_JHM *D_priors)
{
int c,l,m,ll,mm;
 double SUM=0,SUMa=0,SUMb=0;

 D_para->A=0;
 D_para->B=0;

for (c=0;c<2;c++){
    for (l=0;l<D->L;l++){
      ll=c*D->L+l;
      for (m=0;m<D->NoORF[ll];m++){
	mm=D->NoSUM[ll]+m;
	if (D->y[c*D->SHIFTlmn+l*D->M*D->N + m*D->N + D->NoTIME[mm]-1]<=0){D_para->K_clm[mm]=D_priors->P_mu;}
	else{     
		D_para->K_clm[mm]=log(D->y[c*D->SHIFTlmn+l*D->M*D->N + m*D->N + D->NoTIME[mm]-1]);
		}
	if (c==0){SUMa+=exp(D_para->K_clm[mm]);SUM+=D_para->K_clm[mm];}
	  else{ SUMb+=exp(D_para->K_clm[mm]);}
      }
      if(c==0){D_para->K_o_l[l]=SUM/D->NoORF[l];}
      SUM=0;
      }
  }
 D_para->alpha_c[1]=log((SUMb/(2*D->maxy-D->SHIFTlmn))/(SUMa/D->SHIFTlmn));
 D_para->beta_c[1]=D_para->alpha_c[1];
 SUM=0;
 for (l=0;l<(D->L);l++){SUM+=D_para->K_o_l[l];}
 D_para->K_p=SUM/D->L;
				
	for (l=0;l<(2*D->L);l++)          {D_para->tau_K_cl[l]=D_priors->tau_K_mu;}             
       
	D_para->sigma_K_o=D_priors->eta_K_o;              


	for (c=0;c<2;c++){
		for (l=0;l<D->L;l++){
			ll=c*D->L+l;
			for (m=0;m<D->NoORF[l];m++){
				mm=D->NoSUM[ll]+m;
				D_para->r_clm[mm]=D_priors->r_mu;   
			} 
		}  
	}
               

	for (l=0;l<2*D->L;l++)          {D_para->tau_r_cl[l]=D_priors->tau_r_mu;}           

	for (l=0;l<D->L;l++)          {D_para->r_o_l[l]=D_priors->r_mu;}     	   
	D_para->sigma_r_o=D_priors->eta_r_o;          

	D_para->r_p=D_priors->r_mu;    
	
	
	for (l=0;l<D->L;l++)          {D_para->nu_l[l]=D_priors->nu_mu;}                     
	D_para->sigma_nu=D_priors->eta_nu;   

	D_para->nu_p=D_priors->nu_mu; 
	
  D_para->P=D_priors->P_mu;   

	for (l=0;l<D->L;l++)          {D_para->gamma_cl[l]=0;} 

	for (l=0;l<D->L;l++)          {D_para->omega_cl[l]=0;}
	for (l=0;l<D->L;l++)          {D_para->delta_l[l]=1;}/*!*/  
 
	D_para->alpha_c[0]=0;
	D_para->beta_c[0]=0;
	D_para->sigma_gamma=D_priors->eta_gamma;
	D_para->sigma_omega=D_priors->eta_omega;
	D_para->upsilon_c[0]=0; 
	D_para->upsilon_c[1]=D_priors->upsilon_mu;      
        D_para->sigma_upsilon=D_priors->eta_upsilon;
	for (c=0;c<2;c++){
	D_para->tau_K_p[c]=D_priors->tau_K_mu;
	D_para->sigma_tau_K[c]=D_priors->eta_tau_K;
	D_para->tau_r_p[c]=D_priors->tau_r_mu;
	D_para->sigma_tau_r[c]=D_priors->eta_tau_r;
	}
return 0;
}

int fillpriors_JHM(struct_priors_JHM *D_priors)
{
 
return 0;
}
