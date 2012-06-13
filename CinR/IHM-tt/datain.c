#include "headers.h"
#include "datain.h"


/*TEST*/

int testargc(int argc)
{
 	if (argc!=4) {
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


int datadouble(char filename[],char filename2[],struct_data *D,double *QFADyA,double *QFADyB)
{
	int i,t=0;
       	double data,*K_lm,*r_lm,P_a,P_b,size;

	size=D->MAXmn;
	K_lm=malloc(size*sizeof(double)); 
	r_lm=malloc(size*sizeof(double)); 

	/*file=fopen(filename, "r");*/
	t=-1;
	
		for (i=0;i<D->SHIFTmn;i++){ 
			t+=1;
			K_lm[i]=exp(QFADyA[t]); 
		}
		for (i=0;i<2*D->L+3;i++){
			t+=1;
		}
		P_a=exp(QFADyA[t]);  
		for (i=0;i<D->SHIFTmn;i++){
			t+=1;
			r_lm[i]=exp(QFADyA[t]);   
		}

        for (i=0;i<D->SHIFTmn;i++){
	  if(K_lm[i]<=2*P_a){K_lm[i]=2*P_a;}
        }

	for (i=0;i<D->SHIFTmn;i++){
		D->y[i]=(r_lm[i]/log(2*gsl_max(0,K_lm[i]-P_a)/gsl_max(0,K_lm[i]-2*P_a)))*(log(K_lm[i]/P_a)/log(2));
		D->y[i]=D->y[i];
	
	}

		t=-1;

		for (i=D->SHIFTmn;i<D->MAXmn;i++){
			t+=1;
			K_lm[i]=exp(QFADyB[t]);   	
		}
		for (i=0;i<2*D->L+3;i++){
			t+=1;
		}
		P_b=exp(QFADyB[t]);  
		for (i=D->SHIFTmn;i<D->MAXmn;i++){
			t+=1;
			r_lm[i]=exp(QFADyB[t]);   
		}


	for (i=D->SHIFTmn;i<D->MAXmn;i++){
          if(K_lm[i]<=2*P_b){K_lm[i]=2*P_b;}
        }

	for (i=D->SHIFTmn;i<D->MAXmn;i++){
		D->y[i]=(r_lm[i]/log(2*gsl_max(0,K_lm[i]-P_b)/gsl_max(0,K_lm[i]-2*P_b)))*(log(K_lm[i]/P_b)/log(2));
		D->y[i]=(D->y[i]);
	}

return 0;
}

int dataint(char filename[],char filename2[], int datavec[] ,int length)
{

return 0;
}

int dataLMN(char filename[],char filename2[], int *datavecL)
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
    D_priors->Z_mu=PRIORS[0];	

	D_priors->eta_Z_p=PRIORS[1];      

	D_priors->eta_Z=PRIORS[2];	

	D_priors->psi_Z=PRIORS[3];

	D_priors->eta_nu=PRIORS[4];	

	D_priors->psi_nu=PRIORS[5];  

	D_priors->nu_mu=PRIORS[6];    	

	D_priors->eta_nu_p=PRIORS[7]; 

	D_priors->alpha_mu=PRIORS[8];      	

	D_priors->eta_alpha=PRIORS[9];

	D_priors->p=PRIORS[10];    

	D_priors->eta_gamma=PRIORS[11];	

	D_priors->psi_gamma=PRIORS[12];

	D_priors->eta_upsilon=PRIORS[13];	

	D_priors->psi_upsilon=PRIORS[14];

	D_priors->upsilon_mu=PRIORS[15];			
	D_priors->df=3;
	/*fillpriors(priors);*/
return 0;
}

int inzstruct_data(struct_data *data,int *QFAIA,double *QFADyA,int *QFADNoORFA,int *QFAIB,double *QFADyB,int *QFADNoORFB)
{
	int i;
	long size;
		data->L=QFAIA[0];     
		data->L=QFAIB[0];
	size=2*data->L;
	data->NoORF=malloc(size*sizeof(double));    

 for (i=0;i<(data->L);i++){
data->NoORF[i]=QFADNoORFA[i];
data->NoORF[i+data->L]=QFADNoORFB[i];
}


	data->NoSUM=malloc(size*sizeof(double));  
	filldata(data);

	size=data->MAXmn;
	data->y=malloc(size*sizeof(double));  

	datadouble("dataA2.txt","dataB2.txt",data,QFADyA,QFADyB);
return 0;
}

int inzstruct_para(struct_para *para,struct_data *data,struct_priors *priors)
{
	long size;
	size=data->L;
	para->delta_l=malloc(size*sizeof(double));
	para->gamma_cl=malloc(size*sizeof(double));
	para->Z_l=malloc(size*sizeof(double));
	para->nu_l=malloc(size*sizeof(double));
	size=2;
	para->alpha_c=malloc(size*sizeof(double));
	para->upsilon_c=malloc(size*sizeof(double));
	fillpara(para,data,priors);
return 0;
}

/*FILL*/

int fillMH(struct_MH *MH)
{
	MH->halpha_c=0.008;
	MH->hsigma_gamma=0.5;
	MH->hupsilon_c=0.5;
	MH->hsigma_upsilon=0.5;
	MH->hsigma_nu=0.5;
	MH->hsigma_Z=0.5;
	MH->hnu_p=0.5;
	MH->hgamma_cl=0.5;
	MH->hZ_l=0.2;
	MH->hnu_l=0.5;
	MH->accept_Z=0;
	MH->accept_up=0;
	MH->accept_nu=0; /*h sd; accept=0*/
return 0;
}

int filldata(struct_data *D)
{
	int l;
	D->NoSUM[0]=0;
	for (l=1;l<(D->L);l++){
		D->NoSUM[l]=D->NoSUM[l-1]+D->NoORF[l-1];
	}

	D->NoSUM[D->L]=D->SHIFTmn=D->NoSUM[D->L-1]+D->NoORF[D->L-1];

	for (l=(1+D->L);l<(2*D->L);l++){
		D->NoSUM[l]=D->NoSUM[l-1]+D->NoORF[l-1];
	}
	
	D->MAXmn=D->NoSUM[2*D->L-1]+D->NoORF[2*D->L-1];

return 0;
}

int fillpara(struct_para *D_para, struct_data *D,struct_priors *D_priors)
{
  int l,m,mm,ll;
  double SUM=0,SUMa=0,SUMb=0;
	/*initials*/
	for (l=0;l<D->L;l++){
	  for (m=0;m<D->NoORF[l];m++){
	    mm=D->NoSUM[l]+m;
	    SUM += D->y[mm];
	  }
	  /*  if ((SUM/D->NoORF[l])<1){
	    D_para->Z_l[l]=0;
	  }
	  else{*/
	    D_para->Z_l[l]=log(SUM/D->NoORF[l]);
	    /* }*/
	  SUM=0;
	  SUMa+=D_para->Z_l[l];
	}
	D_para->Z_p=(SUMa/D->L);

	D_para->sigma_Z=D_priors->eta_Z;     


	for (l=0;l<D->L;l++)          {D_para->nu_l[l]=D_priors->nu_mu;}          
	D_para->sigma_nu=D_priors->eta_nu;   
  
	D_para->nu_p=D_priors->nu_mu;    

	for (l=0;l<D->L;l++)          {D_para->gamma_cl[l]=0;} 

	for (l=0;l<D->L;l++)          {D_para->delta_l[l]=1;}/*!*/  

	D_para->alpha_c[0]=0;/**/

	for (l=0;l<D->L;l++){
	  SUMa+=exp(D_para->Z_l[l]);
	  ll=l+D->L;
	  for (m=0;m<D->NoORF[ll];m++){
	    mm=D->NoSUM[ll]+m;
	      SUM+=D->y[mm];
	  }
	  if ((SUM/D->NoORF[ll])<1){
	    SUMb+=0;
	  }
	  else{
	    SUMb+=SUM/D->NoORF[ll];
	  }
	  SUM=0;
	}

	D_para->alpha_c[1]=log(SUMb/SUMa);	/**/
	D_para->sigma_gamma=D_priors->eta_gamma;
	D_para->upsilon_c[0]=0; /**/
	D_para->upsilon_c[1]=0;	
	D_para->sigma_upsilon=D_priors->eta_upsilon;

return 0;
}

int fillpriors(struct_priors *D_priors)
{

return 0;
}

/* eof  */
