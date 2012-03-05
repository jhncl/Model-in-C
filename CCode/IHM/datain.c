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


int datadouble(char filename[],char filename2[],struct_data *D)
{
	int i;
	char number[20];
	double data,*K_lm,*r_lm,P_a,P_b,size;
	FILE *file;
	FILE *file2;
	size=D->MAXmn;
	K_lm=malloc(size*sizeof(double)); 
	r_lm=malloc(size*sizeof(double)); 

	file=fopen(filename, "r");
	i=0;
	if ( file != NULL ){
		for (i=0;i<2*D->SHIFTmn+5*D->L+7+4;i++){
			fscanf(file, "%c %s",number,number);
		}
		for (i=0;i<D->SHIFTmn;i++){
			fscanf(file, "%c %lf",number,&data);
			K_lm[i]=exp(data);   
		}
		for (i=0;i<2*D->L+3;i++){
			fscanf(file, "%c %lf",number,&data);
		}
		P_a=exp(data);  
		for (i=0;i<D->SHIFTmn;i++){
			fscanf(file, "%c %lf",number,&data);
			r_lm[i]=exp(data);   
		}
	}
	else{perror(filename);}
	fclose(file);

	for (i=0;i<D->SHIFTmn;i++){
		D->y[i]=(r_lm[i]/log(2*gsl_max(0,K_lm[i]-P_a)/gsl_max(0,K_lm[i]-2*P_a)))*(log(K_lm[i]/P_a)/log(2));
	}

	file2=fopen(filename2, "r");
	i=0;
	if ( file != NULL ){
		for (i=0;i<2*D->SHIFTmn+5*D->L+7+4;i++){
			fscanf(file, "%c %s",number,number);
		}
		for (i=D->SHIFTmn;i<D->MAXmn;i++){
			fscanf(file, "%c %lf",number,&data);
			K_lm[i]=exp(data);   
		}
		for (i=0;i<2*D->L+3;i++){
			fscanf(file, "%c %lf",number,&data);
		}
		P_b=exp(data);  
		for (i=D->SHIFTmn;i<D->MAXmn;i++){
			fscanf(file, "%c %lf",number,&data);
			r_lm[i]=exp(data);   
		}
	}
	else{perror(filename2);}
	fclose(file2);

	for (i=D->SHIFTmn;i<D->MAXmn;i++){
		D->y[i]=(r_lm[i]/log(2*gsl_max(0,K_lm[i]-P_b)/gsl_max(0,K_lm[i]-2*P_b)))*(log(K_lm[i]/P_b)/log(2));
	}
return 0;
}

int dataint(char filename[],char filename2[], int datavec[] ,int length)
{
	int i=0;
	char number[20];
	double data;
	FILE *file;
	FILE *file2;
	file = fopen(filename, "r");
	if ( file != NULL ){
		fscanf(file, "%s %lf",number,&data);
		while (fscanf(file, "%s %lf",number,&data)!=-1){
			datavec[i]=data;
			i++;
        	}
	testsame(length,i);
    	}
	else{perror(filename);}
	fclose(file);
	file2 = fopen(filename2, "r");
  	if ( file2 != NULL ){
		fscanf(file2, "%s %lf",number,&data);
		while ( fscanf(file2, "%s %lf",number,&data)!=-1)
	        {
			datavec[i]=data;
			i++;
	        }
	testsame(2*length,i);
	}
	else{perror(filename2);}
	fclose(file2);

return 0;
}

int dataLMN(char filename[],char filename2[], int *datavecL)
{
	char number[20];
	double data;
	FILE *file;
	FILE *file2;
	file=fopen(filename, "r");
	if ( file != NULL ){
		fscanf(file, "%s %lf",number,&data);
		fscanf(file, "%s %lf",number,&data);
		*datavecL=data;     
	}
	else{perror(filename);}
	fclose(file);
	file2=fopen(filename2, "r");
	if ( file2 != NULL ){
		fscanf(file2, "%s %lf",number,&data);
		fscanf(file2, "%s %lf",number,&data);

		testsame(*datavecL,data);
		*datavecL=data;     
	}
	else{perror(filename2);}
	fclose(file2);
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
	dataLMN("LMNmaxdataA1.txt","LMNmaxdataB1.txt",&data->L);      
	size=2*data->L;
	data->NoORF=malloc(size*sizeof(double));    
        dataint("NoORFdataA1.txt","NoORFdataB1.txt",data->NoORF,data->L);

	data->NoSUM=malloc(size*sizeof(double));  
	filldata(data);

	size=data->MAXmn;
	data->y=malloc(size*sizeof(double));  
	datadouble("dataA1.txt","dataB1.txt",data);
	/*if (data->y==NULL||data->x==NULL||data->NoORF==NULL||data->NoSUM==NULL||data->NoTIME==NULL) {
		perror("malloc failed");
    		exit(EXIT_FAILURE);
  	};*/
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
	MH->halpha_c=0.5;
	MH->hsigma_gamma=0.5;
	MH->hupsilon_c=0.5;
	MH->hsigma_upsilon=0.5;
	MH->hsigma_nu=0.5;
	MH->hsigma_Z=0.5;
	MH->hnu_p=0.5;
	MH->hgamma_cl=0.5;
	MH->hZ_l=0.5;
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
	int l;

	/*initials*/
	for (l=0;l<D->L;l++){D_para->Z_l[l]=D_priors->Z_mu;}
	D_para->sigma_Z=D_priors->eta_Z;     
	D_para->Z_p=D_priors->Z_mu;  

	for (l=0;l<D->L;l++)          {D_para->nu_l[l]=D_priors->nu_mu;}          
	D_para->sigma_nu=D_priors->eta_nu;   
  
	D_para->nu_p=D_priors->nu_mu;    

	for (l=0;l<D->L;l++)          {D_para->gamma_cl[l]=0;} 

	for (l=0;l<D->L;l++)          {D_para->delta_l[l]=1;}/*!*/  

	D_para->alpha_c[0]=0;/**/
	D_para->alpha_c[1]=D_priors->alpha_mu;
	D_para->sigma_gamma=D_priors->eta_gamma;
	D_para->upsilon_c[0]=0; /**/
	D_para->upsilon_c[1]=0;	
	D_para->sigma_upsilon=D_priors->eta_upsilon;

return 0;
}

int fillpriors(struct_priors *D_priors)
{
/*priors*/
  char number[20];
  double data;
  FILE *file = fopen("priors.txt", "r");
  if ( file != NULL ){  
    fscanf(file, "%s %lf",number,&data);
    fscanf(file, "%s %lf",number,&data);
	D_priors->Z_mu=data;	
    fscanf(file, "%s %lf",number,&data);
	D_priors->eta_Z_p=data;      
    fscanf(file, "%s %lf",number,&data);
	D_priors->eta_Z=data;	
    fscanf(file, "%s %lf",number,&data);
	D_priors->psi_Z=data;
    fscanf(file, "%s %lf",number,&data);
	D_priors->eta_nu=data;	
    fscanf(file, "%s %lf",number,&data);
	D_priors->psi_nu=data;  
    fscanf(file, "%s %lf",number,&data);
	D_priors->nu_mu=data;    	
    fscanf(file, "%s %lf",number,&data);
	D_priors->eta_nu_p=data; 
    fscanf(file, "%s %lf",number,&data);
	D_priors->alpha_mu=data;      	
    fscanf(file, "%s %lf",number,&data);
	D_priors->eta_alpha=data;
    fscanf(file, "%s %lf",number,&data);
	D_priors->p=data;    
    fscanf(file, "%s %lf",number,&data);
	D_priors->eta_gamma=data;	
    fscanf(file, "%s %lf",number,&data);
	D_priors->psi_gamma=data;
    fscanf(file, "%s %lf",number,&data);
	D_priors->eta_upsilon=data;	
    fscanf(file, "%s %lf",number,&data);
	D_priors->psi_upsilon=data;
    fscanf(file, "%s %lf",number,&data);
	D_priors->upsilon_mu=data;			
  }
  else{perror("Priors");}  
  fclose(file);
return 0;
}

/* eof  */
