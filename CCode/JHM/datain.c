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

/*READ IN*/


int datadouble(char filename[], char filename2[], double datavec[],int length)
{
	int i;
	char number[20];
	double data;
	FILE *file;
	FILE *file2;
	file=fopen(filename, "r");
	i=0;
	if ( file != NULL ){
		fscanf(file, "%s %lf",number,&data);
		while ( fscanf(file, "%s %lf",number,&data)!=-1){
			datavec[i]=data;   
			i++;
		}
	testsame(length,i);
	}
	else{perror(filename);}
	fclose(file);

	file2=fopen(filename2, "r");
	if ( file != NULL ){
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

int dataint(char filename[],char filename2[], int datavec[] ,int lengtha,int lengthb )
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
	testsame(lengtha,i);
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
	testsame(lengtha+lengthb,i);
	}
	else{perror(filename2);}
	fclose(file2);

return 0;
}

int dataLMN(char filename[],char filename2[], int *datavecL,int *datavecM,int *datavecN,int *datavecmaxy,int *datavecmaxTIMEa,int *datavecmaxTIMEb)
{
	char number[20];
	double data;
	FILE *file;
	FILE *file2;
	file = fopen(filename, "r");
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
		*datavecmaxTIMEa=data;
	}
	else{perror(filename);}
	fclose(file);
	file2 = fopen(filename2, "r");
	if ( file2 != NULL ){
		fscanf(file2, "%s %lf",number,&data);
		fscanf(file2, "%s %lf",number,&data);
		testsame(*datavecL,data); 	
		*datavecL=data;     

		fscanf(file2, "%s %lf",number,&data);
		*datavecM=gsl_max(*datavecM,data);

		fscanf(file2, "%s %lf",number,&data);
		*datavecN=gsl_max(*datavecN,data);

		fscanf(file2, "%s %lf",number,&data);
		testsame(*datavecmaxy,data); 
		*datavecmaxy=data;
		testsame(*datavecmaxy,(*datavecL)*(*datavecM)*(*datavecN));

		fscanf(file2, "%s %lf",number,&data);
		*datavecmaxTIMEb=data;
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
		dataLMN("LMNmaxdataA1.txt","LMNmaxdataB1.txt",
		&data->L,&data->M,&data->N,&data->maxy,&data->maxTIMEa,&data->maxTIMEb);     
		
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

	datadouble("ydataA1.txt","ydataB1.txt",data->y,data->maxy);
        datadouble("xdataA1.txt","xdataB1.txt",data->x,data->maxy);
        dataint("NoORFdataA1.txt","NoORFdataB1.txt",data->NoORF,data->L,data->L);
        dataint("NoTIMEdataA1.txt","NoTIMEdataB1.txt",data->NoTIME,data->maxTIMEa,data->maxTIMEb);

	filldata(data);
return 0;
}

int inzstruct_para(struct_para *para,struct_data *data,struct_priors *priors)
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

	fillpara(para,data,priors);
return 0;
}

/*FILL*/

int fillMH(struct_MH *MH)
{
	MH->hK=0.0001;	MH->accept_K=0;
	MH->hr=0.001;	MH->accept_r=0;
	MH->hnu=0.3;	MH->accept_nu=0;
	MH->hP=0.1;	MH->accept_P=0;  /*h sd; accept=0*/
return 0;
}

int filldata(struct_data *D)
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

int fillpara(struct_para *D_para, struct_data *D,struct_priors *D_priors)
{
int c,l,m,ll,mm;
	/*initials*/

	/*K*/
	for (c=0;c<2;c++){
		for (l=0;l<D->L;l++){
			ll=c*D->L+l;
			for (m=0;m<D->NoORF[ll];m++){
				mm=D->NoSUM[ll]+m;
				D_para->K_clm[mm]=D->y[c*D->SHIFTlmn+l*D->M*D->N + m*D->N + D->NoTIME[mm]-1];
				if(D_para->K_clm[mm]>0){
					D_para->K_clm[mm]=D_priors->K_mu;
				}
			}
		}
	}

				
	for (l=0;l<(2*D->L);l++)          {D_para->tau_K_cl[l]=D_priors->sigma_K;}                  /*Precision*/

	for (l=0;l<D->L;l++)          {D_para->K_o_l[l]=D_priors->K_mu;}        /*LMean*/
	D_para->sigma_K_o=D_priors->eta_K_o;               /*Precision*/
	D_para->K_p=D_priors->K_mu;       /*LMean*/

	/*r*/
	for (c=0;c<2;c++){
		for (l=0;l<D->L;l++){
			ll=c*D->L+l;
			for (m=0;m<D->NoORF[l];m++){
				mm=D->NoSUM[ll]+m;
				D_para->r_clm[mm]=D_priors->r_mu;   
			} 
		}  
	}
/*
	for (c=0;c<2;c++){
		for (l=0;l<D->L;l++){
			ll=c*D->L+l;
			for (m=0;m<D->NoORF[ll];m++){
				mm=D->NoSUM[ll]+m;
				nn=c*D->SHIFTlmn+l*D->M*D->N + m*D->N + D->NoTIME[mm]/2-1;
				if(exp(D_para->K_clm[mm])>D->y[nn]){
					D_para->r_clm[mm]=((D->y[nn]*(exp(D_para->K_clm[mm])-exp(D_para->P)))/(exp(D_para->P)*
						(exp(D_para->K_clm[mm])-D->y[nn]) ) )/D->x[nn];
				}
				else{
					D_para->r_clm[mm]=0.2;
				}
				if(D_para->r_clm[mm]!=0){
					D_para->r_clm[mm]=gsl_sf_log(D_para->r_clm[mm]);
				}
			}
		}
	} */  

               

	for (l=0;l<2*D->L;l++)          {D_para->tau_r_cl[l]=D_priors->sigma_r;}                  /*Precision*/

	for (l=0;l<D->L;l++)          {D_para->r_o_l[l]=D_priors->r_mu;}     	   /*LMean*/
	D_para->sigma_r_o=D_priors->eta_r_o;               /*Precision*/

	D_para->r_p=D_priors->r_mu;       /*LMean*/

	/*nu*/
	for (l=0;l<D->L;l++)          {D_para->nu_l[l]=D_priors->nu_mu;}                      /*LMean*/
	D_para->sigma_nu=D_priors->eta_nu;   /*Precision for lMean*/

	D_para->nu_p=D_priors->nu_mu;   /*LMean*/
	/*P*/
  D_para->P=D_priors->P_mu;      /*LMean*/

	for (l=0;l<D->L;l++)          {D_para->gamma_cl[l]=0;} 

	for (l=0;l<D->L;l++)          {D_para->omega_cl[l]=0;}
	for (l=0;l<D->L;l++)          {D_para->delta_l[l]=1;}/*!*/  
 
	D_para->alpha_c[0]=0;
	D_para->beta_c[0]=0;
	D_para->alpha_c[1]=D_priors->alpha_mu;
	D_para->beta_c[1]=D_priors->beta_mu;  
	D_para->sigma_gamma=D_priors->eta_gamma;
	D_para->sigma_omega=D_priors->eta_omega;
	D_para->upsilon_c[0]=0; 
	D_para->upsilon_c[1]=D_priors->upsilon_mu;      
        D_para->sigma_upsilon=D_priors->eta_upsilon;
return 0;
}

int fillpriors(struct_priors *D_priors)
{
 /*Priors*/
  char number[20];
  double data;
  FILE *file = fopen("priors.txt", "r");
  if ( file != NULL ){  
    fscanf(file, "%s %lf",number,&data);
    fscanf(file, "%s %lf",number,&data);
	/*K*/
	D_priors->sigma_K=data;        
fscanf(file, "%s %lf",number,&data);
       D_priors->phi_K=data;             
fscanf(file, "%s %lf",number,&data);
	D_priors->eta_K_o=data;   
fscanf(file, "%s %lf",number,&data);
            D_priors->psi_K_o=data;        
	/*r*/
fscanf(file, "%s %lf",number,&data);
	D_priors->sigma_r=data;           
fscanf(file, "%s %lf",number,&data);
    D_priors->phi_r=data;              
fscanf(file, "%s %lf",number,&data);
	D_priors->eta_r_o=data;            
fscanf(file, "%s %lf",number,&data);
   D_priors->psi_r_o=data;         
	/*nu*/
fscanf(file, "%s %lf",number,&data);
	D_priors->eta_nu=data;    
fscanf(file, "%s %lf",number,&data);
          D_priors->psi_nu=data;           

	/*K*//*r*//*nu*//*P*/
fscanf(file, "%s %lf",number,&data);
	D_priors->K_mu=data;   
fscanf(file, "%s %lf",number,&data);
   D_priors->eta_K_p=data;     
fscanf(file, "%s %lf",number,&data);
	D_priors->r_mu=data;       
fscanf(file, "%s %lf",number,&data);
     D_priors->eta_r_p=data;      
fscanf(file, "%s %lf",number,&data);
	D_priors->nu_mu=data;           
fscanf(file, "%s %lf",number,&data);
 D_priors->eta_nu_p=data;    
fscanf(file, "%s %lf",number,&data);
	D_priors->P_mu=data;     
fscanf(file, "%s %lf",number,&data); 
   D_priors->eta_P=data;   
	/*data2.c*/       
fscanf(file, "%s %lf",number,&data);
	D_priors->alpha_mu=data;     
fscanf(file, "%s %lf",number,&data);
     D_priors->eta_alpha=data;
fscanf(file, "%s %lf",number,&data);
	D_priors->beta_mu=data;        
fscanf(file, "%s %lf",number,&data);
   D_priors->eta_beta=data;
fscanf(file, "%s %lf",number,&data);
	D_priors->p=data;    
fscanf(file, "%s %lf",number,&data);
	D_priors->eta_gamma=data;
fscanf(file, "%s %lf",number,&data);
	D_priors->psi_gamma=data;
fscanf(file, "%s %lf",number,&data);
	D_priors->eta_omega=data;	
fscanf(file, "%s %lf",number,&data);
D_priors->psi_omega=data;
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
