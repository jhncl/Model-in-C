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

gsl_matrix *getParticleMatrix(char *filename)
{
    int line_length = 10000000;    
    int nrows, ncols;
       
    FILE* f;  
    char *pch;
    char line[10000000];
    gsl_matrix *particles;
    
    
    f=fopen(filename, "r");
    if(NULL==f)
    {
        fprintf(stderr, "Cannot open file %s\n", filename);
        exit(1);
    }

    nrows = 0, ncols = 0;
    /*Scan once to get the dimensions
        there doesn't seem to be a realloc matrix function
    */

    while(fgets(line, line_length, f) != NULL){
        pch = strtok(line," ");
        while(nrows == 0 && pch != NULL )
        {
            ncols++;
            pch = strtok(NULL," ");
        }
        nrows++;
    }
    
    fclose(f);
        
    /*Create matrix and fill up*/
    particles = gsl_matrix_alloc(nrows, ncols);
    nrows = 0; ncols = 0;
    f=fopen(filename, "r");
      
    while(fgets(line, line_length, f) != NULL){
        pch = strtok(line," ");
        while(pch != NULL )
        {
            gsl_matrix_set(particles, nrows, ncols, atof(pch));
            ncols++;
            pch = strtok(NULL," ");
        }
        ncols = 0;
        nrows++;
    }
    fclose(f);
    
    return(particles);    
}    

double meanParticleMatrix(gsl_matrix* matrix,int col,int start,int end){
	int i;
	double mean;
	for (i=start;i<=end;i++){
		mean=mean+gsl_matrix_get(matrix,i,col);
	}
	mean=mean/(end-start+1);
	return(mean);
}

int datadouble(char filename[], char filename2[], double datavec[],int length,struct_data *D)
{
int i;
double P,K,r;
gsl_matrix* matrixA=getParticleMatrix(filename);
gsl_matrix* matrixB=getParticleMatrix(filename2);
i=D->SHIFTmn+2*D->L+3-1;
P=exp(meanParticleMatrix(matrixA,i,1,2));
for (i=0;i<D->SHIFTmn;i++){
K=exp(meanParticleMatrix(matrixA,i,1,2));
r=exp(meanParticleMatrix(matrixA,i+D->SHIFTmn+2*D->L+4-1,1,2));
D->y[i]=(r/log(2*(K-P)/(K-2*P)))*(log(K/P)/log(2));
}
i=D->SHIFTmn+2*D->L+3-1;
P=exp(meanParticleMatrix(matrixB,i,1,2));
for (i=0;i<(D->MAXmn-D->SHIFTmn);i++){
K=exp(meanParticleMatrix(matrixB,i,1,2));
r=exp(meanParticleMatrix(matrixB,i+(D->MAXmn-D->SHIFTmn)+2*D->L+4-1,1,2));
D->y[D->SHIFTmn+i]=(r/log(2*(K-P)/(K-2*P)))*(log(K/P)/log(2));
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
	datadouble("dd.txt","ee.txt",data->y,size,data);
	/*if (data->y==NULL||data->x==NULL||data->NoORF==NULL||data->NoSUM==NULL||data->NoTIME==NULL) {
		perror("malloc failed");
    		exit(EXIT_FAILURE);
  	};*/
return 0;
}

int inzstruct_para(struct_para *para,struct_data *data)
{
	long size;
	size=data->L;
	para->delta=malloc(size*sizeof(double));
	para->gamma=malloc(size*sizeof(double));
	para->Z_l=malloc(size*sizeof(double));
	para->nu_l=malloc(size*sizeof(double));
	size=2;
	para->alpha=malloc(size*sizeof(double));
	para->upsilon_c=malloc(size*sizeof(double));
	fillpara(para,data);
return 0;
}

/*FILL*/

int fillMH(struct_MH *MH)
{
	MH->hZ=0.1;	MH->accept_Z=0;
	MH->hup=0.1;	MH->accept_up=0;
	MH->hnu=0.1;	MH->accept_nu=0; /*h sd; accept=0*/
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

int fillpara(struct_para *D_para, struct_data *D)
{
	int l;

	/*initials*/
	for (l=0;l<D->L;l++){D_para->Z_l[l]=gsl_sf_log(50);}
	D_para->sigma_Z=1/100;     
	D_para->Z_p=1;  

	for (l=0;l<D->L;l++)          {D_para->nu_l[l]=18;}          
	D_para->sigma_nu=0.0025;   
  
	D_para->nu_p=18;    

	for (l=0;l<D->L;l++)          {D_para->gamma[l]=0;} 

	for (l=0;l<D->L;l++)          {D_para->delta[l]=1;}  

	D_para->alpha[0]=gsl_sf_log(1);
	D_para->alpha[1]=gsl_sf_log(1);
	D_para->sigma_gamma=1/4;
	D_para->upsilon_c[0]=1; 
	D_para->upsilon_c[1]=1;       D_para->sigma_upsilon=1;
	D_para->upsilon_p=1;

return 0;
}

int fillpriors(struct_priors *D_priors)
{
	/*Priors*/
	/*K*/


	D_priors->Z_mu=gsl_sf_log(50);   	   D_priors->eta_Z_mu=1;      /*Normal  LMean; Precisions */
	D_priors->eta_K=1;           		   D_priors->psi_Z=200000;               /*Gamma  Shape; Scale */

	/*nu*/
	D_priors->eta_nu_p=2;			   D_priors->psi_nu=200000;              /*Gamma  Shape; Scale */
	D_priors->nu_mu=gsl_sf_log(31);            D_priors->eta_nu_mu=1;     /*Normal  LMean; Precisions */
	D_priors->alpha_mu=gsl_sf_log(1);          D_priors->eta_alpha=1/3^2;
	D_priors->p=0.01;    
	D_priors->eta_gamma=1;           	   D_priors->psi_gamma=1;
	D_priors->eta_upsilon=1;	 	   D_priors->phi_upsilon=1;	    
	D_priors->upsilon_mu=1;			   D_priors->eta_upsilon_mu=1;
return 0;
}

/* eof  */
