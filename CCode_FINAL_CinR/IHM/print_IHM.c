#include "headers_IHM.h"
#include "print_IHM.h"

void printheader_IHM(struct struct_data_IHM *D,char **HEADER){
  int l,t=0;
	for (l=0;l<D->L;l++){
	sprintf(HEADER[t],"Z_l[%i]",l);t+=1;
	}

	sprintf(HEADER[t],"sigma_Z");t+=1;
	sprintf(HEADER[t],"Z_p");t+=1;

	for (l=0;l<D->L;l++){
		sprintf(HEADER[t],"nu_l[%i]",l);t+=1;
	}

	sprintf(HEADER[t],"sigma_nu");t+=1;

	sprintf(HEADER[t],"nu_p");t+=1;

	sprintf(HEADER[t],"alpha_c");t+=1;

	for (l=0;l<D->L;l++){
		sprintf(HEADER[t],"delta_l[%i]",l);t+=1;
	}

	for (l=0;l<D->L;l++){
		sprintf(HEADER[t],"gamma_cl[%i]",l);t+=1;
	}

	sprintf(HEADER[t],"sigma_gamma");t+=1;

}

void printdata_IHM(struct struct_data_IHM *D,struct struct_para_IHM *D_para,struct struct_MH_IHM *D_MH,double *OUT,int *T){
  int l;
	for (l=0;l<D->L;l++){
		OUT[*T]=D_para->Z_l[l];*T+=1;
	}
	OUT[*T]=D_para->sigma_Z;*T+=1;
	OUT[*T]=D_para->Z_p;*T+=1;

	for (l=0;l<D->L;l++){
	OUT[*T]=D_para->nu_l[l];*T+=1;
	}

	OUT[*T]=D_para->sigma_nu;*T+=1;

	OUT[*T]=D_para->nu_p;*T+=1;

	OUT[*T]=D_para->alpha_c[1];*T+=1;

	for (l=0;l<D->L;l++){
		OUT[*T]=D_para->delta_l[l];*T+=1;
	}

	for (l=0;l<D->L;l++){
		OUT[*T]=D_para->gamma_cl[l];*T+=1;
	}

	OUT[*T]=D_para->sigma_gamma;*T+=1;
}



