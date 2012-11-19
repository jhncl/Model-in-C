#include "headers_IHM.h"
#include "print_IHM.h"


void printheader_IHM(struct struct_data_IHM *D,char **HEADER){
  int l,t=0;
	for (l=0;l<D->L;l++){
		/*printf ("Z_l[%i] ",l);*/
sprintf(HEADER[t],"Z_l[%i]",l);t+=1;
	}
	/*printf ("sigma_Z ");*/sprintf(HEADER[t],"sigma_Z");t+=1;
	/*printf ("Z_p ");*/sprintf(HEADER[t],"Z_p");t+=1;


	for (l=0;l<D->L;l++){
		/*printf ("nu_l[%i] ",l);*/sprintf(HEADER[t],"nu_l[%i]",l);t+=1;
	}

	/*printf ("sigma_nu ");*/sprintf(HEADER[t],"sigma_nu");t+=1;

	/*printf ("nu_p ");*/sprintf(HEADER[t],"nu_p");t+=1;

	/*printf("alpha_c ");*/sprintf(HEADER[t],"alpha_c");t+=1;

	for (l=0;l<D->L;l++){
		/*printf ("delta_l[%i] ",l);*/sprintf(HEADER[t],"delta_l[%i]",l);t+=1;
	}

	for (l=0;l<D->L;l++){
		/*printf ("gamma_cl[%i] ",l);*/sprintf(HEADER[t],"gamma_cl[%i]",l);t+=1;
	}

	/*printf("sigma_gamma ");*/sprintf(HEADER[t],"sigma_gamma");t+=1;

	/*printf("upsilon_c ");*/sprintf(HEADER[t],"upsilon_c");t+=1;

	/*printf("sigma_upsilon ");*/sprintf(HEADER[t],"sigma_upsilon");t+=1;
/*
	printf("accept_Z ");
	printf("accept_up ");
	printf("accept_nu ");

	printf ("\n");*/
}

void printdata_IHM(struct struct_data_IHM *D,struct struct_para_IHM *D_para,struct struct_MH_IHM *D_MH,double *OUT,int *T){
  int l;
	for (l=0;l<D->L;l++){
		/*printf ("%g ",D_para->Z_l[l]);*/ OUT[*T]=D_para->Z_l[l];*T+=1;
	}
	/*printf ("%g ",D_para->sigma_Z);*/ OUT[*T]=D_para->sigma_Z;*T+=1;
	/*printf ("%g ",D_para->Z_p);*/ OUT[*T]=D_para->Z_p;*T+=1;


	for (l=0;l<D->L;l++){
		/*printf ("%g ",D_para->nu_l[l]);*/ OUT[*T]=D_para->nu_l[l];*T+=1;
	}

	/*printf ("%g ",D_para->sigma_nu);*/ OUT[*T]=D_para->sigma_nu;*T+=1;

	/*printf ("%g ",D_para->nu_p);*/ OUT[*T]=D_para->nu_p;*T+=1;

	/*printf ("%g ",D_para->alpha_c[1]);*/ OUT[*T]=D_para->alpha_c[1];*T+=1;

	for (l=0;l<D->L;l++){
		/*printf ("%g ",D_para->delta_l[l]);*/ OUT[*T]=D_para->delta_l[l];*T+=1;
	}

	for (l=0;l<D->L;l++){
		/*printf ("%g ",D_para->gamma_cl[l]);*/ OUT[*T]=D_para->gamma_cl[l];*T+=1;
	}

	/*printf ("%g ",D_para->sigma_gamma);*/ OUT[*T]=D_para->sigma_gamma;*T+=1;

	/*printf ("%g ",D_para->upsilon_c[1]);*/ OUT[*T]=D_para->upsilon_c[1];*T+=1;

	/*printf ("%g ",D_para->sigma_upsilon);*/  OUT[*T]=D_para->sigma_upsilon;*T+=1;

/*
	printf ("%g ",D_MH->accept_Z);
	printf ("%g ",D_MH->accept_up);
	printf ("%g ",D_MH->accept_nu);
	printf ("\n");*/
}



