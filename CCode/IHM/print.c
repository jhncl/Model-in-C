#include "headers.h"
#include "print.h"


void printheader(struct struct_data *D){
  int l;
	for (l=0;l<D->L;l++){
		printf ("Z_l[%i] ",l);
	}
	printf ("sigma_z ");
	printf ("Z_p ");


	for (l=0;l<D->L;l++){
		printf ("nu_l[%i] ",l);
	}

	printf ("sigma_nu ");

	printf ("nu_p ");

	printf("alpha ");

	printf("beta ");

	for (l=0;l<D->L;l++){
		printf ("delta[%i] ",l);
	}

	for (l=0;l<D->L;l++){
		printf ("gamma[%i] ",l);
	}

	printf("sigma_gamma ");

	for (l=0;l<D->L;l++){
		printf ("omega[%i] ",l);
	}

	printf("sigma_omega ");

	printf("upsilon_c[0] ");
	printf("upsilon_c[1] ");

	printf("sigma_upsilon ");

	printf("upsilon_p ");
	printf ("\n");
}

void printdata(struct struct_data *D,struct struct_para *D_para,struct struct_MH *D_MH){
  int l;
	for (l=0;l<D->L;l++){
		printf ("%g ",D_para->Z_l[l]);
	}
	printf ("%g ",D_para->sigma_Z);
	printf ("%g ",D_para->Z_p);


	for (l=0;l<D->L;l++){
		printf ("%g ",D_para->nu_l[l]);
	}

	printf ("%g ",D_para->sigma_nu);

	printf ("%g ",D_para->nu_p);

	printf ("%g ",D_para->alpha[1]);

	printf ("%g ",D_para->beta[1]);

	for (l=0;l<D->L;l++){
		printf ("%g ",D_para->delta[l]);
	}

	for (l=0;l<D->L;l++){
		printf ("%g ",D_para->gamma[l]);
	}

	printf ("%g ",D_para->sigma_gamma);

	for (l=0;l<D->L;l++){
		printf ("%g ",D_para->omega[l]);
	}

	printf ("%g ",D_para->sigma_omega);

	printf ("%g ",D_para->upsilon_c[0]);
	printf ("%g ",D_para->upsilon_c[1]);

	printf ("%g ",D_para->sigma_upsilon);

	printf ("%g ",D_para->upsilon_p);  

	printf ("%g ",D_MH->accept_Z);
	printf ("%g ",D_MH->accept_up);
	printf ("%g ",D_MH->accept_nu);
	printf ("\n");
  printf ("\n");
}



