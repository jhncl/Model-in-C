#include "headers.h"
#include "print.h"


void printheader(struct struct_data *D){
  int l;
	for (l=0;l<D->L;l++){
		printf ("Z_l[%i] ",l);
	}
	printf ("sigma_Z ");
	printf ("Z_p ");


	for (l=0;l<(2*D->L);l++){
		printf ("nu_cl[%i] ",l);
	}

	printf ("sigma_nu ");

	printf ("nu_p ");

	printf("alpha_c ");

	for (l=0;l<D->L;l++){
		printf ("delta_l[%i] ",l);
	}

	for (l=0;l<D->L;l++){
		printf ("gamma_cl[%i] ",l);
	}

	printf("sigma_gamma ");

	printf("accept_Z ");
	printf("accept_up ");
	printf("accept_nu ");

	printf ("\n");
}

void printdata(struct struct_data *D,struct struct_para *D_para,struct struct_MH *D_MH){
  int l;
	for (l=0;l<D->L;l++){
		printf ("%g ",D_para->Z_l[l]);
	}
	printf ("%g ",D_para->sigma_Z);
	printf ("%g ",D_para->Z_p);


	for (l=0;l<(2*D->L);l++){
		printf ("%g ",D_para->nu_l[l]);
	}

	printf ("%g ",D_para->sigma_nu);

	printf ("%g ",D_para->nu_p);

	printf ("%g ",D_para->alpha_c[1]);

	for (l=0;l<D->L;l++){
		printf ("%g ",D_para->delta_l[l]);
	}

	for (l=0;l<D->L;l++){
		printf ("%g ",D_para->gamma_cl[l]);
	}

	printf ("%g ",D_para->sigma_gamma);

	printf ("%g ",D_MH->accept_Z);
	printf ("%g ",D_MH->accept_up);
	printf ("%g ",D_MH->accept_nu);
	printf ("\n");
}



