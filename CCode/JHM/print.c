#include "headers.h"
#include "print.h"

void printheader(struct struct_data *D){
	int c,l,m,ll,mm;
	for (c=0;c<2;c++){
		for (l=0;l<D->L;l++){
			ll=c*D->L+l;
			for (m=0;m<D->NoORF[l];m++){ 
				mm=D->NoSUM[ll]+m;
				printf ("K_clm[%i] ",mm);
			}
		}
	}

	for (l=0;l<2*D->L;l++){
		printf ("tau_K_cl[%i] ",l);
	}

	for (l=0;l<D->L;l++){
		printf ("K_o_l[%i] ",l);
	}

	printf ("sigma_K_o ");

	printf ("K_p ");

	printf ("P ");

	for (c=0;c<2;c++){
		for (l=0;l<D->L;l++){
			ll=c*D->L+l;
			for (m=0;m<D->NoORF[l];m++){ 
				mm=D->NoSUM[ll]+m;
				printf ("r_clm[%i] ",mm);
			}
		}
	}

	for (l=0;l<2*D->L;l++){
		printf ("tau_r_cl[%i] ",l);
	}

	for (l=0;l<D->L;l++){
		printf ("r_o_l[%i] ",l);
	}

	printf ("sigma_r_o ");

	printf ("r_p ");


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
  


	printf ("accept_K ");
	printf ("accept_r ");
	printf ("accept_nu ");
	printf ("accept_P");
	printf ("\n");
}

void printdata(struct struct_data *D,struct struct_para *D_para,struct struct_MH *D_MH){
  int c,l,m,ll,mm;

	for (c=0;c<2;c++){
		for (l=0;l<D->L;l++){
			ll=c*D->L+l;
			for (m=0;m<D->NoORF[l];m++){ 
				mm=D->NoSUM[ll]+m;
				printf ("%g ",D_para->K_clm[mm]);
			}
		}
	}

	for (l=0;l<2*D->L;l++){
		printf ("%g ",D_para->tau_K_cl[l]);
	}

	for (l=0;l<D->L;l++){
		printf ("%g ",D_para->K_o_l[l]);
	}

	printf ("%g ",D_para->sigma_K_o);

	printf ("%g ",D_para->K_p);

	printf ("%g ",D_para->P);


	for (c=0;c<2;c++){
		for (l=0;l<D->L;l++){
			ll=c*D->L+l;
			for (m=0;m<D->NoORF[l];m++){ 
				mm=D->NoSUM[ll]+m;
				printf ("%g ",D_para->r_clm[mm]);
			}
		}
	}

	for (l=0;l<2*D->L;l++){
		printf ("%g ",D_para->tau_r_cl[l]);
	}

	for (l=0;l<D->L;l++){
		printf ("%g ",D_para->r_o_l[l]);
	}

	printf ("%g ",D_para->sigma_r_o);

	printf ("%g ",D_para->r_p);

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

	printf ("%g ",D_MH->accept_K);
	printf ("%g ",D_MH->accept_r);
	printf ("%g ",D_MH->accept_nu);
	printf ("%g",D_MH->accept_P);
	printf ("\n");
}
