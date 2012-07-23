#include "headers_JHM.h"
#include "print_JHM.h"

void printheader_JHM(struct struct_data_JHM *D,char **HEADER){
	int c,l,m,ll,mm,t=0;
	for (c=0;c<2;c++){
		for (l=0;l<D->L;l++){
			ll=c*D->L+l;
			for (m=0;m<D->NoORF[l];m++){ 
				mm=D->NoSUM[ll]+m;
				/*printf ("K_clm[%i] ",mm);*/
sprintf(HEADER[t],"K_lm[%i]",mm);t+=1;
			}
		}
	}

	for (l=0;l<2*D->L;l++){
		/*printf ("tau_K_cl[%i] ",l);*/
sprintf(HEADER[t],"tau_K_cl[%i]",l);t+=1;
	}

	for (l=0;l<D->L;l++){
		/*printf ("K_o_l[%i] ",l);*/
sprintf(HEADER[t],"K_o_l[%i]",l);t+=1;
	}

	/*printf ("sigma_K_o ");*/sprintf(HEADER[t],"sigma_K_o");t+=1;

	/*printf ("K_p ");*/sprintf(HEADER[t],"K_p");t+=1;

	/*printf ("P ");*/sprintf(HEADER[t],"P");t+=1;

	for (c=0;c<2;c++){
		for (l=0;l<D->L;l++){
			ll=c*D->L+l;
			for (m=0;m<D->NoORF[l];m++){ 
				mm=D->NoSUM[ll]+m;
				/*printf ("r_clm[%i] ",mm);*/
sprintf(HEADER[t],"r_clm[%i]",mm);t+=1;
			}
		}
	}

	for (l=0;l<2*D->L;l++){
		/*printf ("tau_r_cl[%i] ",l);*/
sprintf(HEADER[t],"tau_r_cl[%i]",l);t+=1;
	}

	for (l=0;l<D->L;l++){
		/*printf ("r_o_l[%i] ",l);*/
sprintf(HEADER[t],"r_o_l[%i]",l);t+=1;
	}

	/*printf ("sigma_r_o ");*/sprintf(HEADER[t],"sigma_r_o");t+=1;

	/*printf ("r_p ");*/sprintf(HEADER[t],"r_p");t+=1;


	for (l=0;l<D->L;l++){
		/*printf ("nu_l[%i] ",l);*/
sprintf(HEADER[t],"nu_l[%i]",l);t+=1;
	}

	/*printf ("sigma_nu ");*/sprintf(HEADER[t],"sigma_nu");t+=1;

	/*printf ("nu_p ");*/sprintf(HEADER[t],"nu_p");t+=1;

	/*printf("alpha_c ");*/sprintf(HEADER[t],"alpha_c");t+=1;

	/*printf("beta_c ");*/sprintf(HEADER[t],"beta_c");t+=1;

	for (l=0;l<D->L;l++){
		/*printf ("delta_l[%i] ",l);*/
sprintf(HEADER[t],"delta_l[%i]",l);t+=1;
	}

	for (l=0;l<D->L;l++){
		/*printf ("gamma_cl[%i] ",l);*/
sprintf(HEADER[t],"gamma_cl[%i]",l);t+=1;
	}

	/*printf("sigma_gamma ");*/sprintf(HEADER[t],"sigma_gamma");t+=1;

	for (l=0;l<D->L;l++){
		/*printf ("omega_cl[%i] ",l);*/
sprintf(HEADER[t],"omega_cl[%i]",l);t+=1;
	}

	/*printf("sigma_omega ");*/sprintf(HEADER[t],"sigma_omega");t+=1;

	/*printf("upsilon_c ");*/sprintf(HEADER[t],"upsilon_c");t+=1;

	/*printf("sigma_upsilon ");*/sprintf(HEADER[t],"sigma_upsilon");t+=1;

	for (c=0;c<2;c++){
	/*printf ("tau_K_p[%i] ",c);*/
sprintf(HEADER[t],"tau_K_p[%i]",c);t+=1;
        /*printf ("sigma_tau_K[%i] ",c);*/
sprintf(HEADER[t],"sigma_tau_K[%i]",c);t+=1;
	}
        for (c=0;c<2;c++){
	  /*printf ("tau_r_p[%i] ",c);*/
sprintf(HEADER[t],"tau_r_p[%i]",c);t+=1;
	  /*printf ("sigma_tau_r[%i] ",c);*/
sprintf(HEADER[t],"sigma_tau_r[%i]",c);t+=1;
	}


/*
	printf ("A ");
        printf ("B ");


	printf ("A ");
	printf ("B ");
	printf ("accept_K ");
	printf ("accept_r ");
	printf ("accept_nu ");
	printf ("accept_P");
	printf ("\n");*/
}

void printdata_JHM(struct struct_data_JHM *D,struct struct_para_JHM *D_para,struct struct_MH_JHM *D_MH,double *OUT,int *T){
  int c,l,m,ll,mm;

	for (c=0;c<2;c++){
		for (l=0;l<D->L;l++){
			ll=c*D->L+l;
			for (m=0;m<D->NoORF[l];m++){ 
				mm=D->NoSUM[ll]+m;
				/*printf ("%g ",D_para->K_clm[mm]);*/
OUT[*T]=D_para->K_clm[mm];*T+=1;
			}
		}
	}

	for (l=0;l<2*D->L;l++){
		/*printf ("%g ",D_para->tau_K_cl[l]);*/
OUT[*T]=D_para->tau_K_cl[l];*T+=1;
	}

	for (l=0;l<D->L;l++){
		/*printf ("%g ",D_para->K_o_l[l]);*/
OUT[*T]=D_para->K_o_l[l];*T+=1;
	}

	/*printf ("%g ",D_para->sigma_K_o);*/OUT[*T]=D_para->sigma_K_o;*T+=1;

	/*printf ("%g ",D_para->K_p);*/OUT[*T]=D_para->K_p;*T+=1;

	/*printf ("%g ",D_para->P);*/OUT[*T]=D_para->P;*T+=1;


	for (c=0;c<2;c++){
		for (l=0;l<D->L;l++){
			ll=c*D->L+l;
			for (m=0;m<D->NoORF[l];m++){ 
				mm=D->NoSUM[ll]+m;
				/*printf ("%g ",D_para->r_clm[mm]);*/
OUT[*T]=D_para->r_clm[mm];*T+=1;
			}
		}
	}

	for (l=0;l<2*D->L;l++){
		/*printf ("%g ",D_para->tau_r_cl[l]);*/
OUT[*T]=D_para->tau_r_cl[l];*T+=1;
	}

	for (l=0;l<D->L;l++){
		/*printf ("%g ",D_para->r_o_l[l]);*/
OUT[*T]=D_para->r_o_l[l];*T+=1;
	}

	/*printf ("%g ",D_para->sigma_r_o);*/OUT[*T]=D_para->sigma_r_o;*T+=1;

	/*printf ("%g ",D_para->r_p);*/OUT[*T]=D_para->r_p;*T+=1;

	for (l=0;l<D->L;l++){
		/*printf ("%g ",D_para->nu_l[l]);*/
OUT[*T]=D_para->nu_l[l];*T+=1;
	}

	/*printf ("%g ",D_para->sigma_nu);*/OUT[*T]=D_para->sigma_nu;*T+=1;

	/*printf ("%g ",D_para->nu_p);*/OUT[*T]=D_para->nu_p;*T+=1;

	/*printf ("%g ",D_para->alpha_c[1]);*/OUT[*T]=D_para->alpha_c[1];*T+=1;

	/*printf ("%g ",D_para->beta_c[1]);*/OUT[*T]=D_para->beta_c[1];*T+=1;

	for (l=0;l<D->L;l++){
		/*printf ("%g ",D_para->delta_l[l]);*/
OUT[*T]=D_para->delta_l[l];*T+=1;
	}

	for (l=0;l<D->L;l++){
		/*printf ("%g ",D_para->gamma_cl[l]);*/
OUT[*T]=D_para->gamma_cl[l];*T+=1;
	}

	/*printf ("%g ",D_para->sigma_gamma);*/OUT[*T]=D_para->sigma_gamma;*T+=1;

	for (l=0;l<D->L;l++){
		/*printf ("%g ",D_para->omega_cl[l]);*/
OUT[*T]=D_para->omega_cl[l];*T+=1;
	}

	/*printf ("%g ",D_para->sigma_omega);*/OUT[*T]=D_para->sigma_omega;*T+=1;

	/*printf ("%g ",D_para->upsilon_c[1]);*/OUT[*T]=D_para->upsilon_c[1];*T+=1;

	/*printf ("%g ",D_para->sigma_upsilon);*/OUT[*T]=D_para->sigma_upsilon;*T+=1;

	for (c=0;c<2;c++){
	  /*printf ("%g ", D_para->tau_K_p[c]);*/
OUT[*T]=D_para->tau_K_p[c];*T+=1;
	  /*printf ("%g ",D_para->sigma_tau_K[c]);*/
OUT[*T]=D_para->sigma_tau_K[c];*T+=1;
        }
        for (c=0;c<2;c++){
          /*printf ("%g ",D_para->tau_r_p[c]);*/
OUT[*T]=D_para->tau_r_p[c];*T+=1;
          /*printf ("%g ",D_para->sigma_tau_r[c]);*/
OUT[*T]=D_para->sigma_tau_r[c];*T+=1;
        }
 /*
	printf ("%g ",D_para->A);
	printf ("%g ",D_para->B);

	printf ("%g ",D_MH->accept_K);
	printf ("%g ",D_MH->accept_r);
	printf ("%g ",D_MH->accept_nu);
	printf ("%g",D_MH->accept_P);
	printf ("\n");*/
}
