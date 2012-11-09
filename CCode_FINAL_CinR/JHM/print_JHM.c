#include "headers_JHM.h"
#include "print_JHM.h"

void printheader_JHM(struct struct_data_JHM *D,char **HEADER){
	int c,l,m,ll,mm,t=0;
	for (c=0;c<2;c++){
		for (l=0;l<D->L;l++){
			ll=c*D->L+l;
			for (m=0;m<D->NoORF[l];m++){ 
				mm=D->NoSUM[ll]+m;
				sprintf(HEADER[t],"K_lm[%i]",mm);t+=1;
			}
		}
	}

	for (l=0;l<2*D->L;l++){
		sprintf(HEADER[t],"tau_K_cl[%i]",l);t+=1;
	}

	for (l=0;l<D->L;l++){
		sprintf(HEADER[t],"K_o_l[%i]",l);t+=1;
	}

sprintf(HEADER[t],"sigma_K_o");t+=1;

sprintf(HEADER[t],"K_p");t+=1;

sprintf(HEADER[t],"P");t+=1;

	for (c=0;c<2;c++){
		for (l=0;l<D->L;l++){
			ll=c*D->L+l;
			for (m=0;m<D->NoORF[l];m++){ 
				mm=D->NoSUM[ll]+m;
				sprintf(HEADER[t],"r_clm[%i]",mm);t+=1;
			}
		}
	}

	for (l=0;l<2*D->L;l++){
		sprintf(HEADER[t],"tau_r_cl[%i]",l);t+=1;
	}

	for (l=0;l<D->L;l++){
		sprintf(HEADER[t],"r_o_l[%i]",l);t+=1;
	}

sprintf(HEADER[t],"sigma_r_o");t+=1;

sprintf(HEADER[t],"r_p");t+=1;


	for (l=0;l<D->L;l++){
sprintf(HEADER[t],"nu_l[%i]",l);t+=1;
	}

sprintf(HEADER[t],"sigma_nu");t+=1;

sprintf(HEADER[t],"nu_p");t+=1;

sprintf(HEADER[t],"alpha_c");t+=1;

sprintf(HEADER[t],"beta_c");t+=1;

	for (l=0;l<D->L;l++){
sprintf(HEADER[t],"delta_l[%i]",l);t+=1;
	}

	for (l=0;l<D->L;l++){
sprintf(HEADER[t],"gamma_cl[%i]",l);t+=1;
	}

sprintf(HEADER[t],"sigma_gamma");t+=1;


	for (l=0;l<D->L;l++){
sprintf(HEADER[t],"omega_cl[%i]",l);t+=1;
	}

sprintf(HEADER[t],"sigma_omega");t+=1;

	for (c=0;c<2;c++){
sprintf(HEADER[t],"tau_K_p[%i]",c);t+=1;
sprintf(HEADER[t],"sigma_tau_K[%i]",c);t+=1;
	}
        for (c=0;c<2;c++){
sprintf(HEADER[t],"tau_r_p[%i]",c);t+=1;
sprintf(HEADER[t],"sigma_tau_r[%i]",c);t+=1;
	}
}

void printdata_JHM(struct struct_data_JHM *D,struct struct_para_JHM *D_para,struct struct_MH_JHM *D_MH,double *OUT,int *T){
  int c,l,m,ll,mm;

	for (c=0;c<2;c++){
		for (l=0;l<D->L;l++){
			ll=c*D->L+l;
			for (m=0;m<D->NoORF[l];m++){ 
				mm=D->NoSUM[ll]+m;
OUT[*T]=D_para->K_clm[mm];*T+=1;
			}
		}
	}

	for (l=0;l<2*D->L;l++){
OUT[*T]=D_para->tau_K_cl[l];*T+=1;
	}

	for (l=0;l<D->L;l++){
OUT[*T]=D_para->K_o_l[l];*T+=1;
	}

OUT[*T]=D_para->sigma_K_o;*T+=1;

OUT[*T]=D_para->K_p;*T+=1;

OUT[*T]=D_para->P;*T+=1;


	for (c=0;c<2;c++){
		for (l=0;l<D->L;l++){
			ll=c*D->L+l;
			for (m=0;m<D->NoORF[l];m++){ 
				mm=D->NoSUM[ll]+m;
OUT[*T]=D_para->r_clm[mm];*T+=1;
			}
		}
	}

	for (l=0;l<2*D->L;l++){
OUT[*T]=D_para->tau_r_cl[l];*T+=1;
	}

	for (l=0;l<D->L;l++){
OUT[*T]=D_para->r_o_l[l];*T+=1;
	}

OUT[*T]=D_para->sigma_r_o;*T+=1;

OUT[*T]=D_para->r_p;*T+=1;

	for (l=0;l<D->L;l++){
OUT[*T]=D_para->nu_l[l];*T+=1;
	}

OUT[*T]=D_para->sigma_nu;*T+=1;

OUT[*T]=D_para->nu_p;*T+=1;

OUT[*T]=D_para->alpha_c[1];*T+=1;

OUT[*T]=D_para->beta_c[1];*T+=1;


	for (l=0;l<D->L;l++){
OUT[*T]=D_para->delta_l[l];*T+=1;
	}

	for (l=0;l<D->L;l++){
OUT[*T]=D_para->gamma_cl[l];*T+=1;
	}

OUT[*T]=D_para->sigma_gamma;*T+=1;

	for (l=0;l<D->L;l++){
		OUT[*T]=D_para->omega_cl[l];*T+=1;
	}

OUT[*T]=D_para->sigma_gamma;*T+=1;

	for (c=0;c<2;c++){
OUT[*T]=D_para->tau_K_p[c];*T+=1;
OUT[*T]=D_para->sigma_tau_r[c];*T+=1;
        }
        for (c=0;c<2;c++){
OUT[*T]=D_para->tau_r_p[c];*T+=1;
OUT[*T]=D_para->sigma_tau_r[c];*T+=1;
        }
}
