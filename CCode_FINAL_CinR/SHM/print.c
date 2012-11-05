#include "headers.h"
#include "print.h"


void printheader(struct struct_data *D,char **HEADER){
  int l,m,mm,t=0;

for (l=0;l<D->L;l++){
	for (m=0;m<D->NoORF[l];m++){
		mm=D->NoSUM[l]+m;
sprintf(HEADER[t],"K_lm[%i]",mm);t+=1;
	}
}

for (l=0;l<D->L;l++){
sprintf(HEADER[t],"tau_K_l[%i]",l);t+=1;
}

for (l=0;l<D->L;l++){
sprintf(HEADER[t],"K_o_l[%i]",l);t+=1;
}

sprintf(HEADER[t],"sigma_K_o");t+=1;
sprintf(HEADER[t],"K_p");t+=1;
sprintf(HEADER[t],"P");t+=1;

for (l=0;l<D->L;l++){
	for (m=0;m<D->NoORF[l];m++){
		mm=D->NoSUM[l]+m;
		sprintf(HEADER[t],"r_lm[%i]",mm);t+=1;

	}
}


for (l=0;l<D->L;l++){
	sprintf(HEADER[t],"tau_r_l[%i]",l);t+=1;
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
sprintf(HEADER[t],"tau_K_p");t+=1;
sprintf(HEADER[t],"sigma_tau_K");t+=1;
sprintf(HEADER[t],"tau_r_p");t+=1;
sprintf(HEADER[t],"sigma_tau_r");t+=1;
}

void printdata(struct struct_data *D,struct struct_para *D_para,struct struct_MH *D_MH,double *OUT,int *T){
	int l,m,mm;
	for (l=0;l<D->L;l++){
		for (m=0;m<D->NoORF[l];m++){
			mm=D->NoSUM[l]+m;
OUT[*T]=D_para->K_lm[mm];*T+=1;
		}
	}

  for (l=0;l<D->L;l++)            { OUT[*T]=D_para->tau_K_l[l];*T+=1;}

  for (l=0;l<D->L;l++)            { OUT[*T]=D_para->K_o_l[l];*T+=1;}

 OUT[*T]=D_para->sigma_K_o;*T+=1;

  OUT[*T]=D_para->K_p;*T+=1;

 OUT[*T]=D_para->P;*T+=1;

  for (l=0;l<D->L;l++){
	for (m=0;m<D->NoORF[l];m++){
		mm=D->NoSUM[l]+m;
		OUT[*T]=D_para->r_lm[mm];*T+=1;	
	}
}
 
  for (l=0;l<D->L;l++)            { OUT[*T]=D_para->tau_r_l[l];*T+=1;}

  for (l=0;l<D->L;l++)            { OUT[*T]=D_para->r_o_l[l];*T+=1;}

  OUT[*T]=D_para->sigma_r_o;*T+=1;

  OUT[*T]=D_para->r_p;*T+=1;

  for (l=0;l<D->L;l++){ OUT[*T]=D_para->nu_l[l];*T+=1;}

OUT[*T]=D_para->sigma_nu;*T+=1;
OUT[*T]=D_para->nu_p;*T+=1;
OUT[*T]=D_para->tau_K_p;*T+=1;
OUT[*T]=D_para->sigma_tau_K;*T+=1;
OUT[*T]=D_para->tau_r_p;*T+=1;
OUT[*T]=D_para->sigma_tau_r;*T+=1;

}



