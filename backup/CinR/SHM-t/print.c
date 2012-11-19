#include "headers.h"
#include "print.h"


void printheader(struct struct_data *D,char **HEADER){
  int l,m,mm,t=0;


for (l=0;l<D->L;l++){
	for (m=0;m<D->NoORF[l];m++){
		mm=D->NoSUM[l]+m;/*printf ("K_lm[%i] ",mm);*/
sprintf(HEADER[t],"K_lm[%i]",mm);t+=1;
	}
}

for (l=0;l<D->L;l++){
	/*printf ("tau_K_l[%i] ",l);*/
sprintf(HEADER[t],"tau_K_l[%i]",l);t+=1;
}

for (l=0;l<D->L;l++){
	/*printf ("K_o_l[%i] ",l);*/
sprintf(HEADER[t],"K_o_l[%i]",l);t+=1;
}

/*printf ("sigma_K_o ");*/sprintf(HEADER[t],"sigma_K_o");t+=1;
/*printf ("K_p ");*/sprintf(HEADER[t],"K_p");t+=1;
/*printf ("P ");*/sprintf(HEADER[t],"P");t+=1;
for (l=0;l<D->L;l++){
	for (m=0;m<D->NoORF[l];m++){
		mm=D->NoSUM[l]+m;
		/*printf ("r_lm[%i] ",mm);*/sprintf(HEADER[t],"r_lm[%i]",mm);t+=1;

	}
}

for (l=0;l<D->L;l++){
	/*printf ("tau_r_l[%i] ",l);*/sprintf(HEADER[t],"tau_r_l[%i]",l);t+=1;
}

for (l=0;l<D->L;l++){
	/*printf ("r_o_l[%i] ",l);*/sprintf(HEADER[t],"r_o_l[%i]",l);t+=1;
}

/*printf ("sigma_r_o ");*/sprintf(HEADER[t],"sigma_r_o");t+=1;

/*printf ("r_p ");*/sprintf(HEADER[t],"r_p");t+=1;

for (l=0;l<D->L;l++){
	/*printf ("nu_l[%i] ",l);*/sprintf(HEADER[t],"nu_l[%i]",l);t+=1;
}

/*printf ("sigma_nu ");*/sprintf(HEADER[t],"sigma_nu");t+=1;

/*printf ("nu_p ");*/sprintf(HEADER[t],"nu_p");t+=1;
 /*printf ("tau_K_p ");*/sprintf(HEADER[t],"tau_K_p");t+=1;
 /*printf ("sigma_tau_K ");*/sprintf(HEADER[t],"sigma_tau_K");t+=1;
 /*printf ("tau_r_p ");*/sprintf(HEADER[t],"tau_r_p");t+=1;
 /*printf ("sigma_tau_r ");*/sprintf(HEADER[t],"sigma_tau_r");t+=1;
 /*printf ("A ");
 printf ("B ");


printf ("accept_K ");
printf ("accept_r ");
printf ("accept_nu ");
printf ("accept_P");
printf ("\n");*/
}

void printdata(struct struct_data *D,struct struct_para *D_para,struct struct_MH *D_MH,double *OUT,int *T){
	int l,m,mm;
	for (l=0;l<D->L;l++){
		for (m=0;m<D->NoORF[l];m++){
			mm=D->NoSUM[l]+m;
			 /* printf ("%g ",D_para->K_lm[mm]);*/
OUT[*T]=D_para->K_lm[mm];*T+=1;
		}
	}

  for (l=0;l<D->L;l++)            { /* printf ("%g ",D_para->tau_K_l[l]);*/OUT[*T]=D_para->tau_K_l[l];*T+=1;}

  for (l=0;l<D->L;l++)            { /* printf ("%g ",D_para->K_o_l[l]);*/OUT[*T]=D_para->K_o_l[l];*T+=1;}

  /*  printf ("%g ",D_para->sigma_K_o);*/OUT[*T]=D_para->sigma_K_o;*T+=1;

  /*  printf ("%g ",D_para->K_p);*/OUT[*T]=D_para->K_p;*T+=1;

  /*  printf ("%g ",D_para->P);*/OUT[*T]=D_para->P;*T+=1;

  for (l=0;l<D->L;l++){
	for (m=0;m<D->NoORF[l];m++){
		mm=D->NoSUM[l]+m;
		 /* printf ("%g ",D_para->r_lm[mm]);*/OUT[*T]=D_para->r_lm[mm];*T+=1;	
	}
}
 
  for (l=0;l<D->L;l++)            { /* printf ("%g ",D_para->tau_r_l[l]);*/OUT[*T]=D_para->tau_r_l[l];*T+=1;}

  for (l=0;l<D->L;l++)            { /* printf ("%g ",D_para->r_o_l[l]);*/OUT[*T]=D_para->r_o_l[l];*T+=1;}

  /*  printf ("%g ",D_para->sigma_r_o);*/OUT[*T]=D_para->sigma_r_o;*T+=1;

  /*  printf ("%g ",D_para->r_p);*/OUT[*T]=D_para->r_p;*T+=1;

  for (l=0;l<D->L;l++){ /* printf ("%g ",D_para->nu_l[l]);*/OUT[*T]=D_para->nu_l[l];*T+=1;}

  /*  printf ("%g ",D_para->sigma_nu);*/OUT[*T]=D_para->sigma_nu;*T+=1;

  /*  printf ("%g ",D_para->nu_p);*/OUT[*T]=D_para->nu_p;*T+=1;

  /*  printf ("%g ",D_para->tau_K_p);*/OUT[*T]=D_para->tau_K_p;*T+=1;
  /*  printf ("%g ",D_para->sigma_tau_K);*/OUT[*T]=D_para->sigma_tau_K;*T+=1;
  /*  printf ("%g ",D_para->tau_r_p);*/OUT[*T]=D_para->tau_r_p;*T+=1;
  /*  printf ("%g ",D_para->sigma_tau_r);*/OUT[*T]=D_para->sigma_tau_r;*T+=1;

 /* printf ("%g ",D_para->A);
  printf ("%g ",D_para->B);


  printf ("%g ",D_MH->accept_K);
  printf ("%g ",D_MH->accept_r);
  printf ("%g ",D_MH->accept_nu);
  printf ("%g",D_MH->accept_P);
  printf ("\n");*/
}



