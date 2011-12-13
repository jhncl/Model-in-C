#include "headers.h"
#include "print.h"


void printheader(struct struct_data *D){
  int l,m,mm;

for (l=0;l<D->L;l++){
	for (m=0;m<D->NoORF[l];m++){
		mm=D->NoSUM[l]+m;printf ("K_lm[%i] ",mm);
	}
}

for (l=0;l<D->L;l++){
	printf ("tau_K_l[%i] ",l);
}

for (l=0;l<D->L;l++){
	printf ("K_o_l[%i] ",l);
}

printf ("sigma_K_o ");

printf ("K_p ");

printf ("P ");

for (l=0;l<D->L;l++){
	for (m=0;m<D->NoORF[l];m++){
		mm=D->NoSUM[l]+m;
		printf ("r_lm[%i] ",mm);
	}
}

for (l=0;l<D->L;l++){
	printf ("tau_r_l[%i] ",l);
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

printf ("accept_K ");
printf ("accept_r ");
printf ("accept_nu ");
printf ("accept_P");
printf ("\n");
}

void printdata(struct struct_data *D,struct struct_para *D_para,struct struct_MH *D_MH){
	int l,m,mm;
	for (l=0;l<D->L;l++){
		for (m=0;m<D->NoORF[l];m++){
			mm=D->NoSUM[l]+m;
			printf ("%g ",D_para->K_lm[mm]);
		}
	}

  for (l=0;l<D->L;l++)            {printf ("%g ",D_para->tau_K_l[l]);}

  for (l=0;l<D->L;l++)            {printf ("%g ",D_para->K_o_l[l]);}

  printf ("%g ",D_para->sigma_K_o);

  printf ("%g ",D_para->K_p);

  printf ("%g ",D_para->P);

  for (l=0;l<D->L;l++){for (m=0;m<D->NoORF[l];m++){mm=D->NoSUM[l]+m;printf ("%g ",D_para->r_lm[mm]);}}
 
  for (l=0;l<D->L;l++)            {printf ("%g ",D_para->tau_r_l[l]);}

  for (l=0;l<D->L;l++)            {printf ("%g ",D_para->r_o_l[l]);}

  printf ("%g ",D_para->sigma_r_o);

  printf ("%g ",D_para->r_p);

  for (l=0;l<D->L;l++)            {printf ("%g ",D_para->nu_l[l]);}

  printf ("%g ",D_para->sigma_nu);

  printf ("%g ",D_para->nu_p);

  printf ("%g ",D_MH->accept_K);
  printf ("%g ",D_MH->accept_r);
  printf ("%g ",D_MH->accept_nu);
  printf ("%g",D_MH->accept_P);
  printf ("\n");
}



