#include "headers_JHM.h"
#include "functions_JHM.h"
#include "print_JHM.h"
/*Logistic Growth*/

double logistic_function_JHM(double t,double K, double r, double P){
	double output;
	output=(K*P*exp(r*t))/(K+P*(exp(r*t)-1));
	return(output);
}

double logistic_function_E_JHM(double t,double K, double r, double P){
	double output;
	output=logistic_function_JHM(t,exp(K), exp(r), exp(P));
	return(output);
}

/*Gibbs*/

double gauss_sample_JHM(struct_data_JHM *D ,int start, int N,double x[],double tau,double mu_0,double tau_0){
	double vec,Ndouble,SUM=0;
	int i;
	Ndouble=N;
	for (i=start;i<(start+N);i++){SUM=SUM+x[i];}
	vec=rnorm((tau_0*mu_0+tau*SUM)/(tau_0+Ndouble*tau),1/sqrt(tau_0+Ndouble*tau));
	return(vec);
}

/*MCMC MH*/

double MCMC_base_JHM(struct_data_JHM *D,struct_para_JHM *D_para,struct_priors_JHM *D_priors,double *accept,double *h,double para,double (*foo)(struct struct_data_JHM *D,struct struct_para_JHM *D_para,struct struct_priors_JHM *D_priors,double,int,int,int),int c,int l, int m){
	double logu,logaprob,can;
	can=rnorm(para,*h);
	logaprob=(*foo)(D,D_para,D_priors,can,c,l,m)-(*foo)(D,D_para,D_priors,para,c,l,m);

  logu=log(1-runif(0,1));
	if (logaprob>logu){para=can;*accept=*accept+1;}
	return(para); 
	}

double MCMC_base_truncate_low_JHM(double truncate,struct_data_JHM *D,struct_para_JHM *D_para,struct_priors_JHM *D_priors,double *accept,double *h,double para,double (*foo)(struct struct_data_JHM *D,struct struct_para_JHM *D_para,struct struct_priors_JHM *D_priors,double,int,int,int),int c,int l, int m){
  double logu,logaprob,can;

	can=rnorm(para,*h);
    if(can<(truncate)){
      can=para;
    }
 
  logaprob=(*foo)(D,D_para,D_priors,can,c,l,m)-(*foo)(D,D_para,D_priors,para,c,l,m);

  logu=log(1-runif(0,1));
  if (logaprob>logu){para=can;*accept=*accept+1;}
  return(para);
}
double MCMC_base_truncate_high_JHM(double truncate,struct_data_JHM *D,struct_para_JHM *D_para,struct_priors_JHM *D_priors,double *accept,double *h,double para,double (*foo)(struct struct_data_JHM *D,struct struct_para_JHM *D_para,struct struct_priors_JHM *D_priors,double,int,int,int),int c,int l, int m){
  double logu,logaprob,can;
 
	can=rnorm(para,*h);
    if(can>(truncate)){
      can=para;
    }
  logaprob=(*foo)(D,D_para,D_priors,can,c,l,m)-(*foo)(D,D_para,D_priors,para,c,l,m);

  logu=log(1-runif(0,1));
  if (logaprob>logu){para=can;*accept=*accept+1;}
  return(para);
}



double MCMC_K_clm_JHM(struct_data_JHM *D,struct_para_JHM *D_para,struct_priors_JHM *D_priors,double para,int c,int l, int m){
	double density,F,SUM=0;
	int n,mm,nn,ll=c*D->L+l;
	mm=D->NoSUM[ll]+m;
	for (n=0;n<D->NoTIME[mm];n++){
		nn=c*D->SHIFTlmn+l*D->M*D->N + m*D->N+n;
		F=D->y[nn]-logistic_function_E_JHM(D->x[nn],para, D_para->r_clm[mm],D_para->P);
		SUM=(1+D_priors->df2)*log(1+pow(F,2)*exp(D_para->nu_l[l]+D_para->upsilon_c[c])/D_priors->df2)+SUM;
/*F*F*exp(D_para->nu_l[l]+D_para->upsilon_c[c])+SUM;*/
	}	
	F=para-(D_para->alpha_c[c]+D_para->K_o_l[l]+c*D_para->delta_l[l]*D_para->gamma_cl[l]);
	density=F*F*exp(D_para->tau_K_cl[ll])+SUM; 
	return(-0.5*density);
}

double MCMC_r_clm_JHM(struct_data_JHM *D,struct_para_JHM *D_para,struct_priors_JHM *D_priors,double para,int c,int l, int m){
	double density,F,SUM=0;
	int n,mm,nn,ll=c*D->L+l;
	mm=D->NoSUM[ll]+m;
	for (n=0;n<D->NoTIME[mm];n++){
		nn=c*D->SHIFTlmn+l*D->M*D->N + m*D->N+n;
		F=D->y[nn]-logistic_function_E_JHM(D->x[nn], D_para->K_clm[mm],para,D_para->P);	
		SUM=(1+D_priors->df2)*log(1+pow(F,2)*exp(D_para->nu_l[l]+D_para->upsilon_c[c])/D_priors->df2)+SUM;
/*F*F*exp(D_para->nu_l[l]+D_para->upsilon_c[c])+SUM;*/
	}	
	F=para-(D_para->beta_c[c]+D_para->r_o_l[l]+c*D_para->delta_l[l]*D_para->omega_cl[l]);
	density=F*F*exp(D_para->tau_r_cl[ll])+SUM; 
	return(-0.5*density); 
}


double MCMC_tau_K_cl_JHM(struct_data_JHM *D,struct_para_JHM *D_para,struct_priors_JHM *D_priors,double para,int c,int l, int m){
	double density,SUM=0,F;
	int mm,ll=c*D->L+l;
	for (m=0;m<D->NoORF[ll];m++){
		mm=D->NoSUM[ll]+m;
		F=D_para->K_clm[mm]-(D_para->alpha_c[c]+D_para->K_o_l[l]+c*D_para->delta_l[l]*D_para->gamma_cl[l]);
		SUM=-para+F*F*exp(para)+SUM;
	}	
	F=para-D_para->tau_K_p[c]-exp(D_para->alpha_c[c]+D_para->K_o_l[l])*D_para->A;
	density=F*F*exp(D_para->sigma_tau_K[c])+SUM; 
	return(-0.5*density); 
}

double MCMC_tau_r_cl_JHM(struct_data_JHM *D,struct_para_JHM *D_para,struct_priors_JHM *D_priors,double para,int c,int l, int m){
	double density,SUM=0,F;
	int mm,ll=c*D->L+l;
	for (m=0;m<D->NoORF[ll];m++){
		mm=D->NoSUM[ll]+m;
		F=D_para->r_clm[mm]-(D_para->beta_c[c]+D_para->r_o_l[l]+c*D_para->delta_l[l]*D_para->omega_cl[l]);
		SUM=-para+F*F*exp(para)+SUM;
	}	
	F=para-D_para->tau_r_p[c]-exp(D_para->beta_c[c]+D_para->r_o_l[l])*D_para->B;
	density=F*F*exp(D_para->sigma_tau_r[c])+SUM; 
	return(-0.5*density); 
}



double MCMC_nu_l_JHM(struct_data_JHM *D,struct_para_JHM *D_para,struct_priors_JHM *D_priors,double para,int c,int l, int m){
	double density,F,SUM=0;
	int n,mm,nn,ll;
	for (c=0;c<2;c++){
		ll=c*D->L+l;
		for (m=0;m<D->NoORF[ll];m++){
			mm=D->NoSUM[ll]+m;
			for (n=0;n<D->NoTIME[mm];n++){
				nn=c*D->SHIFTlmn+l*D->M*D->N + m*D->N+n;
				F=D->y[nn]-logistic_function_E_JHM(D->x[nn],D_para->K_clm[mm], D_para->r_clm[mm],D_para->P);
				SUM=-para+(1+D_priors->df2)*log(1+pow(F,2)*exp(para+D_para->upsilon_c[c])/D_priors->df2)+SUM;
/*-para+F*F*exp(para+D_para->upsilon_c[c])+SUM;*/
			}
		}
	}
	F=para-D_para->nu_p;
	density=F*F*exp(D_para->sigma_nu)+SUM; 
	return(-0.5*density); 
}

double MCMC_sigma_nu_JHM(struct_data_JHM *D,struct_para_JHM *D_para,struct_priors_JHM *D_priors,double para,int c,int l, int m){
	double density,F,SUM=0;
	for (l=0;l<D->L;l++){
		F=D_para->nu_l[l]-D_para->nu_p;
		SUM=-para+F*F*exp(para)+SUM;
	}	
	F=para-D_priors->eta_nu;
	density=F*F*D_priors->psi_nu+SUM; 
	return(-0.5*density); 
}
/*
double MCMC_gamma_cl(struct_data_JHM *D,struct_para_JHM *D_para,struct_priors_JHM *D_priors,double para,int c,int l, int m){
	double density,F,SUM=0;
	int ll=D->L+l,mm;
	for (m=0;m<D->NoORF[ll];m++){
		mm=D->NoSUM[ll]+m;
		F=D_para->K_clm[mm]-(D_para->alpha_c[1]+D_para->K_o_l[l]+D_para->delta_l[l]*para);
		SUM=F*F*exp(D_para->tau_K_cl[ll])+SUM;
	}	
	F=para;
	density=F*F*exp(D_para->sigma_gamma)+SUM; 
	return(-0.5*density); 
}

double MCMC_omega_cl(struct_data_JHM *D,struct_para_JHM *D_para,struct_priors_JHM *D_priors,double para,int c,int l, int m){
	double density,F,SUM=0;
	int ll=D->L+l,mm;
	for (m=0;m<D->NoORF[ll];m++){
		mm=D->NoSUM[ll]+m;
		F=D_para->r_clm[mm]-(D_para->beta_c[1]+D_para->r_o_l[l]+D_para->delta_l[l]*para);
		SUM=F*F*exp(D_para->tau_r_cl[ll])+SUM;
	}	
	F=para;
	density=F*F*exp(D_para->sigma_omega)+SUM; 
	return(-0.5*density); 
}

double MCMC_sigma_gamma(struct_data_JHM *D,struct_para_JHM *D_para,struct_priors_JHM *D_priors,double para,int c,int l, int m){
	double density,F,SUM=0;
	for (l=0;l<D->L;l++){
		F=D_para->gamma_cl[l];
		SUM=-para+F*F*exp(para)+SUM;
	}	
	F=para-D_priors->eta_gamma;
	density=F*F*D_priors->psi_gamma+SUM; 
	return(-0.5*density); 
}

double MCMC_sigma_omega(struct_data_JHM *D,struct_para_JHM *D_para,struct_priors_JHM *D_priors,double para,int c,int l, int m){
	double density,F,SUM=0;
	for (l=0;l<D->L;l++){
		F=D_para->omega_cl[l];
		SUM=-para+F*F*exp(para)+SUM;
	}	
	F=para-D_priors->eta_omega;
	density=F*F*D_priors->psi_omega+SUM; 
	return(-0.5*density); 
	}*/


double MCMC_alpha_c_JHM(struct_data_JHM *D,struct_para_JHM *D_para,struct_priors_JHM *D_priors,double para,int c,int l, int m){
	double density,F,SUM=0;
	int ll,mm;
	for (l=0;l<D->L;l++){
		ll=D->L+l;
 		for (m=0;m<D->NoORF[ll];m++){
			mm=D->NoSUM[ll]+m;
			F=D_para->K_clm[mm]-(para+D_para->K_o_l[l]+D_para->delta_l[l]*D_para->gamma_cl[l]);
			SUM=F*F*exp(D_para->tau_K_cl[ll])+SUM;
		}
	}	
	F=para-D_priors->alpha_mu;
	density=F*F*D_priors->eta_alpha+SUM; 

	return(-0.5*density); 
}

double MCMC_beta_c_JHM(struct_data_JHM *D,struct_para_JHM *D_para,struct_priors_JHM *D_priors,double para,int c,int l, int m){
	double density,F,SUM=0;
	int ll,mm;
	for (l=0;l<D->L;l++){
		ll=D->L+l;
 		for (m=0;m<D->NoORF[ll];m++){
			mm=D->NoSUM[ll]+m;
			F=D_para->r_clm[mm]-(para+D_para->r_o_l[l]+D_para->delta_l[l]*D_para->omega_cl[l]);
			SUM=F*F*exp(D_para->tau_r_cl[ll])+SUM;
		}
	}	
	F=para-D_priors->beta_mu;	
	density=F*F*D_priors->eta_beta+SUM; 
	return(-0.5*density); 
}

double MCMC_upsilon_c_JHM(struct_data_JHM *D,struct_para_JHM *D_para,struct_priors_JHM *D_priors,double para,int c,int l, int m){
	double density,F,SUM=0;
	int n,mm,nn,ll;
	for (l=0;l<D->L;l++){
		ll=D->L+l;
		for (m=0;m<D->NoORF[ll];m++){
			mm=D->NoSUM[ll]+m;
			for (n=0;n<D->NoTIME[mm];n++){
				nn=D->SHIFTlmn+l*D->M*D->N + m*D->N + n;
				F=D->y[nn]-logistic_function_E_JHM(D->x[nn],D_para->K_clm[mm], D_para->r_clm[mm],D_para->P);
				SUM=-para+(1+D_priors->df2)*log(1+pow(F,2)*exp(D_para->nu_l[l]+para)/D_priors->df2)+SUM;				
/*SUM=-para+F*F*exp(D_para->nu_l[l]+para)+SUM;*/
			}
		}
	}
	F=para-D_priors->upsilon_mu;
	density=F*F*exp(D_para->sigma_upsilon)+SUM; 
	return(-0.5*density); 
}

double MCMC_sigma_upsilon_JHM(struct_data_JHM *D,struct_para_JHM *D_para,struct_priors_JHM *D_priors,double para,int c,int l, int m){
	double density,F,SUM=0;
	F=D_para->upsilon_c[1]-D_priors->upsilon_mu;
	SUM=-para+F*F*exp(para)+SUM;
	F=para-D_priors->eta_upsilon;
	density=F*F*D_priors->psi_upsilon+SUM; 
	return(-0.5*density); 
}

double MCMC_P_JHM(struct_data_JHM *D,struct_para_JHM *D_para,struct_priors_JHM *D_priors,double para,int c,int l, int m){
	double density,F,SUM=0;
	int n,mm,nn,ll;
	for (c=0;c<2;c++){
		for (l=0;l<D->L;l++){
		ll=c*D->L+l;
        		for (m=0;m<D->NoORF[ll];m++){
				mm=D->NoSUM[ll]+m;
				for (n=0;n<D->NoTIME[mm];n++){
					nn=c*D->SHIFTlmn+l*D->M*D->N + m*D->N + n;
					F=D->y[nn]-logistic_function_E_JHM(D->x[nn],D_para->K_clm[mm], D_para->r_clm[mm],para);
					SUM=(1+D_priors->df2)*log(1+pow(F,2)*exp(D_para->nu_l[l]+D_para->upsilon_c[c])/D_priors->df2)+SUM;
					/*SUM=F*F*exp(D_para->nu_l[l]+D_para->upsilon_c[c])+SUM;*/
				}
			}
		}
	}
	F=para-D_priors->P_mu;
	density=F*F*D_priors->eta_P+SUM;
	return(-0.5*density); 
}

/*Data Aug*/


double aug_delta_l_JHM(struct_data_JHM *D,struct_para_JHM *D_para,struct_priors_JHM *D_priors,int l){
	double u,para,PRODa=1,PRODb=1,F,F1;
	int m,ll=D->L+l,mm;
        for (m=0;m<D->NoORF[ll];m++){
		mm=D->NoSUM[ll]+m;
		F=D_para->K_clm[mm]-(D_para->alpha_c[1]+D_para->K_o_l[l]);
		F1=D_para->r_clm[mm]-(D_para->beta_c[1]+D_para->r_o_l[l]);
		PRODa=PRODa*(
			sqrt(exp(D_para->tau_K_cl[ll])) *  exp(-0.5*exp(D_para->tau_K_cl[ll])*F*F)*
			sqrt(exp(D_para->tau_r_cl[ll])) *  exp(-0.5*exp(D_para->tau_r_cl[ll])*F1*F1) );
	}	
	PRODa=(1-D_priors->p)*PRODa;

        for (m=0;m<D->NoORF[ll];m++){
		mm=D->NoSUM[ll]+m;
		F=D_para->K_clm[mm]-(D_para->alpha_c[1]+D_para->K_o_l[l]+D_para->gamma_cl[l]);
		F1=D_para->r_clm[mm]-(D_para->beta_c[1]+D_para->r_o_l[l]+D_para->omega_cl[l]);
		PRODb=PRODb*(
			 sqrt(exp(D_para->tau_K_cl[ll])) *  exp(-0.5*exp(D_para->tau_K_cl[ll])*F*F)*
			 sqrt(exp(D_para->tau_r_cl[ll])) *  exp(-0.5*exp(D_para->tau_r_cl[ll])*F1*F1) );
		}
	PRODb=D_priors->p*PRODb;

	u=runif(0,1);
	if (u<(PRODb/(PRODa+PRODb))){para=1;}
	else{para=0;}
	return(para); 
}



double MCMC_K_p_JHM(struct_data_JHM *D,struct_para_JHM *D_para,struct_priors_JHM *D_priors,double para,int c,int l, int m){
  double density,F,SUM=0;
  for (l=0;l<D->L;l++){
    F=exp(D_para->K_o_l[l])-exp(para);
    SUM+=(1+D_priors->df)*log(1+pow(F,2)*exp(D_para->sigma_K_o)/D_priors->df);
  }
  F=para-(D_priors->K_mu);
  density=F*F*D_priors->eta_K_p+SUM; 
  return(-0.5*density); 
}

double MCMC_K_o_l_JHM(struct_data_JHM *D,struct_para_JHM *D_para,struct_priors_JHM *D_priors,double para,int c,int l, int m){
  double density,F,SUM=0;
	int ll,mm;
	for (c=0;c<2;c++){
		ll=c*D->L+l;
		for (m=0;m<D->NoORF[ll];m++){
			mm=D->NoSUM[ll]+m;
			F=D_para->K_clm[mm]-(D_para->alpha_c[c]+log(para)+c*D_para->delta_l[l]*D_para->gamma_cl[l]);
      			SUM+=F*F*exp(D_para->tau_K_cl[ll]);
    		}
	}
	F=(para)-exp(D_para->K_p);
	density=(1+D_priors->df)*log(1+pow(F,2)*exp(D_para->sigma_K_o)/D_priors->df)+SUM;
  return(-0.5*density); 
}

double MCMC_r_p_JHM(struct_data_JHM *D,struct_para_JHM *D_para,struct_priors_JHM *D_priors,double para,int c,int l, int m){
  double density,F,SUM=0;
  for (l=0;l<D->L;l++){
    F=exp(D_para->r_o_l[l])-exp(para);
    SUM+=(1+D_priors->df)*log(1+pow(F,2)*exp(D_para->sigma_r_o)/D_priors->df);
  }
  F=para-(D_priors->r_mu);
  density=F*F*D_priors->eta_r_p+SUM; 
  return(-0.5*density); 
}


double MCMC_r_o_l_JHM(struct_data_JHM *D,struct_para_JHM *D_para,struct_priors_JHM *D_priors,double para,int c,int l, int m){
  double density,F,SUM=0;
	int ll,mm;
	for (c=0;c<2;c++){
		ll=c*D->L+l;
		for (m=0;m<D->NoORF[ll];m++){
			mm=D->NoSUM[ll]+m;
			F=D_para->r_clm[mm]-(D_para->beta_c[c]+log(para)+c*D_para->delta_l[l]*D_para->omega_cl[l]);
      			SUM+=F*F*exp(D_para->tau_r_cl[ll]);
    		}
	}
	F=(para)-exp(D_para->r_p);
  density=(1+D_priors->df)*log(1+pow(F,2)*exp(D_para->sigma_r_o)/D_priors->df)+SUM;
  return(-0.5*density); 
}


double MCMC_sigma_K_o_JHM(struct_data_JHM *D,struct_para_JHM *D_para,struct_priors_JHM *D_priors,double para,int c,int l, int m){
  double density,SUM=0,F;
  for (l=0;l<D->L;l++){
    F=exp(D_para->K_o_l[l])-exp(D_para->K_p);
    SUM+=(1+D_priors->df)*log(1+pow(F,2)*exp(para)/D_priors->df)-(para);
  }
  F=(para)-(D_priors->eta_K_o);
  density=F*F*D_priors->psi_K_o+SUM; 
  return(-0.5*density); 
}

double MCMC_sigma_r_o_JHM(struct_data_JHM *D,struct_para_JHM *D_para,struct_priors_JHM *D_priors,double para,int c,int l, int m){
  double density,SUM=0,F;
  for (l=0;l<D->L;l++){
    F=exp(D_para->r_o_l[l])-exp(D_para->r_p);
    SUM+=(1+D_priors->df)*log(1+pow(F,2)*exp(para)/D_priors->df)-(para);
  }
  F=(para)-(D_priors->eta_r_o);
  density=F*F*D_priors->psi_r_o+SUM; 
  return(-0.5*density); 
}

double MCMC_gamma_cl_JHM(struct_data_JHM *D,struct_para_JHM *D_para,struct_priors_JHM *D_priors,double para,int c,int l, int m){
  double density,F,SUM=0;
  int ll=D->L+l,mm;
  for (m=0;m<D->NoORF[ll];m++){
    mm=D->NoSUM[ll]+m;
    F=D_para->K_clm[mm]-(D_para->alpha_c[1]+D_para->K_o_l[l]+D_para->delta_l[l]*log(para));
    SUM=F*F*exp(D_para->tau_K_cl[ll])+SUM;
  }
  F=para-1;
  density=(1+D_priors->df)*log(1+pow(F,2)*exp(D_para->sigma_gamma)/D_priors->df)+SUM;
  return(-0.5*density);
}

double MCMC_omega_cl_JHM(struct_data_JHM *D,struct_para_JHM *D_para,struct_priors_JHM *D_priors,double para,int c,int l, int m){
  double density,F,SUM=0;
  int ll=D->L+l,mm;
  for (m=0;m<D->NoORF[ll];m++){
    mm=D->NoSUM[ll]+m;
    F=D_para->r_clm[mm]-(D_para->beta_c[1]+D_para->r_o_l[l]+D_para->delta_l[l]*log(para));
    SUM=F*F*exp(D_para->tau_r_cl[ll])+SUM;
  }
  F=para-1;
  density=(1+D_priors->df)*log(1+pow(F,2)*exp(D_para->sigma_omega)/D_priors->df)+SUM;
  return(-0.5*density);
}

double MCMC_sigma_gamma_JHM(struct_data_JHM *D,struct_para_JHM *D_para,struct_priors_JHM *D_priors,double para,int c,int l, int m){
  double density,F,SUM=0;
  for (l=0;l<D->L;l++){
    F=exp(D_para->gamma_cl[l])-1;
    SUM+=(1+D_priors->df)*log(1+pow(F,2)*exp(para)/D_priors->df)-(para);
  }
  F=para-D_priors->eta_gamma;
  density=F*F*D_priors->psi_gamma+SUM;
  return(-0.5*density);
}

double MCMC_sigma_omega_JHM(struct_data_JHM *D,struct_para_JHM *D_para,struct_priors_JHM *D_priors,double para,int c,int l, int m){
  double density,F,SUM=0;
  for (l=0;l<D->L;l++){
    F=exp(D_para->omega_cl[l])-1;
    SUM+=(1+D_priors->df)*log(1+pow(F,2)*exp(para)/D_priors->df)-(para);
  }
  F=para-D_priors->eta_omega;
  density=F*F*D_priors->psi_omega+SUM;
  return(-0.5*density);
}

double MCMC_tau_K_p_JHM(struct_data_JHM *D,struct_para_JHM *D_para,struct_priors_JHM *D_priors,double para,int c,int l, int m){
  double density,F,SUM=0;
  int ll;
  for (l=0;l<D->L;l++){
    ll=l+c*D->L;
    F=D_para->tau_K_cl[ll]-(para)-exp(D_para->alpha_c[c]+D_para->K_o_l[l])*D_para->A;
    SUM+=F*F*exp(D_para->sigma_tau_K[c]);
  }
  F=D_priors->tau_K_mu-para;
  density=F*F*D_priors->eta_tau_K_p+SUM;
  return(-0.5*density);
}

double MCMC_A_JHM(struct_data_JHM *D,struct_para_JHM *D_para,struct_priors_JHM *D_priors,double para,int c,int l, int m){
  double density,F,SUM=0;
  int ll;
  for (c=0;c<2;c++){
  for (l=0;l<D->L;l++){
    ll=l+c*D->L;
    F=D_para->tau_K_cl[ll]-(D_para->tau_K_p[c])-exp(D_para->alpha_c[c]+D_para->K_o_l[l])*para;
    SUM+=F*F*exp(D_para->sigma_tau_K[c]);
  }
  }
  F=para;
  density=F*F/(6*6)+SUM;
  return(-0.5*density);
}


double MCMC_B_JHM(struct_data_JHM *D,struct_para_JHM *D_para,struct_priors_JHM *D_priors,double para,int c,int l, int m){
  double density,F,SUM=0;
  int ll;
  for (c=0;c<2;c++){
    for (l=0;l<D->L;l++){
      ll=l+c*D->L;
      F=D_para->tau_r_cl[ll]-(D_para->tau_r_p[c])-exp(D_para->beta_c[c]+D_para->r_o_l[l])*para;
      SUM+=F*F*exp(D_para->sigma_tau_r[c]);
    }
  }
  F=para;
  density=F*F/(4*4)+SUM;
  return(-0.5*density);
}


double MCMC_sigma_tau_K_JHM(struct_data_JHM *D,struct_para_JHM *D_para,struct_priors_JHM *D_priors,double para,int c,int l, int m){
  double density,SUM=0,F;
  int ll;
  for (l=0;l<D->L;l++){
    ll=l+c*D->L;
    F=D_para->tau_K_cl[ll]-D_para->tau_K_p[c]-D_para->A*exp(D_para->alpha_c[c]+D_para->K_o_l[l]);
    SUM=-para+F*F*exp(para)+SUM;
  }
  F=para-D_priors->eta_tau_K;
  density=F*F*D_priors->psi_tau_K+SUM;
  return(-0.5*density);
}
double MCMC_tau_r_p_JHM(struct_data_JHM *D,struct_para_JHM *D_para,struct_priors_JHM *D_priors,double para,int c,int l, int m){
  double density,F,SUM=0;
  int ll;
 for (l=0;l<D->L;l++){
   ll=l+c*D->L;
   F=D_para->tau_r_cl[ll]-para-exp(D_para->beta_c[c]+D_para->r_o_l[l])*D_para->B;
    SUM+=F*F*exp(D_para->sigma_tau_r[c]);
  }
  F=D_priors->tau_r_mu-para;
  density=F*F*D_priors->eta_tau_r_p+SUM;
  return(-0.5*density);
}

double MCMC_sigma_tau_r_JHM(struct_data_JHM *D,struct_para_JHM *D_para,struct_priors_JHM *D_priors,double para,int c,int l,int m){
  double density,SUM=0,F;
  int ll;
  for (l=0;l<D->L;l++){
    ll=l+c*D->L;
    F=D_para->tau_r_cl[ll]-(D_para->tau_r_p[c])-exp(D_para->beta_c[c]+D_para->r_o_l[l])*D_para->B;
    SUM=-para+F*F*exp(para)+SUM;
  }
  F=para-D_priors->eta_tau_r;
  density=F*F*D_priors->psi_tau_r+SUM;
  return(-0.5*density);
}



/*Gibbs and MH steps*/
int gibbsandMHloop_JHM(int iter,int thin,struct_data_JHM *D,struct_para_JHM *D_para,struct_priors_JHM *D_priors ,struct_MH_JHM *D_MH,int print,double *OUT,char **HEADER){
int i,j,l,m,mm,c,ll,*T,t;
T=&t;
*T=0;
/* printf("%g %g",D_priors->eta_nu_p,D_priors->eta_nu) ;print=3;*/

	if (print==0){
		printheader_JHM(D,HEADER);
	}
	for (i=0;i<iter;i++){
  for (j=0;j<thin;j++){
			D_para->P=MCMC_base_JHM(D,D_para,D_priors,&D_MH->accept_P,&D_MH->hP,D_para->P,MCMC_P_JHM,-999,-999,-999);

			D_para->alpha_c[1]=MCMC_base_JHM(D,D_para,D_priors,&D_MH->accept_K,&D_MH->halpha,D_para->alpha_c[1],MCMC_alpha_c_JHM,-999,-999,-999);

			D_para->beta_c[1]=MCMC_base_JHM(D,D_para,D_priors,&D_MH->accept_K,&D_MH->halpha,D_para->beta_c[1],MCMC_beta_c_JHM,-999,-999,-999);

			/*			printf("%g %g %g %g\n",(D_para->nu_l[0]),exp( D_para->tau_r_cl[0]),exp(D_para->upsilon_c[1]),(D_para->nu_p));*/
			D_para->sigma_gamma=MCMC_base_JHM(D,D_para,D_priors,&D_MH->accept_K,&D_MH->hK,D_para->sigma_gamma,MCMC_sigma_gamma_JHM,-999,-999,-999);

			D_para->sigma_omega=MCMC_base_JHM(D,D_para,D_priors,&D_MH->accept_K,&D_MH->hK,D_para->sigma_omega,MCMC_sigma_omega_JHM,-999,-999,-999);

			D_para->sigma_upsilon=MCMC_base_JHM(
				D,D_para,D_priors,&D_MH->accept_K,&D_MH->hP,D_para->sigma_upsilon,MCMC_sigma_upsilon_JHM,-999,-999,-999);

			D_para->sigma_nu=MCMC_base_JHM(D,D_para,D_priors,&D_MH->accept_K,&D_MH->hP,D_para->sigma_nu,MCMC_sigma_nu_JHM,-999,-999,-999);

			D_para->sigma_K_o=MCMC_base_JHM(D,D_para,D_priors,&D_MH->accept_K,&D_MH->hK,D_para->sigma_K_o,MCMC_sigma_K_o_JHM,-999,-999,-999);

			D_para->sigma_r_o=MCMC_base_JHM(D,D_para,D_priors,&D_MH->accept_K,&D_MH->hK,D_para->sigma_r_o,MCMC_sigma_r_o_JHM,-999,-999,-999);

			D_para->K_p=MCMC_base_JHM(D,D_para,D_priors,&D_MH->accept_K,&D_MH->hK,D_para->K_p,MCMC_K_p_JHM,-999,-999,-999);

			D_para->r_p=MCMC_base_JHM(D,D_para,D_priors,&D_MH->accept_K,&D_MH->hr,D_para->r_p,MCMC_r_p_JHM,-999,-999,-999);

			D_para->nu_p=gauss_sample_JHM(D,0,D->L,D_para->nu_l,exp(D_para->sigma_nu),D_priors->nu_mu,D_priors->eta_nu_p);

			/*				 D_para->A=MCMC_base(D,D_para,D_priors,&D_MH->accept_K,&D_MH->hnu,D_para->A,MCMC_A,-999,-999,-999);						
						      			      																     																      D_MH->hK=0.05;
   D_para->B=MCMC_base(D,D_para,D_priors,&D_MH->accept_K,&D_MH->hK,D_para->B,MCMC_B,-999,-999,-999);
   D_MH->hK=0.1;*/

			D_para->upsilon_c[1]=MCMC_base_JHM(D,D_para,D_priors,&D_MH->accept_K,&D_MH->hK,D_para->upsilon_c[1],MCMC_upsilon_c_JHM,-999,-999,-999);

			for (c=0;c<2;c++){
			  D_para->tau_K_p[c]=MCMC_base_JHM(D,D_para,D_priors,&D_MH->accept_K,&D_MH->hK,D_para->tau_K_p[c],MCMC_tau_K_p_JHM,c,-999,-999);

			  D_para->tau_r_p[c]=MCMC_base_JHM(D,D_para,D_priors,&D_MH->accept_K,&D_MH->hK,D_para->tau_r_p[c],MCMC_tau_r_p_JHM,c,-999,-999);

			  D_para->sigma_tau_K[c]=MCMC_base_JHM(D,D_para,D_priors,&D_MH->accept_K,&D_MH->hK,D_para->sigma_tau_K[c],MCMC_sigma_tau_K_JHM,c,-999,-999);

			  D_para->sigma_tau_r[c]=MCMC_base_JHM(D,D_para,D_priors,&D_MH->accept_K,&D_MH->hr,D_para->sigma_tau_r[c],MCMC_sigma_tau_r_JHM,c,-999,-999);

			}

			for (l=0;l<D->L;l++){
			  D_MH->hK=0.01;			       
			  D_para->K_o_l[l]=MCMC_base_truncate_low_JHM(0,D,D_para,D_priors,&D_MH->accept_nu,&D_MH->hK,exp(D_para->K_o_l[l]),MCMC_K_o_l_JHM,-999,l,-999);
			  D_MH->hK=0.1;
			  D_para->K_o_l[l]=log(D_para->K_o_l[l]);

				D_para->r_o_l[l]=MCMC_base_truncate_low_JHM(0,D,D_para,D_priors,&D_MH->accept_nu,&D_MH->hK,exp(D_para->r_o_l[l]),MCMC_r_o_l_JHM,-999,l,-999);
				D_para->r_o_l[l]=log(D_para->r_o_l[l]);

				D_para->nu_l[l]=MCMC_base_JHM(D,D_para,D_priors,&D_MH->accept_nu,&D_MH->hnu,D_para->nu_l[l],MCMC_nu_l_JHM,-999,l,-999);

				D_para->gamma_cl[l]=MCMC_base_truncate_low_JHM(0,D,D_para,D_priors,&D_MH->accept_nu,&D_MH->hnu,exp(D_para->gamma_cl[l]),MCMC_gamma_cl_JHM,-999,l,-999);
	  D_para->gamma_cl[l]=log(D_para->gamma_cl[l]);

				D_para->omega_cl[l]=MCMC_base_truncate_low_JHM(0,D,D_para,D_priors,&D_MH->accept_nu,&D_MH->hnu,exp(D_para->omega_cl[l]),MCMC_omega_cl_JHM,-999,l,-999);
	  D_para->omega_cl[l]=log(D_para->omega_cl[l]);

	  D_para->delta_l[l]=aug_delta_l_JHM(D,D_para,D_priors,l);

				for (c=0;c<2;c++){

					ll=c*D->L+l;

					D_para->tau_K_cl[ll]=MCMC_base_truncate_low_JHM(0,D,D_para,D_priors,&D_MH->accept_nu,&D_MH->hnu,D_para->tau_K_cl[ll],MCMC_tau_K_cl_JHM,c,l,-999);

					D_para->tau_r_cl[ll]=MCMC_base_JHM(D,D_para,D_priors,&D_MH->accept_nu,&D_MH->hnu,D_para->tau_r_cl[ll],MCMC_tau_r_cl_JHM,c,l,-999);


					for (m=0;m<D->NoORF[l];m++){ 
						mm=D->NoSUM[ll]+m;
						D_MH->hK=0.01;				
	D_para->K_clm[mm]=MCMC_base_truncate_high_JHM(0,
							D,D_para,D_priors,&D_MH->accept_K,&D_MH->hK,D_para->K_clm[mm],MCMC_K_clm_JHM,c,l,m);
	D_MH->hK=0.1;

	D_para->r_clm[mm]=MCMC_base_truncate_high_JHM(3.5,
							D,D_para,D_priors,&D_MH->accept_r,&D_MH->hr,D_para->r_clm[mm],MCMC_r_clm_JHM,c,l,m);

					}
				}
			}
		}	
		if (print==1){
			printdata_JHM(D,D_para,D_MH,OUT,T);
		}
	}
return 0;
}







