#include "headers.h"
#include "functions.h"
#include "print.h"
/*Logistic Growth*/

double logistic_function(double t,double K, double r, double P){
	double output;
	output=(K*P*exp(r*t))/(K+P*(exp(r*t)-1));
	return(output);
}

double logistic_function_E(double t,double K, double r, double P){
	double output;
	output=logistic_function(t,exp(K), exp(r), exp(P));
	return(output);
}

/*Gibbs*/

double gauss_sample(gsl_rng *RNG, struct_data *D ,int start, int N,double x[],double tau,double mu_0,double tau_0){
	double vec,Ndouble,SUM=0;
	int i;
	Ndouble=N;
	for (i=start;i<(start+N);i++){SUM=SUM+x[i];}
	vec=(tau_0*mu_0+tau*SUM)/(tau_0+Ndouble*tau)+gsl_ran_gaussian(RNG,1/sqrt(tau_0+Ndouble*tau));
	return(vec);
}

double gamma_sample(gsl_rng *RNG, struct_data *D ,int start, int N, double x[], double mu,double alpha,double beta){
	double vec,Ndouble,SUM=0;
	int i;
	Ndouble=N;
	for (i=start;i<(start+N);i++){SUM=SUM+(x[i]-mu)*(x[i]-mu);}
	vec=gsl_ran_gamma(RNG,alpha+Ndouble/2,1/(1/beta+SUM/2));
	return(vec);
}

/*MCMC MH*/

double MCMC_base(gsl_rng *RNG, struct_data *D,struct_para *D_para,struct_priors *D_priors,double *accept,double *h,double para,double (*foo)(struct struct_data *D,struct struct_para *D_para,struct struct_priors *D_priors,double,int,int),int l, int m){
	double logu,logaprob,can;
	can=para+gsl_ran_gaussian(RNG,*h);
	logaprob=(*foo)(D,D_para,D_priors,can,l,m)-(*foo)(D,D_para,D_priors,para,l,m);

	logu=gsl_sf_log(1-gsl_rng_uniform(RNG));
	if (logaprob>logu){para=can;*accept=*accept+1;}
	return(para); 
	}

double MCMC_P(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int l, int m){
	double density,F,SUM=0;
	int n,mm,nn;
	for (l=0;l<1;l++){
		for (m=0;m<D->NoORF[l];m++){
			mm=D->NoSUM[l]+m;
			for (n=0;n<D->NoTIME[mm];n++){
				nn=l*D->M*D->N + m*D->N + n;
				F=logistic_function_E(D->x[nn],D_para->K_lm[mm],D_para->r_lm[mm],para);
				SUM=SUM+F*(F-2*D->y[nn])*exp(D_para->nu_l[l]);
			}
		}
	}
	density=para*( para -2*D_priors->P_mu )*D_priors->eta_P_mu+SUM;

	return(-0.5*density); 

}

double MCMC_nu_l(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int l, int m){
	double density,F,SUM=0;
	int n,mm,nn;
	for (m=0;m<D->NoORF[l];m++){
		mm=D->NoSUM[l]+m;
		for (n=0;n<D->NoTIME[mm];n++){
			nn=l*D->M*D->N + m*D->N + n;
			F=logistic_function_E(D->x[nn],D_para->K_lm[mm],D_para->r_lm[mm],D_para->P);
			SUM=SUM+(D->y[nn]-F)*(D->y[nn]-F)*exp(para)-para;
		}
	}
	density=para*( para -2*D_para->nu_p )*exp(D_para->sigma_nu)+SUM; 
	return(-0.5*density); 
}

double MCMC_K_lm(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int l, int m){
	double density,F,SUM=0;
	int n,mm,nn;
	mm=D->NoSUM[l]+m;
	for (n=0;n<D->NoTIME[mm];n++){
		nn=l*D->M*D->N + m*D->N + n;
		F=logistic_function_E(D->x[nn],para, D_para->r_lm[mm],D_para->P);
		SUM=SUM+F*(F-2*D->y[nn])*exp(D_para->nu_l[l]);
	}	
	density=para*(para-2*D_para->K_o_l[l])*exp(D_para->tau_K_l[l])+SUM; 
	return(-0.5*density);
}

double MCMC_r_lm(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int l, int m){
	double density,F,SUM=0;
	int n,mm,nn;
	mm=D->NoSUM[l]+m;
	for (n=0;n<D->NoTIME[mm];n++){
		nn=l*D->M*D->N + m*D->N + n;
		F=logistic_function_E(D->x[nn], D_para->K_lm[mm],para,D_para->P);	
		SUM=SUM+F*(F-2*D->y[nn])*exp(D_para->nu_l[l]);
	}	
	density=para*(para-2*D_para->r_o_l[l])*exp(D_para->tau_r_l[l])+SUM; 
	return(-0.5*density); 
}

/*Gibbs and MH steps*/



int gibbsandMHloop(int iter,int thin,gsl_rng *RNG,struct_data *D,struct_para *D_para,struct_priors *D_priors ,struct_MH *D_MH,int CAPL,int print){
int i,j,l,m,mm;
D->L=gsl_min(D->L,CAPL);

	if (print==0){printheader(D);}	
	for (i=0;i<iter;i++){
	for (j=0;j<thin;j++){
		D_para->P=MCMC_base(RNG,D,D_para,D_priors,&D_MH->accept_P,&D_MH->hP,D_para->P,MCMC_P,-999,-999);
		D_para->sigma_K_o=MCMC_base(RNG,D,D_para,D_priors,&D_MH->accept_P,&D_MH->hP,D_para->sigma_K_o,MCMC_sigma_K_o,-999,-999);
		D_para->sigma_r_o=MCMC_base(RNG,D,D_para,D_priors,&D_MH->accept_P,&D_MH->hP,D_para->sigma_r_o,MCMC_sigma_r_o,-999,-999);
		D_para->sigma_nu=MCMC_base(RNG,D,D_para,D_priors,&D_MH->accept_P,&D_MH->hP,D_para->sigma_nu,MCMC_sigma_nu,-999,-999);
		D_para->K_p=gauss_sample(RNG,D,0,D->L,D_para->K_o_l,exp(D_para->sigma_K_o),D_priors->K_mu,D_priors->eta_K_mu);
		D_para->r_p=gauss_sample(RNG,D,0,D->L,D_para->r_o_l,exp(D_para->sigma_r_o),D_priors->r_mu,D_priors->eta_r_mu);
		D_para->nu_p=gauss_sample(RNG,D,0,D->L,D_para->nu_l,exp(D_para->sigma_nu),D_priors->nu_mu,D_priors->eta_nu_mu);

		for (l=0;l<D->L;l++){
			D_para->tau_K_l[l]=MCMC_base(RNG,D,D_para,D_priors,&D_MH->accept_P,&D_MH->hP,D_para->tau_K_l[l],MCMC_tau_K_l,l,-999);
			D_para->tau_r_l[l]=MCMC_base(RNG,D,D_para,D_priors,&D_MH->accept_P,&D_MH->hP,D_para->tau_r_l[l],MCMC_tau_r_l,l,-999);

			D_para->K_o_l[l]=gauss_sample(RNG,D,D->NoSUM[l],D->NoORF[l],D_para->K_lm,exp(D_para->tau_K_l[l]),D_para->K_p,exp(D_para->sigma_K_o));
			D_para->r_o_l[l]=gauss_sample(RNG,D,D->NoSUM[l],D->NoORF[l],D_para->r_lm,exp(D_para->tau_r_l[l]),D_para->r_p,exp(D_para->sigma_r_o));
			D_para->nu_l[l]=MCMC_base(RNG,D,D_para,D_priors,&D_MH->accept_nu,&D_MH->hnu,D_para->nu_l[l],MCMC_nu_l,l,-999);
			for (m=0;m<D->NoORF[l];m++){ 
				mm=D->NoSUM[l]+m;
				D_para->K_lm[mm]=MCMC_base(RNG,D,D_para,D_priors,&D_MH->accept_K,&D_MH->hK,D_para->K_lm[mm],MCMC_K_lm,l,m);
				D_para->r_lm[mm]=gsl_min(3.5,MCMC_base(RNG,D,D_para,D_priors,&D_MH->accept_r,&D_MH->hr,D_para->r_lm[mm],MCMC_r_lm,l,m));
			}
		}
	}
	if (print==1){printdata(D,D_para,D_MH);}
	}
return 0;
}


double MCMC_sigma_K_o(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int l, int m){
	double density,SUM=0;
	for (l=0;l<D->L;l++){
		SUM=SUM-para+(D_para->K_o_l[l]-D_para->K_p)*(D_para->K_o_l[l]-D_para->K_p)*exp(para);
	}	
	density=para*(para-2*D_priors->sigma_K)*D_priors->phi_K+SUM; 
	return(-0.5*density); 
}


double MCMC_sigma_r_o(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int l, int m){
	double density,SUM=0;
	for (l=0;l<D->L;l++){
		SUM=SUM-para+(D_para->r_o_l[l]-D_para->r_p)*(D_para->r_o_l[l]-D_para->r_p)*exp(para);
	}	
	density=para*(para-2*D_priors->sigma_r)*D_priors->phi_r+SUM; 
	return(-0.5*density); 
}

double MCMC_tau_K_l(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int l, int m){
	double density,SUM=0,F;
	int mm;
	for (m=0;m<D->NoORF[l];m++){
		mm=D->NoSUM[l]+m;
		F=D_para->K_lm[mm]-D_para->K_o_l[l];
		SUM=SUM-para+F*F*exp(para);
	}	
	density=para*(para-2*D_priors->eta_K_p)*D_priors->psi_K_o+SUM; 
	return(-0.5*density); 
}

double MCMC_tau_r_l(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int l, int m){
	double density,SUM=0,F;
	int mm;
	for (m=0;m<D->NoORF[l];m++){
		mm=D->NoSUM[l]+m;
		F=D_para->r_lm[mm]-D_para->r_o_l[l];
		SUM=SUM-para+F*F*exp(para);
	}	
	density=para*(para-2*D_priors->eta_r_p)*D_priors->psi_r_o+SUM; 
	return(-0.5*density); 
}

double MCMC_sigma_nu(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int l, int m){
	double density,SUM=0,F;
	for (l=0;l<D->L;l++){
		F=D_para->nu_l[l]-D_para->nu_p;
		SUM=SUM-para+F*F*exp(para);
	}	
	density=para*(para-2*D_priors->eta_nu)*D_priors->psi_nu+SUM; 
	return(-0.5*density); 
}
