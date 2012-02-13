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

/*MCMC MH*/

double MCMC_base(gsl_rng *RNG, struct_data *D,struct_para *D_para,struct_priors *D_priors,double *accept,double *h,double para,double (*foo)(struct struct_data *D,struct struct_para *D_para,struct struct_priors *D_priors,double,int,int),int l, int m){
	double logu,logaprob,can;
	can=para+gsl_ran_gaussian(RNG,*h);
	logaprob=(*foo)(D,D_para,D_priors,can,l,m)-(*foo)(D,D_para,D_priors,para,l,m);

	logu=gsl_sf_log(1-gsl_rng_uniform(RNG));
	if (logaprob>logu){para=can;*accept=*accept+1;}
	return(para); 
	}


double MCMC_K_lm(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int l, int m){
	double density,F,SUM=0;
	int n,mm,nn;
	mm=D->NoSUM[l]+m;
	for (n=0;n<D->NoTIME[mm];n++){
		nn=l*D->M*D->N + m*D->N + n;
		F=D->y[nn]-logistic_function_E(D->x[nn],para, D_para->r_lm[mm],D_para->P);
		SUM=F*F*exp(D_para->nu_l[l])+SUM;
	}	
	F=para-D_para->K_o_l[l];
	density=F*F*exp(D_para->tau_K_l[l])+SUM; 
	return(-0.5*density);
}

double MCMC_r_lm(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int l, int m){
	double density,F,SUM=0;
	int n,mm,nn;
	mm=D->NoSUM[l]+m;
	for (n=0;n<D->NoTIME[mm];n++){
		nn=l*D->M*D->N + m*D->N + n;
		F=D->y[nn]-logistic_function_E(D->x[nn], D_para->K_lm[mm],para,D_para->P);	
		SUM=F*F*exp(D_para->nu_l[l])+SUM;
	}	
	F=para-D_para->r_o_l[l];
	density=F*F*exp(D_para->tau_r_l[l])+SUM; 
	return(-0.5*density); 
}


double MCMC_tau_K_l(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int l, int m){
	double density,SUM=0,F;
	int mm;
	for (m=0;m<D->NoORF[l];m++){
		mm=D->NoSUM[l]+m;
		F=D_para->K_lm[mm]-D_para->K_o_l[l];
		SUM=-para+F*F*exp(para)+SUM;
	}	
	F=para-D_priors->sigma_K;
	density=F*F*D_priors->phi_K+SUM; 
	return(-0.5*density); 
}

double MCMC_tau_r_l(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int l, int m){
	double density,SUM=0,F;
	int mm;
	for (m=0;m<D->NoORF[l];m++){
		mm=D->NoSUM[l]+m;
		F=D_para->r_lm[mm]-D_para->r_o_l[l];
		SUM=-para+F*F*exp(para)+SUM;
	}	
	F=para-D_priors->sigma_r;
	density=F*F*D_priors->phi_r+SUM; 
	return(-0.5*density); 
}

double MCMC_nu_l(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int l, int m){
	double density,F,SUM=0;
	int n,mm,nn;
	for (m=0;m<D->NoORF[l];m++){
		mm=D->NoSUM[l]+m;
		for (n=0;n<D->NoTIME[mm];n++){
			nn=l*D->M*D->N + m*D->N + n;
			F=D->y[nn]-logistic_function_E(D->x[nn],D_para->K_lm[mm],D_para->r_lm[mm],D_para->P);
			SUM=-para+F*F*exp(para)+SUM;
		}
	}
	F=para-D_para->nu_p;
	density=F*F*exp(D_para->sigma_nu)+SUM; 
	return(-0.5*density); 
}


double MCMC_sigma_K_o(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int l, int m){
	double density,SUM=0,F;
	for (l=0;l<D->L;l++){
		F=D_para->K_o_l[l]-D_para->K_p;
		SUM=-para+F*F*exp(para)+SUM;
	}	
	F=para-D_priors->eta_K_o;
	density=F*F*D_priors->psi_K_o+SUM; 
	return(-0.5*density); 
}


double MCMC_sigma_r_o(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int l, int m){
	double density,SUM=0,F;
	for (l=0;l<D->L;l++){
		F=D_para->r_o_l[l]-D_para->r_p;
		SUM=-para+F*F*exp(para)+SUM;
	}	
	F=para-D_priors->eta_r_o;
	density=F*F*D_priors->psi_r_o+SUM; 
	return(-0.5*density); 
}


double MCMC_sigma_nu(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int l, int m){
	double density,SUM=0,F;
	for (l=0;l<D->L;l++){
		F=D_para->nu_l[l]-D_para->nu_p;
		SUM=-para+F*F*exp(para)+SUM;
	}	
	F=para-D_priors->eta_nu;
	density=F*F*D_priors->psi_nu+SUM; 
	return(-0.5*density); 
}

double MCMC_P(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int l, int m){
	double density,F,SUM=0;
	int n,mm,nn;
	for (l=0;l<D->L;l++){
		for (m=0;m<D->NoORF[l];m++){
			mm=D->NoSUM[l]+m;
			for (n=0;n<D->NoTIME[mm];n++){
				nn=l*D->M*D->N + m*D->N + n;
				F=D->y[nn]-logistic_function_E(D->x[nn],D_para->K_lm[mm],D_para->r_lm[mm],para);
				SUM=F*F*exp(D_para->nu_l[l])+SUM;
			}
		}
	}
	F=para-D_priors->P_mu;
	density=F*F*D_priors->eta_P+SUM;

	return(-0.5*density); 

}



/*Gibbs and MH steps*/

void update_lik_a(double para_mean[],double para_prec[],struct_data *D,struct_para *D_para,double *mu,double *sig)
{
  int m,l,mm;
  double sum=0;

	for (l=0;l<D->L;l++){
		sum=0;
		for (m=0;m<D->NoORF[l];m++){
			mm=D->NoSUM[l]+m;
			sum+=para_mean[mm];
		}  
	mu[l]=sum/D->NoORF[l];
  	sig[l]=1/(exp(para_prec[l])*D->NoORF[l]);
	}
}

void update_lik_b(double pass_mean[],double para_var,double pass_var[],struct_data *D,struct_para *D_para,double *mu,double *sig)
{
  int l;
  double prec=0,sum=0;
  for (l=0;l<D->L;l++) {
    prec+=1/(para_var+pass_var[l]);
    sum+=pass_mean[l]/(para_var+pass_var[l]);
  }
  *mu=sum/prec;
  *sig=1/prec;
}

double sim_para(gsl_rng *RNG,double pass_mean,double pass_var,double para_mean,double para_var)
{
  double mean,var,out_value;
    mean=((pass_var)*para_mean + (para_var)*(pass_mean))
            / (para_var + pass_var);
	var=(pass_var)*(para_var)/(para_var + pass_var);
	out_value=mean+gsl_ran_gaussian(RNG,sqrt(var));
	return(out_value);
}

int gibbsandMHloop(int iter,int thin,gsl_rng *RNG,struct_data *D,struct_para *D_para,struct_priors *D_priors ,struct_MH *D_MH,int CAPL,int print){
int i,j,l,m,mm;
D->L=gsl_min(D->L,CAPL);
	if (print==0){printheader(D);}
	for (i=0;i<iter;i++){
	for (j=0;j<thin;j++){
		D_para->P=MCMC_base(RNG,D,D_para,D_priors,&D_MH->accept_P,&D_MH->hP,D_para->P,MCMC_P,-999,-999);
D_MH->hP=1;		
		D_para->sigma_K_o=MCMC_base(RNG,D,D_para,D_priors,&D_MH->accept_P,&D_MH->hP,D_para->sigma_K_o,MCMC_sigma_K_o,-999,-999);
D_MH->hP=0.2;
		D_para->sigma_r_o=MCMC_base(RNG,D,D_para,D_priors,&D_MH->accept_P,&D_MH->hP,D_para->sigma_r_o,MCMC_sigma_r_o,-999,-999);
		D_para->sigma_nu=MCMC_base(RNG,D,D_para,D_priors,&D_MH->accept_P,&D_MH->hP,D_para->sigma_nu,MCMC_sigma_nu,-999,-999);

		D_para->nu_p=gauss_sample(RNG,D,0,D->L,D_para->nu_l,exp(D_para->sigma_nu),D_priors->nu_mu,D_priors->eta_nu_p);

		for (l=0;l<D->L;l++){
			D_para->tau_K_l[l]=MCMC_base(RNG,D,D_para,D_priors,&D_MH->accept_P,&D_MH->hP,D_para->tau_K_l[l],MCMC_tau_K_l,l,-999);
			D_para->tau_r_l[l]=MCMC_base(RNG,D,D_para,D_priors,&D_MH->accept_P,&D_MH->hP,D_para->tau_r_l[l],MCMC_tau_r_l,l,-999);

			D_para->nu_l[l]=MCMC_base(RNG,D,D_para,D_priors,&D_MH->accept_nu,&D_MH->hnu,D_para->nu_l[l],MCMC_nu_l,l,-999);
			for (m=0;m<D->NoORF[l];m++){ 
				mm=D->NoSUM[l]+m;
				D_para->K_lm[mm]=MCMC_base(RNG,D,D_para,D_priors,&D_MH->accept_K,&D_MH->hK,D_para->K_lm[mm],MCMC_K_lm,l,m);
				D_para->r_lm[mm]=gsl_min(3.5,MCMC_base(RNG,D,D_para,D_priors,&D_MH->accept_r,&D_MH->hr,D_para->r_lm[mm],MCMC_r_lm,l,m));
			}
		}
		
		update_lik_a(D_para->K_lm,D_para->tau_K_l,D,D_para,D_para->star_K_lm,D_para->star_tau_K_l);
		update_lik_b(D_para->star_K_lm,1/exp(D_para->sigma_K_o),D_para->star_tau_K_l,D,D_para,&D_para->star_K_l,&D_para->star_sigma_K_o);
		D_para->K_p=sim_para(RNG,D_para->star_K_l,D_para->star_sigma_K_o,D_priors->K_mu,1/D_priors->eta_K_p);
		for (l=0;l<D->L;l++){	
			D_para->K_o_l[l]=sim_para(RNG,D_para->star_K_lm[l],D_para->star_tau_K_l[l],D_para->K_p,1/exp(D_para->sigma_K_o));
		}
		update_lik_a(D_para->r_lm,D_para->tau_r_l,D,D_para,D_para->star_r_lm,D_para->star_tau_r_l);
		update_lik_b(D_para->star_r_lm,1/exp(D_para->sigma_r_o),D_para->star_tau_r_l,D,D_para,&D_para->star_r_l,&D_para->star_sigma_r_o);
		D_para->r_p=sim_para(RNG,D_para->star_r_l,D_para->star_sigma_r_o,D_priors->r_mu,1/D_priors->eta_r_p);
		for (l=0;l<D->L;l++){	
			D_para->r_o_l[l]=sim_para(RNG,D_para->star_r_lm[l],D_para->star_tau_r_l[l],D_para->r_p,1/exp(D_para->sigma_r_o));
		}

	}
	if (print==1){printdata(D,D_para,D_MH);}
	}
return 0;
}
