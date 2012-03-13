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

double MCMC_base(gsl_rng *RNG, struct_data *D,struct_para *D_para,struct_priors *D_priors,double *accept,double *h,double para,double (*foo)(struct struct_data *D,struct struct_para *D_para,struct struct_priors *D_priors,double,int,int,int),int c,int l, int m){
	double logu,logaprob,can;
	can=para+gsl_ran_gaussian(RNG,*h);
	logaprob=(*foo)(D,D_para,D_priors,can,c,l,m)-(*foo)(D,D_para,D_priors,para,c,l,m);

	logu=gsl_sf_log(1-gsl_rng_uniform(RNG));
	if (logaprob>logu){para=can;*accept=*accept+1;}
	return(para); 
	}


double MCMC_K_clm(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int c,int l, int m){
	double density,F,SUM=0;
	int n,mm,nn,ll=c*D->L+l;
	mm=D->NoSUM[ll]+m;
	for (n=0;n<D->NoTIME[mm];n++){
		nn=c*D->SHIFTlmn+l*D->M*D->N + m*D->N+n;
		F=D->y[nn]-logistic_function_E(D->x[nn],para, D_para->r_clm[mm],D_para->P);
		SUM=F*F*exp(D_para->nu_l[l]+D_para->upsilon_c[c])+SUM;
	}	
	F=para-(D_para->alpha_c[c]+D_para->K_o_l[l]+c*D_para->delta_l[l]*D_para->gamma_cl[l]);
	density=F*F*exp(D_para->tau_K_cl[ll])+SUM; 
	return(-0.5*density);
}

double MCMC_r_clm(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int c,int l, int m){
	double density,F,SUM=0;
	int n,mm,nn,ll=c*D->L+l;
	mm=D->NoSUM[ll]+m;
	for (n=0;n<D->NoTIME[mm];n++){
		nn=c*D->SHIFTlmn+l*D->M*D->N + m*D->N+n;
		F=D->y[nn]-logistic_function_E(D->x[nn], D_para->K_clm[mm],para,D_para->P);	
		SUM=F*F*exp(D_para->nu_l[l]+D_para->upsilon_c[c])+SUM;
	}	
	F=para-(D_para->beta_c[c]+D_para->r_o_l[l]+c*D_para->delta_l[l]*D_para->omega_cl[l]);
	density=F*F*exp(D_para->tau_r_cl[ll])+SUM; 
	return(-0.5*density); 
}


double MCMC_tau_K_cl(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int c,int l, int m){
	double density,SUM=0,F;
	int mm,ll=c*D->L+l;
	for (m=0;m<D->NoORF[ll];m++){
		mm=D->NoSUM[ll]+m;
		F=D_para->K_clm[mm]-(D_para->alpha_c[c]+D_para->K_o_l[l]+c*D_para->delta_l[l]*D_para->gamma_cl[l]);
		SUM=-para+F*F*exp(para)+SUM;
	}	
	F=para-D_priors->sigma_K;
	density=F*F*D_priors->phi_K+SUM; 
	return(-0.5*density); 
}

double MCMC_tau_r_cl(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int c,int l, int m){
	double density,SUM=0,F;
	int mm,ll=c*D->L+l;
	for (m=0;m<D->NoORF[ll];m++){
		mm=D->NoSUM[ll]+m;
		F=D_para->r_clm[mm]-(D_para->beta_c[c]+D_para->r_o_l[l]+c*D_para->delta_l[l]*D_para->omega_cl[l]);
		SUM=-para+F*F*exp(para)+SUM;
	}	
	F=para-D_priors->sigma_r;
	density=F*F*D_priors->phi_r+SUM; 
	return(-0.5*density); 
}



double MCMC_K_o_l(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int c,int l, int m){
	double density,F,SUM=0;
	int ll,mm;
	for (c=0;c<2;c++){
		ll=c*D->L+l;
		for (m=0;m<D->NoORF[ll];m++){
			mm=D->NoSUM[ll]+m;
			F=D_para->K_clm[mm]-(D_para->alpha_c[c]+para+c*D_para->delta_l[l]*D_para->gamma_cl[l]);
			SUM=F*F*exp(D_para->tau_K_cl[ll])+SUM;
		}
	}
	F=para-D_para->K_p;
	density=F*F*exp(D_para->sigma_K_o)+SUM; 
	return(-0.5*density);
}

double MCMC_r_o_l(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int c,int l, int m){
	double density,F,SUM=0;
	int ll,mm;
	for (c=0;c<2;c++){
		ll=c*D->L+l;
		for (m=0;m<D->NoORF[ll];m++){
			mm=D->NoSUM[ll]+m;
			F=D_para->r_clm[mm]-(D_para->beta_c[c]+para+c*D_para->delta_l[l]*D_para->omega_cl[l]);
			SUM=F*F*exp(D_para->tau_r_cl[ll])+SUM;
		}
	}	
	F=para-D_para->r_p;
	density=F*F*exp(D_para->sigma_r_o)+SUM; 
	return(-0.5*density);
}

double MCMC_nu_l(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int c,int l, int m){
	double density,F,SUM=0;
	int n,mm,nn,ll;
	for (c=0;c<2;c++){
		ll=c*D->L+l;
		for (m=0;m<D->NoORF[ll];m++){
			mm=D->NoSUM[ll]+m;
			for (n=0;n<D->NoTIME[mm];n++){
				nn=c*D->SHIFTlmn+l*D->M*D->N + m*D->N+n;
				F=D->y[nn]-logistic_function_E(D->x[nn],D_para->K_clm[mm], D_para->r_clm[mm],D_para->P);
				SUM=-para+F*F*exp(para+D_para->upsilon_c[c])+SUM;
			}
		}
	}
	F=para-D_para->nu_p;
	density=F*F*exp(D_para->sigma_nu)+SUM; 
	return(-0.5*density); 
}

double MCMC_sigma_K_o(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int c,int l, int m){
	double density,F,SUM=0;
	for (l=0;l<D->L;l++){
		F=D_para->K_o_l[l]-D_para->K_p;
		SUM=-para+F*F*exp(para)+SUM;
	}	
	F=para-D_priors->eta_K_o;
	density=F*F*D_priors->psi_K_o+SUM; 
	return(-0.5*density); 
}

double MCMC_sigma_r_o(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int c,int l, int m){
	double density,F,SUM=0;
	for (l=0;l<D->L;l++){
		F=D_para->r_o_l[l]-D_para->r_p;
		SUM=-para+F*F*exp(para)+SUM;
	}	
	F=para-D_priors->eta_r_o;
	density=F*F*D_priors->psi_r_o+SUM; 
	return(-0.5*density); 
}



double MCMC_sigma_nu(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int c,int l, int m){
	double density,F,SUM=0;
	for (l=0;l<D->L;l++){
		F=D_para->nu_l[l]-D_para->nu_p;
		SUM=-para+F*F*exp(para)+SUM;
	}	
	F=para-D_priors->eta_nu;
	density=F*F*D_priors->psi_nu+SUM; 
	return(-0.5*density); 
}

double MCMC_gamma_cl(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int c,int l, int m){
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

double MCMC_omega_cl(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int c,int l, int m){
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

double MCMC_sigma_gamma(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int c,int l, int m){
	double density,F,SUM=0;
	for (l=0;l<D->L;l++){
		F=D_para->gamma_cl[l];
		SUM=-para+F*F*exp(para)+SUM;
	}	
	F=para-D_priors->eta_gamma;
	density=F*F*D_priors->psi_gamma+SUM; 
	return(-0.5*density); 
}

double MCMC_sigma_omega(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int c,int l, int m){
	double density,F,SUM=0;
	for (l=0;l<D->L;l++){
		F=D_para->omega_cl[l];
		SUM=-para+F*F*exp(para)+SUM;
	}	
	F=para-D_priors->eta_omega;
	density=F*F*D_priors->psi_omega+SUM; 
	return(-0.5*density); 
}


double MCMC_alpha_c(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int c,int l, int m){
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

double MCMC_beta_c(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int c,int l, int m){
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

double MCMC_upsilon_c(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int c,int l, int m){
	double density,F,SUM=0;
	int n,mm,nn,ll;
	for (l=0;l<D->L;l++){
		ll=D->L+l;
		for (m=0;m<D->NoORF[ll];m++){
			mm=D->NoSUM[ll]+m;
			for (n=0;n<D->NoTIME[mm];n++){
				nn=D->SHIFTlmn+l*D->M*D->N + m*D->N + n;
				F=D->y[nn]-logistic_function_E(D->x[nn],D_para->K_clm[mm], D_para->r_clm[mm],D_para->P);
				SUM=-para+F*F*exp(D_para->nu_l[l]+para)+SUM;
			}
		}
	}
	F=para-D_priors->upsilon_mu;
	density=F*F*exp(D_para->sigma_upsilon)+SUM; 
	return(-0.5*density); 
}

double MCMC_sigma_upsilon(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int c,int l, int m){
	double density,F,SUM=0;
	F=D_para->upsilon_c[1]-D_priors->upsilon_mu;
	SUM=-para+F*F*exp(para)+SUM;
	F=para-D_priors->eta_upsilon;
	density=F*F*D_priors->psi_upsilon+SUM; 
	return(-0.5*density); 
}

double MCMC_P(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int c,int l, int m){
	double density,F,SUM=0;
	int n,mm,nn,ll;
	for (c=0;c<2;c++){
		for (l=0;l<D->L;l++){
		ll=c*D->L+l;
        		for (m=0;m<D->NoORF[ll];m++){
				mm=D->NoSUM[ll]+m;
				for (n=0;n<D->NoTIME[mm];n++){
					nn=c*D->SHIFTlmn+l*D->M*D->N + m*D->N + n;
					F=D->y[nn]-logistic_function_E(D->x[nn],D_para->K_clm[mm], D_para->r_clm[mm],para);
					SUM=F*F*exp(D_para->nu_l[l]+D_para->upsilon_c[c])+SUM;
				}
			}
		}
	}
	F=para-D_priors->P_mu;
	density=F*F*D_priors->eta_P+SUM;
	return(-0.5*density); 
}

/*Data Aug*/



double aug_delta_l(gsl_rng *RNG,struct_data *D,struct_para *D_para,struct_priors *D_priors,int l){
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

	u=gsl_rng_uniform(RNG);
	if (u<(PRODb/(PRODa+PRODb))){para=1;}
	else{para=0;}
	return(para); 
}

/*Gibbs and MH steps*/

void update_lik_a(double para_mean[],double para_prec[],double scale[],double interaction[],struct_data *D,struct_para *D_para,double *mu,double *sig)
{
  int m,l,mm,c,ll;
  double sum=0;
	for (l=0;l<D->L;l++){
		for (c=0;c<2;c++){
			ll=c*D->L+l;
			sum=0;
			for (m=0;m<D->NoORF[ll];m++){
				mm=D->NoSUM[ll]+m;
				sum+=(para_mean[mm]-scale[c]-interaction[l]*D_para->delta_l[l]);
			}
		} 
	mu[l]=sum/(D->NoORF[l]+D->NoORF[l+D->L]);
  	sig[l]=1/(exp(para_prec[l])*D->NoORF[l]+exp(para_prec[l+D->L])*D->NoORF[l+D->L]);
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

int gibbsandMHloop(int iter,int thin,gsl_rng *RNG,struct_data *D,struct_para *D_para,struct_priors *D_priors ,struct_MH *D_MH,int print){
int i,j,l,m,mm,c,ll;
	if (print==0){
		printheader(D);
	}
	for (i=0;i<iter;i++){
		for (j=0;j<thin;j++){
			D_para->P=MCMC_base(RNG,D,D_para,D_priors,&D_MH->accept_P,&D_MH->hP,D_para->P,MCMC_P,-999,-999,-999);

D_MH->hK=0.00003;
			D_para->alpha_c[1]=MCMC_base(RNG,D,D_para,D_priors,&D_MH->accept_K,&D_MH->hK,D_para->alpha_c[1],MCMC_alpha_c,-999,-999,-999);
			D_para->beta_c[1]=MCMC_base(RNG,D,D_para,D_priors,&D_MH->accept_K,&D_MH->hK,D_para->beta_c[1],MCMC_beta_c,-999,-999,-999);
D_MH->hK=0.0001;
			D_para->sigma_gamma=MCMC_base(RNG,D,D_para,D_priors,&D_MH->accept_K,&D_MH->hK,D_para->sigma_gamma,MCMC_sigma_gamma,-999,-999,-999);
			D_para->sigma_omega=MCMC_base(RNG,D,D_para,D_priors,&D_MH->accept_K,&D_MH->hK,D_para->sigma_omega,MCMC_sigma_omega,-999,-999,-999);
D_MH->hK=1;
			D_para->sigma_upsilon=MCMC_base(
				RNG,D,D_para,D_priors,&D_MH->accept_K,&D_MH->hK,D_para->sigma_upsilon,MCMC_sigma_upsilon,-999,-999,-999);
D_MH->hK=0.0001;
			D_para->sigma_nu=MCMC_base(RNG,D,D_para,D_priors,&D_MH->accept_K,&D_MH->hK,D_para->sigma_nu,MCMC_sigma_nu,-999,-999,-999);
	
D_MH->hK=0.1;
			D_para->sigma_K_o=MCMC_base(RNG,D,D_para,D_priors,&D_MH->accept_K,&D_MH->hK,D_para->sigma_K_o,MCMC_sigma_K_o,-999,-999,-999);
			D_para->sigma_r_o=MCMC_base(RNG,D,D_para,D_priors,&D_MH->accept_K,&D_MH->hK,D_para->sigma_r_o,MCMC_sigma_r_o,-999,-999,-999);
D_MH->hK=0.0001;	


			D_para->nu_p=gauss_sample(RNG,D,0,D->L,D_para->nu_l,exp(D_para->sigma_nu),D_priors->nu_mu,D_priors->eta_nu);
	
			D_para->upsilon_c[1]=MCMC_base(RNG,D,D_para,D_priors,&D_MH->accept_K,&D_MH->hK,D_para->upsilon_c[1],MCMC_upsilon_c,-999,-999,-999);
	
			for (l=0;l<D->L;l++){

				D_para->nu_l[l]=MCMC_base(RNG,D,D_para,D_priors,&D_MH->accept_nu,&D_MH->hnu,D_para->nu_l[l],MCMC_nu_l,-999,l,-999);

				D_para->gamma_cl[l]=MCMC_base(
					RNG,D,D_para,D_priors,&D_MH->accept_nu,&D_MH->hnu,D_para->gamma_cl[l],MCMC_gamma_cl,-999,l,-999);
				D_para->omega_cl[l]=MCMC_base(
					RNG,D,D_para,D_priors,&D_MH->accept_nu,&D_MH->hnu,D_para->omega_cl[l],MCMC_omega_cl,-999,l,-999);
				D_para->delta_l[l]=aug_delta_l(RNG,D,D_para,D_priors,l);

				for (c=0;c<2;c++){
					ll=c*D->L+l;

					D_para->tau_K_cl[ll]=gsl_min(7,MCMC_base(
										 RNG,D,D_para,D_priors,&D_MH->accept_nu,&D_MH->hnu,D_para->tau_K_cl[ll],MCMC_tau_K_cl,c,l,-999));
					D_para->tau_r_cl[ll]=gsl_min(11,MCMC_base(
										  RNG,D,D_para,D_priors,&D_MH->accept_nu,&D_MH->hnu,D_para->tau_r_cl[ll],MCMC_tau_r_cl,c,l,-999));

					for (m=0;m<D->NoORF[l];m++){ 
						mm=D->NoSUM[ll]+m;
D_MH->hK=0.1;

						D_para->K_clm[mm]=gsl_min(0,MCMC_base(
							RNG,D,D_para,D_priors,&D_MH->accept_K,&D_MH->hK,D_para->K_clm[mm],MCMC_K_clm,c,l,m));
D_MH->hK=0.0001;
						D_para->r_clm[mm]=gsl_min(3.5,MCMC_base(
							RNG,D,D_para,D_priors,&D_MH->accept_r,&D_MH->hr,D_para->r_clm[mm],MCMC_r_clm,c,l,m));
					}
				}

			}
			update_lik_a(D_para->K_clm,D_para->tau_K_cl,D_para->gamma_cl,D_para->alpha_c,D,D_para,D_para->star_K_lm,D_para->star_tau_K_l);
			update_lik_b(D_para->star_K_lm,1/exp(D_para->sigma_K_o),D_para->star_tau_K_l,D,D_para,&D_para->star_K_l,&D_para->star_sigma_K_o);
			D_para->K_p=sim_para(RNG,D_para->star_K_l,D_para->star_sigma_K_o,D_priors->K_mu,1/D_priors->eta_K_p);
			for (l=0;l<D->L;l++){	
				D_para->K_o_l[l]=sim_para(RNG,D_para->star_K_lm[l],D_para->star_tau_K_l[l],D_para->K_p,1/exp(D_para->sigma_K_o));
			}
			update_lik_a(D_para->r_clm,D_para->tau_r_cl,D_para->omega_cl,D_para->beta_c,D,D_para,D_para->star_r_lm,D_para->star_tau_r_l);
			update_lik_b(D_para->star_r_lm,1/exp(D_para->sigma_r_o),D_para->star_tau_r_l,D,D_para,&D_para->star_r_l,&D_para->star_sigma_r_o);
			D_para->r_p=sim_para(RNG,D_para->star_r_l,D_para->star_sigma_r_o,D_priors->r_mu,1/D_priors->eta_r_p);
			for (l=0;l<D->L;l++){	
				D_para->r_o_l[l]=sim_para(RNG,D_para->star_r_lm[l],D_para->star_tau_r_l[l],D_para->r_p,1/exp(D_para->sigma_r_o));
			}
		}	
		if (print==1){
			printdata(D,D_para,D_MH);
		}
	}
return 0;
}







