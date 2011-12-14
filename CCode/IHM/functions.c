#include "headers.h"
#include "functions.h"
#include "print.h"

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

double MCMC_base(gsl_rng *RNG, struct_data *D,struct_para *D_para,struct_priors *D_priors,double *accept,double *h,double para,double (*foo)(struct struct_data *D,struct struct_para *D_para,struct struct_priors *D_priors,double,int,int,int),int c,int l, int m){
	double logu,logaprob,can;
	can=para+gsl_ran_gaussian(RNG,*h);
	logaprob=(*foo)(D,D_para,D_priors,can,c,l,m)-(*foo)(D,D_para,D_priors,para,c,l,m);

	logu=gsl_sf_log(1-gsl_rng_uniform(RNG));
	if (logaprob>logu){para=can;*accept=*accept+1;}
	return(para); 
	}

double MCMC_Z_l(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int c,int l, int m){
	double density,F,SUM=0;
	int ll,mm;
	for (c=0;c<2;c++){
		ll=D->L+l;
		for (m=0;m<D->NoORF[ll];m++){
			mm=D->NoSUM[ll]+m;
			F=exp(D_para->alpha[c]+para+c*D_para->delta[l]*D_para->gamma[l]);
			SUM=SUM+F*(F-2*D->y[mm])*exp(D_para->nu_l[l]+D_para->upsilon_c[c]);
		}	
	}
	density=para*( para -2*D_para->Z_p )*exp(D_para->sigma_Z)+SUM;
	return(-0.5*density); 
}


double MCMC_nu_l(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int c,int l, int m){
	double density,F,SUM=0;
	int mm;
	for (c=0;c<2;c++){
		for (m=0;m<D->NoORF[c*D->L+l];m++){
			mm=D->NoSUM[c*D->L+l]+m;
			F=exp(D_para->alpha[c]+D_para->Z_l[l]+c*D_para->delta[l]*D_para->gamma[l]);
			SUM=SUM+(D->y[mm]-F)*(D->y[mm]-F)*exp(para+D_para->upsilon_c[c])-para;
		}
	}
	density=para*( para -2*D_para->nu_p )*exp(D_para->sigma_nu)+SUM; 
	return(-0.5*density); 
}


double aug_delta(gsl_rng *RNG,struct_data *D,struct_para *D_para,struct_priors *D_priors,int l){
	double u,para,SUMa=1,SUMb=1,F,precision;
	int m,ll=D->L+l,mm;
        for (m=0;m<D->NoORF[ll];m++){
		mm=D->NoSUM[ll]+m;
		F=D->y[mm]-exp(D_para->alpha[1]+D_para->Z_l[l]);
		precision=exp(D_para->nu_l[l]+D_para->upsilon_c[1]);
		SUMa=SUMa*exp(-0.5*F*F*precision);
	}	

SUMa=(1-D_priors->p)*SUMa;

        for (m=0;m<D->NoORF[ll];m++){
		mm=D->NoSUM[ll]+m;
		F=D->y[mm]-exp(D_para->alpha[1]+D_para->Z_l[l]+D_para->gamma[l]);
		precision=exp(D_para->nu_l[l]+D_para->upsilon_c[1]);
		SUMb=SUMb*exp(-0.5*F*F*precision);
		}

SUMb=D_priors->p*SUMb;
	u=gsl_rng_uniform(RNG);
	if (u<(SUMb/(SUMa+SUMb))){para=1;}
	else{para=0;}
	return(para); 
}

double MCMC_sigma_Z(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int c,int l, int m){
	double density,SUM=0,F;
	for (l=0;l<D->L;l++){
		F=D_para->Z_l[l]-D_para->Z_p;
		SUM=SUM-para+F*F*exp(para);
	}	
	density=para*(para-2*D_priors->eta_Z)*D_priors->psi_Z+SUM; 
	return(-0.5*density); 
}

double MCMC_sigma_nu(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int c,int l, int m){
	double density,SUM=0,F;
	for (l=0;l<D->L;l++){
		F=D_para->nu_l[l]-D_para->nu_p;
		SUM=SUM-para+F*F*exp(para);
	}	
	density=para*(para-2*D_priors->eta_nu)*D_priors->psi_nu+SUM; 
	return(-0.5*density); 
}

double MCMC_gamma(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int c,int l, int m){
	double density,F,SUM=0;
	int ll=D->L+l,mm;
	for (m=0;m<D->NoORF[ll];m++){
		mm=D->NoSUM[ll]+m;
			F=exp( D_para->alpha[1]+D_para->Z_l[l]+D_para->delta[l]*para );
			SUM=SUM+F*(F-2*D->y[mm])*exp(D_para->nu_l[l]+D_para->upsilon_c[c]);
	}	
	density=para*(para-0)*exp(D_para->sigma_gamma)+SUM; 
	return(-0.5*density); 
}

double MCMC_sigma_gamma(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int c,int l, int m){
	double density,SUM=0;
	for (l=0;l<D->L;l++){
		SUM=SUM-para+D_para->gamma[l]*D_para->gamma[l]*exp(para);
	}	
	density=para*(para-2*D_priors->eta_gamma)*D_priors->psi_gamma+SUM; 
	return(-0.5*density); 
}

double MCMC_alpha(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int c,int l, int m){
	double density,F,SUM=0;
	int ll,mm;
	for (l=0;l<D->L;l++){
		ll=D->L+l;
 		for (m=0;m<D->NoORF[ll];m++){
			mm=D->NoSUM[ll]+m;
			F=exp( para+D_para->Z_l[l]+D_para->delta[l]*D_para->gamma[l] );
			SUM=SUM+F*(F-2*D->y[mm])*exp(D_para->nu_l[l]+D_para->upsilon_c[1]);
		}
	}	
	density=para*(para-2*D_priors->alpha_mu)*D_priors->eta_alpha+SUM; 

	return(-0.5*density); 
}


double MCMC_upsilon_c(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int c,int l, int m){
	double density,F,SUM=0;
	int mm,ll;
	for (l=0;l<D->L;l++){
		ll=c*D->L+l;
		for (m=0;m<D->NoORF[ll];m++){
			mm=D->NoSUM[ll]+m;
			F=exp( D_para->alpha[c]+D_para->Z_l[l]+c*D_para->delta[l]*D_para->gamma[l] );
			SUM=SUM-para+(D->y[mm]-F)*(D->y[mm]-F)*exp(para+D_para->nu_l[l]);
		}
	}
	density=para*( para -2*D_para->upsilon_p )*exp(D_para->sigma_upsilon)+SUM; 
	return(-0.5*density); 
}


double MCMC_sigma_upsilon(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int c,int l, int m){
	double density,SUM=0,F;
	for (c=0;c<2;c++){
		F=D_para->upsilon_c[c]-D_para->upsilon_p;
		SUM=SUM-para+F*F*exp(para);
	}	
	density=para*(para-2*D_priors->eta_upsilon)*D_priors->phi_upsilon+SUM; 
	return(-0.5*density); 
}



/*Gibbs and MH steps*/

int gibbsandMHloop(int iter,int thin,gsl_rng *RNG,struct_data *D,struct_para *D_para,struct_priors *D_priors ,struct_MH *D_MH,int print){
int i,j,l;
	if (print==0){printheader(D);}	
	for (i=0;i<iter;i++){
		for (j=0;j<thin;j++){
			D_para->alpha[1]=MCMC_base(RNG,D,D_para,D_priors,&D_MH->accept_Z,&D_MH->hZ,D_para->alpha[1],MCMC_alpha,-999,-999,-999);
			D_para->sigma_gamma=MCMC_base(RNG,D,D_para,D_priors,&D_MH->accept_Z,&D_MH->hZ,D_para->sigma_gamma,MCMC_sigma_gamma,-999,-999,-999);
			D_para->upsilon_c[0]=MCMC_base(RNG,D,D_para,D_priors,&D_MH->accept_up,&D_MH->hup,D_para->upsilon_c[0],MCMC_upsilon_c,0,-999,-999);
			D_para->upsilon_c[1]=MCMC_base(RNG,D,D_para,D_priors,&D_MH->accept_up,&D_MH->hup,D_para->upsilon_c[1],MCMC_upsilon_c,1,-999,-999);
			D_para->sigma_upsilon=MCMC_base(
				RNG,D,D_para,D_priors,&D_MH->accept_up,&D_MH->hup,D_para->sigma_upsilon,MCMC_sigma_upsilon,-999,-999,-999);

			D_para->upsilon_p=gauss_sample(RNG,D,0,2,D_para->upsilon_c,exp(D_para->sigma_upsilon),D_priors->upsilon_mu,D_priors->eta_upsilon_p);
			D_para->sigma_nu=MCMC_base(RNG,D,D_para,D_priors,&D_MH->accept_nu,&D_MH->hnu,D_para->sigma_nu,MCMC_sigma_nu,-999,-999,-999);
			D_para->Z_p=gauss_sample(RNG,D,0,D->L,D_para->Z_l,exp(D_para->sigma_Z),D_priors->Z_mu,D_priors->eta_Z_p);
			D_para->sigma_Z=MCMC_base(
				RNG,D,D_para,D_priors,&D_MH->accept_Z,&D_MH->hZ,D_para->sigma_Z,MCMC_sigma_Z,-999,-999,-999);         

			D_para->nu_p=gauss_sample(RNG,D,0,D->L,D_para->nu_l,exp(D_para->sigma_nu),D_priors->nu_mu,D_priors->eta_nu_p);

			for (l=0;l<D->L;l++){
				D_para->delta[l]=aug_delta(RNG,D,D_para,D_priors,l);
				D_para->gamma[l]=MCMC_base(RNG,D,D_para,D_priors,&D_MH->accept_nu,&D_MH->hnu,D_para->gamma[l],MCMC_gamma,-999,l,-999);
				D_para->Z_l[l]=MCMC_base(RNG,D,D_para,D_priors,&D_MH->accept_Z,&D_MH->hZ,D_para->Z_l[l],MCMC_Z_l,-999,l,-999);
				D_para->nu_l[l]=MCMC_base(RNG,D,D_para,D_priors,&D_MH->accept_nu,&D_MH->hnu,D_para->nu_l[l],MCMC_nu_l,-999,l,-999);
			}
		}
		if (print==1){printdata(D,D_para,D_MH);}
	}
return 0;
}

