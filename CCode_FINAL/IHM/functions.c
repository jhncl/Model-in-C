#include "headers.h"
#include "functions.h"
#include "print.h"


/*trun*/
double trun_const_low(double x,double m, double l)
{
  return(1-0.5-(  (pow(3*l,0.5))*(x-m)/(((l*pow(x-m,2))+3))         -atan((pow(l/3,0.5))*(x-m)))/M_PI);
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

  logu=log(1-gsl_rng_uniform(RNG));
  if (logaprob>logu){para=can;*accept=*accept+1;}
  return(para); 
}

double MCMC_base_truncate_low(double truncate, gsl_rng *RNG, struct_data *D,struct_para *D_para,struct_priors *D_priors,double *accept,double *h,double para,double (*foo)(struct struct_data *D,struct struct_para *D_para,struct struct_priors *D_priors,double,int,int,int),int c,int l, int m){
  double logu,logaprob,can;
    can=para+gsl_ran_gaussian(RNG,*h);
    if(can<=(truncate)){
      can=para;
    }
 
  logaprob=(*foo)(D,D_para,D_priors,can,c,l,m)-(*foo)(D,D_para,D_priors,para,c,l,m);

  logu=log(1-gsl_rng_uniform(RNG));
  if (logaprob>logu){para=can;*accept=*accept+1;}
  return(para);
}


double MCMC_nu_l(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int c,int l, int m){
  double density,F,SUM=0;
  int ll,mm;
    ll=c*D->L+l;
    for (m=0;m<D->NoORF[ll];m++){
      mm=D->NoSUM[ll]+m;
      F=D->y[mm]-exp(D_para->alpha_c[c]+D_para->Z_l[l]+c*D_para->delta_l[l]*D_para->gamma_cl[l]);
      SUM+=-para+F*F*exp(para);
    }
  F=para-D_para->nu_p;
  density=F*F*exp(D_para->sigma_nu)+SUM; 
  return(-0.5*density); 
}


double MCMC_sigma_nu(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int c,int l, int m){
  double density,SUM=0,F;
  for (l=0;l<2*D->L;l++){
    F=D_para->nu_l[l]-D_para->nu_p;
    SUM=-para+F*F*exp(para)+SUM;
  }
  F=para-D_priors->eta_nu;
  density=F*F*D_priors->psi_nu+SUM; 
  return(-0.5*density); 
}

double MCMC_alpha_c(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int c,int l, int m){
  double density,F,SUM=0;
  int ll,mm;
  for (l=0;l<D->L;l++){
    ll=D->L+l;
    for (m=0;m<D->NoORF[ll];m++){
      mm=D->NoSUM[ll]+m;
      F=D->y[mm]-exp( para+D_para->Z_l[l]+D_para->delta_l[l]*D_para->gamma_cl[l]);
      SUM+=F*F*exp(D_para->nu_l[ll]);
    }
  }
  F=para-D_priors->alpha_mu;
  density=F*F*D_priors->eta_alpha+SUM; 

  return(-0.5*density); 
}

/*
double aug_delta_l(gsl_rng *RNG,struct_data *D,struct_para *D_para,struct_priors *D_priors,int l){
  double u,para,PRODa=1,PRODb=1,F,precision;
  int m,ll=D->L+l,mm;
  for (m=0;m<D->NoORF[ll];m++){
    mm=D->NoSUM[ll]+m;
    F=D->y[mm]-exp(D_para->alpha_c[1]+D_para->Z_l[l]);
    precision=exp(D_para->nu_l[ll]);
    PRODa=PRODa*exp(-0.5*F*F*precision);
  }
  PRODa=(1-D_priors->p)*PRODa;

  for (m=0;m<D->NoORF[ll];m++){
    mm=D->NoSUM[ll]+m;
    F=D->y[mm]-exp(D_para->alpha_c[1]+D_para->Z_l[l]+D_para->gamma_cl[l]);
    precision=exp(D_para->nu_l[ll]);
    PRODb=PRODb*exp(-0.5*F*F*precision);
  }
  PRODb=D_priors->p*PRODb;

  u=gsl_rng_uniform(RNG);
  if (u<(PRODb/(PRODa+PRODb))){para=1;}
  else{para=0;}
  return(para); 
}
*/
double aug_delta_l(gsl_rng *RNG,struct_data *D,struct_para *D_para,struct_priors *D_priors,int l){
  double u,para,PRODa=0,PRODb=0,F;
  int m,ll=D->L+l,mm;
  for (m=0;m<D->NoORF[ll];m++){
    mm=D->NoSUM[ll]+m;
    F=D->y[mm]-exp(D_para->alpha_c[1]+D_para->Z_l[l]);
    PRODa+=-0.5*F*F*exp(D_para->nu_l[ll]);
  }
  PRODa+=log(1-D_priors->p);

  for (m=0;m<D->NoORF[ll];m++){
    mm=D->NoSUM[ll]+m;
    F=D->y[mm]-exp(D_para->alpha_c[1]+D_para->Z_l[l]+D_para->gamma_cl[l]);
    PRODb+=-0.5*F*F*exp(D_para->nu_l[ll]);
  }
  PRODb+=log(D_priors->p);

  u=gsl_rng_uniform(RNG);
  if (log(u/(1-u))<(PRODa-PRODb) ){para=0;}
  else{para=1;}
  return(para);
}

double MCMC_Z_p(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int c,int l, int m){
  double density,F,SUM=0;
  for (l=0;l<D->L;l++){
    F=exp(D_para->Z_l[l])-exp(para);
    SUM+=(1+D_priors->df)*log(1+pow(F,2)*exp(D_para->sigma_Z)/D_priors->df)/**/+2*log(trun_const_low(0,exp(para),exp(D_para->sigma_Z)));
  }
  F=para-D_priors->Z_mu;
  density=F*F*D_priors->eta_Z_p+SUM; 
  return(-0.5*density); 
}

double MCMC_Z_l(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int c,int l, int m){
  double density,F,SUM=0;
  int ll,mm;
  for (c=0;c<2;c++){
    ll=c*D->L+l;
    for (m=0;m<D->NoORF[ll];m++){
      mm=D->NoSUM[ll]+m;
      F=D->y[mm]-para*exp(D_para->alpha_c[c]+c*D_para->delta_l[l]*D_para->gamma_cl[l]);
      SUM+=F*F*exp(D_para->nu_l[ll]);
    }
  }
  F=para-exp(D_para->Z_p);
  density=(1+D_priors->df)*log(1+pow(F,2)*exp(D_para->sigma_Z)/D_priors->df)+SUM;
  return(-0.5*density); 
}

double MCMC_sigma_Z(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int c,int l, int m){
  double density,SUM=0,F;
  for (l=0;l<D->L;l++){
    F=exp(D_para->Z_l[l])-exp(D_para->Z_p);
    SUM+=(1+D_priors->df)*log(1+pow(F,2)*exp(para)/D_priors->df)-(para)/**/+2*log(trun_const_low(0,exp(D_para->Z_p),exp(para)));
  }
  F=para-D_priors->eta_Z;
  density=F*F*D_priors->psi_Z+SUM; 
  return(-0.5*density); 
}

double MCMC_gamma_cl(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int c,int l, int m){
  double density,F,SUM=0;
  int ll=D->L+l,mm;
  for (m=0;m<D->NoORF[ll];m++){
    mm=D->NoSUM[ll]+m;
    F=D->y[mm]-exp( D_para->alpha_c[1]+D_para->Z_l[l]+D_para->delta_l[l]*log(para));
    SUM+=F*F*exp(D_para->nu_l[ll]);
  }
  F=para-1;
  density=(1+D_priors->df)*log(1+pow(F,2)*exp(D_para->sigma_gamma)/D_priors->df)+SUM; 
  return(-0.5*density); 
}

double MCMC_sigma_gamma(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int c,int l, int m){
  double density,F,SUM=0;
  for (l=0;l<D->L;l++){
    F=exp(D_para->gamma_cl[l])-1;
    SUM+=(1+D_priors->df)*log(1+pow(F,2)*exp(para)/D_priors->df)-(para)/**/+2*log(trun_const_low(0,1,exp(para)));
  }
  F=para-D_priors->eta_gamma;
  density=F*F*D_priors->psi_gamma+SUM; 
  return(-0.5*density); 
}



/*Gibbs and MH steps*/

int gibbsandMHloop(int iter,int thin,gsl_rng *RNG,struct_data *D,struct_para *D_para,struct_priors *D_priors ,struct_MH *D_MH,int print){
  int i,j,l,c,ll;
  /* print=3;*/
    if (print==0){printheader(D);}
  for (i=0;i<iter;i++){
    for (j=0;j<thin;j++){
      /*       printf("%g\n",D_para->sigma_Z);*/
      /*      printf("%g %g %g\n",D_para->delta_l[2412],D_para->delta_l[3593],D_para->delta_l[4211]);*/
      D_para->alpha_c[1]=MCMC_base(RNG,D,D_para,D_priors,&D_MH->accept_nu,&D_MH->halpha_c,D_para->alpha_c[1],MCMC_alpha_c,-999,-999,-999);
      D_para->sigma_gamma=MCMC_base(RNG,D,D_para,D_priors,&D_MH->accept_Z,&D_MH->hsigma_gamma,D_para->sigma_gamma,MCMC_sigma_gamma,-999,-999,-999);      
      D_para->sigma_nu=MCMC_base(RNG,D,D_para,D_priors,&D_MH->accept_Z,&D_MH->hsigma_nu,D_para->sigma_nu,MCMC_sigma_nu,-999,-999,-999);
      D_para->Z_p=MCMC_base(RNG,D,D_para,D_priors,&D_MH->accept_Z,&D_MH->hZ_l,D_para->Z_p,MCMC_Z_p,-999,-999,-999);
      D_para->sigma_Z=MCMC_base(RNG,D,D_para,D_priors,&D_MH->accept_Z,&D_MH->hsigma_Z,D_para->sigma_Z,MCMC_sigma_Z,-999,-999,-999);   
      D_para->nu_p=gauss_sample(RNG,D,0,2*D->L,D_para->nu_l,exp(D_para->sigma_nu),D_priors->nu_mu,D_priors->eta_nu_p);

      for (l=0;l<D->L;l++){
	D_para->delta_l[l]=aug_delta_l(RNG,D,D_para,D_priors,l);
	D_para->gamma_cl[l]=MCMC_base_truncate_low(0,RNG,D,D_para,D_priors,&D_MH->accept_Z,&D_MH->hgamma_cl,exp(D_para->gamma_cl[l]),MCMC_gamma_cl,-999,l,-999);
	D_para->gamma_cl[l]=log( D_para->gamma_cl[l]);
	D_para->Z_l[l]=MCMC_base_truncate_low(0,RNG,D,D_para,D_priors,&D_MH->accept_Z,&D_MH->hZ_l,exp(D_para->Z_l[l]),MCMC_Z_l,-999,l,-999);
	D_para->Z_l[l]=log( D_para->Z_l[l]);
	for (c=0;c<2;c++){
	  ll=c*D->L+l;
	  D_para->nu_l[ll]=MCMC_base(RNG,D,D_para,D_priors,&D_MH->accept_Z,&D_MH->hnu_l,D_para->nu_l[ll],MCMC_nu_l,c,l,-999);
	}
      }
    }
    if (print==1){printdata(D,D_para,D_MH);}
  }

  return 0;
}
