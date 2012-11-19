#include "headers_IHM.h"
#include "functions_IHM.h"
#include "print_IHM.h"

/*Gibbs*/

double gauss_sample_IHM(struct_data_IHM *D ,int start, int N,double x[],double tau,double mu_0,double tau_0){
  double vec,Ndouble,SUM=0;
  int i;
  Ndouble=N;
  for (i=start;i<(start+N);i++){SUM=SUM+x[i];}
  vec=rnorm((tau_0*mu_0+tau*SUM)/(tau_0+Ndouble*tau),1/sqrt(tau_0+Ndouble*tau));
  return(vec);
}

/*MCMC MH*/

double MCMC_base_adaptive_IHM(struct_data_IHM *D,struct_para_IHM *D_para,struct_priors_IHM *D_priors,double *accept,double *h,double para,double (*foo)(struct struct_data_IHM *D,struct struct_para_IHM *D_para,struct struct_priors_IHM *D_priors,double,int,int,int),int c,int l, int m){
  double logu,logaprob,can;
  can=rnorm(para,*h);
  logaprob=(*foo)(D,D_para,D_priors,can,c,l,m)-(*foo)(D,D_para,D_priors,para,c,l,m);

  logu=log(1-runif(0,1));
  if (logaprob>logu){para=can;*accept=*accept+1;*h=*h+0.6;}
  else {
    if (*h>=0.3){*h=*h-0.05;}
  }
  return(para); 
}


double MCMC_base_IHM(struct_data_IHM *D,struct_para_IHM *D_para,struct_priors_IHM *D_priors,double *accept,double *h,double para,double (*foo)(struct struct_data_IHM *D,struct struct_para_IHM *D_para,struct struct_priors_IHM *D_priors,double,int,int,int),int c,int l, int m){
  double logu,logaprob,can;
  can=rnorm(para,*h);
  logaprob=(*foo)(D,D_para,D_priors,can,c,l,m)-(*foo)(D,D_para,D_priors,para,c,l,m);

  logu=log(1-runif(0,1));
  if (logaprob>logu){para=can;*accept=*accept+1;}
  return(para); 
}




double MCMC_base_truncate_low_IHM(double truncate, struct_data_IHM *D,struct_para_IHM *D_para,struct_priors_IHM *D_priors,double *accept,double *h,double para,double (*foo)(struct struct_data_IHM *D,struct struct_para_IHM *D_para,struct struct_priors_IHM *D_priors,double,int,int,int),int c,int l, int m){
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


double MCMC_nu_l_IHM(struct_data_IHM *D,struct_para_IHM *D_para,struct_priors_IHM *D_priors,double para,int c,int l, int m){
  double density,F,SUM=0;
  int ll,mm;
  for (c=0;c<2;c++){
    ll=c*D->L+l;
    for (m=0;m<D->NoORF[ll];m++){
      mm=D->NoSUM[ll]+m;
      F=D->y[mm]-exp(D_para->alpha_c[c]+D_para->Z_l[l]+c*D_para->delta_l[l]*D_para->gamma_cl[l]);
      SUM=-para+F*F*exp(para+D_para->upsilon_c[c])+SUM;
    }
  }
  F=para-D_para->nu_p;
  density=F*F*exp(D_para->sigma_nu)+SUM; 
  return(-0.5*density); 
}


double MCMC_sigma_nu_IHM(struct_data_IHM *D,struct_para_IHM *D_para,struct_priors_IHM *D_priors,double para,int c,int l, int m){
  double density,SUM=0,F;
  for (l=0;l<D->L;l++){
    F=D_para->nu_l[l]-D_para->nu_p;
    SUM=-para+F*F*exp(para)+SUM;
  }
  F=para-D_priors->eta_nu;
  density=F*F*D_priors->psi_nu+SUM; 
  return(-0.5*density); 
}

/*double MCMC_gamma_cl(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int c,int l, int m){
  double density,F,SUM=0;
  int ll=D->L+l,mm;
  for (m=0;m<D->NoORF[ll];m++){
    mm=D->NoSUM[ll]+m;
    F=D->y[mm]-1-exp( D_para->alpha_c[1]+gsl_max(0,gsl_max(0,D_para->Z_l[l])+D_para->delta_l[l]*para ));
    SUM=F*F*exp(D_para->nu_l[l]+D_para->upsilon_c[1])+SUM;
  }
  F=para;
  density=F*F*exp(D_para->sigma_gamma)+SUM; 
  return(-0.5*density); 
}

double MCMC_sigma_gamma(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int c,int l, int m){
  double density,F,SUM=0;
  for (l=0;l<D->L;l++){
    F=D_para->gamma_cl[l];
    SUM=-para+F*F*exp(para)+SUM;
  }x
  F=para-D_priors->eta_gamma;
  density=F*F*D_priors->psi_gamma+SUM; 
  return(-0.5*density); 
  }*/

double MCMC_alpha_c_IHM(struct_data_IHM *D,struct_para_IHM *D_para,struct_priors_IHM *D_priors,double para,int c,int l, int m){
  double density,F,SUM=0;
  int ll,mm;
  for (l=0;l<D->L;l++){
    ll=D->L+l;
    for (m=0;m<D->NoORF[ll];m++){
      mm=D->NoSUM[ll]+m;
      F=D->y[mm]-exp( para+D_para->Z_l[l]+D_para->delta_l[l]*D_para->gamma_cl[l]);
      SUM=F*F*exp(D_para->nu_l[l]+D_para->upsilon_c[1])+SUM;
    }
  }
  F=para-D_priors->alpha_mu;
  density=F*F*D_priors->eta_alpha+SUM; 

  return(-0.5*density); 
}


double MCMC_upsilon_c_IHM(struct_data_IHM *D,struct_para_IHM *D_para,struct_priors_IHM *D_priors,double para,int c,int l, int m){
  double density,F,SUM=0;
  int mm,ll;
  for (l=0;l<D->L;l++){
    ll=D->L+l;
    for (m=0;m<D->NoORF[ll];m++){
      mm=D->NoSUM[ll]+m;
      F=D->y[mm]-exp(D_para->alpha_c[1]+D_para->Z_l[l]+D_para->delta_l[l]*D_para->gamma_cl[l] );
      SUM=-para+F*F*exp(D_para->nu_l[l]+para)+SUM;
    }
  }
  F=para-D_priors->upsilon_mu;
  density=F*F*exp(D_para->sigma_upsilon)+SUM; 
  return(-0.5*density); 
}


double MCMC_sigma_upsilon_IHM(struct_data_IHM *D,struct_para_IHM *D_para,struct_priors_IHM *D_priors,double para,int c,int l, int m){
  double density,F,SUM=0;
  F=D_para->upsilon_c[1]-D_priors->upsilon_mu;
  SUM=-para+F*F*exp(para)+SUM;
  F=para-D_priors->eta_upsilon;
  density=F*F*D_priors->psi_upsilon+SUM; 
  return(-0.5*density); 
}

/*Data Aug*/
double aug_delta_l_IHM(struct_data_IHM *D,struct_para_IHM *D_para,struct_priors_IHM *D_priors,int l){
  double u,para,PRODa=1,PRODb=1,F,precision;
  int m,ll=D->L+l,mm;
  for (m=0;m<D->NoORF[ll];m++){
    mm=D->NoSUM[ll]+m;
    F=D->y[mm]-exp(D_para->alpha_c[1]+D_para->Z_l[l]);
    precision=exp(D_para->nu_l[l]+D_para->upsilon_c[1]);
    PRODa=PRODa*exp(-0.5*F*F*precision);
  }
  PRODa=(1-D_priors->p)*PRODa;

  for (m=0;m<D->NoORF[ll];m++){
    mm=D->NoSUM[ll]+m;
    F=D->y[mm]-exp(D_para->alpha_c[1]+D_para->Z_l[l]+D_para->gamma_cl[l]);
    precision=exp(D_para->nu_l[l]+D_para->upsilon_c[1]);
    PRODb=PRODb*exp(-0.5*F*F*precision);
  }
  PRODb=D_priors->p*PRODb;

  u=runif(0,1);
  if (u<(PRODb/(PRODa+PRODb))){para=1;}
  else{para=0;}
  return(para); 
}

double MCMC_Z_p_IHM(struct_data_IHM *D,struct_para_IHM *D_para,struct_priors_IHM *D_priors,double para,int c,int l, int m){
  double density,F,SUM=0;
  for (l=0;l<D->L;l++){
    F=exp(D_para->Z_l[l])-para;
    SUM+=(1+D_priors->df)*log(1+pow(F,2)*exp(D_para->sigma_Z)/D_priors->df);
  }
  F=D_para->Z_p-para;
  density=F*F*D_priors->eta_Z_p+SUM; 
  return(-0.5*density); 
}

double MCMC_Z_l_IHM(struct_data_IHM *D,struct_para_IHM *D_para,struct_priors_IHM *D_priors,double para,int c,int l, int m){
  double density,F,SUM=0;
  int ll,mm;
  for (c=0;c<2;c++){
    ll=c*D->L+l;
    for (m=0;m<D->NoORF[ll];m++){
      mm=D->NoSUM[ll]+m;
      F=D->y[mm]-/**/para*exp(D_para->alpha_c[c]/*+log(para)*/+c*D_para->delta_l[l]*D_para->gamma_cl[l]);
      SUM+=F*F*exp(D_para->nu_l[l]+D_para->upsilon_c[c]);
    }
  }
  F=para-D_para->Z_p;
  density=(1+D_priors->df)*log(1+pow(F,2)*exp(D_para->sigma_Z)/D_priors->df)+SUM;
  return(-0.5*density); 
}

double MCMC_sigma_Z_IHM(struct_data_IHM *D,struct_para_IHM *D_para,struct_priors_IHM *D_priors,double para,int c,int l, int m){
  double density,SUM=0,F;
  for (l=0;l<D->L;l++){
    F=exp(D_para->Z_l[l])-D_para->Z_p;
    SUM+=(1+D_priors->df)*log(1+pow(F,2)*exp(para)/D_priors->df)-(para);
  }
  F=para-D_priors->eta_Z;
  density=F*F*D_priors->psi_Z+SUM; 
  return(-0.5*density); 
}

double MCMC_gamma_cl_IHM(struct_data_IHM *D,struct_para_IHM *D_para,struct_priors_IHM *D_priors,double para,int c,int l, int m){
  double density,F,SUM=0;
  int ll=D->L+l,mm;
  for (m=0;m<D->NoORF[ll];m++){
    mm=D->NoSUM[ll]+m;
    F=D->y[mm]-exp( D_para->alpha_c[1]+D_para->Z_l[l]+D_para->delta_l[l]*log(para));
    SUM=F*F*exp(D_para->nu_l[l]+D_para->upsilon_c[1])+SUM;
  }
  F=para-1;
  density=(1+D_priors->df)*log(1+pow(F,2)*exp(D_para->sigma_gamma)/D_priors->df)+SUM; 
  return(-0.5*density); 
}

double MCMC_sigma_gamma_IHM(struct_data_IHM *D,struct_para_IHM *D_para,struct_priors_IHM *D_priors,double para,int c,int l, int m){
  double density,F,SUM=0;
  for (l=0;l<D->L;l++){
    F=exp(D_para->gamma_cl[l])-1;
    SUM+=(1+D_priors->df)*log(1+pow(F,2)*exp(para)/D_priors->df)-(para);
  }
  F=para-D_priors->eta_gamma;
  density=F*F*D_priors->psi_gamma+SUM; 
  return(-0.5*density); 
}



/*Gibbs and MH steps*/

int gibbsandMHloop_IHM(int iter,int thin,struct_data_IHM *D,struct_para_IHM *D_para,struct_priors_IHM *D_priors ,struct_MH_IHM *D_MH,int print,double *OUT,char **HEADER){
  int i,j,l,*T,t;
T=&t;
*T=0;
  if (print==0){printheader_IHM(D,HEADER);}
  /*for (i=0;i<iter;i++){
    D_para->alpha_c[1]=MCMC_base_adaptive(
    D,D_para,D_priors,&D_MH->accept_nu,&D_MH->halpha_c,D_para->alpha_c[1],MCMC_alpha_c,-999,-999,-999);

    D_para->sigma_gamma=MCMC_base_adaptive(
    D,D_para,D_priors,&D_MH->accept_Z,&D_MH->hsigma_gamma,D_para->sigma_gamma,MCMC_sigma_gamma,-999,-999,-999);
    D_para->upsilon_c[1]=MCMC_base_adaptive(
    D,D_para,D_priors,&D_MH->accept_up,&D_MH->hupsilon_c,D_para->upsilon_c[1],MCMC_upsilon_c,-999,-999,-999);
    D_para->sigma_upsilon=MCMC_base_adaptive(
    D,D_para,D_priors,&D_MH->accept_Z,&D_MH->hsigma_upsilon,D_para->sigma_upsilon,MCMC_sigma_upsilon,-999,-999,-999);
    D_para->sigma_nu=MCMC_base_adaptive(D,D_para,D_priors,&D_MH->accept_Z,&D_MH->hsigma_nu,D_para->sigma_nu,MCMC_sigma_nu,-999,-999,-999);
    D_para->Z_p=gauss_sample(D,0,D->L,D_para->Z_l,exp(D_para->sigma_Z),D_priors->Z_mu,D_priors->eta_Z_p);
    D_para->sigma_Z=gsl_min(-4,MCMC_base(
    D,D_para,D_priors,&D_MH->accept_Z,&D_MH->hsigma_Z,D_para->sigma_Z,MCMC_sigma_Z,-999,-999,-999));        
    D_para->nu_p=gauss_sample(D,0,D->L,D_para->nu_l,exp(D_para->sigma_nu),D_priors->nu_mu,D_priors->eta_nu_p);

    for (l=0;l<D->L;l++){
    D_para->delta_l[l]=aug_delta_l(D,D_para,D_priors,l);
    D_para->gamma_cl[l]=MCMC_base_adaptive(
    D,D_para,D_priors,&D_MH->accept_Z,&D_MH->hgamma_cl,D_para->gamma_cl[l],MCMC_gamma_cl,-999,l,-999);
    D_para->Z_l[l]=MCMC_base_adaptive(D,D_para,D_priors,&D_MH->accept_Z,&D_MH->hZ_l,D_para->Z_l[l],MCMC_Z_l,-999,l,-999);
    D_para->nu_l[l]=MCMC_base_adaptive(D,D_para,D_priors,&D_MH->accept_Z,&D_MH->hnu_l,D_para->nu_l[l],MCMC_nu_l,-999,l,-999);
    }
    }
}
  */
  /*if (print==1){*/
  for (i=0;i<iter;i++){
    for (j=0;j<thin;j++){
      D_para->alpha_c[1]=MCMC_base_IHM(D,D_para,D_priors,&D_MH->accept_nu,&D_MH->halpha_c,D_para->alpha_c[1],MCMC_alpha_c_IHM,-999,-999,-999);

      D_para->sigma_gamma=MCMC_base_IHM(D,D_para,D_priors,&D_MH->accept_Z,&D_MH->hsigma_gamma,D_para->sigma_gamma,MCMC_sigma_gamma_IHM,-999,-999,-999);

      D_para->upsilon_c[1]=MCMC_base_IHM(D,D_para,D_priors,&D_MH->accept_up,&D_MH->hupsilon_c,D_para->upsilon_c[1],MCMC_upsilon_c_IHM,-999,-999,-999);

      D_para->sigma_upsilon=MCMC_base_IHM(
				      D,D_para,D_priors,&D_MH->accept_Z,&D_MH->hsigma_upsilon,D_para->sigma_upsilon,MCMC_sigma_upsilon_IHM,-999,-999,-999);

      D_para->sigma_nu=MCMC_base_IHM(D,D_para,D_priors,&D_MH->accept_Z,&D_MH->hsigma_nu,D_para->sigma_nu,MCMC_sigma_nu_IHM,-999,-999,-999);

      D_para->Z_p=MCMC_base_IHM(D,D_para,D_priors,&D_MH->accept_Z,&D_MH->hZ_l,D_para->Z_p,MCMC_Z_p_IHM,-999,-999,-999);

      D_para->sigma_Z=MCMC_base_IHM(
				D,D_para,D_priors,&D_MH->accept_Z,&D_MH->hsigma_Z,D_para->sigma_Z,MCMC_sigma_Z_IHM,-999,-999,-999);

      D_para->nu_p=gauss_sample_IHM(D,0,D->L,D_para->nu_l,exp(D_para->sigma_nu),D_priors->nu_mu,D_priors->eta_nu_p);


      for (l=0;l<D->L;l++){
	D_para->delta_l[l]=aug_delta_l_IHM(D,D_para,D_priors,l);

	D_para->gamma_cl[l]=MCMC_base_truncate_low_IHM(0,
						   D,D_para,D_priors,&D_MH->accept_Z,&D_MH->hgamma_cl,exp(D_para->gamma_cl[l]),MCMC_gamma_cl_IHM,-999,l,-999);
D_para->gamma_cl[l]=log( D_para->gamma_cl[l]);

	D_para->Z_l[l]=MCMC_base_truncate_low_IHM(0,D,D_para,D_priors,&D_MH->accept_Z,&D_MH->hZ_l,exp(D_para->Z_l[l]),MCMC_Z_l_IHM,-999,l,-999);
	D_para->Z_l[l]=log( D_para->Z_l[l]);

	D_para->nu_l[l]=MCMC_base_IHM(D,D_para,D_priors,&D_MH->accept_Z,&D_MH->hnu_l,D_para->nu_l[l],MCMC_nu_l_IHM,-999,l,-999);

      }
    }
  if (print==1){printdata_IHM(D,D_para,D_MH,OUT,T);}
  }
  /*}*/


  return 0;
}
