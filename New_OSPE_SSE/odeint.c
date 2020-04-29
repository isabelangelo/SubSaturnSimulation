#include <stdio.h>
#include <math.h>
#define NRANSI
#include "nrutil.h"
#define MAXSTP 2000000000//10000000//1000000//2000000000
#define TINY 1.0e-30
#include "HTI.h"


double Roche(double r,double q);

extern int kmax,kount,flgout;
extern double *xp,**yp, dxsav,**vp,*ve;


void odeint(double ystart[], int nvar, double x1, double x2, double eps, double h1,
	double hmin, int *nok, int *nbad,
	void (*derivs)(double, double [], double []),
	void (*rkqs)(double [], double [], int, double *, double, double, double [],
	double *, double *, void (*)(double, double [], double [])))
{
	int nstp,i;
	double xsav=0,x,hnext,hdid,h;
	double *yscal,*y,*dydx;
	double e1,e2,htot,G2,a1,a2,P1,R1,R2,q,Naozval;
	yscal=vector(1,nvar);
	y=vector(1,nvar);
	dydx=vector(1,nvar);
	x=x1;
	h=SIGN(h1,x2-x1);
	*nok = (*nbad) = kount = 0;
       	
	for (i=1;i<=nvar;i++) y[i]=ystart[i];
	if (kmax > 0) xsav=x-dxsav*2.0;

	for (nstp=1;nstp<=MAXSTP;nstp++) { 
	  (*derivs)(x,y,dydx);
	  /*if (flgout==1){///finshed to integarting
	    break;
	    }*/
	  for (i=1;i<=nvar;i++) yscal[i]=fabs(y[i])+fabs(dydx[i]*h)+TINY;
	  if (kmax > 0 && kount < kmax-1 && fabs(x-xsav) > fabs(dxsav)) { 
	    xp[++kount]=x;  
	    for (i=1;i<=nvar;i++) {
	      yp[i][kount]=y[i];
	    }
	    for (i=1;i<=n_orb_param;i++){
	      vp[i][kount]=ve[i];
	    }
	    xsav=x;
	    if (flgout==1){
	      fprintf(stderr, " Oh no, you got nan's, something is wrong, \n");
	      fprintf(stderr, " terminating the intergation now....\n");
              sur=-1;
              return;
	    }

	    //smadar:: more conditions to quit:
	    e1 = sqrt(yp[1][kount]*yp[1][kount] + yp[2][kount]*yp[2][kount]);
	    htot=yp[8][kount];
	    a1= SQR(htot)/ (k2 * (yp[16][kount]+yp[20][kount]) *(1-SQR(e1)));
	    R1=yp[15][kount];
	    R2=yp[19][kount];
	    q=yp[16][kount]/(yp[20][kount]+yp[16][kount]);
	    Roche1=Roche(R1,q);
	    q=yp[20][kount]/(yp[16][kount]+yp[20][kount]);
	    Roche2=Roche(R2,q);
	    if (a1*(1-e1)<Roche1 || a1*(1-e1)<Roche2){

	      //    if ( vp[7][kount]>=a1*(1-e1)*Roche1 || vp[8][kount]>=a1*(1-e1)*Roche2 ){

	      fprintf(stderr, "   ******************* Roche Lobe overflow!!! *******************\n");
	      fprintf(stderr, " star 1: ::::: a1*(1-e1)=%le AU >= Roche1 (=%le AU) ::::: \n",a1*(1-e1),Roche1);
	      fprintf(stderr, " star 2: ::::: a1*(1-e1)=%le AU >= Roche2 (=%le AU) ::::: \n",a1*(1-e1),Roche2);
	      fprintf(stderr, "  the time is: %le years \n",x);
              fprintf(stderr, " Terminating now...\n");
              sur=0;
	      if (e1<1e-4) sur=5;
	      flgout=1;
              return;
            }
	    // fprintf(stderr, "a1*(1-e1)=%le Roche=%le\n",a1*(1-e1),Roche);
	    /* if (a1*(1-e1)<=Roche){
              fprintf(stderr, " ::::: a1*(1-e1)(=%le AU) < Roche (=%le AU) ::::: \n",a1*(1-e1),Roche);
              fprintf(stderr, " The planet just died....poor thing....\n");
              fprintf(stderr, " Terminating now...\n");
	       sur=0;
	       return;
	       }*/
	    /* if (a1<=0.03 && e1<0.005) {
	      fprintf(stderr, " The planet got captured....\n");
	      sur=1;
	      return;
	      }*/
	  }
	  if ((x+h-x2)*(x+h-x1) > 0.0) h=x2-x;
	  if ((h-dxsav)*(x+h-x1) > 0.0) h=dxsav;
	  (*rkqs)(y,dydx,nvar,&x,h,eps,yscal,&hdid,&hnext,derivs);
	  // if (hdid == h) ++(*nok); else ++(*nbad);
	  //smadar: change to have another condition
	  if (flgout==1){
	    fprintf(stderr, " Oh no, you got nan's, something is wrong, \n");
	    fprintf(stderr, " terminating the intergation now....\n");
	        
	    if (kmax) {
              xp[++kount]=x;
	      for (i=1;i<=nvar;i++) {
                yp[i][kount]=y[i];
              }
	      for (i=1;i<=n_orb_param;i++){
                vp[i][kount]=ve[i];
              }
	    }
	    sur=-1;
	        
	    return;
	  }

	  if (hdid == h) ++(*nok); else ++(*nbad);

	  if ((x-x2)*(x2-x1) >= 0.0) {
	    for (i=1;i<=nvar;i++) {
	      ystart[i]=y[i];
	    }
	    if (kmax) {
	      xp[++kount]=x;
	      for (i=1;i<=nvar;i++) {
		yp[i][kount]=y[i];
	      }
	      for (i=1;i<=n_orb_param;i++){
		vp[i][kount]=ve[i];
	      }
	   
	    }
	    free_vector(dydx,1,nvar);
	    free_vector(y,1,nvar);
	    free_vector(yscal,1,nvar);
	    //fprintf(stderr, "Integration Successful \n");
	    sur=1;
	    //checking BH Roche limit
	    e2 = sqrt(yp[3][kount]*yp[3][kount] + yp[4][kount]*yp[4][kount]);
	    G2=yp[18][kount];
	    a2= SQR(G2 / ( (yp[16][kount]+yp[20][kount])*vp[4][kount] ) ) *(yp[16][kount]+yp[20][kount]+vp[4][kount]) / (k2 *(1-SQR(e2)));
	    q=vp[4][kount]/(yp[16][kount]+yp[20][kount]);
	    Naozval=pow(3*q,1./3.)*(1+e1)/(1-e2);
	    if (a2/a1 < Naozval){
	      sur=3;
	      sur2=3;
	      return;
	    } 
	    return;
	  }
	  if (fabs(hnext) <= hmin) {
	    fprintf(stderr,"Warnning:  Step size too small in odeint  - n");
	    break;
	  }//nrerror("Step size too small in odeint");
	  h=hnext;
	  //if (nstp==3) exit(1);
	}
	fprintf(stderr,"Too many steps in routine odeint -  - assuming that  TF had e->0\n");
	fprintf(stderr," !!!! Warnning:::: if TF is not on - something is wrong!!!!\n");
	//	nrerror("Too many steps in routine odeint");
}
#undef MAXSTP
#undef TINY
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software ?421.1-9. */