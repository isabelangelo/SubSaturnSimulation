#include <math.h>
#include "HTI.h"
#include <stdio.h>
#include <stdlib.h>
#include "nrutil.h"

double *tvec,*mvec,*Rvec,*dmdtvec,*tofdmdtvec,*tvec2,*mvec2,*Rvec2,*dmdtvec2,*tofdmdtvec2;
int length;
extern  int Ndis,Nlin,Ndis2,Nlin2;
extern double *tdis,*tdis2;
extern double *timeSSE,*MtSSE,*logRSSE,*timeSSE2,*MtSSE2,*logRSSE2,*spinSSE,*spinSSE2; //*typeSSE,*MoSSE,*MtSSE,*logLSSE,*logRSSE,*logTSSE,*McSSE,*MenvSSE,*epochSSE,*spinSSE;
extern double *dmdtvec,*tofdmdtvec,*dRdtvec,*mvec,*Rvec,*dmdtvec2,*tofdmdtvec2,*dRdtvec2,*mvec2,*Rvec2;
extern int *indexdis,*indexdis2; //this is the index of discontinutites
extern double R1f,R2f;
extern int Nvar;
extern double t_inter;
extern double *spinvec,*spinvec2,*dspindtvec,*dspindtvec2;

void linint(double *xa,double *ya,int n,double x,double *y);
double dfridr(double (*func)(double), double x, double h, double *err);
void linintMinMax(double *xa,double *ya,int minv,int n,double x,double *y);


/*******************************************/
double dmdt(double t, double m, double m0){
  double mf;
  double dm=0.;
  double err=1.e-8;
  double Deltat;
  int indexval,min,max,length,count,j,iSSE,minIn,maxIn;
  t=t+tMS;

  indexval=indexdis[Ndis];
  minIn=1;
  maxIn=indexdis[Ndis];
  mf=MtSSE[Nlin];


  if (!(SSE1)){
    if (t<=tofdmdtvec[2]){
      dm=dmdtvec[1];
    }
    else if (t>tofdmdtvec[2] && t<=tofdmdtvec[(Ndis+1)*Nvar]){
      linint(tofdmdtvec,dmdtvec,(Ndis+1)*Nvar,t,&dm);
      if (t>tofdmdtvec[(Ndis+1)*Nvar-1]){
	dm=dmdtvec[(Ndis+1)*Nvar-1];
      }
    }
    else if (t>tofdmdtvec[(Ndis+1)*Nvar]){
      if (m>mf){    
	dm=(mf-m)/t_inter;    
      }
      else {
	dm=0;
      }
    }    
    else {
      dm=0.;
    }        
  }
  else{
    dm=0;
  }

  if (isinf(dm)) {
    fprintf(stderr,"dm=%le m=%le\n",dm,m);
    fprintf(stderr,"t=%le tdis[1]=%le minIn=%d maxIn=%d Ndis=%d \n",t,tdis[1],minIn,maxIn,Ndis);
    for (count=minIn;count<=maxIn;count++)  fprintf(stderr,"dmdtvec[%d]=%le tofdmdtvec[%d]=%le\n",count,dmdtvec[count],count,tofdmdtvec[count]);
    fprintf(stderr,"terminating now...\n");
    exit(6);
  }

  if (m<=mf) dm=0;
  if (dm>=0) dm=0;

  return dm;
}
/*******************************************/
double dspindt(double t,double spin){
  double dspin=0.;
  double spinf;
  spinf=spinSSE[Nlin];
  t=t+tMS;

  if (!(SSE1)){
    if (t<=tofdmdtvec[2]){
      dspin=dspindtvec[1]*spin;
    }
    else if (t>tofdmdtvec[2] && t<=tofdmdtvec[(Ndis+1)*Nvar]){
      linint(tofdmdtvec,dspindtvec,(Ndis+1)*Nvar,t,&dspin);
      dspin=dspin*spin;
      if (t>tofdmdtvec[(Ndis+1)*Nvar-1]){
	dspin=dspindtvec[(Ndis+1)*Nvar-1]*spin;
      }
    }
    /*else if (t>tofdmdtvec[(Ndis+1)*Nvar]){
      if (fabs(spin-spinf)>1.e-6){    
	dspin=(spinf-spin)/t_inter;    
      }
      else {
	dspin=0;
      }
    }*/    
    else {
      dspin=0.;
    }        
  }
  else{
    dspin=0;
  }

  return dspin;
}


/*******************************************/

double dRdt(double t,double R){
  double dR,Rf1;
  int count,count2,minIn,maxIn;
  dR=0;
  minIn=1;
  maxIn=indexdis[Ndis]-1;
  t=t+tMS; 

  Rf1=pow(10.,logRSSE[Nlin])*Rsun;

  if (!(SSE1)){

    if (t<=tofdmdtvec[2]){
      dR=dRdtvec[1]*Rsun;
    }

    else if (t>tofdmdtvec[2] && t<=tofdmdtvec[(Ndis+1)*Nvar]){
      linint(tofdmdtvec,dRdtvec,(Ndis+1)*Nvar,t,&dR);
      dR=dR*Rsun;
      if (t>tofdmdtvec[(Ndis+1)*Nvar-1]){
	dR=dRdtvec[(Ndis+1)*Nvar-1]*Rsun;
      }
    }
    
    else if (t>tofdmdtvec[(Ndis+1)*Nvar]){
      if (R>Rf1){
	dR=-(R-Rf1)/(fabs(tofdmdtvec[(Ndis)*Nvar-1]-tofdmdtvec[(Ndis)*Nvar]));
      }
      else {
	dR=0;
      }
    }

  }
  else{
    dR=0;
  }

  if (isinf(dR)) {
    fprintf(stderr,"dR=%le R=%le\n",dR,R);
    fprintf(stderr,"t=%le tdis[1]=%le minIn=%d maxIn=%d Ndis=%d \n",t,tdis[1],minIn,maxIn,Ndis);
    for (count=minIn;count<=maxIn;count++)  fprintf(stderr,"dRdtvec[%d]=%le tofdmdtvec[%d]=%le\n",count,dRdtvec[count],count,tofdmdtvec[count]);
    fprintf(stderr,"terminating now...\n");
    exit(6);
  }

  if (R<=Rf1) dR=0;

  return dR;
}

/*******************************************/
double dmdt2(double t,double m, double m0){
  double m2f;
  double dm2=0.;
  double err=1.e-8;
  double Deltat;
  int indexval,min,max,length,count,j,iSSE,minIn,maxIn;
  t=t+tMS;

  indexval=indexdis2[Ndis2];
  minIn=1;
  maxIn=indexdis2[Ndis2];
  m2f=MtSSE2[Nlin2];

  if (!(SSE2)){
    if (t<=tofdmdtvec2[2]){
      dm2=dmdtvec2[1];
    }
    else if (t>tofdmdtvec2[2] && t<=tofdmdtvec2[(Ndis2+1)*Nvar]){
      linint(tofdmdtvec2,dmdtvec2,(Ndis2+1)*Nvar,t,&dm2);
      if (t>tofdmdtvec2[(Ndis2+1)*Nvar-1]){
	dm2=dmdtvec2[(Ndis2+1)*Nvar-1];
      }
    }
    else if (t>tofdmdtvec2[(Ndis2+1)*Nvar]){
      if (m>m2f){    
	dm2=(m2f-m)/t_inter;    
      }
      else {
	dm2=0;
      }
    }    
    else {
      dm2=0.;
    }        
  }
  else{
    dm2=0;
  }

  if (isinf(dm2)) {
    fprintf(stderr,"dm2=%le m2=%le\n",dm2,m);
    fprintf(stderr,"t=%le tdis2[1]=%le minIn=%d maxIn=%d Ndis2=%d \n",t,tdis2[1],minIn,maxIn,Ndis2);
    for (count=minIn;count<=maxIn;count++)  fprintf(stderr,"dmdtvec[%d]=%le tofdmdtvec[%d]=%le\n",count,dmdtvec2[count],count,tofdmdtvec2[count]);
    fprintf(stderr,"terminating now...\n");
    exit(6);
  }

  if (m<=m2f) dm2=0;
  if (dm2>=0) dm2=0;

  return dm2;
}

/*******************************************/
double dspindt2(double t,double spin2){
  double dspin2=0.;
  double spinf2;
  spinf2=spinSSE2[Nlin2];
  t=t+tMS;

  if (!(SSE2)){
    if (t<=tofdmdtvec2[2]){
      dspin2=dspindtvec2[1]*spin2;
    }
    else if (t>tofdmdtvec2[2] && t<=tofdmdtvec2[(Ndis2+1)*Nvar]){
      linint(tofdmdtvec2,dspindtvec2,(Ndis2+1)*Nvar,t,&dspin2);
      dspin2=dspin2*spin2;
      if (t>tofdmdtvec2[(Ndis2+1)*Nvar-1]){
	dspin2=dspindtvec2[(Ndis2+1)*Nvar-1]*spin2;
      }
    }
    /*else if (t>tofdmdtvec2[(Ndis2+1)*Nvar]){
      if (fabs(spin2-spinf2)>1.e-6){    
	dspin2=(spinf2-spin2)/t_inter;    
      }
      else {
	dspin2=0;
      }
      }*/    
    else {
      dspin2=0.;
    }        
  }
  else{
    dspin2=0;
  }

  return dspin2;
}


/*******************************************/

double dRdt2(double t,double R){
  double dR2,Rf2;
  int count,count2,minIn,maxIn;
  dR2=0;
  minIn=1;
  maxIn=indexdis2[Ndis2]-1;
  t=t+tMS; 

  Rf2=pow(10.,logRSSE2[Nlin2])*Rsun;

  if (!(SSE2)){

    if (t<=tofdmdtvec2[2]){
      dR2=dRdtvec2[1]*Rsun;
    }

    else if (t>tofdmdtvec2[2] && t<=tofdmdtvec2[(Ndis2+1)*Nvar]){
      linint(tofdmdtvec2,dRdtvec2,(Ndis2+1)*Nvar,t,&dR2);
      dR2=dR2*Rsun;
      if (t>tofdmdtvec2[(Ndis2+1)*Nvar-1]){
	dR2=dRdtvec2[(Ndis2+1)*Nvar-1]*Rsun;
      }
    }
    
    else if (t>tofdmdtvec2[(Ndis2+1)*Nvar]){
      if (R>Rf2){
	dR2=-(R-Rf2)/(fabs(tofdmdtvec2[(Ndis2)*Nvar-1]-tofdmdtvec2[(Ndis2)*Nvar]));
      }
      else {
	dR2=0;
      }
    }

  }
  else{
    dR2=0;
  }

  if (isinf(dR2)) {
    fprintf(stderr,"dR2=%le R2=%le\n",dR2,R);
    fprintf(stderr,"t=%le tdis2[1]=%le minIn=%d maxIn=%d Ndis2=%d \n",t,tdis2[1],minIn,maxIn,Ndis2);
    for (count=minIn;count<=maxIn;count++)  fprintf(stderr,"dRdtvec2[%d]=%le tofdmdtvec2[%d]=%le\n",count,dRdtvec2[count],count,tofdmdtvec2[count]);
    fprintf(stderr,"terminating now...\n");
    exit(6);
  }

  if (R<=Rf2) dR2=0;

  return dR2;
}