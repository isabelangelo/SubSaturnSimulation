/*==============================================*/
/* generate multiple files to be read later    */
/*using perl into octopule.                  */
/*==============================================*/
// complie gcc -O ICgen_StarsRL_uni.c ran1.c gasdev.c ran2.c nrutil.c -lm
#include <stdio.h> 
#include <stdlib.h> 
#include <string.h> 
#include <math.h>
#include "nrutil.h"
#define daytoyr (1./365.25)
#define Be (1./(0.33*0.33))
#define Rsun (695500*6.68459e-9)                                   /* sun radius au       */


int Ndismax = 100; /*maximum number of type changes allowed */
int Nlinmax = 5000; /*maximum number of lines from sse output */

void readssecl(double mass, int *Ndis, int *Nlin, double *tdis, double *timeSSE, double *typeSSE, double *MoSSE, double *MtSSE, double *logLSSE, double *logRSSE, double *logTSSE, double *McSSE, double *MenvSSE, double *epochSSE, double *spinSSE)
{

  FILE *fp;
  int status;
  char command[100], output[200];
  char *pch;
  int j,pos,iSSE;
  fprintf(stderr,"running SSE with mass=%le M_sun ...\n",mass);
  sprintf(command, "SSE/./sse_cl %g", mass);
  //  sprintf(command, "./sse %g", mass);
  //sprintf(command, "./bse %g", mass);

  /* run command */
  fp = popen(command, "r");
  if (fp == NULL) {
    printf("Failed to run command\n" );
    exit(1);
  }
 

  /* Read the output line-by-line and save the output in the arrays. */
 j=0;
  *Ndis = -1;
  *Nlin = 0;
  while (fgets(output, sizeof(output)-1, fp) != NULL) {
    if (j == 0) *Ndis = atof(output);
    if (j > 0 && j < *Ndis+1) tdis[j] = atof(output);
    if (*Ndis != -1 && j > *Ndis+1) {
      (*Nlin)++;
      pos = j-*Ndis-1;
      //  fprintf(stderr,"pos=%d *Ndis=%i j=%d size=%lu\n",pos,*Ndis,j,sizeof(output));
      //  exit(7);
      pch = strtok(output," ");
      timeSSE[pos] = atof(pch);
      //  fprintf(stderr,"timeSSE[%d]=%le\n",pos,timeSSE[pos]);
      pch = strtok(NULL," ");
      typeSSE[pos] = atof(pch);
      pch = strtok(NULL," ");
      MoSSE[pos] = atof(pch);
      pch = strtok(NULL," ");
      MtSSE[pos] = atof(pch);
      pch = strtok(NULL," ");
      logLSSE[pos] = atof(pch);
      pch = strtok(NULL," ");
      logRSSE[pos] = atof(pch);
      pch = strtok(NULL," ");
      logTSSE[pos] = atof(pch);
      pch = strtok(NULL," ");
      McSSE[pos] = atof(pch);
      pch = strtok(NULL," ");
      MenvSSE[pos] = atof(pch);
      pch = strtok(NULL," ");
      epochSSE[pos] = atof(pch);
      pch = strtok(NULL," ");
      spinSSE[pos] = atof(pch);
    }
    j++;
    /* safety checks */
    if (*Ndis > Ndismax) {
      fprintf(stderr,"WARNING: Ndis=%i > Ndismax=%i !!!\n",*Ndis,Ndismax);
      exit(0);
    }
   if (*Nlin > Nlinmax) {
      fprintf(stderr,"WARNING: Nlin=%i > Nlinmax=%i !!!\n",*Nlin,Nlinmax);
      exit(0);
    }
  }

  pclose(fp);
  
}

int countMAX=30;

double rL(double q){
  // Eggleton 83 equation                                                                                                                                                       
  // q is mdonner/mreciver     
  /// This Roche may be funcky use the one from my code - done                                                                                                                                                 
  double num1=0.49,num2=0.6;
  return (num1*pow(q,2./3.)/(num2* pow(q,2./3.)+log(1+pow(q,1./3.))));
}


double Mar(double q,double eout,double imu){
  //Mardling & Aarseth 2001 stability
  //qout=m3/(m1+m2)
  double val;
  val=2.8*pow((1.+q),(2./5.))*pow(1.+eout,2./5.)*pow(1.-eout,-6./5.)*(1.-0.3*imu/180);
  //  printf("q=%le eout=%le imu=%le\n",q,eout,imu);
 return val;
}

double Naoz(double q, double ein, double eout){
  double val;
  val = pow(3*q,1./3.)*(1+ein)/(1-eout);
  return val;
}

int Num=4000;//7000; //this is the max number of steps taken in order to generate the files
//int Num2=1000; //max number of desired files
double logPBAR=4.8;
double siglogP=2.3;
double q12mean=0.23;
double q12sig=0.42;
int Mmax=500;
int MoreN=1;//580;

double *tdis;
double *timeSSE,*typeSSE,*MoSSE,*MtSSE,*logLSSE,*logRSSE,*logTSSE,*McSSE,*MenvSSE,*epochSSE,*spinSSE,*dmdtvec,*tofdmdtvec,*dRdtvec;
double *tdis2;
double *timeSSE2,*typeSSE2,*MoSSE2,*MtSSE2,*logLSSE2,*logRSSE2,*logTSSE2,*McSSE2,*MenvSSE2,*epochSSE2,*spinSSE2,*dmdtvec2,*tofdmdtvec2,*dRdtvec2;
double *tdis3;
double *timeSSE3,*typeSSE3,*MoSSE3,*MtSSE3,*logLSSE3,*logRSSE3,*logTSSE3,*McSSE3,*MenvSSE3,*epochSSE3,*spinSSE3,*dmdtvec3,*tofdmdtvec3,*dRdtvec3;
int *indexdis,*indexdis2,*indexdis3; //this is the index of discontinutites
double type1,type2,type3;
int Ndis,Nlin;
int Ndis2,Nlin2;
int Ndis3,Nlin3;

main(){
  int i,j,counter=0,countm,o;
  FILE *out,*outF,*outTest;
  double m1,m2,m3,R1,R2,spin1,beta,beta2,gamma,gamma2,spin2,a1,a2,e1,e2,g1,g2,inc,age,t_ev,t_rr;
  double cosi,LittleOmega;
  double mine2,maxe2,minm3,maxm3,mincosi,maxcosi;  
  double Marval,Naozval,Roche1,Roche2;
  double logP,logPval,P,P2,P1,PT1,PT2;
  double q12,qout;
  double alpha,kappa;
  double epsilon,a1GR,a1high;
  double k2=4*M_PI*M_PI;
  double c2=3.993953869e9;
  char pathname[500],name[500],failname[500];
  long idum=-13;
  /* flags:  0 yes 1 no  */
  int QUADRUPOLE=0,OCTUPOLE=0,GR=0,TF=0,ML1=0,MB1=0,ML2=1,MB2=1,SSE1=0,SSE2=1,SSE3=0,tMSMyr=0.0; 
  float ran1(long *);  
  float ran2(long *);  
  float gasdev(long *);  
  float alpha2, m1max, m1min, a2max, a2min, P2min, P2max;
  int iSSE,count;
  double tSN1,tSN2,tSN3;

  double m_store;

  tdis = vector(1,Ndismax);
  indexdis = ivector(1,Ndismax);

  timeSSE = vector(1,Nlinmax);
  typeSSE = vector(1,Nlinmax);
  MoSSE = vector(1,Nlinmax);
  MtSSE = vector(1,Nlinmax);
  logLSSE = vector(1,Nlinmax);
  logRSSE = vector(1,Nlinmax);
  logTSSE = vector(1,Nlinmax);
  McSSE = vector(1,Nlinmax);
  MenvSSE = vector(1,Nlinmax);
  epochSSE = vector(1,Nlinmax);
  spinSSE = vector(1,Nlinmax);

  tdis2 = vector(1,Ndismax);
  indexdis2 = ivector(1,Ndismax);

  timeSSE2 = vector(1,Nlinmax);
  typeSSE2 = vector(1,Nlinmax);
  MoSSE2 = vector(1,Nlinmax);
  MtSSE2 = vector(1,Nlinmax);
  logLSSE2 = vector(1,Nlinmax);
  logRSSE2 = vector(1,Nlinmax);
  logTSSE2 = vector(1,Nlinmax);
  McSSE2 = vector(1,Nlinmax);
  MenvSSE2 = vector(1,Nlinmax);
  epochSSE2 = vector(1,Nlinmax);
  spinSSE2 = vector(1,Nlinmax);

  tdis3 = vector(1,Ndismax);
  indexdis3 = ivector(1,Ndismax);

  timeSSE3 = vector(1,Nlinmax);
  typeSSE3 = vector(1,Nlinmax);
  MoSSE3 = vector(1,Nlinmax);
  MtSSE3 = vector(1,Nlinmax);
  logLSSE3 = vector(1,Nlinmax);
  logRSSE3 = vector(1,Nlinmax);
  logTSSE3 = vector(1,Nlinmax);
  McSSE3 = vector(1,Nlinmax);
  MenvSSE3 = vector(1,Nlinmax);
  epochSSE3 = vector(1,Nlinmax);
  spinSSE3 = vector(1,Nlinmax);


  printf(" ====== Shalom, this code generate initial conditions files For Stars! ======\n");
  printf(" ====== The file names are triple.in##                      ======\n");

  alpha = 2.35;//1.7;//2.35;
  alpha2 = 0.4;
  m1min = 1.;//1.6;//1.;
  m1max = 1.6;//2.4;//150.;
  a2min = 700.;//700.;
  a2max = 20650.;
  kappa = -0.55;//power law parameter
  /******/
  /** set some parameters  **/
  //  m1=1.;//Msun
  //R1=1.;
  // m2=0.001; //for now, set it to be jpiter
  //R2=0.10012;// RJ
  //spin1=10.;
  //spin2=10.;
  /*
  g1=45.;
  g2=0.;
  */
  beta=0.;
  //a1=5.;//6.;
  //a2=500.; //randomize - see below
  //  e1=0.01; //start in situ of the disk
  
  //age=10000;//Myr
  /******/
  /** min max for the random varibels **/
  mine2=0.;
  maxe2=1.;
  mincosi=-1.;//cos(0.*M_PI/180.);
  maxcosi=1.;//cos(180.*M_PI/180.);
  mincosi=cos(40.*M_PI/180.);
  maxcosi=cos(140.*M_PI/180.);

sprintf(pathname,"/u/home/a/alexpste/IniConWD2020_focus/");
  i=1;
  j=1;
  o=1;
  i=MoreN;
  j=MoreN;
  Num=Num+MoreN;
  sprintf(failname,"%sfailedStarsuni.txt",pathname);
  outF=fopen(failname,"w");
  outTest=fopen("ICsStarsRL_uni.txt","w"); 

  //for (i=1;i<=Num;i++){
  printf(" generating ...\n");
  while (j<=Num){
    //printf("You are here 1 \n");
    // sprintf(name,"%striple.in%3.3d",pathname,i);
    for (countm=0;countm<=Mmax;countm++){
      //m1 = pow((pow(m1max,1-alpha)-pow(m1min,1-alpha))*ran2(&idum)+pow(m1min,1-alpha),1./(1-alpha)); // Salpeter dist. with alpha = 1.7//2.35
      m1 = m1min+(m1max-m1min)*ran2(&idum);
      if (m1>=1.0) break;
      else if (m1<1.0 && countm==Mmax){
	printf("oh no, you have tried %d times to find a posotive mass (m1) and faild\n");
        printf("Some thing is wrong...\n");
        printf("Terminating now\n");
        exit(5);
      }
    }
    // Choose paramters for the inner orbit:
    for (countm=0;countm<=Mmax;countm++){
      //q12=gasdev(&idum);
      //m2=pow((pow(m1max,1-alpha)-pow(m1min,1-alpha))*ran2(&idum)+pow(m1min,1-alpha),1./(1-alpha)); // Salpeter dist. with alpha = 1.7//2.35  
      // m2=m1*(q12*q12sig+q12mean);
      q12=0.1+0.9*ran2(&idum);//pow((1.-pow(0.01,1.-alpha2))*ran2(&idum)+pow(0.01,1.-alpha2),1./(1.-alpha2));
      // m2 =m1/q12;
      //m2 = m1min+(m1max-m1min)*ran2(&idum);
      if (ran2(&idum)>=0.5){
	m2 = 0.0009546;//0.0009543;//0.00005149;//8.345e-9;
	R2 = 0.1;
	spin2 = 9.9/24.;
      }
      else{
	m2 = 0.00005149;
	R2 = 0.0354;
	spin2 = 16.11/24.;
      }
      //printf("You are here 2 \n");
      if (m2>0.) break;//(m2>=1. && m2<=150.) break;
      else if (m2<=0. && countm==Mmax){
	printf("oh no, you have tried %d times to find a posotive mass (m2) and faild\n");
      //printf("Some thing is wrong...\n");
      //printf("Terminating now\n");
      //exit(5);
      }
    }

    for(countm=0;countm<=Mmax;countm++){
      qout=gasdev(&idum);
      //m3 = 1e-8;
      m3=(m1)*(qout*q12sig+q12mean);                                                                                                                        
      //m3 = 2.*m2;
      //m3=0.75;
      /*if(m3>=1.){
	q12 = ran2(&idum);
	if(q12 >= 0.5){
	  m_store = m1;
	  m1 = m3;
	  m3 = m_store;
	  }
      }*/                             
      //m3=4.0e6;                                                                                                                                                                               
      //m3 = pow((pow(8.,1-alpha)-pow(0.1,1-alpha))*ran2(&idum)+pow(0.1,1-alpha),1./(1-alpha)); // Salpeter dist. with alpha = 1.7//2.35                                                        
      if (m3>0.1) break;
      else if (m3<=0.1 &&countm==Mmax){
        printf("oh no, you have tried %d times to find a posotive mass (m3)  and faild\n");
        printf("Some thing is wrong...\n");
        printf("Terminating now\n");
        exit(5);
      }
    }

  
  for (iSSE=1; iSSE<=Ndismax; iSSE++) tdis[iSSE] = -1.;
  for (iSSE=1; iSSE<=Nlinmax; iSSE++){
    timeSSE[iSSE] = -1;
    typeSSE[iSSE] = -1;
    MoSSE[iSSE] = -1;
    MtSSE[iSSE] = -1;
    logLSSE[iSSE] = -1;
    logRSSE[iSSE] = -1;
    logTSSE[iSSE] = -1;
    McSSE[iSSE] = -1;
    MenvSSE[iSSE] = -1;
    epochSSE[iSSE] = -1;
    spinSSE[iSSE] = -1;
  }
  readssecl(m1, &Ndis, &Nlin, tdis, timeSSE, typeSSE, MoSSE, MtSSE, logLSSE, logRSSE, logTSSE, McSSE, MenvSSE, epochSSE, spinSSE);
  spin1=spinSSE[1];
  R1=pow(10.,logRSSE[1]);

  for (iSSE=1; iSSE<=Ndismax; iSSE++) tdis2[iSSE] = -1.;
  for (iSSE=1; iSSE<=Nlinmax; iSSE++){
    timeSSE2[iSSE] = -1;
    typeSSE2[iSSE] = -1;
    MoSSE2[iSSE] = -1;
    MtSSE2[iSSE] = -1;
    logLSSE2[iSSE] = -1;
    logRSSE2[iSSE] = -1;
    logTSSE2[iSSE] = -1;
    McSSE2[iSSE] = -1;
    MenvSSE2[iSSE] = -1;
    epochSSE2[iSSE] = -1;
    spinSSE2[iSSE] = -1;
  }
  //readssecl(m2, &Ndis2, &Nlin2, tdis2, timeSSE2, typeSSE2, MoSSE2, MtSSE2, logLSSE2, logRSSE2, logTSSE2, McSSE2, MenvSSE2, epochSSE2, spinSSE2);
  //spin2=spinSSE2[1];
  //R2=pow(10.,logRSSE2[1]);
  //spin2 = 10./24.;//1.;//10./24.;
  //R2 = 0.1;//0.00167;//0.1;

  for (iSSE=1; iSSE<=Ndismax; iSSE++) tdis3[iSSE] = -1.;
  for (iSSE=1; iSSE<=Nlinmax; iSSE++){
    timeSSE3[iSSE] = -1;
    typeSSE3[iSSE] = -1;
    MoSSE3[iSSE] = -1;
    MtSSE3[iSSE] = -1;
    logLSSE3[iSSE] = -1;
    logRSSE3[iSSE] = -1;
    logTSSE3[iSSE] = -1;
    McSSE3[iSSE] = -1;
    MenvSSE3[iSSE] = -1;
    epochSSE3[iSSE] = -1;
    spinSSE3[iSSE] = -1;
  }
  readssecl(m3, &Ndis3, &Nlin3, tdis3, timeSSE3, typeSSE3, MoSSE3, MtSSE3, logLSSE3, logRSSE3, logTSSE3, McSSE3, MenvSSE3, epochSSE3, spinSSE3);

  count=1;
  iSSE=1;
  while (timeSSE[count]<=tdis[Ndis]){
    // fprintf(stderr,"timesSSE[%d]=%le tdis[%d]=%le timeSSE[%d +1]=%le \n",count,timeSSE[count],iSSE,tdis[iSSE],count,timeSSE[count+1]);
    // if (fabs((timeSSE[count]-tdis[iSSE])/tdis[iSSE])<6.2e-7){
    if (timeSSE[count]<=tdis[iSSE] && timeSSE[count+1]>tdis[iSSE]){
      //  fprintf(stderr,"==times[%d]=%le tdis[%d]=%le  \n",count,timeSSE[count],iSSE,tdis[iSSE]);
      indexdis[iSSE]=count;
      iSSE++;
      count++;
    }
    else{
      count++;
    }
  }

  /*count=1;
  iSSE=1;
  while (timeSSE2[count]<=tdis2[Ndis2]){
    // fprintf(stderr,"timesSSE[%d]=%le tdis[%d]=%le timeSSE[%d +1]=%le \n",count,timeSSE[count],iSSE,tdis[iSSE],count,timeSSE[count+1]);                                                       
    // if (fabs((timeSSE[count]-tdis[iSSE])/tdis[iSSE])<6.2e-7){                                                                                                                                
    if (timeSSE2[count]<=tdis2[iSSE] && timeSSE2[count+1]>tdis2[iSSE]){
      //  fprintf(stderr,"==times[%d]=%le tdis[%d]=%le  \n",count,timeSSE[count],iSSE,tdis[iSSE]);                                                                                              
      indexdis2[iSSE]=count;
      iSSE++;
      count++;
    }
    else{
      count++;
    }
    }*/

  count=1;
  iSSE=1;
  while (timeSSE3[count]<=tdis3[Ndis3]){
    // fprintf(stderr,"timesSSE[%d]=%le tdis[%d]=%le timeSSE[%d +1]=%le \n",count,timeSSE[count],iSSE,tdis[iSSE],count,timeSSE[count+1]);                                                       
    // if (fabs((timeSSE[count]-tdis[iSSE])/tdis[iSSE])<6.2e-7){                                                                                                                                
    if (timeSSE3[count]<=tdis3[iSSE] && timeSSE3[count+1]>tdis3[iSSE]){
      //  fprintf(stderr,"==times[%d]=%le tdis[%d]=%le  \n",count,timeSSE[count],iSSE,tdis[iSSE]);                                                                                              
      indexdis3[iSSE]=count;
      iSSE++;
      count++;
    }
    else{
      count++;
    }
    }
  //printf("I am here(1)\n");
  tSN1=timeSSE[indexdis[Ndis]-1];
  //printf("runtime=%le\n",tSN1);
  //printf("I am here(2)\n");
  tSN2=timeSSE2[indexdis2[Ndis2]-1];
  //printf("I am here(3)\n");
  tSN3=timeSSE3[indexdis3[Ndis3]-1];
  // printf("I am here(1)\n");

    //R2=R1*pow(m2/m1,0.8);


    //assuming uniform in log dis since we only start from tail of  Duquennoy & Mayor log-normal, so we approximate as const.
  /*for(countm=0;countm<=Mmax;countm++){

    logPval=gasdev(&idum);
    //logP=logPval*siglogP+logPBAR; //this is the log of the period of the inner orbit as a result of normal dis with mean logPBAR and sigma=siglogP
    if (m2 <= 15. && m1 <= 15.){
      PT1=pow(10,0.15+(3.5-0.15)*ran2(&idum))*daytoyr;
    }
    else if (m2 > 15. || m1 > 15.){
      PT1=pow(10,pow((pow(3.5,1+kappa)-pow(0.15,1+kappa))*ran2(&idum)+pow(0.15,1+kappa),1/(1+kappa)))*daytoyr;
    }
    if (log10(PT1/daytoyr)>=0.15 && log10(PT1/daytoyr)<=3.5) break;
    else if (log10(PT1/daytoyr)<0.15 || log10(PT1/daytoyr)>3.5 && countm==Mmax){
      printf("oh no, you have tried %d times to find a valid inner Period (P1)  and faild\n");
      printf("Some thing is wrong...\n");
      printf("Terminating now\n");
      exit(5);
    }
    }*/
  for(countm=0;countm<=Mmax;countm++){

                    
    /*if ( (m1/m2 > 0.1 && m1/m2 < 10.) || (m1<15. && m2>15.)){
      //PT1=pow(10,0.15+(3.5-0.15)*ran2(&idum))*daytoyr;
      a1=pow(10,-1.+(3.+1.)*ran2(&idum));
    }
    else if ((m1/m2 < 0.1 || m1/m2 > 10.) && (m1>15. || m2>15.)){
      PT1=pow(10,pow((pow(3.5,1+kappa)-pow(0.15,1+kappa))*ran2(&idum)+pow(0.15,1+kappa),1/(1+kappa)))*daytoyr;
      P1 = PT1;
      a1=pow(SQR(P1)*(m1+m2),1./3.);
      }*/

    //PT1 = pow(10.,3.65*ran2(&idum));
    //P1 = PT1*daytoyr;
    //a1=pow(SQR(P1)*(m1+m2),1./3.);  

    //a1 = 0.02+(1.0-0.02)*ran2(&idum);
    e1 = 0.01;//0.15*ran2(&idum);
    a1 = 3.+(100.-3.)*ran2(&idum);//20.0+(50.0-20.0)*ran2(&idum);//40.0+(100.0-40.0)*ran2(&idum);
    Roche1=rL(m1/m2);
    Roche2=rL(m2/m1);

    if (R1*Rsun<a1*(1-e1)*Roche1 && R2*Rsun<a1*(1-e1)*Roche2){
      break;
    }
  }
  
    //logPval=gasdev(&idum);
    //logP=logPval*siglogP+logPBAR; //this is the log of the period of the outer orbit as a result of normal dis with mean logP2BAR and sigma=siglogP2
    // PT2=exp(logP)*daytoyr;
    //PT2=pow(10.,logP)*daytoyr;

    /*
    P1=PT1;//DMIN(PT1,PT2);
    P2=PT2;//DMAX(PT1,PT2);
    */

    //P1=DMIN(PT1,PT2);
    //P2=DMAX(PT1,PT2);
    //printf("%lf \n", PT1);
    //P1 = PT1;
    
    //a1=pow(SQR(P1)*(m1+m2),1./3.); 
    //printf("a1=%le\n", a1);
    
    // printf("P1=%le yr, P1=%le day, a1=%le\n",P1,P1*365.25,a1);
    /*
    if (P1*365.25>1000){

      e1=sqrt(ran2(&idum));//thermal
      // printf("thermal : e1=%le\n",e1);
    }
    else{//Rayleigh
      // e1=sqrt(-log(2.*Be*ran2(&idum))/Be);
      e1=sqrt(-log(ran2(&idum))/Be);
      //printf(" Rayleigh : 2.*Be*ran2(&idum))/Be=%le e1=%le\n",ran2(&idum),e1);
      //printf(" Be=%le\n",Be);
    }
    */
    //e1=(maxe2-mine2)*ran2(&idum)+mine2; //uniform
  //e1 = 0.01;
    /* draw a random nuber for the relevent variables  */

    // see Payne et al 2011 (cof. procceding) why we choose a uniform dis. 
    cosi = (maxcosi-mincosi)*ran2(&idum)+mincosi ;//uniform
    inc = acos(cosi)*180/M_PI;


    // now for the outer orbit
   

    /*for(countm=0;countm<=Mmax;countm++){
      qout=gasdev(&idum);   
      //m3=(m1+m2)*(qout*q12sig+q12mean);
      //m3=4.0e6;
      m3 = pow((pow(100.,1-alpha)-pow(0.1,1-alpha))*ran2(&idum)+pow(0.1,1-alpha),1./(1-alpha)); // Salpeter dist. with alpha = 1.7//2.35
      if (m3>0) break;
      else if (m3<=0 &&countm==Mmax){
        printf("oh no, you have tried %d times to find a posotive mass (m3)  and faild\n");
        printf("Some thing is wrong...\n");
        printf("Terminating now\n");
        exit(5);
      }
      }*/

    //printf("You are here 4 \n");
    //P2min = sqrt(pow(a2min,3.)/(m1+m2+m3));
    //P2max = sqrt(pow(a2max,3.)/(m1+m2+m3));
    //P2=pow(10,log10(P2min)+(log10(P2max)-log10(P2min))*ran2(&idum)); // uniform in log   
    for(countm=0;countm<=Mmax;countm++){
      /*if ((m3/(m1+m2)>0.1 && m3/(m1+m2)<10.) || (m3 < 15. && m2 < 15. && m3 < 15.)){
	//P2=pow(10,0.15+(3.5-0.15)*ran2(&idum))*daytoyr;
	a2=pow(10,log10(a1)+(4.-log10(a1))*ran2(&idum));
      }
	else if ((m3/(m1+m2)<0.1 || m3/(m1+m2)>10. ) && (m3 > 15. || m1 > 15. || m2 > 15.)){
	P2=pow(10,pow((pow(3.5,1+kappa)-pow(0.15,1+kappa))*ran2(&idum)+pow(0.15,1+kappa),1/(1+kappa)))*daytoyr;
	a2=pow(SQR(P2)*(m1+m2+m3),1./3.);
	}*/
      logPval=gasdev(&idum);
      logP=logPval*siglogP+logPBAR;
      P2=pow(10.,logP)*daytoyr;
      a2=pow(SQR(P2)*(m1+m2+m3),1./3.);
      //a2=(100.-10.)*ran2(&idum)+10.;
      e2=(maxe2-mine2)*ran2(&idum)+mine2;
      epsilon = (a1/a2) *e2/(1.-(e2*e2));
      //Marval=Mar(m3/(m1+m2),e2,inc);
      if (epsilon<0.1 && a2<10000. && (a1/a2)<0.1 ) break;
      else if ((epsilon >= 0.1 || a2 > 10000 || (a1/a2)>=0.1) && countm==Mmax){
	printf("couldn't find stable system\n");
	//printf("a1=%le R1=%le a2=%le R2=%le\n",a1,R1,a2,R2);
      }
    }
    //a2=pow(SQR(P2)*(m1+m2+m3),1./3.); 
    //printf("a2=%le\n", a2);
    //printf("You are here 5 \n");
    //e2=(maxe2-mine2)*ran2(&idum)+mine2; //uniform 
    //e2 = sqrt(ran2(&idum));//thermal

    //a1GR=pow(pow(a2,3.)*3*k2*m1*m2*sqrt(1-e2*e2)/(c2*m3),1./4.);

    //a1high=0.1*a2*(1-e2*e2)/e2;

    //a1high=a2/Naoz(m3/(m1+m2),e1,e2);

    //a1=(a1high-a1GR)*ran2(&idum)+a1GR;
    /*
    if (P2*365.25>1000){

      e2=sqrt(ran2(&idum));//thermal

    }
    else{
      e2=sqrt(-log(ran2(&idum))/Be);
      // e2=sqrt(-log(2.*Be*ran2(&idum))/Be);
    }
    */

    //e2=(maxe2-mine2)*ran2(&idum)+mine2; //uniform
   
    //g1, g2 uniform in [0,2pi]
    g1=360.*ran2(&idum);
    g2=360.*ran2(&idum);
    beta=0.;//360.*ran2(&idum);
    beta2=0.;
    gamma=45.;
    gamma2=45.;

    Marval=Mar(m3/(m1+m2),e2,inc);
    Naozval=Naoz(m3/(m1+m2),e1,e2);

    Roche1=rL(m1/m2);
    Roche2=rL(m2/m1);

    //age = 2.0*sqrt(k2*m3/a2)/(16.0*sqrt(M_PI)*log(10.0)*a1*1.4*pow(10.0,-7.0)*pow(a2/82500.0,-1.4))*1e-6;
   
    //evaporation time scale derived by Alexander from Binney and Tremaine
    t_ev=sqrt(3.0/M_PI)*(280.*0.211*sqrt(20650.0/a2)/(32*k2*a1*log(15.)*5.2e5*pow(a2/20650.0*0.5,-1.8)/pow(206500,3)))*(m1+m2);

    //resonant relaxation time scale from Smadar

    /*if (t_ev < t_rr){
        age = t_ev*1e-6;
    }
    else {
	age = t_rr*1e-6;
    }*/
    //age = t_ev*1e-6;
    //if (age > 10){ // for high mass sample runs
    //  age=10;
    //}
    
    /*if (m1 >= m2 && m1 >= m3){
      age=tSN1;
    }
    else if (m2 >= m1 && m2 >= m3){
      age=tSN2;
    }
    else if (m1 < m3 && m2 < m3){
      age=tSN3;
      }*/
    age = 13000.;
    // printf("runtime=%le\n",age);
    //  printf("Marval=%le e2=%le inc=%le a2/a1=%le\n",Marval,e2,inc,a2/a1);
    // Add this condtion epsilon= a1/a2 *e2/(1-e2^2)<0.1
    epsilon = (a1/a2) *e2/(1-e2*e2);
    //printf("You are here 6 \n");
    if (a2/a1<Marval){
      o++;
    }
     if (R1*Rsun<a1*(1-e1)*Roche1 && R2*Rsun<a1*(1-e1)*Roche2 && epsilon<0.1){
      //printf("You are here 7 \n");
    //if (a2/a1>Marval  ){
  fprintf(outTest,"%le %le %le %le  %le %le %le %le %le %le %le  %le %le %le\n",a1,a2,P1,P2,e1,e2,m1,m2,m3,R1,R2,g1,g2,inc);	
    
	 
      sprintf(name,"%striple.in%d",pathname,i);
      out=fopen(name,"w");
      
      fprintf(out,"#########################################################################\n");
      fprintf(out,"## Read in the initial orbital parameters of the triple system         ##\n");
      fprintf(out,"## Mass in M_sol, semimajor axis in AU                                 ##\n");
      fprintf(out,"#########################################################################\n");
      fprintf(out,"\n");
      fprintf(out,"__m1_______m2___________m3________R1(Rsun)____R2(Rsun)____Spin1P(day)__Spin2P(day)__beta(s_O_deg)__beta2(s_O_deg)__gamma(s_O_deg)__gamma2(s_O_deg)__a1______a2____e1____e2____g1(deg)__g2(deg)__i(deg)__age(Myr)___ :::\n");
      fprintf(out,"\n");
      
      
      fprintf(out,"%le  %le  %le  %le  %le %le  %le  %le  %le  %le %le %le %le %le  %le  %le  %le  %le %le\n",m1,m2,m3,R1,R2,spin1,spin2,beta,beta2,gamma,gamma2,a1,a2,e1,e2,g1,g2,inc,age);
      fprintf(out,"\n");
      fprintf(out,"#########################################################################\n");
      fprintf(out,"## Flags                                                               ##\n");
      fprintf(out,"#########################################################################\n");
      fprintf(out,"\n");
      fprintf(out,"__QUADRUPOLE(0_yes_1_no)___ :::\n");
      fprintf(out,"%d\n",QUADRUPOLE);
      fprintf(out,"\n");
      fprintf(out,"#########################################################################\n");
      fprintf(out,"__OCTUPOLE(0_yes_1_no)___ :::\n");
      fprintf(out,"\n");
      fprintf(out,"%d\n",OCTUPOLE);
      fprintf(out,"\n");
      fprintf(out,"#########################################################################\n");
      fprintf(out,"__GR(0_yes_1_no)___ :::\n");
      fprintf(out,"\n");
      fprintf(out,"%d\n",GR);
      fprintf(out,"\n");
      fprintf(out,"#########################################################################\n");
      fprintf(out,"__TF(0_yes_1_no)___ :::\n");
      fprintf(out,"\n");
      fprintf(out,"%d\n",TF);
      fprintf(out,"\n");
      fprintf(out,"########################################################################\n");
      fprintf(out,"__ML1(0_yes_1_no)___ :::\n");
      fprintf(out,"\n");
      fprintf(out,"%d\n",ML1);
      fprintf(out,"\n");
      fprintf(out,"########################################################################\n");
      fprintf(out,"__MB1(0_yes_1_no)___ :::\n");
      fprintf(out,"\n");
      fprintf(out,"%d\n",MB1);
      fprintf(out,"\n");
      fprintf(out,"########################################################################\n");
      fprintf(out,"__ML2(0_yes_1_no)___ :::\n");
      fprintf(out,"\n");
      fprintf(out,"%d\n",ML2);
      fprintf(out,"\n");
      fprintf(out,"########################################################################\n");
      fprintf(out,"__MB2(0_yes_1_no)___ :::\n");
      fprintf(out,"\n");
      fprintf(out,"%d\n",MB2);
      fprintf(out,"\n");
      fprintf(out,"########################################################################\n");
      fprintf(out,"__SSE1(0_yes_1_no)___ :::\n");
      fprintf(out,"\n");
      fprintf(out,"%d\n",SSE1);
      fprintf(out,"\n");
      fprintf(out,"########################################################################\n");
      fprintf(out,"__SSE2(0_yes_1_no)___ :::\n");
      fprintf(out,"\n");
      fprintf(out,"%d\n",SSE2);
      fprintf(out,"\n");
      fprintf(out,"########################################################################\n");
      fprintf(out,"__SSE3(0_yes_1_no)___ :::\n");
      fprintf(out,"\n");
      fprintf(out,"%d\n",SSE3);
      fprintf(out,"\n");
      fprintf(out,"########################################################################\n");
      fprintf(out,"__tMSMyr(0_yes_1_no)___ :::\n");
      fprintf(out,"\n");
      fprintf(out,"%d\n",tMSMyr);
      fprintf(out,"\n");
      fprintf(out,"#########################################################################\n");
      fprintf(out,"## Control Parameters                                                  ##\n");
      fprintf(out,"#########################################################################\n");
      fprintf(out,"\n");
      fprintf(out,"__eps___ :::\n");
      fprintf(out,"1e-16\n");
      
      fclose(out);
      i++;
      }
      else {
	fprintf(outF,"a2/a1=%le Mar=%le %le  %le  %le  %le  %le %le  %le  %le  %le  %le %le  %le  %le  %le  %le %le\n",a2/a1,Marval,m1,m2,m3,R1,R2,spin1,spin2,beta,a1,a2,e1,e2,g1,g2,inc,age);
	counter++;
      }
    j++;  
  }
 
  fprintf(outF,"out of %d  there where %d un-stabel ones\n",j-1,counter);
  fprintf(outF,"%d violated Mard. stability criterion\n",o-1);
  printf("   The end\n");

}
