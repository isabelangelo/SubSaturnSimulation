/****************************************************************************
 * Integrate octupole perturbation equations for heirarchical triple systems 
 * Takes input from "triple.in" 
 * * smadar: adding TF May 2010   ** smadar adding MAss loss and magnetic breaking Aug- sept 2014   
 ****************************************************************************/


#include <stdio.h> 
#include <stdlib.h> 
#include <string.h> 
#include <math.h>
#include "nrutil.h"
#include "HTI.h"
#include <time.h>

#define TINY 1e-30


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
  /*
 	fprintf(stderr,"%d %d\n",*Ndis,*Nlin);
	//	exit(3);
for (iSSE=0; iSSE<*Ndis; iSSE++) fprintf(stderr,"%f\n",tdis[iSSE]); 
 for (iSSE=0; iSSE<*Nlin; iSSE++) {
   fprintf(stderr,"%f %f %f %f %f %f %f %f %f %f %le\n",timeSSE[iSSE],typeSSE[iSSE],MoSSE[iSSE],MtSSE[iSSE],logLSSE[iSSE],logRSSE[iSSE],logTSSE[iSSE],McSSE[iSSE],MenvSSE[iSSE],epochSSE[iSSE],spinSSE[iSSE]); 
   exit(4);
 }
  exit(2);
  */
}

/*********************************/
double Roche(double r,double q){
  //Roche=R_1 /(0.462 q^(1/3)  ) , where q=m_1/m2
  double RL,factor;
  //factor = 0.49*pow(q, 2./3.)/(0.6* pow(q,2./3.) + log(1. + pow(q,1./3.))); // Eggleton for 2 stars
  factor=0.6*pow(q,1./3.);  // check Roche after running
  RL = r/factor;
  
 //  RL=r/(0.462 * pow( q,1./3.  ));
  return RL;
  //r/(0.462 * pow( q, 1./3.  )); 
}

/*********************************/
/* odeint control parameters */
int kmax = 10000000,    kount = 0;
double *v, *xp, **yp, dxsav,**vp,*ve;						/* storage arrays */
//smadar: adding a storage matrix vp for v parameters 
int flgout;

double *tdis;
double *timeSSE,*typeSSE,*MoSSE,*MtSSE,*logLSSE,*logRSSE,*logTSSE,*McSSE,*MenvSSE,*epochSSE,*spinSSE,*dmdtvec,*tofdmdtvec,*dRdtvec,*mvec,*Rvec,*spinvec,*dspindtvec;
double *tdis2;
double *timeSSE2,*typeSSE2,*MoSSE2,*MtSSE2,*logLSSE2,*logRSSE2,*logTSSE2,*McSSE2,*MenvSSE2,*epochSSE2,*spinSSE2,*dmdtvec2,*tofdmdtvec2,*dRdtvec2,*mvec2,*Rvec2,*spinvec2,*dspindtvec2;
double *tdis3;
double *timeSSE3,*typeSSE3,*MoSSE3,*MtSSE3,*logLSSE3,*logRSSE3,*logTSSE3,*McSSE3,*MenvSSE3,*epochSSE3,*spinSSE3;
int *indexdis,*indexdis2,*indexdis3; //this is the index of discontinutites
double type1,type2,type3;
int Ndis,Nlin;
int Ndis2,Nlin2;
int Ndis3,Nlin3;
int Nvar;
double R1f,R2f;
double dR1dt,dR2dt;
double t_inter;
double *tmidDis,*tmidDis2;
double *toutput;
double *detotdtTF,*dhtotdtTF;
void linintMinMax(double *xa,double *ya,int minv,int n,double x,double *y);
/*********************************/


void findt1(double t,int *minIn,int *maxIn){

  int count;
  *minIn=1;
  *maxIn=Nlin;
  if (t<tdis[1]){
       *minIn=1;
       *maxIn=indexdis[1]-2;
  }
  else if (t>=tdis[Ndis]){
    *minIn=indexdis[Ndis-1];
    *maxIn=Nlin;
       
  }
  else {
     
    for (count=2;count<=Ndis;count++){
      if (t>=tdis[count-1] && t<tdis[count]){
	*minIn=indexdis[count-1]+2;
	*maxIn=indexdis[count]-2;
	break;
      }
    }
  }
  //fprintf(stderr, " ->:: minIn=%d maxIn=%d\n",*minIn,*maxIn);
  // return;
}

/*********************************/

void findt2(double t,int *minIn,int *maxIn){

  int count;
  *minIn=1;
  *maxIn=Nlin2;
  if (t<tdis2[1]){
       *minIn=1;
       *maxIn=indexdis2[1]-2;
  }
  else if (t>=tdis2[Ndis2]){
    *minIn=indexdis2[Ndis2-1];
    *maxIn=Nlin2;
       
  }
  else {
     
    for (count=2;count<=Ndis2;count++){
      if (t>=tdis2[count-1] && t<tdis2[count]){
	*minIn=indexdis2[count-1]+2;
	*maxIn=indexdis2[count]-2;
	break;
      }
    }
  }
  //fprintf(stderr, " ->:: minIn=%d maxIn=%d\n",*minIn,*maxIn);
  // return;
}

/*********************************/


void findt3(double t,int *minIn,int *maxIn){

  int count;
 *minIn=1;
  *maxIn=Nlin3;
  if (t<tdis3[1]){
       *minIn=1;
       *maxIn=indexdis3[1]-2;
  }
  else if (t>=tdis3[Ndis3]){
    *minIn=indexdis3[Ndis3-1];
    *maxIn=Nlin3;
       
  }
  else {
     
    for (count=2;count<=Ndis3;count++){
      if (t>=tdis3[count-1] && t<tdis3[count]){
	*minIn=indexdis3[count-1]+2;
	*maxIn=indexdis3[count]-2;
	break;
      }
    }
  }
  // return;
}

/*********************************/
/* main program, calls for odeint.c */
//int main (void) {
int main (int argc,char **argv) {	
	double 	*x, *x_i, *dxdt;
	double 	t_i, t_f, eps=0, h1=0, hmin=0, t_ini=0.;
	int		nok=0, nbad=0, j=0, k_f=0;
	int		scnret;
	int iSSE,count,min,max,length,count2;
	int minIn,maxIn;
	char	dummy[1024];
	double	m1,m2,m3,a1,a2,e1,e2,i,g1,g2,e1max_calc;

	double	P1,P2,P_koz,P_GR,range;
	FILE*	input;
	FILE *out;//smadar - save the max and min values
	float ideg;  //smadar
	double imax=0.,imin=360.,e1max=0,e1min=1.,e2max=0,e2min=1.,imumax=0.,imumin=360.,tempimin=360.,tempimax=0.;//smadar: what is the min and max inclinations reached? 
	double betamax=0,betamin=360;
	double Sp1=365.25,Sp2=365.25;//smadar:spin period in days!
	double i1,i2,imu,spin1,spin2,theta1,theta2;//smadar: adding the spin, for now - change later!
	double L1,L2,G1,G2,m12,m123,Hsq,t,theta,H1,H2;
	char outputname[700];
	double Omega1,Omega2,spin1x,spin1y,spin1z,spin2x,spin2y,spin2z; //smadar adding 
	double hz,hy,hx,ex,ey,ez,etot,htot;
	double ex_u,ey_u,ez_u,qx_u,qy_u,qz_u,hx_u,hy_u,hz_u;//smadar: unit vector componants with respect to (xyz) plane
	double sing1, cosg1,cosOmega1,sinOmega1;
	double spin1e,spin1q,spin1h,spin2e,spin2q,spin2h;
	double tempi,temp;
	double R1,R2;
	double beta,beta2,spintot,spintot2,spinDh,spinDh2,gamma,gamma2;
	double t_b,tloop_f,Deltat;//,t_inter;
	double QoftV1,QoftV2;
	double kL1=0.5*Q1/(1.-Q1);
	double kL2=0.5*Q2/(1.-Q2);
	char name[500],pathname[500];
        int nsnap;

	double m1_c,m2_c,a1_c,a2_c,e1_c,e2_c,i_c,i1_c,i2_c,g1_c,g2_c,m3_c,R1_c,R2_c,spin1e_c,spin1h_c,spin1q_c,spin2e_c,spin2h_c,spin2q_c,htot_c,type1_c,type2_c,type3_c,beta_c,vp_c,spintot_c,t_i_c,Roche1_c,beta2_c,gamma_c,gamma2_c;
	int sur_c,sur2_c;

	double t_quad,t_GR,t_TF1,t_TF2,tF1,tF2,t_rot1,t_rot2;

	setbuf(stdout, NULL);

	void linint(double *xa,double *ya,int n,double x,double *y);
	// initilizing some paramters:
	minIn=1;
	maxIn=Nlin;
	i1=i2=0;
	sur=2;
	sur2=2;
	Nvar=2000;


	type1=-1;
	type2=-1;
	type3=-1;

	if(argc!=3)    {

	  printf(" Error: please specifay : <path> <number of snap>!\n");
	  return 0;
	}

	strcpy(pathname, argv[1]);
	nsnap=atoi(argv[2]);

	if (nsnap<=0){
	  printf(" Error: nsnap should be larger then 1\n");
	  return 0;
	}

	//	spin1=2*M_PI*365.25/Sp1;
	//	spin2=2*M_PI*365.25/Sp2;

	//sprintf(outputname,"Angular.txt");
	//	out=fopen(outputname,"w");

	flgout=0;
	/* takes in input from input.txt */
	sprintf(name,"%striple.in%d",pathname,nsnap);
	fprintf(stderr, " Reading: %s \n",name);
	input = fopen(name, "r");
	sur=2;
	flgout=0;
	if (ferror(input)!=0) {
	    printf("Error from input file \n");
	    return 0;
	}

	for(;;) {	/* read in until :::, end of header */
		scnret = fscanf(input, "%s", dummy);
		if (!(strcmp(dummy, ":::"))) break;
	}

	fscanf(input, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf ",\
	       &m1,&m2,&m3,&R1,&R2,&spin1,&spin2,&beta,&beta2,&gamma,&gamma2,&a1,&a2,&e1,&e2,&g1,&g2,&i,&t_f);


	for(;;) {	/* read in until :::, end of header */
		scnret = fscanf(input, "%s", dummy);
		if (!(strcmp(dummy, ":::"))) break;
	}

	fscanf(input, "%d",&QUADRUPOLE);

	//	fprintf(stderr," %d\n",QUADRUPOLE);
	for(;;) {	/* read in until :::, end of header */
		scnret = fscanf(input, "%s", dummy);
		if (!(strcmp(dummy, ":::"))) break;
	}

	fscanf(input, "%d",&OCTUPOLE);

	//	fprintf(stderr," %d\n",OCTUPOLE);

	for(;;) {	/* read in until :::, end of header */
		scnret = fscanf(input, "%s", dummy);
		if (!(strcmp(dummy, ":::"))) break;
	}

	fscanf(input, "%d",&GR);

	for(;;) {	/* read in until :::, end of header */
		scnret = fscanf(input, "%s", dummy);
		if (!(strcmp(dummy, ":::"))) break;
	}

	fscanf(input, "%d",&TF);
	//	fprintf(stderr," %d\n",GR);


	for(;;) {	/* read in until :::, end of header */
		scnret = fscanf(input, "%s", dummy);
		if (!(strcmp(dummy, ":::"))) break;
	}

	fscanf(input, "%d",&ML1);

	for(;;) {	/* read in until :::, end of header */
		scnret = fscanf(input, "%s", dummy);
		if (!(strcmp(dummy, ":::"))) break;
	}

	fscanf(input, "%d",&MB1);

	for(;;) {	/* read in until :::, end of header */
		scnret = fscanf(input, "%s", dummy);
		if (!(strcmp(dummy, ":::"))) break;
	}

	fscanf(input, "%d",&ML2);

	for(;;) {	/* read in until :::, end of header */
		scnret = fscanf(input, "%s", dummy);
		if (!(strcmp(dummy, ":::"))) break;
	}

	fscanf(input, "%d",&MB2);

	for(;;) {	/* read in until :::, end of header */
		scnret = fscanf(input, "%s", dummy);
		if (!(strcmp(dummy, ":::"))) break;
	}

	fscanf(input, "%d",&SSE1);

	for(;;) {	/* read in until :::, end of header */
		scnret = fscanf(input, "%s", dummy);
		if (!(strcmp(dummy, ":::"))) break;
	}

	fscanf(input, "%d",&SSE2);

	for(;;) {	/* read in until :::, end of header */
		scnret = fscanf(input, "%s", dummy);
		if (!(strcmp(dummy, ":::"))) break;
	}

	fscanf(input, "%d",&SSE3);

	for(;;) {	/* read in until :::, end of header */
		scnret = fscanf(input, "%s", dummy);
		if (!(strcmp(dummy, ":::"))) break;
	}

	fscanf(input, "%lf",&tMS);

	for(;;) {	/* read in until :::, end of header */
		scnret = fscanf(input, "%s", dummy);
		if (!(strcmp(dummy, ":::"))) break;
	}
	
/* read in control parameters */

	fscanf(input, "%lf",&eps);

      	fclose(input);

	if (TF==1 && QUADRUPOLE==1 && OCTUPOLE==1 && GR==1 && ML1==1 && ML2==1){
	  fprintf(stderr, "   =====  oops....you did not choose any physical effect to explor =====\n");
	  fprintf(stderr, "   ===== now terminating.... bye bye \n");
	  exit(30);
	}

	if (e1==0){
	  e1=TINY;
	}

	
	/* memory allocation */
	x = 	vector(1,n_eq);
	x_i = 	vector(1,n_eq);
	dxdt = 	vector(1,n_eq);
	v = 	vector(1,n_orb_param);
	xp = 	vector(1,kmax);
	yp 	=	matrix(1,n_eq,1,kmax);
	vp = matrix(1,n_orb_param,1,kmax); //smadar: define vp
	ve =	vector(1,n_orb_param);//smadar allocating vexternal saving
	// allocation of vecotrs that use SSE 
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

	/*detotdtTF and dhtotdtTF*/
	detotdtTF = vector(1,1);
	dhtotdtTF = vector(1,1);
	/* End of memory allocation */ 
	ideg=i;
	i *=M_PI/180.0;
	beta *=M_PI/180.0;
	beta2 *=M_PI/180.0;
	gamma *=M_PI/180.0;
	gamma2 *=M_PI/180.0;
	g1*=M_PI/180.0;
	g2*=M_PI/180.0;
	R1*=Rsun;
	R2*=Rsun;
	e1max_calc = sqrt(1-5/3.0*cos(i)*cos(i));
	tMS *=1.e6;
	t_f *= 1e6;
	t_i =  0;
	t_b=t_i-Tprint-1;
	//	t_inter=fabs(t_f-t_i)/400000;//4000;//10000;//100.;// /100000;//100000.;
	//	t_inter=fabs(t_f-t_i)/400000;//10000;//100.;// /100000;//100000.;
	//	t_inter=fabs(t_f-t_i)/4000;//10000;//100.;// /100000;//100000.;

	fprintf(stderr,"spin=%le\n",spin1);
 	spin1=2.*M_PI*365.25/spin1;
	spin2=2.*M_PI*365.25/spin2;
	fprintf(stderr,"spin=%le\n",spin1);
	R1f=R1;



	if (!(SSE1)){
	/*******************************************/
	/*******************READ SSe and reset paratmers **************/
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
  /* test the output */
		
	/*fprintf(stderr,"%i %i\n",Ndis,Nlin);
	for (iSSE=1; iSSE<=Ndis; iSSE++) fprintf(stderr,"%f\n",tdis[iSSE]); 

	for (iSSE=1; iSSE<=Nlin; iSSE++) {
	  fprintf(stderr,"%f %f %f %f %f %f %f %f %f %f %le\n",timeSSE[iSSE],typeSSE[iSSE],MoSSE[iSSE],MtSSE[iSSE],logLSSE[iSSE],logRSSE[iSSE],logTSSE[iSSE],McSSE[iSSE],MenvSSE[iSSE],epochSSE[iSSE],spinSSE[iSSE]); 
	  //	exit(4);
	  //	  break;
	}*/
	//	exit(4);
	
  /*
  exit(3);
*/
	//get the index assosiated with the discontinuities 
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
	  /*	  for (iSSE=1;iSSE<=Ndis;iSSE++){
	    type1[iSSE]=typeSSE[indexdis[iSSE]-1];
	    }*/
	  //  exit(2);
	fprintf(stderr,"in SSE: spin=%le\n",spin1);
	//	for  (iSSE=1; iSSE<=Ndis; iSSE++) fprintf(stderr,"%d tdis[%d]=%le Myr type %d\n",indexdis[iSSE],iSSE,tdis[iSSE], type1[iSSE]);
	for  (iSSE=1; iSSE<=Ndis; iSSE++) fprintf(stderr,"%d tdis[%d]=%le Myr \n",indexdis[iSSE],iSSE,tdis[iSSE]);
	//	exit(2);
	// resetting values according to SSE:
	  for (iSSE=1; iSSE<=Nlin; iSSE++) timeSSE[iSSE]=timeSSE[iSSE]*1.e6;
	  for (iSSE=1; iSSE<=Ndis; iSSE++) tdis[iSSE]=tdis[iSSE]*1.e6;

	  toutput=vector(1,100*Ndis+101);
	  for (count=1; count<100; count+=1){
	    toutput[count]=(count*tdis[1]/100.);
	  }
	  for (count=100; (count/100)<Ndis; count+=100){
	    toutput[count]=tdis[count/100];
	    for (count2=1; count2<100; count2+=1){
	      toutput[count+count2]=(tdis[(count/100)]+count2*(tdis[(count/100)+1]-tdis[(count/100)])/100.);
	    }
	  }	  
	  toutput[100*Ndis]=tdis[Ndis];
	  for (count=1; count<=100; count+=1){
	    toutput[100*Ndis+count]=(tdis[Ndis]+count*(t_f+tMS-tdis[Ndis])/100.);
	  }
	  toutput[100*Ndis+101]=t_f+tMS+1e9;

	  mvec=vector(1,(Ndis+1)*Nvar);
	  Rvec=vector(1,(Ndis+1)*Nvar);
	  spinvec=vector(1,(Ndis+1)*Nvar);
	  dmdtvec=vector(1,(Ndis+1)*Nvar);
	  tofdmdtvec=vector(1,(Ndis+1)*Nvar);
	  dRdtvec=vector(1,(Ndis+1)*Nvar);
	  dspindtvec=vector(1,(Ndis+1)*Nvar);
	  tmidDis=vector(1,(Ndis+2));
	  j=1;
	  tofdmdtvec[1]=timeSSE[1];
	  mvec[1]=MtSSE[1];
	  Rvec[1]=logRSSE[1];
	  for (count=1;count<=Nlin;count++){
	    spinSSE[count]=log(spinSSE[count]);
	  }
	  spinvec[1]=spinSSE[1];

	  for (count=1;count<=Ndis+2;count++){
	    if (count == 1){
	      tmidDis[count]=tofdmdtvec[1]+1e3;
	    }
	    else if (count==2){
	      tmidDis[count]=(tdis[count-1]-(tdis[count]-tdis[count-1])/2.);
	    }
	    else if (count>2 && count<Ndis+2){
	      tmidDis[count]=(tdis[count-1]-(tdis[count-1]-tdis[count-2])/2.);
	    }
	    else if (count==Ndis+2){
	      tmidDis[count]=(tdis[count-2]+(tdis[count-2]-tdis[count-3])/2.);
	    }  
	  }

	  //for (count=1;count<=Ndis+2;count++) fprintf(stderr, "%le\n",tmidDis[count]);

	  for (count2=1;count2<=Ndis+1;count2++){  
	    for (count=1;count<=Nvar;count++){
	      if (count!=1 && count2==1){
		tofdmdtvec[count]=tmidDis[count2]+(count)*(tmidDis[count2+1]-tmidDis[count2])/(Nvar+1);
		linintMinMax(timeSSE,MtSSE,1,Nlin,tofdmdtvec[count],&mvec[count]);
		linintMinMax(timeSSE,spinSSE,1,Nlin,tofdmdtvec[count],&spinvec[count]);
		linintMinMax(timeSSE,logRSSE,1,Nlin,tofdmdtvec[count],&Rvec[count]);
	      }
	      else if (count2>1){
		tofdmdtvec[count+(count2-1)*Nvar]=tmidDis[count2]+(count-1)*(tmidDis[count2+1]-tmidDis[count2])/Nvar;
		linintMinMax(timeSSE,MtSSE,1,Nlin,tofdmdtvec[count+(count2-1)*Nvar],&mvec[count+(count2-1)*Nvar]);
		linintMinMax(timeSSE,spinSSE,1,Nlin,tofdmdtvec[count+(count2-1)*Nvar],&spinvec[count+(count2-1)*Nvar]);
		linintMinMax(timeSSE,logRSSE,1,Nlin,tofdmdtvec[count+(count2-1)*Nvar],&Rvec[count+(count2-1)*Nvar]);
	      }
	    }
	  }
	    

	  for (count=2;count<=(Ndis+1)*Nvar;count++){
	    if (fabs(tofdmdtvec[count]-tofdmdtvec[count-1]) != 0.){
	      dmdtvec[j]=( (mvec[count]) - (mvec[count-1]))/(fabs(tofdmdtvec[count])-(tofdmdtvec[count-1]));
	      dspindtvec[j]=(( (spinvec[count]) - (spinvec[count-1]))/(fabs(tofdmdtvec[count])-(tofdmdtvec[count-1])));
	      dRdtvec[j]=( pow(10.,Rvec[count]) - pow(10.,Rvec[count-1]))/(fabs(tofdmdtvec[count]-tofdmdtvec[count-1]));
	    }
	    else {
	      dmdtvec[j]=0.;
	      dspindtvec[j]=0.;
	      dRdtvec[j]=0.;
	      fprintf(stderr, "%d %le \n",j,fabs(tofdmdtvec[j+1]-tofdmdtvec[j]));
	    }
	    j++;
	  }
	  
	  dmdtvec[(Ndis+1)*Nvar]=((MtSSE[Nlin])-(mvec[(Ndis+1)*Nvar]))/((timeSSE[indexdis[Ndis]+1])-(tofdmdtvec[(Ndis+1)*Nvar]));
	  dspindtvec[(Ndis+1)*Nvar]=(((spinSSE[Nlin])-(spinvec[(Ndis+1)*Nvar]))/((timeSSE[indexdis[Ndis]+1])-(tofdmdtvec[(Ndis+1)*Nvar])));
	  dRdtvec[(Ndis+1)*Nvar]=(pow(10.,logRSSE[Nlin])-pow(10.,Rvec[(Ndis+1)*Nvar]))/((timeSSE[indexdis[Ndis]+1])-(tofdmdtvec[(Ndis+1)*Nvar]));
	  j++;
	    

	  if (tMS>0.1){
	    //   linint(timeSSE,MtSSE,Nlin,tdis[5],&m1);
	    linint(timeSSE,MtSSE,Nlin,tMS,&m1);
	    linint(timeSSE,logRSSE,Nlin,tMS,&R1);
	    R1= pow(10.,R1);//*Rsun;
	    //linint(timeSSE,spinSSE,Nlin,tMS,&spin1);
	    //spin1=exp(spin1);
	    fprintf(stderr,"   SSE says that at the age of tMS=%le yrs the mass of the star is actualy m1=%le M_sun\n",tMS,m1);
	    fprintf(stderr,"                                                           and R1=%le Rsun and spin=%le\n",R1,spin1);
	    R1=R1*Rsun;
	  }
	  else {
	    spin1=exp(spinSSE[1]);
	R1=pow(10.,logRSSE[1])*Rsun;
	}
	linintMinMax(timeSSE,logRSSE,indexdis[Ndis-1],Nlin,tdis[Ndis],&R1f);
	R1f= pow(10.,R1f)*Rsun;
	fprintf(stderr,"   final raduis is: %le AU which is %le Rsun \n",R1f,R1f/Rsun);
	//	exit(4);
	/*******************************************/
	}
	/* release some memory */


	//	free_vector(typeSSE,1,Nlinmax);
	free_vector(MoSSE,1,Nlinmax);
	free_vector(logTSSE,1,Nlinmax);
	free_vector(McSSE,1,Nlinmax);
	free_vector(epochSSE,1,Nlinmax);
	free_vector(logLSSE,1,Nlinmax);

	/* End of: release some memory */

	R2f = R2;

	if (!(SSE2)){
	/*******************************************/
	/*******************READ SSe and reset paratmers **************/
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

	readssecl(m2, &Ndis2, &Nlin2, tdis2, timeSSE2, typeSSE2, MoSSE2, MtSSE2, logLSSE2, logRSSE2, logTSSE2, McSSE2, MenvSSE2, epochSSE2, spinSSE2);
  /* test the output */
	/*	
	fprintf(stderr,"%i %i\n",Ndis,Nlin);
	for (iSSE=1; iSSE<=Ndis; iSSE++) fprintf(stderr,"%f\n",tdis[iSSE]); 

	for (iSSE=1; iSSE<=Nlin; iSSE++) {
	  fprintf(stderr,"%f %f %f %f %f %f %f %f %f %f %le\n",timeSSE[iSSE],typeSSE[iSSE],MoSSE[iSSE],MtSSE[iSSE],logLSSE[iSSE],logRSSE[iSSE],logTSSE[iSSE],McSSE[iSSE],MenvSSE[iSSE],epochSSE[iSSE],spinSSE[iSSE]); 
	  //	exit(4);
	  //	  break;
	}
	//	exit(4);
	*/
  /*
  exit(3);
*/
	//get the index assosiated with the discontinuities 
	count=1;
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
	  }
	  /*	  for (iSSE=1;iSSE<=Ndis;iSSE++){
	    type1[iSSE]=typeSSE[indexdis[iSSE]-1];
	    }*/
	  //  exit(2);
	fprintf(stderr,"in SSE: spin=%le\n",spin2);
	//	for  (iSSE=1; iSSE<=Ndis; iSSE++) fprintf(stderr,"%d tdis[%d]=%le Myr type %d\n",indexdis[iSSE],iSSE,tdis[iSSE], type1[iSSE]);
	for  (iSSE=1; iSSE<=Ndis2; iSSE++) fprintf(stderr,"%d tdis[%d]=%le Myr \n",indexdis2[iSSE],iSSE,tdis2[iSSE]);
	//	exit(2);
	// resetting values according to SSE:
	  for (iSSE=1; iSSE<=Nlin2; iSSE++) timeSSE2[iSSE]=timeSSE2[iSSE]*1.e6;
	  for (iSSE=1; iSSE<=Ndis2; iSSE++) tdis2[iSSE]=tdis2[iSSE]*1.e6;

	  mvec2=vector(1,(Ndis2+1)*Nvar);
	  Rvec2=vector(1,(Ndis2+1)*Nvar);
	  spinvec2=vector(1,(Ndis2+1)*Nvar);
	  dmdtvec2=vector(1,(Ndis2+1)*Nvar);
	  tofdmdtvec2=vector(1,(Ndis2+1)*Nvar);
	  dRdtvec2=vector(1,(Ndis2+1)*Nvar);
	  dspindtvec2=vector(1,(Ndis+1)*Nvar);
	  tmidDis2=vector(1,(Ndis2+2));
	  j=1;
	  tofdmdtvec2[1]=timeSSE2[1];
	  mvec2[1]=MtSSE2[1];
	  Rvec2[1]=logRSSE2[1];
	  for (count=1;count<=Nlin2;count++){
	    spinSSE2[count]=log(spinSSE2[count]);
	  }
	  spinvec2[1]=spinSSE2[1];

	  for (count=1;count<=Ndis2+2;count++){
	    if (count == 1){
	      tmidDis2[count]=tofdmdtvec2[1]+1e3;
	    }
	    else if (count==2){
	      tmidDis2[count]=(tdis2[count-1]-(tdis2[count]-tdis2[count-1])/2.);
	    }
	    else if (count>2 && count<Ndis2+2){
	      tmidDis2[count]=(tdis2[count-1]-(tdis2[count-1]-tdis2[count-2])/2.);
	    }
	    else if (count==Ndis2+2){
	      tmidDis2[count]=(tdis2[count-2]+(tdis2[count-2]-tdis2[count-3])/2.);
	    }  
	  }

	  for (count2=1;count2<=Ndis2+1;count2++){  
	    for (count=1;count<=Nvar;count++){
	      if (count!=1 && count2==1){
		tofdmdtvec2[count]=tmidDis2[count2]+(count)*(tmidDis2[count2+1]-tmidDis2[count2])/(Nvar+1);
		linintMinMax(timeSSE2,MtSSE2,1,Nlin2,tofdmdtvec2[count],&mvec2[count]);
		linintMinMax(timeSSE2,spinSSE2,1,Nlin2,tofdmdtvec2[count],&spinvec2[count]);
		linintMinMax(timeSSE2,logRSSE2,1,Nlin2,tofdmdtvec2[count],&Rvec2[count]);
	      }
	      else if (count2>1){
		tofdmdtvec2[count+(count2-1)*Nvar]=tmidDis2[count2]+(count-1)*(tmidDis2[count2+1]-tmidDis2[count2])/Nvar;
		linintMinMax(timeSSE2,MtSSE2,1,Nlin2,tofdmdtvec2[count+(count2-1)*Nvar],&mvec2[count+(count2-1)*Nvar]);
		linintMinMax(timeSSE2,spinSSE2,1,Nlin2,tofdmdtvec2[count+(count2-1)*Nvar],&spinvec2[count+(count2-1)*Nvar]);
		linintMinMax(timeSSE2,logRSSE2,1,Nlin2,tofdmdtvec2[count+(count2-1)*Nvar],&Rvec2[count+(count2-1)*Nvar]);
	      }
	    }
	  }


	  for (count=2;count<=(Ndis2+1)*Nvar;count++){
	    if (fabs(tofdmdtvec2[count]-tofdmdtvec2[count-1]) != 0.){
	      dmdtvec2[j]=( (mvec2[count]) - (mvec2[count-1]))/(fabs(tofdmdtvec2[count])-(tofdmdtvec2[count-1]));
	      dspindtvec2[j]=( (spinvec2[count]) - (spinvec2[count-1]))/(fabs(tofdmdtvec2[count])-(tofdmdtvec2[count-1]));
	      dRdtvec2[j]=( pow(10.,Rvec2[count]) - pow(10.,Rvec2[count-1]))/(fabs(tofdmdtvec2[count]-tofdmdtvec2[count-1]));
	    }
	    else {
	      dmdtvec2[j]=0.;
	      dspindtvec2[j]=0.;
	      dRdtvec2[j]=0.;
	      fprintf(stderr, "%d %le \n",j,fabs(tofdmdtvec2[j+1]-tofdmdtvec2[j]));
	    }
	    j++;
	  }
	  
	  dmdtvec2[(Ndis2+1)*Nvar]=((MtSSE2[Nlin])-(mvec2[(Ndis2+1)*Nvar]))/((timeSSE2[indexdis2[Ndis2]+1])-(tofdmdtvec2[(Ndis2+1)*Nvar]));
	  dspindtvec2[(Ndis2+1)*Nvar]=((spinSSE2[Nlin])-(spinvec2[(Ndis2+1)*Nvar]))/((timeSSE2[indexdis2[Ndis2]+1])-(tofdmdtvec2[(Ndis2+1)*Nvar]));
	  dRdtvec2[(Ndis2+1)*Nvar]=(pow(10.,logRSSE2[Nlin2])-pow(10.,Rvec2[(Ndis2+1)*Nvar]))/((timeSSE2[indexdis2[Ndis2]+1])-(tofdmdtvec2[(Ndis2+1)*Nvar]));
	  j++;


	  if (tMS>0.1){
	    //   linint(timeSSE,MtSSE,Nlin,tdis[5],&m1);
	    linint(timeSSE2,MtSSE2,Nlin2,tMS,&m2);
	    linint(timeSSE2,logRSSE2,Nlin2,tMS,&R2);
	    R2= pow(10.,R2);//*Rsun;
	    //linint(timeSSE2,spinSSE2,Nlin2,tMS,&spin2);
	    //spin2=exp(spin2);
	    fprintf(stderr,"   SSE says that at the age of tMS=%le yrs the mass of the *second* star is actualy m2=%le M_sun\n",tMS,m2);
	    fprintf(stderr,"                                                           and R2=%le Rsun and spin=%le\n",R2,spin2);
	    R2=R2*Rsun;
	  }
	  else {
	    spin2=exp(spinSSE2[1]);
	R2=pow(10.,logRSSE2[1])*Rsun;
	}
	linintMinMax(timeSSE2,logRSSE2,indexdis2[Ndis2-1],Nlin2,tdis2[Ndis2],&R2f);
	R2f= pow(10.,R2f)*Rsun;
	fprintf(stderr,"   final raduis is: %le AU which is %le Rsun \n",R2f,R2f/Rsun);
	//	exit(4);
	/*******************************************/
	}
	/* release some memory */


	//	free_vector(typeSSE,1,Nlinmax);
	free_vector(MoSSE2,1,Nlinmax);
	free_vector(logTSSE2,1,Nlinmax);
	free_vector(McSSE2,1,Nlinmax);
	free_vector(epochSSE2,1,Nlinmax);
	free_vector(logLSSE2,1,Nlinmax);

	/* End of: release some memory */



	if (!(SSE3)){
	/*******************************************/
	/*******************READ SSe and reset paratmers **************/
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
	  //  exit(2);
	fprintf(stderr,"in SSE: spin=%le\n",spin1);
	for  (iSSE=1; iSSE<=Ndis3; iSSE++) fprintf(stderr,"%d tdis[%d]=%le Myr\n",indexdis3[iSSE],iSSE,tdis3[iSSE]); 
	// resetting values according to SSE:
	  for (iSSE=1; iSSE<=Nlin3; iSSE++) timeSSE3[iSSE]=timeSSE3[iSSE]*1.e6;
	    for (iSSE=1; iSSE<=Ndis3; iSSE++) tdis3[iSSE]=tdis3[iSSE]*1.e6;



	if (tMS>0.1){
	  //   linint(timeSSE,MtSSE,Nlin,tdis[5],&m1);
	    linint(timeSSE3,MtSSE3,Nlin3,tMS,&m3);
	
	  fprintf(stderr,"   SSE says that at the age of tMS=%le yrs the mass of the *third* star is actualy m3=%le M_sun\n",tMS,m3);

	 
	}
	free_vector(logRSSE3,1,Nlinmax);
	free_vector(spinSSE3,1,Nlinmax);
	//	free_vector(typeSSE3,1,Nlinmax);
	free_vector(MoSSE3,1,Nlinmax);
	free_vector(MenvSSE3,1,Nlinmax);
	free_vector(logTSSE3,1,Nlinmax);
	free_vector(McSSE3,1,Nlinmax);
	free_vector(epochSSE3,1,Nlinmax);
	free_vector(logLSSE3,1,Nlinmax);
	//	exit(4);
	/*******************************************/
	}

	fprintf(stderr,"spin=%le\n",spin1);
	//		exit(2);
	/* various timescales */
	P1 = 2*M_PI*sqrt(a1*a1*a1/(k2*(m1+m2)));
	P2 = 2*M_PI*sqrt(a2*a2*a2/(k2*(m1+m2+m3)));
	P_koz = (8./(15*M_PI))*((m1+m2+m3)/m3)*(P2*P2/P1)*sqrt(1-e2*e2)*(1-e2*e2);
	P_GR  = 3.36e7*(1-e1*e1)*P1*a1/m2;
	range = 5/3.0*cos(i)*cos(i);	//to find minimium P_GR 
	P_GR*=range;


	/* some usful parameters */
	m12 = m1 + m2;
	m123 = m1 + m2 + m3;

	/* i1 and i2 */
	L1 = m1 * m2  / m12  * sqrt( k2 * m12  * a1 );
	L2 = m3 * m12 / m123 * sqrt( k2 * m123 * a2 );
	G1 = L1 * sqrt(1 - e1*e1);
	G2 = L2 * sqrt(1 - e2*e2);
        Hsq = G1*G1 + 2.*G1*G2*cos(i) + G2*G2;

	fprintf(stderr,"Hsq=%le L1=%le L2=%le L1art=%le\n",Hsq,L1,L2,m1 * m2  / m12  * sqrt( k2 * m12 ));

	theta1=(Hsq-G2*G2+G1*G1)/(2.*sqrt(Hsq)*G1);
	theta2=(Hsq+G2*G2-G1*G1)/(2.*sqrt(Hsq)*G2);

	if(fabs(theta1)<=1)    i1=acos(theta1);
	else if(theta1>1)      i1 = 0.0;
	else if(theta1<-1)      i1 = M_PI;

	if(fabs(theta2)<=1)    i2=acos(theta2);
	else if(theta2>1)      i2 = 0.0;
	else if(theta2<-1)      i2 = M_PI;


	fprintf(stderr,"H1=%le H2=%le\n",G1*theta1,G2*theta2);
	fprintf(stderr,"G1=%le G2=%le\n",G1,G2);
fprintf(stderr,"theta1=%le theta2=%le\n",theta1,theta2);
//exit(1);
	////

	///c
	htot=sqrt(k2 * (m1+m2) *a1 *(1-SQR(e1)));

	spin1h=spin1*cos(beta);
	// damping the rest in the two other directions.
	spin1e=spin1*sin(beta)*cos(gamma);//spin1q =sqrt(( SQR(spin1)-SQR(spin1h))/2 );
	spin1q=spin1*sin(beta)*sin(gamma);

	//Not assuming that the planet is allined
	spin2h=spin2*cos(beta2);
	spin2e=spin2*sin(beta2)*cos(gamma2);//spin2q =sqrt(( SQR(spin2)-SQR(spin2h))/2 );
	spin2q=spin2*sin(beta2)*sin(gamma2);

	spin1x=2*M_PI*365.25/spin1h;
	spin1y=2*M_PI*365.25/spin1e;
	spin1z=2*M_PI*365.25/spin1q;
	spin2x=2*M_PI*365.25/spin2h;
	spin2y=2*M_PI*365.25/spin2e;
	spin2z=2*M_PI*365.25/spin2q;

	/* control parameters */

	h1  	= eps * P2;
	hmin = 0;
	dxsav 	= 10  * P2;

	/* initial condition */
	x_i[1] = e1*sin(g1);
	x_i[2] = e1*cos(g1);
	x_i[3] = e2*sin(g2);
	x_i[4] = e2*cos(g2);

	//smadar: adding parameters for TF


	x_i[5]=cos(i1);//theta1;
	x_i[6]=cos(i2);//theta2;
	x_i[7]=cos(i);//theta;
	x_i[8]=htot; 
	x_i[9]=spin1e;
	x_i[10]=spin1q;
	x_i[11]=spin1h;
	x_i[12]=spin2e;
	x_i[13]=spin2q;
	x_i[14]=spin2h;
	x_i[15]=R1;
	x_i[16]=m1;
	x_i[17]=G1;
	x_i[18]=G2;
	x_i[19]=R2;
	x_i[20]=m2;

	//x_i[15] = 0;//Omega1 - no need - only need it's derivative and calcualting it in octupole.c
	
	// parameters:
	v[1]	= i;
	v[2] 	= m1;
	v[3] 	= m2;
	v[4] 	= m3;
	v[5]	= a1;
	v[6]	= a2;
	v[7]    = R1;
	v[8]    = R2;

	/* print out the initial conditions */
	fprintf(stderr,"\n============================== Initial Conditions ==============================\n");
	fprintf(stderr,"	m1=%1.3f m2=%1.3f m3=%1.3f [Msun] a1=%1.3f a2=%1.3f [AU] k2=%1.3f\n"\
				,m1,m2,m3,a1,a2,k2);
	
	fprintf(stderr,"	R1=%le R2=%le  [AU] Spin period: %le %le [rad/yr] \n"\
		,R1,R2,spin1,spin2);
       
	fprintf(stderr,"	e1=%1.3f g1=%1.5f e2=%1.3f g2=%1.5f i=%1.3f\n"\
				,e1,g1/M_PI*180,e2,g2/M_PI*180,i/M_PI*180);

	fprintf(stderr,"        Initial incl: :: i1=%f i2=%f i1+i2=%f i=%f\n\n",i1*180/M_PI,i2*180/M_PI,(i1+i2)*180/M_PI,i*180/M_PI);
	fprintf(stderr,"	Integration stopped at %.0f years, Inner Period=%1.0g years, Outer Period=%e years\n"\
				,t_f,P1, P2);
	if(i>i_koz)  	{
	  fprintf(stderr,"	e1_max calc'd from i= %lf\n",e1max_calc);
	  fprintf(stderr,"	This means that r_p= %le AU\n",a1*(1-e1max_calc));
	  //	exit(1);
	}
	else if(i<i_koz)	fprintf(stderr, "	Inclination too low\n");
	fprintf(stderr, "	Kozai Period= %1.5e years\n	GR Period= %.5e years---%.5e years\n\n"\
				,P_koz,P_GR*range, P_GR);

	fprintf(stderr, "       Running OSPE with the following options on: \n");
	if(!(OCTUPOLE)){
	fprintf(stderr, "       OCTUPOLE \n");
	}
	if(!(QUADRUPOLE)){
	fprintf(stderr, "       QUADRUPOLE \n");
	}
	if(!(GR)){
	fprintf(stderr, "       GR \n");
	}
	if(!(TF)){
	fprintf(stderr, "       TIDAL FRICTION \n");
	}
	if(!(ML1)){
	fprintf(stderr, "       MASS LOSS for m1 \n");
	}
	if(!(MB1)){
	fprintf(stderr, "       Magnetic Braking for m1 \n");
	}
	if(!(ML2)){
	fprintf(stderr, "       MASS LOSS for m2 \n");
	}
	if(!(MB2)){
	fprintf(stderr, "       Magnetic Braking for m2 \n");
	}
	if(!(SSE1)){
	fprintf(stderr, "       SSE for m1\n");
	}
	if(!(SSE2)){
	fprintf(stderr, "       SSE for m2\n");
	}
	if(!(SSE3)){
	fprintf(stderr, "       SSE for m3\n");
	}

	fprintf(stderr, "   \n");
	fprintf(stderr, "        TV1=%lfyr TV2=%lfyr\n",TV1,TV2);
	fprintf(stderr, "        kL1=%lfyr kL2=%lfyr\n",kL1,kL2);
	QoftV1=4.*kL1*4.*M_PI*M_PI*m1*(1./365.)*TV1/(3.*SQR(1.+2.*kL1)*CBE(R1)*2.*M_PI);
	QoftV2=4.*kL2*4.*M_PI*M_PI*m2*(3./365.)*TV2/(3.*SQR(1.+2.*kL2)*CBE(R2)*2.*M_PI);
	fprintf(stderr, "        QoftV1=%le QoftV2=%le\n",QoftV1,QoftV2);
	fprintf(stderr, "   \n epsilon=%le epsilonM=%le\n",a1/a2*e2/(1-e2*e2),a1/a2*e2/(1-e2*e2)*fabs(m1-m2)/(m1+m2));
	//double Roche=2.*R2*pow(m1 / m2,1./3.);
   //  fprintf(stderr, " \n The 2*Roche limit is %le AU\n",Roche);
	//	double Roche;//=3.*R2*pow(m2/(m2+m1),-1./3.)/2.; 
	Roche2=Roche(R2,m2/(m1+m2));//R2*pow(m2/(m2+m1),-1./3.)/0.6;
	//	fprintf(stderr, " \n The 2*Roche limit is %le AU also note that 0.462*q^1/3*a =%le\n",Roche,0.462*pow(m2/m1,1./3.)*R1);
	fprintf(stderr, " \n The 2*Roche limit is %le AU\n",Roche2);
	//	exit(1);
	// fprintf(stderr, " \n The 2*Roche limit is %le AU\n",3.*R2*pow(m2/(m2+m1),-1./3.)/2.);
	/*
	double sq=m2/m1;
	double val=pow(sq,1./3.);
	fprintf(stderr, " \n The 2*Roche limit is %le AU\n",0.49*val*val/(0.6*val*val+log(1.0+val))*a1);
	exit(1);
	*/
	fprintf(stderr,"\n================================================================================\n\n");
	fprintf(stderr, "    Some control paramters:\n");
	//	h1 = eps;
	fprintf(stderr, "    first step: %le\n",h1);
 	fprintf(stderr, "    eps  %le\n",eps);
	fprintf(stderr, "Integrating ....\n\n");
	//	exit(1);
	/******************************************************************/
	
	//print the initial values
	//	printf("%lf\t",xp[j]);
	//printf("%1.6lf \t %1.6lf \t %1.6lf \t %1.6lf  %1.6lf \t  %1.6lf \t  %1.6lf   \t  %1.6lf \t  %1.6lf \t  %1.6lf \t  %1.6lf \t  %1.6lf\n",e1,e2,g1,g2,a1,i1*180./M_PI,i2*180./M_PI,i*180./M_PI,spin1h,spintot,beta*180./M_PI,vp[1][j]*180/M_PI);


	
	//t_i=TSTART;
	

  fprintf(stderr,"before integrations\n");
	  fprintf(stderr," Hsquar=%le t=%le yr\n", Hsq,xp[j]);
	  fprintf(stderr,"G1=%le G2=%le cosi=%le htot=%le\n",G1,G2,cos(i),htot);
	  fprintf(stderr,"m1=%le m2=%le m3=%le e2=%le a2=%le a1=%le e1=%le\n",m1,m2,m3,e2,a2,a1,e1);
	  /******************/
	  /* first find where tMS is among tdis */
	  /* 
	  min=1;
	  max=Ndis;
	  for (count=1; count<=Ndis; count++){
    
	    if (tMS>=tdis[Ndis]){

	      min=1;
	      max=1;

	      break;
	    }
	    if (tMS>tdis[count]){
	      min=count;
	      count++;
	    }
	    else {
	      max=Ndis;
	      break;
	    }
	    }*/
	  /******************/

	  //	  fprintf(stderr, "min=%d max=%d\n",min,max);
	/* integrate the octupole equations (use bsstep rather than rkqs) */

	  if (!(SSE1)){
	    if (xp[j]+tMS>tdis[Ndis]){
	      type1=typeSSE[indexdis[Ndis]+1];
	            
	    }
	    else {
	      if (xp[j]+tMS<timeSSE[2]) type1=typeSSE[1];
	      else {
		//fprintf(stderr, " Are you here type 1\n");
		linint(timeSSE,typeSSE,Nlin,xp[j]+tMS,&type1);
	      }
	      //    fprintf(stderr, " :: type1=%le\n",type1);
	    }
	  }

	  if (!(SSE2)){
	    if (xp[j]+tMS>tdis2[Ndis2]){
	      type2=typeSSE2[indexdis2[Ndis2]+1];
	            
	    }
	    else {
	      if (xp[j]+tMS<timeSSE2[2]) type2=typeSSE2[1];
	      else {
		//fprintf(stderr, " Are you here type 1\n");
		linint(timeSSE2,typeSSE2,Nlin2,xp[j]+tMS,&type2);
	      }
	      //    fprintf(stderr, " :: type1=%le\n",type1);
	    }
	  }
	  //  fprintf(stderr, " :: tdis[Ndis]=%le xp[j]=%le type1=%e\n",tdis[Ndis],xp[j],type1);
	  //  exit(5);
	  if (!(SSE3)){
	    if (xp[j]+tMS>tdis3[Ndis3]){
	      type3=typeSSE3[indexdis3[Ndis3]+1];
	    }
	    else {
	      if (xp[j]+tMS<timeSSE3[2]) type3=typeSSE3[1];
	      else {
		//fprintf(stderr, " Are you here type 3\n");
		linint(timeSSE3,typeSSE3,Nlin3,xp[j]+tMS,&type3);
	      }
	      //    findt3(xp[j]+tMS,&minIn,&maxIn);
	      //    linintMinMax(timeSSE,typeSSE,minIn,maxIn,xp[j]+tMS,&type3);
	    }
	        
	  }
	  //  fprintf(stderr, " :: tdis3[Ndis]=%le xp[j]=%le type3=%e\n",tdis3[Ndis3],xp[j],type3);
	  //  exit(4);
	  //  fprintf(stderr, " :: tdis[Ndis]=%le xp[j]=%le type1=%e\n",tdis[Ndis],xp[j],type1);
	  // if (t-t_b >=Tprint){
	  printf("%d \t %d \t %lf \t",1,1,xp[j]);
	  printf("%le \t %le \t %le \t %le \t %le \t  %le \t  %le   \t  %le \t  %le \t  %le \t  %le \t  %le \t  %le   \t  %le \t  %le \t  %le \t  %le \t  %le \t  %le \t  %le  \t  %le \t %le \t %le \t %le \t %le \t %le \t %lf \t %lf \t %lf \t %le \t %le \t %le \t \n",e1,e2,g1,g2,a1,i1*180./M_PI,i2*180./M_PI,i*180./M_PI,spin1h,spintot,beta*180./M_PI,vp[1][j]*180/M_PI,spin1e,spin1q,spin2e,spin2q,spin2h,htot,m1,R1,m2,R2,a2,m3,Roche1,R1/(0.462 * pow( m1/m2,1./3.  )),type1,type2,type3,beta2*180./M_PI,gamma*180./M_PI,gamma2*180./M_PI);

	  t_ini = t_i;

	  P1 = 2*M_PI*sqrt(a1*a1*a1/(k2*(m1+m2)));
	  P2 = 2*M_PI*sqrt(a2*a2*a2/(k2*(m1+m2+m3)));
	  P_koz = (8./(15*M_PI))*((m1+m2+m3)/m3)*(P2*P2/P1)*sqrt(1-e2*e2)*(1-e2*e2);
	  t_quad = 2.*M_PI*pow(a2,3)*pow(1.0-pow(e2,2),1.5)*sqrt(m1+m2)/(pow(a1,1.5)*m3*k2);
	  t_GR = 2*M_PI*pow(a1,2.5)*c2*(1-pow(e1,2))/(3*pow(k2*sqrt(m1+m2),3));
	  t_TF1 = pow(a1,13./2.)*m2*pow(1-(e1*e1),5)/(sqrt(k2)*kL2*(1.+1.5*(e1*e1)+pow(e1,4)/8.)*m1*sqrt(m1+m2)*pow(R2,5));
	  t_TF2 = pow(a1,13./2.)*m1*pow(1-(e1*e1),5)/(sqrt(k2)*kL1*(1.+1.5*(e1*e1)+pow(e1,4)/8.)*m2*sqrt(m1+m2)*pow(R2,5));
	  tF1 = TV1*pow(a1/R1,8)*SQR(m1)/ ( 9.*(m1+m2)*m2*SQR(1+2.*kL1));
	  tF2 = TV2*pow(a1/R2,8)*SQR(m2)/ ( 9.*(m1+m2)*m1*SQR(1+2.*kL2));
	  t_rot1 = sqrt(k2)*pow(a1,7./2.)*pow(1.-e1*e1,2.)/(kL1*pow(spintot,2)*sqrt(m1+m2)*pow(R1,5));
	  t_rot2 = sqrt(k2)*pow(a1,7./2.)*pow(1.-e1*e1,2.)/(kL2*pow(spintot2,2)*sqrt(m1+m2)*pow(R2,5));

	  if (P_koz/100.<=fabs(t_f-t_ini)/400000){
	    t_inter=P_koz/100.;
	  }
	  else t_inter=fabs(t_f-t_ini)/400000;

	  /*if (t_TF1 < t_inter || t_TF2 < t_inter || tF1 < t_inter || tF2 < t_inter || t_rot1 < t_inter || t_rot2 < t_inter){
	    if (t_rot1 < t_TF1 && t_rot1 < t_TF2 && t_rot1 < t_rot2 && t_rot1 < tF1 && t_rot1 < tF2) t_inter = fabs(t_rot1);
	    else if (t_rot2 < t_TF1 && t_rot2 < t_TF2 && t_rot2 < t_rot1 && t_rot2 < tF1 && t_rot2 < tF2) t_inter = fabs(t_rot2);
	    else if (t_TF1 < t_rot1 && t_TF1 < t_TF2 && t_TF1 < t_rot2 && t_TF1 < tF1 && t_TF1 < tF2) t_inter = fabs(t_TF1);
	    else if (t_TF2 < t_rot1 && t_TF2 < t_TF1 && t_TF2 < t_rot2 && t_TF2 < tF1 && t_TF2 < tF2) t_inter = fabs(t_TF2);
	    else if (tF1 < t_rot1 && tF1 < t_TF2 && tF1 < t_rot2 && tF1 < t_TF1 && tF1 < tF2) t_inter = fabs(tF1);
	    else if (tF2 < t_rot1 && tF2 < t_TF2 && tF2 < t_rot2 && tF2 < t_TF1 && tF2 < tF1) t_inter = fabs(tF2);

	    if (t_inter < 100.) t_inter = 100.;
	  }
	  else{
	    if (P_koz/250.<=fabs(t_f-t_ini)/400000){
	      t_inter=P_koz/250.;
	    }
	    else t_inter=fabs(t_f-t_ini)/400000;
	    }*/

	  // t_inter=fabs(t_f-t_i)/40000;//10000;//100.;// /100000;//100000.;
	  //if (P_koz/250. <= fabs(t_f-t_i)/400000){
	  //  t_inter=P_koz/250.;
	  //}
	  //else t_inter=fabs(t_f-t_i)/400000;//400000;//10000;//100.;// /100000;//100000.;
	  //	  t_inter=fabs(t_f-t_i)/4000000;//10000;//100.;// /100000;//100000.;
	  count = 1;
	  count2 = 0;
	  while (tMS>=toutput[count+1]) count=count+1;

	  sur_c = 1;
	  sur2_c = 1;
	  t_i_c = 0.;
	  a1_c = a1;
	  a2_c = a2;
	  e1_c = e1;
	  e2_c = e2;
	  m1_c = m1;
	  m2_c = m2;
	  m3_c = m3;
	  g1_c = g1;
	  g2_c = g2;
	  i1_c = i1;
	  i2_c = i2;
	  i_c = i;
	  beta_c = beta;
	  spin1h_c = spin1h;
	  spintot_c = spintot;
	  vp_c = 0.;
	  spin1e_c = spin1e;
	  spin1q_c = spin1q;
	  spin2e_c = spin2e;
	  spin2h_c = spin2h;
	  R1_c = R1;
	  R2_c = R2;
	  Roche1_c = Roche1;
	  htot_c = htot;
	  type1_c = type1;
	  type2_c = type2;
	  type3_c = type3;
	  beta2_c = beta2;
	  gamma_c = gamma;
	  gamma2_c = gamma2;

	  time_t start_time, check_time, duration;
	  start_time = time(0);
	  
	  /* initialize time flags */
	  flagtime1 = 0;
	  flagtime2 = 0;
	  flagtime3 = 0;


	while (t_i<=t_f){// && sur==2){
	
	  /*if (e1>= 0.95){
	    if (P_koz/500.<= 10.) t_inter=10.;
	    else t_inter=P_koz/500.;
	    }*/
	 
	  tloop_f=t_i+t_inter;
	  //fprintf(stderr, " h1=%le eps=%le\n",h1,eps);
	  //fprintf(stderr, " Integrating on: tstart=%le tend=%le\n",t_i,tloop_f);
	  h1=fabs(tloop_f-t_i)/1000000.;//100000.;///10000.;//*eps;///10;
	  //	  h1=P2/100.;
	  clock_t start, finish;
	  start = clock();
	  odeint(x_i, n_eq, t_i, tloop_f, eps, h1, hmin, &nok, &nbad, octupole, integrator);
	  finish = clock();
	  //fprintf(stderr, "Clock time for calculation step = %le\n",((double)(finish - start))/CLOCKS_PER_SEC );
	k_f = kount;
	//	fprintf(stderr, "====\n");
	for(j=1;j<=k_f;j++) {
	  
	  e1 = sqrt(yp[1][j]*yp[1][j] + yp[2][j]*yp[2][j]);
	  e2 = sqrt(yp[3][j]*yp[3][j] + yp[4][j]*yp[4][j]);
	  g1 = atan2(yp[1][j],yp[2][j]);
	  g2 = atan2(yp[3][j],yp[4][j]);
	  t=xp[j];

	  htot=yp[8][j];

	  if (xp[j]+tMS > 0.1 && SSE1 == 0){
	    linint(timeSSE,logRSSE,Nlin,tMS+xp[j],&R1);
	    R1= pow(10.,R1);
	    R1*=Rsun;
	  }
	  else R1=yp[15][j];
	  m1=yp[16][j];

	  //	  m3=vp[4][j];
	  m3=ve[4];

	  G1=yp[17][j];
	  G2=yp[18][j];
	  
	  if (xp[j]+tMS > 0.1 && SSE2 == 0){
	    linint(timeSSE2,logRSSE2,Nlin2,tMS+xp[j],&R2);
	    R2= pow(10.,R2);
	    R2*=Rsun;
	  }
	  else R2=yp[19][j];
	  m2=yp[20][j];

	  m12 = m1 + m2;
	  m123 = m1 + m2 + m3;

	  a1= SQR(htot)/ (k2 * (m1+m2) *(1-SQR(e1)));

	  //	  G1=htot*(m1*m2)/(m1+m2);
	  /*
	  if (abs(G1-htot*(m1*m2)/(m1+m2))>1.e-6){
	    fprintf(stderr, "====\n");
	    fprintf(stderr, " Oh no G1=%le while from htot caluclaiton you get: %le\n",G1,htot*(m1*m2)/(m1+m2));
	    fprintf(stderr, "   you may also want to know m1=%le m2=%le htot=%le\n",m1,m2,htot);
	    fprintf(stderr, "      Terminating now....\n");
	    exit(88);
	    }*/
	  L1= G1/sqrt(1-SQR(e1));
	  L2= G2/sqrt(1-SQR(e2));

	  a2 = SQR(G2 / ( (m1+m2)*m3 ) ) *(m1+m2+m3) / (k2 *(1-SQR(e2)));
	  /*
	  if (t+tMS> 309.4e6){
    fprintf(stderr," m3=%le t=%le tMS=%le a2=%le\n",m3,xp[j],tMS,a2);
    exit(7);
    }*/
	  Roche1=Roche(R1,m1/(m1+m2));
	  Roche2=Roche(R2,m2/(m1+m2));

	  if(fabs(yp[5][j])<=1)    i1=acos(yp[5][j]);
	  else if(yp[5][j]>1)      i1 = 0.0;
	  else if(yp[5][j]<-1)      i1 = M_PI;


	  if(fabs(yp[6][j])<=1)    i2=acos(yp[6][j]);
	  else if(yp[6][j]>1)      i2 = 0.0;
	  else if(yp[6][j]<-1)      i2 = M_PI;

	  if(fabs(yp[7][j])<=1)    i=acos(yp[7][j]);
	  else if(yp[7][j]>1)      i = 0.0;
	  else if(yp[7][j]<-1)      i = M_PI;

	  imu=i1+i2;


	  //L1 = m1 * m2  / m12  * sqrt( k2 * m12  * a1 );
	  // L2 = m3 * m12 / m123 * sqrt( k2 * m123 * a2 );
	  //	  G1 = L1 * sqrt(1 - e1*e1);
	  //	  G2 = L2 * sqrt(1 - e2*e2);


	 Hsq = G1*G1 + 2.*G1*G2*cos(i) + G2*G2;
	  H1 = G1*yp[5][j];
	  H2 = G2*yp[6][j];
	  /*
	  fprintf(stderr,"After integraitons\n");
	  fprintf(stderr," Hsquar=%le t=%le yr\n", Hsq,xp[j]);
	  fprintf(stderr,"G1=%le G2=%le cosi=%le htot=%le\n",G1,G2,cos(i),htot);
	  fprintf(stderr,"m1=%le m2=%le m3=%le e2=%le a2=%le a1=%le e1=%le\n",m1,m2,m3,e2,a2,a1,e1);
	  // exit(1);
	  */
	  // calc the spins:
	  spin1h=yp[11][j];
	  spin2h=yp[14][j];

	  //spin1h = spin*hat(h)*cos(beta)

	  spintot=sqrt(SQR(spin1h)+SQR(yp[10][j])+SQR(yp[9][j])  );
	  spintot2=sqrt(SQR(spin2h)+SQR(yp[13][j])+SQR(yp[12][j])  );
	  spinDh =  spin1h / (spintot);//spin1h / (spintot * htot);
	  spinDh2 = spin2h /spintot2;

	  if(fabs(spinDh)<=1)    beta=acos(spinDh);
	  else if(spinDh>1)      beta = 0.0;
	  else if(spinDh<-1)      beta = M_PI;

	  if(fabs(spinDh2)<=1)    beta2=acos(spinDh2);
	  else if(spinDh2>1)      beta2 = 0.0;
	  else if(spinDh2<-1)      beta2 = M_PI;

	  spin1e=yp[9][j];
	  spin1q=yp[10][j];
	  spin2e=yp[12][j];
	  spin2q=yp[13][j];

	  if(fabs(spin1e/(spintot*sin(beta)))<=1)    gamma=acos(spin1e/(spintot*sin(beta)));
	  else if(spin1e/(spintot*sin(beta))>1)      gamma = 0.0;
	  else if(spin1e/(spintot*sin(beta))<-1)      gamma = M_PI;

	  if(fabs(spin2e/(spintot2*sin(beta2)))<=1)    gamma2=acos(spin2e/(spintot2*sin(beta2)));
	  else if(spin2e/(spintot2*sin(beta2))>1)      gamma2 = 0.0;
	  else if(spin2e/(spintot2*sin(beta2))<-1)      gamma2 = M_PI;

	  
	  if (!(SSE1)){
	    if (xp[j]+tMS>tdis[Ndis]){
	      type1=typeSSE[indexdis[Ndis]+1];
	      
	    }
	    else {
	      if (xp[j]+tMS<timeSSE[2]) type1=typeSSE[1];
	      else {
		//		fprintf(stderr, " Are you here type 1\n");
		linint(timeSSE,typeSSE,Nlin,xp[j]+tMS,&type1);
	      }
	      //	    fprintf(stderr, " :: type1=%le\n",type1);
	    }
	  }

	  if (!(SSE2)){
	    if (xp[j]+tMS>tdis2[Ndis2]){
	      type2=typeSSE2[indexdis2[Ndis2]+1];
	      
	    }
	    else {
	      if (xp[j]+tMS<timeSSE2[2]) type2=typeSSE2[1];
	      else {
		//		fprintf(stderr, " Are you here type 1\n");
		linint(timeSSE2,typeSSE2,Nlin2,xp[j]+tMS,&type2);
	      }
	      //	    fprintf(stderr, " :: type1=%le\n",type1);
	    }
	  }
	  //	  fprintf(stderr, " :: tdis[Ndis]=%le xp[j]=%le type1=%e\n",tdis[Ndis],xp[j],type1);
	  //	  exit(5);
	  if (!(SSE3)){
	    if (xp[j]+tMS>tdis3[Ndis3]){
	      type3=typeSSE3[indexdis3[Ndis3]+1];
	    }
	    else {
	      if (xp[j]+tMS<timeSSE3[2]) type3=typeSSE3[1];
	      else {
		//		fprintf(stderr, " Are you here type 3\n");
		linint(timeSSE3,typeSSE3,Nlin3,xp[j]+tMS,&type3);
	      }
	      //	    findt3(xp[j]+tMS,&minIn,&maxIn);
	      //	    linintMinMax(timeSSE,typeSSE,minIn,maxIn,xp[j]+tMS,&type3);
	    }
	    
	  }
	  //	  fprintf(stderr, " :: tdis3[Ndis]=%le xp[j]=%le type3=%e\n",tdis3[Ndis3],xp[j],type3);
	  //	  exit(4);
	  //	  fprintf(stderr, " :: tdis[Ndis]=%le xp[j]=%le type1=%e\n",tdis[Ndis],xp[j],type1);
	  
	  /*if ( a1*(1.-e1) < a1_c*(1.-e1_c) ){
	    sur_c = sur;
	    sur2_c = sur2;
	    t_i_c = t_i;
	    a1_c = a1;
	    a2_c = a2;
	    e1_c = e1;
	    e2_c = e2;
	    m1_c = m1;
	    m2_c = m2;
	    m3_c = m3;
	    g1_c = g1;
	    g2_c = g2;
	    i1_c = i1;
	    i2_c = i2;
	    i_c = i;
	    beta_c = beta;
	    spin1h_c = spin1h;
	    spintot_c = spintot;
	    vp_c = vp[1][j];
	    spin1e_c = spin1e;
	    spin1q_c = spin1q;
	    spin2e_c = spin2e;
	    spin2h_c = spin2h;
	    R1_c = R1;
	    R2_c = R2;
	    Roche1_c = Roche1;
	    htot_c = htot;
	    type1_c = type1;
	    type2_c = type2;
	    type3_c = type3;
	    beta2_c = beta2;
	    gamma_c = gamma;
	    gamma2_c = gamma2;
	    }*/

	  // if (t-t_b >=Tprint){
	  /*if ( xp[j]+tMS >= toutput[count] ){
	    printf("%d \t %d \t %lf \t",sur_c,sur2_c,t_i_c);
	      printf("%le \t %le \t %le \t %le \t %le \t  %le \t  %le   \t  %le \t  %le \t  %le \t  %le \t  %le \t  %le   \t  %le \t  %le \t  %le \t  %le \t  %le \t  %le \t  %le  \t  %le \t %le \t %le \t %le \t %le \t %le \t %lf \t %lf \t %lf \t %le \t  %le \t  %le \t \n",e1_c,e2_c,g1_c,g2_c,a1_c,i1_c*180./M_PI,i2_c*180./M_PI,i_c*180./M_PI,spin1h_c,spintot_c,beta_c*180./M_PI,vp_c*180/M_PI,spin1e_c,spin1q_c,spin2e_c,spin2q_c,spin2h_c,htot_c,m1_c,R1_c,m2_c,R2_c,a2_c,m3_c,Roche1_c,R1_c/(pow( (m1_c+m2_c)/m2_c,1./3.  )),type1_c,type2_c,type3_c,beta2_c*180./M_PI,gamma_c*180./M_PI,gamma2_c*180./M_PI);

	  printf("%d \t %d \t %lf \t",sur,sur2,t_i);
	  printf("%le \t %le \t %le \t %le \t %le \t  %le \t  %le   \t  %le \t  %le \t  %le \t  %le \t  %le \t  %le   \t  %le \t  %le \t  %le \t  %le \t  %le \t  %le \t  %le  \t  %le \t %le \t %le \t %le \t %le \t %le \t %lf \t %lf \t %lf \t %le \t  %le \t  %le \t \n",e1,e2,g1,g2,a1,i1*180./M_PI,i2*180./M_PI,i*180./M_PI,spin1h,spintot,beta*180./M_PI,vp[1][j]*180/M_PI,spin1e,spin1q,spin2e,spin2q,spin2h,htot,m1,R1,m2,R2,a2,m3,Roche1,R1/(pow( (m1+m2)/m2,1./3.  )),type1,type2,type3,beta2*180./M_PI,gamma*180./M_PI,gamma2*180./M_PI);

	    count++;
	  }*/
	    /*sur_c = sur;
	    sur2_c = sur2;
	    t_i_c = t_i;
	    a1_c = a1;
	    a2_c = a2;
	    e1_c = e1;
	    e2_c = e2;
	    m1_c = m1;
	    m2_c = m2;
	    m3_c = m3;
	    g1_c = g1;
	    g2_c = g2;
	    i1_c = i1;
	    i2_c = i2;
	    i_c = i;
	    beta_c = beta;
	    spin1h_c = spin1h;
	    spintot_c = spintot;
	    vp_c = vp[1][j];
	    spin1e_c = spin1e;
	    spin1q_c = spin1q;
	    spin2e_c = spin2e;
	    spin2h_c = spin2h;
	    R1_c = R1;
	    R2_c = R2;
	    Roche1_c = Roche1;
	    htot_c = htot;
	    type1_c = type1;
	    type2_c = type2;
	    type3_c = type3;
	    beta2_c = beta2;
            gamma_c = gamma;
            gamma2_c = gamma2;
	    }

	  check_time = time(0);
	  duration = check_time-start_time;

	  if ( duration >= 167.0*3600. && count2 <= 2 ){	//prints out parameters after 23 hours and 00 minutes (three times max)  
	    printf("%d \t %d \t %lf \t",sur_c,sur2_c,t_i_c);
	    printf("%le \t %le \t %le \t %le \t %le \t  %le \t  %le   \t  %le \t  %le \t  %le \t  %le \t  %le \t  %le   \t  %le \t  %le \t  %le \t  %le \t  %le \t  %le \t  %le  \t  %le \t %le \t %le \t %le \t %le \t %le \t %lf \t %lf \t %lf \t %le \t %le \t %le \t \n",e1_c,e2_c,g1_c,g2_c,a1_c,i1_c*180./M_PI,i2_c*180./M_PI,i_c*180./M_PI,spin1h_c,spintot_c,beta_c*180./M_PI,vp_c*180/M_PI,spin1e_c,spin1q_c,spin2e_c,spin2q_c,spin2h_c,htot_c,m1_c,R1_c,m2_c,R2_c,a2_c,m3_c,Roche1_c,R1_c/(pow( (m1_c+m2_c)/m2_c,1./3.  )),type1_c,type2_c,type3_c,beta2_c*180./M_PI,gamma_c*180./M_PI,gamma2_c*180./M_PI);

	    printf("%d \t %d \t %lf \t",sur,sur2,t_i);
	    printf("%le \t %le \t %le \t %le \t %le \t  %le \t  %le   \t  %le \t  %le \t  %le \t  %le \t  %le \t  %le   \t  %le \t  %le \t  %le \t  %le \t  %le \t  %le \t  %le  \t  %le \t %le \t %le \t %le \t %le \t %le \t %lf \t %lf \t %lf \t %le \t %le \t %le \t \n",e1,e2,g1,g2,a1,i1*180./M_PI,i2*180./M_PI,i*180./M_PI,spin1h,spintot,beta*180./M_PI,vp[1][j]*180/M_PI,spin1e,spin1q,spin2e,spin2q,spin2h,htot,m1,R1,m2,R2,a2,m3,Roche1,R1/(pow( (m1+m2)/m2,1./3.  )),type1,type2,type3,beta2*180./M_PI,gamma*180./M_PI,gamma2*180./M_PI);

	    start_time = time(0);
	    count2 += 1;
	    }*/

	  if ( a1*(1.-e1) <= 0.02  && count2 <= 2 ){        //prints out parameters if planet periastron <0.02 AU (max three times)
            printf("%d \t %d \t %lf \t",sur,sur2,t_i);
	    printf("%le \t %le \t %le \t %le \t %le \t  %le \t  %le   \t  %le \t  %le \t  %le \t  %le \t  %le \t  %le   \t  %le \t  %le \t  %le \t  %le \t  %le \t  %le \t  %le  \t  %le \t %le \t %le \t %le \t %le \t %le \t %lf \t %lf \t %lf \t %le \t %le \t %le \t \n",e1,e2,g1,g2,a1,i1*180./M_PI,i2*180./M_PI,i*180./M_PI,spin1h,spintot,beta*180./M_PI,vp[1][j]*180/M_PI,spin1e,spin1q,spin2e,spin2q,spin2h,htot,m1,R1,m2,R2,a2,m3,Roche1,R1/(pow((m1+m2)/m2,1./3.)),type1,type2,type3,beta2*180./M_PI,gamma*180./M_PI,gamma2*180./M_PI);    
            count2 += 1; 
	  }

	  /* Smadar :: example for Isable::*/

	  /* print statement for time=0.8*t_sys*/
	  if ( xp[j]<=6.3*0.8e9 && flagtime1==0){

            printf("%d \t %d \t %lf \t",sur,sur2,t_i);
	    printf("%le \t %le \t %le \t %le \t %le \t  %le \t  %le   \t  %le \t  %le \t  %le \t  %le \t  %le \t  %le   \t  %le \t  %le \t  %le \t  %le \t  %le \t  %le \t  %le  \t  %le \t %le \t %le \t %le \t %le \t %le \t %lf \t %lf \t %lf \t %le \t %le \t %le \t \n",e1,e2,g1,g2,a1,i1*180./M_PI,i2*180./M_PI,i*180./M_PI,spin1h,spintot,beta*180./M_PI,vp[1][j]*180/M_PI,spin1e,spin1q,spin2e,spin2q,spin2h,htot,m1,R1,m2,R2,a2,m3,Roche1,R1/(pow((m1+m2)/m2,1./3.)),type1,type2,type3,beta2*180./M_PI,gamma*180./M_PI,gamma2*180./M_PI);    
 
flagtime1 =1;          
}

	  /* print statement for time=0.8*t_sys*/
          if ( xp[j]<=6.3*1e9 && flagtime2==0){

            printf("%d \t %d \t %lf \t",sur,sur2,t_i);
            printf("%le \t %le \t %le \t %le \t %le \t  %le \t  %le   \t  %le \t  %le \t  %le \t  %le \t  %le \t  %le   \t  %le \t  %le \t  %le \t  %le \t  %le \t  %le \t  %le\
  \t  %le \t %le \t %le \t %le \t %le \t %le \t %lf \t %lf \t %lf \t %le \t %le \t %le \t \n",e1,e2,g1,g2,a1,i1*180./M_PI,i2*180./M_PI,i*180./M_PI,spin1h,spintot,beta*180./M_P\
		   I,vp[1][j]*180/M_PI,spin1e,spin1q,spin2e,spin2q,spin2h,htot,m1,R1,m2,R2,a2,m3,Roche1,R1/(pow((m1+m2)/m2,1./3.)),type1,type2,type3,beta2*180./M_PI,gamma*180./M_PI,gamma2*180./M\
		   _PI);

	    flagtime2 =1;
	  }


	  /* print statement for time=0.8*t_sys*/
          if ( xp[j]<=6.3*1.2e9 && flagtime3==0){

            printf("%d \t %d \t %lf \t",sur,sur2,t_i);
            printf("%le \t %le \t %le \t %le \t %le \t  %le \t  %le   \t  %le \t  %le \t  %le \t  %le \t  %le \t  %le   \t  %le \t  %le \t  %le \t  %le \t  %le \t  %le \t  %le\
  \t  %le \t %le \t %le \t %le \t %le \t %le \t %lf \t %lf \t %lf \t %le \t %le \t %le \t \n",e1,e2,g1,g2,a1,i1*180./M_PI,i2*180./M_PI,i*180./M_PI,spin1h,spintot,beta*180./M_P\
		   I,vp[1][j]*180/M_PI,spin1e,spin1q,spin2e,spin2q,spin2h,htot,m1,R1,m2,R2,a2,m3,Roche1,R1/(pow((m1+m2)/m2,1./3.)),type1,type2,type3,beta2*180./M_PI,gamma*180./M_PI,gamma2*180./M\
		   _PI);

	    flagtime3 =1;
	  }





	  if (	  i1*180./M_PI <imin){
	    imin=i1*180./M_PI;
	  }
	  if( i1*180./M_PI>imax){
	    imax=i1*180./M_PI;
	  }


	  if (	  beta*180./M_PI <betamin){
	    betamin=beta*180./M_PI;
	  }
	  if( beta*180./M_PI>betamax){
	    betamax=beta*180./M_PI;
	  }

	  if (	  i*180./M_PI <imumin){
	    imumin=i*180./M_PI;
	  }
	  if( i*180./M_PI>imumax){
	    imumax=i*180./M_PI;
	  }
	  if (e1>e1max){
	    e1max=e1;
	  }
	  if(e1<e1min){
	    e1min=e1;
	  }
	  if (e2>e2max){
	    e2max=e2;
	  }
	  if(e1<e2min){
	    e2min=e2;
	  }

	if( ((t_f-t_i)/t_inter)*((double)(finish - start))/CLOCKS_PER_SEC >= 7200 ){
	    t_quad = 2.*M_PI*pow(a2,3)*pow(1.0-pow(e2,2),1.5)*sqrt(m1+m2)/(pow(a1,1.5)*m3*k2);
	    t_GR = 2*M_PI*pow(a1,2.4)*c2*(1-pow(e1,2))/(3*pow(k2*sqrt(m1+m2),3));
	    //t_TF = (1/(9*365.0))*pow(a1/R1,8)*pow(m1,2)/((m1+m2)*m2);
	    t_TF1 = pow(a1,13./2.)*m2*pow(1-(e1*e1),5)/(sqrt(k2)*kL2*(1.+1.5*(e1*e1)+pow(e1,4)/8.)*m1*sqrt(m1+m2)*pow(R2,5));
	    t_TF2 = pow(a1,13./2.)*m1*pow(1-(e1*e1),5)/(sqrt(k2)*kL1*(1.+1.5*(e1*e1)+pow(e1,4)/8.)*m2*sqrt(m1+m2)*pow(R2,5));
	    //fprintf(stderr, "Needed time is %le more days!\n",(((t_f-t_i)/t_inter)*((double)(finish - start))/CLOCKS_PER_SEC)/86400. );
	    if(t_TF1 < t_quad || t_TF2 < t_quad || tF1 < t_quad || tF2 < t_quad){
	      //fprintf(stderr, "Tidal, GR, Kozai time scales = %le, %le, %le years\n",t_TF,t_GR,t_quad );
	      //if(e1<=1.e-5){
		//if(a1/Roche1 < 3. && a1/Roche1 > 1.5){
		P1 = 2*M_PI*sqrt(a1*a1*a1/(k2*(m1+m2)));
		if( e1 <= 1e-5  ){
		//if( (a1/Roche1 <= 5 || a1/Roche2 <= 5 || a1 <= 0.1) && e1 <= 1e-2 && ( e1 <= 1e-4  || ((fabs((2*2*M_PI/(P1*(cos(beta)+pow(cos(beta),-1))))/spintot -1.) <= 1e-3 || fabs((2*2*M_PI/(P1*(cos(acos(spin2h/spintot2))+pow(cos(acos(spin2h/spintot2)),-1))))/spintot2 -1.) <= 1e-3)) ) || (a1 <= 0.1 && e1 <= 1e-2) ){
		  //if( (e1 <= 1e-3 && (a1/Roche1 <= 2. || a1/Roche2 <= 2.)) || fabs((2*2*M_PI/(P1*(cos(beta)+pow(cos(beta),-1))))/spintot -1.) <= 1e-3 || fabs((2*2*M_PI/(P1*(cos(acos(spin2h/spintot2))+pow(cos(acos(spin2h/spintot2)),-1))))/spintot2 -1.) <= 1e-3 ){
		  flgout=1;
		  sur=4;
		  //fprintf(stderr, "Tidal locking at t = %le years\n",t_i );
		}
		/*else if(a1/Roche2 < 3. && a1/Roche2 > 1.5){
		  flgout=1;
		  sur=4;
		  fprintf(stderr, "Tidal locking at t = %le years\n",t_i );
		  }*/
		//}
	    }
	  }

	if ( j == k_f && xp[k_f]+tMS >= toutput[count] ){
          printf("%d \t %d \t %lf \t",sur,sur2,t_i);
          printf("%le \t %le \t %le \t %le \t %le \t  %le \t  %le   \t  %le \t  %le \t  %le \t  %le \t  %le \t  %le   \t  %le \t  %le \t  %le \t  %le \t  %le \t  %le \t  %le  \t  %le \t %le \t %le \t %le \t %le \t %le \t %lf \t %lf \t %lf \t %le \t  %le \t  %le \t \n",e1,e2,g1,g2,a1,i1*180./M_PI,i2*180./M_PI,i*180./M_PI,spin1h,spintot,beta*180./M_PI,vp[1][j]*180/M_PI,spin1e,spin1q,spin2e,spin2q,spin2h,htot,m1,R1,m2,R2,a2,m3,Roche1,R1/(pow( (m1+m2)/m2,1./3.  )),type1,type2,type3,beta2*180./M_PI,gamma*180./M_PI,gamma2*180./M_PI);

          count++;
        }


	}// end of  priting loop on the steps that were saved for odeint

	t_b=t;
	if (flgout==1) break;
	//	exit(5);

	//initialize again for the next step of integration::
	x_i[1] = yp[1][k_f];
	x_i[2] =  yp[2][k_f];
	x_i[3] =  yp[3][k_f];
	x_i[4] = yp[4][k_f];
	x_i[5]=yp[5][k_f];
	x_i[6]=yp[6][k_f];
	x_i[7]=yp[7][k_f];
	x_i[8]=yp[8][k_f]; 
	x_i[9]=yp[9][k_f];
	x_i[10]=yp[10][k_f];
	x_i[11]=yp[11][k_f];
	x_i[12]=yp[12][k_f];
	x_i[13]=yp[13][k_f];
	x_i[14]=yp[14][k_f];
	x_i[15]=R1;//yp[15][k_f];
	x_i[16]=yp[16][k_f];
	x_i[17]=yp[17][k_f];
	x_i[18]=yp[18][k_f];
	x_i[19]=R2;//yp[19][k_f];
	x_i[20]=yp[20][k_f];

	t_i=xp[k_f];
	t=t_i;
	// parameters:
	v[1]	= i;
	//	v[2] 	= m1;
	//  v[3] 	= m2;
	v[4] 	= m3;
	v[5]	= a1;
	v[6]	= a2;
	//	v[7]    = R1;
	//  v[8]    = R2;

	P1 = 2*M_PI*sqrt(a1*a1*a1/(k2*(m1+m2)));
	P2 = 2*M_PI*sqrt(a2*a2*a2/(k2*(m1+m2+m3)));
	P_koz = (8./(15*M_PI))*((m1+m2+m3)/m3)*(P2*P2/P1)*sqrt(1-e2*e2)*(1-e2*e2);
	t_quad = 2.*M_PI*pow(a2,3)*pow(1.0-pow(e2,2),1.5)*sqrt(m1+m2)/(pow(a1,1.5)*m3*k2);
	t_GR = 2*M_PI*pow(a1,2.5)*c2*(1-pow(e1,2))/(3*pow(k2*sqrt(m1+m2),3));
	t_TF1 = pow(a1,13./2.)*m2*pow(1-(e1*e1),5)/(sqrt(k2)*kL2*(1.+1.5*(e1*e1)+pow(e1,4)/8.)*m1*sqrt(m1+m2)*pow(R2,5));
	t_TF2 = pow(a1,13./2.)*m1*pow(1-(e1*e1),5)/(sqrt(k2)*kL1*(1.+1.5*(e1*e1)+pow(e1,4)/8.)*m2*sqrt(m1+m2)*pow(R2,5));
	tF1 = TV1*pow(a1/R1,8)*SQR(m1)/ ( 9.*(m1+m2)*m2*SQR(1+2.*kL1));
	tF2 = TV2*pow(a1/R2,8)*SQR(m2)/ ( 9.*(m1+m2)*m1*SQR(1+2.*kL2));
	t_rot1 = sqrt(k2)*pow(a1,7./2.)*pow(1.-e1*e1,2.)/(kL1*pow(spintot,2)*sqrt(m1+m2)*pow(R1,5));
	t_rot2 = sqrt(k2)*pow(a1,7./2.)*pow(1.-e1*e1,2.)/(kL2*pow(spintot2,2)*sqrt(m1+m2)*pow(R2,5));


	/*if (3.*R1 >= a1 || 3.*R2 >= a1){
	  if (dhtotdtTF[1] != 0. && detotdtTF[1] != 0.){
	    //fprintf(stderr," Are you in here?\n");
	    if (!(SSE1)) linint(tofdmdtvec,dRdtvec,(Ndis+1)*Nvar,xp[j]+tMS,&dR1dt); 
	    else dR1dt = 0.;
	    //fprintf(stderr," Are you in here2?\n");
	    if (!(SSE2)) linint(tofdmdtvec2,dRdtvec2,(Ndis2+1)*Nvar,xp[j]+tMS,&dR2dt); 
	    else dR2dt = 0.;
	    //fprintf(stderr," Are you in here4?\n");
	    if (fabs(dR1dt/R1)/2. <= fabs((2.*dhtotdtTF[1]/htot) + (2.*detotdtTF[1]*e1/(1-e1*e1)))  ||  (fabs(dR2dt/R2)/2. <= fabs((2*dhtotdtTF[1]/htot) + (2*detotdtTF[1]*e1/(1-e1*e1))))) {
	      if (P2 <= 100.) t_inter = 100.;
	      else t_inter = P2;
	      //t_inter = 25.;
	      //fprintf(stderr," Are you in here?\n");
	    }
	    else{ 
	      if (P_koz/250.<=fabs(t_f-t_ini)/400000){
		t_inter=P_koz/250.;
	      }
	      else t_inter=fabs(t_f-t_ini)/400000;
	    }
	  }
	}
	else{ 
	  if (P_koz/250.<=fabs(t_f-t_ini)/400000){
	    t_inter=P_koz/250.;
	  }
	  else t_inter=fabs(t_f-t_ini)/400000;
	}

	if (tF1/100. <= t_quad || tF2/100. <= t_quad){
	  if (t_inter > tF1/250.){
	    if (P1 < tF1/250.) t_inter=tF1/250.;
	    else t_inter = P1;
	  } 
	  else if (t_inter > tF2/250. && tF2 < tF1){
	    if (P1 < tF2/250.) t_inter=tF2/250.;
	    else t_inter = P1;
	  } 
	  }*/

	if (P_koz/100.<=fabs(t_f-t_ini)/400000){
	  t_inter=P_koz/100.;
	}
	else t_inter=fabs(t_f-t_ini)/400000;

	/*if (t_TF1 < t_inter || t_TF2 < t_inter || tF1 < t_inter || tF2 < t_inter || t_rot1 < t_inter || t_rot2 < t_inter){
	  if (t_rot1 < t_TF1 && t_rot1 < t_TF2 && t_rot1 < t_rot2 && t_rot1 < tF1 && t_rot1 < tF2) t_inter = fabs(t_rot1);
	  else if (t_rot2 < t_TF1 && t_rot2 < t_TF2 && t_rot2 < t_rot1 && t_rot2 < tF1 && t_rot2 < tF2) t_inter = fabs(t_rot2);
	  else if (t_TF1 < t_rot1 && t_TF1 < t_TF2 && t_TF1 < t_rot2 && t_TF1 < tF1 && t_TF1 < tF2) t_inter = fabs(t_TF1);
	  else if (t_TF2 < t_rot1 && t_TF2 < t_TF1 && t_TF2 < t_rot2 && t_TF2 < tF1 && t_TF2 < tF2) t_inter = fabs(t_TF2);
	  else if (tF1 < t_rot1 && tF1 < t_TF2 && tF1 < t_rot2 && tF1 < t_TF1 && tF1 < tF2) t_inter = fabs(tF1);
	  else if (tF2 < t_rot1 && tF2 < t_TF2 && tF2 < t_rot2 && tF2 < t_TF1 && tF2 < tF1) t_inter = fabs(tF2);

	  if (t_inter < 100.) t_inter = 100.;
	}
	else{
	  if (P_koz/250.<=fabs(t_f-t_ini)/400000){
	    t_inter=P_koz/250.;
	  }
	  else t_inter=fabs(t_f-t_ini)/400000;
	  }*/

	//if (t_TF1 < t_quad/5. || t_TF2 < t_quad/5. || tF1 < t_quad/5. || tF2 < t_quad/5.){
	//  OCTUPOLE = QUADRUPOLE = 1;
	//}
	//else OCTUPOLE = QUADRUPOLE = 0;

	}//end of integrating loop
	j=k_f;

	if (!(SSE1)){
	  if (xp[j]+tMS>tdis[Ndis]){
	    type1=typeSSE[indexdis[Ndis]+1];
	          
	  }
	  else {
	    if (xp[j]+tMS<timeSSE[2]) type1=typeSSE[1];
	    else {
	      //fprintf(stderr, " Are you here type 1\n");
	      linint(timeSSE,typeSSE,Nlin,xp[j]+tMS,&type1);
	    }
	    //    fprintf(stderr, " :: type1=%le\n",type1);
	  }
	}

	if (!(SSE2)){
	  if (xp[j]+tMS>tdis2[Ndis2]){
	    type2=typeSSE2[indexdis2[Ndis2]+1];
	          
	  }
	  else {
	    if (xp[j]+tMS<timeSSE2[2]) type2=typeSSE2[1];
	    else {
	      //fprintf(stderr, " Are you here type 1\n");
	      linint(timeSSE2,typeSSE2,Nlin2,xp[j]+tMS,&type2);
	    }
	    //    fprintf(stderr, " :: type1=%le\n",type1);
	  }
	}
	//  fprintf(stderr, " :: tdis[Ndis]=%le xp[j]=%le type1=%e\n",tdis[Ndis],xp[j],type1);
	//  exit(5);
	if (!(SSE3)){
	  if (xp[j]+tMS>tdis3[Ndis3]){
	    type3=typeSSE3[indexdis3[Ndis3]+1];
	  }
	  else {
	    if (xp[j]+tMS<timeSSE3[2]) type3=typeSSE3[1];
	    else {
	      //fprintf(stderr, " Are you here type 3\n");
	      linint(timeSSE3,typeSSE3,Nlin3,xp[j]+tMS,&type3);
	    }
	    //    findt3(xp[j]+tMS,&minIn,&maxIn);
	    //    linintMinMax(timeSSE,typeSSE,minIn,maxIn,xp[j]+tMS,&type3);
	  }
	      
	}
	//  fprintf(stderr, " :: tdis3[Ndis]=%le xp[j]=%le type3=%e\n",tdis3[Ndis3],xp[j],type3);
	//  exit(4);
	//  fprintf(stderr, " :: tdis[Ndis]=%le xp[j]=%le type1=%e\n",tdis[Ndis],xp[j],type1);
	// if (t-t_b >=Tprint){
	printf("%d \t %d \t %lf \t",sur,sur2,t);
	printf("%le \t %le \t %le \t %le \t %le \t  %le \t  %le   \t  %le \t  %le \t  %le \t  %le \t  %le \t  %le   \t  %le \t  %le \t  %le \t  %le \t  %le \t  %le \t  %le  \t  %le \t %le \t %le \t %le \t %le \t %le \t %lf \t %lf \t %lf \t %le \t %le \t %le \t \n",e1,e2,g1,g2,a1,i1*180./M_PI,i2*180./M_PI,i*180./M_PI,spin1h,spintot,beta*180./M_PI,vp[1][j]*180/M_PI,spin1e,spin1q,spin2e,spin2q,spin2h,htot,m1,R1,m2,R2,a2,m3,Roche1,R1/(0.462 * pow( m1/m2,1./3.  )),type1,type2,type3,beta2*180./M_PI,gamma*180./M_PI,gamma2*180./M_PI);


	//	}
	fprintf(stderr,"	    -------\n");

	fprintf(stderr,"	    i1max=%f i1min=%f \n",imax,imin);
	fprintf(stderr,"	    mutual: imax=%f imin=%f \n",imumax,imumin);
	fprintf(stderr,"	            betamax=%f betamin=%f \n",betamax,betamin);


	fprintf(stderr,"	    emin=%f emax=%f\n",e1min,e1max);

	//fclose(out);

	/* free memory */
	free_vector(x,1,2);
	free_vector(x_i,1,2);
	free_vector(dxdt,1,2);
	free_vector(v,1,n_orb_param);
	free_vector(ve,1,n_orb_param);
	free_vector(xp,1,kmax);
	free_matrix(yp,1,kmax,1,kmax);
	free_matrix(vp,1,n_orb_param,1,kmax);
	free_vector(dmdtvec,1,Nlinmax);
	free_vector(tofdmdtvec,1,Nlinmax);
	free_vector(dmdtvec2,1,Nlinmax);
	free_vector(tofdmdtvec2,1,Nlinmax);
	free_vector(dRdtvec,1,Nlinmax);
	free_vector(dRdtvec2,1,Nlinmax);
	free_vector(MtSSE,1,Nlinmax);
	free_vector(timeSSE,1,Nlinmax);
	free_vector(MtSSE2,1,Nlinmax);
	free_vector(timeSSE2,1,Nlinmax);
	free_vector(MtSSE3,1,Nlinmax);
	free_vector(timeSSE3,1,Nlinmax);
	free_vector(mvec,1,Nlinmax);
	free_vector(Rvec,1,Nlinmax);

	free_vector(detotdtTF,1,1);
	free_vector(dhtotdtTF,1,1);

	return 0;
}

