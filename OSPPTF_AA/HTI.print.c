/****************************************************************************
 * Integrate octupole perturbation equations for heirarchical triple systems 
 * Takes input from "triple.in" 
 * * smadar: adding TF May 2010                                                    
 ****************************************************************************/


#include <stdio.h> 
#include <stdlib.h> 
#include <string.h> 
#include <math.h>
#include "nrutil.h"
#include "HTI.h"

#define TINY 1e-30
/* odeint control parameters */
int kmax = 10000000,    kount = 0;
double *v, *xp, **yp, dxsav,**vp,*ve;						/* storage arrays */
//smadar: adding a storage matrix vp for v parameters 
int flgout;
/* main program, calls for odeint.c */
int main (void) {
	
	double 	*x, *x_i, *dxdt;
	double 	t_i, t_f, eps=0, h1=0, hmin=0;
	int		nok=0, nbad=0, j=0, k_f=0;
	int		scnret;
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
	double beta,spintot,spinDh;
	double t_b,tloop_f,t_inter;
	double QoftV1,QoftV2;
	double kL1=0.5*Q1/(1.-Q1);
	
	double kL2=0.5*Q2/(1.-Q2);

	sur=2;
	//	spin1=2*M_PI*365.25/Sp1;
	//	spin2=2*M_PI*365.25/Sp2;

	//sprintf(outputname,"Angular.txt");
	//	out=fopen(outputname,"w");

	flgout=0;
	/* takes in input from input.txt */
	input = fopen("triple.in", "r");
	if (ferror(input)!=0) {
	    printf("Error from input file \n");
	    return 0;
	}

	for(;;) {	/* read in until :::, end of header */
		scnret = fscanf(input, "%s", dummy);
		if (!(strcmp(dummy, ":::"))) break;
	}

	fscanf(input, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf ",\
	       &m1,&m2,&m3,&R1,&R2,&spin1,&spin2,&beta,&a1,&a2,&e1,&e2,&g1,&g2,&i,&t_f);


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
	//	fprintf(stderr," %d\n",GR);


	for(;;) {	/* read in until :::, end of header */
		scnret = fscanf(input, "%s", dummy);
		if (!(strcmp(dummy, ":::"))) break;
	}

	fscanf(input, "%d",&TF);
	//	fprintf(stderr," %d\n",TF);

	for(;;) {	/* read in until :::, end of header */
		scnret = fscanf(input, "%s", dummy);
		if (!(strcmp(dummy, ":::"))) break;
	}
	
/* read in control parameters */

	fscanf(input, "%lf",&eps);

      	fclose(input);

	if (TF==1 && QUADRUPOLE==1 && OCTUPOLE==1 && GR==1){
	  fprintf(stderr, "   =====  oops....you did not choose any physical effect to explor =====\n");
	  fprintf(stderr, "   ===== now terminating.... bye bye \n");
	  exit(30);
	}

	if (e1==0){
	  e1=TINY;
	}

	ideg=i;
	i *=M_PI/180.0;
	beta *=M_PI/180.0;
	g1*=M_PI/180.0;
	g2*=M_PI/180.0;
	R1*=Rsun;
	R2*=Rsun;
	e1max_calc = sqrt(1-5/3.0*cos(i)*cos(i));
	t_f *= 1e6;
	t_i = 0;
	t_b=t_i-Tprint-1;
	//	t_inter=fabs(t_f-t_i)/400000;//4000;//10000;//100.;// /100000;//100000.;
	//	t_inter=fabs(t_f-t_i)/400000;//10000;//100.;// /100000;//100000.;
	//	t_inter=fabs(t_f-t_i)/4000;//10000;//100.;// /100000;//100000.;
	t_inter=fabs(t_f-t_i)/40000;//10000;//100.;// /100000;//100000.;
 	spin1=2.*M_PI*365.25/spin1;
	spin2=2.*M_PI*365.25/spin2;


	/* various timescales */
	P1 = 2*M_PI*sqrt(a1*a1*a1/(k2*(m1+m2)));
	P2 = 2*M_PI*sqrt(a2*a2*a2/(k2*(m1+m2+m3)));
	P_koz = P1*(m1+m2)/m3*pow(a2/a1,3)*sqrt(1-e2*e2)*(1-e2*e2);
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
	else if(theta1<1)      i1 = M_PI;

	if(fabs(theta2)<=1)    i2=acos(theta2);
	else if(theta2>1)      i2 = 0.0;
	else if(theta2<1)      i2 = M_PI;


	fprintf(stderr,"H1=%le H2=%le\n",G1*theta1,G2*theta2);
	fprintf(stderr,"G1=%le G2=%le\n",G1,G2);
fprintf(stderr,"theta1=%le theta2=%le\n",theta1,theta2);
//exit(1);
	////

	///c
	htot=sqrt(k2 * (m1+m2) *a1 *(1-SQR(e1)));
	spin1h=spin1*cos(beta);
	// damping the rest in the two other directions.
	spin1e=spin1q =sqrt(( SQR(spin1)-SQR(spin1h))/2 );

	//assuming theat the planet is alline
	spin2h=spin2;
	spin2e=spin2q =0;

	spin1x=2*M_PI*365.25/spin1;
	spin1y=spin1z=0;
	spin2x=spin1x;
	spin2y=spin2z=0;
	
	/* memory allocation */
	x = 	vector(1,n_eq);
	x_i = 	vector(1,n_eq);
	dxdt = 	vector(1,n_eq);
	v = 	vector(1,n_orb_param);
	xp = 	vector(1,kmax);
	yp 	=	matrix(1,n_eq,1,kmax);
	vp = matrix(1,n_orb_param,1,kmax); //smadar: define vp
	ve =	vector(1,n_orb_param);//smadar allocating vexternal saving
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
	Roche=R2*pow(m2/(m2+m1),-1./3.)/0.6;
	//	fprintf(stderr, " \n The 2*Roche limit is %le AU also note that 0.462*q^1/3*a =%le\n",Roche,0.462*pow(m2/m1,1./3.)*R1);
	fprintf(stderr, " \n The 2*Roche limit is %le AU\n",Roche);
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

	/* integrate the octupole equations (use bsstep rather than rkqs) */
	
	//t_i=TSTART;

	while (t_i<=t_f && sur==2){
	  tloop_f=t_i+t_inter;
	  //fprintf(stderr, " h1=%le eps=%le\n",h1,eps);
	  //fprintf(stderr, " Integrating on: tstart=%le tend=%le\n",t_i,tloop_f);
	  	  h1=fabs(tloop_f-t_i)/100000.;///10000.;//*eps;///10;
	  //	  h1=P2/100.;
	odeint(x_i, n_eq, t_i, tloop_f, eps, h1, hmin, &nok, &nbad, octupole, integrator);
	k_f = kount;
	//	fprintf(stderr, "====\n");
	for(j=1;j<=k_f;j++) {
	  
	  e1 = sqrt(yp[1][j]*yp[1][j] + yp[2][j]*yp[2][j]);
	  e2 = sqrt(yp[3][j]*yp[3][j] + yp[4][j]*yp[4][j]);
	  g1 = atan2(yp[1][j],yp[2][j]);
	  g2 = atan2(yp[3][j],yp[4][j]);
	  t=xp[j];

	  htot=yp[8][j];

	  a1= SQR(htot)/ (k2 * (m1+m2) *(1-SQR(e1)));

	  spin1h=yp[11][j];

	  //spin1h = spin*hat(h)*cos(beta)

	  spintot=sqrt(SQR(spin1h)+SQR(yp[10][j])+SQR(yp[9][j])  );
	  spinDh =  spin1h / (spintot);//spin1h / (spintot * htot);


	  if(fabs(spinDh)<=1)    beta=acos(spinDh);
	  else if(spinDh>1)      beta = 0.0;
	  else if(spinDh<1)      beta = M_PI;

	  spin1e=yp[9][j];
	  spin1q=yp[10][j];
	  spin2e=yp[12][j];
	  spin2q=yp[13][j];

	  spin2h=yp[14][j];
	  
	  if(fabs(yp[5][j])<=1)    i1=acos(yp[5][j]);
	  else if(yp[5][j]>1)      i1 = 0.0;
	  else if(yp[5][j]<1)      i1 = M_PI;


	  if(fabs(yp[6][j])<=1)    i2=acos(yp[6][j]);
	  else if(yp[6][j]>1)      i2 = 0.0;
	  else if(yp[6][j]<1)      i2 = M_PI;

	  if(fabs(yp[7][j])<=1)    i=acos(yp[7][j]);
	  else if(yp[7][j]>1)      i = 0.0;
	  else if(yp[7][j]<1)      i = M_PI;

	  imu=i1+i2;


	  L1 = m1 * m2  / m12  * sqrt( k2 * m12  * a1 );
	  L2 = m3 * m12 / m123 * sqrt( k2 * m123 * a2 );
	  G1 = L1 * sqrt(1 - e1*e1);
	  G2 = L2 * sqrt(1 - e2*e2);
	  Hsq = G1*G1 + 2.*G1*G2*cos(i) + G2*G2;
	  H1 = G1*yp[5][j];
	  H2 = G2*yp[6][j];

	  // if (t-t_b >=Tprint){
	    printf("%lf\t",xp[j]);
	    printf("%1.6lf \t %1.6lf \t %1.6lf \t %1.6lf  %1.6lf \t  %1.6lf \t  %1.6lf   \t  %1.6lf \t  %1.6lf \t  %1.6lf \t  %1.6lf \t  %1.6lf \t  %1.6lf   \t  %1.6lf \t  %1.6lf \t  %1.6lf \t  %1.6lf \t  %1.6lf\n",e1,e2,g1,g2,a1,i1*180./M_PI,i2*180./M_PI,i*180./M_PI,spin1h,spintot,beta*180./M_PI,vp[1][j]*180/M_PI,spin1e,spin1q,spin2e,spin2q,spin2h,htot);
	    if (a1<=Roche){
              fprintf(stderr, " ::::: a1(=%le AU) < Roche (=%le AU) ::::: \n",a1,Roche);
              fprintf(stderr, " The planet just died....poor thing....\n");
              fprintf(stderr, " Terminating now...\n");
              exit(7);
            }

	    /*
	    if (a1*(1-e1)<=Roche){
              fprintf(stderr, " ::::: a1*(1-e1)(=%le AU) < Roche (=%le AU) ::::: \n",a1*(1-e1),Roche);
              fprintf(stderr, " The planet just died....poor thing....\n");
              fprintf(stderr, " Terminating now...\n");
              exit(7);
            }
	    */
	    //fprintf(out,"%le %le %le %le %le %le %le %le\n",xp[j],H1,H2,Hsq,G1,G2,L1,L2);

	    // fprintf(stderr,"%le %le %le %le %le %le %le %le\n",xp[j],H1,H2,Hsq,G1,G2,L1,L2);
 // exit(1);

	    // }
	    //fprintf(stderr, " xp[%d]=%le\n",j,xp[j]);

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


	}// end of  priting loop on the steps that were saved for odeint
	t_b=t;
	 

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

	t_i=xp[k_f];
	// parameters:
	v[1]	= i;
	v[2] 	= m1;
	v[3] 	= m2;
	v[4] 	= m3;
	v[5]	= a1;
	v[6]	= a2;
	v[7]    = R1;
	v[8]    = R2;


	}//end of integrating loop
	j=k_f;
	//fprintf(stderr,"%lf\t",xp[j]);
	//	    fprintf(stderr,"%1.6lf \t %1.6lf \t %1.6lf \t %1.6lf  %1.6lf \t  %1.6lf \t  %1.6lf   \t  %1.6lf \t  %1.6lf \t  %1.6lf \t  %1.6lf \t  %1.6lf\n",e1,e2,g1,g2,a1,i1*180./M_PI,i2*180./M_PI,i*180./M_PI,spin1h,spintot,beta*180./M_PI,vp[1][j]*180/M_PI);

	/* output result */
	/*
	fprintf(stderr,"	\n");

	for(j=1;j<=k_f;j++) {
	  
	  printf("%lf\t",xp[j]);
	  
	  e1 = sqrt(yp[1][j]*yp[1][j] + yp[2][j]*yp[2][j]);
	  e2 = sqrt(yp[3][j]*yp[3][j] + yp[4][j]*yp[4][j]);
	  g1 = atan2(yp[1][j],yp[2][j]);
	  g2 = atan2(yp[3][j],yp[4][j]);
	  t=xp[j];



	

	L1 = m1 * m2  / m12  * sqrt( k2 * m12  * a1 );
	L2 = m3 * m12 / m123 * sqrt( k2 * m123 * a2 );
	G1 = L1 * sqrt(1 - e1*e1);
	G2 = L2 * sqrt(1 - e2*e2);
       	if(t==0)    Hsq = G1*G1 + 2.*G1*G2*cos(i) + G2*G2;
	theta = (Hsq - G1*G1 - G2*G2)/(2.0*G1*G2);
	if(fabs(theta)<=1)    tempi=acos(theta);
	else if(theta>1)      tempi = 0.0;
	else if(theta<1)      tempi = M_PI;



	htot=yp[8][j];

	a1= SQR(htot)/ (k2 * (m1+m2) *(1-SQR(e1)));

	spin1h=yp[11][j];

	//spin1h = spin*hat(h)*cos(beta)

	spintot=sqrt(SQR(spin1h)+SQR(yp[10][j])+SQR(yp[9][j])  );
	spinDh =  spin1h / (spintot);//spin1h / (spintot * htot);
	//	fprintf(stderr,"%1.6lf \t %1.6lf \t %1.6lf \t  %1.6lf \t  %1.6lf \t  %1.6lf\n",yp[9][j],yp[10][j],spin1h,htot,spintot,spinDh);
	if(fabs(spinDh)<=1)    beta=acos(spinDh);
	else if(spinDh>1)      beta = 0.0;
	else if(spinDh<1)      beta = M_PI;


	spin2h=yp[14][j];

	if(fabs(yp[5][j])<=1)    i1=acos(yp[5][j]);
	else if(yp[5][j]>1)      i1 = 0.0;
	else if(yp[5][j]<1)      i1 = M_PI;


	if(fabs(yp[6][j])<=1)    i2=acos(yp[6][j]);
	else if(yp[6][j]>1)      i2 = 0.0;
	else if(yp[6][j]<1)      i2 = M_PI;



	if(fabs(yp[7][j])<=1)    i=acos(yp[7][j]);
	else if(yp[7][j]>1)      i = 0.0;
	else if(yp[7][j]<1)      i = M_PI;

	imu=i1+i2;

	printf("%1.6lf \t %1.6lf \t %1.6lf \t %1.6lf  %1.6lf \t  %1.6lf \t  %1.6lf   \t  %1.6lf \t  %1.6lf \t  %1.6lf \t  %1.6lf \t  %1.6lf\n",e1,e2,g1,g2,a1,i1*180./M_PI,i2*180./M_PI,i*180./M_PI,spin1h,spintot,beta*180./M_PI,vp[1][j]*180/M_PI);


	//	fprintf(stderr,"%1.6lf \t %1.6lf \t %1.6lf \t %1.6lf  %1.6lf \t  %1.6lf \t  %1.6lf \t  %1.6lf \t  %1.6lf \t  %1.6lf \n",e1,e2,g1,g2,yp[5][j],acos(yp[6][j])*180./M_PI,acos(yp[7][j])*180./M_PI,acos(yp[6][j])*180./M_PI+acos(yp[7][j])*180./M_PI,yp[8][j],yp[9][j]);

	//]]if (j==12)	exit(1);

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

	}
	*/
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
	return 0;
}

