#include <math.h>
#include "HTI.h"
#include <stdio.h>
#include <stdlib.h>
#include "nrutil.h"


extern double *v,*ve;
extern int flgout;
extern double *tdis, *tdis2, *tdis3;
extern double R1f,R2f;

extern 	int Ndis,Nlin,Ndis2,Nlin2,Ndis3,Nlin3;
extern double *timeSSE,*MenvSSE,*timeSSE2,*MenvSSE2,*timeSSE3,*MtSSE3,*indexdis3,*spinSSE,*spinSSE2;
extern double type1, type2, type3;
extern double t_inter;
extern double *detotdtTF,*dhtotdtTF;
/* ODE for a damped harmonic oscillator to be integrated */
void octupole (double t, double *x, double *dxdt) {

  //static double  Hsq;
  double L1,L2, C3G23, C4G25, m12, m123,m1,a1,a2,P1,P2,R1,m3,m2,R2,Hsq;
	static int	  flg=0;
	int j;
	double e1sing1, e1cosg1, e2sing2, e2cosg2;
	double de1sing1dt, de1cosg1dt, de2sing2dt, de2cosg2dt;
	double g1,g2,e1,e2,i=0.0;                          /* varrying quantities */
	double dg1dt, dg2dt, de1dt, de2dt;                 /* derivatives         */
	double G1,G2,theta,A,B,C3,C4,cosphi;               /* ang.mom., etc.      */
	double sing1, cosg1, sing2, cosg2, sin2g1, cos2g1; /* trig funcs          */
	//smadar adding varible for calc the inclinations:
	double tF2,tF1,In1,In2,mu,detotdt=0.,didtTF=0;
	double dFdg2,dFdg1,di1dt,di2dt,didt,dBigG1dt=0.,dBigG2dt=0.,di1_onlydt,Ztot,Xtot,Ytot;
	double i1,i2,di1dtTF,theta1,H1;
	double dBigG1dtML=0,dBigG2dtML=0,BigH1dotML=0,BigH2dotML=0,BigHtotdotML=0,di1dtML=0,di2dtML=0,BigHtotdot=0,dJspindtMB=0,dJspin2dtMB=0;
	double BigH1dot=0,BigH2dot=0;
	double Omega1,Omega2;
	double dOmega1dt,dOmega2dt,dOmega1dtTF,htot,dhtotdt,dg1dtTF;
	double n;
	//spins:
	double spin1h,spin2h,spin1q,spin2q,spin1e,spin2e,spin1,spin2,spin1tot,spin2tot,rootval; //spin parameters
	double dspin1edt,dspin1qdt,dspin1hdt,dspin2edt,dspin2qdt,dspin2hdt;
	double X1,X2,Y1,Y2,V1,V2,W1,W2,Z1,Z2,ZGR; //TF cooeficents
	double e1_2,e1_4,e1_6,e1_8; //for TF 
	double sini1,sini;
	double dm1dt,dR1dt,dm2dt,dR2dt;
	double menv=0.,menv2=0.,val;
	void linint(double *xa,double *ya,int n,double x,double *y);
	void linintMinMax(double *xa,double *ya,int minv,int n,double x,double *y);
	double dmdt(double t, double m, double m0);
	double dRdt(double t, double R);
	double dmdt2(double t, double m, double m0);
	double dRdt2(double t, double R);
	double dspindt(double t, double spin);
	double dspindt2(double t, double spin2);
	dspin1edt=dspin1qdt=dspin1hdt=dspin2edt=dspin2qdt=dspin2hdt=0.;
	X1=X2=Y1=Y2=Z1=Z2=W1=W2=V1=V2=0;
	dm1dt=dR1dt=dm2dt=dR2dt=0;
	i1=i2=0.;
	double spin1SSE,spin1SSEold,spin2SSE,spin2SSEold,dspintot1dt,dspintot2dt;

	double kL1=0.5*Q1/(1.-Q1);
	
	double kL2=0.5*Q2/(1.-Q2);
	//	fprintf(stderr,"t=%le\n",t);
	/*
	R1=2.13*695500*6.68459e-9;//
	R2=1.95*71492*6.68459e-9;//au
	*/
	/* copy values from the vector to variables */
	e1sing1 = x[1];
	e1cosg1 = x[2];
	e2sing2 = x[3];
	e2cosg2 = x[4];

	e1 = sqrt(e1sing1*e1sing1+e1cosg1*e1cosg1);
	e2 = sqrt(e2sing2*e2sing2+e2cosg2*e2cosg2);
	g1 = atan2(e1sing1,e1cosg1);
	g2 = atan2(e2sing2,e2cosg2);



	htot = x[8];

	spin1e = x[9];
	spin1q = x[10];
	spin1h = x[11];
	
	spin2e = x[12];
	spin2q = x[13];
	spin2h = x[14];

	spin1=sqrt(spin1e*spin1e+spin1q*spin1q+spin1h*spin1h);
	spin2=sqrt(spin2e*spin2e+spin2q*spin2q+spin2h*spin2h);


	//	fprintf(stderr," spin1e=%le spin1q=%le spin1h=%le spin2e=%le spin2q=%le spin2h=%le\n",spin1e,spin1q,spin1h,spin2e,spin2q,spin2h);
	if(!(flg))		{
		i =v[1];
		m1=v[2];
		m2=v[3];
		m3=v[4];
		a1=v[5];
		a2=v[6];
		R1=v[7];
		R2=v[8];
		m12 = m1 + m2;
		m123 = m1 + m2 + m3;
		
		P1 = 2*M_PI*sqrt(a1*a1*a1/(k2*(m1+m2)));
		P2 = 2*M_PI*sqrt(a2*a2*a2/(k2*(m1+m2+m3)));
    
		L1 = m1 * m2  / m12  * sqrt( k2 * m12  * a1 );
		L2 = m3 * m12 / m123 * sqrt( k2 * m123 * a2 );
    
		C3G23 = k2*k2/16.0 * pow(m12*m3,7)/pow(m123*m1*m2,3) * pow(L1/L2,3) * L1;
		C4G25 = -k2*k2*(15/64.0) * pow(m12*m3,9)*(m1-m2)/pow(m123,4)/pow(m1*m2,5) * pow(L1*L1/L2,3);
		
		/* debug */
		DebugPrint2("L1=%g L2=%g  ",L1, L2);
		DebugPrint2("C3G23=%g C4G25=%g\n",C3G23,C4G25);
	
		flg++;
	}

	//smadar adding mass loss:
	/*	if (!(SSE) && R1<R1f) {
	  fprintf(stderr," octupole: input x[15]=%le\n",x[15]);
	  x[15]=R1f;
	  fprintf(stderr," octupole: R1=%le AU which is %le Rsun\n",x[15],x[15]/Rsun);
	  }*/
	R1=x[15];

	m1=x[16];
	G1=x[17];
	G2=x[18];

	R2=x[19];
	m2=x[20];	
	
	if (!(SSE3)){
	  // minIndex=indexdis3[Ndis3-1];
	  // maxIndex=indexdis3[Ndis3];
	  //  linintMinMax(timeSSE3,MtSSE3,Nlin3,t+tMS,&val);
	  if (t+tMS<timeSSE3[2]) val=MtSSE3[1];
	  linint(timeSSE3,MtSSE3,Nlin3,t+tMS,&val);
	  //	  fprintf(stderr," Are you in SSE3?\n");
	  if (isinf(val) || val!=val) {
	    fprintf(stderr," got val=%le at t=%le tMS=%le\n",val,t,tMS);
	  }
	  else {
	    m3=val;
	  
	  }
	  v[4]=m3;
	}
	else {
	  m3=v[4];
	}

	
	m12 = m1 + m2;
	m123 = m1 + m2 + m3;
	//smadar: in case of tidal frictions these parameters vary
	
	a1 = SQR(htot)/ (k2 * (m1+m2) *(1-SQR(e1)));

	a2 = SQR(G2 / ( (m1+m2)*m3 ) ) *(m1+m2+m3) / (k2 *(1-SQR(e2)));
	/*
  if (t+tMS> 309.4e6){
    fprintf(stderr," m3=%le t=%le tMS=%le a2=%le\n",m3,t,tMS,a2);
    }*/
	/*
	if (t>309.4*1e6+tMS){
	  fprintf(stderr,"t=%le m3=%le a2=%le\n",t,m3,a2);
	  }*/
	/*
	  fprintf(stderr,"in octupole \n");
	  fprintf(stderr,"t=%le Hsq=%le\n",t,Hsq);
	  fprintf(stderr,"G1=%le G2=%le cosi=%le\n",G1,G2,x[7]);
	  fprintf(stderr,"m1=%le m2=%le m3=%le e2=%le a2=%le a1=%le e1=%le\n",m1,m2,m3,e2,a2,a1,e1);
	*/
	  //	  exit(7);
	// now that we have the new paramters, a1, a2, mass etc we can re-calcualte the relavent paramters below:

	P1 = 2*M_PI*sqrt(a1*a1*a1/(k2*(m1+m2)));
	P2 = 2*M_PI*sqrt(a2*a2*a2/(k2*(m1+m2+m3)));
		
	//	L1 = m1 * m2  / m12  * sqrt( k2 * m12  * a1 );
	//	L2 = m3 * m12 / m123 * sqrt( k2 * m123 * a2 );
		/* Compute convenient variables */
	//	G1 = L1 * sqrt(1 - e1*e1);
	//	G2 = L2 * sqrt(1 - e2*e2);

	L1 = G1 / sqrt(1 - e1*e1);
	L2= G2 / sqrt(1 - e2*e2);

	C3G23 = k2*k2/16.0 * pow(m12*m3,7)/pow(m123*m1*m2,3) * pow(L1/L2,3) * L1;
	C4G25 = -k2*k2*(15/64.0) * pow(m12*m3,9)*(m1-m2)/pow(m123,4)/pow(m1*m2,5) * pow(L1*L1/L2,3);


	
	n = 2.*M_PI/P1;
	//paramters for the spins
	In1=alpha1*m1*R1*R1;
	In2=alpha2*m2*R2*R2;
	mu=m1*m2/(m1+m2);


	cosg1 = e1cosg1/e1;
	sing1 = e1sing1/e1;
	cosg2 = e2cosg2/e2;
	sing2 = e2sing2/e2;
	cos2g1 = cosg1*cosg1-sing1*sing1;
	sin2g1 = 2*cosg1*sing1;

	/* debug */
	DebugPrint2("cos(g1)=%g sin(g1)=%g ",cosg1,sing1);
	DebugPrint2("cos(g2)=%g sin(g2)=%g ",cosg2,sing2);
	DebugPrint2("cos(2g1)=%g sin(2g1)=%g  \n",cos2g1,sin2g1);


	/*
		  fprintf(stderr,"in octupole \n");
	  fprintf(stderr,"t=%le Hsq=%le\n",t, G1*G1 + 2.*G1*G2*x[7] + G2*G2);
	  fprintf(stderr,"G1=%le G2=%le cosi=%le\n",G1,G2,x[7]);
	  fprintf(stderr,"m1=%le m2=%le m3=%le e2=%le a2=%le a1=%le e1=%le\n",m1,m2,m3,e2,a2,a1,e1);
	  //  exit(7);
	  */
	
	//	if(t==0)    Hsq = G1*G1 + 2.*G1*G2*cos(i) + G2*G2;
		  Hsq = G1*G1 + 2.*G1*G2*x[7] + G2*G2;
	if(t==0){
	  Hsq = G1*G1 + 2.*G1*G2*x[7] + G2*G2;
	  fprintf(stderr,"in octupole t=0\n");
	  fprintf(stderr,"t=%le Hsq=%le\n",t,Hsq);
	  fprintf(stderr,"G1=%le G2=%le cosi=%le htot=%le\n",G1,G2,x[7],htot);
	  fprintf(stderr,"m1=%le m2=%le m3=%le e2=%le a2=%le a1=%le e1=%le\n",m1,m2,m3,e2,a2,a1,e1);
	}
	/*if (t>1e9){
  fprintf(stderr,"t=%le Hsq=%le\n",t,Hsq);
  }*/
	/*
	if (t<1009917.6 && t>1009916){
	  //fprintf(stderr,"t=%le Hsq=%le\n",t,Hsq);
	  }*/
	//	exit(3);

       	theta = (Hsq - G1*G1 - G2*G2)/(2.0*G1*G2);

	//fprintf(stderr," theta=%le x[7]=%le\n",theta,x[7]);

	if(fabs(theta)<=1)    v[1]=acos(theta);
	else if(theta>1)      v[1] = 0.0;
	else if(theta<1)      v[1] = M_PI;

	theta=x[7];
	//finding the inclinations angles

	if(fabs(theta)<=1)    i=acos(theta);
	else if(theta>1)      i = 0.0;
	else if(theta<1)      i = M_PI;
	sini=sin(i);

	if(fabs(x[5])<=1)    	i1 = acos(x[5]);
	else if(x[5]>1)      i1 = 0.0;
	else if(x[5]<1)       i1 = M_PI;
	theta1=cos(i1);
	H1= G1*theta1;

	sini1=sin(i1);
	if(fabs(x[6])<=1)    	i2 = acos(x[6]);
	else if(x[6]>1)      i2 = 0.0;
	else if(x[6]<1)       i2 = M_PI;

	//fprintf(stderr,"t=%le H1=%le H2=%le val=(%le,%le)\n",H1,G2*cos(i2),sqrt(1-e1*e1)*cos(i1),sqrt(1-e2*e2)*cos(i2),t);

	B = 2 + 5*e1*e1 - 7*e1*e1*cos2g1;
	A = 4 + 3*e1*e1 - 2.5*(1-theta*theta)*B;
	C3 = C3G23 / pow(G2,3);
	C4 = C4G25 / pow(G2,5);

	cosphi = - cosg1 * cosg2 - theta * sing1 * sing2;

	/* debug */
	DebugPrint2("G1=%g G2=%g  "		,G1,G2); 
	DebugPrint2("Hsq=%g theta=%g  "	,Hsq,theta);
	DebugPrint2("B=%g A=%g "		,B,A);
	DebugPrint2("C3=%g C4=%g "		,C3,C4);
	DebugPrint1("cosphi=%g  \n"		,cosphi);
	
	dg1dt = 0;
	de1dt = 0;
	dg2dt = 0;
	de2dt = 0;

	dhtotdt =0;
	dspin1edt = dspin1qdt  = dspin1hdt = dspin2edt = dspin2qdt  = dspin2hdt =0; 
	didt=di1dt=di2dt=0;
	di1dtTF=dg1dtTF=dOmega1dtTF=dOmega1dt=0;
	dFdg1=dFdg2=0.;
	dBigG1dt=dBigG2dt=0;


	//#ifdef QUADRUPOLE                  /* quadrupole approximation */
	if(!(QUADRUPOLE)){	
		dg1dt += C3 * 6 * ((1/G1)*(4*theta*theta+(5*cos2g1-1)*(1-e1*e1-theta*theta) ) \
       		+ (theta/G2)*(2+e1*e1*(3-5*cos2g1)));
	


		de1dt += C3 * ((1-e1*e1)/G1) * (30*e1*(1-theta*theta)*sin2g1);
		
		dg2dt += C3 * 3 * ((2*theta/G1)*(2+e1*e1*(3-5*cos2g1) ) + \
				(1/G2)*(4+6*e1*e1 + (5*theta*theta-3)*(2+e1*e1*(3-5*cos2g1))));
		
		de2dt += 0.;

		dFdg1 += C3*15. * SQR(e1)* (1 - SQR(theta))* ( -2.*sin2g1) ; 

		dFdg2 += 0.;
	
		dOmega1dt += -C3*(6.*theta*(2.+3.*e1*e1)-30.*theta*e1*e1*cos2g1)*sini/(G1*sini1);

}
	
	  //#ifdef OCTUPOLE                   /* octupole approximation */
	if(!(OCTUPOLE)){

	dg1dt += -C4 * e2 * (e1*(1/G2 + theta/G1) * (sing1*sing2 \
			* (A+10*(3*theta*theta-1)*(1-e1*e1)) - 5*theta*B*cosphi ) \
			- (1-e1*e1)/(e1*G1) * ( sing1*sing2*10*theta*(1-theta*theta)\
			* (1-3*e1*e1) + cosphi * (3*A - 10*theta*theta+2)));

	de1dt += C4 * e2 * (1-e1*e1) / G1 * (35*cosphi*(1-theta*theta)*e1*e1*sin2g1\
			- 10*theta*(1-e1*e1)*(1-theta*theta)*cosg1*sing2 \
			- A*(sing1*cosg2-theta*cosg1*sing2));

	dg2dt += C4 * e1 * (sing1*sing2*((4*e2*e2+1)/(e2*G2)*10*theta*(1-theta*theta)\
			* (1 - e1*e1) - e2*(1/G1+theta/G2) * (A + 10*(3*theta*theta-1)\
			* (1 - e1*e1))) + cosphi*(5*B*theta*e2* (1/G1 + theta/G2) \
			+ (4*e2*e2 + 1)/(e2*G2)*A));

	de2dt += -C4 * e1 * (1-e2*e2)/G2 * (10*theta*(1-theta*theta)*(1-e1*e1)\
			* sing1*cosg2 + A*(cosg1*sing2-theta*sing1*cosg2));


	dFdg1 += C4 * e1 * e2 *( -35. *SQR(e1) *(1.-SQR(theta) )*sin2g1 * cosphi + A * (sing1 * cosg2 - theta* cosg1 * sing2) + 10.*theta *(1. -SQR(theta)) * (1.-SQR(e1))*cosg1 * sing2 );  

	dFdg2 += C4*e1*e2 *( A* (cosg1 * sing2 - theta * sing1 * cosg2 ) +  10.*theta *(1. -SQR(theta)) * (1.-SQR(e1))* sing1 * cosg2 );


	dOmega1dt += -C4*e1*e2*(5*theta*B*cosphi - A*sing1*sing2 + (10.-30.*theta*theta)*(1-e1*e1)*sing1*sing2)*sini/(G1*sini1);
	//#endif
	}

//#ifdef GR						/* include general reletavistic precession */
	if(!(GR)){
	dg1dt += 6*M_PI*k2*m12/(P1*a1*(1.0-e1*e1)*c2);
	dg2dt += 6*M_PI*k2*m123/(P2*a2*(1.0-e2*e2)*c2);
	}

	//#endif





	/* smadar: general derivatives:  a result only from non-disipative effects! */

	dhtotdt += dFdg1 * (m1+m2)/(m1*m2) ;

	dBigG1dt += dFdg1;
	dBigG2dt += dFdg2;

	//	dhtotdt += dBigG1dt * (m1+m2)/(m1*m2) ;

	//note: This is actualy dcos(i)dt!!!
	//Note that these are actualy dcosidt
	/*	di1dt=-dBigG1dt * (Hsq-G2*G2+G1*G1)/(2.*sqrt(Hsq)*G1*G1) +(-2.*G2*dBigG2dt+2.*G1*dBigG1dt )/(2.*sqrt(Hsq)*G1);

	di2dt=-dBigG2dt * (Hsq+G2*G2-G1*G1)/(2.*sqrt(Hsq)*G2*G2) +(+2.*G2*dBigG2dt-2.*G1*dBigG1dt )/(2.*sqrt(Hsq)*G2);
	*/

	//	di1dt=BigH1dot/G1-dBigG1dt/G1 *cos(i1);
	//	di2dt=BigH2dot/G2-dBigG2dt/G2 *cos(i2);


	if(!(ML1)){

	  /// mass loss happen with htot constant. Here we calcl the derivative only due to mass loss.
	  dm1dt= dmdt(t, m1, v[2]);
	 
	  

	  //	  BigH1dotML=BigHtotdotML + G1/ sqrt(Hsq) * dBigG1dtML  - G2 /sqrt(Hsq) * dBigG2dtML - BigHtotdotML / sqrt(Hsq) * G1 *cos(i1);

	  //	  BigH2dotML=BigHtotdotML + G2/ sqrt(Hsq) * dBigG2dtML  - G1 /sqrt(Hsq) * dBigG1dtML - BigHtotdotML / sqrt(Hsq) * G2 *cos(i2);

	  // note this is actualy d cos(i1) dt, also - 
	  //	  di1dtML=BigH1dotML/G1-dBigG1dtML/G1 *cos(i1);
	  //	  di2dtML=BigH2dotML/G2-dBigG2dtML/G2 *cos(i2);
	  
	  //	  dR1dt=0;
	  	  dR1dt=dRdt(t,R1);
	}

	if(!(ML2)){

	  dm2dt= dmdt2(t,m2,v[3]);
	  dR2dt= dRdt2(t,R2);
	}


	dBigG1dtML += dm1dt* m2/ m1* G1 /(m1+m2)  + dm2dt* m1/ m2* G1/(m1+m2);
	  //	  dhtotdt +=   dm1dt* m2/ m1* G1 / (m1*m2);

	
	dBigG2dtML += (dm1dt+dm2dt) *m3 /( (m1+m2)*(m1+m2+m3)) * G2 ;

	  
	BigHtotdotML=dBigG1dtML*cos(i1)+dBigG2dtML*cos(i2);

	//		fprintf(stderr,"dm1dt=%le\n",dm1dt);

	//	BigH1dot=G1/sqrt(Hsq) *dBigG1dt-G2/ sqrt(Hsq) *dBigG2dt;

	//BigH1dot=sin(i2)/sin(i)*dBigG1dt-sin(i1)/sin(i)*dBigG2dt;
	//	BigH2dot=-BigH1dot;

	dBigG1dt += dBigG1dtML;
	dBigG2dt += dBigG2dtML;



	BigHtotdot=BigHtotdotML;

	BigH1dot=BigHtotdot  - BigHtotdot / sqrt(Hsq) * G1 *cos(i1) + G1/ sqrt(Hsq) * dBigG1dt  - G2 /sqrt(Hsq) * dBigG2dt;

	BigH2dot=BigHtotdot  - BigHtotdot / sqrt(Hsq) * G2 *cos(i2) + G2/ sqrt(Hsq) * dBigG2dt  - G1 /sqrt(Hsq) * dBigG1dt;

	  //	BigH1dot +=BigH1dotML;
	  //	BigH2dot +=BigH2dotML;

	di1dt=BigH1dot/G1-dBigG1dt/G1 *cos(i1);

	di2dt=BigH2dot/G2-dBigG2dt/G2 *cos(i2);

	//if (!(MB1)){
	//  if (!(SSE1)){
	//    if (t+tMS<timeSSE[2]) {
	//       menv=MenvSSE[1];
	//     }
	//    else {
	//      linint(timeSSE,MenvSSE,Nlin,t+tMS,&menv);
	//     }
	//  }
	//  if (menv<=1.e-10) menv=0.;

	//  dJspindtMB =  -MBfactor* pow(R1*spin1,3)*SQR(Rsun) * menv/m1;
	//  dspin1edt += dJspindtMB/In1 ; //               +Ztot * spin1q - Ytot * spin1h;
	//  dspin1qdt += dJspindtMB/In1 ;// -Ztot * spin1e                 + Xtot * spin1h ;
	//  dspin1hdt += dJspindtMB/In1 ;//+Ytot * spin1e -Xtot * spin1q;
	//}

	//if (!(MB2)){
	//  if (!(SSE2)){
	//  	if (t+tMS<timeSSE2[2]) {
	//       menv2=MenvSSE2[1];
	//     }
	//    else {
	//      linint(timeSSE2,MenvSSE2,Nlin2,t+tMS,&menv2);
	//     }
	//  }
	//  if (menv2<=1.e-10) menv2=0.;

	//  dJspin2dtMB =  -MBfactor* pow(R2*spin2,3)*SQR(Rsun) * menv2/m2;
	//  dspin2edt += dJspin2dtMB/In2 ; //               +Ztot * spin1q - Ytot * spin1h;
	//  dspin2qdt += dJspin2dtMB/In2 ;// -Ztot * spin1e                 + Xtot * spin1h ;
	//  dspin2hdt += dJspin2dtMB/In2 ;//+Ytot * spin1e -Xtot * spin1q;
	//}


	//	di1dt +=di1dtML;
	//	di2dt +=di2dtML;

	//	dhtotdt += dBigG1dt * (m1+m2)/(m1*m2) ;

	if(!(TF)){
	  //	fprintf(stderr,"TF=%d\n",TF);
	  //	exit(1);
	//This is from Dan's paper
	  
	//tF1=TV1  *pow(a1/R1,8)*SQR(m1)/ ( 9.*(m1+m2)*m2*SQR(1+2.*kL1));
	//tF2=TV2  *pow(a1/R2,8)*SQR(m2)/ ( 9.*(m1+m2)*m1*SQR(1+2.*kL2));
	
	if((m1 < 1.4)||(type1!=1)||(type1!=11)){ 
	tF1=TV1  *pow(a1/R1,8)*SQR(m1)/ ( 9.*(m1+m2)*m2*SQR(1+2.*kL1)); //convective stars
	}
	else if((m1 > 1.4)&&(type1==1)||(type1==11)){ 
	tF1=pow(a1/R1,9)*pow((pow(a1,3)/(k2*m1)),0.5)*m1/ (9.*1.592e-9*pow(m1,2.84)*m2*pow(1.+(m2/m1),11./6.)); //radiative stars
	//tF1= (2./21.) * pow(a1/R1,9) * pow(SQR(m1)/ ( (m1+m2)*m2 ), 6./11.) * sqrt(pow(a1,3)/(k2*m1)) * 6.2814e8 * pow(m1, -2.84); //radiative stars - test
	}

	if((m2 < 1.4)||(type2!=1)||(type2!=11)){ 
	tF2=TV2  *pow(a1/R2,8)*SQR(m2)/ ( 9.*(m1+m2)*m1*SQR(1+2.*kL2)); //convective stars
	}
	else if((m2 > 1.4)&&(type2==1)||(type2==11)){ 
	tF2=pow(a1/R2,9)*pow((pow(a1,3)/(k2*m2)),0.5)*m2/ (9.*1.592e-9*pow(m2,2.84)*m1*pow(1.+(m1/m2),11./6.)); //radiative stars
	//tF1= (2./21.) * pow(a1/R1,9) * pow(SQR(m1)/ ( (m1+m2)*m2 ), 6./11.) * sqrt(pow(a1,3)/(k2*m1)) * 6.2814e8 * pow(m1, -2.84); //radiative stars - test
	}

	//This is from Dan's code
	/*
	tF1=TV1  *pow(a1/R1,8)*SQR(m1)* SQR(1.-Q1)/ (9* (m1+m2)*m2);
      	tF2=TV2  *pow(a1/R2,8)*SQR(m2) *SQR(1.-Q2) / (9* (m1+m2)*m1);
	*/

       	e1_2=e1*e1;
	e1_4=e1*e1*e1*e1;
	e1_6=e1_2*e1_2*e1_2;
	e1_8=e1_4*e1_4;

	//TF coefficents
	V1= 9.*( (1. +15./4. *e1_2 +15./ 8. *e1_4+ 5./64. *e1_6)/pow(1.-SQR(e1),13./2.) - 11.* spin1h / (18. * n )*( 1. + 3./2. * e1_2 + e1_4/8.)/pow(1-e1_2,5) ) / tF1  ;

       	V2= 9.*( (1. +15./4. *e1_2 +15./ 8. *e1_4+ 5./64. *e1_6)/pow(1.-SQR(e1),13./2.) - 11.* spin2h / (18. * n )*( 1. + 3./2.* e1_2 + e1_4/8.)/pow(1-e1_2,5) ) / tF2  ;
	//V2=0;

	W1= ( (1. +15./2. *e1_2 +45./ 8. *e1_4+ 5./16. *e1_6)/pow(1.-SQR(e1),13./2.) - spin1h / ( n )*( 1. + 3.* e1_2 + 3. *e1_4/8.)/pow(1-e1_2,5) ) / tF1  ;

	W2= ( (1. +15./2. *e1_2 +45./ 8. *e1_4+ 5./16. *e1_6)/pow(1.-SQR(e1),13./2.) - spin2h / ( n )*( 1. + 3.* e1_2 + 3. *e1_4/8.)/pow(1-e1_2,5) ) / tF2  ;
	//	W2=0.;

	X1= -pow(R1/a1,5) *kL1*m2/(mu *n ) *spin1h*spin1e/(SQR(1-e1_2)) - spin1q / (2. *n *tF1 ) *( 1. +9. /2. * e1_2 + 5./ 8. *e1_4 )/pow(1.-e1_2,5);

	Y1= -pow(R1/a1,5) *kL1*m2/(mu *n ) *spin1h*spin1q/(SQR(1-e1_2)) + spin1e / (2. *n *tF1 ) *( 1. +3. /2. * e1_2 + 1./ 8. *e1_4 )/pow(1.-e1_2,5);

	Z1= pow(R1/a1,5) *kL1*m2/(mu *n ) * ( (2* SQR(spin1h)- SQR(spin1q)- SQR(spin1q) )/ (2.* SQR(1-e1_2) )+ 15. * k2 *m2 / pow(a1,3)*( 1. +3. /2. * e1_2 + 1./ 8. *e1_4 )/pow(1.-e1_2,5));


	X2= -pow(R2/a1,5) *kL2*m1/(mu *n ) *spin2h*spin2e/(SQR(1-e1_2)) - spin2q / (2. *n *tF2 ) *( 1. +9. /2. * e1_2 + 5./ 8. *e1_4 )/pow(1.-e1_2,5);
	//X2=0.;

	Y2= -pow(R2/a1,5) *kL2*m1/(mu *n ) *spin2h*spin2q/(SQR(1-e1_2)) + spin2e / (2. *n *tF2 ) *( 1. +3. /2. * e1_2 + 1./ 8. *e1_4 )/pow(1.-e1_2,5);
	//Y2=0.;			     
	
	Z2= pow(R2/a1,5) *kL2*m1/(mu *n ) * ( (2* SQR(spin2h)- SQR(spin2q)- SQR(spin2q) )/ (2.* SQR(1-e1_2))+ 15. * k2 *m1 / pow(a1,3)*( 1. +3. /2. * e1_2 + 1./ 8. *e1_4 )/pow(1.-e1_2,5));
	//Z2=0.;

	//	fprintf(stderr,"X1=%f X2=%f Y1=%f Y2=%f Z1=%f Z2=%f V1=%f V2=%f  W1=%le W2=%le\n",X1,X2,Y1,Y2,Z1,Z2,V1,V2,W1,W2);

	//	ZGR = 3.*k2*(m1+m2)* n/(a1 * SQR(1-e1*e1) *c2); 
	
	detotdt= ( - (V1+ V2) )*(e1);
	dhtotdt +=(  -(W1 +W2)  ) * (htot);

	detotdtTF[1] = detotdt;
	dhtotdtTF[1] = dhtotdt;
	
	if (fabs(sini1)<1.e-12){
	  dOmega1dtTF = 0.;
	}
	else{
	  dOmega1dtTF = ( (X1+X2)*sing1 + (Y1+Y2)*cosg1 )/(sini1);
	}

	di1dtTF = (X1+X2)*cosg1 - (Y1+Y2)*sing1 ;


	
	dg1dtTF = Z1+Z2 - dOmega1dtTF *theta1;
	



	//#endif
	}

    

	//addeing the TF if flaged:


	

	/* 
	   From non-desepetive only secular dynamics we clauclated dcos(i)dt, the TF needed to be added to the 
	   actual didt. Note dcos(i)dt=-sin(i)didt
	*/

	if (fabs(sin(i1))<1.e-12){    	//dcos(i1)dt = -sin(i1)di1dt
	  di1_onlydt = 0;
	}
	else {
	  di1_onlydt= -di1dt/sin(i1);
	}



	de1dt += detotdt;
	dg1dt += dg1dtTF;
	/* now adding the TF fraction which is actualy di1dt and NOT dcos(i1)dt */

	di1_onlydt = di1_onlydt + di1dtTF;
	
	/* now transform back to dcos(i1)dt */
	di1dt = -di1_onlydt * sin(i1);



	dOmega1dt += dOmega1dtTF;

	//calcualting d(cos(i)/ dt

	
	if (fabs(sin(i1))>1.e-12 && fabs(sin(i2))>1.e-12){
	  didt = sin(i)* (di1dt / sin(i1) + di2dt/sin(i2));
	}
	else if (fabs(sin(i1))>1.e-12 && fabs(sin(i2))<1.e-12){
	  didt = sin(i)* (di1dt / sin(i1) );
	}
	else if (fabs(sin(i1))<1.e-12 && fabs(sin(i2))>1.e-12){
	  didt = sin(i)* (di2dt/sin(i2));
	}
	else {
	  didt = 0;
	}

	//di1dtTF=-sin(i1)*di1dtTF;  //since we are clacualting actualy dcos(i)dt=-sin(i)didt
	
	Ztot = dg1dt + dOmega1dt*theta1;
	Ytot = -di1_onlydt * sing1 + dOmega1dt * cosg1 * sini1;
	Xtot = di1_onlydt * cosg1 +  dOmega1dt * sing1 * sini1;

	dspin1edt +=-Y1*htot*mu/In1                +Ztot * spin1q - Ytot * spin1h;
	dspin1qdt +=X1*htot*mu/In1 -Ztot * spin1e                 + Xtot * spin1h ;
	dspin1hdt +=W1*htot*mu/In1+Ytot * spin1e -Xtot * spin1q;
	
	
	dspin2edt +=-Y2*htot*mu/In2                +Ztot * spin2q - Ytot * spin2h;
	dspin2qdt +=X2*htot*mu/In2 -Ztot * spin2e                 +Xtot * spin2h ;
	dspin2hdt +=W2*htot*mu/In2+Ytot * spin2e -Xtot * spin2q;

	if (!(MB1)){
	  if ((SSE1)==1){
	    linint(timeSSE,MenvSSE,Nlin,t,&menv);
	     
	    if (menv<=1.e-10) menv=0.;

	    if (m1<=1.2) dJspindtMB =  -MBfactor* pow(R1*spin1/Rsun,3)*SQR(Rsun) * menv/m1;
	    else dJspindtMB =  -0.1*MBfactor* pow(R1*spin1/Rsun,3)*SQR(Rsun) * menv/m1;

	    dspin1edt += spin1e * dJspindtMB/(In1*spin1);//              +Ztot *spin1q - Ytot * spin1h;
	    dspin1qdt += spin1q * dJspindtMB/(In1*spin1);//  -Ztot * spin1e      + Xtot * spin1h ;
	    dspin1hdt += spin1h * dJspindtMB/(In1*spin1);// +Ytot * spin1e -Xtot * spin1q;
	  }
	  else if (!(SSE1)){
	    dspintot1dt = dspindt(t,spin1);
	    dspin1edt += spin1e * dspintot1dt/(spin1);//              +Ztot *spin1q - Ytot * spin1h;
	    dspin1qdt += spin1q * dspintot1dt/(spin1);//  -Ztot * spin1e      + Xtot * spin1h ;
	    dspin1hdt += spin1h * dspintot1dt/(spin1);// +Ytot * spin1e -Xtot * spin1q;
	  }

	}


	if (!(MB2)){
	  if ((SSE2)==1){
	    linint(timeSSE2,MenvSSE2,Nlin2,t+tMS,&menv2);
	     
	    if (menv2<=1.e-10) menv2=0.;
	    
	    if (m1<=1.2) dJspindtMB =  -MBfactor* pow(R2*spin2/Rsun,3)*SQR(Rsun) * menv2/m2;
	    else dJspindtMB =  -0.1*MBfactor* pow(R2*spin2/Rsun,3)*SQR(Rsun) * menv2/m2;

	    dspin2edt += spin2e * dJspindtMB/(In2*spin2);//              +Ztot *spin2q - Ytot * spin2h;
	    dspin2qdt += spin2q * dJspindtMB/(In2*spin2);//  -Ztot * spin2e       + Xtot * spin2h ;
	    dspin2hdt += spin2h * dJspindtMB/(In2*spin2);// +Ytot * spin2e -Xtot * spin2q;
	  }
	  else if (!(SSE2)){
	    dspintot2dt = dspindt2(t,spin2);
	    dspin2edt += spin2e * dspintot2dt/(spin2);//              +Ztot *spin2q - Ytot * spin2h;
	    dspin2qdt += spin2q * dspintot2dt/(spin2);//  -Ztot * spin2e       + Xtot * spin2h ;
	    dspin2hdt += spin2h * dspintot2dt/(spin2);// +Ytot * spin2e -Xtot * spin2q;
	  }

	}

	/*
	dspin2edt=0.;
	dspin2qdt=0.;
	dspin2hdt=0.;
	*/
	//		fprintf(stderr," == t=%le i=%le i1=%le i2=%le m1=%le e1=%le e2=%le \n",t,i,i1,i2,m1,e1,e2);
		//	fprintf(stderr,"    didt=%le di1dt=%le di2dt=%le \n",didt,di1dt,di2dt);


	//	fprintf(stderr," t=%le d1spin(e,q,h)=(%le %le %le) d1spin(e,q,h)=(%le %le %le)\n",t,dspin1edt,dspin1qdt,dspin1hdt,dspin2edt,dspin2qdt,dspin2hdt);
	
	dxdt[1] =  e1cosg1 * dg1dt + sing1*de1dt;
	dxdt[2] = -e1sing1 * dg1dt + cosg1*de1dt;
	dxdt[3] =  e2cosg2 * dg2dt + sing2*de2dt;
	dxdt[4] = -e2sing2 * dg2dt + cosg2*de2dt;
	
	dxdt[5] = di1dt;
	dxdt[6] = di2dt;
	dxdt[7] = didt; 

	dxdt[8] = dhtotdt;

	//	dspin1edt=dspin1qdt=dspin1hdt=dspin2edt=dspin2qdt=dspin2hdt=0.;

	dxdt[9] = dspin1edt;
	dxdt[10] = dspin1qdt;
	dxdt[11] = dspin1hdt;
	
	dxdt[12] = dspin2edt;
	dxdt[13] = dspin2qdt;
	dxdt[14] = dspin2hdt;

	dxdt[15]=0.;//dR1dt;
	dxdt[16]=dm1dt;

	dxdt[17]=dBigG1dt;

	dxdt[18]=dBigG2dt;

	dxdt[19]=0.;//dR2dt;
	dxdt[20]=dm2dt;
	/*
	if (t>3577360000){
	fprintf(stderr," t=%le d1spin(e,q,h)=(%le %le %le) d1spin(e,q,h)=(%le %le %le)\n",t,dspin1edt,dspin1qdt,dspin1hdt,dspin2edt,dspin2qdt,dspin2hdt);
	fprintf(stderr," == t=%le i=%le i1=%le i2=%le m1=%le m2=%le m3=%le e1=%le e2=%le \n",t,i,i1*180/M_PI,i2*180/M_PI,m1,m2,m3,e1,e2);
		fprintf(stderr,"    didt=%le di1dt=%le di2dt=%le \n",didt,di1dt,di2dt);
		fprintf(stderr, "  R1=%le dR1dt=%le\n",R1,dR1dt);
		//	exit(1);	
		}*/
	
	//dxdt[15] = dOmega1dt;	
	





	for (j=1;j<=n_orb_param;j++){
	  ve[j]=v[j];
	}



	for (j=1;j<=n_eq;j++){

	  // if (fabs(dxdt[j])<SmallNumber) dxdt[j]=0.;

	  if (dxdt[j]!=dxdt[j]){
	    fprintf(stderr,"t=%le dxdt[%d]=%le\n",t,j,dxdt[j]);
	    flgout=1;
	   
	    // exit(1);
	  }
	}
	int m;

	v[1] = i;
	v[4]=m3;
	ve[4]=m3;
	/*	fprintf(stderr,"   t=%le e2=%le\n",t,e2);
	fprintf(stderr," de2dt = %le dg2dt=%le \n",de2dt,dg2dt);
       	fprintf(stderr," dxdt[3]=%g dxdt[4]=%g\n",dxdt[3],dxdt[4]);
	*/
	/* debug */
	DebugPrint2("dg1dt=%g de1dt=%g ",dxdt[1],dxdt[2]);
	DebugPrint2("dg2dt=%g de2dt%g\n",dxdt[3],dxdt[4]);

}
