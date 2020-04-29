#include <math.h>
#include "HTI.h"
#include <stdio.h>
#include "nrutil.h"


extern double *v,*ve;
extern int flgout;

/* ODE for a damped harmonic oscillator to be integrated */
void octupole (double t, double *x, double *dxdt) {

  static double m1,m2,m3,a1,a2,P1,P2,R1,R2;
	static double L1,L2, C3G23, C4G25, m12, m123, Hsq;
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
	double Omega1,Omega2;
	double dOmega1dt,dOmega2dt,dOmega1dtTF,htot,dhtotdt,dg1dtTF;
	double n;
	//spins:
	double spin1h,spin2h,spin1q,spin2q,spin1e,spin2e,spin1,spin2,spin1tot,spin2tot,rootval; //spin parameters
	double dspin1edt,dspin1qdt,dspin1hdt,dspin2edt,dspin2qdt,dspin2hdt;
	double X1,X2,Y1,Y2,V1,V2,W1,W2,Z1,Z2,ZGR; //TF cooeficents
	double e1_2,e1_4,e1_6,e1_8; //for TF 
	double sini1,sini;
	dspin1edt=dspin1qdt=dspin1hdt=dspin2edt=dspin2qdt=dspin2hdt=0.;
	X1=X2=Y1=Y2=Z1=Z2=W1=W2=V1=V2=0;

	double kL1=0.5*Q1/(1.-Q1);
	
	double kL2=0.5*Q2/(1.-Q2);

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

	
	//smadar: in case of tidal frictions these parameters vary
	
	a1= SQR(htot)/ (k2 * (m1+m2) *(1-SQR(e1)));

	P1 = 2*M_PI*sqrt(a1*a1*a1/(k2*(m1+m2)));
	P2 = 2*M_PI*sqrt(a2*a2*a2/(k2*(m1+m2+m3)));
		
	L1 = m1 * m2  / m12  * sqrt( k2 * m12  * a1 );
	L2 = m3 * m12 / m123 * sqrt( k2 * m123 * a2 );
	
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

	/* Compute convenient variables */
	G1 = L1 * sqrt(1 - e1*e1);
	G2 = L2 * sqrt(1 - e2*e2);


	//	if(t==0)    Hsq = G1*G1 + 2.*G1*G2*cos(i) + G2*G2;
	
	if(t==0){
	  Hsq = G1*G1 + 2.*G1*G2*x[7] + G2*G2;
	  fprintf(stderr,"t=%le Hsq=%le\n",t,Hsq);
	}
	/*if (t>1e9){
  fprintf(stderr,"t=%le Hsq=%le\n",t,Hsq);
  }*/
	/*
	if (t<1009917.6 && t>1009916){
	  //fprintf(stderr,"t=%le Hsq=%le\n",t,Hsq);
	  }*/
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
		/*
		dBigG1dt += dFdg1;
		dBigG2dt += 0.;
		*/

	
		//dOmega1dt += C3 * 4. * H1 * ( 2.+3.*e1*e1) *(3.*SQR(theta1)-1) + 15.*SQR(e1)*(1-SQR(theta1) *cos2g1);
		//	dOmega1dt += -C3*(6.*theta*(2.+3.*e1*e1)-30.*theta*e1*e1*cos2g1)*sini/(G1*sini1);
		//		dOmega1dt += -C3*(6.*theta*(2.+3.*e1*e1)-30.*theta*e1*e1*cos2g1)*sqrt(Hsq)/(G1*G2);
		dOmega1dt += -C3*(6.*theta*(2.+3.*e1*e1)-30.*theta*e1*e1*cos2g1)*sini/(G1*sini1);

}
		//#endif
	//	fprintf(stderr,"4: t=%le de1dt=%le dg1dt=%le\n", t,de1dt,dg1dt);
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

	/*
	dBigG1dt += dFdg1;
	dBigG2dt += dFdg2;
	*/

	//dOmega1dt += C4 * 6. * H1 * e1*e2* (A * cosphi +10.*theta1* (1-SQR(theta1))*(1-SQR(e1))*sing1*sing2 );
	//	dOmega1dt += -C4*e1*e2*(5*theta*B*cosphi - A*sing1*sing2 + 10.+30.*theta*theta)*sini/(G1*sini1);
	//	dOmega1dt += -C4*e1*e2*(5*theta*B*cosphi - A*sing1*sing2 + 10.+30.*theta*theta)*sini/(G1*sini1);
	//	dOmega1dt += -C4*e1*e2*(5*theta*B*cosphi - A*sing1*sing2 + (10.-30.*theta*theta)*SQR(1-e1*e1)*sing1*sing2)*sini/(G1*sini1);
	//dOmega1dt += -C4*e1*e2*(5*theta*B*cosphi - A*sing1*sing2 + (10.-30.*theta*theta)*(1-e1*e1)*sing1*sing2)*sqrt(Hsq)/(G1*G2);
	dOmega1dt += -C4*e1*e2*(5*theta*B*cosphi - A*sing1*sing2 + (10.-30.*theta*theta)*(1-e1*e1)*sing1*sing2)*sini/(G1*sini1);
	//#endif
	}
	/*
	fprintf(stderr,"8: t=%le de1dt=%le dg1dt=%le\n", t, C4 * e2 * (1-e1*e1) / G1 * (35*cosphi*(1-theta*theta)*e1*e1*sin2g1 \
			- 10*theta*(1-e1*e1)*(1-theta*theta)*cosg1*sing2 \
											- A*(sing1*cosg2-theta*cosg1*sing2)), C4 * e2 * (1-e1*e1) / G1 * (35*cosphi*(1-theta*theta)*e1*e1*sin2g1\
			- 10*theta*(1-e1*e1)*(1-theta*theta)*cosg1*sing2 \
																			  - A*(sing1*cosg2-theta*cosg1*sing2)));
	*/	
//#ifdef GR						/* include general reletavistic precession */
	if(!(GR)){
	dg1dt += 6*M_PI*k2*m12/(P1*a1*(1.0-e1*e1)*c2);
	dg2dt += 6*M_PI*k2*m123/(P2*a2*(1.0-e2*e2)*c2);
	}

	//#endif



	/* smadar: general derivatives:  a result only from non-disipative effects! */

	dhtotdt += dFdg1 * (m1+m2)/(m1*m2) ;

	//note: This is actualy dcos(i)dt!!!


	dBigG1dt += dFdg1;
	dBigG2dt += dFdg2;
	//Note that these are actualy dcosidt
	di1dt=-dBigG1dt * (Hsq-G2*G2+G1*G1)/(2.*sqrt(Hsq)*G1*G1) +(-2.*G2*dBigG2dt+2.*G1*dBigG1dt )/(2.*sqrt(Hsq)*G1);

	di2dt=-dBigG2dt * (Hsq+G2*G2-G1*G1)/(2.*sqrt(Hsq)*G2*G2) +(+2.*G2*dBigG2dt-2.*G1*dBigG1dt )/(2.*sqrt(Hsq)*G2);


	//	dOmega1dt += dg1dt * theta1 ;
	//dOmega1dt +=dg1dt / theta1 ;
	//	dOmega1dt += dg1dt *G1 ;
	//	dOmega1dt += dg1dt *sqrt(Hsq)/(G1);



	if(!(TF)){
	  //	fprintf(stderr,"TF=%d\n",TF);
	  //	exit(1);
	//This is from Dan's paper
	  
	tF1=TV1  *pow(a1/R1,8)*SQR(m1)/ ( 9.*(m1+m2)*m2*SQR(1+2.*kL1));
	tF2=TV2  *pow(a1/R2,8)*SQR(m2)/ ( 9.*(m1+m2)*m1*SQR(1+2.*kL2));
	  

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


	W1= ( (1. +15./2. *e1_2 +45./ 8. *e1_4+ 5./16. *e1_6)/pow(1.-SQR(e1),13./2.) - spin1h / ( n )*( 1. + 3.* e1_2 + 3. *e1_4/8.)/pow(1-e1_2,5) ) / tF1  ;

	W2= ( (1. +15./2. *e1_2 +45./ 8. *e1_4+ 5./16. *e1_6)/pow(1.-SQR(e1),13./2.) - spin2h / ( n )*( 1. + 3.* e1_2 + 3. *e1_4/8.)/pow(1-e1_2,5) ) / tF2  ;


	X1= -pow(R1/a1,5) *kL1*m2/(mu *n ) *spin1h*spin1e/(SQR(1-e1_2)) - spin1q / (2. *n *tF1 ) *( 1. +9. /2. * e1_2 + 5./ 8. *e1_4 )/pow(1.-e1_2,5);

	Y1= -pow(R1/a1,5) *kL1*m2/(mu *n ) *spin1h*spin1q/(SQR(1-e1_2)) + spin1e / (2. *n *tF1 ) *( 1. +3. /2. * e1_2 + 1./ 8. *e1_4 )/pow(1.-e1_2,5);

	Z1= pow(R1/a1,5) *kL1*m2/(mu *n ) * ( (2* SQR(spin1h)- SQR(spin1q)- SQR(spin1q) )/ (2.* SQR(1-e1_2) )+ 15. * k2 *m2 / pow(a1,3)*( 1. +3. /2. * e1_2 + 1./ 8. *e1_4 )/pow(1.-e1_2,5));


	X2= -pow(R2/a1,5) *kL2*m1/(mu *n ) *spin2h*spin2e/(SQR(1-e1_2)) - spin2q / (2. *n *tF2 ) *( 1. +9. /2. * e1_2 + 5./ 8. *e1_4 )/pow(1.-e1_2,5);
	
	Y2= -pow(R2/a1,5) *kL2*m1/(mu *n ) *spin2h*spin2q/(SQR(1-e1_2)) + spin2e / (2. *n *tF2 ) *( 1. +3. /2. * e1_2 + 1./ 8. *e1_4 )/pow(1.-e1_2,5);
					     
	Z2= pow(R2/a1,5) *kL2*m1/(mu *n ) * ( (2* SQR(spin2h)- SQR(spin2q)- SQR(spin2q) )/ (2.* SQR(1-e1_2))+ 15. * k2 *m1 / pow(a1,3)*( 1. +3. /2. * e1_2 + 1./ 8. *e1_4 )/pow(1.-e1_2,5));

	//	fprintf(stderr,"X1=%f X2=%f Y1=%f Y2=%f Z1=%f Z2=%f V1=%f V2=%f  W1=%le W2=%le\n",X1,X2,Y1,Y2,Z1,Z2,V1,V2,W1,W2);

	//	ZGR = 3.*k2*(m1+m2)* n/(a1 * SQR(1-e1*e1) *c2); 
	
	detotdt= ( - (V1+ V2) )*(e1);
	dhtotdt +=(  -(W1 +W2)  ) * (htot);


	
	if (fabs(sini1)<1.e-12){
	  dOmega1dtTF = 0.;
	}
	else{
	  dOmega1dtTF = ( (X1+X2)*sing1 + (Y1+Y2)*cosg1 )/(sini1);
	}

	di1dtTF = (X1+X2)*cosg1 - (Y1+Y2)*sing1 ;

	
	//	fprintf(stderr," t=%le di1dtTF=%le i1=%le \n",t,di1dtTF,i1);
	
//	di1dtTF=-sin(i1)*di1dtTF;  //since we are clacualting actualy dcos(i)dt=-sin(i)didt
	
	dg1dtTF = Z1+Z2 - dOmega1dtTF *theta1;
	
	//	de1dt += detotdt;

	//the spin derivative are below!!



	//		fprintf(stderr," t=%le d1spin(e,q,h)=(%le %le %le) d1spin(e,q,h)=(%le %le %le)\n",t,dspin1edt,dspin1qdt,dspin1hdt,dspin2edt,dspin2qdt,dspin2hdt);
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
	dg1dt +=dg1dtTF;
	/* now adding the TF fraction which is actualy di1dt and NOT dcos(i1)dt */

	di1_onlydt = di1_onlydt + di1dtTF;
	
	/* now transform back to dcos(i1)dt */
	di1dt = -di1_onlydt * sin(i1);

	//di1dt +=di1dtTF;

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

	//dxdt[15] = dOmega1dt;	
	//	exit(1);

	for (j=1;j<=n_orb_param;j++){
	  ve[j]=v[j];
	}


	if (t>143819.074087){
	  //	    fprintf(stderr,"t=%le e1=%le de1dt=%le dspin1(%le,%le,%le), did1=%le di1d1=%le dg1dt=%le\n",t,e1,de1dt,dspin1edt,dspin1qdt,dspin1hdt,didt,di1dt,dg1dt);

	}
	
	//if (t>161833.689873 || e1<0.01){
	if (t>143819.074087 && e1<0.1){
	  //  exit(1);
	}

	for (j=1;j<=n_eq;j++){

	  /*	if (t>143819.074087){
	    fprintf(stderr,"t=%le x[%d]=%le dxdt[%d]=%le\n",t,j,dxdt[j]);
	    }
	  */
	  //  dxdt[j]= t* dxdt[j];
	  //  fprintf(stderr," dxdt[%d]=%le\n",j,dxdt[j]);
	  if (dxdt[j]!=dxdt[j]){
	    fprintf(stderr,"t=%le dxdt[%d]=%le\n",t,j,dxdt[j]);
	    flgout=1;
	   
	    exit(1);
	  }
	}
	int m;
	/*	if (flgout==1){
	  	  fprintf(stderr,"t=%le flgout=%d\n",t,flgout);
	  for (m=1;m<=n_eq;m++){
	    dxdt[m]=0;
	  }
	  }*/
	/*	if (t>1e9){
	  fprintf(stderr,"dOmegadt=%le\n",dOmega1dt);
	fprintf(stderr," t=%le d1spin(e,q,h)=(%le %le %le) d1spin(e,q,h)=(%le %le %le)\n",t,dspin1edt,dspin1qdt,dspin1hdt,dspin2edt,dspin2qdt,dspin2hdt);
	}*/
	v[1] = i;
	/*	fprintf(stderr,"   t=%le e2=%le\n",t,e2);
	fprintf(stderr," de2dt = %le dg2dt=%le \n",de2dt,dg2dt);
       	fprintf(stderr," dxdt[3]=%g dxdt[4]=%g\n",dxdt[3],dxdt[4]);
	*/
	/* debug */
	DebugPrint2("dg1dt=%g de1dt=%g ",dxdt[1],dxdt[2]);
	DebugPrint2("dg2dt=%g de2dt%g\n",dxdt[3],dxdt[4]);

}
