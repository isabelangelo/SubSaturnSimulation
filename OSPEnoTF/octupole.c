#include <math.h>
#include "HTI.h"
#include <stdio.h>
#include "nrutil.h"


extern double *v;


/* ODE for a damped harmonic oscillator to be integrated */
void octupole (double t, double *x, double *dxdt) {

	static double m1,m2,m3,a1,a2,P1,P2;
	static double L1,L2, C3G23, C4G25, m12, m123, Hsq;
	static int	  flg=0;
	
	double e1sing1, e1cosg1, e2sing2, e2cosg2;
	double de1sing1dt, de1cosg1dt, de2sing2dt, de2cosg2dt;
	double g1,g2,e1,e2,i=0.0;                          /* varrying quantities */
	double dg1dt, dg2dt, de1dt, de2dt;                 /* derivatives         */
	double G1,G2,theta,A,B,C3,C4,cosphi;               /* ang.mom., etc.      */
	double sing1, cosg1, sing2, cosg2, sin2g1, cos2g1; /* trig funcs          */

	double dFdg2=0.,dFdg1=0.,di1dt=0.,di2dt=0.,didt=0.,dBigG1dt=0.,dBigG2dt=0.;
	double i1,i2;
	double n;
	int j;
	/* copy values from the vector to variables */
	e1sing1 = x[1];
	e1cosg1 = x[2];
	e2sing2 = x[3];
	e2cosg2 = x[4];

	e1 = sqrt(e1sing1*e1sing1+e1cosg1*e1cosg1);
	e2 = sqrt(e2sing2*e2sing2+e2cosg2*e2cosg2);
	g1 = atan2(e1sing1,e1cosg1);
	g2 = atan2(e2sing2,e2cosg2);


	if(!(flg))		{
		i =v[1];
		m1=v[2];
		m2=v[3];
		m3=v[4];
		a1=v[5];
		a2=v[6];
		
		m12 = m1 + m2;
		m123 = m1 + m2 + m3;
		
		P1 = 2*M_PI*sqrt(a1*a1*a1/(k2*(m1+m2)));
		P2 = 2*M_PI*sqrt(a2*a2*a2/(k2*(m1+m2+m3)));
    
		L1 = m1 * m2  / m12  * sqrt( k2 * m12  * a1 );
		L2 = m3 * m12 / m123 * sqrt( k2 * m123 * a2 );
    
		//	C3G23 = k2*k2/16.0 * pow(m12*m3,7)/pow(m123*m1*m2,3) * pow(L1/L2,3) * L1;
		//wrong sign for the quadrapole:
		C3G23 =  k2*k2/16.0 * pow(m12*m3,7)/pow(m123*m1*m2,3) * pow(L1/L2,3) * L1;
		C4G25 = -k2*k2*(15/64.0) * pow(m12*m3,9)*(m1-m2)/pow(m123,4)/pow(m1*m2,5) * pow(L1*L1/L2,3);
		
		/* debug */
		DebugPrint2("L1=%g L2=%g  ",L1, L2);
		DebugPrint2("C3G23=%g C4G25=%g\n",C3G23,C4G25);
	
		flg++;
	}

	n = 2.*M_PI/P1;

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
	
	if(t==0)    Hsq = G1*G1 + 2.*G1*G2*cos(i) + G2*G2;
	
	theta = (Hsq - G1*G1 - G2*G2)/(2.0*G1*G2);

	if(fabs(theta)<=1)    v[1]=acos(theta);
	else if(theta>1)      v[1] = 0.0;
	else if(theta<1)      v[1] = M_PI;
	
	theta=x[7];
	
	if(fabs(theta)<=1)    i=acos(theta);
	else if(theta>1)      i = 0.0;
	else if(theta<1)      i = M_PI;


	if(fabs(x[5])<=1)    	i1 = acos(x[5]);
	else if(x[5]>1)      i1 = 0.0;
	else if(x[5]<1)       i1 = M_PI;


	if(fabs(x[6])<=1)    	i2 = acos(x[6]);
	else if(x[6]>1)      i2 = 0.0;
	else if(x[6]<1)       i2 = M_PI;
    
	B = 2 + 5*e1*e1 - 7*e1*e1*cos2g1;
	A = 4 + 3*e1*e1 - 2.5*(1-theta*theta)*B;
	C3 = C3G23 / pow(G2,3);
	C4 =  C4G25 / pow(G2,5);

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



#ifdef QUADRUPOLE                  /* quadrupole approximation */
	
	dg1dt += C3 * 6 * ((1/G1)*(4*theta*theta+(5*cos2g1-1)*(1-e1*e1-theta*theta) ) \
       		+ (theta/G2)*(2+e1*e1*(3-5*cos2g1)));
			
	de1dt += C3 * ((1-e1*e1)/G1) * (30*e1*(1-theta*theta)*sin2g1);

	dg2dt += C3 * 3 * ((2*theta/G1)*(2+e1*e1*(3-5*cos2g1) ) +\
				(1/G2)*(4+6*e1*e1 + (5*theta*theta-3)*(2+e1*e1*(3-5*cos2g1))));

	de2dt += 0.;

	dFdg1 += C3*15. * SQR(e1)* (1 - SQR(theta))* ( -2.*sin2g1) ; 

	dFdg2 += 0.;
	

#endif


#ifdef OCTUPOLE                   /* octupole approximation */


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
#endif

#ifdef GR						/* include general reletavistic precession */
	dg1dt += 6*M_PI*k2*m12/(P1*a1*(1.0-e1*e1)*c2);
	dg2dt += 6*M_PI*k2*m123/(P2*a2*(1.0-e2*e2)*c2);
	//fprintf(stderr,"  IN GR \n");

#endif


	dBigG1dt += dFdg1;
	dBigG2dt += dFdg2;

	


	//note: This is actualy dcos(i)dt!!!
	di1dt=-dBigG1dt * (Hsq-G2*G2+G1*G1)/(2.*sqrt(Hsq)*G1*G1) +(-2.*G2*dBigG2dt+2.*G1*dBigG1dt )/(2.*sqrt(Hsq)*G1);
	di2dt=-dBigG2dt * (Hsq+G2*G2-G1*G1)/(2.*sqrt(Hsq)*G2*G2) +(+2.*G2*dBigG2dt-2.*G1*dBigG1dt )/(2.*sqrt(Hsq)*G2);

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




	dxdt[1] =  e1cosg1 * dg1dt + sing1*de1dt;
	dxdt[2] = -e1sing1 * dg1dt + cosg1*de1dt;
	dxdt[3] =  e2cosg2 * dg2dt + sing2*de2dt;
	dxdt[4] = -e2sing2 * dg2dt + cosg2*de2dt;
	

	dxdt[5] = di1dt;
	dxdt[6] = di2dt;
	dxdt[7] = didt;

	for (j=1;j<=n_eq;j++){
	  if (dxdt[j]!=dxdt[j]){
	    fprintf(stderr,"t=%le dxdt[%d]=%le\n",t,j,dxdt[j]);
	    fprintf(stderr," Oh-oh you go NaNs, maybe try to have a large resultion run....\n");
	    fprintf(stderr,"   --- Terminating now...\n");
	    exit(1);
	  }
	}

	v[1] = i;
	/*	fprintf(stderr,"   t=%le e2=%le\n",t,e2);
	fprintf(stderr," de2dt = %le dg2dt=%le \n",de2dt,dg2dt);
       	fprintf(stderr," dxdt[3]=%g dxdt[4]=%g\n",dxdt[3],dxdt[4]);
	*/
	/* debug */
	DebugPrint2("dg1dt=%g de1dt=%g ",dxdt[1],dxdt[2]);
	DebugPrint2("dg2dt=%g de2dt%g\n",dxdt[3],dxdt[4]);

}
