/****************************************************************************
 * Integrate octupole perturbation equations for heirarchical triple systems 
 * Takes input from "triple.in" 
 * Update: 29 June, 2006                                                    
 ****************************************************************************/

#include <stdio.h> 
#include <stdlib.h> 
#include <string.h> 
#include <math.h>
#include "nrutil.h"
#include "HTI.h"

/* odeint control parameters */
int kmax = 10000000,    kount = 0;
double *v, *xp, **yp, dxsav;						/* storage arrays */


/* main program, calls for odeint.c */
int main (void) {
	
	double 	*x, *x_i, *dxdt;
	double 	t_i, t_f, eps=0, h1=0, hmin=0;
	int		nok=0, nbad=0, j=0, k_f=0;
	int		scnret;
	char	dummy[1024];
	double	m1,m2,m3,a1,a2,e1,e2,i,g1,g2,emax,e1max_calc;
	double	P1,P2,P_koz,P_GR,range;
	FILE*	input;  
	double imax=0.,imin=360.,e1max=0,e1min=1.,e2max=0,e2min=1.,i1max=0.;//smadar: what is the min and max inclinations reached? 
	double i1,i2;
	double L1,L2,G1,G2,Hsq,theta,m12,m123,im,theta1,theta2;
	double t_inter,tloop_f;
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

	fscanf(input, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf ",\
					&m1,&m2,&m3,&a1,&a2,&e1,&e2,&g1,&g2,&i,&t_f);
	for(;;) {	/* read in until :::, end of header */
		scnret = fscanf(input, "%s", dummy);
		if (!(strcmp(dummy, ":::"))) break;
	}
	
/* read in control parameters */

	fscanf(input, "%lf",&eps);

	fclose(input);
#ifdef QUADRUPOLE
	fprintf(stderr,"QUADRUPOLE\n");
#endif
#ifdef OCTUPOLE
	fprintf(stderr,"OCTUPOLE\n");
#endif
#ifdef GR
	fprintf(stderr,"GR\n");
#endif

	i *=M_PI/180.0;
	g1*=M_PI/180.0;
	g2*=M_PI/180.0;
	e1max_calc = sqrt(1-5/3.0*cos(i)*cos(i));
	t_f *= 1e6;
	t_i = 0;

	

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
	fprintf(stderr," G1=%le G2=%le\n",G1,G2);

	theta1=(Hsq-G2*G2+G1*G1)/(2.*sqrt(Hsq)*G1);
	theta2=(Hsq+G2*G2-G1*G1)/(2.*sqrt(Hsq)*G2);

	if(fabs(theta1)<=1)    i1=acos(theta1);
	else if(theta1>1)      i1 = 0.0;
	else if(theta1<1)      i1 = M_PI;

	if(fabs(theta2)<=1)    i2=acos(theta2);
	else if(theta2>1)      i2 = 0.0;
	else if(theta2<1)      i2 = M_PI;

	im=i; //this is a redandancy - just to check if all is OK
	/* various timescales */
	P1 = 2*M_PI*sqrt(a1*a1*a1/(k2*(m1+m2)));
	P2 = 2*M_PI*sqrt(a2*a2*a2/(k2*(m1+m2+m3)));
	P_koz = P1*(m1+m2)/m3*pow(a2/a1,3)*sqrt(1-e2*e2)*(1-e2*e2);
	P_GR  = 3.36e7*(1-e1*e1)*P1*a1/m2;
	range = 5/3.0*cos(i)*cos(i);	//to find minimium P_GR 
	P_GR*=range;

	/* memory allocation */
	x = 	vector(1,n_eq);
	x_i = 	vector(1,n_eq);
	dxdt = 	vector(1,n_eq);
	v = 	vector(1,n_orb_param);
	xp = 	vector(1,kmax);
	yp 	=	matrix(1,n_eq,1,kmax);
	
	/* control parameters */

	h1  	= eps * P2;
	hmin = 0;
	dxsav 	= 10  * P2;

	/* initial condition */
	x_i[1] = e1*sin(g1);
	x_i[2] = e1*cos(g1);
	x_i[3] = e2*sin(g2);
	x_i[4] = e2*cos(g2);
	x_i[5] = cos(i1);
	x_i[6] = cos(i2);
	x_i[7] = cos(i);
	v[1]	= i;
	v[2] 	= m1;
	v[3] 	= m2;
	v[4] 	= m3;
	v[5]	= a1;
	v[6]	= a2;
		

	/* print out the initial conditions */
	fprintf(stderr,"\n============================== Initial Conditions ==============================\n");
	fprintf(stderr,"	m1=%1.3f m2=%1.3f m3=%1.3f a1=%1.3f a2=%1.3f k2=%1.3f\n"\
				,m1,m2,m3,a1,a2,k2);
	fprintf(stderr,"	e1=%1.3f g1=%1.5f e2=%1.3f g2=%1.5f i=%1.3f\n"\
				,e1,g1/M_PI*180,e2,g2/M_PI*180,i/M_PI*180);
	fprintf(stderr,"	i1=%1.3f i2=%1.5f \n"\
				,i1/M_PI*180,i2/M_PI*180);
	fprintf(stderr,"	Integration stopped at %.0f years, Inner Period=%1.0g years, Outer Period=%e years\n"\
				,t_f,P1, P2);
	if(i>i_koz)  		fprintf(stderr,"	e1_max calc'd from i= %lf\n",e1max_calc);
	else if(i<i_koz)	fprintf(stderr, "	Inclination too low\n");
	fprintf(stderr, "	Kozai Period= %1.5e years\n	GR Period= %.5e years---%.5e years\n\n"\
				,P_koz,P_GR*range, P_GR);
	

	fprintf(stderr,"\n epsilon=%lg\n",a1*e2/(a2*(1.-e2*e2)));
	fprintf(stderr,"\n epsilonM=%lg\n",a1*e2/(a2*(1.-e2*e2))*(m1-m2)/(m1+m2));	
	/******************************************************************/
	/* integrate the octupole equations (use bsstep rather than rkqs) */
	/* if ETrigG is on: odeint will output in this order
		time, e1sing1, e1cosg1, e2sing2, e2cosg2, i , e1, e2  cos(i1) cos(i2) cos(i)*/
	/* deviding the  integration to smaller steps to prevent overshoting. */
	/* in very long integration times (i.e., t_f is very large) try to devide to bigger numbers  */ 
	//	t_inter=fabs(t_f-t_i)/5000;
		t_inter=fabs(t_f-t_i)/90000;
	//	t_inter=fabs(t_f-t_i)/800000;
	while (t_i<=t_f){
	  tloop_f=t_i+t_inter;
	  h1=fabs(tloop_f-t_i)/1000000.; //initiial step size - decreas in NaN
	  odeint(x_i, n_eq, t_i, tloop_f, eps, h1, hmin, &nok, &nbad, octupole, integrator);
	  k_f = kount;

	/* output result */
       
	for(j=1;j<=k_f;j++) {
	  
	  printf("%lf\t",xp[j]);
	  
	  e1 = sqrt(yp[1][j]*yp[1][j] + yp[2][j]*yp[2][j]);
	  e2 = sqrt(yp[3][j]*yp[3][j] + yp[4][j]*yp[4][j]);
	  g1 = atan2(yp[1][j],yp[2][j]);
	  g2 = atan2(yp[3][j],yp[4][j]);

	
	  m12 = m1 + m2;
	  m123 = m1 + m2 + m3;
	L1 = m1 * m2  / m12  * sqrt( k2 * m12  * a1 );
	L2 = m3 * m12 / m123 * sqrt( k2 * m123 * a2 );
	G1 = L1 * sqrt(1 - e1*e1);
	G2 = L2 * sqrt(1 - e2*e2);
	if(xp[j]==0)    Hsq = G1*G1 + 2.*G1*G2*cos(im) + G2*G2;
	theta = (Hsq - G1*G1 - G2*G2)/(2.0*G1*G2);
	if(fabs(theta)<=1)    im=acos(theta);
	else if(theta>1)      im = 0.0;
	else if(theta<1)      im = M_PI;


	if(fabs(yp[5][j])<=1)    i1=acos(yp[5][j]);
	else if(yp[5][j]>1)      i1 = 0.0;
	else if(yp[5][j]<1)      i1 = M_PI;

	if(fabs(yp[6][j])<=1)    i2=acos(yp[6][j]);
	else if(yp[6][j]>1)      i2 = 0.0;
	else if(yp[6][j]<1)      i2 = M_PI;



	if(fabs(yp[7][j])<=1)    i=acos(yp[7][j]);
	else if(yp[7][j]>1)      i = 0.0;
	else if(yp[7][j]<1)      i = M_PI;


	printf("%1.6lf \t %1.6lf \t %1.6lf \t %1.6lf  \t %1.6lf \t %1.6lf \t %1.6lf  \t %1.6lf \t %1.6lf \t %1.6lf \t %1.6lf \t %1.6lf\n",e1,e2,g1,g2,i1*180./M_PI,i2*180./M_PI,i*180./M_PI,im*180./M_PI,G1*yp[5][j],G1,G2,G1*G1 + 2.*G1*G2*cos(i) + G2*G2);
	

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
	  if( i1*180./M_PI>i1max){
	    i1max=i1*180./M_PI;
	  }

       
	}// end of  printing loop on the steps that were saved for odeint

	//initialize again for the next step of integration::
	x_i[1] = yp[1][k_f];
	x_i[2] =  yp[2][k_f];
	x_i[3] =  yp[3][k_f];
	x_i[4] = yp[4][k_f];
	x_i[5]=yp[5][k_f];
	x_i[6]=yp[6][k_f];
	x_i[7]=yp[7][k_f];
	t_i=xp[k_f];
	// parameters:
	v[1]	= i;
	v[2] 	= m1;
	v[3] 	= m2;
	v[4] 	= m3;
	v[5]	= a1;
	v[6]	= a2;
	}//end of integrating loop
	fprintf(stderr,"	    -------\n");

	

	fprintf(stderr,"	    emin=%f emax=%f i1max=%le\n",e1min,e1max,i1max);
  
	
	/* free memory */
	free_vector(x,1,2);
	free_vector(x_i,1,2);
	free_vector(dxdt,1,2);
	free_vector(v,1,n_orb_param);
	free_vector(xp,1,kmax);
	free_matrix(yp,1,kmax,1,kmax);

	return 0;
}

