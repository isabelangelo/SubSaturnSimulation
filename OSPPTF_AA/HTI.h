///last update: April 23 - smadar ading tidal friction

//#define TINY 1e-30

#define TSTART 0.//1
#define Tprint 5. // printing interval

/* debugging */
//#define DEBUG 
#ifdef DEBUG 
#define DebugPrint1( f, a )		fprintf( stderr, f, a )
#define DebugPrint2( f, a, b )	fprintf( stderr, f, a, b )
#else
#define DebugPrint1( f, a )				/* print out one value */
#define DebugPrint2( f, a, b )			/* print out two values */
#endif

/* control parameters */
#define n_eq		14					/* number of ODE */
//smadar: adding 4 more integrating param.
#define n_orb_param	8					/* number of orbital parameters */

//#define integrator bsstep    						/* Burlisch-Stoer            */
#define integrator rkqs      						/* Runge-Kutta 4th Order     */
//a paramter that define if the object survivied
double sur;

/*
#define QUADRUPOLE                                 
#define OCTUPOLE
#define GR
//#define TF                                   
*/
int QUADRUPOLE,OCTUPOLE,GR,TF;

//int	  flgout;

/* physical constants */
#define i_koz 0.684719203					/* critical Kozai angle */
#define c2  3.993953869e9					/* c ^ 2 */
#define	k2  (4*M_PI*M_PI)						/* gravitational constant */
#define Rsun (695500*6.68459e-9)                                   /* sun radius au       */


/* tidal parameters*/

#define Q1 0.028 //0.0272   // Beware: Eggleton's Q, not dissipative quality factor.
#define Q2 0.028 //0.3333 
//n=3 politrop - star gives k=0.014 - Q=0.038 and n=1 politrop gas gaiant k=0.25 - Q=0.3333
//"Gyroradii":  I = alpha * M * R^2  
#define alpha1 0.08 
#define alpha2 0.08 //0.26 


//#define TV1 100. //50. //10000000.//100. //50. //1000.//50.//9.3e4 //1000.//1000000000000. //50. //5. //yr Viscous times, which go into tidal dissipation rate.  
//#define TV2 12. //10000000. //12. //0.001 //12 //0.1//100000000000.//0.01 //0.03 //0.001 //5. %%% 0.03 with 
//#define TV1 1.e9///50. //5. //50. //0.01 //5. //15.//100000000000.//2. //50. //10000000.//100. //50. //1000.//50.//9.3e4 //1000.//1000000000000. //50. //5. //yr Viscous times, which go into tidal dissipation rate.  
#define TV1 5.//0.1 //50. //0.001 //50. ///50. //5. //50. //0.01 //5. //15.//100000000000.//2. //50. //10000000.//100. //50. //1000.//50.//9.3e4 //1000.//1000000000000. //50. //5. //yr Viscous times, which go into tidal dissipation rate.  
#define TV2 5.//0.000005 //1. //0.001//1. //1.5 //0.001 //1.5 //5.//100000000000.//0.5 //10000000. //12. //0.001 //12 //0.1//100000000000.//0.01 //0.03 //0.001 //5. %%% 0.03 with 

double Roche;

/* NR built-in functions */
void bsstep(double y[], double dydx[], int nv, double *xx, double htry, double eps,
	double yscal[], double *hdid, double *hnext,
	void (*derivs)(double, double [], double []));
void odeint(double ystart[], int nvar, double x1, double x2, double eps, double h1,
	double hmin, int *nok, int *nbad,
	void (*derivs)(double, double [], double []),
	void (*rkqs)(double [], double [], int, double *, double, double, double [],
	double *, double *, void (*)(double, double [], double [])));
void rkqs(double y[], double dydx[], int n, double *x, double htry, double eps,
	double yscal[], double *hdid, double *hnext,
	void (*derivs)(double, double [], double []));
/* differential equations for a damped harmonic oscillator */
void octupole (double t, double *x, double *dxdt);
