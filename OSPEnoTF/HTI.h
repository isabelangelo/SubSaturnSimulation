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
#define n_eq		7					/* number of ODE */
#define n_orb_param	6					/* number of orbital parameters */

#define integrator bsstep    						/* Burlisch-Stoer            */
//#define integrator rkqs      						/* Runge-Kutta 4th Order     */

#define QUADRUPOLE                                 
#define OCTUPOLE
#define GR

/* physical constants */
#define i_koz 0.684719203					/* critical Kozai angle */
#define c2  3.993953869e9					/* c ^ 2
											 * 2.599358e8 for earlier PSR1620 runs
											 * CHECK!!!!                           */
#define	k2  (4*M_PI*M_PI)						/* gravitational constant */



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

