/* CAUTION: This is the traditional K&R C (only) version of the Numerical
   Recipes utility file nrutil.c.  Do not confuse this file with the
   same-named file nrutil.c that is supplied in the same subdirectory or
   archive as the header file nrutil.h.  *That* file contains both ANSI and
   traditional K&R versions, along with #ifdef macros to select the
   correct version.  *This* file contains only traditional K&R.           */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define NR_END 1
#define FREE_ARG char*

void nrerror(error_text)
char error_text[];
/* Numerical Recipes standard error handler */
{
	void exit();

	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}

double *vector(nl,nh)
long nh,nl;
/* allocate a double vector with subscript range v[nl..nh] */
{
	double *v;

	v=(double *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(double)));
	if (!v) nrerror("allocation failure in vector()");
	return v-nl+NR_END;
}

int *ivector(nl,nh)
long nh,nl;
/* allocate an int vector with subscript range v[nl..nh] */
{
	int *v;

	v=(int *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(int)));
	if (!v) nrerror("allocation failure in ivector()");
	return v-nl+NR_END;
}

unsigned char *cvector(nl,nh)
long nh,nl;
/* allocate an unsigned char vector with subscript range v[nl..nh] */
{
	unsigned char *v;

	v=(unsigned char *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(unsigned char)));
	if (!v) nrerror("allocation failure in cvector()");
	return v-nl+NR_END;
}

unsigned long *lvector(nl,nh)
long nh,nl;
/* allocate an unsigned long vector with subscript range v[nl..nh] */
{
	unsigned long *v;

	v=(unsigned long *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(long)));
	if (!v) nrerror("allocation failure in lvector()");
	return v-nl+NR_END;
}

double *dvector(nl,nh)
long nh,nl;
/* allocate a double vector with subscript range v[nl..nh] */
{
	double *v;

	v=(double *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(double)));
	if (!v) nrerror("allocation failure in dvector()");
	return v-nl+NR_END;
}

double **matrix(nrl,nrh,ncl,nch)
long nch,ncl,nrh,nrl;
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	double **m;

	/* allocate pointers to rows */
	m=(double **) malloc((unsigned int)((nrow+NR_END)*sizeof(double*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(double *) malloc((unsigned int)((nrow*ncol+NR_END)*sizeof(double)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

double **dmatrix(nrl,nrh,ncl,nch)
long nch,ncl,nrh,nrl;
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	double **m;

	/* allocate pointers to rows */
	m=(double **) malloc((unsigned int)((nrow+NR_END)*sizeof(double*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(double *) malloc((unsigned int)((nrow*ncol+NR_END)*sizeof(double)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

int **imatrix(nrl,nrh,ncl,nch)
long nch,ncl,nrh,nrl;
/* allocate a int matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	int **m;

	/* allocate pointers to rows */
	m=(int **) malloc((unsigned int)((nrow+NR_END)*sizeof(int*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;


	/* allocate rows and set pointers to them */
	m[nrl]=(int *) malloc((unsigned int)((nrow*ncol+NR_END)*sizeof(int)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

double **submatrix(a,oldrl,oldrh,oldcl,oldch,newrl,newcl)
double **a;
long newcl,newrl,oldch,oldcl,oldrh,oldrl;
/* point a submatrix [newrl..][newcl..] to a[oldrl..oldrh][oldcl..oldch] */
{
	long i,j,nrow=oldrh-oldrl+1,ncol=oldcl-newcl;
	double **m;

	/* allocate array of pointers to rows */
	m=(double **) malloc((unsigned int) ((nrow+NR_END)*sizeof(double*)));
	if (!m) nrerror("allocation failure in submatrix()");
	m += NR_END;
	m -= newrl;

	/* set pointers to rows */
	for(i=oldrl,j=newrl;i<=oldrh;i++,j++) m[j]=a[i]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

double **convert_matrix(a,nrl,nrh,ncl,nch)
double *a;
long nch,ncl,nrh,nrl;
/* allocate a double matrix m[nrl..nrh][ncl..nch] that points to the matrix
declared in the standard C manner as a[nrow][ncol], where nrow=nrh-nrl+1
and ncol=nch-ncl+1. The routine should be called with the address
&a[0][0] as the first argument. */
{
	long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1;
	double **m;

	/* allocate pointers to rows */
	m=(double **) malloc((unsigned int) ((nrow+NR_END)*sizeof(double*)));
	if (!m) nrerror("allocation failure in convert_matrix()");
	m += NR_END;
	m -= nrl;

	/* set pointers to rows */
	m[nrl]=a-ncl;
	for(i=1,j=nrl+1;i<nrow;i++,j++) m[j]=m[j-1]+ncol;
	/* return pointer to array of pointers to rows */
	return m;
}

double ***f3tensor(nrl,nrh,ncl,nch,ndl,ndh)
long nch,ncl,ndh,ndl,nrh,nrl;
/* allocate a double 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
	long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
	double ***t;

	/* allocate pointers to pointers to rows */
	t=(double ***) malloc((unsigned int)((nrow+NR_END)*sizeof(double**)));
	if (!t) nrerror("allocation failure 1 in f3tensor()");
	t += NR_END;
	t -= nrl;

	/* allocate pointers to rows and set pointers to them */
	t[nrl]=(double **) malloc((unsigned int)((nrow*ncol+NR_END)*sizeof(double*)));
	if (!t[nrl]) nrerror("allocation failure 2 in f3tensor()");
	t[nrl] += NR_END;
	t[nrl] -= ncl;

	/* allocate rows and set pointers to them */
	t[nrl][ncl]=(double *) malloc((unsigned int)((nrow*ncol*ndep+NR_END)*sizeof(double)));
	if (!t[nrl][ncl]) nrerror("allocation failure 3 in f3tensor()");
	t[nrl][ncl] += NR_END;
	t[nrl][ncl] -= ndl;

	for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
	for(i=nrl+1;i<=nrh;i++) {
		t[i]=t[i-1]+ncol;
		t[i][ncl]=t[i-1][ncl]+ncol*ndep;
		for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
	}

	/* return pointer to array of pointers to rows */
	return t;
}

void free_vector(v,nl,nh)
double *v;
long nh,nl;
/* free a double vector allocated with vector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_ivector(v,nl,nh)
int *v;
long nh,nl;
/* free an int vector allocated with ivector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_cvector(v,nl,nh)
long nh,nl;
unsigned char *v;
/* free an unsigned char vector allocated with cvector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_lvector(v,nl,nh)
long nh,nl;
unsigned long *v;
/* free an unsigned long vector allocated with lvector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_dvector(v,nl,nh)
double *v;
long nh,nl;
/* free a double vector allocated with dvector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_matrix(m,nrl,nrh,ncl,nch)
double **m;
long nch,ncl,nrh,nrl;
/* free a double matrix allocated by matrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

void free_dmatrix(m,nrl,nrh,ncl,nch)
double **m;
long nch,ncl,nrh,nrl;
/* free a double matrix allocated by dmatrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

void free_imatrix(m,nrl,nrh,ncl,nch)
int **m;
long nch,ncl,nrh,nrl;
/* free an int matrix allocated by imatrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

void free_submatrix(b,nrl,nrh,ncl,nch)
double **b;
long nch,ncl,nrh,nrl;
/* free a submatrix allocated by submatrix() */
{
	free((FREE_ARG) (b+nrl-NR_END));
}

void free_convert_matrix(b,nrl,nrh,ncl,nch)
double **b;
long nch,ncl,nrh,nrl;
/* free a matrix allocated by convert_matrix() */
{
	free((FREE_ARG) (b+nrl-NR_END));
}

void free_f3tensor(t,nrl,nrh,ncl,nch,ndl,ndh)
double ***t;
long nch,ncl,ndh,ndl,nrh,nrl;
/* free a double f3tensor allocated by f3tensor() */
{
	free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
	free((FREE_ARG) (t[nrl]+ncl-NR_END));
	free((FREE_ARG) (t+nrl-NR_END));
}
/*****************************************************/
   
void linint(double *xa,double *ya,int n,double x,double *y)
     /*
       From Numerical Recipes, but modified to do linear interpolation 
       (and extrapolation) instead of cubic spline, and to deal with 
       < or > ordering.
     */
{
  int k,khi=n/*-1*/,klo=1/*0*/;
  double a,b,h;
  /*
   printf("in linint\n"); 
   printf("x=%le\n",exp(x));
  printf("n=%d\n",n);
  printf("khi=%d klo=%d\n",khi,klo);  
 printf("xa[khi]=%le xa[klo]=%le\n",xa[khi],xa[klo]);
   printf("y[khi]=%le y[klo]=%le\n",ya[khi],ya[klo]);
  */  
  //  fprintf(stderr,"x=%le xa[1]=%le xa[n]=%le\n",pow(10,x),pow(10,xa[klo]),pow(10,xa[khi]));

if (xa[khi] >= xa[klo])
    {
      while (khi-klo > 1) 
	{
	  k=(khi+klo)/2;
	  if (xa[k] > x) khi=k;
	  else klo=k;
	}
      //  printf("k=%d\n",k);
    }
  else
    {
      while (khi-klo > 1) 
	{
	  k=(khi+klo)/2;
	  if (xa[k] < x) khi=k;
	  else klo=k;
	}
    }
  h=xa[khi]-xa[klo];
  //  printf("h=%le\n",h);
  if (h == 0.) 
    {
       fprintf(stderr,"======\n");
       fprintf(stderr,"%g %g %g\n",xa[1],xa[n],x);
      nrerror("bad xa input in linint B");
    }
  a=(xa[khi]-x)/h;
  b=(x-xa[klo])/h;
  //printf("a=%le b=%le\n",a,b);
  *y=a*ya[klo]+b*ya[khi];
  //  printf("y=%le\n",*y);
  // return;
}


/****************************************************/
/*****************************************************/
   
void linintMinMax(double *xa,double *ya,int minv,int n,double x,double *y)
     /*
       From Numerical Recipes, but modified to do linear interpolation 
       (and extrapolation) instead of cubic spline, and to deal with 
       < or > ordering.
     */
{
  int k,khi=n/*-1*/,klo=minv/*0*/;
  double a,b,h;
  if (klo>khi){
    fprintf(stderr,"======\n");
    fprintf(stderr," oh no, the min index of the vecor %d is larger from the max one %d\n",klo,khi);
    nrerror("fix index - exsiting ... \n");
  }

  /*
   printf("in linint\n"); 
   printf("x=%le\n",exp(x));
  printf("n=%d\n",n);
  printf("khi=%d klo=%d\n",khi,klo);  
 printf("xa[khi]=%le xa[klo]=%le\n",xa[khi],xa[klo]);
   printf("y[khi]=%le y[klo]=%le\n",ya[khi],ya[klo]);
  */  
  //  fprintf(stderr,"x=%le xa[1]=%le xa[n]=%le\n",pow(10,x),pow(10,xa[klo]),pow(10,xa[khi]));

if (xa[khi] >= xa[klo])
    {
      while (khi-klo > 1) 
	{
	  k=(khi+klo)/2;
	  if (xa[k] > x) khi=k;
	  else klo=k;
	}
      //  printf("k=%d\n",k);
    }
  else
    {
      while (khi-klo > 1) 
	{
	  k=(khi+klo)/2;
	  if (xa[k] < x) khi=k;
	  else klo=k;
	}
    }
  h=xa[khi]-xa[klo];
  //  printf("h=%le\n",h);
  if (h == 0.) 
    {
       fprintf(stderr,"======\n");
       fprintf(stderr,"%g %g %g\n",xa[minv],xa[n],x);
      nrerror("bad xa input in linint B");
    }
  a=(xa[khi]-x)/h;
  b=(x-xa[klo])/h;
  //printf("a=%le b=%le\n",a,b);
  *y=a*ya[klo]+b*ya[khi];
  //  printf("y=%le\n",*y);
  // return;
}


/****************************************************/
