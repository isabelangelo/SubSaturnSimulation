#include <math.h>
#include "HTI.h"
#include <stdio.h>
#include <stdlib.h>
#include "nrutil.h"

double *tvec,*mvec,*Rvec,*dmdtvec,*tofdmdtvec;
int length;
extern 	int Ndis,Nlin;
extern double *tdis;
extern double *timeSSE,*MtSSE,*logRSSE; //*typeSSE,*MoSSE,*MtSSE,*logLSSE,*logRSSE,*logTSSE,*McSSE,*MenvSSE,*epochSSE,*spinSSE;
extern int *indexdis; //this is the index of discontinutites

void linint(double *xa,double *ya,int n,double x,double *y);
double dfridr(double (*func)(double), double x, double h, double *err);

double moft(double t){
  int count,min,max,length,indexval;
  double moftval;
 
  min=1;
  max=Nlin;

  // fprintf(stderr,"=================\n");
  if (t<timeSSE[1]){
    moftval=MtSSE[1];
 
    return moftval;
  }
  else{ 
    for (count=1; count<=Ndis; count++){
    
      if (t>=tdis[Ndis]){

	min=indexdis[Ndis];
	max=Nlin;

	break;
      }
  
      if (t>tdis[count]){
	min=indexdis[count];
	count++;

      }
      else {
	max=indexdis[count-1];
	break;
      }
    }

  if (max==min) {
    fprintf(stderr,"Something funny max=min=%d you did something wrong!\n",max);
  }
 
 
  length=max-min;

  indexval=indexdis[Ndis];
  
   if (length==1 || length==2){ 
     moftval=MtSSE[min];
   }
  
   else{
  tvec=vector(1,length);
  mvec=vector(1,length);
 
  //fill new vercotrs
  
  for (count=1;count<=length;count++){
    mvec[count]=log10(MtSSE[count+min]);
    tvec[count]=log10(timeSSE[count+min]);
   
  }
  // fprintf(stderr,"t=%le tvec[1]=%le tvec[length]=%le\n",t,pow(10,tvec[1]),pow(10,tvec[length]));
  linint(tvec,mvec,length,log10(t),&moftval);
  moftval=pow(10.,moftval);

  free_vector(mvec,1,length);
  free_vector(tvec,1,length);
  
   }
  
  
  
  return moftval;
  }
}
/*******************************************/
double moftNew(double t){
  double mval;
  int j,count;
  double *mvec,*tvec;

  if (t<=tdis[5]){
    linint(timeSSE,MtSSE,Nlin,tdis[5],&mval);
  }
  if (t>tdis[5] && t<tdis[6]){
    tvec=vector(1,Nlin);
    mvec=vector(1,Nlin);
    j=1;
    for (count=1;count<=Nlin;count++){
      if (timeSSE[count]>tdis[5] && timeSSE[count]<tdis[6]) {
	mvec[j]=(MtSSE[count]);
	tvec[j]=(timeSSE[count]);
	
	j++;
      }
      if (timeSSE[count]>tdis[6]) break;
    }
    linint(tvec,mvec,j-1,t,&mval);
    free_vector(mvec,1,Nlin);
    free_vector(tvec,1,Nlin);
  }
  else {
	   linint(timeSSE,MtSSE,Nlin,tdis[6],&mval);
	  }
  return mval;
}

/*******************************************/
double dmdt(double t, double m, double m0){
  double mf;
  double dm=0.;
  double err=1.e-8;
  double Deltat;
  int indexval,min,max,length,count,j;
  t=t+tMS;
  indexval=indexdis[Ndis];

  if (!(SSE)){

    min=1;
    max=Nlin;
    /*
    // fprintf(stderr,"=================\n");
    if (t<timeSSE[1]){
      dm=0;
 
    }
    else{ 
      for (count=1; count<=Ndis; count++){
    
	if (t>=tdis[Ndis]){

	  min=indexdis[Ndis];
	  max=Nlin;

	  break;
	}
        if (t>tdis[count]){
	  min=indexdis[count];
	  count++;
	}
	else {
	  max=indexdis[count-1];
	  break;
	}
      }

      if (max==min) {
	fprintf(stderr,"Something funny max=min=%d you did something wrong!\n",max);
      }
    
 
      length=max-min;
      dmdtvec=vector(1,length-1);
      // tvec=vector(1,length);
      // mvec=vector(1,length);
      tofdmdtvec=vector(1,length-1);
      //fill new vercotrs
      
      for (count=1;count<=length;count++){
	//  mvec[count]=(MtSSE[count+min]);
	//  tvec[count]=(timeSSE[count+min]);
	if (count>=2){
	  dmdtvec[count-1]=( MtSSE[count] - MtSSE[count-1])/(timeSSE[count]-timeSSE[count-1]);
	  tofdmdtvec[count-1]=(timeSSE[count]+timeSSE[count-1])/2.;
	}
      }
      
      linint(tofdmdtvec,dmdtvec,length-1,(t),&dm);
      
      //  free_vector(mvec,1,length);
      //  free_vector(tvec,1,length);
      free_vector(dmdtvec,1,length-1);
      free_vector(tofdmdtvec,1,length-1);
    }
    */

if (t<=tdis[5]){
             dm=0;
	}
	if (t>tdis[5] && t<tdis[6]){
	  /*
	 === put here derivative from SSE
	 use t instead of the index, i.e., if t >timeSSE[iSSE] and t< time[iSSE+1} etc - also use the mass astimated at tdis[5] so it will work well with the discontinurty you impose!*/
	  linint(timeSSE,MtSSE,Nlin,tdis[6],&mf);
	  dm=(mf-m0)/(tdis[6]-tdis[5]); //mf=1.280100e+00 m0=6.745018e+00 dm=-1.063028e-05
	  /*
	  dmdtvec=vector(1,Nlin);
	  //    tvec=vector(1,Nlin);
	  //       mvec=vector(1,Nlin);
      tofdmdtvec=vector(1,Nlin);
      j=1;
	  for (count=1;count<=Nlin;count++){
	    if (timeSSE[count]>tdis[5] && timeSSE[count]<tdis[6]) {
		mvec[j]=(MtSSE[count]);
		tvec[j]=(timeSSE[count]);

		j++;
	      }
	  */
	      /*
	    if (timeSSE[count]>tdis[5] && timeSSE[count]<tdis[6] && count>2){
	      dmdtvec[j]=( MtSSE[count] - MtSSE[count-1])/(timeSSE[count]-timeSSE[count-1]);
	      tofdmdtvec[j]=(timeSSE[count]+timeSSE[count-1])/2.;
	      j++;
	    }
	      */
	  /*
	    if (timeSSE[count]>tdis[6]) break;
	  }
	    //	  linint(tofdmdtvec,dmdtvec,j-1,(t),&dm);// dm=-6.281882e-07

	  //	  fprintf(stderr," dm=%le\n",dm);
	    //  free_vector(mvec,1,Nlin);
	    //        free_vector(tvec,1,Nlin);
	  free_vector(dmdtvec,1,Nlin);
      free_vector(tofdmdtvec,1,Nlin);
*/

	    Deltat=tdis[6]-tdis[5];
	    err=1.e-8;
	    dm=dfridr(&moftNew,t,Deltat,&err);


	}
	else {
	  dm = 0.;
	  }



    /*
   if (t<=tdis[5]){
             dm=0;
	}
	if (t>tdis[5] && t<tdis[6]){
	  mf=0.109*m0+0.394;
	  dm=(mf-m0)/(tdis[6]-tdis[5]);
	  fprintf(stderr,"mf=%le m0=%le dm=%le\n",mf,m0,dm);
	}
	else {
	  dm = 0.;
	  }*/

  }
    else{
       
      /*if (t<=3e6){
             dm=0;
	}
	if (t>3e6 && t<4e6){
	  mf=0.109*m0+0.394;
	  dm=(mf-m0)/(1.e6);
	}
	else {
	  dm = 0.;
	}
	}*/
      if (t<=55.77067e6){
             dm=0;
	}
	if (t>55.77067e6 && t<56.28476e6){
	  mf=0.109*m0+0.394;
	  dm=(mf-m0)/(56.28476e6-55.77067e6);
	}
	else {
	  dm = 0.;
	}
	}
    //fprintf(stderr,"dmdt mass/yera =%le \n",dm);
    //   exit(1);
	return dm;
 }

double dmdtOLD(double t, double m, double m0){
  double mf;
  double dm=0.;
  double err=1.e-8;
  double Deltat;
  int indexval;
  t=t+tMS;
  indexval=indexdis[Ndis];
  
    if (!(SSE)){


      if (t<tdis[Ndis]){
   
	if (t<tdis[Ndis-1]){
	  Deltat=(tdis[Ndis-1]-tdis[Ndis-2])*1.2;
	}
	else{
	  Deltat=(tdis[Ndis]-tdis[Ndis-1])*1.2;	
	}
     

	dm=dfridr(&moft,t,Deltat,&err);
      }
      else{
	dm=0.;
      }


  }
    else{
        if (t<=3e6){
             dm=0;
	}
	if (t>3e6 && t<4e6){
	  mf=0.109*m0+0.394;
	  dm=(mf-m0)/(1.e6);
	}
	else {
	  dm = 0.;
	}
	 }

    //fprintf(stderr,"dmdt mass/yera =%le \n",dm);
    //   exit(1);
	return dm;
}

/*******************************************/

double dRdt(double t, double m, double dmdt){
double dR;
  dR=0;
	return dR;
}
