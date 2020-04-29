/*gcc -o read_sse_cl read_sse_cl.c */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int Ndismax = 100; /*maximum number of type changes allowed */
int Nlinmax = 5000; /*maximum number of lines from sse output */

void readssecl(double mass, int *Ndis, int *Nlin, double *tdis, double *time, double *type, double *Mo, double *Mt, double *logL, double *logR, double *logT, double *Mc, double *Menv, double *epoch, double *spin)
{

  FILE *fp;
  int status;
  char command[100], output[200];
  char *pch;
  int i,pos;

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
  i=0;
  *Ndis = -1;
  *Nlin = 0;
  while (fgets(output, sizeof(output)-1, fp) != NULL) {
    if (i == 0) *Ndis = atof(output);
    if (i > 0 && i < *Ndis+1) tdis[i-1] = atof(output);
    if (*Ndis != -1 && i > *Ndis+1) {
      (*Nlin)++;
      pos = i-*Ndis-2;
      pch = strtok(output," ");
      time[pos] = atof(pch);
      pch = strtok(NULL," ");
      type[pos] = atof(pch);
      pch = strtok(NULL," ");
      Mo[pos] = atof(pch);
      pch = strtok(NULL," ");
      Mt[pos] = atof(pch);
      pch = strtok(NULL," ");
      logL[pos] = atof(pch);
      pch = strtok(NULL," ");
      logR[pos] = atof(pch);
      pch = strtok(NULL," ");
      logT[pos] = atof(pch);
      pch = strtok(NULL," ");
      Mc[pos] = atof(pch);
      pch = strtok(NULL," ");
      Menv[pos] = atof(pch);
      pch = strtok(NULL," ");
      epoch[pos] = atof(pch);
      pch = strtok(NULL," ");
      spin[pos] = atof(pch);
    }
    i++;
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

 
}

int main( int argc, char *argv[] )
{
  double mass;

  mass = 7.;

  /* define the arrays and initialize */
  double tdis[Ndismax];
  double time[Nlinmax],type[Nlinmax],Mo[Nlinmax],Mt[Nlinmax],logL[Nlinmax],logR[Nlinmax],logT[Nlinmax],Mc[Nlinmax],Menv[Nlinmax],epoch[Nlinmax],spin[Nlinmax];
  int Ndis,Nlin;
  int i;
  for (i=0; i<Ndismax; i++) tdis[i] = -1.;
  for (i=0; i<Nlinmax; i++){
    time[i] = -1;
    type[i] = -1;
    Mo[i] = -1;
    Mt[i] = -1;
    logL[i] = -1;
    logR[i] = -1;
    logT[i] = -1;
    Mc[i] = -1;
    Menv[i] = -1;
    epoch[i] = -1;
    spin[i] = -1;
  }

  //printf("==before func\n");
  readssecl(mass, &Ndis, &Nlin, tdis, time, type, Mo, Mt, logL, logR, logT, Mc, Menv, epoch, spin);
  // printf("==before func\n");

  /* test the output */
  fprintf(stdout,"%i %i\n",Ndis,Nlin);
  for (i=0; i<Ndis; i++) fprintf(stdout,"%f\n",tdis[i]); 
  for (i=0; i<Nlin; i++) fprintf(stdout,"%f %f %f %f %f %f %f %f %f %f %le\n",time[i],type[i],Mo[i],Mt[i],logL[i],logR[i],logT[i],Mc[i],Menv[i],epoch[i],spin[i]); 


  return 0;
}
