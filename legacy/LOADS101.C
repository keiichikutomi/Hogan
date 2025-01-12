/*LOADS001.C SINCE 1996.5.26.JUNSATO.*/
/*LAST CHANGE:1996.9.6.*/

/*CMQ OF GIRDER FOR KYOTO ENTRANCE.*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <malloc.h>
#include <math.h>
#include <limits.h>
#include <time.h>

#ifndef PI
#define PI 3.1415926535897932384
#endif

struct structcmq{double Ci,Cj,Mc0,Qi0,Qj0;};

int cmqconcentration(double W0,
                     double L0,
                     double L1,double L2,
                     struct structcmq *cmq);
int cmqrightangledtriangle(double w0,
                           double L0,
                           double L1,double L2,double L3,
                           struct structcmq *cmq);

void main()
{
  char non[20];
  double W0,w0,L0,L1,L2,L3;
  double H0;
  struct structcmq cmq,sum={0.0,0.0,0.0,0.0,0.0};

  /*       .              */
  /*     .:::.            */
  /*   .:::::::.  H0=L0/2 */
  /* .:::::::::::.        */
  /*       L0             */

  L0=6.300;
  H0=L0/2.0;

  w0=2.86*H0;
  L1=0.0;
  L2=H0;
  L3=H0;
  cmqrightangledtriangle(w0,L0,L1,L2,L3,&cmq);
  printf("L0=%.3f[m] PEAK LOAD=%.5f[tf/m]\n",L0,w0);
  printf("Ci=%.5f[tfm] Cj=%.5f[tfm] Mc=%.5f[tf] ",
         cmq.Ci,cmq.Cj,cmq.Mc0);
  printf("Qi=%.5f[tf] Qj=%.5f[tf]\n",cmq.Qi0,cmq.Qj0);

  sum.Ci=cmq.Ci-cmq.Cj; sum.Cj=-cmq.Ci+cmq.Cj;
  sum.Mc0+=2.0*cmq.Mc0;
  sum.Qi0=cmq.Qi0+cmq.Qj0; sum.Qj0=cmq.Qi0+cmq.Qj0;

  printf("SUM:");
  printf("Ci=%.5f[tfm] Cj=%.5f[tfm] Mc=%.5f[tf] ",
         sum.Ci,sum.Cj,sum.Mc0);
  printf("Qi=%.5f[tf] Qj=%.5f[tf]\n",sum.Qi0,sum.Qj0);

  gets(non);
  return;
}/*main*/

int cmqconcentration(double W0,
                     double L0,
                     double L1,double L2,
                     struct structcmq *cmq)
/*LOAD:CONCENTRATED.*/
/*L0:LENGTH OF ELEMENT = L1+L2*/
/*W0:LOAD AT z=L1*/
{
  if(L0!=(L1+L2)) return 0;

  cmq->Ci=W0*L1*L2*L2/(L0*L0);
  cmq->Cj=-W0*L1*L1*L2/(L0*L0);

  cmq->Qi0=W0*L2/L0;
  cmq->Qj0=W0*L1/L0;

  if(L1>=L2) cmq->Mc0=0.5*W0*L2;
  if(L1<L2)  cmq->Mc0=0.5*W0*L1;

  return 1;
}/*cmqconcentration.......................................*/

int cmqrightangledtriangle(double w0,
                           double L0,
                           double L1,double L2,double L3,
                           struct structcmq *cmq)
/*LOAD:TRIANGLE BETWEEN z=L1 AND z=L1+L2 WITH RIGHT ANGLE.*/
/*L0:LENGTH OF ELEMENT = L1+L2+L3*/
/*w0:PEAK OF TRIANGLE LOAD AT z=L1+L2*/
{
  if(L0!=(L1+L2+L3)) return 0;

  cmq->Ci=1.0/60.0*w0*L2/(L0*L0)
          *(5.0*(L2*L2+4.0*L2*L3+ 6.0*L3*L3)*L1
           +2.0*(L2*L2+5.0*L2*L3+10.0*L3*L3)*L2);
  cmq->Cj=-1.0/60.0*w0*L2/(L0*L0)
          *(5.0*( 6.0*L1*L1+ 8.0*L1*L2+3.0*L2*L2)*L3
               +(10.0*L1*L1+10.0*L1*L2+3.0*L2*L2)*L2);

  cmq->Qi0=1.0/6.0*w0*L2/L0*(L2+3.0*L3);
  cmq->Qj0=1.0/6.0*w0*L2/L0*(3.0*L1+2.0*L2);

  if(0.5*L0<=L1) cmq->Mc0=0.5*L0*cmq->Qi0;
  if(L1<0.5*L0 && 0.5*L0<(L1+L2))
  {
    cmq->Mc0=1.0/48.0*w0/L2
             *(3.0*L2*L2*(L1+L2+3.0*L3)
               +(L1-L3)*(L1-L3)*(L1-3.0*L2-L3));
  }
  if((L1+L2)<=0.5*L0) cmq->Mc0=0.5*L0*cmq->Qj0;

  return 1;
}/*cmqrightangledtriangle.................................*/

