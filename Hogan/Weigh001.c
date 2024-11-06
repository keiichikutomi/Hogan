/*WEIGH001.C SINCE 1996.5.26.JUNSATO.*/
/*LAST CHANGE:1998.8.14.*/

/*CMQ OF GIRDER.*/

/*#include <stdio.h>*/
#include <stdlib.h>
#include <string.h>
#include <malloc.h>
#include <math.h>
#include <limits.h>

/*
#ifndef PI
#define PI 3.1415926535897932384
#endif
*/

int cmqpolygon(struct obans cban,double wban,
               struct cmqelem *csum);
int cmqconcentration(double W0,
                     double L0,
                     double L1,double L2,
                     struct cmqelem *cmq);
int cmqrightangledtriangle(double w0,
                           double L0,
                           double L1,double L2,double L3,
                           struct cmqelem *cmq);

int cmqpolygon(struct obans cban,double wban,
               struct cmqelem *csum)
/*CMQ BY POLYGONAL LOAD.INPUT BAN AS POLYGON FORM.*/
/*GIRDER ON 1ST LINE OF BAN.*/
/*wban=WEIGHT OF BAN*/
{
  char non[256];
  int idot,jdot;
  int ndot;
  double w0,L0,L1,L2,L3;
  double H0;
  double **drccos;
  struct cmqelem cmq;
  struct onode *lnods;
  struct obans ban;

  csum->Ci=0.0;
  csum->Cj=0.0;
  csum->Mc0=0.0;
  csum->Qi0=0.0;
  csum->Qj0=0.0;

  ndot=cban.nnod;

  /*PROJECTED BAN*/
  ban.nods=(struct onode **)malloc(ndot*sizeof(struct onode *));
  lnods=(struct onode *)malloc(ndot*sizeof(struct onode));

  drccos=filmdrccos(*(*(cban.nods+0)),
                    *(*(cban.nods+1)),
                    *(*(cban.nods+ndot-1)));    /*DIRECTION COSINE.*/

  for(idot=0;idot<ndot;idot++)       /*TRANSFORM GLOBAL INTO LOCAL.*/
  {
    *(lnods+idot)=transgtol(drccos,*(*(cban.nods+0)),
                                   *(*(cban.nods+idot)));
    *(ban.nods+idot)=(lnods+idot);
  }

  L0=((*(ban.nods+1))->d[EX])-((*(ban.nods+0))->d[EX]);   /*LENGTH.*/

  for(idot=ndot-1;idot>1;idot--)
  {
    H0=(*(ban.nods+idot))->d[EY]; /*Yi=HEIGHT*/
    w0=wban*H0;

    /*          wi          */
    /*          .           */
    /*        .::           */
    /*      .::::           */
    /* ....::::::.......... */
    /*    Xj    Xi          */
    if(idot==ndot-1) jdot=0;
    else             jdot=idot+1;

    L1=(*(ban.nods+jdot))->d[EX];
    L2=(*(ban.nods+idot))->d[EX]-L1;
    L3=L0-L1-L2;
    if(!cmqrightangledtriangle(w0,L0,L1,L2,L3,&cmq)) return 0;

    csum->Ci+=cmq.Ci;
    csum->Cj+=cmq.Cj;
    csum->Mc0+=cmq.Mc0;
    csum->Qi0=cmq.Qi0;
    csum->Qj0=cmq.Qj0;

    /*          wi          */
    /*          .           */
    /*          :.          */
    /*          ::.         */
    /* .........:::........ */
    /*          Xi Xj       */
    jdot=idot-1;

    L1=L0-((*(ban.nods+jdot))->d[EX]);
    L2=L1-((*(ban.nods+idot))->d[EX]);
    L3=L0-L1-L2;
    if(!cmqrightangledtriangle(w0,L0,L1,L2,L3,&cmq)) return 0;

    csum->Ci-=cmq.Cj;
    csum->Cj-=cmq.Ci;
    csum->Mc0+=cmq.Mc0;
    csum->Qi0+=cmq.Qj0;
    csum->Qj0+=cmq.Qi0;
  }


  sprintf(non,"SUM:Ci=%.3f[tfm] Cj=%.3f[tfm] Mc=%.3f[tfm]",
          csum->Ci,csum->Cj,csum->Mc0);
  MessageBox(NULL,non,"CMQ",MB_OK);

  sprintf(non,"SUM:Qi=%.5f[tf] Qj=%.5f[tf]",csum->Qi0,csum->Qj0);
  MessageBox(NULL,non,"CMQ",MB_OK);

  free(lnods);
  free(ban.nods);
  for(idot=2;idot>=0;idot--) free(*(drccos+idot));
  free(drccos);

  return 1;
}/*cmqpolygon*/

int cmqconcentration(double W0,
                     double L0,
                     double L1,double L2,
                     struct cmqelem *cmq)
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
}/*cmqconcentration*/

int cmqrightangledtriangle(double w0,
                           double L0,
                           double L1,double L2,double L3,
                           struct cmqelem *cmq)
/*LOAD:TRIANGLE BETWEEN z=L1 AND z=L1+L2 WITH RIGHT ANGLE.*/
/*L0:LENGTH OF ELEMENT = L1+L2+L3*/
/*w0:PEAK OF TRIANGLE LOAD AT z=L1+L2*/
/*            w0          */
/*            .           */
/*          .::           */
/*        .::::           */
/* ......::::::.......... */
/*   L1    L2      L3     */
/*          L0            */
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
}/*cmqrightangledtriangle*/

