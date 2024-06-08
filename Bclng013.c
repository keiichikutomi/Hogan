/*=================================================================*/
/*     PROGRAM BCLNG001                                            */
/*     ELASTIC BUCKLING FOR 3D FRAME                               */
/*     CODED BY JUN SATO                                           */
/*     DATE 1993.7.28 ... 1994.1.28                                */
/*     BASED ON 'Y6FRAMEB' BY YOSHINOBU FUJITANI                   */
/*=================================================================*/
/*                                                                 */
/* DEIGAB:GENERALIZED EIGENVALUE PROBLEM.                          */
/* DEIGRS:STANDARD EIGENVALUE PROBLEM.                             */
/* LAST DEBUG:1997.11.15.                                          */
/*                                                                 */
/* MODIFY HOUSEHOLDER. 1997.10.18.                                 */
/* MODIFY BISECTION. 1997.10.18.                                   */
/*                                                                 */
/* [A]{x}=L{x} CORRECT SOLUTION FOR FOLLOWING [A].                 */
/*                                                                 */
/* |[A]-L[I]|=|1-L  2 | BY ARRAY.1997.10.12.                       */
/*            | 2  5-L|                                            */
/*                                                                 */
/*            |1-L  2   3 |                                        */
/* |[A]-L[I]|=| 2  5-L  2 | BY ARRAY.1997.10.18.                   */
/*            | 3   2  1-L|                                        */
/*                                                                 */
/* [A]{x}=L[B]{x} CORRECT SOLUTION FOR FOLLOWING [A],[B].          */
/* SYMMETRIC BY FULL ARRAY WITHOUT CONFS:1997.10.19.               */
/* SYMMETRIC BY LOWER TRIANGLE GCOMP WITHOUT CONFS:1997.11.09.     */
/* BY GCOMP WITH ONE CONF:CONF=1 FOR ALL CASE.1997.11.15.          */
/*                                                                 */
/* 4x4 MATRIX                                                      */
/* [A]=  2             [B]= 2             {CONFS}= 0               */
/*       4 10              -1  4                   0               */
/*       0 12 14            2 -3  9                0               */
/*       8  0 16 18        -4  1 -2 11             0               */
/*       1  1  1  1  1      1  1  1  1  1          1               */
/*                                                                 */
/* FOR GCOMP STYLE WITHOUT CONFS                                   */
/*   CHOLESKI DECOMPOSITION COMPLETED. 1997.10.19.                 */
/*   WITHOUT CONFS COMPLETED. 1997.11.09.                          */
/*   4x4 MATRIX WITH ALL CONFS=FREE COMPLETED. 1997.11.15.         */
/*                                                                 */
/* CORRECT SOLUTION BY GCOMP FOR FOLLOWING [A],[B].1997.11.16.     */
/*                                                                 */
/* [A]= 2                    [B]= 5                  {CONFS}= 0    */
/*      4 10                     -1  8                        0    */
/*      1  0  1                   1  0  1                     1    */
/*      0  0  1  1                0  0  1  1                  1    */
/*      0 12  0  0 14             0 -6  1  0  9               0    */
/*      8  0  0  1  0 18         -9  0  0  1 -1 19            0    */
/*      3  0  1  0 16  5 15       3 -2  1  0  0 -5 21         0    */
/*      1  1  0  0  1  0  1  1    1  1  0  0  1  1  1  1      1    */
/*                                                                 */

#include <windows.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <malloc.h>
#include <math.h>
#include <limits.h>
#include <time.h>

#define MSIZE   8 /*MATRIX SIZE.*/
#define NEIGEN  5 /*NUMBERS OF EIGEN VALUES TO BE OBTAINED.*/

struct gcomponent{unsigned short int m,n;                /*LINE,ROW*/
                  double value;
                  struct gcomponent *down,*left;};  /*SIZE 20 BYTES*/
struct confinement{signed char iconf;
                   double value;};                  /*CONFINEMENTS.*/
/*-----------------------------------------------------------------*/
struct gcomponent *gdefine(unsigned short int m,
                           unsigned short int n,
                           double value,
                           struct gcomponent *down,
                           struct gcomponent *left);
int gread(struct gcomponent *gmtx,
          long int i,long int j,double *data);
int gwrite(struct gcomponent *gmtx,
           long int i,long int j,double data);
struct gcomponent *copygcompmatrix(struct gcomponent *gmtx,
                                   long int msize);
double eigensubstitution(struct gcomponent *amtx,
                         struct gcomponent *bmtx,
                         struct confinement *confs,
                         long int msize,
                         double eig,double *vct);
/*-----------------------------------------------------------------*/
void currentvalues(char *str,
                   long int n,long int ne,
                   double A[][MSIZE],
                   double W[][MSIZE],
                   double E[],double V[][MSIZE]);
void deigab(double A[][MSIZE],double B[][MSIZE],
            long int N,long int NSIZE,long int NE,long int NV,
            double EPS,double W[][MSIZE],
            double E[],double V[][MSIZE]);
void deigrs(double A[][MSIZE],
            long int N,long int N1,long int NE,long int NV,
            double EPS,double W[][MSIZE],double LW[],
            double E[],double V[][MSIZE]);
/*-----------------------------------------------------------------*/
void currentvalue(char *str,
                  long int n,long int ne,
                  struct gcomponent *A,
                  double **W,
                  double *E,double **V);
void deigabgeneral(struct gcomponent *A,
                   struct gcomponent *B,
                   struct confinement *confs,
                   long int N,long int NE,long int NV,
                   double EPS,
                   double *E,double **V);
void deigrsstandard(struct gcomponent *A,
                    struct confinement *confs,
                    long int N,long int NE,long int NV,
                    double EPS,
                    double *E,double **V);
/*-----------------------------------------------------------------*/

void main(void)
{
  char non[10],str[256],s[80];
  int i,j,k;
  int ii,jj;
  long int lc,nko;
  long int loff;
  double eps=1.0E-16;
  double ai,bi,data;
  double GK[MSIZE][MSIZE],GG[MSIZE][MSIZE];
  double GKmem[MSIZE][MSIZE],GGmem[MSIZE][MSIZE];
  double V[NEIGEN][MSIZE],W[7][MSIZE],XL[NEIGEN];

  double **v,*vv,*xl;
  struct gcomponent *kmtx,*gmtx; /*GLOBAL MATRIX*/
  struct gcomponent *kmem,*gmem;
  struct gcomponent ginit={0,0,0.0,NULL,NULL};

  struct confinement *confs;

  lc=MSIZE;
  nko=NEIGEN;

  loff=MSIZE;                         /*DIAGONALS OF GLOBAL MATRIX.*/
  kmtx=(struct gcomponent *)malloc(loff*sizeof(struct gcomponent));
  gmtx=(struct gcomponent *)malloc(loff*sizeof(struct gcomponent));
  confs=(struct confinement *)
        malloc(loff*sizeof(struct confinement));
  if(kmtx==NULL || gmtx==NULL) exit(1);
  for(i=0;i<=loff-1;i++)            /*GLOBAL MATRIX INITIALIZATION.*/
  {
    (kmtx+i)->down=NULL;
    /*(kmtx+i)->left=NULL;*/
    (gmtx+i)->down=NULL;
    /*(gmtx+i)->left=NULL;*/

    (confs+i)->iconf=0;
    (confs+i)->value=0.0;
  }

  loff=MSIZE;                       /*GLOBAL MATRIX INITIALIZATION.*/
  for(i=1;i<=loff;i++)
  {
    ginit.m=(unsigned short int)i;
    ginit.n=(unsigned short int)i;

    *(kmtx+(i-1))=ginit;
    *(gmtx+(i-1))=ginit;
  }

  /* [K]=  2                   [G]= 5                  {CONFS}= 0 */
  /*       4 10                    -1  8                        0 */
  /*       1  0  1                  1  0  1                     1 */
  /*       0  0  1  1               0  0  1  1                  1 */
  /*       0 12  0  0 14            0 -6  1  0  9               0 */
  /*       8  0  0  1  0 18        -9  0  0  1 -1 19            0 */
  /*       3  0  1  0 16  5 15      3 -2  1  0  0 -5 21         0 */
  /*       1  1  0  0  1  0  1  1   1  1  0  0  1  1  1  1      1 */

  (confs+2)->iconf=1;
  (confs+3)->iconf=1;
  (confs+7)->iconf=1;

  gwrite(kmtx,1,1, 2.0);   gwrite(gmtx,1,1, 5.0);
  gwrite(kmtx,2,1, 4.0);   gwrite(gmtx,2,1,-1.0);
  gwrite(kmtx,3,1, 1.0);   gwrite(gmtx,3,1, 1.0);
  gwrite(kmtx,6,1, 8.0);   gwrite(gmtx,6,1,-9.0);
  gwrite(kmtx,7,1, 3.0);   gwrite(gmtx,7,1, 3.0);
  gwrite(kmtx,8,1, 1.0);   gwrite(gmtx,8,1, 1.0);

  gwrite(kmtx,2,2,10.0);   gwrite(gmtx,2,2, 8.0);
  gwrite(kmtx,5,2,12.0);   gwrite(gmtx,5,2,-6.0);
                           gwrite(gmtx,7,2,-2.0);
  gwrite(kmtx,8,2, 1.0);   gwrite(gmtx,8,2, 1.0);

  gwrite(kmtx,3,3, 1.0);   gwrite(gmtx,3,3, 1.0);
  gwrite(kmtx,4,3, 1.0);   gwrite(gmtx,4,3, 1.0);
                           gwrite(gmtx,5,3, 1.0);
  gwrite(kmtx,7,3, 1.0);   gwrite(gmtx,7,3, 1.0);

  gwrite(kmtx,4,4, 1.0);   gwrite(gmtx,4,4, 1.0);
  gwrite(kmtx,6,4, 1.0);   gwrite(gmtx,6,4, 1.0);

  gwrite(kmtx,5,5,14.0);   gwrite(gmtx,5,5, 9.0);
                           gwrite(gmtx,6,5,-1.0);
  gwrite(kmtx,7,5,16.0);
  gwrite(kmtx,8,5, 1.0);   gwrite(gmtx,8,5, 1.0);

  gwrite(kmtx,6,6,18.0);   gwrite(gmtx,6,6,19.0);
  gwrite(kmtx,7,6, 5.0);   gwrite(gmtx,7,6,-5.0);
                           gwrite(gmtx,8,6, 1.0);

  gwrite(kmtx,7,7,15.0);   gwrite(gmtx,7,7,21.0);
  gwrite(kmtx,8,7, 1.0);   gwrite(gmtx,8,7, 1.0);

  gwrite(kmtx,8,8, 1.0);   gwrite(gmtx,8,8, 1.0);

  fprintf(stderr,"[GK]\n");
  for(ii=1;ii<=MSIZE;ii++) /*CONFIRM GLOBAL MATRIX*/
  {
    sprintf(str,"\0");
    for(jj=1;jj<=ii;jj++)
    {
      gread(kmtx,ii,jj,&data);
      sprintf(s," %9.3f",data);
      strcat(str,s);
    }
    fprintf(stderr,"%s\n",str);
  }

  fprintf(stderr,"[GG]\n");
  for(ii=1;ii<=MSIZE;ii++) /*CONFIRM GLOBAL MATRIX*/
  {
    sprintf(str,"\0");
    for(jj=1;jj<=ii;jj++)
    {
      gread(gmtx,ii,jj,&data);
      sprintf(s," %9.3f",data);
      strcat(str,s);
    }
    fprintf(stderr,"%s\n",str);
  }

  sprintf(str,"CONFS\n");
  for(ii=0;ii<MSIZE;ii++)
  {
    fprintf(stderr," %9d",(confs+ii)->iconf);
  }
  fprintf(stderr,"\n");

  gets(non);

  xl=(double *)malloc(NEIGEN*sizeof(double));
  v=(double **)malloc(NEIGEN*sizeof(double *));
  for(i=0;i<NEIGEN;i++)
  {
    vv=(double *)malloc(MSIZE*sizeof(double));
    *(v+i)=vv;
  }

  kmem=copygcompmatrix(kmtx,MSIZE);
  gmem=copygcompmatrix(gmtx,MSIZE);

  deigabgeneral(kmtx,gmtx,confs,MSIZE,nko,nko,eps,xl,v);

  for(i=0;i<NEIGEN;i++)
  {
    fprintf(stderr,"SOLUTION %ld\n",i+1);
    eigensubstitution(kmem,gmem,confs,MSIZE,*(xl+i),*(v+i));
    gets(non);
  }

  /*BY ARRAY.*/
  /* [K]=  2                  [G]= 2             */
  /*       4 10                   -1  4          */
  /*       0 12 14                 2 -3  9       */
  /*       8  0  0 18             -4  1 -2 11    */
  /*       3  0 16  5 15           3 -2  0 -1 21 */

  GK[0][0]= 2.0; GK[0][1]= 4.0; GK[0][2]= 0.0; GK[0][3]= 8.0;
  GK[1][0]= 4.0; GK[1][1]=10.0; GK[1][2]=12.0; GK[1][3]= 0.0;
  GK[2][0]= 0.0; GK[2][1]=12.0; GK[2][2]=14.0; GK[2][3]= 0.0;
  GK[3][0]= 8.0; GK[3][1]= 0.0; GK[3][2]= 0.0; GK[3][3]=18.0;
  GK[4][0]= 3.0; GK[4][1]= 0.0; GK[4][2]=16.0; GK[4][3]= 5.0;

  GK[0][4]= 3.0;
  GK[1][4]= 0.0;
  GK[2][4]=16.0;
  GK[3][4]= 5.0;
  GK[4][4]=15.0;

  GG[0][0]= 2.0; GG[0][1]=-1.0; GG[0][2]= 2.0; GG[0][3]=-4.0;
  GG[1][0]=-1.0; GG[1][1]= 4.0; GG[1][2]=-3.0; GG[1][3]= 1.0;
  GG[2][0]= 2.0; GG[2][1]=-3.0; GG[2][2]= 9.0; GG[2][3]=-2.0;
  GG[3][0]=-4.0; GG[3][1]= 1.0; GG[3][2]=-2.0; GG[3][3]=11.0;
  GG[4][0]= 3.0; GG[4][1]=-2.0; GG[4][2]= 0.0; GG[4][3]=-1.0;

  GG[0][4]= 3.0;
  GG[1][4]=-2.0;
  GG[2][4]= 0.0;
  GG[3][4]=-1.0;
  GG[4][4]=21.0;

  for(i=0;i<5;i++)
  {
    for(j=0;j<5;j++)
    {
      GKmem[i][j]=GK[i][j];
      GGmem[i][j]=GG[i][j];
    }
  }

  /*deigab(GK,GG,lc,MSIZE,-nko,nko,eps,W,XL,V);*/
  deigab(GK,GG,5,5,5,5,eps,W,XL,V);

  /*CHECK SOLUTIONS.*/

  for(i=0;i<5;i++)
  {
    fprintf(stdout,"[A]{v}=e{v}\n");
    for(j=0;j<5;j++)
    {
      ai=0.0;
      bi=0.0;
      for(k=0;k<5;k++)
      {
        ai+=GKmem[j][k]*V[i][k];
        bi+=XL[i]*GGmem[j][k]*V[i][k];
      }

      fprintf(stdout," {A%dT}{v}=%9.5f ",(j+1),ai);
      fprintf(stdout,"L{B%dT}{v}=%9.5f\n",(j+1),bi);
    }
    gets(non);
  }

  return;
}/*main*/

/*-----------------------------------------------------------------*/
void currentvalues(char *str,
                   long int n,long int ne,
                   double A[][MSIZE],
                   double W[][MSIZE],
                   double E[],double V[][MSIZE])
/*CHECK CURRENT VALUES FOR DEBUG.*/
{
  char non[10];
  long int i,j;

  ne=labs(ne);
/*
  for(i=0;i<n;i++)
  {
    fprintf(stdout,"A%d:",i+1);
    for(j=0;j<n;j++)
    {
      fprintf(stdout," %12.5f",A[i][j]);
    }
    fprintf(stdout,"\n");
  }
*/
  fprintf(stdout,"\n");
  for(i=1;i<=ne;i++) fprintf(stdout,"           E%ld",i);
  fprintf(stdout,"\n");
  for(j=0;j<ne;j++)
  {
    fprintf(stdout," %12.5f",E[j]);
  }
  fprintf(stdout,"\n");

  fprintf(stdout,"\n");
  for(i=1;i<=ne;i++) fprintf(stdout,"           V%ld",i);
  fprintf(stdout,"\n");
  for(i=0;i<n;i++)
  {
    for(j=0;j<ne;j++)
    {
      fprintf(stdout," %12.5f",V[j][i]);
    }
    fprintf(stdout,"\n");
  }

/*
  fprintf(stdout,"\n");
  for(i=0;i<7;i++)
  {
    fprintf(stdout,"W%d:",i);
    for(j=0;j<n;j++)
    {
      fprintf(stdout," %12.5E",W[i][j]);
    }
    fprintf(stdout,"\n");
  }
*/

  fprintf(stdout,"%s\n",str);
  gets(non);
  if(strcmp(non,"\0")) exit(1);
  return;
}/*currentvalues*/

/*-----------------------------------------------------------------*/
void deigab(double A[][MSIZE],double B[][MSIZE],
            long int N,long int NSIZE,long int NE,long int NV,
            double EPS,double W[][MSIZE],
            double E[],double V[][MSIZE])
{
  /* SUBROUTINE FOR GENERALIZED EIGENVALUE PROBLEM                */
  /*  [A]{X} = LAMBDA [B]{X}                                      */
  /* FOR REAL SYMMETRIC MATRICES [A] & [B], THE LATTER BEING      */
  /* POSITIVE DEFINITE.                                           */
  /*                                                              */
  /* USAGE:                                                       */
  /* CALL DEIGAB( A, B, N, NSIZE, NE, NV, EPS, W, E, V )          */
  /*                                                              */
  /* INPUT:                                                       */
  /*   A[N][N]  : 2-DIM. ARRAY CONTAINING REAL SYMMETRIC MATRIX.  */
  /*   B[N][N]  : 2-DIM. ARRAY CONTAINING REAL SYMMETRIC          */
  /*              POSITIVE DEFINITE MATRIX.                       */
  /*   N        : ORDER OF MATRIX.                                */
  /*   NSIZE    : SIZE OF THE 2-DIM. ARRAYS  A, B, W  & V         */
  /*              DEFINED IN 'DIMENSION' STATEMENT (FIRST INDEX). */
  /*   NE       : NUMBER OF EIGENVALUES TO BE OBTAINED.           */
  /*              IN DESCENDING ORDER  WHEN  NE > 0,              */
  /*              IN ASCENDING  ORDER  WHEN  NE < 0.              */
  /*   NV       : NUMBER OF EIGENVECTORS TO BE OBTAINED.          */
  /*   EPS      : ACCURACY  ( STANDARD VALUE  1.0D-16 ).          */
  /*                                                              */
  /* OUTPUT:                                                      */
  /*   E[NE]    : 1-DIM. ARRAY CONTAINING THE OUTPUT EIGENVALUES  */
  /*   V[NV][N] : 2-DIM. ARRAY CONTAINING THE OUTPUT EIGENVECTORS */
  /*              THE VECTOR ( V(1,K), V(2,K),..., V(N,K) )       */
  /*              BELONGS TO THE EIGENVALUE  E(K).                */
  /*                                                              */
  /* WORKING SPACE:                                               */
  /*   W[7][N]  : 2-DIM. ARRAY  (7,N)  USED AS THE WORKING AREA.  */
  /*                                                              */
  /* NOTE:                                                        */
  /*  THE MATRICES [A] AND [B] ARE DESTROYED.                     */
  /* METHOD:                                                      */
  /*  CHOLESKY REDUCTION FOLLOWED BY HOUSEHOLDER DIAGONALIZATION. */
  /* SUBROUTINE USED : DEIGRS                                     */

  /*char non[256];*/
  long int I,J,K,K1,nev,neva,nvec,R,RSUB1;
  double SUM,T1,T;

  /* CHECK INPUT DATA. */
  nev  = NE;
  neva = labs(nev);
  nvec = NV;

  if(N<=0 || NSIZE-N<0 || neva==0 || N-neva<0 ||
     nvec<0 || neva-nvec<0)
  {
    fprintf(stdout,"DEIGAB:INVALID ARGUMENT.");
    fprintf(stdout," N,NSIZE,NE,NV = %ld,%ld,%ld,%ld",
            N,NSIZE,nev,nvec);
    return;
  }

  if(N==1)
  {
    T = B[0][0];
    if(T<=0.0)
    {
      fprintf(stdout,"DEIGAB:MATRIX [B] IS NOT POSITIVE DEFINITE.");
    }
    else
    {
      B[0][0] = sqrt(1.0/T);
      E[0] = A[0][0]/T;
      V[0][0] = 1.0;
    }
    return;
  }

  /* CHOLESKI DECOMPOSITION OF THE POSITIVE DEFINITE */
  /* MATRIX [B] INTO A PRODUCT OF A LOWER TRIANGULAR */
  /* MATRIX [L] WITH ITS TRANSPOSED MATRIX.          */
  /* DIAGONALS OF [L] ARE REMAINED INVERSE.          */

  T = B[0][0];
  if(T<=0.0)
  {
    fprintf(stdout,"DEIGAB:MATRIX [B] IS NOT POSITIVE DEFINITE.");
    return;
  }
  T = sqrt(1.0/T);
  B[0][0] = T;

  for(I=1;I<N;I++)
  {
    B[0][I] = B[I][0] * T;
  }

  for(R=1;R<N;R++)
  {
    SUM = 0.0;
    for(K=0;K<=R-1;K++)
    {
      SUM = B[K][R]*B[K][R]  + SUM;
    }
    T = B[R][R] - SUM;
    if(T<=0.0)
    {
      fprintf(stdout,"DEIGAB:MATRIX [B] IS NOT POSITIVE DEFINITE.\n");
      fprintf(stdout,"T%ld=%9.5E-%9.5E=%9.5E\n",R,B[R][R],SUM,T);
      return;
    }
    T = sqrt(1.0/T);
    B[R][R] = T;

    if(R>=N-1) break;

    for(I=R+1;I<N;I++)
    {
      SUM = 0.0;
      for(K=0;K<=R-1;K++)
      {
        SUM = B[K][I] * B[K][R] + SUM;
      }
      B[R][I] = ( B[I][R] - SUM ) * T;
    }
  }

  /*currentvalues("DEIGAB:[B] DECOMPOSED.",N,NE,B,W,E,V);*/

  /* CHOLESKI DECOMPOSITION IS COMPLETED.                */
  /* [B] = [L][LT],  THE UPPER TRIANGULAR MATRIX [LT] IS */
  /* STORED IN THE ARRAY [B].                            */
  /* NOW, PROCEED TO EVALUATE  [L]**(-1) [A] [LT]**(-1)  */
  /* FIRST, PREMULTIPLY  [L]**(-1).                      */

  for(J=0;J<N;J++)
  {
    A[0][J] = A[J][0] * B[0][0];
  }
  for(J=0;J<N;J++)
  {
    for(R=1;R<N;R++)
    {
      RSUB1 = R - 1;
      SUM = 0.0;
      for(K=0;K<=RSUB1;K++)
      {
        SUM = B[K][R] * A[K][J]  +  SUM;
      }
      A[R][J] = ( A[R][J] - SUM ) * B[R][R];
    }
  }

  /* NEXT, POSTMULTIPLY  [LT]**(-1). */

  for(J=0;J<N;J++)
  {
    A[J][0] = A[J][0] * B[0][0];
  }
  for(R=1;R<N;R++)
  {
    RSUB1 = R - 1;
    T1 = B[R][R];
    for(K=0;K<=RSUB1;K++)
    {
      T = - B[K][R];
      for(J=R;J<N;J++)
      {
        A[J][R] += A[J][K] * T;
      }
    }
    for(J=R;J<N;J++)
    {
      A[J][R] = A[J][R] * T1;
    }
  }

  /*currentvalues("DEIGAB:TRANSFORMED.",N,NE,A,W,E,V);*/

  /*FIND EIGENVALUES AND EIGENVECTORS OF THE TRANSFORMED MATRIX.*/
  deigrs( A, N, NSIZE, nev, nvec, EPS, W, W[6], E, V );

  /*BACK TRANSFORMATION OF EIGENVECTORS.*/
  if(nvec==0) return;

  for(J=0;J<nvec;J++)
  {
    K = N-1;
    while(1)
    {
      T = V[J][K] * B[K][K];
      V[J][K] = T;
      K1 = K - 1;
      if(K1<0) break;
      for(R=0;R<=K1;R++)
      {
        V[J][R] -= B[R][K] * T;
      }
      K = K1;
    }
  }
  currentvalues("DEIGAB END.",N,NE,A,W,E,V);

  return;
}/*deigab*/

/*-----------------------------------------------------------------*/
void deigrs(double A[][MSIZE],
            long int N,long int N1,long int NE,long int NV,
            double EPS,double W[][MSIZE],double LW[],
            double E[],double V[][MSIZE])
{
  /* SUBROUTINE FOR STANDARD EIGENVALUE PROBLEM                   */
  /* HOUSEHOLDER'S TRIDIAGONAL REDUCTION.                         */
  /* EIGENVALUES BY BISECTION.                                    */
  /* EIGENVECTORS BY INVERSE ITERATION.                           */
  /*                                                              */
  /*  [A]{X} = LAMBDA {X}                                         */
  /*                                                              */
  /* INPUT:                                                       */
  /*   A[N][N]  : 2-DIM. ARRAY CONTAINING REAL SYMMETRIC MATRIX.  */
  /*   N        : ORDER OF MATRIX.                                */
  /*   N1       : SIZE OF THE 2-DIM. ARRAYS  A, B, W  & V         */
  /*              DEFINED IN 'DIMENSION' STATEMENT (FIRST INDEX). */
  /*   NE       : NUMBER OF EIGENVALUES TO BE OBTAINED.           */
  /*              IN DESCENDING ORDER  WHEN  NE > 0,              */
  /*              IN ASCENDING  ORDER  WHEN  NE < 0.              */
  /*   NV       : NUMBER OF EIGENVECTORS TO BE OBTAINED. NV<=|NE| */
  /*              ONLY EIGENVALUES IF NV=0.                       */
  /*   EPS      : ACCURACY  ( STANDARD VALUE  1.0D-16 ).          */
  /*                                                              */
  /*   W[6][N]  : WORK SPACE FOR TRIDIAGONALS,ITERATION etc.      */
  /*   LW[N]    :                                                 */
  /*                                                              */
  /* OUTPUT:                                                      */
  /*   E[NE]    : 1-DIM. ARRAY CONTAINING THE OUTPUT EIGENVALUES  */
  /*   V[NV][N] : 2-DIM. ARRAY CONTAINING THE OUTPUT EIGENVECTORS */
  /*              THE VECTOR { V(K,1), V(K,2),..., V(K,N) }       */
  /*              BELONGS TO THE EIGENVALUE  E(K).                */

  /*char non[256];*/
  int SW,of;
  long int I,I1,J,K,K1,MM,NEA,NVA,nm1,nm2,M;
  double T,R,RR1,RR2,D,F,Q,S,SR,A1,EPS1,EPS2;

  NEA=labs(NE);
  if(NEA==0)
  {
    fprintf(stdout,"DEIGRS: NE = %ld\n",NE);
    fprintf(stdout,
            "NE SHOULD BE NON-ZERO. RETURN WITH NO CALCULATION.\n");
    return;
  }
  NVA=labs(NV);
  if(NVA>NEA || NEA>N || N>N1)
  {
    fprintf(stdout,"DEIGRS: NV,NE,N,N1 = %ld,%ld,%ld,%ld\n",
            NV,NE,N,N1);
    fprintf(stdout,"NV,NE,N,N1");
    fprintf(stdout,
            " SHOULD SATISFY THE FOLLOWING INEQUALITIES.\n");
    fprintf(stdout,
            "|NV|<=|NE|<=N<=N1 RETURN WITH NO CALCULATION.\n");

    E[0]=0.0;
    return;
  }

  nm1=N-1;
  nm2=N-2;
  if(EPS<0.0) EPS=1.0E-16;

  /* CASE N=1 : ONLY ONE NUMBER. */
  if(nm2<0)
  {
    E[0]=A[0][0];
    if(NV!=0) V[0][0] = 1.0;
    return;
  }

  /* CASE N=2 : COMPUTE EIGENVALUES OF 2x2 MATRIX. */
  else if(nm2==0)
  {
    W[0][0]=A[1][0];
    T = 0.5*(A[0][0]+A[1][1]);
    R=A[0][0]*A[1][1]-A[1][0]*A[1][0];
    D=T*T-R;
    Q=fabs(T)+sqrt(D);
    if(T<0.0) Q=-Q;
    T = T*(double)NE;
    if(T>=0.0)
    {
      E[0]=Q;
      if(NEA==2) E[1]=R/Q;
    }
    else
    {
      E[0]=R/Q;
      if(NEA==2) E[1]=Q;
    }
  }

  /* CASE N=3,4,... */
  /* REDUCE TO TRIDIAGONAL FORM BY HOUSEHOLDER'S METHOD */
  else if(nm2>0)
  {
    for(K=0;K<nm2;K++) /* N-2 TIMES REDUCTION. */
    {
      K1=K+1;
      S=0.0;
      for(I=K1;I<N;I++) S=S+A[I][K]*A[I][K]; /* S={aT}{a} */

      W[0][K]=0.0;
      if(S!=0.0)
      {
        SR=sqrt(S); /* SR=SQRT( {aT}{a} ) */

        A1=A[K1][K];
        if(A1<0.0) SR=-SR; /* SGN(SR)=SGN(A1) */
        W[0][K]=-SR; /* a1=-SR */

        R = 1.0/(S+A1*SR);

        A[K1][K]=A1+SR; /*{V}={a(K+1) a(K+2)...a(N)}+{s 0 0...0}*/

        for(I=K1;I<N;I++)
        {
          S=0.0;
          for(J=K1;J<I;J++)
          {
            S+=A[J][K]*A[I][J]; /*Si={VT}{Ai}*/
          }
          for(J=I;J<N;J++)
          {
            S+=A[J][K]*A[J][I];
          }
          W[0][I]=S*R; /*{Woi}=R{Si}*/
        }

        S=0.0;
        for(I=K1;I<N;I++) S+=W[0][I]*A[I][K]; /*S'={VT}{Wo}*/

        T = 0.5*R*S;
        for(I=K1;I<N;I++)
        {
          W[0][I]-=T*A[I][K]; /*{Wo'}={Wo}-0.5RS'{V}*/
        }

        for(J=K1;J<N;J++) /*[Aij']=[Aij-Wj'Vi-Wi'Vj]*/
        {
          for(I=J;I<N;I++)
          {
            A[I][J]-=W[0][J]*A[I][K]+W[0][I]*A[J][K];
          }
        }
      }
    }
    W[0][nm1-1]=A[N-1][nm1-1];

    /*currentvalues("DEIGRS:REDUCED.",N,NE,A,W,E,V);*/

    /*COMPUTE EIGENVALUES BY BISECTION METHOD*/

    for(I=0;I<N;I++) W[5][I]=A[I][I]; /*{W5}={Aii}*/

    RR1=fabs(W[5][0])+fabs(W[0][0]);
    RR2=fabs(W[0][nm1-1])+fabs(W[5][N-1]);
    if(RR1>=RR2) R=RR1;
    else         R=RR2;

    for(I=1;I<nm1;I++)
    {
      T=fabs(W[0][I-1])+fabs(W[5][I])+fabs(W[0][I]);
      if(T>R) R=T;
    }
    EPS1=R*1.0E-16;
    EPS2=R*EPS;
    for(I=0;I<nm1;I++) W[1][I]=W[0][I]*W[0][I]; /*{W1}={AiiAii}*/
    if(NE<0) R=-R;
    F=R; /*F:UPPER BOUND FOR DECENDING,LOWER FOR ACENDING.*/
    for(I=0;I<NEA;I++) E[I]=-R;

    for(K=0;K<NEA;K++)
    {
      D=E[K]; /*D:LOWER BOUND FOR DECENDING,UPPER FOR ACENDING.*/

      while(1)
      {
        T = 0.5*(D+F); /*DIVIDE SECTION.*/

        if(fabs(D-F)<=EPS2 || T==D || T==F) break;

        J=0;
        I=0;
        while(1)
        {
          Q=W[5][I]-T;

          while(1)
          {
            if(Q>=0.0) J=J+1;
            if(Q!=0.0)
            {
              I=I+1;
              if(I>=N) break;

              Q=W[5][I]-T-W[1][I-1]/Q;
            }
            else /*Q==0.0*/
            {
              I=I+2;
              break;
            }
          }
          if(I>=N) break;
        }
        if(NE<0) J=N-J;

        if(J<=K)
        {
          F=T; /*F:UPPER BOUND FOR DECENDING,LOWER FOR ACENDING.*/
        }
        else
        {
          D=T; /*D:LOWER BOUND FOR DECENDING,UPPER FOR ACENDING.*/
          if(J<=NEA) M=J; /*M=MIN{J,NEA}*/
          else       M=NEA;
          for(I=K;I<M;I++) E[I]=T;
        }
      }
      E[K]=T;
    }
  }

  /* COMPUTE EIGENVECTORS BY INVERSE ITERATION */
  if(NV==0) return;
  if(N==2)
  {
    W[5][0]=A[0][0];
    W[5][1]=A[1][1];
  }

  W[0][N-1]=0.0;
  MM=584287;

  for(I=0;I<NVA;I++)
  {
    for(J=0;J<N;J++)
    {
      W[1][J]=W[5][J]-E[I];
      W[2][J]=W[0][J];
      V[I][J] = 1.0;
    }
    SW=FALSE;

    /*REDUCE TO TRIANGULAR FORM*/
    for(J=0;J<nm1;J++)
    {
      if(fabs(W[1][J])>=fabs(W[0][J]))
      {
        if(W[1][J]==0.0) W[1][J]=1.0E-30;
        W[4][J]=W[0][J]/W[1][J];
        LW[J]=0.0; /*FALSE*/
        W[1][J+1]=W[1][J+1]-W[4][J]*W[2][J];
        W[3][J]=0.0;
      }
      else
      {
        W[4][J]=W[1][J]/W[0][J];
        LW[J]=1.0; /*TRUE*/
        W[1][J]=W[0][J];
        T=W[2][J];
        W[2][J]=W[1][J+1];
        W[3][J]=W[2][J+1];
        W[1][J+1]=T-W[4][J]*W[2][J];
        W[2][J+1]=-W[4][J]*W[3][J];
      }
    }
    if(W[1][N-1]==0.0) W[1][N-1]=1.0E-30;

    /*BEGIN BACK SUBSTITUTION*/
    if(I!=0 && fabs(E[I]-E[I-1])<EPS1)
    {
      /*GENERATE RANDOM NUMBERS*/
      for(J=0;J<N;J++)
      {
        MM=MM*48828125;
        V[I][J]=(double)MM*0.4656613E-9;
      }
    }

    while(1)
    {
      T=V[I][N-1];
      R=V[I][N-2];

      while(1)
      {
        V[I][N-1]=T/W[1][N-1];
        V[I][nm1-1]=(R-W[2][nm1-1]*V[I][N-1])/W[1][nm1-1];

        of=0; /*OVERFLOW FLAG.*/
        if(T>1.0E+05) of=1;
        if(R>1.0E+05) of=1;
        for(J=0;J<nm2;J++)
        {
          if(V[I][J]>1.0E+05) of=1;
        }

        if(of) /*IF POSITIVE OVERFLOW.*/
        {
          for(J=0;J<nm2;J++) V[I][J]*=1.0E-5;
          T=T*1.0E-5;
          R=R*1.0E-5;
        }
        else break;
      }

      if(N!=2)
      {
        K=nm2-1;
        while(1)
        {
          T=V[I][K];
          while(1)
          {
            V[I][K]=(T-W[2][K]*V[I][K+1]-W[3][K]*V[I][K+2])
                    /W[1][K];

            of=0; /*OVERFLOW FLAG.*/
            if(T>1.0E+5) of=1;
            for(J=0;J<N;J++)
            {
              if(V[I][J]>1.0E+5) of=1;
            }

            if(of) /*IF POSITIVE OVERFLOW.*/
            {
              for(J=0;J<N;J++) V[I][J]*=1.0E-5;
              T=T*1.0E-5;
            }
            else break;
          }
          K=K-1;
          if(K<0) break;
        }
      }

      if(SW) break;
      SW=TRUE;
      for(J=0;J<nm1;J++)
      {
        if(!(int)(LW[J]))
        {
          V[I][J+1]=V[I][J+1]-W[4][J]*V[I][J];
        }
        else
        {
          T=V[I][J];
          V[I][J]=V[I][J+1];
          V[I][J+1]=T-W[4][J]*V[I][J+1];
        }
      }
    }
  }

  /* BEGIN BACK TRANSFORMATION */
  if(N!=2)
  {
    for(I=0;I<nm2;I++) W[0][I]=-W[0][I]*A[I+1][I];
    for(I=0;I<NVA;I++)
    {
      K=nm2-1;
      while(1)
      {
        R=W[0][K];
        if(R!=0.0)
        {
          R = 1.0/R;
          S=0.0;
          K1=K+1;
          for(J=K1;J<N;J++)
          {
            S+=A[J][K]*V[I][J];
          }
          R=R*S;
          for(J=K1;J<N;J++)
          {
            V[I][J]-=R*A[J][K];
          }
        }
        K=K-1;
        if(K<0) break;
      }
    }
  }

  /*NORMALIZE EIGENVECTORS          */
  /*NORMALIZE AS MAXIMUM ELEMENT = 1*/
  for(I=0;I<NVA;I++)
  {
    T=fabs(V[I][0]);
    K=0;
    for(J=1;J<N;J++)
    {
      R=fabs(V[I][J]);
      if(T<R)
      {
        T=R;
        K=J;
      }
    }
    T = 1.0/V[I][K];
    for(J=0;J<N;J++) V[I][J]*=T;
  }
  if(NV<0) return;

  /*ORTHONORMALIZE AS NORM = 1*/
  for(I=0;I<NVA;I++)
  {
    if(I!=0 && fabs(E[I]-E[I-1])<EPS1)
    {
      /* ORTHONORMALIZE EIGENVECTORS FOR DEGENERATED EIGENVALUES */
      I1=I-1;
      for(J=M;J<I1;J++)
      {
        S=0.0;
        for(K=0;K<N;K++) S+=V[J][K]*V[I][K];
        for(K=0;K<N;K++) V[I][K]-=S*V[J][K];
      }
    }
    else
    {
      M=I;
    }

    /*NORMALIZE AS NORM = 1*/
    S=0.0;
    for(J=0;J<N;J++) S+=V[I][J]*V[I][J];
    T=0.0;
    if(S!=0.0) T = sqrt(1.0/S);
    for(J=0;J<N;J++) V[I][J]*=T;
  }

  return;
}/*deigrs*/

/*=================================================================*/
void currentvalue(char *str,
                  long int n,long int ne,
                  struct gcomponent *A,
                  double **W,
                  double *E,double **V)
/*CHECK CURRENT VALUES FOR DEBUG.*/
{
  char non[10];
  long int i,j;
  double data;

  ne=labs(ne);

  if(A!=NULL)
  {
    for(i=1;i<=n;i++)
    {
      fprintf(stderr,"A%d:",i);
      for(j=1;j<=i;j++)
      {
        gread(A,i,j,&data);
        fprintf(stderr," %11.5f",data);
      }
      fprintf(stderr,"\n");
    }
  }

  if(E!=NULL)
  {
    fprintf(stdout,"\n");
    for(i=1;i<=ne;i++) fprintf(stdout,"           E%ld",i);
    fprintf(stdout,"\n");
    for(j=0;j<ne;j++)
    {
      fprintf(stdout," %12.5f",*(E+j));
    }
    fprintf(stdout,"\n");
  }

  if(V!=NULL)
  {
    fprintf(stdout,"\n");
    for(i=1;i<=ne;i++) fprintf(stdout,"           V%ld",i);
    fprintf(stdout,"\n");
    for(i=0;i<n;i++)
    {
      for(j=0;j<ne;j++)
      {
        fprintf(stdout," %12.5f",*(*(V+j)+i));
      }
      fprintf(stdout,"\n");
    }
  }

  if(W!=NULL)
  {
    fprintf(stdout,"\n");
    for(i=0;i<6;i++)
    {
      fprintf(stdout,"W%d:",i);
      for(j=0;j<n;j++)
      {
        fprintf(stdout," %12.5E",*(*(W+i)+j));
      }
      fprintf(stdout,"\n");
    }
  }

  fprintf(stdout,"%s\n",str);
  gets(non);
  if(strcmp(non,"\0")) exit(1);
  return;
}/*currentvalue*/

/*-----------------------------------------------------------------*/
void deigabgeneral(struct gcomponent *A,
                   struct gcomponent *B,
                   struct confinement *confs,
                   long int N,long int NE,long int NV,
                   double EPS,
                   double *E,double **V)
{
  /* SUBROUTINE FOR GENERALIZED EIGENVALUE PROBLEM                */
  /*  [A]{X} = LAMDA [B]{X}                                       */
  /* FOR REAL SYMMETRIC MATRICES [A] & [B], THE LATTER BEING      */
  /* POSITIVE DEFINITE.                                           */
  /*                                                              */
  /* INPUT:                                                       */
  /*   A[N][N]  : 2-DIM. ARRAY CONTAINING REAL SYMMETRIC MATRIX.  */
  /*   B[N][N]  : 2-DIM. ARRAY CONTAINING REAL SYMMETRIC          */
  /*              POSITIVE DEFINITE MATRIX.                       */
  /*   N        : ORDER OF MATRIX.                                */
  /*              SIZE OF THE 2-DIM. ARRAYS  A, B, W, V           */
  /*              DEFINED IN 'DIMENSION' STATEMENT (FIRST INDEX). */
  /*   NE       : NUMBER OF EIGENVALUES TO BE OBTAINED.           */
  /*              IN DESCENDING ORDER  WHEN  NE > 0,              */
  /*              IN ASCENDING  ORDER  WHEN  NE < 0.              */
  /*   NV       : NUMBER OF EIGENVECTORS TO BE OBTAINED.          */
  /*   EPS      : ACCURACY  ( STANDARD VALUE  1.0D-16 ).          */
  /*                                                              */
  /* OUTPUT:                                                      */
  /*   E[NE]    : 1-DIM. ARRAY CONTAINING THE OUTPUT EIGENVALUES  */
  /*   V[NV][N] : 2-DIM. ARRAY CONTAINING THE OUTPUT EIGENVECTORS */
  /*              THE VECTOR ( V(1,K), V(2,K),..., V(N,K) )       */
  /*              BELONGS TO THE EIGENVALUE  E(K).                */
  /*                                                              */
  /* WORKING SPACE:                                               */
  /*   W[7][N]  : 2-DIM. ARRAY  (7,N)  USED AS THE WORKING AREA.  */
  /*                                                              */
  /* NOTE:                                                        */
  /*  THE MATRICES [A] AND [B] ARE DESTROYED.                     */
  /* METHOD:                                                      */
  /*  CHOLESKY REDUCTION FOLLOWED BY HOUSEHOLDER DIAGONALIZATION. */
  /* SUBROUTINE USED : DEIGRS                                     */

  /*char non[256];*/
  long int i,ii,kk;
  long int J,K,nev,neva,nvec,R;
  double T;

  unsigned short int mm,mmm;
  double lkj,lim;
  struct gcomponent *gi,*gj,*gr,*gk,*gp,*gg;
  struct gcomponent *bj,*bm;

  long int ipiv; /*PIVOT LINE.*/
  signed char ic; /*CONFINEMENT ID.*/
  long int ndim; /*MATRIX DIMENSIONS.*/

  /*CHECK INPUT DATA.*/
  nev  = NE;
  neva = labs(nev);
  nvec = NV;

  if(N<=0 || neva==0 || N-neva<0 || nvec<0 || neva-nvec<0)
  {
    fprintf(stdout,"DEIGAB:INVALID ARGUMENT.");
    fprintf(stdout," N,NE,NV = %ld,%ld,%ld",N,nev,nvec);
    return;
  }

  ndim=0;
  for(i=0;i<N;i++) /*COUNT DIMENSIONS.*/
  {
    if(!(confs+i)->iconf) ndim++;
  }
  if(ndim<=0) return;

  ipiv=0; /*SEARCH HEAD LINE.*/
  while((confs+ipiv)->iconf && ipiv<N) ipiv++;
  if(ipiv>=N) return;

  if(ndim==1)
  {
    T = (B+ipiv)->value;
    if(T<=0.0)
    {
      fprintf(stdout,"DEIGAB:MATRIX [B] IS NOT POSITIVE DEFINITE.");
    }
    else
    {
      (B+ipiv)->value = sqrt(1.0/T);
      *(E+0) = ((A+ipiv)->value)/T;
      *(*(V+0)+ipiv) = 1.0;
    }
    return;
  }

  /* CHOLESKI DECOMPOSITION OF THE POSITIVE DEFINITE */
  /* MATRIX [B] INTO A PRODUCT OF A LOWER TRIANGULAR */
  /* MATRIX [L] WITH ITS TRANSPOSED MATRIX.          */
  /* DIAGONALS OF [L] ARE REMAINED INVERSE.          */

  T = (B+ipiv)->value;
  if(T<=0.0)
  {
    fprintf(stdout,"DEIGAB:MATRIX [B] IS NOT POSITIVE DEFINITE.");
    return;
  }
  T = sqrt(1.0/T);
  (B+ipiv)->value = T;

  gi=(B+ipiv);
  while(gi->down != NULL) /*FIRST ROW.*/
  {
    gi = gi->down;
    ic=(confs+(gi->m)-1)->iconf;
    if(!ic) gi->value *= T;
  }

  for(R=(ipiv+1);R<N;R++)
  {
    gr=(B+R-1);
    ic=(confs+R-1)->iconf;
    if(!ic)
    {
      while(gr->down != NULL)
      {
        gr=gr->down;
        mm=gr->m;
        ic=(confs+mm-1)->iconf;
        if(!ic)
        {
          (B+mm-1)->value -= (gr->value)*(gr->value); /*DIAGONALS.*/

          gk=(B+mm-1);
          gi=gr;
          while(gi->down != NULL)
          {
            gi=gi->down;
            mmm=gi->m;
            ic=(confs+mmm-1)->iconf;
            if(!ic)
            {
              while((gk->m) < mmm  &&  gk->down!=NULL) /*SEARCH.*/
              {
                gp=gk;
                gk=gk->down;
              }

              if(gk->m ==mmm)
              {
                gk->value -= (gr->value)*(gi->value);
              }
              else if(gk->m < mmm)
              {
                gg=gdefine(mmm,mm,-(gr->value)*(gi->value),
                           NULL,NULL);
                gk->down=gg;
                gp=gk;
                gk=gg;
              }
              else if(gk->m > mmm)
              {
                gg=gdefine(mmm,mm,-(gr->value)*(gi->value),
                           gk,NULL);
                gp->down=gg;
                gk=gg;
              }
            }
          }
        }
      }
    }

    ic=(confs+R)->iconf;
    if(!ic)
    {
      if((B+R)->value <= 0.0)
      {
        fprintf(stdout,
                "DEIGAB:MATRIX [B] IS NOT POSITIVE DEFINITE.\n");
        fprintf(stdout,"T%ld=%9.5E\n",R,(B+R)->value);
        return;
      }
      T = sqrt(1.0/(B+R)->value);
      (B+R)->value = T;

      gi=(B+R);
      while(gi->down != NULL)
      {
        gi=gi->down;
        ic=(confs+(gi->m)-1)->iconf;
        if(!ic) gi->value *= T;
      }
    }
  }

  /*currentvalue("DEIGABGENERAL:[B] DECOMPOSED.",N,NE,B,NULL,NULL,NULL);*/

  /* CHOLESKI DECOMPOSITION IS COMPLETED.                */
  /* [B] = [L][LT],  THE UPPER TRIANGULAR MATRIX [LT] IS */
  /* STORED IN THE ARRAY [B].                            */
  /* NOW, PROCEED TO EVALUATE  [L]**(-1) [A] [LT]**(-1)  */

  for(J=ipiv;J<N;J++)
  {
    gj=(A+J); /*DIAGONAL.*/
    ic=(confs+J)->iconf;
    if(!ic)
    {
      while(1)
      {
        mmm=gj->m;
        ic=(confs+mmm-1)->iconf;
        if(!ic)
        {
          gj->value *= ((B+J)->value)*((B+mmm-1)->value);

          bj=(B+J);
          while(1)
          {
            kk=bj->m;
            ic=(confs+kk-1)->iconf;
            if(!ic)
            {
              if(kk==J+1) lkj=1/(bj->value); /*DIAGONAL.*/
              else        lkj=bj->value;

              bm=(B+mmm-1);
              if(kk==J+1 && (bm->m)<=mmm) bm=bm->down;

              if(bm!=NULL)
              {
                gk=(A+kk-1);

                while(1) /*FOR ROW M.*/
                {
                  ii=bm->m;
                  ic=(confs+ii-1)->iconf;
                  if(!ic)
                  {
                    if(ii==mmm) lim=1.0/(bm->value);
                    else        lim=bm->value;

                    if(mmm!=J+1 && ii!=J+1 && ii<=kk)
                    {
                      gi=(A+ii-1);

                      while((gi->m)<kk &&
                            gi->down!=NULL) /*SEEK LINE I.*/
                      {
                        gp=gi;
                        gi=gi->down;
                      }
                      if(gi->m == kk)
                      {
                        gi->value -= lkj*lim*(gj->value);
                      }
                      else if(gi->m < kk) /*ADD.*/
                      {
                        gg=gdefine((unsigned short int)kk,
                                   (unsigned short int)ii,
                                   -lkj*lim*(gj->value),NULL,NULL);
                        gi->down=gg;
                        gp=gi;
                        gi=gg;
                      }
                      else if(gi->m > kk) /*INSERT.*/
                      {
                        gg=gdefine((unsigned short int)kk,
                                   (unsigned short int)ii,
                                   -lkj*lim*(gj->value),gi,NULL);
                        gp->down=gg;
                        gp=gg;
                      }
                    }

                    if(ii>=kk)
                    {
                      while((gk->m)<ii &&
                            gk->down!=NULL) /*SEEK LINE M.*/
                      {
                        gp=gk;
                        gk=gk->down;
                      }
                      if(gk->m == ii)
                      {
                        gk->value -= lkj*lim*(gj->value);
                      }
                      else if(gk->m < ii) /*ADD.*/
                      {
                        gg=gdefine((unsigned short int)ii,
                                   (unsigned short int)kk,
                                   -lkj*lim*(gj->value),NULL,NULL);
                        gk->down=gg;
                        gp=gk;
                        gk=gg;
                      }
                      else if(gk->m > ii) /*INSERT.*/
                      {
                        gg=gdefine((unsigned short int)ii,
                                   (unsigned short int)kk,
                                   -lkj*lim*(gj->value),gk,NULL);
                        gp->down=gg;
                        gp=gg;
                      }
                    }
                  }
                  if(bm->down!=NULL) bm=bm->down;
                  else break;
                }
              }
            }
            if(bj->down!=NULL) bj=bj->down;
            else break;
          }
        }
        if(gj->down!=NULL) gj=gj->down;
        else break;
      }
    }
  }

  /*currentvalue("DEIGABGENERAL:TRANSFORMED.",N,NE,A,NULL,NULL,NULL);*/

  /*FIND EIGENVALUES AND EIGENVECTORS OF THE TRANSFORMED MATRIX.*/
  deigrsstandard(A,confs,N,nev,nvec,EPS,E,V);

  /*BACK TRANSFORMATION OF EIGENVECTORS.*/
  if(nvec==0) return;

  for(J=0;J<nvec;J++)
  {
    ipiv=N-1;
    while((confs+ipiv)->iconf && ipiv>=0) ipiv--;
    if(ipiv<0) return;

    *(*(V+J)+ipiv) *= ((B+ipiv)->value);
    for(K=(ipiv-1);K>=0;K--)
    {
      gk=(B+K);
      ic=(confs+K)->iconf;
      if(!ic)
      {
        while(gk->down!=NULL)
        {
          gk=gk->down;
          R=(gk->m)-1;
          ic=(confs+R)->iconf;
          if(!ic) *(*(V+J)+K) -= (gk->value)*(*(*(V+J)+R));
        }
        *(*(V+J)+K) *= ((B+K)->value);
      }
    }
  }
  currentvalue("DEIGABGENERAL END.",N,NE,A,NULL,E,V);

  return;
}/*deigabgeneral*/

/*-----------------------------------------------------------------*/
void deigrsstandard(struct gcomponent *A,
                    struct confinement *confs,
                    long int N,long int NE,long int NV,
                    double EPS,
                    double *E,double **V)
{
  /* SUBROUTINE FOR STANDARD EIGENVALUE PROBLEM                   */
  /* HOUSEHOLDER'S TRIDIAGONAL REDUCTION.                         */
  /* EIGENVALUES BY BISECTION.                                    */
  /* EIGENVECTORS BY INVERSE ITERATION.                           */
  /*                                                              */
  /*  [A]{V} = LAMBDA {V}                                         */
  /*                                                              */
  /* INPUT:                                                       */
  /*   A[N][N]  : 2-DIM. ARRAY CONTAINING REAL SYMMETRIC MATRIX.  */
  /*   N        : ORDER OF MATRIX.                                */
  /*              SIZE OF THE 2-DIM. ARRAYS  A, B, W, V           */
  /*              DEFINED IN 'DIMENSION' STATEMENT (FIRST INDEX). */
  /*   NE       : NUMBER OF EIGENVALUES TO BE OBTAINED.           */
  /*              IN DESCENDING ORDER  WHEN  NE > 0,              */
  /*              IN ASCENDING  ORDER  WHEN  NE < 0.              */
  /*   NV       : NUMBER OF EIGENVECTORS TO BE OBTAINED. NV<=|NE| */
  /*              ONLY EIGENVALUES IF NV=0.                       */
  /*   EPS      : ACCURACY  ( STANDARD VALUE  1.0D-16 ).          */
  /*                                                              */
  /*   W[6][N]  : WORK SPACE FOR TRIDIAGONALS,ITERATION etc.      */
  /*   LW[N]    :                                                 */
  /*                                                              */
  /* OUTPUT:                                                      */
  /*   E[NE]    : 1-DIM. ARRAY CONTAINING THE OUTPUT EIGENVALUES  */
  /*   V[NV][N] : 2-DIM. ARRAY CONTAINING THE OUTPUT EIGENVECTORS */
  /*              THE VECTOR { V(K,1), V(K,2),..., V(K,N) }       */
  /*              BELONGS TO THE EIGENVALUE  E(K).                */

  /*char non[256];*/
  int SW,of;
  long int I,I1,J,J1,K,K1,MM,NEA,NVA,nm1,nm2,M;
  double T,R,RR1,RR2,D,F,Q,S,SR,A1,EPS1,EPS2;

  long int i,j,k,k1,k2,ii,jj;
  double *ww;

  struct gcomponent *gi,*gj,*gk,*gp,*gg;
  struct gcomponent **gpp;

  int *LW;
  double **W;

  long int ipiv,ipiv1,mpiv; /*PIVOT LINE.*/
  signed char ic;
  long int ndim; /*DIMENSIONS.*/
  long int *lines;

  NEA=labs(NE);
  if(NEA==0)
  {
    fprintf(stdout,"DEIGRS: NE = %ld\n",NE);
    fprintf(stdout,
            "NE SHOULD BE NON-ZERO. RETURN WITH NO CALCULATION.\n");
    return;
  }
  NVA=labs(NV);
  if(NVA>NEA || NEA>N)
  {
    fprintf(stdout,"DEIGRS: NV,NE,N = %ld,%ld,%ld\n",NV,NE,N);
    fprintf(stdout,"NV,NE,N");
    fprintf(stdout," SHOULD SATISFY THE FOLLOWING INEQUALITIES.\n");
    fprintf(stdout,"|NV|<=|NE|<=N RETURN WITH NO CALCULATION.\n");
    return;
  }

  W=(double **)malloc(6*sizeof(double *));
  for(i=0;i<6;i++)
  {
    ww=(double *)malloc(N*sizeof(double));
    *(W+i)=ww;
  }
  LW=(int *)malloc(N*sizeof(int));
  gpp=(struct gcomponent **)malloc(N*sizeof(struct gcomponent *));

  ndim=0;
  for(i=0;i<N;i++) /*COUNT DIMENSIONS.*/
  {
    if(!(confs+i)->iconf) ndim++;
  }
  if(ndim<=0) return;
  if(NEA>ndim)
  {
    fprintf(stderr,"ORDERING TOO MANY EIGEN VALUES.\n");
    fprintf(stdout,"RETURN WITH NO CALCULATION.\n");
    return;
  }

  lines=(long int *)malloc(ndim*sizeof(long int));
  j=0;
  for(i=0;i<N;i++) /*PICK UP LINES.*/
  {
    if(!(confs+i)->iconf)
    {
      *(lines+j)=i;
      j++;
    }
  }

  ipiv = *(lines+0); /*HEAD LINE.*/
  if(ndim>=2) ipiv1 = *(lines+1); /*SECOND LINE.*/

  if(ndim>=2) nm2 = *(lines+ndim-2); /*SECOND TAIL LINE.*/
  nm1 = *(lines+ndim-1); /*TAIL LINE.*/

  if(EPS<0.0) EPS=1.0E-16;

  /* CASE N=1 : ONLY ONE NUMBER. */
  if(ndim==1)
  {
    *(E+0)=(A+ipiv)->value;
    if(NV!=0) *(*(V+0)+ipiv) = 1.0;
    return;
  }

  /* CASE N=2 : COMPUTE EIGENVALUES OF 2x2 MATRIX. */
  else if(ndim==2)
  {
    gi=(A+ipiv);
    while((gi->m)<(ipiv1+1) && gi->down!=NULL) gi=gi->down;

    if((gi->m)==(ipiv1+1))
    {
      *(*(W+0)+ipiv)=gi->value; /*W[0]=A[1][0]*/
    }
    else *(*(W+0)+ipiv)=0.0;

    T = 0.5*(((A+ipiv)->value)+((A+ipiv1)->value));
    R=((A+ipiv)->value)*((A+ipiv1)->value)
     -(*(*(W+0)+ipiv))*(*(*(W+0)+ipiv));
    D=T*T-R;
    Q=fabs(T)+sqrt(D);
    if(T<0.0) Q=-Q;
    T = T*(double)NE;
    if(T>=0.0)
    {
      *(E+0)=Q;
      if(NEA==2) *(E+1)=R/Q;
    }
    else
    {
      *(E+0)=R/Q;
      if(NEA==2) *(E+1)=Q;
    }
  }

  /* CASE N=3,4,... */
  /* REDUCE TO TRIDIAGONAL FORM BY HOUSEHOLDER'S METHOD */
  else if(ndim>=3)
  {
    for(K=ipiv;K<nm2;K++) /* N-2 TIMES REDUCTION. */
    {
      K1=K+1;

      gi=(A+K);
      ic=(confs+(gi->m)-1)->iconf;
      if(!ic)
      {
        S=0.0;
        while(gi->down!=NULL)
        {
          gi=gi->down;
          ic=(confs+(gi->m)-1)->iconf;
          if(!ic) S += (gi->value)*(gi->value); /* S={aT}{a} */
        }

        *(*(W+0)+K)=0.0;

        gi=(A+K);
        if(S!=0.0)
        {
          SR=sqrt(S); /* SR=SQRT( {aT}{a} ) */

          mpiv=K+1; /*SEARCH NEXT OF DIAGONAL LINE.*/
          while((confs+mpiv)->iconf) mpiv++;

          gj=gi; /*SEARCH NEXT OF DIAGONAL.*/
          while((gj->m)<(mpiv+1) && gj->down!=NULL)
          {
            gp=gj;
            gj=gj->down;
          }

          if(gj->m==mpiv+1) A1=gj->value;
          else              A1=0.0;

          if(A1<0.0) SR=-SR; /* SGN(SR)=SGN(A1) */

          *(*(W+0)+K)=-SR; /* a1=-SR */

          R = 1.0/(S+A1*SR);

          /*{V}={a(K+1) a(K+2)...a(N)}+{s 0 0...0}*/
          if((gj->m)==(mpiv+1))
          {
            gj->value=A1+SR;
          }
          else if((gj->m)<(mpiv+1)) /*ADD.*/
          {
            gg=gdefine((unsigned short int)(mpiv+1),
                       (unsigned short int)(K+1),
                       (A1+SR),
                       gj->down,NULL);
            gj->down=gg;
          }
          else if((gj->m)>(mpiv+1)) /*FILL.*/
          {
            gg=gdefine((unsigned short int)(mpiv+1),
                       (unsigned short int)(K+1),
                       (A1+SR),
                       gj,NULL);
            gp->down=gg;
          }

          gi=(A+K);
          for(J=K1;J<N;J++) /*INITIAL.*/
          {
            *(gpp+J)=(A+J);
            *(*(W+0)+J)=0.0;
          }
          while(gi->down!=NULL)
          {
            gi=gi->down; /*A[J][K]*/
            jj=gi->m;
            ic=(confs+jj-1)->iconf;
            if(!ic)
            {
              for(J=K1;J<jj;J++)
              {
                ic=(confs+J)->iconf;
                if(!ic)
                {
                  while(((*(gpp+J))->m)<jj &&
                        ((*(gpp+J))->down)!=NULL)
                  {
                    *(gpp+J)=(*(gpp+J))->down;
                  }

                  if(((*(gpp+J))->m)==jj) /*{Woi}=R{Si}*/
                  {
                    S=(gi->value)
                     *((*(gpp+J))->value); /*Si={VT}{Ai}*/

                     *(*(W+0)+J)+=R*S;
                  }
                }
              }

              gj=(A+jj-1);
              while(gj->down!=NULL)  /*Si={VT}{Ai}*/
              {
                gj=gj->down;
                J=(gj->m)-1;
                ic=(confs+J)->iconf;
                if(!ic) *(*(W+0)+J)+=R*(gi->value)*(gj->value);
              }
            }
          }

          S=0.0;
          gi=(A+K);
          while(gi->down!=NULL) /*S'={VT}{Wo}*/
          {
            gi=gi->down;
            I=(gi->m)-1;
            ic=(confs+I)->iconf;
            if(!ic) S += (*(*(W+0)+I))*(gi->value);
          }

          T = 0.5*R*S;
          gi=(A+K);
          while(gi->down!=NULL) /*{Wo'}={Wo}-0.5RS'{V}*/
          {
            gi=gi->down;
            I=(gi->m)-1;
            ic=(confs+I)->iconf;
            if(!ic) *(*(W+0)+I) -= T*(gi->value);
          }

          gk=(A+K);
          while(gk->down!=NULL) /*[Aij']=[Aij-Wj'Vi-Wi'Vj]*/
          {
            gk=gk->down;
            J=(gk->m)-1; /*Ajk*/
            ic=(confs+J)->iconf;
            if(!ic)
            {
              gi=gk;
              gj=(A+J);
              while(1)
              {
                ii=(gi->m); /*Aik*/
                I=(gi->m)-1;
                ic=(confs+I)->iconf;
                if(!ic)
                {
                  while((gj->m)<ii && gj->down!=NULL)
                  {
                    gp=gj;
                    gj=gj->down;
                  }

                  if(gj->m == ii)
                  {
                    gj->value-=(*(*(W+0)+J))*(gi->value)  /*WjAik*/
                              +(*(*(W+0)+I))*(gk->value); /*WiAjk*/
                  }
                  else if(gj->m < ii) /*ADD.*/
                  {
                    gg=gdefine((unsigned short int)ii,
                               (unsigned short int)(J+1),
                               -((*(*(W+0)+J))*(gi->value)
                                +(*(*(W+0)+I))*(gk->value)),
                               NULL,NULL);
                    gj->down=gg;
                    gp=gj;
                    gj=gg;
                  }
                  else if(gj->m > ii) /*FILL IN.*/
                  {
                    gg=gdefine((unsigned short int)ii,
                               (unsigned short int)(J+1),
                               -((*(*(W+0)+J))*(gi->value)
                                +(*(*(W+0)+I))*(gk->value)),
                               gj,NULL);
                    gp->down=gg;
                    gp=gg;
                  }
                }
                if(gi->down!=NULL) gi=gi->down;
                else break;
              }
            }
          }
        }
      }
    }

    gi=(A+nm2);
    while((gi->m)<(nm1+1) && gi->down!=NULL)
    {
      gi=gi->down;
    }

    if((gi->m)==(nm1+1)) *(*(W+0)+nm2)=gi->value;
    else                 *(*(W+0)+nm2)=0.0;

/*currentvalue("DEIGRSSTANDARD:REDUCED.",N,NE,A,NULL,NULL,NULL);*/

    /*COMPUTE EIGENVALUES BY BISECTION METHOD*/

    for(i=0;i<ndim;i++)
    {
      I=*(lines+i);
      *(*(W+5)+I)=(A+I)->value; /*{W5}={Aii}*/
    }

    RR1=fabs(*(*(W+5)+ipiv))+fabs(*(*(W+0)+ipiv));
    RR2=fabs(*(*(W+0)+nm2))+fabs(*(*(W+5)+nm1));
    if(RR1>=RR2) R=RR1;
    else         R=RR2;

    for(I=1;I<(ndim-1);I++)
    {
      i=*(lines+I);
      j=*(lines+I-1);

      T=fabs(*(*(W+0)+j))+fabs(*(*(W+5)+i))+fabs(*(*(W+0)+i));
      if(T>R) R=T;
    }

    EPS1=R*1.0E-16;
    EPS2=R*EPS;

    for(I=0;I<(ndim-1);I++)
    {
      i=*(lines+I);
      *(*(W+1)+i)=(*(*(W+0)+i))*(*(*(W+0)+i)); /*{W1}={AiiAii}*/
    }
    if(NE<0) R=-R;
    F=R; /*F:UPPER BOUND FOR DECENDING,LOWER FOR ACENDING.*/
    for(I=0;I<NEA;I++) *(E+I)=-R;

    for(K=0;K<NEA;K++)
    {
      D=*(E+K); /*D:LOWER BOUND FOR DECENDING,UPPER FOR ACENDING.*/

      while(1)
      {
        T = 0.5*(D+F); /*DIVIDE SECTION.*/

        if(fabs(D-F)<=EPS2 || T==D || T==F) break;

        J=0;
        I=0;
        while(1)
        {
          i=*(lines+I);
          Q=(*(*(W+5)+i))-T;

          while(1)
          {
            if(Q>=0.0) J=J+1;
            if(Q!=0.0)
            {
              I=I+1;
              if(I>=ndim) break;

              i=*(lines+I);
              j=*(lines+I-1);
              Q=(*(*(W+5)+i))-T-(*(*(W+1)+j))/Q;
            }
            else /*Q==0.0*/
            {
              I=I+2;
              break;
            }
          }
          i=*(lines+I);
          if(i>=ndim) break;
        }
        if(NE<0) J=ndim-J;

        if(J<=K)
        {
          F=T; /*F:UPPER BOUND FOR DECENDING,LOWER FOR ACENDING.*/
        }
        else
        {
          D=T; /*D:LOWER BOUND FOR DECENDING,UPPER FOR ACENDING.*/
          if(J<=NEA) M=J; /*M=MIN{J,NEA}*/
          else       M=NEA;
          for(I=K;I<M;I++) *(E+I)=T;
        }
      }
      *(E+K)=T;
    }
  }

  /* COMPUTE EIGENVECTORS BY INVERSE ITERATION */
  if(NV==0) return;
  if(N==2)
  {
    *(*(W+5)+ipiv)=(A+ipiv)->value;
    *(*(W+5)+ipiv1)=(A+ipiv1)->value;
  }

  *(*(W+0)+nm1)=0.0;
  MM=584287;

  for(I=0;I<NVA;I++)
  {
    for(j=0;j<ndim;j++)
    {
      J=*(lines+j);
      *(*(W+1)+J)=*(*(W+5)+J)-*(E+I);
      *(*(W+2)+J)=*(*(W+0)+J);
      *(*(V+I)+J)=1.0;
    }
    SW=FALSE;

    /*REDUCE TO TRIANGULAR FORM*/
    for(j=0;j<(ndim-1);j++)
    {
      J=*(lines+j);
      J1=*(lines+j+1);

      if(fabs(*(*(W+1)+J))>=fabs(*(*(W+0)+J)))
      {
        if((*(*(W+1)+J))==0.0) *(*(W+1)+J)=(1.0E-30);
        *(*(W+4)+J)=(*(*(W+0)+J))/(*(*(W+1)+J));
        *(LW+J)=0; /*FALSE*/
        *(*(W+1)+J1) -= (*(*(W+4)+J))*(*(*(W+2)+J));
        *(*(W+3)+J)=0.0;
      }
      else
      {
        *(*(W+4)+J)=(*(*(W+1)+J))/(*(*(W+0)+J));
        *(LW+J)=1; /*TRUE*/
        *(*(W+1)+J)=(*(*(W+0)+J));
        T=(*(*(W+2)+J));
        *(*(W+2)+J)=(*(*(W+1)+J1));
        *(*(W+3)+J)=(*(*(W+2)+J1));
        *(*(W+1)+J1)=T-(*(*(W+4)+J))*(*(*(W+2)+J));
        *(*(W+2)+J1)=-(*(*(W+4)+J))*(*(*(W+3)+J));
      }
    }
    if(*(*(W+1)+nm1)==0.0) *(*(W+1)+nm1)=(1.0E-30);

    /*BEGIN BACK SUBSTITUTION*/
    if(I!=0 && fabs((*(E+I))-(*(E+I-1)))<EPS1)
    {
      /*GENERATE RANDOM NUMBERS*/
      for(j=0;j<ndim;j++)
      {
        J=*(lines+j);
        MM=MM*48828125;
        *(*(V+I)+J)=(double)MM*(0.4656613E-9);
      }
    }

    while(1)
    {
      T=(*(*(V+I)+nm1));
      R=(*(*(V+I)+nm2));

      while(1)
      {
        *(*(V+I)+nm1)=T/(*(*(W+1)+nm1));
        *(*(V+I)+nm2)=(R - (*(*(W+2)+nm2)) * (*(*(V+I)+nm1)))
                       /(*(*(W+1)+nm2));

        of=0; /*OVERFLOW FLAG.*/
        if(T>(1.0E+05)) of=1;
        if(R>(1.0E+05)) of=1;
        for(j=0;j<(ndim-2);j++)
        {
          J=*(lines+j);
          if((*(*(V+I)+J))>(1.0E+05)) of=1;
        }

        if(of) /*IF POSITIVE OVERFLOW.*/
        {
          for(j=0;j<(ndim-2);j++)
          {
            J=*(lines+j);
            *(*(V+I)+J) *= (1.0E-5);
          }
          T=T*(1.0E-5);
          R=R*(1.0E-5);
        }
        else break;
      }

      if(ndim!=2)
      {
        k=ndim-3;
        while(1)
        {
          K=*(lines+k);
          k1=*(lines+k+1);
          k2=*(lines+k+2);

          T=*(*(V+I)+K);
          while(1)
          {
            *(*(V+I)+K)=(T-(*(*(W+2)+K))*(*(*(V+I)+k1))
                          -(*(*(W+3)+K))*(*(*(V+I)+k2)))
                       /(*(*(W+1)+K));

            of=0; /*OVERFLOW FLAG.*/
            if(T>(1.0E+5)) of=1;
            for(j=0;j<ndim;j++)
            {
              J=*(lines+j);
              if((*(*(V+I)+J))>(1.0E+5)) of=1;
            }

            if(of) /*IF POSITIVE OVERFLOW.*/
            {
              for(j=0;j<ndim;j++)
              {
                J=*(lines+j);
                *(*(V+I)+J) *= (1.0E-5);
              }
              T=T*(1.0E-5);
            }
            else break;
          }
          k=k-1;
          if(k<0) break;
        }
      }

      if(SW) break;
      SW=TRUE;
      for(j=0;j<ndim-1;j++)
      {
        J=*(lines+j);
        J1=*(lines+j+1);
        if(!(int)(*(LW+J)))
        {
          *(*(V+I)+J1) -= (*(*(W+4)+J))*(*(*(V+I)+J));
        }
        else
        {
          T=(*(*(V+I)+J));
          *(*(V+I)+J)=(*(*(V+I)+J1));
          *(*(V+I)+J1)=T-(*(*(W+4)+J))*(*(*(V+I)+J1));
        }
      }
    }
  }

  /*BEGIN BACK TRANSFORMATION*/
  if(ndim!=2)
  {
    for(i=0;i<(ndim-2);i++)
    {
      I=*(lines+i);
      I1=*(lines+i+1);

      gi=(A+I);
      while((gi->m)<(I1+1) && gi->down!=NULL) gi=gi->down;

      if((gi->m)==(I1+1)) *(*(W+0)+I) *= -(gi->value);
      else                *(*(W+0)+I) = 0.0;
    }

    for(I=0;I<NVA;I++)
    {
      k=ndim-3;
      while(1)
      {
        K=*(lines+k);

        R=(*(*(W+0)+K));
        if(R!=0.0)
        {
          R=1.0/R;
          S=0.0;

          gk=(A+K);
          while(gk->down!=NULL)
          {
            gk=gk->down;
            J=(gk->m)-1;
            ic=(confs+J)->iconf;
            if(!ic) S += (gk->value)*(*(*(V+I)+J));
          }

          R=R*S;
          gk=(A+K);
          while(gk->down!=NULL)
          {
            gk=gk->down;
            J=(gk->m)-1;
            ic=(confs+J)->iconf;
            if(!ic) *(*(V+I)+J) -= R*(gk->value);
          }
        }
        k=k-1;
        if(k<0) break;
      }
    }
  }

  /*NORMALIZE EIGENVECTORS          */
  /*NORMALIZE AS MAXIMUM ELEMENT = 1*/
  for(I=0;I<NVA;I++)
  {
    T=fabs(*(*(V+I)+ipiv));
    K=ipiv;
    for(j=1;j<ndim;j++)
    {
      J=*(lines+j);

      R=fabs(*(*(V+I)+J));
      if(T<R)
      {
        T=R;
        K=J;
      }
    }
    T = 1.0/(*(*(V+I)+K));
    for(j=0;j<ndim;j++)
    {
      J=*(lines+j);

      *(*(V+I)+J) *= T;
    }
  }
  if(NV<0) return;

  /*ORTHONORMALIZE AS NORM = 1*/
  for(I=0;I<NVA;I++)
  {
    if(I!=0 && fabs(*(E+I)-*(E+I-1))<EPS1)
    {
      /* ORTHONORMALIZE EIGENVECTORS FOR DEGENERATED EIGENVALUES */
      I1=I-1;
      for(J=M;J<I1;J++)
      {
        S=0.0;
        for(k=0;k<ndim;k++)
        {
          K=*(lines+k);
          S+=(*(*(V+J)+K))*(*(*(V+I)+K));
        }
        for(k=0;k<ndim;k++)
        {
          K=*(lines+k);
          *(*(V+I)+K) -= S*(*(*(V+J)+K));
        }
      }
    }
    else
    {
      M=I;
    }

    /* NORMALIZE AS NORM = 1 */
    S=0.0;
    for(j=0;j<ndim;j++)
    {
      J=*(lines+j);
      S+=(*(*(V+I)+J))*(*(*(V+I)+J));
    }
    T=0.0;
    if(S!=0.0) T = sqrt(1.0/S);
    for(j=0;j<ndim;j++)
    {
      J=*(lines+j);
      *(*(V+I)+J) *= T;
    }
  }

  free(lines);
  free(LW);
  free(gpp);
  for(i=0;i<6;i++) free(*(W+i));
  free(W);

  return;
}/*deigrsstandard*/

/*-----------------------------------------------------------------*/
struct gcomponent *gdefine(unsigned short int m,
                           unsigned short int n,
                           double value,
                           struct gcomponent *down,
                           struct gcomponent *left)
/*DEFINE GCOMP PARAMETERS.*/
{
  struct gcomponent *g;

  g=(struct gcomponent *)malloc(sizeof(struct gcomponent));
  if(g==NULL) return NULL;

  g->m=m;
  g->n=n;
  g->value=value;
  g->down=down;
  g->left=left;

  return g;
}/*gdefine*/

int gread(struct gcomponent *gmtx,
          long int i,long int j,double *data)
/*READ DATA FROM FILE GLOBAL MATRIX.*/
{
  long int k;
  struct gcomponent *g;

  if(i<j) {k=i; i=j; j=k;}                       /*MATRIX SYMMETRIC*/

  g=(gmtx+(j-1));                                        /*DIAGONAL*/

  while((g->m)<i && (g->down)!=NULL) g=g->down;

  if((g->m)==i) *data=g->value; /*EXISTENT.*/
  else          *data=0.0;      /*EMPTY.*/

  return 1;
}/*gread*/

int gwrite(struct gcomponent *gmtx,
           long int i,long int j,double data)
/*WRITE DATA INTO FILE GLOBAL MATRIX.*/
{
  int n=0;
  long int k;
  struct gcomponent *gcomp,*gdown,*gleft,*pdown,*pleft;

  if(i<j) {k=i; i=j; j=k;}                       /*MATRIX SYMMETRIC*/

  gdown=(gmtx+(j-1));                                      /*DIAGONAL*/
  while((gdown->m)<i && (gdown->down)!=NULL)
  {
    pdown=gdown;
    gdown=gdown->down;
  }

  gleft=(gmtx+(i-1));                                      /*DIAGONAL*/
  while((gleft->n)>j && (gleft->left)!=NULL)
  {
    pleft=gleft;
    gleft=gleft->left;
  }

  if(gdown->m==i /*&& gleft->n==j*/)          /*COMPONENT EXISTENT.*/
  {
    if(data==0.0 && i!=j)     /*COMPONENT VANISHED.EXCEPT DIAGONAL.*/
    {
      pdown->down=gdown->down;                    /*"gdown"="gleft"*/
      pleft->left=gleft->left;
      free(gdown);

      n=1;
    }
    else
    {
      gdown->value=data;

      n=2;
    }
  }

  else if(data!=0.0)                             /*COMPONENT EMPTY.*/
  {
    gcomp=(struct gcomponent *)malloc(sizeof(struct gcomponent));
    if(gcomp==NULL) return 0;

    gcomp->m=(unsigned short int)i;
    gcomp->n=(unsigned short int)j;
    gcomp->value=data;

    if((gdown->m)<i)                                  /*ADD TO ROW.*/
    {
      gcomp->down=NULL;
      gdown->down=gcomp;

      n=2;
    }
    else if((gdown->m)>i)                          /*FILL INTO ROW.*/
    {
      gcomp->down=pdown->down;
      pdown->down=gcomp;

      n=4;
    }

    if((gleft->n)>j)                                 /*ADD TO LINE.*/
    {
      gcomp->left=NULL;
      gleft->left=gcomp;

      n++;
    }
    else if((gleft->n)<j)                         /*FILL INTO LINE.*/
    {
      gcomp->left=pleft->left;
      pleft->left=gcomp;

      n+=2;
    }
  }

  return n;
}/*gwrite*/

struct gcomponent *copygcompmatrix(struct gcomponent *gmtx,
                                   long int msize)
/*CREATE COPY OF GCOMP MATRIX.*/
{
  long int k,ncomp;
  struct gcomponent *gcpy,*go,*gc;

  gcpy=(struct gcomponent *)malloc(msize*sizeof(struct gcomponent));

  for(k=0;k<msize;k++)
  {
    (gcpy+k)->m=(gmtx+k)->m; /*DIAGONAL.*/
    (gcpy+k)->n=(gmtx+k)->n;
    (gcpy+k)->value=(gmtx+k)->value;

    go=(gmtx+k);
    ncomp=0;
    while(go->down!=NULL) /*COUNT COMPS IN COLUMN.*/
    {
      go=go->down;
      ncomp++;
    }

    gc=(struct gcomponent *)malloc(ncomp*sizeof(struct òˇ?ı*ë ˝¸ã • û	
˛\Rˇ« J £ R˝Áˇ9 ôˇÀ˛ò˛' ÿ˛3ı¡˛J¸Z˝y˚1˘€Pˇ	 , ª˝E ' ¢ F†ˇI˜ˇ˙uˆóˇŒˇtÈˇí ˇÍˇŒˇÄˇ  t )ˇ. Ÿ Ã ô ∑ˇdˇúˇ·ˇH
Bˇ.˜« #˘¬ 5 & •ÈB © `Ÿˇjˇˇˇ2˛) ûˇíˇÿˇ}ÈBˇbˇú˙Õ ˇn˛Q UˇóÉˇ¿˛8˛√˛È∫Ã˛ïTˇíƒ 7˛*˚@ˇª˝¡ˇˇô ∫˚Àˇ”ˇ™ˇ @ˇúı4	ô ƒˇ! ∆ «ˇ≈XÙπ˛∫˛–ˇ= ] ﬁˇYˇ…˛à\ü9 –≥ˇÓ (  Ù‹ˇ ˇq˝ˇÍˇ«˝∏˝Óˇ›ˇ∏ˇ‰ˇP Œ˛” º kˇÙ1t_ÆˇHˇ±ˇøˇÄ
n
 ˇØˇVˇ ®ˇz ãˇ! 3 „ˇèˇcÙ. "ÛZ  Lˇˇ¬ˇ+ Áˇª Ñ˛„ˇßˇä †˛Á H
û …Òˇ^k *¸ˇ†˚à Uˇ÷ˇ-
Fˇn ß ât„ˇO˝Q ç9˛£,  ¸˛  ¶ˇ»˝ÖˇQˇË ⁄ˇD	o 8øˇKd Ë ¨˝˛	∑˚⁄¢ˇ˜bˇ•ı ˇ Àˇ≠ˇ˚ˇÃˇÕ òˇ⁄˛Ÿ <ˇ€ˇ˘ˇ  ‰ˇ”ˇ<  ) o é ]„ÿˇ8ˇ Zˇ˚—ˇˇ∞˝§˝! S 4ˇï˚ZˇD«˙∏ 2 Æÿ  Hˇâ 5   W Âˇó d 
 —ˇ vˇ _Ùà˛: ® jê q˛Ÿ¸ M˘ rˇ≥ˇŒ˛ˆ¢Æˇêè˛< ‡Y˛Ó˝© / Î % Q ™ˇ‰ˇ»ˇ≈ˇ3ˇÇˇ ﬂâ Êˇg˛Ìˇj ]ˇ¥ˇ«ˇÎ  É Øˇ;ˇ¬˛6äÙÃ…˙# Ã˛£ˇÄˇ±?ˇÈ˙¢ ﬁˇ©ˇcˇá …˛£˛^ ˙˛  %ˇ9  ¸ˇ
 Ôˇ& ¿˛0 ú P ¯0 -  ¸ „ )˝yˇÎˇP˝3|5sﬁ˝Ì ¡ˆﬂˇV©ˇ∫7ˇ–Â˝‡˛¸ç |˛∏ˇÏˇ≈ˇ¯Àˇõˇ– ˇΩ¸ˇ 6^ ñ˛]ˇÄˇ Íˇoˇ§ˇˇBÙ˝∂
Ï˛Œ Uˆ?s˝™˛Ñ K Úˇ6 ™ B V˚˛≈µˇ¿zˇ  ô´ zˇm 99® ú˛5 ÎΩˇ+ x Ñˇâù	Ôˇˇ£˛YˇbˇÌˇT $˝4¸W ˚`ˇ ç  6≈˛ùˇÄD ] È¯? ™ˇy €˛“ s ? Rˇh»¯Ñˇ˚*π ı˛ºóˇ^ˇˇN˛^˝¶ˇH˜& ® ˚˛!F Ç Øˇ 8 9 D Ïˇnˇ hˇÔˇüˇ}ˇX
G ·Ú( Î˛9U Êˇ) 6ˇL  ∞ˇ`ˇöˇÔˇõ Ω VÚB˛ôˇˇ3ˇ¬˝|ˇ© Œ ô ˝Nˇ£ˇgÙ ˇlˇã F ª˛u ˇÛ	 + ıˇÑ˛Bˇ∆˛<Ù2 ∂˛Êˇ∆ô  ÅÙ: Ù»ˇ ( 	 A€ˇ3˛Íˇä˛9 5ˇ ÒˇÁˇˇ€ˇ¯ˇ◊ˇåˇF ≈ˇM »   ’ˇv ∏ˇΩˇh ˙ˇN ´  ôˇb ~Ì4ˇ  ∏ˇ! ¬ˇ/ ù Ï mˇ`˛*˝B ∆¢ ˛ˇ?	vˇGˇ÷1F π Íˇ3ˇ% ¡ˇ±˛s˛+ 
 cÊêˇˇj ˇˇ