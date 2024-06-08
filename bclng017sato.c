/*=================================================================*/
/*     SUBROUTINE BCLNG001                                         */
/*     ELASTIC BUCKLING FOR 3D FRAME                               */
/*     CODED BY JUN SATO                                           */
/*     DATE 1993.7.28 ... 1994.1.28                                */
/*     BASED ON 'Y6FRAMEB' BY YOSHINOBU FUJITANI                   */
/*=================================================================*/
/*                                                                 */
/* DEIGAB:GENERALIZED EIGENVALUE PROBLEM.                          */
/* DEIGRS:STANDARD EIGENVALUE PROBLEM.                             */
/* LAST DEBUG:1997.11.29.                                          */
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

#define MSIZE  24 /*MAX MATRIX SIZE BY ARRAY.*/
#define KSIZE  12 /*MATRIX SIZE FOR CONDENSATION.*/
#define NEIGEN  1 /*NUMBERS OF EIGEN VALUES TO BE OBTAINED.*/

int bclng001(struct arclmframe *af);
int bclng002(struct arclmframe *af); /*BUCKLING CONDENSATION*/
void currentvalues(char *str,
				   long int n,long int ne,
				   double A[][MSIZE],
				   double W[][MSIZE],
				   double E[],double V[][MSIZE]);
void deigqr(double A[][KSIZE],double B[][KSIZE],
			long int N,long int NSIZE,long int NE,long int NV,
			double EPS,double W[][KSIZE],
			double E[],double V[][KSIZE],
			signed char CF[]);
void deigqrcf(double A[][KSIZE],double B[][KSIZE],
			  long int N,long int NSIZE,long int NE,long int NV,
			  double EPS,double W[][KSIZE],
			  double E[],double V[][KSIZE],
			signed char CF[]);
void deigab(double A[][MSIZE],double B[][MSIZE],
			long int N,long int NSIZE,long int NE,long int NV,
			double EPS,double W[][MSIZE],
			double E[],double V[][MSIZE]);
void deigrs(double A[][MSIZE],
			long int N,long int N1,long int NE,long int NV,
			double EPS,double W[][MSIZE],double LW[],
			double E[],double V[][MSIZE]);

void currentvalue(char *string,
				  long int n,long int ne,
				  struct gcomponent *A,
				  double **W,
				  double *E,double **V);
void deigabgeneral(struct gcomponent *A,
				   struct gcomponent *B,
				   struct oconf *confs,
				   long int N,long int NE,long int NV,
				   double EPS,
				   double *E,double **V);
void deigrsstandard(struct gcomponent *A,
					struct oconf *confs,
					long int N,long int NE,long int NV,
					double EPS,
					double *E,double **V);

struct gcomponent *gdefine(unsigned short int m,
						   unsigned short int n,
						   double value,
						   struct gcomponent *down,
						   struct gcomponent *left);
struct gcomponent *copygcompmatrix(struct gcomponent *gmtx,
								   long int msize);
double eigensubstitution(struct gcomponent *amtx,
						 struct gcomponent *bmtx,
						 struct oconf *confs,
						 long int msize,
						 double eig,double *vct);
void outputmode(double *gvct,FILE *fout,int nnode,
				struct onode *nodes);
void updatemode(struct arclmframe *af,double *gvct);

/*EXTERNAL PARAMETERS*/
extern FILE *globalfile; /*GLOBAL FILE.*/

int bclng001(struct arclmframe *af)
/*ELASTIC BUCKLING FOR ARCLM FRAME.*/
/*AFTER "ARCLM001".*/
{
  DWORD memory0,memory1;

  FILE /**fin,*/*fout;                               /*FILE 8 BYTES*/
  /*FILE *felem,*fdisp,*freact;*/
  /*char dir[]=DIRECTORY;*/                        /*DATA DIRECTORY*/
  char s[80],string[1024];
  int i,j,ii,jj;
  int nnode,nelem,nsect/*,nreact*/;
  long int /*fsize,*/loff,msize;
  /*long int time;*/
  struct gcomponent ginit={0,0,0.0,NULL};
  struct gcomponent *kmtx,*gmtx;                    /*GLOBAL MATRIX*/
  double **gvct;                                    /*GLOBAL VECTOR*/
  double **drccos,**tmatrix,**estiff,*estress;
  /*double determinant,data;*/
  double gdata;
  clock_t t0/*,t1,t2*/;

  /*struct osect *sects;*/
  struct onode /**nodes,*/*ninit;
  struct owire elem;
  struct owire *elems;
  struct oconf *confs;

  long int neig;
  double eps=1.0E-16,*eigen;

  memory0=availablephysicalmemory("INITIAL:");   /*MEMORY AVAILABLE*/

  fout=fgetstofopen("\0","w",ID_OUTPUTFILE);          /*OUTPUT FILE*/
  globalfile=fout;

  t0=clock();                                        /*CLOCK BEGIN.*/

  nnode=af->nnode;                                 /*INPUT INITIAL.*/
  nelem=af->nelem;
  nsect=af->nsect;
  /*nreact=af->nreact;*/
  sprintf(string,"NODES=%d ELEMS=%d SECTS=%d",nnode,nelem,nsect);
  errormessage(string);
  fprintf(fout,"%s\n",string);

  msize=6*nnode;                           /*SIZE OF GLOBAL MATRIX.*/
  neig=NEIGEN;

  kmtx=(struct gcomponent *)          /*DIAGONALS OF GLOBAL MATRIX.*/
		malloc(msize*sizeof(struct gcomponent));
  gmtx=(struct gcomponent *)
		malloc(msize*sizeof(struct gcomponent));
  if(kmtx==NULL || gmtx==NULL) return 0;
  for(i=0;i<msize;i++)
  {
	(kmtx+i)->down=NULL;            /*GLOBAL MATRIX INITIALIZATION.*/
    (gmtx+i)->down=NULL;
  }

  gvct=(double **)malloc(neig*sizeof(double *));   /*GLOBAL VECTORS*/
  eigen=(double *)malloc(neig*sizeof(double));       /*EIGEN VALUES*/
  if(gvct==NULL || eigen==NULL) return 0;
  for(i=0;i<neig;i++)
  {
    *(gvct+i)=(double *)malloc(msize*sizeof(double));
    for(j=0;j<msize;j++) *(*(gvct+i)+j)=0.0;
  }

  /*fdisp=af->fdisp;*/                 /*DISPLACEMENT:6 DIRECTIONS.*/
  /*felem=af->felem;*/              /*CODE,12 BOUNDARIES,12 STRESS.*/
  /*freact=af->freact;*/                           /*REACTION FILE.*/

  /*sects=af->sects;*/
  /*nodes=af->nodes;*/
  ninit=af->ninit;
  elems=af->elems;
  confs=af->confs;

  GetAsyncKeyState(VK_LBUTTON);                   /*CLEAR KEY LEFT.*/
  GetAsyncKeyState(VK_RBUTTON);                  /*CLEAR KEY RIGHT.*/

  errormessage("BCLNG001:BUCKLING LINEAR.");
  availablephysicalmemory("REMAIN:");            /*MEMORY AVAILABLE*/

  for(i=1;i<=msize;i++)             /*GLOBAL MATRIX INITIALIZATION.*/
  {
	ginit.m=(unsigned short int)i;
    /*ginit.n=(unsigned short int)i;*/
	*(kmtx+(i-1))=ginit;
    *(gmtx+(i-1))=ginit;
  }
  /*comps=msize;*/ /*INITIAL COMPONENTS=DIAGONALS.*/

  laptime("ASSEMBLING GLOBAL MATRIX.",t0);

  estress=(double *)malloc(12*sizeof(double));

  for(i=1;i<=nelem;i++)                 /*ASSEMBLAGE GLOBAL MATRIX.*/
  {
    /*elem=*(elems+i-1);*/                     /*READ ELEMENT DATA.*/
    inputelem(elems,af->melem,i-1,&elem);

    for(ii=0;ii<=1;ii++) /*INITIALIZE COORDINATION.*/
    {
	  loff=elem.node[ii]->loff;
      for(jj=0;jj<3;jj++)
      {
        elem.node[ii]->d[jj]=(ninit+loff)->d[jj];
      }
    }

	drccos=directioncosine(elem.node[0]->d[0],
                           elem.node[0]->d[1],
                           elem.node[0]->d[2],
                           elem.node[1]->d[0],
                           elem.node[1]->d[1],
						   elem.node[1]->d[2],
                           elem.cangle);                 /*[DRCCOS]*/

    tmatrix=transmatrix(drccos);           /*TRANSFORMATION MATRIX.*/
    estiff=assememtx(elem);            /*ELASTIC MATRIX OF ELEMENT.*/
    estiff=modifyhinge(elem,estiff);               /*MODIFY MATRIX.*/
    estiff=transformation(estiff,tmatrix);         /*[K]=[Tt][k][T]*/
/*
for(ii=0;ii<12;ii++)
{
  sprintf(string,"\0");
  for(jj=0;jj<=ii;jj++)
  {
    sprintf(s," %12.5E",*(*(estiff+ii)+jj));
	strcat(string,s);
  }
  errormessage(string);
}
*/
    assemgstiffness(kmtx,estiff,&elem);       /*ASSEMBLAGE ELASTIC.*/

    for(ii=0;ii<=11;ii++) free(*(estiff+ii));
    free(estiff);

    for(ii=0;ii<=1;ii++)                                /*STRESSES.*/
    {
      for(jj=0;jj<6;jj++) *(estress+6*ii+jj)=elem.stress[ii][jj];
    }
/*
sprintf(string,"\0");
for(ii=0;ii<12;ii++)
{
  sprintf(s," %12.5E",*(estress+ii));
  strcat(string,s);
}
errormessage(string);
*/
	estiff=assemgmtx(elem,estress);
    estiff=modifyhinge(elem,estiff);               /*MODIFY MATRIX.*/
    estiff=transformation(estiff,tmatrix);         /*[K]=[Tt][k][T]*/
    for(ii=0;ii<12;ii++)
    {
      for(jj=0;jj<12;jj++) *(*(estiff+ii)+jj)=-*(*(estiff+ii)+jj);
    }

    assemgstiffness(gmtx,estiff,&elem);     /*ASSEMBLAGE GEOMETRIC.*/

    for(ii=0;ii<=11;ii++) free(*(estiff+ii));
    free(estiff);

    for(ii=0;ii<=2;ii++) free(*(drccos+ii));
	free(drccos);
	for(ii=0;ii<=11;ii++) free(*(tmatrix+ii));
    free(tmatrix);
  }
  free(estress);
  laptime("GLOBAL MATRIX ASSEMBLED.",t0);

  /*currentvalue("GLOBAL MATRIX:[K]",msize,neig,kmtx,NULL,NULL,NULL);*/
  /*currentvalue("GLOBAL MATRIX:[G]",msize,neig,gmtx,NULL,NULL,NULL);*/

/*
for(ii=1;ii<=msize;ii++)
{
  sprintf(string,"%3d",ii);
  for(jj=1;jj<=msize;jj++)
  {
	gread(gmtx,ii,jj,&gdata);
	sprintf(s," %12.3f",gdata);
	strcat(string,s);
jj++;
  }
  if(fout!=NULL) fprintf(fout,"%s\n",string);
ii++;
}
for(ii=1;ii<=msize;ii++)
{
  sprintf(string,"%3d",ii);
  for(jj=1;jj<=msize;jj++)
  {
	gread(kmtx,ii,jj,&gdata);
	sprintf(s," %12.3f",gdata);
	strcat(string,s);
jj++;
  }
  if(fout!=NULL) fprintf(fout,"%s\n",string);
ii++;
}
*/

  deigabgeneral(gmtx,kmtx,confs,msize,neig,neig,eps,eigen,gvct);
  laptime("EIGEN COMPLETED.",t0);

  for(i=0;i<neig;i++)
  {
    *(eigen+i)=1.0/(*(eigen+i));
	sprintf(string,"EIGEN VALUE %ld=%.5E",(i+1),*(eigen+i));
	fprintf(fout,"%s\n",string);
    errormessage(string);

    outputmode(*(gvct+i),fout,nnode,ninit);
  }

  updatemode(af,*(gvct+0)); /*FORMATION UPDATE.*/

  af->nlaps=neig;
  af->eigenval=eigen;
  af->eigenvec=gvct;

  gfree(kmtx,nnode); /*FREE GLOBAL MATRIX.*/
  gfree(gmtx,nnode); /*FREE GLOBAL MATRIX.*/

  errormessage(" ");
  errormessage("COMPLETED.");
  fprintf(fout,"COMPLETED.\n");

  fclose(fout);

  memory1=availablephysicalmemory("REMAIN:");
  sprintf(string,"CONSUMPTION:%ld[BYTES]",(memory0-memory1));
  errormessage(string);

  return 1;
}/*bclng001*/

int bclng002(struct arclmframe *af)
/*ELASTIC BUCKLING FOR ARCLM FRAME.*/
/*AFTER "ARCLM001".*/
{
  DWORD memory0,memory1;

  FILE /**fin,*/*fout,*feig,*frat;                   /*FILE 8 BYTES*/
  /*FILE *felem,*fdisp,*freact;*/
  /*char dir[]=DIRECTORY;*/                        /*DATA DIRECTORY*/
  char s[80],string[1024];
  int i,j,ii,jj,mm,nn;
  int nnode,nelem,nsect/*,nreact*/;
  long int /*fsize,*/loff,ioff,joff,msize;
  /*long int time;*/
  struct gcomponent ginit={0,0,0.0,NULL};
  struct gcomponent *kmtx,*gmtx;                    /*GLOBAL MATRIX*/
  double **kcmtx,**gcmtx,**ktmtx,**gtmtx;        /*CONDENSED MATRIX*/
  double **gvct;                                    /*GLOBAL VECTOR*/
  double **drccos,**tmatrix,**estiff,*estress;
  /*double determinant,data;*/
  double gdata;
  clock_t t0/*,t1,t2*/;

  /*struct osect *sects;*/
  struct onode /**nodes,*/*ninit;
  struct owire elem;
  struct owire *elems;
  struct oconf *confs;

  long int neig;
  double AA[12][12],BB[12][12],WW[12][12],EE[12],VV[12][12];
  signed char CF[12];
  double eps=1.0E-16,*eigen;

  memory0=availablephysicalmemory("INITIAL:");   /*MEMORY AVAILABLE*/

  fout=fgetstofopen("\0","w",ID_OUTPUTFILE);          /*OUTPUT FILE*/
  globalfile=fout;

  GetDlgItemText((wmenu.childs+2)->hwnd,ID_OUTPUTFILE,s,80);
  /*GetDlgItemText((wmenu.childs+2)->hwnd,ID_SECTIONFILE,s,80);*/
  nn=strcspn(s,".");
  strncpy(string,s,nn);
  string[nn]='\0';
  strcat(string,".eig");
  feig=fopen(string,"w");

  sprintf(string,"\0");
  strncpy(string,s,nn);
  string[nn]='\0';
  strcat(string,".rat");
  frat=fopen(string,"w");

  fprintf(feig,"ELEMENT EIGEN VALUES\n");

  t0=clock();                                        /*CLOCK BEGIN.*/

  nnode=af->nnode;                                 /*INPUT INITIAL.*/
  nelem=af->nelem;
  nsect=af->nsect;
  /*nreact=af->nreact;*/
  sprintf(string,"NODES=%d ELEMS=%d SECTS=%d",nnode,nelem,nsect);
  errormessage(string);
  fprintf(fout,"%s\n",string);

  msize=6*nnode;                           /*SIZE OF GLOBAL MATRIX.*/
  neig=NEIGEN;

  kmtx=(struct gcomponent *)          /*DIAGONALS OF GLOBAL MATRIX.*/
		malloc(msize*sizeof(struct gcomponent));
  gmtx=(struct gcomponent *)
		malloc(msize*sizeof(struct gcomponent));
  /*kcmtx=(struct gcomponent *)
		 malloc(msize*sizeof(struct gcomponent));*/
  /*gcmtx=(struct gcomponent *)
		 malloc(msize*sizeof(struct gcomponent));*/
  kcmtx=(double **)malloc(msize*sizeof(double *));
  gcmtx=(double **)malloc(msize*sizeof(double *));
  ktmtx=(double **)malloc(msize*sizeof(double *));
  gtmtx=(double **)malloc(msize*sizeof(double *));
  if(kmtx==NULL || gmtx==NULL) return 0;
  if(kcmtx==NULL || gcmtx==NULL) return 0;
  if(ktmtx==NULL || gtmtx==NULL) return 0;
  for(i=0;i<msize;i++)
  {
	(kmtx+i)->down=NULL;            /*GLOBAL MATRIX INITIALIZATION.*/
	(gmtx+i)->down=NULL;
	/*(kcmtx+i)->down=NULL;*/
	/*(gcmtx+i)->down=NULL;*/
	*(kcmtx+i)=(double *)malloc(msize*sizeof(double));
	*(gcmtx+i)=(double *)malloc(msize*sizeof(double));
	*(ktmtx+i)=(double *)malloc(msize*sizeof(double));
	*(gtmtx+i)=(double *)malloc(msize*sizeof(double));
  }

  gvct=(double **)malloc(neig*sizeof(double *));   /*GLOBAL VECTORS*/
  eigen=(double *)malloc(neig*sizeof(double));       /*EIGEN VALUES*/
  if(gvct==NULL || eigen==NULL) return 0;
  for(i=0;i<neig;i++)
  {
	*(gvct+i)=(double *)malloc(msize*sizeof(double));
	for(j=0;j<msize;j++) *(*(gvct+i)+j)=0.0;
  }

  /*fdisp=af->fdisp;*/                 /*DISPLACEMENT:6 DIRECTIONS.*/
  /*felem=af->felem;*/              /*CODE,12 BOUNDARIES,12 STRESS.*/
  /*freact=af->freact;*/                           /*REACTION FILE.*/

  /*sects=af->sects;*/
  /*nodes=af->nodes;*/
  ninit=af->ninit;
  elems=af->elems;
  confs=af->confs;

  GetAsyncKeyState(VK_LBUTTON);                   /*CLEAR KEY LEFT.*/
  GetAsyncKeyState(VK_RBUTTON);                  /*CLEAR KEY RIGHT.*/

  errormessage("BCLNG001:BUCKLING LINEAR.");
  availablephysicalmemory("REMAIN:");            /*MEMORY AVAILABLE*/

  for(i=1;i<=msize;i++)             /*GLOBAL MATRIX INITIALIZATION.*/
  {
	ginit.m=(unsigned short int)i;
	/*ginit.n=(unsigned short int)i;*/
	*(kmtx+(i-1))=ginit;
	*(gmtx+(i-1))=ginit;
  }
  /*comps=msize;*/ /*INITIAL COMPONENTS=DIAGONALS.*/

  laptime("ASSEMBLING GLOBAL MATRIX.",t0);

  estress=(double *)malloc(12*sizeof(double));

  for(i=1;i<=nelem;i++)                 /*ASSEMBLAGE GLOBAL MATRIX.*/
  {
	/*elem=*(elems+i-1);*/                     /*READ ELEMENT DATA.*/
	inputelem(elems,af->melem,i-1,&elem);

	for(ii=0;ii<=1;ii++) /*INITIALIZE COORDINATION.*/
	{
	  loff=elem.node[ii]->loff;
	  for(jj=0;jj<3;jj++)
	  {
		elem.node[ii]->d[jj]=(ninit+loff)->d[jj];
	  }
	}

	drccos=directioncosine(elem.node[0]->d[0],
						   elem.node[0]->d[1],
						   elem.node[0]->d[2],
						   elem.node[1]->d[0],
						   elem.node[1]->d[1],
						   elem.node[1]->d[2],
						   elem.cangle);                 /*[DRCCOS]*/

	tmatrix=transmatrix(drccos);           /*TRANSFORMATION MATRIX.*/
	estiff=assememtx(elem);            /*ELASTIC MATRIX OF ELEMENT.*/
	estiff=modifyhinge(elem,estiff);               /*MODIFY MATRIX.*/
	estiff=transformation(estiff,tmatrix);         /*[K]=[Tt][k][T]*/

/*
fprintf(fout,"Element Matrix %d [ke%d]\n",i,i);
for(ii=0;ii<12;ii++)
{
  sprintf(string,"\0");
  for(jj=0;jj<=ii;jj++)
  {
	sprintf(s," %12.5E",*(*(estiff+ii)+jj));
	strcat(string,s);
  }
  if(fout!=NULL) fprintf(fout,"%s\n",string);
  errormessage(string);
}
*/

	assemgstiffness(kmtx,estiff,&elem);       /*ASSEMBLAGE ELASTIC.*/

	for(ii=0;ii<=11;ii++) free(*(estiff+ii));
	free(estiff);

	for(ii=0;ii<=1;ii++)                                /*STRESSES.*/
	{
	  for(jj=0;jj<6;jj++) *(estress+6*ii+jj)=elem.stress[ii][jj];
	}
/*
sprintf(string,"\0");
for(ii=0;ii<12;ii++)
{
  sprintf(s," %12.5E",*(estress+ii));
  strcat(string,s);
}
errormessage(string);
*/
	estiff=assemgmtx(elem,estress);
	estiff=modifyhinge(elem,estiff);               /*MODIFY MATRIX.*/
	estiff=transformation(estiff,tmatrix);         /*[K]=[Tt][k][T]*/
	for(ii=0;ii<12;ii++)
	{
	  for(jj=0;jj<12;jj++) *(*(estiff+ii)+jj)=-*(*(estiff+ii)+jj);
	}

	assemgstiffness(gmtx,estiff,&elem);     /*ASSEMBLAGE GEOMETRIC.*/

	for(ii=0;ii<=11;ii++) free(*(estiff+ii));
	free(estiff);

	for(ii=0;ii<=2;ii++) free(*(drccos+ii));
	free(drccos);
	for(ii=0;ii<=11;ii++) free(*(tmatrix+ii));
	free(tmatrix);
  }
  free(estress);
  laptime("GLOBAL MATRIX ASSEMBLED.",t0);

  /*currentvalue("GLOBAL MATRIX:[Ke]",msize,neig,kmtx,NULL,NULL,NULL);*/
  /*currentvalue("GLOBAL MATRIX:[Kg]",msize,neig,gmtx,NULL,NULL,NULL);*/

/*
fprintf(fout,"GLOBAL ELASTIC MATRIX [Ke]\n");
for(ii=1;ii<=msize;ii++)
{
  sprintf(string,"%3d",ii);
  for(jj=1;jj<=msize;jj++)
  {
	gread(kmtx,ii,jj,&gdata);
	sprintf(s," %12.3f",gdata);
	strcat(string,s);
  }
  if(fout!=NULL) fprintf(fout,"%s\n",string);
}
fprintf(fout,"GLOBAL GEOMETRIC MATRIX [Kg]\n");
for(ii=1;ii<=msize;ii++)
{
  sprintf(string,"%3d",ii);
  for(jj=1;jj<=msize;jj++)
  {
	gread(gmtx,ii,jj,&gdata);
	sprintf(s," %12.3f",gdata);
	strcat(string,s);
  }
  if(fout!=NULL) fprintf(fout,"%s\n",string);
}
*/

  /*CONDENSATION*/
fprintf(fout,"CONDENSATION BEGIN\n");

  for(i=0;i<nelem;i++)
  {

fprintf(fout,"\nELEM %d OFFSET=%d\n",(af->elems+i)->code,i+1);

	for(ii=0;ii<msize;ii++) /*DUPLICATE [Ke],[Kg]*/
	{
	  for(jj=0;jj<msize;jj++)
	  {
		gread(kmtx,ii+1,jj+1,&gdata);
		*(*(kcmtx+ii)+jj)=gdata;
		if(ii!=jj) *(*(kcmtx+jj)+ii)=gdata;
		gread(gmtx,ii+1,jj+1,&gdata);
		*(*(gcmtx+ii)+jj)=gdata;
		if(ii!=jj) *(*(gcmtx+jj)+ii)=gdata;
	  }
	}

/*
fprintf(fout,"[Ke]\n");
for(ii=0;ii<msize;ii++)
{
  sprintf(string,"\0");
  for(jj=0;jj<msize;jj++)
  {
	sprintf(s," %12.3f",*(*(kcmtx+ii)+jj));
	strcat(string,s);
  }
  fprintf(fout,"%s\n",string);
}
fprintf(fout,"[Kg]\n");
for(ii=0;ii<msize;ii++)
{
  sprintf(string,"\0");
  for(jj=0;jj<msize;jj++)
  {
	sprintf(s," %12.3f",*(*(gcmtx+ii)+jj));
	strcat(string,s);
  }
  fprintf(fout,"%s\n",string);
}
*/

	ioff=6*((af->elems+i)->node[0]->loff);
	joff=6*((af->elems+i)->node[1]->loff);
	if(ioff>joff) /*SWAP*/
	{
	  loff=ioff;
	  ioff=joff;
	  joff=loff;
	}

sprintf(string,"ELEM ORDER=%d",i+1);
MessageBox(NULL,string,"CONDENSE",MB_OK);

	/*CONDENSATION*/ /*WITH CONSIDERING CONF*/
	for(j=0;j<msize;j++)
	{
	  if(j==ioff) j+=6;
	  if(j==joff) j+=6;
	  if(j>=msize) break;

	  if(!(confs+j)->iconf)
	  {
		if(*(*(kcmtx+j)+j)==0.0)
		{
		  sprintf(string,"INSTABLE TERMINATION K%d%d=%9.3f",j+1,j+1,*(*(kcmtx+j)+j));
		  MessageBox(NULL,string,"CONDENSE",MB_OK);
		  return 0;
		}
/*
sprintf(string,"i=%d j=%d",i,j);
MessageBox(NULL,string,"CONDENSE",MB_OK);
*/
		if(j<ioff) ii=j;
		else ii=ioff;
		while(ii<msize)
		{
		  if(ioff+6<=j && j<joff && ii==ioff+6) ii=j;
		  if(j>joff && ii==ioff+6) ii=joff;
		  if(j>=joff+6 && ii==joff+6) ii=j;
/*
sprintf(string,"i=%d j=%d ii=%d",i,j,ii);
MessageBox(NULL,string,"CONDENSE",MB_OK);
*/
		  if(!(confs+ii)->iconf)
		  /*if(ii!=j)*/
		  {
			if(j<ioff) jj=j;
			else jj=ioff;
			while(jj<msize)
			{
			  if(ioff+6<=j && j<joff && jj==ioff+6) jj=j;
			  if(j>joff && jj==ioff+6) jj=joff;
			  if(j>=joff+6 && jj==joff+6) jj=j;
/*
sprintf(string,"i=%d j=%d ii=%d jj=%d",i,j,ii,jj);
MessageBox(NULL,string,"CONDENSE",MB_OK);
*/
			  if(!(confs+jj)->iconf)
			  /*if(jj!=j)*/
			  {
				*(*(ktmtx+ii)+jj)=(*(*(kcmtx+ii)+jj))-(*(*(kcmtx+ii)+j))*(*(*(kcmtx+j)+jj))/(*(*(kcmtx+j)+j));
				/*
				*(*(gtmtx+ii)+jj)=(*(*(gcmtx+ii)+jj))-(*(*(kcmtx+ii)+j))*(*(*(gcmtx+j)+jj))/(*(*(kcmtx+j)+j));
				*/
				if(*(*(gcmtx+j)+j)==0.0)
				{
				  /*
				  *(*(gtmtx+ii)+jj)=0.0;
				  */
				}
				else
				{
				  *(*(gtmtx+ii)+jj)=(*(*(gcmtx+ii)+jj))-(*(*(gcmtx+ii)+j))*(*(*(gcmtx+j)+jj))/(*(*(gcmtx+j)+j));
				}

			  }
			  jj++;
			}
		  }
		  ii++;
		}

		for(ii=0;ii<msize;ii++)
		{
		  if(!(confs+ii)->iconf)
		  {
			for(jj=0;jj<msize;jj++)
			{
			  if(!(confs+jj)->iconf)
			  {
				*(*(kcmtx+ii)+jj)=*(*(ktmtx+ii)+jj);
				*(*(gcmtx+ii)+jj)=*(*(gtmtx+ii)+jj);
			  }
			}
		  }
		}

/*
fprintf(fout,"CONDENSED LINE=%d\n",j+1);
fprintf(fout,"ELEM %d ORDER=%d\n",(af->elems+i)->code,i+1);
fprintf(fout,"[Ke']\n");
for(ii=0;ii<msize;ii++)
{
  sprintf(string,"%3d",ii+1);
  for(jj=0;jj<msize;jj++)
  {
	sprintf(s," %12.3f",*(*(kcmtx+ii)+jj));
	strcat(string,s);
  }
  fprintf(fout,"%s\n",string);
}
fprintf(fout,"[Kg']\n");
for(ii=0;ii<msize;ii++)
{
  sprintf(string,"%3d",ii+1);
  for(jj=0;jj<msize;jj++)
  {
	sprintf(s," %12.3f",*(*(gcmtx+ii)+jj));
	strcat(string,s);
  }
  fprintf(fout,"%s\n",string);
}
*/
	  }
	} /*END CONDENSATION*/


fprintf(fout,"[ke']\n");
for(ii=ioff;ii<joff+6;ii++)
{
  if(ii==ioff+6) ii=joff;
  sprintf(string,"%3d",ii+1);
  for(jj=ioff;jj<joff+6;jj++)
  {
	if(jj==ioff+6) jj=joff;
	sprintf(s," %18.8f",*(*(kcmtx+ii)+jj));
	strcat(string,s);
  }
  fprintf(fout,"%s\n",string);
}
fprintf(fout,"[kg']\n");
for(ii=ioff;ii<joff+6;ii++)
{
  if(ii==ioff+6) ii=joff;
  sprintf(string,"%3d",ii+1);
  for(jj=ioff;jj<joff+6;jj++)
  {
	if(jj==ioff+6) jj=joff;
	sprintf(s," %18.8f",*(*(gcmtx+ii)+jj));
	strcat(string,s);
  }
  fprintf(fout,"%s\n",string);
}

/*
fprintf(fout,"2D Part of [ke'][kg']\n");
fprintf(fout,"[ke']\n");
for(ii=ioff;ii<joff+6;ii++)
{
  if(ii==ioff+6) ii=joff;
  sprintf(string,"\0");
  for(jj=ioff;jj<joff+6;jj++)
  {
	if(jj==ioff+6) jj=joff;
	sprintf(s," %12.5f",*(*(kcmtx+ii)+jj));
	strcat(string,s);
	jj++;
  }
  fprintf(fout,"%s\n",string);
  ii++;
}
fprintf(fout,"[kg']\n");
for(ii=ioff;ii<joff+6;ii++)
{
  if(ii==ioff+6) ii=joff;
  sprintf(string,"\0");
  for(jj=ioff;jj<joff+6;jj++)
  {
	if(jj==ioff+6) jj=joff;
	sprintf(s," %12.5f",*(*(gcmtx+ii)+jj));
	strcat(string,s);
	jj++;
  }
  fprintf(fout,"%s\n",string);
  ii++;
}
*/

	/*EXTRACT CONCERNING ELEMENT LINES*/
	mm=0;
	for(ii=ioff;ii<joff+6;ii++)
	{
	  if(ii==ioff+6) ii=joff;

	  CF[mm]=(af->confs+ii)->iconf;

	  nn=0;
	  for(jj=ioff;jj<joff+6;jj++)
	  {
		if(jj==ioff+6) jj=joff;
		BB[mm][nn]=*(*(kcmtx+ii)+jj);
		nn++;
	  }
	  mm++;
	}
	mm=0;
	for(ii=ioff;ii<joff+6;ii++)
	{
	  if(ii==ioff+6) ii=joff;
	  nn=0;
	  for(jj=ioff;jj<joff+6;jj++)
	  {
		if(jj==ioff+6) jj=joff;
		AA[mm][nn]=*(*(gcmtx+ii)+jj);
		nn++;
	  }
	  mm++;
	}


fprintf(fout,"{CONF} :");
for(ii=0;ii<12;ii++) fprintf(fout," %3d",CF[ii]);
fprintf(fout,"\n");


	/*deigqr(AA,BB,12,12,12,12,eps,WW,EE,VV,CF);*/
	deigqrcf(AA,BB,12,12,12,12,eps,WW,EE,VV,CF);


fprintf(globalfile,"EIGEN VECTOR [V]\n");
for(ii=0;ii<12;ii++)
{
  if(!CF[ii])
  {
	for(jj=0;jj<12;jj++)
	{
	  if(!CF[jj]) fprintf(globalfile," %18.8f",VV[ii][jj]);
	}
	fprintf(globalfile,"\n");
  }
}
fprintf(globalfile,"\n");


fprintf(globalfile,"EIGEN VALUE INVERSE {E}\n");
for(ii=0;ii<12;ii++)
{
  if(!CF[ii]) fprintf(globalfile," %18.8f\n",EE[ii]);
}
fprintf(globalfile,"\n");


	for(ii=0;ii<12;ii++)
	{
	  if(!CF[ii])
	  {
		mm=ii;
		break;
	  }
	}
	if(EE[mm]==0.0) MessageBox(NULL,"EIGEN VALUE=0.0","CONDENSE",MB_OK);
	fprintf(fout,"ELEM %d EIGEN VALUE=%12.5f",(af->elems+i)->code,1/fabs(EE[mm]));
	/*fprintf(fout," LINE=%d",mm);*/
	fprintf(fout,"\n");

	fprintf(feig,"ELEM %d EIGEN VALUE=%12.5f SAFETY=%12.5f\n",(af->elems+i)->code,1/fabs(EE[mm]),fabs(EE[mm]));
	fprintf(frat,"ELEM: %5d SSECT: %4d %12.5f 0.00000 0.00000 0.00000\n",(af->elems+i)->code,(af->elems+i)->sect->code,fabs(EE[mm]));
  }

  laptime("EIGEN COMPLETED.",t0);

  af->nlaps=neig;
  af->eigenval=eigen;
  af->eigenvec=gvct;

  gfree(kmtx,nnode); /*FREE GLOBAL MATRIX.*/
  gfree(gmtx,nnode); /*FREE GLOBAL MATRIX.*/

  errormessage(" ");
  errormessage("COMPLETED.");
  fprintf(fout,"COMPLETED.\n");

  fclose(fout);

  memory1=availablephysicalmemory("REMAIN:");
  sprintf(string,"CONSUMPTION:%ld[BYTES]",(memory0-memory1));
  errormessage(string);

  return 1;
}/*bclng002*/

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

  if(A!=NULL)
  {
	for(i=0;i<n;i++)
	{
	  fprintf(stdout,"A%d:",i+1);
	  for(j=0;j<n;j++)
	  {
		fprintf(stdout," %12.5f",A[i][j]);
	  }
	  fprintf(stdout,"\n");
	}
  }

  if(E!=NULL)
  {
	fprintf(stdout,"\n");
	for(i=1;i<=ne;i++) fprintf(stdout,"           E%ld",i);
    fprintf(stdout,"\n");
    for(j=0;j<ne;j++)
    {
      fprintf(stdout," %12.5f",E[j]);
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
		fprintf(stdout," %12.5f",V[j][i]);
      }
	  fprintf(stdout,"\n");
    }
  }

  if(W!=NULL)
  {
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
  }

  fprintf(stdout,"%s\n",str);
  gets(non);
  if(strcmp(non,"\0")) exit(1);
  return;
}/*currentvalues*/

void deigqr(double A[][KSIZE],double B[][KSIZE],
			long int N,long int NSIZE,long int NE,long int NV,
			double EPS,double W[][KSIZE],
			double E[],double V[][KSIZE],
			signed char CF[]) /*CONF UNDER CONSTRUCTION*/
{
  /* SUBROUTINE FOR GENERALIZED EIGENVALUE PROBLEM                */
  /*  [A]{X} = LAMBDA [B]{X}                                      */
  /* FOR REAL ASYMMETRIC MATRICES [A] AND SYMMETRIC MATRICES [B], */
  /* THE LATTER BEING POSITIVE DEFINITE.                          */
  /*                                                              */
  /* USAGE:                                                       */
  /* CALL DEIGQR( A, B, N, NSIZE, NE, NV, EPS, W, E, V )          */
  /*                                                              */
  /* INPUT:                                                       */
  /*   A[N][N]  : 2-DIM. ARRAY CONTAINING REAL ASYMMETRIC MATRIX. */
  /*   B[N][N]  : 2-DIM. ARRAY CONTAINING REAL SYMMETRIC          */
  /*              POSITIVE DEFINITE MATRIX.                       */
  /*              BECAUSE OF CHOLESKI DECOMPOSITION.              */
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
  /*              IN ORDER OF ABSOLUTE VALUES                     */
  /*   V[NV][N] : 2-DIM. ARRAY CONTAINING THE OUTPUT EIGENVECTORS */
  /*              THE VECTOR {V[1,K], V[2,K],..., V[N,K]}         */
  /*              BELONGS TO THE EIGENVALUE E[K].                 */
  /*                                                              */
  /* WORKING SPACE:                                               */
  /*   W[7][N]  : 2-DIM. ARRAY  (7,N)  USED AS THE WORKING AREA.  */
  /*                                                              */
  /* NOTE:                                                        */
  /*  THE MATRICES [A] AND [B] WILL BE DESTROYED.                 */
  /* METHOD:                                                      */
  /*  CHOLESKY REDUCTION FOLLOWED BY HOUSEHOLDER DIAGONALIZATION. */

  char non[256];
  long int I,J,K,K1,nev,neva,nvec,R,RSUB1;
  double SUM,T1,T,Ep,En,value;
  long int ii,jj,kk,nn;
  double AA[KSIZE][KSIZE],LI[KSIZE][KSIZE];
  double QQ[KSIZE][KSIZE],RR[KSIZE][KSIZE],QR[KSIZE][KSIZE];
  double qi[KSIZE][KSIZE],qj[KSIZE],VV[KSIZE][KSIZE];


fprintf(globalfile,"[A]\n");
for(ii=0;ii<N;ii++)
{
  for(jj=0;jj<N;jj++) fprintf(globalfile," %12.5f",A[ii][jj]);
  fprintf(globalfile,"\n");
}
fprintf(globalfile,"\n");
fprintf(globalfile,"[B]\n");
for(ii=0;ii<N;ii++)
{
  for(jj=0;jj<N;jj++) fprintf(globalfile," %12.5f",B[ii][jj]);
  fprintf(globalfile,"\n");
}
fprintf(globalfile,"\n");


  /* CHECK INPUT DATA. */
  nev =NE;
  neva=labs(nev);
  nvec=NV;

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
	if(B[0][0]<=0.0)
	{
	  fprintf(stdout,"DEIGAB 1:MATRIX [B] IS NOT POSITIVE DEFINITE.");
	  gets(non);
	}
	else
	{
	  E[0]=A[0][0]/B[0][0];
	  B[0][0]=sqrt(1.0/B[0][0]);
	  V[0][0]=1.0;
	}
	return;
  }

  /* CHOLESKI DECOMPOSITION OF [B] INTO [L][Lt]. */
  /* NOTE DIAGONALS OF [L] ARE REMAINED INVERSE. */

  if(B[0][0]<=0.0)
  {
	fprintf(stdout,"DEIGAB 2:MATRIX [B] IS NOT POSITIVE DEFINITE.");
	gets(non);
	return;
  }
  T=sqrt(1.0/B[0][0]);
  B[0][0]=T;

  for(ii=1;ii<N;ii++) B[0][ii]=B[ii][0]*T;

  for(R=1;R<N;R++)
  {
	T=B[R][R];
	for(K=0;K<R;K++) T-=B[K][R]*B[K][R];
	if(T<=0.0)
	{
	  fprintf(stdout,"DEIGAB 3:MATRIX [B] IS NOT POSITIVE DEFINITE.\n");
	  fprintf(stdout,"T%ld=%9.5E-%9.5E=%9.5E\n",R,B[R][R],SUM,T);
	  gets(non);
	  return;
	}
	T=sqrt(1.0/T);
	B[R][R]=T;

	if(R>=N-1) break;

	for(I=R+1;I<N;I++)
	{
	  SUM=B[I][R];
	  for(K=0;K<=R-1;K++) SUM-=B[K][I]*B[K][R];
	  B[R][I]=SUM*T;
	}
  }

  /* CHOLESKI DECOMPOSITION IS COMPLETED.                */
  /* [B] = [L][LT],  THE UPPER TRIANGULAR MATRIX [LT] IS */
  /* STORED IN THE ARRAY [B].                            */
  /* NOW, PROCEED TO EVALUATE  [L]**(-1) [A] [LT]**(-1)  */
  /* FIRST, PREMULTIPLY  [L]**(-1).                      */

/*
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
*/

  for(ii=0;ii<N;ii++)
  {
	for(jj=0;jj<N;jj++)
	{
	  if(ii==jj) LI[ii][jj]=1.0;
	  else       LI[ii][jj]=0.0;
	}
  }
  for(J=0;J<nvec;J++)
  {
	K=N-1;
	while(1)
	{
	  T=LI[J][K]*B[K][K];
	  LI[J][K]=T;
	  K1=K-1;
	  if(K1<0) break;
	  for(R=0;R<=K1;R++)
	  {
		LI[J][R]-=B[R][K]*T;
	  }
	  K=K1;
	}
  }


fprintf(globalfile,"[L-1]\n");
for(ii=0;ii<N;ii++)
{
  for(jj=0;jj<N;jj++) fprintf(globalfile," %9.5f",LI[ii][jj]);
  fprintf(globalfile,"\n");
}
fprintf(globalfile,"\n");


  for(ii=0;ii<N;ii++)
  {
	for(jj=0;jj<N;jj++)
	{
	  AA[ii][jj]=0.0;
	  for(kk=0;kk<N;kk++) AA[ii][jj]+=LI[ii][kk]*A[kk][jj];
	}
  }

  /* NEXT, POSTMULTIPLY  [LT]**(-1). */

/*
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
*/

  for(ii=0;ii<N;ii++)
  {
	for(jj=0;jj<N;jj++)
	{
	  A[ii][jj]=0.0;
	  for(kk=0;kk<N;kk++) A[ii][jj]+=AA[ii][kk]*LI[jj][kk];
	}
  }

fprintf(globalfile,"TRANSFORMED [A']\n");
for(ii=0;ii<N;ii++)
{
  for(jj=0;jj<N;jj++) fprintf(globalfile," %12.8f",A[ii][jj]);
  fprintf(globalfile,"\n");
}
fprintf(globalfile,"\n");


fprintf(globalfile,"TRANSFORMED [B']\n");
for(ii=0;ii<N;ii++)
{
  for(jj=0;jj<N;jj++) fprintf(globalfile," %12.8f",B[ii][jj]);
  fprintf(globalfile,"\n");
}
fprintf(globalfile,"\n");


for(ii=0;ii<N;ii++)
{
  for(jj=ii;jj<N;jj++)
  {
	QR[ii][jj]=0.0;
	for(kk=0;kk<=ii;kk++)
	{
	  if(kk==ii && kk==jj) QR[ii][jj]+=1.0/(B[kk][ii]*B[kk][jj]);
	  else if(kk==ii)      QR[ii][jj]+=(1.0/B[kk][ii])*B[kk][jj];
	  else if(kk==jj)      QR[ii][jj]+=B[kk][ii]*(1.0/B[kk][jj]);
	  else                 QR[ii][jj]+=B[kk][ii]*B[kk][jj];
	}
  }
}
fprintf(globalfile,"[L][Lt]\n");
for(ii=0;ii<N;ii++)
{
  for(jj=0;jj<ii;jj++) fprintf(globalfile,"             ");
  for(jj=ii;jj<N;jj++) fprintf(globalfile," %12.8f",QR[ii][jj]);
  fprintf(globalfile,"\n");
}
fprintf(globalfile,"\n");

  /*
  for(ii=0;ii<N;ii++)
  {
	for(jj=ii+1;jj<N;jj++) A[ii][jj]=A[jj][ii];
  }
  */

fprintf(globalfile,"TRANSFORMED FULL [A']\n");
for(ii=0;ii<N;ii++)
{
  for(jj=0;jj<N;jj++) fprintf(globalfile," %12.8f",A[ii][jj]);
  fprintf(globalfile,"\n");
}
fprintf(globalfile,"\n");


for(ii=0;ii<N;ii++)
{
  for(jj=0;jj<N;jj++) AA[ii][jj]=A[ii][jj];
}


  for(ii=0;ii<N;ii++)
  {
	for(jj=0;jj<N;jj++) RR[ii][jj]=0.0;
  }

  nn=0;
  while(nn<1000)
  /*for(nn=0;nn<100;nn++)*/
  {
	/*[A]=[Q][R] DECOMPOSITION*/
	for(ii=0;ii<N;ii++)
	{
	  RR[ii][ii]=0.0;
	  for(jj=0;jj<ii;jj++)
	  {
		RR[jj][ii]=0.0;
		for(kk=0;kk<N;kk++) RR[jj][ii]+=A[kk][ii]*QQ[kk][jj];
	  }
	  for(jj=0;jj<N;jj++)
	  {
		QQ[jj][ii]=A[jj][ii];
		for(kk=0;kk<ii;kk++) QQ[jj][ii]-=RR[kk][ii]*QQ[jj][kk];
		RR[ii][ii]+=QQ[jj][ii]*QQ[jj][ii];
	  }
	  if(QQ[ii][ii]>0.0) RR[ii][ii]= sqrt(RR[ii][ii]);
	  else               RR[ii][ii]=-sqrt(RR[ii][ii]); /*CASE SIGN OF EIGEN VALUE = Qii*/

	  if(RR[ii][ii]==0.0)
	  {
		fprintf(globalfile,"QR INSTABLE AT LINE=%d Rii=%.5f\n",ii+1,RR[ii][ii]);
		return;
	  }
	  for(jj=0;jj<N;jj++) QQ[jj][ii]/=RR[ii][ii];
	}


fprintf(globalfile,"[Q]\n");
for(ii=0;ii<N;ii++)
{
  for(jj=0;jj<N;jj++) fprintf(globalfile," %12.8f",QQ[ii][jj]);
  fprintf(globalfile,"\n");
}
fprintf(globalfile,"\n");


fprintf(globalfile,"[R]\n");
for(ii=0;ii<N;ii++)
{
  for(jj=0;jj<ii;jj++) fprintf(globalfile,"             ");
  for(jj=ii;jj<N;jj++) fprintf(globalfile," %12.8f",RR[ii][jj]);
  fprintf(globalfile,"\n");
}
fprintf(globalfile,"\n");


for(ii=0;ii<N;ii++)
{
  for(jj=0;jj<N;jj++)
  {
	QR[ii][jj]=0.0;
	for(kk=0;kk<=jj;kk++) QR[ii][jj]+=QQ[ii][kk]*RR[kk][jj];
  }
}
fprintf(globalfile,"[Q][R]\n");
for(ii=0;ii<N;ii++)
{
  for(jj=0;jj<N;jj++) fprintf(globalfile," %12.8f",QR[ii][jj]);
  fprintf(globalfile,"\n");
}
fprintf(globalfile,"\n");


	/*EIGENVECTOR [Q]=[Q][Qi]*/
	if(nn==0)
	{
	  for(ii=0;ii<N;ii++)
	  {
		for(jj=0;jj<N;jj++) V[ii][jj]=QQ[ii][jj];
	  }
	}
	else
	{
	  /*for(ii=0;ii<N;ii++)
	  {
		for(jj=0;jj<N;jj++)
		{
		  qi[jj]=0.0;
		  for(kk=0;kk<N;kk++) qi[jj]+=V[ii][kk]*QQ[kk][jj];
		  V[ii][jj]=qi[jj];
		}
	  }*/
	  for(ii=0;ii<N;ii++)
	  {
		for(jj=0;jj<N;jj++)
		{
		  qi[ii][jj]=0.0;
		  for(kk=0;kk<N;kk++) qi[ii][jj]+=V[ii][kk]*QQ[kk][jj];
		}
	  }
	  for(ii=0;ii<N;ii++)
	  {
		for(jj=0;jj<N;jj++) V[ii][jj]=qi[ii][jj];
	  }
	}


fprintf(globalfile,"[V%d]\n",nn+1);
for(ii=0;ii<N;ii++)
{
  for(jj=0;jj<N;jj++) fprintf(globalfile," %12.8f",V[ii][jj]);
  fprintf(globalfile,"\n");
}
fprintf(globalfile,"\n");


	/*ITERATION CHECK*/
	/*En=RR[0][0];*/ /*TEST BY 1ST EIGENVALUE.*/
	En=RR[KSIZE-1][KSIZE-1]; /*TEST BY LAST EIGENVALUE.*/
	if(nn==0) Ep=En;
	else
	{
	  value=fabs(Ep-En);
	  if(value<EPS)
	  {
		fprintf(stdout,"ITERATION COMPLETE N=%d EPS=%.3E",nn+1,value);
		gets(non);
		break;
	  }
	  Ep=En;
	}

	/*[A]=[R][Q] RECOMPOSITION*/
	for(ii=0;ii<N;ii++)
	{
	  for(jj=0;jj<N;jj++)
	  {
		qj[jj]=0.0;
		for(kk=ii;kk<N;kk++) qj[jj]+=RR[ii][kk]*QQ[kk][jj];
		A[ii][jj]=qj[jj];
	  }
	}


	fprintf(globalfile,"[A']=[R][Q]\n");
	for(ii=0;ii<N;ii++)
	{
	  for(jj=0;jj<N;jj++) fprintf(globalfile," %12.8f",A[ii][jj]);
	  fprintf(globalfile,"\n");
	}
	fprintf(globalfile,"\n");

/*
for(ii=0;ii<N;ii++)
{
  for(jj=0;jj<N;jj++)
  {
	QR[ii][jj]=0.0;
	for(kk=0;kk<N;kk++) QR[ii][jj]+=RR[ii][kk]*QQ[kk][jj];
  }
}
fprintf(globalfile,"[R][Q]\n");
for(ii=0;ii<N;ii++)
{
  for(jj=0;jj<N;jj++) fprintf(globalfile," %12.8f",QR[ii][jj]);
  fprintf(globalfile,"\n");
}
fprintf(globalfile,"\n");
*/

	nn++;
	if(nn>=1000)
	{
	  fprintf(stdout,"ABORT : ITERATION OVER 1000 TIMES.");
	  gets(non);
	}
  }

  for(ii=0;ii<N;ii++) E[ii]=RR[ii][ii];

/*
fprintf(globalfile,"[R]\n");
for(ii=0;ii<N;ii++)
{
  for(jj=0;jj<N;jj++) fprintf(globalfile," %9.5f",RR[ii][jj]);
  fprintf(globalfile,"\n");
}
fprintf(globalfile,"\n");
*/

  for(ii=N-1;ii>=0;ii--)
  {
	for(jj=0;jj<ii;jj++) VV[ii][jj]=0.0;
	VV[ii][ii]=1.0;
	for(jj=ii+1;jj<N;jj++)
	{
	  VV[ii][jj]=0.0;
	  for(kk=ii+1;kk<=jj;kk++) VV[ii][jj]+=RR[ii][kk]*VV[kk][jj];
	  if(RR[jj][jj]==RR[ii][ii])
	  {
		fprintf(stdout,"R%d%d-R%d%d OVERFLOW.",jj+1,jj+1,ii+1,ii+1);
		gets(non);
		return;
	  }
	  VV[ii][jj]/=(RR[jj][jj]-RR[ii][ii]);
	}
  }


fprintf(globalfile,"[V'']\n");
for(ii=0;ii<N;ii++)
{
  for(jj=0;jj<N;jj++) fprintf(globalfile," %12.8f",VV[ii][jj]);
  fprintf(globalfile,"\n");
}
fprintf(globalfile,"\n");


for(ii=0;ii<N;ii++)
{
  for(jj=0;jj<N;jj++)
  {
	QR[ii][jj]=0.0;
	for(kk=0;kk<N;kk++) QR[ii][jj]+=RR[ii][kk]*VV[kk][jj];
  }
}
fprintf(globalfile,"[R][V'']\n");
for(ii=0;ii<N;ii++)
{
  for(jj=0;jj<N;jj++) fprintf(globalfile," %12.8f",QR[ii][jj]);
  fprintf(globalfile,"\n");
}
fprintf(globalfile,"\n");


fprintf(globalfile,"[V''][QRdiag]\n");
for(ii=0;ii<N;ii++)
{
  for(jj=0;jj<N;jj++) fprintf(globalfile," %12.8f",RR[jj][jj]*VV[ii][jj]);
  fprintf(globalfile,"\n");
}
fprintf(globalfile,"\n");


  for(ii=0;ii<N;ii++)
  {
	for(jj=0;jj<N;jj++)
	{
	  qi[ii][jj]=0.0;
	  for(kk=0;kk<N;kk++) qi[ii][jj]+=V[ii][kk]*VV[kk][jj];
	}
  }
  for(ii=0;ii<N;ii++)
  {
	for(jj=0;jj<N;jj++) V[ii][jj]=qi[ii][jj];
  }


fprintf(globalfile,"[A']\n");
for(ii=0;ii<N;ii++)
{
  for(jj=0;jj<N;jj++) fprintf(globalfile," %9.5f",AA[ii][jj]);
  fprintf(globalfile,"\n");
}
fprintf(globalfile,"\n");


fprintf(globalfile,"[V']\n");
for(ii=0;ii<N;ii++)
{
  for(jj=0;jj<N;jj++) fprintf(globalfile," %9.5f",V[ii][jj]);
  fprintf(globalfile,"\n");
}
fprintf(globalfile,"\n");


fprintf(globalfile,"{E}\n");
for(ii=0;ii<N;ii++) fprintf(globalfile," %9.5f\n",E[ii]);
fprintf(globalfile,"\n");


for(ii=0;ii<N;ii++)
{
  for(jj=0;jj<N;jj++)
  {
	QR[ii][jj]=0.0;
	for(kk=0;kk<N;kk++) QR[ii][jj]+=AA[ii][kk]*V[kk][jj];
  }
}
fprintf(globalfile,"[A'][V']\n");
for(ii=0;ii<N;ii++)
{
  for(jj=0;jj<N;jj++) fprintf(globalfile," %9.5f",QR[ii][jj]);
  fprintf(globalfile,"\n");
}
fprintf(globalfile,"\n");


fprintf(globalfile,"{E}[V']\n");
for(ii=0;ii<N;ii++)
{
  for(jj=0;jj<N;jj++) fprintf(globalfile," %9.5f",E[jj]*V[ii][jj]);
  fprintf(globalfile,"\n");
}
fprintf(globalfile,"\n");


  /*BACK TRANSFORMATION OF EIGENVECTORS [L-1][V].*/

/*
fprintf(globalfile,"[L-1]\n");
for(ii=0;ii<N;ii++)
{
  for(jj=0;jj<N;jj++) fprintf(globalfile," %9.5f",LI[ii][jj]);
  fprintf(globalfile,"\n");
}
fprintf(globalfile,"\n");
*/

  for(ii=0;ii<N;ii++)
  {
	for(jj=N-1;jj>=0;jj--)
	{
	  V[jj][ii]*=B[jj][jj];
	  for(kk=0;kk<jj;kk++) V[kk][ii]-=B[kk][jj]*V[jj][ii];
	}
  }

  /*
  for(ii=0;ii<N;ii++)
  {
	for(jj=0;jj<N;jj++)
	{
	  qi[ii][jj]=0.0;
	  for(kk=0;kk<N;kk++) qi[ii][jj]+=LI[kk][ii]*V[kk][jj];
	}
  }
  for(ii=0;ii<N;ii++)
  {
	for(jj=0;jj<N;jj++) V[ii][jj]=qi[ii][jj];
  }
  */

  return;
}/*deigqr*/

void deigqrcf(double A[][KSIZE],double B[][KSIZE],
			  long int N,long int NSIZE,long int NE,long int NV,
			  double EPS,double W[][KSIZE],
			  double E[],double V[][KSIZE],
			  signed char CF[]) /*CONF UNDER CONSTRUCTION*/
{
  /* SUBROUTINE FOR GENERALIZED EIGENVALUE PROBLEM                */
  /*  [A]{X} = LAMBDA [B]{X}                                      */
  /* FOR REAL ASYMMETRIC MATRICES [A] AND SYMMETRIC MATRICES [B], */
  /* THE LATTER BEING POSITIVE DEFINITE.                          */
  /*                                                              */
  /* USAGE:                                                       */
  /* CALL DEIGQR( A, B, N, NSIZE, NE, NV, EPS, W, E, V )          */
  /*                                                              */
  /* INPUT:                                                       */
  /*   A[N][N]  : 2-DIM. ARRAY CONTAINING REAL ASYMMETRIC MATRIX. */
  /*   B[N][N]  : 2-DIM. ARRAY CONTAINING REAL SYMMETRIC          */
  /*              POSITIVE DEFINITE MATRIX.                       */
  /*              BECAUSE OF CHOLESKI DECOMPOSITION.              */
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
  /*              IN ORDER OF ABSOLUTE VALUES                     */
  /*   V[NV][N] : 2-DIM. ARRAY CONTAINING THE OUTPUT EIGENVECTORS */
  /*              THE VECTOR {V[1,K], V[2,K],..., V[N,K]}         */
  /*              BELONGS TO THE EIGENVALUE E[K].                 */
  /*                                                              */
  /* WORKING SPACE:                                               */
  /*   W[7][N]  : 2-DIM. ARRAY  (7,N)  USED AS THE WORKING AREA.  */
  /*                                                              */
  /* NOTE:                                                        */
  /*  THE MATRICES [A] AND [B] WILL BE DESTROYED.                 */
  /* METHOD:                                                      */
  /*  CHOLESKY REDUCTION FOLLOWED BY HOUSEHOLDER DIAGONALIZATION. */

  char non[256],iflag;
  long int I,J,K,K1,nev,neva,nvec,R,RSUB1;
  double SUM,T1,T,Ep,En,value;
  long int ii,jj,kk,nn,n1,n2;
  double AA[KSIZE][KSIZE],LI[KSIZE][KSIZE];
  double QQ[KSIZE][KSIZE],RR[KSIZE][KSIZE],QR[KSIZE][KSIZE];
  double qi[KSIZE][KSIZE],qj[KSIZE],VV[KSIZE][KSIZE];
  double eps=1.0E-08;


fprintf(globalfile,"[A]\n");
for(ii=0;ii<N;ii++)
{
  for(jj=0;jj<N;jj++) fprintf(globalfile," %18.8f",A[ii][jj]);
  fprintf(globalfile,"\n");
}
fprintf(globalfile,"\n");
fprintf(globalfile,"[B]\n");
for(ii=0;ii<N;ii++)
{
  for(jj=0;jj<N;jj++) fprintf(globalfile," %18.8f",B[ii][jj]);
  fprintf(globalfile,"\n");
}
fprintf(globalfile,"\n");


  /*TEST CONF=1 FOR 0 LINE & ROW*/
  /*
  for(ii=0;ii<N;ii++)
  {
	if(!CF[ii])
	{
	  iflag=0;
	  for(jj=0;jj<N;jj++)
	  {
		if(!CF[jj])
		{
		  if(A[ii][jj]!=0.0) iflag=1;
		  if(A[jj][ii]!=0.0) iflag=1;
		}
	  }
	  if(iflag==0) CF[ii]=1;
	}
  }
  */


fprintf(globalfile,"[A] FREE CONF PART\n");
for(ii=0;ii<N;ii++)
{
  if(!CF[ii])
  {
	for(jj=0;jj<N;jj++)
	{
	  if(!CF[jj]) fprintf(globalfile," %18.8f",A[ii][jj]);
	}
	fprintf(globalfile,"\n");
  }
}
fprintf(globalfile,"\n");
fprintf(globalfile,"[B] FREE CONF PART\n");
for(ii=0;ii<N;ii++)
{
  if(!CF[ii])
  {
	for(jj=0;jj<N;jj++)
	{
	  if(!CF[jj]) fprintf(globalfile," %18.8f",B[ii][jj]);
	}
	fprintf(globalfile,"\n");
  }
}
fprintf(globalfile,"\n");


  /* CHECK INPUT DATA. */
  nev =NE;
  neva=labs(nev);
  nvec=NV;

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
	if(!CF[0])
	{
	  fprintf(stdout,"DEIGAB 1 : MSIZE=1 CONF=1 NON SOLUTION");
	  gets(non);
	}
	else if(B[0][0]<=0.0)
	{
	  fprintf(stdout,"DEIGAB 1 : MATRIX [B] IS NOT POSITIVE DEFINITE.");
	  gets(non);
	}
	else
	{
	  E[0]=A[0][0]/B[0][0];
	  B[0][0]=sqrt(1.0/B[0][0]);
	  V[0][0]=1.0;
	}
	return;
  }

  /* CHOLESKI DECOMPOSITION OF [B] INTO [L][Lt]. */
  /* NOTE DIAGONALS OF [L] ARE REMAINED INVERSE. */

  n1=-1;
  for(ii=0;ii<N;ii++)
  {
	if(!CF[ii])
	{
	  if(n1==-1) n1=ii;
	  n2=ii;
	}
  }
/*fprintf(globalfile,"n1=%d n2=%d\n",n1,n2);*/

  if(B[n1][n1]<=0.0)
  {
	fprintf(stdout,"DEIGAB 2:MATRIX [B] IS NOT POSITIVE DEFINITE.");
	gets(non);
	return;
  }
  T=sqrt(1.0/B[n1][n1]);
  B[n1][n1]=T;

  for(ii=n1+1;ii<N;ii++)
  {
	if(!CF[ii]) B[n1][ii]=B[ii][n1]*T;
  }

  for(R=n1+1;R<N;R++)
  {
	if(!CF[R])
	{
	  T=B[R][R];
	  for(K=0;K<R;K++)
	  {
		if(!CF[K])
		{
		  T-=B[K][R]*B[K][R];
/*
fprintf(globalfile,"R=%d K=%d Brr=%9.5f Bkr=%9.5f T=%9.5e 1/T=%9.5f\n",R,K,B[R][R],B[K][R],T,1.0/T);
*/
		}
	  }
	  /*if(T<=0.0)*/
	  if(T<=eps)
	  {
		/*fprintf(stdout,"DEIGQR 3:MATRIX [B] IS NOT POSITIVE DEFINITE.\n");*/
		/*fprintf(stdout,"T%ld=%9.5E-%9.5E=%9.5E\n",R,B[R][R],SUM,T);*/
		fprintf(globalfile,"DEIGQR 3:MATRIX [B] IS NOT POSITIVE DEFINITE.\n");
		fprintf(globalfile,"T%ld=%9.5E-%9.5E=%9.5E\n",R,B[R][R],SUM,T);
		/*gets(non);*/
		/*return;*/
		/*T=0.0;*/
		/*T=1.0;*/
		T=1/eps; /*?*/
	  }
	  else T=sqrt(1.0/T);

	  B[R][R]=T;
/*
if(T==0.0) fprintf(globalfile,"R=%d Brr=%9.5f T=%9.5f\n\n",R,B[R][R],T);
else       fprintf(globalfile,"R=%d Brr=%9.5f SQRT(1/T)=%9.5f\n\n",R,B[R][R],T);
*/
	  /*if(R>=N-1) break;*/
	  if(R>=n2) break;

	  for(I=R+1;I<N;I++)
	  {
		if(!CF[I])
		{
		  SUM=B[I][R];
		  for(K=0;K<=R-1;K++)
		  {
			if(!CF[K])
			{
			  SUM-=B[K][I]*B[K][R];
/*
fprintf(globalfile,"R=%d K=%d I=%d Bir= %12.8f Bki x Bkr = %12.8f x %12.8f SUM=%9.5e T=%12.8f\n",
		R,K,I,B[I][R],B[K][I],B[K][R],SUM,T);
*/
			  if(-eps<=SUM && SUM<=eps) SUM=0.0; /*?*/
			}
		  }
		  B[R][I]=SUM*T;
/*
fprintf(globalfile,"R=%d I=%d Bri=%12.8f\n\n",R,I,B[R][I]);
*/
		}
	  }
	}
  }

fprintf(globalfile,"DECOMPOSED [B] FREE CONF\n");
for(ii=0;ii<N;ii++)
{
  if(!CF[ii])
  {
	for(jj=0;jj<N;jj++)
	{
	  if(!CF[jj]) fprintf(globalfile," %18.8f",B[ii][jj]);
	}
	fprintf(globalfile,"\n");
  }
}
fprintf(globalfile,"\n");


  /* CHOLESKI DECOMPOSITION IS COMPLETED.                */
  /* [B] = [L][LT],  THE UPPER TRIANGULAR MATRIX [LT] IS */
  /* STORED IN THE ARRAY [B].                            */
  /* NOW, PROCEED TO EVALUATE  [L]**(-1) [A] [LT]**(-1)  */
  /* FIRST, PREMULTIPLY  [L]**(-1).                      */

/*
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
*/

  for(ii=0;ii<N;ii++)
  {
	for(jj=0;jj<N;jj++)
	{
	  if(ii==jj) LI[ii][jj]=1.0;
	  else       LI[ii][jj]=0.0;
	}
  }
  for(J=0;J<nvec;J++)
  {
	if(!CF[J])
	{
	  K=N-1;
	  while(1)
	  {
		if(K<0) break;
		if(!CF[K])
		{

/*fprintf(globalfile," LI%d%d x B%d%d = %18.8f x %18.8f = %18.8f / %9.5e\n",
		J,K,K,K,LI[J][K],B[K][K],LI[J][K]*B[K][K],LI[J][K]*B[K][K]);*/

		  T=LI[J][K]*B[K][K];
		  LI[J][K]=T;
		  /*if(K<0) break;*/
		  for(R=0;R<K;R++)
		  {
			if(!CF[R]) LI[J][R]-=B[R][K]*T;
		  }
		}
		K--;
	  }

/*
fprintf(globalfile,"[L-1] J=%d\n",J);
for(ii=0;ii<N;ii++)
{
  if(!CF[ii])
  {
	for(jj=0;jj<N;jj++)
	{
	  if(!CF[jj]) fprintf(globalfile," %9.5f",LI[ii][jj]);
	}
	fprintf(globalfile,"\n");
  }
}
fprintf(globalfile,"\n");
*/

	}
  }

/*
fprintf(globalfile,"[L-1]\n");
for(ii=0;ii<N;ii++)
{
  if(!CF[ii])
  {
	for(jj=0;jj<N;jj++)
	{
	  if(!CF[jj]) fprintf(globalfile," %18.8f",LI[ii][jj]);
	}
	fprintf(globalfile,"\n");
  }
}
fprintf(globalfile,"\n");
*/

  for(ii=0;ii<N;ii++)
  {
	if(!CF[ii])
	{
	  for(jj=0;jj<N;jj++)
	  {
		if(!CF[jj])
		{
		  AA[ii][jj]=0.0;
		  for(kk=0;kk<N;kk++)
		  {
			if(!CF[kk]) AA[ii][jj]+=LI[ii][kk]*A[kk][jj];
		  }
		  /*if(ii==jj && -eps<=AA[ii][jj] && AA[ii][jj]<=eps)
		  {
			fprintf(globalfile,"AA%d%d=%9.5e\n",ii,jj,AA[ii][jj]);
			AA[ii][jj]=eps;
		  }*/
		}
	  }
	}
  }

  /* NEXT, POSTMULTIPLY  [LT]**(-1). */

/*
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
*/

  for(ii=0;ii<N;ii++)
  {
	if(!CF[ii])
	{
	  for(jj=0;jj<N;jj++)
	  {
		if(!CF[jj])
		{
		  A[ii][jj]=0.0;
		  for(kk=0;kk<N;kk++)
		  {
			if(!CF[kk]) A[ii][jj]+=AA[ii][kk]*LI[jj][kk];
		  }
		}
	  }
	}
  }
/*
fprintf(globalfile,"TRANSFORMED [A']\n");
for(ii=0;ii<N;ii++)
{
  if(!CF[ii])
  {
	for(jj=0;jj<N;jj++)
	{
	  if(!CF[jj]) fprintf(globalfile," %18.8f",A[ii][jj]);
	}
	fprintf(globalfile,"\n");
  }
}
fprintf(globalfile,"\n");
*/
/*
fprintf(globalfile,"TRANSFORMED [B']\n");
for(ii=0;ii<N;ii++)
{
  if(!CF[ii])
  {
	for(jj=0;jj<N;jj++)
	{
	  if(!CF[jj]) fprintf(globalfile," %18.8f",B[ii][jj]);
	}
	fprintf(globalfile,"\n");
  }
}
fprintf(globalfile,"\n");
*/

/*CHECK [L][Lt]=[B]*/

for(ii=0;ii<N;ii++)
{
  if(!CF[ii])
  {
	for(jj=ii;jj<N;jj++)
	{
	  if(!CF[jj])
	  {
		QR[ii][jj]=0.0;
		for(kk=0;kk<=ii;kk++)
		{
		  if(!CF[kk])
		  {
			if(kk==ii && kk==jj && B[kk][ii]!=0.0 && B[kk][jj]!=0.0)
			{
			  QR[ii][jj]+=1.0/(B[kk][ii]*B[kk][jj]);
			}
			else if(kk==ii && B[kk][ii]!=0.0) QR[ii][jj]+=(1.0/B[kk][ii])*B[kk][jj];
			else if(kk==jj && B[kk][jj]!=0.0) QR[ii][jj]+=B[kk][ii]*(1.0/B[kk][jj]);
			else                              QR[ii][jj]+=B[kk][ii]*B[kk][jj];
		  }
		}
	  }
	}
  }
}
fprintf(globalfile,"CHECK:[L][Lt]=[B]\n");
for(ii=0;ii<N;ii++)
{
  if(!CF[ii])
  {
	for(jj=0;jj<ii;jj++)
	{
	  if(!CF[jj]) fprintf(globalfile,"                   ");
	}
	for(jj=ii;jj<N;jj++)
	{
	  if(!CF[jj]) fprintf(globalfile," %18.8f",QR[ii][jj]);
	}
	fprintf(globalfile,"\n");
  }
}
fprintf(globalfile,"\n");


  /*
  for(ii=0;ii<N;ii++)
  {
	for(jj=ii+1;jj<N;jj++) A[ii][jj]=A[jj][ii];
  }
  */


fprintf(globalfile,"TRANSFORMED FULL [A']\n");
for(ii=0;ii<N;ii++)
{
  if(!CF[ii])
  {
	fprintf(globalfile,"%3d",ii+1);
	for(jj=0;jj<N;jj++)
	{
	  if(!CF[jj]) fprintf(globalfile," %18.8f",A[ii][jj]);
	}
	fprintf(globalfile,"\n");
  }
}
fprintf(globalfile,"\n");


  for(ii=0;ii<N;ii++)
  {
	for(jj=0;jj<N;jj++)
	{
	  AA[ii][jj]=A[ii][jj];
	  RR[ii][jj]=0.0;
	}

	if(!CF[ii] && A[ii][ii]==0.0)
	{
	  CF[ii]=1;
	  fprintf(globalfile,"NOTE : SKIPPED LINE=%d DIAGONAL=0.0\n",ii+1);
	}
  }

  nn=0;
  while(nn<1000)
  /*for(nn=0;nn<100;nn++)*/
  {
	/*[A]=[Q][R] DECOMPOSITION*/
	for(ii=0;ii<N;ii++)
	{
	  if(!CF[ii])
	  {
		RR[ii][ii]=0.0;
		for(jj=0;jj<ii;jj++)
		{
		  if(!CF[jj])
		  {
			RR[jj][ii]=0.0;
			for(kk=0;kk<N;kk++)
			{
			  if(!CF[kk]) RR[jj][ii]+=A[kk][ii]*QQ[kk][jj];
			}
		  }
		}
		for(jj=0;jj<N;jj++)
		{
		  if(!CF[jj])
		  {
			QQ[jj][ii]=A[jj][ii];
			for(kk=0;kk<ii;kk++)
			{
			  if(!CF[kk]) QQ[jj][ii]-=RR[kk][ii]*QQ[jj][kk];
			}
			RR[ii][ii]+=QQ[jj][ii]*QQ[jj][ii];
		  }
		}
		if(QQ[ii][ii]>0.0) RR[ii][ii]= sqrt(RR[ii][ii]);
		else               RR[ii][ii]=-sqrt(RR[ii][ii]); /*CASE SIGN OF EIGEN VALUE = Qii*/

		if(RR[ii][ii]==0.0)
		{
		  fprintf(globalfile,"QR INSTABLE AT LINE=%d Rii=%.5f\n",ii+1,RR[ii][ii]);
		  return;
		}
		for(jj=0;jj<N;jj++)
		{
		  if(!CF[jj]) QQ[jj][ii]/=RR[ii][ii];
		}
	  }
	}

/*
fprintf(globalfile,"[Q]\n");
for(ii=0;ii<N;ii++)
{
  if(!CF[ii])
  {
	for(jj=0;jj<N;jj++)
	{
	  if(!CF[jj]) fprintf(globalfile," %18.8f",QQ[ii][jj]);
	}
	fprintf(globalfile,"\n");
  }
}
fprintf(globalfile,"\n");
*/
/*
fprintf(globalfile,"[R]\n");
for(ii=0;ii<N;ii++)
{
  if(!CF[ii])
  {
	for(jj=0;jj<ii;jj++)
	{
	  if(!CF[jj]) fprintf(globalfile,"                   ");
	}
	for(jj=ii;jj<N;jj++)
	{
	  if(!CF[jj]) fprintf(globalfile," %18.8f",RR[ii][jj]);
	}
	fprintf(globalfile,"\n");
  }
}
fprintf(globalfile,"\n");
*/
/*
for(ii=0;ii<N;ii++)
{
  if(!CF[ii])
  {
	for(jj=0;jj<N;jj++)
	{
	  if(!CF[jj])
	  {
		QR[ii][jj]=0.0;
		for(kk=0;kk<=jj;kk++)
		{
		  if(!CF[kk]) QR[ii][jj]+=QQ[ii][kk]*RR[kk][jj];
		}
	  }
	}
  }
}
fprintf(globalfile,"[Q][R]\n");
for(ii=0;ii<N;ii++)
{
  if(!CF[ii])
  {
	for(jj=0;jj<N;jj++)
	{
	  if(!CF[jj]) fprintf(globalfile," %18.8f",QR[ii][jj]);
	}
	fprintf(globalfile,"\n");
  }
}
fprintf(globalfile,"\n");
*/

	/*EIGENVECTOR [Q]=[Q][Qi]*/
	if(nn==0)
	{
	  for(ii=0;ii<N;ii++)
	  {
		if(!CF[ii])
		{
		  for(jj=0;jj<N;jj++)
		  {
			if(!CF[jj]) V[ii][jj]=QQ[ii][jj];
		  }
		}
	  }
	}
	else
	{
	  /*for(ii=0;ii<N;ii++)
	  {
		for(jj=0;jj<N;jj++)
		{
		  qi[jj]=0.0;
		  for(kk=0;kk<N;kk++) qi[jj]+=V[ii][kk]*QQ[kk][jj];
		  V[ii][jj]=qi[jj];
		}
	  }*/
	  for(ii=0;ii<N;ii++)
	  {
		if(!CF[ii])
		{
		  for(jj=0;jj<N;jj++)
		  {
			if(!CF[jj])
			{
			  qi[ii][jj]=0.0;
			  for(kk=0;kk<N;kk++)
			  {
				if(!CF[kk]) qi[ii][jj]+=V[ii][kk]*QQ[kk][jj];
			  }
			}
		  }
		}
	  }
	  for(ii=0;ii<N;ii++)
	  {
		if(!CF[ii])
		{
		  for(jj=0;jj<N;jj++)
		  {
			if(!CF[jj])	V[ii][jj]=qi[ii][jj];
		  }
		}
	  }
	}

/*
fprintf(globalfile,"[V%d]\n",nn+1);
for(ii=0;ii<N;ii++)
{
  if(!CF[ii])
  {
	for(jj=0;jj<N;jj++)
	{
	  if(!CF[jj]) fprintf(globalfile," %18.8f",V[ii][jj]);
	}
	fprintf(globalfile,"\n");
  }
}
fprintf(globalfile,"\n");
*/

	/*ITERATION CHECK*/
	/*TEST BY 1ST EIGENVALUE.*/
	for(ii=0;ii<KSIZE;ii++)
	{
	  if(!CF[ii]) break;
	}
	/*TEST BY LAST EIGENVALUE.*/
	/*for(ii=KSIZE-1;ii>=0;ii--)
	{
	  if(!CF[ii]) break;
	}*/

	En=RR[ii][ii];
	if(nn==0) Ep=En;
	else
	{
	  value=fabs(Ep-En);
	  if(value<EPS)
	  {
		fprintf(stdout,"ITERATION COMPLETE N=%d EPS=%.3E",nn+1,value);
		gets(non);
		break;
	  }
	  Ep=En;
	}

	/*[A]=[R][Q] RECOMPOSITION*/
	for(ii=0;ii<N;ii++)
	{
	  if(!CF[ii])
	  {
		for(jj=0;jj<N;jj++)
		{
		  if(!CF[jj])
		  {
			qj[jj]=0.0;
			for(kk=ii;kk<N;kk++)
			{
			  if(!CF[kk]) qj[jj]+=RR[ii][kk]*QQ[kk][jj];
			}
			A[ii][jj]=qj[jj];
		  }
		}
	  }
	}

	/*
	fprintf(globalfile,"[A']=[R][Q]\n");
	for(ii=0;ii<N;ii++)
	{
	  if(!CF[ii])
	  {
		for(jj=0;jj<N;jj++)
		{
		  if(!CF[jj]) fprintf(globalfile," %18.8f",A[ii][jj]);
		}
		fprintf(globalfile,"\n");
	  }
	}
	fprintf(globalfile,"\n");
	*/
/*
for(ii=0;ii<N;ii++)
{
  for(jj=0;jj<N;jj++)
  {
	QR[ii][jj]=0.0;
	for(kk=0;kk<N;kk++) QR[ii][jj]+=RR[ii][kk]*QQ[kk][jj];
  }
}
fprintf(globalfile,"[R][Q]\n");
for(ii=0;ii<N;ii++)
{
  for(jj=0;jj<N;jj++) fprintf(globalfile," %12.8f",QR[ii][jj]);
  fprintf(globalfile,"\n");
}
fprintf(globalfile,"\n");
*/

	nn++;
	if(nn>=1000)
	{
	  fprintf(stdout,"ABORT : ITERATION OVER 1000 TIMES.");
	  gets(non);
	}
  }

  for(ii=0;ii<N;ii++) E[ii]=RR[ii][ii];

/*
fprintf(globalfile,"[R]\n");
for(ii=0;ii<N;ii++)
{
  for(jj=0;jj<N;jj++) fprintf(globalfile," %9.5f",RR[ii][jj]);
  fprintf(globalfile,"\n");
}
fprintf(globalfile,"\n");
*/

  for(ii=N-1;ii>=0;ii--)
  {
	if(!CF[ii])
	{
	  for(jj=0;jj<ii;jj++)
	  {
		if(!CF[jj]) VV[ii][jj]=0.0;
	  }
	  VV[ii][ii]=1.0;
	  for(jj=ii+1;jj<N;jj++)
	  {
		if(!CF[jj])
		{
		  VV[ii][jj]=0.0;
		  for(kk=ii+1;kk<=jj;kk++)
		  {
			if(!CF[kk]) VV[ii][jj]+=RR[ii][kk]*VV[kk][jj];
		  }
		  if(RR[jj][jj]==RR[ii][ii])
		  {
			fprintf(stdout,"R%d%d-R%d%d OVERFLOW.",jj+1,jj+1,ii+1,ii+1);
			gets(non);
			return;
		  }
		  VV[ii][jj]/=(RR[jj][jj]-RR[ii][ii]);
		}
	  }
	}
  }

/*
fprintf(globalfile,"[V'']\n");
for(ii=0;ii<N;ii++)
{
  if(!CF[ii])
  {
	for(jj=0;jj<N;jj++)
	{
	  if(!CF[jj]) fprintf(globalfile," %18.8f",VV[ii][jj]);
	}
	fprintf(globalfile,"\n");
  }
}
fprintf(globalfile,"\n");
*/
/*
for(ii=0;ii<N;ii++)
{
  if(!CF[ii])
  {
	for(jj=0;jj<N;jj++)
	{
	  if(!CF[jj])
	  {
		QR[ii][jj]=0.0;
		for(kk=0;kk<N;kk++)
		{
		  if(!CF[kk]) QR[ii][jj]+=RR[ii][kk]*VV[kk][jj];
		}
	  }
	}
  }
}
fprintf(globalfile,"[R][V'']\n");
for(ii=0;ii<N;ii++)
{
  if(!CF[ii])
  {
	for(jj=0;jj<N;jj++)
	{
	  if(!CF[jj]) fprintf(globalfile," %18.8f",QR[ii][jj]);
	}
	fprintf(globalfile,"\n");
  }
}
fprintf(globalfile,"\n");
*/
/*
fprintf(globalfile,"[V''][QRdiag]\n");
for(ii=0;ii<N;ii++)
{
  if(!CF[ii])
  {
	for(jj=0;jj<N;jj++)
	{
	  if(!CF[jj]) fprintf(globalfile," %18.8f",RR[jj][jj]*VV[ii][jj]);
	}
	fprintf(globalfile,"\n");
  }
}
fprintf(globalfile,"\n");
*/

  for(ii=0;ii<N;ii++)
  {
	if(!CF[ii])
	{
	  for(jj=0;jj<N;jj++)
	  {
		if(!CF[jj])
		{
		  qi[ii][jj]=0.0;
		  for(kk=0;kk<N;kk++)
		  {
			if(!CF[kk]) qi[ii][jj]+=V[ii][kk]*VV[kk][jj];
		  }
		}
	  }
	}
  }
  for(ii=0;ii<N;ii++)
  {
	if(!CF[ii])
	{
	  for(jj=0;jj<N;jj++)
	  {
		if(!CF[jj]) V[ii][jj]=qi[ii][jj];
	  }
	}
  }


fprintf(globalfile,"TRANSFORMED [A']\n");
for(ii=0;ii<N;ii++)
{
  if(!CF[ii])
  {
	for(jj=0;jj<N;jj++)
	{
	  if(!CF[jj]) fprintf(globalfile," %18.8f",AA[ii][jj]);
	}
	fprintf(globalfile,"\n");
  }
}
fprintf(globalfile,"\n");


fprintf(globalfile,"EIGEN VECTOR [V']\n");
for(ii=0;ii<N;ii++)
{
  if(!CF[ii])
  {
	for(jj=0;jj<N;jj++)
	{
	  if(!CF[jj]) fprintf(globalfile," %18.8f",V[ii][jj]);
	}
	fprintf(globalfile,"\n");
  }
}
fprintf(globalfile,"\n");


fprintf(globalfile,"EIGEN VALUE {E}\n");
for(ii=0;ii<N;ii++)
{
  if(!CF[ii]) fprintf(globalfile," %18.8f\n",E[ii]);
}
fprintf(globalfile,"\n");


for(ii=0;ii<N;ii++)
{
  if(!CF[ii])
  {
	for(jj=0;jj<N;jj++)
	{
	  if(!CF[jj])
	  {
		QR[ii][jj]=0.0;
		for(kk=0;kk<N;kk++)
		{
		  if(!CF[kk]) QR[ii][jj]+=AA[ii][kk]*V[kk][jj];
		}
	  }
	}
  }
}
fprintf(globalfile,"CHECK [A'][V']\n");
for(ii=0;ii<N;ii++)
{
  if(!CF[ii])
  {
	for(jj=0;jj<N;jj++)
	{
	  if(!CF[jj]) fprintf(globalfile," %18.8f",QR[ii][jj]);
	}
	fprintf(globalfile,"\n");
  }
}
fprintf(globalfile,"\n");


fprintf(globalfile,"CHECK {E}[V']\n");
for(ii=0;ii<N;ii++)
{
  if(!CF[ii])
  {
	for(jj=0;jj<N;jj++)
	{
	  if(!CF[jj]) fprintf(globalfile," %18.8f",E[jj]*V[ii][jj]);
	}
	fprintf(globalfile,"\n");
  }
}
fprintf(globalfile,"\n");


  /*BACK TRANSFORMATION OF EIGENVECTORS [L-1][V].*/

/*
fprintf(globalfile,"[L-1]\n");
for(ii=0;ii<N;ii++)
{
  for(jj=0;jj<N;jj++) fprintf(globalfile," %9.5f",LI[ii][jj]);
  fprintf(globalfile,"\n");
}
fprintf(globalfile,"\n");
*/

  for(ii=0;ii<N;ii++)
  {
	if(!CF[ii])
	{
	  for(jj=N-1;jj>=0;jj--)
	  {
		if(!CF[jj])
		{
		  V[jj][ii]*=B[jj][jj];
		  for(kk=0;kk<jj;kk++)
		  {
			if(!CF[kk]) V[kk][ii]-=B[kk][jj]*V[jj][ii];
		  }
		}
	  }
	}
  }

  /*
  for(ii=0;ii<N;ii++)
  {
	for(jj=0;jj<N;jj++)
	{
	  qi[ii][jj]=0.0;
	  for(kk=0;kk<N;kk++) qi[ii][jj]+=LI[kk][ii]*V[kk][jj];
	}
  }
  for(ii=0;ii<N;ii++)
  {
	for(jj=0;jj<N;jj++) V[ii][jj]=qi[ii][jj];
  }
  */

  return;
}/*deigqrcf*/

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

  long int ii,jj;

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

  /*currentvalues("DEIGAB:[B] DECOMPOSED.",N,NE,B,NULL,NULL,NULL);*/

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

/*
if(globalfile!=NULL) fprintf(globalfile,"Transformed [A]\n");
for(ii=0;ii<24;ii++)
{
  for(jj=0;jj<24;jj++)
  {
    if(globalfile!=NULL) fprintf(globalfile," %12.5E",A[ii][jj]);
  }
  if(globalfile!=NULL) fprintf(globalfile,"\n");
}
*/

  /*currentvalues("DEIGAB:TRANSFORMED.",N,NE,A,NULL,NULL,NULL);*/

  /*FIND EIGENVALUES AND EIGENVECTORS OF THE TRANSFORMED MATRIX.*/
  deigrs( A, N, NSIZE, nev, nvec, EPS, W, W[6], E, V );

/*
if(globalfile!=NULL) fprintf(globalfile,"Transformed [A] after DEIGRS\n");
for(ii=0;ii<24;ii++)
{
  for(jj=0;jj<24;jj++)
  {
    if(globalfile!=NULL) fprintf(globalfile," %12.5E",A[ii][jj]);
  }
  if(globalfile!=NULL) fprintf(globalfile,"\n");
}
*/

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
  /*currentvalues("DEIGAB END.",N,NE,A,NULL,E,V);*/

  return;
}/*deigab*/

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

/*if(globalfile!=NULL && A[J][K]*A[I][J]!=0.0)
fprintf(globalfile,"1 A[%d][%d]*A[%d][%d] = %12.5E * %12.5E = %12.5E\n",
        J,K,I,J,A[J][K],A[I][J],A[J][K]*A[I][J]);*/
          }
          for(J=I;J<N;J++)
		  {
            S+=A[J][K]*A[J][I];

/*if(globalfile!=NULL && A[J][K]*A[J][I]!=0.0)
fprintf(globalfile,"2 A[%d][%d]*A[%d][%d] = %12.5E * %12.5E = %12.5E\n",
        J,K,J,I,A[J][K],A[J][I],A[J][K]*A[J][I]);*/
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
/*if(globalfile!=NULL)
fprintf(globalfile," %d %d %d : %12.5E\n",
        K,J,I,W[0][J]*A[I][K]+W[0][I]*A[J][K]);*/
          }
        }
      }
    }
    W[0][nm1-1]=A[N-1][nm1-1];

    /*currentvalues("DEIGRS:REDUCED.",N,NE,A,W,NULL,NULL);*/

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

void currentvalue(char *string,
				  long int n,long int ne,
				  struct gcomponent *A,
				  double **W,
				  double *E,double **V)
/*CHECK CURRENT VALUES FOR DEBUG.*/
{
  char /*non[10],*/s[80],str[400];
  long int i,j;
  double data;

  ne=labs(ne);

  if(A!=NULL)
  {
	for(i=1;i<=n;i++)
    {
      sprintf(str,"A%d:",i);
      for(j=1;j<=i;j++)
      {
        gread(A,i,j,&data);
        sprintf(s," %11.5f",data);
        strcat(str,s);
      }
      errormessage(str);
    }
  }

  if(E!=NULL)
  {
    errormessage("\0");
    sprintf(str,"\0");
    for(i=1;i<=ne;i++)
    {
      sprintf(s,"           E%ld",i);
      strcat(str,s);
    }
    errormessage(str);

    sprintf(str,"\0");
    for(j=0;j<ne;j++)
    {
	  sprintf(s," %12.5E",*(E+j));
      strcat(str,s);
    }
    errormessage(str);
  }

  if(V!=NULL)
  {
    errormessage("\0");
    sprintf(str,"\0");
    for(i=1;i<=ne;i++)
    {
      sprintf(s,"           V%ld",i);
      strcat(str,s);
    }
	errormessage(str);

    for(i=0;i<n;i++)
    {
      sprintf(str,"\0");
      for(j=0;j<ne;j++)
      {
        sprintf(s," %12.5E",*(*(V+j)+i));
        strcat(str,s);
      }
      errormessage(str);
    }
  }

  if(W!=NULL)
  {
    errormessage("\0");
    for(i=0;i<6;i++)
    {
      sprintf(str,"W%d:",i);
      for(j=0;j<n;j++)
      {
        sprintf(s," %12.5E",*(*(W+i)+j));
        strcat(str,s);
      }
      errormessage(str);
    }
  }

  errormessage(string);
  return;
}/*currentvalue*/

void deigabgeneral(struct gcomponent *A,
                   struct gcomponent *B,
                   struct oconf *confs,
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

  char str[256];
  long int i,ii,jj,kk;
  long int J,K,nev,neva,nvec,R;
  double T;

  unsigned short int mm,mmm;
  double lkj,lim;
  struct gcomponent *gi,*gj,*gr,*gk,*gp,*gg;
  struct gcomponent *bj,*bm;

  long int ipiv; /*PIVOT LINE.*/
  signed char ic; /*CONFINEMENT ID.*/
  long int ndim; /*MATRIX DIMENSIONS.*/

  long int iii,jjj;
  double data;

  errormessage("DEIGAB:BEGIN.");

/*
if(globalfile!=NULL) fprintf(globalfile,"[A] in DEIGAB\n");
for(ii=0;ii<N;ii++)
{
  for(jj=0;jj<N;jj++)
  {
    gread(A,ii+1,jj+1,&data);
    if(globalfile!=NULL) fprintf(globalfile," %12.5E",data);
  }
  if(globalfile!=NULL) fprintf(globalfile,"\n");
}
if(globalfile!=NULL) fprintf(globalfile,"\n");
if(globalfile!=NULL) fprintf(globalfile,"[B] in DEIGAB\n");
for(ii=0;ii<N;ii++)
{
  for(jj=0;jj<N;jj++)
  {
    gread(B,ii+1,jj+1,&data);
    if(globalfile!=NULL) fprintf(globalfile," %12.5E",data);
  }
  if(globalfile!=NULL) fprintf(globalfile,"\n");
}
if(globalfile!=NULL) fprintf(globalfile,"\n");
*/

  /*CHECK INPUT DATA.*/
  nev  = NE;
  neva = labs(nev);
  nvec = NV;

  if(N<=0 || neva==0 || N-neva<0 || nvec<0 || neva-nvec<0)
  {
    errormessage("DEIGAB:INVALID ARGUMENT.");
    sprintf(str," N,NE,NV = %ld,%ld,%ld",N,nev,nvec);
	errormessage(str);
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
      errormessage("DEIGAB 1:MATRIX [B] IS NOT POSITIVE DEFINITE.");
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

  /*currentvalue("DEIGABGENERAL:[A]",N,NE,A,NULL,NULL,NULL);*/
  /*currentvalue("DEIGABGENERAL:[B]",N,NE,B,NULL,NULL,NULL);*/

  T = (B+ipiv)->value;
  if(T<=0.0)
  {
	errormessage("DEIGAB 2:MATRIX [B] IS NOT POSITIVE DEFINITE.");
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

/*
if(globalfile!=NULL) fprintf(globalfile,"[B %d]\n",R);
for(iii=1;iii<=18;iii++)
{
  for(jjj=1;jjj<=18;jjj++)
  {
	gread(B,iii,jjj,&data);
    if(globalfile!=NULL) fprintf(globalfile," %12.5E",data);
  }
  if(globalfile!=NULL) fprintf(globalfile,"\n");
}
if(globalfile!=NULL) fprintf(globalfile,"\n");
*/

    ic=(confs+R)->iconf;
	if(!ic)
    {
      if((B+R)->value <= 0.0)
      {
        errormessage("DEIGAB 3:MATRIX [B] IS NOT POSITIVE DEFINITE.");
        sprintf(str,"T%ld=%9.5E",R,(B+R)->value);
        errormessage(str);
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

    currentpivot(R+1,N);
  }

  errormessage("DEIGAB:[B] DECOMPOSED.");
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

    currentpivot(J+1,N);
  }

/*
if(globalfile!=NULL) fprintf(globalfile,"Transformed [A]\n");
for(iii=int(N/2)+1;iii<=N;iii++)
{
  for(jjj=int(N/2)+1;jjj<=N;jjj++)
  {
	gread(A,iii,jjj,&data);
    if(globalfile!=NULL) fprintf(globalfile," %12.5E",data);
  }
  if(globalfile!=NULL) fprintf(globalfile,"\n");
}
*/

  errormessage("DEIGAB:[A] TRANSFORMED.");
  /*currentvalue("DEIGABGENERAL:TRANSFORMED.",N,NE,A,NULL,NULL,NULL);*/

  /*FIND EIGENVALUES AND EIGENVECTORS OF THE TRANSFORMED MATRIX.*/
  deigrsstandard(A,confs,N,nev,nvec,EPS,E,V);

/*
if(globalfile!=NULL) fprintf(globalfile,"Transformed [A] after DEIGRSSTANDARD\n");
for(iii=int(N/2)+1;iii<=N;iii++)
{
  for(jjj=int(N/2)+1;jjj<=N;jjj++)
  {
    gread(A,iii,jjj,&data);
    if(globalfile!=NULL) fprintf(globalfile," %12.5E",data);
  }
  if(globalfile!=NULL) fprintf(globalfile,"\n");
}
*/

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

    currentpivot(J+1,nvec);
  }

/*
if(globalfile!=NULL) fprintf(globalfile,"E1=%12.5E\n",*(E+0));
if(globalfile!=NULL) fprintf(globalfile,"{V} in DEIGAB\n");
for(ii=0;ii<N;ii++)
{
  if(!((confs+ii)->iconf) && globalfile!=NULL) fprintf(globalfile," %5d %12.5E\n",ii+1,*(*(V+0)+ii));
}
*/

  errormessage("DEIGAB:END.");
  /*currentvalue("DEIGABGENERAL END.",N,NE,NULL,NULL,E,V);*/

  return;
}/*deigabgeneral*/

void deigrsstandard(struct gcomponent *A,
					struct oconf *confs,
					long int N,long int NE,long int NV,
					double EPS,
					double *E,double **V)
{
  /* SUBROUTINE FOR STANDARD EIGENVALUE PROBLEM                   */
  /* HOUSEHOLDER'S TRIDIAGONAL REDUCTION.                         */
  /* EIGENVALUES BY BISECTION.                                    */
  /* EIGENVECTORS BY INVERSE ITERATION.                           */
  /*                                                              */
  /*  [A]{V} = LAMDA {V}                                         */
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

  char str[256];
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

  double data;

  errormessage("DEIGRS:BEGIN.");

/*
if(globalfile!=NULL) fprintf(globalfile,"[A] in DEIGRS\n");
for(ii=0;ii<N;ii++)
{
  if((confs+ii)->iconf!=1)
  {
    for(jj=0;jj<N;jj++)
	{
      if((confs+jj)->iconf!=1)
      {
        gread(A,ii+1,jj+1,&data);
        if(globalfile!=NULL) fprintf(globalfile," %18.14E",data);
      }
    }
    if(globalfile!=NULL) fprintf(globalfile,"\n");
  }
}
if(globalfile!=NULL) fprintf(globalfile,"\n");
*/

  NEA=labs(NE);
  if(NEA==0)
  {
    sprintf(str,"DEIGRS: NE = %ld",NE);
    errormessage(str);
    errormessage("NE SHOULD BE NON-ZERO.");
    errormessage("RETURN WITH NO CALCULATION.");
	return;
  }
  NVA=labs(NV);
  if(NVA>NEA || NEA>N)
  {
    sprintf(str,"DEIGRS: NV,NE,N = %ld,%ld,%ld",NV,NE,N);
    errormessage(str);
    sprintf(str,"NV,NE,N");
    strcat(str," SHOULD SATISFY THE FOLLOWING INEQUALITIES.");
    errormessage(str);
    errormessage("|NV|<=|NE|<=N RETURN WITH NO CALCULATION.");
    return;
  }

  W=(double **)malloc(6*sizeof(double *));
  for(i=0;i<6;i++)
  {
    ww=(double *)malloc(N*sizeof(double));
    *(W+i)=ww;

    /*for(j=0;j<N;j++)
    {
	  *(*(W+i)+j)=0.0;
    }*/
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
    errormessage("DEIGRS:ORDERING TOO MANY EIGEN VALUES.");
    errormessage("RETURN WITH NO CALCULATION.");
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
          for(I=K1;I<N;I++)
          {
/*S=0.0;*/
            ic=(confs+I)->iconf;
/*if(globalfile!=NULL) fprintf(globalfile,"conf %d=%d\n",I,ic);*/
            if(!ic)
            {
              gj=(A+K);
              while(gj->down!=NULL)
              {
                gj=gj->down;
                J=gj->m-1;
                ic=(confs+J)->iconf;
                if(!ic)
                {
                  if(J<I)
                  {
*(gpp+J)=(A+J);
                    while(((*(gpp+J))->m)<(I+1) &&
                          ((*(gpp+J))->down)!=NULL)
                    {
                      *(gpp+J)=(*(gpp+J))->down;
                    }

                    if(((*(gpp+J))->m)==(I+1)) /*{Woi}=R{Si}*/
                    {
                      S=(gj->value)
                       *((*(gpp+J))->value); /*Si={VT}{Ai}*/

                       *(*(W+0)+I)+=R*S;

/*if(globalfile!=NULL)
fprintf(globalfile,"1 A[%d][%d]*A[%d][%d] = %12.5E * %12.5E = %12.5E\n",
        J-24,K-24,I-24,J-24,(gj->value),(*(gpp+J))->value,S);*/
                    }
                  }
                  if(J>=I)
				  {
*(gpp+I)=(A+I);
                    while(((*(gpp+I))->m)<(J+1) &&
                          ((*(gpp+I))->down)!=NULL)
                    {
                      *(gpp+I)=(*(gpp+I))->down;
                    }

                    if(((*(gpp+I))->m)==(J+1)) /*{Woi}=R{Si}*/
                    {
                      S=(gj->value)
                       *((*(gpp+I))->value); /*Si={VT}{Ai}*/

                       *(*(W+0)+I)+=R*S;

/*if(globalfile!=NULL)
fprintf(globalfile,"2 A[%d][%d]*A[%d][%d] = %12.5E * %12.5E = %12.5E\n",
        J-24,K-24,J-24,I-24,(gj->value),(*(gpp+I))->value,S);*/
                    }
                  }
                }
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
          for(J=K1;J<N;J++)
          {
            jj=J+1;
            if((gk->m)<jj && gk->down!=NULL) gk=gk->down;

            ic=(confs+J)->iconf;
            if(!ic)
            {
              gi=gk;
              gj=(A+J);
              for(I=J;I<N;I++)
              {
                ii=I+1;
                if((gi->m)<ii && gi->down!=NULL) gi=gi->down;
                if((gj->m)<ii && gj->down!=NULL)
                {
                  gp=gj;
                  gj=gj->down;
				}

                ic=(confs+I)->iconf;
                if(!ic)
                {
                  S=0.0;
                  if(gi->m==ii) S+=(*(*(W+0)+J))*(gi->value); /*WjAik*/
                  if(gk->m==jj) S+=(*(*(W+0)+I))*(gk->value); /*WiAjk*/
/*if(globalfile!=NULL)
fprintf(globalfile," %d %d %d : %12.5E\n",
        K-24,J-24,I-24,S);*/

                  if(S!=0.0)
                  {
                    if(gj->m == ii)
                    {
                      gj->value-=S;
                    }
                    else if(gj->m < ii) /*ADD.*/
                    {
					  gg=gdefine((unsigned short int)ii,
                                 (unsigned short int)jj,
                                 -S,NULL,NULL);
                      gj->down=gg;
                      gp=gj;
                      gj=gg;
                    }
                    else if(gj->m > ii) /*FILL IN.*/
                    {
                      gg=gdefine((unsigned short int)ii,
                                 (unsigned short int)jj,
                                 -S,gj,NULL);
                      gp->down=gg;
                      gp=gg;
                    }
                  }
                }
              }
            }
          }
        }
      }
	  currentpivot(K,nm2);
    }

    gi=(A+nm2);
    while((gi->m)<(nm1+1) && gi->down!=NULL)
    {
      gi=gi->down;
    }

    if((gi->m)==(nm1+1)) *(*(W+0)+nm2)=gi->value;
    else                 *(*(W+0)+nm2)=0.0;

    errormessage("DEIGRS:HOUSEHOLDER REDUCED.");
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
/*if(NE<0) R=0;*/

    F=R; /*F:UPPER BOUND FOR DECENDING,LOWER FOR ACENDING.*/
    for(I=0;I<NEA;I++) *(E+I)=-R;

    for(K=0;K<NEA;K++)
    {
      currentpivot(K,NEA);

      D=*(E+K); /*D:LOWER BOUND FOR DECENDING,UPPER FOR ACENDING.*/
/*
if(globalfile!=NULL) fprintf(globalfile,"D%d = %12.5f\n",K+1,D);
*/

      while(1)
	  {
        T = 0.5*(D+F); /*DIVIDE SECTION.*/
/*
if(globalfile!=NULL) fprintf(globalfile,"T%d = 0.5 ( %12.5f + %12.5f ) = %12.5f\n",K+1,D,F,T);
*/

        if(fabs(D-F)<=EPS2 || T==D || T==F) break;

        J=0;
        I=0;
        while(1)
        {
          i=*(lines+I);

/*sprintf(str,"K=%d W[5][%d]=%.5E T=%.5E",K,i,*(*(W+5)+i),T);
MessageBox(NULL,str,"DEIGRS",MB_OK);*/

/*ERROR*/ Q=(*(*(W+5)+i))-T; /*ERROR AVOIDED WITHOUT ARCLM001 BEFORE GNSHN101.*/

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
/*ERROR*/ if(I>=ndim) break; /*?*/
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

if(globalfile!=NULL) fprintf(globalfile,"E%d = %18.14E\n",K+1,*(E+K));

    }
  }
  errormessage("DEIGRS:BISECTION END.");

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
    currentpivot(I,NVA);

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
  errormessage("DEIGRS:INVERSE ITERATION END.");

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

if(globalfile!=NULL) fprintf(globalfile,"E1=%12.5E\n",*(E+0));
if(globalfile!=NULL) fprintf(globalfile,"{V} in DEIGRS\n");
for(ii=0;ii<ndim;ii++)
{
  J=*(lines+ii);
  if(globalfile!=NULL) fprintf(globalfile," %5d %12.5E\n",J+1,*(*(V+0)+J));
}

  free(lines);
  free(LW);
  free(gpp);
  for(i=0;i<6;i++) free(*(W+i));
  free(W);

  errormessage("DEIGRS:END.");
  return;
}/*deigrsstandard*/

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

    gc=(struct gcomponent *)malloc(ncomp*sizeof(struct gcomponent));
    (gcpy+k)->down=(gc+0);

    go=(gmtx+k);
    while(go->down!=NULL) /*COPY COMPS IN COLUMN.*/
    {
      go=go->down;

      gc->m=go->m;
      gc->n=go->n;
      gc->value=go->value;

      if(go->down!=NULL)
	  {
        gc->down=(gc+1);
        gc++;
      }
      else gc->down=NULL;
    }
  }

  return gcpy;
}/*copygcompmatrix*/

double eigensubstitution(struct gcomponent *amtx,
                         struct gcomponent *bmtx,
                         struct oconf *confs,
                         long int msize,
                         double eig,double *vct)
{
  char str[256];
  long int i,j;
  double adata,bdata;
  double ai,bi,gosa;

  errormessage("[A]{v}=E[B]{v}");
  sprintf(str,"EIGEN VALUE=%.5E=%.5f",eig,eig);
  errormessage(str);

  gosa=0.0;
  for(i=1;i<=msize;i++)
  {
    ai=0.0;
    bi=0.0;
    if(!(confs+i-1)->iconf)
    {
      for(j=1;j<=msize;j++)
      {
        if(!(confs+j-1)->iconf)
        {
          gread(amtx,i,j,&adata);
          gread(bmtx,i,j,&bdata);

          ai+=adata*(*(vct+j-1));
          bi+=eig*bdata*(*(vct+j-1));
        }
      }
      sprintf(str," {A%ld}{V}=%9.5f  E{B%ld}{V}=%9.5f",i,ai,i,bi);
      errormessage(str);

      gosa+=(ai-bi)*(ai-bi);
    }
  }
  gosa=sqrt(gosa);
  sprintf(str," GOSA=%.5E\n",gosa);
  errormessage(str);

  return gosa;
}/*eigensubstitution*/

void outputmode(double *gvct,FILE *fout,int nnode,
                struct onode *nodes)
/*OUTPUT MODE DISPLACEMENT.*/
{
  char string[256];
  int i,j;
  double data[6];

  for(i=1;i<=nnode;i++)
  {
    for(j=0;j<=5;j++) data[j]=*(gvct+6*(i-1)+j);
    fprintf(fout,"NODE:%5ld {dU}=",(nodes+i-1)->code);
	sprintf(string," %12.5E %12.5E %12.5E %12.5E %12.5E %12.5E",
            data[0],data[1],data[2],data[3],data[4],data[5]);
    if(fout!=NULL) fprintf(fout,"%s\n",string);
  }
}/*outputmode*/

void updatemode(struct arclmframe *af,double *gvct)
/*FORMATION UPDATE.*/
{
  int i,j;
  long int nnode,loff;
  double data,ddata;

  nnode=af->nnode;

  for(i=0;i<nnode;i++)
  {
    for(j=0;j<=2;j++)
    {
      loff=6*i+j;
      ddata=*(gvct+loff);

      data=((af->ninit+i)->d[j])+ddata; /*{U}+{dU}*/
      (af->nodes+i)->d[j]=data;

      *(af->ddisp+loff)=data;
    }
    for(j=3;j<=5;j++)
    {
      loff=6*i+j;
      ddata=*(gvct+loff);
      *(af->ddisp+loff)=ddata;
    }
  }
  return;
}/*updatemode*/

void bclngoutputtomemory(FILE *ftext,struct arclmframe *af)
/*TRANSLATE BCLNG OUTPUTFILE TEXT INTO MEMORY.*/
{
  char **data,str[256]="\0";
  int i,j,n;
  long int ncode;
  double ddata;

  fseek(ftext,0L,SEEK_SET);

  af->ddisp=(double *)malloc(6*(af->nnode)*sizeof(double));
                                       /*DISPLACEMENT:6 DIRECTIONS.*/

  while(strncmp(str,"EIGEN",5)) fgets(str,256,ftext);


  for(i=0;i<(af->nnode);i++) /*DISPLACEMENTS.*/
  {
    data=fgetsbrk(ftext,&n);
    if(n!=9) return;

    ncode=strtol(*(data+1),NULL,10);
    if(ncode!=(af->nodes+i)->code) return;

    for(j=0;j<3;j++)
    {
      ddata=strtod(*(data+j+3),NULL);
      ddata+=(af->ninit+i)->d[j];
      (af->nodes+i)->d[j]=ddata;
      *(af->ddisp+6*i+j)=ddata;
    }
    for(j=3;j<6;j++)
    {
      ddata=strtod(*(data+j+3),NULL);
      *(af->ddisp+6*i+j)=ddata;
    }

    for(;n>0;n--) free(*(data+n-1));
    free(data);
  }

  return;
}/*bclngoutputtomemory*/
