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

#define NEIGEN 1  /*NUMBERS OF EIGEN VALUES TO BE OBTAINED.*/
#define SOLVER 1 /* 0:DEIGABGENERAL, 1:BISECSYLVESTER */
#define BISECEPS 1e-8/* BISECEPS: ACCURACY OF EIGENVALUE FOR BISECSYLVESTER */
/* BISECRIGHT: INITIAL UPPER BOUND FOR BISECSYLVESTER */
/*             IT IS ASSUMED THAT EIGENVALUE IS HIGHER THAN INVERSE OF BISECRIGHT */
/*             i.e. IF BISECRIGHT=1000.0, EIGENVALUE > 0.001 */
#define BISECRIGHT 1.0

int bclng001(struct arclmframe *af);/*ELASTIC BUCKLING ANALYSIS FOR ARCLM FRAME*/
/*-----------------------------------------------------------------*/
int bclng002(struct arclmframe *af);
                           /*BUCKLING CONDENSATION INTO INDIVIDUAL ELEMENTS*/
int bclng003(struct arclmframe *af,struct owire **multiwire,int nmultiwire);
                                     /*BUCKLING CONDENSATION INTO MULTIWIRE*/
/*-----------------------------------------------------------------*/
int bclng011(struct arclmframe *af,struct arclmframe *af0);
                                /*ELASTIC BUCKLING ANALYSIS WITH PRE-STRESS*/
int bclng101(struct arclmframe *af);
                                         /*ELASTO-PLASTIC BUCKLING ANALYSIS*/
/*-----------------------------------------------------------------*/
struct gcomponent *gdefine(unsigned int m,
						   unsigned int n,
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
void outputmodeII(double *gvct,FILE *fout,int nnode,struct onode *nodes,
                long int *loffs,int nmultinode);
void updatemode(struct arclmframe *af,double *gvct);
/*-----------------------------------------------------------------*/
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
/*-----------------------------------------------------------------*/
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

void bclngoutputtomemory(FILE *ftext,struct arclmframe *af);
double bclngoutputtomemoryII(FILE *ftext,struct arclmframe *af);        //ujok

void bisecsylvester(struct gcomponent *A,
					struct gcomponent *B,
					struct oconf *confs,
					long int N,long int NE,long int NV,
					double EPS,
					double *E,double **V);
void bisecgeneral(struct gcomponent *A,double factorA,
				  struct gcomponent *B,double factorB,
				  struct oconf *confs,
				  long int N,long int NE,double defsign,
				  double EPS,
				  double *E,double **V,
				  double BL, double BR);
double inversemethod(struct gcomponent *gmtx, struct oconf *confs, double *evct, int msize);
struct gcomponent *gcomponentadd2(struct gcomponent *mtx1,
								  struct gcomponent *mtx2,
                                  double factor,
								  int msize);
struct gcomponent *gcomponentadd3(struct gcomponent *mtx1,double factor1,
								  struct gcomponent *mtx2,double factor2,
								  int msize);
void definencr(struct arclmframe *af,double *ncr);

/*EXTERNAL PARAMETERS*/
extern FILE *globalfile; /*GLOBAL FILE.*/


int bclng001(struct arclmframe *af)
/*ELASTIC BUCKLING FOR ARCLM FRAME.*/
/*AFTER "ARCLM001".*/
{
  DWORDLONG memory0,memory1;
  FILE *fin,*fout;                              		 /*FILE 8 BYTES*/
  char dir[]=DIRECTORY;                        /*DATA DIRECTORY*/
  char s[800],string[400],fname[256];
  clock_t t0/*,t1,t2*/;

  int i,j,ii,jj;

  int nnode,nelem,nshell,nsect,nconstraint;
  long int msize;

  /*ARCLMFRAME*/
  struct osect *sects;
  struct onode *nodes;
  struct onode *ninit;
  struct owire *elems;
  struct oshell *shells;
  struct oconf *confs;
  struct memoryelem* melem;
  struct memoryshell* mshell;
  long int* constraintmain;
  double* ddisp, * iform;

  struct gcomponent ginit={0,0,0.0,NULL};
  struct gcomponent *gmtx1,*gmtx2,*epsgmtx1,*gmtx1cpy,*gcomp1;                    /*GLOBAL MATRIX*/
  double sign, determinant;
  long int nline;

	/*FOR READING ANALYSISDATA*/
	FILE *fdata;
	int nstr, pstr, readflag;
	char **data;
	char filename[256];
	char* dot;
	int neig = 1;
	int solver = 0;
	int method = 0;
	/*
	neig=NEIGEN;
	solver=SOLVER;
	*/

#if 0
	fin = fgetstofopenII(dir, "r", (wdraw.childs+1)->inpfile);
	if(fin==NULL)
	{
	  errormessage("ACCESS IMPOSSIBLE.");
	  return 0;
	}
	inputinitII(fin, &nnode, &nelem, &nshell, &nsect, &nconstraint); /*INPUT INITIAL.*/
#endif
#if 1
	nnode=af->nnode;
	nelem=af->nelem;
	nshell=af->nshell;
	nsect=af->nsect;
	nconstraint=af->nconstraint;
#endif

	msize=6*nnode;/*SIZE OF GLOBAL MATRIX.*/

	/*OUTPUT FILE NAME IS SAME AS INPUT*/
	strcpy(filename, (wdraw.childs+1)->inpfile);
	dot = strrchr(filename, '.');
	*dot = '\0';

	fdata = fgetstofopenII(dir, "r", "analysisdata.txt");

	if (fdata == NULL)
	{
		printf("couldn't open analysisdata.txt\n");
		if(MessageBox(NULL,"DEIGABGENERAL","SOLVER",MB_OKCANCEL)==IDOK)solver=0;
		else if(MessageBox(NULL,"BISECSYLVESTER","SOLVER",MB_OKCANCEL)==IDOK)solver=1;

		/*READ ARCLENGTH PARAMS*/
		getincrement((wmenu.childs+2)->hwnd,&neig,NULL);
	}
	else
	{
		readflag = 1;
		while (readflag)
		{
			data = fgetsbrk(fdata, &nstr);
			if (nstr == 0)
			{
				readflag = 0;
			}
			else
			{
				pstr = 0; /*POSITION IN "DATA".*/
				while ((nstr - pstr) > 0)
				{
					if (nstr - pstr == 1) /*POINTING LAST STRING.*/
					{
						pstr++;
					}
					else
					{
						if (!strcmp(*(data + pstr), "FILENAME"))
						{
							pstr++;
							if(strstr(*(data + pstr),filename)==NULL)
							{
								readflag = 0;
								if(MessageBox(NULL,"DEIGABGENERAL","SOLVER",MB_OKCANCEL)==IDOK)solver=0;
								else if(MessageBox(NULL,"BISECSYLVESTER","SOLVER",MB_OKCANCEL)==IDOK)solver=1;
								getincrement((wmenu.childs+2)->hwnd,&neig,NULL);
							}
							else
							{
								strcpy(filename, *(data + pstr));
							}
						}
						if (!strcmp(*(data + pstr), "NEIG"))
						{
							pstr++;
							neig = (int)strtol(*(data + pstr), NULL, 10);
						}
						if (!strcmp(*(data + pstr), "METHOD"))
						{
							pstr++;
							method = (int)strtol(*(data + pstr), NULL, 10);
						}
						else
						{
							pstr++;
						}
					}
				}
			}
		}
	}

  memory0=availablephysicalmemoryEx("INITIAL:");   /*MEMORY AVAILABLE*/

  snprintf(fname, sizeof(fname), "%s.%s", filename, "otp");
  fout = fopen(fname, "w");

  t0=clock();                                        /*CLOCK BEGIN.*/

#if 0
	/*MEMORY NOT ALLOCATED*/
	free(af->sects);
	free(af->nodes);
	free(af->elems);
	free(af->shells);
	free(af->confs);
	free(af->iform);
	free(af->ddisp);
	free(af->melem);
	free(af->mshell);
	free(af->constraintmain);

	sects = (struct osect*)malloc(nsect * sizeof(struct osect));
	nodes = (struct onode*)malloc(nnode * sizeof(struct onode));
	ninit = (struct onode*)malloc(nnode * sizeof(struct onode));
	elems = (struct owire*)malloc(nelem * sizeof(struct owire));
	shells = (struct oshell*)malloc(nshell * sizeof(struct oshell));
	confs = (struct oconf*)malloc(msize * sizeof(struct oconf));
	iform = (double*)malloc(msize * sizeof(double));
	ddisp = (double*)malloc(msize * sizeof(double));
	melem = (struct memoryelem*)malloc(nelem * sizeof(struct memoryelem));
	mshell = (struct memoryshell*)malloc(nshell * sizeof(struct memoryshell));
	constraintmain = (long int*)malloc(msize * sizeof(long int));

	af->sects = sects;
	af->nodes = nodes;
	af->ninit = ninit;
	af->elems = elems;
	af->shells = shells;
	af->confs = confs;
	af->iform = iform;
	af->ddisp = ddisp;
	af->melem = melem;
	af->mshell = mshell;
	af->constraintmain = constraintmain;

	inputtexttomemory(fin, af);
	fclose(fin);
#endif
#if 1
	/*MEMORY ALREADY ALLOCATED*/
	sects = af->sects;
	nodes = af->nodes;
	ninit = af->ninit;
	elems = af->elems;
	shells = af->shells;
	confs = af->confs;
	iform = af->iform;
	ddisp = af->ddisp;
	melem = af->melem;
	mshell = af->mshell;
	constraintmain = af->constraintmain;
#endif

  /*DIAGONALS OF GLOBAL MATRIX.*/
  gmtx1=(struct gcomponent *)malloc(msize*sizeof(struct gcomponent));
  gmtx2=(struct gcomponent *)malloc(msize*sizeof(struct gcomponent));
  if(gmtx1==NULL || gmtx2==NULL) return 0;
  for(i=0;i<msize;i++)
  {
	(gmtx1+i)->down=NULL;            /*GLOBAL MATRIX INITIALIZATION.*/
	(gmtx2+i)->down=NULL;
  }

  for(i=1;i<=msize;i++)             /*GLOBAL MATRIX INITIALIZATION.*/
  {
	ginit.m=(unsigned int)i;
	/*ginit.n=(unsigned int)i;*/
	*(gmtx1+(i-1))=ginit;
	*(gmtx2+(i-1))=ginit;
  }
  comps=msize; /*INITIAL COMPONENTS=DIAGONALS.*/

  GetAsyncKeyState(VK_LBUTTON);                  /*CLEAR KEY LEFT.*/
  GetAsyncKeyState(VK_RBUTTON);                  /*CLEAR KEY RIGHT.*/

  errormessage("BCLNG001:BUCKLING LINEAR.");
  //availablephysicalmemoryEx("REMAIN:");            /*MEMORY AVAILABLE*/
  //laptime("ASSEMBLING GLOBAL MATRIX.",t0);


	struct owire elem;
	struct oshell shell;
	int nnod;
	double* gforminit, * eforminit;                      /*DEFORMED COORDINATION OF ELEMENT*/
	double* gform, * eform;                      /*INITIAL COORDINATION OF ELEMENT*/
	double* edisp;
	double* einternal;                        /*INTERNAL FORCE OF ELEMENT*/

	double** Ke,** Kg,** drccos,** drccosinit;                           /*MATRIX*/
	double** T,** HPT;
	int* loffset;
	double area;


  /*ELEMENT STIFFNESS ASSEMBLAGE.*/
  if(method==0)
  {
	for (i = 1; i <= nelem; i++)
	{
		inputelem(elems,melem,i-1,&elem);
		nnod=elem.nnod;
		elem.sect=(elems+i-1)->sect;
		loffset = (int*)malloc(6 * nnod * sizeof(int));
		for (ii = 0; ii < nnod; ii++)
		{
			for (jj = 0; jj < 6; jj++)
			{
				*(loffset + (6 * ii + jj)) = 6 * (elem.node[ii]->loff) + jj;
			}
		}

		/*INITIAL CONFIGURATION*/
		for (ii = 0; ii < nnod; ii++)
		{
			inputnode(iform, elem.node[ii]);
		}
		drccosinit = directioncosine(elem.node[0]->d[0],
							elem.node[0]->d[1],
							elem.node[0]->d[2],
							elem.node[1]->d[0],
							elem.node[1]->d[1],
							elem.node[1]->d[2],
							elem.cangle);
		gforminit = extractdisplacement(elem, iform);
		eforminit = extractlocalcoord(gforminit,drccosinit,nnod);

		Ke = assememtx(elem);
		//Ke = modifyhinge(elem,Ke);             /*MODIFY MATRIX.*/

		/*DEFORMED CONFIDURATION*/
		for (ii = 0; ii < nnod; ii++)
		{
			inputnode(ddisp, elem.node[ii]);
		}
		gform = extractdisplacement(elem, ddisp);
		drccos = updatedrccos(drccosinit, gforminit, gform);
		eform = extractlocalcoord(gform,drccos,nnod);

		T = transmatrixIII(drccos, nnod);										/*[T].*/

		edisp = extractdeformation(eforminit, eform, nnod);                		/*{Ue}*/
		//einternal = matrixvector(Ke, edisp, 6 * nnod);          				/*{Fe}=[Ke]{Ue}.*/
		einternal=(double *)malloc(6*nnod*sizeof(double));
		for(ii=0;ii<nnod;ii++)
		{
		  for(jj=0;jj<6;jj++) *(einternal+6*ii+jj)=(melem+i)->stress[ii][jj];
		}

		HPT = (double**)malloc(6 * nnod * sizeof(double*));
		for (ii = 0; ii <6 * nnod; ii++)
		{
			*(HPT + ii) = (double*)malloc(6 * nnod * sizeof(double));
		}
		Kg = assemgmtxCR(eform, edisp, einternal, T, nnod);		/*[Kg]=[Kgr]+[Kgp]+[Kgm]*/
		symmetricmtx(Kg, 6*nnod);											/*SYMMETRIC TANGENTIAL MATRIX[Ksym].*/
		Ke = transformationIII(Ke, HPT, 6*nnod);							/*[Ke]=[Tt][Pt][Ht][K][H][P][T]*/

		assemgstiffnesswithDOFelimination(gmtx2, Ke, &elem, constraintmain); /*ASSEMBLAGE TANGENTIAL STIFFNESS MATRIX.*/
		assemgstiffnesswithDOFelimination(gmtx1, Kg, &elem, constraintmain); /*ASSEMBLAGE TANGENTIAL STIFFNESS MATRIX.*/

		freematrix(drccosinit, 3);
		freematrix(drccos, 3);
		freematrix(T, 6 * nnod);
		freematrix(Ke, 6 * nnod);
		freematrix(Kg, 6 * nnod);
		freematrix(HPT, 6 * nnod);

		free(einternal);

		free(eforminit);
		free(gforminit);
		free(eform);
		free(gform);
		free(edisp);
		free(loffset);

	}
	for (i = 1; i <= nshell; i++)
	{
		inputshell(shells, mshell, i - 1, &shell);
		nnod = shell.nnod;
		shell.sect = (shells + i - 1)->sect;                      /*READ SECTION DATA.*/
		loffset = (int*)malloc(6 * nnod * sizeof(int));
		for (ii = 0; ii < nnod; ii++)
		{
			for (jj = 0; jj < 6; jj++)
			{
				*(loffset + (6 * ii + jj)) = 6 * (shell.node[ii]->loff) + jj;
			}
		}

		/*INITIAL CONFIGURATION*/
		for (ii = 0; ii < nnod; ii++)
		{
			inputnode(iform, shell.node[ii]);
		}
		drccosinit = shelldrccos(shell);
		gforminit = extractshelldisplacement(shell, iform);                 /*{Xg}*/
		eforminit = extractlocalcoord(gforminit,drccosinit,nnod);        	/*{Xe}*/

		Ke = assemshellemtx(shell);                        /*[Ke].*/

		/*DEFORMED CONFIGFURATION*/
		for (ii = 0; ii < nnod; ii++)
		{
			inputnode(ddisp, shell.node[ii]);
		}
		drccos = shelldrccos(shell);
		gform = extractshelldisplacement(shell, ddisp);                     /*{Xg+Ug}*/
		eform = extractlocalcoord(gform,drccos,nnod); 			       	    /*{Xe+Ue}*/

		T = transmatrixIII(drccos, nnod);         							/*[T].*/

		edisp = extractdeformation(eforminit, eform, nnod);           		/*{Ue}*/
		//einternal = matrixvector(Ke, edisp, 6 * nnod);      				/*{Fe}=[Ke]{Ue}.*/
		einternal=(double *)malloc(6*nnod*sizeof(double));
		for(ii=0;ii<nnod;ii++)
		{
		  for(jj=0;jj<6;jj++) *(einternal+6*ii+jj)=(mshell+i)->stress[ii][jj];
		}

		HPT = (double**)malloc(6 * nnod * sizeof(double*));
		for (ii = 0; ii <6 * nnod; ii++)
		{
			*(HPT + ii) = (double*)malloc(6 * nnod * sizeof(double));
		}
		Kg = assemgmtxCR(eform, edisp, einternal, T, nnod);		/*[Kg]=[Kgr]+[Kgp]+[Kgm]*/
		symmetricmtx(Kg, 6*nnod);											/*SYMMETRIC TANGENTIAL MATRIX[Ksym].*/
		Ke = transformationIII(Ke, HPT, 6*nnod);							/*[Ke]=[Tt][Pt][Ht][K][H][P][T]*/

		assemgstiffnessIIwithDOFelimination(gmtx2, Ke, &shell, constraintmain); 	/*ASSEMBLAGE TANGENTIAL STIFFNESS MATRIX.*/
		assemgstiffnessIIwithDOFelimination(gmtx1, Kg, &shell, constraintmain); 	/*ASSEMBLAGE TANGENTIAL STIFFNESS MATRIX.*/

		freematrix(drccos, 3);
		freematrix(drccosinit, 3);
		freematrix(T, 6 * nnod);
		freematrix(Ke, 6 * nnod);
		freematrix(Kg, 6 * nnod);
		freematrix(HPT, 6 * nnod);

		free(einternal);

		free(eforminit);
		free(gforminit);
		free(eform);
		free(gform);
		free(edisp);
		free(loffset);
	}
  }
#if 0
  if(0)
  {
  	/*ELEMENT STIFFNESS & FORCE ASSEMBLAGE*/
	assemelem(elems, NULL, nelem, constraintmain, NULL, gmtx1, iform, ddisp);
	assemshell(shells, NULL, nshell, constraintmain, NULL, gmtx1, iform, ddisp);

	gmtx1cpy=copygcompmatrix(gmtx1,msize);

	nline = croutlu(gmtx1, confs, msize, &determinant, &sign, gcomp1);
	nline = forwardbackward(gmtx1, dup, confs, msize, gcomp1);

	for (ii = 0; ii < msize; ii++)/*UPDATE FORM FOR DIRECTIONAL DERIVATIVE*/
	{
		*(epsgvct + ii) = eps * (*(dup + ii)+*(givendisp + ii));
		*(epsddisp + ii) = *(ddisp + ii);
	}
	for (ii = 0; ii < msize; ii++)
	{
		*(epsgvct + ii) = *(epsgvct + *(constraintmain + ii));
	}
	updateform(epsddisp, epsgvct, nnode);

	/*ELEMENT STIFFNESS & FORCE ASSEMBLAGE*/
	assemelem(elems, NULL, nelem, constraintmain, NULL, epsgmtx1, iform, epsddisp);
	assemshell(shells, NULL, nshell, constraintmain, NULL, epsgmtx1, iform, epsddisp);

	gmtx2=gcomponentadd3(epsgmtx,1.0/eps,gmtx1,-1.0/eps,msize);
  }
#endif
#if 0

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
	T=transmatrix(drccos);           /*TRANSFORMATION MATRIX.*/

	einternal=(double *)malloc(12*sizeof(double));
	for(ii=0;ii<=1;ii++)                                /*STRESSES.*/
	{
	  for(jj=0;jj<6;jj++) *(einternal+6*ii+jj)=elem.stress[ii][jj];
	}


	Ke=assememtx(elem);            /*ELASTIC MATRIX OF ELEMENT.*/
	Ke=modifyhinge(elem,Ke);               /*MODIFY MATRIX.*/
	Ke=transformation(Ke,T);         /*[K]=[Tt][k][T]*/

	assemgstiffness(gmtx2,Ke,&elem);       /*ASSEMBLAGE ELASTIC.*/


	Kg=assemgmtx(elem,einternal);
	Kg=modifyhinge(elem,Kg);               /*MODIFY MATRIX.*/
	Kg=transformation(Kg,T);         /*[K]=[Tt][k][T]*/
	for(ii=0;ii<12;ii++)
	{
	  for(jj=0;jj<12;jj++) *(*(Kg+ii)+jj)=-*(*(Kg+ii)+jj);
	}

	assemgstiffness(gmtx1,Kg,&elem);     /*ASSEMBLAGE GEOMETRIC.*/

	freematrix(Ke, 12);
	freematrix(Kg, 12);
	free(einternal);
	freematrix(drccos, 3);
	freematrix(T, 12);
  }

#endif


  //laptime("GLOBAL MATRIX ASSEMBLED.",t0);

  //currentvalue("GLOBAL MATRIX:[K]",msize,neig,kmtx,NULL,NULL,NULL);
  //currentvalue("GLOBAL MATRIX:[G]",msize,neig,gmtx,NULL,NULL,NULL);

  /*EIGEN ANALYSIS START.       	*/
  /*([A]-eigen*[B])*{evct}=0    	*/
  /*                            	*/
  /*ONE-POINT BUCKLING ANALYSIS		*/
  /*([-Kg]-1/lambda*[Ke])*{evct}=0	*/
  /*A:[-Kg]							*/
  /*B:[Ke]							*/
  /*eigen=1/lambda     				*/
  /*lambda=1/eigen              	*/
  /*0<eigen                     	*/
  /*			                   	*/
  /*TWO-POINTS BUCKLING ANALYSIS	*/
  /*([Kt]-lambda*[-dKt])*{evct}=0	*/
  /*A:-([Kt]+[dKt])					*/
  /*B:[Kt]							*/
  /*eigen=(1-lambda)/lambda     	*/
  /*lambda=1/(1+eigen)     			*/
  /*-1<eigen                    	*/

  double **evct, *eigen;
  double eps=1.0E-16;
  double biseceps=1.0E-10;/*BISECEPS*/

  evct=(double **)malloc(neig*sizeof(double *));   /*GLOBAL VECTORS*/
  eigen=(double *)malloc(neig*sizeof(double));       /*EIGEN VALUES*/
  if(evct==NULL || eigen==NULL) return 0;
  for(i=0;i<neig;i++)
  {
	*(evct+i)=(double *)malloc(msize*sizeof(double));
	for(j=0;j<msize;j++) *(*(evct+i)+j)=0.0;
  }

  /*
  if(globalmessageflag==0||MessageBox(NULL,"DEIGABGENERAL","SOLVER",MB_OKCANCEL)==IDOK)
		deigabgeneral(gmtx,kmtx,confs,msize,neig,neig,eps,eigen,evct);
  else if(MessageBox(NULL,"BISECSYLVESTER","SOLVER",MB_OKCANCEL)==IDOK)
		bisecsylvester(gmtx,kmtx,confs,msize,neig,neig,biseceps,eigen,evct);
  */

  /*CALCULATE GENERAL EIGENVALUE PROBLEM FROM MAX.*/
  /*(factorA*[A]+eigen*factorB*[B]){evct}=0.*/
  bisecgeneral(gmtx1,-1.0,gmtx2,-1.0,confs,msize,neig,0,biseceps,eigen,evct,0.0,1.0);

  laptime("EIGEN COMPLETED.",t0);


  /*OUTPUT*/
  if(fout!=NULL)
  {
	  for(i=0;i<neig;i++)
	  {
		//outputmode(*(evct+i),fout,nnode,ninit);

		fprintf(fout,"LAP = %d MODE = %d GENERALIZED EIGENVALUE = %e STANDARD EIGENVALUE = %e ", af->nlaps, (i+1), *(eigen + i), 0.0);

		if (*(eigen + i) > 0.0)
		{
			*(eigen+i) = 1.0/(*(eigen + i));
			fprintf(fout, "LAMBDA = %e\n",*(eigen+i));
		}
		else
		{
			fprintf(fout, "ERROR:EIGEN VALUE NEGATIVE.\n");
		}

		for (ii = 0; ii < msize; ii++)
		{
			*(*(evct + i) + ii) = *(*(evct + i) + *(constraintmain + ii));
		}
		for (ii = 0; ii < nnode; ii++)
		{
			fprintf(fout,
			"%4ld %e %e %e %e %e %e\n",
			(nodes + ii)->code, *(*(evct + i) + 6*ii + 0),
			*(*(evct + i) + 6*ii + 1),
			*(*(evct + i) + 6*ii + 2),
			*(*(evct + i) + 6*ii + 3),
			*(*(evct + i) + 6*ii + 4),
			*(*(evct + i) + 6*ii + 5));
		}
	  }
  }
  fclose(fout);

  //af->nlaps=neig;
  af->eigenval=eigen;
  af->eigenvec=evct;

  /*EIGEN ANALYSIS END.*/

  //updatemode(af,*(evct+0)); /*FORMATION UPDATE.*/

  gfree(gmtx1,nnode); /*FREE GLOBAL MATRIX.*/
  gfree(gmtx2,nnode); /*FREE GLOBAL MATRIX.*/

  errormessage(" ");
  errormessage("COMPLETED.");

  memory1=availablephysicalmemoryEx("REMAIN:");
  sprintf(string,"CONSUMPTION:%ld[BYTES]",(memory0-memory1));
  errormessage(string);

#if 0
	/*STRESS UPDATE WITH BUCKLING STRESS.*/ //20210917.
  if(globalmessageflag && MessageBox(NULL,"Stress Update with Buckling Stress.","BCLNG011",MB_OKCANCEL)==IDOK)
  {
	for(i=0;i<af->nelem;i++)
	{
	  for(ii=0;ii<=1;ii++)
	  {
		for(jj=0;jj<6;jj++)
		{
		  (af->melem+((af->elems+i)->loff))->stress[ii][jj]*=*(eigen+0);
/*
		  if(ii==0 && jj==0)
		  {
			sprintf(string,"ELEM %ld Nz=%.5f",(af->elems+i)->code,(af->elems+i)->stress[ii][jj]);
			errormessage(string);
		  }
*/
		}
	  }
	}
   stressintofile(af); //extract
  }
#endif

  return 1;
}/*bclng001*/



/*-----------------------------------------------------------------*/
#if 1
int bclng002(struct arclmframe *af)
/*ELASTIC BUCKLING FOR ARCLM FRAME.*/
/*AFTER "ARCLM001".*/
{
  DWORD memory0,memory1;

  FILE /**fin,*/*fout,*feig,*frat;             /*FILE 8 BYTES*/
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
  double factor;
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
  fout=NULL; feig=NULL; frat=NULL;
  if(globalmessageflag)
  {
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
  strcat(string,".brat");
  frat=fopen(string,"w");
  }

  if(feig!=NULL) fprintf(feig,"ELEMENT EIGEN VALUES\n");

  t0=clock();                                        /*CLOCK BEGIN.*/

  nnode=af->nnode;                                 /*INPUT INITIAL.*/
  nelem=af->nelem;
  nsect=af->nsect;
  /*nreact=af->nreact;*/
  sprintf(string,"NODES=%d ELEMS=%d SECTS=%d",nnode,nelem,nsect);
  errormessage(string);
  if(fout!=NULL) fprintf(fout,"%s\n",string);

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
/*
  kcmtx=(double **)malloc(msize*sizeof(double *));
  gcmtx=(double **)malloc(msize*sizeof(double *));
  ktmtx=(double **)malloc(msize*sizeof(double *));
  gtmtx=(double **)malloc(msize*sizeof(double *));
*/
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
/*
	*(kcmtx+i)=(double *)malloc(msize*sizeof(double));
	*(gcmtx+i)=(double *)malloc(msize*sizeof(double));
	*(ktmtx+i)=(double *)malloc(msize*sizeof(double));
	*(gtmtx+i)=(double *)malloc(msize*sizeof(double));
*/
	*(kcmtx+i)=(double *)malloc(msize*sizeof(double));
	*(gcmtx+i)=(double *)malloc(msize*sizeof(double));
	*(ktmtx+i)=(double *)malloc(msize*sizeof(double));
	*(gtmtx+i)=(double *)malloc(msize*sizeof(double));
  }
//  gvct=(double **)malloc(neig*sizeof(double *));   /*GLOBAL VECTORS*/
  gvct=(double **)malloc(neig*sizeof(double *));   /*GLOBAL VECTORS*/
  eigen=(double *)malloc(neig*sizeof(double));       /*EIGEN VALUES*/
  if(gvct==NULL || eigen==NULL) return 0;
  for(i=0;i<neig;i++)
  {
//	*(gvct+i)=(double *)malloc(msize*sizeof(double));
 	*(gvct+i)=(double *)malloc(msize*sizeof(double));
 	for(j=0;j<msize;j++)  *(*(gvct+i)+j)=0.0;
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
if(fout!=NULL) fprintf(fout,"CONDENSATION BEGIN\n");

  for(i=0;i<nelem;i++)
  {
    if(fout!=NULL) fprintf(fout,"\nELEM %d OFFSET=%d\n",(af->elems+i)->code,i+1);

#if 1   /*FOR SYMMETRIC STRUCTURE:BEZIERSURFACEII*/
    if((af->elems+i)->node[0]->d[GX]<0 || (af->elems+i)->node[0]->d[GY]<0 ||
       (af->elems+i)->node[1]->d[GX]<0 || (af->elems+i)->node[1]->d[GY]<0 )
    {
      for(ii=0;ii<i;ii++)
      {
         if((abs(abs((af->elems+i)->node[0]->d[GX])-abs((af->elems+ii)->node[0]->d[GX]))<0.1) &&
            (abs(abs((af->elems+i)->node[0]->d[GY])-abs((af->elems+ii)->node[0]->d[GY]))<0.1) &&
            (abs(abs((af->elems+i)->node[1]->d[GX])-abs((af->elems+ii)->node[1]->d[GX]))<0.1) &&
            (abs(abs((af->elems+i)->node[1]->d[GY])-abs((af->elems+ii)->node[1]->d[GY]))<0.1)
           )
         {
            (af->elems+i)->srate[0]=(af->elems+ii)->srate[0];  /*Buckling safety ratio*/  //ujioka

            if(fout!=NULL)
            {
              fprintf(fout,"ELEM %d EIGEN VALUE=%12.5f",(af->elems+i)->code,1/((af->elems+i)->srate[0]));         //tsutsumi20171218
 	          fprintf(fout,"\n");
            }
            if(feig!=NULL) fprintf(feig,"ELEM %d EIGEN VALUE=%12.5f SAFETY=%12.5f\n",(af->elems+i)->code,1/((af->elems+i)->srate[0]),(af->elems+i)->srate[0]);                 //tsutsumi20171218
            if(frat!=NULL) fprintf(frat,"ELEM: %5d SECT: %4d %12.5f 0.00000 0.00000 0.00000\n",(af->elems+i)->code,(af->elems+i)->sect->code,(af->elems+i)->srate[0]);//ujioka negative-safety-ratio
            break;
         }
      }
    }
    else
#endif

    {
	inputelem(elems,af->melem,i,&elem);

    if(elem.stress[0][0]>=0)  factor=1.0;         /*COMPRESSION*/
    if(elem.stress[0][0]<0)   factor=-1.0;        /*TENSION*/
/*
    sprintf(string,"stress=%12.5f,factor=%12.0f",elem.stress[0][0],factor);
    MessageBox(NULL,string,"CONDENSE",MB_OK);
*/

	for(ii=0;ii<msize;ii++) /*DUPLICATE [Ke],[Kg]*/
	{
	  for(jj=0;jj<msize;jj++)
	  {
		gread(kmtx,ii+1,jj+1,&gdata);
		*(*(kcmtx+ii)+jj)=gdata;
		if(ii!=jj) *(*(kcmtx+jj)+ii)=gdata;
		gread(gmtx,ii+1,jj+1,&gdata);
		*(*(gcmtx+ii)+jj)=factor*gdata;
		if(ii!=jj) *(*(gcmtx+jj)+ii)=factor*gdata;
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
//MessageBox(NULL,string,"CONDENSE",MB_OK);               //tsutsumi20171218

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
				  *(*(gtmtx+ii)+jj)=(*(*(gcmtx+ii)+jj));    //tsutsumi
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
    if(globalmessageflag) currentpivot(j+1,msize);
    else currentpivot(i+1,nelem);
	} /*END CONDENSATION*/


if(fout!=NULL)
{
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
if(fout!=NULL)
{
fprintf(fout,"{CONF} :");
for(ii=0;ii<12;ii++) fprintf(fout," %3d",CF[ii]);
fprintf(fout,"\n");
}

	/*deigqr(AA,BB,12,12,12,12,eps,WW,EE,VV,CF);*/
	deigqrcf(AA,BB,12,12,12,12,eps,WW,EE,VV,CF);

if(globalmessageflag)
{
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
  if(!CF[ii]) fprintf(globalfile," %18.8f\n",factor*EE[ii]);
}
fprintf(globalfile,"\n");
}

	for(ii=0;ii<12;ii++)
	{
	  if(!CF[ii])
	  {
		mm=ii;
		break;
	  }
	}
    if(fout!=NULL)
    {
//	if(EE[mm]==0.0) MessageBox(NULL,"EIGEN VALUE=0.0","CONDENSE",MB_OK);
//	fprintf(fout,"ELEM %d EIGEN VALUE=%12.5f",(af->elems+i)->code,1/fabs(EE[mm]));
	fprintf(fout,"ELEM %d EIGEN VALUE=%12.5f",(af->elems+i)->code,factor/EE[mm]);         //tsutsumi20171218
	/*fprintf(fout," LINE=%d",mm);*/
	fprintf(fout,"\n");
    if(feig!=NULL && frat!=NULL)
    {
//	fprintf(feig,"ELEM %d EIGEN VALUE=%12.5f SAFETY=%12.5f\n",(af->elems+i)->code,1/fabs(EE[mm]),fabs(EE[mm]));
	fprintf(feig,"ELEM %d EIGEN VALUE=%12.5f SAFETY=%12.5f\n",(af->elems+i)->code,factor/EE[mm],factor*EE[mm]);                 //tsutsumi20171218
//	fprintf(frat,"ELEM: %5d SECT: %4d %12.5f 0.00000 0.00000 0.00000\n",(af->elems+i)->code,(af->elems+i)->sect->code,fabs(EE[mm]));
	fprintf(frat,"ELEM: %5d SECT: %4d %12.5f 0.00000 0.00000 0.00000\n",(af->elems+i)->code,(af->elems+i)->sect->code,factor*EE[mm]);//ujioka negative-safety-ratio
    }
    }

    (af->elems+i)->srate[0]=factor*EE[mm];  /*Buckling safety ratio*/  //ujioka
    }
    if(globalmessageflag)
    {
      setlaps((wmenu.childs+2)->hwnd,i+1,nelem);
    }
  }
  laptime("EIGEN COMPLETED.",t0);

  af->nlaps=neig;
  af->eigenval=eigen;
//  af->eigenvec=gvct;

  gfree(kmtx,nnode); /*FREE GLOBAL MATRIX.*/
  gfree(gmtx,nnode); /*FREE GLOBAL MATRIX.*/

  for(ii=0;ii<msize;ii++) free(*(kcmtx+ii));
  free(kcmtx);
  for(ii=0;ii<msize;ii++) free(*(gcmtx+ii));
  free(gcmtx);
  for(ii=0;ii<msize;ii++) free(*(ktmtx+ii));
  free(ktmtx);
  for(ii=0;ii<msize;ii++) free(*(gtmtx+ii));
  free(gtmtx);
  for(ii=0;ii<neig;ii++) free(*(gvct+ii));
  free(gvct);

  errormessage(" ");
  errormessage("COMPLETED.");
  if(globalmessageflag) fprintf(fout,"COMPLETED.\n");

  if(fout!=NULL) fclose(fout);
  if(feig!=NULL) fclose(feig);
  if(frat!=NULL) fclose(frat);

  memory1=availablephysicalmemory("REMAIN:");
  sprintf(string,"CONSUMPTION:%ld[BYTES]",(memory0-memory1));
  errormessage(string);

  /*Buckling safety ratio*/  //ujioka
  (wdraw.childs+1)->vparam.vflag.ev.srcancolor=1;
  (wdraw.childs+1)->vparam.vflag.ev.srcanrate=1;

  return 1;
}/*bclng002*/
#endif


/*-----------------------------------------------------------------*/
#if 0  /*float*/
int bclng002(struct arclmframe *af)
/*ELASTIC BUCKLING FOR ARCLM FRAME.*/
/*AFTER "ARCLM001".*/
{
  DWORD memory0,memory1;

  FILE /**fin,*/*fout,*feig,*frat;             /*FILE 8 BYTES*/
  /*FILE *felem,*fdisp,*freact;*/
  /*char dir[]=DIRECTORY;*/                        /*DATA DIRECTORY*/
  char s[80],string[1024];
  int i,j,ii,jj,mm,nn;
  int nnode,nelem,nsect/*,nreact*/;
  long int /*fsize,*/loff,ioff,joff,msize;
  /*long int time;*/
  struct gcomponent ginit={0,0,0.0,NULL};
  struct gcomponent *kmtx,*gmtx;                    /*GLOBAL MATRIX*/
  float **kcmtx,**gcmtx,**ktmtx,**gtmtx;        /*CONDENSED MATRIX*/
  float **gvct;                                    /*GLOBAL VECTOR*/
  double **drccos,**tmatrix,**estiff,*estress;
  /*double determinant,data;*/
  double gdata;
  double factor;

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
  strcat(string,".brat");
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
/*
  kcmtx=(double **)malloc(msize*sizeof(double *));
  gcmtx=(double **)malloc(msize*sizeof(double *));
  ktmtx=(double **)malloc(msize*sizeof(double *));
  gtmtx=(double **)malloc(msize*sizeof(double *));
*/
  kcmtx=(float **)malloc(msize*sizeof(float *));
  gcmtx=(float **)malloc(msize*sizeof(float *));
  ktmtx=(float **)malloc(msize*sizeof(float *));
  gtmtx=(float **)malloc(msize*sizeof(float *));
  if(kmtx==NULL || gmtx==NULL) return 0;
  if(kcmtx==NULL || gcmtx==NULL) return 0;
  if(ktmtx==NULL || gtmtx==NULL) return 0;
  for(i=0;i<msize;i++)
  {
	(kmtx+i)->down=NULL;            /*GLOBAL MATRIX INITIALIZATION.*/
	(gmtx+i)->down=NULL;
	/*(kcmtx+i)->down=NULL;*/
	/*(gcmtx+i)->down=NULL;*/
/*
	*(kcmtx+i)=(double *)malloc(msize*sizeof(double));
	*(gcmtx+i)=(double *)malloc(msize*sizeof(double));
	*(ktmtx+i)=(double *)malloc(msize*sizeof(double));
	*(gtmtx+i)=(double *)malloc(msize*sizeof(double));
*/
	*(kcmtx+i)=(float *)malloc(msize*sizeof(float));
	*(gcmtx+i)=(float *)malloc(msize*sizeof(float));
	*(ktmtx+i)=(float *)malloc(msize*sizeof(float));
	*(gtmtx+i)=(float *)malloc(msize*sizeof(float));
  }
//  gvct=(double **)malloc(neig*sizeof(double *));   /*GLOBAL VECTORS*/
  gvct=(float **)malloc(neig*sizeof(float *));   /*GLOBAL VECTORS*/
  eigen=(double *)malloc(neig*sizeof(double));       /*EIGEN VALUES*/
  if(gvct==NULL || eigen==NULL) return 0;
  for(i=0;i<neig;i++)
  {
//	*(gvct+i)=(double *)malloc(msize*sizeof(double));
 	*(gvct+i)=(float *)malloc(msize*sizeof(float));
 	for(j=0;j<msize;j++)  *(*(gvct+i)+j)=0.0;
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
#if 0   /*FOR SYMMETRIC STRUCTURE:BEZIERSURFACEII*/
    if((af->elems+i)->node[0]->d[GX]<0 || (af->elems+i)->node[0]->d[GY]<0 ||
       (af->elems+i)->node[1]->d[GX]<0 || (af->elems+i)->node[1]->d[GY]<0 )
    {
      for(ii=0;ii<i;ii++)
      {
         if((abs((af->elems+i)->node[0]->d[GX])==abs((af->elems+ii)->node[0]->d[GX])) &&
            (abs((af->elems+i)->node[0]->d[GY])==abs((af->elems+ii)->node[0]->d[GY])) &&
            (abs((af->elems+i)->node[1]->d[GX])==abs((af->elems+ii)->node[1]->d[GX])) &&
            (abs((af->elems+i)->node[1]->d[GY])==abs((af->elems+ii)->node[1]->d[GY]))
           )
         {
            (af->elems+i)->srate[0]=(af->elems+ii)->srate[0];  /*Buckling safety ratio*/  //ujioka

            fprintf(fout,"ELEM %d EIGEN VALUE=%12.5f",(af->elems+i)->code,1/((af->elems+i)->srate[0]));         //tsutsumi20171218
 	        fprintf(fout,"\n");
            fprintf(feig,"ELEM %d EIGEN VALUE=%12.5f SAFETY=%12.5f\n",(af->elems+i)->code,1/((af->elems+i)->srate[0]),(af->elems+i)->srate[0]);                 //tsutsumi20171218
            fprintf(frat,"ELEM: %5d SECT: %4d %12.5f 0.00000 0.00000 0.00000\n",(af->elems+i)->code,(af->elems+i)->sect->code,(af->elems+i)->srate[0]);//ujioka negative-safety-ratio
            break;
         }
      }
    }
    else
#endif
    {
	inputelem(elems,af->melem,i,&elem);

    if(elem.stress[0][0]>=0)  factor=1.0;         /*COMPRESSION*/
    if(elem.stress[0][0]<0)   factor=-1.0;        /*TENSION*/
/*
    sprintf(string,"stress=%12.5f,factor=%12.0f",elem.stress[0][0],factor);
    MessageBox(NULL,string,"CONDENSE",MB_OK);
*/

	for(ii=0;ii<msize;ii++) /*DUPLICATE [Ke],[Kg]*/
	{
	  for(jj=0;jj<msize;jj++)
	  {
		gread(kmtx,ii+1,jj+1,&gdata);
		*(*(kcmtx+ii)+jj)=gdata;
		if(ii!=jj) *(*(kcmtx+jj)+ii)=gdata;
		gread(gmtx,ii+1,jj+1,&gdata);
		*(*(gcmtx+ii)+jj)=factor*gdata;
		if(ii!=jj) *(*(gcmtx+jj)+ii)=factor*gdata;
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
//MessageBox(NULL,string,"CONDENSE",MB_OK);               //tsutsumi20171218

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
				  *(*(gtmtx+ii)+jj)=(*(*(gcmtx+ii)+jj));    //tsutsumi
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
    if(globalmessageflag) currentpivot(j+1,msize);
    else currentpivot(i+1,nelem);
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
//	if(EE[mm]==0.0) MessageBox(NULL,"EIGEN VALUE=0.0","CONDENSE",MB_OK);
//	fprintf(fout,"ELEM %d EIGEN VALUE=%12.5f",(af->elems+i)->code,1/fabs(EE[mm]));
	fprintf(fout,"ELEM %d EIGEN VALUE=%12.5f",(af->elems+i)->code,factor/EE[mm]);         //tsutsumi20171218
	/*fprintf(fout," LINE=%d",mm);*/
	fprintf(fout,"\n");

//	fprintf(feig,"ELEM %d EIGEN VALUE=%12.5f SAFETY=%12.5f\n",(af->elems+i)->code,1/fabs(EE[mm]),fabs(EE[mm]));
	fprintf(feig,"ELEM %d EIGEN VALUE=%12.5f SAFETY=%12.5f\n",(af->elems+i)->code,factor/EE[mm],factor*EE[mm]);                 //tsutsumi20171218
//	fprintf(frat,"ELEM: %5d SECT: %4d %12.5f 0.00000 0.00000 0.00000\n",(af->elems+i)->code,(af->elems+i)->sect->code,fabs(EE[mm]));
	fprintf(frat,"ELEM: %5d SECT: %4d %12.5f 0.00000 0.00000 0.00000\n",(af->elems+i)->code,(af->elems+i)->sect->code,factor*EE[mm]);//ujioka negative-safety-ratio

    (af->elems+i)->srate[0]=factor*EE[mm];  /*Buckling safety ratio*/  //ujioka
    }
    if(globalmessageflag)
    {
      setlaps((wmenu.childs+2)->hwnd,i+1,nelem);
    }
  }
  laptime("EIGEN COMPLETED.",t0);

  af->nlaps=neig;
  af->eigenval=eigen;
//  af->eigenvec=gvct;

  gfree(kmtx,nnode); /*FREE GLOBAL MATRIX.*/
  gfree(gmtx,nnode); /*FREE GLOBAL MATRIX.*/

  errormessage(" ");
  errormessage("COMPLETED.");
  fprintf(fout,"COMPLETED.\n");

  fclose(fout);
  fclose(feig);
  fclose(frat);

  memory1=availablephysicalmemory("REMAIN:");
  sprintf(string,"CONSUMPTION:%ld[BYTES]",(memory0-memory1));
  errormessage(string);

  /*Buckling safety ratio*/  //ujioka
  (wdraw.childs+1)->vparam.vflag.ev.srcancolor=1;
  (wdraw.childs+1)->vparam.vflag.ev.srcanrate=1;

  return 1;
}/*bclng002*/
#endif


/*-----------------------------------------------------------------*/
#if 0
int bclng003(struct arclmframe *af,struct owire **multiwire,int nmultiwire)
/*ELASTIC BUCKLING FOR ARCLM FRAME.*/
/*AFTER "ARCLM001".*/
/*BUCKLING CONDENSATION INTO MULTIWIRE*/  /*UJIOKA*/
{
  DWORD memory0,memory1;

  FILE /**fin,*/*fout,*feig/*,*frat*/;                   /*FILE 8 BYTES*/
  /*FILE *felem,*fdisp,*freact;*/
  /*char dir[]=DIRECTORY;*/                        /*DATA DIRECTORY*/
  char s[4096],string[4096];
  int i,j,ii,jj,mm,nn;
  int nnode,nelem,nsect/*,nreact*/;
  long int /*fsize,*/loff,ioff,joff,msize;
  long int *loffs,*moffs,*noffs,k,kk;
  int nmultinode;
  /*long int time;*/
  struct gcomponent ginit={0,0,0.0,NULL};
  struct gcomponent *kmtx,*gmtx;                    /*GLOBAL MATRIX*/
  struct gcomponent *kmtx2,*gmtx2;                  /*GLOBAL MATRIX*/
  struct gcomponent *ge,*g,*p;                      /*GLOBAL MATRIX*/
//  double **kcmtx,**gcmtx,**ktmtx,**gtmtx;        /*CONDENSED MATRIX*/
//  double **kcmtx2,**gcmtx2;
  float **kcmtx,**gcmtx,**ktmtx,**gtmtx;        /*CONDENSED MATRIX*/
  float **kcmtx2,**gcmtx2;
  double **gvct,**gvct2;                            /*GLOBAL VECTOR*/
  double **drccos,**tmatrix,**estiff,*estress;
  /*double determinant,data;*/
  double gdata,factor;
  clock_t t0/*,t1,t2*/;

  /*struct osect *sects;*/
  struct onode /**nodes,*/*ninit;
  struct owire elem;
  struct owire *elems;
  struct oconf *confs,*confs2;
  struct oconf *initialconfs;

  long int neig;
//  double AA[12][12],BB[12][12],WW[12][12],EE[12],VV[12][12];
//  signed char CF[12];
//  double **AA,**BB,**WW,*EE,**VV;
//  signed char *CF;
  double eps=1.0E-16,*eigen,biseceps;

  /*FOR SORT NODES*/
  int flag,tmp;
  struct owire *pe;
  char str[80];
  nmultinode=0;
  /*FOR SORT NODES*/

  memory0=availablephysicalmemory("INITIAL:");   /*MEMORY AVAILABLE*/

  fout=fgetstofopen("\0","w",ID_OUTPUTFILE);          /*OUTPUT FILE*/
  //fout=fgetstofopen("\0","a",ID_OUTPUTFILE);          /*OUTPUT FILE*/
  globalfile=fout;

  GetDlgItemText((wmenu.childs+2)->hwnd,ID_OUTPUTFILE,s,80);
  /*GetDlgItemText((wmenu.childs+2)->hwnd,ID_SECTIONFILE,s,80);*/

  biseceps=BISECEPS;

  nn=strcspn(s,".");
  strncpy(string,s,nn);
  string[nn]='\0';
  strcat(string,".eig");
  feig=fopen(string,"w");
/*
  sprintf(string,"\0");
  strncpy(string,s,nn);
  string[nn]='\0';
  strcat(string,".rat");
  frat=fopen(string,"w");
*/
  fprintf(fout,"CONDENSED MULTIELEM EIGEN VALUES\n");
//  fprintf(feig,"CONDENSED MULTIELEM EIGEN VALUES\n");

  t0=clock();                                        /*CLOCK BEGIN.*/

  nnode=af->nnode;                                 /*INPUT INITIAL.*/
  nelem=af->nelem;
  nsect=af->nsect;
  /*nreact=af->nreact;*/
  sprintf(string,"NODES=%d ELEMS=%d SECTS=%d",nnode,nelem,nsect);
  errormessage(string);
  fprintf(fout,"%s\n",string);
  fprintf(feig,"%s\n",string);

  msize=6*nnode;                           /*SIZE OF GLOBAL MATRIX.*/
  neig=NEIGEN;
//  neig=12;

  kmtx=(struct gcomponent *)          /*DIAGONALS OF GLOBAL MATRIX.*/
		malloc(msize*sizeof(struct gcomponent));
  gmtx=(struct gcomponent *)
		malloc(msize*sizeof(struct gcomponent));
  /*kcmtx=(struct gcomponent *)
		 malloc(msize*sizeof(struct gcomponent));*/
  /*gcmtx=(struct gcomponent *)
		 malloc(msize*sizeof(struct gcomponent));*/
/*
  kcmtx=(double **)malloc(msize*sizeof(double *));
  gcmtx=(double **)malloc(msize*sizeof(double *));
  kcmtx2=(double **)malloc(msize*sizeof(double *));
  gcmtx2=(double **)malloc(msize*sizeof(double *));
  ktmtx=(double **)malloc(msize*sizeof(double *));
  gtmtx=(double **)malloc(msize*sizeof(double *));
*/
  kcmtx=(float **)malloc(msize*sizeof(float *));
  gcmtx=(float **)malloc(msize*sizeof(float *));
  kcmtx2=(float **)malloc(msize*sizeof(float *));
  gcmtx2=(float **)malloc(msize*sizeof(float *));
  ktmtx=(float **)malloc(msize*sizeof(float *));
  gtmtx=(float **)malloc(msize*sizeof(float *));

  if(kmtx==NULL || gmtx==NULL) return 0;
  if(kcmtx==NULL || gcmtx==NULL) return 0;
  if(ktmtx==NULL || gtmtx==NULL) return 0;
  for(i=0;i<msize;i++)
  {
	(kmtx+i)->down=NULL;            /*GLOBAL MATRIX INITIALIZATION.*/
	(gmtx+i)->down=NULL;
	/*(kcmtx+i)->down=NULL;*/
	/*(gcmtx+i)->down=NULL;*/
/*
	*(kcmtx+i)=(double *)malloc(msize*sizeof(double));
	*(gcmtx+i)=(double *)malloc(msize*sizeof(double));
	*(kcmtx2+i)=(double *)malloc(msize*sizeof(double));
	*(gcmtx2+i)=(double *)malloc(msize*sizeof(double));
	*(ktmtx+i)=(double *)malloc(msize*sizeof(double));
	*(gtmtx+i)=(double *)malloc(msize*sizeof(double));
*/
	*(kcmtx+i)=(float *)malloc(msize*sizeof(float));
	*(gcmtx+i)=(float *)malloc(msize*sizeof(float));
	*(kcmtx2+i)=(float *)malloc(msize*sizeof(float));
	*(gcmtx2+i)=(float *)malloc(msize*sizeof(float));
	*(ktmtx+i)=(float *)malloc(msize*sizeof(float));
	*(gtmtx+i)=(float *)malloc(msize*sizeof(float));
  }

  gvct=(double **)malloc(neig*sizeof(double *));   /*GLOBAL VECTORS*/
  gvct2=(double **)malloc(neig*sizeof(double *));   /*GLOBAL VECTORS*/
  eigen=(double *)malloc(neig*sizeof(double));       /*EIGEN VALUES*/
  if(gvct==NULL || gvct2==NULL || eigen==NULL) return 0;
  for(i=0;i<neig;i++)
  {
	*(gvct+i)=(double *)malloc(msize*sizeof(double));
    *(gvct2+i)=(double *)malloc(msize*sizeof(double));
	for(j=0;j<msize;j++)
    {
       *(*(gvct+i)+j)=0.0;
       *(*(gvct2+i)+j)=0.0;
    }
  }

  /*fdisp=af->fdisp;*/                 /*DISPLACEMENT:6 DIRECTIONS.*/
  /*felem=af->felem;*/              /*CODE,12 BOUNDARIES,12 STRESS.*/
  /*freact=af->freact;*/                           /*REACTION FILE.*/

  /*sects=af->sects;*/
  /*nodes=af->nodes;*/
  ninit=af->ninit;
  elems=af->elems;
  confs=af->confs;
//  confs2=af->confs;        /*WRONG*/
  confs2=(struct oconf *)malloc(msize*sizeof(struct oconf));

  initialconfs=(struct oconf *)malloc(msize*sizeof(struct oconf));
  for(i=0;i<msize;i++)
  {
    (initialconfs+i)->iconf=(confs+i)->iconf;
    (initialconfs+i)->value=(confs+i)->value;
  }

  GetAsyncKeyState(VK_LBUTTON);                   /*CLEAR KEY LEFT.*/
  GetAsyncKeyState(VK_RBUTTON);                  /*CLEAR KEY RIGHT.*/

  errormessage("BCLNG003:BUCKLING LINEAR.");
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


	assemgstiffness(kmtx,estiff,&elem);       /*ASSEMBLAGE ELASTIC.*/

	for(ii=0;ii<=11;ii++) free(*(estiff+ii));
	free(estiff);

	for(ii=0;ii<=1;ii++)                                /*STRESSES.*/
	{
	  for(jj=0;jj<6;jj++) *(estress+6*ii+jj)=elem.stress[ii][jj];
	}

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


  /*SORT NODES IN MULTIWIRE*/
  loffs=(long int *)malloc(sizeof(long int));/*loff of multiwire*/

  for(i=0;i<nnode;i++)
  {
    flag=0;
    for(j=0;j<nmultiwire;j++)
    {
      pe=*(multiwire+j);

      for(k=0;k<2;k++)
      {
        if((af->nodes+i)->code==((pe->node[k])->code)) flag=1;
      }
      if(flag) break;
    }

    if(flag)
    {
      if(nmultinode)
      {
        nmultinode+=1;
        loffs=(long int *)realloc(loffs,nmultinode*sizeof(long int));
        *(loffs+nmultinode-1)=(af->nodes+i)->loff;
      }
      else     //initial
      {
        *(loffs+0)=(af->nodes+i)->loff;
        nmultinode=1;
      }
    }
  }
/*
  for(i=0;i<nmultinode;i++)
  {
   sprintf(string,"loffs[%ld]=%ld\n",i,*(loffs+i));
   errormessage(string);
  }
*/

  moffs=(long int *)malloc(nnode*sizeof(long int));
  noffs=(long int *)malloc(nnode*sizeof(long int));

  /*CONDENSATION*/
  sprintf(string,"\nBCLNG003:CONDENSATION BEGIN\n");
  errormessage(string);

  sprintf(string,"CONDENSED NODES=%ld ELEMS=%ld",nmultinode,nmultiwire);
  errormessage(string);
  fprintf(fout,"%s\n",string);

  for(i=0;i<nmultinode;i++)
  {
    sprintf(string,"NODE  %ld",(af->nodes+(*(loffs+i)))->code);
//    errormessage(string);
//    fprintf(fout,"%s\n",string);
//    fprintf(feig,"%s\n",string);
  }

  for(i=0;i<nmultiwire;i++)
  {
    sprintf(string,"ELEM  %ld",(*(multiwire+i))->code);
//    errormessage(string);
    fprintf(fout,"%s\n",string);
//    fprintf(feig,"%s\n",string);
  }
  fprintf(fout,"\n");
//  fprintf(feig,"\n");

    factor=1.0;
#if 1       /*FOR TENSION COLUMN BUCKLING LOAD*/
	inputelem(elems,af->melem,(*(multiwire+0))->loff,&elem);
    if(elem.stress[0][0]<0)   factor=-1.0;        /*TENSION*/
#endif

  /* for(i=0;i<nelem;i++) */
  {
	for(ii=0;ii<msize;ii++) /*DUPLICATE [Ke],[Kg]*/
	{
	  for(jj=0;jj<msize;jj++)
	  {
		gread(kmtx,ii+1,jj+1,&gdata);
		*(*(kcmtx+ii)+jj)=gdata;
		if(ii!=jj) *(*(kcmtx+jj)+ii)=gdata;
		gread(gmtx,ii+1,jj+1,&gdata);
		*(*(gcmtx+ii)+jj)=factor*gdata;
		if(ii!=jj) *(*(gcmtx+jj)+ii)=factor*gdata;
	  }
	}


/*
    for(i=0;i<6*nnode;i++)
    {
      sprintf(string,"befor:kcmtx[%ld][%ld]=%ld\n",i+1,i+1,*(*(kcmtx+i)+i));
      errormessage(string);
    }
*/
/*
    for(i=0;i<6;i++)
    {

      sprintf(string,"before:kcmtx[%ld[%ld]=%ld\n",6**(loffs+0)+i+1,6**(loffs+0)+i+1,
      				*(*(kcmtx+6**(loffs+0)+i)+6**(loffs+0)+i));
      errormessage(string);
      fprintf(fout,"%s\n",string);
    }
*/

    /*EXCHANGING LINES FOR CONDENSATION*/
    i=0;
    ii=0;
    for(i=0;i<nnode;i++)
    {
      *(moffs+i)=i;
      *(noffs+i)=i;
    }

    for(i=0;i<nnode;i++)
    {
      for(ii=0;ii<nmultinode;ii++)
      {
        if(i==*(loffs+ii))
        {
        tmp=*(noffs+i);
        *(noffs+i)=*(noffs+ii);
        *(noffs+ii)=tmp;
        }
      }
    }

/*
    for(i=0;i<nnode;i++)
    {
      sprintf(string,"moffs[%ld]=%ld\n",i,*(moffs+i));
      errormessage(string);
      fprintf(fout,"%s\n",string);
    }
    for(i=0;i<nnode;i++)
    {
      sprintf(string,"noffs[%ld]=%ld\n",i,*(noffs+i));
      errormessage(string);
      fprintf(fout,"%s\n",string);
    }
*/
    for(i=0;i<msize;i++)
    {
      (confs2+i)->iconf=0;
      (confs2+i)->value=0.0;
    }

    exchangelinesIIfloat(kcmtx,gcmtx,confs,kcmtx2,gcmtx2,confs2,
                    moffs,noffs,af->nnode,0);
/*
    for(i=0;i<6*nnode;i++)
    {
      sprintf(string,"iconf[%ld]=%ld\n",i+1,(confs+i)->iconf);
      errormessage(string);
      fprintf(fout,"%s\n",string);
    }
*/
	for(ii=0;ii<msize;ii++)
	{
	  for(jj=0;jj<msize;jj++)
	  {
		*(*(kcmtx+ii)+jj)=*(*(kcmtx2+ii)+jj);
		if(ii!=jj) *(*(kcmtx+jj)+ii)=*(*(kcmtx2+ii)+jj);
		*(*(gcmtx+ii)+jj)=*(*(gcmtx2+ii)+jj);
		if(ii!=jj) *(*(gcmtx+jj)+ii)=*(*(gcmtx2+ii)+jj);
	  }
      *(confs+ii)=*(confs2+ii);
	}

    laptime("LINES EXCHANGED.",t0);

/*
    for(i=0;i<nnode*6;i++)
    {
      sprintf(string,"kcmtx[%ld][%ld]=%ld\n",i+1,i+1,*(*(kcmtx+i)+i));
      errormessage(string);
      fprintf(fout,"%s\n",string);
    }

    for(i=0;i<6*nnode;i++)
    {
      sprintf(string,"iconf[%ld]=%ld\n",i+1,(confs+i)->iconf);
      errormessage(string);
    }

    for(i=0;i<6*nnode;i++)
    {
      sprintf(string,"iconf[%ld]=%ld\n",i+1,(confs+i)->iconf);
      errormessage(string);
      fprintf(fout,"%s\n",string);
    }
*/

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

//sprintf(string,"ELEM ORDER=%d",i+1);
//MessageBox(NULL,string,"CONDENSE",MB_OK);               //tsutsumi20171218

    laptime("MATRIX ELIMINATION BEGIN.",t0);

	/*CONDENSATION*/ /*WITH CONSIDERING CONF*/
	for(j=/*0*/6*nmultinode;j<msize;j++)/*j:縮約で削除する行、列*/
	{
	  if(j>=msize) break;
/*
      sprintf(string,"kcmtx[%ld][%ld]=%ld\n",j+1,j+1,
      				*(*(kcmtx+j)+j));
      errormessage(string);
      fprintf(fout,"%s\n",string);
*/
	  if(!(confs+j)->iconf)
	  {
		if(*(*(kcmtx+j)+j)==0.0)
		{
		  sprintf(string,"INSTABLE TERMINATION K%d%d=%9.3f",j+1,j+1,*(*(kcmtx+j)+j));
		  MessageBox(NULL,string,"CONDENSE",MB_OK);
          fclose(fout);
		  return 0;
		}
/*
sprintf(string,"i=%d j=%d",i,j);
MessageBox(NULL,string,"CONDENSE",MB_OK);
*/
        ii=0;

		while(ii<msize)
		{
          if(ii==6*(nmultinode)) ii=j;
/*
sprintf(string,"i=%d j=%d ii=%d",i,j,ii);
MessageBox(NULL,string,"CONDENSE",MB_OK);
*/
		  if(!(confs+ii)->iconf)
		  {
            jj=0;
			while(jj<msize)
			{
              if(jj==6*(nmultinode)) jj=j;
/*
sprintf(string,"i=%d j=%d ii=%d jj=%d",i,j,ii,jj);
MessageBox(NULL,string,"CONDENSE",MB_OK);
*/
			  if(!(confs+jj)->iconf)
			  /*if(jj!=j)*/
			  {
              /*縮約（弾性剛性マトリクス）*/
				*(*(ktmtx+ii)+jj)=(*(*(kcmtx+ii)+jj))-(*(*(kcmtx+ii)+j))*(*(*(kcmtx+j)+jj))/(*(*(kcmtx+j)+j));

              /*縮約（幾何剛性マトリクス）*/
				/*
				*(*(gtmtx+ii)+jj)=(*(*(gcmtx+ii)+jj))-(*(*(kcmtx+ii)+j))*(*(*(gcmtx+j)+jj))/(*(*(kcmtx+j)+j));
				*/
				if(*(*(gcmtx+j)+j)==0.0)
				{
				  *(*(gtmtx+ii)+jj)=(*(*(gcmtx+ii)+jj));    //tsutsumi
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
      currentpivot(j+1,6*nnode);
	  }
	} /*END CONDENSATION*/
    laptime("MATRIX ELIMINATION COMPLETED.",t0);

/*
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
*/

/*縮約された剛性マトリクスの出力*/
#if 0
fprintf(fout,"[ke']\n");
for(ii=0;ii<6*nmultinode;ii++)
{
  sprintf(string,"%3d",ii+1);
  for(jj=0;jj<6*nmultinode;jj++)
  {
	sprintf(s," %18.8f",*(*(kcmtx+ii)+jj));
	strcat(string,s);
  }
  fprintf(fout,"%s\n",string);
}
fprintf(fout,"[kg']\n");
for(ii=0;ii<6*nmultinode;ii++)
{
  sprintf(string,"%3d",ii+1);
  for(jj=0;jj<6*nmultinode;jj++)
  {
	sprintf(s," %18.8f",*(*(gcmtx+ii)+jj));
	strcat(string,s);
  }
  fprintf(fout,"%s\n",string);
}
#endif

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
/*
    AA=(double **)malloc(6*nmultinode*sizeof(double *));
    BB=(double **)malloc(6*nmultinode*sizeof(double *));
    WW=(double **)malloc(6*nmultinode*sizeof(double *));
    EE=(double *)malloc(6*nmultinode*sizeof(double));
    VV=(double **)malloc(6*nmultinode*sizeof(double *));
    CF=(signed char *)malloc(6*nmultinode*sizeof(signed char));
    for(i=0;i<6*nmultinode;i++)
    {
	  *(AA+i)=(double *)malloc(6*nmultinode*sizeof(double));
	  *(BB+i)=(double *)malloc(6*nmultinode*sizeof(double));
	  *(WW+i)=(double *)malloc(6*nmultinode*sizeof(double));
	  *(VV+i)=(double *)malloc(6*nmultinode*sizeof(double));
    }
*/

    kmtx2=(struct gcomponent *)
          malloc(6*nmultinode*sizeof(struct gcomponent));
    for(i=1;i<=6*nmultinode;i++)           /*GLOBAL MATRIX INITIALIZATION.*/
    {
      g=(kmtx2+(i-1))->down;   /*NEXT OF DIAGONAL.*/

      while(g!=NULL) /*CLEAR ROW.*/
      {
        p=g;
        g=g->down;
        free(p);
      }

      ginit.m=(unsigned short int)i;
      /*ginit.n=(unsigned short int)i;*/

      *(kmtx2+(i-1))=ginit;
    }

    gmtx2=(struct gcomponent *)
          malloc(6*nmultinode*sizeof(struct gcomponent));
    for(i=1;i<=6*nmultinode;i++)           /*GLOBAL MATRIX INITIALIZATION.*/
    {
      g=(gmtx2+(i-1))->down;   /*NEXT OF DIAGONAL.*/

      while(g!=NULL) /*CLEAR ROW.*/
      {
        p=g;
        g=g->down;
        free(p);
      }

      ginit.m=(unsigned short int)i;
      /*ginit.n=(unsigned short int)i;*/

      *(gmtx2+(i-1))=ginit;
    }

    confs2=(struct oconf *)malloc(6*nmultinode*sizeof(struct oconf));
    for(j=0;j<6*nmultinode;j++)
    {
      (confs2+j)->iconf=0.0;
    }

    eigen=(double *)malloc(neig*sizeof(double));       /*EIGEN VALUES*/

    gvct2=(double **)malloc(neig*sizeof(double *));   /*GLOBAL VECTORS*/
    for(i=0;i<neig;i++)
    {
      *(gvct2+i)=(double *)malloc(6*nmultinode*sizeof(double));
	  for(j=0;j<6*nmultinode;j++)
      {
        *(*(gvct2+i)+j)=0.0;
      }
    }

	mm=0;
	for(ii=1;ii<=6*nmultinode;ii++)
	{
//      *(CF+mm)=(confs+ii)->iconf;
      (confs2+ii-1)->iconf=(confs+ii-1)->iconf;
	  for(jj=1;jj<=6*nmultinode;jj++)
	  {
//		*(*(BB+mm-1)+nn-1)=*(*(kcmtx+ii-1)+jj-1);
        gdata=*(*(kcmtx+ii-1)+jj-1);
        gwrite(kmtx2,ii,jj,gdata);
	  }
	}
/*
    for(ii=0;ii<6*nmultinode;ii++)
    {
      for(jj=0;jj<6*nmultinode;jj++)
      {
        gread(kmtx2,ii,jj,&gdata);
        sprintf(str,"gcomp[%ld][%ld]=%12.5f\n",ii,jj,gdata);
        if(ii<10 && gdata) errormessage(str);
      }
    }
    for(ii=0;ii<6*nmultinode;ii++)
    {
        sprintf(str,"gcomp[%ld]=%12.5f\n",ii,(kmtx2+ii)->value);
        errormessage(str);
    }
*/
	for(ii=1;ii<=6*nmultinode-1;ii++)
	{
	  for(jj=1;jj<=6*nmultinode-1;jj++)
	  {
//		*(*(AA+mm-1)+nn-1)=*(*(gcmtx+ii-1)+jj-1);
        gdata=*(*(gcmtx+ii-1)+jj-1);
        gwrite(gmtx2,ii,jj,gdata);
	  }
	}

#if 0
    fprintf(fout,"{CONF} :");
    for(ii=0;ii<6*nmultinode;ii++) fprintf(fout," %3d",(confs2+ii)->iconf);
    fprintf(fout,"\n");
#endif

    laptime("CONDENSED MATRIX EXTRACTED.",t0);

	/*deigqrcf(AA,BB,12,12,12,12,eps,WW,EE,VV,CF);*/
    deigabgeneral(gmtx2,kmtx2,confs2,6*nmultinode,neig,neig,eps,eigen,gvct2);
//    bisecsylvester(gmtx2,kmtx2,confs2,6*nmultinode,neig,neig,biseceps,eigen,gvct2);

/*
fprintf(globalfile,"EIGEN VECTOR [V]\n");
for(ii=0;ii<neig;ii++)
{
  if(!*(CF+ii))
  {
	for(jj=0;jj<6*nmultinode;jj++)
	{
	  if(!*(CF+jj))
      {
        fprintf(globalfile," %18.8f",*(*(VV+ii)+jj));
        *(*(gvct2+ii)+jj)=*(*(VV+ii)+jj);
      }
	}
	fprintf(globalfile,"\n");
  }
}
fprintf(globalfile,"\n");


fprintf(globalfile,"EIGEN VALUE INVERSE {E}\n");
for(ii=0;ii<neig;ii++)
{
  if(!*(CF+ii)) fprintf(globalfile," %18.8f\n",*(EE+ii));
}
fprintf(globalfile,"\n");


	for(ii=0;ii<neig;ii++)
	{
	  if(!*(CF+ii))
	  {
		mm=ii;
		break;
	  }
	}
*/

/*
    fprintf(globalfile,"EIGEN VECTOR [V]\n");
    for(jj=0;jj<6*nmultinode;jj++)
    {
      for(ii=0;ii<neig;ii++)
      {
        fprintf(globalfile," %18.8f\n",*(*(gvct2+ii)+jj));
      }
//      fprintf(globalfile,"\n");
    }
*/

/*
    for(ii=0;ii<neig;ii++)
    {
      for(jj=0;jj<6*nmultinode;jj++)
      {
        fprintf(globalfile," %18.8f\n",*(*(gvct2+ii)+jj));
      }
//      fprintf(globalfile,"\n");
    }
*/

/*
    fprintf(globalfile,"EIGEN VALUE INVERSE {E}\n");
    for(ii=0;ii<neig;ii++)
    {
      fprintf(globalfile," %18.8f\n",*(eigen+ii));
    }
    fprintf(globalfile,"\n");
*/

/*
	for(ii=0;ii<neig;ii++)
	{
	  if(*(eigen+ii))
	  {
		mm=ii;
		break;
	  }
	}
*/
	if(*(eigen+mm)==0.0)
    {
      MessageBox(NULL,"EIGEN VALUE=1/0.0","CONDENSE",MB_OK);
      return 0;
    }

//	fprintf(fout,"EIGEN VALUE=%12.5f",1/(*(eigen+0)));
	/*fprintf(fout," LINE=%d",mm);*/
//	fprintf(fout,"\n");

//	fprintf(feig,"EIGEN VALUE=%12.5f SAFETY=%12.5f\n",1/(*(eigen+0)),(*(eigen+0)));                 //tsutsumi20171218
  }

  laptime("EIGEN COMPLETED.",t0);

  af->nlaps=neig;
  af->eigenval=eigen;
  af->eigenvec=gvct;

  for(i=0;i<neig;i++)
  {
    *(eigen+i)*=factor;
    for(ii=0;ii<nmultinode;ii++)
    {
      for(jj=0;jj<6;jj++)
      {
//        sprintf(str,"gvct2[%ld][%d]=%9.8f\n",i,6*ii+jj,*(*(gvct2+i)+6*ii+jj));
//        errormessage(str);
        *(*(gvct+i)+6*(*(loffs+ii))+jj)=*(*(gvct2+i)+6*ii+jj);
      }
    }
  }

  for(i=0;i<neig;i++)
  {
    sprintf(string,"EIGEN VALUE %ld=%.5E",(i+1),1/(*(eigen+i)));
    fprintf(fout,"%s\n",string);
    fprintf(feig,"%s\n",string);
    errormessage(string);
    outputmode(*(gvct+i),feig,nnode,ninit);
    outputmodeII(*(gvct+i),fout,nnode,ninit,loffs,nmultinode);
//    outputmodeII(*(gvct+i),feig,nnode,ninit,loffs,nmultinode);
  }

  for(i=0;i<nmultiwire;i++)
  {
    for(j=0;j<nelem;j++)
    {
      if(((*(multiwire+i))->code)==(elems+j)->code)
      {
        (af->elems+j)->srate[0]=*(eigen+0);  /*Buckling safety ratio*/
        break;
      }
    }
  }

  for(i=0;i<msize;i++)
  {
    (af->confs+i)->iconf=(initialconfs+i)->iconf;
    (af->confs+i)->value=(initialconfs+i)->value;
  }

  free(kcmtx);  free(gcmtx);  free(ktmtx);  free(gtmtx);
  free(kcmtx2);  free(gcmtx2);

  gfree(kmtx,nnode); /*FREE GLOBAL MATRIX.*/
  gfree(gmtx,nnode); /*FREE GLOBAL MATRIX.*/
  gfree(kmtx2,nmultinode);
  gfree(gmtx2,nmultinode);

  updatemode(af,*(gvct+0)); /*FORMATION UPDATE.*/

/*
  for(ii=0;ii<msize;ii++)
  {
    sprintf(str,"gvct[%ld][%d]=%9.8f\n",i,ii,*(*(gvct+i)+ii));
    errormessage(str);
  }
*/

  errormessage(" ");
  errormessage("COMPLETED.");
//  fprintf(fout,"COMPLETED.\n");

  fclose(fout);
  fclose(feig);
//  fclose(frat);

  memory1=availablephysicalmemory("REMAIN:");
  sprintf(string,"CONSUMPTION:%ld[BYTES]",(memory0-memory1));
  errormessage(string);

  /*Buckling safety ratio*/  //ujioka
//  (wdraw.childs+1)->vparam.vflag.ev.srcancolor=1;
//  (wdraw.childs+1)->vparam.vflag.ev.srcanrate=1;

  return 1;
}/*bclng003*/
#endif

#if 1
int bclng003(struct arclmframe *af,struct owire **multiwire,int nmultiwire)
/*ELASTIC BUCKLING FOR ARCLM FRAME.*/
/*AFTER "ARCLM001".*/
/*BUCKLING CONDENSATION INTO MULTIWIRE*/  /*UJIOKA*/
{
  DWORD memory0,memory1;

  FILE /**fin,*/*fout,*feig/*,*frat*/;                   /*FILE 8 BYTES*/
  /*FILE *felem,*fdisp,*freact;*/
  /*char dir[]=DIRECTORY;*/                        /*DATA DIRECTORY*/
  char s[4096],string[4096];
  int i,j,ii,jj,mm,nn;
  int nnode,nelem,nsect/*,nreact*/;
  long int /*fsize,*/loff,ioff,joff,msize;
  long int *loffs,*moffs,*noffs,k,kk;
  int nmultinode;
  /*long int time;*/
  struct gcomponent ginit={0,0,0.0,NULL};
  struct gcomponent *kmtx,*gmtx;                    /*GLOBAL MATRIX*/
  struct gcomponent *kmtx2,*gmtx2;                  /*GLOBAL MATRIX*/
  struct gcomponent *ge,*g,*p;                      /*GLOBAL MATRIX*/
  double **kcmtx,**gcmtx,**ktmtx,**gtmtx;        /*CONDENSED MATRIX*/
  double **kcmtx2,**gcmtx2;
  double **gvct,**gvct2;                            /*GLOBAL VECTOR*/
  double **drccos,**tmatrix,**estiff,*estress;
  /*double determinant,data;*/
  double gdata,factor;
  clock_t t0/*,t1,t2*/;

  /*struct osect *sects;*/
  struct onode /**nodes,*/*ninit;
  struct owire elem;
  struct owire *elems;
  struct oconf *confs,*confs2;
  struct oconf *initialconfs;

  long int neig;
//  double AA[12][12],BB[12][12],WW[12][12],EE[12],VV[12][12];
//  signed char CF[12];
//  double **AA,**BB,**WW,*EE,**VV;
//  signed char *CF;
  double eps=1.0E-16,*eigen,biseceps;

  /*FOR SORT NODES*/
  int flag,tmp;
  struct owire *pe;
  char str[80];
  nmultinode=0;
  /*FOR SORT NODES*/

  memory0=availablephysicalmemory("INITIAL:");   /*MEMORY AVAILABLE*/

  fout=fgetstofopen("\0","w",ID_OUTPUTFILE);          /*OUTPUT FILE*/
  //fout=fgetstofopen("\0","a",ID_OUTPUTFILE);          /*OUTPUT FILE*/
  globalfile=fout;

  GetDlgItemText((wmenu.childs+2)->hwnd,ID_OUTPUTFILE,s,80);
  /*GetDlgItemText((wmenu.childs+2)->hwnd,ID_SECTIONFILE,s,80);*/

  biseceps=BISECEPS;

  nn=strcspn(s,".");
  strncpy(string,s,nn);
  string[nn]='\0';
  strcat(string,".eig");
  feig=fopen(string,"w");
/*
  sprintf(string,"\0");
  strncpy(string,s,nn);
  string[nn]='\0';
  strcat(string,".rat");
  frat=fopen(string,"w");
*/
  fprintf(fout,"CONDENSED MULTIELEM EIGEN VALUES\n");
//  fprintf(feig,"CONDENSED MULTIELEM EIGEN VALUES\n");

  t0=clock();                                        /*CLOCK BEGIN.*/

  nnode=af->nnode;                                 /*INPUT INITIAL.*/
  nelem=af->nelem;
  nsect=af->nsect;
  /*nreact=af->nreact;*/
  sprintf(string,"NODES=%d ELEMS=%d SECTS=%d",nnode,nelem,nsect);
  errormessage(string);
  fprintf(fout,"%s\n",string);
  fprintf(feig,"%s\n",string);

  msize=6*nnode;                           /*SIZE OF GLOBAL MATRIX.*/
  neig=NEIGEN;
//  neig=12;

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
  kcmtx2=(double **)malloc(msize*sizeof(double *));
  gcmtx2=(double **)malloc(msize*sizeof(double *));
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
	*(kcmtx2+i)=(double *)malloc(msize*sizeof(double));
	*(gcmtx2+i)=(double *)malloc(msize*sizeof(double));
	*(ktmtx+i)=(double *)malloc(msize*sizeof(double));
	*(gtmtx+i)=(double *)malloc(msize*sizeof(double));
  }

  gvct=(double **)malloc(neig*sizeof(double *));   /*GLOBAL VECTORS*/
  gvct2=(double **)malloc(neig*sizeof(double *));   /*GLOBAL VECTORS*/
  eigen=(double *)malloc(neig*sizeof(double));       /*EIGEN VALUES*/
  if(gvct==NULL || gvct2==NULL || eigen==NULL) return 0;
  for(i=0;i<neig;i++)
  {
	*(gvct+i)=(double *)malloc(msize*sizeof(double));
    *(gvct2+i)=(double *)malloc(msize*sizeof(double));
	for(j=0;j<msize;j++)
    {
       *(*(gvct+i)+j)=0.0;
       *(*(gvct2+i)+j)=0.0;
    }
  }

  /*fdisp=af->fdisp;*/                 /*DISPLACEMENT:6 DIRECTIONS.*/
  /*felem=af->felem;*/              /*CODE,12 BOUNDARIES,12 STRESS.*/
  /*freact=af->freact;*/                           /*REACTION FILE.*/

  /*sects=af->sects;*/
  /*nodes=af->nodes;*/
  ninit=af->ninit;
  elems=af->elems;
  confs=af->confs;
//  confs2=af->confs;        /*WRONG*/
  confs2=(struct oconf *)malloc(msize*sizeof(struct oconf));

  initialconfs=(struct oconf *)malloc(msize*sizeof(struct oconf));
  for(i=0;i<msize;i++)
  {
    (initialconfs+i)->iconf=(confs+i)->iconf;
    (initialconfs+i)->value=(confs+i)->value;
  }

  GetAsyncKeyState(VK_LBUTTON);                   /*CLEAR KEY LEFT.*/
  GetAsyncKeyState(VK_RBUTTON);                  /*CLEAR KEY RIGHT.*/

  errormessage("BCLNG003:BUCKLING LINEAR.");
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

  /*SORT NODES IN MULTIWIRE*/
  loffs=(long int *)malloc(sizeof(long int));/*loff of multiwire*/

  for(i=0;i<nnode;i++)
  {
    flag=0;
    for(j=0;j<nmultiwire;j++)
    {
      pe=*(multiwire+j);

      for(k=0;k<2;k++)
      {
        if((af->nodes+i)->code==((pe->node[k])->code)) flag=1;
      }
      if(flag) break;
    }

    if(flag)
    {
      if(nmultinode)
      {
        nmultinode+=1;
        loffs=(long int *)realloc(loffs,nmultinode*sizeof(long int));
        *(loffs+nmultinode-1)=(af->nodes+i)->loff;
      }
      else     //initial
      {
        *(loffs+0)=(af->nodes+i)->loff;
        nmultinode=1;
      }
    }
  }
/*
  for(i=0;i<nmultinode;i++)
  {
   sprintf(string,"loffs[%ld]=%ld\n",i,*(loffs+i));
   errormessage(string);
  }
*/

  moffs=(long int *)malloc(nnode*sizeof(long int));
  noffs=(long int *)malloc(nnode*sizeof(long int));

  /*CONDENSATION*/
//  fprintf(fout,"CONDENSATION BEGIN\n");
  sprintf(string,"\nBCLNG003:CONDENSATION BEGIN\n");
  errormessage(string);

  sprintf(string,"CONDENSED NODES=%ld ELEMS=%ld",nmultinode,nmultiwire);
  errormessage(string);
  fprintf(fout,"%s\n",string);
//  fprintf(feig,"%s\n",string);

  for(i=0;i<nmultinode;i++)
  {
    sprintf(string,"NODE  %ld",(af->nodes+(*(loffs+i)))->code);
//    errormessage(string);
//    fprintf(fout,"%s\n",string);
//    fprintf(feig,"%s\n",string);
  }

  for(i=0;i<nmultiwire;i++)
  {
    sprintf(string,"ELEM  %ld",(*(multiwire+i))->code);
//    errormessage(string);
    fprintf(fout,"%s\n",string);
//    fprintf(feig,"%s\n",string);
  }
  fprintf(fout,"\n");
//  fprintf(feig,"\n");

    factor=1.0;
#if 1       /*FOR TENSION COLUMN BUCKLING LOAD*/
	inputelem(elems,af->melem,(*(multiwire+0))->loff,&elem);
    if(elem.stress[0][0]<0)   factor=-1.0;        /*TENSION*/
#endif

  /* for(i=0;i<nelem;i++) */
  {
	for(ii=0;ii<msize;ii++) /*DUPLICATE [Ke],[Kg]*/
	{
	  for(jj=0;jj<msize;jj++)
	  {
		gread(kmtx,ii+1,jj+1,&gdata);
		*(*(kcmtx+ii)+jj)=gdata;
		if(ii!=jj) *(*(kcmtx+jj)+ii)=gdata;
		gread(gmtx,ii+1,jj+1,&gdata);
		*(*(gcmtx+ii)+jj)=factor*gdata;
		if(ii!=jj) *(*(gcmtx+jj)+ii)=factor*gdata;
	  }
	}


/*
    for(i=0;i<6*nnode;i++)
    {
      sprintf(string,"befor:kcmtx[%ld][%ld]=%ld\n",i+1,i+1,*(*(kcmtx+i)+i));
      errormessage(string);
    }
*/
/*
    for(i=0;i<6;i++)
    {

      sprintf(string,"before:kcmtx[%ld[%ld]=%ld\n",6**(loffs+0)+i+1,6**(loffs+0)+i+1,
      				*(*(kcmtx+6**(loffs+0)+i)+6**(loffs+0)+i));
      errormessage(string);
      fprintf(fout,"%s\n",string);
    }
*/

    /*EXCHANGING LINES FOR CONDENSATION*/
    i=0;
    ii=0;
    for(i=0;i<nnode;i++)
    {
      *(moffs+i)=i;
      *(noffs+i)=i;
    }

    for(i=0;i<nnode;i++)
    {
      for(ii=0;ii<nmultinode;ii++)
      {
        if(i==*(loffs+ii))
        {
        tmp=*(noffs+i);
        *(noffs+i)=*(noffs+ii);
        *(noffs+ii)=tmp;
        }
      }
    }

/*
    for(i=0;i<nnode;i++)
    {
      sprintf(string,"moffs[%ld]=%ld\n",i,*(moffs+i));
      errormessage(string);
      fprintf(fout,"%s\n",string);
    }
    for(i=0;i<nnode;i++)
    {
      sprintf(string,"noffs[%ld]=%ld\n",i,*(noffs+i));
      errormessage(string);
      fprintf(fout,"%s\n",string);
    }
*/
    for(i=0;i<msize;i++)
    {
      (confs2+i)->iconf=0;
      (confs2+i)->value=0.0;
    }

    exchangelinesII(kcmtx,gcmtx,confs,kcmtx2,gcmtx2,confs2,
                    moffs,noffs,af->nnode,0);
/*
    for(i=0;i<6*nnode;i++)
    {
      sprintf(string,"iconf[%ld]=%ld\n",i+1,(confs+i)->iconf);
      errormessage(string);
      fprintf(fout,"%s\n",string);
    }
*/
	for(ii=0;ii<msize;ii++)
	{
	  for(jj=0;jj<msize;jj++)
	  {
		*(*(kcmtx+ii)+jj)=*(*(kcmtx2+ii)+jj);
		if(ii!=jj) *(*(kcmtx+jj)+ii)=*(*(kcmtx2+ii)+jj);
		*(*(gcmtx+ii)+jj)=*(*(gcmtx2+ii)+jj);
		if(ii!=jj) *(*(gcmtx+jj)+ii)=*(*(gcmtx2+ii)+jj);
	  }
      *(confs+ii)=*(confs2+ii);
	}

    laptime("LINES EXCHANGED.",t0);

/*
    for(i=0;i<nnode*6;i++)
    {
      sprintf(string,"kcmtx[%ld][%ld]=%ld\n",i+1,i+1,*(*(kcmtx+i)+i));
      errormessage(string);
      fprintf(fout,"%s\n",string);
    }

    for(i=0;i<6*nnode;i++)
    {
      sprintf(string,"iconf[%ld]=%ld\n",i+1,(confs+i)->iconf);
      errormessage(string);
    }

    for(i=0;i<6*nnode;i++)
    {
      sprintf(string,"iconf[%ld]=%ld\n",i+1,(confs+i)->iconf);
      errormessage(string);
      fprintf(fout,"%s\n",string);
    }
*/

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

//sprintf(string,"ELEM ORDER=%d",i+1);
//MessageBox(NULL,string,"CONDENSE",MB_OK);               //tsutsumi20171218

    laptime("MATRIX ELIMINATION BEGIN.",t0);

	/*CONDENSATION*/ /*WITH CONSIDERING CONF*/
	for(j=/*0*/6*nmultinode;j<msize;j++)/*j:縮約で削除する行、列*/
	{
	  if(j>=msize) break;
/*
      sprintf(string,"kcmtx[%ld][%ld]=%ld\n",j+1,j+1,
      				*(*(kcmtx+j)+j));
      errormessage(string);
      fprintf(fout,"%s\n",string);
*/
	  if(!(confs+j)->iconf)
	  {
		if(*(*(kcmtx+j)+j)==0.0)
		{
		  sprintf(string,"INSTABLE TERMINATION K%d%d=%9.3f",j+1,j+1,*(*(kcmtx+j)+j));
		  MessageBox(NULL,string,"CONDENSE",MB_OK);
          fclose(fout);
		  return 0;
		}
/*
sprintf(string,"i=%d j=%d",i,j);
MessageBox(NULL,string,"CONDENSE",MB_OK);
*/
        ii=0;

		while(ii<msize)
		{
          if(ii==6*(nmultinode)) ii=j;
/*
sprintf(string,"i=%d j=%d ii=%d",i,j,ii);
MessageBox(NULL,string,"CONDENSE",MB_OK);
*/
		  if(!(confs+ii)->iconf)
		  {
            jj=0;
			while(jj<msize)
			{
              if(jj==6*(nmultinode)) jj=j;
/*
sprintf(string,"i=%d j=%d ii=%d jj=%d",i,j,ii,jj);
MessageBox(NULL,string,"CONDENSE",MB_OK);
*/
			  if(!(confs+jj)->iconf)
			  /*if(jj!=j)*/
			  {
              /*縮約（弾性剛性マトリクス）*/
				*(*(ktmtx+ii)+jj)=(*(*(kcmtx+ii)+jj))-(*(*(kcmtx+ii)+j))*(*(*(kcmtx+j)+jj))/(*(*(kcmtx+j)+j));

              /*縮約（幾何剛性マトリクス）*/
				/*
				*(*(gtmtx+ii)+jj)=(*(*(gcmtx+ii)+jj))-(*(*(kcmtx+ii)+j))*(*(*(gcmtx+j)+jj))/(*(*(kcmtx+j)+j));
				*/
				if(*(*(gcmtx+j)+j)==0.0)
				{
				  *(*(gtmtx+ii)+jj)=(*(*(gcmtx+ii)+jj));    //tsutsumi
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
      currentpivot(j+1,6*nnode);
	  }
	} /*END CONDENSATION*/
    laptime("MATRIX ELIMINATION COMPLETED.",t0);

/*
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
*/

/*縮約された剛性マトリクスの出力*/
#if 0
fprintf(fout,"[ke']\n");
for(ii=0;ii<6*nmultinode;ii++)
{
  sprintf(string,"%3d",ii+1);
  for(jj=0;jj<6*nmultinode;jj++)
  {
	sprintf(s," %18.8f",*(*(kcmtx+ii)+jj));
	strcat(string,s);
  }
  fprintf(fout,"%s\n",string);
}
fprintf(fout,"[kg']\n");
for(ii=0;ii<6*nmultinode;ii++)
{
  sprintf(string,"%3d",ii+1);
  for(jj=0;jj<6*nmultinode;jj++)
  {
	sprintf(s," %18.8f",*(*(gcmtx+ii)+jj));
	strcat(string,s);
  }
  fprintf(fout,"%s\n",string);
}
#endif

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
/*
    AA=(double **)malloc(6*nmultinode*sizeof(double *));
    BB=(double **)malloc(6*nmultinode*sizeof(double *));
    WW=(double **)malloc(6*nmultinode*sizeof(double *));
    EE=(double *)malloc(6*nmultinode*sizeof(double));
    VV=(double **)malloc(6*nmultinode*sizeof(double *));
    CF=(signed char *)malloc(6*nmultinode*sizeof(signed char));
    for(i=0;i<6*nmultinode;i++)
    {
	  *(AA+i)=(double *)malloc(6*nmultinode*sizeof(double));
	  *(BB+i)=(double *)malloc(6*nmultinode*sizeof(double));
	  *(WW+i)=(double *)malloc(6*nmultinode*sizeof(double));
	  *(VV+i)=(double *)malloc(6*nmultinode*sizeof(double));
    }
*/

    kmtx2=(struct gcomponent *)
          malloc(6*nmultinode*sizeof(struct gcomponent));
    for(i=1;i<=6*nmultinode;i++)           /*GLOBAL MATRIX INITIALIZATION.*/
    {
//      g=(kmtx2+(i-1))->down;   /*NEXT OF DIAGONAL.*/

//      while(g!=NULL) /*CLEAR ROW.*/
/*
      {
        p=g;
        g=g->down;
        free(p);
      }
*/
      ginit.m=(unsigned short int)i;
      /*ginit.n=(unsigned short int)i;*/

      *(kmtx2+(i-1))=ginit;
    }

    gmtx2=(struct gcomponent *)
          malloc(6*nmultinode*sizeof(struct gcomponent));
    for(i=1;i<=6*nmultinode;i++)           /*GLOBAL MATRIX INITIALIZATION.*/
    {
//      g=(gmtx2+(i-1))->down;   /*NEXT OF DIAGONAL.*/

//      while(g!=NULL) /*CLEAR ROW.*/
/*
      {
        p=g;
        g=g->down;
        free(p);
      }
*/
      ginit.m=(unsigned short int)i;
      /*ginit.n=(unsigned short int)i;*/

      *(gmtx2+(i-1))=ginit;
    }

    confs2=(struct oconf *)malloc(6*nmultinode*sizeof(struct oconf));
    for(j=0;j<6*nmultinode;j++)
    {
      (confs2+j)->iconf=0.0;
    }

    eigen=(double *)malloc(neig*sizeof(double));       /*EIGEN VALUES*/

    gvct2=(double **)malloc(neig*sizeof(double *));   /*GLOBAL VECTORS*/
    for(i=0;i<neig;i++)
    {
      *(gvct2+i)=(double *)malloc(6*nmultinode*sizeof(double));
	  for(j=0;j<6*nmultinode;j++)
      {
        *(*(gvct2+i)+j)=0.0;
      }
    }

	mm=0;
	for(ii=1;ii<=6*nmultinode;ii++)
	{
//      *(CF+mm)=(confs+ii)->iconf;
      (confs2+ii-1)->iconf=(confs+ii-1)->iconf;
	  for(jj=1;jj<=6*nmultinode;jj++)
	  {
//		*(*(BB+mm-1)+nn-1)=*(*(kcmtx+ii-1)+jj-1);
        gdata=*(*(kcmtx+ii-1)+jj-1);
        gwrite(kmtx2,ii,jj,gdata);
	  }
	}
/*
    for(ii=0;ii<6*nmultinode;ii++)
    {
      for(jj=0;jj<6*nmultinode;jj++)
      {
        gread(kmtx2,ii,jj,&gdata);
        sprintf(str,"gcomp[%ld][%ld]=%12.5f\n",ii,jj,gdata);
        if(ii<10 && gdata) errormessage(str);
      }
    }
    for(ii=0;ii<6*nmultinode;ii++)
    {
        sprintf(str,"gcomp[%ld]=%12.5f\n",ii,(kmtx2+ii)->value);
        errormessage(str);
    }
*/
	for(ii=1;ii<=6*nmultinode-1;ii++)
	{
	  for(jj=1;jj<=6*nmultinode-1;jj++)
	  {
//		*(*(AA+mm-1)+nn-1)=*(*(gcmtx+ii-1)+jj-1);
        gdata=*(*(gcmtx+ii-1)+jj-1);
        gwrite(gmtx2,ii,jj,gdata);
	  }
	}

#if 0
    fprintf(fout,"{CONF} :");
    for(ii=0;ii<6*nmultinode;ii++) fprintf(fout," %3d",(confs2+ii)->iconf);
    fprintf(fout,"\n");
#endif

    laptime("CONDENSED MATRIX EXTRACTED.",t0);

	/*deigqrcf(AA,BB,12,12,12,12,eps,WW,EE,VV,CF);*/
    deigabgeneral(gmtx2,kmtx2,confs2,6*nmultinode,neig,neig,eps,eigen,gvct2);
//    bisecsylvester(gmtx2,kmtx2,confs2,6*nmultinode,neig,neig,biseceps,eigen,gvct2);

/*
fprintf(globalfile,"EIGEN VECTOR [V]\n");
for(ii=0;ii<neig;ii++)
{
  if(!*(CF+ii))
  {
	for(jj=0;jj<6*nmultinode;jj++)
	{
	  if(!*(CF+jj))
      {
        fprintf(globalfile," %18.8f",*(*(VV+ii)+jj));
        *(*(gvct2+ii)+jj)=*(*(VV+ii)+jj);
      }
	}
	fprintf(globalfile,"\n");
  }
}
fprintf(globalfile,"\n");


fprintf(globalfile,"EIGEN VALUE INVERSE {E}\n");
for(ii=0;ii<neig;ii++)
{
  if(!*(CF+ii)) fprintf(globalfile," %18.8f\n",*(EE+ii));
}
fprintf(globalfile,"\n");


	for(ii=0;ii<neig;ii++)
	{
	  if(!*(CF+ii))
	  {
		mm=ii;
		break;
	  }
	}
*/

/*
    fprintf(globalfile,"EIGEN VECTOR [V]\n");
    for(jj=0;jj<6*nmultinode;jj++)
    {
      for(ii=0;ii<neig;ii++)
      {
        fprintf(globalfile," %18.8f\n",*(*(gvct2+ii)+jj));
      }
//      fprintf(globalfile,"\n");
    }
*/

/*
    for(ii=0;ii<neig;ii++)
    {
      for(jj=0;jj<6*nmultinode;jj++)
      {
        fprintf(globalfile," %18.8f\n",*(*(gvct2+ii)+jj));
      }
//      fprintf(globalfile,"\n");
    }
*/

/*
    fprintf(globalfile,"EIGEN VALUE INVERSE {E}\n");
    for(ii=0;ii<neig;ii++)
    {
      fprintf(globalfile," %18.8f\n",*(eigen+ii));
    }
    fprintf(globalfile,"\n");
*/

/*
	for(ii=0;ii<neig;ii++)
	{
	  if(*(eigen+ii))
	  {
		mm=ii;
		break;
	  }
	}
*/
	if(*(eigen+mm)==0.0)
    {
      MessageBox(NULL,"EIGEN VALUE=1/0.0","CONDENSE",MB_OK);
      return 0;
    }

//	fprintf(fout,"EIGEN VALUE=%12.5f",1/(*(eigen+0)));
	/*fprintf(fout," LINE=%d",mm);*/
//	fprintf(fout,"\n");

//	fprintf(feig,"EIGEN VALUE=%12.5f SAFETY=%12.5f\n",1/(*(eigen+0)),(*(eigen+0)));                 //tsutsumi20171218
  }

  laptime("EIGEN COMPLETED.",t0);

  af->nlaps=neig;
  af->eigenval=eigen;
  af->eigenvec=gvct;

  for(i=0;i<neig;i++)
  {
    *(eigen+i)*=factor;
    for(ii=0;ii<nmultinode;ii++)
    {
      for(jj=0;jj<6;jj++)
      {
//        sprintf(str,"gvct2[%ld][%d]=%9.8f\n",i,6*ii+jj,*(*(gvct2+i)+6*ii+jj));
//        errormessage(str);
        *(*(gvct+i)+6*(*(loffs+ii))+jj)=*(*(gvct2+i)+6*ii+jj);
      }
    }
  }

  for(i=0;i<neig;i++)
  {
    sprintf(string,"EIGEN VALUE %ld=%.5E",(i+1),1/(*(eigen+i)));
    fprintf(fout,"%s\n",string);
    fprintf(feig,"%s\n",string);
    errormessage(string);
    outputmode(*(gvct+i),feig,nnode,ninit);
    outputmodeII(*(gvct+i),fout,nnode,ninit,loffs,nmultinode);
//    outputmodeII(*(gvct+i),feig,nnode,ninit,loffs,nmultinode);
  }

  for(i=0;i<nmultiwire;i++)
  {
    for(j=0;j<nelem;j++)
    {
      if(((*(multiwire+i))->code)==(elems+j)->code)
      {
        (af->elems+j)->srate[0]=*(eigen+0);  /*Buckling safety ratio*/
        break;
      }
    }
  }

  for(i=0;i<msize;i++)
  {
    (af->confs+i)->iconf=(initialconfs+i)->iconf;
    (af->confs+i)->value=(initialconfs+i)->value;
  }

  free(kcmtx);  free(gcmtx);  free(ktmtx);  free(gtmtx);
  free(kcmtx2);  free(gcmtx2);

  gfree(kmtx,nnode); /*FREE GLOBAL MATRIX.*/
  gfree(gmtx,nnode); /*FREE GLOBAL MATRIX.*/
  gfree(kmtx2,nmultinode);
  gfree(gmtx2,nmultinode);

  updatemode(af,*(gvct+0)); /*FORMATION UPDATE.*/

/*
  for(ii=0;ii<msize;ii++)
  {
    sprintf(str,"gvct[%ld][%d]=%9.8f\n",i,ii,*(*(gvct+i)+ii));
    errormessage(str);
  }
*/

  errormessage(" ");
  errormessage("COMPLETED.");
//  fprintf(fout,"COMPLETED.\n");

  fclose(fout);
  fclose(feig);
//  fclose(frat);

  memory1=availablephysicalmemory("REMAIN:");
  sprintf(string,"CONSUMPTION:%ld[BYTES]",(memory0-memory1));
  errormessage(string);

  /*Buckling safety ratio*/  //ujioka
//  (wdraw.childs+1)->vparam.vflag.ev.srcancolor=1;
//  (wdraw.childs+1)->vparam.vflag.ev.srcanrate=1;

  return 1;
}/*bclng003*/
#endif
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

/*-----------------------------------------------------------------*/

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
  /*for(nn=0;nn<100;nn+column10_1+)*/
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


if(globalmessageflag)
{
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
}

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


if(globalmessageflag)
{
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
}

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

if(globalmessageflag)
{
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
if(globalmessageflag)
{fprintf(globalfile,"CHECK:[L][Lt]=[B]\n");
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
}

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
	  if(globalmessageflag) fprintf(globalfile,"NOTE : SKIPPED LINE=%d DIAGONAL=0.0\n",ii+1);
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


if(globalmessageflag)
{
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
}

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
if(globalmessageflag)
{
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
}

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
/*=================================================================*/
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

/*-----------------------------------------------------------------*/
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
  long int i,ii,kk;
  long int J,K,nev,neva,nvec,R;
  double T;

  unsigned int mm,mmm;
  double lkj,lim;
  struct gcomponent *gi,*gj,*gr,*gk,*gp,*gg;
  struct gcomponent *bj,*bm;

  long int ipiv; /*PIVOT LINE.*/
  signed char ic; /*CONFINEMENT ID.*/
  long int ndim; /*MATRIX DIMENSIONS.*/

  long int iii,jjj;
  double data;

  long int jj;

  errormessage("DEIGAB:BEGIN.");

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
/*
  sprintf(str,"LINE %ld T=%f",ipiv+1,T);
  errormessage(str);
*/
  if(T<=0.0)
  {
    errormessage("DEIGAB 2:MATRIX [B] IS NOT POSITIVE DEFINITE.");
      sprintf(str,"LINE %ld",ipiv+1);
      errormessage(str);
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
                        gg=gdefine((unsigned int)kk,
                                   (unsigned int)ii,
                                   -lkj*lim*(gj->value),NULL,NULL);
                        gi->down=gg;
                        gp=gi;
                        gi=gg;
                      }
                      else if(gi->m > kk) /*INSERT.*/
                      {
                        gg=gdefine((unsigned int)kk,
                                   (unsigned int)ii,
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
                        gg=gdefine((unsigned int)ii,
                                   (unsigned int)kk,
                                   -lkj*lim*(gj->value),NULL,NULL);
                        gk->down=gg;
                        gp=gk;
                        gk=gg;
                      }
                      else if(gk->m > ii) /*INSERT.*/
                      {
                        gg=gdefine((unsigned int)ii,
                                   (unsigned int)kk,
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
for(iii=25;iii<=N;iii++)
{
  for(jjj=25;jjj<=N;jjj++)
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
for(iii=25;iii<=N;iii++)
{
  for(jjj=25;jjj<=N;jjj++)
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
  errormessage("DEIGAB:END.");
  /*currentvalue("DEIGABGENERAL END.",N,NE,NULL,NULL,E,V);*/

  return;
}/*deigabgeneral*/

/*-----------------------------------------------------------------*/
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

  errormessage("DEIGRS:BEGIN.");

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
            gg=gdefine((unsigned int)(mpiv+1),
                       (unsigned int)(K+1),
                       (A1+SR),
                       gj->down,NULL);
            gj->down=gg;
          }
          else if((gj->m)>(mpiv+1)) /*FILL.*/
          {
            gg=gdefine((unsigned int)(mpiv+1),
                       (unsigned int)(K+1),
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
                      gg=gdefine((unsigned int)ii,
                                 (unsigned int)jj,
                                 -S,NULL,NULL);
                      gj->down=gg;
                      gp=gj;
                      gj=gg;
                    }
                    else if(gj->m > ii) /*FILL IN.*/
                    {
                      gg=gdefine((unsigned int)ii,
                                 (unsigned int)jj,
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
    F=R; /*F:UPPER BOUND FOR DECENDING,LOWER FOR ACENDING.*/
    for(I=0;I<NEA;I++) *(E+I)=-R;

    for(K=0;K<NEA;K++)
    {
      currentpivot(K,NEA);

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

  free(lines);
  free(LW);
  free(gpp);
  for(i=0;i<6;i++) free(*(W+i));
  free(W);

  errormessage("DEIGRS:END.");
  return;
}/*deigrsstandard*/

/*-----------------------------------------------------------------*/
struct gcomponent *gdefine(unsigned int m,
                           unsigned int n,
                           double value,
						   struct gcomponent *down,
                           struct gcomponent *left)
/*DEFINE GCOMP PARAMETERS.*/
{
  struct gcomponent *g;

  g=(struct gcomponent *)malloc(sizeof(struct gcomponent));
  if(g==NULL) {
    errormessage("gdefine: memory error");
    return NULL;
  }

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
  long int k;
  struct gcomponent *gcpy,*go,*gc,*gi;
  struct gcomponent ginit={0,0,0.0,NULL};

  gcpy=(struct gcomponent *)malloc(msize*sizeof(struct gcomponent));
  if (gcpy==NULL) {
	  errormessage("copygconpmatrix memory error");
      return NULL;
  }
  for(k=0;k<msize;k++)
  {
    (gcpy+k)->down=NULL;
    ginit.m=(unsigned int)(k+1);
    *(gcpy+k)=ginit;
  }

  for(k=0;k<msize;k++)
  {
    (gcpy+k)->m=(gmtx+k)->m; /*DIAGONAL.*/
    (gcpy+k)->n=(gmtx+k)->n;
    (gcpy+k)->value=(gmtx+k)->value;

    go=(gmtx+k);
    gc=(gcpy+k);
    while(go->down!=NULL) /*COPY COMPS IN COLUMN.*/
    {
      go=go->down;

      gi=gdefine((unsigned int)go->m,
                 (unsigned int)go->n,
                 go->value,NULL,NULL);
      gc->down=gi;

      gc=gc->down;
    }
    gc->down=NULL;
  }

  return gcpy;
}/*copygcompmatrix*/

struct gcomponent *copygcompmatrix2(struct gcomponent *gmtx,
                                    long int msize,double factor)
/*CREATE COPY OF GCOMP MATRIX.*/
{
  long int k;
  struct gcomponent *gcpy,*go,*gc,*gi;
  struct gcomponent ginit={0,0,0.0,NULL};

  gcpy=(struct gcomponent *)malloc(msize*sizeof(struct gcomponent));
  if (gcpy==NULL) {
      errormessage("copygconpmatrix memory error");
      return NULL;
  }
  for(k=0;k<msize;k++)
  {
    (gcpy+k)->down=NULL;
    ginit.m=(unsigned int)(k+1);
    *(gcpy+k)=ginit;
  }

  for(k=0;k<msize;k++)
  {
    (gcpy+k)->m=(gmtx+k)->m; /*DIAGONAL.*/
    (gcpy+k)->n=(gmtx+k)->n;
    (gcpy+k)->value=factor*((gmtx+k)->value);

	go=(gmtx+k);
    gc=(gcpy+k);
    while(go->down!=NULL) /*COPY COMPS IN COLUMN.*/
    {
      go=go->down;

      gi=gdefine((unsigned int)go->m,
                 (unsigned int)go->n,
                 factor*(go->value),NULL,NULL);
      gc->down=gi;

      gc=gc->down;
    }
    gc->down=NULL;
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

void outputmodeII(double *gvct,FILE *fout,int nnode,struct onode *nodes,
                long int *loffs,int nmultinode)
/*OUTPUT MODE DISPLACEMENT.*/
{
  char string[256];
  int i,j;
  double data[6];

  for(i=0;i<nmultinode;i++)
  {
    for(j=0;j<6;j++) data[j]=*(gvct+6*(*(loffs+i))+j);
    fprintf(fout,"NODE:%5ld {dU}=",(nodes+(*(loffs+i)))->code);
    sprintf(string," %12.5E %12.5E %12.5E %12.5E %12.5E %12.5E",
            data[0],data[1],data[2],data[3],data[4],data[5]);
    if(fout!=NULL) fprintf(fout,"%s\n",string);
  }
}/*outputmodeII*/

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

double bclngoutputtomemoryII(FILE *ftext,struct arclmframe *af)
/*TRANSLATE BCLNG OUTPUTFILE TEXT INTO MEMORY.*/
/*FOR HIGHER BCLNG MODES. 2021.02.15 UJIOKA*/
/*
  RETURN:BUCKLING EIGEN VALUE
  2021.08.25 UJIOKA

  RETURN:THE INNER PRODUCT OF LOAD VECTOR AND BUCKLING MODE VECTOR.
  2021.03.01 UJIOKA
*/
{
  char **data,str[256]="\0";
  char *s;
  int i,ii,j,jj,m,n;
  long int ncode;
  double ddata;
  int laps;
  double dsafety;

  double eigen;               /*Buckling Eigen Value*/
  double inproduct;           /*inner product*/
  double loadnorm,modenorm;   /*norm of {P},{Ue}*/

  getincrement((wmenu.childs+2)->hwnd,&laps,&dsafety);
  sprintf(str,"LAP:%d",laps);
//  if(MessageBox(NULL,str,"OPEN RESULT",MB_OKCANCEL)==IDCANCEL)  return;

  fseek(ftext,0L,SEEK_SET);

  af->ddisp=(double *)malloc(6*(af->nnode)*sizeof(double));
									   /*DISPLACEMENT:6 DIRECTIONS.*/

//  while(strncmp(str,"EIGEN",5)) fgets(str,256,ftext);
  m=0;
  while(1)
  {
    data=fgetsbrk(ftext,&n);
  //    fgets(str,256,ftext);
//    if(strncmp(str,"LAP:",4)==0 && MessageBox(NULL,str,"OPEN RESULT",MB_OKCANCEL)==IDOK) break;
    if(strncmp(*data,"EIGEN",5)==0) m++;
    if(m==laps) break;
    free(data);
  }
//  fgets(str,256,ftext);

  sprintf(str,"%s %s %s",*(data),*(data+1),*(data+2));
  if(MessageBox(NULL,str,"OPEN RESULT",MB_OKCANCEL)==IDCANCEL) return -1;
  errormessage(str);
  SetDlgItemText((wmenu.childs+2)->hwnd,ID_LAPS,*(data+2));
  sprintf((wdraw.childs+1)->title,"BUCKLING %s",str);

  s = strtok(*(data+2),"=");
  s = strtok(NULL,"=");
  eigen = strtod(s,NULL);

  free(data);

  for(i=0;i<(af->nnode);i++) /*DISPLACEMENTS.*/
  {
	data=fgetsbrk(ftext,&n);
	if(n!=9) return -1;

	ncode=strtol(*(data+1),NULL,10);
	if(ncode!=(af->nodes+i)->code) return -1;

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

  inproduct=0;
  loadnorm=0;
  modenorm=0;

  for(ii=0;ii<(af->nnode);ii++)
  {
    for(jj=0;jj<6;jj++)
    {
    inproduct+=((af->confs+6*ii+jj)->value)*(*(af->ddisp+6*ii+jj)-((af->ninit+ii)->d[jj]));
    /*<{P},{Ue}>(inner product)*/

    loadnorm+=((af->confs+6*ii+jj)->value)*((af->confs+6*ii+jj)->value);
    /*|{P}|^2*/

    modenorm+=(*(af->ddisp+6*ii+jj)-((af->ninit+ii)->d[jj]))*(*(af->ddisp+6*ii+jj)-((af->ninit+ii)->d[jj]));
    /*|{Ue}|^2*/

    /*FOR TEST*/
    /*
    if((af->confs+6*ii+jj)->value!=0)
     {
      sprintf(str,"%.5f",inproduct);
      errormessage(str);
     }
    */
  	}
  }

  loadnorm=sqrt(loadnorm);         /*|{P}|*/
  modenorm=sqrt(modenorm);         /*|{Ue}|*/
  inproduct/=(loadnorm*modenorm);  /*<{P},{Ue}>/(|{P}||{Ue}|)(normalization)*/
  inproduct=abs(inproduct);

  sprintf(str,"%.5f",inproduct);
  SetDlgItemText((wmenu.childs+2)->hwnd,ID_SAFETY,str);

  sprintf(str,"COSINE: <{P},{Ue}>/(|{P}||{Ue}|) = %.5f",inproduct);
  errormessage(str);

//  return inproduct;
  return eigen;
}/*bclngoutputtomemoryII*/

double bclngoutputtomemoryIII(FILE *ftext,struct arclmframe *af,struct owire ***multiwire,int *nmultiwire)
/*TRANSLATE BCLNG OUTPUTFILE TEXT INTO MEMORY.*/
/*FOR BCLNG003. 2021.08.25 UJIOKA*/
/*
  RETURN:BUCKLING EIGEN VALUE
*/
{
  char **data,str[256]="\0";
  char *s;
  int i,ii,j,jj,m,n;
  long int ncode;
  int code;
  int nmultinode;
  double ddata;
  int laps;
  double dsafety;

  double eigen;               /*Buckling Eigen Value*/
  double inproduct;           /*inner product*/
  double loadnorm,modenorm;   /*norm of {P},{Ue}*/

  getincrement((wmenu.childs+2)->hwnd,&laps,&dsafety);
  sprintf(str,"LAP:%d",laps);
//  if(MessageBox(NULL,str,"OPEN RESULT",MB_OKCANCEL)==IDCANCEL)  return;

  fseek(ftext,0L,SEEK_SET);

  af->ddisp=(double *)malloc(6*(af->nnode)*sizeof(double));
									   /*DISPLACEMENT:6 DIRECTIONS.*/

  while(strncmp(str,"NODES",4)) fgets(str,256,ftext);
  data=fgetsbrk(ftext,&n);
  sprintf(str,"%s %s %s",*data,*(data+1),*(data+2));
  errormessage(str);

  s = strtok(*(data+1),"=");
  s = strtok(NULL,"=");

  nmultinode = strtol(s,NULL,10);

  s = strtok(*(data+2),"=");
  s = strtok(NULL,"=");

  *nmultiwire = strtol(s,NULL,10);

  *multiwire=(struct owire **)malloc((*nmultiwire)*sizeof(struct owire *));
  free(data);

  for(i=0;i<*nmultiwire;i++)
  {
    data=fgetsbrk(ftext,&n);
    code=strtol(*(data+1),NULL,10);
    for(j=0;j<af->nelem;j++)
    {
      if(code==(af->elems+j)->code)
      {
        *(*multiwire+i)=af->elems+j;
        break;
      }
    }
    free(data);
  }

/*
  for(i=0;i<*nmultiwire;i++)
  {
      sprintf(str,"MULTIWIRE[%ld]:CODE=%ld",i+1,(*(*multiwire+i))->code);
      errormessage(str);
  }
*/

  m=0;
  while(1)
  {
    data=fgetsbrk(ftext,&n);
//    fgets(str,256,ftext);
//    if(strncmp(str,"LAP:",4)==0 && MessageBox(NULL,str,"OPEN RESULT",MB_OKCANCEL)==IDOK) break;
    if(strncmp(*data,"EIGEN",5)==0) m++;
    if(m==laps) break;
    free(data);
  }
//  fgets(str,256,ftext);

  sprintf(str,"%s %s %s",*(data),*(data+1),*(data+2));
  if(MessageBox(NULL,str,"OPEN RESULT",MB_OKCANCEL)==IDCANCEL) return 0;
  errormessage(str);
  SetDlgItemText((wmenu.childs+2)->hwnd,ID_LAPS,*(data+2));
  sprintf((wdraw.childs+1)->title,"BUCKLING CONDENSATION: %s",str);

  s = strtok(*(data+2),"=");
  s = strtok(NULL,"=");
  eigen = strtod(s,NULL);
/*
  sprintf(str,"%f",eigen);
  errormessage(str);
*/
  free(data);

  for(i=0;i<nmultinode;i++) /*DISPLACEMENTS.*/
  {
	data=fgetsbrk(ftext,&n);
	if(n!=9) return 0;

	ncode=strtol(*(data+1),NULL,10);
//	if(ncode!=(af->nodes+i)->code) return;
	for(jj=0;jj<af->nnode;jj++)
    {
      if(ncode==(af->nodes+jj)->code)
      {
        ii=(af->nodes+jj)->loff;
        break;
      }
    }
	for(j=0;j<3;j++)
	{
	  ddata=strtod(*(data+j+3),NULL);
	  ddata+=(af->ninit+ii)->d[j];
	  (af->nodes+ii)->d[j]=ddata;
	  *(af->ddisp+6*ii+j)=ddata;
    }

    for(j=3;j<6;j++)
    {
      ddata=strtod(*(data+j+3),NULL);
      *(af->ddisp+6*ii+j)=ddata;
	}

    for(;n>0;n--) free(*(data+n-1));
    free(data);
  }

  inproduct=0;
  loadnorm=0;
  modenorm=0;

  for(ii=0;ii<(af->nnode);ii++)
  {
    for(jj=0;jj<6;jj++)
    {
    inproduct+=((af->confs+6*ii+jj)->value)*(*(af->ddisp+6*ii+jj)-((af->ninit+ii)->d[jj]));
    /*<{P},{Ue}>(inner product)*/

    loadnorm+=((af->confs+6*ii+jj)->value)*((af->confs+6*ii+jj)->value);
    /*|{P}|^2*/

    modenorm+=(*(af->ddisp+6*ii+jj)-((af->ninit+ii)->d[jj]))*(*(af->ddisp+6*ii+jj)-((af->ninit+ii)->d[jj]));
    /*|{Ue}|^2*/

    /*FOR TEST*/
    /*
    if((af->confs+6*ii+jj)->value!=0)
     {
      sprintf(str,"%.5f",inproduct);
      errormessage(str);
     }
    */
  	}
  }

  loadnorm=sqrt(loadnorm);         /*|{P}|*/
  modenorm=sqrt(modenorm);         /*|{Ue}|*/
  inproduct/=(loadnorm*modenorm);  /*<{P},{Ue}>/(|{P}||{Ue}|)(normalization)*/
  inproduct=abs(inproduct);

  sprintf(str,"%.5f",inproduct);
  SetDlgItemText((wmenu.childs+2)->hwnd,ID_SAFETY,str);

  sprintf(str,"COSINE: <{P},{Ue}>/(|{P}||{Ue}|) = %.5f",inproduct);
  errormessage(str);

  return eigen;
}/*bclngoutputtomemoryIII*/

void freqencyoutputtomemory(FILE *ftext,double **freq)
/*TRANSLATE BCLNG OUTPUTFILE TEXT INTO MEMORY.*/
{
  char **data,str[256]="\0",non[256];
  int i,j,n;
  long int ncode;
  double ddata;

  /*fseek(ftext,0L,SEEK_SET);*/

  *(freq)=(double *)malloc(200*sizeof(double));

  while(strncmp(str,"ACCELERATION=",13)) fgets(str,256,ftext);

  for(i=0;i<200;i++) /*FREQENCY.*/
  {
	data=fgetsbrk(ftext,&n);
	if(n!=3) return;

	ddata=strtod(*(data+2),NULL);

	*(*(freq)+i)=sqrt(ddata)/2.0/PI;

//sprintf(non,"ddata=%10f freq=%10f",ddata,*(*(freq)+i));
//MessageBox(NULL,non, "Eigen Mode", MB_OK);

	for(;n>0;n--) free(*(data+n-1));
	free(data);
  }

  return;
}/*freqencyoutputtomemory*/

void eigenmodeoutputtomemory(FILE *ftext,struct arclmframe *af,double *dfact)
/*TRANSLATE BCLNG OUTPUTFILE TEXT INTO MEMORY.*/
{
  char **data,str[256]="\0",non[256];
  int i,j,n;
  long int ncode;
  double ddata;
  double dx,dy,dz,ddist,ddistmax;

  /*fseek(ftext,0L,SEEK_SET);*/

  *dfact=0.0,ddistmax=0.0;

  af->ddisp=(double *)malloc(6*(af->nnode)*sizeof(double));
									   /*DISPLACEMENT:6 DIRECTIONS.*/

  while(strncmp(str,"EIGEN",5)) fgets(str,256,ftext);

  for(i=0;i<(af->nnode);i++) /*DISPLACEMENTS.*/
  {
    dx=0.0,dy=0.0,dz=0.0,ddist=0.0;

	data=fgetsbrk(ftext,&n);
	if(n!=9) return;

	ncode=strtol(*(data+1),NULL,10);
	if(ncode!=(af->nodes+i)->code) return;

	for(j=0;j<3;j++)
	{
	  ddata=strtod(*(data+j+3),NULL);

	  if(j==0)      dx=ddata;
	  else if(j==1) dy=ddata;
	  else if(j==2) dz=ddata;

	  ddata+=(af->ninit+i)->d[j];
	  (af->nodes+i)->d[j]=ddata;
	  *(af->ddisp+6*i+j)=ddata;
	}
	for(j=3;j<6;j++)
	{
	  ddata=strtod(*(data+j+3),NULL);
	  *(af->ddisp+6*i+j)=ddata;
	}

	ddist=sqrt(dx*dx+dy*dy+dz*dz);
	if(ddist>ddistmax)
	{
	  ddistmax=ddist;
	  *dfact=1.0/ddistmax;
	}
//sprintf(non,"dx=%10f dy=%10f dz=%10f ddist=%10f dfact=%3f",dx,dy,dz,ddist,dfact);
//MessageBox(NULL,non, "Eigen Mode", MB_OK);

	for(;n>0;n--) free(*(data+n-1));
	free(data);
  }

  return;
}/*eigenmodeoutputtomemory*/


void bisecsylvester(struct gcomponent *A,
                    struct gcomponent *B,
                    struct oconf *confs,
                    long int N,long int NE,long int NV,
					double EPS,
					double *E,double **V)
{
  char err[256],str[256];
  int i,j,k;
  double lambda,lastlambda,sum;
  double determinant,sign;
  struct gcomponent *gmtx,*gcomp1;
  double *lastvec,*evec;
  int neg;
  int nnode=N/6;
  double LL,LR;
  LL=0.0;
  LR=BISECRIGHT;

  lastvec=(double *)malloc(N*sizeof(double));
  evec=(double *)malloc(N*sizeof(double));
  if(lastvec==NULL || evec==NULL) return;

  for(i=0;i<N;i++)
  {
	*(lastvec+i)=0.0;
  }

  for(i=0;i<NE;i++)
  {
	neg=0;
    sprintf(err,"EIG=%d",i+1);
	errormessage(err);
    lambda=0.5*(LL+LR);
    lastlambda=lambda;

	/* BISECTION METHOD */
    while(1)
	{
bisecstart:
	  /* CALCULATE (A-λB)x=b */
      gmtx=gcomponentadd2(A,B,-1.0*lambda,N);
	  croutlu(gmtx,confs,N,&determinant,&sign,gcomp1);

      for(j=0;j<N;j++)
      {
        /* CHECK IF QUADRATIC FORM xAx<0.0 */
        if((confs+j)->iconf==0 && (gmtx+j)->value>0.0)
        {
          neg++;
          if(neg>i)
          {
            sprintf(err,"LR-LL=%.3E",LR-LL);
            errormessage(err);
            if(LR-LL<EPS)
            {
              gfree(gmtx,nnode);
              goto bisecend;
            }
			sprintf(err,"LAMBDA<%.14f",1.0/lambda);
			errormessage(err);
			LL=lambda;
            lambda=0.5*(LL+LR);
			gfree(gmtx,nnode);
			goto bisecstart;
          }
        }
      }
	  sprintf(err,"LAMBDA>%.14f",1.0/lambda);
      errormessage(err);

      /* CALCULATE CANDIDATE FOR EIGENVECTOR */
      for(j=0;j<N;j++)
      {
        *(evec+j)=((double)rand()+1.0)/((double)RAND_MAX+2.0);
      }
      forwardbackward(gmtx,evec,confs,N,gcomp1);

      /* NORMALIZE */
      sum=0.0;
      for(j=0;j<N;j++)
      {
        if((confs+j)->iconf==0)
        {
          sum+=evec[j]*evec[j];
        }
      }
      sum=sqrt(sum);
      for(j=0;j<N;j++)
	  {
        if((confs+j)->iconf==0)
        {
          if(sum>0.0)
          {
            evec[j]/=sum;
          }
          lastvec[j]=evec[j];
        }
        else
        {
          lastvec[j]=0.0;
        }
	  }
	  FILE *ftmp=fopen("mode.tmp","w");
	  if(ftmp!=NULL)
      {
		char string[256];
        int ii;
        fprintf(ftmp,"DEIGABGENERAL EIGEN VALUE %ld=%.5f\n",(i+1),1.0/lastlambda);
		if(lastlambda>0.0)
        {
		  double Ti=2.0*PI*sqrt(lastlambda);
		  sprintf(string,"PERIOD T%ld=%.5f [sec]",(i+1),Ti);
          fprintf(ftmp,"%s\n",string);
        }

        fprintf(ftmp,"\nDEIGABGENERAL EIGEN VECTOR %ld\n",(i+1));
        for(ii=0;ii<nnode;ii++)
        {
          fprintf(ftmp,"%4ld %11.7f %11.7f %11.7f %11.7f %11.7f %11.7f\n",
                       ii,
                       lastvec[6*ii+0],
                       lastvec[6*ii+1],
                       lastvec[6*ii+2],
                       lastvec[6*ii+3],
					   lastvec[6*ii+4],
                       lastvec[6*ii+5]);
		}
		fclose(ftmp);
	  }

      lastlambda=lambda;
	  LR=lambda;
      lambda=0.5*(LL+LR);
      gfree(gmtx,nnode);
	}
bisecend:

    /* SET EIGENVALUE & EIGENVECTOR */
    E[i]=lastlambda;
    for(j=0;j<N;j++)
    {
      V[i][j]=lastvec[j];
    }

	/* SHIFT FOR NEXT EIGENVALUE */
	LR=LL;
    LL=0.0;
  }

  free(lastvec);
  return;
}/*bisecsylvester*/



void bisecgeneral(struct gcomponent *A,double factorA,
				  struct gcomponent *B,double factorB,
				  struct oconf *confs,
				  long int N,long int NE,double defsign,
				  double EPS,
				  double *E,double **V,
				  double BL, double BR)
{
  char err[256],str[256];
  int i,j,k,ii;
  double lambda;
  double determinant,sign;
  struct gcomponent *gmtx,*gcomp1;
  double *evct;
  double eigen;
  double neg;
  int nnode=N/6;
  MSG msg;
  double LL=0.0;
  double LR=BISECRIGHT;
  double LM;

  if(BL!=NULL)LL=BL;
  if(BR!=NULL)LR=BR;
  if(defsign==NULL)defsign=0;


  /*EIGEN VALUE BOUND CHECKING.*/
  while(1)
  {
	while (PeekMessage(&msg, NULL, 0, 0, PM_REMOVE))
	{
	  TranslateMessage(&msg);
	  DispatchMessage(&msg);
	}
	neg=0;
	//gmtx=gcomponentadd2(A,B,-1.0*LR,N);

	gmtx=gcomponentadd3(A,factorA,B,factorB*LR,N);

	croutlu(gmtx,confs,N,&determinant,&sign,gcomp1);
	for(j=0;j<N;j++)
	{
	  if((confs+j)->iconf==0 && (gmtx+j)->value>0.0)
	  {
		neg++;
		if(neg>defsign)
		{
		  LR*=2.0;
		  gfree(gmtx,nnode);
		  break;
		}
	  }
	}
	if(neg==defsign)
	{
	  gfree(gmtx,nnode);
	  break;
	}
  }
  sprintf(err,"LL=%f LR=%f\n",LL,LR);
  errormessage(err);

  evct=(double *)malloc(N*sizeof(double));

  for(i=0;i<NE;i++)
  {
	sprintf(err,"NEIG=%ld",i);
	errormessage(err);
	neg=0;

	LM=0.5*(LL+LR);
	lambda=LR;/*SAFETY SIDE*/

	/* BISECTION METHOD */
	while(1)
	{
	  //MESSAGE FOR UPDATE UI
	  //AVOID FREEZE ON LONG RUNNING TASK
	  bisecstart:

	  while (PeekMessage(&msg, NULL, 0, 0, PM_REMOVE))
	  {
		TranslateMessage(&msg);
		DispatchMessage(&msg);
	  }

      while(GetAsyncKeyState(VK_LBUTTON))  /*LEFT CLICK TO CONTINUE.*/
	  {
		if(GetAsyncKeyState(VK_RBUTTON))      /*RIGHT CLICK TO ABORT.*/
		{
		  return;
		}
	  }

	  //gmtx=gcomponentadd2(A,B,-1.0*LM,N);
	  gmtx=gcomponentadd3(A,factorA,B,factorB*LM,N);

	  croutlu(gmtx,confs,N,&determinant,&sign,gcomp1);

	  if(LR-LL<EPS)break;

	  neg=0;
	  for(j=0;j<N;j++)
	  {
		/* CHECK IF QUADRATIC FORM xAx<0.0 */
		if((confs+j)->iconf==0 && (gmtx+j)->value>0.0)
		{
		  neg++;
		  if(neg>defsign+i)
		  {
			sprintf(err,"EIGENVALUE<%.14f NEG=%f",LM, neg);
			errormessage(err);
			LL=LM;
			LM=0.5*(LL+LR);
			gfree(gmtx,nnode);
			goto bisecstart;
		  }
		}
	  }
	  sprintf(err,"EIGENVALUE>%.14f NEG=%f",LM, neg);
	  errormessage(err);
	  lambda=LM;
	  LR=LM;
	  LM=0.5*(LL+LR);
	  gfree(gmtx,nnode);
	}

	/* INITIAL EIGEN VECTOR */
	for(j=0;j<N;j++)
	{
	  if((confs+j)->iconf==0)
	  {
		*(evct+j)=((double)rand()+1.0)/((double)RAND_MAX+2.0);
	  }
	  else
	  {
		*(evct+j)=0.0;
	  }
	}

	/*INVERS METHOD FOR EIGEN VECTOR*/
	eigen = inversemethod(gmtx,confs,evct,N);

	/* SET EIGENVALUE & EIGENVECTOR */
	E[i]=LM/*lambda*/;
	for(j=0;j<N;j++)
	{
	  V[i][j]=evct[j];
	}
	/* SHIFT FOR NEXT EIGENVALUE */
	LL=0.0;
	if(BL!=NULL)LL=BL;
  }
  free(evct);
  return;
}/*bisecgeneral*/

double inversemethod(struct gcomponent *gmtx, struct oconf *confs, double *evct, int msize)
{
	char string[500];
	int i,nline;
	int iteration = 0;
	double len = 1.0;

	double* lastevct;
	double eigen,evctevct,evctlastevct;
	struct gcomponent *gcomp1;

	lastevct=(double *)malloc(msize*sizeof(double));

	vectornormalize(evct, msize);

	while (len > 1.0e-12 && iteration<20)/*INVERSE METHOD*/
	{
		for (i = 0; i < msize; i++)
		{
			*(lastevct + i) = *(evct + i);
		}
		nline = forwardbackward(gmtx, evct, confs, msize, gcomp1);

		/*QUADRATIC FORM*/
		/*[K]{x}=λ{x} <=> λ={x}^t[K]{x}/{x}^t{x}*/
		evctlastevct = dotproduct(evct, lastevct, msize);/*{x}^t[K]{x} evct:{x] lastevct:[K]{x}*/
		evctevct = dotproduct(evct, evct, msize);
		eigen = evctlastevct / evctevct;

		vectornormalize(evct, msize);

		for (i = 0; i < msize; i++)
		{
			*(lastevct + i) = abs(*(lastevct + i)) - abs(*(evct + i));
		}
		len = vectorlength(lastevct, msize);/*FOR CHECKING CONVERGENCE.*/

		//sprintf(string, "INVERSE ITERATION : STANDARD EIGENVALUE= %e LEN= %e\n", eigen, len);
		//errormessage(string);
		iteration++;
	}
	free(lastevct);

	return eigen;
}


struct gcomponent *gcomponentadd2(struct gcomponent *mtx1,
								  struct gcomponent *mtx2,
								  double factor,
								  int msize)
{
  long int i,j;
  struct gcomponent *gcomp,*gdown1,*gdown2,*pdown1;
  struct gcomponent *rtn;

  rtn=copygcompmatrix(mtx1,msize);

  if(factor==0.0)
  {
	return rtn;
  }

  for(j=1;j<=msize;j++)
  {
	gdown1=(rtn+(j-1));
	gdown2=(mtx2+(j-1));
	gdown1->value+=factor*gdown2->value;
	while(gdown2->down!=NULL) /*DOWNWARD.*/
	{
	  gdown2=gdown2->down;

	  i=gdown2->m;
	  while((gdown1->m)<i && (gdown1->down)!=NULL)/*(gdown1->m) >= (gdown2->m) or (gdown1->down) == NULL*/
	  {
		pdown1=gdown1;
		gdown1=gdown1->down;
	  }

	  if(gdown1->m ==i)/*(gdown1->m) == (gdown2->m)*/
	  {
		gdown1->value+=factor*gdown2->value;/*UPDATE NEW GCOMPONENT.*/
	  }
	  else if(gdown1->m < i)/*(gdown1->down) == NULL*/
      {
		gcomp=gdefine(i,j,factor*(gdown2->value),
					  NULL,NULL);/*ADD NEW GCOMPONENT.*/
		gdown1->down=gcomp;
		pdown1=gdown1;
		gdown1=gcomp;
	  }
	  else if(gdown1->m > i)/*(gdown1->m) > (gdown2->m)*/
	  {
		gcomp=gdefine(i,j,factor*(gdown2->value),
					  gdown1,NULL);/*INSERT NEW GCOMPONENT.*/
		pdown1->down=gcomp;
		gdown1=gcomp;
	  }
	}
  }
  return rtn;
}/*gcomponentadd2*/



struct gcomponent *gcomponentadd3(struct gcomponent *mtx1,double factor1,
								  struct gcomponent *mtx2,double factor2,
								  int msize)
{
  long int i,j;
  struct gcomponent *gcomp,*gdown1,*gdown2,*pdown1;
  struct gcomponent *rtn;

  rtn=copygcompmatrix2(mtx1,msize,factor1);

  if(factor2==0.0)
  {
	return rtn;
  }

  for(j=1;j<=msize;j++)
  {
	gdown1=(rtn+(j-1));
	gdown2=(mtx2+(j-1));
	gdown1->value+=factor2*gdown2->value;
	while(gdown2->down!=NULL) /*DOWNWARD.*/
	{
	  gdown2=gdown2->down;

	  i=gdown2->m;
	  while((gdown1->m)<i && (gdown1->down)!=NULL)/*(gdown1->m) >= (gdown2->m) or (gdown1->down) == NULL*/
	  {
		pdown1=gdown1;
		gdown1=gdown1->down;
	  }

	  if(gdown1->m ==i)/*(gdown1->m) == (gdown2->m)*/
	  {
		gdown1->value+=factor2*(gdown2->value);/*UPDATE NEW GCOMPONENT.*/
	  }
	  else if(gdown1->m < i)/*(gdown1->down) == NULL*/
	  {
		gcomp=gdefine(i,j,factor2*(gdown2->value),
					  NULL,NULL);/*ADD NEW GCOMPONENT.*/
		gdown1->down=gcomp;
        pdown1=gdown1;
		gdown1=gcomp;
	  }
	  else if(gdown1->m > i)/*(gdown1->m) > (gdown2->m)*/
	  {
		gcomp=gdefine(i,j,factor2*(gdown2->value),
					  gdown1,NULL);/*INSERT NEW GCOMPONENT.*/
		pdown1->down=gcomp;
		gdown1=gcomp;
	  }
	}
  }
  return rtn;
}/*gcomponentadd3*/


void definencr(struct arclmframe *af,double *ncr)       /*UJIOKA*/
{
  FILE *fout,*fload;                                    /*FILE 8 BYTES*/
  FILE *frat = NULL;
  char fname[1024],str[100];
  int n,i;

  char s[80],string[1024];
  int ii,nn;

  double *nz,*bsafety;                                 /*Nz,buckling safety*/
  double nu=1.00;
//  double nu=2.17/1.5;

  nz=mallocdoublevector(af->nelem);
  bsafety=mallocdoublevector(af->nelem);

/*OPEN FILES*/
  fout=fgetstofopen("\0","r",ID_OUTPUTFILEZ);          /*OTL FILE*/

  n=strcspn((wdraw.childs+1)->otpfile,".");            /*RATE FILE*/
  strncpy(fname,(wdraw.childs+1)->otpfile,n);
  fname[n]='\0';
//  strcat(fname,".rat");
  strcat(fname,".brat");
  frat=fopen(fname,"r");

  fload=fopen("bclngloads.txt","w");

  if(fout==NULL)
  {
	  errormessage("No OTL File");
	  MessageBox(NULL,"No OTL File","DefineNcr",MB_OK);
	  return;
  }

  if(frat==NULL)
  {
	  errormessage("No Rate File");
	  MessageBox(NULL,"No Rate File","DefineNcr",MB_OK);
	  MessageBox(NULL,fname,"DefineNcr",MB_OK);
	  return;
  }

/*READ FILES*/
  frameoutputtomemory(fout,&*af);
  readsrcanrate(frat,&*af);

  for(i=0;i<af->nelem;i++)
  {
	  nz[i]=(af->elems+i)->stress[0][0];              /*Nz*/
	  bsafety[i]=(af->elems+i)->srate[0];             /*buckling safety*/
//      if(nz[i]>0)                                     /*compression*/
	  if(nz[i]*bsafety[i]>0.0)                            /*Ncr:compression*/  //20210816.Ujioka
	  {
		  ncr[i]=nz[i]/bsafety[i]/nu;                    /*Ncr*/
/*
		  sprintf(str,"ELEM %d :Nz=%5.8f,1/λ'=%12.5f,Ncr=%5.8f\n",
					(af->elems+i)->code,nz[i],bsafety[i],ncr[i]);
		  errormessage(str);
*/
		  sprintf(str,"ELEM %d SECT %d :Nz= %5.8f 1/λ'=%12.5f Ncr'= %5.8f\n",
					(af->elems+i)->code,(af->elems+i)->sect->code,
					nz[i],bsafety[i],ncr[i]);
//     	  errormessage(str);
		  fprintf(fload,str);
	  }
	  else ncr[i]=0.0;                                  /*tension*/
  }

/*CLOSE FILES*/
  fclose(fout);
  fclose(frat);
  fclose(fload);

return;                                            /*Return Ncr*/
}
/*definencr*/

int bclng011(struct arclmframe *af,struct arclmframe *af0)
/*ELASTIC BUCKLING FOR ARCLM FRAME WITH PRE-STRESS.

  ([KE]-[KG0]){Ue}=λ[KG]{Ue}
  [KE]:ELASTIC
  [KG0]:GEOMETRIC FROM PRE-STRESS
  [KG]:GEOMETRIC FROM LOAD
*/
{
  DWORD memory0,memory1;

  FILE *fin,*fout;                               /*FILE 8 BYTES*/
  /*FILE *felem,*fdisp,*freact;*/
  /*char dir[]=DIRECTORY;*/                        /*DATA DIRECTORY*/
  char /*s[80],*/string[400];
  int i,j,ii,jj;
  int nnode,nelem,nsect/*,nreact*/;
  long int /*fsize,*/loff,msize;
  /*long int time;*/
  struct gcomponent ginit={0,0,0.0,NULL};
  struct gcomponent *kmtx,*gmtx;          /*GLOBAL MATRIX:[KE],[KG]*/
  struct gcomponent *gmtx0;                   /*GLOBAL MATRIX:[KG0]*/
  double **gvct;                                    /*GLOBAL VECTOR*/
  double **drccos,**tmatrix,**estiff,*estress;
  double *estress0;                                    /*PRE-STRESS*/
  /*double determinant,data;*/
  clock_t t0/*,t1,t2*/;

  /*struct osect *sects;*/
  struct onode /**nodes,*/*ninit;
  struct owire elem;
  struct owire *elems;
  struct oconf *confs;
  struct memoryelem *melem;

  long int neig;
  double eps=1.0E-16,*eigen,biseceps;

  memory0=availablephysicalmemory("INITIAL:");   /*MEMORY AVAILABLE*/

  fout=fgetstofopen("\0","w",ID_OUTPUTFILE);          /*OUTPUT FILE*/

  t0=clock();                                        /*CLOCK BEGIN.*/

  nnode=af->nnode;                                 /*INPUT INITIAL.*/
  nelem=af->nelem;
  nsect=af->nsect;
  /*nreact=af->nreact;*/
  sprintf(string,"NODES=%d ELEMS=%d SECTS=%d",nnode,nelem,nsect);
  errormessage(string);
  fprintf(fout,"%s\n",string);

  msize=6*nnode;                           /*SIZE OF GLOBAL MATRIX.*/
  neig=NEIGEN;                             /*NUMBERS OF EIGEN VALUES TO BE OBTAINED.*/

  //if(SOLVER==1)
  //{
  //  eps=BISECEPS;
    /* neig=1; */
  //}

  biseceps=BISECEPS;   //ujok

  kmtx=(struct gcomponent *)          /*DIAGONALS OF GLOBAL MATRIX.*/
        malloc(msize*sizeof(struct gcomponent));
  gmtx=(struct gcomponent *)
        malloc(msize*sizeof(struct gcomponent));
  gmtx0=(struct gcomponent *)
        malloc(msize*sizeof(struct gcomponent));
  if(kmtx==NULL || gmtx==NULL || gmtx0==NULL) return 0;
  for(i=0;i<msize;i++)
  {
    (kmtx+i)->down=NULL;              /*GLOBAL MATRIX INITIALIZATION.*/
    (gmtx+i)->down=NULL;
    (gmtx0+i)->down=NULL;
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
  ninit=af->ninit;                                /*struct onode*/
  elems=af->elems;
  confs=af->confs;

  GetAsyncKeyState(VK_LBUTTON);                  /*CLEAR KEY LEFT.*/
  GetAsyncKeyState(VK_RBUTTON);                  /*CLEAR KEY RIGHT.*/

  errormessage("BCLNG001:BUCKLING LINEAR.");
  availablephysicalmemory("REMAIN:");            /*MEMORY AVAILABLE*/

  for(i=1;i<=msize;i++)             /*GLOBAL MATRIX INITIALIZATION.*/
  {
	ginit.m=(unsigned int)i;
	/*ginit.n=(unsigned int)i;*/
	*(kmtx+(i-1))=ginit;
    *(gmtx+(i-1))=ginit;
    *(gmtx0+(i-1))=ginit;
  }
  /*comps=msize;*/ /*INITIAL COMPONENTS=DIAGONALS.*/

  laptime("ASSEMBLING GLOBAL MATRIX.",t0);

  estress=(double *)malloc(12*sizeof(double));
  estress0=(double *)malloc(12*sizeof(double));

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

    /*-----------------------------------------------------------------*/
    estiff=assememtx(elem);            /*ELASTIC MATRIX OF ELEMENT.*/
    estiff=modifyhinge(elem,estiff);               /*MODIFY MATRIX.*/
    estiff=transformation(estiff,tmatrix);         /*[K]=[Tt][k][T]*/

    assemgstiffness(kmtx,estiff,&elem);       /*ASSEMBLAGE ELASTIC.*/

    /*-----------------------------------------------------------------*/
    /*[KG0]*/
    for(ii=0;ii<=11;ii++) free(*(estiff+ii));
    free(estiff);

    for(ii=0;ii<=1;ii++)                                /*STRESSES.*/
    {
      for(jj=0;jj<6;jj++)
      {
        *(estress0+6*ii+jj)=(af0->elems+i-1)->stress[ii][jj];
        /*
        if(ii==0 && jj==0)
        {
          sprintf(string,"ELEM %ld Nz=%.5E",((af->elems+i-1)->code),*(estress0+6*ii+jj));
          errormessage(string);           //for check
        }
        */
      }
    }

    estiff=assemgmtx(elem,estress0);
    estiff=modifyhinge(elem,estiff);               /*MODIFY MATRIX.*/
    estiff=transformation(estiff,tmatrix);         /*[K]=[Tt][k][T]*/
    for(ii=0;ii<12;ii++)
    {
      for(jj=0;jj<12;jj++) *(*(estiff+ii)+jj)=-*(*(estiff+ii)+jj);
    }

    assemgstiffness(gmtx0,estiff,&elem);     /*ASSEMBLAGE GEOMETRIC.*/

    /*-----------------------------------------------------------------*/
    for(ii=0;ii<=11;ii++) free(*(estiff+ii));
    free(estiff);

    for(ii=0;ii<=1;ii++)                                /*STRESSES.*/
    {
      for(jj=0;jj<6;jj++)
      {
        *(estress+6*ii+jj)=elem.stress[ii][jj];
        /*
        if(ii==0 && jj==0)
        {
          sprintf(string,"ELEM %ld Nz=%.5E",((af->elems+i-1)->code),*(estress+6*ii+jj));
          errormessage(string);           //for check
        }
        */
      }

      }
    estiff=assemgmtx(elem,estress);
    estiff=modifyhinge(elem,estiff);               /*MODIFY MATRIX.*/
    estiff=transformation(estiff,tmatrix);         /*[K]=[Tt][k][T]*/
    for(ii=0;ii<12;ii++)
    {
      for(jj=0;jj<12;jj++) *(*(estiff+ii)+jj)=-*(*(estiff+ii)+jj);
    }

    assemgstiffness(gmtx,estiff,&elem);     /*ASSEMBLAGE GEOMETRIC.*/

    /*-----------------------------------------------------------------*/
    for(ii=0;ii<=11;ii++) free(*(estiff+ii));
	free(estiff);

    for(ii=0;ii<=2;ii++) free(*(drccos+ii));
    free(drccos);
    for(ii=0;ii<=11;ii++) free(*(tmatrix+ii));
    free(tmatrix);
  }
  free(estress);
  laptime("GLOBAL MATRIX ASSEMBLED.",t0);

  kmtx=gcomponentadd2(kmtx,gmtx0,-1.0,msize);    /*[K]=[KE]-[KG0]*/
  //NEED TO CHECK

  /*currentvalue("GLOBAL MATRIX:[K]",msize,neig,kmtx,NULL,NULL,NULL);*/
  /*currentvalue("GLOBAL MATRIX:[G]",msize,neig,gmtx,NULL,NULL,NULL);*/

  //SOLVER  ujioka
  if(globalmessageflag==0||MessageBox(NULL,"DEIGABGENERAL","SOLVER",MB_OKCANCEL)==IDOK)
		deigabgeneral(gmtx,kmtx,confs,msize,neig,neig,eps,eigen,gvct);
  else if(MessageBox(NULL,"BISECSYLVESTER","SOLVER",MB_OKCANCEL)==IDOK)
		bisecsylvester(gmtx,kmtx,confs,msize,neig,neig,biseceps,eigen,gvct);

  /*
  if(SOLVER==0)
  {
    deigabgeneral(gmtx,kmtx,confs,msize,neig,neig,eps,eigen,gvct);
  }
  else
  {
	bisecsylvester(gmtx,kmtx,confs,msize,neig,neig,biseceps,eigen,gvct);
  }
  */

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

  /*STRESS UPDATE WITH BUCKLING STRESS.*/ //20210917.
  if(MessageBox(NULL,"Stress Update with Buckling Stress.","BCLNG011",MB_OKCANCEL)==IDOK)
  {
    for(i=0;i<af->nelem;i++)
    {
      for(ii=0;ii<=1;ii++)
      {
        for(jj=0;jj<6;jj++)
        {
          *(estress0+6*ii+jj)=(af0->elems+i)->stress[ii][jj];

          (af->melem+((af->elems+i)->loff))->stress[ii][jj]*=*(eigen+0);
          (af->melem+((af->elems+i)->loff))->stress[ii][jj]+=*(estress0+6*ii+jj);
/*
          if(ii==0 && jj==0)
          {
            sprintf(string,"ELEM %ld Nz=%.5f",(af->elems+i)->code,(af->elems+i)->stress[ii][jj]);
            errormessage(string);
          }
*/
        }
      }
    }

   stressintofile(af); //extract

   for(i=0;i<=7;i++)   //view setting
   {
     if(i==1||i==4)    //column and brace
     {
     (wdraw.childs+1)->vparam.vflag.ev.etype[i]=1;

     (wdraw.childs+1)->vparam.vflag.ev.stress[i][0]=1;
     }
   }
  }

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
}/*bclng011*/

int bclng101(struct arclmframe *af,int idinput)
/*ELASTO-PLASTIC BUCKLING ANALYSIS FOR ARCLM FRAME.*/
/*INCREMENTAL ANALYSIS:CALCULATING EIGEN-VALUES AND MODES AT EACH STEP(LAP)*/
/*NEWTON-RAPTHON METHOD IS NOT USED...*/
{
  DWORD memory0,memory1,memory2;

  FILE *fin,*fout,*ffig,*ferr,*fsrf;                       /*FILE 8 BYTES*/
  FILE *feig;
  double *ddisp,*dreact;
  struct memoryelem *melem;
  /*FILE *felem,*fdisp,*freact;*/
  char dir[]=DIRECTORY;                            /*DATA DIRECTORY*/
  char s[80],string[400]/*,fname[256]*/;
  int i,j,ii,jj,nn;
  int nnode,nelem,nsect,nlap,laps,nreact;
  long int loffset,msize;
  long int time;
  struct gcomponent ginit={0,0,0.0,NULL};
  struct gcomponent *kmtx,*gmtx;                    /*GLOBAL MATRIX*/
  struct gcomponent *g,*p;                          /*GLOBAL MATRIX*/
  double *gvct;                                     /*GLOBAL VECTOR*/
  double *gvctcmq;         /*****CMQ ZOBUN*****/    /*GLOBAL VECTOR*/
  double **eigenvct,*eigen;
  double **drccos,**tmatrix,**estiff,*estress;
  double determinant,sign,safety,dsafety;
  double func[2];
  clock_t t0,t1,t2;

  struct osect *sects;
  struct onode *nodes,*ninit;
  struct owire elem;
  struct owire *elems;
  struct oconf *confs;

  double initarea;
  long int neig;
  double eps=1.0E-16,biseceps;

  memory0=availablephysicalmemory("INITIAL:");   /*MEMORY AVAILABLE*/

  fin=fgetstofopen(dir,"r",idinput);              /*OPEN INPUT FILE*/
  if(fin==NULL)
  {
    errormessage("ACCESS IMPOSSIBLE.");
    return 0;
  }
  fout=fgetstofopen("\0","w",ID_OUTPUTFILE);          /*OUTPUT FILE*/
  ffig=fopen("hogan.fig","w");                    /*HYSTERISIS FILE*/
  ferr=fopen("hogan.txt","w");                 /*ERROR MESSAGE FILE*/
  fsrf=fopen("surface.txt","w");                     /*SURFACE FILE*/

  //  feig=fopen("eigen.txt","w");                         /*EIGEN FILE*/
  GetDlgItemText((wmenu.childs+2)->hwnd,ID_OUTPUTFILE,s,80);
  /*GetDlgItemText((wmenu.childs+2)->hwnd,ID_SECTIONFILE,s,80);*/
  nn=strcspn(s,".");
  strncpy(string,s,nn);
  string[nn]='\0';
  strcat(string,".eig");
  feig=fopen(string,"w");


  /*fgets(string,256,fin);*/                    /*INPUT APPELATION.*/
  /*errormessage(string);*/

  t0=clock();                                        /*CLOCK BEGIN.*/

  getincrement((wmenu.childs+2)->hwnd,&laps,&dsafety);

  inputinit(fin,&nnode,&nelem,&nsect);             /*INPUT INITIAL.*/
  sprintf(string,"NODES=%d ELEMS=%d SECTS=%d",nnode,nelem,nsect);
  errormessage(string);
  if(fout!=NULL) fprintf(fout,"%s\n",string);

  msize=6*nnode;                           /*SIZE OF GLOBAL MATRIX.*/

  biseceps=BISECEPS;
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

  eigenvct=(double **)malloc(neig*sizeof(double *));   /*GLOBAL VECTORS*/
  eigen=(double *)malloc(neig*sizeof(double));       /*EIGEN VALUES*/
  if(eigenvct==NULL || eigen==NULL) return 0;
  for(i=0;i<neig;i++)
  {
    *(eigenvct+i)=(double *)malloc(msize*sizeof(double));
    for(j=0;j<msize;j++) *(*(eigenvct+i)+j)=0.0;
  }

  gvct=(double *)malloc(msize*sizeof(double));      /*GLOBAL VECTOR*/

  /*****CMQ ZOBUN*****/
  gvctcmq=(double *)malloc(msize*sizeof(double));   /*GLOBAL VECTOR*/

  if(gvct==NULL) return 0;
  for(i=0;i<=msize-1;i++)
  {
    *(gvct+i)=0.0;                  /*GLOBAL VECTOR INITIALIZATION.*/
    *(gvctcmq+i)=0.0;               /*GLOBAL VECTOR INITIALIZATION.*/
  }

  if(idinput==ID_INPUTFILE)
  {
    free(af->sects);
    free(af->nodes);
    free(af->ninit);
    free(af->elems);
    free(af->confs);
    free(af->ddisp);
    free(af->melem);

    sects=(struct osect *)malloc(nsect*sizeof(struct osect));
    if(sects==NULL) return 0;
    nodes=(struct onode *)malloc(nnode*sizeof(struct onode));
    if(nodes==NULL) return 0;
    ninit=(struct onode *)malloc(nnode*sizeof(struct onode));
    if(ninit==NULL) return 0;
    elems=(struct owire *)malloc(nelem*sizeof(struct owire));
    if(elems==NULL) return 0;
    confs=(struct oconf *)malloc(msize*sizeof(struct oconf));
    if(confs==NULL) return 0;
    ddisp=(double *)malloc(6*nnode*sizeof(double));
    if(ddisp==NULL) return 0;
    melem=(struct memoryelem *)
          malloc(nelem*sizeof(struct memoryelem));
    if(melem==NULL) return 0;

    af->sects=sects;
    af->nodes=nodes;
    af->ninit=ninit;
    af->elems=elems;
    af->confs=confs;
    af->ddisp=ddisp;                   /*DISPLACEMENT:6 DIRECTIONS.*/
    af->melem=melem;

    inputtexttomemory(fin,af);      /*READY TO READ LONG REACTIONS.*/
    nnode=af->nnode;
    nelem=af->nelem;
    nsect=af->nsect;
    nreact=af->nreact;

    initialform(nodes,ddisp,nnode);         /*ASSEMBLAGE FORMATION.*/
    initialelem(elems,melem,nelem);          /*ASSEMBLAGE ELEMENTS.*/

    dreact=(double *)malloc(nreact*sizeof(double));     /*REACTION.*/
    af->dreact=dreact;
    initialreact(fin,dreact,nreact);   /*ASSEMBLAGE LONG REACTIONS.*/
  }
  else
  {
    ddisp=af->ddisp;
    melem=af->melem;
    dreact=af->dreact;

    sects=af->sects;
    nodes=af->nodes;
    ninit=af->ninit;
    elems=af->elems;
    confs=af->confs;

    nnode=af->nnode;
    nelem=af->nelem;
    nsect=af->nsect;
    nreact=af->nreact;

    inputloadtomemory(fin,af); /*INPUT LOAD.*/

    if(nreact!=af->nreact)
    {
      MessageBox(NULL,"Input Error.","Arclm101",MB_OK);
      return 0;
    }
  }

  GetAsyncKeyState(VK_LBUTTON);                   /*CLEAR KEY LEFT.*/
  GetAsyncKeyState(VK_RBUTTON);                  /*CLEAR KEY RIGHT.*/

  drawglobalaxis((wdraw.childs+1)->hdcC,(wdraw.childs+1)->vparam,
                 0,0,255);                      /*DRAW GLOBAL AXIS.*/

  if(wsurf.hwnd!=NULL)
  {
    drawyieldsurface((wsurf.childs+1)->hdcC,
                     (wsurf.childs+1)->vparam,2,4,0,NULL);
    overlayhdc(*(wsurf.childs+1),SRCPAINT);       /*UPDATE DISPLAY.*/
  }
  af->fsurface=fopen("canbin.sfc","wb+");      /*STRESS ON SURFACE.*/

  /*fclose(fout);*/

  /*****CMQ ZOBUN*****/
  for(i=1;i<=nelem;i++)               /*ASSEMBLAGE GLOBAL MATRIX.*/
  {
    elem=*(elems+i-1);                       /*READ ELEMENT DATA.*/

    inputnode(ddisp,elem.node[0]);                         /*HEAD*/
    inputnode(ddisp,elem.node[1]);                         /*TAIL*/

    drccos=directioncosine(elem.node[0]->d[0],
                           elem.node[0]->d[1],
                           elem.node[0]->d[2],
                           elem.node[1]->d[0],
                           elem.node[1]->d[1],
                           elem.node[1]->d[2],
                           elem.cangle);               /*[DRCCOS]*/

    tmatrix=transmatrix(drccos);         /*TRANSFORMATION MATRIX.*/

  	modifycmq(melem,&elem);
  	assemcmq101(elem,tmatrix,confs,gvctcmq,dsafety);    /*ASSEMBLAGE CMQ AS LOADS.*/

    for(ii=0;ii<=2;ii++) free(*(drccos+ii));
    free(drccos);
    for(ii=0;ii<=11;ii++) free(*(tmatrix+ii));
    free(tmatrix);
  }
  /*****CMQ ZOBUN*****/

  for(nlap=1;nlap<=laps;nlap++)
  {
    /*sprintf(fname,"arclm%d.lap",nlap);*/
    /*fout=fopen(fname,"w");*/                          /*LAP FILE*/

    /*af->nlaps=nlap;*/
    af->nlaps=1;

    sprintf(string,"LAP:%d/%d",nlap,laps);
    errormessage(string);
    if(fout!=NULL) fprintf(fout,"%s\n",string);
    if(ffig!=NULL) fprintf(ffig,"%s",string);

    setincrement((wmenu.childs+2)->hwnd,
                 laps,nlap,dsafety,(nlap*dsafety));

    memory1=availablephysicalmemory("REMAIN:");  /*MEMORY AVAILABLE*/
/*ARCLM101 PROCESS:CALCULATE DEFORMATION and STRESS*/
    for(i=1;i<=msize;i++)           /*GLOBAL MATRIX INITIALIZATION.*/
    {
      g=(kmtx+(i-1))->down; /*NEXT OF DIAGONAL.*/
      while(g!=NULL) /*CLEAR ROW.*/
      {
        p=g;
        g=g->down;
        free(p);
      }

      ginit.m=(unsigned short int)i;
      /*ginit.n=(unsigned short int)i;*/

      *(kmtx+(i-1))=ginit;
    }
    comps=msize; /*INITIAL COMPONENTS=DIAGONALS.*/

    laptime("ASSEMBLING GLOBAL MATRIX.",t0);

    for(i=1;i<=nelem;i++)               /*ASSEMBLAGE GLOBAL MATRIX.*/
    {
     inputelem(elems,melem,i-1,&elem);        /*READ ELEMENT DATA.*/

     for(ii=0;ii<=1;ii++)
      {
        for(jj=0;jj<=5;jj++)
        {
          (elems+i-1)->iconf[ii][jj]=elem.iconf[ii][jj];
        }
      }
      inputnode(ddisp,elem.node[0]);                         /*HEAD*/
      inputnode(ddisp,elem.node[1]);                         /*TAIL*/

      elem.sect=(elems+i-1)->sect;             /*READ SECTION DATA.*/

      if((wdraw.childs+1)->hdcC!=NULL)     /*DRAW DEFORMED ELEMENT.*/
      {
        drawglobalwire((wdraw.childs+1)->hdcC,
                       (wdraw.childs+1)->vparam,
                       *af,elem,255,255,255,
                                255,255,255,0,ONSCREEN);
      }

      drccos=directioncosine(elem.node[0]->d[0],
                             elem.node[0]->d[1],
                             elem.node[0]->d[2],
                             elem.node[1]->d[0],
                             elem.node[1]->d[1],
                             elem.node[1]->d[2],
                             elem.cangle);               /*[DRCCOS]*/

      tmatrix=transmatrix(drccos);         /*TRANSFORMATION MATRIX.*/
      estiff=assememtx(elem);          /*ELASTIC MATRIX OF ELEMENT.*/
      estiff=modifyhinge(elem,estiff);             /*MODIFY MATRIX.*/
      estiff=assempmtx(elem,estiff);          /*ADD PLASTIC MATRIX.*/
      estiff=transformation(estiff,tmatrix);       /*[K]=[Tt][k][T]*/

      assemgstiffness(kmtx,estiff,&elem);             /*ASSEMBLAGE.*/

      for(ii=0;ii<=2;ii++) free(*(drccos+ii));
      free(drccos);
      for(ii=0;ii<=11;ii++) free(*(tmatrix+ii));
      free(tmatrix);
      for(ii=0;ii<=11;ii++) free(*(estiff+ii));
      free(estiff);
    }
    sprintf(string,"GLOBAL MATRIX %ld COMPS ASSEMBLED.",comps);
    laptime(string,t0);

    overlayhdc(*(wdraw.childs+1),SRCPAINT);       /*UPDATE DISPLAY.*/

    /*errormessage("ASSEMBLAGE GLOBAL VECTOR.");*/
#if 1
    /*****CMQ ZOBUN*****/
    for(i=0;i<=msize-1;i++)
    {
      *(gvct+i)=*(gvctcmq+i);                      /*GLOBAL VECTOR.*/
    }
  	/*****CMQ ZOBUN*****/
#endif

//    assemconf(confs,gvct,dsafety,nnode);           /*GLOBAL VECTOR.*/
    assemconf201(confs,gvct,dsafety,nnode);        /*GLOBAL VECTOR.*/

    modifygivend(kmtx,gvct,confs,nnode);   /*0:LOAD 1:DISPLACEMENT.*/

    laptime("CROUT LU DECOMPOSITION.",t0);
    croutludecomposition(kmtx,
                         gvct,confs,
                         6*nnode,
                         &determinant,&sign);        /*[K]{dU}={dF}*/

    sprintf(string,"LOG{DETERMINANT}=%.5E SIGN=%.1f COMPS=%ld",
            determinant,sign,comps);
    errormessage(string);
    if(fout!=NULL) fprintf(fout,"%s\n",string);

    safety=nlap*dsafety;
    sprintf(string,"SAFETY FACTOR=%.5f",safety);
    errormessage(string);
    if(fout!=NULL) fprintf(fout,"%s\n",string);
    if(ffig!=NULL) fprintf(ffig," SAFETY= %.5f",safety);

    if(sign<=0.0)
    {
      errormessage(" ");
      errormessage("INSTABLE TERMINATION.");
      if(fout!=NULL) fprintf(fout,"INSTABLE TERMINATION.\n");

      laptime("\0",t0);

      fclose(fin);
      fclose(fout);
      fclose(ffig);
      /*fclose(felem);*/
      /*fclose(fdisp);*/
      /*fclose(freact);*/

      gfree(kmtx,nnode); /*FREE GLOBAL MATRIX.*/
      free(gvct);
      /*free(confs);*/

      memory2=availablephysicalmemory("REMAIN:");
      sprintf(string,"CONSUMPTION:%ld[BYTES]",(memory0-memory2));
      errormessage(string);

      return 1;
    }
/*ARCLM101 PROCESS:CALCULATE DEFORMATION and STRESS*/

//    laptime("OUTPUT INTO FILE.",t0);

if(ferr!=NULL) fprintf(ferr,"LAP %3d / %3d",nlap,laps);
if(feig!=NULL) fprintf(feig,"LAP %3d / %3d\n",nlap,laps);

    if(fout!=NULL) fprintf(fout,"\"DISPLACEMENT\"\n");
    outputdisp(gvct,fout,nnode,nodes);  /*INCREMENTAL DISPLACEMENT.*/
    /*while(!GetAsyncKeyState(VK_LBUTTON))
    ;*/                                   /*LEFT CLICK TO CONTINUE.*/

    if(fout!=NULL) fprintf(fout,"\"STRESS\"\n");
    for(i=1;i<=nelem;i++)                   /*STRESS OUTPUT,UPDATE.*/
    {
      inputelem(elems,melem,i-1,&elem);

      inputnode(ddisp,elem.node[0]);
      inputnode(ddisp,elem.node[1]);

      elem.sect=(elems+i-1)->sect;             /*READ SECTION DATA.*/

      estress=elemstress(&elem,gvct,melem,fout,func);   /*STRESS UPDATE.*/

      outputstress(elem,estress,fout,func);
      free(estress);

    }
    if(wsurf.hwnd!=NULL)
    {

      drawyieldsurface((wsurf.childs+1)->hdcC,
                       (wsurf.childs+1)->vparam,
                       SURFACEX,SURFACEY,SURFACEZ,
                       af->fsurface);
      overlayhdc(*(wsurf.childs+1),SRCPAINT);     /*UPDATE DISPLAY.*/
    }
    /*while(!GetAsyncKeyState(VK_LBUTTON))
	;*/                                   /*LEFT CLICK TO CONTINUE.*/


/***UJIOKA:BUCKLING ANALYSIS***/
    laptime("ASSEMBLING GLOBAL MATRIX.",t0);

    estress=(double *)malloc(12*sizeof(double));

    for(i=1;i<=msize;i++)             /*GLOBAL MATRIX INITIALIZATION.*/
    {
      ginit.m=(unsigned int)i;
	  /*ginit.n=(unsigned int)i;*/
   	  *(kmtx+(i-1))=ginit;
      *(gmtx+(i-1))=ginit;
    }
//    comps=msize; /*INITIAL COMPONENTS=DIAGONALS.*/

    for(i=1;i<=nelem;i++)               /*ASSEMBLAGE GLOBAL MATRIX.*/
    {
      inputelem(elems,melem,i-1,&elem);        /*READ ELEMENT DATA.*/
//      elem=*(elems+i-1);                       /*READ ELEMENT DATA.*/
      for(ii=0;ii<=1;ii++)
      {
        for(jj=0;jj<=5;jj++)
        {
          (elems+i-1)->iconf[ii][jj]=elem.iconf[ii][jj];
        }
      }
      inputnode(ddisp,elem.node[0]);                         /*HEAD*/
      inputnode(ddisp,elem.node[1]);                         /*TAIL*/

      elem.sect=(elems+i-1)->sect;             /*READ SECTION DATA.*/

      if((wdraw.childs+1)->hdcC!=NULL)     /*DRAW DEFORMED ELEMENT.*/
      {
        drawglobalwire((wdraw.childs+1)->hdcC,
                       (wdraw.childs+1)->vparam,
                       *af,elem,255,255,255,
                                255,255,255,0,ONSCREEN);
      }

      drccos=directioncosine(elem.node[0]->d[0],
                             elem.node[0]->d[1],
                             elem.node[0]->d[2],
                             elem.node[1]->d[0],
                             elem.node[1]->d[1],
                             elem.node[1]->d[2],
                             elem.cangle);               /*[DRCCOS]*/

      tmatrix=transmatrix(drccos);         /*TRANSFORMATION MATRIX.*/
      estiff=assememtx(elem);          /*ELASTIC MATRIX OF ELEMENT.*/
      estiff=modifyhinge(elem,estiff);             /*MODIFY MATRIX.*/
      estiff=assempmtx(elem,estiff);          /*ADD PLASTIC MATRIX.*/
      estiff=transformation(estiff,tmatrix);       /*[K]=[Tt][k][T]*/

      assemgstiffness(kmtx,estiff,&elem);             /*ASSEMBLAGE.*/

      for(ii=0;ii<=11;ii++) free(*(estiff+ii));
      free(estiff);

      /*GEOMETRIC*/
      //estress=elemstress(&elem,gvct,melem,fout,func);
      for(ii=0;ii<=1;ii++)                                /*STRESSES.*/
      {
        for(jj=0;jj<6;jj++) *(estress+6*ii+jj)=elem.stress[ii][jj];
      }

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
      /*GEOMETRIC*/

      for(ii=0;ii<=2;ii++) free(*(drccos+ii));
      free(drccos);
      for(ii=0;ii<=11;ii++) free(*(tmatrix+ii));
      free(tmatrix);
    }
    sprintf(string,"GLOBAL MATRIX %ld COMPS ASSEMBLED.",comps);
    laptime(string,t0);

	bisecsylvester(gmtx,kmtx,confs,msize,neig,neig,biseceps,eigen,eigenvct);/*EIGEN VALUE PROBLEM*/
/***UJIOKA:BUCKLING ANALYSIS***/

	if(fout!=NULL) fprintf(fout,"\"REACTION\"\n");
	outputreaction(kmtx,gvct,nodes,confs,dreact,fout,nnode);

	updateform(ddisp,gvct,nnode);               /*FORMATION UPDATE.*/
	if(fout!=NULL) fprintf(fout,"\"CURRENT FORM\"\n");

	for(ii=0;ii<nnode;ii++)
	{
	  sprintf(string,"NODE:%5ld {U}=",(nodes+ii)->code);
	  for(jj=0;jj<6;jj++)
	  {
		loffset=6*ii+jj;
		sprintf(s," %14.5f",*(ddisp+loffset));
		strcat(string,s);
	  }
	  if(fout!=NULL) fprintf(fout,"%s\n",string);

	  if(ffig!=NULL && (nodes+ii)->code==108)
	  {
		fprintf(ffig," NODE:%5ld %s\n",(nodes+ii)->code,string);
	  }
	}

    for(i=0;i<neig;i++)
    {
      *(eigen+i)=1.0/(*(eigen+i));
      sprintf(string,"EIGEN VALUE %ld=%.5E",(i+1),*(eigen+i));
      fprintf(feig,"%s\n",string);
      errormessage(string);

      outputmode(*(eigenvct+i),feig,nnode,ninit);
	}

    if (*(eigen+0)<1.0)
    {
       sprintf(string,"BUCKLING at LAP:%ld",nlap);

       MessageBox(NULL,string,"BCLNG101",MB_OK);

       updatemode(af,*(eigenvct+0)); /*FORMATION UPDATE.*/
       break;
    }
    af->nlaps=neig;
    af->eigenval=eigen;
    af->eigenvec=eigenvct;


	t1=laptime("\0",t0);

	memory2=availablephysicalmemory(NULL);
	sprintf(string,"CONSUMPTION:%ld[BYTES]",(memory1-memory2));
	errormessage(string);

	errormessage("L:CONTINUE R:ABORT");            /*L=LEFT R=RIGHT*/
	while(!GetAsyncKeyState(VK_LBUTTON))  /*LEFT CLICK TO CONTINUE.*/
	{
	  if(GetAsyncKeyState(VK_RBUTTON))      /*RIGHT CLICK TO ABORT.*/
	  {
		fclose(fin);


		gfree(kmtx,nnode); /*FREE GLOBAL MATRIX.*/
		gfree(gmtx,nnode); /*FREE GLOBAL MATRIX.*/
		free(gvct);
        /*free(confs);*/

        errormessage(" ");
        errormessage("ABORTED.");
        if(fout!=NULL) fprintf(fout,"ABORTED.\n");

        fclose(fout);
        fclose(ffig);

        laptime("\0",t0);
        return 1;
      }
      t2=clock();
      time=(t2-t1)/CLK_TCK;
      if(time>=WAIT) break;               /*CONTINUE AFTER WAITING.*/
    }

    /*fclose(fout);*/
  }                                        /*REPEAT UNTIL INSTABLE.*/

  if((wdraw.childs+1)->hdcC!=NULL &&
     melem!=NULL && ddisp!=NULL)                 /*DRAW LAST FRAME.*/
  {
    for(i=1;i<=nelem;i++)
    {
      inputelem(elems,melem,i-1,&elem);
      for(ii=0;ii<=1;ii++) /*COPY HINGE DATA.*/
      {
        for(jj=0;jj<=5;jj++)
        {
          (elems+i-1)->iconf[ii][jj]=elem.iconf[ii][jj];
        }
      }

      inputnode(ddisp,elem.node[0]);
      inputnode(ddisp,elem.node[1]);

      drawglobalwire((wdraw.childs+1)->hdcC,
                     (wdraw.childs+1)->vparam,
                     *af,elem,255,255,255,
                              255,255,255,0,ONSCREEN);
	}
    overlayhdc(*(wdraw.childs+1),SRCPAINT);       /*UPDATE DISPLAY.*/
  }

  fclose(fin);


  gfree(kmtx,nnode); /*FREE GLOBAL MATRIX.*/
  gfree(gmtx,nnode); /*FREE GLOBAL MATRIX.*/  /*free(gvct);*/
  /*free(confs);*/

  af->eigenvec=(double **)malloc(1*sizeof(double *));
  *((af->eigenvec)+0)=gvct;

  errormessage(" ");
  errormessage("COMPLETED.");
  if(fout!=NULL) fprintf(fout,"COMPLETED.\n");

  if(fout!=NULL) fclose(fout);
  if(ffig!=NULL) fclose(ffig);
  if(ferr!=NULL) fclose(ferr);
  if(fsrf!=NULL) fclose(fsrf);
  if(feig!=NULL) fclose(feig);

  memory2=availablephysicalmemory("REMAIN:");
  sprintf(string,"CONSUMPTION:%ld[BYTES]",(memory0-memory2));
  errormessage(string);
  errormessage(" ");

  return 0;
}/*bclng101*/

int matrixdeterminant(struct gcomponent *gmtx,
            struct oconf *confs,
            long int msize,
            double *det,double *sign)
{
  /*char iconf;*/ /*0:FREE 1:FIXED*/
  char str[256];
  long int i,j,k;
  /*double det=1.0;*/
  double data1;
  struct gcomponent *pivot,*pcomp;

  *det=0.0;
  *sign=1.0;

  for(j=1;j<=msize;j++)                             /*DECOMPOSITION.*/
  {
    if((confs+j-1)->iconf==0) /*FREE*/
    {
      pivot=(gmtx+(j-1)); /*PIVOT.*/

      if((pivot->value)==0.0){*sign=0.0; return (j-1);}  /*INSTABLE.*/
      /*det*=pivot->value;*/
      *det+=log10(fabs(pivot->value));   /*LOG BY 10 OF DETERMINANT.*/
      *sign*=pivot->value/fabs(pivot->value); /*SIGN OF DETERMINANT.*/
      currentpivot(j,msize);
    }
  }
  if(*sign<=0.0) return 0;
  else return 1;
}/* croutlu */

