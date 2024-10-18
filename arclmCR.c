int arclmCR(struct arclmframe* af);

/*FOR VECTOR CALCURATION*/
double dotproduct(double* vct1, double* vct2, int vsize);
double vectorlength(double* vct, int vsize);
void vectornormalize(double* vct, int vsize);

/*FOR ROTATION CALCURATION*/
double* rotationvct(double** rmtx);
double** rotationmtx(double* rvct);
double** spinmtx(double* rvct);
double** jacobimtx(double* rvct);

double** spinfittermtx(double* eform, int nnod);/*G*/
double** projectionmtx(double* eform, double** G, int nnod);/*P*/
double** blockjacobimtx(double* edisp, double* estress, double** M, int nnod);/*H&M*/
double** transmatrixHPT(double* eform, double* edisp, double** T, int nnod);
double** assemtmtxCR(double** Ke, double* eform, double* edisp, double* estress, double* gstress, double** T, double** HPT, int nnod);
double** assemgmtxCR(double* eform, double* edisp, double* estress, double* gstress, double** T, double** HPT, int nnod);
void symmetricmtx(double** estiff, int msize);

void assemelem(struct owire* elems, struct memoryelem* melem, int nelem, long int* constraintmain,
			   struct gcomponent* mmtx, struct gcomponent* gmtx,
			   double* iform, double* ddisp, double* finternal, double* fexternal);
void assemshell(struct oshell* shells, struct memoryshell* mshell, int nshell, long int* constraintmain,
			   struct gcomponent* mmtx, struct gcomponent* gmtx,
			   double* iform, double* ddisp, double* finternal, double* fexternal);

void updaterotation(double* ddisp, double* gvct, int nnode);

double* quaternionvct(double* rvct);
double** updatedrccos(double** drccosinit, double* gforminit, double* gform);
double** interpolatermtx(double* rvct1, double* rvct2, double alpha);
//void updateforce(double* pvct, double* ddisp, int nnode);

double* extractlocalcoord(double* gform, double** drccos, double nnod);
double* extractshelldisplacement(struct oshell shell, double* ddisp);
double* extractdeformation(double* eforminit, double* eform, int nnod);

//void initialformCR(struct onode* nodes, double* ddisp, int nnode);
void LDLmode(struct gcomponent* gmtx, struct oconf* confs, int* m, int nmode, double** mode, double* norm, double* dm, long int msize);

double shellvolume(struct oshell shell, double** drccos, double area);
void assemshellvolume(struct oshell* shells, int nshell, double* ddisp, double* volume);

double equilibriumcurvature(double* weight, double* lapddisp, double laploadfactor, double* dup, int msize);

void dbgvct(double* vct, int size, int linesize, const char* str);
void dbgmtx(double** mtx, int rows, int cols, const char* str);

extern struct gcomponent *copygcompmatrix(struct gcomponent *gmtx,long int msize);
extern struct gcomponent *gcomponentadd3(struct gcomponent *mtx1,double factor1,
										 struct gcomponent *mtx2,double factor2,
										 int msize);
extern void bisecgeneral(struct gcomponent *A,double factorA,
						 struct gcomponent *B,double factorB,
						 struct oconf *confs,
						 long int N,long int NE,double defsign,
						 double EPS,
						 double *E,double **V,
						 double BL, double BR);
extern double inversemethod(struct gcomponent *gmtx, struct oconf *confs, double *evct, int msize);

/*
PINPOINTMODE
0:DO NOT DETECT BUCKLING
1:CALCULATE EIGEN VECTOR
*/



int arclmCR(struct arclmframe* af)
{
	DWORDLONG memory0,memory1;

	char dir[]=DIRECTORY;
	char string[400];
	char fname[256];
	clock_t t0;
	/*INPUT & OUTPUT FILE*/
	FILE * fin, * fonl, * fdsp, * fexf, * finf, * fubf, * frct, * fstr, * fene, * ffig, * fbcl, * feig, * fout;         /*FILE 8 BYTES*/

	int i, ii, jj;


	int nnode, nelem, nshell, nsect, nconstraint;
	long int msize;

	/*ARCLMFRAME*/
	struct osect* sects;
	struct onode* nodes;
	struct onode* ninit;
	struct owire* elems;
	struct oshell* shells;
	struct oconf* confs;
	struct memoryelem* melem;
	struct memoryshell* mshell;
	long int* constraintmain;

	/*GLOBAL MATRIX*/
	struct gcomponent ginit = { 0,0,0.0,NULL };
	struct gcomponent *epsgmtx, *gmtx, *gmtxcpy, *dgmtx, * g, * p, * gcomp1;
	double gg;
	double sign, determinant;
	long int nline;

	/*GLOBAL VECTOR*/
	double* ddisp, * iform;
	double* gvct;
	double* funbalance, * fexternal, * finternal, * freaction, * fpressure, * fbaseload, * fgivendisp, * fswitching;
    double* givendisp;

	///FOR ARC-LENGTH PARAMETER///
	int nlap, laps;
	int iteration;
	int maxiteration = 20;
	double tolerance = 1.0e-4;
	double residual;
	double loadfactor=0.0;
	double lambda;
	/*ARC-LENGTH METHOD*/
	double arclength;
	double arcsum, predictorsign;/*FOR PREDICTOR*/
	double dupdue, dupdup;       /*FOR CORRECTOR*/
	double* due, * dup;
	/*ARCLENGTH CONTROL*/
	int SCALINGARCFLAG = 0;
	double k1, k, scaledarclength;
	/*BIGINING ARC-LENGTH*/
	int BIGININGARCFLAG = 0;
	double biginingarcratio = 1.0;
	int biginingarclap = 0;
	/*WEIGHT*/
	double *weight;
	int node;

	/*OUTPUT*/
	int outputmode = 0;/*0:ConvergedLaps.1:AllLaps.*/
	double volume = 0.0;


	///FOR PINPOINTING///
	int pinpointmode = 0;/*0:NotPinpointing.1:BisecPinpointing.2:ExtendedSystemPinpointing.*/
	int BCLFLAG = 0;/*FLAG FOR BUCKLING DETECTION*/
	/*LAST LAP*/
	double* lastddisp,* lastgvct,* lastpivot;
	double lastloadfactor,lastlambda,lastsign;
	/*INCREMENT OF CURRENT LAP*/
	double* lapddisp;
	double laploadfactor;
	/*FOR BISECSYLVESTER*/
	double LL,LR,LM;
	double* nextddisp;
	double nextloadfactor;
	double biseceps = 1.0e-5;
	/*FOR INVERSE ITERATION*/
	double biseceigen;
	double* bisecevct;
	double evctdot;/*FOR BIFURCATION DETECTION*/

	/*FOR EXTENDED SYSTEM (LDL MODE & PATH SWITCHING)*/
	int ldlneig = 0;
	double** ldlevct;   /*NORMALIZED EIGEN MODE*/
	double* norm, * dm;	/*NORM & PIVOT OF LDL-MODE*/
	int* m;            	/*LINE OF LDL-MODE*/
	double* epsddisp, * epsgvct;
	double* re, * rp;
	double gradeps = 1.0e-8;
	double gamma;
	double evctfunbalance, epsevctfunbalance;
	double * epsfunbalance, * epsfinternal, * epsfexternal, * epsfpressure;

	int neig = 0;
	double eps=1.0e-3;
	double **evct, *eigen;

	/*ANALYSIS TERMINATION*/
	int ENDFLAG = 0;
	int fnode=NULL,fnodedrc=NULL;
	double fnodemin, fnodemax;
	/*UNLOADING*/
	int UNLOADFLAG = 0;


	/*FOR READING ANALYSISDATA*/
	FILE *fdata;
	int nstr, pstr, readflag;
	char** data;
	char filename[256];
	char* dot;

	memory0 = availablephysicalmemoryEx("INITIAL:");   /*MEMORY AVAILABLE*/


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

	msize = 6 * nnode;                                      /*SIZE OF GLOBAL MATRIX.*/

	weight = (double *)malloc((6*nnode + 1) * sizeof(double));
	for (i=0;i<6*nnode;i++)
	{
		*(weight + i) = 0.0;
	}
	*(weight + 6*nnode) = 1.0;



	strcpy(filename, (wdraw.childs+1)->inpfile);
	dot = strrchr(filename, '.');
	*dot = '\0';

	fdata = fgetstofopenII(dir, "r", "analysisdata.txt");

	if (fdata == NULL)
	{
		errormessage("couldn't open analysisdata.txt\n");
		getincrement((wmenu.childs+2)->hwnd,&laps,&arclength);
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
								getincrement((wmenu.childs+2)->hwnd,&laps,&arclength);
							}
							else
							{
								strcpy(filename, *(data + pstr));
                            }
						}
						if (!strcmp(*(data + pstr), "LAPS"))
						{
							pstr++;
							laps = (int)strtol(*(data + pstr), NULL, 10);
						}
						if (!strcmp(*(data + pstr), "ITERMAX"))
						{
							pstr++;
							maxiteration = (int)strtol(*(data + pstr), NULL, 10);
						}
						if (!strcmp(*(data + pstr), "ARCLENGTH"))
						{
							pstr++;
							arclength = (double)strtod(*(data + pstr), NULL);
						}
						if (!strcmp(*(data + pstr), "SCALINGARC"))
						{
							SCALINGARCFLAG=1;
						}
						if (!strcmp(*(data + pstr), "BIGININGARC"))
						{
							BIGININGARCFLAG=1;
							pstr++;
							biginingarcratio = (double)strtod(*(data + pstr), NULL);
							pstr++;
							biginingarclap = (int)strtol(*(data + pstr), NULL, 10);
						}
						if (!strcmp(*(data + pstr), "WEIGHT"))
						{
							pstr++;
							if (!strcmp(*(data + pstr), "LOAD"))
							{
								pstr++;
								*(weight + 6 * nnode) = (double)strtod(*(data + pstr), NULL);
							}
							if (!strcmp(*(data + pstr), "NODE"))
							{
								pstr++;
								node= (int)strtol(*(data + pstr), NULL, 10);
								node -= 101;
								pstr++;
								i = (int)strtol(*(data + pstr), NULL, 10);
								pstr++;
								*(weight+6*node+i) = (double)strtod(*(data + pstr), NULL);
							}
						}
						if (!strcmp(*(data + pstr), "OUTPUTMODE"))
						{
							pstr++;
							outputmode = (int)strtol(*(data + pstr), NULL, 10);
						}
						if (!strcmp(*(data + pstr), "PINPOINTMODE"))
						{
							pstr++;
							pinpointmode = (int)strtol(*(data + pstr), NULL, 10);
						}
						if (!strcmp(*(data + pstr), "FNODE"))
						{
							pstr++;
							fnode = (int)strtol(*(data + pstr), NULL, 10);
							pstr++;
							fnodedrc = (int)strtol(*(data + pstr), NULL, 10);
							pstr++;
							fnodemin = (double)strtod(*(data + pstr), NULL);
							pstr++;
							fnodemax = (double)strtod(*(data + pstr), NULL);

						}
						if (!strcmp(*(data + pstr), "LOADFACTOR"))
						{
							pstr++;
							loadfactor = (double)strtod(*(data + pstr), NULL);
						}
                        if (!strcmp(*(data + pstr), "NEIG"))
						{
							pstr++;
							neig = (int)strtol(*(data + pstr), NULL, 10);
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

	sprintf(string, "FILENAME : %s\n LAPS = %d\n MAX ITERATION= %d\n ARCLENGTH = %lf\n", filename, laps, maxiteration, arclength);
	errormessage(string);

	if (outputmode == 0)errormessage("OUTPUT CONVERGED RESULT\n");
	if (outputmode == 1)errormessage("OUTPUT ALL RESULT\n");


	snprintf(fname, sizeof(fname), "%s.%s", filename, "dsp");
	fdsp = fopen(fname, "w");
	snprintf(fname, sizeof(fname), "%s.%s", filename, "inf");
	finf = fopen(fname, "w");
	snprintf(fname, sizeof(fname), "%s.%s", filename, "exf");
	fexf = fopen(fname, "w");
	snprintf(fname, sizeof(fname), "%s.%s", filename, "ubf");
	fubf = fopen(fname, "w");
	snprintf(fname, sizeof(fname), "%s.%s", filename, "rct");
	frct = fopen(fname, "w");
	snprintf(fname, sizeof(fname), "%s.%s", filename, "str");
	fstr = fopen(fname, "w");
	snprintf(fname, sizeof(fname), "%s.%s", filename, "onl");
	fonl = fopen(fname, "w");
	snprintf(fname, sizeof(fname), "%s.%s", filename, "fig");
	ffig = fopen(fname, "w");
	snprintf(fname, sizeof(fname), "%s.%s", filename, "bcl");
	fbcl = fopen(fname, "w");
	snprintf(fname, sizeof(fname), "%s.%s", filename, "ene");
	fene = fopen(fname, "w");
	snprintf(fname, sizeof(fname), "%s.%s", filename, "eig");
	feig = fopen(fname, "w");
	snprintf(fname, sizeof(fname), "%s.%s", filename, "otl");
	fout = fopen(fname, "w");

	t0 = clock();                                                   /*CLOCK BEGIN.*/


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

	/***GLOBAL MATRIX***/
	gmtx = (struct gcomponent*)malloc(msize * sizeof(struct gcomponent));/*DIAGONALS OF GLOBAL MATRIX.*/
	epsgmtx = (struct gcomponent*)malloc(msize * sizeof(struct gcomponent));/*DIAGONALS OF GLOBAL MATRIX.*/

	/***GLOBAL VECTOR***/
	gvct = (double *)malloc(msize * sizeof(double));/*INCREMENTAL GLOBAL VECTOR.*/

	funbalance = (double*)malloc(msize * sizeof(double));          /*UNBALANCED INTERNAL FORCE VECTOR.*/
	freaction = (double*)malloc(msize * sizeof(double));           /*EXTERNAL FORCE VECTOR.*/
	fexternal = (double*)malloc(msize * sizeof(double));           /*EXTERNAL FORCE VECTOR.*/
	finternal = (double*)malloc(msize * sizeof(double));           /*INTERNAL FORCE VECTOR.*/
	fpressure = (double*)malloc(msize * sizeof(double));           /*PRESSURE VECTOR.*/
	fbaseload = (double*)malloc(msize * sizeof(double));           /*BASE LOAD VECTOR.*/
	fgivendisp = (double*)malloc(msize * sizeof(double));           /*BASE LOAD VECTOR.*/
	fswitching = (double*)malloc(msize * sizeof(double));


	due = (double*)malloc(msize * sizeof(double));
	dup = (double*)malloc(msize * sizeof(double));

	lastddisp = (double*)malloc(msize * sizeof(double));
	lastgvct = (double*)malloc(msize * sizeof(double));
	lastpivot = (double*)malloc(msize * sizeof(double));    /*PIVOT SIGN OF TANGENTIAL STIFFNESS.*/
	lapddisp = (double*)malloc(msize * sizeof(double));		/*INCREMENT IN THE LAP*/
	nextddisp = (double*)malloc(msize * sizeof(double));
	givendisp = (double*)malloc(msize * sizeof(double));

	bisecevct = (double*)malloc(msize * sizeof(double));

	epsddisp = (double*)malloc(6 * nnode * sizeof(double));
	epsgvct = (double*)malloc(6 * nnode * sizeof(double));

	epsfunbalance = (double*)malloc(6 * nnode * sizeof(double));
	epsfinternal = (double*)malloc(6 * nnode * sizeof(double));
	epsfexternal = (double*)malloc(6 * nnode * sizeof(double));
	epsfpressure = (double*)malloc(6 * nnode * sizeof(double));
	re = (double*)malloc(msize * sizeof(double));
	rp = (double*)malloc(msize * sizeof(double));

	evct=(double **)malloc(neig*sizeof(double *));   /*GLOBAL VECTORS*/
	eigen=(double *)malloc(neig*sizeof(double));       /*EIGEN VALUES*/
	for(i=0;i<neig;i++)
	{
	  *(evct+i)=(double *)malloc(msize*sizeof(double));
	  for(ii=0;ii<msize;ii++)
	  {
		*(*(evct+i)+ii)=0.0;
	  }
	}

	for (i = 0; i < msize; i++)
	{
		(gmtx + i)->down = NULL;
        (epsgmtx + i)->down = NULL;
		*(gvct + i) = 0.0;
	}


	initialelem(elems, melem, nelem);             /*ASSEMBLAGE ELEMENTS.*/
	initialshell(shells, mshell, nshell);         /*ASSEMBLAGE ELEMENTS.*/

	setviewpoint((wdraw.childs+0)->hwnd,*af,&((wdraw.childs+1)->vparam));
	setviewparam((wmenu.childs+2)->hwnd,(wdraw.childs+1)->vparam);

	GetAsyncKeyState(VK_LBUTTON);                   /*CLEAR KEY LEFT.*/
	GetAsyncKeyState(VK_RBUTTON);                  /*CLEAR KEY RIGHT.*/
	if(globaldrawflag==1)
	{
	  drawglobalaxis((wdraw.childs+1)->hdcC,(wdraw.childs+1)->vparam,0,0,255);                     /*DRAW GLOBAL AXIS.*/
	}

	nlap = 1;
	iteration = 1;
	residual = 0.0;

	assemconf(confs,fbaseload,1.0,nnode);               /*GLOBAL VECTOR.*/
	assemgivend(confs,givendisp,1.0,nnode);


	while (nlap <= laps && ENDFLAG == 0)
	{
		af->nlaps = nlap;

		if ((outputmode == 0 && iteration == 1) || outputmode == 1)
		{
			sprintf(string, "LAP: %5ld / %5ld ITERATION: %5ld\n", nlap, laps, iteration);
			fprintf(fdsp, string);
			fprintf(finf, string);
			fprintf(fexf, string);
			fprintf(fubf, string);
			fprintf(frct, string);
			fprintf(fstr, string);
			fprintf(fene, string);
		}

		for (i = 1; i <= msize; i++)/*MATRIX & VECTOR INITIALIZATION.*/
		{
			g = (gmtx + (i - 1))->down;   /*NEXT OF DIAGONAL.*/
			while (g != NULL) 	      /*CLEAR ROW.*/
			{
				p = g;
				g = g->down;
				free(p);
			}
			ginit.m = (unsigned short int)i;
			*(gmtx + (i - 1)) = ginit;
		}
		for (i = 1; i <= msize; i++)/*MATRIX & VECTOR INITIALIZATION.*/
		{
			g = (epsgmtx + (i - 1))->down;   /*NEXT OF DIAGONAL.*/
			while (g != NULL) 	      /*CLEAR ROW.*/
			{
				p = g;
				g = g->down;
				free(p);
			}
			ginit.m = (unsigned short int)i;
			*(epsgmtx + (i - 1)) = ginit;
		}

		for (i = 0; i < msize; i++)/*MATRIX & VECTOR INITIALIZATION.*/
		{
			*(finternal + i) = 0.0;			 /*GLOBAL VECTOR INITIALIZATION.*/
			*(fexternal + i) = 0.0;
			*(fpressure + i) = 0.0;
			*(funbalance + i) = 0.0;
			*(freaction + i) = 0.0;
		}
		comps = msize; /*INITIAL COMPONENTS = DIAGONALS.*/

		assemshellvolume(shells, nshell, ddisp, &volume);



		/*ELEMENT STIFFNESS & FORCE ASSEMBLAGE*/
		assemelem(elems, melem, nelem, constraintmain, NULL, gmtx, iform, ddisp, finternal, fpressure);
		assemshell(shells, mshell, nshell, constraintmain, NULL, gmtx, iform, ddisp, finternal, fpressure);

		if(iteration==1)
		{
		  clearwindow(*(wdraw.childs+1));
		  drawarclmframe((wdraw.childs+1)->hdcC,(wdraw.childs+1)->vparam,*af,0,ONSCREEN);
		  overlayhdc(*(wdraw.childs + 1), SRCPAINT);                  /*UPDATE DISPLAY.*/
		}

		/*EXTERNAL FORCE & UNBALANCED FORCE*/
		if(iteration==1)
		{
			for (i = 0; i < msize; i++)/*IMPOSED DISP.*/
			{
				*(fgivendisp + i) = *(givendisp + i);
			}
			modifygivend(gmtx,fgivendisp,confs,nnode);      /*0:LOAD 1:DISPLACEMENT.*/
		}
		residual = 0.0;
		if(/*UNLOADFLAG==1*/0)
		{
			for (i = 0; i < msize; i++)
			{
				if(nlap==1 && iteration==1)
				{
					*(fbaseload + i) = -*(finternal + i);
				}
				*(dup + i) = *(fbaseload + i)+*(fpressure + i);
				*(fexternal + i) = - *(fbaseload + i) + loadfactor * *(dup + i);
				*(funbalance + i) = *(fexternal + i) - *(finternal + i);            /*funbalance:UNBALANCED FORCE -{E}.*/
				if ((confs + i)->iconf == 1)
				{
					*(freaction + i) = *(funbalance + i);
					*(funbalance + i) = 0.0;
					*(dup + i) = 0.0;
				}
				residual += *(funbalance + i) * *(funbalance + i);
				*(due + i) = *(funbalance + i);
			}

		}
		else
		{
			for (i = 0; i < msize; i++)
			{
				*(dup + i) = *(fbaseload + i)+*(fpressure + i);
				*(fexternal + i) = loadfactor * *(dup + i);
				*(funbalance + i) = *(fexternal + i) - *(finternal + i);            /*funbalance:UNBALANCED FORCE -{E}.*/
				if ((confs + i)->iconf == 1)
				{
					*(freaction + i) = *(funbalance + i);
					*(funbalance + i) = 0.0;
					*(dup + i) = 0.0;
				}
				residual += *(funbalance + i) * *(funbalance + i);
				*(due + i) = *(funbalance + i);
				if(iteration==1)*(dup + i) += *(fgivendisp + i);
			}
		}

		/*OUTPUT*/
		if ((outputmode == 0 && iteration == 1) || outputmode == 1)
		{
			outputdisp(ddisp, fdsp, nnode, nodes);                    /*FORMATION OUTPUT.*/
			outputdisp(finternal, finf, nnode, nodes);                /*FORMATION OUTPUT.*/
			outputdisp(fexternal, fexf, nnode, nodes);                /*FORMATION OUTPUT.*/
			outputdisp(funbalance, fubf, nnode, nodes);               /*FORMATION OUTPUT.*/
			outputdisp(freaction, frct, nnode, nodes);                /*FORMATION OUTPUT.*/

			/*STRESS & ENERGY OUTPUT*/
			for(i = 0; i < nshell; i++)
			{
				fprintf(fene, "%5ld %e %e %e\n", (shells+i)->code, (mshell+i)->SEp, (mshell+i)->SEb, (mshell+i)->SE);
				fprintf(fstr, "%5ld %e %e %e %e %e %e\n", (shells+i)->code,
				(mshell+i)->stress[0][0],
				(mshell+i)->stress[0][1],
				(mshell+i)->stress[0][2],
				(mshell+i)->stress[0][3],
				(mshell+i)->stress[0][4],
				(mshell+i)->stress[0][5]);
			}
		}


		/*CHECK BEFORE FORWARDBACKWARD : dup & due == 0.0 IF ICONF == 1*/



		if (BCLFLAG < 1)/*REGULAR*/
		{
			if(iteration==1 && neig!=0)gmtxcpy=copygcompmatrix(gmtx,msize);

			nline = croutlu(gmtx, confs, msize, &determinant, &sign, gcomp1);
			//ldlneig = 0;
			//for (ii = 0; ii < msize; ii++)
			//{
			//	if ((confs + ii)->iconf == 0 && *(lastpivot + ii) * ((gmtx + ii)->value) < 0)
			//	{
			//		ldlneig++;
			//	}
			//	if ((confs + ii)->iconf == 0 && ((gmtx + ii)->value) < 0 && iteration==1)
			//	{
			//		errormessage("%d\n ", ii);
			//	}
			//}

			sprintf(string, "LAP: %4d ITER: %2d {LOAD}= % 5.8f {RESD}= %1.6e {DET}= %8.5f {SIGN}= %2.0f {BCL}= %1d {EPS}= %1.5e {V}= %8.5f\n",
				nlap, iteration, loadfactor, residual, determinant, sign, BCLFLAG, 0.0, volume);
			fprintf(ffig, "%s", string);
			errormessage(string);
			/*LDL DECOMPOSITION FAILED*/
			if (sign < 0.0)
			{
				for (ii = 1; ii <= msize; ii++)
				{
					gg = 0.0;
					gread(gmtx, ii, ii, &gg);

					if (gg < 0.0)
					{
						sprintf(string, "INSTABLE TERMINATION AT NODE %ld.",
							(nodes + int((ii - 1) / 6))->code);
						errormessage(" ");
						errormessage(string);
						if (fonl != NULL) fprintf(fonl, "%s\n", string);
					}
				}

				fclose(fonl);
				fclose(fdsp);
				fclose(fexf);
				fclose(finf);
				fclose(fubf);
				fclose(frct);
				fclose(fstr);
				fclose(fene);
				fclose(ffig);
				fclose(fbcl);

				gfree(gmtx,nnode);
				gfree(epsgmtx,nnode);
				free(weight);
				free(gvct);
				free(funbalance);
				free(freaction);
				free(fexternal);
				free(finternal);
				free(fpressure);
				free(fbaseload);
				free(fgivendisp);
				free(fswitching);
				free(due);
				free(dup);
				free(lastddisp);
				free(lastgvct);
				free(lastpivot);
				free(lapddisp);
				free(nextddisp);
				free(bisecevct);
				free(epsddisp);
				free(epsgvct);
				free(epsfunbalance);
				free(epsfinternal);
				free(epsfexternal);
				free(epsfpressure);
				free(re);
				free(rp);

				return 1;
			}


			if (iteration == 1 && (sign - lastsign) != 0 && nlap != 1 && pinpointmode != 0 && BCLFLAG == 0)/*BUCKLING DETECTED*/
			{
				BCLFLAG = 1;
				sprintf(string, "BUCKLING DITECTED LAP: %4d ITER: %2d\n", nlap, iteration);
				fprintf(fbcl, "%s", string);

				if (pinpointmode == 1)/*GO TO BISEC*/
				{
					LR = 1.0;
					LL = 0.0;
					LM = 0.5;
					/*NEW TARGET*/
					nextloadfactor = loadfactor;
					loadfactor = 0.5 * (nextloadfactor + lastloadfactor);
					for (ii = 0; ii < msize; ii++)
					{
						*(nextddisp + ii) = *(ddisp + ii);
						*(ddisp + ii) = 0.5 * (*(nextddisp + ii) + *(lastddisp + ii));
					}
				}
				if(pinpointmode == 2)/*EXTENDED SYSTEM*/
				{
					ldlneig = 0;
					for (ii = 0; ii < msize; ii++)
					{
						if ((confs + ii)->iconf == 0 && *(lastpivot + ii) * ((gmtx + ii)->value) < 0)
						{
							ldlneig++;
						}
					}

					norm = (double *)malloc(ldlneig * sizeof(double));/*NORM OF LDL MODE*/
					dm = (double *)malloc(ldlneig * sizeof(double));/*DIAGONAL PIVOT VALUE*/
					m = (int*)malloc(ldlneig * sizeof(int));/*BUCKLING DETECTED LINE*/

					ldlevct =  (double **)malloc(ldlneig * sizeof(double *));
					for (ii = 0; ii < ldlneig; ii++)
					{
					   *(ldlevct + ii) = (double *)malloc(msize * sizeof(double));
					}

					jj=0;
					for (ii = 0; ii < msize; ii++)
					{
						if ((confs + ii)->iconf == 0 && *(lastpivot + ii) * ((gmtx + ii)->value) < 0)
						{
							*(m + jj) = ii;
							jj++;
							sprintf(string, "BUCKLING DITECTED LAP: %4d ITER: %2d LINE: %5ld ", nlap, iteration, ii);
							fprintf(fbcl, "%s\n", string);
							fprintf(feig, "%s", string);
							if (jj == ldlneig)
							{
								break;
							}
						}
					}
					/*PIVOT ZERO LINE CHECK*/
				}

				iteration++;
			}
			else/*NORMAL CALCULATION IN CASE OF REGULAR*/
			{
				if (iteration == 1)/*PREDICTOR CALCULATION*/
				{
					nline = forwardbackward(gmtx, dup, confs, msize, gcomp1);
					/*ARCLENGTH CONTROL*/
					if (nlap < biginingarclap+1)
					{
						scaledarclength = biginingarcratio * arclength;
					}
					else if (nlap == biginingarclap+1)
					{
						scaledarclength = arclength;
					}
					else if (nlap == biginingarclap+2)
					{
						if(SCALINGARCFLAG==1)
						{
							k1 = equilibriumcurvature(weight, lapddisp, laploadfactor, dup, msize);
						}
					}
					else
					{
						if(SCALINGARCFLAG==1)
						{
							k = equilibriumcurvature(weight, lapddisp, laploadfactor, dup, msize);
							scaledarclength = arclength * sqrt(k1 / k);
							if (scaledarclength > arclength)scaledarclength = arclength;
						}
					}

					/*SIGN OF PREDICTOR VECTOR & ARC-LENGTH SCALING FACTOR*/
					if(nlap == 1)
					{
						for (ii = 0; ii < msize; ii++)
						{
							if ((confs + ii)->iconf != 1 && *(weight + ii) != 0.0 && *(dup + ii)!=0.0)
							{
								*(weight + ii) = 1.0 / abs(*(dup + ii));
							}
							else
							{
								*(weight + ii) = 0.0;
                            }
						}
					}

					if (nlap==1)
					{
						predictorsign = 1.0;
					}
					else
					{
						predictorsign = *(weight + msize) * *(weight + msize) * laploadfactor * 1.0;
						for (ii = 0; ii < msize; ii++)
						{
							predictorsign += *(weight + ii) * *(weight + ii) * *(lapddisp + ii) * *(dup + ii);
						}

						if(predictorsign != 0.0)
						{
							predictorsign /= abs(predictorsign);
						}
						else
						{
							predictorsign = 1.0;
						}
					}

					/*INCREMENTAL CALCULATION*/
					arcsum = *(weight + msize) * *(weight + msize);
					for (ii = 0; ii < msize; ii++)
					{
						arcsum += *(weight + ii) * *(weight + ii) * *(dup + ii) * *(dup + ii);/*SCALED DISPLACEMENT.*/
					}
					lambda = predictorsign * scaledarclength / sqrt(arcsum);/*INCREMANTAL LOAD FACTOR.*/
					for (ii = 0; ii < msize; ii++)
					{
						*(gvct + ii) = lambda * (*(dup + ii)+*(givendisp + ii));/*INCREMANTAL DISPLACEMENT.*/
					}
					fprintf(fonl, "LAP: %4d ITER: %2d {LAMBDA}= % 5.8f {TOP}= % 5.8f {BOTTOM}= % 5.8f\n", nlap, iteration, lambda, scaledarclength, sqrt(arcsum));


					/*MEMORY OF PREDICTOR*/
					for (ii = 0; ii < msize; ii++)
					{
						*(lastddisp + ii) = *(ddisp + ii);/*NODE COORDINATION AT CONVERGED POINT*/
						*(lastgvct + ii) = *(gvct + ii);/*INCREMENTAL DISPLACEMENT OF PREDICTOR*/
						*(lastpivot + ii) = (gmtx + ii)->value;/*PIVOT SIGN OF TANGENTIAL STIFFNESS*/
					}
					lastloadfactor = loadfactor;/*LOAD FACTOR AT CONVERGED POINT*/
					lastlambda = lambda;/*INCREMENTAL LOAD FACTOR OF PREDICTOR*/
					lastsign = sign;/*SIGN AT CONVERGED POINT*/
					BCLFLAG = 0;

				}
				else if (BCLFLAG == -2)/*CORRECTOR CALCULATION*/
				{
					nline = forwardbackward(gmtx, due, confs, msize, gcomp1);
					lambda = 0.0;
					for (ii = 0; ii < msize; ii++)
					{
						*(gvct + ii) = *(due + ii);/*gvct:{ƒÂU_e+ƒÂƒ©U_p}*/
					}
					fprintf(fonl, "LAP: %4d ITER: %2d {LAMBDA}= % 5.8f {TOP}= % 5.8f {BOTTOM}= % 5.8f\n", nlap, iteration, lambda, 0.0, 0.0);
				}
				else/*CORRECTOR CALCULATION*/
				{
					nline = forwardbackward(gmtx, dup, confs, msize, gcomp1);
					nline = forwardbackward(gmtx, due, confs, msize, gcomp1);

					#if 1
					/*MINIMUM RESIDUAL QUANTITIES METHOD*/
					dupdue = 0.0;
					dupdup = *(weight + msize) * *(weight + msize);
					for (ii = 0; ii < msize; ii++)
					{
						dupdue += *(weight + ii) * *(weight + ii) * *(dup + ii) * *(due + ii);
						dupdup += *(weight + ii) * *(weight + ii) * *(dup + ii) * *(dup + ii);
						//dupdue += *(dup + ii) * *(due + ii);
						//dupdup += *(dup + ii) * *(dup + ii);
					}
					lambda = -dupdue / dupdup;
					for (ii = 0; ii < msize; ii++)
					{
						*(gvct + ii) = lambda * *(dup + ii) + *(due + ii);
					}
					fprintf(fonl, "LAP: %4d ITER: %2d {LAMBDA}= % 5.8f {TOP}= % 5.8f {BOTTOM}= % 5.8f\n", nlap, iteration, lambda, dupdue, dupdup);
					#endif

					#if 0
					/*HYPERSPHERE CONSTRAINT WITH RADIAL RETURN*/
					/*USING PREDICTOR DATA*/
					dupdue = 0.0;
					dupdup = *(weight + msize) * *(weight + msize) * lambda;
					for (ii = 0; ii < msize; ii++)
					{
						dupdue += *(weight + ii) * *(weight + ii) * *(gvct + ii) * *(due + ii);
						dupdup += *(weight + ii) * *(weight + ii) * *(gvct + ii) * *(dup + ii);
					}
					lambda += - dupdue / dupdup;
					fprintf(fonl, "LAP: %4d ITER: %2d {LAMBDA}= %e {TOP}= %e {BOTTOM}= %e\n", nlap, iteration, lambda, dupdue, dupdup);
					for (ii = 0; ii < msize; ii++)
					{
						*(gvct + ii) += -(dupdue / dupdup) * *(dup + ii) + *(due + ii);
					}
					arcsum = *(weight + msize) * *(weight + msize) * lambda * lambda;
					for (ii = 0; ii < msize; ii++)
					{
						arcsum += *(weight + ii) * *(weight + ii) * *(gvct + ii) * *(gvct + ii);
					}
					lambda *= scaledarclength / sqrt(arcsum);
					loadfactor = lastloadfactor;
					for (ii = 0; ii < msize; ii++)
					{
						*(gvct + ii) *= scaledarclength / sqrt(arcsum);
						*(ddisp + ii) = *(lastddisp + ii);
					}
					#endif
				}

				if(iteration == 1  && neig!=0)/*TWO-POINTS BUCKLING ANALYSIS.*/
				{
					eps=1e-8;
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
					assemelem(elems, NULL, nelem, constraintmain, NULL, epsgmtx, iform, epsddisp, NULL, NULL);
					assemshell(shells, NULL, nshell, constraintmain, NULL, epsgmtx, iform, epsddisp, NULL, NULL);

					dgmtx=gcomponentadd3(epsgmtx,1.0/eps,gmtxcpy,-1.0/eps,msize);

					bisecgeneral(dgmtx,-1.0,gmtxcpy,-1.0,confs,msize,neig,sign,1.0E-10,eigen,evct,0.0,1.0E+3);

					/*OUTPUT*/
					if(feig!=NULL)
					{
						for(i=0;i<neig;i++)
						{
							fprintf(feig,"LAP = %d MODE = %d GENERALIZED EIGENVALUE = %e STANDARD EIGENVALUE = %e ", af->nlaps, (i+1), *(eigen + i), 0.0);

							if (*(eigen + i) > 0.0)
							{
								*(eigen+i) = 1.0/(*(eigen + i));
								fprintf(feig, "LAMBDA = %e BUCKLINGLOAD = %e\n",*(eigen+i),loadfactor+*(eigen+i));
							}
							else
							{
								fprintf(feig, "ERROR:EIGEN VALUE NEGATIVE.\n");
							}

							for (ii = 0; ii < msize; ii++)
							{
								*(*(evct + i) + ii) = *(*(evct + i) + *(constraintmain + ii));
							}
							for (ii = 0; ii < nnode; ii++)
							{
								fprintf(feig,
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
					gfree(gmtxcpy,nnode);
					//gfree(epsgmtx,nnode);
					gfree(dgmtx,nnode);
				}


				for (ii = 0; ii < msize; ii++)
				{
					*(gvct + ii) = *(gvct + *(constraintmain + ii));
				}

				/*lapddisp : INCREMENTAL TRANSITION & ROTATION IN THIS LAP.*/
				if(iteration == 1)
				{
					laploadfactor=0.0;
					for(i=0;i<msize;i++)
					{
						*(lapddisp+i)=0.0;
					}
				}
				laploadfactor += lambda;
				updaterotation(lapddisp, gvct, nnode);

				loadfactor += lambda;
				updaterotation(ddisp, gvct, nnode);


				if ((residual < tolerance || iteration > maxiteration) && iteration != 1)
				{
					nlap++;
					iteration = 0;
				}
				else if (maxiteration == 1)
				{
					nlap++;
					iteration = 0;
				}
				iteration++;
			}
		}
		else if (BCLFLAG == 1)/*BUCKLING DETECTED.PIN-POINTING EXCUTION.*/
		{
			nline = croutlu(gmtx, confs, msize, &determinant, &sign, gcomp1);/*FOT COUNTING NEGATIVE PIVOT*/
			sprintf(string, "LAP: %4d ITER: %2d {LOAD}= % 5.8f {RESD}= %1.6e {DET}= %8.5f {SIGN}= %2.0f {BCL}= %1d {EPS}=%1.5e {V}= %8.5f\n",
				nlap, iteration, loadfactor, residual, determinant, sign, BCLFLAG, LM, volume);
			fprintf(ffig, "%s", string);
			fprintf(fbcl, "%s", string);
			errormessage(string);

			/*LDL DECOMPOSITION FAILED*/
			if (sign < 0.0)
			{
				for (ii = 1; ii <= msize; ii++)
				{
					gg = 0.0;
					gread(gmtx, ii, ii, &gg);

					if (gg < 0.0)
					{
						sprintf(string, "INSTABLE TERMINATION AT NODE %ld.",
							(nodes + int((ii - 1) / 6))->code);
						errormessage(" ");
						errormessage(string);
						if (fonl != NULL) fprintf(fonl, "%s\n", string);
					}
				}

				fclose(fin);
				fclose(fonl);
				fclose(fdsp);
				fclose(fexf);
				fclose(finf);
				fclose(fubf);
				fclose(frct);
				fclose(fstr);
				fclose(fene);
				fclose(ffig);
				fclose(fbcl);

				gfree(gmtx,nnode);
                gfree(epsgmtx,nnode);
				free(weight);
				free(gvct);
				free(funbalance);
				free(freaction);
				free(fexternal);
				free(finternal);
				free(fpressure);
				free(fbaseload);
				free(fgivendisp);
				free(fswitching);
				free(due);
				free(dup);
				free(lastddisp);
				free(lastgvct);
				free(lastpivot);
				free(lapddisp);
				free(nextddisp);
				free(bisecevct);
				free(epsddisp);
				free(epsgvct);
				free(epsfunbalance);
				free(epsfinternal);
				free(epsfexternal);
				free(epsfpressure);
				free(re);
				free(rp);

				return 1;
			}


			if (pinpointmode == 1)/*PIN-POINTING BASED ON SYLVESTER'S LAW OF INERTIA (BISECTION METHOD).*/
			{
				/*BISECTION PIN-POINTING*/
				if (lastsign == sign)
				{
					LL = LM;
				}
				else
				{
					LR = LM;
				}
				LM = 0.5 * (LL + LR);
				loadfactor = (1 - LM) * lastloadfactor + LM * nextloadfactor;
				for (ii = 0; ii < msize; ii++)
				{
					*(ddisp + ii) = (1 - LM) * *(lastddisp + ii) + LM * *(nextddisp + ii);
				}

				if (LR - LL < biseceps)/*CONVERGED*/
				{
					/*INITIAL EIGEN VECTOR*/
					for (ii = 0; ii < msize; ii++)
					{
						*(bisecevct + ii) = *(due + ii);
					}
					/*INVERS METHOD FOR EIGEN VECTOR*/
					biseceigen=inversemethod(gmtx,confs,bisecevct,msize);

					/*BIFURCATION JUDGE*/
					vectornormalize(fexternal, msize);
					evctdot = dotproduct(bisecevct, fexternal, msize);

					fprintf(feig,"LAP = %d MODE = %d GENERALIZED EIGENVALUE = %e STANDARD EIGENVALUE = %e LAMBDA = %e DOT = %e\n", nlap, 1, LM, biseceigen, loadfactor-lastloadfactor, evctdot);
					for (ii = 0; ii < msize; ii++)
					{
						*(bisecevct + ii) = *(bisecevct + *(constraintmain + ii));
					}
					for (ii = 0; ii < nnode; ii++)
					{
						fprintf(feig,
						"%4ld %e %e %e %e %e %e\n",
						(nodes + ii)->code,
						*(bisecevct + 6*ii + 0),
						*(bisecevct + 6*ii + 1),
						*(bisecevct + 6*ii + 2),
						*(bisecevct + 6*ii + 3),
						*(bisecevct + 6*ii + 4),
						*(bisecevct + 6*ii + 5));
					}

					if (1)//(abs(evctdot) > 1.0e-3)/*LIMIT POINT*/
					{
						BCLFLAG = -1;
						loadfactor = nextloadfactor;
						for (ii = 0; ii < msize; ii++)
						{
							*(ddisp + ii) = *(nextddisp + ii);
						}
					}
					else/*PATH=SWITCHING FOR BIFURCATION*/
					{
						BCLFLAG = 2;
						loadfactor = lastloadfactor;
						for (ii = 0; ii < msize; ii++)
						{
							*(ddisp + ii) = *(lastddisp + ii);
						}
					}
					nlap++;
					iteration = 0;
				}
			}
			if (pinpointmode == 2)/*PIN-POINTING BASED ON EXTENDED SYSTEM.*/
			{
#if 0
				nline = forwardbackward(gmtx, dup, confs, msize, gcomp1);
				nline = forwardbackward(gmtx, due, confs, msize, gcomp1);

				/*FOR DIRECTIONAL DERIVATIVE OF TANGENTIAL STIFFNESS MATRIX*/
				LDLmode(gmtx, confs, m, ldlneig, ldlevct, norm, dm, msize);
				for (i = 0; i < ldlneig; i++)
				{
					fprintf(fbcl, "CRITICAL EIGEN VECTOR : LINE = %5ld NORM = %12.9f dm = %12.9f\n", *(m + i), *(norm + i), *(dm + i));
					for (ii = 0; ii < nnode; ii++)
					{
						fprintf(fbcl, "%5ld %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f\n", (nodes + ii)->code,
							*(*(ldlevct + i) + (6 * ii + 0)), *(*(ldlevct + i) + (6 * ii + 1)), *(*(ldlevct + i) + (6 * ii + 2)),
							*(*(ldlevct + i) + (6 * ii + 3)), *(*(ldlevct + i) + (6 * ii + 4)), *(*(ldlevct + i) + (6 * ii + 5)));
					}
				}

				for (ii = 0; ii < msize; ii++)/*UPDATE FORM FOR DIRECTIONAL DERIVATIVE*/
				{
					*(epsgvct + ii) = eps * *(*(ldlevct + 0) + ii);
					*(epsddisp + ii) = *(ddisp + ii);
				}
				for (ii = 0; ii < msize; ii++)
				{
					*(epsgvct + ii) = *(epsgvct + *(constraintmain + ii));
				}
				updaterotation(epsddisp, epsgvct, nnode);





				for (ii = 0; ii < msize; ii++)
				{
					*(re + ii) = 0.0;
					*(rp + ii) = 0.0;
				}
				volume = 0.0;
				for (i = 1; i <= nshell; i++)
				{
					inputshell(shells, mshell, i - 1, &shell);
					nnod = shell.nnod;
					shell.sect = (shells + i - 1)->sect;                      /*READ SECTION DATA.*/
					loffset = (int*)malloc(6*nnod * sizeof(int));
					for (ii = 0; ii < nnod; ii++)
					{
						for (jj = 0; jj < 6; jj++)
						{
							*(loffset + (6 * ii + jj)) = 6 * (shell.node[ii]->loff) + jj;
						}
					}

					for (ii = 0; ii < nnod; ii++)
					{
						inputnode(iform, shell.node[ii]);
					}
					drccosinit = shelldrccos(shell, &area);
					gforminit = extractshelldisplacement(shell, iform);
					eforminit = extractlocalcoord(gforminit,drccosinit,nnod);

					Ke = assemshellemtx(shell, drccosinit, NULL); /*ELASTIC MATRIX OF SHELL[K].*/

					for (ii = 0; ii < nnod; ii++)
					{
						inputnode(ddisp, shell.node[ii]);
					}
					drccos = shelldrccos(shell, &area);
					gform = extractshelldisplacement(shell, ddisp);
					eform = extractlocalcoord(gform,drccos,nnod);

					edisp = extractdeformation(eforminit, eform, nnod);                   /*{Ue}*/
					einternal = matrixvector(Ke, edisp, 6 * nnod);

					T = transmatrixIII(drccos, nnod);     /*TRANSFORMATION MATRIX[T].*/

					volume += shellvolume(shell, drccos, area);                    /*VOLUME*/
					epressure = assemshellpvct(shell, drccos);/*ELEMENT EXTERNAL FORCE{Fe}.*/

					ginternal = (double*)malloc(6 * nnod * sizeof(double));
					HPT = (double**)malloc(6 * nnod * sizeof(double*));
					for (ii = 0; ii <6 * nnod; ii++)
					{
						*(HPT + ii) = (double*)malloc(6 * nnod * sizeof(double));
					}
					Kt = assemtmtxCR(Ke, eform, edisp, einternal, ginternal, T, HPT, nnod);


					for (ii = 0; ii < 18; ii++)
					{
						for (jj = 0; jj < 18; jj++)
						{
							*(re + *(loffset + ii)) += *(*(Kt + ii) + jj) * *(due + *(loffset + jj));
							*(rp + *(loffset + ii)) += *(*(Kt + ii) + jj) * *(dup + *(loffset + jj));
						}
					}

					freematrix(drccosinit, 3);
					freematrix(drccos, 3);
					freematrix(T, 6*nnod);
					freematrix(Ke, 6*nnod);
					freematrix(Kt, 6*nnod);
					freematrix(HPT, 6*nnod);

					free(loffset);

					free(einternal);
					free(ginternal);

					free(epressure);
					free(eforminit);
					free(gforminit);
					free(eform);
					free(gform);
					free(edisp);
				}

				for (i = 0; i < msize; i++)/*UPDATE FORM FOR DIRECTIONAL DERIVATIVE*/
				{
					*(re + i) -= *(funbalance + i);
					*(rp + i) -= *(fpressure + i);
					if ((confs + i)->iconf != 0)
					{
						*(re + i) = 0.0;
						*(rp + i) = 0.0;
					}
				}

				/*dupdue=*(dm+0)/(*(norm+0)**(norm+0)); */
				fprintf(fbcl, "EIGENVALUE=%12.9f\n", dupdue);

				/*dupdue*=eps;*/
				dupdue = 0;
				dupdup = 0;
				fprintf(fonl, "dupdue=%12.15f dupdup=%12.15f\n", dupdue, dupdup);
				for (i = 0; i < msize; i++)
				{
					if ((confs + i)->iconf == 0)
					{
						dupdup += *(rp + i) * *(*(ldlevct + 0) + i);
						dupdue += *(re + i) * *(*(ldlevct + 0) + i);
					}
				}

				fprintf(fonl, "dupdue=%12.15f dupdup=%12.15f\n", dupdue, dupdup);
				dupdue += eps * *(dm + 0) / (*(norm + 0) * *(norm + 0));

				lambda = -dupdue / dupdup;
				fprintf(fonl, "LAMBDA=%12.9f\n", lambda);

				/*NEW TARGET*/
				loadfactor += lambda;
				for (ii = 0; ii < msize; ii++)
				{
					*(gvct + ii) = *(dup + ii) * lambda + *(due + ii);/*gvct:{ƒÂU_e+ƒÂƒ©U_p}*/
				}


				/*PINPOINTING END*/
				if (iteration != 1 && fabs(*(dm + 0)) < tolerance)
				{
					nlap++;
					iteration = 0;
					ldlneig = 0;
					freematrix(ldlevct, ldlneig);
					free(norm);
					free(dm);
					free(m);
				}
				/*FORMATION UPDATE FOR PREDICTOR/CORRECTOR*/
				for (ii = 0; ii < msize; ii++)
				{
					*(gvct + ii) = *(gvct + *(constraintmain + ii));
				}
				updaterotation(ddisp, gvct, nnode);
#endif
			}
			iteration++;
		}
		else if (BCLFLAG == 2)
		{
#if 0
			eps = 0.0;
			gamma = 0.0;

			residual = 0.0;
			for (i = 0; i < msize; i++)
			{
				*(dup + i) = *(*(evct + 0) + i);
				*(fswitching + i) = gamma * *(*(evct + 0) + i) + *(funbalance + i);
				if ((confs + i)->iconf == 1)
				{
					*(fswitching + i) = 0.0;
					*(dup + i) = 0.0;
				}
				residual += *(fswitching + i) * *(fswitching + i);
				*(due + i) = *(fswitching + i);
			}


			/*GLOBALLY CONVERGENT NONLINEAR SOLUTION METHOD*/
			nline = croutlu(gmtx, confs, msize, &determinant, &sign, gcomp1);

			sprintf(string, "LAP: %4d ITER: %2d {LOAD}= % 5.8f {RESD}= %1.6e {DET}= %8.5f {SIGN}= %2.0f {BCL}= %1d {EPS}=%1.5e {V}= %8.5f\n",
				nlap, iteration, loadfactor, residual, determinant, sign, BCLFLAG, gamma, volume);
			fprintf(ffig, "%s", string);
			errormessage(string);


			if (iteration == 1)/*PREDICTOR CALCULATION*/
			{
				nline = forwardbackward(gmtx, dup, confs, msize, gcomp1);
				scaledarclength = 1e-3 * arclength;


				/*INCREMENTAL CALCULATION*/
				arcsum = *(weight + msize) * *(weight + msize);
				for (ii = 0; ii < msize; ii++)
				{
					arcsum += *(weight + ii) * *(weight + ii) * *(dup + ii) * *(dup + ii);/*SCALED DISPLACEMENT.*/
				}
				lambda = scaledarclength / sqrt(arcsum);/*INCREMANTAL LOAD FACTOR.*/
				for (ii = 0; ii < msize; ii++)
				{
					if ((confs + ii)->iconf != 1)
					{
						*(gvct + ii) = gamma * *(dup + ii);/*INCREMANTAL DISPLACEMENT.*/
					}
				}
				fprintf(fonl, "LAP: %4d ITER: %2d {LAMBDA}= % 5.8f {TOP}= % 5.8f {BOTTOM}= % 5.8f\n", nlap, iteration, lambda, scaledarclength, sqrt(arcsum));
			}
			else/*CORRECTOR CALCULATION*/
			{
				nline = forwardbackward(gmtx, dup, confs, msize, gcomp1);
				nline = forwardbackward(gmtx, due, confs, msize, gcomp1);
				/*MINIMUM RESIDUAL QUANTITIES METHOD*/
				dupdue = 0.0;
				dupdup = *(weight + msize) * *(weight + msize);
				for (ii = 0; ii < msize; ii++)
				{
					if ((confs + ii)->iconf != 1)
					{
						dupdue += *(weight + ii) * *(weight + ii) * *(dup + ii) * *(due + ii);
						dupdup += *(weight + ii) * *(weight + ii) * *(dup + ii) * *(dup + ii);
					}
				}
				lambda = -dupdue / dupdup;
				for (ii = 0; ii < msize; ii++)
				{
					if ((confs + ii)->iconf != 1)
					{
						*(gvct + ii) = gamma * *(dup + ii) + *(due + ii);
					}
				}
				fprintf(fonl, "LAP: %4d ITER: %2d {LAMBDA}= % 5.8f {TOP}= % 5.8f {BOTTOM}= % 5.8f\n", nlap, iteration, lambda, dupdue, dupdup);
			}

			gamma += lambda;
			for (ii = 0; ii < msize; ii++)
			{
				*(gvct + ii) = *(gvct + *(constraintmain + ii));
			}
			updaterotation(ddisp, gvct, nnode);
			if ((residual < tolerance || iteration > maxiteration) && iteration != 1)
			{
				nlap++;
				iteration = 0;
			}
			iteration++;
#endif
			/*LINE SEARCH*/
			evctfunbalance = dotproduct(bisecevct, funbalance, msize);
			fprintf(fbcl, " %5.8e %5.8e \n", eps, evctfunbalance);


			for (ii = 0; ii < msize; ii++)/*UPDATE FORM FOR DIRECTIONAL DERIVATIVE*/
			{
				*(epsddisp + ii) = *(ddisp + ii);
				*(epsgvct + ii) = gradeps * *(bisecevct + ii);
			}
			for (ii = 0; ii < msize; ii++)
			{
				*(epsgvct + ii) = *(epsgvct + *(constraintmain + ii));
			}
			updaterotation(epsddisp, epsgvct, nnode);

			for (ii = 0; ii < msize; ii++)
			{
				*(epsfinternal + ii) = 0.0;
				*(epsfpressure + ii) = 0.0;
			}


			assemshell(shells, mshell, nshell, constraintmain,
					   NULL, gmtx,
					   iform, epsddisp, epsfinternal, epsfexternal);

			for (ii = 0; ii < msize; ii++)
			{
				*(epsfunbalance + ii) = loadfactor * *(epsfpressure + ii) - *(epsfinternal + ii);
				if ((confs + ii)->iconf == 1) *(epsfunbalance + ii) = 0.0;
			}
			fprintf(fbcl, " %5.8e %5.8e \n", eps + gradeps, epsevctfunbalance);
			epsevctfunbalance = (dotproduct(bisecevct, epsfunbalance, msize) - evctfunbalance) / gradeps;/*GRADIENT*/
			eps -= evctfunbalance / epsevctfunbalance;

			for (ii = 0; ii < msize; ii++)/*UPDATE FORM FOR DIRECTIONAL DERIVATIVE*/
			{
				*(ddisp + ii) = *(lastddisp + ii);
				*(gvct + ii) = eps * *(bisecevct + ii);
			}
			for (ii = 0; ii < msize; ii++)/*FORMATION UPDATE FOR PREDICTOR/CORRECTOR*/
			{
					*(gvct + ii) = *(gvct + *(constraintmain + ii));
			}
			updaterotation(ddisp, gvct, nnode);
			iteration++;

			/*if (abs(evctfunbalance) < 1e-8)
			{
				BCLFLAG = -2;
				ldlneig = 0;
			}*/

		}


		/*TERMINATION FLAG*/
		if (iteration == 1)
		{
			if (fnode != NULL && fnodedrc != NULL)
			{
				for (ii = 0; ii < nnode; ii++)
				{
					if ((nodes + ii)->code == fnode)
					{
						if (*(ddisp + 6 * ii + fnodedrc) - *(iform + 6 * ii + fnodedrc) < fnodemin || *(ddisp + 6 * ii + fnodedrc) - *(iform + 6 * ii + fnodedrc) > fnodemax)
						{
							ENDFLAG = 1;
							sprintf(string,"DEFORMATION EXCEED THRESHOLD(DEFORMATION = %e). ANALYSIS TERMINATED.\n", *(ddisp + 6 * ii + fnodedrc) - *(iform + 6 * ii + fnodedrc));
							errormessage(string);
						}
					}
				}
			}

		}
		/*GO TO NEXT LAP&ITERATION*/
		if (BCLFLAG==-1 && sign > 20)
		{
			ENDFLAG = 1;
			sprintf(string,"DIVERGENCE DITECTED(SIGN = %5ld). ANALYSIS TERMINATED.\n", sign);
			errormessage(string);
		}
		if (0/*nlap>20 && loadfactor < 0.0*/)
		{
			ENDFLAG = 1;
			sprintf(string,"NEGATIVE LOAD DITECTED(LOADFACTOR = %8.5f). ANALYSIS TERMINATED.\n", loadfactor);
			errormessage(string);
		}



		//MESSAGE FOR UPDATE UI
		//AVOID FREEZE ON LONG RUNNING TASK
		MSG msg;
		while (PeekMessage(&msg, NULL, 0, 0, PM_REMOVE))
		{
			TranslateMessage(&msg);
			DispatchMessage(&msg);
		}
#if 1
		while(GetAsyncKeyState(VK_LBUTTON))  /*LEFT CLICK TO CONTINUE.*/
		{
		  if(GetAsyncKeyState(VK_RBUTTON))      /*RIGHT CLICK TO ABORT.*/
		  {
			if(fonl!=NULL) fprintf(fonl,"\nABORTED.\n");

			fclose(fonl);
			fclose(fdsp);
			fclose(fexf);
			fclose(finf);
			fclose(fubf);
			fclose(frct);
			fclose(fstr);
			fclose(fene);
			fclose(ffig);
			fclose(fbcl);
			fclose(feig);


			gfree(gmtx,nnode);
			gfree(epsgmtx,nnode);
			free(weight);
			free(gvct);
			free(funbalance);
			free(freaction);
			free(fexternal);
			free(finternal);
			free(fpressure);
			free(fbaseload);
			free(fgivendisp);
			free(fswitching);
			free(due);
			free(dup);
			free(lastddisp);
			free(lastgvct);
			free(lastpivot);
			free(lapddisp);
			free(nextddisp);
			free(bisecevct);
			free(epsddisp);
			free(epsgvct);
			free(epsfunbalance);
			free(epsfinternal);
			free(epsfexternal);
			free(epsfpressure);
			free(re);
			free(rp);

			laptime("\0",t0);
			return 1;
		  }

		  //t2=clock();
		  //if((t2-t1)/CLK_TCK>=WAIT) break;               /*CONTINUE AFTER WAITING.*/
		}
#endif
	}





	assemelem(elems, melem, nelem, constraintmain, NULL, NULL, iform, ddisp, NULL, NULL);
	assemshell(shells, mshell, nshell, constraintmain, NULL, NULL, iform, ddisp, NULL, NULL);

	////////// OUTPUT RESULT FOR SRCAN//////////
	if(fout!=NULL)
	{
		laptime("OUTPUT INTO FILE.",t0);
		errormessage("STRESS.");
		fprintf(fout,"\n\n");
		fprintf(fout,"** FORCES OF MEMBER\n\n");
		fprintf(fout,"   NO   KT NODE            N           Q1           Q2");
		fprintf(fout,"           MT           M1           M2\n\n");
		fprintf(fout,"\n\n");
        for (i = 0; i < msize; i++)
		{
		  *(gvct+i)=*(ddisp+i)-*(iform+i);
		}
		outputdisp001(gvct,fout,nnode,nodes);
		fprintf(fout,"\n\n");
	}

	clearwindow(*(wdraw.childs+1));
	drawarclmframe((wdraw.childs+1)->hdcC,(wdraw.childs+1)->vparam,*af,0,ONSCREEN);
	overlayhdc(*(wdraw.childs + 1), SRCPAINT);                  /*UPDATE DISPLAY.*/

	fclose(fonl);
	fclose(fdsp);
	fclose(fexf);
	fclose(finf);
	fclose(fubf);
	fclose(frct);
	fclose(fstr);
	fclose(fene);
	fclose(ffig);
	fclose(fbcl);
	fclose(fout);
	fclose(feig);

	gfree(gmtx,nnode);
	gfree(epsgmtx,nnode);
	free(weight);
	free(gvct);
	free(funbalance);
	free(freaction);
	free(fexternal);
	free(finternal);
	free(fpressure);
	free(fbaseload);
	free(fgivendisp);
	free(fswitching);
	free(due);
	free(dup);
	free(lastddisp);
	free(lastgvct);
	free(lastpivot);
	free(lapddisp);
	free(nextddisp);
	free(bisecevct);
	free(epsddisp);
	free(epsgvct);
	free(epsfunbalance);
	free(epsfinternal);
	free(epsfexternal);
	free(epsfpressure);
	free(re);
	free(rp);

	errormessage(" ");
	errormessage("COMPLETED.");

	memory1=availablephysicalmemoryEx("REMAIN:");
	sprintf(string,"CONSUMPTION:%ld[BYTES]",(memory0-memory1));
	errormessage(string);

	return 0;
}
/*arclmCR*/




double dotproduct(double* vct1, double* vct2, int vsize)
{
	int i;
	double dot;
	dot = 0.0;
	for (i = 0; i < vsize; i++) dot += (*(vct1 + i)) * (*(vct2 + i));
	return dot;
}/*dotproduct*/

double vectorlength(double* vct, int vsize)
{
	int i;
	double len;
	len = 0.0;
	for (i = 0; i < vsize; i++) len += (*(vct + i)) * (*(vct + i));
	len = sqrt(len);
    return len;
}/*vectorlength*/

void vectornormalize(double* vct, int vsize)
{
	int i;
	double len;
	len = vectorlength(vct,vsize);
	if(len!=0.0)
	{
		for (i = 0; i < vsize; i++) *(vct + i)/=len;
    }
	return;
}/*vectorlength*/




double* rotationvct(double** rmtx)
{
	//char str[50];
	double c, s;
	double* rvct;
	double theta;
	rvct = (double*)malloc(3 * sizeof(double));
	c = 0.5*(  *(*(rmtx + 0) + 0) + *(*(rmtx + 1) + 1) + *(*(rmtx + 2) + 2) - 1  );/*cos(theta)*/
	s = 0.5*(sqrt(pow( *(*(rmtx + 2) + 1) - *(*(rmtx + 1) + 2) ,2)
				+ pow( *(*(rmtx + 0) + 2) - *(*(rmtx + 2) + 0) ,2)
				+ pow( *(*(rmtx + 1) + 0) - *(*(rmtx + 0) + 1) ,2))); /*|sin(theta)|>=0*/



	theta = atan2(s , c);/*theta>=0*/

	if(s > 1e-15)
	{
		*(rvct + 0) = 0.5 * ((*(*(rmtx + 2) + 1) - *(*(rmtx + 1) + 2)))* theta / s;
		*(rvct + 1) = 0.5 * ((*(*(rmtx + 0) + 2) - *(*(rmtx + 2) + 0)))* theta / s;
		*(rvct + 2) = 0.5 * ((*(*(rmtx + 1) + 0) - *(*(rmtx + 0) + 1)))* theta / s;
	}
	else if(c < 0.0)/*theta=PI*/
	{
		*(rvct + 0) = sqrt( ( *(*(rmtx + 0) + 0) + 1.0 ) / 2.0 );
		if (*(rvct + 0) != 0.0)
		{
			*(rvct + 1) = 0.5 * *(*(rmtx + 0) + 1) / (*(rvct + 0));
			*(rvct + 2) = 0.5 * *(*(rmtx + 0) + 2) / (*(rvct + 0));
		}
		else
		{
			*(rvct + 1) = sqrt( ( *(*(rmtx + 1) + 1) + 1.0 ) / 2.0 );
			if (*(rvct + 1) != 0.0)
			{
				*(rvct + 2) = 0.5 * *(*(rmtx + 1) + 2) / (*(rvct + 1));
			}
			else
			{
				*(rvct + 2) = sqrt( ( *(*(rmtx + 2) + 2) + 1.0 ) / 2.0 );
			}
		}
		*(rvct + 0) *= PI;
		*(rvct + 1) *= PI;
		*(rvct + 2) *= PI;
		//sprintf(str,"%e %e %e\n",s,c,theta);
		//errormessage(str);

	}
	else/*theta=0*/
	{
		*(rvct + 0) = 0.0;
		*(rvct + 1) = 0.0;
		*(rvct + 2) = 0.0;
	}

	return rvct;
}


double** rotationmtx(double* rvct)
{
	int i;
	double** rmtx;
	double* n;
	double theta;

	n = (double*)malloc(3 * sizeof(double));
	rmtx = (double**)malloc(3 * sizeof(double*));
	for (i = 0; i < 3; i++)
	{
		*(rmtx + i) = (double*)malloc(3 * sizeof(double));
	}
	theta = sqrt(*(rvct + 0) * *(rvct + 0) + *(rvct + 1) * *(rvct + 1) + *(rvct + 2) * *(rvct + 2));
	if (theta > 1e-15)
	{
		for (i = 0; i < 3; i++)
		{
			*(n + i) = *(rvct + i) / theta;
		}

		*(*(rmtx + 0) + 0) =  cos(theta)            + (1 - cos(theta)) * *(n + 0) * *(n + 0);
		*(*(rmtx + 0) + 1) = -sin(theta) * *(n + 2) + (1 - cos(theta)) * *(n + 0) * *(n + 1);
		*(*(rmtx + 0) + 2) =  sin(theta) * *(n + 1) + (1 - cos(theta)) * *(n + 0) * *(n + 2);
		*(*(rmtx + 1) + 0) =  sin(theta) * *(n + 2) + (1 - cos(theta)) * *(n + 1) * *(n + 0);
		*(*(rmtx + 1) + 1) =  cos(theta)            + (1 - cos(theta)) * *(n + 1) * *(n + 1);
		*(*(rmtx + 1) + 2) = -sin(theta) * *(n + 0) + (1 - cos(theta)) * *(n + 1) * *(n + 2);
		*(*(rmtx + 2) + 0) = -sin(theta) * *(n + 1) + (1 - cos(theta)) * *(n + 2) * *(n + 0);
		*(*(rmtx + 2) + 1) =  sin(theta) * *(n + 0) + (1 - cos(theta)) * *(n + 2) * *(n + 1);
		*(*(rmtx + 2) + 2) =  cos(theta)            + (1 - cos(theta)) * *(n + 2) * *(n + 2);
	}
	else
	{
		*(*(rmtx + 0) + 0) = 1.0;
		*(*(rmtx + 0) + 1) = 0.0;
		*(*(rmtx + 0) + 2) = 0.0;
		*(*(rmtx + 1) + 0) = 0.0;
		*(*(rmtx + 1) + 1) = 1.0;
		*(*(rmtx + 1) + 2) = 0.0;
		*(*(rmtx + 2) + 0) = 0.0;
		*(*(rmtx + 2) + 1) = 0.0;
		*(*(rmtx + 2) + 2) = 1.0;
	}
	free(n);
	return rmtx;
}

double** spinmtx(double* rvct)
{
	int i;
	double** smtx;

	smtx = (double**)malloc(3 * sizeof(double*));
	for (i = 0; i < 3; i++)
	{
		*(smtx + i) = (double*)malloc(3 * sizeof(double));
	}

	*(*(smtx + 0) + 0) = 0.0;
	*(*(smtx + 0) + 1) = -*(rvct + 2);
	*(*(smtx + 0) + 2) = *(rvct + 1);
	*(*(smtx + 1) + 0) = *(rvct + 2);
	*(*(smtx + 1) + 1) = 0.0;
	*(*(smtx + 1) + 2) = -*(rvct + 0);
	*(*(smtx + 2) + 0) = -*(rvct + 1);
	*(*(smtx + 2) + 1) = *(rvct + 0);
	*(*(smtx + 2) + 2) = 0.0;

	return smtx;
}

double** spinfittermtx(double* eform, int nnod)
/*spin-fitter matrix for beam & shell*/
/*G     :3*6nnod matrix(rotational variation by translational motion)*/
/*input :6nnod vector(variation of nodal displacement in local which inclouds noneffective rotational DOF)*/
/*output:3 vector(variation of local coord psuedo rotation)*/
{
	int i, j;
	double A, len;
	double x1, x2, y1, y2, z1;
	double** G;


	G = (double**)malloc(3 * sizeof(double*));
	for (i = 0; i < 3; i++)
	{
		*(G + i) = (double*)malloc(6*nnod * sizeof(double));
		for (j = 0; j < 6*nnod; j++)
		{
			*(*(G + i) + j) = 0.0;
		}
	}


	if(nnod==2)/*[G1,0,G2,0]*/
	{

		x1 = *(eform + 6) - *(eform + 0);
		y1 = *(eform + 7) - *(eform + 1);
		z1 = *(eform + 8) - *(eform + 2);
		len = sqrt(x1 * x1 + y1 * y1 + z1 * z1);

		/*G1*/
		*(*(G + 1) + 2) =  1 / len;
		*(*(G + 2) + 1) = -1 / len;
		/*G2*/
		*(*(G + 1) + 8) = -1 / len;
		*(*(G + 2) + 7) =  1 / len;
	}
	if(nnod==3)/*[G1,0,G2,0,G3,0]*/
	{

		x1 = *(eform + 6) - *(eform + 0);
		y1 = *(eform + 7) - *(eform + 1);
		x2 = *(eform + 12) - *(eform + 0);
		y2 = *(eform + 13) - *(eform + 1);
		A = 0.5 * (x1 * y2 - x2 * y1);
		len = sqrt(x1 * x1 + y1 * y1);

		/*G1*/
		*(*(G + 0) + 2) = 0.5 * (*(eform + 12) - *(eform + 6)) / A;
		*(*(G + 1) + 2) = 0.5 * (*(eform + 13) - *(eform + 7)) / A;
		*(*(G + 2) + 1) = -1 / len;
		/*G2*/
		*(*(G + 0) + 8) = 0.5 * (*(eform + 0) - *(eform + 12)) / A;
		*(*(G + 1) + 8) = 0.5 * (*(eform + 1) - *(eform + 13)) / A;
		*(*(G + 2) + 7) = 1 / len;
		/*G3*/
		*(*(G + 0) + 14) = 0.5 * (*(eform + 6) - *(eform + 0)) / A;
		*(*(G + 1) + 14) = 0.5 * (*(eform + 7) - *(eform + 1)) / A;
	}
	return G;
}

double** projectionmtx(double* eform, double** G,int nnod)
/*element projection matrix for beam & shell*/
/*G     :6nnod*6nnod matrix(variation correction)*/
/*input :6nnod vector(initial variation of displacement in local coord)*/
/*output:6nnod vector(corrected variation of displacement in local coord)*/
{
	int i, j, a, b;
	double* node;
	double** P, ** Gu, ** S, ** SGu;

	P = (double**)malloc(6*nnod * sizeof(double*));
	for (i = 0; i < 6*nnod; i++)
	{
		*(P + i) = (double*)malloc(6*nnod * sizeof(double));
	}

	/*eform : latest coordination of nodes in local*/
	node = (double*)malloc(3 * sizeof(double));

	/*for extracting Gu_b(b=1,2,3) from G*/
	Gu = (double**)malloc(3 * sizeof(double*));
	for (i = 0; i < 3; i++)
	{
		*(Gu + i) = (double*)malloc(3 * sizeof(double));
	}
	//Pab = [ delta*I-1/3*I+S_a*Gu_b    0       ]
	//      [         -Gu_b             delta*I ]

	for (a = 0; a < nnod; a++)/*row 6*a+0~6*a+5 of P*/
	{
		for (i = 0; i < 3; i++)
		{
			*(node + i) = *(eform + 6 * a + i);
		}
		S = spinmtx(node);
		for (b = 0; b < nnod; b++)/*column 6*b+0~6*b+5 of P*/
		{
			for (i = 0; i < 3; i++)
			{
				for (j = 0; j < 3; j++)
				{
					*(*(Gu + i) + j) = *(*(G + i) + 6 * b + j);
				}
			}
			SGu = matrixmatrix(S, Gu, 3);
			for (i = 0; i < 3; i++)
			{
				for (j = 0; j < 3; j++)
				{
					*(*(P + 6 * a + 3 + i) + 6 * b + j) = -*(*(Gu + i) + j);/*-Gu_b*/
					*(*(P + 6 * a + i) + 6 * b + 3 + j) = 0.0;				/*0*/
					if (a == b && i == j)
					{
						*(*(P + 6 * a + i) + 6 * b + j) = 1.0 - 1.0/(double)nnod;	/*delta*I-1/3*I*/
						*(*(P + 6 * a + 3 + i) + 6 * b + 3 + j) = 1.0;		/*delta*I*/
					}
					else if (i == j)
					{
						*(*(P + 6 * a + i) + 6 * b + j) = - 1.0/(double)nnod;		/*delta*I-1/3*I*/
						*(*(P + 6 * a + 3 + i) + 6 * b + 3 + j) = 0.0;		/*delta*I*/
					}
					else
					{
						*(*(P + 6 * a + i) + 6 * b + j) = 0.0;				/*delta*I-1/3*I*/
						*(*(P + 6 * a + 3 + i) + 6 * b + 3 + j) = 0.0;		/*delta*I*/
					}
					*(*(P + 6 * a + i) + 6 * b + j) += *(*(SGu + i) + j);	/*S_a*Gu_b*/
				}
			}
			freematrix(SGu, 3);
		}
		freematrix(S, 3);
	}
	free(node);
	freematrix(Gu, 3);
	return P;
}

double** jacobimtx(double* rvct)
/*transformation from additional infinitesimal incremental rotation to increment of total rotational pseudo-vector*/
{
	int n, i, j;
	double theta, eta;
	double** thetaspin, ** thetaspin2;
	double** H;


	H = (double**)malloc(3 * sizeof(double*));
	for (i = 0; i < 3; i++)
	{
		*(H + i) = (double*)malloc(3 * sizeof(double));
		for (j = 0; j < 3; j++)
		{
			*(*(H + i) + j) = 0.0;
		}
	}

	theta = vectorlength(rvct,3);

	if (theta < PI / 30.0)
	{
		eta = 1.0 / 12.0 + pow(theta, 2) / 720.0 + pow(theta, 4) / 30240.0 + pow(theta, 6) / 1209600.0;
	}
	else
	{
		eta = (1.0 - 0.5 * theta * (1.0 / tan(0.5 * theta))) / pow(theta, 2);
	}
	thetaspin = spinmtx(rvct);
	thetaspin2 = matrixmatrix(thetaspin, thetaspin, 3);


	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < 3; j++)
		{
			*(*(H + i) + j) = -0.5 * *(*(thetaspin + i) + j) + eta * *(*(thetaspin2 + i) + j);
			if (i == j)*(*(H + i) + j) += 1.0;
		}
	}
	freematrix(thetaspin, 3);
	freematrix(thetaspin2, 3);
	return H;
}

double** blockjacobimtx(double* edisp, double* estress, double** M, int nnod)
{
	int n, i, j;
	double theta, dot, eta, mu;
	double* rvct, * mvct;
	double** thetaspin, ** thetaspin2, ** mspin, ** mtheta, ** mtheta2;
	double** H, ** Ha, ** Ma;

	/*H=diag[I, H_1, I, H_2(, I, H_3)]*/
	H = (double**)malloc(6*nnod * sizeof(double*));
	for (i = 0; i < 6*nnod; i++)
	{
		*(H + i) = (double*)malloc(6*nnod * sizeof(double));
		for (j = 0; j < 6*nnod; j++)
		{
			*(*(H + i) + j) = 0.0;
		}
	}
	rvct = (double*)malloc(3 * sizeof(double));

	if(estress!=NULL && M!=NULL)
	{
		mtheta = (double**)malloc(3 * sizeof(double*));                         /*{ma}{ƒÆat}*/
		for (i = 0; i < 3; i++)
		{
			*(mtheta + i) = (double*)malloc(3 * sizeof(double));
		}
		mvct = (double*)malloc(3 * sizeof(double));
	}

	for (n = 0; n < nnod; n++)
	{
		for (i = 0; i < 3; i++)
		{
			*(rvct + i) = *(edisp + 6 * n + 3 + i);
		}
		theta = vectorlength(rvct,3);
        thetaspin = spinmtx(rvct);
		thetaspin2 = matrixmatrix(thetaspin, thetaspin, 3);

		if (theta < PI / 30.0)
		{
			eta = 1.0 / 12.0 + pow(theta, 2) / 720.0 + pow(theta, 4) / 30240.0 + pow(theta, 6) / 1209600.0;
			mu = 1.0 / 360.0 + pow(theta, 2) / 7560.0 + pow(theta, 4) / 201600.0 + pow(theta, 6) / 5987520.0;
		}
		else
		{
			eta = (1.0 - 0.5 * theta / tan(0.5 * theta)) / pow(theta, 2);
			mu = (theta * theta + 4.0 * cos(theta) + theta * sin(theta) - 4.0)
				/ (4.0 * pow(theta, 4) * sin(0.5 * theta) * sin(0.5 * theta));
		}

		Ha = (double**)malloc(3 * sizeof(double*));
		for (i = 0; i < 3; i++)
		{
			*(Ha + i) = (double*)malloc(3 * sizeof(double));
		}
		for (i = 0; i < 3; i++)
		{
			for (j = 0; j < 3; j++)
			{
				*(*(Ha + i) + j) = -0.5 * *(*(thetaspin + i) + j) + eta * *(*(thetaspin2 + i) + j);
				if (i == j)*(*(Ha + i) + j) += 1.0;
			}
		}
		for (i = 0; i < 3; i++)
		{
			for (j = 0; j < 3; j++)
			{
				*(*(H + 6 * n + 3 + i) + 6 * n + 3 + j) = *(*(Ha + i) + j);
				if (i == j)*(*(H + 6 * n + i) + 6 * n + j) = 1.0;
			}
		}


		if(estress!=NULL && M!=NULL)
		{
			for (i = 0; i < 3; i++)
			{
				*(mvct + i) = *(estress + 6 * n + 3 + i);
			}
			dot = dotproduct(rvct,mvct,3);
			for (i = 0; i < 3; i++)
			{
				for (j = 0; j < 3; j++)
				{
					*(*(mtheta + i) + j) = *(mvct + i) * *(rvct + j);
				}
			}
			mspin = spinmtx(mvct);
			mtheta2 = matrixmatrix(thetaspin2, mtheta, 3);
			for (i = 0; i < 3; i++)
			{
				for (j = 0; j < 3; j++)
				{
					*(*(mtheta2 + i) + j) *= mu;
					*(*(mtheta2 + i) + j) += eta * ( *(*(mtheta + j) + i) - 2.0 * *(*(mtheta + i) + j) ) - 0.5 * *(*(mspin + i) + j);
					if (i == j)*(*(mtheta2 + i) + j) += eta * dot;
				}
			}
			Ma = matrixmatrix(mtheta2, Ha, 3);
			for (i = 0; i < 3; i++)
			{
				for (j = 0; j < 3; j++)
				{
					*(*(M + 6 * n + 3 + i) + 6 * n + 3 + j) = *(*(Ma + i) + j);
				}
			}

			freematrix(mspin, 3);
			freematrix(mtheta2, 3);
			freematrix(Ma, 3);
		}

		freematrix(thetaspin, 3);
		freematrix(thetaspin2, 3);
		freematrix(Ha, 3);
	}

	free(rvct);
	if(estress!=NULL && M!=NULL)
	{
		free(mvct);
		freematrix(mtheta, 3);
	}
	return H;
}

double** transmatrixHPT(double* eform, double* edisp, double** T, int nnod)
{
	double** G, ** P, ** H;
	double** HP, ** HPT;

	G = spinfittermtx(eform, nnod);                     						/*SPIN-FITTER MATRIX[G].*/
	P = projectionmtx(eform, G, nnod);    										/*PROJECTION MATRIX[P].*/
	H = blockjacobimtx(edisp, NULL, NULL, nnod);								/*JACOBIAN MATRIX OF ROTATION[H].*/

	HP = matrixmatrix(H, P, 6*nnod);
	HPT = matrixmatrix(HP, T, 6*nnod);

	freematrix(G, 3);
	freematrix(P, 6*nnod);
	freematrix(H, 6*nnod);
	freematrix(HP, 6*nnod);

	return HPT;
}


double** assemgmtxCR(double* eform, double* edisp, double* estress, double* gstress, double** T, double** HPT, int nnod)
{
	int i, j, n;
	double** G, ** P, ** H;
	double** HP, ** PtHt;
	double** TtPtHt;
	double* pstress;
	double* nm, ** spinnm;
	double** Fnm, ** Fn, ** FnG, ** GtFnt;
	double** Kg, ** Kgr, ** Kgp, ** Kgm;


	Kg = (double**)malloc(6*nnod * sizeof(double*));                            /*[Kg]*/
	for (i = 0; i < 6*nnod; i++)
	{
		*(Kg + i) = (double*)malloc(6*nnod * sizeof(double));
	}

	Kgm = (double**)malloc(6*nnod * sizeof(double*));                           /*[M]&[Kgm]*/
	for (i = 0; i < 6*nnod; i++)
	{
		*(Kgm + i) = (double*)malloc(6*nnod * sizeof(double));
		for (j = 0; j < 6*nnod; j++)
		{
			*(*(Kgm + i) + j) = 0.0;
		}
	}

	G = spinfittermtx(eform, nnod);                     						/*SPIN-FITTER MATRIX[G].*/
	P = projectionmtx(eform, G, nnod);    										/*PROJECTION MATRIX[P].*/
	H = blockjacobimtx(edisp, estress, Kgm, nnod);								/*JACOBIAN MATRIX OF ROTATION[H].*/

	HP = matrixmatrix(H, P, 6*nnod);                    						/*[H][P]*/
	PtHt = matrixtranspose(HP, 6*nnod);
	pstress = matrixvector(PtHt, estress, 6*nnod);     							/*projected estress {Fp}*/

	if(HPT!=NULL)
	{
		matrixmatrixII(HPT, HP, T, 6*nnod);                 						/*[H][P][T]*/
		if(gstress!=NULL)
		{
		TtPtHt = matrixtranspose(HPT, 6*nnod);
		matrixvectorII(gstress, TtPtHt, estress, 6*nnod);      						/*global estress {Fg}*/
		freematrix(TtPtHt, 6*nnod);
		}
	}


	nm = (double*)malloc(3 * sizeof(double));           						/*projected estress of each node {n}&{m}*/

	Fnm = (double**)malloc(6*nnod * sizeof(double*));                          	/*[Fnm]*/
	for (i = 0; i < 6*nnod; i++)
	{
		*(Fnm + i) = (double*)malloc(3 * sizeof(double));
	}
	Fn = (double**)malloc(6*nnod * sizeof(double*));                            /*[Fn]*/
	for (i = 0; i < 6*nnod; i++)
	{
		*(Fn + i) = (double*)malloc(3 * sizeof(double));
	}

	for (n = 0; n < nnod; n++)
	{
		for (i = 0; i < 3; i++)
		{
			*(nm + i) = *(pstress + 6 * n + i);
		}
		spinnm = spinmtx(nm);
		for (i = 0; i < 3; i++)
		{
			for (j = 0; j < 3; j++)
			{
				*(*(Fnm + 6 * n + i) + j) = *(*(spinnm + i) + j);
				*(*(Fn + 6 * n + i) + j) = *(*(spinnm + i) + j);
			}
		}
		freematrix(spinnm, 3);

		for (i = 0; i < 3; i++)
		{
			*(nm + i) = *(pstress + 6 * n + 3 + i);
		}
		spinnm = spinmtx(nm);
		for (i = 0; i < 3; i++)
		{
			for (j = 0; j < 3; j++)
			{
				*(*(Fnm + 6 * n + 3 + i) + j) = *(*(spinnm + i) + j);
				*(*(Fn + 6 * n + 3 + i) + j) = 0.0;
			}
		}
		freematrix(spinnm, 3);
	}

	Kgr = matrixmatrixIII(Fnm, G, 6*nnod, 3, 6*nnod);/*[Fnm][G]*/

	FnG = matrixmatrixIII(Fn, G, 6*nnod, 3, 6*nnod);/*[Fn][G]*/
	GtFnt = matrixtranspose(FnG, 6*nnod);/*[Gt][Fnt]*/
	Kgp = matrixmatrix(GtFnt, P, 6*nnod);/*[Gt][Fnt][P]*/

	Kgm = transformationIII(Kgm, P, 6*nnod);/*[Pt][M][P]*/

	for (i = 0; i < 6*nnod; i++)
	{
		for (j = 0; j < 6*nnod; j++)
		{
			*(*(Kg + i) + j) = - *(*(Kgr + i) + j) - *(*(Kgp + i) + j) + *(*(Kgm + i) + j);
		}
	}

	Kg = transformationIII(Kg, T, 6*nnod);

	free(pstress);
	free(nm);

	freematrix(G, 3);
	freematrix(P, 6*nnod);
	freematrix(H, 6*nnod);
	freematrix(HP, 6*nnod);
	freematrix(PtHt, 6*nnod);
	freematrix(Fnm, 6*nnod);
	freematrix(Fn, 6*nnod);

	freematrix(Kgr, 6*nnod);
	freematrix(FnG, 6*nnod);
	freematrix(GtFnt, 6*nnod);
	freematrix(Kgp, 6*nnod);
	freematrix(Kgm, 6*nnod);
	return Kg;
}


double** assemtmtxCR(double** Ke, double* eform, double* edisp, double* estress, double* gstress, double** T, double** HPT, int nnod)
{
	int i, j;
	double** Kt;

	Kt = assemgmtxCR(eform, edisp, estress, gstress, T, HPT, nnod);/*[Kg]=[Kgr]+[Kgp]+[Kgm]*/
	Ke = transformationIII(Ke, HPT, 6*nnod);/*[Ke]=[Tt][Pt][Ht][K][H][P][T]*/

	for (i = 0; i < 6*nnod; i++)
	{
		for (j = 0; j < 6*nnod; j++)
		{
			*(*(Kt + i) + j) += *(*(Ke + i) + j);/*[Kt]=[Ke]+[Kg]*/
		}
	}
	return Kt;
}


void symmetricmtx(double** estiff, int msize)
{
	int i, j;
	for (i = 0; i < msize; i++)
	{
		for (j = 0; j < i; j++)
		{
			*(*(estiff + i) + j) = 0.5 * (*(*(estiff + i) + j) + *(*(estiff + j) + i));
			*(*(estiff + j) + i) = *(*(estiff + i) + j);
		}
	}
	return;
}


void assemelem(struct owire* elems, struct memoryelem* melem, int nelem, long int* constraintmain,
			   struct gcomponent* mmtx, struct gcomponent* gmtx,
			   double* iform, double* ddisp, double* finternal, double* fexternal)
{
	struct owire elem;
	int i,j,ii,jj;
	int nnod;
	double* gforminit, * eforminit;                      /*DEFORMED COORDINATION OF ELEMENT*/
	double* gform, * eform;                      /*INITIAL COORDINATION OF ELEMENT*/
	double* edisp;
	double* ginternal, * einternal;                        /*INTERNAL FORCE OF ELEMENT*/
	//double* gexternal, * eexternal;                          /*EXTERNAL FORCE OF ELEMENT*/
	double** Me,** Ke,** Kt,** DBe,** B,** drccos,** drccosinit;                           /*MATRIX*/
	double** T,** Tt, **HPT,**TtPtHt;
	int* loffset;
	double area;


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
		Ke = modifyhinge(elem,Ke);             /*MODIFY MATRIX.*/

		/*DEFORMED CONFIDURATION*/
		for (ii = 0; ii < nnod; ii++)
		{
			inputnode(ddisp, elem.node[ii]);
		}
		gform = extractdisplacement(elem, ddisp);
		drccos = updatedrccos(drccosinit, gforminit, gform);
		eform = extractlocalcoord(gform,drccos,nnod);

		T = transmatrixIII(drccos, nnod);									/*[T].*/
		Tt = matrixtranspose(T, 6 * nnod);                    				/*[Tt].*/

		edisp = extractdeformation(eforminit, eform, nnod);                	/*{Ue}*/
		einternal = matrixvector(Ke, edisp, 6 * nnod);          			/*{Fe}=[Ke]{Ue}.*/

		ginternal = (double*)malloc(6 * nnod * sizeof(double));
		HPT = (double**)malloc(6 * nnod * sizeof(double*));
		for (ii = 0; ii <6 * nnod; ii++)
		{
			*(HPT + ii) = (double*)malloc(6 * nnod * sizeof(double));
		}
		Kt = assemtmtxCR(Ke, eform, edisp, einternal, ginternal, T, HPT, nnod);	/*TANGENTIAL MATRIX[Kt].*/
		symmetricmtx(Kt, 6*nnod);											/*SYMMETRIC TANGENTIAL MATRIX[Ksym].*/

		if(gmtx!=NULL)assemgstiffnesswithDOFelimination(gmtx, Kt, &elem, constraintmain); /*ASSEMBLAGE TANGENTIAL STIFFNESS MATRIX.*/



		for (ii = 0; ii < nnod; ii++)
		{
			for (jj = 0; jj < 6; jj++)
			{
				if(finternal!=NULL)
				{
					*(finternal + *(constraintmain + *(loffset + (6 * ii + jj)))) += *(ginternal + 6 * ii + jj);
				}
				(elems+i-1)->stress[ii][jj]=*(einternal+6*ii+jj);
				if(melem!=NULL)
				{
					(melem+i-1)->stress[ii][jj]=*(einternal+6*ii+jj);
				}
			}
		}




		freematrix(drccosinit, 3);
		freematrix(drccos, 3);
		freematrix(T, 6 * nnod);
		freematrix(Tt, 6 * nnod);
		freematrix(HPT, 6 * nnod);
		freematrix(Ke, 6 * nnod);
		freematrix(Kt, 6 * nnod);

		free(einternal);
		free(ginternal);

		free(eforminit);
		free(gforminit);
		free(eform);
		free(gform);
		free(edisp);
		free(loffset);
	}
	return;
}

void assemshell(struct oshell* shells, struct memoryshell* mshell, int nshell, long int* constraintmain,
			   struct gcomponent* mmtx, struct gcomponent* gmtx,
			   double* iform, double* ddisp, double* finternal, double* fexternal)
{
	struct oshell shell;
	int i,j,ii,jj;
	int nnod;
	double* gforminit, * eforminit;                      /*DEFORMED COORDINATION OF ELEMENT*/
	double* gform, * eform;                      /*INITIAL COORDINATION OF ELEMENT*/
	double* edisp;
	double* ginternal, * einternal;                        /*INTERNAL FORCE OF ELEMENT*/
	double* gexternal, * eexternal;                          /*EXTERNAL FORCE OF ELEMENT*/
	double** Me,** Ke,** Kt,** DBe,** B,** drccos,** drccosinit;                           /*MATRIX*/
	double** T,** Tt,/*** HP,**PtHt,*/**HPT,**TtPtHt;
	int* loffset;
	double area;


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
		drccosinit = shelldrccos(shell, &area);
		gforminit = extractshelldisplacement(shell, iform);                 /*{Xg}*/
		eforminit = extractlocalcoord(gforminit,drccosinit,nnod);        	/*{Xe}*/

		B = (double***)malloc(7 * sizeof(double**));/*FOR 3-NODES ELEMENT*/
		for (ii = 0; ii < 7; ii++)
		{
			*(B + ii) = (double**)malloc(6 * sizeof(double*));
			for (jj = 0; jj < 6; jj++)
			{
				*(*(B + ii) + jj) = (double*)malloc(6*nnod * sizeof(double));
				for (kk = 0; kk < 6*nnod; kk++)
				{
					*(*(*(B + ii) + jj) + kk) = 0.0;
				}
			}
		}
		B = shellshape(shell, drccosinit);

		Ke = assemshellemtx(shell, drccosinit, B);                        /*[Ke]*/

		Me = assemshellmmtx(shell, drccosinit);         				  /*[Me]*/



		/*DEFORMED CONFIGURATION OF LAST LAP.*/
		for (ii = 0; ii < nnod; ii++)
		{
			inputnode(lastddisp, shell.node[ii]);
		}
		lastdrccos = shelldrccos(shell, &area);
		lastgform = extractshelldisplacement(shell, lastddisp);             /*{Xg+Ug}*/
		lasteform = extractlocalcoord(lastgform,lastdrccos,nnod);           /*{Xe+Ue}*/

		lastT = transmatrixIII(lastdrccos, nnod);         					/*[T]*/
		lastTt = matrixtranspose(lastT, 6 * nnod);                  		/*[Tt]*/

		lastedisp = extractdeformation(eforminit, lasteform, nnod);         /*{Ue}*/

		/*DEFORMED CONFIGFURATION*/
		for (ii = 0; ii < nnod; ii++)
		{
			inputnode(ddisp, shell.node[ii]);
		}
		drccos = shelldrccos(shell, &area);
		gform = extractshelldisplacement(shell, ddisp);                     /*{Xg+Ug}*/
		eform = extractlocalcoord(gform,drccos,nnod); 			       	    /*{Xe+Ue}*/

		T = transmatrixIII(drccos, nnod);         							/*[T].*/
		Tt = matrixtranspose(T, 6 * nnod);                  				/*[Tt].*/

		edisp = extractdeformation(eforminit, eform, nnod);           		/*{Ue}*/





		for(ii = 0; ii < 6*nnod; ii++)
		{
		  *(edisp+ii) -= *(lastedisp+ii);
		}

		for (j = 0; j < 7; j++)/*FOR EACH INTEGRATION POINT*/
		{

		  estrain = matrixmatrixIII(*(B+j),edisp,6,6*nnod,1); /*INCREMENTAL STRAIN d{ƒÃx ƒÃy ƒÃxy ƒÁx ƒÁy ƒÁxy} FOR NODE-ii*/
		  estress = matrixvector(C,estrain,6);/*ELASTIC PREDICTOR*/

		  for(ii = 0; ii < 6; ii++)
		  {
			*(estress+ii) += (mshell+i-1)->laststress[j][ii];/*TRIAL STRESS{ƒÐtry}*/

			//*(backstress+ii) = (mshell+i-1)->backstress[j][ii];
			//*(estress+ii) -= *(backstress+ii);
		  }

		  lambda = returnmapIlyushin(shell, estress);

		  /*STRESS RESULTANT OF INTEGRATION POINTS*/
		  for (ii = 0; ii < 6; ii++)
		  {
			(mshell+i-1)->stress[j][ii] = *(estress+ii);
		  }
		  free(lambda);

		  Bt=matrixtransposeIII(B,6,6*nnod);
		  einternal=matrixmatrixIII(Bt,estress,6*nnod,6,1);
		}
		/*STRESS RESULTANT UPDATED.*/
		/*INTERNAL FORCE UPDATED*/



		//einternal = matrixvector(Ke, edisp, 6 * nnod);      				/*{Fe}=[Ke]{Ue}.*/
		eexternal= assemshellpvct(shell, drccos);                			/*{Pe}.*/
		//volume += shellvolume(shell, drccos, area);                   		/*VOLUME*/

		ginternal = (double*)malloc(6 * nnod * sizeof(double));
		HPT = (double**)malloc(6 * nnod * sizeof(double*));
		for (ii = 0; ii < 6 * nnod; ii++)
		{
			*(HPT + ii) = (double*)malloc(6 * nnod * sizeof(double));
		}
		Kt = assemtmtxCR(Ke, eform, edisp, einternal, ginternal, T, HPT, nnod);	/*TANGENTIAL MATRIX[Kt].*/
		symmetricmtx(Kt, 6*nnod);											/*SYMMETRIC TANGENTIAL MATRIX[Ksym].*/
		if(gmtx!=NULL)assemgstiffnessIIwithDOFelimination(gmtx, Kt, &shell, constraintmain); /*ASSEMBLAGE TANGENTIAL STIFFNESS MATRIX.*/

		Me = transformationIII(Me, T, 6*nnod);/*[Ke]=[Tt][Pt][Ht][K][H][P][T]*/
		if(mmtx!=NULL)assemgstiffnessIIwithDOFelimination(mmtx, Me, &shell, constraintmain); /*ASSEMBLAGE TANGENTIAL STIFFNESS MATRIX.*/

		gexternal = matrixvector(Tt, eexternal, 6 * nnod);  /*GLOBAL EXTERNAL FORCE{Pg}.*/




		for (ii = 0; ii < nnod; ii++)
		{
			for (jj = 0; jj < 6; jj++)
			{
				if(finternal!=NULL) *(finternal + *(constraintmain + *(loffset + (6 * ii + jj)))) += *(ginternal + 6 * ii + jj);
				if(fexternal!=NULL) *(fexternal + *(constraintmain + *(loffset + (6 * ii + jj)))) += *(gexternal + 6 * ii + jj);
				(shells+i-1)->stress[ii][jj]=*(einternal+6*ii+jj);
			}
		}



		if(mshell!=NULL)
		{
			(mshell+i-1)->SE  = 0.0;
			(mshell+i-1)->SEp = 0.0;
			(mshell+i-1)->SEb = 0.0;
			for (ii = 0; ii < nnod; ii++)
			{
				for (jj = 0; jj < 2; jj++)
				{
					(mshell+i-1)->SEp += 0.5 * *(edisp + 6 * ii + jj) * *(einternal + 6 * ii + jj);
				}
				for (jj = 2; jj < 5; jj++)
				{
					(mshell+i-1)->SEb += 0.5 * *(edisp + 6 * ii + jj) * *(einternal + 6 * ii + jj);
				}
				(mshell+i-1)->SE += 0.5 * *(edisp + 6 * ii + 5) * *(einternal + 6 * ii + 5);
			}
			(mshell+i-1)->SE += (mshell+i-1)->SEp + (mshell+i-1)->SEb;
		}

		freematrix(drccos, 3);
		freematrix(drccosinit, 3);
		freematrix(T, 6 * nnod);
		freematrix(Tt, 6 * nnod);
		freematrix(HPT, 6 * nnod);
		freematrix(Ke, 6 * nnod);
		freematrix(Kt, 6 * nnod);

		freematrix(B, 6 * (2*nnod+1));
		freematrix(Me, 6 * nnod);

		free(einternal);
		free(ginternal);

		free(eexternal);
		free(gexternal);

		free(eforminit);
		free(gforminit);
		free(eform);
		free(gform);
		free(edisp);

		free(loffset);
	}
	return;
}


void updaterotation(double* ddisp, double* gvct, int nnode)
/*FORMATION UPDATE IF ROTATION IS FINITE.*/
{
	int i, j;
	long int loff;
	double* rvctR, * rvctL, * rvct;
	double** rmtxR, ** rmtxL, ** rmtx;
	rvctR = (double*)malloc(3 * sizeof(double));
	rvctL = (double*)malloc(3 * sizeof(double));
	for (i = 0; i < nnode; i++)
	{
		for (j = 0; j < 3; j++)
		{
			loff = 6 * i + j;
			*(ddisp + loff) += *(gvct + loff);
		}
		for (j = 0; j < 3; j++)
		{
			loff = 6 * i + 3 + j;
			*(rvctR + j) = *(ddisp + loff);
			*(rvctL + j) = *(gvct + loff);
		}
		rmtxR = rotationmtx(rvctR);
		rmtxL = rotationmtx(rvctL);
		rmtx = matrixmatrix(rmtxL, rmtxR, 3);
		rvct = rotationvct(rmtx);
		for (j = 0; j < 3; j++)
		{
			loff = 6 * i + 3 + j;
			*(ddisp + loff) = *(rvct + j);
		}
		freematrix(rmtxR, 3);
		freematrix(rmtxL, 3);
		freematrix(rmtx, 3);
		free(rvct);
	}
	free(rvctR);
	free(rvctL);
	return;
}/*updaterotation*/



double* quaternionvct(double* rvct)
{
	int i;
	double* qvct;
	double theta;

	qvct = (double*)malloc(4 * sizeof(double));
	theta = sqrt(*(rvct + 0) * *(rvct + 0) + *(rvct + 1) * *(rvct + 1) + *(rvct + 2) * *(rvct + 2));
	if(theta!=0)
	{
		for (i = 0; i < 3; i++)
		{
			*(qvct + i) = sin(0.5*theta)**(rvct + i)/theta;
		}
	}
	else
	{
        for (i = 0; i < 3; i++)
		{
			*(qvct + i) = 0.0;
		}
    }
	*(qvct + 3) = cos(0.5*theta);
	return qvct;
}


double** updatedrccos(double** drccosinit, double* gforminit, double* gform)
{
	int i;
	double** drccos;
	double* rvct1, * rvct2;
	double** rmtxinit, ** rmtxinitt, **rmtx, ** drccosrmtx;

	rvct1 = (double*)malloc(3 * sizeof(double));
	rvct2 = (double*)malloc(3 * sizeof(double));

	for (i = 0; i < 3; i++)
	{
		*(rvct1 + i) = *(gforminit + 3 + i);
		*(rvct2 + i) = *(gforminit + 9 + i);
	}
	rmtxinit = interpolatermtx(rvct1, rvct2, 0.5);
	rmtxinitt = matrixtranspose(rmtxinit, 3);

	for (i = 0; i < 3; i++)
	{
		*(rvct1 + i) = *(gform + 3 + i);
		*(rvct2 + i) = *(gform + 9 + i);
	}
	rmtx = interpolatermtx(rvct1, rvct2, 0.5);

	drccosrmtx = matrixmatrix(rmtx, rmtxinitt, 3);
	drccos = matrixmatrix(drccosrmtx, drccosinit, 3);

	free(rvct1);
	free(rvct2);

	freematrix(rmtxinit, 3);
	freematrix(rmtxinitt, 3);
	freematrix(rmtx, 3);
	freematrix(drccosrmtx, 3);

	return drccos;
}


double** interpolatermtx(double* rvct1, double* rvct2, double alpha)
{
	int i;
	double* rvct;
	double** rmtx1, ** rmtx2, **trmtx1, ** rmtx;
	double** alpharmtx, ** midrmtx;


	/*
	qvct1 = quaternionvct(rvct1);
	qvct2 = quaternionvct(rvct2);
	dot = dotproduct(qvct1,qvct2,4);
	theta = acos(dot);

	if(theta!= 0)
	{
		for (i = 0; i < 4; i++)
		{
			*(midqvct + i) = (*(qvct1 + i)+*(qvct2 + i))/(2*cos(0.5*theta));
		}
	}
	*/



	rmtx1 = rotationmtx(rvct1);
	rmtx2 = rotationmtx(rvct2);
	trmtx1 = matrixtranspose(rmtx1, 3);
	rmtx = matrixmatrix(rmtx2, trmtx1, 3);
	rvct = rotationvct(rmtx);

	for (i = 0; i < 3; i++)
	{
		*(rvct+i)*=alpha;
	}

	alpharmtx = rotationmtx(rvct);
	midrmtx = matrixmatrix(alpharmtx, rmtx1, 3);


	freematrix(rmtx1, 3);
	freematrix(rmtx2, 3);
	freematrix(trmtx1, 3);
	freematrix(rmtx, 3);
	free(rvct);
	freematrix(alpharmtx, 3);

	return midrmtx;
}



double* extractlocalcoord(double* gform, double** drccos, double nnod)
/*EXTRACT LOCAL ELEMENT DEFORMATION FROM GLOBAL.*/
/*UPDATE PSUEDO-ROTATION VECTOR*/
{
	long int i, n;
	double* d, * r, * c, * td, * tr, * eform;
	double** trmtx, ** rmtx;

	eform = (double*)malloc(6 * nnod * sizeof(double));

	c = (double*)malloc(3 * sizeof(double));
	d = (double*)malloc(3 * sizeof(double));
	r = (double*)malloc(3 * sizeof(double));

	for (i = 0; i < 3; i++)
	{
		*(c + i) = 0.0;
		for (n = 0; n < nnod; n++)
		{
			*(c + i) += *(gform + 6 * n + i) / nnod;
		}
	}
	/*CENTER*/
	for (n = 0; n < nnod; n++)
	{
		for (i = 0; i < 3; i++)
		{
			*(d + i) = *(gform + 6 * n + i) - *(c + i);
			*(r + i) = *(gform + 6 * n + 3 + i);
		}

		td = matrixvector(drccos, d, 3);
		/*EACH NODE FROM CENTER*/

		rmtx = rotationmtx(r);

		/*rmtx:Ra*/
		trmtx = matrixmatrix(drccos, rmtx, 3);/*TRIAD DIRECTION MATRIX(3 VECTOR) IN LOCAL*/
		tr = rotationvct(trmtx);

		for (i = 0; i < 3; i++)
		{
			*(eform + 6 * n + i) = *(td + i);
			*(eform + 6 * n + 3 + i) = *(tr + i);
		}
		free(td);
		free(tr);
		freematrix(rmtx, 3);
		freematrix(trmtx, 3);
	}
	free(c);
	free(d);
	free(r);

	return eform;
}/*extractlocalcoord*/


double* extractshelldisplacement(struct oshell shell, double* ddisp)
/*EXTRACT ELEMENT DEFORMATION{dU} FROM GLOBAL VECTOR.*/
{
	long int i, loffset;
	int n, nnod;
	double* d;

	nnod = shell.nnod;
	d = (double*)malloc(6 * nnod * sizeof(double));
	for (n = 0; n < nnod; n++)
	{
		for (i = 0; i < 6; i++)
		{
			loffset = 6 * (shell.node[n]->loff) + i;
			*(d + 6 * n + i) = *(ddisp + loffset);
		}
	}
	return d;
}/*extractshelldisplacement*/


double*  extractdeformation(double* eforminit, double* eform, int nnod)
/*EXTRACT LOCAL ELEMENT DEFORMATION FROM GLOBAL.*/
/*UPDATE PSUEDO-ROTATION VECTOR*/
{
	int n, i;
	double* edisp;
	double* r, * rinit, * rvct;
	double** rmtx, ** rh, ** rt, ** rtt;

	edisp = (double*)malloc(6*nnod * sizeof(double));
	r     = (double*)malloc(3 * sizeof(double));
	rinit = (double*)malloc(3 * sizeof(double));

	for (n = 0; n < nnod; n++)
	{
		for (i = 0; i < 3; i++)
		{
			*(r + i)     = *(eform     + 6 * n + 3 + i);
			*(rinit + i) = *(eforminit + 6 * n + 3 + i);
		}

		rh = rotationmtx(r);
		rt = rotationmtx(rinit);
		rtt = matrixtranspose(rt, 3);
		rmtx = matrixmatrix(rh, rtt, 3);
		rvct = rotationvct(rmtx);

		for (i = 0; i < 3; i++)
		{
			*(edisp + 6 * n + i)     = *(eform + 6 * n + i) - *(eforminit + 6 * n + i);
			*(edisp + 6 * n + 3 + i) = *(rvct + i);
		}

		freematrix(rh, 3);
		freematrix(rt, 3);
		freematrix(rtt, 3);
		freematrix(rmtx, 3);
		free(rvct);
	}

	free(r);
	free(rinit);
	return edisp;
}

#if 0
void initialformCR(struct onode* nodes, double* ddisp, int nnode)
/*INITIAL FORMATION INTO DISPLACEMENT.WITHOUT NODE CODE.*/
{
	int i, j;

	for (i = 0; i < nnode; i++)
	{
		for (j = 0; j < 3; j++)
		{
			*(ddisp + 6 * i + j) = (nodes + i)->d[j];
			if ((nodes + i)->r[j] != NULL)
			{
				*(ddisp + 6 * i + 3 + j) = (nodes + i)->r[j];
			}
			else
			{
				*(ddisp + 6 * i + 3 + j) = 0.0;
			}
		}
	}


	return;
}/*initialformCR*/
#endif

void LDLmode(struct gcomponent* gmtx, struct oconf* confs, int*  m, int nmode, double** mode, double* norm, double* dm, long int msize)
{
	double data;
	double LDLnorm;
	struct gcomponent* gcomp1;
	long int i, j, k;
	long int ii;
	for (i = 0; i < nmode; i++)
	{
		ii = *(m + i);
		*(dm + i) = (gmtx + ii)->value;
		LDLnorm = 0.0;
		for (j = msize - 1; j >= 0; j--)                                 /*BACKWARD.*/
		{
			if (j == ii)
			{
				*(*(mode + i) + j) = 1.0;
			}
			else
			{
				*(*(mode + i) + j) = 0.0;
			}

			if ((confs + j)->iconf == 0)
			{
				data = *(*(mode + i) + j);
				gcomp1 = (gmtx + j); /*DIAGONAL.*/
				while (gcomp1->down != NULL) /*DOWNWARD.*/
				{
					gcomp1 = gcomp1->down;
					k = gcomp1->m;
					if ((confs + k - 1)->iconf == 0) /*FREE*/
					{
						data -= (gcomp1->value) * (*(*(mode + i) + k - 1));
					}
				}
				*(*(mode + i) + j) = data;
				LDLnorm += data * data;
			}
		}
		//*(norm + i) = sqrt(LDLnorm);
		LDLnorm = sqrt(LDLnorm);
		for (j = 0; j < msize; j++)                                 /*BACKWARD.*/
		{
			//*(*(mode + i) + j) /= *(norm + i); ]
			*(*(mode + i) + j) /= LDLnorm;
		}

	}
	return;
}

double shellvolume(struct oshell shell, double** drccos, double area)
{
	double volume;
	struct onode node = *(shell.node[0]);

	volume = (*(*(drccos + 2) + 0) * node.d[GX]
		+ *(*(drccos + 2) + 1) * node.d[GY]
		+ *(*(drccos + 2) + 2) * node.d[GZ]) * area / 3.0;
	return volume;
}/*shellvolume*/

void assemshellvolume(struct oshell* shells, int nshell, double* ddisp, double* volume)
{
	struct oshell shell;
	int i,j,ii,jj;
	int nnod;
	double* gforminit, * eforminit;                      /*DEFORMED COORDINATION OF ELEMENT*/
	double* gform, * eform;                      /*INITIAL COORDINATION OF ELEMENT*/
	double** drccos,** drccosinit;                           /*MATRIX*/
	double area;

	for (i = 1; i <= nshell; i++)
	{
		inputshell(shells, NULL, i - 1, &shell);
		nnod = shell.nnod;


		/*INITIAL CONFIGURATION*/
		/*
		for (ii = 0; ii < nnod; ii++)
		{
			inputnode(iform, shell.node[ii]);
		}
		drccosinit = shelldrccos(shell, &area);
		gforminit = extractshelldisplacement(shell, iform);
		eforminit = extractlocalcoord(gforminit,drccosinit,nnod);
		*/

		/*DEFORMED CONFIGFURATION*/
		for (ii = 0; ii < nnod; ii++)
		{
			inputnode(ddisp, shell.node[ii]);
		}
		drccos = shelldrccos(shell, &area);
		//gform = extractshelldisplacement(shell, ddisp);                     /*{Xg+Ug}*/
		//eform = extractlocalcoord(gform,drccos,nnod); 			       	    /*{Xe+Ue}*/

		*volume += shellvolume(shell, drccos, area);                   		/*VOLUME*/

		freematrix(drccos, 3);
		/*
		freematrix(drccosinit, 3);

		free(eforminit);
		free(gforminit);
		free(eform);
		free(gform);
		*/

	}
	return;
}
/*
void pressureupdate(double volume, double mol, double loadfactor)
{

	 loadfactor = mol*8.3*(20+273)/volume;
	 return;
}
*/

double equilibriumcurvature(double* weight, double* lapddisp, double laploadfactor, double* dup, int msize)
{
	int i;
	double cos, dot, len, lastlen, curvature;
	char string[40];

	dot = laploadfactor;
	len = 1.0;
	lastlen = laploadfactor * laploadfactor;
	for (i = 0; i < msize; i++)
	{
		dot += *(weight + i) **(dup + i) **(weight + i) **(lapddisp + i);
		len += *(weight + i) **(dup + i) **(weight + i) **(dup + i);
		lastlen += *(weight + i) **(lapddisp + i) **(weight + i) **(lapddisp + i);
	}
	cos = abs(dot / sqrt(len * lastlen));
	curvature = acos(cos) / sqrt(lastlen);
	return curvature;
}



void dbgvct(double* vct, int size, int linesize, const char* str)
{
	FILE *fout = fopen("hogan.dbg", "a");
	if (fout == NULL)return;

	if (str != NULL)
	{
		fprintf(fout, "%s\n", str);
	}
	if (linesize == NULL)
	{
		linesize = 1;
	}
	for (int i = 0; i < size; i++)
	{
		fprintf(fout, "%e ", *(vct+i));
		if ((i + 1) % linesize == 0)fprintf(fout, "\n");
	}
	fclose(fout);
}

void dbgmtx(double** mtx, int rows, int cols, const char* str)
{
	FILE *fout = fopen("hogan.dbg", "a");
	if (fout == NULL)return;

	if (str != NULL)
	{
		fprintf(fout, "%s\n", str);
	}
	if (cols== NULL)
	{
		cols = rows;
	}
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			fprintf(fout, "%e ", *(*(mtx+i)+j));
		}
		fprintf(fout, "\n");
	}
	fclose(fout);
}


