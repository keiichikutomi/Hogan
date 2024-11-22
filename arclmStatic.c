int arclmStatic(struct arclmframe* af);

void LDLmode(struct gcomponent* gmtx, struct oconf* confs, int* m, int nmode, double** mode, double* norm, double* dm, long int msize);
double equilibriumcurvature(double* weight, double* lapddisp, double laploadfactor, double* dup, int msize);


void dbgstr(const char* str);
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


int arclmStatic(struct arclmframe* af)
{
	DWORDLONG memory0,memory1;
	char dir[]=DIRECTORY;
	char string[400],fname[100];
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
	struct gcomponent *gmtx, * g, * p, * gcomp1, *epsgmtx, *gmtxcpy, *dgmtx;
	double gg;
	double sign, determinant;
	long int nline;

	int USINGEIGENFLAG = 0;

	/*GLOBAL VECTOR*/
	double* ddisp, * iform;
	double* givendisp;
	double* gvct;
	double* funbalance, * fexternal, * finternal, * freaction;
	double* fbaseload, * fdeadload, * fpressure,* fgivendisp;
	double* fswitching;



	///FOR ARC-LENGTH PARAMETER///
	int nlap, laps;
	int iteration;
	int maxiteration = 20;
	double tolerance = 1.0e-4;
	double residual;
	double gvctlen;
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

	double volume = 0.0;


	/*ANALYSIS MODE*/
	int outputmode = 0;/*0:ConvergedLaps.1:AllLaps.*/
	int pinpointmode = 0;/*0:NotPinpointing.1:BisecPinpointing.2:ExtendedSystemPinpointing.*/
	int BCLFLAG = 0;/*FLAG FOR BUCKLING DETECTION*/
	int UNLOADFLAG = 0;
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

	/***GLOBAL MATRIX*/
	gmtx = (struct gcomponent*)malloc(msize * sizeof(struct gcomponent));/*DIAGONALS OF GLOBAL MATRIX.*/
	epsgmtx = (struct gcomponent*)malloc(msize * sizeof(struct gcomponent));/*DIAGONALS OF GLOBAL MATRIX.*/

	/***GLOBAL VECTOR*/
	gvct = (double *)malloc(msize * sizeof(double));/*INCREMENTAL GLOBAL VECTOR.*/
	due = (double*)malloc(msize * sizeof(double));
	dup = (double*)malloc(msize * sizeof(double));

	/*FORCE VECTOR INITIALIZATION*/
	funbalance = (double*)malloc(msize * sizeof(double));          /*UNBALANCED INTERNAL FORCE VECTOR.*/
	freaction = (double*)malloc(msize * sizeof(double));           /*EXTERNAL FORCE VECTOR.*/
	fexternal = (double*)malloc(msize * sizeof(double));           /*EXTERNAL FORCE VECTOR.*/
	finternal = (double*)malloc(msize * sizeof(double));           /*INTERNAL FORCE VECTOR.*/
	fbaseload = (double*)malloc(msize * sizeof(double));           /*BASE LOAD VECTOR.*/
	fdeadload = (double*)malloc(msize * sizeof(double));           /*DEAD LOAD VECTOR.*/
	fpressure = (double*)malloc(msize * sizeof(double));           /*PRESSURE VECTOR.*/
	fgivendisp = (double*)malloc(msize * sizeof(double));           /*BASE LOAD VECTOR.*/
	fswitching = (double*)malloc(msize * sizeof(double));

	/*POSITION VECTOR INITIALIZATION*/
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

	nlap = 0;
	iteration = 0;
	residual = 0.0;

	assemconf(confs,fdeadload,1.0,nnode);
	assemgivend(confs,givendisp,1.0,nnode);


	while (nlap <= laps && ENDFLAG == 0)
	{
		std::vector<Triplet> Ktriplet;
		std::vector<Triplet> Mtriplet;

		SparseMatrix Kglobal(msize, msize);
		/*EXECUTE BY WIN64. AMD ORDERING IS AVAILABLE ONLY BY WIN64*/
		//Eigen::SimplicialLDLT<SparseMatrix,Eigen::Lower,Eigen::NaturalOrdering<int>> solver;
		Eigen::SimplicialLDLT<SparseMatrix> solver;

		for (i = 1; i <= msize; i++)
		{
			g = (gmtx + (i - 1))->down;
			while (g != NULL)
			{
				p = g;
				g = g->down;
				free(p);
			}
			ginit.m = (unsigned short int)i;
			*(gmtx + (i - 1)) = ginit;
		}
		for (i = 1; i <= msize; i++)
		{
			g = (epsgmtx + (i - 1))->down;
			while (g != NULL)
			{
				p = g;
				g = g->down;
				free(p);
			}
			ginit.m = (unsigned short int)i;
			*(epsgmtx + (i - 1)) = ginit;
		}


		for (i = 0; i < msize; i++)	 /*GLOBAL VECTOR INITIALIZATION.*/
		{
			*(finternal + i) = 0.0;
			*(fexternal + i) = 0.0;
			*(funbalance + i) = 0.0;
			*(freaction + i) = 0.0;
			*(fbaseload + i) = 0.0;
			*(fpressure + i) = 0.0;
		}
		comps = msize; /*INITIAL COMPONENTS = DIAGONALS.*/

		volume = 0.0;
		assemshellvolume(shells, nshell, ddisp, &volume);

		/*POST PROCESS*/
		if(USINGEIGENFLAG==1)
		{
		  //assemelemEx(elems, melem, nelem, constraintmain, confs, NULL, NULL, iform, ddisp, finternal, fpressure);
		  //assemshellEx(shells, mshell, nshell, constraintmain, confs, NULL, NULL, iform, ddisp, finternal, fpressure);
		}
		else
		{
		  assemelem(elems, melem, nelem, constraintmain, NULL, NULL, iform, ddisp, finternal, fpressure);
		  assemshell(shells, mshell, nshell, constraintmain, NULL, NULL, iform, ddisp, finternal, fpressure);
		}


		if(/*UNLOADFLAG==1*/0)
		{
			for (i = 0; i < msize; i++)
			{
				if(nlap==1 && iteration==1)
				{
					*(fdeadload + i) = -*(finternal + i);
				}
				*(dup + i) = *(fdeadload + i)+*(fpressure + i);
				*(fexternal + i) = - *(fdeadload + i) + loadfactor * *(dup + i);
				*(funbalance + i) = *(fexternal + i) - *(finternal + i);            /*funbalance:UNBALANCED FORCE -{E}.*/
				if ((confs + i)->iconf == 1)
				{
					*(freaction + i) = *(funbalance + i);
					*(funbalance + i) = 0.0;
					*(dup + i) = 0.0;
				}

				*(due + i) = *(funbalance + i);
			}


		}
		else
		{
			for (i = 0; i < msize; i++)
			{
				*(fbaseload + i) = *(fdeadload + i)+*(fpressure + i);
				*(fexternal + i) = loadfactor * *(fbaseload + i);
				*(funbalance + i) = *(fexternal + i) - *(finternal + i);            /*funbalance:UNBALANCED FORCE -{E}.*/
				if ((confs + i)->iconf == 1)
				{
					*(freaction + i) = *(funbalance + i);
					*(funbalance + i) = 0.0;
					*(fbaseload + i) = 0.0;
				}
				*(dup + i) = *(fbaseload + i);
				*(due + i) = *(funbalance + i);
			}

		}

		/*CONVERGENCE JUDGEMENT.*/
		residual = vectorlength(funbalance,msize);

		if (residual < tolerance || iteration > maxiteration-1)
		{
		  nlap++;
		  iteration = 0;
		}
		iteration++;
		setincrement((wmenu.childs+2)->hwnd,laps,nlap,maxiteration,iteration);

		if(iteration==1)
		{
			clearwindow(*(wdraw.childs+1));
			drawarclmframe((wdraw.childs+1)->hdcC,(wdraw.childs+1)->vparam,*af,0,ONSCREEN);
			overlayhdc(*(wdraw.childs + 1), SRCPAINT);                  /*UPDATE DISPLAY.*/

			af->nlaps = nlap;

			/*LOCAL INCREMENTAL IN THIS LAP.*/
			for(i = 0; i < nshell; i++)
			{
			  outputmemoryshell(shells,mshell,i);
			}
		}

		/*STIFFNESS ASSEMBLAGE*/
		if(USINGEIGENFLAG==1)
		{
		  //assemelemEx(elems, melem, nelem, constraintmain, confs, NULL, Ktriplet, iform, ddisp, NULL, NULL);
		  //assemshellEx(shells, mshell, nshell, constraintmain, confs, NULL, Ktriplet, iform, ddisp, NULL, NULL);

		  /*GLOBAL MATRIX USING EIGEN*/
		  Kglobal.reserve(msize*msize);
		  Kglobal.setFromTriplets(Ktriplet.begin(), Ktriplet.end());//SPARSE MATRIX FROM TRIPLET
		}
		else
		{
		  assemelem(elems, melem, nelem, constraintmain, NULL, gmtx, iform, ddisp, NULL, NULL);
		  assemshell(shells, mshell, nshell, constraintmain, NULL, gmtx, iform, ddisp, NULL, NULL);
		}

		/*OUTPUT.*/
		if ((outputmode == 0 && iteration == 1) || outputmode == 1)
		{
			sprintf(string, "LAP: %5ld / %5ld ITERATION: %5ld\n", nlap, laps, iteration);

			/*OUTPUT FOR GLOBAL*/
			//dbgstr(string);

			fprintf(fdsp, string);
			fprintf(finf, string);
			fprintf(fexf, string);
			fprintf(fubf, string);
			fprintf(frct, string);
			outputdisp(ddisp, fdsp, nnode, nodes);                    /*FORMATION OUTPUT.*/
			outputdisp(finternal, finf, nnode, nodes);                /*FORMATION OUTPUT.*/
			outputdisp(fexternal, fexf, nnode, nodes);                /*FORMATION OUTPUT.*/
			outputdisp(funbalance, fubf, nnode, nodes);               /*FORMATION OUTPUT.*/
			outputdisp(freaction, frct, nnode, nodes);                /*FORMATION OUTPUT.*/

			/*OUTPUT FOR LOCAL*/
			fprintf(fstr, string);
			fprintf(fene, string);
			for(i = 0; i < nshell; i++)
			{
				//fprintf(fene, "%5ld %e %e %e\n", (shells+i)->code, (mshell+i)->SEp, (mshell+i)->SEb, (mshell+i)->SE);
				fprintf(fstr, "%5ld %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n", (shells+i)->code,
				((shells+i)->gp[0]). stress[0],((shells+i)->gp[0]). stress[1],((shells+i)->gp[0]). stress[2],((shells+i)->gp[0]). stress[3],((shells+i)->gp[0]). stress[4],((shells+i)->gp[0]). stress[5],
				((shells+i)->gp[0]).estrain[0],((shells+i)->gp[0]).estrain[1],((shells+i)->gp[0]).estrain[2],((shells+i)->gp[0]).estrain[3],((shells+i)->gp[0]).estrain[4],((shells+i)->gp[0]).estrain[5],
				((shells+i)->gp[0]).pstrain[0],((shells+i)->gp[0]).pstrain[1],((shells+i)->gp[0]).pstrain[2],((shells+i)->gp[0]).pstrain[3],((shells+i)->gp[0]).pstrain[4],((shells+i)->gp[0]).pstrain[5],
				((shells+i)->gp[0]).qn,((shells+i)->gp[0]).qm,((shells+i)->gp[0]).qnm,((shells+i)->gp[0]).y,((shells+i)->gp[0]).yinit,((shells+i)->gp[0]).f[0],((shells+i)->gp[0]).f[1]
				);
			}
		}

		if(iteration==1 && USINGEIGENFLAG==0)
		{
			for (i = 0; i < msize; i++)
			{
				*(fgivendisp + i) = *(givendisp + i);
			}
			modifygivend(gmtx,fgivendisp,confs,nnode);

			for (i = 0; i < msize; i++)
			{
				*(dup + i) += *(fgivendisp + i);
			}
		}


		/*SOLVE [K]*/
		if(USINGEIGENFLAG==1)
		{
			solver.compute(Kglobal);

			Eigen::VectorXd D = solver.vectorD();
			determinant= 0.0;
			sign = 0;
			for (i = 0; i < msize; i++)
			{
			  if((confs+i)->iconf!=1)
			  {
				if(D(i)<0.0) sign += 1;
				determinant += log10(fabs(D(i)));
			  }
			}

			if (solver.info() != Eigen::Success)sign=-1;
		}
		else
		{
			if(iteration==1 && neig!=0)gmtxcpy=copygcompmatrix(gmtx,msize);

			nline = croutlu(gmtx, confs, msize, &determinant, &sign, gcomp1);
			if (sign < 0.0)
			{
				for (ii = 1; ii <= msize; ii++)
				{
					gg = 0.0;
					gread(gmtx, ii, ii, &gg);

					if (gg < 0.0)
					{
						sprintf(string, "INSTABLE TERMINATION AT NODE %ld.\n",(nodes + int((ii - 1) / 6))->code);
						fprintf(ffig, "%s", string);
						errormessage(string);
					}
				}
			}

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
		}

		/*COULD NOT SOLVE*/
		if (sign < 0.0)
		{
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
			free(fdeadload);
			free(fgivendisp);
			free(fswitching);
			free(due);
			free(dup);
			free(lastddisp);
			free(lastgvct);
			free(lastpivot);
			free(lapddisp);

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

		sprintf(string, "LAP: %4d ITER: %2d {LOAD}= % 5.8f {RESD}= %1.6e {DET}= %8.5f {SIGN}= %2.0f {BCL}= %1d {EPS}= %1.5e {V}= %8.5f\n",
				nlap, iteration, loadfactor, residual, determinant, sign, BCLFLAG, 0.0, volume);
		fprintf(ffig, "%s", string);
		errormessage(string);

		if (iteration == 1 && (sign - lastsign) != 0 && nlap != 1 && pinpointmode == 1)/*BUCKLING DETECTED*/
		{
			sprintf(string, "BUCKLING DITECTED LAP: %4d ITER: %2d\n", nlap, iteration);
			fprintf(fbcl, "%s", string);

			/*BISEC INITIAL*/
			nextloadfactor = loadfactor;
			for (ii = 0; ii < msize; ii++)
			{
				*(nextddisp + ii) = *(ddisp + ii);
			}

			LR = 1.0;
			LL = 0.0;

			while(1)
			{
				LM = 0.5 * (LL + LR);

				loadfactor = (1 - LM) * lastloadfactor + LM * nextloadfactor;
				for (ii = 0; ii < msize; ii++)
				{
					*(ddisp + ii) = (1 - LM) * *(lastddisp + ii) + LM * *(nextddisp + ii);
				}

				for (i = 1; i <= msize; i++)
				{
					g = (gmtx + (i - 1))->down;
					while (g != NULL)
					{
						p = g;
						g = g->down;
						free(p);
					}
					ginit.m = (unsigned short int)i;
					*(gmtx + (i - 1)) = ginit;
				}

				assemelem(elems, melem, nelem, constraintmain, NULL, gmtx, iform, ddisp, NULL, NULL);
				assemshell(shells, mshell, nshell, constraintmain, NULL, gmtx, iform, ddisp, NULL, NULL);

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
							sprintf(string, "INSTABLE TERMINATION AT NODE %ld.",(nodes + int((ii - 1) / 6))->code);
							errormessage(" ");
							errormessage(string);
							if (fonl != NULL) fprintf(fonl, "%s\n", string);
						}
					}
					return 1;
				}


				/*BISECTION PIN-POINTING*/
				if (lastsign == sign)
				{
					LL = LM;
				}
				else
				{
					LR = LM;
				}
				if(LR - LL < biseceps)break;/*CONVERGED*/
			}

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
				loadfactor = nextloadfactor;
				for (ii = 0; ii < msize; ii++)
				{
					*(ddisp + ii) = *(nextddisp + ii);
				}
			}
			else/*PATH=SWITCHING FOR BIFURCATION*/
			{
				/*UNDER CONSTRUCTION*/
				loadfactor = lastloadfactor;
				for (ii = 0; ii < msize; ii++)
				{
					*(ddisp + ii) = *(lastddisp + ii);
				}
			}

		}
		/*NORMAL CALCULATION IN CASE OF REGULAR*/

		if (iteration == 1)/*PREDICTOR CALCULATION*/
		{
			if(USINGEIGENFLAG==1)
			{
			  Vector P = Vector::Zero(msize);
			  for(i=0;i<msize;i++)P(i)=*(fbaseload+i);
			  Vector Up = solver.solve(P);
			  for(i=0;i<msize;i++)*(dup+i)=Up(i);
			  if (solver.info() != Eigen::Success)return -1;
			}
			else
			{
			  nline = forwardbackward(gmtx, dup, confs, msize, gcomp1);
			}

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
				*(lapddisp + ii) = 0.0;
			}
			lastloadfactor = loadfactor;/*LOAD FACTOR AT CONVERGED POINT*/
			lastlambda = lambda;/*INCREMENTAL LOAD FACTOR OF PREDICTOR*/
			laploadfactor = 0.0;
			lastsign = sign;/*SIGN AT CONVERGED POINT*/
			BCLFLAG = 0;

		}
		else/*CORRECTOR CALCULATION*/
		{
			if(USINGEIGENFLAG==1)
			{
				Vector P = Vector::Zero(msize);
				for(i=0;i<msize;i++)P(i)=*(fbaseload+i);
				Vector Up = solver.solve(P);
				for(i=0;i<msize;i++)*(dup+i)=Up(i);
				if (solver.info() != Eigen::Success)return -1;

				Vector E = Vector::Zero(msize);
				for(i=0;i<msize;i++)E(i)=*(funbalance+i);
				Vector Ue = solver.solve(E);
				for(i=0;i<msize;i++)*(due+i)=Ue(i);
				if (solver.info() != Eigen::Success)return -1;
			}
			else
			{
				nline = forwardbackward(gmtx, dup, confs, msize, gcomp1);
				nline = forwardbackward(gmtx, due, confs, msize, gcomp1);
			}



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
			/*HYPERSPHERE CONSTRAINT : CRISFIELD'S METHOD*/
			/*USING PREDICTOR DATA*/

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


		/*TWO-POINTS BUCKLING ANALYSIS.*/
		if(iteration == 1 && neig!=0 && nlap%10==0)
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
			assemelem(elems, melem, nelem, constraintmain, NULL, epsgmtx, iform, epsddisp, NULL, NULL);
			assemshell(shells, mshell, nshell, constraintmain, NULL, epsgmtx, iform, epsddisp, NULL, NULL);

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

		laploadfactor += lambda;
		updaterotation(lapddisp, gvct, nnode);

		loadfactor += lambda;
		updaterotation(ddisp, gvct, nnode);


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
#if 0
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
            free(fdeadload);
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
	free(fdeadload);
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



		#if 0

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
		else if (BCLFLAG == 1)/*BUCKLING DETECTED.PIN-POINTING EXCUTION.*/
		{




			if (pinpointmode == 2)/*PIN-POINTING BASED ON EXTENDED SYSTEM.*/
			{

				BCLFLAG = 1;
				sprintf(string, "BUCKLING DITECTED LAP: %4d ITER: %2d\n", nlap, iteration);
				fprintf(fbcl, "%s", string);
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

			}
			iteration++;
		}
		else if (BCLFLAG == 2)
		{

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


			assemshell(shells, mshell, nshell, constraintmain,NULL, gmtx,iform, epsddisp, epsfinternal, epsfexternal);

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


void dbgstr(const char* str)
{
	FILE *fout = fopen("hogan.dbg", "a");
	if (fout == NULL)return;
	fprintf(fout, "%s", str);
	fclose(fout);
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


