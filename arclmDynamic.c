/* ========================================================= */
/* PROGRAM GNSHN CR FOR OUTER SURFACE OF LUNAR MARS PROJECT  */
/* DYNAMIC ANALYSIS FOR LARGE DEFORMATION                    */
/* USING CR FORMULATION                                      */
/* CODED BY KEIICHI KUTOMI SINSE 2024.05.26                  */
/* ========================================================= */

int arclmStatDyna(struct arclmframe* af);
int arclmDynamic(struct arclmframe* af);
double timestepcontrol(double ddt, double* lapddisp_m, double* udd_m, double* lastudd_m, double* lastlastudd_m, double beta, int msize);
double loadfactormap(double time);

int arclmStatDyna(struct arclmframe* af)
{
	double targetload;
	double ENDFLAG = 0;
	targetload = 1e+10;
	while(1)
	{
		errormessage("arclmStatic START.");
		ENDFLAG = arclmStatic(af);
		if(af->loadfactor > targetload || ENDFLAG == 1)break;
		errormessage("arclmDynamic START.");
		ENDFLAG = arclmDynamic(af);
		if(af->loadfactor > targetload || ENDFLAG == 1)break;
	}
	errormessage("arclmStatDyna END.");
	return 0;
}


int arclmDynamic(struct arclmframe* af)
{
	DWORDLONG memory0,memory1;
	char dir[] = DIRECTORY;
	char string[400], fname[100];
	clock_t t0;
	FILE *fin, *fonl, * fdsp, * fexf, * finf, * fubf, * frct, * finr, * fstr, * fene, * ffig, * fbcl, * feig, *fout, *fplst;         /*FILE 8 BYTES*/
	FILE *fvel,*facc;

	int i, j, ii, jj;

	/*MODEL SIZE*/
	int nnode, nelem, nshell, nsect, nconstraint;
	long int msize,csize;

	///ARCLMFRAME///
	struct osect* sects;
	struct onode* nodes;
	struct onode* ninit;
	struct owire* elems;
	struct oshell* shells;
	struct oconf* confs;
	struct memoryelem* melem;
	struct memoryshell* mshell;
	long int* constraintmain;
	struct oconstraint* constraints;
    double *nmass;

	/*GLOBAL MATRIX*/
	struct gcomponent ginit = {0,0,0.0,NULL};
	struct gcomponent* gmtx, * g, * p, * gcomp1;/*GLOBAL MATRIX*/
	struct gcomponent* gmtx2;
	double gg;
	double sign ,determinant;
	double sign2,determinant2;
	long int nline,nline2;

	int USINGEIGENFLAG = 1;

	/*GLOBAL VECTOR*/
	double* iform,* ddisp,* lastddisp,* lapddisp,* lapddisp_m,*givendisp;
	double* lambda, * lastlambda;
	double* funbalance, * fexternal, * finternal, * freaction;
	double* fbaseload, * fdeadload, * fpressure,* finertial,* fdamping, * fgivendisp;
	double* fconstraint;
	double* constraintvct;

	double* funbalance_m;

	double* gvct;

	double* lastlastudd_m;
	double* lastud_m,* lastudd_m;
	double* ud_m,* udd_m;

	double* rightud_m;/*VELOCITY AT t+dt/2 FOR EXPLICIT*/
	double* leftud_m;/*VELOCITY AT t-dt/2 FOR EXPLICIT*/
	double* mvct, * cvct;/*MASS & DAMPING VECTOR(DIAGONAL MATRIX)*/

	///FOR ITERATIVE CALCULATION///
	int nlap = 1;
	int beginlap = 1;
    int targetlap;
	int laps = 1000;
	int iteration = 1;
	int maxiteration = 20;
	double tolerance = 1.0E-2;
	double residual = 0.0;
	double constraintresidual = 0.0;
	double gvctlen = 0.0;

	double ddt = 0.0001;/*TIME INCREMENT[sec]*/
	double initialddt = ddt;
	double time = 0.0;/*TOTAL TIME[sec]*/

	double initialloadfactor = 0.0;
	double loadfactor = 0.0;
	double lastloadfactor = 0.0;
	double loadlambda = 0.0;

	double volume = 0.0;

	/* ENERGY */
	double Wkt = 0.0;
	double Wet = 0.0;
	double Wpt = 0.0;
	double Wot = 0.0;
	double lastWkt = 0.0;

	/*ANALYSIS MODE*/
	int outputmode   = 0;/*0:ConvergedLaps.1:AllLaps.*/
	int pinpointmode = 0;/*0:NotPinpointing.1:BisecPinpointing.2:ExtendedSystemPinpointing.*/
	int solver = 0;
	int method = 0;
	double rho,alpham,alphaf,xi,beta,gamma;
	/*
	[solver]
	0:INPLICIT
	1:EXPLICIT

	[method]
	solver  INPLICIT         EXPLICIT
	0      | Newmark-beta    | Centered Difference(CD)
	1      | HTT-alpha       | Runge-Kutta(RK4)
	2      | Energy Momentum | Dynamic Relaxation(CD)
	*/

	int STATDYNAFLAG = 0;
	int UNLOADFLAG = 0;
	int RELAXATION = 1;
	int PEAKFLAG = 0;
	int DIVFLAG = 0;

	/*ANALYSIS TERMINATION*/
	int ENDFLAG = 0;
	int fnode=NULL,fnodedrc=NULL;
	double fnodemin, fnodemax;

	/*FOR READING ANALYSISDATA*/
	FILE *fdata;
	int nstr, pstr, readflag;
	char **data;
	char filename[256],inputfilename[256],outputfilename[256];
	char description[256];
	char* dot;

	memory0 = availablephysicalmemoryEx("INITIAL:");   /*MEMORY AVAILABLE*/

#if 0
	fin = fgetstofopenII(dir, "r", (wdraw.childs+1)->inpfile);
	if(fin==NULL)
	{
	  errormessage("ACCESS IMPOSSIBLE.");
	  return 1;
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



	strcpy(filename, (wdraw.childs+1)->inpfile);
	dot = strrchr(filename, '.');
	*dot = '\0';

	fdata = fgetstofopenII(dir, "r", "analysisdata.txt");

	if (fdata == NULL)
	{
		errormessage("couldn't open analysisdata.txt\n");
		getincrement((wmenu.childs+2)->hwnd,&laps,&ddt);
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
								readflag=0;
								getincrement((wmenu.childs+2)->hwnd,&laps,&ddt);
							}
							else
							{
								strcpy(filename, *(data + pstr));
							}
						}
						if (!strcmp(*(data + pstr), "INPUTFILENAME"))
						{
							pstr++;
							strcpy(inputfilename, *(data + pstr));
						}
						if (!strcmp(*(data + pstr), "OUTPUTFILENAME"))
						{
							pstr++;
							strcpy(outputfilename, *(data + pstr));
						}
						if (!strcmp(*(data + pstr), "DESCRIPTION"))
						{
							pstr++;
							strcpy(description, *(data + pstr));
							strcat(filename, "_");
							strcat(filename, description);
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
						if (!strcmp(*(data + pstr), "TIMEINCREMENT"))
						{
							pstr++;
							ddt = (double)strtod(*(data + pstr), NULL);
							initialddt = ddt;
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
							initialloadfactor = (double)strtod(*(data + pstr), NULL);
						}
						if (!strcmp(*(data + pstr), "TARGETLAP"))
						{
							pstr++;
							targetlap = (int)strtol(*(data + pstr), NULL, 10);
						}
						if (!strcmp(*(data + pstr), "SOLVER"))
						{
							pstr++;
							solver = (int)strtol(*(data + pstr), NULL, 10);
						}
						if (!strcmp(*(data + pstr), "METHOD"))
						{
							pstr++;
							method = (int)strtol(*(data + pstr), NULL, 10);
						}
						if (!strcmp(*(data + pstr), "STATDYNA"))
						{
							STATDYNAFLAG = 1;
						}
						if (!strcmp(*(data + pstr), "UNLOAD"))
						{
							UNLOADFLAG = 1;
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
    fclose(fdata);



	sprintf(string,"FILENAME : %s\n LAPS = %d\n TIME INCREMENT = %lf\n", filename, laps, ddt);
	errormessage(string);

	if (outputmode == 0)errormessage("OUTPUT CONVERGED RESULT\n");
	if (outputmode == 1)errormessage("OUTPUT ALL RESULT\n");


	if(method==1)/*HTT-alpha'S PARAMETER.*/
	{
		rho = 0.904762;
		alpham = 0.0;  /*MID-POINT USED TO EVALUATE INERTIAL FORCE*/
		alphaf = (1.0-rho)/(1.0+rho);        /*MID-POINT USED TO EVALUATE INTERNAL FORCE*/
		xi = 0.0;   	 				/*NUMERICAL DISSIPATION COEFFICIENT*/
		beta = pow(1+alphaf,2)/4.0;
		gamma = 0.5+alphaf;
	}
	else if(method==2)/*ENERGY MOMENTUM METHOD'S PARAMETER.*/
	{
		rho = 1.0;
		//rho = 0.8;
		alpham = (2.0*rho-1)/(rho+1);  /*MID-POINT USED TO EVALUATE INERTIAL FORCE*/
		alphaf = rho/(rho+1);        /*MID-POINT USED TO EVALUATE INTERNAL FORCE*/
		//xi = 0.0;   	 				/*NUMERICAL DISSIPATION COEFFICIENT*/
		xi = 0.056;   	 				/*NUMERICAL DISSIPATION COEFFICIENT*/
		beta = 0.25*pow(1-alpham+alphaf,2);
		gamma = 0.5-alpham+alphaf;
	}
	else/*NEWMARK METHOD'S PARAMETER.*/
	{
		rho = 1.0;
		alpham = 0.0;  /*MID-POINT USED TO EVALUATE INERTIAL FORCE*/
		alphaf = 0.0;        /*MID-POINT USED TO EVALUATE INTERNAL FORCE*/
		xi = 0.0;   	 				/*NUMERICAL DISSIPATION COEFFICIENT*/
		beta = 1.0/pow(rho+1,2);
		gamma = (3.0-rho)/(2.0*(rho+1.0));
	}
	sprintf(string,"alpham=%.3f alphaf=%.3f beta=%.3f gamma=%.3f ",alpham,alphaf,beta,gamma);
	errormessage(string);



	///INPUT FILE OPEN & OUTPUT FILE SETTING///

	if(STATDYNAFLAG == 1 && af->nlaps != NULL)
	{
		snprintf(fname, sizeof(fname), "%s.%s", filename, "dsp");
		fdsp = fopen(fname, "a");
		snprintf(fname, sizeof(fname), "%s.%s", filename, "inf");
		finf = fopen(fname, "a");
		snprintf(fname, sizeof(fname), "%s.%s", filename, "exf");
		fexf = fopen(fname, "a");
		snprintf(fname, sizeof(fname), "%s.%s", filename, "ubf");
		fubf = fopen(fname, "a");
		snprintf(fname, sizeof(fname), "%s.%s", filename, "rct");
		frct = fopen(fname, "a");
		snprintf(fname, sizeof(fname), "%s.%s", filename, "inr");
		finr = fopen(fname, "a");
		snprintf(fname, sizeof(fname), "%s.%s", filename, "str");
		fstr = fopen(fname, "a");
		snprintf(fname, sizeof(fname), "%s.%s", filename, "ene");
		fene = fopen(fname, "a");
		snprintf(fname, sizeof(fname), "%s.%s", filename, "onl");
		fonl = fopen(fname, "a");
		snprintf(fname, sizeof(fname), "%s.%s", filename, "fig");
		ffig = fopen(fname, "a");
		snprintf(fname, sizeof(fname), "%s.%s", filename, "bcl");
		fbcl = fopen(fname, "a");
		snprintf(fname, sizeof(fname), "%s.%s", filename, "eig");
		feig = fopen(fname, "a");
		snprintf(fname, sizeof(fname), "%s.%s", filename, "otl");
		fout = fopen(fname, "a");
		snprintf(fname, sizeof(fname), "%s.%s", filename, "vel");
		fvel = fopen(fname, "a");
		snprintf(fname, sizeof(fname), "%s.%s", filename, "acc");
		facc = fopen(fname, "a");
		snprintf(fname, sizeof(fname), "%s.%s", filename, "plst");
		fplst = fopen(fname, "a");
	}
	else
	{
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
		snprintf(fname, sizeof(fname), "%s.%s", filename, "inr");
		finr = fopen(fname, "w");
		snprintf(fname, sizeof(fname), "%s.%s", filename, "str");
		fstr = fopen(fname, "w");
		snprintf(fname, sizeof(fname), "%s.%s", filename, "ene");
		fene = fopen(fname, "w");
		snprintf(fname, sizeof(fname), "%s.%s", filename, "onl");
		fonl = fopen(fname, "w");
		snprintf(fname, sizeof(fname), "%s.%s", filename, "fig");
		ffig = fopen(fname, "w");
		snprintf(fname, sizeof(fname), "%s.%s", filename, "bcl");
		fbcl = fopen(fname, "w");
		snprintf(fname, sizeof(fname), "%s.%s", filename, "eig");
		feig = fopen(fname, "w");
		snprintf(fname, sizeof(fname), "%s.%s", filename, "otl");
		fout = fopen(fname, "w");
		snprintf(fname, sizeof(fname), "%s.%s", filename, "vel");
		fvel = fopen(fname, "w");
		snprintf(fname, sizeof(fname), "%s.%s", filename, "acc");
		facc = fopen(fname, "w");
		snprintf(fname, sizeof(fname), "%s.%s", filename, "plst");
		fplst = fopen(fname, "w");
	}

	t0 = clock();                                                   /*CLOCK BEGIN.*/




#if 0
	/*MEMORY NOT ALLOCATED*/
	free(af->sects);
	free(af->nodes);
	free(af->elems);
	free(af->shells);
	free(af->confs);
	free(af->constraintmain);
	//free(af->iform);
	//free(af->ddisp);
	//free(af->melem);
	//free(af->mshell);

	sects = (struct osect*)malloc(nsect * sizeof(struct osect));
	nodes = (struct onode*)malloc(nnode * sizeof(struct onode));
	ninit = (struct onode*)malloc(nnode * sizeof(struct onode));
	elems = (struct owire*)malloc(nelem * sizeof(struct owire));
	shells = (struct oshell*)malloc(nshell * sizeof(struct oshell));
	confs = (struct oconf*)malloc(msize * sizeof(struct oconf));
	constraintmain = (long int*)malloc(msize * sizeof(long int));
	constraints = (struct oconstraint*)malloc(nconstraint * sizeof(struct oconstraint));
	//iform = (double*)malloc(msize * sizeof(double));		/*INITIAL*/
	//ddisp = (double*)malloc(msize * sizeof(double));		/*LATEST ITERATION*/
	//melem = (struct memoryelem*)malloc(nelem * sizeof(struct memoryelem));
	//mshell = (struct memoryshell*)malloc(nshell * sizeof(struct memoryshell));

	af->sects = sects;
	af->nodes = nodes;
	af->ninit = ninit;
	af->elems = elems;
	af->shells = shells;
	af->confs = confs;
	af->constraintmain = constraintmain;
	af->constraints = constraints;
	af->nmass = nmass;
	//af->iform = iform;
	//af->ddisp = ddisp;
	//af->melem = melem;
	//af->mshell = mshell;

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
	constraintmain = af->constraintmain;
	constraints = af->constraints;
	nmass = af->nmass;
	//iform = af->iform;
	//ddisp = af->ddisp;
	//melem = af->melem;
	//mshell = af->mshell;
#endif


	iform = (double*)malloc(msize * sizeof(double));
	ddisp = (double*)malloc(msize * sizeof(double));
	melem = (struct memoryelem*)malloc(nelem * sizeof(struct memoryelem));
	mshell = (struct memoryshell*)malloc(nshell * sizeof(struct memoryshell));

	af->iform = iform;
	af->ddisp = ddisp;
	af->melem = melem;
	af->mshell = mshell;


	csize = 0;
	for (i = 0; i < nconstraint; i++)
	{
		csize += (constraints+i)->neq;
	}

	constraintvct = (double*)malloc(csize * sizeof(double));

	if(af->lambda==NULL)
	{
		lambda = (double*)malloc(csize * sizeof(double));
		for (i = 0; i < csize; i++)
		{
			*(lambda + i) = 0.0;
		}
		af->lambda = lambda;
	}
	else
	{
		lambda = af->lambda;
	}
	lastlambda = (double*)malloc(csize * sizeof(double));




	/*GLOBAL MATRIX*/
	gmtx = (struct gcomponent*)malloc((msize+csize) * sizeof(struct gcomponent));/*DIAGONALS OF GLOBAL MATRIX.*/
	gmtx2 = (struct gcomponent*)malloc((msize+csize) * sizeof(struct gcomponent));/*DIAGONALS OF GLOBAL MATRIX.*/
	for (i = 0; i < (msize+csize); i++)
	{
		(gmtx + i)->down = NULL;
		(gmtx2 + i)->down = NULL;
	}

	/*GLOBAL VECTOR*/
	gvct = (double *)malloc((msize+csize) * sizeof(double));/*INCREMENTAL GLOBAL VECTOR.*/


	/*FORCE VECTOR INITIALIZATION*/
	funbalance = (double*)malloc(msize * sizeof(double));        /*UNBALANCED FORCE VECTOR.*/
	freaction = (double*)malloc(msize * sizeof(double));         /*INERTIAL FORCE VECTOR.*/
	fexternal = (double*)malloc(msize * sizeof(double));         /*BASE EXTERNAL FORCE VECTOR.*/
	finternal = (double*)malloc(msize * sizeof(double));         /*INTERNAL FORCE VECTOR.*/
	finertial = (double*)malloc(msize * sizeof(double));         /*INERTIAL FORCE VECTOR.*/
	fdamping = (double*)malloc(msize * sizeof(double));         /*DAMPING FORCE VECTOR.*/
	fbaseload = (double*)malloc(msize * sizeof(double));           /*BASE LOAD VECTOR.*/
	fdeadload = (double*)malloc(msize * sizeof(double));           /*DEAD LOAD VECTOR.*/
	fpressure = (double*)malloc(msize * sizeof(double));         /*BASE PRESSURE FORCE VECTOR.*/
	fgivendisp = (double*)malloc(msize * sizeof(double));
	fconstraint = (double*)malloc(msize * sizeof(double));

	/*POSITION VECTOR INITIALIZATION*/
	/*IN SPATIAL FORM*/
	lastddisp = (double*)malloc(msize * sizeof(double));	/*LATEST LAP*/
	lapddisp  = (double*)malloc(msize * sizeof(double));	/*INCREMENT IN THE LAP*/
	givendisp = (double*)malloc(msize * sizeof(double));


	/*VELOSITY & ACCELERATION VECTOR INITIALIZATION*/
	//ud : VELOSITY, udd : ACCELERATION
	//_m : ROTATIONAL DOFs ARE REPRESENTED IN SPATIAL & MATERIAL FORM
	lastlastudd_m = (double*)malloc(msize * sizeof(double));
	lastud_m = (double*)malloc(msize * sizeof(double));     /*LATEST LAP IN MATERIAL*/
	lastudd_m = (double*)malloc(msize * sizeof(double));
	ud_m = (double*)malloc(msize * sizeof(double));  		/*LATEST ITERATION IN MATERIAL*/
	udd_m = (double*)malloc(msize * sizeof(double));

	/*
	if(inputfilename!=NULL && targetlap!=NULL)
	{
	  inputdsp(af, inputfilename, targetlap, 1);
	  inputplst(af, inputfilename, targetlap, 1);
	}
	*/

	initialform(ninit, iform, nnode);           /*ASSEMBLAGE FORMATION.*/
	initialform(nodes, ddisp, nnode);           /*ASSEMBLAGE FORMATION.*/
	initialelem(elems, melem, nelem);             /*ASSEMBLAGE ELEMENTS.*/
	initialshell(shells, mshell, nshell);         /*ASSEMBLAGE ELEMENTS.*/

	assemconf(confs,fdeadload,1.0,nnode);
	assemgivend(confs,givendisp,1.0,nnode);

	setviewpoint((wdraw.childs+0)->hwnd,arc,&((wdraw.childs+1)->vparam));
	setviewparam((wmenu.childs+2)->hwnd,(wdraw.childs+1)->vparam);

	GetAsyncKeyState(VK_LBUTTON);                   /*CLEAR KEY LEFT.*/
	GetAsyncKeyState(VK_RBUTTON);                  /*CLEAR KEY RIGHT.*/
	if(globaldrawflag==1)
	{
	  drawglobalaxis((wdraw.childs+1)->hdcC,(wdraw.childs+1)->vparam,0,0,255);  /*DRAW GLOBAL AXIS.*/
	}


	if(af->nlaps != NULL)
	{
		nlap = af->nlaps;/*INITIAL FROM ARCLMFRAME*/
	}
	beginlap = nlap;/*BEGINING OF THIS ANALYSIS*/
	setincrement((wmenu.childs+2)->hwnd,laps,nlap,maxiteration,iteration);

	/*INITIAL*/
	if(nlap==1)/*INITIAL FROM ANALYSISDATA*/
	{
		loadfactor = initialloadfactor + loadfactormap(time);
	}
	else if(STATDYNAFLAG == 1)/*INITIAL FROM ARCLMFRAME*/
	{
		initialloadfactor = af->loadfactor;
		loadfactor = initialloadfactor + loadfactormap(time);
	}
	else
	{
		loadfactor = loadfactormap(time);
	}
	lastloadfactor = loadfactor;

	for(i = 0; i < msize; i++)
	{
		*(lapddisp  + i) = 0.0;/*lapddisp : INCREMENTAL TRANSITION & ROTATION IN THIS LAP.*/
		*(lastddisp + i) = *(ddisp + i);

		*(lastlastudd_m + i) = 0.0;

		*(lastud_m + i) = 0.0;
		*(lastudd_m + i) = 0.0;

		*(ud_m + i) = 0.0;
		*(udd_m + i) = 0.0;
	}
	for (i = 0; i < csize; i++)
	{
		*(lastlambda + i) = *(lambda + i);
	}

	for(i = 0; i < msize; i++)
	{
		*(finertial  + i) = 0.0;
		*(finternal  + i) = 0.0;
		*(fexternal  + i) = 0.0;
		*(funbalance + i) = 0.0;
		*(freaction  + i) = 0.0;
		*(fbaseload  + i) = 0.0;
		*(fpressure  + i) = 0.0;
		*(fconstraint + i) = 0.0;
	}
	for(i = 0; i < csize; i++)
	{
		*(constraintvct + i) = 0.0;
	}

	volume = 0.0;
	assemshellvolume(shells, nshell, ddisp, &volume);

	/*POST PROCESS*/
	/*elemstress_DYNA(elems, melem, nelem, constraintmain,
					iform, lastddisp, ddisp,
					lastudd_m, udd_m,
					lastud_m, ud_m,
					finertial, fdamping, finternal, fexternal,
					alpham, alphaf, xi);*/
	shellstress_DYNA(shells, mshell, nshell, constraintmain,
					 iform, lastddisp, ddisp,
					 lastudd_m, udd_m,
					 lastud_m, ud_m,
					 finertial, fdamping, finternal, fpressure,
					 alpham, alphaf, xi);
	/*constraintstress_DYNA(constraints, nconstraint, constraintmain,
						  iform, lastddisp, ddisp, lastlambda, lambda,
						  fconstraint,constraintvct,
						  alphaf);*/


	strainenergy(af, &Wet, &Wpt);
	kineticenergy(af, ud_m, &Wkt);


	for (i = 0; i < msize; i++)
	{
		*(fbaseload + i) = *(fdeadload + i)+*(fpressure + i);
		*(fexternal + i) = (alphaf * lastloadfactor + (1-alphaf) * loadfactor) * *(fbaseload + i);
		*(funbalance + i) = *(fexternal + i) - *(finertial + i) - *(finternal + i) - *(fconstraint + i);
		if ((confs + i)->iconf == 1)
		{
			*(freaction + i) = *(funbalance + i);
			*(funbalance + i) = 0.0;
			*(fbaseload + i) = 0.0;
		}
		*(gvct + i) = *(funbalance + i);
	}
	for (i = 0; i < csize; i++)
	{
		*(gvct + msize + i) = *(constraintvct + i);
	}


	/*OUTPUT(INITIAL)*/
	sprintf(string, "LAP: %5ld / %5ld ITERATION: %5ld\n", nlap, laps, iteration);

	fprintf(fdsp, string);
	fprintf(fvel, string);
	fprintf(facc, string);
	outputdisp(ddisp, fdsp, nnode, nodes);
	outputdisp(ud_m, fvel, nnode, nodes);
	outputdisp(udd_m, facc, nnode, nodes);

	fprintf(finr, string);
	fprintf(finf, string);
	fprintf(fexf, string);
	fprintf(fubf, string);
	fprintf(frct, string);
	outputdisp(finertial, finr, nnode, nodes);
	outputdisp(finternal, finf, nnode, nodes);
	outputdisp(fexternal, fexf, nnode, nodes);
	outputdisp(funbalance, fubf, nnode, nodes);
	outputdisp(freaction, frct, nnode, nodes);

	fprintf(fstr, string);
	fprintf(fene, string);
	fprintf(fplst,string);
	for(i = 0; i < nshell; i++)
	{
		//fprintf(fene, "%5ld %e %e %e\n", (shells+i)->code, (mshell+i)->SEp, (mshell+i)->SEb, (mshell+i)->SE);
		fprintf(fstr, "%5ld %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n", (shells+i)->code,
		((shells+i)->gp[0]). stress[0],((shells+i)->gp[0]). stress[1],((shells+i)->gp[0]). stress[2],((shells+i)->gp[0]). stress[3],((shells+i)->gp[0]). stress[4],((shells+i)->gp[0]). stress[5],
		((shells+i)->gp[0]).estrain[0],((shells+i)->gp[0]).estrain[1],((shells+i)->gp[0]).estrain[2],((shells+i)->gp[0]).estrain[3],((shells+i)->gp[0]).estrain[4],((shells+i)->gp[0]).estrain[5],
		((shells+i)->gp[0]).pstrain[0],((shells+i)->gp[0]).pstrain[1],((shells+i)->gp[0]).pstrain[2],((shells+i)->gp[0]).pstrain[3],((shells+i)->gp[0]).pstrain[4],((shells+i)->gp[0]).pstrain[5],
		((shells+i)->gp[0]).alpha,((shells+i)->gp[0]).f[0],((shells+i)->gp[0]).f[1]
		);
	}
	for(i = 0; i < nshell; i++)
	{
		fprintf(fplst, "%5ld %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n", (shells+i)->code,
		((shells+i)->gp[0]).pstrain[0],((shells+i)->gp[0]).pstrain[1],((shells+i)->gp[0]).pstrain[2],((shells+i)->gp[0]).pstrain[3],((shells+i)->gp[0]).pstrain[4],((shells+i)->gp[0]).pstrain[5],((shells+i)->gp[0]).alpha,
		((shells+i)->gp[1]).pstrain[0],((shells+i)->gp[1]).pstrain[1],((shells+i)->gp[1]).pstrain[2],((shells+i)->gp[1]).pstrain[3],((shells+i)->gp[1]).pstrain[4],((shells+i)->gp[1]).pstrain[5],((shells+i)->gp[1]).alpha,
		((shells+i)->gp[2]).pstrain[0],((shells+i)->gp[2]).pstrain[1],((shells+i)->gp[2]).pstrain[2],((shells+i)->gp[2]).pstrain[3],((shells+i)->gp[2]).pstrain[4],((shells+i)->gp[2]).pstrain[5],((shells+i)->gp[2]).alpha,
		((shells+i)->gp[3]).pstrain[0],((shells+i)->gp[3]).pstrain[1],((shells+i)->gp[3]).pstrain[2],((shells+i)->gp[3]).pstrain[3],((shells+i)->gp[3]).pstrain[4],((shells+i)->gp[3]).pstrain[5],((shells+i)->gp[3]).alpha
		);
	}
	fprintf(fene, "%e %e %e %e\n", Wet, Wpt, Wkt, Wot);

	clearwindow(*(wdraw.childs+1));
	drawarclmframe((wdraw.childs+1)->hdcC,(wdraw.childs+1)->vparam,*af,0,ONSCREEN);
	overlayhdc(*(wdraw.childs + 1), SRCPAINT);                  /*UPDATE DISPLAY.*/


	while (nlap <= laps && ENDFLAG == 0)
	{
		/*MATRIX INITIALIZE*/
		std::vector<Triplet> Ktriplet;/*FOR DYNAMIC TANGENTIAL STIFFNESS*/
		std::vector<Triplet> Ktriplet2;/*FOR STATIC TANGENTIAL STIFFNESS*/

		SparseMatrix Kglobal((msize+csize), (msize+csize));
		SparseMatrix Kglobal2((msize+csize), (msize+csize));
		/*EXECUTE BY WIN64. AMD ORDERING IS AVAILABLE ONLY BY WIN64*/
		//Eigen::SimplicialLDLT<SparseMatrix,Eigen::Lower,Eigen::NaturalOrdering<int>> solver;
		//Eigen::SimplicialLDLT<SparseMatrix> solver;
		Eigen::SparseLU<SparseMatrix> solver;


		//Eigen::BiCGSTAB<SparseMatrix> solver;
		//solver.setTolerance(1e-9); // 許容誤差の設定
		//solver.setMaxIterations(10000); // 最大反復回数

		for (i = 1; i <= (msize+csize); i++)/*FOR DYNAMIC TANGENTIAL STIFFNESS*/
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
		for (i = 1; i <= (msize+csize); i++)/*FOR STATIC TANGENTIAL STIFFNESS*/
		{
			g = (gmtx2 + (i - 1))->down;
			while (g != NULL)
			{
				p = g;
				g = g->down;
				free(p);
			}
			ginit.m = (unsigned short int)i;
			*(gmtx2 + (i - 1)) = ginit;
		}
		//comps = msize;

		/*MATRIX ASSEMBLAGE*/
		if(USINGEIGENFLAG==1)
		{
		  assemshellEigen_DYNA(shells, mshell, nshell, constraintmain, confs,
							   Ktriplet,Ktriplet2,
							   iform, lastddisp, ddisp,
							   ud_m, udd_m,
							   alpham, alphaf, xi, beta, ddt);
		  //assemconstraintEigen_DYNA;
		  /*GLOBAL MATRIX USING EIGEN*/
		  modifytriplet(Ktriplet,confs,msize);
		  modifytriplet(Ktriplet2,confs,msize);

		  Kglobal.reserve((msize+csize)*(msize+csize));
		  Kglobal.setFromTriplets(Ktriplet.begin(), Ktriplet.end());//SPARSE MATRIX FROM TRIPLET
		  Kglobal2.reserve((msize+csize)*(msize+csize));
		  Kglobal2.setFromTriplets(Ktriplet2.begin(), Ktriplet2.end());//SPARSE MATRIX FROM TRIPLET
		}
		else
		{
		  assemshell_DYNA(shells, mshell, nshell, constraintmain,
						  gmtx,gmtx2,
						  iform, lastddisp, ddisp,
						  ud_m, udd_m,
						  alpham, alphaf, xi, beta, ddt);

		  assemconstraint_DYNA(constraints, nconstraint, constraintmain,
							   gmtx,
							   iform, lastddisp, ddisp, lastlambda, lambda,
							   alphaf, msize);
		}

		/*SOLVE [K]*/
		if(USINGEIGENFLAG==1)
		{
			solver.compute(Kglobal);

			/*
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
			*/

			if (solver.info() != Eigen::Success)sign=-1;
			/*DECOMPOSITION FAILED*/
			if (sign < 0.0)
			{
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
                fclose(fplst);
				return 1;
			}

			sprintf(string, "LAP: %4d ITER: %2d {LOAD}= %5.8f {RESD}= %1.5e {DET}= %5.5f {SIGN}= %2.0f {BCL}= %1d {EPS}= %1.5e {V}= %5.5f {TIME}= %5.5f\n",
					nlap, iteration, loadfactor,residual, Wkt, 0, 0, ddt, volume, time);
			fprintf(ffig, "%s", string);
			errormessage(string);

			Vector P = Vector::Zero(msize + csize);
			for(i=0;i<msize + csize;i++)P(i)=*(gvct+i);
			Vector U = solver.solve(P);

			/*solver.iteration();
			solver.error();*/


			for(i=0;i<msize + csize;i++)*(gvct+i)=U(i);
			if (solver.info() != Eigen::Success)return 1;
		}
		else
		{
			/*CROUT LU DECOMPOSITION.*/
			nline = croutluII(gmtx, confs, msize, csize, &determinant, &sign, gcomp1);
			nline2 = croutluII(gmtx2, confs, msize, csize, &determinant2, &sign2, gcomp1);

			/*DECOMPOSITION FAILED*/
			if (sign < 0.0)
			{

				sprintf(string, "INSTABLE TERMINATION AT NODE %ld.\n",nline);
				fprintf(ffig, "%s", string);
				errormessage(string);

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
                fclose(fplst);
				return 1;
			}
			sprintf(string, "LAP: %4d ITER: %2d {LOAD}= %5.8f {RESD}= %1.5e {DET}= %5.5f {SIGN}= %2.0f {BCL}= %1d {EPS}= %1.5e {V}= %5.5f {TIME}= %5.5f\n",
					nlap, iteration, loadfactor,residual, determinant, sign2, 0, ddt, volume, time);
			fprintf(ffig, "%s", string);
			errormessage(string);


			nline = forwardbackwardII(gmtx, gvct, confs, msize, csize, gcomp1);
        }

		for (ii = 0; ii < msize; ii++)
		{
			*(gvct + ii) = *(gvct + *(constraintmain + ii));
		}

		/*COORD UPDATE*/
		updaterotation(lapddisp, gvct, nnode);
		updaterotation(ddisp, gvct, nnode);
		for(ii = 0; ii < nnode;ii++)
		{
			inputnode(ddisp,nodes+ii);
		}

		/*CONSTRAINT UPDATE*/
		for (ii = 0; ii < csize; ii++)
		{
		  *(lambda+ii)+=*(gvct+msize+ii);
		}

		/*ACC & VEL UPDATE*/
		lapddisp_m = pullback(lastddisp,lapddisp,nnode);
		for(ii = 0; ii < msize; ii++)
		{
			*(ud_m + ii)  =  (gamma / (beta * ddt)) **(lapddisp_m + ii)
						   - (gamma / beta - 1.0) **(lastud_m + ii)
						   - (gamma / (2.0 * beta) - 1.0) * ddt **(lastudd_m + ii);
			*(udd_m + ii) =  (1.0 / (beta * ddt * ddt)) **(lapddisp_m + ii)
						   - (1.0 / (beta * ddt)) **(lastud_m + ii)
						   - (1.0 / (2.0 * beta) - 1.0) **(lastudd_m + ii);
		}

		/*FORCE INITIALIZATION*/
		for(ii = 0; ii < msize; ii++)
		{
			*(finertial  + ii) = 0.0;
			*(finternal  + ii) = 0.0;
			*(fexternal  + ii) = 0.0;
			*(funbalance + ii) = 0.0;
			*(freaction  + ii) = 0.0;
			*(fbaseload  + ii) = 0.0;
			*(fpressure  + ii) = 0.0;
			*(fconstraint + ii) = 0.0;
		}
		for (i = 0; i < csize; i++)
		{
			*(constraintvct  + i) = 0.0;
		}

		volume = 0.0;
		assemshellvolume(shells, nshell, ddisp, &volume);

		/*POST PROCESS*/
		/*elemstress_DYNA(elems, melem, nelem, constraintmain,
						iform, lastddisp, ddisp,
						lastudd_m, udd_m,
						lastud_m, ud_m,
						finertial, fdamping, finternal, fexternal,
						alpham, alphaf, xi);*/
		shellstress_DYNA(shells, mshell, nshell, constraintmain,
						 iform, lastddisp, ddisp,
						 lastudd_m, udd_m,
						 lastud_m, ud_m,
						 finertial, fdamping, finternal, fpressure,
					     alpham, alphaf, xi);
		/*constraintstress_DYNA(constraints, nconstraint, constraintmain,
							  iform, lastddisp, ddisp, lastlambda, lambda,
							  fconstraint,constraintvct,
							  alphaf);*/

		strainenergy(af, &Wet, &Wpt);
		kineticenergy(af, ud_m, &Wkt);

		for (i = 0; i < msize; i++)
		{
			*(fbaseload + i) = *(fdeadload + i)+*(fpressure + i);
			*(fexternal + i) = (alphaf * lastloadfactor + (1-alphaf) * loadfactor) * *(fbaseload + i);
			*(funbalance + i) = *(fexternal + i) - *(finertial + i) - *(finternal + i) - *(fconstraint + i);
			if ((confs + i)->iconf == 1)
			{
				*(freaction + i) = *(funbalance + i);
				*(funbalance + i) = 0.0;
				*(fbaseload + i) = 0.0;
			}
			*(gvct + i) = *(funbalance + i);
		}
		for (i = 0; i < csize; i++)	 /*GLOBAL VECTOR INITIALIZATION.*/
		{
			*(gvct + msize + i) = *(constraintvct  + i);
		}


		/*CONVERGENCE JUDGEMENT.*/
		residual = vectorlength(funbalance,msize);
		constraintresidual = vectorlength(constraintvct,csize);
		gvctlen = vectorlength(gvct,msize);

		if (!isfinite(sign) || sign > 20 || !isfinite(residual) || residual > 1e+10)
		{
			DIVFLAG = 1;
			//ENDFLAG = 1;
			sprintf(string,"DIVERGENCE DITECTED(SIGN = %f). ANALYSIS TERMINATED.\n", sign);
			errormessage(string);
		}
		else
		{
			DIVFLAG = 0;
		}

		//if ((residual < tolerance || iteration > maxiteration-1)) && iteration!=1)
		if ((gvctlen < tolerance || iteration > maxiteration-1) && iteration!=1)
		{
			nlap += 1;
			iteration = 0;
			if(RELAXATION == 1 && lastWkt > Wkt && time>0.01)
			{
				PEAKFLAG = 1;
				ddt = initialddt;
				if(Wkt<0.01)
				{
					ENDFLAG = 1;
					sprintf(string,"KINEMATIC ENERGY CONVERGED LAP: %4d ITER: %2d\n", nlap, iteration);
				}
				time += ddt;
			}
			else if(DIVFLAG == 1)
			{
				ddt *= 0.1;
			}
			else
			{
				ddt = timestepcontrol(ddt, lapddisp_m, udd_m, lastudd_m, lastlastudd_m, beta, msize);
				time += ddt;
			}

			clearwindow(*(wdraw.childs+1));
			drawarclmframe((wdraw.childs+1)->hdcC,(wdraw.childs+1)->vparam,*af,0,ONSCREEN);
			overlayhdc(*(wdraw.childs + 1), SRCPAINT);                  /*UPDATE DISPLAY.*/
		}
		iteration++;
		setincrement((wmenu.childs+2)->hwnd,laps,nlap,maxiteration,iteration);

		free(lapddisp_m);

		if (((outputmode == 0 && iteration == 1) || outputmode == 1) && DIVFLAG == 0)
		{
			sprintf(string, "LAP: %5ld / %5ld ITERATION: %5ld\n", nlap, laps, iteration);

			fprintf(fdsp, string);
			fprintf(fvel, string);
			fprintf(facc, string);
			outputdisp(ddisp, fdsp, nnode, nodes);
			outputdisp(ud_m, fvel, nnode, nodes);
			outputdisp(udd_m, facc, nnode, nodes);

			fprintf(finr, string);
			fprintf(finf, string);
			fprintf(fexf, string);
			fprintf(fubf, string);
			fprintf(frct, string);
			outputdisp(finertial, finr, nnode, nodes);
			outputdisp(finternal, finf, nnode, nodes);
			outputdisp(fexternal, fexf, nnode, nodes);
			outputdisp(funbalance, fubf, nnode, nodes);
			outputdisp(freaction, frct, nnode, nodes);

			fprintf(fstr, string);
			fprintf(fene, string);
			fprintf(fplst,string);
			for(i = 0; i < nshell; i++)
			{
				//fprintf(fene, "%5ld %e %e %e\n", (shells+i)->code, (mshell+i)->SEp, (mshell+i)->SEb, (mshell+i)->SE);
				fprintf(fstr, "%5ld %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n", (shells+i)->code,
				((shells+i)->gp[0]). stress[0],((shells+i)->gp[0]). stress[1],((shells+i)->gp[0]). stress[2],((shells+i)->gp[0]). stress[3],((shells+i)->gp[0]). stress[4],((shells+i)->gp[0]). stress[5],
				((shells+i)->gp[0]).estrain[0],((shells+i)->gp[0]).estrain[1],((shells+i)->gp[0]).estrain[2],((shells+i)->gp[0]).estrain[3],((shells+i)->gp[0]).estrain[4],((shells+i)->gp[0]).estrain[5],
				((shells+i)->gp[0]).pstrain[0],((shells+i)->gp[0]).pstrain[1],((shells+i)->gp[0]).pstrain[2],((shells+i)->gp[0]).pstrain[3],((shells+i)->gp[0]).pstrain[4],((shells+i)->gp[0]).pstrain[5],
				((shells+i)->gp[0]).alpha,((shells+i)->gp[0]).f[0],((shells+i)->gp[0]).f[1]
				);
			}
			for(i = 0; i < nshell; i++)
			{
			    fprintf(fplst, "%5ld %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n", (shells+i)->code,
				((shells+i)->gp[0]).pstrain[0],((shells+i)->gp[0]).pstrain[1],((shells+i)->gp[0]).pstrain[2],((shells+i)->gp[0]).pstrain[3],((shells+i)->gp[0]).pstrain[4],((shells+i)->gp[0]).pstrain[5],((shells+i)->gp[0]).alpha,
				((shells+i)->gp[1]).pstrain[0],((shells+i)->gp[1]).pstrain[1],((shells+i)->gp[1]).pstrain[2],((shells+i)->gp[1]).pstrain[3],((shells+i)->gp[1]).pstrain[4],((shells+i)->gp[1]).pstrain[5],((shells+i)->gp[1]).alpha,
				((shells+i)->gp[2]).pstrain[0],((shells+i)->gp[2]).pstrain[1],((shells+i)->gp[2]).pstrain[2],((shells+i)->gp[2]).pstrain[3],((shells+i)->gp[2]).pstrain[4],((shells+i)->gp[2]).pstrain[5],((shells+i)->gp[2]).alpha,
				((shells+i)->gp[3]).pstrain[0],((shells+i)->gp[3]).pstrain[1],((shells+i)->gp[3]).pstrain[2],((shells+i)->gp[3]).pstrain[3],((shells+i)->gp[3]).pstrain[4],((shells+i)->gp[3]).pstrain[5],((shells+i)->gp[3]).alpha
				);
			}
			fprintf(fene, "%e %e %e %e\n", Wet, Wpt, Wkt, Wot);
		}



		if(iteration==1)
		{
			if(DIVFLAG == 1)/*TERMINATION AT LAST LAP*/
			{
				loadfactor = lastloadfactor;
				for (ii = 0; ii < msize; ii++)
				{
					*(lapddisp + ii) = 0.0;
					*(ddisp + ii) = *(lastddisp + ii);

					*(ud_m + ii)  =- (gamma / beta - 1.0) **(lastud_m + ii)
								   - (gamma / (2.0 * beta) - 1.0) * ddt **(lastudd_m + ii);
					*(udd_m + ii) =- (1.0 / (beta * ddt)) **(lastud_m + ii)
								   - (1.0 / (2.0 * beta) - 1.0) **(lastudd_m + ii);
				}
				for(ii = 0; ii < nnode;ii++)
				{
					inputnode(ddisp,nodes+ii);
				}
				for (ii = 0; ii < csize; ii++)
				{
					*(lambda+ii) = *(lastlambda+ii) ;
				}
				Wkt = lastWkt;
			}
			else
			{
				initialelem(elems,melem,nelem);    /*melem (LAST LAP DATA) UPDATE.*/
				initialshell(shells,mshell,nshell);/*mshell(LAST LAP DATA) UPDATE.*/

				lastloadfactor = loadfactor;
				if(STATDYNAFLAG == 1)/*INITIAL FROM ARCLMFRAME*/
				{
					loadfactor = initialloadfactor + loadfactormap(time);
				}
				else
				{
					loadfactor = loadfactormap(time);
				}

				if(PEAKFLAG == 1)
				{
					for (ii = 0; ii < msize; ii++)
					{
						*(lapddisp  + ii) = 0.0;/*lapddisp : INCREMENTAL TRANSITION & ROTATION IN THIS LAP.*/
						*(lastddisp + ii) = *(ddisp + ii);

						*(lastlastudd_m  + ii) = 0.0;

						*(lastud_m  + ii) = 0.0;
						*(lastudd_m + ii) = 0.0;

						*(ud_m + ii)  = 0.0;
						*(udd_m + ii) = 0.0;
					}
					for (ii = 0; ii < csize; ii++)
					{
						*(lastlambda + ii) = *(lambda + ii);
					}
					PEAKFLAG = 0;
					lastWkt = 0.0;
				}
				else
				{
					for (ii = 0; ii < msize; ii++)
					{
						*(lapddisp  + ii) = 0.0;/*lapddisp : INCREMENTAL TRANSITION & ROTATION IN THIS LAP.*/
						*(lastddisp + ii) = *(ddisp + ii);

						*(lastlastudd_m  + ii) = *(lastud_m  + ii);

						*(lastud_m  + ii) = *(ud_m  + ii);
						*(lastudd_m + ii) = *(udd_m + ii);

						*(ud_m + ii)  =- (gamma / beta - 1.0) **(lastud_m + ii)
									   - (gamma / (2.0 * beta) - 1.0) * ddt **(lastudd_m + ii);
						*(udd_m + ii) =- (1.0 / (beta * ddt)) **(lastud_m + ii)
									   - (1.0 / (2.0 * beta) - 1.0) **(lastudd_m + ii);
					}
					for (ii = 0; ii < csize; ii++)
					{
						*(lastlambda + ii) = *(lambda + ii);
					}
					lastWkt = Wkt;
				}
			}

			for(i = 0; i < msize; i++)
			{
				*(finertial  + i) = 0.0;
				*(finternal  + i) = 0.0;
				*(fexternal  + i) = 0.0;
				*(funbalance + i) = 0.0;
				*(freaction  + i) = 0.0;
				*(fbaseload  + i) = 0.0;
				*(fpressure  + i) = 0.0;
				*(fconstraint + i) = 0.0;
			}
			for(i = 0; i < csize; i++)
			{
				*(constraintvct + i) = 0.0;
			}


			//volume = 0.0;
			//assemshellvolume(shells, nshell, ddisp, &volume);

			/*STRESS*/
			/*elemstress_DYNA(elems, melem, nelem, constraintmain,
							  iform, lastddisp, ddisp,
							  lastudd_m, udd_m,
							  lastud_m, ud_m,
							  finertial, fdamping, finternal, fexternal,
							  alpham, alphaf, xi);*/
			shellstress_DYNA(shells, mshell, nshell, constraintmain,
							 iform, lastddisp, ddisp,
							 lastudd_m, udd_m,
							 lastud_m, ud_m,
							 finertial, fdamping, finternal, fpressure,
							 alpham, alphaf, xi);
			/*constraintstress_DYNA(constraints, nconstraint, constraintmain,
								  iform, lastddisp, ddisp, lastlambda, lambda,
								  fconstraint,constraintvct,
								  alphaf);*/


			//strainenergy(af, &Wet, &Wpt);
			//kineticenergy(af, ud_m, &Wkt);

			for (i = 0; i < msize; i++)
			{
				*(fbaseload + i) = *(fdeadload + i)+*(fpressure + i);
				*(fexternal + i) = (alphaf * lastloadfactor + (1-alphaf) * loadfactor) * *(fbaseload + i);
				*(funbalance + i) = *(fexternal + i) - *(finertial + i) - *(finternal + i) - *(fconstraint + i);
				if ((confs + i)->iconf == 1)
				{
					*(freaction + i) = *(funbalance + i);
					*(funbalance + i) = 0.0;
					*(fbaseload + i) = 0.0;
				}
				*(gvct + i) = *(funbalance + i);
			}
			for (i = 0; i < csize; i++)
			{
				*(gvct + msize + i) = *(constraintvct + i);
			}

		}




		//MESSAGE FOR UPDATE UI
		//AVOID FREEZE ON LONG RUNNING TASK
		MSG msg;
		if(PeekMessage(&msg, NULL, 0, 0, PM_REMOVE))
		{
			TranslateMessage(&msg);
			DispatchMessage(&msg);
		}
#if 0
		while(GetAsyncKeyState(VK_LBUTTON))  /*LEFT CLICK TO CONTINUE.*/
		{
		  if(GetAsyncKeyState(VK_RBUTTON))      /*RIGHT CLICK TO ABORT.*/
		  {
			gfree(gmtx, nnode);  /*FREE GLOBAL MATRIX.*/
			free(gvct);
			free(funbalance);
			free(finternal);
			free(fexternal);
			free(fpressure);

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
            fclose(fplst);

			laptime("\0",t0);
			return 1;
		  }

		  //t2=clock();
		  //if((t2-t1)/CLK_TCK>=WAIT) break;               /*CONTINUE AFTER WAITING.*/
		}
#endif
	}

	af->nlaps = nlap;
	af->loadfactor = loadfactor;

	clearwindow(*(wdraw.childs+1));
	drawarclmframe((wdraw.childs+1)->hdcC,(wdraw.childs+1)->vparam,*af,0,ONSCREEN);
	overlayhdc(*(wdraw.childs + 1), SRCPAINT);                  /*UPDATE DISPLAY.*/


	fclose(fonl);
	fclose(fdsp);
	fclose(fexf);
	fclose(finf);
	fclose(fubf);
	fclose(frct);
	fclose(finr);
	fclose(fstr);
	fclose(fene);
	fclose(ffig);
	fclose(fbcl);
	fclose(feig);
	fclose(fout);
    fclose(fvel);
	fclose(facc);
	fclose(fplst);

	gfree(gmtx, nnode);  	/*FREE GLOBAL MATRIX.*/
	gfree(gmtx2, nnode);

	free(gvct);

	free(lastlambda);
	free(lastddisp);
	free(lapddisp);
	free(givendisp);

	free(lastlastudd_m);
	free(lastud_m);
	free(lastudd_m);
	free(ud_m);
	free(udd_m);

	free(funbalance);
	free(freaction);
	free(fexternal);
	free(finternal);
	free(finertial);
	free(fdamping);
	free(fbaseload);
	free(fdeadload);
	free(fpressure);
	free(fgivendisp);
	free(fconstraint);

	free(iform);
	free(ddisp);
	free(melem);
	free(mshell);
	free(constraintvct);

	errormessage(" ");
	errormessage("COMPLETED.");

	memory1=availablephysicalmemoryEx("REMAIN:");
	sprintf(string,"CONSUMPTION:%ld[BYTES]",(memory0-memory1));
	errormessage(string);

	if(!isfinite(sign) || sign > 20 || !isfinite(residual) || residual > 1e+10)return 1;
	else return 0;
}




double timestepcontrol(double ddt, double* lapddisp_m, double* udd_m, double* lastudd_m, double* lastlastudd_m, double beta, int msize)
{
	int i;
	char string[100];
	double eta_upper  = 0.006;
	double eta_target = 0.005;
	double eta_lower  = 0.004;
	double eta, errornorm, dispnorm;
	double* error, *disp;

	error = (double *)malloc(msize * sizeof(double));
	disp = (double *)malloc(msize * sizeof(double));
	for (i = 0; i < msize; i++)
	{
	  if(i%6<3)
	  {
#if 0
		  /*Time Step Control Algorithm by Schweizerhof & Riccius*/
		  *(error + i) = ddt * ddt * ( (beta-(1.0/8.0))**(udd_m + i) + ((1.0/12.0)-beta)**(lastudd_m + i) + (1.0/24.0)**(lastlastudd_m + i));
#endif
#if 1
		  /*Time Step Control Algorithm by Zienkiewics & Xie*/
		  *(error + i) = ddt * ddt * (beta-(1.0/6.0)) * (*(udd_m + i) - *(lastudd_m + i));
#endif
		  *(disp + i) =  *(lapddisp_m + i);
	  }
	  else
	  {
		  *(error + i) =  0.0;
		  *(disp + i) =  0.0;
	  }
	}

	errornorm = vectorlength(error, msize);
	dispnorm  = vectorlength(lapddisp_m, msize);

	eta = errornorm / dispnorm;

	//sprintf(string, "eta: %e\n",eta);
	//errormessage(string);

	if(eta > eta_upper || eta < eta_lower)
	{
	  ddt *= sqrt(eta_target/eta);
	}

	free(error);
	free(disp);


	return ddt;
}


double loadfactormap(double time)
{
	double loadfactor;

#if 0/*FOR SHELL*/

	  loadfactor = 0.0 * time;

#endif
#if 0/*FOR SHELL*/
	if(time<=0.2)
	{
	  loadfactor = 2.5e+8 * time;
	}
	else
	{
	  loadfactor = 5.0e+7;
	}
#endif
#if 1/*FOR HINGE*/
/*0.1Mpa = 100000Pa*/
	if(time<=100)
	{
	  loadfactor = 100* time;
	}
	else
	{
	  loadfactor = 10000;
	}
#endif
#if 0/*FOR CYLINDER*/
	if(time<=0.5)
	{
	  loadfactor = 10.0 * time;
	}
	else if(time<=1.0 && time>0.5)
	{
	  loadfactor = 10.0-10.0*time;
	}
	else
	{
	  loadfactor = 0.0;
	}
#endif
	return loadfactor;
}


