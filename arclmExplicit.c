/* ========================================================= */
/* PROGRAM GNSHN CR FOR OUTER SURFACE OF LUNAR MARS PROJECT  */
/* DYNAMIC ANALYSIS FOR LARGE DEFORMATION                    */
/* USING CR FORMULATION                                      */
/* CODED BY KEIICHI KUTOMI SINSE 2024.05.26                  */
/* ========================================================= */


int arclmExplicit(struct arclmframe* af)
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
	int solver = 1;
	int method = 0;
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

	if(inputfilename!=NULL && targetlap!=NULL)
	{
	  inputdsp(af, inputfilename, targetlap, 1);
	  inputplst(af, inputfilename, targetlap, 1);
	}

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



		/*SOLVE*/


		/*DYNAMIC RELAXATION.*/

		double leftleftWkt,leftWkt,Wkt,rightWkt;/*KINETIC ENERGY AT t-3*dt/2, t-dt/2, t, t+dt/2*/
		double q;

		rightud_m = (double*)malloc(msize * sizeof(double));
		leftud_m = (double*)malloc(msize * sizeof(double));
		for (i = 0; i < msize; i++)
		{
			*(rightud_m + i) = 0.0;
			*(leftud_m + i) = 0.0;
		}

		for (i = 0; i < msize; i++)
		{
			if (time == 0.0 || PEAKFLAG)
			{
				*(rightud_m + i) = *(funbalance + i) * ddt / (2.0 * *(mvct + i));
				*(leftud_m  + i) = -*(rightud_m  + i);
				leftleftWkt = 0.0;
				leftWkt = 0.0;

				PEAKFLAG = 0;
			}
			else
			{
				*(leftud_m  + i) = *(rightud_m  + i);
				*(rightud_m  + i) = (*(funbalance + i) * ddt / *(mvct + i)) + *(leftud_m  + i);
				//*(rightud_m  + i) = ( *(funbalance_m + i) + (*(mvct + i) / ddt - *(cvct + i) / 2.0) * *(leftud_m  + i) ) / ( *(mvct + i) / ddt + *(cvct + i) / 2.0 );
				leftleftWkt = leftWkt;
				leftWkt = rightWkt;
			}
			Wkt = 0.0;
			rightWkt = 0.0;

			*(ud_m  + i) = (*(rightud_m  + i) + *(leftud_m  + i)) / 2.0;
			*(udd_m  + i) = (*(rightud_m  + i) - *(leftud_m  + i)) / ddt;

			Wkt += 0.5 * *(mvct + i) * *(ud_m  + i) * *(ud_m  + i);
			rightWkt += 0.5 * *(mvct + i) * *(rightud_m  + i) * *(rightud_m  + i);
		}

		/*FOR DINAMIC RELAXATION*/
		if (Wkt < lastWkt)
		{
			q = (leftWkt - rightWkt) / (2 * leftWkt - rightWkt - leftleftWkt);
			for (i = 0; i < msize; i++)
			{
				*(gvct + i) = -q * *(leftud_m  + i) * ddt;
			}
			PEAKFLAG = 1;
		}
		else
		{
			for (i = 0; i < msize; i++)
			{
				*(gvct + i) = *(rightud_m  + i) * ddt;
			}
		}




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
				return 1;
			}
			sprintf(string, "LAP: %4d ITER: %2d {LOAD}= %5.8f {RESD}= %1.5e {DET}= %5.5f {SIGN}= %2.0f {BCL}= %1d {EPS}= %1.5e {V}= %5.5f {TIME}= %5.5f\n",
					nlap, iteration, loadfactor,residual, determinant, sign2, 0, ddt, volume, time);
			fprintf(ffig, "%s", string);
			errormessage(string);


			nline = forwardbackwardII(gmtx, gvct, confs, msize, csize, gcomp1);




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

		//if ((residual < tolerance || iteration > maxiteration-1)) && iteration!=1)
		if ((gvctlen < tolerance || iteration > maxiteration-1) && iteration!=1)
		{
			nlap += 1;;
			iteration = 0;
			if(RELAXATION == 1 && lastWkt > Wkt)
			{
				PEAKFLAG = 1;
			}
			time += ddt;

			clearwindow(*(wdraw.childs+1));
			drawarclmframe((wdraw.childs+1)->hdcC,(wdraw.childs+1)->vparam,*af,0,ONSCREEN);
			overlayhdc(*(wdraw.childs + 1), SRCPAINT);                  /*UPDATE DISPLAY.*/
		}
		iteration++;
		setincrement((wmenu.childs+2)->hwnd,laps,nlap,maxiteration,iteration);


		free(lapddisp_m);



		if ((outputmode == 0 && iteration == 1) || outputmode == 1)
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
			if (!isfinite(sign) || sign > 20 || !isfinite(residual) || residual > 1e+10)
			{
				ENDFLAG = 1;
				sprintf(string,"DIVERGENCE DITECTED(SIGN = %f). ANALYSIS TERMINATED.\n", sign);
				errormessage(string);
			}
			if(PEAKFLAG==1 && Wkt<0.01)
			{
				ENDFLAG = 1;
				sprintf(string,"KINEMATIC ENERGY CONVERGED LAP: %4d ITER: %2d\n", nlap, iteration);
				errormessage(string);
			}
		}

		if(iteration==1)
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


			if(0)
			{
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







