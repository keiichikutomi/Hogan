/* ========================================================= */
/* PROGRAM GNSHN CR FOR OUTER SURFACE OF LUNAR MARS PROJECT  */
/* DYNAMIC ANALYSIS FOR LARGE DEFORMATION                    */
/* USING CR FORMULATION                                      */
/* CODED BY KEIICHI KUTOMI SINSE 2024.05.26                  */
/* ========================================================= */

int arclmDynamic(struct arclmframe* af);
double timestepcontrol(double ddt, double* lapddisp_m, double* udd_m, double* lastudd_m, double* lastlastudd_m, double beta, int msize);

int arclmDynamic(struct arclmframe* af)
{
	DWORDLONG memory0,memory1;
	char dir[] = DIRECTORY;
	char string[400], fname[100];
	clock_t t0;
	FILE *fin, *fonl, * fdsp, * fexf, * finf, * fubf, * frct, * finr, * fstr, * fene, * ffig, * fbcl, * feig, *fout;         /*FILE 8 BYTES*/
	FILE *fvel,*facc;

	int i, ii, jj;

	/*MODEL SIZE*/
	int nnode, nelem, nshell, nsect, nconstraint;
	long int msize;

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

	/*GLOBAL MATRIX*/
	struct gcomponent ginit = {0,0,0.0,NULL};
	struct gcomponent* gmtx, * g, * p, * gcomp1;/*GLOBAL MATRIX*/
	double gg;
	double sign ,determinant;
	long int nline;

	int USINGEIGENFLAG = 0;

	/*GLOBAL VECTOR*/
	double* iform,* ddisp,* lastddisp,* lapddisp,* lapddisp_m;
	double* funbalance, * fexternal, * finternal, * freaction;
	double* fbaseload, * fdeadload, * fpressure,* finertial,* fdamping;
	double* funbalance_m;
	double* gvct,* gvct_s;
	double* udinit_m,* uddinit_m;
	double* lastlastudd_m;
	double* lastud_m,* lastudd_m;
	double* ud_m,* udd_m;

	double* rightud_m;/*VELOCITY AT t+dt/2 FOR EXPLICIT*/
	double* leftud_m;/*VELOCITY AT t-dt/2 FOR EXPLICIT*/
	double* mvct, * cvct;/*MASS & DAMPING VECTOR(DIAGONAL MATRIX)*/

	///FOR ITERATIVE CALCULATION///
	int nlap, laps;
	int iteration;
	int maxiteration = 20;
	double tolerance = 1.0E-4;
	double residual;
	double gvctlen;

	double ddt = 0.001;/*TIME INCREMENT[sec]*/
	double time = 0.0;/*TOTAL TIME[sec]*/

	double loadfactor = 0.0;
	double lastloadfactor = 0.0;
	//double lambda = 0.0;
	double volume = 0.0;

	double leftleftKE,leftKE,KE,rightKE;/*KINETIC ENERGY AT t-3*dt/2, t-dt/2, t, t+dt/2*/
	int PEAKFLAG = 0;
	double q;






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

	/*ANALYSIS TERMINATION*/
	int ENDFLAG = 0;
	int fnode=NULL,fnodedrc=NULL;
	double fnodemin, fnodemax;

	/*FOR READING ANALYSISDATA*/
	FILE *fdata;
	int nstr, pstr, readflag;
	char **data;
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
						else
						{
							pstr++;
						}
					}
				}
			}
		}
	}



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
		alpham = (2.0*rho-1)/(rho+1);  /*MID-POINT USED TO EVALUATE INERTIAL FORCE*/
		alphaf = rho/(rho+1);        /*MID-POINT USED TO EVALUATE INTERNAL FORCE*/
		xi = 0.0;   	 				/*NUMERICAL DISSIPATION COEFFICIENT*/
		beta = pow(1-alpham+alphaf,2)/4.0;
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
	sprintf(string,"alpham=%e alphaf=%e beta=%e gamma=%e ",alpham,alphaf,beta,gamma);
	errormessage(string);



	///INPUT FILE OPEN & OUTPUT FILE SETTING///
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
	iform = (double*)malloc(msize * sizeof(double));		/*INITIAL*/
	ddisp = (double*)malloc(msize * sizeof(double));		/*LATEST ITERATION*/
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

	/*GLOBAL MATRIX*/
	gmtx = (struct gcomponent*)malloc(msize * sizeof(struct gcomponent));/*DIAGONALS OF GLOBAL MATRIX.*/

	/*GLOBAL VECTOR*/
	gvct = (double *)malloc(msize * sizeof(double));/*INCREMENTAL GLOBAL VECTOR.*/

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


	/*POSITION VECTOR INITIALIZATION*/
	/*IN SPATIAL FORM*/
	lastddisp = (double*)malloc(msize * sizeof(double));	/*LATEST LAP*/
	lapddisp  = (double*)malloc(msize * sizeof(double));		/*INCREMENT IN THE LAP*/
	//lapddisp_m = (double*)malloc(msize * sizeof(double));		/*INCREMENT IN THE LAP*/


	/*VELOSITY & ACCELERATION VECTOR INITIALIZATION*/
	//ud : VELOSITY, udd : ACCELERATION
	//_m : ROTATIONAL DOFs ARE REPRESENTED IN SPATIAL & MATERIAL FORM
	udinit_m = (double*)malloc(msize * sizeof(double));  	   	/*NEWMARK INITIAL IN LAP*/
	uddinit_m = (double*)malloc(msize * sizeof(double));
	lastlastudd_m = (double*)malloc(msize * sizeof(double));
	lastud_m = (double*)malloc(msize * sizeof(double));     /*LATEST LAP IN MATERIAL*/
	lastudd_m = (double*)malloc(msize * sizeof(double));
	ud_m = (double*)malloc(msize * sizeof(double));  		/*LATEST ITERATION IN MATERIAL*/
	udd_m = (double*)malloc(msize * sizeof(double));

	for (i = 0; i < msize; i++)
	{
		*(lastlastudd_m + i) = 0.0;

		*(udinit_m + i) = 0.0;
		*(uddinit_m + i) = 0.0;

		*(lastud_m + i) = 0.0;
		*(lastudd_m + i) = 0.0;

		*(ud_m + i) = 0.0;
		*(udd_m + i) = 0.0;
	}




	initialelem(elems, melem, nelem);             /*ASSEMBLAGE ELEMENTS.*/
	initialshell(shells, mshell, nshell);         /*ASSEMBLAGE ELEMENTS.*/

	setviewpoint((wdraw.childs+0)->hwnd,arc,&((wdraw.childs+1)->vparam));
	setviewparam((wmenu.childs+2)->hwnd,(wdraw.childs+1)->vparam);

	GetAsyncKeyState(VK_LBUTTON);                   /*CLEAR KEY LEFT.*/
	GetAsyncKeyState(VK_RBUTTON);                  /*CLEAR KEY RIGHT.*/
	if(globaldrawflag==1)
	{
	  drawglobalaxis((wdraw.childs+1)->hdcC,(wdraw.childs+1)->vparam,0,0,255);  /*DRAW GLOBAL AXIS.*/
	}

	nlap = 0;
	iteration = 0;
	residual = 0.0;

	assemconf(confs,fdeadload,1.0,nnode);               /*GLOBAL VECTOR.*/
	//assemgivend(confs,givendisp,1.0,nnode);

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

		for (i = 0; i < msize; i++)	 /*GLOBAL VECTOR INITIALIZATION.*/
		{
			*(finertial  + i) = 0.0;
			*(finternal  + i) = 0.0;
			*(fexternal  + i) = 0.0;
			*(fpressure  + i) = 0.0;
			*(freaction  + i) = 0.0;
			*(funbalance + i) = 0.0;
		}
		comps = msize;

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
					   finertial, fdamping, finternal, fexternal,
					   alpham, alphaf, xi);

		for (i = 0; i < msize; i++)
		{
			*(fbaseload + i) = *(fdeadload + i)+*(fpressure + i);
			*(fexternal + i) = (alphaf * lastloadfactor + (1-alphaf) * loadfactor) * *(fbaseload + i);
			*(funbalance + i) = *(fexternal + i) - *(finertial + i) - *(finternal + i);
			if ((confs + i)->iconf == 1)
			{
				*(freaction + i) = *(funbalance + i);
				*(funbalance + i) = 0.0;
				*(fbaseload + i) = 0.0;
			}
			*(gvct + i) = *(funbalance + i);
		}

		/*CONVERGENCE JUDGEMENT.*/
		residual = vectorlength(funbalance,msize);
		//gvctlen = vectorlength(gvct,msize);

		if (residual < tolerance || iteration > maxiteration-1)
		{
		  nlap++;
		  time += ddt;
		  iteration = 0;
		}
		/*
		if (gvctlen < tolerance || iteration > maxiteration-1)
		{
		  nlap++;
		  time += ddt;
		  iteration = 0;
		}
		*/
		iteration++;
		setincrement((wmenu.childs+2)->hwnd,laps,nlap,maxiteration,iteration);


		if(iteration==1)
		{
			clearwindow(*(wdraw.childs+1));
			drawarclmframe((wdraw.childs+1)->hdcC,(wdraw.childs+1)->vparam,*af,0,ONSCREEN);
			overlayhdc(*(wdraw.childs + 1), SRCPAINT);                  /*UPDATE DISPLAY.*/

			af->nlaps = nlap;


			/*
			if(nlap>3)
			{
				ddt = timestepcontrol(ddt, lapddisp_m, udd_m, lastudd_m, lastlastudd_m, beta, msize);
				sprintf(string, "ddt: %e time: %e\n",ddt, time);
				errormessage(string);
			}
			*/

			/*LOCAL INCREMENTAL IN THIS LAP.*/
			for(i = 0; i < nshell; i++)
			{
			  outputmemoryshell(shells,mshell,i);
			}

			lastloadfactor = loadfactor;

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
#if 0/*FOR HINGE*/
			if(time<=1.0)
			{
			  loadfactor = 0.012 * time;
			}
			else
			{
			  loadfactor = 0.012;
			}
#endif
#if 1/*FOR HINGE*/
			if(time<=10.0)
			{
			  loadfactor = 0.0012* time;
			}
			else
			{
			  loadfactor = 0.012;
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

			for (i = 0; i < msize; i++)
			{
				*(lapddisp  + i) = 0.0;/*lapddisp : INCREMENTAL TRANSITION & ROTATION IN THIS LAP.*/
				*(lastddisp + i) = *(ddisp + i);

				*(lastlastudd_m  + i) = *(lastud_m  + i);

				*(lastud_m  + i) = *(ud_m  + i);
				*(lastudd_m + i) = *(udd_m + i);

				*(udinit_m + i)  =- (gamma / beta - 1.0) **(lastud_m + i)
								  - (gamma / (2.0 * beta) - 1.0) * ddt **(lastudd_m + i);
				*(uddinit_m + i) =- (1.0 / (beta * ddt)) **(lastud_m + i)
								  - (1.0 / (2.0 * beta) - 1.0) **(lastudd_m + i);
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
		  //assemelem_DYNA(elems, melem, nelem, constraintmain, NULL, gmtx, iform, ddisp);
		  assemshell_DYNA(shells, mshell, nshell, constraintmain,
						  gmtx,
						  iform, lastddisp, ddisp, lapddisp,
						  ud_m, udd_m,
						  alpham, alphaf, xi, beta, ddt);
		}


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

		/*CROUT LU DECOMPOSITION.*/
		nline = croutlu(gmtx, confs, msize, &determinant, &sign, gcomp1);
		/*DECOMPOSITION FAILED*/
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

			gfree(gmtx, nnode);  /*FREE GLOBAL MATRIX.*/
			free(gvct);
			free(funbalance);
			free(finternal);
			free(fexternal);
			return 1;
		}
		nline = forwardbackward(gmtx, gvct, confs, msize, gcomp1);

		sprintf(string, "LAP: %4d ITER: %2d {LOAD}= %5.8f {RESD}= %1.5e {DET}= %5.5f {SIGN}= %2.0f {BCL}= %1d {EPS}= %1.5e {V}= %5.5f {TIME}= %5.5f\n",
			nlap, iteration, loadfactor,/*residual*/gvctlen, determinant, sign, 0, 0.0, 0.0, time);
		fprintf(ffig, "%s", string);
		errormessage(string);

		for (ii = 0; ii < msize; ii++)
		{
			*(gvct + ii) = *(gvct + *(constraintmain + ii));
		}

		updaterotation(lapddisp, gvct, nnode);
		updaterotation(ddisp, gvct, nnode);

		lapddisp_m = pullback(lastddisp,lapddisp,nnode);
		for(ii = 0; ii < msize; ii++)
		{
			*(ud_m + ii)  = (gamma / (beta * ddt)) **(lapddisp_m + ii) + *(udinit_m + ii);
			*(udd_m + ii) = (1.0 / (beta * ddt * ddt)) **(lapddisp_m + ii) + *(uddinit_m + ii);
		}
		free(lapddisp_m);

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

			laptime("\0",t0);
			return 1;
		  }

		  //t2=clock();
		  //if((t2-t1)/CLK_TCK>=WAIT) break;               /*CONTINUE AFTER WAITING.*/
		}
#endif
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

	gfree(gmtx, nnode);  	/*FREE GLOBAL MATRIX.*/
	free(fexternal);		/*FREE VECTOR*/
	free(finternal);

	errormessage(" ");
	errormessage("COMPLETED.");

	memory1=availablephysicalmemoryEx("REMAIN:");
	sprintf(string,"CONSUMPTION:%ld[BYTES]",(memory0-memory1));
	errormessage(string);

	return 0;
}













#if 0
		if(solver==1)
		{

			for (i = 0; i < msize; i++)
			{
				*(mvct + i) = 0.0;
				*(cvct + i) = 0.0;
			}
			/*ELEMENTS ASSEMBLAGE.*/
			for (i = 1; i <= nshell; i++)
			{
				inputshell(shells, mshell, i - 1, &shell);
				nnod = shell.nnod;
				shell.sect = (shells + i - 1)->sect;
				loffset = (int*)malloc(6 * nnod * sizeof(int));
				for (ii = 0; ii < nnod; ii++)
				{
					for (jj = 0; jj < 6; jj++)
					{
						*(loffset + (6 * ii + jj)) = 6 * (shell.node[ii]->loff) + jj;
					}
				}

				/*INITIAL CONFIGURATION.*/
				for (ii = 0; ii < nnod; ii++)
				{
					inputnode(iform, shell.node[ii]);
				}
				drccosinit = shelldrccos(shell);
				gforminit = extractshelldisplacement(shell, iform);                 /*{Xg}*/
				eforminit = extractlocalcoord(gforminit,drccosinit,nnod);     		/*{Xe}*/


				//Kp = assemshellemtx(shell, C, B);
				Ke = assemshellemtx(shell);   						/*[Ke]*/
				Me = assemshellmmtx(shell, drccosinit);          					/*[Me]*/

				/*DEFORMED CONFIGURATION OF LAST ITERATION.*/
				for (ii = 0; ii < nnod; ii++)
				{
					inputnode(ddisp, shell.node[ii]);
				}
				drccos = shelldrccos(shell);
				gform = extractshelldisplacement(shell, ddisp);                     /*{Xg+Ug}*/
				eform = extractlocalcoord(gform,drccos,nnod);                 		/*{Xe+Ue}*/

				T = transmatrixIII(drccos, nnod);         					        /*[T]*/
				Tt = matrixtranspose(T, 6 * nnod);                  				/*[Tt]*/

				Me = transformationIII(Me, T, 6*nnod);/*[Tt][M][T]*/

				masstotal = 0.0;
				massdiag = 0.0;
				for (ii = 0; ii < 18; ii++)
				{
					for (jj = 0; jj < 18; jj++)
					{
						masstotal += Me[ii][jj];
						if (ii == jj && ii % 6 < 3)massdiag += Me[ii][jj];
					}
				}
				for (ii = 0; ii < nnod; ii++)
				{
					for (jj = 0; jj < 6; jj++)
					{
						*(mvct + *(constraintmain + *(loffset + (6 * ii + jj)))) += *(*(Me + 6 * ii + jj) + 6 * ii + jj) * masstotal / massdiag;
					}
				}



				edisp = extractdeformation(eforminit, eform, nnod);           		/*{Ue}*/
				einternal = matrixvector(Ke, edisp, 6 * nnod);      				/*{Fe}=[Ke]{Ue}*/


				midHPT = transmatrixHPT(eform, edisp, T, nnod);
				midTtPtHt = matrixtranspose(midHPT, 6 * nnod);

				midginternal = matrixvector(midTtPtHt, einternal, 6 * nnod);

				epressure = assemshellpvct(shell, drccos);                		    /*{Pe}*/
				midgpressure = matrixvector(Tt, epressure, 6 * nnod);



				shellstress = matrixvector(DBe, edisp, 6 * nnod);
				gvel_m = extractshelldisplacement(shell, ud_m);
				gmomentum_m = matrixvector(Me, gvel_m, 6 * nnod);

				/*GLOBAL VECTOR & MATRIX ASSEMBLY.*/
				for (ii = 0; ii < nnod; ii++)
				{
					for (jj = 0; jj < 6; jj++)
					{
						*(finternal + *(constraintmain + *(loffset + (6 * ii + jj)))) += *(midginternal + 6 * ii + jj);
						*(fpressure + *(constraintmain + *(loffset + (6 * ii + jj)))) += *(midgpressure + 6 * ii + jj);
						(shells+i-1)->stress[ii][jj]=*(einternal+6*ii+jj);
						(mshell+i-1)->stress[ii][jj]=*(einternal+6*ii+jj);

						(mshell+i-1)->stress[ii][jj]=*(shellstress+6*ii+jj);
					}
				}





				(mshell+i-1)->SE = 0.0;
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

				(mshell+i-1)->KE = 0.0;
				(mshell+i-1)->KEt = 0.0;
				(mshell+i-1)->KEr = 0.0;
				for (ii = 0; ii < nnod; ii++)
				{
					for (jj = 0; jj < 3; jj++)
					{
						(mshell+i-1)->KEt += 0.5 * *(gvel_m + 6 * ii + jj) * *(gmomentum_m + 6 * ii + jj);
					}
					for (jj = 3; jj < 6; jj++)
					{
						(mshell+i-1)->KEr += 0.5 * *(gvel_m + 6 * ii + jj) * *(gmomentum_m + 6 * ii + jj);
					}
				}
				(mshell+i-1)->KE = (mshell+i-1)->KEt + (mshell+i-1)->KEr;

				free(gvel_m);
				free(gmomentum_m);
				free(shellstress);




				/*MEMORY FREE : INITIAL CONFIGURATION.*/
				free(loffset);
				freematrix(drccosinit, 3);
				free(gforminit);
				free(eforminit);
				freematrix(Ke, 6 * nnod);
				freematrix(Me, 6 * nnod);


				/*MEMORY FREE : DEFORMED CONFIGURATION OF LAST ITERATION.*/
				freematrix(drccos, 3);
				freematrix(T, 6 * nnod);
				freematrix(Tt, 6 * nnod);
				freematrix(midHPT, 6 * nnod);
				freematrix(midTtPtHt, 6 * nnod);
				free(gform);
				free(eform);
				free(edisp);
				free(einternal);
				free(midginternal);


				free(epressure);
				free(midgpressure);
			}


		}


		if(solver==1)
		{



		rightud_m = (double*)malloc(msize * sizeof(double));
		leftud_m = (double*)malloc(msize * sizeof(double));
		mvct = (double *)malloc(msize * sizeof(double));
		cvct = (double *)malloc(msize * sizeof(double));

		for (i = 0; i < msize; i++)
		{
			*(rightud_m + i) = 0.0;
			*(leftud_m + i) = 0.0;
		}




			for (i = 0; i < msize; i++)
			{
				*(fbaseload + i) = *(fdeadload + i)+*(fpressure + i);
				*(fexternal + i) = loadfactor * *(fbaseload + i);
				*(funbalance + i) = *(fexternal + i) - *(finternal + i);
				if ((confs + i)->iconf == 1)
				{
					*(freaction + i) = *(funbalance + i);
					*(funbalance + i) = 0.0;
					*(fbaseload + i) = 0.0;
				}
			}



			leftleftKE = leftKE;
			leftKE = rightKE;
			KE = 0.0;
			rightKE = 0.0;

			//funbalance_m = pullback(ddisp,funbalance,nnode);
			for (i = 0; i < msize; i++)
			{

				if (time == 0.0 || PEAKFLAG)
				{
					*(rightud_m + i) = *(funbalance + i) * ddt / (2.0 * *(mvct + i));
					*(leftud_m  + i) = -*(rightud_m  + i);

					*(ud_m  + i) = 0.0;
					*(udd_m  + i) = 2.0 * *(rightud_m  + i) / ddt;
				}
				else
				{
					*(leftud_m  + i) = *(rightud_m  + i);
					*(rightud_m  + i) = (*(funbalance + i) * ddt / *(mvct + i)) + *(leftud_m  + i);
					//*(rightud_m  + i) = ( *(funbalance_m + i) + (*(mvct + i) / ddt - *(cvct + i) / 2.0) * *(leftud_m  + i) ) / ( *(mvct + i) / ddt + *(cvct + i) / 2.0 );

					*(ud_m  + i) = (*(rightud_m  + i) + *(leftud_m  + i)) / 2.0;
					*(udd_m  + i) = (*(rightud_m  + i) - *(leftud_m  + i)) / ddt;
				}
				KE += 0.5 * *(mvct + i) * *(ud_m  + i) * *(ud_m  + i);
				rightKE += 0.5 * *(mvct + i) * *(rightud_m  + i) * *(rightud_m  + i);

			}
			sprintf(string, "LAP: %4d ITER: %2d {LOAD}= %5.8f {RESD}= %1.5e {DET}= %5.5f {SIGN}= %2.0f {BCL}= %1d {EPS}= %1.5e {V}= %5.5f {TIME}= %5.5f\n",
				nlap, 1, loadfactor,/*residual*/0, 0, 0, 0, 0.0, 0.0, time);
			fprintf(ffig, "%s", string);
			errormessage(string);


			/*FOR DINAMIC RELAXATION*/
			PEAKFLAG = 0;
			if (method==2 && rightKE < leftKE)
			{
				q = (leftKE - rightKE) / (2 * leftKE - rightKE - leftleftKE);
				for (i = 0; i < msize; i++)
				{
					*(gvct + i) = -q * *(leftud_m  + i) * ddt;
				}
				leftKE = 0.0;
				rightKE = 0.0;
				PEAKFLAG = 1;
			}
			else
			{
				for (i = 0; i < msize; i++)
				{
					*(gvct + i) = *(rightud_m  + i) * ddt;
				}
			}
			//gvct_s = pushforward(ddisp,gvct,nnode);

			//gvctlen = vectorlength(gvct,msize);

			for (ii = 0; ii < msize; ii++)
			{
				if (*(constraintmain + ii) != ii)
				{
					*(gvct + ii) = *(gvct + *(constraintmain + ii));
				}
			}
			updaterotation(ddisp, gvct, nnode);


			//free(funbalance_m);
			//free(gvct_s);


			nlap++;
		}




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
				//Me = assemmmtx(elem, drccosinit);          					/*[Me]*/

				/*DEFORMED CONFIGURATION OF LAST LAP.*/
				for (ii = 0; ii < nnod; ii++)
				{
					inputnode(lastddisp, elem.node[ii]);
				}
				lastgform = extractdisplacement(elem, lastddisp);
				lastdrccos = updatedrccos(drccosinit, gforminit, lastgform);
				lasteform = extractlocalcoord(lastgform,lastdrccos,nnod);

				lastT = transmatrixIII(lastdrccos, nnod);         					/*[T]*/
				lastTt = matrixtranspose(lastT, 6 * nnod);                  		/*[Tt]*/

				lastedisp = extractdeformation(eforminit, lasteform, nnod);         /*{Ue}*/
				lasteinternal = matrixvector(Ke, lastedisp, 6 * nnod);      		/*{Fe}=[Ke]{Ue}*/
				lastHPT = transmatrixHPT(lasteform, lastedisp, lastT, nnod);

				lastgacc_m = extractdisplacement(elem, lastudd_m);
				lastR = pushforwardmtx(lastgform, nnod);
				lastRt = matrixtranspose(lastR, 6 * nnod);
				lastginertial_m = matrixvector(Me, lastgacc_m, 6 * nnod);
				lastginertial = pushforward(lastgform, lastginertial_m, nnod);

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
				HPT = transmatrixHPT(eform, edisp, T, nnod);

				gacc_m = extractdisplacement(elem, udd_m);
				R=pushforwardmtx(gform,nnod);
				ginertial_m = matrixvector(Me, gacc_m, 6 * nnod);
				ginertial = pushforward(gform, ginertial_m, nnod);

				lapgform = extractdisplacement(elem, lapddisp);
				lapH = blockjacobimtx(lapgform, NULL, NULL, nnod);


				/*(24):MID-POINT TRANSFORMATION MATRIX.*/
				midHPT = midpointmtx(HPT, lastHPT, alphaf, 6 * nnod);
				midTtPtHt = matrixtranspose(midHPT, 6 * nnod);

				midT = midpointmtx(T, lastT, alphaf, 6 * nnod);
				midTt = matrixtranspose(midT, 6 * nnod);

				/*(21)&(22):MID-POINT FORCE VECTOR.*/
				mideinternal = midpointvct(einternal, lasteinternal, alphaf-xi, 6*nnod);
				midginternal = matrixvector(midTtPtHt, mideinternal, 6 * nnod);

				midginertial = midpointvct(ginertial, lastginertial, alpham, 6*nnod);

				/*MID-POINT MASS & STIFFNESS MATRIX.*/
				Keff=assemtmtxCR_DYNA(eform, edisp, mideinternal, T, Ke, midTtPtHt, HPT,
									  ginertial, Me, R, lastRt, lapH,
									  alphaf, alpham, xi, beta, ddt, nnod);
				symmetricmtx(Keff, 6*nnod);
				assemgstiffnesswithDOFelimination(gmtx, Keff, &elem, constraintmain);

				/*GLOBAL VECTOR & MATRIX ASSEMBLY.*/
				for (ii = 0; ii < nnod; ii++)
				{
					for (jj = 0; jj < 6; jj++)
					{
						*(finertial + *(loffset + (6 * ii + jj))) += *(midginertial + 6 * ii + jj);
						*(finternal + *(loffset + (6 * ii + jj))) += *(midginternal + 6 * ii + jj);
					}
				}
			}


#endif



double timestepcontrol(double ddt, double* lapddisp_m, double* udd_m, double* lastudd_m, double* lastlastudd_m, double beta, int msize)
{
	int i;
	char string[100];
	double eta_upper  = 0.06;
	double eta_target = 0.05;
	double eta_lower  = 0.04;
	double eta, errornorm, dispnorm;
	double* error;

	error = (double *)malloc(msize * sizeof(double));
	for (i = 0; i < msize; i++)
	{
#if 0
	  /*Time Step Control Algorithm by Schweizerhof & Riccius*/
	  *(error + i) = ddt * ddt * ((beta-1.0/8.0)**(udd_m + i) + (1.0/12.0 - beta)**(lastudd_m + i) + 1.0/24.0**(lastlastudd_m + i));
#endif
#if 1
	  /*Time Step Control Algorithm by Zienkiewics & Xie*/
	  *(error + i) = ddt * ddt * (beta-1.0/6.0) * (*(udd_m + i) - *(lastudd_m + i));
#endif
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
	return ddt;
}


