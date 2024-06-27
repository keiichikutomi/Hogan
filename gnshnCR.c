/* ========================================================= */
/* PROGRAM GNSHN CR FOR OUTER SURFACE OF LUNAR MARS PROJECT  */
/* DYNAMIC ANALYSIS FOR LARGE DEFORMATION                    */
/* USING CR FORMULATION & EMM ALGORITHM                      */
/* CODED BY KEIICHI KUTOMI SINSE 2024.05.26                  */
/* ========================================================= */

/*MATERIAL & SPATIAL FORM VARIABLES.*/
double* pullback(double* ddisp, double* gvct_s, int nnode)
{
	int i,n;
	double* rvct, * vct_s, * vct_m, * gvct_m;
	double** rmtx, ** trmtx;

	gvct_m = (double*)malloc(6*nnode * sizeof(double));
	vct_s = (double*)malloc(3 * sizeof(double));
	rvct = (double*)malloc(3 * sizeof(double));
	for (n = 0; n < nnode; n++)
	{
		for (i = 0; i < 3; i++)
		{
			*(rvct + i) = *(ddisp + 6*n+3+i);
			*(vct_s + i) = *(gvct_s + 6*n+3+i);
		}
		rmtx = rotationmtx(rvct);
		trmtx = matrixtranspose(rmtx, 3);
		vct_m = matrixvector(trmtx, vct_s, 3);
		for (i = 0; i < 3; i++)
		{
			*(gvct_m + 6*n+i) = *(gvct_s + 6*n+i);
			*(gvct_m + 6*n+3+i) = *(vct_m + i);
		}
	}
	free(rvct);
	free(vct_s);
	free(vct_m);
	freematrix(rmtx,3);
	freematrix(trmtx,3);
	return gvct_m;
}

double** pullbackmtx(double* gform, int nnod)
{
	int i,j,n;
	double* rvct;
	double** rmtx,** trmtx, **Rt;

	Rt = (double**)malloc(6*nnod * sizeof(double*));
	for (i = 0; i < 6*nnod; i++)
	{
		*(Rt + i) = (double*)malloc(6*nnod * sizeof(double));
		for (j = 0; j < 6*nnod; j++)
		{
			*(*(Rt + i) + j) = 0.0;
		}
	}
	rvct = (double*)malloc(3 * sizeof(double));
	for (n = 0; n < nnod; n++)
	{
		for (i = 0; i < 3; i++)
		{
			*(rvct  + i) = *(gform + 6*n+3+i);
		}
		rmtx = rotationmtx(rvct);
		trmtx = matrixtranspose(rmtx, 3);
		for (i = 0; i < 3; i++)
		{
			for (j = 0; j < 3; j++)
			{
				*(*(Rt + 6 * n + 3 + i) + 6 * n + 3 + j) = *(*(trmtx + i) + j);
				if (i == j)*(*(Rt + 6 * n + i) + 6 * n + j) = 1.0;
			}
		}
	}

	free(rvct);
	freematrix(rmtx,3);
	freematrix(trmtx,3);
	return Rt;
}

double* pushforward(double* ddisp, double* gvct_m, int nnode)
{
	int i,n;
	double* rvct, * vct_s, * vct_m, * gvct_s;
	double** rmtx;

	gvct_s = (double*)malloc(6*nnode * sizeof(double));
	vct_m = (double*)malloc(3 * sizeof(double));
	rvct = (double*)malloc(3 * sizeof(double));
	for (n = 0; n < nnode; n++)
	{
		for (i = 0; i < 3; i++)
		{
			*(rvct  + i) = *(ddisp + 6*n+3+i);
			*(vct_m + i) = *(gvct_m + 6*n+3+i);
		}
		rmtx = rotationmtx(rvct);
		vct_s = matrixvector(rmtx, vct_m, 3);
		for (i = 0; i < 3; i++)
		{
			*(gvct_s + 6*n+i) = *(gvct_m + 6*n+i);
			*(gvct_s + 6*n+i+3) = *(vct_s + i);
		}
	}
	free(rvct);
	free(vct_s);
	free(vct_m);
	freematrix(rmtx,3);
	return gvct_s;
}


double** pushforwardmtx(double* gform, int nnod)
{
	int i,j,n;
	double* rvct;
	double** rmtx, **R;

	R = (double**)malloc(6*nnod * sizeof(double*));
	for (i = 0; i < 6*nnod; i++)
	{
		*(R + i) = (double*)malloc(6*nnod * sizeof(double));
		for (j = 0; j < 6*nnod; j++)
		{
			*(*(R + i) + j) = 0.0;
		}
	}
	rvct = (double*)malloc(3 * sizeof(double));
	for (n = 0; n < nnod; n++)
	{
		for (i = 0; i < 3; i++)
		{
			*(rvct  + i) = *(gform + 6*n+3+i);
		}
		rmtx = rotationmtx(rvct);
		for (i = 0; i < 3; i++)
		{
			for (j = 0; j < 3; j++)
			{
				*(*(R + 6 * n + 3 + i) + 6 * n + 3 + j) = *(*(rmtx + i) + j);
				if (i == j)*(*(R + 6 * n + i) + 6 * n + j) = 1.0;
			}
		}
	}

	free(rvct);
	freematrix(rmtx,3);
	return R;
}

/*MID-POINT VARIABLES.*/
double* midpointvct(double* vct,double* lastvct,double alpha,int size)
{
	int i;
	double* midvct;

	midvct = (double*)malloc(size * sizeof(double));
	for (i=0;i<size;i++)
	{
		*(midvct+i)=(1.0-alpha)**(vct+i)+alpha**(lastvct+i);
	}
	return midvct;
}

double** midpointmtx(double** mtx,double** lastmtx,double alpha,int size)
{
	int i,j;
	double** midmtx;

	midmtx= (double**)malloc(size * sizeof(double*));
	for (i=0;i<size;i++)
	{
		*(midmtx+i) = (double*)malloc(size * sizeof(double));
		for(j=0;j<size;j++)
		{
			*(*(midmtx+i)+j)=(1.0-alpha)**(*(mtx+i)+j)+alpha**(*(lastmtx+i)+j);
		}
	}
	return midmtx;
}


double **assemtmtxCR_DYNA(double* eform, double* edisp, double* mideinternal, double** T, double** Ke, double** midTtPtHt, double** HPT,
						  double* ginertial, double** Me, double** R, double** lastRt, double** lapH,
						  double alphaf, double alpham, double xi, double beta, double ddt, int nnod)
{
	int i,j,n;
	double** Kint1, ** Kint2;
	double** Kmas1, ** Kmas2;
	double* mvct;
	double** mspin;
	double** TPHK;
	double** RM, ** RMR;
	double** Keff;

	Kint1 = assemgmtxCR(eform, edisp, mideinternal, NULL, T, NULL, nnod);

	TPHK = matrixmatrix(midTtPtHt, Ke, 6*nnod);
	Kint2 = matrixmatrix(TPHK, HPT, 6*nnod);

	mvct = (double*)malloc(3 * sizeof(double));
	Kmas1 = (double**)malloc(6*nnod * sizeof(double*));
	for (i = 0; i < 6*nnod; i++)
	{
		*(Kmas1 + i) = (double*)malloc(6*nnod * sizeof(double));
		for (j = 0; j < 6*nnod; j++)
		{
			*(*(Kmas1 + i) + j)=0.0;
		}
	}
	for (n = 0; n < nnod; n++)
	{
		for (i = 0; i < 3; i++)
		{
			*(mvct + i) = *(ginertial + 6 * n + 3 + i);
		}
		mspin = spinmtx(mvct);
		for (i = 0; i < 3; i++)
		{
			for (j = 0; j < 3; j++)
			{
				*(*(Kmas1 + 6 * n + 3 + i) + 6 * n + 3 + j) = - *(*(mspin + i) + j);
			}
		}
	}

	RM = matrixmatrix(R, Me, 6*nnod);
	RMR = matrixmatrix(RM, lastRt, 6*nnod);
	Kmas2 = matrixmatrix(RMR, lapH, 6*nnod);

	Keff = (double**)malloc(6*nnod * sizeof(double*));
	for (i = 0; i < 6*nnod; i++)
	{
		*(Keff + i) = (double*)malloc(6*nnod * sizeof(double));
		for (j = 0; j < 6*nnod; j++)
		{
			*(*(Keff + i) + j) = (1-alphaf)     **(*(Kint1 + i) + j)
							   + (1-alphaf + xi)**(*(Kint2 + i) + j)
							   + (1-alpham)     **(*(Kmas1 + i) + j)
							   + (1-alpham)     **(*(Kmas2 + i) + j) / (beta * ddt * ddt)
							   ;
		}
	}
	freematrix(Kint1, 6 * nnod);
	freematrix(TPHK,  6 * nnod);
	freematrix(Kint2, 6 * nnod);
	free(mvct);
	freematrix(mspin, 3);
	freematrix(Kmas1, 6 * nnod);
	freematrix(RM,    6 * nnod);
	freematrix(RMR,   6 * nnod);
	freematrix(Kmas2, 6 * nnod);
	return Keff;
}

int gnshnCR(struct arclmframe* af)
{
	DWORDLONG memory;
	FILE *fdata, *fin, *fonl, * fdsp, * fexf, * finf, * fubf, * frct, * finr, * fstr, * fene, * ffig, * fbcl, * feig, * flog;         /*FILE 8 BYTES*/
	char dir[] = DIRECTORY;
	char s[80], string[400], inpname[50], fname[50];

	int i, ii, jj;
	int nnode, nelem, nshell, nsect, nreact, nconstraint;
	long int msize;

	int nlap, laps;  /*LAP COUNT*/
	double ddt = 0.001;/*TIME INCREMENT[sec]*/
	double time = 0.0;/*TOTAL TIME[sec]*/

	/***FOR ELEMENT***/
	double volume;
	int nnod;
	int* loffset;
	double** Me, ** Ke, ** DBe;
	double* gforminit, * lastgform, * gform;                      	/*GLOBAL COORDINATION OF ELEMENT*/
	double* eforminit, * lasteform, * eform;                      	/*LOCAL COORDINATION OF ELEMENT*/
	double* lastedisp, * edisp;                     			 	/*LOCAL DEFORMATION OF ELEMENT*/
	double* lastgacc, * gacc;
	double* lastgacc_m, * gacc_m;
	double** drccosinit, ** lastdrccos, ** drccos;  				/*ELEMENT COORDINATION MATRIX*/
	double** lastT, ** T;
	double** lastTt, ** Tt;
	double** lastHPT, ** HPT;
	double** midHPT, **midTtPtHt;
	double** lastR, ** R;
	double** lastRt, ** Rt;
	double** lapH;
	double* lapgform;

	double* lasteinternal, * einternal;
	double* mideinternal;
	double* midginternal;

	double* lastginertial_m, * ginertial_m;
	double* lastginertial, * ginertial;
	double* midginertial;

	double* lastgpressure, * gpressure;
	double* midgpressure;
	double** Keff;

	double* shellstress;                           /*σx,σy,τxy,Mx,My,Mxy OF ELEMENT*/
	double Ep, Eb, Ee;                             /*STRAIN ENERGY OF ELEMENT*/

	double area, volumetotal;
	double masstotal,massdiag;

	/***FOR INCREMENTAL***/
	int iteration;
	int maxiteration = 20;
	double residual;
	double tolerance = 1.0E-8;

	double momentumlinear,momentumangular;

	long int nline;
	double det, sign;
	double loadfactor = 0.0;
	double lastloadfactor = 0.0;
	//double lambda = 0.0;

	/*ENERGY MOMENTUM METHOD'S PARAMETER.*/
	double rho = 1.0;
	double alpham = (2.0*rho-1)/(rho+1);  /*MID-POINT USED TO EVALUATE INERTIAL FORCE*/
	double alphaf = rho/(rho+1);        /*MID-POINT USED TO EVALUATE INTERNAL FORCE*/
	double xi = 0.0;   	 				/*NUMERICAL DISSIPATION COEFFICIENT*/

	/*NEWMARK-BETA'S PARAMETER.*/
	double beta = pow(1-alpham+alphaf,2)/4.0;
	double gamma = (0.5-alpham+alphaf)/4.0;

	/*FOR READING ANALYSISDATA*/
	/*ANALYSIS SETTING*/
	int nstr, pstr, readflag;
	char **data, *filename;

	/*NODE DEFORMATION AT ANALYSIS ENDING*/
	int fnode=NULL,fnodedrc=NULL;
	double fnodemin, fnodemax;

	int outputmode   = 0;/*0:ConvergedLaps.1:AllLaps.*/
	int pinpointmode = 0;/*0:NotPinpointing.1:BisecPinpointing.2:ExtendedSystemPinpointing.*/


	fdata = fgetstofopenII(dir, "r", "analysisdata.txt");
	if (fdata == NULL)
	{
		errormessage("couldn't open analysisdata.txt\n");
		getchar();
		exit(EXIT_FAILURE);
	}
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
						filename = *(data + pstr);
					}
					if (!strcmp(*(data + pstr), "LAPS"))
					{
						pstr++;
						laps = (int)strtol(*(data + pstr), NULL, 10);
					}
					if (!strcmp(*(data + pstr), "TIMEINCREMENT"))
					{
						pstr++;
						ddt = (int)strtol(*(data + pstr), NULL, 10);
					}
					if (!strcmp(*(data + pstr), "NNODE"))
					{
						pstr++;
						nnode = (int)strtol(*(data + pstr), NULL, 10);
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
					else
					{
						pstr++;
					}
				}
			}
		}
	}
	sprintf(string,"FILENAME : %s\n", filename);
	errormessage(string);
	if (outputmode == 0)errormessage("OUTPUT CONVERGED RESULT\n");
	if (outputmode == 1)errormessage("OUTPUT ALL RESULT\n");

	sprintf(string, "INITIAL:");
	memory = availablephysicalmemoryEx(string);   /*MEMORY AVAILABLE*/


	///INPUT FILE OPEN & OUTPUT FILE SETTING///

	strcat(filename,".inl");
	fin = fgetstofopenII(dir, "r", filename);              /*OPEN INPUT FILE*/

	strcpy(inpname, filename);
	char* dot = strrchr(inpname, '.');
	*dot = '\0';
	snprintf(fname, sizeof(fname), "%s.%s", inpname, "dsp");
	fdsp = fopen(fname, "w");
	snprintf(fname, sizeof(fname), "%s.%s", inpname, "inf");
	finf = fopen(fname, "w");
	snprintf(fname, sizeof(fname), "%s.%s", inpname, "exf");
	fexf = fopen(fname, "w");
	snprintf(fname, sizeof(fname), "%s.%s", inpname, "ubf");
	fubf = fopen(fname, "w");
	snprintf(fname, sizeof(fname), "%s.%s", inpname, "rct");
	frct = fopen(fname, "w");
	snprintf(fname, sizeof(fname), "%s.%s", inpname, "inr");
	finr = fopen(fname, "w");

	snprintf(fname, sizeof(fname), "%s.%s", inpname, "str");
	fstr = fopen(fname, "w");
	snprintf(fname, sizeof(fname), "%s.%s", inpname, "ene");
	fene = fopen(fname, "w");

	snprintf(fname, sizeof(fname), "%s.%s", inpname, "onl");
	fonl = fopen(fname, "w");

	snprintf(fname, sizeof(fname), "%s.%s", inpname, "fig");
	ffig = fopen(fname, "w");
	snprintf(fname, sizeof(fname), "%s.%s", inpname, "bcl");
	fbcl = fopen(fname, "w");
	snprintf(fname, sizeof(fname), "%s.%s", inpname, "eig");
	feig = fopen(fname, "w");
	snprintf(fname, sizeof(fname), "%s.%s", inpname, "log");
	flog = fopen(fname, "w");


	inputinitII(fin, &nnode, &nelem, &nshell, &nsect, &nconstraint);	/*INPUT INITIAL.*/
	msize = 6*nnode;                           							/*SIZE OF GLOBAL MATRIX.*/

	///MEMORY///
	struct osect* sects;
	struct onode* nodes;
	struct onode* ninit;
	struct owire elem;
	struct owire* elems;
	struct oshell shell;
	struct oshell* shells;
	struct oconf* confs;
	struct memoryelem* melem;
	struct memoryshell* mshell;
	long int mainoff;
	long int* constraintmain;

	sects = (struct osect*)malloc(nsect * sizeof(struct osect));
	nodes = (struct onode*)malloc(nnode * sizeof(struct onode));
	ninit = (struct onode*)malloc(nnode * sizeof(struct onode));
	elems = (struct owire*)malloc(nelem * sizeof(struct owire));
	shells = (struct oshell*)malloc(nshell * sizeof(struct oshell));
	confs = (struct oconf*)malloc(msize * sizeof(struct oconf));
	melem = (struct memoryelem*)malloc(nelem * sizeof(struct memoryelem));
	mshell = (struct memoryshell*)malloc(nshell * sizeof(struct memoryshell));
	constraintmain = (long int*)malloc(msize * sizeof(long int));
    for (i = 0; i < msize; i++)
	{
		*(constraintmain + i) = i;
	}

	///POSITION VECTOR INITIALIZATION///
	/*IN SPATIAL FORM*/
	double* iform,* ddisp,* lastddisp,* lapddisp,* lapddisp_m;

	iform = (double*)malloc(msize * sizeof(double));		/*INITIAL*/
	ddisp = (double*)malloc(msize * sizeof(double));		/*LATEST ITERATION*/
	lastddisp = (double*)malloc(msize * sizeof(double));	/*LATEST LAP*/
	lapddisp  = (double*)malloc(msize * sizeof(double));		/*INCREMENT IN THE LAP*/
	lapddisp_m = (double*)malloc(msize * sizeof(double));		/*INCREMENT IN THE LAP*/


	///VELOSITY & ACCELERATION VECTOR INITIALIZATION///
	/*ud : VELOSITY, udd : ACCELERATION*/
	/*_m : ROTATIONAL DOFs ARE REPRESENTED IN SPATIAL & MATERIAL FORM*/
	double* udinit_m,* uddinit_m;
	//double* lastud,* lastudd;
	double* lastud_m,* lastudd_m;
	//double* ud,* udd;
	double* ud_m,* udd_m;

	udinit_m = (double*)malloc(msize * sizeof(double));  	   	/*NEWMARK INITIAL IN LAP*/
	uddinit_m = (double*)malloc(msize * sizeof(double));

	//lastud = (double*)malloc(msize * sizeof(double));  		/*LATEST LAP IN SPATIAL*/
	//lastudd = (double*)malloc(msize * sizeof(double));

	lastud_m = (double*)malloc(msize * sizeof(double));     /*LATEST LAP IN MATERIAL*/
	lastudd_m = (double*)malloc(msize * sizeof(double));

	//ud = (double*)malloc(msize * sizeof(double)); 			/*LATEST ITERATION IN SPATIAL*/
	//udd = (double*)malloc(msize * sizeof(double));

	ud_m = (double*)malloc(msize * sizeof(double));  		/*LATEST ITERATION IN MATERIAL*/
	udd_m = (double*)malloc(msize * sizeof(double));

	for (i = 0; i < msize; i++)
	{
		*(udinit_m + i) = 0.0;
		*(uddinit_m + i) = 0.0;

		//*(lastud + i) = 0.0;
		//*(lastudd + i) = 0.0;

		*(lastud_m + i) = 0.0;
		*(lastudd_m + i) = 0.0;

		//*(ud + i) = 0.0;
		//*(udd + i) = 0.0;

		*(ud_m + i) = 0.0;
		*(udd_m + i) = 0.0;
	}

	///FORCE VECTOR INITIALIZATION///
	double* fexternal, * fpressure, * finternal,* finertial,* fdamping,* funbalance, * freaction;

	fexternal = (double*)malloc(msize * sizeof(double));         /*BASE EXTERNAL FORCE VECTOR.*/
	fpressure = (double*)malloc(msize * sizeof(double));         /*BASE PRESSURE FORCE VECTOR.*/
	finternal = (double*)malloc(msize * sizeof(double));         /*INTERNAL FORCE VECTOR.*/
	finertial = (double*)malloc(msize * sizeof(double));         /*INERTIAL FORCE VECTOR.*/
	//fdamping= (double*)malloc(msize * sizeof(double));         /*DAMPING FORCE VECTOR.*/
	funbalance = (double*)malloc(msize * sizeof(double));        /*UNBALANCED FORCE VECTOR.*/
	freaction = (double*)malloc(msize * sizeof(double));         /*INERTIAL FORCE VECTOR.*/

	///GLOBAL VECTOR & MATRIX INITIALIZATION FOR SOLVING EQUILIBRIUM EQUATION///
	double* gvct;
	gvct = (double*)malloc(msize * sizeof(double));

	struct gcomponent ginit = {0,0,0.0,NULL};
	struct gcomponent* gmtx, * g, * p, * gcomp1;/*GLOBAL MATRIX*/
	double gg;
	gmtx = (struct gcomponent*)malloc(msize * sizeof(struct gcomponent));  /* DIAGONALS OF MATRIX [K].*/

	for (i = 1; i <= msize; i++)
	{
		ginit.m = (unsigned short int)i;
		*(gmtx + (i - 1)) = ginit;
	}

	///ARCLMFRAME INITIALIZATION///
	free(af->sects);
	free(af->nodes);
	free(af->ninit);
	free(af->elems);
	free(af->shells);
	free(af->confs);
	free(af->ddisp);
	free(af->melem);
	free(af->mshell);
	free(af->constraintmain);

	af->sects = sects;
	af->nodes = nodes;
	af->ninit = ninit;
	af->elems = elems;
	af->shells = shells;
	af->confs = confs;
	af->ddisp = ddisp;
	af->melem = melem;
	af->mshell = mshell;
	af->constraintmain = constraintmain;

	///INITIALIZATION FROM INPUT FILE///
	inputtexttomemory(fin, af);
	fclose(fin);

	nnode = af->nnode;
	ninit = af->ninit;
	nelem = af->nelem;
	nshell = af->nshell;
	nsect = af->nsect;
	nreact = af->nreact;
	nconstraint = af->nconstraint;

	initialformCR(ninit, iform, nnode);           /*ASSEMBLAGE FORMATION.*/
	initialformCR(nodes, lastddisp, nnode);       /*ASSEMBLAGE FORMATION.*/
	initialformCR(nodes, ddisp, nnode);           /*ASSEMBLAGE FORMATION.*/

	initialelem(elems, melem, nelem);             /*ASSEMBLAGE ELEMENTS.*/
	initialshell(shells, mshell, nshell);         /*ASSEMBLAGE ELEMENTS.*/


	for (i = 0; i < msize; i++)
	{
		if (*(constraintmain + i) != i)
		{
			(confs + i)->iconf = (signed char)1;
		}
	}

	setviewpoint((wdraw.childs+0)->hwnd,arc,
						 &((wdraw.childs+1)->vparam));
	setviewparam((wmenu.childs+2)->hwnd,
						 (wdraw.childs+1)->vparam);
	clearwindow(*(wdraw.childs+1));
	drawarclmframe((wdraw.childs+1)->hdcC,
						   (wdraw.childs+1)->vparam,arc,0,ONSCREEN);

	//getincrement((wmenu.childs+2)->hwnd,&laps,&dsafety);

	///FOR DRAWING 1///
	GetAsyncKeyState(VK_LBUTTON);                   /*CLEAR KEY LEFT.*/
	GetAsyncKeyState(VK_RBUTTON);                  /*CLEAR KEY RIGHT.*/
	if(globaldrawflag==1)
	{
	  drawglobalaxis((wdraw.childs+1)->hdcC,(wdraw.childs+1)->vparam,0,0,255);  /*DRAW GLOBAL AXIS.*/
	}
	///FOR DRAWING 1///


	nlap = 1;
	iteration = 1;
	residual = 0.0;

	while (nlap <= laps)
	{
		af->nlaps = nlap;

		///FOR DRAWING 2///
		//setincrement((wmenu.childs+2)->hwnd,laps,nlap,maxiteration,iteration);
		if(iteration==1)clearwindow(*(wdraw.childs+1));
		///FOR DRAWING 2///


		if((outputmode == 0 && iteration == 1) || outputmode == 1)
		{
			sprintf(string, "LAP: %5ld / %5ld ITERATION: %5ld\n", nlap, laps, iteration);
			fprintf(fdsp, string);
			fprintf(finr, string);
			fprintf(finf, string);
			fprintf(fexf, string);
			fprintf(fubf, string);
			fprintf(frct, string);
			fprintf(fstr, string);
			fprintf(fene, string);
		}

		if (iteration == 1)
		{
			time+=ddt;
			lastloadfactor = loadfactor;
			loadfactor = 1.0e+5 * time;

			for (i = 0; i < msize; i++)
			{
				*(lapddisp  + i) = 0.0;
				*(lastddisp + i) = *(ddisp + i);
				*(lastud_m  + i) = *(ud_m  + i);
				*(lastudd_m + i) = *(udd_m + i);
			}
		}

		/*GLOBAL MATRIX INITIALIZATION.*/
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
		comps = msize;

		/*GLOBAL VECTOR INITIALIZATION.*/
		for (i = 0; i < msize; i++)
		{
			*(finertial  + i) = 0.0;
			*(finternal + i) = 0.0;
			*(fexternal + i) = 0.0;
			*(fpressure + i) = 0.0;
			*(freaction  + i) = 0.0;
			*(funbalance + i) = 0.0;
		}

		volume = 0.0;

		assemconf(confs, fexternal, 1.0, nnode);

		clearwindow(*(wdraw.childs + 1));

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
			drccosinit = shelldrccos(shell, &area);
			gforminit = extractshelldisplacement(shell, iform);                 /*{Xg}*/
			eforminit = extractlocalcoord(gforminit,drccosinit,nnod);     		/*{Xe}*/

			DBe = (double**)malloc(18 * sizeof(double*));
			for (ii = 0; ii < 18; ii++)
			{
				*(DBe + ii) = (double*)malloc(18 * sizeof(double));
				for (jj = 0; jj < 18; jj++)
				{
					*(*(DBe + ii) + jj) = 0.0;
				}
			}
			Ke = assemshellemtx(shell, drccosinit, DBe);   						/*[Ke]*/
			Me = assemshellmmtx(shell, drccosinit);          					/*[Me]*/

			/*DEFORMED CONFIGURATION OF LAST LAP.*/
			for (ii = 0; ii < nnod; ii++)
			{
				inputnode(lastddisp, shell.node[ii]);
			}
			lastdrccos = shelldrccos(shell, &area);
			lastT = transmatrixIII(lastdrccos, nnod);         					/*[T]*/
			lastTt = matrixtranspose(lastT, 6 * nnod);                  		/*[Tt]*/

			lastgform = extractshelldisplacement(shell, lastddisp);             /*{Xg+Ug}*/
			lasteform = extractlocalcoord(lastgform,lastdrccos,nnod);           /*{Xe+Ue}*/

			lastedisp = extractdeformation(eforminit, lasteform, nnod);         /*{Ue}*/
			lasteinternal = matrixvector(Ke, lastedisp, 6 * nnod);      		/*{Fe}=[Ke]{Ue}*/
			lastHPT = transmatrixHPT(lasteform, lastedisp, lastT, nnod);

			lastgacc_m = extractshelldisplacement(shell, lastudd_m);
			lastR = pushforwardmtx(lastgform, nnod);
			lastRt = matrixtranspose(lastR, 6 * nnod);
			lastginertial_m = matrixvector(Me, lastgacc_m, 6 * nnod);
			lastginertial = pushforward(lastgform, lastginertial_m, nnod);

			/*DEFORMED CONFIGURATION OF LAST ITERATION.*/
			for (ii = 0; ii < nnod; ii++)
			{
				inputnode(ddisp, shell.node[ii]);
			}
			drccos = shelldrccos(shell, &area);
			T = transmatrixIII(drccos, nnod);         					        /*[T]*/
			Tt = matrixtranspose(T, 6 * nnod);                  				/*[Tt]*/

			gform = extractshelldisplacement(shell, ddisp);                     /*{Xg+Ug}*/
			eform = extractlocalcoord(gform,drccos,nnod);                 		/*{Xe+Ue}*/

			edisp = extractdeformation(eforminit, eform, nnod);           		/*{Ue}*/
			einternal = matrixvector(Ke, edisp, 6 * nnod);      				/*{Fe}=[Ke]{Ue}*/
			HPT = transmatrixHPT(eform, edisp, T, nnod);

			gacc_m = extractshelldisplacement(shell, udd_m);
			R=pushforwardmtx(gform,nnod);
			ginertial_m = matrixvector(Me, gacc_m, 6 * nnod);
			ginertial = pushforward(gform, ginertial_m, nnod);

			lapgform = extractshelldisplacement(shell, lapddisp);                     /*{Xg+Ug}*/
			lapH = blockjacobimtx(lapgform, NULL, NULL, nnod);

			//epressure = assemshellpvct(shell, drccos);                		  /*{Pe}*/
			//gpressure = matrixvector(Tt, epressure, 6 * nnod);  				  /*{Pg}*/
			//volume += shellvolume(shell, drccos, area);                         /*VOLUME*/


			/*MID-POINT TRANSFORMATION MATRIX.*/
			midHPT = midpointmtx(HPT, lastHPT, alphaf, 6 * nnod);
			midTtPtHt = matrixtranspose(midHPT, 6 * nnod);

			/*MID-POINT FORCE VECTOR.*/
			mideinternal = midpointvct(einternal, lasteinternal, alphaf-xi, 6*nnod);/*xi : NUMERICAL DAMPING DISSIPATION*/
			midginternal = matrixvector(midTtPtHt, mideinternal, 6 * nnod);

			midginertial = midpointvct(ginertial, lastginertial, alpham, 6*nnod);

			//midgpressure = midpointvct(gpressure, lastgpressure, alphaf, 6*nnod);

			/*MID-POINT MASS & STIFFNESS MATRIX.*/

			Keff=assemtmtxCR_DYNA(eform, edisp, mideinternal, T, Ke, midTtPtHt, HPT,
								  ginertial, Me, R, lastRt, lapH,
								  alphaf, alpham, xi, beta, ddt, nnod);
			symmetricmtx(Keff, 6*nnod);
			assemgstiffnessIIwithDOFelimination(gmtx, Keff, &shell, constraintmain);

			/*GLOBAL VECTOR & MATRIX ASSEMBLY.*/
			for (ii = 0; ii < nnod; ii++)
			{
				for (jj = 0; jj < 6; jj++)
				{
					*(finertial + *(loffset + (6 * ii + jj))) += *(midginertial + 6 * ii + jj);
					*(finternal + *(loffset + (6 * ii + jj))) += *(midginternal + 6 * ii + jj);
					//*(fpressure + *(loffset + (6 * ii + jj))) += *(midgpressure + 6 * ii + jj);
				}
			}

			/*OUTPUT STRAIN ENERGY & STRESS.*/
			if (iteration == 1)
			{
				Ee = 0.0;
				Ep = 0.0;
				Eb = 0.0;
				for (ii = 0; ii < nnod; ii++)
				{
					for (jj = 0; jj < 2; jj++)
					{
						Ep += 0.5 * *(edisp + 6 * ii + jj) * *(einternal + 6 * ii + jj);
					}
					for (jj = 2; jj < 5; jj++)
					{
						Eb += 0.5 * *(edisp + 6 * ii + jj) * *(einternal + 6 * ii + jj);
					}
					Ee += 0.5 * *(edisp + 6 * ii + 5) * *(einternal + 6 * ii + 5);
				}
				Ee += Ep + Eb;
				fprintf(fene, "%5ld %e %e %e\n", shell.code, Ep, Eb, Ee);

				shellstress = matrixvector(DBe, edisp, 6 * nnod);
				for (ii = 0; ii < nnod; ii++)
				{
					for (jj = 0; jj < 6; jj++)
					{
						shell.stress[ii][jj] = *(shellstress + 6 * ii + jj);
					}
				}
				outputshellstress(shell, shellstress, fstr);
				free(shellstress);
			}

			/*MEMORY FREE : INITIAL CONFIGURATION.*/
			free(loffset);
			freematrix(drccosinit, 3);
			free(gforminit);
			free(eforminit);
			freematrix(DBe, 6 * nnod);
			freematrix(Ke, 6 * nnod);
			freematrix(Me, 6 * nnod);

			/*MEMORY FREE : DEFORMED CONFIGURATION OF LAST LAP.*/
			freematrix(lastdrccos, 3);
			freematrix(lastT, 6 * nnod);
			freematrix(lastTt, 6 * nnod);
			free(lastgform);
			free(lasteform);
			free(lastedisp);
			free(lasteinternal);
			freematrix(lastHPT, 6 * nnod);
			free(lastgacc_m);
			freematrix(lastR, 6 * nnod);
			freematrix(lastRt, 6 * nnod);
			free(lastginertial_m);
			free(lastginertial);

			/*MEMORY FREE : DEFORMED CONFIGURATION OF LAST ITERATION.*/
			freematrix(drccos, 3);
			freematrix(T, 6 * nnod);
			freematrix(Tt, 6 * nnod);
			free(gform);
			free(eform);
			free(edisp);
			free(einternal);
			freematrix(HPT, 6 * nnod);
			free(gacc_m);
			freematrix(R, 6 * nnod);
			free(ginertial_m);
			free(ginertial);

			free(lapgform);
			freematrix(lapH, 6 * nnod);

			/*MEMORY FREE : MID-POINT.*/
			freematrix(midHPT, 6 * nnod);
			freematrix(midTtPtHt, 6 * nnod);
			free(mideinternal);
			free(midginternal);
			free(midginertial);
			//free(midgpressure);
			freematrix(Keff, 6 * nnod);


			if (/*iteration==1 &&*/ (wdraw.childs + 1)->hdcC != NULL)/*DRAW DEFORMED ELEMENT OF LAST ITERATION.*/
			{
				drawglobalshell((wdraw.childs + 1)->hdcC,
					(wdraw.childs + 1)->vparam,
					*af, shell, 255, 255, 255,
					255, 255, 255, 0, ONSCREEN/*,i*/);
			}

		}
		/*DOF ELIMINATION.*/

		for (i = 0; i < msize; i++)
		{
			if (*(constraintmain + i) != i)
			{
				*(finertial + *(constraintmain + i)) += *(finertial + i);
				*(finertial + i) = 0.0;
				*(finternal + *(constraintmain + i)) += *(finternal + i);
				*(finternal + i) = 0.0;
				*(fexternal + *(constraintmain + i)) += *(fexternal + i);
				*(fexternal + i) = 0.0;
				//*(fpressure + *(constraintmain + i)) += *(fpressure + i);
				//*(fpressure + i) = 0.0;
			}
		}

		overlayhdc(*(wdraw.childs + 1), SRCPAINT);                  /*UPDATE DISPLAY.*/


		/*UNBALANCED FORCE & RESIDUAL AT MID-POINT.*/
		residual = 0.0;
		for (i = 0; i < msize; i++)
		{
			*(fexternal + i) = (alphaf*lastloadfactor + (1-alphaf)*loadfactor) * (*(fexternal + i)+*(fpressure + i));
			*(funbalance + i) = *(fexternal + i) - *(finertial + i) - *(finternal + i);
			/*funbalance : UNBALANCED FORCE -{E}.*/
			/*SIGN OF UNBALANCED FORCE IS INVERTED FROM DEFINITION.*/
			if ((confs + i)->iconf == 1)
			{
				*(freaction + i) = *(funbalance + i);
				*(funbalance + i) = 0.0;
			}
			residual += *(funbalance + i) * *(funbalance + i);
			*(gvct + i) = *(funbalance + i);
		}

		if ((outputmode == 0 && iteration == 1) || outputmode == 1)
		{
			outputdisp(ddisp, fdsp, nnode, nodes);                        /*FORMATION OUTPUT.*/
			outputdisp(finertial, finr, nnode, nodes);
			outputdisp(finternal, finf, nnode, nodes);
			outputdisp(fexternal, fexf, nnode, nodes);
			outputdisp(funbalance, fubf, nnode, nodes);
			outputdisp(freaction, frct, nnode, nodes);
		}

		/*CROUT LU DECOMPOSITION.*/
		nline = croutlu(gmtx, confs, msize, &det, &sign, gcomp1);
		sprintf(string, "LAP: %4d ITER: %2d {LOAD}= % 5.8f {RESD}= %1.6e {DET}= %8.5f {SIGN}= %2.0f {BCL}= %1d {EPS}=%1.5e {V}= %8.5f\n",
			nlap, iteration, loadfactor, residual, det, sign, 0, 0.0, volume);
		fprintf(ffig, "%s", string);
		errormessage(string);
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

		/*lapddisp : INCREMENTAL TRANSITION & ROTATION IN THIS LAP.*/
		if(iteration==1)
		{
			for(ii = 0; ii < msize; ii++)
			{
				*(lapddisp + ii)  =   0.0;
				*(udinit_m + ii)  =- (gamma / beta - 1.0) **(lastud_m + ii)
								   - (gamma / (2.0 * beta) - 1.0) * ddt **(lastudd_m + ii);
				*(uddinit_m + ii) =- (1.0 / (beta * ddt)) **(lastud_m + ii)
								   - (1.0 / (2.0 * beta) - 1.0) **(lastudd_m + ii);
			}
		}
		/*INCREMENT DISPLACEMENT OF THIS LAP IN SPATIAL FORM.*/
		updaterotation(lapddisp,gvct,nnode);

		/*INCREMENTAL DISPLACEMENT OF THIS LAP IN MATERIAL FORM.*/
		lapddisp_m = pullback(lastddisp,lapddisp,nnode);

		/*UPDATE VELOCITY & ACCELERATION IN MATERIAL FORM.*/
		for(ii = 0; ii < msize; ii++)
		{
			*(ud_m + ii)  = (gamma / (beta * ddt)) **(lapddisp_m + ii) + *(udinit_m + ii);
			*(udd_m + ii) = (1.0 / (beta * ddt * ddt)) **(lapddisp_m + ii) + *(uddinit_m + ii);
		}

		/*UPDATE POSITION & VELOCITY & ACCELERATION IN SPATIAL FORM.*/
		updaterotation(ddisp, gvct, nnode);

		//ud  = pushforward(ddisp,ud_m,nnode);
		//udd = pushforward(ddisp,udd_m,nnode);

		for (ii = 0; ii < msize; ii++)
		{
			if (*(constraintmain + ii) != ii)
			{
				mainoff = *(constraintmain + ii);
				*(ddisp + ii) = *(ddisp + mainoff);
				*(ud_m + ii)  = *(ud_m + mainoff);
				*(udd_m + ii) = *(udd_m + mainoff);
			}
		}
		free(lapddisp_m);

		if ((residual<tolerance || iteration>maxiteration) && iteration != 1)
		{
			nlap++;
			iteration = 0;
		}
		iteration++;
	}






	if ((wdraw.childs + 1)->hdcC != NULL && ddisp != NULL)	/*DRAW LAST FRAME.*/
	{
		for (i = 1; i <= nelem; i++)
		{
			inputelem(elems, melem, i - 1, &elem);

			for (ii = 0; ii < elem.nnod; ii++)
			{
				inputnode(ddisp, elem.node[ii]);
			}
			if (globaldrawflag == 1)
			{
				drawglobalwire((wdraw.childs + 1)->hdcC,
					(wdraw.childs + 1)->vparam,
					*af, elem, 255, 255, 255,
					255, 255, 255, 0, ONSCREEN/*,i*/);
			}
		}
		for (i = 1; i <= nshell; i++)
		{
			shell = *(shells + i - 1);

			for (ii = 0; ii < shell.nnod; ii++)
			{
				inputnode(ddisp, shell.node[ii]);
			}
			if (globaldrawflag == 1)
			{
				drawglobalshell((wdraw.childs + 1)->hdcC,
					(wdraw.childs + 1)->vparam,
					*af, shell, 255, 255, 255,
					255, 255, 255, 0, ONSCREEN);
			}
		}
		overlayhdc(*(wdraw.childs + 1), SRCPAINT);     	  	/*UPDATE DISPLAY.*/
	}

	fclose(fin);			/*FILE CLOSE.*/
	fclose(fout);
	fclose(fonl);
	fclose(ffig);
	fclose(fbcl);
	fclose(fene);

	gfree(gmtx, nnode);  	/*FREE GLOBAL MATRIX.*/
	free(fexternal);		/*FREE VECTOR*/
	free(finternal);

	errormessage(" ");
	errormessage("COMPLETED.");
	return 0;
}



