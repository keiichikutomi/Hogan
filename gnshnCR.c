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
	rvct = (double*)malloc(3 * sizeof(double));
	for (n = 0; n < nnode; i++)
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

double* pushforward(double* ddisp, double* gvct_m, int nnode)
{
	int i,n;
	double* rvct, * vct_s, * vct_m, * gvct_s;
	double** rmtx;

	gvct_s = (double*)malloc(6*nnode * sizeof(double));
	rvct = (double*)malloc(3 * sizeof(double));
	for (n = 0; n < nnode; i++)
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


/*MID-POINT VARIABLES.*/
double* midpointvct(double* vct,double* lastvct,double alpha,int size)
{
	int i;
	double* midvct;

	midvct= (double*)malloc(size * sizeof(double));
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

int gnshnCR(struct arclmframe* af)
{
	DWORDLONG memory;
	FILE *fdata, *fin, *fonl, * fdsp, * fexf, * finf, * fubf, * fstr, * fene, * ffig, * fbcl, * feig, * flog;         /*FILE 8 BYTES*/
	char dir[] = DIRECTORY;
	char s[80], string[400], inpname[50], fname[50];


	int i, ii, jj;
	int nnode, nelem, nshell, nsect, nreact, nconstraint;
	int nlap, laps;  /*LAP COUNT*/

	long int msize;

	double ddt = 0.00020;/*TIME INCREMENT[sec]*/
	double time = 0;/*TOTAL TIME[sec]*/


	/***FOR ELEMENT***/
	double* gforminit, * gform;                      /*GLOBAL COORDINATION OF ELEMENT*/
	double* eforminit, * eform;                      /*LOCAL COORDINATION OF ELEMENT*/
	double* edisp;                                   /*LOCAL DEFORMATION OF ELEMENT*/

	double* midgdisp, * midedisp;
	double* lastgdisp, * lastedisp;

	double* ginternal, * einternal;
	double* midginternal, * mideinternal;
	double* lastginternal, *lasteinternal;

	double* gexternal, * eexternal;
	double* midgexternal, * mideexternal;
	double* lastgexternal, *lasteexternal;

	double* gacc,* ginertial;
	double* gvel,* gdamping;
	double* shellstress;                           /*σx,σy,τxy,Mx,My,Mxy OF ELEMENT*/
	double Ep, Eb, Ee;                             /*STRAIN ENERGY OF ELEMENT*/
	double** Me, ** Ke, ** Kt, ** DBe;

	double** drccos,** drccosinit, ** T, ** Tt;                           /*ELEMENT COORDINATION MATRIX*/
	double** lastdrccos, ** lastT, ** lastTt;
	double** middrccos, ** midT, ** midTt;

	int* loffset;

	double area, volumetotal;
	double masstotal,massdiag;

	/***FOR INCREMENTAL***/
	int iteration;
	int maxiteration = 20;
	double residual;
	double tolerance = 1.0E-8;

	double momentumlinear,momentumangular;

	long int nline;
	double determinant, sign;
	double loadfactor = 0.0;
	double lambda = 0.0;

	/*ENERGY MOMENTUM METHOD'S PARAMETER.*/
	double rho = 1.0;
	double alpham = 2.0*rho-1/(rho+1);  /*MID-POINT USED TO EVALUATE INERTIAL FORCE*/
	double alphaf = rho/(rho+1);        /*MID-POINT USED TO EVALUATE INTERNAL FORCE*/
	double xi = 0.0;   	 				/*NUMERICAL DISSIPATION COEFFICIENT*/

	/*NEWMARK-BETA'S PARAMETER.*/
	double beta = pow(1-alpham+alphaf,2)/4.0;
	double gamma = (0.5-alpham+alphaf)/4.0;

	/*FOR READING ANALYSISDATA*/
	/*ANALYSIS SETTING*/
	int nstr, pstr, readflag;
	char **data, *filename;

	int fnode=NULL,fnodedrc=NULL;/*NODE DEFORMATION AT ANALYSIS ENDING*/
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
	snprintf(fname, sizeof(fname), "%s.%s", inpname, "str");
	fstr = fopen(fname, "w");
	snprintf(fname, sizeof(fname), "%s.%s", inpname, "onl");
	fonl = fopen(fname, "w");
	snprintf(fname, sizeof(fname), "%s.%s", inpname, "fig");
	ffig = fopen(fname, "w");
	snprintf(fname, sizeof(fname), "%s.%s", inpname, "bcl");
	fbcl = fopen(fname, "w");
	snprintf(fname, sizeof(fname), "%s.%s", inpname, "ene");
	fene = fopen(fname, "w");
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

	///POSITION VECTOR INITIALIZATION///
	/*IN SPATIAL FORM*/
	double* ddisp,* iform,* lastddisp,* lapddisp;

	iform = (double*)malloc(msize * sizeof(double));		/*INITIAL*/
	ddisp = (double*)malloc(msize * sizeof(double));		/*LATEST ITERATION*/
	lastddisp = (double*)malloc(msize * sizeof(double));	/*LATEST LAP*/
	lapddisp = (double*)malloc(msize * sizeof(double));		/*INCREMENT IN THE LAP*/


	///VELOSITY & ACCELERATION VECTOR INITIALIZATION///
	/*ud:VELOSITY, udd:ACCELERATION*/
	/*ROTATIONAL DOFs ARE REPRESENTED IN SPATIAL & MATERIAL FORM*/
	double* udinit,* uddinit,* ud,* udd,* ud_m,* udd_m,* lastud,* lastudd,* lastud_m,* lastudd_m;

	udinit = (double*)malloc(msize * sizeof(double));  	   	/*NEWMARK INITIAL IN LAP*/
	uddinit = (double*)malloc(msize * sizeof(double));

	ud = (double*)malloc(msize * sizeof(double)); 			/*LATEST ITERATION IN SPATIAL*/
	udd = (double*)malloc(msize * sizeof(double));

	ud_m = (double*)malloc(msize * sizeof(double));  		/*LATEST ITERATION IN MATERIAL*/
	udd_m = (double*)malloc(msize * sizeof(double));

	lastud = (double*)malloc(msize * sizeof(double));  		/*LATEST LAP IN SPATIAL*/
	lastudd = (double*)malloc(msize * sizeof(double));

	lastud_m = (double*)malloc(msize * sizeof(double));     /*LATEST LAP IN MATERIAL*/
	lastudd_m = (double*)malloc(msize * sizeof(double));

	/*ONLY FOR GENERALIZED ALPHA ALGOLITHM*/
	double* midud,* midudd,* midud_m,* midudd_m;

	midud = (double*)malloc(msize * sizeof(double));     	 /*MID-POINT IN SPATIAL*/
	midudd = (double*)malloc(msize * sizeof(double));

	midud_m = (double*)malloc(msize * sizeof(double)); 	 /*MID-POINT IN MATERIAL*/
	midudd_m = (double*)malloc(msize * sizeof(double));


	///FORCE VECTOR INITIALIZATION///
	double* fexternal, * finternal,* finertial,* fdamping,* funbalance;

	fexternal = (double*)malloc(msize * sizeof(double));         /*BASE EXTERNAL FORCE VECTOR.*/
	finternal = (double*)malloc(msize * sizeof(double));         /*INTERNAL FORCE VECTOR.*/
	finertial = (double*)malloc(msize * sizeof(double));         /*INERTIAL FORCE VECTOR.*/
	//fdamping= (double*)malloc(msize * sizeof(double));         /*DAMPING FORCE VECTOR.*/
	funbalance = (double*)malloc(msize * sizeof(double));        /*UNBALANCED FORCE VECTOR.*/

	for (i = 0; i < msize; i++)
	{
		/*POSITION INCREMENT*/
		*(lapddisp + i) = 0.0;

		/*VELOCITY & ACCELERATION*/
		*(udinit + i) = 0.0;
		*(uddinit + i) = 0.0;

		*(ud + i) = 0.0;
		*(udd + i) = 0.0;

		*(ud_m + i) = 0.0;
		*(udd_m + i) = 0.0;

		*(lastud + i) = 0.0;
		*(lastudd + i) = 0.0;

		*(lastud_m + i) = 0.0;
		*(lastudd_m + i) = 0.0;

		*(midud + i) = 0.0;
		*(midudd + i) = 0.0;

		*(midud_m + i) = 0.0;
		*(midudd_m + i) = 0.0;

		/*FORCE*/
		*(finertial + i) = 0.0;
		*(finternal + i) = 0.0;
		*(fexternal + i) = 0.0;
		*(funbalance + i) = 0.0;

		*(constraintmain + i) = i;
	}


	///MATRIX INITIALIZATION///
	struct gcomponent ginit = { 0,0,0.0,NULL };
	struct gcomponent* gmtx, * g, * p;/*GLOBAL MATRIX*/
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
	initialformCR(ninit, lastddisp, nnode);           /*ASSEMBLAGE FORMATION.*/
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


	///FOR DRAWING 1///
	GetAsyncKeyState(VK_LBUTTON);                   /*CLEAR KEY LEFT.*/
	GetAsyncKeyState(VK_RBUTTON);                  /*CLEAR KEY RIGHT.*/
	if(globaldrawflag==1)
	{
	  drawglobalaxis((wdraw.childs+1)->hdcC,(wdraw.childs+1)->vparam,0,0,255);                     /*DRAW GLOBAL AXIS.*/
	}
	///FOR DRAWING 1///


	nlap = 1;
	iteration = 1;
	residual = 0.0;

	while (nlap <= laps)
	{
		af->nlaps = nlap;

		///FOR DRAWING 2///
		setincrement((wmenu.childs+2)->hwnd,laps,nlap,maxiteration,iteration);
		if(iteration==1)clearwindow(*(wdraw.childs+1));
		///FOR DRAWING 2///


		if (iteration == 1)
		{
			sprintf(string, "LAP:%5ld/%5ld", nlap, laps);
			errormessage(string);


			for (i = 0; i < msize; i++)
			{
				*(lastddisp + i) = *(ddisp + i);
				*(lapddisp + i) = 0.0;
				*(lastud + i) = *(ud + i);
				*(lastudd + i) = *(udd + i);
			}


		}
		clearwindow(*(wdraw.childs + 1));




		///////////*INITIALIZATION.*///////////
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
		/*FORCE VECTOR INITIALIZATION.*/
		for (i = 0; i < msize; i++)
		{
			*(finertial  + i) = 0.0;
			*(finternal + i) = 0.0;
			*(fexternal + i) = 0.0;
			*(funbalance + i) = 0.0;
		}

		loadfactor = lambda * nlap;

		volume = 0.0;

		assemconf(confs, fexternal, 1.0, nnode);               /*GLOBAL VECTOR.*/



		///////////*ASSEMBLAGE MATRIX.*///////////

		for (i = 1; i <= nshell; i++)
		{
			inputshell(shells, mshell, i - 1, &shell);
			shell.sect = (shells + i - 1)->sect;                      /*READ SECTION DATA.*/
			loffset = (int*)malloc(6 * shell.nnod * sizeof(int));
			for (ii = 0; ii < shell.nnod; ii++)
			{
				for (jj = 0; jj < 6; jj++)
				{
					*(loffset + (6 * ii + jj)) = 6 * (shell.node[ii]->loff) + jj;
				}
			}

			//////////*INITIAL*//////////

			for (ii = 0; ii < shell.nnod; ii++)
			{
				inputnode(iform, shell.node[ii]);
			}
			drccosinit = shelldrccos(shell, &area);
			gforminit = extractshelldisplacement(shell, iform);                 /*{Xg}*/
			eforminit = extractlocalcoord(gforminit,drccosinit,shell.nnod);     /*{Xe}*/

			if (iteration == 1)
			{
				DBe = (double**)malloc(18 * sizeof(double*));
				for (ii = 0; ii < 18; ii++)
				{
					*(DBe + ii) = (double*)malloc(18 * sizeof(double));
					for (jj = 0; jj < 18; jj++)
					{
						*(*(DBe + ii) + jj) = 0.0;
					}
				}
			}
			else
			{
				DBe = NULL;
			}
			Ke = assemshellemtx(shell, drccosinit, DBe);   						/*ELASTIC MATRIX[Ke]*/
			Me = assemshellmmtx(shell, drccosinit);          					/*INERTIAL MATRIX[Me]*/


			//////////*LATEST ITERATION*//////////

			for (ii = 0; ii < shell.nnod; ii++)
			{
				inputnode(ddisp, shell.node[ii]);
			}
			drccos = shelldrccos(shell, &area);
			gform = extractshelldisplacement(shell, ddisp);                     /*{Xg+Ug}*/
			eform = extractlocalcoord(gform,drccos,shell.nnod);                 /*{Xe+Ue}*/

			edisp = extractdeformation(eforminit, eform, shell.nnod);           /*{Ue}*/
			einternal = matrixvector(Ke, edisp, 6 * shell.nnod);      			/*ELEMENT INTERNAL FORCE{Fe}=[Ke]{Ue}.*/

			T = transmatrixIII(drccos, shell.nnod);         					/*TRANSFORMATION MATRIX[T]*/
			Kt = assemtmtx(Ke, eform, edisp, einternal, T, shell.nnod);         /*TANGENTIAL STIFFNESS MATRIX[Kt].*//*PROJECTION of estress[Pt][Ht]{Fe}.*/

			gacc_m = extractshelldisplacement(shell, udd_m);
			ginertial_m = matrixvector(Me, gacc, 6 * shell.nnod);
			ginertial = pushforward(edisp, ginertial_m, 3);

			epressure = assemshellpvct(shell, drccos);          /*ELEMENT EXTERNAL FORCE{Fe}.*/
			gpressure = matrixvector(Tt, epressure, 6 * shell.nnod); /*GLOBAL EXTERNAL FORCE.*/

			volume += shellvolume(shell, drccos, area);                         /*VOLUME*/


			if ((outputmode == 0 && (iteration == 1 || BCLFLAG == 2)) || outputmode == 1)
			{
				Ee = 0.0;
				Ep = 0.0;
				Eb = 0.0;
				for (ii = 0; ii < shell.nnod; ii++)                   /*UPDATE STRAIN ENERGY.*/
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

				shellstress = matrixvector(DBe, edisp, 6 * shell.nnod);
				for (ii = 0; ii < shell.nnod; ii++)                          /*UPDATE STRESS.*/
				{
					for (jj = 0; jj < 6; jj++)
					{
						shell.stress[ii][jj] = *(shellstress + 6 * ii + jj);
					}
				}
				outputshellstress(shell, shellstress, fstr);
				free(shellstress);
				freematrix(DBe, 6 * shell.nnod);
			}

			//////////*MID-POINT*//////////
			midginertial = midpointvct(ginertial, lastginertial, alpham   , 6*shell.nnod);

			mideinternal = midpointvct(einternal, lasteinternal, alphaf-xi, 6*shell.nnod);/*xi : NUMERICAL DAMPING DISSIPATION*/
			midA = midpointmtx(A, lastA, alpha - xi, 6 * shell.nnod);
			midginternal = matrixvector(midA, mideinternal, 6 * shell.nnod);

			midgexternal = midpointvct(gexternal, lastgexternal, alphaf   , 6*shell.nnod);


			//////////*GLOBAL VECTOR & MATRIX ASSEMBLY*//////////
			for (ii = 0; ii < shell.nnod; ii++)
			{
				for (jj = 0; jj < 6; jj++)
				{
					*(finertial + *(loffset + (6 * ii + jj))) += *(midginertial + 6 * ii + jj);
					*(finternal + *(loffset + (6 * ii + jj))) += *(midginternal + 6 * ii + jj);
					*(fexternal + *(loffset + (6 * ii + jj))) += *(midgexternal + 6 * ii + jj);
				}
			}
			assemgstiffnessIIwithDOFelimination(gmtx, K, &shell, constraintmain);

			freematrix(drccosinit, 3);
			freematrix(lastdrccos, 3);
			freematrix(drccos, 3);

			freematrix(T, 18);
			freematrix(Tt, 18);
			freematrix(Ke, 18);
			freematrix(Kt, 18);
			freematrix(Me, 18);

			free(einternal);
			free(ginternal);
			free(eexternal);
			free(gexternal);

			free(eform);
			free(gform);
			free(lastedisp);
			free(lastgdisp);
			free(edisp);
			free(gdisp);


			if (/*iteration==1 &&*/ (wdraw.childs + 1)->hdcC != NULL)/*DRAW DEFORMED ELEMENT.*/
			{
				drawglobalshell((wdraw.childs + 1)->hdcC,
					(wdraw.childs + 1)->vparam,
					*af, shell, 255, 255, 255,
					255, 255, 255, 0, ONSCREEN/*,i*/);
			}

		}

		///////////*DOF ELIMINATION.*///////////

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
			}
		}

		overlayhdc(*(wdraw.childs + 1), SRCPAINT);                  /*UPDATE DISPLAY.*/


		//////////*UNBALANCED FORCE & RESIDUAL AT MID-POINT.*///////////
		residual = 0.0;
		for (i = 0; i < msize; i++)
		{
			*(funbalance + i) = loadfactor * *(fexternal + i) - *(finternal + i) - *(finertial + i);
			/*SIGN OF UNBALANCED FORCE IS INVERTED FROM DEFINITION.*/
			if ((confs + i)->iconf == 1) *(funbalance + i) = 0.0;
			residual += *(funbalance + i) * *(funbalance + i);
		}

		/*CROUT LU DECOMPOSITION.*/
		nline = croutludecomposition(gmtx, gvct, confs, msize, &det, &sign);/*INCREMENT OF ITERATION IN SPATIAL FORM*/
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
			/*FILE CLOSE*/
			/*MEMORY FREE*/
			return 1;
		}


		/*lapddisp : INCREMENTAL TRANSITION & ROTATION IN THIS LAP.*/
		if(iteration==1)
		{
			for(ii = 0; ii < msize; ii++)
			{
				*(lapddisp + ii)=0.0;
				*(udinit_m + ii)  =- (gamma / beta - 1.0) **(lastud_m + ii)
								   - (gamma / (2.0 * beta) - 1.0) * ddt **(lastudd_m + ii);
				*(uddinit_m + ii) =- (1.0 / (beta * ddt)) **(lastud_m + ii)
								   - (1.0 / (2.0 * beta) - 1.0) **(lastudd_m + ii);
			}
		}
		/*INCREMENT DISPLACEMENT OF THIS LAP IN SPATIAL FORM*/
		updaterotation(lapddisp,gvct,nnode);

		/*INCREMENTAL DISPLACEMENT OF THIS LAP IN MATERIAL FORM*/
		lapddisp_m = pullback(lapddisp,lastddisp,nnode);

		/*UPDATE VELOCITY & ACCELERATION IN MATERIAL FORM*/
		for(ii = 0; ii < msize; ii++)
		{
			*(ud_m + ii)  = (gamma / (beta * ddt)) **(lapddisp_m + ii) + *(udinit_m + ii);
			*(udd_m + ii) = (1.0 / (beta * ddt * ddt)) **(lapddisp_m + ii) + *(uddinit_m + ii);
		}

		/*UPDATE POSITION & VELOCITY & ACCELERATION IN SPATIAL FORM*/
		updaterotation(ddisp, lapdisp, nnode);
		ud  = pushforward(ud_m,ddisp,nnode);
		udd = pushforward(udd_m,ddisp,nnode);


		for (ii = 0; ii < msize; ii++)
		{
			if (*(constraintmain + ii) != ii)
			{
				mainoff = *(constraintmain + ii);
				*(ddisp + ii) = *(ddisp + mainoff);
				*(ud + ii) = *(ud + mainoff);
				*(udd + ii) = *(udd + mainoff);
			}
		}

		if ((residual<tolerance || iteration>maxiteration) && iteration != 1)
		{
			nlap++;
			iteration = 0;
			time += ddt;
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



