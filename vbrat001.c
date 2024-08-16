int vbrat001(struct arclmframe* af)
{
	DWORDLONG memory;
	char dir[]=DIRECTORY;
	FILE *fdata, *fin, *fonl, * fdsp, * fexf, * finf, * fubf, * frct, * fstr, * fene, * ffig, * fbcl, * feig, * fout;         /*FILE 8 BYTES*/
	char s[80], string[400], inpname[50], fname[50];

	int i, ii, jj;
	int nnode, nelem, nshell, nsect, nreact, nconstraint;
	int nlap, laps;

	long int msize, nline;
	double time;

	/*FOR READING ANALYSISDATA*/
	int nstr, pstr, readflag, node;
	char **data, *filename;

	/***GLOBAL MATRIX***/
	struct gcomponent ginit = { 0,0,0.0,NULL };
	struct gcomponent* gmtx, * g, * p, * gcomp1;/*GLOBAL MATRIX*/
	double gg;


	/***GLOBAL VECTOR*/
	double* gvct;

	/***FOR EACH ELEMENT***/
	int nnod;
	double* gforminit, * eforminit;                      /*DEFORMED COORDINATION OF ELEMENT*/
	double* gform, * eform;                      /*INITIAL COORDINATION OF ELEMENT*/
	double* edisp;
	double* ginternal, * einternal;                        /*INTERNAL FORCE OF ELEMENT*/
	double* gexternal, * eexternal;                        /*EXTERNAL FORCE OF ELEMENT*/
	double* gpressure, * epressure;                          /*EXTERNAL FORCE OF ELEMENT*/
	double* shellstress;                        /*σx,σy,τxy,Mx,My,Mxy OF ELEMENT*/
	double Ep, Eb, Ee;                           /*STRAIN ENERGY OF ELEMENT*/
	double** Me,** Ke,** Kt,** DBe,** drccos,** drccosinit;                           /*MATRIX*/
	double** T,** Tt,** HP,**PtHt,**HPT,**TtPtHt;
	int* loffset;
	double area;
	/***FOR ARC-LENGTH INCREMENTAL***/
	int iteration;
	int maxiteration = 20;
	int maxtime = 20;
	double ddt = 1.0E-8;
	double residual;
	double tolerance = 1.0e-4;

	double determinant, sign;
	double arcsum, predictorsign;
	double lambda, loadfactor = 0.0;
	double gamma;
	double dupdue, dupdup;
	double volume, masstotal, massdiag;
	/*ARCLENGTH CONTROL*/
	double k1, k, scaledarclength;
	/***FOR BISECSYLVESTER & EXTENDED SYSTEM***/
	double biseceps = 1e-5;
	double gradeps = 1e-8;
	double LL;
	double LR;
	double LM;
	/*FOR BUCKLING DETECTION*/
	double lastsign;
	double* lastpivot;
	double nextloadfactor;
	double lastloadfactor, lastlambda, bisecloadfactor, biseclambda;
	double* lastddisp, * lastgvct, * bisecddisp, * bisecgvct;

	double laploadfactor;
	double* lapddisp;
	double* nextddisp;

	double* epsddisp, * epsgvct, * epsfunbalance, * epsfinternal, * epsfexternal, * epsfpressure, * re, * rp;



	double eigen, lasteigen;
	double* lastevct, * evct;
	int inverseiter;
	double eps;

	double evctdot, len, evctlastevct, evctevct, evctfunbalance, epsevctfunbalance;


	int nmode = 0;
	double* norm, * dm;	/*NORM & PIVOT OF LDL-MODE*/
	int* m;            	/*LINE OF LDL-MODE*/
	double** ldlevct;/*NORMALIZED EIGEN MODE*/


	int BCLFLAG = 0;/*FLAG FOR BUCKLING DETECTION*/
	int ENDFLAG = 0;/*FLAG FOR ANLYSIS TERMINATION*/
	int UNLOADFLAG = 0;/*FLAG FOR ANLYSIS TERMINATION*/

	int SCALINGARCFLAG = 0;
	int BIGININGARCFLAG = 0;
	double biginingarcratio = 1.0;
	int biginingarclap = 0;









	 /*
	FILE *fgetstofopenII(const char *directory,const char *mode,const char *filename)
	{
	  FILE *f=NULL;
	  char fname[256],dandf[256];

		strcpy(fname,filename);
		strcpy(dandf,directory);
		strcat(dandf,fname);
		f=fopen(dandf,mode);

	  return f;
	}*/


	fdata = fgetstofopenII(dir, "r", "analysisdata.txt");
	if (fdata == NULL)
	{
		printf("couldn't open analysisdata.txt\n");
		getchar();
		exit(EXIT_FAILURE);
	}
	readflag = 1;



	sprintf(string, "FILENAME:%s\n LAPS = %d\n MAX ITERATION= %d\n ARCLENGTH = %lf\n", filename, laps, maxiteration, arclength);
	errormessage(string);
	if (outputmode == 0)printf("OUTPUT CONVERGED RESULT\n");
	if (outputmode == 1)printf("OUTPUT ALL RESULT\n");

	strcat(filename,".inl");




	sprintf(string, "INITIAL:");
	memory = availablephysicalmemoryEx(string);   /*MEMORY AVAILABLE*/

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
	snprintf(fname, sizeof(fname), "%s.%s", inpname, "otl");
	fout = fopen(fname, "w");


	clock_t t0;
	t0 = clock();                                                   /*CLOCK BEGIN.*/

	inputinitII(fin, &nnode, &nelem, &nshell, &nsect, &nconstraint); /*INPUT INITIAL.*/

	msize = 6 * nnode;                                      /*SIZE OF GLOBAL MATRIX.*/

	struct osect* sects;
	struct onode* nodes;
	struct onode* ninit;
	struct owire elem;
	struct owire* elems;
	struct oshell shell;
	struct oshell* shells;
	struct oconf* confs;
	sects = (struct osect*)malloc(nsect * sizeof(struct osect));
	nodes = (struct onode*)malloc(nnode * sizeof(struct onode));
	ninit = (struct onode*)malloc(nnode * sizeof(struct onode));
	elems = (struct owire*)malloc(nelem * sizeof(struct owire));
	shells = (struct oshell*)malloc(nshell * sizeof(struct oshell));
	confs = (struct oconf*)malloc(msize * sizeof(struct oconf));

	double* ddisp, * iform;
	iform = (double*)malloc(msize * sizeof(double));          /*INITIAL GLOBAL VECTOR.*/
	ddisp = (double*)malloc(msize * sizeof(double));

    struct memoryelem* melem;
	struct memoryshell* mshell;
	melem = (struct memoryelem*)malloc(nelem * sizeof(struct memoryelem));
	mshell = (struct memoryshell*)malloc(nshell * sizeof(struct memoryshell));

    long int mainoff;
	long int* constraintmain;
	constraintmain = (long int*)malloc(msize * sizeof(long int));


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

	gmtx = (struct gcomponent*)malloc(msize * sizeof(struct gcomponent));/*DIAGONALS OF GLOBAL MATRIX.*/
	gvct = (double *)malloc(msize * sizeof(double));/*INCREMENTAL GLOBAL VECTOR.*/

	double* funbalance, * fexternal, * finternal, * freaction, * fpressure, * fbaseload, * fgivendisp, * fswitching;
	funbalance = (double*)malloc(msize * sizeof(double));          /*UNBALANCED INTERNAL FORCE VECTOR.*/
	freaction = (double*)malloc(msize * sizeof(double));           /*EXTERNAL FORCE VECTOR.*/
	fexternal = (double*)malloc(msize * sizeof(double));           /*EXTERNAL FORCE VECTOR.*/
	finternal = (double*)malloc(msize * sizeof(double));           /*INTERNAL FORCE VECTOR.*/
	fpressure = (double*)malloc(msize * sizeof(double));           /*PRESSURE VECTOR.*/
	fbaseload = (double*)malloc(msize * sizeof(double));           /*BASE LOAD VECTOR.*/
	fgivendisp = (double*)malloc(msize * sizeof(double));           /*BASE LOAD VECTOR.*/
	fswitching = (double*)malloc(msize * sizeof(double));





	evct = (double*)malloc(msize * sizeof(double));
	lastevct = (double*)malloc(msize * sizeof(double));

	nextddisp = (double*)malloc(msize * sizeof(double));

	bisecddisp = (double*)malloc(6 * nnode * sizeof(double));
	bisecgvct = (double*)malloc(6 * nnode * sizeof(double));



	for (i = 0; i < msize; i++)
	{
		(gmtx + i)->down = NULL;
		*(gvct + i) = 0.0;
		*(constraintmain + i) = i;
	}

	inputtexttomemory(fin, af);
	fclose(fin);

	nnode = af->nnode;
	ninit = af->ninit;
	nelem = af->nelem;
	nshell = af->nshell;
	nsect = af->nsect;
	nreact = af->nreact;
	nconstraint = af->nconstraint;

	initialformCR(nodes, ddisp, nnode);           /*ASSEMBLAGE FORMATION.*/
	initialformCR(ninit, iform, nnode);           /*ASSEMBLAGE FORMATION.*/
	initialelem(elems, melem, nelem);             /*ASSEMBLAGE ELEMENTS.*/
	initialshell(shells, mshell, nshell);         /*ASSEMBLAGE ELEMENTS.*/


	for (ii = 0; ii < msize; ii++)
	{
		if (*(constraintmain + ii) != ii)
		{
			(confs + ii)->iconf = (signed char)1;
		}
	}

	setviewpoint((wdraw.childs+0)->hwnd,arc,
						 &((wdraw.childs+1)->vparam));
	setviewparam((wmenu.childs+2)->hwnd,
						 (wdraw.childs+1)->vparam);
	clearwindow(*(wdraw.childs+1));
	drawarclmframe((wdraw.childs+1)->hdcC,
						   (wdraw.childs+1)->vparam,arc,0,ONSCREEN);


	if(globaldrawflag==1)
	{
	  drawglobalaxis((wdraw.childs+1)->hdcC,(wdraw.childs+1)->vparam,0,0,255);                     /*DRAW GLOBAL AXIS.*/
	}
	///FOR DRAWING 1///


	nlap = 1;
	iteration = 1;
	residual = 0.0;

	assemconf(confs,fbaseload,1.0,nnode);               /*GLOBAL VECTOR.*/
	assemgivend(confs,fbaseload,1.0,nnode);


	while (nlap <= laps && ENDFLAG == 0)
	{

		sprintf(string, "LAP: %5ld / %5ld ITERATION: %5ld\n", nlap, laps, iteration);
		af->nlaps = nlap;

		///FOR DRAWING 2///
		//setincrement((wmenu.childs+2)->hwnd,laps,nlap,maxiteration,iteration);
		if(iteration==1)clearwindow(*(wdraw.childs+1));
		///FOR DRAWING 2///

		if ((outputmode == 0 && (iteration == 1 || BCLFLAG == 2)) || outputmode == 1)
		{
			fprintf(fdsp, string);
			fprintf(finf, string);
			fprintf(fexf, string);
			fprintf(fubf, string);
			fprintf(frct, string);
			fprintf(fstr, string);
			fprintf(fene, string);
		}

		for (i = 1; i <= msize; i++)/*MATRIX INITIALIZATION.*/
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

			*(finternal + (i - 1)) = 0.0;			 /*GLOBAL VECTOR INITIALIZATION.*/
			*(fexternal + (i - 1)) = 0.0;
			*(fpressure + (i - 1)) = 0.0;
			*(funbalance + (i - 1)) = 0.0;
			*(freaction+ (i - 1)) = 0.0;
		}
		comps = msize; /*INITIAL COMPONENTS=DIAGONALS.*/

		volume = 0.0;


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
			assemgstiffnesswithDOFelimination(gmtx, Kt, &elem, constraintmain); /*ASSEMBLAGE TANGENTIAL STIFFNESS MATRIX.*/

			for (ii = 0; ii < nnod; ii++)
			{
				for (jj = 0; jj < 6; jj++)
				{
					*(finternal + *(loffset + (6 * ii + jj))) += *(ginternal + 6 * ii + jj);
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

			///FOR DRAWING 3///
			if(iteration==1 && (wdraw.childs+1)->hdcC!=NULL)
			{
			  drawglobalwire((wdraw.childs+1)->hdcC,
							 (wdraw.childs+1)->vparam,
							  *af,elem,255,255,255,
									   255,255,255,0,ONSCREEN);
			}

			///FOR DRAWING 3///
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
			drccosinit = shelldrccos(shell, &area);
			gforminit = extractshelldisplacement(shell, iform);                 /*{Xg}*/
			eforminit = extractlocalcoord(gforminit,drccosinit,nnod);        	/*{Xe}*/

			DBe = (double**)malloc(6 * nnod * sizeof(double*));
			for (ii = 0; ii < 6 * nnod; ii++)
			{
				*(DBe + ii) = (double*)malloc(6 * nnod * sizeof(double));
				for (jj = 0; jj < 6 * nnod; jj++)
				{
					*(*(DBe + ii) + jj) = 0.0;
				}
			}
			Ke = assemshellemtx(shell, drccosinit, DBe);                        /*[Ke].*/

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
			einternal = matrixvector(Ke, edisp, 6 * nnod);      				/*{Fe}=[Ke]{Ue}.*/
			epressure = assemshellpvct(shell, drccos);                			/*{Pe}.*/
			volume += shellvolume(shell, drccos, area);                   		/*VOLUME*/

			ginternal = (double*)malloc(6 * nnod * sizeof(double));
			HPT = (double**)malloc(6 * nnod * sizeof(double*));
			for (ii = 0; ii < 6 * nnod; ii++)
			{
				*(HPT + ii) = (double*)malloc(6 * nnod * sizeof(double));
			}
			Kt = assemtmtxCR(Ke, eform, edisp, einternal, ginternal, T, HPT, nnod);	/*TANGENTIAL MATRIX[Kt].*/
			symmetricmtx(Kt, 6*nnod);											/*SYMMETRIC TANGENTIAL MATRIX[Ksym].*/
			assemgstiffnessIIwithDOFelimination(gmtx, Kt, &shell, constraintmain); /*ASSEMBLAGE TANGENTIAL STIFFNESS MATRIX.*/

			gpressure = matrixvector(Tt, epressure, 6 * nnod);  /*GLOBAL EXTERNAL FORCE{Pg}.*/
			for (ii = 0; ii < nnod; ii++)
			{
				for (jj = 0; jj < 6; jj++)
				{
					*(finternal + *(loffset + (6 * ii + jj))) += *(ginternal + 6 * ii + jj);
					*(fpressure + *(loffset + (6 * ii + jj))) += *(gpressure + 6 * ii + jj);
				}
			}


			/*OUTPUT STRAIN ENERGY & STRESS*/
			if ((outputmode == 0 && (iteration == 1 || BCLFLAG == 2)) || outputmode == 1)
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

			freematrix(drccos, 3);
			freematrix(drccosinit, 3);
			freematrix(T, 6 * nnod);
			freematrix(Tt, 6 * nnod);
			freematrix(HPT, 6 * nnod);
			freematrix(Ke, 6 * nnod);
			freematrix(Kt, 6 * nnod);
			freematrix(DBe, 6 * nnod);

			free(einternal);
			free(ginternal);

			free(epressure);
			free(gpressure);

			free(eforminit);
			free(gforminit);
			free(eform);
			free(gform);
			free(edisp);

			free(loffset);

			///FOR DRAWING 3///
			if(iteration==1 && (wdraw.childs+1)->hdcC!=NULL)   /*DRAW DEFORMED ELEMENT.*/
			{
			  drawglobalshell((wdraw.childs+1)->hdcC,
							  (wdraw.childs+1)->vparam,
							  *af,shell,255,255,255,
										255,255,255,0,ONSCREEN/*,i*/);
			}
			///FOR DRAWING 3///
		}


    }
}




int eigenperiod(struct arclmframe *af, int neig, int msize, int FLAG)
{

	double **evct, *eigen; /* EIGEN VECTORS,EIGEN VALUES */
	double eps = 1.0E-8;

	evct = (double **)malloc(neig*sizeof(double *)); /* EIGEN VECTORS */
	eigen = (double *)malloc(neig*sizeof(double)); /* EIGEN VALUES */
	for (i = 0; i < neig; i++)
	{
		*(evct + i) = (double *)malloc(msize*sizeof(double));
		for (j = 0; j < msize; j++)
		{
			* (*(evct + i) + j) = 0.0;
		}
	}

	mtx1 = copygcompmatrix(gmtx, msize);
	mtx2 = copygcompmatrix(mmtx, msize);

	if(FLAG)
	{
		deigabgeneral(mtx2, mtx1, confs, msize, neig, neig, eps, eigen, evct);
	}
	else
	{
		bisecsylvester(mtx2, mtx1, confs, msize, neig, neig, eps, eigen, evct);
	}

	for (i = 0; i < neig; i++)
	{
		if (FLAG)
		{
			fprintf(fout,"DEIGABGENERAL EIGEN VALUE %ld=%.8f\n", (i + 1),*(eigen + i));

		}
		else // bisecsylvester
		{
			fprintf(fout,"BISECSYLVESTER EIGEN VALUE %ld=%.8f\n", (i + 1),1.0 / (*(eigen + i))); // bisecsylvester
		}


		if (*(eigen + i) > 0.0)
		{
			Ti = 2.0 * PI * sqrt(*(eigen + i));
			sprintf(string, "PERIOD T%ld=%.8f [sec]",(i + 1), Ti);
			fprintf(fout, "%s\n", string);
		}
		else
		{
			fprintf(fout, "ERROR:EIGEN VALUE NEGATIVE.\n");
		}


		fprintf(fout, "\nEIGEN VECTOR %ld\n",(i + 1));
		for (ii = 0; ii < nnode; ii++)
		{
			fprintf(fout,"%4ld %11.7f %11.7f %11.7f %11.7f %11.7f %11.7f\n",
			(nodes + ii)->code,
			*(*(evct + i) + 6*ii + 0),
			*(*(evct + i) + 6*ii + 1),
			*(*(evct + i) + 6*ii + 2),
			*(*(evct + i) + 6*ii + 3),
			*(*(evct + i) + 6*ii + 4),
			*(*(evct + i) + 6*ii + 5));
		}
	}
	w1 = 2.0 * PI / T1;
}
