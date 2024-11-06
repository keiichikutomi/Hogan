int arclmStaticEx(struct arclmframe* af);

void updaterotationEx(double* ddisp, double* gvct, int nnode);
double equilibriumcurvature(double* weight, double* lapddisp, double laploadfactor, double* dup, int msize);

/*
PINPOINTMODE
0:DO NOT DETECT BUCKLING
1:CALCULATE EIGEN VECTOR
*/



int arclmStaticEx(struct arclmframe* af)
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
	double sign, determinant;
	long int nline;

	/*GLOBAL VECTOR*/
	double* ddisp, * iform;

	double* gvct;
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

	/***GLOBAL VECTOR***/
	gvct = (double *)malloc(msize * sizeof(double));/*INCREMENTAL GLOBAL VECTOR.*/



	lastddisp = (double*)malloc(msize * sizeof(double));
	lastgvct = (double*)malloc(msize * sizeof(double));
	lastpivot = (double*)malloc(msize * sizeof(double));    /*PIVOT SIGN OF TANGENTIAL STIFFNESS.*/
	lapddisp = (double*)malloc(msize * sizeof(double));		/*INCREMENT IN THE LAP*/
	nextddisp = (double*)malloc(msize * sizeof(double));

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

	Vector fdeadload  = Vector::Zero(msize);
	Vector givendisp  = Vector::Zero(msize);

	assemconf(confs,fdeadload,1.0,nnode);               /*GLOBAL VECTOR.*/
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



		SparseMatrix gmtx(msize, msize);
		/*GLOBAL MATRIX USING EIGEN*/
		std::vector<Triplet> Ktriplets;

		/*GLOBAL VECTOR INITIALIZATION.*/
		/*
		std::vector<double> finternal(msize,0.0);
		*/
		Vector finternal  = Vector::Zero(msize);
		Vector fexternal  = Vector::Zero(msize);
		Vector funbalance = Vector::Zero(msize);
		Vector freaction  = Vector::Zero(msize);

		Vector fbaseload  = Vector::Zero(msize);
		Vector fpressure  = Vector::Zero(msize);

		Vector fgivendisp = Vector::Zero(msize);


		assemshellvolume(shells, nshell, ddisp, &volume);

		/*ELEMENT STIFFNESS & FORCE ASSEMBLAGE*/
		assemelemEx (elems,  melem,  nelem,  constraintmain, NULL, gmtx, iform, ddisp, finternal, fpressure);
		assemshellEx(shells, mshell, nshell, constraintmain, NULL, gmtx, iform, ddisp, finternal, fpressure);

		if(iteration==1)
		{
		  clearwindow(*(wdraw.childs+1));
		  drawarclmframe((wdraw.childs+1)->hdcC,(wdraw.childs+1)->vparam,*af,0,ONSCREEN);
		  overlayhdc(*(wdraw.childs + 1), SRCPAINT);                  /*UPDATE DISPLAY.*/
		}


		/*EXTERNAL FORCE & UNBALANCED FORCE*/
		/*
		if(iteration==1)
		{
			for (i = 0; i < msize; i++)
			{
				*(fgivendisp + i) = *(givendisp + i);
			}
			modifygivend(gmtx,fgivendisp,confs,nnode);
		}
		*/

		if(iteration==1)
		{
			for(i = 0; i < nshell; i++)
			{
				outputmemoryshell(shells,mshell,i);/*Initialshell‚Å‚Í*/
            }
		}
		residual = 0.0;
		/*
		if(UNLOADFLAG==1)
		{
			for (i = 0; i < msize; i++)
			{
				if(nlap==1 && iteration==1)
				{
					*(fbaseload + i) = -*(finternal + i);
				}
				*(dup + i) = *(fbaseload + i)+*(fpressure + i);
				*(fexternal + i) = - *(fbaseload + i) + loadfactor * *(dup + i);
				*(funbalance + i) = *(fexternal + i) - *(finternal + i);
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
		*/

		for (i = 0; i < msize; i++)
		{
			fbaseload = fdeadload + fpressure;
			fexternal = loadfactor * P;


			funbalance = fexternal - finternal;
			if ((confs + i)->iconf == 1)
			{
				freaction[i] = funbalance[i];
				funbalance[i] = 0.0;
			}
			residual += *(funbalance + i) * *(funbalance + i);

			if(iteration==1)fbaseload += fgivendisp;
		}

		/*OUTPUT*/
		if ((outputmode == 0 && iteration == 1) || outputmode == 1)
		{
			outputdispEx(ddisp,      fdsp, nnode, nodes);
			outputdispEx(finternal,  finf, nnode, nodes);
			outputdispEx(fexternal,  fexf, nnode, nodes);
			outputdispEx(funbalance, fubf, nnode, nodes);
			outputdispEx(freaction,  frct, nnode, nodes);

			/*STRESS & ENERGY OUTPUT*/
			for(i = 0; i < nshell; i++)
			{
				fprintf(fstr, "%5ld %e %e %e %e %e %e\n",
				 (shells+i)->code,
				((shells+i)->gp[0]).stress[0],
				((shells+i)->gp[0]).stress[1],
				((shells+i)->gp[0]).stress[2],
				((shells+i)->gp[0]).stress[3],
				((shells+i)->gp[0]).stress[4],
				((shells+i)->gp[0]).stress[5]);
			}
		}




		if (BCLFLAG < 1)/*REGULAR*/
		{
			gmtx.setFromTriplets(triplets.begin(), triplets.end());

			boundary(gmtx,confs);


			Eigen::SimplicialLDLT<SparseMatrix> solver;
			solver.compute(gmtx);

			determinant = solver.vectorD().array().log().sum();
			sign = 0.0;
			for (i = 0; i < msize; i++)
			{
			  if(solver.vectorD()<0.0) sign += 1;
			}



			if (sign < 0.0 || solver.info() != Eigen::Success)
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
				return -1;
			}
			/*
			Eigen::SparseLU<SparseMatrix> solver;
			solver.analyzePattern(K_global);
			solver.factorize(K_global);
			*/
			sprintf(string, "LAP: %4d ITER: %2d {LOAD}= % 5.8f {RESD}= %1.6e {DET}= %8.5f {SIGN}= %2.0f {BCL}= %1d {EPS}= %1.5e {V}= %8.5f\n",
					nlap, iteration, loadfactor, residual, determinant, sign, BCLFLAG, 0.0, volume);
			fprintf(ffig, "%s", string);
			errormessage(string);







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
					//nline = forwardbackward(gmtx, dup, confs, msize, gcomp1);



					Vector dup = solver.solve(P);
					if (solver.info() != Eigen::Success)return -1;


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
				else/*CORRECTOR CALCULATION*/
				{
					Vector dup = solver.solve(fbaseload);
					Vector due = solver.solve(funbalance);
					if (solver.info() != Eigen::Success)return -1;

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

			iteration++;
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


			free(weight);
			free(gvct);


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
		}
#endif
	}





	assemelemEx(elems, melem, nelem, constraintmain, NULL, NULL, iform, ddisp, NULL, NULL);
	assemshellEx(shells, mshell, nshell, constraintmain, NULL, NULL, iform, ddisp, NULL, NULL);

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


	free(weight);
	free(gvct);

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


