void dbgstr(const char* str);
void dbgvct(double* vct, int size, int linesize, const char* str);
void dbgmtx(double** mtx, int rows, int cols, const char* str);
void dbggcomp(struct gcomponent* gmtx, int size, const char* str);

void inputdsp(struct arclmframe *af, const char *fname, int targetlap, int targetiteration);
void inputplst(struct arclmframe *af, const char *fname, int targetlap, int targetiteration);

double equilibriumcurvature(double* weight, double* lapddisp, double laploadfactor, double* dup, int msize);

int arclmStatic(struct arclmframe* af);
int arclmStaticExtended(struct arclmframe* af);


extern void clearwindow(struct windowparams wp);
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




void dbgmtxEigen(const Eigen::SparseMatrix<double>& mtx, const char* str)
{
    FILE* fout = fopen("hogan.dbg", "a");
    if (fout == nullptr)
        return;

    if (str != nullptr)
    {
        fprintf(fout, "%s\n", str);
    }

	int rows = mtx.rows();
    int cols = mtx.cols();


    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
		{
			fprintf(fout, "%e ", mtx.coeff(i, j));
        }
        fprintf(fout, "\n");
    }

    fclose(fout);
}

void dbggcomp(struct gcomponent* gmtx, int size, const char* str)
{
	double gdata;
	FILE *fout = fopen("hogan.dbg", "a");
	if (fout == NULL)return;

	if (str != NULL)
	{
		fprintf(fout, "%s\n", str);
	}

	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			gread(gmtx,i+1,j+1,&gdata);

			if(fabs(gdata)>1e-5)fprintf(fout, "1 ");
			else fprintf(fout, "0 ");
		}
		fprintf(fout, "\n");
	}
	fclose(fout);
}


void inputdsp(struct arclmframe *af, const char *fname, int targetlap, int targetiteration)
{
	char fullfname[256],line[1024];
	int lap, iteration, laps;
	FILE *fdsp;
	char **data;
	int n,i,j;

	snprintf(fullfname, sizeof(fullfname), "%s.%s", fname, "dsp");
	fdsp = fopen(fullfname, "r");
	if (fdsp == NULL)return;

	while (fgets(line, sizeof(line), fdsp) != NULL)
	{
		if (strstr(line, "LAP:") != NULL)
		{
			if (sscanf(line, "LAP: %d / %d ITERATION: %d", &lap, &laps, &iteration) != 3)continue;
			if (lap == targetlap && iteration == targetiteration)
			{
				for (int i = 0; i < af->nnode; i++)
				{
					data=fgetsbrk(fdsp,&n);
					for (j = 0; j < 3; j++)
					{
						af->nodes[i].d[j] = strtod(*(data + 1 + j), NULL);
						af->nodes[i].r[j] = strtod(*(data + 4 + j), NULL);
					}
					for(;n>0;n--) free(*(data+n-1));
					free(data);
				}
				break;
			}
			else
			{
				for (int i = 0; i < af->nnode; i++)
				{
					fgets(line, sizeof(line), fdsp);
				}
			}
		}
	}
	fclose(fdsp);
}

void inputplst(struct arclmframe *af, const char *fname, int targetlap, int targetiteration)
{
	char fullfname[256], line[1024];
    int lap, iteration, laps;
    FILE *fplst;
    char **data;
    int n, i, j, k;

    snprintf(fullfname, sizeof(fullfname), "%s.%s", fname, "plst");
    fplst = fopen(fullfname, "r");
    if (fplst == NULL) return;

	while (fgets(line, sizeof(line), fplst) != NULL)
	{
		if (strstr(line, "LAP:") != NULL)
		{
			if (sscanf(line, "LAP: %d / %d ITERATION: %d", &lap, &laps, &iteration) != 3)continue;
			if (lap == targetlap && iteration == targetiteration)
			{
				for (i = 0; i < af->nshell; i++)
				{
					data = fgetsbrk(fplst, &n);
					for (j = 0; j < af->shells[i].ngp; j++)
					{
						for (k = 0; k < 6; k++)
						{
							af->shells[i].gp[j].pstrain[k] = strtod(data[7 * j + k + 1], NULL);
						}
						af->shells[i].gp[j].alpha = strtod(data[7 * j + 7], NULL);
					}
					for(j = 0; j < n; j++) free(data[j]);
					free(data);
				}
				break;
            }
			else
			{
				for (i = 0; i < af->nshell; i++)
				{
					fgets(line, sizeof(line), fplst);
				}
			}
		}
    }
    fclose(fplst);
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


int arclmStatic(struct arclmframe* af)
{
	DWORDLONG memory0,memory1;
	char dir[]=DIRECTORY;
	char string[400],fname[100];
	clock_t t0;
	/*INPUT & OUTPUT FILE*/
	FILE * fin, * fonl, * fdsp, * fexf, * finf, * fubf, * frct, * fstr, * fene, * ffig, * fbcl, * feig, * fout, *fplst;         /*FILE 8 BYTES*/

	int i, j, ii, jj, n;

	int nnode, nelem, nshell, nsect, nconstraint;
	long int msize,csize;


	/*ARCLMFRAME*/
	struct osect* sects;
	struct onode* nodes;
	struct onode* ninit;
	struct owire* elems;
	struct oshell* shells;
	struct oconf* confs;
	long int* constraintmain;
	struct oconstraint* constraints;

	/*GLOBAL MATRIX*/
	struct gcomponent ginit = { 0,0,0.0,NULL };
	struct gcomponent *gmtx, * g, * p, * gcomp1;
	double gg;
	double sign, determinant;
	long int nline;
	
	/*COPY FOR TWO-POINT BUCKLING ANALYSIS*/
	struct gcomponent *gmtxcpy;

	double* givendisp;
	double* gvct;

	/*GLOBAL FORCE VECTOR*/
	double* funbalance, * fexternal, * finternal, * freaction;
	double* fbaseload, * fdeadload, * fpressure,* fgivendisp;
	double* fconstraint;
	double* constraintvct;


	///FOR ARC-LENGTH PARAMETER///
	int nlap = 1;
	int beginlap = 1;
	int targetlap = 0;
	int laps = 1000;
	int iteration = 1;
	int maxiteration = 20;
	double tolerance = 1.0e-2;

	double residual = 0.0;
	double constraintresidual = 0.0;
	double gvctlen = 0.0;

	double volume = 0.0;

	/* ENERGY */
	double Wkt = 0.0;
	double Wet = 0.0;
	double Wpt = 0.0;
	double Wot = 0.0;
	//double lastWkt = 0.0;


	/*ARC-LENGTH METHOD*/
	double arclength;
	double arcsum, predictorsign;/*FOR PREDICTOR*/
	double dupdue, dupdup;       /*FOR CORRECTOR*/
	double* due, * dup;
	double loadlambda = 0.0;

	int SCALINGARCFLAG = 0;/*ARCLENGTH CONTROL*/
	double k1, k, scaledarclength;

	int BEGININGARCFLAG = 0;/*BEGINING ARC-LENGTH*/
	double beginingarcratio = 1.0;
	int beginingarclap = 0;

	double *weight;/*WEIGHT*/
	int node;

	/*ANALYSIS MODE*/
	int neig = 0;/*TWO-POINT BUCKLING ANALYSIS*/
	int outputmode = 0;/*0:ConvergedLaps.1:AllLaps.*/
	int pinpointmode = 0;/*0:NotPinpointing.1:BisecPinpointing.2:ExtendedSystemPinpointing.*/
	int UNLOADFLAG = 0;
	int STATDYNAFLAG = 0;
	int USINGEIGENFLAG = 0;
	

	/*INITIAL CONFIG*/
	double* iform;
	double initialloadfactor = 0.0;

	/*CURRENT CONFIG*/
	double* ddisp;
	double* lambda;
	double loadfactor = 0.0;

	/*LAST CONFIG*/
	struct memoryelem* melem;
	struct memoryshell* mshell;

	double* lastddisp;
	double* lastlambda;
	double* lastpivot;
	double lastloadfactor;
	double lastsign;
	double* lapddisp;
	double laploadfactor;

	/*LAST CONFIG*/
	struct memoryelem* lastmelem;
	struct memoryshell* lastmshell;

	double* lastlastddisp;
	double* lastlastlambda;
	double* lastlastpivot;
	double lastlastloadfactor;
	double lastlastsign;
	double* lastlapddisp;
	double lastlaploadfactor;


	/*EXTENDED SYSTEM*/
	int m;/*BUCKLING DETECTED LINE*/
	double dm;/*DIAGONAL PIVOT VALUE*/
	double* evct_L, * evct_R;
	double dotLR, lenR, eigen;/*NORM OF LDL MODE*/
	int EXTENDEDFLAG = 0;


	double eps = 1e-5;
	double* epsddisp, * epsgvct, * epslambda;
	double* re, * rp;
	double evctre, evctrp;

	struct owire elem;
	struct oshell shell;

	int nnod;
	int ngp;
	int nstress;
	int* loffset;

	double* gforminit, * eforminit;
	double* gform, * eform, * edisp;

	double** drccosinit;
	double** drccos,** T,** HPT;

	double* einternal;

	double*** C, ***B;
	double** Ke, ** Kp, ** Kt;
	/*EXTENDED SYSTEM*/

	/*ANALYSIS TERMINATION*/
	int ENDFLAG = 0;
	int DIVFLAG = 0;
	int BISECTIONFLAG = 0;
	int fnode=NULL,fnodedrc=NULL;
	double fnodemin, fnodemax;

	/*FOR READING ANALYSISDATA*/
	FILE *fdata;
	int nstr, pstr, readflag;
	char** data;
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

	weight = (double *)malloc((msize + 1) * sizeof(double));
	for (i=0;i<6*nnode;i++)
	{
		*(weight + i) = 0.0;
	}
	*(weight + msize) = 1.0;

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
						if (!strcmp(*(data + pstr), "ARCLENGTH"))
						{
							pstr++;
							arclength = (double)strtod(*(data + pstr), NULL);
						}
						if (!strcmp(*(data + pstr), "SCALINGARC"))
						{
							SCALINGARCFLAG = 1;
						}
						if (!strcmp(*(data + pstr), "BEGININGARC"))
						{
							BEGININGARCFLAG=1;
							pstr++;
							beginingarcratio = (double)strtod(*(data + pstr), NULL);
							pstr++;
							beginingarclap = (int)strtol(*(data + pstr), NULL, 10);
						}
						if (!strcmp(*(data + pstr), "WEIGHT"))
						{
							pstr++;
							if (!strcmp(*(data + pstr), "LOAD"))
							{
								pstr++;
								*(weight + 6*nnode) = (double)strtod(*(data + pstr), NULL);
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
							initialloadfactor = (double)strtod(*(data + pstr), NULL);
						}
						if (!strcmp(*(data + pstr), "TARGETLAP"))
						{
							pstr++;
							targetlap = (int)strtol(*(data + pstr), NULL, 10);
						}
						if (!strcmp(*(data + pstr), "NEIG"))
						{
							pstr++;
							neig = (int)strtol(*(data + pstr), NULL, 10);
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

	sprintf(string, "FILENAME : %s\n LAPS = %d\n MAX ITERATION= %d\n ARCLENGTH = %lf\n", filename, laps, maxiteration, arclength);
	errormessage(string);

	if (outputmode == 0)errormessage("OUTPUT CONVERGED RESULT\n");
	if (outputmode == 1)errormessage("OUTPUT ALL RESULT\n");

	Eigen::setNbThreads(8);
	//n = Eigen::nbThreads( );
	//sprintf(string, "%d THREADS USING.\n", n);
	//errormessage(string);



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
	free(af->constraints);

	sects = (struct osect*)malloc(nsect * sizeof(struct osect));
	nodes = (struct onode*)malloc(nnode * sizeof(struct onode));
	ninit = (struct onode*)malloc(nnode * sizeof(struct onode));
	elems = (struct owire*)malloc(nelem * sizeof(struct owire));
	shells = (struct oshell*)malloc(nshell * sizeof(struct oshell));
	confs = (struct oconf*)malloc(msize * sizeof(struct oconf));
	constraintmain = (long int*)malloc(msize * sizeof(long int));
	constraints = (struct oconstraint*)malloc(nconstraint * sizeof(struct oconstraint));

	af->sects = sects;
	af->nodes = nodes;
	af->ninit = ninit;
	af->elems = elems;
	af->shells = shells;
	af->confs = confs;
	af->constraintmain = constraintmain;
	af->constraints = constraints;

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
#endif



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


	/*GLOBAL MATRIX FOR SOLVER*/
	gmtx = (struct gcomponent*)malloc((msize+csize) * sizeof(struct gcomponent));/*DIAGONALS OF GLOBAL MATRIX.*/
	for (i = 0; i < (msize+csize); i++)
	{
		(gmtx + i)->down = NULL;
	}

	/*GLOBAL VECTOR FOR SOLVER*/
	gvct = (double *)malloc((msize+csize) * sizeof(double));/*INCREMENTAL GLOBAL VECTOR.*/
	due = (double*)malloc((msize+csize) * sizeof(double));
	dup = (double*)malloc((msize+csize) * sizeof(double));


	/*FORCE VECTOR INITIALIZATION*/
	funbalance = (double*)malloc(msize * sizeof(double));          /*UNBALANCED INTERNAL FORCE VECTOR.*/
	freaction = (double*)malloc(msize * sizeof(double));           /*EXTERNAL FORCE VECTOR.*/
	fexternal = (double*)malloc(msize * sizeof(double));           /*EXTERNAL FORCE VECTOR.*/
	finternal = (double*)malloc(msize * sizeof(double));           /*INTERNAL FORCE VECTOR.*/
	fbaseload = (double*)malloc(msize * sizeof(double));           /*BASE LOAD VECTOR.*/
	fdeadload = (double*)malloc(msize * sizeof(double));           /*DEAD LOAD VECTOR.*/
	fpressure = (double*)malloc(msize * sizeof(double));           /*PRESSURE VECTOR.*/
	fgivendisp = (double*)malloc(msize * sizeof(double));
	fconstraint = (double*)malloc(msize * sizeof(double));


	/*POSITION VECTOR INITIALIZATION*/

	/*INITIAL CONFIGULATION*/
	iform = (double*)malloc(msize * sizeof(double));
	/*CURRENT CONFIGULATION*/
	ddisp = (double*)malloc(msize * sizeof(double));
	givendisp = (double*)malloc(msize * sizeof(double));

	/*LAST CONFIGULATION*/
	melem = (struct memoryelem*)malloc(nelem * sizeof(struct memoryelem));
	mshell = (struct memoryshell*)malloc(nshell * sizeof(struct memoryshell));
	lastddisp = (double*)malloc(msize * sizeof(double));
	lapddisp = (double*)malloc(msize * sizeof(double));		/*INCREMENT IN THE LAP*/
	lastpivot = (double*)malloc(msize * sizeof(double));    /*PIVOT SIGN OF TANGENTIAL STIFFNESS.*/
	lastlambda = (double*)malloc(csize * sizeof(double));
	/*LAST CONFIGULATION*/
	lastmelem = (struct memoryelem*)malloc(nelem * sizeof(struct memoryelem));
	lastmshell = (struct memoryshell*)malloc(nshell * sizeof(struct memoryshell));
	lastlastddisp = (double*)malloc(msize * sizeof(double));
	lastlapddisp = (double*)malloc(msize * sizeof(double));		/*INCREMENT IN THE LAP*/
	lastlastpivot = (double*)malloc(msize * sizeof(double));    /*PIVOT SIGN OF TANGENTIAL STIFFNESS.*/
	lastlastlambda = (double*)malloc(csize * sizeof(double));




	af->iform = iform;
	af->ddisp = ddisp;
	af->melem = melem;
	af->mshell = mshell;

	if(af->nlaps == NULL && inputfilename[0] != '\0' && targetlap!=NULL)
	{
	  errormessage("OUTPUT DATA READING...");
	  inputdsp(af, inputfilename, targetlap, 1);
	  inputplst(af, inputfilename, targetlap, 1);
	}

	initialform(ninit, iform, nnode);           /*ASSEMBLAGE FORMATION.*/
	initialform(nodes, ddisp, nnode);           /*ASSEMBLAGE FORMATION.*/
	initialelem(elems, melem, nelem);             /*ASSEMBLAGE ELEMENTS.*/
	initialshell(shells, mshell, nshell);         /*ASSEMBLAGE ELEMENTS.*/


	assemconf(confs,fdeadload,1.0,nnode);
	assemgivend(confs,givendisp,1.0,nnode);

	setviewpoint((wdraw.childs+0)->hwnd,*af,&((wdraw.childs+1)->vparam));
	setviewparam((wmenu.childs+2)->hwnd,(wdraw.childs+1)->vparam);

	GetAsyncKeyState(VK_LBUTTON);                   /*CLEAR KEY LEFT.*/
	GetAsyncKeyState(VK_RBUTTON);                  /*CLEAR KEY RIGHT.*/
	if(globaldrawflag==1)
	{
	  drawglobalaxis((wdraw.childs+1)->hdcC,(wdraw.childs+1)->vparam,0,0,255);                     /*DRAW GLOBAL AXIS.*/
	}


	if(af->nlaps != NULL)
	{
		nlap = af->nlaps;
	}
	beginlap = nlap;
	setincrement((wmenu.childs+2)->hwnd,laps,nlap,maxiteration,iteration);

	/*INITIAL*/
	if(nlap==1)/*INITIAL FROM ANALYSISDATA*/
	{
		loadfactor = initialloadfactor;
	}
	else if(STATDYNAFLAG == 1)/*INITIAL FROM ARCLMFRAME*/
	{

		initialloadfactor = af->loadfactor;
		loadfactor = initialloadfactor;
	}
	else
	{
		loadfactor = 0.0;
    }
	lastloadfactor = loadfactor;


	for (i = 0; i < msize; i++)
	{
		*(finternal + i) = 0.0;
		*(fexternal + i) = 0.0;
		*(funbalance + i) = 0.0;
		*(freaction + i) = 0.0;
		*(fbaseload + i) = 0.0;
		*(fpressure + i) = 0.0;
		*(fconstraint + i) = 0.0;
	}
	for(i = 0; i < csize; i++)
	{
		*(constraintvct + i) = 0.0;
	}

	volume = 0.0;
	assemshellvolume(shells, nshell, ddisp, &volume);

	/*POST PROCESS*/
	elemstress(elems, melem, nelem, constraintmain, iform, ddisp, finternal, fpressure);
	shellstress(shells, mshell, nshell, constraintmain, iform, ddisp, finternal, fpressure);
	constraintstress(constraints, nconstraint, constraintmain, iform, ddisp, lambda, fconstraint, constraintvct);

	strainenergy(af, &Wet, &Wpt);
	//kineticenergy(af, ud_m, &Wkt);

	if(UNLOADFLAG==1)
	{
		for (i = 0; i < msize; i++)
		{
			*(fdeadload + i) = -*(finternal + i);
		}
		loadfactor = -1.0;
	}

	for (i = 0; i < msize; i++)
	{
		*(fbaseload + i) = *(fdeadload + i)+*(fpressure + i);
		*(fexternal + i) = loadfactor * *(fbaseload + i);
		*(funbalance + i) = *(fexternal + i) - *(finternal + i) - *(fconstraint + i);            /*funbalance:UNBALANCED FORCE -{E}.*/
		if ((confs + i)->iconf == 1)
		{
			*(freaction + i) = *(funbalance + i);
			*(funbalance + i) = 0.0;
			*(fbaseload + i) = 0.0;
		}
		*(dup + i) = *(fbaseload + i);
		*(due + i) = *(funbalance + i);
	}
	for (i = 0; i < csize; i++)
	{
		*(dup + msize + i) = 0.0;
		*(due + msize + i) = *(constraintvct + i);
	}
	residual = vectorlength(funbalance,msize);
	constraintresidual = vectorlength(constraintvct,csize);

	/*OUTPUT(INITIAL)*/
	sprintf(string, "LAP: %5ld / %5ld ITERATION: %5ld\n", nlap, laps, iteration);

	fprintf(fdsp, string);
	outputdisp(ddisp, fdsp, nnode, nodes);

	fprintf(finf, string);
	fprintf(fexf, string);
	fprintf(fubf, string);
	fprintf(frct, string);
	outputdisp(finternal, finf, nnode, nodes);
	outputdisp(fexternal, fexf, nnode, nodes);
	outputdisp(funbalance, fubf, nnode, nodes);
	outputdisp(freaction, frct, nnode, nodes);

	/*OUTPUT FOR LOCAL*/
	fprintf(fstr, string);
	fprintf(fene, string);
	fprintf(fplst, string);
	for(i = 0; i < nshell; i++)
	{
		//fprintf(fene, "%5ld %e %e %e\n", (shells+i)->code, (mshell+i)->SEp, (mshell+i)->SEb, (mshell+i)->SE);
		fprintf(fstr, "%5ld %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n", (shells+i)->code,
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
		std::vector<Triplet> Ktriplet;
		std::vector<Triplet> Mtriplet;


		//Eigen::SparseMatrix<double,Eigen::ColMajor> Kglobal((msize+csize), (msize+csize));
		Eigen::SparseMatrix<double,Eigen::RowMajor> Kglobal((msize+csize), (msize+csize));/*FOR BiCGstab MPC*/

		/*EXECUTE BY WIN64. AMD ORDERING IS AVAILABLE ONLY BY WIN64*/
		//Eigen::SimplicialLDLT<SparseMatrix,Eigen::Lower,Eigen::NaturalOrdering<int>> solver;
		//Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
		Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;

		//Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Eigen::Upper, Eigen::DiagonalPreconditioner<double>> solver;
		//Eigen::BiCGSTAB<Eigen::SparseMatrix<double>,Eigen::IncompleteLUT<double>> solver;
		//solver.setTolerance(1e-5);
		//solver.setMaxIterations(10000);



		for (i = 1; i <= (msize+csize); i++)
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
		//comps = msize; /*INITIAL COMPONENTS = DIAGONALS.*/


		/*STIFFNESS ASSEMBLAGE*/
		if(USINGEIGENFLAG==1)
		{
		  assemelemEigen(elems, melem, nelem, constraintmain, confs, Mtriplet, Ktriplet, iform, ddisp, 1);
		  assemshellEigen(shells, mshell, nshell, constraintmain, confs, Mtriplet, Ktriplet, iform, ddisp, 1);

		  /*GLOBAL MATRIX USING EIGEN*/
		  modifytriplet(Ktriplet,confs,msize);

		  /*GLOBAL MATRIX USING EIGEN*/
		  Kglobal.reserve((msize+csize)*(msize+csize));
		  Kglobal.setFromTriplets(Ktriplet.begin(), Ktriplet.end());//SPARSE MATRIX FROM TRIPLET
		}
		else
		{
		  assemelem(elems, melem, nelem, constraintmain, NULL, gmtx, iform, ddisp);
		  assemshell(shells, mshell, nshell, constraintmain, NULL, gmtx, iform, ddisp);
		  assemconstraint(constraints, nconstraint, constraintmain, gmtx, iform, ddisp, lambda, msize);
		}
		//dbggcomp(gmtx,msize+csize,"GMTX");

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
#if 0
			solver.analyzePattern(Kglobal);
			solver.factorize(Kglobal);
#endif



			solver.compute(Kglobal);

#if 0
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
#endif
			if (solver.info() != Eigen::Success)sign=-1;
		}
		else
		{
		/*
			if(iteration==1 && neig!=0)
			{
				gmtxcpy=copygcompmatrix(gmtx,(msize+csize));
				gfree(gmtxcpy,nnode);
			}
		*/
			nline = croutluII(gmtx, confs, msize, csize, &determinant, &sign, gcomp1);
			if (sign < 0.0)
			{
				sprintf(string, "INSTABLE TERMINATION AT NODE %ld.\n",nline);
				fprintf(ffig, "%s", string);
				errormessage(string);
			}


			//for (ii = 0; ii < msize; ii++)
			//{
			//	if ((confs + ii)->iconf == 0 && ((gmtx + ii)->value) < 0 && iteration==1)
			//	{
			//		errormessage("%d\n ", ii);
			//	}
			//}

			/*
			for (ii = 0; ii < msize+csize; ii++)
			{
				if (((confs + ii)->iconf == 0 && ii<msize) || ii >= msize)
				{
					if(fabs((gmtx + ii)->value)<1e-3)
					{
						sprintf(string,"%e ", (gmtx + ii)->value);
						dbgstr(string);
					}
				}

			}
			dbgstr("\n");
			*/


			//dbggcomp(gmtx,msize+csize,"GMTX");
		}

		/*COULD NOT SOLVE*/
		if (sign < 0.0)
		{
			break;
		}

		sprintf(string, "LAP: %4d ITER: %2d {LOAD}= %5.8f {RESD}= %1.5e {DET}= %5.5f {SIGN}= %2.0f {BCL}= %1d {EPS}= %1.5e {V}= %5.5f {TIME}= %5.5f\n",
					nlap, iteration, loadfactor, residual, determinant, sign, 0, eigen, volume, 0.0);
		fprintf(ffig, "%s", string);
		errormessage(string);

#if 1
		if(iteration == 1 && (sign - lastsign) != 0 && nlap != beginlap)
		{
			/*PIVOT ZERO LINE CHECK*/
			for (ii = 0; ii < msize; ii++)
			{
				if ((confs + ii)->iconf == 0 && *(lastpivot + ii) * ((gmtx + ii)->value) < 0)
				{
					m = ii;
					sprintf(string, "BUCKLING DITECTED LAP: %4d ITER: %2d LINE: %5ld ", nlap, iteration, ii);
					fprintf(fbcl, "%s\n", string);
					errormessage(string);
				}
			}


			if(pinpointmode == 1)
			{
			  BISECTIONFLAG=1;
			}
			if(pinpointmode == 2)
			{
			  EXTENDEDFLAG=1;
			}
		}

		/*TWO-POINTS BUCKLING ANALYSIS.*/
		if(iteration == 1 && neig!=0 && (nlap-beginlap)%10==0)
		{

		}
		/*NORMAL CALCULATION IN CASE OF REGULAR*/
#endif

		if(EXTENDEDFLAG)
		{
					/*CORRECTOR CALCULATION*/
					if(USINGEIGENFLAG==1)
					{
						Vector P = Vector::Zero((msize+csize));
						for(i=0;i<msize;i++)P(i)=*(dup+i);
						Vector Up = solver.solve(P);
						for(i=0;i<(msize+csize);i++)*(dup+i)=Up(i);
						if (solver.info() != Eigen::Success)return 1;

						Vector E = Vector::Zero((msize+csize));
						for(i=0;i<msize;i++)E(i)=*(due+i);
						Vector Ue = solver.solve(E);
						for(i=0;i<(msize+csize);i++)*(due+i)=Ue(i);
						if (solver.info() != Eigen::Success)return 1;
					}
					else
					{
						nline = forwardbackwardII(gmtx, dup, confs, msize, csize, gcomp1);
						nline = forwardbackwardII(gmtx, due, confs, msize, csize, gcomp1);
					}

					evct_L = (double*)malloc(msize * sizeof(double));
					evct_R = (double*)malloc(msize * sizeof(double));


						if(1)
						{

							/*LDL^T MODE*/
							for (ii = 0; ii < msize; ii++)
							{
								*(evct_L + ii) = 0.0;
								*(evct_R + ii) = 0.0;
							}

							*(evct_L + m) = 1.0;
							*(evct_R + m) = 1.0;
							dm = (gmtx + m)->value;

							backward(gmtx, evct_L, confs, msize, csize, gcomp1);
							backward(gmtx, evct_R, confs, msize, csize, gcomp1);


							dotLR = dotproduct(evct_L, evct_R, msize);
							lenR = vectorlength(evct_R, msize);
							eigen = dm/dotLR;

							for (ii = 0; ii < msize; ii++)
							{
								*(evct_L + ii) *= lenR/dotLR;
								*(evct_R + ii) *=   1.0/lenR;
							}
							/*LDL^T MODE*/
						}


					fprintf(fbcl, "CRITICAL EIGEN VECTOR : LINE = %5ld dm = %12.9f EIGENVALUE=%12.9f\n", m, dm, eigen);

					for (ii = 0; ii < nnode; ii++)
					{
						fprintf(fbcl, "%5ld %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f\n", (nodes + ii)->code,
							*(evct_R + (6 * ii + 0)), *(evct_R + (6 * ii + 1)), *(evct_R + (6 * ii + 2)),
							*(evct_R + (6 * ii + 3)), *(evct_R + (6 * ii + 4)), *(evct_R + (6 * ii + 5)));
					}
					fprintf(fbcl, "CRITICAL EIGEN VECTOR : LINE = %5ld  dm = %12.9f\n", m, dm);
					for (ii = 0; ii < nnode; ii++)
					{
						fprintf(fbcl, "%5ld %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f\n", (nodes + ii)->code,
							*(evct_L + (6 * ii + 0)), *(evct_L + (6 * ii + 1)), *(evct_L + (6 * ii + 2)),
							*(evct_L + (6 * ii + 3)), *(evct_L + (6 * ii + 4)), *(evct_L + (6 * ii + 5)));
					}


					epsgvct = (double*)malloc(msize * sizeof(double));
					epsddisp = (double*)malloc(msize * sizeof(double));
					epslambda = (double*)malloc(csize * sizeof(double));

					for (ii = 0; ii < msize; ii++)/*UPDATE FORM FOR DIRECTIONAL DERIVATIVE*/
					{
						*(epsgvct + ii) = eps * *(evct_R + ii);
						*(epsddisp + ii) = *(ddisp + ii);
					}
					for (ii = 0; ii < msize; ii++)
					{
						*(epsgvct + ii) = *(epsgvct + *(constraintmain + ii));
					}

					/*
					for (ii = 0; ii < msize; ii++)
					{
						*(epsddisp + ii) += *(epsgvct + ii);
					}
					*/
					updaterotation(epsddisp, epsgvct, nnode);

					for (ii = 0; ii < csize; ii++)/*UPDATE FORM FOR DIRECTIONAL DERIVATIVE*/
					{
						*(epslambda + ii) = eps * *(evct_R + msize + ii);
					}

					/*ELEMENT STIFFNESS & FORCE ASSEMBLAGE*/
					elemstress(elems, melem, nelem, constraintmain, iform, epsddisp, NULL, NULL);
					shellstress(shells,mshell, nshell, constraintmain, iform, epsddisp, NULL, NULL);
					//constraintstress(constraints, nconstraint, constraintmain, iform, epsddisp, epslambda, NULL, NULL);

					re = (double*)malloc(msize * sizeof(double));
					rp = (double*)malloc(msize * sizeof(double));
					for (ii = 0; ii < msize; ii++)
					{
						*(re + ii) = 0.0;
						*(rp + ii) = 0.0;
					}

					for (i = 1; i <= nelem; i++)
					{
						inputelem(elems,melem,i-1,&elem);
						nnod=elem.nnod;
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

						edisp = extractdeformation(eforminit, eform, nnod);                	/*{Ue}*/

						T = transmatrixIII(drccos, nnod);									/*[T].*/
						HPT = transmatrixHPT(eform, edisp, T, nnod);

						einternal = matrixvector(Ke, edisp, 6 * nnod);          			/*{Fe}=[Ke]{Ue}.*/

						Kt = assemtmtxCR(Ke, eform, edisp, einternal, T, HPT, nnod);	/*TANGENTIAL MATRIX[Kt].*/
						symmetricmtx(Kt, 6*nnod);											/*SYMMETRIC TANGENTIAL MATRIX[Ksym].*/

						free(loffset);
						free(eforminit);
						free(gforminit);
						free(eform);
						free(gform);
						free(edisp);

						free(einternal);

						freematrix(drccosinit, 3);

						freematrix(drccos, 3);
						freematrix(T, 6 * nnod);
						freematrix(HPT, 6 * nnod);

						freematrix(Ke, 6 * nnod);
						//freematrix(Kp, 6 * nnod);
						freematrix(Kt, 6 * nnod);

					}

					for (i = 1; i <= nshell; i++)
					{
						inputshell(shells, mshell, i - 1, &shell);

						nnod = shell.nnod;
						ngp = shell.ngp;
						nstress = shell.nstress;

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

						B = shellB(shell);

						/*DEFORMED CONFIGFURATION*/
						for (ii = 0; ii < nnod; ii++)
						{
							inputnode(ddisp, shell.node[ii]);
						}
						drccos = shelldrccos(shell);
						gform = extractshelldisplacement(shell, ddisp);                     /*{Xg+Ug}*/
						eform = extractlocalcoord(gform,drccos,nnod); 			       	    /*{Xe+Ue}*/

						edisp = extractdeformation(eforminit, eform, nnod);           		/*{Ue}*/

						T = transmatrixIII(drccos, nnod);         							/*[T].*/
						HPT = transmatrixHPT(eform, edisp, T, nnod);

						inputshell(shells, NULL, i - 1, &shell);

						C = shellCconsistent(shell);
						einternal = assemshelleinternal(&shell, B);

						Kp = assemshellpmtx(shell,C,B);
						Kp = transformationIII(Kp, HPT, 6*nnod);/*[Ke]=[Tt][Pt][Ht][K][H][P][T]*/
						Kt = assemgmtxCR(eform, edisp, einternal, T, nnod);/*[Kg]=[Kgr]+[Kgp]+[Kgm]*/
						for (ii = 0; ii < 6*nnod; ii++)
						{
							for (jj = 0; jj < 6*nnod; jj++)
							{
								*(*(Kt + ii) + jj) += *(*(Kp + ii) + jj);/*[Kt]=[Ke]+[Kg]*/
							}
						}
						symmetricmtx(Kt, 6*nnod);											/*SYMMETRIC TANGENTIAL MATRIX[Ksym].*/

						for (ii = 0; ii < 6*nnod; ii++)
						{
							for (jj = 0; jj < 6*nnod; jj++)
							{
								*(re + *(constraintmain + *(loffset + ii))) += *(*(Kt + ii) + jj) * *(due + *(constraintmain + *(loffset + jj)));
								*(rp + *(constraintmain + *(loffset + ii))) += *(*(Kt + ii) + jj) * *(dup + *(constraintmain + *(loffset + jj)));
							}
						}

						for(ii=0;ii<ngp;ii++)freematrix(*(C+ii), nstress);
						free(C);
						for(ii=0;ii<ngp;ii++)freematrix(*(B+ii), nstress);
						free(B);

						freematrix(Kp, 6 * nnod);
						freematrix(Kt, 6 * nnod);

						free(loffset);
						free(eforminit);
						free(gforminit);
						free(eform);
						free(gform);
						free(edisp);
						free(einternal);
						freematrix(drccos, 3);
						freematrix(drccosinit, 3);
						freematrix(T, 6 * nnod);
						freematrix(HPT, 6 * nnod);
					}


					fprintf(fbcl, "RE\n", m, dm);
					for (ii = 0; ii < nnode; ii++)
					{
						fprintf(fbcl, "%5ld %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f\n", (nodes + ii)->code,
							*(re + (6 * ii + 0)), *(re + (6 * ii + 1)), *(re + (6 * ii + 2)),
							*(re + (6 * ii + 3)), *(re + (6 * ii + 4)), *(re + (6 * ii + 5)));
					}

					fprintf(fbcl, "RP\n", m, dm);
					for (ii = 0; ii < nnode; ii++)
					{
						fprintf(fbcl, "%5ld %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f\n", (nodes + ii)->code,
							*(rp + (6 * ii + 0)), *(rp + (6 * ii + 1)), *(rp + (6 * ii + 2)),
							*(rp + (6 * ii + 3)), *(rp + (6 * ii + 4)), *(rp + (6 * ii + 5)));
					}

					fprintf(fbcl, "funbalance\n", m, dm);
					for (ii = 0; ii < nnode; ii++)
					{
						fprintf(fbcl, "%5ld %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f\n", (nodes + ii)->code,
							*(funbalance + (6 * ii + 0)), *(funbalance + (6 * ii + 1)), *(funbalance + (6 * ii + 2)),
							*(funbalance + (6 * ii + 3)), *(funbalance + (6 * ii + 4)), *(funbalance + (6 * ii + 5)));
					}

					fprintf(fbcl, "fpressure\n", m, dm);
					for (ii = 0; ii < nnode; ii++)
					{
						fprintf(fbcl, "%5ld %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f\n", (nodes + ii)->code,
							*(fpressure + (6 * ii + 0)), *(fpressure + (6 * ii + 1)), *(fpressure + (6 * ii + 2)),
							*(fpressure + (6 * ii + 3)), *(fpressure + (6 * ii + 4)), *(fpressure + (6 * ii + 5)));
					}



					for (ii = 0; ii < msize; ii++)/*UPDATE FORM FOR DIRECTIONAL DERIVATIVE*/
					{
						if ((confs + ii)->iconf != 0)
						{
							*(re + ii) = 0.0;
							*(rp + ii) = 0.0;
						}
						else
						{
							*(re + ii) = -(*(re + ii)-*(funbalance + ii))/eps;
							*(rp + ii) = -(*(rp + ii)-*(fpressure + ii))/eps;
						}
					}

					fprintf(fbcl, "RE NEW\n", m, dm);
					for (ii = 0; ii < nnode; ii++)
					{
						fprintf(fbcl, "%5ld %12.9e %12.9e %12.9e %12.9e %12.9e %12.9e\n", (nodes + ii)->code,
							*(re + (6 * ii + 0)), *(re + (6 * ii + 1)), *(re + (6 * ii + 2)),
							*(re + (6 * ii + 3)), *(re + (6 * ii + 4)), *(re + (6 * ii + 5)));
					}

					fprintf(fbcl, "RP NEW\n", m, dm);
					for (ii = 0; ii < nnode; ii++)
					{
						fprintf(fbcl, "%5ld %12.9e %12.9e %12.9e %12.9e %12.9e %12.9e\n", (nodes + ii)->code,
							*(rp + (6 * ii + 0)), *(rp + (6 * ii + 1)), *(rp + (6 * ii + 2)),
							*(rp + (6 * ii + 3)), *(rp + (6 * ii + 4)), *(rp + (6 * ii + 5)));
					}





					if(0)
					{
							/*FOR NORMAL EXTENDED SYSTEM*/
							if(USINGEIGENFLAG==1)
							{
								Vector P = Vector::Zero((msize+csize));
								for(i=0;i<msize;i++)P(i)=*(rp+i);
								Vector Up = solver.solve(P);
								for(i=0;i<(msize+csize);i++)*(rp+i)=Up(i);
								if (solver.info() != Eigen::Success)return 1;

								Vector E = Vector::Zero((msize+csize));
								for(i=0;i<msize;i++)E(i)=*(re+i);
								Vector Ue = solver.solve(E);
								for(i=0;i<(msize+csize);i++)*(re+i)=Ue(i);
								if (solver.info() != Eigen::Success)return 1;
							}
							else
							{
								nline = forwardbackwardII(gmtx, rp, confs, msize, csize, gcomp1);
								nline = forwardbackwardII(gmtx, re, confs, msize, csize, gcomp1);
							}
							lenR = vectorlength(evct_R, msize);
                            loadlambda = (lenR-evctre)/evctrp;

					}

					evctre = dotproduct(re,evct_L,msize);
					evctrp = dotproduct(rp,evct_L,msize);
					fprintf(fbcl, "evctre=%12.15f evctrp=%12.15f\n", evctre, evctrp);



					loadlambda = (eigen-evctre)/evctrp;
					fprintf(fbcl, "LOADLAMBDA=%12.9f\n", loadlambda);


					for (ii = 0; ii < msize; ii++)
					{
						*(gvct + ii) = *(dup + ii) * loadlambda + *(due + ii);
					}
					/*
					for (ii = 0; ii < msize; ii++)
					{
						*(evct_R + ii) = (*(rp + ii) * loadlambda + *(re + ii));
						*(evct_L + ii) = (*(rp + ii) * loadlambda + *(re + ii));
					}
					*/
					loadfactor += loadlambda;

					free(rp);
					free(re);

					free(evct_L);
					free(evct_R);

					free(epsgvct);
					free(epsddisp);
					free(epslambda);
		}
		else
		{
			if (iteration == 1)/*PREDICTOR CALCULATION*/
			{
				if(USINGEIGENFLAG==1)
				{
				  Vector P = Vector::Zero((msize+csize));
				  for(i=0;i<(msize+csize);i++)P(i)=*(dup+i);
				  Vector Up = solver.solve(P);
				  for(i=0;i<(msize+csize);i++)*(dup+i)=Up(i);
				  if (solver.info() != Eigen::Success)return 1;
				}
				else
				{
				  nline = forwardbackwardII(gmtx, dup, confs, msize, csize, gcomp1);
				}

				/*ARCLENGTH CONTROL*/
				if (nlap < beginlap + beginingarclap)
				{
					scaledarclength = beginingarcratio * arclength;
				}
				else if (nlap == beginlap + beginingarclap || nlap == beginlap + beginingarclap + 1)
				{
					scaledarclength = arclength;
				}
				else if (nlap == beginlap + beginingarclap + 2)
				{
					if(SCALINGARCFLAG==1)
					{
						k1 = equilibriumcurvature(weight, lastlapddisp, lastlaploadfactor, dup, msize);
					}
				}
				else
				{
					if(SCALINGARCFLAG==1)
					{
						k = equilibriumcurvature(weight, lastlapddisp, lastlaploadfactor, dup, msize);
						if(k == 0.0 || k1 > k)
						{
						  scaledarclength = arclength;
						}
						else
						{
						  scaledarclength = arclength * sqrt(k1 / k);
						}
					}
				}

				/*SIGN OF PREDICTOR VECTOR & ARC-LENGTH SCALING FACTOR*/
				if(nlap == beginlap)
				{
					for (ii = 0; ii < msize; ii++)
					{
						if ((confs + ii)->iconf != 1 && *(weight + ii) != 0.0 && *(dup + ii) != 0.0)
						{
							*(weight + ii) = 1.0 / fabs(*(dup + ii));/*IF DUP OF THE TARGET DOF IS NEAR ZERO, WEIGHT IS INFINITE NUM & BAD FOR CONVERGENCE*/
						}
						else
						{
							*(weight + ii) = 0.0;
						}
					}
				}

				if (nlap == beginlap || STATDYNAFLAG == 1)
				{
					predictorsign = 1.0;
				}
				else
				{
					predictorsign = *(weight + msize) * *(weight + msize) * lastlaploadfactor * 1.0;
					for (ii = 0; ii < msize; ii++)
					{
						predictorsign += *(weight + ii) * *(weight + ii) * *(lastlapddisp + ii) * *(dup + ii);
					}

					if(predictorsign != 0.0)
					{
						predictorsign /= fabs(predictorsign);
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
				loadlambda = predictorsign * scaledarclength / sqrt(arcsum);/*INCREMANTAL LOAD FACTOR.*/
				for (ii = 0; ii < (msize+csize); ii++)
				{
					*(gvct + ii) = loadlambda * (*(dup + ii)+*(givendisp + ii));/*INCREMANTAL DISPLACEMENT.*/
				}
				fprintf(fonl, "LAP: %4d ITER: %2d {LOADLAMBDA}= % 5.8f {TOP}= % 5.8f {BOTTOM}= % 5.8f\n", nlap, iteration, loadlambda, scaledarclength, sqrt(arcsum));

				lastsign = sign;
				for (ii = 0; ii < msize; ii++)
				{
					*(lastpivot + ii) = (gmtx + ii)->value;
				}
			}
			else/*CORRECTOR CALCULATION*/
			{
				if(USINGEIGENFLAG==1)
				{
					Vector P = Vector::Zero((msize+csize));
					for(i=0;i<msize;i++)P(i)=*(dup+i);
					Vector Up = solver.solve(P);
					for(i=0;i<(msize+csize);i++)*(dup+i)=Up(i);
					if (solver.info() != Eigen::Success)return 1;

					Vector E = Vector::Zero((msize+csize));
					for(i=0;i<msize;i++)E(i)=*(due+i);
					Vector Ue = solver.solve(E);
					for(i=0;i<(msize+csize);i++)*(due+i)=Ue(i);
					if (solver.info() != Eigen::Success)return 1;
				}
				else
				{
					nline = forwardbackwardII(gmtx, dup, confs, msize, csize, gcomp1);
					nline = forwardbackwardII(gmtx, due, confs, msize, csize, gcomp1);
				}



				if (nlap == beginlap && STATDYNAFLAG == 1)
				{
					/*LOAD INCREMENT*/
					dupdue = 0.0;
					dupdup = 0.0;
					loadlambda = 0.0;
					for (ii = 0; ii < (msize+csize); ii++)
					{
						*(gvct + ii) = *(due + ii);
					}
					fprintf(fonl, "LAP: %4d ITER: %2d {LOADLAMBDA}= % 5.8f {TOP}= % 5.8f {BOTTOM}= % 5.8f\n", nlap, iteration, loadlambda, dupdue, dupdup);
				}
				else
				{
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
					loadlambda = -dupdue / dupdup;
					for (ii = 0; ii < (msize+csize); ii++)
					{
						*(gvct + ii) = loadlambda * *(dup + ii) + *(due + ii);
					}
					fprintf(fonl, "LAP: %4d ITER: %2d {LOADLAMBDA}= % 5.8f {TOP}= % 5.8f {BOTTOM}= % 5.8f\n", nlap, iteration, loadlambda, dupdue, dupdup);
					#endif


					#if 0
					/*HYPERSPHERE CONSTRAINT WITH RADIAL RETURN*/
					/*USING PREDICTOR DATA*/
					dupdue = 0.0;
					dupdup = *(weight + msize) * *(weight + msize) * loadlambda;
					for (ii = 0; ii < msize; ii++)
					{
						dupdue += *(weight + ii) * *(weight + ii) * *(gvct + ii) * *(due + ii);
						dupdup += *(weight + ii) * *(weight + ii) * *(gvct + ii) * *(dup + ii);
					}
					loadlambda += - dupdue / dupdup;
					fprintf(fonl, "LAP: %4d ITER: %2d {LOADLAMBDA}= %e {TOP}= %e {BOTTOM}= %e\n", nlap, iteration, loadlambda, dupdue, dupdup);
					for (ii = 0; ii < msize; ii++)
					{
						*(gvct + ii) += -(dupdue / dupdup) * *(dup + ii) + *(due + ii);
					}
					arcsum = *(weight + msize) * *(weight + msize) * loadlambda * loadlambda;
					for (ii = 0; ii < msize; ii++)
					{
						arcsum += *(weight + ii) * *(weight + ii) * *(gvct + ii) * *(gvct + ii);
					}
					loadlambda *= scaledarclength / sqrt(arcsum);
					loadfactor = lastloadfactor;
					for (ii = 0; ii < msize; ii++)
					{
						*(gvct + ii) *= scaledarclength / sqrt(arcsum);
						*(ddisp + ii) = *(lastddisp + ii);
					}
					#endif
				}
			}
		}




		/*UPDATE.*/
		laploadfactor += loadlambda;
		loadfactor += loadlambda;

		for (ii = 0; ii < msize; ii++)
		{
			*(gvct + ii) = *(gvct + *(constraintmain + ii));
		}
		updaterotation(lapddisp, gvct, nnode);
		updaterotation(ddisp, gvct, nnode);

		for(ii = 0; ii < nnode;ii++)
		{
			inputnode(ddisp,nodes+ii);
		}

		for (ii = 0; ii < csize; ii++)
		{
		  *(lambda+ii)+=*(gvct+msize+ii);
		}


		for (i = 0; i < msize; i++)	 /*GLOBAL VECTOR INITIALIZATION.*/
		{
			*(finternal  + i) = 0.0;
			*(fexternal  + i) = 0.0;
			*(funbalance + i) = 0.0;
			*(freaction  + i) = 0.0;
			*(fbaseload  + i) = 0.0;
			*(fpressure  + i) = 0.0;
			*(fconstraint + i) = 0.0;
		}
		for (i = 0; i < csize; i++)	 /*GLOBAL VECTOR INITIALIZATION.*/
		{
			*(constraintvct  + i) = 0.0;
		}

		volume = 0.0;
		assemshellvolume(shells, nshell, ddisp, &volume);

		/*POST PROCESS*/
		elemstress(elems, melem, nelem, constraintmain, iform, ddisp, finternal, fpressure);
		shellstress(shells, mshell, nshell, constraintmain, iform, ddisp, finternal, fpressure);
		constraintstress(constraints, nconstraint, constraintmain, iform, ddisp, lambda, fconstraint, constraintvct);

		strainenergy(af, &Wet, &Wpt);
		//kineticenergy(af, ud_m, &Wkt);


		for (i = 0; i < msize; i++)
		{
			*(fbaseload + i) = *(fdeadload + i)+*(fpressure + i);
			*(fexternal + i) = loadfactor * *(fbaseload + i);
			*(funbalance + i) = *(fexternal + i) - *(finternal + i) - *(fconstraint + i);            /*funbalance:UNBALANCED FORCE -{E}.*/
			if ((confs + i)->iconf == 1)
			{
				*(freaction + i) = *(funbalance + i);
				*(funbalance + i) = 0.0;
				*(fbaseload + i) = 0.0;
			}
			*(dup + i) = *(fbaseload + i);
			*(due + i) = *(funbalance + i);
		}
		for (i = 0; i < csize; i++)	 /*GLOBAL VECTOR INITIALIZATION.*/
		{
			*(dup + msize + i) = 0.0;
			*(due + msize + i) = *(constraintvct  + i);
		}

		/*CONVERGENCE JUDGEMENT.*/
		residual = vectorlength(funbalance,msize);
		constraintresidual = vectorlength(constraintvct,csize);
		gvctlen = vectorlength(gvct,msize+csize);

		if (!isfinite(sign) || sign > 10 || !isfinite(residual) || residual > 1e+10)
		{
			DIVFLAG = 1;
			ENDFLAG = 1;

			nlap++;
			iteration = 0;
			sprintf(string,"DIVERGENCE DITECTED(SIGN = %f RESIDUAL = %f). ANALYSIS TERMINATED.\n", sign, residual);
			errormessage(string);
		}
		else
		{
			if(EXTENDEDFLAG == 1)
			{
				if((residual < tolerance && fabs(eigen) < 1e-8)  && iteration != 1)
				{
					EXTENDEDFLAG = 0;
					eigen=0.0;


					nlap++;
					iteration = 0;
				}
				clearwindow(*(wdraw.childs+1));
				drawarclmframe((wdraw.childs+1)->hdcC,(wdraw.childs+1)->vparam,*af,0,ONSCREEN);
				overlayhdc(*(wdraw.childs + 1), SRCPAINT);                  /*UPDATE DISPLAY.*/
			}
			else
			{
				DIVFLAG = 0;
				ENDFLAG = 0;
				if ((residual < tolerance || iteration > maxiteration-1) && iteration != 1)
				{
					nlap++;
					iteration = 0;

					if(residual>1e+2)
					{
					  DIVFLAG=1;
					  ENDFLAG=1;
					}

					clearwindow(*(wdraw.childs+1));
					drawarclmframe((wdraw.childs+1)->hdcC,(wdraw.childs+1)->vparam,*af,0,ONSCREEN);
					overlayhdc(*(wdraw.childs + 1), SRCPAINT);                  /*UPDATE DISPLAY.*/
				}
			}
		}
		iteration++;
		setincrement((wmenu.childs+2)->hwnd,laps,nlap,maxiteration,iteration);

		/*OUTPUT.*/
		if (((outputmode == 0 && iteration == 1) || outputmode == 1) && DIVFLAG == 0)
		{
			sprintf(string, "LAP: %5ld / %5ld ITERATION: %5ld\n", nlap, laps, iteration);

			/*OUTPUT FOR GLOBAL*/
			//dbgstr(string);

			fprintf(fdsp, string);
			outputdisp(ddisp, fdsp, nnode, nodes);

			fprintf(finf, string);
			fprintf(fexf, string);
			fprintf(fubf, string);
			fprintf(frct, string);
			outputdisp(finternal, finf, nnode, nodes);
			outputdisp(fexternal, fexf, nnode, nodes);
			outputdisp(funbalance, fubf, nnode, nodes);
			outputdisp(freaction, frct, nnode, nodes);

			/*OUTPUT FOR LOCAL*/
			fprintf(fstr, string);
			fprintf(fene, string);
			fprintf(fplst, string);
			for(i = 0; i < nshell; i++)
			{
				fprintf(fstr, "%5ld %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n", (shells+i)->code,
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
#if 0
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
			if (nlap>20 && loadfactor < 0.0)
			{
				ENDFLAG = 1;
				sprintf(string,"NEGATIVE LOAD DITECTED(LOADFACTOR = %8.5f). ANALYSIS TERMINATED.\n", loadfactor);
				errormessage(string);
			}
		}
#endif

		if (iteration == 1)
		{
			if (DIVFLAG == 1)/*TERMINATION AT LAST LAP*/
			{
				laploadfactor = 0.0;
				loadfactor = lastloadfactor;
				for (ii = 0; ii < msize; ii++)
				{
					*(lapddisp + ii) = 0.0;
					*(ddisp + ii) = *(lastddisp + ii);
				}
				for(ii = 0; ii < nnode;ii++)
				{
					inputnode(ddisp,nodes+ii);
				}
				for (ii = 0; ii < csize; ii++)
				{
					*(lambda+ii) = *(lastlambda+ii) ;
				}

			}
			else/*LAST LAP DATA UPDATE.*/
			{
				/*copy data from last-lap to last-last-lap*/
				memcpy(lastmelem, melem, nelem * sizeof(struct memoryelem));
				memcpy(lastmshell, mshell, nshell * sizeof(struct memoryshell));
				lastlaploadfactor = laploadfactor;
				lastlastloadfactor = lastloadfactor;
				lastlastsign = lastsign;
				for (ii = 0; ii < msize; ii++)
				{
					*(lastlapddisp  + ii) = *(lapddisp  + ii);
					*(lastlastddisp + ii) = *(lastddisp + ii);
					*(lastlastpivot + ii) = *(lastpivot + ii);
				}
				for (ii = 0; ii < csize; ii++)
				{
					*(lastlastlambda + ii) = *(lastlambda + ii);
				}

				/*copy data from current-lap to last-lap*/
				initialelem(elems,melem,nelem);    /*melem (LAST LAP DATA) UPDATE.*/
				initialshell(shells,mshell,nshell);/*mshell(LAST LAP DATA) UPDATE.*/
				laploadfactor = 0.0;
				lastloadfactor = loadfactor;
				//lastsign = sign;
				for (ii = 0; ii < msize; ii++)
				{
					*(lapddisp  + ii) = 0.0;
					*(lastddisp + ii) = *(ddisp + ii);
					//*(lastpivot + ii) = (gmtx + ii)->value;
				}
				for (ii = 0; ii < csize; ii++)
				{
					*(lastlambda + ii) = *(lambda + ii);
				}
			}


			for(i = 0; i < msize; i++)
			{
				//*(finertial  + i) = 0.0;
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

			/*POST PROCESS*/
			elemstress(elems, melem, nelem, constraintmain, iform, ddisp, finternal, fpressure);
			shellstress(shells, mshell, nshell, constraintmain, iform, ddisp, finternal, fpressure);
			constraintstress(constraints, nconstraint, constraintmain, iform, ddisp, lambda, fconstraint, constraintvct);

			//strainenergy(af, &Wet, &Wpt);
			//kineticenergy(af, ud_m, &Wkt);


			for (i = 0; i < msize; i++)
			{
				*(fbaseload + i) = *(fdeadload + i)+*(fpressure + i);
				*(fexternal + i) = loadfactor * *(fbaseload + i);
				*(funbalance + i) = *(fexternal + i) - *(finternal + i) - *(fconstraint + i);            /*funbalance:UNBALANCED FORCE -{E}.*/
				if ((confs + i)->iconf == 1)
				{
					*(freaction + i) = *(funbalance + i);
					*(funbalance + i) = 0.0;
					*(fbaseload + i) = 0.0;
				}
				*(dup + i) = *(fbaseload + i);
				*(due + i) = *(funbalance + i);
			}

			for (i = 0; i < csize; i++)	 /*GLOBAL VECTOR INITIALIZATION.*/
			{
				*(dup + msize + i) = 0.0;
				*(due + msize + i) = *(constraintvct  + i);
			}

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
			fclose(fout);
			fclose(feig);
            fclose(fplst);

			gfree(gmtx,nnode);

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
			free(fconstraint);

			free(due);
			free(dup);

			free(lastddisp);
			free(lastpivot);
			free(lapddisp);


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

	////////// OUTPUT RESULT FOR SRCAN//////////
	/*
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
	*/


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
	fclose(fplst);

	gfree(gmtx,nnode);

	free(weight);
	free(gvct);

	free(due);
	free(dup);

	free(funbalance);
	free(freaction);
	free(fexternal);
	free(finternal);
	free(fpressure);
	free(fbaseload);
	free(fdeadload);
	free(fgivendisp);
	free(fconstraint);
	free(constraintvct);

	free(givendisp);

	free(lastddisp);
	free(lastpivot);
	free(lapddisp);

	free(iform);
	free(ddisp);
	free(melem);
	free(mshell);

	errormessage(" ");
	errormessage("COMPLETED.");

	memory1=availablephysicalmemoryEx("REMAIN:");
	sprintf(string,"CONSUMPTION:%ld[BYTES]",(memory0-memory1));
	errormessage(string);

	return 0;
}
/*arclmCR*/

#if 0
int twopointsBuckling(struct arclmframe* af, struct gcomponent* gmtx, double* ddisp, double* lambda, double* gvct, int neig, int msize, int csize)
{
	int i,j,ii,jj;
	char string[400]
	int nnode, nelem, nshell, nsect, nconstraint;
	struct gcomponent *epsgmtx, *dgmtx, *p,*g;;
	struct gcomponent ginit = { 0,0,0.0,NULL };

	double **evct, *eigen;
	double* epsddisp, * epsgvct, * epslambda;
	double eps=1e-3;

	nnode=af->nnode;
	nelem=af->nelem;
	nshell=af->nshell;
	nsect=af->nsect;
	nconstraint=af->nconstraint;

	/*BUCKLING ANALYSIS INITIALIZATION.*/
	epsgmtx = (struct gcomponent*)malloc(msize * sizeof(struct gcomponent));/*DIAGONALS OF GLOBAL MATRIX.*/
	for (i = 0; i < msize; i++)
	{
		(epsgmtx + i)->down = NULL;
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

	epsddisp = (double*)malloc(msize * sizeof(double));
	epsgvct = (double*)malloc(msize * sizeof(double));
	epslambda = (double*)malloc(csize * sizeof(double));
	for (ii = 0; ii < msize; ii++)/*UPDATE FORM FOR DIRECTIONAL DERIVATIVE*/
	{
		*(epsgvct + ii) = eps * (*(gvct + ii));
		*(epsddisp + ii) = *(ddisp + ii);
	}
	for (ii = 0; ii < msize; ii++)
	{
		*(epsgvct + ii) = *(epsgvct + *(constraintmain + ii));
	}
	updateform(epsddisp, epsgvct, nnode);

	for (ii = 0; ii < csize; ii++)/*UPDATE FORM FOR DIRECTIONAL DERIVATIVE*/
	{
		*(epslambda + ii) = eps * (*(gvct + msize + ii));
	}

	/*ELEMENT STIFFNESS & FORCE ASSEMBLAGE*/

	elemstress( af->elems, af->melem, af->nelem, af->constraintmain, af->iform, epsddisp, NULL, NULL);
	shellstress( af->shells, af->mshell, af->nshell, af->constraintmain, af->iform, epsddisp, NULL, NULL);
	//constraintstress( af->constraints,  af->nconstraint,  af->constraintmain,  af->iform, epsddisp, epslambda, NULL, NULL);

	assemelem(af->elems, af->melem, af->nelem, af->constraintmain, NULL, epsgmtx, af->iform, epsddisp);
	assemshell(af->shells, af->mshell, af->nshell, af->constraintmain, NULL, epsgmtx, af->iform, epsddisp);
	//assemconstraint(af->constraints, af->nconstraint, af->constraintmain, epsgmtx, af->iform, epsddisp, epslambda, msize);

	dgmtx=gcomponentadd3(epsgmtx,1.0/eps,gmtxcpy,-1.0/eps,msize);
	bisecgeneral(dgmtx,-1.0,gmtxcpy,-1.0,confs,msize,neig,sign,1.0E-10,eigen,evct,0.0,1.0E+3);

	/*OUTPUT*/
	if(feig!=NULL)
	{
		for(i=0;i<neig;i++)
		{
			fprintf(feig,"LAP = %d MODE = %d GENERALIZED EIGENVALUE = %e STANDARD EIGENVALUE = %e ", nlap, (i+1), *(eigen + i), 0.0);

			if (*(eigen + i) > 0.0)
			{
				*(eigen+i) = 1.0/(*(eigen + i));
				fprintf(feig, "LOADLAMBDA = %e BUCKLINGLOAD = %e\n",*(eigen+i),loadfactor+*(eigen+i));
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
				(nodes + ii)->code,
				*(*(evct + i) + 6*ii + 0),
				*(*(evct + i) + 6*ii + 1),
				*(*(evct + i) + 6*ii + 2),
				*(*(evct + i) + 6*ii + 3),
				*(*(evct + i) + 6*ii + 4),
				*(*(evct + i) + 6*ii + 5));
			}
		}
	}

	gfree(epsgmtx,nnode);
	gfree(dgmtx,nnode);

	free(eigen);
	freematrix(evct,neig);

    free(epsddisp);
	free(epsgvct);
	free(epslambda);

	return 1;
}

#endif

int bisecPinpoint(struct arclmframe* af, struct memoryelem* melem, struct memoryshell* mshell,
				  double loadfactor_R, double loadfactor_L,
				  double* ddisp_R, double* ddisp_L,
				  double* lambda_R, double* lambda_L,
				  double *evctinit,
				  double *fexternal,
				  int msize, int csize,
				  int nlap, double lastsign,
				  FILE* fbcl, FILE* feig)
{
	int i,j,ii,jj;
	char string[400];
	int nnode, nelem, nshell, nsect, nconstraint;
	/*FOR BISECSYLVESTER*/
	double LL,LR,LM;

	double loadfactor;
	double* ddisp;

	double eigen;
	double* evct;
	
	/*GLOBAL MATRIX*/
	struct gcomponent ginit = { 0,0,0.0,NULL };
	struct gcomponent *gmtx, * g, * p, * gcomp1;
	double gg;
	double sign, determinant;
	long int nline;
	
	double eps = 1.0e-5;/*FOR CONVERGENT JUDGEMENT*/
	double evctdot;/*FOR BIFURCATION JUDGEMENT*/

	int iteration = 1;

	nnode=af->nnode;
	nelem=af->nelem;
	nshell=af->nshell;
	nsect=af->nsect;
	nconstraint=af->nconstraint;

	ddisp = (double*)malloc(msize * sizeof(double));
	evct = (double*)malloc(msize * sizeof(double));
	gmtx = (struct gcomponent*)malloc(msize * sizeof(struct gcomponent));/*DIAGONALS OF GLOBAL MATRIX.*/
	for (i = 0; i < msize; i++)
	{
		(gmtx + i)->down = NULL;
	}


	LR = 1.0;
	LL = 0.0;
	while(LR - LL > eps)
	{
		LM = 0.5 * (LL + LR);
		loadfactor = (1 - LM) * loadfactor_L + LM * loadfactor_R;
		for (ii = 0; ii < msize; ii++)
		{
			*(ddisp + ii) = (1 - LM) * *(ddisp_L + ii) + LM * *(ddisp_R + ii);
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


		elemstress(af->elems, melem, af->nelem,  af->constraintmain,  af->iform, ddisp, NULL, NULL);
		shellstress(af->shells, mshell, af->nshell,  af->constraintmain,  af->iform, ddisp, NULL, NULL);
		//constraintstress( af->constraints,  af->nconstraint,  af->constraintmain,  af->iform, ddisp, lambda, NULL,NULL);

		assemelem(af->elems, melem, af->nelem,  af->constraintmain, NULL, gmtx,  af->iform, ddisp);
		assemshell(af->shells, mshell, af->nshell,  af->constraintmain, NULL, gmtx,  af->iform, ddisp);
		//assemconstraint( af->constraints,  af->nconstraint,  af->constraintmain, gmtx,  af->iform, ddisp, lambda, msize);

		/*SOLVE*/
		nline = croutluII(gmtx, af->confs, msize, csize, &determinant, &sign, gcomp1);/*FOT COUNTING NEGATIVE PIVOT*/
		sprintf(string, "LAP: %4d ITER: %2d {LOAD}= %5.8f {RESD}= %1.5e {DET}= %5.5f {SIGN}= %2.0f {BCL}= %1d {EPS}= %1.5e {V}= %5.5f {TIME}= %5.5f\n",
			nlap, iteration, loadfactor, 0.0, determinant, sign, 0, 0.0, 0.0, 0.0);
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
					sprintf(string, "INSTABLE TERMINATION AT NODE %ld.",(af->nodes + int((ii - 1) / 6))->code);
					fprintf(fbcl, "%s\n", string);
					errormessage(string);
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
		
		iteration++;
	}

	/*INITIAL EIGEN VECTOR*/
	for (ii = 0; ii < msize; ii++)
	{
		*(evct + ii) = *(evctinit + ii);
	}
	/*INVERS METHOD FOR EIGEN VECTOR*/
	eigen=inversemethod(gmtx,af->confs,evct,msize);


	/*BIFURCATION JUDGEMANT*/
	if (fexternal!=NULL)
	{
		/*BIFURCATION JUDGE*/
		vectornormalize(fexternal, msize);
		evctdot = dotproduct(evct, fexternal, msize);
	}


	/*OUTPUT*/
	if(feig!=NULL)
	{
		fprintf(feig,"LAP = %d MODE = %d GENERALIZED EIGENVALUE = %e STANDARD EIGENVALUE = %e LOADLAMBDA = %e DOT = %e\n",
		nlap, 1, LM, eigen, loadfactor-loadfactor_L, evctdot);

		for (ii = 0; ii < msize; ii++)
		{
			*(evct + ii) = *(evct + *(af->constraintmain + ii));
		}
		for (ii = 0; ii < nnode; ii++)
		{
			fprintf(feig,
			"%4ld %e %e %e %e %e %e\n",
			(af->nodes + ii)->code,
			*(evct + 6*ii + 0),
			*(evct + 6*ii + 1),
			*(evct + 6*ii + 2),
			*(evct + 6*ii + 3),
			*(evct + 6*ii + 4),
			*(evct + 6*ii + 5));
		}
	}
	free(ddisp);
	free(evct);
	gfree(gmtx,nnode);
	return 1;

}

#if 0

/*FOR PATH SWITCH*/
	double gradeps = 1.0e-8;
	double gamma;
	double evctfunbalance, epsevctfunbalance;
	double* fswitching;
	double * epsfunbalance, * epsfinternal, * epsfexternal, * epsfpressure;
	epsfunbalance = (double*)malloc(6 * nnode * sizeof(double));
	epsfinternal = (double*)malloc(6 * nnode * sizeof(double));
	epsfexternal = (double*)malloc(6 * nnode * sizeof(double));
	epsfpressure = (double*)malloc(6 * nnode * sizeof(double));


	fswitching = (double*)malloc(msize * sizeof(double));
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
	nline = croutluII(gmtx, confs, msize, csize, &determinant, &sign, gcomp1);

	sprintf(string, "LAP: %4d ITER: %2d {LOAD}= % 5.8f {RESD}= %1.6e {DET}= %8.5f {SIGN}= %2.0f {BCL}= %1d {EPS}=%1.5e {V}= %8.5f\n",
		nlap, iteration, loadfactor, residual, determinant, sign, BCLFLAG, gamma, volume);
	fprintf(ffig, "%s", string);
	errormessage(string);


	if (iteration == 1)/*PREDICTOR CALCULATION*/
	{
		nline = forwardbackwardII(gmtx, dup, confs, msize, csize, gcomp1);
		scaledarclength = 1e-3 * arclength;


		/*INCREMENTAL CALCULATION*/
		arcsum = *(weight + msize) * *(weight + msize);
		for (ii = 0; ii < msize; ii++)
		{
			arcsum += *(weight + ii) * *(weight + ii) * *(dup + ii) * *(dup + ii);/*SCALED DISPLACEMENT.*/
		}
		loadlambda = scaledarclength / sqrt(arcsum);/*INCREMANTAL LOAD FACTOR.*/
		for (ii = 0; ii < msize; ii++)
		{
			if ((confs + ii)->iconf != 1)
			{
				*(gvct + ii) = gamma * *(dup + ii);/*INCREMANTAL DISPLACEMENT.*/
			}
		}
		fprintf(fonl, "LAP: %4d ITER: %2d {LOADLAMBDA}= % 5.8f {TOP}= % 5.8f {BOTTOM}= % 5.8f\n", nlap, iteration, loadlambda, scaledarclength, sqrt(arcsum));
	}
	else/*CORRECTOR CALCULATION*/
	{
		nline = forwardbackwardII(gmtx, dup, confs, msize, csize, gcomp1);
		nline = forwardbackwardII(gmtx, due, confs, msize, csize, gcomp1);
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
		loadlambda = -dupdue / dupdup;
		for (ii = 0; ii < msize; ii++)
		{
			if ((confs + ii)->iconf != 1)
			{
				*(gvct + ii) = gamma * *(dup + ii) + *(due + ii);
			}
		}
		fprintf(fonl, "LAP: %4d ITER: %2d {LOADLAMBDA}= % 5.8f {TOP}= % 5.8f {BOTTOM}= % 5.8f\n", nlap, iteration, loadlambda, dupdue, dupdup);
	}

	gamma += loadlambda;
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


	shellstress(shells, mshell, nshell, constraintmain,iform, epsddisp, epsfinternal, epsfexternal);

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



#endif


