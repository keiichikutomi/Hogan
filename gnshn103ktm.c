                                                                           /* ========================================================= */
/* PROGRAM GUNBASHIRA SHINDOU 002 */
/* 6 D.O.F. 3D FRAME WITH POINT MASS */
/* GRAVITY LOADED. */
/* HISENKEI OUTOU KAISEKI */
/* CODED BY JUN SATO */
/* DATE:FORTRAN SINSE 1993. 9.17 */
/* C       SINSE 1998.10.10 */
/* */
/* BASSUI FROM: */
/* 'ZAKU001'  ,BASED ON 'Y6FRAMEB' BY YOSHINOBU FUJITANI */
/* 'D24PE3..9',BY TSUNETA OCHI */
/* ========================================================= */

#define MAXLAP 1000
/* #define MAXLAP 6005 */ /* 6000x0.010 FOR HAKONE,KOKUJI HACHINOHE UD */
/* #define MAXLAP 2005 */ /* 2000x0.010 FOR KOKUJI UD */
/* #define MAXLAP 1205 */ /* 1200x0.010 FOR KANSOKU UD */
/* #define MAXLAP 116 */  /* 116x0.010 FOR EL CENTRO UD */
/* #define MAXLAP 2005 *//* FOR BUCHAREST */
/* #define MAXLAP 4005 */ /* FOR E-DEFENSE */
/* #define MAXLAP 1 */
#define DAMPINGTYPE 1 /*1:PROPORTION TO INITIAL 2:TO OCCASIONAL 3:TO OCCASIONAL+AXIS=0*/
/* #define DAMPINGTYPE 3 */ /* FOR HAKUSHIMA */
#define FLOORS 2 /*INCLUDING GL,RFL*/
/* #define FLOORS 9 */  /* c200-c800,e100 */
/* #define FLOORS 13 */ /* c900,d100,d200,d300 */
/* #define FLOORS 10 */ /* INCLUDING GL,RFL */
#define DAMPING 0.020                                /*DAMPING RATE*/
/* #define DAMPING 0.012 */              /* DAMPING RATE FOR E-DEFENSE */
/* #define DAMPING 0.010 */               /* DAMPING RATE FOR KAMISATO */
/* #define DAMPING 0.000 */                            /* DAMPING RATE */

#define BCLNGCONDENSATION 1 //ujioka
/*
  Yield Surface
  0:Normal , 1:Revised(Buckling Condensation)
*/


struct accdata
{
	long int ndata;
	double dt;
	double *a;
}; /* ACCELERATION DATA. */

struct structrigid
{
	int ielem, floor;
	double x, y;
	double K[2];
}; /* SIZE 50 BYTES */

int gnshn101(struct arclmframe *af);
int gnshn102(struct arclmframe *af);
double *inputacceleration(FILE *fin, long int *ndata, double *dt);
double assemmass(struct gcomponent *mmtx, long int nnode, struct oconf *confs,
	double *nmass);
double accelerationincrement(struct accdata acc, double ddt, long int lap);
void assemaccel(double *gvct, double dacc[], long int nnode,
	struct arclmframe *af);
double newmarkbeta(struct gcomponent *gmtx, struct gcomponent *cmtx,
	struct gcomponent *mmtx, double *gacc, struct oconf *confs, int nnode,
	double *u, double *ud, double *udd, double *du, double *dud, double *dudd,
	double ddt, double beta);
void floorvalues(FILE *fout, FILE *ftxt, struct arclmframe *af, double *ddisp,
	double dacc[]);
double *stresstransform(struct owire *elem, double *estress);

/* EXTERNAL PARAMETERS */
extern FILE *globalfile; /* GLOBAL FILE. */
extern void clearwindow(struct windowparams wp);
extern void drawtexts(HDC hdc, struct snode *strset, int nstr,
	struct viewparam vp);

int gnshn101(struct arclmframe *af)
	/* gosei:first, t1:first*/
	/* DYNAMIC RESPONSE ANALYSIS. */
	/* NODE MASSES ARE ALREADY SUMMED IN EXTRACTION. */
	/* GRAVITY LOADED ANALYSIS IS ALREADY DONE BY ARCLM. */
	{
	DWORD memory0, memory1, memory2;

	FILE *fin, *fout, *ftxt, *ftxt0, *ftxt1, *ferr,*ftest; /* FILE 8 BYTES */
	double *ddisp, *dreact, emax, ratio0, ratio1;
	struct memoryelem *melem;
	char dir[] = DIRECTORY; /* DATA DIRECTORY */
	char s[80], string[4000], name1[20], name2[20];
	int i, j, ii, jj, iii, jjj;
	int n11, n22, n33, n44, n55, n66, n77, n88, n99, n100;
	int nnode, nelem, nsect, nreact,nlen;
	long int nlap, laps;
	long int loffset, msize;
	long int time;

	struct gcomponent ginit = {0, 0, 0.0, NULL};

	struct gcomponent *gmtx, *g, *p; /* GLOBAL MATRIX [K] */

	double **drccos, **tmatrix, **estiff, *estress, **e, **t;
	double *gdisp, *edisp;
	double determinant, func[2];
	clock_t t0, t1, t2;
	struct onode *nodes, *mnode;
	struct owire elem, *elems;
	struct oconf *confs = NULL;
	struct osect esect;
	double *nmass;

	FILE *fxacc, *fyacc, *fzacc;
	long int ndata[3];
	double ddt, dacc[3], tacc[3], tvel[3], tdis[3], afact = 1.0;
	double pacc[3], pvel[3], pdis[3];
	double *gacc;
	double *u, *ud, *udd, *du, *dud, *dudd; /* VECTORS */
	double Wkt, Wet, Wpt, Wot, *wk[3], *wo[3]; /* ENERGY */
	struct accdata acc[3];
	struct gcomponent *mmtx; /* MASS MATRIX [M] */

	/* int neig=NEIGEN,ncount; */
	int neig = 1, ncount; /* T1 */
	/* int neig=5,ncount; */                  /* T1,T2,---,Tn */


	double **evct, *eigen; /* EIGEN VECTORS,EIGEN VALUES */
	double eps = 1.0E-8, data, d, xdisp1, xdisp2, xdispgap, maxxdisp = 0;
	double h1 = DAMPING, w1, T1, Ti, T11,T12, ene1, ene2;
	struct gcomponent *cmtx,*mtx1,*mtx2; /* MATRIX [C] */
	double mtotal;
	int dtype = DAMPINGTYPE;

	struct snode *sn;



    double *ncr; /*BCLNG CONDENSATION RESULT*/      //ujioka
	int bclngcondensation;

    double gdata;

    bclngcondensation=BCLNGCONDENSATION;

	memory0 = availablephysicalmemory("INITIAL:"); /* MEMORY AVAILABLE */

	/* INPUT FILE */
	fin = fgetstofopen(dir, "r", ID_INPUTFILE);
	/* fin=NULL; */

	/* GET YIELD SURFACE */
	if (fin != NULL) {
		inputinit(fin, &(af->nnode), &(af->nelem), &(af->nsect));
		for (i = 0; i < (af->nsect); i++)
			readsect(fin, (af->sects + i));
		fclose(fin);
		fin = NULL;
	}


	fout = fgetstofopen("\0", "w", ID_OUTPUTFILE);

	globalfile = fout;

	ferr = fopen("gnshn.txt", "w"); /* ERROR FILE */
	if (ferr != NULL) {
		fprintf(ferr, "   LAP     TIME");
		fprintf(ferr,"      ELEM          Nzi          Qyi          Qyj          Mxi          Mxj");
		fprintf(ferr,"      ELEM          Nzi          Qyi          Qyj          Mxi          Mxj");
		fprintf(ferr,"      ELEM          Nzi          Qyi          Qyj          Mxi          Mxj");
		fprintf(ferr,"      ELEM          Nzi          Qyi          Qyj          Mxi          Mxj");
		fprintf(ferr,"      ELEM          Nzi          Qyi          Qyj          Mxi          Mxj");
		fprintf(ferr,"      ELEM          Nzi          Qyi          Qyj          Mxi          Mxj");
		fprintf(ferr, "      NODE           Dx");
		fprintf(ferr, "      NODE           Dx");
		fprintf(ferr, "      NODE           Dx");
		fprintf(ferr, "\n");
	}

	ftxt = fopen("gnshn.huk", "w"); /* HUKUGENRYOKU GRAPH */
	if (ftxt != NULL) {
		fprintf(ftxt, "   LAP");
		for (i = 0; i < FLOORS; i++)
			fprintf(ftxt, "         Dx%2d         Dy%2d", i + 1, i + 1);
		for (i = 0; i < FLOORS - 1; i++)
			fprintf(ftxt, "         Qx%2d         Qy%2d", i + 1, i + 1);
		fprintf(ftxt, "      Accx      Accy");
		fprintf(ftxt, "        Vx        Vy");
		fprintf(ftxt, "        Dx        Dy");
		fprintf(ftxt, "          Wk       Wk+We    Wk+We+Wp");
		fprintf(ftxt, "\n");
	}

	/*edit*/
	/* ENERGY INFO by KANAZAWA */
	nlen = strlen((wdraw.childs+1)->inpfile);
	strcpy(name1, (wdraw.childs+1)->inpfile);
	name1[nlen-3] = 'e';
	name1[nlen-2] = 'e';
	name1[nlen-1] = 'n';
	name1[nlen-0] = '\0';
	ftxt0 = fopen(name1, "w");

	strcpy(name2, (wdraw.childs+1)->inpfile);
	name2[nlen-3] = 't';
	name2[nlen-2] = 'e';
	name2[nlen-1] = 'n';
	name2[nlen-0] = '\0';
	ftxt1 = fopen(name2, "w");
	if (ftxt1 != NULL) {
	  fprintf(ftxt1, "Nstep Energy Check\n");
	}



	t0 = clock();

	/* ACCCELERATION DATA INPUT */
	/* change x wave*/
	fxacc = fopen("el1220w04.txt","r");   /*ñÕã[ínêkîgÅFEL CENTRO*/
	fyacc = NULL;
	fzacc = NULL;

	acc[0].a = inputacceleration(fxacc, &(acc[0].ndata), &(acc[0].dt));
	acc[1].a = inputacceleration(fyacc, &(acc[1].ndata), &(acc[1].dt));
	acc[2].a = inputacceleration(fzacc, &(acc[2].ndata), &(acc[2].dt));


	/* change time increment*/
	ddt = 0.020; /* TIME INCREMENT[sec] */

	/* change strongth of earthquake*/
	afact = 125.0/100.061; /*ñÕã[ínêkîgÅFEL CENTRO*/ //BezierDome

	ene1 = 0.0;
	ene2 = 0.0;

	afact *= 1.0; /* LOCATION FACTOR */
	/* afact*=1.25; */ /* LOCATION FACTOR */

	afact *= 0.01; /* ACCELERATION FACTOR [cm/sec2] INTO [m/sec2] */

	/* COUNT TOTAL LAPS */
	for (i = 0; i < 3; i++) {
		if (acc[i].a == NULL)
		{
		  ndata[i] = 0;
		}
		else
		{
		  ndata[i] = (long int)((acc[i].ndata - 1) * (acc[i].dt) / ddt) + 1;
		}
	}
	laps = 0;
	if (ndata[0] > 0)
		laps = ndata[0];
	else if (ndata[1] > 0)
		laps = ndata[1];
	else if (ndata[2] > 0)
		laps = ndata[2];
	if (ndata[1] > 0 && laps > 0 && ndata[1] < laps)
		laps = ndata[1];
	if (ndata[2] > 0 && laps > 0 && ndata[2] < laps)
		laps = ndata[2];

	if (laps > MAXLAP)
		laps = MAXLAP;
	sprintf(string, "LAPS=%d", laps);

	nnode = af->nnode;
	nelem = af->nelem;
	nsect = af->nsect;
	nreact = af->nreact;

	nodes = af->nodes;
	elems = af->elems;
	confs = af->confs;
	nmass = af->nmass;

	/* x num of 1node */
	xdisp1 = (af->nodes+71)->d[0];


	sprintf(string, "NODES=%d ELEMS=%d SECTS=%d", nnode, nelem, nsect);
	errormessage(string);
	if (fout != NULL)
		fprintf(fout, "%s\n", string);

	msize = 6 * nnode; /* SIZE OF GLOBAL MATRIX. */

	ncount = 0;
	for (iii = 0; iii < msize; iii++) {
		if ((confs + iii)->iconf != 1)
			ncount++;
	}
	sprintf(string, "MATRIX DIMENSION=%d", ncount);

	gmtx = (struct gcomponent*)malloc(msize*sizeof(struct gcomponent));
	for (i = 0; i < msize; i++)
	{
	  (gmtx+i)->down = NULL; /* GLOBAL MATRIX INITIALIZATION. */
	}
	mmtx = (struct gcomponent*)malloc(msize*sizeof(struct gcomponent));
	for (i = 0; i < msize; i++) {
		(mmtx + i)->m = (unsigned short int)(i + 1);
		(mmtx + i)->n = (unsigned short int)(i + 1);
		(mmtx + i)->value = 0.0;
		(mmtx + i)->down = NULL;
	}
#if 1 /*FOR PRORORTION TO INITIAL STIFFNESS, OCCASIONAL STIFFNESS II.*/
	if (dtype == 1 || dtype == 3)
	{
		cmtx = (struct gcomponent*)malloc(msize*sizeof(struct gcomponent));
		for (i = 0; i < msize; i++)
		{
		  (cmtx + i)->down = NULL; /* MATRIX [C] INITIALIZATION. */
		}
	}
#endif

	evct = (double **)malloc(neig*sizeof(double *)); /* EIGEN VECTORS */
	eigen = (double *)malloc(neig*sizeof(double)); /* EIGEN VALUES */
	if (evct == NULL || eigen == NULL)return 0;
	for (i = 0; i < neig; i++)
	{
	  *(evct + i) = (double *)malloc(msize*sizeof(double));
	  for (j = 0; j < msize; j++)
	  {
		*(*(evct + i) + j) = 0.0;
	  }
	}

	/* DISPLACEMENT VECTORS. */
	u = mallocdoublevector(msize);
	ud = mallocdoublevector(msize);
	udd = mallocdoublevector(msize);
	du = mallocdoublevector(msize);
	dud = mallocdoublevector(msize);
	dudd = mallocdoublevector(msize);

	gacc = mallocdoublevector(msize); /* ACCEL VECTOR */

	for (i = 0; i < msize; i++)
	{
		*(u + i)    = 0.0;
		*(ud + i)   = 0.0;
		*(udd + i)  = 0.0;
		*(du + i)   = 0.0;
		*(dud + i)  = 0.0;
		*(dudd + i) = 0.0;
	}


	ddisp = (double *)malloc(6 * nnode*sizeof(double));
	melem = (struct memoryelem*)malloc(nelem*sizeof(struct memoryelem));

	af->ddisp = ddisp; /* DISPLACEMENT:6 DIRECTIONS. */
	af->melem = melem; /* CODE,12 BOUNDARIES,12 STRESS. */


	initialform(nodes, ddisp, nnode); /* ASSEMBLAGE FORMATION. */

	mnode = (struct onode*)malloc(nnode*sizeof(struct onode));
	for (i = 0; i < nnode; i++) {
		*(mnode + i) = *(nodes + i);
	}

	initialelem(elems, melem, nelem); /* ASSEMBLAGE ELEMENTS. */


	dreact = (double *)malloc(nreact*sizeof(double)); /* REACTION FILE. */
	af->dreact = dreact;
	initialreact(fin, dreact, nreact); /* ASSEMBLAGE LONG REACTIONS. */



	drawglobalaxis((wdraw.childs + 1)->hdcC, (wdraw.childs + 1)->vparam, 0, 0,
		255); /* DRAW GLOBAL AXIS. */

	if (wsurf.hwnd != NULL) {
		drawyieldsurface((wsurf.childs + 1)->hdcC, (wsurf.childs + 1)->vparam,
			SURFACEX, SURFACEY, SURFACEZ, NULL);
		overlayhdc(*(wsurf.childs + 1), SRCPAINT); /* UPDATE DISPLAY. */
	}
	af->fsurface = fopen("canbin.sfc", "wb+"); /* STRESS ON SURFACE. */

	/* ELEMENT ENERGY WITHOUT GRAVITY. */
	for (i = 0; i < nelem; i++) {
		(elems + i)->Ee[0] = 0.0;
		(elems + i)->Ee[1] = 0.0;
		(elems + i)->Ep[0] = 0.0;
		(elems + i)->Ep[1] = 0.0;
	}

	/* GLOBAL ENERGY. */
	wk[0] = (double *)malloc(nnode*sizeof(double));
	wk[1] = (double *)malloc(nnode*sizeof(double));
	wk[2] = (double *)malloc(nnode*sizeof(double));
	wo[0] = (double *)malloc(nnode*sizeof(double));
	wo[1] = (double *)malloc(nnode*sizeof(double));
	wo[2] = (double *)malloc(nnode*sizeof(double));
	for (ii = 0; ii < 3; ii++) {
		for (jj = 0; jj < nnode; jj++) {
			*(wk[ii] + jj) = 0.0;
			*(wo[ii] + jj) = 0.0;
		}
	}

	for (ii = 0; ii < 3; ii++) {
		tacc[ii] = 0.0;
		tvel[ii] = 0.0;
		tdis[ii] = 0.0;
	}

	/* ASSEMBLAGE MASS MATRIX. */
	mtotal = assemmass(mmtx, nnode, confs, nmass);
	if (fout != NULL)fprintf(fout, "TOTAL MASS = %12.5f\n\n", mtotal);


	/* MEMORIES. */
	drccos = mallocdoublematrix(3);
	estiff = mallocdoublematrix(12);
	tmatrix = mallocdoublematrix(12);
	e = mallocdoublematrix(12);
	t = mallocdoublematrix(12);
	estress = mallocdoublevector(12);
	gdisp = mallocdoublevector(12);
	edisp = mallocdoublevector(12);

    /* BUCKLING LOADS */
    ncr=mallocdoublevector(nelem);

	definencr(&arc,&*ncr);          //ujioka

#if 0 /*FOR PRORORTION TO INITIAL STIFFNESS. (DAMPING MATRIX (dtype=1))*/
	if (dtype == 1)
	{
		for (i = 1; i <= msize; i++) /* INITIALIZATION. */
		{
			g = (cmtx + (i - 1))->down; /* NEXT OF DIAGONAL. */
			while (g != NULL) /* CLEAR ROW. */
			{
				p = g;
				g = g->down;
				free(p);
			}
			ginit.m = (unsigned short int)i;
			*(cmtx + (i - 1)) = ginit;
		}
		for (i = 0; i < nelem;i++) /* TEST */      /* ASSEMBLAGE DAMPING MATRIX. */
		{
			inputelem(elems, melem, i, &elem); /* READ ELEMENT DATA. */

			for (ii = 0; ii <= 1; ii++)
			{
				for (jj = 0; jj < 6; jj++)
				{
					(elems + i)->iconf[ii][jj] = elem.iconf[ii][jj];
				}
			}

			inputnode(ddisp, elem.node[0]); /* HEAD */
			inputnode(ddisp, elem.node[1]); /* TAIL */

			/* elem.sect=(elems+i)->sect; */
			/* READ SECTION DATA. */
			esect = *((elems + i)->sect);
			elem.sect = &esect;
			/* elem.sect->area=0.0; */

			directioncosineII  (elem.node[0]->d[0],
								elem.node[0]->d[1],
								elem.node[0]->d[2],
								elem.node[1]->d[0],
								elem.node[1]->d[1],
								elem.node[1]->d[2],
								elem.cangle,
								drccos); /* [DRCCOS] */

			transmatrixII(drccos, tmatrix); /* TRANSFORMATION MATRIX. */
			assememtxII(elem, estiff); /* ELASTIC MATRIX OF ELEMENT. */
			estiff = modifyhinge(elem, estiff); /* MODIFY MATRIX. */
			/* estiff=assempmtx(elem,estiff); */      /* ADD PLASTIC MATRIX. */
			transformationII(estiff, tmatrix, e, t); /* [K]=[Tt][k][T] */
			assemgstiffness(cmtx, estiff, &elem); /* ASSEMBLAGE. */
		}

		/*CHANGE T1*/
		T1 = 1.77364112 ;/* FUNDAMENTAL NATURAL PERIOD 1éüå≈óLé¸ä˙ */


		w1 = 2.0 * PI / T1;
		for (ii = 1; ii <= msize; ii++) {
			for (jj = 1; jj <= ii; jj++) {
				gread(cmtx, ii, jj, &data);
				data *= (2.0*h1 / w1);
				gwrite(cmtx, ii, jj, data);
			}
		}
		sprintf(string, "DAMPING MATRIX %ld COMPS ASSEMBLED.", comps);
		laptime(string, t0);
	}
#endif


	int n101 = 750;
	/*changeanalysis*/
	bool FLAG2 = 0; /* 0:bisecsylvester  , 1:deigabgeneral */



	/* Ç±Ç±Ç©ÇÁêUìÆâêÕÇÃåJÇËï‘ÇµââéZ */
	for (nlap = 1; nlap <= n101 ; nlap++)
	{ /*Å@n11Ç≈ã≠êßèIóπ or lapsÇ‹Ç≈Å@*/

		/* af->nlaps=nlap; */
		af->nlaps = 1;

		sprintf(string, "LAP:%d/%d", nlap, laps);
		errormessage(string);
		if (fout != NULL)
			fprintf(fout, "%s\n", string);
		if (ferr != NULL)
			fprintf(ferr, " %5d %8.3f", nlap, nlap*ddt);
		if (ftxt != NULL)
			fprintf(ftxt, " %5d", nlap);

		setlaps((wmenu.childs + 2)->hwnd, nlap, laps);

		memory1 = availablephysicalmemory("REMAIN:"); /* MEMORY AVAILABLE */

		for (i = 1; i <= msize; i++) /* GLOBAL MATRIX INITIALIZATION. */
		{
			g = (gmtx + (i - 1))->down; /* NEXT OF DIAGONAL. */
			while (g != NULL) /* CLEAR ROW. */ {
				p = g;
				g = g->down;
				free(p);
			}

			ginit.m = (unsigned short int)i;
			/* ginit.n=(unsigned short int)i; */

			*(gmtx + (i - 1)) = ginit;
		}
		comps = msize; /* INITIAL COMPONENTS=DIAGONALS. */

		laptime("ASSEMBLING GLOBAL MATRIX.", t0);

		clearwindow(*(wdraw.childs + 1));
		/* drawarclmtitle((wdraw.childs+1)->hdcC,
		 (wdraw.childs+1)->vparam,*af,0,ONSCREEN); */
		for (i = 0; i < (wdraw.childs + 1)->org.ntext; i++) {
			sn = (wdraw.childs + 1)->org.texts + i;
			wdraw.strset = addtext(wdraw.strset, &(wdraw.nstring), sn->str,
				sn->n.d[0], sn->n.d[1]);
		}
		SetTextColor((wdraw.childs + 1)->hdcC, RGB(255, 255, 255));
		drawtexts((wdraw.childs + 1)->hdcC, wdraw.strset, wdraw.nstring,
			(wdraw.childs + 1)->vparam);

		for (i = 0; i < nelem; i++) /* ASSEMBLAGE GLOBAL MATRIX. */ {
			inputelem(elems, melem, i, &elem); /* READ ELEMENT DATA. */

			for (ii = 0; ii <= 1; ii++) {
				for (jj = 0; jj <= 5; jj++) {
					(elems + i)->iconf[ii][jj] = elem.iconf[ii][jj];
				}
			}

			inputnode(ddisp, elem.node[0]); /* HEAD */
			inputnode(ddisp, elem.node[1]); /* TAIL */

			elem.sect = (elems + i)->sect; /* READ SECTION DATA. */

			if ((wdraw.childs + 1)->hdcC != NULL) /* DRAW DEFORMED ELEMENT. */ {
				drawglobalwire((wdraw.childs + 1)->hdcC,
					(wdraw.childs + 1)->vparam, *af, elem, 255, 255, 255, 255,
					255, 255, 0, ONSCREEN);
			}

			directioncosineII  (elem.node[0]->d[0],
								elem.node[0]->d[1],
								elem.node[0]->d[2],
								elem.node[1]->d[0],
								elem.node[1]->d[1],
								elem.node[1]->d[2],
								elem.cangle,
								drccos); /* [DRCCOS] */

			transmatrixII(drccos, tmatrix); /* TRANSFORMATION MATRIX. */

			assememtxII(elem, estiff); /* ELASTIC MATRIX OF ELEMENT. */

			estiff = modifyhinge(elem, estiff); /* MODIFY MATRIX. */


			if(bclngcondensation) /***UJIOKA***/
			{
			  estiff = assempmtxbc(elem,estiff,ncr[i]);/* ADD PLASTIC MATRIX. */
			}
			else
			{
			  estiff = assempmtx(elem, estiff); /* ADD PLASTIC MATRIX. */
			}

			transformationII(estiff, tmatrix, e, t); /* [K]=[Tt][k][T] */
			assemgstiffness(gmtx, estiff, &elem); /* ASSEMBLAGE. */
		}



		sprintf(string, "GLOBAL MATRIX %ld COMPS ASSEMBLED.", comps);
		laptime(string, t0);

		overlayhdc(*(wdraw.childs + 1), SRCPAINT); /* UPDATE DISPLAY. */

		for (i = 0; i < 3; i++)
		{
			if (acc[i].a == NULL)
			{
				dacc[i] = 0.0;
			}
			else
			{
				dacc[i] = accelerationincrement(acc[i], ddt, nlap);
				dacc[i] *= afact;
			}

			pacc[i] = tacc[i];
			pvel[i] = tvel[i];
			pdis[i] = tdis[i];

			tacc[i] += dacc[i]; /* GROUND ACCELERATION. */
			tvel[i] += ddt * (pacc[i] + tacc[i]) / 2.0;
			/* TOTAL GROUND VELOCITY. */
			/* tdis[i]+=ddt*(pdis[i]+tvel[i])/2.0; */ /* ERROR */
			/* GROUND DISPLACEMENT. */
			tdis[i] += ddt * (pvel[i] + tvel[i]) / 2.0; /* GROUND DISPLACEMENT. */
		}
		if (fout != NULL)
			fprintf(fout, "ACCELERATION=%.8f,%.8f,%.8f\n", dacc[0], dacc[1],
			dacc[2]);

		assemaccel(gacc, dacc, nnode, af); /* ACCELERATION VECTOR. */


		/* CALCULATE EIGEN PERIOD å≈óLé¸ä˙âêÕ */
		if (nlap == 1)
        {
			if (MessageBox(NULL,"Calculate Eigen Period.","GNSHN",MB_YESNO) == IDYES)
            {
				mtx1 = copygcompmatrix(gmtx, msize);
				mtx2 = copygcompmatrix(mmtx, msize);

                for(ii=1;ii<=msize;ii++)
                {
                   sprintf(string,"%3d",ii);
                   for(jj=1;jj<=msize;jj++)
                   {
					 gread(mtx1,ii,jj,&gdata);
				   }

				}
				for(ii=1;ii<=msize;ii++)
                {
                   sprintf(string,"%3d",ii);
                   for(jj=1;jj<=msize;jj++)
                   {
					 gread(mtx2,ii,jj,&gdata);
				   }
				}


				if (FLAG2) { // deigabgeneral
					MessageBox(NULL, "DEIGABGENERAL BEGIN.", "Gnshn", MB_OK);
					deigabgeneral(mtx2, mtx1, confs, msize, neig, neig, eps,
						eigen, evct);
					laptime("EIGEN COMPLETED.", t0);
					MessageBox(NULL, "DEIGABGENERAL END.", "Gnshn", MB_OK);
				}
				else { // bisecsylvester
					MessageBox(NULL, "BISECSYLVESTER BEGIN.", "Gnshn", MB_OK);
					bisecsylvester(mtx2, mtx1, confs, msize, neig, neig, eps, eigen, evct);
					laptime("EIGEN COMPLETED.", t0);
					MessageBox(NULL, "BISECSYLVESTER END.", "Gnshn", MB_OK);
				}

				if (fout != NULL)
					fprintf(fout, "\n");
				for (i = 0; i < neig; i++)
                {
					if (FLAG2)
                    {
						if (fout != NULL)
							fprintf(fout,
							"DEIGABGENERAL EIGEN VALUE %ld=%.8f\n", (i + 1),
							*(eigen + i));
						if (*(eigen + i) > 0.0)
						{
							Ti = 2.0 * PI * sqrt(*(eigen + i));
							sprintf(string, "PERIOD T%ld=%.8f [sec]",
								(i + 1), Ti);
							fprintf(fout, "%s\n", string);
						}
						else
						{
							fprintf(fout, "ERROR:EIGEN VALUE NEGATIVE.\n");
						}
					}
					else // bisecsylvester
					{
						if (fout != NULL)
							fprintf(fout,
							"BISECSYLVESTER EIGEN VALUE %ld=%.8f\n", (i + 1),
							1.0 / (*(eigen + i))); // bisecsylvester
						if (*(eigen + i) > 0.0)
                        {
							Ti = 2.0 * PI * sqrt(*(eigen + i));
							sprintf(string, "PERIOD T%ld=%.8f [sec]",
								(i + 1), Ti);
							fprintf(fout, "%s\n", string);
						}
						else
                        {
							fprintf(fout, "ERROR:EIGEN VALUE NEGATIVE.\n");
						}
					}

					fprintf(fout, "\nEIGEN VECTOR %ld\n",(i + 1));
					for (ii = 0; ii < nnode; ii++)
					{
							fprintf(fout,
							"%4ld %11.7f %11.7f %11.7f %11.7f %11.7f %11.7f\n",
							(nodes + ii)->code, *(*(evct + i) + 6*ii + 0),
							*(*(evct + i) + 6*ii + 1),
							*(*(evct + i) + 6*ii + 2),
							*(*(evct + i) + 6*ii + 3),
							*(*(evct + i) + 6*ii + 4),
							*(*(evct + i) + 6*ii + 5));
					}

                    if(i==0) T1=Ti;

				}
				//fclose(fout);

//				if (*(eigen + 0) <= 0.0)
//					return 1;
//				w1 = sqrt(*(eigen + 0));
				/* w1=sqrt(E[0]); */
//				T1 = 2.0 * PI / w1;
//				fclose(fout);
//				return 1;
			} /* #endif */    /* å≈óLé¸ä˙âêÕÇ±Ç±Ç‹Ç≈ */

			if (dtype != 1)
            {
				/* T1=0.24422; */ /* [sec] */
				/* T1=0.78; */ /* FOR E-DEFENSE*[sec] */
				/* T1=0.479; */ /* FOR HAKONE[sec] */
				/* T1=0.8; */ /* FOR HAKONE UD[sec] */
				/* T1=7.464; */ /* FOR HAKONE[sec] */
				/* T1=3.96572; */ /* FOR KAMISATO[sec] */
				/* T1=0.87492; */ /* FOR KYOBASHI LINEAR[sec] */
				/* T1=0.70795; */ /* FOR KYOBASHI NONLINEAR[sec] */
				/* T1=0.74210; */ /* FOR HAKUSHIMA Y */
				/* T1=0.19261; */ /* FOR HAKUSHIMA X */
				T1 = 0.823; /* FOR TOSHIMA*[sec] */
			}
			w1 = 2.0 * PI / T1;

#if 1
			if (dtype == 1)
			{
				for (i = 1; i <= msize; i++) /* INITIALIZATION. */
				{
					g = (cmtx + (i - 1))->down; /* NEXT OF DIAGONAL. */
					while (g != NULL) /* CLEAR ROW. */
					{
						p = g;
						g = g->down;
						free(p);
					}
					ginit.m = (unsigned short int)i;
					*(cmtx + (i - 1)) = ginit;
				}
				for (i = 0; i < nelem;i++) /* TEST */      /* ASSEMBLAGE DAMPING MATRIX. */
				{
					inputelem(elems, melem, i, &elem); /* READ ELEMENT DATA. */

					for (ii = 0; ii <= 1; ii++)
					{
						for (jj = 0; jj <= 5; jj++)
						{
							(elems + i)->iconf[ii][jj] = elem.iconf[ii][jj];
						}
					}

					inputnode(ddisp, elem.node[0]); /* HEAD */
					inputnode(ddisp, elem.node[1]); /* TAIL */

					/* elem.sect=(elems+i)->sect; */
					/* READ SECTION DATA. */
					esect = *((elems + i)->sect);
					elem.sect = &esect;
					/* elem.sect->area=0.0; */

					directioncosineII(elem.node[0]->d[0], elem.node[0]->d[1],
						elem.node[0]->d[2], elem.node[1]->d[0], elem.node[1]->d[1],
						elem.node[1]->d[2], elem.cangle, drccos); /* [DRCCOS] */

					transmatrixII(drccos, tmatrix); /* TRANSFORMATION MATRIX. */
					assememtxII(elem, estiff); /* ELASTIC MATRIX OF ELEMENT. */
					estiff = modifyhinge(elem, estiff); /* MODIFY MATRIX. */
					/* estiff=assempmtx(elem,estiff); */      /* ADD PLASTIC MATRIX. */
					transformationII(estiff, tmatrix, e, t); /* [K]=[Tt][k][T] */
					assemgstiffness(cmtx, estiff, &elem); /* ASSEMBLAGE. */
				}


				w1 = 2.0 * PI / T1;
				for (ii = 1; ii <= msize; ii++)
				{
					for (jj = 1; jj <= ii; jj++)
					{
						gread(cmtx, ii, jj, &data);
						data *= (2.0*h1 / w1);
						gwrite(cmtx, ii, jj, data);
					}
				}
				sprintf(string, "DAMPING MATRIX %ld COMPS ASSEMBLED.", comps);
				laptime(string, t0);
            }
        }
#endif



			fprintf(fout, "T1=%.5f w1=%.5f h1=%.5f\n", T1, w1, h1);

		/* NEWMARK'S BETA PROCESS. */
		laptime("NEWMARK'S BETA.", t0);
		determinant = newmarkbeta(gmtx, cmtx, mmtx, gacc, confs, nnode, u, ud,
			udd, du, dud, dudd, ddt, (1.0 / 4.0));

		sprintf(string, "DETERMINANT=%.5E COMPS=%ld", determinant, comps);
		errormessage(string);
		fprintf(fout, "%s\n", string);


		if (determinant <= 0.0)
		{
			errormessage(" ");
			errormessage("INSTABLE TERMINATION.");
			if (fout != NULL)
				fprintf(fout, "INSTABLE TERMINATION.\n");

			laptime("\0", t0);

			fclose(fin);

			gfree(gmtx, nnode); /* FREE GLOBAL MATRIX. */
			free(gacc);
			/* free(confs); */
			freematrix(drccos, 3);
			freematrix(estiff, 12);
			freematrix(tmatrix, 12);
			freematrix(e, 12);
			freematrix(t, 12);
			free(estress);
			free(gdisp);
			free(edisp);

			memory2 = availablephysicalmemory("REMAIN:");
			sprintf(string, "CONSUMPTION:%ld[BYTES]", (memory0 - memory2));
			errormessage(string);

			return 1;
		}

		laptime("OUTPUT INTO FILE.", t0);

		if (fout != NULL)
			fprintf(fout, "\"DISPLACEMENT\"\n");
		outputdisp(du, fout, nnode, nodes); /* INCREMENTAL DISPLACEMENT. */


		for (i = 0; i < nelem; i++) /* STRESS OUTPUT,UPDATE. */
		{
			inputelem(elems, melem, i, &elem);

			inputnode(ddisp, elem.node[0]);
			inputnode(ddisp, elem.node[1]);

			elem.sect = (elems + i)->sect; /* READ SECTION DATA. */

          if (bclngcondensation) /***UJIOKA***/
          {
            if(elem.stress[0][0]>0)  //compression
            {
  	          sprintf(string,"ELEM %d :Ncr=%5.8f\n",
						(af->elems+i-1)->code,ncr[i]);

              elemstressIIbc(estress, &elem, du, melem, fout, drccos, tmatrix,
				estiff, gdisp, edisp, func,ftxt,ncr[i]);
            }
            else                     //tension
            {
              elemstressIIbc(estress, &elem, du, melem, fout, drccos, tmatrix,
				estiff, gdisp, edisp, func,ftxt,0);
              /*
				ï÷ãXè„Ncr=0Ç∆ÇµÅAÇ±ÇÃèÍçáç~ïöã»ñ ÇÕèCê≥ÇµÇ»Ç¢ÅB
				(updatestressbc,coefficientsbcÇ≈ÇÃèåèï™äÚÇ…ÇÊÇÈ)
			  */
            }

          }
          else
          {
			elemstressII(estress, &elem, du, melem, fout, drccos, tmatrix,
				estiff, gdisp, edisp, func, ftxt);
          }
			(elems + i)->Ee[0] = elem.Ee[0];
			(elems + i)->Ee[1] = elem.Ee[1];
			(elems + i)->Ep[0] = elem.Ep[0];
			(elems + i)->Ep[1] = elem.Ep[1];

		}

		/* UPDATE DISPLAY. */
		if (wsurf.hwnd != NULL)
		{
			drawyieldsurface((wsurf.childs + 1)->hdcC,
				(wsurf.childs + 1)->vparam, SURFACEX, SURFACEY, SURFACEZ,
				af->fsurface);
			overlayhdc(*(wsurf.childs + 1), SRCPAINT);

		}





		fprintf(fout, "\"REACTION\"\n");
		outputreaction(gmtx, du, nodes, confs, dreact, fout, nnode);

		updateform(ddisp, du, nnode); /* FORMATION UPDATE. */

		if (fout != NULL)fprintf(fout, "\"CURRENT FORM\"\n\n");
		for (ii = 0; ii < nnode; ii++) {


			/*
			if (0) {

				sprintf(string, "NODE:%5ld(%5ld) {U}=", (nodes + ii)->code, ii);
				for (jj = 0; jj < 6; jj++) {
					loffset = 6 * ii + jj;
					sprintf(s, " %14.5f", *(ddisp + loffset));
					strcat(string, s);
				}
				if (fout != NULL)
					fprintf(fout, "%s\n", string);

				fprintf(ferr,
					" NODE %9d %12.5f %12.5f %12.5f %12.8f %12.8f %12.8f",
					(nodes + ii)->code,
					*(ddisp + 6*ii + 0) - (mnode + ii)->d[0],
					*(ddisp + 6*ii + 1) - (mnode + ii)->d[1],
					*(ddisp + 6*ii + 2) - (mnode + ii)->d[2], *(u + 6*ii + 0),
					*(u + 6*ii + 1), *(ud + 6*ii + 0), *(ud + 6*ii + 1),
					*(udd + 6*ii + 0), *(udd + 6*ii + 1));



			}
			*/
		}

		fprintf(ferr, "\n");

		/* CALCULATE FLOOR DISPLACEMENT,FLOOR HUKUGENRYOKU. */

		/* floorvalues(fout,ftxt,af,ddisp,dacc); */
		floorvalues(NULL, ftxt, af, ddisp, dacc);
		/* if(ftxt!=NULL) fprintf(ftxt,"\n"); */

		/* SUM GLOBAL ENERGY. */
		Wkt = 0.0;
		Wet = 0.0;
		Wpt = 0.0;
		Wot = 0.0;
		for (ii = 0; ii < nnode; ii++) {
			for (jj = 0; jj < 3; jj++) {
				*(wk[jj] + ii) = 0.5 * (*(nmass + ii)) * (*(ud + 6 * ii + jj)) *
					(*(ud + 6 * ii + jj));
				Wkt += *(wk[jj] + ii);

				*(wo[jj] + ii) -=
					0.5 * (*(nmass + ii)) * (2 * tacc[jj] - dacc[jj]) *
					(*(du + 6 * ii + jj));
				Wot += *(wo[jj] + ii);
			}
		}

		for (ii = 0; ii < nelem; ii++) {
			Wet += (elems + ii)->Ee[0];
			Wet += (elems + ii)->Ee[1];
			Wpt += (elems + ii)->Ep[0];
			Wpt += (elems + ii)->Ep[1];
		}

		if (ftxt != NULL) {
			fprintf(ftxt, " %9.3f %9.3f", tacc[0], tacc[1]); /* Acc */
			fprintf(ftxt, " %9.3f %9.3f", tvel[0], tvel[1]); /* V */
			fprintf(ftxt, " %9.3f %9.3f", tdis[0], tdis[1]); /* D */
			fprintf(ftxt, " %11.7f %11.7f %11.7f", Wkt, Wkt + Wet,
				Wkt + Wet + Wpt); /* W */
			fprintf(ftxt, "\n");
		}

		xdisp2 = (af->nodes+71)->d[0];
		xdispgap = abs(xdisp2 - xdisp1);
		if (xdispgap > maxxdisp)
		{
			maxxdisp = xdispgap;
		}

		if (ftxt1 != NULL) {
			fprintf(ftxt1, "%11.7f,%11.7f\n", Wpt,Wot);
		}


		t1 = laptime("\0", t0);

		memory2 = availablephysicalmemory(NULL);
		sprintf(string, "CONSUMPTION:%ld[BYTES]", (memory1 - memory2));
		errormessage(string);

		errormessage("L:CONTINUE R:ABORT"); /* L=LEFT R=RIGHT */
		while (!GetAsyncKeyState(VK_LBUTTON)) /* LEFT CLICK TO CONTINUE. */ {
			if (GetAsyncKeyState(VK_RBUTTON)) /* RIGHT CLICK TO ABORT. */ {
				fclose(fin);

				gfree(gmtx, nnode); /* FREE GLOBAL MATRIX. */
				gfree(mmtx, nnode); /* FREE MASS   MATRIX. */
				free(gacc);
				/* free(confs); */
				freematrix(drccos, 3);
				freematrix(estiff, 12);
				freematrix(tmatrix, 12);
				freematrix(e, 12);
				freematrix(t, 12);
				free(estress);
				free(gdisp);
				free(edisp);

				errormessage(" ");
				errormessage("ABORTED.");
				if (fout != NULL)
					fprintf(fout, "ABORTED.\n");

				laptime("\0", t0);
				return 1;
			}
			t2 = clock();
			time = (t2 - t1) / CLK_TCK;
			if (time >= WAIT)
				break; /* CONTINUE AFTER WAITING. */
		}

	} /* END OF LAP.REPEAT UNTIL INSTABLE. êUìÆåJÇËï‘ÇµÇ±Ç±Ç‹Ç≈ */


	/*edit*/
		if (ftxt0 != NULL) {
			emax = 0.0;
			for (ii = 0; ii < nelem; ii++) {
				if (emax < ((elems + ii)->Ee[0] + (elems + ii)->Ep[0])) {
					emax = ((elems + ii)->Ee[0] + (elems + ii)->Ep[0]);
				}
				if (emax < ((elems + ii)->Ee[1] + (elems + ii)->Ep[1])) {
					emax = ((elems + ii)->Ee[1] + (elems + ii)->Ep[1]);
				}
			}
			fprintf(ftxt0, "//Check for Energy Distribution//\n\n");
			fprintf(ftxt0, "Wpt Total = %-13.6f\n",Wpt);//ó›êœëYê´òcÉGÉlÉãÉMÅ[
			fprintf(ftxt0, "Wkt Total = %-13.6f\n",Wkt);//íeê´êUìÆÉGÉlÉãÉMÅ[
			fprintf(ftxt0, "Wet Total = %-13.6f\n",Wet);//íeê´òcÉGÉlÉãÉMÅ[
			fprintf(ftxt0, "Wot Total = %-13.6f\n",Wot);//ì¸óÕÉGÉlÉãÉMÅ[
			fprintf(ftxt0, "Wkt+Wet+Wpt Total = %-13.6f\n",Wkt+Wet+Wpt);
			fprintf(ftxt0, "VE = %-13.6f\n\n",sqrt(2*Wot*10000/mtotal));
			sprintf(string, "VE = %-13.6f\n", sqrt(2*Wot*10000/mtotal));
			fprintf(ftxt0, "max xdisp of node172 = %-13.6f\n\n", maxxdisp);
			/*MessageBox(NULL, "BISECSYLVESTER START.", "Gnshn", MB_OK);*/
			/*deigabgeneral(mtx3, mtx4, confs, msize, -neig, neig, eps,
						eigen, evct);
			T11 = 2.0 * PI / sqrt(*(eigen + 0));
			bisecsylvester(mtx4, mtx3, confs, msize, neig, neig, eps,
						eigen, evct);
			T11 = 2.0 * PI * sqrt(*(eigen + 0));*/
            /*
			deigabgeneral(mtx3, mtx4, confs, msize, -neig, neig, eps,
						eigen, evct);
			T11 = 2.0 * PI / sqrt(*(eigen + 0));
			fprintf(ftxt0, "LAP %3d  :  T1  = %-13.6f\n\n",n11,T11);

			deigabgeneral(mtx5, mtx6, confs, msize, -neig, neig, eps,
						eigen, evct);
			T11 = 2.0 * PI / sqrt(*(eigen + 0));
			fprintf(ftxt0, "LAP %3d  :  T1  = %-13.6f\n\n",n22,T11);

			deigabgeneral(mtx7, mtx8, confs, msize, -neig, neig, eps,
						eigen, evct);
			T11 = 2.0 * PI / sqrt(*(eigen + 0));
			fprintf(ftxt0, "LAP %3d  :  T1  = %-13.6f\n\n",n33,T11);

			deigabgeneral(mtx9, mtx10, confs, msize, -neig, neig, eps,
						eigen, evct);
			T11 = 2.0 * PI / sqrt(*(eigen + 0));
			fprintf(ftxt0, "LAP %3d  :  T1  = %-13.6f\n\n",n44,T11);

			deigabgeneral(mtx11, mtx12, confs, msize, -neig, neig, eps,
						eigen, evct);
			T11 = 2.0 * PI / sqrt(*(eigen + 0));
			fprintf(ftxt0, "LAP %3d  :  T1  = %-13.6f\n\n",n55,T11);

			deigabgeneral(mtx13, mtx14, confs, msize, -neig, neig, eps,
						eigen, evct);
			T11 = 2.0 * PI / sqrt(*(eigen + 0));
			fprintf(ftxt0, "LAP %3d  :  T1  = %-13.6f\n\n",n66,T11);

			deigabgeneral(mtx15, mtx16, confs, msize, -neig, neig, eps,
						eigen, evct);
			T11 = 2.0 * PI / sqrt(*(eigen + 0));
			fprintf(ftxt0, "LAP %3d  :  T1  = %-13.6f\n\n",n77,T11);

			deigabgeneral(mtx17, mtx18, confs, msize, -neig, neig, eps,
						eigen, evct);
			T11 = 2.0 * PI / sqrt(*(eigen + 0));
			fprintf(ftxt0, "LAP %3d  :  T1  = %-13.6f\n\n",n88,T11);

			deigabgeneral(mtx19, mtx20, confs, msize, -neig, neig, eps,
						eigen, evct);
			T11 = 2.0 * PI / sqrt(*(eigen + 0));
			fprintf(ftxt0, "LAP %3d  :  T1  = %-13.6f\n\n",n99,T11);

			deigabgeneral(mtx21, mtx22, confs, msize, -neig, neig, eps,
						eigen, evct);
			T11 = 2.0 * PI / sqrt(*(eigen + 0));
			fprintf(ftxt0, "LAP %3d  :  T1  = %-13.6f\n\n",n100,T11);
			*/
			MessageBox(NULL, string, "Gnshn", MB_OK);
			/*for (ii = 0; ii < nelem; ii++) {
				ratio0 = (elems + ii)->Ep[0] / (Wpt);
				ratio1 = (elems + ii)->Ep[1] / (Wpt);
				fprintf(ftxt0, "ELEM %d", (elems + ii)->code);
				fprintf(ftxt0, " SPT %d", (elems + ii)->node[0]->code);
				fprintf(ftxt0, " ENE %-8.6f", (elems + ii)->Ep[0]);
				fprintf(ftxt0, " RATIO %-8.6f", ratio0);
				fprintf(ftxt0, " EPT %d", (elems + ii)->node[1]->code);
				fprintf(ftxt0, " ENE %-8.6f", (elems + ii)->Ep[1]);
				fprintf(ftxt0, " RATIO %-8.6f", ratio1);
				fprintf(ftxt0, "\n");

			}*/

/****Modified by UJIOKA****/
            // ãzé˚ÉGÉlÉãÉMÅ[ÇÃílÇèoóÕ
            for (ii = 0; ii < nelem; ii++)
            {
				ratio0 = ((elems + ii)->Ee[0] + (elems + ii)->Ep[0]) / (Wet+Wpt);
				ratio1 = ((elems + ii)->Ee[1] + (elems + ii)->Ep[1]) / (Wet+Wpt);
				fprintf(ftxt0, "ELEM %d", (elems + ii)->code);
				fprintf(ftxt0, " NODE %d Ee+Ep= %-8.12f + %-8.12f = %-8.12f RATIO: %-8.6f",
                                 (elems + ii)->node[0]->code,
                                 (elems + ii)->Ee[0],
                                 (elems + ii)->Ep[0],
                                 (elems + ii)->Ee[0]+(elems + ii)->Ep[0],
                                 ratio0);
				fprintf(ftxt0, " NODE %d Ee+Ep= %-8.12f + %-8.12f = %-8.12f RATIO: %-8.6f",
                                 (elems + ii)->node[1]->code,
                                 (elems + ii)->Ee[1],
                                 (elems + ii)->Ep[1],
                                 (elems + ii)->Ee[1]+(elems + ii)->Ep[1],
                                 ratio1);
				fprintf(ftxt0, "\n");

			}
/****Modified by UJIOKA****/
	fclose(ftxt0);
		}

	if ((wdraw.childs + 1)->hdcC != NULL && melem != NULL && ddisp != NULL)
		/* DRAW LAST FRAME. */ {
		for (i = 1; i <= nelem; i++) {
			inputelem(elems, melem, i - 1, &elem);
			for (ii = 0; ii <= 1; ii++) /* COPY HINGE DATA. */ {
				for (jj = 0; jj <= 5; jj++) {
					(elems + i - 1)->iconf[ii][jj] = elem.iconf[ii][jj];
				}
			}

			inputnode(ddisp, elem.node[0]);
			inputnode(ddisp, elem.node[1]);

			drawglobalwire((wdraw.childs + 1)->hdcC, (wdraw.childs + 1)->vparam,
				*af, elem, 255, 255, 255, 255, 255, 255, 0, ONSCREEN);
		}
		overlayhdc(*(wdraw.childs + 1), SRCPAINT); /* UPDATE DISPLAY. */
	}

	fclose(fin);
	if (ferr != NULL) {
		fclose(ferr);
		globalfile = NULL;
	}
	/*fclose(ftest);*/

	gfree(gmtx, nnode); /* FREE GLOBAL  MATRIX. */
	/* gfree(cmtx,nnode); */ /* ERROR. *//* FREE DAMPING MATRIX. */
	gfree(mmtx, nnode); /* FREE MASS    MATRIX. */

	/* free(confs); */

	freematrix(drccos, 3);
	freematrix(estiff, 12);
	freematrix(tmatrix, 12);
	freematrix(e, 12);
	freematrix(t, 12);
	free(estress);
	free(gdisp);
	free(edisp);

	af->eigenvec = (double **)malloc(1 * sizeof(double *));
	*((af->eigenvec) + 0) = du;

	errormessage(" ");
	errormessage("COMPLETED.");
	if (fout != NULL)
		fprintf(fout, "COMPLETED.\n");
	if (fout != NULL)
		fclose(fout);

	memory2 = availablephysicalmemory("REMAIN:");
	sprintf(string, "CONSUMPTION:%ld[BYTES]", (memory0 - memory2));
	errormessage(string);
	errormessage(" ");

	return 0;
} /* gnshn101 */

int gnshn102(struct arclmframe *af)
    /* gosei:timely, t1:first */
	/* DYNAMIC RESPONSE ANALYSIS. */
	/* NODE MASSES ARE ALREADY SUMMED IN EXTRACTION. */
	/* GRAVITY LOADED ANALYSIS IS ALREADY DONE BY ARCLM. */
    {
	DWORD memory0, memory1, memory2;

	FILE *fin, *fout, *ftxt, *ftxt0, *ftxt1, *ferr; /* FILE 8 BYTES */
	double *ddisp, *dreact, emax, ratio0, ratio1;
	struct memoryelem *melem;
	char dir[] = DIRECTORY; /* DATA DIRECTORY */
	char s[80], string[4000], name1[20], name2[20];
	int i, j, ii, jj, iii, jjj;
	int n11, n22, n33, n44, n55, n66, n77, n88, n99, n100;
	int nnode, nelem, nsect, nreact,nlen;
	long int nlap, laps;
	long int loffset, msize;
	long int time;

	struct gcomponent ginit = {
		0, 0, 0.0, NULL
	};

	struct gcomponent *gmtx, *g, *p; /* GLOBAL MATRIX [K] */
	/* double *gvct; */                                 /* GLOBAL VECTOR */
	double **drccos, **tmatrix, **estiff, *estress, **e, **t;
	double *gdisp, *edisp;
	double determinant, func[2];
	clock_t t0, t1, t2;
	struct onode *nodes, *mnode;
	struct owire elem, *elems;
	struct oconf *confs = NULL;
	struct osect esect;
	double *nmass;

	FILE *fxacc, *fyacc, *fzacc;
	long int ndata[3];
	double ddt, dacc[3], tacc[3], tvel[3], tdis[3], afact;
	double pacc[3], pvel[3], pdis[3];
	double *gacc;
	double *u, *ud, *udd, *du, *dud, *dudd; /* VECTORS */
	double Wkt, Wet, Wpt, Wot, *wk[3], *wo[3]; /* ENERGY */
	struct accdata acc[3];
	struct gcomponent *mmtx; /* MASS MATRIX [M] */

	/* int neig=NEIGEN,ncount; */
	int neig = 1, ncount; /* T1 */
	/* int neig=5,ncount; */                  /* T1,T2,---,Tn */
	double **evct, *eigen; /* EIGEN VECTORS,EIGEN VALUES */
	double eps = 1.0E-16, data, d;
	double h1 = DAMPING, w1, T1, Ti, T11, ene1, ene2;
	struct gcomponent *cmtx, *mtx1, *mtx2, *mtx3, *mtx4, *mtx5, *mtx6, *mtx7, *mtx8, *mtx9, *mtx10, *mtx11, *mtx12, *mtx13, *mtx14, *mtx15, *mtx16, *mtx17, *mtx18, *mtx19, *mtx20, *mtx21, *mtx22; /* MATRIX [C] */
    struct gcomponent *mtx01, *mtx02;
	double mtotal;
	int dtype = DAMPINGTYPE;

	struct snode *sn;



    double *ncr; /*BCLNG CONDENSATION RESULT*/      //ujioka

	memory0 = availablephysicalmemory("INITIAL:"); /* MEMORY AVAILABLE */

	/* INPUT FILE */
	fin = fgetstofopen(dir, "r", ID_INPUTFILE);
	/* fin=NULL; */

	/* GET YIELD SURFACE */
	if (fin != NULL) {
		inputinit(fin, &(af->nnode), &(af->nelem), &(af->nsect));
		for (i = 0; i < (af->nsect); i++)
			readsect(fin, (af->sects + i));

		fclose(fin);
		fin = NULL;
	}

	/* OUTPUT FILE */
	fout = fgetstofopen("\0", "w", ID_OUTPUTFILE);

	/* fout=NULL; */
	globalfile = fout;

	ferr = fopen("gnshn.txt", "w"); /* ERROR FILE */
	/* ferr=NULL; */

	if (ferr != NULL) {
		fprintf(ferr, "   LAP     TIME");
		fprintf(ferr,
			"      ELEM          Nzi          Qyi          Qyj          Mxi          Mxj");
		fprintf(ferr,
			"      ELEM          Nzi          Qyi          Qyj          Mxi          Mxj");
		fprintf(ferr,
			"      ELEM          Nzi          Qyi          Qyj          Mxi          Mxj");
		fprintf(ferr,
			"      ELEM          Nzi          Qyi          Qyj          Mxi          Mxj");
		fprintf(ferr,
			"      ELEM          Nzi          Qyi          Qyj          Mxi          Mxj");
		fprintf(ferr,
			"      ELEM          Nzi          Qyi          Qyj          Mxi          Mxj");
		fprintf(ferr, "      NODE           Dx");
		fprintf(ferr, "      NODE           Dx");
		fprintf(ferr, "      NODE           Dx");
		fprintf(ferr, "\n");
	}

	ftxt = fopen("gnshn.huk", "w"); /* HUKUGENRYOKU GRAPH */
	/* ftxt=NULL; */

	if (ftxt != NULL) {
		fprintf(ftxt, "   LAP");
		for (i = 0; i < FLOORS; i++)
			fprintf(ftxt, "         Dx%2d         Dy%2d", i + 1, i + 1);
		for (i = 0; i < FLOORS - 1; i++)
			fprintf(ftxt, "         Qx%2d         Qy%2d", i + 1, i + 1);
		fprintf(ftxt, "      Accx      Accy");
		fprintf(ftxt, "        Vx        Vy");
		fprintf(ftxt, "        Dx        Dy");
		fprintf(ftxt, "          Wk       Wk+We    Wk+We+Wp");
		fprintf(ftxt, "\n");
	}

	/*edit*/
	/* ENERGY INFO by KANAZAWA */
	nlen = strlen((wdraw.childs+1)->inpfile);
	strcpy(name1, (wdraw.childs+1)->inpfile);
	name1[nlen-3] = 'e';
	name1[nlen-2] = 'e';
	name1[nlen-1] = 'n';
	name1[nlen-0] = '\0';
	ftxt0 = fopen(name1, "w");

	strcpy(name2, (wdraw.childs+1)->inpfile);
	name2[nlen-3] = 't';
	name2[nlen-2] = 'e';
	name2[nlen-1] = 'n';
	name2[nlen-0] = '\0';
	ftxt1 = fopen(name2, "w");
	if (ftxt1 != NULL) {
	  fprintf(ftxt1, "Nstep Energy Check\n");
	}

	/* ENERGY INFO by KANAZAWA */

	/* step number og getting matrix */

    /*
	n11 = 261;
	n22 = 262;
	n33 = 256;
	n44 = 257;
	n55 = 258;
	n66 = 259;
	n77 = 260;
	n88 = 261;
	n99 = 255;
	n100 = 256;
	*/

	/* step number og getting matrix */

	t0 = clock(); /* CLOCK BEGIN. */

	/* ACCCELERATION DATA INPUT */

	/* change x wave*/
	/* fxacc=fopen("elcnns.dat","r"); */ /* X */
	/* fxacc=fopen("phs00_564.txt","r"); */ /* X */
	/* fxacc=fopen("JmaKobeNS.txt","r"); */ /* X **** JMA KOBE N-S **** */
	/* fxacc=fopen("JMA_Kobe_E_NS.txt","r"); */ /* X E-defense */
	fxacc = fopen("TAT95000NS.txt", "r");
	/* X E-defense **** HANSHIN TAKATORI N-S **** */
	/* fxacc=fopen("993_Ko_Hachinohe_EW.dat","r"); */ /* Hakushima */
	/* fxacc=fopen("626_Si_Takatori_NS.dat","r"); */ /* Hakushima */
	/* fxacc=fopen("Zohuku-Hcnhew-Level2.dat","r"); */ /* Toshima */
	/* fxacc=fopen("Zohuku-Kobens-Level2.dat","r"); */ /* Toshima */
	/* fxacc=fopen("HcnhEW.txt","r"); */ /* Toshima  **** HACHINOHE E-W **** */
	/* fxacc=NULL;*/  /* X */

	/* change y wave*/
	/* fyacc=fopen("elcnew.dat","r"); */ /* Y */
	/* fyacc=fopen("tancho.dat","r"); */ /* Y */
	/* fyacc=fopen("ohuku.dat","r"); */ /* Y */
	/* fyacc=fopen("JmaKobeEW.txt","r"); */ /* Y */
	/* fyacc=fopen("JMA_Kobe_E_EW.txt","r"); */ /* Y E-defense */
	/*fyacc=fopen("TAT95000EW.txt","r");*/ /* Y E-defense */
	/* fyacc=fopen("993_Ko_Hachinohe_EW.dat","r"); */ /* Hakushima */
	/* fyacc=fopen("Zohuku-Hcnhew-Level2.dat","r"); */ /* Toshima */
	/* fyacc=fopen("Zohuku-Kobens-Level2.dat","r"); */ /* Toshima */
	/* fyacc=fopen("HcnhEW.txt","r"); */ /* Toshima */
	fyacc = NULL; /* Y */

	/* change z wave*/
	/* fzacc=fopen("94-ElcnUD.txt","r"); */ /* Z */
	/* fzacc=fopen("JmaKobeUD.txt","r"); */ /* Z */
	/* fzacc=fopen("94-TaftUD.txt","r"); */ /* Z */
	/* fzacc=fopen("kokuji_hachinohe_ud.txt","r"); */ /* Z */
	/* fzacc=fopen("kokuji_tohoku_ud.txt","r"); /*Z */
	/* fzacc=fopen("kokuji_jmakobe_ud.txt","r"); */ /* Z */
	/* fzacc=fopen("JmaKobeUD.txt","r"); */ /* Z */
	/* fzacc=fopen("JMA_Kobe_E_UD.txt","r"); */ /* Z E-defense */
	/* fzacc=fopen("TAT95000UD.txt","r"); */ /* Z E-defense */
	fzacc = NULL; /* Z */

	acc[0].a = inputacceleration(fxacc, &(acc[0].ndata), &(acc[0].dt));
	acc[1].a = inputacceleration(fyacc, &(acc[1].ndata), &(acc[1].dt));
	acc[2].a = inputacceleration(fzacc, &(acc[2].ndata), &(acc[2].dt));

	/* change time increment*/
	ddt = 0.020; /* TIME INCREMENT[sec] */
	/* ddt=0.010; */  /* TIME INCREMENT FOR HAKONE UD[sec] */
	/* ddt=0.005; */  /* TIME INCREMENT FOR E-DEFENSE[sec] */
	/* EL CENTRO:dt=0.02[sec],ddt=0.005[sec],4000[laps]=20[sec] */
	if (fout != NULL)
		fprintf(fout, "TIME INCREMENT ddt = %.5f\n", ddt);

	/* change strongth of earthquake*/
	/* äÓèÄâª */
	/* afact=0.2; */ /* KOKUJI LEVEL 1 */
	/* afact=1.0; */ /* KOKUJI LEVEL 2 */
	/* afact=1.0; */ /* LEVEL 3 */
	/* afact=25.0/ 92.065; */ /* JMA KOBE NS LEVEL 1 */
	/* afact=25.0/170.592; */ /* TAKATORI NS+EW E-defense LEVEL 1 */
	/* afact=450.0/341.0; */  /* EL CENTRO NS LEVEL 2 */
	/* afact=50.0/ 92.065; */ /* JMA KOBE NS LEVEL 2 */
	/* afact=150.0/ 92.065; */  /***** JMA KOBE NS**** */
	/* afact=50.0/106.583; */ /* JMA KOBE NS+EW LEVEL 2 */
	/* afact=50.0/103.801; */ /* JMA KOBE NS+EW E-defense LEVEL 2 */
	/* afact=25.0/149.843; */ /* TAKATORI NS LEVEL 1 */
	/* afact=50.0/149.843; */ /* TAKATORI NS LEVEL 2 */
	afact = 90 / 149.843; /***** TAKATORI NS**** */
	/* afact=50.0/170.592; */ /* TAKATORI NS+EW E-defense LEVEL 2 */
	/* afact=60.0/170.592; */ /* TAKATORI NS+EW E-defense LEVEL 3.1 */
	/* afact=75.0/170.592; */ /* TAKATORI NS+EW E-defense LEVEL 3.1 */
	/* afact=2.879; */ /* TAFT EW LEVEL 2 */
	/* afact=50.0/ 36.99; */   /* HACHINOHE EW LEVEL 2 */
	/* afact=150.0/ 36.99; */  /***** HACHINOHE EW**** */

	if (fout != NULL)
		fprintf(fout, "ACCELERATION FACTOR = %.5f\n", afact);

	ene1 = 0.0;
	ene2 = 0.0;

	afact *= 1.0; /* LOCATION FACTOR */
	/* afact*=1.25; */ /* LOCATION FACTOR */

	afact *= 0.01; /* ACCELERATION FACTOR [cm/sec2] INTO [m/sec2] */

	/* COUNT TOTAL LAPS */
	for (i = 0; i < 3; i++) {
		if (acc[i].a == NULL)
			ndata[i] = 0;
		else {
			ndata[i] = (long int)((acc[i].ndata - 1) * (acc[i].dt) / ddt) + 1;
		}
	}
	laps = 0;
	if (ndata[0] > 0)
		laps = ndata[0];
	else if (ndata[1] > 0)
		laps = ndata[1];
	else if (ndata[2] > 0)
		laps = ndata[2];
	if (ndata[1] > 0 && laps > 0 && ndata[1] < laps)
		laps = ndata[1];
	if (ndata[2] > 0 && laps > 0 && ndata[2] < laps)
		laps = ndata[2];

	if (laps > MAXLAP)
		laps = MAXLAP;
	sprintf(string, "LAPS=%d", laps);

	nnode = af->nnode;
	nelem = af->nelem;
	nsect = af->nsect;
	nreact = af->nreact;

	nodes = af->nodes;
	elems = af->elems;
	confs = af->confs;
	nmass = af->nmass;


	sprintf(string, "NODES=%d ELEMS=%d SECTS=%d", nnode, nelem, nsect);
	errormessage(string);
	if (fout != NULL)
		fprintf(fout, "%s\n", string);

	msize = 6 * nnode; /* SIZE OF GLOBAL MATRIX. */

	ncount = 0;
	for (iii = 0; iii < msize; iii++) {
		if ((confs + iii)->iconf != 1)
			ncount++;
	}
	sprintf(string, "MATRIX DIMENSION=%d", ncount);

	gmtx = (struct gcomponent*) /* DIAGONALS OF GLOBAL MATRIX. */
		malloc(msize*sizeof(struct gcomponent));

	gacc = mallocdoublevector(msize); /* ACCEL VECTOR */
	if (gmtx == NULL || gacc == NULL)
		return 0;
	for (i = 0; i < msize; i++) {
		(gmtx + i)->down = NULL; /* GLOBAL MATRIX INITIALIZATION. */
	}

#if 1 /*FOR PRORORTION TO INITIAL STIFFNESS, OCCASIONAL STIFFNESS II.*/
	if (dtype == 1 || dtype == 3) {
		cmtx = (struct gcomponent*)malloc(msize*sizeof(struct gcomponent)); /* DIAGONALS OF MATRIX [C]. */
		if (cmtx == NULL)return 0;
		for (i = 0; i < msize; i++)
		{
			(cmtx + i)->down = NULL; /* MATRIX [C] INITIALIZATION. */
		}
	}
#endif

	/* DISPLACEMENT VECTORS. */
	u = mallocdoublevector(msize);
	ud = mallocdoublevector(msize);
	udd = mallocdoublevector(msize);
	du = mallocdoublevector(msize);
	dud = mallocdoublevector(msize);
	dudd = mallocdoublevector(msize);
	if (u == NULL || ud == NULL || udd == NULL || du == NULL || dud == NULL ||
		dudd == NULL)
		return 0;
	for (i = 0; i < msize; i++) {
		*(u + i) = 0.0;
		*(ud + i) = 0.0;
		*(udd + i) = 0.0;
		*(du + i) = 0.0;
		*(dud + i) = 0.0;
		*(dudd + i) = 0.0;
	}

	mmtx = (struct gcomponent*) /* DIAGONALS OF MASS MATRIX. */
		malloc(msize*sizeof(struct gcomponent));
	if (mmtx == NULL)
		return 0;
	for (i = 0; i < msize; i++) {
		(mmtx + i)->m = (unsigned short int)(i + 1);
		(mmtx + i)->n = (unsigned short int)(i + 1);
		(mmtx + i)->value = 0.0;
		(mmtx + i)->down = NULL;
	}

	ddisp = (double *)malloc(6 * nnode*sizeof(double));
	af->ddisp = ddisp; /* DISPLACEMENT:6 DIRECTIONS. */

	melem = (struct memoryelem*) malloc(nelem*sizeof(struct memoryelem));
	af->melem = melem; /* CODE,12 BOUNDARIES,12 STRESS. */
	/* melem=af->melem; */ /* ? */

	initialform(nodes, ddisp, nnode); /* ASSEMBLAGE FORMATION. */

	mnode = (struct onode*)malloc(nnode*sizeof(struct onode));
	for (i = 0; i < nnode; i++) {
		*(mnode + i) = *(nodes + i);
	}

	initialelem(elems, melem, nelem); /* ASSEMBLAGE ELEMENTS. */
	/* for(i=0;i<nelem;i++) inputelem(elems,melem,i,(elems+i)); */ /* ? */

	dreact = (double *)malloc(nreact*sizeof(double)); /* REACTION FILE. */
	af->dreact = dreact;
	initialreact(fin, dreact, nreact); /* ASSEMBLAGE LONG REACTIONS. */

	GetAsyncKeyState(VK_LBUTTON); /* CLEAR KEY LEFT. */
	GetAsyncKeyState(VK_RBUTTON); /* CLEAR KEY RIGHT. */

	drawglobalaxis((wdraw.childs + 1)->hdcC, (wdraw.childs + 1)->vparam, 0, 0,
		255); /* DRAW GLOBAL AXIS. */

	if (wsurf.hwnd != NULL) {
		drawyieldsurface((wsurf.childs + 1)->hdcC, (wsurf.childs + 1)->vparam,
			SURFACEX, SURFACEY, SURFACEZ, NULL);
		overlayhdc(*(wsurf.childs + 1), SRCPAINT); /* UPDATE DISPLAY. */
	}
	af->fsurface = fopen("canbin.sfc", "wb+"); /* STRESS ON SURFACE. */

	/* ELEMENT ENERGY WITHOUT GRAVITY. */
	for (i = 0; i < nelem; i++) {
		(elems + i)->Ee[0] = 0.0;
		(elems + i)->Ee[1] = 0.0;
		(elems + i)->Ep[0] = 0.0;
		(elems + i)->Ep[1] = 0.0;
	}

	/* GLOBAL ENERGY. */
	wk[0] = (double *)malloc(nnode*sizeof(double));
	wk[1] = (double *)malloc(nnode*sizeof(double));
	wk[2] = (double *)malloc(nnode*sizeof(double));
	wo[0] = (double *)malloc(nnode*sizeof(double));
	wo[1] = (double *)malloc(nnode*sizeof(double));
	wo[2] = (double *)malloc(nnode*sizeof(double));
	for (ii = 0; ii < 3; ii++) {
		for (jj = 0; jj < nnode; jj++) {
			*(wk[ii] + jj) = 0.0;
			*(wo[ii] + jj) = 0.0;
		}
	}

	for (ii = 0; ii < 3; ii++) {
		tacc[ii] = 0.0;
		tvel[ii] = 0.0;
		tdis[ii] = 0.0;
	}

	/* ASSEMBLAGE MASS MATRIX. */
	mtotal = assemmass(mmtx, nnode, confs, nmass);
	if (fout != NULL)fprintf(fout, "TOTAL MASS = %12.5f\n\n", mtotal);


	/* MEMORIES. */
	drccos = mallocdoublematrix(3);
	estiff = mallocdoublematrix(12);
	tmatrix = mallocdoublematrix(12);
	e = mallocdoublematrix(12);
	t = mallocdoublematrix(12);
	estress = mallocdoublevector(12);
	gdisp = mallocdoublevector(12);
	edisp = mallocdoublevector(12);


#if 1 /*FOR PRORORTION TO INITIAL STIFFNESS. (DAMPING MATRIX (dtype=1))*/
	if (dtype == 1) {
		for (i = 1; i <= msize; i++) /* INITIALIZATION. */ {
			g = (cmtx + (i - 1))->down; /* NEXT OF DIAGONAL. */
			while (g != NULL) /* CLEAR ROW. */ {
				p = g;
				g = g->down;
				free(p);
			}
			ginit.m = (unsigned short int)i;
			/* ginit.n=(unsigned short int)i; */
			*(cmtx + (i - 1)) = ginit;
		}
		for (i = 0; i < nelem;
		i++) /* TEST */      /* ASSEMBLAGE DAMPING MATRIX. */ {
			inputelem(elems, melem, i, &elem); /* READ ELEMENT DATA. */

			for (ii = 0; ii <= 1; ii++) {
				for (jj = 0; jj <= 5; jj++) {
					(elems + i)->iconf[ii][jj] = elem.iconf[ii][jj];
				}
			}

			inputnode(ddisp, elem.node[0]); /* HEAD */
			inputnode(ddisp, elem.node[1]); /* TAIL */

			/* elem.sect=(elems+i)->sect; */
			/* READ SECTION DATA. */
			esect = *((elems + i)->sect);
			elem.sect = &esect;
			/* elem.sect->area=0.0; */

			directioncosineII(elem.node[0]->d[0], elem.node[0]->d[1],
				elem.node[0]->d[2], elem.node[1]->d[0], elem.node[1]->d[1],
				elem.node[1]->d[2], elem.cangle, drccos); /* [DRCCOS] */

			transmatrixII(drccos, tmatrix); /* TRANSFORMATION MATRIX. */
			assememtxII(elem, estiff); /* ELASTIC MATRIX OF ELEMENT. */
			estiff = modifyhinge(elem, estiff); /* MODIFY MATRIX. */
			/* estiff=assempmtx(elem,estiff); */      /* ADD PLASTIC MATRIX. */
			transformationII(estiff, tmatrix, e, t); /* [K]=[Tt][k][T] */
			assemgstiffness(cmtx, estiff, &elem); /* ASSEMBLAGE. */
		}

		/*change eigen period*/
		T1 = 1.58; /* FUNDAMENTAL NATURAL PERIOD 1éüå≈óLé¸ä˙ */
		/* T1=0.78; */ /* FOR E-DEFENSE*[sec] */
		/* T1=0.74210; */ /* FOR HAKUSHIMA Y */
		/* T1=0.19261; */ /* FOR HAKUSHIMA X */
		w1 = 2.0 * PI / T1;

		for (ii = 1; ii <= msize; ii++) {
			for (jj = 1; jj <= ii; jj++) {
				gread(cmtx, ii, jj, &data);
				data *= (2.0*h1 / w1);
				gwrite(cmtx, ii, jj, data);
			}
		}

		sprintf(string, "DAMPING MATRIX %ld COMPS ASSEMBLED.", comps);
		laptime(string, t0);
	}
#endif
	ncr=mallocdoublevector(nelem);

	definencr(&arc,&*ncr);          //ujioka

	/*change analysis*/
	bool FLAG2 = 0; /* 0:bisecsylvester  , 1:deigabgeneral */
	/* Ç±Ç±Ç©ÇÁêUìÆâêÕÇÃåJÇËï‘ÇµââéZ */
	for (nlap = 1; nlap <= laps; nlap++) { /*Å@n11Ç≈ã≠êßèIóπ or lapsÇ‹Ç≈Å@*/

		/* af->nlaps=nlap; */
		af->nlaps = 1;

		sprintf(string, "LAP:%d/%d", nlap, laps);
		errormessage(string);
		if (fout != NULL)
			fprintf(fout, "%s\n", string);
		if (ferr != NULL)
			fprintf(ferr, " %5d %8.3f", nlap, nlap*ddt);
		if (ftxt != NULL)
			fprintf(ftxt, " %5d", nlap);

		setlaps((wmenu.childs + 2)->hwnd, nlap, laps);

		memory1 = availablephysicalmemory("REMAIN:"); /* MEMORY AVAILABLE */


		for (i = 0; i < neig; i++)
		{
			*(evct + i) = (double *)malloc(msize*sizeof(double));
			for (j = 0; j < msize; j++)
				* (*(evct + i) + j) = 0.0;
		}

		for (i = 1; i <= msize; i++) /* GLOBAL MATRIX INITIALIZATION. */ {
			g = (gmtx + (i - 1))->down; /* NEXT OF DIAGONAL. */
			while (g != NULL) /* CLEAR ROW. */ {
				p = g;
				g = g->down;
				free(p);
			}

			ginit.m = (unsigned short int)i;
			/* ginit.n=(unsigned short int)i; */

			*(gmtx + (i - 1)) = ginit;
		}

		for (i = 1; i <= msize; i++) /* INITIALIZATION. */
		{
			g = (cmtx + (i - 1))->down; /* NEXT OF DIAGONAL. */
			while (g != NULL) /* CLEAR ROW. */ {
				p = g;
				g = g->down;
				free(p);
			}
			ginit.m = (unsigned short int)i;
			/* ginit.n=(unsigned short int)i; */
			*(cmtx + (i - 1)) = ginit;
		}

		comps = msize; /* INITIAL COMPONENTS=DIAGONALS. */

		laptime("ASSEMBLING GLOBAL MATRIX.", t0);

		clearwindow(*(wdraw.childs + 1));
		/* drawarclmtitle((wdraw.childs+1)->hdcC,
		 (wdraw.childs+1)->vparam,*af,0,ONSCREEN); */
		for (i = 0; i < (wdraw.childs + 1)->org.ntext; i++) {
			sn = (wdraw.childs + 1)->org.texts + i;
			wdraw.strset = addtext(wdraw.strset, &(wdraw.nstring), sn->str,
				sn->n.d[0], sn->n.d[1]);
		}
		SetTextColor((wdraw.childs + 1)->hdcC, RGB(255, 255, 255));
		drawtexts((wdraw.childs + 1)->hdcC, wdraw.strset, wdraw.nstring,
			(wdraw.childs + 1)->vparam);

		for (i = 0; i < nelem; i++) /* ASSEMBLAGE GLOBAL MATRIX. */ {
			inputelem(elems, melem, i, &elem); /* READ ELEMENT DATA. */

			for (ii = 0; ii <= 1; ii++) {
				for (jj = 0; jj <= 5; jj++) {
					(elems + i)->iconf[ii][jj] = elem.iconf[ii][jj];
				}
			}

			inputnode(ddisp, elem.node[0]); /* HEAD */
			inputnode(ddisp, elem.node[1]); /* TAIL */

			elem.sect = (elems + i)->sect; /* READ SECTION DATA. */

			if ((wdraw.childs + 1)->hdcC != NULL) /* DRAW DEFORMED ELEMENT. */ {
				drawglobalwire((wdraw.childs + 1)->hdcC,
					(wdraw.childs + 1)->vparam, *af, elem, 255, 255, 255, 255,
					255, 255, 0, ONSCREEN);
			}

			directioncosineII(elem.node[0]->d[0], elem.node[0]->d[1],
				elem.node[0]->d[2], elem.node[1]->d[0], elem.node[1]->d[1],
				elem.node[1]->d[2], elem.cangle, drccos); /* [DRCCOS] */

			transmatrixII(drccos, tmatrix); /* TRANSFORMATION MATRIX. */


			assememtxII(elem, estiff); /* ELASTIC MATRIX OF ELEMENT. */


			estiff = modifyhinge(elem, estiff); /* MODIFY MATRIX. */


			estiff = assempmtx(elem, estiff); /* ADD PLASTIC MATRIX. */


			transformationII(estiff, tmatrix, e, t); /* [K]=[Tt][k][T] */


			assemgstiffness(gmtx, estiff, &elem); /* ASSEMBLAGE. */

			assemgstiffness(cmtx, estiff, &elem); /* ASSEMBLAGE. */

		}

		/*

		for (i = 1; i <= msize; i++)
		{
			g = (mtx01 + (i - 1))->down;
			while (g != NULL)
			{
				p = g;
				g = g->down;
				free(p);
			}

			ginit.m = (unsigned short int)i;
			*(mtx01 + (i - 1)) = ginit;
		}

        for (i = 0; i < msize; i++) {
		(mtx02 + i)->m = (unsigned short int)(i + 1);
		(mtx02 + i)->n = (unsigned short int)(i + 1);
		(mtx02 + i)->value = 0.0;
		(mtx02 + i)->down = NULL;
		}

		mtx01 = copygcompmatrix(gmtx, msize);
		mtx02 = copygcompmatrix(mmtx, msize);

		deigabgeneral(mtx01, mtx02, confs, msize, -neig, neig, eps,
						eigen, evct);

		T1 = 2.0 * PI / sqrt(*(eigen + 0));
		w1 = 2.0 * PI / T1;

		*/

		for (ii = 1; ii <= msize; ii++)
		{
			for (jj = 1; jj <= ii; jj++)
			{
				gread(cmtx, ii, jj, &data);
				data *= (2.0*h1 / w1);
				gwrite(cmtx, ii, jj, data);
			}
		}





		sprintf(string, "GLOBAL MATRIX %ld COMPS ASSEMBLED.", comps);
		laptime(string, t0);

		overlayhdc(*(wdraw.childs + 1), SRCPAINT); /* UPDATE DISPLAY. */

			for (i = 0; i < 3; i++) {
				if (acc[i].a == NULL)
					dacc[i] = 0.0;
				else {
					dacc[i] = accelerationincrement(acc[i], ddt, nlap);
					dacc[i] *= afact;
				}

				pacc[i] = tacc[i];
				pvel[i] = tvel[i];
				pdis[i] = tdis[i];

				tacc[i] += dacc[i]; /* GROUND ACCELERATION. */
				tvel[i] += ddt * (pacc[i] + tacc[i]) / 2.0;
				/* TOTAL GROUND VELOCITY. */
				/* tdis[i]+=ddt*(pdis[i]+tvel[i])/2.0; */ /* ERROR */
				/* GROUND DISPLACEMENT. */
				tdis[i] += ddt * (pvel[i] + tvel[i]) / 2.0;
				/* GROUND DISPLACEMENT. */
			}

		if (fout != NULL)
			fprintf(fout, "ACCELERATION=%.8f,%.8f,%.8f\n", dacc[0], dacc[1],
			dacc[2]);

		assemaccel(gacc, dacc, nnode, af); /* ACCELERATION VECTOR. */



		if (fout != NULL)
			fprintf(fout, "T1=%.5f w1=%.5f h1=%.5f\n", T1, w1, h1);

		/* NEWMARK'S BETA PROCESS. */
		laptime("NEWMARK'S BETA.", t0);
		determinant = newmarkbeta(gmtx, cmtx, mmtx, gacc, confs, nnode, u, ud,
			udd, du, dud, dudd, ddt, (1.0 / 4.0));

		sprintf(string, "DETERMINANT=%.5E COMPS=%ld", determinant, comps);
		errormessage(string);
		if (fout != NULL)
			fprintf(fout, "%s\n", string);


		if (determinant <= 0.0) {
			errormessage(" ");
			errormessage("INSTABLE TERMINATION.");
			if (fout != NULL)
				fprintf(fout, "INSTABLE TERMINATION.\n");

			laptime("\0", t0);

			fclose(fin);

			gfree(gmtx, nnode); /* FREE GLOBAL MATRIX. */
			gfree(cmtx, nnode);
			/*
			gfree(mtx01, nnode);
			gfree(mtx02, nnode);
			*/
			free(gacc);
			/* free(confs); */
			freematrix(drccos, 3);
			freematrix(estiff, 12);
			freematrix(tmatrix, 12);
			freematrix(e, 12);
			freematrix(t, 12);
			free(estress);
			free(gdisp);
			free(edisp);
			free(evct);
			free(eigen);

			memory2 = availablephysicalmemory("REMAIN:");
			sprintf(string, "CONSUMPTION:%ld[BYTES]", (memory0 - memory2));
			errormessage(string);

			return 1;
		}

		laptime("OUTPUT INTO FILE.", t0);

		if (fout != NULL)
			fprintf(fout, "\"DISPLACEMENT\"\n");
		outputdisp(du, fout, nnode, nodes); /* INCREMENTAL DISPLACEMENT. */
		/* while(!GetAsyncKeyState(VK_LBUTTON))
		 ; */                                   /* LEFT CLICK TO CONTINUE. */

		/* if(ferr!=NULL) fprintf(ferr," %7.5f",((double)nlap)*ddt); */

		if (fout != NULL)
			fprintf(fout, "\"STRESS\"\n");
		for (i = 0; i < nelem; i++) /* STRESS OUTPUT,UPDATE. */ {
			inputelem(elems, melem, i, &elem);

			inputnode(ddisp, elem.node[0]);
			inputnode(ddisp, elem.node[1]);

			elem.sect = (elems + i)->sect; /* READ SECTION DATA. */

			elemstressII(estress, &elem, du, melem, fout, drccos, tmatrix,
				estiff, gdisp, edisp, func, ftxt);

			(elems + i)->Ee[0] = elem.Ee[0];
			(elems + i)->Ee[1] = elem.Ee[1];
			(elems + i)->Ep[0] = elem.Ep[0];
			(elems + i)->Ep[1] = elem.Ep[1];
		}

		/* UPDATE DISPLAY. */
		if (wsurf.hwnd != NULL)
		{
			drawyieldsurface((wsurf.childs + 1)->hdcC,
				(wsurf.childs + 1)->vparam, SURFACEX, SURFACEY, SURFACEZ,
				af->fsurface);
			overlayhdc(*(wsurf.childs + 1), SRCPAINT);

		}




		if (fout != NULL)
			fprintf(fout, "\"REACTION\"\n");
		outputreaction(gmtx, du, nodes, confs, dreact, fout, nnode);

		updateform(ddisp, du, nnode); /* FORMATION UPDATE. */
		if (fout != NULL)
			fprintf(fout, "\"CURRENT FORM\"\n\n");


		fprintf(ferr, "\n");

		/* CALCULATE FLOOR DISPLACEMENT,FLOOR HUKUGENRYOKU. */

		/* floorvalues(fout,ftxt,af,ddisp,dacc); */
		floorvalues(NULL, ftxt, af, ddisp, dacc);
		/* if(ftxt!=NULL) fprintf(ftxt,"\n"); */

		/* SUM GLOBAL ENERGY. */
		Wkt = 0.0;
		Wet = 0.0;
		Wpt = 0.0;
		Wot = 0.0;
		for (ii = 0; ii < nnode; ii++) {
			for (jj = 0; jj < 3; jj++) {
				*(wk[jj] + ii) = 0.5 * (*(nmass + ii)) * (*(ud + 6 * ii + jj)) *
					(*(ud + 6 * ii + jj));
				Wkt += *(wk[jj] + ii);

				*(wo[jj] + ii) -=
					0.5 * (*(nmass + ii)) * (2 * tacc[jj] - dacc[jj]) *
					(*(du + 6 * ii + jj));
				Wot += *(wo[jj] + ii);
			}
		}

		for (ii = 0; ii < nelem; ii++) {
			Wet += (elems + ii)->Ee[0];
			Wet += (elems + ii)->Ee[1];
			Wpt += (elems + ii)->Ep[0];
			Wpt += (elems + ii)->Ep[1];
		}

		if (ftxt != NULL) {
			fprintf(ftxt, " %9.3f %9.3f", tacc[0], tacc[1]); /* Acc */
			fprintf(ftxt, " %9.3f %9.3f", tvel[0], tvel[1]); /* V */
			fprintf(ftxt, " %9.3f %9.3f", tdis[0], tdis[1]); /* D */
			fprintf(ftxt, " %11.7f %11.7f %11.7f", Wkt, Wkt + Wet,
				Wkt + Wet + Wpt); /* W */
			fprintf(ftxt, "\n");
		}

		if (ftxt1 != NULL) {
			fprintf(ftxt1, "%11.7f\n", Wpt); /* W */
		}

		t1 = laptime("\0", t0);

		memory2 = availablephysicalmemory(NULL);
		sprintf(string, "CONSUMPTION:%ld[BYTES]", (memory1 - memory2));
		errormessage(string);

		errormessage("L:CONTINUE R:ABORT"); /* L=LEFT R=RIGHT */
		while (!GetAsyncKeyState(VK_LBUTTON)) /* LEFT CLICK TO CONTINUE. */ {
			if (GetAsyncKeyState(VK_RBUTTON)) /* RIGHT CLICK TO ABORT. */ {
				fclose(fin);

				gfree(gmtx, nnode); /* FREE GLOBAL MATRIX. */
				gfree(mmtx, nnode); /* FREE MASS   MATRIX. */
				gfree(cmtx, nnode);
				/*
				gfree(mtx01, nnode);
				gfree(mtx02, nnode);
				*/
				free(gacc);
				/* free(confs); */
				freematrix(drccos, 3);
				freematrix(estiff, 12);
				freematrix(tmatrix, 12);
				freematrix(e, 12);
				freematrix(t, 12);
				free(estress);
				free(gdisp);
				free(edisp);
				free(evct);
				free(eigen);

				errormessage(" ");
				errormessage("ABORTED.");
				if (fout != NULL)
					fprintf(fout, "ABORTED.\n");

				laptime("\0", t0);
				return 1;
			}
			t2 = clock();
			time = (t2 - t1) / CLK_TCK;
			if (time >= WAIT)
				break; /* CONTINUE AFTER WAITING. */
		}

	} /* END OF LAP.REPEAT UNTIL INSTABLE. êUìÆåJÇËï‘ÇµÇ±Ç±Ç‹Ç≈ */

	/*edit*/
		if (ftxt0 != NULL) {
			emax = 0.0;
			for (ii = 0; ii < nelem; ii++) {
				if (emax < ((elems + ii)->Ee[0] + (elems + ii)->Ep[0])) {
					emax = ((elems + ii)->Ee[0] + (elems + ii)->Ep[0]);
				}
				if (emax < ((elems + ii)->Ee[1] + (elems + ii)->Ep[1])) {
					emax = ((elems + ii)->Ee[1] + (elems + ii)->Ep[1]);
				}
			}
			fprintf(ftxt0, "//Check for Energy Distribution//\n\n");
			fprintf(ftxt0, "Total = %-13.6f\n",Wpt);
			fprintf(ftxt0, "Max   = %-13.6f\n\n",emax);
			/*MessageBox(NULL, "BISECSYLVESTER START.", "Gnshn", MB_OK);*/
			/*deigabgeneral(mtx3, mtx4, confs, msize, -neig, neig, eps,
						eigen, evct);
			T11 = 2.0 * PI / sqrt(*(eigen + 0));
			bisecsylvester(mtx4, mtx3, confs, msize, neig, neig, eps,
						eigen, evct);
			T11 = 2.0 * PI * sqrt(*(eigen + 0));*/
			/*
			deigabgeneral(mtx3, mtx4, confs, msize, -neig, neig, eps,
						eigen, evct);
			T11 = 2.0 * PI / sqrt(*(eigen + 0));
			fprintf(ftxt0, "LAP %3d  :  T1  = %-13.6f\n\n",n11,T11);

			deigabgeneral(mtx5, mtx6, confs, msize, -neig, neig, eps,
						eigen, evct);
			T11 = 2.0 * PI / sqrt(*(eigen + 0));
			fprintf(ftxt0, "LAP %3d  :  T1  = %-13.6f\n\n",n22,T11);

			deigabgeneral(mtx7, mtx8, confs, msize, -neig, neig, eps,
						eigen, evct);
			T11 = 2.0 * PI / sqrt(*(eigen + 0));
			fprintf(ftxt0, "LAP %3d  :  T1  = %-13.6f\n\n",n33,T11);

			deigabgeneral(mtx9, mtx10, confs, msize, -neig, neig, eps,
						eigen, evct);
			T11 = 2.0 * PI / sqrt(*(eigen + 0));
			fprintf(ftxt0, "LAP %3d  :  T1  = %-13.6f\n\n",n44,T11);

			deigabgeneral(mtx11, mtx12, confs, msize, -neig, neig, eps,
						eigen, evct);
			T11 = 2.0 * PI / sqrt(*(eigen + 0));
			fprintf(ftxt0, "LAP %3d  :  T1  = %-13.6f\n\n",n55,T11);

			deigabgeneral(mtx13, mtx14, confs, msize, -neig, neig, eps,
						eigen, evct);
			T11 = 2.0 * PI / sqrt(*(eigen + 0));
			fprintf(ftxt0, "LAP %3d  :  T1  = %-13.6f\n\n",n66,T11);

			deigabgeneral(mtx15, mtx16, confs, msize, -neig, neig, eps,
						eigen, evct);
			T11 = 2.0 * PI / sqrt(*(eigen + 0));
			fprintf(ftxt0, "LAP %3d  :  T1  = %-13.6f\n\n",n77,T11);

			deigabgeneral(mtx17, mtx18, confs, msize, -neig, neig, eps,
						eigen, evct);
			T11 = 2.0 * PI / sqrt(*(eigen + 0));
			fprintf(ftxt0, "LAP %3d  :  T1  = %-13.6f\n\n",n88,T11);

			deigabgeneral(mtx19, mtx20, confs, msize, -neig, neig, eps,
						eigen, evct);
			T11 = 2.0 * PI / sqrt(*(eigen + 0));
			fprintf(ftxt0, "LAP %3d  :  T1  = %-13.6f\n\n",n99,T11);

			deigabgeneral(mtx21, mtx22, confs, msize, -neig, neig, eps,
						eigen, evct);
			T11 = 2.0 * PI / sqrt(*(eigen + 0));
			fprintf(ftxt0, "LAP %3d  :  T1  = %-13.6f\n\n",n100,T11);
			*/
			MessageBox(NULL, "BISECSYLVESTER END.", "Gnshn", MB_OK);
			for (ii = 0; ii < nelem; ii++) {
				ratio0 = ((elems + ii)->Ee[0] + (elems + ii)->Ep[0]) / (Wet+Wpt);
				ratio1 = ((elems + ii)->Ee[1] + (elems + ii)->Ep[1]) / (Wet+Wpt);
				fprintf(ftxt0, "ELEM %d", (elems + ii)->code);
				fprintf(ftxt0, " SPT %d", (elems + ii)->node[0]->code);
				fprintf(ftxt0, " ENE %-8.6f", (elems + ii)->Ee[0] + (elems + ii)->Ep[0]);
				fprintf(ftxt0, " RATIO %-8.6f", ratio0);
				fprintf(ftxt0, " EPT %d", (elems + ii)->node[1]->code);
				fprintf(ftxt0, " ENE %-8.6f", (elems + ii)->Ee[1] + (elems + ii)->Ep[1]);
				fprintf(ftxt0, " RATIO %-8.6f", ratio1);
				fprintf(ftxt0, "\n");

			}
			fclose(ftxt0);
		}

	if ((wdraw.childs + 1)->hdcC != NULL && melem != NULL && ddisp != NULL)
		/* DRAW LAST FRAME. */ {
		for (i = 1; i <= nelem; i++) {
			inputelem(elems, melem, i - 1, &elem);
			for (ii = 0; ii <= 1; ii++) /* COPY HINGE DATA. */ {
				for (jj = 0; jj <= 5; jj++) {
					(elems + i - 1)->iconf[ii][jj] = elem.iconf[ii][jj];
				}
			}

			inputnode(ddisp, elem.node[0]);
			inputnode(ddisp, elem.node[1]);

			drawglobalwire((wdraw.childs + 1)->hdcC, (wdraw.childs + 1)->vparam,
				*af, elem, 255, 255, 255, 255, 255, 255, 0, ONSCREEN);
		}
		overlayhdc(*(wdraw.childs + 1), SRCPAINT); /* UPDATE DISPLAY. */
	}

	fclose(fin);
	if (ferr != NULL) {
		fclose(ferr);
		globalfile = NULL;
	}

	gfree(gmtx, nnode); /* FREE GLOBAL  MATRIX. */
	/* gfree(cmtx,nnode); */ /* ERROR. *//* FREE DAMPING MATRIX. */
	gfree(mmtx, nnode); /* FREE MASS    MATRIX. */

	/* free(confs); */

	freematrix(drccos, 3);
	freematrix(estiff, 12);
	freematrix(tmatrix, 12);
	freematrix(e, 12);
	freematrix(t, 12);
	free(estress);
	free(gdisp);
	free(edisp);

	af->eigenvec = (double **)malloc(1 * sizeof(double *));
	*((af->eigenvec) + 0) = du;

	errormessage(" ");
	errormessage("COMPLETED.");
	if (fout != NULL)
		fprintf(fout, "COMPLETED.\n");
	if (fout != NULL)
		fclose(fout);

	memory2 = availablephysicalmemory("REMAIN:");
	sprintf(string, "CONSUMPTION:%ld[BYTES]", (memory0 - memory2));
	errormessage(string);
	errormessage(" ");

	return 0;
} /* gnshn102 */




double *inputacceleration(FILE *fin, long int *ndata, double *dt)
	/* INPUT ACCELERATION FROM TEXTFILE. */ {
	/* CAPTION */
	/* NUMBER OF DATA */
	/* INCREMENT OF TIME [s] */
	/* ACCELERATION DATA [cm/s2] */

	char **data, str[256];
	int i, nstr;
	long int idata;
	double *a;

	*ndata = 0;
	*dt = 0.0;

	if (fin == NULL)
		return NULL;

	fseek(fin, 0L, SEEK_SET);

	fgets(str, 256, fin); /* CAPTION. */

	data = fgetsbrk(fin, &nstr);
	if (nstr < 2)
		return NULL;
	if (!strcmp(*(data + 0), "DATA")) {
		*ndata = strtol(*(data + 1), NULL, 10); /* NUMBER OF DATA. */
	}

	data = fgetsbrk(fin, &nstr);
	if (nstr < 2)
		return NULL;
	if (!strcmp(*(data + 0), "DT")) {
		*dt = strtod(*(data + 1), NULL); /* INCREMENT OF TIME. */
	}

	a = (double *)malloc((*ndata)*sizeof(double));

	idata = 0;
	while (1) {
		data = fgetsbrk(fin, &nstr);
		if (nstr == 0)
			return a;

		for (i = 0; i < nstr; i++) {
			*(a + idata) = strtod(*(data + i), NULL);

			idata++;
			if (idata >= (*ndata))
				return a;
		}

		freestr(data, nstr);
	}

} /* inputacceleration */

double assemmass(struct gcomponent *mmtx, long int nnode, struct oconf *confs,double *nmass)/* ASSEMBLAGE MASS MATRIX. */
{
	long int i, j, ii;
	double gdata, nm, ge, mtotal;
	double eps = 1.0E-8;
//	double eps = 0.0;

	mtotal = 0.0;
	for (i = 0; i < nnode; i++)
	{
		nm = *(nmass + i);
		mtotal += nm;
		for(j=0;j<3; j++) /* FOR X,Y,Z LINE. */
		{
			ii=6*i+j;
			if ((confs+ii)->iconf==0)
			{
				if (nm != 0.0)
					gwrite(mmtx, ii+1, ii+1, nm);
			}
		}
		for(j=3;j<6;j++) /* FOR TX,TY,TZ LINE. */
		{
			ii=6*i+j;

			if ((confs+ii)->iconf == 0)
			{
				if (nm != 0.0)
				{
				  gdata = nm;
				}
				else
				{
				  gdata = eps;
				}
				if (gdata != 0.0)
				{
				  gwrite(mmtx, ii+1, ii+1, gdata);
				}

			}
		}
	}
	return mtotal;
} /* assemmass */

double accelerationincrement(struct accdata acc, double ddt, long int lap)
	/* INCREMENT OF ACCELERATION. */ {
	long int n;
	double t, a1, a2, da1, da2;
	double dacc;

	t = ddt * (double)(lap - 1);

	n = (long int)(t / (acc.dt));
	a1 = *(acc.a + n);
	a2 = *(acc.a + n + 1);
	da1 = a1 + (a2 - a1) * (t - ((double)n) * (acc.dt)) / (acc.dt);

	t += ddt;
	n = (long int)(t / (acc.dt));
	a1 = *(acc.a + n);
	a2 = *(acc.a + n + 1);
	da2 = a1 + (a2 - a1) * (t - ((double)n) * (acc.dt)) / (acc.dt);

	dacc = da2 - da1;

	return dacc;
} /* accelerationincrement */

void assemaccel(double *gvct, double dacc[], long int nnode, struct arclmframe *af)
	/* ASSEMBLAGE ACCELERATION INTO GLOBAL VECTOR. */ {
	long int i, j, ii;

	for (i = 0; i < (6 * nnode); i++) {
		*(gvct + i) = 0.0; /* INITIALIZATION. */
	}

	for (i = 0; i < nnode; i++) {
		for (j = 0; j < 3; j++) /* FOR X,Y,Z LINE. */ {
			ii = 6 * i + j;

			if ((af->confs + ii)->iconf == 0) {
				*(gvct + ii) = dacc[j];
			}
		}
	}
	return;
} /* assemaccel */

double newmarkbeta(struct gcomponent *gmtx, struct gcomponent *cmtx,
	struct gcomponent *mmtx, double *gacc, struct oconf *confs, int nnode,
	double *u, double *ud, double *udd, double *du, double *dud, double *dudd,
	double ddt, double beta)
	/* NEWMARK'S BETA BY INCREMENT. */
{
	long int i, j, msize;
	int nline;
	double gdata, cdata, mdata, det, sign;

	msize = 6 * nnode;

	for (i = 0; i < msize; i++)
	{
		if ((confs + i)->iconf == 0)
		{
			*(du + i) = 0.0;

			/* for(j=0;j<=i;j++) */
			for (j = 0; j < msize; j++)
			{
				if ((confs + j)->iconf == 0)
				{
					gread(gmtx, (i + 1), (j + 1), &gdata);
					if (cmtx != NULL)
					{
						gread(cmtx, (i + 1), (j + 1), &cdata);
					}
					else
					{
						cdata = 0.0;
					}
					gread(mmtx, (i + 1), (j + 1), &mdata);

					*(du + i) += -mdata * (*(gacc + j)) + mdata * ((*(ud + j)) / beta / ddt +(*(udd + j)) / 2.0 / beta) +
						cdata * ((*(ud + j)) / 2.0 / beta + (*(udd + j)) * ddt *
						(1.0 / 4.0 / beta - 1.0)); /* {dP'}:(3.59) */

					if (j <= i)
					{
						gdata += mdata / beta / (ddt * ddt)
							   + cdata / 2.0 / beta / ddt; /* [K']:(3.59) */

						gwrite(gmtx, (i + 1), (j + 1), gdata);
					}
				}
			}
		}
	}
	/* CROUT LU DECOMPOSITION. [K']{dU}={dP'} */
	nline=croutludecomposition(gmtx, du, confs, msize, &det, &sign); /* (3.58) */
	if (sign < 0.0)
	{
		return sign;
	}

	for (i = 0; i < msize; i++)
	{
		if ((confs + i)->iconf == 0)
		{

			*(dud + i) = (*(du + i)) / 2.0 / beta / ddt -
				(*(ud + i)) / 2.0 / beta - (1.0 / 4.0 / beta - 1.0) *
				(*(udd + i)) * ddt; /* (3.61) */

			(*(dudd + i)) = (*(du + i)) / beta / (ddt * ddt) -
				(*(ud + i)) / beta / ddt - (*(udd + i)) / 2.0 / beta;
			/* (3.62) */

			*(u + i) += *(du + i);
			*(ud + i) += *(dud + i);
			*(udd + i) += *(dudd + i);
		}
	}

	return sign;
} /* newmarkbeta */

void floorvalues(FILE *fout, FILE *ftxt, struct arclmframe *af, double *ddisp,
	double dacc[]) {
	/* char str[256]; */
	int i, ii, jj, d;
	long int nnode, nelem, nsect;
	long int loff;
	double disp;
	double bunbo;
	double intersect;

	double flevel[FLOORS + 1] = {-1.0, 0.1, 100};

	int fnodes[FLOORS][2];
	double meanlevel[FLOORS][2], meandisp[FLOORS][2], maxdisp[FLOORS][2];
	double r[FLOORS - 1][2], R[FLOORS - 1][2]; /* FOR RIGID RATE. */
	double rmean[2]; /* MEAN OF "r". */

	double Gx[FLOORS - 1], Gy[FLOORS - 1]; /* CENTER OF SHEAR LOAD. */
	double summass[FLOORS], sumload[FLOORS - 1][2], summoment[FLOORS - 1][2];

	double Lx[FLOORS - 1], Ly[FLOORS - 1]; /* CENTER OF RIGID. */
	int felems[FLOORS - 1];
	double *estress, *gstress;
	double K, Ex, Ey; /* RIGID,COORDINATION OF ELEMENT. */
	double sumshear[FLOORS - 1][2];
	double sumrigid[FLOORS - 1][2], sumrigidmoment[FLOORS - 1][2];
	double ex[FLOORS - 1], ey[FLOORS - 1]; /* DISTANCE FROM G TO L. */
	double Kr[FLOORS - 1];
	double rex[FLOORS - 1], rey[FLOORS - 1], Rex[FLOORS - 1], Rey[FLOORS - 1];

	double *srate[FLOORS - 1][2]; /* SHEAR  OF EACH SECTIONS. */
	double trate[FLOORS - 1][2][MAXTYPE]; /* SHEAR OF EACH SECTION TYPES. */

	struct onode *node;
	struct owire *elem;

	long int nrigid;
	struct structrigid rigid, *rigids = NULL;

	nnode = af->nnode;
	nelem = af->nelem;
	nsect = af->nsect;

	/* RATE OF FLOOR ANGLE,CENTER OF HORIZONTAL LOAD.................. */
	for (ii = 0; ii < FLOORS; ii++) /* INITIALIZATION */ {
		summass[ii] = 0.0;

		for (d = 0; d < 2; d++) {
			fnodes[ii][d] = 0;
			meanlevel[ii][d] = 0.0;
			meandisp[ii][d] = 0.0;
			maxdisp[ii][d] = 0.0;

			if (ii < (FLOORS - 1)) {
				sumload[ii][d] = 0.0;
				summoment[ii][d] = 0.0;
			}
		}
	}

	/* LEVEL,DISPLACEMENT OF NODE */
	for (i = 1; i <= nnode; i++) {
		node = af->nodes + i - 1;

		for (d = 1; d <= 2; d++) /* DIRECTION X,Y */ {
			loff = 6 * (i - 1) + d - 1; /* OFFSET */

			for (ii = 1; ii <= FLOORS; ii++) {
				if ((flevel[ii - 1] <= (node->d[GZ])) && ((node->d[GZ]) <
					flevel[ii])) {
					fnodes[ii - 1][d - 1]++; /* COUNT FLOOR NODES */
					meanlevel[ii - 1][d - 1] += (node->d[GZ]);
					if (d == 1 /* && (af->confs+6*(i-1)+2)->iconf==0 */) {
						summass[ii - 1] += (*(af->nmass + i - 1));
					}

					/* disp=*(ddisp+loff)-(node->d[d-1]); */
					disp = *(ddisp + loff) - ((af->ninit + i - 1)->d[d - 1]);
					meandisp[ii - 1][d - 1] += disp;

					if (maxdisp[ii - 1][d - 1] < fabs(disp)) {
						maxdisp[ii - 1][d - 1] = fabs(disp);
					}
				}
				if (ii >= 2 && flevel[ii - 1] <= (node->d[GZ])) {
					if ((af->confs + loff)->iconf == 0) {
						sumload[ii - 2][d - 1] +=
							dacc[d - 1] * (*(af->nmass + i - 1));
						if (d == 1) {
							summoment[ii - 2][d - 1] +=
								dacc[d - 1] * (*(af->nmass + i - 1))
								* node->d[GY];
						}
						if (d == 2) {
							summoment[ii - 2][d - 1] +=
								dacc[d - 1] * (*(af->nmass + i - 1))
								* node->d[GX];
						}
					}
				}
			}
		}
	}

	/* "DEFORMATION ANGLE" */
	for (d = 0; d < 2; d++) {
		rmean[d] = 0.0;
		for (ii = 0; ii < FLOORS; ii++) {
			if (fnodes[ii][d] > 0) {
				meanlevel[ii][d] /= fnodes[ii][d];
				meandisp[ii][d] /= fnodes[ii][d];
			}
			if (ii > 0 && meandisp[ii][d] != meandisp[ii - 1][d]) {
				r[ii - 1][d] = (meanlevel[ii][d] - meanlevel[ii - 1][d]) /
					(meandisp[ii][d] - meandisp[ii - 1][d]);
				rmean[d] += r[ii - 1][d];
			}
		}
		if (FLOORS > 1)
			rmean[d] /= (double)(FLOORS - 1);
	}

	if (fout != NULL)
		fprintf(fout, "\n");
	for (ii = 0; ii < FLOORS; ii++) {
		if (fout != NULL) {
			fprintf(fout, "%2dFL:MASS=%8.5f\n", ii + 1, summass[ii]);
		}
	}

	if (fout != NULL)
		fprintf(fout, "\n");
	for (ii = 0; ii < FLOORS; ii++) {
		if (fout != NULL) {
			fprintf(fout, "%2dFL:LEVEL=%8.5f[m] Dx=%8.5f[m] Dy=%8.5f[m]\n",
				ii + 1, meanlevel[ii][0], meandisp[ii][0], meandisp[ii][1]);
		}
		if (ftxt != NULL) {
			/* fprintf(ftxt," dx%d %8.5E dy%d %8.5E",
			 ii+1,meandisp[ii][0],
			 ii+1,meandisp[ii][1]); */
			fprintf(ftxt, " %12.5E %12.5E", meandisp[ii][0], meandisp[ii][1]);
		}
	}

	if (fout != NULL)
		fprintf(fout, "\n");
	for (ii = 0; ii < FLOORS; ii++) {
		if (fout != NULL) {
			fprintf(fout, "%2dFL:Dxmax=%8.5f[m] Dymax=%8.5f[m]\n", ii + 1,
				maxdisp[ii][0], maxdisp[ii][1]);
		}
	}

	if (fout != NULL)
		fprintf(fout, "\n");
	for (ii = 0; ii < (FLOORS - 1); ii++) {
		if (rmean[0] != 0.0)
			R[ii][0] = r[ii][0] / rmean[0];
		else
			R[ii][0] = 1.0;
		if (rmean[1] != 0.0)
			R[ii][1] = r[ii][1] / rmean[1];
		else
			R[ii][1] = 1.0;

		if (fout != NULL) {
			fprintf(fout, "%2dF:Dx/Hi=1/%4.0f Dy/Hi=1/%4.0f", ii + 1, r[ii][0],
				r[ii][1]);
			fprintf(fout, " Rx=%7.5f Ry=%7.5f\n", R[ii][0], R[ii][1]);
		}
	}

	/* CENTER OF HORIZONTAL LOAD...................................... */
	if (fout != NULL)
		fprintf(fout, "\n");
	if (fout != NULL)
		fprintf(fout, "CENTER OF LOAD.\n");
	for (ii = 1; ii <= (FLOORS - 1); ii++) {
		if (sumload[ii - 1][1] != 0.0) {
			Gx[ii - 1] = summoment[ii - 1][1] / sumload[ii - 1][1];
		}
		else
			Gx[ii - 1] = 0.0;

		if (sumload[ii - 1][0] != 0.0) {
			Gy[ii - 1] = summoment[ii - 1][0] / sumload[ii - 1][0];
		}
		else
			Gy[ii - 1] = 0.0;

		if (fout != NULL) {
			fprintf(fout, "%2dFL:Qx=%10.5f[tf] Qy=%10.5f[tf]", (ii + 1),
				sumload[ii - 1][0], sumload[ii - 1][1]);
			fprintf(fout, " Gx=%8.5f[m] Gy=%8.5f[m]\n", Gx[ii - 1], Gy[ii - 1]);
		}
	}

	/* CENTER OF RIGID................................................ */
	for (d = 0; d < 2; d++) /* INITIALIZATION */ {
		for (ii = 0; ii < (FLOORS - 1); ii++) {
			sumshear[ii][d] = 0.0;
			sumrigid[ii][d] = 0.0;
			sumrigidmoment[ii][d] = 0.0;

			srate[ii][d] = (double *)malloc(nsect*sizeof(double));
			for (jj = 0; jj < nsect; jj++)
				* (srate[ii][d] + jj) = 0.0;
			for (jj = 0; jj < MAXTYPE; jj++)
				trate[ii][d][jj] = 0.0;
		}
	}
	for (ii = 0; ii < (FLOORS - 1); ii++)
		felems[ii] = 0;

	estress = (double *)malloc(12 * sizeof(double));
	nrigid = 0;
	for (i = 1; i <= nelem; i++) {
		elem = af->elems + i - 1;

		/* GET STRESS */
		for (ii = 0; ii < 2; ii++) {
			for (jj = 0; jj < 6; jj++) {
				*(estress + 6 * ii + jj) = (af->melem + i - 1)->stress[ii][jj];
			}
		}

		for (ii = 1; ii <= (FLOORS - 1); ii++) {
			intersect = (flevel[ii] - (elem->node[0])->d[GZ]) *
				(flevel[ii] - (elem->node[1])->d[GZ]);

			/* IF ACROSSING JUDGING LEVEL */
			if (intersect < 0.0) {
				felems[ii - 1]++; /* COUNT ELEM */

				/* TRANSFORM STRESS INTO GLOBAL */
				gstress = stresstransform(elem, estress);

				for (d = 0; d < 2; d++) {
					/* SUM OF SHEAR */
					sumshear[ii - 1][d] -= *(gstress + d);
					*(srate[ii - 1][d] + (elem->sect->loff)) -= *(gstress + d);
					trate[ii - 1][d][elem->sect->type] -= *(gstress + d);

					/* SUM OF RIGID */
					bunbo = *(ddisp + (elem->node[1]->loff) * 6 + d) -
						*(ddisp + (elem->node[0]->loff) * 6 + d);

					if (bunbo != 0.0) {
						K = *(gstress + d) / bunbo;
						Ex = 0.5 * ((elem->node[1])->d[GX] + (elem->node[0])
							->d[GX]);
						Ey = 0.5 * ((elem->node[1])->d[GY] + (elem->node[0])
							->d[GY]);

						sumrigid[ii - 1][d] += fabs(K);
						if (d == 0)
							sumrigidmoment[ii - 1][d] += fabs(K) * Ey;
						if (d == 1)
							sumrigidmoment[ii - 1][d] += fabs(K) * Ex;

						rigid.K[d] = fabs(K);
						if (d == 0) {
							rigid.ielem = elem->code;
							rigid.floor = ii;
							rigid.x = Ex;
							rigid.y = Ey;
						}
						if (d == 1) {
							nrigid++;
							rigids = (struct structrigid*) realloc(rigids,
								nrigid*sizeof(struct structrigid));
							*(rigids + nrigid - 1) = rigid;
						}
					}
				}
				free(gstress);
			}
		}
	}
	free(estress);

	if (fout != NULL)
		fprintf(fout, "\n");
	if (fout != NULL)
		fprintf(fout, "CENTER OF RIGID.\n");
	for (ii = 1; ii <= (FLOORS - 1); ii++) {
		if (sumrigid[ii - 1][1] != 0.0) {
			Lx[ii - 1] = sumrigidmoment[ii - 1][1] / sumrigid[ii - 1][1];
		}
		else
			Lx[ii - 1] = 0.0;

		if (sumrigid[ii - 1][0] != 0.0) {
			Ly[ii - 1] = sumrigidmoment[ii - 1][0] / sumrigid[ii - 1][0];
		}
		else
			Ly[ii - 1] = 0.0;

		ex[ii - 1] = Lx[ii - 1] - Gx[ii - 1];
		ey[ii - 1] = Ly[ii - 1] - Gy[ii - 1];

		if (fout != NULL) {
			fprintf(fout, "%2dF:Lx=%8.5f[m] Ly=%8.5f[m]", ii, Lx[ii - 1],
				Ly[ii - 1]);
			fprintf(fout, " ex=%8.5f[m] ey=%8.5f[m]\n", ex[ii - 1], ey[ii - 1]);
		}
	}

	/* RIGID OF TORTION */
	for (ii = 1; ii <= (FLOORS - 1); ii++) {
		Kr[ii - 1] = 0.0;
	}
	for (jj = 0; jj < nrigid; jj++) {
		rigid = *(rigids + jj);

		ii = rigid.floor;
		Kr[ii - 1] += (rigid.K[0]*(rigid.y - Ly[ii - 1])*(rigid.y - Ly[ii - 1]))
			+ (rigid.K[1]*(rigid.x - Lx[ii - 1])*(rigid.x - Lx[ii - 1]));
	}

	/* RATE OF UNEVENESS */
	if (fout != NULL)
		fprintf(fout, "\n");
	if (fout != NULL)
		fprintf(fout, "RATE OF UNEVENESS.\n");
	for (ii = 1; ii <= (FLOORS - 1); ii++) {
		if (sumrigid[ii - 1][0] != 0.0 && Kr[ii - 1] != 0.0) {
			rex[ii - 1] = sqrt(Kr[ii - 1] / sumrigid[ii - 1][0]);
			Rex[ii - 1] = ey[ii - 1] / rex[ii - 1];
		}
		else
			Rex[ii - 1] = 0.0;

		if (sumrigid[ii - 1][1] != 0.0 && Kr[ii - 1] != 0.0) {
			rey[ii - 1] = sqrt(Kr[ii - 1] / sumrigid[ii - 1][1]);
			Rey[ii - 1] = ex[ii - 1] / rey[ii - 1];
		}
		else
			Rey[ii - 1] = 0.0;

		if (fout != NULL)
			fprintf(fout, "%2dF:Rex=%8.5f Rey=%8.5f\n", ii, Rex[ii - 1],
			Rey[ii - 1]);
	}

	/* SHEAR OF FLOOR */
	if (fout != NULL)
		fprintf(fout, "\n");
	if (fout != NULL)
		fprintf(fout, "SHEAR OF FLOOR.\n");
	for (ii = 1; ii <= (FLOORS - 1); ii++) {
		if (fout != NULL) {
			fprintf(fout, "%2dF:Qx=%9.5f Qy=%9.5f\n", ii, sumshear[ii - 1][0],
				sumshear[ii - 1][1]);
		}
		if (ftxt != NULL) {
			/* fprintf(ftxt," Qx%d %9.5E Qy%d %9.5E",
			 ii,sumshear[ii-1][0],ii,sumshear[ii-1][1]); */
			fprintf(ftxt, " %12.5E %12.5E", sumshear[ii - 1][0],
				sumshear[ii - 1][1]);
		}
	}

	/* SHEAR OF EACH SECTION */
	if (fout != NULL)
		fprintf(fout, "\n");
	if (fout != NULL)
		fprintf(fout, "SHEAR OF SECTION.\n");
	for (ii = 1; ii <= (FLOORS - 1); ii++) {
		if (fout != NULL) {
			fprintf(fout, "%2dF:Sect%4d Qx=%9.5f Qy=%9.5f\n", ii,
				(af->sects + 0)->code, *(srate[ii - 1][0] + 0),
				*(srate[ii - 1][1] + 0));
			for (jj = 1; jj < nsect; jj++) {
				fprintf(fout, "    Sect%4d Qx=%9.5f Qy=%9.5f\n",
					(af->sects + jj)->code, *(srate[ii - 1][0] + jj),
					*(srate[ii - 1][1] + jj));
			}
		}
	}

	/* SHEAR OF EACH SECTION TYPE */
	if (fout != NULL)
		fprintf(fout, "\n");
	if (fout != NULL)
		fprintf(fout, "SHEAR OF TYPE.\n");
	for (ii = 1; ii <= (FLOORS - 1); ii++) {
		if (fout != NULL) {
			fprintf(fout, "%2dF:Type%2d Qx=%9.5f Qy=%9.5f\n", ii, 0,
				trate[ii - 1][0][0], trate[ii - 1][1][0]);

			for (jj = 1; jj < MAXTYPE; jj++) {
				fprintf(fout, "    Type%2d Qx=%9.5f Qy=%9.5f\n", jj,
					trate[ii - 1][0][jj], trate[ii - 1][1][jj]);
			}
		}
	}

	return;
} /* floorvalues */

double *stresstransform(struct owire *elem, double *estress)
	/* ELEMENT STRESS INTO GLOBAL STRESS. */ {
	int ii;
	double **drccos, **t, **tt, *gstress;

	drccos = directioncosine(elem->node[0]->d[GX], elem->node[0]->d[GY],
		elem->node[0]->d[GZ], elem->node[1]->d[GX], elem->node[1]->d[GY],
		elem->node[1]->d[GZ], elem->cangle); /* [DRCCOS] */
	t = transmatrix(drccos); /* [T] */
	tt = matrixtranspose(t, 12); /* [Tt] */

	gstress = matrixvector(tt, estress, 12); /* {F}=[Tt]{f} */

	for (ii = 0; ii <= 2; ii++)
		free(*(drccos + ii));
	free(drccos);
	for (ii = 0; ii <= 11; ii++)
		free(*(t + ii));
	free(t);
	for (ii = 0; ii <= 11; ii++)
		free(*(tt + ii));
	free(tt);

	return gstress;
} /* stresstransform */

void energyoutputtomemory(FILE *ftext,struct arclmframe *af)
/*TRANSLATE ENERGY OUTPUTFILE TEXT INTO MEMORY.*/
{
  char **data,str[256]="\0";
  char *s;
  int i,ii,j,jj,m,n;
  long int ncode;
  double ddata;

  fseek(ftext,0L,SEEK_SET);

  af->ddisp=(double *)malloc(6*(af->nnode)*sizeof(double));
									   /*DISPLACEMENT:6 DIRECTIONS.*/

  while(1)
  {
    data=fgetsbrk(ftext,&n);
    if(!strncmp(*data,"ELEM",4)) break;
    free(data);
  }

  for(i=0;i<(af->nelem);i++) /*DISPLACEMENTS.*/
  {
	if(n!=22) return;

    ddata=strtod(*(data+5),NULL);
    (af->elems+i)->Ee[0]=ddata;

	ddata=strtod(*(data+7),NULL);
    (af->elems+i)->Ep[0]=ddata;

    ddata=strtod(*(data+15),NULL);
    (af->elems+i)->Ee[1]=ddata;

    ddata=strtod(*(data+17),NULL);
    (af->elems+i)->Ep[1]=ddata;

	for(;n>0;n--) free(*(data+n-1));
	free(data);


    data=fgetsbrk(ftext,&n);
  }
  return;
}



