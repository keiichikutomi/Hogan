


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

extern int beziersurfaceII(int m,int n,double *x,double *y,double *z,
					double u,double v,struct onode *node);



struct outputdata
{

};

void conjugategradientaf(struct arclmframe *af)     /*CONJUGATE*/
/*OPTIMIZE ORGAN WITH CONJUGATE GRADIENT.*/
{
  FILE *ftxt;
  char str[500];
  int i,j,k,m,n;
  int ii,jj;

  int nnode,nelem,nshell;
  double df,f1,f2,fa,ftarget;
  double c1,c2,c3,alpha,beta,gamma,vsize,eps;

  double dfact; /*COORDINATES*/
  double *u1,*u2,*ua,*fgrad1,*fgrad2; /*GRADIENT VECTOR*/

  double xi,tau,test; /*FOR LINE SEARCH(ARMIJO RULE)*/

  /*BEZIER SURFACE*/
  FILE *fbezier;
  char **data/*,str[256]="\0"*/;
  double ddata;
  int ndata;
  long int ncode;

  /*CONTROL POINTS*/
  int ncontrol;
  int n1,n2;
  double *x,*y,*z;              /*CURRENT CONTROLE POINTS*/
  double *xini,*yini,*zini;     /*INITIAL CONTROLE POINTS*/
  int *fp;                      /*ICONF OF CONTROLE POINTS*/

  /*ALL NODES*/
  struct onode node;
  double *u,*v;                 /*U-V COORDINATES*/


  ftarget=5.0;
  gamma=0.01;
  dfact=0.1;
  eps=0.001;
  /*PARAMETERS FOR LINE SEARCH(ARMIJO RULE)*/
  xi=0.001;/*0.0<=xi<=1.0*/
  tau=0.9; /*0<tau<1*/


  nnode=af->nnode;
  nelem=af->nelem;
  nshell=af->nshell;

  /*CONTROL SURFACE FOR DESIGN VARIABLES*/

  /*CREATE INITIAL BEZIER SURFACE*/
  fbezier=fopen("bezier.txt","r");   /*BEZIER SURFACE DATA*/
  fseek(fbezier,0L,SEEK_SET);

  data=fgetsbrk(fbezier,&ndata);
  n1=strtol(*(data+0),NULL,10);
  n2=strtol(*(data+1),NULL,10);
  if(nnode!=strtol(*(data+2),NULL,10)) return;

  /*CONTROLE POINTS*/
  ncontrol=(n1+1)*(n2+1);
  fp=(int *)malloc(ncontrol*sizeof(int));
  x=(double *)malloc(ncontrol*sizeof(double));
  y=(double *)malloc(ncontrol*sizeof(double));
  z=(double *)malloc(ncontrol*sizeof(double));
  zini=(double *)malloc(ncontrol*sizeof(double));

  fgrad1=(double *)malloc(ncontrol*sizeof(double));
  fgrad2=(double *)malloc(ncontrol*sizeof(double));
  u1=(double *)malloc(ncontrol*sizeof(double));
  u2=(double *)malloc(ncontrol*sizeof(double));
  ua=(double *)malloc(ncontrol*sizeof(double));


  for(i=0;i<ncontrol;i++)
  {
	data=fgetsbrk(fbezier,&ndata);
	if(ndata!=3) return;
	fp[i]=0;
	x[i]=strtod(*(data+0),NULL);
    y[i]=strtod(*(data+1),NULL);
	z[i]=strtod(*(data+2),NULL);
	zini[i]=z[i];

    for(;ndata>0;ndata--) free(*(data+ndata-1));
	free(data);
  }

  /*U-V COORDINATE OF EACH NODE*/
  u=(double *)malloc(nnode*sizeof(double));
  v=(double *)malloc(nnode*sizeof(double));
  for(i=0;i<nnode;i++)
  {
	data=fgetsbrk(fbezier,&ndata);
	if(ndata!=5) return;

	ncode=strtol(*(data+1),NULL,10);
	if(ncode!=(af->nodes+i)->code) return;

	*(u+i)=strtod(*(data+3),NULL);
	*(v+i)=strtod(*(data+4),NULL);

    for(;ndata>0;ndata--) free(*(data+ndata-1));
    free(data);
  }
  fclose(fbezier);




  ///*INITIAL*///
  for(ii=0;ii<nnode;ii++)
  {
	node=*(af->nodes+ii);
	beziersurfaceII(n1,n2,x,y,z,*(u+ii),*(v+ii),&node);

	(af->nodes+ii)->d[0]=node.d[0];
	(af->nodes+ii)->d[1]=node.d[1];
	(af->nodes+ii)->d[2]=node.d[2];
  }
  f1=arclmStatic(af);
  ///*INITIAL*///

  ///*SENSITIVITY*///
  for(i=0;i<ncontrol;i++)
  {
	if(*(fp+i)==0)
	{
	  for(j=0;j<ncontrol;j++)
	  {
		if(i==j)
		{
		  z[j]=zini[j]+dfact;
		}
		else
		{
		  z[j]=zini[j];
		}
	  }

	  for(ii=0;ii<nnode;ii++)
      {
		node=*(af->nodes+ii);
		beziersurfaceII(n1,n2,x,y,z,*(u+ii),*(v+ii),&node);

		(af->nodes+ii)->d[0]=node.d[0];
		(af->nodes+ii)->d[1]=node.d[1];
		(af->nodes+ii)->d[2]=node.d[2];
	  }
	  df=arclmStatic(af);
	}
	else
	{
	  df=f1;
	}
	*(fgrad1+i)=gamma*(f1-df)/dfact;/*NEGATIVE GRADIENT*/
	fprintf(ftxt,"CONTROLE POINT[%d] NEGATIVE GRADIENT=%9.5f\n",i,(f1-df)/dfact);
  }
  ///*SENSITIVITY*///
  for(i=0;i<ncontrol;i++)
  {
	*(u1+i)=*(fgrad1+i);
  }


  k=0;
  while(1)/*OPTIMIZATION BEGIN.*/
  {
	k++;
	sprintf(str,"%d",k);
	SetDlgItemText((wmenu.childs+2)->hwnd,ID_LAPS,str);
	/*BACKTRACKING LINE SEARCH : TEST STEP*/

	///*UPDATE*///
	for(i=0;i<ncontrol;i++)
	{
	  z[i]=zini[i]+*(u1+i);
	}
	for(ii=0;ii<nnode;ii++)
	{
	  node=*(af->nodes+ii);
	  beziersurfaceII(n1,n2,x,y,z,*(u+ii),*(v+ii),&node);

	  (af->nodes+ii)->d[0]=node.d[0];
	  (af->nodes+ii)->d[1]=node.d[1];
	  (af->nodes+ii)->d[2]=node.d[2];
	}
	fa=arclmStatic(af);
	///*UPDATE*///

	///*SENSITIVITY*///
	for(i=0;i<ncontrol;i++)
	{
	  if(*(fp+i)==0)
	  {
		for(j=0;j<ncontrol;j++)
		{
		  if(i==j)
		  {
			z[j]=zini[j]+*(u1+i)+dfact;
		  }
		  else
		  {
			z[j]=zini[j]+*(u1+i);
		  }
		}
		for(ii=0;ii<nnode;ii++)
		{
		  node=*(af->nodes+ii);
		  beziersurfaceII(n1,n2,x,y,z,*(u+ii),*(v+ii),&node);

		  (af->nodes+ii)->d[0]=node.d[0];
		  (af->nodes+ii)->d[1]=node.d[1];
		  (af->nodes+ii)->d[2]=node.d[2];
		}
		df=arclmStatic(af);
	  }
	  else
	  {
		df=fa;
	  }
	  *(fgrad2+i)=gamma*(fa-df)/dfact;/*NEGATIVE GRADIENT*/
	  fprintf(ftxt,"CONTROLE POINT[%d] NEGATIVE GRADIENT=%9.5f\n",i,(fa-df)/dfact);
	}
	///*SENSITIVITY*///
	for(i=0;i<ncontrol;i++)
	{
	  *(ua+i)=-*(fgrad2+i)+*(fgrad1+i); /*DIRECTION DERIVATIVE OF GRADIENT*/
	}


	/*BACKTRACKING LINE SEARCH : INITIAL STEP SIZE*/
	c1=0.0;
	c2=0.0;
	for(i=0; i<m; i++)
	{
	  c1+=(*(fgrad1+i))*(*(fgrad1+i));
	  c2+=(*(u1+i))*(*(ua+i));
	}
	alpha=c1/c2;

	/*ARMIJO RULE*/
	while(1)
	{

	  ///*UPDATE*///
	  for(i=0;i<ncontrol;i++)
	  {
		z[i]=zini[i]+alpha**(u1+i);
	  }
	  for(ii=0;ii<nnode;ii++)
	  {
		node=*(af->nodes+ii);
		beziersurfaceII(n1,n2,x,y,z,*(u+ii),*(v+ii),&node);

		(af->nodes+ii)->d[0]=node.d[0];
		(af->nodes+ii)->d[1]=node.d[1];
		(af->nodes+ii)->d[2]=node.d[2];
	  }
	  f2=arclmStatic(af);
	  ///*UPDATE*///

	  /*ARMIJO END.*/
	  test=0.0;
	  for(i=0;i<ncontrol;i++)
	  {
		*(fgrad2+i)=*(fgrad1+i)-alpha**(ua+i);
		test+=*(fgrad2+i)**(u1+i);
	  }
	  test=f2-f1+abs(xi*alpha*test);
	  if(test<=0.0 || (alpha<=0.01 && f1>f2))
	  {
		sprintf(str,"LINE SEARCH : ALPHA= %e f1= %e f2= %e TEST= %8.3f\n",alpha,f1,f2,test);
		fprintf(ftxt,"%s",str);
		break;
	  }
	  else
	  {
		alpha*=tau;
		sprintf(str,"LINE SEARCH : ALPHA= %e f1= %e f2= %e TEST= %8.3f\n",alpha,f1,f2,test);
		fprintf(ftxt,"%s",str);
	  }
	}

	///*SENSITIVITY*///
	for(i=0;i<ncontrol;i++)
	{
	  if(*(fp+i)==0)
	  {
		for(j=0;j<ncontrol;j++)
		{
		  if(i==j)
		  {
			z[j]=zini[j]+alpha**(u1+i)+dfact;
		  }
		  else
		  {
			z[j]=zini[j]+alpha**(u1+i);
		  }
		}
		for(ii=0;ii<nnode;ii++)
		{
		  node=*(af->nodes+ii);
		  beziersurfaceII(n1,n2,x,y,z,*(u+ii),*(v+ii),&node);

		  (af->nodes+ii)->d[0]=node.d[0];
		  (af->nodes+ii)->d[1]=node.d[1];
		  (af->nodes+ii)->d[2]=node.d[2];
		}
		df=arclmStatic(af);
	  }
	  else
	  {
		df=f2;
	  }
	  *(fgrad2+i)=gamma*(f2-df)/dfact;/*NEGATIVE GRADIENT*/
	  fprintf(ftxt,"CONTROLE POINT[%d] NEGATIVE GRADIENT=%9.5f\n",i,(f2-df)/dfact);
	}
	///*SENSITIVITY*///



	sprintf(str,"STEP= %d OBJECTIVE FUNCTION= %.5f GRADIENT SIZE= %.5f\n",k,f2,vsize);
	fprintf(ftxt,"%s",str);

	/*OPTIMIZATION END.*/
	vsize=vectorlength(fgrad2,ncontrol);
	if(vsize<eps || f2<=ftarget)
	{
	  sprintf(str,"COMPLETED\n");
	  fprintf(ftxt,"%s",str);
	  break;
	}

	/*UPDATE CONJUGATE GRADIENT DIRECTION*/
	c3=0.0;
	for(i=0; i<m; i++)
	{
	  c3+=*(fgrad2+i)**(fgrad2+i);
	}
	beta=c3/c1;/*FLETCHER-REEVES*/
	for(i=0;i<ncontrol;i++)
	{
	  *(u2+i)=*(fgrad2+i)+beta**(u1+i);
	}

	/*GO TO NEXT LAP*/
	for(i=0;i<ncontrol;i++)
	{
	  zini[i]=z[i];
	  *(u1+i)=*(u2+i);
	  *(fgrad1+i)=*(fgrad2+i);
	}
	f1=f2;

	test=0.0;
	for(i=0;i<ncontrol;i++)
	{
	  test+=*(u1+i)**(fgrad1+i);/*CHECK DIRECTION : COMPARE CONJUGATE GRADIENT WITH STEEPEST GRADIENT*/
	}
	if(test<=0)/*CONJUGATE GRADIENT IS NOT DESCENT DIRECTION*/
	{
	  for(i=0;i<ncontrol;i++)
	  {
		 *(u1+i)=*(fgrad1+i);/*THE STEEPEST DESCENT METHOD*/
	  }
	  sprintf(str,"THE STEEPEST DESCENT METHOD AT STEP %d\n",k+1);
	  fprintf(ftxt,"%s",str);
	}

	/*OUTPUT*/
	/*
	sprintf(str,"hogtxt_opt%d.inp",k);
	fresult=fopen(str,"w");
	if(fresult==NULL) break;
	saveorganization(fresult,&((wdraw.childs+1)->org),
					&((wdraw.childs+1)->vparam));
	fclose(fresult);
	*/
  }
  fclose(ftxt);

  return;
}/*conjugategradientaf*/





int arclmStaticOptimization(struct arclmframe* af, struct outputdata* op)
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
	int USINGEIGENFLAG = 1;


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
	int ESMODE = 0;


	double eps = 1e-8;
	double* epsddisp, * epsgvct, * epslambda;
	double* re, * rp;
	double* Kinvre, * Kinvrp, *Kevct;
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
	constraintvct = (double*)malloc(csize * sizeof(double));

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
	initialform(nodes, lastddisp, nnode);           /*ASSEMBLAGE FORMATION.*/


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
	for (i = 0; i < csize; i++)
	{
		*(lastlambda + i) = *(lambda + i);
	}


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

		Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
		//Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;

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

#if 1
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
					nlap, iteration, loadfactor, residual, /*determinant*/lenR, sign, EXTENDEDFLAG, eigen, volume, 0.0);
		fprintf(ffig, "%s", string);
		errormessage(string);

#if 1
		   /*if(nlap == 20)
			{
			  EXTENDEDFLAG=1;
			}*/


		if(iteration == 1 && (sign - lastsign) != 0 && nlap != beginlap)
		{
			/*PIVOT ZERO LINE CHECK*/
		/*	for (ii = 0; ii < msize; ii++)
			{
				if ((confs + ii)->iconf == 0 && *(lastpivot + ii) * ((gmtx + ii)->value) < 0)
				{
					m = ii;
					sprintf(string, "BUCKLING DITECTED LAP: %4d ITER: %2d LINE: %5ld ", nlap, iteration, ii);
					fprintf(fbcl, "%s\n", string);
					errormessage(string);
				}
			}      */


			if(pinpointmode == 1)
			{
			  BISECTIONFLAG=1;
			}
			if(pinpointmode == 2 && nlap >70)
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


					if(iteration==1)
					{
						evct_R = (double*)malloc(msize * sizeof(double));
						for (ii = 0; ii < msize; ii++)
						{
							*(evct_R + ii) = *(dup + ii);
						}
						lenR = vectorlength(evct_R, msize);

						for (ii = 0; ii < msize; ii++)
						{
							*(evct_R + ii) *= 1.0/lenR;
						}
						/*LDL^T MODE*/
					}


					dbgvct(evct_R,msize,6,"EVCT");

					fprintf(fbcl, "CRITICAL EIGEN VECTOR : LINE = %5ld dm = %12.9f EIGENVALUE=%12.9f\n", m, dm, eigen);
					for (ii = 0; ii < nnode; ii++)
					{
						fprintf(fbcl, "%5ld %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f\n", (nodes + ii)->code,
							*(evct_R + (6 * ii + 0)), *(evct_R + (6 * ii + 1)), *(evct_R + (6 * ii + 2)),
							*(evct_R + (6 * ii + 3)), *(evct_R + (6 * ii + 4)), *(evct_R + (6 * ii + 5)));
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
					Kinvre = (double*)malloc(msize * sizeof(double));
					Kinvrp = (double*)malloc(msize * sizeof(double));
					Kevct = (double*)malloc(msize * sizeof(double));
					for (ii = 0; ii < msize; ii++)
					{
						*(re + ii) = 0.0;
						*(rp + ii) = 0.0;
						*(Kinvre + ii) = 0.0;
						*(Kinvrp + ii) = 0.0;
						*(Kevct + ii) = 0.0;
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
						gform = extractdisplacement(elem, epsddisp);
						drccos = updatedrccos(drccosinit, gforminit, gform);
						eform = extractlocalcoord(gform,drccos,nnod);

						edisp = extractdeformation(eforminit, eform, nnod);                	/*{Ue}*/

						T = transmatrixIII(drccos, nnod);									/*[T].*/
						HPT = transmatrixHPT(eform, edisp, T, nnod);

						einternal = matrixvector(Ke, edisp, 6 * nnod);          			/*{Fe}=[Ke]{Ue}.*/
						//einternal = assemshelleinternal(&elem);

						Kt = assemtmtxCR(Ke, eform, edisp, einternal, T, HPT, nnod);	/*TANGENTIAL MATRIX[Kt].*/
						symmetricmtx(Kt, 6*nnod);											/*SYMMETRIC TANGENTIAL MATRIX[Ksym].*/

						for (ii = 0; ii < 6*nnod; ii++)
						{
							for (jj = 0; jj < 6*nnod; jj++)
							{
								*(re + *(constraintmain + *(loffset + ii))) += *(*(Kt + ii) + jj) * *(due + *(constraintmain + *(loffset + jj)));
								*(rp + *(constraintmain + *(loffset + ii))) += *(*(Kt + ii) + jj) * *(dup + *(constraintmain + *(loffset + jj)));
								*(Kevct + *(constraintmain + *(loffset + ii))) += *(*(Kt + ii) + jj) * *(evct_R + *(constraintmain + *(loffset + jj)));
							}
						}

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
							inputnode(epsddisp, shell.node[ii]);
						}
						drccos = shelldrccos(shell);
						gform = extractshelldisplacement(shell, epsddisp);                     /*{Xg+Ug}*/
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
								*(Kevct + *(constraintmain + *(loffset + ii))) += *(*(Kt + ii) + jj) * *(evct_R + *(constraintmain + *(loffset + jj)));
							}
						}

						for(ii=0;ii<ngp;ii++)freematrix(*(C+ii), nstress);
						free(C);
						for(ii=0;ii<ngp;ii++)freematrix(*(B+ii), nstress);
						free(B);

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

						freematrix(Kp, 6 * nnod);
						freematrix(Kt, 6 * nnod);
					}

					dbgvct(re,msize,6,"RE");
					dbgvct(funbalance,msize,6,"FUNB");


					for (ii = 0; ii < msize; ii++)/*UPDATE FORM FOR DIRECTIONAL DERIVATIVE*/
					{
						if ((confs + ii)->iconf != 0)
						{
							*(re + ii) = 0.0;
						}
						else
						{
							*(re + ii) = -(*(re + ii)-*(funbalance + ii))/eps;
						}
					}
					dbgvct(re,msize,6,"epsRE");

					Vector eigenRe = Vector::Zero((msize+csize));
					for(i=0;i<msize;i++)eigenRe(i)=*(re+i);
					Vector eigenKinvre = solver.solve(eigenRe);
					for(i=0;i<(msize+csize);i++)*(Kinvre+i)=eigenKinvre(i);
					if (solver.info() != Eigen::Success)return 1;

					dbgvct(re,msize,6,"Kinvre");



					dbgvct(rp,msize,6,"RP");
					dbgvct(fbaseload,msize,6,"FBASE");

					for (ii = 0; ii < msize; ii++)/*UPDATE FORM FOR DIRECTIONAL DERIVATIVE*/
					{
						if ((confs + ii)->iconf != 0)
						{
							*(rp + ii) = 0.0;
						}
						else
						{
							*(rp + ii) = -(*(rp + ii)-*(fbaseload + ii))/eps;
						}
					}
					dbgvct(rp,msize,6,"epsRP");

					Vector eigenRp = Vector::Zero((msize+csize));
					for(i=0;i<msize;i++)eigenRp(i)=*(rp+i);
					Vector eigenKinvrp = solver.solve(eigenRp);
					for(i=0;i<(msize+csize);i++)*(Kinvrp+i)=eigenKinvrp(i);
					if (solver.info() != Eigen::Success)return 1;

					dbgvct(re,msize,6,"Kinvrp");


					lenR = vectorlength(evct_R, msize);

					evctre = dotproduct(Kinvre,evct_R,msize);
					evctrp = dotproduct(Kinvrp,evct_R,msize);
					fprintf(fbcl, "evctre=%12.15f evctrp=%12.15f\n", evctre, evctrp);

					loadlambda = (lenR-evctre)/evctrp;
					eigen = dotproduct(evct_R,Kevct,msize);
					eigen *= 1/(lenR*lenR);
					fprintf(fbcl, "LOADLAMBDA=%12.9f\n", loadlambda);



					for (ii = 0; ii < msize; ii++)
					{
						*(gvct + ii) = *(dup + ii) * loadlambda + *(due + ii);
					}

					for (ii = 0; ii < msize; ii++)
					{
						*(evct_R + ii) = (*(Kinvrp + ii) * loadlambda + *(Kinvre + ii)) ;
					}


					free(rp);
					free(re);
					free(Kinvrp);
					free(Kinvre);
					free(Kevct);

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
			/*
				if((residual < tolerance && fabs(eigen) < 1e-8)  && iteration != 1)
				{
					EXTENDEDFLAG = 0;
					//eigen=0.0;
					if(ESMODE==0)
					{

                    }

					nlap++;
					iteration = 0;
				}
			*/
				clearwindow(*(wdraw.childs+1));
				drawarclmframe((wdraw.childs+1)->hdcC,(wdraw.childs+1)->vparam,*af,0,ONSCREEN);
				overlayhdc(*(wdraw.childs + 1), SRCPAINT);
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
					  sprintf(string,"CONVERGENCE ERROR(RESIDUAL = %f). ANALYSIS TERMINATED.\n", sign, residual);
					  errormessage(string);
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
			dbgstr(string);

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
	free(lastmelem);
	free(lastmshell);

	errormessage(" ");
	errormessage("COMPLETED.");

	memory1=availablephysicalmemoryEx("REMAIN:");
	sprintf(string,"CONSUMPTION:%ld[BYTES]",(memory0-memory1));
	errormessage(string);

	return 0;
}
/*arclmCR*/

