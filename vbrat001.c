int vbrat001(struct arclmframe* af)
{
  DWORDLONG memory0,memory1;
  FILE *fin,*fout;                               /*FILE 8 BYTES*/
  char dir[]=DIRECTORY;                        /*DATA DIRECTORY*/
  char s[800],string[40000];
  clock_t t0/*,t1,t2*/;

  int i,j,ii,jj;

  int nnode,nelem,nshell,nsect,nconstraint;
  long int msize;

  /*ARCLMFRAME*/
  struct osect *sects;
  struct onode *nodes;
  struct onode *ninit;
  struct owire *elems;
  struct oshell *shells;
  struct oconf *confs;
  struct memoryelem* melem;
  struct memoryshell* mshell;
  long int mainoff;
  long int* constraintmain;

  struct gcomponent ginit={0,0,0.0,NULL};
  struct gcomponent *kmtx,*mmtx;                    /*GLOBAL MATRIX*/
  double gdata;

  double* ddisp, * iform;
  double *eigen;
  double **evct;
  double T;

	/*FOR READING ANALYSISDATA*/
	FILE *fdata;
	int nstr, pstr, readflag;
	char **data;
	char *filename;
	int neig = 1;
	int solver = 0;
	/*
	neig=NEIGEN;
	solver=SOLVER;
	*/


	fdata = fgetstofopenII(dir, "r", "analysisdata.txt");
	if (fdata == NULL)
	{
		printf("couldn't open analysisdata.txt\n");

		/*OUTPUT FILE NAME IS SAME AS INPUT*/
		strcpy(filename, (wdraw.childs+1)->inpfile);
		char* dot = strrchr(filename, '.');
		*dot = '\0';

		if(MessageBox(NULL,"DEIGABGENERAL","SOLVER",MB_OKCANCEL)==IDOK)solver=0;
		else if(MessageBox(NULL,"BISECSYLVESTER","SOLVER",MB_OKCANCEL)==IDOK)solver=1;

		/*READ ARCLENGTH PARAMS*/
		getincrement((wmenu.childs+2)->hwnd,&neig,NULL);

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
							filename = *(data + pstr);
						}
						if (!strcmp(*(data + pstr), "NEIG"))
						{
							pstr++;
							neig = (int)strtol(*(data + pstr), NULL, 10);
						}
						if (!strcmp(*(data + pstr), "SOLVER"))
						{
							pstr++;
							solver = (int)strtol(*(data + pstr), NULL, 10);
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



  memory0=availablephysicalmemoryEx("INITIAL:");   /*MEMORY AVAILABLE*/

  fout=fgetstofopen("\0","w",ID_OUTPUTFILE);          /*OUTPUT FILE*/

  //snprintf(fname, sizeof(fname), "%s.%s", filename, "otp");
  //fout = fopen(fname, "w");

  t0=clock();                                        /*CLOCK BEGIN.*/

  double eps=1.0E-16;
  double biseceps=BISECEPS;

  nnode=af->nnode;
  nelem=af->nelem;
  nshell=af->nshell;
  nsect=af->nsect;
  nconstraint=af->nconstraint;

  msize=6*nnode;                           /*SIZE OF GLOBAL MATRIX.*/

  sprintf(string,"NODES=%d ELEMS=%d SHELLS=%d SECTS=%d",nnode,nelem,nshell,nsect);
  errormessage(string);


	/*MEMORY NOT ALLOCATED*/
	//sects = (struct osect*)malloc(nsect * sizeof(struct osect));
	//nodes = (struct onode*)malloc(nnode * sizeof(struct onode));
	//ninit = (struct onode*)malloc(nnode * sizeof(struct onode));
	//elems = (struct owire*)malloc(nelem * sizeof(struct owire));
	//shells = (struct oshell*)malloc(nshell * sizeof(struct oshell));
	//confs = (struct oconf*)malloc(msize * sizeof(struct oconf));
	iform = (double*)malloc(msize * sizeof(double));          /*INITIAL GLOBAL VECTOR.*/
	ddisp = (double*)malloc(msize * sizeof(double));
	melem = (struct memoryelem*)malloc(nelem * sizeof(struct memoryelem));
	mshell = (struct memoryshell*)malloc(nshell * sizeof(struct memoryshell));
	//constraintmain = (long int*)malloc(msize * sizeof(long int));

	//af->sects = sects;
	//af->nodes = nodes;
	//af->ninit = ninit;
	//af->elems = elems;
	//af->shells = shells;
	//af->confs = confs;
	af->ddisp = ddisp;
	af->melem = melem;
	af->mshell = mshell;
	//af->constraintmain = constraintmain;

	/*MEMORY ALREADY ALLOCATED*/
	sects = af->sects;
	nodes = af->nodes;
	ninit = af->ninit;
	elems = af->elems;
	shells = af->shells;
	confs = af->confs;
	//ddisp = af->ddisp;
	//melem = af->melem;
	//mshell = af->mshell;
	constraintmain = af->constraintmain;


  /***GLOBAL MATRIX***/
  kmtx=(struct gcomponent *)malloc(msize*sizeof(struct gcomponent));
  mmtx=(struct gcomponent *)malloc(msize*sizeof(struct gcomponent));
  for(i=0;i<msize;i++)
  {
    (kmtx+i)->down=NULL;            /*GLOBAL MATRIX INITIALIZATION.*/
    (mmtx+i)->down=NULL;
  }

  for(i=1;i<=msize;i++)
  {
	ginit.m=(unsigned short int)i;
	/*ginit.n=(unsigned int)i;*/
	*(kmtx+(i-1))=ginit;
	*(mmtx+(i-1))=ginit;
  }
  comps=msize;

  /***GLOBAL VECTOR***/
  evct=(double **)malloc(neig*sizeof(double *));    /*EIGEN VECTORS*/
  eigen=(double *)malloc(neig*sizeof(double));       /*EIGEN VALUES*/
  if(evct==NULL || eigen==NULL) return 0;
  for(i=0;i<neig;i++)
  {
	*(evct+i)=(double *)malloc(msize*sizeof(double));
	for(j=0;j<msize;j++) *(*(evct+i)+j)=0.0;
  }

  iform = (double*)malloc(msize * sizeof(double));          /*INITIAL GLOBAL VECTOR.*/
  ddisp = (double*)malloc(msize * sizeof(double));

  initialformCR(ninit, iform, nnode);           /*ASSEMBLAGE FORMATION.*/
  initialformCR(nodes, ddisp, nnode);           /*ASSEMBLAGE FORMATION.*/

  GetAsyncKeyState(VK_LBUTTON);                  /*CLEAR KEY LEFT.*/
  GetAsyncKeyState(VK_RBUTTON);                  /*CLEAR KEY RIGHT.*/

  errormessage("VBRAT001:EIGEN ANALYSIS OF VIBRATION.");
  availablephysicalmemoryEx("REMAIN:");            /*MEMORY AVAILABLE*/
  laptime("ASSEMBLING GLOBAL MATRIX.",t0);

  assemelem(elems, melem, nelem, constraintmain, mmtx, kmtx, iform, ddisp, NULL, NULL);
  assemshell(shells, mshell, nshell, constraintmain, mmtx, kmtx, iform, ddisp, NULL, NULL);

  laptime("GLOBAL MATRIX ASSEMBLED.",t0);


  if(solver==0)
  {
	deigabgeneral(mmtx,kmtx,confs,msize,neig,neig,eps,eigen,evct);
  }
  else
  {
	bisecsylvester(mmtx,kmtx,confs,msize,neig,neig,biseceps,eigen,evct);
  }
  laptime("EIGEN COMPLETED.",t0);

  for(i=0;i<neig;i++)
  {
	//outputmode(*(evct+i),fout,nnode,ninit);

	if (solver==0)
	{
		if (fout != NULL)fprintf(fout,"DEIGABGENERAL EIGEN VALUE = %.8f ", *(eigen + i));
		if (*(eigen + i) > 0.0)
		{
			T = 2.0 * PI * sqrt(*(eigen + i));
			sprintf(string, "PERIOD = %.8f [sec] FREQUENCY = %.8f [Hz]",T,1/T);
			if (fout != NULL)fprintf(fout, "%s\n", string);
		}
		else
		{
			if (fout != NULL)
				fprintf(fout, "ERROR:EIGEN VALUE NEGATIVE.\n");
		}
	}
	else // bisecsylvester
	{
		if (fout != NULL)fprintf(fout,"BISECSYLVESTER EIGEN VALUE = %.8f ", 1.0 / (*(eigen + i))); // bisecsylvester
		if (*(eigen + i) > 0.0)
		{
			T = 2.0 * PI * sqrt(*(eigen + i));
			sprintf(string, "PERIOD = %.8f [sec] FREQUENCY = %.8f [Hz]",T,1/T);
			if (fout != NULL)fprintf(fout, "%s\n", string);

		}
		else
		{
			if (fout != NULL)
				fprintf(fout, "ERROR:EIGEN VALUE NEGATIVE.\n");
		}
	}

	//if (fout != NULL)fprintf(fout, "EIGEN VECTOR %ld\n",(i + 1));
	for (ii = 0; ii < nnode; ii++)
	{
		if (fout != NULL)
			fprintf(fout,
			"%4ld %11.7f %11.7f %11.7f %11.7f %11.7f %11.7f\n",
			(nodes + ii)->code, *(*(evct + i) + 6*ii + 0),
			*(*(evct + i) + 6*ii + 1),
			*(*(evct + i) + 6*ii + 2),
			*(*(evct + i) + 6*ii + 3),
			*(*(evct + i) + 6*ii + 4),
			*(*(evct + i) + 6*ii + 5));
	}

  }

  updatemode(af,*(evct+0)); /*FORMATION UPDATE.*/

  af->nlaps=neig;
  af->eigenval=eigen;
  af->eigenvec=evct;


  gfree(kmtx,nnode); /*FREE GLOBAL MATRIX.*/
  gfree(mmtx,nnode); /*FREE GLOBAL MATRIX.*/

  errormessage(" ");
  errormessage("COMPLETED.");
  fprintf(fout,"COMPLETED.\n");
  fclose(fout);

  memory1=availablephysicalmemoryEx("REMAIN:");
  sprintf(string,"CONSUMPTION:%ld[BYTES]",(memory0-memory1));
  errormessage(string);

  return 1;
}/*bclng001*/






#if 0
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
#endif
