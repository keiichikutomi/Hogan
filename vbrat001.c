int vbrat001(struct arclmframe* af);

int vbrat001(struct arclmframe* af)
{
  DWORDLONG memory0,memory1;
  FILE *fin,*fout;                               /*FILE 8 BYTES*/
  char dir[]=DIRECTORY;                        /*DATA DIRECTORY*/
  char s[800],string[400],fname[256];
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
  long int* constraintmain;
  double* ddisp, * iform;

  struct gcomponent ginit={0,0,0.0,NULL};
  struct gcomponent *kmtx,*mmtx;                    /*GLOBAL MATRIX*/

	/*FOR READING ANALYSISDATA*/
	FILE *fdata;
	int nstr, pstr, readflag;
	char **data;
	char filename[256];
	char* dot;
	int neig = 1;
	int solver = 0;
	/*
	neig=NEIGEN;
	solver=SOLVER;
	*/

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

	msize=6*nnode;/*SIZE OF GLOBAL MATRIX.*/


	strcpy(filename, (wdraw.childs+1)->inpfile);
	dot = strrchr(filename, '.');
	*dot = '\0';

	fdata = fgetstofopenII(dir, "r", "analysisdata.txt");

	if (fdata == NULL)
	{
		errormessage("couldn't open analysisdata.txt\n");
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
							if(strstr(*(data + pstr),filename)==NULL)
							{
								readflag = 0;
								if(MessageBox(NULL,"DEIGABGENERAL","SOLVER",MB_OKCANCEL)==IDOK)solver=0;
								else if(MessageBox(NULL,"BISECSYLVESTER","SOLVER",MB_OKCANCEL)==IDOK)solver=1;
								getincrement((wmenu.childs+2)->hwnd,&neig,NULL);
							}
							else
							{
								strcpy(filename, *(data + pstr));
							}
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

	snprintf(fname, sizeof(fname), "%s.%s", filename, "otp");
	fout = fopen(fname, "w");

	t0=clock();                                        /*CLOCK BEGIN.*/



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



  /***GLOBAL MATRIX***/
  kmtx=(struct gcomponent *)malloc(msize*sizeof(struct gcomponent));
  mmtx=(struct gcomponent *)malloc(msize*sizeof(struct gcomponent));
  if(kmtx==NULL || mmtx==NULL) return 0;
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

  GetAsyncKeyState(VK_LBUTTON);                  /*CLEAR KEY LEFT.*/
  GetAsyncKeyState(VK_RBUTTON);                  /*CLEAR KEY RIGHT.*/

  errormessage("VBRAT001:EIGEN ANALYSIS OF VIBRATION.");
  availablephysicalmemoryEx("REMAIN:");            /*MEMORY AVAILABLE*/
  //laptime("ASSEMBLING GLOBAL MATRIX.",t0);

  assemelem(elems, melem, nelem, constraintmain, mmtx, kmtx, iform, ddisp, NULL, NULL);
  assemshell(shells, mshell, nshell, constraintmain, mmtx, kmtx, iform, ddisp, NULL, NULL);

  //laptime("GLOBAL MATRIX ASSEMBLED.",t0);


  /*EIGEN ANALYSIS START.       	*/
  /*([A]-eigen*[B])*{evct}=0    	*/
  /*                            	*/
  /*VIBRATION ANALYSIS				*/
  /*([M]-1/(omega^2)*[K])*{evct}=0		*/
  /*A:[M]							*/
  /*B:[K]							*/
  /*eigen=1/(omega^2)  				*/
  /*omega=1/(eigen^0.5)             */
  /*T=2*PI/omega=2*PI*eigen^0.5     */
  /*f=1/(2*PI*eigen^0.5)		    */
  /*0<eigen                     	*/
  /*			                   	*/


  double **evct, *eigen; /* EIGEN VECTORS,EIGEN VALUES */
  double T;
  double eps=1.0E-10;
  double biseceps=1.0E-10;

  evct=(double **)malloc(neig*sizeof(double *));    /*EIGEN VECTORS*/
  eigen=(double *)malloc(neig*sizeof(double));       /*EIGEN VALUES*/
  if(evct==NULL || eigen==NULL) return 0;
  for(i=0;i<neig;i++)
  {
	*(evct+i)=(double *)malloc(msize*sizeof(double));
	for(j=0;j<msize;j++) *(*(evct+i)+j)=0.0;
  }

  bisecgeneral(mmtx,1.0,kmtx,-1.0,confs,msize,neig,biseceps,eigen,evct,0.0,1.0);

  laptime("EIGEN COMPLETED.",t0);


  /*OUTPUT*/
  if(fout!=NULL)
  {
	  for(i=0;i<neig;i++)
	  {
		//outputmode(*(evct+i),fout,nnode,ninit);
		//ninit FOR NODE ID
		//NODE:%5ld {dU}= %12.5E %12.5E %12.5E %12.5E %12.5E %12.5E\n

		fprintf(fout,"LAP = %d MODE = %d GENERALIZED EIGENVALUE = %e STANDARD EIGENVALUE = %e ", af->nlaps, (i+1), *(eigen + i), 0.0);

		if (*(eigen + i) > 0.0)
		{
			T = 2.0 * PI * sqrt(*(eigen + i));
			fprintf(fout, "PERIOD = %.8f FREQUENCY = %.8f\n",T,1/T);
		}
		else
		{
			fprintf(fout, "ERROR:EIGEN VALUE NEGATIVE.\n");
		}

		for (ii = 0; ii < msize; ii++)
		{
			*(*(evct + i) + ii) = *(*(evct + i) + *(constraintmain + ii));
		}

		//if (fout != NULL)fprintf(fout, "EIGEN VECTOR %ld\n",(i + 1));
		for (ii = 0; ii < nnode; ii++)
		{
			fprintf(fout,
			"%4ld %e %e %e %e %e %e\n",
			(nodes + ii)->code, *(*(evct + i) + 6*ii + 0),
			*(*(evct + i) + 6*ii + 1),
			*(*(evct + i) + 6*ii + 2),
			*(*(evct + i) + 6*ii + 3),
			*(*(evct + i) + 6*ii + 4),
			*(*(evct + i) + 6*ii + 5));
		}
	  }
  }
  fclose(fout);

  //af->nlaps=neig;
  af->eigenval=eigen;
  af->eigenvec=evct;

  /*EIGEN ANALYSIS END.*/

  //updatemode(af,*(evct+0)); /*FORMATION UPDATE.*/

  gfree(kmtx,nnode); /*FREE GLOBAL MATRIX.*/
  gfree(mmtx,nnode); /*FREE GLOBAL MATRIX.*/

  errormessage(" ");
  errormessage("COMPLETED.");

  memory1=availablephysicalmemoryEx("REMAIN:");
  sprintf(string,"CONSUMPTION:%ld[BYTES]",(memory0-memory1));
  errormessage(string);

  return 1;
}/*vbrat001*/
