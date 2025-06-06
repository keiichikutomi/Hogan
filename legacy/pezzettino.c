/*PEZZETTINO STATIC NONLINEAR ANALYSIS*/  /*UNDER CONSTRUCTION*/
int pezzettino_staticnonlinear(struct arclmframe *af,int idinput);
double **calculateDeformationStiffnessMatrixKd_CR(struct owire elem);
double **calculateTransformationMatrixS_CR(struct owire elem,double **R);
double **calculateRotationMatrixKr_CR(struct owire elem,struct corotational *corot);
double **assemElementTangentStiffness_CR(struct owire elem,struct corotational *corot,double **Kd,double **S,double **Kr);
void updateElementProperty_CR(struct owire *elem,struct corotational *corot,double *gvct);
double *calculateNodalForce_CR(struct owire *elem,struct corotational *corot,struct memoryelem *melem);

/*CR PROP*/
struct corotational{double mode[6];   //DEFORMATION MODE VECTOR AND PhiAsymX LENGTH0
					double phiAsymX;
					double length0;
					double force[6];  //INTERNAL ELEMENT FORCE
					double basis[9];  //ELEMENT BASIS OF COROTATIONAL BEAM
				   };

/*CR PROP*/
void initializecorotational(struct arclmframe *af,struct corotational *corot);
void corotationaloutputtomemory(FILE *ftext,struct arclmframe *af,struct corotational *corot);
void outputcorotational(FILE *fout,struct arclmframe *af,struct corotational *corot);


/*fstart*/
int saveloadasfstart(char *fname,struct arclmframe *af);

/***FOR PEZZETTINO 003 FILE TYPE (-inl2 style)***/
int saveasarclm2(char *fname,struct arclmframe *af);
void inputtexttomemory2(FILE *ftext,struct arclmframe *af);
void inputinit2(FILE *fin,int *nnode,int *nelem,int *nsect);
void readsect2(FILE *fin,struct osect *sect);

/***FOR PEZZETTINO 003 INCREMENTAL ANALYSIS***/
int updateorganization(struct organ *org);

extern struct arclmframe arc;   /*GLOBAL ARCLM FRAME.*/
extern struct arclmframe arci;  /*NULL ARCLM FRAME.*/
extern struct arclmframe arcx;  /*ARCLM FRAME FOR X LOAD.*/
extern struct arclmframe arcy;  /*ARCLM FRAME FOR Y LOAD.*/

#define FILEFORCHANGESECT  "changesect_test.inp"  /*file for test*/

int pezzettino_staticnonlinear(struct arclmframe *af,int idinput)
/*
  PEZZETTINO STATIC NONLINEAR ANALYSIS
  UNDER CONSTRUCTION
*/
{
  DWORD memory0,memory1,memory2;

  FILE *fin,*fout,*fonl;                           /*FILE 8 BYTES*/
  FILE *fiteration,*fcr,*fresidual;
  FILE *ffstart,*fnew;
  double *ddisp,*dreact,*iform;
  struct memoryelem *melem;
  /*FILE *felem,*fdisp,*freact;*/
  char dir[]=DIRECTORY;                            /*DATA DIRECTORY*/
  char s[80],string[400]/*,fname[256]*/;
  int i,j,ii,jj,iii;
  int nnode,nelem,nsect,nreact;
  int code;
  long int loffset,msize,nline;
  long int time;
  struct gcomponent ginit={0,0,0.0,NULL};
  struct gcomponent *gmtx,*k,*ge,*g,*p;/*GLOBAL MATRIX*/
  double gg,kk;
  double *gvct,*gvct0,*gvct1,*gvct2;                       /*GLOBAL VECTOR*/
  double **drccos,**tmatrix,**estiff/*,*estress*/,*qvct,**tt;
  double **Kd,**S,**Kr,**I;
  double determinant,sign,safety,dsafety;
  double func[2],*tforce1,*tforce2;
  double Ny;

  /**** setting NONLINEAR ****/
  int laps=100;                /*nIStep*/
  int maxiteration=20;         /*nJStep*/
  int firstnlaps=20;           /*firstNStep*/
  double firstsafety=0.0001;   /*firstIncrementRatio*/
  double tolerance=0.000001;   /*epsilon*/
  int nfnodes=3;
  int fnodes[3];               /*node*/
  fnodes[0]=243,fnodes[1]=233,fnodes[2]=196;
  int fdirection=2;            /*direction*/
  /**** setting NONLINEAR ****/

  tforce1=(double *)malloc(nfnodes*sizeof(double));
  tforce2=(double *)malloc(nfnodes*sizeof(double));

  int nlap;
  int iteration;
  double residual;

  clock_t t0,t1,t2;

  struct osect *sects;
  struct onode *nodes,*ninit;
  struct owire elem;
  struct owire *elems;
  struct oconf *confs;

  struct corotational *corot;   /*CR PROP*/

  memory0=availablephysicalmemory("INITIAL:");   /*MEMORY AVAILABLE*/

  fin=fgetstofopen(dir,"r",idinput);              /*OPEN INPUT FILE*/
  if(fin==NULL)
  {
	errormessage("ACCESS IMPOSSIBLE.");
	return 0;
  }
  fout=fgetstofopen("\0","w",ID_OUTPUTFILE);          /*OUTPUT FILE*/
  fonl=fopen("hognon.onl","w");             /*ITELATION OUTPUT FILE*/

  fiteration=fopen("iteration.txt","w");
//  ffstart=fgetstofopenII(dir,"r","fstart.txt");

  t0=clock();                                        /*CLOCK BEGIN.*/

  /*getincrement((wmenu.childs+2)->hwnd,&laps,&dsafety);*/
  dsafety=(1.0-firstnlaps*firstsafety)/(laps-firstnlaps);

  inputinit(fin,&nnode,&nelem,&nsect);             /*INPUT INITIAL.*/
  sprintf(string,"NODES=%d ELEMS=%d SECTS=%d",nnode,nelem,nsect);
  errormessage(string);
  if(fout!=NULL) fprintf(fonl,"%s\n",string);

  msize=6*nnode;                           /*SIZE OF GLOBAL MATRIX.*/

  gmtx =(struct gcomponent *)                 /*DIAGONALS OF GLOBAL MATRIX.*/
		malloc(msize*sizeof(struct gcomponent));
  gvct =(double *)malloc(msize*sizeof(double));            /*GLOBAL VECTOR.*/
  gvct0=(double *)malloc(msize*sizeof(double)); /*GLOBAL VECTOR FOR FSTART.*/
  gvct1=(double *)malloc(msize*sizeof(double));   /*GLOBAL VECTOR ONLY CMQ.*/
  gvct2=(double *)malloc(msize*sizeof(double));       /*TRUE GLOBAL VECTOR.*/

  if(gmtx==NULL || gvct==NULL || gvct1==NULL || gvct2==NULL) return 0;
  for(i=0;i<=msize-1;i++)
  {
	(gmtx+i)->down=NULL;            /*GLOBAL MATRIX INITIALIZATION.*/
	*(gvct+i)=0.0;                  /*GLOBAL VECTOR INITIALIZATION.*/
	*(gvct1+i)=0.0;                 /*GLOBAL VECTOR INITIALIZATION.*/
	*(gvct2+i)=0.0;                 /*GLOBAL VECTOR INITIALIZATION.*/
  }

  I=(double **)malloc(3*sizeof(double *));
  for(i=0;i<3;i++)
  {
	*(I+i)=(double *)malloc(3*sizeof(double));
	for(j=0;j<3;j++)
	{
	  if(i==j) *(*(I+i)+j)=1.0;  /*Identity.*/
	  else     *(*(I+i)+j)=0.0;
	}
  }

  free(af->sects);
  free(af->nodes);
  free(af->ninit);
  free(af->elems);
  free(af->confs);
  free(af->ddisp);
  free(af->melem);

  sects=(struct osect *)malloc(nsect*sizeof(struct osect));
  if(sects==NULL) return 0;
  nodes=(struct onode *)malloc(nnode*sizeof(struct onode));
  if(nodes==NULL) return 0;
  ninit=(struct onode *)malloc(nnode*sizeof(struct onode));
  if(ninit==NULL) return 0;
  elems=(struct owire *)malloc(nelem*sizeof(struct owire));
  if(elems==NULL) return 0;
  confs=(struct oconf *)malloc(msize*sizeof(struct oconf));
  if(confs==NULL) return 0;
  ddisp=(double *)malloc(6*nnode*sizeof(double));
  if(ddisp==NULL) return 0;
  iform=(double *)malloc(6*nnode*sizeof(double));
  if(iform==NULL) return 0;
  melem=(struct memoryelem *)
		malloc(nelem*sizeof(struct memoryelem));
  if(melem==NULL) return 0;

  af->sects=sects;
  af->nodes=nodes;
  af->ninit=ninit;
  af->elems=elems;
  af->confs=confs;
  af->ddisp=ddisp;                     /*DISPLACEMENT:6 DIRECTIONS.*/
  af->melem=melem;

  inputtexttomemory(fin,af);        /*READY TO READ LONG REACTIONS.*/
  nnode=af->nnode;
  nelem=af->nelem;
  nsect=af->nsect;
  nreact=af->nreact;

  initialform(nodes,ddisp,nnode);           /*ASSEMBLAGE FORMATION.*/
  initialelem(elems,melem,nelem);            /*ASSEMBLAGE ELEMENTS.*/

  initialform(nodes,iform,nnode);              /*INITIAL FORMATION.*/

  dreact=(double *)malloc(nreact*sizeof(double));       /*REACTION.*/
  af->dreact=dreact;
  initialreact(fin,dreact,nreact);     /*ASSEMBLAGE LONG REACTIONS.*/

  corot=(struct corotational *)
		 malloc((af->nelem)*sizeof(struct corotational));

  if(MessageBox(NULL,"INITIALIZE CR PROP.","PEZZETTINO",MB_OKCANCEL)==IDOK)
  {
	 initializecorotational(af,corot);
  }
  else
  {
	 MessageBox(NULL,"LOAD CR PROP.","PEZZETTINO",MB_OK);
	 fcr=fgetstofopenII(dir,"r","corotational.txt");  /*OPEN FILE.*/
	 if(fcr!=NULL)
	 {
		corotationaloutputtomemory(fcr,af,corot);
		fclose(fcr);
	 }
	 else
	 {
		MessageBox(NULL,"ACCESS IMPOSSIBLE.","PEZZETTINO",MB_OK);
		MessageBox(NULL,"Could not open corotational.txt.","PEZZETTINO",MB_OK);
		return 0;
	 }
  }
/******FSTART: UNDER CONSTRUCTION******/

  GetAsyncKeyState(VK_LBUTTON);                   /*CLEAR KEY LEFT.*/
  GetAsyncKeyState(VK_RBUTTON);                  /*CLEAR KEY RIGHT.*/

  drawglobalaxis((wdraw.childs+1)->hdcC,(wdraw.childs+1)->vparam,
				 0,0,255);                      /*DRAW GLOBAL AXIS.*/

  nlap=1;
  iteration=1;
  residual=0.0;


  double norm=0.0;
  for(i=0;i<6*nnode;i++)
  {
	norm+=((confs+i)->value)*((confs+i)->value);
  }
  norm=pow(norm,0.5);


  while(nlap<=laps)      /***UJIOKA FOR LOAD INCREMENTAL***/
  {
	fcr=fgetstofopenII(dir,"w","corotational.txt");  /*OPEN FILE.*/
	if(fcr!=NULL)
	{
	  outputcorotational(fcr,af,corot);
	  fclose(fcr);
	}
	af->nlaps=nlap;

	if(iteration==1)
	{
	  sprintf(string,"\nLAP:%d/%d",nlap,laps);
	  errormessage(string);
	  if(fonl!=NULL) fprintf(fonl,"%s\n",string);
	}

	sprintf(string,"\nITERATION:%d",iteration);
	errormessage(string);
	if(fonl!=NULL) fprintf(fonl,"%s\n",string);

	if(nlap<=firstnlaps)  safety=nlap*firstsafety;
	else                  safety=firstnlaps*firstsafety
								 +(nlap-firstnlaps)*dsafety;

	setincrementII((wmenu.childs+2)->hwnd,
				   laps,nlap,dsafety,safety);

	memory1=availablephysicalmemory("REMAIN:");  /*MEMORY AVAILABLE*/

	for(i=1;i<=msize;i++)           /*GLOBAL MATRIX INITIALIZATION.*/
	{
	  g=(gmtx+(i-1))->down;   /*NEXT OF DIAGONAL.*/
	  while(g!=NULL) /*CLEAR ROW.*/
	  {
		p=g;
		g=g->down;
		free(p);
	  }
	  ginit.m=(unsigned short int)i;
	  /*ginit.n=(unsigned short int)i;*/
	  *(gmtx+(i-1))=ginit;
	  *(gvct+(i-1))=0.0;            /*GLOBAL VECTOR INITIALIZATION.*/
	}
	comps=msize; /*INITIAL COMPONENTS=DIAGONALS.*/

	laptime("ASSEMBLING GLOBAL MATRIX.",t0);

	for(i=0;i<nelem;i++)               /*ASSEMBLAGE GLOBAL MATRIX.*/
	{
	  inputelem(elems,melem,i,&elem);        /*READ ELEMENT DATA.*/
	  for(ii=0;ii<2;ii++)
	  {
		for(jj=0;jj<6;jj++)
		{
		  (elems+i)->iconf[ii][jj]=elem.iconf[ii][jj];
		}
	  }
	  inputnode(ddisp,elem.node[0]);                         /*HEAD*/
	  inputnode(ddisp,elem.node[1]);                         /*TAIL*/

	  elem.sect=(elems+i)->sect;               /*READ SECTION DATA.*/

	  if((wdraw.childs+1)->hdcC!=NULL)     /*DRAW DEFORMED ELEMENT.*/
	  {
		drawglobalwire((wdraw.childs+1)->hdcC,
					   (wdraw.childs+1)->vparam,
					   *af,elem,255,255,255,
								255,255,255,0,ONSCREEN/*,i*/);
	  }

	  Kd = calculateDeformationStiffnessMatrixKd_CR(elem);
	  S  = calculateTransformationMatrixS_CR(elem,I);  //to local (R = Identity)
	  Kr = calculateRotationMatrixKr_CR(elem,corot+i); //Symmetric form of Kr (Krenk Book P.125)

	  estiff = assemElementTangentStiffness_CR(elem,corot+i,Kd,S,Kr);  /*ELEMENT TANGENT STIFFNESS.*/  //to global
	  assemgstiffness(gmtx,estiff,&elem); /*ASSEMBLAGE TANGENT STIFFNESS.*/
	  /******CMQ UNAVAILABLE******/

	  for(ii=0;ii<6;ii++) free(*(Kd+ii));
	  free(Kd);
	  for(ii=0;ii<12;ii++) free(*(S+ii));
	  free(S);
	  for(ii=0;ii<12;ii++) free(*(Kr+ii));
	  free(Kr);
	  for(ii=0;ii<12;ii++) free(*(estiff+ii));
	  free(estiff);
	}
	sprintf(string,"GLOBAL MATRIX %ld COMPS ASSEMBLED.",comps);
	laptime(string,t0);
	overlayhdc(*(wdraw.childs+1),SRCPAINT);             /*UPDATE DISPLAY.*/

	/*CALCULATE UNBALANCED FORCE.*/
	/*errormessage("ASSEMBLAGE GLOBAL VECTOR.");*/
	assemconf(confs,gvct,safety,nnode);                  /*GLOBAL VECTOR.*/
	if(fonl!=NULL && /*nlap==1*/iteration==1) fprintf(fonl,"\"TARGET FORCE {F}\"\n");
	for(i=0;i<nnode;i++)
	{
	  if(fonl!=NULL && /*nlap==1*/iteration==1)
	  {
		fprintf(fonl,"NODE:%5ld %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f\n",(nodes+i)->code,
				*(gvct+(6*i+0)),*(gvct+(6*i+1)),*(gvct+(6*i+2)),
				*(gvct+(6*i+3)),*(gvct+(6*i+4)),*(gvct+(6*i+5)));
	  }
	  for(ii=0;ii<nfnodes;ii++)
	  {
		if((nodes+i)->code==fnodes[ii] && fiteration!=NULL)
		{
		  tforce1[ii]=*(gvct+(6*i+fdirection));
		  tforce2[ii]=*(gvct2+(6*i+fdirection));
		}
	  }
	}
	residual=0.0;
	fresidual=fopen("residual.txt","w");
	for(i=0;i<msize;i++)
	{
	  if((confs+i)->iconf==1) *(gvct2+i)=0.0;
	  else
	  {
		*(gvct+i)-=*(gvct2+i);   /*UNBALANCED FORCE.*/
		residual+=*(gvct+i)**(gvct+i);

		sprintf(string,"  %ld   %ld   %.10f",
				   (nodes+int(i/6))->code,i-(int(i/6))*6,
				   *(gvct+i));
		if(fresidual!=NULL) fprintf(fresidual,"%s\n",string);
	  }
	}
	residual=sqrt(residual);
	sprintf(string,"Norm: %.10f",residual);
	if(fresidual!=NULL) fprintf(fresidual,"%s\n",string);
	if(fresidual!=NULL) fclose(fresidual);

	modifygivend(gmtx,gvct,confs,nnode);   /*0:LOAD 1:DISPLACEMENT.*/

	/*SOLVE {dF}=[K]{dU}*/
	laptime("CROUT LU DECOMPOSITION.",t0);
	nline=croutludecompositionnl(gmtx,
								 gvct,confs,
								 6*nnode,
								 &determinant,&sign);        /*[K]{dU}={dF}*/

	sprintf(string,"LOG{DETERMINANT}=%.5E SIGN=%.1f COMPS=%ld",
			determinant,sign,comps);
	errormessage(string);
	if(fonl!=NULL) fprintf(fonl,"%s\n",string);

	sprintf(string,"SAFETY FACTOR=%.5f",safety);
	errormessage(string);
	if(fonl!=NULL) fprintf(fonl,"%s\n",string);

	if (sign < 0.0)
	{
	  for(ii=1;ii<=msize;ii++)
	  {
		gg=0.0;
		gread(gmtx,ii,ii,&gg);

		if(gg<0.0)
		{
		  sprintf(string,"INSTABLE TERMINATION AT NODE %ld.",
				  (nodes+int((ii-1)/6))->code);
		  errormessage(" ");
		  errormessage(string);
		  if(fonl!=NULL) fprintf(fonl,"%s\n",string);
		}
	  }

	  laptime("\0",t0);

	  fclose(fin);
	  fclose(fonl);
	  fclose(fout);
	  fclose(fiteration);

	  gfree(gmtx,nnode);  /*FREE GLOBAL MATRIX.*/
	  free(gvct);
	  free(gvct1);
	  free(gvct2);

	  memory2=availablephysicalmemory("REMAIN:");
	  sprintf(string,"CONSUMPTION:%ld[BYTES]",(memory0-memory2));
	  errormessage(string);

	  return 1;
	}

	/*UPDATE PARAMETERS*/
	laptime("OUTPUT INTO FILE.",t0);
	if(nlap==laps)
	{
	  if(fout!=NULL)
	  {
		fprintf(fout,"\n\n");
		fprintf(fout,"** FORCES OF MEMBER\n\n");
		fprintf(fout,"  NO   KT NODE         N        Q1        Q2");
		fprintf(fout,"        MT        M1        M2\n\n");
	  }
	}

	if(fonl!=NULL) fprintf(fonl,"\"DISPLACEMENT\"\n");
	outputdisp(gvct,fonl,nnode,nodes); /*INCREMENTAL DISPLACEMENT.*/
	if(fonl!=NULL) fprintf(fonl,"\"STRESS\"\n");

	updateform(ddisp,gvct,nnode);              /*FORMATION UPDATE.*/

	fnew=fopen("newcoord.txt","w");
	for(i=0;i<nnode;i++)
	{
	  {
		sprintf(string,"  %ld    %6f     %.6f     %.6f",
				(nodes+i)->code,(nodes+i)->d[0],(nodes+i)->d[1],(nodes+i)->d[2]);
		if(fnew!=NULL) fprintf(fnew,"%s\n",string);
	  }
	}
	if(fnew!=NULL) fclose(fnew);

	for(i=0;i<msize;i++)
	{
	  *(gvct2+i)=0.0;        /*Nodal Element Force in global frame*/
	}

	for(i=0;i<nelem;i++)                   /*STRESS OUTPUT,UPDATE.*/
	{
	  inputelem(elems,melem,i,&elem);

	  inputnode(ddisp,elem.node[0]);
	  inputnode(ddisp,elem.node[1]);

	  elem.sect=(elems+i)->sect;             /*READ SECTION DATA.*/

	  updateElementProperty_CR(&elem,corot+i,gvct);
	  qvct=calculateNodalForce_CR(&elem,corot+i,melem);

	  for(ii=0;ii<6;ii++)                /*CALCULATE NODAL FORCE.*/
	  {
		*(gvct2+((6*elem.node[0]->loff)+ii))+=*(qvct+ii);
		*(gvct2+((6*elem.node[1]->loff)+ii))+=*(qvct+(6+ii));
	  }

	  //if(nlap==laps)
	  {
//		outputstress(elem,estress,fonl,func);
		outputstressnl(elem,NULL,fonl);
	  }
	  /*free(estress);*/
	  free(qvct);
	}

	/***UJIOKA FOR LOAD INCREMENTAL***/
//	if( (residual<tolerance*norm*safety || iteration>=maxiteration
//		 || (residual<tolerance*norm && nlap>laps) /*&& iteration!=1*/))
	if( residual<tolerance*norm || iteration>=maxiteration )
	{
	  nlap++;
	  iteration=0;
#if 0
	  for(i=0;i<nelem;i++)
	  {
		inputelem(elems,melem,i,&elem);

		elem.sect=(elems+i)->sect;             /*READ SECTION DATA.*/

		if ((elem.sect)->code == 501||
			(elem.sect)->code == 502||
			(elem.sect)->code == 503||
			(elem.sect)->code == 504
		)
		{
		  code = (elem.sect)->code + 10;
		  if ((elem.sect)->code == 501)  Ny=-0.035;
//          if ((elem.sect)->code == 501)  Ny=-0.105;
		  if ((elem.sect)->code == 502)  Ny=-0.033;
		  if ((elem.sect)->code == 503)  Ny=-0.057;
		  if ((elem.sect)->code == 504)  Ny=-0.115;
		  if(elem.stress[0][0]<Ny && nlap<=laps)     /*yield*/
		  {
//     		sprintf(str,"ELEM %d :ESEC %d \n",elem.code,elem.sect->code);
//    	    MessageBox(NULL,str,"yield",MB_OK);
			for(j=0;j<af->nsect;j++)
			{
			  if(code==(af->sects+j)->code)
			  {
				elem.sect=af->sects+j;
				break;
			  }
			}
//     	   	sprintf(str,"ELEM %d :ESEC %d \n",elem.code,(elem.sect)->code);


//          errormessage(str);
//        	MessageBox(NULL,str,"changed(compression)",MB_OK);
		  }
		}
		if ((elem.sect)->code == 511||
			(elem.sect)->code == 512||
			(elem.sect)->code == 513||
			(elem.sect)->code == 514
		)
		{
		  code = (elem.sect)->code - 10;
		  if ((elem.sect)->code == 511)  Ny=-0.035;
//          if ((elem.sect)->code == 511)  Ny=-0.105;
		  if ((elem.sect)->code == 512)  Ny=-0.033;
		  if ((elem.sect)->code == 513)  Ny=-0.057;
		  if ((elem.sect)->code == 514)  Ny=-0.115;
		  if(elem.stress[0][0]>Ny || nlap>laps)     /*elastic*/
		  {
//     		sprintf(str,"ELEM %d :ESEC %d \n",elem.code,elem.sect->code);
//    	    MessageBox(NULL,str,"yield",MB_OK);
			for(j=0;j<af->nsect;j++)
			{
			  if(code==(af->sects+j)->code)
			  {
				elem.sect=af->sects+j;
				break;
			  }
			}
//     	   	sprintf(str,"ELEM %d :ESEC %d \n",elem.code,(elem.sect)->code);


//          errormessage(str);
//        	MessageBox(NULL,str,"changed(compression)",MB_OK);
		  }
		}

		(af->elems+i)->sect=elem.sect;    /*section data update*/

	  }
#endif
	  clearwindow(*(wdraw.childs+1));
	  (wdraw.childs+1)->vparam.vflag.ev.deformation=1;
	  drawarclmframe((wdraw.childs+1)->hdcC,
					 (wdraw.childs+1)->vparam,arc,0,ONSCREEN);
	  SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);
	}
	iteration++;
	sprintf(string,"\nITERATION:%d",iteration);
	errormessage(string);
	/***UJIOKA FOR LOAD INCREMENTAL***/

	if(fonl!=NULL)
	{
	  fprintf(fonl,"\"CURRENT FORM\"\n");
	  for(ii=0;ii<nnode;ii++)
	  {
		sprintf(string,"NODE:%5ld {U}=",(nodes+ii)->code);
		for(jj=0;jj<6;jj++)
		{
		  loffset=6*ii+jj;
		  sprintf(s," %16.12f",*(ddisp+loffset));
		  strcat(string,s);
		}
        fprintf(fonl,"%s\n",string);

		for(iii=0;iii<nfnodes;iii++)
		{
		  if((nodes+ii)->code==fnodes[iii] && fiteration!=NULL)
		  {
			if(safety!=0.0)
			{
			  fprintf(fiteration,"LAP: %3d / %3d NODE %3d {Fz}= %18.12f {U}= %16.12f test= %16.12f\n",
					  nlap,laps,(nodes+ii)->code,tforce1[iii],*(ddisp+6*ii+fdirection),
					  residual/tolerance/norm/safety);
			}
			else
			{
			  fprintf(fiteration,"LAP: %3d / %3d NODE %3d {Fz}= %18.12f {U}= %16.12f\n",
					  nlap,laps,(nodes+ii)->code,tforce1[iii],*(ddisp+6*ii+fdirection));
			}
			fprintf(fiteration,"LAP: %3d / %3d NODE %3d {Fz}= %18.12f {U}= %16.12f\n",
					nlap,laps,(nodes+ii)->code,*(gvct2+(6*ii+fdirection)),*(ddisp+6*ii+fdirection));
		  }
		}
	  }
	  fprintf(fonl,"\"REACTION\"\n");
	  outputreaction(gmtx,gvct,nodes,confs,dreact,fonl,nnode);
	}

	if(nlap==laps && fout!=NULL)
	{
	  fprintf(fout,"\n\n");
	  fprintf(fout,"** DISPLACEMENT OF NODE\n\n");
	  fprintf(fout,"  NO          U          V          W");
	  fprintf(fout,"        KSI        ETA      OMEGA\n\n");

	  for(ii=0;ii<nnode;ii++)
	  {
		sprintf(string,"%4d",(nodes+ii)->code);
		for(jj=0;jj<6;jj++)
		{
		  loffset=6*ii+jj;
		  sprintf(s," %10.6f",*(ddisp+loffset)-*(iform+loffset));
		  strcat(string,s);
		}
		fprintf(fout,"%s\n",string);
	  }

	  fprintf(fout,"\n\n");
	  fprintf(fout,"** INITIAL FORM\n\n");
	  fprintf(fout,"  NO          U          V          W");
	  fprintf(fout,"        KSI        ETA      OMEGA\n\n");

	  for(ii=0;ii<nnode;ii++)
	  {
		sprintf(string,"%4d",(nodes+ii)->code);
		for(jj=0;jj<6;jj++)
		{
		  loffset=6*ii+jj;
		  sprintf(s," %10.6f",*(iform+loffset));
		  strcat(string,s);
		}
		fprintf(fout,"%s\n",string);
	  }

	  fprintf(fout,"\n\n");
	  fprintf(fout,"** REACTION\n\n");
	  fprintf(fout,"  NO  DIRECTION              R    NC\n\n");

	  outputreactionnl(gmtx,gvct,nodes,confs,dreact,fout,nnode);/*REACTION.*/
	}

	t1=laptime("\0",t0);

	memory2=availablephysicalmemory(NULL);
	sprintf(string,"CONSUMPTION:%ld[BYTES]",(memory1-memory2));
	errormessage(string);

	errormessage("L:CONTINUE R:ABORT");            /*L=LEFT R=RIGHT*/
	while(!GetAsyncKeyState(VK_LBUTTON))  /*LEFT CLICK TO CONTINUE.*/
	{
	  if(GetAsyncKeyState(VK_RBUTTON))      /*RIGHT CLICK TO ABORT.*/
	  {
		fclose(fin);
		gfree(gmtx,nnode);  /*FREE GLOBAL MATRIX.*/
		free(gvct);
		free(gvct1);
		free(gvct2);
		free(qvct);

		errormessage(" ");
		errormessage("ABORTED.");
		if(fonl!=NULL) fprintf(fonl,"ABORTED.\n");

		fclose(fonl);

		laptime("\0",t0);
		return 1;
	  }
	  t2=clock();
	  time=(t2-t1)/CLK_TCK;
	  if(time>=WAIT) break;               /*CONTINUE AFTER WAITING.*/
	}
  }                                        /*REPEAT UNTIL INSTABLE.*/

  if((wdraw.childs+1)->hdcC!=NULL &&
	 melem!=NULL && ddisp!=NULL)                 /*DRAW LAST FRAME.*/
  {
	for(i=1;i<=nelem;i++)
	{
	  inputelem(elems,melem,i-1,&elem);
	  for(ii=0;ii<=1;ii++) /*COPY HINGE DATA.*/
	  {
		for(jj=0;jj<=5;jj++)
		{
		  (elems+i-1)->iconf[ii][jj]=elem.iconf[ii][jj];
		}
	  }

	  inputnode(ddisp,elem.node[0]);
	  inputnode(ddisp,elem.node[1]);

	  drawglobalwire((wdraw.childs+1)->hdcC,
					 (wdraw.childs+1)->vparam,
					 *af,elem,255,255,255,
							  255,255,255,0,ONSCREEN/*,i*/);
	}
	overlayhdc(*(wdraw.childs+1),SRCPAINT);       /*UPDATE DISPLAY.*/
  }

  fcr=fgetstofopenII(dir,"w","corotational.txt");  /*OPEN FILE.*/
  if(fcr!=NULL)
  {
	outputcorotational(fcr,af,corot);
	fclose(fcr);
  }
  else
  {
	MessageBox(NULL,"ACCESS IMPOSSIBLE.","PEZZETTINO",MB_OK);
	MessageBox(NULL,"Could not open corotational.txt.","PEZZETTINO",MB_OK);
  }

  fclose(fin);
  fclose(fout);
  fclose(fiteration);

  gfree(gmtx,nnode);  /*FREE GLOBAL MATRIX.*/

  af->eigenvec=(double **)malloc(1*sizeof(double *));
  *((af->eigenvec)+0)=gvct;

  errormessage(" ");
  errormessage("COMPLETED.");
  if(fonl!=NULL) fprintf(fonl,"COMPLETED.\n");

  if(fonl!=NULL) fclose(fonl);

  memory2=availablephysicalmemory("REMAIN:");
  sprintf(string,"CONSUMPTION:%ld[BYTES]",(memory0-memory2));
  errormessage(string);
  errormessage(" ");

  return 0;
}/*pezzettino_staticnonlinear*/

double **calculateDeformationStiffnessMatrixKd_CR(struct owire elem)
/*Kd: 6x6 DEFORMATION STIFFNESS MATRIX.*/
{
  int i,j;
  double **Kd;
  double dx,dy,dz,dl;
  double E,poi,A,Ixx,Iyy,J,G;
  double PsiX = 1.0; // shear flexibility parameter
  double PsiY = 1.0; // shear flexibitily parameter

  dx=elem.node[1]->d[0]-elem.node[0]->d[0];
  dy=elem.node[1]->d[1]-elem.node[0]->d[1];
  dz=elem.node[1]->d[2]-elem.node[0]->d[2];
  dl=sqrt(dx*dx+dy*dy+dz*dz);
  E=elem.sect->E;
  poi=elem.sect->poi;
  A=elem.sect->area;
  Ixx=elem.sect->Ixx;
  Iyy=elem.sect->Iyy;
  J=elem.sect->Jzz;
  G=0.5*E/(1.0+poi);

  Kd=(double **)malloc(6*sizeof(double *));
  for(i=0;i<6;i++)
  {
	*(Kd+i)=(double *)malloc(6*sizeof(double));
	for(j=0;j<6;j++)
	{
	  *(*(Kd+i)+j)=0.0;                               /*INITIAL.*/
	}
  }

  // consititutive stiffness
  *(*(Kd+0)+0)=G*J/dl;
  *(*(Kd+1)+1)=E*Ixx/dl;
  *(*(Kd+2)+2)=E*Iyy/dl;
  *(*(Kd+3)+3)=E*A/dl;
  *(*(Kd+4)+4)=3.0 * PsiX * E * Ixx/dl;
  *(*(Kd+5)+5)=3.0 * PsiY * E * Iyy/dl;

  return Kd;
}/*calculateDeformationStiffnessMatrixKd_CR*/

double **calculateTransformationMatrixS_CR(struct owire elem,double **R)
/*S: 12�~6 TRANSFORMATION MATRIX.*/
{
  int i,j;
  double **S;
  double dx,dy,dz,dl;

  dx=elem.node[1]->d[0]-elem.node[0]->d[0];
  dy=elem.node[1]->d[1]-elem.node[0]->d[1];
  dz=elem.node[1]->d[2]-elem.node[0]->d[2];
  dl=sqrt(dx*dx+dy*dy+dz*dz);

  S=(double **)malloc(12*sizeof(double *));
  for(i=0;i<12;i++)
  {
	*(S+i)=(double *)malloc(6*sizeof(double));
	for(j=0;j<6;j++)
	{
	  *(*(S+i)+j)=0.0;                               /*INITIAL.*/
	}
  }

  for(i=0;i<3;i++)
  {
	for(j=0;j<1;j++)
	{
	  *(*(S+i+0)+j+3)=-1.0*(*(*(R+i)+j+0));  //S.block(0,3,3,1) = -R_.col(0);
	}
  }
  for(i=0;i<3;i++)
  {
	for(j=0;j<1;j++)
	{
	  *(*(S+i+0)+j+4)=-2.0*(*(*(R+i)+j+2))/dl;  //S.block(0,4,3,1) = -2.0 * R_.col(2) / calculateLength();
	}
  }
  for(i=0;i<3;i++)
  {
	for(j=0;j<1;j++)
	{
	  *(*(S+i+0)+j+5)=2.0*(*(*(R+i)+j+1))/dl;  //S.block(0,5,3,1) = 2.0 * R_.col(1) / calculateLength();
	}
  }
  for(i=0;i<3;i++)
  {
	for(j=0;j<1;j++)
	{
	  *(*(S+i+3)+j+0)=-1.0*(*(*(R+i)+j+0));  //S.block(3,0,3,1) = -R_.col(0);
	}
  }
  for(i=0;i<3;i++)
  {
	for(j=0;j<1;j++)
	{
	  *(*(S+i+3)+j+1)=-1.0*(*(*(R+i)+j+1));  //S.block(3,1,3,1) = -R_.col(1);
	}
  }
  for(i=0;i<3;i++)
  {
	for(j=0;j<1;j++)
	{
	  *(*(S+i+3)+j+2)=-1.0*(*(*(R+i)+j+2));  //S.block(3,2,3,1) = -R_.col(2);
	}
  }
  for(i=0;i<3;i++)
  {
	for(j=0;j<1;j++)
	{
	  *(*(S+i+3)+j+4)=1.0*(*(*(R+i)+j+1));  //S.block(3,4,3,1) = R_.col(1);
	}
  }
  for(i=0;i<3;i++)
  {
	for(j=0;j<1;j++)
	{
	  *(*(S+i+3)+j+5)=1.0*(*(*(R+i)+j+2));  //S.block(3,5,3,1) = R_.col(2);
	}
  }
  for(i=0;i<3;i++)
  {
	for(j=0;j<1;j++)
	{
	  *(*(S+i+6)+j+3)=1.0*(*(*(R+i)+j+0));  //S.block(6,3,3,1) = R_.col(0);
	}
  }
  for(i=0;i<3;i++)
  {
	for(j=0;j<1;j++)
	{
	  *(*(S+i+6)+j+4)=2.0*(*(*(R+i)+j+2))/dl;  //S.block(6,4,3,1) = 2.0 * R_.col(2) / calculateLength();
	}
  }
  for(i=0;i<3;i++)
  {
	for(j=0;j<1;j++)
	{
	  *(*(S+i+6)+j+5)=-2.0*(*(*(R+i)+j+1))/dl;  //S.block(6,5,3,1) = -2.0 * R_.col(1) / calculateLength();
	}
  }
  for(i=0;i<3;i++)
  {
	for(j=0;j<1;j++)
	{
	  *(*(S+i+9)+j+0)=1.0*(*(*(R+i)+j+0));  //S.block(9,0,3,1) = R_.col(0);
	}
  }
  for(i=0;i<3;i++)
  {
	for(j=0;j<1;j++)
	{
	  *(*(S+i+9)+j+1)=1.0*(*(*(R+i)+j+1));  //S.block(9,1,3,1) = R_.col(1);
	}
  }
  for(i=0;i<3;i++)
  {
	for(j=0;j<1;j++)
	{
	  *(*(S+i+9)+j+2)=1.0*(*(*(R+i)+j+2));  //S.block(9,2,3,1) = R_.col(2);
	}
  }
  for(i=0;i<3;i++)
  {
	for(j=0;j<1;j++)
	{
	  *(*(S+i+9)+j+4)=1.0*(*(*(R+i)+j+1));  //S.block(9,4,3,1) = R_.col(1);
	}
  }
  for(i=0;i<3;i++)
  {
	for(j=0;j<1;j++)
	{
	  *(*(S+i+9)+j+5)=1.0*(*(*(R+i)+j+2));  //S.block(9,5,3,1) = R_.col(2);
	}
  }

  return S;
}/*calculateTransformationMatrixS_CR*/

double **calculateRotationMatrixKr_CR(struct owire elem,struct corotational *corot)
/*Kr: ROTATION STIFFNESS MATRIX.*/
/*cf. Symmetric form of Kr (Krenk Book P.125)*/
{
  int i,j;
  double **Kr,**S,**R;
  double **Kr11,**Kr12,**Kr14,**Kr22,**Kr44,**Kr24;;
  double *qe;
  double dx,dy,dz,dl;
  double N, Qx, Qy, M, mAx, mAy, mBx, mBy; // x,y are y,z respectively in Krenk book

  dx=elem.node[1]->d[0]-elem.node[0]->d[0];
  dy=elem.node[1]->d[1]-elem.node[0]->d[1];
  dz=elem.node[1]->d[2]-elem.node[0]->d[2];
  dl=sqrt(dx*dx+dy*dy+dz*dz);

  Kr=(double **)malloc(12*sizeof(double *));
  for(i=0;i<12;i++)
  {
	*(Kr+i)=(double *)malloc(12*sizeof(double));
	for(j=0;j<12;j++)
	{
	  *(*(Kr+i)+j)=0.0;                               /*INITIAL.*/
	}
  }

  qe=(double *)malloc(12*sizeof(double));

  Kr11=(double **)malloc(3*sizeof(double *));
  Kr12=(double **)malloc(3*sizeof(double *));
  Kr14=(double **)malloc(3*sizeof(double *));
  Kr22=(double **)malloc(3*sizeof(double *));
  Kr44=(double **)malloc(3*sizeof(double *));
  Kr24=(double **)malloc(3*sizeof(double *));
  for(i=0;i<3;i++)
  {
	*(Kr11+i)=(double *)malloc(3*sizeof(double));
	*(Kr12+i)=(double *)malloc(3*sizeof(double));
	*(Kr14+i)=(double *)malloc(3*sizeof(double));
	*(Kr22+i)=(double *)malloc(3*sizeof(double));
	*(Kr44+i)=(double *)malloc(3*sizeof(double));
	*(Kr24+i)=(double *)malloc(3*sizeof(double));
  }

  R=(double **)malloc(3*sizeof(double *));
  for(i=0;i<3;i++)
  {
	*(R+i)=(double *)malloc(3*sizeof(double));
	for(j=0;j<3;j++)
	{
	  if(i==j) *(*(R+i)+j)=1.0;  /*Identity.*/
	  else     *(*(R+i)+j)=0.0;
	}
  }

  S=calculateTransformationMatrixS_CR(elem,R); //to local (R = Identity)

  for(i=0;i<12;i++)
  {
	*(qe+i)=0.0;
	for(j=0;j<6;j++)
	{
	  *(qe+i) += (*(*(S+i)+j))*(corot->force[j]);
	}
  }

  M = corot->force[0];
  N = corot->force[3];
  mAx = qe[4];
  mAy = qe[5];
  mBx = qe[10];
  mBy = qe[11];
  Qx = -2.0 * corot->force[5] / dl;
  Qy = 2.0 * corot->force[4] / dl;

  /*** not including local geometric stiffness ***/
  //	Kr11 << 0.0, -Qx,  -Qy,
  //	        -Qx,   N,  0.0,
  //	        -Qy, 0.0,    N;
  //	Kr11 /= l;
  *(*(Kr11+0)+0)=0.0;  *(*(Kr11+0)+1)=-Qx;  *(*(Kr11+0)+2)=-Qy;
  *(*(Kr11+1)+0)=-Qx;  *(*(Kr11+1)+1)=  N;  *(*(Kr11+1)+2)=0.0;
  *(*(Kr11+2)+0)=-Qy;  *(*(Kr11+2)+1)=0.0;  *(*(Kr11+2)+2)=  N;
  for(i=0;i<3;i++)
  {
	for(j=0;j<3;j++)
	{
	  *(*(Kr11+i)+j) /= dl;
	}
  }
  //	Kr12 << 0.0, 0.0, 0.0,
  //	        mAx,   M, 0.0,
  //	        mAy, 0.0,   M;
  //	Kr12 /= l;
  *(*(Kr12+0)+0)=0.0;  *(*(Kr12+0)+1)=0.0;  *(*(Kr12+0)+2)=0.0;
  *(*(Kr12+1)+0)=mAx;  *(*(Kr12+1)+1)=  M;  *(*(Kr12+1)+2)=0.0;
  *(*(Kr12+2)+0)=mAy;  *(*(Kr12+2)+1)=0.0;  *(*(Kr12+2)+2)=  M;
  for(i=0;i<3;i++)
  {
	for(j=0;j<3;j++)
	{
	  *(*(Kr12+i)+j) /= dl;
	}
  }
  //	Kr14 << 0.0, 0.0, 0.0,
  //	        mBx, -M, 0.0,
  //	        mBy, 0.0, -M;
  //	Kr14 /= l;
  *(*(Kr14+0)+0)=0.0;  *(*(Kr14+0)+1)=0.0;  *(*(Kr14+0)+2)=0.0;
  *(*(Kr14+1)+0)=mBx;  *(*(Kr14+1)+1)= -M;  *(*(Kr14+1)+2)=0.0;
  *(*(Kr14+2)+0)=mBy;  *(*(Kr14+2)+1)=0.0;  *(*(Kr14+2)+2)= -M;
  for(i=0;i<3;i++)
  {
	for(j=0;j<3;j++)
	{
	  *(*(Kr14+i)+j) /= dl;
	}
  }
  //	Kr22 << 0.0, -mAy, mAx,
  //	        -mAy, 0.0, 0.0,
  //			 mAx, 0.0, 0.0;
  //	Kr22 /= 2.0;
  *(*(Kr22+0)+0)=0.0;  *(*(Kr22+0)+1)=-mAy; *(*(Kr22+0)+2)=mAx;
  *(*(Kr22+1)+0)=-mAy; *(*(Kr22+1)+1)=0.0;  *(*(Kr22+1)+2)=0.0;
  *(*(Kr22+2)+0)=mAx;  *(*(Kr22+2)+1)=0.0;  *(*(Kr22+2)+2)=0.0;
  for(i=0;i<3;i++)
  {
	for(j=0;j<3;j++)
	{
	  *(*(Kr22+i)+j) /= 2.0;
	}
  }
  //	Kr44 << 0.0, -mBy, mBx,
  //	        -mBy, 0.0, 0.0,
  //		     mBx, 0.0, 0.0;
  //	Kr44 /= 2.0;
  *(*(Kr44+0)+0)=0.0;  *(*(Kr44+0)+1)=-mBy; *(*(Kr44+0)+2)=mBx;
  *(*(Kr44+1)+0)=-mBy; *(*(Kr44+1)+1)=0.0;  *(*(Kr44+1)+2)=0.0;
  *(*(Kr44+2)+0)=mBx;  *(*(Kr44+2)+1)=0.0;  *(*(Kr44+2)+2)=0.0;
  for(i=0;i<3;i++)
  {
	for(j=0;j<3;j++)
	{
	  *(*(Kr44+i)+j) /= 2.0;
	}
  }
  //	Kr24 << 0.0, 0.0, 0.0,
  //			0.0, 0.0,   M,
  //			0.0,  -M, 0.0;
  //	Kr24 /= 2.0;
  *(*(Kr24+0)+0)=0.0;  *(*(Kr24+0)+1)=0.0;  *(*(Kr24+0)+2)=0.0;
  *(*(Kr24+1)+0)=0.0;  *(*(Kr24+1)+1)=0.0;  *(*(Kr24+1)+2)=  M;
  *(*(Kr24+2)+0)=0.0;  *(*(Kr24+2)+1)= -M;  *(*(Kr24+2)+2)=0.0;
  for(i=0;i<3;i++)
  {
	for(j=0;j<3;j++)
	{
	  *(*(Kr24+i)+j) /= 2.0;
	}
  }


  for(i=0;i<3;i++)
  {
	for(j=0;j<3;j++)
	{
	  *(*(Kr+i+0)+j+0)=*(*(Kr11+i)+j);  //Kr.block(0,0,3,3) = Kr11;
	}
  }
  for(i=0;i<3;i++)
  {
	for(j=0;j<3;j++)
	{
	  *(*(Kr+i+0)+j+3)=*(*(Kr12+i)+j);  //Kr.block(0,3,3,3) = Kr12;
	}
  }
  for(i=0;i<3;i++)
  {
	for(j=0;j<3;j++)
	{
	  *(*(Kr+i+0)+j+6)=-1.0*(*(*(Kr11+i)+j));  //Kr.block(0,6,3,3) = -Kr11;
	}
  }
  for(i=0;i<3;i++)
  {
	for(j=0;j<3;j++)
	{
	  *(*(Kr+i+0)+j+9)=*(*(Kr14+i)+j);  //Kr.block(0,9,3,3) = Kr14;
	}
  }
  for(i=0;i<3;i++)
  {
	for(j=0;j<3;j++)
	{
	  *(*(Kr+i+3)+j+0)=*(*(Kr12+j)+i);  //Kr.block(3,0,3,3) = Kr12.transpose();
	}
  }
  for(i=0;i<3;i++)
  {
	for(j=0;j<3;j++)
	{
	  *(*(Kr+i+3)+j+3)=*(*(Kr22+i)+j);  //Kr.block(3,3,3,3) = Kr22;
	}
  }
  for(i=0;i<3;i++)
  {
	for(j=0;j<3;j++)
	{
	  *(*(Kr+i+3)+j+6)=-1.0*(*(*(Kr12+j)+i));  //Kr.block(3,6,3,3) = -Kr12.transpose();
	}
  }
  for(i=0;i<3;i++)
  {
	for(j=0;j<3;j++)
	{
	  *(*(Kr+i+3)+j+9)=*(*(Kr24+i)+j);  //Kr.block(3,9,3,3) = Kr24;
	}
  }
  for(i=0;i<3;i++)
  {
	for(j=0;j<3;j++)
	{
	  *(*(Kr+i+6)+j+0)=-1.0*(*(*(Kr11+i)+j));  //Kr.block(6,0,3,3) = -Kr11;
	}
  }
  for(i=0;i<3;i++)
  {
	for(j=0;j<3;j++)
	{
	  *(*(Kr+i+6)+j+3)=-1.0*(*(*(Kr12+i)+j));  //Kr.block(6,3,3,3) = -Kr12;
	}
  }
  for(i=0;i<3;i++)
  {
	for(j=0;j<3;j++)
	{
	  *(*(Kr+i+6)+j+6)=*(*(Kr11+i)+j);  //Kr.block(6,6,3,3) = Kr11;
	}
  }
  for(i=0;i<3;i++)
  {
	for(j=0;j<3;j++)
	{
	  *(*(Kr+i+6)+j+9)=-1.0*(*(*(Kr14+i)+j));  //Kr.block(6,9,3,3) = -Kr14;
	}
  }
  for(i=0;i<3;i++)
  {
	for(j=0;j<3;j++)
	{
	  *(*(Kr+i+9)+j+0)=*(*(Kr14+j)+i);  //Kr.block(9,0,3,3) = Kr14.transpose();
	}
  }
  for(i=0;i<3;i++)
  {
	for(j=0;j<3;j++)
	{
	  *(*(Kr+i+9)+j+3)=*(*(Kr24+j)+i);  //Kr.block(9,3,3,3) = Kr24.transpose();
	}
  }
  for(i=0;i<3;i++)
  {
	for(j=0;j<3;j++)
	{
	  *(*(Kr+i+9)+j+6)=-1.0*(*(*(Kr14+j)+i));  //Kr.block(9,6,3,3) = -Kr14.transpose();
	}
  }
  for(i=0;i<3;i++)
  {
	for(j=0;j<3;j++)
	{
	  *(*(Kr+i+9)+j+9)=*(*(Kr44+i)+j);  //Kr.block(9,9,3,3) = Kr44;
	}
  }

  freematrix(S,12); freematrix(R,3);
  freematrix(Kr11,3); freematrix(Kr12,3);
  freematrix(Kr14,3); freematrix(Kr22,3);
  freematrix(Kr44,3); freematrix(Kr24,3);
  free(qe);

  return Kr;
}/*calculateRotationMatrixKr_CR*/

double **assemElementTangentStiffness_CR(struct owire elem,struct corotational *corot,double **Kd,double **S,double **Kr)
/*Ke: ELEMENT TANGENT STIFFNESS MATRIX.*/
{
  int i,j,ii,jj;
  double **Ke;
  double **R,**tmatrix;
  double **St,**SKd;
  double dx,dy,dz,dl;

  R=(double **)malloc(3*sizeof(double *));
  for(i=0;i<3;i++)
  {
	*(R+i)=(double *)malloc(3*sizeof(double));
	for(j=0;j<3;j++)
	{
	  *(*(R+i)+j)=corot->basis[i+3*j];
	}
  }
  tmatrix=transmatrix(R);  /*TRANSFORMATION MATRIX.*/
  tmatrix=matrixtranspose(tmatrix,12);     /*[T]=[Rt]*/

  /*[Ke]=[S][Kd][St]+[Kr]*/ //to local
  //not in block matrix format
  SKd=matrixmatrixIII(S,Kd,12,6,6);   /*[S][Kd]*/
  St=matrixtransposeIII(S,12,6);      /*[St]*/
  Ke=matrixmatrixIII(SKd,St,12,6,12); /*[S][Kd][St]*/
  for(i=0;i<12;i++)
  {
	for(j=0;j<12;j++)
	{
	  *(*(Ke+i)+j) += *(*(Kr+i)+j);   /*[Ke]=[S][Kd][St]+[Kr]*/
	}
  }

  Ke=transformation(Ke,tmatrix);      /*[K]=[Tt][Ke][T]*/  //to global

  freematrix(R,3);
  freematrix(tmatrix,12);
  freematrix(St,6);
  freematrix(SKd,12);

  return Ke;
}/*assemElementTangentStiffness_CR*/

void updateElementProperty_CR(struct owire *elem,struct corotational *corot,double *gvct)
{
  int i,j;
  double dx, dy, dz,dl;
  double norm=0.0;
  double **R,**Rt,**Kd,**I,*gdisp;
  double *modevector,*t;
  double *deltax, *nx;
  double *v3d1, *v3d2;    // VECTOR 3d.
  double **m2d;           // MATRIX 2d.
  double **m3d1, **m3d2;  // MATRIX 3d.
  double **nynz, **nynz0; // MATRIX 3x2.
  double **nv, **nvt, **nnt;  // nv = vector n in matrix form.
  double *dPhiSym,*dPhiAsym;

  dx=elem->node[1]->d[0]-elem->node[0]->d[0];
  dy=elem->node[1]->d[1]-elem->node[0]->d[1];
  dz=elem->node[1]->d[2]-elem->node[0]->d[2];
  dl=sqrt(dx*dx+dy*dy+dz*dz);

  modevector = (double *)malloc(6*sizeof(double));
  deltax = (double *)malloc(3*sizeof(double));
  nx = (double *)malloc(3*sizeof(double));

  for(i=0;i<3;i++)
  {
//	*(nx+i)=corot->basis[3*i];    /*WRONG.*/
	*(nx+i)=corot->basis[i];   /*Debugged.*/
  }
  R=(double **)malloc(3*sizeof(double *));
  for(i=0;i<3;i++)
  {
	*(R+i)=(double *)malloc(3*sizeof(double));
	for(j=0;j<3;j++)
	{
	  *(*(R+i)+j)=corot->basis[i+3*j];
	}
  }
  Rt=matrixtranspose(R,3);  /*MATRIX [Rt].*/

  m2d=(double **)malloc(2*sizeof(double *));
  for(i=0;i<2;i++)
  {
	*(m2d+i)=(double *)malloc(2*sizeof(double));
	for(j=0;j<2;j++)
	{
	  *(*(m2d+i)+j)=0.0;  /*INITIAL.*/
	}
  }
  I=(double **)malloc(3*sizeof(double *));
  for(i=0;i<3;i++)
  {
	*(I+i)=(double *)malloc(3*sizeof(double));
	for(j=0;j<3;j++)
	{
	   if(i==j) *(*(I+i)+j)=1.0;  /*Identity.*/
	   else     *(*(I+i)+j)=0.0;
	}
  }

  gdisp=extractdisplacement(*elem,gvct);

/*calculateRotationIncrement*/
  v3d1=(double *)malloc(3*sizeof(double));
  for(i=0;i<3;i++)
  {
	*(v3d1+i) = (*(gdisp+9+i))-(*(gdisp+3+i));
  }
  dPhiSym = matrixvector(Rt,v3d1,3);
  for(i=0;i<3;i++)
  {
	*(v3d1+i) = (*(gdisp+6+i))-(*(gdisp+i));
  }
  v3d2 = outerproductII(nx,v3d1);
  for(i=0;i<3;i++)
  {
	*(v3d2+i) *= 2.0/dl;
	*(v3d1+i) = (*(gdisp+9+i)) + (*(gdisp+3+i)) - *(v3d2+i);
  }
  dPhiAsym = matrixvector(Rt,v3d1,3);
  for(i=0;i<3;i++)
  {
	corot->mode[i] += dPhiSym[i];
  }
  corot->phiAsymX += dPhiAsym[0];
  corot->mode[3]   = dl - corot->length0;
  corot->mode[4]  += dPhiAsym[1];
  corot->mode[5]  += dPhiAsym[2];
  for(i=0;i<6;i++)
  {
	*(modevector+i)=corot->mode[i];
  }

/*updateElementBasisR*/
// rotate element basis around axis
  *(*(m2d+0)+0) = cos(dPhiAsym[0]);
  *(*(m2d+0)+1) =-sin(dPhiAsym[0]);
  *(*(m2d+1)+0) = sin(dPhiAsym[0]);
  *(*(m2d+1)+1) = cos(dPhiAsym[0]);
  nynz0=(double **)malloc(3*sizeof(double *));
  for(i=0;i<3;i++)
  {
	*(nynz0+i)=(double *)malloc(2*sizeof(double));
	for(j=0;j<2;j++)
	{
	  *(*(nynz0+i)+j)=*(*(R+i)+j+1);
	}
  }
  nynz = matrixmatrixIII(nynz0,m2d,3,2,2);  //R_.rightCols<2>() = R_.rightCols<2>() * m2d
  for(i=0;i<3;i++)
  {
	for(j=0;j<2;j++)
	{
	  *(*(R+i)+j+1)=*(*(nynz+i)+j);
	}
  }
// rotate basis to new element axis and update R
  *(deltax+0) = dx;
  *(deltax+1) = dy;
  *(deltax+2) = dz;
  nv =(double **)malloc(3*sizeof(double *));
  nvt=(double **)malloc(1*sizeof(double *));
  *(nvt+0) =(double *)malloc(3*sizeof(double));
  for(i=0;i<3;i++)
  {
	*(nv+i)  =(double *)malloc(1*sizeof(double));
	*(*(nv+i)+0) = (*(*(R+i)+0)) + (*(deltax+i))/dl;  //{n}={nx}+deltax/dl
  }
  for(i=0;i<3;i++)
  {
	norm += (*(*(nv+i))+0)*(*(*(nv+i)+0));   //norm=|{n}|^2
  }
  norm=sqrt(norm);  //norm=|{n}|
  for(i=0;i<3;i++)
  {
	*(*(nv+i)+0) /= norm;         //{n}={n}/|{n}|
	*(*(nvt+0)+i) = *(*(nv+i)+0); //{nt}={n}T
  }
  m3d1 =(double **)malloc(3*sizeof(double *));
  m3d2 =(double **)malloc(3*sizeof(double *));
  for(i=0;i<3;i++)
  {
	*(m3d1+i)=(double *)malloc(3*sizeof(double));
	*(m3d2+i)=(double *)malloc(3*sizeof(double));
  }
  for(i=0;i<3;i++)
  {
	*(*(m3d1+i)+0)=-1.0*(*(*(R+i)+0));
	*(*(m3d1+i)+1)= 1.0*(*(*(R+i)+1));
	*(*(m3d1+i)+2)= 1.0*(*(*(R+i)+2));      //[m3d1]=[{-nx},{ny},{nz}]
  }
  for(i=0;i<3;i++)
  {
	for(j=0;j<3;j++)
	{
	  if(i==j) *(*(m3d2+i)+j) = 1.0;  /*Identity*/
	  else     *(*(m3d2+i)+j) = 0.0;
	}
  }
  nnt= matrixmatrixIII(nv,nvt,3,1,3);
  for(i=0;i<3;i++)
  {
	for(j=0;j<3;j++)
	{
	  *(*(m3d2+i)+j) -= 2.0*(*(*(nnt+i)+j));    //[m3d2]=[I]-2{n}{n}T
	}
  }
  matrixmatrixII(R,m3d2,m3d1,3);     //R=[m3d2][m3d1]=([I]-2{n}{n}T)[{-nx},{ny},{nz}]

  for(i=0;i<3;i++)
  {
	for(j=0;j<3;j++)
	{
	  corot->basis[i+3*j] = *(*(R+i)+j);
	}
  }

/*calculateInternalElementForce*/
  Kd = calculateDeformationStiffnessMatrixKd_CR(*elem);
  t = matrixvector(Kd,modevector,6); // {t} = [Kd]{v} in case Kd includes only constitutive stiffness
  for(i=0;i<6;i++)
  {
	corot->force[i] = *(t+i);
  }

  freematrix(R,3);
  freematrix(Rt,3);
  freematrix(Kd,6);
  freematrix(I,3);
  free(gdisp);
  free(t);
  free(deltax);
  free(nx);
  free(modevector);
  freematrix(nynz,3); freematrix(nynz0,3); freematrix(m2d,2);
  freematrix(m3d1,3); freematrix(m3d2,3);
  freematrix(nv,3); freematrix(nvt,1); freematrix(nnt,3);
  free(dPhiSym); free(dPhiAsym);

  return;
}/*updateElementProperty_CR*/

double *calculateNodalForce_CR(struct owire *elem,struct corotational *corot,struct memoryelem *melem)
{
  int i,j;
  double **tmatrix;
  double **S,**R,*q,*qe;
  double **I;
  double dx,dy,dz,dl;

  dx=elem->node[1]->d[0]-elem->node[0]->d[0];
  dy=elem->node[1]->d[1]-elem->node[0]->d[1];
  dz=elem->node[1]->d[2]-elem->node[0]->d[2];
  dl=sqrt(dx*dx+dy*dy+dz*dz);

  I=(double **)malloc(3*sizeof(double *));
  for(i=0;i<3;i++)
  {
	*(I+i)=(double *)malloc(3*sizeof(double));
	for(j=0;j<3;j++)
	{
	   if(i==j) *(*(I+i)+j)=1.0;  /*Identity.*/
	   else     *(*(I+i)+j)=0.0;
	}
  }

  qe=(double *)malloc(12*sizeof(double));
  q=(double *)malloc(12*sizeof(double));
  for(i=0;i<12;i++)
  {
	 *(qe+i)=0.0;
	 *(q+i)=0.0;
  }

  R=(double **)malloc(3*sizeof(double *));
  for(i=0;i<3;i++)
  {
	*(R+i)=(double *)malloc(3*sizeof(double));
	for(j=0;j<3;j++)
	{
	  *(*(R+i)+j)=corot->basis[i+3*j];
	}
  }
  tmatrix=transmatrix(R);                /*TRANSFORMATION MATRIX.*/

  S = calculateTransformationMatrixS_CR(*elem,I);  //to local (R = Identity)

  //in local frame
  for(i=0;i<12;i++)
  {
	 for(j=0;j<6;j++)
	 {
	   *(qe+i)+=(*(*(S+i)+j))*(corot->force[j]);  //  qe = S * t
	 }
  }

  elem->stress[0][0] =-(corot->force[3]);          //Nzi
  elem->stress[0][1] = 2.0*(corot->force[5])/dl;   //Qxi
  elem->stress[0][2] =-2.0*(corot->force[4])/dl;   //Qyi
  elem->stress[0][3] =-(corot->force[0]);          //Mzi
  elem->stress[0][4] = qe[4];                      //Mxi
  elem->stress[0][5] = qe[5];                      //Myi
  elem->stress[1][0] = (corot->force[3]);          //Nzj
  elem->stress[1][1] =-2.0*(corot->force[5])/dl;   //Qxj
  elem->stress[1][2] = 2.0*(corot->force[4])/dl;   //Qyj
  elem->stress[1][3] = (corot->force[0]);          //Mzj
  elem->stress[1][4] = qe[10];                     //Mxi
  elem->stress[1][5] = qe[11];                     //Myi

  for(i=0;i<2;i++)                           /*UPDATE STRESS.*/
  {
	for(j=0;j<6;j++)
	{
	  (melem+(elem->loff))->stress[i][j]=elem->stress[i][j];
	}
  }

  //in global frame
  for(i=0;i<12;i++)
  {
	 for(j=0;j<12;j++)
	 {
	   *(q+i)+=(*(*(tmatrix+i)+j))*(*(qe+j));  //  q = Re * S * t

	 }
  }

  free(qe);
  freematrix(R,3);
  freematrix(I,3);
  freematrix(S,12);
  freematrix(tmatrix,12);

  return q;
}/*calculateNodalForce_CR*/

void initializecorotational(struct arclmframe *af,struct corotational *corot)
/*INITIALIZE CR PROP.*/
{
  int i,j,k;
  double dx,dy,dz;
  double **drccos,**R;

  for(i=0;i<(af->nelem);i++)
  {
	dx=(af->elems+i)->node[1]->d[0]-(af->elems+i)->node[0]->d[0];
	dy=(af->elems+i)->node[1]->d[1]-(af->elems+i)->node[0]->d[1];
	dz=(af->elems+i)->node[1]->d[2]-(af->elems+i)->node[0]->d[2];

	drccos=directioncosine((af->elems+i)->node[0]->d[0],
						   (af->elems+i)->node[0]->d[1],
						   (af->elems+i)->node[0]->d[2],
						   (af->elems+i)->node[1]->d[0],
						   (af->elems+i)->node[1]->d[1],
						   (af->elems+i)->node[1]->d[2],
						   (af->elems+i)->cangle);  /*[DRCCOS]*/

	R=matrixtranspose(drccos,3);  /*[R]=[DRCCOS]T*/

	for(j=0;j<6;j++)
	{
	  (corot+i)->mode[j]=0.0;
	}
	(corot+i)->phiAsymX=0.0;
	(corot+i)->length0=sqrt(dx*dx+dy*dy+dz*dz);

	for(j=0;j<6;j++)
	{
	  (corot+i)->force[j]=0.0;
	}

	for(j=0;j<3;j++)
	{
	  for(k=0;k<3;k++)
	  {
		(corot+i)->basis[j+3*k]=*(*(R+j)+k);
	  }
	}
	freematrix(drccos,3);
	freematrix(R,3);
  }
  return;
}/*corotationaloutputtomemory*/

void corotationaloutputtomemory(FILE *ftext,struct arclmframe *af,struct corotational *corot)
/*TRANSLATE CR PROP TEXT INTO MEMORY.*/
{
  char **data,str[256]="\0";
  int i,j,k,n;
  long int ecode,hcode;
  double ddata;

  fseek(ftext,0L,SEEK_SET);


/*  corot=(struct corotational *)
		 malloc((af->nelem)*sizeof(struct corotational));
*/
  while(strncmp(str,"** DEFORMATION",14)) fgets(str,256,ftext);
  for(i=1;i<=3;i++) fgets(str,256,ftext);

  for(i=1;i<=(af->nelem);i++)
  {
	data=fgetsbrk(ftext,&n);
	if(n!=10) return;

	ecode=strtol(*(data+0),NULL,10);
	if(ecode!=(af->elems+i-1)->code) return;

	/*scode=strtol(*(data+1),NULL,10);*/ /*SECTION.*/
	/*hcode=strtol(*(data+2),NULL,10);*/ /*HEAD NODE.*/

	k=2;
	for(j=0;j<6;j++)
	{
	  (corot+i-1)->mode[j]=strtod(*(data+k),NULL);
	  k++;
	}
	  (corot+i-1)->phiAsymX=strtod(*(data+k),NULL);
	  k++;
	  (corot+i-1)->length0=strtod(*(data+k),NULL);

	for(;n>0;n--) free(*(data+n-1));
	free(data);
  }

  while(strncmp(str,"** INTERNAL",11)) fgets(str,256,ftext);
  for(i=1;i<=3;i++) fgets(str,256,ftext);

  for(i=1;i<=(af->nelem);i++)
  {
	data=fgetsbrk(ftext,&n);
	if(n!=8) return;

	/*hcode=strtol(*(data+0),NULL,10);*/
	/*if(hcode!=(af->nodes+i-1)->code) return;*/

	k=2;
	for(j=0;j<6;j++)
	{
	  (corot+i-1)->force[j]=strtod(*(data+k),NULL);
	  k++;
	}

	for(;n>0;n--) free(*(data+n-1));
	free(data);
  }

  while(strncmp(str,"** ELEMENT BASIS",16)) fgets(str,256,ftext);
  for(i=1;i<=3;i++) fgets(str,256,ftext);

  for(i=1;i<=(af->nelem);i++)
  {
	data=fgetsbrk(ftext,&n);
	if(n!=11) return;

	*(af->dreact+i)=ddata;

	k=2;
	for(j=0;j<9;j++)
	{
	  (corot+i-1)->basis[j]=strtod(*(data+k),NULL);
	  k++;
	}

	for(;n>0;n--) free(*(data+n-1));
	free(data);
  }

  return;
}/*corotationaloutputtomemory*/

void outputcorotational(FILE *fout,struct arclmframe *af,struct corotational *corot)
/*OUTPUT CR PROP.*/
{
  char **data,str[256]="\0";
  int i,j,k,n;
  long int ecode,scode;
  double ddata;

  fprintf(fout,"\n\n** DEFORMATION MODE VECTOR AND PhiAsymX LENGTH0\n");
  fprintf(fout,"  NO  SECT      PhiSymX          PhiSymY          PhiSymZ          U         PhiAsymY         PhiAsymZ        PhiAsymX         LENGTH0\n\n");
  fprintf(fout,"\n");
  for(i=1;i<=(af->nelem);i++)
  {
	ecode=(af->elems+i-1)->code;
	scode=((af->elems+i-1)->sect)->code;

	sprintf(str," %d  %d   %.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f\n",
			ecode,scode,(corot+i-1)->mode[0],(corot+i-1)->mode[1],(corot+i-1)->mode[2],
						(corot+i-1)->mode[3],(corot+i-1)->mode[4],(corot+i-1)->mode[5],
						(corot+i-1)->phiAsymX,(corot+i-1)->length0);
	fprintf(fout,"%s",str);
  }

  fprintf(fout,"\n\n** INTERNAL ELEMENT FORCE\n");
  fprintf(fout,"  NO  SECT          Mz          Mxs          Mys         N         Mxa       Mya\n\n");
  fprintf(fout,"\n");
  for(i=1;i<=(af->nelem);i++)
  {
	ecode=(af->elems+i-1)->code;
	scode=((af->elems+i-1)->sect)->code;

    sprintf(str," %d  %d   %.12f %.12f %.12f %.12f %.12f %.12f\n",
            ecode,scode,(corot+i-1)->force[0],(corot+i-1)->force[1],(corot+i-1)->force[2],
                        (corot+i-1)->force[3],(corot+i-1)->force[4],(corot+i-1)->force[5]);
    fprintf(fout,"%s",str);
  }


  fprintf(fout,"\n\n** ELEMENT BASIS OF COROTATIONAL BEAM\n");
  fprintf(fout,"   NO  SECT      NX1      NX2      NX3      NY1      NY2      NY3      NZ1      NZ2      NZ3\n\n");
  fprintf(fout,"\n");
  for(i=1;i<=(af->nelem);i++)
  {
    ecode=(af->elems+i-1)->code;
    scode=((af->elems+i-1)->sect)->code;

    sprintf(str," %d  %d   %.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f\n",
            ecode,scode,(corot+i-1)->basis[0],(corot+i-1)->basis[1],(corot+i-1)->basis[2],
                        (corot+i-1)->basis[3],(corot+i-1)->basis[4],(corot+i-1)->basis[5],
						(corot+i-1)->basis[6],(corot+i-1)->basis[7],(corot+i-1)->basis[8]);
    fprintf(fout,"%s",str);
  }

  return;
}/*outputcorotational*/

int saveloadasfstart(char *fname,struct arclmframe *af)
/*SAVE AS FSTART INPUTFILE.*/
{
  FILE *fin;
  char str[256],dandf[256];
  char dir[]=DIRECTORY;
  /*char dir[]="\0";*/
  int i,j,ii,jj,offset;

  /*"INPUTFILE SPECIFICATION"*/
  /*APPELATION*/
  /*NNODE NELEM NSECT*/
  /*ISECT E POI alpha A Ixx Iyy J QxMAX.....MzMIN*/
  /*INODE X Y Z*/
  /*IELEM ISECT NODEI J COORDANGLE BOUNDARY... LONGSTRESS...*/
  /*INODE CONFINEMENTTYPE... CONFINEMENTVALUE...*/
  /*INODE DIRECTION LONGREACTION*/

  strcpy(dandf,dir);
  /*strcat(dandf,"cansav.inp");*/
  strcat(dandf,fname);
  fin=fopen(dandf,"w");                                /*SAVE FILE.*/
  if(fin==NULL)
  {
	errormessage("ACCESS IMPOSSIBLE.");
	return 0;
  }

  for(i=0;i<(af->nnode);i++) /*CONF VECTOR:CONFINEMENT,VALUE.*/
  {
	fprintf(fin,"%5d",(af->nodes+i)->code);

	for(j=0;j<6;j++)
	{
	  offset=6*i+j;
	  fprintf(fin," %1d",(af->confs+offset)->iconf);
	}
	for(j=0;j<6;j++)
	{
	  offset=6*i+j;
	  fprintf(fin," %11.8f",(af->confs+offset)->value);
	}
	fprintf(fin,"\n");
  }
  /*LONGREACTION:UNDER CONSTRUCTION.*/
  fclose(fin);

  if(globalmessageflag==1)
  {
	sprintf(str,"Saved As %s",dandf);
	MessageBox(NULL,str,"Save",MB_OK);
  }

  return 1;
}/*saveloadasfstart*/

/***FOR PEZZETTINO 003 FILE TYPE (-inl2 style)***/
int saveasarclm2(char *fname,struct arclmframe *af)
/*SAVE AS ARCLM2 INPUTFILE (-inl2 style).*/
{
  FILE *fin;
  char str[256],dandf[256];
  char dir[]=DIRECTORY;
  /*char dir[]="\0";*/
  int i,j,ii,jj,offset;

  /*"INPUTFILE SPECIFICATION"*/
  /*APPELATION*/
  /*NNODE NELEM NSECT*/
  /*ISECT E POI alpha A Ixx Iyy J QxMAX.....MzMIN*/
  /*INODE X Y Z*/
  /*IELEM ISECT NODEI J COORDANGLE BOUNDARY... LONGSTRESS...*/
  /*INODE CONFINEMENTTYPE... CONFINEMENTVALUE...*/
  /*INODE DIRECTION LONGREACTION*/

  strcpy(dandf,dir);
  /*strcat(dandf,"cansav.inp");*/
  strcat(dandf,fname);
  fin=fopen(dandf,"w");                                /*SAVE FILE.*/
  if(fin==NULL)
  {
	errormessage("ACCESS IMPOSSIBLE.");
	return 0;
  }

  fprintf(fin,"%5d %5d %5d %5d\n",af->nnode,0,af->nelem,af->nsect);

  for(i=0;i<(af->nsect);i++)
  {
	fprintf(fin,"%5d %.5E %.5f %.4f %.4f %.8f %.8f %.8f",
			(af->sects+i)->code,
			(af->sects+i)->E,
			(af->sects+i)->poi,
            0.0/*(af->sects+i)->alpha*/, /*���c����*/
			(af->sects+i)->area,
			(af->sects+i)->Ixx,
			(af->sects+i)->Iyy,
			(af->sects+i)->Jzz);
	for(j=0;j<6;j++)
	{
	  fprintf(fin," %12.6f %12.6f",(af->sects+i)->fmax[j],
								 (af->sects+i)->fmin[j]);
	}
	if((af->sects+i)->type!=TYPENULL)
	{
//	  fprintf(fin," %5d",(af->sects+i)->type);
	  fprintf(fin," %5d",7);     //temporary
	}
    else
    {
      fprintf(fin," %5d",0);
    }

/*    fprintf(fin," %5d",(af->sects+i)->ocode);    */

    fprintf(fin,"\n");
  }
  for(i=0;i<(af->nnode);i++)
  {
//	fprintf(fin,"%5d %7.3f %7.3f %7.3f\n",
	fprintf(fin,"%5d %7.6f %7.6f %7.6f\n",
			(af->nodes+i)->code,
			(af->nodes+i)->d[GX],
			(af->nodes+i)->d[GY],
			(af->nodes+i)->d[GZ]);
  }
  for(i=0;i<(af->nelem);i++)
  {
    fprintf(fin,"%5d  %5d  %5d %5d  %.5f",
            (af->elems+i)->code,
            (af->elems+i)->sect->code,
            (af->elems+i)->node[0]->code,
            (af->elems+i)->node[1]->code,
            (af->elems+i)->cangle);

    for(ii=0;ii<2;ii++)                  /*BOUNDARY.0:RIGID 1:HINGE*/
    {
      fprintf(fin," %1d %1d %1d",
              (af->elems+i)->iconf[ii][3],
              (af->elems+i)->iconf[ii][4],
              (af->elems+i)->iconf[ii][5]);
    }
    for(ii=0;ii<2;ii++)                                   /*STRESS.*/
    {
      for(jj=0;jj<6;jj++)
      {
		fprintf(fin," %8.5f",(af->elems+i)->stress[ii][jj]);
      }
    }
    fprintf(fin,"\n");
  }

  for(i=0;i<(af->nnode);i++) /*CONF VECTOR:CONFINEMENT,VALUE.*/
  {
    fprintf(fin,"%5d",(af->nodes+i)->code);

    for(j=0;j<6;j++)
    {
      offset=6*i+j;
      fprintf(fin," %1d",(af->confs+offset)->iconf);
    }
    for(j=0;j<6;j++)
    {
      offset=6*i+j;
	  fprintf(fin," %11.8f",(af->confs+offset)->value);
	}
    fprintf(fin,"\n");
  }
  /*LONGREACTION:UNDER CONSTRUCTION.*/
  fclose(fin);

  if(globalmessageflag==1)
  {
    sprintf(str,"Saved As %s",dandf);
    MessageBox(NULL,str,"Save",MB_OK);
  }

  return 1;
}/*saveasarclm2*/

void inputtexttomemory2(FILE *ftext,struct arclmframe *af)
/*TRANSLATE ARCLM INPUTFILE TEXT INTO MEMORY. (-inl2 style)*/
{
  char **data;
  int i,j,ii,jj,k,n;
  long int offset;
  long int scode,hcode,tcode;

  fseek(ftext,0L,SEEK_SET);
  inputinit2(ftext,&(af->nnode),&(af->nelem),&(af->nsect));

  /*"INPUTFILE SPECIFICATION"*/
  /*APPELATION*/
  /*NNODE NELEM NSECT*/
  /*ISECT E POI alpha A Ixx Iyy J QxMAX.....MzMIN*/
  /*INODE X Y Z*/
  /*IELEM ISECT NODEI J COORDANGLE BOUNDARY... LONGSTRESS...*/
  /*INODE CONFINEMENTTYPE... CONFINEMENTVALUE...*/
  /*INODE DIRECTION LONGREACTION*/

  for(i=0;i<(af->nsect);i++)
  {
    (af->sects+i)->loff=i;
    readsect(ftext,(af->sects+i));

    if((af->sects+i)->area==0.0 &&
       (af->sects+i)->Ixx ==0.0 &&
       (af->sects+i)->Iyy ==0.0 &&
       (af->sects+i)->Jzz ==0.0)
    {
      (af->sects+i)->dflag=0;
    }
    else (af->sects+i)->dflag=1;

    (af->sects+i)->hiju[0]=0.0;
    (af->sects+i)->hiju[1]=0.0;
    (af->sects+i)->hiju[2]=0.0;
	(af->sects+i)->lload[0]=0.0;
    (af->sects+i)->lload[1]=0.0;
    (af->sects+i)->lload[2]=0.0;
	(af->sects+i)->perpl[0]=0.0;
    (af->sects+i)->perpl[1]=0.0;
    (af->sects+i)->perpl[2]=0.0;
    (af->sects+i)->dcolor.r=255;
    (af->sects+i)->dcolor.g=255;
    (af->sects+i)->dcolor.b=255;
    (af->sects+i)->ppc.npcurve=0;
  }
  for(i=0;i<(af->nnode);i++)
  {
    (af->nodes+i)->loff=i;

    data=fgetsbrk(ftext,&n);

    (af->nodes+i)->code=strtol(*(data+0),NULL,10);
	(af->nodes+i)->d[0]=strtod(*(data+1),NULL);
    (af->nodes+i)->d[1]=strtod(*(data+2),NULL);
    (af->nodes+i)->d[2]=strtod(*(data+3),NULL);

	*(af->ninit+i)=*(af->nodes+i);

    for(;n>0;n--) free(*(data+n-1));
    free(data);
  }
  for(i=1;i<=(af->nelem);i++)
  {
	(af->elems+i-1)->loff=i-1;

    data=fgetsbrk(ftext,&n);

    (af->elems+i-1)->code=strtol(*(data+0),NULL,10);
    scode=strtol(*(data+1),NULL,10); /*SECTION.*/
    hcode=strtol(*(data+2),NULL,10); /*HEAD NODE.*/
    tcode=strtol(*(data+3),NULL,10); /*TAIL NODE.*/
    (af->elems+i-1)->cangle=strtod(*(data+4),NULL);

    k=5;
    for(ii=0;ii<=1;ii++)                 /*BOUNDARY.0:RIGID 1:HINGE*/
    {
      for(jj=0;jj<=2;jj++)
      {
        (af->elems+i-1)->iconf[ii][jj]=0;
      }
      for(jj=3;jj<=5;jj++)
      {
        (af->elems+i-1)->iconf[ii][jj]
        =(signed char)strtol(*(data+k),NULL,10);
        k++;
      }
    }
    for(ii=0;ii<=1;ii++)                                  /*STRESS.*/
    {
      for(jj=0;jj<=5;jj++)
      {
        (af->elems+i-1)->stress[ii][jj]
        =strtod(*(data+k),NULL);
        k++;
      }
     }

    for(;n>0;n--) free(*(data+n-1));
    free(data);

    offset=0;
    for(ii=0;ii<1;)                                      /*SECTION.*/
    {
      if((af->sects+offset)->code==scode)
      {
        (af->elems+i-1)->sect=af->sects+offset;
        ii++;
      }
      offset++;
    }
    offset=0;
    for(ii=0;ii<2;)                                        /*NODES.*/
    {
      if((af->nodes+offset)->code==hcode)
      {
        (af->elems+i-1)->node[0]=af->nodes+offset;
        ii++;
      }
      if((af->nodes+offset)->code==tcode)
      {
        (af->elems+i-1)->node[1]=af->nodes+offset;
        ii++;
      }
      offset++;
    }
  }

  af->nreact=0;
  for(i=1;i<=(af->nnode);i++) /*CONF VECTOR:CONFINEMENT,VALUE.*/
  {
	data=fgetsbrk(ftext,&n);

    for(j=1;j<=6;j++)
    {
      offset=6*(i-1)+(j-1);

      (af->confs+offset)->iconf
	  =(signed char)strtol(*(data+j),NULL,10);
      (af->confs+offset)->value
      =strtod(*(data+j+6),NULL);

      if((af->confs+offset)->iconf==1) (af->nreact)++;
    }

    for(;n>0;n--) free(*(data+n-1));
    free(data);
  }
  return;
}/*inputtexttomemory2*/

void inputinit2(FILE *fin,int *nnode,int *nelem,int *nsect)
{
  char **data;
  int n;

  data=fgetsbrk(fin,&n);
  *nnode=strtol(*(data+0),NULL,10);
  *nelem=strtol(*(data+2),NULL,10);
  *nsect=strtol(*(data+3),NULL,10);
  for(;n>0;n--) free(*(data+n-1));
  free(data);

  return;
}/*inputinit2*/

void readsect2(FILE *fin,struct osect *sect)
/*READ SECTION FROM INPUT FILE.*/
{
  char **data;
  int i,k,n;

  data=fgetsbrk(fin,&n);
  sect->code=strtol(*(data+0),NULL,10);             /*SECTION CODE.*/
  sect->E=strtod(*(data+1),NULL);                /*YOUNG'S MODULUS.*/
  sect->poi=strtod(*(data+2),NULL);              /*POISSON'S RATIO.*/
//  sect->alpha=strtod(*(data+3),NULL);
  sect->area=strtod(*(data+4),NULL);                        /*AREA.*/
  sect->Ixx=strtod(*(data+5),NULL);
  sect->Iyy=strtod(*(data+6),NULL);
  sect->Jzz=strtod(*(data+7),NULL);          /*ST.VENANT'S TORTION.*/
  k=7;
  for(i=0;i<=5;i++)           /*UPPER,LOWER LIMIT OF YIELD SURFACE.*/
  {
	sect->fmax[i]=strtod(*(data+k),NULL); k++;
	sect->fmin[i]=strtod(*(data+k),NULL); k++;
  }
  if(n>19) sect->type=strtol(*(data+19),NULL,10);   /*SECTION TYPE.*/
  else     sect->type=TYPENULL;

  for(;n>0;n--) free(*(data+n-1));
  free(data);

  return;
}/*readsect2*/

/*UPDATE ORGAN FRAME FOR KIRIGAMI*/ /*updateorganization*/
int updateorganization(struct organ *org)
{
  char str[256],non[80],dir[]=DIRECTORY;
  int i,j,k,n,ii,jj,code;
  int nelem;
  struct oelem *einit;
  FILE *fin,*fout,*frat;

  /*FOR CR PROP UPDATE*/
  FILE *ftxt;
  struct corotational *corot;
  double dl,dx,dy,dz;

  /*TURN OFF MESSAGES*/
  globalmessageflag=0;
  globaldrawflag=0;

  arc  = arci;
  arcx = arci;
  arcy = arci;
  /*free((wdraw.childs+1)->org.loads);*/

  for(i=0;i<org->nelem;i++)
  {
	for(j=0;j<2;j++)
	{
	  for(k=0;k<6;k++) (org->elems+i)->initial[j][k]=0.0;
	}
  } /*INITIAL CMQ UNAVAILABLE.*/

  einit=(struct oelem *)malloc(org->nelem*sizeof(struct oelem));
  for(k=0; k<org->nelem; k++)  *(einit+k)=*(org->elems+k);

  /*EXTRACT ARCLM*/
  extractarclmfromorgan(org,&arc,&arcx,&arcy);

  (wmenu.childs+2)->vparam.vflag.mv.ftype=F_ARCLM;

  /*SAVE AS ARCLM*/
//  saveasarclm((wdraw.childs+1)->inpfilez,&arc);
  saveloadasfstart("fstart.txt",&arc);

  getviewparam((wmenu.childs+2)->hwnd,
			   &((wdraw.childs+1)->vparam));
  (wdraw.childs+1)->vparam.vflag.ev.deformation=0;

  /*OPEN OTL-FILE*/
  fout=fgetstofopen("\0","r",ID_OUTPUTFILEZ);  /*OPEN FILE.*/
  frameoutputtomemory(fout,&arc);
  fclose(fout);

  /*OPEN RAT-FILE*/
  frat=fgetstofopen("\0","r",ID_OUTPUTFILE);
  if(frat==NULL) return -1;
  readsrcanrate(frat,&arc);
  fclose(frat);

/*
  clearwindow(*(wdraw.childs+1));
  drawarclmframe((wdraw.childs+1)->hdcC,
				 (wdraw.childs+1)->vparam,arc,0,ONSCREEN);
  SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);
*/
  //SetDlgItemText(hdwnd,IDV_MODENUM,"1");

  /*LOAD CR PROP*/
  ftxt=fgetstofopenII(dir,"r","corotational.txt");  /*OPEN FILE.*/
  if(ftxt==NULL) return -1;
  corot=(struct corotational *)
		 malloc((arc.nelem)*sizeof(struct corotational));

  corotationaloutputtomemory(ftxt,&arc,corot);
  fclose(ftxt);

  (wdraw.childs+1)->lstatus=ROTATE;
  (wmenu.childs+2)->vparam.vflag.mv.ftype=F_ORGAN;

  /*Initialize Organ Frame and Reopen the inp-file*/ /*20201202 ujok*/
  freeorganization(&((wdraw.childs+1)->org));
  fin=fopen((wdraw.childs+1)->inpfile,"r");

  sprintf(str,"OPENED=%s",(wdraw.childs+1)->inpfile);
  errormessage(str);

  inputorganization(fin,&((wdraw.childs+1)->org),
						&((wdraw.childs+1)->vparam));

  fclose(fin);

  for(i=0;i<org->nnode;i++)
  {
	for(j=0;j<3;j++)
	{
	  (org->nodes+i)->d[j]=*(arc.ddisp+6*i+j);
	}
  }
  /*CHANGE SECTION OF ORGAN FRAME*/ /*20201207 ujok*/
  jj=0;      /*number of walls and slabs*/
  for(i=0;i<org->nelem;i++)
  {
	if((org->elems+i)->type==WALL||(org->elems+i)->type==SLAB)
	{
	  jj = jj+1;;
	}
	ii=i-jj; /*index of arclm frame*/

	dx=(org->elems+i)->nods[1]->d[0]-(org->elems+i)->nods[0]->d[0];
	dy=(org->elems+i)->nods[1]->d[1]-(org->elems+i)->nods[0]->d[1];
	dz=(org->elems+i)->nods[1]->d[2]-(org->elems+i)->nods[0]->d[2];
	dl=sqrt(dx*dx+dy*dy+dz*dz);
  #if 1
	if (((org->elems+i)->sect)->code == 501||
		((org->elems+i)->sect)->code == 502||
		((org->elems+i)->sect)->code == 503||
		((org->elems+i)->sect)->code == 504
	   )
	{
		code = ((org->elems+i)->sect)->code + 10;
		if((arc.elems+ii)->srate[0]>1.0)          /*yielded*/
		{
			sprintf(str,"ELEM %d :ESEC %d \n",(org->elems+i)->code,((org->elems+i)->sect)->code);
			for(j=0;j<((wdraw.childs+1)->org.nsect);j++)
			{
			  if(code==((wdraw.childs+1)->org.sects+j)->code)
			  {
				/*UPDATE CR PROP*/
				(corot+ii)->mode[3]/=(
									 ((((wdraw.childs+1)->org.sects+j)->figs+0)->prop->E)
									 /(((org->elems+i)->sect->figs+0)->prop->E)
									 );
				(corot+ii)->length0 =dl-
									((dl-(corot+ii)->length0)*
									 (((org->elems+i)->sect->figs+0)->prop->E)
									 /((((wdraw.childs+1)->org.sects+j)->figs+0)->prop->E)
									);
				(org->elems+i)->sect=(wdraw.childs+1)->org.sects+j;
				break;
			  }
			}
			sprintf(str,"ELEM %d :ESEC %d \n",(org->elems+i)->code,((org->elems+i)->sect)->code);
			errormessage(str);
		}
	}
  #endif
  #if 0 /*���חp*/
	if (((org->elems+i)->sect)->code == 511||
		((org->elems+i)->sect)->code == 512||
		((org->elems+i)->sect)->code == 513||
		((org->elems+i)->sect)->code == 514
	   )
	{
		code = ((org->elems+i)->sect)->code - 10;
		/*if((arc.elems+ii)->srate[0]<1.0)*/
		{
			sprintf(str,"ELEM %d :ESEC %d \n",(org->elems+i)->code,((org->elems+i)->sect)->code);
			for(j=0;j<((wdraw.childs+1)->org.nsect);j++)
			{
			  if(code==((wdraw.childs+1)->org.sects+j)->code)
			  {
				/*UPDATE CR PROP*/
				(corot+ii)->mode[3]/=(
									 ((((wdraw.childs+1)->org.sects+j)->figs+0)->prop->E)
									 /(((org->elems+i)->sect->figs+0)->prop->E)
									 );
				(corot+ii)->length0 =dl-
									((dl-(corot+ii)->length0)*
									 (((org->elems+i)->sect->figs+0)->prop->E)
									 /((((wdraw.childs+1)->org.sects+j)->figs+0)->prop->E)
									);
				(org->elems+i)->sect=(wdraw.childs+1)->org.sects+j;
				break;
			  }
			}
			sprintf(str,"ELEM %d :ESEC %d \n",(org->elems+i)->code,((org->elems+i)->sect)->code);
			errormessage(str);
		}
	}
  #endif
  }

  fin=fopen(FILEFORCHANGESECT,"w");
  saveorganization(fout, &((wdraw.childs + 1)->org),&((wdraw.childs + 1)->vparam));
  fclose(fin);

  ftxt=fgetstofopenII(dir,"w","corotational_new.txt");  /*OPEN FILE.*/
  if(ftxt==NULL) return -1;
  outputcorotational(ftxt,&arc,corot);
  fclose(ftxt);

  sprintf(str,"Saved as %s \n",FILEFORCHANGESECT);
  MessageBox(NULL,str,"change section automatic",MB_OK);

  globalmessageflag=1;
  globaldrawflag=1;

return 0;
}/*updateorganization*/

int arclm202(struct arclmframe *af,int idinput)
/*
  Non-Linear Shape Analysis for KIRIGAMI
*/
{
  DWORD memory0,memory1,memory2;

  FILE *fin,*fout,*fonl;                           /*FILE 8 BYTES*/
  FILE *fiteration;
  double *ddisp,*dreact,*iform;
  struct memoryelem *melem;
  /*FILE *felem,*fdisp,*freact;*/
  char dir[]=DIRECTORY;                            /*DATA DIRECTORY*/
  char s[80],string[400],str[256]/*,fname[256]*/;
  int i,ii,iii,j,jj;
  int nnode,nelem,nsect,nreact;
  int code;
  long int loffset,msize,nline;
  long int time;
  struct gcomponent ginit={0,0,0.0,NULL};
  struct gcomponent /**kmtx,*/*gemtx,*gmtx,*k,*ge,*g,*p;/*GLOBAL MATRIX*/
  double gg,kk;
  double *gvct,*gvct1,*gvct2;                       /*GLOBAL VECTOR*/
  double **drccos,**tmatrix,**estiff,*estress,*estress2,**tt;
  double determinant,sign,safety,dsafety;
  double func[2],*tforce1,*tforce2;
  double Ny;

  /**** setting NONLINEAR ****/
  int laps=300;                /*nIStep*/
  int maxiteration=20;         /*nJStep*/
  int firstnlaps=100;           /*firstNStep*/
  double firstsafety=0.00002;  /*firstIncrementRatio*/
  double tolerance=0.05;    /*epsilon*/
  int nfnodes=3;
  int fnodes[3];               /*node*/
  fnodes[0]=243,fnodes[1]=233,fnodes[2]=196;
//  int fdirection=2;            /*direction*/
  /**** setting NONLINEAR ****/

  tforce1=(double *)malloc(nfnodes*sizeof(double));
  tforce2=(double *)malloc(nfnodes*sizeof(double));

  int nlap;
  int iteration;
  double residual;

  clock_t t0,t1,t2;

  struct osect *sects;
  struct onode *nodes,*ninit;
  struct owire elem;
  struct owire *elems;
  struct oconf *confs;

  memory0=availablephysicalmemory("INITIAL:");   /*MEMORY AVAILABLE*/

  fin=fgetstofopen(dir,"r",idinput);              /*OPEN INPUT FILE*/
  if(fin==NULL)
  {
    errormessage("ACCESS IMPOSSIBLE.");
    return 0;
  }
  fout=fgetstofopen("\0","w",ID_OUTPUTFILE);          /*OUTPUT FILE*/
  fonl=fopen("hognon.onl","w");             /*ITELATION OUTPUT FILE*/

  fiteration=fopen("iteration.txt","w");

  t0=clock();                                        /*CLOCK BEGIN.*/

  /*getincrement((wmenu.childs+2)->hwnd,&laps,&dsafety);*/
  dsafety=(1.0-firstnlaps*firstsafety)/(laps-firstnlaps);

  inputinit(fin,&nnode,&nelem,&nsect);             /*INPUT INITIAL.*/
  sprintf(string,"NODES=%d ELEMS=%d SECTS=%d",nnode,nelem,nsect);
  errormessage(string);
//  if(fout!=NULL) fprintf(fonl,"%s\n",string);

  msize=6*nnode;                           /*SIZE OF GLOBAL MATRIX.*/

  gemtx=(struct gcomponent *)      /*DIAGONALS OF GEOMETRIC MATRIX.*/
        malloc(msize*sizeof(struct gcomponent));
  gmtx=(struct gcomponent *)          /*DIAGONALS OF GLOBAL MATRIX.*/
        malloc(msize*sizeof(struct gcomponent));
  gvct=(double *)malloc(msize*sizeof(double));           /*GLOBAL VECTOR.*/
  gvct1=(double *)malloc(msize*sizeof(double)); /*GLOBAL VECTOR ONLY CMQ.*/
  gvct2=(double *)malloc(msize*sizeof(double));     /*TRUE GLOBAL VECTOR.*/

  if(gemtx==NULL || gmtx==NULL || gvct==NULL || gvct1==NULL || gvct2==NULL) return 0;
  for(i=0;i<=msize-1;i++)
  {
    (gemtx+i)->down=NULL;        /*GEOMETRIC MATRIX INITIALIZATION.*/
    (gmtx+i)->down=NULL;            /*GLOBAL MATRIX INITIALIZATION.*/
    *(gvct+i)=0.0;                  /*GLOBAL VECTOR INITIALIZATION.*/
    *(gvct1+i)=0.0;                 /*GLOBAL VECTOR INITIALIZATION.*/
    *(gvct2+i)=0.0;                 /*GLOBAL VECTOR INITIALIZATION.*/
  }

  free(af->sects);
  free(af->nodes);
  free(af->ninit);
  free(af->elems);
  free(af->confs);
  free(af->ddisp);
  free(af->melem);

  sects=(struct osect *)malloc(nsect*sizeof(struct osect));
  if(sects==NULL) return 0;
  nodes=(struct onode *)malloc(nnode*sizeof(struct onode));
  if(nodes==NULL) return 0;
  ninit=(struct onode *)malloc(nnode*sizeof(struct onode));
  if(ninit==NULL) return 0;
  elems=(struct owire *)malloc(nelem*sizeof(struct owire));
  if(elems==NULL) return 0;
  confs=(struct oconf *)malloc(msize*sizeof(struct oconf));
  if(confs==NULL) return 0;
  ddisp=(double *)malloc(6*nnode*sizeof(double));
  if(ddisp==NULL) return 0;
  iform=(double *)malloc(6*nnode*sizeof(double));
  if(iform==NULL) return 0;
  melem=(struct memoryelem *)
        malloc(nelem*sizeof(struct memoryelem));
  if(melem==NULL) return 0;

  af->sects=sects;
  af->nodes=nodes;
  af->ninit=ninit;
  af->elems=elems;
  af->confs=confs;
  af->ddisp=ddisp;                     /*DISPLACEMENT:6 DIRECTIONS.*/
  af->melem=melem;

  inputtexttomemory(fin,af);        /*READY TO READ LONG REACTIONS.*/
  nnode=af->nnode;
  nelem=af->nelem;
  nsect=af->nsect;
  nreact=af->nreact;

  initialform(nodes,ddisp,nnode);           /*ASSEMBLAGE FORMATION.*/
  initialelem(elems,melem,nelem);            /*ASSEMBLAGE ELEMENTS.*/

  initialform(nodes,iform,nnode);              /*INITIAL FORMATION.*/

  dreact=(double *)malloc(nreact*sizeof(double));       /*REACTION.*/
  af->dreact=dreact;
  initialreact(fin,dreact,nreact);     /*ASSEMBLAGE LONG REACTIONS.*/

  GetAsyncKeyState(VK_LBUTTON);                   /*CLEAR KEY LEFT.*/
  GetAsyncKeyState(VK_RBUTTON);                  /*CLEAR KEY RIGHT.*/

  drawglobalaxis((wdraw.childs+1)->hdcC,(wdraw.childs+1)->vparam,
                 0,0,255);                      /*DRAW GLOBAL AXIS.*/

  nlap=1;
  iteration=1;
  residual=0.0;
  double norm=0.0;
  for(i=0;i<6*nnode;i++)
  {
    norm+=((confs+i)->value)*((confs+i)->value);
  }
  norm=pow(norm,0.5);

//  for(nlap=1;nlap<=laps;nlap++)
//  while(nlap<=laps)      /***UJIOKA FOR LOAD INCREMENTAL***/
  while(nlap<=laps*1.1)      /***UJIOKA FOR LOAD INCREMENTAL***/
  {
    /***UJIOKA FOR LOAD INCREMENTAL***/
    af->nlaps=nlap;
    /*af->nlaps=1;*/
//    sprintf(string,"\nRESIDUAL:%.5f",residual);
//    errormessage(string);

    if(iteration==1)
    {
      sprintf(string,"\nLAP:%d/%d",nlap,laps);
      errormessage(string);
      if(fonl!=NULL) fprintf(fonl,"%s\n",string);
    }

    sprintf(string,"\nITERATION:%d",iteration);
    errormessage(string);
    if(fonl!=NULL) fprintf(fonl,"%s\n",string);

    if(nlap<=firstnlaps)  safety=nlap*firstsafety;
	else if(nlap<=laps)   safety=firstnlaps*firstsafety+(nlap-firstnlaps)*dsafety;
	else                  safety=1.0-10.0*double((nlap-laps))/double(laps);
//	else if(nlap<=laps+5)   safety=1.0-50.0*double((nlap-laps))/double(laps);
//	else                  safety=0.1*(1-double((nlap-laps-5))/50.0);

//	if(nlap==120) maxiteration=0;

	setincrementII((wmenu.childs+2)->hwnd,
                   laps,nlap,dsafety,safety);

    memory1=availablephysicalmemory("REMAIN:");  /*MEMORY AVAILABLE*/

    for(i=1;i<=msize;i++)           /*GLOBAL MATRIX INITIALIZATION.*/
    {
      g=(gmtx+(i-1))->down;   /*NEXT OF DIAGONAL.*/

      while(g!=NULL) /*CLEAR ROW.*/
      {
        p=g;
        g=g->down;
        free(p);
      }

      ginit.m=(unsigned short int)i;
      /*ginit.n=(unsigned short int)i;*/

      *(gmtx+(i-1))=ginit;

      *(gvct+(i-1))=0.0;            /*GLOBAL VECTOR INITIALIZATION.*/
    }
    comps=msize; /*INITIAL COMPONENTS=DIAGONALS.*/

    laptime("ASSEMBLING GLOBAL MATRIX.",t0);

    for(i=1;i<=nelem;i++)               /*ASSEMBLAGE GLOBAL MATRIX.*/
    {
      inputelem(elems,melem,i-1,&elem);        /*READ ELEMENT DATA.*/
      for(ii=0;ii<2;ii++)
      {
        for(jj=0;jj<6;jj++)
        {
          (elems+i-1)->iconf[ii][jj]=elem.iconf[ii][jj];
        }
      }
      inputnode(ddisp,elem.node[0]);                         /*HEAD*/
      inputnode(ddisp,elem.node[1]);                         /*TAIL*/

      elem.sect=(elems+i-1)->sect;             /*READ SECTION DATA.*/

      if((wdraw.childs+1)->hdcC!=NULL)     /*DRAW DEFORMED ELEMENT.*/
      {
		drawglobalwire((wdraw.childs+1)->hdcC,
                       (wdraw.childs+1)->vparam,
                       *af,elem,255,255,255,
                                255,255,255,0,ONSCREEN/*,i*/);
      }

	  drccos=directioncosine(elem.node[0]->d[0],
							 elem.node[0]->d[1],
							 elem.node[0]->d[2],
							 elem.node[1]->d[0],
							 elem.node[1]->d[1],
							 elem.node[1]->d[2],
							 elem.cangle);                     /*[DRCCOS]*/

      tmatrix=transmatrix(drccos);               /*TRANSFORMATION MATRIX.*/
      estiff=assememtx(elem);                /*ELASTIC MATRIX OF ELEMENT.*/

      estiff=modifyhinge(elem,estiff,af);               /*MODIFY MATRIX.*/
      estiff=transformation(estiff,tmatrix);             /*[K]=[Tt][k][T]*/
      assemgstiffness(gmtx,estiff,&elem);           /*ASSEMBLAGE ELASTIC.*/
#if 0
      if(nlap==1 && iteration==1)
      {
        modifycmq201(melem,&elem);
        assemcmq201(elem,tmatrix,confs,gvct1); /*ASSEMBLAGE CMQ AS LOADS.*/
        for(ii=0;ii<2;ii++)                          /*ELEM STRESS RESET.*/
        {
          for(jj=0;jj<6;jj++)
          {
            elem.stress[ii][jj]=0.0;
            (melem+(elem.loff))->stress[ii][jj]=elem.stress[ii][jj];
          }
        }
      }
#endif

      for(ii=0;ii<=2;ii++) free(*(drccos+ii));
      free(drccos);
      for(ii=0;ii<=11;ii++) free(*(tmatrix+ii));
      free(tmatrix);
      for(ii=0;ii<=11;ii++) free(*(estiff+ii));
      free(estiff);
    }

	for(ii=1;ii<=msize;ii++)
	{
	  for(jj=1;jj<=ii;jj++)
	  {
		gg=0.0;
		kk=0.0;

		gread(gemtx,ii,jj,&gg);

		if(gg!=0)
		{
		  gread(gmtx,ii,jj,&kk);
		  kk+=gg;
		  gwrite(gmtx,ii,jj,kk);          /*ASSEMBLAGE GLOBAL MATRIX.*/
		}
	  }
	}
//	gcomponentadd(gmtx,gemtx,msize);

    sprintf(string,"GLOBAL MATRIX %ld COMPS ASSEMBLED.",comps);
    laptime(string,t0);

    overlayhdc(*(wdraw.childs+1),SRCPAINT);             /*UPDATE DISPLAY.*/

    /*errormessage("ASSEMBLAGE GLOBAL VECTOR.");*/
    assemconf(confs,gvct,safety,nnode);                  /*GLOBAL VECTOR.*/

	if(fonl!=NULL && /*nlap==1*/iteration==1) fprintf(fonl,"\"TARGET FORCE {F}\"\n");
	for(i=0;i<nnode;i++)
	{
	  if(fonl!=NULL && /*nlap==1*/iteration==1)
	  {
/*
		fprintf(fonl,"NODE:%5ld %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f\n",(nodes+i)->code,
				*(gvct+(6*i+0)),*(gvct+(6*i+1)),*(gvct+(6*i+2)),
				*(gvct+(6*i+3)),*(gvct+(6*i+4)),*(gvct+(6*i+5)));
*/
	  }

	  for(ii=0;ii<nfnodes;ii++)
	  {
		if((nodes+i)->code==fnodes[ii] && fiteration!=NULL)
		{
		  tforce1[ii]=*(gvct+(6*i+2));
		  tforce2[ii]=*(gvct2+(6*i+2));
	    }
	  }
	}

    if(fonl!=NULL && /*nlap!=1*/iteration!=1)
    {
      fprintf(fonl,"\"TRUE FORCE {F-dF}:BEFORE MODIFY CONF\"\n");
      for(i=0;i<nnode;i++)
      {
/*
        fprintf(fonl,"NODE:%5ld %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f\n",(nodes+i)->code,
                *(gvct2+(6*i+0)),*(gvct2+(6*i+1)),*(gvct2+(6*i+2)),
                *(gvct2+(6*i+3)),*(gvct2+(6*i+4)),*(gvct2+(6*i+5)));
*/
      }
    }
    residual=0.0;  /***UJIOKA***/
    for(i=0;i<msize;i++)
    {
      if((confs+i)->iconf==1) *(gvct2+i)=0.0;
      *(gvct+i)-=*(gvct2+i);   /*UNBALANCED FORCE.*/
      residual+=*(gvct+i)**(gvct+i);  /***UJIOKA***/
    }
    residual=pow(residual,0.5);

    modifygivend(gmtx,gvct,confs,nnode);   /*0:LOAD 1:DISPLACEMENT.*/

    if(fonl!=NULL && /*nlap!=1*/iteration!=1)
    {
      fprintf(fonl,"\"TRUE FORCE {F-dF}:AFTER MODIFY CONF\"\n");
      for(i=0;i<nnode;i++)
      {
/*
        fprintf(fonl,"NODE:%5ld %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f\n",(nodes+i)->code,
                *(gvct2+(6*i+0)),*(gvct2+(6*i+1)),*(gvct2+(6*i+2)),
                *(gvct2+(6*i+3)),*(gvct2+(6*i+4)),*(gvct2+(6*i+5)));
*/
      }
      fprintf(fonl,"\"UNBALANCED FORCE {dF}\"\n");
      for(i=0;i<nnode;i++)
      {
/*
        fprintf(fonl,"NODE:%5ld %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f\n",(nodes+i)->code,
                *(gvct+(6*i+0)),*(gvct+(6*i+1)),*(gvct+(6*i+2)),
                *(gvct+(6*i+3)),*(gvct+(6*i+4)),*(gvct+(6*i+5)));
*/
      }
    }

    laptime("CROUT LU DECOMPOSITION.",t0);
	nline=croutludecomposition(gmtx,
                               gvct,confs,
                               6*nnode,
                               &determinant,&sign);        /*[K]{dU}={dF}*/

    sprintf(string,"LOG{DETERMINANT}=%.5E SIGN=%.1f COMPS=%ld",
            determinant,sign,comps);
    errormessage(string);
    if(fonl!=NULL) fprintf(fonl,"%s\n",string);

    sprintf(string,"SAFETY FACTOR=%.5f",safety);
    errormessage(string);
    if(fonl!=NULL) fprintf(fonl,"%s\n",string);

    if(sign<=0.0)
    {
      for(ii=1;ii<=msize;ii++)
      {
        gg=0.0;
        gread(gmtx,ii,ii,&gg);

        if(gg<0.0)
        {
          sprintf(string,"INSTABLE TERMINATION AT NODE %ld.",
                  (nodes+int((ii-1)/6))->code);
          errormessage(" ");
          errormessage(string);
          if(fonl!=NULL) fprintf(fonl,"%s\n",string);
        }
      }

//      sprintf(string,"INSTABLE TERMINATION AT NODE %ld.",
//              (af->nodes+(int)(nline/6))->code);

//      errormessage(" ");
//      errormessage(string);

//      if(fonl!=NULL) fprintf(fonl,"%s\n",string);

      laptime("\0",t0);

      fclose(fin);
      fclose(fonl);
      fclose(fout);
      fclose(fiteration);
      /*fclose(felem);*/
      /*fclose(fdisp);*/
      /*fclose(freact);*/

      /*gfree(kmtx,nnode);*/  /*FREE ELASTIC MATRIX.*/
      gfree(gemtx,nnode); /*FREE GEOMETRIC MATRIX.*/
      gfree(gmtx,nnode);  /*FREE GLOBAL MATRIX.*/
      free(gvct);
      free(gvct1);
      free(gvct2);
      /*free(estress2);*/
      /*free(confs);*/

      memory2=availablephysicalmemory("REMAIN:");
      sprintf(string,"CONSUMPTION:%ld[BYTES]",(memory0-memory2));
      errormessage(string);

      return 1;
    }

    laptime("OUTPUT INTO FILE.",t0);

    if(nlap==laps || nlap==2*laps)
    {
      if(fout!=NULL)
      {
        fprintf(fout,"\n\n");
        fprintf(fout,"** FORCES OF MEMBER\n\n");
        fprintf(fout,"  NO   KT NODE         N        Q1        Q2");
        fprintf(fout,"        MT        M1        M2\n\n");
      }
    }

    if(fonl!=NULL) fprintf(fonl,"\"DISPLACEMENT\"\n");
    outputdisp(gvct,fonl,nnode,nodes);  /*INCREMENTAL DISPLACEMENT.*/
    /*while(!GetAsyncKeyState(VK_LBUTTON))
    ;*/                                   /*LEFT CLICK TO CONTINUE.*/

    if(fonl!=NULL) fprintf(fonl,"\"STRESS\"\n");

    for(i=1;i<=msize;i++)        /*GEOMETRIC MATRIX INITIALIZATION.*/
    {
      ge=(gemtx+(i-1))->down; /*NEXT OF DIAGONAL.*/

      while(ge!=NULL) /*CLEAR ROW.*/
      {
        p=ge;
        ge=ge->down;
        free(p);
      }

      ginit.m=(unsigned short int)i;
      /*ginit.n=(unsigned short int)i;*/

      *(gemtx+(i-1))=ginit;
	  *(gvct2+(i-1))=0.0;           /*GLOBAL VECTOR INITIALIZATION.*/
    }

    updateform(ddisp,gvct,nnode);               /*FORMATION UPDATE.*/

    for(i=1;i<=nelem;i++)                   /*STRESS OUTPUT,UPDATE.*/
    {
      inputelem(elems,melem,i-1,&elem);

      inputnode(ddisp,elem.node[0]);
      inputnode(ddisp,elem.node[1]);

      elem.sect=(elems+i-1)->sect;             /*READ SECTION DATA.*/

///      estress=elemstress(&elem,gvct,melem,fonl,func);
      estress=elemstressnl(&elem,gvct,melem,af);                /*{f}.*/

      estress2=(double *)malloc(12*sizeof(double));  /*ACCUMULATION STRESS*/
      for(ii=0;ii<=1;ii++)                                /*{f+df}.*/
      {
        for(jj=0;jj<6;jj++)
        {
          *(estress2+(6*ii+jj))=0.0;
          *(estress2+(6*ii+jj))=elem.stress[ii][jj];

/*if(jj>=3) *(estress2+(6*ii+jj))=0.0;*/ /*NONLINEAR TEST*/
        }
      }

      /*outputstress(elem,estress,fonl,func);
      if(nlap==laps) outputstressnl(elem,estress,fout);*/
      if(nlap==laps || nlap==2*laps)
      {
        outputstress(elem,estress,fonl,func);
        outputstressnl(elem,estress,fout);
      }

	  drccos=directioncosine(elem.node[0]->d[0],
                             elem.node[0]->d[1],
                             elem.node[0]->d[2],
                             elem.node[1]->d[0],
                             elem.node[1]->d[1],
                             elem.node[1]->d[2],
                             elem.cangle);               /*[DRCCOS]*/

      tmatrix=transmatrix(drccos);         /*TRANSFORMATION MATRIX.*/

      estiff=assemgmtx(elem,estress2);/*GEOMETRIC MATRIX OF ELEMENT.*/

#if 0
for(ii=0;ii<12;ii++)
{
  fprintf(fonl,"[kg](%2d)",ii);
  for(jj=0;jj<=ii;jj++)
  {
    fprintf(fonl," %14.5f",estiff[ii][jj]);
  }
  fprintf(fonl,"\n");
}
#endif

//      estiff=assemgmtx(elem,estress);/*GEOMETRIC MATRIX OF ELEMENT.*/
//      estiff=modifyhinge(elem,estiff);             /*MODIFY MATRIX.*/
#if 0
for(ii=0;ii<12;ii++)
{
  fprintf(fonl,"[kg](%2d)",ii);
  for(jj=0;jj<=ii;jj++)
  {
    fprintf(fonl," %14.5f",estiff[ii][jj]);
  }
  fprintf(fonl,"\n");
}
#endif

      estiff=transformation(estiff,tmatrix);       /*[K]=[Tt][k][T]*/
      assemgstiffness(gemtx,estiff,&elem);  /*ASSEMBLAGE GEOMETRIC.*/

      tt=matrixtranspose(tmatrix,12);                       /*[Tt].*/

//Modify by fukushima///////////////////////////////////////////////////////////
	  estress2=matrixvectorIII(tt,estress2,12); /*TRUE FORCE OF ELEMENT {F}=[Tt]{f}.*/
	  /*estress2=matrixvector(tt,estress2,12);*/ /*TRUE FORCE OF ELEMENT {F}=[Tt]{f}.*/
////////////////////////////////////////////////////////////////////////////////
      /*estress2=matrixvector(tmatrix,estress2,12);*/ /*TRUE FORCE OF ELEMENT {F}=[T]{f}.*/
      /*estress=matrixvector(tmatrix,estress,12);*/   /*{F}=[T]{f}.*/

	  for(ii=0;ii<6;ii++)                             /*TRUE FORCE.*/
      {
        /**(gvct2+((6*elem.node[0]->loff)+ii))-=*(estress2+ii);*/
        /**(gvct2+((6*elem.node[1]->loff)+ii))-=*(estress2+(6+ii));*/
        *(gvct2+((6*elem.node[0]->loff)+ii))+=*(estress2+ii);
        *(gvct2+((6*elem.node[1]->loff)+ii))+=*(estress2+(6+ii));
      }

//150220 fukushima for arclm201/////////////////////////////////////////////////
      for(ii=0;ii<=2;ii++) free(*(drccos+ii));
      free(drccos);
      for(ii=0;ii<=11;ii++) free(*(tmatrix+ii));
      free(tmatrix);
      for(ii=0;ii<=11;ii++) free(*(tt+ii));
      free(tt);
////////////////////////////////////////////////////////////////////////////////

	  free(estress);
      free(estress2);
      for(ii=0;ii<=11;ii++) free(*(estiff+ii));
      free(estiff);
    }

    /*if(wsurf.hwnd!=NULL)
    {
      drawyieldsurface((wsurf.childs+1)->hdcC,
                       (wsurf.childs+1)->vparam,
                       SURFACEX,SURFACEY,SURFACEZ,
                       af->fsurface);
      overlayhdc(*(wsurf.childs+1),SRCPAINT);*/   /*UPDATE DISPLAY.*/
    /*}*/
    /*while(!GetAsyncKeyState(VK_LBUTTON))
    ;*/                                   /*LEFT CLICK TO CONTINUE.*/

    /*updateform(ddisp,gvct,nnode);*/           /*FORMATION UPDATE.*/

    if(fonl!=NULL)
    {
      fprintf(fonl,"\"CURRENT FORM\"\n");
      for(ii=0;ii<nnode;ii++)
      {
        sprintf(string,"NODE:%5ld {U}=",(nodes+ii)->code);
        for(jj=0;jj<6;jj++)
        {
          loffset=6*ii+jj;
          sprintf(s," %16.12f",*(ddisp+loffset));
          strcat(string,s);
        }
//        fprintf(fonl,"%s\n",string);

		for(iii=0;iii<nfnodes;iii++)
		{
		  if((nodes+ii)->code==fnodes[iii] && fiteration!=NULL)
		  {
			if(safety!=0.0)
			{
			  fprintf(fiteration,"LAP: %3d / %3d NODE %3d {Fz}= %18.12f {U}= %16.12f test= %16.12f\n",
					  nlap,laps,(nodes+ii)->code,tforce1[iii],*(ddisp+6*ii+2),
					  residual/tolerance/norm/safety);
			}
			else
			{
			  fprintf(fiteration,"LAP: %3d / %3d NODE %3d {Fz}= %18.12f {U}= %16.12f\n",
					  nlap,laps,(nodes+ii)->code,tforce1[iii],*(ddisp+6*ii+2));
			}
			fprintf(fiteration,"LAP: %3d / %3d NODE %3d {Fz}= %18.12f {U}= %16.12f\n",
					nlap,laps,(nodes+ii)->code,*(gvct2+(6*ii+2)),*(ddisp+6*ii+2));
		  }
		}
      }
      fprintf(fonl,"\"REACTION\"\n");
      outputreaction(gmtx,gvct,nodes,confs,dreact,fonl,nnode);
    }

    if((nlap==laps || nlap==2*laps) && fout!=NULL)
    {
      fprintf(fout,"\n\n");
      fprintf(fout,"** DISPLACEMENT OF NODE\n\n");
      fprintf(fout,"  NO          U          V          W");
      fprintf(fout,"        KSI        ETA      OMEGA\n\n");

      for(ii=0;ii<nnode;ii++)
      {
        sprintf(string,"%4d",(nodes+ii)->code);
        for(jj=0;jj<6;jj++)
        {
          loffset=6*ii+jj;
		  sprintf(s," %10.6f",*(ddisp+loffset)-*(iform+loffset));
		  strcat(string,s);
        }
        fprintf(fout,"%s\n",string);
      }

      fprintf(fout,"\n\n");
      fprintf(fout,"** INITIAL FORM\n\n");
      fprintf(fout,"  NO          U          V          W");
      fprintf(fout,"        KSI        ETA      OMEGA\n\n");

      for(ii=0;ii<nnode;ii++)
      {
        sprintf(string,"%4d",(nodes+ii)->code);
        for(jj=0;jj<6;jj++)
        {
          loffset=6*ii+jj;
          sprintf(s," %10.6f",*(iform+loffset));
          strcat(string,s);
        }
        fprintf(fout,"%s\n",string);
      }

      fprintf(fout,"\n\n");
      fprintf(fout,"** REACTION\n\n");
      fprintf(fout,"  NO  DIRECTION              R    NC\n\n");

      outputreactionnl(gmtx,gvct,nodes,confs,dreact,fout,nnode);/*REACTION.*/
    }

    /***UJIOKA FOR LOAD INCREMENTAL***/
	if( (residual<tolerance*norm*safety || iteration>=maxiteration
	     || (residual<tolerance*norm && nlap>laps) /*&& iteration!=1*/))
    {
      nlap++;
      iteration=0;
	  for(i=0;i<nelem;i++)
	  {
		inputelem(elems,melem,i,&elem);

		elem.sect=(elems+i)->sect;             /*READ SECTION DATA.*/

		if ((elem.sect)->code == 501||
			(elem.sect)->code == 502||
			(elem.sect)->code == 503||
			(elem.sect)->code == 504
		)
		{
		  code = (elem.sect)->code + 10;
		  if ((elem.sect)->code == 501)  Ny=-0.035;
//          if ((elem.sect)->code == 501)  Ny=-0.105;
		  if ((elem.sect)->code == 502)  Ny=-0.033;
		  if ((elem.sect)->code == 503)  Ny=-0.057;
		  if ((elem.sect)->code == 504)  Ny=-0.115;
		  if(elem.stress[0][0]<Ny && nlap<=laps)     /*yield*/
		  {
//     		sprintf(str,"ELEM %d :ESEC %d \n",elem.code,elem.sect->code);
//    	    MessageBox(NULL,str,"yield",MB_OK);
			for(j=0;j<af->nsect;j++)
			{
			  if(code==(af->sects+j)->code)
			  {
				elem.sect=af->sects+j;
				break;
			  }
			}
//     	   	sprintf(str,"ELEM %d :ESEC %d \n",elem.code,(elem.sect)->code);


//          errormessage(str);
//        	MessageBox(NULL,str,"changed(compression)",MB_OK);
		  }
		}
		if ((elem.sect)->code == 511||
			(elem.sect)->code == 512||
			(elem.sect)->code == 513||
			(elem.sect)->code == 514
		)
		{
		  code = (elem.sect)->code - 10;
		  if ((elem.sect)->code == 511)  Ny=-0.035;
//          if ((elem.sect)->code == 511)  Ny=-0.105;
		  if ((elem.sect)->code == 512)  Ny=-0.033;
		  if ((elem.sect)->code == 513)  Ny=-0.057;
		  if ((elem.sect)->code == 514)  Ny=-0.115;
		  if(elem.stress[0][0]>Ny || nlap>laps)     /*elastic*/
		  {
//     		sprintf(str,"ELEM %d :ESEC %d \n",elem.code,elem.sect->code);
//    	    MessageBox(NULL,str,"yield",MB_OK);
			for(j=0;j<af->nsect;j++)
			{
			  if(code==(af->sects+j)->code)
			  {
				elem.sect=af->sects+j;
				break;
			  }
			}
//     	   	sprintf(str,"ELEM %d :ESEC %d \n",elem.code,(elem.sect)->code);


//          errormessage(str);
//        	MessageBox(NULL,str,"changed(compression)",MB_OK);
		  }
		}

		(af->elems+i)->sect=elem.sect;    /*section data update*/

	  }
	  clearwindow(*(wdraw.childs+1));
	  (wdraw.childs+1)->vparam.vflag.ev.deformation=1;
	  drawarclmframe((wdraw.childs+1)->hdcC,
					 (wdraw.childs+1)->vparam,arc,0,ONSCREEN);
	  SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);
	}
	iteration++;
    sprintf(string,"\nITERATION:%d",iteration);
    errormessage(string);
    /***UJIOKA FOR LOAD INCREMENTAL***/


    t1=laptime("\0",t0);

    memory2=availablephysicalmemory(NULL);
    sprintf(string,"CONSUMPTION:%ld[BYTES]",(memory1-memory2));
    errormessage(string);

    errormessage("L:CONTINUE R:ABORT");            /*L=LEFT R=RIGHT*/
    while(!GetAsyncKeyState(VK_LBUTTON))  /*LEFT CLICK TO CONTINUE.*/
    {
      if(GetAsyncKeyState(VK_RBUTTON))      /*RIGHT CLICK TO ABORT.*/
      {
        fclose(fin);
        /*fclose(felem);*/
        /*fclose(fdisp);*/
        /*fclose(freact);*/

        /*gfree(kmtx,nnode);*/  /*FREE ELASTIC MATRIX.*/
        gfree(gemtx,nnode); /*FREE GEOMETRIC MATRIX.*/
        gfree(gmtx,nnode);  /*FREE GLOBAL MATRIX.*/
        free(gvct);
        free(gvct1);
        free(gvct2);
        free(estress2);
        /*free(confs);*/

        errormessage(" ");
        errormessage("ABORTED.");
        if(fonl!=NULL) fprintf(fonl,"ABORTED.\n");

        fclose(fonl);

        laptime("\0",t0);
        return 1;
      }
      t2=clock();
      time=(t2-t1)/CLK_TCK;
      if(time>=WAIT) break;               /*CONTINUE AFTER WAITING.*/
    }

    /*fclose(fonl);*/
  }                                        /*REPEAT UNTIL INSTABLE.*/

  if((wdraw.childs+1)->hdcC!=NULL &&
     melem!=NULL && ddisp!=NULL)                 /*DRAW LAST FRAME.*/
  {
    for(i=1;i<=nelem;i++)
    {
      inputelem(elems,melem,i-1,&elem);
      for(ii=0;ii<=1;ii++) /*COPY HINGE DATA.*/
      {
        for(jj=0;jj<=5;jj++)
        {
          (elems+i-1)->iconf[ii][jj]=elem.iconf[ii][jj];
        }
      }

      inputnode(ddisp,elem.node[0]);
      inputnode(ddisp,elem.node[1]);

	  drawglobalwire((wdraw.childs+1)->hdcC,
                     (wdraw.childs+1)->vparam,
                     *af,elem,255,255,255,
							  255,255,255,0,ONSCREEN/*,i*/);
    }
    overlayhdc(*(wdraw.childs+1),SRCPAINT);       /*UPDATE DISPLAY.*/
  }

  fclose(fin);
  fclose(fout);
  fclose(fiteration);
  /*fclose(felem);*/
  /*fclose(fdisp);*/
  /*fclose(freact);*/

  /*gfree(kmtx,nnode);*/  /*FREE ELASTIC MATRIX.*/
  gfree(gemtx,nnode); /*FREE GEOMETRIC MATRIX.*/
  gfree(gmtx,nnode);  /*FREE GLOBAL MATRIX.*/
  /*free(gvct);*/
  /*free(confs);*/
  /*free(estress2);*/

  af->eigenvec=(double **)malloc(1*sizeof(double *));
  *((af->eigenvec)+0)=gvct;

  errormessage(" ");
  errormessage("COMPLETED.");
  if(fonl!=NULL) fprintf(fonl,"COMPLETED.\n");

  if(fonl!=NULL) fclose(fonl);

  memory2=availablephysicalmemory("REMAIN:");
  sprintf(string,"CONSUMPTION:%ld[BYTES]",(memory0-memory2));
  errormessage(string);
  errormessage(" ");

  return 0;
}/*arclm202*/

