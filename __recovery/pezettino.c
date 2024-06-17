/*CR PROP*/
struct corotational{double mode[8];   //DEFORMATION MODE VECTOR AND PhiAsymX LENGTH0
                    double force[6];  //INTERNAL ELEMENT FORCE
                    double basis[9];  //ELEMENT BASIS OF COROTATIONAL BEAM
                   };

/*PEZETTINO STATIC NONLINEAR ANALYSIS*/  /*UNDER CONSTRUCTION*/
int pezettino_staticnonlinear(struct arclmframe *af,int idinput);

/***FOR PEZETTINO 003 FILE TYPE (-inl2 style)***/
int saveasarclm2(char *fname,struct arclmframe *af);
void inputtexttomemory2(FILE *ftext,struct arclmframe *af);
void inputinit2(FILE *fin,int *nnode,int *nelem,int *nsect);
void readsect2(FILE *fin,struct osect *sect);

/***FOR PEZETTINO 003 INCREMENTAL ANALYSIS***/
int updateorganization(struct organ *org);
/*CR PROP*/
void corotationaloutputtomemory(FILE *ftext,struct arclmframe *af,struct corotational *corot);
void outputcorotational(FILE *fout,struct arclmframe *af,struct corotational *corot);
/*fstart*/
int saveloadasfstart(char *fname,struct arclmframe *af);

int pezettino_staticnonlinear(struct arclmframe *af,int idinput)
/*
  PEZETTINO STATIC NONLINEAR ANALYSIS
  UNDER CONSTRUCTION
*/
{
  DWORD memory0,memory1,memory2;

  FILE *fin,*fout,*fonl;                           /*FILE 8 BYTES*/
  FILE *fiteration,*fcr,*fresidual;
  FILE *fstart,*fnew;
  double *ddisp,*dreact,*iform;
  struct memoryelem *melem;
  /*FILE *felem,*fdisp,*freact;*/
  char dir[]=DIRECTORY;                            /*DATA DIRECTORY*/
  char s[80],string[400]/*,fname[256]*/;
  int i,ii,jj;
  int nnode,nelem,nsect,nreact;
  long int loffset,msize,nline;
  long int time;
  struct gcomponent ginit={0,0,0.0,NULL};
  struct gcomponent /**kmtx,*/*gemtx,*gmtx,*k,*ge,*g,*p;/*GLOBAL MATRIX*/
  double gg,kk;
  double *gvct,*gvct0,*gvct1,*gvct2;                       /*GLOBAL VECTOR*/
  double **drccos,**tmatrix,**estiff,*estress,*estress2,**tt;
  double determinant,sign,safety,dsafety;
  double func[2],tforce1,tforce2;

  /**** setting NONLINEAR ****/
  int laps=500;                 /*nIStep*/
  int maxiteration=20;         /*nJStep*/
  int firstnlaps=200;           /*firstNStep*/
  double firstsafety=0.0001;    /*firstIncrementRatio*/
  double tolerance=0.00001;   /*epsilon*/
  int fnode=243;               /*node*/
  int fdirection=2;            /*direction*/
  /**** setting NONLINEAR ****/

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
  fresidual=fopen("residual.txt","w");
  fstart=fgetstofopenII(dir,"r","fstart.txt");
  fnew=fgetstofopenII(dir,"w","newcoord.txt");

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
  gvct0=(double *)malloc(msize*sizeof(double)); /*fstart*/          /*GLOBAL VECTOR.*/
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

  corot=(struct corotational *)
         malloc((arc.nelem)*sizeof(struct corotational));

  fcr=fgetstofopenII(dir,"r","corotational.txt");  /*OPEN FILE.*/
  if(fcr!=NULL &&
     MessageBox(NULL,"LOAD CR PROP.","PEZETTINO",MB_OKCANCEL)==IDOK);
  {
     corotationaloutputtomemory(fcr,&arc,&*corot);
     fclose(fcr);
  }

  GetAsyncKeyState(VK_LBUTTON);                   /*CLEAR KEY LEFT.*/
  GetAsyncKeyState(VK_RBUTTON);                  /*CLEAR KEY RIGHT.*/

  drawglobalaxis((wdraw.childs+1)->hdcC,(wdraw.childs+1)->vparam,
                 0,0,255);                      /*DRAW GLOBAL AXIS.*/

  nlap=1;
  iteration=1;
  residual=0.0;
//  for(nlap=1;nlap<=laps;nlap++)
  while(nlap<=laps)      /***UJIOKA FOR LOAD INCREMENTAL***/
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

      /*estiff=modifyhinge(elem,estiff);*/               /*MODIFY MATRIX.*/
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

      if((nodes+i)->code==fnode && fiteration!=NULL)
      {
        tforce1=*(gvct+(6*i+2));
        tforce2=*(gvct2+(6*i+2));
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
    residual=0;  /***UJIOKA***/
    for(i=0;i<msize;i++)
    {
      if((confs+i)->iconf==1) *(gvct2+i)=0.0;
      *(gvct+i)-=*(gvct2+i);   /*UNBALANCED FORCE.*/
      residual+=*(gvct+i)**(gvct+i);  /***UJIOKA***/
    }

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
	  estress=elemstressnl(&elem,gvct,melem,fonl);                /*{f}.*/

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
      if(nlap==laps)
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


    /***UJIOKA FOR LOAD INCREMENTAL***/
    if( (residual<tolerance || iteration>=maxiteration)&& iteration!=1)
    {
      nlap++;
      iteration=0;
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
//        fprintf(fonl,"%s\n",string);

        if((nodes+ii)->code==fnode && fiteration!=NULL)
        {
          fprintf(fiteration,"LAP: %3d / %3d NODE %3d {Fz}= %18.12f {U}= %16.12f\n",
                  nlap,laps,(nodes+ii)->code,tforce1,*(ddisp+6*ii+2));
          fprintf(fiteration,"LAP: %3d / %3d NODE %3d {Fz}= %18.12f {U}= %16.12f\n",
                  nlap,laps,(nodes+ii)->code,*(gvct2+(6*ii+2)),*(ddisp+6*ii+2));
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

  fcr=fgetstofopenII(dir,"w","corotational_new.txt");  /*OPEN FILE.*/
  if(fcr==NULL) return -1;
  outputcorotational(fcr,&arc,&*corot);
  fclose(fcr);

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
}/*arclm201*/

/***FOR PEZETTINO 003 FILE TYPE (-inl2 style)***/
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
            0.0/*(af->sects+i)->alpha*/, /*線膨張率*/
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
  /*LONGREACTION:UNDER CONSTRUCTION.*/

  return;
}/*inputtexttomemory*/

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
}/*inputinit*/

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
}/*readsect*/

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

//  arc.nlaps=1;

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

  corotationaloutputtomemory(ftxt,&arc,&*corot);
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
  ii=0;      /*index of arclm frame*/
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
//    	    MessageBox(NULL,str,"compression",MB_OK);
            for(j=0;j<((wdraw.childs+1)->org.nsect);j++)
            {
              if(code==((wdraw.childs+1)->org.sects+j)->code)
              {
                /*UPDATE CR PROP*/
                (corot+ii)->mode[3]/=(
                                     ((((wdraw.childs+1)->org.sects+j)->figs+0)->prop->E)
                                     /(((org->elems+i)->sect->figs+0)->prop->E)
                                     );
                (corot+ii)->mode[7]=dl-
                                    ((dl-(corot+ii)->mode[7])*
                                     (((org->elems+i)->sect->figs+0)->prop->E)
                                     /((((wdraw.childs+1)->org.sects+j)->figs+0)->prop->E)
                                    );
                (org->elems+i)->sect=(wdraw.childs+1)->org.sects+j;
                break;
              }
            }
     	   	sprintf(str,"ELEM %d :ESEC %d \n",(org->elems+i)->code,((org->elems+i)->sect)->code);
            errormessage(str);
//        	MessageBox(NULL,str,"changed(compression)",MB_OK);
     	}
    }
  #if 0
  	else if (((org->elems+i)->sect)->code == 511||
        ((org->elems+i)->sect)->code == 512||
        ((org->elems+i)->sect)->code == 513||
        ((org->elems+i)->sect)->code == 614
       )
    {
        code = ((org->elems+i)->sect)->code - 10;
        if((arc.elems+ii)->srate[0]<1.0)
     	{
           sprintf(str,"ELEM %d :ESECT %d \n",(org->elems+i)->code,((org->elems+i)->sect)->code);
//         MessageBox(NULL,str,"tension",MB_OK);
            for(j=0;j<((wdraw.childs+1)->org.nsect);j++)
            {
              if(code==((wdraw.childs+1)->org.sects+j)->code)
              {
                (org->elems+i)->sect=(wdraw.childs+1)->org.sects+j;
                break;
              }
            }
           sprintf(str,"ELEM %d :ESEC %d \n",(org->elems+i)->code,((org->elems+i)->sect)->code);
           errormessage(str);
//         MessageBox(NULL,str,"changed(tension)",MB_OK);
     	}
    }
  #endif
  }

#define FILEFORCHANGESECT  "changesect_test.inp"  /*file for test*/
  fin=fopen(FILEFORCHANGESECT,"w");
  saveorganizationforchangsectautomatic(fin,&((wdraw.childs+1)->org),&arc,
                                  &((wdraw.childs+1)->vparam));
  fclose(fin);


  ftxt=fgetstofopenII(dir,"w","corotational_new.txt");  /*OPEN FILE.*/
  if(ftxt==NULL) return -1;
  outputcorotational(ftxt,&arc,&*corot);
  fclose(ftxt);

  sprintf(str,"Saved as %s \n",FILEFORCHANGESECT);
  MessageBox(NULL,str,"change section automatic",MB_OK);

  globalmessageflag=1;
  globaldrawflag=1;

return 0;
}/*updateorganization*/

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
}/*saveasarclm2*/

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
    for(j=0;j<8;j++)
    {
      (corot+i-1)->mode[j]=strtod(*(data+k),NULL);
      k++;
    }

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
                        (corot+i-1)->mode[6],(corot+i-1)->mode[7]);
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
}/*corotationaloutputtomemory*/

int arclm202(struct arclmframe *af,int idinput)
/*
  Non-Linear Shape Analysis for KIRIGAMI
*/
{
  DWORD memory0,memory1,memory2;

  FILE *fin,*fout,*fonl;                           /*FILE 8 BYTES*/
  FILE *fiteration;
  FILE *fstart,*fnew;
  double *ddisp,*dreact,*iform;
  struct memoryelem *melem;
  /*FILE *felem,*fdisp,*freact;*/
  char dir[]=DIRECTORY;                            /*DATA DIRECTORY*/
  char s[80],string[400],str[256]/*,fname[256]*/;
  int i,ii,j,jj;
  int nnode,nelem,nsect,nreact;
  int code;
  long int loffset,msize,nline;
  long int time;
  struct gcomponent ginit={0,0,0.0,NULL};
  struct gcomponent /**kmtx,*/*gemtx,*gmtx,*k,*ge,*g,*p;/*GLOBAL MATRIX*/
  double gg,kk;
  double *gvct,*gvct0,*gvct1,*gvct2;                       /*GLOBAL VECTOR*/
  double **drccos,**tmatrix,**estiff,*estress,*estress2,**tt;
  double determinant,sign,safety,dsafety;
  double func[2],tforce1,tforce2;
  double Ny;

  /**** setting NONLINEAR ****/
  int laps=200;                /*nIStep*/
  int maxiteration=10;         /*nJStep*/
  int firstnlaps=10;           /*firstNStep*/
  double firstsafety=0.00005;   /*firstIncrementRatio*/
  double tolerance=0.0001;     /*epsilon*/
  int fnode=243;               /*node*/
  int fdirection=2;            /*direction*/
  /**** setting NONLINEAR ****/

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
  gvct0=(double *)malloc(msize*sizeof(double)); /*fstart*/          /*GLOBAL VECTOR.*/
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
//  for(nlap=1;nlap<=laps;nlap++)
  while(nlap<=laps)      /***UJIOKA FOR LOAD INCREMENTAL***/
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

      estiff=modifyhinge(elem,estiff);               /*MODIFY MATRIX.*/
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

      if((nodes+i)->code==fnode && fiteration!=NULL)
      {
        tforce1=*(gvct+(6*i+2));
        tforce2=*(gvct2+(6*i+2));
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
    residual=0;  /***UJIOKA***/
    for(i=0;i<msize;i++)
    {
      if((confs+i)->iconf==1) *(gvct2+i)=0.0;
      *(gvct+i)-=*(gvct2+i);   /*UNBALANCED FORCE.*/
      residual+=*(gvct+i)**(gvct+i);  /***UJIOKA***/
    }

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
	  estress=elemstressnl(&elem,gvct,melem,fonl);                /*{f}.*/

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
      if(nlap==laps)
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


    /***UJIOKA FOR LOAD INCREMENTAL***/
    if( (residual<tolerance || iteration>=maxiteration)&& iteration!=1)
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
          if ((elem.sect)->code == 502)  Ny=-0.033;
          if ((elem.sect)->code == 503)  Ny=-0.057;
          if ((elem.sect)->code == 504)  Ny=-0.115;
  	      if(elem.stress[0][0]<Ny)     /*yield*/
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

        if((nodes+ii)->code==fnode && fiteration!=NULL)
        {
          fprintf(fiteration,"LAP: %3d / %3d NODE %3d {Fz}= %18.12f {U}= %16.12f\n",
                  nlap,laps,(nodes+ii)->code,tforce1,*(ddisp+6*ii+2));
          fprintf(fiteration,"LAP: %3d / %3d NODE %3d {Fz}= %18.12f {U}= %16.12f\n",
                  nlap,laps,(nodes+ii)->code,*(gvct2+(6*ii+2)),*(ddisp+6*ii+2));
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

