/*ARCLM001.C FOR WIN32 SINCE 1995.11.24.JUNSATO.*/
/*LAST CHANGE:09-Aug-2018.*/

/*ANALYSIS STATIC LINEAR.*/
/*OUTPUT WITH FORMAT "FRAME3".*/

/*"INPUTFILE SPECIFICATION"*/
/*APPELATION*/
/*NNODE NELEM NSECT*/
/*SECTCODE E POI A Ixx Iyy J (QxMAX.....MzMIN)*/
/*NODECODE X Y Z*/
/*ELEMCODE SECT NODEI NODEJ COORDANGLE BOUNDARY... CMQ...*/
/*NODECODE CONFINEMENTTYPE... CONFINEMENTVALUE...*/

/*"OUTPUTFILE SPECIFICATION"*/
/*APPELATION*/
/*ELEMCODE SECT HEAD N Q1 Q2 MT M1 M2*/
/*              TAIL N Q1 Q2 MT M1 M2*/
/*NODECODE U V W TX TY TZ*/
/*NODECODE DIRECTION REACTIONVALUE CONFINEMENTTYPE*/

/*#include "canhead.h"*/                /*DEFINITION OF COMMAND ID.*/

int arclm001(struct arclmframe *af,int idinput,int idoutput);
int arclm001_lxy(struct arclmframe *af,int idinputs[],int idoutputs[]); // 150515 fukushima for all
int arclmtest(struct arclmframe *af,int idinput,int idoutput);
void initialelem001(struct owire *elems,
                    struct memoryelem *melem,int nelem);
double *elemstress001(struct owire *elem,
                      double *gvct,struct memoryelem *melem,struct arclmframe *af);
double *elemstress002(struct owire *elem,
                      double *gvct,struct memoryelem *melem,
                      long int *moff,struct arclmframe *af);
double *elemstress003(FILE *fout,struct owire *elem,
                      double *gvct,struct memoryelem *melem,struct arclmframe *af);
void updatestress001(struct memoryelem *melem,double *dstress,
                     struct owire *elem);
void outputdisp001(double *gvct,FILE *fout,int nnode,
                   struct onode *nodes);
void outputdisp002(double *gvct,FILE *fout,int nnode,
                   struct onode *nodes,long int *moff);
void outputstress001(struct owire elem,
                     double *estress,FILE *fout);
void outputreaction001(struct gcomponent *gmtx,
                       double *gvct,
                       struct onode *nodes,
                       struct oconf *confs,
                       double *dreact,FILE *fout,int nnode);
void outputreaction002(struct gcomponent *gmtx,
                       double *gvct,
                       struct onode *nodes,
                       struct oconf *confs,
                       double *dreact,FILE *fout,int nnode,
                       long int *moff);
void modifycmq(struct memoryelem *melem,struct owire *elem);
void assemcmq(struct owire elem,double **tmatrix,
              struct oconf *confs,double *gvct);
void assemconf001(struct oconf *confs,double *gvct,
                  double dsafety,int nnode);

void frameoutputtomemory(FILE *ftext,struct arclmframe *af);
void openarclmlastfile(FILE *ftext,struct arclmframe *af);

int arclm001(struct arclmframe *af,int idinput,int idoutput)
{
  DWORD memory0,memory1,memory2;

  FILE *fin,*fout;                                   /*FILE 8 BYTES*/
  FILE *fmtx;
  /* FILE *fmtx; */
  double *ddisp,*dreact;
  struct memoryelem *melem;
  /*FILE *felem,*fdisp,*freact;*/
  char dir[]=DIRECTORY;                            /*DATA DIRECTORY*/
  char s[20],string[1024];
  int i,ii;
  int nnode,nelem,nsect,nreact;
  long int msize,nline,*moff,*noff;
  struct gcomponent ginit={0,0,0.0,NULL};
  struct gcomponent *gmtx,*gmtx2;                   /*GLOBAL MATRIX*/
  double *gvct,*gvct2;                              /*GLOBAL VECTOR*/
  double **drccos,**tmatrix,**estiff,*estress;
  double determinant,sign;
  /*** araki for instable termination ***/
  int m,jj;
  double data;
  /**************************************/
  clock_t t0;

  struct osect *sects;
  struct onode *nodes,*ninit;
  struct owire elem;
  struct owire *elems;
  struct oconf *confs,*confs2;

  memory0=availablephysicalmemory("INITIAL:");   /*MEMORY AVAILABLE*/

  fin=fgetstofopen(dir,"r",idinput);              /*OPEN INPUT FILE*/
  if(fin==NULL)
  {
    errormessage("ACCESS IMPOSSIBLE.");
    return 0;
  }
  fout=fgetstofopen("\0","w",idoutput);               /*OUTPUT FILE*/
  /* fmtx=fopen("matrix001.txt","w"); */

  /*fgets(string,256,fin);*/                    /*INPUT APPELATION.*/
  /*errormessage(string);*/

  t0=clock();                                        /*CLOCK BEGIN.*/

  inputinit(fin,&nnode,&nelem,&nsect);             /*INPUT INITIAL.*/
  sprintf(string,"NODES=%d ELEMS=%d SECTS=%d",nnode,nelem,nsect);
  errormessage(string);
  /*if(fout!=NULL) fprintf(fout,"%s\n",string);*/

  msize=6*nnode; /*SIZE OF GLOBAL MATRIX.*/

  gmtx=(struct gcomponent *)          /*DIAGONALS OF GLOBAL MATRIX.*/
        malloc(msize*sizeof(struct gcomponent));
  gvct=(double *)malloc(msize*sizeof(double));      /*GLOBAL VECTOR*/
  if(gmtx==NULL || gvct==NULL) return 0;
  for(i=0;i<msize;i++)
  {
    (gmtx+i)->down=NULL;            /*GLOBAL MATRIX INITIALIZATION.*/
    *(gvct+i)=0.0;                  /*GLOBAL VECTOR INITIALIZATION.*/
  }

  free(af->sects);
  free(af->nodes);
  free(af->ninit);
  free(af->elems);
  free(af->confs);
  free(af->ddisp);                     /*DISPLACEMENT:6 DIRECTIONS.*/
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
  melem=(struct memoryelem *)malloc(nelem*sizeof(struct memoryelem));
  if(melem==NULL) return 0;

  af->sects=sects;
  af->nodes=nodes;
  af->ninit=ninit;
  af->elems=elems;
  af->confs=confs;
  af->ddisp=ddisp;
  af->melem=melem;

  inputtexttomemory(fin,af);        /*READY TO READ LONG REACTIONS.*/
  nnode=af->nnode;
  nelem=af->nelem;
  nsect=af->nsect;
  nreact=af->nreact;

  initialform(nodes,ddisp,nnode);           /*ASSEMBLAGE FORMATION.*/

  initialelem001(elems,melem,nelem);         /*ASSEMBLAGE ELEMENTS.*/

  dreact=(double *)malloc(nreact*sizeof(double));       /*REACTION.*/
  af->dreact=dreact;
  initialreact(fin,dreact,nreact);     /*ASSEMBLAGE LONG REACTIONS.*/

  GetAsyncKeyState(VK_LBUTTON);                   /*CLEAR KEY LEFT.*/
  GetAsyncKeyState(VK_RBUTTON);                  /*CLEAR KEY RIGHT.*/

  drawglobalaxis((wdraw.childs+1)->hdcC,(wdraw.childs+1)->vparam,
                 0,0,255);                      /*DRAW GLOBAL AXIS.*/

  errormessage("ARCLM001:ANALYSIS LINEAR.");
  memory1=availablephysicalmemory("REMAIN:");    /*MEMORY AVAILABLE*/

  for(i=1;i<=msize;i++)             /*GLOBAL MATRIX INITIALIZATION.*/
  {
    ginit.m=(unsigned short int)i;
    /*ginit.n=(unsigned short int)i;*/
    *(gmtx+(i-1))=ginit;
  }
  comps=msize; /*INITIAL COMPONENTS=DIAGONALS.*/

  fmtx = fopen("matrix001.txt","w");
  laptime("ASSEMBLING GLOBAL MATRIX.",t0);

  for(i=1;i<=nelem;i++)               /*ASSEMBLAGE GLOBAL MATRIX.*/
  {
    elem=*(elems+i-1);                       /*READ ELEMENT DATA.*/

    inputnode(ddisp,elem.node[0]);                         /*HEAD*/
    inputnode(ddisp,elem.node[1]);                         /*TAIL*/

    elem.sect=(elems+i-1)->sect;             /*READ SECTION DATA.*/

    if((wdraw.childs+1)->hdcC!=NULL)     /*DRAW DEFORMED ELEMENT.*/
    {
      drawglobalwire((wdraw.childs+1)->hdcC,
                     (wdraw.childs+1)->vparam,
                     *af,elem,255,255,255,
                              255,255,255,0,ONSCREEN);
    }

    drccos=directioncosine(elem.node[0]->d[0],
                           elem.node[0]->d[1],
                           elem.node[0]->d[2],
                           elem.node[1]->d[0],
                           elem.node[1]->d[1],
                           elem.node[1]->d[2],
                           elem.cangle);               /*[DRCCOS]*/

    tmatrix=transmatrix(drccos);         /*TRANSFORMATION MATRIX.*/
    estiff=assememtx(elem);          /*ELASTIC MATRIX OF ELEMENT.*/
    estiff=modifyhinge(elem,estiff,af);          /*MODIFY MATRIX.*/
    estiff=transformation(estiff,tmatrix);       /*[K]=[Tt][k][T]*/

    assemgstiffness(gmtx,estiff,&elem);             /*ASSEMBLAGE.*/

    modifycmq(melem,&elem);
    assemcmq(elem,tmatrix,confs,gvct);    /*ASSEMBLAGE CMQ AS LOADS.*/

    for(ii=0;ii<=2;ii++) free(*(drccos+ii));
    free(drccos);
    for(ii=0;ii<=11;ii++) free(*(tmatrix+ii));
    free(tmatrix);
    for(ii=0;ii<=11;ii++) free(*(estiff+ii));
    free(estiff);
  }
  for(m=0;m<msize;m++)
  {
    fprintf(fmtx,"%d %25.18f\n",m,gvct[m]);
  }
  sprintf(string,"GLOBAL MATRIX %ld COMPS ASSEMBLED.",comps);
  laptime(string,t0);

  /* if(fmtx!=NULL) */
  /* { */
  /*   outputmtx_conf(fmtx,gmtx,confs,msize); */
  /*   fclose(fmtx); */
  /* } */

  overlayhdc(*(wdraw.childs+1),SRCPAINT);          /*UPDATE DISPLAY.*/

  /*errormessage("ASSEMBLAGE GLOBAL VECTOR.");*/
  assemconf001(confs,gvct,1.0,nnode);               /*GLOBAL VECTOR.*/
  modifygivend(gmtx,gvct,confs,nnode);      /*0:LOAD 1:DISPLACEMENT.*/

#if 0
  if(globalmessageflag==1 &&
     MessageBox(NULL,"Decrease Band Width ?","ARCLM001",MB_OKCANCEL)
     ==IDOK)
  {
    laptime("DECREASE BAND WIDTH.",t0);
    moff=(long int *)malloc(nnode*sizeof(long int));
    noff=(long int *)malloc(nnode*sizeof(long int));
    if(moff==NULL || noff==NULL) return 0;
    for(i=0;i<nnode;i++)
    {
      *(moff+i)=i;
      *(noff+i)=i;
    }
    decreaseband(moff,noff,af);              /*DECREASE BAND WIDTH.*/

    gmtx2=(struct gcomponent *)malloc(msize*sizeof(struct gcomponent));
    gvct2=(double *)malloc(msize*sizeof(double));
    confs2=(struct oconf *)malloc(msize*sizeof(struct oconf));
    if(gmtx2==NULL || gvct2==NULL || confs2==NULL) return 0;
    for(i=0;i<msize;i++)
    {
      ginit.m=(unsigned short int)(i+1);
      *(gmtx2+i)=ginit;

      *(gvct2+i)=0.0;

      (confs2+i)->iconf=0;
      (confs2+i)->value=0.0;
    }

    exchangelines(gmtx,gvct,confs,
                  gmtx2,gvct2,confs2,
                  moff,noff,nnode,0);
  }
  else
  {
    moff=NULL;
    noff=NULL;
    gmtx2=gmtx;
    gvct2=gvct;
    confs2=confs;
  }
#endif

#if 1
  moff=NULL;
  noff=NULL;
  gmtx2=gmtx;
  gvct2=gvct;
  confs2=confs;
#endif

  laptime("CROUT LU DECOMPOSITION.",t0);
  fclose(fmtx);
  nline=croutludecomposition(gmtx2,
                             gvct2,confs2,
                             6*nnode,
                             &determinant,&sign);  /*[K]{dU}={dF}*/

  sprintf(string,"LOG{DETERMINANT}=%.5E SIGN=%.1f COMPS=%ld",
          determinant,sign,comps);
  errormessage(string);
  /*if(fout!=NULL) fprintf(fout,"%s\n",string);*/

  /*** 10.11.26 araki for instable termination ****************************/
    for(m=1;m<=msize;m++)
    {
      gread(gmtx2,m,m,&data);
      /*if(fout!=NULL)
            fprintf(fout,"NODE %ld. DIAGNAL %20.5f\n",
                    (af->nodes+int((m-1)/6))->code,data);*/
      if(data<=0.0)
      {
        sprintf(string,"INSTABLE TERMINATION AT NODE %ld.\n",
                (af->nodes+int((m-1)/6))->code);
        errormessage(string);
      }
    }
  /************************************************************************/

  if(sign<=0.0)
  {
    if(noff==NULL)
    {
      sprintf(string,"INSTABLE TERMINATION AT NODE %ld.",
              (af->nodes+(int)(nline/6))->code);
    }
    else
    {
      sprintf(string,"INSTABLE TERMINATION AT NODE %ld.",
              (af->nodes+(*(noff+(int)(nline/6))))->code);
    }

    errormessage(" ");
    errormessage(string);
    if(fout!=NULL) fprintf(fout,"%s\n",string);

    laptime("\0",t0);

    fclose(fin);

    if(gmtx2==gmtx)
    {
      gfree(gmtx,nnode); /*FREE GLOBAL MATRIX.*/
      free(gvct);
      /*free(confs);*/
    }
    else
    {
      gfree(gmtx,nnode); /*FREE GLOBAL MATRIX.*/
      gfree(gmtx2,nnode);
      free(gvct);
      free(gvct2);
      /*free(confs);*/
      /*free(confs2);*/
    }

    memory2=availablephysicalmemory("REMAIN:");
    sprintf(string,"CONSUMPTION:%ld[BYTES]",(memory0-memory2));
    errormessage(string);

    return 1;
  }

  /*laptime("REEXCHANGE MATRIX.",t0);*/
  /*exchangelines(gmtx,gvct,confs,moff,noff,nnode,1);*/

  laptime("OUTPUT INTO FILE.",t0);

  errormessage("STRESS.");
  if(fout!=NULL)
  {
    fprintf(fout,"\n\n");
    fprintf(fout,"** FORCES OF MEMBER\n\n");
    fprintf(fout,"  NO   KT NODE         N        Q1        Q2");
    fprintf(fout,"        MT        M1        M2\n\n");
  }
  for(i=1;i<=nelem;i++)                   /*STRESS OUTPUT,UPDATE.*/
  {
    inputelem(elems,melem,i-1,&elem);

    inputnode(ddisp,elem.node[0]);
    inputnode(ddisp,elem.node[1]);

    elem.sect=(elems+i-1)->sect;             /*READ SECTION DATA.*/

    /*estress=elemstress001(&elem,gvct,melem);*/
    estress=elemstress002(&elem,gvct2,melem,moff,af);

    outputstress001(elem,estress,fout);
    free(estress);
  }
  if(fout!=NULL) fprintf(fout,"\n\n");

  errormessage("DISPLACEMENT.");
  /*outputdisp001(gvct,fout,nnode,nodes);*/           /*DISPLACEMENT.*/
  outputdisp002(gvct2,fout,nnode,nodes,moff);     /*DISPLACEMENT.*/
  if(fout!=NULL) fprintf(fout,"\n\n");
  /*while(!GetAsyncKeyState(VK_LBUTTON))
  ;*/                                   /*LEFT CLICK TO CONTINUE.*/

  errormessage("REACTION.");
  if(fout!=NULL)
  {
    fprintf(fout,"** REACTION\n\n");
    fprintf(fout,"  NO  DIRECTION              R    NC\n\n");
  }
  /*outputreaction001(gmtx,gvct,nodes,confs,dreact,fout,nnode);*/
  outputreaction002(gmtx2,gvct2,nodes,confs2,dreact,fout,nnode,moff);

  fclose(fin);
  fclose(fout);
  errormessage("FILES CLOSED.");

  /*updateform(ddisp,gvct,nnode);*/               /*FORMATION UPDATE.*/
  updateform2(ddisp,gvct2,nnode,moff);        /*FORMATION UPDATE.*/

  laptime("\0",t0);

  memory2=availablephysicalmemory(NULL);
  sprintf(string,"CONSUMPTION:%ld[BYTES]",(memory1-memory2));
  errormessage(string);

  if((wdraw.childs+1)->hdcC!=NULL &&
     melem!=NULL && ddisp!=NULL)               /*DRAW LAST FRAME.*/
  {
    for(i=1;i<=nelem;i++)
    {
      elem=*(elems+i-1);                     /*READ ELEMENT DATA.*/

      inputnode(ddisp,elem.node[0]);
      inputnode(ddisp,elem.node[1]);

      drawglobalwire((wdraw.childs+1)->hdcC,
                     (wdraw.childs+1)->vparam,
                     *af,elem,255,255,255,
                              255,255,255,0,ONSCREEN);
    }
    overlayhdc(*(wdraw.childs+1),SRCPAINT);     /*UPDATE DISPLAY.*/
  }

  if(gmtx2==gmtx)
  {
    gfree(gmtx,nnode); /*FREE GLOBAL MATRIX.*/
    /*free(gvct);*/
    /*free(confs);*/
  }
  else
  {
    gfree(gmtx,nnode); /*FREE GLOBAL MATRIX.*/
    gfree(gmtx2,nnode);
    /*free(gvct);*/
    /*free(gvct2);*/
    /*free(confs);*/
    /*free(confs2);*/
  }

  af->nlaps=1;
  /*af->eigenvec=(double **)malloc(1*sizeof(double *));*/
  /**((af->eigenvec)+0)=gvct2;*/

  errormessage(" ");
  errormessage("COMPLETED.");

  memory2=availablephysicalmemory("REMAIN:");
  sprintf(string,"CONSUMPTION:%ld[BYTES]",(memory0-memory2));
  errormessage(string);

  return 0;
}/*arclm001*/

/* 150515 fukushima for all */
int arclm001_lxy(struct arclmframe *afs[],int idinputs[],int idoutputs[])
{
  DWORD memory0,memory1,memory2;

  FILE *fin,*fout;                                   /*FILE 8 BYTES*/
  double *ddisp,*dreact;
  struct memoryelem *melem;
  /*FILE *felem,*fdisp,*freact;*/
  char dir[]=DIRECTORY;                            /*DATA DIRECTORY*/
  char s[20],string[1024];
  int i,ii;
  int nnode,nelem,nsect,nreact;
  long int msize,nline,*moff,*noff;
  struct gcomponent ginit={0,0,0.0,NULL};
  struct gcomponent *gmtx,*gmtx2;                   /*GLOBAL MATRIX*/
  double *gvct,*gvct2;                              /*GLOBAL VECTOR*/
  double **drccos,**tmatrix,**estiff,*estress;
  double determinant,sign;
  /*** araki for instable termination ***/
  int m,jj;
  int lxy;
  struct gcomponent *lu;
  double data;
  /**************************************/
  clock_t t0;

  struct osect *sects;
  struct onode *nodes,*ninit;
  struct owire elem;
  struct owire *elems;
  struct oconf *confs,*confs2;

  memory0=availablephysicalmemory("INITIAL:");   /*MEMORY AVAILABLE*/

  for(lxy=0;lxy<3;lxy++)
  {
    fin=fgetstofopen(dir,"r",idinputs[lxy]);              /*OPEN INPUT FILE*/
    if(fin==NULL)
    {
      errormessage("ACCESS IMPOSSIBLE.");
      return 0;
    }
    fout=fgetstofopen("\0","w",idoutputs[lxy]);               /*OUTPUT FILE*/

    /*fgets(string,256,fin);*/                    /*INPUT APPELATION.*/
    /*errormessage(string);*/

    t0=clock();                                        /*CLOCK BEGIN.*/

    inputinit(fin,&nnode,&nelem,&nsect);             /*INPUT INITIAL.*/
    sprintf(string,"NODES=%d ELEMS=%d SECTS=%d",nnode,nelem,nsect);
    errormessage(string);
    /*if(fout!=NULL) fprintf(fout,"%s\n",string);*/

    msize=6*nnode; /*SIZE OF GLOBAL MATRIX.*/

    if(lxy==0) /* Initialise global matrix and vector only at LONG perild */
    {
      gmtx=(struct gcomponent *)          /*DIAGONALS OF GLOBAL MATRIX.*/
            malloc(msize*sizeof(struct gcomponent));
      gvct=(double *)malloc(msize*sizeof(double));      /*GLOBAL VECTOR*/
      if(gmtx==NULL || gvct==NULL) return 0;
      for(i=0;i<msize;i++)
      {
        (gmtx+i)->down=NULL;            /*GLOBAL MATRIX INITIALIZATION.*/
        *(gvct+i)=0.0;                  /*GLOBAL VECTOR INITIALIZATION.*/
      }
    } /* Initialise global matrix and vector only at LONG perild */
    else /* Initialise global vector at X, Y perild */
    {
      gvct=(double *)malloc(msize*sizeof(double));      /*GLOBAL VECTOR*/
      if(gvct==NULL) return 0;
      for(i=0;i<msize;i++)
      {
        *(gvct+i)=0.0;                  /*GLOBAL VECTOR INITIALIZATION.*/
      }
    }

    free(afs[lxy]->sects);
    free(afs[lxy]->nodes);
    free(afs[lxy]->ninit);
    free(afs[lxy]->elems);
    free(afs[lxy]->confs);
    free(afs[lxy]->ddisp);                     /*DISPLACEMENT:6 DIRECTIONS.*/
    free(afs[lxy]->melem);

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
    melem=(struct memoryelem *)malloc(nelem*sizeof(struct memoryelem));
    if(melem==NULL) return 0;

    afs[lxy]->sects=sects;
    afs[lxy]->nodes=nodes;
    afs[lxy]->ninit=ninit;
    afs[lxy]->elems=elems;
    afs[lxy]->confs=confs;
    afs[lxy]->ddisp=ddisp;
    afs[lxy]->melem=melem;

    inputtexttomemory(fin,afs[lxy]);        /*READY TO READ LONG REACTIONS.*/
    nnode=afs[lxy]->nnode;
    ninit=afs[lxy]->ninit;
    nelem=afs[lxy]->nelem;
    nsect=afs[lxy]->nsect;
    nreact=afs[lxy]->nreact;

    initialform(nodes,ddisp,nnode);           /*ASSEMBLAGE FORMATION.*/

    initialelem001(elems,melem,nelem);         /*ASSEMBLAGE ELEMENTS.*/

    dreact=(double *)malloc(nreact*sizeof(double));       /*REACTION.*/
    afs[lxy]->dreact=dreact;
    initialreact(fin,dreact,nreact);     /*ASSEMBLAGE LONG REACTIONS.*/

    GetAsyncKeyState(VK_LBUTTON);                   /*CLEAR KEY LEFT.*/
    GetAsyncKeyState(VK_RBUTTON);                  /*CLEAR KEY RIGHT.*/

    drawglobalaxis((wdraw.childs+1)->hdcC,(wdraw.childs+1)->vparam,
                   0,0,255);                      /*DRAW GLOBAL AXIS.*/

    errormessage("ARCLM001:ANALYSIS LINEAR.");
    memory1=availablephysicalmemory("REMAIN:");    /*MEMORY AVAILABLE*/

    if(lxy==0) /* Assemble global matrix only at LONG period */
    {
      for(i=1;i<=msize;i++)             /*GLOBAL MATRIX INITIALIZATION.*/
      {
        ginit.m=(unsigned short int)i;
        /*ginit.n=(unsigned short int)i;*/
        *(gmtx+(i-1))=ginit;
      }
      comps=msize; /*INITIAL COMPONENTS=DIAGONALS.*/

      laptime("ASSEMBLING GLOBAL MATRIX.",t0);

      for(i=1;i<=nelem;i++)               /*ASSEMBLAGE GLOBAL MATRIX.*/
      {
        elem=*(elems+i-1);                       /*READ ELEMENT DATA.*/

        inputnode(ddisp,elem.node[0]);                         /*HEAD*/
        inputnode(ddisp,elem.node[1]);                         /*TAIL*/

        elem.sect=(elems+i-1)->sect;             /*READ SECTION DATA.*/

        if((wdraw.childs+1)->hdcC!=NULL)     /*DRAW DEFORMED ELEMENT.*/
        {
          drawglobalwire((wdraw.childs+1)->hdcC,
                         (wdraw.childs+1)->vparam,
                         *afs[lxy],elem,255,255,255,
                                  255,255,255,0,ONSCREEN);
        }

        drccos=directioncosine(elem.node[0]->d[0],
                               elem.node[0]->d[1],
                               elem.node[0]->d[2],
                               elem.node[1]->d[0],
                               elem.node[1]->d[1],
                               elem.node[1]->d[2],
                               elem.cangle);               /*[DRCCOS]*/

        tmatrix=transmatrix(drccos);         /*TRANSFORMATION MATRIX.*/
        estiff=assememtx(elem);          /*ELASTIC MATRIX OF ELEMENT.*/
        estiff=modifyhinge(elem,estiff,afs[lxy]);    /*MODIFY MATRIX.*/
        estiff=transformation(estiff,tmatrix);       /*[K]=[Tt][k][T]*/

        assemgstiffness(gmtx,estiff,&elem);             /*ASSEMBLAGE.*/

        modifycmq(melem,&elem);
        assemcmq(elem,tmatrix,confs,gvct);    /*ASSEMBLAGE CMQ AS LOADS.*/

        for(ii=0;ii<=2;ii++) free(*(drccos+ii));
        free(drccos);
        for(ii=0;ii<=11;ii++) free(*(tmatrix+ii));
        free(tmatrix);
        for(ii=0;ii<=11;ii++) free(*(estiff+ii));
        free(estiff);
      }
      sprintf(string,"GLOBAL MATRIX %ld COMPS ASSEMBLED.",comps);
      laptime(string,t0);
    } /* Assemble global matrix only at LONG period */
    else /* Assemble global vector at X, Y period */
    {
      for(i=1;i<=nelem;i++)               /*ASSEMBLAGE GLOBAL MATRIX.*/
      {
        elem=*(elems+i-1);                       /*READ ELEMENT DATA.*/

        inputnode(ddisp,elem.node[0]);                         /*HEAD*/
        inputnode(ddisp,elem.node[1]);                         /*TAIL*/

        elem.sect=(elems+i-1)->sect;             /*READ SECTION DATA.*/

        if((wdraw.childs+1)->hdcC!=NULL)     /*DRAW DEFORMED ELEMENT.*/
        {
          drawglobalwire((wdraw.childs+1)->hdcC,
                         (wdraw.childs+1)->vparam,
                         *afs[lxy],elem,255,255,255,
                                  255,255,255,0,ONSCREEN);
        }

        drccos=directioncosine(elem.node[0]->d[0],
                               elem.node[0]->d[1],
                               elem.node[0]->d[2],
                               elem.node[1]->d[0],
                               elem.node[1]->d[1],
                               elem.node[1]->d[2],
                               elem.cangle);               /*[DRCCOS]*/

        tmatrix=transmatrix(drccos);         /*TRANSFORMATION MATRIX.*/
        //estiff=assememtx(elem);          /*ELASTIC MATRIX OF ELEMENT.*/
        //estiff=modifyhinge(elem,estiff);             /*MODIFY MATRIX.*/
        //estiff=transformation(estiff,tmatrix);       /*[K]=[Tt][k][T]*/

        //assemgstiffness(gmtx,estiff,&elem);             /*ASSEMBLAGE.*/

        modifycmq(melem,&elem);
        assemcmq(elem,tmatrix,confs,gvct);    /*ASSEMBLAGE CMQ AS LOADS.*/

        for(ii=0;ii<=2;ii++) free(*(drccos+ii));
        free(drccos);
        for(ii=0;ii<=11;ii++) free(*(tmatrix+ii));
        free(tmatrix);
        //for(ii=0;ii<=11;ii++) free(*(estiff+ii));
        //free(estiff);
      }
    } /* Assemble global vector at X, Y period */

    overlayhdc(*(wdraw.childs+1),SRCPAINT);          /*UPDATE DISPLAY.*/

    /*errormessage("ASSEMBLAGE GLOBAL VECTOR.");*/
    assemconf001(confs,gvct,1.0,nnode);               /*GLOBAL VECTOR.*/
    modifygivend(gmtx,gvct,confs,nnode);      /*0:LOAD 1:DISPLACEMENT.*/

    if(lxy==0) /* Decrease band width only at LONG period */
    {
      if(globalmessageflag==1 &&
         MessageBox(NULL,"Decrease Band Width ?","ARCLM001",MB_OKCANCEL)
         ==IDOK)
      {
        laptime("DECREASE BAND WIDTH.",t0);
        moff=(long int *)malloc(nnode*sizeof(long int));
        noff=(long int *)malloc(nnode*sizeof(long int));
        if(moff==NULL || noff==NULL) return 0;
        for(i=0;i<nnode;i++)
        {
          *(moff+i)=i;
          *(noff+i)=i;
        }
        decreaseband(moff,noff,afs[lxy]);              /*DECREASE BAND WIDTH.*/

        gmtx2=(struct gcomponent *)malloc(msize*sizeof(struct gcomponent));
        gvct2=(double *)malloc(msize*sizeof(double));
        confs2=(struct oconf *)malloc(msize*sizeof(struct oconf));
        if(gmtx2==NULL || gvct2==NULL || confs2==NULL) return 0;
        for(i=0;i<msize;i++)
        {
          ginit.m=(unsigned short int)(i+1);
          *(gmtx2+i)=ginit;

          *(gvct2+i)=0.0;

          (confs2+i)->iconf=0;
          (confs2+i)->value=0.0;
        }

        exchangelines(gmtx,gvct,confs,
                      gmtx2,gvct2,confs2,
                      moff,noff,nnode,0);
      }
      else
      {
        moff=NULL;
        noff=NULL;
        gmtx2=gmtx;
        gvct2=gvct;
        confs2=confs;
      }
    } /* Decrease band width only at LONG period */
    else /* period: X, Y */
    {
      gvct2=gvct;
    }

    laptime("CROUT LU DECOMPOSITION.",t0);
    /* croutludecomposition */
    if(lxy==0) /* LU decomposition only at LONG period */
    {
      nline=croutlu(gmtx2,
                    confs2,
                    6*nnode,
                    &determinant,&sign,
                    lu);
    } /* LU decomposition only at LONG period */

    forwardbackward(gmtx2,
                    gvct2,confs2,
                    6*nnode,
                    lu);
    
    //nline=croutludecomposition(gmtx2,
    //                           gvct2,confs2,
    //                           6*nnode,
    //                           &determinant,&sign);  /*[K]{dU}={dF}*/

    /* croutludecomposition */
    

    sprintf(string,"LOG{DETERMINANT}=%.5E SIGN=%.1f COMPS=%ld",
            determinant,sign,comps);
    errormessage(string);
    /*if(fout!=NULL) fprintf(fout,"%s\n",string);*/

    /*** 10.11.26 araki for instable termination ****************************/
      for(m=1;m<=msize;m++)
      {
        gread(gmtx2,m,m,&data);
        /*if(fout!=NULL)
              fprintf(fout,"NODE %ld. DIAGNAL %20.5f\n",
                      (af->nodes+int((m-1)/6))->code,data);*/
        if(data<=0.0)
        {
          sprintf(string,"INSTABLE TERMINATION AT NODE %ld.\n",
                  (afs[lxy]->nodes+int((m-1)/6))->code);
          errormessage(string);
        }
      }
    /************************************************************************/

    if(sign<=0.0)
    {
      if(noff==NULL)
      {
        sprintf(string,"INSTABLE TERMINATION AT NODE %ld.",
                (afs[lxy]->nodes+(int)(nline/6))->code);
      }
      else
      {
        sprintf(string,"INSTABLE TERMINATION AT NODE %ld.",
                (afs[lxy]->nodes+(*(noff+(int)(nline/6))))->code);
      }

      errormessage(" ");
      errormessage(string);
      if(fout!=NULL) fprintf(fout,"%s\n",string);

      laptime("\0",t0);

      fclose(fin);

      if(gmtx2==gmtx)
      {
        gfree(gmtx,nnode); /*FREE GLOBAL MATRIX.*/
        free(gvct);
        /*free(confs);*/
      }
      else
      {
        gfree(gmtx,nnode); /*FREE GLOBAL MATRIX.*/
        gfree(gmtx2,nnode);
        free(gvct);
        free(gvct2);
        /*free(confs);*/
        /*free(confs2);*/
      }

      memory2=availablephysicalmemory("REMAIN:");
      sprintf(string,"CONSUMPTION:%ld[BYTES]",(memory0-memory2));
      errormessage(string);

      return 1;
    }

    /*laptime("REEXCHANGE MATRIX.",t0);*/
    /*exchangelines(gmtx,gvct,confs,moff,noff,nnode,1);*/

    laptime("OUTPUT INTO FILE.",t0);

    errormessage("STRESS.");
    if(fout!=NULL)
    {
      fprintf(fout,"\n\n");
      fprintf(fout,"** FORCES OF MEMBER\n\n");
      fprintf(fout,"  NO   KT NODE         N        Q1        Q2");
      fprintf(fout,"        MT        M1        M2\n\n");
    }
    for(i=1;i<=nelem;i++)                   /*STRESS OUTPUT,UPDATE.*/
    {
      inputelem(elems,melem,i-1,&elem);

      inputnode(ddisp,elem.node[0]);
      inputnode(ddisp,elem.node[1]);

      elem.sect=(elems+i-1)->sect;             /*READ SECTION DATA.*/

      /*estress=elemstress001(&elem,gvct,melem);*/
      estress=elemstress002(&elem,gvct2,melem,moff,afs[lxy]);

      outputstress001(elem,estress,fout);
      free(estress);
    }
    if(fout!=NULL) fprintf(fout,"\n\n");

    errormessage("DISPLACEMENT.");
    /*outputdisp001(gvct,fout,nnode,nodes);*/           /*DISPLACEMENT.*/
    outputdisp002(gvct2,fout,nnode,nodes,moff);     /*DISPLACEMENT.*/
    if(fout!=NULL) fprintf(fout,"\n\n");
    /*while(!GetAsyncKeyState(VK_LBUTTON))
    ;*/                                   /*LEFT CLICK TO CONTINUE.*/

    errormessage("REACTION.");
    if(fout!=NULL)
    {
      fprintf(fout,"** REACTION\n\n");
      fprintf(fout,"  NO  DIRECTION              R    NC\n\n");
    }
    /*outputreaction001(gmtx,gvct,nodes,confs,dreact,fout,nnode);*/
    outputreaction002(gmtx2,gvct2,nodes,confs2,dreact,fout,nnode,moff);

    fclose(fin);
    fclose(fout);
    errormessage("FILES CLOSED.");

    /*updateform(ddisp,gvct,nnode);*/               /*FORMATION UPDATE.*/
    updateform2(ddisp,gvct2,nnode,moff);        /*FORMATION UPDATE.*/

    laptime("\0",t0);

    memory2=availablephysicalmemory(NULL);
    sprintf(string,"CONSUMPTION:%ld[BYTES]",(memory1-memory2));
    errormessage(string);

    if((wdraw.childs+1)->hdcC!=NULL &&
       melem!=NULL && ddisp!=NULL)               /*DRAW LAST FRAME.*/
    {
      for(i=1;i<=nelem;i++)
      {
        elem=*(elems+i-1);                     /*READ ELEMENT DATA.*/

        inputnode(ddisp,elem.node[0]);
        inputnode(ddisp,elem.node[1]);

        drawglobalwire((wdraw.childs+1)->hdcC,
                       (wdraw.childs+1)->vparam,
                       *afs[lxy],elem,255,255,255,
                                255,255,255,0,ONSCREEN);
      }
      overlayhdc(*(wdraw.childs+1),SRCPAINT);     /*UPDATE DISPLAY.*/
    }

    if(lxy==2) /* Free global matrix only at Y period */
    {
      if(gmtx2==gmtx)
      {
        gfree(gmtx,nnode); /*FREE GLOBAL MATRIX.*/
        /*free(gvct);*/
        /*free(confs);*/
      }
      else
      {
        gfree(gmtx,nnode); /*FREE GLOBAL MATRIX.*/
        gfree(gmtx2,nnode);
        /*free(gvct);*/
        /*free(gvct2);*/
        /*free(confs);*/
        /*free(confs2);*/
      }
    } /* Free global matrix only at Y period */

    afs[lxy]->nlaps=1;
    /*af->eigenvec=(double **)malloc(1*sizeof(double *));*/
    /**((af->eigenvec)+0)=gvct2;*/

    errormessage(" ");
    errormessage("COMPLETED.");

    memory2=availablephysicalmemory("REMAIN:");
    sprintf(string,"CONSUMPTION:%ld[BYTES]",(memory0-memory2));
    errormessage(string);
  }

  return 0;
}/*arclm001_lxy*/
/***/

int arclm401(struct arclmframe *af,int idinput)
{
  /* Analysis for Prestressed Structure */

  DWORD memory0,memory1,memory2;

  FILE *fin,*fout,*fins;                             /*FILE 8 BYTES*/
  FILE *fmtx1, *fmtx2, *fmtx3;
  double *ddisp,*dreact;
  struct memoryelem *melem;
  char dir[]=DIRECTORY;                            /*DATA DIRECTORY*/
  char s[120],string[400]/*,fname[256]*/;
  int i,ii,jj;
  int nnode,nelem,nsect,nreact;
  long int msize;
  long int loff;
  struct gcomponent ginit={0,0,0.0,NULL};
  struct gcomponent *kmtx,*gmtx,*k,*g,*p;                    /*GLOBAL MATRIX*/
  double *gvct;                                     /*GLOBAL VECTOR*/
  double **drccos,**tmatrix,**estiff,*estress;
  double determinant,sign;
  clock_t t0;
  char **data;
  int n,times,nlap;

  struct osect *sects;
  struct onode *nodes,*ninit;
  struct owire elem;
  struct owire *elems;
  struct oconf *confs;

  double gg;

  memory0=availablephysicalmemory("INITIAL:");   /*MEMORY AVAILABLE*/

  fins=fgetstofopen(dir,"r",idinput);              /*OPEN INPUT FILE*/
  if(fins==NULL)
  {
    errormessage("ACCESS IMPOSSIBLE.");
    return 0;
  }
  //fpre=fgetstofopen(dir,"r",idpre);
  fout=fgetstofopen("\0","w",ID_OUTPUTFILE);          /*OUTPUT FILE*/

  fmtx1=fopen("matrix401_1.txt","w");
  fmtx2=fopen("matrix401_2.txt","w");
  fmtx3=fopen("matrix401_3.txt","w");

  /*fgets(string,256,fin);*/                    /*INPUT APPELATION.*/
  /*errormessage(string);*/

  t0=clock();                                        /*CLOCK BEGIN.*/

  while (1)
  {
      data=fgetsbrk(fins,&n);
      if (n==0) break;

      fin=fopen(*(data+0),"r");
      /* if(fin==NULL) */
      /* { */
      /*   errormessage("ACCESS IMPOSSIBLE."); */
      /*   return 0; */
      /* } */
      times=(int)strtol(*(data+1),NULL,10);

      sprintf(string,"INPUT: %s TIMES: %d",
              *(data+0),times);
      errormessage(string);

      inputinit(fin,&nnode,&nelem,&nsect);             /*INPUT INITIAL.*/
      sprintf(string,"NODES=%d ELEMS=%d SECTS=%d",nnode,nelem,nsect);
      errormessage(string);
      //if(fout!=NULL) fprintf(fout,"%s\n",string);

      for (nlap=1;nlap<=times;nlap++)
      {
          msize=6*nnode;                           /*SIZE OF GLOBAL MATRIX.*/

          kmtx=(struct gcomponent *)          /*DIAGONALS OF GLOBAL MATRIX.*/
                malloc(msize*sizeof(struct gcomponent));
          gmtx=(struct gcomponent *)          /*DIAGONALS OF GLOBAL MATRIX.*/
                malloc(msize*sizeof(struct gcomponent));
          gvct=(double *)malloc(msize*sizeof(double));      /*GLOBAL VECTOR*/
          if(kmtx==NULL || gmtx==NULL || gvct==NULL) return 0;
          for(i=0;i<=msize-1;i++)
          {
            (kmtx+i)->down=NULL;            /*GLOBAL MATRIX INITIALIZATION.*/
            (gmtx+i)->down=NULL;            /*GLOBAL MATRIX INITIALIZATION.*/
            *(gvct+i)=0.0;                  /*GLOBAL VECTOR INITIALIZATION.*/
          }

          //errormessage("Read Prestress.");
          //frameoutputtomemory(fpre,&arc);
          //errormessage("Read Prestress Succeeded.");

          free(af->sects);
          free(af->nodes);
          free(af->ninit);
          free(af->elems);
          free(af->confs);
          //free(af->ddisp);
          //free(af->melem);

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
          melem=(struct memoryelem *)
                malloc(nelem*sizeof(struct memoryelem));
          if(melem==NULL) return 0;

          af->sects=sects;
          af->nodes=nodes;
          af->ninit=ninit;
          af->elems=elems;
          af->confs=confs;
          //af->ddisp=ddisp;                   /*DISPLACEMENT:6 DIRECTIONS.*/
          //af->melem=melem;
          melem=af->melem;
          ddisp=af->ddisp;

          inputtexttomemory(fin,af);      /*READY TO READ LONG REACTIONS.*/
          nnode=af->nnode;
          nelem=af->nelem;
          nsect=af->nsect;
          nreact=af->nreact;

          /* initialform(nodes,ddisp,nnode);         #<{(|ASSEMBLAGE FORMATION.|)}># */
          initialelem401(elems,melem,nelem);          /*ASSEMBLAGE ELEMENTS.*/

          dreact=(double *)malloc(nreact*sizeof(double));     /*REACTION.*/
          af->dreact=dreact;
          initialreact(fin,dreact,nreact);   /*ASSEMBLAGE LONG REACTIONS.*/

          GetAsyncKeyState(VK_LBUTTON);                   /*CLEAR KEY LEFT.*/
          GetAsyncKeyState(VK_RBUTTON);                  /*CLEAR KEY RIGHT.*/

          drawglobalaxis((wdraw.childs+1)->hdcC,(wdraw.childs+1)->vparam,
                         0,0,255);                      /*DRAW GLOBAL AXIS.*/

          errormessage("ARCLM401:ANALYSIS LINEAR.");
          memory1=availablephysicalmemory("REMAIN:");    /*MEMORY AVAILABLE*/

          for(i=1;i<=msize;i++)           /*GLOBAL MATRIX INITIALIZATION.*/
          {
            /* only for incremental analysis? */
            k=(kmtx+(i-1))->down; /*NEXT OF DIAGONAL.*/
            g=(gmtx+(i-1))->down; /*NEXT OF DIAGONAL.*/
            while(k!=NULL) /*CLEAR ROW.*/
            {
              p=k;
              k=k->down;
              free(p);
            }
            while(g!=NULL) /*CLEAR ROW.*/
            {
              p=g;
              g=g->down;
              free(p);
            }

            ginit.m=(unsigned short int)i;
            /*ginit.n=(unsigned short int)i;*/

            *(kmtx+(i-1))=ginit;
            *(gmtx+(i-1))=ginit;
          }
          comps=msize; /*INITIAL COMPONENTS=DIAGONALS.*/

          laptime("ASSEMBLING GLOBAL MATRIX.",t0);

          estress=(double *)malloc(12*sizeof(double));

          for(i=1;i<=nelem;i++)               /*ASSEMBLAGE GLOBAL MATRIX.*/
          {
            /* melem or af->melem? */
            inputelem(elems,melem,i-1,&elem);        /*READ ELEMENT DATA.*/

            /* only for incremental analysis? */
            for(ii=0;ii<=1;ii++)
            {
              for(jj=0;jj<=5;jj++)
              {
                (elems+i-1)->iconf[ii][jj]=elem.iconf[ii][jj];
              }
            }

            /* taken from bclng001 */
            //for(ii=0;ii<=1;ii++) /*INITIALIZE COORDINATION.*/
            //{
              //loff=elem.node[ii]->loff;
              //for(jj=0;jj<3;jj++)
              //{
                //elem.node[ii]->d[jj]=(ninit+loff)->d[jj];
              //}
            //}
            inputnode(ddisp,elem.node[0]);                         /*HEAD*/
            inputnode(ddisp,elem.node[1]);                         /*TAIL*/

            elem.sect=(elems+i-1)->sect;             /*READ SECTION DATA.*/

            if((wdraw.childs+1)->hdcC!=NULL)     /*DRAW DEFORMED ELEMENT.*/
            {
              drawglobalwire((wdraw.childs+1)->hdcC,
                             (wdraw.childs+1)->vparam,
                             *af,elem,255,255,255,
                                      255,255,255,0,ONSCREEN);
            }

            drccos=directioncosine(elem.node[0]->d[0],
                                   elem.node[0]->d[1],
                                   elem.node[0]->d[2],
                                   elem.node[1]->d[0],
                                   elem.node[1]->d[1],
                                   elem.node[1]->d[2],
                                   elem.cangle);               /*[DRCCOS]*/

            tmatrix=transmatrix(drccos);         /*TRANSFORMATION MATRIX.*/
            estiff=assememtx(elem);          /*ELASTIC MATRIX OF ELEMENT.*/
            //errormessage("Modify hinge for e");
            estiff=modifyhinge(elem,estiff,af);          /*MODIFY MATRIX.*/
            estiff=transformation(estiff,tmatrix);       /*[K]=[Tt][k][T]*/

            assemgstiffness(kmtx,estiff,&elem);             /*ASSEMBLAGE.*/

            for(ii=0;ii<=11;ii++) free(*(estiff+ii));
            free(estiff);

            for(ii=0;ii<=1;ii++)                                /*STRESSES.*/
            {
              for(jj=0;jj<6;jj++) *(estress+6*ii+jj)=elem.stress[ii][jj];
            }

            // add GG to gmtx

            estiff=assemgmtx(elem,estress);
            //errormessage("Modify hinge for g");
            estiff=modifyhinge401(elem,estiff);               /*MODIFY MATRIX.*/
            estiff=transformation(estiff,tmatrix);         /*[K]=[Tt][k][T]*/
            //for(ii=0;ii<12;ii++)
            //{
              //for(jj=0;jj<12;jj++) *(*(estiff+ii)+jj)=-*(*(estiff+ii)+jj);
            //}

            assemgstiffness(gmtx,estiff,&elem);     /*ASSEMBLAGE GEOMETRIC.*/

            // modifycmq?
            //modifycmq(melem,&elem);
            //assemcmq(elem,tmatrix,confs,gvct);    /*ASSEMBLAGE CMQ AS LOADS.*/

            for(ii=0;ii<=2;ii++) free(*(drccos+ii));
            free(drccos);
            for(ii=0;ii<=11;ii++) free(*(tmatrix+ii));
            free(tmatrix);
            for(ii=0;ii<=11;ii++) free(*(estiff+ii));
            free(estiff);
          }
          sprintf(string,"GLOBAL MATRIX %ld COMPS ASSEMBLED.",comps);
          laptime(string,t0);

          if(fmtx1!=NULL)
          {
            outputmtx_conf(fmtx1,kmtx,confs,msize);
            fclose(fmtx1);
          }
          if(fmtx2!=NULL)
          {
            outputmtx_conf(fmtx2,gmtx,confs,msize);
            fclose(fmtx2);
          }

          gcomponentadd(gmtx,kmtx,msize);
          //gcompsub(gmtx,kmtx,msize);

          if(fmtx3!=NULL)
          {
            outputmtx_conf(fmtx3,gmtx,confs,msize);
            fclose(fmtx3);
          }

          overlayhdc(*(wdraw.childs+1),SRCPAINT);       /*UPDATE DISPLAY.*/

          /*errormessage("ASSEMBLAGE GLOBAL VECTOR.");*/
          assemconf(confs,gvct,1.0,nnode);           /*GLOBAL VECTOR.*/
          modifygivend(gmtx,gvct,confs,nnode);   /*0:LOAD 1:DISPLACEMENT.*/
          //for(i=0;i<msize;i++)
          //{
            //if(*(gvct+i)!=0.0)
            //{
              //sprintf(string,"%ld %12.5f",i,*(gvct+i));
              //errormessage(string);
            //}
          //}
          laptime("CROUT LU DECOMPOSITION.",t0);
          croutludecomposition(gmtx,
                               gvct,confs,
                               6*nnode,
                               &determinant,&sign);        /*[K]{dU}={dF}*/

          sprintf(string,"LOG{DETERMINANT}=%.5E SIGN=%.1f COMPS=%ld",
                  determinant,sign,comps);
          errormessage(string);

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
              }
            }
            errormessage(" ");
            errormessage(string);
            if(fout!=NULL) fprintf(fout,"%s\n",string);

            laptime("\0",t0);

            fclose(fin);
            //fclose(fout);
            //fclose(ffig);
            /*fclose(felem);*/
            /*fclose(fdisp);*/
            /*fclose(freact);*/

            gfree(gmtx,nnode); /*FREE GLOBAL MATRIX.*/
            free(gvct);
            /*free(confs);*/

            memory2=availablephysicalmemory("REMAIN:");
            sprintf(string,"CONSUMPTION:%ld[BYTES]",(memory0-memory2));
            errormessage(string);

            return 1;
          }

          laptime("OUTPUT INTO FILE.",t0);

          //if(fout!=NULL) fprintf(fout,"\"DISPLACEMENT\"\n");
          //outputdisp(gvct,fout,nnode,nodes);  /*DISPLACEMENT(this don't include those by prestress.*/

          errormessage("STRESS.");
          if(fout!=NULL)
          {
            fprintf(fout,"\n\n");
            fprintf(fout,"** FORCES OF MEMBER\n\n");
            fprintf(fout,"  NO   KT NODE         N        Q1        Q2");
            fprintf(fout,"        MT        M1        M2\n\n");
          }
          for(i=1;i<=nelem;i++)                   /*STRESS OUTPUT,UPDATE.*/
          {
            inputelem(elems,melem,i-1,&elem);

            inputnode(ddisp,elem.node[0]);
            inputnode(ddisp,elem.node[1]);

            elem.sect=(elems+i-1)->sect;             /*READ SECTION DATA.*/

            estress=elemstress001(&elem,gvct,melem,af);

            outputstress001(elem,estress,fout);
            free(estress);
          }
          if(fout!=NULL) fprintf(fout,"\n\n");

          errormessage("DISPLACEMENT.");
          outputdisp001(gvct,fout,nnode,nodes);           /*DISPLACEMENT.*/
          if(fout!=NULL) fprintf(fout,"\n\n");

          errormessage("REACTION.");
          if(fout!=NULL)
          {
            fprintf(fout,"** REACTION\n\n");
            fprintf(fout,"  NO  DIRECTION              R    NC\n\n");
          }
          outputreaction001(gmtx,gvct,nodes,confs,dreact,fout,nnode);

          updateform(ddisp,gvct,nnode);               /*FORMATION UPDATE.*/

          laptime("\0",t0);
      }
  }

  memory2=availablephysicalmemory(NULL);
  sprintf(string,"CONSUMPTION:%ld[BYTES]",(memory1-memory2));
  errormessage(string);

  if((wdraw.childs+1)->hdcC!=NULL &&
     melem!=NULL && ddisp!=NULL)                 /*DRAW LAST FRAME.*/
  {
    for(i=1;i<=nelem;i++)
    {
      inputelem(elems,melem,i-1,&elem);

      inputnode(ddisp,elem.node[0]);
      inputnode(ddisp,elem.node[1]);

      drawglobalwire((wdraw.childs+1)->hdcC,
                     (wdraw.childs+1)->vparam,
                     *af,elem,255,255,255,
                              255,255,255,0,ONSCREEN);
    }
    overlayhdc(*(wdraw.childs+1),SRCPAINT);       /*UPDATE DISPLAY.*/
  }

  fclose(fin);
  /*fclose(felem);*/
  /*fclose(fdisp);*/
  /*fclose(freact);*/

  gfree(gmtx,nnode); /*FREE GLOBAL MATRIX.*/
  /*free(gvct);*/
  /*free(confs);*/

  af->nlaps=1;
  //af->eigenvec=(double **)malloc(1*sizeof(double *));
  //*((af->eigenvec)+0)=gvct;

  errormessage(" ");
  errormessage("COMPLETED.");

  if(fout!=NULL) fclose(fout);

  memory2=availablephysicalmemory("REMAIN:");
  sprintf(string,"CONSUMPTION:%ld[BYTES]",(memory0-memory2));
  errormessage(string);
  errormessage(" ");

  return 0;
}/*arclm401*/

int arclmtest(struct arclmframe *af,int idinput,int idoutput)
{
  FILE *fin,*fout;                                   /*FILE 8 BYTES*/
  double *ddisp,*dreact;
  struct memoryelem *melem;
  char dir[]=DIRECTORY;                            /*DATA DIRECTORY*/
  char string[400];
  int i,ii,jj;
  int nnode,nelem,nsect,nreact;
  long int msize;
  struct gcomponent ginit={0,0,0.0,NULL};
  struct gcomponent *gmtx,*mmtx;                    /*GLOBAL MATRIX*/
  double *gvct;                                     /*GLOBAL VECTOR*/
  double **drccos,**tmatrix,**estiff,*estress;
  double determinant,sign,data,value;

  struct osect *sects;
  struct onode *nodes;
  struct owire elem;
  struct owire *elems;
  struct oconf *confs;

  /*fin=fgetstofopen(dir,"r",idinput);*/              /*OPEN INPUT FILE*/
  fin=fopen("c:\\cdocs\\data\\model01.inp","r");
  /*fout=fgetstofopen("\0","w",idoutput);*/               /*OUTPUT FILE*/
  fout=fopen("model01.otp","w");
  if(fin==NULL || fout==NULL) return 0;

  fprintf(fout,"INPUT :MODEL01.INP\n");
  fprintf(fout,"OUTPUT:MODEL01.OTP\n");

  inputinit(fin,&nnode,&nelem,&nsect);             /*INPUT INITIAL.*/
  fprintf(fout,"êﬂì_êî=%d óvëfêî=%d ífñ É^ÉCÉvêî=%d\n",nnode,nelem,nsect);
  fprintf(fout,"\n");

  msize=6*nnode; /*SIZE OF GLOBAL MATRIX.*/

  gmtx=(struct gcomponent *)          /*DIAGONALS OF GLOBAL MATRIX.*/
        malloc(msize*sizeof(struct gcomponent));
  gvct=(double *)malloc(msize*sizeof(double));      /*GLOBAL VECTOR*/
  if(gmtx==NULL || gvct==NULL) return 0;
  for(i=0;i<=msize-1;i++)
  {
    (gmtx+i)->down=NULL;            /*GLOBAL MATRIX INITIALIZATION.*/
    *(gvct+i)=0.0;                  /*GLOBAL VECTOR INITIALIZATION.*/
  }

  mmtx=(struct gcomponent *)                 /*MEMORY GLOBAL MATRIX*/
        malloc(msize*sizeof(struct gcomponent));
  if(mmtx==NULL) return 0;
  for(i=0;i<=msize-1;i++) (mmtx+i)->down=NULL;

  free(af->sects);
  free(af->nodes);
  free(af->elems);
  free(af->confs);
  free(af->ddisp);                     /*DISPLACEMENT:6 DIRECTIONS.*/
  free(af->melem);

  sects=(struct osect *)malloc(nsect*sizeof(struct osect));
  if(sects==NULL) return 0;
  nodes=(struct onode *)malloc(nnode*sizeof(struct onode));
  if(nodes==NULL) return 0;
  elems=(struct owire *)malloc(nelem*sizeof(struct owire));
  if(elems==NULL) return 0;
  confs=(struct oconf *)malloc(msize*sizeof(struct oconf));
  if(confs==NULL) return 0;
  ddisp=(double *)malloc(6*nnode*sizeof(double));
  if(ddisp==NULL) return 0;
  melem=(struct memoryelem *)malloc(nelem*sizeof(struct memoryelem));
  if(melem==NULL) return 0;

  af->sects=sects;
  af->nodes=nodes;
  af->elems=elems;
  af->confs=confs;
  af->ddisp=ddisp;
  af->melem=melem;

  inputtexttomemory(fin,af);        /*READY TO READ LONG REACTIONS.*/
  nnode=af->nnode;
  nelem=af->nelem;
  nsect=af->nsect;
  nreact=af->nreact;

  initialform(nodes,ddisp,nnode);           /*ASSEMBLAGE FORMATION.*/

  initialelem001(elems,melem,nelem);         /*ASSEMBLAGE ELEMENTS.*/

  dreact=(double *)malloc(nreact*sizeof(double));       /*REACTION.*/
  af->dreact=dreact;
  initialreact(fin,dreact,nreact);     /*ASSEMBLAGE LONG REACTIONS.*/

  for(i=1;i<=msize;i++)             /*GLOBAL MATRIX INITIALIZATION.*/
  {
    ginit.m=(unsigned short int)i;
    *(gmtx+(i-1))=ginit;
    *(mmtx+(i-1))=ginit;
  }
  comps=msize; /*INITIAL COMPONENTS=DIAGONALS.*/

  for(i=1;i<=nelem;i++)               /*ASSEMBLAGE GLOBAL MATRIX.*/
  {
    elem=*(elems+i-1);                       /*READ ELEMENT DATA.*/

    inputnode(ddisp,elem.node[0]);                         /*HEAD*/
    inputnode(ddisp,elem.node[1]);                         /*TAIL*/
    fprintf(fout,"óvëf %4d éní[=%3d ç¿ïW={%6.3f %6.3f %6.3f}\n",
            elem.code,elem.node[0]->code,elem.node[0]->d[0]
                                        ,elem.node[0]->d[1]
                                        ,elem.node[0]->d[2]);
    fprintf(fout,"          èIí[=%3d ç¿ïW={%6.3f %6.3f %6.3f}\n",
            elem.node[1]->code,elem.node[1]->d[0]
                              ,elem.node[1]->d[1]
                              ,elem.node[1]->d[2]);

    elem.sect=(elems+i-1)->sect;             /*READ SECTION DATA.*/
    fprintf(fout,"          ífñ =%3d\n",elem.sect->code);
    fprintf(fout,"          Çk=%.3f",sqrt((elem.node[1]->d[0]-elem.node[0]->d[0])
                                     *(elem.node[1]->d[0]-elem.node[0]->d[0])
                                     +(elem.node[1]->d[1]-elem.node[0]->d[1])
                                     *(elem.node[1]->d[1]-elem.node[0]->d[1])
                                     +(elem.node[1]->d[2]-elem.node[0]->d[2])
                                     *(elem.node[1]->d[2]-elem.node[0]->d[2])));
    fprintf(fout," Çd=%.5E"    ,elem.sect->E);
    fprintf(fout," ÉÀ=%.5f"    ,elem.sect->poi);
    fprintf(fout," Ç`=%.4f"    ,elem.sect->area);
    fprintf(fout," Çhxx=%.8f"  ,elem.sect->Ixx);
    fprintf(fout," Çhyy=%.8f"  ,elem.sect->Iyy);
    fprintf(fout," Çizz=%.8f\n",elem.sect->Jzz);

    drccos=directioncosine(elem.node[0]->d[0],
                           elem.node[0]->d[1],
                           elem.node[0]->d[2],
                           elem.node[1]->d[0],
                           elem.node[1]->d[1],
                           elem.node[1]->d[2],
                           elem.cangle);               /*[DRCCOS]*/
    fprintf(fout,"ï˚å¸ó]å∑ [Çc] =");
    for(ii=0;ii<3;ii++)
    {
      if(ii>0) fprintf(fout,"               ");
      for(jj=0;jj<3;jj++)
      {
        fprintf(fout," %6.3f",*((*(drccos+ii))+jj));
      }
      fprintf(fout,"\n");
    }

    tmatrix=transmatrix(drccos);         /*TRANSFORMATION MATRIX.*/
    fprintf(fout,"ç¿ïWïœä∑çsóÒ [Çs] =");
    for(ii=0;ii<12;ii++)
    {
      if(ii>0) fprintf(fout,"                   ");
      for(jj=0;jj<12;jj++)
      {
        fprintf(fout," %6.3f",*((*(tmatrix+ii))+jj));
      }
      fprintf(fout,"\n");
    }

    estiff=assememtx(elem);          /*ELASTIC MATRIX OF ELEMENT.*/
    fprintf(fout,"ïîçﬁçÑê´çsóÒ [Çãe] =");
    for(ii=0;ii<12;ii++)
    {
      if(ii>0) fprintf(fout,"                    ");
      for(jj=0;jj<12;jj++)
      {
        fprintf(fout," %5.0f",*((*(estiff+ii))+jj));
      }
      fprintf(fout,"\n");
    }

    estiff=modifyhinge(elem,estiff,af);          /*MODIFY MATRIX.*/

    estiff=transformation(estiff,tmatrix);       /*[K]=[Tt][k][T]*/
    fprintf(fout,"ç¿ïWïœä∑å„ÇÃïîçﬁçÑê´çsóÒ [Çje] =");
    for(ii=0;ii<12;ii++)
    {
      if(ii>0) fprintf(fout,"                                ");
      for(jj=0;jj<12;jj++)
      {
        fprintf(fout," %5.0f",*((*(estiff+ii))+jj));
      }
      fprintf(fout,"\n");
    }

    assemgstiffness(gmtx,estiff,&elem);             /*ASSEMBLAGE.*/
    fprintf(fout,"ëSëÃçÑê´çsóÒ [Çjg] =");
    for(ii=0;ii<6*nnode;ii++)
    {
      if(ii>0) fprintf(fout,"                    ");
      for(jj=0;jj<3*nnode;jj++)
      {
        gread(gmtx,ii+1,jj+1,&data);
        fprintf(fout," %5.0f",data);
      }
      fprintf(fout,"\n");
    }
    fprintf(fout,"å„îº\n");
    for(ii=0;ii<6*nnode;ii++)
    {
      for(jj=3*nnode;jj<6*nnode;jj++)
      {
        gread(gmtx,ii+1,jj+1,&data);
        fprintf(fout," %5.0f",data);
      }
      fprintf(fout,"\n");
    }

    modifycmq(melem,&elem);
    assemcmq(elem,tmatrix,confs,gvct); /*ASSEMBLAGE CMQ AS LOADS.*/

    for(ii=0;ii<=2;ii++) free(*(drccos+ii));
    free(drccos);
    for(ii=0;ii<=11;ii++) free(*(tmatrix+ii));
    free(tmatrix);
    for(ii=0;ii<=11;ii++) free(*(estiff+ii));
    free(estiff);
  }

  fprintf(fout,"ëSëÃçÑê´çsóÒç≈èIíl\n");
  fprintf(fout,"[Çjg] =");
  for(ii=0;ii<6*nnode;ii++)
  {
    if(ii>0) fprintf(fout,"       ");
    for(jj=0;jj<3*nnode;jj++)
    {
      gread(gmtx,ii+1,jj+1,&data);
      fprintf(fout," %5.0f",data);
    }
    fprintf(fout,"\n");
  }
  fprintf(fout,"å„îº\n");
  for(ii=0;ii<6*nnode;ii++)
  {
    for(jj=3*nnode;jj<6*nnode;jj++)
    {
      gread(gmtx,ii+1,jj+1,&data);
      fprintf(fout," %5.0f",data);
    }
    fprintf(fout,"\n");
  }

  assemconf001(confs,gvct,1.0,nnode);            /*GLOBAL VECTOR.*/
  modifygivend(gmtx,gvct,confs,nnode);   /*0:LOAD 1:DISPLACEMENT.*/

  fprintf(fout,"çSë©èåèÇ…ÇÊÇËèCê≥Ç≥ÇÍÇΩëSëÃçÑê´çsóÒ\n");
  fprintf(fout,"{Çe} = [Çjg]{Çt}\n");
  for(ii=0;ii<6*nnode;ii++)
  {
    if((confs+ii)->iconf==1) fprintf(fout,"   R%d%d",ii/6+1,ii%6+1);
    else                     fprintf(fout," %5.1f",*(gvct+ii));

    if(ii==0) fprintf(fout," =");
    else      fprintf(fout,"  ");

    for(jj=0;jj<3*nnode;jj++)
    {
      gread(gmtx,ii+1,jj+1,&data);
      fprintf(fout," %5.0f",data);
    }
    fprintf(fout,"\n");
  }
  fprintf(fout,"å„îº\n");
  for(ii=0;ii<6*nnode;ii++)
  {
    for(jj=3*nnode;jj<6*nnode;jj++)
    {
      gread(gmtx,ii+1,jj+1,&data);
      fprintf(fout," %5.0f",data);
    }

    if((confs+ii)->iconf==1) fprintf(fout,"  %.3f",(confs+ii)->value);
    else                     fprintf(fout,"   U%d%d",ii/6+1,ii%6+1);

    if((confs+ii)->iconf==1) fprintf(fout,"  çSë©");

    fprintf(fout,"\n");
  }

  /*MEMORY GLOBAL MATRIX*/
  for(ii=0;ii<6*nnode;ii++)
  {
    for(jj=0;jj<6*nnode;jj++)
    {
      gread(gmtx,ii+1,jj+1,&data);
      gwrite(mmtx,ii+1,jj+1,data);
    }
  }

  croutludecomposition(gmtx,
                       gvct,confs,
                       6*nnode,
                       &determinant,&sign);        /*[K]{dU}={dF}*/

  sprintf(string,"LOG{DETERMINANT}=%.5E SIGN=%.1f COMPS=%ld",
          determinant,sign,comps);
  errormessage(string);

  if(sign<=0.0)
  {
    if(fout!=NULL) fprintf(fout,"INSTABLE TERMINATION.\n");

    fclose(fin);

    gfree(gmtx,nnode); /*FREE GLOBAL MATRIX.*/
    free(gvct);

    return 1;
  }

  fprintf(fout,"ìæÇÁÇÍÇΩâÇÃåüéZ\n");
  for(ii=0;ii<6*nnode;ii++)
  {
    if(ii==0) fprintf(fout,"[Çjg]{Çt} =");
    else      fprintf(fout,"           ");

    for(jj=0;jj<3*nnode;jj++)
    {
      gread(mmtx,ii+1,jj+1,&data);
      fprintf(fout," %5.0f",data);
    }
    fprintf(fout,"\n");
  }
  fprintf(fout,"å„îº\n");
  for(ii=0;ii<6*nnode;ii++)
  {
    for(jj=3*nnode;jj<6*nnode;jj++)
    {
      gread(mmtx,ii+1,jj+1,&data);
      fprintf(fout," %5.0f",data);
    }

    fprintf(fout,"   %6.3f",*(gvct+ii));

    fprintf(fout,"\n");
  }

  fprintf(fout,"åÎç∑\n");
  for(ii=0;ii<6*nnode;ii++)
  {
    if(ii==0) fprintf(fout,"[Çjg]{Çt} =");
    else      fprintf(fout,"           ");

    value=0.0;

    for(jj=0;jj<6*nnode;jj++)
    {
      gread(mmtx,ii+1,jj+1,&data);
      value+=data*(*(gvct+jj));
    }
    fprintf(fout," %4.1f",value);

    if(ii==0) fprintf(fout," {Çe} =");
    else      fprintf(fout,"       ");

    if((confs+ii)->iconf==1) fprintf(fout,"  R%d%d",ii/6+1,ii%6+1);
    else                     fprintf(fout," %4.1f",(confs+ii)->value);

    if(ii==0) fprintf(fout," åÎç∑ [Çjg]{Çt}-{Çe} =");
    else      fprintf(fout,"                      ");

    if((confs+ii)->iconf==1) fprintf(fout,"  çSë©");
    else                     fprintf(fout," %10.3E",((confs+ii)->value)-value);

    fprintf(fout,"\n");
  }

  fprintf(fout,"ïîçﬁâûóÕ\n");
  /*if(fout!=NULL)
  {
    fprintf(fout,"\n\n");
    fprintf(fout,"** FORCES OF MEMBER\n\n");
    fprintf(fout,"  NO   KT NODE         N        Qx        Qy");
    fprintf(fout,"        Mz        Mx        My\n\n");
  }*/
  for(i=1;i<=nelem;i++)                   /*STRESS OUTPUT,UPDATE.*/
  {
    inputelem(elems,melem,i-1,&elem);

    inputnode(ddisp,elem.node[0]);
    inputnode(ddisp,elem.node[1]);

    elem.sect=(elems+i-1)->sect;             /*READ SECTION DATA.*/

    fprintf(fout,"óvëf %4d éní[=%3d èIí[=%3d\n",
            elem.code,elem.node[0]->code,elem.node[1]->code);
    estress=elemstress003(fout,&elem,gvct,melem,af);

    /*outputstress001(elem,estress,fout);*/
    free(estress);
  }
  if(fout!=NULL) fprintf(fout,"\n\n");

  errormessage("DISPLACEMENT.");
  outputdisp001(gvct,fout,nnode,nodes);           /*DISPLACEMENT.*/
  if(fout!=NULL) fprintf(fout,"\n\n");
  /*while(!GetAsyncKeyState(VK_LBUTTON))
  ;*/                                   /*LEFT CLICK TO CONTINUE.*/

  errormessage("REACTION.");
  if(fout!=NULL)
  {
    fprintf(fout,"** REACTION\n\n");
    fprintf(fout,"  NO  DIRECTION              R    NC\n\n");
  }
  outputreaction001(gmtx,gvct,nodes,confs,dreact,fout,nnode);

  fclose(fin);
  fclose(fout);

  updateform(ddisp,gvct,nnode);               /*FORMATION UPDATE.*/

  gfree(gmtx,nnode); /*FREE GLOBAL MATRIX.*/

  return 0;
}/*arclmtest*/

void initialelem001(struct owire *elems,
                    struct memoryelem *melem,int nelem)
/*ASSEMBLAGE CONFINEMENT OF ELEMENT,CLEAR STRESS.*/
{
  int i,j;

  for(i=0;i<nelem;i++)
  {
    (melem+i)->code=(elems+i)->code;                        /*CODE.*/

    for(j=0;j<6;j++)                                 /*ICONF[2][6].*/
    {
      (melem+i)->bond[0][j]=(elems+i)->iconf[0][j];
      (melem+i)->bond[1][j]=(elems+i)->iconf[1][j];
    }
    for(j=0;j<6;j++)                           /*LONG STRESS[2][6].*/
    {
      (melem+i)->stress[0][j]=0.0;
      (melem+i)->stress[1][j]=0.0;
    }
  }

  return;
}/*initialelem001*/

double *elemstress001(struct owire *elem,
                      double *gvct,struct memoryelem *melem,struct arclmframe *af)
/*ELEMENT STRESS INCREMENTAL.*/
{
  int ii;
  double **drccos,**tmatrix,**estiff,*gdisp,*edisp,*estress;

  drccos=directioncosine(elem->node[0]->d[0],
                         elem->node[0]->d[1],
                         elem->node[0]->d[2],
                         elem->node[1]->d[0],
                         elem->node[1]->d[1],
                         elem->node[1]->d[2],
                         elem->cangle);                  /*[DRCCOS]*/
  tmatrix=transmatrix(/* *elem,*/drccos);                     /*[T]*/
  estiff=assememtx(*elem);                                   /*[ke]*/
  estiff=modifyhinge(*elem,estiff,af);

  gdisp=extractdisplacement(*elem,gvct);                     /*{dU}*/
  edisp=matrixvector(tmatrix,gdisp,12);              /*{du}=[T]{dU}*/
  estress=matrixvector(estiff,edisp,12);             /*{df}=[k]{du}*/

  updatestress001(melem,estress,elem); /*{f}+{df}*/

  free(gdisp);
  free(edisp);
  for(ii=0;ii<=2;ii++) free(*(drccos+ii));
  free(drccos);
  for(ii=0;ii<=11;ii++) free(*(tmatrix+ii));
  free(tmatrix);
  for(ii=0;ii<=11;ii++) free(*(estiff+ii));
  free(estiff);

  return estress;
}/*elemstress001*/

double *elemstress002(struct owire *elem,
                      double *gvct,struct memoryelem *melem,
                      long int *moff,struct arclmframe *af)
/*ELEMENT STRESS INCREMENTAL.*/
{
  int ii;
  double **drccos,**tmatrix,**estiff,*gdisp,*edisp,*estress;

  drccos=directioncosine(elem->node[0]->d[0],
                         elem->node[0]->d[1],
                         elem->node[0]->d[2],
                         elem->node[1]->d[0],
                         elem->node[1]->d[1],
                         elem->node[1]->d[2],
                         elem->cangle);                  /*[DRCCOS]*/
  tmatrix=transmatrix(/* *elem,*/drccos);                     /*[T]*/
  estiff=assememtx(*elem);                                   /*[ke]*/
  estiff=modifyhinge(*elem,estiff,af);

  gdisp=extractdisplacement2(*elem,gvct,moff);               /*{dU}*/
  edisp=matrixvector(tmatrix,gdisp,12);              /*{du}=[T]{dU}*/
  estress=matrixvector(estiff,edisp,12);             /*{df}=[k]{du}*/

  updatestress001(melem,estress,elem); /*{f}+{df}*/

  free(gdisp);
  free(edisp);
  for(ii=0;ii<=2;ii++) free(*(drccos+ii));
  free(drccos);
  for(ii=0;ii<=11;ii++) free(*(tmatrix+ii));
  free(tmatrix);
  for(ii=0;ii<=11;ii++) free(*(estiff+ii));
  free(estiff);

  return estress;
}/*elemstress002*/

double *elemstress003(FILE *fout,struct owire *elem,
                      double *gvct,struct memoryelem *melem,struct arclmframe *af)
/*ELEMENT STRESS INCREMENTAL FOR ARCLMTEST.*/
{
  int ii,jj;
  double **drccos,**tmatrix,**estiff,*gdisp,*edisp,*estress;

  drccos=directioncosine(elem->node[0]->d[0],
                         elem->node[0]->d[1],
                         elem->node[0]->d[2],
                         elem->node[1]->d[0],
                         elem->node[1]->d[1],
                         elem->node[1]->d[2],
                         elem->cangle);                  /*[DRCCOS]*/

  /*fprintf(fout,"ï˚å¸ó]å∑ [Çc] =");
  for(ii=0;ii<3;ii++)
  {
    if(ii>0) fprintf(fout,"               ");
    for(jj=0;jj<3;jj++)
    {
      fprintf(fout," %6.3f",*((*(drccos+ii))+jj));
    }
    fprintf(fout,"\n");
  }*/

  tmatrix=transmatrix(/* *elem,*/drccos);                     /*[T]*/

  /*fprintf(fout,"ç¿ïWïœä∑çsóÒ [Çs] =");
  for(ii=0;ii<12;ii++)
  {
    if(ii>0) fprintf(fout,"                   ");
    for(jj=0;jj<12;jj++)
    {
      fprintf(fout," %6.3f",*((*(tmatrix+ii))+jj));
    }
    fprintf(fout,"\n");
  }*/

  estiff=assememtx(*elem);                                   /*[ke]*/
  estiff=modifyhinge(*elem,estiff,af);

  fprintf(fout,"ïîçﬁçÑê´çsóÒ [Çãe] =");
  for(ii=0;ii<12;ii++)
  {
    if(ii>0) fprintf(fout,"                    ");
    for(jj=0;jj<12;jj++)
    {
      fprintf(fout," %5.0f",*((*(estiff+ii))+jj));
    }
    fprintf(fout,"\n");
  }

  gdisp=extractdisplacement(*elem,gvct);                     /*{dU}*/

  fprintf(fout,"ëSëÃç¿ïWÇ≈ÇÃïîçﬁïœà  {Çt} =");
  for(ii=0;ii<12;ii++)
  {
    if(ii>0) fprintf(fout,"                           ");
    fprintf(fout," %7.3f\n",*(gdisp+ii));
  }

  edisp=matrixvector(tmatrix,gdisp,12);              /*{du}=[T]{dU}*/

  fprintf(fout,"ïîçﬁç¿ïWÇ≈ÇÃïîçﬁïœà  {Çï} = [Çs]{Çt} =");
  for(ii=0;ii<12;ii++)
  {
    if(ii>0) fprintf(fout,"                                    ");
    fprintf(fout," %7.3f\n",*(edisp+ii));
  }

  estress=matrixvector(estiff,edisp,12);             /*{df}=[k]{du}*/

  fprintf(fout,"ïîçﬁâûóÕ {ÇÜ} = [Çãe]{Çï} =");
  for(ii=0;ii<12;ii++)
  {
    if(ii>0) fprintf(fout,"                           ");
    for(jj=0;jj<12;jj++)
    {
      fprintf(fout," %5.0f",*((*(estiff+ii))+jj));
    }
    fprintf(fout,"   %7.3f",*(edisp+ii));
    fprintf(fout,"   %7.3f",*(estress+ii));
    fprintf(fout,"\n");
  }

  updatestress001(melem,estress,elem); /*{f}+{df}*/

  free(gdisp);
  free(edisp);
  for(ii=0;ii<=2;ii++) free(*(drccos+ii));
  free(drccos);
  for(ii=0;ii<=11;ii++) free(*(tmatrix+ii));
  free(tmatrix);
  for(ii=0;ii<=11;ii++) free(*(estiff+ii));
  free(estiff);

  return estress;
}/*elemstress003*/

void updatestress001(struct memoryelem *melem,double *dstress,
                     struct owire *elem)
/*ELEMENT STRESS UPDATE.*/
{
  int i,j;

  for(i=0;i<2;i++)                                 /*UPDATE STRESS.*/
  {
    for(j=0;j<6;j++)
    {
      elem->stress[i][j]+=*(dstress+6*i+j);
      (melem+(elem->loff))->stress[i][j]=elem->stress[i][j];
    }
  }

  return;
}/*updatestress001*/

void outputdisp001(double *gvct,FILE *fout,int nnode,
                   struct onode *nodes)
/*OUTPUT NODE DISPLACEMENT AS FRAME3.*/
{
  char string[256];
  int i,j;
  double data[6];

  if(fout!=NULL)
  {
    fprintf(fout,"** DISPLACEMENT OF NODE\n\n");
    fprintf(fout,"  NO          U          V          W");
    fprintf(fout,"         KSI         ETA       OMEGA\n\n");
  }
  for(i=1;i<=nnode;i++)
  {
    for(j=0;j<=5;j++) data[j]=*(gvct+6*(i-1)+j);
    sprintf(string,
            "%4d %10.6f %10.6f %10.6f %11.7f %11.7f %11.7f",
            (nodes+i-1)->code,
            data[0],data[1],data[2],data[3],data[4],data[5]);
    if(fout!=NULL) fprintf(fout,"%s\n",string);
  }
}/*outputdisp001*/

void outputdisp002(double *gvct,FILE *fout,int nnode,
                   struct onode *nodes,long int *moff)
/*OUTPUT NODE DISPLACEMENT AS FRAME3.*/
{
  char string[256];
  int i,j,n;
  double data[6];

  if(fout!=NULL)
  {
    fprintf(fout,"** DISPLACEMENT OF NODE\n\n");
    fprintf(fout,"  NO          U          V          W");
    fprintf(fout,"         KSI         ETA       OMEGA\n\n");
  }
  for(i=0;i<nnode;i++)
  {
    if(moff==NULL) n=i;
    else           n=*(moff+i);

    for(j=0;j<6;j++) data[j]=*(gvct+6*n+j);
    sprintf(string,
            "%4d %10.6f %10.6f %10.6f %11.7f %11.7f %11.7f",
            (nodes+i)->code,
            data[0],data[1],data[2],data[3],data[4],data[5]);
    if(fout!=NULL) fprintf(fout,"%s\n",string);
  }
}/*outputdisp002*/

void outputstress001(struct owire elem,
                     double *estress,FILE *fout)
/*ELEMENT STRESS OUTPUT.*/
{
  char s[80],string[256];
  int i,n,nn[2];

  for(n=0;n<=1;n++)
  {
    nn[n]=elem.node[n]->code;
    if(n==0) sprintf(string,"%5d %4d %4d",
                     elem.code,elem.sect->code,nn[n]);
    if(n==1) sprintf(string,"           %4d",nn[n]);
    for(i=0;i<=5;i++)
    {
      /*sprintf(s," %9.3f",*(estress+6*n+i));*/       /*WITHOUT CMQ*/
      sprintf(s," %15.12f",elem.stress[n][i]);             /*WITH CMQ*/
      strcat(string,s);
    }
    if(fout!=NULL) fprintf(fout,"%s\n",string);
  }
  return;
}/*outputstress001*/

void outputreaction001(struct gcomponent *gmtx,
                       double *gvct,
                       struct onode *nodes,
                       struct oconf *confs,
                       double *dreact,FILE *fout,int nnode)
/*REACTIONS UPDATE,OUTPUT AS FRAME3.*/
{
  char string[256];
  char iconf;
  int offset;
  long int i,j,msize,nreact=0;
  double gstiff,dreaction;

  msize=6*nnode;
  for(i=1;i<=msize;i++)
  {
    iconf=(confs+i-1)->iconf;
    offset=(int)((i-1)/6);

    if(iconf==1)
    {
      dreaction=0.0;
      for(j=1;j<=msize;j++)
      {
        gread(gmtx,i,j,&gstiff);
        dreaction+=gstiff*(*(gvct+j-1));             /*{dR}=[K]{dU}*/
      }

      *(dreact+nreact)+=dreaction;                       /*{R}+{dR}*/

      sprintf(string,"%4ld %10ld %14.6f     1",
              (nodes+offset)->code,(i-1)%6+1,dreaction);
                                          /*1:ONLY FIXED AVAILABLE.*/
      if(fout!=NULL) fprintf(fout,"%s\n",string);

      nreact++;
    }

    currentpivot(i,msize);
  }
  return;
}/*outputreaction001*/

void outputreaction002(struct gcomponent *gmtx,
                       double *gvct,
                       struct onode *nodes,
                       struct oconf *confs,
                       double *dreact,FILE *fout,int nnode,
                       long int *moff)
/*REACTIONS UPDATE,OUTPUT AS FRAME3 USING EXCHANGED MATRIX.*/
/*GVCT,CONFS REEXCHANGED.*/
{
  char string[256];
  char iconf;
  long int i,ii,j,jj,m,n,msize,nreact=0;
  double gstiff,dreaction;

  msize=6*nnode;
  for(i=0;i<nnode;i++)
  {
    for(ii=0;ii<6;ii++)
    {
      if(moff==NULL) m=6*i+ii;
      else           m=6**(moff+i)+ii;

      /*iconf=(confs+6*i+ii)->iconf;*/
      iconf=(confs+m)->iconf;

      if(iconf==1)
      {
        dreaction=0.0;
        for(j=0;j<nnode;j++)
        {
          for(jj=0;jj<6;jj++)
          {
            if(moff==NULL) n=6*j+jj;
            else           n=6**(moff+j)+jj;
            gread(gmtx,m+1,n+1,&gstiff);
            /*dreaction+=gstiff*(*(gvct+6*j+jj));*/      /*{dR}=[K]{dU}*/
            dreaction+=gstiff*(*(gvct+n));           /*{dR}=[K]{dU}*/
          }
        }

        *(dreact+nreact)+=dreaction;                     /*{R}+{dR}*/

        sprintf(string,"%4ld %10ld %14.6f     1",
                (nodes+i)->code,ii+1,dreaction);
                                          /*1:ONLY FIXED AVAILABLE.*/
        if(fout!=NULL) fprintf(fout,"%s\n",string);

        nreact++;
      }
      currentpivot(6*i+ii+1,msize);
    }
  }
  return;
}/*outputreaction002*/

/* TODO: 2018-08-09 UNDERCONSTRUCTION FOR ROTATIONAL SPRING */
void modifycmq(struct memoryelem *melem,struct owire *elem)
/*MODIFY CMQ BY HINGE.*/
{
  int i,j;
  double dl,dx,dy,dz;

  dx=(elem->node[1]->d[0])-(elem->node[0]->d[0]);
  dy=(elem->node[1]->d[1])-(elem->node[0]->d[1]);
  dz=(elem->node[1]->d[2])-(elem->node[0]->d[2]);
  dl=sqrt(dx*dx+dy*dy+dz*dz);

  if((elem->iconf[0][4]==1)&&(elem->iconf[1][4]==1))
  {
    elem->stress[0][4]=0.0;
    elem->stress[1][4]=0.0;
  }
  if((elem->iconf[0][5]==1)&&(elem->iconf[1][5]==1))
  {
    elem->stress[0][5]=0.0;
    elem->stress[1][5]=0.0;
  }
  if((elem->iconf[0][4]==1)&&(elem->iconf[1][4]==0))
  {
    elem->stress[1][4]-=elem->stress[0][4]/2.0;
    elem->stress[0][2]+=elem->stress[0][4]*1.5/dl;
    elem->stress[1][2]-=elem->stress[0][4]*1.5/dl;
    elem->stress[0][4]=0.0;
  }
  if((elem->iconf[0][4]==0)&&(elem->iconf[1][4]==1))
  {
    elem->stress[0][4]-=elem->stress[1][4]/2.0;
    elem->stress[0][2]+=elem->stress[1][4]*1.5/dl;
    elem->stress[1][2]-=elem->stress[1][4]*1.5/dl;
    elem->stress[1][4]=0.0;
  }
  if((elem->iconf[0][5]==1)&&(elem->iconf[1][5]==0))
  {
    elem->stress[1][5]-=elem->stress[0][5]/2.0;
    elem->stress[0][1]-=elem->stress[0][5]*1.5/dl;
    elem->stress[1][1]+=elem->stress[0][5]*1.5/dl;
    elem->stress[0][5]=0.0;
  }
  if((elem->iconf[0][5]==1)&&(elem->iconf[1][5]==0))
  {
    elem->stress[0][5]-=elem->stress[1][5]/2.0;
    elem->stress[0][1]-=elem->stress[1][5]*1.5/dl;
    elem->stress[1][1]+=elem->stress[1][5]*1.5/dl;
    elem->stress[1][5]=0.0;
  }

  for(i=0;i<2;i++)                                 /*UPDATE STRESS.*/
  {
    for(j=0;j<6;j++)
    {
      (melem+(elem->loff))->stress[i][j]=elem->stress[i][j];
    }
  }

  return;
}/*modifycmq*/

void assemcmq(struct owire elem,double **tmatrix,
              struct oconf *confs,double *gvct)
/*ASSEMBLAGE CMQ INTO GLOBAL VECTOR.*/
{
  long int i,loffset;
  int n,node[2];
  double *load,*c,**t;

  c=(double *)malloc(12*sizeof(double));
  for(n=1;n<=2;n++)
  {
    for(i=1;i<=6;i++)
    {
      *(c+6*(n-1)+i-1)=elem.stress[n-1][i-1];
    }
  }
  t=matrixtranspose(tmatrix,12);
  load=matrixvector(t,c,12);
  node[0]=elem.node[0]->loff;
  node[1]=elem.node[1]->loff;
  for(n=0;n<=1;n++)
  {
    for(i=0;i<=5;i++)
    {
      loffset=6*node[n]+i;
      if((confs+loffset)->iconf==0)
      {
        *(gvct+loffset)-=*(load+6*n+i);
      }
    }
  }

  for(i=0;i<=11;i++) free(*(t+i));
  free(t);
  free(load);
  free(c);

  return;
}/*assemcmq*/

void assemconf001(struct oconf *confs,double *gvct,
                  double dsafety,int nnode)
/*ASSEMBLAGE CONFINEMENT VALUE INTO GLOBAL VECTOR.*/
{
  int i;
  double conf;

  for(i=1;i<=6*nnode;i++)     /*LOADS OR DISPS GIVEN INCREMENTALLY.*/
  {
    conf=(confs+i-1)->value;

    *(gvct+i-1)+=dsafety*conf;                 /*"+=":ADD WITH CMQ.*/
  }

  return;
}/*assemconf001*/

void frameoutputtomemory(FILE *ftext,struct arclmframe *af)
/*TRANSLATE FRAME OUTPUTFILE TEXT INTO MEMORY.*/
{
  char **data,str[256]="\0";
  int i,j,k,n;
  long int ecode,hcode;
  double ddata;

  fseek(ftext,0L,SEEK_SET);

  /*"OUTPUTFILE SPECIFICATION"*/
  /*APPELATION*/
  /*ELEMCODE SECT HEAD N Q1 Q2 MT M1 M2*/
  /*              TAIL N Q1 Q2 MT M1 M2*/
  /*NODECODE U V W TX TY TZ*/
  /*NODECODE DIRECTION REACTIONVALUE CONFINEMENTTYPE*/

  af->ddisp=(double *)malloc(6*(af->nnode)*sizeof(double));
                                       /*DISPLACEMENT:6 DIRECTIONS.*/
  af->melem=(struct memoryelem *)
            malloc((af->nelem)*sizeof(struct memoryelem));
                                        /*CODE,12 BOUNDS,12 STRESS.*/
  af->dreact=(double *)malloc((af->nreact)*sizeof(double));
                                                        /*REACTION.*/

  while(strncmp(str,"** FORCES",9)) fgets(str,256,ftext);
  for(i=1;i<=3;i++) fgets(str,256,ftext);

  for(i=1;i<=(af->nelem);i++) /*STRESSES.*/
  {
    data=fgetsbrk(ftext,&n);
    if(n!=9) return;

    ecode=strtol(*(data+0),NULL,10);
    if(ecode!=(af->elems+i-1)->code) return;

    /*scode=strtol(*(data+1),NULL,10);*/ /*SECTION.*/
    /*hcode=strtol(*(data+2),NULL,10);*/ /*HEAD NODE.*/

    k=3;
    for(j=0;j<6;j++)
    {
      (af->elems+i-1)->stress[0][j]=strtod(*(data+k),NULL);
      k++;
    }

    for(;n>0;n--) free(*(data+n-1));
    free(data);

    data=fgetsbrk(ftext,&n);
    if(n!=7) return;

    /*tcode=strtol(*(data+0),NULL,10);*/ /*TAIL NODE.*/

    k=1;
    for(j=0;j<6;j++)
    {
      (af->elems+i-1)->stress[1][j]=strtod(*(data+k),NULL);
      k++;
    }

    for(;n>0;n--) free(*(data+n-1));
    free(data);

    (af->melem+i-1)->code=(af->elems+i-1)->code;
    for(j=0;j<6;j++)
    {
      (af->melem+i-1)->bond[0][j]=(af->elems+i-1)->iconf[0][j];
      (af->melem+i-1)->bond[1][j]=(af->elems+i-1)->iconf[1][j];
      (af->melem+i-1)->stress[0][j]=(af->elems+i-1)->stress[0][j];
      (af->melem+i-1)->stress[1][j]=(af->elems+i-1)->stress[1][j];
    }
  }

  while(strncmp(str,"** DISPLACEMENT",15)) fgets(str,256,ftext);
  for(i=1;i<=3;i++) fgets(str,256,ftext);

  for(i=1;i<=(af->nnode);i++) /*DISPLACEMENTS.*/
  {
    data=fgetsbrk(ftext,&n);
    if(n!=7) return;

    hcode=strtol(*(data+0),NULL,10);
    if(hcode!=(af->nodes+i-1)->code) return;

    for(j=1;j<=3;j++)
    {
      ddata=strtod(*(data+j),NULL);
      ddata+=(af->ninit+i-1)->d[j-1];
      (af->nodes+i-1)->d[j-1]=ddata;
      *(af->ddisp+6*(i-1)+(j-1))=ddata;
    }
    for(j=4;j<=6;j++)
    {
      ddata=strtod(*(data+j),NULL);
      *(af->ddisp+6*(i-1)+(j-1))=ddata;
    }

    for(;n>0;n--) free(*(data+n-1));
    free(data);
  }

  while(strncmp(str,"** REACTION",11)) fgets(str,256,ftext);
  for(i=1;i<=3;i++) fgets(str,256,ftext);

  for(i=0;i<(af->nreact);i++) /*REACTIONS.*/
  {
    data=fgetsbrk(ftext,&n);
    if(n!=4) return;

    ddata=strtod(*(data+2),NULL);
    *(af->dreact+i)=ddata;

    for(;n>0;n--) free(*(data+n-1));
    free(data);
  }

  return;
}/*frameoutputtomemory*/

void openarclmlastfile(FILE *ftext,struct arclmframe *af)
/*TRANSLATE ARCLM101 OUTPUTFILE TEXT INTO MEMORY.*/
{
  char **data,str[256]="\0";
  int i,j,k,n;
  long int ecode,hcode;
  double ddata,func;

  fseek(ftext,0L,SEEK_SET);

  /*"OUTPUTFILE SPECIFICATION"*/
  /*APPELATION*/
  /*ELEMCODE SECT HEAD N Q1 Q2 MT M1 M2*/
  /*              TAIL N Q1 Q2 MT M1 M2*/
  /*NODECODE U V W TX TY TZ*/
  /*NODECODE DIRECTION REACTIONVALUE CONFINEMENTTYPE*/

  af->ddisp=(double *)malloc(6*(af->nnode)*sizeof(double));
                                       /*DISPLACEMENT:6 DIRECTIONS.*/
  af->melem=(struct memoryelem *)
            malloc((af->nelem)*sizeof(struct memoryelem));
                                        /*CODE,12 BOUNDS,12 STRESS.*/
  af->dreact=(double *)malloc((af->nreact)*sizeof(double));
                                                        /*REACTION.*/

  while(strncmp(str,"\"STRESS\"",8)) fgets(str,256,ftext);

  for(i=1;i<=(af->nelem);i++) /*STRESSES.*/
  {
    /*HEAD*/
    data=fgetsbrk(ftext,&n);
    while(strcmp(*(data+0),"ELEM"))
    {
      for(;n>0;n--) free(*(data+n-1));
      free(data);
      data=fgetsbrk(ftext,&n);
    }
    if(n!=20) return;

    ecode=strtol(*(data+1),NULL,10);
    if(ecode!=(af->elems+i-1)->code) return;

    /*scode=strtol(*(data+3),NULL,10);*/ /*SECTION.*/
    /*hcode=strtol(*(data+5),NULL,10);*/ /*HEAD NODE.*/

    k=7;
    for(j=0;j<6;j++)
    {
      (af->elems+i-1)->stress[0][j]=strtod(*(data+k),NULL);
      k+=2;
    }

    func=strtod(*(data+19),NULL);
    if(func>pow(RADIUS,EXPONENT))
    {
      for(j=0;j<6;j++)
      {
        (af->elems+i-1)->iconf[0][j]=-1; /*-1:PLASTIC*/
      }
    }

    for(;n>0;n--) free(*(data+n-1));
    free(data);

    /*TAIL*/
    data=fgetsbrk(ftext,&n);
    while(strcmp(*(data+0),"ELEM"))
    {
      for(;n>0;n--) free(*(data+n-1));
      free(data);
      data=fgetsbrk(ftext,&n);
    }
    if(n!=20) return;

    ecode=strtol(*(data+1),NULL,10);
    if(ecode!=(af->elems+i-1)->code) return;

    /*tcode=strtol(*(data+5),NULL,10);*/ /*TAIL NODE.*/

    k=7;
    for(j=0;j<6;j++)
    {
      (af->elems+i-1)->stress[1][j]=strtod(*(data+k),NULL);
      k+=2;
    }

    func=strtod(*(data+19),NULL);
    if(func>pow(RADIUS,EXPONENT))
    {
      for(j=0;j<6;j++)
      {
        (af->elems+i-1)->iconf[1][j]=-1; /*-1:PLASTIC*/
      }
    }

    for(;n>0;n--) free(*(data+n-1));
    free(data);

    (af->melem+i-1)->code=(af->elems+i-1)->code;
    for(j=0;j<6;j++)
    {
      (af->melem+i-1)->bond[0][j]=(af->elems+i-1)->iconf[0][j];
      (af->melem+i-1)->bond[1][j]=(af->elems+i-1)->iconf[1][j];
      (af->melem+i-1)->stress[0][j]=(af->elems+i-1)->stress[0][j];
      (af->melem+i-1)->stress[1][j]=(af->elems+i-1)->stress[1][j];
    }
  }

  while(strncmp(str,"\"REACTION\"",10)) fgets(str,256,ftext);

  for(i=0;i<(af->nreact);i++) /*REACTIONS.*/
  {
    data=fgetsbrk(ftext,&n);
    if(n!=4) return;

    ddata=strtod(*(data+4),NULL);
    *(af->dreact+i)=ddata;

    for(;n>0;n--) free(*(data+n-1));
    free(data);
  }

  while(strncmp(str,"\"CURRENT FORM\"",14)) fgets(str,256,ftext);

  for(i=1;i<=(af->nnode);i++) /*DISPLACEMENTS.*/
  {
    data=fgetsbrk(ftext,&n);
    if(n!=7) return;

    hcode=strtol(*(data+1),NULL,10);
    if(hcode!=(af->nodes+i-1)->code) return;

    for(j=2;j<=4;j++)
    {
      ddata=strtod(*(data+j),NULL);
      (af->nodes+i-1)->d[j-1]=ddata;
      *(af->ddisp+6*(i-1)+(j-1))=ddata;
    }
    for(j=5;j<=7;j++)
    {
      ddata=strtod(*(data+j),NULL);
      *(af->ddisp+6*(i-1)+(j-1))=ddata;
    }

    for(;n>0;n--) free(*(data+n-1));
    free(data);
  }

  return;
}/*openarclmlastfile*/

