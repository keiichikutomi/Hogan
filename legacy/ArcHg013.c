/*ARCLM001.C FOR WIN32 SINCE 1995.11.24.JUNSATO.*/
/*LAST CHANGE:1997.11.30.*/
 
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
void initialelem001(struct owire *elems,FILE *felem,int nelem);
double *elemstress001(struct owire *elem,
                      double *gvct,FILE *felem);
void updatestress001(FILE *felem,double *dstress,
                     struct owire *elem);
void outputdisp001(double *gvct,FILE *fout,int nnode,
                   struct onode *nodes);
void outputstress001(struct owire elem,
                     double *estress,FILE *fout);
void outputreaction001(struct gcomponent *gmtx,
                       double *gvct,
                       struct onode *nodes,
                       struct oconf *confs,
                       FILE *freact,FILE *fout,int nnode);
void modifycmq(FILE *felem,struct owire *elem);
void assemcmq(struct owire elem,double **tmatrix,
              struct oconf *confs,double *gvct);
void assemconf001(struct oconf *confs,double *gvct,
                  double dsafety,int nnode);

void frameoutputtomemory(FILE *ftext,struct arclmframe *af);

int arclm001(struct arclmframe *af,int idinput,int idoutput)
{
  DWORD memory0,memory1,memory2;

  FILE *fin,*fout;                                   /*FILE 8 BYTES*/
  FILE *felem,*fdisp,*freact;
  char dir[]=DIRECTORY;                            /*DATA DIRECTORY*/
  char /*s[20],*/string[400];
  int i,ii;
  int nnode,nelem,nsect,nreact;
  long int fsize,msize;
  struct gcomponent ginit={0,0,0.0,NULL};
  struct gcomponent *gmtx;                          /*GLOBAL MATRIX*/
  double *gvct;                                     /*GLOBAL VECTOR*/
  double **drccos,**tmatrix,**estiff,*estress;
  double determinant,sign;
  clock_t t0;

  struct osect *sects;
  struct onode *nodes;
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
  fout=fgetstofopen("\0","w",idoutput);               /*OUTPUT FILE*/

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
  for(i=0;i<=msize-1;i++)
  {
    (gmtx+i)->down=NULL;            /*GLOBAL MATRIX INITIALIZATION.*/
    *(gvct+i)=0.0;                  /*GLOBAL VECTOR INITIALIZATION.*/
  }

  fdisp=fopen("canbin.dsp","wb+");     /*DISPLACEMENT:6 DIRECTIONS.*/
  af->fdisp=fdisp;
  fsize=sizeof(double);
  setvbuf(fdisp,NULL,_IOFBF,fsize);

  felem=fopen("canbin.elm","wb+"); /*IELEM,12 BOUNDARIES,12 STRESS.*/
  af->felem=felem;
  fsize=sizeof(long int)+12*sizeof(signed char)+12*sizeof(double);
  setvbuf(felem,NULL,_IOFBF,fsize);

  free(af->sects);
  free(af->nodes);
  free(af->elems);
  free(af->confs);

  sects=(struct osect *)malloc(nsect*sizeof(struct osect));
  if(sects==NULL) return 0;
  nodes=(struct onode *)malloc(nnode*sizeof(struct onode));
  if(nodes==NULL) return 0;
  elems=(struct owire *)malloc(nelem*sizeof(struct owire));
  if(elems==NULL) return 0;
  confs=(struct oconf *)malloc(msize*sizeof(struct oconf));
  if(confs==NULL) return 0;

  af->sects=sects;
  af->nodes=nodes;
  af->elems=elems;
  af->confs=confs;

  inputtexttomemory(fin,af);        /*READY TO READ LONG REACTIONS.*/
  nnode=af->nnode;
  nelem=af->nelem;
  nsect=af->nsect;
  nreact=af->nreact;

  initialform(nodes,fdisp,nnode);           /*ASSEMBLAGE FORMATION.*/

  initialelem001(elems,felem,nelem);         /*ASSEMBLAGE ELEMENTS.*/

  freact=fopen("canbin.rct","wb+");                /*REACTION FILE.*/
  af->freact=freact;
  fsize=sizeof(double);
  setvbuf(freact,NULL,_IOFBF,fsize);
  initialreact(fin,freact,nreact);     /*ASSEMBLAGE LONG REACTIONS.*/

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

  laptime("ASSEMBLING GLOBAL MATRIX.",t0);

  for(i=1;i<=nelem;i++)               /*ASSEMBLAGE GLOBAL MATRIX.*/
  {
    elem=*(elems+i-1);                       /*READ ELEMENT DATA.*/

    inputnode(fdisp,elem.node[0]);                         /*HEAD*/
    inputnode(fdisp,elem.node[1]);                         /*TAIL*/

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
    estiff=modifyhinge(elem,estiff);             /*MODIFY MATRIX.*/
    estiff=transformation(estiff,tmatrix);       /*[K]=[Tt][k][T]*/

    assemgstiffness(gmtx,estiff,&elem);             /*ASSEMBLAGE.*/

    modifycmq(felem,&elem);
    assemcmq(elem,tmatrix,confs,gvct); /*ASSEMBLAGE CMQ AS LOADS.*/

    for(ii=0;ii<=2;ii++) free(*(drccos+ii));
    free(drccos);
    for(ii=0;ii<=11;ii++) free(*(tmatrix+ii));
    free(tmatrix);
    for(ii=0;ii<=11;ii++) free(*(estiff+ii));
    free(estiff);
  }
  sprintf(string,"GLOBAL MATRIX %ld COMPS ASSEMBLED.",comps);
  laptime(string,t0);

  overlayhdc(*(wdraw.childs+1),SRCPAINT);       /*UPDATE DISPLAY.*/

  /*errormessage("ASSEMBLAGE GLOBAL VECTOR.");*/
  assemconf001(confs,gvct,1.0,nnode);            /*GLOBAL VECTOR.*/
  modifygivend(gmtx,gvct,confs,nnode);   /*0:LOAD 1:DISPLACEMENT.*/

  laptime("CROUT LU DECOMPOSITION.",t0);
  croutludecomposition(gmtx,
                       gvct,confs,
                       6*nnode,
                       &determinant,&sign);        /*[K]{dU}={dF}*/

  sprintf(string,"LOG{DETERMINANT}=%.5E SIGN=%.1f COMPS=%ld",
          determinant,sign,comps);
  errormessage(string);
  /*if(fout!=NULL) fprintf(fout,"%s\n",string);*/

  if(sign<=0.0)
  {
    errormessage(" ");
    errormessage("INSTABLE TERMINATION.");
    if(fout!=NULL) fprintf(fout,"INSTABLE TERMINATION.\n");

    laptime("\0",t0);

    fclose(fin);

    gfree(gmtx,nnode); /*FREE GLOBAL MATRIX.*/
    free(gvct);
    /*free(confs);*/

    memory2=availablephysicalmemory("REMAIN:");
    sprintf(string,"CONSUMPTION:%ld[BYTES]",(memory0-memory2));
    errormessage(string);

    return 1;
  }

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
    inputelem(elems,felem,i-1,&elem);

    inputnode(fdisp,elem.node[0]);
    inputnode(fdisp,elem.node[1]);

    elem.sect=(elems+i-1)->sect;             /*READ SECTION DATA.*/

    estress=elemstress001(&elem,gvct,felem);

    outputstress001(elem,estress,fout);
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
  outputreaction001(gmtx,gvct,nodes,confs,freact,fout,nnode);

  fclose(fin);
  fclose(fout);
  errormessage("FILES CLOSED.");

  updateform(fdisp,gvct,nnode);               /*FORMATION UPDATE.*/

  laptime("\0",t0);

  memory2=availablephysicalmemory(NULL);
  sprintf(string,"CONSUMPTION:%ld[BYTES]",(memory1-memory2));
  errormessage(string);

  if((wdraw.childs+1)->hdcC!=NULL &&
     felem!=NULL && fdisp!=NULL)               /*DRAW LAST FRAME.*/
  {
    for(i=1;i<=nelem;i++)
    {
      elem=*(elems+i-1);                     /*READ ELEMENT DATA.*/

      inputnode(fdisp,elem.node[0]);
      inputnode(fdisp,elem.node[1]);

      drawglobalwire((wdraw.childs+1)->hdcC,
                     (wdraw.childs+1)->vparam,
                     *af,elem,255,255,255,
                              255,255,255,0,ONSCREEN);
    }
    overlayhdc(*(wdraw.childs+1),SRCPAINT);     /*UPDATE DISPLAY.*/
  }

  gfree(gmtx,nnode); /*FREE GLOBAL MATRIX.*/
  /*free(gvct);*/
  /*free(confs);*/

  af->nlaps=1;
  af->eigenvec=(double **)malloc(1*sizeof(double *));
  *((af->eigenvec)+0)=gvct;

  errormessage(" ");
  errormessage("COMPLETED.");

  memory2=availablephysicalmemory("REMAIN:");
  sprintf(string,"CONSUMPTION:%ld[BYTES]",(memory0-memory2));
  errormessage(string);

  return 0;
}/*arclm001*/

void initialelem001(struct owire *elems,FILE *felem,int nelem)
/*ASSEMBLAGE CONFINEMENT OF ELEMENT,CLEAR STRESS.*/
{
  int i,j;
  double zero=0.0;

  for(i=0;i<nelem;i++)
  {
    fwrite(&((elems+i)->code),sizeof(long int),1,felem);    /*CODE.*/

    for(j=0;j<=1;j++)                                /*ICONF[2][6].*/
    {
      fwrite(&((elems+i)->iconf[j]),sizeof(signed char),6,felem);
    }
    for(j=0;j<=11;j++)                              /*STRESS[2][6].*/
    {
      fwrite(&zero,sizeof(double),1,felem);
    }
  }
  fflush(felem);

  return;
}/*initialelem001*/

double *elemstress001(struct owire *elem,
                      double *gvct,FILE *felem)
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
  estiff=modifyhinge(*elem,estiff);

  gdisp=extractdisplacement(*elem,gvct);                     /*{dU}*/
  edisp=matrixvector(tmatrix,gdisp,12);              /*{du}=[T]{dU}*/
  estress=matrixvector(estiff,edisp,12);             /*{df}=[k]{du}*/

  updatestress001(felem,estress,elem); /*{f}+{df}*/

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

void updatestress001(FILE *felem,double *dstress,
                     struct owire *elem)
/*ELEMENT STRESS UPDATE.*/
{
  int i,j;
  long int loffset;

  loffset=(elem->loff)
         *(sizeof(long int)+12*sizeof(char)+12*sizeof(double))
         +sizeof(long int)+12*sizeof(char);
  fseek(felem,loffset,SEEK_SET);
  for(i=1;i<=2;i++)                                /*UPDATE STRESS.*/
  {
    for(j=1;j<=6;j++)
    {
      elem->stress[i-1][j-1]+=*(dstress+6*(i-1)+j-1);
    }
    fwrite(elem->stress[i-1],sizeof(double),6,felem);
  }
  fflush(felem);

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
      sprintf(s," %9.3f",elem.stress[n][i]);             /*WITH CMQ*/
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
                       FILE *freact,FILE *fout,int nnode)
/*REACTIONS UPDATE,OUTPUT AS FRAME3.*/
{
  char string[256];
  char iconf;
  int offset;
  long int i,j,nreact=0;
  double gstiff,reaction,dreaction;

  for(i=1;i<=6*nnode;i++)
  {
    iconf=(confs+i-1)->iconf;
    offset=(i-1)/6;

    if(iconf==1)
    {
      nreact++;
      dreaction=0.0;
      for(j=1;j<=6*nnode;j++)
      {
        gread(gmtx,i,j,&gstiff);
        dreaction+=gstiff*(*(gvct+j-1));             /*{dR}=[K]{dU}*/
      }
      vread(freact,nreact,&reaction);            /*UPDATE REACTION.*/
      reaction+=dreaction;                               /*{R}+{dR}*/
      vwrite(freact,nreact,&reaction);

      sprintf(string,"%4d %10ld %14.6f     1",
              (nodes+offset)->code,(i-1)%6+1,dreaction);
                                          /*1:ONLY FIXED AVAILABLE.*/
      if(fout!=NULL) fprintf(fout,"%s\n",string);
    }
  }
  return;
}/*outputreaction001*/

void modifycmq(FILE *felem,struct owire *elem)
/*MODIFY CMQ BY HINGE.*/
{
  int i;
  long int loffset;
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

  loffset=(elem->loff)
         *(sizeof(long int)+12*sizeof(signed char)+12*sizeof(double))
         +sizeof(long int)+12*sizeof(signed char);
  fseek(felem,loffset,SEEK_SET);
  for(i=0;i<=1;i++)                                /*UPDATE STRESS.*/
  {
    fwrite(elem->stress[i],sizeof(double),6,felem);
  }
  fflush(felem);

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

  af->fdisp=fopen("canbin.dsp","wb+"); /*DISPLACEMENT:6 DIRECTIONS.*/
  af->felem=fopen("canbin.elm","wb+");  /*CODE,12 BOUNDS,12 STRESS.*/
  af->freact=fopen("canbin.rct","wb+");            /*REACTION FILE.*/

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

    fwrite(&((af->elems+i-1)->code),sizeof(long int),1,af->felem);
    for(j=0;j<=1;j++) fwrite((af->elems+i-1)->iconf[j],
                             sizeof(signed char),6,af->felem);
    for(j=0;j<=1;j++) fwrite((af->elems+i-1)->stress[j],
                             sizeof(double),6,af->felem);
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
      fwrite(&ddata,sizeof(double),1,af->fdisp);
    }
    for(j=4;j<=6;j++)
    {
      ddata=strtod(*(data+j),NULL);
      fwrite(&ddata,sizeof(double),1,af->fdisp);
    }

    for(;n>0;n--) free(*(data+n-1));
    free(data);
  }

  while(strncmp(str,"** REACTION",11)) fgets(str,256,ftext);
  for(i=1;i<=3;i++) fgets(str,256,ftext);

  for(i=1;i<=(af->nreact);i++) /*REACTIONS.*/
  {
    data=fgetsbrk(ftext,&n);
    if(n!=4) return;

    ddata=strtod(*(data+2),NULL);
    fwrite(&ddata,sizeof(double),1,af->freact);

    for(;n>0;n--) free(*(data+n-1));
    free(data);
  }

  return;
}/*frameoutputtomemory*/

