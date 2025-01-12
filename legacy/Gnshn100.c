/*=========================================================*/
/*   PROGRAM GUNBASHIRA SHINDOU 002                        */
/*     6 D.O.F. 3D FRAME WITH POINT MASS                   */
/*     GRAVITY LOADED.                                     */
/*     HISENKEI OUTOU KAISEKI                              */
/*   CODED BY JUN SATO                                     */
/*   DATE:FORTRAN SINSE 1993. 9.17                         */
/*        C       SINSE 1998.10.10                         */
/*                                                         */
/*   BASSUI FROM:                                          */
/*   'ZAKU001'  ,BASED ON 'Y6FRAMEB' BY YOSHINOBU FUJITANI */
/*   'D24PE3..9',BY TSUNETA OCHI                           */
/*=========================================================*/

/*#define MAXLAP 65*/
#define MAXLAP 200

struct accdata{long int ndata;
               double dt;
               double *a;}; /*ACCELERATION DATA.*/

int gnshn101(struct arclmframe *af);
double *inputacceleration(FILE *fin,
                          long int *ndata,double *dt);
void assemmass(struct gcomponent *mmtx,
               struct arclmframe *af);
double accelerationincrement(struct accdata acc,
                             double ddt,long int lap);
void assemaccel(double *gvct,
                double dacc[],struct arclmframe *af);
double newmarkbeta(struct gcomponent *gmtx,
                   struct gcomponent *mmtx,
                   double *gacc,
                   struct oconf *confs,
                   int nnode,
                   double *u, double *ud, double *udd,
                   double *du,double *dud,double *dudd,
                   double ddt,double beta);

int gnshn101(struct arclmframe *af)
/*DYNAMIC RESPONSE ANALYSIS.*/
/*GRAVITY LOADED ANALYSIS IS ALREADY DONE BY ARCLM.*/
{
  DWORD memory0,memory1,memory2;

  FILE *fin,*fout;                                   /*FILE 8 BYTES*/
  FILE *felem,*fdisp,*freact;
  /*char dir[]=DIRECTORY;*/                        /*DATA DIRECTORY*/
  char s[80],string[400];
  int i,ii,jj;
  int nnode,nelem,nsect,nreact;
  long int nlap,laps;
  long int loffset,msize;
  long int time;
  struct gcomponent ginit={0,0,0.0,NULL};
  struct gcomponent *gmtx,*g,*p;                    /*GLOBAL MATRIX*/
  /*double *gvct;*/                                 /*GLOBAL VECTOR*/
  double **drccos,**tmatrix,**estiff,*estress,**e,**t;
  double *gdisp,*edisp;
  double determinant,data,func[2];
  clock_t t0,t1,t2;
  struct onode *nodes;
  struct owire elem,*elems;
  struct oconf *confs;

  FILE *fxacc,*fyacc,*fzacc;
  long int ndata[3];
  double ddt,dacc[3],afact;
  double *gacc;
  double *u,*ud,*udd,*du,*dud,*dudd;                      /*VECTORS*/
  struct accdata acc[3];
  struct gcomponent *mmtx;                            /*MASS MATRIX*/

  memory0=availablephysicalmemory("INITIAL:");   /*MEMORY AVAILABLE*/

  /*INPUT FILE*/
  fin=NULL;

  /*OUTPUT FILE*/
  /*fout=fgetstofopen("\0","w",ID_OUTPUTFILE);*/
  fout=NULL;

  t0=clock();                                        /*CLOCK BEGIN.*/

  /*ACCCELERATION DATA INPUT.*/
  fxacc=fopen("elcnew.dat","r");
  fyacc=fopen("elcnns.dat","r");
  fzacc=NULL;

  acc[0].a=inputacceleration(fxacc,&(acc[0].ndata),&(acc[0].dt));
  acc[1].a=inputacceleration(fyacc,&(acc[1].ndata),&(acc[1].dt));
  acc[2].a=inputacceleration(fzacc,&(acc[2].ndata),&(acc[2].dt));
  /*ddt=0.005;*/
  ddt=0.005;
  /*afact=0.05;*/ /*ACCELERATION FACTOR*/
  afact=0.05;
  for(i=0;i<3;i++)
  {
    if(acc[i].a==NULL) ndata[i]=0;
    else
    {
      ndata[i]=(long int)((acc[i].ndata-1)*(acc[i].dt)/ddt)+1;
    }
  }
  if(ndata[0]>0 && ndata[0]<ndata[1]) laps=ndata[0];
  else                                laps=ndata[1];
  if(ndata[2]>0 && ndata[2]<laps)     laps=ndata[2];

  if(laps>MAXLAP) laps=MAXLAP;

  nnode=af->nnode;
  nelem=af->nelem;
  nsect=af->nsect;
  nreact=af->nreact;

  nodes=af->nodes;
  elems=af->elems;
  confs=af->confs;

  sprintf(string,"NODES=%d ELEMS=%d SECTS=%d",nnode,nelem,nsect);
  errormessage(string);
  if(fout!=NULL) fprintf(fout,"%s\n",string);

  msize=6*nnode;                           /*SIZE OF GLOBAL MATRIX.*/

  gmtx=(struct gcomponent *)          /*DIAGONALS OF GLOBAL MATRIX.*/
        malloc(msize*sizeof(struct gcomponent));
  gacc=mallocdoublevector(msize);                    /*ACCEL VECTOR*/
  if(gmtx==NULL || gacc==NULL) return 0;
  for(i=0;i<msize;i++)
  {
    (gmtx+i)->down=NULL;            /*GLOBAL MATRIX INITIALIZATION.*/
  }

  /*DISPLACEMENT VECTORS.*/
  u   =mallocdoublevector(msize);
  ud  =mallocdoublevector(msize);
  udd =mallocdoublevector(msize);
  du  =mallocdoublevector(msize);
  dud =mallocdoublevector(msize);
  dudd=mallocdoublevector(msize);
  if(u ==NULL || ud ==NULL || udd ==NULL ||
     du==NULL || dud==NULL || dudd==NULL) return 0;
  for(i=0;i<msize;i++)
  {
    *(u   +i)=0.0;
    *(ud  +i)=0.0;
    *(udd +i)=0.0;
    *(du  +i)=0.0;
    *(dud +i)=0.0;
    *(dudd+i)=0.0;
  }

  mmtx=(struct gcomponent *)            /*DIAGONALS OF MASS MATRIX.*/
        malloc(msize*sizeof(struct gcomponent));
  if(mmtx==NULL) return 0;
  for(i=0;i<=msize-1;i++) (mmtx+i)->down=NULL;

  fdisp=fopen("canbin.dsp","wb+");     /*DISPLACEMENT:6 DIRECTIONS.*/
  af->fdisp=fdisp;

  felem=fopen("canbin.elm","wb+");  /*CODE,12 BOUNDARIES,12 STRESS.*/
  af->felem=felem;

  initialform(nodes,fdisp,nnode);           /*ASSEMBLAGE FORMATION.*/

  initialelem(elems,felem,nelem);            /*ASSEMBLAGE ELEMENTS.*/

  freact=fopen("canbin.rct","wb+");                /*REACTION FILE.*/
  af->freact=freact;
  initialreact(fin,freact,nreact);     /*ASSEMBLAGE LONG REACTIONS.*/

  GetAsyncKeyState(VK_LBUTTON);                   /*CLEAR KEY LEFT.*/
  GetAsyncKeyState(VK_RBUTTON);                  /*CLEAR KEY RIGHT.*/

  drawglobalaxis((wdraw.childs+1)->hdcC,(wdraw.childs+1)->vparam,
                 0,0,255);                      /*DRAW GLOBAL AXIS.*/

  if(wsurf.hwnd!=NULL)
  {
    drawyieldsurface((wsurf.childs+1)->hdcC,
                     (wsurf.childs+1)->vparam,
                     SURFACEX,SURFACEY,SURFACEZ,NULL);
    overlayhdc(*(wsurf.childs+1),SRCPAINT);       /*UPDATE DISPLAY.*/
  }
  af->fsurface=fopen("canbin.sfc","wb+");      /*STRESS ON SURFACE.*/

  /*ASSEMBLAGE MASS MATRIX.*/
  assemmass(mmtx,af);

  /*MEMORIES.*/
  drccos =mallocdoublematrix(3);
  estiff =mallocdoublematrix(12);
  tmatrix=mallocdoublematrix(12);
  e      =mallocdoublematrix(12);
  t      =mallocdoublematrix(12);
  estress=mallocdoublevector(12);
  gdisp  =mallocdoublevector(12);
  edisp  =mallocdoublevector(12);

  for(nlap=1;nlap<=laps;nlap++)
  {
    /*af->nlaps=nlap;*/
    af->nlaps=1;

    sprintf(string,"LAP:%d/%d",nlap,laps);
    errormessage(string);
    if(fout!=NULL) fprintf(fout,"%s\n",string);
    setlaps((wmenu.childs+2)->hwnd,nlap,laps);

    memory1=availablephysicalmemory("REMAIN:");  /*MEMORY AVAILABLE*/

    for(i=1;i<=msize;i++)           /*GLOBAL MATRIX INITIALIZATION.*/
    {
      g=(gmtx+(i-1))->down; /*NEXT OF DIAGONAL.*/
      while(g!=NULL) /*CLEAR ROW.*/
      {
        p=g;
        g=g->down;
        free(p);
      }

      ginit.m=(unsigned short int)i;
      /*ginit.n=(unsigned short int)i;*/

      *(gmtx+(i-1))=ginit;
    }
    comps=msize; /*INITIAL COMPONENTS=DIAGONALS.*/

    laptime("ASSEMBLING GLOBAL MATRIX.",t0);

    for(i=1;i<=nelem;i++)               /*ASSEMBLAGE GLOBAL MATRIX.*/
    {
      inputelem(elems,felem,i-1,&elem);        /*READ ELEMENT DATA.*/
      for(ii=0;ii<=1;ii++)
      {
        for(jj=0;jj<=5;jj++)
        {
          (elems+i-1)->iconf[ii][jj]=elem.iconf[ii][jj];
        }
      }
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

      directioncosineII(elem.node[0]->d[0],
                        elem.node[0]->d[1],
                        elem.node[0]->d[2],
                        elem.node[1]->d[0],
                        elem.node[1]->d[1],
                        elem.node[1]->d[2],
                        elem.cangle,
                        drccos);                         /*[DRCCOS]*/

      transmatrixII(drccos,tmatrix);       /*TRANSFORMATION MATRIX.*/
      assememtxII(elem,estiff);        /*ELASTIC MATRIX OF ELEMENT.*/
      estiff=modifyhinge(elem,estiff);             /*MODIFY MATRIX.*/
      estiff=assempmtx(elem,estiff);          /*ADD PLASTIC MATRIX.*/
      transformationII(estiff,tmatrix,e,t);        /*[K]=[Tt][k][T]*/

      assemgstiffness(gmtx,estiff,&elem);             /*ASSEMBLAGE.*/
    }
    sprintf(string,"GLOBAL MATRIX %ld COMPS ASSEMBLED.",comps);
    laptime(string,t0);

    overlayhdc(*(wdraw.childs+1),SRCPAINT);       /*UPDATE DISPLAY.*/

    for(i=0;i<3;i++)
    {
      if(acc[i].a==NULL) dacc[i]=0.0;
      else
      {
        dacc[i]=accelerationincrement(acc[i],ddt,nlap);
        dacc[i]*=afact;
      }
    }
    assemaccel(gacc,dacc,af);                /*ACCELERATION VECTOR.*/

    /*NEWMARK'S BETA PROCESS.*/
    laptime("NEWMARK'S BETA.",t0);
    determinant=newmarkbeta(gmtx,mmtx,gacc,confs,nnode,
                            u,ud,udd,du,dud,dudd,
                            ddt,(1.0/4.0));

    sprintf(string,"DETERMINANT=%.5E COMPS=%ld",determinant,comps);
    errormessage(string);
    if(fout!=NULL) fprintf(fout,"%s\n",string);

    /*
    safety=nlap*dsafety;
    sprintf(string,"SAFETY FACTOR=%.5f",safety);
    errormessage(string);
    if(fout!=NULL) fprintf(fout,"%s\n",string);
    */

    if(determinant<=0.0)
    {
      errormessage(" ");
      errormessage("INSTABLE TERMINATION.");
      if(fout!=NULL) fprintf(fout,"INSTABLE TERMINATION.\n");

      laptime("\0",t0);

      fclose(fin);

      gfree(gmtx,nnode); /*FREE GLOBAL MATRIX.*/
      free(gacc);
      /*free(confs);*/
      freematrix(drccos,3);
      freematrix(estiff, 12);
      freematrix(tmatrix,12);
      freematrix(e,12);
      freematrix(t,12);
      free(estress);
      free(gdisp);
      free(edisp);

      memory2=availablephysicalmemory("REMAIN:");
      sprintf(string,"CONSUMPTION:%ld[BYTES]",(memory0-memory2));
      errormessage(string);

      return 1;
    }

    laptime("OUTPUT INTO FILE.",t0);

    if(fout!=NULL) fprintf(fout,"\"DISPLACEMENT\"\n");
    outputdisp(du,fout,nnode,nodes);  /*INCREMENTAL DISPLACEMENT.*/
    /*while(!GetAsyncKeyState(VK_LBUTTON))
    ;*/                                   /*LEFT CLICK TO CONTINUE.*/

    if(fout!=NULL) fprintf(fout,"\"STRESS\"\n");
    for(i=1;i<=nelem;i++)                   /*STRESS OUTPUT,UPDATE.*/
    {
      inputelem(elems,felem,i-1,&elem);

      inputnode(fdisp,elem.node[0]);
      inputnode(fdisp,elem.node[1]);

      elem.sect=(elems+i-1)->sect;             /*READ SECTION DATA.*/

      elemstressII(estress,
                   &elem,du,felem,fout,
                   drccos,tmatrix,estiff,gdisp,edisp,func);

      outputstress(elem,estress,fout,func);
    }
    if(wsurf.hwnd!=NULL)
    {
      drawyieldsurface((wsurf.childs+1)->hdcC,
                       (wsurf.childs+1)->vparam,
                       SURFACEX,SURFACEY,SURFACEZ,
                       af->fsurface);
      overlayhdc(*(wsurf.childs+1),SRCPAINT);     /*UPDATE DISPLAY.*/
    }
    /*while(!GetAsyncKeyState(VK_LBUTTON))
    ;*/                                   /*LEFT CLICK TO CONTINUE.*/

    if(fout!=NULL) fprintf(fout,"\"REACTION\"\n");
    outputreaction(gmtx,du,nodes,confs,freact,fout,nnode);

    updateform(fdisp,du,nnode);                 /*FORMATION UPDATE.*/
    if(fout!=NULL) fprintf(fout,"\"CURRENT FORM\"\n");
    for(ii=0;ii<nnode;ii++)
    {
      sprintf(string,"NODE:%5ld {U}=",(nodes+ii)->code);
      for(jj=1;jj<=6;jj++)
      {
        loffset=6*ii+jj;
        vread(fdisp,loffset,&data);
        sprintf(s," %14.5f",data);
        strcat(string,s);
      }
      if(fout!=NULL) fprintf(fout,"%s\n",string);
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

        gfree(gmtx,nnode); /*FREE GLOBAL MATRIX.*/
        gfree(mmtx,nnode); /*FREE MASS   MATRIX.*/
        free(gacc);
        /*free(confs);*/
        freematrix(drccos,3);
        freematrix(estiff, 12);
        freematrix(tmatrix,12);
        freematrix(e,12);
        freematrix(t,12);
        free(estress);
        free(gdisp);
        free(edisp);

        errormessage(" ");
        errormessage("ABORTED.");
        if(fout!=NULL) fprintf(fout,"ABORTED.\n");

        laptime("\0",t0);
        return 1;
      }
      t2=clock();
      time=(t2-t1)/CLK_TCK;
      if(time>=WAIT) break;               /*CONTINUE AFTER WAITING.*/
    }
  }/*END OF LAP.REPEAT UNTIL INSTABLE.*/

  if((wdraw.childs+1)->hdcC!=NULL &&
     felem!=NULL && fdisp!=NULL)                 /*DRAW LAST FRAME.*/
  {
    for(i=1;i<=nelem;i++)
    {
      inputelem(elems,felem,i-1,&elem);
      for(ii=0;ii<=1;ii++) /*COPY HINGE DATA.*/
      {
        for(jj=0;jj<=5;jj++)
        {
          (elems+i-1)->iconf[ii][jj]=elem.iconf[ii][jj];
        }
      }

      inputnode(fdisp,elem.node[0]);
      inputnode(fdisp,elem.node[1]);

      drawglobalwire((wdraw.childs+1)->hdcC,
                     (wdraw.childs+1)->vparam,
                     *af,elem,255,255,255,
                              255,255,255,0,ONSCREEN);
    }
    overlayhdc(*(wdraw.childs+1),SRCPAINT);       /*UPDATE DISPLAY.*/
  }

  fclose(fin);

  gfree(gmtx,nnode); /*FREE GLOBAL MATRIX.*/
  gfree(mmtx,nnode); /*FREE MASS   MATRIX.*/
  /*free(confs);*/
  freematrix(drccos,3);
  freematrix(estiff, 12);
  freematrix(tmatrix,12);
  freematrix(e,12);
  freematrix(t,12);
  free(estress);
  free(gdisp);
  free(edisp);

  af->eigenvec=(double **)malloc(1*sizeof(double *));
  *((af->eigenvec)+0)=du;

  errormessage(" ");
  errormessage("COMPLETED.");
  if(fout!=NULL) fprintf(fout,"COMPLETED.\n");

  memory2=availablephysicalmemory("REMAIN:");
  sprintf(string,"CONSUMPTION:%ld[BYTES]",(memory0-memory2));
  errormessage(string);
  errormessage(" ");

  return 0;
}/*gnshn101*/

double *inputacceleration(FILE *fin,
                          long int *ndata,double *dt)
/*INPUT ACCELERATION FROM TEXTFILE.*/
{
  /*CAPTION*/
  /*NUMBER OF DATA*/
  /*INCREMENT OF TIME [s]*/
  /*ACCELERATION DATA [cm/s2]*/

  char **data,str[256];
  int i,nstr;
  long int idata;
  double *a;

  if(fin==NULL) return NULL;

  fseek(fin,0L,SEEK_SET);

  fgets(str,256,fin); /*CAPTION.*/

  data=fgetsbrk(fin,&nstr);
  if(nstr<2) return NULL;
  if(!strcmp(*(data+0),"DATA"))
  {
    *ndata=strtol(*(data+1),NULL,10); /*NUMBER OF DATA.*/
  }

  data=fgetsbrk(fin,&nstr);
  if(nstr<2) return NULL;
  if(!strcmp(*(data+0),"DT"))
  {
    *dt=strtod(*(data+1),NULL); /*INCREMENT OF TIME.*/
  }

  a=(double *)malloc((*ndata)*sizeof(double));

  idata=0;
  while(1)
  {
    data=fgetsbrk(fin,&nstr);
    if(nstr==0) return a;

    for(i=0;i<nstr;i++)
    {
      *(a+idata)=strtod(*(data+i),NULL);

      idata++;
      if(idata>=(*ndata)) return a;
    }

    freestr(data,nstr);
  }

}/*inputacceleration*/

void assemmass(struct gcomponent *mmtx,
               struct arclmframe *af)
/*ASSEMBLAGE MASS MATRIX.*/
{
  long int i,j,ii;
  double gdata;

  for(i=0;i<(af->nnode);i++)
  {
    for(j=0;j<=2;j++) /*FOR X,Y,Z LINE.*/
    {
      ii=6*i+j+1;

      if((af->confs+ii-1)->iconf==0)
      {
        gdata=*(af->nmass+i);
        if(gdata!=0.0) gwrite(mmtx,ii,ii,gdata);
      }
    }
  }
  return;
}/*assemmass*/

double accelerationincrement(struct accdata acc,
                             double ddt,long int lap)
/*INCREMENT OF ACCELERATION.*/
{
  long int n;
  double t,a1,a2,da1,da2;
  double dacc;

  t=ddt*(double)(lap-1);

  n=(long int)(t/(acc.dt));
  a1=*(acc.a+n);
  a2=*(acc.a+n+1);
  da1=a1+(a2-a1)*(t-((double)n)*(acc.dt))/(acc.dt);

  t+=ddt;
  n=(long int)(t/(acc.dt));
  a1=*(acc.a+n);
  a2=*(acc.a+n+1);
  da2=a1+(a2-a1)*(t-((double)n)*(acc.dt))/(acc.dt);

  dacc=da2-da1;

  return dacc;
}/*accelerationincrement*/

void assemaccel(double *gvct,
                double dacc[],struct arclmframe *af)
/*ASSEMBLAGE ACCELERATION INTO GLOBAL VECTOR.*/
{
  int i,j,ii;

  for(i=0;i<6*(af->nnode);i++)
  {
    *(gvct+i)=0.0; /*INITIALIZATION.*/
  }

  for(i=0;i<(af->nnode);i++)
  {
    for(j=0;j<=2;j++) /*FOR X,Y,Z LINE.*/
    {
      ii=6*i+j;

      if((af->confs+ii)->iconf==0)
      {
        *(gvct+ii)=dacc[j];
      }
    }
  }
  return;
}/*assemaccel*/

double newmarkbeta(struct gcomponent *gmtx,
                   struct gcomponent *mmtx,
                   double *gacc,
                   struct oconf *confs,
                   int nnode,
                   double *u, double *ud, double *udd,
                   double *du,double *dud,double *dudd,
                   double ddt,double beta)
/*NEWMARK'S BETA BY INCREMENT.*/
{
  long int i,j,msize;
  double gdata,mdata,det,sign;

  msize=6*nnode;

  for(i=0;i<msize;i++)
  {
    if((confs+i)->iconf==0)
    {
      *(du+i)=0.0;

      for(j=0;j<=i;j++)
      {
        if((confs+j)->iconf==0)
        {
          gread(gmtx,(i+1),(j+1),&gdata);
          gread(mmtx,(i+1),(j+1),&mdata);

          *(du+i)+=-mdata*(*(gacc+j))
                   +mdata*((*(ud +j))/beta/ddt
                          +(*(udd+j))/2.0/beta); /*{dP'}:(3.59)*/

          gdata+=mdata/beta/(ddt*ddt); /*[K']:(3.59)*/

          gwrite(gmtx,(i+1),(j+1),gdata);
        }
      }
    }
  }

  /*CROUT LU DECOMPOSITION. [K']{dU}={dP'}*/
  croutludecomposition(gmtx,du,confs,msize,&det,&sign); /*(3.58)*/
  if(sign<=0.0) return sign;

  for(i=0;i<msize;i++)
  {
    if((confs+i)->iconf==0)
    {
      *(dud+i)=(*(du+i))/2.0/beta/ddt-(*(ud+i))/2.0/beta
              -(1.0/4.0/beta-1.0)*(*(udd+i))*ddt; /*(3.61)*/

      (*(dudd+i))=(*(du+i))/beta/(ddt*ddt)-(*(ud+i))/beta/ddt
                  -(*(udd+i))/2.0/beta; /*(3.62)*/

      (*(u  +i))+=(*(du  +i));
      (*(ud +i))+=(*(dud +i));
      (*(udd+i))+=(*(dudd+i));
    }
  }

  return det;
}/*newmarkbeta*/


