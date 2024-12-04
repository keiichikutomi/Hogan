/*QADHG001.C FOR WIN32 SINCE 1995.11.24.JUNSATO.*/
/*LAST CHANGE:1997.6.28.*/

/*NODE:12 DOF.*/
/*LOAD:STATIC.*/

#define DU    0
#define DUDX  1
#define DUDY  2
#define DUDZ  3
#define DV    4
#define DVDX  5
#define DVDY  6
#define DVDZ  7
#define DW    8
#define DWDX  9
#define DWDY 10
#define DWDZ 11

#define FMATRIX 0
#define SMATRIX 1

/*FOR BIQUADRATIC WIRE,FILM========================================*/
struct node12{long int code,loff;  /*12DOF/NODE:U,DU/DX,DU/DY,DU/DZ*/
              double c[3];         /*           V,DV/DX,DV/DY,DV/DZ*/
              double d[12];};      /*           W,DW/DX,DW/DY,DW/DZ*/

struct node09{long int code,loff;     /*9DOF/NODE:DU/DX,DU/DY,DU/DZ*/
              double c[3];            /*          DV/DX,DV/DY,DV/DZ*/
              double d[9];};          /*          DW/DX,DW/DY,DW/DZ*/

struct qsect{long int code,loff;
             double E,poi;                              /*MATERIAL.*/
             double area,Ixx,Iyy,Jzz;          /*WIRE COEFFICIENTS.*/
             double thick;};
struct qwire{long int code,loff;
             double cangle;
             signed char iconf[3][9]; /*iconf[3][12] FOR NONLINEAR.*/
             double stress[3][12];
             struct node12 *(n12[2]);
             struct node09 *n09;
             struct qsect *sect;};    /*WIRE ELEM FOR ARCLM001,101.*/
struct qfilm{long int code,loff;
             signed char iconf[6][9];
             double stress[6][12];
             struct node12 *(n12[3]);
             struct node09 *(n09[3]);
             struct qsect *sect;
             double ndisp[63];};           /*FILM ELEM FOR ARCLM301.*/

struct biquadframe{long int code,loff;
                   char *appelation;
                   int nnode,nnmid,nwire,nfilm,nwsec,nfsec,nreact;
                   FILE *fdisp,*fwire,*ffilm,*freact,*fsurface;
                   double *ddisp;
                   struct node12 *nodes,*ninit;
                   struct node09 *nmids;
                   struct qsect *sects;
                   struct qwire *wires;
                   struct qfilm *films;
                   struct oconf *confs;
                  }; /*BIQUADRATIC FRAME.*/
/*-----------------------------------------------------------------*/
int arclm301(struct biquadframe *bqf);
int factorial(int n);
double *polypoly1(int n,double *c1,double *c2);
double *polypoly2(int n,double *c1,double *c2);
void bquadgfree(struct gcomponent *gmtx,int nnode,int nnmid);
void inputbiquadframe(FILE *ftext,struct biquadframe *bqf);
long int countmidnodes(struct biquadframe *bqf);
void createmidnodes(struct biquadframe *bqf);
void inputbquadinit(FILE *fin,int *nnode,
							  int *nwire,int *nfilm,
                              int *nwsec,int *nfsec);
void inputbquadnode12(FILE *fdisp,struct node12 *n12);
void inputbquadnode09(FILE *fdisp,int nnode,struct node09 *n09);
void inputwire(struct qwire *wires,FILE *fwire,int offset,
               struct qwire *wire);
void inputfilm(struct qfilm *films,FILE *ffilm,int offset,
               struct qfilm *film);
void initialbquadform(FILE *fdisp,struct node12 *nodes,int nnode,
                                  struct node09 *nmids,int nnmid);
void initialwire(struct qwire *wires,FILE *fwire,int nwire);
void initialfilm(struct qfilm *films,FILE *ffilm,int nfilm);
void initialbquadreact(FILE *fin,FILE *freact,int nreact);
double **bquadwiredrccos(double x1,double y1,double z1,
                         double x2,double y2,double z2,
                         double cangle);
double **bquadfilmdrccos(double x1,double y1,double z1,
                         double x2,double y2,double z2,
                         double x3,double y3,double z3);
double **assemwiremtx(struct qwire wire);
double **assemfilmmtx(struct qfilm film);
double **assemtrifilmmtx1(struct qfilm film,int mmode);
double **transmatrix12(double **drccos);
double **transmatrix09(double **drccos);
double **wiretransmatrix(double **drccos);
double **filmtransmatrix(double **drccos);
double **bquadtransformation(double **estiff,double **tmatrix,
                             int msize);
double **modifywirehinge(struct qwire wire,double **wstiff);
double **modifyfilmhinge(struct qfilm film,double **fstiff);
void assemwtog(struct gcomponent *gmtx,
               double **wstiff,struct qwire *wire,int nnode);
void assemftog(struct gcomponent *gmtx,
               double **fstiff,struct qfilm *film,int nnode);
double *extractwirevct(struct qwire wire,double *gvct,int nnode);
double *extractfilmvct(struct qfilm film,double *gvct,int nnode);
double *wirestress(struct qwire *wire,double *gvct,int nnode);
double *filmstress(struct qfilm *film,double *gvct,int nnode);
void bquadassemconf(struct oconf *confs,double *gvct,double dsafety,
                    int nnode,int nnmid);
void bquadmodifygivend(struct gcomponent *gmtx,double *gvct,
                       struct oconf *confs,int nnode,int nnmid);

void bquadoutputdisp(double *gvct,FILE *fout,int nnode,int nnmid,
                     struct node12 *nodes,struct node09 *nmids);
void bquadupdateform(FILE *fdisp,double *gvct,int nnode,int nnmid);
/*void bquadupdateform(double *ddisp,double *gvct,int nnode,int nnmid);*/
void outputwirestress(struct qwire wire,double *wstress,FILE *fout);
void outputfilmstress(struct qfilm film,double *fstress,FILE *fout);
void bquadoutputreaction(struct gcomponent *gmtx,
                         double *gvct,
                         struct node12 *nodes,
                         struct node09 *nmids,
                         struct oconf *confs,
                         FILE *freact,FILE *fout,
                         int nnode,int nnmid);

void drawbiquadframe(HDC hdc,struct viewparam vp,
					 struct biquadframe bf,
					 int mode);
void drawbiquadnodes(HDC hdc,struct viewparam vp,
					 struct biquadframe bf);
void drawbiquadnode12(HDC hdc,struct viewparam vp,struct node12 bn);
void drawbiquadnode09(HDC hdc,struct viewparam vp,struct node09 bn);
void drawbiquadwire(HDC hdc,struct viewparam vp,
                    struct qwire qw,
                    int mode);
void drawbiquadfilm(HDC hdc,struct viewparam vp,
                    struct qfilm qf,
                    int mode);
/*-----------------------------------------------------------------*/
FILE *fout;

int arclm301(struct biquadframe *bqf)
{
  DWORD memory0,memory1,memory2;

  FILE *fin;                                         /*FILE 8 BYTES*/
  FILE *fwire,*ffilm,*fdisp,*freact;
  char dir[]="\0";                                 /*DATA DIRECTORY*/
  char string[1000],str[80];
  int i,ii,jj;
  int nnode,nnmid,nwire,nfilm,nwsec,nfsec,nlap,nreact;
  long int msize;
  long int time;
  struct gcomponent ginit={0,0,0.0,NULL};
  struct gcomponent *gmtx,*g,*p;                    /*GLOBAL MATRIX*/
  double *gvct,data;                                     /*GLOBAL VECTOR*/
  double **drccos,**tmatrix,**wstiff,**fstiff,*wstress,*fstress;
  double determinant,sign,safety,dsafety=1.0;
  clock_t t0,t1,t2;

  struct qsect *sects;
  struct node12 node,*nodes,*ninit;
  struct node09 *nmids;
  struct qwire wire,*wires;
  struct qfilm film,*films;
  struct oconf *confs;

  int il,ir;
  double **mtxtest1,**mtxtest2;

  memory0=availablephysicalmemory("INITIAL:");   /*MEMORY AVAILABLE*/

  fin=fgetstofopen(dir,"r",ID_INPUTFILE);         /*OPEN INPUT FILE*/
  if(fin==NULL)
  {
    errormessage("ACCESS IMPOSSIBLE.");
    return 0;
  }
  fout=fgetstofopen("\0","w",ID_OUTPUTFILE);          /*OUTPUT FILE*/

  /*fgets(string,256,fin);*/                    /*INPUT APPELATION.*/
  /*errormessage(string);*/

  t0=clock();                                        /*CLOCK BEGIN.*/

  inputbquadinit(fin,&nnode,&nwire,&nfilm,&nwsec,&nfsec);/*INITIAL.*/
  sprintf(string,"NODES=%d WIRES=%d FILMS=%d SECTS=%d",
          nnode,nwire,nfilm,(nwsec+nfsec));
  errormessage(string);
  fprintf(fout,"%s\n",string);

  fdisp=fopen("canbin.dsp","wb+");  /*DISPLACEMENT:12,9 DIRECTIONS.*/
  bqf->fdisp=fdisp;

  fwire=fopen("canbin.wir","wb+");                    /*WIRES FILE.*/
  bqf->fwire=fwire;

  ffilm=fopen("canbin.fim","wb+");                    /*FILMS FILE.*/
  bqf->ffilm=ffilm;

  free(bqf->sects);
  free(bqf->nodes);
  free(bqf->ninit);
  free(bqf->nmids);
  free(bqf->wires);
  free(bqf->films);
  free(bqf->confs);

  sects=(struct qsect *)malloc((nwsec+nfsec)*sizeof(struct qsect));
  if(sects==NULL) return 0;
  nodes=(struct node12 *)malloc(nnode*sizeof(struct node12));
  if(nodes==NULL) return 0;
  ninit=(struct node12 *)malloc(nnode*sizeof(struct node12));
  if(ninit==NULL) return 0;
  if(nwire>0)
  {
    wires=(struct qwire *)malloc(nwire*sizeof(struct qwire));
    if(wires==NULL) return 0;
  }
  if(nfilm>0)
  {
    films=(struct qfilm *)malloc(nfilm*sizeof(struct qfilm));
    if(films==NULL) return 0;
  }
  confs=(struct oconf *)malloc(12*(nnode)*sizeof(struct oconf));
  if(confs==NULL) return 0;
  errormessage("MALLOC COMPLETED.");

  bqf->sects=sects;
  bqf->nodes=nodes;
  bqf->ninit=ninit;
  if(nwire>0) bqf->wires=wires;
  if(nfilm>0) bqf->films=films;
  bqf->confs=confs;

  inputbiquadframe(fin,bqf);        /*READY TO READ LONG REACTIONS.*/

/*
for(ii=0;ii<nnode;ii++)
{
  fprintf(fout,"Node %d CONF=",(bqf->nodes+ii)->code);
  for(jj=0;jj<12;jj++)
  {
    fprintf(fout," %d",(confs+12*ii+jj)->iconf);
  }
  fprintf(fout,"\n");
}
*/

  nnode=bqf->nnode;
  nwire=bqf->nwire;
  nfilm=bqf->nfilm;
  nwsec=bqf->nwsec;
  nfsec=bqf->nfsec;
  nreact=bqf->nreact;

  errormessage("COUNT REACTIONS.");
  freact=fopen("canbin.rct","wb+");                /*REACTION FILE.*/
  if(freact==NULL) return 0;
  bqf->freact=freact;

  initialbquadreact(fin,freact,nreact);     /*ASSEMBLAGE REACTIONS.*/

  errormessage("CREATE MID NODES.");
  bqf->nnmid=countmidnodes(bqf);
  nnmid=bqf->nnmid;
  bqf->nmids=(struct node09 *)malloc(nnmid*sizeof(struct node09));
  bqf->confs=(struct oconf *)realloc(confs,(12*nnode+9*nnmid)*sizeof(struct oconf));
  confs=bqf->confs;

  createmidnodes(bqf);
  nnmid=bqf->nnmid;
  nmids=bqf->nmids;

/*
for(ii=0;ii<nnode;ii++)
{
  fprintf(fout,"Node %d CONF=",(bqf->nodes+ii)->code);
  for(jj=0;jj<12;jj++)
  {
    fprintf(fout," %d",(confs+12*ii+jj)->iconf);
  }
  fprintf(fout,"\n");
}
for(ii=0;ii<nnmid;ii++)
{
  fprintf(fout,"Nmid %d CONF=",(bqf->nmids+ii)->code);
  for(jj=0;jj<9;jj++)
  {
    fprintf(fout," %d",(confs+12*nnode+9*ii+jj)->iconf);
  }
  fprintf(fout,"\n");
}
*/

  msize=12*nnode+9*nnmid;                  /*SIZE OF GLOBAL MATRIX.*/

  errormessage("MALLOC MATRIX.");
  gmtx=(struct gcomponent *)          /*DIAGONALS OF GLOBAL MATRIX.*/
        malloc(msize*sizeof(struct gcomponent));
  gvct=(double *)malloc(msize*sizeof(double));      /*GLOBAL VECTOR*/
  if(gmtx==NULL || gvct==NULL) return 0;
  for(i=0;i<=msize-1;i++)
  {
    (gmtx+i)->down=NULL;            /*GLOBAL MATRIX INITIALIZATION.*/
    *(gvct+i)=0.0;                  /*GLOBAL VECTOR INITIALIZATION.*/
  }

  bqf->ddisp=(double *)malloc(msize*sizeof(double)); /*DISPLACEMENTS.*/
  if(bqf->ddisp==NULL) return 0;
  for(i=0;i<msize;i++) *(bqf->ddisp+i)=0.0;

  initialbquadform(fdisp,nodes,nnode,nmids,nnmid);  /*INITIAL FORM.*/
  if(nwire>0) initialwire(wires,fwire,nwire);   /*ASSEMBLAGE WIRES.*/
  if(nfilm>0) initialfilm(films,ffilm,nfilm);   /*ASSEMBLAGE FILMS.*/

/*
for(i=0;i<nnode;i++)
{
  node.loff=i;
  inputbquadnode12(fdisp,&node);
  fprintf(fout,"Node%d (%3.1f,%3.1f,%3.1f)\n",
          i+1,
          node.d[0],node.d[4],node.d[8]);
}
for(i=0;i<nfilm;i++)
{
  inputfilm(films,ffilm,i,&film);
  fprintf(fout,"FILM%d %ld(%ld,%ld,%ld)\n",
          i+1,film.code,
          film.n12[0]->code,
          film.n12[1]->code,
          film.n12[2]->code);
}
*/

  GetAsyncKeyState(VK_LBUTTON);                   /*CLEAR KEY LEFT.*/
  GetAsyncKeyState(VK_RBUTTON);                  /*CLEAR KEY RIGHT.*/

  drawglobalaxis((wdraw.childs+1)->hdcC,(wdraw.childs+1)->vparam,
                 0,0,255);                      /*DRAW GLOBAL AXIS.*/
  /*drawbiquadframe((wdraw.childs+1)->hdcC,
                  (wdraw.childs+1)->vparam,*bqf,ONSCREEN);*/
  /*overlayhdc(*(wdraw.childs+1),SRCPAINT);*/         /*UPDATE DISPLAY.*/

  for(nlap=1;nlap<=1 /*LAPS*/;nlap++)
  {
    sprintf(string,"LAP:%d",nlap);
    errormessage(string);
    fprintf(fout,"%s\n",string);

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

    laptime("ASSEMBLING WIRE MATRIX.",t0);
    for(i=0;i<nwire;i++)     /*ASSEMBLAGE WIRES INTO GLOBAL MATRIX.*/
    {
      inputwire(wires,fwire,i,&wire);          /*READ ELEMENT DATA.*/
      for(ii=0;ii<3;ii++)
      {
        for(jj=0;jj<9;jj++)
        {
          (wires+i)->iconf[ii][jj]=wire.iconf[ii][jj];
        }
      }

      inputbquadnode12(fdisp,wire.n12[0]);                   /*HEAD*/
      inputbquadnode12(fdisp,wire.n12[1]);                   /*TAIL*/
      inputbquadnode09(fdisp,nnode,wire.n09);                 /*MID*/

      /*wire.sect=(wires+i)->sect;*/           /*READ SECTION DATA.*/

      drccos=bquadwiredrccos(wire.n12[0]->d[0],
                             wire.n12[0]->d[4],
                             wire.n12[0]->d[8],
                             wire.n12[1]->d[0],
                             wire.n12[1]->d[4],
                             wire.n12[1]->d[8],
                             wire.cangle);      /*DIRECTION COSINE.*/

      tmatrix=wiretransmatrix(drccos);     /*TRANSFORMATION MATRIX.*/
      wstiff=assemwiremtx(wire);             /*ELASTIC WIRE MATRIX.*/

/*
 for(ii=0;ii<2;ii++)
 {
   fprintf(fout,"Node%d %ld(%3.1f,%3.1f,%3.1f)\n",
           ii+1,
           wire.n12[ii]->code,
           wire.n12[ii]->d[0],
           wire.n12[ii]->d[4],
           wire.n12[ii]->d[8]);
   fprintf(fout,"         (%3.1f,%3.1f,%3.1f)\n",
           wire.n12[ii]->d[1],
           wire.n12[ii]->d[5],
           wire.n12[ii]->d[9]);
   fprintf(fout,"         (%3.1f,%3.1f,%3.1f)\n",
           wire.n12[ii]->d[2],
           wire.n12[ii]->d[6],
           wire.n12[ii]->d[10]);
   fprintf(fout,"         (%3.1f,%3.1f,%3.1f)\n",
           wire.n12[ii]->d[3],
           wire.n12[ii]->d[7],
           wire.n12[ii]->d[11]);
 }
*/
/*
 for(ii=0;ii<3;ii++)
 {
   sprintf(string,"\0");
   for(jj=0;jj<3;jj++)
   {
     sprintf(str," %5.3f",*(*(drccos+ii)+jj));
     strcat(string,str);
   }
   fprintf(fout,"%s\n",string);
 }
 for(ii=0;ii<33;ii++)
 {
   sprintf(string,"\0");
   for(jj=0;jj<=ii;jj++)
   {
     sprintf(str," %1.0f",(*(*(wstiff+ii)+jj))/(*(*(wstiff+ii)+jj)));
     strcat(string,str);
   }
   fprintf(fout,"%s\n",string);
 }
*/
      wstiff=modifywirehinge(wire,wstiff);
      wstiff=bquadtransformation(wstiff,tmatrix,
                                 33);              /*[K]=[Tt][k][T]*/
      assemwtog(gmtx,wstiff,&wire,nnode); /*ASSEMBLAGE.*/

      /*modifycmq(&wire);*/
      /*assemcmq(wire,tmatrix,gvct);*/ /*ASSEMBLAGE CMQ AS LOAD.*/

      for(ii=0;ii<3;ii++) free(*(drccos+ii));
      free(drccos);
      for(ii=0;ii<33;ii++)
      {
        free(*(tmatrix+ii));
        free(*(wstiff+ii));
      }
      free(tmatrix);
      free(wstiff);
    }

    laptime("ASSEMBLING FILM MATRIX.",t0);
    for(i=0;i<nfilm;i++)     /*ASSEMBLAGE FILMS INTO GLOBAL MATRIX.*/
    {
      inputfilm(films,ffilm,i,&film);          /*READ ELEMENT DATA.*/

/*
fprintf(fout,"FILM%d %ld(%ld,%ld,%ld)\n",
        i+1,film.code,
        film.n12[0]->code,
        film.n12[1]->code,
        film.n12[2]->code);
*/

      for(ii=0;ii<6;ii++)
      {
        for(jj=0;jj<9;jj++)
        {
          (films+i)->iconf[ii][jj]=film.iconf[ii][jj];
        }
      }

      inputbquadnode12(fdisp,film.n12[0]);                  /*NODE1*/
      inputbquadnode12(fdisp,film.n12[1]);                  /*NODE2*/
      inputbquadnode12(fdisp,film.n12[2]);                  /*NODE3*/
      inputbquadnode09(fdisp,nnode,film.n09[0]);             /*MID1*/
      inputbquadnode09(fdisp,nnode,film.n09[1]);             /*MID2*/
      inputbquadnode09(fdisp,nnode,film.n09[2]);             /*MID3*/

      /*film.sect=(films+i)->sect;*/           /*READ SECTION DATA.*/

      drccos=bquadfilmdrccos(film.n12[0]->d[0],
                             film.n12[0]->d[4],
                             film.n12[0]->d[8],
                             film.n12[1]->d[0],
                             film.n12[1]->d[4],
                             film.n12[1]->d[8],
                             film.n12[2]->d[0],
                             film.n12[2]->d[4],
                             film.n12[2]->d[8]);         /*[DRCCOS]*/
/*
 for(ii=0;ii<3;ii++)
 {
   fprintf(fout,"Node%d %ld(%3.1f,%3.1f,%3.1f)\n",
           ii+1,
           film.n12[ii]->code,
           film.n12[ii]->d[0],
           film.n12[ii]->d[4],
           film.n12[ii]->d[8]);
   fprintf(fout,"         (%3.1f,%3.1f,%3.1f)\n",
           film.n12[ii]->d[1],
           film.n12[ii]->d[5],
           film.n12[ii]->d[9]);
   fprintf(fout,"         (%3.1f,%3.1f,%3.1f)\n",
           film.n12[ii]->d[2],
           film.n12[ii]->d[6],
           film.n12[ii]->d[10]);
   fprintf(fout,"         (%3.1f,%3.1f,%3.1f)\n",
           film.n12[ii]->d[3],
           film.n12[ii]->d[7],
           film.n12[ii]->d[11]);
 }
*/

      tmatrix=filmtransmatrix(drccos);     /*TRANSFORMATION MATRIX.*/
      /*fstiff=assemfilmmtx(film);*/             /*ELASTIC FILM MATRIX.*/
fstiff=assemtrifilmmtx1(film,FMATRIX); /*TRIANGLE LINEAR FILM MATRIX.*/
      if(fstiff==NULL) return 0;

 /*for(ii=0;ii<3;ii++)
 {
   sprintf(string,"\0");
   for(jj=0;jj<3;jj++)
   {
     sprintf(str," %5.3f",*(*(drccos+ii)+jj));
     strcat(string,str);
   }
   fprintf(fout,"%s\n",string);
 }*/
 /*for(ii=0;ii<63;ii++)
 {
   sprintf(string,"\0");
   for(jj=0;jj<=ii;jj++)
   {
     sprintf(str," %8.1E",*(*(fstiff+ii)+jj));
     strcat(string,str);
   }
   fprintf(fout,"%s\n",string);
 }*/

      fstiff=modifyfilmhinge(film,fstiff);
      fstiff=bquadtransformation(fstiff,tmatrix,
                                 63);              /*[K]=[Tt][k][T]*/
/*
 for(ii=0;ii<63;ii++)
 {
   sprintf(string,"\0");
   for(jj=0;jj<=ii;jj++)
   {
     if(*(*(fstiff+ii)+jj)==0.0) sprintf(str," %8.1f",*(*(fstiff+ii)+jj));
     else sprintf(str," %8.1E",*(*(fstiff+ii)+jj));
     strcat(string,str);
   }
   fprintf(fout,"%s\n",string);
 }
*/

      assemftog(gmtx,fstiff,&film,nnode); /*ASSEMBLAGE.*/

      /*modifycmq(&film);*/
      /*assemcmq(film,tmatrix,gvct);*/ /*ASSEMBLAGE CMQ AS LOAD.*/

      for(ii=0;ii<3;ii++) free(*(drccos+ii));
      free(drccos);
      for(ii=0;ii<63;ii++)
      {
        free(*(tmatrix+ii));
        free(*(fstiff+ii));
      }
      free(tmatrix);
      free(fstiff);
    }
    sprintf(string,"GLOBAL MATRIX %ld COMPS ASSEMBLED.",comps);
    laptime(string,t0);

    /*overlayhdc(*(wdraw.childs+1));*/            /*UPDATE DISPLAY.*/

    /*errormessage("ASSEMBLAGE GLOBAL VECTOR.");*/
    bquadassemconf(confs,gvct,dsafety,nnode,nnmid);/*GLOBAL VECTOR.*/
	bquadmodifygivend(gmtx,gvct,confs,nnode,nnmid);   /*COMPULSORY.*/

for(ii=12*nnode;ii<msize;ii++) /*PASS MID NODE.*/
{
  (confs+ii)->iconf=1;
}

/*
for(ii=0;ii<nnode;ii++)
{
  fprintf(fout,"Node %d CONF=",(bqf->nodes+ii)->code);
  for(jj=0;jj<12;jj++)
  {
    fprintf(fout," %d",(confs+12*ii+jj)->iconf);
  }
  fprintf(fout,"\n");
}
for(ii=0;ii<nnmid;ii++)
{
  fprintf(fout,"Nmid %d CONF=",(bqf->nmids+ii)->code);
  for(jj=0;jj<9;jj++)
  {
    fprintf(fout," %d",(confs+12*nnode+9*ii+jj)->iconf);
  }
  fprintf(fout,"\n");
}
*/

/*
fprintf(fout,"MSIZE=%d\n",msize);
for(ii=0;ii<msize;ii++)
{
  sprintf(string,"\0");
  for(jj=0;jj<=ii;jj++)
  {
    gread(gmtx,ii+1,jj+1,&data);
    if(data==0.0) sprintf(str," %8.1f",data);
    else sprintf(str," %8.1E",data);
    strcat(string,str);
  }
  fprintf(fout,"%s ICONF=%d VCONF=%7.3f\n",
          string,(confs+ii)->iconf,(confs+ii)->value);
}
*/

/*
fprintf(fout,"[M],{F}\n");
mtxtest1=(double **)malloc(5*sizeof(double *));
for(ii=0;ii<5;ii++) *(mtxtest1+ii)=(double *)malloc(5*sizeof(double));
il=0;
for(ii=0;ii<msize;ii++)
{
  if((confs+ii)->iconf==0)
  {
    sprintf(string,"\0");
    ir=0;
    for(jj=0;jj<msize;jj++)
    {
      if((confs+jj)->iconf==0)
      {
        gread(gmtx,ii+1,jj+1,&data);
        *(*(mtxtest1+il)+ir)=data;

        if(data==0.0) sprintf(str," %8.1f",data);
        else sprintf(str," %8.1E",data);
        strcat(string,str);
        ir++;
      }
    }
    fprintf(fout,"%s VCONF=%7.3f GVCT=%7.3f\n",string,(confs+ii)->value,*(gvct+ii));
    il++;
  }
}
fprintf(fout,"[M]\n");
for(ii=0;ii<5;ii++)
{
  sprintf(string,"\0");
  for(jj=0;jj<5;jj++)
  {
    sprintf(str," %8.1E",*(*(mtxtest1+ii)+jj));
    strcat(string,str);
  }
  fprintf(fout,"%s\n",string);
}
*/

    laptime("CROUT LU DECOMPOSITION.",t0);
    croutludecomposition(gmtx,
                         gvct,confs,
                         (12*nnode+9*nnmid),
                         &determinant,&sign);        /*[K]{dU}={dF}*/

    sprintf(string,"LOG{DETERMINANT}=%.5E SIGN=%.1f COMPS=%ld",
            determinant,sign,comps);
    errormessage(string);
    if(fout!=NULL) fprintf(fout,"%s\n",string);

    safety=nlap*dsafety;
    sprintf(string,"SAFETY FACTOR=%.5f",safety);
    errormessage(string);
    if(fout!=NULL) fprintf(fout,"%s\n",string);

/*
mtxtest2=matrixinverse(mtxtest1,5);
fprintf(fout,"[M-1]\n");
for(ii=0;ii<5;ii++)
{
  sprintf(string,"\0");
  for(jj=0;jj<5;jj++)
  {
    sprintf(str," %11.3E",*(*(mtxtest2+ii)+jj));
    strcat(string,str);
  }
  fprintf(fout,"%s\n",string);
}

fprintf(fout,"{U}\n");
for(ii=0;ii<msize;ii++)
{
  if((confs+ii)->iconf==0)
  {
    fprintf(fout," %12.3E\n",*(gvct+ii));
  }
}
*/

    if(determinant<=0.0)
    {
      errormessage(" ");
      errormessage("INSTABLE TERMINATION.");
      if(fout!=NULL) fprintf(fout,"INSTABLE TERMINATION.\n");

      laptime("\0",t0);

      fclose(fin);

      bquadgfree(gmtx,nnode,nnmid); /*FREE GLOBAL MATRIX.*/
      free(gvct);
      /*free(confs);*/

      memory2=availablephysicalmemory("REMAIN:");
      sprintf(string,"CONSUMPTION:%ld[BYTES]",(memory0-memory2));
      errormessage(string);

      return 1;
    }

	laptime("OUTPUT INTO FILE.",t0);

	if(fout!=NULL) fprintf(fout,"\"DISPLACEMENT\"\n");
	bquadoutputdisp(gvct,fout,nnode,nnmid,nodes,nmids);    /*DISPS.*/
	/*while(!GetAsyncKeyState(VK_LBUTTON))
	;*/                                   /*LEFT CLICK TO CONTINUE.*/

	if(fout!=NULL) fprintf(fout,"\"STRESS\"\n");
	for(i=0;i<nwire;i++)                    /*STRESS OUTPUT,UPDATE.*/
	{
	  inputwire(wires,fwire,i,&wire);          /*READ ELEMENT DATA.*/

	  inputbquadnode12(fdisp,wire.n12[0]);                   /*HEAD*/
	  inputbquadnode12(fdisp,wire.n12[1]);                   /*TAIL*/
	  inputbquadnode09(fdisp,nnode,wire.n09);                 /*MID*/

	  wire.sect=(wires+i)->sect;               /*READ SECTION DATA.*/

	  wstress=wirestress(&wire,gvct,nnode);

	  outputwirestress(wire,wstress,fout);
	  free(wstress);
	}
	for(i=0;i<nfilm;i++)                    /*STRESS OUTPUT,UPDATE.*/
	{
	  inputfilm(films,ffilm,i,&film);          /*READ ELEMENT DATA.*/

	  inputbquadnode12(fdisp,film.n12[0]);                  /*NODE1*/
	  inputbquadnode12(fdisp,film.n12[1]);                  /*NODE2*/
	  inputbquadnode12(fdisp,film.n12[2]);                  /*NODE3*/
	  inputbquadnode09(fdisp,nnode,film.n09[0]);             /*MID1*/
	  inputbquadnode09(fdisp,nnode,film.n09[1]);             /*MID2*/
	  inputbquadnode09(fdisp,nnode,film.n09[2]);             /*MID3*/

	  film.sect=(films+i)->sect;               /*READ SECTION DATA.*/

	  fstress=filmstress(&film,gvct,nnode);

	  outputfilmstress(film,fstress,fout);
	  free(fstress);
	}
	/*while(!GetAsyncKeyState(VK_LBUTTON))
	;*/                                   /*LEFT CLICK TO CONTINUE.*/

	if(fout!=NULL) fprintf(fout,"\"REACTION\"\n");
	bquadoutputreaction(gmtx,gvct,nodes,nmids,confs,
						freact,fout,nnode,nnmid);

	bquadupdateform(fdisp,gvct,nnode,nnmid); /*FORMATION UPDATE.*/
	/*bquadupdateform(bqf->ddisp,gvct,nnode,nnmid);*/
	/*if(fout!=NULL) fprintf(fout,"\"CURRENT FORM\"\n");
	for(ii=0;ii<nnode;ii++)
	{
	  sprintf(string,"NODE:%5ld {U}=",(nodes+ii)->code);
	  for(jj=1;jj<=12;jj++)
	  {
		loffset=12*ii+jj;
		vread(fdisp,loffset,&data);
		sprintf(s," %14.5f",data);
		strcat(string,s);
	  }
	  if(fout!=NULL) fprintf(fout,"%s\n",string);
	}*/
	for(i=0;i<nfilm;i++)                    /*STRESS OUTPUT,UPDATE.*/
	{
	  inputbquadnode12(fdisp,(films+i)->n12[0]);            /*NODE1*/
	  inputbquadnode12(fdisp,(films+i)->n12[1]);            /*NODE2*/
	  inputbquadnode12(fdisp,(films+i)->n12[2]);            /*NODE3*/
	  inputbquadnode09(fdisp,nnode,(films+i)->n09[0]);       /*MID1*/
	  inputbquadnode09(fdisp,nnode,(films+i)->n09[1]);       /*MID2*/
	  inputbquadnode09(fdisp,nnode,(films+i)->n09[2]);       /*MID3*/

	  for(ii=0;ii<3;ii++)
	  {
		for(jj=0;jj<12;jj++)
		{
		  (films+i)->ndisp[12*ii+jj]=*(gvct+12*(((films+i)->n12[ii])->loff)+jj);
		}
	  }
	  for(ii=0;ii<3;ii++)
	  {
		for(jj=0;jj<9;jj++)
		{
		  (films+i)->ndisp[36+9*ii+jj]=*(gvct+12*nnode+9*(((films+i)->n09[ii])->loff)+jj);
		}
	  }

	  drawbiquadfilm((wdraw.childs+1)->hdcC,
					 (wdraw.childs+1)->vparam,*(films+i),ONSCREEN);
	}
	overlayhdc(*(wdraw.childs+1),SRCPAINT);       /*UPDATE DISPLAY.*/

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

        bquadgfree(gmtx,nnode,nnmid); /*FREE GLOBAL MATRIX.*/
        free(gvct);
        /*free(confs);*/

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
  }                                        /*REPEAT UNTIL INSTABLE.*/

  fclose(fin);

  bquadgfree(gmtx,nnode,nnmid); /*FREE GLOBAL MATRIX.*/
  free(gvct);
  /*free(confs);*/

  errormessage(" ");
  errormessage("COMPLETED.");
  if(fout!=NULL) fprintf(fout,"COMPLETED.\n");

  memory2=availablephysicalmemory("REMAIN:");
  sprintf(string,"CONSUMPTION:%ld[BYTES]",(memory0-memory2));
  errormessage(string);

  return 0;
}/*arclm301*/

int factorial(int n)
/*RETURN:FACTORIAL OF N.*/
{
  int f=1;

  for(;n>=1;n--) f*=n;
  return f;
}/*factorial*/

double *polypoly1(int n,double *c1,double *c2)
/*MULTIPLICATION OF PERFECT n DIMENSIONAL POLYNOMIALS WITH 1 */
/*VARIABLES. RETURN:COEFFICIENTS OF THE ANSWER 2n DIMENSIONAL*/
/*POLYNOMIAL.*/
{
  int i,j,ncoeff;
  double *c;

  ncoeff=2*n+1;
  c=(double *)malloc(ncoeff*sizeof(double));
  for(i=1;i<=ncoeff;i++) *(c+i-1)=0.0;
  for(i=0;i<=n;i++)
  {
    for(j=0;j<=n;j++)
    {
      *(c+i+j)+=(*(c1+i))*(*(c2+j));
    }
  }
  return c;
}/*polypoly1*/

double *polypoly2(int n,double *c1,double *c2)
/*MULTIPLICATION OF PERFECT n DIMENSIONAL POLYNOMIALS WITH 2 */
/*VARIABLES. RETURN:COEFFICIENTS OF THE ANSWER 2n DIMENSIONAL*/
/*POLYNOMIAL.*/
{
  int i,ii,j,k,kk,iorder,ncoff1,ncoff2;
  double *c;

  ncoff1=(n+1)*(2*n+1);
  ncoff2=(n+1)*(n+2)/2;
  c=(double *)malloc(ncoff1*sizeof(double));
  for(i=1;i<=ncoff1;i++)
  {
    *(c+i-1)=0.0;
  }
  for(i=1;i<=ncoff2;i++)
  {
    j=i;
    iorder=sqrt(2.0*(i-1)+0.25)-0.5;
    k=1;
    for(ii=1;ii<ncoff2;)
    {
      for(kk=1;kk<=k;kk++)
      {
        *(c+j-1)+=(*(c1+i-1))*(*(c2+ii-1));
        ii++; j++;
      }
      j+=iorder;
      k++;
    }
  }
  return c;
}/*polypoly2*/

void bquadgfree(struct gcomponent *gmtx,int nnode,int nnmid)
/*FREE GLOBAL MATRIX.*/
{
  long int n,i;
  struct gcomponent *g,*p;

  n=12*nnode+9*nnmid;
  for(i=0;i<n;i++)
  {
    g=(gmtx+i)->down; /*NEXT OF DIAGONAL.*/
    while(g!=NULL) /*CLEAR ROW.*/
    {
      p=g;
      g=g->down;
      free(p);
    }
  }
  free(gmtx);

  return;
}/*bquadgfree*/

void inputbiquadframe(FILE *ftext,struct biquadframe *bqf)
/*TRANSLATE BQUAD INPUTFILE TEXT INTO MEMORY.*/
{
  char **data;
  int i,j,ii,jj,kk,k,n;
  long int loff;
  long int scode,code1,code2,code3;

  fseek(ftext,0L,SEEK_SET);
  inputbquadinit(ftext,&(bqf->nnode),
                       &(bqf->nwire),&(bqf->nfilm),
                       &(bqf->nwsec),&(bqf->nfsec));

  /*"INPUTFILE SPECIFICATION"*/
  /*APPELATION*/
  /*NNODE NWIRE NFILM NWSEC NFSEC*/
  /*ISECT E POI A Ixx Iyy Jzz*/
  /*ISECT E POI THICKNESS*/
  /*INODE X Y Z*/
  /*IWIRE ISECT NODEI J COORDANGLE BOUNDARY:27 LONGSTRESS:33*/
  /*IFILM ISECT NODEI J K BOUNDARY:54 LONGSTRESS:63*/
  /*INODE CONFINEMENTTYPE:12 CONFINEMENTVALUE:12*/
  /*INODE DIRECTION LONGREACTION*/

  for(i=0;i<(bqf->nwsec);i++) /*WIRE SECTION.*/
  {
    (bqf->sects+i)->loff=i;

    data=fgetsbrk(ftext,&n);

    (bqf->sects+i)->code=strtol(*(data+0),NULL,10);
    (bqf->sects+i)->E=strtod(*(data+1),NULL);
    (bqf->sects+i)->poi=strtod(*(data+2),NULL);
    (bqf->sects+i)->area=strtod(*(data+3),NULL);
    (bqf->sects+i)->Ixx=strtod(*(data+4),NULL);
    (bqf->sects+i)->Iyy=strtod(*(data+5),NULL);
    (bqf->sects+i)->Jzz=strtod(*(data+6),NULL);

    freestr(data,n);
  }
  for(i=bqf->nwsec;i<(bqf->nwsec+bqf->nfsec);i++) /*FILM SECTION.*/
  {
    (bqf->sects+i)->loff=i;

    data=fgetsbrk(ftext,&n);

    (bqf->sects+i)->code=strtol(*(data+0),NULL,10);
    (bqf->sects+i)->E=strtod(*(data+1),NULL);
    (bqf->sects+i)->poi=strtod(*(data+2),NULL);
    (bqf->sects+i)->thick=strtod(*(data+3),NULL);

    freestr(data,n);
  }
  for(i=0;i<(bqf->nnode);i++) /*NODES.*/
  {
    (bqf->nodes+i)->loff=i;

    data=fgetsbrk(ftext,&n);

    (bqf->nodes+i)->code=strtol(*(data+0),NULL,10);
    (bqf->nodes+i)->d[0]=strtod(*(data+1),NULL);
    (bqf->nodes+i)->d[4]=strtod(*(data+2),NULL);
    (bqf->nodes+i)->d[8]=strtod(*(data+3),NULL);

    *(bqf->ninit+i)=*(bqf->nodes+i);

    freestr(data,n);
  }
  for(i=0;i<(bqf->nwire);i++) /*WIRES.*/
  {
    (bqf->wires+i)->loff=i;

    data=fgetsbrk(ftext,&n);

    (bqf->wires+i)->code=strtol(*(data+0),NULL,10);
    scode=strtol(*(data+1),NULL,10); /*SECTION.*/
    code1=strtol(*(data+2),NULL,10); /*HEAD NODE.*/
    code2=strtol(*(data+3),NULL,10); /*TAIL NODE.*/
    (bqf->wires+i)->cangle=strtod(*(data+4),NULL);

    k=5;
    for(ii=0;ii<=2;ii++)     /*BOUNDARY OF ENDS,MID.0:RIGID 1:HINGE*/
    {
      for(kk=0;kk<=2;kk++) /*U,V,W*/
      {
        loff=4*kk;
        (bqf->wires+i)->iconf[ii][loff]=0; /*U,V,W*/
        for(jj=(loff+1);jj<=(loff+3);jj++) /*d/dX,d/dY,d/dZ*/
        {
          (bqf->wires+i)->iconf[ii][jj]
          =(signed char)strtol(*(data+k),NULL,10);
          k++;
        }
      }
    }
    for(ii=0;ii<2;ii++)                              /*ENDS STRESS.*/
    {
      for(jj=0;jj<12;jj++)
      {
        (bqf->wires+i)->stress[ii][jj]
        =strtod(*(data+k),NULL);
        k++;
      }
    }
    for(kk=0;kk<3;kk++)                               /*MID STRESS.*/
    {
      for(jj=1;jj<=3;jj++)
      {
        (bqf->wires+i)->stress[2][4*kk+jj]
        =strtod(*(data+k),NULL);
        k++;
      }
    }
    freestr(data,n);

    loff=0;
    for(ii=0;ii<1;)                                      /*SECTION.*/
    {
      if((bqf->sects+loff)->code==scode)
      {
        (bqf->wires+i)->sect=bqf->sects+loff;
        ii++;
      }
      loff++;
    }
	loff=0;
    for(ii=0;ii<2;)                                        /*NODES.*/
    {
      if((bqf->nodes+loff)->code==code1)
      {
        (bqf->wires+i)->n12[0]=bqf->nodes+loff;
        ii++;
      }
      if((bqf->nodes+loff)->code==code2)
      {
        (bqf->wires+i)->n12[1]=bqf->nodes+loff;
        ii++;
      }
      loff++;
    }
  }
  for(i=0;i<(bqf->nfilm);i++) /*FILMS.*/
  {
    (bqf->films+i)->loff=i;

    data=fgetsbrk(ftext,&n);

    (bqf->films+i)->code=strtol(*(data+0),NULL,10);
    scode=strtol(*(data+1),NULL,10); /*SECTION.*/
    code1=strtol(*(data+2),NULL,10); /*NODE1.*/
    code2=strtol(*(data+3),NULL,10); /*NODE2.*/
    code3=strtol(*(data+4),NULL,10); /*NODE3.*/

    k=5;
    for(ii=0;ii<=5;ii++)    /*BOUNDARY OF ENDS,MIDS.0:RIGID 1:HINGE*/
    {
      for(kk=0;kk<=2;kk++) /*U,V,W*/
      {
        loff=4*kk;
        (bqf->films+i)->iconf[ii][loff]=0; /*U,V,W*/
        for(jj=(loff+1);jj<=(loff+3);jj++) /*d/dX,d/dY,d/dZ*/
        {
          (bqf->films+i)->iconf[ii][jj]
          =(signed char)strtol(*(data+k),NULL,10);
          k++;
        }
      }
    }
    for(ii=0;ii<3;ii++)                              /*ENDS STRESS.*/
    {
      for(jj=0;jj<12;jj++)
      {
        (bqf->films+i)->stress[ii][jj]
        =strtod(*(data+k),NULL);
        k++;
      }
    }
    for(ii=0;ii<3;ii++)                              /*MIDS STRESS.*/
    {
      for(kk=0;kk<3;kk++) /*U,V,W*/
      {
        for(jj=1;jj<=3;jj++) /*d/dX,d/dY,d/dZ*/
        {
          (bqf->films+i)->stress[2][4*kk+jj]
          =strtod(*(data+k),NULL);
          k++;
        }
      }
    }
    freestr(data,n);

    loff=0;
    for(ii=0;ii<1;)                                      /*SECTION.*/
    {
      if((bqf->sects+loff)->code==scode)
      {
        (bqf->films+i)->sect=bqf->sects+loff;
        ii++;
      }
      loff++;
    }
    loff=0;
    for(ii=0;ii<3;)                                        /*NODES.*/
    {
      if((bqf->nodes+loff)->code==code1)
      {
        (bqf->films+i)->n12[0]=bqf->nodes+loff;
        ii++;
      }
      if((bqf->nodes+loff)->code==code2)
      {
        (bqf->films+i)->n12[1]=bqf->nodes+loff;
        ii++;
      }
      if((bqf->nodes+loff)->code==code3)
      {
        (bqf->films+i)->n12[2]=bqf->nodes+loff;
        ii++;
      }
      loff++;
    }

    for(ii=0;ii<63;ii++) (bqf->films+i)->ndisp[ii]=0.0;
  }

  bqf->nreact=0;
  for(i=0;i<(bqf->nnode);i++) /*CONF VECTOR:CONFINEMENT,VALUE.*/
  {
	data=fgetsbrk(ftext,&n);

    for(j=0;j<12;j++)
    {
      loff=12*i+j;

      (bqf->confs+loff)->iconf
      =(signed char)strtol(*(data+j+1),NULL,10);
      (bqf->confs+loff)->value
      =strtod(*(data+j+13),NULL);

      if((bqf->confs+loff)->iconf==1) (bqf->nreact)++;
    }

    freestr(data,n);
  }

  /*LONGREACTION:UNDER CONSTRUCTION.*/

  return;
}/*inputbiquadframe*/

long int countmidnodes(struct biquadframe *bqf)
/*COUNT MIDPOINT NODES.*/
{
  int i,j,ii,jj,kk,ll,soe[4]={0,1,2,0},count[3];
  int nnmid;
  long int loff;
  struct node09 mid,*nmids;
  struct oconf *confs;

  nnmid=bqf->nwire;

  for(i=0;i<(bqf->nfilm);i++) /*COUNT NMID.*/
  {
    count[0]=0;
    count[1]=0;
    count[2]=0;

    for(j=0;j<(bqf->nwire);j++)
    {
      for(ii=0;ii<3;ii++) /*FIND WIRE ON FILM LINE.*/
      {
        if(((bqf->films+i)->n12[soe[ii]]
            ==(bqf->wires+j)->n12[0] &&
            (bqf->films+i)->n12[soe[ii+1]]
            ==(bqf->wires+j)->n12[1]) ||
           ((bqf->films+i)->n12[soe[ii]]
            ==(bqf->wires+j)->n12[1] &&
            (bqf->films+i)->n12[soe[ii+1]]
            ==(bqf->wires+j)->n12[0]))
        {
          count[ii]=1;
        }
      }
    }

    for(j=0;j<i;j++)
    {
      for(ii=0;ii<3;ii++) /*FIND MID NODE ALREADY CREATED ON FILM.*/
      {
        if(count[ii]==0)
        {
          for(jj=0;jj<3;jj++)
          {
            if(((bqf->films+i)->n12[soe[ii]]
                ==(bqf->films+j)->n12[soe[jj]] &&
                (bqf->films+i)->n12[soe[ii+1]]
                ==(bqf->films+j)->n12[soe[jj+1]]) ||
               ((bqf->films+i)->n12[soe[ii]]
                ==(bqf->films+j)->n12[soe[jj+1]] &&
                (bqf->films+i)->n12[soe[ii+1]]
				==(bqf->films+j)->n12[soe[jj]]))
            {
              count[ii]=1;
            }
		  }
        }
      }
    }

    for(ii=0;ii<3;ii++) /*CASE NO WIRE ON FILM LINE.*/
    {
      if(count[ii]==0)
      {
        nnmid++;
      }
    }
  }

  return nnmid;
}/*countmidnodes*/

void createmidnodes(struct biquadframe *bqf)
/*CREATE MIDPOINT NODES.*/
{
  int i,j,ii,jj,kk,ll,soe[4]={0,1,2,0},count[3];
  int nnmid=0;
  long int loff;
  struct node09 mid,*nmids;
  struct oconf *confs;

  for(i=0;i<(bqf->nwire);i++)
  {
    mid.loff=nnmid;
    mid.code=nnmid+1; /*UNDER CONSTRUCTION.*/

    mid.c[GX]=0.5*((bqf->wires+i)->n12[0]->d[0]
                  +(bqf->wires+i)->n12[1]->d[0]);
    mid.c[GY]=0.5*((bqf->wires+i)->n12[0]->d[4]
                  +(bqf->wires+i)->n12[1]->d[4]);
    mid.c[GZ]=0.5*((bqf->wires+i)->n12[0]->d[8]
                  +(bqf->wires+i)->n12[1]->d[8]);

    *(bqf->nmids+nnmid)=mid;
    (bqf->wires+i)->n09=bqf->nmids+nnmid; /*MIDNODE OF WIRE.*/

    nnmid++;
  }

  for(i=0;i<(bqf->nfilm);i++)
  {
    count[0]=0;
    count[1]=0;
    count[2]=0;

    for(j=0;j<(bqf->nwire);j++)
    {
      for(ii=0;ii<3;ii++) /*FIND WIRE ON FILM LINE.*/
      {
        if(((bqf->films+i)->n12[soe[ii]]
            ==(bqf->wires+j)->n12[0] &&
            (bqf->films+i)->n12[soe[ii+1]]
            ==(bqf->wires+j)->n12[1]) ||
           ((bqf->films+i)->n12[soe[ii]]
            ==(bqf->wires+j)->n12[1] &&
            (bqf->films+i)->n12[soe[ii+1]]
            ==(bqf->wires+j)->n12[0]))
        {
          count[ii]=1;
          (bqf->films+i)->n09[ii]=(bqf->wires+j)->n09;
        }
      }
    }

    for(j=0;j<i;j++)
    {
      for(ii=0;ii<3;ii++) /*FIND MID NODE ALREADY CREATED ON FILM.*/
      {
        if(count[ii]==0)
        {
          for(jj=0;jj<3;jj++)
          {
            if(((bqf->films+i)->n12[soe[ii]]
                ==(bqf->films+j)->n12[soe[jj]] &&
                (bqf->films+i)->n12[soe[ii+1]]
                ==(bqf->films+j)->n12[soe[jj+1]]) ||
               ((bqf->films+i)->n12[soe[ii]]
                ==(bqf->films+j)->n12[soe[jj+1]] &&
                (bqf->films+i)->n12[soe[ii+1]]
                ==(bqf->films+j)->n12[soe[jj]]))
            {
              count[ii]=1;
			  (bqf->films+i)->n09[ii]=(bqf->films+j)->n09[jj];
			}
          }
        }
      }
    }

    for(ii=0;ii<3;ii++) /*CASE NO WIRE ON FILM LINE.*/
    {
      if(count[ii]==0)
      {
        /*bqf->nmids=(struct node09 *)
                   realloc(bqf->nmids,
                           (nnmid+1)*sizeof(struct node09));*/
        /*
        nmids=(struct node09 *)malloc((nnmid+1)*sizeof(struct node09));
        for(jj=0;jj<nnmid;jj++)
        {
          *(nmids+jj)=*(bqf->nmids+jj);
          for(kk=0;kk<(bqf->nwire);kk++)
          {
            if((bqf->wires+kk)->n09==bqf->nmids+jj)
            {
              (bqf->wires+kk)->n09=(nmids+jj);
            }
          }
          for(kk=0;kk<i;kk++)
          {
            for(ll=0;ll<3;ll++)
            {
              if((bqf->films+kk)->n09[ll]==(bqf->nmids+jj))
              {
                (bqf->films+kk)->n09[ll]=(nmids+jj);
              }
            }
          }
        }
        free(bqf->nmids);
        bqf->nmids=nmids;
        */

        mid.loff=nnmid;
        mid.code=nnmid+1; /*UNDER CONSTRUCTION.*/

        mid.c[GX]=0.5*((bqf->films+i)->n12[soe[ii]]->d[0]
                      +(bqf->films+i)->n12[soe[ii+1]]->d[0]);
        mid.c[GY]=0.5*((bqf->films+i)->n12[soe[ii]]->d[4]
                      +(bqf->films+i)->n12[soe[ii+1]]->d[4]);
        mid.c[GZ]=0.5*((bqf->films+i)->n12[soe[ii]]->d[8]
                      +(bqf->films+i)->n12[soe[ii+1]]->d[8]);

        *(bqf->nmids+nnmid)=mid;                 /*MIDNODE OF FILM.*/
        (bqf->films+i)->n09[soe[ii]]=bqf->nmids+nnmid;

        nnmid++;
      }
    }
  }
  if(bqf->nnmid!=nnmid) MessageBox(NULL,"Nnmid Incorrect.","Create",MB_OK);

  /*bqf->confs=(struct oconf *)
             realloc(bqf->confs,
                     (12*(bqf->nnode)+9*nnmid)*sizeof(struct oconf));*/
  /*
  confs=(struct oconf *)malloc((12*(bqf->nnode)+9*nnmid)*sizeof(struct oconf));
  for(i=0;i<12*(bqf->nnode);i++)
  {
    *(confs+i)=*(bqf->confs+i);
  }
  */
  for(i=0;i<nnmid;i++) /*CONFINEMENTS OF MIDNODES.*/
  {
    for(j=0;j<9;j++)
    {
      loff=12*(bqf->nnode)+9*i+j;

      (bqf->confs+loff)->iconf=0;
      (bqf->confs+loff)->value=0.0;
      /*(confs+loff)->iconf=0;*/
      /*(confs+loff)->value=0.0;*/
    }
  }
  /*free(bqf->confs);*/
  /*bqf->confs=confs;*/

  return;
}/*createmidnodes*/

void inputbquadinit(FILE *fin,int *nnode,
                              int *nwire,int *nfilm,
                              int *nwsec,int *nfsec)
/*INPUT INITIAL DATA.*/
{
  char **data;
  int n;

  data=fgetsbrk(fin,&n);
  *nnode=strtol(*(data+0),NULL,10);
  *nwire=strtol(*(data+1),NULL,10);
  *nfilm=strtol(*(data+2),NULL,10);

  *nwsec=strtol(*(data+3),NULL,10);
  *nfsec=strtol(*(data+4),NULL,10);

  freestr(data,n);

  return;
}/*inputbquadinit*/

void inputbquadnode12(FILE *fdisp,struct node12 *n12)
{
  int i;
  long int loffset;
  double d[12];

  loffset=(n12->loff)*12*sizeof(double);
  fseek(fdisp,loffset,SEEK_SET);
  fread(d,sizeof(double),12,fdisp);
  for(i=0;i<12;i++) n12->d[i]=d[i];

  return;
}/*inputbquadnode12*/

void inputbquadnode09(FILE *fdisp,int nnode,struct node09 *n09)
{
  int i;
  long int loffset;
  double d[9];

  loffset=(12*nnode+9*(n09->loff))*sizeof(double);
  fseek(fdisp,loffset,SEEK_SET);
  fread(d,sizeof(double),9,fdisp);
  for(i=0;i<9;i++) n09->d[i]=d[i];

  return;
}/*inputbquadnode09*/

void inputwire(struct qwire *wires,FILE *fwire,int offset,
			   struct qwire *wire)
{
  long int loffset;

  if(fwire!=NULL)
  {
	loffset=offset
			*(sizeof(struct qwire));
	fseek(fwire,loffset,SEEK_SET);
	fread(wire,sizeof(struct qwire),1,fwire);
  }
  else
  {
	*wire=*(wires+offset);
  }

  return;
}/*inputwire*/

void inputfilm(struct qfilm *films,FILE *ffilm,int offset,
			   struct qfilm *film)
{
  long int loffset;

  if(ffilm!=NULL)
  {
	loffset=offset
			*(sizeof(struct qfilm));
	fseek(ffilm,loffset,SEEK_SET);
	fread(film,sizeof(struct qfilm),1,ffilm);
  }
  else
  {
	*film=*(films+offset);
  }

  return;
}/*inputfilm*/

void initialbquadform(FILE *fdisp,struct node12 *nodes,int nnode,
                                  struct node09 *nmids,int nnmid)
/*INITIAL FORMATION INTO DISPLACEMENT FILE.WITHOUT NODE CODE.*/
{
  int i,ii,jj;
  double zero=0.0;

  fseek(fdisp,0L,SEEK_SET);
  for(i=0;i<nnode;i++)
  {
    for(ii=0;ii<3;ii++) /*U,V,W.*/
    {
      fwrite(&((nodes+i)->d[4*ii]),sizeof(double),1,fdisp);

      for(jj=0;jj<3;jj++) /*d/dX,d/dY,d/dZ.*/
      {fwrite(&zero,sizeof(double),1,fdisp);}
    }
  }
  for(i=0;i<nnmid;i++)
  {
    for(ii=0;ii<9;ii++) /*dU/dX,dU/dY,dU/dZ,.....,dW/dZ.*/
    {fwrite(&zero,sizeof(double),1,fdisp);}
  }
  fflush(fdisp);

  return;
}/*initialbquadform*/

void initialwire(struct qwire *wires,FILE *fwire,int nwire)
/*ASSEMBLAGE LONG STRESS OF WIRE.*/
{
  int i;

  for(i=0;i<nwire;i++)
  {
    fwrite((wires+i),sizeof(struct qwire),1,fwire);
  }
  fflush(fwire);

  return;
}/*initialwire*/

void initialfilm(struct qfilm *films,FILE *ffilm,int nfilm)
/*ASSEMBLAGE LONG STRESS OF FILM.*/
{
  int i;

  for(i=0;i<nfilm;i++)
  {
    fwrite((films+i),sizeof(struct qfilm),1,ffilm);
  }
  fflush(ffilm);

  return;
}/*initialfilm*/

void initialbquadreact(FILE *fin,FILE *freact,int nreact)
/*ASSEMBLAGE LONG REACTIONS.*/
{
  char **data;
  int i,n;
  double ddata;

  for(i=1;i<=nreact;i++)
  {
    data=fgetsbrk(fin,&n);
    if(data==NULL || n<=2) ddata=0.0;
    else                   ddata=strtod(*(data+2),NULL);
    fwrite(&ddata,sizeof(double),1,freact);

    for(;n>0;n--) free(*(data+n-1));
    free(data);
  }
  fflush(freact);

  return;
}/*initialbquadreact*/

double **bquadwiredrccos(double x1,double y1,double z1,
                         double x2,double y2,double z2,
                         double cangle)
/*RETURN:WIRE DIRECTION COSINE.*/
{
  int i;
  double dl,Xx,Yx,Zx,Xy,Yy,Zy,Xz,Yz,Zz;
  double **drccos,**rotate;

  drccos=(double **)malloc(3*sizeof(double *));
  for(i=0;i<=2;i++) *(drccos+i)=(double *)malloc(3*sizeof(double));

  Xz=x2-x1; /*z:WIRE DIRECTION.*/
  Yz=y2-y1;
  Zz=z2-z1;
  dl=sqrt(Xz*Xz+Yz*Yz+Zz*Zz);
  *(*(drccos+2)+0)=Xz/dl;
  *(*(drccos+2)+1)=Yz/dl;
  *(*(drccos+2)+2)=Zz/dl;

  Xx=Yz*1.0-Zz*0.0; /*x:OUTER PRODUCT z BY Z.*/
  Yx=Zz*0.0-Xz*1.0;
  if((Xx==0.0)&&(Yx==0.0)) Xx=1.0;
  Zx=Xz*0.0-Yz*0.0;
  dl=sqrt(Xx*Xx+Yx*Yx+Zx*Zx);
  *(*(drccos+0)+0)=Xx/dl;
  *(*(drccos+0)+1)=Yx/dl;
  *(*(drccos+0)+2)=Zx/dl;

  Xy=Yz*Zx-Zz*Yx; /*y:OUTER PRODUCT z BY x.*/
  Yy=Zz*Xx-Xz*Zx;
  Zy=Xz*Yx-Yz*Xx;
  dl=sqrt(Xy*Xy+Yy*Yy+Zy*Zy);
  *(*(drccos+1)+0)=Xy/dl;
  *(*(drccos+1)+1)=Yy/dl;
  *(*(drccos+1)+2)=Zy/dl;

  rotate=(double **)malloc(3*sizeof(double *));
  for(i=0;i<=2;i++) *(rotate+i)=(double *)malloc(3*sizeof(double));

  *(*(rotate+0)+0)=cos(cangle); /*ROTATION MATRIX AROUND Z.*/
  *(*(rotate+0)+1)=sin(cangle);
  *(*(rotate+0)+2)=0.0;
  *(*(rotate+1)+0)=-sin(cangle);
  *(*(rotate+1)+1)=cos(cangle);
  *(*(rotate+1)+2)=0.0;
  *(*(rotate+2)+0)=0.0;
  *(*(rotate+2)+1)=0.0;
  *(*(rotate+2)+2)=1.0;

  drccos=matrixmatrix(drccos,rotate,3);

  return drccos;
}/*bquadwiredrccos*/

double **bquadfilmdrccos(double x1,double y1,double z1,
                         double x2,double y2,double z2,
                         double x3,double y3,double z3)
/*RETURN:FILM DIRECTION COSINE.*/
{
  int i;
  double dl,Xx,Yx,Zx,Xy,Yy,Zy,Xz,Yz,Zz;
  double **drccos;

  drccos=(double **)malloc(3*sizeof(double *));
  for(i=0;i<=2;i++) *(drccos+i)=(double *)malloc(3*sizeof(double));

  Xx=x2-x1;
  Yx=y2-y1;
  Zx=z2-z1;
  dl=sqrt(Xx*Xx+Yx*Yx+Zx*Zx);
  *(*(drccos+0)+0)=Xx/dl;
  *(*(drccos+0)+1)=Yx/dl;
  *(*(drccos+0)+2)=Zx/dl;

  Xz=Yx*(z3-z1)-Zx*(y3-y1); /*OUTER PRODUCT.*/
  Yz=Zx*(x3-x1)-Xx*(z3-z1);
  Zz=Xx*(y3-y1)-Yx*(x3-x1);
  dl=sqrt(Xz*Xz+Yz*Yz+Zz*Zz);
  *(*(drccos+2)+0)=Xz/dl;
  *(*(drccos+2)+1)=Yz/dl;
  *(*(drccos+2)+2)=Zz/dl;

  Xy=Yz*Zx-Zz*Yx; /*OUTER PRODUCT.*/
  Yy=Zz*Xx-Xz*Zx;
  Zy=Xz*Yx-Yz*Xx;
  dl=sqrt(Xy*Xy+Yy*Yy+Zy*Zy);
  *(*(drccos+1)+0)=Xy/dl;
  *(*(drccos+1)+1)=Yy/dl;
  *(*(drccos+1)+2)=Zy/dl;

  return drccos;
}/*bquadfilmdrccos*/


double **assemwiremtx(struct qwire wire)
/*ASSEMBLAGE WIRE ELEMENT MATRIX.*/
{
  char str[256];
  int i,j,k;
  int utow[5]={ 0,12, 3,15,26};
  int vtow[5]={ 4,16, 7,19,29};
  int wtow[5]={ 8,20,11,23,32};
  int txtow[3]={2,14,25};
  int tytow[3]={5,17,27};
  double **wmtx,**C,**CC,**Bcoeff,*BB,**U,**V,**W,**Tx,**Ty;
  double integral[7];
  double dx,dy,dz,dl[3];
  double E,poi,A,Ixx,Iyy,Jx,Jy,G;

  dl[0]=0.0;
  dx=wire.n09->c[0]-wire.n12[0]->d[0];
  dy=wire.n09->c[1]-wire.n12[0]->d[4];
  dz=wire.n09->c[2]-wire.n12[0]->d[8];
  dl[1]=sqrt(dx*dx+dy*dy+dz*dz);
  dx=wire.n12[1]->d[0]-wire.n12[0]->d[0];
  dy=wire.n12[1]->d[4]-wire.n12[0]->d[4];
  dz=wire.n12[1]->d[8]-wire.n12[0]->d[8];
  dl[2]=sqrt(dx*dx+dy*dy+dz*dz);

  E=wire.sect->E;
  poi=wire.sect->poi;
  A=wire.sect->area;
  Ixx=wire.sect->Ixx;
  Iyy=wire.sect->Iyy;
  Jx=wire.sect->Ixx;
  Jy=wire.sect->Iyy;
  G=0.5*E/(1.0+poi);

  for(i=0;i<=6;i++) integral[i]=1.0/((i+1)*1.0)*pow(dl[1],i*1.0);

  wmtx=(double **)malloc(33*sizeof(double *));
  for(i=0;i<33;i++)
  {
    *(wmtx+i)=(double *)malloc(33*sizeof(double));
    for(j=0;j<33;j++) *(*(wmtx+i)+j)=0.0;
  }

  Bcoeff=(double **)malloc(5*sizeof(double *));
  for(i=0;i<5;i++)
  {
    *(Bcoeff+i)=(double *)malloc(4*sizeof(double));
    if(*(Bcoeff+i)==NULL) errormessage("Overflow.");
  }

  /*FOR U,V,W.*/
  C=(double **)malloc(5*sizeof(double *));
  for(i=0;i<5;i++)
  {
    *(C+i)=(double *)malloc(5*sizeof(double));
    if(*(C+i)==NULL) errormessage("Overflow.");
  }
  for(i=0;i<=1;i++) /*ASSEMBLAGE C MATRIX.*/
  {
    *(*(C+i)+0)=1.0;
    for(j=1;j<5;j++) *(*(C+i)+j)=pow(dl[i],j*1.0);
  }
  for(i=2;i<5;i++)
  {
    *(*(C+i)+0)=0.0;
    *(*(C+i)+1)=1.0;
    *(*(C+i)+2)=2.0*dl[i-2];
    *(*(C+i)+3)=3.0*dl[i-2]*dl[i-2];
    *(*(C+i)+4)=4.0*dl[i-2]*dl[i-2]*dl[i-2];
  }
  CC=fullmatrixcroutlu(C,5);
  if(CC==NULL) errormessage("Crout LU Failure.");

  /*FOR W.*/
  for(i=0;i<5;i++) /*ASSEMBLAGE B MATRIX.*/
  {
    *(*(Bcoeff+i)+0)=(*(*(CC+1)+i));
    *(*(Bcoeff+i)+1)=(*(*(CC+2)+i))*2.0;
    *(*(Bcoeff+i)+2)=(*(*(CC+3)+i))*3.0;
    *(*(Bcoeff+i)+3)=(*(*(CC+4)+i))*4.0;
  }

  W=(double **)malloc(5*sizeof(double *));
  for(i=0;i<=4;i++)
  {
    *(W+i)=(double *)malloc(5*sizeof(double));
    for(j=0;j<=4;j++) *(*(W+i)+j)=0.0;
  }
  for(i=0;i<=4;i++) /*ASSEMBLAGE W STIFFNESS MATRIX.*/
  {
    for(j=0;j<=4;j++)
    {
      BB=polypoly1(3,*(Bcoeff+i),*(Bcoeff+j));
      for(k=0;k<=6;k++) *(*(W+i)+j)+=E*A*(*(BB+k))*integral[k];
    }
  }
  for(i=0;i<=4;i++) /*W MATRIX INTO WIRE MATRIX.*/
  {
    for(j=0;j<=4;j++)
    {*(*(wmtx+wtow[i])+wtow[j])=*(*(W+i)+j);}
  }

  /*FOR U,V.*/
  for(i=0;i<5;i++) /*ASSEMBLAGE B MATRIX.*/
  {
    *(*(Bcoeff+i)+0)=*(*(CC+2)+i)*2.0;
    *(*(Bcoeff+i)+1)=*(*(CC+3)+i)*6.0;
    *(*(Bcoeff+i)+2)=*(*(CC+4)+i)*12.0;
  }

  U=(double **)malloc(5*sizeof(double *));
  V=(double **)malloc(5*sizeof(double *));
  for(i=0;i<=4;i++)
  {
    *(U+i)=(double *)malloc(5*sizeof(double));
    *(V+i)=(double *)malloc(5*sizeof(double));
    for(j=0;j<=4;j++)
    {
      *(*(U+i)+j)=0.0;
      *(*(V+i)+j)=0.0;
    }
  }
  for(i=0;i<=4;i++) /*ASSEMBLAGE U,V STIFFNESS MATRIX.*/
  {
    for(j=0;j<=4;j++)
    {
      BB=polypoly1(2,*(Bcoeff+i),*(Bcoeff+j));
      for(k=0;k<=4;k++)
      {
        *(*(U+i)+j)+=E*Ixx*(*(BB+k))*integral[k];
        *(*(V+i)+j)+=E*Iyy*(*(BB+k))*integral[k];
      }
    }
  }
  for(i=0;i<=4;i++) /*U,V MATRIX INTO WIRE MATRIX.*/
  {
    for(j=0;j<=4;j++)
    {
      *(*(wmtx+utow[i])+utow[j])=*(*(U+i)+j);
      *(*(wmtx+vtow[i])+vtow[j])=*(*(V+i)+j);
    }
  }

  for(i=0;i<5;i++) free(*(C+i));
  free(C);
  for(i=0;i<5;i++) free(*(CC+i));
  free(CC);

  /*FOR THETA Z.*/
  C=(double **)malloc(3*sizeof(double *));
  for(i=1;i<=3;i++) *(C+i-1)=(double *)malloc(3*sizeof(double));
  for(i=0;i<=2;i++) /*ASSEMBLAGE C MATRIX.*/
  {
    *(*(C+i)+0)=1.0;
    for(j=1;j<=2;j++) *(*(C+i)+j)=pow(dl[i],j*1.0);
  }
  CC=fullmatrixcroutlu(C,3);
  if(CC==NULL) errormessage("Crout LU Failure.");

  for(i=0;i<3;i++) /*ASSEMBLAGE B MATRIX.*/
  {
    *(*(Bcoeff+i)+0)=*(*(CC+1)+i);
    *(*(Bcoeff+i)+1)=*(*(CC+2)+i)*2.0;
  }

  Tx=(double **)malloc(3*sizeof(double *));
  Ty=(double **)malloc(3*sizeof(double *));
  for(i=0;i<=2;i++)
  {
    *(Tx+i)=(double *)malloc(3*sizeof(double));
    *(Ty+i)=(double *)malloc(3*sizeof(double));
    for(j=0;j<=2;j++)
    {
      *(*(Tx+i)+j)=0.0;
      *(*(Ty+i)+j)=0.0;
    }
  }
  for(i=0;i<=2;i++) /*ASSEMBLAGE THETA Z STIFFNESS MATRIX.*/
  {
    for(j=0;j<=2;j++)
    {
      BB=polypoly1(1,*(Bcoeff+i),*(Bcoeff+j));
      for(k=0;k<=2;k++)
      {
        *(*(Tx+i)+j)+=G*Jx*(*(BB+k))*integral[k];
        *(*(Ty+i)+j)+=G*Jy*(*(BB+k))*integral[k];
      }
    }
  }
  for(i=0;i<=2;i++) /*THETA Z MATRIX INTO WIRE MATRIX.*/
  {
    for(j=0;j<=2;j++)
    {
      *(*(wmtx+txtow[i])+txtow[j])=*(*(Tx+i)+j);
      *(*(wmtx+tytow[i])+tytow[j])=*(*(Ty+i)+j);
    }
  }

  for(i=0;i<3;i++) free(*(C+i));
  free(C);
  for(i=0;i<3;i++) free(*(CC+i));
  free(CC);

  for(i=0;i<5;i++) free(*(Bcoeff+i));
  free(Bcoeff);

  return wmtx;
}/*assemwiremtx*/

double **assemfilmmtx(struct qfilm film)
/*ASSEMBLAGE FILM MATRIX.*/
/*DEFORMATION BIQUADRATIC:INNER 30 DOF,OUTER 15 DOF.*/
{
  /*char string[300],str[80];*/
  int i,j,k,l,m,n,ii,jj,nn;
  int itof[30]={ 0, 1, 2, 4, 5, 6,  /*U1,V1*/
                12,13,14,16,17,18,  /*U2,V2*/
                24,25,26,28,29,30,  /*U3,V3*/
                   36,37,   39,40,  /*U4,V4*/
                   45,46,   48,49,  /*U5,V5*/
                   54,55,   57,58}; /*U6,V6*/
  int otof[15]={ 8, 9,10,  /*W1*/
                20,21,22,  /*W2*/
                32,33,34,  /*W3*/
                   42,43,  /*W4*/
                   51,52,  /*W5*/
                   60,61}; /*W6*/
  int ttof[18]={ 3, 7,11,  /*dU1/dZ,dV1/dZ,dW1/dZ*/
                15,19,23,  /*dU2/dZ,dV2/dZ,dW2/dZ*/
                27,31,35,  /*dU3/dZ,dV3/dZ,dW3/dZ*/
                38,41,44,  /*dU4/dZ,dV4/dZ,dW4/dZ*/
                47,50,53,  /*dU5/dZ,dV5/dZ,dW5/dZ*/
                56,59,62}; /*dU6/dZ,dV6/dZ,dW6/dZ*/
  double **fmtx,**inner,**outer,**drccos,*coord,*projected;
  double **Cb,**Ci,*CC,***BIcoeff,*BOcoeff[3][15];
  double *BB,integral[28],delta,fact,x1,x2,x3,y1,y2,y3;
  double d1,d2,d3;
  double thick,E,poi,x[6],y[6],z[6],org[3];

  E=film.sect->E;
  poi=film.sect->poi;
  thick=film.sect->thick;

  drccos=bquadfilmdrccos(film.n12[0]->d[0],
                         film.n12[0]->d[4],
                         film.n12[0]->d[8],
                         film.n12[1]->d[0],
                         film.n12[1]->d[4],
                         film.n12[1]->d[8],
                         film.n12[2]->d[0],
                         film.n12[2]->d[4],
                         film.n12[2]->d[8]);

  /*ORIGIN POINT*/
  org[0]=film.n12[0]->d[0];
  org[1]=film.n12[0]->d[4];
  org[2]=film.n12[0]->d[8];

  coord=(double *)malloc(3*sizeof(double));
  for(i=0;i<3;i++)
  {
    *(coord+0)=film.n12[i]->d[0]-org[0];
    *(coord+1)=film.n12[i]->d[4]-org[1];
    *(coord+2)=film.n12[i]->d[8]-org[2];
    projected=matrixvector(drccos,coord,3);
    x[i]=*(projected+0);
    y[i]=*(projected+1);
    z[i]=*(projected+2);

    free(projected);
  }
  for(i=3;i<6;i++)
  {
    *(coord+0)=film.n09[i-3]->c[0]-org[0];
    *(coord+1)=film.n09[i-3]->c[1]-org[1];
    *(coord+2)=film.n09[i-3]->c[2]-org[2];
    projected=matrixvector(drccos,coord,3);
    x[i]=*(projected+0);
    y[i]=*(projected+1);
    z[i]=*(projected+2);

    free(projected);
  }
  free(coord);

  for(i=0;i<=2;i++) free(*(drccos+i));
  free(drccos);

for(ii=0;ii<6;ii++)
{
  fprintf(fout,"Node%d (%.1f,%.1f,%.1f)\n",
          ii+1,x[ii],y[ii],z[ii]);
}

  x1=x[0]; x2=x[1]; x3=x[2];
  y1=y[0]; y2=y[1]; y3=y[2];
  delta=(x2*y3-x3*y2)+(x3*y1-x1*y3)+(x1*y2-x2*y1);
  for(i=0;i<=27;i++) integral[i]=0.0;
  for(i=0;i<=6;i++) /*INTEGRAL OF POW(X,M)POW(Y,N).*/
  {
    for(j=0;j<=6;j++)
    {
      for(k=0;k<=6;k++)
      {
        for(l=0;l<=6;l++)
        {
          for(m=0;m<=6;m++)
          {
            for(n=0;n<=6;n++)
            {
              if(i+j+k+l+m+n<=6)
              {
                ii=i+j+k; jj=l+m+n; nn=0.5*(ii+jj)*(ii+jj+1)+jj;
                fact=delta
                    *factorial(i+l)*factorial(j+m)*factorial(k+n)
                    *factorial(i+j+k)*factorial(l+m+n)
                    /factorial(i)/factorial(j)/factorial(k)
                    /factorial(l)/factorial(m)/factorial(n)
                    /factorial(i+j+k+l+m+n+2);

                if(i!=0) fact*=pow(x1,i);
                if(j!=0) fact*=pow(x2,j);
                if(k!=0) fact*=pow(x3,k);
                if(l!=0) fact*=pow(y1,l);
                if(m!=0) fact*=pow(y2,m);
                if(n!=0) fact*=pow(y3,n);

                integral[nn]+=fact;
              }
            }
          }
        }
      }
    }
  }

  Cb=(double **)malloc(15*sizeof(double *));
  for(n=0;n<=2;n++) /*ASSEMBLAGE C MATRIX.*/
  {
    CC=(double *)malloc(15*sizeof(double));
    *(CC+0)=1.0;
    *(CC+1)=x[n];
    *(CC+2)=y[n];
    *(CC+3)=pow(x[n],2.0);
    *(CC+4)=x[n]*y[n];
    *(CC+5)=pow(y[n],2.0);
    *(CC+6)=pow(x[n],3.0);
    *(CC+7)=pow(x[n],2.0)*y[n];
    *(CC+8)=x[n]*pow(y[n],2.0);
    *(CC+9)=pow(y[n],3.0);
    *(CC+10)=pow(x[n],4.0);
    *(CC+11)=pow(x[n],3.0)*y[n];
    *(CC+12)=pow(x[n],2.0)*pow(y[n],2.0);
    *(CC+13)=x[n]*pow(y[n],3.0);
    *(CC+14)=pow(y[n],4.0);

    *(Cb+3*n)=CC;
  }
  for(n=0;n<=5;n++)
  {
    CC=(double *)malloc(15*sizeof(double));
    *(CC+0)=0.0;
    *(CC+1)=1.0;
    *(CC+2)=0.0;
    *(CC+3)=2.0*x[n];
    *(CC+4)=y[n];
    *(CC+5)=0.0;
    *(CC+6)=3.0*pow(x[n],2.0);
    *(CC+7)=2.0*x[n]*y[n];
    *(CC+8)=pow(y[n],2.0);
    *(CC+9)=0.0;
    *(CC+10)=4.0*pow(x[n],3.0);
    *(CC+11)=3.0*pow(x[n],2.0)*y[n];
    *(CC+12)=2.0*x[n]*pow(y[n],2.0);
    *(CC+13)=pow(y[n],3.0);
    *(CC+14)=0.0;

    *(Cb+3*n+1+(2-n)*(int)(n/3))=CC;

    CC=(double *)malloc(15*sizeof(double));
    *(CC+0)=0.0;
    *(CC+1)=0.0;
    *(CC+2)=1.0;
    *(CC+3)=0.0;
    *(CC+4)=x[n];
    *(CC+5)=2.0*y[n];
    *(CC+6)=0.0;
    *(CC+7)=pow(x[n],2.0);
    *(CC+8)=2.0*x[n]*y[n];
    *(CC+9)=3.0*pow(y[n],2.0);
    *(CC+10)=0.0;
    *(CC+11)=pow(x[n],3.0);
    *(CC+12)=2.0*pow(x[n],2.0)*y[n];
    *(CC+13)=3.0*x[n]*pow(y[n],2.0);
    *(CC+14)=4.0*pow(y[n],3.0);

    *(Cb+3*n+2+(2-n)*(int)(n/3))=CC;
  }

fprintf(fout,"[C]\n");
for(ii=0;ii<15;ii++)
{
  for(jj=0;jj<15;jj++)
  {
    fprintf(fout," %5.1f",*(*(Cb+ii)+jj));
  }
  fprintf(fout,"\n");
}

  Ci=matrixinverse(Cb,15); /*NOT SYMMETRIC.*/
  /*Ci=fullmatrixcroutlu(Cb,15);*/
  if(Ci==NULL)
  {
    errormessage("Crout LU Failure.");
    return NULL;
  }

MessageBox(NULL,"Inversed.","[C]",MB_OK);
fprintf(fout,"[C]Inverse\n");
for(ii=0;ii<15;ii++)
{
  for(jj=0;jj<15;jj++)
  {
    fprintf(fout," %8.1E",*(*(Ci+ii)+jj));
  }
  fprintf(fout,"\n");
}
fmtx=matrixmatrix(Cb,Ci,15);
fprintf(fout,"[Cb][Ci]\n");
for(ii=0;ii<15;ii++)
{
  for(jj=0;jj<15;jj++)
  {
    fprintf(fout," %5.1f",*(*(fmtx+ii)+jj));
  }
  fprintf(fout,"\n");
}

  for(i=0;i<15;i++) free(*(Cb+i));
  free(Cb);

  fmtx=(double **)malloc(63*sizeof(double *));
  if(fmtx==NULL)
  {
    errormessage("INNER MATRIX OVERFLOW.");
    return NULL;
  }
  for(i=0;i<63;i++)
  {
    *(fmtx+i)=(double *)malloc(63*sizeof(double));
    if(*(fmtx+i)==NULL)
    {
      errormessage("INNER MATRIX OVERFLOW.");
      return NULL;
    }
    for(j=0;j<63;j++) *(*(fmtx+i)+j)=0.0; /*INITIAL.*/
  }

  /*INNER MATRIX.*/
  d1=E/(1.0-poi*poi);
  d2=E*poi/(1.0-poi*poi);
  d3=E/2.0/(1.0+poi);

  BIcoeff=(double ***)malloc(3*sizeof(double **));
  if(BIcoeff==NULL)
  {
    errormessage("MATRIX [B] OVERFLOW.");
    return NULL;
  }
  for(i=0;i<3;i++)
  {
    *(BIcoeff+i)=(double **)malloc(30*sizeof(double *));
    if(*(BIcoeff+i)==NULL)
    {
      errormessage("MATRIX [Bi] OVERFLOW.");
      return NULL;
    }
    for(j=0;j<30;j++)
    {
      *(*(BIcoeff+i)+j)=(double *)malloc(10*sizeof(double));
      if(*(*(BIcoeff+i)+j)==NULL)
      {
        errormessage("MATRIX [Bij] OVERFLOW.");
        return NULL;
      }
    }
  }
  for(i=0;i<15;i++) /*ASSEMBLAGE B MATRIX.*/
  {
    *(*(*(BIcoeff+0)+i)+0)=*(*(Ci+1)+i);
    *(*(*(BIcoeff+0)+i)+1)=*(*(Ci+3)+i)*2.0;
    *(*(*(BIcoeff+0)+i)+2)=*(*(Ci+4)+i);
    *(*(*(BIcoeff+0)+i)+3)=*(*(Ci+6)+i)*3.0;
    *(*(*(BIcoeff+0)+i)+4)=*(*(Ci+7)+i)*2.0;
    *(*(*(BIcoeff+0)+i)+5)=*(*(Ci+8)+i);
    *(*(*(BIcoeff+0)+i)+6)=*(*(Ci+10)+i)*4.0;
    *(*(*(BIcoeff+0)+i)+7)=*(*(Ci+11)+i)*3.0;
    *(*(*(BIcoeff+0)+i)+8)=*(*(Ci+12)+i)*2.0;
    *(*(*(BIcoeff+0)+i)+9)=*(*(Ci+13)+i);

    *(*(*(BIcoeff+1)+15+i)+0)=*(*(Ci+2)+i);
    *(*(*(BIcoeff+1)+15+i)+1)=*(*(Ci+4)+i);
    *(*(*(BIcoeff+1)+15+i)+2)=*(*(Ci+5)+i)*2.0;
    *(*(*(BIcoeff+1)+15+i)+3)=*(*(Ci+7)+i);
    *(*(*(BIcoeff+1)+15+i)+4)=*(*(Ci+8)+i)*2.0;
    *(*(*(BIcoeff+1)+15+i)+5)=*(*(Ci+9)+i)*3.0;
    *(*(*(BIcoeff+1)+15+i)+6)=*(*(Ci+11)+i);
    *(*(*(BIcoeff+1)+15+i)+7)=*(*(Ci+12)+i)*2.0;
    *(*(*(BIcoeff+1)+15+i)+8)=*(*(Ci+13)+i)*3.0;
    *(*(*(BIcoeff+1)+15+i)+9)=*(*(Ci+14)+i)*4.0;

    for(j=0;j<10;j++)
    {
      *(*(*(BIcoeff+0)+15+i)+j)=0.0;
      *(*(*(BIcoeff+1)+i)+j)   =0.0;
      *(*(*(BIcoeff+2)+i)+j)   =*(*(*(BIcoeff+1)+15+i)+j);
      *(*(*(BIcoeff+2)+15+i)+j)=*(*(*(BIcoeff+0)+i)+j);
    }
  }

/*for(ii=0;ii<3;ii++)
{
  fprintf(fout,"[Binner %d]\n",ii+1);
  for(jj=0;jj<30;jj++)
  {
    sprintf(string,"\0");
    for(nn=0;nn<10;nn++)
    {
      sprintf(str," %12.3E",*(*(*(BIcoeff+ii)+jj)+nn));
      strcat(string,str);
    }
    fprintf(fout,"%s\n",string);
  }
}*/

  inner=(double **)malloc(30*sizeof(double *));
  for(i=0;i<30;i++)
  {
    *(inner+i)=(double *)malloc(30*sizeof(double));
    for(j=0;j<30;j++) *(*(inner+i)+j)=0.0;
  }
  for(i=0;i<30;i++) /*ASSEMBLAGE INNER STIFFNESS MATRIX.*/
  {
    for(j=0;j<30;j++)
    {
      BB=polypoly2(3,*(*(BIcoeff+0)+i),*(*(BIcoeff+0)+j));
      for(k=0;k<=27;k++)
      {*(*(inner+i)+j)+=thick*d1*(*(BB+k))*integral[k];}
      BB=polypoly2(3,*(*(BIcoeff+0)+i),*(*(BIcoeff+1)+j));
      for(k=0;k<=27;k++)
      {*(*(inner+i)+j)+=thick*d2*(*(BB+k))*integral[k];}
      BB=polypoly2(3,*(*(BIcoeff+1)+i),*(*(BIcoeff+0)+j));
      for(k=0;k<=27;k++)
      {*(*(inner+i)+j)+=thick*d2*(*(BB+k))*integral[k];}
      BB=polypoly2(3,*(*(BIcoeff+1)+i),*(*(BIcoeff+1)+j));
      for(k=0;k<=27;k++)
      {*(*(inner+i)+j)+=thick*d1*(*(BB+k))*integral[k];}
      BB=polypoly2(3,*(*(BIcoeff+2)+i),*(*(BIcoeff+2)+j));
      for(k=0;k<=27;k++)
      {*(*(inner+i)+j)+=thick*d3*(*(BB+k))*integral[k];}
    }
  }
  for(i=0;i<30;i++) /*INNER MATRIX INTO PLATE MATRIX.*/
  {
    for(j=0;j<30;j++)
    {*(*(fmtx+itof[i])+itof[j])=*(*(inner+i)+j);}
  }

/*fprintf(fout,"[Inner]\n");
for(ii=0;ii<30;ii++)
{
  sprintf(string,"\0");
  for(jj=0;jj<30;jj++)
  {
    sprintf(str," %1.0f",*(*(inner+ii)+jj)
                       /(*(*(inner+ii)+jj)));
    strcat(string,str);
  }
  fprintf(fout,"%s\n",string);
}*/

MessageBox(NULL,"LU Begin","Inner",MB_OK);
Cb=fullmatrixcroutlu(inner,30);
if(Cb==NULL) MessageBox(NULL,"Instable.","Inner",MB_OK);
if(Cb!=NULL) MessageBox(NULL,"Stable.",  "Inner",MB_OK);

  for(i=0;i<30;i++) free(*(inner+i));
  free(inner);

  /*OUTER MATRIX.*/
  d1=E*pow(thick,3.0)/12.0/(1.0-poi*poi);
  d2=E*poi*pow(thick,3.0)/12.0/(1.0-poi*poi);
  d3=E*pow(thick,3.0)/24.0/(1.0+poi);

  for(i=0;i<15;i++) /*ASSEMBLAGE B MATRIX.*/
  {
    BOcoeff[0][i]=(double *)malloc(6*sizeof(double));
    BOcoeff[1][i]=(double *)malloc(6*sizeof(double));
    BOcoeff[2][i]=(double *)malloc(6*sizeof(double));

    *(BOcoeff[0][i]+0)=*(*(Ci+3)+i)*2.0;
    *(BOcoeff[0][i]+1)=*(*(Ci+6)+i)*6.0;
    *(BOcoeff[0][i]+2)=*(*(Ci+7)+i)*2.0;
    *(BOcoeff[0][i]+3)=*(*(Ci+10)+i)*12.0;
    *(BOcoeff[0][i]+4)=*(*(Ci+11)+i)*6.0;
    *(BOcoeff[0][i]+5)=*(*(Ci+12)+i)*2.0;

    *(BOcoeff[1][i]+0)=*(*(Ci+5)+i)*2.0;
    *(BOcoeff[1][i]+1)=*(*(Ci+8)+i)*2.0;
    *(BOcoeff[1][i]+2)=*(*(Ci+9)+i)*6.0;
    *(BOcoeff[1][i]+3)=*(*(Ci+12)+i)*2.0;
    *(BOcoeff[1][i]+4)=*(*(Ci+13)+i)*6.0;
    *(BOcoeff[1][i]+5)=*(*(Ci+14)+i)*12.0;

    *(BOcoeff[2][i]+0)=*(*(Ci+4)+i)*2.0;
    *(BOcoeff[2][i]+1)=*(*(Ci+7)+i)*4.0;
    *(BOcoeff[2][i]+2)=*(*(Ci+8)+i)*4.0;
    *(BOcoeff[2][i]+3)=*(*(Ci+11)+i)*6.0;
    *(BOcoeff[2][i]+4)=*(*(Ci+12)+i)*8.0;
    *(BOcoeff[2][i]+5)=*(*(Ci+13)+i)*6.0;
  }

  outer=(double **)malloc(15*sizeof(double *));
  for(i=0;i<15;i++)
  {
    *(outer+i)=(double *)malloc(15*sizeof(double));
    for(j=0;j<15;j++) *(*(outer+i)+j)=0.0;
  }
  for(i=0;i<15;i++) /*ASSEMBLAGE OUTER STIFFNESS MATRIX.*/
  {
    for(j=0;j<15;j++)
    {
      BB=polypoly2(3,*(*(BIcoeff+0)+i),*(*(BIcoeff+0)+j));
      for(k=0;k<=14;k++)
      {*(*(outer+i)+j)+=d1*(*(BB+k))*integral[k];}
      BB=polypoly2(3,*(*(BIcoeff+0)+i),*(*(BIcoeff+1)+j));
      for(k=0;k<=14;k++)
      {*(*(outer+i)+j)+=d2*(*(BB+k))*integral[k];}
      BB=polypoly2(3,*(*(BIcoeff+1)+i),*(*(BIcoeff+0)+j));
      for(k=0;k<=14;k++)
      {*(*(outer+i)+j)+=d2*(*(BB+k))*integral[k];}
      BB=polypoly2(3,*(*(BIcoeff+1)+i),*(*(BIcoeff+1)+j));
      for(k=0;k<=14;k++)
      {*(*(outer+i)+j)+=d1*(*(BB+k))*integral[k];}
      BB=polypoly2(3,*(*(BIcoeff+2)+i),*(*(BIcoeff+2)+j));
      for(k=0;k<=14;k++)
      {*(*(outer+i)+j)+=d3*(*(BB+k))*integral[k];}
    }
  }
  for(i=0;i<=14;i++) /*OUTER MATRIX INTO PLATE MATRIX.*/
  {
    for(j=0;j<=14;j++)
    {*(*(fmtx+otof[i])+otof[j])=*(*(outer+i)+j);}
  }

  for(i=0;i<15;i++) free(*(outer+i));
  free(outer);
  for(i=0;i<15;i++) free(*(Ci+i));
  free(Ci);

  for(i=0;i<18;i++) /*LINES OUT OF USE.UNDER CONSTRUCTION.*/
  {
    *(*(fmtx+ttof[i])+ttof[i])=1.0;
  }

  for(i=0;i<3;i++)
  {
    for(j=0;j<30;j++) free(*(*(BIcoeff+i)+j));
    free(*(BIcoeff+i));
  }
  free(BIcoeff);

  return fmtx;
}/*assemfilmmtx*/

double **assemtrifilmmtx1(struct qfilm film,int mmode)
/*ASSEMBLAGE TRIANGLE LINEAR FILM MATRIX.*/
/*DEFORMATION LINEAR:INNER 2 DOF.*/
{
  char string[300],str[80];
  int i,j,k,l,m,n,ii,jj;
  int itof[6]={  0,12,24,  /*U*/
                 4,16,28}; /*V*/
  double **fmtx,**inner,**drccos,*coord,*projected;
  double **Cb,**Cbt,**Ci,*CC,BIcoeff[3][6][1],**DD,**BD;
  double *BB,delta,fact,x1,x2,x3,y1,y2,y3;
  double d1,d2,d3;
  double thick,E,poi,x[6],y[6],z[6],org[3];
  double **DB;

  E=film.sect->E;
  poi=film.sect->poi;
  thick=film.sect->thick;

  drccos=bquadfilmdrccos(film.n12[0]->d[0],
                         film.n12[0]->d[4],
                         film.n12[0]->d[8],
                         film.n12[1]->d[0],
                         film.n12[1]->d[4],
                         film.n12[1]->d[8],
                         film.n12[2]->d[0],
                         film.n12[2]->d[4],
                         film.n12[2]->d[8]);

  /*ORIGIN POINT*/
  org[0]=film.n12[0]->d[0];
  org[1]=film.n12[0]->d[4];
  org[2]=film.n12[0]->d[8];

  coord=(double *)malloc(3*sizeof(double));
  for(i=0;i<3;i++)
  {
    *(coord+0)=film.n12[i]->d[0]-org[0];
    *(coord+1)=film.n12[i]->d[4]-org[1];
    *(coord+2)=film.n12[i]->d[8]-org[2];
    projected=matrixvector(drccos,coord,3);
    x[i]=*(projected+0);
    y[i]=*(projected+1);
    z[i]=*(projected+2);

    free(projected);
  }
  for(i=3;i<6;i++)
  {
    *(coord+0)=film.n09[i-3]->c[0]-org[0];
    *(coord+1)=film.n09[i-3]->c[1]-org[1];
    *(coord+2)=film.n09[i-3]->c[2]-org[2];
    projected=matrixvector(drccos,coord,3);
    x[i]=*(projected+0);
    y[i]=*(projected+1);
    z[i]=*(projected+2);

    free(projected);
  }
  free(coord);

  for(i=0;i<=2;i++) free(*(drccos+i));
  free(drccos);

/*
if(fout!=NULL)
{
  for(ii=0;ii<6;ii++)
  {
    fprintf(fout,"Node%d (%.1f,%.1f,%.1f)\n",
            ii+1,x[ii],y[ii],z[ii]);
  }
}
*/

  x1=x[0]; x2=x[1]; x3=x[2];
  y1=y[0]; y2=y[1]; y3=y[2];
  delta=(x2*y3-x3*y2)+(x3*y1-x1*y3)+(x1*y2-x2*y1);

  Cb=(double **)malloc(3*sizeof(double *));
  for(n=0;n<3;n++) /*ASSEMBLAGE C MATRIX.*/
  {
    CC=(double *)malloc(3*sizeof(double));

    *(CC+0)=1.0;
    *(CC+1)=x[n];
    *(CC+2)=y[n];

    *(Cb+n)=CC;
  }

/*
if(fout!=NULL)
{
  fprintf(fout,"[Cb]\n");
  for(ii=0;ii<3;ii++)
  {
    for(jj=0;jj<3;jj++)
    {
      fprintf(fout," %5.1f",*(*(Cb+ii)+jj));
    }
    fprintf(fout,"\n");
  }
}
*/

  Ci=mallocdoublematrix(3);
  Cbt=matrixtranspose(Cb,3); /*MATRIX [Ct].*/
  Ci=matrixinverse(Cbt,3); /*NOT SYMMETRIC.*/
  /*matrixinverseII(Ci,Cbt,3);*/ /*NOT SYMMETRIC.*/
  if(Ci==NULL)
  {
    errormessage("Inverse Failure.");
    return NULL;
  }

/*MessageBox(NULL,"Inversed.","[Ci]",MB_OK);*/
/*INVERSE CHECK*/
/*
if(fout!=NULL)
{
  fprintf(fout,"[Ci]=Inversed [Cbt]\n");
  Cbt=matrixtranspose(Cb,3);
  for(ii=0;ii<3;ii++)
  {
    for(jj=0;jj<3;jj++)
    {
      fprintf(fout," %8.1E",*(*(Ci+ii)+jj));
    }
    fprintf(fout,"    ");
    for(jj=0;jj<3;jj++)
    {
      fprintf(fout," %8.1E",*(*(Cbt+ii)+jj));
    }
    fprintf(fout,"\n");
  }
  fmtx=matrixmatrix(Ci,Cbt,3);
  fprintf(fout,"CHECK [Ci][Cbt]\n");
  for(ii=0;ii<3;ii++)
  {
    for(jj=0;jj<3;jj++)
    {
      fprintf(fout," %5.1f",*(*(fmtx+ii)+jj));
    }
    fprintf(fout,"\n");
  }
}
*/

  for(i=0;i<3;i++) free(*(Cb+i));
  free(Cb);
  for(i=0;i<3;i++) free(*(Cbt+i));
  free(Cbt);

  for(i=0;i<3;i++) /*ASSEMBLAGE B MATRIX.*/
  {
    BIcoeff[0][i][0]=*(*(Ci+i)+1);
    BIcoeff[1][3+i][0]=*(*(Ci+i)+2);

    BIcoeff[0][3+i][0]=0.0;
    BIcoeff[1][i][0]  =0.0;
    BIcoeff[2][i][0]  =BIcoeff[1][3+i][0];
    BIcoeff[2][3+i][0]=BIcoeff[0][i][0];
  }

/*for(ii=0;ii<3;ii++)
{
  fprintf(fout,"[Binner %d]\n",ii+1);
  for(jj=0;jj<6;jj++)
  {
    sprintf(string,"\0");
    sprintf(str," %12.3E",BIcoeff[ii][jj][0]);
    strcat(string,str);
    fprintf(fout,"%s\n",string);
  }
}*/
/*
if(fout!=NULL)
{
  fprintf(fout,"[Binner]\n");
  for(ii=0;ii<3;ii++)
  {
    for(jj=0;jj<6;jj++)
    {
      fprintf(fout," %12.3E",BIcoeff[ii][jj][0]);
    }
    fprintf(fout,"\n");
  }
}
*/

  fmtx=(double **)malloc(63*sizeof(double *));
  if(fmtx==NULL)
  {
    errormessage("INNER MATRIX OVERFLOW.");
    return NULL;
  }
  for(i=0;i<63;i++)
  {
    *(fmtx+i)=(double *)malloc(63*sizeof(double));
    if(*(fmtx+i)==NULL)
    {
      errormessage("INNER MATRIX OVERFLOW.");
      return NULL;
    }
    for(j=0;j<63;j++) /*INITIAL.*/
    {
      *(*(fmtx+i)+j)=0.0;
    }
  }

  /*INNER MATRIX.*/
  DD=(double **)malloc(3*sizeof(double *));
  for(i=0;i<3;i++)
  {
    *(DD+i)=(double *)malloc(3*sizeof(double));
    for(j=0;j<3;j++) *(*(DD+i)+j)=0.0;
  }

  d1=E/(1.0-poi*poi);
  d2=E*poi/(1.0-poi*poi);
  d3=E/2.0/(1.0+poi);

  /*
  d1=E*(1.0-poi)/(1.0+poi)/(1.0-2.0*poi);
  d2=E*poi/(1.0+poi)/(1.0-2.0*poi);
  d3=E/2.0/(1.0+poi);
  */

  *(*(DD+0)+0)=d1;
  *(*(DD+0)+1)=d2;
  *(*(DD+0)+2)=0.0;
  *(*(DD+1)+0)=d2;
  *(*(DD+1)+1)=d1;
  *(*(DD+1)+2)=0.0;
  *(*(DD+2)+0)=0.0;
  *(*(DD+2)+1)=0.0;
  *(*(DD+2)+2)=d3;

/*
if(fout!=NULL)
{
  fprintf(fout,"[D]\n");
  for(ii=0;ii<3;ii++)
  {
    for(jj=0;jj<3;jj++)
    {
      fprintf(fout," %12.3E",*(*(DD+ii)+jj));
    }
    fprintf(fout,"\n");
  }
}
*/

if(mmode==SMATRIX) /*RETURN [D][B]*/
{
  DB=(double **)malloc(3*sizeof(double *));
  for(i=0;i<3;i++)
  {
    *(DB+i)=(double *)malloc(6*sizeof(double));
    for(j=0;j<6;j++) *(*(DB+i)+j)=0.0;
  }
  for(i=0;i<3;i++)
  {
    for(j=0;j<6;j++)
    {
      for(k=0;k<3;k++)
      {
        *(*(DB+i)+j)+=(*(*(DD+i)+k))*BIcoeff[k][j][0];
      }
    }
  }
  for(i=0;i<3;i++) free(*(DD+i));
  free(DD);

  return DB;
}

  BD=(double **)malloc(6*sizeof(double *));
  for(i=0;i<6;i++)
  {
    *(BD+i)=(double *)malloc(3*sizeof(double));
    for(j=0;j<3;j++) *(*(BD+i)+j)=0.0;
  }
  inner=(double **)malloc(6*sizeof(double *));
  for(i=0;i<6;i++)
  {
    *(inner+i)=(double *)malloc(6*sizeof(double));
    for(j=0;j<6;j++) *(*(inner+i)+j)=0.0;
  }

  for(i=0;i<6;i++)
  {
    for(j=0;j<3;j++)
    {
      for(k=0;k<3;k++)
      {
        *(*(BD+i)+j)+=BIcoeff[k][i][0]*(*(*(DD+k)+j));
      }
    }
  }

/*
if(fout!=NULL)
{
  fprintf(fout,"[Bt][D]\n");
  for(ii=0;ii<6;ii++)
  {
    for(jj=0;jj<3;jj++)
    {
      fprintf(fout," %12.3E",*(*(BD+ii)+jj));
    }
    fprintf(fout,"\n");
  }
}
*/

  for(i=0;i<6;i++)
  {
    for(j=0;j<6;j++)
    {
      for(k=0;k<3;k++)
      {
        *(*(inner+i)+j)+=thick*delta
                         *(*(*(BD+i)+k))*BIcoeff[k][j][0];
      }
    }
  }

/*
if(fout!=NULL)
{
  fprintf(fout,"[Inner]\n");
  for(ii=0;ii<6;ii++)
  {
    sprintf(string,"\0");
    for(jj=0;jj<6;jj++)
    {
      sprintf(str," %11.3E",*(*(inner+ii)+jj));
      strcat(string,str);
    }
    fprintf(fout,"%s\n",string);
  }
}
*/

/*MessageBox(NULL,"Stability Check","Inner",MB_OK);*/
/*Cb=matrixinverse(inner,6);*/
/*
Cb=fullmatrixcroutlu(inner,6);
if(Cb==NULL) MessageBox(NULL,"Instable.","Inner",MB_OK);
if(Cb!=NULL) MessageBox(NULL,"Stable.",  "Inner",MB_OK);
*/
  for(i=0;i<6;i++) /*INNER MATRIX INTO PLATE MATRIX.*/
  {
    for(j=0;j<6;j++)
    {*(*(fmtx+itof[i])+itof[j])=*(*(inner+i)+j);}
  }

  for(i=0;i<3;i++) free(*(DD+i));
  free(DD);
  for(i=0;i<3;i++) free(*(BD+i));
  free(BD);
  for(i=0;i<6;i++) free(*(inner+i));
  free(inner);
  for(i=0;i<3;i++) free(*(Ci+i));
  free(Ci);

  return fmtx;
}/*assemtrifilmmtx1*/

double **transmatrix12(double **drccos)
/*ASSEMBLAGE TRANSFORMATION MATRIX FOR 12 DIMENSIONS.*/
{
  int i,j,ii,jj;
  double **t;

  t=(double **)malloc(12*sizeof(double *));
  for(i=0;i<12;i++)
  {
    *(t+i)=(double *)malloc(12*sizeof(double));
    for(j=0;j<12;j++) *(*(t+i)+j)=0.0; /*INITIAL.*/
  }

  for(i=0;i<3;i++) /*FOR U V W.*/
  {
    for(j=0;j<3;j++)
    {*(*(t+4*i)+4*j)=*(*(drccos+i)+j);}

    for(j=0;j<3;j++) /*FOR d/dX,d/dY,d/dZ.*/
    {
      for(ii=0;ii<3;ii++)
      {
        for(jj=0;jj<3;jj++)
        {
          *(*(t+4*i+1+ii)+4*j+1+jj)=(*(*(drccos+i)+j))
                                   *(*(*(drccos+ii)+jj));
        }
      }
    }
  }

  return t;
}/*transmatrix12*/

double **transmatrix09(double **drccos)
/*ASSEMBLAGE TRANSFORMATION MATRIX FOR 9 DIMENSIONS.*/
{
  int i,j,ii,jj;
  double **t;

  t=(double **)malloc(9*sizeof(double *));
  for(i=0;i<9;i++)
  {
    *(t+i)=(double *)malloc(9*sizeof(double));
    for(j=0;j<9;j++) *(*(t+i)+j)=0.0; /*INITIAL.*/
  }

  for(i=0;i<3;i++) /*FOR U V W.*/
  {
    for(j=0;j<3;j++)
    {
      for(ii=0;ii<3;ii++)
      {
        for(jj=0;jj<3;jj++)
        {
          *(*(t+3*i+ii)+3*j+jj)=(*(*(drccos+i)+j))
                               *(*(*(drccos+ii)+jj));
        }
      }
    }
  }

  return t;
}/*transmatrix09*/

double **wiretransmatrix(double **drccos)
/*ASSEMBLAGE TRANSFORMATION MATRIX FOR WIRE.*/
{
  int i,j,n;
  double **t,**t12,**t09;

  t=(double **)malloc(33*sizeof(double *));
  for(i=0;i<33;i++)
  {
    *(t+i)=(double *)malloc(33*sizeof(double));
    for(j=0;j<33;j++) *(*(t+i)+j)=0.0; /*INITIAL.*/
  }

  t12=transmatrix12(drccos);
  t09=transmatrix09(drccos);

  for(n=0;n<2;n++) /*FOR NODE 1,2.*/
  {
    for(i=0;i<12;i++)
    {
      for(j=0;j<12;j++)
      {
        *(*(t+12*n+i)+12*n+j)=*(*(t12+i)+j);
      }
    }
  }
  for(i=0;i<9;i++) /*FOR MID NODE.*/
  {
    for(j=0;j<9;j++)
    {
      *(*(t+24+i)+24+j)=*(*(t09+i)+j);
    }
  }

  for(i=0;i<12;i++) free(*(t12+i));
  free(t12);
  for(i=0;i<9;i++) free(*(t09+i));
  free(t09);

  return t;
}/*wiretransmatrix*/

double **filmtransmatrix(double **drccos)
/*ASSEMBLAGE TRANSFORMATION MATRIX FOR FILM.*/
{
  int i,j,n;
  double **t,**t12,**t09;

  t=(double **)malloc(63*sizeof(double *));
  for(i=0;i<63;i++)
  {
    *(t+i)=(double *)malloc(63*sizeof(double));
    for(j=0;j<63;j++) *(*(t+i)+j)=0.0; /*INITIAL.*/
  }

  t12=transmatrix12(drccos);
  t09=transmatrix09(drccos);

  for(n=0;n<3;n++) /*FOR NODE 1,2,3.*/
  {
    for(i=0;i<12;i++)
    {
      for(j=0;j<12;j++)
      {
        *(*(t+12*n+i)+12*n+j)=*(*(t12+i)+j);
      }
    }
  }
  for(n=0;n<3;n++) /*FOR MID NODE 1,2,3.*/
  {
    for(i=0;i<9;i++)
    {
      for(j=0;j<9;j++)
      {
        *(*(t+36+9*n+i)+36+9*n+j)=*(*(t09+i)+j);
      }
    }
  }

  for(i=0;i<12;i++) free(*(t12+i));
  free(t12);
  for(i=0;i<9;i++) free(*(t09+i));
  free(t09);

  return t;
}/*filmtransmatrix*/

double **bquadtransformation(double **estiff,double **tmatrix,
                             int msize)
/*ELEMENT MATRIX TRANSFORMATION.*/
{
  double **ematrix;

  ematrix=matrixmatrix(estiff,tmatrix,msize);
  tmatrix=matrixtranspose(tmatrix,msize);
  estiff=matrixmatrix(tmatrix,ematrix,msize);
  return estiff;
}/*bquadtransformation*/

double **modifywirehinge(struct qwire wire,double **wstiff)
/*MODIFY WIRE ELEMENT MATRIX BY HINGE.*/
{
  int n,i,ii,jj,kk;
  double h[33][33],**w=wstiff;

  for(n=0;n<3;n++)
  {
    for(i=0;i<9;i++)
    {
      if(wire.iconf[n][i]==1)
      {
        kk=12*n+3*(1-(int)(n/2))+i; /*HINGE LINE.*/
        for(ii=0;ii<33;ii++)
        {
          for(jj=0;jj<33;jj++)
          {
            h[ii][jj]=- *(*(w+ii)+kk)
                      / *(*(w+kk)+kk)
                      * *(*(w+kk)+jj);
          }
        }
        for(ii=0;ii<33;ii++)
        {
          for(jj=0;jj<33;jj++) *(*(w+ii)+jj)+=h[ii][jj];
        }
      }
    }
  }
  return w;
}/*modifywirehinge*/

double **modifyfilmhinge(struct qfilm film,double **fstiff)
/*MODIFY FILM ELEMENT MATRIX BY HINGE.*/
{
  int n,i,ii,jj,kk;
  double h[63][63],**f=fstiff;

  for(n=0;n<6;n++)
  {
    for(i=0;i<9;i++)
    {
      if(film.iconf[n][i]==1)
      {
        kk=12*n+3-(3*n-6)*(int)(n/3)+i;
        for(ii=0;ii<63;ii++)
        {
          for(jj=0;jj<63;jj++)
          {
            h[ii][jj]=- *(*(f+ii)+kk)
                      / *(*(f+kk)+kk)
                      * *(*(f+kk)+jj);
          }
        }
        for(ii=0;ii<63;ii++)
        {
          for(jj=0;jj<63;jj++) *(*(f+ii)+jj)+=h[ii][jj];
        }
      }
    }
  }
  return f;
}/*modifyfilmhinge*/

void assemwtog(struct gcomponent *gmtx,
               double **wstiff,struct qwire *wire,int nnode)
/*ASSEMBLAGE WIRE MATRIX INTO GLOBAL MATRIX.*/
{
  long int loff[3];
  long int i,j,ii,jj,n1,n2;
  double wdata,gdata;

  for(i=0;i<2;i++) loff[i]=wire->n12[i]->loff;
  loff[2]=wire->n09->loff;

  for(n1=0;n1<2;n1++) /*END NODES.*/
  {
    for(i=0;i<12;i++)
    {
      ii=12*loff[n1]+i;
      for(n2=0;n2<2;n2++)
      {
        for(j=0;j<12;j++)
        {
          jj=12*loff[n2]+j;
          if(ii>=jj)
          {
            wdata=*(*(wstiff+12*n1+i)+12*n2+j);
            if(wdata!=0.0)
            {
              gread(gmtx,ii+1,jj+1,&gdata);
              gdata+=wdata;
              gwrite(gmtx,ii+1,jj+1,gdata);
            }
          }
        }
      }
      for(j=0;j<9;j++)
      {
        jj=12*nnode+9*(loff[2]-nnode)+j;
        if(ii>=jj)
        {
          wdata=*(*(wstiff+12*n1+i)+24+j);
          if(wdata!=0.0)
          {
            gread(gmtx,ii+1,jj+1,&gdata);
            gdata+=wdata;
            gwrite(gmtx,ii+1,jj+1,gdata);
          }
        }
      }
    }
  }
  for(i=0;i<9;i++) /*MID NODE.*/
  {
    ii=12*nnode+9*(loff[2]-nnode)+i;
    for(n2=0;n2<2;n2++)
    {
      for(j=0;j<12;j++)
      {
        jj=12*loff[n2]+j;
        if(ii>=jj)
        {
          wdata=*(*(wstiff+24+i)+12*n2+j);
          if(wdata!=0.0)
          {
            gread(gmtx,ii+1,jj+1,&gdata);
            gdata+=wdata;
            gwrite(gmtx,ii+1,jj+1,gdata);
          }
        }
      }
    }
    for(j=0;j<9;j++)
    {
      jj=12*nnode+9*(loff[2]-nnode)+j;
      if(ii>=jj)
      {
        wdata=*(*(wstiff+24+i)+24+j);
        if(wdata!=0.0)
        {
          gread(gmtx,ii+1,jj+1,&gdata);
          gdata+=wdata;
          gwrite(gmtx,ii+1,jj+1,gdata);
        }
      }
    }
  }
  return;
}/*assemwtog*/

void assemftog(struct gcomponent *gmtx,
               double **fstiff,struct qfilm *film,int nnode)
/*ASSEMBLAGE FILM MATRIX INTO GLOBAL MATRIX.*/
{
  long int loff[6];
  long int i,j,ii,jj,n1,n2;
  double fdata,gdata;

  for(i=0;i<3;i++) loff[i]=film->n12[i]->loff;
  for(i=0;i<3;i++) loff[i+3]=film->n09[i]->loff;

  for(n1=0;n1<3;n1++) /*END NODES.*/
  {
    for(i=0;i<12;i++)
    {
      ii=12*loff[n1]+i;
      for(n2=0;n2<3;n2++)
      {
        for(j=0;j<12;j++)
        {
          jj=12*loff[n2]+j;
          if(ii>=jj)
          {
            fdata=*(*(fstiff+12*n1+i)+12*n2+j);
            if(fdata!=0.0)
            {
              gread(gmtx,ii+1,jj+1,&gdata);
              gdata+=fdata;
              gwrite(gmtx,ii+1,jj+1,gdata);
            }
          }
        }
      }
      for(n2=0;n2<3;n2++)
      {
        for(j=0;j<9;j++)
        {
          jj=12*nnode+9*loff[n2+3]+j;
          if(ii>=jj)
          {
            fdata=*(*(fstiff+12*n1+i)+36+9*n2+j);
            if(fdata!=0.0)
            {
              gread(gmtx,ii+1,jj+1,&gdata);
              gdata+=fdata;
              gwrite(gmtx,ii+1,jj+1,gdata);
            }
          }
        }
      }
    }
  }
  for(n1=0;n1<3;n1++) /*MID NODES.*/
  {
    for(i=0;i<9;i++)
    {
      ii=12*nnode+9*loff[n1+3]+i;
      for(n2=0;n2<3;n2++)
      {
        for(j=0;j<12;j++)
        {
          jj=12*loff[n2]+j;
          if(ii>=jj)
          {
            fdata=*(*(fstiff+36+9*n1+i)+12*n2+j);
            if(fdata!=0.0)
            {
              gread(gmtx,ii+1,jj+1,&gdata);
              gdata+=fdata;
              gwrite(gmtx,ii+1,jj+1,gdata);
            }
          }
        }
      }
      for(n2=0;n2<3;n2++)
      {
        for(j=0;j<9;j++)
        {
          jj=12*nnode+9*loff[n2+3]+j;
          if(ii>=jj)
          {
            fdata=*(*(fstiff+36+9*n1+i)+36+9*n2+j);
            if(fdata!=0.0)
            {
              gread(gmtx,ii+1,jj+1,&gdata);
              gdata+=fdata;
              gwrite(gmtx,ii+1,jj+1,gdata);
            }
          }
        }
      }
    }
  }
  return;
}/*assemftog*/

double *extractwirevct(struct qwire wire,double *gvct,int nnode)
/*EXTRACT WIRE DISPLACEMENT FROM GLOBAL VECTOR.*/
{
  long int i,k,loff;
  int n;
  double *d;

  d=(double *)malloc(33*sizeof(double));
  if(d==NULL)
  {
    errormessage("MEMORY INSUFFICIENT.");
    return NULL;
  }
  k=0;
  for(n=0;n<2;n++)
  {
    for(i=0;i<12;i++)
    {
      loff=12*(wire.n12[n]->loff)+i;
      *(d+k)=*(gvct+loff);
      k++;
    }
  }
  for(i=0;i<9;i++)
  {
    loff=12*nnode+9*((wire.n09->loff)-nnode)+i;
    *(d+k)=*(gvct+loff);
    k++;
  }
  return d;
}/*extractwirevct*/

double *extractfilmvct(struct qfilm film,double *gvct,int nnode)
/*EXTRACT FILM DISPLACEMENT FROM GLOBAL VECTOR.*/
{
  long int i,k,loff;
  int n;
  double *d;

  d=(double *)malloc(63*sizeof(double));
  if(d==NULL)
  {
    errormessage("MEMORY INSUFFICIENT.");
    return NULL;
  }
  k=0;
  for(n=0;n<3;n++)
  {
    for(i=0;i<12;i++)
    {
      loff=12*(film.n12[n]->loff)+i;
      *(d+k)=*(gvct+loff);
      k++;
    }
  }
  for(n=0;n<3;n++)
  {
    for(i=0;i<9;i++)
    {
      loff=12*nnode+9*((film.n09[n]->loff)-nnode)+i;
      *(d+k)=*(gvct+loff);
      k++;
    }
  }
  return d;
}/*extractfilmvct*/

double *wirestress(struct qwire *wire,double *gvct,int nnode)
/*WIRE ELEMENT STRESS.*/
{
  int ii;
  double **drccos,**tmatrix,**wstiff,*wstress;

  drccos=bquadwiredrccos(wire->n12[0]->d[0],
                         wire->n12[0]->d[4],
                         wire->n12[0]->d[8],
                         wire->n12[1]->d[0],
                         wire->n12[1]->d[4],
                         wire->n12[1]->d[8],
                         wire->cangle);         /*DIRECTION COSINE.*/

  tmatrix=wiretransmatrix(drccos);         /*TRANSFORMATION MATRIX.*/
  wstiff=assemwiremtx(*wire);                /*ELASTIC WIRE MATRIX.*/
  wstiff=modifywirehinge(*wire,wstiff);

  wstress=extractwirevct(*wire,gvct,nnode);
  wstress=matrixvector(tmatrix,wstress,33);
  wstress=matrixvector(wstiff,wstress,33);

  for(ii=0;ii<=2;ii++) free(*(drccos+ii));
  free(drccos);
  for(ii=0;ii<=33;ii++)
  {
    free(*(tmatrix+ii));
    free(*(wstiff+ii));
  }
  free(tmatrix);
  free(wstiff);

  return wstress;
}/*wirestress*/

double *filmstress(struct qfilm *film,double *gvct,int nnode)
/*FILM ELEMENT STRESS.*/
{
  int ii;
  double **drccos,**tmatrix,**fstiff,*fstress;

  drccos=bquadfilmdrccos(film->n12[0]->d[0],
                         film->n12[0]->d[4],
                         film->n12[0]->d[8],
                         film->n12[1]->d[0],
                         film->n12[1]->d[4],
                         film->n12[1]->d[8],
                         film->n12[2]->d[0],
                         film->n12[2]->d[4],
                         film->n12[2]->d[8]);            /*[DRCCOS]*/

  tmatrix=filmtransmatrix(drccos);         /*TRANSFORMATION MATRIX.*/
  /*fstiff=assemfilmmtx(*film);*/                /*ELASTIC FILM MATRIX.*/
fstiff=assemtrifilmmtx1(*film,FMATRIX); /*TRIANGLE LINEAR FILM MATRIX.*/
  fstiff=modifyfilmhinge(*film,fstiff);

  fstress=extractfilmvct(*film,gvct,nnode);
  fstress=matrixvector(tmatrix,fstress,63);
  fstress=matrixvector(fstiff,fstress,63);

  for(ii=0;ii<=2;ii++) free(*(drccos+ii));
  free(drccos);
  /*
  for(ii=0;ii<=63;ii++)
  {
    free(*(tmatrix+ii));
    free(*(fstiff+ii));
  }
  free(tmatrix);
  free(fstiff);
  */
  
  return fstress;
}/*filmstress*/

void bquadassemconf(struct oconf *confs,double *gvct,double dsafety,
                    int nnode,int nnmid)
/*ASSEMBLAGE CONFINEMENT VALUE INTO GLOBAL VECTOR.*/
{
  int i;
  double conf;

  for(i=0;i<12*nnode;i++)     /*LOADS OR DISPS GIVEN INCREMENTALLY.*/
  {
    conf=(confs+i)->value;

    /**(gvct+i-1)+=dsafety*conf;*/   /*"+=":IF GVECTOR INITIALIZED.*/
    *(gvct+i)=dsafety*conf;         /*"=":IF GVECTOR UNINITIALIZED.*/
  }
  for(i=12*nnode;i<(12*nnode+9*nnmid);i++)         /*FOR MIDPOINTS.*/
  {
    *(gvct+i)=0.0;
  }

  return;
}/*bquadassemconf*/

void bquadmodifygivend(struct gcomponent *gmtx,double *gvct,
                       struct oconf *confs,int nnode,int nnmid)
/*MODIFY GLOBAL VECTOR BY GIVEN DISPLACEMENT.*/
{
  signed char iconf;
  long int ii,jj,msize;
  double gstiff,disp;

  msize=12*nnode+9*nnmid;

  for(ii=1;ii<=msize;ii++)
  {
    disp=*(gvct+ii-1);
    iconf=(confs+ii-1)->iconf;

    if(iconf==1 && disp!=0.0)
    {
      for(jj=1;jj<=msize;jj++)
	  {
        iconf=(confs+jj-1)->iconf;

        if(iconf!=1)
        {
          gread(gmtx,jj,ii,&gstiff);
          if(gstiff!=0.0)
          {
            *(gvct+jj-1)-=gstiff*disp;
          }
        }
      }
    }
  }
  return;
}/*bquadmodifygivend*/

void bquadoutputdisp(double *gvct,FILE *fout,int nnode,int nnmid,
                     struct node12 *nodes,struct node09 *nmids)
/*OUTPUT NODE DISPLACEMENT.*/
{
  int i,j;

  if(fout==NULL) return;

  for(i=0;i<nnode;i++) /*END NODES.*/
  {
    fprintf(fout,"NODE:%5ld {dU}=",(nodes+i)->code);
    for(j=0;j<12;j++)
    {
      fprintf(fout," %9.5f",*(gvct+12*i+j));
    }
    fprintf(fout,"\n");
  }
}/*bquadoutputdisp*/

void bquadupdateform(FILE *fdisp,double *gvct,int nnode,int nnmid)
/*void bquadupdateform(double *ddisp,double *gvct,int nnode,int nnmid)*/
/*FORMATION UPDATE.*/
{
  int i,j;
  long int loff,msize;
  double *ddisp;

  msize=12*nnode+9*nnmid;
  ddisp=(double *)malloc(msize*sizeof(double));

  fseek(fdisp,0L,SEEK_SET);
  fread(ddisp,sizeof(double),msize,fdisp);

  for(i=0;i<msize;i++) *(ddisp+i)+=*(gvct+i); /*{U}+{dU}*/

  fseek(fdisp,0L,SEEK_SET);
  fwrite(ddisp,sizeof(double),msize,fdisp);

  return;
}/*bquadupdateform*/

void outputwirestress(struct qwire wire,double *wstress,FILE *fout)
/*OUTPUT WIRE ELEMENT STRESS.*/
{
  int i,n;

  for(n=0;n<2;n++)
  {
    fprintf(fout,"WIRE:%5ld NODE:%5ld",
            wire.code,(wire.n12[n]->code));
    for(i=0;i<12;i++) fprintf(fout," %5.2E",*(wstress+12*n+i));
    fprintf(fout,"\n");
  }
  fprintf(fout,"WIRE:%5ld NMID:%5ld",(wire.n09->code));
  for(i=0;i<9;i++) fprintf(fout," %5.2E",*(wstress+24+i));
  fprintf(fout,"\n");

  return;
}/*outputwirestress*/

void outputfilmstress(struct qfilm film,double *fstress,FILE *fout)
/*OUTPUT FILM ELEMENT STRESS.*/
{
  int i,n;

  for(n=0;n<3;n++) /*END NODES.*/
  {
	fprintf(fout,"FILM:%5ld NODE:%5ld",
            film.code,(film.n12[n]->code));
    for(i=0;i<12;i++) fprintf(fout," %5.2E",*(fstress+12*n+i));
    fprintf(fout,"\n");
  }
  for(n=0;n<3;n++) /*MID NODES.*/
  {
    fprintf(fout,"FILM:%5ld NMID:%5ld",
            film.code,(film.n09[n]->code));
    for(i=0;i<9;i++) fprintf(fout," %5.2E",*(fstress+36+i));
    fprintf(fout,"\n");
  }
  return;
}/*outputfilmstress*/

void bquadoutputreaction(struct gcomponent *gmtx,
                         double *gvct,
                         struct node12 *nodes,
                         struct node09 *nmids,
                         struct oconf *confs,
                         FILE *freact,FILE *fout,
                         int nnode,int nnmid)
/*REACTIONS UPDATE,OUTPUT.*/
{
  char string[256];
  char iconf;
  int offset;
  long int i,j,nreact=0;
  double gstiff,reaction,dreaction;

  for(i=1;i<=12*nnode;i++)
  {
    iconf=(confs+i-1)->iconf;
    offset=(i-1)/12;

    if(iconf==1)
    {
      nreact++;
      dreaction=0.0;
      for(j=1;j<=(12*nnode+9*nnmid);j++)
      {
        gread(gmtx,i,j,&gstiff);
		dreaction+=gstiff*(*(gvct+j-1));             /*{dR}=[K]{dU}*/
      }
      vread(freact,nreact,&reaction);            /*UPDATE REACTION.*/
	  reaction+=dreaction;                               /*{R}+{dR}*/
	  vwrite(freact,nreact,&reaction);

	  sprintf(string,"NODE:%5ld %ld %14.5f",
			  (nodes+offset)->code,(i-1)%12+1,dreaction);
	  if(fout!=NULL) fprintf(fout,"%s\n",string);
	}
  }
  return;
}/*bquadoutputreaction*/

void drawbiquadframe(HDC hdc,struct viewparam vp,
					 struct biquadframe bf,
					 int mode)
{
  HPEN hpen,ppen;
  int i;
  long int loff;
  struct onode on[3];

  drawglobalaxis(hdc,vp,0,0,255);               /*DRAW GLOBAL AXIS.*/

  for(i=0;i<bf.nwire;i++) /*WIRES.*/
  {
	drawbiquadwire(hdc,vp,*(bf.wires+i),mode);
  }
  for(i=0;i<bf.nfilm;i++) /*FILMS.*/
  {
	/*DRAW INITIAL FILM*/
	if(mode==ONPRINTER) hpen=CreatePen(PS_SOLID,1,RGB(0,0,0));
	else hpen=CreatePen(PS_SOLID,1,RGB(150,150,150));
	ppen=(HPEN)SelectObject(hdc,hpen);

	loff=((bf.films+i)->n12[0])->loff;
	on[0].d[0]=(bf.ninit+loff)->d[0];
	on[0].d[1]=(bf.ninit+loff)->d[4];
	on[0].d[2]=(bf.ninit+loff)->d[8];
	loff=((bf.films+i)->n12[1])->loff;
	on[1].d[0]=(bf.ninit+loff)->d[0];
	on[1].d[1]=(bf.ninit+loff)->d[4];
	on[1].d[2]=(bf.ninit+loff)->d[8];
	loff=((bf.films+i)->n12[2])->loff;
	on[2].d[0]=(bf.ninit+loff)->d[0];
	on[2].d[1]=(bf.ninit+loff)->d[4];
	on[2].d[2]=(bf.ninit+loff)->d[8];
	drawgloballine(hdc,vp,on[0],on[1]);
	drawgloballine(hdc,vp,on[1],on[2]);
	drawgloballine(hdc,vp,on[2],on[0]);

	SelectObject(hdc,ppen);
	DeleteObject(hpen);

	/*DRAW DEFORMED FILM*/
	if(vp.vflag.ev.deformation)
	{
/*
	  for(ii=0;ii<3;ii++)
	  {
		for(jj=0;jj<12;jj++)
		{
		  film.ndisp[12*ii+jj]=*(gvct+12*(film.n12[ii]->loff)+jj);
		}
	  }
	  for(ii=0;ii<3;ii++)
	  {
		for(jj=0;jj<9;jj++)
		{
		  film.ndisp[36+9*ii+jj]=*(gvct+12*nnode+9*(film.n09[ii]->loff)+jj);
		}
	  }
*/
	  drawbiquadfilm(hdc,vp,*(bf.films+i),mode);
	}
  }

  SetTextColor(hdc,RGB(150,150,255));
  for(i=0;i<bf.nnode;i++) /*NODES.*/
  {
	drawbiquadnode12(hdc,vp,*(bf.ninit+i));
	/*drawbiquadnode12(hdc,vp,*(bf.nodes+i));*/
  }
  SetTextColor(hdc,RGB(255,150,255));
  if(bf.nmids!=NULL)
  {
	for(i=0;i<bf.nnmid;i++) /*MID NODES.*/
	{
	  drawbiquadnode09(hdc,vp,*(bf.nmids+i));
	}
  }

  return;
}/*drawbiquadframe*/

void drawbiquadnodes(HDC hdc,struct viewparam vp,
					 struct biquadframe bf)
{
  int i;

  drawglobalaxis(hdc,vp,0,0,255);               /*DRAW GLOBAL AXIS.*/

  SetTextColor(hdc,RGB(150,150,255));
  /*NODES.*/
  for(i=0;i<bf.nnode;i++)
  {
	drawbiquadnode12(hdc,vp,*(bf.ninit+i));
  }

 /*MID NODES.*/
/*
  SetTextColor(hdc,RGB(255,150,255));
  if(bf.nmids!=NULL)
  {
	for(i=0;i<bf.nnmid;i++)
	{
	  drawbiquadnode09(hdc,vp,*(bf.nmids+i));
	}
  }
*/
  return;
}/*drawbiquadnodes*/

void drawbiquadnode12(HDC hdc,struct viewparam vp,struct node12 bn)
/*DRAW NODE GLOBAL ON SCREEN BY AXONOMETRICS,PERSPECTIVE.*/
{
  char str[20];
  struct onode gn;

  sprintf(str,"%ld",bn.code);

  gn.d[0]=bn.d[0];
  gn.d[1]=bn.d[4];
  gn.d[2]=bn.d[8];

  drawglobaltext(hdc,vp,gn,str);

  return;
}/*drawbiquadnode12*/

void drawbiquadnode09(HDC hdc,struct viewparam vp,struct node09 bn)
/*DRAW NODE GLOBAL ON SCREEN BY AXONOMETRICS,PERSPECTIVE.*/
{
  char str[20];
  struct onode gn;

  sprintf(str,"%ld",bn.code);

  gn.d[0]=bn.c[0];
  gn.d[1]=bn.c[1];
  gn.d[2]=bn.c[2];

  drawglobaltext(hdc,vp,gn,str);

  return;
}/*drawbiquadnode09*/


void drawbiquadwire(HDC hdc,struct viewparam vp,
					struct qwire qw,
					int mode)
/*DRAW WIRE ON SCREEN BY AXONOMETRICS,PERSPECTIVE.*/
{
  HPEN hpen,ppen;
  struct onode on1,on2;

  if(mode==ONPRINTER) hpen=CreatePen(PS_SOLID,1,RGB(0,0,0));
  else hpen=CreatePen(PS_SOLID,1,RGB(255,255,150));

  ppen=(HPEN)SelectObject(hdc,hpen);

  on1.d[0]=qw.n12[0]->d[0];
  on1.d[1]=qw.n12[0]->d[4];
  on1.d[2]=qw.n12[0]->d[8];
  on2.d[0]=qw.n12[1]->d[0];
  on2.d[1]=qw.n12[1]->d[4];
  on2.d[2]=qw.n12[1]->d[8];
  drawgloballine(hdc,vp,on1,on2);

  SelectObject(hdc,ppen);
  DeleteObject(hpen);

  return;
}/*drawbiquadwire*/

void drawbiquadfilm(HDC hdc,struct viewparam vp,
					struct qfilm qf,
					int mode)
/*DRAW FILM ON SCREEN BY AXONOMETRICS,PERSPECTIVE.*/
{
  HPEN hpen,ppen;
  DWORD rop; /*RASTER OPERATION.*/
  char str[256],s[64];
  struct onode on[3];
  struct obans ob={0,0,3,NULL};
  int i,j,sr,sg,sb;
  int itof[6]={0,12,24,  /*U*/
			   4,16,28}; /*V*/
  double **drccos,**tmatrix,*uvct,**DB,*svct;
  double smax,smin,value,rate;

  if(mode==ONPRINTER) hpen=CreatePen(PS_SOLID,1,RGB(0,0,0));
  else hpen=CreatePen(PS_SOLID,1,RGB(255,150,150));
  ppen=(HPEN)SelectObject(hdc,hpen);

  if(mode==ONPRINTER) rop=SRCAND;
  else                rop=SRCPAINT;

  on[0].d[0]=qf.n12[0]->d[0];
  on[0].d[1]=qf.n12[0]->d[4];
  on[0].d[2]=qf.n12[0]->d[8];
  on[1].d[0]=qf.n12[1]->d[0];
  on[1].d[1]=qf.n12[1]->d[4];
  on[1].d[2]=qf.n12[1]->d[8];
  on[2].d[0]=qf.n12[2]->d[0];
  on[2].d[1]=qf.n12[2]->d[4];
  on[2].d[2]=qf.n12[2]->d[8];

  ob.nods=(struct onode **)malloc(3*sizeof(struct onode *));
  *(ob.nods+0)=&(on[0]);
  *(ob.nods+1)=&(on[1]);
  *(ob.nods+2)=&(on[2]);

/*DRAW STRESS*/
uvct=(double *)malloc(63*sizeof(double));
for(i=0;i<63;i++) *(uvct+i)=qf.ndisp[i];
drccos=bquadfilmdrccos(qf.n12[0]->d[0],
					   qf.n12[0]->d[4],
					   qf.n12[0]->d[8],
					   qf.n12[1]->d[0],
					   qf.n12[1]->d[4],
					   qf.n12[1]->d[8],
					   qf.n12[2]->d[0],
					   qf.n12[2]->d[4],
					   qf.n12[2]->d[8]); /*[DRCCOS]*/
tmatrix=filmtransmatrix(drccos);         /*TRANSFORMATION MATRIX.*/
uvct=matrixvector(tmatrix,uvct,63);
DB=assemtrifilmmtx1(qf,SMATRIX);
smax= 500.000;
smin=-500.000;
svct=(double *)malloc(3*sizeof(double));
for(i=0;i<3;i++)
{
  *(svct+i)=0.0;
  for(j=0;j<6;j++)
  {
	*(svct+i)+=*(*(DB+i)+j)*(*(uvct+itof[j]));
  }
}

if(fout!=NULL)
{
  fprintf(fout,"{s}=[DB][u]\n");
  for(i=0;i<3;i++)
  {
	sprintf(str,"{ %9.3f}",*(svct+i));
	for(j=0;j<6;j++)
	{
	  sprintf(s," %12.3f",*(*(DB+i)+j));
	  strcat(str,s);
	}
	sprintf(s," { %9.3f}",*(uvct+itof[i]));
	strcat(str,s);
	fprintf(fout,"%s\n",str);
  }
  for(i=3;i<6;i++)
  {
	for(j=0;j<7;j++) fprintf(fout,"             ");
	fprintf(fout,"{ %9.3f}\n",*(uvct+itof[i]));
  }
}

value=0.5*(*(svct+0)+*(svct+1))
	  +sqrt(0.25*(*(svct+0)-*(svct+1))*(*(svct+0)-*(svct+1))
			+(*(svct+2))*(*(svct+2))); /*MAXIMUM PRINCIPLE STRESS.*/
rate=(value-smin)/(smax-smin);

/*
sprintf(str,"RATE=%.5f",rate);
MessageBox(NULL,str,"DRAWBIQUADFILM",MB_OK);
*/

sr=(int)(255.0*rate);
sg=(int)0;
sb=(int)(255.0*(1.0-rate*0.5));
SelectObject(hdc,ppen);
DeleteObject(hpen);
hpen=CreatePen(PS_SOLID,1,RGB(sr,sg,sb));
ppen=(HPEN)SelectObject(hdc,hpen);

sprintf(str,"FILM %d STRESS=%.3f RGB={%3d %3d %3d}",
		qf.code,value,sr,sg,sb);
if(fout!=NULL) fprintf(fout,"%s\n",str);
/*MessageBox(NULL,str,"DRAW",MB_OK);*/

  drawglobalban(hdc,vp,ob,sr,sg,sb,rop);

  free(ob.nods);
  free(uvct);
  free(svct);
  for(i=0;i<3;i++) free(*(drccos+i));
  free(drccos);
  for(i=0;i<63;i++) free(*(tmatrix+i));
  free(tmatrix);

  SelectObject(hdc,ppen);
  DeleteObject(hpen);
  return;
}








