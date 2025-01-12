/*SRCAL004.H:SRC SUBROUTINES SINCE 1995.12.21.JUNSATO.*/
/*HEADER FOR HOGAN001.C.*/
/*BASED ON SRCAM005.H.*/
/*LAST CHANGE:2000.8.24.*/

/*DRAW S,RC,SRC,PC SECTIONS.*/
/*SECTION FROM SECTION LIST FILE.*/

/*
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <malloc.h>
#include <math.h>
#include <limits.h>
*/

#define S   1                                        /*SECTION TYPE*/
#define RC  2
#define SRC 3
#define PC  4

#define PLONG     0 /*PERIOD*/
#define PSHORT    1
#define PULTIMATE 2

#define SX 0 /*AXIS OF SECTION*/
#define SY 1

#define HSTRONG 0 /*AXIS OF H*/
#define HWEAK   1

#define HEAD 0 /*END OF ELEMENT*/
#define TAIL 1

#define END 0 /*POSITION OF ELEMENT FOR PC*/
#define MID 1

#define UPPER 0 /*TENSED SIDE FOR BENDING OF PC GIRDER*/
#define LOWER 1

#define STRNDA 0 /*PC BAR TYPE:A,B,C,WIRE STRAND*/
#define STRNDB 1
#define STRNDC 2
#define STRNDW 3

#define MAXSECT 300 /*MAX OF SECTIONS*/
#define MAXCMQ  100 /*MAX OF CMQ LIST*/
#define MAXCOMMENT 100 /*MAX OF COMMENTS PER ONE SECTION*/

#define MAXSRECT    5                     /*MAX OF STEEL RECTANGLES*/
#define MAXREINS   50                       /*MAX OF REINFORCEMENTS*/
#define MAXCRECT    3                  /*MAX OF CONCRETE RECTANGLES*/
#define MAXSTRND   50                           /*MAX OF PC STRANDS*/

struct materials{
                 double sE,sF,sft,sfc,sfb,sfs,sftu,sfcu;        /*S*/
                 double rE,rft,rfc,wft,rftu,rfcu,wfp;           /*R*/
                 double cE,Fc,cfc,cfs,srcfc,srcfs,cfcu;         /*C*/
                 double pcE,pcF,pcfc,pcft,pcfs,pcfcu,pcfsu;    /*PC*/
                 double pcfco,pcfto;                  /*PC FOR DEAD*/
                 double stfact;    /*EFFECTIVE FACTOR FOR PC STRAND*/
                 double stE,stft[4],stftu[4];              /*STRAND*/
                };
struct materialrect{double top,bottom,left,right;}; /*RECTANGLE[mm]*/
struct reinforcement{double area,x,y;};         /*REIN AREA,X,Y[mm]*/
struct pcstrand{
                int type;                  /*STRAND TYPE:A,B,C,WIRE*/
                double area[2],x[2],y[2];       /*REIN AREA,X,Y[mm]*/
                double Ni[2];                 /*INITIAL TENSION[tf]*/
               };

struct section{
               long int code;                     /*CODE OF SECTION*/
               int soff;                        /*OFFSET OF SECTION*/
               int stype,etype;         /*SECTION TYPE,ELEMENT TYPE*/
               int nsteel,nrein,nconc;     /*STEELS,REINS,CONCRETES*/
               int nstrnd;                             /*PC STRANDS*/
               struct materialrect srect[MAXSRECT];   /*STEEL RECTS*/
               struct reinforcement rein[MAXREINS];         /*REINS*/
               struct materialrect crect[MAXCRECT]; /*CONCRETE RECT*/
               struct pcstrand strnd[MAXSTRND];        /*PC STRANDS*/
               double shearrein[2];        /*RATE OF REIN FOR SHEAR*/
               double thick,wlength,wheight,windowrate;
               double face[2][2];          /*FACE LENGTH[AXIS][END]*/
               double safetya;       /*MAX RATE OF SAFETY ALLOWABLE*/
               double safetyu;        /*MAX RATE OF SAFETY ULTIMATE*/

               int ncomment;
               char comment[MAXCOMMENT][80];
              };

struct structnode{
                  long int code;
                  double x,y,z;
                 };

struct stress{double N,Q[2],Mt,M[2];};     /*STRESS ON 6 DIRECTION.*/
struct estress{struct stress x,y,z;};  /*x,y:HORIZONTAL z:VERTICAL.*/
struct element{
               long int code;                        /*ELEMENT CODE*/
               struct section *sect;                      /*SECTION*/
               struct structnode *node[2];             /*NODES[END]*/
               struct estress head,tail;                 /*STRESSES*/
               int cmqcode;
               double Mo;                                     /*CMQ*/
              };

/*SUBROUTINES*/
int srcan001(char fname[]);
int createyieldsurface(struct arclmframe *af);

/*SUBROUTINES FOR STRING.*/
char **fgetscut(FILE *fin,int *n);

/*SUBROUTINES FROM SRCAN.*/
int getsectionform(FILE *flist,int code,struct section *sect);

/*SUBROUTINES FOR DRAWING.*/
void drawsectionlist(HDC hdc,int nsect,struct section *slist,
                     long int ox,long int oy,
                     long int pagewidth,long int pageheight,
                     struct viewparam vp,int mode);
void getsectionspace(HDC hdc,
                     struct section *sect,
                     double *Xmin,double *Xmax,
                     double *Ymin,double *Ymax,
                     int *commentwidth,int *commentheight);
void drawsection(HDC hdc,struct section *sect,
                 int sr,int sg,int sb,
                 int rr,int rg,int rb,
                 int cr,int cg,int cb,
                 int pr,int pg,int pb,
                 int tr,int tg,int tb,
                 int lw,
                 long int ox,long int oy,
                 struct viewparam vp);
void drawmaterialrect(HDC hdc,
                      struct materialrect *mr,
                      long int ox,long int oy,
                      struct viewparam vp);
void drawreincircle(HDC hdc,
                    struct reinforcement *r,
                    long int ox,long int oy,
                    struct viewparam vp);
void drawstrandcircle(HDC hdc,
                      struct pcstrand *ps,int position,
                      long int ox,long int oy,
                      struct viewparam vp);

/*SUBROUTINES FROM SRCAM*/
void comments(FILE *f);
void inputfiletomemory(FILE *ftext,
                       int *nnode,int *nelem,int *nsect,
                       struct section *sects,
                       struct structnode *nodes,
                       struct element *elems);
void inputfiletomemory2(FILE *ftext,
                       int *nnode,int *nelem,int *nsect,
                       struct section *sects,
                       struct structnode *nodes,
                       struct element *elems);
void fgetinitial(FILE *fin,int *nnode,int *nelem,int *nsect,
                 double *E,double *G);
void fgetinitial2(FILE *fin,int *nnode,int *nelem,int *nsect);
void initializestress(struct stress *s);
int fgetstress(FILE *fin,int *ielem,int *isect,int *inode,
               struct stress *s);
int getsrcansection(FILE *flist,long int code,struct section *sect);
int getcodelist(FILE *flist,long int codelist[]);
int getcmqlist(FILE *fcmq,int cmqcode[],double Mo[]);

double elementlength(struct element elem);
double walllength(struct element elem);
double wallheight(struct element elem);

void translatesection(struct section sect,
                      struct materials m,
                      double *As,double *Ac,double *Ar,
                      double *Yg,
                      int axis);
double allowablebendingofrc(struct element elem,
                            struct materials m,
                            int axis,
                            double Nd);
double allowablebendingofsrc(struct element elem,
                             struct materials m,
                             int axis,
                             double Nd);
double ultimatebendingofsrc(struct element elem,
                            struct materials m,
                            int axis,
                            double Nd,
                            double *sN,double *sM,
                            double *rN,double *rM,
                            double *cN,double *cM);

double boundaryvalue(struct element elem,struct materials m,
                     double Yn,
                     double sBi[],double sDi[],
                     double cBi[],double cDi[],
                     double sYi[],double sYj[],
                     double cYi[],double cYj[],
                     double sfi[],double sfj[],
                     double cfi[],double cfj[],
                     double rAi[],double rYi[],double rfi[],
                     double *gN,double *gM);
double boundaryvalueultimate(struct element elem,struct materials m,
                             double Yn,
                             double sBi[],double sDi[],
                             double cBi[],double cDi[],
                             double sYi[],double sYj[],
                             double cYi[],double cYj[],
                             double sftu,double sfcu,
                             double rftu,double rfcu,
                             double cftu,double cfcu,
                             double rAi[],double rYi[],
                             double *sN,double *sM,
                             double *rN,double *rM,
                             double *cN,double *cM,
                             double *gN,double *gM);

double allowableshearofrclong(struct element elem,
                              struct materials m,
                              int axis,
                              double Qd,
                              double Md);
double allowableshearofrcshort(struct element elem,
                               struct materials m,
                               int axis,
                               double Qd,
                               double Md);
double ultimateshearofrc(struct element elem,
                         struct materials m,
                         int axis,
                         double Nd,double Qd,double Md);

double allowableshearofsrclong(struct element elem,
                               struct materials m,
                               int axis,
                               int haxis,
                               double srcQd,
                               double srcMd);
double allowableshearofsrcshort(struct element elem,
                                struct materials m,
                                int axis,
                                int haxis,
                                double srcQd,
                                double srcMd);
double ultimateshearofsrc(struct element elem,
                          struct materials m,
                          int axis,
                          int haxis,
                          double rcQd,
                          double rcMd);

double allowultimshearofrcwall(struct element elem,
                               struct materials m,
                               double l, /*LENGTH[cm]*/
                               int period);

double allowabletensionofs(double sF,
                           struct section sect);
double allowablecompressionofs(double E,double F,
                               double lk,
                               struct section sect,int axis);
double allowablebendingofhstrong(int period,
                                 double Nd,
                                 double lb,
                                 struct materials m,
                                 struct section sect,int axis);
double allowablebendingofsweak(int period,
                               double Nd,
                               double lk,
                               struct materials m,
                               struct section sect,int axis);
double allowultimshearofhstrong(int period,
                                double F,
                                struct section sect,int axis);
double allowultimshearofhweak(int period,
                              double F,
                              struct section sect,int axis);

double coeffA(int nrect,struct materialrect rect[]);
double coeffI(int nrect,struct materialrect rect[],int axis,
              double *yg,double *yt,double *yc);
double coeffZ(int nrect,struct materialrect rect[],int axis);
double coeffi(int nrect,struct materialrect rect[],int axis);

double allowablebendingofpccolumn(struct element elem,
                                  struct materials m,int axis,
                                  double Nd);
double allowablebendingofpcgirder(struct element elem,
                                  struct materials m,int axis);
double allowablebendingofpcgirderonmid(struct element elem,
                                       struct materials m,
                                       double length,int axis);
double ultimatebendingofpccolumn(struct element elem,
                                 struct materials m,int axis,
                                 double Nd);
double ultimatebendingofpcgirder(struct element elem,
                                 struct materials m,int axis,
                                 int tensedside);
double allowableshearofpccolumn(struct element elem,
                                struct materials m,int axis,
                                double Nd);
double ultimateshearofpccolumn(struct element elem,
                               struct materials m,int axis,
                               double Nd,double Qd,double Md);
double ultimateshearofpcgirder(struct element elem,
                               struct materials m,int axis,
                               int tensedside,
                               double Nd,double Qd,double Md);

/*GLOBAL VARIABLES*/
FILE *fout0;
struct materials gmaterial;                      /*GLOBAL MATERIAL.*/

/**********/
/*SRCAN001.C:SRC SKELETON FOR WIN32 SINCE 1995.12.21.JUNSATO.*/
/*HEADER:SRCAM001.H*/
/*LAST CHANGE:1997.10.16.*/

/*FOR HAKODATE UNIV....PC GIRDER,COLUMN.....UNDER CONSTRUCTION.*/

/*FOR HOSHO CHURCH.....STRONG AXIS OF H IN COLUMN IS ALWAYS y=X.*/
/*                     STRONG AXIS OF H IN GIRDER IS ALWAYS y=Z.*/
/*                     ONLY H IN S,SRC.*/
/*FOR RC COLUMN,RC GIRDER,RC WALL,SRC COLUMN.*/

/*FILES:  SRCAN004.C*/
/*        SRCAM005.H*/
/*       SRCAN02.LST*/
/*        INC08C.INL*/
/*        INC08C.OTL*/
/*        INC08C.OHX*/
/*        INC08C.OHY*/
/*OUTPUT:SRCAN02.TST*/

/*SECTION FROM SECTION LIST FILE.*/
/*STRESS FROM "INC08C.OTL,OHX,OHY" OUTPUT OF "FRAME3EM.BAS".*/

/*SECTION LIST FILE SPECIFICATION:*/
/*"CODE" CODENUMBER MATERIALTYPE ELEMENTTYPE*/
/*"SRECT" BOTTOMLEFTX BOTTOMLEFTY TOPRIGHTX TOPRIGHTY*/
/*"REINS" AREA X Y*/
/*"CRECT" BOTTOMLEFTX BOTTOMLEFTY TOPRIGHTX TOPRIGHTY*/
/*"HOOPS" HOOPRATEX HOOPRATEY*/
/*"THICK" THICKNESS*/
/*"SREIN" SHEARREINRATE*/
/*"WRECT" WINDOWWIDTH WINDOWHEIGHT*/

/*
#define PROJECT "test"
#define HENSHINX 0.15
#define HENSHINY 0.15
#define CMQFILE     "c:\\cdocs\\srcan\\data\\srcan.cmq"
#define OUTPUTFILE  "srcan.tst"
#define RESULTLIST  "srcan.rlt"
#define RATELIST    "srcan.rat"
*/
/*
#define PROJECT "ken"
#define HENSHINX 0.15
#define HENSHINY 0.15
#define INPUTFILE   "e:\\kenyu\\kenyu26.inl"
#define CMQFILE     "e:\\kenyu\\kenyu26.cmq"
#define RESULTFILEX "e:\\kenyu\\kenyu26.ohx"
#define RESULTFILEY "e:\\kenyu\\kenyu26.ohy"
#define RESULTFILEZ "e:\\kenyu\\kenyu26.otl"
#define OUTPUTFILE  "e:\\kenyu\\ken26.tst"
#define RESULTLIST  "e:\\kenyu\\ken26.rlt"
#define SECTIONLIST "e:\\kenyu\\ken26.lst"
#define RATELIST    "e:\\kenyu\\ken26.rat"
*/

#define PROJECT "aki"
#define HENSHINX 0.15
#define HENSHINY 0.15
#define CMQFILE     "c:\\cdocs\\srcan\\data\\srcan.cmq"

int srcan001(char fname[])
{
  int ii,jj,ie,is,ic,ia,in,na,n1,n2,n3,ielem,isect,ni,nj;
  int hoko;                 /*DIRECTION OF HORIZONTAL LOAD 0:X 1:Y.*/
  int sign;        /*SIGN OF HORIZONTAL LOAD 0:POSITIVE 1:NEGATIVE.*/
  double dsign;
  double lN,eN; /*LONG,EARTHQUAKE STRESS.*/
  double sN,sQ[2][2],sMt,sM[2][2];  /*s:SHORT [END][AXIS]*/
  double uN,uQ[2][2],uQe[2][2],uMt,uM[2][2],Qp[2];  /*u:ULTIMATE*/
  double Ma,Mu,Muhead,Mutail,Qa,Qu; /*a:ALLOABLE u:ULTIMATE*/
  double Ns,Ms,Nr,Mr,Nc,Mc;            /*s:STEEL r:REIN c:CONCRETE.*/
  double nfact=1.0,mfact=1.0,qfact=2.0,Co=0.2,Ds=0.3,Fes;
  double wafact=2.0,wufact=3.0; /*FACTORS FOR WALL.*/
  double h,l; /*WALL HEIGHT,LENGTH*/
  double h0[2],l0; /*INNER HEIGHT,LENGTH*/
  double facei,facej,face[2][2];
  struct section sect1;
  long int codelist[MAXSECT];
  int cmqcode[MAXCMQ];
  double Mo[MAXCMQ],Mm;
  double rate,ratema,ratemu,rateqa,ratequ;
  double marate[MAXSECT],murate[MAXSECT];
  double qarate[MAXSECT],qurate[MAXSECT];
  double wrate1,wrate2; /*RATE OF WINDOW.*/
  /*int check;*/

  FILE *fin,*fcmq=NULL,*flist,*fsafe,*frate;
  FILE *fz=NULL,*fx=NULL,*fy=NULL;
  char txt[400],non[80],str[256],prj[256];
  char strhoko[2][10]={"X","Y"};
  /*char straxis[2][10]={"x","y"};*/
  /*char strend[2][10]={"i","j"};*/
  /*char strsign[2][10]={"+","-"};*/
  char strhugo[2][10]={"正","負"};
  char strQ[5][2][256],strM[3][2][256]; /*STRINGS FOR OUTPUT.*/
  int nnode,nelem,nsect,ncmq=0,soffset;
  double As,Ar,Ac,Yg/*,E,G*/;
  double jis=1.1;                                  /*1.1=JIS STEEL.*/
  struct element elem;
  struct stress *pstress[2]; /*HEAD,TAIL*/

  struct structnode *nodes;
  struct element *elems;
  struct section *sects;

  int tensedside;

  char dir[]=DIRECTORY;

  strcpy(prj,PROJECT);
  sprintf(txt,"PROJECT:%s\n",prj);

  sprintf(non,"FILES  :%s.TST,RLT,RAT\n",fname);
  strcat(txt,non);

  sprintf(non,"Continue or Cancel"); strcat(txt,non);

  if(MessageBox(NULL,txt,"SRCan",MB_OKCANCEL)==IDCANCEL)
  {
    return 0;
  }

  strcpy(txt,fname);
  strcat(txt,".tst");
  fout0=fopen(txt,"w"); if(fout0==NULL) return 0;

  strcpy(txt,fname);
  strcat(txt,".rlt");
  fsafe=fopen(txt,"w"); if(fsafe==NULL) return 0;

  strcpy(txt,fname);
  strcat(txt,".rat");
  frate=fopen(txt,"w"); if(frate==NULL) return 0;

  /*LIST UP FILES*/
  GetDlgItemText((wmenu.childs+2)->hwnd,ID_INPUTFILEZ,non,80);
  fprintf(fout0,"断面算定 \"S,RC,SRC\" 長期,短期,終局\n");
  fprintf(fout0,"入力        =%s\n",non);
  fprintf(fsafe,"%s\n",non);

  GetDlgItemText((wmenu.childs+2)->hwnd,ID_OUTPUTFILEZ,non,80);
  fprintf(fout0,"鉛直時結果  =%s\n",non);
  fprintf(fsafe,"%s\n",non);
  GetDlgItemText((wmenu.childs+2)->hwnd,ID_OUTPUTFILEX,non,80);
  fprintf(fout0,"水平時結果 X=%s\n",non);
  fprintf(fsafe,"%s\n",non);
  GetDlgItemText((wmenu.childs+2)->hwnd,ID_OUTPUTFILEY,non,80);
  fprintf(fout0,"           Y=%s\n",non);
  fprintf(fsafe,"%s\n",non);

  GetDlgItemText((wmenu.childs+2)->hwnd,ID_SECTIONFILE,non,80);
  fprintf(fout0,"断面入力    =%s\n",non);
  fprintf(fsafe,"%s\n",non);

  fin  =fgetstofopen(dir,"r",ID_INPUTFILEZ);
  flist=fgetstofopen(dir,"r",ID_SECTIONFILE);
  fz   =fgetstofopen("\0","r",ID_OUTPUTFILEZ);
  fx   =fgetstofopen("\0","r",ID_OUTPUTFILEX);
  fy   =fgetstofopen("\0","r",ID_OUTPUTFILEY);
  if(fin==NULL)   return 0;
  if(flist==NULL) return 0;
  if(fz==NULL)    return 0;
  if(fx==NULL)    return 0;
  if(fy==NULL)    return 0;

  fcmq=fopen(CMQFILE,"r");      /*if(fcmq==NULL)  return 0;*/

  /*INITIAL SETTING*/
  fgetinitial2(fin,&nnode,&nelem,&nsect);  /*INPUT INITIAL.*/

  sects=(struct section *)malloc(nsect*sizeof(struct section));
  if(sects==NULL) return 0;
  nodes=(struct structnode *)malloc(nnode*sizeof(struct structnode));
  if(nodes==NULL) return 0;
  elems=(struct element *)malloc(nelem*sizeof(struct element));
  if(elems==NULL) return 0;

  inputfiletomemory2(fin,&nnode,&nelem,&nsect,sects,nodes,elems);

  comments(fout0);
  /*if(MessageBox(NULL,"Continue or Cancel","SRCan",
                MB_OKCANCEL)==IDCANCEL) return 0;*/

  for(is=0;is<nsect;is++) /*SECTIONS INTO MEMORY.*/
  {
    getsrcansection(flist,(sects+is)->code,(sects+is));
  }
  nsect=getcodelist(flist,codelist); /*nsect CHANGED.*/

  if(nsect>MAXSECT)
  {
    MessageBox(NULL,"MAIN:SECTIONS OVERFLOW.","SRCan",MB_OK);
    return 0;
  }
  for(is=0;is<nsect;is++)
  {
    qarate[is]=0.0;
    qurate[is]=0.0;
    marate[is]=0.0;
    murate[is]=0.0;
  }

  if(fcmq!=NULL)
  {
    ncmq=getcmqlist(fcmq,cmqcode,Mo);
    if(ncmq>MAXCMQ)
    {
      MessageBox(NULL,"MAIN:CMQ OVERFLOW.","SRCan",MB_OK);
      return 0;
    }
  }
/*
if(MessageBox(NULL,"Pass 1","SRCan",MB_OKCANCEL)==IDCANCEL) return 0;
*/
  ie=0;
  while(1) /*CALCULATION*/
  {
    if(ie>=nelem) break;
    currentpivot(ie+1,nelem);

    elem=*(elems+ie);            /*ALL ELEMENTS POINTING EACH SECT.*/

    initializestress(&(elem.head.x));           /*INITIALIZE STRESS*/
    initializestress(&(elem.tail.x));
    initializestress(&(elem.head.y));
    initializestress(&(elem.tail.y));
    initializestress(&(elem.head.z));
    initializestress(&(elem.tail.z));

    n1=fgetstress(fz,&ielem,&isect,&ni,&(elem.head.z));
    n2=fgetstress(fx,&ielem,&isect,&ni,&(elem.head.x));
    n3=fgetstress(fy,&ielem,&isect,&ni,&(elem.head.y));
    if(n1==0 || n2==0 || n3==0) break; /*END OF FILE.*/

    if(n1==9 && n2==9 && n3==9)
    {
      if(ielem)
      {
        n1=fgetstress(fz,&ielem,&isect,&nj,&(elem.tail.z));
        n2=fgetstress(fx,&ielem,&isect,&nj,&(elem.tail.x));
        n3=fgetstress(fy,&ielem,&isect,&nj,&(elem.tail.y));

        if(n1<7 || n2<7 || n3<7) break;
        else ie++; /*COUNT ELEMENTS.*/

        rateqa=0.0;
        ratequ=0.0;
        ratema=0.0;
        ratemu=0.0;

        /*fprintf(stderr,"ELEM:%d SECT:%d\r",ielem,isect);*/

        /*soffset=getsrcansection(flist,isect,&sect1);*/
        soffset=elem.sect->soff;
        sect1=*(elem.sect);

        elem.sect->code=isect;

        if(elem.cmqcode!=0) /*CMQ DATA.*/
        {
          for(ic=0;ic<ncmq;ic++)
          {
            if(elem.cmqcode==cmqcode[ic]) elem.Mo=Mo[ic];
          }
        }

        if(soffset && sect1.etype!=WALL && sect1.etype!=SLAB)
        {
          h=100.0*elementlength(elem); /*[cm]*/

          if(sect1.etype==COLUMN)
          {
            h0[SX]=h-elem.sect->face[SX][HEAD]
                    -elem.sect->face[SX][TAIL]; /*INNER LENGTH.*/
            h0[SY]=h-elem.sect->face[SY][HEAD]
                    -elem.sect->face[SY][TAIL]; /*INNER LENGTH.*/
          }
          else
          {
            h0[SX]=h-elem.sect->face[SX][HEAD]
                    -elem.sect->face[SX][TAIL]; /*INNER LENGTH.*/
          }

          fprintf(fout0,"-------------------------------------");
          fprintf(fout0,"-------------------------------------");
          fprintf(fout0,"-----------------------------------\n");
          fprintf(fout0,"部材:%4d ",ielem);
          fprintf(fout0,"始端:%3d 終端:%3d ",ni,nj);

          fprintf(fout0,"断面:%3d",isect);
          if(sect1.stype==S)        fprintf(fout0,"=Ｓ");
          else if(sect1.stype==RC)  fprintf(fout0,"=ＲＣ");
          else if(sect1.stype==SRC) fprintf(fout0,"=ＳＲＣ");
          else if(sect1.stype==PC)  fprintf(fout0,"=ＰＣ");
          else                      fprintf(fout0,"=  ");
          if(sect1.etype==COLUMN)      fprintf(fout0,"柱 ");
          else if(sect1.etype==GIRDER) fprintf(fout0,"大梁 ");
          else if(sect1.etype==BEAM)   fprintf(fout0,"小梁 ");
          else if(sect1.etype==BRACE)  fprintf(fout0,"筋違 ");
          else                         fprintf(fout0,"不明 ");

          fprintf(fout0,"材長=%.1f[cm] Mx内法=%.1f[cm]",h,h0[SX]);
          if(sect1.etype==COLUMN)
          {
            fprintf(fout0," My内法=%.1f[cm]\n",h0[SY]);
          }
          else
          {
            fprintf(fout0,"\n");
          }

          /*MATERIAL INDEPENDENT ON PERIOD.*/
          gmaterial.sE=2100000.0; /*[kgf/cm2]*/
          gmaterial.rE=2100000.0;
          /*gmaterial.cE=210000.0;*/
          gmaterial.cE=gmaterial.rE/15.0;
          gmaterial.sF=3300.0*jis;                   /*STEEL SN490B*/
          gmaterial.Fc=240.0;                      /*CONCRETE Fc240*/

          translatesection(sect1,gmaterial,&As,&Ac,&Ar,&Yg,SX);
          /*fprintf(fout0,"As=%7.2f[cm2] ",As);
          fprintf(fout0,"Ar=%7.2f[cm2] ",Ar);
          fprintf(fout0,"Ac=%7.2f[cm2]\n",Ac);*/

          fprintf(fout0,"応力[tf,tfm]         N");
          fprintf(fout0,"      Qxi      Qxj      Qyi      Qyj");
          fprintf(fout0,"       Mt");
          fprintf(fout0,"      Mxi      Mxj      Myi      Myj\n");

          fprintf(fout0,"鉛直時Z     : %8.3f",elem.head.z.N);
          fprintf(fout0," %8.3f %8.3f",
                  elem.head.z.Q[0],elem.tail.z.Q[0]);
          fprintf(fout0," %8.3f %8.3f",
                  elem.head.z.Q[1],elem.tail.z.Q[1]);
          fprintf(fout0," %8.3f",elem.head.z.Mt);
          fprintf(fout0," %8.3f %8.3f",
                  elem.head.z.M[0],elem.tail.z.M[0]);
          fprintf(fout0," %8.3f %8.3f\n",
                  elem.head.z.M[1],elem.tail.z.M[1]);

          fprintf(fout0,"水平時X     : %8.3f",elem.head.x.N);
          fprintf(fout0," %8.3f %8.3f",
                  elem.head.x.Q[0],elem.tail.x.Q[0]);
          fprintf(fout0," %8.3f %8.3f",
                  elem.head.x.Q[1],elem.tail.x.Q[1]);
          fprintf(fout0," %8.3f",elem.head.x.Mt);
          fprintf(fout0," %8.3f %8.3f",
                  elem.head.x.M[0],elem.tail.x.M[0]);
          fprintf(fout0," %8.3f %8.3f\n",
                  elem.head.x.M[1],elem.tail.x.M[1]);

          fprintf(fout0,"水平時Y     : %8.3f",elem.head.y.N);
          fprintf(fout0," %8.3f %8.3f",
                  elem.head.y.Q[0],elem.tail.y.Q[0]);
          fprintf(fout0," %8.3f %8.3f",
                  elem.head.y.Q[1],elem.tail.y.Q[1]);
          fprintf(fout0," %8.3f",elem.head.y.Mt);
          fprintf(fout0," %8.3f %8.3f",
                  elem.head.y.M[0],elem.tail.y.M[0]);
          fprintf(fout0," %8.3f %8.3f\n",
                  elem.head.y.M[1],elem.tail.y.M[1]);
          fprintf(fout0,"\n");

          /*if(check){gets(non); if(non[0]!='\0') return 0;}*/

          /*LONG...................................................*/
          /*MATERIAL FOR LONG*/
          gmaterial.sft=3300.0/1.5*jis;               /*STEEL SM490*/
          gmaterial.sfc=-gmaterial.sft; /*FOR SRC*/
          gmaterial.sfb=3300.0/1.5*jis; /*FOR SRC*/     /*b:BENDING*/
          gmaterial.sfs=3300.0/1.5/sqrt(3.0)*jis;         /*s:SHEAR*/

if(!strcmp(prj,"ken"))
{
  if(isect==141 || isect==341)
  {
    gmaterial.sft=3000.0/1.5*jis;               /*STEEL SM490*/
    gmaterial.sfc=-gmaterial.sft; /*FOR SRC*/
    gmaterial.sfb=3000.0/1.5*jis; /*FOR SRC*/     /*b:BENDING*/
    gmaterial.sfs=3000.0/1.5/sqrt(3.0)*jis;         /*s:SHEAR*/
  }
}

if(!strcmp(prj,"saka") || !strcmp(prj,"kuro"))
{
  gmaterial.rft=2000.0*jis;     /*REINFORCEMENT SD295*/
  gmaterial.wft=2000.0;
}
else{  gmaterial.rft=2200.0*jis; /*REINFORCEMENT SD345*/
       gmaterial.wft=2000.0;}
          gmaterial.rfc=-gmaterial.rft;
          /*n=15.0;*/                                  /*RATE Er/Ec*/
          gmaterial.cfc=-160.0/2.0;                /*CONCRETE Fc240*/
          gmaterial.cfs=11.1/1.5;

          gmaterial.pcF=500.0;                  /*PC CONCRETE Fc500*/
          gmaterial.pcfc=gmaterial.pcF/3.0;         /*+:COMPRESSION*/
          if(gmaterial.pcfc>210.0) gmaterial.pcfc=210.0;
          gmaterial.pcft=0.0;
          gmaterial.pcfsu=7.5+1.5/100.0*gmaterial.pcF;
          if(gmaterial.pcfsu>16.5) gmaterial.pcfsu=16.5;

          gmaterial.pcfco=gmaterial.pcF*0.45;         /*+:COMPRESSION*/
          if(gmaterial.pcfco>210.0) gmaterial.pcfco=210.0;
          gmaterial.pcfto=-0.07*gmaterial.pcfco;
          if(gmaterial.pcfto<-21.0) gmaterial.pcfto=-21.0;

          gmaterial.stfact=0.8; /*EFFECTIVE RATE OF PRESTRESS.*/
          gmaterial.stft[STRNDA]=0.8* 8000.0; /*PC STRAND A [kgf/cm2]*/
          gmaterial.stft[STRNDB]=0.8* 9500.0; /*PC STRAND B*/
          gmaterial.stft[STRNDC]=0.8*11000.0; /*PC STRAND C No.2*/
          gmaterial.stft[STRNDW]=0.8*16000.0; /*PC STRAND WIRE SWPR7B*/
          gmaterial.stftu[STRNDA]= 8000.0;
          gmaterial.stftu[STRNDB]= 9500.0;
          gmaterial.stftu[STRNDC]=11000.0;
          gmaterial.stftu[STRNDW]=16000.0;

          face[0][0]=1.0;
          face[0][1]=1.0;
          face[1][0]=1.0;
          face[1][1]=1.0;

  if(!strcmp(prj,"ken")) /*FACE RATES.*/
  {
    if(elem.sect->etype==COLUMN) na=1;
    else                         na=0;
    for(ia=0;ia<=na;ia++) /*FOR AXIS*/
    {
      if(elem.head.z.M[ia]==0.0) face[0][ia]=1.0;
      else
      {
        face[0][ia]=(fabs(elem.head.z.M[ia]/elem.head.z.Q[1-ia])
                    -elem.sect->face[ia][HEAD]/100.0)
                    /fabs(elem.head.z.M[ia]/elem.head.z.Q[1-ia]);
      }
      if(elem.tail.z.M[ia]==0.0) face[1][ia]=1.0;
      else
      {
        face[1][ia]=(fabs(elem.tail.z.M[ia]/elem.tail.z.Q[1-ia])
                    -elem.sect->face[ia][TAIL]/100.0)
                    /fabs(elem.tail.z.M[ia]/elem.tail.z.Q[1-ia]);
      }
    }
  }
          fprintf(fout0,"長期        : %8.3f",elem.head.z.N);
          fprintf(fout0," %8.3f %8.3f",
                  elem.head.z.Q[0],elem.tail.z.Q[0]);
          fprintf(fout0," %8.3f %8.3f",
                  elem.head.z.Q[1],elem.tail.z.Q[1]);
          fprintf(fout0," %8.3f",elem.head.z.Mt);
          fprintf(fout0," %8.3f %8.3f",
                  face[0][0]*elem.head.z.M[0],
                  face[1][0]*elem.tail.z.M[0]);
          fprintf(fout0," %8.3f %8.3f\n",
                  face[0][1]*elem.head.z.M[1],
                  face[1][1]*elem.tail.z.M[1]);

          for(ii=0;ii<=4;ii++)
          {
            for(jj=0;jj<=1;jj++) sprintf(strQ[ii][jj],"\0");
          }
          for(ii=0;ii<=2;ii++)
          {
            for(jj=0;jj<=1;jj++) sprintf(strM[ii][jj],"\0");
          }

          if(elem.sect->etype==COLUMN) na=1;
          else                         na=0;
          for(ia=0;ia<=na;ia++) /*FOR AXIS*/
          {
            pstress[0]=&(elem.head.z); pstress[1]=&(elem.tail.z);

            lN=1000.0*elem.head.z.N; /*LONG N[kgf].*/

            Ma=0.0;
if(!strcmp(prj,"ken")) /*KENYU GIRDERS.*/
{
if(sect1.stype==RC && sect1.etype==GIRDER) lN=0.0;
}
            if(sect1.stype==S)                                  /*S*/
            {
              if(sect1.etype==GIRDER || sect1.etype==BEAM)
              {
                Ma=allowablebendingofhstrong(PLONG,lN,h,
                                             gmaterial,sect1,SX);
              }
              if(sect1.etype==COLUMN && ia==SX)
              {
                Ma=allowablebendingofhstrong(PLONG,lN,h,
                                             gmaterial,sect1,SX);
              }
              if(sect1.etype==COLUMN && ia==SY)
              {
                Ma=allowablebendingofsweak(PLONG,lN,h,
                                           gmaterial,sect1,SY);
              }
            }
            else if(sect1.stype==RC)                           /*RC*/
            {
/*for(lN=-39600.0;lN<=158400.0;lN=lN+19800.0)*/
{
              Ma=allowablebendingofrc(elem,gmaterial,ia,(-lN));
}
            }
            else if(sect1.stype==SRC)                         /*SRC*/
            {
/*for(lN=-394.4;lN<=524.7;lN=lN+9.191)*/
{
              Ma=allowablebendingofsrc(elem,gmaterial,ia,(-lN));
}
            }
            else if(sect1.stype==PC)                           /*PC*/
            {
              if(sect1.etype==COLUMN)
              {
                Ma=allowablebendingofpccolumn(elem,gmaterial,ia,lN);
              }
              else if(sect1.etype==GIRDER || sect1.etype==BEAM)
              {
                Ma=allowablebendingofpcgirder(elem,gmaterial,ia);
              }
            }
            if(Ma<0.1) Ma=0.1;

            sprintf(str," %8.3f %8.3f",Ma/100000.0,Ma/100000.0);
            strcat(strM[1][ia],str);
            sprintf(str," %8.3f %8.3f",
                    fabs(face[0][ia]*pstress[0]->M[ia]*100000.0/Ma),
                    fabs(face[1][ia]*pstress[1]->M[ia]*100000.0/Ma));
            strcat(strM[2][ia],str);

            /*
            fprintf(fout0,"Ma%s=%.3f[tfm] ",straxis[ia],Ma/100000.0);
            fprintf(fout0," Mi/Ma=%.3f Mj/Ma=%.3f\n",
                    fabs(pstress[0]->M[ia]*100000.0/Ma),
                    fabs(pstress[1]->M[ia]*100000.0/Ma));
            */
            rate=fabs(face[0][ia]*pstress[0]->M[ia]*100000.0/Ma);
            if(rate>ratema) ratema=rate;
            if(rate>marate[soffset-1]) marate[soffset-1]=rate;
            rate=fabs(face[1][ia]*pstress[1]->M[ia]*100000.0/Ma);
            if(rate>ratema) ratema=rate;
            if(rate>marate[soffset-1]) marate[soffset-1]=rate;

            for(in=HEAD;in<=TAIL;in++) /*FOR END*/
            {
              Qa=0.0;
              if(sect1.stype==S)                                /*S*/
              {
                if(sect1.etype==GIRDER || sect1.etype==BEAM)
                {
                  Qa=allowultimshearofhstrong(PLONG,gmaterial.sF,
                                              sect1,SX);
                }
                if(sect1.etype==COLUMN && ia==SX)
                {
                  Qa=allowultimshearofhstrong(PLONG,gmaterial.sF,
                                              sect1,SX);
                }
                if(sect1.etype==COLUMN && ia==SY)
                {
                  Qa=allowultimshearofhweak(PLONG,gmaterial.sF,
                                            sect1,SY);
                }
              }
              else if(sect1.stype==RC)                         /*RC*/
              {
                Qa=allowableshearofrclong(elem,gmaterial,ia,
                                          pstress[in]->Q[1-ia]
                                          *1000.0,
                                          pstress[in]->M[ia]
                                          *100000.0);
              }
              else if(sect1.stype==SRC)                       /*SRC*/
              {
                if(sect1.etype==GIRDER || sect1.etype==BEAM)
                {
                  Qa=allowableshearofsrclong(elem,gmaterial,
                                             SX,HSTRONG,
                                             pstress[in]->Q[1-ia]
                                             *1000.0,
                                             pstress[in]->M[ia]
                                             *100000.0);
                }
                if(sect1.etype==COLUMN && ia==SX)
                {
                  Qa=allowableshearofsrclong(elem,gmaterial,
                                             SX,HSTRONG,
                                             pstress[in]->Q[1-ia]
                                             *1000.0,
                                             pstress[in]->M[ia]
                                             *100000.0);
                }
                if(sect1.etype==COLUMN && ia==SY)
                {
                  Qa=allowableshearofsrclong(elem,gmaterial,
                                             SY,HWEAK,
                                             pstress[in]->Q[1-ia]
                                             *1000.0,
                                             pstress[in]->M[ia]
                                             *100000.0);
                }
              }
              else if(sect1.stype==PC)                         /*PC*/
              {
                Qa=allowableshearofpccolumn(elem,gmaterial,ia,lN);
              }
              sprintf(str," %8.3f",Qa/1000.0);
              strcat(strQ[1][1-ia],str);

              /*
              fprintf(fout0,"Qa%s=%.3f[tf]",straxis[1-ia],Qa/1000.0);
              */
              if(Qa<0.1) Qa=0.1;

              rate=fabs(pstress[in]->Q[1-ia]*1000.0/Qa);
              if(rate>rateqa) rateqa=rate;
              if(rate>qarate[soffset-1]) qarate[soffset-1]=rate;

              sprintf(str," %8.3f",rate);
              strcat(strQ[2][1-ia],str);
              /*fprintf(fout0," Q%s/Qa=%.3f\n",strend[in],rate);*/
            }
          }
          if(sect1.etype==GIRDER || sect1.etype==BEAM)
          {
            sprintf(strQ[1][0],"                  ");
            sprintf(strQ[2][0],"                  ");
            sprintf(strM[1][1],"                  ");
            sprintf(strM[2][1],"                  ");
          }

          fprintf(fout0,"      許容値:         %s%s",
                  strQ[1][0],strQ[1][1]);
          fprintf(fout0,"         %s%s\n",
                  strM[1][0],strM[1][1]);
          fprintf(fout0,"      到達率:         %s%s",
                  strQ[2][0],strQ[2][1]);
          fprintf(fout0,"         %s%s\n",
                  strM[2][0],strM[2][1]);

          /*FOR MIDPOINT BENDING Mm.*/
          if(sect1.stype==PC &&
             elem.sect->etype==GIRDER &&
             elem.cmqcode!=0)
          {
            Mm=elem.Mo-0.5*(elem.tail.z.M[0]-elem.head.z.M[0]);
            fprintf(fout0,"        中央:Mo=%8.3f Mm=%8.3f",
                    elem.Mo,Mm);

            if(Mm<=0.0)
            {
              fprintf(fout0," 上端引張\n\n\n");
            }
            else
            {
              Ma=allowablebendingofpcgirderonmid(elem,gmaterial,
                                                 h0[SX],0);
              if(Ma<0.1) Ma=0.1;

              rate=fabs(Mm*100000.0/Ma);
              if(rate>ratema) ratema=rate;
              if(rate>marate[soffset-1]) marate[soffset-1]=rate;

              fprintf(fout0,"\n");
              fprintf(fout0,"      許容値:               %8.3f\n",
                      Ma/100000.0);
              fprintf(fout0,"      到達率:               %8.3f\n",
                      rate);
            }
          }
          else if(elem.sect->etype==GIRDER)
          {
            /*fprintf(fout0,"\n\n\n");*/
          }

          fprintf(fout0,"\n");
          /*if(check){gets(non); if(non[0]!='\0') return 0;}*/

          /*SHORT..................................................*/
          /*MATERIAL FOR SHORT*/
          gmaterial.sft=3300.0;                       /*STEEL SM490*/
          gmaterial.sfc=-gmaterial.sft;
          gmaterial.sfb=3300.0; /*FOR SRC*/             /*b:BENDING*/
          gmaterial.sfs=3300.0/sqrt(3.0);                 /*s:SHEAR*/
if(!strcmp(prj,"saka") || !strcmp(prj,"kuro"))
{
  gmaterial.rft=3000.0*jis;     /*REINFORCEMENT SD295*/
  gmaterial.wft=3000.0;
}
else{  gmaterial.rft=3500.0*jis; /*REINFORCEMENT SD345*/
       gmaterial.wft=3000.0;}
          gmaterial.rfc=-gmaterial.rft;
          /*n=15.0;*/                                  /*RATE Er/Ec*/
          gmaterial.cfc=-160.0;                          /*CONCRETE*/
          gmaterial.cfs=11.1;

          for(hoko=0;hoko<=1;hoko++) /*LOAD DIRECTION X,Y.*/
          {
            if(sect1.stype==PC) break;  /*PC ROUTE 3a : SKIP SHORT.*/

            if(hoko==0)
            {
              pstress[0]=&(elem.head.x); pstress[1]=&(elem.tail.x);
            }
            if(hoko==1)
            {
              pstress[0]=&(elem.head.y); pstress[1]=&(elem.tail.y);
            }

            for(sign=0;sign<=1;sign++) /*LOAD POSITIVE,NEGATIVE.*/
            {
              for(ii=0;ii<=4;ii++)
              {
                for(jj=0;jj<=1;jj++) sprintf(strQ[ii][jj],"\0");
              }
              for(ii=0;ii<=2;ii++)
              {
                for(jj=0;jj<=1;jj++) sprintf(strM[ii][jj],"\0");
              }

              sprintf(strQ[0][0],"短期%s%s方向 :",
                      strhoko[hoko],strhugo[sign]);

              if(elem.sect->etype==COLUMN) na=1;
              else                         na=0;
              for(ia=0;ia<=na;ia++) /*FOR AXIS*/
              {
                if(pstress[0]->M[ia]==0.0) facei=1.0;       /*RATE.*/
                else
                {
                  facei=(fabs(pstress[0]->M[ia]/pstress[0]->Q[1-ia])
                         -elem.sect->face[ia][HEAD]/100.0)
                        /fabs(pstress[0]->M[ia]/pstress[0]->Q[1-ia]);
                }
                if(pstress[1]->M[ia]==0.0) facej=1.0;
                else
                {
                  facej=(fabs(pstress[1]->M[ia]/pstress[1]->Q[1-ia])
                         -elem.sect->face[ia][TAIL]/100.0)
                        /fabs(pstress[1]->M[ia]/pstress[1]->Q[1-ia]);
                }

                if(sign==0) dsign=1.0; /*POSITIVE*/
                if(sign==1) dsign=-1.0;/*NEGATIVE*/

                sN=1000.0*(elem.head.z.N
                           +dsign*nfact*(pstress[HEAD]->N));
                sQ[HEAD][1-ia]=(elem.head.z.Q[1-ia]
                                +dsign
                                *qfact*(pstress[HEAD]->Q[1-ia]))
                               *1000.0;
                sQ[TAIL][1-ia]=(elem.tail.z.Q[1-ia]
                                +dsign
                                *qfact*(pstress[TAIL]->Q[1-ia]))
                               *1000.0;
                sMt=(elem.head.z.Mt+dsign*mfact*(pstress[HEAD]->Mt))
                    *100000.0;
                sM[HEAD][ia]=(elem.head.z.M[ia]
                              +dsign*facei
                              *mfact*(pstress[HEAD]->M[ia]))
                             *100000.0;
                sM[TAIL][ia]=(elem.tail.z.M[ia]
                              +dsign*facej
                              *mfact*(pstress[TAIL]->M[ia]))
                             *100000.0;

                if(ia==0)
                {
                  sprintf(str," %8.3f",sN/1000.0);
                  strcat(strQ[0][0],str);
                }
                if(sect1.etype==GIRDER || sect1.etype==BEAM)
                {
                  strcat(strQ[0][0],"                  ");
                }
                sprintf(str," %8.3f %8.3f",
                        sQ[HEAD][1-ia]/1000.0,
                        sQ[TAIL][1-ia]/1000.0);
                strcat(strQ[0][1-ia],str);

                if(ia==0) sprintf(strM[0][0]," %8.3f",sMt/100000.0);
                sprintf(str," %8.3f %8.3f",
                        sM[HEAD][ia]/100000.0,
                        sM[TAIL][ia]/100000.0);
                strcat(strM[0][ia],str);

                Ma=0.0;

if(!strcmp(prj,"ken")) /*KENYU GIRDERS.*/
{
if(sect1.stype==RC && sect1.etype==GIRDER) sN=0.0;
}
                if(sect1.stype==S)                              /*S*/
                {
                  if(sect1.etype==GIRDER || sect1.etype==BEAM)
                  {
                    Ma=allowablebendingofhstrong(PSHORT,sN,h,
                                                 gmaterial,sect1,SX);
                  }
                  if(sect1.etype==COLUMN && ia==SX)
                  {
                    Ma=allowablebendingofhstrong(PSHORT,sN,h,
                                                 gmaterial,sect1,SX);
                  }
                  if(sect1.etype==COLUMN && ia==SY)
                  {
                    Ma=allowablebendingofsweak(PSHORT,sN,h,
                                               gmaterial,sect1,SY);
                  }
                }
                else if(sect1.stype==RC)                       /*RC*/
                {
if(0)
{
if(isect==3011) sN=0.0; /*HAKODATE FG.*/
if(isect==3012) sN=0.0;
if(isect==302) sN=0.0;
if(isect==303) sN=0.0;
if(isect==304) sN=0.0;
if(isect==305) sN=0.0;
}
                  Ma=allowablebendingofrc(elem,gmaterial,ia,(-sN));
                }
                else if(sect1.stype==SRC)                     /*SRC*/
                {
                  Ma=allowablebendingofsrc(elem,gmaterial,ia,(-sN));
                }
                if(Ma<0.1) Ma=0.1;

                sprintf(str," %8.3f %8.3f",Ma/100000.0,Ma/100000.0);
                strcat(strM[1][ia],str);
                sprintf(str," %8.3f %8.3f",fabs(sM[HEAD][ia]/Ma),
                                           fabs(sM[TAIL][ia]/Ma));
                strcat(strM[2][ia],str);
                /*
                fprintf(fout0,"Ma%s=%.3f[tfm]",
                        straxis[ia],Ma/100000.0);
                fprintf(fout0," Mi/Ma=%.3f Mj/Ma=%.3f\n",
                        fabs(sM[HEAD][ia]/Ma),
                        fabs(sM[TAIL][ia]/Ma));
                */
                rate=fabs(sM[HEAD][ia]/Ma);
                if(rate>ratema) ratema=rate;
                if(rate>marate[soffset-1]) marate[soffset-1]=rate;
                rate=fabs(sM[TAIL][ia]/Ma);
                if(rate>ratema) ratema=rate;
                if(rate>marate[soffset-1]) marate[soffset-1]=rate;

                for(in=HEAD;in<=TAIL;in++) /*FOR END*/
                {
                  Qa=0.0;

                  if(sect1.stype==S)                            /*S*/
                  {
                    if(sect1.etype==GIRDER || sect1.etype==BEAM)
                    {
                      Qa=allowultimshearofhstrong(PSHORT,gmaterial.sF,
                                                  sect1,SX);
                    }
                    if(sect1.etype==COLUMN && ia==SX)
                    {
                      Qa=allowultimshearofhstrong(PSHORT,gmaterial.sF,
                                                  sect1,SX);
                    }
                    if(sect1.etype==COLUMN && ia==SY)
                    {
                      Qa=allowultimshearofhweak(PSHORT,gmaterial.sF,
                                                sect1,SY);
                    }
                  }
                  else if(sect1.stype==RC)                     /*RC*/
                  {
                    Qa=allowableshearofrcshort(elem,gmaterial,ia,
                                               sQ[in][1-ia],
                                               sM[in][ia]);
                  }
                  else if(sect1.stype==SRC)                   /*SRC*/
                  {
                    if(sect1.etype==GIRDER || sect1.etype==BEAM)
                    {
                      Qa=allowableshearofsrcshort(elem,gmaterial,
                                                  SX,HSTRONG,
                                                  sQ[in][1-ia],
                                                  sM[in][ia]);
                    }
                    if(sect1.etype==COLUMN && ia==SX)
                    {
                      Qa=allowableshearofsrcshort(elem,gmaterial,
                                                  SX,HSTRONG,
                                                  sQ[in][1-ia],
                                                  sM[in][ia]);
                    }
                    if(sect1.etype==COLUMN && ia==SY)
                    {
                      Qa=allowableshearofsrcshort(elem,gmaterial,
                                                  SY,HWEAK,
                                                  sQ[in][1-ia],
                                                  sM[in][ia]);
                    }
                  }
                  sprintf(str," %8.3f",Qa/1000.0);
                  strcat(strQ[1][1-ia],str);

                  /*fprintf(fout0,"Qa%s=%8.3f[tf]",
                          straxis[1-ia],Qa/1000.0);*/
                  if(Qa<0.1) Qa=0.1;

                  rate=fabs(sQ[in][1-ia]/Qa);
                  if(rate>rateqa) rateqa=rate;
                  if(rate>qarate[soffset-1])
                  {
                    qarate[soffset-1]=rate;
                  }
                  sprintf(str," %8.3f",rate);
                  strcat(strQ[2][1-ia],str);

                  /*
                  fprintf(fout0," Q%s/Qa=%.3f\n",strend[in],rate);
                  */
                }
              }

              if(sect1.etype==GIRDER || sect1.etype==BEAM)
              {
                sprintf(strQ[1][0],"                  ");
                sprintf(strQ[2][0],"                  ");
                sprintf(strM[1][1],"                  ");
                sprintf(strM[2][1],"                  ");
              }
              fprintf(fout0,"%s%s%s%s\n",
                      strQ[0][0],strQ[0][1],strM[0][0],strM[0][1]);
              fprintf(fout0,"      許容値:         %s%s",
                      strQ[1][0],strQ[1][1]);
              fprintf(fout0,"         %s%s\n",
                      strM[1][0],strM[1][1]);
              fprintf(fout0,"      到達率:         %s%s",
                      strQ[2][0],strQ[2][1]);
              fprintf(fout0,"         %s%s\n",
                      strM[2][0],strM[2][1]);
            }
            fprintf(fout0,"\n");
            /*if(check){gets(non); if(non[0]!='\0') return 0;}*/
          }

          /*ULTIMATE...............................................*/
          /*fprintf(stderr,"ELEM:%d SECT:%d\n",ielem,isect);*/

          /*MATERIAL FOR ULTIMATE*/
          gmaterial.sftu=3300.0*jis;                  /*STEEL SM490*/
          gmaterial.sfcu=-gmaterial.sftu;
          gmaterial.rftu=3500.0*jis;          /*REINFORCEMENT SD345*/
          gmaterial.rfcu=-gmaterial.rftu;
          gmaterial.wfp=3000.0*jis;           /*REINFORCEMENT SD295*/

          for(hoko=0;hoko<=1;hoko++) /*LOAD DIRECTION X,Y.*/
          {
            if(hoko==0)
            {
              Fes=1.0+(HENSHINX-0.15)/0.15*0.5;
              pstress[0]=&(elem.head.x); pstress[1]=&(elem.tail.x);
            }
            if(hoko==1)
            {
              Fes=1.0+(HENSHINY-0.15)/0.15*0.5;
              pstress[0]=&(elem.head.y); pstress[1]=&(elem.tail.y);
            }

            for(sign=0;sign<=1;sign++) /*LOAD POSITIVE,NEGATIVE.*/
            {
              for(ii=0;ii<=4;ii++)
              {
                for(jj=0;jj<=1;jj++) sprintf(strQ[ii][jj],"\0");
              }
              for(ii=0;ii<=2;ii++)
              {
                for(jj=0;jj<=1;jj++) sprintf(strM[ii][jj],"\0");
              }

              sprintf(strQ[0][0],"終局%s%s方向 :",
                      strhoko[hoko],strhugo[sign]);

              if(elem.sect->etype==COLUMN) na=1;
              else                         na=0;
              for(ia=0;ia<=na;ia++) /*FOR AXIS*/
              {
                if(pstress[0]->M[ia]==0.0) facei=1.0;       /*RATE.*/
                else
                {
                  facei=(fabs(pstress[0]->M[ia]/pstress[0]->Q[1-ia])
                         -elem.sect->face[ia][HEAD]/100.0)
                        /fabs(pstress[0]->M[ia]/pstress[0]->Q[1-ia]);
                }
                if(pstress[1]->M[ia]==0.0) facej=1.0;
                else
                {
                  facej=(fabs(pstress[1]->M[ia]/pstress[1]->Q[1-ia])
                         -elem.sect->face[ia][TAIL]/100.0)
                        /fabs(pstress[1]->M[ia]/pstress[1]->Q[1-ia]);
                }

                if(sign==0) dsign=1.0; /*POSITIVE*/
                if(sign==1) dsign=-1.0;/*NEGATIVE*/

                uN=1000.0*(elem.head.z.N
                           +dsign*(Fes*Ds/Co)*(pstress[HEAD]->N));
                uQ[HEAD][1-ia]=(elem.head.z.Q[1-ia]
                                +dsign
                                *(Fes*Ds/Co)*(pstress[HEAD]->Q[1-ia]))
                               *1000.0;
                uQ[TAIL][1-ia]=(elem.tail.z.Q[1-ia]
                                +dsign
                                *(Fes*Ds/Co)*(pstress[TAIL]->Q[1-ia]))
                               *1000.0;
                uMt=100000.0*(elem.head.z.Mt
                              +dsign*(Fes*Ds/Co)*(pstress[HEAD]->Mt));
                uM[HEAD][ia]=(elem.head.z.M[ia]
                              +dsign*facei
                              *(Fes*Ds/Co)*(pstress[HEAD]->M[ia]))
                             *100000.0;
                uM[TAIL][ia]=(elem.tail.z.M[ia]
                              +dsign*facej
                              *(Fes*Ds/Co)*(pstress[TAIL]->M[ia]))
                             *100000.0;

                if(ia==0)
                {
                  sprintf(str," %8.3f",uN/1000.0);
                  strcat(strQ[0][0],str);
                }
                if(sect1.etype==GIRDER || sect1.etype==BEAM)
                {
                  strcat(strQ[0][0],"                  ");
                }
                sprintf(str," %8.3f %8.3f",
                        uQ[HEAD][1-ia]/1000.0,
                        uQ[TAIL][1-ia]/1000.0);
                strcat(strQ[0][1-ia],str);

                if(ia==0) sprintf(strM[0][0]," %8.3f",uMt/100000.0);
                sprintf(str," %8.3f %8.3f",
                        uM[HEAD][ia]/100000.0,
                        uM[TAIL][ia]/100000.0);
                strcat(strM[0][ia],str);

                Mu=0.0;
                Qp[HEAD]=0.0;
                Qp[TAIL]=0.0;

                if(sect1.stype==S)                              /*S*/
                {
                  Mu=ultimatebendingofsrc(elem,gmaterial,ia,(-uN),
                                          &Ns,&Ms,&Nr,&Mr,&Nc,&Mc);
                }
                else if(sect1.stype==RC)                       /*RC*/
                {
if(0)
{
if(isect==3011) uN=0.0; /*HAKODATE FG.*/
if(isect==3012) uN=0.0;
if(isect==302) uN=0.0;
if(isect==303) uN=0.0;
if(isect==304) uN=0.0;
if(isect==305) uN=0.0;
if(isect==331) uN=0.0;
}
                  Mu=ultimatebendingofsrc(elem,gmaterial,ia,(-uN),
                                          &Ns,&Ms,&Nr,&Mr,&Nc,&Mc);
                }
                else if(sect1.stype==SRC)                     /*SRC*/
                {
                  Mu=ultimatebendingofsrc(elem,gmaterial,ia,(-uN),
                                          &Ns,&Ms,&Nr,&Mr,&Nc,&Mc);
                }
                else if(sect1.stype==PC)                       /*PC*/
                {
                  if(sect1.etype==COLUMN)
                  {
                    Mu=ultimatebendingofpccolumn(elem,gmaterial,ia,uN);
                  }
                  if(sect1.etype==GIRDER)
                  {
                    if(uM[HEAD][ia]>=0)     tensedside=LOWER;
                    else if(uM[HEAD][ia]<0) tensedside=UPPER;
                    Muhead=ultimatebendingofpcgirder(elem,gmaterial,ia,
                                                     tensedside);

                    if(uM[TAIL][ia]>=0)     tensedside=UPPER;
                    else if(uM[TAIL][ia]<0) tensedside=LOWER;
                    Mutail=ultimatebendingofpcgirder(elem,gmaterial,ia,
                                                     tensedside);
                  }
                }

                if(sect1.stype==PC && sect1.etype==GIRDER)
                {
                  if(h0[ia]<=0.0) /*ALL LENGTH IN FACE.*/
                  {
                    Muhead=0.1;
                    Mutail=0.1;
                  }

                  Qp[HEAD]=2.0*fabs(Muhead)/h0[ia];
                  Qp[TAIL]=2.0*fabs(Mutail)/h0[ia];

                  if(Muhead<0.1) Muhead=0.1;
                  if(Mutail<0.1) Mutail=0.1;

                  sprintf(str," %8.3f %8.3f",Muhead/100000.0,
                                             Mutail/100000.0);
                  strcat(strM[1][ia],str);
                  sprintf(str," %8.3f %8.3f",fabs(uM[HEAD][ia]/Muhead),
                                             fabs(uM[TAIL][ia]/Mutail));
                  strcat(strM[2][ia],str);

                  rate=fabs(uM[HEAD][ia]/Muhead);
                  if(rate>ratemu) ratemu=rate;
                  if(rate>murate[soffset-1]) murate[soffset-1]=rate;
                  rate=fabs(uM[TAIL][ia]/Mutail);
                  if(rate>ratemu) ratemu=rate;
                  if(rate>murate[soffset-1]) murate[soffset-1]=rate;
                }
                else
                {
                  Qp[HEAD]=2.0*fabs(Mu)/h0[ia];
                  Qp[TAIL]=2.0*fabs(Mu)/h0[ia];

                  if(Mu<0.1) Mu=0.1;

                  sprintf(str," %8.3f %8.3f",Mu/100000.0,Mu/100000.0);
                  strcat(strM[1][ia],str);
                  sprintf(str," %8.3f %8.3f",fabs(uM[HEAD][ia]/Mu),
                                             fabs(uM[TAIL][ia]/Mu));
                  strcat(strM[2][ia],str);
                  /*
                  fprintf(fout0,"Mu%s=%.3f[tfm]",
                          straxis[ia],Mu/100000.0);
                  fprintf(fout0," Mi/Mu=%.3f Mj/Mu=%.3f\n",
                          fabs(uM[HEAD][ia]/Mu),
                          fabs(uM[TAIL][ia]/Mu));
                  */
                  rate=fabs(uM[HEAD][ia]/Mu);
                  if(rate>ratemu) ratemu=rate;
                  if(rate>murate[soffset-1]) murate[soffset-1]=rate;
                  rate=fabs(uM[TAIL][ia]/Mu);
                  if(rate>ratemu) ratemu=rate;
                  if(rate>murate[soffset-1]) murate[soffset-1]=rate;
                }

                for(in=HEAD;in<=TAIL;in++) /*FOR END*/
                {
                  Qu=0.0;

                  if(sect1.stype==S)                            /*S*/
                  {
                    if(sect1.etype==GIRDER || sect1.etype==BEAM)
                    {
                      Qu=allowultimshearofhstrong(PULTIMATE,
                                                  gmaterial.sF,
                                                  sect1,SX);
                    }
                    if(sect1.etype==COLUMN && ia==SX)
                    {
                      Qu=allowultimshearofhstrong(PULTIMATE,
                                                  gmaterial.sF,
                                                  sect1,SX);
                    }
                    if(sect1.etype==COLUMN && ia==SY)
                    {
                      Qu=allowultimshearofhweak(PULTIMATE,
                                                gmaterial.sF,
                                                sect1,SY);
                    }
                  }
                  else if(sect1.stype==RC)                     /*RC*/
                  {
                    Qu=ultimateshearofrc(elem,gmaterial,ia,
                                         (-uN),
                                         uQ[in][1-ia],
                                         uM[in][ia]);
                  }
                  else if(sect1.stype==SRC)                   /*SRC*/
                  {
                    if(sect1.etype==GIRDER || sect1.etype==BEAM)
                    {
                      Qu=ultimateshearofsrc(elem,gmaterial,
                                            SX,HSTRONG,
                                            uQ[in][1-ia],
                                            uM[in][ia]);
                      /*APPROXIMATION.Q,M MUST BE OF RC.*/
                    }
                    if(sect1.etype==COLUMN && ia==SX)
                    {
                      Qu=ultimateshearofsrc(elem,gmaterial,
                                            SX,HSTRONG,
                                            uQ[in][1-ia],
                                            uM[in][ia]);
                      /*APPROXIMATION.Q,M MUST BE OF RC.*/
                    }
                    if(sect1.etype==COLUMN && ia==SY)
                    {
                      Qu=ultimateshearofsrc(elem,gmaterial,
                                            SY,HWEAK,
                                            uQ[in][1-ia],
                                            uM[in][ia]);
                      /*APPROXIMATION.Q,M MUST BE OF RC.*/
                    }
                  }
                  else if(sect1.stype==PC)                     /*PC*/
                  {
                    if(sect1.etype==COLUMN)
                    {
                      Qu=ultimateshearofpccolumn(elem,gmaterial,ia,
                                                 uN,
                                                 uQ[in][1-ia],
                                                 uM[in][ia]);
                    }
                    if(sect1.etype==GIRDER)
                    {
                      if(in==HEAD)
                      {
                        if(uM[HEAD][ia]>=0)     tensedside=LOWER;
                        else if(uM[HEAD][ia]<0) tensedside=UPPER;
                      }
                      else if(in==TAIL)
                      {
                        if(uM[TAIL][ia]>=0)     tensedside=UPPER;
                        else if(uM[TAIL][ia]<0) tensedside=LOWER;
                      }
                      Qu=ultimateshearofpcgirder(elem,gmaterial,ia,
                                                 tensedside,
                                                 uN,
                                                 uQ[in][1-ia],
                                                 uM[in][ia]);
                    }
                  }

                  sprintf(str," %8.3f",Qu/1000.0);
                  strcat(strQ[1][1-ia],str);

                  /*fprintf(fout0,"Qu%s=%8.3f[tf]",
                          straxis[1-ia],Qu/1000.0);*/
                  if(Qu<0.1) Qu=0.1;

                  rate=fabs(uQ[in][1-ia]/Qu);
                  if(rate>ratequ) ratequ=rate;
                  if(rate>qurate[soffset-1])
                  {
                    qurate[soffset-1]=rate;
                  }
                  sprintf(str," %8.3f",rate);
                  strcat(strQ[2][1-ia],str);

                  /*
                  fprintf(fout0," Q%s/Qu=%.3f\n",strend[in],rate);
                  */

                  if(sect1.stype==RC)
                  {
                    Qp[in]*=1.1; /*"CENTER STANDARD" P.288*/
                  }

                  if(in==HEAD)
                  {
                    Qp[in]+=1000.0*fabs(elem.head.z.Q[1-ia]);
                  }
                  if(in==TAIL)
                  {
                    Qp[in]+=1000.0*fabs(elem.tail.z.Q[1-ia]);
                  }

                  if(sect1.stype==PC || sect1.stype==RC)
                  {
                    if(in==HEAD)
                    {
                      uQe[HEAD][1-ia]=(elem.head.z.Q[1-ia]
                                       +dsign
                                       *(Fes*2.0)*(pstress[HEAD]->Q[1-ia]))
                                      *1000.0;
                    }
                    if(in==TAIL)
                    {
                      uQe[TAIL][1-ia]=(elem.tail.z.Q[1-ia]
                                       +dsign
                                       *(Fes*2.0)*(pstress[TAIL]->Q[1-ia]))
                                      *1000.0;
                    }

                    if(Qp[in]>fabs(uQe[in][1-ia])) Qp[in]=uQe[in][1-ia];
                  }

                  sprintf(str," %8.3f",Qp[in]/1000.0);
                  strcat(strQ[3][1-ia],str);

                  rate=fabs(Qp[in]/Qu);
                  if(rate>ratequ) ratequ=rate;
                  if(rate>qurate[soffset-1]) qurate[soffset-1]=rate;
                  sprintf(str," %8.3f",rate);
                  strcat(strQ[4][1-ia],str);
                }
              }

              if(sect1.etype==GIRDER || sect1.etype==BEAM)
              {
                sprintf(strQ[1][0],"                  ");
                sprintf(strQ[2][0],"                  ");
                sprintf(strQ[3][0],"                  ");
                sprintf(strQ[4][0],"                  ");
                sprintf(strM[1][1],"                  ");
                sprintf(strM[2][1],"                  ");
              }
              fprintf(fout0,"%s%s%s%s\n",
                      strQ[0][0],strQ[0][1],strM[0][0],strM[0][1]);
              fprintf(fout0,"      終局値:         %s%s",
                      strQ[1][0],strQ[1][1]);
              fprintf(fout0,"         %s%s\n",
                      strM[1][0],strM[1][1]);
              fprintf(fout0,"      到達率:         %s%s",
                      strQ[2][0],strQ[2][1]);
              fprintf(fout0,"         %s%s\n",
                      strM[2][0],strM[2][1]);
              fprintf(fout0,"  機構形成時:         %s%s\n",
                      strQ[3][0],strQ[3][1]);
              fprintf(fout0,"      到達率:         %s%s\n",
                      strQ[4][0],strQ[4][1]);
            }
            fprintf(fout0,"\n");
            /*if(check){gets(non); if(non[0]!='\0') return 0;}*/
          }

          fprintf(fout0,"MAX:Q/Qa=%.5f",rateqa);
          fprintf(fout0," Q/Qu=%.5f",ratequ);
          fprintf(fout0," M/Ma=%.5f",ratema);
          fprintf(fout0," M/Mu=%.5f",ratemu);
          fprintf(fout0,"\n");

          if(rateqa>1.0 || ratequ>1.0 ||
             ratema>1.0 || ratemu>1.0)
          {
            fprintf(fsafe,"NG ELEM:%4d SECT:%4d ",ielem,isect);
            fprintf(fsafe,"MAX:Q/Qa=%.5f Q/Qu=%.5f ",rateqa,ratequ);
            fprintf(fsafe,"M/Ma=%.5f M/Mu=%.5f\n",ratema,ratemu);
          }
          /*if(rateqa<=1.0 && ratequ<=1.0 &&
             ratema<=1.0 && ratemu<=1.0)
          {
            fprintf(fsafe,"ELEM:%4d SECT:%4d ",ielem,isect);
            fprintf(fsafe,"MAX:Q/Qa=%.5f Q/Qu=%.5f ",rateqa,ratequ);
            fprintf(fsafe,"M/Ma=%.5f M/Mu=%.5f\n",ratema,ratemu);
          }*/

          if(frate!=NULL)
          {
            fprintf(frate,"ELEM: %5d SECT: %5d  %.5f %.5f %.5f %.5f\n",
                    ielem,isect,rateqa,ratequ,ratema,ratemu);
          }

          /*if(check){gets(non); if(non[0]!='\0') return 0;}*/
        }

        if(soffset && sect1.etype==WALL)/*...........SHEAR OF WALL.*/
        {
          Qa=0.0;
          Qu=0.0;

          /*fprintf(stderr,"ELEM:%d SECT:%d\n",ielem,isect);*/

          fprintf(fout0,"-------------------------------------");
          fprintf(fout0,"-------------------------------------");
          fprintf(fout0,"-----------------------------------\n");
          fprintf(fout0,"部材:%4d ",ielem);
          fprintf(fout0,"始端:%3d 終端:%3d ",ni,nj);

          fprintf(fout0,"断面:%3d",isect);
          if(sect1.stype==RC)  fprintf(fout0,"=ＲＣ壁 ");
          else                 fprintf(fout0,"=不明壁 ");

          fprintf(fout0,"ＲＣ壁厚:%.0f[mm] ",sect1.thick*10.0);

          /*MATERIAL SHORT*/
          gmaterial.Fc=240.0;                      /*CONCRETE Fc240*/
          gmaterial.wft=3000.0;
          gmaterial.wfp=3300.0;
          gmaterial.cfs=11.1;

          h=100.0*wallheight(elem); /*[cm]*/
          h0[1]=h-elem.sect->face[1][HEAD]-elem.sect->face[1][TAIL];
          if(h0[1]<0.0) h0[1]=0.0;
          l=100.0*walllength(elem); /*[cm]*/
          l0=l-elem.sect->face[0][HEAD]-elem.sect->face[0][TAIL];
          if(l0<0.0) l0=0.0;
          fprintf(fout0,"内法:L=%.1f[cm] H=%.1f[cm] ",l0,h0[1]);

          wrate1=(elem.sect->wlength)/l0; /*RATE OF WINDOW.*/
          wrate2=sqrt((elem.sect->wlength)*(elem.sect->wheight)
                      /(l0*h0[1]));
          if(wrate1>=wrate2) elem.sect->windowrate=wrate1;
          else               elem.sect->windowrate=wrate2;
          fprintf(fout0,"開口率:1-r=%.3f\n",elem.sect->windowrate);

          l0/=2.0; /*1 OF 2 BRACES.*/

          Qa=allowultimshearofrcwall(elem,gmaterial,l0,PSHORT);
          /*Qu=allowultimshearofrcwall(elem,gmaterial,l0,PULTIMATE);*/

          for(hoko=0;hoko<=1;hoko++)
          {
            if(hoko==0) eN=elem.head.x.N;
            if(hoko==1) eN=elem.head.y.N;

            sQ[0][0]=eN*l/sqrt(l*l+h*h)*1000.0;
            fprintf(fout0,"水平時%s:",strhoko[hoko]);
            fprintf(fout0,"N=%7.3f[tf] 水平成分:Qh=%7.3f[tf] ",
                    eN,sQ[0][0]/1000.0);

            fprintf(fout0,"短期=%.1f×水平:",wafact);
            fprintf(fout0,"Qs=%7.3f[tf]",(wafact*sQ[0][0])/1000.0);

if(0){ /*FOR ROUTE 1*/
            fprintf(fout0," 終局=%.1f×水平:",wufact);
            fprintf(fout0,"Qu=%7.3f[tf]",(wufact*sQ[0][0])/1000.0);
}/*FOR ROUTE 1*/

            /*fprintf(fout0,"\n");*/

            rate=fabs(wafact*sQ[0][0]/Qa);
            if(rate>rateqa) rateqa=rate;
            if(rate>qarate[soffset-1]) qarate[soffset-1]=rate;
            fprintf(fout0," Qs/Qa=%7.5f",rate);

if(0){ /*FOR ROUTE 1*/
            rate=fabs(wufact*sQ[0][0]/Qu);
            if(rate>ratequ) ratequ=rate;
            if(rate>qurate[soffset-1]) qurate[soffset-1]=rate;
            fprintf(fout0," Qs/Qu=%7.5f",rate);
}/*FOR ROUTE 1*/

            fprintf(fout0,"\n");
          }

          if(rateqa>1.0 || ratequ>1.0)
          {
            fprintf(fsafe,"NG ELEM:%4d SECT:%4d ",ielem,isect);
            fprintf(fsafe,"MAX:Q/Qa=%.5f Q/Qu=%.5f\n",rateqa,ratequ);
          }
          /*if(rateqa<=1.0 && ratequ<=1.0)
          {
            fprintf(fsafe,"ELEM:%4d SECT:%4d ",ielem,isect);
            fprintf(fsafe,"MAX:Q/Qa=%.5f Q/Qu=%.5f\n",rateqa,ratequ);
          }*/

          if(frate!=NULL)
          {
            fprintf(frate,"ELEM: %5d SECT: %5d  %.5f %.5f\n",
                    ielem,isect,rateqa,ratequ);
          }
        }

        if(soffset && sect1.etype==SLAB)/*...........SHEAR OF SLAB.*/
        {
          /*UNDER CONSTRUCTION*/
        }
      }
    }
  }
  fprintf(fout0,"=====================================");
  fprintf(fout0,"=====================================");
  fprintf(fout0,"===================================\n");
  fprintf(fout0,"各断面種別の許容,終局曲げ到達率の最大値\n\n");

  rateqa=0.0; ratequ=0.0;
  ratema=0.0; ratemu=0.0;

  for(is=1;is<=nsect;is++)
  {
    soffset=getsrcansection(flist,codelist[is-1],&sect1);
    if(soffset)
    {
      if(sect1.etype!=WALL && sect1.etype!=SLAB)
      {
        fprintf(fout0,"断面記号:%4d",sect1.code);

        translatesection(sect1,gmaterial,&As,&Ac,&Ar,&Yg,SX);
        fprintf(fout0," As=%7.2f[cm2]",As);
        fprintf(fout0," Ar=%7.2f[cm2]",Ar);
        fprintf(fout0," Ac=%8.2f[cm2]",Ac);

        fprintf(fout0," MAX:Q/Qa=%7.5f",qarate[soffset-1]);
        /*fprintf(fout0," Q/Qu=%7.5f",qurate[soffset-1]);*/
        fprintf(fout0," M/Ma=%7.5f",marate[soffset-1]);
        /*fprintf(fout0," M/Mu=%7.5f",murate[soffset-1]);*/
        fprintf(fout0,"\n");
      }
      if(sect1.etype==WALL)
      {
        fprintf(fout0,"断面記号:%4d",sect1.code);

        fprintf(fout0," t=%4.1f[cm]",sect1.thick);
        fprintf(fout0," ps=%7.5f",sect1.shearrein[0]);

        fprintf(fout0,"                          ");
        fprintf(fout0," MAX:Q/Qa=%7.5f",qarate[soffset-1]);
        /*fprintf(fout0," Q/Qu=%7.5f",qurate[soffset-1]);*/
        fprintf(fout0,"\n");
      }
      if(qarate[soffset-1]>rateqa) rateqa=qarate[soffset-1];
      if(qurate[soffset-1]>ratequ) ratequ=qurate[soffset-1];
      if(marate[soffset-1]>ratema) ratema=marate[soffset-1];
      if(murate[soffset-1]>ratemu) ratemu=murate[soffset-1];
    }
  }

  sprintf(txt,"安全率の最大値");
  sprintf(non," Q/Qa=%7.5f Q/Qu=%7.5f",rateqa,ratequ);
  strcat(txt,non);
  sprintf(non," M/Ma=%7.5f M/Mu=%7.5f",ratema,ratemu);
  strcat(txt,non);
  fprintf(fout0,"\n%s\n",txt);

  fclose(fin);
  fclose(flist);
  fclose(fx);
  fclose(fy);
  fclose(fz);
  fclose(fout0);
  fclose(fsafe);
  fclose(frate);

  sprintf(non,"\nCompleted."); strcat(txt,non);
  MessageBox(NULL,txt,"SRCan",MB_OK);
  return 1;
}/*srcan001*/

int createyieldsurface(struct arclmframe *af)
{
  int ii,ie,is,ia,na;
  /*int hoko;*/                 /*DIRECTION OF HORIZONTAL LOAD 0:X 1:Y.*/
  int sign;        /*SIGN OF HORIZONTAL LOAD 0:POSITIVE 1:NEGATIVE.*/
  /*double lN,eN;*/ /*LONG,EARTHQUAKE STRESS.*/
  /*double sN,sQ[2][2],sMt,sM[2][2];*/  /*s:SHORT [END][AXIS]*/
  double uN,/*uQ[2][2],uQe[2][2],uMt,uM[2][2],*/Qp[2];  /*u:ULTIMATE*/
  double Nu,Nmax,Nmin,Mtu;
  double Mu,Qa,Qu; /*a:ALLOABLE u:ULTIMATE*/
  double Ns,Ms,Nr,Mr,Nc,Mc;            /*s:STEEL r:REIN c:CONCRETE.*/
  double h,l; /*WALL HEIGHT,LENGTH*/
  double h0[2],l0; /*INNER HEIGHT,LENGTH*/
  /*double facei,facej,face[2][2];*/
  struct section sect1;
  long int codelist[MAXSECT];
  double wrate1,wrate2; /*RATE OF WINDOW.*/
  /*int check;*/

  FILE *fin,*flist;
  char txt[400],non[80],/*str[256],*/prj[256];
  int nnode,nelem,nsect,soffset;
  /*double E,G;*/
  double jis=1.1;                                  /*1.1=JIS STEEL.*/
  struct element elem;
  struct stress *pstress[2]; /*HEAD,TAIL*/
  struct structnode *nodes;
  struct element *elems;
  struct section *sects;
  int tensedside;
  char dir[]=DIRECTORY;

  struct owire *aelem;
  struct surface ys;

  /*BEGINNING MESSAGE*/
  strcpy(prj,PROJECT);
  sprintf(txt,"Project Code:%s\n",prj);
  sprintf(non,"Continue or Cancel"); strcat(txt,non);
  if(MessageBox(NULL,txt,"Create Yield Surface",MB_OKCANCEL)
     ==IDCANCEL) return 0;

  /*LIST UP FILES*/
  fin  =fgetstofopen(dir,"r",ID_INPUTFILE);
  flist=fgetstofopen(dir,"r",ID_SECTIONFILE);
  if(fin==NULL)   return 0;
  if(flist==NULL) return 0;

  /*INITIAL SETTING*/
  fgetinitial2(fin,&nnode,&nelem,&nsect);  /*INPUT INITIAL.*/
  if(nnode!=af->nnode ||
     nelem!=af->nelem ||
     nsect!=af->nsect)
  {
    MessageBox(NULL,"FILE ERROR.","SURFACE",MB_OK);
    return 0;
  }

  sects=(struct section *)malloc(nsect*sizeof(struct section));
  if(sects==NULL) return 0;
  nodes=(struct structnode *)malloc(nnode*sizeof(struct structnode));
  if(nodes==NULL) return 0;
  elems=(struct element *)malloc(nelem*sizeof(struct element));
  if(elems==NULL) return 0;

  inputfiletomemory2(fin,&nnode,&nelem,&nsect,sects,nodes,elems);

  /*if(MessageBox(NULL,"Continue or Cancel","SRCan",
                MB_OKCANCEL)==IDCANCEL) return 0;*/

  for(is=0;is<nsect;is++) /*SECTIONS INTO MEMORY.*/
  {
    getsrcansection(flist,(sects+is)->code,(sects+is));
  }
  nsect=getcodelist(flist,codelist); /*nsect CHANGED.*/

  if(nsect>MAXSECT)
  {
    MessageBox(NULL,"MAIN:SECTIONS OVERFLOW.","SRCan",MB_OK);
    return 0;
  }

  ie=0;
  while(1) /*CALCULATION*/
  {
    if(ie>=nelem) break;
    currentpivot(ie+1,nelem);

    elem=*(elems+ie);            /*ALL ELEMENTS POINTING EACH SECT.*/
    aelem=af->elems+ie;

    /*INITIALIZE SURFACE*/
    ys.exp=0.0;
    for(ii=0;ii<6;ii++)
    {
      ys.fmax[ii]=0.0;
      ys.fmin[ii]=0.0;
    }

    initializestress(&(elem.head.x));           /*INITIALIZE STRESS*/
    initializestress(&(elem.tail.x));
    initializestress(&(elem.head.y));
    initializestress(&(elem.tail.y));
    initializestress(&(elem.head.z));
    initializestress(&(elem.tail.z));

        /*soffset=getsrcansection(flist,isect,&sect1);*/
        soffset=elem.sect->soff;
        sect1=*(elem.sect);

        /*elem.sect->code=isect;*/

        if(soffset && sect1.etype!=WALL && sect1.etype!=SLAB)
        {
          /*ELEMENT LENGTH*/
          h=100.0*elementlength(elem); /*[cm]*/

          if(sect1.etype==COLUMN)
          {
            h0[SX]=h-elem.sect->face[SX][HEAD]
                    -elem.sect->face[SX][TAIL]; /*INNER LENGTH.*/
            h0[SY]=h-elem.sect->face[SY][HEAD]
                    -elem.sect->face[SY][TAIL]; /*INNER LENGTH.*/
          }
          else
          {
            h0[SX]=h-elem.sect->face[SX][HEAD]
                    -elem.sect->face[SX][TAIL]; /*INNER LENGTH.*/
          }

          /*MATERIAL INDEPENDENT ON PERIOD.*/
          gmaterial.sE=2100000.0; /*[kgf/cm2]*/
          gmaterial.rE=2100000.0;
          /*gmaterial.cE=210000.0;*/
          gmaterial.cE=gmaterial.rE/15.0;
          gmaterial.sF=3300.0*jis;                   /*STEEL SN490B*/
          gmaterial.Fc=240.0;                      /*CONCRETE Fc240*/

          /*LONG...................................................*/
          /*MATERIAL FOR LONG*/
          gmaterial.pcF=500.0;                  /*PC CONCRETE Fc500*/
          gmaterial.pcfc=gmaterial.pcF/3.0;         /*+:COMPRESSION*/
          if(gmaterial.pcfc>210.0) gmaterial.pcfc=210.0;
          gmaterial.pcft=0.0;
          gmaterial.pcfsu=7.5+1.5/100.0*gmaterial.pcF;
          if(gmaterial.pcfsu>16.5) gmaterial.pcfsu=16.5;

          gmaterial.pcfco=gmaterial.pcF*0.45;         /*+:COMPRESSION*/
          if(gmaterial.pcfco>210.0) gmaterial.pcfco=210.0;
          gmaterial.pcfto=-0.07*gmaterial.pcfco;
          if(gmaterial.pcfto<-21.0) gmaterial.pcfto=-21.0;

          gmaterial.stfact=0.8; /*EFFECTIVE RATE OF PRESTRESS.*/
          gmaterial.stft[STRNDA]=0.8* 8000.0; /*PC STRAND A [kgf/cm2]*/
          gmaterial.stft[STRNDB]=0.8* 9500.0; /*PC STRAND B*/
          gmaterial.stft[STRNDC]=0.8*11000.0; /*PC STRAND C No.2*/
          gmaterial.stft[STRNDW]=0.8*16000.0; /*PC STRAND WIRE SWPR7B*/
          gmaterial.stftu[STRNDA]= 8000.0;
          gmaterial.stftu[STRNDB]= 9500.0;
          gmaterial.stftu[STRNDC]=11000.0;
          gmaterial.stftu[STRNDW]=16000.0;

          /*SHORT..................................................*/
          /*MATERIAL FOR SHORT*/
          gmaterial.sft=3300.0;                       /*STEEL SM490*/
          gmaterial.sfc=-gmaterial.sft;
          gmaterial.sfb=3300.0; /*FOR SRC*/             /*b:BENDING*/
          gmaterial.sfs=3300.0/sqrt(3.0);                 /*s:SHEAR*/
          if(!strcmp(prj,"saka") || !strcmp(prj,"kuro"))
          {
            gmaterial.rft=3000.0*jis;     /*REINFORCEMENT SD295*/
            gmaterial.wft=3000.0;
          }
          else{  gmaterial.rft=3500.0*jis; /*REINFORCEMENT SD345*/
                 gmaterial.wft=3000.0;}
          gmaterial.rfc=-gmaterial.rft;
          /*n=15.0;*/                                  /*RATE Er/Ec*/
          gmaterial.cfc=-160.0;                          /*CONCRETE*/
          gmaterial.cfs=11.1;

          /*ULTIMATE...............................................*/
          /*MATERIAL FOR ULTIMATE*/
          gmaterial.sftu=3300.0*jis;                  /*STEEL SM490*/
          gmaterial.sfcu=-gmaterial.sftu;
          gmaterial.rftu=3500.0*jis;          /*REINFORCEMENT SD345*/
          gmaterial.rfcu=-gmaterial.rftu;
          gmaterial.wfp=3000.0*jis;           /*REINFORCEMENT SD295*/

          ys.exp=1.5; /*1.0<exp<2.0*/

          /*AXIAL FORCE*/
          if(sect1.stype==S) /*S*/
          {
            /*Nmax=ultimatecompressionofh();*/
            /*Nmin=ultimatetensionofh();*/
          }
          else if(sect1.stype==RC) /*RC*/
          {
            /*Nmax=ultimatecompressionofrc();*/
            /*Nmin=ultimatetensionofrc();*/
          }
          else if(sect1.stype==SRC) /*SRC*/
          {
            /*Nmax=ultimatecompressionofsrc();*/
            /*Nmin=ultimatetensionofsrc();*/
          }
          else if(sect1.stype==PC) /*PC*/
          {
            if(sect1.etype==COLUMN)
            {
              /*Nmax=ultimatecompressionofpccolumn();*/
              /*Nmin=ultimatetensionofpccolumn();*/
            }
            else if(sect1.etype==GIRDER || sect1.etype==BEAM)
            {
              /*Nmax=ultimatecompressionofpcgirder();*/
              /*Nmin=ultimatetensionofpcgirder();*/
            }
          }

          /*Nmax,Nmin*/
          ys.fmax[0]=Nmax;
          ys.fmin[0]=Nmin;

          /*CENTER OF SURFACE : N=0.5(Nmax+Nmin)*/
          elem.head.x.N=0.5*(ys.fmax[0]+ys.fmin[0]);
          elem.tail.x.N=-elem.head.x.N;

          pstress[0]=&(elem.head.x); pstress[1]=&(elem.tail.x);

            for(sign=0;sign<=1;sign++) /*LOAD POSITIVE,NEGATIVE.*/
            {
              /*sprintf(strQ[0][0],"終局%s%s方向 :",
                      strhoko[hoko],strhugo[sign]);*/

              if(elem.sect->etype==COLUMN) na=1;
              else                         na=0;
              for(ia=0;ia<=na;ia++) /*FOR AXIS*/
              {
                uN=1000.0*(pstress[HEAD]->N);

                Mu=0.0;
                Qp[HEAD]=0.0;
                Qp[TAIL]=0.0;

                if(sect1.stype==S)                              /*S*/
                {
                  Mu=ultimatebendingofsrc(elem,gmaterial,ia,(-uN),
                                          &Ns,&Ms,&Nr,&Mr,&Nc,&Mc);
                }
                else if(sect1.stype==RC)                       /*RC*/
                {
                  Mu=ultimatebendingofsrc(elem,gmaterial,ia,(-uN),
                                          &Ns,&Ms,&Nr,&Mr,&Nc,&Mc);
                }
                else if(sect1.stype==SRC)                     /*SRC*/
                {
                  Mu=ultimatebendingofsrc(elem,gmaterial,ia,(-uN),
                                          &Ns,&Ms,&Nr,&Mr,&Nc,&Mc);
                }
                else if(sect1.stype==PC)                       /*PC*/
                {
                  if(sect1.etype==COLUMN)
                  {
                    Mu=ultimatebendingofpccolumn(elem,gmaterial,ia,uN);
                  }
                  if(sect1.etype==GIRDER)
                  {
                    if(sign==0)      tensedside=LOWER;
                    else if(sign==1) tensedside=UPPER;
                    Mu=ultimatebendingofpcgirder(elem,gmaterial,ia,
                                                 tensedside);
                  }
                }

                if(Mu<0.1) Mu=0.1;

                /*Mmax,Mmin*/
                if(sign==0) ys.fmax[ia+4]=Mu;
                if(sign==1) ys.fmin[ia+4]=-Mu;

                  Qu=0.0;
                  if(sect1.stype==S)                            /*S*/
                  {
                    if(sect1.etype==GIRDER || sect1.etype==BEAM)
                    {
                      Qu=allowultimshearofhstrong(PULTIMATE,
                                                  gmaterial.sF,
                                                  sect1,SX);
                    }
                    if(sect1.etype==COLUMN && ia==SX)
                    {
                      Qu=allowultimshearofhstrong(PULTIMATE,
                                                  gmaterial.sF,
                                                  sect1,SX);
                    }
                    if(sect1.etype==COLUMN && ia==SY)
                    {
                      Qu=allowultimshearofhweak(PULTIMATE,
                                                gmaterial.sF,
                                                sect1,SY);
                    }
                  }
                  else if(sect1.stype==RC)                     /*RC*/
                  {
                    Qu=ultimateshearofrc(elem,gmaterial,ia,
                                         (-uN),0.0,0.0);
                  }
                  else if(sect1.stype==SRC)                   /*SRC*/
                  {
                    if(sect1.etype==GIRDER || sect1.etype==BEAM)
                    {
                      Qu=ultimateshearofsrc(elem,gmaterial,
                                            SX,HSTRONG,
                                            0.0,0.0);
                      /*APPROXIMATION.Q,M MUST BE OF RC.*/
                    }
                    if(sect1.etype==COLUMN && ia==SX)
                    {
                      Qu=ultimateshearofsrc(elem,gmaterial,
                                            SX,HSTRONG,
                                            0.0,0.0);
                      /*APPROXIMATION.Q,M MUST BE OF RC.*/
                    }
                    if(sect1.etype==COLUMN && ia==SY)
                    {
                      Qu=ultimateshearofsrc(elem,gmaterial,
                                            SY,HWEAK,
                                            0.0,0.0);
                      /*APPROXIMATION.Q,M MUST BE OF RC.*/
                    }
                  }
                  else if(sect1.stype==PC)                     /*PC*/
                  {
                    if(sect1.etype==COLUMN)
                    {
                      Qu=ultimateshearofpccolumn(elem,gmaterial,ia,
                                                 uN,0.0,0.0);
                    }
                    if(sect1.etype==GIRDER)
                    {
                      if(sign==0) tensedside=LOWER;
                      if(sign==1) tensedside=UPPER;
                      Qu=ultimateshearofpcgirder(elem,gmaterial,ia,
                                                 tensedside,
                                                 uN,0.0,0.0);
                    }
                  }

                  /*fprintf(fout0,"Qu%s=%8.3f[tf]",
                          straxis[1-ia],Qu/1000.0);*/
                  if(Qu<0.1) Qu=0.1;

                /*Qmax,Qmin*/
                if(sign==0) ys.fmax[2-ia]=Mu;
                if(sign==1) ys.fmin[2-ia]=-Mu;

              }/*ia*/
            }/*sign*/

          /*TORTION*/
          if(sect1.stype==S) /*S*/
          {
            /*Mtu=ultimatetortionofh();*/
          }
          else if(sect1.stype==RC) /*RC*/
          {
            /*Mtu=ultimatetortionofrc();*/
          }
          else if(sect1.stype==SRC) /*SRC*/
          {
            /*Mtu=ultimatetortionofsrc();*/
          }
          else if(sect1.stype==PC) /*PC*/
          {
            if(sect1.etype==COLUMN)
            {
              /*Mtu=ultimatetortionofpccolumn();*/
            }
            else if(sect1.etype==GIRDER || sect1.etype==BEAM)
            {
              /*Mtu=ultimatetortionofpcgirder();*/
            }
          }

          /*Mtmax,Mtmin*/
          ys.fmax[3]=Mtu;
          ys.fmin[3]=-Mtu;

          /*if(check){gets(non); if(non[0]!='\0') return 0;}*/
        }

        if(soffset && sect1.etype==WALL)/*...........SHEAR OF WALL.*/
        {
          Qa=0.0;
          Qu=0.0;

          /*MATERIAL SHORT*/
          gmaterial.Fc=240.0;                      /*CONCRETE Fc240*/
          gmaterial.wft=3000.0;
          gmaterial.wfp=3300.0;
          gmaterial.cfs=11.1;

          h=100.0*wallheight(elem); /*[cm]*/
          h0[1]=h-elem.sect->face[1][HEAD]-elem.sect->face[1][TAIL];
          if(h0[1]<0.0) h0[1]=0.0;
          l=100.0*walllength(elem); /*[cm]*/
          l0=l-elem.sect->face[0][HEAD]-elem.sect->face[0][TAIL];
          if(l0<0.0) l0=0.0;
          /*fprintf(fout0,"内法:L=%.1f[cm] H=%.1f[cm] ",l0,h0[1]);*/

          wrate1=(elem.sect->wlength)/l0; /*RATE OF WINDOW.*/
          wrate2=sqrt((elem.sect->wlength)*(elem.sect->wheight)
                      /(l0*h0[1]));
          if(wrate1>=wrate2) elem.sect->windowrate=wrate1;
          else               elem.sect->windowrate=wrate2;
          /*
          fprintf(fout0,"開口率:1-r=%.3f\n",elem.sect->windowrate);
          */

          l0/=2.0; /*1 OF 2 BRACES.*/

          /*Qa=allowultimshearofrcwall(elem,gmaterial,l0,PSHORT);*/
          Qu=allowultimshearofrcwall(elem,gmaterial,l0,PULTIMATE);
          Nu=Qu*sqrt(l*l+h*h)/l;

          /*Nmax,Nmin*/
          ys.fmax[0]=Nu;
          ys.fmin[0]=-Nu;
        }

        if(soffset && sect1.etype==SLAB)/*...........SHEAR OF SLAB.*/
        {
          /*UNDER CONSTRUCTION*/
        }

    aelem->yield=ys;
  }

  fclose(fin);
  fclose(flist);

  sprintf(non,"\nCompleted."); strcat(txt,non);
  MessageBox(NULL,txt,"Surface",MB_OK);
  return 0;
}/*createyieldsurface*/






/*****DEBUGGING*****/

/**********/

char **fgetscut(FILE *fin,int *n)
/*ONE LINE FROM FILE BREAK INTO SEVERAL STRINGS.*/
/*SPACES INCLUDED.*/
/*RETURN:HEAD OF STRINGS.n=0 IF END OF FILE.*/
{
  FILE *f=fin;
  char **s; /*ADDRESS SET*/
  char *ss; /*STRING*/
  int ichr; /*CHARACTER*/
  int i,ii;

  *n=0;
  i=0;
  s=NULL;
  while((ichr=getc(f))==' ')
  ;
  while(1)
  {
    ii=0;
    ss=NULL;
    while(1)
    {
      if(ichr==EOF) return NULL;

      ss=(char *)realloc(ss,(ii+1)*sizeof(char));
      if(ss==NULL)
      {
        /*errormessage("FGETSBRK:MEMORY INSUFFICIENT.");*/
        return NULL;
      }
      if(ichr=='\n') {*(ss+ii)='\0'; break;}
      else if(ichr==' ') break;
      else {*(ss+ii)=(char)ichr;}

      ii++;
      ichr=getc(f);
    }
    while(ichr==' ')
    {
      *(ss+ii)=' ';
      ii++;
      ss=(char *)realloc(ss,(ii+1)*sizeof(char));
      if(ss==NULL) return NULL;
      *(ss+ii)='\0';

      ichr=getc(f);
    }

    s=(char **)realloc(s,(i+1)*sizeof(char *));
    if(s==NULL)
    {
      /*errormessage("FGETSBRK:MEMORY INSUFFICIENT.");*/
      return NULL;
    }
    *(s+i)=ss;
    i++;
    if(ichr=='\n') break;
  }
  *n=i;

  return s;
}/*fgetscut*/

int getsectionform(FILE *flist,int code,struct section *sect)
/*GET SECTION FROM SECTION LIST.*/
{
  char **data;
  int i,j,n,ns=0;
  int type;
  int nsteel=0,nrein=0,nconc=0,nstrnd=0,ncomment=0;

  fseek(flist,0L,SEEK_SET);

  sect->code=0;                                   /*INITIALIZATION.*/
  sect->soff=0;
  sect->stype=0;
  sect->etype=0;

  sect->nsteel=0;
  sect->nrein=0;
  sect->nconc=0;
  sect->nstrnd=0;
  sect->ncomment=0;

  sect->face[0][0]=0.0;
  sect->face[0][1]=0.0;
  sect->face[1][0]=0.0;
  sect->face[1][1]=0.0;

  for(i=0;i<MAXSTRND;i++)
  {
    sect->strnd[i].area[END]=0.0;
    sect->strnd[i].area[MID]=0.0;
    sect->strnd[i].Ni[END]=0.0;
    sect->strnd[i].Ni[MID]=0.0;
  }

  while(1)
  {
    /*data=fgetsbrk(flist,&n);*/
    data=fgetscut(flist,&n);
    if(n==0) return 0;

/*for(i=0;i<n;i++)
{
  sprintf(non,"/%s/",*(data+i));
  MessageBox(NULL,non,"Test",MB_OK);
}*/
    if(!strncmp(*(data+0),"CODE",4))
    {
      ns++;                                         /*OFFSET COUNT.*/
      sect->code=(int)strtol(*(data+1),NULL,10);

      if(sect->code==code)
      {
        if(!strncmp(*(data+2),"S ",2))   sect->stype=S;
        if(!strncmp(*(data+2),"RC",2))  sect->stype=RC;
        if(!strncmp(*(data+2),"SRC",3)) sect->stype=SRC;
        if(!strncmp(*(data+2),"PC",2))  sect->stype=PC;

        if(!strncmp(*(data+3),"COLUMN",6)) sect->etype=COLUMN;
        if(!strncmp(*(data+3),"GIRDER",6)) sect->etype=GIRDER;
        if(!strncmp(*(data+3),"BEAM",4))   sect->etype=BEAM;
        if(!strncmp(*(data+3),"BRACE",5))  sect->etype=BRACE;
        if(!strncmp(*(data+3),"WALL",4))   sect->etype=WALL;
        if(!strncmp(*(data+3),"SLAB",4))   sect->etype=SLAB;

        sprintf(sect->comment[ncomment],"CODE:%d %s%s",
                sect->code,*(data+2),*(data+3));
        ncomment++;

        for(i=0;i<n;i++)
        {
          if(!strncmp(*(data+i),"\"",1))
          {
            sprintf(sect->comment[ncomment],"\0");
            for(j=i;j<n;j++)
            {
              strcat(sect->comment[ncomment],*(data+j));
            }
            ncomment++;
            break;
          }
        }
        freestr(data,n);

        if(sect->etype==COLUMN ||
           sect->etype==GIRDER ||
           sect->etype==BEAM)
        {
          while(1)
          {
            /*data=fgetsbrk(flist,&n);*/
            data=fgetscut(flist,&n);
            if(n==0) break;
            if(!strncmp(*(data+0),"CODE",4)) break;

            if(!strncmp(*(data+0),"SRECT",5))
            {
              sect->srect[nsteel].left  =strtod(*(data+1),NULL);
              sect->srect[nsteel].bottom=strtod(*(data+2),NULL);
              sect->srect[nsteel].right =strtod(*(data+3),NULL);
              sect->srect[nsteel].top   =strtod(*(data+4),NULL);
              nsteel++;
            }
            if(!strncmp(*(data+0),"CRECT",4))
            {
              sect->crect[nconc].left  =strtod(*(data+1),NULL);
              sect->crect[nconc].bottom=strtod(*(data+2),NULL);
              sect->crect[nconc].right =strtod(*(data+3),NULL);
              sect->crect[nconc].top   =strtod(*(data+4),NULL);
              nconc++;
            }
            if(!strncmp(*(data+0),"REINS",5))
            {
              sect->rein[nrein].area=strtod(*(data+1),NULL);
              sect->rein[nrein].x=strtod(*(data+2),NULL);
              sect->rein[nrein].y=strtod(*(data+3),NULL);
              nrein++;
            }
            if(!strncmp(*(data+0),"HOOPS",5))
            {
              sect->shearrein[0]=strtod(*(data+1),NULL);
              sect->shearrein[1]=strtod(*(data+2),NULL);
            }
            if(!strncmp(*(data+0),"XFACE",5))
            {
              sect->face[SX][HEAD]=strtod(*(data+1),NULL);
              sect->face[SX][TAIL]=strtod(*(data+2),NULL);
            }
            if(!strncmp(*(data+0),"YFACE",5))
            {
              sect->face[SY][HEAD]=strtod(*(data+1),NULL);
              sect->face[SY][TAIL]=strtod(*(data+2),NULL);
            }
            if(!strncmp(*(data+0),"STRND",5))
            {
              if(!strncmp(*(data+2),"A",1))      type=STRNDA;
              else if(!strncmp(*(data+2),"B",1)) type=STRNDB;
              else if(!strncmp(*(data+2),"C",1)) type=STRNDC;
              else if(!strncmp(*(data+2),"W",1)) type=STRNDW;
              else                               type=4;

              sect->strnd[nstrnd].type=type;

              if(!strncmp(*(data+1),"ALL",3) ||
                 !strncmp(*(data+1),"END",3))
              {
                sect->strnd[nstrnd].area[END]=strtod(*(data+3),NULL);
                sect->strnd[nstrnd].x[END]=strtod(*(data+4),NULL);
                sect->strnd[nstrnd].y[END]=strtod(*(data+5),NULL);
                sect->strnd[nstrnd].Ni[END]=strtod(*(data+6),NULL);
              }
              if(!strncmp(*(data+1),"ALL",3) ||
                 !strncmp(*(data+1),"MID",3))
              {
                sect->strnd[nstrnd].area[MID]=strtod(*(data+3),NULL);
                sect->strnd[nstrnd].x[MID]=strtod(*(data+4),NULL);
                sect->strnd[nstrnd].y[MID]=strtod(*(data+5),NULL);
                sect->strnd[nstrnd].Ni[MID]=strtod(*(data+6),NULL);
              }
              nstrnd++;
            }

            for(i=0;i<n;i++)
            {
              if(!strncmp(*(data+i),"\"",1))
              {
                sprintf(sect->comment[ncomment],"\0");
                for(j=i;j<n;j++)
                {
                  strcat(sect->comment[ncomment],*(data+j));
                }
                ncomment++;
                break;
              }
            }
            freestr(data,n);
          }
        }
        if(sect->etype==WALL || sect->etype==SLAB)
        {
          while(1)
          {
            /*data=fgetsbrk(flist,&n);*/
            data=fgetscut(flist,&n);
            if(n==0) break;
            if(!strncmp(*(data+0),"CODE",4)) break;

            if(!strncmp(*(data+0),"THICK",5))
            {
              sect->thick=strtod(*(data+1),NULL);
            }
            if(!strncmp(*(data+0),"SREIN",5))
            {
              sect->shearrein[0]=strtod(*(data+1),NULL);
            }
            if(!strncmp(*(data+0),"WRECT",5))
            {
              sect->wlength=strtod(*(data+1),NULL);
              sect->wheight=strtod(*(data+2),NULL);
            }
            if(!strncmp(*(data+0),"XFACE",5)) /*FACE FOR LENGTH.*/
            {
              sect->face[0][HEAD]=strtod(*(data+1),NULL);
              sect->face[0][TAIL]=strtod(*(data+2),NULL);
            }
            if(!strncmp(*(data+0),"YFACE",5)) /*FACE FOR HEIGHT.*/
            {
              sect->face[1][HEAD]=strtod(*(data+1),NULL);
              sect->face[1][TAIL]=strtod(*(data+2),NULL);
            }

            for(i=0;i<n;i++)
            {
              if(!strncmp(*(data+i),"\"",1))
              {
                sprintf(sect->comment[ncomment],"\0");
                for(j=i;j<n;j++)
                {
                  strcat(sect->comment[ncomment],*(data+j));
                }
                ncomment++;
                break;
              }
            }
            freestr(data,n);
          }
        }

        sect->soff=ns;
        sect->nsteel=nsteel;
        sect->nrein=nrein;
        sect->nconc=nconc;
        sect->nstrnd=nstrnd;
        sect->ncomment=ncomment;

        return ns;
      }
    }
    freestr(data,n);
  }
}/*getsectionform*/

void drawsectionlist(HDC hdc,int nsect,struct section *slist,
                     long int ox,long int oy,
                     long int pagewidth,long int pageheight,
                     struct viewparam vp,int mode)
{
  int i;
  int commentwidth,commentheight;
  long int top,height;
  double blank=20;
  double Xmin,Xmax,Ymin,Ymax;

  if(oy>pageheight) return;
  top=oy;

  if(mode==ONPRINTER) StartPage(hdc);

  for(i=0;i<nsect;i++)
  {
    getsectionspace(hdc,(slist+i),
                    &Xmin,&Xmax,&Ymin,&Ymax,
                    &commentwidth,&commentheight);

    height=(long int)(vp.gfactor*(Ymax-Ymin));
    if(height<commentheight) height=commentheight;

    if(fabs(Xmax-Xmin)>150.0)
    {
      height=(long int)(vp.gfactor*(Ymax-Ymin))
            +commentheight
            +10;
    }

    if(mode==ONSCREEN)
    {
      if((top<pageheight && (top+height)>0) ||
         (top<0 && (top+height)>0))
      {
        drawsection(hdc,(slist+i),
                    255,0,0,
                    255,255,0,
                    0,255,255,
                    0,255,0,
                    255,255,255,
                    1,
                    ox,top,
                    vp);
      }
    }
    if(mode==ONPRINTER)
    {
      if((top+height)>pageheight)
      {
        top=oy;

        EndPage(hdc);
        StartPage(hdc);
      }

      drawsection(hdc,(slist+i),
                  0,0,0,
                  0,0,0,
                  0,0,0,
                  0,0,0,
                  0,0,0,
                  3,
                  ox,top,
                  vp);
    }

    top+=height+(long int)(vp.gfactor*blank);
  }

  if(mode==ONPRINTER) EndPage(hdc);

  return;
}/*drawsectionlist*/

void getsectionspace(HDC hdc,
                     struct section *sect,
                     double *Xmin,double *Xmax,
                     double *Ymin,double *Ymax,
                     int *commentwidth,int *commentheight)
{
  SIZE size;
  int i;

  *Xmin=0.0; *Xmax=0.0;
  *Ymin=0.0; *Ymax=0.0;

  *commentwidth =0;
  *commentheight=0;

  for(i=0;i<sect->nsteel;i++)
  {
    if(*Xmin > sect->srect[i].left)  *Xmin=sect->srect[i].left;
    if(*Xmin > sect->srect[i].right) *Xmin=sect->srect[i].right;
    if(*Xmax < sect->srect[i].left)  *Xmax=sect->srect[i].left;
    if(*Xmax < sect->srect[i].right) *Xmax=sect->srect[i].right;

    if(*Ymin > sect->srect[i].top)    *Ymin=sect->srect[i].top;
    if(*Ymin > sect->srect[i].bottom) *Ymin=sect->srect[i].bottom;
    if(*Ymax < sect->srect[i].top)    *Ymax=sect->srect[i].top;
    if(*Ymax < sect->srect[i].bottom) *Ymax=sect->srect[i].bottom;
  }

  for(i=0;i<sect->nrein;i++)
  {
    if(*Xmin > sect->rein[i].x) *Xmin=sect->rein[i].x;
    if(*Xmax < sect->rein[i].x) *Xmax=sect->rein[i].x;

    if(*Ymin > sect->rein[i].y) *Ymin=sect->rein[i].y;
    if(*Ymax < sect->rein[i].y) *Ymax=sect->rein[i].y;
  }

  for(i=0;i<sect->nconc;i++)
  {
    if(*Xmin > sect->crect[i].left)  *Xmin=sect->crect[i].left;
    if(*Xmin > sect->crect[i].right) *Xmin=sect->crect[i].right;
    if(*Xmax < sect->crect[i].left)  *Xmax=sect->crect[i].left;
    if(*Xmax < sect->crect[i].right) *Xmax=sect->crect[i].right;

    if(*Ymin > sect->crect[i].top)    *Ymin=sect->crect[i].top;
    if(*Ymin > sect->crect[i].bottom) *Ymin=sect->crect[i].bottom;
    if(*Ymax < sect->crect[i].top)    *Ymax=sect->crect[i].top;
    if(*Ymax < sect->crect[i].bottom) *Ymax=sect->crect[i].bottom;
  }

  for(i=0;i<sect->nstrnd;i++)
  {
    if(*Xmin > sect->strnd[i].x[0]) *Xmin=sect->strnd[i].x[0];
    if(*Xmax < sect->strnd[i].x[0]) *Xmax=sect->strnd[i].x[0];

    if(*Ymin > sect->strnd[i].y[0]) *Ymin=sect->strnd[i].y[0];
    if(*Ymax < sect->strnd[i].y[0]) *Ymax=sect->strnd[i].y[0];

    if(*Xmin > sect->strnd[i].x[1]) *Xmin=sect->strnd[i].x[1];
    if(*Xmax < sect->strnd[i].x[1]) *Xmax=sect->strnd[i].x[1];

    if(*Ymin > sect->strnd[i].y[1]) *Ymin=sect->strnd[i].y[1];
    if(*Ymax < sect->strnd[i].y[1]) *Ymax=sect->strnd[i].y[1];
  }

  for(i=0;i<sect->ncomment;i++)
  {
    GetTextExtentPoint32(hdc,
                         sect->comment[i],
                         strlen(sect->comment[i]),&size);

    if(*commentwidth < size.cx) *commentwidth=size.cx;
    *commentheight+=size.cy;
  }

  return;
}/*getsectionspace*/

void drawsection(HDC hdc,struct section *sect,
                 int sr,int sg,int sb, /*STEEL COLOR*/
                 int rr,int rg,int rb, /*REIN COLOR*/
                 int cr,int cg,int cb, /*CONCRETE COLOR*/
                 int pr,int pg,int pb, /*STRAND COLOR*/
                 int tr,int tg,int tb, /*TEXT COLOR*/
                 int lw,               /*LINE WIDTH*/
                 long int ox,long int oy,
                 struct viewparam vp)
{
  SIZE size;
  HPEN hpenS,hpenR,hpenC,hpenP,ppen;
  int i;
  int commentwidth,commentheight,commentleft;
  double textleft=200.0;
  double Xmin,Xmax,Ymin,Ymax;

  commentleft=ox+(int)(vp.gfactor*textleft);

  getsectionspace(hdc,sect,
                  &Xmin,&Xmax,&Ymin,&Ymax,
                  &commentwidth,&commentheight);

  ox-=(long int)(vp.gfactor*Xmin);
  oy+=(long int)(vp.gfactor*Ymax);

  hpenS=CreatePen(PS_SOLID,lw,RGB(sr,sg,sb)); /*STEEL*/
  hpenR=CreatePen(PS_SOLID,lw,RGB(rr,rg,rb)); /*REIN*/
  hpenC=CreatePen(PS_SOLID,lw,RGB(cr,cg,cb)); /*CONCRETE*/
  hpenP=CreatePen(PS_SOLID,lw,RGB(pr,pg,pb)); /*STRAND*/

  ppen=SelectObject(hdc,hpenS);
  for(i=0;i<sect->nsteel;i++)
  {
    drawmaterialrect(hdc,&(sect->srect[i]),ox,oy,vp);
  }

  SelectObject(hdc,hpenR);
  for(i=0;i<sect->nrein;i++)
  {
    drawreincircle(hdc,&(sect->rein[i]),ox,oy,vp);
  }

  SelectObject(hdc,hpenC);
  for(i=0;i<sect->nconc;i++)
  {
    drawmaterialrect(hdc,&(sect->crect[i]),ox,oy,vp);
  }

  SelectObject(hdc,hpenP);
  for(i=0;i<sect->nstrnd;i++)
  {
    drawstrandcircle(hdc,&(sect->strnd[i]),END,ox,oy,vp);
  }

  SelectObject(hdc,ppen);
  DeleteObject(hpenS);
  DeleteObject(hpenR);
  DeleteObject(hpenC);
  DeleteObject(hpenP);

  if(fabs(Xmax-Xmin)>textleft)
  {
    oy+=(long int)(vp.gfactor*(Ymax-Ymin))+10;
  }

  SetTextColor(hdc,RGB(tr,tg,tb));
  for(i=0;i<sect->ncomment;i++)
  {
    GetTextExtentPoint32(hdc,
                         sect->comment[i],
                         strlen(sect->comment[i]),&size);
    TextOut(hdc,
            commentleft,
            oy-(int)(vp.gfactor*Ymax)+i*(size.cy),
            sect->comment[i],strlen(sect->comment[i]));
  }

  return;
}/*drawsection*/

void drawmaterialrect(HDC hdc,
                      struct materialrect *mr,
                      long int ox,long int oy,
                      struct viewparam vp)
{
  MoveToEx(hdc,ox+(int)(vp.gfactor*mr->left),
               oy-(int)(vp.gfactor*mr->top) ,NULL);
  LineTo(hdc,ox+(int)(vp.gfactor*mr->left),
             oy-(int)(vp.gfactor*mr->bottom));
  LineTo(hdc,ox+(int)(vp.gfactor*mr->right),
             oy-(int)(vp.gfactor*mr->bottom));
  LineTo(hdc,ox+(int)(vp.gfactor*mr->right),
             oy-(int)(vp.gfactor*mr->top));
  LineTo(hdc,ox+(int)(vp.gfactor*mr->left),
             oy-(int)(vp.gfactor*mr->top));

  return;
}/*drawmaterialrect*/

void drawreincircle(HDC hdc,
                    struct reinforcement *r,
                    long int ox,long int oy,
                    struct viewparam vp)
{
  double radius,dradius;
  int iradius;

  radius=sqrt(fabs(r->area)/PI);

  dradius=vp.gfactor*radius;
  iradius=(int)dradius;

  if(dradius<1.0)
  {
    MoveToEx(hdc,
             (int)ox+(int)(vp.gfactor*(r->x)),
             (int)oy-(int)(vp.gfactor*(r->y)),
             NULL);
    LineTo(hdc,
           (int)ox+(int)(vp.gfactor*(r->x)),
           (int)oy-(int)(vp.gfactor*(r->y))+1);
  }
  else if(dradius<2.0)
  {
    Arc(hdc,
        (int)ox+(int)(vp.gfactor*(r->x))-1,
        (int)oy-(int)(vp.gfactor*(r->y))-1,
        (int)ox+(int)(vp.gfactor*(r->x))+2,
        (int)oy-(int)(vp.gfactor*(r->y))+2,
        0,0,0,0); /*HOLLOW CIRCLE.*/
  }
  else
  {
    Arc(hdc,
        (int)ox+(int)(vp.gfactor*(r->x))-iradius,
        (int)oy-(int)(vp.gfactor*(r->y))-iradius,
        (int)ox+(int)(vp.gfactor*(r->x))+iradius,
        (int)oy-(int)(vp.gfactor*(r->y))+iradius,
        0,0,0,0); /*HOLLOW CIRCLE.*/
  }
  return;
}/*drawreincircle*/

void drawstrandcircle(HDC hdc,
                      struct pcstrand *ps,int position,
                      long int ox,long int oy,
                      struct viewparam vp)
{
  double radius,dradius;
  int iradius1,iradius2;

  radius=sqrt((ps->area[position])/PI);

  dradius=vp.gfactor*radius;

  iradius1=(int)dradius;
  iradius2=(int)dradius*0.7;

  if(dradius<1.0)
  {
    MoveToEx(hdc,
             (int)ox+(int)(vp.gfactor*(ps->x[position])),
             (int)oy-(int)(vp.gfactor*(ps->y[position])),
             NULL);
    LineTo(hdc,
           (int)ox+(int)(vp.gfactor*(ps->x[position])),
           (int)oy-(int)(vp.gfactor*(ps->y[position]))+1);
  }
  else if(dradius<2.0)
  {
    Arc(hdc,
        (int)ox+(int)(vp.gfactor*(ps->x[position]))-1,
        (int)oy-(int)(vp.gfactor*(ps->y[position]))-1,
        (int)ox+(int)(vp.gfactor*(ps->x[position]))+2,
        (int)oy-(int)(vp.gfactor*(ps->y[position]))+2,
        0,0,0,0); /*HOLLOW CIRCLE.*/
  }
  else
  {
    Arc(hdc,
        (int)ox+(int)(vp.gfactor*(ps->x[position]))-iradius1,
        (int)oy-(int)(vp.gfactor*(ps->y[position]))-iradius1,
        (int)ox+(int)(vp.gfactor*(ps->x[position]))+iradius1,
        (int)oy-(int)(vp.gfactor*(ps->y[position]))+iradius1,
        0,0,0,0); /*HOLLOW CIRCLE.*/
  }

  Arc(hdc,
      ox+(int)(vp.gfactor*(ps->x[position]))-iradius2,
      oy-(int)(vp.gfactor*(ps->y[position]))-iradius2,
      ox+(int)(vp.gfactor*(ps->x[position]))+iradius2,
      oy-(int)(vp.gfactor*(ps->y[position]))+iradius2,
      0,0,0,0); /*HOLLOW CIRCLE.*/

  MoveToEx(hdc,ox+(int)(vp.gfactor*(ps->x[position])),
               oy-(int)(vp.gfactor*(ps->y[position]+1.5*radius)),
               NULL);
  LineTo(hdc,ox+(int)(vp.gfactor*(ps->x[position])),
             oy-(int)(vp.gfactor*(ps->y[position]-1.5*radius)));

  MoveToEx(hdc,ox+(int)(vp.gfactor*(ps->x[position]-1.5*radius)),
               oy-(int)(vp.gfactor*(ps->y[position])),
               NULL);
  LineTo(hdc,ox+(int)(vp.gfactor*(ps->x[position]+1.5*radius)),
             oy-(int)(vp.gfactor*(ps->y[position])));
  return;
}/*drawstrandcircle*/

/*********/
/*SRCAM001.H:SRC SUBROUTINES SINCE 1995.12.21.JUNSATO.*/
/*HEADER FOR SRCAN001.C.*/
/*LAST CHANGE:1997.11.17.*/

/*S,RC REMAINING DIFFERENT PLAIN WITH SAME NEUTRAL AXIS.*/

/*S COLUMN...................................................*/
/*ALLOWABLE TENSION OF S RECTANGLES................COMPLETED.*/
/*ALLOWABLE COMPRESSION OF S RECTANGLES............COMPLETED.*/
/*ALLOWABLE BENDING OF S COLUMN WEAK...............COMPLETED.*/
/*ALLOWABLE BENDING OF S COLUMN STRONG..........ONLY FOR H,[.*/
/*ALLOWABLE SHEAR OF S RECTANGLES...............ONLY FOR H,[.*/
/*PLASTIC BENDING OF S COLUMN......................COMPLETED.*/
/*ULTIMATE BENDING OF S COLUMN......................*/
/*ULTIMATE SHEAR OF S RECTANGLES................ONLY FOR H,[.*/

/*RC COLUMN..................................................*/
/*ALLOWABLE BENDING OF RC COLUMN...................COMPLETED.*/
/*ALLOWABLE SHEAR OF RC COLUMN,GIRDER FOR LONG.....COMPLETED.*/
/*ALLOWABLE SHEAR OF RC COLUMN,GIRDER FOR SHORT....COMPLETED.*/
/*ULTIMATE BENDING OF RC COLUMN....................COMPLETED.*/
/*ULTIMATE SHEAR OF RC COLUMN......................COMPLETED.*/

/*RC WALL....................................................*/
/*ALLOWABLE SHEAR OF RC WALL.......................COMPLETED.*/
/*ULTIMATE SHEAR OF RC WALL.........................*/

/*SRC COLUMN.................................................*/
/*ALLOWABLE BENDING OF SRC COLUMN..................COMPLETED.*/
/*ALLOWABLE SHEAR OF SRC COLUMN......ONLY ONE CRECT WITH H,[.*/
/*ULTIMATE BENDING OF SRC COLUMN...................COMPLETED.*/
/*ULTIMATE SHEAR OF SRC COLUMN.......ONLY ONE CRECT WITH H,[.*/

/*PC COLUMN..................................................*/
/*ALLOWABLE BENDING OF PC COLUMN..............ONLY ONE CRECT.*/
/*ALLOWABLE SHEAR OF PC COLUMN................ONLY ONE CRECT.*/
/*ULTIMATE BENDING OF PC COLUMN...............ONLY ONE CRECT.*/
/*ULTIMATE SHEAR OF PC COLUMN.................ONLY ONE CRECT.*/

/*PC GIRDER..................................................*/
/*ALLOWABLE BENDING OF PC GIRDER..............ONLY ONE CRECT.*/
/*ALLOWABLE SHEAR OF PC GIRDER................SAME AS COLUMN.*/
/*ULTIMATE BENDING OF PC GIRDER...............ONLY ONE CRECT.*/
/*ULTIMATE SHEAR OF PC GIRDER.................ONLY ONE CRECT.*/

/*FOR ROUTE 2.3.*/
/*SECTION FROM SECTION LIST FILE.*/

void comments(FILE *f)
{
  fprintf(f,"\n");
  fprintf(f,"As:鉄骨 Ar:主筋 Ac:コンクリート Ap:ＰＣストランド\n");
  fprintf(f,"N:軸力 Q:せん断力 ");
  fprintf(f,"Mt:ねじりモーメント M:曲げモーメント\n");
  fprintf(f,"添字 i:始端 j:終端 c:中央\n");
  fprintf(f,"     a:許容 u:終局\n");
  return;
}/*comments*/

void inputfiletomemory(FILE *ftext,
                       int *nnode,int *nelem,int *nsect,
                       struct section *sects,
                       struct structnode *nodes,
                       struct element *elems)
/*TRANSLATE FRAME3 INPUTFILE INTO MEMORY.*/
{
  char **data;
  int i,ii,n,count;
  long int hcode,tcode,scode;
  double E,G;

  fseek(ftext,0L,SEEK_SET);
  fgetinitial(ftext,nnode,nelem,nsect,&E,&G);      /*INPUT INITIAL.*/

  /*"INPUTFILE SPECIFICATION"*/
  /*NNODE NELEM MSIZE BAND 0 E G 1.0E-04 1.0E+8 0 1 0 NSECT*/
  /*ISECT TYPE A Ixx Iyy J*/
  /*INODE X Y Z*/
  /*IELEM NODEI J SECT BOUNDARY... COORD EFACT GFACT CMQTYPE CMQ...*/
  /*INODE CONFINEMENTVALUE... CONFINEMENTTYPE...*/

  for(i=0;i<*nsect;i++)
  {
    data=fgetsbrk(ftext,&n);

    (sects+i)->code=strtol(*(data+0),NULL,10);

    freestr(data,n);
  }
  for(i=0;i<*nnode;i++)
  {
    data=fgetsbrk(ftext,&n);

    (nodes+i)->code=strtol(*(data+0),NULL,10);
    (nodes+i)->x=strtod(*(data+1),NULL);
    (nodes+i)->y=strtod(*(data+2),NULL);
    (nodes+i)->z=strtod(*(data+3),NULL);

    freestr(data,n);
  }
  for(i=0;i<*nelem;i++)
  {
    data=fgetsbrk(ftext,&n);

    (elems+i)->code=strtol(*(data+0),NULL,10);

    hcode=strtol(*(data+1),NULL,10); /*HEAD*/
    tcode=strtol(*(data+2),NULL,10); /*TAIL*/

    scode=strtol(*(data+3),NULL,10); /*SECTION*/

    (elems+i)->cmqcode=strtol(*(data+13),NULL,10); /*CMQ*/

    freestr(data,n);

    count=0;
    for(ii=0;(ii<(*nnode) && count<2);ii++)
    {
      if((nodes+ii)->code==hcode)
      {
        (elems+i)->node[0]=(nodes+ii); /*POINT HEAD.*/
        count++;
      }
      if((nodes+ii)->code==tcode)
      {
        (elems+i)->node[1]=(nodes+ii); /*POINT TAIL.*/
        count++;
      }
    }
    for(ii=0;ii<(*nsect);ii++) /*POINT SECTION.*/
    {
      if((sects+ii)->code==scode)
      {
        (elems+i)->sect=(sects+ii);
        break;
      }
    }
  }

  return;
}/*inputfiletomemory*/

void inputfiletomemory2(FILE *ftext,
                       int *nnode,int *nelem,int *nsect,
                       struct section *sects,
                       struct structnode *nodes,
                       struct element *elems)
/*TRANSLATE ARCLM INPUTFILE INTO MEMORY.*/
{
  char **data;
  int i,ii,n,count;
  long int hcode,tcode,scode;
  /*double E,G;*/

  fseek(ftext,0L,SEEK_SET);
  fgetinitial2(ftext,nnode,nelem,nsect);      /*INPUT INITIAL.*/

  /*"INPUTFILE SPECIFICATION"*/
  /*NNODE NELEM MSIZE BAND 0 E G 1.0E-04 1.0E+8 0 1 0 NSECT*/
  /*ISECT TYPE A Ixx Iyy J*/
  /*INODE X Y Z*/
  /*IELEM NODEI J SECT BOUNDARY... COORD EFACT GFACT CMQTYPE CMQ...*/
  /*INODE CONFINEMENTVALUE... CONFINEMENTTYPE...*/

  for(i=0;i<*nsect;i++)
  {
    data=fgetsbrk(ftext,&n);

    (sects+i)->code=strtol(*(data+0),NULL,10);

    freestr(data,n);
  }
  for(i=0;i<*nnode;i++)
  {
    data=fgetsbrk(ftext,&n);

    (nodes+i)->code=strtol(*(data+0),NULL,10);
    (nodes+i)->x=strtod(*(data+1),NULL);
    (nodes+i)->y=strtod(*(data+2),NULL);
    (nodes+i)->z=strtod(*(data+3),NULL);

    freestr(data,n);
  }
  for(i=0;i<*nelem;i++)
  {
    data=fgetsbrk(ftext,&n);

    (elems+i)->code=strtol(*(data+0),NULL,10);

    hcode=strtol(*(data+2),NULL,10); /*HEAD*/
    tcode=strtol(*(data+3),NULL,10); /*TAIL*/

    scode=strtol(*(data+1),NULL,10); /*SECTION*/

    /*(elems+i)->cmqcode=strtol(*(data+13),NULL,10);*/ /*CMQ*/

    freestr(data,n);

    count=0;
    for(ii=0;(ii<(*nnode) && count<2);ii++)
    {
      if((nodes+ii)->code==hcode)
      {
        (elems+i)->node[0]=(nodes+ii); /*POINT HEAD.*/
        count++;
      }
      if((nodes+ii)->code==tcode)
      {
        (elems+i)->node[1]=(nodes+ii); /*POINT TAIL.*/
        count++;
      }
    }
    for(ii=0;ii<(*nsect);ii++) /*POINT SECTION.*/
    {
      if((sects+ii)->code==scode)
      {
        (elems+i)->sect=(sects+ii);
        break;
      }
    }
  }

  return;
}/*inputfiletomemory2*/

void fgetinitial(FILE *fin,int *nnode,int *nelem,int *nsect,
                 double *E,double *G)
/*GET INITIAL DATA FROM FRAME3 INPUTFILE.*/
{
  char **data;
  int n;

  data=fgetsbrk(fin,&n);

  *nnode=strtol(*(data+0),NULL,10);
  *nelem=strtol(*(data+1),NULL,10);
  *E=strtod(*(data+5),NULL);
  *G=strtod(*(data+6),NULL);
  *nsect=strtol(*(data+12),NULL,10);

  freestr(data,n);

  return;
}/*fgetinitial*/

void fgetinitial2(FILE *fin,int *nnode,int *nelem,int *nsect)
/*GET INITIAL DATA FROM ARCLM INPUTFILE.*/
{
  char **data;
  int n;

  data=fgetsbrk(fin,&n);

  *nnode=strtol(*(data+0),NULL,10);
  *nelem=strtol(*(data+1),NULL,10);
  *nsect=strtol(*(data+2),NULL,10);

  freestr(data,n);

  return;
}/*fgetinitial2*/

void initializestress(struct stress *s)
/*INITIALIZE STRESSES.*/
{
  s->N=0.0;
  s->Q[SX]=0.0;
  s->Q[SY]=0.0;
  s->Mt=0.0;
  s->M[SX]=0.0;
  s->M[SY]=0.0;

  return;
}/*initializestress*/

int fgetstress(FILE *fin,int *ielem,int *isect,int *inode,
               struct stress *s)
/*6 STRESSES FROM OUTPUT FILE.*/
{
  char **data;
  int i,nstr;

  data=fgetsbrk(fin,&nstr);
  if(nstr==9) /*HEAD*/
  {
    *ielem=(int)strtol(*(data+0),NULL,10);
    *isect=(int)strtol(*(data+1),NULL,10);
    *inode=(int)strtol(*(data+2),NULL,10);
    s->N=strtod(*(data+3),NULL);
    s->Q[0]=strtod(*(data+4),NULL);
    s->Q[1]=strtod(*(data+5),NULL);
    s->Mt=strtod(*(data+6),NULL);
    s->M[0]=strtod(*(data+7),NULL);
    s->M[1]=strtod(*(data+8),NULL);
  }
  else if(nstr==7) /*TAIL*/
  {
    *inode=(int)strtol(*(data+0),NULL,10);
    s->Q[0]=strtod(*(data+2),NULL);
    s->Q[1]=strtod(*(data+3),NULL);
    s->M[0]=strtod(*(data+5),NULL);
    s->M[1]=strtod(*(data+6),NULL);
  }
  else if(nstr==8) /*FOR KUROSAWA TAIL*/
  {
    *inode=(int)strtol(*(data+1),NULL,10);
    s->Q[0]=strtod(*(data+3),NULL);
    s->Q[1]=strtod(*(data+4),NULL);
    s->M[0]=strtod(*(data+6),NULL);
    s->M[1]=strtod(*(data+7),NULL);
  }
  for(i=1;i<=nstr;i++) free(*(data+i-1));
  free(data);

  return nstr;
}/*fgetstress*/

int getsrcansection(FILE *flist,long int code,struct section *sect)
/*GET SECTION FROM SECTION LIST.*/
{
  char **data;
  int i,n,ns=0;
  int type;
  int nsteel=0,nrein=0,nconc=0,nstrnd=0;

  fseek(flist,0L,SEEK_SET);

  sect->code=0;                                   /*INITIALIZATION.*/
  sect->soff=0;
  sect->stype=0;
  sect->etype=0;

  sect->nsteel=0;
  sect->nrein=0;
  sect->nconc=0;
  sect->nstrnd=0;

  sect->face[0][0]=0.0;
  sect->face[0][1]=0.0;
  sect->face[1][0]=0.0;
  sect->face[1][1]=0.0;

  for(i=0;i<MAXSTRND;i++)
  {
    sect->strnd[i].area[END]=0.0;
    sect->strnd[i].area[MID]=0.0;
    sect->strnd[i].Ni[END]=0.0;
    sect->strnd[i].Ni[MID]=0.0;
  }

  while(1)
  {
    data=fgetsbrk(flist,&n);
    if(n==0) return 0;

    if(!strcmp(*(data+0),"CODE"))
    {
      ns++;                                         /*OFFSET COUNT.*/
      sect->code=strtol(*(data+1),NULL,10);

      if(sect->code==code)
      {
        if(!strcmp(*(data+2),"S"))   sect->stype=S;
        if(!strcmp(*(data+2),"RC"))  sect->stype=RC;
        if(!strcmp(*(data+2),"SRC")) sect->stype=SRC;
        if(!strcmp(*(data+2),"PC"))  sect->stype=PC;

        if(!strcmp(*(data+3),"COLUMN")) sect->etype=COLUMN;
        if(!strcmp(*(data+3),"GIRDER")) sect->etype=GIRDER;
        if(!strcmp(*(data+3),"BEAM"))   sect->etype=BEAM;
        if(!strcmp(*(data+3),"BRACE"))  sect->etype=BRACE;
        if(!strcmp(*(data+3),"WALL"))   sect->etype=WALL;
        if(!strcmp(*(data+3),"SLAB"))   sect->etype=SLAB;

        freestr(data,n);

        if(sect->etype==COLUMN ||
           sect->etype==GIRDER ||
           sect->etype==BEAM)
        {
          while(1)
          {
            data=fgetsbrk(flist,&n);
            if(n==0) break;
            if(!strcmp(*(data+0),"CODE")) break;

            if(!strcmp(*(data+0),"SRECT"))
            {
              sect->srect[nsteel].left  =strtod(*(data+1),NULL);
              sect->srect[nsteel].bottom=strtod(*(data+2),NULL);
              sect->srect[nsteel].right =strtod(*(data+3),NULL);
              sect->srect[nsteel].top   =strtod(*(data+4),NULL);
              nsteel++;
            }
            if(!strcmp(*(data+0),"CRECT"))
            {
              sect->crect[nconc].left  =strtod(*(data+1),NULL);
              sect->crect[nconc].bottom=strtod(*(data+2),NULL);
              sect->crect[nconc].right =strtod(*(data+3),NULL);
              sect->crect[nconc].top   =strtod(*(data+4),NULL);
              nconc++;
            }
            if(!strcmp(*(data+0),"REINS"))
            {
              sect->rein[nrein].area=strtod(*(data+1),NULL);
              sect->rein[nrein].x=strtod(*(data+2),NULL);
              sect->rein[nrein].y=strtod(*(data+3),NULL);
              nrein++;
            }
            if(!strcmp(*(data+0),"HOOPS"))
            {
              sect->shearrein[0]=strtod(*(data+1),NULL);
              sect->shearrein[1]=strtod(*(data+2),NULL);
            }
            if(!strcmp(*(data+0),"XFACE"))
            {
              sect->face[SX][HEAD]=strtod(*(data+1),NULL);
              sect->face[SX][TAIL]=strtod(*(data+2),NULL);
            }
            if(!strcmp(*(data+0),"YFACE"))
            {
              sect->face[SY][HEAD]=strtod(*(data+1),NULL);
              sect->face[SY][TAIL]=strtod(*(data+2),NULL);
            }
            if(!strcmp(*(data+0),"STRND"))
            {
              if(!strcmp(*(data+2),"A"))      type=STRNDA;
              else if(!strcmp(*(data+2),"B")) type=STRNDB;
              else if(!strcmp(*(data+2),"C")) type=STRNDC;
              else if(!strcmp(*(data+2),"W")) type=STRNDW;
              else                            type=4;

              sect->strnd[nstrnd].type=type;

              if(!strcmp(*(data+1),"ALL") ||
                 !strcmp(*(data+1),"END"))
              {
                sect->strnd[nstrnd].area[END]=strtod(*(data+3),NULL);
                sect->strnd[nstrnd].x[END]=strtod(*(data+4),NULL);
                sect->strnd[nstrnd].y[END]=strtod(*(data+5),NULL);
                sect->strnd[nstrnd].Ni[END]=strtod(*(data+6),NULL);
              }
              if(!strcmp(*(data+1),"ALL") ||
                 !strcmp(*(data+1),"MID"))
              {
                sect->strnd[nstrnd].area[MID]=strtod(*(data+3),NULL);
                sect->strnd[nstrnd].x[MID]=strtod(*(data+4),NULL);
                sect->strnd[nstrnd].y[MID]=strtod(*(data+5),NULL);
                sect->strnd[nstrnd].Ni[MID]=strtod(*(data+6),NULL);
              }
              nstrnd++;
            }
            freestr(data,n);
          }
        }
        if(sect->etype==WALL || sect->etype==SLAB)
        {
          while(1)
          {
            data=fgetsbrk(flist,&n);
            if(n==0) break;
            if(!strcmp(*(data+0),"CODE")) break;

            if(!strcmp(*(data+0),"THICK"))
            {
              sect->thick=strtod(*(data+1),NULL);
            }
            if(!strcmp(*(data+0),"SREIN"))
            {
              sect->shearrein[0]=strtod(*(data+1),NULL);
            }
            if(!strcmp(*(data+0),"WRECT"))
            {
              sect->wlength=strtod(*(data+1),NULL);
              sect->wheight=strtod(*(data+2),NULL);
            }
            if(!strcmp(*(data+0),"XFACE")) /*FACE FOR LENGTH.*/
            {
              sect->face[0][HEAD]=strtod(*(data+1),NULL);
              sect->face[0][TAIL]=strtod(*(data+2),NULL);
            }
            if(!strcmp(*(data+0),"YFACE")) /*FACE FOR HEIGHT.*/
            {
              sect->face[1][HEAD]=strtod(*(data+1),NULL);
              sect->face[1][TAIL]=strtod(*(data+2),NULL);
            }

            freestr(data,n);
          }
        }

        sect->soff=ns;
        sect->nsteel=nsteel;
        sect->nrein=nrein;
        sect->nconc=nconc;
        sect->nstrnd=nstrnd;

        return ns;
      }
    }
    freestr(data,n);
  }
}/*getsrcansection*/

int getcodelist(FILE *flist,long int codelist[])
/*GET SECTION CODE LIST FROM SECTION LIST.*/
{
  char **data;
  int n,ns=0;

  fseek(flist,0L,SEEK_SET);

  while(1)
  {
    data=fgetsbrk(flist,&n);

    if(n==0) return ns;

    if(!strcmp(*(data+0),"CODE"))
    {
      codelist[ns]=strtol(*(data+1),NULL,10);
      ns++;
    }
    freestr(data,n);
  }
}/*getcodelist*/

int getcmqlist(FILE *fcmq,int cmqcode[],double Mo[])
/*GET CMQ FROM CMQ LIST.*/
{
  char **data;
  int n,ns=0;

  fseek(fcmq,0L,SEEK_SET);

  while(1)
  {
    data=fgetsbrk(fcmq,&n);

    if(n==0) return ns;

    if(!strcmp(*(data+0),"CMQ"))
    {
      cmqcode[ns]=(int)strtol(*(data+1),NULL,10);
      Mo[ns]=strtod(*(data+4),NULL);
      ns++;
    }
    freestr(data,n);
  }
}/*getcmqlist*/

double elementlength(struct element elem)
/*LENGTH OF ELEMENT.*/
{
  double dx,dy,dz,length;

  dx=(elem.node[1]->x)-(elem.node[0]->x);
  dy=(elem.node[1]->y)-(elem.node[0]->y);
  dz=(elem.node[1]->z)-(elem.node[0]->z);

  length=sqrt(dx*dx+dy*dy+dz*dz);

  return length;
}/*elementlength*/

double walllength(struct element elem)
/*LENGTH OF WALL ELEMENT.*/
{
  double dx,dy,length;

  dx=(elem.node[1]->x)-(elem.node[0]->x);
  dy=(elem.node[1]->y)-(elem.node[0]->y);

  length=sqrt(dx*dx+dy*dy);

  return length;
}/*walllength*/

double wallheight(struct element elem)
/*HIGHT OF WALL ELEMENT.*/
{
  double height;

  height=fabs((elem.node[1]->z)-(elem.node[0]->z));

  return height;
}/*wallheight*/

void translatesection(struct section sect,
                      struct materials m,
                      double *As,double *Ac,double *Ar,     /*AREAS*/
                      double *Yg,               /*CENTER OF GRAVITY*/
                      int axis)
/*TRANSLATE SECTION DATA INTO EACH AREA OF MATERIAL etc.*/
{
  int i;
  double sBi[MAXSRECT],sDi[MAXSRECT];       /*WIDTH,DEPTH OF STEELS*/
  double cBi[MAXCRECT],cDi[MAXCRECT];    /*WIDTH,DEPTH OF CONCRETES*/
  double rAi[MAXREINS];
  double sYi[MAXSRECT],sYj[MAXSRECT];         /*UPPER,LOWER OF RECT*/
  double rYi[MAXREINS];
  double cYi[MAXCRECT],cYj[MAXCRECT];
  double sYt,sYc,rYt,rYc,cYt,cYc;
  double dA,Ag;

  Ag=0.0; *Yg=0.0;

  *As=0.0; /*STEEL AREA TOTAL[cm2]*/
  sYc=-1000.0; sYt=1000.0; /*[cm]*/
  for(i=0;i<(sect.nsteel);i++)
  {
    if(axis==SX) /*AROUND x.*/
    {
      sBi[i]=fabs(sect.srect[i].right
                 -sect.srect[i].left);
      sYi[i]=sect.srect[i].top;
      sYj[i]=sect.srect[i].bottom;
    }
    if(axis==SY) /*AROUND y.*/
    {
      sBi[i]=fabs(sect.srect[i].top
                 -sect.srect[i].bottom);
      sYi[i]=sect.srect[i].right;
      sYj[i]=sect.srect[i].left;
    }
    sDi[i]=sYi[i]-sYj[i];

    *As+=(sBi[i]*sDi[i]); /*[cm2]*/

    dA=(sBi[i]*sDi[i])*(m.sE);
    Ag+=dA;
    *Yg+=dA*0.5*(sYi[i]+sYj[i]);

    if(sYc<sYi[i]) sYc=sYi[i];
    if(sYt>sYj[i]) sYt=sYj[i];
  }

  *Ac=0.0; /*CONCRETE AREA TOTAL[cm2]*/
  cYc=-1000.0; cYt=1000.0; /*[cm]*/
  for(i=0;i<(sect.nconc);i++)
  {
    if(axis==SX)
    {
      cBi[i]=fabs(sect.crect[i].right
                 -sect.crect[i].left);
      cYi[i]=sect.crect[i].top;
      cYj[i]=sect.crect[i].bottom;
    }
    if(axis==SY)
    {
      cBi[i]=fabs(sect.crect[i].top
                 -sect.crect[i].bottom);
      cYi[i]=sect.crect[i].right;
      cYj[i]=sect.crect[i].left;
    }
    cDi[i]=cYi[i]-cYj[i];

    *Ac+=(cBi[i]*cDi[i]); /*[cm2]*/

    dA=(cBi[i]*cDi[i])*(m.cE);
    Ag+=dA;
    *Yg+=dA*0.5*(cYi[i]+cYj[i]);

    if(cYc<cYi[i]) cYc=cYi[i];
    if(cYt>cYj[i]) cYt=cYj[i];
  }

  *Ar=0.0;
  rYc=-1000.0; rYt=1000.0; /*[cm]*/
  for(i=0;i<(sect.nrein);i++)
  {
    rAi[i]=sect.rein[i].area;
    if(axis==SX) rYi[i]=sect.rein[i].y;
    if(axis==SY) rYi[i]=sect.rein[i].x;

    *Ar+=rAi[i];

    dA=rAi[i]*(m.rE);
    Ag+=dA;
    *Yg+=dA*rYi[i];

    if(rYc<rYi[i]) rYc=rYi[i];
    if(rYt>rYi[i]) rYt=rYi[i];
  }

  *Yg/=Ag;

  return;
}/*translatesection*/

double allowablebendingofrc(struct element elem,
                            struct materials m,
                            int axis,                     /*0:x 1:y*/
                            double Nd /*[kgf]*/)
/*RETURN:ALLOWABLE BENDING OF RC COLUMN*/
/*N=+:TENSION -:COMPRESSION [kgf]*/
{
  char code[20]="\0";
  int i,ii;
  double sBi[MAXSRECT],sDi[MAXSRECT];       /*WIDTH,DEPTH OF STEELS*/
  double cBi[MAXCRECT],cDi[MAXCRECT];    /*WIDTH,DEPTH OF CONCRETES*/
  double rAi[MAXREINS];
  double Yg,Yn[16],rcYn[6],yn;     /*CENTER OF GRAVITY,NEUTRAL AXIS*/
  double sYi[MAXSRECT],sYj[MAXSRECT];         /*UPPER,LOWER OF RECT*/
  double rYi[MAXREINS];
  double cYi[MAXCRECT],cYj[MAXCRECT],cYb[2*MAXCRECT];
  int nbound;                               /*SUPPLEMENTARY BOUNDS.*/
  double rYt,rYc,cYt,cYc;
  double ret,rec,cec;                         /*MAX STRAINS*/

  double C0,C1,C2,C3;
  double sfi[MAXSRECT],sfj[MAXSRECT];  /*STRESS UPPER,LOWER OF RECT*/
  double rfi[MAXREINS];
  double cfi[MAXCRECT],cfj[MAXCRECT];
  double rft,rfc,cfc;                  /*ALLOWABLE STRESSES*/
  double N[16],rcN[6],Nover;                              /*TENSION*/
  double M[16];
  double Ng,Mg,Ma=0.0;                                /*a:ALLOWABLE*/

  double cB,cBYi,cBYYi;
  double Ac,Ar;
  double rAY,cAupper,cAY;

  int nanswer,found; /*ANSWERS OF CUBIC EQUATION.*/
  double answer[3];

  Yg=0.0;

  Ac=0.0; /*CONCRETE AREA TOTAL[cm2]*/
  cYc=-1000.0; cYt=1000.0; /*[cm]*/
  for(i=0;i<(elem.sect->nconc);i++)
  {
    if(axis==0)
    {
      cBi[i]=fabs(elem.sect->crect[i].right
                 -elem.sect->crect[i].left);
      cYi[i]=elem.sect->crect[i].top;
      cYj[i]=elem.sect->crect[i].bottom;
    }
    if(axis==1)
    {
      cBi[i]=fabs(elem.sect->crect[i].top
                 -elem.sect->crect[i].bottom);
      cYi[i]=elem.sect->crect[i].right;
      cYj[i]=elem.sect->crect[i].left;
    }
    cDi[i]=cYi[i]-cYj[i];

    Ac+=(cBi[i]*cDi[i]); /*[cm2]*/

    if(cYc<cYi[i]) cYc=cYi[i];
    if(cYt>cYj[i]) cYt=cYj[i];
  }

  Ar=0.0;
  rAY=0.0;
  rYc=-1000.0; rYt=1000.0; /*[cm]*/
  for(i=0;i<(elem.sect->nrein);i++)
  {
    rAi[i]=elem.sect->rein[i].area;
    if(axis==0) rYi[i]=elem.sect->rein[i].y;
    if(axis==1) rYi[i]=elem.sect->rein[i].x;

    Ar+=rAi[i];
    rAY+=rAi[i]*rYi[i];

    if(rYc<rYi[i]) rYc=rYi[i];
    if(rYt>rYi[i]) rYt=rYi[i];
  }

  rft=m.rft; ret=rft/m.rE;
  rfc=m.rfc; rec=rfc/m.rE;

  if(elem.sect->stype==RC) cfc=m.cfc;
  else                     return 0.0;
  cec=cfc/m.cE;

  /*
  fprintf(fout0,"Ar=%.3f rYc=%.3f rYt=%.3f[cm2,cm]\n",Ar,rYc,rYt);
  fprintf(fout0,"Ac=%.3f cYc=%.3f cYt=%.3f[cm2,cm]\n",Ac,cYc,cYt);
  */

  /*BOUNDARY (iii).*/
  Yn[3]=cYc;

  for(i=0;i<(elem.sect->nconc);i++)
  {
    cfi[i]=0.0;
    cfj[i]=0.0;
  }
  for(i=0;i<(elem.sect->nrein);i++)
  {
    rfi[i]=rft;
  }
  boundaryvalue(elem,m,Yn[3],sBi,sDi,cBi,cDi,
                             sYi,sYj,cYi,cYj,
                             sfi,sfj,cfi,cfj,
                             rAi,rYi,rfi,
                             &N[3],&M[3]);

  /*BOUNDARY (v).*/
  Yn[7]=cYc;

  for(i=0;i<(elem.sect->nconc);i++)
  {
    cfi[i]=0.0;
    cfj[i]=0.0;
  }
  for(i=0;i<(elem.sect->nrein);i++)
  {
    rfi[i]=(Yn[7]-rYi[i])/(Yn[7]-rYt)*rft;
  }
  boundaryvalue(elem,m,Yn[7],sBi,sDi,cBi,cDi,
                             sYi,sYj,cYi,cYj,
                             sfi,sfj,cfi,cfj,
                             rAi,rYi,rfi,
                             &N[7],&M[7]);

  /*BOUNDARY (vii).*/
  Yn[10]=(ret*cYc-cec*rYt)/(ret-cec);

  for(i=0;i<(elem.sect->nconc);i++)
  {
    cfi[i]=(cYi[i]-Yn[10])/(cYc-Yn[10])*cfc;
    cfj[i]=(cYj[i]-Yn[10])/(cYc-Yn[10])*cfc;
  }
  for(i=0;i<(elem.sect->nrein);i++)
  {
    rfi[i]=(Yn[10]-rYi[i])/(Yn[10]-rYt)*rft;
  }
  boundaryvalue(elem,m,Yn[10],sBi,sDi,cBi,cDi,
                              sYi,sYj,cYi,cYj,
                              sfi,sfj,cfi,cfj,
                              rAi,rYi,rfi,
                              &N[10],&M[10]);

  /*BOUNDARY (xi).*/
  Yn[13]=cYt;

  for(i=0;i<(elem.sect->nconc);i++)
  {
    cfi[i]=(cYi[i]-Yn[13])/(cYc-Yn[13])*cfc;
    cfj[i]=(cYj[i]-Yn[13])/(cYc-Yn[13])*cfc;
  }
  for(i=0;i<(elem.sect->nrein);i++)
  {
    rfi[i]=(Yn[13]-rYi[i])/(Yn[13]-cYc)*cec*(m.rE);
  }
  boundaryvalue(elem,m,Yn[13],sBi,sDi,cBi,cDi,
                              sYi,sYj,cYi,cYj,
                              sfi,sfj,cfi,cfj,
                              rAi,rYi,rfi,
                              &N[13],&M[13]);

  /*BOUNDARY (xiii).*/
  Yn[14]=cYt;

  for(i=0;i<(elem.sect->nconc);i++)
  {
    cfi[i]=cfc;
    cfj[i]=cfc;
  }
  for(i=0;i<(elem.sect->nrein);i++)
  {
    rfi[i]=cec*(m.rE);
  }

  boundaryvalue(elem,m,Yn[14],sBi,sDi,cBi,cDi,
                              sYi,sYj,cYi,cYj,
                              sfi,sfj,cfi,cfj,
                              rAi,rYi,rfi,
                              &N[14],&M[14]);
  /*
  fprintf(fout0,"BOUND 3:Yn=%11.3f N=%12.3f M=%12.3f\n",
          Yn[3],N[3],M[3]);
  fprintf(fout0,"BOUND 7:Yn=%11.3f N=%12.3f M=%12.3f\n",
          Yn[7],N[7],M[7]);
  fprintf(fout0,"BOUND10:Yn=%11.3f N=%12.3f M=%12.3f\n",
          Yn[10],N[10],M[10]);
  fprintf(fout0,"BOUND13:Yn=%11.3f N=%12.3f M=%12.3f\n",
          Yn[13],N[13],M[13]);
  fprintf(fout0,"BOUND14:Yn=%11.3f N=%12.3f M=%12.3f\n",
          Yn[14],N[14],M[14]);
  */

  /*RANGE REDEFINITION FOR RC.*/
  rcN[1]=N[3];  rcYn[1]=Yn[3];
  rcN[2]=N[7];  rcYn[2]=Yn[7];
  rcN[3]=N[10]; rcYn[3]=Yn[10];
  rcN[4]=N[13]; rcYn[4]=Yn[13];
  rcN[5]=N[14]; rcYn[5]=Yn[14];

  /*RANGE DETERMINATION.*/
  if(rcN[1]<=Nd)
  {
    Nover=Nd-rcN[1];
    return 0.0;
  }
  else if(Nd<=rcN[5])
  {
    Nover=Nd-rcN[5];
    return 0.0;
  }
  else Nover=0.0;

  if(rcN[2]<Nd && Nd<=rcN[1]) /*I:(iii)...(v)*/
  {
    C3=0.0;
    C2=0.0;
    C1=Nd-rft*Ar;
    C0=-Nd*rYt+rft*rAY;

    if(C1==0.0) return 0.0;
    yn=-C0/C1;
    if(yn<=rcYn[2])
    {
      fprintf(fout0,"I:NO ANSWER.Yn=%.5f\n",yn);
      return 0.0;
    }

    for(i=0;i<(elem.sect->nconc);i++)
    {
      cfi[i]=0.0;
      cfj[i]=0.0;
    }
    for(i=0;i<(elem.sect->nrein);i++)
    {
      rfi[i]=(yn-rYi[i])/(yn-rYt)*rft;
    }

    sprintf(code,"(I)  ");
    Yg=boundaryvalue(elem,m,yn,sBi,sDi,cBi,cDi,
                               sYi,sYj,cYi,cYj,
                               sfi,sfj,cfi,cfj,
                               rAi,rYi,rfi,
                               &Ng,&Mg);
    Ma=Mg;
  }

  else if(rcN[3]<=Nd && Nd<=rcN[2]) /*II:(v)...(vii)*/
  {
    nbound=0;
    cYb[0]=rcYn[2];
    for(i=0;i<(elem.sect->nconc);i++)   /*GET SUPPLEMENTARY BOUNDS.*/
    {
      if(rcYn[3]<cYi[i] && cYi[i]<rcYn[2])
      {
        nbound++;
        cYb[nbound]=cYi[i];
      }
      if(rcYn[3]<cYj[i] && cYj[i]<rcYn[2])
      {
        nbound++;
        cYb[nbound]=cYj[i];
      }
    }
    nbound++;
    cYb[nbound]=rcYn[3];

    sortdouble(cYb,(nbound+1));

    found=0;
    for(i=0;(!found && i<nbound);i++)
    {
      cB=0.0; cBYi=0.0; cBYYi=0.0;
      cAupper=0.0; cAY=0.0;
      for(ii=0;ii<(elem.sect->nconc);ii++) /*GET SUMS.*/
      {
        if(cYj[ii]>=cYb[i])
        {
          cAupper+=(cBi[i]*cDi[i]); /*[cm2]*/
          cAY+=(cBi[i]*cDi[i])*(cYi[i]+cYj[i]); /*[cm3]*/
        }
        else if(cYi[ii]>cYb[i+1])
        {
          cB+=cBi[i]; /*[cm]*/
          cBYi+=cBi[i]*cYi[i]; /*[cm2]*/
          cBYYi+=cBi[i]*cYi[i]*cYi[i]; /*[cm3]*/
        }
      }

      C3=0.0;
      C2=0.5*ret*(m.cE)*cB;
      C1=Nd-rft*Ar-ret*(m.cE)*cAupper-ret*(m.cE)*cBYi;
      C0=-Nd*rYt+rft*rAY+0.5*ret*(m.cE)*cAY+0.5*ret*(m.cE)*cBYYi;

      nanswer=quadraticequation(C2,C1,C0,answer);
      for(ii=0;ii<nanswer;ii++)
      {
        if(cYb[i]>=answer[ii] && answer[ii]>=cYb[i+1])
        {
          /*fprintf(fout0,"CHOSEN ANSWER=%.5f\n",answer[ii]);*/
          yn=answer[ii];
          found=1;
        }
      }
    }
    if(!found)
    {
      fprintf(fout0,"II:NO ANSWER.\n");
      return 0.0;
    }

    for(i=0;i<(elem.sect->nconc);i++)
    {
      cfi[i]=(cYi[i]-yn)/(rYt-yn)*ret*(m.cE);
      cfj[i]=(cYj[i]-yn)/(rYt-yn)*ret*(m.cE);
    }
    for(i=0;i<(elem.sect->nrein);i++)
    {
      rfi[i]=(yn-rYi[i])/(yn-rYt)*rft;
    }

    sprintf(code,"(II) ");
    Yg=boundaryvalue(elem,m,yn,sBi,sDi,cBi,cDi,
                               sYi,sYj,cYi,cYj,
                               sfi,sfj,cfi,cfj,
                               rAi,rYi,rfi,
                               &Ng,&Mg);
    Ma=Mg;
  }

  else if(rcN[4]<=Nd && Nd<=rcN[3]) /*III:(vii)...(xi)*/
  {
    nbound=0;
    cYb[0]=rcYn[3];
    for(i=0;i<(elem.sect->nconc);i++)   /*GET SUPPLEMENTARY BOUNDS.*/
    {
      if(rcYn[4]<cYi[i] && cYi[i]<rcYn[3])
      {
        nbound++;
        cYb[nbound]=cYi[i];
      }
      if(rcYn[4]<cYj[i] && cYj[i]<rcYn[3])
      {
        nbound++;
        cYb[nbound]=cYj[i];
      }
    }
    nbound++;
    cYb[nbound]=rcYn[4];

    sortdouble(cYb,(nbound+1));

    found=0;
    for(i=0;(!found && i<nbound);i++)
    {
      cB=0.0; cBYi=0.0; cBYYi=0.0;
      cAupper=0.0; cAY=0.0;
      for(ii=0;ii<(elem.sect->nconc);ii++) /*GET SUMS.*/
      {
        if(cYj[ii]>=cYb[i])
        {
          cAupper+=(cBi[i]*cDi[i]);             /*[cm2]*/
          cAY+=(cBi[i]*cDi[i])*(cYi[i]+cYj[i]); /*[cm3]*/
        }
        else if(cYi[ii]>cYb[i+1])
        {
          cB+=cBi[i];                  /*[cm]*/
          cBYi+=cBi[i]*cYi[i];         /*[cm2]*/
          cBYYi+=cBi[i]*cYi[i]*cYi[i]; /*[cm3]*/
        }
      }

      C3=0.0;
      C2=0.5*cfc*cB;
      C1=Nd-cec*(m.rE)*Ar-cfc*cAupper-cfc*cBYi;
      C0=-Nd*cYc+cec*(m.rE)*rAY+0.5*cfc*cAY+0.5*cfc*cBYYi;

      nanswer=quadraticequation(C2,C1,C0,answer);
      for(ii=0;ii<nanswer;ii++)
      {
        if(cYb[i]>=answer[ii] && answer[ii]>=cYb[i+1])
        {
          /*fprintf(fout0,"CHOSEN ANSWER=%.5f\n",answer[ii]);*/
          yn=answer[ii];
          found=1;
        }
      }
    }
    if(!found)
    {
      fprintf(fout0,"III:NO ANSWER.\n");
      return 0.0;
    }

    for(i=0;i<(elem.sect->nconc);i++)
    {
      cfi[i]=(cYi[i]-yn)/(cYc-yn)*cfc;
      cfj[i]=(cYj[i]-yn)/(cYc-yn)*cfc;
    }
    for(i=0;i<(elem.sect->nrein);i++)
    {
      rfi[i]=(yn-rYi[i])/(yn-cYc)*cec*(m.rE);
    }

    sprintf(code,"(III)");
    Yg=boundaryvalue(elem,m,yn,sBi,sDi,cBi,cDi,
                               sYi,sYj,cYi,cYj,
                               sfi,sfj,cfi,cfj,
                               rAi,rYi,rfi,
                               &Ng,&Mg);
    Ma=Mg;
  }

  else if(rcN[5]<=Nd && Nd<=rcN[4]) /*IV:(xi)...(xiii)*/
  {
    nbound=0;
    cYb[0]=rcYn[4];
    nbound++;
    cYb[nbound]=rcYn[5];

    sortdouble(cYb,(nbound+1));

    found=0;
    for(i=0;(!found && i<nbound);i++)
    {
      cB=0.0; cBYi=0.0; cBYYi=0.0;
      cAupper=0.0; cAY=0.0;
      for(ii=0;ii<(elem.sect->nconc);ii++) /*GET SUMS.*/
      {
        if(cYj[ii]>=cYb[i])
        {
          cAupper+=(cBi[i]*cDi[i]);             /*[cm2]*/
          cAY+=(cBi[i]*cDi[i])*(cYi[i]+cYj[i]); /*[cm3]*/
        }
        else if(cYi[ii]>cYb[i+1])
        {
          cB+=cBi[i];                  /*[cm]*/
          cBYi+=cBi[i]*cYi[i];         /*[cm2]*/
          cBYYi+=cBi[i]*cYi[i]*cYi[i]; /*[cm3]*/
        }
      }

      C3=0.0;
      C2=0.0;
      C1=Nd-cec*(m.rE)*Ar-cfc*cAupper;
      C0=-Nd*cYc+cec*(m.rE)*rAY+0.5*cfc*cAY;

      if(C1!=0.0)
      {
        yn=-C0/C1;
        if(yn<=rcYn[4])
        {
          /*fprintf(fout0,"CHOSEN ANSWER=%.5f\n",answer[ii]);*/
          found=1;
        }
      }
    }
    if(!found)
    {
      fprintf(fout0,"IV:NO ANSWER.\n");
      return 0.0;
    }

    for(i=0;i<(elem.sect->nconc);i++)
    {
      cfi[i]=(cYi[i]-yn)/(cYc-yn)*cfc;
      cfj[i]=(cYj[i]-yn)/(cYc-yn)*cfc;
    }
    for(i=0;i<(elem.sect->nrein);i++)
    {
      rfi[i]=(yn-rYi[i])/(yn-cYc)*cec*(m.rE);
    }

    sprintf(code,"(IV) ");
    Yg=boundaryvalue(elem,m,yn,sBi,sDi,cBi,cDi,
                               sYi,sYj,cYi,cYj,
                               sfi,sfj,cfi,cfj,
                               rAi,rYi,rfi,
                               &Ng,&Mg);
    Ma=Mg;
  }

  /*else
  {
    fprintf(fout0,"OUT OF RANGE.\n");
    return 0.0;
  }*/

  /*fprintf(fout0,"RC軸力 rcN=%11.3f[tf]",-Nd/1000.0);
  fprintf(fout0," 重心:%11.3f",Yg);
  fprintf(fout0," 中立軸範囲:%s 中立軸=%11.3f",code,yn);
  fprintf(fout0," 許容値 rcMa=%11.3f[tfm]\n",Ma/100000.0);*/

  /*fprintf(fout0,"B=%.3f D=%.3f ",cBi[0],(cYi[0]-cYj[0]));
  fprintf(fout0,"N/bD=%.3f Ma/bD2=%.3f\n",
          Nd/cBi[0]/(cYi[0]-cYj[0]),
          Ma/cBi[0]/(cYi[0]-cYj[0])/(cYi[0]-cYj[0]));*/

  /*fprintf(fout0," %11.3f %11.3f 0.0\n",Ma/100000.0,-Nd/1000.0);
  fprintf(fout0,"%11.3f %11.3f 0.0",Ma/100000.0,-Nd/1000.0);*/

  return Ma;
}/*allowablebendingofrc*/

double allowablebendingofsrc(struct element elem,
                             struct materials m,
                             int axis,                    /*0:x 1:y*/
                             double Nd /*[kgf]*/)
/*RETURN:ALLOWABLE BENDING OF SRC COLUMN*/
/*N=+:TENSION -:COMPRESSION*/
{
  char code[20]="\0";
  int i,ii;
  double sBi[MAXSRECT],sDi[MAXSRECT];       /*WIDTH,DEPTH OF STEELS*/
  double cBi[MAXCRECT],cDi[MAXCRECT];    /*WIDTH,DEPTH OF CONCRETES*/
  double rAi[MAXREINS];
  double Yg,Yn[16],srcYn[7],yn;    /*CENTER OF GRAVITY,NEUTRAL AXIS*/
  double sYi[MAXSRECT],sYj[MAXSRECT];         /*UPPER,LOWER OF RECT*/
  double rYi[MAXREINS];
  double cYi[MAXCRECT],cYj[MAXCRECT],cYb[2*MAXCRECT];
  int nbound;                               /*SUPPLEMENTARY BOUNDS.*/
  double sYt,sYc,rYt,rYc,cYt,cYc;
  double set,sec,ret,rec,cec;                         /*MAX STRAINS*/

  double C0,C1,C2,C3;
  double sfi[MAXSRECT],sfj[MAXSRECT];  /*STRESS UPPER,LOWER OF RECT*/
  double rfi[MAXREINS];
  double cfi[MAXCRECT],cfj[MAXCRECT];
  double sft,sfc,rft,rfc,cfc;                  /*ALLOWABLE STRESSES*/
  double N[16],srcN[7],Nover;                             /*TENSION*/
  double M[16];
  /*double Ns,Ms,Nr,Mr,Nc,Mc;*//*s:STEEL r:REINFORCEMENT c:CONCRETE*/
  double Ng,Mg,Ma=0.0;                                /*a:ALLOWABLE*/

  double phi5,phi7,phi9,phi10; /*CURVATURE*/

  double cB,cBYi,cBYYi;
  double As,Ac,Ar;
  double sAY,rAY,cAupper,cAY;

  int nanswer,found; /*ANSWERS OF CUBIC EQUATION.*/
  double answer[3];

  Yg=0.0;

  As=0.0; /*STEEL AREA TOTAL[cm2]*/
  sAY=0.0;
  sYc=-1000.0; sYt=1000.0; /*[cm]*/
  for(i=0;i<(elem.sect->nsteel);i++)
  {
    if(axis==0)
    {
      sBi[i]=fabs(elem.sect->srect[i].right
                 -elem.sect->srect[i].left);
      sYi[i]=elem.sect->srect[i].top;
      sYj[i]=elem.sect->srect[i].bottom;
    }
    if(axis==1)
    {
      sBi[i]=fabs(elem.sect->srect[i].top
                 -elem.sect->srect[i].bottom);
      sYi[i]=elem.sect->srect[i].right;
      sYj[i]=elem.sect->srect[i].left;
    }
    sDi[i]=sYi[i]-sYj[i];

    As+=(sBi[i]*sDi[i]); /*[cm2]*/

    sAY+=(sBi[i]*sDi[i])*(sYi[i]+sYj[i]); /*[cm3]*/

    if(sYc<sYi[i]) sYc=sYi[i];
    if(sYt>sYj[i]) sYt=sYj[i];
  }

  Ac=0.0; /*CONCRETE AREA TOTAL[cm2]*/
  cYc=-1000.0; cYt=1000.0; /*[cm]*/
  for(i=0;i<(elem.sect->nconc);i++)
  {
    if(axis==0)
    {
      cBi[i]=fabs(elem.sect->crect[i].right
                 -elem.sect->crect[i].left);
      cYi[i]=elem.sect->crect[i].top;
      cYj[i]=elem.sect->crect[i].bottom;
    }
    if(axis==1)
    {
      cBi[i]=fabs(elem.sect->crect[i].top
                 -elem.sect->crect[i].bottom);
      cYi[i]=elem.sect->crect[i].right;
      cYj[i]=elem.sect->crect[i].left;
    }
    cDi[i]=cYi[i]-cYj[i];

    Ac+=(cBi[i]*cDi[i]); /*[cm2]*/

    if(cYc<cYi[i]) cYc=cYi[i];
    if(cYt>cYj[i]) cYt=cYj[i];
  }

  Ar=0.0;
  rAY=0.0;
  rYc=-1000.0; rYt=1000.0; /*[cm]*/
  for(i=0;i<(elem.sect->nrein);i++)
  {
    rAi[i]=elem.sect->rein[i].area;
    if(axis==0) rYi[i]=elem.sect->rein[i].y;
    if(axis==1) rYi[i]=elem.sect->rein[i].x;

    Ar+=rAi[i];
    rAY+=rAi[i]*rYi[i];

    if(rYc<rYi[i]) rYc=rYi[i];
    if(rYt>rYi[i]) rYt=rYi[i];
  }

  sft=m.sft; set=sft/m.sE;
  sfc=m.sfc; sec=sfc/m.sE;
  rft=m.rft; ret=rft/m.rE;
  rfc=m.rfc; rec=rfc/m.rE;
  if(elem.sect->stype==RC)
  {
    cfc=m.cfc;
  }
  else if(elem.sect->stype==SRC)
  {
    cfc=m.cfc*(1.0-15.0*(As/2.0/Ac)); /*"SRC STANDARD"(29)*/
  }
  cec=cfc/m.cE;
/*
  Ns=0.0; Ms=0.0;
  Nr=0.0; Mr=0.0;
  Nc=0.0; Mc=0.0;
*/
  phi5 =set/(cYc-sYt);
  phi7 =ret/(cYc-rYt);
  phi9 =(set-cec)/(cYc-sYt);
  phi10=(ret-cec)/(cYc-rYt);

  /*
  fprintf(fout0,"As=%.3f sYc=%.3f sYt=%.3f[cm2,cm]\n",As,sYc,sYt);
  fprintf(fout0,"Ar=%.3f rYc=%.3f rYt=%.3f[cm2,cm]\n",Ar,rYc,rYt);
  fprintf(fout0,"Ac=%.3f cYc=%.3f cYt=%.3f[cm2,cm]\n",Ac,cYc,cYt);
  */

  /*BOUNDARY (i):UPPER LIMIT OF TENSION.*/
  Yn[1]=cYc;
  for(i=0;i<(elem.sect->nsteel);i++)
  {
    sfi[i]=sft;
    sfj[i]=sft;
  }
  for(i=0;i<(elem.sect->nconc);i++)
  {
    cfi[i]=0.0;
    cfj[i]=0.0;
  }
  for(i=0;i<(elem.sect->nrein);i++)
  {
    if(set>=ret) rfi[i]=rft;
    else         rfi[i]=set*(m.rE);
  }

  /*fprintf(fout0,"BOUND 1:");*/
  boundaryvalue(elem,m,Yn[1],sBi,sDi,cBi,cDi,
                             sYi,sYj,cYi,cYj,
                             sfi,sfj,cfi,cfj,
                             rAi,rYi,rfi,
                             &N[1],&M[1]);

  /*BOUNDARY (i').*/
  if(set==ret)
  {
    Yn[2]=Yn[1];
    N[2]=N[1];
    M[2]=M[1];
  }
  else
  {
    Yn[2]=(set*rYt-ret*sYt)/(set-ret); /*NEUTRAL AXIS*/

    for(i=0;i<(elem.sect->nsteel);i++)
    {
      sfi[i]=(Yn[2]-sYi[i])/(Yn[2]-sYt)*sft;
      sfj[i]=(Yn[2]-sYj[i])/(Yn[2]-sYt)*sft;
    }
    for(i=0;i<(elem.sect->nconc);i++)
    {
      cfi[i]=0.0;
      cfj[i]=0.0;
    }
    for(i=0;i<(elem.sect->nrein);i++)
    {
      rfi[i]=(Yn[2]-rYi[i])/(Yn[2]-rYt)*rft;
    }

    /*fprintf(fout0,"BOUND 2:");*/
    boundaryvalue(elem,m,Yn[2],sBi,sDi,cBi,cDi,
                               sYi,sYj,cYi,cYj,
                               sfi,sfj,cfi,cfj,
                               rAi,rYi,rfi,
                               &N[2],&M[2]);
  }

  /*BOUNDARY (iii).*/
  Yn[3]=cYc;
  for(i=0;i<(elem.sect->nsteel);i++)
  {
    if(set>=ret)
    {
      sfi[i]=ret*(m.sE);
      sfj[i]=ret*(m.sE);
    }
    else
    {
      sfi[i]=sft;
      sfj[i]=sft;
    }
  }
  for(i=0;i<(elem.sect->nconc);i++)
  {
    cfi[i]=0.0;
    cfj[i]=0.0;
  }
  for(i=0;i<(elem.sect->nrein);i++)
  {
    rfi[i]=rft;
  }

  /*fprintf(fout0,"BOUND 3:");*/
  boundaryvalue(elem,m,Yn[3],sBi,sDi,cBi,cDi,
                             sYi,sYj,cYi,cYj,
                             sfi,sfj,cfi,cfj,
                             rAi,rYi,rfi,
                             &N[3],&M[3]);

  /*BOUNDARY (iii').*/
  Yn[4]=Yn[2];
  N[4]=N[2];
  M[4]=M[2];

  /*BOUNDARY (iv).*/
  Yn[5]=cYc; /*NEUTRAL AXIS*/

  for(i=0;i<(elem.sect->nsteel);i++)
  {
    sfi[i]=(Yn[5]-sYi[i])/(Yn[5]-sYt)*sft;
    sfj[i]=(Yn[5]-sYj[i])/(Yn[5]-sYt)*sft;
  }
  for(i=0;i<(elem.sect->nconc);i++)
  {
    cfi[i]=0.0;
    cfj[i]=0.0;
  }
  for(i=0;i<(elem.sect->nrein);i++)
  {
    if((set/(cYc-sYt))>=(ret/(cYc-rYt)))
    {
      rfi[i]=(Yn[5]-rYi[i])/(Yn[5]-rYt)*rft;
    }
    else
    {
      rfi[i]=(Yn[5]-rYi[i])/(Yn[5]-sYt)*set*(m.rE);
    }
  }

  /*fprintf(fout0,"BOUND 5:");*/
  boundaryvalue(elem,m,Yn[5],sBi,sDi,cBi,cDi,
                             sYi,sYj,cYi,cYj,
                             sfi,sfj,cfi,cfj,
                             rAi,rYi,rfi,
                             &N[5],&M[5]);

  /*BOUNDARY (iv').*/
  if(set==ret)
  {
    Yn[6]=Yn[1];
    N[6]=N[1];
    M[6]=M[1];
  }
  else
  {
    Yn[6]=(set*rYt-ret*sYt)/(set-ret); /*NEUTRAL AXIS*/

    for(i=0;i<(elem.sect->nsteel);i++)
    {
      sfi[i]=(Yn[6]-sYi[i])/(Yn[6]-sYt)*sft;
      sfj[i]=(Yn[6]-sYj[i])/(Yn[6]-sYt)*sft;
    }
    for(i=0;i<(elem.sect->nconc);i++)
    {
      cfi[i]=(Yn[6]-cYi[i])/(Yn[6]-rYt)*ret*(m.cE);
      cfj[i]=(Yn[6]-cYj[i])/(Yn[6]-rYt)*ret*(m.cE);
    }
    for(i=0;i<(elem.sect->nrein);i++)
    {
      rfi[i]=(Yn[6]-rYi[i])/(Yn[6]-rYt)*rft;
    }

    /*fprintf(fout0,"BOUND 6:");*/
    boundaryvalue(elem,m,Yn[6],sBi,sDi,cBi,cDi,
                               sYi,sYj,cYi,cYj,
                               sfi,sfj,cfi,cfj,
                               rAi,rYi,rfi,
                               &N[6],&M[6]);
  }

  /*BOUNDARY (v).*/
  Yn[7]=cYc; /*NEUTRAL AXIS*/

  for(i=0;i<(elem.sect->nsteel);i++)
  {
    if((set/(cYc-sYt))>=(ret/(cYc-rYt)))
    {
      sfi[i]=(Yn[7]-sYi[i])/(Yn[7]-rYt)*ret*(m.sE);
      sfj[i]=(Yn[7]-sYj[i])/(Yn[7]-rYt)*ret*(m.sE);
    }
    else
    {
      sfi[i]=(Yn[7]-sYi[i])/(Yn[7]-sYt)*sft;
      sfj[i]=(Yn[7]-sYj[i])/(Yn[7]-sYt)*sft;
    }
  }
  for(i=0;i<(elem.sect->nconc);i++)
  {
    cfi[i]=0.0;
    cfj[i]=0.0;
  }
  for(i=0;i<(elem.sect->nrein);i++)
  {
    rfi[i]=(Yn[7]-rYi[i])/(Yn[7]-rYt)*rft;
  }

  /*fprintf(fout0,"BOUND 7:");*/
  boundaryvalue(elem,m,Yn[7],sBi,sDi,cBi,cDi,
                             sYi,sYj,cYi,cYj,
                             sfi,sfj,cfi,cfj,
                             rAi,rYi,rfi,
                             &N[7],&M[7]);

  /*BOUNDARY (v').*/
  Yn[8]=Yn[6];
  N[8]=N[6];
  M[8]=M[6];

  /*BOUNDARY (vi),(vi').*/
  if(((set-cec)/(cYc-sYt))>=((ret-cec)/(cYc-rYt)))
  {
    Yn[9]=(ret*cYc-cec*rYt)/(ret-cec); /*NEUTRAL AXIS*/
  }
  else
  {
    Yn[9]=(set*cYc-cec*sYt)/(set-cec); /*NEUTRAL AXIS*/
  }

  for(i=0;i<(elem.sect->nsteel);i++)
  {
    sfi[i]=(Yn[9]-sYi[i])/(Yn[9]-sYt)*sft;
    sfj[i]=(Yn[9]-sYj[i])/(Yn[9]-sYt)*sft;
  }
  for(i=0;i<(elem.sect->nconc);i++)
  {
    cfi[i]=(cYi[i]-Yn[9])/(cYc-Yn[9])*cfc;
    cfj[i]=(cYj[i]-Yn[9])/(cYc-Yn[9])*cfc;
  }
  for(i=0;i<(elem.sect->nrein);i++)
  {
    if(((set-cec)/(cYc-sYt))>=((ret-cec)/(cYc-rYt)))
    {
      rfi[i]=(Yn[9]-rYi[i])/(Yn[9]-rYt)*rft;
    }
    else
    {
      rfi[i]=(Yn[9]-rYi[i])/(Yn[9]-sYt)*set*(m.rE);
    }
  }

  /*fprintf(fout0,"BOUND 9:");*/
  boundaryvalue(elem,m,Yn[9],sBi,sDi,cBi,cDi,
                             sYi,sYj,cYi,cYj,
                             sfi,sfj,cfi,cfj,
                             rAi,rYi,rfi,
                             &N[9],&M[9]);

  /*BOUNDARY (vii).*/
  Yn[10]=(ret*cYc-cec*rYt)/(ret-cec); /*NEUTRAL AXIS*/

  for(i=0;i<(elem.sect->nsteel);i++)
  {
    if(((set-cec)/(cYc-sYt))>=((ret-cec)/(cYc-rYt)))
    {
      sfi[i]=(Yn[10]-sYi[i])/(Yn[10]-rYt)*ret*(m.sE);
      sfj[i]=(Yn[10]-sYj[i])/(Yn[10]-rYt)*ret*(m.sE);
    }
    else
    {
      sfi[i]=(Yn[10]-sYi[i])/(Yn[10]-sYt)*sft;
      sfj[i]=(Yn[10]-sYj[i])/(Yn[10]-sYt)*sft;
    }
  }
  for(i=0;i<(elem.sect->nconc);i++)
  {
    cfi[i]=(cYi[i]-Yn[10])/(cYc-Yn[10])*cfc;
    cfj[i]=(cYj[i]-Yn[10])/(cYc-Yn[10])*cfc;
  }
  for(i=0;i<(elem.sect->nrein);i++)
  {
    rfi[i]=(Yn[10]-rYi[i])/(Yn[10]-rYt)*rft;
  }

  /*fprintf(fout0,"BOUND 10:");*/
  boundaryvalue(elem,m,Yn[10],sBi,sDi,cBi,cDi,
                              sYi,sYj,cYi,cYj,
                              sfi,sfj,cfi,cfj,
                              rAi,rYi,rfi,
                              &N[10],&M[10]);

  /*BOUNDARY (viii).*/
  Yn[11]=(set*sYc-sec*sYt)/(set-sec); /*NEUTRAL AXIS*/

  for(i=0;i<(elem.sect->nsteel);i++)
  {
    sfi[i]=(Yn[11]-sYi[i])/(Yn[11]-sYt)*sft;
    sfj[i]=(Yn[11]-sYj[i])/(Yn[11]-sYt)*sft;
  }
  for(i=0;i<(elem.sect->nconc);i++)
  {
    cfi[i]=(cYi[i]-Yn[11])/(cYc-Yn[11])*cfc;
    cfj[i]=(cYj[i]-Yn[11])/(cYc-Yn[11])*cfc;
  }
  for(i=0;i<(elem.sect->nrein);i++)
  {
    rfi[i]=(Yn[11]-rYi[i])/(Yn[11]-cYc)*cec*(m.rE);
  }

  /*fprintf(fout0,"BOUND 11:");*/
  boundaryvalue(elem,m,Yn[11],sBi,sDi,cBi,cDi,
                              sYi,sYj,cYi,cYj,
                              sfi,sfj,cfi,cfj,
                              rAi,rYi,rfi,
                              &N[11],&M[11]);

  /*BOUNDARY (x).*/
  Yn[12]=cYt; /*NEUTRAL AXIS*/

  for(i=0;i<(elem.sect->nsteel);i++)
  {
    sfi[i]=(Yn[12]-sYi[i])/(Yn[12]-sYc)*sfc;
    sfj[i]=(Yn[12]-sYj[i])/(Yn[12]-sYc)*sfc;
  }
  for(i=0;i<(elem.sect->nconc);i++)
  {
    cfi[i]=(cYi[i]-Yn[12])/(cYc-Yn[12])*cfc;
    cfj[i]=(cYj[i]-Yn[12])/(cYc-Yn[12])*cfc;
  }
  for(i=0;i<(elem.sect->nrein);i++)
  {
    rfi[i]=(Yn[12]-rYi[i])/(Yn[12]-cYc)*cec*(m.rE);
  }

  /*fprintf(fout0,"BOUND 12:");*/
  boundaryvalue(elem,m,Yn[12],sBi,sDi,cBi,cDi,
                              sYi,sYj,cYi,cYj,
                              sfi,sfj,cfi,cfj,
                              rAi,rYi,rfi,
                              &N[12],&M[12]);

  /*BOUNDARY (xi).*/
  Yn[13]=cYt; /*NEUTRAL AXIS*/

  for(i=0;i<(elem.sect->nsteel);i++)
  {
    sfi[i]=(Yn[13]-sYi[i])/(Yn[13]-cYc)*cec*(m.sE);
    sfj[i]=(Yn[13]-sYj[i])/(Yn[13]-cYc)*cec*(m.sE);
  }
  for(i=0;i<(elem.sect->nconc);i++)
  {
    cfi[i]=(cYi[i]-Yn[13])/(cYc-Yn[13])*cfc;
    cfj[i]=(cYj[i]-Yn[13])/(cYc-Yn[13])*cfc;
  }
  for(i=0;i<(elem.sect->nrein);i++)
  {
    rfi[i]=(Yn[13]-rYi[i])/(Yn[13]-cYc)*cec*(m.rE);
  }

  /*fprintf(fout0,"BOUND 13:");*/
  boundaryvalue(elem,m,Yn[13],sBi,sDi,cBi,cDi,
                              sYi,sYj,cYi,cYj,
                              sfi,sfj,cfi,cfj,
                              rAi,rYi,rfi,
                              &N[13],&M[13]);

  /*BOUNDARY (xiii).*/
  Yn[14]=cYt; /*NEUTRAL AXIS*/

  for(i=0;i<(elem.sect->nsteel);i++)
  {
    sfi[i]=cec*(m.sE);
    sfj[i]=cec*(m.sE);
  }
  for(i=0;i<(elem.sect->nconc);i++)
  {
    cfi[i]=cfc;
    cfj[i]=cfc;
  }
  for(i=0;i<(elem.sect->nrein);i++)
  {
    rfi[i]=cec*(m.rE);
  }

  /*fprintf(fout0,"BOUND 14:");*/
  boundaryvalue(elem,m,Yn[14],sBi,sDi,cBi,cDi,
                              sYi,sYj,cYi,cYj,
                              sfi,sfj,cfi,cfj,
                              rAi,rYi,rfi,
                              &N[14],&M[14]);

  /*BOUNDARY (xiv).*/
  Yn[15]=cYt; /*NEUTRAL AXIS*/

  for(i=0;i<(elem.sect->nsteel);i++)
  {
    sfi[i]=sfc;
    sfj[i]=sfc;
  }
  for(i=0;i<(elem.sect->nconc);i++)
  {
    cfi[i]=cfc;
    cfj[i]=cfc;
  }
  for(i=0;i<(elem.sect->nrein);i++)
  {
    rfi[i]=cec*(m.rE);
  }

  /*fprintf(fout0,"BOUND 15:");*/
  boundaryvalue(elem,m,Yn[15],sBi,sDi,cBi,cDi,
                              sYi,sYj,cYi,cYj,
                              sfi,sfj,cfi,cfj,
                              rAi,rYi,rfi,
                              &N[15],&M[15]);

  /*for(i=1;i<=15;i++)
  {
    fprintf(fout0,"BOUND%2d:Yn=%11.3f N=%12.3f M=%12.3f\n",
            i,Yn[i],N[i],M[i]);
  }*/

  /*RANGE REDEFINITION FOR SRC.*/
  if(set>=ret){srcN[1]=N[1]; srcYn[1]=Yn[1];}
  else        {srcN[1]=N[3]; srcYn[1]=Yn[3];}

  if(phi5<=phi7){srcN[2]=N[7]; srcYn[2]=Yn[7];}
  else          {srcN[2]=N[5]; srcYn[2]=Yn[5];}

  if(phi9<=phi10){srcN[3]=N[10]; srcYn[3]=Yn[10];}
  else           {srcN[3]=N[9];  srcYn[3]=Yn[9];}

  srcN[4]=N[11]; srcYn[4]=Yn[11];
  srcN[5]=N[12]; srcYn[5]=Yn[12];
  srcN[6]=N[15]; srcYn[6]=Yn[15];

  /*RANGE DETERMINATION.*/
  if(srcN[1]<=Nd)
  {
    Nover=Nd-srcN[1];
    return 0.0;
  }
  else if(Nd<=srcN[6])
  {
    Nover=Nd-srcN[6];
    return 0.0;
  }
  else Nover=0.0;

  if(srcN[2]<=Nd && Nd<=srcN[1]) /*I:(iii)...(v) etc.*/
  {
    found=0;

    C3=0.0;
    C2=Nd-sft*As-rft*Ar;
    C1=-Nd*(sYt+rYt)+0.5*sft*(2.0*rYt*As+sAY)
      +rft*(sYt*Ar+rAY);
    C0=Nd*sYt*rYt-0.5*sft*(rYt*sAY)-rft*sYt*rAY;

    nanswer=quadraticequation(C2,C1,C0,answer);
    for(ii=0;ii<nanswer;ii++)
    {
      if(answer[ii]>=srcYn[2])
      {
        /*fprintf(fout0,"CHOSEN ANSWER=%.5f\n",answer[ii]);*/
        yn=answer[ii];
        found=1;
      }
    }
    if(!found)
    {
      fprintf(fout0,"I:NO ANSWER.\n");
      return 0.0;
    }

    for(i=0;i<(elem.sect->nsteel);i++)
    {
      sfi[i]=(yn-sYi[i])/(yn-sYt)*sft;
      sfj[i]=(yn-sYj[i])/(yn-sYt)*sft;
    }
    for(i=0;i<(elem.sect->nconc);i++)
    {
      cfi[i]=0.0;
      cfj[i]=0.0;
    }
    for(i=0;i<(elem.sect->nrein);i++)
    {
      rfi[i]=(yn-rYi[i])/(yn-rYt)*rft;
    }

    sprintf(code,"(I)  ");
    Yg=boundaryvalue(elem,m,yn,sBi,sDi,cBi,cDi,
                               sYi,sYj,cYi,cYj,
                               sfi,sfj,cfi,cfj,
                               rAi,rYi,rfi,
                               &Ng,&Mg);
    Ma=Mg;
  }

  else if(srcN[3]<=Nd && Nd<=srcN[2]) /*II:(v)...(vii) etc.*/
  {
    nbound=0;
    cYb[0]=srcYn[2];
    for(i=0;i<(elem.sect->nconc);i++)   /*GET SUPPLEMENTARY BOUNDS.*/
    {
      if(srcYn[3]<cYi[i] && cYi[i]<srcYn[2])
      {
        nbound++;
        cYb[nbound]=cYi[i];
      }
      if(srcYn[3]<cYj[i] && cYj[i]<srcYn[2])
      {
        nbound++;
        cYb[nbound]=cYj[i];
      }
    }
    nbound++;
    cYb[nbound]=srcYn[3];

    sortdouble(cYb,(nbound+1));

    found=0;
    for(i=0;(!found && i<nbound);i++)
    {
      cB=0.0; cBYi=0.0; cBYYi=0.0;
      cAupper=0.0; cAY=0.0;
      for(ii=0;ii<(elem.sect->nconc);ii++) /*GET SUMS.*/
      {
        if(cYj[ii]>=cYb[i])
        {
          cAupper+=(cBi[i]*cDi[i]); /*[cm2]*/

          cAY+=(cBi[i]*cDi[i])*(cYi[i]+cYj[i]); /*[cm3]*/
        }
        else if(cYi[ii]>cYb[i+1])
        {
          cB+=cBi[i]; /*[cm]*/

          cBYi+=cBi[i]*cYi[i]; /*[cm2]*/
          cBYYi+=cBi[i]*cYi[i]*cYi[i]; /*[cm3]*/
        }
      }

      C3=0.5*ret*(m.cE)*cB;
      C2=Nd-sft*As-rft*Ar-ret*(m.cE)*cAupper
        -0.5*ret*(m.cE)*(sYt*cB+2.0*cBYi);
      C1=-Nd*(sYt+rYt)+0.5*sft*(2.0*rYt*As+sAY)
        +rft*(sYt*Ar+rAY)+0.5*ret*(m.cE)*(2.0*sYt*cAupper+cAY)
        +0.5*ret*(m.cE)*(2.0*sYt*cBYi+cBYYi);
      C0=Nd*sYt*rYt-0.5*sft*(rYt*sAY)-rft*sYt*rAY
        -0.5*ret*(m.cE)*(sYt*cAY)-0.5*ret*(m.cE)*sYt*cBYYi;

      nanswer=cubicequation(C3,C2,C1,C0,answer);
      for(ii=0;ii<nanswer;ii++)
      {
        if(cYb[i]>=answer[ii] && answer[ii]>=cYb[i+1])
        {
          /*fprintf(fout0,"CHOSEN ANSWER=%.5f\n",answer[ii]);*/
          yn=answer[ii];
          found=1;
        }
      }
    }
    if(!found)
    {
      fprintf(fout0,"II:NO ANSWER.\n");
      return 0.0;
    }

    for(i=0;i<(elem.sect->nsteel);i++)
    {
      sfi[i]=(yn-sYi[i])/(yn-sYt)*sft;
      sfj[i]=(yn-sYj[i])/(yn-sYt)*sft;
    }
    for(i=0;i<(elem.sect->nconc);i++)
    {
      cfi[i]=(cYi[i]-yn)/(rYt-yn)*ret*(m.cE);
      cfj[i]=(cYj[i]-yn)/(rYt-yn)*ret*(m.cE);
    }
    for(i=0;i<(elem.sect->nrein);i++)
    {
      rfi[i]=(yn-rYi[i])/(yn-rYt)*rft;
    }

    sprintf(code,"(II) ");
    Yg=boundaryvalue(elem,m,yn,sBi,sDi,cBi,cDi,
                               sYi,sYj,cYi,cYj,
                               sfi,sfj,cfi,cfj,
                               rAi,rYi,rfi,
                               &Ng,&Mg);
    Ma=Mg;
  }

  else if(srcN[4]<=Nd && Nd<=srcN[3]) /*III:(vii)...(viii) etc.*/
  {
    nbound=0;
    cYb[0]=srcYn[3];
    for(i=0;i<(elem.sect->nconc);i++)   /*GET SUPPLEMENTARY BOUNDS.*/
    {
      if(srcYn[4]<cYi[i] && cYi[i]<srcYn[3])
      {
        nbound++;
        cYb[nbound]=cYi[i];
      }
      if(srcYn[4]<cYj[i] && cYj[i]<srcYn[3])
      {
        nbound++;
        cYb[nbound]=cYj[i];
      }
    }
    nbound++;
    cYb[nbound]=srcYn[4];

    sortdouble(cYb,(nbound+1));

    found=0;
    for(i=0;(!found && i<nbound);i++)
    {
      cB=0.0; cBYi=0.0; cBYYi=0.0;
      cAupper=0.0; cAY=0.0;
      for(ii=0;ii<(elem.sect->nconc);ii++) /*GET SUMS.*/
      {
        if(cYj[ii]>=cYb[i])
        {
          cAupper+=(cBi[i]*cDi[i]); /*[cm2]*/

          cAY+=(cBi[i]*cDi[i])*(cYi[i]+cYj[i]); /*[cm3]*/
        }
        else if(cYi[ii]>cYb[i+1])
        {
          cB+=cBi[i]; /*[cm]*/

          cBYi+=cBi[i]*cYi[i]; /*[cm2]*/
          cBYYi+=cBi[i]*cYi[i]*cYi[i]; /*[cm3]*/
        }
      }

      C3=0.5*cfc*cB;
      C2=Nd-sft*As-cec*(m.rE)*Ar-cfc*cAupper
        -0.5*cfc*(sYt*cB+2.0*cBYi);
      C1=-Nd*(sYt+cYc)+0.5*sft*(2.0*cYc*As+sAY)
        +cec*(m.rE)*(sYt*Ar+rAY)+0.5*cfc*(2.0*sYt*cAupper+cAY)
        +0.5*cfc*(2.0*sYt*cBYi+cBYYi);
      C0=Nd*sYt*cYc-0.5*sft*(cYc*sAY)-cec*(m.rE)*sYt*rAY
        -0.5*cfc*(sYt*cAY)-0.5*cfc*sYt*cBYYi;

      nanswer=cubicequation(C3,C2,C1,C0,answer);
      for(ii=0;ii<nanswer;ii++)
      {
        if(cYb[i]>=answer[ii] && answer[ii]>=cYb[i+1])
        {
          /*fprintf(fout0,"CHOSEN ANSWER=%.5f\n",answer[ii]);*/
          yn=answer[ii];
          found=1;
        }
      }
    }
    if(!found)
    {
      fprintf(fout0,"III:NO ANSWER.\n");
      return 0.0;
    }

    for(i=0;i<(elem.sect->nsteel);i++)
    {
      sfi[i]=(yn-sYi[i])/(yn-sYt)*sft;
      sfj[i]=(yn-sYj[i])/(yn-sYt)*sft;
    }
    for(i=0;i<(elem.sect->nconc);i++)
    {
      cfi[i]=(cYi[i]-yn)/(cYc-yn)*cfc;
      cfj[i]=(cYj[i]-yn)/(cYc-yn)*cfc;
    }
    for(i=0;i<(elem.sect->nrein);i++)
    {
      rfi[i]=(yn-rYi[i])/(yn-cYc)*cec*(m.rE);
    }

    sprintf(code,"(III)");
    Yg=boundaryvalue(elem,m,yn,sBi,sDi,cBi,cDi,
                               sYi,sYj,cYi,cYj,
                               sfi,sfj,cfi,cfj,
                               rAi,rYi,rfi,
                               &Ng,&Mg);
    Ma=Mg;
  }

  else if(srcN[5]<=Nd && Nd<=srcN[4]) /*IV:(viii)...(x) etc.*/
  {
    nbound=0;
    cYb[0]=srcYn[4];
    for(i=0;i<(elem.sect->nconc);i++)   /*GET SUPPLEMENTARY BOUNDS.*/
    {
      if(srcYn[5]<cYi[i] && cYi[i]<srcYn[4])
      {
        nbound++;
        cYb[nbound]=cYi[i];
      }
      if(srcYn[5]<cYj[i] && cYj[i]<srcYn[4])
      {
        nbound++;
        cYb[nbound]=cYj[i];
      }
    }
    nbound++;
    cYb[nbound]=srcYn[5];

    sortdouble(cYb,(nbound+1));

    found=0;
    for(i=0;(!found && i<nbound);i++)
    {
      cB=0.0; cBYi=0.0; cBYYi=0.0;
      cAupper=0.0; cAY=0.0;
      for(ii=0;ii<(elem.sect->nconc);ii++) /*GET SUMS.*/
      {
        if(cYj[ii]>=cYb[i])
        {
          cAupper+=(cBi[i]*cDi[i]); /*[cm2]*/

          cAY+=(cBi[i]*cDi[i])*(cYi[i]+cYj[i]); /*[cm3]*/
        }
        else if(cYi[ii]>cYb[i+1])
        {
          cB+=cBi[i]; /*[cm]*/

          cBYi+=cBi[i]*cYi[i]; /*[cm2]*/
          cBYYi+=cBi[i]*cYi[i]*cYi[i]; /*[cm3]*/
        }
      }

      C3=0.5*cfc*cB;
      C2=Nd-sfc*As-cec*(m.rE)*Ar-cfc*cAupper
        -0.5*cfc*(sYc*cB+2.0*cBYi);
      C1=-Nd*(sYc+cYc)+0.5*sfc*(2.0*cYc*As+sAY)
        +cec*(m.rE)*(sYc*Ar+rAY)+0.5*cfc*(2.0*sYc*cAupper+cAY)
        +0.5*cfc*(2.0*sYc*cBYi+cBYYi);
      C0=Nd*sYc*cYc-0.5*sfc*(cYc*sAY)-cec*(m.rE)*sYc*rAY
        -0.5*cfc*(sYc*cAY)-0.5*cfc*sYc*cBYYi;

      nanswer=cubicequation(C3,C2,C1,C0,answer);
      for(ii=0;ii<nanswer;ii++)
      {
        if(cYb[i]>=answer[ii] && answer[ii]>=cYb[i+1])
        {
          /*fprintf(fout0,"CHOSEN ANSWER=%.5f\n",answer[ii]);*/
          yn=answer[ii];
          found=1;
        }
      }
    }
    if(!found)
    {
      fprintf(fout0,"IV:NO ANSWER.\n");
      return 0.0;
    }

    for(i=0;i<(elem.sect->nsteel);i++)
    {
      sfi[i]=(yn-sYi[i])/(yn-sYc)*sfc;
      sfj[i]=(yn-sYj[i])/(yn-sYc)*sfc;
    }
    for(i=0;i<(elem.sect->nconc);i++)
    {
      cfi[i]=(cYi[i]-yn)/(cYc-yn)*cfc;
      cfj[i]=(cYj[i]-yn)/(cYc-yn)*cfc;
    }
    for(i=0;i<(elem.sect->nrein);i++)
    {
      rfi[i]=(yn-rYi[i])/(yn-cYc)*cec*(m.rE);
    }

    sprintf(code,"(IV) ");
    Yg=boundaryvalue(elem,m,yn,sBi,sDi,cBi,cDi,
                               sYi,sYj,cYi,cYj,
                               sfi,sfj,cfi,cfj,
                               rAi,rYi,rfi,
                               &Ng,&Mg);
    Ma=Mg;
  }

  else if(srcN[6]<=Nd && Nd<=srcN[5]) /*V:(x)...(xiv) etc.*/
  {
    nbound=0;
    cYb[0]=srcYn[5];
    nbound++;
    cYb[nbound]=srcYn[6];

    found=0;
    for(i=0;(!found && i<nbound);i++)
    {
      cB=0.0; cBYi=0.0; cBYYi=0.0;
      cAupper=0.0; cAY=0.0;
      for(ii=0;ii<(elem.sect->nconc);ii++) /*GET SUMS.*/
      {
        if(cYj[ii]>=cYb[i])
        {
          cAupper+=(cBi[i]*cDi[i]); /*[cm2]*/

          cAY+=(cBi[i]*cDi[i])*(cYi[i]+cYj[i]); /*[cm3]*/
        }
      }

      C3=0.0;
      C2=Nd-sfc*As-cec*(m.rE)*Ar-cfc*cAupper;
      C1=-Nd*(sYc+cYc)+0.5*sfc*(2.0*cYc*As+sAY)
        +cec*(m.rE)*(sYc*Ar+rAY)+0.5*cfc*(2.0*sYc*cAupper+cAY);
      C0=Nd*sYc*cYc-0.5*sfc*(cYc*sAY)-cec*(m.rE)*sYc*rAY
        -0.5*cfc*(sYc*cAY);

      nanswer=quadraticequation(C2,C1,C0,answer);
      for(ii=0;ii<nanswer;ii++)
      {
        if(srcYn[5]>=answer[ii])
        {
          /*fprintf(fout0,"CHOSEN ANSWER=%.5f\n",answer[ii]);*/
          yn=answer[ii];
          found=1;
        }
      }
    }
    if(!found)
    {
      fprintf(fout0,"V:NO ANSWER.\n");
      return 0.0;
    }

    for(i=0;i<(elem.sect->nsteel);i++)
    {
      sfi[i]=(yn-sYi[i])/(yn-sYc)*sfc;
      sfj[i]=(yn-sYj[i])/(yn-sYc)*sfc;
    }
    for(i=0;i<(elem.sect->nconc);i++)
    {
      cfi[i]=(cYi[i]-yn)/(cYc-yn)*cfc;
      cfj[i]=(cYj[i]-yn)/(cYc-yn)*cfc;
    }
    for(i=0;i<(elem.sect->nrein);i++)
    {
      rfi[i]=(yn-rYi[i])/(yn-cYc)*cec*(m.rE);
    }

    sprintf(code,"(V)  ");
    Yg=boundaryvalue(elem,m,yn,sBi,sDi,cBi,cDi,
                               sYi,sYj,cYi,cYj,
                               sfi,sfj,cfi,cfj,
                               rAi,rYi,rfi,
                               &Ng,&Mg);
    Ma=Mg;
  }

  /*fprintf(fout0,"SRC軸力 srcN=%11.3f[tf]",-Nd/1000.0);
  fprintf(fout0," 重心:%11.3f",Yg);
  fprintf(fout0," 中立軸範囲:%s 中立軸=%11.3f",code,yn);
  fprintf(fout0," 許容値 srcMa=%11.3f[tfm]\n",Ma/100000.0);*/

  /*fprintf(fout0," %11.3f %11.3f 0.0\n",Ma/100000.0,-Nd/1000.0);
  fprintf(fout0,"%11.3f %11.3f 0.0",Ma/100000.0,-Nd/1000.0);*/

  return Ma;
}/*allowablebendingofsrc*/

double ultimatebendingofsrc(struct element elem,
                            struct materials m,
                            int axis,                     /*0:x 1:y*/
                            double Nd,                      /*[kgf]*/
                            double *sN,double *sM,
                            double *rN,double *rM,
                            double *cN,double *cM)
/*RETURN:ULTIMATE BENDING OF RC,SRC COLUMN.*/
/*       PLASTIC BENDING OF S COLUMN.*/
/*N=+:TENSION -:COMPRESSION*/
{
  char code[20]="\0";
  int i,j;
  double sBi[MAXSRECT],sDi[MAXSRECT];       /*WIDTH,DEPTH OF STEELS*/
  double cBi[MAXCRECT],cDi[MAXCRECT];    /*WIDTH,DEPTH OF CONCRETES*/
  double rAi[MAXREINS];
  double Yg,Yn,yn;                 /*CENTER OF GRAVITY,NEUTRAL AXIS*/
  double sYi[MAXSRECT],sYj[MAXSRECT];         /*UPPER,LOWER OF RECT*/
  double rYi[MAXREINS];
  double cYi[MAXCRECT],cYj[MAXCRECT];
  int nbound,ibound;                               /*SUPPLEMENTARY BOUNDS.*/
  double bounds[2*MAXSRECT+MAXREINS+2*MAXCRECT];
  double sYt,sYc,rYt,rYc,cYt,cYc;

  /*double sfi[MAXSRECT],sfj[MAXSRECT];*/  /*STRESS UPPER,LOWER OF RECT*/
  /*double rfi[MAXREINS];*/
  /*double cfi[MAXCRECT],cfj[MAXCRECT];*/
  double sftu,sfcu,rftu,rfcu,cfcu;              /*ULTIMATE STRESSES*/
  double N[2*MAXSRECT+MAXREINS+2*MAXCRECT],Nover;         /*TENSION*/
  double M[2*MAXSRECT+MAXREINS+2*MAXCRECT];
  double Ng,Mg,Mu=0.0;                                 /*u:ULTIMATE*/

  double As,Ac,Ar;
  double sAY,rAY;
  double bunbo,bunsi; /*BUNBO,BUNSI OF Yn.*/

  double rAt,rYtg;

  nbound=0;
  Yg=0.0;

  As=0.0; /*STEEL AREA TOTAL[cm2]*/
  sAY=0.0;
  sYc=-1000.0; sYt=1000.0; /*[cm]*/
  for(i=0;i<(elem.sect->nsteel);i++)
  {
    if(axis==0)
    {
      sBi[i]=fabs(elem.sect->srect[i].right
                 -elem.sect->srect[i].left);
      sYi[i]=elem.sect->srect[i].top;
      sYj[i]=elem.sect->srect[i].bottom;
    }
    if(axis==1)
    {
      sBi[i]=fabs(elem.sect->srect[i].top
                 -elem.sect->srect[i].bottom);
      sYi[i]=elem.sect->srect[i].right;
      sYj[i]=elem.sect->srect[i].left;
    }
    sDi[i]=sYi[i]-sYj[i];

    As+=(sBi[i]*sDi[i]); /*[cm2]*/

    sAY+=(sBi[i]*sDi[i])*(sYi[i]+sYj[i]); /*[cm3]*/

    bounds[nbound]=sYi[i]; nbound++;
    bounds[nbound]=sYj[i]; nbound++;

    if(sYc<sYi[i]) sYc=sYi[i];
    if(sYt>sYj[i]) sYt=sYj[i];
  }

  Ac=0.0; /*CONCRETE AREA TOTAL[cm2]*/
  cYc=-1000.0; cYt=1000.0; /*[cm]*/
  for(i=0;i<(elem.sect->nconc);i++)
  {
    if(axis==0)
    {
      cBi[i]=fabs(elem.sect->crect[i].right
                 -elem.sect->crect[i].left);
      cYi[i]=elem.sect->crect[i].top;
      cYj[i]=elem.sect->crect[i].bottom;
    }
    if(axis==1)
    {
      cBi[i]=fabs(elem.sect->crect[i].top
                 -elem.sect->crect[i].bottom);
      cYi[i]=elem.sect->crect[i].right;
      cYj[i]=elem.sect->crect[i].left;
    }
    cDi[i]=cYi[i]-cYj[i];

    Ac+=(cBi[i]*cDi[i]); /*[cm2]*/

    bounds[nbound]=cYi[i]; nbound++;
    bounds[nbound]=cYj[i]; nbound++;

    if(cYc<cYi[i]) cYc=cYi[i];
    if(cYt>cYj[i]) cYt=cYj[i];
  }

  Ar=0.0;
  rAY=0.0;
  rAt=0.0; rYtg=0.0; /*TENSED REIN AREA.*/
  rYc=-1000.0; rYt=1000.0; /*[cm]*/
  for(i=0;i<(elem.sect->nrein);i++)
  {
    rAi[i]=elem.sect->rein[i].area;

    if(axis==0) rYi[i]=elem.sect->rein[i].y;
    if(axis==1) rYi[i]=elem.sect->rein[i].x;

/*if(elem.sect->code==101  ||
   elem.sect->code==1021 ||
   elem.sect->code==1022 ||
   elem.sect->code==103  ||
   elem.sect->code==1041 ||
   elem.sect->code==1042 ||
   elem.sect->code==105  ||
   elem.sect->code==106)*/ /*FOR HAKODATE.*/
{
  if(rYi[i]<=0.0) rAi[i]=elem.sect->rein[i].area;
  else            rAi[i]=0.0;
}

    Ar +=rAi[i];
    rAY+=rAi[i]*rYi[i];

    if(rYi[i]<0.0)
    {
      rAt +=rAi[i];
      rYtg+=rAi[i]*rYi[i];
    }

    bounds[nbound]=rYi[i]; nbound++;

    if(rYc<rYi[i]) rYc=rYi[i];
    if(rYt>rYi[i]) rYt=rYi[i];
  }
  rYtg/=rAt;

  sortdouble(bounds,nbound); /*SORT BOUNDARIES.*/

  sftu=m.sftu;
  sfcu=m.sfcu;
  rftu=m.rftu;
  rfcu=m.rfcu;
  if(elem.sect->stype==RC)
  {
    cfcu=-0.85*m.Fc;
  }
  else if(elem.sect->stype==SRC)
  {
    cfcu=-m.Fc*(0.85-2.5*(As/2.0/Ac)); /*"SRC STANDARD"*/
  }

  /*
  fprintf(fout0,"As=%.3f sYc=%.3f sYt=%.3f[cm2,cm]\n",As,sYc,sYt);
  fprintf(fout0,"Ar=%.3f rYc=%.3f rYt=%.3f[cm2,cm]\n",Ar,rYc,rYt);
  fprintf(fout0,"Ac=%.3f cYc=%.3f cYt=%.3f[cm2,cm]\n",Ac,cYc,cYt);
  */

  for(i=0;i<nbound;i++) /*BOUNDARY DEFINITION.*/
  {
    Yn=bounds[i];

    boundaryvalueultimate(elem,m,Yn,sBi,sDi,cBi,cDi,
                                    sYi,sYj,cYi,cYj,
                                    sftu,sfcu,
                                    rftu,rfcu,
                                    0.0,cfcu,
                                    rAi,rYi,
                                    sN,sM,rN,rM,cN,cM,
                                    &N[i],&M[i]);
  }

  /*for(i=0;i<nbound;i++)
  {
    fprintf(fout0,"BOUND%2d:Yn=%11.3f N=%12.3f M=%12.3f\n",
            i,bounds[i],N[i],M[i]);
  }*/

  /*RANGE DETERMINATION.*/
  if(N[0]<=Nd)
  {
    Nover=Nd-N[0];
    return 0.0;
  }
  else if(Nd<=N[nbound-1])
  {
    Nover=Nd-N[nbound-1];
    return 0.0;
  }
  else Nover=0.0;

  ibound=0;
  while(N[ibound+1]>=Nd) ibound++;

  if(Nd==N[ibound])
  {
    yn=bounds[ibound];
  }
  else
  {
    bunsi=Nd;
    bunbo=0.0;

    for(j=0;j<(elem.sect->nsteel);j++)
    {
      if(sYi[j]<=bounds[ibound+1]) /*TENSED.*/
      {
        bunsi-=sftu*sBi[j]*sDi[j];
      }
      else if(sYj[j]>=bounds[ibound]) /*COMPRESSED.*/
      {
        bunsi-=sfcu*sBi[j]*sDi[j];
      }
      else /*ACROSSING.*/
      {
        bunsi+=sYj[j]*sBi[j]*sftu-sYi[j]*sBi[j]*sfcu;
        bunbo+=sBi[j]*sftu-sBi[j]*sfcu;
      }
    }
    for(j=0;j<(elem.sect->nconc);j++)
    {
      if(cYi[j]<=bounds[ibound+1]) /*TENSED.*/
      {
        bunsi-=0.0;
      }
      else if(cYj[j]>=bounds[ibound]) /*COMPRESSED.*/
      {
        bunsi-=cfcu*sBi[j]*sDi[j];
      }
      else /*ACROSSING.*/
      {
        bunsi-=cYi[j]*cBi[j]*cfcu;
        bunbo-=cBi[j]*cfcu;
      }
    }
    for(j=0;j<(elem.sect->nrein);j++)
    {
      if(rYi[j]<=bounds[ibound+1]) /*TENSED.*/
      {
        bunsi-=rftu*rAi[j];
      }
      else if(rYi[j]>=bounds[ibound]) /*COMPRESSED.*/
      {
        bunsi-=rfcu*rAi[j];
      }
    }
    yn=bunsi/bunbo;
  }

  if(yn<bounds[ibound+1] || bounds[ibound]<yn)
  {
    /*fprintf(fout0,"NEUTRAL AXIS OVERFLOW.%.3f<%.3f<%.3f[cm]\n",
            bounds[ibound+1],yn,bounds[ibound]);*/

    if(elem.sect->stype==RC)
    {
      Mu=0.9*rAt*rftu*(cYc-rYtg); /*CENTER STANDARD*/
      return Mu;
    }
    else return 0.0;
  }

  Yg=boundaryvalueultimate(elem,m,yn,sBi,sDi,cBi,cDi,
                                     sYi,sYj,cYi,cYj,
                                     sftu,sfcu,
                                     rftu,rfcu,
                                     0.0,cfcu,
                                     rAi,rYi,
                                     sN,sM,rN,rM,cN,cM,
                                     &Ng,&Mg);
  Mu=Mg;

  /*fprintf(fout0,"SRC軸力 srcN=%11.3f[tf]",-Ng/1000.0);
  fprintf(fout0," 重心:%11.3f",Yg);
  fprintf(fout0," 中立軸=%11.3f",yn);
  fprintf(fout0," 終局値 srcMu=%11.3f[tfm]\n",Mu/100000.0);*/

  /*fprintf(fout0,"%11.3f %11.3f 0.0",Mu/100000.0,-Nd/1000.0);*/

  return Mu;
}/*ultimatebendingofsrc*/

double boundaryvalue(struct element elem,struct materials m,
                     double Yn,
                     double sBi[],double sDi[],
                     double cBi[],double cDi[],
                     double sYi[],double sYj[],
                     double cYi[],double cYj[],
                     double sfi[],double sfj[],
                     double cfi[],double cfj[],
                     double rAi[],double rYi[],double rfi[],
                     double *gN,double *gM)
/*BOUNDARY VALUE gN,gM.*/
{
  int i;
  double gA,gY,Ygi,tA,tY;         /*Ygi:CENTER OF GRAVITY FOR RECT.*/
  double sNi,sMi,rNi,cNi;
  double sN,sM,rN,rM,cN,cM;

  gY=0.0; gA=0.0;
  for(i=0;i<(elem.sect->nsteel);i++)                       /*STEELS*/
  {
    gY+=(m.sE)*sBi[i]*sDi[i]*0.5*(sYi[i]+sYj[i]);
    gA+=(m.sE)*sBi[i]*sDi[i];
  }
  for(i=0;i<(elem.sect->nrein);i++)                /*REINFORCEMENTS*/
  {
    gY+=(m.rE)*rAi[i]*rYi[i];
    gA+=(m.rE)*rAi[i];
  }
  tA=gA; tY=gY;
  for(i=0;i<(elem.sect->nconc);i++)                     /*CONCRETES*/
  {
    tY+=(m.cE)*cBi[i]*cDi[i]*0.5*(cYi[i]+cYj[i]);
    tA+=(m.cE)*cBi[i]*cDi[i];

    if(cYj[i]>=Yn) /*RECT COMPRESSED.*/
    {
      gY+=(m.cE)*cBi[i]*cDi[i]*0.5*(cYi[i]+cYj[i]);
      gA+=(m.cE)*cBi[i]*cDi[i];
    }
    else if(cYi[i]>Yn) /*RECT ACROSSING NEUTRAL.*/
    {
      gY+=(m.cE)*cBi[i]*(cYi[i]-Yn)*0.5*(cYi[i]+Yn);
      gA+=(m.cE)*cBi[i]*(cYi[i]-Yn);
    }
  }
  gY/=gA; /*CENTER OF GRAVITY ON REMAINING SECTION.*/
  tY/=tA; /*CENTER OF GRAVITY ON WHOLE SECTION.*/

  gY=tY;  /*UNDER CONSTRUCTION.ALMOST CENTER OF ANALYSIS.*/

  sN=0.0; sM=0.0;
  for(i=0;i<(elem.sect->nsteel);i++)             /*STRESS OF STEELS*/
  {
    if(sYj[i]>=Yn || sYi[i]<=Yn)     /*ALL COMPRESSED OR ALL TENSED*/
    {
      Ygi=((2.0*sfj[i]+sfi[i])*sYj[i]+(2.0*sfi[i]+sfj[i])*sYi[i])
         /(3.0*(sfi[i]+sfj[i]));

      sNi=0.5*(sfi[i]+sfj[i])*sBi[i]*sDi[i];
      sMi=sNi*(gY-Ygi);
    }
    else /*ACROSSING*/
    {
      Ygi=Yn+2.0/3.0*(sYi[i]-Yn);
      sMi=0.5*sfi[i]*sBi[i]*(sYi[i]-Yn)*(gY-Ygi);

      Ygi=Yn+2.0/3.0*(sYj[i]-Yn);
      sMi+=-0.5*sfj[i]*sBi[i]*(sYj[i]-Yn)*(gY-Ygi);

      sNi=0.5*(sfi[i]+sfj[i])*sBi[i]*sDi[i];
    }
    sN+=sNi;
    sM+=sMi;
  }
  rN=0.0; rM=0.0;
  for(i=0;i<(elem.sect->nrein);i++)      /*STRESS OF REINFORCEMENTS*/
  {
    rNi=rfi[i]*rAi[i];
    rN+=rNi;
    rM+=rNi*(gY-rYi[i]);
  }
  cN=0.0; cM=0.0;
  for(i=0;i<(elem.sect->nconc);i++)           /*STRESS OF CONCRETES*/
  {
    if(cYj[i]>=Yn) /*COMPRESSED*/
    {
      Ygi=((2.0*cfj[i]+cfi[i])*cYj[i]+(2.0*cfi[i]+cfj[i])*cYi[i])
         /(3.0*(cfi[i]+cfj[i]));
      cNi=0.5*(cfi[i]+cfj[i])*cBi[i]*cDi[i];
    }
    else if(cYi[i]>Yn) /*ACROSSING*/
    {
      Ygi=Yn+2.0/3.0*(cYi[i]-Yn);
      cNi=0.5*cfi[i]*cBi[i]*(cYi[i]-Yn);
    }
    else /*TENSED*/
    {
      Ygi=0.0;
      cNi=0.0;
    }
    cN+=cNi;
    cM+=cNi*(gY-Ygi);
  }

  *gN=sN+rN+cN;
  *gM=sM+rM+cM;

  /*fprintf(fout0,"Yn=%.5f[cm]\n",Yn);*/
  /*fprintf(fout0,"gN=S%15.3f+R%15.3f+C%15.3f=%15.3f[kgf]\n",
          sN,rN,cN,*gN);
  fprintf(fout0,"gM=S%15.3f+R%15.3f+C%15.3f=%15.3f[kgfcm]\n",
          sM,rM,cM,*gM);*/

  return gY;
}/*boundaryvalue*/

double boundaryvalueultimate(struct element elem,struct materials m,
                             double Yn,
                             double sBi[],double sDi[],
                             double cBi[],double cDi[],
                             double sYi[],double sYj[],
                             double cYi[],double cYj[],
                             double sftu,double sfcu,
                             double rftu,double rfcu,
                             double cftu,double cfcu,
                             double rAi[],double rYi[],
                             double *sN,double *sM,
                             double *rN,double *rM,
                             double *cN,double *cM,
                             double *gN,double *gM)
/*BOUNDARY VALUE gN,gM.*/
{
  int i;
  double gA,gY,Ygi,tA,tY;         /*Ygi:CENTER OF GRAVITY FOR RECT.*/
  double sNi,sMi,rNi,cNi;

  gY=0.0; gA=0.0;
  for(i=0;i<(elem.sect->nsteel);i++)                       /*STEELS*/
  {
    gY+=sfcu*sBi[i]*sDi[i]*0.5*(sYi[i]+sYj[i]);
    gA+=sfcu*sBi[i]*sDi[i];
  }
  for(i=0;i<(elem.sect->nrein);i++)                /*REINFORCEMENTS*/
  {
    gY+=rfcu*rAi[i]*rYi[i];
    gA+=rfcu*rAi[i];
  }
  tA=gA; tY=gY;
  for(i=0;i<(elem.sect->nconc);i++)                     /*CONCRETES*/
  {
    tY+=cfcu*cBi[i]*cDi[i]*0.5*(cYi[i]+cYj[i]);
    tA+=cfcu*cBi[i]*cDi[i];

    if(cYj[i]>=Yn) /*RECT COMPRESSED.*/
    {
      gY+=cfcu*cBi[i]*cDi[i]*0.5*(cYi[i]+cYj[i]);
      gA+=cfcu*cBi[i]*cDi[i];
    }
    else if(cYi[i]>Yn) /*RECT ACROSSING NEUTRAL.*/
    {
      gY+=cfcu*cBi[i]*(cYi[i]-Yn)*0.5*(cYi[i]+Yn);
      gA+=cfcu*cBi[i]*(cYi[i]-Yn);
    }
  }
  gY/=gA; /*PLASTIC CENTER OF GRAVITY ON REMAINING SECTION.*/
  tY/=tA; /*PLASTIC CENTER OF GRAVITY ON WHOLE SECTION.*/

  gY=tY;  /*UNDER CONSTRUCTION.ALMOST CENTER OF ANALYSIS.*/

  *sN=0.0; *sM=0.0;
  for(i=0;i<(elem.sect->nsteel);i++)             /*STRESS OF STEELS*/
  {
    if(sYj[i]>=Yn)     /*ALL COMPRESSED*/
    {
      Ygi=0.5*(sYi[i]+sYj[i]);

      sNi=sfcu*sBi[i]*sDi[i];
      sMi=sNi*(gY-Ygi);
    }
    else if(sYi[i]<=Yn) /*ALL TENSED*/
    {
      Ygi=0.5*(sYi[i]+sYj[i]);

      sNi=sftu*sBi[i]*sDi[i];
      sMi=sNi*(gY-Ygi);
    }
    else /*ACROSSING*/
    {
      Ygi=0.5*(sYi[i]+Yn);
      sNi=sfcu*sBi[i]*(sYi[i]-Yn);
      sMi=sfcu*sBi[i]*(sYi[i]-Yn)*(gY-Ygi);

      Ygi=0.5*(sYj[i]+Yn);
      sNi+=sftu*sBi[i]*(Yn-sYi[i]);
      sMi+=sftu*sBi[i]*(Yn-sYj[i])*(gY-Ygi);
    }
    *sN+=sNi;
    *sM+=sMi;
  }
  *rN=0.0; *rM=0.0;
  for(i=0;i<(elem.sect->nrein);i++)      /*STRESS OF REINFORCEMENTS*/
  {
    if(rYi[i]>Yn)     /*ALL COMPRESSED*/
    {
      rNi=rfcu*rAi[i];
      *rN+=rNi;
      *rM+=rNi*(gY-rYi[i]);
    }
    else if(rYi[i]<Yn)     /*ALL TENSED*/
    {
      rNi=rftu*rAi[i];
      *rN+=rNi;
      *rM+=rNi*(gY-rYi[i]);
    }
  }
  *cN=0.0; *cM=0.0;
  for(i=0;i<(elem.sect->nconc);i++)           /*STRESS OF CONCRETES*/
  {
    if(cYj[i]>=Yn) /*COMPRESSED*/
    {
      Ygi=0.5*(cYi[i]+cYj[i]);
      cNi=cfcu*cBi[i]*cDi[i];
    }
    else if(cYi[i]>Yn) /*ACROSSING*/
    {
      Ygi=0.5*(cYi[i]+Yn);
      cNi=cfcu*cBi[i]*(cYi[i]-Yn);
    }
    else /*TENSED*/
    {
      Ygi=0.0;
      cNi=0.0;
    }
    *cN+=cNi;
    *cM+=cNi*(gY-Ygi);
  }

  *gN=*sN+*rN+*cN;
  *gM=*sM+*rM+*cM;

  /*fprintf(fout0,"Yn=%.5f[cm] ",Yn);*/
  /*fprintf(fout0,"gN=S%15.3f+R%15.3f+C%15.3f=%15.3f[kgf]\n",
          *sN,*rN,*cN,*gN);*/
  /*fprintf(fout0,"gM=S%15.3f+R%15.3f+C%15.3f=%15.3f[kgfcm]\n",
          *sM,*rM,*cM,*gM);*/

  return gY;
}/*boundaryvalueultimate*/

double allowableshearofrclong(struct element elem,
                              struct materials m,
                              int axis,                   /*0:x 1:y*/
                              double Qd,                    /*[kgf]*/
                              double Md)                  /*[kgfcm]*/
/*RETURN:ALLOWABLE SHEAR OF RC COLUMN,GIRDER FOR LONG.*/
{
  int i;
  double cBi[MAXCRECT],cDi[MAXCRECT];    /*WIDTH,DEPTH OF CONCRETES*/
  double cYi[MAXCRECT],cYj[MAXCRECT];
  double rYi[MAXREINS];
  double rYt,cYt,cYc;
  double wft,cfs; /*ALLOWABLE STRESSES*/
  double pw;

  double dc,ac,Ac,eAc;
  double Qa;

  rYt=1000.0; /*[cm]*/
  for(i=0;i<(elem.sect->nrein);i++)
  {
    if(axis==SX) rYi[i]=elem.sect->rein[i].y;
    if(axis==SY) rYi[i]=elem.sect->rein[i].x;

    if(rYt>rYi[i]) rYt=rYi[i];
  }

  Ac=0.0; /*CONCRETE AREA TOTAL[cm2]*/
  eAc=0.0; /*EFFECTIVE CONCRETE AREA=ABOVE BOTTOM REIN.*/
  cYc=-1000.0; cYt=1000.0; /*[cm]*/
  for(i=0;i<(elem.sect->nconc);i++)
  {
    if(axis==SX)
    {
      cBi[i]=fabs(elem.sect->crect[i].right
                 -elem.sect->crect[i].left);
      cYi[i]=elem.sect->crect[i].top;
      cYj[i]=elem.sect->crect[i].bottom;
    }
    if(axis==SY)
    {
      cBi[i]=fabs(elem.sect->crect[i].top
                 -elem.sect->crect[i].bottom);
      cYi[i]=elem.sect->crect[i].right;
      cYj[i]=elem.sect->crect[i].left;
    }
    cDi[i]=cYi[i]-cYj[i];

    Ac+=(cBi[i]*cDi[i]); /*[cm2]*/

    if(cYc<cYi[i]) cYc=cYi[i];
    if(cYt>cYj[i]) cYt=cYj[i];

    if(cYj[i]>=rYt)     eAc+=(cBi[i]*cDi[i]);       /*[cm2]*/
    else if(cYi[i]>rYt) eAc+=(cBi[i]*(cYi[i]-rYt)); /*[cm2]*/
  }

  dc=fabs(rYt-cYc);

  ac=4.0/(fabs(Md/Qd/dc)+1.0);
  if(ac<1.0) ac=1.0;
  if(ac>2.0) ac=2.0;

  cfs=fabs(m.cfs);
  wft=fabs(m.wft);

  pw=elem.sect->shearrein[1-axis];
  if(pw>0.012) pw=0.012;

  if(elem.sect->etype==COLUMN) /*"RC STANDARD"16.(25)*/
  {
    Qa=7.0/8.0*ac*eAc*cfs;
  }
  else if(elem.sect->etype==GIRDER || elem.sect->etype==BEAM)
  {
    Qa=7.0/8.0*eAc*(ac*cfs+0.5*wft*(pw-0.002)); /*16.(22)*/
  }
  else Qa=0.0;

  /*
  fprintf(fout0,"Ac=%.3f cYc=%.3f cYt=%.3f[cm2,cm]\n",Ac,cYc,cYt);
  fprintf(fout0,"rYt=%.3f[cm] pw=%.3f\n",rYt,pw);
  fprintf(fout0,"rcQa=%11.3f[tf]",fabs(Qa)/1000.0);
  */

  return Qa;
}/*allowableshearofrclong*/

double allowableshearofrcshort(struct element elem,
                               struct materials m,
                               int axis, /*0:x 1:y*/
                               double Qd, /*[kgf]*/
                               double Md)
/*RETURN:ALLOWABLE SHEAR OF RC COLUMN,GIRDER FOR SHORT.*/
{
  int i;
  double cBi[MAXCRECT],cDi[MAXCRECT];    /*WIDTH,DEPTH OF CONCRETES*/
  double cYi[MAXCRECT],cYj[MAXCRECT];
  double rYi[MAXREINS];
  double rYt,cYt,cYc;
  double cfs,wft; /*ALLOWABLE STRESSES*/
  double pw;

  double dc,ac,Ac,eAc;
  double Qa;

  rYt=1000.0; /*[cm]*/
  for(i=0;i<(elem.sect->nrein);i++)
  {
    if(axis==SX) rYi[i]=elem.sect->rein[i].y;
    if(axis==SY) rYi[i]=elem.sect->rein[i].x;

    if(rYt>rYi[i]) rYt=rYi[i];
  }

  Ac=0.0; /*CONCRETE AREA TOTAL[cm2]*/
  eAc=0.0; /*EFFECTIVE CONCRETE AREA=ABOVE BOTTOM REIN.*/
  cYc=-1000.0; cYt=1000.0; /*[cm]*/
  for(i=0;i<(elem.sect->nconc);i++)
  {
    if(axis==SX)
    {
      cBi[i]=fabs(elem.sect->crect[i].right
                 -elem.sect->crect[i].left);
      cYi[i]=elem.sect->crect[i].top;
      cYj[i]=elem.sect->crect[i].bottom;
    }
    if(axis==SY)
    {
      cBi[i]=fabs(elem.sect->crect[i].top
                 -elem.sect->crect[i].bottom);
      cYi[i]=elem.sect->crect[i].right;
      cYj[i]=elem.sect->crect[i].left;
    }
    cDi[i]=cYi[i]-cYj[i];

    Ac+=(cBi[i]*cDi[i]); /*[cm2]*/

    if(cYc<cYi[i]) cYc=cYi[i];
    if(cYt>cYj[i]) cYt=cYj[i];

    if(cYj[i]>=rYt)     eAc+=(cBi[i]*cDi[i]);       /*[cm2]*/
    else if(cYi[i]>rYt) eAc+=(cBi[i]*(cYi[i]-rYt)); /*[cm2]*/
  }

  dc=fabs(rYt-cYc);

  ac=4.0/(fabs(Md/Qd/dc)+1.0);
  if(ac<1.0) ac=1.0;
  if(ac>2.0) ac=2.0;

  cfs=fabs(m.cfs);
  wft=fabs(m.wft);

  pw=elem.sect->shearrein[1-axis];
  if(pw>0.012) pw=0.012;

  /*fprintf(fout0,"dc=%.3f cfs=%.3f wft=%.3f pw=%.5f\n",
          dc,cfs,wft,pw);*/

  if(elem.sect->etype==COLUMN)               /*"RC STANDARD"16.(25)*/
  {
    Qa=7.0/8.0*eAc*(cfs+0.5*wft*(pw-0.002));
  }
  else if(elem.sect->etype==GIRDER || elem.sect->etype==BEAM)
  {
    Qa=7.0/8.0*eAc*(ac*cfs+0.5*wft*(pw-0.002)); /*16.(22)*/
  }
  else Qa=0.0;

  /*
  fprintf(fout0,"Ac=%.3f cYc=%.3f cYt=%.3f[cm2,cm]\n",Ac,cYc,cYt);
  fprintf(fout0,"rYt=%.3f[cm] pw=%.3f\n",rYt,pw);
  fprintf(fout0,"rcQa=%11.3f[tf]",fabs(Qa)/1000.0);
  */

  return Qa;
}/*allowableshearofrcshort*/

double ultimateshearofrc(struct element elem,
                         struct materials m,
                         int axis,
                         double Nd,double Qd,double Md)
/*ULTIMATE SHEAR OF RC COLUMN.*/
{
  int i;
  double at,pt,pw;
  double s0; /*AXIAL STRESS.*/
  double cBi[MAXCRECT],cDi[MAXCRECT];    /*WIDTH,DEPTH OF CONCRETES*/
  double cYi[MAXCRECT],cYj[MAXCRECT];
  double rYi[MAXREINS];
  double rYt,cYt,cYc;
  double Fc,wfp; /*PLASTIC STRESSES*/
  double dc,ac,Ac,eAc;
  double gA,Yg,dA;
  double bQu,rcQu;

  gA=0.0; Yg=0.0;

  rYt=1000.0; /*[cm]*/
  for(i=0;i<(elem.sect->nrein);i++)
  {
    if(axis==0) rYi[i]=elem.sect->rein[i].y;
    if(axis==1) rYi[i]=elem.sect->rein[i].x;

    dA=(elem.sect->rein[i].area)*(m.rE);
    gA+=dA;
    Yg+=dA*rYi[i];

    if(rYt>rYi[i]) rYt=rYi[i];
  }

  Ac=0.0; /*CONCRETE AREA TOTAL[cm2]*/
  eAc=0.0; /*EFFECTIVE CONCRETE AREA=ABOVE BOTTOM REIN.*/
  cYc=-1000.0; cYt=1000.0; /*[cm]*/
  for(i=0;i<(elem.sect->nconc);i++)
  {
    if(axis==0)
    {
      cBi[i]=fabs(elem.sect->crect[i].right
                 -elem.sect->crect[i].left);
      cYi[i]=elem.sect->crect[i].top;
      cYj[i]=elem.sect->crect[i].bottom;
    }
    if(axis==1)
    {
      cBi[i]=fabs(elem.sect->crect[i].top
                 -elem.sect->crect[i].bottom);
      cYi[i]=elem.sect->crect[i].right;
      cYj[i]=elem.sect->crect[i].left;
    }
    cDi[i]=cYi[i]-cYj[i];

    dA=(cBi[i]*cDi[i])*(m.cE);
    gA+=dA;
    Yg+=dA*0.5*(cYi[i]+cYj[i]);

    Ac+=(cBi[i]*cDi[i]); /*[cm2]*/

    if(cYc<cYi[i]) cYc=cYi[i];
    if(cYt>cYj[i]) cYt=cYj[i];

    if(cYj[i]>=rYt)     eAc+=(cBi[i]*cDi[i]);       /*[cm2]*/
    else if(cYi[i]>rYt) eAc+=(cBi[i]*(cYi[i]-rYt)); /*[cm2]*/
  }

  dc=fabs(rYt-cYc);

  s0=Nd/Ac; /*AXIAL STRESS[kgf/cm2] +:COMPRESSION.*/
  ac=fabs(Md/Qd/dc);
  if(ac<1.0) ac=1.0;
  if(ac>3.0) ac=3.0;

  Yg/=gA; /*CENTER OF GRAVITY.*/
  at=0.0;
  for(i=0;i<(elem.sect->nrein);i++)
  {
    if(rYi[i]<Yg) at+=(elem.sect->rein[i].area); /*TENSED REIN.*/
  }
  pt=at/eAc*100.0; /*RATE OF TENSED REIN [%].*/

  Fc=fabs(m.Fc);
  wfp=fabs(m.wfp);

  pw=elem.sect->shearrein[1-axis]; /*RATE OF HOOP [RATE].*/
  if(pw>0.012) pw=0.012;

  if((0.9+s0/250.0)<=0.0) return 0.0;

  bQu=((0.068*pow(pt,0.23)*(Fc+180.0))/(ac+0.12)+2.7*sqrt(pw*wfp))
     *7.0/8.0*eAc;                    /*"CENTER STANDARD"(付1.7.2b)*/
  rcQu=(0.9+s0/250.0)*bQu;            /*"CENTER STANDARD"(付1.7.4b)*/

  /*fprintf(fout0,"Nd=%.3f Qd=%.3f Md=%.3f\n",Nd,Qd,Md);
  fprintf(fout0,"d=%.3f[cm] pw=%.5f\n",dc,pw);
  fprintf(fout0,"rcQu=%11.3f[tf]\n",rcQu/1000.0);*/

  return rcQu; /*[kgf]*/
}/*ultimateshearofrc*/

double allowableshearofsrclong(struct element elem,
                               struct materials m,
                               int axis,                  /*0:x 1:y*/
                               int haxis,  /*H AXIS 0:STRONG 1:WEAK*/
                               double srcQd,                /*[kgf]*/
                               double srcMd)              /*[kgfcm]*/
/*RETURN:ALLOWABLE SHEAR OF SRC COLUMN WITH H FOR LONG.*/
{
  int i,ii;
  double eBrate,rate,sBi,sYi,sYj;
  double cBi[MAXCRECT],cDi[MAXCRECT];    /*WIDTH,DEPTH OF CONCRETES*/
  double cYi[MAXCRECT],cYj[MAXCRECT];
  double rYi[MAXREINS];
  double rYt,cYt,cYc;
  double wft,cfs; /*ALLOWABLE STRESSES*/
  /*double pw;*/

  double dc,ac,Ac,eAc;
  double sQa,rcQa,srcQa;

  rYt=1000.0; /*[cm]*/
  for(i=0;i<(elem.sect->nrein);i++)
  {
    if(axis==SX) rYi[i]=elem.sect->rein[i].y;
    if(axis==SY) rYi[i]=elem.sect->rein[i].x;

    if(rYt>rYi[i]) rYt=rYi[i];
  }

  Ac=0.0; /*CONCRETE AREA TOTAL[cm2]*/
  eAc=0.0; /*EFFECTIVE CONCRETE AREA=ABOVE BOTTOM REIN.*/
  cYc=-1000.0; cYt=1000.0; /*[cm]*/
  for(i=0;i<(elem.sect->nconc);i++)
  {
    if(axis==SX)
    {
      cBi[i]=fabs(elem.sect->crect[i].right
                 -elem.sect->crect[i].left);
      cYi[i]=elem.sect->crect[i].top;
      cYj[i]=elem.sect->crect[i].bottom;
    }
    if(axis==SY)
    {
      cBi[i]=fabs(elem.sect->crect[i].top
                 -elem.sect->crect[i].bottom);
      cYi[i]=elem.sect->crect[i].right;
      cYj[i]=elem.sect->crect[i].left;
    }
    cDi[i]=cYi[i]-cYj[i];

    Ac+=(cBi[i]*cDi[i]); /*[cm2]*/

    if(cYc<cYi[i]) cYc=cYi[i];
    if(cYt>cYj[i]) cYt=cYj[i];

    if(cYj[i]>=rYt)     eAc+=(cBi[i]*cDi[i]);       /*[cm2]*/
    else if(cYi[i]>rYt) eAc+=(cBi[i]*(cYi[i]-rYt)); /*[cm2]*/
  }

  eBrate=1.0;
  for(i=0;i<(elem.sect->nsteel);i++)
  {
    if(axis==SX)
    {
      sBi=fabs(elem.sect->srect[i].right
              -elem.sect->srect[i].left);
      sYi=elem.sect->srect[i].top;
      sYj=elem.sect->srect[i].bottom;
    }
    if(axis==SY)
    {
      sBi=fabs(elem.sect->srect[i].top
              -elem.sect->srect[i].bottom);
      sYi=elem.sect->srect[i].right;
      sYj=elem.sect->srect[i].left;
    }

    for(ii=0;ii<(elem.sect->nconc);ii++)
    {
      if(cYi[ii]>sYj && cYj[ii]<sYi)
      {
        rate=(cBi[ii]-sBi)/cBi[ii];          /*UNDER CONSIDERATION.*/
                                      /*INCORRECT IF cB IS DIVIDED.*/
        if(rate<eBrate) eBrate=rate; /*EFFECTIVE WIDTH OF CONCRETE.*/
      }
    }
  }

  dc=fabs(cYc-rYt); /*[cm]*/

  eBrate*=3.0; /*"SRC STANDARD"(59)*/
  ac=4.0/(fabs(srcMd/srcQd/dc)+1.0); /*"SRC STANDARD"(53)*/
  if(ac<1.0) ac=1.0;
  if(ac>2.0) ac=2.0;
  if(ac>eBrate) ac=eBrate;

  cfs=fabs(m.cfs);
  wft=fabs(m.wft);

  /*pw=elem.sect->shearrein[1-axis];
  if(pw>0.006) pw=0.006;*/

  rcQa=7.0/8.0*ac*eAc*cfs; /*"SRC STANDARD"(58)*/

  if(haxis==HSTRONG)
  {
    sQa=allowultimshearofhstrong(PLONG,m.sF,*(elem.sect),axis);
  }
  if(haxis==HWEAK)
  {
    sQa=allowultimshearofhweak(PLONG,m.sF,*(elem.sect),axis);
  }
  srcQa=sQa+rcQa;

  /*
  fprintf(fout0,"Ac=%.3f cYc=%.3f cYt=%.3f[cm2,cm]\n",Ac,cYc,cYt);
  fprintf(fout0,"rYt=%.3f[cm]\n",rYt);
  fprintf(fout0,"sQa=%.3f[tf] rcQa=%.3f[tf]\n",
          sQa/1000.0,rcQa/1000.0);
  fprintf(fout0,"srcQa=%11.3f[tf]",srcQa/1000.0);
  */

  return srcQa;
}/*allowableshearofsrclong*/

double allowableshearofsrcshort(struct element elem,
                                struct materials m,
                                int axis,                 /*0:x 1:y*/
                                int haxis, /*H AXIS 0:STRONG 1:WEAK*/
                                double srcQd,               /*[kgf]*/
                                double srcMd)             /*[kgfcm]*/
/*RETURN:ALLOWABLE SHEAR OF SRC COLUMN WITH H FOR SHORT.*/
{
  int i,ii;
  double eBrate,rate,sBi,sYi,sYj;
  double cBi[MAXCRECT],cDi[MAXCRECT];    /*WIDTH,DEPTH OF CONCRETES*/
  double cYi[MAXCRECT],cYj[MAXCRECT];
  double rYi[MAXREINS];
  double rYt,cYt,cYc;
  double wft,cfs; /*ALLOWABLE STRESSES*/
  double pw;

  double dc,Ac,eAc;
  /*double ac;*/
  double sQa,rcQa,rcQa1,rcQa2,srcQa;

  rYt=1000.0; /*[cm]*/
  for(i=0;i<(elem.sect->nrein);i++)
  {
    if(axis==SX) rYi[i]=elem.sect->rein[i].y;
    if(axis==SY) rYi[i]=elem.sect->rein[i].x;

    if(rYt>rYi[i]) rYt=rYi[i];
  }

  Ac=0.0; /*CONCRETE AREA TOTAL[cm2]*/
  eAc=0.0; /*EFFECTIVE CONCRETE AREA=ABOVE BOTTOM REIN.*/
  cYc=-1000.0; cYt=1000.0; /*[cm]*/
  for(i=0;i<(elem.sect->nconc);i++)
  {
    if(axis==SX)
    {
      cBi[i]=fabs(elem.sect->crect[i].right
                 -elem.sect->crect[i].left);
      cYi[i]=elem.sect->crect[i].top;
      cYj[i]=elem.sect->crect[i].bottom;
    }
    if(axis==SY)
    {
      cBi[i]=fabs(elem.sect->crect[i].top
                 -elem.sect->crect[i].bottom);
      cYi[i]=elem.sect->crect[i].right;
      cYj[i]=elem.sect->crect[i].left;
    }
    cDi[i]=cYi[i]-cYj[i];

    Ac+=(cBi[i]*cDi[i]); /*[cm2]*/

    if(cYc<cYi[i]) cYc=cYi[i];
    if(cYt>cYj[i]) cYt=cYj[i];

    if(cYj[i]>=rYt)     eAc+=(cBi[i]*cDi[i]);       /*[cm2]*/
    else if(cYi[i]>rYt) eAc+=(cBi[i]*(cYi[i]-rYt)); /*[cm2]*/
  }

  eBrate=1.0;
  for(i=0;i<(elem.sect->nsteel);i++)
  {
    if(axis==SX)
    {
      sBi=fabs(elem.sect->srect[i].right
              -elem.sect->srect[i].left);
      sYi=elem.sect->srect[i].top;
      sYj=elem.sect->srect[i].bottom;
    }
    if(axis==SY)
    {
      sBi=fabs(elem.sect->srect[i].top
              -elem.sect->srect[i].bottom);
      sYi=elem.sect->srect[i].right;
      sYj=elem.sect->srect[i].left;
    }

    for(ii=0;ii<(elem.sect->nconc);ii++)
    {
      if(cYi[ii]>sYj && cYj[ii]<sYi)
      {
        rate=(cBi[ii]-sBi)/cBi[ii];          /*UNDER CONSIDERATION.*/
                                      /*INCORRECT IF cB IS DIVIDED.*/
        if(rate<eBrate) eBrate=rate; /*EFFECTIVE WIDTH OF CONCRETE.*/
      }
    }
  }

  dc=fabs(cYc-rYt); /*[cm]*/

  /*ac=4.0/(fabs(srcMd/srcQd/dc)+1.0);
  if(ac<1.0) ac=1.0;
  if(ac>2.0) ac=2.0;*/

  cfs=fabs(m.cfs);
  wft=fabs(m.wft);

  pw=elem.sect->shearrein[1-axis];
  if(pw>0.006) pw=0.006;

  rcQa1=7.0/8.0*eAc*(cfs+0.5*pw*wft); /*"SRC STANDARD"(63)*/
  rcQa2=7.0/8.0*eAc*(2.0*eBrate*cfs+pw*wft);

  if(rcQa1<=rcQa2) rcQa=rcQa1; /*"SRC STANDARD"(62)*/
  else             rcQa=rcQa2;

  if(haxis==HSTRONG)
  {
    sQa=allowultimshearofhstrong(PSHORT,m.sF,*(elem.sect),axis);
  }
  if(haxis==HWEAK)
  {
    sQa=allowultimshearofhweak(PSHORT,m.sF,*(elem.sect),axis);
  }
  srcQa=sQa+rcQa;

  /*
  fprintf(fout0,"Ac=%.3f cYc=%.3f cYt=%.3f[cm2,cm]\n",Ac,cYc,cYt);
  fprintf(fout0,"rYt=%.3f[cm] pw=%.3f\n",rYt,pw);
  fprintf(fout0,"sQa=%.3f[tf] rcQa=%.3f[tf]\n",
          sQa/1000.0,rcQa/1000.0);
  fprintf(fout0,"srcQa=%11.3f[tf]",srcQa/1000.0);
  */

  return srcQa;
}/*allowableshearofsrcshort*/

double ultimateshearofsrc(struct element elem,
                          struct materials m,
                          int axis,                       /*0:x 1:y*/
                          int haxis,       /*H AXIS 0:STRONG 1:WEAK*/
                          double rcQd,                      /*[kgf]*/
                          double rcMd)                    /*[kgfcm]*/
/*RETURN:ULTIMATE SHEAR OF SRC COLUMN WITH H.*/
{
  int i,ii;
  double eBrate,rate,sBi,sYi,sYj;
  double cBi[MAXCRECT],cDi[MAXCRECT];    /*WIDTH,DEPTH OF CONCRETES*/
  double cYi[MAXCRECT],cYj[MAXCRECT];
  double rYi[MAXREINS];
  double rYt,cYt,cYc;
  double wfp,cfsu,cfs1,cfs2; /*ULTIMATE STRESSES*/
  double pw;

  double dc,ac,Ac,eAc;
  double sQu,sQu1,rcQu,rcQu1,rcQu2,srcQu;

  rYt=1000.0; /*[cm]*/
  for(i=0;i<(elem.sect->nrein);i++)
  {
    if(axis==SX) rYi[i]=elem.sect->rein[i].y;
    if(axis==SY) rYi[i]=elem.sect->rein[i].x;

    if(rYt>rYi[i]) rYt=rYi[i];
  }

  Ac=0.0; /*CONCRETE AREA TOTAL[cm2]*/
  eAc=0.0; /*EFFECTIVE CONCRETE AREA=ABOVE BOTTOM REIN.*/
  cYc=-1000.0; cYt=1000.0; /*[cm]*/
  for(i=0;i<(elem.sect->nconc);i++)
  {
    if(axis==SX)
    {
      cBi[i]=fabs(elem.sect->crect[i].right
                 -elem.sect->crect[i].left);
      cYi[i]=elem.sect->crect[i].top;
      cYj[i]=elem.sect->crect[i].bottom;
    }
    if(axis==SY)
    {
      cBi[i]=fabs(elem.sect->crect[i].top
                 -elem.sect->crect[i].bottom);
      cYi[i]=elem.sect->crect[i].right;
      cYj[i]=elem.sect->crect[i].left;
    }
    cDi[i]=cYi[i]-cYj[i];

    Ac+=(cBi[i]*cDi[i]); /*[cm2]*/

    if(cYc<cYi[i]) cYc=cYi[i];
    if(cYt>cYj[i]) cYt=cYj[i];

    if(cYj[i]>=rYt)     eAc+=(cBi[i]*cDi[i]);       /*[cm2]*/
    else if(cYi[i]>rYt) eAc+=(cBi[i]*(cYi[i]-rYt)); /*[cm2]*/
  }

  eBrate=1.0;
  for(i=0;i<(elem.sect->nsteel);i++)
  {
    if(axis==SX)
    {
      sBi=fabs(elem.sect->srect[i].right
              -elem.sect->srect[i].left);
      sYi=elem.sect->srect[i].top;
      sYj=elem.sect->srect[i].bottom;
    }
    if(axis==SY)
    {
      sBi=fabs(elem.sect->srect[i].top
              -elem.sect->srect[i].bottom);
      sYi=elem.sect->srect[i].right;
      sYj=elem.sect->srect[i].left;
    }

    for(ii=0;ii<(elem.sect->nconc);ii++)
    {
      if(cYi[ii]>sYj && cYj[ii]<sYi)
      {
        rate=(cBi[ii]-sBi)/cBi[ii];          /*UNDER CONSIDERATION.*/
                                      /*INCORRECT IF cB IS DIVIDED.*/
        if(rate<eBrate) eBrate=rate; /*EFFECTIVE WIDTH OF CONCRETE.*/
      }
    }
  }

  dc=fabs(cYc-rYt); /*[cm]*/

  ac=4.0/(fabs(rcMd/rcQd/dc)+1.0); /*"SRC STANDARD"(45)*/
  if(ac<1.0) ac=1.0;
  if(ac>2.0) ac=2.0;

  cfs1=0.15*fabs(m.Fc);              /*"SRC STANDARD"(124)[kgf/cm2]*/
  cfs2=22.5+4.5*fabs(m.Fc)/100;
  if(cfs1<=cfs2) cfsu=cfs1;
  else           cfsu=cfs2;

  wfp=fabs(m.wfp);

  pw=elem.sect->shearrein[1-axis];
  if(pw>0.006) pw=0.006;
  if(pw<0.001) fprintf(fout0,"HOOP TOO LITTLE.Pw=%.7f < 0.001\n",pw);

  /*rcQu0=2*rcMu/inner;*/                     /*"SRC STANDARD"(121)*/

  rcQu1=7.0/8.0*eAc*(0.5*ac*cfsu+0.5*pw*wfp); /*"SRC STANDARD"(123)*/
  rcQu2=7.0/8.0*eAc*(eBrate*cfsu+pw*wfp);

  if(rcQu1<=rcQu2) rcQu=rcQu1;
  else             rcQu=rcQu2;

  /*sQu0=2*sMu/inner;*/                       /*"SRC STANDARD"(126)*/

  if(haxis==HSTRONG)
  {
    sQu1=allowultimshearofhstrong(PULTIMATE,m.sF,*(elem.sect),axis);
  }
  if(haxis==HWEAK)
  {
    sQu1=allowultimshearofhweak(PULTIMATE,m.sF,*(elem.sect),axis);
  }
  sQu=sQu1;

  srcQu=sQu+rcQu;                             /*"SRC STANDARD"(119)*/

  /*fprintf(fout0,"Ac=%.3f cYc=%.3f cYt=%.3f[cm2,cm]\n",Ac,cYc,cYt);
  fprintf(fout0,"rYt=%.3f[cm] pw=%.5f\n",rYt,pw);
  fprintf(fout0,"sMu=%.3f rcMu=%.3f[tfm]\n",
          sMu/100000.0,rcMu/100000.0);
  fprintf(fout0,"sQu0=%.3f sQu1=%.3f[tf]\n",
          sQu0/1000.0,sQu1/1000.0);
  fprintf(fout0,"rcQu0=%.3f rcQu1=%.3f rcQu2=%.3f[tf]\n",
          rcQu0/1000.0,rcQu1/1000.0,rcQu2/1000.0);
  fprintf(fout0,"sQu=%.3f[tf] rcQu=%.3f[tf]\n",
          sQu/1000.0,rcQu/1000.0);
  fprintf(fout0,"srcQu=%11.3f[tf]\n",srcQu/1000.0);*/

  return srcQu;
}/*ultimateshearofsrc*/

double allowultimshearofrcwall(struct element elem,
                               struct materials m,
                               double l, /*LENGTH[cm]*/
                               int period)
/*ALLOWABLE,ULTIMATE SHEAR OF RC WALL.*/
{
  char code[20];
  double t,r;                           /*t:THICKNESS r:WINDOW RATE*/
  double ps;
  double cfs,wft;
  double Qc,Qw,Qau;

  if(period==PULTIMATE)
  {
    cfs=2.4*sqrt(fabs(m.Fc)); /*"RC STANDARD" APP21.(20)*/
    wft=fabs(m.wfp);
  }
  else
  {
    cfs=fabs(m.cfs);
    wft=fabs(m.wft);
  }

  t=elem.sect->thick;
  ps=elem.sect->shearrein[0];
  r=1.0-elem.sect->windowrate;

  /*fprintf(fout0,"THICK=%.3f Ps=%.3f WINDOW=%.3f",t,ps,r);
  fprintf(fout0," cfs=%.3f wft=%.3f\n",cfs,wft);*/

  Qc=r*t*l*cfs;                             /*"RC STANDARD" 18.(30)*/
  Qw=r*ps*t*l*wft;                          /*"RC STANDARD" 18.(32)*/

  if(Qc>=Qw)
  {
    Qau=Qc;
    strcpy(code,"CONCRETE");
  }
  else
  {
    Qau=Qw;
    strcpy(code,"REINFORCE");
  }
  fprintf(fout0,"CONCRETE:Qc=%7.3f[tf]",Qc/1000.0);
  fprintf(fout0," REINFORCE:Qw=%7.3f[tf]",Qw/1000.0);

  if(period==PULTIMATE) fprintf(fout0," Qu=");
  else                 fprintf(fout0," Qa=");

  fprintf(fout0,"%7.3f[tf] BY %s.\n",Qau/1000.0,code);

  return Qau;
}/*allowableshearofrcwall*/

double allowabletensionofs(double sF,
                           struct section sect)
/*ALLOWABLE TENSION OF STEEL LONG. =SHORT/1.5*/
{
  double A;
  double sfta,sNta;                 /*s:STEEL t:TENSION a:ALLOWABLE*/

  sfta=sF/1.5;                                  /*"S STANDARD"(5.1)*/

  A=coeffA(sect.nsteel,sect.srect);                         /*[cm2]*/
  sNta=-sfta*A; /*-:TENSION*/                               /*[kgf]*/

  /*fprintf(fout0,"sft=%.3f[kgf/cm2] As=%.3f[cm2]\n",sfta,A);*/

  return sNta;
}/*allowabletensionofs*/

double allowablecompressionofs(double E,double F,
                               double lk, /*BUCKLING LENGTH[cm]*/
                               struct section sect,int axis)
/*ALLOWABLE COMPRESSION OF STEEL LONG. =SHORT/1.5*/
{
  double i,A;
  double lam,LAM,nyu;
  double sfca,sNca;             /*s:STEEL c:COMPRESSION a:ALLOWABLE*/

  LAM=sqrt(PI*PI*E/0.6/F);                      /*"S STANDARD"(5.5)*/

  i=coeffi(sect.nsteel,sect.srect,axis);

  lam=lk/i;                                    /*"S STANDARD"(11.1)*/
  if((sect.etype==COLUMN && lam>200)||(lam>250))
  {
    fprintf(fout0,"ELEMENT TOO LONG.Lk/i=%.5f\n",lam);
    /*return 0.0;*/
  }

  nyu=1.5+(lam/LAM)*(lam/LAM)/1.5;

  if(lam<=LAM)
  {
    sfca=(1.0-0.4*(lam/LAM)*(lam/LAM))*F/nyu;   /*"S STANDARD"(5.3)*/
  }
  else
  {
    sfca=0.277*F/(lam/LAM)/(lam/LAM);           /*"S STANDARD"(5.4)*/
  }

  A=coeffA(sect.nsteel,sect.srect);                         /*[cm2]*/
  sNca=sfca*A; /*+:COMPRESSION*/                            /*[kgf]*/

  /*
  fprintf(fout0,"Lk=%.3f[cm] i=%.3f[cm] lam=%.3f LAM=%.3f NYU=%.3f\n",
          lk,i,lam,LAM,nyu);
  fprintf(fout0,"sfc=%.3f[kgf/cm2] As=%.3f[cm2]\n",sfca,A);
  */

  return sNca;
}/*allowablecompressionofs*/

double allowablebendingofhstrong(int period, /*0:LONG 1:SHORT*/
                                 double Nd, /*[kgf]*/
                                 double lb, /*BUCKLING LENGTH [cm]*/
                                 struct materials m,
                                 struct section sect,int axis)
/*ALLOWABLE BENDING OF H,[ IN STRONG AXIS.*/
/*AXIAL FORCE +:COMPRESSION -:TENSION [tf]*/
{
  double H,B,tf;
  double Af,Z;
  double sfta,sfba;       /*s:STEEL t:TENSION b:BENDING a:ALLOWABLE*/
  double sNa=0.0,sMa;

  if(sect.nsteel!=3)
  {
    fprintf(fout0,"STEEL NOT H,[.\n");
    return 0.0;
  }

  if(Nd>0.0) sNa=allowablecompressionofs(m.sE,m.sF,lb,sect,
                                         (1-axis)); /*FOR WEAK.*/
  if(Nd<0.0) sNa=allowabletensionofs(m.sF,sect);

  if(period==PSHORT) sNa*=1.5;

  if((Nd/sNa)>=1.0) return 0.0;

  if(axis==SX)
  {
    H =(sect.srect[0].top-sect.srect[2].bottom);             /*[cm]*/
    B =(sect.srect[0].right-sect.srect[0].left);             /*[cm]*/
    tf=(sect.srect[0].top-sect.srect[0].bottom);             /*[cm]*/
  }
  if(axis==SY)
  {
    H =(sect.srect[2].right-sect.srect[0].left);             /*[cm]*/
    B =(sect.srect[0].top-sect.srect[0].bottom);             /*[cm]*/
    tf=(sect.srect[0].right-sect.srect[0].left);             /*[cm]*/
  }

  Af=B*tf;

  sfta=m.sF; /*[kgf/cm2]*/
  if(period==PLONG) sfta/=1.5;

  sfba=900.0*Af/lb/H*1000.0;          /*[kgf/cm2] "S STANDARD"(5.8)*/
  if(period==PSHORT) sfba*=1.5;
  if(sfba>sfta) sfba=sfta;

  Z=coeffZ(sect.nsteel,sect.srect,axis);                    /*[cm3]*/

  sMa=(1.0-Nd/sNa)*sfba*Z;                                /*[kgfcm]*/

  /*
  fprintf(fout0,"STRONG. Nd=%.3f[kgf] sNa=%.3f[kgf] Nd/sNa=%.3f\n",
          Nd,sNa,(Nd/sNa));
  fprintf(fout0,"H=%.1f[cm] B=%.1f[cm] tf=%.1f[cm]\n",H,B,tf);
  fprintf(fout0,"sft=%.3f[kgf/cm2] sfb=%.3f[kgf/cm2] lb=%.3f[cm]",
          sfta,sfba,lb);
  fprintf(fout0," Z=%.3f[cm3]\n",Z);
  */

  return sMa;
}/*allowablebendingofhstrong*/

double allowablebendingofsweak(int period, /*0:LONG 1:SHORT*/
                               double Nd,
                               double lk, /*BUCKLING LENGTH [cm]*/
                               struct materials m,
                               struct section sect,int axis)
/*RETURN:ALLOWABLE BENDING [kgfcm] OF S WEAK,CLOSED.*/
{
  double Ze;
  double sfb;
  double sNa=0.0,sMa;                                 /*a:ALLOWABLE*/

  sfb=m.sft; /*"S STANDARD" 5.1.(4).(b)*/

  if(Nd<0.0) sNa=allowabletensionofs(m.sF,sect);
  if(Nd>0.0) sNa=allowablecompressionofs(m.sE,m.sF,lk,sect,axis);

  if(period==PSHORT) sNa*=1.5;

  if((Nd/sNa)>=1.0) return 0.0;

  Ze=coeffZ(sect.nsteel,sect.srect,axis);

  sMa=Ze*sfb*(1.0-Nd/sNa); /*"S STANDARD"(6.1)(6.2)(6.3)(6.4)*/

  if(sMa<0.0) sMa=0.0;

  /*
  fprintf(fout0,"WEAK. Nd=%.3f[kgf] sNa=%.3f[kgf] Nd/sNa=%.3f\n",
          Nd,sNa,(Nd/sNa));
  fprintf(fout0,"sfb=%.3f[kgf/cm2]",sfb);
  fprintf(fout0," Z=%.3f[cm3]\n",Ze);
  fprintf(fout0,"鉄骨軸力 sN=%.3f[tf] 許容値 sMa=%.3f[tfm]",
          Nd/1000.0,sMa/100000.0);
  */

  return sMa;
}/*allowablebendingofsweak*/

double allowultimshearofhstrong(int period,
                                double F,
                                struct section sect,int axis)
/*ALLOWABLE,ULTIMATE SHEAR OF H,[ STRONG.*/
{
  double H,B,tw,tf;
  double sfsa,sQa;                    /*s:STEEL s:SHEAR a:ALLOWABLE*/

  if(sect.nsteel!=3)
  {
    fprintf(fout0,"STEEL NOT H,[.\n");
    return 0.0;
  }

  if(axis==SX) /*SECTION AXIS AROUND X,WITH STRONG H.*/
  {
    H =(sect.srect[0].top-sect.srect[2].bottom);             /*[cm]*/
    B =(sect.srect[0].right-sect.srect[0].left);             /*[cm]*/
    tw=(sect.srect[1].right-sect.srect[1].left);             /*[cm]*/
    tf=(sect.srect[0].top-sect.srect[0].bottom);             /*[cm]*/
  }
  if(axis==SY) /*SECTION AXIS AROUND Y,WITH STRONG H.*/
  {
    H =(sect.srect[2].right-sect.srect[0].left);             /*[cm]*/
    B =(sect.srect[0].top-sect.srect[0].bottom);             /*[cm]*/
    tw=(sect.srect[1].top-sect.srect[1].bottom);             /*[cm]*/
    tf=(sect.srect[0].right-sect.srect[0].left);             /*[cm]*/
  }

  sfsa=F/sqrt(3.0)/1.5;    /*"S STANDARD"(5.1),"SRC STANDARD"(127).*/

  sQa=sfsa*tw*(H-2.0*tf);       /*[kgf] SAME AS "SRC STANDARD"(42).*/
  if(period==PSHORT || period==PULTIMATE) sQa*=1.5;

  /*
  fprintf(fout0,"STRONG. H=%.1f[cm] B=%.1f[cm]",H,B);
  fprintf(fout0," tw=%.1f[cm] tf=%.1f[cm]\n",tw,tf);
  fprintf(fout0,"sfs=%.3f[kgf/cm2]\n",sfsa);
  */

  return sQa;
}/*allowultimshearofhstrong*/

double allowultimshearofhweak(int period,
                              double F,
                              struct section sect,int axis)
/*ALLOWABLE,ULTIMATE SHEAR OF H,[ WEAK.*/
{
  double H,B,tw,tf;
  double sfsa,sQa;                    /*s:STEEL s:SHEAR a:ALLOWABLE*/

  if(sect.nsteel!=3)
  {
    fprintf(fout0,"STEEL NOT H,[.\n");
    return 0.0;
  }

  if(axis==SY) /*SECTION AXIS AROUND Y,WITH WEAK H.*/
  {
    H =(sect.srect[0].top-sect.srect[2].bottom);             /*[cm]*/
    B =(sect.srect[0].right-sect.srect[0].left);             /*[cm]*/
    tw=(sect.srect[1].right-sect.srect[1].left);             /*[cm]*/
    tf=(sect.srect[0].top-sect.srect[0].bottom);             /*[cm]*/
  }
  if(axis==SX) /*SECTION AXIS AROUND X,WITH WEAK H.*/
  {
    H =(sect.srect[2].right-sect.srect[0].left);             /*[cm]*/
    B =(sect.srect[0].top-sect.srect[0].bottom);             /*[cm]*/
    tw=(sect.srect[1].top-sect.srect[1].bottom);             /*[cm]*/
    tf=(sect.srect[0].right-sect.srect[0].left);             /*[cm]*/
  }
  /*fprintf(fout0,"H-%.0fx%.0fx%.0fx%.0f ",
          H*10.0,B*10.0,tw*10.0,tf*10.0);*/

  sfsa=F/sqrt(3.0)/1.5;    /*"S STANDARD"(5.1),"SRC STANDARD"(127).*/

  sQa=sfsa*2.0/3.0*(2.0*tf*B);  /*[kgf] SAME AS "SRC STANDARD"(61).*/
  if(period==PSHORT || period==PULTIMATE) sQa*=1.5;

  /*
  fprintf(fout0,"WEAK. H=%.1f[cm] B=%.1f[cm]",H,B);
  fprintf(fout0," tw=%.1f[cm] tf=%.1f[cm]\n",tw,tf);
  fprintf(fout0,"sfs=%.3f[kgf/cm2]\n",sfsa);
  */

  return sQa;
}/*allowultimshearofhweak*/

double coeffA(int nrect,struct materialrect rect[])
/*COEFFICIENT A OF RECTANGLES.*/
{
  int n;
  double w,h;                                                /*[cm]*/
  double A=0.0;                                             /*[cm2]*/

  for(n=1;n<=nrect;n++)
  {
    w=rect[n-1].right-rect[n-1].left;
    h=rect[n-1].top-rect[n-1].bottom;
    A+=fabs(w*h);                                           /*[cm2]*/
  }

  return A;
}/*coeffA*/

double coeffI(int nrect,struct materialrect rect[],int axis,
              double *yg,double *yt,double *yc)
/*COEFFICIENT I OF RECTANGLES.*/
{
  int n;
  double h[MAXSRECT],w[MAXSRECT];
  double yi[MAXSRECT],yj[MAXSRECT];
  double tA,dA;                                             /*[cm2]*/
  double I=0.0;                                             /*[cm4]*/

  tA=0.0;
  (*yg)=0.0;
  (*yt)=1000.0; (*yc)=-1000.0;

  for(n=0;n<nrect;n++)
  {
    if(axis==SX) /*AROUND x.*/
    {
      w[n]=rect[n].right-rect[n].left; /*WIDTH*/
      h[n]=rect[n].top-rect[n].bottom; /*HEIGHT*/
      yi[n]=rect[n].top;
      yj[n]=rect[n].bottom;
    }
    if(axis==SY) /*AROUND y.*/
    {
      w[n]=rect[n].top-rect[n].bottom; /*WIDTH*/
      h[n]=rect[n].right-rect[n].left; /*HEIGHT*/
      yi[n]=rect[n].right;
      yj[n]=rect[n].left;
    }

    dA=h[n]*w[n];
    tA+=dA; /*AREA OF ALL.*/
    (*yg)+=dA*0.5*(yi[n]+yj[n]);

    if(yi[n]>(*yc)) (*yc)=yi[n]; /*TOP OF ALL.*/
    if(yj[n]<(*yt)) (*yt)=yj[n]; /*BOTTOM OF ALL.*/
  }
  (*yg)/=tA; /*CENTER OF GRAVITY.*/

  for(n=0;n<nrect;n++)
  {
    I=I+1.0/3.0*w[n]*pow((yi[n]-(*yg)),3.0)
       -1.0/3.0*w[n]*pow((yj[n]-(*yg)),3.0);                /*[cm4]*/
  }

  /*fprintf(stderr,"I[cm4]=%.5f\n",I);*/

  return I;
}/*coeffI*/

double coeffZ(int nrect,struct materialrect rect[],int axis)
/*COEFFICIENT Z OF RECTANGLES.*/
{
  double I,yg,yt,yc;
  double Zupper,Zlower;
  double Ze;                                                /*[cm3]*/

  I=coeffI(nrect,rect,axis,&yg,&yt,&yc);

  Zupper=I/fabs(yc-yg);
  Zlower=I/fabs(yg-yt);

  if(Zupper<=Zlower) Ze=Zupper;
  else               Ze=Zlower;

  return Ze;
}/*coeffZ*/

double coeffi(int nrect,struct materialrect rect[],int axis)
/*COEFFICIENT i OF RECTANGLES.*/
{
  double I,A,yg,yt,yc;
  double i;                                                  /*[cm]*/

  A=coeffA(nrect,rect);
  I=coeffI(nrect,rect,axis,&yg,&yt,&yc);

  i=sqrt(I/A);

  return i;
}/*coeffi*/

double allowablebendingofpccolumn(struct element elem,
                                  struct materials m,
                                  int axis, /*0:Mx 1:My*/
                                  double Nd /*[kgf]*/)
/*RETURN:ALLOWABLE BENDING OF PC COLUMN*/
/*N=-:TENSION +:COMPRESSION [kgf]*/
/*FOR ONLY 1 CRECT.*/
{
  int i;
  double cBi[MAXCRECT],cDi[MAXCRECT];    /*WIDTH,DEPTH OF CONCRETES*/
  double pAi[MAXSTRND];
  double pYi[MAXSTRND];
  double cYi[MAXCRECT],cYj[MAXCRECT];
  double cYt,cYc,pYt,pYc;

  double cfc,cft,pft;                          /*ALLOWABLE STRESSES*/
  /*double Nover;*/                             /*REMAINING TENSION*/
  double Ma;                                          /*a:ALLOWABLE*/
  double Ninit;               /*TOTAL INITIAL PRESTRESS IN STRANDS.*/

  double Ac,Ap;
  double pAY;
  double cZ1,cZ2; /*SECTION FACTOR OF CONCRETE.*/
  double Ma1,Ma2;

  Ac=0.0; /*TOTAL AREA OF CONCRETE.[cm2]*/
  cYc=-1000.0; cYt=1000.0; /*[cm]*/
  for(i=0;i<(elem.sect->nconc);i++)
  {
    if(axis==0)
    {
      cBi[i]=fabs(elem.sect->crect[i].right
                 -elem.sect->crect[i].left);
      cYi[i]=elem.sect->crect[i].top;
      cYj[i]=elem.sect->crect[i].bottom;
    }
    if(axis==1)
    {
      cBi[i]=fabs(elem.sect->crect[i].top
                 -elem.sect->crect[i].bottom);
      cYi[i]=elem.sect->crect[i].right;
      cYj[i]=elem.sect->crect[i].left;
    }
    cDi[i]=cYi[i]-cYj[i];

    Ac+=(cBi[i]*cDi[i]); /*[cm2]*/

    if(cYc<cYi[i]) cYc=cYi[i];
    if(cYt>cYj[i]) cYt=cYj[i];
  }

  Ap=0.0; /*TOTAL AREA OF PC STRANDS.[cm2]*/
  pAY=0.0;
  pYc=-1000.0; pYt=1000.0; /*[cm]*/
  Ninit=0.0;
  for(i=0;i<(elem.sect->nstrnd);i++)
  {
    pAi[i]=elem.sect->strnd[i].area[END];
    if(axis==0) pYi[i]=elem.sect->strnd[i].y[END];
    if(axis==1) pYi[i]=elem.sect->strnd[i].x[END];

    Ninit+=elem.sect->strnd[i].Ni[END]*1000.0; /*[kgf]*/

    Ap+=pAi[i];
    pAY+=pAi[i]*pYi[i];

    if(pYc<pYi[i]) pYc=pYi[i];
    if(pYt>pYi[i]) pYt=pYi[i];

    pft=m.stft[elem.sect->strnd[i].type]; /*ALLOWABLE PC STRESS*/
    if(pft<(elem.sect->strnd[i].Ni[END])/pAi[i])
    {
      fprintf(fout0,"TOO MUCH STRESS.\n");
    }
  }

  if(elem.sect->stype==PC)
  {
    cfc=m.pcfc; /*+:COMPRESSION*/
    cft=m.pcft; /*-:TENSION*/
  }
  else return 0.0;

  /*
  fprintf(fout0,"Ac=%.3f Ap=%.3f[cm2]\n",Ac,Ap);
  fprintf(fout0,"cfc=%.3f cft=%.3f [kgf/cm2]\n",cfc,cft);
  fprintf(fout0,"Ninit=%.3f Nd=%.3f[tf]\n",Ninit/1000.0,Nd/1000.0);
  */

  cZ1=1/6.0*cBi[0]*cDi[0]*cDi[0]; /*[cm3]*/
  cZ2=cZ1;

  Ma1=cZ1*(cfc-(Nd+(m.stfact)*Ninit)/Ac);
  Ma2=cZ2*(-cft+(Nd+(m.stfact)*Ninit)/Ac);

  if(Ma1<=Ma2) Ma=Ma1;
  else         Ma=Ma2;

  /*fprintf(fout0,"pcMa=%.3f[kgfcm]\n",Ma);*/

  return Ma;
}/*allowablebendingofpccolumn*/

double allowablebendingofpcgirder(struct element elem,
                                  struct materials m,
                                  int axis /*0:Mx 1:My*/)
/*RETURN:ALLOWABLE BENDING ON END OF PC GIRDER FOR LONG.*/
/*N=-:TENSION +:COMPRESSION [kgf]*/
/*FOR ONLY 1 CRECT.*/
/*FOR ONLY CENTER OF PRESTRESS ABOVE CENTER OF CONCRETE.*/
{
  int i;
  double cBi[MAXCRECT],cDi[MAXCRECT];    /*WIDTH,DEPTH OF CONCRETES*/
  double pAi[MAXSTRND];
  double pYi[MAXSTRND];
  double cYi[MAXCRECT],cYj[MAXCRECT];
  double cYt,cYc,pYt,pYc;

  double cfc,cft,pft;                          /*ALLOWABLE STRESSES*/
  /*double Nover;*/                             /*REMAINING TENSION*/
  double Ma;                                          /*a:ALLOWABLE*/
  double Ninit;               /*TOTAL INITIAL PRESTRESS IN STRANDS.*/

  double Ac,Ap;
  double cZ1,cZ2; /*SECTION FACTOR OF CONCRETE.*/
  double Ma1,Ma2;

  double cYg,pYg; /*CENTER OF CONCRETE,PRESTRESS.*/

  Ac=0.0; /*TOTAL AREA OF CONCRETE.[cm2]*/
  cYc=-1000.0; cYt=1000.0; /*[cm]*/
  cYg=0.0; /*[cm]*/
  for(i=0;i<(elem.sect->nconc);i++)
  {
    if(axis==0)
    {
      cBi[i]=fabs(elem.sect->crect[i].right
                 -elem.sect->crect[i].left);
      cYi[i]=elem.sect->crect[i].top;
      cYj[i]=elem.sect->crect[i].bottom;
    }
    if(axis==1)
    {
      cBi[i]=fabs(elem.sect->crect[i].top
                 -elem.sect->crect[i].bottom);
      cYi[i]=elem.sect->crect[i].right;
      cYj[i]=elem.sect->crect[i].left;
    }
    cDi[i]=cYi[i]-cYj[i];

    Ac+=(cBi[i]*cDi[i]); /*[cm2]*/
    cYg+=(cBi[i]*cDi[i])*0.5*(cYi[i]+cYj[i]);

    if(cYc<cYi[i]) cYc=cYi[i];
    if(cYt>cYj[i]) cYt=cYj[i];
  }
  cYg/=Ac;

  Ap=0.0; /*TOTAL AREA OF PC STRANDS.[cm2]*/
  pYc=-1000.0; pYt=1000.0; /*[cm]*/
  pYg=0.0; /*[cm]*/
  Ninit=0.0;
  for(i=0;i<(elem.sect->nstrnd);i++)
  {
    pAi[i]=elem.sect->strnd[i].area[END];
    if(axis==0) pYi[i]=elem.sect->strnd[i].y[END];
    if(axis==1) pYi[i]=elem.sect->strnd[i].x[END];

    Ninit+=elem.sect->strnd[i].Ni[END]*1000.0; /*[kgf]*/

    Ap+=pAi[i];
    pYg+=pYi[i]*elem.sect->strnd[i].Ni[END]*1000.0;

    if(pYc<pYi[i]) pYc=pYi[i];
    if(pYt>pYi[i]) pYt=pYi[i];

    pft=m.stft[elem.sect->strnd[i].type]; /*ALLOWABLE PC STRESS*/
    if(pft<(elem.sect->strnd[i].Ni[END])/pAi[i])
    {
      fprintf(fout0,"TOO MUCH STRESS.\n");
    }
  }
  pYg/=Ninit;

  if(elem.sect->stype==PC)
  {
    cfc=m.pcfc; /*+:COMPRESSION*/
    cft=m.pcft; /*-:TENSION*/
  }
  else return 0.0;

  /*
  fprintf(fout0,"Ac=%.3f Ap=%.3f[cm2]\n",Ac,Ap);
  fprintf(fout0,"cfc=%.3f cft=%.3f [kgf/cm2]\n",cfc,cft);
  fprintf(fout0,"Ninit=%.3f[tf]\n",Ninit/1000.0);
  */

  cZ1=1/6.0*cBi[0]*cDi[0]*cDi[0]; /*[cm3]*/
  cZ2=cZ1;

  Ma1=cZ1*( cfc-(m.stfact)*Ninit/Ac
               +(m.stfact)*Ninit*(pYg-cYg)/cZ1);
  Ma2=cZ2*(-cft+(m.stfact)*Ninit/Ac
               +(m.stfact)*Ninit*(pYg-cYg)/cZ2);

  if(Ma1<=Ma2) Ma=Ma1;
  else         Ma=Ma2;

  /*fprintf(fout0,"pcMa=%.3f[kgfcm]\n",Ma);*/

  return Ma;
}/*allowablebendingofpcgirder*/

double allowablebendingofpcgirderonmid(struct element elem,
                                       struct materials m,
                                       double length, /*[cm]*/
                                       int axis /*0:Mx 1:My*/)
/*RETURN:ALLOWABLE BENDING ON MID OF PC GIRDER FOR LONG.*/
/*N=-:TENSION +:COMPRESSION [kgf]*/
/*FOR ONLY 1 CRECT.*/
/*FOR ONLY CENTER OF PRESTRESS UNDER CENTER OF CONCRETE.*/
{
  int i;
  double cBi[MAXCRECT],cDi[MAXCRECT];    /*WIDTH,DEPTH OF CONCRETES*/
  double pAi[MAXSTRND];
  double pYi[MAXSTRND];
  double cYi[MAXCRECT],cYj[MAXCRECT];
  double cYt,cYc,pYt,pYc;

  double cfc,cft,cfco,cfto,pft;                /*ALLOWABLE STRESSES*/
  /*double Nover;*/                             /*REMAINING TENSION*/
  double Ma;                                          /*a:ALLOWABLE*/
  double Ninit;               /*TOTAL INITIAL PRESTRESS IN STRANDS.*/

  double Ac,Ap;
  double cZ1,cZ2; /*SECTION FACTOR OF CONCRETE.*/
  double Ma1,Ma2;
  double Wdead,Mdead; /*DEAD LOAD OF GIRDER.*/
  double cYg,pYg; /*CENTER OF CONCRETE,PRESTRESS.*/

  Ac=0.0; /*TOTAL AREA OF CONCRETE.[cm2]*/
  cYc=-1000.0; cYt=1000.0; /*[cm]*/
  cYg=0.0; /*[cm]*/
  for(i=0;i<(elem.sect->nconc);i++)
  {
    if(axis==0)
    {
      cBi[i]=fabs(elem.sect->crect[i].right
                 -elem.sect->crect[i].left);
      cYi[i]=elem.sect->crect[i].top;
      cYj[i]=elem.sect->crect[i].bottom;
    }
    if(axis==1)
    {
      cBi[i]=fabs(elem.sect->crect[i].top
                 -elem.sect->crect[i].bottom);
      cYi[i]=elem.sect->crect[i].right;
      cYj[i]=elem.sect->crect[i].left;
    }
    cDi[i]=cYi[i]-cYj[i];

    Ac+=(cBi[i]*cDi[i]); /*[cm2]*/
    cYg+=(cBi[i]*cDi[i])*0.5*(cYi[i]+cYj[i]);

    if(cYc<cYi[i]) cYc=cYi[i];
    if(cYt>cYj[i]) cYt=cYj[i];
  }
  cYg/=Ac;

  Wdead=Ac*2.4/1000.0; /*[kgf/cm]*/
  Mdead=Wdead*length*length/8.0; /*Mo[kgfcm] BY DEAD LOAD.*/

  Ap=0.0; /*TOTAL AREA OF PC STRANDS.[cm2]*/
  pYc=-1000.0; pYt=1000.0; /*[cm]*/
  pYg=0.0; /*[cm]*/
  Ninit=0.0;
  for(i=0;i<(elem.sect->nstrnd);i++)
  {
    pAi[i]=elem.sect->strnd[i].area[MID];
    if(axis==0) pYi[i]=elem.sect->strnd[i].y[MID];
    if(axis==1) pYi[i]=elem.sect->strnd[i].x[MID];

    Ninit+=elem.sect->strnd[i].Ni[MID]*1000.0; /*[kgf]*/

    Ap+=pAi[i];
    pYg+=pYi[i]*elem.sect->strnd[i].Ni[MID]*1000.0;

    if(pYc<pYi[i]) pYc=pYi[i];
    if(pYt>pYi[i]) pYt=pYi[i];

    pft=m.stft[elem.sect->strnd[i].type]; /*ALLOWABLE PC STRESS*/
    if(pft<(elem.sect->strnd[i].Ni[MID])/pAi[i])
    {
      fprintf(fout0,"TOO MUCH STRESS.\n");
    }
  }
  pYg/=Ninit;

  if(elem.sect->stype==PC)
  {
    cfc =m.pcfc;  /*+:COMPRESSION*/
    cft =m.pcft;  /*-:TENSION*/
    cfco=m.pcfco; /*+:COMPRESSION*/
    cfto=m.pcfto; /*-:TENSION*/
  }
  else return 0.0;
  
  /*
  fprintf(fout0,"cfc =%.3f cft =%.3f [kgf/cm2]\n",cfc,cft);
  fprintf(fout0,"cfco=%.3f cfto=%.3f [kgf/cm2]\n",cfco,cfto);

  fprintf(fout0,"Ac=%.3f Ap=%.3f[cm2]\n",Ac,Ap);
  fprintf(fout0,"cYg=%.3f pYg=%.3f e=%.3f[cm]\n",cYg,pYg,(cYg-pYg));
  fprintf(fout0,"Ninit=%.3f[tf]\n",Ninit/1000.0);
  fprintf(fout0,"Mdead=%.3f[tfm]\n",Mdead/100000.0);
  */

  cZ1=1/6.0*cBi[0]*cDi[0]*cDi[0]; /*[cm3]*/
  cZ2=cZ1;
  /*fprintf(fout0,"Z1=%.0f Z2=%.0f[cm3]\n",cZ1,cZ2);*/

  Ma1=cZ1*( cfto-Ninit/Ac+Ninit*(cYg-pYg)/cZ1);
  Ma2=cZ2*(-cfco+Ninit/Ac+Ninit*(cYg-pYg)/cZ2);

  if(Mdead<Ma1 || Mdead<Ma2)
  {
    fprintf(fout0,"\nTOO MUCH STRESS FOR DEAD LOAD.\n");
    return 0.0;
  }

  Ma1=cZ1*( cfc-(m.stfact)*Ninit/Ac
               +(m.stfact)*Ninit*(cYg-pYg)/cZ1);
  Ma2=cZ2*(-cft+(m.stfact)*Ninit/Ac
               +(m.stfact)*Ninit*(cYg-pYg)/cZ2);

  if(Ma1<=Ma2) Ma=Ma1;
  else         Ma=Ma2;

  /*fprintf(fout0,"pcMa=%.3f[kgfcm]\n",Ma);*/

  return Ma;
}/*allowablebendingofpcgirderonmid*/

double ultimatebendingofpccolumn(struct element elem,
                                 struct materials m,
                                 int axis, /*0:Mx 1:My*/
                                 double Nd /*[kgf]*/)
/*RETURN:ULTIMATE BENDING OF PC COLUMN*/
/*N=-:TENSION +:COMPRESSION [kgf]*/
/*FOR ONLY 1 CRECT.*/
{
  int i;
  double cBi[MAXCRECT],cDi[MAXCRECT];    /*WIDTH,DEPTH OF CONCRETES*/
  double rAi[MAXREINS],pAi[MAXSTRND];
  double rYi[MAXREINS],pYi[MAXSTRND];
  double cYi[MAXCRECT],cYj[MAXCRECT];
  double cYt,cYc,rYt,rYc,pYt,pYc;

  double cfc,/*cft,*/pft;                      /*ALLOWABLE STRESSES*/
  /*double Nover;*/                             /*REMAINING TENSION*/
  double Mu;                                           /*u:ULTIMATE*/
  double Ninit;               /*TOTAL INITIAL PRESTRESS IN STRANDS.*/
  double rpNy;           /*TOTAL YEILD STRESS OF REINS AND STRANDS.*/

  double Ac,Ar,Ap;
  double rAY,pAY;
  double Muc,Mup;

  double k1,k2,dn,pN;

  Ac=0.0; /*TOTAL AREA OF CONCRETE.[cm2]*/
  cYc=-1000.0; cYt=1000.0; /*[cm]*/
  for(i=0;i<(elem.sect->nconc);i++)
  {
    if(axis==0)
    {
      cBi[i]=fabs(elem.sect->crect[i].right
                 -elem.sect->crect[i].left);
      cYi[i]=elem.sect->crect[i].top;
      cYj[i]=elem.sect->crect[i].bottom;
    }
    if(axis==1)
    {
      cBi[i]=fabs(elem.sect->crect[i].top
                 -elem.sect->crect[i].bottom);
      cYi[i]=elem.sect->crect[i].right;
      cYj[i]=elem.sect->crect[i].left;
    }
    cDi[i]=cYi[i]-cYj[i];

    Ac+=(cBi[i]*cDi[i]); /*[cm2]*/

    if(cYc<cYi[i]) cYc=cYi[i];
    if(cYt>cYj[i]) cYt=cYj[i];
  }

  Ar=0.0; /*TOTAL AREA OF REINS.[cm2]*/
  rAY=0.0;
  rYc=-1000.0; rYt=1000.0; /*[cm]*/
  rpNy=0.0;
  for(i=0;i<(elem.sect->nrein);i++)
  {
    rAi[i]=elem.sect->rein[i].area;
    if(axis==0) rYi[i]=elem.sect->rein[i].y;
    if(axis==1) rYi[i]=elem.sect->rein[i].x;

    Ar+=rAi[i];
    rAY+=rAi[i]*rYi[i];
    rpNy+=m.rftu*rAi[i]; /*[kgf]*/

    if(rYc<rYi[i]) rYc=rYi[i];
    if(rYt>rYi[i]) rYt=rYi[i];
  }

  Ap=0.0; /*TOTAL AREA OF PC STRANDS.[cm2]*/
  pAY=0.0;
  pYc=-1000.0; pYt=1000.0; /*[cm]*/
  Ninit=0.0;
  for(i=0;i<(elem.sect->nstrnd);i++)
  {
    pAi[i]=elem.sect->strnd[i].area[END];
    if(axis==0) pYi[i]=elem.sect->strnd[i].y[END];
    if(axis==1) pYi[i]=elem.sect->strnd[i].x[END];

    Ninit+=elem.sect->strnd[i].Ni[END]*1000.0; /*[kgf]*/
    rpNy+=m.stftu[elem.sect->strnd[i].type]*pAi[i]; /*[kgf]*/

    Ap+=pAi[i];
    pAY+=pAi[i]*pYi[i];

    if(pYc<pYi[i]) pYc=pYi[i];
    if(pYt>pYi[i]) pYt=pYi[i];

    pft=m.stft[elem.sect->strnd[i].type]; /*ALLOWABLE PC STRESS*/
    if(pft<(elem.sect->strnd[i].Ni[END])/pAi[i])
    {
      fprintf(fout0,"TOO MUCH STRESS.\n");
    }
  }

  if(elem.sect->stype==PC) cfc=m.pcF; /*+:COMPRESSION*/
  else return 0.0;

  /*
  fprintf(fout0,"Ac=%.3f Ap=%.3f[cm2]\n",Ac,Ap);
  fprintf(fout0,"cfc=%.3f cft=%.3f [kgf/cm2]\n",cfc,cft);
  fprintf(fout0,"Ninit=%.3f Nd=%.3f[tf]\n",Ninit/1000.0,Nd/1000.0);
  */

  k1=1.0; /*EFFECTIVE Fc*/
  k2=0.5; /*CENTER OF STRESS BY COMPRESSED CONCRETE.*/

  /*COMPRESSED CONCRETE DEPTH = NEUTRAL AXIS.*/
  dn=(Nd+(rpNy/2))/(k1*cfc)/cBi[0];
  if(dn>(cDi[0]/2))
  {
    fprintf(fout0,"COMPRESSED CONCRETE OVERLAPPED.\n");
  }

  Muc=(Nd+(rpNy/2))*((1-k2)*dn);
  Mup=0.0;
  for(i=0;i<(elem.sect->nstrnd);i++)
  {
    if(pYi[i]<=0)
    {
      pft=m.stftu[elem.sect->strnd[i].type]; /*+:TENSION*/
      pN=pft*pAi[i];
      Mup+=pN*(cYc-dn-pYi[i]);

      /*fprintf(fout0,"pN=%.3f[tf] dc=%.3f[cm] Mup1=%.3f[tfm]\n",
              pN/1000.0,(cYc-dn-pYi[i]),
              pN*(cYc-dn-pYi[i])/100000.0);*/
    }
  }
  for(i=0;i<(elem.sect->nrein);i++)
  {
    if(rYi[i]<=0)
    {
      pN=m.rftu*rAi[i];
      Mup+=pN*(cYc-dn-rYi[i]);
    }
  }

  Mu=Muc+Mup;
  /*
  fprintf(fout0,"Nd=%.3f rpNy=%.3f[tf]\n",Nd/1000.0,rpNy/1000.0);
  fprintf(fout0,"dn=%.3f[cm]\n",dn);
  fprintf(fout0,"Muc=%.3f Mup=%.3f[kgfcm]\n",Muc,Mup);
  */
  /*fprintf(fout0,"pcMu=%.3f[tfm]\n",Mu/100000.0);*/

  return Mu;
}/*ultimatebendingofpccolumn*/

double ultimatebendingofpcgirder(struct element elem,
                                 struct materials m,
                                 int axis, /*0:Mx 1:My*/
                                 int tensedside)
/*RETURN:ULTIMATE BENDING OF PC GIRDER*/
/*N=-:TENSION +:COMPRESSION [kgf]*/
/*FOR ONLY 1 CRECT.*/
{
  int i;
  double cBi[MAXCRECT],cDi[MAXCRECT];    /*WIDTH,DEPTH OF CONCRETES*/
  double pAi[MAXSTRND];
  double pYi[MAXSTRND];
  double cYi[MAXCRECT],cYj[MAXCRECT];
  double cYt,cYc,pYt,pYc;

  double cfc,/*cft,*/pft;                      /*ALLOWABLE STRESSES*/
  /*double Nover;*/                             /*REMAINING TENSION*/
  double Mu;                                           /*u:ULTIMATE*/
  double Ninit;               /*TOTAL INITIAL PRESTRESS IN STRANDS.*/
  double rpNy,rpNyi;     /*TOTAL YEILD STRESS OF REINS AND STRANDS.*/

  double Ac,Ap;
  double cYg,pYg;
  double Muc,Mup;

  double k1,k2,dn,pN;

  Ac=0.0; /*TOTAL AREA OF CONCRETE.[cm2]*/
  cYc=-1000.0; cYt=1000.0; /*[cm]*/
  cYg=0.0;
  for(i=0;i<(elem.sect->nconc);i++)
  {
    if(axis==0)
    {
      cBi[i]=fabs(elem.sect->crect[i].right
                 -elem.sect->crect[i].left);
      cYi[i]=elem.sect->crect[i].top;
      cYj[i]=elem.sect->crect[i].bottom;
    }
    if(axis==1)
    {
      cBi[i]=fabs(elem.sect->crect[i].top
                 -elem.sect->crect[i].bottom);
      cYi[i]=elem.sect->crect[i].right;
      cYj[i]=elem.sect->crect[i].left;
    }
    cDi[i]=cYi[i]-cYj[i];

    Ac+=(cBi[i]*cDi[i]); /*[cm2]*/
    cYg+=(cBi[i]*cDi[i])*0.5*(cYi[i]+cYj[i]);

    if(cYc<cYi[i]) cYc=cYi[i];
    if(cYt>cYj[i]) cYt=cYj[i];
  }
  cYg/=Ac;

  Ap=0.0; /*TOTAL AREA OF PC STRANDS.[cm2]*/
  pYg=0.0;
  pYc=-1000.0; pYt=1000.0; /*[cm]*/
  Ninit=0.0;
  rpNy=0.0;
  for(i=0;i<(elem.sect->nstrnd);i++)
  {
    pAi[i]=elem.sect->strnd[i].area[END];
    if(axis==0) pYi[i]=elem.sect->strnd[i].y[END];
    if(axis==1) pYi[i]=elem.sect->strnd[i].x[END];

    Ap+=pAi[i];
    Ninit+=elem.sect->strnd[i].Ni[END]*1000.0; /*[kgf]*/

    if((tensedside==UPPER && pYi[i]>cYg) ||
       (tensedside==LOWER && pYi[i]<cYg))
    {
      rpNyi=m.stftu[elem.sect->strnd[i].type]*pAi[i]; /*[kgf]*/
      rpNy+=rpNyi;

      pYg+=rpNyi*pYi[i];
    }

    if(pYc<pYi[i]) pYc=pYi[i];
    if(pYt>pYi[i]) pYt=pYi[i];

    pft=m.stft[elem.sect->strnd[i].type]; /*ALLOWABLE PC STRESS*/
    if(pft<(elem.sect->strnd[i].Ni[END])/pAi[i])
    {
      fprintf(fout0,"TOO MUCH STRESS.\n");
    }
  }
  pYg/=rpNy; /*CENTER OF YIELD STRESS OF STRANDS.*/

  if(elem.sect->stype==PC) cfc=m.pcF; /*+:COMPRESSION*/
  else return 0.0;

  k1=1.0; /*EFFECTIVE Fc*/
  k2=0.5; /*CENTER OF STRESS BY COMPRESSED CONCRETE.*/

  /*COMPRESSED CONCRETE DEPTH = NEUTRAL AXIS.*/
  dn=rpNy/(k1*cfc)/cBi[0];
  if((tensedside==UPPER && dn>(cYg-cYt)) ||
     (tensedside==LOWER && dn>(cYc-cYg)))
  {
    fprintf(fout0,"COMPRESSED CONCRETE OVERLAPPED.\n");
  }

  Muc=rpNy*((1-k2)*dn);
  Mup=0.0;
  for(i=0;i<(elem.sect->nstrnd);i++)
  {
    pft=m.stftu[elem.sect->strnd[i].type]; /*+:TENSION*/
    pN=pft*pAi[i];

    if(tensedside==LOWER && pYi[i]<cYg)
    {
      Mup+=pN*(cYc-dn-pYi[i]);
      /*fprintf(fout0,"pN=%.3f[tf] dc=%.3f[cm] Mup1=%.3f[tfm]\n",
              pN/1000.0,(cYc-pYi[i]),pN*(cYc-dn-pYi[i])/100000.0);*/
    }
    else if(tensedside==UPPER && pYi[i]>cYg)
    {
      Mup+=pN*(pYi[i]-cYt-dn);
      /*fprintf(fout0,"pN=%.3f[tf] dc=%.3f[cm] Mup1=%.3f[tfm]\n",
              pN/1000.0,(pYi[i]-cYt),pN*(pYi[i]-cYt-dn)/100000.0);*/
    }
  }

  Mu=Muc+Mup;
  /*
  fprintf(fout0,"Ac=%.3f Ap=%.3f[cm2]\n",Ac,Ap);
  fprintf(fout0,"cfc=%.3f cft=%.3f [kgf/cm2]\n",cfc,cft);
  fprintf(fout0,"Ninit=%.3f rpNy=%.3f[tf]\n",
          Ninit/1000.0,rpNy/1000.0);

  fprintf(fout0,"dn=%.3f[cm]\n",dn);
  fprintf(fout0,"Muc=%.3f Mup=%.3f[kgfcm]\n",Muc,Mup);
  fprintf(fout0,"pcMu=%.5f[tfm]\n",Mu/100000.0);
  */
  return Mu;
}/*ultimatebendingofpcgirder*/

double allowableshearofpccolumn(struct element elem,
                                struct materials m,
                                int axis, /*0:Qy 1:Qx*/
                                double Nd /*[kgf]*/)
/*RETURN:ALLOWABLE SHEAR OF PC COLUMN,GIRDER.*/
/*N=-:TENSION +:COMPRESSION [kgf]*/
/*FOR ONLY 1 CRECT.*/
/*FOR FULL PRESTRESSING.*/
{
  int i;
  double cBi[MAXCRECT],cDi[MAXCRECT];    /*WIDTH,DEPTH OF CONCRETES*/
  double pAi[MAXSTRND];
  double pYi[MAXSTRND];
  double cYi[MAXCRECT],cYj[MAXCRECT];
  double cYt,cYc,pYt,pYc,cYg;

  double /*cfc,*/cft,pft;                      /*ALLOWABLE STRESSES*/
  double sg;
  /*double Nover;*/                             /*REMAINING TENSION*/
  double Qa;                                          /*a:ALLOWABLE*/
  double Ninit;               /*TOTAL INITIAL PRESTRESS IN STRANDS.*/

  double Ac,Ap;
  double pAY;
  double cSn,cI; /*SECTION FACTOR OF CONCRETE.*/

  Ac=0.0; /*TOTAL AREA OF CONCRETE.[cm2]*/
  cYg=0.0; /*CENTER OF GRAVITY OF CONCRETE.[cm]*/
  cYc=-1000.0; cYt=1000.0; /*[cm]*/
  for(i=0;i<(elem.sect->nconc);i++)
  {
    if(axis==SX) /*FOR Qy*/
    {
      cBi[i]=fabs(elem.sect->crect[i].right
                 -elem.sect->crect[i].left);
      cYi[i]=elem.sect->crect[i].top;
      cYj[i]=elem.sect->crect[i].bottom;
    }
    if(axis==SY) /*FOR Qx*/
    {
      cBi[i]=fabs(elem.sect->crect[i].top
                 -elem.sect->crect[i].bottom);
      cYi[i]=elem.sect->crect[i].right;
      cYj[i]=elem.sect->crect[i].left;
    }
    cDi[i]=cYi[i]-cYj[i];

    Ac+=(cBi[i]*cDi[i]); /*[cm2]*/
    cYg+=(cBi[i]*cDi[i])*0.5*(cYi[i]+cYj[i]);

    if(cYc<cYi[i]) cYc=cYi[i];
    if(cYt>cYj[i]) cYt=cYj[i];
  }
  cYg/=Ac;

  Ap=0.0; /*TOTAL AREA OF PC STRANDS.[cm2]*/
  pAY=0.0;
  pYc=-1000.0; pYt=1000.0; /*[cm]*/
  Ninit=0.0;
  for(i=0;i<(elem.sect->nstrnd);i++)
  {
    pAi[i]=elem.sect->strnd[i].area[END];
    if(axis==SX) pYi[i]=elem.sect->strnd[i].y[END];
    if(axis==SY) pYi[i]=elem.sect->strnd[i].x[END];

    Ninit+=elem.sect->strnd[i].Ni[END]*1000.0; /*[kgf]*/

    Ap+=pAi[i];
    pAY+=pAi[i]*pYi[i];

    if(pYc<pYi[i]) pYc=pYi[i];
    if(pYt>pYi[i]) pYt=pYi[i];

    pft=m.stft[elem.sect->strnd[i].type]; /*ALLOWABLE PC STRESS*/
    if(pft<(elem.sect->strnd[i].Ni[END])/pAi[i])
    {
      fprintf(fout0,"TOO MUCH STRESS.\n");
    }
  }

  if(elem.sect->stype==PC)
  {
    cft=-0.07*(m.pcfc); /*-:TENSION*/
  }
  else return 0.0;

  /*
  fprintf(fout0,"Ac=%.3f Ap=%.3f[cm2]\n",Ac,Ap);
  fprintf(fout0,"cfc=%.3f cft=%.3f [kgf/cm2]\n",cfc,cft);
  fprintf(fout0,"Ninit=%.3f Nd=%.3f[tf]\n",Ninit/1000.0,Nd/1000.0);
  */

  cI=1/12.0*cBi[0]*cDi[0]*cDi[0]*cDi[0]; /*[cm4]*/

  cSn=0.0; /*NEUTRAL AXIS ON CENTER OF CONCRETE.[cm3]*/
  for(i=0;i<(elem.sect->nconc);i++)
  {
    if(cYj[i]>=cYg)
    {
      cSn+=cBi[i]*cDi[i]*(0.5*(cYi[i]+cYj[i])-cYg);
    }
    else if(cYi[i]>=cYg)
    {
      cSn+=cBi[i]*(cYi[i]-cYg)*(0.5*(cYi[i]-cYg));
    }
  }

  sg=(Nd+Ninit)/Ac;
  if(sg<=cft) return 0.0;

  Qa=cBi[0]*cI/cSn*sqrt(cft*(cft-sg));
  /*
  fprintf(fout0,"cSn=%.3f[cm3] cI=%.3f[cm4]\n",cSn,cI);
  fprintf(fout0,"cft=%.3f sg=%.3f[kgf/cm2]\n",cft,sg);
  */
  /*fprintf(fout0,"pcQa=%.3f[kgf]\n",Qa);*/

  return Qa;
}/*allowableshearofpccolumn*/

double ultimateshearofpccolumn(struct element elem,
                               struct materials m,
                               int axis, /*0:Qy 1:Qx*/
                               double Nd,double Qd,double Md)
/*RETURN:ULTIMATE SHEAR OF PC COLUMN.*/
/*N=-:TENSION +:COMPRESSION [kgf]*/
/*FOR ONLY 1 CRECT.*/
/*FOR FULL PRESTRESSING.*/
{
  int i;
  double cBi[MAXCRECT],cDi[MAXCRECT];    /*WIDTH,DEPTH OF CONCRETES*/
  double pAi[MAXSTRND];
  double pYi[MAXSTRND];
  double cYi[MAXCRECT],cYj[MAXCRECT];
  double cYt,cYc,pYt,pYgt,pYc,cYg;
  double pw;

  double cfs,pft,wft;                          /*ALLOWABLE STRESSES*/
  double sg;
  double Qu;                                /*a:ALLOWABLE*/
  double Ninit;               /*TOTAL INITIAL PRESTRESS IN STRANDS.*/

  double Ac,Ap,pAt,cJ;
  double pAY;
  double alpha;

  Ac=0.0; /*TOTAL AREA OF CONCRETE.[cm2]*/
  cYg=0.0; /*CENTER OF GRAVITY OF CONCRETE.[cm]*/
  cYc=-1000.0; cYt=1000.0; /*[cm]*/
  for(i=0;i<(elem.sect->nconc);i++)
  {
    if(axis==SX) /*FOR Qy*/
    {
      cBi[i]=fabs(elem.sect->crect[i].right
                 -elem.sect->crect[i].left);
      cYi[i]=elem.sect->crect[i].top;
      cYj[i]=elem.sect->crect[i].bottom;
    }
    if(axis==SY) /*FOR Qx*/
    {
      cBi[i]=fabs(elem.sect->crect[i].top
                 -elem.sect->crect[i].bottom);
      cYi[i]=elem.sect->crect[i].right;
      cYj[i]=elem.sect->crect[i].left;
    }
    cDi[i]=cYi[i]-cYj[i];

    Ac+=(cBi[i]*cDi[i]); /*[cm2]*/
    cYg+=(cBi[i]*cDi[i])*0.5*(cYi[i]+cYj[i]);

    if(cYc<cYi[i]) cYc=cYi[i];
    if(cYt>cYj[i]) cYt=cYj[i];
  }
  cYg/=Ac;

  Ap=0.0; /*TOTAL AREA OF PC STRANDS.[cm2]*/
  pAt=0.0; /*AREA OF TENSION STRANDS.*/
  pAY=0.0;
  pYgt=0.0; /*CENTER OF TENSION STRANDS.*/
  pYc=-1000.0; pYt=1000.0; /*[cm]*/
  Ninit=0.0;
  for(i=0;i<(elem.sect->nstrnd);i++)
  {
    pAi[i]=elem.sect->strnd[i].area[END];
    if(axis==SX) pYi[i]=elem.sect->strnd[i].y[END];
    if(axis==SY) pYi[i]=elem.sect->strnd[i].x[END];

    Ninit+=elem.sect->strnd[i].Ni[END]*1000.0; /*[kgf]*/

    Ap += pAi[i];
    pAY += pAi[i]*pYi[i];
    if(pYi[i]<cYg)
    {
      pAt += pAi[i];
      pYgt += pAi[i]*pYi[i];
    }
    if(pYc<pYi[i]) pYc=pYi[i];
    if(pYt>pYi[i]) pYt=pYi[i];

    pft=m.stft[elem.sect->strnd[i].type]; /*ALLOWABLE PC STRESS*/
    if(pft<(elem.sect->strnd[i].Ni[END])/pAi[i])
    {
      fprintf(fout0,"TOO MUCH STRESS.\n");
    }
  }
  pYgt/=pAt;

  if(elem.sect->stype==PC) cfs=m.pcfsu;
  else return 0.0;

  wft=m.wfp;
  pw=elem.sect->shearrein[1-axis];
  if(pw>0.012) pw=0.012;

  /*
  fprintf(fout0,"Ac=%.3f Ap=%.3f[cm2]\n",Ac,Ap);
  fprintf(fout0,"cfc=%.3f cft=%.3f [kgf/cm2]\n",cfc,cft);
  fprintf(fout0,"Ninit=%.3f Nd=%.3f[tf]\n",Ninit/1000.0,Nd/1000.0);
  */

  alpha=4.0/(Md/Qd/(cYc-pYgt)+1.0);
  if(alpha<1.0)      alpha=1.0;
  else if(2.0<alpha) alpha=2.0;

  sg=(Nd+Ninit)/Ac;

  cJ=7.0/8.0*(cYc-pYgt);

  Qu=(alpha*(cfs+0.1*sg)+0.5*wft*(pw-0.002))*cBi[0]*cJ;
  /*
  fprintf(fout0,"alpha=%.3f\n",alpha);
  fprintf(fout0,"cfs=%.3f sg=%.3f[kgf/cm2]\n",cfs,sg);
  fprintf(fout0,"wft=%.3f[kgf/cm2] pw=%.3f\n",wft,pw);
  fprintf(fout0,"cB=%.3f dc=%.3f cJ=%.3f[cm]\n",
          cBi[0],(cYc-pYgt),cJ);
  */
  /*fprintf(fout0,"pcQu=%.3f[kgf]\n",Qu);*/

  return Qu;
}/*ultimateshearofpccolumn*/

double ultimateshearofpcgirder(struct element elem,
                               struct materials m,
                               int axis, /*0:Qy 1:Qx*/
                               int tensedside,
                               double Nd,double Qd,double Md)
/*RETURN:ULTIMATE SHEAR OF PC GIRDER.*/
/*N=-:TENSION +:COMPRESSION [kgf]*/
/*FOR ONLY 1 CRECT.*/
/*FOR FULL PRESTRESSING.*/
{
  int i;
  double cBi[MAXCRECT],cDi[MAXCRECT];    /*WIDTH,DEPTH OF CONCRETES*/
  double pAi[MAXSTRND];
  double pYi[MAXSTRND];
  double cYi[MAXCRECT],cYj[MAXCRECT];
  double cYt,cYc,pYt,pYgt,pYc,cYg;
  double pw;

  double cfs,pft,wft;                          /*ALLOWABLE STRESSES*/
  double sg;
  double Qu;                                /*a:ALLOWABLE*/
  double Ninit;               /*TOTAL INITIAL PRESTRESS IN STRANDS.*/

  double Ac,Ap,pAt,cJ,de;
  double pAY;
  double alpha;

  Ac=0.0; /*TOTAL AREA OF CONCRETE.[cm2]*/
  cYg=0.0; /*CENTER OF GRAVITY OF CONCRETE.[cm]*/
  cYc=-1000.0; cYt=1000.0; /*[cm]*/
  for(i=0;i<(elem.sect->nconc);i++)
  {
    if(axis==SX) /*FOR Qy*/
    {
      cBi[i]=fabs(elem.sect->crect[i].right
                 -elem.sect->crect[i].left);
      cYi[i]=elem.sect->crect[i].top;
      cYj[i]=elem.sect->crect[i].bottom;
    }
    if(axis==SY) /*FOR Qx*/
    {
      cBi[i]=fabs(elem.sect->crect[i].top
                 -elem.sect->crect[i].bottom);
      cYi[i]=elem.sect->crect[i].right;
      cYj[i]=elem.sect->crect[i].left;
    }
    cDi[i]=cYi[i]-cYj[i];

    Ac+=(cBi[i]*cDi[i]); /*[cm2]*/
    cYg+=(cBi[i]*cDi[i])*0.5*(cYi[i]+cYj[i]);

    if(cYc<cYi[i]) cYc=cYi[i];
    if(cYt>cYj[i]) cYt=cYj[i];
  }
  cYg/=Ac;

  Ap=0.0; /*TOTAL AREA OF PC STRANDS.[cm2]*/
  pAt=0.0; /*AREA OF TENSION STRANDS.*/
  pAY=0.0;
  pYgt=0.0; /*CENTER OF TENSION STRANDS.*/
  pYc=-1000.0; pYt=1000.0; /*[cm]*/
  Ninit=0.0;
  for(i=0;i<(elem.sect->nstrnd);i++)
  {
    pAi[i]=elem.sect->strnd[i].area[END];
    if(axis==SX) pYi[i]=elem.sect->strnd[i].y[END];
    if(axis==SY) pYi[i]=elem.sect->strnd[i].x[END];

    Ninit+=elem.sect->strnd[i].Ni[END]*1000.0; /*[kgf]*/

    Ap += pAi[i];
    pAY += pAi[i]*pYi[i];

    if((tensedside==UPPER && pYi[i]>cYg) ||
       (tensedside==LOWER && pYi[i]<cYg))
    {
      pAt += pAi[i];
      pYgt += pAi[i]*pYi[i];
    }
    if(pYc<pYi[i]) pYc=pYi[i];
    if(pYt>pYi[i]) pYt=pYi[i];

    pft=m.stft[elem.sect->strnd[i].type]; /*ALLOWABLE PC STRESS*/
    if(pft<(elem.sect->strnd[i].Ni[END])/pAi[i])
    {
      fprintf(fout0,"TOO MUCH STRESS.\n");
    }
  }
  pYgt/=pAt;

  if(elem.sect->stype==PC) cfs=m.pcfsu;
  else return 0.0;

  wft=m.wfp;
  pw=elem.sect->shearrein[1-axis];
  if(pw>0.012) pw=0.012;

  /*
  fprintf(fout0,"Ac=%.3f Ap=%.3f[cm2]\n",Ac,Ap);
  fprintf(fout0,"cfc=%.3f cft=%.3f [kgf/cm2]\n",cfc,cft);
  fprintf(fout0,"Ninit=%.3f Nd=%.3f[tf]\n",Ninit/1000.0,Nd/1000.0);
  */

  if(tensedside==UPPER) de=pYgt-cYt;
  if(tensedside==LOWER) de=cYc-pYgt;

  alpha=4.0/(Md/Qd/de+1.0);
  if(alpha<1.0)      alpha=1.0;
  else if(2.0<alpha) alpha=2.0;

  sg=(Nd+Ninit)/Ac;

  cJ=7.0/8.0*de;

  Qu=(alpha*(cfs+0.1*sg)+0.5*wft*(pw-0.002))*cBi[0]*cJ;
  /*
  fprintf(fout0,"alpha=%.3f\n",alpha);
  fprintf(fout0,"cfs=%.3f sg=%.3f[kgf/cm2]\n",cfs,sg);
  fprintf(fout0,"wft=%.3f[kgf/cm2] pw=%.3f\n",wft,pw);
  fprintf(fout0,"cB=%.3f dc=%.3f cJ=%.3f[cm]\n",
          cBi[0],(cYc-pYgt),cJ);
  */
  /*fprintf(fout0,"pcQu=%.3f[kgf]\n",Qu);*/

  return Qu;
}/*ultimateshearofpcgirder*/

/*********/

