/*ARCLM101.C FOR WIN32 SINCE 1995.11.24.JUNSATO.*/
/*LAST CHANGE:1997.6.28.*/
 
/*ENABLE CREATION OF FRAME "CHANGE","REFER".....NOT YET.*/
/*TRANSLATION OF INPUTFILE ORGAN INTO ARCLM UNAVAILABLE.*/
/*SLAB DIVISION FOR CMQ.....NOT YET.*/
/*SORT COMPONENTS.....NOT YET.*/
/*IMPROVE "CROUTLUDECOMPOSITION".....NOT YET.*/
/*COMPONENT CONSISTS OF LINE,ROW,VALUE,DOWN.*/

/*CODES AT WILL.*/
/*FIXED LINES OF MATRIX INTO BOW.....INEFFECTIVE.*/
/*ALLOCATE ALL MEMORY FOR COMPONENTS FILL IN.....INEFFECTIVE.*/
/*INPUTFILE INTO MEMORY.*/
/*GLOBAL VECTOR INTO MEMORY.*/
/*CONFINEMENTS INTO MEMORY.*/
/*"GAUSSJORDANSWEEP" INTO "CROUTLUDECOMPOSITION".*/
/*"COEFFICIENTS","UPDATESTRESS" EXPONENT>1.0*/
/*GLOBAL MATRIX INTO MEMORY.*/
/*ONLY NONZERO COMPONENT IN GLOBAL MATRIX.*/

/*NODE:6 DOF.*/
/*LOAD:INCREMENTAL.*/
/*PLASTIC HINGE:PERFECT ELASTOPLASTIC HINGE ON EACH END.*/
/*YIELD SURFACE:Nz,Qx,Qy,Mz,Mx,My WITH EXPONENT>1.0.*/

/*"ARCLM INPUTFILE SPECIFICATION"*/
/*APPELATION*/
/*NNODE NELEM NSECT*/
/*ISECT E POI A Ixx Iyy J QxMAX.....MzMIN*/
/*INODE X Y Z*/
/*IELEM ISECT NODEI NODEJ COORDANGLE BOUNDARY... LONGSTRESS...*/
/*INODE CONFINEMENTTYPE... CONFINEMENTVALUE...*/
/*INODE DIRECTION LONGREACTION*/

/*PRINTER SETTINGS FOR MIFES 3.0.*/
/*FONT:TERMINAL 8 PTS*/
/*MARGIN LEFT:12 RIGHT:8 TOP:0 BOTTOM:0*/
/*LINE PITCH:NORMAL*/
/*ROW:DOUBLE*/
/*HEADER:%y.%m.%d %t "%F" PAGE%p*/
/*FOOTER:NONE.*/

#include <windows.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <malloc.h>
#include <math.h>
#include <limits.h>
#include <time.h>

#include "canhead.h"                    /*DEFINITION OF COMMAND ID.*/

/*#ifndef PI*/
#define PI 3.1415926535897932384
/*#endif*/

#define DSAFETY 0.0500                 /*INCREMENT OF SAFETY FACTOR*/
/*TEST FOR 1000 NODES:DSAFETY=0.05*/
/*TEST FOR 8 NODES INCREMENTAL:DSAFETY=0.0005 LAPS=40*/
#define LAPS 3                                           /*ALL LAPS*/
#define WAIT 0                                     /*WAIT TIME[sec]*/
#define EXPONENT 1.500                  /*EXPONENT OF YIELD SURFACE*/

#define DIRECTORY "c:\\cdocs\\frame3\\data\\"    /*DIRECTORY*/

#define GX 0
#define GY 1
#define GZ 2
#define EX 1 /*ANOTHER CASE 0*/
#define EY 2 /*             1*/
#define EZ 0 /*             2*/

#define RED   0
#define GREEN 1
#define BLUE  2

#define ROLEWEIGHT 1
#define ROLERIGID  2
#define ROLESTRESS 3

#define COLUMN 1 /*ELEMENT TYPE.*/
#define GIRDER 2
#define BEAM   3
#define BRACE  4
#define WALL   5
#define SLAB   6

#define AXONOMETRIC 0
#define PERSPECTIVE 1

struct gcomponent{unsigned short int m,n;               /*LINE,ROW.*/
                  double value;
                  struct gcomponent *down,*left;}; /*GLOBAL MATRIX.*/

struct rgbcolor{int r,g,b;};                           /*RGB COLOR.*/
struct nodecolor{struct rgbcolor code;};             /*NODE COLORS.*/
struct elemcolor{struct rgbcolor code,
                                 axis,
                                 line;};          /*ELEMENT COLORS.*/
struct orgcolor{struct rgbcolor code,
                                axis;};             /*ORGAN COLORS.*/

struct onode{long int code,loff;     /*CODE NUMBER,OFFSET POSITION.*/
             double d[3];};
struct oconf{signed char iconf;
             double value;}; /*ONLY FOR ANALYSIS.*/
struct oprop{long int code,loff;
             char *name;
             double hiju,E,poi;
             double rfle[3],rfra[3];}; /*R,G,B.*/
struct ofigs{long int code,loff;
             double area,Ixx,Iyy,Jzz; /*COEFFICIENTS.*/
             double thick;
             struct oprop *prop;}; /*SECTION FIGURES.*/
struct osect{long int code,loff;
             char *name;
             char dflag; /*DRAW FLAG.*/
             int nfig;
             double E,poi,area,Ixx,Iyy,Jzz; /*TOTAL COEFFICIENTS.*/
             double exp,fmax[6],fmin[6]; /*YIELD SURFACE.*/
             struct ofigs *figs;
             struct rgbcolor dcolor;};
struct obans{long int code,loff;
             int nnod;
             struct onode **nods;}; /*PLANE OF ELEMENT.*/
struct owire{long int code,loff;
             double cangle;
             signed char iconf[2][6];
             double stress[2][6];
             struct onode *(node[2]);
             struct osect *sect;}; /*WIRE ELEM FOR ARCLM001,101.*/
struct ofilm{long int code,loff;
             struct onode *(node[3]);
             struct osect *sect;}; /*FILM ELEM FOR ARCLM301.*/

struct oelem{long int code,loff;
             int role,type;
             int nnod,nban;
             double cangle;
             double initial[2][6],current[2][6]; /*STRESSES.*/
             struct onode **nods;
             signed char *bonds; /*6 BONDS OF EACH NODES.*/
             struct obans *bans;
             struct osect *sect;
             struct elemcolor color;};

struct organ{long int code,loff;
             long int nnode,nelem,nprop,nsect;
             struct onode *nodes;
             struct oconf *confs;
             struct oprop *props;
             struct osect *sects;
             struct oelem *elems;
             struct orgcolor color;};

struct snode{int loff;
             struct onode n;
             char str[256];};                     /*STRING ON NODE.*/

struct vector{double dc[3];};                       /*SIZE 24 BYTES*/
struct line{long int code,loff;
            int r,g,b;
            struct onode ends[2];};
struct plane{long int code,loff;
             struct onode nods[3];                  /*CROSSING DOTS*/
             struct vector nvec;                  /*VECTOR VERTICAL*/
             double a,b,c,d;/*ax+by+cz+d=0*/};     /*SIZE 128 BYTES*/
struct cmqdot{struct onode n;
              struct line *(parents[3]);};

struct viewrange{struct onode max,min;};              /*VIEW RANGE.*/
struct menuvisible{int draw,
                       error,
                       surface,
                       sectlist,
                       weight,
                       cmq,
                       horizon,
                       ftype,    /*FILE TYPE.*/
                       savetype; /*SAVE AS TYPE.*/
                  };                /*VISIBLE FLAG FOR OPTION MENU.*/
struct nodevisible{char code,
                        loads[6],
                        confs[6],
                        disps[3],             /*DISPLACEMENT X,Y,Z.*/
                        react[6];
                  };                       /*VISIBLE FLAG FOR NODE.*/
struct elemvisible{char code,
                        axis,
                        hinge,
                        sectioncode,
                        deformation,
                        stress[6],             /*Nz,Qx,Qy,Mz,Mx,My.*/
                        cmqline;               /*CMQ DIVISION LINE.*/
                  };                    /*VISIBLE FLAG FOR ELEMENT.*/
struct globvisible{char axis;
                   struct menuvisible mv;
                   struct nodevisible nv;
                   struct elemvisible ev;
                  };                     /*VISIBLE FLAG FOR GLOBAL.*/
struct drawparam{double gaxis,
                        eaxis,
                        dfact,
                        qfact,
                        mfact,
                        pitch,
                        hsize;
                }; /*PARAMETERS OF DRAWING FRAME.*/

struct viewparam{int type;           /*0:AXONOMETRICS 1:PERSPECTIVE*/
                 int Xo,Yo;                /*ORIGIN OF 2D DRAWINGS.*/
                 double theta,phi,r;
                 double vr,odv,ov[3],axono[3][3];     /*PROJECTION.*/
                 double gfactor;                  /*FACTOR OF SIZE.*/
                 struct onode focus;                 /*FOCUS POINT.*/
                 struct viewrange range;           /*RANGE VISIBLE.*/
                 struct globvisible vflag;         /*VISIBLE FLAGS.*/
                 struct drawparam dparam;        /*DRAW PARAMETERS.*/
                };                /*PARAMETERS OF PERSPECTIVE VIEW.*/

struct hoganparam{double arrowsize,notchsize;
                  char title[3][80];
                  char unit[3][80];
                  double scale[3];
                  double max[3];
                  double min[3];
                  double pitch[3];
                  char label[3][80];
                 };                       /*PARAMETERS OF HOGANSHI.*/

struct windowparams{
                    int code;
                    int nstring,nchilds;
                    int sstatus;                    /*SCROLL STATUS*/
                    char *classname;               /*NAME OF CLASS.*/
                    HWND hwnd;
                    HDC hdcD,hdcB,hdcC;     /*DISPLAY,BACK,COMPATI.*/
                    int lstatus,rstatus;            /*MOUSE STATUS.*/
                    struct windowparams *childs;   /*CHILD WINDOWS.*/

                    int tx,ty;             /*MESSAGE TEXT POSITION.*/
                    struct snode *strset;      /*STRINGS ON WINDOW.*/

                    /*RECT grect,lrect;*/              /*RECT SIZE.*/
                    RECT vbar,hbar;                   /*SCROLL BAR.*/

                    char inpfile[80],otpfile[80];      /*FILE NAME.*/
                    struct viewparam vparam;  /*3D VIEW PARAMETERS.*/
                    struct organ org;
                   };
struct winparamsreg{
                    int code;
                    int nwin;
                    struct windowparams **wp;
                   };                   /*REGISTRY OF WINDOWPARAMS.*/

struct arclmframe{long int code,loff;
                  char *appelation;
                  int nnode,nelem,nsect,nreact,nlaps;
                  double *eigenval,**eigenvec;
                  FILE *fdisp,*felem,*freact,*fsurface;
                  struct onode *nodes,*ninit;
                  struct osect *sects;
                  struct owire *elems;
                  struct oconf *confs;
                 }; /*ARCLM FRAME.*/

extern int globalstatus;
extern struct windowparams wdraw,wmenu,wmesg,wsurf;
extern struct arclmframe arc; /*GLOBAL ARCLM FRAME.*/
long int comps; /*COMPONENTS IN GLOBAL MATRIX.*/
extern struct oelem gelem;
/*=================================================================*/
/*REGULATED FOR ORGAN.*/
/*"INPUT":INPUT FROM TEXT FILE TO MEMORY.*/
/*"CURRENT":CURRENT DATA FROM BINARY FILE TO MEMORY.*/
/*"SELECT":SELECT FROM SCREEN BY MOUSE.*/
/*"GET":GET FROM OPTIONS DIALOG BOX.*/
/*"SET":SET WINDOW OBJECTS,SET TEXT INTO OPTIONS DIALOG BOX.*/

/*MATHEMATICS SUBROUTINES.*/
double *matrixvector(double **mtx,double *vct,int msize);
double **matrixmatrix(double **mtxhead,double **mtxtail,int msize);
double **matrixtranspose(double **mtx,int msize);

struct vector vectoraddition(struct vector vct1,
                             struct vector vct2);
struct vector vectorsubtraction(struct vector vct1,
                                struct vector vct2);
double lengthvector(struct vector vct);
double innerproduct(struct vector vct1,
                    struct vector vct2);
struct vector outerproduct(struct vector vct1,
                           struct vector vct2);

double distancedotdot(struct onode n1,struct onode n2);
double distancedotline(struct onode n,struct line l,
                       struct onode *nearest);
double distancedotplane(struct onode n,struct plane p,
                        struct onode *nearest);
int distancelineline(struct line l1,
                     struct line l2,
                     double *nearest1,
                     double *nearest2,
                     double *distance);
double intersectlineplane(struct line l,struct plane p,
                          struct onode *ncross);
double intersectplaneplane(struct plane p1,struct plane p2,
                           struct line *lcross);
int divisionlineline(struct line l1,
                     struct line l2,
                     struct line *l0);

int intersectlineline(int x1,int y1,int x2,int y2,
                      int x3,int y3,int x4,int y4);
int insideban(struct obans b,struct onode n);
struct plane *bantoplane(struct obans b,struct plane *pl);

/*STRINGS SUBROUTINES.*/
char **fgetsbrk(FILE *fin,int *n);
void freestr(char **str,int nstr);
int fcancelstr(FILE *fin,const int nline);

/*STRUCTURES SUBROUTINES.*/
struct onode setcoord(double x,double y,double z);
struct line setlineends(double x1,double y1,double z1,
                        double x2,double y2,double z2);

/*WINDOW SUBROUTINES.*/
void errormessage(char *str);
void setfontformat(HDC hdc,int h,int w,
                   char *fontname,int r,int g,int b);
void getclientsize(HWND hwnd,long int *w,long int *h);
void getwindowsize(HWND hwnd,long int *w,long int *h);
int overlayhdc(struct windowparams wp,DWORD dwrop);
HDC extendhdc(HDC hdc,long int w,long int h);

/*PROJECTION SUBROUTINES.*/
void getviewparam(HWND hdwnd,struct viewparam *vp);
void setviewparam(HWND hdwnd,struct viewparam vp);
void createviewdata(struct viewparam *vp);
int nodeprojection(struct onode ng,struct onode *np,
                   struct viewparam vp);
int nodeontoscreen(struct onode ng,int *ix,int *iy,
                   struct viewparam vp);
void drawglobaltext(HDC hdc,struct viewparam vp,
                    struct onode tn,char *str);
void drawglobaltextaligned(HDC hdc,struct viewparam vp,
                           struct onode tn,char *str,
                           int wtype,int htype);
void drawglobalnode(HDC hdc,struct viewparam vp,struct onode gn);
void drawgloballine(HDC hdc,struct viewparam vp,struct onode gn1,
                                                struct onode gn2);
void drawlinestruct(HDC hdc,struct viewparam vp,struct line l);
void drawglobalarrow(HDC hdc,struct viewparam vp,
                     struct onode gn1,
                     struct onode gn2,
                     double asize);
void drawcircleonglobaldot(HDC hdc,
                           struct viewparam vp,
                           struct onode center,double radius);
void fillglobalban(HDC hdc,struct viewparam vp,struct obans gb,
                   int r,int g,int b,
                   DWORD rop);
void drawglobalban(HDC hdc,struct viewparam vp,struct obans gb,
                   int r,int g,int b,
                   DWORD rop);
void drawelement(HDC hdc,struct viewparam vp,
                 struct oelem ge,
                 int mode);
void drawglobalaxis(HDC hdc,
                    struct viewparam vp,
                    int r,int g,int b);
void draworganization(HDC hdc,struct viewparam vp,
                      struct organ go,
                      int mode);
void draworgannodes(HDC hdc,struct viewparam vp,struct organ go);
struct onode transltog(double **tdrccos,
                       struct onode on,
                       struct onode ln);
struct onode transgtol(double **drccos,
                       struct onode on,
                       struct onode gn);
void drawtextonlocaldot(HDC hdc,struct viewparam vp,
                        double **tdrccos,
                        struct onode on,
                        struct onode ln,
                        char *str);
void drawlocalline(HDC hdc,struct viewparam vp,
                   double **tdrccos,
                   struct onode on,
                   struct onode ln1,
                   struct onode ln2);
void drawlocalarrow(HDC hdc,struct viewparam vp,
                    double **tdrccos,
                    struct onode on,
                    struct onode ln1,
                    struct onode ln2);
void drawlocalban(HDC hdc,struct viewparam vp,
                  double **tdrccos,
                  struct onode on,
                  struct obans lb);

int insiderange(struct onode n,struct viewrange r);

int getdirection(struct viewparam vpcurrent,
                 long int dx,long int dy,
                 struct onode focus,
                 double *length);

/*ARCLM WIREFRAME SUBROUTINES.*/
void drawarclmnodes(HDC hdc,struct viewparam vp,
                    struct arclmframe af,long int code);
void drawglobalwire(HDC hdc,struct viewparam vp,
                    struct arclmframe af,
                    struct owire elem,
                    int cred,int cgreen,int cblue,
                    int ered,int egreen,int eblue,
                    long int selectcode);
void drawwirestress(HDC hdc,struct viewparam vp,
                    struct arclmframe af,
                    struct owire elem);
void drawwireaxis(HDC hdc,
                  struct viewparam vp,
                  struct onode n1,struct onode n2,double cangle);
void drawarclmframe(HDC hdc,struct viewparam vp,
                    struct arclmframe af,
                    long int code,
                    int mode);
void drawrotatingarclm(HDC hdc,struct viewparam vp,
                       struct arclmframe af);
void drawyieldsurface(HDC hdc,struct viewparam vp,
                      int f1,int f2,int f3,FILE *fsfc);

struct osect *getsection(HWND hdwnd,struct osect *sects,int nsect,
                         long int code);
void setsection(HWND hdwnd,struct osect *sect);
struct onode *getnode(HWND hdwnd,struct arclmframe *af,
                      long int code);
void setnode(HWND hdwnd,struct onode *node,
                        struct oconf *confs);
void setnodeconf(HWND hdwnd,int idiconf,int idvconf,
                 struct oconf *conf);
struct onode *selectnode(struct viewparam vp,
                         struct arclmframe *af,POINT point);
struct owire *getelement(HWND hdwnd,struct arclmframe *af,
                         long int code);
void setelement(HWND hdwnd,struct owire *elem);
void setelemconf(HWND hdwnd,signed char iconf,int idiconf);
struct owire *selectelement(struct viewparam vp,
                            struct arclmframe *af,POINT point);

/*ORGANIZATION SUBROUTINES.*/
void initializeorganization(struct organ *org);
int inputorganization(FILE *fin,struct organ *org);
int saveorganization(FILE *fout,struct organ *org);
void freeorganization(struct organ *org);
void slabdivision(HDC hdc,struct viewparam vp,struct obans gban);

struct rgbcolor setrgbcolor(int r,int g,int b);

/*ANALYSIS SUBROUTINES.*/
clock_t laptime(char *comment,clock_t t0);
void getincrement(HWND hdwnd,int *laps,double *dsafety);
double croutludecomposition(struct gcomponent *gmtx,
                            double *gvct,struct oconf *confs,
                            long int msize);
void currentpivot(long int i,long int msize);
int gread(struct gcomponent *gmtx,
          long int i,long int j,double *data);
int gwrite(struct gcomponent *gmtx,
           long int i,long int j,double data);
void gfree(struct gcomponent *gmtx,long int nnode);
int vread(FILE *fvct,long int i,double *data);
int vwrite(FILE *fvct,long int i,double *data);

int arclm101(struct arclmframe *af);
DWORD availablephysicalmemory(char *comment);
FILE *fgetstofopen(const char *directory,const char *mode,
                   int dlgitem);
void inputtexttomemory(FILE *ftext,struct arclmframe *af);
int saveasarclm(struct arclmframe *af);
void initialform(struct onode *nodes,FILE *fdisp,int nnode);
int initialnode(struct onode *nodes,int nnode,int code,
                struct onode *node);
void initialelem(struct owire *elems,FILE *felem,int nelem);
void initialreact(FILE *fin,FILE *freact,int nreact);
void inputinit(FILE *fin,int *nnode,int *nelem,int *nsect);
void inputnode(FILE *fdisp,struct onode *node);
void inputelem(struct owire *elems,FILE *felem,int offset,
               struct owire *elem);
void readsect(FILE *fin,struct osect *sect);
double **directioncosine(double x1,double y1,double z1,
                         double x2,double y2,double z2,
                         double cangle);
double **filmdrccos(struct onode n1,struct onode n2,struct onode n3);
double **transmatrix(/*struct owire elem,*/double **drccos);
double **assememtx(struct owire elem);
double **assemgmtx(struct owire elem,double *estress);
double **assempmtx(struct owire elem,double **estiff);
void coefficients(struct owire elem,double **estiff,
                  double f[],double q[][2],double a[][2]);
double **modifyhinge(struct owire elem,double **estiff);
double **transformation(double **estiff,double **tmatrix);
void assemgstiffness(struct gcomponent *gmtx,
                     double **estiff,
                     struct owire *elem);
void assemconf(struct oconf *confs,double *gvct,double dsafety,
               int nnode);
void modifygivend(struct gcomponent *gmtx,double *gvct,
                  struct oconf *confs,int nnode);
double *extractdisplacement(struct owire elem,double *gvct);
double *elemstress(struct owire *elem,
                   double *gvct,FILE *felem,FILE *fout);
void updatestress(FILE *felem,FILE *fout,
                  double *edisp,double *dstress,double **estiff,
                  struct owire *elem);
void outputdisp(double *gvct,FILE *fout,int nnode,
                struct onode *nodes);
void updateform(FILE *fdisp,double *gvct,int nnode);
void copyform(struct arclmframe *af,double *gvct);
void outputstress(struct owire elem,double *estress,FILE *fout);
void outputreaction(struct gcomponent *gmtx,
                    double *gvct,
                    struct onode *nodes,
                    struct oconf *confs,
                    FILE *freact,FILE *fout,int nnode);

/*HOGANSHI SUBROUTINES.*/
void gethoganparam(HWND hdwnd,struct hoganparam *hp);
void drawhoganaxis(HDC hdc,
                   struct viewparam vp,
                   struct hoganparam hp,
                   int r,int g,int b);
int drawhoganlines(HDC hdc,HWND hdwnd,struct viewparam vp,
                   FILE *fin);

/*ORGAN CREATE SUBROUTINES.*/
struct oconf *addconf(struct oconf conf[6],
                      long int nodeoffset,
                      struct organ *orgbefore);
struct oconf *changeconf(struct oconf conf[6],
                         long int nodeoffset,
                         struct oconf *cinit);
struct oconf *deleteconf(long int nodeoffset,
                         struct organ *orgbefore);
struct onode *addnode(struct onode node,
                      struct organ *orgbefore);
struct onode *deletenode(long int nodeoffset,
                         struct organ *orgbefore);
struct oelem *addelement(struct oelem *elem,
                         struct organ *orgbefore);
struct oelem *deleteelement(long int elemoffset,
                            struct organ *orgbefore);
void freeelement(struct oelem *elem,
                 char fnods,char fbonds,char fbans,char fbannods);
int copyelement(struct oelem *efrom,struct oelem *eto,
                struct organ *org);
struct onode *createnodeonplane(struct viewparam vp,
                                double mx,double my,
                                struct plane pl,
                                struct onode *ncross);
struct onode *createorgannode(long int code,
                              struct viewparam vp,
                              long int mx,long int my,
                              struct organ *org);
struct onode findlastnode(long int code,
                          struct viewparam vp,
                          long int mx,long int my,
                          struct plane *pe,
                          struct organ *org);
void chaseelement(struct oelem *elem,
                  struct viewparam vp,
                  struct windowparams wp);


double *matrixvector(double **mtx,double *vct,int msize)
/*MULTIPLY VECTOR BY MATRIX.*/
{
  int i,j;
  double *v;

  v=(double *)malloc(msize*sizeof(double));
  if(v==NULL) return NULL;
  for(i=0;i<=msize-1;i++)
  {
    *(v+i)=0.0;
    for(j=0;j<=msize-1;j++)
    {*(v+i)+=(*(*(mtx+i)+j))*(*(vct+j));}
  }

  return v;
}/*matrixvector*/

double **matrixmatrix(double **mtxhead,double **mtxtail,int msize)
/*MULTIPLY MATRIX BY MATRIX.*/
{
  int i,j,jj;
  double **mtx,*mline;

  mtx=(double **)malloc(msize*sizeof(double *));
  if(mtx==NULL) return NULL;
  for(i=0;i<=msize-1;i++)
  {
    mline=(double *)malloc(msize*sizeof(double));
    if(mline==NULL) return NULL;
    for(j=0;j<=msize-1;j++)
    {
      *(mline+j)=0.0;
      for(jj=0;jj<=msize-1;jj++)
      {*(mline+j)+=(*(*(mtxhead+i)+jj))*(*(*(mtxtail+jj)+j));}
    }
    *(mtx+i)=mline;
  }

  return mtx;
}/*matrixmatrix*/

double **matrixtranspose(double **mtx,int msize)
/*MATRIX TRANSPOSITION.*/
{
  double **m;
  double *mm;
  int i,j;

  m=(double **)malloc(msize*sizeof(double *));
  if(m==NULL) return NULL;
  for(i=1;i<=msize;i++)
  {
    mm=(double *)malloc(msize*sizeof(double));
    if(mm==NULL) return NULL;
    for(j=1;j<=msize;j++)
    {*(mm+j-1)=*(*(mtx+j-1)+i-1);}
    *(m+i-1)=mm;
  }

  return m;
}/*matrixtranspose*/

struct vector vectoraddition(struct vector vct1,
                             struct vector vct2)
/*RETURN:ADDITION OF VECTOR.*/
{
  struct vector add;

  add.dc[GX]=vct1.dc[GX]+vct2.dc[GX];
  add.dc[GY]=vct1.dc[GY]+vct2.dc[GY];
  add.dc[GZ]=vct1.dc[GZ]+vct2.dc[GZ];

  return add;
}/*vectoraddition*/

struct vector vectorsubtraction(struct vector vct1,
                                struct vector vct2)
/*RETURN:SUBTRACTION OF VECTOR VCT2-VCT1.*/
{
  struct vector sbt;

  sbt.dc[GX]=vct2.dc[GX]-vct1.dc[GX];
  sbt.dc[GY]=vct2.dc[GY]-vct1.dc[GY];
  sbt.dc[GZ]=vct2.dc[GZ]-vct1.dc[GZ];

  return sbt;
}/*vectorsubtraction*/

double lengthvector(struct vector vct)
/*RETURN:LENGTH OF VECTOR.*/
{
  double length;

  length=sqrt(vct.dc[GX]*vct.dc[GX]
             +vct.dc[GY]*vct.dc[GY]
             +vct.dc[GZ]*vct.dc[GZ]);

  return length;
}/*lengthvector*/

double innerproduct(struct vector vct1,
                    struct vector vct2)
/*RETURN:INNERPRODUCT {vct1} AND {vct2}.*/
{
  double value;

  value=vct1.dc[GX]*vct2.dc[GX]
       +vct1.dc[GY]*vct2.dc[GY]
       +vct1.dc[GZ]*vct2.dc[GZ];

  return value;
}/*innerproduct*/

struct vector outerproduct(struct vector vct1,
                           struct vector vct2)
/*RETURN:OUTERPRODUCT {vct1}X{vct2}.*/
{
  struct vector vct;

  vct.dc[GX]=vct1.dc[GY]*vct2.dc[GZ]-vct1.dc[GZ]*vct2.dc[GY];
  vct.dc[GY]=vct1.dc[GZ]*vct2.dc[GX]-vct1.dc[GX]*vct2.dc[GZ];
  vct.dc[GZ]=vct1.dc[GX]*vct2.dc[GY]-vct1.dc[GY]*vct2.dc[GX];

  return vct;
}/*outerproduct*/

double distancedotdot(struct onode n1,struct onode n2)
/*RETURN:DISTANCE OF 2 DOTS.*/
{
  double dx,dy,dz,distance;

  dx=n2.d[GX]-n1.d[GX];
  dy=n2.d[GY]-n1.d[GY];
  dz=n2.d[GZ]-n1.d[GZ];

  distance=sqrt(dx*dx+dy*dy+dz*dz);

  return distance;
}/*distancedotdot*/

double distancedotline(struct onode n,struct line l,
                       struct onode *nearest)
/*RETURN:DISTANCE FROM DOT TO LINE,NEAREST DOT.*/
{
  double inner,htlength,otlength,fact,distance;
  struct vector vctht;          /*VECTOR FROM HEAD TO TAIL OF LINE.*/
  struct vector vcthd;           /*VECTOR FROM HEAD OF LINE TO DOT.*/
  struct vector outer;                             /*OUTER PRODUCT.*/

  vctht.dc[GX]=l.ends[1].d[GX]-l.ends[0].d[GX];
  vctht.dc[GY]=l.ends[1].d[GY]-l.ends[0].d[GY];
  vctht.dc[GZ]=l.ends[1].d[GZ]-l.ends[0].d[GZ];

  htlength=lengthvector(vctht);
  if(htlength<=0.0)
  {
    nearest->d[GX]=0.0;
    nearest->d[GY]=0.0;
    nearest->d[GZ]=0.0;
    return 0.0;
  }

  vcthd.dc[GX]=n.d[GX]-l.ends[0].d[GX];
  vcthd.dc[GY]=n.d[GY]-l.ends[0].d[GY];
  vcthd.dc[GZ]=n.d[GZ]-l.ends[0].d[GZ];

  inner=innerproduct(vctht,vcthd);
  outer=outerproduct(vctht,vcthd);

  otlength=lengthvector(outer);

  distance=otlength/htlength;

  fact=inner/htlength;
  nearest->d[GX]=l.ends[0].d[GX]+fact*vctht.dc[GX];
  nearest->d[GY]=l.ends[0].d[GY]+fact*vctht.dc[GY];
  nearest->d[GZ]=l.ends[0].d[GZ]+fact*vctht.dc[GZ];

  return distance;
}/*distancedotline*/

double distancedotplane(struct onode n,struct plane p,
                        struct onode *nearest)
/*RETURN:DISTANCE FROM DOT TO PLANE,NEAREST DOT.*/
{
  double func,fact,nx,ny,nz,dn,distance;

  func=p.a*n.d[GX]+p.b*n.d[GY]+p.c*n.d[GZ]+p.d;   /*FUNCTION VALUE.*/

  nx=p.nvec.dc[GX];
  ny=p.nvec.dc[GY];
  nz=p.nvec.dc[GZ];

  dn=sqrt(nx*nx+ny*ny+nz*nz); /*LENGTH OF VECTOR.*/

  distance=func/dn; /*DISTANCE.*/

  fact=distance/dn; /*FACTOR OF VECTOR.*/

  nearest->d[GX]=n.d[GX]-fact*nx; /*NEAREST POINT.*/
  nearest->d[GY]=n.d[GY]-fact*ny;
  nearest->d[GZ]=n.d[GZ]-fact*nz;

  return distance;
}/*distancedotplane*/

int distancelineline(struct line l1,
                     struct line l2,
                     double *nearest1,
                     double *nearest2,
                     double *distance)
/*RETURN:FLAG OF SUCCESS,PARAMETER OF NEAREST DOTS.*/
{
  double inner10,inner20,inner12;
  double htlength1,htlength2,otlength;
  struct vector vctht1,vctht2;  /*VECTOR FROM HEAD TO TAIL OF LINE.*/
  struct vector vcthh;            /*VECTOR FROM HEAD L1 TO HEAD L2.*/
  struct vector outer;                             /*OUTER PRODUCT.*/
  struct onode n1,n2;

  vctht1.dc[GX]=l1.ends[1].d[GX]-l1.ends[0].d[GX];
  vctht1.dc[GY]=l1.ends[1].d[GY]-l1.ends[0].d[GY];
  vctht1.dc[GZ]=l1.ends[1].d[GZ]-l1.ends[0].d[GZ];

  vctht2.dc[GX]=l2.ends[1].d[GX]-l2.ends[0].d[GX];
  vctht2.dc[GY]=l2.ends[1].d[GY]-l2.ends[0].d[GY];
  vctht2.dc[GZ]=l2.ends[1].d[GZ]-l2.ends[0].d[GZ];

  vcthh.dc[GX]=l2.ends[0].d[GX]-l1.ends[0].d[GX];
  vcthh.dc[GY]=l2.ends[0].d[GY]-l1.ends[0].d[GY];
  vcthh.dc[GZ]=l2.ends[0].d[GZ]-l1.ends[0].d[GZ];

  inner10=innerproduct(vctht1,vcthh);
  inner20=innerproduct(vctht2,vcthh);
  inner12=innerproduct(vctht1,vctht2);
  outer=outerproduct(vctht1,vctht2);

  if(innerproduct(vcthh,outer)==0.0)
  {
    if(outer.dc[GZ]!=0.0)
    {
      otlength=outer.dc[GZ];
      htlength1=vcthh.dc[GX]*vctht2.dc[GY]
               -vcthh.dc[GY]*vctht2.dc[GX];
      htlength2=vcthh.dc[GX]*vctht1.dc[GY]
               -vcthh.dc[GY]*vctht1.dc[GX];
    }
    else if(outer.dc[GX]!=0.0)
    {
      otlength=outer.dc[GX];
      htlength1=vcthh.dc[GY]*vctht2.dc[GZ]
               -vcthh.dc[GZ]*vctht2.dc[GY];
      htlength2=vcthh.dc[GY]*vctht1.dc[GZ]
               -vcthh.dc[GZ]*vctht1.dc[GY];
    }
    else if(outer.dc[GY]!=0.0)
    {
      otlength=outer.dc[GY];
      htlength1=vcthh.dc[GZ]*vctht2.dc[GX]
               -vcthh.dc[GX]*vctht2.dc[GZ];
      htlength2=vcthh.dc[GZ]*vctht1.dc[GX]
               -vcthh.dc[GX]*vctht1.dc[GZ];
    }
    else /*PARALLEL*/
    {
      *distance=distancedotline(l1.ends[0],l2,&n2);
      *nearest1=0.0;
      if(vctht2.dc[GX]!=0.0)
      {
        *nearest2=(n2.d[GX]-l2.ends[0].d[GX])/vctht2.dc[GX];
      }
      else if(vctht2.dc[GY]!=0.0)
      {
        *nearest2=(n2.d[GY]-l2.ends[0].d[GY])/vctht2.dc[GY];
      }
      else if(vctht2.dc[GZ]!=0.0)
      {
        *nearest2=(n2.d[GZ]-l2.ends[0].d[GZ])/vctht2.dc[GZ];
      }
      else
      {
        *nearest2=0.0;
        return 0;
      }
      return 2;
    }

    *nearest1=htlength1/otlength;
    *nearest2=htlength2/otlength;
  }
  else
  {
    otlength=lengthvector(outer);
    htlength1=lengthvector(vctht1);
    htlength2=lengthvector(vctht2);
    if(otlength<=0.0 || htlength1<=0.0 || htlength2<=0.0)
    {
      *nearest1=0.0;
      *nearest2=0.0;
      *distance=0.0;
      return 0;
    }

    *nearest1=(inner12*inner20-htlength2*htlength2*inner10)
             /(otlength*otlength);
    *nearest2=(inner12*inner10-htlength1*htlength1*inner20)
             /(otlength*otlength);
  }

  n1.d[GX]=l1.ends[0].d[GX]+(*nearest1)*vctht1.dc[GX];
  n1.d[GY]=l1.ends[0].d[GY]+(*nearest1)*vctht1.dc[GY];
  n1.d[GZ]=l1.ends[0].d[GZ]+(*nearest1)*vctht1.dc[GZ];

  n2.d[GX]=l2.ends[0].d[GX]+(*nearest2)*vctht2.dc[GX];
  n2.d[GY]=l2.ends[0].d[GY]+(*nearest2)*vctht2.dc[GY];
  n2.d[GZ]=l2.ends[0].d[GZ]+(*nearest2)*vctht2.dc[GZ];

  *distance=distancedotdot(n1,n2);

  return 1;
}/*distancelineline*/

double intersectlineplane(struct line l,struct plane p,
                          struct onode *ncross)
/*RETURN:DISTANCE FROM LINE TO PLANE,INTERSECTION DOT.*/
{
  double inner,distance,x1,y1,z1,x2,y2,z2;
  struct vector vctht;          /*VECTOR FROM HEAD TO TAIL OF LINE.*/
  struct onode node;

  x1=l.ends[0].d[GX];
  y1=l.ends[0].d[GY];
  z1=l.ends[0].d[GZ];
  x2=l.ends[1].d[GX];
  y2=l.ends[1].d[GY];
  z2=l.ends[1].d[GZ];

  vctht.dc[GX]=x2-x1;
  vctht.dc[GY]=y2-y1;
  vctht.dc[GZ]=z2-z1;

  inner=innerproduct(vctht,p.nvec);

  if(inner==0.0)
  {
    distance=distancedotplane(l.ends[0],p,&node);
    ncross->code=0;
    return distance;
  }

  ncross->code=1;
  ncross->d[GX]=(p.b*(y2*x1-y1*x2)+p.c*(z2*x1-z1*x2)-p.d*(x2-x1))
                /inner;
  ncross->d[GY]=(p.c*(z2*y1-z1*y2)+p.a*(x2*y1-x1*y2)-p.d*(y2-y1))
                /inner;
  ncross->d[GZ]=(p.a*(x2*z1-x1*z2)+p.b*(y2*z1-y1*z2)-p.d*(z2-z1))
                /inner;

  return 0.0;
}/*intersectlineplane*/

double intersectplaneplane(struct plane p1,struct plane p2,
                           struct line *lcross)
/*RETURN:DISTANCE FROM PLANE TO PLANE,INTERSECTION LINE.*/
{
  double distance,ad,bd,cd;
  double eps=1.0E-4;
  struct vector vctht;          /*VECTOR FROM HEAD TO TAIL OF LINE.*/
  struct onode node;

  vctht=outerproduct(p1.nvec,p2.nvec);

  ad=p1.a*p2.d-p2.a*p1.d;
  bd=p1.b*p2.d-p2.b*p1.d;
  cd=p1.c*p2.d-p2.c*p1.d;

  if(vctht.dc[GX]>eps || vctht.dc[GX]<(-eps))
  {
    lcross->ends[0].d[GX]=0.0;
    lcross->ends[0].d[GY]= cd/vctht.dc[GX];
    lcross->ends[0].d[GZ]=-bd/vctht.dc[GX];
  }
  else if(vctht.dc[GY]>eps || vctht.dc[GY]<(-eps))
  {
    lcross->ends[0].d[GX]=-cd/vctht.dc[GY];
    lcross->ends[0].d[GY]=0.0;
    lcross->ends[0].d[GZ]= ad/vctht.dc[GY];
  }
  else if(vctht.dc[GZ]>eps || vctht.dc[GZ]<(-eps))
  {
    lcross->ends[0].d[GX]= bd/vctht.dc[GZ];
    lcross->ends[0].d[GY]=-ad/vctht.dc[GZ];
    lcross->ends[0].d[GZ]=0.0;
  }
  else
  {
    distance=distancedotplane(p1.nods[0],p2,&node);
    lcross->code=0;
    return distance;
  }

  lcross->code=1;
  lcross->ends[1].d[GX]=lcross->ends[0].d[GX]+vctht.dc[GX];
  lcross->ends[1].d[GY]=lcross->ends[0].d[GY]+vctht.dc[GY];
  lcross->ends[1].d[GZ]=lcross->ends[0].d[GZ]+vctht.dc[GZ];

  return 0.0;
}/*intersectplaneplane*/

int divisionlineline(struct line l1,
                     struct line l2,
                     struct line *l0)
/*RETURN:FLAG OF SUCCESS,LINE OF DIVISION.*/
{
  double htlength1,htlength2;
  double param1,param2,distance;
  struct vector vctht1,vctht2;  /*VECTOR FROM HEAD TO TAIL OF LINE.*/
  struct onode n1,n2;

  if(!distancelineline(l1,l2,&param1,&param2,&distance))
  {
    *l0=l1;
    return 1;
  }

  vctht1.dc[GX]=l1.ends[1].d[GX]-l1.ends[0].d[GX];
  vctht1.dc[GY]=l1.ends[1].d[GY]-l1.ends[0].d[GY];
  vctht1.dc[GZ]=l1.ends[1].d[GZ]-l1.ends[0].d[GZ];

  vctht2.dc[GX]=l2.ends[1].d[GX]-l2.ends[0].d[GX];
  vctht2.dc[GY]=l2.ends[1].d[GY]-l2.ends[0].d[GY];
  vctht2.dc[GZ]=l2.ends[1].d[GZ]-l2.ends[0].d[GZ];

  n1.d[GX]=l1.ends[0].d[GX]+param1*vctht1.dc[GX];
  n1.d[GY]=l1.ends[0].d[GY]+param1*vctht1.dc[GY];
  n1.d[GZ]=l1.ends[0].d[GZ]+param1*vctht1.dc[GZ];

  n2.d[GX]=l2.ends[0].d[GX]+param2*vctht2.dc[GX];
  n2.d[GY]=l2.ends[0].d[GY]+param2*vctht2.dc[GY];
  n2.d[GZ]=l2.ends[0].d[GZ]+param2*vctht2.dc[GZ];

  htlength1=lengthvector(vctht1);
  htlength2=lengthvector(vctht2);

  vctht1.dc[GX]/=htlength1; /*UNIT.*/
  vctht1.dc[GY]/=htlength1;
  vctht1.dc[GZ]/=htlength1;

  vctht2.dc[GX]/=htlength2;
  vctht2.dc[GY]/=htlength2;
  vctht2.dc[GZ]/=htlength2;

  l0->ends[0].d[GX]=0.5*(n1.d[GX]+n2.d[GX]); /*HEAD.*/
  l0->ends[0].d[GY]=0.5*(n1.d[GY]+n2.d[GY]);
  l0->ends[0].d[GZ]=0.5*(n1.d[GZ]+n2.d[GZ]);

  l0->ends[1].d[GX]=(l0->ends[0].d[GX])            /*TAIL.*/
                   +(vctht1.dc[GX]+vctht2.dc[GX]);
  l0->ends[1].d[GY]=(l0->ends[0].d[GY])
                   +(vctht1.dc[GY]+vctht2.dc[GY]);
  l0->ends[1].d[GZ]=(l0->ends[0].d[GZ])
                   +(vctht1.dc[GZ]+vctht2.dc[GZ]);

  return 1;
}/*divisionlineline*/

int intersectlineline(int x1,int y1,int x2,int y2,
                      int x3,int y3,int x4,int y4)
/*INTERSECTION OF 2 LINES.*/
{
  double n1,n2,n3,n4;

  n1=(double)((x2-x1)*(y3-y1)-(y2-y1)*(x3-x1));
  n2=(double)((x2-x1)*(y4-y1)-(y2-y1)*(x4-x1));
  n3=(double)((x4-x3)*(y1-y3)-(y4-y3)*(x1-x3));
  n4=(double)((x4-x3)*(y2-y3)-(y4-y3)*(x2-x3));

  if(n1==0.0 && n2==0.0) /*ON LINE.*/
  {
    n1=(double)((x3-x1)*(x4-x1));
    n2=(double)((y3-y1)*(y4-y1));
    n3=(double)((x3-x2)*(x4-x2));
    n4=(double)((y3-y2)*(y4-y2));

    if((n1*n2)<0.0 || (n3*n4)<0.0) return 1;
    else return 0;
  }
  else if((n1*n2)>0.0) return 0;                 /*NOT INTERSECTED.*/
  else if((n3*n4)>0.0) return 0;                 /*NOT INTERSECTED.*/
  else                 return 1;                     /*INTERSECTED.*/

}/*intersection*/

int insideban(struct obans b,struct onode n)
/*RETURN:1=NODE INSIDE BAN COLUMN,0=OUTSIDE.*/
{
  /*char non[256];*/
  int i,count;
  double fact,xo,yo,xi,yi,xj,yj;

  xo=n.d[GX];
  yo=n.d[GY];

  count=0;
  for(i=0;i<b.nnod;i++)
  {
    xi=(*(b.nods+i))->d[GX];
    yi=(*(b.nods+i))->d[GY];

    if(i==b.nnod-1)
    {
      xj=(*(b.nods+0))->d[GX];
      yj=(*(b.nods+0))->d[GY];
    }
    else
    {
      xj=(*(b.nods+i+1))->d[GX];
      yj=(*(b.nods+i+1))->d[GY];
    }

    if(yi!=yj)
    {
/*sprintf(non,"(%.0f %.0f)(%.0f %.0f)(%.0f %.0f)",xo,yo,xi,yi,xj,yj);
MessageBox(NULL,non,"Inside",MB_OK);*/
      if(!(yo<yi && yo<yj) && !(yo>yi && yo>yj))
      {
/*MessageBox(NULL,"GetIn 2.","Inside",MB_OK);*/
        /*INTERSECTION WITH HORIZONTAL LINE.*/
        fact=(xj-xi)*(yo-yi)/(yj-yi)-(xo-xi);
        if(fact>=0.0) count++;
      }
    }
  }

  if((count%2)==1) return 1; /*KISUU*/
  else             return 0; /*GUSUU*/
}/*insideban*/

struct plane *bantoplane(struct obans b,struct plane *pl)
{
  struct vector v01,v02;

  if(b.nnod<3) return NULL;

  pl->nods[0].d[GX]=(*(b.nods+0))->d[GX];
  pl->nods[0].d[GY]=(*(b.nods+0))->d[GY];
  pl->nods[0].d[GZ]=(*(b.nods+0))->d[GZ];

  pl->nods[1].d[GX]=(*(b.nods+1))->d[GX];
  pl->nods[1].d[GY]=(*(b.nods+1))->d[GY];
  pl->nods[1].d[GZ]=(*(b.nods+1))->d[GZ];

  pl->nods[2].d[GX]=(*(b.nods+b.nnod-1))->d[GX];
  pl->nods[2].d[GY]=(*(b.nods+b.nnod-1))->d[GY];
  pl->nods[2].d[GZ]=(*(b.nods+b.nnod-1))->d[GZ];

  v01.dc[GX]=(pl->nods[1].d[GX])-(pl->nods[0].d[GX]);
  v01.dc[GY]=(pl->nods[1].d[GY])-(pl->nods[0].d[GY]);
  v01.dc[GZ]=(pl->nods[1].d[GZ])-(pl->nods[0].d[GZ]);

  v02.dc[GX]=(pl->nods[2].d[GX])-(pl->nods[0].d[GX]);
  v02.dc[GY]=(pl->nods[2].d[GY])-(pl->nods[0].d[GY]);
  v02.dc[GZ]=(pl->nods[2].d[GZ])-(pl->nods[0].d[GZ]);

  pl->nvec=outerproduct(v01,v02);
  /*pl->nvec=outerproduct(v02,v01);*/

  pl->a=pl->nvec.dc[GX];
  pl->b=pl->nvec.dc[GY];
  pl->c=pl->nvec.dc[GZ];
  pl->d=-(pl->nvec.dc[GX])*(pl->nods[0].d[GX])
        -(pl->nvec.dc[GY])*(pl->nods[0].d[GY])
        -(pl->nvec.dc[GZ])*(pl->nods[0].d[GZ]);

  return pl;
}/*bantoplane*/

char **fgetsbrk(FILE *fin,int *n)
/*ONE LINE FROM FILE BREAK INTO SEVERAL STRINGS.*/
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
  while(1)
  {
    ii=0;
    ss=NULL;
    while((ichr=getc(f))==' ')
    ;
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
      else if(ichr==' ') {*(ss+ii)='\0'; break;}
      else {*(ss+ii)=(char)ichr;}

      ii++;
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
}/*fgetsbrk*/

void freestr(char **str,int nstr)
/*FREE STRINGS ALLOCATED BY "fgetsbrk".*/
{
  for(;nstr>0;nstr--) free(*(str+nstr-1));
  free(str);

  return;
}/*freestr*/

int fcancelstr(FILE *fin,const int nline)
/*SKIP LINES FROM FILE.*/
{
  FILE *f=fin;
  int i,chr;

  for(i=0;i<nline;)
  {
    chr=getc(f);
    if(chr=='\n') i++;
    if(chr==EOF){i++; break;}
  }
  return i;
}/*fcancelstr*/

struct onode setcoord(double x,double y,double z)
{
  struct onode n;

  n.d[0]=x;
  n.d[1]=y;
  n.d[2]=z;

  return n;
}/*setcoord*/

struct line setlineends(double x1,double y1,double z1,
                        double x2,double y2,double z2)
{
  struct line l;

  l.ends[0].d[0]=x1;
  l.ends[0].d[1]=y1;
  l.ends[0].d[2]=z1;

  l.ends[1].d[0]=x2;
  l.ends[1].d[1]=y2;
  l.ends[1].d[2]=z2;

  return l;
}/*setlineends*/

void errormessage(char *str)
/*OUTPUT ERROR MESSAGE.*/
{
  HWND hoya;
  SIZE size;
  long int x,y,tx,ty,cw,ch,pw,ph;
  POINT pc,pp;

  if(wmesg.nchilds>=2 && (wmesg.childs+1)->hwnd!=NULL)
  {
    tx=(wmesg.childs+1)->tx;
    ty=(wmesg.childs+1)->ty;

    GetTextExtentPoint32((wmesg.childs+1)->hdcC,
                         str,strlen(str),&size);

    pc.x=0;
    pc.y=0;
    ClientToScreen((wmesg.childs+1)->hwnd,&pc); /*TOPLEFT OF CHILD.*/

    hoya=GetParent((wmesg.childs+1)->hwnd);
    pp.x=0;
    pp.y=0;
    ClientToScreen(hoya,&pp); /*TOPLEFT OF PARENT.*/

    x=pp.x-pc.x;
    y=pp.y-pc.y;

    getclientsize(hoya,&pw,&ph);
    getclientsize((wmesg.childs+1)->hwnd,&cw,&ch);

    if(ch<(ty+size.cy-2)) /*EXTEND HEIGHT IF FILLED UP.*/
    {
      ch=ty+size.cy-2;
      extendhdc((wmesg.childs+1)->hdcC,cw,ch);
      extendhdc((wmesg.childs+1)->hdcB,cw,ch);
      MoveWindow((wmesg.childs+1)->hwnd,-x,-y,cw,ch,TRUE);
    }
    if(cw<size.cx) /*EXTEND WIDTH FOR LONG TEXT.*/
    {
      cw=size.cx;
      extendhdc((wmesg.childs+1)->hdcC,cw,ch);
      extendhdc((wmesg.childs+1)->hdcB,cw,ch);
      MoveWindow((wmesg.childs+1)->hwnd,-x,-y,cw,ch,TRUE);
    }
    if((y+ph)<(ty+size.cy-2)) /*SCROLL IF OVERFLOW FROM PARENT.*/
    {
      MoveWindow((wmesg.childs+1)->hwnd,
                 (-x),(ph-(ty+size.cy-2)),cw,ch,TRUE);
    }
    SendMessage(wmesg.hwnd,WM_PAINT,0,0);

    SetTextColor((wmesg.childs+1)->hdcC,RGB(255,255,255));
    TextOut((wmesg.childs+1)->hdcC,tx,ty,str,strlen(str));

    overlayhdc(*(wmesg.childs+1),SRCPAINT);

    (wmesg.childs+1)->ty+=size.cy-2;
  }
  return;
}/*errormessage*/

void setfontformat(HDC hdc,int h,int w,
                   char *fontname,int r,int g,int b)
/*SET FORMAT OF FONT.*/
{
  HFONT hfont,pfont;

  if(hdc!=NULL)
  {
    hfont=CreateFont(h,w,                            /*HEIGHT,WIDTH*/
                     0,0,      /*TANGENT OF LINE,LETTER[0.1 DEGREE]*/
                     100,                                  /*WEIGHT*/
                     0,0,0,            /*ITALIC,UNDERLINE,STRIKEOUT*/
                     ANSI_CHARSET,                       /*CHAR SET*/
                     OUT_DEFAULT_PRECIS,
                     CLIP_DEFAULT_PRECIS,
                     DEFAULT_QUALITY,
                     DEFAULT_PITCH | FF_DONTCARE,
                     fontname);                         /*FONT NAME*/
    pfont=SelectObject(hdc,hfont);
    DeleteObject(pfont);

    SetBkMode(hdc,TRANSPARENT);               /*BACKGROUND OF TEXT.*/
    SetTextColor(hdc,RGB(r,g,b));                     /*TEXT COLOR.*/
  }
  return;
}/*setfontformat*/

void getclientsize(HWND hwnd,long int *w,long int *h)
/*GET WIDTH,HEIGHT OF CLIENT OF WINDOW.*/
{
  RECT rect;

  GetClientRect(hwnd,&rect);
  *w = rect.right - rect.left;
  *h = rect.bottom - rect.top;

  return;
}/*getclientsize*/

void getwindowsize(HWND hwnd,long int *w,long int *h)
/*GET WIDTH,HEIGHT OF WINDOW.*/
{
  RECT rect;

  GetWindowRect(hwnd,&rect);
  *w = rect.right - rect.left;
  *h = rect.bottom - rect.top;

  return;
}/*getwindowsize*/

int overlayhdc(struct windowparams wp,DWORD dwrop)
/*OVERLAY HDC COMPATI,BACK.*/
{
  HWND hoya;
  HDC hdc;
  POINT pc,pp;
  long int x,y,maxX,maxY;

  if(wp.hdcB != NULL && wp.hdcC != NULL)
  {
    pc.x=0;
    pc.y=0;
    ClientToScreen(wp.hwnd,&pc);

    hoya=GetParent(wp.hwnd);
    pp.x=0;
    pp.y=0;
    ClientToScreen(hoya,&pp);

    x=pp.x-pc.x;
    y=pp.y-pc.y;

    getwindowsize(hoya,&maxX,&maxY);

    PatBlt(wp.hdcB,x,y,maxX,maxY,PATCOPY);
    BitBlt(wp.hdcB,x,y,maxX,maxY,wp.hdcC,x,y,
           dwrop); /*GROUND BLACK:SRCPAINT WHITE:SRCAND.*/

    hdc = GetDC(wp.hwnd);
    BitBlt(hdc,x,y,maxX,maxY,wp.hdcB,x,y,SRCCOPY);
    ReleaseDC(wp.hwnd,hdc);

    return 1;
  }
  return 0;
}/*overlayhdc*/

HDC extendhdc(HDC hdc,long int w,long int h)
/*EXTEND HDC.*/
{
  HDC hdcC;
  HBITMAP hbit,pbit;
  BITMAP bmp;

  if(hdc==NULL) return NULL;

  hdcC=CreateCompatibleDC(hdc);
  hbit=CreateCompatibleBitmap(hdc,w,h);
  pbit=SelectObject(hdc,hbit);
  SelectObject(hdcC,pbit);

  GetObject(pbit,sizeof(BITMAP),(LPSTR)&bmp);

  PatBlt(hdc,0,0,w,h,PATCOPY);
  BitBlt(hdc,0,0,bmp.bmWidth,bmp.bmHeight,hdcC,0,0,SRCCOPY);

  DeleteObject(pbit);
  DeleteDC(hdcC);

  return hdc;
}/*extendhdc*/

void getviewparam(HWND hdwnd,struct viewparam *vp)
/*GET VIEW PARAMETERS FROM DIALOG.*/
{
  char data[80];

  /*vp->type=PERSPECTIVE;*/
  GetDlgItemText(hdwnd,IDV_GFACTOR,data,80); /*PROJECTION PARAMS.*/
  vp->gfactor=strtod(data,NULL);

  GetDlgItemText(hdwnd,IDV_X,data,80);
  vp->focus.d[0]=strtod(data,NULL);
  GetDlgItemText(hdwnd,IDV_Y,data,80);
  vp->focus.d[1]=strtod(data,NULL);
  GetDlgItemText(hdwnd,IDV_Z,data,80);
  vp->focus.d[2]=strtod(data,NULL);

  GetDlgItemText(hdwnd,IDV_THETA,data,80);
  vp->theta=strtod(data,NULL);
  GetDlgItemText(hdwnd,IDV_PHI,data,80);
  vp->phi=strtod(data,NULL);
  GetDlgItemText(hdwnd,IDV_R,data,80);
  vp->r=strtod(data,NULL);
  GetDlgItemText(hdwnd,IDV_L,data,80);
  vp->odv=strtod(data,NULL);

  GetDlgItemText(hdwnd,IDV_ORIGINX,data,80);
  vp->Xo=(int)strtol(data,NULL,10);
  GetDlgItemText(hdwnd,IDV_ORIGINY,data,80);
  vp->Yo=(int)strtol(data,NULL,10);

  GetDlgItemText(hdwnd,IDR_XMAX,data,80); /*RANGE DATA.*/
  vp->range.max.d[GX]=strtod(data,NULL);
  GetDlgItemText(hdwnd,IDR_XMIN,data,80);
  vp->range.min.d[GX]=strtod(data,NULL);
  GetDlgItemText(hdwnd,IDR_YMAX,data,80);
  vp->range.max.d[GY]=strtod(data,NULL);
  GetDlgItemText(hdwnd,IDR_YMIN,data,80);
  vp->range.min.d[GY]=strtod(data,NULL);
  GetDlgItemText(hdwnd,IDR_ZMAX,data,80);
  vp->range.max.d[GZ]=strtod(data,NULL);
  GetDlgItemText(hdwnd,IDR_ZMIN,data,80);
  vp->range.min.d[GZ]=strtod(data,NULL);

  GetDlgItemText(hdwnd,IDV_GAXISLENGTH,data,80);
  vp->dparam.gaxis=strtod(data,NULL);
  GetDlgItemText(hdwnd,IDV_EAXISLENGTH,data,80);
  vp->dparam.eaxis=strtod(data,NULL);
  GetDlgItemText(hdwnd,IDV_DFACTOR,data,80);
  vp->dparam.dfact=strtod(data,NULL);
  GetDlgItemText(hdwnd,IDV_QFACTOR,data,80);
  vp->dparam.qfact=strtod(data,NULL);
  GetDlgItemText(hdwnd,IDV_MFACTOR,data,80);
  vp->dparam.mfact=strtod(data,NULL);
  GetDlgItemText(hdwnd,IDV_GYOPITCH,data,80);
  vp->dparam.pitch=strtod(data,NULL);
  GetDlgItemText(hdwnd,IDV_HINGESIZE,data,80);
  vp->dparam.hsize=strtod(data,NULL);

  createviewdata(vp);
  return;
}/*getviewparam*/

void setviewparam(HWND hdwnd,struct viewparam vp)
/*SET VIEW PARAMETERS INTO DIALOG.*/
{
  char str[80];

  /*vp->type=PERSPECTIVE;*/

  sprintf(str,"%.1f",vp.gfactor);
  SetDlgItemText(hdwnd,IDV_GFACTOR,str);

  sprintf(str,"%.1f",vp.focus.d[0]);
  SetDlgItemText(hdwnd,IDV_X,str);
  sprintf(str,"%.1f",vp.focus.d[1]);
  SetDlgItemText(hdwnd,IDV_Y,str);
  sprintf(str,"%.1f",vp.focus.d[2]);
  SetDlgItemText(hdwnd,IDV_Z,str);

  sprintf(str,"%.1f",vp.theta);
  SetDlgItemText(hdwnd,IDV_THETA,str);
  sprintf(str,"%.1f",vp.phi);
  SetDlgItemText(hdwnd,IDV_PHI,str);
  sprintf(str,"%.1f",vp.r);
  SetDlgItemText(hdwnd,IDV_R,str);
  sprintf(str,"%.1f",vp.odv);
  SetDlgItemText(hdwnd,IDV_L,str);

  sprintf(str,"%d",(vp.Xo));
  SetDlgItemText(hdwnd,IDV_ORIGINX,str);
  sprintf(str,"%d",(vp.Yo));
  SetDlgItemText(hdwnd,IDV_ORIGINY,str);

  return;
}/*setviewparam*/

void createviewdata(struct viewparam *vp)
/*CREATE VIEW DATA.*/
{
  double vth,vph;

  vth=vp->theta*PI/180.0;                      /*DEGREE INTO RADIAN*/
  vph=vp->phi*PI/180.0;
  vp->vr=vp->r;
  vp->ov[0]=(vp->vr)*cos(vph)*cos(vth); /*VIEW POINT X.*/
  vp->ov[1]=(vp->vr)*cos(vph)*sin(vth); /*VIEW POINT Y.*/
  vp->ov[2]=(vp->vr)*sin(vph);         /*VIEW POINT Z.*/
  vp->axono[0][0]=cos(vph)*cos(vth);
  vp->axono[0][1]=cos(vph)*sin(vth);
  vp->axono[0][2]=sin(vph);
  vp->axono[1][0]=-sin(vth);
  vp->axono[1][1]=cos(vth);
  vp->axono[1][2]=0.0;
  vp->axono[2][0]=-sin(vph)*cos(vth);
  vp->axono[2][1]=-sin(vph)*sin(vth);
  vp->axono[2][2]=cos(vph);

  return;
}/*createviewdata*/

int nodeprojection(struct onode ng,struct onode *np,
                   struct viewparam vp)
/*PROJECT GLOBAL NODE BY AXONOMETRICS,PERSPECTIVE.*/
/*NP:PROJECTED NODE.*/
{
  int i;
  double ppers[3],pv[3],vnai,fact;

  ng.d[GX]-=vp.focus.d[GX];
  ng.d[GY]-=vp.focus.d[GY];
  ng.d[GZ]-=vp.focus.d[GZ];

  for(i=0;i<=2;i++)
  {
    ppers[i]=vp.axono[i][0]*ng.d[GX]
            +vp.axono[i][1]*ng.d[GY]
            +vp.axono[i][2]*ng.d[GZ];
    pv[i]=vp.ov[i]-ng.d[i];
  }

  if(vp.type==PERSPECTIVE)      /*TYPE 0:AXONOMETRICS 1:PERSPECTIVE*/
  {
    vnai=vp.ov[0]*pv[0]+vp.ov[1]*pv[1]+vp.ov[2]*pv[2];
    if(vnai==0.0) return 0;

    fact=vp.vr*vp.odv/vnai;
    ppers[1]*=fact;
    ppers[2]*=fact;
  }

  np->d[GX]=ppers[1];
  np->d[GY]=ppers[2];
  np->d[GZ]=0.0;

  return 1;
}/*nodeprojection*/

int nodeontoscreen(struct onode ng,int *ix,int *iy,
                   struct viewparam vp)
/*PROJECT GLOBAL NODE ONTO SCREEN AS INT.*/
{
  struct onode np;

  if(nodeprojection(ng,&np,vp))
  {
    *ix=vp.Xo+(int)(vp.gfactor*np.d[GX]);
    *iy=vp.Yo-(int)(vp.gfactor*np.d[GY]);
    return 1;
  }
  else
  {
    *ix=0;
    *iy=0;
    return 0;
  }
}/*nodeontoscreen*/

void drawglobaltext(HDC hdc,struct viewparam vp,
                    struct onode tn,char *str)
/*DRAW TEXT GLOBAL ON SCREEN BY AXONOMETRICS,PERSPECTIVE.*/
/*TEXTCOLOR DEFINITION MUST BE ALREADY DONE FOR FAST DRAWING.*/
{
  int Ox,Oy;                                   /*COORDINATION ARROW*/

  if(!nodeontoscreen(tn,&Ox,&Oy,vp)) return;           /*PROJECTION*/
  TextOut(hdc,Ox,Oy,str,strlen(str));

  return;
}/*drawglobaltext*/

void drawglobaltextaligned(HDC hdc,struct viewparam vp,
                           struct onode tn,char *str,
                           int wtype,int htype)
/*DRAW TEXT GLOBAL ON SCREEN ALIGNED.*/
{
  SIZE tsize;
  int Ox,Oy;                                   /*COORDINATION ARROW*/

  if(!nodeontoscreen(tn,&Ox,&Oy,vp)) return;           /*PROJECTION*/

  GetTextExtentPoint32(hdc,str,strlen(str),&tsize);

  if(wtype==TEXT_CENTER) Ox-=(int)(tsize.cx/2);
  if(wtype==TEXT_RIGHT)  Ox-=(int)(tsize.cx);

  if(htype==TEXT_CENTER) Oy-=(int)(tsize.cy/2);
  if(htype==TEXT_BOTTOM) Oy-=(int)(tsize.cy);

  TextOut(hdc,Ox,Oy,str,strlen(str));

  return;
}/*drawglobaltextaligned*/

void drawglobalnode(HDC hdc,struct viewparam vp,struct onode gn)
/*DRAW NODE GLOBAL ON SCREEN BY AXONOMETRICS,PERSPECTIVE.*/
/*TEXTCOLOR DEFINITION MUST BE ALREADY DONE FOR FAST DRAWING.*/
{
  char str[20];

  sprintf(str,"%ld",gn.code);
  drawglobaltext(hdc,vp,gn,str);

  return;
}/*drawglobalnode*/

void drawgloballine(HDC hdc,struct viewparam vp,struct onode gn1,
                                                struct onode gn2)
/*DRAW LINE ON SCREEN BY AXONOMETRICS,PERSPECTIVE.*/
{
  int Hx,Hy,Tx,Ty;

  if(!nodeontoscreen(gn1,&Hx,&Hy,vp)) return;          /*PROJECTION*/
  if(!nodeontoscreen(gn2,&Tx,&Ty,vp)) return;          /*PROJECTION*/

  MoveToEx(hdc,Hx,Hy,NULL);
  LineTo(hdc,Tx,Ty);

  return;
}/*drawgloballine*/

void drawlinestruct(HDC hdc,struct viewparam vp,struct line l)
/*DRAW LINE STRUCT ON SCREEN BY AXONOMETRICS,PERSPECTIVE.*/
{
  HPEN hpen,ppen;

  hpen=CreatePen(PS_SOLID,1,RGB(l.r,l.g,l.b));
  ppen=SelectObject(hdc,hpen);

  drawgloballine(hdc,vp,l.ends[0],l.ends[1]);

  SelectObject(hdc,ppen);
  DeleteObject(hpen);
  return;
}/*drawlinestruct*/

void drawglobalarrow(HDC hdc,struct viewparam vp,
                     struct onode gn1,
                     struct onode gn2,
                     double asize)
/*DRAW ARROW ON 2D SCREEN.*/
{
  int Hx,Hy,Tx,Ty;
  double dx,dy,arrowlength,arrowsize,arrowtan;

  if(!nodeontoscreen(gn1,&Hx,&Hy,vp)) return;         /*PROJECTION.*/
  if(!nodeontoscreen(gn2,&Tx,&Ty,vp)) return;         /*PROJECTION.*/

  MoveToEx(hdc,Hx,Hy,NULL);
  LineTo(hdc,Tx,Ty);                                        /*AXIS.*/

  dx=(double)(Tx-Hx);
  dy=(double)(Ty-Hy);
  arrowlength=sqrt(dx*dx+dy*dy);

  if(asize<=0.0) arrowsize=0.2*arrowlength;
  else           arrowsize=asize;

  arrowtan=0.1*arrowsize;

  if(arrowlength!=0.0)  /*ARROW.*/
  {
    Hx=Tx
       -(int)((arrowsize/arrowlength)*dx+(dy/arrowlength)*arrowtan);
    Hy=Ty
       -(int)((arrowsize/arrowlength)*dy-(dx/arrowlength)*arrowtan);
    MoveToEx(hdc,Hx,Hy,NULL);
    LineTo(hdc,Tx,Ty);

    Hx=Tx
       -(int)((arrowsize/arrowlength)*dx-(dy/arrowlength)*arrowtan);
    Hy=Ty
       -(int)((arrowsize/arrowlength)*dy+(dx/arrowlength)*arrowtan);
    LineTo(hdc,Hx,Hy);
  }
  return;
}/*drawglobalarrow*/

void drawcircleonglobaldot(HDC hdc,
                           struct viewparam vp,
                           struct onode center,double radius)
/*DRAW CIRCLE ON 3D POINT DRAWN BY AXONOMETRICS,PERSPECTIVE.*/
{
  int Ox,Oy,R;                                 /*POINT ON DRAWINGS.*/
  double inner,fact;
  struct onode pinit;
  struct vector gov,gpv;

  pinit.d[0]=center.d[GX];
  pinit.d[1]=center.d[GY];
  pinit.d[2]=center.d[GZ];
  if(!nodeontoscreen(pinit,&Ox,&Oy,vp)) return;        /*PROJECTION*/

  radius*=vp.gfactor;                           /*PROJECTED RADIUS.*/
  if(vp.type)
  {
    gov.dc[GX]=vp.ov[GX];
    gov.dc[GY]=vp.ov[GY];
    gov.dc[GZ]=vp.ov[GZ];

    gpv.dc[GX]=vp.ov[GX]-(center.d[GX]-vp.focus.d[GX]);
    gpv.dc[GY]=vp.ov[GY]-(center.d[GY]-vp.focus.d[GY]);
    gpv.dc[GZ]=vp.ov[GZ]-(center.d[GZ]-vp.focus.d[GZ]);

    inner=innerproduct(gov,gpv);

    fact=vp.vr*vp.odv/inner;
    radius*=fabs(fact);
  }
  R=(int)radius;

  Arc(hdc,(Ox-R),(Oy-R),(Ox+R),(Oy+R),0,0,0,0);     /*HOLLOW CIRCLE*/

  return;
}/*drawcircleonglobaldot*/

void fillglobalban(HDC hdc,struct viewparam vp,struct obans gb,
                   int r,int g,int b,
                   DWORD rop) /*RASTER OPERATION.*/
/*FILL BAN ON SCREEN WITH SEMITRANSPARENT BRUSH.*/
{
  HDC hdcC;
  HBITMAP hbitC,pbit;
  HBRUSH hbrush,pbrush;
  POINT *points;
  int i,ix,iy,maxX,minX,maxY,minY;

  points=(POINT *)malloc(gb.nnod*sizeof(POINT));
  if(points==NULL) return;

  for(i=0;i<gb.nnod;i++) /*PROJECTION*/
  {
    if(!nodeontoscreen(*(*(gb.nods+i)),&ix,&iy,vp)) return;
    (points+i)->x=ix;
    (points+i)->y=iy;

    if(i==0)
    {
      maxX=(points+i)->x;
      minX=(points+i)->x;
      maxY=(points+i)->y;
      minY=(points+i)->y;
    }
    else
    {
      if(maxX<((points+i)->x)) maxX=(points+i)->x;
      if(minX>((points+i)->x)) minX=(points+i)->x;
      if(maxY<((points+i)->y)) maxY=(points+i)->y;
      if(minY>((points+i)->y)) minY=(points+i)->y;
    }
  }

  hdcC = CreateCompatibleDC(hdc);

  hbitC = CreateCompatibleBitmap(hdc,maxX,maxY);
  pbit=SelectObject(hdcC,hbitC);

  hbrush = (HBRUSH)GetCurrentObject(hdc,OBJ_BRUSH);
  pbrush = SelectObject(hdcC,hbrush);

  PatBlt(hdcC,minX,minY,(maxX-minX),(maxY-minY),PATCOPY);

  hbrush = (HBRUSH)CreateSolidBrush(RGB(r,g,b));
  SelectObject(hdcC,hbrush);

  Polygon(hdcC,points,gb.nnod);

  BitBlt(hdc,minX,minY,(maxX-minX),(maxY-minY),
         hdcC,minX,minY,
         rop); /*GROUND BLACK:SRCPAINT WHITE:SRCAND.*/

  SelectObject(hdcC,pbrush);
  SelectObject(hdcC,pbit);
  DeleteObject(hbrush);
  DeleteObject(hbitC);
  DeleteDC(hdcC);
  free(points);

  return;
}/*fillglobalban*/

void drawglobalban(HDC hdc,struct viewparam vp,struct obans gb,
                   int r,int g,int b,
                   DWORD rop) /*RASTER OPERATION.*/
/*DRAW BAN ON SCREEN BY AXONOMETRICS,PERSPECTIVE.*/
{
  HPEN hpen,ppen;
  int i;
  struct onode gn1,gn2;

  if(rop==SRCAND) hpen=CreatePen(PS_SOLID,1,RGB(0.5*r,0.5*g,0.5*b));
  else            hpen=CreatePen(PS_SOLID,1,RGB(r,g,b));
  ppen=SelectObject(hdc,hpen);
  for(i=0;i<gb.nnod;i++)
  {
    gn1=*(*(gb.nods+i));
    if(i<gb.nnod-1) gn2=*(*(gb.nods+i+1));
    else            gn2=*(*(gb.nods+0));

    drawgloballine(hdc,vp,gn1,gn2);
  }
  SelectObject(hdc,ppen);
  DeleteObject(hpen);

  if(rop==SRCAND)
  {
    fillglobalban(hdc,vp,gb,r,g,b,rop);
  }
  else
  {
    fillglobalban(hdc,vp,gb,(int)(0.5*r),(int)(0.5*g),(int)(0.5*b),
                  rop);
  }

  return;
}/*drawglobalban*/

void drawelement(HDC hdc,struct viewparam vp,
                 struct oelem ge,
                 int mode)
/*DRAW ELEMENT ON SCREEN BY AXONOMETRICS,PERSPECTIVE.*/
{
  HPEN hpen,ppen;
  DWORD rop; /*RASTER OPERATION.*/
  int i;

  if(mode==ONPRINTER) hpen=CreatePen(PS_SOLID,1,RGB(0,0,0));
  else hpen=CreatePen(PS_SOLID,1,RGB(ge.color.line.r,
                                     ge.color.line.g,
                                     ge.color.line.b));
  ppen=SelectObject(hdc,hpen);
  if(ge.nnod==2)
  {
    drawgloballine(hdc,vp,*(*(ge.nods+0)),*(*(ge.nods+1)));
  }
  SelectObject(hdc,ppen);
  DeleteObject(hpen);

  if(mode==ONPRINTER) rop=SRCAND;
  else                rop=SRCPAINT;

  for(i=0;i<ge.nban;i++)
  {
    drawglobalban(hdc,vp,*(ge.bans+i),ge.color.line.r,
                                      ge.color.line.g,
                                      ge.color.line.b,rop);
    if(ge.type==SLAB)
    {
      if(mode==ONPRINTER)
      {
        hpen=CreatePen(PS_SOLID,1,RGB((int)(0.5*ge.color.line.r),
                                      (int)(0.5*ge.color.line.g),
                                      (int)(0.5*ge.color.line.b)));
      }
      else
      {
        hpen=CreatePen(PS_SOLID,1,RGB(ge.color.line.r,
                                      ge.color.line.g,
                                      ge.color.line.b));
      }

      ppen=SelectObject(hdc,hpen);

      slabdivision(hdc,vp,*(ge.bans+i));

      SelectObject(hdc,ppen);
      DeleteObject(hpen);
    }
  }
  return;
}/*drawelement*/

void drawglobalaxis(HDC hdc,
                    struct viewparam vp,
                    int r,int g,int b)
/*DRAW GLOBAL AXIS BY AXONOMETRICS,PERSPECTIVE.*/
{
  HPEN hpen,ppen;
  double arrowlength;
  struct onode hn,on;

  SetTextColor(hdc,RGB(100,100,100));                  /*TEXT GRAY.*/
  hpen=CreatePen(PS_SOLID,1,RGB(r,g,b));             /*ARROW COLOR.*/
  ppen=SelectObject(hdc,hpen);

  arrowlength=vp.dparam.gaxis;

  on=setcoord(0.0,0.0,0.0);
  drawglobaltext(hdc,vp,on,"O");                           /*ORIGIN*/

  hn=setcoord(arrowlength,0.0,0.0);
  drawglobaltext(hdc,vp,hn,"X");                           /*TEXT X*/
  drawglobalarrow(hdc,vp,on,hn,0.0);                      /*ARROW X*/

  hn=setcoord(0.0,arrowlength,0.0);
  drawglobaltext(hdc,vp,hn,"Y");                           /*TEXT Y*/
  drawglobalarrow(hdc,vp,on,hn,0.0);                      /*ARROW Y*/

  hn=setcoord(0.0,0.0,arrowlength);
  drawglobaltext(hdc,vp,hn,"Z");                           /*TEXT Z*/
  drawglobalarrow(hdc,vp,on,hn,0.0);                      /*ARROW Z*/

  SelectObject(hdc,ppen);
  DeleteObject(hpen);

  return;
}/*drawglobalaxis*/

void draworganization(HDC hdc,struct viewparam vp,
                      struct organ go,
                      int mode)
/*DRAW ORGANIZATION ON SCREEN BY AXONOMETRICS,PERSPECTIVE.*/
{
  int i;

  if(vp.vflag.axis==1)
  {
    drawglobalaxis(hdc,vp,0,0,255);             /*DRAW GLOBAL AXIS.*/
  }

  for(i=0;i<go.nelem;i++) /*ELEMENTS.*/
  {
    if((go.elems+i)->sect->dflag==1)
    {
      (go.elems+i)->color.line.r=(go.elems+i)->sect->dcolor.r;
      (go.elems+i)->color.line.g=(go.elems+i)->sect->dcolor.g;
      (go.elems+i)->color.line.b=(go.elems+i)->sect->dcolor.b;

      drawelement(hdc,vp,*(go.elems+i),mode);
    }
  }

  SetTextColor(hdc,RGB(150,150,255));
  for(i=0;i<go.nnode;i++) /*NODES.*/
  {
    drawglobalnode(hdc,vp,*(go.nodes+i));
  }

  return;
}/*draworganization*/

void draworgannodes(HDC hdc,struct viewparam vp,struct organ go)
/*DRAW NODES OF ORGANIZATION ON SCREEN.*/
{
  int i;

  drawglobalaxis(hdc,vp,0,0,255);               /*DRAW GLOBAL AXIS.*/

  SetTextColor(hdc,RGB(0,255,255));
  for(i=0;i<go.nnode;i++)
  {
    drawglobalnode(hdc,vp,*(go.nodes+i));
  }

  return;
}/*draworgannodes*/

struct onode transltog(double **tdrccos,
                       struct onode on, /*ORIGIN GLOBAL.*/
                       struct onode ln) /*LOCAL NODE.*/
/*TRANSFORM POINT LOCAL INTO GLOBAL.*/
{
  struct onode gn; /*TRANSFORMED NODE.*/
  double *ldot,*tdot;

  ldot=(double *)malloc(3*sizeof(double));
  *(ldot+EX)=ln.d[EX];
  *(ldot+EY)=ln.d[EY];
  *(ldot+EZ)=ln.d[EZ];

  tdot=matrixvector(tdrccos,ldot,3);

  gn.d[GX]=*(tdot+GX)+on.d[GX];
  gn.d[GY]=*(tdot+GY)+on.d[GY];
  gn.d[GZ]=*(tdot+GZ)+on.d[GZ];

  free(ldot);
  free(tdot);

  return gn;
}/*transltog*/

struct onode transgtol(double **drccos,
                       struct onode on, /*ORIGIN LOCAL.*/
                       struct onode gn) /*GLOBAL NODE.*/
/*TRANSFORM POINT GLOBAL INTO LOCAL.*/
{
  struct onode ln; /*TRANSFORMED NODE.*/
  double *gdot,*tdot;

  gdot=(double *)malloc(3*sizeof(double));
  *(gdot+GX)=gn.d[GX]-on.d[GX];
  *(gdot+GY)=gn.d[GY]-on.d[GY];
  *(gdot+GZ)=gn.d[GZ]-on.d[GZ];

  tdot=matrixvector(drccos,gdot,3);

  ln.d[EX]=*(tdot+EX);
  ln.d[EY]=*(tdot+EY);
  ln.d[EZ]=*(tdot+EZ);

  free(gdot);
  free(tdot);

  return ln;
}/*transgtol*/

void drawtextonlocaldot(HDC hdc,struct viewparam vp,
                        double **tdrccos,
                        struct onode on, /*ORIGIN*/
                        struct onode ln, /*DOT*/
                        char *str)
/*DRAW GLOBAL TEXT ON LOCAL NODE PROJECTED.*/
{
  struct onode gn;

  gn=transltog(tdrccos,on,ln);
  drawglobaltext(hdc,vp,gn,str);

  return;
}/*drawtextonlocaldot*/

void drawlocalline(HDC hdc,struct viewparam vp,
                   double **tdrccos,
                   struct onode on, /*ORIGIN.*/
                   struct onode ln1, /*DOT HEAD.*/
                   struct onode ln2) /*DOT TAIL.*/
/*DRAW LOCAL LINE PROJECTED.*/
{
  struct onode gn1,gn2;

  gn1=transltog(tdrccos,on,ln1);
  gn2=transltog(tdrccos,on,ln2);
  drawgloballine(hdc,vp,gn1,gn2);

  return;
}/*drawlocalline*/

void drawlocalarrow(HDC hdc,struct viewparam vp,
                    double **tdrccos,
                    struct onode on, /*ORIGIN.*/
                    struct onode ln1, /*DOT HEAD.*/
                    struct onode ln2) /*DOT TAIL.*/
/*DRAW LOCAL LINE WITH GLOBAL ARROW.*/
{
  struct onode gn1,gn2;

  gn1=transltog(tdrccos,on,ln1);
  gn2=transltog(tdrccos,on,ln2);
  drawglobalarrow(hdc,vp,gn1,gn2,0.0);

  return;
}/*drawlocalarrow*/

void drawlocalban(HDC hdc,struct viewparam vp,
                  double **tdrccos,
                  struct onode on, /*ORIGIN.*/
                  struct obans lb) /*LOCAL BAN.*/
/*DRAW LOCAL BAN PROJECTED.*/
{
  /*HPEN hpen,ppen;*/
  int i;
  struct onode ln1,ln2;

  /*hpen=CreatePen(PS_SOLID,1,RGB(r,g,b));*/
  /*ppen=SelectObject(hdc,hpen);*/
  for(i=0;i<lb.nnod;i++)
  {
    ln1=*(*(lb.nods+i));
    if(i<lb.nnod-1) ln2=*(*(lb.nods+i+1));
    else            ln2=*(*(lb.nods+0));

    drawlocalline(hdc,vp,tdrccos,on,ln1,ln2);
  }
  /*SelectObject(hdc,ppen);*/
  /*DeleteObject(hpen);*/

  return;
}/*drawlocalban*/

int insiderange(struct onode n,struct viewrange r)
/*NODE INSIDE RANGE.*/
{
  if(n.d[GX]<r.min.d[GX] || r.max.d[GX]<n.d[GX]) return 0;
  if(n.d[GY]<r.min.d[GY] || r.max.d[GY]<n.d[GY]) return 0;
  if(n.d[GZ]<r.min.d[GZ] || r.max.d[GZ]<n.d[GZ]) return 0;

  return 1;
}/*insiderange*/

int getdirection(struct viewparam vpcurrent,
                 long int dx,long int dy,
                 struct onode focus, /*INITIAL FOCUS.*/
                 double *length)
/*GET DIRECTION OF MOVEMENT BY AXONOMETRICS,PERSPECTIVE.*/
/*{dx,dy}:RELATIVE MOVEMENT OF MOUSE.*/
{
  int Ox,Oy,Xx,Xy,Yx,Yy,Zx,Zy,Dx,Dy;           /*COORDINATION ARROW*/
  double arrowlength;
  double coordfactor[3]={0.0,0.0,0.0},coordlength;
  double distance[3]={0.0,0.0,0.0};
  struct onode pinit;
  struct viewparam vp;

  vp=vpcurrent;

  vp.focus.d[GX]=focus.d[GX]; /*FOCUS CONSTANT.*/
  vp.focus.d[GY]=focus.d[GY];
  vp.focus.d[GZ]=focus.d[GZ];

  arrowlength=1.0;

  pinit=setcoord(0.0,0.0,0.0);                             /*ORIGIN*/
  nodeontoscreen(pinit,&Ox,&Oy,vp);

  pinit=setcoord(arrowlength,0.0,0.0);                     /*AXIS X*/
  nodeontoscreen(pinit,&Xx,&Xy,vp);

  pinit=setcoord(0.0,arrowlength,0.0);                     /*AXIS Y*/
  nodeontoscreen(pinit,&Yx,&Yy,vp);

  pinit=setcoord(0.0,0.0,arrowlength);                     /*AXIS Z*/
  nodeontoscreen(pinit,&Zx,&Zy,vp);

  Dx=Xx-Ox;                                     /*VECTOR OF AXIS X.*/
  Dy=Xy-Oy;
  coordlength=Dx*Dx+Dy*Dy;            /*SQUARE OF LENGTH OF AXIS X.*/
  if(coordlength!=0.0)    /*INNER PRODUCT {AXIS X}{PERPENDICULAR}=0*/
  {
    coordfactor[0]=(dx*Dx+dy*Dy)/coordlength;
    distance[0]=coordfactor[0]*coordfactor[0]*coordlength;
  }
  Dx=Yx-Ox;                                     /*VECTOR OF AXIS Y.*/
  Dy=Yy-Oy;
  coordlength=Dx*Dx+Dy*Dy;            /*SQUARE OF LENGTH OF AXIS Y.*/
  if(coordlength!=0.0)    /*INNER PRODUCT {AXIS Y}{PERPENDICULAR}=0*/
  {
    coordfactor[1]=(dx*Dx+dy*Dy)/coordlength;
    distance[1]=coordfactor[1]*coordfactor[1]*coordlength;
  }
  Dx=Zx-Ox;                                     /*VECTOR OF AXIS Z.*/
  Dy=Zy-Oy;
  coordlength=Dx*Dx+Dy*Dy;            /*SQUARE OF LENGTH OF AXIS Z.*/
  if(coordlength!=0.0)    /*INNER PRODUCT {AXIS Z}{PERPENDICULAR}=0*/
  {
    coordfactor[2]=(dx*Dx+dy*Dy)/coordlength;
    distance[2]=coordfactor[2]*coordfactor[2]*coordlength;
  }

  if(distance[2]>=distance[0] && distance[2]>=distance[1])
  {
    *length=coordfactor[2]*arrowlength;
    return DIRECTIONZ;
  }
  else if(distance[1]>=distance[2] && distance[1]>=distance[0])
  {
    *length=coordfactor[1]*arrowlength;
    return DIRECTIONY;
  }
  else if(distance[0]>=distance[1] && distance[0]>=distance[2])
  {
    *length=coordfactor[0]*arrowlength;
    return DIRECTIONX;
  }
  return 0;
}/*getdirection*/

void drawarclmnodes(HDC hdc,struct viewparam vp,
                    struct arclmframe af,long int code)
/*DRAW NODES OF ARCLM FRAME.*/
/*CODE,CONFINEMENT,DISPLACEMENT,REACTION.*/
{
  char str[20];
  char flabel[6][10]={"Fx","Fy","Fz","Mx","My","Mz"};       /*LOAD.*/
  char clabel[6][10]={"Dx","Dy","Dz","Rx","Ry","Rz"}; /*COMPULSORY.*/
  char dlabel[3][10]={"Dx","Dy","Dz"};              /*DISPLACEMENT.*/
  char rlabel[6][10]={"Fx","Fy","Fz","Mx","My","Mz"};   /*REACTION.*/
  long int i,j,loff;
  long int nreact=0;
  double pitch=vp.dparam.pitch;
  double value;
  struct onode tn; /*TEXT POSITION.*/

  for(i=0;i<af.nnode;i++)
  {
    if(insiderange(*(af.ninit+i),vp.range))
    {
      if(vp.vflag.nv.code) /*CODES.*/
      {
        if(globalstatus==SELECTNODE && (af.nodes+i)->code==code)
        {
          SetTextColor(hdc,RGB(0,255,255));
        }
        else SetTextColor(hdc,RGB(150,150,255));

        drawglobalnode(hdc,vp,*(af.ninit+i));
      }

      loff=6*((af.ninit+i)->loff);
      tn=*(af.ninit+i);

      for(j=5;j>=0;j--)
      {
        if(vp.vflag.nv.loads[j]) /*LOADS.*/
        {
          if((!(af.confs+loff+j)->iconf) &&
             ((af.confs+loff+j)->value!=0.0))
          {
            sprintf(str,"%s %.5f",flabel[j],
                    (af.confs+loff+j)->value);
            tn.d[GZ]+=pitch;
            SetTextColor(hdc,RGB(255,0,255));       /*LOAD MAGENTA.*/
            drawglobaltext(hdc,vp,tn,str);                  /*TEXT.*/
          }
        }
        if(vp.vflag.nv.confs[j]) /*CONFINEMENTS.*/
        {
          if((af.confs+loff+j)->iconf)
          {
            sprintf(str,"%s %.5f",clabel[j],
                    (af.confs+loff+j)->value);
            tn.d[GZ]+=pitch;
            SetTextColor(hdc,RGB(255,100,0));        /*CONF ORANGE.*/
            drawglobaltext(hdc,vp,tn,str);                  /*TEXT.*/
          }
        }
      }

      for(j=2;j>=0;j--) /*DISPLACEMETS.*/
      {
        if(vp.vflag.nv.disps[j] &&
           (!(af.confs+loff+j)->iconf) &&
           af.fdisp!=NULL)
        {
          vread(af.fdisp,(loff+j+1),&value);
          value-=(af.ninit+i)->d[j];

          sprintf(str,"%s %.5f",dlabel[j],value);

          tn.d[GZ]+=pitch;
          SetTextColor(hdc,RGB(255,255,0));  /*DISPLACEMENT YELLOW.*/
          drawglobaltext(hdc,vp,tn,str);                    /*TEXT.*/
        }
      }

      if(af.freact!=NULL)
      {
        tn=*(af.ninit+i);
        for(j=0;j<=5;j++) /*REACTIONS.*/
        {
          if((af.confs+loff+j)->iconf==1)
          {
            fseek(af.freact,nreact*sizeof(double),SEEK_SET);
            fread(&value,sizeof(double),1,af.freact);
            nreact++;

            if(vp.vflag.nv.react[j])
            {
              sprintf(str,"%s %.5f",rlabel[j],value);

              tn.d[GZ]-=pitch;
              SetTextColor(hdc,RGB(255,0,0));       /*REACTION RED.*/
              drawglobaltext(hdc,vp,tn,str);                /*TEXT.*/
            }
          }
        }
      }
    }
    else
    {
      for(j=0;j<6;j++)
      {
        if((af.confs+loff+j)->iconf==1) nreact++;
      }
    }
  }
  return;
}/*drawarclmnodes*/

void drawglobalwire(HDC hdc,struct viewparam vp,
                    struct arclmframe af,
                    struct owire elem,
                    int cred,int cgreen,int cblue, /*CODE*/
                    int ered,int egreen,int eblue, /*LINE*/
                    long int selectcode)
/*DRAW ELEMENT,PLASTIC HINGE BY AXONOMETRICS,PERSPECTIVE.*/
{
  HPEN hpen,ppen;                                  /*HANDLE OF PEN.*/
  HBRUSH hbrush,pbrush;                          /*HANDLE OF BRUSH.*/
  char str[10];
  int ix1,iy1,ix2,iy2,dx1,dy1,dx2,dy2;    /*PROJECTED COORDINATION.*/
  double dx,dy;
  int hingeX,hingeY;                             /*CENTER OF HINGE.*/
  double oblique;                                   /*OBLIQUE SIDE.*/
  struct onode inode1,inode2;                       /*INITIAL NODE.*/
  struct onode dnode1,dnode2;                      /*DEFORMED NODE.*/

  initialnode(af.ninit,af.nnode,(elem.node[0]->code),&inode1);
  if(!insiderange(inode1,vp.range)) return;
  if(af.fdisp==NULL || !vp.vflag.ev.deformation)
  {
    dnode1.d[GX]=0.0;
    dnode1.d[GY]=0.0;
    dnode1.d[GZ]=0.0;
  }
  else                                              /*DISPLACEMENT.*/
  {
    /*inputnode(af.fdisp,elem.node[0]);*/ /*DEFORMED NODE.*/
    dnode1.d[GX]=elem.node[0]->d[GX]-inode1.d[GX];
    dnode1.d[GY]=elem.node[0]->d[GY]-inode1.d[GY];
    dnode1.d[GZ]=elem.node[0]->d[GZ]-inode1.d[GZ];

    dnode1.d[GX]=inode1.d[GX]+(vp.dparam.dfact)*dnode1.d[GX];
    dnode1.d[GY]=inode1.d[GY]+(vp.dparam.dfact)*dnode1.d[GY];
    dnode1.d[GZ]=inode1.d[GZ]+(vp.dparam.dfact)*dnode1.d[GZ];

    nodeontoscreen(dnode1,&dx1,&dy1,vp);              /*PROJECTION.*/
  }
  nodeontoscreen(inode1,&ix1,&iy1,vp);                /*PROJECTION.*/

  initialnode(af.ninit,af.nnode,(elem.node[1]->code),&inode2);
  if(!insiderange(inode2,vp.range)) return;
  if(af.fdisp==NULL || !vp.vflag.ev.deformation)
  {
    dnode2.d[GX]=0.0;
    dnode2.d[GY]=0.0;
    dnode2.d[GZ]=0.0;
  }
  else                                              /*DISPLACEMENT.*/
  {
    /*inputnode(af.fdisp,elem.node[0]);*/ /*DEFORMED NODE.*/
    dnode2.d[GX]=elem.node[1]->d[GX]-inode2.d[GX];
    dnode2.d[GY]=elem.node[1]->d[GY]-inode2.d[GY];
    dnode2.d[GZ]=elem.node[1]->d[GZ]-inode2.d[GZ];

    dnode2.d[GX]=inode2.d[GX]+(vp.dparam.dfact)*dnode2.d[GX];
    dnode2.d[GY]=inode2.d[GY]+(vp.dparam.dfact)*dnode2.d[GY];
    dnode2.d[GZ]=inode2.d[GZ]+(vp.dparam.dfact)*dnode2.d[GZ];

    nodeontoscreen(dnode2,&dx2,&dy2,vp);              /*PROJECTION.*/
  }
  nodeontoscreen(inode2,&ix2,&iy2,vp);                /*PROJECTION.*/

  if(vp.vflag.nv.code)
  {
    if(globalstatus==SELECTNODE &&
       elem.node[0]->code==selectcode)
    {
      SetTextColor(hdc,RGB(0,255,255));
    }
    else SetTextColor(hdc,RGB(150,150,255));

    sprintf(str,"%d",elem.node[0]->code);
    TextOut(hdc,ix1,iy1,str,strlen(str));

    if(globalstatus==SELECTNODE &&
       elem.node[1]->code==selectcode)
    {
      SetTextColor(hdc,RGB(0,255,255));
    }
    else SetTextColor(hdc,RGB(150,150,255));

    sprintf(str,"%d",elem.node[1]->code);
    TextOut(hdc,ix2,iy2,str,strlen(str));
  }
  if(vp.vflag.ev.code)
  {
    SetTextColor(hdc,RGB(cred,cgreen,cblue)); /*ELEM CODE.*/
    sprintf(str,"%d",elem.code);
    TextOut(hdc,(0.55*(ix2-ix1)+ix1),
                (0.55*(iy2-iy1)+iy1),str,strlen(str));
  }
  if(vp.vflag.ev.sectioncode)
  {
    SetTextColor(hdc,RGB(elem.sect->dcolor.r,
                         elem.sect->dcolor.g,
                         elem.sect->dcolor.b)); /*SECT CODE.*/
    sprintf(str,"%d",elem.sect->code);
    TextOut(hdc,(0.45*(ix2-ix1)+ix1),
                (0.45*(iy2-iy1)+iy1),str,strlen(str));
  }

  hpen=CreatePen(PS_SOLID,1,RGB(ered,egreen,eblue)); /*ELEMENT PEN.*/
  ppen=SelectObject(hdc,hpen);
  MoveToEx(hdc,ix1,iy1,NULL);
  LineTo(hdc,ix2,iy2);
  if(vp.vflag.ev.deformation)
  {
    MoveToEx(hdc,dx1,dy1,NULL);
    LineTo(hdc,dx2,dy2);
  }
  SelectObject(hdc,ppen);
  DeleteObject(hpen);

  if(vp.vflag.ev.axis)                              /*ELEMENT AXIS.*/
  {
    drawwireaxis(hdc,vp,inode1,inode2,elem.cangle);
  }

  if(vp.vflag.ev.hinge)                            /*INITIAL HINGE.*/
  {
    hpen=CreatePen(PS_SOLID,1,RGB(255,255,255));
    ppen=SelectObject(hdc,hpen);
    if(elem.iconf[0][3]==1 ||
       elem.iconf[0][4]==1 ||
       elem.iconf[0][5]==1)                        /*HINGE ON HEAD.*/
    {
      dx=ix2-ix1; dy=iy2-iy1;
      oblique=sqrt(dx*dx+dy*dy);
      if(oblique==0.0)
      {
        hingeX=ix1;
        hingeY=iy1;
      }
      else
      {
        hingeX=ix1+(int)(vp.dparam.hsize*dx/oblique);
        hingeY=iy1+(int)(vp.dparam.hsize*dy/oblique);
      }
      Arc(hdc,
          (hingeX-vp.dparam.hsize),(hingeY-vp.dparam.hsize),
          (hingeX+vp.dparam.hsize),(hingeY+vp.dparam.hsize),
          0,0,0,0);                                /*HOLLOW CIRCLE.*/
    }
    if(elem.iconf[1][3]==1 ||
       elem.iconf[1][4]==1 ||
       elem.iconf[1][5]==1)                        /*HINGE ON TAIL.*/
    {
      dx=ix1-ix2; dy=iy1-iy2;
      oblique=sqrt(dx*dx+dy*dy);
      if(oblique==0.0)
      {
        hingeX=ix2;
        hingeY=iy2;
      }
      else
      {
        hingeX=ix2+(int)(vp.dparam.hsize*dx/oblique);
        hingeY=iy2+(int)(vp.dparam.hsize*dy/oblique);
      }
      Arc(hdc,
          (hingeX-vp.dparam.hsize),(hingeY-vp.dparam.hsize),
          (hingeX+vp.dparam.hsize),(hingeY+vp.dparam.hsize),
          0,0,0,0);                                /*HOLLOW CIRCLE.*/
    }

    SelectObject(hdc,ppen);
    DeleteObject(hpen);
  }
  if(vp.vflag.ev.hinge && af.fdisp!=NULL) /*INITIAL HINGE DEFORMED.*/
  {
    hpen=CreatePen(PS_SOLID,1,RGB(255,255,255));
    ppen=SelectObject(hdc,hpen);
    if(elem.iconf[0][3]==1 ||
       elem.iconf[0][4]==1 ||
       elem.iconf[0][5]==1)                        /*HINGE ON HEAD.*/
    {
      dx=dx2-dx1; dy=dy2-dy1;
      oblique=sqrt(dx*dx+dy*dy);
      if(oblique==0.0)
      {
        hingeX=dx1;
        hingeY=dy1;
      }
      else
      {
        hingeX=dx1+(int)(vp.dparam.hsize*dx/oblique);
        hingeY=dy1+(int)(vp.dparam.hsize*dy/oblique);
      }
      Arc(hdc,
          (hingeX-vp.dparam.hsize),(hingeY-vp.dparam.hsize),
          (hingeX+vp.dparam.hsize),(hingeY+vp.dparam.hsize),
          0,0,0,0);                                /*HOLLOW CIRCLE.*/
    }
    if(elem.iconf[1][3]==1 ||
       elem.iconf[1][4]==1 ||
       elem.iconf[1][5]==1)                        /*HINGE ON TAIL.*/
    {
      dx=dx1-dx2; dy=dy1-dy2;
      oblique=sqrt(dx*dx+dy*dy);
      if(oblique==0.0)
      {
        hingeX=dx2;
        hingeY=dy2;
      }
      else
      {
        hingeX=dx2+(int)(vp.dparam.hsize*dx/oblique);
        hingeY=dy2+(int)(vp.dparam.hsize*dy/oblique);
      }
      Arc(hdc,
          (hingeX-vp.dparam.hsize),(hingeY-vp.dparam.hsize),
          (hingeX+vp.dparam.hsize),(hingeY+vp.dparam.hsize),
          0,0,0,0);                                /*HOLLOW CIRCLE.*/
    }

    SelectObject(hdc,ppen);
    DeleteObject(hpen);
  }
  if(vp.vflag.ev.hinge && af.fdisp!=NULL)
  {
    hpen=CreatePen(PS_SOLID,1,RGB(0,255,255));     /*PLASTIC HINGE.*/
    ppen=SelectObject(hdc,hpen);
    hbrush=CreateSolidBrush(RGB(0,255,255));      /*BRUSH OF HINGE.*/
    pbrush = SelectObject(hdc,hbrush);
    if(elem.iconf[0][0]==-1)                /*PLASTIC HINGE ON HEAD*/
    {
      dx=dx2-dx1; dy=dy2-dy1;
      oblique=sqrt(dx*dx+dy*dy);
      if(oblique==0.0)
      {
        hingeX=dx1;
        hingeY=dy1;
      }
      else
      {
        hingeX=dx1+(int)(vp.dparam.hsize*dx/oblique);
        hingeY=dy1+(int)(vp.dparam.hsize*dy/oblique);
      }
      Ellipse(hdc,                                  /*FILLED CIRCLE*/
              (hingeX-vp.dparam.hsize),(hingeY-vp.dparam.hsize),
              (hingeX+vp.dparam.hsize),(hingeY+vp.dparam.hsize));

    }
    if(elem.iconf[1][0]==-1)                /*PLASTIC HINGE ON TAIL*/
    {
      dx=dx1-dx2; dy=dy1-dy2;
      oblique=sqrt(dx*dx+dy*dy);
      if(oblique==0.0)
      {
        hingeX=dx2;
        hingeY=dy2;
      }
      else
      {
        hingeX=dx2+(int)(vp.dparam.hsize*dx/oblique);
        hingeY=dy2+(int)(vp.dparam.hsize*dy/oblique);
      }
      Ellipse(hdc,                                  /*FILLED CIRCLE*/
              (hingeX-vp.dparam.hsize),(hingeY-vp.dparam.hsize),
              (hingeX+vp.dparam.hsize),(hingeY+vp.dparam.hsize));

    }

    SelectObject(hdc,ppen);
    DeleteObject(hpen);
    SelectObject(hdc,pbrush);
    DeleteObject(hbrush);
  }

  return;
}/*drawglobalwire*/

void drawwirestress(HDC hdc,struct viewparam vp,
                    struct arclmframe af,
                    struct owire elem)
/*DRAW STRESS OF ONE ELEMENT BY AXONOMETRICS,PERSPECTIVE.*/
{

  HPEN hpenM,hpenQ,ppen;                            /*HANDLE OF PEN*/
  char str[20];
  char slabel[6][10]={"Nz","Qx","Qy","Mz","Mx","My"};     /*LABELS.*/
  int i;
  int qred=255,qgreen=100,qblue=100; /*SHEAR COLOR*/
  int mred=150,mgreen=  0,mblue=255; /*MOMENT COLOR*/
  double value,Qi,Qj,Mi,Mj;
  double elength;
  double **drccos,**tdrccos;
  double mfact,qfact;
  struct onode ninit[2];
  struct onode hn,tn;
  struct owire einit=elem;

  qfact=vp.dparam.qfact;
  mfact=vp.dparam.mfact;

  initialnode(af.ninit,af.nnode,
              (elem.node[0]->code),&(ninit[0]));
  if(!insiderange(ninit[0],vp.range)) return;
  initialnode(af.ninit,af.nnode,
              (elem.node[1]->code),&(ninit[1]));
  if(!insiderange(ninit[1],vp.range)) return;

  einit.node[0]=&(ninit[0]); /*INITIAL HEAD.*/
  einit.node[1]=&(ninit[1]); /*INITIAL TAIL.*/

  elength=distancedotdot(ninit[0],ninit[1]);

  drccos=directioncosine(einit.node[0]->d[GX],
                         einit.node[0]->d[GY],
                         einit.node[0]->d[GZ],
                         einit.node[1]->d[GX],
                         einit.node[1]->d[GY],
                         einit.node[1]->d[GZ],
                         einit.cangle);                  /*[DRCCOS]*/
  tdrccos=matrixtranspose(drccos,3);                          /*[T]*/

  for(i=0;i<=5;i++) /*0:Nz 1:Qx 2:Qy 3:Mz 4:Mx 5:My*/
  {
    if(i==0) SetTextColor(hdc,RGB(qred,qgreen,qblue));     /*SHEAR.*/
    if(i==3) SetTextColor(hdc,RGB(mred,mgreen,mblue));    /*MOMENT.*/

    if(vp.vflag.ev.stress[i])
    {
      sprintf(str,"%s %.5f",slabel[i],elem.stress[0][i]);   /*HEAD.*/
      tn=setcoord(((0.06*(6-i))*elength),0.0,0.0);
      drawtextonlocaldot(hdc,vp,tdrccos,ninit[0],tn,str);

      sprintf(str,"%s %.5f",slabel[i],elem.stress[1][i]);   /*TAIL.*/
      tn=setcoord(((1-0.06*(i+1))*elength),0.0,0.0);
      drawtextonlocaldot(hdc,vp,tdrccos,ninit[0],tn,str);
    }
  }

  hpenQ=CreatePen(PS_SOLID,1,RGB(qred,qgreen,qblue));  /*SHEAR PEN.*/
  ppen=SelectObject(hdc,hpenQ);

  for(i=0;i<=1;i++) /*HEAD,TAIL.*/
  {
    if(vp.vflag.ev.stress[0] && elem.stress[i][0]!=0.0) /*Nz*/
    {
      value=elem.stress[i][0];

      if(value>0.0) /*COMPRESSION*/
      {
        hn=setcoord(((0.05+0.9*i)*elength),0.0,0.0);
        tn=setcoord(((0.05+0.9*i)*elength+qfact*value),0.0,0.0);
        drawlocalarrow(hdc,vp,tdrccos,ninit[0],hn,tn);
      }
      if(value<0.0) /*TENSION*/
      {
        hn=setcoord(((0.05+0.9*i)*elength-qfact*value),0.0,0.0);
        tn=setcoord(((0.05+0.9*i)*elength),0.0,0.0);
        drawlocalarrow(hdc,vp,tdrccos,ninit[0],hn,tn);
      }
    }
    if(vp.vflag.ev.stress[1] && elem.stress[0][1]!=0.0) /*Qx*/
    {
      value=elem.stress[0][1];

      hn=setcoord(((0.05+0.9*i)*elength),(-qfact*0.5*value),0.0);
      tn=setcoord(((0.05+0.9*i)*elength),(qfact*0.5*value),0.0);
      drawlocalarrow(hdc,vp,tdrccos,ninit[0],hn,tn);
    }
    if(vp.vflag.ev.stress[2] && elem.stress[0][2]!=0.0) /*Qy*/
    {
      value=elem.stress[0][2];

      hn=setcoord(((0.05+0.9*i)*elength),0.0,(-qfact*0.5*value));
      tn=setcoord(((0.05+0.9*i)*elength),0.0,(qfact*0.5*value));
      drawlocalarrow(hdc,vp,tdrccos,ninit[0],hn,tn);
    }
  }

  hpenM=CreatePen(PS_SOLID,1,RGB(mred,mgreen,mblue)); /*MOMENT PEN.*/
  SelectObject(hdc,hpenM);
  if(vp.vflag.ev.stress[4]) /*Mx*/
  {
    Qi=elem.stress[0][2];
    Qj=elem.stress[1][2];
    Mi=elem.stress[0][4];
    Mj=elem.stress[1][4];

    hn=setcoord(0.0,0.0,0.0);
    tn=setcoord(0.0,0.0,(-mfact*Mi));
    drawlocalline(hdc,vp,tdrccos,ninit[0],hn,tn);

    if(Qi!=(-Qj))
    {
      hn.d[EZ]=(Qj*elength-Mi-Mj)/(Qi+Qj);
      hn.d[EX]=0.0;
      hn.d[EY]=mfact*(Qi*Mj-Qj*Mi-Qi*Qj*elength)/(Qi+Qj);

      drawlocalline(hdc,vp,tdrccos,ninit[0],hn,tn);
      tn=hn;
    }

    hn=setcoord(elength,0.0,(mfact*Mj));
    drawlocalline(hdc,vp,tdrccos,ninit[0],hn,tn);

    tn=setcoord(elength,0.0,0.0);
    drawlocalline(hdc,vp,tdrccos,ninit[0],hn,tn);
  }
  if(vp.vflag.ev.stress[5]) /*My*/
  {
    Qi=elem.stress[0][1];
    Qj=elem.stress[1][1];
    Mi=elem.stress[0][5];
    Mj=elem.stress[1][5];

    hn=setcoord(0.0,0.0,0.0);
    tn=setcoord(0.0,(mfact*Mi),0.0);
    drawlocalline(hdc,vp,tdrccos,ninit[0],hn,tn);

    if(Qi!=(-Qj))
    {
      hn.d[EZ]=(Qj*elength+Mi+Mj)/(Qi+Qj);
      hn.d[EX]=mfact*(-Qi*Mj+Qj*Mi-Qi*Qj*elength)/(Qi+Qj);
      hn.d[EY]=0.0;

      drawlocalline(hdc,vp,tdrccos,ninit[0],hn,tn);
      tn=hn;
    }

    hn=setcoord(elength,(-mfact*Mj),0.0);
    drawlocalline(hdc,vp,tdrccos,ninit[0],hn,tn);

    tn=setcoord(elength,0.0,0.0);
    drawlocalline(hdc,vp,tdrccos,ninit[0],hn,tn);
  }
  for(i=0;i<=2;i++) free(*(drccos+i));
  free(drccos);
  for(i=0;i<=2;i++) free(*(tdrccos+i));
  free(tdrccos);

  SelectObject(hdc,ppen);
  DeleteObject(hpenQ);
  DeleteObject(hpenM);

  return;
}/*drawwirestress*/

void drawwireaxis(HDC hdc,
                  struct viewparam vp,
                  struct onode n1,struct onode n2,double cangle)
/*DRAW AXIS OF WIRE ELEMENT.*/
{
  HPEN hpen,ppen;                                   /*HANDLE OF PEN*/
  int i;
  double elength;
  double **drccos,**tdrccos;
  double arrowlength=vp.dparam.eaxis;
  struct onode on,hn,tn; /*ORIGIN,HEAD,TAIL.*/

  elength=distancedotdot(n1,n2);

  drccos=directioncosine(n1.d[GX],n1.d[GY],n1.d[GZ],
                         n2.d[GX],n2.d[GY],n2.d[GZ],
                         cangle);                        /*[DRCCOS]*/
  tdrccos=matrixtranspose(drccos,3);                          /*[T]*/

  SetTextColor(hdc,RGB(255,0,255));                 /*AXIS MAGENTA.*/
  hpen=CreatePen(PS_SOLID,1,RGB(255,0,255));   /*AXIS LINE MAGENTA.*/
  ppen=SelectObject(hdc,hpen);

  on=n1;                                            /*LOCAL ORIGIN.*/

  hn=setcoord((0.5*elength-0.5*arrowlength),0.0,0.0);
  tn=setcoord((0.5*elength+0.5*arrowlength),0.0,0.0);
  drawtextonlocaldot(hdc,vp,tdrccos,on,tn,"z");            /*TEXT z*/
  drawlocalarrow(hdc,vp,tdrccos,on,hn,tn);                /*ARROW z*/

  hn=setcoord((0.5*elength-0.5*arrowlength),0.0,0.0);
  tn=setcoord((0.5*elength-0.5*arrowlength),arrowlength,0.0);
  drawtextonlocaldot(hdc,vp,tdrccos,on,tn,"x");            /*TEXT x*/
  drawlocalarrow(hdc,vp,tdrccos,on,hn,tn);                /*ARROW x*/

  hn=setcoord((0.5*elength-0.5*arrowlength),0.0,0.0);
  tn=setcoord((0.5*elength-0.5*arrowlength),0.0,arrowlength);
  drawtextonlocaldot(hdc,vp,tdrccos,on,tn,"y");            /*TEXT y*/
  drawlocalarrow(hdc,vp,tdrccos,on,hn,tn);                /*ARROW y*/

  for(i=0;i<=2;i++) free(*(drccos+i));
  free(drccos);
  for(i=0;i<=2;i++) free(*(tdrccos+i));
  free(tdrccos);

  SelectObject(hdc,ppen);
  DeleteObject(hpen);

  return;
}/*drawwireaxis*/

void drawarclmframe(HDC hdc,struct viewparam vp,
                    struct arclmframe af,
                    long int code,
                    int mode)
{
  int i;
  struct owire elem;

  if(vp.vflag.axis==1)
  {
    drawglobalaxis(hdc,vp,0,0,255);             /*DRAW GLOBAL AXIS.*/
  }

  drawarclmnodes(hdc,vp,af,code);

  for(i=0;i<af.nelem;i++)
  {
    if((af.elems+i)->sect->dflag==1)
    {
      if(globalstatus==SELECTELEMENT && (af.elems+i)->code==code)
      {
        drawglobalwire(hdc,vp,af,*(af.elems+i),255,0,255,
                                               255,0,255,code);
      }
      else if(globalstatus==SELECTSECTION &&
              (af.elems+i)->sect->code==code)
      {
        drawglobalwire(hdc,vp,af,*(af.elems+i),255,0,255,
                                               255,0,255,code);
      }
      else if(mode==ONPRINTER)
      {
        drawglobalwire(hdc,vp,af,*(af.elems+i),0,0,0,
                                               0,0,0,code);
      }
      else
      {
        drawglobalwire(hdc,vp,af,*(af.elems+i),255,255,255,
                       (af.elems+i)->sect->dcolor.r,
                       (af.elems+i)->sect->dcolor.g,
                       (af.elems+i)->sect->dcolor.b,code);
      }

      if(af.felem!=NULL)
      {
        inputelem(af.elems,af.felem,i,&elem);
        drawwirestress(hdc,vp,af,elem);
      }
    }
  }

  return;
}/*drawarclmframe*/

void drawrotatingarclm(HDC hdc,struct viewparam vp,
                       struct arclmframe af)
{
  int i;

  drawglobalaxis(hdc,vp,0,0,255);               /*DRAW GLOBAL AXIS.*/

  SetTextColor(hdc,RGB(0,255,255));
  for(i=0;i<af.nnode;i++)
  {
    drawglobalnode(hdc,vp,*(af.ninit+i));
  }

  return;
}/*drawrotatingarclm*/

void drawyieldsurface(HDC hdc,struct viewparam vp,
                      int f1,int f2,int f3,FILE *fsfc)
/*DRAW YIELD SURFACE BY AXONOMETRICS,PERSPECTIVE.*/
{
  HPEN hpen1,hpen2,ppen;
  char str[10];
  char slabel[6][10]={"Nz","Qx","Qy","Mz","Mx","My"};     /*LABELS.*/
  int n;
  double alength=vp.dparam.gaxis;
  struct line sline;
  struct onode hn,tn;

  SetTextColor(hdc,RGB(100,100,100));            /*AXIS LABEL GRAY.*/
  hpen1=CreatePen(PS_SOLID,1,RGB(0,0,255));       /*AXIS LINE BLUE.*/
  ppen=SelectObject(hdc,hpen1);

  tn=setcoord(0.0,0.0,0.0);
  drawglobaltext(hdc,vp,tn,"O"); /*ORIGIN.*/

  sprintf(str,"%s",slabel[f1]);  /*AXIS 1.*/
  hn=setcoord(-alength,0.0,0.0);
  tn=setcoord(alength,0.0,0.0);
  drawglobaltext(hdc,vp,tn,str);
  drawglobalarrow(hdc,vp,hn,tn,0.0);

  sprintf(str,"%s",slabel[f2]);  /*AXIS 2.*/
  hn=setcoord(0.0,-alength,0.0);
  tn=setcoord(0.0,alength,0.0);
  drawglobaltext(hdc,vp,tn,str);
  drawglobalarrow(hdc,vp,hn,tn,0.0);

  sprintf(str,"%s",slabel[f3]);  /*AXIS 3.*/
  hn=setcoord(0.0,0.0,-alength);
  tn=setcoord(0.0,0.0,alength);
  drawglobaltext(hdc,vp,tn,str);
  drawglobalarrow(hdc,vp,hn,tn,0.0);

  tn=setcoord(1.0,0.0,0.0);
  drawglobaltext(hdc,vp,tn,"1.0");
  tn=setcoord(0.0,1.0,0.0);
  drawglobaltext(hdc,vp,tn,"1.0");
  tn=setcoord(0.0,0.0,1.0);
  drawglobaltext(hdc,vp,tn,"1.0");

  hpen2=CreatePen(PS_SOLID,1,RGB(100,100,100));     /*SURFACE GRAY.*/
  SelectObject(hdc,hpen2);

  sline=setlineends(1.0,0.0,0.0,0.0, 1.0,0.0);
  drawgloballine(hdc,vp,sline.ends[0],sline.ends[1]);
  sline=setlineends(1.0,0.0,0.0,0.0,-1.0,0.0);
  drawgloballine(hdc,vp,sline.ends[0],sline.ends[1]);
  sline=setlineends(1.0,0.0,0.0,0.0, 0.0, 1.0);
  drawgloballine(hdc,vp,sline.ends[0],sline.ends[1]);
  sline=setlineends(1.0,0.0,0.0,0.0, 0.0,-1.0);
  drawgloballine(hdc,vp,sline.ends[0],sline.ends[1]);
  sline=setlineends(-1.0,0.0,0.0,0.0, 1.0,0.0);
  drawgloballine(hdc,vp,sline.ends[0],sline.ends[1]);
  sline=setlineends(-1.0,0.0,0.0,0.0,-1.0,0.0);
  drawgloballine(hdc,vp,sline.ends[0],sline.ends[1]);
  sline=setlineends(-1.0,0.0,0.0,0.0, 0.0, 1.0);
  drawgloballine(hdc,vp,sline.ends[0],sline.ends[1]);
  sline=setlineends(-1.0,0.0,0.0,0.0, 0.0,-1.0);
  drawgloballine(hdc,vp,sline.ends[0],sline.ends[1]);
  sline=setlineends(0.0, 1.0,0.0,0.0, 0.0, 1.0);
  drawgloballine(hdc,vp,sline.ends[0],sline.ends[1]);
  sline=setlineends(0.0, 1.0,0.0,0.0, 0.0,-1.0);
  drawgloballine(hdc,vp,sline.ends[0],sline.ends[1]);
  sline=setlineends(0.0,-1.0,0.0,0.0, 0.0, 1.0);
  drawgloballine(hdc,vp,sline.ends[0],sline.ends[1]);
  sline=setlineends(0.0,-1.0,0.0,0.0, 0.0,-1.0);
  drawgloballine(hdc,vp,sline.ends[0],sline.ends[1]);

  SelectObject(hdc,ppen);
  DeleteObject(hpen1);
  DeleteObject(hpen2);

  if(fsfc!=NULL)
  {
    fseek(fsfc,0L,SEEK_SET);
    while(1)
    {
      n=fread(&sline,sizeof(struct line),1,fsfc);
      if(n==0) break;

      drawlinestruct(hdc,vp,sline);
    }
  }

  return;
}/*drawyieldsurface*/

struct osect *getsection(HWND hdwnd,struct osect *sects,int nsect,
                         long int code)
/*GET SECTION DATA FROM SECTIONS.*/
{
  int i;
  struct osect sect;

  for(i=0;i<nsect;i++)
  {
    if((sects+i)->code==code)
    {
      sect=*(sects+i);                                   /*SECTION.*/
      break;
    }
  }
  if(i>=nsect)
  {
    setsection(hdwnd,NULL);
    return NULL;
  }
  setsection(hdwnd,&sect);

  return sects+i;
}/*getsection*/

void setsection(HWND hdwnd,struct osect *sect)
/*SET SECTION DATA INTO DIALOG.*/
{
  char str[80];

  if(sect==NULL)
  {
    SetDlgItemText(hdwnd,IDS_CODE,    "0");                 /*CODE.*/
    SetDlgItemText(hdwnd,IDS_E,      "\0");      /*YOUNG'S MODULUS.*/
    SetDlgItemText(hdwnd,IDS_POISSON,"\0");      /*POISSON'S RATIO.*/
    SetDlgItemText(hdwnd,IDS_A,      "\0");                 /*AREA.*/
    SetDlgItemText(hdwnd,IDS_IXX,    "\0");
    SetDlgItemText(hdwnd,IDS_IYY,    "\0");
    SetDlgItemText(hdwnd,IDS_J,      "\0");  /*ST.VENANT'S TORTION.*/
    SetDlgItemText(hdwnd,IDS_NUPPER, "\0");     /*UPPER LIMIT OF N.*/
    SetDlgItemText(hdwnd,IDS_NLOWER, "\0");     /*LOWER LIMIT OF N.*/
    SetDlgItemText(hdwnd,IDS_QXUPPER,"\0");    /*UPPER LIMIT OF Qx.*/
    SetDlgItemText(hdwnd,IDS_QXLOWER,"\0");    /*LOWER LIMIT OF Qx.*/
    SetDlgItemText(hdwnd,IDS_QYUPPER,"\0");    /*UPPER LIMIT OF Qy.*/
    SetDlgItemText(hdwnd,IDS_QYLOWER,"\0");    /*LOWER LIMIT OF Qy.*/
    SetDlgItemText(hdwnd,IDS_MXUPPER,"\0");    /*UPPER LIMIT OF Mx.*/
    SetDlgItemText(hdwnd,IDS_MXLOWER,"\0");    /*LOWER LIMIT OF Mx.*/
    SetDlgItemText(hdwnd,IDS_MYUPPER,"\0");    /*UPPER LIMIT OF My.*/
    SetDlgItemText(hdwnd,IDS_MYLOWER,"\0");    /*LOWER LIMIT OF My.*/
    SetDlgItemText(hdwnd,IDS_MZUPPER,"\0");    /*UPPER LIMIT OF Mz.*/
    SetDlgItemText(hdwnd,IDS_MZLOWER,"\0");    /*LOWER LIMIT OF Mz.*/
    return;
  }
  sprintf(str,"%ld",sect->code);
  SetDlgItemText(hdwnd,IDS_CODE,   str);                    /*CODE.*/
  sprintf(str,"%.5f",sect->E);
  SetDlgItemText(hdwnd,IDS_E,      str);         /*YOUNG'S MODULUS.*/
  sprintf(str,"%.5f",sect->poi);
  SetDlgItemText(hdwnd,IDS_POISSON,str);         /*POISSON'S RATIO.*/
  sprintf(str,"%.5f",sect->area);
  SetDlgItemText(hdwnd,IDS_A,      str);                    /*AREA.*/
  sprintf(str,"%.5f",sect->Ixx);
  SetDlgItemText(hdwnd,IDS_IXX,    str);
  sprintf(str,"%.5f",sect->Iyy);
  SetDlgItemText(hdwnd,IDS_IYY,    str);
  sprintf(str,"%.5f",sect->Jzz);
  SetDlgItemText(hdwnd,IDS_J,      str);       /*ST.VENANT TORTION.*/
  sprintf(str,"%.5f",sect->fmax[0]);
  SetDlgItemText(hdwnd,IDS_NUPPER, str);        /*UPPER LIMIT OF N.*/
  sprintf(str,"%.5f",sect->fmin[0]);
  SetDlgItemText(hdwnd,IDS_NLOWER, str);        /*LOWER LIMIT OF N.*/
  sprintf(str,"%.5f",sect->fmax[1]);
  SetDlgItemText(hdwnd,IDS_QXUPPER,str);       /*UPPER LIMIT OF Qx.*/
  sprintf(str,"%.5f",sect->fmin[1]);
  SetDlgItemText(hdwnd,IDS_QXLOWER,str);       /*LOWER LIMIT OF Qx.*/
  sprintf(str,"%.5f",sect->fmax[2]);
  SetDlgItemText(hdwnd,IDS_QYUPPER,str);       /*UPPER LIMIT OF Qy.*/
  sprintf(str,"%.5f",sect->fmin[2]);
  SetDlgItemText(hdwnd,IDS_QYLOWER,str);       /*LOWER LIMIT OF Qy.*/
  sprintf(str,"%.5f",sect->fmax[4]);
  SetDlgItemText(hdwnd,IDS_MXUPPER,str);       /*UPPER LIMIT OF Mx.*/
  sprintf(str,"%.5f",sect->fmin[4]);
  SetDlgItemText(hdwnd,IDS_MXLOWER,str);       /*LOWER LIMIT OF Mx.*/
  sprintf(str,"%.5f",sect->fmax[5]);
  SetDlgItemText(hdwnd,IDS_MYUPPER,str);       /*UPPER LIMIT OF My.*/
  sprintf(str,"%.5f",sect->fmin[5]);
  SetDlgItemText(hdwnd,IDS_MYLOWER,str);       /*LOWER LIMIT OF My.*/
  sprintf(str,"%.5f",sect->fmax[3]);
  SetDlgItemText(hdwnd,IDS_MZUPPER,str);       /*UPPER LIMIT OF Mz.*/
  sprintf(str,"%.5f",sect->fmin[3]);
  SetDlgItemText(hdwnd,IDS_MZLOWER,str);       /*LOWER LIMIT OF Mz.*/

  return;
}/*setsection*/

struct onode *getnode(HWND hdwnd,struct arclmframe *af,
                      long int code)
/*GET NODE DATA FROM NODEPOINTS.*/
{
  int i;
  struct onode node;

  for(i=0;i<(af->nnode);i++)
  {
    if((af->nodes+i)->code==code)
    {
      node=*(af->ninit+i);                                  /*NODE.*/
      break;
    }
  }

  if(i>=af->nnode)
  {
    setnode(hdwnd,NULL,af->confs);
    return NULL;
  };
  setnode(hdwnd,&node,af->confs);

  return (af->ninit+i);
}/*getnode*/

void setnode(HWND hdwnd,struct onode *node,
                        struct oconf *confs)
/*SET NODE DATA INTO DIALOG.*/
{
  char str[80];
  long int loff;

  if(node==NULL)
  {
    SetDlgItemText(hdwnd,IDN_CODE,"0");

    SetDlgItemText(hdwnd,IDN_X,"\0");
    SetDlgItemText(hdwnd,IDN_Y,"\0");
    SetDlgItemText(hdwnd,IDN_Z,"\0");

    SetDlgItemText(hdwnd,IDN_CONFX,  "\0");
    SetDlgItemText(hdwnd,IDN_VALUEX, "\0");
    SetDlgItemText(hdwnd,IDN_CONFY,  "\0");
    SetDlgItemText(hdwnd,IDN_VALUEY, "\0");
    SetDlgItemText(hdwnd,IDN_CONFZ,  "\0");
    SetDlgItemText(hdwnd,IDN_VALUEZ, "\0");
    SetDlgItemText(hdwnd,IDN_CONFTX, "\0");
    SetDlgItemText(hdwnd,IDN_VALUETX,"\0");
    SetDlgItemText(hdwnd,IDN_CONFTY, "\0");
    SetDlgItemText(hdwnd,IDN_VALUETY,"\0");
    SetDlgItemText(hdwnd,IDN_CONFTZ, "\0");
    SetDlgItemText(hdwnd,IDN_VALUETZ,"\0");
    return;
  };
  sprintf(str,"%ld",node->code);
  SetDlgItemText(hdwnd,IDN_CODE,str);
  sprintf(str,"%.5f",node->d[GX]);
  SetDlgItemText(hdwnd,IDN_X,str);
  sprintf(str,"%.5f",node->d[GY]);
  SetDlgItemText(hdwnd,IDN_Y,str);
  sprintf(str,"%.5f",node->d[GZ]);
  SetDlgItemText(hdwnd,IDN_Z,str);

  loff=6*(node->loff); /*NODE CONFINEMENTS.*/
  setnodeconf(hdwnd,IDN_CONFX ,IDN_VALUEX ,(confs+loff+0));
  setnodeconf(hdwnd,IDN_CONFY ,IDN_VALUEY ,(confs+loff+1));
  setnodeconf(hdwnd,IDN_CONFZ ,IDN_VALUEZ ,(confs+loff+2));
  setnodeconf(hdwnd,IDN_CONFTX,IDN_VALUETX,(confs+loff+3));
  setnodeconf(hdwnd,IDN_CONFTY,IDN_VALUETY,(confs+loff+4));
  setnodeconf(hdwnd,IDN_CONFTZ,IDN_VALUETZ,(confs+loff+5));

  return;
}/*setnode*/

void setnodeconf(HWND hdwnd,int idiconf,int idvconf,
                 struct oconf *conf)
/*SET CONFINEMENT INTO DIALOG.*/
{
  char str[80];

  if(conf->iconf==0) SetDlgItemText(hdwnd,idiconf, "FREE");
  if(conf->iconf==1) SetDlgItemText(hdwnd,idiconf, "FIXED");
  sprintf(str,"%.5f",conf->value);
  SetDlgItemText(hdwnd,idvconf,str);

  return;
}/*setnodeconf*/

struct onode *selectnode(struct viewparam vp,
                         struct arclmframe *af,POINT point)
/*SELECT NODE BY CURSOR.*/
/*RETURN:CODE OF NODE.*/
{
  int i;
  int ix1,iy1;                                /*WINDOW COORDINATION*/
  struct onode node;
  int boxsize=5;                                   /*CURSOR REGION.*/

  for(i=0;i<(af->nnode);i++)
  {
    node=*(af->ninit+i);
    if(insiderange(node,vp.range))
    {
      nodeontoscreen(node,&ix1,&iy1,vp);               /*PROJECTION*/

      if((point.x-boxsize)<=ix1 && ix1<=(point.x+boxsize))
      {
        if((point.y-boxsize)<=iy1 && iy1<=(point.y+boxsize))
        {
          return (af->ninit+i);
        }
      }
    }
  }
  return NULL;
}/*selectnode*/

struct owire *getelement(HWND hdwnd,struct arclmframe *af,
                         long int code)
/*GET ELEMENT DATA FROM ELEMENTS.*/
{
  int i;
  struct owire elem;

  for(i=0;i<(af->nelem);i++)
  {
    if((af->elems+i)->code==code)
    {
      elem=*(af->elems+i);                               /*ELEMENT.*/
      break;
    }
  }

  if(i>=(af->nelem))
  {
    setelement(hdwnd,NULL);
    return NULL;
  };
  setelement(hdwnd,&elem);
  setnode(hdwnd,elem.node[0],af->confs);
  setsection(hdwnd,elem.sect);

  return (af->elems+i);
}/*getelement*/

void setelement(HWND hdwnd,struct owire *elem)
/*SET ELEMENT DATA INTO DIALOG.*/
{
  char str[80];

  if(elem==NULL)
  {
    SetDlgItemText(hdwnd,IDE_CODE,       "0");

    SetDlgItemText(hdwnd,IDE_NHEAD,     "\0");
    SetDlgItemText(hdwnd,IDE_NTAIL,     "\0");
    SetDlgItemText(hdwnd,IDE_SECTION,   "\0");
    SetDlgItemText(hdwnd,IDE_COORDANGLE,"\0");
    SetDlgItemText(hdwnd,IDE_BHEADX,    "\0");
    SetDlgItemText(hdwnd,IDE_BHEADY,    "\0");
    SetDlgItemText(hdwnd,IDE_BHEADZ,    "\0");
    SetDlgItemText(hdwnd,IDE_BHEADTX,   "\0");
    SetDlgItemText(hdwnd,IDE_BHEADTY,   "\0");
    SetDlgItemText(hdwnd,IDE_BHEADTZ,   "\0");
    SetDlgItemText(hdwnd,IDE_BTAILX,    "\0");
    SetDlgItemText(hdwnd,IDE_BTAILY,    "\0");
    SetDlgItemText(hdwnd,IDE_BTAILZ,    "\0");
    SetDlgItemText(hdwnd,IDE_BTAILTX,   "\0");
    SetDlgItemText(hdwnd,IDE_BTAILTY,   "\0");
    SetDlgItemText(hdwnd,IDE_BTAILTZ,   "\0");
    return;
  };
  sprintf(str,"%ld",elem->code);
  SetDlgItemText(hdwnd,IDE_CODE,     str);
  sprintf(str,"%ld",elem->node[0]->code);
  SetDlgItemText(hdwnd,IDE_NHEAD,     str);
  sprintf(str,"%ld",elem->node[1]->code);
  SetDlgItemText(hdwnd,IDE_NTAIL,     str);
  sprintf(str,"%ld",elem->sect->code);
  SetDlgItemText(hdwnd,IDE_SECTION,   str);
  sprintf(str,"%.5f",elem->cangle);
  SetDlgItemText(hdwnd,IDE_COORDANGLE,str);

  setelemconf(hdwnd,elem->iconf[0][4],IDE_BHEADTX);
  setelemconf(hdwnd,elem->iconf[0][5],IDE_BHEADTY);
  setelemconf(hdwnd,elem->iconf[0][3],IDE_BHEADTZ);
  setelemconf(hdwnd,elem->iconf[1][4],IDE_BTAILTX);
  setelemconf(hdwnd,elem->iconf[1][5],IDE_BTAILTY);
  setelemconf(hdwnd,elem->iconf[1][3],IDE_BTAILTZ);

  return;
}/*setelement*/

void setelemconf(HWND hdwnd,signed char iconf,int idiconf)
/*SET CONFINEMENT INTO DIALOG.*/
{
  if(iconf==0) SetDlgItemText(hdwnd,idiconf,"RIGID");
  if(iconf==1) SetDlgItemText(hdwnd,idiconf,"HINGE");

  return;
}/*setelemconf*/

struct owire *selectelement(struct viewparam vp,
                            struct arclmframe *af,POINT point)
/*SELECT ELEMENT BY CURSOR.*/
/*RETURN:CODE OF ELEMENT.*/
{
  int i,n1,n2,n3,n4;
  int ix1,iy1,ix2,iy2;                       /*WINDOW COORDINATION.*/
  struct owire elem;
  int boxsize=5;                                   /*CURSOR REGION.*/

  for(i=0;i<(af->nelem);i++)
  {
    elem=*(af->elems+i);

    if(insiderange(*(elem.node[0]),vp.range))
    {
      nodeontoscreen(*(elem.node[0]),&ix1,&iy1,vp);    /*PROJECTION*/

      if(insiderange(*(elem.node[1]),vp.range))
      {
        nodeontoscreen(*(elem.node[1]),&ix2,&iy2,vp);  /*PROJECTION*/

        n1=intersectlineline((point.x-boxsize),(point.y-boxsize),
                             (point.x+boxsize),(point.y-boxsize),
                             ix1,iy1,ix2,iy2);
        n2=intersectlineline((point.x+boxsize),(point.y-boxsize),
                             (point.x+boxsize),(point.y+boxsize),
                             ix1,iy1,ix2,iy2);
        n3=intersectlineline((point.x+boxsize),(point.y+boxsize),
                             (point.x-boxsize),(point.y+boxsize),
                             ix1,iy1,ix2,iy2);
        n4=intersectlineline((point.x-boxsize),(point.y+boxsize),
                             (point.x-boxsize),(point.y-boxsize),
                             ix1,iy1,ix2,iy2);

        if(n1 || n2 || n3 || n4) return (af->elems+i);
      }
    }
  }
  return NULL;
}/*selectelement*/

void initializeorganization(struct organ *org)
/*INITIALIZE ORGAN.*/
{
  org->code=0;
  org->loff=0;

  org->nprop=0;
  org->nsect=0;
  org->nnode=0;
  org->nelem=0;

  return;
}/*initializeorganization*/

int inputorganization(FILE *fin,struct organ *org)
/*INPUT ORGAN FROM TEXTFILE.*/
{
  /*APPELATION*/
  /*NNODE,NELEM,NPROP,NSECT*/
  /*PROP:PNAME,HIJU,E,POI*/
  /*SECT:SNAME*/
  /*     NFIG*/
  /*     FIG:FPROP,AREA,IXX,IYY,VEN*/
  /*     THICK,EXP,NZMAX,...,MYMIN*/
  /*NODE:CORD,ICON,VCON*/
  /*ELEM:ESECT,CANG,TYPE,ROLE*/
  /*     ENODS*/
  /*     ENOD:BONDS*/
  /*     EBANS*/
  /*     EBAN:BNODS,BNOD:*/

  char **data,str[256];
  int i,j,nstr,pstr,ident;
  long int code;
  long int pnode=(-1),pelem=(-1),pprop=(-1),psect=(-1);
  long int psfig,peban;

  fseek(fin,0L,SEEK_SET);

  org->code=0;
  org->loff=0;

  org->nnode=0; /*INITIALIZATION.*/
  org->nelem=0;
  org->nprop=0;
  org->nsect=0;

  fgets(str,256,fin); /*APPELATION.*/

  while(1)
  {
    data=fgetsbrk(fin,&nstr);
    if(nstr==0) return 0;

    pstr=0; /*POSITION IN "DATA".*/

    while((nstr-pstr)>0)
    {
      if(nstr-pstr==1) /*POINTING LAST STRING.*/
      {
        pstr++;
      }
      else
      {
        ident=0;

        if(nstr-pstr>=2 && !strcmp(*(data+pstr),"NNODE"))
        {
          pstr++;
          org->nnode=(int)strtol(*(data+pstr),NULL,10);
          org->nodes=(struct onode *)malloc(org->nnode
                                            *sizeof(struct onode));
          org->confs=(struct oconf *)malloc(6*(org->nnode)
                                            *sizeof(struct oconf));
          pstr++;
          ident++;
        }
        if(nstr-pstr>=2 && !strcmp(*(data+pstr),"NELEM"))
        {
          pstr++;
          org->nelem=strtol(*(data+pstr),NULL,10);
          org->elems=(struct oelem *)malloc(org->nelem
                                            *sizeof(struct oelem));
          pstr++;

          for(i=0;i<(org->nelem);i++)
          {
            (org->elems+i)->role=0;
            (org->elems+i)->type=0;
            (org->elems+i)->nnod=0;
            (org->elems+i)->nban=0;
          }
          ident++;
        }
        if(nstr-pstr>=2 && !strcmp(*(data+pstr),"NPROP"))
        {
          pstr++;
          org->nprop=strtol(*(data+pstr),NULL,10);
          org->props=(struct oprop *)malloc(org->nprop
                                            *sizeof(struct oprop));
          pstr++;
          ident++;
        }
        if(nstr-pstr>=2 && !strcmp(*(data+pstr),"NSECT"))
        {
          pstr++;
          org->nsect=strtol(*(data+pstr),NULL,10);
          org->sects=(struct osect *)malloc(org->nsect
                                            *sizeof(struct osect));
          pstr++;
          ident++;
        }

        /*PROPERTIES:CODE,HIJU,E,POI.*/
        if(nstr-pstr>=2 && (org->nprop)-pprop>1 &&
           !strcmp(*(data+pstr),"PROP"))
        {
          pstr++;
          pprop++;
          (org->props+pprop)->code=strtol(*(data+pstr),NULL,10);
          (org->props+pprop)->loff=pprop;
          (org->props+pprop)->name=NULL;
          pstr++;
          ident++;
        }
        if(nstr-pstr>=2 && (org->nprop)-pprop>0 &&
           !strcmp(*(data+pstr),"PNAME"))
        {
          pstr++;
          (org->props+pprop)->name=(char *)
                                   malloc((strlen(*(data+pstr))+1)
                                          *sizeof(char));
          strcpy((org->props+pprop)->name,*(data+pstr));
          pstr++;
          ident++;
        }
        if(nstr-pstr>=2 && (org->nprop)-pprop>0 &&
           !strcmp(*(data+pstr),"HIJU"))
        {
          pstr++;
          (org->props+pprop)->hiju=strtod(*(data+pstr),NULL);
          pstr++;
          ident++;
        }
        if(nstr-pstr>=2 && (org->nprop)-pprop>0 &&
           !strcmp(*(data+pstr),"E"))
        {
          pstr++;
          (org->props+pprop)->E=strtod(*(data+pstr),NULL);
          pstr++;
          ident++;
        }
        if(nstr-pstr>=2 && (org->nprop)-pprop>0 &&
           !strcmp(*(data+pstr),"POI"))
        {
          pstr++;
          (org->props+pprop)->poi=strtod(*(data+pstr),NULL);
          pstr++;
          ident++;
        }

        /*SECTIONS:CODE,FIGS,SURFACE.*/
        if(nstr-pstr>=2 && (org->nsect)-psect>1 &&
           !strcmp(*(data+pstr),"SECT"))
        {
          pstr++;
          psect++;
          (org->sects+psect)->code=strtol(*(data+pstr),NULL,10);
          (org->sects+psect)->loff=psect;

          (org->sects+psect)->name=NULL;
          (org->sects+psect)->nfig=0;
          (org->sects+psect)->dflag=1;
          (org->sects+psect)->dcolor.r=255;
          (org->sects+psect)->dcolor.g=255;
          (org->sects+psect)->dcolor.b=255;
          pstr++;
          ident++;
        }
        if(nstr-pstr>=2 && (org->nsect)-psect>0 &&
           !strcmp(*(data+pstr),"SNAME"))
        {
          pstr++;
          (org->sects+psect)->name=(char *)
                                   malloc((strlen(*(data+pstr))+1)
                                          *sizeof(char));
          strcpy((org->sects+psect)->name,*(data+pstr));
          pstr++;
          ident++;
        }
        if(nstr-pstr>=2 && (org->nsect)-psect>0 &&
           !strcmp(*(data+pstr),"NFIG"))
        {
          pstr++;
          psfig=(-1); /*POINTING "FIGS".*/
          (org->sects+psect)->nfig=strtol(*(data+pstr),NULL,10);
          (org->sects+psect)->figs=(struct ofigs *)
                                   malloc((org->sects+psect)->nfig
                                          *sizeof(struct ofigs));
          pstr++;
          ident++;
        }
        if(nstr-pstr>=2 &&
           (org->nsect)-psect>0 &&
           !strcmp(*(data+pstr),"FIG") &&
           ((org->sects+psect)->nfig)-psfig>1)
        {
          pstr++;
          psfig++;
          ((org->sects+psect)->figs+psfig)->code
          =strtol(*(data+pstr),NULL,10);
          ((org->sects+psect)->figs+psfig)->loff=psfig;
          pstr++;
          ident++;
        }
        if(nstr-pstr>=2 &&
           (org->nsect)-psect>0 &&
           !strcmp(*(data+pstr),"FPROP") &&
           ((org->sects+psect)->nfig)-psfig>0)
        {
          pstr++;
          code=strtol(*(data+pstr),NULL,10);
          pstr++;

          for(i=0;i<(org->nprop);i++) /*FIND PROP.*/
          {
            if(code==((org->props+i)->code))
            {
              ((org->sects+psect)->figs+psfig)->prop
              =(org->props+i);
              break;
            }
          }
          ident++;
        }
        if(nstr-pstr>=2 &&
           (org->nsect)-psect>0 &&
           !strcmp(*(data+pstr),"AREA") &&
           ((org->sects+psect)->nfig)-psfig>0)
        {
          pstr++;
          ((org->sects+psect)->figs+psfig)->area
          =strtod(*(data+pstr),NULL);
          pstr++;
          ident++;
        }
        if(nstr-pstr>=2 &&
           (org->nsect)-psect>0 &&
           !strcmp(*(data+pstr),"IXX") &&
           ((org->sects+psect)->nfig)-psfig>0)
        {
          pstr++;
          ((org->sects+psect)->figs+psfig)->Ixx
          =strtod(*(data+pstr),NULL);
          pstr++;
          ident++;
        }
        if(nstr-pstr>=2 &&
           (org->nsect)-psect>0 &&
           !strcmp(*(data+pstr),"IYY") &&
           ((org->sects+psect)->nfig)-psfig>0)
        {
          pstr++;
          ((org->sects+psect)->figs+psfig)->Iyy
          =strtod(*(data+pstr),NULL);
          pstr++;
          ident++;
        }
        if(nstr-pstr>=2 &&
           (org->nsect)-psect>0 &&
           !strcmp(*(data+pstr),"VEN") &&
           ((org->sects+psect)->nfig)-psfig>0)
        {
          pstr++;
          ((org->sects+psect)->figs+psfig)->Jzz
          =strtod(*(data+pstr),NULL);
          pstr++;
          ident++;
        }
        if(nstr-pstr>=2 &&
           (org->nsect)-psect>0 &&
           !strcmp(*(data+pstr),"THICK") &&
           ((org->sects+psect)->nfig)-psfig>0)
        {
          pstr++;
          ((org->sects+psect)->figs+psfig)->thick
          =strtod(*(data+pstr),NULL);
          pstr++;
          ident++;
        }

        if(nstr-pstr>=2 && (org->nsect)-psect>0 &&
           !strcmp(*(data+pstr),"EXP"))
        {
          pstr++;
          (org->sects+psect)->exp=strtod(*(data+pstr),NULL);
          pstr++;
          ident++;
        }
        if(nstr-pstr>=2 && (org->nsect)-psect>0 &&
           !strcmp(*(data+pstr),"NZMAX"))
        {
          pstr++;
          (org->sects+psect)->fmax[0]=strtod(*(data+pstr),NULL);
          pstr++;
          ident++;
        }
        if(nstr-pstr>=2 && (org->nsect)-psect>0 &&
           !strcmp(*(data+pstr),"NZMIN"))
        {
          pstr++;
          (org->sects+psect)->fmin[0]=strtod(*(data+pstr),NULL);
          pstr++;
          ident++;
        }
        if(nstr-pstr>=2 && (org->nsect)-psect>0 &&
           !strcmp(*(data+pstr),"QXMAX"))
        {
          pstr++;
          (org->sects+psect)->fmax[1]=strtod(*(data+pstr),NULL);
          pstr++;
          ident++;
        }
        if(nstr-pstr>=2 && (org->nsect)-psect>0 &&
           !strcmp(*(data+pstr),"QXMIN"))
        {
          pstr++;
          (org->sects+psect)->fmin[1]=strtod(*(data+pstr),NULL);
          pstr++;
          ident++;
        }
        if(nstr-pstr>=2 && (org->nsect)-psect>0 &&
           !strcmp(*(data+pstr),"QYMAX"))
        {
          pstr++;
          (org->sects+psect)->fmax[2]=strtod(*(data+pstr),NULL);
          pstr++;
          ident++;
        }
        if(nstr-pstr>=2 && (org->nsect)-psect>0 &&
           !strcmp(*(data+pstr),"QYMIN"))
        {
          pstr++;
          (org->sects+psect)->fmin[2]=strtod(*(data+pstr),NULL);
          pstr++;
          ident++;
        }
        if(nstr-pstr>=2 && (org->nsect)-psect>0 &&
           !strcmp(*(data+pstr),"MZMAX"))
        {
          pstr++;
          (org->sects+psect)->fmax[3]=strtod(*(data+pstr),NULL);
          pstr++;
          ident++;
        }
        if(nstr-pstr>=2 && (org->nsect)-psect>0 &&
           !strcmp(*(data+pstr),"MZMIN"))
        {
          pstr++;
          (org->sects+psect)->fmin[3]=strtod(*(data+pstr),NULL);
          pstr++;
          ident++;
        }
        if(nstr-pstr>=2 && (org->nsect)-psect>0 &&
           !strcmp(*(data+pstr),"MXMAX"))
        {
          pstr++;
          (org->sects+psect)->fmax[4]=strtod(*(data+pstr),NULL);
          pstr++;
          ident++;
        }
        if(nstr-pstr>=2 && (org->nsect)-psect>0 &&
           !strcmp(*(data+pstr),"MXMIN"))
        {
          pstr++;
          (org->sects+psect)->fmin[4]=strtod(*(data+pstr),NULL);
          pstr++;
          ident++;
        }
        if(nstr-pstr>=2 && (org->nsect)-psect>0 &&
           !strcmp(*(data+pstr),"MYMAX"))
        {
          pstr++;
          (org->sects+psect)->fmax[5]=strtod(*(data+pstr),NULL);
          pstr++;
          ident++;
        }
        if(nstr-pstr>=2 && (org->nsect)-psect>0 &&
           !strcmp(*(data+pstr),"MYMIN"))
        {
          pstr++;
          (org->sects+psect)->fmin[5]=strtod(*(data+pstr),NULL);
          pstr++;
          ident++;
        }

        /*NODES:CODE,COORDS,CONFS.*/
        if(nstr-pstr>=2 && (org->nnode)-pnode>1 &&
           !strcmp(*(data+pstr),"NODE"))
        {
          pstr++;
          pnode++;
          (org->nodes+pnode)->code=strtol(*(data+pstr),NULL,10);
          (org->nodes+pnode)->loff=pnode;
          pstr++;
          ident++;
        }
        if(nstr-pstr>=4 && (org->nnode)-pnode>0 &&
           !strcmp(*(data+pstr),"CORD")) /*COORDINATION X,Y,Z.*/
        {
          pstr++;
          (org->nodes+pnode)->d[0]=strtod(*(data+pstr+0),NULL);
          (org->nodes+pnode)->d[1]=strtod(*(data+pstr+1),NULL);
          (org->nodes+pnode)->d[2]=strtod(*(data+pstr+2),NULL);
          pstr+=3;
          ident++;
        }
        if(nstr-pstr>=7 && (org->nnode)-pnode>0 &&
           !strcmp(*(data+pstr),"ICON")) /*CONFINEMENT ID.*/
        {
          pstr++;
          for(i=0;i<=5;i++)
          {
            (org->confs+6*((org->nodes+pnode)->loff)+i)->iconf
            =(signed char)strtol(*(data+pstr+i),NULL,10);
          }
          pstr+=6;
          ident++;
        }
        if(nstr-pstr>=7 && (org->nnode)-pnode>0 &&
           !strcmp(*(data+pstr),"VCON")) /*CONFINEMENT VALUE.*/
        {
          pstr++;
          for(i=0;i<=5;i++)
          {
            (org->confs+6*((org->nodes+pnode)->loff)+i)->value
            =strtod(*(data+pstr+i),NULL);
          }
          pstr+=6;
          ident++;
        }

        /*ELEMENTS:CODE,SECT,NODES,BONDS,BANS,COORDANGLE,TYPE,ROLE.*/
        if(nstr-pstr>=2 && (org->nelem)-pelem>1 &&
           !strcmp(*(data+pstr),"ELEM"))
        {
          pstr++;
          pelem++;
          (org->elems+pelem)->code=strtol(*(data+pstr),NULL,10);
          (org->elems+pelem)->loff=pelem;
          (org->elems+pelem)->type=0;
          (org->elems+pelem)->nnod=0;
          (org->elems+pelem)->nban=0;
          (org->elems+pelem)->color.line=setrgbcolor(255,255,255);
          pstr++;
          ident++;
        }
        if(nstr-pstr>=2 && (org->nelem)-pelem>0 &&
           !strcmp(*(data+pstr),"ESECT")) /*SECTION POINTER.*/
        {
          pstr++;
          code=strtol(*(data+pstr),NULL,10);
          pstr++;

          for(i=0;i<(org->nsect);i++) /*FIND SECT.*/
          {
            if(code==((org->sects+i)->code))
            {
              (org->elems+pelem)->sect=(org->sects+i);
              break;
            }
          }
          ident++;
        }
        if(nstr-pstr>=2 && (org->nelem)-pelem>0 &&
           !strcmp(*(data+pstr),"CANG")) /*COORD ANGLE FOR WIRE.*/
        {
          pstr++;
          (org->elems+pelem)->cangle=strtod(*(data+pstr),NULL);
          pstr++;
          ident++;
        }
        if(nstr-pstr>=13 && (org->nelem)-pelem>0 &&
           !strcmp(*(data+pstr),"CMQ")) /*CMQ AS INITIAL.*/
        {
          pstr++;
          for(i=0;i<=1;i++)
          {
            for(j=0;j<=5;j++)
            {
              (org->elems+pelem)->initial[i][j]
              =strtod(*(data+pstr+6*i+j),NULL);
            }
          }
          pstr+=12;
          ident++;
        }
        if(nstr-pstr>=2 && (org->nelem)-pelem>0 &&
           !strcmp(*(data+pstr),"TYPE")) /*ELEMENT TYPE.*/
        {
          pstr++;
          if(!strcmp(*(data+pstr),"WALL"))
          {
            (org->elems+pelem)->type=WALL;
            (org->elems+pelem)->color.line=setrgbcolor(150,100,255);
          }
          else if(!strcmp(*(data+pstr),"SLAB"))
          {
            (org->elems+pelem)->type=SLAB;
            (org->elems+pelem)->color.line=setrgbcolor(255,0,255);
          }
          pstr++;
          ident++;
        }
        if(nstr-pstr>=2 && (org->nelem)-pelem>0 &&
           !strcmp(*(data+pstr),"ROLE")) /*ELEMENT ROLE.*/
        {
          pstr++;
          if(!strcmp(*(data+pstr),"W"))
          {
            (org->elems+pelem)->role=ROLEWEIGHT;
          }
          else if(!strcmp(*(data+pstr),"R"))
          {
            (org->elems+pelem)->role=ROLERIGID;
          }
          else if(!strcmp(*(data+pstr),"S"))
          {
            (org->elems+pelem)->role=ROLESTRESS;
          }
          else
          {
            (org->elems+pelem)->role=0;
          }
          pstr++;
          ident++;
        }
        if(nstr-pstr>=2 && (org->nelem)-pelem>0 &&
           !strcmp(*(data+pstr),"ENODS"))
        {
          pstr++;
          (org->elems+pelem)->nnod=strtol(*(data+pstr),NULL,10);
          (org->elems+pelem)->nods=(struct onode **)
                                   malloc((org->elems+pelem)->nnod
                                          *sizeof(struct onode *));
          (org->elems+pelem)->bonds=(signed char *)
                                    malloc(6*(org->elems+pelem)->nnod
                                           *sizeof(char));
          pstr++;
          ident++;
        }
        if(nstr-pstr>=2 &&
           (org->nelem)-pelem>0 &&
           !strcmp(*(data+pstr),"ENOD") &&
           nstr-pstr>=((org->elems+pelem)->nnod+1))
        {
          pstr++;
          for(i=0;i<(org->elems+pelem)->nnod;i++)
          {
            code=strtol(*(data+pstr+i),NULL,10);

            for(j=0;j<(org->nnode);j++) /*FIND NODE.*/
            {
              if(code==((org->nodes+j)->code))
              {
                *((org->elems+pelem)->nods+i)=(org->nodes+j);
                break;
              }
            }
          }
          pstr+=(org->elems+pelem)->nnod;
          ident++;
        }
        if(nstr-pstr>=2 &&
           (org->nelem)-pelem>0 &&
           !strcmp(*(data+pstr),"BONDS") &&
           nstr-pstr>=(6*((org->elems+pelem)->nnod)+1))
        {
          pstr++;
          for(i=0;i<6*((org->elems+pelem)->nnod);i++)
          {
            *((org->elems+pelem)->bonds+i)
            =(signed char)strtol(*(data+pstr+i),NULL,10);
          }
          pstr+=6*((org->elems+pelem)->nnod);
          ident++;
        }
        if(nstr-pstr>=2 && (org->nelem)-pelem>0 &&
           !strcmp(*(data+pstr),"EBANS"))
        {
          pstr++;
          peban=(-1);
          (org->elems+pelem)->nban=strtol(*(data+pstr),NULL,10);
          (org->elems+pelem)->bans=(struct obans *)
                                   malloc((org->elems+pelem)->nban
                                          *sizeof(struct obans));
          pstr++;
          ident++;
        }
        if(nstr-pstr>=2 &&
           (org->nelem)-pelem>0 &&
           !strcmp(*(data+pstr),"EBAN") &&
           ((org->elems+pelem)->nban)-peban>1)
        {
          pstr++;
          peban++;
          ((org->elems+pelem)->bans+peban)->code
          =strtol(*(data+pstr),NULL,10);
          pstr++;
          ident++;
        }
        if(nstr-pstr>=2 &&
           (org->nelem)-pelem>0 &&
           !strcmp(*(data+pstr),"BNODS") &&
           ((org->elems+pelem)->nban)-peban>0)
        {
          pstr++;
          ((org->elems+pelem)->bans+peban)->nnod
          =strtol(*(data+pstr),NULL,10);
          ((org->elems+pelem)->bans+peban)->nods
          =(struct onode **)
           malloc(((org->elems+pelem)->bans+peban)->nnod
                  *sizeof(struct onode *));
          pstr++;
          ident++;
        }
        if(nstr-pstr>=2 &&
           (org->nelem)-pelem>0 &&
           !strcmp(*(data+pstr),"BNOD") &&
           ((org->elems+pelem)->nban)-peban>0 &&
           nstr-pstr>=(((org->elems+pelem)->bans+peban)->nnod+1))
        {
          pstr++;
          for(i=0;i<(((org->elems+pelem)->bans+peban)->nnod);i++)
          {
            code=strtol(*(data+pstr+i),NULL,10);

            for(j=0;j<((org->elems+pelem)->nnod);j++) /*FIND NODE.*/
            {
              if(code==((*((org->elems+pelem)->nods+j))->code))
              {
                *(((org->elems+pelem)->bans+peban)->nods+i)
                =*((org->elems+pelem)->nods+j);
                break;
              }
            }
          }
          pstr+=((org->elems+pelem)->bans+peban)->nnod;
          ident++;
        }
        if(ident==0) /*UNIDENTIFIED.*/
        {
          pstr++;
        }
      }
    }
    freestr(data,nstr);
  }
}/*inputorganization*/

int saveorganization(FILE *fout,struct organ *org)
/*SAVE ORGAN INTO OUTPUT FILE.*/
{
  int i,j,k;
  long int nnode,nelem,nprop,nsect;

  fseek(fout,0L,SEEK_SET);

  fprintf(fout,"\"CREATED ORGAN FRAME.\"\n");

  nnode=org->nnode;
  nelem=org->nelem;
  nprop=org->nprop;
  nsect=org->nsect;

  fprintf(fout,"NNODE %ld\n",nnode);
  fprintf(fout,"NELEM %ld\n",nelem);
  fprintf(fout,"NPROP %ld\n",nprop);
  fprintf(fout,"NSECT %ld\n",nsect);
  fprintf(fout,"\n");

  for(i=0;i<nprop;i++)
  {
    fprintf(fout,"PROP %3ld ",(org->props+i)->code);

    if((org->props+i)->name!=NULL)
    {
      fprintf(fout,"PNAME %s\n",(org->props+i)->name);
    }
    else
    {
      fprintf(fout,"\"None\"\n");
    }

    fprintf(fout,"         HIJU %12.3f\n",(org->props+i)->hiju);
    fprintf(fout,"         E    %12.3f\n",(org->props+i)->E);
    fprintf(fout,"         POI  %12.5f\n",(org->props+i)->poi);
  }
  fprintf(fout,"\n");

  for(i=0;i<nsect;i++)
  {
    fprintf(fout,"SECT %3ld ",(org->sects+i)->code);

    if((org->sects+i)->name!=NULL)
    {
      fprintf(fout,"SNAME %s\n",(org->sects+i)->name);
    }
    else
    {
      fprintf(fout,"\"None\"\n");
    }

    fprintf(fout,"         NFIG %d\n",(org->sects+i)->nfig);
    for(j=0;j<(org->sects+i)->nfig;j++)
    {
      fprintf(fout,"         FIG %3d ",
              ((org->sects+i)->figs+j)->code);
      fprintf(fout,"FPROP %d\n",
              ((org->sects+i)->figs+j)->prop->code);

      fprintf(fout,"                 ");
      fprintf(fout,"AREA  %.4f\n",((org->sects+i)->figs+j)->area);
      fprintf(fout,"                 ");
      fprintf(fout,"IXX   %.8f\n",((org->sects+i)->figs+j)->Ixx);
      fprintf(fout,"                 ");
      fprintf(fout,"IYY   %.8f\n",((org->sects+i)->figs+j)->Iyy);
      fprintf(fout,"                 ");
      fprintf(fout,"VEN   %.8f\n",((org->sects+i)->figs+j)->Jzz);
      fprintf(fout,"                 ");
      fprintf(fout,"THICK %.3f\n",((org->sects+i)->figs+j)->thick);
    }

    fprintf(fout,"         EXP %.3f\n",(org->sects+i)->exp);

    fprintf(fout,"         ");
    fprintf(fout,"NZMAX %8.1f ", (org->sects+i)->fmax[0]);
    fprintf(fout,"NZMIN %8.1f\n",(org->sects+i)->fmin[0]);
    fprintf(fout,"         ");
    fprintf(fout,"QXMAX %8.1f ", (org->sects+i)->fmax[1]);
    fprintf(fout,"QXMIN %8.1f\n",(org->sects+i)->fmin[1]);
    fprintf(fout,"         ");
    fprintf(fout,"QYMAX %8.1f ", (org->sects+i)->fmax[2]);
    fprintf(fout,"QYMIN %8.1f\n",(org->sects+i)->fmin[2]);
    fprintf(fout,"         ");
    fprintf(fout,"MZMAX %8.1f ", (org->sects+i)->fmax[3]);
    fprintf(fout,"MZMIN %8.1f\n",(org->sects+i)->fmin[3]);
    fprintf(fout,"         ");
    fprintf(fout,"MXMAX %8.1f ", (org->sects+i)->fmax[4]);
    fprintf(fout,"MXMIN %8.1f\n",(org->sects+i)->fmin[4]);
    fprintf(fout,"         ");
    fprintf(fout,"MYMAX %8.1f ", (org->sects+i)->fmax[5]);
    fprintf(fout,"MYMIN %8.1f\n",(org->sects+i)->fmin[5]);
  }
  fprintf(fout,"\n");

  for(i=0;i<nnode;i++)
  {
    fprintf(fout,"NODE %3ld  ",(org->nodes+i)->code);
    fprintf(fout,"CORD %7.3f %7.3f %7.3f  ",
            (org->nodes+i)->d[0],
            (org->nodes+i)->d[1],
            (org->nodes+i)->d[2]);
    fprintf(fout,"ICON %1d %1d %1d %1d %1d %1d  ",
           (org->confs+6*i+0)->iconf,
           (org->confs+6*i+1)->iconf,
           (org->confs+6*i+2)->iconf,
           (org->confs+6*i+3)->iconf,
           (org->confs+6*i+4)->iconf,
           (org->confs+6*i+5)->iconf);
    fprintf(fout,"VCON %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f\n",
           (org->confs+6*i+0)->value,
           (org->confs+6*i+1)->value,
           (org->confs+6*i+2)->value,
           (org->confs+6*i+3)->value,
           (org->confs+6*i+4)->value,
           (org->confs+6*i+5)->value);
  }
  fprintf(fout,"\n");

  for(i=0;i<nelem;i++)
  {
    fprintf(fout,"ELEM %3ld ",(org->elems+i)->code);
    fprintf(fout,"ESECT %3ld ",(org->elems+i)->sect->code);

    fprintf(fout,"ENODS %d ",(org->elems+i)->nnod);
    fprintf(fout,"ENOD");
    for(j=0;j<(org->elems+i)->nnod;j++)
    {
      fprintf(fout," %d",(*((org->elems+i)->nods+j))->code);
    }
    fprintf(fout," BONDS");
    for(j=0;j<((org->elems+i)->nnod);j++)
    {
      fprintf(fout,"  %1d %1d %1d %1d %1d %1d",
              *((org->elems+i)->bonds+6*j+0),
              *((org->elems+i)->bonds+6*j+1),
              *((org->elems+i)->bonds+6*j+2),
              *((org->elems+i)->bonds+6*j+3),
              *((org->elems+i)->bonds+6*j+4),
              *((org->elems+i)->bonds+6*j+5));
    }
    fprintf(fout,"\n");

    if(((org->elems+i)->nban)>=1)
    {
      fprintf(fout,"                   ");
      fprintf(fout,"EBANS %d ",(org->elems+i)->nban);
      for(j=0;j<(org->elems+i)->nban;j++)
      {
        if(j>=1) fprintf(fout,"                           ");

        fprintf(fout,"EBAN %ld ",((org->elems+i)->bans+j)->code);
        fprintf(fout,"BNODS %d ",((org->elems+i)->bans+j)->nnod);
        fprintf(fout,"BNOD");
        for(k=0;k<((org->elems+i)->bans+j)->nnod;k++)
        {
          fprintf(fout," %ld",
                  (*(((org->elems+i)->bans+j)->nods+k))->code);
        }

        fprintf(fout,"\n");
      }
    }

    if(((org->elems+i)->nnod)==2)
    {
      fprintf(fout,"         CANG %.5f\n",(org->elems+i)->cangle);

      fprintf(fout,"         CMQ ");
      for(j=0;j<2;j++)
      {
        for(k=0;k<6;k++)
        {
          fprintf(fout," %.1f",(org->elems+i)->initial[j][k]);
        }
        fprintf(fout," ");
      }
      fprintf(fout,"\n");
    }

    if((org->elems+i)->type==WALL)
    {
      fprintf(fout,"         TYPE WALL\n");
    }
    if((org->elems+i)->type==SLAB)
    {
      fprintf(fout,"         TYPE SLAB\n");
    }

    if((org->elems+i)->role==ROLEWEIGHT)
    {
      fprintf(fout,"         ROLE W\n");
    }
    if((org->elems+i)->role==ROLERIGID)
    {
      fprintf(fout,"         ROLE R\n");
    }
    if((org->elems+i)->role==ROLESTRESS)
    {
      fprintf(fout,"         ROLE S\n");
    }
  }

  return 1;
}/*saveorganization*/

void freeorganization(struct organ *org)
/*FREE ORGAN FROM MEMORY.*/
{
  int i,j;

  org->code=0;
  org->loff=0;

  for(i=0;i<(org->nelem);i++) /*ELEMENTS.*/
  {
    free((org->elems+i)->nods);
    free((org->elems+i)->bonds);

    for(j=0;j<((org->elems+i)->nban);j++)
    {
      free(((org->elems+i)->bans+j)->nods);
    }
    free((org->elems+i)->bans);
  }
  free(org->elems);
  org->nelem=0;

  free(org->nodes); /*NODES.*/
  free(org->confs);
  org->nnode=0;

  for(i=0;i<(org->nsect);i++) /*SECTIONS.*/
  {
    free((org->sects+i)->figs);
  }
  free(org->sects);
  org->nsect=0;

  free(org->props); /*PROPERTIES.*/
  org->nprop=0;

  return;
}/*freeorganization*/

void slabdivision(HDC hdc,struct viewparam vp,struct obans gban)
/*DIVIDE SLAB BAN INTO CMQ PARTS.*/
{
  int ndot,bdot,hdot,idot,jdot,kdot;
  int iidot;
  int ibegin;
  double param0,param1,param2,distance;
  double **drccos,**tdrccos;
  struct onode *lnods;
  struct obans ban;
  struct obans *dban; /*DIVIDED BANS.*/
  struct line dline,rline,*hen0,*hen1;
  struct cmqdot *ddot;

  if(hdc==NULL || gban.nods==NULL) return;

  ndot=gban.nnod;

  ban.nods=(struct onode **)malloc(ndot*sizeof(struct onode *));
  lnods=(struct onode *)malloc(ndot*sizeof(struct onode));

  hen0=(struct line *)malloc(ndot*sizeof(struct line));
  hen1=(struct line *)malloc(ndot*sizeof(struct line));

  dban=(struct obans *)malloc(ndot*sizeof(struct obans));
  ddot=(struct cmqdot *)malloc(ndot*sizeof(struct cmqdot));

  if(ban.nods==NULL ||
     lnods==NULL ||
     hen0==NULL ||
     hen1==NULL ||
     dban==NULL ||
     ddot==NULL) return;

  drccos=filmdrccos(*(*(gban.nods+0)),
                    *(*(gban.nods+1)),
                    *(*(gban.nods+ndot-1)));    /*DIRECTION COSINE.*/
  tdrccos=matrixtranspose(drccos,3);                  /*MATRIX [T].*/

  for(idot=0;idot<ndot;idot++) /*TRANSFORM GLOBAL INTO LOCAL.*/
  {
    *(lnods+idot)=transgtol(drccos,*(*(gban.nods+0)),
                                   *(*(gban.nods+idot)));
    *(ban.nods+idot)=(lnods+idot);

    if(idot>=1)
    {
      dline.ends[0].d[EX]=(lnods+idot-1)->d[EX];
      dline.ends[0].d[EY]=(lnods+idot-1)->d[EY];
      dline.ends[0].d[EZ]=0.0; /*PROJECT TO xy PLANE.*/
      dline.ends[1].d[EX]=(lnods+idot)->d[EX];
      dline.ends[1].d[EY]=(lnods+idot)->d[EY];
      dline.ends[1].d[EZ]=0.0;

      *(hen0+idot-1)=dline; /*HEN OF GBAN ON LOCAL PLANE.*/
    }
    if(idot==(ndot-1))
    {
      dline.ends[0].d[EX]=(lnods+idot)->d[EX];
      dline.ends[0].d[EY]=(lnods+idot)->d[EY];
      dline.ends[0].d[EZ]=0.0;
      dline.ends[1].d[EX]=(lnods+0)->d[EX];
      dline.ends[1].d[EY]=(lnods+0)->d[EY];
      dline.ends[1].d[EZ]=0.0;

      *(hen0+idot)=dline; /*LAST LINE.*/
    }
  }

  for(idot=0;idot<ndot;idot++) /*FIRST DIVIDING LINES.*/
  {
    if(idot==0) hdot=ndot-1;
    else        hdot=idot-1;

    rline.ends[0].d[EX]=(hen0+hdot)->ends[1].d[EX]; /*REVERSE ENDS.*/
    rline.ends[0].d[EY]=(hen0+hdot)->ends[1].d[EY];
    rline.ends[0].d[EZ]=(hen0+hdot)->ends[1].d[EZ];
    rline.ends[1].d[EX]=(hen0+hdot)->ends[0].d[EX];
    rline.ends[1].d[EY]=(hen0+hdot)->ends[0].d[EY];
    rline.ends[1].d[EZ]=(hen0+hdot)->ends[0].d[EZ];

    divisionlineline(rline,*(hen0+idot),&dline);

    *(hen1+idot)=dline; /*DIVIDING LINES FROM VERTICES.*/
  }

  for(idot=0;idot<ndot;idot++) /*DIVIDING DOTS.*/
  {
    if(idot==0) hdot=ndot-1;
    else        hdot=idot-1;

    (ddot+idot)->parents[0]=(hen0+hdot);
    (ddot+idot)->parents[1]=(hen0+idot);

    ibegin=0;
    for(iidot=0;iidot<ndot;iidot++)
    {
      if(iidot!=idot && iidot!=hdot)
      {
        ibegin++;

        rline.ends[0].d[EX]=(hen0+idot)->ends[1].d[EX];
        rline.ends[0].d[EY]=(hen0+idot)->ends[1].d[EY];
        rline.ends[0].d[EZ]=(hen0+idot)->ends[1].d[EZ];
        rline.ends[1].d[EX]=(hen0+idot)->ends[0].d[EX];
        rline.ends[1].d[EY]=(hen0+idot)->ends[0].d[EY];
        rline.ends[1].d[EZ]=(hen0+idot)->ends[0].d[EZ];

        divisionlineline(rline,*(hen0+iidot),
                         &dline);
        distancelineline(*(hen1+idot),dline,
                         &param1,&param2,&distance);

        if(ibegin==1)
        {
          param0=param1;
          (ddot+idot)->parents[2]=(hen0+iidot);
        }
        else if(fabs(param1)<fabs(param0))
        {
          param0=param1;
          (ddot+idot)->parents[2]=(hen0+iidot);
        }

        rline.ends[0].d[EX]=(hen0+hdot)->ends[1].d[EX];
        rline.ends[0].d[EY]=(hen0+hdot)->ends[1].d[EY];
        rline.ends[0].d[EZ]=(hen0+hdot)->ends[1].d[EZ];
        rline.ends[1].d[EX]=(hen0+hdot)->ends[0].d[EX];
        rline.ends[1].d[EY]=(hen0+hdot)->ends[0].d[EY];
        rline.ends[1].d[EZ]=(hen0+hdot)->ends[0].d[EZ];

        divisionlineline(rline,*(hen0+iidot),
                         &dline);
        distancelineline(*(hen1+idot),dline,
                         &param1,&param2,&distance);

        if(fabs(param1)<fabs(param0))
        {
          param0=param1;
          (ddot+idot)->parents[2]=(hen0+iidot);
        }
      }
    }

    (ddot+idot)->n.d[EX]=(1.0-param0)
                         *((hen1+idot)->ends[0].d[EX])
                         +param0
                         *((hen1+idot)->ends[1].d[EX]);
    (ddot+idot)->n.d[EY]=(1.0-param0)
                         *((hen1+idot)->ends[0].d[EY])
                         +param0
                         *((hen1+idot)->ends[1].d[EY]);
    (ddot+idot)->n.d[EZ]=(1.0-param0)
                         *((hen1+idot)->ends[0].d[EZ])
                         +param0
                         *((hen1+idot)->ends[1].d[EZ]);

    /*drawlocalline(hdc,vp,tdrccos,*(*(gban.nods+0)),
                                 *(lnods+idot),
                                 (ddot+idot)->n);*/
  }

  for(idot=0;idot<ndot;idot++) /*DICIDING DIVIDED BANS.*/
  {
    if(idot==ndot-1) jdot=0;
    else             jdot=idot+1;

    bdot=2; /*1st AND 2nd DOTS=ENDS ON PARENT HEN OF BAN.*/
    (dban+idot)->nnod=bdot;
    (dban+idot)->nods=(struct onode **)
                      malloc(bdot*sizeof(struct onode *));

    if((dban+idot)->nods==NULL)
    {
      MessageBox(NULL,"Malloc Failed.","Division",MB_OK);
      return;
    }

    *((dban+idot)->nods+0)=(lnods+idot); /*1st DOT.*/
    *((dban+idot)->nods+1)=(lnods+jdot); /*2nd DOT.*/

    for(iidot=1;iidot<=ndot;iidot++)
    {
      kdot=idot+iidot;           /*CODE OF DDOT.*/
      if(kdot>=ndot) kdot-=ndot;

      if((ddot+kdot)->parents[0]==(hen0+idot) ||
         (ddot+kdot)->parents[1]==(hen0+idot) ||
         (ddot+kdot)->parents[2]==(hen0+idot))
      {
        bdot++;
        (dban+idot)->nnod=bdot;
        (dban+idot)->nods=(struct onode **)
                          realloc((dban+idot)->nods,
                          bdot*sizeof(struct onode *));
        if((dban+idot)->nods==NULL)
        {
          MessageBox(NULL,"Realloc Failed.","Division",MB_OK);
          return;
        }

        *((dban+idot)->nods+(bdot-1))=&((ddot+kdot)->n);
      }/*DIVIDED BAN OF HEN CONSISTS OF DDOT HAVING HEN AS PARENT.*/
    }

    drawlocalban(hdc,vp,tdrccos,*(*(gban.nods+0)),*(dban+idot));
  }

  free(lnods);
  free(ban.nods);
  free(hen0);
  free(hen1);
  free(ddot);
  for(idot=0;idot<ndot;idot++) free((dban+idot)->nods);
  free(dban);

  free(drccos);
  free(tdrccos);

  return;
}/*slabdivision*/

struct rgbcolor setrgbcolor(int r,int g,int b)
{
  struct rgbcolor c;

  c.r=r;
  c.g=g;
  c.b=b;

  return c;
}/*setrgbcolor*/

clock_t laptime(char *comment,clock_t t0)
/*RETURN:TIME PASSAGE FROM BEGINNING [sec].*/
{
  char s[80];
  clock_t t1;
  long int time;

  t1=clock();
  time=(t1-t0)/CLK_TCK;
  sprintf(s,"%sTIME:%ld[sec]",comment,time);
  errormessage(s);

  return t1;
}/*laptime*/

void getincrement(HWND hdwnd,int *laps,double *dsafety)
/*GET INCREMENT PARAMETERS FROM DIALOG.*/
{
  char data[80];

  GetDlgItemText(hdwnd,ID_LAPS,data,80);
  *laps=(int)strtol(data,NULL,10);
  GetDlgItemText(hdwnd,ID_SAFETY,data,80);
  *dsafety=strtod(data,NULL);

  return;
}/*getincrement*/

double croutludecomposition(struct gcomponent *gmtx,
                            double *gvct,struct oconf *confs,
                            long int msize)
/*SOLVE SIMULTANEOUS LINEAR EQUATIONS.*/
{
  /*char iconf;*/ /*0:FREE 1:FIXED*/
  long int i,j,k;
  double det=1.0;
  double data1;
  struct gcomponent *pivot,*pcomp;
  struct gcomponent *gcomp1,*gcomp2,*gcomp3,*gcomp4;

  for(j=1;j<=msize;j++)                            /*DECOMPOSITION.*/
  {
    if((confs+j-1)->iconf==0) /*FREE*/
    {
      pivot=(gmtx+(j-1)); /*PIVOT.*/

      if((pivot->value)==0.0) return 0.0;               /*INSTABLE.*/
      /*det*=pivot->value;*/
      det*=pivot->value/fabs(pivot->value);  /*SIGN OF DETERMINANT.*/

      gcomp1=pivot;
      while(gcomp1->down!=NULL) /*DOWNWARD.*/
      {
        gcomp1=gcomp1->down;

        i=gcomp1->m; /*i:LINE CODE IN PIVOT ROW.*/
        if((confs+i-1)->iconf==0) /*FREE*/
        {
          gcomp1->value/=(pivot->value);
        }
      }

      gcomp1=pivot;
      while(gcomp1->down!=NULL) /*DOWNWARD.*/
      {
        gcomp1=gcomp1->down; /*Aij*/

        i=gcomp1->m; /*i:LINE CODE IN PIVOT ROW.*/
        if((confs+i-1)->iconf==0) /*FREE*/
        {
          gcomp2=gcomp1; /*Akj*/
          gcomp3=(gmtx+(i-1)); /*Aki*/
          while(1)
          {
            k=gcomp2->m;
            if((confs+k-1)->iconf==0) /*FREE*/
            {
              if((gcomp3->m)<k) /*ADD*/
              {
                gcomp4=(struct gcomponent *)
                       malloc(sizeof(struct gcomponent));
                if(gcomp4==NULL)
                {
                  errormessage("CROUT:MEMORY INSUFFICIENT.");
                  return 0.0;
                }
                gcomp3->down=gcomp4;
                gcomp4->m=(unsigned short int)k;
                /*gcomp4->n=(unsigned short int)i;*/
                gcomp4->down=NULL;
                gcomp4->value=-(pivot->value)
                             *(gcomp1->value)
                             *(gcomp2->value);
                gcomp3=gcomp4;

                comps++;
              }
              else if((gcomp3->m)==k)
              {
                gcomp3->value-=(pivot->value)
                              *(gcomp1->value)
                              *(gcomp2->value);
              }
              else if((gcomp3->m)>k) /*FILL*/
              {
                gcomp4=(struct gcomponent *)
                       malloc(sizeof(struct gcomponent));
                if(gcomp4==NULL)
                {
                  errormessage("CROUT:MEMORY INSUFFICIENT.");
                  return 0.0;
                }

                pcomp->down=gcomp4;
                gcomp4->m=(unsigned short int)k;
                /*gcomp4->n=(unsigned short int)i;*/
                gcomp4->down=gcomp3;
                gcomp4->value=-(pivot->value)
                             *(gcomp1->value)
                             *(gcomp2->value);
                pcomp=gcomp4;

                comps++;
              }
            }

            if(gcomp2->down==NULL) break;
            else gcomp2=gcomp2->down;

            while((gcomp3->m)<(gcomp2->m) && gcomp3->down!=NULL)
            {
              pcomp=gcomp3;
              gcomp3=gcomp3->down;
            }
          }
        }
      }
      currentpivot(j,msize);
    }
  }

  errormessage("FORWARD ELIMINATION.");
  for(j=1;j<=msize;j++)                                  /*FORWARD.*/
  {
    if((confs+j-1)->iconf==0) /*FREE*/
    {
      data1=*(gvct+j-1);
      gcomp1=(gmtx+(j-1)); /*DIAGONAL.*/

      while(gcomp1->down!=NULL) /*DOWNWARD.*/
      {
        gcomp1=gcomp1->down;
        i=gcomp1->m;

        if((confs+i-1)->iconf==0) /*FREE*/
        {
          *(gvct+i-1)-=gcomp1->value*data1;
        }
      }
    }
    currentpivot(j,msize);
  }

  errormessage("BACKWARD SUBSTITUTION.");
  for(j=msize;j>=1;j--)                                 /*BACKWARD.*/
  {
    if((confs+j-1)->iconf==0) /*FREE*/
    {
      data1=*(gvct+j-1);
      gcomp1=(gmtx+(j-1)); /*DIAGONAL.*/
      data1/=gcomp1->value;

      while(gcomp1->down!=NULL) /*DOWNWARD.*/
      {
        gcomp1=gcomp1->down;
        i=gcomp1->m;

        if((confs+i-1)->iconf==0) /*FREE*/
        {
          data1-=gcomp1->value*(*(gvct+i-1));
        }
      }
      *(gvct+j-1)=data1;
    }
    currentpivot((j-1),msize);
  }

  return det;
}/*croutludecomposition*/

void currentpivot(long int i,long int msize)
/*WRITE PIVOT LINE INTO DIALOG BOX.*/
{
  HWND hdwnd=(wmenu.childs+2)->hwnd;
  char str[20];

  sprintf(str,"%4ld/%4ld",i,msize);
  SetDlgItemText(hdwnd,ID_PIVOT,str);
  SendDlgItemMessage(hdwnd,ID_PIVOT,WM_PAINT,0,0);

  return;
}/*currentpivot*/

int gread(struct gcomponent *gmtx,
          long int i,long int j,double *data)
/*READ DATA FROM FILE GLOBAL MATRIX.*/
{
  long int k;
  struct gcomponent *g;

  if(i<j) {k=i; i=j; j=k;}                       /*MATRIX SYMMETRIC*/

  g=(gmtx+(j-1));                                        /*DIAGONAL*/

  while((g->m)<i && (g->down)!=NULL) g=g->down;

  if((g->m)==i) *data=g->value;                         /*EXISTENT.*/
  else          *data=0.0;                                 /*EMPTY.*/

  return 1;
}/*gread*/

int gwrite(struct gcomponent *gmtx,
           long int i,long int j,double data)
/*WRITE DATA INTO FILE GLOBAL MATRIX.*/
{
  int n=0;
  long int k;
  struct gcomponent *gcomp,*gdown,*pdown;

  if(i<j) {k=i; i=j; j=k;}                       /*MATRIX SYMMETRIC*/

  gdown=(gmtx+(j-1));                                    /*DIAGONAL*/
  while((gdown->m)<i && (gdown->down)!=NULL)
  {
    pdown=gdown;
    gdown=gdown->down;
  }

  if(gdown->m==i)                             /*COMPONENT EXISTENT.*/
  {
    if(data==0.0 && i!=j)     /*COMPONENT VANISHED.EXCEPT DIAGONAL.*/
    {
      pdown->down=gdown->down;
      free(gdown);

      comps--;
      n=1;
    }
    else
    {
      gdown->value=data;

      n=2;
    }
  }

  else if(data!=0.0)                             /*COMPONENT EMPTY.*/
  {
    gcomp=(struct gcomponent *)malloc(sizeof(struct gcomponent));
    if(gcomp==NULL)
    {
      errormessage("GWRITE:MEMORY INSUFFICIENT.");
      return 0;
    }

    gcomp->m=(unsigned short int)i;
    /*gcomp->n=(unsigned short int)j;*/
    gcomp->value=data;

    if((gdown->m)<i)                                  /*ADD TO ROW.*/
    {
      gcomp->down=NULL;
      gdown->down=gcomp;

      n=3;
    }
    else if((gdown->m)>i)                          /*FILL INTO ROW.*/
    {
      gcomp->down=pdown->down;
      pdown->down=gcomp;

      n=4;
    }
    comps++;
  }

  return n;
}/*gwrite*/

void gfree(struct gcomponent *gmtx,long int nnode)
/*FREE GLOBAL MATRIX.*/
{
  long int n,i;
  struct gcomponent *g,*p;

  n=6*nnode;
  for(i=0;i<=n-1;i++)
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
}/*gfree*/

int vread(FILE *fvct,long int i,double *data)
/*READ DATA FROM FILE GLOBAL VECTOR.*/
{
  long int loffset;
  int n;

  loffset=(i-1)*sizeof(double);
  fseek(fvct,loffset,SEEK_SET);
  n=fread(data,sizeof(double),1,fvct);

  return n;
}/*vread*/

int vwrite(FILE *fvct,long int i,double *data)
/*WRITE DATA INTO FILE GLOBAL VECTOR.*/
{
  long int loffset;
  int n;

  loffset=(i-1)*sizeof(double);
  fseek(fvct,loffset,SEEK_SET);
  n=fwrite(data,sizeof(double),1,fvct);
  fflush(fvct);

  return n;
}/*vwrite*/

int arclm101(struct arclmframe *af)
{
  DWORD memory0,memory1,memory2;

  FILE *fin,*fout;                                   /*FILE 8 BYTES*/
  FILE *felem,*fdisp,*freact;
  char dir[]=DIRECTORY;                            /*DATA DIRECTORY*/
  char s[80],string[400];
  int i,ii,jj;
  int nnode,nelem,nsect,nlap,laps=LAPS,nreact;
  long int fsize,loffset,msize;
  long int time;
  struct gcomponent ginit={0,0,0.0,NULL};
  struct gcomponent *gmtx,*g,*p;                    /*GLOBAL MATRIX*/
  double *gvct;                                     /*GLOBAL VECTOR*/
  double **drccos,**tmatrix,**estiff,*estress;
  double determinant,data,safety,dsafety=DSAFETY;
  clock_t t0,t1,t2;

  struct osect *sects;
  struct onode *nodes,*ninit;
  struct owire elem;
  struct owire *elems;
  struct oconf *confs;

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

  getincrement((wmenu.childs+2)->hwnd,&laps,&dsafety);

  inputinit(fin,&nnode,&nelem,&nsect);             /*INPUT INITIAL.*/
  sprintf(string,"NODES=%d ELEMS=%d SECTS=%d",nnode,nelem,nsect);
  errormessage(string);
  fprintf(fout,"%s\n",string);

  msize=6*nnode;                           /*SIZE OF GLOBAL MATRIX.*/

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

  felem=fopen("canbin.elm","wb+");  /*CODE,12 BOUNDARIES,12 STRESS.*/
  af->felem=felem;
  fsize=sizeof(long int)+12*sizeof(signed char)+12*sizeof(double);
  setvbuf(felem,NULL,_IOFBF,fsize);

  free(af->sects);
  free(af->nodes);
  free(af->ninit);
  free(af->elems);
  free(af->confs);

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

  af->sects=sects;
  af->nodes=nodes;
  af->ninit=ninit;
  af->elems=elems;
  af->confs=confs;

  inputtexttomemory(fin,af);        /*READY TO READ LONG REACTIONS.*/
  nnode=af->nnode;
  nelem=af->nelem;
  nsect=af->nsect;
  nreact=af->nreact;

  initialform(nodes,fdisp,nnode);           /*ASSEMBLAGE FORMATION.*/

  initialelem(elems,felem,nelem);            /*ASSEMBLAGE ELEMENTS.*/

  freact=fopen("canbin.rct","wb+");                /*REACTION FILE.*/
  af->freact=freact;
  fsize=sizeof(double);
  setvbuf(freact,NULL,_IOFBF,fsize);
  initialreact(fin,freact,nreact);     /*ASSEMBLAGE LONG REACTIONS.*/

  GetAsyncKeyState(VK_LBUTTON);                   /*CLEAR KEY LEFT.*/
  GetAsyncKeyState(VK_RBUTTON);                  /*CLEAR KEY RIGHT.*/

  drawglobalaxis((wdraw.childs+1)->hdcC,(wdraw.childs+1)->vparam,
                 0,0,255);                      /*DRAW GLOBAL AXIS.*/

  if(wsurf.hwnd!=NULL)
  {
    drawyieldsurface((wsurf.childs+1)->hdcC,
                     (wsurf.childs+1)->vparam,2,4,0,NULL);
    overlayhdc(*(wsurf.childs+1),SRCPAINT);       /*UPDATE DISPLAY.*/
  }
  af->fsurface=fopen("canbin.sfc","wb+");      /*STRESS ON SURFACE.*/

  for(nlap=1;nlap<=laps;nlap++)
  {
    /*af->nlaps=nlap;*/
    af->nlaps=1;

    sprintf(string,"LAP:%d/%d",nlap,laps);
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
                                255,255,255,0);
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
      estiff=assempmtx(elem,estiff);          /*ADD PLASTIC MATRIX.*/
      estiff=transformation(estiff,tmatrix);       /*[K]=[Tt][k][T]*/

      assemgstiffness(gmtx,estiff,&elem);             /*ASSEMBLAGE.*/

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
    assemconf(confs,gvct,dsafety,nnode);           /*GLOBAL VECTOR.*/
    modifygivend(gmtx,gvct,confs,nnode);   /*0:LOAD 1:DISPLACEMENT.*/

    laptime("CROUT LU DECOMPOSITION.",t0);
    determinant=croutludecomposition(gmtx,
                                     gvct,confs,
                                     6*nnode);       /*[K]{dU}={dF}*/

    sprintf(string,"DETERMINANT=%.5E COMPS=%ld",determinant,comps);
    errormessage(string);
    if(fout!=NULL) fprintf(fout,"%s\n",string);

    safety=nlap*dsafety;
    sprintf(string,"SAFETY FACTOR=%.5f",safety);
    errormessage(string);
    if(fout!=NULL) fprintf(fout,"%s\n",string);

    if(determinant<=0.0)
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

    if(fout!=NULL) fprintf(fout,"\"DISPLACEMENT\"\n");
    outputdisp(gvct,fout,nnode,nodes);  /*INCREMENTAL DISPLACEMENT.*/
    /*while(!GetAsyncKeyState(VK_LBUTTON))
    ;*/                                   /*LEFT CLICK TO CONTINUE.*/

    if(fout!=NULL) fprintf(fout,"\"STRESS\"\n");
    for(i=1;i<=nelem;i++)                   /*STRESS OUTPUT,UPDATE.*/
    {
      inputelem(elems,felem,i-1,&elem);

      inputnode(fdisp,elem.node[0]);
      inputnode(fdisp,elem.node[1]);

      elem.sect=(elems+i-1)->sect;             /*READ SECTION DATA.*/

      estress=elemstress(&elem,gvct,felem,fout);

      outputstress(elem,estress,fout);
      free(estress);
    }
    if(wsurf.hwnd!=NULL)
    {
      drawyieldsurface((wsurf.childs+1)->hdcC,
                       (wsurf.childs+1)->vparam,2,4,0,
                       af->fsurface);
      overlayhdc(*(wsurf.childs+1),SRCPAINT);     /*UPDATE DISPLAY.*/
    }
    /*while(!GetAsyncKeyState(VK_LBUTTON))
    ;*/                                   /*LEFT CLICK TO CONTINUE.*/

    if(fout!=NULL) fprintf(fout,"\"REACTION\"\n");
    outputreaction(gmtx,gvct,nodes,confs,freact,fout,nnode);

    updateform(fdisp,gvct,nnode);               /*FORMATION UPDATE.*/
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
                              255,255,255,0);
    }
    overlayhdc(*(wdraw.childs+1),SRCPAINT);       /*UPDATE DISPLAY.*/
  }

  fclose(fin);

  gfree(gmtx,nnode); /*FREE GLOBAL MATRIX.*/
  /*free(gvct);*/
  /*free(confs);*/

  af->eigenvec=(double **)malloc(1*sizeof(double *));
  *((af->eigenvec)+0)=gvct;

  errormessage(" ");
  errormessage("COMPLETED.");
  if(fout!=NULL) fprintf(fout,"COMPLETED.\n");

  memory2=availablephysicalmemory("REMAIN:");
  sprintf(string,"CONSUMPTION:%ld[BYTES]",(memory0-memory2));
  errormessage(string);
  errormessage(" ");

  return 0;
}/*arclm101*/

DWORD availablephysicalmemory(char *comment)
/*RETURN:PHYSICAL MEMORY AVAILABLE [BYTES].*/
{
  MEMORYSTATUS memstat;
  char com[256],str[256];

  GlobalMemoryStatus(&memstat);
  if(comment!=NULL)
  {
    sprintf(str,"%ld/%ld[BYTES] AVAILABLE.",
            memstat.dwAvailPhys,memstat.dwTotalPhys);
    strcpy(com,comment);
    strcat(com,str);
    errormessage(com);
  }

  return memstat.dwAvailPhys;
}/*availablephysicalmemory*/

FILE *fgetstofopen(const char *directory,const char *mode,
                   int dlgitem)
/*INPUT FILENAME AND THEN OPEN FILE.*/
{
  HWND hdwnd=(wmenu.childs+2)->hwnd;
  FILE *f=NULL;
  char fname[256],dandf[256];

  if(hdwnd!=NULL)
  {
    GetDlgItemText(hdwnd,dlgitem,fname,80);
    if(fname[0]=='\0' || fname[0]=='\n') return NULL;

    strcpy(dandf,directory);
    strcat(dandf,fname);
    f=fopen(dandf,mode);

    sprintf(fname,"OPENED=%s",dandf);
    errormessage(fname);
  }
  return f;
}/*fgetstofopen*/

void inputtexttomemory(FILE *ftext,struct arclmframe *af)
/*TRANSLATE ARCLM INPUTFILE TEXT INTO MEMORY.*/
{
  char **data;
  int i,j,ii,jj,k,n;
  long int offset;
  long int scode,hcode,tcode;

  fseek(ftext,0L,SEEK_SET);
  inputinit(ftext,&(af->nnode),&(af->nelem),&(af->nsect));

  /*"INPUTFILE SPECIFICATION"*/
  /*APPELATION*/
  /*NNODE NELEM NSECT*/
  /*ISECT E POI A Ixx Iyy J QxMAX.....MzMIN*/
  /*INODE X Y Z*/
  /*IELEM ISECT NODEI J COORDANGLE BOUNDARY... LONGSTRESS...*/
  /*INODE CONFINEMENTTYPE... CONFINEMENTVALUE...*/
  /*INODE DIRECTION LONGREACTION*/

  for(i=0;i<(af->nsect);i++)
  {
    (af->sects+i)->loff=i;
    readsect(ftext,(af->sects+i));

    (af->sects+i)->dflag=1;
    (af->sects+i)->dcolor.r=255;
    (af->sects+i)->dcolor.g=255;
    (af->sects+i)->dcolor.b=255;
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

int saveasarclm(struct arclmframe *af)
/*SAVE AS ARCLM INPUTFILE.*/
{
  FILE *fin;
  char dandf[256],dir[]=DIRECTORY;
  int i,j,ii,jj,offset;

  /*"INPUTFILE SPECIFICATION"*/
  /*APPELATION*/
  /*NNODE NELEM NSECT*/
  /*ISECT E POI A Ixx Iyy J QxMAX.....MzMIN*/
  /*INODE X Y Z*/
  /*IELEM ISECT NODEI J COORDANGLE BOUNDARY... LONGSTRESS...*/
  /*INODE CONFINEMENTTYPE... CONFINEMENTVALUE...*/
  /*INODE DIRECTION LONGREACTION*/

  strcpy(dandf,dir);
  strcat(dandf,"cansav.inp");
  fin=fopen(dandf,"w");                                /*SAVE FILE.*/
  if(fin==NULL)
  {
    errormessage("ACCESS IMPOSSIBLE.");
    return 0;
  }

  fprintf(fin,"%5d %5d %5d\n",af->nnode,af->nelem,af->nsect);

  for(i=0;i<(af->nsect);i++)
  {
    fprintf(fin,"%5d %.5E %.5f %.5E %.5E %.5E %.5E",
            (af->sects+i)->code,
            (af->sects+i)->E,
            (af->sects+i)->poi,
            (af->sects+i)->area,
            (af->sects+i)->Ixx,
            (af->sects+i)->Iyy,
            (af->sects+i)->Jzz);
    for(j=0;j<=5;j++)
    {
      fprintf(fin," %.1f %.1f",(af->sects+i)->fmax[j],
                               (af->sects+i)->fmin[j]);
    }
    fprintf(fin,"\n");
  }
  for(i=0;i<(af->nnode);i++)
  {
    fprintf(fin,"%5d %.5f %.5f %.5f\n",
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

    for(ii=0;ii<=1;ii++)                 /*BOUNDARY.0:RIGID 1:HINGE*/
    {
      fprintf(fin," %1d %1d %1d",
              (af->elems+i)->iconf[ii][3],
              (af->elems+i)->iconf[ii][4],
              (af->elems+i)->iconf[ii][5]);
    }
    for(ii=0;ii<=1;ii++)                                  /*STRESS.*/
    {
      for(jj=0;jj<=5;jj++)
      {
        fprintf(fin," %.5f",(af->elems+i)->stress[ii][jj]);
      }
    }
    fprintf(fin,"\n");
  }

  for(i=0;i<(af->nnode);i++) /*CONF VECTOR:CONFINEMENT,VALUE.*/
  {
    fprintf(fin,"%5d",(af->nodes+i)->code);

    for(j=0;j<=5;j++)
    {
      offset=6*i+j;
      fprintf(fin," %1d",(af->confs+offset)->iconf);
    }
    for(j=0;j<=5;j++)
    {
      offset=6*i+j;
      fprintf(fin," %.5f",(af->confs+offset)->value);
    }
    fprintf(fin,"\n");
  }
  /*LONGREACTION:UNDER CONSTRUCTION.*/
  fclose(fin);

  return 1;
}/*saveasarclm*/

void initialform(struct onode *nodes,FILE *fdisp,int nnode)
/*INITIAL FORMATION INTO DISPLACEMENT FILE.WITHOUT NODE CODE.*/
{
  int i,ii;
  double zero=0.0;

  for(i=0;i<nnode;i++)
  {
    fwrite(&((nodes+i)->d[0]),sizeof(double),1,fdisp);         /*X.*/
    fwrite(&((nodes+i)->d[1]),sizeof(double),1,fdisp);         /*Y.*/
    fwrite(&((nodes+i)->d[2]),sizeof(double),1,fdisp);         /*Z.*/

    for(ii=1;ii<=3;ii++)
    {fwrite(&zero,sizeof(double),1,fdisp);}          /*THETA X,Y,Z.*/
  }
  fflush(fdisp);

  return;
}/*initialform*/

int initialnode(struct onode *nodes,int nnode,int code,
                struct onode *node)
/*INITIAL POSITION OF NODE.....NOT YET.*/
{
  int i;

  for(i=0;i<nnode;i++)
  {
    if((nodes+i)->code==code) break;
  }
  *node=*(nodes+i);

  return i;
}/*initialnode*/

void initialelem(struct owire *elems,FILE *felem,int nelem)
/*ASSEMBLAGE LONG STRESS.*/
{
  int i,j;

  for(i=0;i<nelem;i++)
  {
    fwrite(&((elems+i)->code),sizeof(long int),1,felem);  /*CODE.*/

    for(j=0;j<=1;j++)                                /*ICONF[2][6].*/
    {
      fwrite(&((elems+i)->iconf[j]),sizeof(signed char),6,felem);
    }
    for(j=0;j<=1;j++)                          /*LONG STRESS[2][6].*/
    {
      fwrite(&((elems+i)->stress[j]),sizeof(double),6,felem);
    }
  }
  fflush(felem);

  return;
}/*initialelem*/

void initialreact(FILE *fin,FILE *freact,int nreact)
/*ASSEMBLAGE LONG REACTIONS.*/
{
  char **data;
  int i,n;
  double ddata;

  for(i=1;i<=nreact;i++)
  {
    data=fgetsbrk(fin,&n);
    if(n==0) ddata=0.0;
    else     ddata=strtod(*(data+2),NULL);    /*0:INODE 1:DIRECTION*/
    fwrite(&ddata,sizeof(double),1,freact);

    for(;n>0;n--) free(*(data+n-1));
    free(data);
  }
  fflush(freact);

  return;
}/*initialreact*/

void inputinit(FILE *fin,int *nnode,int *nelem,int *nsect)
{
  char **data;
  int n;

  data=fgetsbrk(fin,&n);
  *nnode=strtol(*(data+0),NULL,10);
  *nelem=strtol(*(data+1),NULL,10);
  *nsect=strtol(*(data+2),NULL,10);
  for(;n>0;n--) free(*(data+n-1));
  free(data);

  return;
}/*inputinit*/

void inputnode(FILE *fdisp,struct onode *node)
{
  long int loffset;

  loffset=(node->loff)*6*sizeof(double);
  fseek(fdisp,loffset,SEEK_SET);
  fread(&(node->d[0]),sizeof(double),1,fdisp);
  fread(&(node->d[1]),sizeof(double),1,fdisp);
  fread(&(node->d[2]),sizeof(double),1,fdisp);

  return;
}/*inputnode*/

void inputelem(struct owire *elems,FILE *felem,int offset,
               struct owire *elem)
{
  int i;
  long int loffset;

  elem->code=(elems+offset)->code;                  /*ELEMENT CODE.*/
  elem->loff=offset;
  elem->sect=(elems+offset)->sect;               /*SECTION POINTER.*/
  elem->node[0]=(elems+offset)->node[0];            /*HEAD POINTER.*/
  elem->node[1]=(elems+offset)->node[1];            /*TAIL POINTER.*/
  elem->cangle=(elems+offset)->cangle;               /*COORD ANGLE.*/

  if(felem!=NULL)
  {
    loffset=offset
           *(sizeof(long int)
             +12*sizeof(signed char)
             +12*sizeof(double))
           +sizeof(long int);
    fseek(felem,loffset,SEEK_SET);
    for(i=0;i<=1;i++)                                   /*BOUNDARY.*/
    {fread(elem->iconf[i],sizeof(signed char),6,felem);}
    for(i=0;i<=1;i++)                                     /*STRESS.*/
    {fread(elem->stress[i],sizeof(double),6,felem);}
  }
  else
  {
    for(i=0;i<=5;i++)
    {
      elem->iconf[0][i]=(elems+offset)->iconf[0][i];
      elem->iconf[1][i]=(elems+offset)->iconf[1][i];
    }
  }
  return;
}/*inputelem*/

void readsect(FILE *fin,struct osect *sect)
/*READ SECTION FROM INPUT FILE.*/
{
  char **data;
  int i,k,n;

  data=fgetsbrk(fin,&n);
  sect->code=strtol(*(data+0),NULL,10);             /*SECTION CODE.*/
  sect->E=strtod(*(data+1),NULL);                /*YOUNG'S MODULUS.*/
  sect->poi=strtod(*(data+2),NULL);              /*POISSON'S RATIO.*/
  sect->area=strtod(*(data+3),NULL);                        /*AREA.*/
  sect->Ixx=strtod(*(data+4),NULL);
  sect->Iyy=strtod(*(data+5),NULL);
  sect->Jzz=strtod(*(data+6),NULL);          /*ST.VENANT'S TORTION.*/
  k=7;
  for(i=0;i<=5;i++)           /*UPPER,LOWER LIMIT OF YIELD SURFACE.*/
  {
    sect->fmax[i]=strtod(*(data+k),NULL); k++;
    sect->fmin[i]=strtod(*(data+k),NULL); k++;
  }
  for(;n>0;n--) free(*(data+n-1));
  free(data);

  return;
}/*readsect*/

double **directioncosine(double x1,double y1,double z1,
                         double x2,double y2,double z2,
                         double cangle)
{
  int i;
  double dl1,dl2,dx,dy,dz;
  double c0,c1,c2;
  double ol2,ox,oy,oz,al1,xx,xy,xz,yx,yy,yz; /*AXIS VECTOR*/
  double **drccos;

  drccos=(double **)malloc(3*sizeof(double *));
  for(i=0;i<=2;i++) *(drccos+i)=(double *)malloc(3*sizeof(double));
  dx=x2-x1;
  dy=y2-y1;
  dz=z2-z1;
  dl1=sqrt(dx*dx+dy*dy+dz*dz);
  dl2=sqrt(dx*dx+dy*dy);

  if(dl2==0.0) /*PERPENDICULAR COLUMN*/
  {
    c0=0.0; c1=1.0; c2=0.0;

    *(*(drccos+0)+0)=dx/dl1;
    *(*(drccos+0)+1)=dy/dl1;
    *(*(drccos+0)+2)=dz/dl1;
    *(*(drccos+1)+0)=-c2*cos(cangle)-c1**(*(drccos+0)+2)*sin(cangle);
    *(*(drccos+1)+1)= c1*cos(cangle)-c2**(*(drccos+0)+2)*sin(cangle);
    *(*(drccos+1)+2)= c0*sin(cangle);
    *(*(drccos+2)+0)= c2*sin(cangle)-c1**(*(drccos+0)+2)*cos(cangle);
    *(*(drccos+2)+1)=-c1*sin(cangle)-c2**(*(drccos+0)+2)*cos(cangle);
    *(*(drccos+2)+2)= c0*cos(cangle);
  }
  else if((dl2/dl1)<0.1) /*DEFORMED COLUMN*/
  {
    ox=dx/dl1;
    oy=dy/dl1;
    oz=dz/dl1;
    ol2=ox*ox+oy*oy;
    *(*(drccos+0)+0)=dx/dl1; /*AXIS EZ*/
    *(*(drccos+0)+1)=dy/dl1;
    *(*(drccos+0)+2)=dz/dl1;

    xx=(ox*oy*(oz-1.0)*cos(cangle)
        -(ox*ox*oz+oy*oy)*sin(cangle))/ol2; /*AXIS EX*/
    xy=((ox*ox+oy*oy*oz)*cos(cangle)
        -ox*oy*(oz-1.0)*sin(cangle))/ol2;
    xz=ox*sin(cangle)-oy*cos(cangle);
    al1=sqrt(xx*xx+xy*xy+xz*xz);
    *(*(drccos+1)+0)=xx/al1;
    *(*(drccos+1)+1)=xy/al1;
    *(*(drccos+1)+2)=xz/al1;

    yx=oy*xz-oz*xy; /*AXIS EY=EZxEX*/
    yy=oz*xx-ox*xz;
    yz=ox*xy-oy*xx;
    al1=sqrt(yx*yx+yy*yy+yz*yz);
    *(*(drccos+2)+0)=yx/al1;
    *(*(drccos+2)+1)=yy/al1;
    *(*(drccos+2)+2)=yz/al1;
  }
  else /*BEAM*/
  {
    c0=dl2/dl1; c1=dx/dl2; c2=dy/dl2;

    *(*(drccos+0)+0)=dx/dl1;
    *(*(drccos+0)+1)=dy/dl1;
    *(*(drccos+0)+2)=dz/dl1;
    *(*(drccos+1)+0)=-c2*cos(cangle)-c1**(*(drccos+0)+2)*sin(cangle);
    *(*(drccos+1)+1)= c1*cos(cangle)-c2**(*(drccos+0)+2)*sin(cangle);
    *(*(drccos+1)+2)= c0*sin(cangle);
    *(*(drccos+2)+0)= c2*sin(cangle)-c1**(*(drccos+0)+2)*cos(cangle);
    *(*(drccos+2)+1)=-c1*sin(cangle)-c2**(*(drccos+0)+2)*cos(cangle);
    *(*(drccos+2)+2)= c0*cos(cangle);
  }

  return drccos;
}/*directioncosine*/

double **filmdrccos(struct onode n1,struct onode n2,struct onode n3)
/*RETURN:FILM DIRECTION COSINE.*/
{
  int i;
  double dl,Xx,Yx,Zx,Xy,Yy,Zy,Xz,Yz,Zz;
  double **drccos;

  drccos=(double **)malloc(3*sizeof(double *));
  for(i=0;i<=2;i++) *(drccos+i)=(double *)malloc(3*sizeof(double));

  Xx=n2.d[GX]-n1.d[GX]; /*AXIS x=FIRST LINE.*/
  Yx=n2.d[GY]-n1.d[GY];
  Zx=n2.d[GZ]-n1.d[GZ];
  dl=sqrt(Xx*Xx+Yx*Yx+Zx*Zx);
  *(*(drccos+EX)+0)=Xx/dl;
  *(*(drccos+EX)+1)=Yx/dl;
  *(*(drccos+EX)+2)=Zx/dl;

  Xz=Yx*(n3.d[GZ]-n1.d[GZ])-Zx*(n3.d[GY]-n1.d[GY]); /*AXIS z=OUTER.*/
  Yz=Zx*(n3.d[GX]-n1.d[GX])-Xx*(n3.d[GZ]-n1.d[GZ]);
  Zz=Xx*(n3.d[GY]-n1.d[GY])-Yx*(n3.d[GX]-n1.d[GX]);
  dl=sqrt(Xz*Xz+Yz*Yz+Zz*Zz);
  *(*(drccos+EZ)+0)=Xz/dl;
  *(*(drccos+EZ)+1)=Yz/dl;
  *(*(drccos+EZ)+2)=Zz/dl;

  Xy=Yz*Zx-Zz*Yx; /*AXIS y=OUTER.*/
  Yy=Zz*Xx-Xz*Zx;
  Zy=Xz*Yx-Yz*Xx;
  dl=sqrt(Xy*Xy+Yy*Yy+Zy*Zy);
  *(*(drccos+EY)+0)=Xy/dl;
  *(*(drccos+EY)+1)=Yy/dl;
  *(*(drccos+EY)+2)=Zy/dl;

  return drccos;
}/*filmdrccos*/

double **transmatrix(/*struct owire elem,*/double **drccos)
/*ASSEMBLAGE TRANSFORMATION MATRIX.*/
{
  int i,j,n;
  double **t;

  t=(double **)malloc(12*sizeof(double *));
  for(n=1;n<=4;n++)
  {
    for(i=1;i<=3;i++)
    {
      *(t+3*(n-1)+i-1)=(double *)malloc(12*sizeof(double));
      for(j=1;j<=12;j++)
      {*(*(t+3*(n-1)+i-1)+j-1)=0.0;}
      for(j=1;j<=3;j++)
      {*(*(t+3*(n-1)+i-1)+3*(n-1)+j-1)=*(*(drccos+i-1)+j-1);}
    }
  }
  return t;
}/*transmatrix*/

double **assememtx(struct owire elem)
/*ASSEMBLAGE ELASTIC MATRIX.*/
{
  int i,j;
  double **e;
  double dx,dy,dz,dl;
  double E,poi,A,Ixx,Iyy,J,G;

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

  e=(double **)malloc(12*sizeof(double *));
  for(i=1;i<=12;i++)
  {
    *(e+i-1)=(double *)malloc(12*sizeof(double));
    for(j=1;j<=12;j++)
    {
      *(*(e+i-1)+j-1)=0.0;                               /*INITIAL.*/
    }
  }
  *(*(e+0)+0)=E*A/dl; *(*(e+0)+6)=-*(*(e+0)+0);
  *(*(e+1)+1)=12.0*E*Iyy/(dl*dl*dl); *(*(e+1)+5)=6.0*E*Iyy/(dl*dl);
  *(*(e+1)+7)=-*(*(e+1)+1); *(*(e+1)+11)=*(*(e+1)+5);
  *(*(e+2)+2)=12.0*E*Ixx/(dl*dl*dl); *(*(e+2)+4)=-6.0*E*Ixx/(dl*dl);
  *(*(e+2)+8)=-*(*(e+2)+2); *(*(e+2)+10)=*(*(e+2)+4);
  *(*(e+3)+3)=G*J/dl; *(*(e+3)+9)=-*(*(e+3)+3);
  *(*(e+4)+2)=*(*(e+2)+4); *(*(e+4)+4)=4.0*E*Ixx/dl;
  *(*(e+4)+8)=-*(*(e+2)+4); *(*(e+4)+10)=2.0*E*Ixx/dl;
  *(*(e+5)+1)=*(*(e+1)+5); *(*(e+5)+5)=4.0*E*Iyy/dl;
  *(*(e+5)+7)=-*(*(e+1)+5); *(*(e+5)+11)=2.0*E*Iyy/dl;
  *(*(e+6)+0)=*(*(e+0)+6); *(*(e+6)+6)=*(*(e+0)+0);
  *(*(e+7)+1)=*(*(e+1)+7); *(*(e+7)+5)=*(*(e+5)+7);
  *(*(e+7)+7)=*(*(e+1)+1); *(*(e+7)+11)=*(*(e+5)+7);
  *(*(e+8)+2)=*(*(e+2)+8); *(*(e+8)+4)=*(*(e+4)+8);
  *(*(e+8)+8)=*(*(e+2)+2); *(*(e+8)+10)=*(*(e+4)+8);
  *(*(e+9)+3)=*(*(e+3)+9); *(*(e+9)+9)=*(*(e+3)+3);
  *(*(e+10)+2)=*(*(e+2)+10); *(*(e+10)+4)=*(*(e+4)+10);
  *(*(e+10)+8)=*(*(e+8)+10); *(*(e+10)+10)=*(*(e+4)+4);
  *(*(e+11)+1)=*(*(e+1)+11); *(*(e+11)+5)=*(*(e+5)+11);
  *(*(e+11)+7)=*(*(e+7)+11); *(*(e+11)+11)=*(*(e+5)+5);
  return e;
}/*assememtx*/

double **assemgmtx(struct owire elem,double *estress)
/*ASSEMBLAGE GEOMETRICAL MATRIX.*/
{
  int i,j;
  double **g;
  double dx,dy,dz,dl;
  double A,Ixx,Iyy;
  /*double E,poi,J,G;*/
  double gam;
  double N,Mxi,Mxj,Myi,Myj,Qx,Qy;

  dx=elem.node[1]->d[0]-elem.node[0]->d[0];
  dy=elem.node[1]->d[1]-elem.node[0]->d[1];
  dz=elem.node[1]->d[2]-elem.node[0]->d[2];
  dl=sqrt(dx*dx+dy*dy+dz*dz);
  /*E=elem.sect->E;*/
  /*poi=elem.sect->poi;*/
  A=elem.sect->area;
  Ixx=elem.sect->Ixx;
  Iyy=elem.sect->Iyy;
  /*J=elem.sect->Jzz;*/
  /*G=0.5*E/(1.0+poi);*/

  gam=-(Ixx+Iyy)/A;

  N  =-*(estress+ 0); /*+:TENSION*/
  Qx = *(estress+ 1);
  Qy = *(estress+ 2);
  Mxi= *(estress+ 4);
  Myi= *(estress+ 5);
  Mxj= *(estress+10);
  Myj= *(estress+11);

  g=(double **)malloc(12*sizeof(double *));
  for(i=1;i<=12;i++)
  {
    *(g+i-1)=(double *)malloc(12*sizeof(double));
    for(j=1;j<=12;j++)
    {
      *(*(g+i-1)+j-1)=0.0;                               /*INITIAL.*/
    }
  }

  *(*(g+ 1)+ 1)=1.2*N/dl;
  *(*(g+ 1)+ 3)=-Mxi/dl;       *(*(g+ 3)+1)= *(*(g+1)+ 3);
  *(*(g+ 1)+ 5)=0.1*N;         *(*(g+ 5)+1)= *(*(g+1)+ 5);
  *(*(g+ 1)+ 7)=-*(*(g+1)+1);  *(*(g+ 7)+1)= *(*(g+1)+ 7);
  *(*(g+ 1)+ 9)=-Mxj/dl;       *(*(g+ 9)+1)= *(*(g+1)+ 9);
  *(*(g+ 1)+11)= *(*(g+1)+5);  *(*(g+11)+1)= *(*(g+1)+11);

  *(*(g+ 2)+ 2)=*(*(g+1)+1);
  *(*(g+ 2)+ 3)= Myi/dl;       *(*(g+ 3)+2)=*(*(g+2)+ 3);
  *(*(g+ 2)+ 4)=-*(*(g+1)+5);  *(*(g+ 4)+2)=*(*(g+2)+ 4);
  *(*(g+ 2)+ 8)=-*(*(g+1)+1);  *(*(g+ 8)+2)=*(*(g+2)+ 8);
  *(*(g+ 2)+ 9)= Myj/dl;       *(*(g+ 9)+2)=*(*(g+2)+ 9);
  *(*(g+ 2)+10)=-*(*(g+1)+5);  *(*(g+10)+2)=*(*(g+2)+10);

  *(*(g+ 3)+ 3)=N*gam/dl;
  *(*(g+ 3)+ 4)=Qx*dl/6.0;     *(*(g+ 4)+3)=*(*(g+3)+ 4);
  *(*(g+ 3)+ 5)=Qy*dl/6.0;     *(*(g+ 5)+3)=*(*(g+3)+ 5);
  *(*(g+ 3)+ 7)=-*(*(g+1)+3);  *(*(g+ 7)+3)=*(*(g+3)+ 7);
  *(*(g+ 3)+ 8)=-*(*(g+2)+3);  *(*(g+ 8)+3)=*(*(g+3)+ 8);
  *(*(g+ 3)+ 9)=-*(*(g+3)+3);  *(*(g+ 9)+3)=*(*(g+3)+ 9);
  *(*(g+ 3)+10)=-*(*(g+3)+4);  *(*(g+10)+3)=*(*(g+3)+10);
  *(*(g+ 3)+11)=-*(*(g+3)+5);  *(*(g+11)+3)=*(*(g+3)+11);

  *(*(g+ 4)+ 4)=N*dl/7.5;
  *(*(g+ 4)+ 8)= *(*(g+1)+5);  *(*(g+ 8)+4)=*(*(g+4)+ 8);
  *(*(g+ 4)+ 9)=-*(*(g+3)+4);  *(*(g+ 9)+4)=*(*(g+4)+ 9);
  *(*(g+ 4)+10)=-N*dl/30.0;    *(*(g+10)+4)=*(*(g+4)+10);

  *(*(g+ 5)+ 5)= *(*(g+4)+4);
  *(*(g+ 5)+ 7)=-*(*(g+1)+5);  *(*(g+ 7)+5)=*(*(g+5)+ 7);
  *(*(g+ 5)+ 9)=-*(*(g+3)+5);  *(*(g+ 9)+5)=*(*(g+5)+ 9);
  *(*(g+ 5)+11)= *(*(g+4)+10); *(*(g+11)+5)=*(*(g+5)+11);

  *(*(g+ 7)+ 7)= *(*(g+1)+1);
  *(*(g+ 7)+ 9)=-*(*(g+1)+9);  *(*(g+ 9)+7)=*(*(g+7)+ 9);
  *(*(g+ 7)+11)=-*(*(g+1)+5);  *(*(g+11)+7)=*(*(g+7)+11);

  *(*(g+ 8)+ 8)= *(*(g+1)+1);
  *(*(g+ 8)+ 9)=-*(*(g+2)+9);  *(*(g+ 9)+8)=*(*(g+8)+ 9);
  *(*(g+ 8)+10)= *(*(g+1)+5);  *(*(g+10)+8)=*(*(g+8)+10);

  *(*(g+ 9)+ 9)= *(*(g+3)+3);
  *(*(g+ 9)+10)= *(*(g+3)+4);  *(*(g+10)+9)=*(*(g+9)+10);
  *(*(g+ 9)+11)= *(*(g+3)+5);  *(*(g+11)+9)=*(*(g+9)+11);

  *(*(g+10)+10)= *(*(g+4)+4);

  *(*(g+11)+11)= *(*(g+4)+4);
  return g;
}/*assemgmtx*/

double **assempmtx(struct owire elem,double **estiff)
/*ASSEMBLAGE PLASTIC MATRIX.*/
{
  /*char str[256],s[80];*/
  signed char iconf[12];
  int i,j;
  double **p,f[2],q[12][2],a[2][2],det;

  for(i=0;i<=1;i++)                   /*ICONF[2][6] INTO ICONF[12].*/
  {
    for(j=0;j<=5;j++) iconf[6*i+j]=elem.iconf[i][j];
  }

  coefficients(elem,estiff,f,q,a);

  p=(double **)malloc(12*sizeof(double *));
  for(i=1;i<=12;i++)
  {
    *(p+i-1)=(double *)malloc(12*sizeof(double));
    for(j=1;j<=12;j++) *(*(p+i-1)+j-1)=0.0;              /*INITIAL.*/
  }

  if(a[0][0]!=0.0 && a[1][1]==0.0)         /*IF I:PLASTIC J:ELASTIC*/
  {
    for(i=1;i<=12;i++)
    {
      if(iconf[i-1]!=1)
      {
        for(j=1;j<=12;j++)
        {
          if(iconf[j-1]!=1)
          {
            *(*(p+i-1)+j-1)=-1.0/a[0][0]*q[i-1][0]*q[j-1][0];
          }
        }
      }
    }
  }
  if(a[0][0]==0.0 && a[1][1]!=0.0)         /*IF I:ELASTIC J:PLASTIC*/
  {
    for(i=1;i<=12;i++)
    {
      if(iconf[i-1]!=1)
      {
        for(j=1;j<=12;j++)
        {
          if(iconf[j-1]!=1)
          {
            *(*(p+i-1)+j-1)=-1.0/a[1][1]*q[i-1][1]*q[j-1][1];
          }
        }
      }
    }
  }
  if((a[0][0]!=0.0)&&(a[1][1]!=0.0))       /*IF I:PLASTIC J:PLASTIC*/
  {
    for(i=1;i<=12;i++)
    {
      if(iconf[i-1]!=1)
      {
        for(j=1;j<=12;j++)
        {
          if(iconf[j-1]!=1)
          {
            if((a[0][0]*a[1][1]-a[0][1]*a[1][0])==0.0)
            {
              errormessage("ASSEMPMTX:UNDER CONSIDERATION.");
            }
            det=-1.0/(a[0][0]*a[1][1]-a[0][1]*a[1][0]);
            *(*(p+i-1)+j-1)=det*(a[1][1]*q[i-1][0]*q[j-1][0]
                                -a[0][1]*q[i-1][0]*q[j-1][1]
                                -a[1][0]*q[i-1][1]*q[j-1][0]
                                +a[0][0]*q[i-1][1]*q[j-1][1]);
          }
        }
      }
    }
  }

  for(i=1;i<=12;i++)                                    /*ADDITION.*/
  {
    if(iconf[i-1]!=1)
    {
      for(j=1;j<=12;j++)
      {
        if(iconf[j-1]!=1)
        {
          *(*(estiff+i-1)+j-1)+=*(*(p+i-1)+j-1);    /*[k]=[ke]+[kp]*/
        }
      }
    }
  }
  for(i=0;i<=11;i++) free(*(p+i));
  free(p);

  return estiff;
}/*assempmtx*/

void coefficients(struct owire elem,double **estiff,
                  double f[],double q[][2],double a[][2])
/*ASSEMBLAGE PLASTIC COEFFICIENTS.*/
{
  /*char str[256];*/
  signed char iconf[12];
  int i,j,ii,jj;
  double fc[6],fu[6],dfdp[2][6],unit,value;

  for(i=0;i<=1;i++)                   /*ICONF[2][6] INTO ICONF[12].*/
  {
    for(j=0;j<=5;j++) iconf[6*i+j]=elem.iconf[i][j];
  }
  for(i=0;i<=5;i++)                /*CENTER,WIDTH OF YIELD SURFACE.*/
  {
    fc[i]=0.5*(elem.sect->fmax[i]+elem.sect->fmin[i]);     /*CENTER*/
    fu[i]=0.5*(elem.sect->fmax[i]-elem.sect->fmin[i]);      /*WIDTH*/
  }
  for(i=0;i<=1;i++)   /*ASSEMBLAGE YIELD FUNCTION:f,VECTOR:{df/dp}.*/
  {       /*f=POW(|Qx/Qxu|,EXP)+POW(|Qy/Qyu|,EXP)+POW(|Nz/Nzu|,EXP)*/
          /* +POW(|Mx/Mxu|,EXP)+POW(|My/Myu|,EXP)+POW(|Mz/Mzu|,EXP)*/
    f[i]=0.0;
    for(j=0;j<=5;j++)
    {
      if(elem.iconf[i][j]!=1)        /*-1:PLASTIC 0:ELASTIC 1:HINGE*/
      {
        if(j==0 || j==3) /*N,Mz*/
        {
          value=fabs(elem.stress[i][j]-pow(-1.0,i)*fc[j])/fu[j];
          f[i]+=pow(value,EXPONENT);
        }
        else /*Qx,Qy,Mx,My*/
        {
          value=fabs(elem.stress[i][j]-fc[j])/fu[j];
          f[i]+=pow(value,EXPONENT);
        }
        if(j==0 || j==3) unit=elem.stress[i][j]-pow(-1.0,i)*fc[j];
        else             unit=elem.stress[i][j]-fc[j];

        if(unit!=0.0) unit/=fabs(unit);                     /*SIGN.*/
        else          unit=0.0;                  /*ON ILLEGAL LINE.*/
        dfdp[i][j]=unit/fu[j]*pow(value,(EXPONENT-1.0));
      }
    }
  }

  for(i=0;i<=11;i++)                        /*ASSEMBLAGE VECTOR:{q}*/
  {
    if(iconf[i]!=1)
    {
      for(j=0;j<=1;j++)
      {
        q[i][j]=0.0;
        for(jj=0;jj<=5;jj++)
        {
          if(iconf[6*j+jj]!=1) /*if(iconf[6*j+jj]==-1)*/
          {
            q[i][j]+=*(*(estiff+i)+6*j+jj)*dfdp[j][jj];
          }
        }
      }
    }
  }
  for(i=0;i<=1;i++)                    /*INNER PRODUCT a={df/dp}{q}*/
  { for(j=0;j<=1;j++)
    {
      a[i][j]=0.0;
      for(ii=0;ii<=5;ii++)
      {
        if((elem.iconf[i][ii]==-1)&&(elem.iconf[j][ii]==-1))
        {
          a[i][j]+=dfdp[i][ii]*q[6*i+ii][j];
        }
      }
    }
  }
  return;
}/*coefficients*/

double **modifyhinge(struct owire elem,double **estiff)
/*MODIFY ELEMENT STIFFNESS MATRX BY HINGE.*/
{
  int n,i,ii,jj,kk;
  double h[12][12],**e=estiff;

  for(n=1;n<=2;n++)
  {
    for(i=1;i<=6;i++)
    {
      if(elem.iconf[n-1][i-1]==1)
      {
        kk=6*(n-1)+i;
        for(ii=1;ii<=12;ii++)
        {
          for(jj=1;jj<=12;jj++)
          {
            h[ii-1][jj-1]=- *(*(e+ii-1)+kk-1)
                          / *(*(e+kk-1)+kk-1)
                          * *(*(e+kk-1)+jj-1);
          }
        }
        for(ii=1;ii<=12;ii++)
        {
          for(jj=1;jj<=12;jj++) *(*(e+ii-1)+jj-1)+=h[ii-1][jj-1];
        }
      }
    }
  }
  return e;
}/*modifyhinge*/

double **transformation(double **estiff,double **tmatrix)
/*TRANSFORM ELEMENT MATRIX INTO GLOBAL.*/
{
  int i;
  double **e,**t;

  e=matrixmatrix(estiff,tmatrix,12);
  t=matrixtranspose(tmatrix,12);                             /*[Tt]*/

  for(i=0;i<=11;i++) free(*(estiff+i));
  free(estiff);

  estiff=matrixmatrix(t,e,12);                     /*[K]=[Tt][k][T]*/

  for(i=0;i<=11;i++) free(*(e+i));
  free(e);
  for(i=0;i<=11;i++) free(*(t+i));
  free(t);

  return estiff;
}/*transformation*/

void assemgstiffness(struct gcomponent *gmtx,
                     double **estiff,
                     struct owire *elem)
/*ASSEMBLAGE ELEMENT MATRIX INTO GLOBAL MATRIX.*/
{
  long int i,j,ii,jj,n1,n2;
  double edata,gdata;

  for(n1=0;n1<=1;n1++)
  {
    for(i=1;i<=6;i++)
    {
      ii=6*(elem->node[n1]->loff)+i;
      for(n2=0;n2<=1;n2++)
      {
        for(j=1;j<=6;j++)
        {
          jj=6*(elem->node[n2]->loff)+j;
          if(ii>=jj)
          {
            edata=*(*(estiff+6*n1+i-1)+6*n2+j-1);
            if(edata!=0.0)
            {
              gread(gmtx,ii,jj,&gdata);
              gdata+=edata;
              gwrite(gmtx,ii,jj,gdata);
            } /*THESE ARE PUT TOGETHER INTO "gadd".NOT YET.*/
          }
        }
      }
    }
  }
  return;
}/*assemgstiffness*/

void assemconf(struct oconf *confs,double *gvct,double dsafety,
               int nnode)
/*ASSEMBLAGE CONFINEMENT VALUE INTO GLOBAL VECTOR.*/
{
  int i;
  double conf;

  for(i=0;i<6*nnode;i++)      /*LOADS OR DISPS GIVEN INCREMENTALLY.*/
  {
    conf=(confs+i)->value;

    /**(gvct+i-1)+=dsafety*conf;*/   /*"+=":IF GVECTOR INITIALIZED.*/
    *(gvct+i)=dsafety*conf;         /*"=":IF GVECTOR UNINITIALIZED.*/
  }

  return;
}/*assemconf*/

void modifygivend(struct gcomponent *gmtx,double *gvct,
                  struct oconf *confs,int nnode)
/*MODIFY GLOBAL VECTOR BY GIVEN DISPLACEMENT.*/
{
  signed char iconf;
  long int ii,jj;
  double gstiff,disp;

  for(ii=1;ii<=6*nnode;ii++)
  {
    disp=*(gvct+ii-1);
    iconf=(confs+ii-1)->iconf;

    if(iconf==1 && disp!=0.0)
    {
      for(jj=1;jj<=(6*nnode);jj++)
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
}/*modifygivend*/

double *extractdisplacement(struct owire elem,double *gvct)
/*EXTRACT ELEMENT DEFORMATION{dU} FROM GLOBAL VECTOR.*/
{
  long int i,loffset;
  int n;
  double *d;

  d=(double *)malloc(12*sizeof(double));
  if(d==NULL)
  {
    errormessage("EXTRACTDISPLACEMENT:MEMORY INSUFFICIENT.");
    return NULL;
  }
  for(n=1;n<=2;n++)
  {
    for(i=1;i<=6;i++)
    {
      loffset=6*(elem.node[n-1]->loff)+i;
      *(d+6*(n-1)+i-1)=*(gvct+loffset-1);
    }
  }
  return d;
}/*extractdisplacement*/

double *elemstress(struct owire *elem,
                   double *gvct,FILE *felem,FILE *fout)
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
  estiff=assempmtx(*elem,estiff);                   /*[k]=[ke]+[kp]*/

  gdisp=extractdisplacement(*elem,gvct);                     /*{dU}*/
  edisp=matrixvector(tmatrix,gdisp,12);              /*{du}=[T]{dU}*/
  estress=matrixvector(estiff,edisp,12);             /*{df}=[k]{du}*/

  updatestress(felem,fout,edisp,estress,estiff,elem);    /*{f}+{df}*/

  free(gdisp);
  free(edisp);
  for(ii=0;ii<=2;ii++) free(*(drccos+ii));
  free(drccos);
  for(ii=0;ii<=11;ii++) free(*(tmatrix+ii));
  free(tmatrix);
  for(ii=0;ii<=11;ii++) free(*(estiff+ii));
  free(estiff);

  return estress;
}/*elemstress*/

void updatestress(FILE *felem,FILE *fout,
                  double *edisp,double *dstress,double **estiff,
                  struct owire *elem)
/*ELEMENT STRESS UPDATE.*/
{
  char s[80],string[256];
  char iconf[12];
  long int i,ii,j,nn[2];
  long int loffset;
  double lamda[2]={0.0,0.0};
  double fc[6],fu[6],f[2],q[12][2],a[2][2],rate;
  double det,detinverse,function;
  double ysX1,ysY1,ysZ1,ysX2,ysY2,ysZ2;
  struct line line;

  nn[0]=elem->node[0]->code;
  nn[1]=elem->node[1]->code;

  loffset=(elem->loff)
         *(sizeof(long int)+12*sizeof(signed char)+12*sizeof(double))
         +sizeof(long int)+12*sizeof(signed char);
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

  for(i=0;i<=1;i++)                   /*ICONF[2][6] INTO ICONF[12].*/
  {
    for(j=0;j<=5;j++) iconf[6*i+j]=elem->iconf[i][j];
  }
  for(i=0;i<=5;i++)                /*CENTER,WIDTH OF YIELD SURFACE.*/
  {
    fc[i]=0.5*(elem->sect->fmax[i]+elem->sect->fmin[i]);
    fu[i]=0.5*(elem->sect->fmax[i]-elem->sect->fmin[i]);
  }

  coefficients(*elem,estiff,f,q,a);

  if(a[0][0]!=0.0 && a[1][1]==0.0)         /*IF I:PLASTIC J:ELASTIC*/
  {
    for(i=0;i<=11;i++)
    {
      if(iconf[i]!=1)
      {
        lamda[0]+=1.0/a[0][0]*q[i][0]*(*(edisp+i));
      }
    }
  }
  if(a[0][0]==0.0 && a[1][1]!=0.0)         /*IF I:ELASTIC J:PLASTIC*/
  {
    for(i=0;i<=11;i++)
    {
      if(iconf[i]!=1)
      {
        lamda[1]+=1.0/a[1][1]*q[i][1]*(*(edisp+i));
      }
    }
  }
  if(a[0][0]!=0.0 && a[1][1]!=0.0)         /*IF I:PLASTIC J:PLASTIC*/
  {
    for(i=0;i<=11;i++)
    {
      if(iconf[i]!=1)
      {
        detinverse=a[0][0]*a[1][1]-a[0][1]*a[1][0];
        if(detinverse==0.0)
        {
          det=0.0;                           /*UNDER CONSIDERATION.*/
          errormessage("UPDATESTRESS:UNDER CONSIDERATION.");
        }
        else det=1.0/detinverse;
        lamda[0]+=det*( a[1][1]*q[i][0]*(*(edisp+i))
                       -a[0][1]*q[i][1]*(*(edisp+i)));
        lamda[1]+=det*(-a[1][0]*q[i][0]*(*(edisp+i))
                       +a[0][0]*q[i][1]*(*(edisp+i)));
      }
    }
  }

  for(i=0;i<=1;i++)
  {
    if(f[i]>=1.0)                                        /*YIELDED.*/
    {
      sprintf(string,"YIELDED:ELEM%d NODE%ld ",elem->code,nn[i]);
      /*errormessage(string);*/
      strcat(string,"FUNCTION:{");
      function=0.0;
      for(ii=0;ii<=5;ii++)
      {
        if(ii==0 || ii==3) /*N,Mz*/
        {
          rate=(elem->stress[i][ii]-pow(-1.0,i)*fc[ii])/fu[ii];
        }
        else /*Qx,Qy,Mx,My*/
        {
          rate=(elem->stress[i][ii]-fc[ii])/fu[ii];
        }
        function+=pow(fabs(rate),EXPONENT);    /*VALUE OF FUNCTION.*/
        sprintf(s," %8.5f",rate);
        strcat(string,s);
        if(elem->iconf[i][ii]==0)
        {
          elem->iconf[i][ii]=-1;     /*-1:PLASTIC 0:ELASTIC 1:HINGE*/
        }
      }
      sprintf(s,"}=%8.5f",function); /*VALUE OF FUNCTION.*/
      strcat(string,s);
      /*errormessage(string);*/
      if(fout!=NULL) fprintf(fout,"%s\n",string);
    }
  }
  for(i=0;i<=1;i++)
  {
    if(lamda[i]<0.0)                                   /*DISLOADED.*/
    {
      sprintf(string,"DISLOAD:ELEM%d NODE%ld",elem->code,nn[i]);
      errormessage(string);
      if(fout!=NULL) fprintf(fout,"%s\n",string);
      for(ii=0;ii<=5;ii++)
      {
        if(elem->iconf[i][ii]==-1)
        {
          elem->iconf[i][ii]=0;      /*-1:PLASTIC 0:ELASTIC 1:HINGE*/
        }
      }
    }
  }

  loffset=(elem->loff)
         *(sizeof(long int)+12*sizeof(signed char)+12*sizeof(double))
         +sizeof(long int);
  fseek(felem,loffset,SEEK_SET);
  for(i=0;i<=1;i++)
  {
    fwrite(elem->iconf[i],sizeof(signed char),6,felem);
  }
  fflush(felem);

  if(wsurf.hwnd!=NULL /*&& (elem->ielem==1 || elem->ielem==5)*/)
  {
    fseek(arc.fsurface,0L,SEEK_END);
    for(i=0;i<=1;i++)                                   /*HEAD,TAIL*/
    {                               /*0:Nz 1:Qx 2:Qy 3:Mz 4:Mx 5:My*/
      ysX1=(elem->stress[i][2]-*(dstress+6*i+2)-fc[2])/fu[2];
      ysY1=(elem->stress[i][4]-*(dstress+6*i+4)-fc[4])/fu[4];
      ysZ1=(elem->stress[i][0]-*(dstress+6*i+0)-pow(-1.0,i)*fc[0])
           /fu[0];
      ysX2=(elem->stress[i][2]-fc[2])/fu[2];
      ysY2=(elem->stress[i][4]-fc[4])/fu[4];
      ysZ2=(elem->stress[i][0]-pow(-1.0,i)*fc[0])/fu[0];

      line.code=elem->code;
      line.ends[0].d[0]=ysX1;
      line.ends[0].d[1]=ysY1;
      line.ends[0].d[2]=ysZ1;
      line.ends[1].d[0]=ysX2;
      line.ends[1].d[1]=ysY2;
      line.ends[1].d[2]=ysZ2;

      if(f[i]>=1.0) /*YIELDED*/
      {
        line.r=255; line.g=100; line.b=255;
      }
      else
      {
        line.r=100; line.g=200; line.b=200;
      }
      fwrite(&line,sizeof(struct line),1,arc.fsurface);
    }
  }

  return;
}/*updatestress*/

void outputdisp(double *gvct,FILE *fout,int nnode,
                struct onode *nodes)
/*OUTPUT NODE DISPLACEMENT.*/
{
  char string[256];
  int i,j;
  double data[6];

  for(i=1;i<=nnode;i++)
  {
    for(j=0;j<=5;j++) data[j]=*(gvct+6*(i-1)+j);
    sprintf(string,
            "NODE:%5ld {dU}= %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f",
            (nodes+i-1)->code,
            data[0],data[1],data[2],data[3],data[4],data[5]);
    if(fout!=NULL) fprintf(fout,"%s\n",string);
  }
}/*outputdisp*/

void updateform(FILE *fdisp,double *gvct,int nnode)
/*FORMATION UPDATE.*/
{
  int i,j;
  long int loffset;
  double data,ddata;

  for(i=1;i<=nnode;i++)
  {
    for(j=1;j<=6;j++)
    {
      loffset=6*(i-1)+j;
      ddata=*(gvct+loffset-1);

      vread(fdisp,loffset,&data);
      data+=ddata;                                       /*{U}+{dU}*/
      vwrite(fdisp,loffset,&data);
    }
  }
  return;
}/*updateform*/

void copyform(struct arclmframe *af,double *gvct)
/*FORMATION COPY.*/
{
  int i,j;
  long int nnode,loff;
  double ddata;

  nnode=af->nnode;

  for(i=0;i<nnode;i++)
  {
    for(j=0;j<3;j++)
    {
      loff=6*i+j;
      ddata=*(gvct+loff);

      (af->nodes+i)->d[j]=(af->ninit+i)->d[j]+ddata;
    }
  }
  return;
}/*copyform*/

void outputstress(struct owire elem,
                  double *estress,FILE *fout)
/*ELEMENT STRESS OUTPUT.*/
{
  char s[80],string[256];
  long int i,n,nn[2];

  for(n=1;n<=2;n++)
  {
    nn[n-1]=elem.node[n-1]->code;
    sprintf(string,"ELEM:%5d NODE%5ld:",elem.code,nn[n-1]);
    for(i=1;i<=6;i++)
    {
      sprintf(s," %11.2E/%11.2E",
              *(estress+6*(n-1)+i-1),elem.stress[n-1][i-1]);
      strcat(string,s);
    }
    if(fout!=NULL) fprintf(fout,"%s\n",string);
  }
  return;
}/*outputstress*/

void outputreaction(struct gcomponent *gmtx,
                    double *gvct,
                    struct onode *nodes,
                    struct oconf *confs,
                    FILE *freact,FILE *fout,int nnode)
/*REACTIONS UPDATE,OUTPUT.*/
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

      sprintf(string,"NODE:%5ld %ld %14.5f",
              (nodes+offset)->code,(i-1)%6+1,dreaction);
      if(fout!=NULL) fprintf(fout,"%s\n",string);
    }
  }
  return;
}/*outputreaction*/

void gethoganparam(HWND hdwnd,struct hoganparam *hp)
/*GET HOGANSHI PARAMETERS FROM DIALOG.*/
{
  char data[80];

  GetDlgItemText(hdwnd,IDH_ASIZE,data,80);
  hp->arrowsize=strtod(data,NULL);
  GetDlgItemText(hdwnd,IDH_NSIZE,data,80);
  hp->notchsize=strtod(data,NULL);

  GetDlgItemText(hdwnd,IDH_TITLEX,(hp->title[0]),80);
  GetDlgItemText(hdwnd,IDH_TITLEY,(hp->title[1]),80);

  GetDlgItemText(hdwnd,IDH_UNITX,(hp->unit[0]),80);
  GetDlgItemText(hdwnd,IDH_UNITY,(hp->unit[1]),80);

  GetDlgItemText(hdwnd,IDH_SCALEX,data,80);
  hp->scale[0]=strtod(data,NULL);
  GetDlgItemText(hdwnd,IDH_SCALEY,data,80);
  hp->scale[1]=strtod(data,NULL);

  GetDlgItemText(hdwnd,IDH_XMAX,data,80);
  hp->max[0]=strtod(data,NULL);
  GetDlgItemText(hdwnd,IDH_YMAX,data,80);
  hp->max[1]=strtod(data,NULL);

  GetDlgItemText(hdwnd,IDH_XMIN,data,80);
  hp->min[0]=strtod(data,NULL);
  GetDlgItemText(hdwnd,IDH_YMIN,data,80);
  hp->min[1]=strtod(data,NULL);

  GetDlgItemText(hdwnd,IDH_PITCHX,data,80);
  hp->pitch[0]=strtod(data,NULL);
  GetDlgItemText(hdwnd,IDH_PITCHY,data,80);
  hp->pitch[1]=strtod(data,NULL);

  return;
}/*gethoganparam*/

void drawhoganaxis(HDC hdc,
                   struct viewparam vp,
                   struct hoganparam hp,
                   int r,int g,int b)
/*DRAW GLOBAL AXIS BY AXONOMETRICS,PERSPECTIVE.*/
{
  HPEN hpen,ppen;
  char str[80];
  int i,nnotch;
  double pos;
  struct onode on,hn,tn;

  SetTextColor(hdc,RGB(255,255,255));                 /*TEXT COLOR.*/
  hpen=CreatePen(PS_SOLID,1,RGB(r,g,b));             /*ARROW COLOR.*/
  ppen=SelectObject(hdc,hpen);

  on=setcoord(0.0,0.0,0.0);
  drawglobaltext(hdc,vp,on,"O");                           /*ORIGIN*/

  hn=setcoord((hp.scale[0])*(hp.min[0]-hp.pitch[0]),0.0,0.0);
  tn=setcoord((hp.scale[0])*(hp.max[0]+2.0*(hp.pitch[0])),0.0,0.0);
  drawglobaltext(hdc,vp,tn,hp.title[0]);                   /*TEXT X*/
  drawglobalarrow(hdc,vp,hn,tn,hp.arrowsize);             /*ARROW X*/

  hn=setcoord(0.0,(hp.scale[1])*(hp.min[1]-hp.pitch[1]),0.0);
  tn=setcoord(0.0,(hp.scale[1])*(hp.max[1]+2.0*(hp.pitch[1])),0.0);
  drawglobaltextaligned(hdc,vp,tn,hp.title[1],
                        TEXT_RIGHT,TEXT_BOTTOM);           /*TEXT Y*/
  drawglobalarrow(hdc,vp,hn,tn,hp.arrowsize);             /*ARROW Y*/

  nnotch=(int)(hp.max[0]/hp.pitch[0]);
  for(i=1;i<=nnotch;i++) /*NOTCH OF AXIS X POSITIVE.*/
  {
    pos=i*(hp.pitch[0]);
    sprintf(str,"%.3f",pos);

    hn=setcoord((hp.scale[0])*pos,-hp.notchsize,0.0);
    tn=setcoord((hp.scale[0])*pos, hp.notchsize,0.0);

    drawglobaltext(hdc,vp,hn,str);
    drawgloballine(hdc,vp,hn,tn);
  }
  nnotch=(int)(hp.min[0]/hp.pitch[0]);
  for(i=-1;i>=nnotch;i--) /*NOTCH OF AXIS X NEGATIVE.*/
  {
    pos=i*(hp.pitch[0]);
    sprintf(str,"%.3f",pos);

    hn=setcoord((hp.scale[0])*pos,-hp.notchsize,0.0);
    tn=setcoord((hp.scale[0])*pos, hp.notchsize,0.0);

    drawglobaltext(hdc,vp,hn,str);
    drawgloballine(hdc,vp,hn,tn);
  }

  nnotch=(int)(hp.max[1]/hp.pitch[1]);
  for(i=1;i<=nnotch;i++) /*NOTCH OF AXIS Y POSITIVE.*/
  {
    pos=i*(hp.pitch[1]);
    sprintf(str,"%.3f",pos);

    hn=setcoord(-hp.notchsize,(hp.scale[1])*pos,0.0);
    tn=setcoord( hp.notchsize,(hp.scale[1])*pos,0.0);

    drawglobaltextaligned(hdc,vp,hn,str,TEXT_RIGHT,TEXT_CENTER);
    drawgloballine(hdc,vp,hn,tn);
  }
  nnotch=(int)(hp.min[1]/hp.pitch[1]);
  for(i=-1;i>=nnotch;i--) /*NOTCH OF AXIS Y NEGATIVE.*/
  {
    pos=i*(hp.pitch[1]);
    sprintf(str,"%.3f",pos);

    hn=setcoord(-hp.notchsize,(hp.scale[1])*pos,0.0);
    tn=setcoord( hp.notchsize,(hp.scale[1])*pos,0.0);

    drawglobaltextaligned(hdc,vp,hn,str,TEXT_RIGHT,TEXT_CENTER);
    drawgloballine(hdc,vp,hn,tn);
  }

  SelectObject(hdc,ppen);
  DeleteObject(hpen);

  return;
}/*drawhoganaxis*/

int drawhoganlines(HDC hdc,HWND hdwnd,struct viewparam vp,
                   FILE *fin)
/*DRAW 2D LINES.*/
{
  /*SAMPLE=C:\CDOCS\HOGAN\HOGAN01.INP*/
  /*SPECIFICATION:6 COORDINATIONS IN EACH LINES.*/
  /*                           SKIP BLANK LINES.*/
  /*Xhead Yhead Zhead Xtail Ytail Ztail*/

  char **data;
  int n;
  struct hoganparam hp;
  struct line l;

  if(hdc==NULL || fin==NULL) return 0;
  fseek(fin,0L,SEEK_SET);

  gethoganparam(hdwnd,&hp);

  drawhoganaxis(hdc,vp,hp,0,0,255);

  while(1)
  {
    data=fgetsbrk(fin,&n);
    if(n==0) break;

    if(n>=6)
    {
      l.ends[0].d[0]=hp.scale[0]*strtod(*(data+0),NULL);
      l.ends[0].d[1]=hp.scale[1]*strtod(*(data+1),NULL);
      l.ends[0].d[2]=0.0;
      l.ends[1].d[0]=hp.scale[0]*strtod(*(data+3),NULL);
      l.ends[1].d[1]=hp.scale[1]*strtod(*(data+4),NULL);
      l.ends[1].d[2]=0.0;

      l.r=0; l.g=255; l.b=255; /*LINE COLOR.*/

      drawlinestruct(hdc,vp,l);
    }

    for(;n>0;n--) free(*(data+n-1));
    free(data);
  }

  return 1;
}/*drawhoganlines*/

struct oconf *addconf(struct oconf conf[6],
                      long int nodeoffset,
                      struct organ *orgbefore)
/*ADD CONFINEMENT OF FRAME.BEFORE ADDITION OF NODE.*/
/*NNODE NOT UPDATE.*/
{
  long int i,j,i1,i2;
  long int loff,nnode;

  nnode=orgbefore->nnode;

  orgbefore->confs=(struct oconf *)
                   realloc(orgbefore->confs,
                           6*(nnode+1)*sizeof(struct oconf));
  if(orgbefore->confs==NULL)
  {
    MessageBox(NULL,"Buffer Null.","AddConf",MB_OK);
    return NULL;
  }

  i1=nnode;
  i2=nodeoffset+1;
  loff=6*(nnode+1)-1;
  for(i=i1;i>=i2;i--)
  {
    for(j=0;j<6;j++)
    {
      *(orgbefore->confs+loff)=*(orgbefore->confs+loff-6);
      loff--;
    }
  }

  loff=6*nodeoffset;
  for(j=0;j<6;j++) *(orgbefore->confs+loff+j)=conf[j];

  return orgbefore->confs;
}/*addconf*/

struct oconf *changeconf(struct oconf conf[6],
                         long int nodeoffset,
                         struct oconf *cinit)
/*CHANGE CONFINEMENT OF FRAME.*/
{
  long int j,loff;

  loff=6*nodeoffset;
  for(j=0;j<6;j++)
  {
    *(cinit+loff)=conf[j];
    loff++;
  }

  return cinit;
}/*changeconf*/

struct oconf *deleteconf(long int nodeoffset,
                         struct organ *orgbefore)
/*DELETE CONFINEMENT OF FRAME.BEFORE DELETING NODE.*/
{
  int i,i1,i2;
  long int nnode;
  /*struct oconf *confs;*/

  nnode=orgbefore->nnode;

  i1=6*nodeoffset;
  i2=6*(nnode-1)-1;
  for(i=i1;i<=i2;i++)
  {
    *(orgbefore->confs+i)=*(orgbefore->confs+i+6);
  }

  orgbefore->confs=(struct oconf *)
                   realloc(orgbefore->confs,
                           6*(nnode-1)*sizeof(struct oconf));
  if(orgbefore->confs==NULL)
  {
    MessageBox(NULL,"Buffer Null.","DeleteConf",MB_OK);
  }
  return orgbefore->confs;
}/*deleteconf*/

struct onode *addnode(struct onode node,
                      struct organ *orgbefore)
/*ADD NODE OF ORGAN.*/
/*RETURN:ADDED NODE POINTER TO NODES OF ORGAN.*/
{
  /*char non[80],none[256];*/
  int i,j,jj,k;
  long int nnode,nelem;
  struct onode *nodes;
  struct oconf c[6];

  nnode=orgbefore->nnode;
  nelem=orgbefore->nelem;

  for(i=0;i<6;i++)
  {
    c[i].iconf=0;
    c[i].value=0.0;
  }

  if(addconf(c,(node.loff),orgbefore)==NULL) return NULL;

  nnode++;
  nodes=(struct onode *)malloc(nnode*sizeof(struct onode));
  if(nodes==NULL)
  {
    MessageBox(NULL,"Buffer Null.","AddNode",MB_OK);
    return NULL;
  }

  orgbefore->nnode=nnode;

  for(i=(nnode-2);i>=node.loff;i--)
  {
    (nodes+i+1)->code=(orgbefore->nodes+i)->code;
    (nodes+i+1)->loff=((orgbefore->nodes+i)->loff)+1;
    (nodes+i+1)->d[0]=(orgbefore->nodes+i)->d[0];
    (nodes+i+1)->d[1]=(orgbefore->nodes+i)->d[1];
    (nodes+i+1)->d[2]=(orgbefore->nodes+i)->d[2];

    for(j=0;j<nelem;j++)
    {
/*sprintf(none,"\0");*/
      for(jj=0;jj<((orgbefore->elems+j)->nnod);jj++)
      {
        if(*((orgbefore->elems+j)->nods+jj)==(orgbefore->nodes+i))
        {
          *((orgbefore->elems+j)->nods+jj)=(nodes+i+1);
        }
/*sprintf(non," %ld",
        *((orgbefore->elems+j)->nods+jj)->code);
strcat(none,non);*/
      }
      for(jj=0;jj<((orgbefore->elems+j)->nban);jj++)
      {
        for(k=0;k<(((orgbefore->elems+j)->bans+jj)->nnod);k++)
        {
          if(*(((orgbefore->elems+j)->bans+jj)->nods+k)
             ==(orgbefore->nodes+i))
          {
            *(((orgbefore->elems+j)->bans+jj)->nods+k)=(nodes+i+1);
          }
        }
      }
/*MessageBox(NULL,none,"AddNode",MB_OK);*/
    }

    for(jj=0;jj<(gelem.nnod-1);jj++)
    {
      if(*(gelem.nods+jj)==(orgbefore->nodes+i))
      {
        *(gelem.nods+jj)=(nodes+i+1);
      }
    }
    for(jj=0;jj<(gelem.nban);jj++)
    {
      for(k=0;k<((gelem.bans+jj)->nnod-1);k++)
      {
        if(*((gelem.bans+jj)->nods+k)==(orgbefore->nodes+i))
        {
          *((gelem.bans+jj)->nods+k)=(nodes+i+1);
        }
      }
    }
  }
  for(i=(node.loff-1);i>=0;i--)
  {
    (nodes+i)->code=(orgbefore->nodes+i)->code;
    (nodes+i)->loff=(orgbefore->nodes+i)->loff;
    (nodes+i)->d[0]=(orgbefore->nodes+i)->d[0];
    (nodes+i)->d[1]=(orgbefore->nodes+i)->d[1];
    (nodes+i)->d[2]=(orgbefore->nodes+i)->d[2];

    for(j=0;j<nelem;j++)
    {
/*sprintf(none,"\0");*/
      for(jj=0;jj<((orgbefore->elems+j)->nnod);jj++)
      {
        if(*((orgbefore->elems+j)->nods+jj)==(orgbefore->nodes+i))
        {
          *((orgbefore->elems+j)->nods+jj)=(nodes+i);
        }
/*sprintf(non," %ld",
        *((orgbefore->elems+j)->nods+jj)->code);
strcat(none,non);*/
      }
      for(jj=0;jj<((orgbefore->elems+j)->nban);jj++)
      {
        for(k=0;k<(((orgbefore->elems+j)->bans+jj)->nnod);k++)
        {
          if(*(((orgbefore->elems+j)->bans+jj)->nods+k)
             ==(orgbefore->nodes+i))
          {
            *(((orgbefore->elems+j)->bans+jj)->nods+k)=(nodes+i);
          }
        }
      }
/*MessageBox(NULL,none,"AddNode",MB_OK);*/
    }

    for(jj=0;jj<(gelem.nnod-1);jj++)
    {
      if(*(gelem.nods+jj)==(orgbefore->nodes+i))
      {
        *(gelem.nods+jj)=(nodes+i);
      }
    }
    for(jj=0;jj<(gelem.nban);jj++)
    {
      for(k=0;k<((gelem.bans+jj)->nnod-1);k++)
      {
        if(*((gelem.bans+jj)->nods+k)==(orgbefore->nodes+i))
        {
          *((gelem.bans+jj)->nods+k)=(nodes+i);
        }
      }
    }
  }

  (nodes+(node.loff))->code=node.code;
  (nodes+(node.loff))->loff=node.loff;
  (nodes+(node.loff))->d[0]=node.d[0];
  (nodes+(node.loff))->d[1]=node.d[1];
  (nodes+(node.loff))->d[2]=node.d[2];

  free(orgbefore->nodes);
  orgbefore->nodes=nodes;

  return (nodes+(node.loff));
}/*addnode*/

struct onode *deletenode(long int nodeoffset,
                         struct organ *orgbefore)
/*DELETE NODE OF ORGAN.*/
/*RETURN:NODES OF ORGAN.*/
{
  int i,j,jj;
  long int nnode,nelem;
  struct onode *nodes;

  nnode=orgbefore->nnode;
  nelem=orgbefore->nelem;

  deleteconf(nodeoffset,orgbefore);

  nnode--;
  nodes=(struct onode *)malloc(nnode*sizeof(struct onode));
  if(nodes==NULL)
  {
    MessageBox(NULL,"Buffer Null.","DeleteNode",MB_OK);
    return orgbefore->nodes;
  }
  orgbefore->nnode=nnode;

  for(i=0;i<=nodeoffset-1;i++)
  {
    *(nodes+i)=*(orgbefore->nodes+i);

    for(j=0;j<nelem;j++)
    {
      for(jj=0;jj<((orgbefore->elems+j)->nnod);jj++)
      {
        if(*((orgbefore->elems+j)->nods+jj)==(orgbefore->nodes+i))
        {
          *((orgbefore->elems+j)->nods+jj)=(nodes+i);
        }
      }
    }
  }
  for(i=nodeoffset;i<=nnode;i++)
  {
    *(nodes+i)=*(orgbefore->nodes+i+1);
    ((nodes+i)->loff)--;

    for(j=0;j<nelem;j++)
    {
      for(jj=0;jj<((orgbefore->elems+j)->nnod);jj++)
      {
        if(*((orgbefore->elems+j)->nods+jj)==(orgbefore->nodes+i+1))
        {
          *((orgbefore->elems+j)->nods+jj)=(nodes+i);
        }
      }
    }
  }

  for(j=0;j<nelem;j++) /*DELETE JOINTED ELEMENT.*/
  {
    for(jj=0;jj<((orgbefore->elems+j)->nnod);jj++)
    {
      if(*((orgbefore->elems+j)->nods+jj)
         ==(orgbefore->nodes+nodeoffset))
      {
        deleteelement(j,orgbefore);
        break;
      }
    }
  }

  free(orgbefore->nodes);
  orgbefore->nodes=nodes;

  return nodes;
}/*deletenode*/

struct oelem *addelement(struct oelem *elem,
                         struct organ *orgbefore)
/*ADD ELEMENT OF ORGAN.OFFSET KNOWN,NODES ALREADY ADDED.*/
{
  /*int i;*/
  long int i,nelem;
  struct oelem *elems;

  nelem=orgbefore->nelem;

  nelem++;
  elems=(struct oelem *)malloc(nelem*sizeof(struct oelem));
  if(elems==NULL)
  {
    MessageBox(NULL,"Buffer Null.","AddElement",MB_OK);
    return orgbefore->elems;
  }

  for(i=0;i<(elem->loff);i++)
  {
    if(!copyelement((orgbefore->elems+i),(elems+i),orgbefore))
    {
      return orgbefore->elems;
    }
  }
  for(i=(elem->loff);i<(nelem-1);i++)
  {
    if(!copyelement((orgbefore->elems+i),(elems+i+1),orgbefore))
    {
      return orgbefore->elems;
    }

    ((elems+i+1)->loff)++;
  }

  if(!copyelement(elem,(elems+(elem->loff)),orgbefore))
  {
    return orgbefore->elems;
  }

  /*for(i=0;i<(orgbefore->nelem);i++)
  {
    freeelement((orgbefore->elems+i),0,0,1,0);
  }*/
  free(orgbefore->elems);

  orgbefore->nelem=nelem;
  orgbefore->elems=elems;

/*MessageBox(NULL,"Completed.","AddElement",MB_OK);*/

  return elems;
}/*addelement*/

struct oelem *deleteelement(long int elemoffset,
                            struct organ *orgbefore)
/*DELETE ELEMENT OF ORGAN.*/
{
  int i,j,jj,count;
  long int nelem;
  struct oelem *elems;

  nelem=orgbefore->nelem;

  nelem--;
  elems=(struct oelem *)malloc(nelem*sizeof(struct oelem));
  if(elems==NULL)
  {
    MessageBox(NULL,"Buffer Null.","DeleteElem",MB_OK);
    return orgbefore->elems;
  }
  orgbefore->nelem=nelem;

  for(i=0;i<=elemoffset-1;i++)
  {
    *(elems+i)=*(orgbefore->elems+i); /* CopyElement ? */
  }
  for(i=elemoffset;i<nelem;i++)
  {
    *(elems+i)=*(orgbefore->elems+i+1); /* CopyElement ? */
    ((elems+i)->loff)--;
  }

  for(i=0;i<((orgbefore->elems+elemoffset)->nnod);i++)
  {
    count=0;

    for(j=0;j<nelem;j++)
    {
      for(jj=0;jj<((elems+j)->nnod);jj++)
      {
        if(*((elems+j)->nods+jj)
           ==*((orgbefore->elems+elemoffset)->nods+i))
        {
          count++;
        }
      }
    }

    if(count==0)
    {
      deletenode((*((orgbefore->elems+elemoffset)->nods+i))->loff,
                 orgbefore);
    }
  }

  for(j=0;j<=nelem;j++) /*FREE EINIT.*/
  {
    if(j==elemoffset) freeelement(orgbefore->elems+j,1,1,1,1);
    /*else              freeelement(orgbefore->elems+j,0,0,1,0);*/
  }
  free(orgbefore->elems);

  orgbefore->elems=elems;

  return elems;
}/*deleteelement*/

void freeelement(struct oelem *elem,
                 char fnods,char fbonds,char fbans,char fbannods)
/*FREE BUFFER OF ELEMENT IF FLAG TRUE.*/
{
  int i;

  if((elem->nban)>0)
  {
    if(fbannods)
    {
      for(i=0;i<(elem->nban);i++)
      {
        free((elem->bans+i)->nods);
      }
    }
    if(fbans) free(elem->bans);
  }
  if(fnods)  free(elem->nods);
  if(fbonds) free(elem->bonds);

  return;
}/*freeelement*/

int copyelement(struct oelem *efrom,struct oelem *eto,
                struct organ *org)
{
  int i,j,k;

  eto->code  =efrom->code;
  eto->loff  =efrom->loff;
  eto->role  =efrom->role;
  eto->type  =efrom->type;
  eto->nnod  =efrom->nnod;
  eto->nban  =efrom->nban;
  eto->cangle=efrom->cangle;

  for(i=0;i<2;i++)
  {
    for(j=0;j<6;j++) eto->initial[i][j]=efrom->initial[i][j];
  }
  for(i=0;i<2;i++)
  {
    for(j=0;j<6;j++) eto->current[i][j]=efrom->current[i][j];
  }

  eto->nods=efrom->nods;
  eto->bonds=efrom->bonds;

  /*eto->nods=(struct onode **)malloc(eto->nnod
                                    *sizeof(struct onode *));
  eto->bonds=(signed char *)malloc(6*(eto->nnod)
                                   *sizeof(signed char));
  if(eto->nods==NULL || eto->bonds==NULL)
  {
    MessageBox(NULL,"Buffer Null.","CopyElement",MB_OK);
    return 0;
  }

  for(i=0;i<eto->nnod;i++)
  {
    *(eto->nods+i)=*(efrom->nods+i);
  }
  for(i=0;i<6*(eto->nnod);i++)
  {
    *(eto->bonds+i)=*(efrom->bonds+i);
  }*/

  if((eto->nban)>0)
  {
    eto->bans=(struct obans *)malloc((eto->nban)
                                     *sizeof(struct obans));
    if(eto->bans==NULL)
    {
      MessageBox(NULL,"Buffer Null.","CopyElement",MB_OK);
      return 0;
    }

    for(i=0;i<(eto->nban);i++)
    {
      /**(eto->bans+i)=*(efrom->bans+i);*/

      (eto->bans+i)->code=(efrom->bans+i)->code;
      (eto->bans+i)->loff=(efrom->bans+i)->loff;
      (eto->bans+i)->nnod=(efrom->bans+i)->nnod;

      (eto->bans+i)->nods=(efrom->bans+i)->nods;
      /*(eto->bans+i)->nods=(struct onode **)
                          malloc(((eto->bans+i)->nnod)
                                 *sizeof(struct onode *));
      if((eto->bans+i)->nods==NULL)
      {
        MessageBox(NULL,"Buffer Null.","CopyElement",MB_OK);
        return 0;
      }

      for(j=0;j<((efrom->bans+i)->nnod);j++)
      {
        *((eto->bans+i)->nods+j)=*((efrom->bans+i)->nods+j);
      }*/
      /*for(j=0;j<((efrom->bans+i)->nnod);j++)
      {
        for(k=0;k<(efrom->nnod);k++)
        {
          if(*((efrom->bans+i)->nods+j)==*(efrom->nods+k))
          {
            *((eto->bans+i)->nods+j)=*(eto->nods+k);
          }
        }
      }*/
    }
  }

  /*eto->sect=org->sects+(efrom->sect->loff);*/
  eto->sect=efrom->sect;

  eto->color.code=efrom->color.code;
  eto->color.axis=efrom->color.axis;
  eto->color.line=efrom->color.line;

  return 1;
}/*copyelement*/

struct onode *createnodeonplane(struct viewparam vp,
                                double mx,double my,
                                struct plane pl,
                                struct onode *ncross)
/*CREATE NODE ON PLANE.PLANE IS ALREADY KNOWN.*/
{
  int i;
  double **drccos,**tdrccos;
  struct plane proj;
  struct line vline;
  struct onode nview,nmouse,norigin;

  /*SCREEN ORIGIN TO GLOBAL.*/
  norigin.d[GX]=(vp.r-vp.odv)/vp.r*vp.ov[GX]+vp.focus.d[GX];
  norigin.d[GY]=(vp.r-vp.odv)/vp.r*vp.ov[GY]+vp.focus.d[GY];
  norigin.d[GZ]=(vp.r-vp.odv)/vp.r*vp.ov[GZ]+vp.focus.d[GZ];

  /*DEFINITION PROJECTION PLANE.*/
  proj.nvec.dc[GX]=vp.ov[GX];
  proj.nvec.dc[GY]=vp.ov[GY];
  proj.nvec.dc[GZ]=vp.ov[GZ];

  proj.nods[0].d[GX]=norigin.d[GX];
  proj.nods[0].d[GY]=norigin.d[GY];
  proj.nods[0].d[GZ]=norigin.d[GZ];
  proj.nods[1].d[GX]=norigin.d[GX]+vp.axono[1][0];    /*AXONO[1]=y.*/
  proj.nods[1].d[GY]=norigin.d[GY]+vp.axono[1][1];
  proj.nods[1].d[GZ]=norigin.d[GZ]+vp.axono[1][2];
  proj.nods[2].d[GX]=norigin.d[GX]+vp.axono[2][0];    /*AXONO[2]=z.*/
  proj.nods[2].d[GY]=norigin.d[GY]+vp.axono[2][1];
  proj.nods[2].d[GZ]=norigin.d[GZ]+vp.axono[2][2];

  /*MOUSE POINT TO GLOBAL.*/
  nmouse.d[EX]=mx;
  nmouse.d[EY]=my;
  nmouse.d[EZ]=0.0;

  drccos=filmdrccos(proj.nods[0],
                    proj.nods[1],
                    proj.nods[2]);              /*DIRECTION COSINE.*/
  tdrccos=matrixtranspose(drccos,3);                  /*MATRIX [T].*/
  nmouse=transltog(tdrccos,proj.nods[0],nmouse);

  for(i=0;i<3;i++) free(*(drccos+i));
  free(drccos);
  for(i=0;i<3;i++) free(*(tdrccos+i));
  free(tdrccos);

  /*VIEWPOINT.*/
  if(vp.type==PERSPECTIVE)
  {
    nview.d[GX]=vp.ov[GX]+vp.focus.d[GX];
    nview.d[GY]=vp.ov[GY]+vp.focus.d[GY];
    nview.d[GZ]=vp.ov[GZ]+vp.focus.d[GZ];
  }
  if(vp.type==AXONOMETRIC)
  {
    nview.d[GX]=nmouse.d[GX]+proj.nvec.dc[GX];
    nview.d[GY]=nmouse.d[GY]+proj.nvec.dc[GY];
    nview.d[GZ]=nmouse.d[GZ]+proj.nvec.dc[GZ];
  }

  /*INTERSECTION.*/
  vline.ends[0]=nview;
  vline.ends[1]=nmouse;

  intersectlineplane(vline,pl,ncross);

  return ncross;
}/*createnodeonplane*/

struct onode *createorgannode(long int code,
                              struct viewparam vp,
                              long int mx,long int my,
                              struct organ *org)
/*CREATE NODE OF ORGAN.*/
{
  /*char non[256];*/
  int found;
  int Nx,Ny;
  long int i,j,k,loff;
  long int codemax;
  double dmx,dmy;
  struct onode nmouse; /*MOUSE POSITION AS NODE.*/
  struct onode *bnods,nadd;
  struct obans ban,*cban;
  struct plane pl;

  nmouse.d[GX]=(double)mx;
  nmouse.d[GY]=(double)my;
  nmouse.d[GZ]=0.0;

  bnods=(struct onode *)malloc((org->nnode)*sizeof(struct onode));
  if(bnods==NULL)
  {
    MessageBox(NULL,"Buffer Null.","CreateOrganNode",MB_OK);
    return NULL;
  }

  codemax=0;
  for(i=0;i<(org->nnode);i++)
  {
    if(!nodeontoscreen(*(org->nodes+i),&Nx,&Ny,vp)) /*PROJECTION*/
    {
      MessageBox(NULL,"Projection Failed.","CreateOrganNode",MB_OK);
      return NULL;
    }

    /*(bnods+i)->d[GX]=(double)( Nx-vp.Xo);
    (bnods+i)->d[GY]=(double)(-Ny+vp.Yo);*/
    (bnods+i)->d[GX]=(double)Nx;
    (bnods+i)->d[GY]=(double)Ny;
    (bnods+i)->d[GZ]=0.0;

    if(codemax<((org->nodes+i)->code))
    {
      codemax=(org->nodes+i)->code;
    }
  }

  found=0;
  for(i=0;i<(org->nelem);i++)
  {
    for(j=0;j<((org->elems+i)->nban);j++)
    {
      cban=(org->elems+i)->bans+j;

      ban.nods=(struct onode **)malloc((cban->nnod)
                                       *sizeof(struct onode *));
      if(ban.nods==NULL)
      {
        MessageBox(NULL,"Buffer Null.","CreateOrganNode",MB_OK);
        return NULL;
      }

      ban.nnod=cban->nnod;

      for(k=0;k<(cban->nnod);k++)
      {
        loff=(*(cban->nods+k))->loff;
        (*(ban.nods+k))=(bnods+loff);
      }

      found=insideban(ban,nmouse);

      free(ban.nods);
      if(found==1) break;
    }
    if(found==1) break;
  }
  free(bnods);

  if(!found) /*CREATE ON GROUND IF NO INTERSECTION.*/
  {
    pl.nods[0].d[GX]=0.0;
    pl.nods[0].d[GY]=0.0;
    pl.nods[0].d[GZ]=0.0;

    pl.nvec.dc[GX]=0.0;
    pl.nvec.dc[GY]=0.0;
    pl.nvec.dc[GZ]=1.0;

    pl.a=0.0; pl.b=0.0; pl.c=1.0; pl.d=0.0; /*PLANE:Z=0*/
  }
  else
  {
    bantoplane(*cban,&pl);
  }

  dmx=(double)( mx-vp.Xo);
  dmy=(double)(-my+vp.Yo);
  createnodeonplane(vp,dmx,dmy,pl,&nadd);

  if(nadd.code==0)
  {
    MessageBox(NULL,"Intersection Failed.","CreateOrganNode",MB_OK);
    return NULL; /*FAILED BY INTERSECTION.*/
  }
  loff=0;
  while(code>((org->nodes+loff)->code) &&
        loff<(org->nnode)) loff++;

  if(code==0 || code>=((org->nodes+loff)->code)) /*ZERO OR SAME.*/
  {
    code=codemax+1;
    loff=org->nnode;
  }

  nadd.code=code;
  nadd.loff=loff;

/*sprintf(non,"ADD NODE:CODE=%ld,LOFF=%ld",
        nadd.code,nadd.loff);
MessageBox(NULL,non,"Create",MB_OK);*/

  return addnode(nadd,org);
}/*createorgannode*/

struct onode findlastnode(long int code,
                          struct viewparam vp,
                          long int mx,long int my,
                          struct plane *pe, /*ELEMENT PLANE.*/
                          struct organ *org)
/*FIND LAST NODE OF ELEMENT WHILE CREATING.*/
{
  /*char non[256];*/
  int found;
  int Nx,Ny;
  long int i,j,k,loff;
  long int codemax;
  double dmx,dmy;
  double eps=1.0E-4;
  struct onode nmouse; /*MOUSE POSITION AS NODE.*/
  struct onode *bnods;
  struct onode nfind={0,0,{0.0,0.0,0.0}};
  struct obans ban,*cban;
  struct plane pb;
  struct line lcross;

  double **drccos,**tdrccos;
  double nearest1,nearest2,distance;
  struct plane proj;
  struct line vline;
  struct onode nview,norigin;

  nmouse.d[GX]=(double)mx;
  nmouse.d[GY]=(double)my;
  nmouse.d[GZ]=0.0;

  bnods=(struct onode *)malloc((org->nnode)*sizeof(struct onode));
  if(bnods==NULL)
  {
    MessageBox(NULL,"Buffer Null.","FindLastNode",MB_OK);
    return nfind;
  }

  codemax=0;
  for(i=0;i<(org->nnode);i++)
  {
    if(!nodeontoscreen(*(org->nodes+i),&Nx,&Ny,vp)) /*PROJECTION*/
    {
      MessageBox(NULL,"Projection Failed.","FindLastNode",MB_OK);
      return nfind;
    }

    /*(bnods+i)->d[GX]=(double)( Nx-vp.Xo);
    (bnods+i)->d[GY]=(double)(-Ny+vp.Yo);*/
    (bnods+i)->d[GX]=(double)Nx;
    (bnods+i)->d[GY]=(double)Ny;
    (bnods+i)->d[GZ]=0.0;

    if(codemax<((org->nodes+i)->code))
    {
      codemax=(org->nodes+i)->code;
    }
  }

  found=0;
  for(i=0;i<(org->nelem);i++)
  {
    for(j=0;j<((org->elems+i)->nban);j++)
    {
      cban=(org->elems+i)->bans+j;

      ban.nods=(struct onode **)malloc((cban->nnod)
                                       *sizeof(struct onode *));
      if(ban.nods==NULL)
      {
        MessageBox(NULL,"Buffer Null.","FindLastNode",MB_OK);
        return nfind;
      }

      ban.nnod=cban->nnod;

      for(k=0;k<(cban->nnod);k++)
      {
        loff=(*(cban->nods+k))->loff;
        (*(ban.nods+k))=(bnods+loff);
      }

      found=insideban(ban,nmouse);

      free(ban.nods);
      if(found==1) break;
    }
    if(found==1) break;
  }
  free(bnods);

  if(!found) /*CREATE ON GROUND IF NO INTERSECTION.*/
  {
    pb.nods[0].d[GX]=0.0;
    pb.nods[0].d[GY]=0.0;
    pb.nods[0].d[GZ]=0.0;

    pb.nvec.dc[GX]=0.0;
    pb.nvec.dc[GY]=0.0;
    pb.nvec.dc[GZ]=1.0;

    pb.a=0.0; pb.b=0.0; pb.c=1.0; pb.d=0.0; /*PLANE:Z=0*/
  }
  else
  {
    bantoplane(*cban,&pb);
  }

  if(pe==NULL) /*WIRE ELEMENT*/
  {
    dmx=(double)( mx-vp.Xo);
    dmy=(double)(-my+vp.Yo);
    createnodeonplane(vp,dmx,dmy,pb,&nfind);
  }
  else
  {
    distance=intersectplaneplane((*pe),pb,&lcross);
    if((-eps)<distance && distance<eps && lcross.code==0)
    {
      dmx=(double)( mx-vp.Xo);
      dmy=(double)(-my+vp.Yo);
      createnodeonplane(vp,dmx,dmy,(*pe),&nfind);
    }
    else if(lcross.code==0)
    {
      MessageBox(NULL,"Intersection Failed.","FindLastNode",MB_OK);
      return nfind; /*FAILED.*/
    }
    else
    {
      /*SCREEN ORIGIN TO GLOBAL.*/
      norigin.d[GX]=(vp.r-vp.odv)/vp.r*vp.ov[GX]+vp.focus.d[GX];
      norigin.d[GY]=(vp.r-vp.odv)/vp.r*vp.ov[GY]+vp.focus.d[GY];
      norigin.d[GZ]=(vp.r-vp.odv)/vp.r*vp.ov[GZ]+vp.focus.d[GZ];

      /*DEFINITION PROJECTION PLANE.*/
      proj.nvec.dc[GX]=vp.ov[GX];
      proj.nvec.dc[GY]=vp.ov[GY];
      proj.nvec.dc[GZ]=vp.ov[GZ];

      proj.nods[0].d[GX]=norigin.d[GX];
      proj.nods[0].d[GY]=norigin.d[GY];
      proj.nods[0].d[GZ]=norigin.d[GZ];
      proj.nods[1].d[GX]=norigin.d[GX]+vp.axono[1][0]; /*AXIS y.*/
      proj.nods[1].d[GY]=norigin.d[GY]+vp.axono[1][1];
      proj.nods[1].d[GZ]=norigin.d[GZ]+vp.axono[1][2];
      proj.nods[2].d[GX]=norigin.d[GX]+vp.axono[2][0]; /*AXIS z.*/
      proj.nods[2].d[GY]=norigin.d[GY]+vp.axono[2][1];
      proj.nods[2].d[GZ]=norigin.d[GZ]+vp.axono[2][2];

      /*MOUSE POINT TO GLOBAL.*/
      nmouse.d[EX]=(double)( mx-vp.Xo);
      nmouse.d[EY]=(double)(-my+vp.Yo);
      nmouse.d[EZ]=0.0;

      drccos=filmdrccos(proj.nods[0],
                        proj.nods[1],
                        proj.nods[2]);          /*DIRECTION COSINE.*/
      tdrccos=matrixtranspose(drccos,3);              /*MATRIX [T].*/
      nmouse=transltog(tdrccos,proj.nods[0],nmouse);

      for(i=0;i<3;i++) free(*(drccos+i));
      free(drccos);
      for(i=0;i<3;i++) free(*(tdrccos+i));
      free(tdrccos);

      /*VIEWPOINT.*/
      if(vp.type==PERSPECTIVE)
      {
        nview.d[GX]=vp.ov[GX]+vp.focus.d[GX];
        nview.d[GY]=vp.ov[GY]+vp.focus.d[GY];
        nview.d[GZ]=vp.ov[GZ]+vp.focus.d[GZ];
      }
      if(vp.type==AXONOMETRIC)
      {
        nview.d[GX]=nmouse.d[GX]+proj.nvec.dc[GX];
        nview.d[GY]=nmouse.d[GY]+proj.nvec.dc[GY];
        nview.d[GZ]=nmouse.d[GZ]+proj.nvec.dc[GZ];
      }

      /*INTERSECTION.*/
      vline.ends[0]=nview;
      vline.ends[1]=nmouse;

      if(!distancelineline(lcross,vline,
                           &nearest1,&nearest2,
                           &distance))
      {
        MessageBox(NULL,"Measurement Failed.","FindLastNode",MB_OK);
        return nfind;
      }

      nfind.code=1;
      nfind.d[GX]=lcross.ends[0].d[GX]
                 -nearest1*(lcross.ends[1].d[GX]
                           -lcross.ends[0].d[GX]);
      nfind.d[GY]=lcross.ends[0].d[GY]
                 -nearest1*(lcross.ends[1].d[GY]
                           -lcross.ends[0].d[GY]);
      nfind.d[GZ]=lcross.ends[0].d[GZ]
                 -nearest1*(lcross.ends[1].d[GZ]
                           -lcross.ends[0].d[GZ]);
    }
  }

  if(nfind.code!=0)
  {
    loff=0;
    while(code>((org->nodes+loff)->code) &&
          loff<(org->nnode)) loff++;

    if(code==0 || code>=((org->nodes+loff)->code)) /*ZERO OR SAME.*/
    {
      code=codemax+1;
      loff=org->nnode;
    }
    nfind.code=code;
    nfind.loff=loff;
  }

  return nfind;
}/*findlastnode*/

void chaseelement(struct oelem *elem,
                  struct viewparam vp,
                  struct windowparams wp)
{
  HWND hoya;
  HDC hdc;
  POINT pc,pp;
  long int x,y,maxX,maxY;

  if(wp.hdcB != NULL && wp.hdcC != NULL)
  {
    pc.x=0;
    pc.y=0;
    ClientToScreen(wp.hwnd,&pc);

    hoya=GetParent(wp.hwnd);
    pp.x=0;
    pp.y=0;
    ClientToScreen(hoya,&pp);

    x=pp.x-pc.x;
    y=pp.y-pc.y;

    getwindowsize(hoya,&maxX,&maxY);

    /*PatBlt(wp.hdcC,x,y,maxX,maxY,PATCOPY);*/
    drawelement(wp.hdcC,vp,*(elem),ONSCREEN);

    /*BitBlt(wp.hdcB,x,y,maxX,maxY,wp.hdcC,x,y,
           dwrop);*/ /*GROUND BLACK:SRCPAINT WHITE:SRCAND.*/

    hdc = GetDC(wp.hwnd);
    BitBlt(hdc,x,y,maxX,maxY,wp.hdcB,x,y,SRCCOPY);
    BitBlt(hdc,x,y,maxX,maxY,wp.hdcC,x,y,SRCPAINT);
    ReleaseDC(wp.hwnd,hdc);

    PatBlt(wp.hdcC,x,y,maxX,maxY,PATCOPY);
  }

  return;
}/*chaseelement*/


