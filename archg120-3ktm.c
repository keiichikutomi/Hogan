
/*ARCHG120-2（白島応答用）との同期未完了, directioncosine と directioncosineII のみコピー*/

/*ARCLM101.C FOR WIN32 SINCE 1995.11.24.JUNSATO.*/
/*LAST CHANGE:22-Jan-2013.*/

/*ENABLE CREATION OF FRAME "CHANGE","REFER".....NOT YET.*/
/*TRANSLATION OF INPUTFILE ORGAN INTO ARCLM UNAVAILABLE.*/
/*SLAB DIVISION FOR CMQ.....NOT YET.*/
/*SORT COMPONENTS.....NOT YET.*/
/*IMPROVE "CROUTL]UDECOMPOSITION".....NOT YET.*/

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

/*#include <windows.h> */
/*#include <stdio.h>*/
#include <stdlib.h>
/*#include <string.h> */
#include <malloc.h>
#include <math.h>
#include <limits.h>
#include <time.h>

/*#include "canhead.h"*/                /*DEFINITION OF COMMAND ID.*/

/*#ifndef PI*/
#define PI 3.1415926535897932384
/*#endif*/

#define SIUNIT 9.80665                             /*SI UNIT FACTOR*/

/*#define DSAFETY 0.0500*/             /*INCREMENT OF SAFETY FACTOR*/
/*TEST FOR 1000 NODES:DSAFETY=0.05*/
/*TEST FOR 8 NODES INCREMENTAL:DSAFETY=0.0005 LAPS=40*/
#define LAPS 3                                           /*ALL LAPS*/
#define WAIT 0                                     /*WAIT TIME[sec]*/
#define EXPONENT 1.500                  /*EXPONENT OF YIELD SURFACE*/
#define RADIUS 0.950                      /*RADIUS OF YIELD SURFACE*/
/*#define RADIUS 0.600*/                      /*RADIUS OF YIELD SURFACE*/  /*for toshima HachnKanYLv2*/
/*#define RADIUS 1.1875*/ /*=0.950x1.250*/        /*RADIUS OF YIELD SURFACE*/
/*#define RADIUS 1.0925*/ /*=1.150x0.75*/     /*RADIUS OF YIELD SURFACE*/

/*AXIS OF YIELD SURFACE.0:Nz 1:Qx 2:Qy 3:Mz 4:Mx 5:My*/
#define SURFACEX 2 // default: 2
#define SURFACEY 4 // default: 4
#define SURFACEZ 0

//#define DIRECTORY "..\\Data\\" /*DATA DIRECTORY*/
#define DIRECTORY "C:\\Users\\keiic\\Desktop\\Hogan\\Data\\" /*DATA DIRECTORY*/
//#define DIRECTORY "\0"               /*CURRENT*/

#define BOXMOUSE 5 /*BOX MOUSE SIZE.*/

#define GX 0
#define GY 1
#define GZ 2
#define EX 1 /*ANOTHER CASE 0*/
#define EY 2 /*             1*/
#define EZ 0 /*             2*/

#define RED   0
#define GREEN 1
#define BLUE  2

#define ROLENULL   0
#define ROLEWEIGHT 1
#define ROLERIGID  2
#define ROLESTRESS 3
#define ROLEHOJO   4

#define MAXTYPE  8 /*ELEMENT TYPES*/
#define TYPENULL 0 /*ELEMENT TYPE*/
#define COLUMN   1
#define GIRDER   2
#define BEAM     3
#define BRACE    4
#define WALL     5
#define SLAB     6
/*#define CURTAIN  7*/
#define NONLINEARBEAM 7 /*UJIOKA FOR PEZETTINO*/

#define CTYPE_DOT     1 /*CURVE TYPE*/
#define CTYPE_LINE    2
#define CTYPE_CIRCLE  3
#define CTYPE_ELLIPSE 4
#define CTYPE_SPLINE  5

#define CTYPE_FB      6
#define CTYPE_L       7
#define CTYPE_H       8
#define CTYPE_CT      9
#define CTYPE_C      10
#define CTYPE_O      11
#define CTYPE_BOX    12

#define WEIGHTSLAB  0 /*WEIGHT TYPE*/
#define WEIGHTFRAME 1
#define WEIGHTEQ    2
#define WEIGHTEQX   3
#define WEIGHTEQY   4


#define LOADZ 0
#define LOADX 1
#define LOADY 2

#define AXONOMETRIC 0
#define PERSPECTIVE 1

#define BOLDRED 0  /*UJIOKA FOR COLOR*/

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

struct onode{long int code,loff;     /*CODE NUMBER,OFFSET POSITION.*/   /*code=節点番号, loff=通し番号*/
			 double d[3];
			 double r[3];};          /*ROTATIONAL DOF BY KUTOMI*/
struct oconf{signed char iconf;
			 double value;}; /*ONLY FOR ANALYSIS.*/
struct oprop{long int code,loff;
			 char *name;
			 double hiju,E,poi;
			 double rfle[3],rfra[3];
			 int r,g,b; /*R,G,B.*/
			 double F,fut,fuc;};                  /*ULTIMATE STRESS*/
struct ofigs{long int code,loff;
			 double area,Ixx,Iyy,Jzz; /*COEFFICIENTS.*/
			 double thick;
			 struct oprop *prop;}; /*SECTION FIGURES.*/

struct curve{long int loff;
			 int type,hugo;
			 double radius[2],cangle;       /*RADIUSES,COORD ANGLE.*/
			 double angle[2];                      /*ANGLE OF ENDS.*/
			 struct onode *center;
			 struct onode *(dots[3]);
			 struct line *(tan[3]);
			};
struct polycurve{long int loff;
				 int ncurve,type;
				 struct curve *curves;
				 /*int r,g,b;*/
				 struct oprop prop;};
struct polypolycurve{int npcurve;
					 struct polycurve *pcurves;
					 char name[256];
					 ICONINFO ici;
					 HICON hico;
					 HCURSOR hcur;};
struct features{double Ax,Ay,Sx,Sy,Ixx,Iyy,Igx,Igy,Zx[2],Zy[2],
					   Hx,Hy,hx[2],hy[2],ix,iy;
				double Jzz,hiju;
				double E,poi;
				struct onode Gx,Gy,Gg;};

struct osect{long int code,loff;
			 long int ocode;
			 char *name;
			 char dflag; /*DRAW FLAG.*/
			 int nfig;
			 double E,poi,area,Ixx,Iyy,Jzz; /*TOTAL COEFFICIENTS.*/
			 double exp,fmax[6],fmin[6]; /*YIELD SURFACE.*/
			 double yieldinit,yieldcoefficient,yieldpower;/*ISOTROPIC MATERIAL HARDENING*/
			 double hiju[3];
			 struct ofigs *figs;
			 struct rgbcolor dcolor;
			 int role,type;
			 struct polypolycurve ppc;
			 double lload[3];
			 double perpl[3];    /*PERPENDICULAR LOAD*/
			 long int smax;};
struct obans{long int code,loff;
             int nnod;
             struct onode **nods;}; /*PLANE OF ELEMENT.*/
struct surface{double exp,fmax[6],fmin[6];};

struct gausspoint{double estrain[7];
				  double pstrain[7];
				  double stress[7];
				  double backstress[7];
				  double qn,qm,qnm;
				  double yinit,y,alpha;
				  double f[2];
				  double lambda[2];
				  double Ee,Ep;
				  };

struct owire{long int code,loff;
			 int nnod;
			 double cangle;
			 signed char iconf[2][6];
			 double stress[2][6];
			 struct onode *(node[2]);
			 struct osect *sect;

			 double srate[4]; /*SAFETY RATE OF SRCAN RESULT.*/

			 struct surface allow,yield; /*YIELD SURFACE.*/

			 double Ee[2],Ep[2]; /*STRAIN ENERGY.*/
			}; /*WIRE ELEM FOR ARCLM001,101.*/

struct memoryelem{
				   long int code;
				   signed char bond[2][6];
				   double stress[2][6];
				 };
							  /*ELEMENT MEMORY FOR ARCLM.*/
struct oshell{long int code,loff;
			  int nnod,ngp,nstress;
			  struct onode *(node[4]);
			  struct osect *sect;
			  double stress[4][6];
			  double area,w[9];
			  struct gausspoint gp[9];/*INTEGRATION POINTS PARAMS*/

			  double prate;/*IN-PLANE STIFFNESS RATIO*/
			  double brate;/*BENDING STIFFNESS RATIO*/
			 }; /*SHELL ELEM*/

struct memoryshell{
					long int code;
					double stress[4][6];
					struct gausspoint gp[9];/*INTEGRATION POINTS PARAMS*/
				  };                     /*SHELL ELEMENT MEMORY FOR ARCLM.*/

struct oconstraint{long int code,loff;
				   int nnod;
				   int type;
				   int neq,leq;
				   //double *lambda;
				   struct onode *(node[2]);
				   double axis[3][3];
				  };


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
             struct elemcolor color;
			 double lface[2],hface[2],wrect[2];};

struct oload{struct onode *nod;
//             double w[3];     /*WEIGHT FOR SLAB,FRAME,EARTHQUAKE.*/
             double w[5];       /*WEIGHT FOR SLAB,FRAME,EARTHQUAKE(XY).*/
//             double aifact;   /*FACTOR FROM Ai LOAD.*/     /*LOADS.*/
             double aifact[2];  /*FACTOR FROM Ai LOAD.*/     /*LOADS.*/
             double mass;};

struct aiparam{int nfloor;
               char **fnames; /*FLOOR NAMES.*/
               int *nnode;
//               double hmax,T1,Tc,Z,Rt,Co,Cf;
               double hmax,T1,Tc,Z,Rt,Cox,Coy,Cfx,Cfy;
               double *lbound,*levels;
//               double *wi,*Wi,*Ai,*Ci,*Qi,*Hi,*facts;
               double *wi,*Wi,*Ai,*Cix,*Ciy,*Qix,*Qiy,*Hix,*Hiy,*factsx,*factsy;    //ujioka for BASE
               double *Ds,*Fes;
              }; /*PARAMETERS FOR Ai LOADS.*/

struct organ{long int code,loff;
             long int nnode,nelem,nprop,nsect;
             struct onode *nodes;
             struct oconf *confs;
             struct oprop *props;
             struct osect *sects;
             struct oelem *elems;
             struct orgcolor color;
             double opaque;
             struct oload *loads;
             int ntext;
             struct snode *texts;                /*TEXTS ON WINDOW.*/
             struct aiparam ai;};

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
struct cmqelem{double Ci,Cj,Mc0,Qi0,Qj0;};

struct viewrange{struct onode max,min;};              /*VIEW RANGE.*/
struct selectrange{int il,ir,it,ib;};               /*SELECT RANGE.*/
struct menuvisible{char draw,
                        error,
                        surface,
                        sectlist,
                        weight,
                        cmq,
                        horizon,
                        ftype,    /*FILE TYPE.*/
                        savetype, /*SAVE AS TYPE.*/
                        inputfile,outputfile,
                        title,pagetitle,pagenum,
                        view;
                  };                /*VISIBLE FLAG FOR OPTION MENU.*/

///////////////////////////////////////////by MIHARA for Print Range Line
struct printvisible{char printrange,
                         a4tate,
                         a4yoko,
                         a3tate,
                         a3yoko;
                  };                /*VISIBLE FLAG FOR PRINT RANGE.*/
/////////////////////////////////////////////////////////////////////////

struct nodevisible{char code,
                        d[3],    //mihara for zahyo 20190329
						loads[6],
                        confs[6],
                        disps[3],             /*DISPLACEMENT X,Y,Z.*/
                        react[6],
                        mcircle,mvalue,                      /*MASS*/
                        conffig;
                  };                       /*VISIBLE FLAG FOR NODE.*/
struct elemvisible{char code,
                        axis,
                        hinge,
						sectioncode,
						sectionshape, //honda
                        deformation,
                        stress[8][6],          /*Nz,Qx,Qy,Mz,Mx,My.*/
						cmqline,               /*CMQ DIVISION LINE.*/
                        cmqcheck,              /*CMQ CHECK*/ //kaza & uji
						etype[8],                   /*ELEMENT TYPE.*/
						srcanrate,srcancolor,srcanmax,
						ecircle,plasticenergy,evalue;                    /*ENERGY*/
				  };                    /*VISIBLE FLAG FOR ELEMENT.*/
struct globvisible{char axis;
                   struct menuvisible mv;
                   struct printvisible pv;      //by MIHARA for Print Range Line
                   struct nodevisible nv;
                   struct elemvisible ev;
                  };                     /*VISIBLE FLAG FOR GLOBAL.*/
struct drawparam{double gaxis,
						eaxis,
                        dfact,
                        qfact,
                        mfact,
                        pitch,
						hsize,
						csize;
                }; /*PARAMETERS OF DRAWING FRAME.*/

struct viewparam{int type;           /*0:AXONOMETRICS 1:PERSPECTIVE*/
                 int Xo,Yo;                /*ORIGIN OF 2D DRAWINGS.*/
                 int chapter,section,subsection;      /*PAGE NUMBER*/
                 double theta,phi,r;
                 double gfactor;                  /*FACTOR OF SIZE.*/
                 struct onode focus;                 /*FOCUS POINT.*/
                 struct viewrange range;           /*RANGE VISIBLE.*/
				 struct drawparam dparam;        /*DRAW PARAMETERS.*/
                 double odv,vr,ov[3],axono[3][3];     /*PROJECTION.*/
                 struct globvisible vflag;         /*VISIBLE FLAGS.*/
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
                    int gstatus;                  /*GENERAL STATUS.*/
                    int lstatus,rstatus;            /*MOUSE STATUS.*/
                    struct windowparams *childs;   /*CHILD WINDOWS.*/

                    int tx,ty;             /*MESSAGE TEXT POSITION.*/
                    struct snode *strset;      /*STRINGS ON WINDOW.*/

                    /*RECT grect,lrect;*/              /*RECT SIZE.*/
                    RECT vbar,hbar;                   /*SCROLL BAR.*/

                    char inpfile[80],otpfile[80];      /*FILE NAME.*/
                    struct viewparam vparam;  /*3D VIEW PARAMETERS.*/
                    struct organ org;
                    char inpfilex[80],inpfiley[80],inpfilez[80];
					char otpfilex[80],otpfiley[80],otpfilez[80];
					char sctfile[80];                  /*FILE NAME.*/
					char title[256],pagetitle[256];
				   };
struct winparamsreg{
                    int code;
                    int nwin;
                    struct windowparams **wp;
				   };                   /*REGISTRY OF WINDOWPARAMS.*/

struct arclmframe{long int code,loff;
				  char *appelation;
				  int nnode,nelem,nshell,nsect,nreact,nlaps,neig,nconstraint;
				  double loadfactor;
				  double *eigenval,**eigenvec;
				  double *iform, *ddisp, *lambda;
				  struct memoryelem *melem;
				  struct memoryshell *mshell;
				  double *dreact;
				  FILE *fsurface;
				  struct onode *nodes,*ninit;
				  struct osect *sects;
				  struct owire *elems;
				  struct oshell *shells;
				  struct oconf *confs;
				  struct oconstraint *constraints;
				  long int *constraintmain;
				  //double *constraintval,**constraintvec;

				  double *nmass;/*NODE MASS FOR GNSHN.*/

				  int nosect; 		/****SRCANMAX****/
				  long int *iosect; /****SRCANMAX****/
				  int *eosect;		/****SRCANMAX****/
				  double *dosect; 	/****SRCANMAX****/
				 };/*ARCLM FRAME.*/



struct print{
			 char pflag; /*0:NOT PRINTING 1:UNDER PRINTING*/
			 char pageflag; /*0:PAGE BEGIN 1:PAGE WAITING*/
			 int margin[4];
             int jiheight,jiwidth,jipitch;
			 int gyopitch;
             int dans,dangap;
             int cWidthPels,cHeightPels,caps;
			 int cWidthDpi,cHeightDpi;
             PRINTDLG pd;
             DOCINFO di;
            };                                     /*PRINT OPTIONS.*/

/*EXTERNAL PARAMETERS*/
extern char prj[20];
extern int globalstatus, globalflag;
extern int globalmessageflag, globaldrawflag;     /*Conjugate*/
extern int globalcondensationflag;                /*CONDENSATION*/
extern struct windowparams wdraw,wmenu,wmesg,wsurf;
extern struct arclmframe arc; /*GLOBAL ARCLM FRAME.*/
long int comps; /*COMPONENTS IN GLOBAL MATRIX.*/
extern struct oelem gelem,*pelem;
extern int nmultinode,nmultielem,nmultiwire;
extern struct oelem **multielem; /*POINTING CURRENT ELEMENT.*/
extern struct onode **multinode; /*POINTING CURRENT ELEMENT.*/
extern struct owire **multiwire;
extern double globalunit; /*UNIT FACTOR*/
extern struct print gprn; /*GLOBAL PRINT PARAMETERS.*/
extern FILE *globalfile; /*GLOBAL FILE.*/
double totalweight =0.0;           //honda for sydney
double totallength =0.0;           //honda for sydney
double Gcx =0.0;           //honda for sydney
double Gcy =0.0;           //honda for sydney
double Gcz =0.0;           //honda for sydney

extern void modifycmq(struct memoryelem *melem,struct owire *elem);
extern void assemcmq(struct owire elem,double **tmatrix,
        	         struct oconf *confs,double *gvct);
extern void assemcmq101(struct owire elem,double **tmatrix,
                     struct oconf *confs,double *gvct,double dsafety);
extern void definencr(struct arclmframe *af,double *ncr);                   //ujioka
extern void drawcontrolepoint(HDC hdc,struct viewparam vp,struct onode gn);

/*=================================================================*/
/*REGULATED FOR ORGAN.*/
/*"INPUT":INPUT FROM TEXT FILE TO MEMORY.*/
/*"CURRENT":CURRENT DATA FROM BINARY FILE TO MEMORY.*/
/*"SELECT":SELECT FROM SCREEN BY MOUSE.*/
/*"GET":GET FROM OPTIONS DIALOG BOX.*/
/*"SET":SET WINDOW OBJECTS,SET TEXT INTO OPTIONS DIALOG BOX.*/

/*MATHEMATICS SUBROUTINES.*/
double *mallocdoublevector(int vsize);
double **mallocdoublematrix(int msize);
double **mallocdoublematrixxyz(int msize,int msize2);
void freematrix(double **mtx,int msize);

/*FOR VECTOR CALCURATION*//*KUTOMI*/
double dotproduct(double* vct1, double* vct2, int vsize);
double* crossproduct(double* vct1, double* vct2);
double vectorlength(double* vct, int vsize);
void vectornormalize(double* vct, int vsize);

double *matrixvector(double **mtx,double *vct,int msize);
void matrixvectorII(double *v,double **mtx,double *vct,int msize);
double *matrixvectorIII(double **mtx,double *vct,int m,int n);

//double *matrixvectorIII(double **mtx,double *vct,int msize); // 150220 fukushima for arclm201



double **matrixmatrix(double **mtxhead,double **mtxtail,int msize);
void matrixmatrixII(double **mtx,double **mtxhead,double **mtxtail,int msize);
double **matrixmatrixIII(double **mtxhead,double **mtxtail,int m,int n,int l);

double **matrixtranspose(double **mtx,int msize);
void matrixtransposeII(double **m,double **mtx,int msize);
double **matrixtransposeIII(double **mtx,int m,int n);

double **matrixinverse(double **mtx,int msize);
double **fullmatrixcroutlu(double **mtx,int msize);

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
					   struct onode *nearest,int *iside);
double distancedotplane(struct onode n,struct plane p,
                        struct onode *nearest);
int distancelineline(struct line l1,
					 struct line l2,
					 double *nearest1,
					 double *nearest2,
					 double *distance,int *iside);
double intersectlineplane(struct line l,struct plane p,
                          struct onode *ncross);
double intersectplaneplane(struct plane p1,struct plane p2,
                           struct line *lcross);
int divisionlineline(struct line l1,
                     struct line l2,
                     struct line *l0);

int intersectlineline(int x1,int y1,int x2,int y2,
                      int x3,int y3,int x4,int y4);
int intersectlinemouse(int ix1,int iy1,int ix2,int iy2,POINT point);
int insideban(struct obans b,struct onode n);
struct plane *bantoplane(struct obans b,struct plane *pl);

/*STRINGS SUBROUTINES.*/
char **fgetsbrk(FILE *fin,int *n);
char **fgetscut(FILE *fin,int *n);
void freestr(char **str,int nstr);
int fcancelstr(FILE *fin,const int nline);
struct snode *addtext(struct snode *strset,int *nstr,char *str,
                      double tx,double ty);

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
void setviewpoint(HWND hdwnd,
                  struct arclmframe af,struct viewparam *vp);
void getdoublefromdialog(HWND hdwnd,int iditem,double *d);
void setdoubleintodialog(HWND hdwnd,int iditem,double d);
void getmasshiju(HWND hdwnd,double *hiju);
void getkatakouparam(HWND hdwnd,
                     double *H,double *B,
                     double *tf1,double *tf2,double *tw,
                     double *ri1,double *ri2,
                     double *ri3,double *ri4,
                     double *ro1,double *ro2,
                     double *ro3,double *ro4);
void getproperty(HWND hdwnd,
                 long int *code,char *name,
                 double *E,double *poi,double *hiju);
void setlaps(HWND hdwnd,int ilap,int nlap);
void createviewdata(struct viewparam *vp);
int nodeprojection(struct onode ng,struct onode *np,
                   struct viewparam vp);
int nodeontoscreen(struct onode ng,int *ix,int *iy,
                   struct viewparam vp);
int nodeinsiderect(struct onode on,
                   int il,int ir,int it,int ib,
                   struct viewparam vp);
void drawglobaltext(HDC hdc,struct viewparam vp,
                    struct onode tn,char *str);
void drawglobaltextaligned(HDC hdc,struct viewparam vp,
                           struct onode tn,char *str,
                           int wtype,int htype);
void drawglobalnode(HDC hdc,struct viewparam vp,
					struct onode gn,struct oload *gl);
void drawgloballine(HDC hdc,struct viewparam vp,struct onode gn1,
                                                struct onode gn2);

//Modifyed by MIHARA 2007.07.11/////////////////////////////////////////////////
void drawprintrange(HDC hdc,struct viewparam vp);
void drawglobalconfig(HDC hdc,struct viewparam vp,struct organ go,
                                struct onode gn,int mode,long int loff);
////////////////////////////////////////////////////////////////////////////////
//081209 araki for conffig/////////////////////////////////////////////////
void drawglobalconfigofaf(HDC hdc,struct viewparam vp,struct arclmframe go,
                                struct onode gn,int mode,long int loff);
///////////////////////////////////////////////////////////////////////////////

double angleofgloballine(struct viewparam vp,struct onode gn1,
                                             struct onode gn2);
void drawlinestruct(HDC hdc,struct viewparam vp,struct line l);
void drawglobalarrow(HDC hdc,struct viewparam vp,
                     struct onode gn1,
                     struct onode gn2,
                     double asize);

/*******/
void sortlines(struct line *pl,int nline,double hugo);

int dividearc(struct curve *cv,
              double *dx,double *dy,double *angle,
              struct curve *dcv);
int dividearcbyline(struct curve *cv,struct line l,
                    struct curve *dcv);
int dividearcbylineext(struct curve *cv,
                       struct line l,double *fa,double *fb,
                       struct curve *dcv);

double definecircleangle(struct curve *cv,int iend,
                         double x,double y);
void chordrange(struct curve *cv,
                double *minx,double *miny,
                double *maxx,double *maxy);

int insidelinecurveint(int x,int y,
                       int x1,int y1,int x2,int y2,
                       int ox,int oy);
int insidelinecurveext(double x,double y,
                       struct curve *cv,
                       struct onode origin);
int dotinchord(struct curve *cv,double x,double y);
int dotinpolycurve(struct polycurve *pc,double x,double y);
struct polycurve *pickpolycurve(struct polypolycurve *ppc,
                                struct viewparam vp,
                                int mx,int my);
int arcinsidechord(struct curve *cvouter,struct curve *cvinner);

struct line tangentoncircle(struct curve *cv,
                            struct onode *on,double *angle);
void copycurve(struct curve *cvto,struct curve *cvfrom);
void copypolycurve(struct polycurve *pcto,
                   struct polycurve *pcfrom);
void copypolypolycurve(struct polypolycurve *ppcto,
                       struct polypolycurve *ppcfrom);

void drawlinecurveext(HDC hdc,
                      struct viewparam vp,
                      struct curve *cv,
                      HPEN hpwaku,HPEN hpfill,HBRUSH hbfill,
                      struct onode *origin);
void drawglobalchord(HDC hdc,
                     struct viewparam vp,
                     struct curve cv,
                     HPEN hpwaku,HPEN hpfill,HPEN hphollow,
                     HBRUSH hbfill,HBRUSH hbhollow);
void drawglobalchordext(HDC hdc,
                        struct viewparam vp,
                        struct curve cv,
                        HPEN hpwaku,
                        COLORREF crwaku,
                        COLORREF crfill,
                        COLORREF crhollow,
                        int fillorhollow,
                        struct onode origin);
void drawglobalpolycurve(HDC hdc,
                         struct viewparam vp,
                         struct polypolycurve *pc);
HICON createpropertyicon(HINSTANCE hinst,HDC hdc,
                         struct oprop *op,int iorc);
HICON createsectionicon(HDC hdc,int ix,int iy,
                        struct viewparam vp,
                        struct polypolycurve *pc);
void drawiconbmp(HDC hdc,int x,int y,
                 ICONINFO *ici,DWORD rop);
void rotatepolycurve(struct polypolycurve *pc,double theta);
/*******/

void drawcircleonglobaldot(HDC hdc,
                           struct viewparam vp,
                           struct onode center,double radius);
void fillglobalcircle(HDC hdc,struct viewparam vp,
                      struct onode center,double radius,
                      int r,int g,int b,double orate,
                      DWORD rop);
void fillglobalban(HDC hdc,struct viewparam vp,struct obans gb,
                   int r,int g,int b,
                   DWORD rop);
void drawglobalban(HDC hdc,struct viewparam vp,struct obans gb,
                   int r,int g,int b,
                   DWORD rop);
void drawelement(HDC hdc,struct viewparam vp,
                 struct oelem ge,
                 int mode,int idraw);
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
void drawtextalignedonlocaldot(HDC hdc,struct viewparam vp,
                               double **tdrccos,
                               struct onode on,
                               struct onode ln,
                               char *str,
                               int wtype,int htype);
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
                    struct arclmframe af,long int code,
                    int mode);
void drawglobalwire(HDC hdc,struct viewparam vp,
					struct arclmframe af,
					struct owire elem,
					int cred,int cgreen,int cblue,
					int ered,int egreen,int eblue,
					long int selectcode,
					int mode);
void drawglobalshell(HDC hdc,struct viewparam vp,
					struct arclmframe af,
					struct oshell shell,
					int cred,int cgreen,int cblue,
					int ered,int egreen,int eblue,
					long int selectcode,
					int mode);
void drawwirestress(HDC hdc,struct viewparam vp,
					struct arclmframe af,
					struct owire elem,int mode);
void drawwireaxis(HDC hdc,
                  struct viewparam vp,
				  struct onode n1,struct onode n2,double cangle);
void drawwireshape(HDC hdc,
				  struct viewparam vp,
				  struct onode n1,struct onode n2,double cangle,
				  int ered,int egreen,int eblue); //honda
void drawarclmframe(HDC hdc,struct viewparam vp,
                    struct arclmframe af,
					long int code,
                    int mode);
void savearclmasdxf(FILE *fout,struct viewparam vp,
                    struct arclmframe af);
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
struct onode *selectorgannode(struct viewparam vp,
                              struct organ *org,POINT point);
int selectelemnode(struct viewparam vp,struct oelem *oe,POINT point);
struct onode *selectpolycurvenode(struct viewparam vp,
                                  struct polypolycurve *pc,
                                  POINT point);
struct owire *getelement(HWND hdwnd,struct arclmframe *af,
                         long int code);
void setelement(HWND hdwnd,struct owire *elem);
void setelemconf(HWND hdwnd,signed char iconf,int idiconf);
struct owire *selectelement(struct viewparam vp,
                            struct arclmframe *af,POINT point);
struct oelem *selectorganelement(struct viewparam vp,
                                 struct organ *org,
                                 POINT point,
                                 struct line **selectedline,
                                 struct obans **selectedban);

/*ORGANIZATION SUBROUTINES.*/
void initializeorganization(struct organ *org);
int inputorganization(FILE *fin,
                      struct organ *org,
                      struct viewparam *vp);
int saveorganization(FILE *fout,
                     struct organ *org,
                     struct viewparam *vp);
int checkorganization(struct organ *org);
void freeorganization(struct organ *org);
int saveorganasdxf(FILE *fout,
                   struct organ go,
                   struct viewparam vp);
void outputtextasdxf(FILE *f,char str[],long int code,
                     double x,double y,double z,double size);
void outputlineasdxf(FILE *f,long int code,
                     char layer[],
                     double x1,double y1,double z1,
                     double x2,double y2,double z2);
void outputglobalarrowasdxf(FILE *f,
                            struct onode gn1,
                            struct onode gn2,
                            double asize,
                            long int *code);
void outputglobalaxisasdxf(FILE *f,long int *code,double size,
                           struct viewparam vp);
void outputglobalbanasdxf(FILE *f,struct viewparam vp,
                          struct obans gb,long int code);
int extractarclmfromorgan(struct organ *org,
                          struct arclmframe *az,
                          struct arclmframe *ax,
						  struct arclmframe *ay);
struct oelem *recttobrace(struct oelem *oe,
                          struct organ *org,
                          int loff,int nbrace,
                          double rfact);
int comparesection(struct osect *s1,struct osect *s2);
int comparefirstfig(struct osect *s1,struct osect *s2);
void slabdivision1(HDC hdc,struct viewparam vp,struct obans gban);
void slabdivision2(HDC hdc,struct viewparam *vp,
                   struct obans gban,struct obans *cmqbans);

struct rgbcolor setrgbcolor(int r,int g,int b);

/*ANALYSIS SUBROUTINES.*/
clock_t laptime(char *comment,clock_t t0);
void getincrement(HWND hdwnd,int *laps,double *dsafety);
void getlaps(HWND hdwnd,int *laps);
void setincrement(HWND hdwnd,
                  int laps,int lap,
                  double dsafes,double dsafe);
int decreaseband(long int *moff,long int *noff,
                 struct arclmframe *af);
int exchangelines(struct gcomponent *gmtx1,
                  double *gvct1,struct oconf *confs1,
                  struct gcomponent *gmtx2,
				  double *gvct2,struct oconf *confs2,
                  long int *moff,long int *noff,long int nnode,
                  char flag);
int exchangelinesII(double **kmtx1,double **gmtx1,struct oconf *confs1,
                  double **kmtx2,double **gmtx2,struct oconf *confs2,
                  long int *moff,long int *noff,long int nnode,
                  char flag);
int exchangelinesIIfloat(float **kmtx1,float **gmtx1,struct oconf *confs1,
                  float **kmtx2,float **gmtx2,struct oconf *confs2,
                  long int *moff,long int *noff,long int nnode,
                  char flag);
/*double croutludecomposition(struct gcomponent *gmtx,
                            double *gvct,struct oconf *confs,
							long int msize);*/

/* 150515 fukushima for all */
int croutlu(struct gcomponent *gmtx,
			struct oconf *confs,
			long int msize,
			double *det,double *sign,
			struct gcomponent *gcomp1);
int forwardbackward(struct gcomponent *gmtx,
					double *gvct, struct oconf *confs,
					long int msize,
					struct gcomponent *gcomp1);

int croutluII(struct gcomponent *gmtx,
			struct oconf *confs,
			long int msize, long int csize,
			double *det,double *sign,
			struct gcomponent *gcomp1);
int forwardbackwardII(struct gcomponent *gmtx,
					double *gvct, struct oconf *confs,
					long int msize, long int csize,
					struct gcomponent *gcomp1);


int croutludecomposition(struct gcomponent *gmtx,
						 double *gvct,struct oconf *confs,
						 long int msize,
						 double *det,double *sign);
int croutludecomposition_arclength(struct gcomponent *gmtx,
								   double *gvct,double *gvct2,struct oconf *confs,
								   long int msize,
								   double *det,double *sign,int iteration);
void currentpivot(long int i,long int msize);
int gread(struct gcomponent *gmtx,
		  long int i,long int j,double *data);
int gwrite(struct gcomponent *gmtx,
		   long int i,long int j,double data);
void gfree(struct gcomponent *gmtx,long int nnode);
int vread(FILE *fvct,long int i,double *data);
int vwrite(FILE *fvct,long int i,double *data);

int arclm101(struct arclmframe *af,int idinput);
int arclm101_bc(struct arclmframe *af,int idinput);  /*by Ujioka for Buckling Condensation Application*/
int arclm201(struct arclmframe *af,int idinput);  /*by Mihara for Non-Linear Shape Analysis 100510*/
DWORD availablephysicalmemory(char *comment);
DWORDLONG availablephysicalmemoryEx(char *comment);
FILE *fgetstofopen(const char *directory,const char *mode,
				   int dlgitem);
FILE *fgetstofopenII(const char *directory,const char *mode, const char *filename);        //


extern double shellarea(struct oshell shell);
void inputtexttomemory(FILE *ftext,struct arclmframe *af);
void inputloadtomemory(FILE *ftext,struct arclmframe *af);
void inputframetomemory(FILE *ftext,struct arclmframe *af);
int saveasarclm(char *fname,struct arclmframe *af);
int savebanddecreasedarclm(struct arclmframe *af);
int savebanddecreasedorgan(FILE *fout,
                           struct organ *org,
                           struct viewparam *vp,
                           struct arclmframe *af); /*Araki*/
void initialform(struct onode *nodes,double *ddisp,int nnode);

int initialnode(struct onode *nodes,int nnode,int code,
				struct onode *node);
void initialelem(struct owire *elems,
				 struct memoryelem *melem,int nelem);
void initialshell(struct oshell *shells,
				 struct memoryshell *mshell,int nshell);
void initialreact(FILE *fin,double *dreact,int nreact);
void inputinit(FILE *fin,int *nnode,int *nelem,int *nsect);
void inputinitII(FILE *fin,int *nnode,int *nelem,int *nshell,int *nsect,int *nconstraint);
void inputframeinit(FILE *fin,int *nnode,int *nelem,int *nsect);
void inputnode(double *ddisp,struct onode *node);
void inputelem(struct owire *elems,
			   struct memoryelem *melem,int offset,
			   struct owire *elem);
void readsect(FILE *fin,struct osect *sect);
void readsrcanrate(FILE *fin,struct arclmframe *af);
double **directioncosine(double x1,double y1,double z1,
						 double x2,double y2,double z2,
                         double cangle);
void directioncosineII(double x1,double y1,double z1,
                       double x2,double y2,double z2,
					   double cangle,
                       double **drccos);
double **filmdrccos(struct onode n1,struct onode n2,struct onode n3);
double **transmatrix(/*struct owire elem,*/double **drccos);
void transmatrixII(double **drccos,double **t);
double **transmatrixIII(double **drccos,int nnod);
double **assememtx(struct owire elem);
void assememtxII(struct owire elem,double **e);
double **assemgmtx(struct owire elem,double *estress);
double **assempmtx(struct owire elem,double **estiff);
double **assempmtxbc(struct owire elem,double **estiff,double ncr);



void coefficients(struct owire elem,double **estiff,
				  double f[],double dfdp[][6],
				  double q[][2],double a[][2]);
void coefficientsbc(struct owire elem,double **estiff,
				  double f[],double dfdp[][6],
				  double q[][2],double a[][2],double ncr);
double **modifyhinge(struct owire elem,double **estiff);

double **transformation(double **estiff,double **tmatrix);
void transformationII(double **estiff,double **tmatrix,
					  double **e,double **t);
double **transformationIII(double **estiff,double **tmatrix,int msize);
double **transformationEx(double **A,double **T,int row,int col);

void assemgstiffness(struct gcomponent *gmtx,
                     double **estiff,
					 struct owire *elem);
void assemgstiffnessII(struct gcomponent *gmtx,
					   double **estiff,
					   struct oshell *shell);
void assemgstiffnesswithDOFelimination(struct gcomponent *gmtx,
					   double **estiff,
					   struct owire *elem,
					   long int *constraintmain);
void assemgstiffnessIIwithDOFelimination(struct gcomponent *gmtx,
					   double **estiff,
					   struct oshell *shell,
					   long int *constraintmain);


void assemconf(struct oconf *confs,double *gvct,double dsafety,
			   int nnode);
void assemnodenorm(struct oshell shell,double *gvct);
void modifygivend(struct gcomponent *gmtx,double *gvct,
				  struct oconf *confs,int nnode);
void modifycmq201(struct memoryelem *melem,struct owire *elem);
void assemcmq201(struct owire elem,double **tmatrix,
			  struct oconf *confs,double *gvct);
void assemconf201(struct oconf *confs,double *gvct,
                  double dsafety,int nnode);
double *extractdisplacement(struct owire elem,double *gvct);
double *extractdisplacement2(struct owire elem,double *gvct,
							 long int *moff);
void extractdisplacementII(double *d,
						   struct owire elem,double *gvct);
double *extractshelldisplacement(struct oshell shell,double *ddisp);

double *elemstress(struct owire *elem,double *gvct,
                   struct memoryelem *melem,FILE *fout,
                   double func[]);
double *elemstressbc(struct owire *elem,double *gvct,
                   struct memoryelem *melem,FILE *fout,FILE *fsrf,
                   double func[],double ncr);
void elemstressII(double *estress,
				  struct owire *elem,
                  double *gvct,struct memoryelem *melem,FILE *fout,
                  double **drccos,
                  double **tmatrix,
                  double **estiff,
                  double *gdisp,double *edisp,
                  double func[],FILE *ftxt);
double *elemstressnl(struct owire *elem,double *gvct,
				   struct memoryelem *melem);
void updatestress(struct memoryelem *melem,FILE *fout,
                  double *edisp,double *dstress,double **estiff,
                  struct owire *elem,double func[],FILE *ftxt);
void updatestressbc(struct memoryelem *melem,FILE *fout,FILE *fsrf,
                  double *edisp,double *dstress,double **estiff,
                  struct owire *elem,double func[],FILE *ftxt,double ncr);
void outputdisp(double *gvct,FILE *fout,int nnode,
                struct onode *nodes);
void updatestressnl(struct memoryelem *melem,double *dstress,
					 struct owire *elem);
/*** 09.08.28 araki for Tsurukawa *********************************************/
void outputdisp02(double *gvct,double *displong,FILE *fout,int nnode,
                  struct onode *nodes);
/******************************************************************************/
void updateform(double *ddisp,double *gvct,int nnode);
void updateform2(double *ddisp,double *gvct,int nnode,long int *moff);

void copyform(struct arclmframe *af,double *gvct);
void outputstress(struct owire elem,double *estress,FILE *fout,
				  double func[]);
void outputshellstress(struct oshell shell,
				  double *estress,FILE *fout);
void outputstressnl(struct owire elem,
					 double *estress,FILE *fout);
void outputreaction(struct gcomponent *gmtx,
                    double *gvct,
                    struct onode *nodes,
					struct oconf *confs,
                    double *dreact,FILE *fout,int nnode);
void outputreactionnl(struct gcomponent *gmtx,
                       double *gvct,
                       struct onode *nodes,
                       struct oconf *confs,
                       double *dreact,FILE *fout,int nnode);
/*** 09.08.28 araki for Tsurukawa *********************************************/
void outputstress02(struct owire elem,double *estress,double *stresslong,int i,FILE *fout,
                    double func[]);
void outputreaction02(struct gcomponent *gmtx,
					  double *gvct,
                      struct onode *nodes,
                      struct oconf *confs,
                      double *dreact,double *reactlong,FILE *fout,int nnode);
/******************************************************************************/
void outputreactionwithoutupdate(struct gcomponent *gmtx,
                                 double *gvct,
                                 struct onode *nodes,
                                 struct oconf *confs,
                                 double *dreact,FILE *fout,FILE *ftxt,
                                 int nnode,int mhoko);
int  gcomponentadd(struct gcomponent *mtx1, // 150219 fukushima for gcomponentadd
                   struct gcomponent *mtx2,
				   int msize);

/*HOGANSHI SUBROUTINES.*/
void gethoganparam(HWND hdwnd,struct hoganparam *hp);
void drawhoganaxis(HDC hdc,
                   struct viewparam vp,
                   struct hoganparam hp,
                   int r,int g,int b);
int drawhoganlines(HDC hdc,HWND hdwnd,struct viewparam vp,
                   FILE *fin);

/*ORGAN CREATE SUBROUTINES.*/
struct oprop *addproperty(struct oprop *op,
                          struct organ *org);
struct osect *addsection(struct osect *os,
                         struct organ *org);
struct oconf *addconf(struct oconf conf[6],
                      long int nodeoffset,
                      struct organ *orgbefore);
struct oconf *changeconf(struct oconf conf[6],
                         long int nodeoffset,
                         struct oconf *cinit);
struct oconf *deleteconf(long int nodeoffset,
                         struct organ *orgbefore);
struct onode *addnode(struct onode node,
                      struct oconf *oc,
                      struct organ *orgbefore);
struct onode *deletenode(long int nodecode,
                         struct organ *orgbefore);
struct onode *deletenodewithelem(long int nodecode,
								 struct organ *orgbefore);
long int singulatenode(struct organ *orgbefore);
struct oelem *addelement(struct oelem *elem,
                         struct organ *orgbefore);
int addelementwithnode(struct oelem *pe,
                       struct organ *org,
                       double dx,double dy,double dz);
struct oelem *deleteelement(long int elemoffset,
                            struct organ *orgbefore,
                            long int pointednodecode);
void freeelement(struct oelem *elem,
                 char fnods,char fbonds,char fbans,char fbannods);
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
/*** araki for multinode*******************************************************/
struct onode **selectmultinode(struct organ *org,
                               struct viewparam vp,
                               int il,int ir,int it,int ib);
/******************************************************************************/
void chasenode(struct onode *node,
			   struct viewparam vp,
               struct windowparams wp);
void chaseelement(struct oelem *elem,
                  struct viewparam vp,
                  struct windowparams wp);
void chaserect(int il,int ir,int it,int ib,
               struct windowparams wp);
struct oelem **selectmultielem(struct organ *org,
                               struct viewparam vp,
                               int il,int ir,int it,int ib);
struct oelem **copymultielem(struct organ *org,
                             double dx,double dy,double dz);
struct oelem **movemultielem(struct organ *org,
                             double dx,double dy,double dz);
struct oelem **deletemultielem(struct organ *org);

/*** ujioka for multiwire*******************************************************/
struct owire **selectmultiwire(struct arclmframe *af,
                               struct viewparam vp,
                               int il,int ir,int it,int ib);
/******************************************************************************/

/*** araki for multinode*******************************************************/
struct onode **deletemultinode(struct organ *org);
struct onode **movemultinode(struct organ *org,
							 double dx,double dy,double dz);  //09.01.29 araki
/******************************************************************************/

/*POLYCURVE FEATURES SUBROUTINES.*/
void chasecurve(struct curve *cv,
                struct viewparam vp,
                struct windowparams wp);
void chasepolycurve(struct polycurve *pcv,
                    struct viewparam vp,
                    struct windowparams wp);
void chasepolypolycurve(struct polypolycurve *pcv,
                        struct viewparam vp,
                        struct windowparams wp);
int linefeatures(double x1,double y1,double x2,double y2,
                 struct features *f,
                 struct onode origin);
int circlefeatures(struct curve *cv,
                   struct features *f,
                   struct onode origin);
int curvefeatures(struct curve *cv,
                  struct features *f,
                  struct onode origin);
struct features polycurvefeatures(HWND hwnd,
                                  struct viewparam *vp,
                                  struct polypolycurve *pc);
double ultimatebendingofpolycurve(HWND hwnd,
                                  struct viewparam *vp,
                                  double Nu,
                                  struct polypolycurve *ppc);
int dividechord(struct polycurve *dpc,
                struct curve *ocv);
void boundaryvaluepolypolycurve(struct polypolycurve *ppc,
                                double Yn,double Yg,
                                double *Nb,double *Mb);
void boundaryvaluecurve(struct curve *cv,
                       double Yn,double Yg,
                       struct oprop op,
                       double *Nb,double *Mb);
struct features sectionfigfeatures(struct osect *os);
struct features sectionfeatures(struct osect *os);
void initializecurve(struct curve *cv);
struct curve *malloccurves(int ncurve);
void createcircledata(struct curve *cv);
void setcurveascircle(struct curve *cv,
                      int loff,int hugo,
                      double *radius,
                      double *angle1,double *angle2,
                      double cx,double cy,
                      double x1,double y1,
                      double x2,double y2);
void setcurveasline(struct curve *cv,
                    int loff,
                    double x1,double y1,
                    double x2,double y2);
void setpolycurveaskatakou(struct polycurve *pc,
                           int type,
                           double H,double B,
                           double tf1,double tf2,double tw,
                           double ri1,double ri2,
                           double ri3,double ri4,
                           double ro1,double ro2,
                           double ro3,double ro4);
void addcurvestopolycurve(struct polycurve *pc,
                          struct polycurve *add);
void deletecurveinpolycurve(struct polycurve *pc,int offset);
void freepolycurve(struct polycurve *pc);
void freepolypolycurve(struct polypolycurve *ppc);

/*ORGAN LOADS SUBROUTINES.*/
int cmqpolygon(struct obans cban,double wban,
               struct cmqelem *csum);
int cmquniform(struct oelem elem,double w0,
               struct cmqelem *csum);
int cmqconcentration(double W0,
                     double L0,
                     double L1,double L2,
                     struct cmqelem *cmq);
int cmqrightangledtriangle(double w0,
                           double L0,
                           double L1,double L2,double L3,
                           struct cmqelem *cmq);
void weightdistribution(HDC hdc,FILE *fout,struct viewparam *vp,
						struct organ *org);
/*SUBROUTINES IN SRCAL004.C*/
void sortdoublepointer(double *value,int nvalue);
void sortdouble(double value[],int nvalue);
int quadraticequation(double c2,double c1,double c0,
                      double answer[]);
int cubicequation(double c3,double c2,double c1,double c0,
                  double answer[]);

/*CADRE SUBROUTINES.*/
void inputcadretomemory(FILE *ftext,struct arclmframe *af);
void cadreoutputtomemory(FILE *ftext,FILE *fin,
                         struct arclmframe *af);

/*OBSERVATION POINTS*/    /*isebo and ujioka for dimple analysis*/

//void obspoints();

/*drawelemcmqcheck*/  /*FOR CMQ CHECK*/ /*kaza & uji for Lunar/MarsBase*/
void drawelemcmqcheck(HDC hdc,
				  struct viewparam vp,
				  struct onode n1,struct onode n2,double cangle,double pi,double pj);



double *mallocdoublevector(int vsize)
/*MALLOC DOUBLE VECTOR.*/
{
  double *v;

  v=(double *)malloc(vsize*sizeof(double));

  return v;
}/*mallocdoublevector*/

double **mallocdoublematrix(int msize)
/*MALLOC DOUBLE MATRIX.*/
{
  int i;
  double **mtx;

  mtx=(double **)malloc(msize*sizeof(double *));
  for(i=0;i<msize;i++)
  {
    *(mtx+i)=(double *)malloc(msize*sizeof(double));
  }

  return mtx;
}/*mallocdoublematrix*/

double **mallocdoublematrixxyz(int msize1,int msize2)/*Coded by Fujimoto*/
/*MALLOC DOUBLE MATRIX XYZ*/
{
	int i;
	double **mtx;

	mtx=(double **)malloc(msize1*sizeof(double *));
	for(i=0;i<msize1;i++)
	{
		*(mtx+i)=(double *)malloc(msize2*sizeof(double));
	}

	return mtx;
}/*mallocdoublematrixxyz*/

void freematrix(double **mtx,int msize)
/*FREE MATRIX.*/
{
  int i;

  for(i=0;i<msize;i++) free(*(mtx+i));
  free(mtx);

  return;
}/*freematrix*/


double dotproduct(double* vct1, double* vct2, int vsize)
{
	int i;
	double dot;
	dot = 0.0;
	for (i = 0; i < vsize; i++) dot += (*(vct1 + i)) * (*(vct2 + i));
	return dot;
}/*dotproduct*/

double* crossproduct(double* vct1, double* vct2)
{
	double* cross;
	cross = (double *)malloc(3*sizeof(double));

	*(cross + 0) =  *(vct1 + 1)**(vct2 + 2) - *(vct1 + 2)**(vct2 + 1);
	*(cross + 1) =  *(vct1 + 2)**(vct2 + 0) - *(vct1 + 0)**(vct2 + 2);
	*(cross + 2) =  *(vct1 + 0)**(vct2 + 1) - *(vct1 + 1)**(vct2 + 0);

	return cross;
}/*crossproduct*/



double vectorlength(double* vct, int vsize)
{
	int i;
	double len;
	len = 0.0;
	for (i = 0; i < vsize; i++) len += (*(vct + i)) * (*(vct + i));
	len = sqrt(len);
	return len;
}/*vectorlength*/

void vectornormalize(double* vct, int vsize)
{
	int i;
	double len;
	len = vectorlength(vct,vsize);
	if(len!=0.0)
	{
		for (i = 0; i < vsize; i++) *(vct + i)/=len;
	}
	return;
}/*vectorlength*/


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

void matrixvectorII(double *v,
                    double **mtx,double *vct,int msize)
/*MULTIPLY VECTOR BY MATRIX WITHOUT MALLOC.*/
{
  int i,j;

  for(i=0;i<msize;i++)
  {
    *(v+i)=0.0;
	for(j=0;j<msize;j++) *(v+i)+=(*(*(mtx+i)+j))*(*(vct+j));
  }

  return;
}/*matrixvectorII*/

#if 0
double *matrixvectorIII(double **mtx,double *vct,int msize) // 150220 fukushima for arclm201
/*MULTIPLY VECTOR BY MATRIX WITH free(vct).*/
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

  free(vct);

  return v;
}/*matrixvectorIII*/
#endif

double *matrixvectorIII(double **mtx,double *vct,int m,int n)
/*MULTIPLY MATRIX BY MATRIX.*/
{
  int i,j,jj;
  double* v;

  v=(double *)malloc(m*sizeof(double));
  if(v==NULL) return NULL;
  for(i=0;i<m;i++)
  {
	*(v+i)=0.0;
	for(j=0;j<n;j++)
	{
	  *(v+i)+=*(*(mtx+i)+j)**(vct+j);
	}
  }

  return v;
}/*matrixmatrixIII*/



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


void matrixmatrixII(double **mtx,
                    double **mtxhead,double **mtxtail,int msize)
/*MULTIPLY MATRIX BY MATRIX WITHOUT MALLOC.*/
{
  int i,j,jj;
  double *mline;

  for(i=0;i<msize;i++)
  {
    mline=*(mtx+i);
    for(j=0;j<msize;j++)
    {
      *(mline+j)=0.0;
	  for(jj=0;jj<msize;jj++)
      {*(mline+j)+=(*(*(mtxhead+i)+jj))*(*(*(mtxtail+jj)+j));}
    }
  }

  return;
}/*matrixmatrixII*/

double **matrixmatrixIII(double **mtxhead,double **mtxtail,int m,int n,int l)
/*MULTIPLY MATRIX BY MATRIX.*/
{
  int i,j,k;
  double **mtx,*mline;

  mtx=(double **)malloc(m*sizeof(double *));
  if(mtx==NULL) return NULL;
  for(i=0;i<m;i++)
  {
	mline=(double *)malloc(l*sizeof(double));
	if(mline==NULL) return NULL;
	for(j=0;j<l;j++)
	{
	  *(mline+j)=0.0;
	  for(k=0;k<n;k++)
	  {
		*(mline+j)+=*(*(mtxhead+i)+k)**(*(mtxtail+k)+j);
	  }
	}
	*(mtx+i)=mline;
  }

  return mtx;
}/*matrixmatrixIII*/

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

void matrixtransposeII(double **m,double **mtx,int msize)
/*MATRIX TRANSPOSITION.*/
{
  int i,j;

  for(i=0;i<msize;i++)
  {
	for(j=0;j<msize;j++) *(*(m+i)+j)=*(*(mtx+j)+i);
  }

  return;
}/*matrixtransposeII*/

double **matrixtransposeIII(double **mtx,int m,int n)
/*MATRIX TRANSPOSITION.*/
{
  double **mm;
  double *mline;
  int i,j;

  mm=(double **)malloc(n*sizeof(double *));
  if(mm==NULL) return NULL;
  for(i=0;i<n;i++)
  {
	mline=(double *)malloc(m*sizeof(double));
	if(mline==NULL) return NULL;
	for(j=0;j<m;j++)
	{
	  *(mline+j)=*(*(mtx+j)+i);
	}
	*(mm+i)=mline;
  }

  return mm;
}/*matrixtransposeIII*/

double **matrixinverse(double **mtx,int msize)
/*RETURN:INVERSE MATRIX BY GAUSSJORDAN SWEEP.*/
{
  char str[400],s[80];
  int i,j,k;
  double **m,det=1.0,pivot,data;

/*MessageBox(NULL,"Begin","Sweep",MB_OK);*/
  m=mallocdoublematrix(msize);
  for(i=0;i<msize;i++) /*ASSEMBLAGE UNIT MATRIX.*/
  {
	for(j=0;j<msize;j++)
	{
	  if(i==j) *(*(m+i)+j)=1.0;
	  else     *(*(m+i)+j)=0.0;
	}
  }

  for(i=0;i<msize;i++) /*PIVOT LINE.*/
  {
	pivot=*(*(mtx+i)+i);
	if(pivot==0.0)
	{
	  MessageBox(NULL,"INSTABLE.","Sweep",MB_OK);
	  return NULL; /*INSTABLE.*/
	}
	det*=pivot/fabs(pivot); /*DETERMINANT.*/
	for(j=0;j<msize;j++) /*PIVOT LINE DIVIDE BY PIVOT.*/
	{
	  *(*(mtx+i)+j)/=pivot;
	  *(*(m+i)+j)/=pivot;
	}
	for(j=0;j<msize;j++) /*FOR EACH LINES.*/
	{
	  if(j!=i)
	  {
        data=*(*(mtx+j)+i);
        for(k=0;k<msize;k++) /*FOR EACH ROWS.*/
		{
		  /*if(k!=i)*/
		  if(1)
		  {
			/**(*(mtx+j)+k)-=(*(*(mtx+j)+i))*(*(*(mtx+i)+k));*/
			/**(*(m+j)+k)-=(*(*(mtx+j)+i))*(*(*(m+i)+k));*/
			*(*(mtx+j)+k)-=data*(*(*(mtx+i)+k));
			*(*(m+j)+k)-=data*(*(*(m+i)+k));
		  }
		}

	  }
	}
  }

/*
sprintf(str,"\0");
for(i=0;i<msize;i++)
{
  for(j=0;j<msize;j++)
  {
	sprintf(s," %8.3f",*(*(m+i)+j));
	strcat(str,s);
  }
  strcat(str,"\n");
}
MessageBox(NULL,str,"Sweep",MB_OK);
*/

/*MessageBox(NULL,"End","Sweep",MB_OK);*/
  return m;
}/*matrixinverse*/

double **fullmatrixcroutlu(double **mtx,int msize)
/*RETURN:INVERSE MATRIX SYMMETRIC BY CROUT'S LU.*/
{
  int i,j,k;
  double **mo,**mi,det=1.0,pivot;

  mo=mallocdoublematrix(msize);
  if(mo==NULL) return NULL;
  for(i=0;i<msize;i++)
  {
	for(j=0;j<msize;j++) /*INITIALIZATION.*/
	{
	  *(*(mo+i)+j)=*(*(mtx+i)+j);
	}
  }
  mi=mallocdoublematrix(msize);
  if(mi==NULL) return NULL;
  for(i=0;i<msize;i++)
  {
    for(j=0;j<msize;j++) /*ASSEMBLAGE UNIT MATRIX.*/
    {
      if(i==j) *(*(mi+i)+j)=1.0;
      else     *(*(mi+i)+j)=0.0;
    }
  }

  for(j=1;j<=(msize-1);j++) /*DECOMPOSITION.*/
  {
    pivot=*(*(mo+j-1)+j-1);
    if(pivot==0.0)
    {
      MessageBox(NULL,"Terminated as Matrix Instable.",
                 "Crout",MB_OK);
      return NULL; /*INSTABLE.*/
    }

    /*det*=pivot;*/
    det*=pivot/fabs(pivot); /*SIGN OF DETERMINANT.*/

    for(i=(j+1);i<=msize;i++)
    {
      *(*(mo+i-1)+j-1)=*(*(mo+i-1)+j-1)/pivot;
    }
    for(i=(j+1);i<=msize;i++)
    {
      for(k=i;k<=msize;k++)
      {
        *(*(mo+k-1)+i-1)=*(*(mo+k-1)+i-1)
                         -pivot
                         *(*(*(mo+k-1)+j-1))
                         *(*(*(mo+i-1)+j-1));
      }
    }
  }
  for(i=2;i<=msize;i++) /*FORWARD ELIMINATION.*/
  {
    for(j=1;j<=(i-1);j++)
    {
      for(k=1;k<=msize;k++)
      {
        *(*(mi+i-1)+k-1)=*(*(mi+i-1)+k-1)-(*(*(mo+i-1)+j-1))
                                         *(*(*(mi+j-1)+k-1));
      }
    }
  }
  for(i=msize;i>=1;i--) /*BACKWARD SUBSTITUTION.*/
  {
    for(k=1;k<=msize;k++)
    {
      *(*(mi+i-1)+k-1)=(*(*(mi+i-1)+k-1))/(*(*(mo+i-1)+i-1));
      for(j=(i+1);j<=msize;j++)
      {
        *(*(mi+i-1)+k-1)=*(*(mi+i-1)+k-1)-(*(*(mo+i-1)+j-1))
                                         *(*(*(mi+j-1)+k-1));
      }
    }
  }

  freematrix(mo,msize);

  if(det<=0.0) errormessage("Completed but Matrix Instable.");

  return mi;
}/*fullmatrixcroutlu*/

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
  double dv,length;

  dv=vct.dc[GX]*vct.dc[GX]
    +vct.dc[GY]*vct.dc[GY]
    +vct.dc[GZ]*vct.dc[GZ];

if(dv<0.0) MessageBox(NULL,"Error.","Length",MB_OK);
  if(dv<0.0) return -1.0;

  length=sqrt(dv);

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
                       struct onode *nearest,int *iside)
/*RETURN:DISTANCE FROM DOT TO LINE,NEAREST DOT.*/
/*RETURN ISIDE:1=NEAREST DOT INSIDE LINE.*/
/*             2=NEAREST DOT OUTSIDE LINE.*/
{
  double inner,htlength,otlength,fact,distance;
  double eps=1.0E-04;
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

  fact=inner/htlength/htlength;
  nearest->d[GX]=l.ends[0].d[GX]+fact*vctht.dc[GX];
  nearest->d[GY]=l.ends[0].d[GY]+fact*vctht.dc[GY];
  nearest->d[GZ]=l.ends[0].d[GZ]+fact*vctht.dc[GZ];

  if(iside!=NULL)
  {
    if(eps<fact && fact<1.0-eps) *iside=1;
    else                         *iside=2;
  }

  return distance;
}/*distancedotline*/

double distancedotplane(struct onode n,struct plane p,
                        struct onode *nearest)
/*RETURN:DISTANCE FROM DOT TO PLANE,NEAREST DOT.*/
{
  double func,fact=0.0,nx,ny,nz,dn,distance=0.0;

  func=p.a*n.d[GX]+p.b*n.d[GY]+p.c*n.d[GZ]+p.d;   /*FUNCTION VALUE.*/

  nx=p.nvec.dc[GX];
  ny=p.nvec.dc[GY];
  nz=p.nvec.dc[GZ];

  dn=sqrt(nx*nx+ny*ny+nz*nz); /*LENGTH OF VECTOR.*/

  if(dn>0.0)
  {
    distance=func/dn; /*DISTANCE.*/
    fact=distance/dn; /*FACTOR OF VECTOR.*/
  }
  nearest->d[GX]=n.d[GX]-fact*nx; /*NEAREST POINT.*/
  nearest->d[GY]=n.d[GY]-fact*ny;
  nearest->d[GZ]=n.d[GZ]-fact*nz;

  return distance;
}/*distancedotplane*/

int distancelineline(struct line l1,
					 struct line l2,
					 double *nearest1,
					 double *nearest2,
					 double *distance,int *iside)
/*RETURN:FLAG OF SUCCESS,PARAMETER OF NEAREST DOTS.*/
/*RETURN:1=CROSSING INSIDE L1,ON EDGE OF L2*/
/*       2=CROSSING INSIDE L2,ON EDGE OF L1*/
/*       3=CROSSING INSIDE L1 AND L2*/
{
  double inner10,inner20,inner12;
  double htlength1,htlength2,otlength;
  double eps1=1.0E-08,eps2=1.0E-08;
  struct vector vctht1,vctht2;  /*VECTOR FROM HEAD TO TAIL OF LINE.*/
  struct vector vcthh;            /*VECTOR FROM HEAD L1 TO HEAD L2.*/
  struct vector outer;                             /*OUTER PRODUCT.*/
  struct onode n1,n2;

  if(iside!=NULL) *iside=0;

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

  if(lengthvector(vcthh)<eps1) /*CASE SAME HEAD.*/
  {
	*nearest1=0.0;
	*nearest2=0.0;
	*distance=0.0;
	return 1;
  }

  if(fabs(innerproduct(vcthh,outer))<eps1)
  /*CASE L1 AND L2 ON SAME PLANE.*/
  {
	if(fabs(outer.dc[GZ])>eps1)
    {
      otlength=outer.dc[GZ];
      htlength1=vcthh.dc[GX]*vctht2.dc[GY]
               -vcthh.dc[GY]*vctht2.dc[GX];
      htlength2=vcthh.dc[GX]*vctht1.dc[GY]
               -vcthh.dc[GY]*vctht1.dc[GX];
    }
	else if(fabs(outer.dc[GX])>eps1)
    {
      otlength=outer.dc[GX];
      htlength1=vcthh.dc[GY]*vctht2.dc[GZ]
               -vcthh.dc[GZ]*vctht2.dc[GY];
      htlength2=vcthh.dc[GY]*vctht1.dc[GZ]
               -vcthh.dc[GZ]*vctht1.dc[GY];
    }
	else if(fabs(outer.dc[GY])>eps1)
    {
      otlength=outer.dc[GY];
      htlength1=vcthh.dc[GZ]*vctht2.dc[GX]
               -vcthh.dc[GX]*vctht2.dc[GZ];
      htlength2=vcthh.dc[GZ]*vctht1.dc[GX]
               -vcthh.dc[GX]*vctht1.dc[GZ];
    }
    else /*PARALLEL CASE OUTER=0*/
    {
      *distance=distancedotline(l1.ends[0],l2,&n2,NULL);
      *nearest1=0.0;
	  if(fabs(vctht2.dc[GX])>eps1)
      {
        *nearest2=(n2.d[GX]-l2.ends[0].d[GX])/vctht2.dc[GX];
      }
	  else if(fabs(vctht2.dc[GY])>eps1)
      {
        *nearest2=(n2.d[GY]-l2.ends[0].d[GY])/vctht2.dc[GY];
      }
	  else if(fabs(vctht2.dc[GZ])>eps1)
      {
        *nearest2=(n2.d[GZ]-l2.ends[0].d[GZ])/vctht2.dc[GZ];
      }
      else /*ERROR CASE VECTOR L2=0.*/
      {
        *nearest2=0.0;
        return 0;
      }
      return 2; /*CASE PARALLEL BUT DISTANCE SUCCEED.*/
    }

    *nearest1=htlength1/otlength;
    *nearest2=htlength2/otlength;
  }
  else /*CASE L1 AND L2 NOT ON SAME PLANE.*/
  {
    otlength=lengthvector(outer);
    htlength1=lengthvector(vctht1);
    htlength2=lengthvector(vctht2);
	if(otlength<eps1 || htlength1<eps1 || htlength2<eps1)
    {
      *nearest1=0.0;
      *nearest2=0.0;
      *distance=0.0;
      return 0;
      /*ERROR CASE PARALLEL OR L1=0 OR L2=0.*/
      /*ERROR IF "NOT ON SAME PLANE" AND "PARALLEL".*/
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

  if(iside!=NULL)
  {
    if(eps2<*nearest1 && *nearest1<1.0-eps2 &&
       eps2<*nearest2 && *nearest2<1.0-eps2)        *iside=3;
    else if(-eps2<=*nearest1 && *nearest1<=eps2 &&
             eps2< *nearest2 && *nearest2<1.0-eps2) *iside=2;
    else if( eps2< *nearest1 && *nearest1<1.0-eps2 &&
            -eps2<=*nearest2 && *nearest2<=eps2)    *iside=1;
    else                                            *iside=0;
  }

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
  double eps=1.0E-04;
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

  if(!distancelineline(l1,l2,&param1,&param2,&distance,NULL))
  {
    /*DISTANCELINELINE RETURN 0 ONLY IN ERROR CASE.*/
    *l0=l1;
    return 0;
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

  /*HEAD=MIDPOINT OF NEAREST DOTS N1,N2.*/
  l0->ends[0].d[GX]=0.5*(n1.d[GX]+n2.d[GX]);
  l0->ends[0].d[GY]=0.5*(n1.d[GY]+n2.d[GY]);
  l0->ends[0].d[GZ]=0.5*(n1.d[GZ]+n2.d[GZ]);

  /*DIRECTION=VECTOR L1+L2.*/
  l0->ends[1].d[GX]=(l0->ends[0].d[GX])
                   +(vctht1.dc[GX]+vctht2.dc[GX]);
  l0->ends[1].d[GY]=(l0->ends[0].d[GY])
                   +(vctht1.dc[GY]+vctht2.dc[GY]);
  l0->ends[1].d[GZ]=(l0->ends[0].d[GZ])
                   +(vctht1.dc[GZ]+vctht2.dc[GZ]);

  /*IF L0 HEAD=TAIL,L1 AND L2 MAKE 180 DEGREE.*/
  /*CAN NOT DECIDE PERPENDICULAR LINE IN 3D SPACE.*/
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

}/*intersectlineline*/

int intersectlinemouse(int ix1,int iy1,int ix2,int iy2,POINT point)
/*INTERSECTION OF LINE ON SCREEN AND BOX MOUSE.*/
{
  int n1,n2,n3,n4;
  int boxsize=BOXMOUSE;

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

  if(n1 || n2 || n3 || n4) return 1;
  else                     return 0;

}/*intersectlinemouse*/

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

  if(fin==NULL) return NULL;

  i=0;
  s=NULL;
  while(1)
  {
	ii=0;
	ss=NULL;
	while((ichr=getc(f))==' ')
	;/*IGNORE SPACE STRINGS AT THE HEAD.*/
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

struct snode *addtext(struct snode *strset,int *nstr,char *str,
                      double tx,double ty)
{
  (*nstr)++;
  strset=(struct snode *)realloc(strset,
                                 (*nstr)*sizeof(struct snode));
  strcpy((strset+(*nstr)-1)->str,str);
  (strset+(*nstr)-1)->n.d[0]=tx;
  (strset+(*nstr)-1)->n.d[1]=ty;
  (strset+(*nstr)-1)->n.d[2]=0.0;
  return strset;
}/*addtext*/

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
    GetAsyncKeyState(VK_RBUTTON);                /*CLEAR KEY RIGHT.*/

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

    /*while(!GetAsyncKeyState(VK_RBUTTON))
    ;*/                               /*RIGHT CLICK TO CONTINUE.*/
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
					 /*ANSI_CHARSET,*/
					 DEFAULT_CHARSET,
                     /*SHIFTJIS_CHARSET,*/    /*CHAR SET*/
					 OUT_DEFAULT_PRECIS,
					 CLIP_DEFAULT_PRECIS,
					 DEFAULT_QUALITY,
					 DEFAULT_PITCH | FF_DONTCARE,
					 fontname);                         /*FONT NAME*/
	pfont=(HFONT)SelectObject(hdc,hfont);
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
  pbit=(HBITMAP)SelectObject(hdc,hbit);
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

  sprintf(str,"%.3f",vp.range.max.d[GX]); /*RANGE DATA.*/
  SetDlgItemText(hdwnd,IDR_XMAX,str);
  sprintf(str,"%.3f",vp.range.min.d[GX]);
  SetDlgItemText(hdwnd,IDR_XMIN,str);
  sprintf(str,"%.3f",vp.range.max.d[GY]);
  SetDlgItemText(hdwnd,IDR_YMAX,str);
  sprintf(str,"%.3f",vp.range.min.d[GY]);
  SetDlgItemText(hdwnd,IDR_YMIN,str);
  sprintf(str,"%.3f",vp.range.max.d[GZ]);
  SetDlgItemText(hdwnd,IDR_ZMAX,str);
  sprintf(str,"%.3f",vp.range.min.d[GZ]);
  SetDlgItemText(hdwnd,IDR_ZMIN,str);

  sprintf(str,"%.1f",vp.dparam.gaxis);
  SetDlgItemText(hdwnd,IDV_GAXISLENGTH,str);
  sprintf(str,"%.1f",vp.dparam.eaxis);
  SetDlgItemText(hdwnd,IDV_EAXISLENGTH,str);
  sprintf(str,"%.1f",vp.dparam.dfact);
  SetDlgItemText(hdwnd,IDV_DFACTOR,str);
  sprintf(str,"%.3f",vp.dparam.qfact);
  SetDlgItemText(hdwnd,IDV_QFACTOR,str);
  sprintf(str,"%.3f",vp.dparam.mfact);
  SetDlgItemText(hdwnd,IDV_MFACTOR,str);
  sprintf(str,"%.1f",vp.dparam.pitch);
  SetDlgItemText(hdwnd,IDV_GYOPITCH,str);
  sprintf(str,"%.0f",vp.dparam.hsize);
  SetDlgItemText(hdwnd,IDV_HINGESIZE,str);

  return;
}/*setviewparam*/

void setviewpoint(HWND hdwnd,
                  struct arclmframe af,struct viewparam *vp)
/*SET VIEW POINT OF ARCLM FRAME.*/
{
  int i,nx,ny,imax[2],imin[2];
  long int wx,wy;
  double dmax[3],dmin[3],fact[4],fmin;
  struct onode on;

  /*RANGE*/
  for(i=0;i<af.nnode;i++)
  {
    if(!nodeontoscreen(*(af.nodes+i),&nx,&ny,*vp)) return;

    if(i==0) /*INITIAL*/
    {
      imax[0]=nx; imin[0]=nx;
      imax[1]=ny; imin[1]=ny;
      dmax[GX]=(af.nodes+i)->d[GX]; dmin[GX]=(af.nodes+i)->d[GX];
      dmax[GY]=(af.nodes+i)->d[GY]; dmin[GY]=(af.nodes+i)->d[GY];
      dmax[GZ]=(af.nodes+i)->d[GZ]; dmin[GZ]=(af.nodes+i)->d[GZ];
    }

    if(dmax[GX]<(af.nodes+i)->d[GX]) dmax[GX]=(af.nodes+i)->d[GX];
    if(dmax[GY]<(af.nodes+i)->d[GY]) dmax[GY]=(af.nodes+i)->d[GY];
    if(dmax[GZ]<(af.nodes+i)->d[GZ]) dmax[GZ]=(af.nodes+i)->d[GZ];
    if(dmin[GX]>(af.nodes+i)->d[GX]) dmin[GX]=(af.nodes+i)->d[GX];
    if(dmin[GY]>(af.nodes+i)->d[GY]) dmin[GY]=(af.nodes+i)->d[GY];
    if(dmin[GZ]>(af.nodes+i)->d[GZ]) dmin[GZ]=(af.nodes+i)->d[GZ];

    if(imax[0]<nx) imax[0]=nx;
    if(imax[1]<ny) imax[1]=ny;
    if(imin[0]>nx) imin[0]=nx;
    if(imin[1]>ny) imin[1]=ny;
  }

  /*SET FOCUS ON CENTER*/
  on.d[0]=0.5*(dmax[0]+dmin[0]);
  on.d[1]=0.5*(dmax[1]+dmin[1]);
  on.d[2]=0.5*(dmax[2]+dmin[2]);

  if(!nodeontoscreen(on,&nx,&ny,*vp)) return;

  vp->focus.d[0]=on.d[0];
  vp->focus.d[1]=on.d[1];
  vp->focus.d[2]=on.d[2];

  /*SET DISTANCE*/
  getclientsize(hdwnd,&wx,&wy);

  if(imax[0]!=nx) fact[0]=((double)wx/2)/(imax[0]-nx);
  if(imin[0]!=nx) fact[1]=((double)wx/2)/(nx-imin[0]);
  if(imax[1]!=ny) fact[2]=((double)wy/2)/(imax[1]-ny);
  if(imin[1]!=ny) fact[3]=((double)wy/2)/(ny-imin[1]);

  if(imax[0]!=nx)      fmin=fact[0];
  else if(imin[0]!=nx) fmin=fact[1];
  else if(imax[0]!=ny) fmin=fact[2];
  else if(imin[0]!=ny) fmin=fact[3];
  else                 fmin=1.0;

  if(imax[0]!=nx && fmin>fact[0]) fmin=fact[0];
  if(imin[0]!=nx && fmin>fact[1]) fmin=fact[1];
  if(imax[1]!=ny && fmin>fact[2]) fmin=fact[2];
  if(imin[1]!=ny && fmin>fact[3]) fmin=fact[3];

  vp->odv*=0.8*fmin;

  vp->range.max.d[GX]=dmax[GX]+100.0;
  vp->range.min.d[GX]=dmin[GX]-100.0;
  vp->range.max.d[GY]=dmax[GY]+100.0;
  vp->range.min.d[GY]=dmin[GY]-100.0;
  vp->range.max.d[GZ]=dmax[GZ]+100.0;
  vp->range.min.d[GZ]=dmin[GZ]-100.0;

  createviewdata(vp);
  return;
}/*setviewpoint*/

void getdoublefromdialog(HWND hdwnd,int iditem,double *d)
/*GET OPAQUE RATE FROM DIALOG.*/
{
  char data[80];

  GetDlgItemText(hdwnd,iditem,data,80);
  *d=strtod(data,NULL);
  return;
}/*getdoublefromdialog*/

void setdoubleintodialog(HWND hdwnd,int iditem,double d)
/*SET OPAQUE RATE INTO DIALOG.*/
{
  char str[80];

  sprintf(str,"%.3f",d);
  SetDlgItemText(hdwnd,iditem,str);
  return;
}/*setdoubleintodialog*/

void getmasshiju(HWND hdwnd,double *hiju)
/*GET HIJU OF MASS FROM DIALOG.*/
{
  char data[80];

  GetDlgItemText(hdwnd,IDVS_MASSHIJU,data,80);
  *hiju=strtod(data,NULL);

  if(*hiju==0.0) *hiju=1.0;
  return;
}/*getmasshiju*/

void getkatakouparam(HWND hdwnd,
                     double *H,double *B,
                     double *tf1,double *tf2,double *tw,
                     double *ri1,double *ri2,
                     double *ri3,double *ri4,
                     double *ro1,double *ro2,
                     double *ro3,double *ro4)
/*GET PARAMETERS OF KATAKOU FROM DIALOG.*/
{
  char data[80];

  GetDlgItemText(hdwnd,IDK_POPHEIGHT,data,80);
  *H=strtod(data,NULL);
  GetDlgItemText(hdwnd,IDK_POPWIDTH,data,80);
  *B=strtod(data,NULL);
  GetDlgItemText(hdwnd,IDK_POPTF1,data,80);
  *tf1=strtod(data,NULL);
  GetDlgItemText(hdwnd,IDK_POPTF2,data,80);
  *tf2=strtod(data,NULL);
  GetDlgItemText(hdwnd,IDK_POPTW,data,80);
  *tw=strtod(data,NULL);

  GetDlgItemText(hdwnd,IDK_POPRI1,data,80);
  *ri1=strtod(data,NULL);
  GetDlgItemText(hdwnd,IDK_POPRI2,data,80);
  *ri2=strtod(data,NULL);
  GetDlgItemText(hdwnd,IDK_POPRI3,data,80);
  *ri3=strtod(data,NULL);
  GetDlgItemText(hdwnd,IDK_POPRI4,data,80);
  *ri4=strtod(data,NULL);

  GetDlgItemText(hdwnd,IDK_POPRO1,data,80);
  *ro1=strtod(data,NULL);
  GetDlgItemText(hdwnd,IDK_POPRO2,data,80);
  *ro2=strtod(data,NULL);
  GetDlgItemText(hdwnd,IDK_POPRO3,data,80);
  *ro3=strtod(data,NULL);
  GetDlgItemText(hdwnd,IDK_POPRO4,data,80);
  *ro4=strtod(data,NULL);

  return;
}/*getkatakouparam*/

void getproperty(HWND hdwnd,
                 long int *code,char *name,
                 double *E,double *poi,double *hiju)
/*GET PROPERTY OF POLYCURVE FROM DIALOG.*/
{
  char data[80];

  GetDlgItemText(hdwnd,IDPP_CODE,data,80);
  *code=(int)strtol(data,NULL,10);

  GetDlgItemText(hdwnd,IDPP_CAPTION,name,80);

  GetDlgItemText(hdwnd,IDPP_E,data,80);
  *E=strtod(data,NULL);
  GetDlgItemText(hdwnd,IDPP_POISSON,data,80);
  *poi=strtod(data,NULL);
  GetDlgItemText(hdwnd,IDPP_HIJU,data,80);
  *hiju=strtod(data,NULL);

  return;
}/*getproperty*/

void setlaps(HWND hdwnd,int ilap,int nlap)
/*SET CURRENT LAP INTO DIALOG.*/
{
  char str[80];

  sprintf(str,"%d/%d",ilap,nlap);
  SetDlgItemText(hdwnd,ID_LAPS,str);
  SendDlgItemMessage(hdwnd,ID_LAPS,WM_PAINT,0,0);

  return;
}/*setlaps*/

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
    *ix=vp.Xo+(int)(vp.gfactor*np.d[GX]+0.5);
    *iy=vp.Yo-(int)(vp.gfactor*np.d[GY]+0.5);
    return 1;
  }
  else
  {
    *ix=0;
    *iy=0;
    return 0;
  }
}/*nodeontoscreen*/

int nodeinsiderect(struct onode on,
                   int il,int ir,int it,int ib,
                   struct viewparam vp)
/*PROJECTED NODE INSIDE OR OUTSIDE RECT.*/
{
  int ix,iy,in;

  if(nodeontoscreen(on,&ix,&iy,vp))
  {
    if(il>ir){in=il; il=ir; ir=in;} /*CORRECT TO LEFT<RIGHT*/
    if(it>ib){in=it; it=ib; ib=in;} /*CORRECT TO TOP<BOTTOM*/

    if(il<=ix && ix<=ir && /*LEFT,RIGHT*/
       it<=iy && iy<=ib)   /*TOP,BOTTOM*/
    {
      return 1;
    }
    else return 0;
  }
  else return 0;

}/*nodeinsiderect*/

void drawglobaltext(HDC hdc,struct viewparam vp,
                    struct onode tn,char *str)
/*DRAW TEXT GLOBAL ON SCREEN BY AXONOMETRICS,PERSPECTIVE.*/
/*TEXTCOLOR DEFINITION MUST BE ALREADY DONE FOR FAST DRAWING.*/
{
  int Ox,Oy;

  if(!nodeontoscreen(tn,&Ox,&Oy,vp)) return;           /*PROJECTION*/
  TextOut(hdc,Ox+2,Oy,str,strlen(str));

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

  TextOut(hdc,Ox+2,Oy,str,strlen(str));

  return;
}/*drawglobaltextaligned*/

void drawglobalnode(HDC hdc,struct viewparam vp,struct onode gn,
                       struct oload *gl)
/*DRAW NODE GLOBAL ON SCREEN BY AXONOMETRICS,PERSPECTIVE.*/
/*TEXTCOLOR DEFINITION MUST BE ALREADY DONE FOR FAST DRAWING.*/
{
  char str[20];
  int Ox,Oy;                                   /*COORDINATION ARROW*/
  SIZE size;

  if(!nodeontoscreen(gn,&Ox,&Oy,vp)) return;           /*PROJECTION*/

  if(vp.vflag.nv.code)
  {
sprintf(str,"%ld",gn.code);
//  sprintf(str,"x=%.3f",gn.d[GX]);
//  sprintf(str,"%.3f",gn.d[GX]);
//TextOut(hdc,Ox+2,Oy,str,strlen(str));
//   sprintf(str,"y=%.3f",gn.d[GY]);
//   sprintf(str,"%.3f",gn.d[GY]);
//TextOut(hdc,Ox+2,Oy-12,str,strlen(str));
//   sprintf(str,"z=%.3f",gn.d[GZ]);                  //mihara
//   sprintf(str,"%.3f",gn.d[GZ]);                  //mihara
//TextOut(hdc,Ox+2,Oy-24,str,strlen(str));

	TextOut(hdc,Ox+2,Oy,str,strlen(str));

//node circle mihara////////////////////////////////////////////////////////////
#if 0
		  HPEN hpen,hpenrange,ppen;
		  HBRUSH hbrush,hbrushrange,pbrush;
		  LOGBRUSH lb;

		  hpen        = (HPEN)GetCurrentObject(hdc,OBJ_PEN);
		  hbrush      = (HBRUSH)GetCurrentObject(hdc,OBJ_BRUSH);
		  hbrushrange = CreateBrushIndirect(&lb);
		  hpenrange   = CreatePen(PS_SOLID,1,RGB(150,150,150));
		  ppen        = (HPEN)SelectObject(hdc,hpenrange);
		  pbrush      = (HBRUSH)SelectObject(hdc,hbrushrange);

		  Ellipse(hdc,(Ox-5),(Oy-5),(Ox+5),(Oy+5));  /*FILLED CIRCLE*/

		  SelectObject(hdc,hbrush);
		  SelectObject(hdc,hpen);
		  DeleteObject(hpenrange);
		  DeleteObject(hbrushrange);
#endif
//node circle mihara////////////////////////////////////////////////////////////

//    sprintf(str,"%.3f",gn.d[GY]);                  //mihara
//    TextOut(hdc,Ox+2,Oy+20,str,strlen(str));

  }
/*mihara for zahyo 20190329///////////////////////////////////////////////////*/
  if(vp.vflag.nv.d[0])
  {
	if(vp.vflag.nv.code)
	{
	  GetTextExtentPoint32(hdc,str,strlen(str),&size);
	  Oy+=size.cy-2;
	}
	sprintf(str,"%.3f",gn.d[GX]);
	TextOut(hdc,Ox+2,Oy,str,strlen(str));
  }
  if(vp.vflag.nv.d[1])
  {
	if(vp.vflag.nv.code||vp.vflag.nv.d[0])
	{
	  GetTextExtentPoint32(hdc,str,strlen(str),&size);
	  Oy+=size.cy-2;
	}
	sprintf(str,"%.3f",gn.d[GY]);
	TextOut(hdc,Ox+2,Oy,str,strlen(str));
  }
  if(vp.vflag.nv.d[2])
  {
	if(vp.vflag.nv.code||vp.vflag.nv.d[0]||vp.vflag.nv.d[1])
	{
	  GetTextExtentPoint32(hdc,str,strlen(str),&size);
	  Oy+=size.cy-2;
	}
	sprintf(str,"%.3f",gn.d[GZ]);
	TextOut(hdc,Ox+2,Oy,str,strlen(str));
  }
/*mihara for zahyo 20190329///////////////////////////////////////////////////*/

  if(gl!=NULL && vp.vflag.nv.mvalue)
  {
	if(vp.vflag.nv.code)
	if(vp.vflag.nv.code||vp.vflag.nv.d[0]||vp.vflag.nv.d[1]||vp.vflag.nv.d[2]) //mihara for zahyo 20190329
	{
	  GetTextExtentPoint32(hdc,str,strlen(str),&size);
	  Oy+=size.cy-2;
	}
	sprintf(str,"%.3f",globalunit*(gl+(gn.loff))->w[WEIGHTSLAB]);
	TextOut(hdc,Ox+2,Oy,str,strlen(str));
  }
  return;
}/*drawglobalnode*/
#if 0
void drawglobalnode(HDC hdc,struct viewparam vp,
					struct onode gn,struct oload *gl)
/*DRAW NODE GLOBAL ON SCREEN BY AXONOMETRICS,PERSPECTIVE.*/
/*TEXTCOLOR DEFINITION MUST BE ALREADY DONE FOR FAST DRAWING.*/
{
  char str[20];
  int Ox,Oy;                                   /*COORDINATION ARROW*/
  SIZE size;

  if(!nodeontoscreen(gn,&Ox,&Oy,vp)) return;           /*PROJECTION*/

  if(vp.vflag.nv.code)
  {
	sprintf(str,"%ld",gn.code);
	TextOut(hdc,Ox+2,Oy,str,strlen(str));

	if(globalflag==1)
	{
	  sprintf(str,"X=%.3f",gn.d[GX]);
	  TextOut(hdc,Ox+2,Oy+14,str,strlen(str));
	  sprintf(str,"Y=%.3f",gn.d[GY]);
	  TextOut(hdc,Ox+2,Oy+28,str,strlen(str));
	  sprintf(str,"Z=%.3f",gn.d[GZ]);
	  TextOut(hdc,Ox+2,Oy+42,str,strlen(str));
	}
  }
/*mihara for zahyo 20190329///////////////////////////////////////////////////*/
  if(vp.vflag.nv.d[0])
  {
	if(vp.vflag.nv.code)
	{
	  GetTextExtentPoint32(hdc,str,strlen(str),&size);
	  Oy+=size.cy-2;
	}
	sprintf(str,"%.3f",gn.d[GX]);
	TextOut(hdc,Ox+2,Oy,str,strlen(str));
  }
  if(vp.vflag.nv.d[1])
  {
	if(vp.vflag.nv.code||vp.vflag.nv.d[0])
	{
	  GetTextExtentPoint32(hdc,str,strlen(str),&size);
	  Oy+=size.cy-2;
	}
	sprintf(str,"%.3f",gn.d[GY]);
	TextOut(hdc,Ox+2,Oy,str,strlen(str));
  }
  if(vp.vflag.nv.d[2])
  {
	if(vp.vflag.nv.code||vp.vflag.nv.d[0]||vp.vflag.nv.d[1])
	{
	  GetTextExtentPoint32(hdc,str,strlen(str),&size);
	  Oy+=size.cy-2;
	}
	sprintf(str,"%.3f",gn.d[GZ]);
	TextOut(hdc,Ox+2,Oy,str,strlen(str));
  }
/*mihara for zahyo 20190329///////////////////////////////////////////////////*/
  if(gl!=NULL && vp.vflag.nv.mvalue)
  {
    if(vp.vflag.nv.code)
	if(vp.vflag.nv.code||vp.vflag.nv.d[0]||vp.vflag.nv.d[1]||vp.vflag.nv.d[2]) //mihara for zahyo 20190329
	{
	  GetTextExtentPoint32(hdc,str,strlen(str),&size);
	  Oy+=size.cy-2;
	}
	sprintf(str,"%.3f",globalunit*(gl+(gn.loff))->w[WEIGHTSLAB]);
	TextOut(hdc,Ox+2,Oy,str,strlen(str));
  }

  return;
}/*drawglobalnode*/
#endif
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

void drawprintrange(HDC hdc,struct viewparam vp)
{
          long int maxX,maxY;
          int cWidthPels,cHeightPels;
          int cWidthPelsb,cHeightPelsb;
          char   str[256];
          HPEN hpen,hpenrange,ppen;
          HBRUSH hbrush,hbrushrange,pbrush;
          LOGBRUSH lb;

          lb.lbStyle = BS_NULL;

		  hpen        = (HPEN)GetCurrentObject(hdc,OBJ_PEN);
          hbrush      = (HBRUSH)GetCurrentObject(hdc,OBJ_BRUSH);
          //hbrush    = CreateSolidBrush(RGB(255,255,255));
          hbrushrange = CreateBrushIndirect(&lb);
          hpenrange   = CreatePen(PS_SOLID,1,RGB(150,150,150));
          ppen        = (HPEN)SelectObject(hdc,hpenrange);
		  pbrush      = (HBRUSH)SelectObject(hdc,hbrushrange);


          getclientsize((wdraw.childs+1)->hwnd,&maxX,&maxY);

          cWidthPelsb  = maxX;
          cHeightPelsb = maxY;

          vp.Xo = 0.5*cWidthPelsb;
          vp.Yo = 0.5*cHeightPelsb;

          if((wdraw.childs+1)->vparam.vflag.pv.a4tate==1)
          {
//MessageBox(NULL,"A4Tate","PrintRange",MB_OK);
          cWidthPels  = 0.5*496.062;     //forA4TATE
          cHeightPels = 0.5*701.15734;   //forA4TATE
          }

          else if((wdraw.childs+1)->vparam.vflag.pv.a4yoko==1)
          {
//MessageBox(NULL,"A4Yoko","PrintRange",MB_OK);
          cWidthPels  = 0.5*701.15734;   //forA4YOKO
          cHeightPels = 0.5*496.062;     //forA4YOKO
          }

          else if((wdraw.childs+1)->vparam.vflag.pv.a3tate==1)
          {
//MessageBox(NULL,"A3Tate","PrintRange",MB_OK);
          cWidthPels  = 0.5*701.15734;   //forA3TATE
          cHeightPels = 0.5*992.124;     //forA3TATE
          }

          else if((wdraw.childs+1)->vparam.vflag.pv.a3yoko==1)
          {
//MessageBox(NULL,"A3Yoko","PrintRange",MB_OK);
          cWidthPels  = 0.5*992.124;     //forA3YOKO
          cHeightPels = 0.5*701.15734;   //forA3YOKO
          }

          else
          {
          MessageBox(NULL,"Select Paper Size","PrintRange",MB_OK);
          cWidthPels  = 0.5*496.062;     //forA4TATE
          cHeightPels = 0.5*701.15734;   //forA4TATE
          }

          Rectangle(hdc,vp.Xo-cWidthPels,vp.Yo-cHeightPels,
                        vp.Xo+cWidthPels,vp.Yo+cHeightPels);

          SelectObject(hdc,hbrush);
          SelectObject(hdc,hpen);
          DeleteObject(hpenrange);
          DeleteObject(hbrushrange);

          return;
}/*drawprintrange*/

void drawglobalconfig(HDC hdc,struct viewparam vp,struct organ go,
                                struct onode gn,int mode,long int loff)
/*DRAW CONFFIG GLOBAL ON SCREEN BY AXONOMETRICS,PERSPECTIVE.*/
{
  int    Ox,Oy;                                   /*COORDINATION ARROW*/
  POINT  ptLine[4];
  HPEN   hpen,hpenpin,hpenrigid,hpenprt,ppen;
  HBRUSH hbrush,hbrushpin,hbrushrigid,hbrushprt,pbrush;
  char   str[256];

  if(!nodeontoscreen(gn,&Ox,&Oy,vp)) return;           /*PROJECTION*/

  hpen   = (HPEN)GetCurrentObject(hdc,OBJ_PEN);
  hbrush = (HBRUSH)GetCurrentObject(hdc,OBJ_BRUSH);

  hpenpin   = CreatePen(PS_SOLID,1,RGB(150,255,0));     /*X,Y,Z FIXED*/
  hbrushpin = CreateSolidBrush(RGB(150,255,0));

  hpenrigid   = CreatePen(PS_SOLID,1,RGB(255,0,150));   /*Tx,Ty,Tz FIXED*/
  hbrushrigid = CreateSolidBrush(RGB(255,0,150));

  hpenprt   = CreatePen(PS_SOLID,1,RGB(0,0,0));         /*BLACK*/
  hbrushprt = CreateSolidBrush(RGB(120,120,120));       /*GRAY*/

  if((go.confs+(loff+3))->iconf==1 ||
     (go.confs+(loff+4))->iconf==1 ||
     (go.confs+(loff+5))->iconf==1) /*Tx,Ty,Tz FIXED*/
  {
    if(mode==ONPRINTER||mode==ONPREVIEW)
    {
	  ppen=(HPEN)SelectObject(hdc,hpenprt);
      pbrush=(HBRUSH)SelectObject(hdc,hbrushprt);
    }
    else
    {
      ppen=(HPEN)SelectObject(hdc,hpenrigid);
      pbrush=(HBRUSH)SelectObject(hdc,hbrushrigid);
    }

	ptLine[0].x = Ox-(vp.dparam.hsize)*sqrt(3.0);  ptLine[0].y = Oy;
	ptLine[1].x = Ox+(vp.dparam.hsize)*sqrt(3.0);  ptLine[1].y = Oy;
	ptLine[2].x = Ox+(vp.dparam.hsize)*sqrt(3.0);  ptLine[2].y = Oy+2.0*(vp.dparam.hsize)*sqrt(3.0);
	ptLine[3].x = Ox-(vp.dparam.hsize)*sqrt(3.0);  ptLine[3].y = Oy+2.0*(vp.dparam.hsize)*sqrt(3.0);

    Polygon(hdc,ptLine,4);

  }
  else if((go.confs+(loff+0))->iconf==1 ||
          (go.confs+(loff+1))->iconf==1 ||
          (go.confs+(loff+2))->iconf==1) /*X,Y,Z FIXED*/
  {
    if(mode==ONPRINTER||mode==ONPREVIEW)
    {
	  ppen=(HPEN)SelectObject(hdc,hpenprt);
      pbrush=(HBRUSH)SelectObject(hdc,hbrushprt);
    }
    else
    {
      ppen=(HPEN)SelectObject(hdc,hpenpin);
      pbrush=(HBRUSH)SelectObject(hdc,hbrushpin);
    }

    ptLine[0].x = Ox;                        ptLine[0].y = Oy;
	ptLine[1].x = Ox+2.0*(vp.dparam.hsize);  ptLine[1].y = Oy+2.0*(vp.dparam.hsize)*sqrt(3.0);
	ptLine[2].x = Ox-2.0*(vp.dparam.hsize);  ptLine[2].y = Oy+2.0*(vp.dparam.hsize)*sqrt(3.0);

//sprintf(str,"ConfSize=%.0f",vp.dparam.csize);
//MessageBox(NULL,str,"Conffig",MB_OK);

    Polygon(hdc,ptLine,3);

  }

//  if((go.confs+(loff+0))->iconf!=1 &&
//          (go.confs+(loff+1))->iconf!=1 &&
//          (go.confs+(loff+2))->iconf==1 &&
//          (gn.d[GZ]!=0.0) && (gn.d[GX]<=20.0) && (gn.d[GX]>=-20.0)) /*X,Y,Z FIXED*/    //Toyota
//  {
//    if(mode==ONPRINTER||mode==ONPREVIEW)
//    {
//      ppen=(HPEN)SelectObject(hdc,hpenprt);
//      pbrush=(HBRUSH)SelectObject(hdc,hbrushprt);
//    }
//    else
//    {
//      ppen=(HPEN)SelectObject(hdc,hpenpin);
//      pbrush=(HBRUSH)SelectObject(hdc,hbrushpin);
//    }

//    ptLine[0].x = Ox;                        ptLine[0].y = Oy;
//    ptLine[1].x = Ox+2.0*(vp.dparam.hsize);  ptLine[1].y = Oy+2.0*(vp.dparam.hsize)*sqrt(3);
//    ptLine[2].x = Ox-2.0*(vp.dparam.hsize);  ptLine[2].y = Oy+2.0*(vp.dparam.hsize)*sqrt(3);

//sprintf(str,"ConfSize=%.0f",vp.dparam.csize);
//MessageBox(NULL,str,"Conffig",MB_OK);

//    Polygon(hdc,ptLine,3);

//  }

  SelectObject(hdc,hpen);
  SelectObject(hdc,hbrush);
  DeleteObject(hpenpin);
  DeleteObject(hbrushpin);
  DeleteObject(hpenrigid);
  DeleteObject(hbrushrigid);
  DeleteObject(hpenprt);
  DeleteObject(hbrushprt);

//  free(ptLine);

  return;
}/*drawglobalconfig*/

////////081209 araki for conffing////////////////////////////////
void drawglobalconfigofaf(HDC hdc,struct viewparam vp,struct arclmframe go,
                                struct onode gn,int mode,long int loff)
/*DRAW CONFFIG GLOBAL ON SCREEN BY AXONOMETRICS,PERSPECTIVE.*/
{
  int    Ox,Oy;                                   /*COORDINATION ARROW*/
  POINT  ptLine[4];
  HPEN   hpen,hpenpin,hpenrigid,hpenprt,ppen;
  HBRUSH hbrush,hbrushpin,hbrushrigid,hbrushprt,pbrush;
  char   str[256];

  if(!nodeontoscreen(gn,&Ox,&Oy,vp)) return;           /*PROJECTION*/

  hpen   = (HPEN)GetCurrentObject(hdc,OBJ_PEN);
  hbrush = (HBRUSH)GetCurrentObject(hdc,OBJ_BRUSH);

  hpenpin   = CreatePen(PS_SOLID,1,RGB(150,255,0));     /*X,Y,Z FIXED*/
  hbrushpin = CreateSolidBrush(RGB(150,255,0));

  hpenrigid   = CreatePen(PS_SOLID,1,RGB(255,0,150));   /*Tx,Ty,Tz FIXED*/
  hbrushrigid = CreateSolidBrush(RGB(255,0,150));

  hpenprt   = CreatePen(PS_SOLID,1,RGB(0,0,0));         /*BLACK*/
  hbrushprt = CreateSolidBrush(RGB(120,120,120));       /*GRAY*/

  if((go.confs+loff+3)->iconf==1 ||
     (go.confs+loff+4)->iconf==1 ||
     (go.confs+loff+5)->iconf==1) /*Tx,Ty,Tz FIXED*/
  {
    if(mode==ONPRINTER||mode==ONPREVIEW)
    {
	  ppen=(HPEN)SelectObject(hdc,hpenprt);
      pbrush=(HBRUSH)SelectObject(hdc,hbrushprt);
    }
    else
    {
	  ppen=(HPEN)SelectObject(hdc,hpenrigid);
      pbrush=(HBRUSH)SelectObject(hdc,hbrushrigid);
    }

    ptLine[0].x = Ox-(vp.dparam.hsize)*sqrt(3.0);  ptLine[0].y = Oy;
    ptLine[1].x = Ox+(vp.dparam.hsize)*sqrt(3.0);  ptLine[1].y = Oy;
    ptLine[2].x = Ox+(vp.dparam.hsize)*sqrt(3.0);  ptLine[2].y = Oy+2.0*(vp.dparam.hsize)*sqrt(3.0);
    ptLine[3].x = Ox-(vp.dparam.hsize)*sqrt(3.0);  ptLine[3].y = Oy+2.0*(vp.dparam.hsize)*sqrt(3.0);

    Polygon(hdc,ptLine,4);

  }
  else if((go.confs+loff+0)->iconf==1 ||
          (go.confs+loff+1)->iconf==1 ||
          (go.confs+loff+2)->iconf==1) /*X,Y,Z FIXED*/
  {
    if(mode==ONPRINTER||mode==ONPREVIEW)
    {
	  ppen=(HPEN)SelectObject(hdc,hpenprt);
	  pbrush=(HBRUSH)SelectObject(hdc,hbrushprt);
    }
    else
    {
	  ppen=(HPEN)SelectObject(hdc,hpenpin);
      pbrush=(HBRUSH)SelectObject(hdc,hbrushpin);
    }

    ptLine[0].x = Ox;                        ptLine[0].y = Oy;
    ptLine[1].x = Ox+2.0*(vp.dparam.hsize);  ptLine[1].y = Oy+2.0*(vp.dparam.hsize)*sqrt(3.0);
    ptLine[2].x = Ox-2.0*(vp.dparam.hsize);  ptLine[2].y = Oy+2.0*(vp.dparam.hsize)*sqrt(3.0);

//sprintf(str,"ConfSize=%.0f",vp.dparam.csize);
//MessageBox(NULL,str,"Conffig",MB_OK);

    Polygon(hdc,ptLine,3);

  }

//  if((go.confs+(loff+0))->iconf!=1 &&
//          (go.confs+(loff+1))->iconf!=1 &&
//          (go.confs+(loff+2))->iconf==1 &&
//          (gn.d[GZ]!=0.0) && (gn.d[GX]<=20.0) && (gn.d[GX]>=-20.0)) /*X,Y,Z FIXED*/    //Toyota
//  {
//    if(mode==ONPRINTER||mode==ONPREVIEW)
//    {
//      ppen=SelectObject(hdc,hpenprt);
//      pbrush=SelectObject(hdc,hbrushprt);
//    }
//    else
//    {
//      ppen=SelectObject(hdc,hpenpin);
//      pbrush=SelectObject(hdc,hbrushpin);
//    }

//    ptLine[0].x = Ox;                        ptLine[0].y = Oy;
//    ptLine[1].x = Ox+2.0*(vp.dparam.hsize);  ptLine[1].y = Oy+2.0*(vp.dparam.hsize)*sqrt(3);
//    ptLine[2].x = Ox-2.0*(vp.dparam.hsize);  ptLine[2].y = Oy+2.0*(vp.dparam.hsize)*sqrt(3);

//sprintf(str,"ConfSize=%.0f",vp.dparam.csize);
//MessageBox(NULL,str,"Conffig",MB_OK);

//    Polygon(hdc,ptLine,3);

//  }

  SelectObject(hdc,hpen);
  SelectObject(hdc,hbrush);
  DeleteObject(hpenpin);
  DeleteObject(hbrushpin);
  DeleteObject(hpenrigid);
  DeleteObject(hbrushrigid);
  DeleteObject(hpenprt);
  DeleteObject(hbrushprt);

//  free(ptLine);

  return;
}/*drawglobalconfigofaf*/
/////////////////////////////////////////////////////////

double angleofgloballine(struct viewparam vp,struct onode gn1,
                                             struct onode gn2)
/*ANGLE OF GLOBAL LINE ON SCREEN BY AXONOMETRICS,PERSPECTIVE.*/
{
  int Hx,Hy,Tx,Ty;
  double angle; /*[DEGREE]*/

  if(!nodeontoscreen(gn1,&Hx,&Hy,vp)) return 0.0;      /*PROJECTION*/
  if(!nodeontoscreen(gn2,&Tx,&Ty,vp)) return 0.0;      /*PROJECTION*/

  if(Hx==Tx && Hy==Ty)     angle=  0.0;
  else if(Hx==Tx && Hy<Ty) angle=-90.0;
  else if(Hx==Tx && Hy>Ty) angle= 90.0;
  else angle=atan(-(double)(Ty-Hy)/(double)(Tx-Hx))*180.0/PI;

  return angle;
}/*angleofgloballine*/

void drawlinestruct(HDC hdc,struct viewparam vp,struct line l)
/*DRAW LINE STRUCT ON SCREEN BY AXONOMETRICS,PERSPECTIVE.*/
{
  HPEN hpen,ppen;

  hpen=CreatePen(PS_SOLID,1,RGB(l.r,l.g,l.b));
  ppen=(HPEN)SelectObject(hdc,hpen);

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

void sortlines(struct line *pl,int nline,double hugo)
/*SORT LINES BY "HEAP SORT".*/
/*HUGO=1.0:KOUJUN HUGO=-1.0:SHOUJUN*/
{
  int i,j,ir,m;
  double rra;
  struct line rral;

  if(nline<2) return;

  m=(int)(nline/2);
  ir=nline-1;

  for(;;)
  {
    if(m>0)
    {
      m--;
      rra=hugo*((pl+m)->ends[0].d[GX]);
      rral=*(pl+m);
    }
    else
    {
      rra=hugo*((pl+ir)->ends[0].d[GX]);
      rral=*(pl+ir);

      *(pl+ir)=*(pl+0);

      if(--ir==0)
      {
        *(pl+0)=rral;
        return;
      }
    }
    i=m;
    j=m*2+1;
    while(j<=ir)
    {
      if(j<ir &&
         (hugo*((pl+j)->ends[0].d[GX]))>
         (hugo*((pl+j+1)->ends[0].d[GX]))) ++j;

      if(rra>(hugo*((pl+j)->ends[0].d[GX])))
      {
        *(pl+i)=*(pl+j);

        i=j;
        j+=j+1;
      }
      else j=ir+1;
    }
    *(pl+i)=rral;
  }
}/*sortlines*/

int dividearc(struct curve *cv,
              double *dx,double *dy,double *angle,
              struct curve *dcv)
/*DIVIDE ARC BY DOT OR X OR Y OR ANGLE.*/
/*RETURN:NUMBER DIVIDED.*/
/*ARRAY DCV MALLOCED AS 3 ELEMENTS.*/
{
  double r,x0,y0,x1,y1,x2,y2,a0,a1,a2,cx,cy;
  double eps=1.0E-08;

  (dcv+0)->type=cv->type;
  (dcv+1)->type=cv->type;
  (dcv+2)->type=cv->type;
  (dcv+0)->hugo=cv->hugo;
  (dcv+1)->hugo=cv->hugo;
  (dcv+2)->hugo=cv->hugo;
  *((dcv+0)->center)=*(cv->center);
  *((dcv+1)->center)=*(cv->center);
  *((dcv+2)->center)=*(cv->center);
  (dcv+0)->radius[0]=cv->radius[0];
  (dcv+1)->radius[0]=cv->radius[0];
  (dcv+2)->radius[0]=cv->radius[0];

  (dcv+0)->angle[0]=cv->angle[0];
  (dcv+0)->angle[1]=cv->angle[1];
  *((dcv+0)->dots[0])=*(cv->dots[0]);
  *((dcv+0)->dots[1])=*(cv->dots[1]);

  if(angle!=NULL)
  {
    if((cv->hugo==1 && (*angle)>(cv->angle[0]+eps) &&
                       (*angle)<(cv->angle[1]-eps)) ||
       (cv->hugo==-1 && (*angle)<(cv->angle[0]-eps) &&
                        (*angle)>(cv->angle[1]+eps))
      )
    {
      (dcv+0)->angle[0]=cv->angle[0];
      (dcv+0)->angle[1]=(*angle);
      (dcv+1)->angle[0]=(*angle);
	  (dcv+1)->angle[1]=cv->angle[1];

      *((dcv+0)->dots[0])=*(cv->dots[0]);
      (dcv+0)->dots[1]->d[GX]
      =(cv->center->d[GX])+(cv->radius[0])*cos(*angle);
      (dcv+0)->dots[1]->d[GY]
	  =(cv->center->d[GY])+(cv->radius[0])*sin(*angle);
      (dcv+0)->dots[1]->d[GZ]=0.0;

      *((dcv+1)->dots[0])=*((dcv+0)->dots[1]);
      *((dcv+1)->dots[1])=*(cv->dots[1]);

      return 2;
    }
  }
  else if(dx!=NULL && dy==NULL)
  {
    r=cv->radius[0];

    cx=cv->center->d[GX];
    cy=cv->center->d[GY];

    if(r<fabs((*dx)-cx)) return 1;

    y1=-sqrt(r*r-((*dx)-cx)*((*dx)-cx))+cy;
    y2= sqrt(r*r-((*dx)-cx)*((*dx)-cx))+cy;

    a1=definecircleangle(cv,1,(*dx),y1);
    a2=definecircleangle(cv,1,(*dx),y2);

    if(a1>a2){a0=a1; a1=a2; a2=a0;
              y0=y1; y1=y2; y2=y0;} /*A1<A2*/

    if((a1!=a2 && cv->hugo==1 &&
        a1>(cv->angle[0]+eps) && a1<(cv->angle[1]-eps) &&
        a2>(cv->angle[0]+eps) && a2<(cv->angle[1]-eps)) ||
       (a1!=a2 && cv->hugo==-1 &&
        a1<(cv->angle[0]-eps) && a1>(cv->angle[1]+eps) &&
        a2<(cv->angle[0]-eps) && a2>(cv->angle[1]+eps)))
    {
      if(cv->hugo==-1){a0=a1; a1=a2; a2=a0;
                       y0=y1; y1=y2; y2=y0;} /*A1>A2*/

      (dcv+0)->angle[0]=cv->angle[0];
      (dcv+0)->angle[1]=a1;
      (dcv+1)->angle[0]=a1;
      (dcv+1)->angle[1]=a2;
      (dcv+2)->angle[0]=a2;
      (dcv+2)->angle[1]=cv->angle[1];

      *((dcv+0)->dots[0])=*(cv->dots[0]);
      (dcv+0)->dots[1]->d[GX]=(*dx);
      (dcv+0)->dots[1]->d[GY]=y1;
      (dcv+0)->dots[1]->d[GZ]=0.0;
      *((dcv+1)->dots[0])=*((dcv+0)->dots[1]);
      (dcv+1)->dots[1]->d[GX]=(*dx);
      (dcv+1)->dots[1]->d[GY]=y2;
      (dcv+1)->dots[1]->d[GZ]=0.0;
      *((dcv+2)->dots[0])=*((dcv+1)->dots[1]);
      *((dcv+2)->dots[1])=*(cv->dots[1]);

      return 3;
    }
    else if((cv->hugo==1 && a1>(cv->angle[0]+eps)
                         && a1<(cv->angle[1]-eps)) ||
            (cv->hugo==-1 && a1<(cv->angle[0]-eps)
                          && a1>(cv->angle[1]+eps)))
    {
      (dcv+0)->angle[0]=cv->angle[0];
      (dcv+0)->angle[1]=a1;
      (dcv+1)->angle[0]=a1;
      (dcv+1)->angle[1]=cv->angle[1];

      *((dcv+0)->dots[0])=*(cv->dots[0]);
      (dcv+0)->dots[1]->d[GX]=(*dx);
      (dcv+0)->dots[1]->d[GY]=y1;
      (dcv+0)->dots[1]->d[GZ]=0.0;
      *((dcv+1)->dots[0])=*((dcv+0)->dots[1]);
      *((dcv+1)->dots[1])=*(cv->dots[1]);

      return 2;
    }
    else if((cv->hugo==1 && a2>(cv->angle[0]+eps)
                         && a2<(cv->angle[1]-eps)) ||
            (cv->hugo==-1 && a2<(cv->angle[0]-eps)
                          && a2>(cv->angle[1]+eps)))
    {
      (dcv+0)->angle[0]=cv->angle[0];
      (dcv+0)->angle[1]=a2;
      (dcv+1)->angle[0]=a2;
      (dcv+1)->angle[1]=cv->angle[1];

      *((dcv+0)->dots[0])=*(cv->dots[0]);
      (dcv+0)->dots[1]->d[GX]=(*dx);
      (dcv+0)->dots[1]->d[GY]=y2;
      (dcv+0)->dots[1]->d[GZ]=0.0;
      *((dcv+1)->dots[0])=*((dcv+0)->dots[1]);
      *((dcv+1)->dots[1])=*(cv->dots[1]);

      return 2;
    }
  }
  else if(dx==NULL && dy!=NULL)
  {
    r=cv->radius[0];

    cx=cv->center->d[GX];
    cy=cv->center->d[GY];

    if(r<fabs((*dy)-cy)) return 1;

    x1=-sqrt(r*r-((*dy)-cy)*((*dy)-cy))+cx;
    x2= sqrt(r*r-((*dy)-cy)*((*dy)-cy))+cx;

    a1=definecircleangle(cv,1,x1,(*dy));
    a2=definecircleangle(cv,1,x2,(*dy));
    if(a1>a2){a0=a1; a1=a2; a2=a0;
              x0=x1; x1=x2; x2=x0;} /*A1<A2*/

    if((a1!=a2 && cv->hugo==1 &&
        a1>(cv->angle[0]+eps) && a1<(cv->angle[1]-eps) &&
        a2>(cv->angle[0]+eps) && a2<(cv->angle[1]-eps)) ||
       (a1!=a2 && cv->hugo==-1 &&
        a1<(cv->angle[0]-eps) && a1>(cv->angle[1]+eps) &&
        a2<(cv->angle[0]-eps) && a2>(cv->angle[1]+eps)))
    {
      if(cv->hugo==-1){a0=a1; a1=a2; a2=a0;
                       x0=x1; x1=x2; x2=x0;} /*A1>A2*/

      (dcv+0)->angle[0]=cv->angle[0];
      (dcv+0)->angle[1]=a1;
      (dcv+1)->angle[1]=a1;
      (dcv+1)->angle[0]=a2;
      (dcv+2)->angle[0]=a2;
      (dcv+2)->angle[1]=cv->angle[1];

      *((dcv+0)->dots[0])=*(cv->dots[0]);
      (dcv+0)->dots[1]->d[GX]=x1;
      (dcv+0)->dots[1]->d[GY]=(*dy);
      (dcv+0)->dots[1]->d[GZ]=0.0;
      *((dcv+1)->dots[0])=*((dcv+0)->dots[1]);
      (dcv+1)->dots[1]->d[GX]=x2;
      (dcv+1)->dots[1]->d[GY]=(*dy);
      (dcv+1)->dots[1]->d[GZ]=0.0;
      *((dcv+2)->dots[0])=*((dcv+1)->dots[1]);
      *((dcv+2)->dots[1])=*(cv->dots[1]);

      return 3;
    }
    else if((cv->hugo==1 && a1>(cv->angle[0]+eps)
                         && a1<(cv->angle[1]-eps)) ||
            (cv->hugo==-1 && a1<(cv->angle[0]-eps)
                          && a1>(cv->angle[1]+eps)))
    {
      (dcv+0)->angle[0]=cv->angle[0];
      (dcv+0)->angle[1]=a1;
      (dcv+1)->angle[0]=a1;
      (dcv+1)->angle[1]=cv->angle[1];

      *((dcv+0)->dots[0])=*(cv->dots[0]);
      (dcv+0)->dots[1]->d[GX]=x1;
      (dcv+0)->dots[1]->d[GY]=(*dy);
      (dcv+0)->dots[1]->d[GZ]=0.0;
      *((dcv+1)->dots[0])=*((dcv+0)->dots[1]);
      *((dcv+1)->dots[1])=*(cv->dots[1]);

      return 2;
    }
    else if((cv->hugo==1 && a2>(cv->angle[0]+eps)
                         && a2<(cv->angle[1]-eps)) ||
            (cv->hugo==-1 && a2<(cv->angle[0]-eps)
                          && a2>(cv->angle[1]+eps)))
    {
      (dcv+0)->angle[0]=cv->angle[0];
      (dcv+0)->angle[1]=a2;
      (dcv+1)->angle[0]=a2;
      (dcv+1)->angle[1]=cv->angle[1];

      *((dcv+0)->dots[0])=*(cv->dots[0]);
      (dcv+0)->dots[1]->d[GX]=x2;
      (dcv+0)->dots[1]->d[GY]=(*dy);
      (dcv+0)->dots[1]->d[GZ]=0.0;
      *((dcv+1)->dots[0])=*((dcv+0)->dots[1]);
      *((dcv+1)->dots[1])=*(cv->dots[1]);

      return 2;
    }
  }
  else if(dx!=NULL && dy!=NULL)
  {
    a1=definecircleangle(cv,1,(*dx),(*dy));

    if((cv->hugo==1 && a1>(cv->angle[0]+eps)
                    && a1<(cv->angle[1]-eps)) ||
       (cv->hugo==-1 && a1<(cv->angle[0]-eps)
                     && a1>(cv->angle[1]+eps)))
    {
      (dcv+0)->angle[0]=cv->angle[0];
      (dcv+0)->angle[1]=a1;
      (dcv+1)->angle[0]=a1;
      (dcv+1)->angle[1]=cv->angle[1];

      *((dcv+0)->dots[0])=*(cv->dots[0]);
      (dcv+0)->dots[1]->d[GX]=(*dx);
      (dcv+0)->dots[1]->d[GY]=(*dy);
      (dcv+0)->dots[1]->d[GZ]=0.0;
      *((dcv+1)->dots[0])=*((dcv+0)->dots[1]);
      *((dcv+1)->dots[1])=*(cv->dots[1]);

      return 2;
    }
  }

  return 1;
}/*dividearc*/

int dividearcbyline(struct curve *cv,struct line l,
                    struct curve *dcv)
/*DIVIDE ARC BY LINE.*/
/*RETURN:NUMBER DIVIDED.*/
/*ARRAY DCV MALLOCED AS 3 ELEMENTS.*/
{
  int i,ndiv,mdiv;
  double r,cx,cy,x1,y1,x2,y2,xa,ya,xb,yb;
  double aa,ab;
  double Do,fa,fb;
  double eps=1.0E-08;
  struct curve *dc;

  r=cv->radius[0];

  cx=cv->center->d[GX];
  cy=cv->center->d[GY];

  x1=l.ends[0].d[GX];
  y1=l.ends[0].d[GY];
  x2=l.ends[1].d[GX];
  y2=l.ends[1].d[GY];

  (dcv+0)->type=cv->type;
  (dcv+1)->type=cv->type;
  (dcv+2)->type=cv->type;
  (dcv+0)->hugo=cv->hugo;
  (dcv+1)->hugo=cv->hugo;
  (dcv+2)->hugo=cv->hugo;
  *((dcv+0)->center)=*(cv->center);
  *((dcv+1)->center)=*(cv->center);
  *((dcv+2)->center)=*(cv->center);
  (dcv+0)->radius[0]=cv->radius[0];
  (dcv+1)->radius[0]=cv->radius[0];
  (dcv+2)->radius[0]=cv->radius[0];

  (dcv+0)->angle[0]=cv->angle[0];
  (dcv+0)->angle[1]=cv->angle[1];
  *((dcv+0)->dots[0])=*(cv->dots[0]);
  *((dcv+0)->dots[1])=*(cv->dots[1]);

  ndiv=1;

  /*DOT1 AS ORIGIN.*/
  cx-=x1;
  cy-=y1;
  x2-=x1;
  y2-=y1;
  if(x2==0.0 && y2==0.0) return 1; /*LINE LENGTH=0*/

  Do=r*r*(x2*x2+y2*y2)-(cx*y2-x2*cy)*(cx*y2-x2*cy);
  if(Do<eps) return 1;

  fa=((cx*x2+cy*y2)-sqrt(Do))/(x2*x2+y2*y2);
  fb=((cx*x2+cy*y2)+sqrt(Do))/(x2*x2+y2*y2);

  dc=(struct curve *)malloc(3*sizeof(struct curve));
  for(i=0;i<3;i++)
  {
    (dc+i)->center=(struct onode *)malloc(sizeof(struct onode));
    (dc+i)->dots[0]=(struct onode *)malloc(sizeof(struct onode));
    (dc+i)->dots[1]=(struct onode *)malloc(sizeof(struct onode));
  }

  /*DIVISION.RELEASE DOT1=ORIGIN.*/
  if(0.0<=fa && fa<=1.0)
  {
    xa=fa*x2+x1;
    ya=fa*y2+y1;
    aa=definecircleangle(cv,1,xa,ya);

    mdiv=dividearc(cv,NULL,NULL,&aa,dc);

    if(mdiv==2)
    {
      for(i=0;i<mdiv;i++)
      {
        (dcv+i)->angle[0]=(dc+i)->angle[0];
        (dcv+i)->angle[1]=(dc+i)->angle[1];
        *((dcv+i)->dots[0])=*((dc+i)->dots[0]);
        *((dcv+i)->dots[1])=*((dc+i)->dots[1]);
      }
      ndiv++;
    }
  }

  if(0.0<=fb && fb<=1.0)
  {
    xb=fb*x2+x1;
    yb=fb*y2+y1;

    for(i=0;i<ndiv;i++)
    {
      ab=definecircleangle((dcv+i),1,xb,yb);

      mdiv=dividearc((dcv+i),NULL,NULL,&ab,dc);

      if(mdiv==2)
      {
        (dcv+i)->angle[0]=(dc+0)->angle[0];
        (dcv+i)->angle[1]=(dc+0)->angle[1];
        *((dcv+i)->dots[0])=*((dc+0)->dots[0]);
        *((dcv+i)->dots[1])=*((dc+0)->dots[1]);
        (dcv+ndiv)->angle[0]=(dc+1)->angle[0];
        (dcv+ndiv)->angle[1]=(dc+1)->angle[1];
        *((dcv+ndiv)->dots[0])=*((dc+1)->dots[0]);
        *((dcv+ndiv)->dots[1])=*((dc+1)->dots[1]);

        ndiv++;
        break;
      }
    }
  }

  for(i=0;i<3;i++)
  {
    free((dc+i)->center);
    free((dc+i)->dots[0]);
    free((dc+i)->dots[1]);
  }
  free(dc);

  return ndiv;
}/*dividearcbyline*/

int dividearcbylineext(struct curve *cv,
                       struct line l,double *fa,double *fb,
                       struct curve *dcv)
/*DIVIDE ARC BY EXTENDED LINE.*/
/*RETURN:NUMBER DIVIDED.Fa,Fb PARAMETERS ON DIVIDING DOTS.*/
/*ARRAY DCV MALLOCED AS 3 ELEMENTS.*/
{
  int i,ndiv,mdiv;
  double r,cx,cy,x1,y1,x2,y2,xa,ya,xb,yb;
  double aa,ab;
  double Do;
  double eps=1.0E-08;
  struct curve *dc;

  r=cv->radius[0];

  cx=cv->center->d[GX];
  cy=cv->center->d[GY];

  x1=l.ends[0].d[GX];
  y1=l.ends[0].d[GY];
  x2=l.ends[1].d[GX];
  y2=l.ends[1].d[GY];

  (dcv+0)->type=cv->type;
  (dcv+1)->type=cv->type;
  (dcv+2)->type=cv->type;
  (dcv+0)->hugo=cv->hugo;
  (dcv+1)->hugo=cv->hugo;
  (dcv+2)->hugo=cv->hugo;
  *((dcv+0)->center)=*(cv->center);
  *((dcv+1)->center)=*(cv->center);
  *((dcv+2)->center)=*(cv->center);
  (dcv+0)->radius[0]=cv->radius[0];
  (dcv+1)->radius[0]=cv->radius[0];
  (dcv+2)->radius[0]=cv->radius[0];

  (dcv+0)->angle[0]=cv->angle[0];
  (dcv+0)->angle[1]=cv->angle[1];
  *((dcv+0)->dots[0])=*(cv->dots[0]);
  *((dcv+0)->dots[1])=*(cv->dots[1]);

  ndiv=1;

  /*DOT1 AS ORIGIN.*/
  cx-=x1;
  cy-=y1;
  x2-=x1;
  y2-=y1;
  if(x2==0.0 && y2==0.0) return 1; /*LINE LENGTH=0*/

  Do=r*r*(x2*x2+y2*y2)-(cx*y2-x2*cy)*(cx*y2-x2*cy);
  if(Do<eps) return 1;

  (*fa)=((cx*x2+cy*y2)-sqrt(Do))/(x2*x2+y2*y2);
  (*fb)=((cx*x2+cy*y2)+sqrt(Do))/(x2*x2+y2*y2);

  dc=(struct curve *)malloc(3*sizeof(struct curve));
  for(i=0;i<3;i++)
  {
    (dc+i)->center=(struct onode *)malloc(sizeof(struct onode));
    (dc+i)->dots[0]=(struct onode *)malloc(sizeof(struct onode));
    (dc+i)->dots[1]=(struct onode *)malloc(sizeof(struct onode));
  }

  /*DIVISION.RELEASE DOT1=ORIGIN.*/
  xa=(*fa)*x2+x1;
  ya=(*fa)*y2+y1;
  aa=definecircleangle(cv,1,xa,ya);

  mdiv=dividearc(cv,NULL,NULL,&aa,dc);

  if(mdiv==2)
  {
    for(i=0;i<mdiv;i++)
    {
      (dcv+i)->angle[0]=(dc+i)->angle[0];
      (dcv+i)->angle[1]=(dc+i)->angle[1];
      *((dcv+i)->dots[0])=*((dc+i)->dots[0]);
      *((dcv+i)->dots[1])=*((dc+i)->dots[1]);
    }
    ndiv++;
  }

  xb=(*fb)*x2+x1;
  yb=(*fb)*y2+y1;

  for(i=0;i<ndiv;i++)
  {
    ab=definecircleangle((dcv+i),1,xb,yb);

    mdiv=dividearc((dcv+i),NULL,NULL,&ab,dc);

    if(mdiv==2)
    {
      (dcv+i)->angle[0]=(dc+0)->angle[0];
      (dcv+i)->angle[1]=(dc+0)->angle[1];
      *((dcv+i)->dots[0])=*((dc+0)->dots[0]);
      *((dcv+i)->dots[1])=*((dc+0)->dots[1]);
      (dcv+ndiv)->angle[0]=(dc+1)->angle[0];
      (dcv+ndiv)->angle[1]=(dc+1)->angle[1];
      *((dcv+ndiv)->dots[0])=*((dc+1)->dots[0]);
      *((dcv+ndiv)->dots[1])=*((dc+1)->dots[1]);

      ndiv++;
      break;
    }
  }

  for(i=0;i<3;i++)
  {
    free((dc+i)->center);
    free((dc+i)->dots[0]);
    free((dc+i)->dots[1]);
  }
  free(dc);

  return ndiv;
}/*dividearcbylineext*/

double definecircleangle(struct curve *cv,int iend,
                         double x,double y)
/*DEFINE ANGLE OF CIRCLE BY DOT.*/
/*RETURN:ANGLE.*/
{
  double x1,y1,x2,y2,a1,a2;

  x1=cv->center->d[GX];
  y1=cv->center->d[GY];
  x2=x;
  y2=y;

  if(iend==0 && cv->hugo==1)
  {
    if((x2-x1)!=0.0)
    {
      a1=atan((y2-y1)/(x2-x1));

      if(x2<x1) a1+=PI;
      if(a1>=0.0) a1-=2.0*PI;
    }
    else if(y1>y2) a1=-0.5*PI;
    else if(y1<y2) a1=-1.5*PI;
    else           a1=0.0;

    return a1;
  }
  else if(iend==0 && cv->hugo==-1)
  {
    if((x2-x1)!=0.0)
    {
      a1=atan((y2-y1)/(x2-x1));

      if(x2<x1) a1+=PI;
      if(a1<0.0) a1+=2.0*PI;
    }
    else if(y1>y2) a1=1.5*PI;
    else if(y1<y2) a1=0.5*PI;
    else           a1=0.0;

    return a1;
  }

  if(iend==1 && cv->hugo==1)
  {
    a1=cv->angle[0];

    if((x2-x1)!=0.0)
    {
      a2=atan((y2-y1)/(x2-x1));

      if(x2<x1) a2+=PI;
      if(a2>=0.0 && (a2-2.0*PI)>a1) a2-=2.0*PI;
      else if(a2<=a1)               a2+=2.0*PI;
    }
    else if(y1>y2)
    {
      if(a1<-0.5*PI) a2=-0.5*PI;
      else           a2= 1.5*PI;
    }
    else if(y1<y2)
    {
      if(a1<-1.5*PI) a2=-1.5*PI;
      else           a2= 0.5*PI;
    }
    else a2=0.0;

    return a2;
  }
  else if(iend==1 && cv->hugo==-1)
  {
    a1=cv->angle[0];

    if((x2-x1)!=0.0)
    {
      a2=atan((y2-y1)/(x2-x1));

      if(x2<x1) a2+=PI;
      if(a2<=0.0 && (a2+2.0*PI)<a1) a2+=2.0*PI;
      else if(a2>=a1)               a2-=2.0*PI;
    }
    else if(y1>y2)
    {
      if(a1>1.5*PI) a2= 1.5*PI;
      else          a2=-0.5*PI;
    }
    else if(y1<y2)
    {
      if(a1>0.5*PI) a2= 0.5*PI;
      else          a2=-1.5*PI;
    }
    else a2=0.0;

    return a2;
  }

  return 0.0;
}/*definecircleangle*/

void chordrange(struct curve *cv,
                double *minx,double *miny,
                double *maxx,double *maxy)
{
  double pi20=2.0*PI, pi15=1.5*PI, pi10=PI, pi05=0.5*PI;

  *minx=cv->dots[0]->d[GX];
  *miny=cv->dots[0]->d[GY];
  *maxx=cv->dots[0]->d[GX];
  *maxy=cv->dots[0]->d[GY];

  if((*minx)>(cv->dots[1]->d[GX])) *minx=cv->dots[1]->d[GX];
  if((*miny)>(cv->dots[1]->d[GY])) *miny=cv->dots[1]->d[GY];
  if((*maxx)<(cv->dots[1]->d[GX])) *maxx=cv->dots[1]->d[GX];
  if((*maxy)<(cv->dots[1]->d[GY])) *maxy=cv->dots[1]->d[GY];

  if(cv->angle[0]==-pi20 || cv->angle[1]==-pi20 ||
     cv->angle[0]== pi20 || cv->angle[1]== pi20)
  *maxx=(cv->center->d[GX])+(cv->radius[0]);

  if((cv->hugo== 1 && cv->angle[0]<-pi15 && cv->angle[1]>-pi15)
  || (cv->hugo==-1 && cv->angle[1]<-pi15 && cv->angle[0]>-pi15))
  *maxy=(cv->center->d[GY])+(cv->radius[0]);

  if((cv->hugo== 1 && cv->angle[0]<-pi10 && cv->angle[1]>-pi10)
  || (cv->hugo==-1 && cv->angle[1]<-pi10 && cv->angle[0]>-pi10))
  *minx=(cv->center->d[GX])-(cv->radius[0]);

  if((cv->hugo== 1 && cv->angle[0]<-pi05 && cv->angle[1]>-pi05)
  || (cv->hugo==-1 && cv->angle[1]<-pi05 && cv->angle[0]>-pi05))
  *miny=(cv->center->d[GY])-(cv->radius[0]);

  if((cv->hugo== 1 && cv->angle[0]<0.0 && cv->angle[1]>0.0)
  || (cv->hugo==-1 && cv->angle[1]<0.0 && cv->angle[0]>0.0))
  *maxx=(cv->center->d[GX])+(cv->radius[0]);

  if((cv->hugo== 1 && cv->angle[0]<pi05 && cv->angle[1]>pi05)
  || (cv->hugo==-1 && cv->angle[1]<pi05 && cv->angle[0]>pi05))
  *maxy=(cv->center->d[GY])+(cv->radius[0]);

  if((cv->hugo== 1 && cv->angle[0]<pi10 && cv->angle[1]>pi10)
  || (cv->hugo==-1 && cv->angle[1]<pi10 && cv->angle[0]>pi10))
  *minx=(cv->center->d[GX])-(cv->radius[0]);

  if((cv->hugo== 1 && cv->angle[0]<pi15 && cv->angle[1]>pi15)
  || (cv->hugo==-1 && cv->angle[1]<pi15 && cv->angle[0]>pi15))
  *miny=(cv->center->d[GY])-(cv->radius[0]);

  return;
}/*chordrange*/

int insidelinecurveint(int x,int y,
                       int x1,int y1,int x2,int y2,
                       int ox,int oy)
/*DOT INSIDE DAIKEI BY LINE AND AXIS X ON ORIGIN.*/
/*ON LINE IS INSIDE.*/
{
  int ymax;

  if((x<x1 && x<x2) || (x1<x && x2<x)) return 0;

  if(y<oy) return 0;

  if(x2-x1==0) return 0;
  else ymax=y1+(y2-y1)/(x2-x1)*(x-x1);

  if(y<=ymax) return 1;
  else return 0;
}/*insidelinecurveint*/

int insidelinecurveext(double x,double y,
                       struct curve *cv,
                       struct onode origin)
/*DOT INSIDE DAIKEI BY LINE AND AXIS X ON ORIGIN.*/
/*ON LINE IS INSIDE.*/
{
  double x1,y1,x2,y2,ymax,ymin;
  double eps=1.0E-08;

  x1=cv->dots[0]->d[GX];
  y1=cv->dots[0]->d[GY];
  x2=cv->dots[1]->d[GX];
  y2=cv->dots[1]->d[GY];

  if((x<x1 && x<x2) || (x1<x && x2<x)) return 0;

  ymin=origin.d[GY];
  if(y<ymin) return 0;

  if(fabs(x2-x1)<eps)
  {
    /*if(fabs(x-x1)<eps && ymin<=y && (y<=y1 || y<=y2)) return 1;
    else*/ return 0;
  }
  else ymax=y1+(y2-y1)/(x2-x1)*(x-x1);

  if(y<=ymax) return 1;
  else return 0;
}/*insidelinecurveext*/

int dotinchord(struct curve *cv,double x,double y)
/*DOT INSIDE OR OUTSIDE CHORD.*/
{
  double r,x1,y1,x2,y2,cx,cy,f;
  double eps=5.0E-02;

  r=cv->radius[0];

  cx=cv->center->d[GX];
  cy=cv->center->d[GY];

  /*if(r<=fabs(x-cx) || r<=fabs(y-cy)) return 0;*/

  x1=cv->dots[0]->d[GX];
  y1=cv->dots[0]->d[GY];
  x2=cv->dots[1]->d[GX];
  y2=cv->dots[1]->d[GY];

  f=(x-cx)*(x-cx)+(y-cy)*(y-cy);
  if(f>r*r) return 0;

  if(fabs(x2-x1)<eps && fabs(y2-y1)<eps) return 1;

  f=(y2-y1)*(x-x1)-(x2-x1)*(y-y1);

  if((cv->hugo== 1 && f>0.0) ||
     (cv->hugo==-1 && f<0.0)) return 1;
  else return 0;

}/*dotinchord*/

int dotinpolycurve(struct polycurve *pc,double x,double y)
/*RETURN:1=INSIDE,ELSE=OUTSIDE.*/
{
  int i,ncurve,count;
  double x1,y1,x2,y2;
  struct onode nmin;

  /*FIND MIN,MAX DOT.*/
  ncurve=pc->ncurve;
  if(ncurve>0)
  {
    nmin.d[GX]=(pc->curves+0)->dots[0]->d[GX];
    nmin.d[GY]=(pc->curves+0)->dots[0]->d[GY];
    nmin.d[GZ]=0.0;
    for(i=0;i<ncurve;i++)
    {
      if((pc->curves+i)->type==CTYPE_LINE)
      {
        x1=(pc->curves+i)->dots[0]->d[GX];
        y1=(pc->curves+i)->dots[0]->d[GY];
        x2=(pc->curves+i)->dots[1]->d[GX];
        y2=(pc->curves+i)->dots[1]->d[GY];
      }
      if((pc->curves+i)->type==CTYPE_CIRCLE)
      {
        chordrange((pc->curves+i),&x1,&y1,&x2,&y2);
      }

      if(nmin.d[GX]>x1) nmin.d[GX]=x1;
      if(nmin.d[GY]>y1) nmin.d[GY]=y1;
      if(nmin.d[GX]>x2) nmin.d[GX]=x2;
      if(nmin.d[GY]>y2) nmin.d[GY]=y2;
    }
  }

  count=0;
  for(i=0;i<ncurve;i++)
  {
    /*IF INSIDE LINE COUNT++,COUNT--*/
    if(insidelinecurveext(x,y,(pc->curves+i),
                          nmin)==1)
    {
      if(((pc->curves+i)->dots[0]->d[GX]) >=
         ((pc->curves+i)->dots[1]->d[GX]))
      {
        count++;
      }
      else count--;
    }

    if((pc->curves+i)->type==CTYPE_CIRCLE)
    {
      /*IF INSIDE CHORD COUNT++,COUNT--*/
      if(dotinchord((pc->curves+i),x,y)==1)
      {
        if((pc->curves+i)->hugo== 1) count++;
        if((pc->curves+i)->hugo==-1) count--;
      }
    }
  }

  return count;
}/*dotinpolycurve*/

struct polycurve *pickpolycurve(struct polypolycurve *ppc,
                                struct viewparam vp,
                                int mx,int my)
{
  int i;
  double dmx,dmy;
  struct plane pl;
  struct onode nadd;

  pl.nods[0].d[GX]=0.0;
  pl.nods[0].d[GY]=0.0;
  pl.nods[0].d[GZ]=0.0;

  pl.nvec.dc[GX]=0.0;
  pl.nvec.dc[GY]=0.0;
  pl.nvec.dc[GZ]=1.0;

  pl.a=0.0; pl.b=0.0; pl.c=1.0; pl.d=0.0; /*PLANE:Z=0*/

  dmx=(double)( mx-(vp.Xo));
  dmy=(double)(-my+(vp.Yo));
  createnodeonplane(vp,dmx,dmy,pl,&nadd);

  if(nadd.code==0)
  {
    /*MessageBox(NULL,"Failed.","Node",MB_OK);*/
    return NULL;
  }

  for(i=0;i<(ppc->npcurve);i++)
  {
    if(dotinpolycurve((ppc->pcurves+i),nadd.d[GX],nadd.d[GY])
       ==1) return (ppc->pcurves+i);
  }

  return NULL;
}/*pickpolycurve*/

int arcinsidechord(struct curve *cvouter,struct curve *cvinner)
/*ARC INSIDE OR OUTSIDE CHORD.*/
/*RETURN:0=OUTSIDE 1=CROSSED 2=INSIDE*/
{
  int ninside=0;

  if(dotinchord(cvouter,
                (cvinner->dots[0]->d[GX]),
                (cvinner->dots[0]->d[GY]))==1) ninside++;

  if(dotinchord(cvouter,
                (cvinner->dots[1]->d[GX]),
                (cvinner->dots[1]->d[GY]))==1) ninside++;

  return ninside;
}/*arcinsidechord*/

struct line tangentoncircle(struct curve *cv,
                            struct onode *on,double *angle)
/*TANGENT OF CIRCLE ON DOT OR ANGLE.*/
{
  double r,cx,cy,ox,oy;
  struct line lt;

  r=cv->radius[0];

  cx=cv->center->d[GX];
  cy=cv->center->d[GY];

  if(on!=NULL)
  {
    ox=on->d[GX];
    oy=on->d[GY];
  }
  else if(angle!=NULL)
  {
    ox=cx+r*cos(*angle);
    oy=cy+r*sin(*angle);
  }

  lt.ends[0].d[GX]=ox;
  lt.ends[0].d[GY]=oy;
  if(cv->hugo==1)
  {
    lt.ends[1].d[GX]=ox-oy+cy;
    lt.ends[1].d[GY]=oy+ox-cx;
  }
  else if(cv->hugo==-1)
  {
    lt.ends[1].d[GX]=ox+oy-cy;
    lt.ends[1].d[GY]=oy-ox+cx;
  }

  return lt;
}/*tangentoncircle*/

void copycurve(struct curve *cvto,struct curve *cvfrom)
{
  cvto->type=cvfrom->type;
  cvto->loff=cvfrom->loff;
  cvto->hugo=cvfrom->hugo;
  *(cvto->dots[0])=*(cvfrom->dots[0]);
  *(cvto->dots[1])=*(cvfrom->dots[1]);

  if(cvfrom->type==CTYPE_CIRCLE)
  {
    *(cvto->center)=*(cvfrom->center);
    cvto->radius[0]=cvfrom->radius[0];
    cvto->angle[0]=cvfrom->angle[0];
    cvto->angle[1]=cvfrom->angle[1];
  }
  return;
}/*copycurve*/

void copypolycurve(struct polycurve *pcto,
                   struct polycurve *pcfrom)
{
  int i;

  if((pcto->ncurve)>0) freepolycurve(pcto);

  pcto->loff  =pcfrom->loff;
  pcto->ncurve=pcfrom->ncurve;
  pcto->type  =pcfrom->type;
  pcto->prop  =pcfrom->prop;

  pcto->curves=(struct curve *)malloc(pcto->ncurve
                                      *sizeof(struct curve));

  for(i=0;i<pcto->ncurve;i++)
  {
    initializecurve(pcto->curves+i);

    (pcto->curves+i)->dots[0]=(struct onode *)
                              malloc(sizeof(struct onode));
    (pcto->curves+i)->dots[1]=(struct onode *)
                              malloc(sizeof(struct onode));
    if((pcfrom->curves+i)->type==CTYPE_CIRCLE)
    {
      (pcto->curves+i)->center=(struct onode *)
                               malloc(sizeof(struct onode));
    }

    copycurve((pcto->curves+i),(pcfrom->curves+i));
  }

  return;
}/*copypolycurve*/

void copypolypolycurve(struct polypolycurve *ppcto,
                       struct polypolycurve *ppcfrom)
{
  int i;

  if((ppcto->npcurve)>0) freepolypolycurve(ppcto);

  ppcto->npcurve=ppcfrom->npcurve;
  ppcto->ici    =ppcfrom->ici;
  ppcto->hico   =ppcfrom->hico;
  ppcto->hcur   =ppcfrom->hcur;

  strcpy(ppcto->name,ppcfrom->name);

  ppcto->pcurves=(struct polycurve *)
                 malloc((ppcto->npcurve)*sizeof(struct polycurve));

  for(i=0;i<(ppcto->npcurve);i++)
  {
    (ppcto->pcurves+i)->ncurve=0;
    copypolycurve((ppcto->pcurves+i),(ppcfrom->pcurves+i));
  }

  return;
}/*copypolypolycurve*/

void drawlinecurveext(HDC hdc,
                      struct viewparam vp,
                      struct curve *cv,
                      HPEN hpwaku,HPEN hpfill,HBRUSH hbfill,
                      struct onode *origin)
/*DRAW LINE OF CURVE.*/
{
  HDC hdcC;
  HBITMAP hbitC,pbit;
  HPEN ppen;
  HBRUSH hbrush,pbrush;
  POINT *points;
  /*char s[80],str[400];*/
  int i,npoint;
  int maxX,minX,maxY,minY;
  int X1,Y1; /*POINT ON DRAWINGS.*/
  double eps=1.0E-08;
  struct onode pn[4];

  if(origin!=NULL)
  {
    /*FOR X,IF NOT VERTICAL.*/
    if(eps<fabs((cv->dots[0]->d[GX])-(cv->dots[1]->d[GX])))
    {
      /*LINES INTO POLYGON.*/
      npoint=4;
      pn[0]=*(cv->dots[0]);
      pn[1]=*(cv->dots[1]);
      pn[2].d[GX]=pn[1].d[GX];
      pn[2].d[GY]=origin->d[GY];
      pn[2].d[GZ]=0.0;
      pn[3].d[GX]=pn[0].d[GX];
      pn[3].d[GY]=origin->d[GY];
      pn[3].d[GZ]=0.0;

      /*PROJECTION.*/
      points=(POINT *)malloc(npoint*sizeof(POINT));
      if(points==NULL) return;
      for(i=0;i<npoint;i++)
      {
        if(!nodeontoscreen(pn[i],&X1,&Y1,vp)) return;
        (points+i)->x=X1;
        (points+i)->y=Y1;

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
	  pbit=(HBITMAP)SelectObject(hdcC,hbitC);

      hbrush = (HBRUSH)GetCurrentObject(hdc,OBJ_BRUSH);
	  pbrush = (HBRUSH)SelectObject(hdcC,hbrush);

      PatBlt(hdcC,minX,minY,(maxX-minX),(maxY-minY),PATCOPY);

      /*PEN,BRUSH.*/
      if(((cv->dots[0])->d[GX])>=((cv->dots[1])->d[GX]))
      {
        ppen=(HPEN)SelectObject(hdcC,hpfill);
        SelectObject(hdcC,hbfill);
      }
      else if(((cv->dots[0])->d[GX])<((cv->dots[1])->d[GX]))
      {
        ppen=(HPEN)SelectObject(hdcC,hpfill);
        SelectObject(hdcC,hbfill);
      }

      Polygon(hdcC,points,npoint);

      BitBlt(hdc,minX,minY,(maxX-minX),(maxY-minY),
             hdcC,minX,minY,
             SRCINVERT);
      /*RASTER OPERATIONS:SRCPAINT    NG*/
      /*                  SRCCOPY     NG*/
      /*                  SRCAND      NG*/
      /*                  SRCERASE    NG*/
      /*                  SRCINVERT   OK*/
      /*                  NOTSRCCOPY  NG*/
      /*                  NOTSRCERASE NG*/

      SelectObject(hdcC,ppen);
      SelectObject(hdcC,pbrush);
      SelectObject(hdcC,pbit);
      DeleteObject(hbrush);
      DeleteObject(hbitC);
      DeleteDC(hdcC);
      free(points);
    }

    /*FOR Y*/
    /*Polygon(hdc,points,npoint);*/
  }
  else if(cv->type==CTYPE_LINE)
  {
    ppen=(HPEN)SelectObject(hdc,hpwaku);
    drawgloballine(hdc,vp,*(cv->dots[0]),*(cv->dots[1]));
    SelectObject(hdc,ppen);
  }

  return;
}/*drawlinecurveext*/

void drawglobalchord(HDC hdc,
                     struct viewparam vp,
                     struct curve cv,
                     HPEN hpwaku,HPEN hpfill,HPEN hphollow,
                     HBRUSH hbfill,HBRUSH hbhollow)
/*DRAW CHORD OF CURVE.*/
{
  int minX,minY,maxX,maxY,X1,Y1,X2,Y2;
  struct onode node;

  if(cv.type!=CTYPE_CIRCLE) return;

  node.d[GX]=(cv.center)->d[GX]-cv.radius[0];
  node.d[GY]=(cv.center)->d[GY]+cv.radius[0];
  node.d[GZ]=0.0;
  if(!nodeontoscreen(node,&minX,&minY,vp)) return; /*MIN DOT*/
  node.d[GX]=(cv.center)->d[GX]+cv.radius[0];
  node.d[GY]=(cv.center)->d[GY]-cv.radius[0];
  node.d[GZ]=0.0;
  if(!nodeontoscreen(node,&maxX,&maxY,vp)) return; /*MAX DOT*/

  if(cv.dots[0]!=NULL) /*START DOT.*/
  {
    node=*(cv.dots[0]);
    if(!nodeontoscreen(node,&X1,&Y1,vp)) return;
  }
  else return;

  if(cv.dots[1]!=NULL) /*END DOT.*/
  {
    node=*(cv.dots[1]);
    if(!nodeontoscreen(node,&X2,&Y2,vp)) return;
  }
  else return;

    if(cv.hugo==1)
    {
      SelectObject(hdc,hpfill);
      SelectObject(hdc,hbfill);
      Chord(hdc,minX,minY,maxX,maxY,X1,Y1,X2,Y2);

      SelectObject(hdc,hpwaku);
      Arc(hdc,minX,minY,maxX,maxY,X1,Y1,X2,Y2);
    }
    if(cv.hugo==-1)
    {
      SelectObject(hdc,hphollow);
      SelectObject(hdc,hbhollow);
      Chord(hdc,minX,minY,maxX,maxY,X2,Y2,X1,Y1);

      SelectObject(hdc,hpwaku);
      Arc(hdc,minX,minY,maxX,maxY,X2,Y2,X1,Y1);
    }

  return;
}/*drawglobalchord*/

void drawglobalchordext(HDC hdc,
                        struct viewparam vp,
                        struct curve cv,
                        HPEN hpwaku,
                        COLORREF crwaku,
                        COLORREF crfill,
                        COLORREF crhollow,
                        int fillorhollow,
                        struct onode origin)
/*DRAW CHORD OF CURVE.*/
{
  HDC hdcC;
  HBITMAP hbitC,pbit;
  HPEN ppen;
  HBRUSH hbrush,pbrush;
  int Ox,Oy,X0,Y0,X1,Y1,X2,Y2;
  int maxX,minX,maxY,minY;
  int i,j,count;
  double r,oy,cx;
  double delta,minx,maxx,x,y;
  struct onode node;

  if(cv.type!=CTYPE_CIRCLE) return;

  node.d[GX]=(cv.center)->d[GX]-cv.radius[0];
  node.d[GY]=(cv.center)->d[GY]+cv.radius[0];
  node.d[GZ]=0.0;
  if(!nodeontoscreen(node,&minX,&minY,vp)) return; /*MIN DOT*/
  node.d[GX]=(cv.center)->d[GX]+cv.radius[0];
  node.d[GY]=(cv.center)->d[GY]-cv.radius[0];
  node.d[GZ]=0.0;
  if(!nodeontoscreen(node,&maxX,&maxY,vp)) return; /*MAX DOT*/

  node=*(cv.dots[0]);
  if(!nodeontoscreen(node,&X1,&Y1,vp)) return; /*1st END*/
  node=*(cv.dots[1]);
  if(!nodeontoscreen(node,&X2,&Y2,vp)) return; /*2nd END*/

  if(cv.hugo==-1) /*REVERSE.*/
  {
    X0=X1; X1=X2; X2=X0;
    Y0=Y1; Y1=Y2; Y2=Y0;
  }

  if(fillorhollow)
  {
    if(!nodeontoscreen(origin,&Ox,&Oy,vp)) return; /*ORIGIN*/
    /*ox=origin.d[GX];*/
    oy=origin.d[GY];
    cx=(cv.center)->d[GX];
    /*cy=(cv.center)->d[GY];*/
    r =cv.radius[0];
    minx=cx-r;
    /*miny=cy-r;*/
    maxx=cx+r;
    /*maxy=cy+r;*/

    delta=(maxx-minx)/(double)(maxX-minX); /*LENGTH=1 ON SCREEN.*/

    hdcC = CreateCompatibleDC(hdc);

    hbitC = CreateCompatibleBitmap(hdc,maxX,maxY);
    pbit=(HBITMAP)SelectObject(hdcC,hbitC);

    hbrush = (HBRUSH)GetCurrentObject(hdc,OBJ_BRUSH);
    pbrush = (HBRUSH)SelectObject(hdcC,hbrush);

    PatBlt(hdcC,minX,minY,(maxX-minX),(maxY-minY),PATCOPY);

    x=minx;
    for(i=minX; i<=maxX; i++) /*SCANNING X*/
    {
      y=oy;
      for(j=Oy; j>=minY; j--) /*SCANNING Y*/
      {
        count=0;

        /*IF INSIDE LINE  COUNT++*/
        if(insidelinecurveext(x,y,&cv,origin)==1) count++;

        /*IF INSIDE CHORD COUNT++*/
        if(dotinchord(&cv,x,y)==1) count++;

        /*IF KISUU,PAINT DOT{i,j}.*/
        if(count%2==1) SetPixelV(hdcC,i,j,crfill);

        y+=delta;
      }
      x+=delta;
    }

    BitBlt(hdc,minX,minY,(maxX-minX),(maxY-minY),
           hdcC,minX,minY,
           SRCINVERT);

    SelectObject(hdcC,pbrush);
    SelectObject(hdcC,pbit);
    DeleteObject(hbrush);
    DeleteObject(hbitC);
    DeleteDC(hdcC);
  }
  else
  {
    ppen=(HPEN)SelectObject(hdc,hpwaku);
    Arc(hdc,minX,minY,maxX,maxY,X1,Y1,X2,Y2);
    SelectObject(hdc,ppen);
  }
  return;
}/*drawglobalchordext*/

void drawglobalpolycurve(HDC hdc,
                         struct viewparam vp,
                         struct polypolycurve *pc)
{
  HPEN hpen,ppen;
  COLORREF crmore,crposi,crzero,crnega,crless;
  /*char str[256];*/
  int i,j,ii,jj,ncurve,count;
  int minX,minY,maxX,maxY,r,g,b;
  double arrowlength,x,y,x1,y1,x2,y2,deltax,deltay;
  struct onode on,hn,nmax,nmin;

  /*FIND MIN,MAX DOT.*/
  count=0;
  for(jj=0;jj<(pc->npcurve);jj++)
  {
    if(((pc->pcurves+jj)->ncurve)>0)
    {
      if(count==0)
      {
        nmin.d[GX]=((pc->pcurves+jj)->curves+0)->dots[0]->d[GX];
        nmin.d[GY]=((pc->pcurves+jj)->curves+0)->dots[0]->d[GY];
        nmin.d[GZ]=0.0;
        nmax.d[GX]=((pc->pcurves+jj)->curves+0)->dots[0]->d[GX];
        nmax.d[GY]=((pc->pcurves+jj)->curves+0)->dots[0]->d[GY];
        nmax.d[GZ]=0.0;
      }
      count=1;
      ncurve=(pc->pcurves+jj)->ncurve;

      for(i=0;i<ncurve;i++)
      {
        if(((pc->pcurves+jj)->curves+i)->type==CTYPE_LINE)
        {
          x1=((pc->pcurves+jj)->curves+i)->dots[0]->d[GX];
          y1=((pc->pcurves+jj)->curves+i)->dots[0]->d[GY];
          x2=((pc->pcurves+jj)->curves+i)->dots[1]->d[GX];
          y2=((pc->pcurves+jj)->curves+i)->dots[1]->d[GY];
        }
        if(((pc->pcurves+jj)->curves+i)->type==CTYPE_CIRCLE)
        {
          chordrange(((pc->pcurves+jj)->curves+i),&x1,&y1,&x2,&y2);
        }

        if(nmin.d[GX]>x1) nmin.d[GX]=x1;
        if(nmin.d[GY]>y1) nmin.d[GY]=y1;
        if(nmin.d[GX]>x2) nmin.d[GX]=x2;
        if(nmin.d[GY]>y2) nmin.d[GY]=y2;

        if(nmax.d[GX]<x1) nmax.d[GX]=x1;
        if(nmax.d[GY]<y1) nmax.d[GY]=y1;
        if(nmax.d[GX]<x2) nmax.d[GX]=x2;
        if(nmax.d[GY]<y2) nmax.d[GY]=y2;
      }
    }
  }
  if(count)
  {
    if(!nodeontoscreen(nmin,&minX,&minY,vp)) return; /*MIN DOT*/
    if(!nodeontoscreen(nmax,&maxX,&maxY,vp)) return; /*MAX DOT*/

    deltax=(nmax.d[GX]-nmin.d[GX])/(double)(maxX-minX);
    deltay=(nmax.d[GY]-nmin.d[GY])/(double)(minY-maxY);

/*sprintf(str,"Min=(%.3f %.3f) Max=(%.3f %.3f)",
        nmin.d[GX],nmin.d[GY],nmax.d[GX],nmax.d[GY]);
MessageBox(NULL,str,"PolyCurve",MB_OK);*/
  }

  /*PAINT DOTS.*/
  for(jj=0;jj<(pc->npcurve);jj++)
  {
    ncurve=(pc->pcurves+jj)->ncurve;

    r=(pc->pcurves+jj)->prop.r;
    g=(pc->pcurves+jj)->prop.g;
    b=(pc->pcurves+jj)->prop.b;

    crmore=RGB(r,      g,      b);       /*+N*/
    crposi=RGB(r/3,    g/3,    b/3);     /*+1*/
    crzero=RGB(  0,      0,      0);     /* 0*/
    crnega=RGB( 85-r/3, 85-g/3, 85-b/3); /*-1*/
    crless=RGB(255-r,  255-g,  255-b);   /*-N*/

    x=nmin.d[GX];
    for(i=minX; i<=maxX; i++) /*SCANNING X*/
    {
      y=nmin.d[GY];
      for(j=minY; j>=maxY; j--) /*SCANNING Y*/
      {
        count=0;

        for(ii=0;ii<ncurve;ii++)
        {
          /*IF INSIDE LINE COUNT++,COUNT--*/
          if(insidelinecurveext(x,y,((pc->pcurves+jj)->curves+ii),
                                nmin)==1)
          {
            if((((pc->pcurves+jj)->curves+ii)->dots[0]->d[GX]) >=
               (((pc->pcurves+jj)->curves+ii)->dots[1]->d[GX]))
            {
              count++;
            }
            else count--;
          }

          if(((pc->pcurves+jj)->curves+ii)->type==CTYPE_CIRCLE)
          {
            /*IF INSIDE CHORD COUNT++,COUNT--*/
            if(dotinchord(((pc->pcurves+jj)->curves+ii),x,y)==1)
            {
              if(((pc->pcurves+jj)->curves+ii)->hugo== 1) count++;
              if(((pc->pcurves+jj)->curves+ii)->hugo==-1) count--;
            }
          }
        }

        /*IF KISUU,PAINT DOT{i,j}.*/
        if(count>= 2) SetPixelV(hdc,i,j,crmore);
        if(count== 1) SetPixelV(hdc,i,j,crposi);
        /*if(count== 0) SetPixelV(hdc,i,j,crzero);*/
        if(count==-1) SetPixelV(hdc,i,j,crnega);
        if(count<=-2) SetPixelV(hdc,i,j,crless);

        y+=deltay;
      }
      x+=deltax;
    }
  }

  /*DRAW FRAMES.*/
  hpen = CreatePen(PS_SOLID,1,RGB(255,255,255)); /*LINE*/
  for(jj=0;jj<(pc->npcurve);jj++)
  {
    ncurve=(pc->pcurves+jj)->ncurve;
    for(i=0;i<ncurve;i++)
    {
      /*DRAW LINES.*/
      if(((pc->pcurves+jj)->curves+i)->type==CTYPE_LINE)
      {
        drawlinecurveext(hdc,vp,((pc->pcurves+jj)->curves+i),
                         hpen,NULL,NULL,NULL);
      }
      /*DRAW ARCS.*/
      if(((pc->pcurves+jj)->curves+i)->type==CTYPE_CIRCLE)
      {
        drawglobalchordext(hdc,vp,*((pc->pcurves+jj)->curves+i),
                           hpen,NULL,NULL,NULL,0,nmin);
      }
    }
  }
  DeleteObject(hpen);

  /*DRAW AXIS.*/
  SetTextColor(hdc,RGB(150,150,150));                  /*TEXT GRAY.*/
  hpen=CreatePen(PS_SOLID,1,RGB(0,0,255));           /*ARROW COLOR.*/
  ppen=(HPEN)SelectObject(hdc,hpen);

  arrowlength=vp.dparam.gaxis;

  on=setcoord(4.0,0.0,0.0);
  drawglobaltext(hdc,vp,on,"O");                           /*ORIGIN*/

  hn=setcoord(arrowlength+4.0,10.0,0.0);
  drawglobaltext(hdc,vp,hn,"X");                           /*TEXT X*/
  on=setcoord(-arrowlength,0.0,0.0);
  hn=setcoord(arrowlength,0.0,0.0);
  drawglobalarrow(hdc,vp,on,hn,0.0);                      /*ARROW X*/

  hn=setcoord(-4.0,arrowlength+20.0,0.0);
  drawglobaltext(hdc,vp,hn,"Y");                           /*TEXT Y*/
  on=setcoord(0.0,-arrowlength,0.0);
  hn=setcoord(0.0,arrowlength,0.0);
  drawglobalarrow(hdc,vp,on,hn,0.0);                      /*ARROW Y*/

  SelectObject(hdc,ppen);
  DeleteObject(hpen);

  return;
}/*drawglobalpolycurve*/

HICON createpropertyicon(HINSTANCE hinst,HDC hdc,
                         struct oprop *op,int iorc)
{
  HDC hdcI,hdcM;
  HBITMAP pbitI,pbitM;
  HBRUSH hbrushM,pbrushM;
  HICON hicon;
  ICONINFO ici;
  COLORREF cr;
  /*char str[256];*/
  int ix,iy;
  int i,j;
  int ir,ig,ib;

  /*COPY TO ICON.*/
  ix=32;
  iy=32;

  ici.fIcon=iorc; /*TRUE=ICON FALSE=CURSOR*/
  ici.xHotspot=11;
  ici.yHotspot=13;

  hdcM = CreateCompatibleDC(hdc);
  hdcI = CreateCompatibleDC(hdc);

  ici.hbmMask =CreateCompatibleBitmap(hdc,ix,iy);
  pbitM = (HBITMAP)SelectObject(hdcM,ici.hbmMask);
  hbrushM = (HBRUSH)CreateSolidBrush(RGB(255,255,255));
  pbrushM = (HBRUSH)SelectObject(hdcM,hbrushM);
  PatBlt(hdcM,0,0,ix,iy,PATCOPY);

  ici.hbmColor = LoadBitmap(hinst,"BITMAPPROP");
  pbitI = (HBITMAP)SelectObject(hdcI,ici.hbmColor);

  for(i=1; i<=ix; i++)
  {
    for(j=1; j<=iy; j++)
    {
      cr=GetPixel(hdcI,i,j);

      if(cr==0x00ff0000)
      {
        ir=op->r;
        ig=op->g;
        ib=op->b;

        SetPixelV(hdcI,i,j,RGB(ir,ig,ib));
        /*SetPixelV(hdcM,i,j,RGB(0,0,0));*/
      }
    }
  }

  SelectObject(hdcM,pbitM);
  SelectObject(hdcM,pbrushM);
  DeleteObject(hbrushM);
  DeleteDC(hdcM);

  SelectObject(hdcI,pbitI);
  DeleteDC(hdcI);

  hicon=CreateIconIndirect(&ici);
  return hicon;
}/*createpropertyicon*/

HICON createsectionicon(HDC hdc,int ix,int iy,
                        struct viewparam vp,
                        struct polypolycurve *pc)
{
  HDC hdcC,hdcI,hdcM;
  HBITMAP hbitC,pbitC,pbitI,pbitM;
  HPEN hpen;
  HBRUSH hbrushC,pbrushC,hbrushI,pbrushI,hbrushM,pbrushM;
  HICON hicon;
  COLORREF crmore,crposi,crzero,crnega,crless,cr;
  /*char str[256];*/
  int i,j,ii,jj,ncurve,count,ndelta,nx,ny;
  int minX,minY,maxX,maxY,r,g,b,ir,ig,ib,xx,yy;
  int iblack;
  double x,y,x1,y1,x2,y2,deltax,deltay;
  struct onode nmax,nmin;

  /*FIND MIN,MAX DOT.*/
  count=0;
  for(jj=0;jj<(pc->npcurve);jj++)
  {
    if(((pc->pcurves+jj)->ncurve)>0)
    {
      if(count==0)
      {
        nmin.d[GX]=((pc->pcurves+jj)->curves+0)->dots[0]->d[GX];
        nmin.d[GY]=((pc->pcurves+jj)->curves+0)->dots[0]->d[GY];
        nmin.d[GZ]=0.0;
        nmax.d[GX]=((pc->pcurves+jj)->curves+0)->dots[0]->d[GX];
        nmax.d[GY]=((pc->pcurves+jj)->curves+0)->dots[0]->d[GY];
        nmax.d[GZ]=0.0;
      }
      count=1;
      ncurve=(pc->pcurves+jj)->ncurve;

      for(i=0;i<ncurve;i++)
      {
        if(((pc->pcurves+jj)->curves+i)->type==CTYPE_LINE)
        {
          x1=((pc->pcurves+jj)->curves+i)->dots[0]->d[GX];
          y1=((pc->pcurves+jj)->curves+i)->dots[0]->d[GY];
          x2=((pc->pcurves+jj)->curves+i)->dots[1]->d[GX];
          y2=((pc->pcurves+jj)->curves+i)->dots[1]->d[GY];
        }
        if(((pc->pcurves+jj)->curves+i)->type==CTYPE_CIRCLE)
        {
          chordrange(((pc->pcurves+jj)->curves+i),&x1,&y1,&x2,&y2);
        }

        if(nmin.d[GX]>x1) nmin.d[GX]=x1;
        if(nmin.d[GY]>y1) nmin.d[GY]=y1;
        if(nmin.d[GX]>x2) nmin.d[GX]=x2;
        if(nmin.d[GY]>y2) nmin.d[GY]=y2;

        if(nmax.d[GX]<x1) nmax.d[GX]=x1;
        if(nmax.d[GY]<y1) nmax.d[GY]=y1;
        if(nmax.d[GX]<x2) nmax.d[GX]=x2;
        if(nmax.d[GY]<y2) nmax.d[GY]=y2;
      }
    }
  }
  if(count)
  {
    if(!nodeontoscreen(nmin,&minX,&minY,vp)) return NULL; /*MIN*/
    if(!nodeontoscreen(nmax,&maxX,&maxY,vp)) return NULL; /*MAX*/

    deltax=(nmax.d[GX]-nmin.d[GX])/(double)(maxX-minX);
    deltay=(nmax.d[GY]-nmin.d[GY])/(double)(minY-maxY);

/*sprintf(str,"Min=(%.3f %.3f) Max=(%.3f %.3f)",
        nmin.d[GX],nmin.d[GY],nmax.d[GX],nmax.d[GY]);
MessageBox(NULL,str,"Icon",MB_OK);*/
  }

  hdcC = CreateCompatibleDC(hdc);

  hbitC = CreateCompatibleBitmap(hdc,maxX,minY);
  pbitC = (HBITMAP)SelectObject(hdcC,hbitC);

  /*hbrushC = (HBRUSH)GetCurrentObject(hdc,OBJ_BRUSH);*/
  /*hbrushC = (HBRUSH)CreateSolidBrush(RGB(0,0,0));*/
  hbrushC = (HBRUSH)CreateSolidBrush(RGB(255,255,255));
  pbrushC = (HBRUSH)SelectObject(hdcC,hbrushC);

  PatBlt(hdcC,minX,maxY,(maxX-minX),(minY-maxY),PATCOPY);

  /*PAINT DOTS.*/
  for(jj=0;jj<(pc->npcurve);jj++)
  {
    ncurve=(pc->pcurves+jj)->ncurve;

    r=(pc->pcurves+jj)->prop.r;
    g=(pc->pcurves+jj)->prop.g;
    b=(pc->pcurves+jj)->prop.b;

    crless=RGB(    r,    g,    b); /*-N*/
    crnega=RGB(    r,    g,    b); /*-1*/
    crzero=RGB(    0,    0,    0); /* 0*/
    crposi=RGB(255-r,255-g,255-b); /*+1*/
    crmore=RGB(255-r,255-g,255-b); /*+N*/

    x=nmin.d[GX];
    for(i=minX; i<=maxX; i++) /*SCANNING X*/
    {
      y=nmin.d[GY];
      for(j=minY; j>=maxY; j--) /*SCANNING Y*/
      {
        count=0;

        for(ii=0;ii<ncurve;ii++)
        {
          /*IF INSIDE LINE COUNT++,COUNT--*/
          if(insidelinecurveext(x,y,((pc->pcurves+jj)->curves+ii),
                                nmin)==1)
          {
            if((((pc->pcurves+jj)->curves+ii)->dots[0]->d[GX]) >=
               (((pc->pcurves+jj)->curves+ii)->dots[1]->d[GX]))
            {
              count++;
            }
            else count--;
          }

          if(((pc->pcurves+jj)->curves+ii)->type==CTYPE_CIRCLE)
          {
            /*IF INSIDE CHORD COUNT++,COUNT--*/
            if(dotinchord(((pc->pcurves+jj)->curves+ii),x,y)==1)
            {
              if(((pc->pcurves+jj)->curves+ii)->hugo== 1) count++;
              if(((pc->pcurves+jj)->curves+ii)->hugo==-1) count--;
            }
          }
        }

        /*IF KISUU,PAINT DOT{i,j}.*/
        if(count>= 2) SetPixelV(hdcC,i,j,crmore);
        if(count== 1) SetPixelV(hdcC,i,j,crposi);
        /*if(count== 0) SetPixelV(hdc,i,j,crzero);*/
        if(count==-1) SetPixelV(hdcC,i,j,crnega);
        if(count<=-2) SetPixelV(hdcC,i,j,crless);

        y+=deltay;
      }
      x+=deltax;
    }
  }

  /*DRAW FRAMES.*/
  /*hpen = CreatePen(PS_SOLID,1,RGB(255,255,255));*/
  hpen = CreatePen(PS_SOLID,1,RGB(0,0,0));
  for(jj=0;jj<(pc->npcurve);jj++)
  {
    ncurve=(pc->pcurves+jj)->ncurve;
    for(i=0;i<ncurve;i++)
    {
      /*DRAW LINES.*/
      if(((pc->pcurves+jj)->curves+i)->type==CTYPE_LINE)
      {
        drawlinecurveext(hdcC,vp,((pc->pcurves+jj)->curves+i),
                         hpen,NULL,NULL,NULL);
      }
      /*DRAW ARCS.*/
      if(((pc->pcurves+jj)->curves+i)->type==CTYPE_CIRCLE)
      {
        drawglobalchordext(hdcC,vp,*((pc->pcurves+jj)->curves+i),
                           hpen,NULL,NULL,NULL,0,nmin);
      }
    }
  }
  DeleteObject(hpen);

/*MessageBox(NULL,"Painted.","Icon",MB_OK);*/

  /*COPY TO ICON.*/
  pc->ici.fIcon=TRUE; /*TRUE=ICON FALSE=CURSOR*/
  pc->ici.xHotspot=ix/2;
  pc->ici.yHotspot=iy/2;

  hdcM = CreateCompatibleDC(hdc);
  hdcI = CreateCompatibleDC(hdc);

  pc->ici.hbmMask =CreateCompatibleBitmap(hdc,ix,iy);
  pc->ici.hbmColor=CreateCompatibleBitmap(hdc,ix,iy);

  pbitM = (HBITMAP)SelectObject(hdcM,(pc->ici.hbmMask));
  hbrushM = (HBRUSH)CreateSolidBrush(RGB(255,255,255));
  pbrushM = (HBRUSH)SelectObject(hdcM,hbrushM);
  PatBlt(hdcM,0,0,ix,iy,PATCOPY);

  pbitI = (HBITMAP)SelectObject(hdcI,(pc->ici.hbmColor));
  hbrushI = (HBRUSH)CreateSolidBrush(RGB(0,0,0));
  pbrushI = (HBRUSH)SelectObject(hdcI,hbrushI);
  PatBlt(hdcI,0,0,ix,iy,PATCOPY);

  nx=(maxX-minX)/ix+1;
  ny=(minY-maxY)/iy+1;

  if(nx>ny) ndelta=nx;
  else      ndelta=ny;

/*sprintf(str,"Icon Max=%dx%d",GetSystemMetrics(SM_CXICON),
                             GetSystemMetrics(SM_CYICON));
MessageBox(NULL,str,"Icon",MB_OK);*/

  for(i=1; i<=ix; i++)
  {
    for(j=iy; j>=1; j--)
    {
      ir=0; ig=0; ib=0;
      count=0;
      iblack=0;

      for(ii=0; ii<ndelta; ii++)
      {
        xx=minX+(i-1)*ndelta+ii;
        if(xx>=maxX) break;

        for(jj=ndelta-1; jj>=0; jj--)
        {
          yy=minY-(iy-j)*ndelta-(ndelta-jj)+1;
          if(yy<maxY) break;

          cr=GetPixel(hdcC,xx,yy);

          if(cr==0x00000000)
          {
            iblack=1;
            break;
          }

          if(cr!=0x00ffffff)
          {
            count++;

            ir+=GetRValue(cr);
            ig+=GetGValue(cr);
            ib+=GetBValue(cr);
          }
        }
        if(iblack) break;
      }

      if(iblack) SetPixelV(hdcI,i,j,RGB(255,255,255));
      else if(count>0)
      {
        r=ir/count; g=ig/count; b=ib/count;

        SetPixelV(hdcI,i,j,RGB(r,g,b));
        /*SetPixelV(hdcM,i,j,RGB(0,0,0));*/
      }
    }
  }

  SelectObject(hdcC,pbrushC);
  SelectObject(hdcC,pbitC);
  DeleteObject(hbrushC);
  DeleteObject(hbitC);
  DeleteDC(hdcC);

  SelectObject(hdcM,pbitM);
  SelectObject(hdcM,pbrushM);
  DeleteObject(hbrushM);
  DeleteDC(hdcM);

  SelectObject(hdcI,pbitI);
  SelectObject(hdcI,pbrushI);
  DeleteObject(hbrushI);
  DeleteDC(hdcI);

  hicon=CreateIconIndirect(&(pc->ici));
  return hicon;
}/*createsectionicon*/

void drawiconbmp(HDC hdc,int x,int y,
                 ICONINFO *ici,DWORD rop)
{
  HDC hdcC;
  HBITMAP pbit;
  HBRUSH hbrush,pbrush;

  hdcC = CreateCompatibleDC(hdc);
  pbit = (HBITMAP)SelectObject(hdcC,(ici->hbmColor));

  /*hbrush = (HBRUSH)GetCurrentObject(hdc,OBJ_BRUSH);
  pbrush = (HBRUSH)SelectObject(hdcC,hbrush);*/
  hbrush = (HBRUSH)CreateSolidBrush(RGB(190,190,190));
  pbrush = (HBRUSH)SelectObject(hdc,hbrush);
  PatBlt(hdc,x,y,32,32,PATCOPY);

  BitBlt(hdc,x,y,32,32,hdcC,0,0,rop);

  SelectObject(hdc,pbrush);
  DeleteObject(hbrush);
  SelectObject(hdcC,pbit);
  DeleteDC(hdcC);
  return;
}/*drawiconext*/

void rotatepolycurve(struct polypolycurve *pc,double theta)
/*ROTATE POLYCURVE.*/
{
  int i,ii;
  struct onode dot;
  struct curve *c;

  for(ii=0;ii<(pc->npcurve);ii++)
  {
    for(i=0;i<((pc->pcurves+ii)->ncurve);i++)
    {
      c=((pc->pcurves+ii)->curves+i);

      dot=*(c->dots[0]);
      c->dots[0]->d[GX]=cos(theta)*(dot.d[GX])
                       -sin(theta)*(dot.d[GY]);
      c->dots[0]->d[GY]=sin(theta)*(dot.d[GX])
                       +cos(theta)*(dot.d[GY]);

      dot=*(c->dots[1]);
      c->dots[1]->d[GX]=cos(theta)*(dot.d[GX])
                       -sin(theta)*(dot.d[GY]);
      c->dots[1]->d[GY]=sin(theta)*(dot.d[GX])
                       +cos(theta)*(dot.d[GY]);

      if(c->type==CTYPE_CIRCLE)
      {
        dot=*(c->center);
        c->center->d[GX]=cos(theta)*(dot.d[GX])
                        -sin(theta)*(dot.d[GY]);
        c->center->d[GY]=sin(theta)*(dot.d[GX])
                        +cos(theta)*(dot.d[GY]);

        c->angle[0]+=theta;
        c->angle[1]+=theta;
        if(c->hugo==1 && (c->angle[0])>0.0)
        {
          c->angle[0]-=2.0*PI;
          c->angle[1]-=2.0*PI;
        }
        if(c->hugo==1 && (c->angle[0])<(-2.0*PI))
        {
          c->angle[0]+=2.0*PI;
          c->angle[1]+=2.0*PI;
        }
        if(c->hugo==-1 && (c->angle[0])<0.0)
        {
          c->angle[0]+=2.0*PI;
          c->angle[1]+=2.0*PI;
        }
        if(c->hugo==-1 && (c->angle[0])>(2.0*PI))
        {
          c->angle[0]-=2.0*PI;
          c->angle[1]-=2.0*PI;
        }
      }
    }
  }

  return;
}/*rotatepolycurve*/

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
  R=(int)(radius+0.5);

  Arc(hdc,(Ox-R),(Oy-R),(Ox+R),(Oy+R),0,0,0,0);     /*HOLLOW CIRCLE*/

  return;
}/*drawcircleonglobaldot*/

void fillglobalcircle(HDC hdc,struct viewparam vp,
                      struct onode center,double radius,
                      int r,int g,int b,double orate,
                      DWORD rop) /*RASTER OPERATION.*/
/*FILL CIRCLE ON SCREEN WITH SEMITRANSPARENT BRUSH.*/
{
  HDC hdcC;
  HBITMAP hbitC,pbit;
  HPEN hpen,ppen;
  HBRUSH hbrush,pbrush;
  int maxX,minX,maxY,minY;
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
  R=(int)(radius+0.5);

  maxX=Ox+R;
  minX=Ox-R;
  maxY=Oy+R;
  minY=Oy-R;

  hdcC = CreateCompatibleDC(hdc);

  hbitC = CreateCompatibleBitmap(hdc,maxX-minX,maxY-minY);
  pbit=(HBITMAP)SelectObject(hdcC,hbitC);

  hbrush = (HBRUSH)GetCurrentObject(hdc,OBJ_BRUSH);
  pbrush = (HBRUSH)SelectObject(hdcC,hbrush);

  PatBlt(hdcC,0,0,(maxX-minX),(maxY-minY),PATCOPY);

  hpen = CreatePen(PS_SOLID,1,RGB(r,g,b));
  ppen = (HPEN)SelectObject(hdcC,hpen);

  if(rop==SRCAND)
  {
    hbrush = (HBRUSH)CreateSolidBrush(RGB((int)((1.0-orate)*(255-r)),
                                          (int)((1.0-orate)*(255-g)),
                                          (int)((1.0-orate)*(255-b))));
  }
  else
  {
    hbrush = (HBRUSH)CreateSolidBrush(RGB((int)(orate*r),
                                          (int)(orate*g),
                                          (int)(orate*b)));
  }
  SelectObject(hdcC,hbrush);

  Ellipse(hdcC,0,0,2*R,2*R);

  BitBlt(hdc,minX,minY,(maxX-minX),(maxY-minY),
         hdcC,0,0,
         rop); /*GROUND BLACK:SRCPAINT WHITE:SRCAND.*/

  SelectObject(hdcC,ppen);
  SelectObject(hdcC,pbrush);
  SelectObject(hdcC,pbit);
  DeleteObject(hpen);
  DeleteObject(hbrush);
  DeleteObject(hbitC);
  DeleteDC(hdcC);

  return;
}/*fillglobalcircle*/

void fillglobalban(HDC hdc,struct viewparam vp,struct obans gb,
				   int r,int g,int b,
				   DWORD rop) /*RASTER OPERATION.*/
/*FILL BAN ON SCREEN WITH SEMITRANSPARENT BRUSH.*/
{
  HDC hdcC;
  HBITMAP hbitC,pbit;
  HBRUSH hbrush,pbrush;
  HPEN hpen,ppen;
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

  /*hbitC = CreateCompatibleBitmap(hdc,maxX,maxY);*/
  hbitC = CreateCompatibleBitmap(hdc,maxX-minX,maxY-minY);
  pbit=(HBITMAP)SelectObject(hdcC,hbitC);

  hbrush = (HBRUSH)GetCurrentObject(hdc,OBJ_BRUSH);
  pbrush = (HBRUSH)SelectObject(hdcC,hbrush);

  /*PatBlt(hdcC,minX,minY,(maxX-minX),(maxY-minY),PATCOPY);*/
  PatBlt(hdcC,0,0,(maxX-minX),(maxY-minY),PATCOPY);

  hbrush = (HBRUSH)CreateSolidBrush(RGB(r,g,b));
  SelectObject(hdcC,hbrush);

  /*hpen=CreatePen(PS_SOLID,1,RGB(0,0,0));*/ /*EDGE BLACK*/
  hpen=CreatePen(PS_SOLID,1,RGB(r,g,b));
  ppen=(HPEN)SelectObject(hdcC,hpen);

  for(i=0;i<gb.nnod;i++)
  {
	(points+i)->x-=minX;
    (points+i)->y-=minY;
  }

  Polygon(hdcC,points,gb.nnod);

  /*BitBlt(hdc,minX,minY,(maxX-minX),(maxY-minY),
         hdcC,minX,minY,
         rop);*/
  BitBlt(hdc,minX,minY,(maxX-minX),(maxY-minY),
		 hdcC,0,0,
         rop); /*GROUND BLACK:SRCPAINT WHITE:SRCAND.*/

  SelectObject(hdcC,ppen);
  SelectObject(hdcC,pbrush);
  SelectObject(hdcC,pbit);
  DeleteObject(hpen);
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
  double orate;
  struct onode gn1,gn2;

  orate=(wdraw.childs+1)->org.opaque;

  /*SRCAND  :ONPRINTER,ONPREVIEW*/
  /*SRCPAINT:ONSCREEN*/
  if(rop==SRCAND)
  {
    /*hpen=CreatePen(PS_SOLID,1,RGB(orate*r,
                                    orate*g,
                                    orate*b));*/ /*FOR EDGELESS*/

    hpen=CreatePen(PS_SOLID,1,RGB(r,g,b));   /*FOR PRESENTATION*/

    /*hpen=CreatePen(PS_SOLID,1,RGB(0,0,0));*/   /*FOR KESANSHO*/
  }
  else hpen=CreatePen(PS_SOLID,1,RGB(r,g,b));

  ppen=(HPEN)SelectObject(hdc,hpen);
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
    fillglobalban(hdc,vp,gb,(int)(orate*r+(1.0-orate)*255.0),
                            (int)(orate*g+(1.0-orate)*255.0),
                            (int)(orate*b+(1.0-orate)*255.0),rop);
  }
  else
  {
    fillglobalban(hdc,vp,gb,(int)(orate*r),
                            (int)(orate*g),
                            (int)(orate*b),rop);
  }

  return;
}/*drawglobalban*/

void drawelement(HDC hdc,struct viewparam vp,
				 struct oelem ge,
				 int mode,int idraw)
/*DRAW ELEMENT ON SCREEN BY AXONOMETRICS,PERSPECTIVE.*/
{
  HPEN hpen,ppen;
  DWORD rop; /*RASTER OPERATION.*/
  SIZE size;
  char str[20]="S";
  int i,j;
  int ix1,iy1,ix2,iy2,Cx,Cy,Hx,Hy; /*PROJECTED COORDINATION.*/
  double dx=0.0,dy=0.0;
  double oblique,orate;
  double Gx,Gy,**drccos,**tdrccos;
  struct onode n0,n1,n2;
  double qyi;

  if(!idraw)
  {
	for(i=0;i<ge.nnod;i++)
	{

	  if(!insiderange(**(ge.nods+i),vp.range)) return;
	}
  }

  /*COLOR OF LINE ELEMENT*/
  if(mode==ONPRINTER)
  {
	hpen=CreatePen(PS_SOLID,2,RGB(0,0,0)); /*FOR KESANSHO*/
  }
  else if(mode==ONPREVIEW)
  {
	/*hpen=CreatePen(PS_SOLID,1,RGB(0,0,0));*/

	hpen=CreatePen(PS_SOLID,1,RGB(255-ge.color.line.r,
								  255-ge.color.line.g,
								  255-ge.color.line.b));
  }
  else hpen=CreatePen(PS_SOLID,1,RGB(ge.color.line.r,
									 ge.color.line.g,
                                     ge.color.line.b));
  ppen=(HPEN)SelectObject(hdc,hpen);

  if(ge.nnod==2)
  {
	/*drawgloballine(hdc,vp,*(*(ge.nods+0)),*(*(ge.nods+1)));*/

	nodeontoscreen(*(*(ge.nods+0)),&ix1,&iy1,vp);
	nodeontoscreen(*(*(ge.nods+1)),&ix2,&iy2,vp);
    MoveToEx(hdc,ix1,iy1,NULL);
    LineTo(hdc,ix2,iy2);

    if(!idraw)
    {
      Cx=(ix1+ix2)/2; /*CENTER*/
	  Cy=(iy1+iy2)/2;

      if(vp.vflag.ev.hinge)                        /*INITIAL HINGE.*/
      {
		SelectObject(hdc,ppen);
        DeleteObject(hpen);

        if(mode==ONPRINTER)
		{
		  hpen=CreatePen(PS_SOLID,2,RGB(0,0,0)); /*BLACK*/
        }
        else if(mode==ONPREVIEW)
		{
          hpen=CreatePen(PS_SOLID,1,RGB(0,0,0)); /*BLACK*/
        }
        else hpen=CreatePen(PS_SOLID,1,RGB(255,255,255)); /*WHITE*/
		ppen=(HPEN)SelectObject(hdc,hpen);

        if(*(ge.bonds+0)==1 ||
           *(ge.bonds+1)==1 ||
		   *(ge.bonds+2)==1 ||
           *(ge.bonds+3)==1 ||
           *(ge.bonds+4)==1 ||
           *(ge.bonds+5)==1)                       /*HINGE ON HEAD.*/
		{
          dx=ix2-ix1; dy=iy2-iy1;
          oblique=sqrt(dx*dx+dy*dy);
          if(oblique==0.0)
		  {
            Hx=ix1;
			Hy=iy1;
          }
		  else
          {
            Hx=ix1+(int)(vp.dparam.hsize*dx/oblique);
			Hy=iy1+(int)(vp.dparam.hsize*dy/oblique);
		  }
          Arc(hdc,
              (Hx-vp.dparam.hsize),(Hy-vp.dparam.hsize),
              (Hx+vp.dparam.hsize),(Hy+vp.dparam.hsize),
			  0,0,0,0);                            /*HOLLOW CIRCLE.*/
        }
        if(*(ge.bonds+ 6)==1 ||
           *(ge.bonds+ 7)==1 ||
		   *(ge.bonds+ 8)==1 ||
           *(ge.bonds+ 9)==1 ||
           *(ge.bonds+10)==1 ||
           *(ge.bonds+11)==1)                      /*HINGE ON TAIL.*/
		{
          dx=ix1-ix2; dy=iy1-iy2;
          oblique=sqrt(dx*dx+dy*dy);
		  if(oblique==0.0)
		  {
            Hx=ix2;
            Hy=iy2;
          }
		  else
          {
            Hx=ix2+(int)(vp.dparam.hsize*dx/oblique);
            Hy=iy2+(int)(vp.dparam.hsize*dy/oblique);
		  }
          Arc(hdc,
			  (Hx-vp.dparam.hsize),(Hy-vp.dparam.hsize),
			  (Hx+vp.dparam.hsize),(Hy+vp.dparam.hsize),
			  0,0,0,0);                            /*HOLLOW CIRCLE.*/
		}

		SelectObject(hdc,ppen);
		DeleteObject(hpen);
	  }
	}
	if(vp.vflag.ev.axis)
	{
	  drawwireaxis(hdc,vp,*(*(ge.nods+0)),*(*(ge.nods+1)),ge.cangle);
	}
	if(vp.vflag.ev.sectionshape)
	{                                                                      //honda
	  drawwireshape(hdc,vp,*(*(ge.nods+0)),*(*(ge.nods+1)),ge.cangle,255,255,255);
	}
	if(!strncmp(prj,"sydney",6))    /*steel weight for sydney*/     //honda
	{
	  SetTextColor(hdc,RGB(255,255,255));
	  char str[256];
	  //setfontformat(hdc,20,8,"MS Mincho",255,255,255);
	  sprintf(str,"\nsteel weight = %.3f tf",totalweight);
	  TextOut(hdc, 300,50,str,strlen(str));
	  sprintf(str,"\narch length  = %.3f m",totallength);
	  TextOut(hdc, 300,65,str,strlen(str));
	  sprintf(str,"\ncenter of gravity = (%.3f,%.3f,%.3f) ",Gcx,Gcy,Gcz);
	  TextOut(hdc, 300,80,str,strlen(str));
	  //setfontformat(hdc,15,6,"MS Mincho",255,255,255);
	} /*steel weight for sydney*/

  }
  else if(!idraw && ge.nnod>2)
  {
	for(i=0;i<ge.nnod;i++)
	{
	  nodeontoscreen(*(*(ge.nods+i)),&ix1,&iy1,vp);
	  dx+=(double)ix1/ge.nnod;
	  dy+=(double)iy1/ge.nnod;
	}

    Cx=(int)dx; /*CENTER*/
	Cy=(int)dy;
  }
  SelectObject(hdc,ppen);
  DeleteObject(hpen);

  if(mode==ONPRINTER)      rop=SRCAND;
  else if(mode==ONPREVIEW) rop=SRCAND;
  else                     rop=SRCPAINT;

  for(i=0;i<ge.nban;i++)
  {
	drawglobalban(hdc,vp,*(ge.bans+i),ge.color.line.r,
                                      ge.color.line.g,
                                      ge.color.line.b,rop);

    /*DRAW WINDOW*/
    if(ge.wrect[0]!=0.0 || ge.wrect[1]!=0.0)
    {
      SelectObject(hdc,ppen);
      DeleteObject(hpen);

	  if(mode==ONPRINTER)
      {
        orate=(wdraw.childs+1)->org.opaque;
        hpen=CreatePen(PS_SOLID,1,
                       RGB((int)(orate*ge.color.line.r),
                           (int)(orate*ge.color.line.g),
                           (int)(orate*ge.color.line.b)));
	  }
      else if(mode==ONPREVIEW)
      {
        orate=(wdraw.childs+1)->org.opaque;
        hpen=CreatePen(PS_SOLID,1,
                       RGB((int)(orate*ge.color.line.r),
                           (int)(orate*ge.color.line.g),
						   (int)(orate*ge.color.line.b)));
      }
      else
      {
        hpen=CreatePen(PS_SOLID,1,RGB(ge.color.line.r,
                                      ge.color.line.g,
                                      ge.color.line.b));
	  }
      ppen=(HPEN)SelectObject(hdc,hpen);

	  drccos=filmdrccos(*(*((ge.bans+i)->nods+0)),
						*(*((ge.bans+i)->nods+1)),
						*(*((ge.bans+i)->nods+(ge.bans+i)->nnod-1)));
                                                /*DIRECTION COSINE.*/
	  tdrccos=matrixtranspose(drccos,3);              /*MATRIX [T].*/

      /*CENTER DOT*/
      Gx=0.0; Gy=0.0;
      for(j=0;j<(ge.bans+i)->nnod;j++)
      {
        n0=transgtol(drccos,*(*((ge.bans+i)->nods+0)),
							*(*((ge.bans+i)->nods+j)));

        Gx+=n0.d[EX];
        Gy+=n0.d[EY];
      }
      Gx/=(ge.bans+i)->nnod;
      Gy/=(ge.bans+i)->nnod;

	  n1=setcoord(0.0,Gx-0.5*ge.wrect[0],Gy-0.5*ge.wrect[1]);
      n2=setcoord(0.0,Gx+0.5*ge.wrect[0],Gy-0.5*ge.wrect[1]);
      drawlocalline(hdc,vp,tdrccos,*(*((ge.bans+i)->nods+0)),n1,n2);
      n1=setcoord(0.0,Gx+0.5*ge.wrect[0],Gy-0.5*ge.wrect[1]);
      n2=setcoord(0.0,Gx+0.5*ge.wrect[0],Gy+0.5*ge.wrect[1]);
      drawlocalline(hdc,vp,tdrccos,*(*((ge.bans+i)->nods+0)),n1,n2);
	  n1=setcoord(0.0,Gx+0.5*ge.wrect[0],Gy+0.5*ge.wrect[1]);
      n2=setcoord(0.0,Gx-0.5*ge.wrect[0],Gy+0.5*ge.wrect[1]);
      drawlocalline(hdc,vp,tdrccos,*(*((ge.bans+i)->nods+0)),n1,n2);
      n1=setcoord(0.0,Gx-0.5*ge.wrect[0],Gy+0.5*ge.wrect[1]);
      n2=setcoord(0.0,Gx-0.5*ge.wrect[0],Gy-0.5*ge.wrect[1]);
      drawlocalline(hdc,vp,tdrccos,*(*((ge.bans+i)->nods+0)),n1,n2);

	  for(j=2;j>=0;j--) free(*(drccos+j));
      free(drccos);
      for(j=2;j>=0;j--) free(*(tdrccos+j));
      free(tdrccos);

      SelectObject(hdc,ppen);
      DeleteObject(hpen);
	}

    if((ge.type==SLAB || ge.type==WALL) &&
       vp.vflag.ev.cmqline)
    {
      SelectObject(hdc,ppen);
	  DeleteObject(hpen);

      if(mode==ONPRINTER)
      {
        orate=(wdraw.childs+1)->org.opaque;
        hpen=CreatePen(PS_SOLID,1,
					   RGB((int)(orate*ge.color.line.r),
                           (int)(orate*ge.color.line.g),
                           (int)(orate*ge.color.line.b)));
      }
      else if(mode==ONPREVIEW)
      {
		orate=(wdraw.childs+1)->org.opaque;
		hpen=CreatePen(PS_SOLID,1,
					   RGB((int)(orate*ge.color.line.r),
						   (int)(orate*ge.color.line.g),
						   (int)(orate*ge.color.line.b)));
	  }
	  else
	  {
        hpen=CreatePen(PS_SOLID,1,RGB(ge.color.line.r,
									  ge.color.line.g,
									  ge.color.line.b));
	  }

	  ppen=(HPEN)SelectObject(hdc,hpen);

      slabdivision2(hdc,&vp,*(ge.bans+i),NULL);

      SelectObject(hdc,ppen);
	  DeleteObject(hpen);
	}
  }

  if(!idraw)
  {
	if(mode==ONPRINTER)      SetTextColor(hdc,RGB(0,0,0)); /*BLACK*/
	else if(mode==ONPREVIEW) SetTextColor(hdc,RGB(0,0,0)); /*BLACK*/
	else SetTextColor(hdc,RGB(ge.color.line.r,
							  ge.color.line.g,
							  ge.color.line.b));
	GetTextExtentPoint32(hdc,str,strlen(str),&size);

	if(vp.vflag.ev.code) /*ELEM CODE*/
	{
	  sprintf(str,"%d",ge.code);
	  TextOut(hdc,Cx,Cy,str,strlen(str));
      Cy+=size.cy;
	}
	if(vp.vflag.ev.sectioncode) /*SECTION CODE.*/
    {
	  sprintf(str,"%d",ge.sect->code);
	  TextOut(hdc,Cx,Cy,str,strlen(str));
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
  ppen=(HPEN)SelectObject(hdc,hpen);

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
  SIZE size;
  char str[256];
  int i,flag,nx,ny;
  long int loff;
  long int imax[2],imin[2];
  double hiju,orate,radius,hsize;

  /*PAGE TITLE*/
  if(vp.vflag.mv.pagetitle)
  {
    if(mode==ONPRINTER)
	{
	  setfontformat(hdc,150,50,"MS Mincho",0,0,0);
    }
    else if(mode==ONPREVIEW)
	{
	  setfontformat(hdc,150,50,"MS Mincho",0,0,0);
    }

    if(mode==ONPRINTER)      SetTextColor(hdc,RGB(0,0,0));
    else if(mode==ONPREVIEW) SetTextColor(hdc,RGB(0,0,0));
    else                     SetTextColor(hdc,RGB(255,255,255));

    sprintf(str,"%s",(wdraw.childs+1)->pagetitle);
//    if(mode==ONPRINTER) TextOut(hdc,1200,500,str,strlen(str));
	if(mode==ONPRINTER) TextOut(hdc, 700,300,str,strlen(str));
	else                TextOut(hdc, 120, 50,str,strlen(str));
  }

  /*PAGE NUMBER*/
  if(vp.vflag.mv.pagenum)
  {
	if(mode==ONPRINTER)
	{
	  setfontformat(hdc,80,30,"MS Mincho",0,0,0);
	}
	else if(mode==ONPREVIEW)
	{
	  setfontformat(hdc,80,30,"MS Mincho",0,0,0);
	}
	if(mode==ONPRINTER)      SetTextColor(hdc,RGB(0,0,0));
	else if(mode==ONPREVIEW) SetTextColor(hdc,RGB(0,0,0));
	else                     SetTextColor(hdc,RGB(255,255,255));

	sprintf(str,"%d.%d.%d",(wdraw.childs+1)->vparam.chapter,
						   (wdraw.childs+1)->vparam.section,
                           (wdraw.childs+1)->vparam.subsection);

	if(mode==ONPRINTER) TextOut(hdc,4400,6600,str,strlen(str));   //要調整
	else if(mode==ONPREVIEW) TextOut(hdc,1200,1000,str,strlen(str));
	else                     TextOut(hdc,1250,1000,str,strlen(str));
  }

  /*GET FRAME RANGE FOR DRAWING TEXTS*/
  if(mode==ONPRINTER || mode==ONPREVIEW)
  {
	setfontformat(hdc,(int)(1.5*gprn.jiheight),
					  (int)(1.5*gprn.jiwidth),"MS Mincho",0,0,0);
  }
  if(mode==ONPRINTER)      SetTextColor(hdc,RGB(0,0,0));
  else if(mode==ONPREVIEW) SetTextColor(hdc,RGB(0,0,0));
  else                     SetTextColor(hdc,RGB(255,255,255));

  imax[0]=0; imin[0]=0;
  imax[1]=0; imin[1]=0;

  flag=0;
  for(i=0;i<go.nnode;i++)
  {
	if(insiderange(*(go.nodes+i),vp.range))
	{
	  if(!nodeontoscreen(*(go.nodes+i),&nx,&ny,vp)) return;

      if(flag==0) /*INITIAL*/
      {
        imax[0]=nx; imin[0]=nx;
        imax[1]=ny; imin[1]=ny;

        flag=1;
      }

	  if(imax[0]<nx) imax[0]=nx;
      if(imax[1]<ny) imax[1]=ny;
      if(imin[0]>nx) imin[0]=nx;
      if(imin[1]>ny) imin[1]=ny;
    }
  }
  sprintf(str,"Text");
  GetTextExtentPoint32(hdc,str,strlen(str),&size);

  /*UPPER TEXTS*/
  imin[1]-=size.cy;
  if(vp.vflag.mv.title)
  {
    imin[1]-=size.cy;
    sprintf(str,"%s",(wdraw.childs+1)->title);
    TextOut(hdc,imin[0],imin[1],str,strlen(str));
  }

  /*LOWER TEXTS*/
  imax[1]+=size.cy;
  if(vp.vflag.nv.code)
  {
    imax[1]+=size.cy;
    sprintf(str,"節点番号");
    TextOut(hdc,imin[0],imax[1],str,strlen(str));
  }
/*mihara for zahyo 20190329///////////////////////////////////////////////////*/
  if(vp.vflag.nv.d[0])
  {
    imax[1]+=size.cy;
	sprintf(str,"X座標[m]");
    TextOut(hdc,imin[0],imax[1],str,strlen(str));
  }
  if(vp.vflag.nv.d[1])
  {
	imax[1]+=size.cy;
	sprintf(str,"Y座標[m]");
    TextOut(hdc,imin[0],imax[1],str,strlen(str));
  }
  if(vp.vflag.nv.d[2])
  {
    imax[1]+=size.cy;
	sprintf(str,"Z座標[m]");
    TextOut(hdc,imin[0],imax[1],str,strlen(str));
  }
/*mihara for zahyo 20190329///////////////////////////////////////////////////*/

  if(vp.vflag.ev.code)
  {
    imax[1]+=size.cy;
	sprintf(str,"部材番号");
	TextOut(hdc,imin[0],imax[1],str,strlen(str));
  }
  if(vp.vflag.ev.sectioncode)
  {
    imax[1]+=size.cy;
	sprintf(str,"断面番号");
    TextOut(hdc,imin[0],imax[1],str,strlen(str));
  }
  if(vp.vflag.nv.mcircle)
  {
    imax[1]+=size.cy;
    sprintf(str,"節点重量図");
    TextOut(hdc,imin[0],imax[1],str,strlen(str));
  }
  if(vp.vflag.nv.mvalue)
  {
    imax[1]+=size.cy;
    if(globalunit==1.0)         sprintf(str,"節点重量値 [tf]");
    else if(globalunit==SIUNIT) sprintf(str,"節点重量値 [kN]");
    else                        sprintf(str,"節点重量値");
    TextOut(hdc,imin[0],imax[1],str,strlen(str));
  }
  if(vp.vflag.ev.axis)
  {
    imax[1]+=size.cy;
    sprintf(str,"部材座標軸");
    TextOut(hdc,imin[0],imax[1],str,strlen(str));
  }
  if(vp.vflag.ev.cmqline)
  {
    imax[1]+=size.cy;
    sprintf(str,"分割線表示");
	TextOut(hdc,imin[0],imax[1],str,strlen(str));
  }
  if(vp.vflag.ev.etype[0] ||
     vp.vflag.ev.etype[1] ||
     vp.vflag.ev.etype[2] ||
     vp.vflag.ev.etype[3] ||
     vp.vflag.ev.etype[4] ||
     vp.vflag.ev.etype[5] ||
     vp.vflag.ev.etype[6] ||
	 vp.vflag.ev.etype[7])
  {
    flag=0;
    imax[1]+=size.cy;
    sprintf(str,"表示部材：");
	/*if(vp.vflag.ev.etype[0])
    {
      strcat(str,"補助");
      flag=1;
	}*/
    if(vp.vflag.ev.etype[1])
    {
      if(flag) strcat(str,", ");
	  strcat(str,"柱");
      flag=1;
    }
    if(vp.vflag.ev.etype[2])
	{
      if(flag) strcat(str,", ");
      strcat(str,"梁");
      flag=1;
    }
    /*if(vp.vflag.ev.etype[3])
    {
	  if(flag) strcat(str,", ");
      strcat(str,"小梁");
      flag=1;
    }*/
    if(vp.vflag.ev.etype[4])
    {
      if(flag) strcat(str,", ");
	  strcat(str,"ブレース");
      flag=1;
    }
    if(vp.vflag.ev.etype[5])
    {
      if(flag) strcat(str,", ");
      strcat(str,"壁");
      flag=1;
    }
	if(vp.vflag.ev.etype[6])
    {
      if(flag) strcat(str,", ");
      strcat(str,"床");
      flag=1;
    }
    /*if(vp.vflag.ev.etype[7])
    {
      if(flag) strcat(str,", ");
	  strcat(str,"張壁");
      flag=1;
    }*/
    TextOut(hdc,imin[0],imax[1],str,strlen(str));
  }
  imax[1]+=size.cy;
  if(vp.vflag.mv.inputfile)
  {
    imax[1]+=size.cy;
	sprintf(str,"Input  : %s",(wdraw.childs+1)->inpfile);
    TextOut(hdc,imin[0],imax[1],str,strlen(str));
  }
  if(vp.vflag.mv.outputfile)
  {
    imax[1]+=size.cy;
    sprintf(str,"Output : %s",(wdraw.childs+1)->otpfile);
    TextOut(hdc,imin[0],imax[1],str,strlen(str));
  }
  if(vp.vflag.mv.view)
  {
    imax[1]+=size.cy;
    sprintf(str,"Focus : %.3f %.3f %.3f",
            vp.focus.d[GX],vp.focus.d[GY],vp.focus.d[GZ]);
    TextOut(hdc,imin[0],imax[1],str,strlen(str));
    imax[1]+=size.cy;
    sprintf(str,"Phi=%.3f Theta=%.3f",vp.phi,vp.theta);
    TextOut(hdc,imin[0],imax[1],str,strlen(str));
    imax[1]+=size.cy;
    sprintf(str,"R=%.3f L=%.3f",vp.r,vp.odv);
    TextOut(hdc,imin[0],imax[1],str,strlen(str));
    imax[1]+=size.cy;
	sprintf(str,"Dfact=%.3f",vp.dparam.dfact);
    TextOut(hdc,imin[0],imax[1],str,strlen(str));
    imax[1]+=size.cy;
    sprintf(str,"Mfact=%.3f",vp.dparam.mfact);
    TextOut(hdc,imin[0],imax[1],str,strlen(str));
  }

  if(vp.vflag.nv.conffig)
  {
    imax[1]+=size.cy;
    imax[1]+=size.cy;
    sprintf(str,"▲印は支点位置を表す");
	TextOut(hdc,imin[0],imax[1],str,strlen(str));
   imax[1]+=size.cy;
    sprintf(str,"(X,Y,Z,θx,θy,θzのうち，X,Y,Zを拘束している)");
    TextOut(hdc,imin[0],imax[1],str,strlen(str));

//    sprintf(str,"■印は支点位置を表す");
//    TextOut(hdc,imin[0],imax[1],str,strlen(str));
//    imax[1]+=size.cy;
//    sprintf(str,"(X,Y,Z,θx,θy,θzのうち，X,Y,Z,θx,θy,θzを拘束している)");
//    TextOut(hdc,imin[0],imax[1],str,strlen(str));

    /*FOR OUTPUT REACTION LONG*/
/*    imax[1]+=size.cy;
    imax[1]+=size.cy;
    sprintf(str,"数値の凡例");
    TextOut(hdc,imin[0],imax[1],str,strlen(str));
    imax[1]+=size.cy;
	sprintf(str,"上段　Ｚ方向反力 [kN]");
    TextOut(hdc,imin[0],imax[1],str,strlen(str));
    imax[1]+=size.cy;
    sprintf(str,"中段　支点重量　 [kN]");
    TextOut(hdc,imin[0],imax[1],str,strlen(str));
    imax[1]+=size.cy;
    sprintf(str,"下段　合計　　　 [kN]");
	TextOut(hdc,imin[0],imax[1],str,strlen(str));
*/

    /*FOR OUTPUT REACTION SHORT*/
/*    imax[1]+=size.cy;
    imax[1]+=size.cy;
    sprintf(str,"数値の凡例");
    TextOut(hdc,imin[0],imax[1],str,strlen(str));
    imax[1]+=size.cy;
    sprintf(str,"Ｚ方向反力　[kN]");
    TextOut(hdc,imin[0],imax[1],str,strlen(str));
*/
  }

  /*FONT*/
  if(mode==ONPRINTER || mode==ONPREVIEW)
  {
	setfontformat(hdc,(int)gprn.jiheight,
                      (int)gprn.jiwidth,"MS Mincho",0,0,0);
  }

  if(vp.vflag.axis==1) /*GLOBAL AXIS.*/
  {
    if(mode==ONPRINTER)      drawglobalaxis(hdc,vp,0,0,0);
	else if(mode==ONPREVIEW) drawglobalaxis(hdc,vp,0,0,0);
    else                     drawglobalaxis(hdc,vp,0,0,255);
  }

  if(mode==ONPRINTER)
  {
    hsize=vp.dparam.hsize;
    vp.dparam.hsize*=5.0;

//    csize=vp.dparam.csize;           /*by MIHARA for conffig*/
    vp.dparam.csize*=5.0;            /*by MIHARA for conffig*/
  }
  else if(mode==ONPREVIEW)
  {
    hsize=vp.dparam.hsize;
	vp.dparam.hsize*=5.0;

//    csize=vp.dparam.csize;           /*by MIHARA for conffig*/
    vp.dparam.csize*=5.0;            /*by MIHARA for conffig*/
  }
  else
  {
    hsize=vp.dparam.hsize;
    vp.dparam.hsize*=0.5; /*FOR INAX*/
  }

//Modifyed by MIHARA for conffig 2007.07.11/////////////////////////////////////

  if(vp.vflag.nv.conffig)
  {
    for(i=0;i<go.nnode;i++) /*CINFFIG*/
   {
        loff=6*((go.nodes+i)->loff);
        drawglobalconfig(hdc,vp,go,*(go.nodes+i),mode,loff);
   }
  }

//printrange 2007.12.04/////////////////////////////////////////////////////////
  if(vp.vflag.pv.printrange)
  {
   if((mode==ONPRINTER)||(mode==ONPREVIEW))
   {}
   else
   {
   drawprintrange(hdc,vp);
   }
  }
////////////////////////////////////////////////////////////////////////////////

  for(i=0;i<go.nelem;i++) /*ELEMENTS.*/
  {
    if((go.elems+i)->sect->dflag==1 &&
       vp.vflag.ev.etype[(go.elems+i)->type]==1)
    {
	  (go.elems+i)->color.line.r=(go.elems+i)->sect->dcolor.r;
      (go.elems+i)->color.line.g=(go.elems+i)->sect->dcolor.g;
      (go.elems+i)->color.line.b=(go.elems+i)->sect->dcolor.b;

      drawelement(hdc,vp,*(go.elems+i),mode,0);
    }
  }
  if(mode==ONPRINTER)      vp.dparam.hsize=hsize;
  else if(mode==ONPREVIEW) vp.dparam.hsize=hsize;

  orate=go.opaque; /*OPAQUE IS ALWAYS UPDATED ON DIALOG INPUT*/
  getmasshiju((wmenu.childs+4)->hwnd,&hiju);

  /*
  if(mode==ONPRINTER)      SetTextColor(hdc,RGB(0,0,0));
  else if(mode==ONPREVIEW) SetTextColor(hdc,RGB(0,0,0));
  else                     SetTextColor(hdc,RGB(150,150,255));
  */

  if(vp.vflag.nv.code ||
     vp.vflag.nv.mcircle ||
     vp.vflag.nv.mvalue)
  {
    for(i=0;i<go.nnode;i++) /*NODES*/
    {
      if(insiderange(*(go.nodes+i),vp.range))
	  {
        loff=6*((go.nodes+i)->loff);

        if(mode==ONPRINTER)      SetTextColor(hdc,RGB(0,0,0));
        else if(mode==ONPREVIEW) SetTextColor(hdc,RGB(0,0,0));
        else if((go.confs+(loff+3))->iconf==1 ||
                (go.confs+(loff+4))->iconf==1 ||
                (go.confs+(loff+5))->iconf==1) /*Tx,Ty,Tz FIXED*/
        {
		  SetTextColor(hdc,RGB(255,0,150));
        }
        else if((go.confs+(loff+0))->iconf==1 ||
                (go.confs+(loff+1))->iconf==1 ||
                (go.confs+(loff+2))->iconf==1) /*X,Y,Z FIXED*/
        {
          SetTextColor(hdc,RGB(150,255,0));
        }
        else SetTextColor(hdc,RGB(150,150,255));

        drawglobalnode(hdc,vp,*(go.nodes+i),go.loads);

        if(go.loads!=NULL) /*DRAW WEIGHT*/
        {
          if(vp.vflag.nv.mcircle)
          {
            radius=pow(3.0/4.0/PI*((go.loads+i)->w[WEIGHTSLAB]/hiju),
                       1.0/3.0);

            if(mode==ONPRINTER)
            {
              fillglobalcircle(hdc,vp,*((go.loads+i)->nod),radius,
                               0,0,0,orate,SRCAND);
            }
            else if(mode==ONPREVIEW)
            {
              fillglobalcircle(hdc,vp,*((go.loads+i)->nod),radius,
							   0,0,0,orate,SRCAND);
            }
            else
            {
              fillglobalcircle(hdc,vp,*((go.loads+i)->nod),radius,
                               150,150,150,orate,SRCPAINT);
            }
          }
        }
	  }
    }
  }

#if 0
  /*Draw Bezier-Surface Controle Points*/
  FILE *ftext;
  char **data/*,str[256]="\0"*/;
  double ddata;
  int ndata;
  long int nnode,ncode;

  double *x,*y,*z;   /*Controle points*/
  int m,n;           /*Degree*/
  int ncontrole;

  double *u,*v;      /*U-V Coordinates*/
  struct onode node;

  /*OPEN FILE*/
  ftext=fopen("bezier.txt","r");   /*bezier-surface data*/
  if(ftext==NULL)
  {
	errormessage("ACCESS IMPOSSIBLE.");
	return;
  }
  fseek(ftext,0L,SEEK_SET);

  data=fgetsbrk(ftext,&ndata);
  m=strtol(*(data+0),NULL,10);
  n=strtol(*(data+1),NULL,10);
  nnode=strtol(*(data+2),NULL,10);

  /*Controle points*/  ncontrole=(m+1)*(n+1);
  x=(double *)malloc((ncontrole)*sizeof(double));
  y=(double *)malloc((ncontrole)*sizeof(double));
  z=(double *)malloc((ncontrole)*sizeof(double));

  for(i=0;i<ncontrole;i++)
  {
	data=fgetsbrk(ftext,&ndata);
	if(ndata!=3) return;
	x[i]=strtod(*(data+0),NULL);
	y[i]=strtod(*(data+1),NULL);
	z[i]=strtod(*(data+2),NULL);

	for(;ndata>0;ndata--) free(*(data+ndata-1));
	free(data);
  }

  fclose(ftext);

  /*DRAW CONTROLE POINTS*/
  for(i=0;i<ncontrole;i++)
  {
	node.d[0]=x[i];
	node.d[1]=y[i];
	node.d[2]=z[i];
	drawcontrolepoint(hdc,vp,node);
  }
#endif

  return;
}/*draworganization*/

void draworgannodes(HDC hdc,struct viewparam vp,struct organ go)
/*DRAW NODES OF ORGANIZATION ON SCREEN.*/
{
  char str[20];
  int i,Ox,Oy;

  drawglobalaxis(hdc,vp,0,0,255);               /*DRAW GLOBAL AXIS.*/

  SetTextColor(hdc,RGB(0,255,255));
  for(i=0;i<go.nnode;i++)
  {
    /*drawglobalnode(hdc,vp,*(go.nodes+i),NULL);*/

    if(!nodeontoscreen(*(go.nodes+i),&Ox,&Oy,vp)) return;
    sprintf(str,"%ld",(go.nodes+i)->code);
    TextOut(hdc,Ox+2,Oy,str,strlen(str));
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

void drawtextalignedonlocaldot(HDC hdc,struct viewparam vp,
                               double **tdrccos,
                               struct onode on, /*ORIGIN*/
                               struct onode ln, /*DOT*/
                               char *str,
                               int wtype,int htype)
/*DRAW GLOBAL TEXT ALIGNED ON LOCAL NODE PROJECTED.*/
{
  struct onode gn;

  gn=transltog(tdrccos,on,ln);
  drawglobaltextaligned(hdc,vp,gn,str,wtype,htype);

  return;
}/*drawtextalignedonlocaldot*/

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
  /*ppen=(HPEN)SelectObject(hdc,hpen);*/
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

  if(distance[2]>=distance[0] && distance[2]>=distance[1]
          && (wdraw.childs+1)->vparam.phi!=90.0
          && ((wdraw.childs+1)->vparam.phi!=270.0
             || (wdraw.childs+1)->vparam.phi!=-90.0))     /*Modified by Ujioka for 2D MODEL*/
  {
    *length=coordfactor[2]*arrowlength;
    return DIRECTIONZ;
  }
  else if(distance[1]>=distance[2] && distance[1]>=distance[0]
          && (wdraw.childs+1)->vparam.phi!=0.0
          && ((wdraw.childs+1)->vparam.phi!=270.0
             || (wdraw.childs+1)->vparam.phi!=-90.0))     /*Modified by Ujioka for 2D MODEL*/
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
                    struct arclmframe af,long int code,
                    int mode)
/*DRAW NODES OF ARCLM FRAME.*/
/*CODE,CONFINEMENT,DISPLACEMENT,REACTION.*/
{
  UINT align;
  char str[20];
  char flabel[6][10]={"Fx","Fy","Fz","Mx","My","Mz"};       /*LOAD.*/
  char clabel[6][10]={"Dx","Dy","Dz","Rx","Ry","Rz"}; /*COMPULSORY.*/
  char dlabel[3][10]={"Dx","Dy","Dz"};              /*DISPLACEMENT.*/
  char rlabel[6][10]={"Fx","Fy","Fz","Mx","My","Mz"};   /*REACTION.*/
  long int i,j,loff;
  long int nreact=0;
  double pitch=vp.dparam.pitch;
  double value,radius,hiju,orate;
  struct onode tn; /*TEXT POSITION.*/

  if(af.nmass!=NULL && vp.vflag.nv.mcircle)
  {
    getdoublefromdialog((wmenu.childs+4)->hwnd,IDVS_OPAQUE,&orate);
    getmasshiju((wmenu.childs+4)->hwnd,&hiju);
  }

  for(i=0;i<af.nnode;i++)
  {
	loff=6*((af.ninit+i)->loff);

    if(insiderange(*(af.ninit+i),vp.range))
    {
      if(vp.vflag.nv.code) /*CODES.*/
      {
        if(mode==ONPRINTER)      SetTextColor(hdc,RGB(0,0,0));
        else if(mode==ONPREVIEW) SetTextColor(hdc,RGB(0,0,0));
        else if(globalstatus==SELECTNODE && (af.nodes+i)->code==code)
        {
          SetTextColor(hdc,RGB(0,255,255));
        }
        else SetTextColor(hdc,RGB(150,150,255));

        drawglobalnode(hdc,vp,*(af.ninit+i),NULL);
      }

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
            if(mode==ONPRINTER)      SetTextColor(hdc,RGB(  0,  0,  0));
            else if(mode==ONPREVIEW) SetTextColor(hdc,RGB(  0,  0,  0));
            else                     SetTextColor(hdc,RGB(255,  0,255));
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
            if(mode==ONPRINTER)      SetTextColor(hdc,RGB(  0,  0,  0));
            else if(mode==ONPREVIEW) SetTextColor(hdc,RGB(  0,  0,  0));
            else                     SetTextColor(hdc,RGB(255,100,  0));
            drawglobaltext(hdc,vp,tn,str);                  /*TEXT.*/
          }
        }
      }

      align = SetTextAlign(hdc,TA_BOTTOM |
                               TA_LEFT |
                               TA_NOUPDATECP);
      for(j=2;j>=0;j--) /*DISPLACEMETS.*/
      {
        if(vp.vflag.nv.disps[j] &&
           (!(af.confs+loff+j)->iconf) &&
           af.ddisp!=NULL)
        {
		  value=*(af.ddisp+loff+j);
          value-=(af.ninit+i)->d[j];

          /*sprintf(str,"%s %.5f",dlabel[j],value);*/ /*[m]*/
          sprintf(str," %6.3f",value*100.0); /*[cm]*/

          tn.d[GZ]+=pitch;
          if(mode==ONPRINTER)      SetTextColor(hdc,RGB(  0,  0,  0));
          else if(mode==ONPREVIEW) SetTextColor(hdc,RGB(  0,  0,  0));
          else                     SetTextColor(hdc,RGB(255,255,  0));
          drawglobaltext(hdc,vp,tn,str);                    /*TEXT.*/
        }
      }
      SetTextAlign(hdc,align);

      if(af.dreact!=NULL)
      {
        tn=*(af.ninit+i);
        for(j=0;j<=5;j++) /*REACTIONS.*/
		{
		  if((af.confs+loff+j)->iconf==1)
		  {
			value=*(af.dreact+nreact);
			nreact++;

			if(vp.vflag.nv.react[j])
			{
			  sprintf(str,"%s %.5f",rlabel[j],value);

			  tn.d[GZ]-=pitch;
			  if(mode==ONPRINTER)      SetTextColor(hdc,
													RGB(  0,  0,  0));
              else if(mode==ONPRINTER) SetTextColor(hdc,
                                                    RGB(  0,  0,  0));
              else                     SetTextColor(hdc,
                                                    RGB(255,  0,  0));
              drawglobaltext(hdc,vp,tn,str);                /*TEXT.*/
            }
          }
        }
	  }

      /*MASS CIRCLE*/
      if(af.nmass!=NULL) /*DRAW WEIGHT*/
      {
        if(vp.vflag.nv.mcircle)
        {
          radius=pow(3.0/4.0/PI*(*(af.nmass+i)/hiju),
                     1.0/3.0);

          if(mode==ONPRINTER)
          {
            fillglobalcircle(hdc,vp,*(af.nodes+i),radius,
                             0,0,0,orate,SRCAND);
          }
          else if(mode==ONPREVIEW)
          {
            fillglobalcircle(hdc,vp,*(af.nodes+i),radius,
                             0,0,0,orate,SRCAND);
          }
          else
          {
            fillglobalcircle(hdc,vp,*(af.nodes+i),radius,
                             150,150,150,orate,SRCPAINT);
          }
        }
        if(vp.vflag.nv.mvalue)
        {
          sprintf(str,"%6.3f",*(af.nmass+i)); /*[]*/
          if(mode==ONPRINTER)      SetTextColor(hdc,RGB(  0,  0,  0));
          else if(mode==ONPREVIEW) SetTextColor(hdc,RGB(  0,  0,  0));
          else                     SetTextColor(hdc,RGB(255,255,255));
          drawglobaltext(hdc,vp,*(af.nodes+i),str);         /*TEXT.*/
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
					long int selectcode,
					int mode)
/*DRAW ELEMENT,PLASTIC HINGE BY AXONOMETRICS,PERSPECTIVE.*/
{
  HPEN hpen,ppen;                                  /*HANDLE OF PEN.*/
  HBRUSH hbrush,pbrush;                          /*HANDLE OF BRUSH.*/
  SIZE size; /****SRCANMAX****/
  char str[10];
  int ix1,iy1,ix2,iy2,dx1,dy1,dx2,dy2;    /*PROJECTED COORDINATION.*/
  long int loff;
  double dx,dy;
  int hingeX,hingeY;                             /*CENTER OF HINGE.*/
  double oblique,radius,hiju,orate,etotal;          /*OBLIQUE SIDE.*/
  struct onode inode1,inode2;                       /*INITIAL NODE.*/
  struct onode dnode1,dnode2;                      /*DEFORMED NODE.*/
  int i,ndiv;
  double rate;

  initialnode(af.ninit,af.nnode,(elem.node[0]->code),&inode1);
  if(!insiderange(inode1,vp.range)) return;
  if(af.ddisp==NULL || !vp.vflag.ev.deformation)
  {
	dnode1.d[GX]=0.0;
	dnode1.d[GY]=0.0;
	dnode1.d[GZ]=0.0;
  }
  else                                              /*DISPLACEMENT.*/
  {
	/*inputnode(af.ddisp,elem.node[0]);*/ /*DEFORMED NODE.*/
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
  if(af.ddisp==NULL || !vp.vflag.ev.deformation)
  {
	dnode2.d[GX]=0.0;
    dnode2.d[GY]=0.0;
    dnode2.d[GZ]=0.0;
  }
  else                                              /*DISPLACEMENT.*/
  {
    /*inputnode(af.ddisp,elem.node[0]);*/ /*DEFORMED NODE.*/
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
    else if(mode==ONPRINTER)
    {
	  SetTextColor(hdc,RGB(0,0,0)); /*BLACK*/
	}
	else if(mode==ONPREVIEW)
	{
	  SetTextColor(hdc,RGB(0,0,0)); /*BLACK*/
	}
	else
	{
	  loff=elem.node[0]->loff;
	  if((af.confs+6*loff+3)->iconf==1 ||
		 (af.confs+6*loff+4)->iconf==1 ||
		 (af.confs+6*loff+5)->iconf==1) /*FIXED*/
	  {
		SetTextColor(hdc,RGB(255,150,150));
      }
      else if((af.confs+6*loff+0)->iconf==1 ||
              (af.confs+6*loff+1)->iconf==1 ||
              (af.confs+6*loff+2)->iconf==1) /*HINGE*/
	  {
        SetTextColor(hdc,RGB(150,255,150));
      }
      else SetTextColor(hdc,RGB(150,150,255)); /*FREE*/
    }
    sprintf(str,"%d",elem.node[0]->code);
	TextOut(hdc,ix1+2,iy1,str,strlen(str));

    if(globalstatus==SELECTNODE &&
       elem.node[1]->code==selectcode)
    {
      SetTextColor(hdc,RGB(0,255,255));
	}
    else if(mode==ONPRINTER)
    {
      SetTextColor(hdc,RGB(0,0,0)); /*BLACK*/
    }
    else if(mode==ONPREVIEW)
	{
      SetTextColor(hdc,RGB(0,0,0)); /*BLACK*/
    }
    else
    {
      loff=elem.node[1]->loff;
	  if((af.confs+6*loff+3)->iconf==1 ||
         (af.confs+6*loff+4)->iconf==1 ||
         (af.confs+6*loff+5)->iconf==1) /*FIXED*/
      {
        SetTextColor(hdc,RGB(255,150,150));
      }
	  else if((af.confs+6*loff+0)->iconf==1 ||
              (af.confs+6*loff+1)->iconf==1 ||
              (af.confs+6*loff+2)->iconf==1) /*HINGE*/
      {
        SetTextColor(hdc,RGB(150,255,150));
      }
	  else SetTextColor(hdc,RGB(150,150,255)); /*FREE*/
    }

    sprintf(str,"%d",elem.node[1]->code);
    TextOut(hdc,ix2+2,iy2,str,strlen(str));
  }
  if(vp.vflag.ev.code)
  {
    SetTextColor(hdc,RGB(cred,cgreen,cblue)); /*ELEM CODE.*/
	/****SRCANMAX****/
	if(vp.vflag.ev.srcanmax) sprintf(str,"部材 %d",elem.code);
	else                     sprintf(str,"%d",elem.code);
	if(vp.vflag.ev.srcanmax && mode==ONSRCANMAX)
	{
	  setfontformat(hdc,(int)(1.2*gprn.jiheight),
						(int)(1.2*gprn.jiwidth),"ＭＳ 明朝",0,0,0);
	}
	//sprintf(str,"%d",elem.code);
	TextOut(hdc,(0.55*(ix2-ix1)+ix1),
                (0.55*(iy2-iy1)+iy1),str,strlen(str));
	/****************/
//    TextOut(hdc,(0.3*(ix2-ix1)+ix1),                           //mihara
//                (0.3*(iy2-iy1)+iy1),str,strlen(str));
  }
  if(vp.vflag.ev.sectioncode) /*SECTION CODE.*/
  {
    /*SetTextColor(hdc,RGB(elem.sect->dcolor.r,
                         elem.sect->dcolor.g,
                         elem.sect->dcolor.b));*/
    if(mode==ONPRINTER)
	{
      SetTextColor(hdc,RGB(0,0,0)); /*BLACK*/
    }
    else if(mode==ONPREVIEW)
    {
      SetTextColor(hdc,RGB(0,0,0)); /*BLACK*/
	}
    else SetTextColor(hdc,RGB(150,150,150));

	/****SRCANMAX****/
	if(vp.vflag.ev.srcanmax) sprintf(str,"断面 %d",elem.sect->ocode);
	else                     sprintf(str,"%d",elem.sect->code);
	if(vp.vflag.ev.srcanmax && mode==ONSRCANMAX)
	{
	  setfontformat(hdc,(int)(1.2*gprn.jiheight),
						(int)(1.2*gprn.jiwidth),"ＭＳ 明朝",0,0,0);
	}
	if(vp.vflag.ev.code)
	{
	  if(vp.vflag.ev.srcanmax)
	  {
		if(mode==ONSRCANMAX) SetTextColor(hdc,RGB(0,0,0));
		else                 SetTextColor(hdc,RGB(255,0,255));
	  }
	  GetTextExtentPoint32(hdc,str,strlen(str),&size);
	  TextOut(hdc,(0.55*(ix2-ix1)+ix1)+4,
				  (0.55*(iy2-iy1)+iy1)+4+size.cy,str,strlen(str));
	}
	else
	{
	  TextOut(hdc,(0.45*(ix2-ix1)+ix1),
				  (0.45*(iy2-iy1)+iy1),str,strlen(str));
	}
	/****************/
	//sprintf(str,"%d",elem.sect->code);
	//TextOut(hdc,(0.45*(ix2-ix1)+ix1),
	//			(0.45*(iy2-iy1)+iy1),str,strlen(str));
//    TextOut(hdc,(0.3*(ix2-ix1)+ix1),                            //mihara
//                (0.3*(iy2-iy1)+iy1+50),str,strlen(str));
  }
  if(vp.vflag.ev.srcanrate) /*SRCAN RATE.*/
  {
	SetTextColor(hdc,RGB(cred,cgreen,cblue));

    rate=0.0;
    for(i=0;i<4;i++)
    {
      if(rate<elem.srate[i]) rate=elem.srate[i];
	}

    sprintf(str,"%.3f",rate);
    TextOut(hdc,(0.45*(ix2-ix1)+ix1),
                (0.45*(iy2-iy1)+iy1),str,strlen(str));
  }
  if(vp.vflag.ev.evalue) /*ENERGY VALUE*/
  {
    SetTextColor(hdc,RGB(cred,cgreen,cblue));
    sprintf(str,"%.3f",elem.Ee[0]+elem.Ep[0]);
    TextOut(hdc,(0.25*(ix2-ix1)+ix1),
                (0.25*(iy2-iy1)+iy1),str,strlen(str));
	sprintf(str,"%.3f",elem.Ee[1]+elem.Ep[1]);
    TextOut(hdc,(0.75*(ix2-ix1)+ix1),
                (0.75*(iy2-iy1)+iy1),str,strlen(str));
  }

  /*INITIAL ELEMENT LINE*/
/*  if(mode==ONPRINTER)
  {
	hpen=CreatePen(PS_SOLID,2,RGB(ered,egreen,eblue));
  }
  else*/ if(mode==ONSRCANMAX)
  {
	hpen=CreatePen(PS_SOLID,4,RGB(ered,egreen,eblue));
  }
  else if(mode==ONPREVIEW)
  {
    hpen=CreatePen(PS_SOLID,1,RGB(ered,egreen,eblue));
  }
  /*UJIOKA FOR COLOR*/
#if BOLDRED
  else if(vp.vflag.ev.srcancolor && ered==255 && egreen==0 && eblue==150)
  {
	  hpen=CreatePen(PS_SOLID,12,RGB(ered,egreen,eblue));
  }
  else if(vp.vflag.ev.srcancolor && ered==255 && egreen==150 && eblue==50)
  {
	  hpen=CreatePen(PS_SOLID,10,RGB(ered,egreen,eblue));
  }
  else if(vp.vflag.ev.srcancolor && ered==150 && egreen==150 && eblue==0)
  {
	  hpen=CreatePen(PS_SOLID,8,RGB(ered,egreen,eblue));
  }
  else if(vp.vflag.ev.srcancolor && ered==255 && egreen==255 && eblue==0)
  {
	  hpen=CreatePen(PS_SOLID,8,RGB(ered,egreen,eblue));
  }
  else if(vp.vflag.ev.srcancolor && ered==150 && egreen==255 && eblue==0)
  {
	  hpen=CreatePen(PS_SOLID,6,RGB(ered,egreen,eblue));
  }
  else if(vp.vflag.ev.srcancolor && ered==0 && egreen==255 && eblue==0)
  {
	  hpen=CreatePen(PS_SOLID,4,RGB(ered,egreen,eblue));
  }
  else if(vp.vflag.ev.srcancolor && ered==0 && egreen==255 && eblue==150)
  {
	  hpen=CreatePen(PS_SOLID,2,RGB(ered,egreen,eblue));
  }
#endif
  else hpen=CreatePen(PS_SOLID,1,RGB(ered,egreen,eblue));

  ppen=(HPEN)SelectObject(hdc,hpen);
  MoveToEx(hdc,ix1,iy1,NULL);
  LineTo(hdc,ix2,iy2);
  SelectObject(hdc,ppen);
  DeleteObject(hpen);

  /*DEFORMED ELEMENT LINE*/
  if(vp.vflag.ev.deformation && af.ddisp!=NULL)
  {
	if(mode==ONPREVIEW /*|| mode==ONPRINTER*/)
	{
	  ndiv=1; /*DEFORMATION DIVISION*/

	  for(i=ndiv;i>0;i--)
	  {
		hpen=CreatePen(PS_SOLID,1,RGB((int)(100.0*i/ndiv),
									  (int)(100.0*i/ndiv),
									  (int)(  0.0*i/ndiv)));/*REVERSE*/
		/*hpen=CreatePen(PS_SOLID,1,RGB((int)(  0.0*i/ndiv),
									  (int)(150.0*i/ndiv),
									  (int)(255.0*i/ndiv)));*//*BLUE*/
		/*hpen=CreatePen(PS_SOLID,1,RGB((int)(255.0*i/ndiv),
									  (int)(150.0*i/ndiv),
									  (int)(255.0*i/ndiv)));*//*MAGENTA*/
		ppen=(HPEN)SelectObject(hdc,hpen);

		MoveToEx(hdc,(int)((double)ix1
						   +(double)(i*(dx1-ix1))/(double)ndiv),
					 (int)((double)iy1
						   +(double)(i*(dy1-iy1))/(double)ndiv),
					 NULL);
		LineTo(hdc,(int)((double)ix2
						 +(double)(i*(dx2-ix2))/(double)ndiv),
				   (int)((double)iy2
						 +(double)(i*(dy2-iy2))/(double)ndiv));

		SelectObject(hdc,ppen);
		DeleteObject(hpen);
	  }
	}
	else
	{
      if(mode==ONPRINTER)
	  {
//        hpen=CreatePen(PS_SOLID,2,RGB(ered,egreen,eblue));
          if((vp.vflag.ev.stress[1][0] && !globalcondensationflag)||
			 (wdraw.childs+1)->vparam.vflag.mv.pagenum)
		  {
			  hpen=CreatePen(PS_DASH,1,RGB(120,120,120));
          }
		  else hpen=CreatePen(PS_SOLID,1,RGB(0,0,255));  /*変形後の部材の色（印刷用）*/
//          else hpen=CreatePen(PS_DASH,1,RGB(120,120,120));
	  }

      //else hpen=CreatePen(PS_DASH,1,RGB(ered,egreen,eblue));
	  //ujioka for deformation
	  else hpen=CreatePen(PS_SOLID,1,RGB(100,100,255));         /*変形後の部材の色*/

	  ppen=(HPEN)SelectObject(hdc,hpen);

      MoveToEx(hdc,dx1,dy1,NULL);
      LineTo(hdc,dx2,dy2);

      SelectObject(hdc,ppen);
      DeleteObject(hpen);
    }
  }

  /*INITIAL ELEMENT LINE*/
  /*
  if(mode==ONPRINTER)
  {
    hpen=CreatePen(PS_SOLID,2,RGB(ered,egreen,eblue));
  }
  else if(mode==ONPREVIEW)
  {
    hpen=CreatePen(PS_SOLID,1,RGB(ered,egreen,eblue));
  }
  else hpen=CreatePen(PS_SOLID,1,RGB(ered,egreen,eblue));
  ppen=(HPEN)SelectObject(hdc,hpen);
  MoveToEx(hdc,ix1,iy1,NULL);
  LineTo(hdc,ix2,iy2);
  SelectObject(hdc,ppen);
  DeleteObject(hpen);
  */

  if(vp.vflag.ev.axis)                              /*ELEMENT AXIS.*/
  {
	drawwireaxis(hdc,vp,inode1,inode2,elem.cangle);
  }
  if(vp.vflag.ev.cmqcheck) /*FOR CMQ CHECK*/ /*kaza & uji for Lunar/MarsBase*/
  {
	drawelemcmqcheck(hdc,vp,inode1,inode2,elem.cangle,-elem.stress[0][2],-elem.stress[1][2]);
  }
  if(vp.vflag.ev.sectionshape)
  {                                                                      //honda
	drawwireshape(hdc,vp,inode1,inode2,elem.cangle,cred,cgreen,cblue);
  }
  if(!strncmp(prj,"sydney",6))    /*steel weight for sydney*/     //honda
  {
	SetTextColor(hdc,RGB(255,255,255));
	char str[256];
	  //setfontformat(hdc,20,8,"MS Mincho",255,255,255);
	  sprintf(str,"\nsteel weight = %.3f tf",totalweight);
	  TextOut(hdc, 300,50,str,strlen(str));
	  sprintf(str,"\narch length  = %.3f m",totallength);
	  TextOut(hdc, 300,65,str,strlen(str));
	  sprintf(str,"\ncenter of gravity = (%.3f,%.3f,%.3f) ",Gcx,Gcy,Gcz);
	  TextOut(hdc, 300,80,str,strlen(str));
	  //setfontformat(hdc,15,6,"MS Mincho",255,255,255);
  } /*steel weight for sydney*/



  if(vp.vflag.ev.hinge)                            /*INITIAL HINGE.*/
  {
	if(mode==ONPRINTER)
	{
	  hpen=CreatePen(PS_SOLID,2,RGB(0,0,0)); /*BLACK*/
	}
	else if(mode==ONPREVIEW)
    {
	  hpen=CreatePen(PS_SOLID,1,RGB(0,0,0)); /*BLACK*/
	}
    else hpen=CreatePen(PS_SOLID,1,RGB(255,255,255));
	ppen=(HPEN)SelectObject(hdc,hpen);

	if(elem.iconf[0][3]==1 ||
	   elem.iconf[0][4]==1 ||
	   elem.iconf[0][5]==1 /*||
	   elem.sect->Ixx==0.0 ||
	   elem.sect->Iyy==0.0 ||
       elem.sect->Jzz==0.0*/)                      /*HINGE ON HEAD.*/
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
	   elem.iconf[1][5]==1 /*||
	   elem.sect->Ixx==0.0 ||
       elem.sect->Iyy==0.0 ||
	   elem.sect->Jzz==0.0*/)                      /*HINGE ON TAIL.*/
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
  if(vp.vflag.ev.hinge && af.ddisp!=NULL
     && (wdraw.childs+1)->vparam.vflag.ev.deformation) /*INITIAL HINGE DEFORMED.*/
  {
	if(mode==ONPRINTER)
	{
	  hpen=CreatePen(PS_SOLID,2,RGB(0,0,0)); /*BLACK*/
	}
	else if(mode==ONPREVIEW)
	{
	  /*hpen=CreatePen(PS_SOLID,1,RGB(0,0,0));*/ /*BLACK*/
	  hpen=CreatePen(PS_SOLID,1,RGB(100,100,0)); /*REVERSE BLUE*/
	  /*hpen=CreatePen(PS_SOLID,1,RGB(0,150,255));*/ /*BLUE*/
	}
	else hpen=CreatePen(PS_SOLID,1,RGB(255,255,255));
	ppen=(HPEN)SelectObject(hdc,hpen);

	if(elem.iconf[0][3]==1 ||
	   elem.iconf[0][4]==1 ||
	   elem.iconf[0][5]==1 /*||
	   elem.sect->Ixx==0.0 ||
	   elem.sect->Iyy==0.0 ||
	   elem.sect->Jzz==0.0*/)                      /*HINGE ON HEAD.*/
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
       elem.iconf[1][5]==1 /*||
       elem.sect->Ixx==0.0 ||
       elem.sect->Iyy==0.0 ||
       elem.sect->Jzz==0.0*/)                      /*HINGE ON TAIL.*/
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
  if(vp.vflag.ev.hinge && af.ddisp!=NULL
     && (wdraw.childs+1)->vparam.vflag.ev.deformation)          /*PLASTIC HINGE.*/
  {
    if(mode==ONPRINTER)
    {
	  hpen=CreatePen(PS_SOLID,1,RGB(0,0,0)); /*BLACK*/
      hbrush=CreateSolidBrush(RGB(0,0,0));   /*BLACK*/
    }
    else if(mode==ONPREVIEW)
    {
      hpen=CreatePen(PS_SOLID,1,RGB(0,0,0)); /*BLACK*/
      hbrush=CreateSolidBrush(RGB(0,0,0));   /*BLACK*/
    }
    else
    {
      hpen=CreatePen(PS_SOLID,1,RGB(0,255,255));    /*PEN OF HINGE.*/
      hbrush=CreateSolidBrush(RGB(0,255,255));    /*BRUSH OF HINGE.*/
    }
	ppen=(HPEN)SelectObject(hdc,hpen);
    pbrush=(HBRUSH)SelectObject(hdc,hbrush);

    if(elem.iconf[0][0]<0)                  /*PLASTIC HINGE ON HEAD*/
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
    if(elem.iconf[1][0]<0)                  /*PLASTIC HINGE ON TAIL*/
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

  if(vp.vflag.ev.ecircle)                                 /*ENERGY.*/
  {
	getdoublefromdialog((wmenu.childs+4)->hwnd,IDVS_OPAQUE,&orate);
	getmasshiju((wmenu.childs+4)->hwnd,&hiju);

	if(mode==ONPRINTER)
	{
      hpen=CreatePen(PS_SOLID,1,RGB(0,0,0));
      hbrush=CreateSolidBrush(RGB(150,255,255));
    }
    else if(mode==ONPREVIEW)
    {
      hpen=CreatePen(PS_SOLID,1,RGB(0,0,0));
      hbrush=CreateSolidBrush(RGB(150,255,255));
    }
    else
    {
      hpen=CreatePen(PS_SOLID,1,RGB(0,255,255));   /*PEN OF ENERGY.*/
      hbrush=CreateSolidBrush(RGB(0,100,100));   /*BRUSH OF ENERGY.*/
    }
    ppen=(HPEN)SelectObject(hdc,hpen);
    pbrush=(HBRUSH)SelectObject(hdc,hbrush);

    /*HEAD*/
if(!vp.vflag.ev.plasticenergy)    etotal=elem.Ee[0]+elem.Ep[0];
else                              etotal=elem.Ep[0]; /*Ujioka*/
    if(etotal>0.0)
    {
      radius=pow(3.0/4.0/PI*(etotal/hiju),1.0/3.0);     /*VOLUME*/
//      radius=pow(1/PI*(etotal/hiju),1.0/2.0);         /*AREA*/
    }
    else radius=0.0;

    if(mode==ONPRINTER) radius*=10.0;
    if(vp.vflag.ev.deformation && af.ddisp!=NULL)
    {
      dx=dx2-dx1; dy=dy2-dy1;
    }
    else                  /*Ujioka*/
    {
      dx=ix2-ix1; dy=iy2-iy1;
    }
    oblique=sqrt(dx*dx+dy*dy);
    if(oblique==0.0)
    {
      if(vp.vflag.ev.deformation && af.ddisp!=NULL)
      {
        hingeX=dx1;
        hingeY=dy1;
      }
      else                /*Ujioka*/
      {
        hingeX=ix1;
        hingeY=iy1;
      }
    }
    else
    {
	  if(vp.vflag.ev.deformation && af.ddisp!=NULL)
      {
        hingeX=dx1+(int)(radius*dx/oblique);
        hingeY=dy1+(int)(radius*dy/oblique);
      }
      else               /*Ujioka*/
      {
        hingeX=ix1+(int)(radius*dx/oblique);
        hingeY=iy1+(int)(radius*dy/oblique);
      }
    }
    Ellipse(hdc,                                  /*FILLED CIRCLE*/
            (hingeX-radius),(hingeY-radius),
            (hingeX+radius),(hingeY+radius));

    /*TAIL*/
if(!vp.vflag.ev.plasticenergy)        etotal=elem.Ee[1]+elem.Ep[1];
else                                  etotal=elem.Ep[1]; /*Ujioka*/
    if(etotal>0.0)
    {
      radius=pow(3.0/4.0/PI*(etotal/hiju),1.0/3.0);     /*VOLUME*/
//      radius=pow(1/PI*(etotal/hiju),1.0/2.0);         /*AREA*/
    }
    else radius=0.0;

    if(mode==ONPRINTER) radius*=10.0;
	if(vp.vflag.ev.deformation && af.ddisp!=NULL)
    {
      dx=dx1-dx2; dy=dy1-dy2;
    }
    else
    {
      dx=ix1-ix2; dy=iy1-iy2;
    }
    oblique=sqrt(dx*dx+dy*dy);
    if(oblique==0.0)
    {
      if(vp.vflag.ev.deformation && af.ddisp!=NULL)
      {
        hingeX=dx2;
        hingeY=dy2;
      }
      else
      {
        hingeX=ix2;
        hingeY=iy2;
      }
    }
    else
    {
      if(vp.vflag.ev.deformation && af.ddisp!=NULL)
      {
		hingeX=dx2+(int)(radius*dx/oblique);
        hingeY=dy2+(int)(radius*dy/oblique);
      }
      else
      {
        hingeX=ix2+(int)(radius*dx/oblique);
        hingeY=iy2+(int)(radius*dy/oblique);
      }
    }
    Ellipse(hdc,                                  /*FILLED CIRCLE*/
            (hingeX-radius),(hingeY-radius),
            (hingeX+radius),(hingeY+radius));

	SelectObject(hdc,ppen);
	DeleteObject(hpen);
	SelectObject(hdc,pbrush);
	DeleteObject(hbrush);
  }

  return;
}/*drawglobalwire*/

void drawglobalshell(HDC hdc,struct viewparam vp,
					struct arclmframe af,
					struct oshell shell,
					int cred,int cgreen,int cblue, /*CODE*/
					int ered,int egreen,int eblue, /*LINE*/
					long int selectcode,
					int mode)
/*DRAW ELEMENT,PLASTIC HINGE BY AXONOMETRICS,PERSPECTIVE.*/
{
  HDC hdcC;
  HBITMAP hbitC,pbit;
  HPEN hpen,ppen;                                  /*HANDLE OF PEN.*/
  HBRUSH hbrush,pbrush;						/*HANDLE OF BRUSH.*/
  POINT *points;
  SIZE size; /****SRCANMAX****/
  DWORD rop; /*RASTER OPERATION.*/
  char str[10];
  int ix1,iy1,ix2,iy2,ix3,iy3,ix4,iy4;
  int dx1,dy1,dx2,dy2,dx3,dy3,dx4,dy4;    /*PROJECTED COORDINATION.*/
  int maxX,minX,maxY,minY;
  long int loff;
  double dx,dy;
  int hingeX,hingeY;                             /*CENTER OF HINGE.*/
  double oblique,radius,hiju,orate,etotal;          /*OBLIQUE SIDE.*/
  struct onode inode1,inode2,inode3;                       /*INITIAL NODE.*/
  struct onode dnode1,dnode2,dnode3;                      /*DEFORMED NODE.*/
  int i,ndiv;
  double rate;


  /*MAKE NODES*/
  initialnode(af.ninit,af.nnode,(shell.node[0]->code),&inode1);
  if(!insiderange(inode1,vp.range)) return;
  if(af.ddisp==NULL || !vp.vflag.ev.deformation)
  {
	dnode1.d[GX]=0.0;
	dnode1.d[GY]=0.0;
	dnode1.d[GZ]=0.0;
  }
  else                                              /*DISPLACEMENT.*/
  {
	/*inputnode(af.ddisp,shell.node[0]);*/ /*DEFORMED NODE.*/
	dnode1.d[GX]=shell.node[0]->d[GX]-inode1.d[GX];
	dnode1.d[GY]=shell.node[0]->d[GY]-inode1.d[GY];
	dnode1.d[GZ]=shell.node[0]->d[GZ]-inode1.d[GZ];

	dnode1.d[GX]=inode1.d[GX]+(vp.dparam.dfact)*dnode1.d[GX];
	dnode1.d[GY]=inode1.d[GY]+(vp.dparam.dfact)*dnode1.d[GY];
	dnode1.d[GZ]=inode1.d[GZ]+(vp.dparam.dfact)*dnode1.d[GZ];

	nodeontoscreen(dnode1,&dx1,&dy1,vp);              /*PROJECTION.*/
  }
  nodeontoscreen(inode1,&ix1,&iy1,vp);                /*PROJECTION.*/




  initialnode(af.ninit,af.nnode,(shell.node[1]->code),&inode2);
  if(!insiderange(inode2,vp.range)) return;
  if(af.ddisp==NULL || !vp.vflag.ev.deformation)
  {
	dnode2.d[GX]=0.0;
    dnode2.d[GY]=0.0;
    dnode2.d[GZ]=0.0;
  }
  else                                              /*DISPLACEMENT.*/
  {
	/*inputnode(af.ddisp,shell.node[1]);*/ /*DEFORMED NODE.*/
	dnode2.d[GX]=shell.node[1]->d[GX]-inode2.d[GX];
	dnode2.d[GY]=shell.node[1]->d[GY]-inode2.d[GY];
	dnode2.d[GZ]=shell.node[1]->d[GZ]-inode2.d[GZ];

    dnode2.d[GX]=inode2.d[GX]+(vp.dparam.dfact)*dnode2.d[GX];
	dnode2.d[GY]=inode2.d[GY]+(vp.dparam.dfact)*dnode2.d[GY];
	dnode2.d[GZ]=inode2.d[GZ]+(vp.dparam.dfact)*dnode2.d[GZ];

	nodeontoscreen(dnode2,&dx2,&dy2,vp);              /*PROJECTION.*/
  }
  nodeontoscreen(inode2,&ix2,&iy2,vp);                /*PROJECTION.*/



  initialnode(af.ninit,af.nnode,(shell.node[2]->code),&inode3);
  if(!insiderange(inode3,vp.range)) return;
  if(af.ddisp==NULL || !vp.vflag.ev.deformation)
  {
	dnode3.d[GX]=0.0;
	dnode3.d[GY]=0.0;
	dnode3.d[GZ]=0.0;
  }
  else                                              /*DISPLACEMENT.*/
  {
	/*inputnode(af.ddisp,shell.node[2]);*/ /*DEFORMED NODE.*/
	dnode3.d[GX]=shell.node[2]->d[GX]-inode3.d[GX];
	dnode3.d[GY]=shell.node[2]->d[GY]-inode3.d[GY];
	dnode3.d[GZ]=shell.node[2]->d[GZ]-inode3.d[GZ];

	dnode3.d[GX]=inode3.d[GX]+(vp.dparam.dfact)*dnode3.d[GX];
	dnode3.d[GY]=inode3.d[GY]+(vp.dparam.dfact)*dnode3.d[GY];
	dnode3.d[GZ]=inode3.d[GZ]+(vp.dparam.dfact)*dnode3.d[GZ];

	nodeontoscreen(dnode3,&dx3,&dy3,vp);              /*PROJECTION.*/
  }
  nodeontoscreen(inode3,&ix3,&iy3,vp);                /*PROJECTION.*/



  /*NODE CODE PRINT*/
  if(vp.vflag.nv.code)
  {
	if(globalstatus==SELECTNODE && shell.node[0]->code==selectcode)
    {
	  SetTextColor(hdc,RGB(0,255,255));
    }
	else if(mode==ONPRINTER)
	{
	  SetTextColor(hdc,RGB(0,0,0)); /*BLACK*/
	}
	else if(mode==ONPREVIEW)
	{
	  SetTextColor(hdc,RGB(0,0,0)); /*BLACK*/
	}
	else
	{
	  loff=shell.node[0]->loff;
	  if((af.confs+6*loff+3)->iconf==1 ||
		 (af.confs+6*loff+4)->iconf==1 ||
		 (af.confs+6*loff+5)->iconf==1) /*FIXED*/
	  {
		SetTextColor(hdc,RGB(255,150,150));
	  }
      else if((af.confs+6*loff+0)->iconf==1 ||
              (af.confs+6*loff+1)->iconf==1 ||
			  (af.confs+6*loff+2)->iconf==1) /*HINGE*/
	  {
		SetTextColor(hdc,RGB(150,255,150));
      }
	  else SetTextColor(hdc,RGB(150,150,255)); /*FREE*/
	}
	sprintf(str,"%d",shell.node[0]->code);
	TextOut(hdc,ix1+2,iy1,str,strlen(str));




	if(globalstatus==SELECTNODE && shell.node[1]->code==selectcode)
	{
	  SetTextColor(hdc,RGB(0,255,255));
	}
	else if(mode==ONPRINTER)
	{
	  SetTextColor(hdc,RGB(0,0,0)); /*BLACK*/
	}
	else if(mode==ONPREVIEW)
	{
	  SetTextColor(hdc,RGB(0,0,0)); /*BLACK*/
	}
	else
	{
	  loff=shell.node[1]->loff;
	  if((af.confs+6*loff+3)->iconf==1 ||
		 (af.confs+6*loff+4)->iconf==1 ||
		 (af.confs+6*loff+5)->iconf==1) /*FIXED*/
	  {
		SetTextColor(hdc,RGB(255,150,150));
	  }
	  else if((af.confs+6*loff+0)->iconf==1 ||
			  (af.confs+6*loff+1)->iconf==1 ||
			  (af.confs+6*loff+2)->iconf==1) /*HINGE*/
	  {
		SetTextColor(hdc,RGB(150,255,150));
	  }
	  else SetTextColor(hdc,RGB(150,150,255)); /*FREE*/
	}

	sprintf(str,"%d",shell.node[1]->code);
	TextOut(hdc,ix2+2,iy2,str,strlen(str));



	if(globalstatus==SELECTNODE && shell.node[2]->code==selectcode)
	{
	  SetTextColor(hdc,RGB(0,255,255));
	}
	else if(mode==ONPRINTER)
	{
	  SetTextColor(hdc,RGB(0,0,0)); /*BLACK*/
	}
	else if(mode==ONPREVIEW)
	{
	  SetTextColor(hdc,RGB(0,0,0)); /*BLACK*/
	}
	else
	{
	  loff=shell.node[2]->loff;
	  if((af.confs+6*loff+3)->iconf==1 ||
		 (af.confs+6*loff+4)->iconf==1 ||
		 (af.confs+6*loff+5)->iconf==1) /*FIXED*/
	  {
		SetTextColor(hdc,RGB(255,150,150));
	  }
	  else if((af.confs+6*loff+0)->iconf==1 ||
			  (af.confs+6*loff+1)->iconf==1 ||
			  (af.confs+6*loff+2)->iconf==1) /*HINGE*/
	  {
		SetTextColor(hdc,RGB(150,255,150));
	  }
	  else SetTextColor(hdc,RGB(150,150,255)); /*FREE*/
	}

	sprintf(str,"%d",shell.node[2]->code);
	TextOut(hdc,ix3+2,iy3,str,strlen(str));



  }
  /*PRINT ELEMENT CODE*/
  if(vp.vflag.ev.code)
  {
	SetTextColor(hdc,RGB(cred,cgreen,cblue)); /*ELEM CODE.*/
	/****SRCANMAX****/
	if(vp.vflag.ev.srcanmax) sprintf(str,"部材 %d",shell.code);
	else                     sprintf(str,"%d",shell.code);
	if(vp.vflag.ev.srcanmax && mode==ONSRCANMAX)
	{
	  setfontformat(hdc,(int)(1.2*gprn.jiheight),
						(int)(1.2*gprn.jiwidth),"ＭＳ 明朝",0,0,0);
	}
	TextOut(hdc,0.33*(ix1+ix2+ix3),
	            0.33*(iy1+iy2+iy3),str,strlen(str));
  }
  /*PRINT SECTION CODE*/
  if(vp.vflag.ev.sectioncode) /*SECTION CODE.*/
  {

	if(mode==ONPRINTER)
	{
	  SetTextColor(hdc,RGB(0,0,0)); /*BLACK*/
    }
    else if(mode==ONPREVIEW)
    {
	  SetTextColor(hdc,RGB(0,0,0)); /*BLACK*/
	}
    else SetTextColor(hdc,RGB(150,150,150));

	/****SRCANMAX****/
	if(vp.vflag.ev.srcanmax) sprintf(str,"断面 %d",shell.sect->ocode);
	else                     sprintf(str,"%d",shell.sect->code);
	if(vp.vflag.ev.srcanmax && mode==ONSRCANMAX)
	{
	  setfontformat(hdc,(int)(1.2*gprn.jiheight),
						(int)(1.2*gprn.jiwidth),"ＭＳ 明朝",0,0,0);
	}
	if(vp.vflag.ev.code)
	{
	  if(vp.vflag.ev.srcanmax)
	  {
		if(mode==ONSRCANMAX) SetTextColor(hdc,RGB(0,0,0));
		else                 SetTextColor(hdc,RGB(255,0,255));
	  }
	  GetTextExtentPoint32(hdc,str,strlen(str),&size);
	  TextOut(hdc,0.33*(ix1+ix2+ix3),
				  0.33*(iy1+iy2+iy3)-size.cy,str,strlen(str));
	}
	else
	{
	  TextOut(hdc,0.33*(ix1+ix2+ix3),
				  0.33*(iy1+iy2+iy3),str,strlen(str));
	}

  }
#if 0
  if(vp.vflag.ev.srcanrate) /*SRCAN RATE.*/
  {
	SetTextColor(hdc,RGB(cred,cgreen,cblue));

    rate=0.0;
    for(i=0;i<4;i++)
    {
      if(rate<elem.srate[i]) rate=elem.srate[i];
	}

    sprintf(str,"%.3f",rate);
    TextOut(hdc,(0.45*(ix2-ix1)+ix1),
                (0.45*(iy2-iy1)+iy1),str,strlen(str));
  }
  if(vp.vflag.ev.evalue) /*ENERGY VALUE*/
  {
    SetTextColor(hdc,RGB(cred,cgreen,cblue));
    sprintf(str,"%.3f",elem.Ee[0]+elem.Ep[0]);
    TextOut(hdc,(0.25*(ix2-ix1)+ix1),
                (0.25*(iy2-iy1)+iy1),str,strlen(str));
	sprintf(str,"%.3f",elem.Ee[1]+elem.Ep[1]);
    TextOut(hdc,(0.75*(ix2-ix1)+ix1),
				(0.75*(iy2-iy1)+iy1),str,strlen(str));
  }
#endif
  /*INITIAL ELEMENT LINE*/
/*  if(mode==ONPRINTER)
  {
	hpen=CreatePen(PS_SOLID,2,RGB(ered,egreen,eblue));
  }
  else*/ if(mode==ONSRCANMAX)
  {
	hpen=CreatePen(PS_SOLID,4,RGB(ered,egreen,eblue));
  }
  else if(mode==ONPREVIEW)
  {
	hpen=CreatePen(PS_SOLID,1,RGB(ered,egreen,eblue));
  }
  /*UJIOKA FOR COLOR*/
#if 0
  else if(vp.vflag.ev.srcancolor && ered==255 && egreen==0 && eblue==150)
  {
	  hpen=CreatePen(PS_SOLID,12,RGB(ered,egreen,eblue));
  }
  else if(vp.vflag.ev.srcancolor && ered==255 && egreen==150 && eblue==50)
  {
	  hpen=CreatePen(PS_SOLID,10,RGB(ered,egreen,eblue));
  }
  else if(vp.vflag.ev.srcancolor && ered==150 && egreen==150 && eblue==0)
  {
	  hpen=CreatePen(PS_SOLID,8,RGB(ered,egreen,eblue));
  }
  else if(vp.vflag.ev.srcancolor && ered==255 && egreen==255 && eblue==0)
  {
	  hpen=CreatePen(PS_SOLID,8,RGB(ered,egreen,eblue));
  }
  else if(vp.vflag.ev.srcancolor && ered==150 && egreen==255 && eblue==0)
  {
	  hpen=CreatePen(PS_SOLID,6,RGB(ered,egreen,eblue));
  }
  else if(vp.vflag.ev.srcancolor && ered==0 && egreen==255 && eblue==0)
  {
	  hpen=CreatePen(PS_SOLID,4,RGB(ered,egreen,eblue));
  }
  else if(vp.vflag.ev.srcancolor && ered==0 && egreen==255 && eblue==150)
  {
	  hpen=CreatePen(PS_SOLID,2,RGB(ered,egreen,eblue));
  }
#endif
  else hpen=CreatePen(PS_SOLID,1,RGB(ered,egreen,eblue));
  //else hpen=CreatePen(PS_SOLID,1,RGB(255,255,255));


  ppen=(HPEN)SelectObject(hdc,hpen);
  MoveToEx(hdc,ix1,iy1,NULL);
  LineTo(hdc,ix2,iy2);
  MoveToEx(hdc,ix2,iy2,NULL);
  LineTo(hdc,ix3,iy3);
  MoveToEx(hdc,ix3,iy3,NULL);
  LineTo(hdc,ix1,iy1);
  SelectObject(hdc,ppen);
  DeleteObject(hpen);


  points=(POINT *)malloc(shell.nnod*sizeof(POINT));
  if(points==NULL) return;
  (points+0)->x=ix1;
  (points+0)->y=iy1;
  (points+1)->x=ix2;
  (points+1)->y=iy2;
  (points+2)->x=ix3;
  (points+2)->y=iy3;
  if(shell.nnod==4)
  {
	(points+3)->x=ix4;
	(points+3)->y=iy4;
  }


  for(i=0;i<shell.nnod;i++) /*PROJECTION*/
  {
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

  if(mode==ONPRINTER)      rop=SRCAND;
  else if(mode==ONPREVIEW) rop=SRCAND;
  else                     rop=SRCPAINT;

  hdcC = CreateCompatibleDC(hdc);


  hbitC = CreateCompatibleBitmap(hdc,maxX-minX,maxY-minY);
  pbit=(HBITMAP)SelectObject(hdcC,hbitC);

  hbrush = (HBRUSH)GetCurrentObject(hdc,OBJ_BRUSH);
  pbrush = (HBRUSH)SelectObject(hdcC,hbrush);

  PatBlt(hdcC,0,0,(maxX-minX),(maxY-minY),PATCOPY);

  hbrush = (HBRUSH)CreateSolidBrush(RGB(0.1*ered,0.1*egreen,0.1*eblue));
  SelectObject(hdcC,hbrush);
  hpen=CreatePen(PS_SOLID,1,RGB(0.1*ered,0.1*egreen,0.1*eblue));
  ppen=(HPEN)SelectObject(hdcC,hpen);

  for(i=0;i<shell.nnod;i++)
  {
	(points+i)->x-=minX;
	(points+i)->y-=minY;
  }

  Polygon(hdcC,points,shell.nnod);

  BitBlt(hdc,minX,minY,(maxX-minX),(maxY-minY),
		 hdcC,0,0,
		 rop); /*GROUND BLACK:SRCPAINT WHITE:SRCAND.*/

  SelectObject(hdcC,ppen);
  SelectObject(hdcC,pbrush);
  SelectObject(hdcC,pbit);
  DeleteObject(hpen);
  DeleteObject(hbrush);
  DeleteObject(hbitC);
  DeleteDC(hdcC);
  free(points);




  /*DEFORMED ELEMENT LINE*/
  if(vp.vflag.ev.deformation && af.ddisp!=NULL)
  {
	if(mode==ONPREVIEW /*|| mode==ONPRINTER*/)
	{
	  ndiv=1; /*DEFORMATION DIVISION*/

	  for(i=ndiv;i>0;i--)
	  {
		hpen=CreatePen(PS_SOLID,1,RGB((int)(100.0*i/ndiv),
									  (int)(100.0*i/ndiv),
									  (int)(  0.0*i/ndiv)));/*REVERSE*/
		/*hpen=CreatePen(PS_SOLID,1,RGB((int)(  0.0*i/ndiv),
									  (int)(150.0*i/ndiv),
									  (int)(255.0*i/ndiv)));*//*BLUE*/
		/*hpen=CreatePen(PS_SOLID,1,RGB((int)(255.0*i/ndiv),
									  (int)(150.0*i/ndiv),
									  (int)(255.0*i/ndiv)));*//*MAGENTA*/
		ppen=(HPEN)SelectObject(hdc,hpen);

		MoveToEx(hdc,(int)((double)ix1
						   +(double)(i*(dx1-ix1))/(double)ndiv),
					 (int)((double)iy1
						   +(double)(i*(dy1-iy1))/(double)ndiv),
					 NULL);
		LineTo(hdc,(int)((double)ix2
						 +(double)(i*(dx2-ix2))/(double)ndiv),
				   (int)((double)iy2
						 +(double)(i*(dy2-iy2))/(double)ndiv));
		MoveToEx(hdc,(int)((double)ix2
						   +(double)(i*(dx2-ix2))/(double)ndiv),
					 (int)((double)iy2
						   +(double)(i*(dy2-iy2))/(double)ndiv),
					 NULL);
		LineTo(hdc,(int)((double)ix3
						 +(double)(i*(dx3-ix3))/(double)ndiv),
				   (int)((double)iy3
						 +(double)(i*(dy3-iy3))/(double)ndiv));
		MoveToEx(hdc,(int)((double)ix3
						   +(double)(i*(dx3-ix3))/(double)ndiv),
					 (int)((double)iy3
						   +(double)(i*(dy3-iy3))/(double)ndiv),
					 NULL);
		LineTo(hdc,(int)((double)ix1
						 +(double)(i*(dx1-ix1))/(double)ndiv),
				   (int)((double)iy1
						 +(double)(i*(dy1-iy1))/(double)ndiv));

		SelectObject(hdc,ppen);
		DeleteObject(hpen);
	  }
	}
	else
	{
      if(mode==ONPRINTER)
	  {
//        hpen=CreatePen(PS_SOLID,2,RGB(ered,egreen,eblue));
          if((vp.vflag.ev.stress[1][0] && !globalcondensationflag)||
			 (wdraw.childs+1)->vparam.vflag.mv.pagenum)
		  {
			  hpen=CreatePen(PS_DASH,1,RGB(120,120,120));
          }
		  else hpen=CreatePen(PS_SOLID,1,RGB(0,0,255));  /*変形後の部材の色（印刷用）*/
//        else hpen=CreatePen(PS_DASH,1,RGB(120,120,120));
	  }

	  //else hpen=CreatePen(PS_DASH,1,RGB(ered,egreen,eblue));
	  //else hpen=CreatePen(PS_SOLID,1,RGB(100,100,255));         /*変形後の部材の色*/
	  else hpen=CreatePen(PS_SOLID,1,RGB(255,255,255));         /*変形後の部材の色*/

	  ppen=(HPEN)SelectObject(hdc,hpen);

      MoveToEx(hdc,dx1,dy1,NULL);
	  LineTo(hdc,dx2,dy2);
	  MoveToEx(hdc,dx2,dy2,NULL);
	  LineTo(hdc,dx3,dy3);
	  MoveToEx(hdc,dx3,dy3,NULL);
	  LineTo(hdc,dx1,dy1);

      SelectObject(hdc,ppen);
      DeleteObject(hpen);
    }
  }

  /*INITIAL ELEMENT LINE*/
  /*
  if(mode==ONPRINTER)
  {
    hpen=CreatePen(PS_SOLID,2,RGB(ered,egreen,eblue));
  }
  else if(mode==ONPREVIEW)
  {
    hpen=CreatePen(PS_SOLID,1,RGB(ered,egreen,eblue));
  }
  else hpen=CreatePen(PS_SOLID,1,RGB(ered,egreen,eblue));
  ppen=(HPEN)SelectObject(hdc,hpen);
  MoveToEx(hdc,ix1,iy1,NULL);
  LineTo(hdc,ix2,iy2);
  SelectObject(hdc,ppen);
  DeleteObject(hpen);
  */






#if 0
  if(vp.vflag.ev.axis)                              /*ELEMENT AXIS.*/
  {
	drawwireaxis(hdc,vp,inode1,inode2,elem.cangle);
  }
  if(vp.vflag.ev.cmqcheck) /*FOR CMQ CHECK*/ /*kaza & uji for Lunar/MarsBase*/
  {
	drawelemcmqcheck(hdc,vp,inode1,inode2,elem.cangle,-elem.stress[0][2],-elem.stress[1][2]);
  }
  if(vp.vflag.ev.sectionshape)
  {                                                                      //honda
	drawwireshape(hdc,vp,inode1,inode2,elem.cangle,cred,cgreen,cblue);
  }




  if(vp.vflag.ev.hinge)                            /*INITIAL HINGE.*/
  {
	if(mode==ONPRINTER)
	{
	  hpen=CreatePen(PS_SOLID,2,RGB(0,0,0)); /*BLACK*/
	}
	else if(mode==ONPREVIEW)
    {
	  hpen=CreatePen(PS_SOLID,1,RGB(0,0,0)); /*BLACK*/
	}
    else hpen=CreatePen(PS_SOLID,1,RGB(255,255,255));
	ppen=(HPEN)SelectObject(hdc,hpen);

	if(elem.iconf[0][3]==1 ||
	   elem.iconf[0][4]==1 ||
	   elem.iconf[0][5]==1 /*||
	   elem.sect->Ixx==0.0 ||
	   elem.sect->Iyy==0.0 ||
       elem.sect->Jzz==0.0*/)                      /*HINGE ON HEAD.*/
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
	   elem.iconf[1][5]==1 /*||
	   elem.sect->Ixx==0.0 ||
       elem.sect->Iyy==0.0 ||
	   elem.sect->Jzz==0.0*/)                      /*HINGE ON TAIL.*/
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

  if(vp.vflag.ev.hinge && af.ddisp!=NULL
	 && (wdraw.childs+1)->vparam.vflag.ev.deformation)          /*PLASTIC HINGE.*/
  {
	if(mode==ONPRINTER)
	{
	  hpen=CreatePen(PS_SOLID,1,RGB(0,0,0)); /*BLACK*/
	  hbrush=CreateSolidBrush(RGB(0,0,0));   /*BLACK*/
	}
	else if(mode==ONPREVIEW)
	{
	  hpen=CreatePen(PS_SOLID,1,RGB(0,0,0)); /*BLACK*/
	  hbrush=CreateSolidBrush(RGB(0,0,0));   /*BLACK*/
    }
    else
    {
      hpen=CreatePen(PS_SOLID,1,RGB(0,255,255));    /*PEN OF HINGE.*/
      hbrush=CreateSolidBrush(RGB(0,255,255));    /*BRUSH OF HINGE.*/
    }
	ppen=(HPEN)SelectObject(hdc,hpen);
    pbrush=(HBRUSH)SelectObject(hdc,hbrush);

    if(elem.iconf[0][0]<0)                  /*PLASTIC HINGE ON HEAD*/
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
    if(elem.iconf[1][0]<0)                  /*PLASTIC HINGE ON TAIL*/
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

  if(vp.vflag.ev.ecircle)                                 /*ENERGY.*/
  {
	getdoublefromdialog((wmenu.childs+4)->hwnd,IDVS_OPAQUE,&orate);
	getmasshiju((wmenu.childs+4)->hwnd,&hiju);

	if(mode==ONPRINTER)
	{
      hpen=CreatePen(PS_SOLID,1,RGB(0,0,0));
      hbrush=CreateSolidBrush(RGB(150,255,255));
    }
    else if(mode==ONPREVIEW)
    {
      hpen=CreatePen(PS_SOLID,1,RGB(0,0,0));
      hbrush=CreateSolidBrush(RGB(150,255,255));
    }
    else
    {
      hpen=CreatePen(PS_SOLID,1,RGB(0,255,255));   /*PEN OF ENERGY.*/
      hbrush=CreateSolidBrush(RGB(0,100,100));   /*BRUSH OF ENERGY.*/
    }
    ppen=(HPEN)SelectObject(hdc,hpen);
    pbrush=(HBRUSH)SelectObject(hdc,hbrush);

    /*HEAD*/
if(!vp.vflag.ev.plasticenergy)    etotal=elem.Ee[0]+elem.Ep[0];
else                              etotal=elem.Ep[0]; /*Ujioka*/
    if(etotal>0.0)
    {
      radius=pow(3.0/4.0/PI*(etotal/hiju),1.0/3.0);     /*VOLUME*/
//      radius=pow(1/PI*(etotal/hiju),1.0/2.0);         /*AREA*/
    }
    else radius=0.0;

    if(mode==ONPRINTER) radius*=10.0;
    if(vp.vflag.ev.deformation && af.ddisp!=NULL)
    {
      dx=dx2-dx1; dy=dy2-dy1;
    }
    else                  /*Ujioka*/
    {
      dx=ix2-ix1; dy=iy2-iy1;
    }
    oblique=sqrt(dx*dx+dy*dy);
    if(oblique==0.0)
    {
      if(vp.vflag.ev.deformation && af.ddisp!=NULL)
      {
        hingeX=dx1;
        hingeY=dy1;
      }
      else                /*Ujioka*/
      {
        hingeX=ix1;
        hingeY=iy1;
      }
    }
    else
    {
	  if(vp.vflag.ev.deformation && af.ddisp!=NULL)
      {
        hingeX=dx1+(int)(radius*dx/oblique);
        hingeY=dy1+(int)(radius*dy/oblique);
      }
      else               /*Ujioka*/
      {
        hingeX=ix1+(int)(radius*dx/oblique);
        hingeY=iy1+(int)(radius*dy/oblique);
      }
    }
    Ellipse(hdc,                                  /*FILLED CIRCLE*/
            (hingeX-radius),(hingeY-radius),
            (hingeX+radius),(hingeY+radius));

    /*TAIL*/
if(!vp.vflag.ev.plasticenergy)        etotal=elem.Ee[1]+elem.Ep[1];
else                                  etotal=elem.Ep[1]; /*Ujioka*/
	if(etotal>0.0)
	{
	  radius=pow(3.0/4.0/PI*(etotal/hiju),1.0/3.0);     /*VOLUME*/
//      radius=pow(1/PI*(etotal/hiju),1.0/2.0);         /*AREA*/
	}
	else radius=0.0;

    if(mode==ONPRINTER) radius*=10.0;
	if(vp.vflag.ev.deformation && af.ddisp!=NULL)
    {
      dx=dx1-dx2; dy=dy1-dy2;
    }
    else
    {
      dx=ix1-ix2; dy=iy1-iy2;
    }
    oblique=sqrt(dx*dx+dy*dy);
    if(oblique==0.0)
    {
      if(vp.vflag.ev.deformation && af.ddisp!=NULL)
      {
        hingeX=dx2;
        hingeY=dy2;
      }
      else
      {
        hingeX=ix2;
        hingeY=iy2;
      }
    }
    else
    {
      if(vp.vflag.ev.deformation && af.ddisp!=NULL)
      {
		hingeX=dx2+(int)(radius*dx/oblique);
        hingeY=dy2+(int)(radius*dy/oblique);
      }
      else
      {
        hingeX=ix2+(int)(radius*dx/oblique);
        hingeY=iy2+(int)(radius*dy/oblique);
      }
    }
    Ellipse(hdc,                                  /*FILLED CIRCLE*/
            (hingeX-radius),(hingeY-radius),
            (hingeX+radius),(hingeY+radius));

	SelectObject(hdc,ppen);
	DeleteObject(hpen);
	SelectObject(hdc,pbrush);
	DeleteObject(hbrush);
  }
#endif
  return;
}/*drawglobalshell*/

void drawglobalwireII(HDC hdc,struct viewparam vp,
					struct arclmframe af,
					struct owire elem,
					int cred,int cgreen,int cblue, /*CODE*/
					int ered,int egreen,int eblue, /*LINE*/
					long int selectcode,
					int mode)
/*DRAW ELEMENT,PLASTIC HINGE BY AXONOMETRICS,PERSPECTIVE.*/
/*UJIOKA FOR BUCKLING CONDENSATION(BCLNG003)*/
{
  HPEN hpen,ppen;                                  /*HANDLE OF PEN.*/
  HBRUSH hbrush,pbrush;                          /*HANDLE OF BRUSH.*/
  char str[10];
  int ix1,iy1,ix2,iy2,dx1,dy1,dx2,dy2;    /*PROJECTED COORDINATION.*/
  long int loff;
  double dx,dy;
  int hingeX,hingeY;                             /*CENTER OF HINGE.*/
  double oblique,radius,hiju,orate,etotal;          /*OBLIQUE SIDE.*/
  struct onode inode1,inode2;                       /*INITIAL NODE.*/
  struct onode dnode1,dnode2;                      /*DEFORMED NODE.*/
  int i,ndiv;
  double rate;

  initialnode(af.ninit,af.nnode,(elem.node[0]->code),&inode1);
  if(!insiderange(inode1,vp.range)) return;
  if(af.ddisp==NULL/* || !vp.vflag.ev.deformation*/)
  {
	dnode1.d[GX]=0.0;
	dnode1.d[GY]=0.0;
	dnode1.d[GZ]=0.0;
  }
  else                                              /*DISPLACEMENT.*/
  {
	/*inputnode(af.ddisp,elem.node[0]);*/ /*DEFORMED NODE.*/
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
  if(af.ddisp==NULL/* || !vp.vflag.ev.deformation*/)
  {
    dnode2.d[GX]=0.0;
    dnode2.d[GY]=0.0;
    dnode2.d[GZ]=0.0;
  }
  else                                              /*DISPLACEMENT.*/
  {
    /*inputnode(af.ddisp,elem.node[0]);*/ /*DEFORMED NODE.*/
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
    else if(mode==ONPRINTER)
    {
      SetTextColor(hdc,RGB(0,0,0)); /*BLACK*/
    }
    else if(mode==ONPREVIEW)
    {
      SetTextColor(hdc,RGB(0,0,0)); /*BLACK*/
    }
    else
    {
      loff=elem.node[0]->loff;
      if((af.confs+6*loff+3)->iconf==1 ||
         (af.confs+6*loff+4)->iconf==1 ||
         (af.confs+6*loff+5)->iconf==1) /*FIXED*/
      {
        SetTextColor(hdc,RGB(255,150,150));
      }
      else if((af.confs+6*loff+0)->iconf==1 ||
              (af.confs+6*loff+1)->iconf==1 ||
              (af.confs+6*loff+2)->iconf==1) /*HINGE*/
      {
        SetTextColor(hdc,RGB(150,255,150));
      }
      else SetTextColor(hdc,RGB(150,150,255)); /*FREE*/
    }
    sprintf(str,"%d",elem.node[0]->code);
    TextOut(hdc,ix1+2,iy1,str,strlen(str));

    if(globalstatus==SELECTNODE &&
       elem.node[1]->code==selectcode)
    {
      SetTextColor(hdc,RGB(0,255,255));
    }
    else if(mode==ONPRINTER)
    {
      SetTextColor(hdc,RGB(0,0,0)); /*BLACK*/
    }
    else if(mode==ONPREVIEW)
    {
      SetTextColor(hdc,RGB(0,0,0)); /*BLACK*/
    }
    else
    {
      loff=elem.node[1]->loff;
      if((af.confs+6*loff+3)->iconf==1 ||
         (af.confs+6*loff+4)->iconf==1 ||
         (af.confs+6*loff+5)->iconf==1) /*FIXED*/
      {
        SetTextColor(hdc,RGB(255,150,150));
      }
      else if((af.confs+6*loff+0)->iconf==1 ||
              (af.confs+6*loff+1)->iconf==1 ||
              (af.confs+6*loff+2)->iconf==1) /*HINGE*/
      {
        SetTextColor(hdc,RGB(150,255,150));
      }
      else SetTextColor(hdc,RGB(150,150,255)); /*FREE*/
    }

    sprintf(str,"%d",elem.node[1]->code);
    TextOut(hdc,ix2+2,iy2,str,strlen(str));
  }
  if(vp.vflag.ev.code)
  {
    SetTextColor(hdc,RGB(cred,cgreen,cblue)); /*ELEM CODE.*/
    sprintf(str,"%d",elem.code);
    TextOut(hdc,(0.55*(ix2-ix1)+ix1),
                (0.55*(iy2-iy1)+iy1),str,strlen(str));
  }
  if(vp.vflag.ev.sectioncode) /*SECTION CODE.*/
  {
    /*SetTextColor(hdc,RGB(elem.sect->dcolor.r,
                         elem.sect->dcolor.g,
                         elem.sect->dcolor.b));*/
    if(mode==ONPRINTER)
    {
      SetTextColor(hdc,RGB(0,0,0)); /*BLACK*/
    }
    else if(mode==ONPREVIEW)
    {
      SetTextColor(hdc,RGB(0,0,0)); /*BLACK*/
    }
    else SetTextColor(hdc,RGB(150,150,150));

    sprintf(str,"%d",elem.sect->code);
    TextOut(hdc,(0.45*(ix2-ix1)+ix1),
                (0.45*(iy2-iy1)+iy1),str,strlen(str));
  }
  if(vp.vflag.ev.srcanrate) /*SRCAN RATE.*/
  {
    SetTextColor(hdc,RGB(cred,cgreen,cblue));

    rate=0.0;
    for(i=0;i<4;i++)
    {
      if(rate<elem.srate[i]) rate=elem.srate[i];
    }

    sprintf(str,"%.3f",rate);
    TextOut(hdc,(0.45*(ix2-ix1)+ix1),
                (0.45*(iy2-iy1)+iy1),str,strlen(str));
  }
  if(vp.vflag.ev.evalue) /*ENERGY VALUE*/
  {
    SetTextColor(hdc,RGB(cred,cgreen,cblue));
    sprintf(str,"%.3f",elem.Ee[0]+elem.Ep[0]);
    TextOut(hdc,(0.25*(ix2-ix1)+ix1),
                (0.25*(iy2-iy1)+iy1),str,strlen(str));
    sprintf(str,"%.3f",elem.Ee[1]+elem.Ep[1]);
    TextOut(hdc,(0.75*(ix2-ix1)+ix1),
                (0.75*(iy2-iy1)+iy1),str,strlen(str));
  }

  /*INITIAL ELEMENT LINE*/
  if(mode==ONPRINTER)
  {
	hpen=CreatePen(PS_DOT,2,RGB(ered,egreen,eblue));
  }
  else if(mode==ONSRCANMAX)
  {
    hpen=CreatePen(PS_DOT,4,RGB(ered,egreen,eblue));
  }
  else if(mode==ONPREVIEW)
  {
    hpen=CreatePen(PS_DOT,1,RGB(ered,egreen,eblue));
  }
  /*UJIOKA FOR COLOR*/
#if BOLDRED
  else if(vp.vflag.ev.srcancolor && ered==255 && egreen==0 && eblue==150)
  {
	  hpen=CreatePen(PS_SOLID,10,RGB(ered,egreen,eblue));
  }
#endif
  else hpen=CreatePen(PS_DOT,1,RGB(ered,egreen,eblue));

  ppen=(HPEN)SelectObject(hdc,hpen);
  MoveToEx(hdc,ix1,iy1,NULL);
  LineTo(hdc,ix2,iy2);
  SelectObject(hdc,ppen);
  DeleteObject(hpen);

  /*DEFORMED ELEMENT LINE*/
//  if(vp.vflag.ev.deformation && af.ddisp!=NULL)
  if(0)
  {
    if(mode==ONPREVIEW /*|| mode==ONPRINTER*/)
    {
      ndiv=1; /*DEFORMATION DIVISION*/

      for(i=ndiv;i>0;i--)
      {
        hpen=CreatePen(PS_SOLID,1,RGB((int)(100.0*i/ndiv),
                                      (int)(100.0*i/ndiv),
                                      (int)(  0.0*i/ndiv)));/*REVERSE*/
        /*hpen=CreatePen(PS_SOLID,1,RGB((int)(  0.0*i/ndiv),
                                      (int)(150.0*i/ndiv),
                                      (int)(255.0*i/ndiv)));*//*BLUE*/
        /*hpen=CreatePen(PS_SOLID,1,RGB((int)(255.0*i/ndiv),
                                      (int)(150.0*i/ndiv),
                                      (int)(255.0*i/ndiv)));*//*MAGENTA*/
        ppen=(HPEN)SelectObject(hdc,hpen);

        MoveToEx(hdc,(int)((double)ix1
                           +(double)(i*(dx1-ix1))/(double)ndiv),
                     (int)((double)iy1
                           +(double)(i*(dy1-iy1))/(double)ndiv),
                     NULL);
        LineTo(hdc,(int)((double)ix2
                         +(double)(i*(dx2-ix2))/(double)ndiv),
                   (int)((double)iy2
                         +(double)(i*(dy2-iy2))/(double)ndiv));

        SelectObject(hdc,ppen);
        DeleteObject(hpen);
      }
    }
    else
    {
      if(mode==ONPRINTER)
      {
//        hpen=CreatePen(PS_SOLID,2,RGB(ered,egreen,eblue));
          hpen=CreatePen(PS_SOLID,2,RGB(0,0,255));  /*変形後の部材の色（印刷用）*/
      }
      //else hpen=CreatePen(PS_DASH,1,RGB(ered,egreen,eblue));
      //ujioka for deformation
      else hpen=CreatePen(PS_SOLID,1,RGB(100,100,255));         /*変形後の部材の色*/

      ppen=(HPEN)SelectObject(hdc,hpen);

      MoveToEx(hdc,dx1,dy1,NULL);
      LineTo(hdc,dx2,dy2);

      SelectObject(hdc,ppen);
      DeleteObject(hpen);
    }
  }

  /*INITIAL ELEMENT LINE*/
  /*
  if(mode==ONPRINTER)
  {
    hpen=CreatePen(PS_SOLID,2,RGB(ered,egreen,eblue));
  }
  else if(mode==ONPREVIEW)
  {
    hpen=CreatePen(PS_SOLID,1,RGB(ered,egreen,eblue));
  }
  else hpen=CreatePen(PS_SOLID,1,RGB(ered,egreen,eblue));
  ppen=(HPEN)SelectObject(hdc,hpen);
  MoveToEx(hdc,ix1,iy1,NULL);
  LineTo(hdc,ix2,iy2);
  SelectObject(hdc,ppen);
  DeleteObject(hpen);
  */

  if(vp.vflag.ev.axis)                              /*ELEMENT AXIS.*/
  {
	drawwireaxis(hdc,vp,inode1,inode2,elem.cangle);
  }
  if(vp.vflag.ev.cmqcheck) /*FOR CMQ CHECK*/ /*kaza & uji for Lunar/MarsBase*/
  {
    drawelemcmqcheck(hdc,vp,inode1,inode2,elem.cangle,-elem.stress[0][2],-elem.stress[1][2]);
  }
  if(vp.vflag.ev.sectionshape)
  {                                                                      //honda
	drawwireshape(hdc,vp,inode1,inode2,elem.cangle,cred,cgreen,cblue);
  }
  if(!strncmp(prj,"sydney",6))    /*steel weight for sydney*/     //honda
  {
	SetTextColor(hdc,RGB(255,255,255));
	char str[256];
	  //setfontformat(hdc,20,8,"MS Mincho",255,255,255);
	  sprintf(str,"\nsteel weight = %.3f tf",totalweight);
	  TextOut(hdc, 300,50,str,strlen(str));
	  sprintf(str,"\narch length  = %.3f m",totallength);
	  TextOut(hdc, 300,65,str,strlen(str));
	  sprintf(str,"\ncenter of gravity = (%.3f,%.3f,%.3f) ",Gcx,Gcy,Gcz);
	  TextOut(hdc, 300,80,str,strlen(str));
	  //setfontformat(hdc,15,6,"MS Mincho",255,255,255);
  } /*steel weight for sydney*/



  if(vp.vflag.ev.hinge)                            /*INITIAL HINGE.*/
  {
    if(mode==ONPRINTER)
    {
      hpen=CreatePen(PS_SOLID,2,RGB(0,0,0)); /*BLACK*/
    }
    else if(mode==ONPREVIEW)
    {
      hpen=CreatePen(PS_SOLID,1,RGB(0,0,0)); /*BLACK*/
    }
    else hpen=CreatePen(PS_SOLID,1,RGB(255,255,255));
    ppen=(HPEN)SelectObject(hdc,hpen);

    if(elem.iconf[0][3]==1 ||
       elem.iconf[0][4]==1 ||
       elem.iconf[0][5]==1 /*||
       elem.sect->Ixx==0.0 ||
       elem.sect->Iyy==0.0 ||
       elem.sect->Jzz==0.0*/)                      /*HINGE ON HEAD.*/
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
       elem.iconf[1][5]==1 /*||
       elem.sect->Ixx==0.0 ||
       elem.sect->Iyy==0.0 ||
       elem.sect->Jzz==0.0*/)                      /*HINGE ON TAIL.*/
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
  if(vp.vflag.ev.hinge && af.ddisp!=NULL) /*INITIAL HINGE DEFORMED.*/
  {
    if(mode==ONPRINTER)
    {
      hpen=CreatePen(PS_SOLID,2,RGB(0,0,0)); /*BLACK*/
    }
    else if(mode==ONPREVIEW)
    {
      /*hpen=CreatePen(PS_SOLID,1,RGB(0,0,0));*/ /*BLACK*/
      hpen=CreatePen(PS_SOLID,1,RGB(100,100,0)); /*REVERSE BLUE*/
      /*hpen=CreatePen(PS_SOLID,1,RGB(0,150,255));*/ /*BLUE*/
    }
    else hpen=CreatePen(PS_SOLID,1,RGB(255,255,255));
    ppen=(HPEN)SelectObject(hdc,hpen);

    if(elem.iconf[0][3]==1 ||
       elem.iconf[0][4]==1 ||
       elem.iconf[0][5]==1 /*||
       elem.sect->Ixx==0.0 ||
       elem.sect->Iyy==0.0 ||
       elem.sect->Jzz==0.0*/)                      /*HINGE ON HEAD.*/
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
       elem.iconf[1][5]==1 /*||
       elem.sect->Ixx==0.0 ||
       elem.sect->Iyy==0.0 ||
       elem.sect->Jzz==0.0*/)                      /*HINGE ON TAIL.*/
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
  if(vp.vflag.ev.hinge && af.ddisp!=NULL)          /*PLASTIC HINGE.*/
  {
    if(mode==ONPRINTER)
    {
      hpen=CreatePen(PS_SOLID,1,RGB(0,0,0)); /*BLACK*/
      hbrush=CreateSolidBrush(RGB(0,0,0));   /*BLACK*/
    }
    else if(mode==ONPREVIEW)
    {
      hpen=CreatePen(PS_SOLID,1,RGB(0,0,0)); /*BLACK*/
      hbrush=CreateSolidBrush(RGB(0,0,0));   /*BLACK*/
    }
    else
    {
      hpen=CreatePen(PS_SOLID,1,RGB(0,255,255));    /*PEN OF HINGE.*/
      hbrush=CreateSolidBrush(RGB(0,255,255));    /*BRUSH OF HINGE.*/
    }
    ppen=(HPEN)SelectObject(hdc,hpen);
    pbrush=(HBRUSH)SelectObject(hdc,hbrush);

    if(elem.iconf[0][0]<0)                  /*PLASTIC HINGE ON HEAD*/
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
    if(elem.iconf[1][0]<0)                  /*PLASTIC HINGE ON TAIL*/
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

  if(vp.vflag.ev.ecircle)                                 /*ENERGY.*/
  {
    getdoublefromdialog((wmenu.childs+4)->hwnd,IDVS_OPAQUE,&orate);
    getmasshiju((wmenu.childs+4)->hwnd,&hiju);

    if(mode==ONPRINTER)
    {
	  hpen=CreatePen(PS_SOLID,1,RGB(0,0,0));
      hbrush=CreateSolidBrush(RGB(150,255,255));
    }
    else if(mode==ONPREVIEW)
    {
      hpen=CreatePen(PS_SOLID,1,RGB(0,0,0));
      hbrush=CreateSolidBrush(RGB(150,255,255));
    }
    else
    {
      hpen=CreatePen(PS_SOLID,1,RGB(0,255,255));   /*PEN OF ENERGY.*/
      hbrush=CreateSolidBrush(RGB(0,100,100));   /*BRUSH OF ENERGY.*/
    }
    ppen=(HPEN)SelectObject(hdc,hpen);
    pbrush=(HBRUSH)SelectObject(hdc,hbrush);

    /*HEAD*/
    etotal=elem.Ee[0]+elem.Ep[0];
    if(etotal>0.0)
    {
      radius=pow(3.0/4.0/PI*(etotal/hiju),1.0/3.0);
    }
    else radius=0.0;

    if(mode==ONPRINTER) radius*=10.0;
    dx=dx2-dx1; dy=dy2-dy1;
    oblique=sqrt(dx*dx+dy*dy);
    if(oblique==0.0)
    {
      hingeX=dx1;
      hingeY=dy1;
    }
    else
    {
      hingeX=dx1+(int)(radius*dx/oblique);
      hingeY=dy1+(int)(radius*dy/oblique);
    }
    Ellipse(hdc,                                  /*FILLED CIRCLE*/
            (hingeX-radius),(hingeY-radius),
            (hingeX+radius),(hingeY+radius));

    /*TAIL*/
    etotal=elem.Ee[1]+elem.Ep[1];
    if(etotal>0.0)
    {
      radius=pow(3.0/4.0/PI*(etotal/hiju),1.0/3.0);
    }
    else radius=0.0;

    if(mode==ONPRINTER) radius*=10.0;
    dx=dx1-dx2; dy=dy1-dy2;
    oblique=sqrt(dx*dx+dy*dy);
    if(oblique==0.0)
    {
      hingeX=dx2;
      hingeY=dy2;
    }
    else
    {
      hingeX=dx2+(int)(radius*dx/oblique);
      hingeY=dy2+(int)(radius*dy/oblique);
    }
    Ellipse(hdc,                                  /*FILLED CIRCLE*/
            (hingeX-radius),(hingeY-radius),
            (hingeX+radius),(hingeY+radius));

    SelectObject(hdc,ppen);
    DeleteObject(hpen);
    SelectObject(hdc,pbrush);
    DeleteObject(hbrush);
  }

  return;
}/*drawglobalwire*/
void drawwirestress(HDC hdc,struct viewparam vp,
                    struct arclmframe af,
                    struct owire elem,int mode)
/*DRAW STRESS OF ONE ELEMENT BY AXONOMETRICS,PERSPECTIVE.*/
{

  HPEN hpenM,hpenQ,ppen;                            /*HANDLE OF PEN*/
  HFONT hfont,pfont;
  LOGFONT lfont;
  UINT align;
  char str[20];
  char slabel[6][10]={"Nz","Qx","Qy","Mz","Mx","My"};     /*LABELS.*/
  int i,j;
  int qred,qgreen,qblue; /*SHEAR COLOR*/
  int mred,mgreen,mblue; /*MOMENT COLOR*/
  double value,Qi,Qj,Mi,Mj;
  double elength,angle;
  double **drccos,**tdrccos;
  double mfact,qfact;
  struct onode ninit[2];
  struct onode hn,tn;
  struct owire einit=elem;
  struct elemvisible ev=vp.vflag.ev;

  double nzmax=0.02;           /*** ujioka for Nz ***/

  if(mode==ONPRINTER)
  {
    qred=  0; qgreen=  0; qblue=  0; /*SHEAR COLOR*/
    mred=  0; mgreen=  0; mblue=  0; /*MOMENT COLOR*/
  }
  else if(mode==ONPREVIEW)
  {
    qred=  0; qgreen=  0; qblue=  0; /*SHEAR COLOR*/
    /*mred=  0; mgreen=  0; mblue=  0;*/ /*MOMENT COLOR*/
    mred=255; mgreen=105; mblue=  0; /*MOMENT COLOR*/
    /*mred=  0; mgreen=150; mblue=255;*/ /*MOMENT COLOR*/
  }
  else
  {
    qred=255; qgreen=150; qblue=150; /*SHEAR COLOR*/
    mred=150; mgreen=  0; mblue=255; /*MOMENT COLOR*/
  }

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

  angle=angleofgloballine(vp,ninit[0],ninit[1]);

  elength=distancedotdot(ninit[0],ninit[1]);

  drccos=directioncosine(einit.node[0]->d[GX],
                         einit.node[0]->d[GY],
                         einit.node[0]->d[GZ],
                         einit.node[1]->d[GX],
                         einit.node[1]->d[GY],
                         einit.node[1]->d[GZ],
                         einit.cangle);                  /*[DRCCOS]*/
  tdrccos=matrixtranspose(drccos,3);                          /*[T]*/

#if 0
if(mode!=ONPREVIEW)
{
#endif
  for(i=0;i<6;i++) /*0:Nz 1:Qx 2:Qy 3:Mz 4:Mx 5:My*/
  {
    if(i==0) SetTextColor(hdc,RGB(qred,qgreen,qblue));     /*SHEAR.*/
    if(i==3) SetTextColor(hdc,RGB(mred,mgreen,mblue));    /*MOMENT.*/

    for(j=1;j<=7;j++) /*ELEMENT TYPE.*/
    {
      if(elem.sect->type==j &&
         ev.stress[j][i])
      {
		hfont=(HFONT)GetCurrentObject(hdc,OBJ_FONT);
        GetObject(hfont,sizeof(LOGFONT),&lfont);
        lfont.lfEscapement = (int)(10.0*angle);
        hfont=CreateFontIndirect(&lfont);
		pfont=(HFONT)SelectObject(hdc,hfont);

        align = SetTextAlign(hdc,TA_BOTTOM |
                                 TA_CENTER |
                                 TA_NOUPDATECP);
#if 1
        /*** araki for Nz ***/
        if(i==0 && !(wdraw.childs+1)->vparam.vflag.mv.pagenum)
          {
            if(elem.stress[0][i]<=0)            /*tension*/
            {
              SetTextColor(hdc,RGB(255,0,150));
              drawglobalwire(hdc,vp,af,elem,255,0,150,
                                            255,0,150,elem.code,mode);
            }
            else                                /*compression*/
            {
              SetTextColor(hdc,RGB(0,150,255));
              drawglobalwire(hdc,vp,af,elem,0,150,255,
                                            0,150,255,elem.code,mode);
            }
          }
        /********************/
#endif

#if 0
        /*** ujioka for Nz ***/
        if(i==0)
          {
            if(elem.stress[0][i]<=-nzmax)            /*tension*/
            {
              SetTextColor(hdc,RGB(255,0,150));
              drawglobalwire(hdc,vp,af,elem,255,0,150,
                                            255,0,150,elem.code,mode);
            }
            else if(elem.stress[0][i]<=0.0)                                 /*tension*/
            {
              SetTextColor(hdc,RGB(
                             255,
                             255+elem.stress[0][i]/nzmax*(255-0),
                             255+elem.stress[0][i]/nzmax*(255-150)
                             ));
              drawglobalwire(hdc,vp,af,elem,
                             255,
                             255+elem.stress[0][i]/nzmax*(255-0),
                             255+elem.stress[0][i]/nzmax*(255-150),
                             255,
                             255+elem.stress[0][i]/nzmax*(255-0),
                             255+elem.stress[0][i]/nzmax*(255-150),
                             elem.code,mode);
            }
            else if(elem.stress[0][i]<=nzmax)                                 /*compression*/
            {
              SetTextColor(hdc,RGB(
                             255-elem.stress[0][i]/nzmax*(255-0),
                             255-elem.stress[0][i]/nzmax*(255-150),
                             255
                             ));
              drawglobalwire(hdc,vp,af,elem,
                             255-elem.stress[0][i]/nzmax*(255-0),
                             255-elem.stress[0][i]/nzmax*(255-150),
                             255,
                             255-elem.stress[0][i]/nzmax*(255-0),
                             255-elem.stress[0][i]/nzmax*(255-150),
                             255,
                             elem.code,mode);
            }
            else                                /*compression*/
            {
              SetTextColor(hdc,RGB(0,150,255));
              drawglobalwire(hdc,vp,af,elem,0,150,255,
                                            0,150,255,elem.code,mode);
            }
          }
        /********************/
#endif

#if 0
        /*** ujioka for Nz (GrayScale)***/
        if(i==0)
          {
            if(elem.stress[0][i]<=0.0)            /*tension*/
            {
              SetTextColor(hdc,RGB(255,0,150));
              drawglobalwire(hdc,vp,af,elem,255,0,150,
                                            255,0,150,elem.code,mode);
            }
            else if(elem.stress[0][i]<=nzmax)                                 /*compression*/
            {
              SetTextColor(hdc,RGB(
                             255-elem.stress[0][i]/nzmax*(255-20),
                             255-elem.stress[0][i]/nzmax*(255-20),
                             255-elem.stress[0][i]/nzmax*(255-20)
                             ));
              drawglobalwire(hdc,vp,af,elem,
                             255-elem.stress[0][i]/nzmax*(255-20),
                             255-elem.stress[0][i]/nzmax*(255-20),
                             255-elem.stress[0][i]/nzmax*(255-20),
                             255-elem.stress[0][i]/nzmax*(255-20),
                             255-elem.stress[0][i]/nzmax*(255-20),
                             255-elem.stress[0][i]/nzmax*(255-20),
                             elem.code,mode);
            }
            else                                /*compression*/
            {
              SetTextColor(hdc,RGB(20,20,20));
              drawglobalwire(hdc,vp,af,elem,20,20,20,
                                            20,20,20,elem.code,mode);
            }
          }
        /********************/
#endif

#if 0
        /*** ujioka for Nz ***/   //for 20211220.kyonan
        if(i==0)
          {
            if(elem.stress[0][i]<-0.9)
            {
              SetTextColor(hdc,RGB(255,0,150));
              drawglobalwire(hdc,vp,af,elem,255,0,150,
                                            255,0,150,elem.code,mode);
            }
            else if(elem.stress[0][i]<-0.7)
            {
              SetTextColor(hdc,RGB(255,150,50));
              drawglobalwire(hdc,vp,af,elem,255,150,50,
                                            255,150,50,elem.code,mode);
            }
            else if(elem.stress[0][i]<-0.5)
            {
              SetTextColor(hdc,RGB(0,255,0));
              drawglobalwire(hdc,vp,af,elem,0,255,0,
                                            0,255,0,elem.code,mode);
            }
            else if(elem.stress[0][i]<-0.3)
            {
              SetTextColor(hdc,RGB(0,255,150));
              drawglobalwire(hdc,vp,af,elem,0,255,150,
                                            0,255,150,elem.code,mode);
            }
            else
            {
              SetTextColor(hdc,RGB(0,150,255));
              drawglobalwire(hdc,vp,af,elem,0,150,255,
                                            0,150,255,elem.code,mode);
            }
          }
        /********************/
#endif
        if (i==0 && !(wdraw.childs+1)->vparam.vflag.mv.pagenum)
        {            /* Nz */
        /*sprintf(str,"%s %.5f",slabel[i],
                  globalunit*elem.stress[0][i]);*/
		sprintf(str,"%.2f",globalunit*elem.stress[0][i]);      /*HEAD.*/
        tn=setcoord(((0.5/**(6-i)*/)*elength),0.0,0.0);
		/*drawtextalignedonlocaldot(hdc,vp,tdrccos,ninit[0],tn,str,
									TEXT_CENTER,TEXT_TOP);*/
		drawtextonlocaldot(hdc,vp,tdrccos,ninit[0],tn,str);

        }

        else
        {
		/*sprintf(str,"%s %.5f",slabel[i],
                  globalunit*elem.stress[0][i]);*/
		sprintf(str,"%.2f",globalunit*elem.stress[0][i]);      /*HEAD.*/
//		sprintf(str,"%.5f",globalunit*elem.stress[0][i]);      /*HEAD.*/
        tn=setcoord(((0.25/**(6-i)*/)*elength),0.0,0.0);
		/*drawtextalignedonlocaldot(hdc,vp,tdrccos,ninit[0],tn,str,
									TEXT_CENTER,TEXT_TOP);*/
		drawtextonlocaldot(hdc,vp,tdrccos,ninit[0],tn,str);

		SetTextAlign(hdc,TA_TOP |
						 TA_CENTER |
						 TA_NOUPDATECP);
		/*sprintf(str,"%s %.5f",slabel[i],
				  globalunit*elem.stress[1][i]);*/
		sprintf(str,"%.2f",globalunit*elem.stress[1][i]);      /*TAIL.*/
//		sprintf(str,"%.5f",globalunit*elem.stress[1][i]);      /*TAIL.*/
		tn=setcoord(((1-0.25/**(i+1)*/)*elength),0.0,0.0);
		/*drawtextalignedonlocaldot(hdc,vp,tdrccos,ninit[0],tn,str,
								  TEXT_CENTER,TEXT_TOP);*/
        drawtextonlocaldot(hdc,vp,tdrccos,ninit[0],tn,str);
        }

        SetTextAlign(hdc,align);
        SelectObject(hdc,pfont);
        DeleteObject(hfont);
      }
    }
  }
#if 0
}
#endif

  hpenQ=CreatePen(PS_SOLID,1,RGB(qred,qgreen,qblue));  /*SHEAR PEN.*/
  ppen=(HPEN)SelectObject(hdc,hpenQ);

  for(i=0;i<=1;i++) /*HEAD,TAIL.*/
  {
    for(j=1;j<=7;j++) /*ELEMENT TYPE.*/
    {
      if(elem.sect->type==j)
      {
        if(ev.stress[j][0] && elem.stress[i][0]!=0.0) /*Nz*/
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
        if(ev.stress[j][1] && elem.stress[0][1]!=0.0) /*Qx*/
        {
          value=elem.stress[0][1];

          hn=setcoord(((0.05+0.9*i)*elength),(-qfact*0.5*value),0.0);
          tn=setcoord(((0.05+0.9*i)*elength),(qfact*0.5*value),0.0);
          drawlocalarrow(hdc,vp,tdrccos,ninit[0],hn,tn);
        }
        if(ev.stress[j][2] && elem.stress[0][2]!=0.0) /*Qy*/
        {
          value=elem.stress[0][2];

          hn=setcoord(((0.05+0.9*i)*elength),0.0,(-qfact*0.5*value));
          tn=setcoord(((0.05+0.9*i)*elength),0.0,(qfact*0.5*value));
          drawlocalarrow(hdc,vp,tdrccos,ninit[0],hn,tn);
        }
      }
    }
  }

  hpenM=CreatePen(PS_SOLID,1,RGB(mred,mgreen,mblue)); /*MOMENT PEN.*/
  SelectObject(hdc,hpenM);
  for(j=1;j<=7;j++) /*ELEMENT TYPE.*/
  {
    if(elem.sect->type==j)
    {
      if(vp.vflag.ev.stress[j][4]) /*Mx*/
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

          if(0.0<hn.d[EZ] && hn.d[EZ]<elength)
          {
            drawlocalline(hdc,vp,tdrccos,ninit[0],hn,tn);
            tn=hn;
          }
        }

        hn=setcoord(elength,0.0,(mfact*Mj));
        drawlocalline(hdc,vp,tdrccos,ninit[0],hn,tn);

        tn=setcoord(elength,0.0,0.0);
        drawlocalline(hdc,vp,tdrccos,ninit[0],hn,tn);
      }
      if(vp.vflag.ev.stress[j][5]) /*My*/
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

          if(0.0<hn.d[EZ] && hn.d[EZ]<elength)
          {
            drawlocalline(hdc,vp,tdrccos,ninit[0],hn,tn);
            tn=hn;
          }
        }

        hn=setcoord(elength,(-mfact*Mj),0.0);
        drawlocalline(hdc,vp,tdrccos,ninit[0],hn,tn);

        tn=setcoord(elength,0.0,0.0);
        drawlocalline(hdc,vp,tdrccos,ninit[0],hn,tn);
      }
    }
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
  ppen=(HPEN)SelectObject(hdc,hpen);

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

void drawwireshape(HDC hdc,
				  struct viewparam vp,
				  struct onode n1,struct onode n2,double cangle,
				  int ered,int egreen,int eblue)				   //honda
/*DRAW AXIS OF WIRE SHAPE.*/
{
  HPEN hpen,ppen;                                   /*HANDLE OF PEN*/
  int i;
  double elength;
  double **drccos,**tdrccos;
  double arrowlength=vp.dparam.eaxis;
  struct onode on,nx1,nx2,ny1,ny2;

  elength=distancedotdot(n1,n2);

  drccos=directioncosine(n1.d[GX],n1.d[GY],n1.d[GZ],
						 n2.d[GX],n2.d[GY],n2.d[GZ],
						 cangle);                        /*[DRCCOS]*/
  tdrccos=matrixtranspose(drccos,3);                          /*[T]*/

  hpen=CreatePen(PS_SOLID,1,RGB(ered,egreen,eblue));

  ppen=(HPEN)SelectObject(hdc,hpen);

  on=n1;

  nx1=setcoord((0.5*elength-0.5*arrowlength),0.5*arrowlength,1.5*arrowlength);
  nx2=setcoord((0.5*elength-0.5*arrowlength),0.5*arrowlength,-1.5*arrowlength);
  ny1=setcoord((0.5*elength-0.5*arrowlength),-0.5*arrowlength,-1.5*arrowlength);
  ny2=setcoord((0.5*elength-0.5*arrowlength),-0.5*arrowlength,1.5*arrowlength);
  drawlocalline(hdc,vp,tdrccos,on,nx1,nx2);
  drawlocalline(hdc,vp,tdrccos,on,nx2,ny1);
  drawlocalline(hdc,vp,tdrccos,on,ny1,ny2);
  drawlocalline(hdc,vp,tdrccos,on,ny2,nx1);
										   /*LOCAL ORIGIN.*/

  for(i=0;i<=2;i++) free(*(drccos+i));
  free(drccos);
  for(i=0;i<=2;i++) free(*(tdrccos+i));
  free(tdrccos);

  SelectObject(hdc,ppen);
  DeleteObject(hpen);

  return;
}/*drawwireshape*/


void drawarclmframe(HDC hdc,struct viewparam vp,
					struct arclmframe af,
					long int code,
					int mode)
{
  SIZE size;
  char str[256];
  char pcode1,pcode2,pcode3; /****SRCANMAX****/
  int i,flag,nx,ny;
  long int imax[2],imin[2];
  double hsize;
  struct owire elem;
  int ii,selected;
  int j; /****SRCANMAX****/


  if(vp.vflag.mv.pagetitle)
  {
    if(mode==ONPRINTER || mode==ONPREVIEW)
    {
      setfontformat(hdc,150,50,"ＭＳ 明朝",0,0,0);
    }
    if(mode==ONPRINTER)      SetTextColor(hdc,RGB(0,0,0));
    else if(mode==ONPREVIEW) SetTextColor(hdc,RGB(0,0,0));
    else                     SetTextColor(hdc,RGB(255,255,255));

	/****SRCANMAX****/
	if(mode==ONPRINTER && vp.vflag.ev.srcanmax)
	{
//	  setfontformat(hdc,120,40,"ＭＳ 明朝",0,0,0);
	  //sprintf(str,"4.3 : 柱，梁，壁，床の断面算定結果");
	  //sprintf(str,"4.3 : 柱，梁，壁の断面算定結果");
	  //sprintf(str,"4.3 : 柱，梁，ブレース，壁，床の断面算定結果");
	  sprintf(str,"4.3 : 柱，梁，ブレースの断面算定結果");
//	  TextOut(hdc,900,500,str,strlen(str));
	  TextOut(hdc,700,300,str,strlen(str));

	  setfontformat(hdc,(int)(1.5*gprn.jiheight),
						(int)(1.5*gprn.jiwidth),"ＭＳ 明朝",0,0,0);
	  sprintf(str,"各断面種別において、検定比が最大の部材は下図の通り。");
	  TextOut(hdc,900,700,str,strlen(str));
	  sprintf(str,"これらについてのみ断面算定の詳細を載せておく。");
	  TextOut(hdc,900,700+(int)(2.0*gprn.jiheight),str,strlen(str));
	}
	else
	{
	  sprintf(str,"%s",(wdraw.childs+1)->pagetitle);
	  if(mode==ONPRINTER) TextOut(hdc,700,300,str,strlen(str));
	  else                TextOut(hdc,120, 50,str,strlen(str));
	}
	/****************/

	//sprintf(str,"%s",(wdraw.childs+1)->pagetitle);
	//if(mode==ONPRINTER) TextOut(hdc,700,300,str,strlen(str));     //1200,500-700,300
    //else                TextOut(hdc, 120, 50,str,strlen(str));

  }


  /*PAGE NUMBER*/
  if(vp.vflag.mv.pagenum)
  {
	if(mode==ONPRINTER)
	{
	  setfontformat(hdc,80,30,"MS Mincho",0,0,0);
	}
	else if(mode==ONPREVIEW)
	{
	  setfontformat(hdc,80,30,"MS Mincho",0,0,0);
	}
	if(mode==ONPRINTER)      SetTextColor(hdc,RGB(0,0,0));
	else if(mode==ONPREVIEW) SetTextColor(hdc,RGB(0,0,0));
	else                     SetTextColor(hdc,RGB(255,255,255));

	sprintf(str,"%d.%d.%d",(wdraw.childs+1)->vparam.chapter,
                           (wdraw.childs+1)->vparam.section,
                           (wdraw.childs+1)->vparam.subsection);

	if(mode==ONPRINTER) TextOut(hdc,4400,6600,str,strlen(str));   //要調整
	else if(mode==ONPREVIEW) TextOut(hdc,1200,1000,str,strlen(str));
	else                     TextOut(hdc,1250,1000,str,strlen(str));
  }


  /*GET FRAME RANGE FOR DRAWING TEXTS*/
  if(mode==ONPRINTER || mode==ONPREVIEW)
  {
	setfontformat(hdc,(int)(1.5*gprn.jiheight),
                      (int)(1.5*gprn.jiwidth),"MS Mincho",0,0,0);
  }
  if(mode==ONPRINTER)      SetTextColor(hdc,RGB(0,0,0));
  else if(mode==ONPREVIEW) SetTextColor(hdc,RGB(0,0,0));
  else                     SetTextColor(hdc,RGB(255,255,255));

  imax[0]=0; imin[0]=0;
  imax[1]=0; imin[1]=0;

  flag=0;
  for(i=0;i<af.nnode;i++)
  {
    if(insiderange(*(af.ninit+i),vp.range))
    {
      if(!nodeontoscreen(*(af.ninit+i),&nx,&ny,vp)) return;

      if(flag==0) /*INITIAL*/
      {
        imax[0]=nx; imin[0]=nx;
        imax[1]=ny; imin[1]=ny;

        flag=1;
      }

      if(imax[0]<nx) imax[0]=nx;
      if(imax[1]<ny) imax[1]=ny;
      if(imin[0]>nx) imin[0]=nx;
      if(imin[1]>ny) imin[1]=ny;
    }
  }
  sprintf(str,"Text");
  GetTextExtentPoint32(hdc,str,strlen(str),&size);

  /*UPPER TEXTS*/
  imin[1]-=size.cy;
  if(vp.vflag.mv.title)
  {
    imin[1]-=size.cy;
    sprintf(str,"%s",(wdraw.childs+1)->title);
    TextOut(hdc,imin[0],imin[1],str,strlen(str));
  }

  /*LOWER TEXTS*/
  imax[1]+=size.cy;
  if(vp.vflag.nv.code)
  {
    imax[1]+=size.cy;
    sprintf(str,"節点番号");
    TextOut(hdc,imin[0],imax[1],str,strlen(str));
  }
  if(vp.vflag.ev.code)
  {
    imax[1]+=size.cy;
    sprintf(str,"部材番号");
    TextOut(hdc,imin[0],imax[1],str,strlen(str));
  }
  if(vp.vflag.ev.sectioncode)
  {
    imax[1]+=size.cy;
    sprintf(str,"断面番号");
    TextOut(hdc,imin[0],imax[1],str,strlen(str));
  }
  if(vp.vflag.nv.mcircle)
  {
    imax[1]+=size.cy;
    sprintf(str,"節点重量図");
    TextOut(hdc,imin[0],imax[1],str,strlen(str));
  }
  if(vp.vflag.nv.mvalue)
  {
    imax[1]+=size.cy;
    sprintf(str,"節点重量値");
    TextOut(hdc,imin[0],imax[1],str,strlen(str));
  }
  if(vp.vflag.ev.ecircle)
  {
    imax[1]+=size.cy;
    sprintf(str,"吸収エネルギー図");
    TextOut(hdc,imin[0],imax[1],str,strlen(str));
  }
  if(vp.vflag.ev.evalue)
  {
    imax[1]+=size.cy;
    sprintf(str,"吸収エネルギー値");
    TextOut(hdc,imin[0],imax[1],str,strlen(str));
  }
  if(vp.vflag.ev.deformation)
  {
    imax[1]+=size.cy;
    sprintf(str,"変形図");
    TextOut(hdc,imin[0],imax[1],str,strlen(str));
  }
  if(vp.vflag.nv.disps[0])
  {
    imax[1]+=size.cy;
    sprintf(str,"Ｘ方向変位 [cm]");
    TextOut(hdc,imin[0],imax[1],str,strlen(str));
  }
  if(vp.vflag.nv.disps[1])
  {
    imax[1]+=size.cy;
    sprintf(str,"Ｙ方向変位 [cm]");
    TextOut(hdc,imin[0],imax[1],str,strlen(str));
  }
  if(vp.vflag.nv.disps[2])
  {
    imax[1]+=size.cy;
    sprintf(str,"Ｚ方向変位 [cm]");
    TextOut(hdc,imin[0],imax[1],str,strlen(str));
  }
  /*if(vp.vflag.nv.mcircle)
  {
    imax[1]+=size.cy;
    sprintf(str,"節点重量図");
    TextOut(hdc,imin[0],imax[1],str,strlen(str));
  }*/
  /*if(vp.vflag.nv.mvalue)
  {
    imax[1]+=size.cy;
    if(globalunit==1.0)         sprintf(str,"節点重量値 [tf]");
    else if(globalunit==SIUNIT) sprintf(str,"節点重量値 [kN]");
    else                        sprintf(str,"節点重量値");
    TextOut(hdc,imin[0],imax[1],str,strlen(str));
  }*/
  if(vp.vflag.ev.axis)
  {
    imax[1]+=size.cy;
    sprintf(str,"部材座標軸");
    TextOut(hdc,imin[0],imax[1],str,strlen(str));
  }
  for(i=0;i<8;i++)
  {
    if(vp.vflag.ev.stress[i][0])
    {
      imax[1]+=size.cy;
      if(globalunit==1.0)         sprintf(str,"軸力 [tf]");
      else if(globalunit==SIUNIT) sprintf(str,"軸力 [kN]");
      else                        sprintf(str,"軸力");
      TextOut(hdc,imin[0],imax[1],str,strlen(str));
      break;
    }
  }
  for(i=0;i<8;i++)
  {
    if(vp.vflag.ev.stress[i][1] || vp.vflag.ev.stress[i][2])
    {
      imax[1]+=size.cy;
      if(globalunit==1.0)         sprintf(str,"せん断力 [tf]");
      else if(globalunit==SIUNIT) sprintf(str,"せん断力 [kN]");
      else                        sprintf(str,"せん断力");
      TextOut(hdc,imin[0],imax[1],str,strlen(str));
      break;
    }
  }
  for(i=0;i<8;i++)
  {
    if(vp.vflag.ev.stress[i][3])
    {
      imax[1]+=size.cy;
      if(globalunit==1.0)         sprintf(str,"ねじりモーメント [tfm]");
      else if(globalunit==SIUNIT) sprintf(str,"ねじりモーメント [kNm]");
      else                        sprintf(str,"ねじりモーメント");
      TextOut(hdc,imin[0],imax[1],str,strlen(str));
      break;
    }
  }
  for(i=0;i<8;i++)
  {
    if(vp.vflag.ev.stress[i][4] || vp.vflag.ev.stress[i][5])
    {
      imax[1]+=size.cy;
      if(globalunit==1.0)         sprintf(str,"曲げモーメント [tfm]");
      else if(globalunit==SIUNIT) sprintf(str,"曲げモーメント [kNm]");
      else                        sprintf(str,"曲げモーメント");
      TextOut(hdc,imin[0],imax[1],str,strlen(str));
      break;
    }
  }
  /*if(vp.vflag.ev.cmqline)
  {
    imax[1]+=size.cy;
    sprintf(str,"分割線表示");
    TextOut(hdc,imin[0],imax[1],str,strlen(str));
  }*/

//////////////////////////////////////////////////////by MIHARA for Safety Level

  /*if(vp.vflag.ev.srcanrate)*/
  if(vp.vflag.ev.srcanrate || vp.vflag.ev.srcancolor)   /*Modified by Ujioka*/
  {
    imax[1]+=size.cy;
	sprintf(str,"安全率 (断面検定比図)");
	TextOut(hdc,imin[0],imax[1],str,strlen(str));

	if(mode==ONPRINTER)
	{
	sprintf(str,"安全率の凡例");/*SRC COLOR MAP*/
#if 0
	TextOut(hdc,3500,imax[1]+size.cy,str,strlen(str));
	sprintf(str,"　 : ≧1.0");
	TextOut(hdc,3500,imax[1]+3*size.cy,str,strlen(str));
	sprintf(str,"　 : 0.9～1.0");
	TextOut(hdc,3500,imax[1]+4*size.cy,str,strlen(str));
	sprintf(str,"　 : 0.8～0.9");
	TextOut(hdc,3500,imax[1]+5*size.cy,str,strlen(str));
	sprintf(str,"　 : 0.7～0.8");
	TextOut(hdc,3500,imax[1]+6*size.cy,str,strlen(str));
	sprintf(str,"　 : 0.6～0.7");
	TextOut(hdc,3500,imax[1]+7*size.cy,str,strlen(str));
	sprintf(str,"　 : 0.5～0.6");
	TextOut(hdc,3500,imax[1]+8*size.cy,str,strlen(str));
	sprintf(str,"　 : ＜0.5");
	TextOut(hdc,3500,imax[1]+9*size.cy,str,strlen(str));
#endif
#if 1
	TextOut(hdc,3500,imax[1]+size.cy,str,strlen(str));
	sprintf(str,"　 : ≧1.0");
	TextOut(hdc,3500,imax[1]+3*size.cy,str,strlen(str));
	sprintf(str,"　 : 0.8～1.0");
	TextOut(hdc,3500,imax[1]+4*size.cy,str,strlen(str));
	sprintf(str,"　 : 0.6～0.8");
	TextOut(hdc,3500,imax[1]+5*size.cy,str,strlen(str));
	sprintf(str,"　 : 0.4～0.6");
	TextOut(hdc,3500,imax[1]+6*size.cy,str,strlen(str));
	sprintf(str,"　 : 0.2～0.4");
	TextOut(hdc,3500,imax[1]+7*size.cy,str,strlen(str));
	sprintf(str,"　 : 0.1～0.2");
	TextOut(hdc,3500,imax[1]+8*size.cy,str,strlen(str));
	sprintf(str,"　 : ＜0.1");
	TextOut(hdc,3500,imax[1]+9*size.cy,str,strlen(str));
#endif
#if 0
/*******FOR BUCKLING SAFETY******/
	sprintf(str,"　 : ≧1.0");
    TextOut(hdc,3500,imax[1]+3*size.cy,str,strlen(str));
    sprintf(str,"　 : 0.5～1.0");
    TextOut(hdc,3500,imax[1]+4*size.cy,str,strlen(str));
    sprintf(str,"　 : 0.2～0.5");
    TextOut(hdc,3500,imax[1]+5*size.cy,str,strlen(str));
    sprintf(str,"　 : 0.1～0.2");
    TextOut(hdc,3500,imax[1]+6*size.cy,str,strlen(str));
    sprintf(str,"　 : 0.05～0.1");
    TextOut(hdc,3500,imax[1]+7*size.cy,str,strlen(str));
    sprintf(str,"　 : 0.02～0.05");
    TextOut(hdc,3500,imax[1]+8*size.cy,str,strlen(str));
    sprintf(str,"　 : ＜0.02");
    TextOut(hdc,3500,imax[1]+9*size.cy,str,strlen(str));
#endif
    SetTextColor(hdc,RGB(255,0,150));
    sprintf(str,"■");
    TextOut(hdc,3500,imax[1]+3*size.cy,str,strlen(str));
    SetTextColor(hdc,RGB(255,150,50));
    sprintf(str,"■");
    TextOut(hdc,3500,imax[1]+4*size.cy,str,strlen(str));
    SetTextColor(hdc,RGB(150,150,0));
    sprintf(str,"■");
    TextOut(hdc,3500,imax[1]+5*size.cy,str,strlen(str));
    SetTextColor(hdc,RGB(150,255,0));
    sprintf(str,"■");
    TextOut(hdc,3500,imax[1]+6*size.cy,str,strlen(str));
    SetTextColor(hdc,RGB(0,255,0));
    sprintf(str,"■");
    TextOut(hdc,3500,imax[1]+7*size.cy,str,strlen(str));
    SetTextColor(hdc,RGB(0,255,150));
    sprintf(str,"■");
	TextOut(hdc,3500,imax[1]+8*size.cy,str,strlen(str));
    SetTextColor(hdc,RGB(0,150,255));
    sprintf(str,"■");
    TextOut(hdc,3500,imax[1]+9*size.cy,str,strlen(str));
    }

    else if (mode==ONPREVIEW)
    {
        sprintf(str,"安全率の凡例");
    TextOut(hdc,imin[0]+1000,imax[1]+size.cy,str,strlen(str));
#if 0
	sprintf(str,"　 : ≧1.0");
	TextOut(hdc,imin[0]+1000,imax[1]+3*size.cy,str,strlen(str));
	sprintf(str,"　 : 0.9～1.0");
	TextOut(hdc,imin[0]+1000,imax[1]+4*size.cy,str,strlen(str));
	sprintf(str,"　 : 0.8～0.9");
	TextOut(hdc,imin[0]+1000,imax[1]+5*size.cy,str,strlen(str));
	sprintf(str,"　 : 0.7～0.8");
	TextOut(hdc,imin[0]+1000,imax[1]+6*size.cy,str,strlen(str));
	sprintf(str,"　 : 0.6～0.7");
	TextOut(hdc,imin[0]+1000,imax[1]+7*size.cy,str,strlen(str));
	sprintf(str,"　 : 0.5～0.6");
	TextOut(hdc,imin[0]+1000,imax[1]+8*size.cy,str,strlen(str));
	sprintf(str,"　 : ＜0.5");
	TextOut(hdc,imin[0]+1000,imax[1]+9*size.cy,str,strlen(str));
#endif
#if 1
	sprintf(str,"　 : ≧1.0");
	TextOut(hdc,imin[0]+1000,imax[1]+3*size.cy,str,strlen(str));
	sprintf(str,"　 : 0.8～1.0");
	TextOut(hdc,imin[0]+1000,imax[1]+4*size.cy,str,strlen(str));
	sprintf(str,"　 : 0.6～0.8");
	TextOut(hdc,imin[0]+1000,imax[1]+5*size.cy,str,strlen(str));
	sprintf(str,"　 : 0.4～0.6");
	TextOut(hdc,imin[0]+1000,imax[1]+6*size.cy,str,strlen(str));
	sprintf(str,"　 : 0.2～0.4");
	TextOut(hdc,imin[0]+1000,imax[1]+7*size.cy,str,strlen(str));
	sprintf(str,"　 : 0.1～0.2");
	TextOut(hdc,imin[0]+1000,imax[1]+8*size.cy,str,strlen(str));
	sprintf(str,"　 : ＜0.1");
	TextOut(hdc,imin[0]+1000,imax[1]+9*size.cy,str,strlen(str));
#endif
#if 0
/*******FOR BUCKLING SAFETY******/
    sprintf(str,"　 : ≧1.0");
    TextOut(hdc,imin[0]+1000,imax[1]+3*size.cy,str,strlen(str));
    sprintf(str,"　 : 0.5～1.0");
    TextOut(hdc,imin[0]+1000,imax[1]+4*size.cy,str,strlen(str));
    sprintf(str,"　 : 0.2～0.5");
    TextOut(hdc,imin[0]+1000,imax[1]+5*size.cy,str,strlen(str));
    sprintf(str,"　 : 0.1～0.2");
    TextOut(hdc,imin[0]+1000,imax[1]+6*size.cy,str,strlen(str));
    sprintf(str,"　 : 0.05～0.1");
    TextOut(hdc,imin[0]+1000,imax[1]+7*size.cy,str,strlen(str));
    sprintf(str,"　 : 0.02～0.5");
    TextOut(hdc,imin[0]+1000,imax[1]+8*size.cy,str,strlen(str));
    sprintf(str,"　 : ＜0.02");
    TextOut(hdc,imin[0]+1000,imax[1]+9*size.cy,str,strlen(str));
#endif

    SetTextColor(hdc,RGB(255,0,150));
    sprintf(str,"■");
    TextOut(hdc,imin[0]+1000,imax[1]+3*size.cy,str,strlen(str));
    SetTextColor(hdc,RGB(255,150,50));
    sprintf(str,"■");
    TextOut(hdc,imin[0]+1000,imax[1]+4*size.cy,str,strlen(str));
    SetTextColor(hdc,RGB(255,255,0));
    sprintf(str,"■");
    TextOut(hdc,imin[0]+1000,imax[1]+5*size.cy,str,strlen(str));
    SetTextColor(hdc,RGB(150,255,0));
    sprintf(str,"■");
    TextOut(hdc,imin[0]+1000,imax[1]+6*size.cy,str,strlen(str));
    SetTextColor(hdc,RGB(0,255,0));
    sprintf(str,"■");
    TextOut(hdc,imin[0]+1000,imax[1]+7*size.cy,str,strlen(str));
    SetTextColor(hdc,RGB(0,255,150));
    sprintf(str,"■");
    TextOut(hdc,imin[0]+1000,imax[1]+8*size.cy,str,strlen(str));
    SetTextColor(hdc,RGB(0,150,255));
    sprintf(str,"■");
    TextOut(hdc,imin[0]+1000,imax[1]+9*size.cy,str,strlen(str));
    }
    
    else
    {
	sprintf(str,"安全率の凡例");
//    TextOut(hdc,imin[0]+250,imax[1]+size.cy,str,strlen(str));        //ujioka
    TextOut(hdc,imin[0]+500,imax[1]+size.cy,str,strlen(str));
#if 0
	sprintf(str,"　 : ≧1.0");
	TextOut(hdc,imin[0]+500,imax[1]+3*size.cy,str,strlen(str));
	sprintf(str,"　 : 0.9～1.0");
	TextOut(hdc,imin[0]+500,imax[1]+4*size.cy,str,strlen(str));
	sprintf(str,"　 : 0.8～0.9");
	TextOut(hdc,imin[0]+500,imax[1]+5*size.cy,str,strlen(str));
	sprintf(str,"　 : 0.7～0.8");
	TextOut(hdc,imin[0]+500,imax[1]+6*size.cy,str,strlen(str));
	sprintf(str,"　 : 0.6～0.7");
	TextOut(hdc,imin[0]+500,imax[1]+7*size.cy,str,strlen(str));
	sprintf(str,"　 : 0.5～0.6");
	TextOut(hdc,imin[0]+500,imax[1]+8*size.cy,str,strlen(str));
	sprintf(str,"　 : ＜0.5");
	TextOut(hdc,imin[0]+500,imax[1]+9*size.cy,str,strlen(str));
#endif
#if 1
	sprintf(str,"　 : ≧1.0");
	TextOut(hdc,imin[0]+500,imax[1]+3*size.cy,str,strlen(str));
	sprintf(str,"　 : 0.8～1.0");
	TextOut(hdc,imin[0]+500,imax[1]+4*size.cy,str,strlen(str));
	sprintf(str,"　 : 0.6～0.8");
	TextOut(hdc,imin[0]+500,imax[1]+5*size.cy,str,strlen(str));
	sprintf(str,"　 : 0.4～0.6");
	TextOut(hdc,imin[0]+500,imax[1]+6*size.cy,str,strlen(str));
	sprintf(str,"　 : 0.2～0.4");
	TextOut(hdc,imin[0]+500,imax[1]+7*size.cy,str,strlen(str));
	sprintf(str,"　 : 0.1～0.2");
	TextOut(hdc,imin[0]+500,imax[1]+8*size.cy,str,strlen(str));
	sprintf(str,"　 : ＜0.1");
	TextOut(hdc,imin[0]+500,imax[1]+9*size.cy,str,strlen(str));
#endif

#if 0
/*******FOR BUCKLING SAFETY******/
    sprintf(str,"　 : ≧1.0");
	TextOut(hdc,imin[0]+500,imax[1]+3*size.cy,str,strlen(str));
    sprintf(str,"　 : 0.5～1.0");
    TextOut(hdc,imin[0]+500,imax[1]+4*size.cy,str,strlen(str));
    sprintf(str,"　 : 0.2～0.5");
    TextOut(hdc,imin[0]+500,imax[1]+5*size.cy,str,strlen(str));
    sprintf(str,"　 : 0.1～0.2");
    TextOut(hdc,imin[0]+500,imax[1]+6*size.cy,str,strlen(str));
    sprintf(str,"　 : 0.05～0.1");
    TextOut(hdc,imin[0]+500,imax[1]+7*size.cy,str,strlen(str));
    sprintf(str,"　 : 0.02～0.05");
    TextOut(hdc,imin[0]+500,imax[1]+8*size.cy,str,strlen(str));
    sprintf(str,"　 : ＜0.02");
    TextOut(hdc,imin[0]+500,imax[1]+9*size.cy,str,strlen(str));
#endif

    SetTextColor(hdc,RGB(255,0,150));
    sprintf(str,"■");
    TextOut(hdc,imin[0]+500,imax[1]+3*size.cy,str,strlen(str));
    SetTextColor(hdc,RGB(255,150,50));
    sprintf(str,"■");
    TextOut(hdc,imin[0]+500,imax[1]+4*size.cy,str,strlen(str));
    SetTextColor(hdc,RGB(255,255,0));
    sprintf(str,"■");
    TextOut(hdc,imin[0]+500,imax[1]+5*size.cy,str,strlen(str));
    SetTextColor(hdc,RGB(150,255,0));
    sprintf(str,"■");
    TextOut(hdc,imin[0]+500,imax[1]+6*size.cy,str,strlen(str));
    SetTextColor(hdc,RGB(0,255,0));
    sprintf(str,"■");
    TextOut(hdc,imin[0]+500,imax[1]+7*size.cy,str,strlen(str));
    SetTextColor(hdc,RGB(0,255,150));
    sprintf(str,"■");
    TextOut(hdc,imin[0]+500,imax[1]+8*size.cy,str,strlen(str));
    SetTextColor(hdc,RGB(0,150,255));
    sprintf(str,"■");
    TextOut(hdc,imin[0]+500,imax[1]+9*size.cy,str,strlen(str));
    }

    if(mode==ONPRINTER)      SetTextColor(hdc,RGB(0,0,0));
    else if(mode==ONPREVIEW) SetTextColor(hdc,RGB(0,0,0));
    else                     SetTextColor(hdc,RGB(255,255,255));
  }

////////////////////////////////////////////////////////////////////////////////

  if(vp.vflag.ev.etype[0] ||
     vp.vflag.ev.etype[1] ||
     vp.vflag.ev.etype[2] ||
     vp.vflag.ev.etype[3] ||
     vp.vflag.ev.etype[4] ||
     vp.vflag.ev.etype[5] ||
     vp.vflag.ev.etype[6] ||
     vp.vflag.ev.etype[7])
  {
    flag=0;
    imax[1]+=size.cy;
    sprintf(str,"表示部材：");
    /*if(vp.vflag.ev.etype[0])
    {
      strcat(str,"補助");
      flag=1;
    }*/
    if(vp.vflag.ev.etype[1])
    {
      if(flag) strcat(str,", ");
      strcat(str,"柱");
      flag=1;
    }
    if(vp.vflag.ev.etype[2])
    {
      if(flag) strcat(str,", ");
      strcat(str,"梁");
      flag=1;
    }
    /*if(vp.vflag.ev.etype[3])
    {
      if(flag) strcat(str,", ");
      strcat(str,"小梁");
      flag=1;
    }*/
    if(vp.vflag.ev.etype[4])
    {
      if(flag) strcat(str,", ");
      strcat(str,"ブレース");
      flag=1;
    }
    if(vp.vflag.ev.etype[5])
    {
      if(flag) strcat(str,", ");
      strcat(str,"壁");
      flag=1;
    }
    if(vp.vflag.ev.etype[6])
    {
      if(flag) strcat(str,", ");
      strcat(str,"床");
      flag=1;
    }
    /*if(vp.vflag.ev.etype[7])
    {
      if(flag) strcat(str,", ");
      strcat(str,"張壁");
      flag=1;
    }*/
    TextOut(hdc,imin[0],imax[1],str,strlen(str));
  }

  imax[1]+=size.cy;
  if(vp.vflag.mv.inputfile)
  {
    imax[1]+=size.cy;
    sprintf(str,"Input  : %s",(wdraw.childs+1)->inpfile);
    TextOut(hdc,imin[0],imax[1],str,strlen(str));
  }
  if(vp.vflag.mv.outputfile)
  {
    imax[1]+=size.cy;
    sprintf(str,"Output : %s",(wdraw.childs+1)->otpfile);
    TextOut(hdc,imin[0],imax[1],str,strlen(str));
  }
  if(vp.vflag.mv.view)
  {
    imax[1]+=size.cy;
    sprintf(str,"Focus : %.3f %.3f %.3f",
            vp.focus.d[GX],vp.focus.d[GY],vp.focus.d[GZ]);
    TextOut(hdc,imin[0],imax[1],str,strlen(str));
    imax[1]+=size.cy;
	sprintf(str,"Phi=%.3f Theta=%.3f",vp.phi,vp.theta);
	TextOut(hdc,imin[0],imax[1],str,strlen(str));
    imax[1]+=size.cy;
    sprintf(str,"R=%.3f L=%.3f",vp.r,vp.odv);
    TextOut(hdc,imin[0],imax[1],str,strlen(str));
    imax[1]+=size.cy;
    sprintf(str,"Dfact=%.3f",vp.dparam.dfact);
    TextOut(hdc,imin[0],imax[1],str,strlen(str));
    imax[1]+=size.cy;
    sprintf(str,"Mfact=%.3f",vp.dparam.mfact);
    TextOut(hdc,imin[0],imax[1],str,strlen(str));
  }
if(vp.vflag.nv.conffig)
  {
    for(i=0;i<(af.nnode);i++) /*CINFFIG*/
    {
      /*loff=6*((af.nodes+i)->loff);*/

      if(mode==ONPRINTER)
      {
        hsize=vp.dparam.hsize;
        vp.dparam.hsize*=5.0;
        drawglobalconfigofaf(hdc,vp,af,*(af.nodes+i),mode,6*((af.nodes+i)->loff));
        vp.dparam.hsize=hsize;
      }
      else if(mode==ONPREVIEW)
      {
        hsize=vp.dparam.hsize;
        vp.dparam.hsize*=5.0;
        drawglobalconfigofaf(hdc,vp,af,*(af.nodes+i),mode,6*((af.nodes+i)->loff));
        vp.dparam.hsize=hsize;
      }
      else drawglobalconfigofaf(hdc,vp,af,*(af.nodes+i),mode,6*((af.nodes+i)->loff));

//      vp.dparam.hsize=hsize;
    }
    imax[1]+=size.cy;
    imax[1]+=size.cy;
    sprintf(str,"▲印は支点位置を表す");
    TextOut(hdc,imin[0],imax[1],str,strlen(str));
    imax[1]+=size.cy;
    sprintf(str,"(X,Y,Z,θx,θy,θzのうち，X,Y,Zを拘束している)");
    TextOut(hdc,imin[0],imax[1],str,strlen(str));
  }

//printrange 2007.12.04/////////////////////////////////////////////////////////
  if(vp.vflag.pv.printrange)
  {
   if((mode==ONPRINTER)||(mode==ONPREVIEW))
   {}
   else
   {
   drawprintrange(hdc,vp);
   }
  }
////////////////////////////////////////////////////////////////////////////////

  /*FONT*/
  if(mode==ONPRINTER || mode==ONPREVIEW)
  {
	setfontformat(hdc,gprn.jiheight,gprn.jiwidth,"MS Mincho",0,0,0);
  }

  if(vp.vflag.axis==1) /*DRAW GLOBAL AXIS.*/
  {
    /*if(mode==ONPREVIEW) drawglobalaxis(hdc,vp,255,255,0);*//*REVERSE*/
	if(mode==ONPREVIEW) drawglobalaxis(hdc,vp,0,0,255);/*BLUE*/
	else                drawglobalaxis(hdc,vp,0,0,255);
  }

  drawarclmnodes(hdc,vp,af,code,mode);

  for(i=0;i<af.nelem;i++)
  {
    /*currentpivot(i+1,af.nelem);*/
    /*currentpivot((af.elems+i)->code,af.nelem);*/

/*
    if((af.elems+i)->sect->dflag==1 &&
       vp.vflag.ev.etype[(af.elems+i)->sect->type]==1)
*/

/*SRC COLOR MAP*/
    if(!globalcondensationflag && ((af.elems+i)->sect->dflag==1 &&
       vp.vflag.ev.etype[(af.elems+i)->sect->type]==1))
    {
	  if(globalstatus==SELECTELEMENT && (af.elems+i)->code==code)
      {
        drawglobalwire(hdc,vp,af,*(af.elems+i),255,0,255,
                                               255,0,255,code,mode);
      }
	  else if(globalstatus==SELECTSECTION &&
              (af.elems+i)->sect->code==code)
      {
        drawglobalwire(hdc,vp,af,*(af.elems+i),255,0,255,
                                               255,0,255,code,mode);
      }
      else if(vp.vflag.ev.srcancolor &&
              ((af.elems+i)->srate[0]>1.0 || /*1.0*//*5.3*/
               (af.elems+i)->srate[1]>1.0 ||
               (af.elems+i)->srate[2]>1.0 ||
               (af.elems+i)->srate[3]>1.0 ))
      {
        drawglobalwire(hdc,vp,af,*(af.elems+i),255,0,150,
                                               255,0,150,code,mode);

/*drawglobalwire(hdc,vp,af,*(af.elems+i),255,150,50,
									   255,150,50,code,mode);*/
      }
#if 0
/*******FOR BUCKLING SAFETY******/
      else if(vp.vflag.ev.srcancolor &&
			  ((af.elems+i)->srate[0]>0.5 || /*0.9*//*0.8*//*1.3*/
               (af.elems+i)->srate[1]>0.5 ||
               (af.elems+i)->srate[2]>0.5 ||
               (af.elems+i)->srate[3]>0.5))
      {
        drawglobalwire(hdc,vp,af,*(af.elems+i),255,150,50,
                                               255,150,50,code,mode);
      }
      else if(vp.vflag.ev.srcancolor &&
			  ((af.elems+i)->srate[0]>0.2 ||
			   (af.elems+i)->srate[1]>0.2 ||
               (af.elems+i)->srate[2]>0.2 ||
               (af.elems+i)->srate[3]>0.2))
      {
        if(mode==ONPRINTER)
        {
          drawglobalwire(hdc,vp,af,*(af.elems+i),150,150,0,
                                                 150,150,0,code,mode);
        }
        else
        {
          drawglobalwire(hdc,vp,af,*(af.elems+i),255,255,0,
                                                 255,255,0,code,mode);
        }
      }
      else if(vp.vflag.ev.srcancolor &&
              ((af.elems+i)->srate[0]>0.1 || /*0.7*//*0.6*//*0.5*/
               (af.elems+i)->srate[1]>0.1 ||
               (af.elems+i)->srate[2]>0.1 ||
               (af.elems+i)->srate[3]>0.1))
      {
        drawglobalwire(hdc,vp,af,*(af.elems+i),150,255,0,
                                               150,255,0,code,mode);
      }
      else if(vp.vflag.ev.srcancolor &&
              ((af.elems+i)->srate[0]>0.05 || /*0.6*//*0.5*//*0.3*/
               (af.elems+i)->srate[1]>0.05 ||
               (af.elems+i)->srate[2]>0.05 ||
               (af.elems+i)->srate[3]>0.05))
      {
        drawglobalwire(hdc,vp,af,*(af.elems+i),0,255,0,
                                               0,255,0,code,mode);
      }
      else if(vp.vflag.ev.srcancolor &&
              ((af.elems+i)->srate[0]>0.02 || /*0.5*//*0.4*//*0.2*/
               (af.elems+i)->srate[1]>0.02 ||
               (af.elems+i)->srate[2]>0.02 ||
               (af.elems+i)->srate[3]>0.02))
      {
        drawglobalwire(hdc,vp,af,*(af.elems+i),0,255,150,
                                               0,255,150,code,mode);
      }
      else if(vp.vflag.ev.srcancolor)
      {
        drawglobalwire(hdc,vp,af,*(af.elems+i),80,50,200,
											   80,50,200,code,mode);
      }
#endif

#if 0
	  else if(vp.vflag.ev.srcancolor &&
			  ((af.elems+i)->srate[0]>0.9 || /*0.9*//*0.8*//*1.3*/
               (af.elems+i)->srate[1]>0.9 ||
			   (af.elems+i)->srate[2]>0.9 ||
			   (af.elems+i)->srate[3]>0.9))
	  {
		drawglobalwire(hdc,vp,af,*(af.elems+i),255,150,50,
											   255,150,50,code,mode);
	  }
	  else if(vp.vflag.ev.srcancolor &&
			  ((af.elems+i)->srate[0]>0.8 || /*0.8*/
			   (af.elems+i)->srate[1]>0.8 ||
			   (af.elems+i)->srate[2]>0.8 ||
			   (af.elems+i)->srate[3]>0.8))
	  {
		if(mode==ONPRINTER)
		{
		  drawglobalwire(hdc,vp,af,*(af.elems+i),150,150,0,
												 150,150,0,code,mode);
		}
        else
		{
		  drawglobalwire(hdc,vp,af,*(af.elems+i),255,255,0,
												 255,255,0,code,mode);
        }
      }
	  else if(vp.vflag.ev.srcancolor &&
              ((af.elems+i)->srate[0]>0.7 || /*0.7*//*0.6*//*0.5*/
			   (af.elems+i)->srate[1]>0.7 ||
               (af.elems+i)->srate[2]>0.7 ||
               (af.elems+i)->srate[3]>0.7))
	  {
        drawglobalwire(hdc,vp,af,*(af.elems+i),150,255,0,
                                               150,255,0,code,mode);
	  }
      else if(vp.vflag.ev.srcancolor &&
              ((af.elems+i)->srate[0]>0.6 || /*0.6*//*0.5*//*0.3*/
			   (af.elems+i)->srate[1]>0.6 ||
               (af.elems+i)->srate[2]>0.6 ||
               (af.elems+i)->srate[3]>0.6))
	  {
        drawglobalwire(hdc,vp,af,*(af.elems+i),0,255,0,
                                               0,255,0,code,mode);
	  }
      else if(vp.vflag.ev.srcancolor &&
              ((af.elems+i)->srate[0]>0.5 || /*0.5*//*0.4*//*0.2*/
			   (af.elems+i)->srate[1]>0.5 ||
               (af.elems+i)->srate[2]>0.5 ||
               (af.elems+i)->srate[3]>0.5))
      {
        drawglobalwire(hdc,vp,af,*(af.elems+i),0,255,150,
											   0,255,150,code,mode);
      }
      else if(vp.vflag.ev.srcancolor)
	  {
        drawglobalwire(hdc,vp,af,*(af.elems+i),80,50,200,
											   80,50,200,code,mode);
	  }
#endif
#if 1
	  else if(vp.vflag.ev.srcancolor &&
			  ((af.elems+i)->srate[0]>0.8 ||
			   (af.elems+i)->srate[1]>0.8 ||
			   (af.elems+i)->srate[2]>0.8 ||
			   (af.elems+i)->srate[3]>0.8))
	  {
		drawglobalwire(hdc,vp,af,*(af.elems+i),255,150,50,
											   255,150,50,code,mode);
	  }
	  else if(vp.vflag.ev.srcancolor &&
			  ((af.elems+i)->srate[0]>0.6 ||
			   (af.elems+i)->srate[1]>0.6 ||
			   (af.elems+i)->srate[2]>0.6 ||
			   (af.elems+i)->srate[3]>0.6))
	  {
		if(mode==ONPRINTER)
		{
		  drawglobalwire(hdc,vp,af,*(af.elems+i),150,150,0,
												 150,150,0,code,mode);
		}
		else
		{
		  drawglobalwire(hdc,vp,af,*(af.elems+i),255,255,0,
												 255,255,0,code,mode);
		}
	  }
	  else if(vp.vflag.ev.srcancolor &&
			  ((af.elems+i)->srate[0]>0.4 ||
			   (af.elems+i)->srate[1]>0.4 ||
			   (af.elems+i)->srate[2]>0.4 ||
			   (af.elems+i)->srate[3]>0.4))
	  {
		drawglobalwire(hdc,vp,af,*(af.elems+i),150,255,0,
											   150,255,0,code,mode);
	  }
	  else if(vp.vflag.ev.srcancolor &&
			  ((af.elems+i)->srate[0]>0.2 ||
			   (af.elems+i)->srate[1]>0.2 ||
			   (af.elems+i)->srate[2]>0.2 ||
			   (af.elems+i)->srate[3]>0.2))
	  {
		drawglobalwire(hdc,vp,af,*(af.elems+i),0,255,0,
											   0,255,0,code,mode);
	  }
	  else if(vp.vflag.ev.srcancolor &&
			  ((af.elems+i)->srate[0]>0.1 ||
			   (af.elems+i)->srate[1]>0.1 ||
			   (af.elems+i)->srate[2]>0.1 ||
			   (af.elems+i)->srate[3]>0.1))
	  {
		drawglobalwire(hdc,vp,af,*(af.elems+i),0,255,150,
											   0,255,150,code,mode);
	  }
	  else if(vp.vflag.ev.srcancolor)
	  {
		drawglobalwire(hdc,vp,af,*(af.elems+i),80,50,200,
											   80,50,200,code,mode);
	  }
#endif

#if 1
	  /****SRCANMAX****/
	  /*else if(vp.vflag.ev.srcanmax &&
			  (af.elems+i)->sect->smax==i)*/
	  else if(vp.vflag.ev.srcanmax)
	  {
		for(j=0;j<(af.nosect);j++)
		{
		  if(*(af.eosect+j)==i)
		  {
			pcode1=vp.vflag.ev.code;
			pcode2=vp.vflag.ev.sectioncode;
			/*pcode3=vp.vflag.mv.pagetitle;*/
			vp.vflag.ev.code=1;
			vp.vflag.ev.sectioncode=1;
			/*vp.vflag.mv.pagetitle=1;*/

			if(mode==ONPRINTER)
			{
			  /*drawglobalwire(hdc,vp,af,*(af.elems+i),0,0,0,
													 0,0,0,code,ONSRCANMAX);*/
			  drawglobalwire(hdc,vp,af,*(af.elems+i),0,0,0,
													 255,0,255,code,ONSRCANMAX);
			  /*drawglobalwire(hdc,vp,af,*(af.elems+i),255,0,255,
													 255,0,255,code,ONSRCANMAX);*/
			}
			else
			{
			  drawglobalwire(hdc,vp,af,*(af.elems+i),255,0,255,
													 255,0,255,code,mode);
			}

			vp.vflag.ev.code=pcode1;
			vp.vflag.ev.sectioncode=pcode2;
			/*vp.vflag.mv.pagetitle=pcode3;*/

			break;
		  }
		}
		if(j>=(af.nosect))
		{
		  if(mode==ONPRINTER)
		  {
			drawglobalwire(hdc,vp,af,*(af.elems+i),0,0,0,
												   0,0,0,code,mode);
		  }
		  else
		  {
			drawglobalwire(hdc,vp,af,*(af.elems+i),255,255,255,
												   255,255,255,code,mode);
		  }
		}
	  }
	  /****************/
#endif
      else if(mode==ONPRINTER)
      {
        hsize=vp.dparam.hsize;
        vp.dparam.hsize*=5.0;
        drawglobalwire(hdc,vp,af,*(af.elems+i),0,0,0,
                                               0,0,0,code,mode);
        vp.dparam.hsize=hsize;
      }
      else if(mode==ONPREVIEW)
      {
        hsize=vp.dparam.hsize;
        vp.dparam.hsize*=5.0;
        drawglobalwire(hdc,vp,af,*(af.elems+i),0,0,0,
                                               0,0,0,code,mode);
        vp.dparam.hsize=hsize;
      }
      else
      {
        drawglobalwire(hdc,vp,af,*(af.elems+i),255,255,255,
                       (af.elems+i)->sect->dcolor.r,
                       (af.elems+i)->sect->dcolor.g,
                       (af.elems+i)->sect->dcolor.b,code,mode);
      }
      if(af.melem!=NULL)
      {
        inputelem(af.elems,af.melem,i,&elem);
        drawwirestress(hdc,vp,af,elem,mode);
      }
    }


    else if(((af.elems+i)->sect->dflag==1 &&
       vp.vflag.ev.etype[(af.elems+i)->sect->type]==1)       )
    {
      selected=0;
      for(ii=0;ii<nmultiwire;ii++)
      {
         if(af.elems+i==*(multiwire+ii))
         {
           selected=1;
//           sprintf(str,"MULTIWIRE[%ld]:CODE=%ld",ii,(*(multiwire+ii))->code);
//           errormessage(str);
         }
      }
      if(selected)
      {
        if(mode==ONPRINTER)
        {
          hsize=vp.dparam.hsize;
          vp.dparam.hsize*=5.0;
          drawglobalwire(hdc,vp,af,*(af.elems+i),0,0,0,
                                                 0,0,0,code,mode);
          vp.dparam.hsize=hsize;
        }
        else if(mode==ONPREVIEW)
        {
          hsize=vp.dparam.hsize;
          vp.dparam.hsize*=5.0;
          drawglobalwire(hdc,vp,af,*(af.elems+i),0,0,0,
                                                 0,0,0,code,mode);
          vp.dparam.hsize=hsize;
        }
        else
        {
          drawglobalwire(hdc,vp,af,*(af.elems+i),255,0,255,
                       (af.elems+i)->sect->dcolor.r,
                       (af.elems+i)->sect->dcolor.g,
                       (af.elems+i)->sect->dcolor.b,code,mode);
        }
        if(af.melem!=NULL)
        {
        inputelem(af.elems,af.melem,i,&elem);
        drawwirestress(hdc,vp,af,elem,mode);
        }
      }
      else if(!selected)
      {
        if(mode==ONPRINTER)
        {
          hsize=vp.dparam.hsize;
          vp.dparam.hsize*=5.0;
          drawglobalwireII(hdc,vp,af,*(af.elems+i),200,200,200,
                                                 200,200,200,code,mode);
          vp.dparam.hsize=hsize;
        }
        else if(mode==ONPREVIEW)
        {
          hsize=vp.dparam.hsize;
          vp.dparam.hsize*=5.0;
          drawglobalwireII(hdc,vp,af,*(af.elems+i),200,200,200,
                                                 200,200,200,code,mode);
          vp.dparam.hsize=hsize;
        }
        else
        {
          drawglobalwireII(hdc,vp,af,*(af.elems+i),255,0,255,
                       100,100,100,code,mode);
		}
/*
        if(af.melem!=NULL)
        {
        inputelem(af.elems,af.melem,i,&elem);
		drawwirestress(hdc,vp,af,elem,mode);
        }
*/
      }
    }
  }

  MSG msg;
  for(i=0;i<af.nshell;i++)
  {
	while (PeekMessage(&msg, NULL, 0, 0, PM_REMOVE))
	{
		TranslateMessage(&msg);
		DispatchMessage(&msg);
	}
	 drawglobalshell(hdc,vp,af,*(af.shells+i),255,0,255,255,0,255,code,mode);
  }
  return;
}/*drawarclmframe*/

void savearclmasdxf(FILE *fout,struct viewparam vp,
					struct arclmframe af)
/*BASED ON SAVEORGANASDXF*/
{
  char str[256];
  int i,flag;
  long int code=11587/*2D43*/;
  double dmax[2],dmin[2];
  double jiheight,jipitch;
  struct onode mp,np,mpd,npd,inode1,inode2,dnode1,dnode2;

  fseek(fout,0L,SEEK_SET);

  fprintf(fout,"\"DXF FILE.\"\n");

  /*TEXT SIZE*/
  jiheight=2.0;
  jipitch =2.5;

  /*GET FRAME RANGE FOR DRAWING TEXTS*/
  dmax[0]=0.0; dmin[0]=0.0;
  dmax[1]=0.0; dmin[1]=0.0;

  for(i=0;i<af.nnode;i++)
  {
    if(insiderange(*(af.ninit+i),vp.range))
    {
      if(!nodeprojection(*(af.ninit+i),&np,vp)) return;

      if(i==0) /*INITIAL*/
      {
        dmax[0]=np.d[GX]; dmin[0]=np.d[GX];
        dmax[1]=np.d[GY]; dmin[1]=np.d[GY];
      }

      if(dmax[0]<np.d[GX]) dmax[0]=np.d[GX];
      if(dmax[1]<np.d[GY]) dmax[1]=np.d[GY];
      if(dmin[0]>np.d[GX]) dmin[0]=np.d[GX];
      if(dmin[1]>np.d[GY]) dmin[1]=np.d[GY];
    }
  }
  fprintf(fout,"RANGE\n");
  fprintf(fout,"X1=%.16f\nY1=%.16f\nX2=%.16f\nY2=%.16f\n",
          dmin[0],dmin[1],dmax[0],dmax[1]);
  fprintf(fout,"\n");

  /*UPPER TEXTS*/
  dmax[1]+=jipitch;
  if(vp.vflag.mv.title)
  {
    code++;
    dmax[1]+=jipitch;
    sprintf(str,"%s",(wdraw.childs+1)->title);
    outputtextasdxf(fout,str,code,dmin[0],dmax[1],0.0,jiheight);
  }

  dmax[1]+=jiheight; /*UPPER EDGE OF TEXT*/

  /*LOWER TEXTS*/
  dmin[1]-=jipitch;
  if(vp.vflag.nv.code)
  {
    code++;
    dmin[1]-=jipitch;
    sprintf(str,"節点番号");
    outputtextasdxf(fout,str,code,dmin[0],dmin[1],0.0,jiheight);
  }
  if(vp.vflag.ev.code)
  {
    code++;
    dmin[1]-=jipitch;
    sprintf(str,"部材番号");
    outputtextasdxf(fout,str,code,dmin[0],dmin[1],0.0,jiheight);
  }
  if(vp.vflag.ev.sectioncode)
  {
    code++;
    dmin[1]-=jipitch;
    sprintf(str,"断面番号");
    outputtextasdxf(fout,str,code,dmin[0],dmin[1],0.0,jiheight);
  }
  if(vp.vflag.ev.deformation)
  {
    code++;
    dmin[1]-=jipitch;
    sprintf(str,"変形図");
    outputtextasdxf(fout,str,code,dmin[0],dmin[1],0.0,jiheight);
  }
  if(vp.vflag.nv.disps[0])
  {
    code++;
    dmin[1]-=jipitch;
    sprintf(str,"Ｘ方向変位 [cm]");
    outputtextasdxf(fout,str,code,dmin[0],dmin[1],0.0,jiheight);
  }
  if(vp.vflag.nv.disps[1])
  {
    code++;
    dmin[1]-=jipitch;
    sprintf(str,"Ｙ方向変位 [cm]");
    outputtextasdxf(fout,str,code,dmin[0],dmin[1],0.0,jiheight);
  }
  if(vp.vflag.nv.disps[2])
  {
    code++;
    dmin[1]-=jipitch;
    sprintf(str,"Ｚ方向変位 [cm]");
    outputtextasdxf(fout,str,code,dmin[0],dmin[1],0.0,jiheight);
  }
  if(vp.vflag.ev.axis)
  {
    code++;
    dmin[1]-=jipitch;
    sprintf(str,"部材座標軸");
    outputtextasdxf(fout,str,code,dmin[0],dmin[1],0.0,jiheight);
  }
  for(i=0;i<8;i++)
  {
    if(vp.vflag.ev.stress[i][0])
    {
      code++;
      dmin[1]-=jipitch;
      if(globalunit==1.0)         sprintf(str,"軸力 [tf]");
      else if(globalunit==SIUNIT) sprintf(str,"軸力 [kN]");
      else                        sprintf(str,"軸力");
      outputtextasdxf(fout,str,code,dmin[0],dmin[1],0.0,jiheight);
      break;
    }
  }
  for(i=0;i<8;i++)
  {
    if(vp.vflag.ev.stress[i][1] || vp.vflag.ev.stress[i][2])
    {
      code++;
      dmin[1]-=jipitch;
      if(globalunit==1.0)         sprintf(str,"せん断力 [tf]");
      else if(globalunit==SIUNIT) sprintf(str,"せん断力 [kN]");
      else                        sprintf(str,"せん断力");
      outputtextasdxf(fout,str,code,dmin[0],dmin[1],0.0,jiheight);
      break;
    }
  }
  for(i=0;i<8;i++)
  {
    if(vp.vflag.ev.stress[i][3])
    {
      code++;
      dmin[1]-=jipitch;
      if(globalunit==1.0)         sprintf(str,"ねじりモーメント [tfm]");
      else if(globalunit==SIUNIT) sprintf(str,"ねじりモーメント [kNm]");
      else                        sprintf(str,"ねじりモーメント");
      outputtextasdxf(fout,str,code,dmin[0],dmin[1],0.0,jiheight);
      break;
    }
  }
  for(i=0;i<8;i++)
  {
    if(vp.vflag.ev.stress[i][4] || vp.vflag.ev.stress[i][5])
    {
      code++;
      dmin[1]-=jipitch;
      if(globalunit==1.0)         sprintf(str,"曲げモーメント [tfm]");
      else if(globalunit==SIUNIT) sprintf(str,"曲げモーメント [kNm]");
      else                        sprintf(str,"曲げモーメント");
      outputtextasdxf(fout,str,code,dmin[0],dmin[1],0.0,jiheight);
      break;
    }
  }
  if(vp.vflag.ev.etype[0] ||
     vp.vflag.ev.etype[1] ||
     vp.vflag.ev.etype[2] ||
     vp.vflag.ev.etype[3] ||
     vp.vflag.ev.etype[4] ||
     vp.vflag.ev.etype[5] ||
     vp.vflag.ev.etype[6] ||
     vp.vflag.ev.etype[7])
  {
    flag=0;
    code++;
    dmin[1]-=jipitch;
    sprintf(str,"表示部材：");
    /*if(vp.vflag.ev.etype[0])
    {
      strcat(str,"補助");
      flag=1;
    }*/
    if(vp.vflag.ev.etype[1])
    {
      if(flag) strcat(str,", ");
      strcat(str,"柱");
      flag=1;
    }
    if(vp.vflag.ev.etype[2])
    {
      if(flag) strcat(str,", ");
      strcat(str,"梁");
      flag=1;
    }
    /*if(vp.vflag.ev.etype[3])
    {
      if(flag) strcat(str,", ");
      strcat(str,"小梁");
      flag=1;
    }*/
    if(vp.vflag.ev.etype[4])
    {
      if(flag) strcat(str,", ");
      strcat(str,"ブレース");
      flag=1;
    }
    if(vp.vflag.ev.etype[5])
    {
      if(flag) strcat(str,", ");
      strcat(str,"壁");
      flag=1;
    }
    if(vp.vflag.ev.etype[6])
    {
      if(flag) strcat(str,", ");
      strcat(str,"床");
      flag=1;
    }
    /*if(vp.vflag.ev.etype[7])
    {
      if(flag) strcat(str,", ");
      strcat(str,"張壁");
      flag=1;
    }*/
    outputtextasdxf(fout,str,code,dmin[0],dmin[1],0.0,jiheight);
  }
  dmin[1]-=jipitch;
  if(vp.vflag.mv.inputfile)
  {
    code++;
    dmin[1]-=jipitch;
    sprintf(str,"Input  : %s",(wdraw.childs+1)->inpfile);
    outputtextasdxf(fout,str,code,dmin[0],dmin[1],0.0,jiheight);
  }
  if(vp.vflag.mv.outputfile)
  {
    code++;
    dmin[1]-=jipitch;
    sprintf(str,"Output : %s",(wdraw.childs+1)->otpfile);
    outputtextasdxf(fout,str,code,dmin[0],dmin[1],0.0,jiheight);
  }
  if(vp.vflag.mv.view)
  {
    code++;
    dmin[1]-=jipitch;
    sprintf(str,"Focus : %.3f %.3f %.3f",
            vp.focus.d[GX],vp.focus.d[GY],vp.focus.d[GZ]);
    outputtextasdxf(fout,str,code,dmin[0],dmin[1],0.0,jiheight);
    code++;
    dmin[1]-=jipitch;
    sprintf(str,"Phi=%.3f Theta=%.3f",vp.phi,vp.theta);
    outputtextasdxf(fout,str,code,dmin[0],dmin[1],0.0,jiheight);
    code++;
    dmin[1]-=jipitch;
    sprintf(str,"R=%.3f L=%.3f",vp.r,vp.odv);
    outputtextasdxf(fout,str,code,dmin[0],dmin[1],0.0,jiheight);
    code++;
    dmin[1]-=jipitch;
    sprintf(str,"Dfact=%.3f",vp.dparam.dfact);
    outputtextasdxf(fout,str,code,dmin[0],dmin[1],0.0,jiheight);
    code++;
    dmin[1]-=jipitch;
    sprintf(str,"Mfact=%.3f",vp.dparam.mfact);
    outputtextasdxf(fout,str,code,dmin[0],dmin[1],0.0,jiheight);
  }

  if(vp.vflag.axis==1) /*GLOBAL AXIS.*/
  {
    outputglobalaxisasdxf(fout,&code,jiheight,vp);
  }

  for(i=0;i<af.nelem;i++) /*ELEMENTS.*/
  {
    if((af.elems+i)->sect->dflag==1 &&
       vp.vflag.ev.etype[(af.elems+i)->sect->type]==1)
    {
      initialnode(af.ninit,af.nnode,
                  ((af.elems+i)->node[0]->code),&inode1);
      if(!insiderange(inode1,vp.range)) return;
      if(af.ddisp==NULL || !vp.vflag.ev.deformation)
      {
        dnode1.d[GX]=0.0;
        dnode1.d[GY]=0.0;
        dnode1.d[GZ]=0.0;
      }
      else                                          /*DISPLACEMENT.*/
      {
        dnode1.d[GX]=(af.elems+i)->node[0]->d[GX]-inode1.d[GX];
        dnode1.d[GY]=(af.elems+i)->node[0]->d[GY]-inode1.d[GY];
        dnode1.d[GZ]=(af.elems+i)->node[0]->d[GZ]-inode1.d[GZ];

        dnode1.d[GX]=inode1.d[GX]+(vp.dparam.dfact)*dnode1.d[GX];
        dnode1.d[GY]=inode1.d[GY]+(vp.dparam.dfact)*dnode1.d[GY];
        dnode1.d[GZ]=inode1.d[GZ]+(vp.dparam.dfact)*dnode1.d[GZ];

        nodeprojection(dnode1,&mpd,vp);               /*PROJECTION.*/
      }
      nodeprojection(inode1,&mp,vp);                  /*PROJECTION.*/

      initialnode(af.ninit,af.nnode,
                  ((af.elems+i)->node[1]->code),&inode2);
      if(!insiderange(inode2,vp.range)) return;
      if(af.ddisp==NULL || !vp.vflag.ev.deformation)
      {
        dnode2.d[GX]=0.0;
        dnode2.d[GY]=0.0;
        dnode2.d[GZ]=0.0;
      }
      else                                          /*DISPLACEMENT.*/
      {
        dnode2.d[GX]=(af.elems+i)->node[1]->d[GX]-inode2.d[GX];
        dnode2.d[GY]=(af.elems+i)->node[1]->d[GY]-inode2.d[GY];
        dnode2.d[GZ]=(af.elems+i)->node[1]->d[GZ]-inode2.d[GZ];

        dnode2.d[GX]=inode2.d[GX]+(vp.dparam.dfact)*dnode2.d[GX];
        dnode2.d[GY]=inode2.d[GY]+(vp.dparam.dfact)*dnode2.d[GY];
        dnode2.d[GZ]=inode2.d[GZ]+(vp.dparam.dfact)*dnode2.d[GZ];

        nodeprojection(dnode2,&npd,vp);               /*PROJECTION.*/
      }
      nodeprojection(inode2,&np,vp);                  /*PROJECTION.*/

      /*ELEMENT LINE*/
      if(vp.vflag.ev.deformation && af.ddisp!=NULL)
      {
        code++;

        if((af.elems+i)->sect->code>900) sprintf(str,"None");
        else sprintf(str,"D%d",(af.elems+i)->sect->code);

        outputlineasdxf(fout,code,str,
                        mpd.d[GX],mpd.d[GY],mpd.d[GZ],
                        npd.d[GX],npd.d[GY],npd.d[GZ]);
      }

      code++;
      if((af.elems+i)->sect->code>900) sprintf(str,"None");
      else sprintf(str,"%d",(af.elems+i)->sect->code);
      outputlineasdxf(fout,code,str,
                      mp.d[GX],mp.d[GY],mp.d[GZ],
                      np.d[GX],np.d[GY],np.d[GZ]);

      /*STRESSES*/
      /*if(af.melem!=NULL)
      {
        inputelem(af.elems,af.melem,i,&elem);
        drawwirestress(hdc,vp,af,elem,mode);
      }*/
    }
  }

  if(vp.vflag.nv.code)
  {
    for(i=0;i<af.nnode;i++) /*NODES*/
    {
      if(insiderange(*(af.nodes+i),vp.range))
      {
        nodeprojection(*(af.nodes+i),&np,vp);
        sprintf(str,"%d",(af.nodes+i)->code);
        code++;
        outputtextasdxf(fout,str,code,np.d[GX],
                                      np.d[GY],
                                      np.d[GZ],jiheight);
      }
    }
  }

  fprintf(fout,"\n");
  fprintf(fout,"HANDSEED=%X\n",code+3);
  fprintf(fout,"VPORT   =%X\n",code+2);

  return;
}/*savearclmasdxf*/

void drawrotatingarclm(HDC hdc,struct viewparam vp,
                       struct arclmframe af)
{
  int i;

  drawglobalaxis(hdc,vp,0,0,255);               /*DRAW GLOBAL AXIS.*/

  SetTextColor(hdc,RGB(0,255,255));
  for(i=0;i<af.nnode;i++)
  {
    drawglobalnode(hdc,vp,*(af.ninit+i),NULL);
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
  ppen=(HPEN)SelectObject(hdc,hpen1);

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
/*RETURN:POINTER OF NODE.*/
{
  int i;
  int ix1,iy1;                                /*WINDOW COORDINATION*/
  struct onode node;
  int boxsize=BOXMOUSE;                            /*CURSOR REGION.*/

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

struct onode *selectorgannode(struct viewparam vp,
                              struct organ *org,POINT point)
/*SELECT NODE OF ORGAN FRAME BY CURSOR.*/
/*RETURN:POINTER OF NODE.*/
{
  /*char non[80];*/
  int i;
  int ix1,iy1;                                /*WINDOW COORDINATION*/
  struct onode node;
  int boxsize=BOXMOUSE;                            /*CURSOR REGION.*/

  for(i=0;i<(org->nnode);i++)
  {
    node=*(org->nodes+i);
    if(insiderange(node,vp.range))
    {
      nodeontoscreen(node,&ix1,&iy1,vp);               /*PROJECTION*/
      /*
      sprintf(non,"MOUSE(%ld,%ld),NODE(%d,%d)",
              point.x,point.y,ix1,iy1);
      MessageBox(NULL,non,"Node",MB_OK);
      */
      if((point.x-boxsize)<=ix1 && ix1<=(point.x+boxsize))
      {
        if((point.y-boxsize)<=iy1 && iy1<=(point.y+boxsize))
        {
          return (org->nodes+i);
        }
      }
    }
  }
  return NULL;
}/*selectorgannode*/

int selectelemnode(struct viewparam vp,struct oelem *oe,POINT point)
/*SELECT NODE OF ORGAN ELEMENT BY CURSOR.*/
/*RETURN:POINTER OF NODE.*/
{
  /*char non[80];*/
  int i;
  int ix1,iy1;                                /*WINDOW COORDINATION*/
  struct onode *node;
  int boxsize=BOXMOUSE;                            /*CURSOR REGION.*/

  for(i=0;i<(oe->nnod);i++)
  {
    node=*(oe->nods+i);
    if(insiderange(*node,vp.range))
    {
      nodeontoscreen(*node,&ix1,&iy1,vp);              /*PROJECTION*/
      /*
      sprintf(non,"MOUSE(%ld,%ld),NODE(%d,%d)",
              point.x,point.y,ix1,iy1);
      MessageBox(NULL,non,"Node",MB_OK);
      */
      if((point.x-boxsize)<=ix1 && ix1<=(point.x+boxsize))
      {
        if((point.y-boxsize)<=iy1 && iy1<=(point.y+boxsize))
        {
          return i+1;
        }
      }
    }
  }
  return 0;
}/*selectelemnode*/

struct onode *selectpolycurvenode(struct viewparam vp,
                                  struct polypolycurve *pc,
                                  POINT point)
/*SELECT NODE OF POLYCURVE BY CURSOR.*/
/*RETURN:POINTER OF NODE.*/
{
  /*char non[80];*/
  int i,j,ii;
  int ix1,iy1;                                /*WINDOW COORDINATION*/
  struct onode *node;
  int boxsize=BOXMOUSE;                            /*CURSOR REGION.*/
  struct curve *c;

  for(ii=0;ii<(pc->npcurve);ii++)
  {
    for(i=0;i<((pc->pcurves+ii)->ncurve);i++)
    {
      c=(pc->pcurves+ii)->curves+i;

      for(j=0;j<=3;j++)
      {
        node=NULL;

        if(j==0) node=c->center;
        if(j==1) node=c->dots[0];
        if(j==2) node=c->dots[1];
        if(j==3) node=c->dots[2];

        if(node!=NULL)
        {
          nodeontoscreen(*node,&ix1,&iy1,vp);         /*PROJECTION*/
          if((point.x-boxsize)<=ix1 && ix1<=(point.x+boxsize))
          {
            if((point.y-boxsize)<=iy1 && iy1<=(point.y+boxsize))
            {
              return node;
            }
          }
        }
      }
    }
  }
  return NULL;
}/*selectpolycurvenode*/

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
  int boxsize=BOXMOUSE;                            /*CURSOR REGION.*/

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

struct oelem *selectorganelement(struct viewparam vp,
                                 struct organ *org,
                                 POINT point,
                                 struct line **selectedline,
                                 struct obans **selectedban)
/*SELECT ELEMENT BY CURSOR.*/
/*RETURN:POINTER OF ELEMENT,SELECTED LINE,SELECTED BAN.*/
{
  int i,j,k;
  int ix1,iy1,ix2,iy2;                       /*WINDOW COORDINATION.*/
  int ihead,itail,inside;
  int Nx,Ny;
  long int loff1,loff2;
  struct onode *bnods,nmouse;
  struct obans ban,*cban;
  struct line *sline;

  *selectedline=NULL;
  *selectedban =NULL;

  sline=(struct line *)malloc(sizeof(struct line));
  if(sline==NULL)
  {
    MessageBox(NULL,"Buffer Null.","SelectOrganElement",MB_OK);
    return NULL;
  }

  bnods=(struct onode *)malloc((org->nnode)*sizeof(struct onode));
  if(bnods==NULL)
  {
    MessageBox(NULL,"Buffer Null.","SelectOrganElement",MB_OK);
    return NULL;
  }

  for(i=0;i<(org->nnode);i++)
  {
    if(!nodeontoscreen(*(org->nodes+i),&Nx,&Ny,vp)) /*PROJECTION*/
    {
      MessageBox(NULL,"Projection Failed.","SelectOrganElement",
                 MB_OK);
      return NULL;
    }

    /*(bnods+i)->d[GX]=(double)( Nx-vp.Xo);
    (bnods+i)->d[GY]=(double)(-Ny+vp.Yo);*/
    (bnods+i)->d[GX]=(double)Nx;
    (bnods+i)->d[GY]=(double)Ny;
    (bnods+i)->d[GZ]=0.0;
  }

  nmouse.d[GX]=(double)(point.x);
  nmouse.d[GY]=(double)(point.y);
  nmouse.d[GZ]=0.0;

  for(i=0;i<(org->nelem);i++)
  {
    /*WHETHER ELEMENT INSIDE RANGE*/
    inside=1;
    if((org->elems+i)->sect->dflag==1 &&
       vp.vflag.ev.etype[(org->elems+i)->type]==1)
    {
      for(j=0;j<(org->elems+i)->nnod;j++)
      {
        if(!insiderange(**((org->elems+i)->nods+j),vp.range))
        {
          inside=0; /*OUTSIDE*/
        }
      }
    }
    else inside=0; /*FLAGS OFF*/

    /*ONLY FOR INSIDE RANGE AND "DRAW FLAG","ELEMENT TYPE FLAG" ON*/
    if(inside)
    {
      if(((org->elems+i)->nban)<=0 && ((org->elems+i)->nnod)==2)
      {
        loff1=(*((org->elems+i)->nods+0))->loff;
        ix1=(int)((bnods+loff1)->d[GX]);
        iy1=(int)((bnods+loff1)->d[GY]);

        loff2=(*((org->elems+i)->nods+1))->loff;
        ix2=(int)((bnods+loff2)->d[GX]);
        iy2=(int)((bnods+loff2)->d[GY]);

        if(intersectlinemouse(ix1,iy1,ix2,iy2,point))
        {
          sline->ends[0]=*(*((org->elems+i)->nods+0));
          sline->ends[1]=*(*((org->elems+i)->nods+1));
          *selectedline=sline; /*FOUND WIRE ELEMENT.*/

          free(bnods);
          return (org->elems+i);
        }
      }
      else if(((org->elems+i)->nban)>0)
      {
        for(j=0;j<((org->elems+i)->nban);j++)
        {
          cban=(org->elems+i)->bans+j;

          ban.nods=(struct onode **)malloc((cban->nnod)
                                           *sizeof(struct onode *));
          if(ban.nods==NULL)
          {
            MessageBox(NULL,"Buffer Null.","SelectOrganElement",MB_OK);
            return NULL;
          }

          ban.nnod=cban->nnod;

          for(k=0;k<(cban->nnod);k++)
          {
            loff1=(*(cban->nods+k))->loff;
            (*(ban.nods+k))=(bnods+loff1);
          }

          for(k=0;k<(cban->nnod);k++)
          {
            ihead=k;
            if(k==(cban->nnod)-1) itail=0;
            else                  itail=k+1;

            loff1=(*(cban->nods+ihead))->loff;
            ix1=(int)((bnods+loff1)->d[GX]);
            iy1=(int)((bnods+loff1)->d[GY]);

            loff2=(*(cban->nods+itail))->loff;
            ix2=(int)((bnods+loff2)->d[GX]);
            iy2=(int)((bnods+loff2)->d[GY]);

            if(intersectlinemouse(ix1,iy1,ix2,iy2,point))
            {
              sline->ends[0]=*(*(cban->nods+ihead));
              sline->ends[1]=*(*(cban->nods+itail));

              *selectedline=sline; /*FOUND LINE OF BAN.*/
              *selectedban=(org->elems+i)->bans+j;

              free(ban.nods);
              free(bnods);
              return (org->elems+i);
            }
          }

          if(insideban(ban,nmouse))
          {
            *selectedban=(org->elems+i)->bans+j; /*FOUND BAN.*/

            free(ban.nods);
            free(bnods);
            return (org->elems+i);
          }

          free(ban.nods);
        }
      }
    }
  }
  free(bnods);

  return NULL;
}/*selectorganelement*/

void initializeorganization(struct organ *org)
/*INITIALIZE ORGAN.*/
{
  org->code=0;
  org->loff=0;

  org->nprop=0;
  org->nsect=0;
  org->nnode=0;
  org->nelem=0;

  org->nodes=NULL;
  org->confs=NULL;
  org->props=NULL;
  org->sects=NULL;
  org->elems=NULL;
  org->loads=NULL;

  org->ai.T1=0.0;
  org->ai.Tc=0.0;
  org->ai.Cox=0.0;
  org->ai.Coy=0.0;
  org->ai.Z =0.0;

  return;
}/*initializeorganization*/

int inputorganization(FILE *fin,
                      struct organ *org,
                      struct viewparam *vp)
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

  char **data,str[256],string[256];
  int i,j,k,ii,nstr,pstr,ident,ncv;
  long int code;
  long int pnode=(-1),pelem=(-1),pprop=(-1),psect=(-1);
  long int psfig,peban,ppcurve,pcurve;

  fseek(fin,0L,SEEK_SET);

  org->code=0;
  org->loff=0;

  org->nnode=0; /*INITIALIZATION.*/
  org->nelem=0;
  org->nprop=0;
  org->nsect=0;

  org->ntext=0;
  org->texts=NULL;

//  org->ai.Co=0.0;
  org->ai.Cox=0.0;          //ujioka for BASE
  org->ai.Coy=0.0;
  org->ai.Z =0.0;
  org->ai.T1=0.0;
  org->ai.Tc=0.0;

//add by fukushima 20140618
  org->ai.nfloor=19;
  org->ai.lbound=(double *)malloc((org->ai.nfloor+1)*sizeof(double));
  org->ai.nfloor=-1;
///////////////////////////

  fgets(str,256,fin); /*APPELATION.*/

  while(1)
  {
    data=fgetsbrk(fin,&nstr);
    if(nstr==0)
    {
      /*ENDING PROCESS.*/
      for(i=0;i<(org->nsect);i++)
      {
        for(j=0;j<((org->sects+i)->ppc.npcurve);j++)
        {
          for(ii=0;ii<(org->nprop);ii++)
          {
            if(((org->sects+i)->ppc.pcurves+j)->prop.code
               ==((org->props+ii)->code))
            {
              ((org->sects+i)->ppc.pcurves+j)->prop
              =*(org->props+ii);

              ((org->sects+i)->ppc.pcurves+j)->prop.name
              =(char *)malloc((strlen((org->props+ii)->name)+1)
                              *sizeof(char));
              strcpy(((org->sects+i)->ppc.pcurves+j)->prop.name,
					 (org->props+ii)->name);
              break;
            }
          }

          for(k=0;k<(((org->sects+i)->ppc.pcurves+j)->ncurve);k++)
          {
            if((((org->sects+i)->ppc.pcurves+j)->curves+k)->type
               ==CTYPE_CIRCLE)
            {
              createcircledata(((org->sects+i)->ppc.pcurves+j)
                               ->curves+k);
            }
          }
        }
      }

      return 0;
    }

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

        if(vp!=NULL && nstr-pstr>=4 && !strcmp(*(data+pstr),"TEXT"))
        {
          pstr++;

          org->texts=(struct snode *)
                     realloc(org->texts,
                             (org->ntext+1)*sizeof(struct snode));
          /*COORDINATION*/
          (org->texts+(org->ntext))->n.d[0]=strtod(*(data+pstr+0),NULL);
          (org->texts+(org->ntext))->n.d[1]=strtod(*(data+pstr+1),NULL);
		  (org->texts+(org->ntext))->n.d[2]=0.0;

          pstr+=2;
          strcpy(str,*(data+pstr));
          while(nstr-pstr>1)
          {
            pstr++;
            strcat(str," "); /*TEXTS DIVIDED BY ONLY ONE SPACE.*/
            strcat(str,*(data+pstr));
          }
/*MessageBox(NULL,str,"Input",MB_OK);*/

		  strcpy((org->texts+(org->ntext))->str,str);

		  org->ntext++;
		  ident++;
		}

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
/*
          sprintf(string,"NSECT=%d",org->nsect);
		  MessageBox(NULL,string,"INPUTORGANIZATION",MB_OK);
*/
		  org->sects=(struct osect *)malloc(org->nsect
                                            *sizeof(struct osect));
          pstr++;
          ident++;
        }

        /*VIEW DATA*/
        if(vp!=NULL && nstr-pstr>=2 && !strcmp(*(data+pstr),"GFACT"))
        {
          pstr++;
          vp->gfactor=strtod(*(data+pstr),NULL);
          pstr++;
          ident++;
        }
        if(vp!=NULL && nstr-pstr>=4 && !strcmp(*(data+pstr),"FOCUS"))
        {
          pstr++;
          vp->focus.d[0]=strtod(*(data+pstr+0),NULL);
          vp->focus.d[1]=strtod(*(data+pstr+1),NULL);
          vp->focus.d[2]=strtod(*(data+pstr+2),NULL);
          pstr+=3;
          ident++;
        }
        if(vp!=NULL && nstr-pstr>=3 && !strcmp(*(data+pstr),"ANGLE"))
        {
          pstr++;
          vp->phi  =strtod(*(data+pstr+0),NULL);
          vp->theta=strtod(*(data+pstr+1),NULL);
          pstr+=2;
          ident++;
        }
        if(vp!=NULL && nstr-pstr>=3 && !strcmp(*(data+pstr),"DISTS"))
        {
          pstr++;
          vp->r  =strtod(*(data+pstr+0),NULL);
          vp->odv=strtod(*(data+pstr+1),NULL);
          pstr+=2;
          ident++;
        }

        /*EARTHQUAKE DATA*/
        if(nstr-pstr>=2 && !strcmp(*(data+pstr),"BASE"))
        {
/*
          pstr++;
          org->ai.Co=strtod(*(data+pstr),NULL);
          pstr++;
          ident++;
*/
		  pstr++;

          org->ai.Cox=strtod(*(data+pstr),NULL);
		  pstr++;

		  if (*(data+pstr)!=NULL)org->ai.Coy=strtod(*(data+pstr),NULL);
          else org->ai.Coy=org->ai.Cox;
		  pstr++;

		  ident++;

        }
		if(nstr-pstr>=2 && !strcmp(*(data+pstr),"LOCATE"))
        {
		  pstr++;
          org->ai.Z=strtod(*(data+pstr),NULL);
          pstr++;
          ident++;
        }
        if(nstr-pstr>=2 && !strcmp(*(data+pstr),"TFACT"))
        {
          pstr++;
          org->ai.T1=strtod(*(data+pstr),NULL);
          pstr++;
          ident++;
        }
        if(nstr-pstr>=2 && !strcmp(*(data+pstr),"GPERIOD"))
        {
          pstr++;
          org->ai.Tc=strtod(*(data+pstr),NULL);
          pstr++;
          ident++;
        }
//add by fukushima 20140618
		if(nstr-pstr>=2 && !strcmp(*(data+pstr),"NFLOOR"))
		{
		  pstr++;
		  org->ai.nfloor=(int)strtol(*(data+pstr),NULL,10);
		  pstr++;
		  ident++;
		}
		if(nstr-pstr>=2+org->ai.nfloor && !strcmp(*(data+pstr),"HEIGHT"))
		{
		  pstr++;
		  for (i=0; i<=org->ai.nfloor; i++)
		  {
			*(org->ai.lbound+ i)= strtod(*(data+pstr),NULL);
			pstr++;
		  }
		  ident++;
		}
////////////////////////////

        /*PROPERTIES:CODE,HIJU,E,POI.*/
        if(nstr-pstr>=2 && (org->nprop)-pprop>1 &&
           !strcmp(*(data+pstr),"PROP"))
        {
          pstr++;
          pprop++;
          (org->props+pprop)->code=strtol(*(data+pstr),NULL,10);
          (org->props+pprop)->loff=pprop;
          /*(org->props+pprop)->name=NULL;*/
		  (org->props+pprop)->name=(char *)malloc(1*sizeof(char));
		  strcpy((org->props+pprop)->name,"\0");
          (org->props+pprop)->r=0;
		  (org->props+pprop)->g=0;
          (org->props+pprop)->b=0;
          pstr++;
          ident++;
        }
        if(nstr-pstr>=2 && (org->nprop)-pprop>0 &&
           !strcmp(*(data+pstr),"PNAME"))
        {
          pstr++;
          (org->props+pprop)->name=(char *)
                                   realloc((org->props+pprop)->name,
                                           (strlen(*(data+pstr))+1)
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
        if(nstr-pstr>=4 && (org->nprop)-pprop>0 &&
           !strcmp(*(data+pstr),"PCOLOR"))
        {
          pstr++;
          (org->props+pprop)->r
          =(int)strtol(*(data+pstr+0),NULL,10);
          (org->props+pprop)->g
          =(int)strtol(*(data+pstr+1),NULL,10);
          (org->props+pprop)->b
          =(int)strtol(*(data+pstr+2),NULL,10);
          pstr+=3;
          ident++;
        }

        /*SECTIONS:CODE,FIGS,SURFACE.*/
        if(nstr-pstr>=2 && (org->nsect)-psect>1 &&
           !strcmp(*(data+pstr),"SECT"))
        {
          pstr++;
          psect++;
          (org->sects+psect)->code=strtol(*(data+pstr),NULL,10);
          (org->sects+psect)->ocode=(org->sects+psect)->code;
          (org->sects+psect)->loff=psect;

          /*(org->sects+psect)->name=NULL;*/
          (org->sects+psect)->name=(char *)malloc(1*sizeof(char));
          strcpy((org->sects+psect)->name,"\0");
          (org->sects+psect)->role=ROLENULL;
          (org->sects+psect)->type=TYPENULL;
          (org->sects+psect)->nfig=0;
          (org->sects+psect)->dflag=1;
          (org->sects+psect)->dcolor.r=255;
          (org->sects+psect)->dcolor.g=255;
          (org->sects+psect)->dcolor.b=255;

          (org->sects+psect)->E   =0.0;
          (org->sects+psect)->poi =0.0;
          (org->sects+psect)->area=0.0;
          (org->sects+psect)->Ixx =0.0;
          (org->sects+psect)->Iyy =0.0;
          (org->sects+psect)->Jzz =0.0;
          (org->sects+psect)->hiju[0]=0.0;
          (org->sects+psect)->hiju[1]=0.0;
          (org->sects+psect)->hiju[2]=0.0;

          (org->sects+psect)->ppc.npcurve=0;

          (org->sects+psect)->lload[0]=0.0;
          (org->sects+psect)->lload[1]=0.0;
          (org->sects+psect)->lload[2]=0.0;
          (org->sects+psect)->perpl[0]=0.0;
          (org->sects+psect)->perpl[1]=0.0;
          (org->sects+psect)->perpl[2]=0.0;

          pstr++;
          ident++;
        }
        if(nstr-pstr>=2 && (org->nsect)-psect>0 &&
           !strcmp(*(data+pstr),"SNAME"))
        {
          pstr++;
          (org->sects+psect)->name=(char *)
                                   realloc((org->sects+psect)->name,
                                           (strlen(*(data+pstr))+1)
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

          ((org->sects+psect)->figs+psfig)->thick=0.0;
          ((org->sects+psect)->figs+psfig)->area =0.0;

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

          ((org->sects+psect)->figs+psfig)->area=0.0;   //ujioka
          ((org->sects+psect)->figs+psfig)->Ixx =0.0;
          ((org->sects+psect)->figs+psfig)->Iyy =0.0;
          ((org->sects+psect)->figs+psfig)->Jzz =0.0;

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

        if(nstr-pstr>=2 && (org->nsect)-psect>0 &&
           !strcmp(*(data+pstr),"SROLE"))
        {
          pstr++;
          if(!strcmp(*(data+pstr),"HOJO"))
          {
            (org->sects+psect)->role=ROLEHOJO;
          }
          pstr++;
          ident++;
        }

        if(nstr-pstr>=4 && (org->nsect)-psect>0 &&
           !strcmp(*(data+pstr),"LLOAD"))
        {
          pstr++;
          (org->sects+psect)->lload[0]=strtod(*(data+pstr+0),NULL);
          (org->sects+psect)->lload[1]=strtod(*(data+pstr+1),NULL);
          (org->sects+psect)->lload[2]=strtod(*(data+pstr+2),NULL);
          pstr+=3;
          ident++;
        }

        if(nstr-pstr>=4 && (org->nsect)-psect>0 &&
           !strcmp(*(data+pstr),"PERPL"))     /*PERPENDICULAR LOAD*/
        {
          pstr++;
          (org->sects+psect)->perpl[0]=strtod(*(data+pstr+0),NULL);
          (org->sects+psect)->perpl[1]=strtod(*(data+pstr+1),NULL);
          (org->sects+psect)->perpl[2]=strtod(*(data+pstr+2),NULL);
          pstr+=3;
          ident++;
        }

        if(nstr-pstr>=4 && (org->nsect)-psect>0 &&
           !strcmp(*(data+pstr),"COLOR"))
        {
          pstr++;
          (org->sects+psect)->dcolor.r
          =(int)strtol(*(data+pstr+0),NULL,10);
          (org->sects+psect)->dcolor.g
          =(int)strtol(*(data+pstr+1),NULL,10);
          (org->sects+psect)->dcolor.b
          =(int)strtol(*(data+pstr+2),NULL,10);
          pstr+=3;
          ident++;
        }

        if(nstr-pstr>=2 && (org->nsect)-psect>0 &&
           !strcmp(*(data+pstr),"NPCURVE"))
        {
          pstr++;
          (org->sects+psect)->ppc.npcurve
          =(int)strtol(*(data+pstr+0),NULL,10);

          (org->sects+psect)->ppc.pcurves
          =(struct polycurve *)
           malloc((org->sects+psect)->ppc.npcurve
                  *sizeof(struct polycurve));

          ppcurve=(-1);
          pstr++;
          ident++;
        }
        if(nstr-pstr>=2 && (org->nsect)-psect>0 &&
           !strcmp(*(data+pstr),"PCURVE"))
        {
          pstr++;
          ppcurve++;
          ((org->sects+psect)->ppc.pcurves+ppcurve)->loff
          =ppcurve;
          /*((org->sects+psect)->ppc.pcurves+ppcurve)->code
          =(int)strtol(*(data+pstr+0),NULL,10);*/
          pstr++;
          ident++;
        }
        if(nstr-pstr>=2 && (org->nsect)-psect>0 &&
           !strcmp(*(data+pstr),"PPROP"))
        {
          pstr++;

          ((org->sects+psect)->ppc.pcurves+ppcurve)->prop.code
          =(int)strtol(*(data+pstr+0),NULL,10);

          pstr++;
          ident++;
        }
        if(nstr-pstr>=2 && (org->nsect)-psect>0 &&
           !strcmp(*(data+pstr),"NCURVE"))
        {
          pstr++;

          ncv=(int)strtol(*(data+pstr+0),NULL,10);

          ((org->sects+psect)->ppc.pcurves+ppcurve)->ncurve=ncv;
          ((org->sects+psect)->ppc.pcurves+ppcurve)->curves
          =malloccurves(ncv);

          pcurve=(-1);
          pstr++;
          ident++;
        }
        if(nstr-pstr>=1 && (org->nsect)-psect>0 &&
           !strcmp(*(data+pstr),"PLINE"))
        {
          pcurve++;

          (((org->sects+psect)->ppc.pcurves+ppcurve)
          ->curves+pcurve)->loff=pcurve;
          (((org->sects+psect)->ppc.pcurves+ppcurve)
          ->curves+pcurve)->type=CTYPE_LINE;

          pstr++;
          ident++;
        }
        if(nstr-pstr>=1 && (org->nsect)-psect>0 &&
           !strcmp(*(data+pstr),"PCIRCLE"))
        {
          pcurve++;

          (((org->sects+psect)->ppc.pcurves+ppcurve)
          ->curves+pcurve)->loff=pcurve;
          (((org->sects+psect)->ppc.pcurves+ppcurve)
          ->curves+pcurve)->type=CTYPE_CIRCLE;
          (((org->sects+psect)->ppc.pcurves+ppcurve)
          ->curves+pcurve)->radius[0]=0.0;

          (((org->sects+psect)->ppc.pcurves+ppcurve)
          ->curves+pcurve)->center
          =(struct onode *)malloc(sizeof(struct onode));

          pstr++;
          ident++;
        }
        if(nstr-pstr>=2 && (org->nsect)-psect>0 &&
           !strcmp(*(data+pstr),"HUGO"))
        {
          pstr++;

          (((org->sects+psect)->ppc.pcurves+ppcurve)
          ->curves+pcurve)->hugo
          =(int)strtol(*(data+pstr+0),NULL,10);

          pstr++;
          ident++;
        }
        if(nstr-pstr>=3 && (org->nsect)-psect>0 &&
           !strcmp(*(data+pstr),"HEAD"))
        {
          pstr++;

          (((org->sects+psect)->ppc.pcurves+ppcurve)
          ->curves+pcurve)->dots[0]
          =(struct onode *)malloc(sizeof(struct onode));
          (((org->sects+psect)->ppc.pcurves+ppcurve)
          ->curves+pcurve)->dots[0]->d[GX]
          =strtod(*(data+pstr+0),NULL);
          (((org->sects+psect)->ppc.pcurves+ppcurve)
          ->curves+pcurve)->dots[0]->d[GY]
          =strtod(*(data+pstr+1),NULL);

          pstr+=2;
          ident++;
        }
        if(nstr-pstr>=3 && (org->nsect)-psect>0 &&
           !strcmp(*(data+pstr),"TAIL"))
        {
          pstr++;

          (((org->sects+psect)->ppc.pcurves+ppcurve)
          ->curves+pcurve)->dots[1]
          =(struct onode *)malloc(sizeof(struct onode));
          (((org->sects+psect)->ppc.pcurves+ppcurve)
          ->curves+pcurve)->dots[1]->d[GX]
          =strtod(*(data+pstr+0),NULL);
          (((org->sects+psect)->ppc.pcurves+ppcurve)
          ->curves+pcurve)->dots[1]->d[GY]
          =strtod(*(data+pstr+1),NULL);

          pstr+=2;
          ident++;
        }
        if(nstr-pstr>=3 && (org->nsect)-psect>0 &&
           !strcmp(*(data+pstr),"CENTER"))
        {
          pstr++;

          (((org->sects+psect)->ppc.pcurves+ppcurve)
          ->curves+pcurve)->center->d[GX]
          =strtod(*(data+pstr+0),NULL);
          (((org->sects+psect)->ppc.pcurves+ppcurve)
          ->curves+pcurve)->center->d[GY]
          =strtod(*(data+pstr+1),NULL);

          pstr+=2;
          ident++;
        }
        if(nstr-pstr>=2 && (org->nsect)-psect>0 &&
           !strcmp(*(data+pstr),"RADIUS"))
        {
          pstr++;

          (((org->sects+psect)->ppc.pcurves+ppcurve)
          ->curves+pcurve)->radius[0]
          =strtod(*(data+pstr+0),NULL);

          pstr++;
          ident++;
        }
        if(nstr-pstr>=2 && (org->nsect)-psect>0 &&
           !strcmp(*(data+pstr),"ANGLE1"))
        {
          pstr++;

          (((org->sects+psect)->ppc.pcurves+ppcurve)
          ->curves+pcurve)->angle[0]
          =strtod(*(data+pstr+0),NULL);

          pstr++;
          ident++;
        }
        if(nstr-pstr>=2 && (org->nsect)-psect>0 &&
           !strcmp(*(data+pstr),"ANGLE2"))
        {
          pstr++;

          (((org->sects+psect)->ppc.pcurves+ppcurve)
          ->curves+pcurve)->angle[1]
          =strtod(*(data+pstr+0),NULL);

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

		  /*CHECK SAME NODE*/
          for(i=0;i<pnode;i++)
          {
			if((org->nodes+pnode)->d[0]==(org->nodes+i)->d[0] &&
               (org->nodes+pnode)->d[1]==(org->nodes+i)->d[1] &&
               (org->nodes+pnode)->d[2]==(org->nodes+i)->d[2])
            {
			  sprintf(str,"NODE SAME %ld=%ld.",
					  (org->nodes+pnode)->code,(org->nodes+i)->code);
              MessageBox(NULL,str,"Input",MB_OK);
            }
          }

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

          (org->elems+pelem)->lface[0]=0.0;
          (org->elems+pelem)->lface[1]=0.0;
          (org->elems+pelem)->hface[0]=0.0;
          (org->elems+pelem)->hface[1]=0.0;
          (org->elems+pelem)->wrect[0]=0.0;
          (org->elems+pelem)->wrect[1]=0.0;

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

          if((org->elems+pelem)->sect->role==ROLEHOJO)
          {
            (org->elems+pelem)->role=ROLEHOJO;
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
          if(!strcmp(*(data+pstr),"COLUMN"))
          {
            (org->elems+pelem)->type=COLUMN;
          }
          else if(!strcmp(*(data+pstr),"GIRDER"))
          {
            (org->elems+pelem)->type=GIRDER;
          }
          else if(!strcmp(*(data+pstr),"BEAM"))
          {
            (org->elems+pelem)->type=BEAM;
          }
          else if(!strcmp(*(data+pstr),"BRACE"))
          {
            (org->elems+pelem)->type=BRACE;
          }
          else if(!strcmp(*(data+pstr),"WALL"))
          {
            (org->elems+pelem)->type=WALL;
            /*(org->elems+pelem)->color.line
            =setrgbcolor(150,100,255);*/
          }
          else if(!strcmp(*(data+pstr),"SLAB"))
          {
            (org->elems+pelem)->type=SLAB;
            /*(org->elems+pelem)->color.line
            =setrgbcolor(255,0,255);*/
          }
/*
          else if(!strcmp(*(data+pstr),"CURTAIN"))
          {
            (org->elems+pelem)->type=CURTAIN;
          }
*/
          else if(!strcmp(*(data+pstr),"NONLINEARBEAM")) /*UJIOKA FOR PEZETTINO*/
          {
            (org->elems+pelem)->type=NONLINEARBEAM;
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
            (org->elems+pelem)->role=ROLENULL;
          }
          pstr++;
          ident++;
        }
        if(nstr-pstr>=3 && (org->nelem)-pelem>0 &&
           !strcmp(*(data+pstr),"LFACE")) /*FACE LENGTH OF WALL.*/
        {
          pstr++;
          (org->elems+pelem)->lface[0]=strtod(*(data+pstr+0),NULL);
          (org->elems+pelem)->lface[1]=strtod(*(data+pstr+1),NULL);
          pstr+=2;
          ident++;
        }
        if(nstr-pstr>=3 && (org->nelem)-pelem>0 &&
           !strcmp(*(data+pstr),"HFACE")) /*FACE HEIGHT OF WALL.*/
        {
          pstr++;
          (org->elems+pelem)->hface[0]=strtod(*(data+pstr+0),NULL);
          (org->elems+pelem)->hface[1]=strtod(*(data+pstr+1),NULL);
          pstr+=2;
          ident++;
        }
        if(nstr-pstr>=3 && (org->nelem)-pelem>0 &&
           !strcmp(*(data+pstr),"WRECT")) /*WINDOW RECT OF WALL.*/
        {
          pstr++;
          (org->elems+pelem)->wrect[0]=strtod(*(data+pstr+0),NULL);
          (org->elems+pelem)->wrect[1]=strtod(*(data+pstr+1),NULL);
          pstr+=2;
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
          for(i=0;i<6*((org->elems+pelem)->nnod);i++)
          {
            *((org->elems+pelem)->bonds+i)=0;
          }

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

          /*CHECK SAME NODE*/
          for(i=0;i<(org->elems+pelem)->nnod-1;i++)
          {
            for(j=i+1;j<(org->elems+pelem)->nnod;j++)
            {
              if(*((org->elems+pelem)->nods+i)
               ==*((org->elems+pelem)->nods+j))
              {
                sprintf(str,"SAME NODE %ld=%ld IN ELEMENT %ld.",
                        (*((org->elems+pelem)->nods+i))->code,
                        (*((org->elems+pelem)->nods+j))->code,
                        (org->elems+pelem)->code);
				MessageBox(NULL,str,"Input",MB_OK);

                /*DELETE SAME NODE*/
                /*UNDER CONSTRUCTION.*/
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

int saveorganization(FILE *fout,
                     struct organ *org,
                     struct viewparam *vp)
/*SAVE ORGAN INTO OUTPUT FILE.*/
{
  int i,j,k;
  long int nnode,nelem,nprop,nsect;
  double eps=1.0E-12;

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

  /*EARTHQUAKE DATA*/
  fprintf(fout,"BASE    %.3f %.3f\n",org->ai.Cox,org->ai.Coy);
  fprintf(fout,"LOCATE  %.3f\n",org->ai.Z);
  fprintf(fout,"TFACT   %.3f\n",org->ai.T1);
  fprintf(fout,"GPERIOD %.3f\n",org->ai.Tc);

//add by fukushima 20140618
  if (org->ai.nfloor > 0)
  {
	fprintf(fout,"NFLOOR %d\n",org->ai.nfloor);
	fprintf(fout,"HEIGHT");
	for (i=0; i<=org->ai.nfloor; i++)
	{
	  fprintf(fout," %.1f",*(org->ai.lbound+i));
	}
	fprintf(fout,"\n");
  }
///////////////////////////

  fprintf(fout,"\n");

  /*VIEW DATA*/
  fprintf(fout,"GFACT %.1f\n",vp->gfactor);
  fprintf(fout,"FOCUS %.1f %.1f %.1f\n",vp->focus.d[0],
                                        vp->focus.d[1],
                                        vp->focus.d[2]);
  fprintf(fout,"ANGLE %.1f %.1f\n",vp->phi,vp->theta);
  fprintf(fout,"DISTS %.1f %.1f\n",vp->r,vp->odv);
  fprintf(fout,"\n");

  for(i=0;i<nprop;i++)
  {
    fprintf(fout,"PROP %3ld ",(org->props+i)->code);

    /*if((org->props+i)->name!=NULL)*/
    if(strlen((org->props+i)->name)>0)
    {
      fprintf(fout,"PNAME %s\n",(org->props+i)->name);
    }
    else
    {
      fprintf(fout,"PNAME \"Noname\"\n");
    }

    fprintf(fout,"         HIJU %13.8f\n",(org->props+i)->hiju);
    fprintf(fout,"         E    %13.3f\n",(org->props+i)->E);
    fprintf(fout,"         POI  %13.5f\n",(org->props+i)->poi);
    fprintf(fout,"         PCOLOR %3d %3d %3d\n",(org->props+i)->r,
                                                 (org->props+i)->g,
                                                 (org->props+i)->b);
  }
  fprintf(fout,"\n");

  for(i=0;i<nsect;i++)
  {
    fprintf(fout,"SECT %3ld ",(org->sects+i)->code);

    /*if((org->sects+i)->name!=NULL)*/
    if(strlen((org->sects+i)->name)>0)
    {
      fprintf(fout,"SNAME %s\n",(org->sects+i)->name);
    }
    else
    {
      fprintf(fout,"SNAME \"Noname\"\n");
    }

    if((org->sects+i)->role==ROLEHOJO)
    {
      fprintf(fout,"         SROLE HOJO\n");
    }
    else
    {
      fprintf(fout,"         NFIG %d\n",(org->sects+i)->nfig);
      for(j=0;j<(org->sects+i)->nfig;j++)
      {
        fprintf(fout,"         FIG %3d ",
                ((org->sects+i)->figs+j)->code);
        fprintf(fout,"FPROP %d\n",
                ((org->sects+i)->figs+j)->prop->code);

        if(((org->sects+i)->figs+j)->area<-eps ||
           ((org->sects+i)->figs+j)->area> eps ||
           ((org->sects+i)->figs+j)->Ixx <-eps ||
           ((org->sects+i)->figs+j)->Ixx > eps ||
           ((org->sects+i)->figs+j)->Iyy <-eps ||
           ((org->sects+i)->figs+j)->Iyy > eps ||
           ((org->sects+i)->figs+j)->Jzz <-eps ||
           ((org->sects+i)->figs+j)->Jzz > eps)
        {
          fprintf(fout,"                 ");
          fprintf(fout,"AREA  %.4f\n",((org->sects+i)->figs+j)->area);
          fprintf(fout,"                 ");
          fprintf(fout,"IXX   %.8f\n",((org->sects+i)->figs+j)->Ixx);
          fprintf(fout,"                 ");
          fprintf(fout,"IYY   %.8f\n",((org->sects+i)->figs+j)->Iyy);
          fprintf(fout,"                 ");
          fprintf(fout,"VEN   %.8f\n",((org->sects+i)->figs+j)->Jzz);
        }
        if(((org->sects+i)->figs+j)->thick<-eps ||
           ((org->sects+i)->figs+j)->thick> eps)
        {
          fprintf(fout,"                 ");
          fprintf(fout,"THICK %.5f\n",((org->sects+i)->figs+j)->thick);
        }
      }

      if((org->sects+i)->exp!=0.0)
      {
        fprintf(fout,"         EXP %.3f\n",(org->sects+i)->exp);

        fprintf(fout,"         ");
        fprintf(fout,"NZMAX %12.6f ", (org->sects+i)->fmax[0]);
        fprintf(fout,"NZMIN %12.6f\n",(org->sects+i)->fmin[0]);
        fprintf(fout,"         ");
        fprintf(fout,"QXMAX %12.6f ", (org->sects+i)->fmax[1]);
        fprintf(fout,"QXMIN %12.6f\n",(org->sects+i)->fmin[1]);
        fprintf(fout,"         ");
        fprintf(fout,"QYMAX %12.6f ", (org->sects+i)->fmax[2]);
        fprintf(fout,"QYMIN %12.6f\n",(org->sects+i)->fmin[2]);
        fprintf(fout,"         ");
        fprintf(fout,"MZMAX %12.6f ", (org->sects+i)->fmax[3]);
        fprintf(fout,"MZMIN %12.6f\n",(org->sects+i)->fmin[3]);
        fprintf(fout,"         ");
        fprintf(fout,"MXMAX %12.6f ", (org->sects+i)->fmax[4]);
        fprintf(fout,"MXMIN %12.6f\n",(org->sects+i)->fmin[4]);
        fprintf(fout,"         ");
        fprintf(fout,"MYMAX %12.6f ", (org->sects+i)->fmax[5]);
        fprintf(fout,"MYMIN %12.6f\n",(org->sects+i)->fmin[5]);
      }
    }

    if((org->sects+i)->lload[0]!=0.0 ||
       (org->sects+i)->lload[1]!=0.0 ||
       (org->sects+i)->lload[2]!=0.0)
    {
      fprintf(fout,"         LLOAD");
      fprintf(fout," %.3f %.3f %.3f\n",(org->sects+i)->lload[0],
                                       (org->sects+i)->lload[1],
                                       (org->sects+i)->lload[2]);
    }

	if((org->sects+i)->perpl[0]!=0.0 ||
       (org->sects+i)->perpl[1]!=0.0 ||
       (org->sects+i)->perpl[2]!=0.0)
    {
      fprintf(fout,"         PERPL");
      fprintf(fout," %.3f %.3f %.3f\n",(org->sects+i)->perpl[0],
                                       (org->sects+i)->perpl[1],
                                       (org->sects+i)->perpl[2]);
    }

    if((org->sects+i)->dcolor.r!=255 ||
       (org->sects+i)->dcolor.g!=255 ||
       (org->sects+i)->dcolor.b!=255)
    {
	  fprintf(fout,"         COLOR");
      fprintf(fout," %d %d %d\n",(org->sects+i)->dcolor.r,
                                 (org->sects+i)->dcolor.g,
                                 (org->sects+i)->dcolor.b);
    }
  }
  fprintf(fout,"\n");

  for(i=0;i<nnode;i++)
  {
	fprintf(fout,"NODE %4ld  ",(org->nodes+i)->code);
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

    fprintf(fout,"VCON  ");
    if((org->confs+6*i+0)->value==0.0) fprintf(fout," 0.0  ");
    else fprintf(fout,"%9.6f ",(org->confs+6*i+0)->value);
    if((org->confs+6*i+1)->value==0.0) fprintf(fout," 0.0  ");
    else fprintf(fout,"%9.6f ",(org->confs+6*i+1)->value);
    if((org->confs+6*i+2)->value==0.0) fprintf(fout," 0.0  ");
    else fprintf(fout,"%9.6f ",(org->confs+6*i+2)->value);
    if((org->confs+6*i+3)->value==0.0) fprintf(fout," 0.0  ");
    else fprintf(fout,"%9.6f ",(org->confs+6*i+3)->value);
    if((org->confs+6*i+4)->value==0.0) fprintf(fout," 0.0  ");
    else fprintf(fout,"%9.6f ",(org->confs+6*i+4)->value);
    if((org->confs+6*i+5)->value==0.0) fprintf(fout," 0.0");
    else fprintf(fout,"%9.6f ",(org->confs+6*i+5)->value);
    fprintf(fout,"\n");
  }
  fprintf(fout,"\n");

  for(i=0;i<nelem;i++)
  {
    fprintf(fout,"ELEM %5ld ",(org->elems+i)->code);
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
      fprintf(fout,"                     ");
      fprintf(fout,"EBANS %d ",(org->elems+i)->nban);
      for(j=0;j<(org->elems+i)->nban;j++)
      {
        if(j>=1) fprintf(fout,"                             ");

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

      fprintf(fout,"           CANG %.5f\n",(org->elems+i)->cangle);


/*FOR SHIZUOKA MODIFY*/
/*
if((org->elems+i)->sect->code==502)
fprintf(fout,"           CANG %.5f\n",-(org->elems+i)->cangle);
*/
/*
if((org->elems+i)->sect->code==401 ||
   (org->elems+i)->sect->code==402 ||
   (org->elems+i)->sect->code==403)
fprintf(fout,"           CANG %.5f\n",0.0);
else fprintf(fout,"           CANG %.5f\n",(org->elems+i)->cangle);
*/
      fprintf(fout,"           CMQ ");
      for(j=0;j<2;j++)
      {
        for(k=0;k<6;k++)
        {
          /*fprintf(fout," %.1f",(org->elems+i)->initial[j][k]);*/
fprintf(fout," %.1f",0.0);
        }
        if(j<1) fprintf(fout," ");
      }
      fprintf(fout,"\n");
    }

    if((org->elems+i)->type==COLUMN)
    {
      fprintf(fout,"           TYPE COLUMN\n");
    }
    if((org->elems+i)->type==GIRDER)
    {
      fprintf(fout,"           TYPE GIRDER\n");
    }
    if((org->elems+i)->type==BEAM)
    {
      fprintf(fout,"           TYPE BEAM\n");
    }
    if((org->elems+i)->type==WALL)
    {
      fprintf(fout,"           TYPE WALL\n");
    }
    if((org->elems+i)->type==SLAB)
    {
      fprintf(fout,"           TYPE SLAB\n");
    }
    if((org->elems+i)->type==BRACE)
    {
      fprintf(fout,"           TYPE BRACE\n");
    }

    if((org->elems+i)->lface[0]!=0.0 ||
       (org->elems+i)->lface[1]!=0.0)
    {
      fprintf(fout,"           LFACE %.3f %.3f\n",
              (org->elems+i)->lface[0],(org->elems+i)->lface[1]);
    }
    if((org->elems+i)->hface[0]!=0.0 ||
       (org->elems+i)->hface[1]!=0.0)
    {
      fprintf(fout,"           HFACE %.3f %.3f\n",
              (org->elems+i)->hface[0],(org->elems+i)->hface[1]);
    }
    if((org->elems+i)->wrect[0]!=0.0 ||
       (org->elems+i)->wrect[1]!=0.0)
    {
      fprintf(fout,"           WRECT %.3f %.3f\n",
              (org->elems+i)->wrect[0],(org->elems+i)->wrect[1]);
    }

    if((org->elems+i)->role==ROLEWEIGHT)
    {
      fprintf(fout,"           ROLE W\n");
    }
    if((org->elems+i)->role==ROLERIGID)
    {
      fprintf(fout,"           ROLE R\n");
    }
    if((org->elems+i)->role==ROLESTRESS)
    {
      fprintf(fout,"           ROLE S\n");
    }
  }

  return 1;
}/*saveorganization*/

int checkorganization(struct organ *org)
/*CHECK ORGAN INPUT DATA*/
{
  FILE *fout;
  int i,j;
  long int nnode,nelem,nprop,nsect;
  double eps1=1.0E-01;

  fout=fopen("hogcheck.txt","w");

  /*CHECK NODE DISTANCE*/
  for(i=0;i<(org->nnode-1);i++)
  {
    for(j=i+1;j<org->nnode;j++)
    {
      if(fabs((org->nodes+i)->d[0]-(org->nodes+j)->d[0])<eps1 &&
         fabs((org->nodes+i)->d[1]-(org->nodes+j)->d[1])<eps1 &&
         fabs((org->nodes+i)->d[2]-(org->nodes+j)->d[2])<eps1)
      {
        fprintf(fout,"DISTANCE SHORT NODE %ld NODE %d\n",
                (org->nodes+i)->code,(org->nodes+j)->code);
      }
    }
  }

  /*CHECK SAME ELEMENT*/
  for(i=0;i<(org->nelem-1);i++)
  {
  }

  fclose(fout);
  return 1;
}/*checkorganization*/

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
  free(org->loads); /*LOADS.*/
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

int saveorganasdxf(FILE *fout,
                   struct organ go,
                   struct viewparam vp)
/*SAVE ORGAN AS DXF FILE.*/
{
  int i,j;
  char str[256];
  int flag;
  long int code=11587/*2D43*/ /*11544=2D18*/;
  double dmax[2],dmin[2];
  double jiheight,jipitch;
  /*double hiju,orate,radius,hsize;*/
  struct onode mp,np;

  fseek(fout,0L,SEEK_SET);

  fprintf(fout,"\"DXF FILE.\"\n");

  /*TEXT SIZE*/
  jiheight=2.0;
  jipitch =2.5;

  /*GET FRAME RANGE FOR DRAWING TEXTS*/
  dmax[0]=0.0; dmin[0]=0.0;
  dmax[1]=0.0; dmin[1]=0.0;

  for(i=0;i<go.nnode;i++)
  {
    if(insiderange(*(go.nodes+i),vp.range))
    {
      if(!nodeprojection(*(go.nodes+i),&np,vp)) return 0;

      if(i==0) /*INITIAL*/
      {
        dmax[0]=np.d[GX]; dmin[0]=np.d[GX];
        dmax[1]=np.d[GY]; dmin[1]=np.d[GY];
      }

      if(dmax[0]<np.d[GX]) dmax[0]=np.d[GX];
      if(dmax[1]<np.d[GY]) dmax[1]=np.d[GY];
      if(dmin[0]>np.d[GX]) dmin[0]=np.d[GX];
      if(dmin[1]>np.d[GY]) dmin[1]=np.d[GY];
    }
  }

  fprintf(fout,"RANGE\n");
  fprintf(fout,"X1=%.16f\nY1=%.16f\nX2=%.16f\nY2=%.16f\n",
          dmin[0],dmin[1],dmax[0],dmax[1]);
  fprintf(fout,"\n");

  /*UPPER TEXTS*/
  dmax[1]+=jipitch;
  if(vp.vflag.mv.title)
  {
    code++;
    dmax[1]+=jipitch;
    sprintf(str,"%s",(wdraw.childs+1)->title);
    outputtextasdxf(fout,str,code,dmin[0],dmax[1],0.0,jiheight);
  }

  dmax[1]+=jiheight; /*UPPER EDGE OF TEXT*/

  /*LOWER TEXTS*/
  dmin[1]-=jipitch;
  if(vp.vflag.nv.code)
  {
    code++;
    dmin[1]-=jipitch;
    sprintf(str,"節点番号");
    outputtextasdxf(fout,str,code,dmin[0],dmin[1],0.0,jiheight);
  }
  if(vp.vflag.ev.code)
  {
    code++;
    dmin[1]-=jipitch;
    sprintf(str,"部材番号");
    outputtextasdxf(fout,str,code,dmin[0],dmin[1],0.0,jiheight);
  }
  if(vp.vflag.ev.sectioncode)
  {
    code++;
    dmin[1]-=jipitch;
    sprintf(str,"断面番号");
    outputtextasdxf(fout,str,code,dmin[0],dmin[1],0.0,jiheight);
  }
  if(vp.vflag.nv.mcircle)
  {
    code++;
    dmin[1]-=jipitch;
    sprintf(str,"節点重量図");
    outputtextasdxf(fout,str,code,dmin[0],dmin[1],0.0,jiheight);
  }
  if(vp.vflag.nv.mvalue)
  {
    code++;
    dmin[1]-=jipitch;
    if(globalunit==1.0)         sprintf(str,"節点重量値 [tf]");
    else if(globalunit==SIUNIT) sprintf(str,"節点重量値 [kN]");
    else                        sprintf(str,"節点重量値");
    outputtextasdxf(fout,str,code,dmin[0],dmin[1],0.0,jiheight);
  }
  if(vp.vflag.ev.axis)
  {
    code++;
    dmin[1]-=jipitch;
    sprintf(str,"部材座標軸");
    outputtextasdxf(fout,str,code,dmin[0],dmin[1],0.0,jiheight);
  }
  if(vp.vflag.ev.cmqline)
  {
    code++;
    dmin[1]-=jipitch;
    sprintf(str,"分割線表示");
    outputtextasdxf(fout,str,code,dmin[0],dmin[1],0.0,jiheight);
  }
  if(vp.vflag.ev.etype[0] ||
     vp.vflag.ev.etype[1] ||
     vp.vflag.ev.etype[2] ||
     vp.vflag.ev.etype[3] ||
     vp.vflag.ev.etype[4] ||
     vp.vflag.ev.etype[5] ||
     vp.vflag.ev.etype[6] ||
     vp.vflag.ev.etype[7])
  {
    flag=0;
    code++;
    dmin[1]-=jipitch;
    sprintf(str,"表示部材：");
    /*if(vp.vflag.ev.etype[0])
    {
      strcat(str,"補助");
      flag=1;
    }*/
    if(vp.vflag.ev.etype[1])
    {
      if(flag) strcat(str,", ");
      strcat(str,"柱");
      flag=1;
    }
    if(vp.vflag.ev.etype[2])
    {
      if(flag) strcat(str,", ");
      strcat(str,"梁");
      flag=1;
    }
    /*if(vp.vflag.ev.etype[3])
    {
      if(flag) strcat(str,", ");
      strcat(str,"小梁");
      flag=1;
    }*/
    if(vp.vflag.ev.etype[4])
    {
      if(flag) strcat(str,", ");
      strcat(str,"ブレース");
      flag=1;
    }
    if(vp.vflag.ev.etype[5])
    {
      if(flag) strcat(str,", ");
      strcat(str,"壁");
      flag=1;
    }
    if(vp.vflag.ev.etype[6])
    {
      if(flag) strcat(str,", ");
      strcat(str,"床");
      flag=1;
    }
    /*if(vp.vflag.ev.etype[7])
    {
      if(flag) strcat(str,", ");
      strcat(str,"張壁");
      flag=1;
    }*/
    outputtextasdxf(fout,str,code,dmin[0],dmin[1],0.0,jiheight);
  }
  dmin[1]-=jipitch;
  if(vp.vflag.mv.inputfile)
  {
    code++;
    dmin[1]-=jipitch;
    sprintf(str,"Input  : %s",(wdraw.childs+1)->inpfile);
    outputtextasdxf(fout,str,code,dmin[0],dmin[1],0.0,jiheight);
  }
  if(vp.vflag.mv.outputfile)
  {
    code++;
    dmin[1]-=jipitch;
    sprintf(str,"Output : %s",(wdraw.childs+1)->otpfile);
    outputtextasdxf(fout,str,code,dmin[0],dmin[1],0.0,jiheight);
  }
  if(vp.vflag.mv.view)
  {
    code++;
    dmin[1]-=jipitch;
    sprintf(str,"Focus : %.3f %.3f %.3f",
            vp.focus.d[GX],vp.focus.d[GY],vp.focus.d[GZ]);
    outputtextasdxf(fout,str,code,dmin[0],dmin[1],0.0,jiheight);
    code++;
    dmin[1]-=jipitch;
    sprintf(str,"Phi=%.3f Theta=%.3f",vp.phi,vp.theta);
    outputtextasdxf(fout,str,code,dmin[0],dmin[1],0.0,jiheight);
    code++;
    dmin[1]-=jipitch;
    sprintf(str,"R=%.3f L=%.3f",vp.r,vp.odv);
    outputtextasdxf(fout,str,code,dmin[0],dmin[1],0.0,jiheight);
    code++;
    dmin[1]-=jipitch;
    sprintf(str,"Dfact=%.3f",vp.dparam.dfact);
    outputtextasdxf(fout,str,code,dmin[0],dmin[1],0.0,jiheight);
    code++;
    dmin[1]-=jipitch;
    sprintf(str,"Mfact=%.3f",vp.dparam.mfact);
    outputtextasdxf(fout,str,code,dmin[0],dmin[1],0.0,jiheight);
  }

  if(vp.vflag.axis==1) /*GLOBAL AXIS.*/
  {
    outputglobalaxisasdxf(fout,&code,jiheight,vp);
  }

  for(i=0;i<go.nelem;i++) /*ELEMENTS.*/
  {
    if((go.elems+i)->sect->dflag==1 &&
       vp.vflag.ev.etype[(go.elems+i)->type]==1)
    {
      if((go.elems+i)->nnod==2)
      {
        nodeprojection(**((go.elems+i)->nods+0),&mp,vp);
        nodeprojection(**((go.elems+i)->nods+1),&np,vp);

        code++;
        sprintf(str,"%d",(go.elems+i)->sect->code);
        outputlineasdxf(fout,code,str,
                        mp.d[GX],mp.d[GY],mp.d[GZ],
                        np.d[GX],np.d[GY],np.d[GZ]);
      }
      else if((go.elems+i)->nnod>2 && (go.elems+i)->nban>0)
      {
        for(j=0;j<(go.elems+i)->nban;j++)
        {
          code++;
          outputglobalbanasdxf(fout,vp,*((go.elems+i)->bans+j),code);
        }
      }
    }
  }

  if(vp.vflag.nv.code)
  {
    for(i=0;i<go.nnode;i++) /*NODES*/
    {
      if(insiderange(*(go.nodes+i),vp.range))
      {
        nodeprojection(*(go.nodes+i),&np,vp);
        sprintf(str,"%d",(go.nodes+i)->code);
        code++;
        outputtextasdxf(fout,str,code,np.d[GX],
                                      np.d[GY],
                                      np.d[GZ],jiheight);
      }
    }
  }

  fprintf(fout,"\n");
  fprintf(fout,"HANDSEED=%X\n",code+3);
  fprintf(fout,"VPORT   =%X\n",code+2);

  return 1;
}/*saveorganasdxf*/

void outputtextasdxf(FILE *f,char str[],long int code,
                     double x,double y,double z,double size)
{
  fprintf(f,"TEXT\n");
  fprintf(f,"  5\n");
  fprintf(f,"%X\n",code);
  fprintf(f,"330\n");
  fprintf(f,"A\n");
  fprintf(f,"100\n");
  fprintf(f,"AcDbEntity\n");
  fprintf(f,"  8\n");
  fprintf(f,"Moji\n"); /*LAYER*/
  fprintf(f,"100\n");
  fprintf(f,"AcDbText\n");
  fprintf(f," 10\n");
  fprintf(f,"%.16f\n",x);
  fprintf(f," 20\n");
  fprintf(f,"%.16f\n",y);
  fprintf(f," 30\n");
  fprintf(f,"%.16f\n",z);
  fprintf(f," 40\n");
  fprintf(f,"%.1f\n",size);
  fprintf(f,"  1\n");
  fprintf(f,"%s\n",str);
  fprintf(f,"100\n");
  fprintf(f,"AcDbText\n");
  fprintf(f,"  0\n");

  return;
}/*outputtextasdxf*/

void outputlineasdxf(FILE *f,long int code,
                     char layer[],
                     double x1,double y1,double z1,
                     double x2,double y2,double z2)
{
  fprintf(f,"LINE\n");
  fprintf(f,"  5\n");
  fprintf(f,"%X\n",code);
  fprintf(f,"330\n");
  fprintf(f,"A\n");
  fprintf(f,"100\n");
  fprintf(f,"AcDbEntity\n");
  fprintf(f,"  8\n");
  fprintf(f,"%s\n",layer); /*LAYER*/
  fprintf(f,"100\n");
  fprintf(f,"AcDbLine\n");
  fprintf(f," 10\n");
  fprintf(f,"%.16f\n",x1);
  fprintf(f," 20\n");
  fprintf(f,"%.16f\n",y1);
  fprintf(f," 30\n");
  fprintf(f,"%.16f\n",z1);
  fprintf(f," 11\n");
  fprintf(f,"%.16f\n",x2);
  fprintf(f," 21\n");
  fprintf(f,"%.16f\n",y2);
  fprintf(f," 31\n");
  fprintf(f,"%.16f\n",z2);
  fprintf(f,"  0\n");

  return;
}/*outputlineasdxf*/

void outputglobalarrowasdxf(FILE *f,
                            struct onode gn1,
                            struct onode gn2,
                            double asize,
                            long int *code)
/*DRAW ARROW AS DXF.*/
{
  double dx,dy,arrowlength,arrowsize,arrowtan;

  (*code)++;
  outputlineasdxf(f,*code,"None",
                  gn1.d[GX],gn1.d[GY],gn1.d[GZ],
                  gn2.d[GX],gn2.d[GY],gn2.d[GZ]);   /*AXIS.*/

  dx=gn2.d[GX]-gn1.d[GX];
  dy=gn2.d[GY]-gn1.d[GY];
  arrowlength=sqrt(dx*dx+dy*dy);

  if(asize<=0.0) arrowsize=0.2*arrowlength;
  else           arrowsize=asize;

  arrowtan=0.1*arrowsize;

  if(arrowlength!=0.0)  /*ARROW.*/
  {
    gn1.d[GX]
   =gn2.d[GX]
    -((arrowsize/arrowlength)*dx+(dy/arrowlength)*arrowtan);
    gn1.d[GY]
   =gn2.d[GY]
    -((arrowsize/arrowlength)*dy-(dx/arrowlength)*arrowtan);

    (*code)++;
    outputlineasdxf(f,*code,"None",
                    gn1.d[GX],gn1.d[GY],gn1.d[GZ],
                    gn2.d[GX],gn2.d[GY],gn2.d[GZ]);

    gn1.d[GX]
   =gn2.d[GX]
    -((arrowsize/arrowlength)*dx-(dy/arrowlength)*arrowtan);
    gn1.d[GY]
   =gn2.d[GY]
    -((arrowsize/arrowlength)*dy+(dx/arrowlength)*arrowtan);

    (*code)++;
    outputlineasdxf(f,*code,"None",
                    gn1.d[GX],gn1.d[GY],gn1.d[GZ],
                    gn2.d[GX],gn2.d[GY],gn2.d[GZ]);
  }
  return;
}/*outputglobalarrowasdxf*/

void outputglobalaxisasdxf(FILE *f,long int *code,double size,
                           struct viewparam vp)
/*DRAW GLOBAL AXIS AS DXF.*/
{
  double arrowlength;
  struct onode hn,on,hp,op;

  arrowlength=vp.dparam.gaxis;

  on=setcoord(0.0,0.0,0.0);
  nodeprojection(on,&op,vp);
  (*code)++;
  outputtextasdxf(f,"O",*code,op.d[GX],
                              op.d[GY],
                              op.d[GZ],size);              /*ORIGIN*/

  hn=setcoord(arrowlength,0.0,0.0);
  nodeprojection(hn,&hp,vp);
  (*code)++;
  outputtextasdxf(f,"X",*code,hp.d[GX],
                              hp.d[GY],
                              hp.d[GZ],size);              /*TEXT X*/
  outputglobalarrowasdxf(f,op,hp,0.0,code);               /*ARROW X*/

  hn=setcoord(0.0,arrowlength,0.0);
  nodeprojection(hn,&hp,vp);
  (*code)++;
  outputtextasdxf(f,"Y",*code,hp.d[GX],
                              hp.d[GY],
                              hp.d[GZ],size);              /*TEXT Y*/
  outputglobalarrowasdxf(f,op,hp,0.0,code);               /*ARROW Y*/

  hn=setcoord(0.0,0.0,arrowlength);
  nodeprojection(hn,&hp,vp);
  (*code)++;
  outputtextasdxf(f,"Z",*code,hp.d[GX],
                              hp.d[GY],
                              hp.d[GZ],size);              /*TEXT Z*/
  outputglobalarrowasdxf(f,op,hp,0.0,code);               /*ARROW Z*/

  return;
}/*outputglobalaxisasdxf*/

void outputglobalbanasdxf(FILE *f,struct viewparam vp,
                          struct obans gb,long int code)
/*DRAW BAN AS DXF.*/
{
  int i;
  struct onode gn,np;

  fprintf(f,"LWPOLYLINE\n");
  fprintf(f,"  5\n");
  fprintf(f,"%X\n",code);
  fprintf(f,"330\n");
  fprintf(f,"A\n");
  fprintf(f,"100\n");
  fprintf(f,"AcDbEntity\n");
  fprintf(f,"  8\n");
  fprintf(f,"None\n"); /*LAYER*/
  fprintf(f,"100\n");
  fprintf(f,"AcDbPolyline\n");
  fprintf(f," 90\n");
  fprintf(f,"        %d\n",gb.nnod+1);
  fprintf(f," 70\n");
  fprintf(f,"   128\n");
  fprintf(f," 43\n");
  fprintf(f,"0.0\n");

  for(i=0;i<gb.nnod;i++)
  {
    gn=*(*(gb.nods+i));
    nodeprojection(gn,&np,vp);

    fprintf(f," 10\n");
    fprintf(f,"%.16f\n",np.d[GX]);
    fprintf(f," 20\n");
    fprintf(f,"%.16f\n",np.d[GY]);
  }

  /*HEAD DOT*/
  gn=*(*(gb.nods+0));
  nodeprojection(gn,&np,vp);

  fprintf(f," 10\n");
  fprintf(f,"%.16f\n",np.d[GX]);
  fprintf(f," 20\n");
  fprintf(f,"%.16f\n",np.d[GY]);

  fprintf(f,"0\n");

  return;
}/*outputglobalbanasdxf*/

int extractarclmfromorgan(struct organ *org,
                          struct arclmframe *az,
						  struct arclmframe *ax,
						  struct arclmframe *ay)
/*EXTRACT ARCLM FROM ORGAN.*/
{
  char str[256],non[80],prj[80];
  int i,j,k,l,m,n,flag,find;
  int count,nsect,pnsect,nelem,pnelem,ielem,nnode,pnnode,inode,nreact;
  long int oloff,aloff,loff;
  double distance,n1,n2,eps=1.0E-02; /*[m]*/
  struct osect *s,*psects;
  struct onode *on,*nodes,*pnodes,*mnodes,*ncross1,*ncross2,nearest,nc1,nc2;
  struct oelem *elems,*pelems,*melems,*penods;
  struct oconf *confs,*mconfs,*pconfs;
  struct oload *loads;
  struct line l1,l2;
  FILE *fout;

  /*PROJECT NAME*/
  GetDlgItemText((wmenu.childs+2)->hwnd,ID_SECTIONFILE,non,80);
  if(!strncmp(non,"vene",4)) strcpy(prj,"vene");
  if(!strncmp(non,"inax",4)) strcpy(prj,"inax");
  if(!strncmp(non,"stack",5)) strcpy(prj,"stack");
  if(!strncmp(non,"LunarBase",9))    strcpy(prj,"lunarbase");

  /*CHECK FACTORS*/
  GetDlgItemText((wmenu.childs+4)->hwnd,IDVS_PFACT,non,80);
  sprintf(str,"Period Factor=%s Base Shear=",non);
  GetDlgItemText((wmenu.childs+4)->hwnd,IDVS_BASESHEAR,non,80);
  strcat(str,non);

  if(globalmessageflag==1)
  {
	if(MessageBox(NULL,str,"Extract",MB_OKCANCEL)==IDCANCEL)
	{
	  return 0;
	}
  }

#if 1
if(!strncmp(prj,"inax",4) || !strncmp(prj,"stack",5))
{
  pnnode=org->nnode;
  /*pnsect=org->nsect;*/
  /*psects=org->sects;*/
  pnelem=org->nelem;
  pelems=org->elems;
  penods=(struct oelem *)malloc(pnelem*sizeof(struct oelem));
  pnodes=(struct onode *)malloc(pnnode*sizeof(struct onode));
  pconfs=(struct oconf *)malloc(6*pnnode*sizeof(struct oconf));
  for(k=0;k<pnnode;k++)
  {
	*(pnodes+k)=*(org->nodes+k);
	for(l=0;l<6;l++)
	{
	  (pconfs+6*k+l)->iconf=(org->confs+6*k+l)->iconf;
	  (pconfs+6*k+l)->value=(org->confs+6*k+l)->value;
	}
  }
  for(k=0;k<pnelem;k++)
  {
	*(penods+k)=*(pelems+k);
	/*(penods+k)->nnod=(pelems+k)->nnod;*/
	/*(penods+k)->nban=(pelems+k)->nban;*/
	(penods+k)->nods
	=(struct onode **)
	 malloc((penods+k)->nnod*sizeof(struct onode *));
	(penods+k)->bonds
	=(signed char *)
	 malloc(6*((pelems+k)->nnod)*sizeof(signed char));
	(penods+k)->bans
	=(struct obans *)
	 malloc((penods+k)->nban*sizeof(struct obans));

	for(l=0;l<(pelems+k)->nnod;l++)
	{
	  *((penods+k)->nods+l)=pnodes+(*((pelems+k)->nods+l))->loff;

	  *((penods+k)->bonds+6*l+0)=*((pelems+k)->bonds+6*l+0);
	  *((penods+k)->bonds+6*l+1)=*((pelems+k)->bonds+6*l+1);
	  *((penods+k)->bonds+6*l+2)=*((pelems+k)->bonds+6*l+2);
	  *((penods+k)->bonds+6*l+3)=*((pelems+k)->bonds+6*l+3);
	  *((penods+k)->bonds+6*l+4)=*((pelems+k)->bonds+6*l+4);
	  *((penods+k)->bonds+6*l+5)=*((pelems+k)->bonds+6*l+5);
	}
	for(l=0;l<(pelems+k)->nban;l++)
	{
	  ((penods+k)->bans+l)->nnod=((pelems+k)->bans+l)->nnod;
	  ((penods+k)->bans+l)->nods
	  =(struct onode **)
	   malloc(((penods+k)->bans+l)->nnod*sizeof(struct onode *));

	  for(m=0;m<((pelems+k)->bans+l)->nnod;m++)
	  {
		*(((penods+k)->bans+l)->nods+m)
		=pnodes+(*(((pelems+k)->bans+l)->nods+m))->loff;
	  }
	}
  }
}
#endif

#if 1
if(!strncmp(prj,"inax",4) || !strncmp(prj,"stack",5)) /*DEVIDE ELEM*/
{
  nnode=org->nnode;
  nodes=org->nodes;
  confs=org->confs;
  nelem=org->nelem;
  elems=org->elems;
  inode=(org->nodes+nnode-1)->code;
  ielem=(org->elems+nelem-1)->code;

sprintf(str,"NNODE=%d NELEM=%d\nELEM%d:HEAD=%d,TAIL=%d\nELEM%d:HEAD=%d,TAIL=%d",
		nnode,nelem,
		(elems+0)->code,
		(*((elems+0)->nods+0))->code,
		(*((elems+0)->nods+1))->code,
		(elems+nelem-1)->code,
		(*((elems+nelem-1)->nods+0))->code,
		(*((elems+nelem-1)->nods+1))->code);
MessageBox(NULL,str,"STACK",MB_OK);

		  /*FIND NODE ON CREATED ELEMENT*/
		  for(j=0;j<nelem;j++)
		  {
/*
sprintf(str,"ELEM %d",(elems+j)->code);
MessageBox(NULL,str,"CROSS",MB_OK);
*/
			if((elems+j)->nnod==2)
			{
			  for(i=0;i<nnode;i++)
			  {
				/*GET NODE COORDINATION*/
				on=*((elems+j)->nods+0);
				l1.ends[0].d[GX]=on->d[GX];
				l1.ends[0].d[GY]=on->d[GY];
				l1.ends[0].d[GZ]=on->d[GZ];

				on=*((elems+j)->nods+1);
				l1.ends[1].d[GX]=on->d[GX];
				l1.ends[1].d[GY]=on->d[GY];
				l1.ends[1].d[GZ]=on->d[GZ];

				distance=distancedotline(*(nodes+i),l1,&nearest,&flag);
				if(flag==1 && fabs(distance)<eps)
				{
				  /*DIVIDE ELEMENT*/
				  nelem++;
				  ielem++;

				  /*elems=(struct oelem *)
						realloc(elems,(nelem+2*count)*sizeof(struct oelem));*/
				  melems=elems;
				  elems=(struct oelem *)
						malloc((nelem/*+2*count*/)*sizeof(struct oelem));
				  for(k=0;k<nelem-1;k++) *(elems+k)=*(melems+k);

				  (elems+nelem-1)->code=ielem;
				  (elems+nelem-1)->sect=(melems+j)->sect;
				  (elems+nelem-1)->type=(melems+j)->type;
				  (elems+nelem-1)->cangle=(melems+j)->cangle;
				  (elems+nelem-1)->role=(melems+j)->role;
				  (elems+nelem-1)->nnod=2;
				  (elems+nelem-1)->nban=0;
				  (elems+nelem-1)->nods
				  =(struct onode **)
				   malloc((elems+nelem-1)->nnod*sizeof(struct onode *));

				  *((elems+nelem-1)->nods+0)=(nodes+i);
				  *((elems+nelem-1)->nods+1)=*((melems+j)->nods+1);
				  *((elems+j)->nods+1)      =(nodes+i);

/*
sprintf(str,"ELEM=%d NELEM=%d",(elems+j)->code,nelem);
MessageBox(NULL,str,"NODE ON ELEM",MB_OK);
*/
				  (elems+nelem-1)->bonds
				  =(signed char *)
				   malloc(6*((elems+nelem-1)->nnod)*sizeof(signed char));

				  *((elems+nelem-1)->bonds+6*0+0)=0;
				  *((elems+nelem-1)->bonds+6*0+1)=0;
				  *((elems+nelem-1)->bonds+6*0+2)=0;
				  *((elems+nelem-1)->bonds+6*0+3)=0;
				  *((elems+nelem-1)->bonds+6*0+4)=0;
				  *((elems+nelem-1)->bonds+6*0+5)=0;

				  *((elems+nelem-1)->bonds+6*1+0)=*((elems+j)->bonds+6*1+0);
				  *((elems+nelem-1)->bonds+6*1+1)=*((elems+j)->bonds+6*1+1);
				  *((elems+nelem-1)->bonds+6*1+2)=*((elems+j)->bonds+6*1+2);
				  *((elems+nelem-1)->bonds+6*1+3)=*((elems+j)->bonds+6*1+3);
				  *((elems+nelem-1)->bonds+6*1+4)=*((elems+j)->bonds+6*1+4);
				  *((elems+nelem-1)->bonds+6*1+5)=*((elems+j)->bonds+6*1+5);

				  *((elems+j)->bonds+6*1+0)=0;
				  *((elems+j)->bonds+6*1+1)=0;
				  *((elems+j)->bonds+6*1+2)=0;
				  *((elems+j)->bonds+6*1+3)=0;
				  *((elems+j)->bonds+6*1+4)=0;
				  *((elems+j)->bonds+6*1+5)=0;
				}
			  }
			}
		  }

/*MessageBox(NULL,"Pass1.2","Stack",MB_OK);*/

		  /*FIND CROSSING POINT ON CREATED ELEMENT*/
		  for(j=0;j<nelem;j++)
		  {
/*
sprintf(str,"NNODE=%d NELEM=%d\nELEM%d:HEAD=%d,TAIL=%d\nELEM%d:HEAD=%d,TAIL=%d\nELEM%d:HEAD=%d,TAIL=%d",
		nnode,nelem,
		(elems+0)->code,
		(*((elems+0)->nods+0))->code,
		(*((elems+0)->nods+1))->code,
		(elems+j)->code,
		(*((elems+j)->nods+0))->code,
		(*((elems+j)->nods+1))->code,
		(elems+nelem-1)->code,
		(*((elems+nelem-1)->nods+0))->code,
		(*((elems+nelem-1)->nods+1))->code);
MessageBox(NULL,str,"STACK",MB_OK);
*/
			if((elems+j)->nnod==2)
			{
			  for(i=j;i<nelem;i++)
			  {
				if((elems+i)->nnod==2)
				{
				  /*GET NODE COORDINATION*/
				  on=*((elems+j)->nods+0);
				  l1.ends[0].d[GX]=on->d[GX];
				  l1.ends[0].d[GY]=on->d[GY];
				  l1.ends[0].d[GZ]=on->d[GZ];

				  on=*((elems+j)->nods+1);
				  l1.ends[1].d[GX]=on->d[GX];
				  l1.ends[1].d[GY]=on->d[GY];
				  l1.ends[1].d[GZ]=on->d[GZ];

				  /*GET NODE COORDINATION*/
				  on=*((elems+i)->nods+0);
				  l2.ends[0].d[GX]=on->d[GX];
				  l2.ends[0].d[GY]=on->d[GY];
				  l2.ends[0].d[GZ]=on->d[GZ];

				  on=*((elems+i)->nods+1);
				  l2.ends[1].d[GX]=on->d[GX];
				  l2.ends[1].d[GY]=on->d[GY];
				  l2.ends[1].d[GZ]=on->d[GZ];

				  distancelineline(l1,l2,&n1,&n2,&distance,&flag);

				  if((flag==1 || flag==3) && fabs(distance)<eps)
				  {
					nc1.loff=0;
					nc1.code=0;
					nc1.d[GX]=l1.ends[0].d[GX]+n1*(l1.ends[1].d[GX]-l1.ends[0].d[GX]);
					nc1.d[GY]=l1.ends[0].d[GY]+n1*(l1.ends[1].d[GY]-l1.ends[0].d[GY]);
					nc1.d[GZ]=l1.ends[0].d[GZ]+n1*(l1.ends[1].d[GZ]-l1.ends[0].d[GZ]);

/*
sprintf(str,"CROSS TYPE1/3 X=%.3f Y=%.3f Z=%.3f\nELEM%d:HEAD=%d,TAIL=%d\nELEM%d:HEAD=%d,TAIL=%d",
		nc1.d[GX],nc1.d[GY],nc1.d[GZ],
		(elems+j)->code,
		(*((elems+j)->nods+0))->code,
		(*((elems+j)->nods+1))->code,
		(elems+i)->code,
		(*((elems+i)->nods+0))->code,
		(*((elems+i)->nods+1))->code);
MessageBox(NULL,str,"STACK",MB_OK);
*/

					/*NODE1 FIND SAME NODE*/
					find=0;
					for(k=0;k<nnode;k++)
					{
					  if((nodes+k)->d[GX]-eps<nc1.d[GX] && nc1.d[GX]<(nodes+k)->d[GX]+eps &&
						 (nodes+k)->d[GY]-eps<nc1.d[GY] && nc1.d[GY]<(nodes+k)->d[GY]+eps &&
						 (nodes+k)->d[GZ]-eps<nc1.d[GZ] && nc1.d[GZ]<(nodes+k)->d[GZ]+eps)
					  {
						ncross1=nodes+k;
						find=1;
						break;
					  }
					}
#if 1
					if(!find)
					{
					  nnode++;
					  inode++;

					  /*nodes=(struct onode *)realloc(nodes,nnode*sizeof(struct onode));*/
					  mnodes=nodes;
					  nodes=(struct onode *)malloc(nnode*sizeof(struct onode));
					  for(k=0;k<nnode-1;k++) *(nodes+k)=*(mnodes+k);

					  /*confs=(struct oconf *)realloc(confs,6*nnode*sizeof(struct oconf));*/
					  mconfs=confs;
					  confs=(struct oconf *)malloc(6*nnode*sizeof(struct oconf));
					  for(k=0;k<6*(nnode-1);k++) *(confs+k)=*(mconfs+k);

					  (nodes+nnode-1)->loff=nnode-1;
					  (nodes+nnode-1)->code=inode;
					  (nodes+nnode-1)->d[GX]=nc1.d[GX];
					  (nodes+nnode-1)->d[GY]=nc1.d[GY];
					  (nodes+nnode-1)->d[GZ]=nc1.d[GZ];

					  (confs+6*(nnode-1)+0)->iconf=0;
					  (confs+6*(nnode-1)+1)->iconf=0;
					  (confs+6*(nnode-1)+2)->iconf=0;
					  (confs+6*(nnode-1)+3)->iconf=0;
					  (confs+6*(nnode-1)+4)->iconf=0;
					  (confs+6*(nnode-1)+5)->iconf=0;
					  (confs+6*(nnode-1)+0)->value=0.0;
					  (confs+6*(nnode-1)+1)->value=0.0;
					  (confs+6*(nnode-1)+2)->value=0.0;
					  (confs+6*(nnode-1)+3)->value=0.0;
					  (confs+6*(nnode-1)+4)->value=0.0;
					  (confs+6*(nnode-1)+5)->value=0.0;

					  ncross1=nodes+nnode-1;

					  for(k=0;k<nelem;k++)
					  {
/*
sprintf(str,"k=%d ELEM=%d nban=%d",k,(elems+k)->code,(elems+k)->nban);
MessageBox(NULL,str,"STACK",MB_OK);
*/
						for(l=0;l<(elems+k)->nnod;l++)
						{
						  loff=(*((elems+k)->nods+l))->loff;
						  *((elems+k)->nods+l)=nodes+loff;
						}
						for(l=0;l<(elems+k)->nban;l++)
						{
						  for(m=0;m<((elems+k)->bans+l)->nnod;m++)
						  {
							loff=(*(((elems+k)->bans+l)->nods+m))->loff; /*ERROR*/
							*(((elems+k)->bans+l)->nods+m)=nodes+loff;
						  }
						}
					  }
					  free(mnodes);
					}
#endif
				  }

				  if((flag==2 || flag==3) && fabs(distance)<eps)
				  {
					nc2.loff=0;
					nc2.code=0;
					nc2.d[GX]=l2.ends[0].d[GX]+n2*(l2.ends[1].d[GX]-l2.ends[0].d[GX]);
					nc2.d[GY]=l2.ends[0].d[GY]+n2*(l2.ends[1].d[GY]-l2.ends[0].d[GY]);
					nc2.d[GZ]=l2.ends[0].d[GZ]+n2*(l2.ends[1].d[GZ]-l2.ends[0].d[GZ]);

/*
sprintf(str,"CROSS TYPE2/3 X=%.3f Y=%.3f Z=%.3f\nELEM%d:HEAD=%d,TAIL=%d\nELEM%d:HEAD=%d,TAIL=%d",
		nc2.d[GX],nc2.d[GY],nc2.d[GZ],
		(elems+j)->code,
		(*((elems+j)->nods+0))->code,
		(*((elems+j)->nods+1))->code,
		(elems+i)->code,
		(*((elems+i)->nods+0))->code,
		(*((elems+i)->nods+1))->code);
MessageBox(NULL,str,"STACK",MB_OK);
*/

					/*NODE2 FIND SAME NODE*/
					find=0;
					for(k=0;k<nnode;k++)
					{
					  if((nodes+k)->d[GX]-eps<nc2.d[GX] && nc2.d[GX]<(nodes+k)->d[GX]+eps &&
						 (nodes+k)->d[GY]-eps<nc2.d[GY] && nc2.d[GY]<(nodes+k)->d[GY]+eps &&
						 (nodes+k)->d[GZ]-eps<nc2.d[GZ] && nc2.d[GZ]<(nodes+k)->d[GZ]+eps)
					  {
						ncross2=nodes+k;
						find=1;
						break;
					  }
					}
#if 1
					if(!find)
					{
					  nnode++;
					  inode++;

					  /*nodes=(struct onode *)realloc(nodes,nnode*sizeof(struct onode));*/
					  mnodes=nodes;
					  nodes=(struct onode *)malloc(nnode*sizeof(struct onode));
					  for(k=0;k<nnode-1;k++) *(nodes+k)=*(mnodes+k);

					  /*confs=(struct oconf *)realloc(confs,6*nnode*sizeof(struct oconf));*/
					  mconfs=confs;
					  confs=(struct oconf *)malloc(6*nnode*sizeof(struct oconf));
					  for(k=0;k<6*(nnode-1);k++) *(confs+k)=*(mconfs+k);

					  (nodes+nnode-1)->loff=nnode-1;
					  (nodes+nnode-1)->code=inode;
					  (nodes+nnode-1)->d[GX]=nc2.d[GX];
					  (nodes+nnode-1)->d[GY]=nc2.d[GY];
					  (nodes+nnode-1)->d[GZ]=nc2.d[GZ];

					  (confs+6*(nnode-1)+0)->iconf=0;
					  (confs+6*(nnode-1)+1)->iconf=0;
					  (confs+6*(nnode-1)+2)->iconf=0;
					  (confs+6*(nnode-1)+3)->iconf=0;
					  (confs+6*(nnode-1)+4)->iconf=0;
					  (confs+6*(nnode-1)+5)->iconf=0;
					  (confs+6*(nnode-1)+0)->value=0.0;
					  (confs+6*(nnode-1)+1)->value=0.0;
					  (confs+6*(nnode-1)+2)->value=0.0;
					  (confs+6*(nnode-1)+3)->value=0.0;
					  (confs+6*(nnode-1)+4)->value=0.0;
					  (confs+6*(nnode-1)+5)->value=0.0;

					  ncross2=nodes+nnode-1;

					  for(k=0;k<nelem;k++)
					  {
						for(l=0;l<(elems+k)->nnod;l++)
						{
						  loff=(*((elems+k)->nods+l))->loff;
						  *((elems+k)->nods+l)=nodes+loff;
						}
						for(l=0;l<(elems+k)->nban;l++)
						{
						  for(m=0;m<((elems+k)->bans+l)->nnod;m++)
						  {
							loff=(*(((elems+k)->bans+l)->nods+m))->loff;
							*(((elems+k)->bans+l)->nods+m)=nodes+loff;
						  }
						}
					  }
					  free(mnodes);
					}
#endif
				  }
#if 1
				  /*DIVIDE ELEMENT 1*/
				  if((flag==1 || flag==3) && fabs(distance)<eps)
				  {
					/*DIVIDE ELEMENT*/
					nelem++;
					ielem++;

					/*elems=(struct oelem *)
						  realloc(elems,(nelem+2*count)*sizeof(struct oelem));*/
					melems=elems;
					elems=(struct oelem *)
						  malloc((nelem/*+2*count*/)*sizeof(struct oelem));
					for(k=0;k<nelem-1;k++) *(elems+k)=*(melems+k);

					(elems+nelem-1)->loff=nelem-1;
					(elems+nelem-1)->code=ielem;
					(elems+nelem-1)->sect=(elems+j)->sect;
					(elems+nelem-1)->type=(elems+j)->type;
					(elems+nelem-1)->cangle=(elems+j)->cangle;
					(elems+nelem-1)->role=(elems+j)->role;
					(elems+nelem-1)->nnod=2;
					(elems+nelem-1)->nban=0;
					(elems+nelem-1)->nods
					=(struct onode **)
					 malloc((elems+nelem-1)->nnod*sizeof(struct onode *));

					*((elems+nelem-1)->nods+0)=ncross1;
					*((elems+nelem-1)->nods+1)=*((elems+j)->nods+1);
					*((elems+j)->nods+1)      =ncross1;

					(elems+nelem-1)->bonds
					=(signed char *)
					 malloc(6*((elems+nelem-1)->nnod)*sizeof(signed char));

					*((elems+nelem-1)->bonds+6*0+0)=0;
					*((elems+nelem-1)->bonds+6*0+1)=0;
					*((elems+nelem-1)->bonds+6*0+2)=0;
					*((elems+nelem-1)->bonds+6*0+3)=0;
					*((elems+nelem-1)->bonds+6*0+4)=0;
					*((elems+nelem-1)->bonds+6*0+5)=0;

					*((elems+nelem-1)->bonds+6*1+0)=*((elems+j)->bonds+6*1+0);
					*((elems+nelem-1)->bonds+6*1+1)=*((elems+j)->bonds+6*1+1);
					*((elems+nelem-1)->bonds+6*1+2)=*((elems+j)->bonds+6*1+2);
					*((elems+nelem-1)->bonds+6*1+3)=*((elems+j)->bonds+6*1+3);
					*((elems+nelem-1)->bonds+6*1+4)=*((elems+j)->bonds+6*1+4);
					*((elems+nelem-1)->bonds+6*1+5)=*((elems+j)->bonds+6*1+5);

					*((elems+j)->bonds+6*1+0)=0;
					*((elems+j)->bonds+6*1+1)=0;
					*((elems+j)->bonds+6*1+2)=0;
					*((elems+j)->bonds+6*1+3)=0;
					*((elems+j)->bonds+6*1+4)=0;
					*((elems+j)->bonds+6*1+5)=0;
/*
sprintf(str,"DIVIDED\nELEM%d:HEAD=%d,TAIL=%d\nELEM%d:HEAD=%d,TAIL=%d",
		(elems+j)->code,
		(*((elems+j)->nods+0))->code,
		(*((elems+j)->nods+1))->code,
		(elems+nelem-1)->code,
		(*((elems+nelem-1)->nods+0))->code,
		(*((elems+nelem-1)->nods+1))->code);
MessageBox(NULL,str,"STACK",MB_OK);
*/
				  }

				  /*DIVIDE ELEMENT 2*/
				  if((flag==2 || flag==3) && fabs(distance)<eps)
				  {
					/*DIVIDE ELEMENT*/
					nelem++;
					ielem++;

					/*elems=(struct oelem *)
						  realloc(elems,(nelem+2*count)*sizeof(struct oelem));*/
					melems=elems;
					elems=(struct oelem *)
						  malloc((nelem/*+2*count*/)*sizeof(struct oelem));
					for(k=0;k<nelem-1;k++) *(elems+k)=*(melems+k);

					(elems+nelem-1)->loff=nelem-1;
					(elems+nelem-1)->code=ielem;
					(elems+nelem-1)->sect=(elems+i)->sect;
					(elems+nelem-1)->type=(elems+i)->type;
					(elems+nelem-1)->cangle=(elems+i)->cangle;
					(elems+nelem-1)->role=(elems+i)->role;
					(elems+nelem-1)->nnod=2;
					(elems+nelem-1)->nban=0;
					(elems+nelem-1)->nods
					=(struct onode **)
					 malloc((elems+nelem-1)->nnod*sizeof(struct onode *));

					*((elems+nelem-1)->nods+0)=ncross2;
					*((elems+nelem-1)->nods+1)=*((elems+i)->nods+1);
					*((elems+i)->nods+1)      =ncross2;

					(elems+nelem-1)->bonds
					=(signed char *)
					 malloc(6*((elems+nelem-1)->nnod)*sizeof(signed char));

					*((elems+nelem-1)->bonds+6*0+0)=0;
					*((elems+nelem-1)->bonds+6*0+1)=0;
					*((elems+nelem-1)->bonds+6*0+2)=0;
					*((elems+nelem-1)->bonds+6*0+3)=0;
					*((elems+nelem-1)->bonds+6*0+4)=0;
					*((elems+nelem-1)->bonds+6*0+5)=0;

					*((elems+nelem-1)->bonds+6*1+0)=*((elems+j)->bonds+6*1+0);
					*((elems+nelem-1)->bonds+6*1+1)=*((elems+j)->bonds+6*1+1);
					*((elems+nelem-1)->bonds+6*1+2)=*((elems+j)->bonds+6*1+2);
					*((elems+nelem-1)->bonds+6*1+3)=*((elems+j)->bonds+6*1+3);
					*((elems+nelem-1)->bonds+6*1+4)=*((elems+j)->bonds+6*1+4);
					*((elems+nelem-1)->bonds+6*1+5)=*((elems+j)->bonds+6*1+5);

					*((elems+j)->bonds+6*1+0)=0;
					*((elems+j)->bonds+6*1+1)=0;
					*((elems+j)->bonds+6*1+2)=0;
					*((elems+j)->bonds+6*1+3)=0;
					*((elems+j)->bonds+6*1+4)=0;
					*((elems+j)->bonds+6*1+5)=0;
/*
sprintf(str,"DIVIDED\nELEM%d:HEAD=%d,TAIL=%d\nELEM%d:HEAD=%d,TAIL=%d",
		(elems+i)->code,
		(*((elems+i)->nods+0))->code,
		(*((elems+i)->nods+1))->code,
		(elems+nelem-1)->code,
		(*((elems+nelem-1)->nods+0))->code,
		(*((elems+nelem-1)->nods+1))->code);
MessageBox(NULL,str,"STACK",MB_OK);
*/
				  }
#endif
				}
			  }
			}
		  }
/*MessageBox(NULL,"Pass1.3","Stack",MB_OK);*/

		  org->nnode=nnode;
		  org->nodes=nodes;
		  org->confs=confs;
		  org->nelem=nelem;
		  org->elems=elems;

/*
sprintf(str,"NNODE=%d NELEM=%d",org->nnode,org->nelem);
MessageBox(NULL,str,"STACK",MB_OK);
*/
/*
for(k=0;k<org->nnode;k++)
{
sprintf(str,"NNODE=%d\nNODE=%d\nICONF=%d,%d,%d,%d,%d,%d\nVCONF=%.3f,%.3f,%.3f,%.3f,%.3f,%.3f",
		org->nnode,
		(org->nodes+k)->code,
        (org->confs+6*k+0)->iconf,
		(org->confs+6*k+1)->iconf,
        (org->confs+6*k+2)->iconf,
        (org->confs+6*k+3)->iconf,
        (org->confs+6*k+4)->iconf,
        (org->confs+6*k+5)->iconf,
        (org->confs+6*k+0)->value,
        (org->confs+6*k+1)->value,
        (org->confs+6*k+2)->value,
        (org->confs+6*k+3)->value,
        (org->confs+6*k+4)->value,
        (org->confs+6*k+5)->value);
MessageBox(NULL,str,"STACK",MB_OK);
}
*/
/*
for(k=0;k<org->nelem;k++)
{
sprintf(str,"NELEM=%d\nELEM%d:HEAD=%d,TAIL=%d",
        org->nelem,
        (org->elems+k)->code,
        (*((org->elems+k)->nods+0))->code,
        (*((org->elems+k)->nods+1))->code);
MessageBox(NULL,str,"STACK",MB_OK);
}
*/
}
#endif



  /*MASS, AI EARTHQUAKE, CMQ*/
  /*荷重集計*/
  fout=fopen("hogtxt.wgt","w");
  weightdistribution(NULL,fout,NULL,org);
  fclose(fout);

  /* 150428 fukushima for original section */
  for(i=0;i<(org->nsect);i++)
  {
	(org->sects+i)->ocode=(org->sects+i)->code;
  }


  /*WALL,SLAB INTO BRACE*/
  nelem=org->nelem;
  count=0;
#if 0
  for(i=0;i<nelem;i++)
  {
    if((org->elems+i)->type==WALL || (org->elems+i)->type==SLAB)
    {
	  if((org->elems+i)->nnod==4 && (((org->elems+i)->sect->figs+0)->prop->E)>0.0) count++;  /*UJIOKA FOR OPTIMIZATION*/
    }
  }
#endif

  elems=(struct oelem *)malloc((nelem+2*count)*sizeof(struct oelem));
  if(elems==NULL)
  {
    MessageBox(NULL,"Buffer Null.","Extract",MB_OK);
    return 0;
  }

  for(i=0;i<nelem;i++)
  {
    *(elems+i)=*(org->elems+i);
  }
  /*free(org->elems);*/

  /*org->nelem=nelem+2*count;*/
  org->elems=elems;

/*globalfile=fopen("hogtxt.brc","w");*/

#if 0
  j=0;
  for(i=0;i<nelem;i++)
  {
    if((org->elems+i)->type==WALL || (org->elems+i)->type==SLAB)
    {
	  if((org->elems+i)->nnod==4 && (((org->elems+i)->sect->figs+0)->prop->E)>0.0)      /*UJIOKA FOR OPTIMIZATION*/
	  {
		recttobrace((org->elems+i),org,(nelem+2*j),2,1.0);
		j++;
	  }
	}
  }
#endif

  org->nelem=nelem+2*count;

/*fclose(globalfile);*/







  /*SECTIONS*/
  nsect=0;
  for(i=0;i<(org->nsect);i++) /*COUNT*/
  {
	s=(org->sects+i);
    sectionfeatures(s);
    if((s->area)>0.0) nsect++;
  }
  az->nsect=nsect;
  ax->nsect=nsect;
  ay->nsect=nsect;
  az->sects=(struct osect *)malloc(nsect*sizeof(struct osect));
  ax->sects=(struct osect *)malloc(nsect*sizeof(struct osect));
  ay->sects=(struct osect *)malloc(nsect*sizeof(struct osect));
  j=0;
  for(i=0;i<(org->nsect);i++)
  {
    s=(org->sects+i);

	/*EXTRACT IF A>0*/
    if((s->area)>1.0E-08)
    {
	  *(az->sects+j)=*s;
      *(ax->sects+j)=*s;
      *(ay->sects+j)=*s;

      (az->sects+j)->loff=j;
	  (ax->sects+j)->loff=j;
      (ay->sects+j)->loff=j;
      (az->sects+j)->dflag=1;
      (ax->sects+j)->dflag=1;
      (ay->sects+j)->dflag=1;

#if 1
if(!strncmp(prj,"vene",4))
{
  if((200<=(az->sects+j)->code)&&((az->sects+j)->code<=299))
  {
    (az->sects+j)->Ixx=0.0;
    (az->sects+j)->Iyy=0.0;
    (az->sects+j)->Jzz=0.0;
  }
  if((600<=(az->sects+j)->code)&&((az->sects+j)->code<=699)) (az->sects+j)->E=0.0000001;
}
#endif
#if 0
if(!strncmp(prj,"vene",4))
{
  if((200<=(az->sects+j)->code)&&((az->sects+j)->code<=299))
  {
	(az->sects+j)->Ixx=0.00000001;
	(az->sects+j)->Iyy=0.00000001;
	(az->sects+j)->Jzz=0.00000001;
  }
  if((600<=(az->sects+j)->code)&&((az->sects+j)->code<=699)) (az->sects+j)->area=0.0;
}
#endif

	  j++;
	}
  }

/*
sprintf(str,"NNODE=%d NELEM=%d",org->nnode,org->nelem);
MessageBox(NULL,str,"STACK",MB_OK);
*/

  /*ELEMENTS*/
  nelem=0;
  for(i=0;i<(org->nelem);i++)
  {
	/*currentpivot((org->elems+i)->code,
				 (org->elems+i)->sect->code);*/

	if((org->elems+i)->role==ROLESTRESS ||
	   (org->elems+i)->role==ROLENULL)
	{
	  for(j=0;j<nsect;j++)
	  {
		if(((org->elems+i)->sect->code)==((az->sects+j)->code))
		{
		  nelem++;
		  break;
		}
	  }
	}
  }
/*
sprintf(str,"NELEM=%d",nelem);
MessageBox(NULL,str,"STACK",MB_OK);
*/

  az->nelem=nelem;
  ax->nelem=nelem;
  ay->nelem=nelem;
  az->elems=(struct owire *)malloc(nelem*sizeof(struct owire));
  ax->elems=(struct owire *)malloc(nelem*sizeof(struct owire));
  ay->elems=(struct owire *)malloc(nelem*sizeof(struct owire));
  k=0;
  for(i=0;i<(org->nelem);i++)
  {
	/*EXTRACT IF ROLE NULL,STRESS AND EXTRACTED SECTION.*/
	if((org->elems+i)->role==ROLESTRESS ||
	   (org->elems+i)->role==ROLENULL)
	{
	  if((org->elems+i)->nnod==2)
	  {
		for(j=0;j<nsect;j++)
		{
		  if(((org->elems+i)->sect->code)==((az->sects+j)->code))
		  {
			(az->elems+k)->code   =(org->elems+i)->code;
			(ax->elems+k)->code   =(org->elems+i)->code;
			(ay->elems+k)->code   =(org->elems+i)->code;

			(az->elems+k)->loff   =(org->elems+i)->loff;
			(ax->elems+k)->loff   =(org->elems+i)->loff;
			(ay->elems+k)->loff   =(org->elems+i)->loff;

			(az->elems+k)->cangle =(org->elems+i)->cangle;
			(ax->elems+k)->cangle =(org->elems+i)->cangle;
			(ay->elems+k)->cangle =(org->elems+i)->cangle;

			(az->elems+k)->node[0]=*((org->elems+i)->nods+0);
			(ax->elems+k)->node[0]=*((org->elems+i)->nods+0);
            (ay->elems+k)->node[0]=*((org->elems+i)->nods+0);

            (az->elems+k)->node[1]=*((org->elems+i)->nods+1);
            (ax->elems+k)->node[1]=*((org->elems+i)->nods+1);
            (ay->elems+k)->node[1]=*((org->elems+i)->nods+1);

			for(m=0;m<6;m++)
			{
              (az->elems+k)->iconf[0][m]=*((org->elems+i)->bonds+m);
              (ax->elems+k)->iconf[0][m]=*((org->elems+i)->bonds+m);
              (ay->elems+k)->iconf[0][m]=*((org->elems+i)->bonds+m);

              (az->elems+k)->iconf[1][m]=*((org->elems+i)->bonds+6+m);
              (ax->elems+k)->iconf[1][m]=*((org->elems+i)->bonds+6+m);
              (ay->elems+k)->iconf[1][m]=*((org->elems+i)->bonds+6+m);

              (az->elems+k)->stress[0][m]=((org->elems+i)->initial[0][m]);
              (az->elems+k)->stress[1][m]=((org->elems+i)->initial[1][m]);

              (ax->elems+k)->stress[0][m]=0.0;
              (ax->elems+k)->stress[1][m]=0.0;

              (ay->elems+k)->stress[0][m]=0.0;
              (ay->elems+k)->stress[1][m]=0.0;
            }

            (az->elems+k)->loff=k;
            (ax->elems+k)->loff=k;
            (ay->elems+k)->loff=k;

            /*SET SECTION.*/
            (az->elems+k)->sect=az->sects+j;
            (ax->elems+k)->sect=ax->sects+j;
            (ay->elems+k)->sect=ay->sects+j;

            (az->elems+k)->sect->type=(org->elems+i)->type;
            (ax->elems+k)->sect->type=(org->elems+i)->type;
            (ay->elems+k)->sect->type=(org->elems+i)->type;

            /*INITIAL ENERGY.*/
            (az->elems+k)->Ee[0]=0.0;
            (ax->elems+k)->Ee[0]=0.0;
            (ay->elems+k)->Ee[0]=0.0;
            (az->elems+k)->Ee[1]=0.0;
            (ax->elems+k)->Ee[1]=0.0;
            (ay->elems+k)->Ee[1]=0.0;

            (az->elems+k)->Ep[0]=0.0;
			(ax->elems+k)->Ep[0]=0.0;
            (ay->elems+k)->Ep[0]=0.0;
            (az->elems+k)->Ep[1]=0.0;
            (ax->elems+k)->Ep[1]=0.0;
            (ay->elems+k)->Ep[1]=0.0;

            /*SAFETY RATE*/
            (az->elems+k)->srate[0]=0.0;
            (az->elems+k)->srate[1]=0.0;
            (az->elems+k)->srate[2]=0.0;
            (az->elems+k)->srate[3]=0.0;
            (ax->elems+k)->srate[0]=0.0;
            (ax->elems+k)->srate[1]=0.0;
            (ax->elems+k)->srate[2]=0.0;
            (ax->elems+k)->srate[3]=0.0;
            (ay->elems+k)->srate[0]=0.0;
            (ay->elems+k)->srate[1]=0.0;
            (ay->elems+k)->srate[2]=0.0;
            (ay->elems+k)->srate[3]=0.0;

            k++;
            break;
          }
        }
      }
    }
  }

  az->nshell=0;
  ax->nshell=0;
  ay->nshell=0;


  /*NODE COORDINATES,LOADS*/
  nnode=0;
  for(i=0;i<(org->nnode);i++)
  {
	/*EXTRACT IF BELONGS TO EXTRACTED ELEMENTS.*/
    for(j=0;j<nelem;j++)
    {
      if((org->nodes+i)->code==(az->elems+j)->node[0]->code ||
         (org->nodes+i)->code==(az->elems+j)->node[1]->code)
      {
        nnode++;
        break;
      }
    }
  }
  az->nnode=nnode;
  ax->nnode=nnode;
  ay->nnode=nnode;

  az->nodes=(struct onode *)malloc(nnode*sizeof(struct onode));
  ax->nodes=(struct onode *)malloc(nnode*sizeof(struct onode));
  ay->nodes=(struct onode *)malloc(nnode*sizeof(struct onode));

  az->ninit=(struct onode *)malloc(nnode*sizeof(struct onode));
  ax->ninit=(struct onode *)malloc(nnode*sizeof(struct onode));
  ay->ninit=(struct onode *)malloc(nnode*sizeof(struct onode));

  az->confs=(struct oconf *)malloc(6*nnode*sizeof(struct oconf));
  ax->confs=(struct oconf *)malloc(6*nnode*sizeof(struct oconf));
  ay->confs=(struct oconf *)malloc(6*nnode*sizeof(struct oconf));

  az->nmass=(double *)malloc(nnode*sizeof(double));
  ax->nmass=(double *)malloc(nnode*sizeof(double));
  ay->nmass=(double *)malloc(nnode*sizeof(double));

  org->ai.lbound=(double *)malloc((org->ai.nfloor+1)*sizeof(double));

  nreact=0;
  k=0;
  for(i=0;i<(org->nnode);i++)
  {
    count=0;

    /*EXTRACT IF BELONGS TO EXTRACTED ELEMENTS.*/
    for(j=0;j<nelem;j++)
    {
      for(m=0;m<2;m++)
      {
		if((org->nodes+i)->code==(az->elems+j)->node[m]->code)
        {
          *(az->nodes+k)=*(org->nodes+i);
          *(ax->nodes+k)=*(org->nodes+i);
          *(ay->nodes+k)=*(org->nodes+i);

		  (az->nodes+k)->loff=k;
		  (ax->nodes+k)->loff=k;
		  (ay->nodes+k)->loff=k;

		  *(az->ninit+k)=*(az->nodes+k);
		  *(ax->ninit+k)=*(ax->nodes+k);
		  *(ay->ninit+k)=*(ay->nodes+k);


		  (az->elems+j)->node[m]=(az->nodes+k); /*SET ELEMENT NODE*/
		  (ax->elems+j)->node[m]=(ax->nodes+k); /*SET ELEMENT NODE*/
		  (ay->elems+j)->node[m]=(ay->nodes+k); /*SET ELEMENT NODE*/

          /*DIRECTLY INPUT LOAD.*/
          for(n=0;n<6;n++)
          {
            oloff=6*i+n;
            aloff=6*k+n;

            (az->confs+aloff)->iconf=(org->confs+oloff)->iconf;
            (ax->confs+aloff)->iconf=(org->confs+oloff)->iconf;
            (ay->confs+aloff)->iconf=(org->confs+oloff)->iconf;

#if 1
if(!strncmp(prj,"vene",4))
{
  if((az->nodes+k)->code==111 && n==0) (az->confs+aloff)->iconf=1;
  if((az->nodes+k)->code==412 && n==0) (az->confs+aloff)->iconf=1;
  if((az->nodes+k)->code==412 && n==1) (az->confs+aloff)->iconf=1;
  if((az->nodes+k)->code==413 && n==1) (az->confs+aloff)->iconf=1;
}
#endif

//            if(n==2) /*ALREADY ADDED IN WEIGHTDISTRIBUTION.*/
            if(n==2 && !((az->confs+aloff)->iconf)) /*MODIFIED BY UJIOKA.*/
            {
              (az->confs+aloff)->value=0.0;
              (ax->confs+aloff)->value=0.0;
			  (ay->confs+aloff)->value=0.0;
            }
            else
            {
              (az->confs+aloff)->value=(org->confs+oloff)->value;
              (ax->confs+aloff)->value=(org->confs+oloff)->value;
              (ay->confs+aloff)->value=(org->confs+oloff)->value;
            }
			/*if(!strncmp(prj,"lunarbase",9))
            {
			  if(((az->confs+aloff)->iconf))    (az->confs+aloff)->value=0.0;
			  if(((ax->confs+aloff)->iconf))    (ax->confs+aloff)->value=0.0;
			  if(((ay->confs+aloff)->iconf))    (ay->confs+aloff)->value=0.0;
			}*/
#if 1
if(!strncmp(prj,"vene",4))
{
  /*(az->confs+aloff)->value=(org->confs+oloff)->value;*/
  (ax->confs+aloff)->value=(org->confs+oloff)->value;
  (ay->confs+aloff)->value=(org->confs+oloff)->value;
}
#endif

            if((az->confs+aloff)->iconf==1) nreact++;
          }

          /*HORIZONTAL X LOAD*/
          aloff=6*k+0;
          if((ax->confs+aloff)->iconf==0)
          {
#if 0
if(!strncmp(prj,"vene",4)) (ax->confs+aloff)->value=0.0;
#endif
            (ax->confs+aloff)->value+=(org->loads+i)->w[WEIGHTEQX];
			if(!strncmp(prj,"lunarbase",9))
            {
			(ax->confs+aloff)->value=0.0;
			(ay->confs+aloff)->value=0.0;
			}
			*(ax->nmass+k)=(org->loads+i)->mass;
          }
          else *(ax->nmass+k)=0.0;

          /*HORIZONTAL Y LOAD*/
          aloff=6*k+1;
          if((ay->confs+aloff)->iconf==0)
          {
            (ay->confs+aloff)->value+=(org->loads+i)->w[WEIGHTEQY];
			if(!strncmp(prj,"lunarbase",9))
            {
			(ax->confs+aloff)->value=0.0;
			(ay->confs+aloff)->value=0.0;
			}
            *(ay->nmass+k)=(org->loads+i)->mass;
          }
          else *(ay->nmass+k)=0.0;

		  /*VERTICAL Z LOAD:FOR LONG PERIOD*/
		  aloff=6*k+2;
          if((az->confs+aloff)->iconf==0)
          {
            (az->confs+aloff)->value-=(org->loads+i)->w[WEIGHTFRAME];
            *(az->nmass+k)=(org->loads+i)->mass;
          }
          else *(az->nmass+k)=0.0;

          k++;
          count++;
          break;
		}
      }
	  if(count>0) break;
    }
  }


  az->nconstraint=0;
  ax->nconstraint=0;
  ay->nconstraint=0;

  az->constraintmain=(long int *)malloc(6*nnode*sizeof(long int));
  ax->constraintmain=(long int *)malloc(6*nnode*sizeof(long int));
  ay->constraintmain=(long int *)malloc(6*nnode*sizeof(long int));

  for (i = 0; i < 6*nnode; i++)
  {
	*((az->constraintmain) + i) = i;
	*((ax->constraintmain) + i) = i;
	*((ay->constraintmain) + i) = i;
  }

  az->nreact=nreact;
  ax->nreact=nreact;
  ay->nreact=nreact;

 if(globalmessageflag==1)
  {
	sprintf(str,"Nodes=%d Elems=%d Sects=%d",nnode,nelem,nsect);
	MessageBox(NULL,str,"Extract",MB_OK);
  }

#if 1
if(!strncmp(prj,"inax",4) || !strncmp(prj,"stack",5))
{
  for(k=0;k<pnelem;k++)
  {
    *(pelems+k)=*(penods+k);
    /*(pelems+k)->nnod=(penods+k)->nnod;*/
    /*(pelems+k)->nban=(penods+k)->nban;*/
    for(l=0;l<(penods+k)->nnod;l++)
    {
      *((pelems+k)->nods+l)=*((penods+k)->nods+l);

      *((pelems+k)->bonds+6*l+0)=*((penods+k)->bonds+6*l+0);
      *((pelems+k)->bonds+6*l+1)=*((penods+k)->bonds+6*l+1);
	  *((pelems+k)->bonds+6*l+2)=*((penods+k)->bonds+6*l+2);
      *((pelems+k)->bonds+6*l+3)=*((penods+k)->bonds+6*l+3);
      *((pelems+k)->bonds+6*l+4)=*((penods+k)->bonds+6*l+4);
      *((pelems+k)->bonds+6*l+5)=*((penods+k)->bonds+6*l+5);
    }
    for(l=0;l<(penods+k)->nban;l++)
    {
      ((pelems+k)->bans+l)->nnod=((penods+k)->bans+l)->nnod;
      for(m=0;m<((penods+k)->bans+l)->nnod;m++)
      {
        *(((pelems+k)->bans+l)->nods+m)
        =*(((penods+k)->bans+l)->nods+m);
      }
    }
  }
  org->nnode=pnnode;
  org->nodes=pnodes;
  org->confs=pconfs;
  /*org->nsect=pnsect;*/
  /*org->sects=psects;*/
  org->nelem=pnelem;
  org->elems=pelems;
}
#endif

  return 1;
}/*extractarclmfromorgan*/

struct oelem *recttobrace(struct oelem *oe,
                          struct organ *org,
                          int loff,int nbrace, /*NUMBER OF BRACE*/
                          double rfact /*RIGID FACTOR*/)
/*RECTANGLE BAN TO 2 BRACE*/
{
  char str[256];
  int i,j,nnod,find;
  double dx,dy,dz;
  double thick,poi,length,height,wrate,wrate1,wrate2,Ae;
  struct osect *ps;
  struct osect os={0,0,0,NULL,1,
                   1,/*NFIG*/
                   0.0, 0.0, 0.0, 0.0, 0.0, 0.0,/*COEFFICIENTS*/
				   0.0,/*EXP*/
				   { 1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0},
				   {-1000.0,-1000.0,-1000.0,-1000.0,-1000.0,-1000.0},
				   0.0,1.0,1.0,
				   {0.0, 0.0, 0.0}/*HIJU*/};
  struct oelem *obraces;
  struct oelem oinit={0,0,ROLENULL,TYPENULL,
					  0,0,/*NNOD,NBAN*/
                      0.0,
                      {{0.0, 0.0, 0.0, 0.0, 0.0, 0.0},/*CMQ HEAD*/
                       {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}}/*CMQ TAIL*/};

  if((oe->nnod)<4 || (oe->nban)<1) return NULL;
  if((oe->nnod)>4 || (oe->nban)>1) return NULL;
  if(nbrace!=2) return NULL;

  obraces=(struct oelem *)malloc(nbrace*sizeof(struct oelem));
  for(i=0;i<nbrace;i++)
  {
    *(obraces+i)=oinit;
    nnod=2;
    (obraces+i)->nnod=nnod;
    (obraces+i)->nods=(struct onode **)
                      malloc(nnod*sizeof(struct onode *));
    (obraces+i)->bonds=(signed char *)
                       malloc(6*nnod*sizeof(signed char));

    for(j=0;j<6*nnod;j++) *((obraces+i)->bonds+j)=0;

    (obraces+i)->type=oe->type;
  }

  *((obraces+0)->nods+0)=*(oe->bans->nods+0);
  *((obraces+0)->nods+1)=*(oe->bans->nods+2);
  *((obraces+1)->nods+0)=*(oe->bans->nods+1);
  *((obraces+1)->nods+1)=*(oe->bans->nods+3);

  thick=(oe->sect->figs+0)->thick;
  poi  =(oe->sect->figs+0)->prop->poi;

  dx=((*(oe->bans->nods+1))->d[0])-((*(oe->bans->nods+0))->d[0]);
  dy=((*(oe->bans->nods+1))->d[1])-((*(oe->bans->nods+0))->d[1]);
  dz=((*(oe->bans->nods+1))->d[2])-((*(oe->bans->nods+0))->d[2]);
  length=sqrt(dx*dx+dy*dy+dz*dz)-(oe->lface[0])-(oe->lface[1]);

  dx=((*(oe->bans->nods+3))->d[0])-((*(oe->bans->nods+0))->d[0]);
  dy=((*(oe->bans->nods+3))->d[1])-((*(oe->bans->nods+0))->d[1]);
  dz=((*(oe->bans->nods+3))->d[2])-((*(oe->bans->nods+0))->d[2]);
  height=sqrt(dx*dx+dy*dy+dz*dz)-(oe->hface[0])-(oe->hface[1]);

  /*RATE OF WINDOW*/
  if(length<=0.0)
  {
    length=0.0;
    Ae=0.0;
  }
  else if(height<=0.0)
  {
    height=0.0;
    Ae=0.0;
  }
  else
  {
    wrate1=(oe->wrect[0])/length;
    wrate2=1.25*sqrt((oe->wrect[0])*(oe->wrect[1])
                     /(length*height));
    if(wrate1>=wrate2) wrate=wrate1;
    else               wrate=wrate2;

    Ae=(1.0-wrate)*rfact/1.2
       *pow((length*length+height*height),1.5)*thick
       /(4.0*(1.0+poi)*length*height); /*1 OF 2 BRACES*/
    Ae*=2.0/(double)nbrace;
  }

  /*sprintf(str,"t =%.5f Lw=%.5f H =%.5f Wl=%.5f Wh=%.5f Ae=%.5f",
          thick,length,height,oe->wrect[0],oe->wrect[1],Ae);*/
  /*if(globalfile!=NULL) fprintf(globalfile,"%s\n",str);*/
  /*MessageBox(NULL,str,"Brace",MB_OK);*/

  os.role=ROLENULL;
  os.type=oe->type;
  os.figs=(struct ofigs *)malloc(sizeof(struct ofigs));
  (os.figs+0)->code =1;
  (os.figs+0)->loff =0;
  (os.figs+0)->area =Ae;
  (os.figs+0)->Ixx  =0.0;
  (os.figs+0)->Iyy  =0.0;
  (os.figs+0)->Jzz  =0.0;
  (os.figs+0)->thick=0.0;
  (os.figs+0)->prop =(oe->sect->figs+0)->prop;
  os.ocode=oe->sect->code; // 200225 furuichi for original section

  find=0;
#if 1
  for(i=0;i<(org->nsect);i++)
  {
    if(comparefirstfig(&os,(org->sects+i)))
    {
      ps=org->sects+i;
      find=1;
      break;
    }
  }
#endif
#if 0
  if(oe->sect->bsect>0)
  {
    for(i=0;i<(org->nsect);i++)
    {
      if((org->sects+i)->code==oe->sect->bsect)
      {
        ps=org->sects+i;
        find=1;
        break;
      }
    }
  }
  else
  {
    for(i=0;i<(org->nsect);i++)
    {
      if((org->sects+i)->ocode==oe->sect->code && comparefirstfig(&os,(org->sects+i)))
      {
        ps=org->sects+i;
        find=1;
        break;
      }
    }
  }
#endif
  if(!find)
  {
    os.code=(org->sects+(org->nsect-1))->code+1;
    os.ocode=oe->sect->code;

    if(oe->type==WALL)
    {
      sprintf(str,"WallBrace");
    }
    else if(oe->type==SLAB)
    {
      sprintf(str,"SlabBrace");
    }
    else
    {
      sprintf(str,"Brace");
    }
    os.dcolor=oe->sect->dcolor;
    os.name=(char *)malloc((strlen(str)+1)*sizeof(char));
    strcpy(os.name,str);

    ps=addsection(&os,org);
  }

  (obraces+0)->sect=ps;
  (obraces+1)->sect=ps;

  for(i=0;i<nbrace;i++)
  {
    (obraces+i)->sect=ps;
    (obraces+i)->code=(org->elems+(loff-1+i))->code+1;
    (obraces+i)->loff=loff+i;

    /*addelement((obraces+i),org);*/
    *(org->elems+loff+i)=*(obraces+i);
  }

  return obraces;
}/*recttobrace*/

int comparesection(struct osect *s1,struct osect *s2)
/*COMPARE SECTION*/
{
  int i;

  if(s1->nfig != s2->nfig) return 0;

  for(i=0;i<(s1->nfig);i++)
  {
    if((s1->figs+i)->area  != (s2->figs+i)->area)  return 0;
    if((s1->figs+i)->Ixx   != (s2->figs+i)->Ixx)   return 0;
    if((s1->figs+i)->Iyy   != (s2->figs+i)->Iyy)   return 0;
    if((s1->figs+i)->Jzz   != (s2->figs+i)->Jzz)   return 0;
    if((s1->figs+i)->thick != (s2->figs+i)->thick) return 0;

    if((s1->figs+i)->prop->hiju != (s2->figs+i)->prop->hiju) return 0;
    if((s1->figs+i)->prop->E    != (s2->figs+i)->prop->E)    return 0;
    if((s1->figs+i)->prop->poi  != (s2->figs+i)->prop->poi)  return 0;
  }

  /*COMPARE POLYPOLYCURVE*/

  return 1;
}/*comparesection*/

int comparefirstfig(struct osect *s1,struct osect *s2)
/*COMPARE FIRST FIGURE*/
{
  double eps=1.0E-12;

  if(s1->type!=s2->type) return 0;
  if(s1->ocode!=s2->ocode) return 0; // 200225 furuichi for original section

  if((s1->nfig)<1 || (s2->nfig)<1) return 0;

  if((s1->figs+0)->area  < (s2->figs+0)->area -eps ||
     (s1->figs+0)->area  > (s2->figs+0)->area +eps)  return 0;
  if((s1->figs+0)->Ixx   < (s2->figs+0)->Ixx  -eps ||
     (s1->figs+0)->Ixx   > (s2->figs+0)->Ixx  +eps)  return 0;
  if((s1->figs+0)->Iyy   < (s2->figs+0)->Iyy  -eps ||
     (s1->figs+0)->Iyy   > (s2->figs+0)->Iyy  +eps)  return 0;
  if((s1->figs+0)->Jzz   < (s2->figs+0)->Jzz  -eps ||
     (s1->figs+0)->Jzz   > (s2->figs+0)->Jzz  +eps)  return 0;
  if((s1->figs+0)->thick < (s2->figs+0)->thick-eps ||
     (s1->figs+0)->thick > (s2->figs+0)->thick+eps)  return 0;

  if((s1->figs+0)->prop->hiju < (s2->figs+0)->prop->hiju-eps ||
	 (s1->figs+0)->prop->hiju > (s2->figs+0)->prop->hiju+eps) return 0;
  if((s1->figs+0)->prop->E    < (s2->figs+0)->prop->E   -eps ||
     (s1->figs+0)->prop->E    > (s2->figs+0)->prop->E   +eps) return 0;
  if((s1->figs+0)->prop->poi  < (s2->figs+0)->prop->poi -eps ||
     (s1->figs+0)->prop->poi  > (s2->figs+0)->prop->poi +eps) return 0;

  /*COMPARE POLYPOLYCURVE*/

  return 1;
}/*comparefirstfig*/

void slabdivision1(HDC hdc,struct viewparam vp,struct obans gban)
/*DIVIDE SLAB BAN INTO CMQ PARTS.....WRONG VERSION.*/
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
     ddot==NULL)
  {
    MessageBox(NULL,"Malloc Failed.","Division",MB_OK);
    return;
  }

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
                         &param1,&param2,&distance,NULL);

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
                         &param1,&param2,&distance,NULL);

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
}/*slabdivision1*/

void slabdivision2(HDC hdc,struct viewparam *vp,
                   struct obans gban,struct obans *cmqbans)
/*DIVIDE SLAB BAN INTO CMQ PARTS.*/
/*cmqbans:DIVIDED BANS OF EACH HEN.*/
{
  /*char non[256];*/
  int ndot,bdot,idot,jdot,kdot;
  int ihen,jhen,khen;
  int lall,dall;
  int jj;
  int ifound;
  double param1,param2,distance,sub0,sub1;
  double x,y,x0,y0,x1,y1;
  double **drccos,**tdrccos;
  double eps=1.0E-08;
  struct onode *lnods;
  struct obans ban;
  struct obans *dban; /*DIVIDED BANS.*/
  struct line dline,rline,*hen0,*hen1;
  struct cmqdot *dot1;
  /*struct cmqelem cmqsum;*/
  struct vector vd,vz,vo;

  if(gban.nods==NULL) return;

  ndot=gban.nnod;

  /*PROJECTED BAN*/
  ban.nods=(struct onode **)malloc(ndot*sizeof(struct onode *));
  lnods=(struct onode *)malloc(ndot*sizeof(struct onode));
  hen0=(struct line *)malloc(ndot*sizeof(struct line));

  /*DIVISION LINES FOR EACH EDE*/
  lall=ndot-1;
  hen1=(struct line *)malloc(lall*sizeof(struct line));

  /*DIVISION DOTS FOR EACH EDGE*/
  dall=(int)(((ndot-1)*(ndot-2))/2);
  dot1=(struct cmqdot *)malloc(dall*sizeof(struct cmqdot));

  /*DIVIDED BANS*/
  dban=(struct obans *)malloc(ndot*sizeof(struct obans));

  if(ban.nods==NULL ||
     lnods==NULL ||
     hen0==NULL ||
     hen1==NULL ||
     dot1==NULL ||
     dban==NULL)
  {
    MessageBox(NULL,"Malloc Failed.","Division",MB_OK);
    return;
  }

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

  for(ihen=0;ihen<ndot;ihen++) /*FOR EACH EDGE.*/
  {
    jhen=ihen;
    for(jj=0;jj<(ndot-1);jj++) /*DIVISION LINE WITH OTHER EDGES.*/
    {
      jhen++;
      if(jhen>=ndot) jhen=0;

      rline.ends[0].d[EX]=(hen0+jhen)->ends[1].d[EX]; /*REVERSE.*/
      rline.ends[0].d[EY]=(hen0+jhen)->ends[1].d[EY];
      rline.ends[0].d[EZ]=(hen0+jhen)->ends[1].d[EZ];
      rline.ends[1].d[EX]=(hen0+jhen)->ends[0].d[EX];
      rline.ends[1].d[EY]=(hen0+jhen)->ends[0].d[EY];
      rline.ends[1].d[EZ]=(hen0+jhen)->ends[0].d[EZ];

      divisionlineline(rline,*(hen0+ihen),&dline);

      vd.dc[EX]=dline.ends[1].d[EX]-dline.ends[0].d[EX];
      vd.dc[EY]=dline.ends[1].d[EY]-dline.ends[0].d[EY];
      vd.dc[EZ]=dline.ends[1].d[EZ]-dline.ends[0].d[EZ];
      if(lengthvector(vd)<eps) /*CASE 2 HENS MAKE 180 DEGREE.*/
      {
        vd.dc[EX]=(hen0+ihen)->ends[1].d[EX]
                 -(hen0+ihen)->ends[0].d[EX];
        vd.dc[EY]=(hen0+ihen)->ends[1].d[EY]
                 -(hen0+ihen)->ends[0].d[EY];
        vd.dc[EZ]=(hen0+ihen)->ends[1].d[EZ]
                 -(hen0+ihen)->ends[0].d[EZ];

        vz.dc[EX]=0.0;
        vz.dc[EX]=0.0;
        vz.dc[EZ]=1.0;

        vo=outerproduct(vz,vd);

        dline.ends[1].d[EX]=dline.ends[0].d[EX]+vo.dc[EX];
        dline.ends[1].d[EY]=dline.ends[0].d[EY]+vo.dc[EY];
        dline.ends[1].d[EZ]=dline.ends[0].d[EZ]+vo.dc[EZ];
      }

      *(hen1+jj)=dline; /*DIVIDING LINES.*//*CASE DEGREE 180 ?*/

      /*drawlocalline(hdc,*vp,tdrccos,*(*(gban.nods+0)),
                    dline.ends[0],dline.ends[1]);
      overlayhdc(*(wdraw.childs+1),SRCPAINT);
      MessageBox(NULL,"Hen1.","Division",MB_OK);*/
    }

    idot=0;
    for(jhen=0;jhen<(ndot-2);jhen++) /*DIVIDING DOTS.*/
    {
      for(khen=(jhen+1);khen<(ndot-1);khen++)
      {
		ifound=distancelineline(*(hen1+jhen),*(hen1+khen),
                                &param1,&param2,&distance,NULL);

        (dot1+idot)->n.code=ifound;          /*0 IF PARALLEL*/
        if(ifound)
        {
          (dot1+idot)->n.d[EX]=(1.0-param1)
                               *((hen1+jhen)->ends[0].d[EX])
                               +param1
                               *((hen1+jhen)->ends[1].d[EX]);
          (dot1+idot)->n.d[EY]=(1.0-param1)
                               *((hen1+jhen)->ends[0].d[EY])
                               +param1
                               *((hen1+jhen)->ends[1].d[EY]);
          (dot1+idot)->n.d[EZ]=(1.0-param1)
                               *((hen1+jhen)->ends[0].d[EZ])
                               +param1
                               *((hen1+jhen)->ends[1].d[EZ]);

            /*drawlocalline(hdc,*vp,tdrccos,*(*(gban.nods+0)),
                          (hen1+jhen)->ends[0],(dot1+idot)->n);
            overlayhdc(*(wdraw.childs+1),SRCPAINT);
            MessageBox(NULL,"Dot1.","Division",MB_OK);*/
        }
        idot++;
      }
    }
    if(idot!=dall)
    {
      MessageBox(NULL,"Dot Missmatch.","Division",MB_OK);
      return;
    }

    /*DICIDING DIVIDED BANS.*/
    idot=ihen;
    if(ihen==ndot-1) jdot=0;
    else             jdot=idot+1;

    bdot=2; /*1st AND 2nd DOTS=ENDS ON PARENT HEN OF BAN.*/
    (dban+ihen)->nnod=bdot;
    (dban+ihen)->nods=(struct onode **)
                      malloc(bdot*sizeof(struct onode *));

    if((dban+ihen)->nods==NULL)
    {
      MessageBox(NULL,"Malloc Failed.","Division",MB_OK);
      return;
    }
    *((dban+ihen)->nods+0)=(lnods+idot); /*1st DOT.*/
    *((dban+ihen)->nods+1)=(lnods+jdot); /*2nd DOT.*/

    /*SEARCH DOTS.*/
    for(kdot=0;kdot<dall;kdot++)
    {
      if((dot1+kdot)->n.code) /*IF EXISTENT CODE=1.*/
      {
        ifound=1;
        for(jhen=0;jhen<(ndot-1);jhen++)
        {
          x0=(hen1+jhen)->ends[0].d[EX];
          y0=(hen1+jhen)->ends[0].d[EY];
          x1=(hen1+jhen)->ends[1].d[EX];
          y1=(hen1+jhen)->ends[1].d[EY];

          x=(lnods+idot)->d[EX];
          y=(lnods+idot)->d[EY];
          sub0=(x-x0)*(y1-y0)-(y-y0)*(x1-x0); /*SUBSTITUTION*/

          if(sub0==0.0)
          {
            x=(lnods+jdot)->d[EX];
            y=(lnods+jdot)->d[EY];
            sub0=(x-x0)*(y1-y0)-(y-y0)*(x1-x0); /*SUBSTITUTION*/
          }
          if(sub0==0.0)
          {
            /*MessageBox(NULL,"Failed.","Division",MB_OK);*/
            ifound=0;
            break;
          }

          x=(dot1+kdot)->n.d[EX];
          y=(dot1+kdot)->n.d[EY];
          sub1=(x-x0)*(y1-y0)-(y-y0)*(x1-x0); /*SUBSTITUTION*/

          if(sub0*sub1<(-eps)) ifound=0; /*DOT ON ANOTHER SIDE*/
        }

        if(ifound)
        {
          bdot++;
          (dban+ihen)->nnod=bdot;
          (dban+ihen)->nods=(struct onode **)
                            realloc((dban+ihen)->nods,
                            bdot*sizeof(struct onode *));
          if((dban+ihen)->nods==NULL)
          {
            MessageBox(NULL,"Realloc Failed.","Division",MB_OK);
            return;
          }

          *((dban+ihen)->nods+(bdot-1))=&((dot1+kdot)->n);

          /*drawlocalline(hdc,*vp,tdrccos,*(*(gban.nods+0)),
                        *(*((dban+ihen)->nods+(bdot-2))),
                        *(*((dban+ihen)->nods+(bdot-1))));
          overlayhdc(*(wdraw.childs+1),SRCPAINT);
          MessageBox(NULL,"Dot1.","Division",MB_OK);*/
        }
      }
    }

    if(cmqbans==NULL) /*ONLY DRAW.*/
    {
      if(hdc!=NULL && vp->vflag.ev.cmqline)
      {
        drawlocalban(hdc,*vp,tdrccos,
                     *(*(gban.nods+0)),*(dban+idot));
      }

      /*cmqpolygon(*(dban+ihen),1.0,&cmqsum);*/
    }
    if(cmqbans!=NULL) /*GIVE DIVIDED BANS.*/
    {
      *(cmqbans+ihen)=*(dban+ihen);

      bdot=(dban+ihen)->nnod;
      (cmqbans+ihen)->nods=(struct onode **)
                           malloc(bdot*sizeof(struct onode *));

      for(jdot=0;jdot<bdot;jdot++)
      {
        *((cmqbans+ihen)->nods+jdot)=(struct onode *)
                                     malloc(sizeof(struct onode));
        *(*((cmqbans+ihen)->nods+jdot))=*(*((dban+ihen)->nods+jdot));
      }
    }
  }

  free(ban.nods);
  free(lnods);
  free(hen0);
  free(hen1);
  free(dot1);
  for(idot=0;idot<ndot;idot++) free((dban+idot)->nods);
  free(dban);

  for(idot=2;idot>=0;idot--) free(*(drccos+idot));
  free(drccos);
  for(idot=2;idot>=0;idot--) free(*(tdrccos+idot));
  free(tdrccos);

  return;
}/*slabdivision2*/

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
  if(laps!=NULL)*laps=(int)strtol(data,NULL,10);

  GetDlgItemText(hdwnd,ID_SAFETY,data,80);
  if(dsafety!=NULL)*dsafety=strtod(data,NULL);

  return;
}/*getincrement*/


void getlaps(HWND hdwnd,int *laps)
/*GET INCREMENT PARAMETERS FROM DIALOG.*/
{
  char data[80];

  GetDlgItemText(hdwnd,ID_LAPS,data,80);
  *laps=(int)strtol(data,NULL,10);

  return;
}/*getincrement*/

void setincrement(HWND hdwnd,
                  int laps,int lap,
                  double dsafes,double dsafe)
/*SET CURRENT INCREMENT PARAMETERS INTO DIALOG.*/
{
  char str[256];

  sprintf(str,"%d:%d",lap,laps);
  SetDlgItemText(hdwnd,ID_LAPS,str);
  SendDlgItemMessage(hdwnd,ID_LAPS,WM_PAINT,0,0);
  if(dsafes!=NULL)         //UJIOKA FOR ONLY LAPS
  {
  sprintf(str,"%.3f:%.3f",dsafe,dsafes);
  SetDlgItemText(hdwnd,ID_SAFETY,str);
  SendDlgItemMessage(hdwnd,ID_SAFETY,WM_PAINT,0,0);
  }
  return;
}/*setincrement*/

void setincrementII(HWND hdwnd,
                  int laps,int lap,
                  double dsafes,double dsafe)
/*SET CURRENT INCREMENT PARAMETERS INTO DIALOG.*/
{
  char str[256];

  sprintf(str,"%d:%d",lap,laps);
  SetDlgItemText(hdwnd,ID_LAPS,str);
  SendDlgItemMessage(hdwnd,ID_LAPS,WM_PAINT,0,0);
  if(dsafes!=NULL)         //UJIOKA FOR ONLY LAPS
  {
  sprintf(str,"%.6f:%.6f",dsafe,dsafes);
  SetDlgItemText(hdwnd,ID_SAFETY,str);
  SendDlgItemMessage(hdwnd,ID_SAFETY,WM_PAINT,0,0);
  }
  return;
}/*setincrement*/

int decreaseband(long int *moff, /*MOVED POSITION OF NODE*/
                 long int *noff, /*CHANGED NODE OFFSETS*/
                 struct arclmframe *af)
/*DECREASE BAND WIDTH.*/
{
  char *cmtx,c; /*CONNECTION MATRIX*/
  char flag,str[512];
  /*char s[80];*/
  long int h,i,j,k,l,m,n,mn,hk,hm,nnode,nelem,off,band;
  long int *ii; /*POSITION OF DIAGONAL*/
  long int mm,nn;

  ///////////////////////////////////////////////////////////araki//
  long int step;
  long int s,t,st;
  long int exchangedband;
  step=0;
  exchangedband=0;
  /*FILE *fband;
  fband=fopen("band.txt","w");  if(fband==NULL) return 0;
  fprintf(fband,"DECREASE BAND WIDTH\n");
  fprintf(fband,"\n");*/
  //////////////////////////////////////////////////////////////////

  nnode=af->nnode;
  nelem=af->nelem;

  /*noff=(long int *)malloc(nnode*sizeof(long int));*/
  /*if(noff==NULL) return 0;*/
  /*for(i=0;i<nnode;i++) *(noff+i)=i;*/

/*
sprintf(str,"NOFF=");
for(mm=0;mm<10;mm++)
{
  sprintf(s," %ld",*(noff+mm));
  strcat(str,s);
}
MessageBox(NULL,str,"Band",MB_OK);
*/

  ii=(long int *)malloc(nnode*sizeof(long int));
  if(ii==NULL) return 0;
  for(i=0;i<nnode;i++) *(ii+i)=((i+1)*(i+2))/2-1;

  mn=(nnode*(nnode+1))/2;
  cmtx=(char *)malloc(mn*sizeof(char));
  if(cmtx==NULL) return 0;
  /*for(i=0;i<mn;i++)
  {
    *(cmtx+mn)=(char)0;
  }*/
  for(mm=0;mm<nnode;mm++)
  {
    for(nn=0;nn<=mm;nn++) *(cmtx+*(ii+mm)-(mm-nn))=(char)0;
  }

  band=0;
  for(i=0;i<nelem;i++) /*ARRANGE CONNECTIONS.*/
  {
    m=(af->elems+i)->node[0]->loff;
    n=(af->elems+i)->node[1]->loff;

    if(m<n){l=m; m=n; n=l;} /*EXCHANGE*/

    mn=(*(ii+m))-(m-n);
    *(cmtx+mn)=(char)1;

	if(band<abs(int(m-n))) band=abs(int(m-n));
  }
  sprintf(str,"Initial Band=%ld",band);
  MessageBox(NULL,str,"Band",MB_OK);

  ///////////////////////////////////////////////////////////araki//
  /*fprintf(fband,"INITIAL\n");*/    /*OUTPUT INITIAL BAND MATRIX*/
  /*for(i=0;i<nnode;i++)
  {
    for(j=0;j<=i;j++)
    {
      mn=*(ii+i)-(i-j);
      fprintf(fband," %d",(int)(*(cmtx+mn)));
    }
    fprintf(fband,"\n");
  }
  fprintf(fband,"\n");
  fprintf(fband,"Initial Band=%ld\n",band);
  fprintf(fband,"\n");*/
  //////////////////////////////////////////////////////////////////

/*
sprintf(str," 0 1 2 3 4 5 6 7 8 9\n");
for(mm=0;mm<10;mm++)
{
  for(nn=0;nn<=mm;nn++)
  {
    sprintf(s," %d",(int)(*(cmtx+(*(ii+mm))-(mm-nn))));
    strcat(str,s);
  }
  strcat(str,"\n");
}
MessageBox(NULL,str,"Band",MB_OK);
*/

  flag=0;
  for(i=nnode-1;i>0;i--) /*SEARCH OUTLAWS.*/
  {
    for(j=i;j<nnode;j++)
    {
      n=j-i;
      m=i+n;
      mn=*(ii+m)-(m-n);

      if(*(cmtx+mn)==1) /*FIND.*/
      {
/*        for(k=n+1;k<m-1;k++) *//*SEARCH CHANGABLE LINE.*/
        for(k=n+1;k<=m-1;k++) /*SEARCH CHANGABLE LINE.*/
        {
          flag=0;
          for(l=0;l<nnode;l++) /*0-n,m-N*/
          {
            if(l<=n) mn=*(ii+k)-(k-l);      //araki
/*            else     mn=*(ii+l)-(l-k); */
            else     mn=*(ii+l)-(l-m);   //araki
            if(*(cmtx+mn)==1)
            {
              flag=1;
              break;
            }
/*            if(l==n) l=m-1;     */
///////////////////////////////////////////////////araki//
            if(l==n)
            {
              /*if(k<=n) l=i+k-1;
              else*/     l=i+k;
            }
//////////////////////////////////////////////////////////
          }
          if(flag==0) /*EXCHANGE k,m.*/
          {
            /*EXCHANGE NODE OFFSETS.*/
/*
sprintf(str,"Exchange k=%ld m=%ld",k,m);
MessageBox(NULL,str,"Band",MB_OK);
*/
            off=*(noff+k);
            *(noff+k)=*(noff+m);
            *(noff+m)=off;

            for(h=0;h<nnode;h++) /*EXCHANGE LINES.*/
            {
              if(h<k)
              {
                hk=*(ii+k)-(k-h);
                hm=*(ii+m)-(m-h);
              }
              else if(h==k)
              {
                hk=*(ii+k);
                hm=*(ii+m);
              }
              else if(h<m)
              {
                hk=*(ii+h)-(h-k);
                hm=*(ii+m)-(m-h);
              }
              else if(m<h)
              {
                hk=*(ii+h)-(h-k);
                hm=*(ii+h)-(h-m);
              }

              if(h!=m)
              {
                /*EXCHANGE.*/
                c=*(cmtx+hk);
                *(cmtx+hk)=*(cmtx+hm);
                *(cmtx+hm)=c;
              }
            }
            ///////////////////////////////////////////////////////////araki//
/*            step++;
            fprintf(fband,"%dSTEP\n",step);
            fprintf(fband,"i=%ld, ",i);
            fprintf(fband,"m=%ld, ",m);
            fprintf(fband,"n=%ld, ",n);
            fprintf(fband,"k=%ld\n",k);
            fprintf(fband,"%ld行目を%ld行目と入れ替えた結果\n",m,k);
            for(s=0;s<nnode;s++)
            {
              for(t=0;t<=s;t++)
              {
                st=*(ii+s)-(s-t);
                fprintf(fband," %d",(int)(*(cmtx+st)));
              }
              fprintf(fband,"\n");
            }
            fprintf(fband,"\n");   */
            //////////////////////////////////////////////////////////////////

/*
sprintf(str,"\0");
for(mm=0;mm<10;mm++)
{
  sprintf(s," %ld",*(noff+mm));
  strcat(str,s);
}
strcat(str,"\n");
for(mm=0;mm<10;mm++)
{
  for(nn=0;nn<=mm;nn++)
  {
    sprintf(s," %d",*(cmtx+*(ii+mm)-(mm-nn)));
    strcat(str,s);
  }
  strcat(str,"\n");
}
MessageBox(NULL,str,"Band",MB_OK);
*/
            break;
          }
        }
/*        if(flag==1) break; */ //araki
          if(flag==1&&exchangedband<i) exchangedband=i;
      }
    }
    currentpivot(i+1,nnode);
/*    if(flag==1) break;  */  //araki
  }

  /*ARRANGE MOVED POSITION OF NODE.*/
  for(mm=0;mm<nnode;mm++) *(moff+*(noff+mm))=mm;

  ///////////////////////////////////////////////////////////araki//
  /*fprintf(fband,"FINAL\n");*/  /*OUTPUT EXCHANGED BAND MATRIX*/
  /*for(s=0;s<nnode;s++)
  {
    for(t=0;t<=s;t++)
    {
      mn=*(ii+s)-(s-t);
      fprintf(fband," %d",(int)(*(cmtx+mn)));
    }
    fprintf(fband,"\n");
  }
  fprintf(fband,"\n");
  fprintf(fband,"Exchanged Band=%ld",exchangedband);
  fclose(fband);*/
  //////////////////////////////////////////////////////////////////

  free(cmtx);
  free(ii);

/*  sprintf(str,"Exchanged Band=%ld",i);  */
  sprintf(str,"Exchanged Band=%ld",exchangedband);   //araki
  MessageBox(NULL,str,"Band",MB_OK);

/*
sprintf(str,"Last NOFF=");
for(mm=0;mm<10;mm++)
{
  sprintf(s," %ld",*(noff+mm));
  strcat(str,s);
}
MessageBox(NULL,str,"Band",MB_OK);
*/

  return i; /*BAND WIDTH.*/
}/*decreaseband*/

#if 0
int decreaseband(long int *moff, /*MOVED POSITION OF NODE*/
                 long int *noff, /*CHANGED NODE OFFSETS*/
                 struct arclmframe *af)
/*DECREASE BAND WIDTH.*/
{
  char *cmtx,c; /*CONNECTION MATRIX*/
  char flag,str[512];
  /*char s[80];*/
  long int h,i,j,k,l,m,n,mn,hk,hm,nnode,nelem,off,band;
  long int *ii; /*POSITION OF DIAGONAL*/
  long int mm,nn;

  nnode=af->nnode;
  nelem=af->nelem;

  /*noff=(long int *)malloc(nnode*sizeof(long int));*/
  /*if(noff==NULL) return 0;*/
  /*for(i=0;i<nnode;i++) *(noff+i)=i;*/

/*
sprintf(str,"NOFF=");
for(mm=0;mm<10;mm++)
{
  sprintf(s," %ld",*(noff+mm));
  strcat(str,s);
}
MessageBox(NULL,str,"Band",MB_OK);
*/

  ii=(long int *)malloc(nnode*sizeof(long int));
  if(ii==NULL) return 0;
  for(i=0;i<nnode;i++) *(ii+i)=((i+1)*(i+2))/2-1;

  mn=(nnode*(nnode+1))/2;
  cmtx=(char *)malloc(mn*sizeof(char));
  if(cmtx==NULL) return 0;
  /*for(i=0;i<mn;i++)
  {
    *(cmtx+mn)=(char)0;
  }*/
  for(mm=0;mm<nnode;mm++)
  {
    for(nn=0;nn<=mm;nn++) *(cmtx+*(ii+mm)-(mm-nn))=(char)0;
  }

  band=0;
  for(i=0;i<nelem;i++) /*ARRANGE CONNECTIONS.*/
  {
    m=(af->elems+i)->node[0]->loff;
    n=(af->elems+i)->node[1]->loff;

    if(m<n){l=m; m=n; n=l;} /*EXCHANGE*/

    mn=(*(ii+m))-(m-n);
    *(cmtx+mn)=(char)1;

    if(band<abs(m-n)) band=abs(m-n);
  }
  sprintf(str,"Initial Band=%ld",band);
  MessageBox(NULL,str,"Band",MB_OK);

/*
sprintf(str," 0 1 2 3 4 5 6 7 8 9\n");
for(mm=0;mm<10;mm++)
{
  for(nn=0;nn<=mm;nn++)
  {
    sprintf(s," %d",(int)(*(cmtx+(*(ii+mm))-(mm-nn))));
    strcat(str,s);
  }
  strcat(str,"\n");
}
MessageBox(NULL,str,"Band",MB_OK);
*/

  flag=0;
  for(i=nnode-1;i>0;i--) /*SEARCH OUTLAWS.*/
  {
    for(j=i;j<nnode;j++)
    {
      n=j-i;
      m=i+n;
      mn=*(ii+m)-(m-n);

      if(*(cmtx+mn)==1) /*FIND.*/
      {
        for(k=n+1;k<m-1;k++) /*SEARCH CHANGABLE LINE.*/
        {
          flag=0;
          for(l=0;l<nnode;l++) /*0-n,m-N*/
          {
            if(l<=n) mn=*(ii+k)-(k-l);
            else     mn=*(ii+l)-(l-k);
            if(*(cmtx+mn)==1)
            {
              flag=1;
              break;
            }
            if(l==n) l=m-1;
          }
          if(flag==0) /*EXCHANGE k,m.*/
          {
            /*EXCHANGE NODE OFFSETS.*/
/*
sprintf(str,"Exchange k=%ld m=%ld",k,m);
MessageBox(NULL,str,"Band",MB_OK);
*/
            off=*(noff+k);
            *(noff+k)=*(noff+m);
            *(noff+m)=off;

            for(h=0;h<nnode;h++) /*EXCHANGE LINES.*/
            {
              if(h<k)
              {
                hk=*(ii+k)-(k-h);
                hm=*(ii+m)-(m-h);
              }
              else if(h==k)
              {
                hk=*(ii+k);
                hm=*(ii+m);
              }
              else if(h<m)
              {
                hk=*(ii+h)-(h-k);
                hm=*(ii+m)-(m-h);
              }
              else if(m<h)
              {
                hk=*(ii+h)-(h-k);
                hm=*(ii+h)-(h-m);
              }

              if(h!=m)
              {
                /*EXCHANGE.*/
                c=*(cmtx+hk);
                *(cmtx+hk)=*(cmtx+hm);
                *(cmtx+hm)=c;
              }
            }
/*
sprintf(str,"\0");
for(mm=0;mm<10;mm++)
{
  sprintf(s," %ld",*(noff+mm));
  strcat(str,s);
}
strcat(str,"\n");
for(mm=0;mm<10;mm++)
{
  for(nn=0;nn<=mm;nn++)
  {
    sprintf(s," %d",*(cmtx+*(ii+mm)-(mm-nn)));
    strcat(str,s);
  }
  strcat(str,"\n");
}
MessageBox(NULL,str,"Band",MB_OK);
*/
            break;
          }
        }
        if(flag==1) break;
      }
    }
    currentpivot(i+1,nnode);
    if(flag==1) break;
  }

  /*ARRANGE MOVED POSITION OF NODE.*/
  for(mm=0;mm<nnode;mm++) *(moff+*(noff+mm))=mm;

  free(cmtx);
  free(ii);

  sprintf(str,"Exchanged Band=%ld",i);
  MessageBox(NULL,str,"Band",MB_OK);

/*
sprintf(str,"Last NOFF=");
for(mm=0;mm<10;mm++)
{
  sprintf(s," %ld",*(noff+mm));
  strcat(str,s);
}
MessageBox(NULL,str,"Band",MB_OK);
*/

  return i; /*BAND WIDTH.*/
}/*decreaseband*/
#endif

int exchangelines(struct gcomponent *gmtx1,
                  double *gvct1,struct oconf *confs1,
                  struct gcomponent *gmtx2,
                  double *gvct2,struct oconf *confs2,
                  long int *moff,long int *noff,long int nnode,
                  char flag)
/*EXCHANGE LINES OF GLOBAL MATRIX,VECTOR,CONFS.*/
/*EXCHANGE:FLAG=0 RETURN:FLAG=1*/
{
  char str[4096],s[80];
  long int i,j,ii,k,l,kk,m,n,msize;
  long int mm,nn;
  double /**gv,*/data;
  /*struct gcomponent *gm;*/
  /*struct gcomponent ginit={0,0,0.0,NULL};*/
  /*struct oconf *oc;*/

  msize=6*nnode;

/*
  gm=(struct gcomponent *)malloc(msize*sizeof(struct gcomponent));
  gv=(double *)malloc(msize*sizeof(double));
  oc=(struct oconf *)malloc(msize*sizeof(struct oconf));
  if(gm==NULL || gv==NULL || oc==NULL) return 0;
  for(i=0;i<msize;i++)
  {
    ginit.m=(unsigned short int)(i+1);
    *(gm+i)=ginit;

    *(gv+i)=0.0;

    (oc+i)->iconf=0;
    (oc+i)->value=0.0;
  }
  comps=msize;
*/
/*
sprintf(str,"NOFF=");
for(mm=0;mm<10;mm++)
{
  sprintf(s," %ld",*(noff+mm));
  strcat(str,s);
}
MessageBox(NULL,str,"Exchange",MB_OK);
*/

  for(i=0;i<nnode;i++)
  {
    if(flag==0) m=*(noff+i);
    if(flag==1) m=*(moff+i);

    for(j=0;j<6;j++)
    {
      ii=6*i+j;
      mm=6*m+j;

      *(gvct2+ii)=*(gvct1+mm);
      *(confs2+ii)=*(confs1+mm);

      for(k=0;k<=i;k++)
      {
        if(flag==0) n=*(noff+k);
        if(flag==1) n=*(moff+k);

        for(l=0;l<6;l++)
        {
          kk=6*k+l;
          nn=6*n+l;

          if(ii>=kk)
          {
            gread(gmtx1,mm+1,nn+1,&data);
            if(data!=0.0) gwrite(gmtx2,ii+1,kk+1,data);
          }
        }
      }
      currentpivot(ii+1,msize);
    }
  }

/*
sprintf(str,"Global Matrix\n");
for(mm=0;mm<12;mm++)
{
  for(nn=0;nn<=mm;nn++)
  {
    gread(gmtx,mm+1,nn+1,&data);
    sprintf(s," %.1E",data);
    strcat(str,s);
  }
  strcat(str,"\n");
}
MessageBox(NULL,str,"Exchange",MB_OK);

sprintf(str,"Exchanged Global Matrix\n");
for(mm=0;mm<12;mm++)
{
  for(nn=0;nn<=mm;nn++)
  {
    gread(gm,mm+1,nn+1,&data);
    sprintf(s," %.1E",data);
    strcat(str,s);
  }
  strcat(str,"\n");
}
MessageBox(NULL,str,"Exchange",MB_OK);
*/

  /*gfree(gmtx,nnode);*/
  /*free(gvct);*/
  /*free(confs);*/

/*
  gmtx=gm;
  gvct=gv;
  confs=oc;
*/

  return 1;
}/*exchangelines*/

int exchangelinesII(double **kmtx1,double **gmtx1,struct oconf *confs1,
                  double **kmtx2,double **gmtx2,struct oconf *confs2,
                  long int *moff,long int *noff,long int nnode,
                  char flag)
/*EXCHANGE LINES OF GLOBAL MATRIX,CONFS.*/
/*EXCHANGE:FLAG=0 RETURN:FLAG=1*/
/*UJIOKA FOR BCLNG003*/
{
  char str[4096],s[80];
  long int i,j,ii,k,l,kk,m,n,msize;
  long int mm,nn;
  double /**gv,*/data;
  /*struct gcomponent *gm;*/
  /*struct gcomponent ginit={0,0,0.0,NULL};*/
  /*struct oconf *oc;*/

  msize=6*nnode;

/*
  gm=(struct gcomponent *)malloc(msize*sizeof(struct gcomponent));
  gv=(double *)malloc(msize*sizeof(double));
  oc=(struct oconf *)malloc(msize*sizeof(struct oconf));
  if(gm==NULL || gv==NULL || oc==NULL) return 0;
  for(i=0;i<msize;i++)
  {
    ginit.m=(unsigned short int)(i+1);
    *(gm+i)=ginit;

    *(gv+i)=0.0;

    (oc+i)->iconf=0;
    (oc+i)->value=0.0;
  }
  comps=msize;
*/
/*
sprintf(str,"NOFF=");
for(mm=0;mm<10;mm++)
{
  sprintf(s," %ld",*(noff+mm));
  strcat(str,s);
}
MessageBox(NULL,str,"Exchange",MB_OK);
*/

  for(i=0;i<nnode;i++)
  {
    if(flag==0) m=*(noff+i);
    if(flag==1) m=*(moff+i);

    for(j=0;j<6;j++)
    {
      ii=6*i+j;
      mm=6*m+j;

/*
      *(gvct2+ii)=*(gvct1+mm);
*/
      *(confs2+ii)=*(confs1+mm);

      for(k=0;k<=i;k++)
      {
        if(flag==0) n=*(noff+k);
        if(flag==1) n=*(moff+k);

        for(l=0;l<6;l++)
        {
          kk=6*k+l;
          nn=6*n+l;

          if(ii>=kk)
          {
            *(*(kmtx2+ii)+kk)=*(*(kmtx1+mm)+nn);
            *(*(gmtx2+ii)+kk)=*(*(gmtx1+mm)+nn);
          /*
            gread(gmtx1,mm+1,nn+1,&data);
            if(data!=0.0) gwrite(gmtx2,ii+1,kk+1,data);
          */
          }
        }
      }
      currentpivot(ii+1,msize);
    }
  }

/*
sprintf(str,"Global Matrix\n");
for(mm=0;mm<12;mm++)
{
  for(nn=0;nn<=mm;nn++)
  {
    gread(gmtx,mm+1,nn+1,&data);
    sprintf(s," %.1E",data);
    strcat(str,s);
  }
  strcat(str,"\n");
}
MessageBox(NULL,str,"Exchange",MB_OK);

sprintf(str,"Exchanged Global Matrix\n");
for(mm=0;mm<12;mm++)
{
  for(nn=0;nn<=mm;nn++)
  {
    gread(gm,mm+1,nn+1,&data);
    sprintf(s," %.1E",data);
    strcat(str,s);
  }
  strcat(str,"\n");
}
MessageBox(NULL,str,"Exchange",MB_OK);
*/

  /*gfree(gmtx,nnode);*/
  /*free(gvct);*/
  /*free(confs);*/

/*
  gmtx=gm;
  gvct=gv;
  confs=oc;
*/

  return 1;
}/*exchangelinesII*/

int exchangelinesIIfloat(float **kmtx1,float **gmtx1,struct oconf *confs1,
                  float **kmtx2,float **gmtx2,struct oconf *confs2,
                  long int *moff,long int *noff,long int nnode,
                  char flag)
/*EXCHANGE LINES OF GLOBAL MATRIX,CONFS.*/
/*EXCHANGE:FLAG=0 RETURN:FLAG=1*/
/*UJIOKA FOR BCLNG003*/
{
  char str[4096],s[80];
  long int i,j,ii,k,l,kk,m,n,msize;
  long int mm,nn;
  double /**gv,*/data;
  /*struct gcomponent *gm;*/
  /*struct gcomponent ginit={0,0,0.0,NULL};*/
  /*struct oconf *oc;*/

  msize=6*nnode;

/*
  gm=(struct gcomponent *)malloc(msize*sizeof(struct gcomponent));
  gv=(double *)malloc(msize*sizeof(double));
  oc=(struct oconf *)malloc(msize*sizeof(struct oconf));
  if(gm==NULL || gv==NULL || oc==NULL) return 0;
  for(i=0;i<msize;i++)
  {
    ginit.m=(unsigned short int)(i+1);
    *(gm+i)=ginit;

    *(gv+i)=0.0;

    (oc+i)->iconf=0;
    (oc+i)->value=0.0;
  }
  comps=msize;
*/
/*
sprintf(str,"NOFF=");
for(mm=0;mm<10;mm++)
{
  sprintf(s," %ld",*(noff+mm));
  strcat(str,s);
}
MessageBox(NULL,str,"Exchange",MB_OK);
*/

  for(i=0;i<nnode;i++)
  {
    if(flag==0) m=*(noff+i);
    if(flag==1) m=*(moff+i);

    for(j=0;j<6;j++)
    {
      ii=6*i+j;
      mm=6*m+j;

/*
      *(gvct2+ii)=*(gvct1+mm);
*/
      *(confs2+ii)=*(confs1+mm);

      for(k=0;k<=i;k++)
      {
        if(flag==0) n=*(noff+k);
        if(flag==1) n=*(moff+k);

        for(l=0;l<6;l++)
        {
          kk=6*k+l;
          nn=6*n+l;

          if(ii>=kk)
          {
            *(*(kmtx2+ii)+kk)=*(*(kmtx1+mm)+nn);
            *(*(gmtx2+ii)+kk)=*(*(gmtx1+mm)+nn);
          /*
            gread(gmtx1,mm+1,nn+1,&data);
            if(data!=0.0) gwrite(gmtx2,ii+1,kk+1,data);
          */
          }
        }
      }
      currentpivot(ii+1,msize);
    }
  }

/*
sprintf(str,"Global Matrix\n");
for(mm=0;mm<12;mm++)
{
  for(nn=0;nn<=mm;nn++)
  {
    gread(gmtx,mm+1,nn+1,&data);
    sprintf(s," %.1E",data);
    strcat(str,s);
  }
  strcat(str,"\n");
}
MessageBox(NULL,str,"Exchange",MB_OK);

sprintf(str,"Exchanged Global Matrix\n");
for(mm=0;mm<12;mm++)
{
  for(nn=0;nn<=mm;nn++)
  {
	gread(gm,mm+1,nn+1,&data);
    sprintf(s," %.1E",data);
    strcat(str,s);
  }
  strcat(str,"\n");
}
MessageBox(NULL,str,"Exchange",MB_OK);
*/

  /*gfree(gmtx,nnode);*/
  /*free(gvct);*/
  /*free(confs);*/

/*
  gmtx=gm;
  gvct=gv;
  confs=oc;
*/

  return 1;
}/*exchangelinesII*/

/* 150515 fukushima for all */
int croutlu(struct gcomponent *gmtx,
			struct oconf *confs,
			long int msize,
			double *det,double *sign,
			struct gcomponent *gcomp1)
{
  /*char iconf;*/ /*0:FREE 1:FIXED*/
  char str[256];
  long int i,j,k;
  /*double det=1.0;*/
  double data1;
  struct gcomponent *pivot,*pcomp;
  struct gcomponent *gcomp2,*gcomp3,*gcomp4;
  double pivotgcomp1;

  *det=0.0;
  *sign=0.0;

  for(j=1;j<=msize;j++)                             /*DECOMPOSITION.*/
  {
	if((confs+j-1)->iconf==0) /*FREE*/
	{
	  pivot=(gmtx+(j-1)); /*PIVOT.*/

	  if((pivot->value) == 0.0){*sign = -1; return (j-1);}  /*INSTABLE.*/
	  if((pivot->value) < 0.0) {*sign += 1; }
	  /*det*=pivot->value;*/
	  *det+=log10(fabs(pivot->value));   /*LOG BY 10 OF DETERMINANT.*/
	  //*sign*=pivot->value/fabs(pivot->value); /*SIGN OF DETERMINANT.*/

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
		pivotgcomp1=(pivot->value)*(gcomp1->value);

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
				  sprintf(str,"GCOMP4-1: %d", i);
				  errormessage(str);
                  errormessage("CROUT:MEMORY INSUFFICIENT.");
				  *sign=-999;
				  return 0;
				}
				gcomp3->down=gcomp4;
				gcomp4->m=(unsigned short int)k;
                /*gcomp4->n=(unsigned short int)i;*/
                gcomp4->down=NULL;
				gcomp4->value=-pivotgcomp1
                             *(gcomp2->value);
                gcomp3=gcomp4;

                comps++;
              }
              else if((gcomp3->m)==k)
              {
				gcomp3->value-=pivotgcomp1
                              *(gcomp2->value);
              }
			  else if((gcomp3->m)>k) /*FILL*/
              {
                gcomp4=(struct gcomponent *)
                       malloc(sizeof(struct gcomponent));
				if(gcomp4==NULL)
				{
                  sprintf(str,"GCOMP4-2: %d", i);
                  errormessage(str);
                  errormessage("CROUT:MEMORY INSUFFICIENT.");
				  *sign=-999;
                  return 0;
                }

                pcomp->down=gcomp4;
                gcomp4->m=(unsigned short int)k;
                /*gcomp4->n=(unsigned short int)i;*/
                gcomp4->down=gcomp3;
				gcomp4->value=-pivotgcomp1
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
	  //currentpivot(j,msize);
    }
  }
  if(*sign<0.0) return 0;
  else return 1;
}/* croutlu */
/***/

/* 150515 fukushima for all */
int forwardbackward(struct gcomponent *gmtx,
					double *gvct, struct oconf *confs,
					long int msize,
					struct gcomponent *gcomp1)
{
  double data1;
  long int i,j,k;

  //errormessage("FORWARD ELIMINATION.");
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
	//currentpivot(j,msize);
  }

  //errormessage("BACKWARD SUBSTITUTION.");
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
	//currentpivot((j-1),msize);
  }

  return 1;
}/* forwardbackward */
/***/



/* 150515 fukushima for all */
int croutluII(struct gcomponent *gmtx,
			struct oconf *confs,
			long int msize, long int csize,
			double *det,double *sign,
			struct gcomponent *gcomp1)
{
  /*char iconf;*/ /*0:FREE 1:FIXED*/
  char str[256];
  long int i,j,k;
  /*double det=1.0;*/
  double data1;
  struct gcomponent *pivot,*pcomp;
  struct gcomponent *gcomp2,*gcomp3,*gcomp4;
  double pivotgcomp1;

  *det=0.0;
  *sign=0.0;

  for(j=1;j<=msize+csize;j++)                             /*DECOMPOSITION.*/
  {
	if((j<=msize && (confs+j-1)->iconf==0) || j>msize) /*FREE*/
	{
	  pivot=(gmtx+(j-1)); /*PIVOT.*/

	  if(fabs(pivot->value) == 0.0){*sign = -1; return (j-1);}  /*INSTABLE.*/
	  if((pivot->value) < 0.0 && j<=msize) {*sign += 1; }
	  /*det*=pivot->value;*/
	  *det+=log10(fabs(pivot->value));   /*LOG BY 10 OF DETERMINANT.*/
	  //*sign*=pivot->value/fabs(pivot->value); /*SIGN OF DETERMINANT.*/

	  gcomp1=pivot;
	  while(gcomp1->down!=NULL) /*DOWNWARD.*/
	  {
		gcomp1=gcomp1->down;

		i=gcomp1->m; /*i:LINE CODE IN PIVOT ROW.*/
		if((i<=msize && (confs+i-1)->iconf==0) || i>msize) /*FREE*/
		{
          gcomp1->value/=(pivot->value);
		}
	  }

	  gcomp1=pivot;
	  while(gcomp1->down!=NULL) /*DOWNWARD.*/
	  {
		gcomp1=gcomp1->down; /*Aij*/
		pivotgcomp1=(pivot->value)*(gcomp1->value);

		i=gcomp1->m; /*i:LINE CODE IN PIVOT ROW.*/
		if((i<=msize && (confs+i-1)->iconf==0) || i>msize) /*FREE*/
        {
		  gcomp2=gcomp1; /*Akj*/
		  gcomp3=(gmtx+(i-1)); /*Aki*/

		  while(1)
		  {
			k=gcomp2->m;
			if((k<=msize && (confs+k-1)->iconf==0) || k>msize) /*FREE*/
			{
			  if((gcomp3->m)<k) /*ADD*/
			  {
                gcomp4=(struct gcomponent *)
                       malloc(sizeof(struct gcomponent));
				if(gcomp4==NULL)
				{
				  sprintf(str,"GCOMP4-1: %d", i);
				  errormessage(str);
				  errormessage("CROUT:MEMORY INSUFFICIENT.");
				  *sign=-999;
				  return 0;
				}
				gcomp3->down=gcomp4;
				gcomp4->m=(unsigned short int)k;
				/*gcomp4->n=(unsigned short int)i;*/
				gcomp4->down=NULL;
				gcomp4->value=-pivotgcomp1
							 *(gcomp2->value);
				gcomp3=gcomp4;

				comps++;
			  }
			  else if((gcomp3->m)==k)
			  {
				gcomp3->value-=pivotgcomp1
							  *(gcomp2->value);
			  }
			  else if((gcomp3->m)>k) /*FILL*/
              {
                gcomp4=(struct gcomponent *)
                       malloc(sizeof(struct gcomponent));
				if(gcomp4==NULL)
				{
                  sprintf(str,"GCOMP4-2: %d", i);
                  errormessage(str);
                  errormessage("CROUT:MEMORY INSUFFICIENT.");
				  *sign=-999;
                  return 0;
                }

                pcomp->down=gcomp4;
                gcomp4->m=(unsigned short int)k;
                /*gcomp4->n=(unsigned short int)i;*/
                gcomp4->down=gcomp3;
				gcomp4->value=-pivotgcomp1
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
	  //currentpivot(j,msize);
	}
  }
  if(*sign<0.0) return 0;
  else return 1;
}/* croutluII */

/* 150515 fukushima for all */
int forwardbackwardII(struct gcomponent *gmtx,
					double *gvct, struct oconf *confs,
					long int msize, long int csize,
					struct gcomponent *gcomp1)
{
  double data1;
  long int i,j,k;

  //errormessage("FORWARD ELIMINATION.");
  for(j=1;j<=msize+csize;j++)                                  /*FORWARD.*/
  {
	if((j<=msize && (confs+j-1)->iconf==0) || j>msize) /*FREE*/
	{
	  data1=*(gvct+j-1);
	  gcomp1=(gmtx+(j-1)); /*DIAGONAL.*/

	  while(gcomp1->down!=NULL) /*DOWNWARD.*/
	  {
		gcomp1=gcomp1->down;
		i=gcomp1->m;

		if((i<=msize && (confs+i-1)->iconf==0) || i>msize) /*FREE*/
		{
		  *(gvct+i-1)-=gcomp1->value*data1;
		}
	  }
	}
	//currentpivot(j,msize);
  }

  //errormessage("BACKWARD SUBSTITUTION.");
  for(j=msize+csize;j>=1;j--)                                 /*BACKWARD.*/
  {
	if((j<=msize && (confs+j-1)->iconf==0) || j>msize) /*FREE*/
	{
	  data1=*(gvct+j-1);
	  gcomp1=(gmtx+(j-1)); /*DIAGONAL.*/
	  data1/=gcomp1->value;

	  while(gcomp1->down!=NULL) /*DOWNWARD.*/
	  {
        gcomp1=gcomp1->down;
        i=gcomp1->m;

		if((i<=msize && (confs+i-1)->iconf==0) || i>msize) /*FREE*/
		{
          data1-=gcomp1->value*(*(gvct+i-1));
        }
      }
	  *(gvct+j-1)=data1;
	}
	//currentpivot((j-1),msize);
  }

  return 1;
}/* forwardbackwardII */
/***/

int croutludecomposition(struct gcomponent *gmtx,
						 double *gvct,struct oconf *confs,
						 long int msize,
						 double *det,double *sign)
	/*SOLVE SIMULTANEOUS LINEAR EQUATIONS.*/
{
  /*char iconf;*/ /*0:FREE 1:FIXED*/
  long int i,j,k;
  /*double det=1.0;*/
  double data1;
  struct gcomponent *pivot,*pcomp;
  struct gcomponent *gcomp1,*gcomp2,*gcomp3,*gcomp4;
  double pivotgcomp1;

  *det=0.0;
  *sign=0.0;

  for(j=1;j<=msize;j++)                             /*DECOMPOSITION.*/
  {
	if((confs+j-1)->iconf==0) /*FREE*/
	{
	  pivot=(gmtx+(j-1)); /*PIVOT.*/
	  if((pivot->value)==0.0){*sign=-1; return (j-1);}  /*INSTABLE.*/
	  if((pivot->value)<0.0){*sign+=1;}
	  *det+=log10(fabs(pivot->value));   /*LOG BY 10 OF DETERMINANT.*/
	  /**sign*=pivot->value/fabs(pivot->value);*/ /*SIGN OF DETERMINANT.*/


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
		pivotgcomp1=(pivot->value)*(gcomp1->value);

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
				  *sign=-999;
                  return 0;
				}
                gcomp3->down=gcomp4;
                gcomp4->m=(unsigned short int)k;
                /*gcomp4->n=(unsigned short int)i;*/
                gcomp4->down=NULL;
				gcomp4->value=-pivotgcomp1
							 *(gcomp2->value);
                gcomp3=gcomp4;

                comps++;
              }
              else if((gcomp3->m)==k)
              {
				gcomp3->value-=pivotgcomp1
                              *(gcomp2->value);
              }
              else if((gcomp3->m)>k) /*FILL*/
              {
                gcomp4=(struct gcomponent *)
					   malloc(sizeof(struct gcomponent));
                if(gcomp4==NULL)
                {
                  errormessage("CROUT:MEMORY INSUFFICIENT.");
				  *sign=-999;
                  return 0;
                }

                pcomp->down=gcomp4;
                gcomp4->m=(unsigned short int)k;
                /*gcomp4->n=(unsigned short int)i;*/
                gcomp4->down=gcomp3;
				gcomp4->value=-pivotgcomp1
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
  if(*sign<0.0) return 0;

  /*errormessage("FORWARD ELIMINATION.");*/
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

  /*errormessage("BACKWARD SUBSTITUTION.");*/
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

  /*errormessage("CROUT LU COMPLETED.");*/
  return 1;
}/*croutludecomposition*/

int croutludecomposition_arclength(struct gcomponent *gmtx,
						 double *gvct,double *gvct2,struct oconf *confs,
						 long int msize,
						 double *det,double *sign,int iteration)
	/*SOLVE SIMULTANEOUS LINEAR EQUATIONS FOR ARC-LENGTH METHOD.*/
{
  /*char iconf;*/ /*0:FREE 1:FIXED*/
  long int i,j,k;
  /*double det=1.0;*/
  double data1,data2;
  struct gcomponent *pivot,*pcomp;
  struct gcomponent *gcomp1,*gcomp2,*gcomp3,*gcomp4;
  double pivotgcomp1;

  *det=0.0;
  *sign=0.0;

  for(j=1;j<=msize;j++)                             /*DECOMPOSITION.*/
  {
	if((confs+j-1)->iconf==0) /*FREE*/
	{
	  pivot=(gmtx+(j-1)); /*PIVOT.*/
	  if((pivot->value)==0.0){*sign=-1; return (j-1);}  /*INSTABLE.*/
	  if((pivot->value)<0.0){*sign+=1;}
	  /*det*=pivot->value;*/
	  *det+=log10(fabs(pivot->value));   /*LOG BY 10 OF DETERMINANT.*/
	  /**sign*=pivot->value/fabs(pivot->value);*/ /*SIGN OF DETERMINANT.*/


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
		pivotgcomp1=(pivot->value)*(gcomp1->value);

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
				  *sign=-999;
                  return 0;
                }
				gcomp3->down=gcomp4;
                gcomp4->m=(unsigned short int)k;
				/*gcomp4->n=(unsigned short int)i;*/
				gcomp4->down=NULL;
				gcomp4->value=-pivotgcomp1
                             *(gcomp2->value);
				gcomp3=gcomp4;

                comps++;
			  }
              else if((gcomp3->m)==k)
              {
				gcomp3->value-=pivotgcomp1
							  *(gcomp2->value);
              }
              else if((gcomp3->m)>k) /*FILL*/
			  {
                gcomp4=(struct gcomponent *)
                       malloc(sizeof(struct gcomponent));
				if(gcomp4==NULL)
                {
				  errormessage("CROUT:MEMORY INSUFFICIENT.");
				  *sign=-999;
				  return 0;
				}

				pcomp->down=gcomp4;
                gcomp4->m=(unsigned short int)k;
                /*gcomp4->n=(unsigned short int)i;*/
				gcomp4->down=gcomp3;
				gcomp4->value=-pivotgcomp1
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
	  /*currentpivot(j,msize);*/
	}
  }
  if(*sign<0) return 0;

  /*errormessage("FORWARD ELIMINATION.");*/
  for(j=1;j<=msize;j++)                                  /*FORWARD.*/
  {
	if((confs+j-1)->iconf==0) /*FREE*/
	{
	  data1=*(gvct+j-1);
	  if(iteration!=1)data2=*(gvct2+j-1);
	  gcomp1=(gmtx+(j-1)); /*DIAGONAL.*/

	  while(gcomp1->down!=NULL) /*DOWNWARD.*/
	  {
		gcomp1=gcomp1->down;
		i=gcomp1->m;

		if((confs+i-1)->iconf==0) /*FREE*/
		{
		  *(gvct+i-1)-=gcomp1->value*data1;
		  if(iteration!=1)*(gvct2+i-1)-=gcomp1->value*data2;
		}
	  }
	}
	/*currentpivot(j,msize);*/
  }

  /*errormessage("BACKWARD SUBSTITUTION.");*/
  for(j=msize;j>=1;j--)                                 /*BACKWARD.*/
  {
	if((confs+j-1)->iconf==0) /*FREE*/
	{
	  data1=*(gvct+j-1);
	  if(iteration!=1)data2=*(gvct2+j-1);
	  gcomp1=(gmtx+(j-1)); /*DIAGONAL.*/
	  data1/=gcomp1->value;
	  if(iteration!=1)data2/=gcomp1->value;
	  while(gcomp1->down!=NULL) /*DOWNWARD.*/
	  {
		gcomp1=gcomp1->down;
		i=gcomp1->m;

		if((confs+i-1)->iconf==0) /*FREE*/
		{
		  data1-=gcomp1->value*(*(gvct+i-1));
		  if(iteration!=1)data2-=gcomp1->value*(*(gvct2+i-1));
		}
	  }
	  *(gvct+j-1)=data1;
	  if(iteration!=1)*(gvct2+j-1)=data2;
	}
	/*currentpivot((j-1),msize);*/
  }

  /*errormessage("CROUT LU COMPLETED.");*/
  return 1;
}/*croutludecomposition_arclength*/



void currentpivot(long int i,long int msize)
/*WRITE PIVOT LINE INTO DIALOG BOX.*/
{
#if 0
  HWND hdwnd=(wmenu.childs+2)->hwnd;
  char str[20];

  sprintf(str,"%4ld/%4ld",i,msize);
  SetDlgItemText(hdwnd,ID_PIVOT,str);
  SendDlgItemMessage(hdwnd,ID_PIVOT,WM_PAINT,0,0);
#endif
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


int arclm101(struct arclmframe *af,int idinput)
/*STATIC INCREMENTAL ANALYSIS.*/
{
  DWORDLONG memory0,memory1,memory2;
  char dir[]=DIRECTORY;                            /*DATA DIRECTORY*/
  char s[80],string[400],fname[256];
  clock_t t0,t1,t2;

  FILE *fin,*fout,*ffig,*ferr,*fsrf;                       /*FILE 8 BYTES*/
  /*FILE *felem,*fdisp,*freact;*/

  int i,ii,jj;

  int nnode,nelem,nshell,nsect,nconstraint,nreact;
  long int loffset,msize;

  /*ARCLMFRAME*/
  struct osect *sects;
  struct onode *nodes;
  struct onode *ninit;
  struct owire *elems;
  struct oconf *confs;
  struct memoryelem *melem;

  /*GLOBAL MATRIX*/
  struct gcomponent ginit={0,0,0.0,NULL};
  struct gcomponent *gmtx,*g,*p;
  double determinant,sign;

  /*GLOBAL VECTOR*/
  double *ddisp,*dreact;
  double *gvct;
  double *gvctcmq;         /*****CMQ ZOBUN*****/

  /*ELEMENT*/
  struct owire elem;
  double **drccos,**tmatrix,**estiff,*estress;

  /*FOR INCREMENTAL*/
  long int time;
  int nlap,laps;
  double safety,dsafety;

  double func[2];
  double initarea;
  double *ncr; /*BCLNG CONDENSATION RESULT*/      //ujioka





  memory0=availablephysicalmemory("INITIAL:");   /*MEMORY AVAILABLE*/

  /*OPEN INPUT FILE*/
#if 1
	fin=fgetstofopen(dir,"r",idinput);              /*OPEN INPUT FILE*/
	//fin = fgetstofopenII(dir, "r", (wdraw.childs+1)->inpfile);
	if(fin==NULL)
	{
	  errormessage("ACCESS IMPOSSIBLE.");
	  return 0;
	}
	inputinitII(fin, &nnode, &nelem, &nshell, &nsect, &nconstraint); /*INPUT INITIAL.*/
#endif
#if 0
	nnode=af->nnode;
	nelem=af->nelem;
	nshell=af->nshell;
	nsect=af->nsect;
	nconstraint=af->nconstraint;
#endif
	msize=6*nnode;                           /*SIZE OF GLOBAL MATRIX.*/

  //sprintf(string,"NODES=%d ELEMS=%d SECTS=%d",nnode,nelem,nsect);
  //errormessage(string);
  //if(fout!=NULL) fprintf(fout,"%s\n",string);


  /*READ ANALYSIS DATA*/
  getincrement((wmenu.childs+2)->hwnd,&laps,&dsafety);

  /*OPEN OUTPUT FILE*/
  fout=fgetstofopen("\0","w",ID_OUTPUTFILE);          /*OUTPUT FILE*/
  ffig=fopen("hogan.fig","w");                    /*HYSTERISIS FILE*/
  ferr=fopen("hogan.txt","w");                 /*ERROR MESSAGE FILE*/
  fsrf=fopen("surface.txt","w");                 /*ERROR MESSAGE FILE*/

  t0=clock();                                        /*CLOCK BEGIN.*/

  if(idinput==ID_INPUTFILE)
  {
  #if 1
	/*MEMORY NOT ALLOCATED*/
	free(af->sects);
	free(af->nodes);
	free(af->ninit);
	free(af->elems);
	free(af->shells);
	free(af->confs);
	free(af->iform);
	free(af->ddisp);
	free(af->melem);
	free(af->mshell);
	free(af->constraints);

	sects = (struct osect*)malloc(nsect * sizeof(struct osect));
	nodes = (struct onode*)malloc(nnode * sizeof(struct onode));
	ninit = (struct onode*)malloc(nnode * sizeof(struct onode));
	elems = (struct owire*)malloc(nelem * sizeof(struct owire));
	//shells = (struct oshell*)malloc(nshell * sizeof(struct oshell));
	confs = (struct oconf*)malloc(msize * sizeof(struct oconf));
	//iform = (double*)malloc(msize * sizeof(double));
	ddisp = (double*)malloc(msize * sizeof(double));
	melem = (struct memoryelem*)malloc(nelem * sizeof(struct memoryelem));
	//mshell = (struct memoryshell*)malloc(nshell * sizeof(struct memoryshell));
	//constraintmain = (long int*)malloc(msize * sizeof(long int));

	af->sects = sects;
	af->nodes = nodes;
	af->ninit = ninit;
	af->elems = elems;
	//af->shells = shells;
	af->confs = confs;
	//af->iform = iform;
	af->ddisp = ddisp;
	af->melem = melem;
	//af->mshell = mshell;
	//af->constraintmain = constraintmain;

	inputtexttomemory(fin, af);
	fclose(fin);
  #endif
	/*
	nnode=af->nnode;
    nelem=af->nelem;
    nsect=af->nsect;
	nreact=af->nreact;
	*/
    initialform(nodes,ddisp,nnode);         /*ASSEMBLAGE FORMATION.*/
    initialelem(elems,melem,nelem);          /*ASSEMBLAGE ELEMENTS.*/

    dreact=(double *)malloc(nreact*sizeof(double));     /*REACTION.*/
    af->dreact=dreact;
    initialreact(fin,dreact,nreact);   /*ASSEMBLAGE LONG REACTIONS.*/
  }
  else
  {

  #if 1
	/*MEMORY ALREADY ALLOCATED*/
	sects = af->sects;
	nodes = af->nodes;
	ninit = af->ninit;
	elems = af->elems;
	//shells = af->shells;
	confs = af->confs;
	//iform = af->iform;
	ddisp = af->ddisp;
	melem = af->melem;
	//mshell = af->mshell;
	//constraintmain = af->constraintmain;

	dreact=af->dreact;
  #endif
	/*
	nnode=af->nnode;
	nelem=af->nelem;
    nsect=af->nsect;
    nreact=af->nreact;
	*/
	inputloadtomemory(fin,af); /*INPUT LOAD.*/

    if(nreact!=af->nreact)
	{
	  MessageBox(NULL,"Input Error.","Arclm101",MB_OK);
	  return 0;
    }
  }

  /***GLOBAL MATRIX***/
  gmtx=(struct gcomponent *)malloc(msize*sizeof(struct gcomponent)); /*DIAGONALS OF GLOBAL MATRIX.*/
  /***GLOBAL VECTOR***/
  gvct=(double *)malloc(msize*sizeof(double));
  gvctcmq=(double *)malloc(msize*sizeof(double));/*CMQ ZOBUN*/

  for(i=0;i<=msize-1;i++)
  {
	(gmtx+i)->down=NULL;            /*GLOBAL MATRIX INITIALIZATION.*/
	*(gvct+i)=0.0;                  /*GLOBAL VECTOR INITIALIZATION.*/
	*(gvctcmq+i)=0.0;               /*GLOBAL VECTOR INITIALIZATION.*/
  }


  GetAsyncKeyState(VK_LBUTTON);                   /*CLEAR KEY LEFT.*/
  GetAsyncKeyState(VK_RBUTTON);                  /*CLEAR KEY RIGHT.*/
  if(globaldrawflag==1)
  {
	drawglobalaxis((wdraw.childs+1)->hdcC,(wdraw.childs+1)->vparam,0,0,255);  /*DRAW GLOBAL AXIS.*/
  }





  if(wsurf.hwnd!=NULL)
  {
    drawyieldsurface((wsurf.childs+1)->hdcC,
                     (wsurf.childs+1)->vparam,2,4,0,NULL);
    overlayhdc(*(wsurf.childs+1),SRCPAINT);       /*UPDATE DISPLAY.*/
  }
  af->fsurface=fopen("canbin.sfc","wb+");      /*STRESS ON SURFACE.*/


  /*****CMQ ZOBUN*****/
  for(i=1;i<=nelem;i++)               /*ASSEMBLAGE GLOBAL MATRIX.*/
  {
    elem=*(elems+i-1);                       /*READ ELEMENT DATA.*/

	inputnode(ddisp,elem.node[0]);                         /*HEAD*/
    inputnode(ddisp,elem.node[1]);                         /*TAIL*/

	drccos=directioncosine(elem.node[0]->d[0],
                           elem.node[0]->d[1],
                           elem.node[0]->d[2],
                           elem.node[1]->d[0],
						   elem.node[1]->d[1],
                           elem.node[1]->d[2],
                           elem.cangle);               /*[DRCCOS]*/

	tmatrix=transmatrix(drccos);         /*TRANSFORMATION MATRIX.*/

	modifycmq(melem,&elem);
	assemcmq101(elem,tmatrix,confs,gvctcmq,dsafety);    /*ASSEMBLAGE CMQ AS LOADS.*/

    for(ii=0;ii<=2;ii++) free(*(drccos+ii));
    free(drccos);
    for(ii=0;ii<=11;ii++) free(*(tmatrix+ii));
    free(tmatrix);
  }
  /*****CMQ ZOBUN*****/

  //definencr(&arc,&*ncr);          //ujioka

  for(nlap=1;nlap<=laps;nlap++)
  {
	af->nlaps=1/*nlap*/;

    sprintf(string,"LAP:%d/%d",nlap,laps);
	errormessage(string);
	if(fout!=NULL) fprintf(fout,"%s\n",string);
    if(ffig!=NULL) fprintf(ffig,"%s",string);

	setincrement((wmenu.childs+2)->hwnd,laps,nlap,dsafety,(nlap*dsafety));

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
	  *(gmtx+(i-1))=ginit;
	}

    comps=msize; /*INITIAL COMPONENTS=DIAGONALS.*/

    laptime("ASSEMBLING GLOBAL MATRIX.",t0);

    for(i=1;i<=nelem;i++)               /*ASSEMBLAGE GLOBAL MATRIX.*/
    {
      inputelem(elems,melem,i-1,&elem);        /*READ ELEMENT DATA.*/
	  for(ii=0;ii<=1;ii++)
      {
        for(jj=0;jj<=5;jj++)
        {
          (elems+i-1)->iconf[ii][jj]=elem.iconf[ii][jj];
		}
      }
      inputnode(ddisp,elem.node[0]);                         /*HEAD*/
      inputnode(ddisp,elem.node[1]);                         /*TAIL*/

	  elem.sect=(elems+i-1)->sect;             /*READ SECTION DATA.*/

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

	//clearwindow(*(wdraw.childs+1));
	drawarclmframe((wdraw.childs+1)->hdcC,(wdraw.childs+1)->vparam,*af,0,ONSCREEN);
	overlayhdc(*(wdraw.childs + 1), SRCPAINT);                  /*UPDATE DISPLAY.*/

	/*****CMQ ZOBUN*****/
    for(i=0;i<=msize-1;i++)
	{
      *(gvct+i)=*(gvctcmq+i);                      /*GLOBAL VECTOR.*/
	}
	assemconf(confs,gvct,dsafety,nnode);           /*GLOBAL VECTOR.*/

	modifygivend(gmtx,gvct,confs,nnode);   /*0:LOAD 1:DISPLACEMENT.*/

    laptime("CROUT LU DECOMPOSITION.",t0);
    croutludecomposition(gmtx,
                         gvct,confs,
                         6*nnode,
                         &determinant,&sign);        /*[K]{dU}={dF}*/

	//sprintf(string,"LOG{DETERMINANT}=%.5E SIGN=%.1f COMPS=%ld",determinant,sign,comps);
	//errormessage(string);
	//if(fout!=NULL) fprintf(fout,"%s\n",string);

	//safety=nlap*dsafety;
	//sprintf(string,"SAFETY FACTOR=%.5f",safety);
	//errormessage(string);
	//if(fout!=NULL) fprintf(fout,"%s\n",string);
	//if(ffig!=NULL) fprintf(ffig," SAFETY= %.5f",safety);

	if(sign<=0.0)
    {
      errormessage(" ");
      errormessage("INSTABLE TERMINATION.");
      if(fout!=NULL) fprintf(fout,"INSTABLE TERMINATION.\n");

      laptime("\0",t0);

      fclose(fin);
      fclose(fout);
	  fclose(ffig);

      gfree(gmtx,nnode); /*FREE GLOBAL MATRIX.*/
	  free(gvct);

      memory2=availablephysicalmemory("REMAIN:");
      sprintf(string,"CONSUMPTION:%ld[BYTES]",(memory0-memory2));
      errormessage(string);

      return 1;
    }

	//laptime("OUTPUT INTO FILE.",t0);

	//if(ferr!=NULL) fprintf(ferr,"LAP %3d / %3d",nlap,laps);

	//if(fout!=NULL) fprintf(fout,"\"DISPLACEMENT\"\n");
	//outputdisp(gvct,fout,nnode,nodes);  /*INCREMENTAL DISPLACEMENT.*/

	//if(fout!=NULL) fprintf(fout,"\"STRESS\"\n");
	for(i=1;i<=nelem;i++)                   /*STRESS OUTPUT,UPDATE.*/
	{
	  inputelem(elems,melem,i-1,&elem);

      inputnode(ddisp,elem.node[0]);
      inputnode(ddisp,elem.node[1]);

      elem.sect=(elems+i-1)->sect;             /*READ SECTION DATA.*/

	  estress=elemstress(&elem,gvct,melem,fout,func);

	  outputstress(elem,estress,fout,func);
      free(estress);

	}
	if(wsurf.hwnd!=NULL)
	{
	  drawyieldsurface((wsurf.childs+1)->hdcC,
                       (wsurf.childs+1)->vparam,
					   SURFACEX,SURFACEY,SURFACEZ,
					   af->fsurface);
	  overlayhdc(*(wsurf.childs+1),SRCPAINT);     /*UPDATE DISPLAY.*/
	}

	if(fout!=NULL) fprintf(fout,"\"REACTION\"\n");
	outputreaction(gmtx,gvct,nodes,confs,dreact,fout,nnode);

	updateform(ddisp,gvct,nnode);               /*FORMATION UPDATE.*/
	//if(fout!=NULL) fprintf(fout,"\"CURRENT FORM\"\n");

	for(ii=0;ii<nnode;ii++)
	{
	  sprintf(string,"NODE:%5ld {U}=",(nodes+ii)->code);
	  for(jj=0;jj<6;jj++)
	  {
		loffset=6*ii+jj;
		sprintf(s," %14.5f",*(ddisp+loffset));
		strcat(string,s);
	  }
	  if(fout!=NULL) fprintf(fout,"%s\n",string);
	}

	//if(ferr!=NULL) fprintf(ferr,"\n");

	t1=laptime("\0",t0);

	memory2=availablephysicalmemory(NULL);
	sprintf(string,"CONSUMPTION:%ld[BYTES]",(memory1-memory2));
	errormessage(string);


	//MESSAGE FOR UPDATE UI
	//AVOID FREEZE ON LONG RUNNING TASK
	MSG msg;
	while (PeekMessage(&msg, NULL, 0, 0, PM_REMOVE))
	{
		TranslateMessage(&msg);
		DispatchMessage(&msg);
	}

	errormessage("L:CONTINUE R:ABORT");            /*L=LEFT R=RIGHT*/
	while(!GetAsyncKeyState(VK_LBUTTON))  /*LEFT CLICK TO CONTINUE.*/
	{
	  if(GetAsyncKeyState(VK_RBUTTON))      /*RIGHT CLICK TO ABORT.*/
	  {
		if(fout!=NULL) fprintf(fout,"ABORTED.\n");
		errormessage(" ");
		errormessage("ABORTED.");

		fclose(fin);
		fclose(fout);
		fclose(ffig);

		gfree(gmtx,nnode); /*FREE GLOBAL MATRIX.*/
		free(gvct);

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
                              255,255,255,0,ONSCREEN);
	}
    overlayhdc(*(wdraw.childs+1),SRCPAINT);       /*UPDATE DISPLAY.*/
  }

  fclose(fin);
  fclose(fout);
  fclose(ffig);
  fclose(ferr);
  gfree(gmtx,nnode);

  errormessage(" ");
  errormessage("COMPLETED.");

  memory2=availablephysicalmemory("REMAIN:");
  sprintf(string,"CONSUMPTION:%ld[BYTES]",(memory0-memory2));
  errormessage(string);

  return 0;
}/*arclm101*/

int arclm101_bc(struct arclmframe *af,int idinput)
/*STATIC INCREMENTAL ANALYSIS.*/
/*UJIOKA :BUCKLING CONDENSATION APPLICATION*/
{
  DWORD memory0,memory1,memory2;

  FILE *fin,*fout,*ffig,*ferr,*fsrf/*,*ftest*/;                 /*FILE 8 BYTES*/
  double *ddisp,*dreact;
  struct memoryelem *melem;
  /*FILE *felem,*fdisp,*freact;*/
  char dir[]=DIRECTORY;                            /*DATA DIRECTORY*/
  char s[80],string[400]/*,fname[256]*/;
  int i,ii,jj;
  int nnode,nelem,nsect,nlap,laps,nreact;
  long int loffset,msize;
  long int time;
  struct gcomponent ginit={0,0,0.0,NULL};
  struct gcomponent *gmtx,*g,*p;                    /*GLOBAL MATRIX*/
  double *gvct;                                     /*GLOBAL VECTOR*/
  double *gvctcmq;         /*****CMQ ZOBUN*****/    /*GLOBAL VECTOR*/
  double **drccos,**tmatrix,**estiff,*estress;
  double determinant,sign,safety,dsafety;
  double func[2];
  clock_t t0,t1,t2;

  struct osect *sects;
  struct onode *nodes,*ninit;
  struct owire elem;
  struct owire *elems;
  struct oconf *confs;

  double initarea;

  double *ncr; /*BCLNG CONDENSATION RESULT*/      //ujioka

  /*****FOR VIERENDEEL ARCH BUCKLING EXPERIMENT*****/
  FILE *fdata=NULL;
  char **data;
  int ndata;
  double *thrust;
  /*****FOR VIERENDEEL ARCH BUCKLING EXPERIMENT*****/

  memory0=availablephysicalmemory("INITIAL:");   /*MEMORY AVAILABLE*/

  fin=fgetstofopen(dir,"r",idinput);              /*OPEN INPUT FILE*/
  if(fin==NULL)
  {
    errormessage("ACCESS IMPOSSIBLE.");
    return 0;
  }
  fout=NULL; ffig=NULL; fsrf=NULL;
  if(globalmessageflag)
  {
  fout=fgetstofopen("\0","w",ID_OUTPUTFILE);          /*OUTPUT FILE*/
  ffig=fopen("hogan.fig","w");                    /*HYSTERISIS FILE*/
  fsrf=fopen("surface.txt","w");                 /*FOR SURFACE FILE*/  //UJIOKA
  /*ftest=fopen("energytest.txt","w");*/       /*FOR STRAIN ENERGY TEST*/  //UJIOKA
  }
  ferr=fopen("hogan.txt","w");                 /*ERROR MESSAGE FILE*/
  /*****FOR VIERENDEEL ARCH BUCKLING EXPERIMENT*****/
//  if(globalmessageflag) fdata=fopen("thrust.txt","r");
  if(globalmessageflag) fdata=fopen("thrust2.txt","r");

  /*fgets(string,256,fin);*/                    /*INPUT APPELATION.*/
  /*errormessage(string);*/

  t0=clock();                                        /*CLOCK BEGIN.*/

  getincrement((wmenu.childs+2)->hwnd,&laps,&dsafety);

  inputinit(fin,&nnode,&nelem,&nsect);             /*INPUT INITIAL.*/
  sprintf(string,"NODES=%d ELEMS=%d SECTS=%d",nnode,nelem,nsect);
  errormessage(string);
  if(fout!=NULL) fprintf(fout,"%s\n",string);
  if(fsrf!=NULL) fprintf(fsrf,"%s\n",string);

  msize=6*nnode;                           /*SIZE OF GLOBAL MATRIX.*/

  gmtx=(struct gcomponent *)          /*DIAGONALS OF GLOBAL MATRIX.*/
        malloc(msize*sizeof(struct gcomponent));
  gvct=(double *)malloc(msize*sizeof(double));      /*GLOBAL VECTOR*/

  /*****CMQ ZOBUN*****/
  gvctcmq=(double *)malloc(msize*sizeof(double));   /*GLOBAL VECTOR*/

  ncr=mallocdoublevector(nelem);

  if(gmtx==NULL || gvct==NULL) return 0;
  for(i=0;i<=msize-1;i++)
  {
    (gmtx+i)->down=NULL;            /*GLOBAL MATRIX INITIALIZATION.*/
    *(gvct+i)=0.0;                  /*GLOBAL VECTOR INITIALIZATION.*/
    *(gvctcmq+i)=0.0;               /*GLOBAL VECTOR INITIALIZATION.*/
  }

  if(idinput==ID_INPUTFILE)
  {
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
    melem=(struct memoryelem *)
          malloc(nelem*sizeof(struct memoryelem));
    if(melem==NULL) return 0;

    af->sects=sects;
    af->nodes=nodes;
    af->ninit=ninit;
    af->elems=elems;
    af->confs=confs;
    af->ddisp=ddisp;                   /*DISPLACEMENT:6 DIRECTIONS.*/
    af->melem=melem;

    inputtexttomemory(fin,af);      /*READY TO READ LONG REACTIONS.*/
    nnode=af->nnode;
    nelem=af->nelem;
    nsect=af->nsect;
    nreact=af->nreact;

    initialform(nodes,ddisp,nnode);         /*ASSEMBLAGE FORMATION.*/
    initialelem(elems,melem,nelem);          /*ASSEMBLAGE ELEMENTS.*/

    dreact=(double *)malloc(nreact*sizeof(double));     /*REACTION.*/
    af->dreact=dreact;
    initialreact(fin,dreact,nreact);   /*ASSEMBLAGE LONG REACTIONS.*/
  }
  else
  {
    ddisp=af->ddisp;
    melem=af->melem;
    dreact=af->dreact;

    sects=af->sects;
    nodes=af->nodes;
    ninit=af->ninit;
    elems=af->elems;
    confs=af->confs;

    nnode=af->nnode;
    nelem=af->nelem;
    nsect=af->nsect;
    nreact=af->nreact;

    inputloadtomemory(fin,af); /*INPUT LOAD.*/

    if(nreact!=af->nreact)
    {
      MessageBox(NULL,"Input Error.","Arclm101",MB_OK);
      return 0;
    }
  }

  GetAsyncKeyState(VK_LBUTTON);                   /*CLEAR KEY LEFT.*/
  GetAsyncKeyState(VK_RBUTTON);                  /*CLEAR KEY RIGHT.*/
  if(globaldrawflag==1)
  {
  drawglobalaxis((wdraw.childs+1)->hdcC,(wdraw.childs+1)->vparam,
                 0,0,255);                      /*DRAW GLOBAL AXIS.*/
  }

  if(wsurf.hwnd!=NULL)
  {
    drawyieldsurface((wsurf.childs+1)->hdcC,
                     (wsurf.childs+1)->vparam,2,4,0,NULL);
    overlayhdc(*(wsurf.childs+1),SRCPAINT);       /*UPDATE DISPLAY.*/
  }
  af->fsurface=fopen("canbin.sfc","wb+");      /*STRESS ON SURFACE.*/

  /* ELEMENT ENERGY WITHOUT GRAVITY. */
  for (i = 0; i < nelem; i++)
  {
      //elem=*(elems+i);                       /*READ ELEMENT DATA.*/
      inputelem(elems,melem,i,&elem);        /*READ ELEMENT DATA.*/

      elem.Ee[0] = 0.0;
      elem.Ee[1] = 0.0;
      elem.Ep[0] = 0.0;
      elem.Ep[1] = 0.0;
  }

  /*fclose(fout);*/
#if 0
  /*****CMQ ZOBUN*****/
  for(i=1;i<=nelem;i++)               /*ASSEMBLAGE GLOBAL MATRIX.*/
  {
    elem=*(elems+i-1);                       /*READ ELEMENT DATA.*/

    inputnode(ddisp,elem.node[0]);                         /*HEAD*/
    inputnode(ddisp,elem.node[1]);                         /*TAIL*/

    drccos=directioncosine(elem.node[0]->d[0],
                           elem.node[0]->d[1],
                           elem.node[0]->d[2],
                           elem.node[1]->d[0],
                           elem.node[1]->d[1],
                           elem.node[1]->d[2],
                           elem.cangle);               /*[DRCCOS]*/

    tmatrix=transmatrix(drccos);         /*TRANSFORMATION MATRIX.*/

  	modifycmq(melem,&elem);
  	assemcmq101(elem,tmatrix,confs,gvctcmq,dsafety);    /*ASSEMBLAGE CMQ AS LOADS.*/

    for(ii=0;ii<=2;ii++) free(*(drccos+ii));
    free(drccos);
    for(ii=0;ii<=11;ii++) free(*(tmatrix+ii));
    free(tmatrix);
  }
  /*****CMQ ZOBUN*****/
#endif

  if(globalmessageflag) definencr(af,ncr);          //ujioka
  else
  {
     for(i=0;i<nelem;i++)
     {
       inputelem(elems,melem,i,&elem);        /*READ ELEMENT DATA.*/
/*
       sprintf(string,"ELEM %d :Nz=%5.8f,safety=%5.8f\n",(af->elems+i)->code,
              elem.stress[0][0],(elems+i)->srate[0]);
       errormessage(string);
*/
       if((elem.stress[0][0])*((elems+i)->srate[0])>1.0e-16)
       {
          ncr[i]=(elem.stress[0][0])/((elems+i)->srate[0]);
       }
       else  ncr[i]=0.0;
     }

  }
  if(ferr!=NULL) fprintf(ferr,"Buckling Loads(Ncr')\n");
  for(i=0;i<nelem;i++)
  {
      if(ferr!=NULL) fprintf(ferr,"ELEM %d :Ncr=%5.8f\n",(af->elems+i)->code,ncr[i]);
  }
  if(ferr!=NULL) fprintf(ferr,"\n");

  /*****FOR VIERENDEEL ARCH BUCKLING EXPERIMENT*****/
  if(fdata!=NULL)
  {
  thrust=(double *)malloc(laps*sizeof(double));     /*GLOBAL VECTOR*/
  for(i=0;i<laps;i++)
  {
	data=fgetsbrk(fdata,&ndata);
	if(!ndata) break;

    *(thrust+i)=strtod(*(data+0),NULL);
    *(thrust+i)/=1000;                             /*UNIT:[mm]->[m]*/

    for(;ndata>0;ndata--) free(*(data+ndata-1));
    free(data);
  }
  }
  /*****FOR VIERENDEEL ARCH BUCKLING EXPERIMENT*****/
  if(fdata!=NULL) fclose(fdata);
/*  for(nlap=1;nlap<=laps;nlap++) */
  nlap=0;
  while(1)
  {
    /*modified by UJIOKA*/
    nlap++;
    if(globalmessageflag)
    {
//      if(nlap>laps) break;
      if(nlap>2*laps) break;
    }
    else
    {
      if(nlap>2*laps) break;      /*FOR OPTIMIZATION*/
    }

    /*sprintf(fname,"arclm%d.lap",nlap);*/
    /*fout=fopen(fname,"w");*/                          /*LAP FILE*/


    /*af->nlaps=nlap;*/
    af->nlaps=1;

    sprintf(string,"LAP:%d/%d",nlap,laps);
    errormessage(string);
    if(globalmessageflag && fout!=NULL) fprintf(fout,"%s\n",string);
    if(globalmessageflag && ffig!=NULL) fprintf(ffig,"%s",string);
    if(globalmessageflag && fsrf!=NULL) fprintf(fsrf,"%s\n",string);

    if(globalmessageflag)
    {
    setincrement((wmenu.childs+2)->hwnd,
                 laps,nlap,dsafety,(nlap*dsafety));
    }
    else currentpivot(nlap,laps);

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
      inputelem(elems,melem,i-1,&elem);        /*READ ELEMENT DATA.*/
      for(ii=0;ii<=1;ii++)
      {
        for(jj=0;jj<=5;jj++)
        {
          (elems+i-1)->iconf[ii][jj]=elem.iconf[ii][jj];
        }
      }
      inputnode(ddisp,elem.node[0]);                         /*HEAD*/
      inputnode(ddisp,elem.node[1]);                         /*TAIL*/

      elem.sect=(elems+i-1)->sect;             /*READ SECTION DATA.*/

/*FOR STAINEDGLASS TEST CASE.*/
/*
if(nlap==1 && i==1) MessageBox(NULL,"Clearance for Staindglass","Arclm101",MB_OK);
if(nlap*dsafety<0.07)
{
  if(elem.sect->code==611 ||
     elem.sect->code==612 ||
     elem.sect->code==613 ||
     elem.sect->code==614 ||
     elem.sect->code==615)
  {
    initarea = elem.sect->area;
    elem.sect->area=0.0;
  }
}*/

      if((wdraw.childs+1)->hdcC!=NULL && globaldrawflag)     /*DRAW DEFORMED ELEMENT.*/
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
      if(elem.stress[0][0]>0)  //compression
      {
//   	    sprintf(string,"ELEM %d :Ncr=%5.8f\n",
//      			     	(af->elems+i-1)->code,ncr[i-1]);
//      errormessage(string);/*for check*/
//        estress=elemstressbc(&elem,gvct,melem,fout,fsrf,func,ncr[i-1]);
        estiff=assempmtxbc(elem,estiff,ncr[i-1]);          /*ADD PLASTIC MATRIX.*/
      }
      else estiff=assempmtxbc(elem,estiff,0);

      estiff=transformation(estiff,tmatrix);       /*[K]=[Tt][k][T]*/

      assemgstiffness(gmtx,estiff,&elem);             /*ASSEMBLAGE.*/

      for(ii=0;ii<=2;ii++) free(*(drccos+ii));
      free(drccos);
      for(ii=0;ii<=11;ii++) free(*(tmatrix+ii));
      free(tmatrix);
      for(ii=0;ii<=11;ii++) free(*(estiff+ii));
      free(estiff);

/*FOR STAINEDGLASS TEST CASE.*/
/*if(nlap*dsafety<0.1)
{
  if(elem.sect->code==611 ||
     elem.sect->code==612 ||
     elem.sect->code==613 ||
     elem.sect->code==614 ||
     elem.sect->code==615)
  {
    elem.sect->area=initarea;
  }
}*/
    }
    sprintf(string,"GLOBAL MATRIX %ld COMPS ASSEMBLED.",comps);
    laptime(string,t0);
    if(globaldrawflag)
    {
    overlayhdc(*(wdraw.childs+1),SRCPAINT);       /*UPDATE DISPLAY.*/
    }
    /*errormessage("ASSEMBLAGE GLOBAL VECTOR.");*/

    /*modified by ujioka for assemconf201*/
    for(i=0;i<=msize-1;i++)
    {
      *(gvct+i)=0.0;                          /*GLOBAL VECTOR.*/
    }

#if 0
    /*****CMQ ZOBUN*****/
    for(i=0;i<=msize-1;i++)
    {
      *(gvct+i)=*(gvctcmq+i);                      /*GLOBAL VECTOR.*/
    }
  	/*****CMQ ZOBUN*****/
#endif

#if 0
    /*****FOR VIERENDEEL ARCH BUCKLING EXPERIMENT*****/
    if(nlap==1)
    {
      *(gvct+(6*(1-1)+0))=-*(thrust+nlap)/2;          /*GLOBAL VECTOR. NODE 101 X*/
      *(gvct+(6*(10-1)+0))=*(thrust+nlap)/2;          /*GLOBAL VECTOR. NODE 110 X*/
    }
    if(nlap>1)
    {
      *(gvct+(6*(1-1)+0))=-(*(thrust+nlap)-*(thrust+nlap-1))/2;          /*GLOBAL VECTOR. NODE 101 X*/
      *(gvct+(6*(10-1)+0))=(*(thrust+nlap)-*(thrust+nlap-1))/2;          /*GLOBAL VECTOR. NODE 110 X*/
    }
/*
    sprintf(string,"thrust[%d]=%.12f",
            nlap,*(thrust+nlap));
    if(MessageBox(NULL,string,"ARCLM101",
			 MB_OKCANCEL)==IDCANCEL) break;
*/
//    if(nlap<=25)
//    {
//      *(gvct+(6*(1-1)+0))=-0.0000100;          /*GLOBAL VECTOR. NODE 101 X*/
//      *(gvct+(6*(10-1)+0))=0.0000100;          /*GLOBAL VECTOR. NODE 110 X*/
//    }
//    if(nlap>25)
//    {
//      *(gvct+(6*(1-1)+0))=-0.0000325;          /*GLOBAL VECTOR. NODE 101 X*/
//      *(gvct+(6*(10-1)+0))=0.0000325;          /*GLOBAL VECTOR. NODE 110 X*/
//    }
    /*****FOR VIERENDEEL ARCH BUCKLING EXPERIMENT*****/
#endif

//    assemconf(confs,gvct,dsafety,nnode);           /*GLOBAL VECTOR.*/
    if(nlap<=laps) assemconf201(confs,gvct,dsafety,nnode);        /*GLOBAL VECTOR.*/
    else           assemconf201(confs,gvct,-1.0*dsafety,nnode);
    modifygivend(gmtx,gvct,confs,nnode);   /*0:LOAD 1:DISPLACEMENT.*/

    laptime("CROUT LU DECOMPOSITION.",t0);
    croutludecomposition(gmtx,
                         gvct,confs,
                         6*nnode,
                         &determinant,&sign);        /*[K]{dU}={dF}*/

    sprintf(string,"LOG{DETERMINANT}=%.5E SIGN=%.1f COMPS=%ld",
            determinant,sign,comps);
    errormessage(string);
    if(globalmessageflag && fout!=NULL) fprintf(fout,"%s\n",string);

if(nlap<=laps)    safety=nlap*dsafety;
else              safety-=dsafety;
    sprintf(string,"SAFETY FACTOR=%.5f",safety);
    errormessage(string);
    if(globalmessageflag && fout!=NULL) fprintf(fout,"%s\n",string);
    if(globalmessageflag && ffig!=NULL) fprintf(ffig," SAFETY= %.5f",safety);

    if(sign<=0.0)
    {
      errormessage(" ");
      errormessage("INSTABLE TERMINATION.");
      if(globalmessageflag && fout!=NULL) fprintf(fout,"INSTABLE TERMINATION.\n");

      laptime("\0",t0);

      fclose(fin);
      fclose(fout);
      fclose(ffig);
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

    if(globalmessageflag && ferr!=NULL) fprintf(ferr,"LAP %3d / %3d",nlap,laps);

    if(globalmessageflag && fout!=NULL) fprintf(fout,"\"DISPLACEMENT\"\n");
    outputdisp(gvct,fout,nnode,nodes);  /*INCREMENTAL DISPLACEMENT.*/
    /*while(!GetAsyncKeyState(VK_LBUTTON))
    ;*/                                   /*LEFT CLICK TO CONTINUE.*/

    if(globalmessageflag && fout!=NULL) fprintf(fout,"\"STRESS\"\n");
    for(i=1;i<=nelem;i++)                   /*STRESS OUTPUT,UPDATE.*/
    {
      inputelem(elems,melem,i-1,&elem);

      inputnode(ddisp,elem.node[0]);
      inputnode(ddisp,elem.node[1]);

      elem.sect=(elems+i-1)->sect;             /*READ SECTION DATA.*/

/*FOR STAINEDGLASS TEST CASE.*/
/*if(nlap*dsafety<0.1)
{
  if(elem.sect->code==611 ||
     elem.sect->code==612 ||
     elem.sect->code==613 ||
     elem.sect->code==614 ||
     elem.sect->code==615)
  {
    initarea = elem.sect->area;
    elem.sect->area=0.0;
  }
}*/

	/***UJIOKA***/
      if(elem.stress[0][0]>0)  //compression
      {
  	  sprintf(string,"ELEM %d :Ncr=%5.8f\n",
					(af->elems+i-1)->code,ncr[i-1]);
//      errormessage(string);/*for check*/
      estress=elemstressbc(&elem,gvct,melem,fout,fsrf,func,ncr[i-1]);
      }
      else                     //tension
      {
	  estress=elemstressbc(&elem,gvct,melem,fout,fsrf,func,0);
      /*
        便宜上Ncr=0とし、この場合降伏曲面は修正しない。
        (updatestressbc,coefficientsbcでの条件分岐による)
      */
      }
	/***UJIOKA***/

    /***UJIOKA FOR STRAIN ENERGY***/
      (elems + i-1)->Ee[0] = elem.Ee[0];
      (elems + i-1)->Ee[1] = elem.Ee[1];
      (elems + i-1)->Ep[0] = elem.Ep[0];
      (elems + i-1)->Ep[1] = elem.Ep[1];
/*
      if(elem.code==1001)
      {
           fprintf(ftest,"ELEM %5ld %-6.12f %-6.12f %-6.12f %-6.12f\n",
                   elem.code,elem.Ee[0],elem.Ep[0],elem.Ee[1],elem.Ep[1]);
      }
*/
    /***UJIOKA FOR STRAIN ENERGY***/


/*FOR STAINEDGLASS TEST CASE.*/
/*if(elem.code==1004 ||
   elem.code==1007 ||
   elem.code==1011 ||
   elem.code==1013 ||
   elem.code==1016)
{
  fprintf(ferr," ELEM %5ld Nz %12.8f Mx %12.8f",
		  elem.code,elem.stress[0][0],elem.stress[0][4]);
}*/

      outputstress(elem,estress,fout,func);
      free(estress);

/*FOR STAINEDGLASS TEST CASE.*/
/*if(nlap*dsafety<0.1)
{
  if(elem.sect->code==611 ||
     elem.sect->code==612 ||
     elem.sect->code==613 ||
     elem.sect->code==614 ||
     elem.sect->code==615)
  {
    elem.sect->area=initarea;
  }
}*/
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

	if(globalmessageflag && fout!=NULL) fprintf(fout,"\"REACTION\"\n");
	outputreaction(gmtx,gvct,nodes,confs,dreact,fout,nnode);

	updateform(ddisp,gvct,nnode);               /*FORMATION UPDATE.*/
	if(globalmessageflag && fout!=NULL) fprintf(fout,"\"CURRENT FORM\"\n");
if(globalmessageflag)
{
	for(ii=0;ii<nnode;ii++)
	{
if(nlap<10)	          sprintf(string,"  NODE:%5ld {U}=",(nodes+ii)->code);
else if(nlap<100)	  sprintf(string," NODE:%5ld {U}=",(nodes+ii)->code);
else	              sprintf(string,"NODE:%5ld {U}=",(nodes+ii)->code);
	  for(jj=0;jj<6;jj++)
	  {
		loffset=6*ii+jj;
		sprintf(s," %14.5f",*(ddisp+loffset));
		strcat(string,s);
	  }
	  if(fout!=NULL) fprintf(fout,"%s\n",string);

      if(ffig!=NULL && (nodes+ii)->code==101)
	  /*if(ffig!=NULL && (nodes+ii)->code==164)*/ /*2FL*/
	  /*if(ffig!=NULL && (nodes+ii)->code==162)*/ /*3FL*/
	  /*if(ffig!=NULL && (nodes+ii)->code==194)*/ /*4FL*/
	  /*if(ffig!=NULL && (nodes+ii)->code==226)*/ /*5FL*/
	  /*if(ffig!=NULL && (nodes+ii)->code==258)*/ /*RFL*/
	  {
		fprintf(ffig," %s\n",string);
	  }
/*
	  if(ffig!=NULL && (nodes+ii)->code==107)
	  {
if(nlap<10)		        fprintf(ffig,"                          %s\n",string);
else if(nlap<100)		fprintf(ffig,"                           %s\n",string);
else	             	fprintf(ffig,"                            %s\n",string);
	  }
*/
/*FOR STAINEDGLASS TEST CASE.*/
/*if(ferr!=NULL &&
   (nodes+ii)->code==119 ||
   (nodes+ii)->code==128 ||
   (nodes+ii)->code==136 ||
   (nodes+ii)->code==145)
{
  fprintf(ferr," NODE %5ld Ux %12.8f",(nodes+ii)->code,*(ddisp+6*ii+0));
}*/
if(ferr!=NULL &&
   (nodes+ii)->code==101 ||
   (nodes+ii)->code==116/*||
   (nodes+ii)->code==122 ||
   (nodes+ii)->code==155 ||
   (nodes+ii)->code==192)
*/
  )
{
  fprintf(ferr," NODE %5ld Ux %12.8f",(nodes+ii)->code,*(ddisp+6*ii+0));
  /*fprintf(ferr," NODE %5ld Uy %12.8f",(nodes+ii)->code,*(ddisp+6*ii+1));*/
  /*fprintf(ferr," NODE %5ld Uz %12.8f",(nodes+ii)->code,*(ddisp+6*ii+2));*/
}
	}

if(ferr!=NULL) fprintf(ferr,"\n");
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

		gfree(gmtx,nnode); /*FREE GLOBAL MATRIX.*/
		free(gvct);
        /*free(confs);*/

        errormessage(" ");
        errormessage("ABORTED.");
        if(fout!=NULL) fprintf(fout,"ABORTED.\n");

        fclose(fout);
        fclose(ffig);

        laptime("\0",t0);
        return 1;
      }
      t2=clock();
      time=(t2-t1)/CLK_TCK;
      if(time>=WAIT) break;               /*CONTINUE AFTER WAITING.*/
    }

    /*fclose(fout);*/
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
      if(globaldrawflag)
      {
      drawglobalwire((wdraw.childs+1)->hdcC,
                     (wdraw.childs+1)->vparam,
                     *af,elem,255,255,255,
                              255,255,255,0,ONSCREEN);
	  }
    }
    if(globaldrawflag)
    {
    overlayhdc(*(wdraw.childs+1),SRCPAINT);       /*UPDATE DISPLAY.*/
    }
  }

  fclose(fin);
  /*fclose(felem);*/
  /*fclose(fdisp);*/
  /*fclose(freact);*/

  gfree(gmtx,nnode); /*FREE GLOBAL MATRIX.*/
  /*free(gvct);*/
  /*free(confs);*/
  free(ncr);         /*FREE Ncr'*/


  af->eigenvec=(double **)malloc(1*sizeof(double *));
  *((af->eigenvec)+0)=gvct;

  errormessage(" ");
  errormessage("COMPLETED.");
  if(fout!=NULL) fprintf(fout,"COMPLETED.\n");

  if(fout!=NULL) fclose(fout);
  if(ffig!=NULL) fclose(ffig);
  if(ferr!=NULL) fclose(ferr);
  if(fsrf!=NULL) fclose(fsrf);
  /*if(ftest!=NULL) fclose(ftest);*/

  memory2=availablephysicalmemory("REMAIN:");
  sprintf(string,"CONSUMPTION:%ld[BYTES]",(memory0-memory2));
  errormessage(string);
  errormessage(" ");

  return 0;
}/*arclm101_bc*/

int arclm201(struct arclmframe *af,int idinput)
{
  DWORD memory0,memory1,memory2;

  FILE *fin,*fout,*fonl,*ffig;                       /*FILE 8 BYTES*/
  double *ddisp,*dreact,*iform;
  struct memoryelem *melem;
  /*FILE *felem,*fdisp,*freact;*/
  char dir[]=DIRECTORY;                            /*DATA DIRECTORY*/
  char s[80],string[400]/*,fname[256]*/;
  int i,ii,jj;
  int nnode,nelem,nsect,nlap,laps,nreact;
  long int loffset,msize,fnode,nline;
  long int time;
  struct gcomponent ginit={0,0,0.0,NULL};
  struct gcomponent /**kmtx,*/*gemtx,*gmtx,*k,*ge,*g,*p;/*GLOBAL MATRIX*/
  double gg,kk;
  double *gvct,*gvct1,*gvct2;                       /*GLOBAL VECTOR*/
  double **drccos,**tmatrix,**estiff,*estress,*estress2,**tt;
  double determinant,sign,safety,dsafety;
  double func[2],tforce1,tforce2;

  /***UJIOKA FOR LOAD INCREMENTAL***/
  int iteration;
  double residual;

  int maxiteration=20;
  double tolerance=0.01;
  /***UJIOKA FOR LOAD INCREMENTAL***/


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

//  fnode=174;
  fnode=1961;
  ffig=fopen("hognon.fig","w");                       /*FIGURE FILE*/

  /*fgets(string,256,fin);*/                    /*INPUT APPELATION.*/
  /*errormessage(string);*/

  t0=clock();                                        /*CLOCK BEGIN.*/

  getincrement((wmenu.childs+2)->hwnd,&laps,&dsafety);

  inputinit(fin,&nnode,&nelem,&nsect);             /*INPUT INITIAL.*/
  sprintf(string,"NODES=%d ELEMS=%d SECTS=%d",nnode,nelem,nsect);
  errormessage(string);
  if(fout!=NULL) fprintf(fonl,"%s\n",string);

  msize=6*nnode;                           /*SIZE OF GLOBAL MATRIX.*/

  /*kmtx=(struct gcomponent *)*/         /*DIAGONALS OF ERASTIC MATRIX.*/
        /*malloc(msize*sizeof(struct gcomponent));*/
  gemtx=(struct gcomponent *)      /*DIAGONALS OF GEOMETRIC MATRIX.*/
		malloc(msize*sizeof(struct gcomponent));
  gmtx=(struct gcomponent *)          /*DIAGONALS OF GLOBAL MATRIX.*/
        malloc(msize*sizeof(struct gcomponent));
  gvct=(double *)malloc(msize*sizeof(double));           /*GLOBAL VECTOR.*/
  gvct1=(double *)malloc(msize*sizeof(double)); /*GLOBAL VECTOR ONLY CMQ.*/
  gvct2=(double *)malloc(msize*sizeof(double));     /*TRUE GLOBAL VECTOR.*/

  if(/*kmtx==NULL ||*/gemtx==NULL || gmtx==NULL || gvct==NULL || gvct1==NULL || gvct2==NULL) return 0;
  for(i=0;i<=msize-1;i++)
  {
    /*(kmtx+i)->down=NULL;*/           /*ERASTIC MATRIX INITIALIZATION.*/
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

  /*fclose(fout);*/
  nlap=1;
  iteration=1;
  residual=0.0;
//  for(nlap=1;nlap<=laps;nlap++)
  while(nlap<=laps)      /***UJIOKA FOR LOAD INCREMENTAL***/
  {
    /*sprintf(fname,"arclm%d.lap",nlap);*/
    /*fout=fopen(fname,"w");*/                          /*LAP FILE*/

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

	safety=nlap*dsafety;
	//if(safety>1.0) safety=1.0;
	//if(safety>0.1) safety=0.1;
	/***UJIOKA FOR LOAD INCREMENTAL***/

	setincrement((wmenu.childs+2)->hwnd,
                 laps,nlap,dsafety,safety);
    /*setincrement((wmenu.childs+2)->hwnd,
                 laps,nlap,dsafety,(nlap*dsafety));*/

    memory1=availablephysicalmemory("REMAIN:");  /*MEMORY AVAILABLE*/

	for(i=1;i<=msize;i++)           /*GLOBAL MATRIX INITIALIZATION.*/
	{
	  /*k=(kmtx+(i-1))->down;*/   /*NEXT OF DIAGONAL.*/
	  g=(gmtx+(i-1))->down;   /*NEXT OF DIAGONAL.*/

	  /*while(k!=NULL)*/ /*CLEAR ROW.*/
	  /*{
		p=k;
		k=k->down;
		free(p);
	  }*/

	  while(g!=NULL) /*CLEAR ROW.*/
      {
        p=g;
        g=g->down;
        free(p);
      }

      ginit.m=(unsigned short int)i;
      /*ginit.n=(unsigned short int)i;*/

      /**(kmtx+(i-1))=ginit;*/
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
							 elem.cangle);               /*[DRCCOS]*/

	  tmatrix=transmatrix(drccos);         /*TRANSFORMATION MATRIX.*/
	  estiff=assememtx(elem);          /*ELASTIC MATRIX OF ELEMENT.*/

#if 0
for(ii=0;ii<12;ii++)
{
  fprintf(fonl,"[ke](%2d)",ii);
  for(jj=0;jj<=ii;jj++)
  {
	fprintf(fonl," %14.5f",estiff[ii][jj]);
  }
  fprintf(fonl,"\n");
}
#endif

//      estiff=modifyhinge(elem,estiff);             /*MODIFY MATRIX.*/
      estiff=transformation(estiff,tmatrix);       /*[K]=[Tt][k][T]*/
      /*assemgstiffness(kmtx,estiff,&elem);*/     /*ASSEMBLAGE ELASTIC.*/
	  assemgstiffness(gmtx,estiff,&elem);     /*ASSEMBLAGE ELASTIC.*/

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


      for(ii=0;ii<=2;ii++) free(*(drccos+ii));
      free(drccos);
      for(ii=0;ii<=11;ii++) free(*(tmatrix+ii));
	  free(tmatrix);
      for(ii=0;ii<=11;ii++) free(*(estiff+ii));
      free(estiff);
    }

//modify by fukushima////////////////////////////////////////////////////////////

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


    #if 0
    for(ii=1;ii<=msize;ii++)
    {
      for(jj=1;jj<=ii;jj++)
      {
        gg=0.0;
        kk=0.0;

        gread(kmtx,ii,jj,&kk);
        gread(gemtx,ii,jj,&gg);

        kk+=gg;
        /*kk-=gg;*/

        /*if(ii==jj && kk<=0.0) kk=0.1;*/ /*TEST*/
        /*if(ii==jj) fprintf(fonl,"NODE:%ld / k(%ld,%ld)=%20.5f kg(%ld,%ld)=%20.5f\n"
                               ,(nodes+int((ii-1)/6))->code,ii,jj,kk,ii,jj,gg);*/

        gwrite(gmtx,ii,jj,kk);          /*ASSEMBLAGE GLOBAL MATRIX.*/
      }
    }
    #endif

    sprintf(string,"GLOBAL MATRIX %ld COMPS ASSEMBLED.",comps);
    laptime(string,t0);

    overlayhdc(*(wdraw.childs+1),SRCPAINT);       /*UPDATE DISPLAY.*/

	/*errormessage("ASSEMBLAGE GLOBAL VECTOR.");*/
	assemconf(confs,gvct,safety,nnode);           /*GLOBAL VECTOR.*/
	/*assemconf(confs,gvct,dsafety,nnode);*/           /*GLOBAL VECTOR.*/
///    assemconf201(confs,gvct,dsafety,nnode);        /*GLOBAL VECTOR.*/

//    for(i=0;i<msize;i++) *(gvct+i)+=*(gvct1+i)*safety; /*TARGET FORCE.?*/

    if(fonl!=NULL && /*nlap==1*/iteration==1) fprintf(fonl,"\"TARGET FORCE {F}\"\n");
	for(i=0;i<nnode;i++)
	{
	  if(fonl!=NULL && /*nlap==1*/iteration==1)
	  {
		fprintf(fonl,"NODE:%5ld %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f\n",(nodes+i)->code,
				*(gvct+(6*i+0)),*(gvct+(6*i+1)),*(gvct+(6*i+2)),
				*(gvct+(6*i+3)),*(gvct+(6*i+4)),*(gvct+(6*i+5)));
	  }

	  if((nodes+i)->code==fnode && ffig!=NULL)
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
        fprintf(fonl,"NODE:%5ld %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f\n",(nodes+i)->code,
                *(gvct2+(6*i+0)),*(gvct2+(6*i+1)),*(gvct2+(6*i+2)),
                *(gvct2+(6*i+3)),*(gvct2+(6*i+4)),*(gvct2+(6*i+5)));
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
        fprintf(fonl,"NODE:%5ld %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f\n",(nodes+i)->code,
                *(gvct2+(6*i+0)),*(gvct2+(6*i+1)),*(gvct2+(6*i+2)),
                *(gvct2+(6*i+3)),*(gvct2+(6*i+4)),*(gvct2+(6*i+5)));
      }
      fprintf(fonl,"\"UNBALANCED FORCE {dF}\"\n");
      for(i=0;i<nnode;i++)
      {
        fprintf(fonl,"NODE:%5ld %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f\n",(nodes+i)->code,
                *(gvct+(6*i+0)),*(gvct+(6*i+1)),*(gvct+(6*i+2)),
                *(gvct+(6*i+3)),*(gvct+(6*i+4)),*(gvct+(6*i+5)));
      }
    }



    laptime("CROUT LU DECOMPOSITION.",t0);
    nline=croutludecomposition(gmtx,
                               gvct,confs,
                               6*nnode,
                               &determinant,&sign);        /*[K]{dU}={dF}*/


	sprintf(string,"LOG{DETERMINANT}=%.5E SIGN=%.1f COMPS=%ld",determinant,sign,comps);
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

      laptime("\0",t0);

      fclose(fin);
      fclose(fonl);
      fclose(fout);
      fclose(ffig);

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
	  estress=elemstressnl(&elem,gvct,melem);                /*{f}.*/

	  estress2=(double *)malloc(12*sizeof(double));  /*ACCUMULATION STRESS*/
      for(ii=0;ii<=1;ii++)                                /*{f+df}.*/
      {
		for(jj=0;jj<6;jj++)
		{
		  *(estress2+(6*ii+jj))=0.0;
		  *(estress2+(6*ii+jj))=elem.stress[ii][jj];
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



//      estiff=assemgmtx(elem,estress);/*GEOMETRIC MATRIX OF ELEMENT.*/
//      estiff=modifyhinge(elem,estiff);             /*MODIFY MATRIX.*/


      estiff=transformation(estiff,tmatrix);       /*[K]=[Tt][k][T]*/
      assemgstiffness(gemtx,estiff,&elem);  /*ASSEMBLAGE GEOMETRIC.*/

      tt=matrixtranspose(tmatrix,12);                       /*[Tt].*/

//Modify by fukushima///////////////////////////////////////////////////////////
	  //estress2=matrixvectorIII(tt,estress2,12); /*TRUE FORCE OF ELEMENT {F}=[Tt]{f}.*/ /*COMMENT OUT BY KUTOMI*/
	  //estress2=matrixvector(tt,estress2,12); /*TRUE FORCE OF ELEMENT {F}=[Tt]{f}.*/
////////////////////////////////////////////////////////////////////////////////

//      estress2=matrixvector(tmatrix,estress2,12); /*TRUE FORCE OF ELEMENT {F}=[T]{f}.*/
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
        fprintf(fonl,"%s\n",string);

        if((nodes+ii)->code==fnode && ffig!=NULL)
        {
          fprintf(ffig,"LAP: %3d / %3d NODE %3d {Fz}= %18.12f {U}= %16.12f\n",
                  nlap,laps,(nodes+ii)->code,tforce1,*(ddisp+6*ii+2));
          fprintf(ffig,"LAP: %3d / %3d NODE %3d {Fz}= %18.12f {U}= %16.12f\n",
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
  fclose(ffig);


  gfree(gemtx,nnode); /*FREE GEOMETRIC MATRIX.*/
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
}/*arclm201*/




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

#include <windows.h>

DWORDLONG availablephysicalmemoryEx(char *comment)
/*RETURN: PHYSICAL MEMORY AVAILABLE [BYTES].*/
{
    MEMORYSTATUSEX memstat;
	memstat.dwLength = sizeof(memstat); // This is important for MEMORYSTATUSEX!
    char com[256], str[256];

    if (!GlobalMemoryStatusEx(&memstat)) {
        // Handle the error, maybe return 0 or another indicator
		return 0;
    }

    if (comment != NULL) {
        sprintf(str, "%llu/%llu[BYTES] AVAILABLE.", memstat.ullAvailPhys, memstat.ullTotalPhys);
		strcpy(com, comment);
		strcat(com, str);
		errormessage(com); // Assuming errormessage handles the message appropriately
	}

	return memstat.ullAvailPhys;
}/*availablephysicalmemoryEx*/

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



FILE *fgetstofopenII(const char *directory,const char *mode,const char *filename)
/*INPUT FILENAME AND THEN OPEN FILE.*/
{
  FILE *f=NULL;
  char fname[256],dandf[256];

	strcpy(fname,filename);
	strcpy(dandf,directory);
	strcat(dandf,fname);
	f=fopen(dandf,mode);

  return f;
}/*fgetstofopenII*/


void inputtexttomemory(FILE *ftext,struct arclmframe *af)
/*TRANSLATE ARCLM INPUTFILE TEXT INTO MEMORY.*/
{
  char **data;
  int i,j,ii,jj,k,l,m,n;
  long int offset,mainoff,suboff;
  long int scode,ncode,dcode,code1,code2,code3,code4;
  char str[50];/*for debug*/

  fseek(ftext,0L,SEEK_SET);
  inputinitII(ftext,&(af->nnode),&(af->nelem),&(af->nshell),&(af->nsect),&(af->nconstraint));

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

	if((af->sects+i)->area==0.0 &&
	   (af->sects+i)->Ixx ==0.0 &&
	   (af->sects+i)->Iyy ==0.0 &&
	   (af->sects+i)->Jzz ==0.0)
	{
	  (af->sects+i)->dflag=0;
	}
	else (af->sects+i)->dflag=1;

	/*
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
	*/
	(af->sects+i)->ppc.npcurve=0;
  }


  for(i=0;i<(af->nnode);i++)
  {
	(af->nodes+i)->loff=i;
	(af->ninit+i)->loff=i;

	data=fgetsbrk(ftext,&n);
	(af->nodes + i)->code = strtol(*(data + 0), NULL, 10);
	(af->ninit + i)->code = (af->nodes + i)->code;

	if (n == 13)
	{
		for (ii = 0; ii < 3; ii++)
		{
			(af->nodes + i)->d[ii] = strtod(*(data + 1 + ii), NULL);
			(af->nodes + i)->r[ii] = strtod(*(data + 4 + ii), NULL);
			(af->ninit + i)->d[ii] = strtod(*(data + 7 + ii), NULL);
			(af->ninit + i)->r[ii] = strtod(*(data + 10 + ii), NULL);
		}
	}
	if (n == 4)
	{
		for (ii = 0; ii < 3; ii++)
		{
			(af->nodes + i)->d[ii] = strtod(*(data + 1 + ii), NULL);
			(af->nodes + i)->r[ii] = 0.0;
			(af->ninit + i)->d[ii] = (af->nodes + i)->d[ii];
			(af->ninit + i)->r[ii] = (af->nodes + i)->r[ii];
		}
	}
	if (n == 7)
	{
		for (ii = 0; ii < 3; ii++)
		{
			(af->nodes + i)->d[ii] = strtod(*(data + 1 + ii), NULL);
			(af->nodes + i)->r[ii] = strtod(*(data + 4 + ii), NULL);
			(af->ninit + i)->d[ii] = (af->nodes + i)->d[ii];
			(af->ninit + i)->r[ii] = (af->nodes + i)->r[ii];
		}
	}

	for(;n>0;n--) free(*(data+n-1));
	free(data);
  }

  //initialform(af->ninit, af->iform, af->nnode);           /*ASSEMBLAGE FORMATION.*/
  //initialform(af->nodes, af->ddisp, af->nnode);           /*ASSEMBLAGE FORMATION.*/

  for(i=1;i<=(af->nelem);i++)
  {
	(af->elems+i-1)->loff=i-1;
	data=fgetsbrk(ftext,&n);
	(af->elems+i-1)->code=strtol(*(data+0),NULL,10);
	scode=strtol(*(data+1),NULL,10); /*SECTION.*/

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


	code1=strtol(*(data+2),NULL,10); /*HEAD NODE.*/
	code2=strtol(*(data+3),NULL,10); /*TAIL NODE.*/
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
	for(ii=0;ii<2;)                                        /*NODES.*/
	{
	  if((af->nodes+offset)->code==code1)
	  {
		(af->elems+i-1)->node[0]=af->nodes+offset;
		ii++;
	  }
	  if((af->nodes+offset)->code==code2)
	  {
		(af->elems+i-1)->node[1]=af->nodes+offset;
		ii++;
	  }
	  offset++;
	}
	for(;n>0;n--) free(*(data+n-1));
	free(data);
  }

  for(i=1;i<=(af->nshell);i++)
  {
	(af->shells+i-1)->loff=i-1;
	data=fgetsbrk(ftext,&n);
	(af->shells+i-1)->code=strtol(*(data+0),NULL,10);
	scode=strtol(*(data+1),NULL,10); /*SECTION.*/

	offset=0;
	for(ii=0;ii<1;)                                      /*SECTION.*/
	{
	  if((af->sects+offset)->code==scode)
	  {
		(af->shells+i-1)->sect=af->sects+offset;
		ii++;
	  }
	  offset++;
	}
	offset=0;
	if(((af->shells+i-1)->sect)->type==7)
	{
	  (af->shells+i-1)->nnod=3;
	  (af->shells+i-1)->ngp=4;
	  //(af->shells+i-1)->ngp=7;
	  (af->shells+i-1)->nstress=6;
	  code1=strtol(*(data+2),NULL,10); /*HEAD NODE.*/
	  code2=strtol(*(data+3),NULL,10); /*TAIL NODE.*/
	  code3=strtol(*(data+4),NULL,10); /*NODE.*/
	  k=5;
	  for(ii=0;ii<3;)                                        /*NODES.*/
	  {
		if((af->nodes+offset)->code==code1)
		{
		  (af->shells+i-1)->node[0]=af->nodes+offset;
		  ii++;
		}

		if((af->nodes+offset)->code==code2)
		{
		  (af->shells+i-1)->node[1]=af->nodes+offset;
		  ii++;
		}

		if((af->nodes+offset)->code==code3)
		{
		  (af->shells+i-1)->node[2]=af->nodes+offset;
		  ii++;
		}
		offset++;
	  }
	}
	else
	{
	  (af->shells+i-1)->nnod=4;
	  (af->shells+i-1)->ngp=9;
	  (af->shells+i-1)->nstress=6;
	  code1=strtol(*(data+2),NULL,10);
	  code2=strtol(*(data+3),NULL,10);
	  code3=strtol(*(data+4),NULL,10);
	  code4=strtol(*(data+5),NULL,10);
	  k=6;
	  for(ii=0;ii<4;)
	  {
		if((af->nodes+offset)->code==code1)
		{
		  (af->shells+i-1)->node[0]=af->nodes+offset;
		  ii++;
		}

		if((af->nodes+offset)->code==code2)
		{
		  (af->shells+i-1)->node[1]=af->nodes+offset;
		  ii++;
		}

		if((af->nodes+offset)->code==code3)
		{
		  (af->shells+i-1)->node[2]=af->nodes+offset;
		  ii++;
		}

		if((af->nodes+offset)->code==code3)
		{
		  (af->shells+i-1)->node[3]=af->nodes+offset;
		  ii++;
		}
		offset++;
	  }
	}
	for(ii=0;ii<(af->shells+i-1)->nnod;ii++)                                  /*STRESS.*/
	{
	  for(jj=0;jj<6;jj++)
	  {
		(af->shells+i-1)->stress[ii][jj]=strtod(*(data+k),NULL);
		k++;
	  }
	}
	for(ii=0;ii<(af->shells+i-1)->ngp;ii++)                                  /*STRESS.*/
	{
	  for(jj=0;jj<(af->shells+i-1)->nstress;jj++)
	  {
		((af->shells+i-1)->gp[ii]).estrain[jj]=0.0;
		((af->shells+i-1)->gp[ii]).pstrain[jj]=0.0;
		((af->shells+i-1)->gp[ii]).stress[jj]=0.0;
		((af->shells+i-1)->gp[ii]).backstress[jj]=0.0;
	  }
	  ((af->shells+i-1)->gp[ii]).qn=0.0;
	  ((af->shells+i-1)->gp[ii]).qm=0.0;
	  ((af->shells+i-1)->gp[ii]).qnm=0.0;
	  ((af->shells+i-1)->gp[ii]).yinit=0.0;
	  ((af->shells+i-1)->gp[ii]).y=0.0;
	  ((af->shells+i-1)->gp[ii]).f[0]=-1.0;
	  ((af->shells+i-1)->gp[ii]).f[1]=-1.0;
	  ((af->shells+i-1)->gp[ii]).lambda[0]=0.0;
	  ((af->shells+i-1)->gp[ii]).lambda[1]=0.0;
	  ((af->shells+i-1)->gp[ii]).alpha=0.0;
	}

	/*INITIAL AREA OF SHELL*/
	(af->shells+i-1)->area = shellarea(*(af->shells+i-1));
	/*WEIGHT OF GAUSS INTEGRATION POINT*/
	/*
	(af->shells+i-1)->w[0]=27.0/60.0;
	(af->shells+i-1)->w[1]=8.0/60.0;
	(af->shells+i-1)->w[2]=(af->shells+i-1)->w[1];
	(af->shells+i-1)->w[3]=(af->shells+i-1)->w[1];
	(af->shells+i-1)->w[4]=3.0/60.0;
	(af->shells+i-1)->w[5]=(af->shells+i-1)->w[4];
	(af->shells+i-1)->w[6]=(af->shells+i-1)->w[4];
	*/

	(af->shells+i-1)->w[0]=0.0;
	(af->shells+i-1)->w[1]=1.0/3.0;
	(af->shells+i-1)->w[2]=(af->shells+i-1)->w[1];
	(af->shells+i-1)->w[3]=(af->shells+i-1)->w[1];

	if(k==n-2)
	{
	  (af->shells+i-1)->prate=strtod(*(data+n-2),NULL);
	  (af->shells+i-1)->brate=strtod(*(data+n-1),NULL);
	  k+=2;
	}
	else
	{
	  (af->shells+i-1)->prate=1.0;
	  (af->shells+i-1)->brate=1.0;
	}

	for(;n>0;n--) free(*(data+n-1));
	free(data);
  }

  af->nreact=0;
  for(i=1;i<=(af->nnode);i++) /*CONF VECTOR:CONFINEMENT,VALUE.*/
  {
	data=fgetsbrk(ftext,&n);

	for(j=1;j<=6;j++)
	{
	  offset=6*(i-1)+(j-1);

	  (af->confs+offset)->iconf = (signed char)strtol(*(data+j),NULL,10);
	  (af->confs+offset)->value =              strtod(*(data+j+6),NULL);

	  if((af->confs+offset)->iconf==1) (af->nreact)++;

	  if(n==19 && af->nmass!=NULL)
	  {
		*(af->nmass+offset) = strtod(*(data+j+12),NULL);
	  }
	}


	for(;n>0;n--) free(*(data+n-1));
	free(data);
  }

  for (i = 0; i < 6*(af->nnode); i++)
  {
	*((af->constraintmain) + i) = i;
  }
  offset = 0;/*low offset of constraint equation matrix*/
  for(i=1;i<=(af->nconstraint);i++) /*CONF VECTOR:CONFINEMENT,VALUE.*/
  {
	(af->constraints+i-1)->loff=i-1;
	(af->constraints+i-1)->leq=offset;
	data=fgetsbrk(ftext,&n);
	(af->constraints+i-1)->code=strtol(*(data+0),NULL,10);
	(af->constraints+i-1)->type=strtol(*(data+1),NULL,10);

	if((af->constraints+i-1)->type==0)/*MPC.*/
	{
		(af->constraints+i-1)->nnod=2;
		(af->constraints+i-1)->neq=0;
		code1=strtol(*(data+2),NULL,10);
		code2=strtol(*(data+3),NULL,10);
		j=0;
		for(k=0;k<2;)
		{
			if((af->nodes+j)->code==code1)
			{
				(af->constraints+i-1)->node[0]=af->nodes+j;
				mainoff = 6*j;
				k++;
			}
			if((af->nodes+j)->code==code2)
			{
				(af->constraints+i-1)->node[1]=af->nodes+j;
				suboff = 6*j;
				for(l=0;l<6;l++)
				{
					*((af->constraintmain)+suboff+l)=mainoff+l;
				}
				k++;
			}
			j++;
		}
	}
	else if((af->constraints+i-1)->type==1)/*REVOLUTE JOINT.*/
	{
		(af->constraints+i-1)->nnod=2;
		(af->constraints+i-1)->neq=2;
		code1=strtol(*(data+2),NULL,10);
		code2=strtol(*(data+3),NULL,10);
		j=0;
		for(k=0;k<2;)
		{
			if((af->nodes+j)->code==code1)
			{
				(af->constraints+i-1)->node[0]=af->nodes+j;
				mainoff = 6*j;
				k++;
			}
			if((af->nodes+j)->code==code2)
			{
				(af->constraints+i-1)->node[1]=af->nodes+j;
				suboff = 6*j;
				for(l=0;l<3;l++)
				{
					*((af->constraintmain)+suboff+l)=mainoff+l;
				}
				k++;
			}
			j++;
		}
		for(j=0;j<3;j++)
		{
			(af->constraints+i-1)->axis[2][j]=strtod(*(data+4+j),NULL);
		}
		if((af->constraints+i-1)->axis[2][0] != 0.0 || (af->constraints+i-1)->axis[2][1] != 0.0)
		{
			(af->constraints+i-1)->axis[0][0] =  -(af->constraints+i-1)->axis[2][1];
			(af->constraints+i-1)->axis[0][1] =   (af->constraints+i-1)->axis[2][0];
			(af->constraints+i-1)->axis[0][2] =   0.0;
		}
		else
		{
			(af->constraints+i-1)->axis[0][0] =  1.0;
			(af->constraints+i-1)->axis[0][1] =  0.0;
			(af->constraints+i-1)->axis[0][2] =  0.0;
		}
		double* cross = crossproduct((af->constraints+i-1)->axis[2],(af->constraints+i-1)->axis[0]);
		for(j=0;j<3;j++)
		{
			(af->constraints+i-1)->axis[1][j]=cross[j];
		}
		free(cross);
		vectornormalize((af->constraints+i-1)->axis[0],3);
		vectornormalize((af->constraints+i-1)->axis[1],3);
		vectornormalize((af->constraints+i-1)->axis[2],3);
	}
	else if((af->constraints+i-1)->type==2)/*SPHERICAL JOINT.*/
	{
		(af->constraints+i-1)->nnod=2;
		(af->constraints+i-1)->neq=0;
		code1=strtol(*(data+2),NULL,10);
		code2=strtol(*(data+3),NULL,10);
        j=0;
		for(k=0;k<2;)
		{
			if((af->nodes+j)->code==code1)
			{
				(af->constraints+i-1)->node[0]=af->nodes+j;
				mainoff = 6*j;
				k++;
			}
			if((af->nodes+j)->code==code2)
			{
				(af->constraints+i-1)->node[1]=af->nodes+j;
				suboff = 6*j;
				for(l=0;l<3;l++)
				{
					*((af->constraintmain)+suboff+l)=mainoff+l;
				}
				k++;
			}
			j++;
		}
	}
	else if((af->constraints+i-1)->type==3)/*PRISMATIC JOINT.*/
	{
		(af->constraints+i-1)->nnod=2;
		(af->constraints+i-1)->neq=3;
		code1=strtol(*(data+2),NULL,10);
		code2=strtol(*(data+3),NULL,10);
	}
	else if((af->constraints+i-1)->type==4)/*CYLINDRICAL JOINT.*/
	{
		(af->constraints+i-1)->nnod=2;
		(af->constraints+i-1)->neq=3;
		code1=strtol(*(data+2),NULL,10);
		code2=strtol(*(data+3),NULL,10);
	}
	else if((af->constraints+i-1)->type==5)/*UNIVERSAL JOINT.*/
	{
		(af->constraints+i-1)->nnod=2;
		(af->constraints+i-1)->neq=3;
		code1=strtol(*(data+2),NULL,10);
		code2=strtol(*(data+3),NULL,10);
	}
	else/*MPC.*/
	{
		(af->constraints+i-1)->nnod=0;
		(af->constraints+i-1)->neq=0;
		for(j=0;j<n/(int)3;j++)
		{
		  ncode=strtol(*(data+3*j+2),NULL,10);
		  dcode=strtol(*(data+3*j+3),NULL,10);
		  if(j==0)
		  {
			for(k=0;k<af->nnode;k++)
			{
			  if((af->nodes+k)->code==ncode)
			  {
				mainoff=6*k+dcode;
				/* *(*((af->constraintvec)+i-1)+mainoff)=strtod(*(data+3*j+2),NULL);*/
				break;
			  }
			}
		  }
		  else
		  {
			for(k=0;k<af->nnode;k++)
			{
			  if((af->nodes+k)->code==ncode)
			  {
				suboff=6*k+dcode;
				/* *(*((af->constraintvec)+i-1)+suboff)=strtod(*(data+3*j+2),NULL);*/
				*((af->constraintmain)+suboff)=mainoff;
				break;
			  }
			}
		  }

		}
		/* *((af->constraintval)+i-1)=strtod(*(data+n-1),NULL);*/
	}

	offset += (af->constraints+i-1)->neq;

	for(;n>0;n--) free(*(data+n-1));
	free(data);

  }


  /*LONGREACTION:UNDER CONSTRUCTION.*/
  for (i = 0; i < 6*(af->nnode); i++)
  {
	if (*((af->constraintmain) + i) != i)
	{
	  (af->confs + i)->iconf = (signed char)1;
	}
  }

  return;
}/*inputtexttomemory*/

void inputloadtomemory(FILE *ftext,struct arclmframe *af)
/*READ ONLY LOAD FROM ARCLM INPUTFILE INTO MEMORY.*/
{
  char **data,str[256];
  int i,j,n;
  long int offset;

  fseek(ftext,0L,SEEK_SET);
  fgets(str,256,ftext); /*INITIAL.*/

  for(i=0;i<(af->nsect);i++)
  {
	fgets(str,256,ftext); /*SECTIONS.*/
  }
  for(i=0;i<(af->nnode);i++)
  {
	fgets(str,256,ftext); /*NODES.*/
  }
  for(i=0;i<(af->nelem);i++)
  {
	fgets(str,256,ftext); /*ELEMENTS.*/
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
}/*inputloadtomemory*/

void inputframetomemory(FILE *ftext,struct arclmframe *af)
/*TRANSLATE FRAME INPUTFILE INTO MEMORY.*/
{
  char **data,msg[256];
  int i,j,ii,jj,k,n;
  long int offset;
  long int scode,hcode,tcode;
  double E,G,efact,gfact;

  fseek(ftext,0L,SEEK_SET);

  data=fgetsbrk(ftext,&n);
  af->nnode=strtol(*(data+0),NULL,10);
  af->nelem=strtol(*(data+1),NULL,10);
  E=strtod(*(data+5),NULL);
  G=strtod(*(data+6),NULL);
  af->nsect=strtol(*(data+12),NULL,10);
  for(;n>0;n--) free(*(data+n-1));
  free(data);

  /*"INPUTFILE SPECIFICATION"*/
  /*NNODE NELEM MSIZE BAND 0 E G 1.0E-04 1.0E+8 0 1 0 NSECT*/
  /*ISECT TYPE A Ixx Iyy J*/
  /*INODE X Y Z*/
  /*IELEM NODEI J SECT BOUNDARY... COORD EFACT GFACT CMQTYPE CMQ...*/
  /*INODE CONFINEMENTVALUE... CONFINEMENTTYPE...*/

  for(i=0;i<(af->nsect);i++)
  {
    (af->sects+i)->loff=i;

    data=fgetsbrk(ftext,&n);

    if(n!=6)
    {
      sprintf(msg,"SECT IN OFFSET=%d DATA ERROR.",i);
      errormessage(msg);
    }

    (af->sects+i)->code=strtol(*(data+0),NULL,10);
    (af->sects+i)->area=strtod(*(data+2),NULL);
    (af->sects+i)->Ixx =strtod(*(data+3),NULL);
    (af->sects+i)->Iyy =strtod(*(data+4),NULL);
    (af->sects+i)->Jzz =strtod(*(data+5),NULL);

    for(;n>0;n--) free(*(data+n-1));
    free(data);

    (af->sects+i)->E=E;
    (af->sects+i)->poi=0.5*E/G-1.0;

    for(ii=0;ii<=5;ii++)      /*UPPER,LOWER LIMIT OF YIELD SURFACE.*/
    {
      (af->sects+i)->fmax[ii]=1.0;
      (af->sects+i)->fmin[ii]=1.0;
    }

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
    if(n!=4)
    {
      sprintf(msg,"NODE IN OFFSET=%d DATA ERROR.",i);
      errormessage(msg);
    }

    (af->nodes+i)->code=strtol(*(data+0),NULL,10);
    (af->nodes+i)->d[0]=strtod(*(data+1),NULL);
    (af->nodes+i)->d[1]=strtod(*(data+2),NULL);
    (af->nodes+i)->d[2]=strtod(*(data+3),NULL);

    *(af->ninit+i)=*(af->nodes+i);

    for(;n>0;n--) free(*(data+n-1));
    free(data);
  }
  for(i=0;i<(af->nelem);i++)
  {
    (af->elems+i)->loff=i;

    data=fgetsbrk(ftext,&n);
    if(n!=14 && n!=26)
    {
      sprintf(msg,"ELEM IN OFFSET=%d DATA ERROR.",i);
      errormessage(msg);
    }

    (af->elems+i)->code=strtol(*(data+0),NULL,10);
    hcode=strtol(*(data+1),NULL,10); /*HEAD NODE.*/
    tcode=strtol(*(data+2),NULL,10); /*TAIL NODE.*/
    scode=strtol(*(data+3),NULL,10); /*SECTION.*/
	(af->elems+i)->cangle=strtod(*(data+10),NULL);

    k=4;
    for(ii=0;ii<=1;ii++)                 /*BOUNDARY.0:RIGID 1:HINGE*/
    {
      for(jj=0;jj<=2;jj++)
      {
        (af->elems+i)->iconf[ii][jj]=0;
      }
      for(jj=3;jj<=5;jj++)
      {
        (af->elems+i)->iconf[ii][jj]
        =(signed char)strtol(*(data+k),NULL,10);
        k++;
      }
    }

    efact=strtod(*(data+11),NULL);
    gfact=strtod(*(data+12),NULL);
    if(efact!=1.0 || gfact!=1.0)
    {
      MessageBox(NULL,"FACTOR E,G NOT 1.0.","Input Frame",MB_OK);
    }

    k=14;
    for(ii=0;ii<=1;ii++)                                  /*STRESS.*/
    {
      for(jj=0;jj<=5;jj++)
      {
        if(n<=14) (af->elems+i)->stress[ii][jj]=0.0;
        else (af->elems+i)->stress[ii][jj]=strtod(*(data+k),NULL);
        k++;
      }
    }

    for(;n>0;n--) free(*(data+n-1));
    free(data);

    offset=0;
    for(ii=0;ii<1;)                                      /*SECTION.*/
    {
      if(offset>=(af->nsect))
      {
        sprintf(msg,"ELEM %ld:SECTION %ld NOT FOUND.",
                (af->elems+i)->code,scode);
        MessageBox(NULL,msg,"Input Frame",MB_OK);
        break;
      }

      if((af->sects+offset)->code==scode)
      {
        (af->elems+i)->sect=af->sects+offset;
        ii++;
      }
      offset++;
    }

    offset=0;
    for(ii=0;ii<2;)                                        /*NODES.*/
    {
      if(offset>=(af->nnode))
      {
        sprintf(msg,"ELEM %ld:NODE %ld OR %ld NOT FOUND.",
                (af->elems+i)->code,hcode,tcode);
        MessageBox(NULL,msg,"Input Frame",MB_OK);
        break;
      }

      if((af->nodes+offset)->code==hcode)
      {
        (af->elems+i)->node[0]=af->nodes+offset;
        ii++;
      }
      if((af->nodes+offset)->code==tcode)
      {
        (af->elems+i)->node[1]=af->nodes+offset;
        ii++;
      }
      offset++;
    }
  }

  af->nreact=0;
  for(i=0;i<(af->nnode);i++) /*CONF VECTOR:CONFINEMENT,VALUE.*/
  {
    data=fgetsbrk(ftext,&n);
    if(n!=13)
    {
      sprintf(msg,"CONF IN OFFSET=%d DATA ERROR.",i);
      errormessage(msg);
    }

    for(j=1;j<=6;j++)
    {
      offset=6*i+(j-1);

      (af->confs+offset)->iconf
      =(signed char)strtol(*(data+j+6),NULL,10);
      (af->confs+offset)->value
      =strtod(*(data+j),NULL);

      if((af->confs+offset)->iconf==1) (af->nreact)++;
    }

    for(;n>0;n--) free(*(data+n-1));
    free(data);
  }
  /*LONGREACTION:UNDER CONSTRUCTION.*/

  return;
}/*inputframetomemory*/

int saveasarclm(char *fname,struct arclmframe *af)
/*SAVE AS ARCLM INPUTFILE.*/
{
  FILE *fin;
  char str[256],dandf[256];
  char dir[]=DIRECTORY;
  /*char dir[]="\0";*/
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
  /*strcat(dandf,"cansav.inp");*/
  strcat(dandf,fname);
  fin=fopen(dandf,"w");                                /*SAVE FILE.*/
  if(fin==NULL)
  {
	errormessage("ACCESS IMPOSSIBLE.");
	return 0;
  }

  fprintf(fin,"%5d %5d %5d\n",af->nnode,af->nelem,af->nsect);

  for(i=0;i<(af->nsect);i++)
  {
	fprintf(fin,"%5d %.5E %.5f %.4f %.8f %.8f %.8f",
			(af->sects+i)->code,
			(af->sects+i)->E,
			(af->sects+i)->poi,
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
	  fprintf(fin," %5d",(af->sects+i)->type);
	}
    else
    {
      fprintf(fin," %5d",0);
    }

    fprintf(fin," %5d",(af->sects+i)->ocode);

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
}/*saveasarclm*/

int savebanddecreasedarclm(struct arclmframe *af)
/*SAVE BAND DECREASED ARCLM.*/
{
  FILE *fin;
  char str[256],dandf[256];
  char dir[]=DIRECTORY;
  /*char dir[]="\0";*/
  int i,j,ii,jj,offset,off;
  long int *moff,*noff;

  strcpy(dandf,dir);
  strcat(dandf,"hogtxt.bda");
  fin=fopen(dandf,"w");                                /*SAVE FILE.*/
  if(fin==NULL)
  {
	errormessage("ACCESS IMPOSSIBLE.");
    return 0;
  }

  moff=(long int *)malloc(af->nnode*sizeof(long int));
  noff=(long int *)malloc(af->nnode*sizeof(long int));
  if(moff==NULL || noff==NULL) return 0;
  for(i=0;i<af->nnode;i++)
  {
    *(moff+i)=i;
    *(noff+i)=i;
  }
  decreaseband(moff,noff,af);                /*DECREASE BAND WIDTH.*/

  fprintf(fin,"%5d %5d %5d\n",af->nnode,af->nelem,af->nsect);

  for(i=0;i<(af->nsect);i++)
  {
    fprintf(fin,"%5d %.5E %.5f %.4f %.8f %.8f %.8f",
            (af->sects+i)->code,
            (af->sects+i)->E,
            (af->sects+i)->poi,
            (af->sects+i)->area,
            (af->sects+i)->Ixx,
            (af->sects+i)->Iyy,
            (af->sects+i)->Jzz);
    for(j=0;j<6;j++)
    {
      fprintf(fin," %9.1f %9.1f",(af->sects+i)->fmax[j],
                                 (af->sects+i)->fmin[j]);
    }
    if((af->sects+i)->type!=TYPENULL)
    {
      fprintf(fin," %5d",(af->sects+i)->type);
    }

    fprintf(fin,"\n");
  }
  for(i=0;i<(af->nnode);i++)
  {
    off=*(noff+i);

    fprintf(fin,"%5d %7.3f %7.3f %7.3f\n",
            (af->nodes+off)->code,
            (af->nodes+off)->d[GX],
            (af->nodes+off)->d[GY],
            (af->nodes+off)->d[GZ]);
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
        fprintf(fin," %.3f",(af->elems+i)->stress[ii][jj]);
      }
    }
    fprintf(fin,"\n");
  }

  for(i=0;i<(af->nnode);i++) /*CONF VECTOR:CONFINEMENT,VALUE.*/
  {
    off=*(noff+i);

    fprintf(fin,"%5d",(af->nodes+off)->code);

    for(j=0;j<6;j++)
    {
      offset=6*off+j;
      fprintf(fin," %1d",(af->confs+offset)->iconf);
    }
    for(j=0;j<6;j++)
    {
      offset=6*off+j;
      fprintf(fin," %8.3f",(af->confs+offset)->value);
    }
    fprintf(fin,"\n");
  }
  /*LONGREACTION:UNDER CONSTRUCTION.*/
  fclose(fin);

  sprintf(str,"Saved As %s",dandf);
  MessageBox(NULL,str,"Save",MB_OK);

  return 1;
}/*savebanddecreasedarclm*/

int savebanddecreasedorgan(FILE *fout,
                           struct organ *org,
                           struct viewparam *vp,
                           struct arclmframe *af)
/*DECREASE BAND ARCLM AND SAVE ORGAN INTO OUTPUT FILE.*/
/*CODED BY ARAKI 2008.08.26.*/
{
  int i,j,k,offset,off;
  long int *moff,*noff;
  long int nnode,nelem,nprop,nsect;
  long int count;
  double eps=1.0E-12;

  moff=(long int *)malloc(af/*org*/->nnode*sizeof(long int));
  noff=(long int *)malloc(af/*org*/->nnode*sizeof(long int));
  if(moff==NULL || noff==NULL) return 0;
  for(i=0;i</*af*/org->nnode;i++)
  {
    *(moff+i)=i;
    *(noff+i)=i;
  }
/*  decreasebandorgan(moff,noff,org); */   /*DECREASE BAND WIDTH OF ORGAN.*/
  decreaseband(moff,noff,af);            /*DECREASE BAND WIDTH.*/

  fseek(fout,0L,SEEK_SET);

  fprintf(fout,"\"CREATED ORGAN FRAME.\"\n");

  nnode=org->nnode;
  nelem=org->nelem;
  nprop=org->nprop;
  nsect=org->nsect;

  count=nelem;
  for(i=0;i<nelem;i++)                       /*DELETE BRACE FROM WALL,SLAB*/
  {
    if((org->elems+i)->type==WALL || (org->elems+i)->type==SLAB)
    {
      if(((org->elems+i)->nban)==1)
      {
        count=count-2;
      }
    }
  }
  nelem=count;

  count=nsect;
  for(i=0;i<nsect;i++)
  {
    if((org->sects+i)->code>900) count--;
  }
  nsect=count;

  fprintf(fout,"NNODE %ld\n",nnode);
  fprintf(fout,"NELEM %ld\n",nelem);
  fprintf(fout,"NPROP %ld\n",nprop);
  fprintf(fout,"NSECT %ld\n",nsect);
  fprintf(fout,"\n");

  /*EARTHQUAKE DATA*/
  fprintf(fout,"BASE    %.3f %.3f\n",org->ai.Cox,org->ai.Coy);
  fprintf(fout,"LOCATE  %.3f\n",org->ai.Z);
  fprintf(fout,"TFACT   %.3f\n",org->ai.T1);
  fprintf(fout,"GPERIOD %.3f\n",org->ai.Tc);
  fprintf(fout,"\n");

  /*VIEW DATA*/
  fprintf(fout,"GFACT %.1f\n",vp->gfactor);
  fprintf(fout,"FOCUS %.1f %.1f %.1f\n",vp->focus.d[0],
                                        vp->focus.d[1],
                                        vp->focus.d[2]);
  fprintf(fout,"ANGLE %.1f %.1f\n",vp->phi,vp->theta);
  fprintf(fout,"DISTS %.1f %.1f\n",vp->r,vp->odv);
  fprintf(fout,"\n");

  for(i=0;i<nprop;i++)
  {
    fprintf(fout,"PROP %3ld ",(org->props+i)->code);

    /*if((org->props+i)->name!=NULL)*/
    if(strlen((org->props+i)->name)>0)
    {
      fprintf(fout,"PNAME %s\n",(org->props+i)->name);
    }
    else
    {
      fprintf(fout,"PNAME \"Noname\"\n");
    }

    fprintf(fout,"         HIJU %13.3f\n",(org->props+i)->hiju);
    fprintf(fout,"         E    %13.3f\n",(org->props+i)->E);
    fprintf(fout,"         POI  %13.5f\n",(org->props+i)->poi);
    fprintf(fout,"         PCOLOR %3d %3d %3d\n",(org->props+i)->r,
                                                 (org->props+i)->g,
                                                 (org->props+i)->b);
  }
  fprintf(fout,"\n");

  for(i=0;i<nsect;i++)
  {
    fprintf(fout,"SECT %3ld ",(org->sects+i)->code);

    /*if((org->sects+i)->name!=NULL)*/
    if(strlen((org->sects+i)->name)>0)
    {
      fprintf(fout,"SNAME %s\n",(org->sects+i)->name);
    }
    else
    {
      fprintf(fout,"SNAME \"Noname\"\n");
    }

    if((org->sects+i)->role==ROLEHOJO)
    {
      fprintf(fout,"         SROLE HOJO\n");
    }
    else
    {
      fprintf(fout,"         NFIG %d\n",(org->sects+i)->nfig);
      for(j=0;j<(org->sects+i)->nfig;j++)
      {
        fprintf(fout,"         FIG %3d ",
                ((org->sects+i)->figs+j)->code);
        fprintf(fout,"FPROP %d\n",
                ((org->sects+i)->figs+j)->prop->code);

        if((((org->sects+i)->figs+j)->area<-eps ||
            ((org->sects+i)->figs+j)->area> eps ||
            ((org->sects+i)->figs+j)->Ixx <-eps ||
            ((org->sects+i)->figs+j)->Ixx > eps ||
            ((org->sects+i)->figs+j)->Iyy <-eps ||
            ((org->sects+i)->figs+j)->Iyy > eps ||
            ((org->sects+i)->figs+j)->Jzz <-eps ||
            ((org->sects+i)->figs+j)->Jzz > eps) &&
            ((org->sects+i)->code)<700)
        {
          fprintf(fout,"                 ");
          fprintf(fout,"AREA  %.4f\n",((org->sects+i)->figs+j)->area);
          fprintf(fout,"                 ");
          fprintf(fout,"IXX   %.8f\n",((org->sects+i)->figs+j)->Ixx);
          fprintf(fout,"                 ");
          fprintf(fout,"IYY   %.8f\n",((org->sects+i)->figs+j)->Iyy);
          fprintf(fout,"                 ");
          fprintf(fout,"VEN   %.8f\n",((org->sects+i)->figs+j)->Jzz);
        }
        if((((org->sects+i)->figs+j)->thick<-eps ||
            ((org->sects+i)->figs+j)->thick> eps) &&
            ((org->sects+i)->code)>700)
        {
          fprintf(fout,"                 ");
          fprintf(fout,"THICK %.5f\n",((org->sects+i)->figs+j)->thick);
        }
      }

      if((org->sects+i)->exp!=0.0)
      {
        fprintf(fout,"         EXP %.3f\n",(org->sects+i)->exp);

        fprintf(fout,"         ");
        fprintf(fout,"NZMAX %8.3f ", (org->sects+i)->fmax[0]);
        fprintf(fout,"NZMIN %8.3f\n",(org->sects+i)->fmin[0]);
        fprintf(fout,"         ");
        fprintf(fout,"QXMAX %8.3f ", (org->sects+i)->fmax[1]);
        fprintf(fout,"QXMIN %8.3f\n",(org->sects+i)->fmin[1]);
        fprintf(fout,"         ");
        fprintf(fout,"QYMAX %8.3f ", (org->sects+i)->fmax[2]);
        fprintf(fout,"QYMIN %8.3f\n",(org->sects+i)->fmin[2]);
        fprintf(fout,"         ");
        fprintf(fout,"MZMAX %8.3f ", (org->sects+i)->fmax[3]);
        fprintf(fout,"MZMIN %8.3f\n",(org->sects+i)->fmin[3]);
        fprintf(fout,"         ");
        fprintf(fout,"MXMAX %8.3f ", (org->sects+i)->fmax[4]);
        fprintf(fout,"MXMIN %8.3f\n",(org->sects+i)->fmin[4]);
        fprintf(fout,"         ");
        fprintf(fout,"MYMAX %8.3f ", (org->sects+i)->fmax[5]);
        fprintf(fout,"MYMIN %8.3f\n",(org->sects+i)->fmin[5]);
      }
    }

    if((org->sects+i)->lload[0]!=0.0 ||
       (org->sects+i)->lload[1]!=0.0 ||
       (org->sects+i)->lload[2]!=0.0)
    {
      fprintf(fout,"         LLOAD");
      fprintf(fout," %.3f %.3f %.3f\n",(org->sects+i)->lload[0],
                                       (org->sects+i)->lload[1],
                                       (org->sects+i)->lload[2]);
    }

    if((org->sects+i)->lload[0]!=0.0 ||
       (org->sects+i)->lload[1]!=0.0 ||
       (org->sects+i)->lload[2]!=0.0)
    {
      fprintf(fout,"         PERPL");
      fprintf(fout," %.3f %.3f %.3f\n",(org->sects+i)->perpl[0],
                                       (org->sects+i)->perpl[1],
                                       (org->sects+i)->perpl[2]);
    }

    if((org->sects+i)->dcolor.r!=255 ||
       (org->sects+i)->dcolor.g!=255 ||
       (org->sects+i)->dcolor.b!=255)
    {
      fprintf(fout,"         COLOR");
      fprintf(fout," %d %d %d\n",(org->sects+i)->dcolor.r,
                                 (org->sects+i)->dcolor.g,
                                 (org->sects+i)->dcolor.b);
    }
  }
  fprintf(fout,"\n");

  for(i=0;i<nnode;i++)
  {
    off=*(moff+i);
    fprintf(fout,"NODE %4ld  ",(org->nodes+i)->code);
    fprintf(fout,"CORD %7.3f %7.3f %7.3f  ",
            (org->nodes+off)->d[0],
            (org->nodes+off)->d[1],
            (org->nodes+off)->d[2]);
    fprintf(fout,"ICON %1d %1d %1d %1d %1d %1d  ",
           (org->confs+6*off+0)->iconf,
           (org->confs+6*off+1)->iconf,
           (org->confs+6*off+2)->iconf,
           (org->confs+6*off+3)->iconf,
           (org->confs+6*off+4)->iconf,
           (org->confs+6*off+5)->iconf);

    fprintf(fout,"VCON  ");
    if((org->confs+6*off+0)->value==0.0) fprintf(fout," 0.0  ");
    else fprintf(fout,"%6.3f",(org->confs+6*off+0)->value);
    if((org->confs+6*off+1)->value==0.0) fprintf(fout," 0.0  ");
    else fprintf(fout,"%6.3f",(org->confs+6*off+1)->value);
    if((org->confs+6*off+2)->value==0.0) fprintf(fout," 0.0  ");
    else fprintf(fout,"%6.3f",(org->confs+6*off+2)->value);
    if((org->confs+6*off+3)->value==0.0) fprintf(fout," 0.0  ");
    else fprintf(fout,"%6.3f",(org->confs+6*off+3)->value);
    if((org->confs+6*off+4)->value==0.0) fprintf(fout," 0.0  ");
    else fprintf(fout,"%6.3f",(org->confs+6*off+4)->value);
    if((org->confs+6*off+5)->value==0.0) fprintf(fout," 0.0");
    else fprintf(fout,"%6.3f",(org->confs+6*off+5)->value);
    fprintf(fout,"\n");
  }
  fprintf(fout,"\n");

  for(i=0;i<nelem;i++)
  {
    fprintf(fout,"ELEM %5ld ",(org->elems+i)->code);
    fprintf(fout,"ESECT %3ld ",(org->elems+i)->sect->code);

    if(((org->elems+i)->nnod)==2)
    {
      fprintf(fout,"CANG %.5f ",(org->elems+i)->cangle);
    }

    fprintf(fout,"ENODS %d ",(org->elems+i)->nnod);
    fprintf(fout,"ENOD");
    for(j=0;j<(org->elems+i)->nnod;j++)
    {
      off=*(noff+(*((org->elems+i)->nods+j))->loff);
/*      fprintf(fout," %d",(*((org->elems+i)->nods+j))->code); */
      fprintf(fout," %d",(org->nodes+off)->code);
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
      fprintf(fout,"                     ");
      fprintf(fout,"EBANS %d ",(org->elems+i)->nban);
      for(j=0;j<(org->elems+i)->nban;j++)
      {
        if(j>=1) fprintf(fout,"                             ");

        fprintf(fout,"EBAN %ld ",((org->elems+i)->bans+j)->code);
        fprintf(fout,"BNODS %d ",((org->elems+i)->bans+j)->nnod);
        fprintf(fout,"BNOD");
        for(k=0;k<((org->elems+i)->bans+j)->nnod;k++)
        {
          off=*(noff+(*(((org->elems+i)->bans+j)->nods+k))->loff);
/*          fprintf(fout," %ld",
                  (*(((org->elems+i)->bans+j)->nods+k))->code);   */
          fprintf(fout," %ld",(org->nodes+off)->code);
        }

        fprintf(fout,"\n");
      }
    }

    if(((org->elems+i)->nnod)==2)
    {
//      fprintf(fout,"           CANG %.5f\n",(org->elems+i)->cangle);

      fprintf(fout,"           CMQ ");
      for(j=0;j<2;j++)
      {
        for(k=0;k<6;k++)
        {
     /*     fprintf(fout," %.1f",(org->elems+i)->initial[j][k]);   */
          fprintf(fout," %.1f",0.0);
        }
        if(j<1) fprintf(fout," ");
      }
      fprintf(fout,"\n");
    }

    if((org->elems+i)->type==COLUMN)
    {
      fprintf(fout,"           TYPE COLUMN\n");
    }
    if((org->elems+i)->type==GIRDER)
    {
      fprintf(fout,"           TYPE GIRDER\n");
    }
    if((org->elems+i)->type==BEAM)
    {
      fprintf(fout,"           TYPE BEAM\n");
    }
    if((org->elems+i)->type==WALL)
    {
      fprintf(fout,"           TYPE WALL\n");
    }
    if((org->elems+i)->type==SLAB)
    {
      fprintf(fout,"           TYPE SLAB\n");
    }
    if((org->elems+i)->type==BRACE)
    {
      fprintf(fout,"           TYPE BRACE\n");
    }
    if((org->elems+i)->type==NULL)
    {
      fprintf(fout,"           TYPE NULL\n");
    }

    if((org->elems+i)->lface[0]!=0.0 ||
       (org->elems+i)->lface[1]!=0.0)
    {
      fprintf(fout,"           LFACE %.3f %.3f\n",
              (org->elems+i)->lface[0],(org->elems+i)->lface[1]);
    }
    if((org->elems+i)->hface[0]!=0.0 ||
       (org->elems+i)->hface[1]!=0.0)
    {
      fprintf(fout,"           HFACE %.3f %.3f\n",
              (org->elems+i)->hface[0],(org->elems+i)->hface[1]);
    }
    if((org->elems+i)->wrect[0]!=0.0 ||
       (org->elems+i)->wrect[1]!=0.0)
    {
      fprintf(fout,"           WRECT %.3f %.3f\n",
              (org->elems+i)->wrect[0],(org->elems+i)->wrect[1]);
    }

    if((org->elems+i)->role==ROLEWEIGHT)
    {
      fprintf(fout,"           ROLE W\n");
    }
    if((org->elems+i)->role==ROLERIGID)
    {
      fprintf(fout,"           ROLE R\n");
    }
    if((org->elems+i)->role==ROLESTRESS)
    {
      fprintf(fout,"           ROLE S\n");
    }
  }

  return 1;
}/*savebanddecreasedorgan*/

void initialform(struct onode *nodes,double *ddisp,int nnode)
/*INITIAL FORMATION INTO DISPLACEMENT.WITHOUT NODE CODE.*/
{
  int i;

  for(i=0;i<nnode;i++)
  {
	*(ddisp+6*i+0)=(nodes+i)->d[0];
	*(ddisp+6*i+1)=(nodes+i)->d[1];
	*(ddisp+6*i+2)=(nodes+i)->d[2];
	*(ddisp+6*i+3)=(nodes+i)->r[0];
	*(ddisp+6*i+4)=(nodes+i)->r[1];
	*(ddisp+6*i+5)=(nodes+i)->r[2];
  }

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




void initialelem(struct owire *elems,
				 struct memoryelem *melem,int nelem)
/*ASSEMBLAGE LONG STRESS.*/
{
  int i,j;

  for(i=0;i<nelem;i++)
  {
	(melem+i)->code=(elems+i)->code;                       /*CODE.*/

	for(j=0;j<6;j++)                                 /*ICONF[2][6].*/
	{
	  (melem+i)->bond[0][j]=(elems+i)->iconf[0][j];
	  (melem+i)->bond[1][j]=(elems+i)->iconf[1][j];
	}
	for(j=0;j<6;j++)                           /*LONG STRESS[2][6].*/
	{
	  (melem+i)->stress[0][j]=(elems+i)->stress[0][j];
	  (melem+i)->stress[1][j]=(elems+i)->stress[1][j];
	}
  }

  return;
}/*initialelem*/





void initialreact(FILE *fin,double *dreact,int nreact)
/*ASSEMBLAGE LONG REACTIONS.*/
{
  char **data;
  int i,n;
  double ddata;

  for(i=0;i<nreact;i++)
  {
    if(fin==NULL) *(dreact+i)=0.0;
    else
    {
      data=fgetsbrk(fin,&n);
      if(n<=2) ddata=0.0;
      else     ddata=strtod(*(data+2),NULL);  /*0:INODE 1:DIRECTION*/

      *(dreact+i)=ddata;

      for(;n>0;n--) free(*(data+n-1));
      free(data);
    }
  }

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

void inputinitII(FILE *fin,int *nnode,int *nelem,int *nshell,int *nsect,int *nconstraint)
{
  char **data;
  int n;

  data=fgetsbrk(fin,&n);
  if(n==5)
  {
	*nnode=strtol(*(data+0),NULL,10);
	*nelem=strtol(*(data+1),NULL,10);
	*nshell=strtol(*(data+2),NULL,10);
	*nsect=strtol(*(data+3),NULL,10);
	*nconstraint=strtol(*(data+4),NULL,10);
  }
  else if(n==4)
  {
	*nnode=strtol(*(data+0),NULL,10);
	*nelem=strtol(*(data+1),NULL,10);
	*nshell=strtol(*(data+2),NULL,10);
	*nsect=strtol(*(data+3),NULL,10);
	*nconstraint=0;
  }
  else
  {
	*nnode=strtol(*(data+0),NULL,10);
	*nelem=strtol(*(data+1),NULL,10);
	*nshell=0;
	*nsect=strtol(*(data+2),NULL,10);
	*nconstraint=0;
  }
  for(;n>0;n--) free(*(data+n-1));
  free(data);

  return;
}/*inputinitII*/

void inputframeinit(FILE *fin,int *nnode,int *nelem,int *nsect)
/*INPUT INITIAL DATA FROM FRAME3 INPUTFILE.*/
{
  char **data;
  int n;

  data=fgetsbrk(fin,&n);
  *nnode=strtol(*(data+0),NULL,10);
  *nelem=strtol(*(data+1),NULL,10);
  *nsect=strtol(*(data+12),NULL,10);
  for(;n>0;n--) free(*(data+n-1));
  free(data);

  return;
}/*inputframeinit*/

void inputnode(double *ddisp,struct onode *node)
{
  long int loff;

  loff=6*(node->loff);

  node->d[0]=*(ddisp+loff+0);
  node->d[1]=*(ddisp+loff+1);
  node->d[2]=*(ddisp+loff+2);
  node->r[0]=*(ddisp+loff+3);
  node->r[1]=*(ddisp+loff+4);
  node->r[2]=*(ddisp+loff+5);

  return;
}/*inputnode*/

void inputelem(struct owire *elems,
			   struct memoryelem *melem,int offset,
			   struct owire *elem)
{
  int i;

  elem->nnod=2;                                      /*KUTOMI*/
  elem->code=(elems+offset)->code;                  /*ELEMENT CODE.*/
  elem->loff=offset;
  elem->sect=(elems+offset)->sect;               /*SECTION POINTER.*/

  elem->node[0]=(elems+offset)->node[0];            /*HEAD POINTER.*/
  elem->node[1]=(elems+offset)->node[1];            /*TAIL POINTER.*/
  elem->cangle=(elems+offset)->cangle;               /*COORD ANGLE.*/

  elem->Ee[0]=(elems+offset)->Ee[0];              /*ELASTIC ENERGY.*/
  elem->Ee[1]=(elems+offset)->Ee[1];
  elem->Ep[0]=(elems+offset)->Ep[0];              /*PLASTIC ENERGY.*/
  elem->Ep[1]=(elems+offset)->Ep[1];

  if(melem!=NULL)
  {
	for(i=0;i<6;i++)                                    /*BOUNDARY.*/
    {
	  elem->iconf[0][i]=(melem+offset)->bond[0][i];
	  elem->iconf[1][i]=(melem+offset)->bond[1][i];
    }
	for(i=0;i<6;i++)                                      /*STRESS.*/
	{
	  elem->stress[0][i]=(melem+offset)->stress[0][i];
	  elem->stress[1][i]=(melem+offset)->stress[1][i];
	}
  }
  else
  {
	for(i=0;i<6;i++)
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
  char str[256];

  data=fgetsbrk(fin,&n);

  sect->code=strtol(*(data+0),NULL,10);             /*SECTION CODE.*/
  sect->E=strtod(*(data+1),NULL);                /*YOUNG'S MODULUS.*/
  sect->poi=strtod(*(data+2),NULL);              /*POISSON'S RATIO.*/
  sect->area=strtod(*(data+3),NULL);                        /*AREA.*/
  sect->Ixx=strtod(*(data+4),NULL);
  sect->Iyy=strtod(*(data+5),NULL);
  sect->Jzz=strtod(*(data+6),NULL);          /*ST.VENANT'S TORTION.*/
  k=7;
  for(i=0;i<6;i++)           /*UPPER,LOWER LIMIT OF YIELD SURFACE.*/
  {
	sect->fmax[i]=strtod(*(data+k),NULL); k++;
	sect->fmin[i]=strtod(*(data+k),NULL); k++;
  }
  if(n==20 || n==21)
  {
	  sect->hiju[0] = 0.0;
	  sect->hiju[1] = 0.0;
	  sect->hiju[2] = 0.0;
	  sect->lload[0] = 0.0;
	  sect->lload[1] = 0.0;
	  sect->lload[2] = 0.0;
	  sect->perpl[0] = 0.0;
	  sect->perpl[1] = 0.0;
	  sect->perpl[2] = 0.0;
	  sect->dcolor.r = 255;
	  sect->dcolor.g = 255;
	  sect->dcolor.b = 255;
	  sect->type=strtol(*(data+19),NULL,10);   /*SECTION TYPE.*/
  }
  else if(n==32)
  {
	  sect->hiju[0] = strtod(*(data+19),NULL);
	  sect->hiju[1] = strtod(*(data+20),NULL);
	  sect->hiju[2] = strtod(*(data+21),NULL);
	  sect->lload[0] = strtod(*(data+22),NULL);
	  sect->lload[1] = strtod(*(data+23),NULL);
	  sect->lload[2] = strtod(*(data+24),NULL);
	  sect->perpl[0] = strtod(*(data+25),NULL);
	  sect->perpl[1] = strtod(*(data+26),NULL);
	  sect->perpl[2] = strtod(*(data+27),NULL);
	  sect->dcolor.r = strtol(*(data+28),NULL,10);
	  sect->dcolor.g = strtol(*(data+29),NULL,10);
	  sect->dcolor.b = strtol(*(data+30),NULL,10);
	  sect->type=strtol(*(data+31),NULL,10);   /*SECTION TYPE.*/
  }
  else if(n==35)
  {
	  sect->hiju[0] = strtod(*(data+19),NULL);
	  sect->hiju[1] = strtod(*(data+20),NULL);
	  sect->hiju[2] = strtod(*(data+21),NULL);
	  sect->lload[0] = strtod(*(data+22),NULL);
	  sect->lload[1] = strtod(*(data+23),NULL);
	  sect->lload[2] = strtod(*(data+24),NULL);
	  sect->perpl[0] = strtod(*(data+25),NULL);
	  sect->perpl[1] = strtod(*(data+26),NULL);
	  sect->perpl[2] = strtod(*(data+27),NULL);
	  sect->dcolor.r = strtol(*(data+28),NULL,10);
	  sect->dcolor.g = strtol(*(data+29),NULL,10);
	  sect->dcolor.b = strtol(*(data+30),NULL,10);
	  sect->type=strtol(*(data+31),NULL,10);   /*SECTION TYPE.*/
	  sect->yieldinit = strtod(*(data+32),NULL);
	  sect->yieldcoefficient = strtod(*(data+33),NULL);
	  sect->yieldpower = strtod(*(data+34),NULL);
  }
  else
  {
	  sect->hiju[0] = 0.0;
	  sect->hiju[1] = 0.0;
	  sect->hiju[2] = 0.0;
	  sect->lload[0] = 0.0;
	  sect->lload[1] = 0.0;
	  sect->lload[2] = 0.0;
	  sect->perpl[0] = 0.0;
	  sect->perpl[1] = 0.0;
	  sect->perpl[2] = 0.0;
	  sect->dcolor.r = 255;
	  sect->dcolor.g = 255;
	  sect->dcolor.b = 255;
	  sect->type=TYPENULL;
  }
  //sprintf(str,"%d",n);
  //errormessage(str);

  for(;n>0;n--) free(*(data+n-1));
  free(data);



  return;
}/*readsect*/

void readsrcanrate(FILE *fin,struct arclmframe *af)
{
  char **data;
  int nstr;
  long int i,j,code,loff;
  double *srate;
  long int countsect=0; /****SRCANMAX****/

  /****SRCANMAX****/
  af->nosect=0;
  af->iosect=(long int *)malloc((af->nelem)*sizeof(long int));
  af->eosect=(int *)malloc((af->nelem)*sizeof(int));
  af->dosect=(double *)malloc((af->nelem)*sizeof(double));
  for(i=0;i<(af->nelem);i++)
  {
	*(af->iosect+i)=0;
	*(af->eosect+i)=(af->nelem)+1;
	*(af->dosect+i)=0.0;
  }
  /****************/

  srate=(double *)malloc(af->nsect*sizeof(double));

  for(i=0;i<(af->nsect);i++) /*INITIAL*/
  {
    *(srate+i)=0.0;
    (af->sects+i)->smax=0;
  }

  for(i=0;i<(af->nelem);i++) /*INITIAL*/
  {
    (af->elems+i)->srate[0]=0.0;
    (af->elems+i)->srate[1]=0.0;
    (af->elems+i)->srate[2]=0.0;
    (af->elems+i)->srate[3]=0.0;
  }

  while(1)
  {
    data=fgetsbrk(fin,&nstr);
    if(nstr<2)
    {
      for(;nstr>0;nstr--) free(*(data+nstr-1));
      free(data);
      free(srate);
      break;
    }

    code=strtol(*(data+1),NULL,10);

    for(i=0;i<(af->nelem);i++) /*SEARCH*/
    {
      if(code==(af->elems+i)->code)
      {
        (af->elems+i)->srate[0]=strtod(*(data+4),NULL);
        (af->elems+i)->srate[1]=strtod(*(data+5),NULL);

        if(nstr>=7)
        {
          (af->elems+i)->srate[2]=strtod(*(data+6),NULL);
          (af->elems+i)->srate[3]=strtod(*(data+7),NULL);
        }

		/****SRCANMAX****/
		if(nstr>=9)
		{
		  (af->elems+i)->sect->ocode=strtod(*(data+9),NULL);
		}

		for(j=0;j<(af->nosect);j++)
		{
		  if(((af->elems+i)->sect->ocode)==*(af->iosect+j))
		  {
			break;
		  }
		}
		if(j==(af->nosect))
		{
		  af->nosect++;
		  *(af->iosect+j)=(af->elems+i)->sect->ocode;
		  *(af->eosect+j)=i;
		}
		loff=j;
		for(j=0;j<4;j++)
		{
		  if(*(af->dosect+loff)<(af->elems+i)->srate[j])
		  {
			*(af->dosect+loff)=(af->elems+i)->srate[j];
			*(af->eosect+loff)=i; /*OFFSET OF MAX ELEMENT*/
		  }
		}

		/*loff=(af->elems+i)->sect->loff;*/
//		for(j=0;j<4;j++)
//		{
//		  if(*(srate+loff)<(af->elems+i)->srate[j])
//		  {
//			*(srate+loff)=(af->elems+i)->srate[j];
//			(af->sects+loff)->smax=i; /*OFFSET OF MAX ELEMENT*/
//		  }
//		}
		/****************/


		break;
	  }
	}

	for(;nstr>0;nstr--) free(*(data+nstr-1));
	free(data);
  }
  return;
}/*readsrcanrate*/

double **directioncosine(double x1,double y1,double z1,
						 double x2,double y2,double z2,
						 double cangle)
{
  int i;
  double dl1,dl2,dx,dy,dz;
  double c0,c1,c2;
  double ol2,ox,oy,oz,al1,xx,xy,xz,yx,yy,yz; /*AXIS VECTOR*/
  double **drccos;

  double theta; /*BY ARAKI FOR HAKUSHIMA,STICK*/

  if(!strncmp(prj,"sunny",5))       theta=0.5;
  else if(!strncmp(prj,"stick",5))  theta=0.9;
  else                              theta=0.1;

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
  /*else if((dl2/dl1)<0.0001)*/
  /*else if((dl2/dl1)<0.1)*/ /*DEFORMED COLUMN*//*UNDER CONSTRUCTION*/
#if 0
  else if((dl2/dl1)<theta) /*DEFORMED COLUMN*/
  {
//    MessageBox(NULL,"There is a vertical member.","NOTE",MB_OK); //honda
	ox=dx/dl1;
	oy=dy/dl1;
	oz=dz/dl1;
	ol2=ox*ox+oy*oy;
	*(*(drccos+0)+0)=dx/dl1; /*AXIS EZ*/
	*(*(drccos+0)+1)=dy/dl1;
	*(*(drccos+0)+2)=dz/dl1;

	if(oz>=0.0)
	{
	  xx=(ox*oy*(oz-1.0)*cos(cangle)
		  -(ox*ox*oz+oy*oy)*sin(cangle))/ol2; /*AXIS EX*/
	  xy=((ox*ox+oy*oy*oz)*cos(cangle)
		  -ox*oy*(oz-1.0)*sin(cangle))/ol2;
	}
	else
	{
	  xx=(ox*oy*(-oz-1.0)*cos(cangle)
		  -(-ox*ox*oz+oy*oy)*sin(cangle))/ol2; /*AXIS EX*/
	  xy=((ox*ox-oy*oy*oz)*cos(cangle)
		  -ox*oy*(-oz-1.0)*sin(cangle))/ol2;
	}
	xz=ox*sin(cangle)-oy*cos(cangle);
	al1=sqrt(xx*xx+xy*xy+xz*xz);
	*(*(drccos+1)+0)=xx/al1;
	*(*(drccos+1)+1)=xy/al1;
	*(*(drccos+1)+2)=xz/al1;

	yx=oy*xz-oz*xy; /*AXIS EY=EZxEX*/
	if(oz>=0.0) yy=oz*xx-ox*xz;
	else        yy=oz*xx+ox*xz;
	yz=ox*xy-oy*xx;
	al1=sqrt(yx*yx+yy*yy+yz*yz);
	*(*(drccos+2)+0)=yx/al1;
	*(*(drccos+2)+1)=yy/al1;
	*(*(drccos+2)+2)=yz/al1;
  }
#endif
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

void directioncosineII(double x1,double y1,double z1,
                       double x2,double y2,double z2,
                       double cangle,
                       double **drccos)
/*DIRECTION COSINE WITHOUT MALLOC.*/
{
  double dl1,dl2,dx,dy,dz;
  double c0,c1,c2;
  double ol2,ox,oy,oz,al1,xx,xy,xz,yx,yy,yz; /*AXIS VECTOR*/

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
  /* else if((dl2/dl1)<0.0001) #<{(|DEFORMED COLUMN|)}>##<{(|UNDER CONSTRUCTION|)}># */

//  else if((dl2/dl1)<0.0001) /*DEFORMED COLUMN*//*UNDER CONSTRUCTION*/    //original:0.1 for sydney:0.0001
  else if((dl2/dl1)<0.1) /*DEFORMED COLUMN*//*UNDER CONSTRUCTION*/    //original:0.1 for sydney:0.0001
  {
	ox=dx/dl1;
	oy=dy/dl1;
	oz=dz/dl1;
	ol2=ox*ox+oy*oy;
	*(*(drccos+0)+0)=dx/dl1; /*AXIS EZ*/
	*(*(drccos+0)+1)=dy/dl1;
	*(*(drccos+0)+2)=dz/dl1;

	//MessageBox(NULL,"There is a vertical member.","NOTE",MB_OK);//honda


	if(oz>=0.0)
    {
      xx=(ox*oy*(oz-1.0)*cos(cangle)
          -(ox*ox*oz+oy*oy)*sin(cangle))/ol2; /*AXIS EX*/
      xy=((ox*ox+oy*oy*oz)*cos(cangle)
          -ox*oy*(oz-1.0)*sin(cangle))/ol2;
    }
    else
    {
      xx=(ox*oy*(-oz-1.0)*cos(cangle)
          -(-ox*ox*oz+oy*oy)*sin(cangle))/ol2; /*AXIS EX*/
	  xy=((ox*ox-oy*oy*oz)*cos(cangle)
          -ox*oy*(-oz-1.0)*sin(cangle))/ol2;
    }
    xz=ox*sin(cangle)-oy*cos(cangle);
    al1=sqrt(xx*xx+xy*xy+xz*xz);
    *(*(drccos+1)+0)=xx/al1;
    *(*(drccos+1)+1)=xy/al1;
    *(*(drccos+1)+2)=xz/al1;

    yx=oy*xz-oz*xy; /*AXIS EY=EZxEX*/
    if(oz>=0.0) yy=oz*xx-ox*xz;
    else        yy=oz*xx+ox*xz;
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

  return;
}/*directioncosineII*/

double **filmdrccos
(struct onode n1,struct onode n2,struct onode n3)
/*RETURN:FILM DIRECTION COSINE.*/
{
  char str[256];
  int i;
  double dl,Xx,Yx,Zx,Xy,Yy,Zy,Xz,Yz,Zz;
  double **drccos;

  drccos=(double **)malloc(3*sizeof(double *));
  for(i=0;i<=2;i++) *(drccos+i)=(double *)malloc(3*sizeof(double));

  Xx=n2.d[GX]-n1.d[GX]; /*AXIS x=FIRST LINE.*/
  Yx=n2.d[GY]-n1.d[GY];
  Zx=n2.d[GZ]-n1.d[GZ];
  dl=sqrt(Xx*Xx+Yx*Yx+Zx*Zx);
if(dl==0.0)
{
  sprintf(str,"N1=%ld N2=%ld N3=%ld",n1.code,n2.code,n3.code);
  MessageBox(NULL,str,"Filmdrccos",MB_OK);
}
  *(*(drccos+EX)+0)=Xx/dl;
  *(*(drccos+EX)+1)=Yx/dl;
  *(*(drccos+EX)+2)=Zx/dl;

  Xz=Yx*(n3.d[GZ]-n1.d[GZ])-Zx*(n3.d[GY]-n1.d[GY]); /*AXIS z=OUTER.*/
  Yz=Zx*(n3.d[GX]-n1.d[GX])-Xx*(n3.d[GZ]-n1.d[GZ]);
  Zz=Xx*(n3.d[GY]-n1.d[GY])-Yx*(n3.d[GX]-n1.d[GX]);
  dl=sqrt(Xz*Xz+Yz*Yz+Zz*Zz);
if(dl==0.0)
{
  sprintf(str,"N1=%ld N2=%ld N3=%ld",n1.code,n2.code,n3.code);
  MessageBox(NULL,str,"Filmdrccos",MB_OK);
}
  *(*(drccos+EZ)+0)=Xz/dl;
  *(*(drccos+EZ)+1)=Yz/dl;
  *(*(drccos+EZ)+2)=Zz/dl;

  Xy=Yz*Zx-Zz*Yx; /*AXIS y=OUTER.*/
  Yy=Zz*Xx-Xz*Zx;
  Zy=Xz*Yx-Yz*Xx;
  dl=sqrt(Xy*Xy+Yy*Yy+Zy*Zy);
if(dl==0.0)
{
  sprintf(str,"N1=%ld N2=%ld N3=%ld",n1.code,n2.code,n3.code);
  MessageBox(NULL,str,"Filmdrccos",MB_OK);
}
  *(*(drccos+EY)+0)=Xy/dl;
  *(*(drccos+EY)+1)=Yy/dl;
  *(*(drccos+EY)+2)=Zy/dl;

  return drccos;
}/*filmdrccos*/

double **transmatrix(double **drccos)
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

void transmatrixII(double **drccos,double **t)
/*ASSEMBLAGE TRANSFORMATION MATRIX WITHOUT MALLOC.*/
{
  int i,j,n;

  for(n=0;n<4;n++)
  {
    for(i=0;i<3;i++)
    {
      for(j=0;j<12;j++) *(*(t+3*n+i)+j)=0.0;
	  for(j=0;j<3;j++)  *(*(t+3*n+i)+3*n+j)=*(*(drccos+i)+j);
    }
  }
  return;
}/*transmatrixII*/

double **transmatrixIII(double **drccos,int nnod)
/*ASSEMBLAGE TRANSFORMATION MATRIX.*/
{
  int i,j,n;
  double **t;

  t=(double **)malloc(6*nnod*sizeof(double *));
  for(n=1;n<=2*nnod;n++)
  {
	for(i=1;i<=3;i++)
	{
	  *(t+3*(n-1)+i-1)=(double *)malloc(6*nnod*sizeof(double));
	  for(j=1;j<=6*nnod;j++)
	  {*(*(t+3*(n-1)+i-1)+j-1)=0.0;}
	  for(j=1;j<=3;j++)
	  {*(*(t+3*(n-1)+i-1)+3*(n-1)+j-1)=*(*(drccos+i-1)+j-1);}
	}
  }
  return t;
}/*transmatrixIII*/


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

void assememtxII(struct owire elem,double **e)
/*ASSEMBLAGE ELASTIC MATRIX WITHOUT MALLOC.*/
{
  int i,j;
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

  for(i=0;i<12;i++)
  {
	for(j=0;j<12;j++) *(*(e+i)+j)=0.0; /*INITIAL.*/
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
  return;
}/*assememtxII*/

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

  if(A!=0.0) gam=-(Ixx+Iyy)/A;

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
  char str[80];
  signed char iconf[12];
  int i,j;
  double p[12][12],f[2],dfdp[2][6],q[12][2],a[2][2],det;

  for(i=0;i<2;i++)                    /*ICONF[2][6] INTO ICONF[12].*/
  {
	for(j=0;j<6;j++) iconf[6*i+j]=elem.iconf[i][j];
  }

  coefficients(elem,estiff,f,dfdp,q,a);

  for(i=0;i<12;i++)
  {
	for(j=0;j<12;j++) p[i][j]=0.0;                       /*INITIAL.*/
  }


  if(a[0][0]!=0.0 && a[1][1]==0.0)         /*IF I:PLASTIC J:ELASTIC*/
  {
	for(i=0;i<12;i++)
	{
	  if(iconf[i]!=1 && iconf[i]!=-2 && iconf[i]!=-3)
	  {
		for(j=0;j<12;j++)
		{
		  if(iconf[j]!=1 && iconf[j]!=-2 && iconf[j]!=-3)
		  {
			p[i][j]=-1.0/a[0][0]*q[i][0]*q[j][0];
		  }
		}
	  }
	}
  }


  if(a[0][0]==0.0 && a[1][1]!=0.0)         /*IF I:ELASTIC J:PLASTIC*/
  {
	for(i=0;i<12;i++)
    {
      if(iconf[i]!=1 && iconf[i]!=-2 && iconf[i]!=-3)
      {
		for(j=0;j<12;j++)
        {
          if(iconf[j]!=1 && iconf[j]!=-2 && iconf[j]!=-3)
          {
			p[i][j]=-1.0/a[1][1]*q[i][1]*q[j][1];
		  }
		}
	  }
	}
  }


  if((a[0][0]!=0.0)&&(a[1][1]!=0.0))       /*IF I:PLASTIC J:PLASTIC*/
  {
    for(i=0;i<12;i++)
    {
	  if(iconf[i]!=1 && iconf[i]!=-2 && iconf[i]!=-3)
	  {
		for(j=0;j<12;j++)
		{
		  if(iconf[j]!=1 && iconf[j]!=-2 && iconf[j]!=-3)
		  {
			if((a[0][0]*a[1][1]-a[0][1]*a[1][0])==0.0)
			{

				if(globalfile!=NULL)
				{
				  fprintf(globalfile,"Assem : Matrix Singular.ELEM=%d SECT=%d i=%d j=%d\n",
						  elem.code,elem.sect->code,i,j);
				  fprintf(globalfile,"%10.5f x %10.5f - %10.5f x %10.5f = 0.0\n",
						  a[0][0],a[1][1],a[0][1],a[1][0]);
				}

              errormessage("ASSEMPMTX:UNDER CONSIDERATION.");
              sprintf(str,"Matrix Singular.ELEM=%d SECT=%d i=%d j=%d",
                      elem.code,elem.sect->code,i,j);

              /*DELETE THIS ELEMENT*/
			  *(*(estiff+i)+j)=0.0;
              p[i][j]=0.0;
            }
            else
			{
			  det=-1.0/(a[0][0]*a[1][1]-a[0][1]*a[1][0]);
              p[i][j]=det*(a[1][1]*q[i][0]*q[j][0]
                          -a[0][1]*q[i][0]*q[j][1]
                          -a[1][0]*q[i][1]*q[j][0]
                          +a[0][0]*q[i][1]*q[j][1]);
			}
          }
        }
	  }
    }
  }

  for(i=0;i<12;i++)                                     /*ADDITION.*/
  {
    if(iconf[i]!=1 && iconf[i]!=-2 && iconf[i]!=-3)
    {
	  for(j=0;j<12;j++)
	  {
		if(iconf[j]!=1 && iconf[j]!=-2 && iconf[j]!=-3)
        {
		  *(*(estiff+i)+j)+=p[i][j];                /*[k]=[ke]+[kp]*/
        }
      }
	}
  }

  return estiff;
}/*assempmtx*/

double **assempmtxbc(struct owire elem,double **estiff,double ncr)
/*ASSEMBLAGE PLASTIC MATRIX.*/
{
  char str[80];
  signed char iconf[12];
  int i,j;
  double p[12][12],f[2],dfdp[2][6],q[12][2],a[2][2],det;

  for(i=0;i<2;i++)                    /*ICONF[2][6] INTO ICONF[12].*/
  {
	for(j=0;j<6;j++) iconf[6*i+j]=elem.iconf[i][j];
  }

  coefficientsbc(elem,estiff,f,dfdp,q,a,ncr);

  for(i=0;i<12;i++)
  {
	for(j=0;j<12;j++) p[i][j]=0.0;                       /*INITIAL.*/
  }

  if(a[0][0]!=0.0 && a[1][1]==0.0)         /*IF I:PLASTIC J:ELASTIC*/
  {
	for(i=0;i<12;i++)
    {
      if(iconf[i]!=1 && iconf[i]!=-2 && iconf[i]!=-3)
	  {
        for(j=0;j<12;j++)
        {
		  if(iconf[j]!=1 && iconf[j]!=-2 && iconf[j]!=-3)
          {
            p[i][j]=-1.0/a[0][0]*q[i][0]*q[j][0];
		  }
        }
      }
	}
  }
  if(a[0][0]==0.0 && a[1][1]!=0.0)         /*IF I:ELASTIC J:PLASTIC*/
  {
    for(i=0;i<12;i++)
    {
	  if(iconf[i]!=1 && iconf[i]!=-2 && iconf[i]!=-3)
      {
        for(j=0;j<12;j++)
		{
          if(iconf[j]!=1 && iconf[j]!=-2 && iconf[j]!=-3)
		  {
			p[i][j]=-1.0/a[1][1]*q[i][1]*q[j][1];
          }
		}
	  }
	}
  }
  if((a[0][0]!=0.0)&&(a[1][1]!=0.0))       /*IF I:PLASTIC J:PLASTIC*/
  {
	for(i=0;i<12;i++)
	{
	  if(iconf[i]!=1 && iconf[i]!=-2 && iconf[i]!=-3)
	  {
		for(j=0;j<12;j++)
		{
		  if(iconf[j]!=1 && iconf[j]!=-2 && iconf[j]!=-3)
		  {
            if((a[0][0]*a[1][1]-a[0][1]*a[1][0])==0.0)
			{

if(globalfile!=NULL)
{
  fprintf(globalfile,"Assem : Matrix Singular.ELEM=%d SECT=%d i=%d j=%d\n",
		  elem.code,elem.sect->code,i,j);
  fprintf(globalfile,"%10.5f x %10.5f - %10.5f x %10.5f = 0.0\n",
		  a[0][0],a[1][1],a[0][1],a[1][0]);
}

			  errormessage("ASSEMPMTX:UNDER CONSIDERATION.");
              sprintf(str,"Matrix Singular.ELEM=%d SECT=%d i=%d j=%d",
					  elem.code,elem.sect->code,i,j);
/*
MessageBox(NULL,str,"Assempmtx",MB_OK);
*/
			  /*DELETE THIS ELEMENT*/
			  *(*(estiff+i)+j)=0.0;
			  p[i][j]=0.0;
			}
            else
			{
			  det=-1.0/(a[0][0]*a[1][1]-a[0][1]*a[1][0]);
			  p[i][j]=det*(a[1][1]*q[i][0]*q[j][0]
						  -a[0][1]*q[i][0]*q[j][1]
						  -a[1][0]*q[i][1]*q[j][0]
						  +a[0][0]*q[i][1]*q[j][1]);
			}
		  }
		}
	  }
	}
  }

  for(i=0;i<12;i++)                                     /*ADDITION.*/
  {
	if(iconf[i]!=1 && iconf[i]!=-2 && iconf[i]!=-3)
	{
	  for(j=0;j<12;j++)
	  {
		if(iconf[j]!=1 && iconf[j]!=-2 && iconf[j]!=-3)
		{
		  *(*(estiff+i)+j)+=p[i][j];                /*[k]=[ke]+[kp]*/
		}
	  }
	}
  }

  return estiff;
}/*assempmtxbc*/

void coefficients(struct owire elem,double **estiff,
				  double f[],double dfdp[][6],
				  double q[][2],double a[][2])
/*ASSEMBLAGE PLASTIC COEFFICIENTS.*/
{
  signed char iconf[12];
  int i,j,ii,jj;
  double fc[6],fu[6];
  double unit,value;

  for(i=0;i<2;i++)                    /*ICONF[2][6] INTO ICONF[12].*/
  {
	for(j=0;j<6;j++) iconf[6*i+j]=elem.iconf[i][j];
  }

  for(i=0;i<6;i++)                 /*CENTER,WIDTH OF YIELD SURFACE.*/
  {
	fc[i]=0.5*(elem.sect->fmax[i]+elem.sect->fmin[i]);     /*CENTER*/
	fu[i]=0.5*(elem.sect->fmax[i]-elem.sect->fmin[i]);      /*WIDTH*/
  }


  for(i=0;i<2;i++)    /*ASSEMBLAGE YIELD FUNCTION:f,VECTOR:{df/dp}.*/
  {
  /*f=POW(|Qx/Qxu|,EXP)+POW(|Qy/Qyu|,EXP)+POW(|Nz/Nzu|,EXP)*/
  /* +POW(|Mx/Mxu|,EXP)+POW(|My/Myu|,EXP)+POW(|Mz/Mzu|,EXP)*/

	f[i]=0.0;
	for(j=0;j<6;j++)
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

  for(i=0;i<12;i++)                         /*ASSEMBLAGE VECTOR:{q}*/
  {
	if(iconf[i]!=1)
	{
	  for(j=0;j<2;j++)
	  {
		q[i][j]=0.0;

		for(jj=0;jj<6;jj++)
		{
		  if(iconf[6*j+jj]!=1) /*if(iconf[6*j+jj]==-1)*/
		  {
			q[i][j]+=*(*(estiff+i)+6*j+jj)*dfdp[j][jj];
		  }
		}
	  }
	}
  }
  for(i=0;i<2;i++)                     /*INNER PRODUCT a={df/dp}{q}*/
  {
	for(j=0;j<2;j++)
	{
	  a[i][j]=0.0;

	  for(ii=0;ii<6;ii++)
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


#if 0
void coefficients(struct owire elem,double **estiff,
				  double f[],double dfdp[][6],
				  double q[][2],double a[][2])
/*ASSEMBLAGE PLASTIC COEFFICIENTS.*/
{
  signed char iconf[12];
  int i,j,ii,jj;
  double fc[6],fu[6],unit,value;

  for(i=0;i<2;i++)                    /*ICONF[2][6] INTO ICONF[12].*/
  {
	for(j=0;j<6;j++) iconf[6*i+j]=elem.iconf[i][j];
  }

/*
if(globalfile!=NULL) fprintf(globalfile,"\nELEM %d\n",elem.code);
if(globalfile!=NULL) fprintf(globalfile,"   CENTER fc   WIDTH fu\n");
*/

  for(i=0;i<6;i++)                 /*CENTER,WIDTH OF YIELD SURFACE.*/
  {
	fc[i]=0.5*(elem.sect->fmax[i]+elem.sect->fmin[i]);     /*CENTER*/
	fu[i]=0.5*(elem.sect->fmax[i]-elem.sect->fmin[i]);      /*WIDTH*/

/*if(globalfile!=NULL) fprintf(globalfile,"%d %10.5f %10.5f\n",
							 i,fc[i],fu[i]);*/
  }

/*
if(globalfile!=NULL)
fprintf(globalfile,"        STRESS       RATE       dfdp\n");
*/

  for(i=0;i<2;i++)    /*ASSEMBLAGE YIELD FUNCTION:f,VECTOR:{df/dp}.*/
  {       /*f=POW(|Qx/Qxu|,EXP)+POW(|Qy/Qyu|,EXP)+POW(|Nz/Nzu|,EXP)*/
		  /* +POW(|Mx/Mxu|,EXP)+POW(|My/Myu|,EXP)+POW(|Mz/Mzu|,EXP)*/
	f[i]=0.0;
	for(j=0;j<6;j++)
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

/*
if(globalfile!=NULL)
fprintf(globalfile,"%d %d %10.5f %10.5f %10.5f\n",
		i,j,elem.stress[i][j],value,dfdp[i][j]);
*/

      }
	}
  }

  for(i=0;i<12;i++)                         /*ASSEMBLAGE VECTOR:{q}*/
  {
	if(iconf[i]!=1)
    {
      for(j=0;j<2;j++)
      {
        q[i][j]=0.0;

/*
if(globalfile!=NULL)
fprintf(globalfile,"q%d%d =\n",i,j);
*/

		for(jj=0;jj<6;jj++)
		{
		  if(iconf[6*j+jj]!=1) /*if(iconf[6*j+jj]==-1)*/
		  {
			q[i][j]+=*(*(estiff+i)+6*j+jj)*dfdp[j][jj];

/*
if(globalfile!=NULL)
fprintf(globalfile,"    + %10.5f x %10.5f\n",
        *(*(estiff+i)+6*j+jj),dfdp[j][jj]);
*/

		  }
        }

/*
if(globalfile!=NULL)
fprintf(globalfile,"    = %10.5f\n",q[i][j]);
*/

	  }
	}
  }
  for(i=0;i<2;i++)                     /*INNER PRODUCT a={df/dp}{q}*/
  {
	for(j=0;j<2;j++)
	{
	  a[i][j]=0.0;

/*
if(globalfile!=NULL)
fprintf(globalfile,"a%d%d =\n",i,j);
*/

	  for(ii=0;ii<6;ii++)
	  {
        if((elem.iconf[i][ii]==-1)&&(elem.iconf[j][ii]==-1))
		{
		  a[i][j]+=dfdp[i][ii]*q[6*i+ii][j];

/*
if(globalfile!=NULL)
fprintf(globalfile,"    + %10.5f x %10.5f\n",
		dfdp[i][ii],q[6*i+ii][j]);
*/

		}
	  }

/*
if(globalfile!=NULL)
fprintf(globalfile,"    = %10.5f\n",a[i][j]);
*/

	}
  }
  return;
}/*coefficients*/
#endif

void coefficientsbc(struct owire elem,double **estiff,
                  double f[],double dfdp[][6],
                  double q[][2],double a[][2],double ncr)
/*ASSEMBLAGE PLASTIC COEFFICIENTS.*/
/*UJIOKA:REFLECTING BUCKLING CONDENSATION RESULT*/
{
  char string[256];
  signed char iconf[12];
  int i,j,ii,jj;
  double fc[6],fu[6],unit,value;

  for(i=0;i<2;i++)                    /*ICONF[2][6] INTO ICONF[12].*/
  {
    for(j=0;j<6;j++) iconf[6*i+j]=elem.iconf[i][j];
  }

/*
if(globalfile!=NULL) fprintf(globalfile,"\nELEM %d\n",elem.code);
if(globalfile!=NULL) fprintf(globalfile,"   CENTER fc   WIDTH fu\n");
*/

  for(i=0;i<6;i++)                 /*CENTER,WIDTH OF YIELD SURFACE.*/
  {
	fc[i]=0.5*(elem.sect->fmax[i]+elem.sect->fmin[i]);     /*CENTER*/
	fu[i]=0.5*(elem.sect->fmax[i]-elem.sect->fmin[i]);      /*WIDTH*/

/*if(globalfile!=NULL) fprintf(globalfile,"%d %10.5f %10.5f\n",
                             i,fc[i],fu[i]);*/
  }

  /***UJIOKA:圧縮部材について、降伏曲面を縮約による座屈荷重(Ncr)に修正***/
  /*
引張部材について、ここでは便宜上Ncr=0になっており、降伏曲面の修正をしない（圧縮側と引張側で降伏曲面の式を場合分けしていないため）
  */
  if(ncr>0 && ncr<elem.sect->fmax[0])
  {
    fc[0]=0;
    fu[0]=ncr;       /*圧縮・引張の両側を修正している*/
//    sprintf(string,"fu[0]=%.5E\n",fu[0]);
//    errormessage(string);

  /*
    fc[0]=0.5*(ncr+elem->sect->fmin[i]);
    fu[0]=0.5*(ncr-elem->sect->fmin[i]);
  */
  }
  /*
　　圧縮部材の座屈荷重が負の値となる場合：未検討、元の降伏曲面のまま
  */
  else if (ncr<0) /*negative buckling load in compression element:yet*/
  {
//     sprintf(string,"Negative Buckling Load in Compression Element\n");
//     errormessage(string);
  }
  /***UJIOKA***/

/*
if(globalfile!=NULL)
fprintf(globalfile,"        STRESS       RATE       dfdp\n");
*/

  for(i=0;i<2;i++)    /*ASSEMBLAGE YIELD FUNCTION:f,VECTOR:{df/dp}.*/
  {       /*f=POW(|Qx/Qxu|,EXP)+POW(|Qy/Qyu|,EXP)+POW(|Nz/Nzu|,EXP)*/
          /* +POW(|Mx/Mxu|,EXP)+POW(|My/Myu|,EXP)+POW(|Mz/Mzu|,EXP)*/
    f[i]=0.0;
	for(j=0;j<6;j++)
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

/*
if(globalfile!=NULL)
fprintf(globalfile,"%d %d %10.5f %10.5f %10.5f\n",
        i,j,elem.stress[i][j],value,dfdp[i][j]);
*/

      }
    }
  }

  for(i=0;i<12;i++)                         /*ASSEMBLAGE VECTOR:{q}*/
  {
    if(iconf[i]!=1)
    {
	  for(j=0;j<2;j++)
	  {
		q[i][j]=0.0;

/*
if(globalfile!=NULL)
fprintf(globalfile,"q%d%d =\n",i,j);
*/

		for(jj=0;jj<6;jj++)
		{
		  if(iconf[6*j+jj]!=1) /*if(iconf[6*j+jj]==-1)*/
		  {
            q[i][j]+=*(*(estiff+i)+6*j+jj)*dfdp[j][jj];

/*
if(globalfile!=NULL)
fprintf(globalfile,"    + %10.5f x %10.5f\n",
        *(*(estiff+i)+6*j+jj),dfdp[j][jj]);
*/

		  }
		}

/*
if(globalfile!=NULL)
fprintf(globalfile,"    = %10.5f\n",q[i][j]);
*/

	  }
	}
  }
  for(i=0;i<2;i++)                     /*INNER PRODUCT a={df/dp}{q}*/
  {
	for(j=0;j<2;j++)
	{
	  a[i][j]=0.0;

/*
if(globalfile!=NULL)
fprintf(globalfile,"a%d%d =\n",i,j);
*/

      for(ii=0;ii<6;ii++)
      {
        if((elem.iconf[i][ii]==-1)&&(elem.iconf[j][ii]==-1))
        {
          a[i][j]+=dfdp[i][ii]*q[6*i+ii][j];

/*
if(globalfile!=NULL)
fprintf(globalfile,"    + %10.5f x %10.5f\n",
        dfdp[i][ii],q[6*i+ii][j]);
*/

		}
      }

/*
if(globalfile!=NULL)
fprintf(globalfile,"    = %10.5f\n",a[i][j]);
*/

    }
  }
  return;
}/*coefficients*/

double **modifyhinge(struct owire elem,double **estiff)
/*MODIFY ELEMENT STIFFNESS MATRX BY HINGE.*/
{
  char s[80];
  int n,i,ii,jj,kk;
  double h[12][12],**e=estiff;


  /*BOTH ENDS CAN NOT BE HINGE FOR N,Mz*/
  if(elem.iconf[0][0]== 1 && elem.iconf[1][0]== 1) elem.iconf[1][0]=0;
  if(elem.iconf[0][0]== 1 && elem.iconf[1][0]==-2) elem.iconf[1][0]=0;
  if(elem.iconf[0][0]== 1 && elem.iconf[1][0]==-3) elem.iconf[1][0]=0;
  if(elem.iconf[0][0]==-2 && elem.iconf[1][0]==-2) elem.iconf[1][0]=0;
  if(elem.iconf[0][0]==-2 && elem.iconf[1][0]==-3) elem.iconf[1][0]=0;
  if(elem.iconf[0][0]==-3 && elem.iconf[1][0]==-3) elem.iconf[1][0]=0;
  if(elem.iconf[0][3]== 1 && elem.iconf[1][3]== 1) elem.iconf[1][3]=0;

  for(n=0;n<=1;n++)
  {
	for(i=0;i<6;i++)
	{
	  if(elem.iconf[n][i]== 1 ||
		 elem.iconf[n][i]==-2 ||
		 elem.iconf[n][i]==-3)
	  {
		kk=6*n+i;

if(*(*(e+kk)+kk)!=0.0) /*UNDER CONSIDERATION.*/
{
#if 0
		if(*(*(e+kk)+kk)==0.0) /*UNDER CONSIDERATION.*/
		{
		  sprintf(s,"Modify : Matrix Singular.ELEM=%ld SECT=%d",elem.code,elem.sect->code);
		  MessageBox(NULL,s,"Modifyhinge",MB_OK);
		}
#endif

		for(ii=0;ii<12;ii++)
		{
		  for(jj=0;jj<12;jj++)
		  {
			h[ii][jj]=- *(*(e+ii)+kk)
					  / *(*(e+kk)+kk)
					  * *(*(e+kk)+jj);
		  }
		}
		for(ii=0;ii<12;ii++)
		{
          for(jj=0;jj<12;jj++) *(*(e+ii)+jj)+=h[ii][jj];
        }
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

void transformationII(double **estiff,double **tmatrix,
					  double **e,double **t)
/*TRANSFORM ELEMENT MATRIX INTO GLOBAL WITHOUT MALLOC.*/
{
  matrixmatrixII(e,estiff,tmatrix,12);
  matrixtransposeII(t,tmatrix,12);                           /*[Tt]*/
  matrixmatrixII(estiff,t,e,12);                   /*[K]=[Tt][k][T]*/
  return;
}/*transformationII*/

double **transformationIII(double **estiff,double **tmatrix,int msize)
/*TRANSFORM ELEMENT MATRIX INTO GLOBAL.*/
{
  int i;
  double **e,**t;

  e=matrixmatrix(estiff,tmatrix,msize);
  t=matrixtranspose(tmatrix,msize);                             /*[Tt]*/

  for(i=0;i<msize;i++) free(*(estiff+i));
  free(estiff);

  estiff=matrixmatrix(t,e,msize);                     /*[K]=[Tt][k][T]*/

  for(i=0;i<msize;i++) free(*(e+i));
  free(e);
  for(i=0;i<msize;i++) free(*(t+i));
  free(t);

  return estiff;
}/*transformationIII*/

double **transformationEx(double **A,double **T,int row, int col)
/*TRANSFORM ELEMENT MATRIX INTO GLOBAL.*/
{
  int i;
  double **AT,**Tt,**TtAT;

  AT=matrixmatrixIII(A,T,row,row,col);
  Tt=matrixtransposeIII(T,row,col);
  TtAT=matrixmatrixIII(Tt,AT,col,row,col);

  freematrix(AT,row);
  freematrix(Tt,col);

  return TtAT;
}/*transformationEx*/

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
			} /*THESE ARE PUT TOGETHER INTO "gadd".UNDER CONSTRUCTION.*/
		  }
		}
	  }
	}
  }
  return;
}/*assemgstiffness*/

void assemgstiffnessII(struct gcomponent *gmtx,
					   double **estiff,
					   struct oshell *shell)
/*ASSEMBLAGE ELEMENT MATRIX INTO GLOBAL MATRIX.*/
{
  long int i,j,ii,jj,n1,n2;
  double edata,gdata;

  for(n1=0;n1<3;n1++)
  {
	for(i=1;i<=6;i++)
	{
	  ii=6*(shell->node[n1]->loff)+i;
	  for(n2=0;n2<3;n2++)
	  {
		for(j=1;j<=6;j++)
		{
		  jj=6*(shell->node[n2]->loff)+j;
		  if(ii>=jj)
		  {
			edata=*(*(estiff+6*n1+i-1)+6*n2+j-1);
			if(edata!=0.0)
			{
			  gread(gmtx,ii,jj,&gdata);
			  gdata+=edata;
			  gwrite(gmtx,ii,jj,gdata);
			} /*THESE ARE PUT TOGETHER INTO "gadd".UNDER CONSTRUCTION.*/
		  }
		}
	  }
	}
  }
  return;
}/*assemgstiffnessII*/

void assemgstiffnesswithDOFelimination(struct gcomponent *gmtx,
					   double **estiff,
					   struct owire *elem,
					   long int *constraintmain)
/*ASSEMBLAGE ELEMENT MATRIX INTO GLOBAL MATRIX.*/
{
  long int i,j,ii,jj,n1,n2;
  double edata,gdata;

  for(n1=0;n1<2;n1++)
  {
	for(i=1;i<=6;i++)
	{
	  ii=6*(elem->node[n1]->loff)+i;
	  ii=*(constraintmain+ii-1)+1;
	  for(n2=0;n2<2;n2++)
	  {
		for(j=1;j<=6;j++)
		{
		  jj=6*(elem->node[n2]->loff)+j;
		  jj=*(constraintmain+jj-1)+1;
		  if(ii>=jj)
		  {
			edata=*(*(estiff+6*n1+i-1)+6*n2+j-1);
			if(edata!=0.0)
			{
			  gread(gmtx,ii,jj,&gdata);
			  gdata+=edata;
			  gwrite(gmtx,ii,jj,gdata);
			} /*THESE ARE PUT TOGETHER INTO "gadd".UNDER CONSTRUCTION.*/
		  }
		}
	  }
	}
  }
  return;
}/*assemgstiffnesswithDOFelimination*/

void assemgstiffnessIIwithDOFelimination(struct gcomponent *gmtx,
					   double **estiff,
					   struct oshell *shell,
					   long int *constraintmain)
/*ASSEMBLAGE ELEMENT MATRIX INTO GLOBAL MATRIX.*/
{
  long int i,j,ii,jj,n1,n2;
  double edata,gdata;

  for(n1=0;n1<3;n1++)
  {
	for(i=1;i<=6;i++)
	{
	  ii=6*(shell->node[n1]->loff)+i;
	  ii=*(constraintmain+ii-1)+1;
	  for(n2=0;n2<3;n2++)
	  {
		for(j=1;j<=6;j++)
		{
		  jj=6*(shell->node[n2]->loff)+j;
		  jj=*(constraintmain+jj-1)+1;
		  if(ii>=jj)
		  {
			edata=*(*(estiff+6*n1+i-1)+6*n2+j-1);
			if(edata!=0.0)
			{
			  gread(gmtx,ii,jj,&gdata);
			  gdata+=edata;
			  gwrite(gmtx,ii,jj,gdata);
			} /*THESE ARE PUT TOGETHER INTO "gadd".UNDER CONSTRUCTION.*/
		  }
		}
	  }
	}
  }
  return;
}/*assemgstiffnessIIwithDOFelimination*/

void assemconf(struct oconf *confs,double *gvct,double dsafety,
			   int nnode)
/*ASSEMBLAGE CONFINEMENT VALUE INTO GLOBAL VECTOR.*/
{
  int i;
  double conf;
  signed char iconf;

  for(i=0;i<6*nnode;i++)      /*LOADS OR DISPS GIVEN INCREMENTALLY.*/
  {
	conf=(confs+i)->value;
	iconf=(confs+i)->iconf;

	//*(gvct+i-1)+=dsafety*conf;   /*"+=":IF GVECTOR INITIALIZED.*/
	//*(gvct+i)=dsafety*conf;         /*"=":IF GVECTOR UNINITIALIZED.*/

	if(iconf!=1 && conf!=0.0)
	{
	   *(gvct+i)=dsafety*conf;         /*"=":IF GVECTOR UNINITIALIZED.*/
	}
	else
	{
	   *(gvct+i)=0.0;
	}

  }

  return;
}/*assemconf*/

void assemgivend(struct oconf *confs,double *gvct,double dsafety,
			   int nnode)
{
  int i;
  double conf;
  signed char iconf;


  for(i=0;i<6*nnode;i++)
  {
	conf=(confs+i)->value;
	iconf=(confs+i)->iconf;

	if(iconf==1 && conf!=0.0)
	{
	   *(gvct+i)=dsafety*conf;         /*"=":IF GVECTOR UNINITIALIZED.*/
	}
	else
	{
	   *(gvct+i)=0.0;
	}
  }
  return;
}


void assemnodenorm(struct oshell shell,double *gvct)
/*ASSEMBLAGE CONFINEMENT VALUE INTO GLOBAL VECTOR.*/
{
  int i,j;
  long int loff;
  struct onode n1,n2,n3,n4;
  double *nvct;
  nvct=(double *)malloc(3*sizeof(double));
  n1=*(shell.node[0]);
  n2=*(shell.node[1]);
  n3=*(shell.node[2]);
  *(nvct+0)=(n2.d[GY]-n1.d[GY])*(n3.d[GZ]-n1.d[GZ])-(n2.d[GZ]-n1.d[GZ])*(n3.d[GY]-n1.d[GY]);
  *(nvct+1)=(n2.d[GZ]-n1.d[GZ])*(n3.d[GX]-n1.d[GX])-(n2.d[GX]-n1.d[GX])*(n3.d[GZ]-n1.d[GZ]);
  *(nvct+2)=(n2.d[GX]-n1.d[GX])*(n3.d[GY]-n1.d[GY])-(n2.d[GY]-n1.d[GY])*(n3.d[GX]-n1.d[GX]);

  if(shell.nnod==4)
  {
	n4=*(shell.node[3]);
	*(nvct+0)+=(n4.d[GY]-n3.d[GY])*(n1.d[GZ]-n3.d[GZ])-(n4.d[GZ]-n3.d[GZ])*(n1.d[GY]-n3.d[GY]);
	*(nvct+1)+=(n4.d[GZ]-n3.d[GZ])*(n1.d[GX]-n3.d[GX])-(n4.d[GX]-n3.d[GX])*(n1.d[GZ]-n3.d[GZ]);
	*(nvct+2)+=(n4.d[GX]-n3.d[GX])*(n1.d[GY]-n3.d[GY])-(n4.d[GY]-n3.d[GY])*(n1.d[GX]-n3.d[GX]);
  }


  for(i=0;i<shell.nnod;i++)      /*LOADS OR DISPS GIVEN INCREMENTALLY.*/
  {
	loff=shell.node[i]->loff;

	for(j=0;j<3;j++)
	{
	  *(gvct+6*loff+j)+=0.5**(nvct+j)/shell.nnod;
	}
  }
  free(nvct);
  return;
}/*assemnodenorm*/


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

void modifycmq201(struct memoryelem *melem,struct owire *elem)
/*MODIFY CMQ BY HINGE.*/
/*NOT SAME AS "modifycmq".*/
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
  if((elem->iconf[0][5]==0)&&(elem->iconf[1][5]==1))
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
/*if(j>=3) elem->stress[i][j]=0.0;*/ /*NONLINEAR TEST*/

	  (melem+(elem->loff))->stress[i][j]=elem->stress[i][j];
    }
  }

  return;
}/*modifycmq201*/

void assemcmq201(struct owire elem,double **tmatrix,
			  struct oconf *confs,double *gvct)
/*ASSEMBLAGE CMQ INTO GLOBAL VECTOR.*/
/*SAME AS "assemcmq".*/
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
}/*assemcmq201*/

void assemconf201(struct oconf *confs,double *gvct,
				  double dsafety,int nnode)
/*ASSEMBLAGE CONFINEMENT VALUE INTO GLOBAL VECTOR.*/
/*SAME AS "assemconf001".*/
{
  int i;
  double conf;

  for(i=1;i<=6*nnode;i++)     /*LOADS OR DISPS GIVEN INCREMENTALLY.*/
  {
	conf=(confs+i-1)->value;

	*(gvct+i-1)+=dsafety*conf;                 /*"+=":ADD WITH CMQ.*/
  }

  return;
}/*assemconf201*/

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
  for(n=0;n<2;n++)
  {
	for(i=0;i<6;i++)
	{
	  loffset=6*(elem.node[n]->loff)+i;
	  *(d+6*n+i)=*(gvct+loffset);
	}
  }
  return d;
}/*extractdisplacement*/

double *extractdisplacement2(struct owire elem,double *gvct,
							 long int *moff)
/*EXTRACT ELEMENT DEFORMATION{dU} FROM GLOBAL VECTOR.*/
{
  long int i,m,n;
  double *d;

  d=(double *)malloc(12*sizeof(double));
  if(d==NULL)
  {
	errormessage("EXTRACTDISPLACEMENT:MEMORY INSUFFICIENT.");
	return NULL;
  }
  for(n=0;n<2;n++)
  {
	if(moff==NULL) m=elem.node[n]->loff;
	else           m=*(moff+(elem.node[n]->loff));

	for(i=0;i<6;i++) *(d+6*n+i)=*(gvct+6*m+i);
  }
  return d;
}/*extractdisplacement2*/

void extractdisplacementII(double *d,
						   struct owire elem,double *gvct)
/*EXTRACT ELEMENT DEFORMATION{dU} FROM GLOBAL VECTOR.*/
{
  long int i,loffset;
  int n;

  for(n=0;n<2;n++)
  {
	for(i=0;i<6;i++)
	{
	  loffset=6*(elem.node[n]->loff)+i;
	  *(d+6*n+i)=*(gvct+loffset);
	}
  }
  return;
}/*extractdisplacementII*/


double *elemstress(struct owire *elem,double *gvct,
				   struct memoryelem *melem,FILE *fout,
                   double func[])
/*ELEMENT STRESS INCREMENTAL.*/
{
  int i,j;
  double **drccos,**tmatrix,**estiff,**ee,*gdisp,*edisp,*estress;

  ee=mallocdoublematrix(12);

  drccos=directioncosine(elem->node[0]->d[0],
						 elem->node[0]->d[1],
						 elem->node[0]->d[2],
						 elem->node[1]->d[0],
						 elem->node[1]->d[1],
						 elem->node[1]->d[2],
						 elem->cangle);                  /*[DRCCOS]*/
  tmatrix=transmatrix(drccos);                     /*[T]*/
  estiff=assememtx(*elem);                                   /*[ke]*/
  estiff=modifyhinge(*elem,estiff);

  for(i=0;i<12;i++)
  {
	for(j=0;j<12;j++) *(*(ee+i)+j)=*(*(estiff+i)+j);
  }

  estiff=assempmtx(*elem,estiff);                   /*[k]=[ke]+[kp]*/

  gdisp=extractdisplacement(*elem,gvct);                     /*{dU}*/
  edisp=matrixvector(tmatrix,gdisp,12);              /*{du}=[T]{dU}*/
  estress=matrixvector(estiff,edisp,12);             /*{df}=[k]{du}*/

  updatestress(melem,fout,edisp,estress,ee,elem,
			   func,NULL);                               /*{f}+{df}*/

  free(gdisp);
  free(edisp);
  freematrix(drccos,3);
  freematrix(tmatrix,12);
  freematrix(estiff,12);
  freematrix(ee,12);

  return estress;
}/*elemstress*/


double *elemstressbc(struct owire *elem,double *gvct,
                   struct memoryelem *melem,FILE *fout,FILE *fsrf,
                   double func[],double ncr)
/*ELEMENT STRESS INCREMENTAL.*/
/*UJIOKA:REFLECTING BUCKLING CONDENSATION RESULT*/
{
  int i,j;
  double **drccos,**tmatrix,**estiff,**ee,*gdisp,*edisp,*estress;
  FILE  *ftxt;

  ee=mallocdoublematrix(12);

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

  for(i=0;i<12;i++)
  {
    for(j=0;j<12;j++) *(*(ee+i)+j)=*(*(estiff+i)+j);
  }

  estiff=assempmtxbc(*elem,estiff,ncr);                   /*[k]=[ke]+[kp]*/

  gdisp=extractdisplacement(*elem,gvct);                     /*{dU}*/
  edisp=matrixvector(tmatrix,gdisp,12);              /*{du}=[T]{dU}*/
  estress=matrixvector(estiff,edisp,12);             /*{df}=[k]{du}*/

//  ftxt=fopen("energy.txt","w");

  updatestressbc(melem,fout,fsrf,edisp,estress,ee,elem,
               func,ftxt/*NULL*/,ncr);                               /*{f}+{df}*/

  free(gdisp);
  free(edisp);
  freematrix(drccos,3);
  freematrix(tmatrix,12);
  freematrix(estiff,12);
  freematrix(ee,12);

  return estress;
}/*elemstressbc*/

void elemstressII(double *estress,
                  struct owire *elem,
                  double *gvct,struct memoryelem *melem,FILE *fout,
                  double **drccos,
                  double **tmatrix,
                  double **estiff,
                  double *gdisp,double *edisp,
                  double func[],FILE *ftxt)
/*ELEMENT STRESS INCREMENTAL WITHOUT MALLOC.*/
{
  int i,j;
  double **ee;

  ee=mallocdoublematrix(12);

  directioncosineII(elem->node[0]->d[0],
                    elem->node[0]->d[1],
                    elem->node[0]->d[2],
                    elem->node[1]->d[0],
                    elem->node[1]->d[1],
                    elem->node[1]->d[2],
                    elem->cangle,
                    drccos);                             /*[DRCCOS]*/
  transmatrixII(drccos,tmatrix);                              /*[T]*/
  assememtxII(*elem,estiff);                                 /*[ke]*/
  estiff=modifyhinge(*elem,estiff);

  for(i=0;i<12;i++)
  {
    for(j=0;j<12;j++) *(*(ee+i)+j)=*(*(estiff+i)+j);
  }

  estiff=assempmtx(*elem,estiff);                   /*[k]=[ke]+[kp]*/

  extractdisplacementII(gdisp,*elem,gvct);                   /*{dU}*/
  matrixvectorII(edisp,tmatrix,gdisp,12);            /*{du}=[T]{dU}*/
  matrixvectorII(estress,estiff,edisp,12);           /*{df}=[k]{du}*/


  updatestress(melem,fout,edisp,estress,ee,elem,
			   func,ftxt);                               /*{f}+{df}*/

  freematrix(ee,12);

  return;
}/*elemstressII*/

void elemstressIIbc(double *estress,
                  struct owire *elem,
                  double *gvct,struct memoryelem *melem,FILE *fout,
                  double **drccos,
                  double **tmatrix,
                  double **estiff,
                  double *gdisp,double *edisp,
                  double func[],FILE *ftxt,
                  double ncr)
/*ELEMENT STRESS INCREMENTAL WITHOUT MALLOC.*/
{
  int i,j;
  double **ee;

  ee=mallocdoublematrix(12);

  directioncosineII(elem->node[0]->d[0],
                    elem->node[0]->d[1],
                    elem->node[0]->d[2],
                    elem->node[1]->d[0],
                    elem->node[1]->d[1],
                    elem->node[1]->d[2],
                    elem->cangle,
                    drccos);                             /*[DRCCOS]*/
  transmatrixII(drccos,tmatrix);                              /*[T]*/
  assememtxII(*elem,estiff);                                 /*[ke]*/
  estiff=modifyhinge(*elem,estiff);

  for(i=0;i<12;i++)
  {
    for(j=0;j<12;j++) *(*(ee+i)+j)=*(*(estiff+i)+j);
  }

//  estiff=assempmtx(*elem,estiff);                   /*[k]=[ke]+[kp]*/
  estiff=assempmtxbc(*elem,estiff,ncr);                   /*[k]=[ke]+[kp]*/

  extractdisplacementII(gdisp,*elem,gvct);                   /*{dU}*/
  matrixvectorII(edisp,tmatrix,gdisp,12);            /*{du}=[T]{dU}*/
  matrixvectorII(estress,estiff,edisp,12);           /*{df}=[k]{du}*/

  /*fprintf(fout,"ELEM %d [Kp]{du}\n",elem->code);
  for(i=0;i<12;i++)
  {
    data=0.0;
    for(j=0;j<12;j++)
    {
      data+=(*(*(ee+i)+j)-*(*(estiff+i)+j))*(*(edisp+j));
    }
    fprintf(fout,"%8.5f\n",data);
  }*/

/*  updatestress(melem,fout,edisp,estress,ee,elem,
               func,ftxt);   */                            /*{f}+{df}*/
  updatestressbc(melem,fout,NULL,edisp,estress,ee,elem,
               func,ftxt,ncr);                               /*{f}+{df}*/

  freematrix(ee,12);

  return;
}/*elemstressIIbc*/

double *elemstressnl(struct owire *elem,double *gvct,
				   struct memoryelem *melem)
/*ELEMENT STRESS INCREMENTAL FOR NON-LINEAR ANALYSIS.*/
{
  int i,j;
  double **drccos,**tmatrix,**estiff,*gdisp,*edisp,*estress;

  drccos=directioncosine(elem->node[0]->d[0],
						 elem->node[0]->d[1],
						 elem->node[0]->d[2],
						 elem->node[1]->d[0],
						 elem->node[1]->d[1],
						 elem->node[1]->d[2],
						 elem->cangle);                  /*[DRCCOS]*/
  tmatrix=transmatrix(drccos);                                /*[T]*/
  estiff=assememtx(*elem);                                   /*[ke]*/
  estiff=modifyhinge(*elem,estiff);

  gdisp=extractdisplacement(*elem,gvct);                     /*{dU}*/
  edisp=matrixvector(tmatrix,gdisp,12);              /*{du}=[T]{dU}*/
  estress=matrixvector(estiff,edisp,12);             /*{df}=[k]{du}*/

  updatestressnl(melem,estress,elem);                    /*{f}+{df}*/

  free(gdisp);
  free(edisp);
  freematrix(drccos,3);
  freematrix(tmatrix,12);
  freematrix(estiff,12);

  return estress;
}/*elemstressnl*/

void updatestress(struct memoryelem *melem,FILE *fout,
				  double *edisp,double *dstress,double **estiff,
				  struct owire *elem,double func[],FILE *ftxt)
/*ELEMENT STRESS UPDATE.*/
{
  char s[80],string[256];
  char iconf[12];
  long int i,ii,j/*,jj*/,nn[2];
  double dL,lamda[2]={0.0,0.0};
  double fc[6],fu[6],f[2],dfdp[2][6],q[12][2],a[2][2],rate;
  double fe[2],dfdpe[2][6],qe[12][2],ae[2][2];     /*1 STEP BEFORE.*/
  double det,detinverse,function;
  double ys[2][6];
  double due[2][6],dup[2][6]/*,dp[2][6]*/;
  struct line line;

  nn[0]=elem->node[0]->code;
  nn[1]=elem->node[1]->code;

  /*FUNCTION 1 STEP BEFORE.*/
  coefficients(*elem,estiff,fe,dfdpe,qe,ae);

  /*UPDATE STRESS.*/
  for(i=0;i<2;i++)
  {
    for(j=0;j<6;j++)
	{
      elem->stress[i][j]+=*(dstress+6*i+j);
	  (melem+(elem->loff))->stress[i][j]=elem->stress[i][j];
    }
  }

  coefficients(*elem,estiff,f,dfdp,q,a);        /*UPDATED FUNCTION.*/

  for(i=0;i<2;i++)                    /*ICONF[2][6] INTO ICONF[12].*/
  {
	for(j=0;j<6;j++) iconf[6*i+j]=elem->iconf[i][j];
  }

  for(i=0;i<6;i++)                 /*CENTER,WIDTH OF YIELD SURFACE.*/
  {
    fc[i]=0.5*(elem->sect->fmax[i]+elem->sect->fmin[i]);
    fu[i]=0.5*(elem->sect->fmax[i]-elem->sect->fmin[i]);
  }

  for(i=0;i<2;i++)                                       /*INITIAL.*/
  {
    for(j=0;j<6;j++) dup[i][j]=0.0;
  }

  /*LAMDA OF PLASTIC END.*/
  if(ae[0][0]!=0.0 && ae[1][1]==0.0)       /*IF I:PLASTIC J:ELASTIC*/
  {
    for(i=0;i<12;i++)
    {
      if(iconf[i]!=1 && iconf[i]!=-2 && iconf[i]!=-3)
      {
        lamda[0]+=1.0/ae[0][0]*qe[i][0]*(*(edisp+i));
      }
    }
    for(j=0;j<6;j++)
    {
      dup[0][j]=lamda[0]*dfdpe[0][j];
    }
  }
  else if(ae[0][0]==0.0 && ae[1][1]!=0.0)  /*IF I:ELASTIC J:PLASTIC*/
  {
    for(i=0;i<12;i++)
    {
      if(iconf[i]!=1 && iconf[i]!=-2 && iconf[i]!=-3)
      {
        lamda[1]+=1.0/ae[1][1]*qe[i][1]*(*(edisp+i));
      }
    }
    for(j=0;j<6;j++)
    {
      dup[1][j]=lamda[1]*dfdpe[1][j];
	}
  }
  else if(ae[0][0]!=0.0 && ae[1][1]!=0.0)  /*IF I:PLASTIC J:PLASTIC*/
  {
    for(i=0;i<12;i++)
    {
      if(iconf[i]!=1 && iconf[i]!=-2 && iconf[i]!=-3)
      {
        detinverse=ae[0][0]*ae[1][1]-ae[0][1]*ae[1][0];
        if(detinverse==0.0)
        {
          det=0.0;                           /*UNDER CONSIDERATION.*/
          errormessage("UPDATESTRESS:UNDER CONSIDERATION.");
		  sprintf(string,"Update : Matrix Singular.SECT=%d",elem->sect->code);
        }
		else det=1.0/detinverse;

        lamda[0]+=det*( ae[1][1]*qe[i][0]*(*(edisp+i))
                       -ae[0][1]*qe[i][1]*(*(edisp+i)));
        lamda[1]+=det*(-ae[1][0]*qe[i][0]*(*(edisp+i))
                       +ae[0][0]*qe[i][1]*(*(edisp+i)));
      }
    }
    for(j=0;j<6;j++)
    {
      dup[0][j]=lamda[0]*dfdpe[0][j];
      dup[1][j]=lamda[1]*dfdpe[1][j];
    }
  }

  /*ELASTIC END.*/
  if(ae[0][0]==0.0 || ae[1][1]==0.0)    /*IF I:ELASTIC OR J:ELASTIC*/
  {
    if(elem->sect->Ixx==0.0 &&
       elem->sect->Iyy==0.0 &&
       elem->sect->Jzz==0.0)                            /*FOR BRACE*/
	{
      dL=*(edisp+6)-*(edisp+0); /*EXTENSION OF LENGTH*/

      if(dL< 0.0 && elem->iconf[0][0]==-2 && ae[0][0]==0.0)
      {
        lamda[0]=-1.0;                         /*LONGING DISLOADED.*/
		if(fout!=NULL) fprintf(fout,"ELEM%d HEAD dL=%.5E LONGING DISLOADED.\n",elem->code,dL);
      }
	  if(dL< 0.0 && elem->iconf[1][0]==-2 && ae[1][1]==0.0)
      {
		lamda[1]=-1.0;                         /*LONGING DISLOADED.*/
		if(fout!=NULL) fprintf(fout,"ELEM%d TAIL dL=%.5E LONGING DISLOADED.\n",elem->code,dL);
	  }
	  if(dL>=0.0 && elem->iconf[0][0]==-3 && ae[0][0]==0.0)
      {
        if(ae[0][0]==0.0) lamda[0]=-1.0;      /*SHORTING DISLOADED.*/
		if(fout!=NULL) fprintf(fout,"ELEM%d HEAD dL=%.5E SHORTING DISLOADED.\n",elem->code,dL);
      }
      if(dL>=0.0 && elem->iconf[1][0]==-3 && ae[1][1]==0.0)
      {
        if(ae[1][1]==0.0) lamda[1]=-1.0;      /*SHORTING DISLOADED.*/
		if(fout!=NULL) fprintf(fout,"ELEM%d TAIL dL=%.5E SHORTING DISLOADED.\n",elem->code,dL);
      }
	}
	else if(ae[0][0]==0.0 && fe[0]>=f[0])
    {
      lamda[0]=-1.0;                                 /*I:DISLOADED.*/
	}
    else if(ae[1][1]==0.0 && fe[1]>=f[1])
    {
      lamda[1]=-1.0;                                 /*J:DISLOADED.*/
    }
  }

  /*dUe,dUp DECISION.*/
  for(i=0;i<2;i++)
  {
    for(j=0;j<6;j++)
    {
      due[i][j]=*(edisp+6*i+j)-dup[i][j];        /*{due}={du}-{dup}*/
	}
  }

  /*STRAIN ENERGY.*/
  for(j=0;j<6;j++)
  {
    for(i=0;i<2;i++)
    {
      elem->Ee[0]+=0.25*(2*elem->stress[i][j]-(*(dstress+6*i+j)))
					   *due[i][j];
      elem->Ee[1]+=0.25*(2*elem->stress[i][j]-(*(dstress+6*i+j)))
                       *due[i][j];
    }
    elem->Ep[0]+=(elem->stress[0][j])*dup[0][j];
    elem->Ep[1]+=(elem->stress[1][j])*dup[1][j];
  }


  /*YIELD JUDGEMENT.*/
  for(i=0;i<2;i++)
  {
    func[i]=f[i];

	if(f[i]>=pow(RADIUS,EXPONENT))                       /*YIELDED.*/
    {
      sprintf(string,"YIELDED:ELEM%d NODE%ld SECT%d",
              elem->code,nn[i],elem->sect->code);
      /*errormessage(string);*/
      strcat(string," FUNCTION:{");
      function=0.0;
      for(ii=0;ii<6;ii++)
      {
        if(ii==0 || ii==3) /*N,Mz*/
        {
		  rate=(elem->stress[i][ii]-pow(-1.0,i*1.0)*fc[ii])/fu[ii];
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
          if(elem->sect->Ixx==0.0 &&
             elem->sect->Iyy==0.0 &&
             elem->sect->Jzz==0.0) /*FOR BRACE*/
          {
            /*elem->sect->area=0.0;*/    /*WRONG.ALL ELEM WITH THIS*/
										 /*SECT WILL BE SET 0.     */

            dL=*(edisp+6)-*(edisp+0); /*EXTENSION OF LENGTH*/

            if(ii==0)
            {
			  if(dL>=0.0 && ((i==0 && rate<0.0) || (i==1 && rate>0.0)))
              {
                elem->iconf[i][0]=-2;
				if(fout!=NULL) fprintf(fout,"ELEM%d dL=%.5E LONGING.\n",elem->code,dL);
              }
			  else if(dL< 0.0 && ((i==0 && rate>0.0) || (i==1 && rate<0.0)))
              {
				elem->iconf[i][0]=-3;
				if(fout!=NULL) fprintf(fout,"ELEM%d dL=%.5E SHORTING.\n",elem->code,dL);
			  }
            }
		  }
          else
          {
            elem->iconf[i][ii]=-1;   /*-1:PLASTIC 0:ELASTIC 1:HINGE*/
          }
		}
	  }
	  sprintf(s,"}=%8.5f",function); /*VALUE OF FUNCTION.*/
	  strcat(string,s);
	  if(fout!=NULL) fprintf(fout,"%s\n",string);
    }
  }

  /*DISLOAD JUDGEMENT.*/
  for(i=0;i<2;i++)
  {
    if(lamda[i]<0.0)                                   /*DISLOADED.*/
    {
      sprintf(string,"DISLOAD:ELEM%d NODE%ld SECT%d",
              elem->code,nn[i],elem->sect->code);

      for(ii=0;ii<6;ii++)
      {
        if(elem->iconf[i][ii]==-1)
        {
          elem->iconf[i][ii]=0;      /*-1:PLASTIC 0:ELASTIC 1:HINGE*/
        }
      }
      if((elem->iconf[i][0]== 1 ||
          elem->iconf[i][0]==-2 ||
          elem->iconf[i][0]==-3) &&
         elem->sect->Ixx==0.0 &&
         elem->sect->Iyy==0.0 &&
         elem->sect->Jzz==0.0) /*FOR BRACE*/
      {
        elem->iconf[i][0]=0;         /*-1:PLASTIC 0:ELASTIC 1:HINGE*/

if(fout!=NULL) fprintf(fout,"ELEM%d DISLOADED.\n",elem->code);
      }
    }
  }

  for(i=0;i<6;i++)
  {
    (melem+(elem->loff))->bond[0][i]=elem->iconf[0][i];
    (melem+(elem->loff))->bond[1][i]=elem->iconf[1][i];
  }


  if(wsurf.hwnd!=NULL)
  {
    fseek(arc.fsurface,0L,SEEK_END);
    for(i=0;i<2;i++)                                    /*HEAD,TAIL*/
    {                               /*0:Nz 1:Qx 2:Qy 3:Mz 4:Mx 5:My*/
      ys[0][0]=(elem->stress[i][0]-*(dstress+6*i+0)
				-pow(-1.0,i*1.0)*fc[0])/fu[0];
	  ys[0][1]=(elem->stress[i][1]-*(dstress+6*i+1)-fc[1])/fu[1];
	  ys[0][2]=(elem->stress[i][2]-*(dstress+6*i+2)-fc[2])/fu[2];
	  ys[0][3]=(elem->stress[i][3]-*(dstress+6*i+3)
				-pow(-1.0,i*1.0)*fc[3])/fu[3];
      ys[0][4]=(elem->stress[i][4]-*(dstress+6*i+4)-fc[4])/fu[4];
      ys[0][5]=(elem->stress[i][5]-*(dstress+6*i+5)-fc[5])/fu[5];

	  ys[1][0]=(elem->stress[i][0]-pow(-1.0,i*1.0)*fc[0])/fu[0];
      ys[1][1]=(elem->stress[i][1]-fc[1])/fu[1];
      ys[1][2]=(elem->stress[i][2]-fc[2])/fu[2];
	  ys[1][3]=(elem->stress[i][3]-pow(-1.0,i*1.0)*fc[3])/fu[3];
      ys[1][4]=(elem->stress[i][4]-fc[4])/fu[4];
      ys[1][5]=(elem->stress[i][5]-fc[5])/fu[5];

      line.code=elem->code;
      line.ends[0].d[0]=ys[0][SURFACEX];
      line.ends[0].d[1]=ys[0][SURFACEY];
      line.ends[0].d[2]=ys[0][SURFACEZ];
      line.ends[1].d[0]=ys[1][SURFACEX];
      line.ends[1].d[1]=ys[1][SURFACEY];
      line.ends[1].d[2]=ys[1][SURFACEZ];

      if(f[i]>=pow(RADIUS,EXPONENT)) /*YIELDED*/
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





void updatestressbc(struct memoryelem *melem,FILE *fout,FILE *fsrf,
                  double *edisp,double *dstress,double **estiff,
                  struct owire *elem,double func[],FILE *ftxt,double ncr)
/*ELEMENT STRESS UPDATE.*/
/*UJIOKA:REFLECTING BUCKLING CONDENSATION RESULT*/
{
  char s[80],string[256];
  char iconf[12];
  long int i,ii,j/*,jj*/,nn[2];
  double dL,lamda[2]={0.0,0.0};
  double fc[6],fu[6],f[2],dfdp[2][6],q[12][2],a[2][2],rate;
  double fe[2],dfdpe[2][6],qe[12][2],ae[2][2];     /*1 STEP BEFORE.*/
  double det,detinverse,function;
  double ys[2][6];
  double due[2][6],dup[2][6]/*,dp[2][6]*/;
  struct line line;


  nn[0]=elem->node[0]->code;
  nn[1]=elem->node[1]->code;

  /*FUNCTION 1 STEP BEFORE.*/
  coefficientsbc(*elem,estiff,fe,dfdpe,qe,ae,ncr);      //UJIOKA

  /*UPDATE STRESS.*/
  for(i=0;i<2;i++)
  {
    for(j=0;j<6;j++)
    {
      elem->stress[i][j]+=*(dstress+6*i+j);
      (melem+(elem->loff))->stress[i][j]=elem->stress[i][j];
    }
  }

  if(elem->stress[0][0]<=0)  ncr=0;
  coefficientsbc(*elem,estiff,f,dfdp,q,a,ncr);        /*UPDATED FUNCTION.*/

  for(i=0;i<2;i++)                 /*ICONF[2][6] INTO ICONF[12].*/
  {
    for(j=0;j<6;j++) iconf[6*i+j]=elem->iconf[i][j];
  }
  for(i=0;i<6;i++)                 /*CENTER,WIDTH OF YIELD SURFACE.*/
  {
    fc[i]=0.5*(elem->sect->fmax[i]+elem->sect->fmin[i]);
    fu[i]=0.5*(elem->sect->fmax[i]-elem->sect->fmin[i]);
  }

  /***UJIOKA:圧縮部材について、降伏曲面を縮約による座屈荷重(Ncr)に修正***/
  /*
  　引張部材について、ここでは便宜上Ncr=0になっており、降伏曲面の修正をしない
  　（圧縮側と引張側で降伏曲面の式を場合分けしていないため）　
  */
  if(ncr>0 && ncr<elem->sect->fmax[0])
  {
    fc[0]=0;
    fu[0]=ncr;
//    sprintf(string,"fu[0]=%.5E\n",fu[0]);
//    errormessage(string);

  /*
    fc[0]=0.5*(ncr+elem->sect->fmin[i]);
    fu[0]=0.5*(ncr-elem->sect->fmin[i]);
  */
  }
  /*
　　圧縮部材の座屈荷重が負の値となる場合：未検討、元の降伏曲面のまま
  */
  else if (ncr<0) /*negative buckling load in compression element:yet*/
  {
//     sprintf(string,"Negative Buckling Load in Compression Element\n");
//     errormessage(string);
  }
  /***UJIOKA***/

  for(i=0;i<2;i++)                                       /*INITIAL.*/
  {
    for(j=0;j<6;j++) dup[i][j]=0.0;
  }

  /*LAMDA OF PLASTIC END.*/
  if(ae[0][0]!=0.0 && ae[1][1]==0.0)       /*IF I:PLASTIC J:ELASTIC*/
  {
    for(i=0;i<12;i++)
    {
      if(iconf[i]!=1 && iconf[i]!=-2 && iconf[i]!=-3)
      {
        lamda[0]+=1.0/ae[0][0]*qe[i][0]*(*(edisp+i));
      }
    }
    for(j=0;j<6;j++)
    {
      dup[0][j]=lamda[0]*dfdpe[0][j];
	}
  }
  else if(ae[0][0]==0.0 && ae[1][1]!=0.0)  /*IF I:ELASTIC J:PLASTIC*/
  {
    for(i=0;i<12;i++)
    {
      if(iconf[i]!=1 && iconf[i]!=-2 && iconf[i]!=-3)
      {
        lamda[1]+=1.0/ae[1][1]*qe[i][1]*(*(edisp+i));
      }
    }
    for(j=0;j<6;j++)
    {
      dup[1][j]=lamda[1]*dfdpe[1][j];
    }
  }
  else if(ae[0][0]!=0.0 && ae[1][1]!=0.0)  /*IF I:PLASTIC J:PLASTIC*/
  {
    for(i=0;i<12;i++)
    {
      if(iconf[i]!=1 && iconf[i]!=-2 && iconf[i]!=-3)
      {
        detinverse=ae[0][0]*ae[1][1]-ae[0][1]*ae[1][0];
        if(detinverse==0.0)
        {
          det=0.0;                           /*UNDER CONSIDERATION.*/
		  errormessage("UPDATESTRESS:UNDER CONSIDERATION.");
          sprintf(string,"Update : Matrix Singular.SECT=%d",elem->sect->code);
/*
MessageBox(NULL,string,"Updatestress",MB_OK);
*/
        }
        else det=1.0/detinverse;

        lamda[0]+=det*( ae[1][1]*qe[i][0]*(*(edisp+i))
                       -ae[0][1]*qe[i][1]*(*(edisp+i)));
        lamda[1]+=det*(-ae[1][0]*qe[i][0]*(*(edisp+i))
                       +ae[0][0]*qe[i][1]*(*(edisp+i)));
      }
    }
    for(j=0;j<6;j++)
    {
      dup[0][j]=lamda[0]*dfdpe[0][j];
      dup[1][j]=lamda[1]*dfdpe[1][j];
    }
  }

  /*ELASTIC END.*/
  if(ae[0][0]==0.0 || ae[1][1]==0.0)    /*IF I:ELASTIC OR J:ELASTIC*/
  {
    if(elem->sect->Ixx==0.0 &&
       elem->sect->Iyy==0.0 &&
       elem->sect->Jzz==0.0)                            /*FOR BRACE*/
    {
      dL=*(edisp+6)-*(edisp+0); /*EXTENSION OF LENGTH*/

      if(dL< 0.0 && elem->iconf[0][0]==-2 && ae[0][0]==0.0)
      {
        lamda[0]=-1.0;                         /*LONGING DISLOADED.*/
if(fout!=NULL) fprintf(fout,"ELEM%d HEAD dL=%.5E LONGING DISLOADED.\n",
elem->code,dL);
      }
      if(dL< 0.0 && elem->iconf[1][0]==-2 && ae[1][1]==0.0)
      {
        lamda[1]=-1.0;                         /*LONGING DISLOADED.*/
if(fout!=NULL) fprintf(fout,"ELEM%d TAIL dL=%.5E LONGING DISLOADED.\n",
elem->code,dL);
      }
      if(dL>=0.0 && elem->iconf[0][0]==-3 && ae[0][0]==0.0)
      {
        if(ae[0][0]==0.0) lamda[0]=-1.0;      /*SHORTING DISLOADED.*/
if(fout!=NULL) fprintf(fout,"ELEM%d HEAD dL=%.5E SHORTING DISLOADED.\n",
elem->code,dL);
      }
      if(dL>=0.0 && elem->iconf[1][0]==-3 && ae[1][1]==0.0)
      {
        if(ae[1][1]==0.0) lamda[1]=-1.0;      /*SHORTING DISLOADED.*/
if(fout!=NULL) fprintf(fout,"ELEM%d TAIL dL=%.5E SHORTING DISLOADED.\n",
elem->code,dL);
      }
    }
    else if(ae[0][0]==0.0 && fe[0]>=f[0])
    {
      lamda[0]=-1.0;                                 /*I:DISLOADED.*/
    }
    else if(ae[1][1]==0.0 && fe[1]>=f[1])
    {
      lamda[1]=-1.0;                                 /*J:DISLOADED.*/
    }
  }

  /*dUe,dUp DECISION.*/
  for(i=0;i<2;i++)
  {
    for(j=0;j<6;j++)
    {
      due[i][j]=*(edisp+6*i+j)-dup[i][j];        /*{due}={du}-{dup}*/
    }
  }
  /*CHECK STRESS {dp}=[ke]{due}.*/
  /*for(i=0;i<2;i++)
  {
    for(j=0;j<6;j++)
    {
      dp[i][j]=0.0;
      for(ii=0;ii<2;ii++)
      {
        for(jj=0;jj<6;jj++)
        {
          dp[i][j]+=(*(*(estiff+6*i+j)+6*ii+jj))*due[ii][jj];
        }
      }
    }
  }*/

  /*fprintf(fout,"ELEM %d [Ke]{dup}\n",elem->code);
  for(i=0;i<12;i++)
  {
    data=0.0;
    for(ii=0;ii<2;ii++)
    {
      for(jj=0;jj<6;jj++)
      {
        data+=(*(*(estiff+i)+6*ii+jj))*dup[ii][jj];
      }
    }
    fprintf(fout,"%8.5f\n",data);
  }*/

  /*STRAIN ENERGY.*/
  for(j=0;j<6;j++)
  {
    for(i=0;i<2;i++)
    {
      elem->Ee[0]+=0.25*(2*elem->stress[i][j]-(*(dstress+6*i+j)))
                       *due[i][j];
      elem->Ee[1]+=0.25*(2*elem->stress[i][j]-(*(dstress+6*i+j)))
                       *due[i][j];
    }
    elem->Ep[0]+=(elem->stress[0][j])*dup[0][j];
    elem->Ep[1]+=(elem->stress[1][j])*dup[1][j];
  }
#if 0
  if(fout!=NULL)
  {
    fprintf(fout,"ELEM %3d",elem->code);
    fprintf(fout,"       du=         due+     dup");
    fprintf(fout,"         dp= [K]{due}");
    fprintf(fout,"          pm          pn");
    fprintf(fout,"           dEei           dEej");
    fprintf(fout,"       dEpi       dEpj\n");
    for(i=0;i<2;i++)
    {
      for(j=0;j<6;j++)
      {
        fprintf(fout,"         %8.5f=%12.5E+%8.5f/*  %9.5f=%9.5f*/",
                *(edisp+6*i+j),due[i][j],dup[i][j]/*,
                *(dstress+6*i+j),dp[i][j]*/);
        fprintf(fout," %11.7f %11.7f",
                (elem->stress[i][j]-(*(dstress+6*i+j))),
                elem->stress[i][j]);
        fprintf(fout," %14.7E %14.7E %10.7f %10.7f\n",
                0.25*(2*elem->stress[i][j]-(*(dstress+6*i+j)))
                       *due[i][j],
                0.25*(2*elem->stress[i][j]-(*(dstress+6*i+j)))
                       *due[i][j],
                (elem->stress[0][j])*dup[0][j],
                (elem->stress[1][j])*dup[1][j]);
      }
    }
    fprintf(fout,"Eei=%10.7E Eej=%10.7E Epi=%10.7E Epj=%10.7E\n",
			elem->Ee[0],elem->Ee[1],elem->Ep[0],elem->Ep[1]);
  }
#endif
  /*YIELD JUDGEMENT.*/
  for(i=0;i<2;i++)
  {
    func[i]=f[i];

    if(f[i]>=pow(RADIUS,EXPONENT))                       /*YIELDED.*/
    {
      sprintf(string,"YIELDED:ELEM%d NODE%ld SECT%d",
              elem->code,nn[i],elem->sect->code);
      /*errormessage(string);*/
      strcat(string," FUNCTION:{");
      function=0.0;
      for(ii=0;ii<6;ii++)
      {
        if(ii==0 || ii==3) /*N,Mz*/
        {
		  rate=(elem->stress[i][ii]-pow(-1.0,i*1.0)*fc[ii])/fu[ii];
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
          if(elem->sect->Ixx==0.0 &&
             elem->sect->Iyy==0.0 &&
             elem->sect->Jzz==0.0) /*FOR BRACE*/
          {
            /*elem->sect->area=0.0;*/    /*WRONG.ALL ELEM WITH THIS*/
                                         /*SECT WILL BE SET 0.     */

            dL=*(edisp+6)-*(edisp+0); /*EXTENSION OF LENGTH*/

            if(ii==0)
            {
              if(dL>=0.0 && ((i==0 && rate<0.0) ||
                             (i==1 && rate>0.0)))
              {
                elem->iconf[i][0]=-2;
if(fout!=NULL) fprintf(fout,"ELEM%d dL=%.5E LONGING.\n",
elem->code,dL);
              }
              else if(dL< 0.0 && ((i==0 && rate>0.0) ||
                                  (i==1 && rate<0.0)))
              {
                elem->iconf[i][0]=-3;
if(fout!=NULL) fprintf(fout,"ELEM%d dL=%.5E SHORTING.\n",
elem->code,dL);
              }
              /*elem->iconf[i][0]=1;*/       /*-2:LONGING -3:SHORTING*/
            }
          }
          else
          {
            elem->iconf[i][ii]=-1;   /*-1:PLASTIC 0:ELASTIC 1:HINGE*/
          }
        }
      }
      sprintf(s,"}=%8.5f",function); /*VALUE OF FUNCTION.*/
      strcat(string,s);
      /*errormessage(string);*/
      if(fout!=NULL) fprintf(fout,"%s\n",string);
      if(fsrf!=NULL) fprintf(fsrf,"%s\n",string);    /*for surface.txt*/
    }
    else        /*for surface.txt*/
    {
      sprintf(string,"ELEM%d NODE%ld SECT%d",
              elem->code,nn[i],elem->sect->code);
      /*errormessage(string);*/
      strcat(string," FUNCTION:{");
      function=0.0;
      for(ii=0;ii<6;ii++)
      {
        if(ii==0 || ii==3) /*N,Mz*/
        {
		  rate=(elem->stress[i][ii]-pow(-1.0,i*1.0)*fc[ii])/fu[ii];
		}
        else /*Qx,Qy,Mx,My*/
        {
          rate=(elem->stress[i][ii]-fc[ii])/fu[ii];
        }
        function+=pow(fabs(rate),EXPONENT);    /*VALUE OF FUNCTION.*/
        sprintf(s," %8.5f",rate);
        strcat(string,s);
      }
      sprintf(s,"}=%8.5f",function); /*VALUE OF FUNCTION.*/
      strcat(string,s);
      /*errormessage(string);*/
//      if(fout!=NULL) fprintf(fout,"%s\n",string);
      if(fsrf!=NULL) fprintf(fsrf,"%s\n",string);    /*for surface.txt*/
    }
  }

/*if(fout!=NULL) fprintf(fout,"ELEM%d Lamda1=%.3E Lamda2=%.3E\n",
elem->code,lamda[0],lamda[1]);*/

  /*DISLOAD JUDGEMENT.*/
  for(i=0;i<2;i++)
  {
    if(lamda[i]<0.0)                                   /*DISLOADED.*/
    {
      sprintf(string,"DISLOAD:ELEM%d NODE%ld SECT%d",
              elem->code,nn[i],elem->sect->code);
      /*errormessage(string);*/
      /*if(fout!=NULL) fprintf(fout,"%s\n",string);*/
      for(ii=0;ii<6;ii++)
      {
        if(elem->iconf[i][ii]==-1)
        {
          elem->iconf[i][ii]=0;      /*-1:PLASTIC 0:ELASTIC 1:HINGE*/
        }
      }
      if((elem->iconf[i][0]== 1 ||
          elem->iconf[i][0]==-2 ||
          elem->iconf[i][0]==-3) &&
         elem->sect->Ixx==0.0 &&
         elem->sect->Iyy==0.0 &&
         elem->sect->Jzz==0.0) /*FOR BRACE*/
      {
        elem->iconf[i][0]=0;         /*-1:PLASTIC 0:ELASTIC 1:HINGE*/

if(fout!=NULL) fprintf(fout,"ELEM%d DISLOADED.\n",elem->code);
      }
    }
  }

  for(i=0;i<6;i++)
  {
    (melem+(elem->loff))->bond[0][i]=elem->iconf[0][i];
    (melem+(elem->loff))->bond[1][i]=elem->iconf[1][i];
  }

  /*if(ftxt!=NULL)
  {
    for(i=0;i<2;i++)
    {
	  ys[1][0]=(elem->stress[i][0]-pow(-1.0,i)*fc[0])/fu[0];
      ys[1][1]=(elem->stress[i][1]-fc[1])/fu[1];
      ys[1][2]=(elem->stress[i][2]-fc[2])/fu[2];
      ys[1][3]=(elem->stress[i][3]-pow(-1.0,i)*fc[3])/fu[3];
      ys[1][4]=(elem->stress[i][4]-fc[4])/fu[4];
      ys[1][5]=(elem->stress[i][5]-fc[5])/fu[5];

      fprintf(ftxt," ELEM%d-%d %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f",
              elem->code,i+1,
              ys[1][0],ys[1][1],ys[1][2],ys[1][3],ys[1][4],ys[1][5]);
    }
    fprintf(ftxt,"\n");
  }*/

/*if(elem->code==1054){*/
  if(wsurf.hwnd!=NULL)
  {
    fseek(arc.fsurface,0L,SEEK_END);
    for(i=0;i<2;i++)                                    /*HEAD,TAIL*/
    {                               /*0:Nz 1:Qx 2:Qy 3:Mz 4:Mx 5:My*/
      ys[0][0]=(elem->stress[i][0]-*(dstress+6*i+0)
				-pow(-1.0,i*1.0)*fc[0])/fu[0];
	  ys[0][1]=(elem->stress[i][1]-*(dstress+6*i+1)-fc[1])/fu[1];
	  ys[0][2]=(elem->stress[i][2]-*(dstress+6*i+2)-fc[2])/fu[2];
	  ys[0][3]=(elem->stress[i][3]-*(dstress+6*i+3)
				-pow(-1.0,i*1.0)*fc[3])/fu[3];
      ys[0][4]=(elem->stress[i][4]-*(dstress+6*i+4)-fc[4])/fu[4];
      ys[0][5]=(elem->stress[i][5]-*(dstress+6*i+5)-fc[5])/fu[5];

	  ys[1][0]=(elem->stress[i][0]-pow(-1.0,i*1.0)*fc[0])/fu[0];
      ys[1][1]=(elem->stress[i][1]-fc[1])/fu[1];
      ys[1][2]=(elem->stress[i][2]-fc[2])/fu[2];
	  ys[1][3]=(elem->stress[i][3]-pow(-1.0,i*1.0)*fc[3])/fu[3];
      ys[1][4]=(elem->stress[i][4]-fc[4])/fu[4];
      ys[1][5]=(elem->stress[i][5]-fc[5])/fu[5];

      line.code=elem->code;
      line.ends[0].d[0]=ys[0][SURFACEX];
      line.ends[0].d[1]=ys[0][SURFACEY];
      line.ends[0].d[2]=ys[0][SURFACEZ];
      line.ends[1].d[0]=ys[1][SURFACEX];
      line.ends[1].d[1]=ys[1][SURFACEY];
      line.ends[1].d[2]=ys[1][SURFACEZ];

      if(f[i]>=pow(RADIUS,EXPONENT)) /*YIELDED*/
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

/*sprintf(string,"ELEM%d\np=%.3f %.3f %.3f %.3f %.3f %.3f\nf=%.3f %.3f %.3f %.3f %.3f %.3f",
        elem->code,
        elem->stress[0][0],
        elem->stress[0][1],
        elem->stress[0][2],
        elem->stress[0][3],
        elem->stress[0][4],
        elem->stress[0][5],
        ys[0][0],ys[0][1],ys[0][2],ys[0][3],ys[0][4],ys[0][5]);
MessageBox(NULL,string,"TEST",MB_OK);*/
/*}*//*TEST*/

  return;
}/*updatestressbc*/

void updatestressnl(struct memoryelem *melem,double *dstress,
					 struct owire *elem)
/*ELEMENT STRESS UPDATE FOR NON-LINEAR ANALYSIS.*/
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
}/*updatestressnl*/

void outputdisp(double *gvct,FILE *fout,int nnode,
				struct onode *nodes)
/*OUTPUT NODE DISPLACEMENT.*/
{
  char string[256];
  int i,j;
  double data[6];

  for(i=0;i<nnode;i++)
  {
	for(j=0;j<6;j++) data[j]=*(gvct+6*i+j);
	sprintf(string,
			"%4d %10.6e %10.6e %10.6e %11.7e %11.7e %11.7e",
			(nodes+i)->code,
			data[0],data[1],data[2],data[3],data[4],data[5]);
	if(fout!=NULL) fprintf(fout,"%s\n",string);
  }
}/*outputdisp*/

/*** 09.08.28 araki for Tsurukawa *********************************************/
void outputdisp02(double *gvct,double *displong,FILE *fout,int nnode,
                  struct onode *nodes)
/*OUTPUT NODE DISPLACEMENT IN OHX,OHY STYLE.*/
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
    for(j=0;j<=5;j++) data[j]=*(gvct+6*(i-1)+j)-*(displong+6*(i-1)+j);
    sprintf(string,
			"%4d %10.6f %10.6f %10.6f %11.7f %11.7f %11.7f",
			(nodes+i-1)->code,
			data[0],data[1],data[2],data[3],data[4],data[5]);
	if(fout!=NULL) fprintf(fout,"%s\n",string);
  }
}/*outputdisp02*/
/******************************************************************************/

void updateform(double *ddisp,double *gvct,int nnode)
/*FORMATION UPDATE.*/
{
  int i,j;
  long int loff;

  for(i=0;i<nnode;i++)
  {
	for(j=0;j<6;j++)
	{
	  loff=6*i+j;
	  *(ddisp+loff)+=*(gvct+loff);                       /*{U}+{dU}*/
	}
  }
  return;
}/*updateform*/

void updateform2(double *ddisp,double *gvct,int nnode,long int *moff)
/*FORMATION UPDATE.*/
{
  int i,j,n;

  for(i=0;i<nnode;i++)                                   /*{U}+{dU}*/
  {
	if(moff==NULL) n=i;
	else           n=*(moff+i);
	for(j=0;j<6;j++) *(ddisp+6*i+j)+=*(gvct+6*n+j);
  }
  return;
}/*updateform2*/


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
				  double *estress,FILE *fout,double func[])
/*ELEMENT STRESS OUTPUT.*/
{
  char s[80],string[256];
  long int i,n,nn[2];

  for(n=1;n<=2;n++)
  {
	nn[n-1]=elem.node[n-1]->code;
	sprintf(string,"ELEM %5d SECT %4d NODE %4ld:",
			elem.code,elem.sect->code,nn[n-1]);
	for(i=0;i<6;i++)
	{
/*      sprintf(s," %9.3f/%9.3f",
			  *(estress+6*(n-1)+i),elem.stress[n-1][i]);
*/
	 sprintf(s," %12.8f/%12.8f",
			  *(estress+6*(n-1)+i),elem.stress[n-1][i]);
	 strcat(string,s);
	}
	if(fout!=NULL) fprintf(fout,"%s f=%6.3f\n",string,func[n-1]);
  }
  return;
}/*outputstress*/

void outputshellstress(struct oshell shell,
    double* estress, FILE* fout)
	/*ELEMENT STRESS OUTPUT.*/
{
    char s[80], string[256];
    long int i, n;
    /*
    for (n = 0; n < shell.nnod; n++)
    {
        if (n == 0)
        {
            sprintf(string, "%5d %4d %4d",
                shell.code, shell.sect->code, shell.node[n]->code);
        }
        else
		{
            sprintf(string, "           %4d", shell.node[n]->code);
        }

        for (i = 0; i < 6; i++)
		{
            sprintf(s, " %12.8f", shell.stress[n][i]);
            strcat(string, s);
        }
        if (fout != NULL) fprintf(fout, "%s\n", string);
	}
    */

	sprintf(string, "%5d",shell.code);
    for (i = 0; i < 6; i++)
    {
		sprintf(s, " %e", shell.stress[0][i]);
        strcat(string, s);
    }
	if (fout != NULL) fprintf(fout, "%s\n", string);
    return;
}/*outputshellstress*/


void outputstress02(struct owire elem,
                    double *estress,double *stresslong,int i,FILE *fout,double func[])
/*ELEMENT STRESS OUTPUT IN OHX,OHY STYLE.*/
{
  char s[80],string[256];
  long int j,n,nn[2];
  double stress;

  for(n=1;n<=2;n++)
  {
    nn[n-1]=elem.node[n-1]->code;
    if(n==1) sprintf(string,"%5d %4d %4d",
                     elem.code,elem.sect->code,nn[n-1]);
    if(n==2) sprintf(string,"           %4d",nn[n-1]);
    for(j=0;j<6;j++)
    {
	if(*stresslong!=NULL)
      stress=elem.stress[n-1][j]-*(stresslong+12*i+j+6*(n-1));                   //araki?
    else
      stress=elem.stress[n-1][j];                   //araki?
      /*sprintf(s," %9.8f",*(estress+6*(n-1)+i));*/       /*WITHOUT CMQ*/
      sprintf(s," %9.8f",stress);                            /*WITH CMQ*/
      strcat(string,s);
    }
    if(fout!=NULL) fprintf(fout,"%s\n",string);
  }
  return;
}/*outputstress02*/

void outputstressnl(struct owire elem,                           //Check MIHARA
                     double *estress,FILE *fout)
/*ELEMENT STRESS OUTPUT.*/
{
  char s[80],string[256];
  int i,n,nn[2];

  for(n=0;n<2;n++)
  {
    nn[n]=elem.node[n]->code;
    if(n==0) sprintf(string,"%5d %4d %4d",
                     elem.code,elem.sect->code,nn[n]);
    if(n==1) sprintf(string,"           %4d",nn[n]);
	for(i=0;i<6;i++)
    {
      /*sprintf(s," %9.3f",*(estress+6*n+i));*/       /*WITHOUT CMQ*/
      sprintf(s," %3.9f",elem.stress[n][i]);             /*WITH CMQ*/      //9.3
      strcat(string,s);
    }
    if(fout!=NULL) fprintf(fout,"%s\n",string);
  }
  return;
}/*outputstressnl*/

void outputreaction(struct gcomponent *gmtx,
                    double *gvct,
					struct onode *nodes,
                    struct oconf *confs,
                    double *dreact,FILE *fout,int nnode)
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
      dreaction=0.0;
      for(j=1;j<=6*nnode;j++)
	  {
        gread(gmtx,i,j,&gstiff);
		dreaction+=gstiff*(*(gvct+j-1));             /*{dR}=[K]{dU}*/
      }

      *(dreact+nreact)+=dreaction;               /*UPDATE REACTION.*/
	  /*reaction=*(dreact+nreact-1);*//*ERROR*/          /*{R}+{dR}*/
	  reaction=*(dreact+nreact);                         /*{R}+{dR}*/

	  sprintf(string,"NODE:%5ld %ld %12.5f / %12.5f",
			  (nodes+offset)->code,(i-1)%6+1,dreaction,reaction);
	  if(fout!=NULL) fprintf(fout,"%s\n",string);


	  nreact++;
	}
  }
  return;
}/*outputreaction*/

void outputreactionnl(struct gcomponent *gmtx,
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

	//if(iconf==1)
	//{
	//  dreaction=0.0;
	//  for(j=1;j<=msize;j++)
	//  {
	//    gread(gmtx,i,j,&gstiff);
	//    dreaction+=gstiff*(*(gvct+j-1));             /*{dR}=[K]{dU}*/
	//  }

	//  *(dreact+nreact)+=dreaction;                       /*{R}+{dR}*/

	//  sprintf(string,"%4ld %10ld %14.6f     1",
	//          (nodes+offset)->code,(i-1)%6+1,dreaction);
										  /*1:ONLY FIXED AVAILABLE.*/
	//  if(fout!=NULL) fprintf(fout,"%s\n",string);

	//  nreact++;
	//}

    if(iconf==1)    //modyfied by fukushima 150224
    {
	  sprintf(string,"%4ld %10ld %14.6f     1",
			  (nodes+offset)->code,(i-1)%6+1,*(dreact+nreact));
										  /*1:ONLY FIXED AVAILABLE.*/
	  if(fout!=NULL) fprintf(fout,"%s\n",string);

      nreact++;
    }

	currentpivot(i,msize);
  }
  return;
}/*outputreactionnl*/

void outputreaction02(struct gcomponent *gmtx,
                      double *gvct,
                      struct onode *nodes,
                      struct oconf *confs,
                      double *dreact,double *reactlong,FILE *fout,int nnode)
/*REACTIONS OUTPUT IN OHX,OHY STYLE.*/
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
      reaction=*(dreact+nreact)-*(reactlong+nreact);     /*Σ{dR}={R}-{RL}*/

      sprintf(string,"%4ld %10ld %14.6f     1",
              (nodes+offset)->code,(i-1)%6+1,reaction);
                                          /*1:ONLY FIXED AVAILABLE.*/
      if(fout!=NULL) fprintf(fout,"%s\n",string);

      nreact++;
    }
  }
  return;
}/*outputreaction02*/
/******************************************************************************/

/*** 12.05.29 araki for Hakushima *********************************************/
void outputreactionwithoutupdate(struct gcomponent *gmtx,
                                 double *gvct,
                                 struct onode *nodes,
								 struct oconf *confs,
								 double *dreact,FILE *fout,FILE *ftxt,
                                 int nnode,int mhoko)
/*REACTIONS OUTPUT WITHOUT UPDATING.*/
{
  char string[256];
  char iconf;
  int offset;
  long int i,j,nreact=0;
  double gstiff,reaction,dreaction;
  double sreaction=0;

  for(i=1;i<=6*nnode;i++)
  {
    iconf=(confs+i-1)->iconf;
    offset=(i-1)/6;

    if(iconf==1)
	{
      dreaction=0.0;
      for(j=1;j<=6*nnode;j++)
      {
		gread(gmtx,i,j,&gstiff);
        dreaction+=gstiff*(*(gvct+j-1));             /*{dR}=[K]{dU}*/
      }
      /*reaction=*(dreact+nreact-1);*/                       /*{R}+{dR}*/
      reaction=*(dreact+nreact);                       /*{R}+{dR}*/             //araki

      sprintf(string,"NODE:%5ld %ld %12.5f / %12.5f",
              (nodes+offset)->code,(i-1)%6+1,dreaction,reaction);
      if(fout!=NULL) fprintf(fout,"%s\n",string);
      if((i-1)%6+1==mhoko+1) sreaction-=reaction;
      nreact++;
    }
  }
  if(ftxt!=NULL) fprintf(ftxt," %10.5f",sreaction);
  return;
}/*outputreactionwithoutupdate*/
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
  ppen=(HPEN)SelectObject(hdc,hpen);

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

struct oprop *addproperty(struct oprop *op,
                          struct organ *org)
/*ADD PROPERTY OF ORGAN FRAME.*/
{
  int i,j;
  long int pprop,loff;
  struct oprop *p,*props;

  pprop=org->nprop;
  org->nprop++;

  /*COPY PROPERTIES.*/
  props=(struct oprop *)malloc((org->nprop)*sizeof(struct oprop));
  for(i=0;i<pprop;i++)
  {
    *(props+i)=*(org->props+i);
    (props+i)->name=(char *)
                    malloc((strlen((org->props+i)->name)+1)
                           *sizeof(char));
    strcpy((props+i)->name,(org->props+i)->name);
  }

  i=org->nprop-1;
  *(props+i)=*op;
  (props+i)->loff=i;
  (props+i)->name=(char *)
                  malloc((strlen(op->name)+1)*sizeof(char));
  strcpy((props+i)->name,op->name);

  /*MODIFY PROPERTY POINTERS.*/
  for(i=0;i<(org->nsect);i++)
  {
    for(j=0;j<(org->sects+i)->nfig;j++)
    {
      loff=((org->sects+i)->figs+j)->loff;
      ((org->sects+i)->figs+j)->prop=(props+loff);
    }
  }

  /*FREE PREVIOUS PROPERTIES.*/
  for(i=0;i<pprop;i++)
  {
    free((org->props+i)->name);
  }
  free(org->props);

  org->props=props;
  p=org->props+(org->nprop-1); /*RETURN LAST.*/

  return p;
}/*addproperty*/

struct osect *addsection(struct osect *os,
                         struct organ *org)
/*ADD SECTION OF ORGAN FRAME.*/
{
  int i;
  long int psect,loff;
  struct osect *s,*sects;

  psect=org->nsect;
  org->nsect++;

  /*COPY SECTIONS.*/
  sects=(struct osect *)malloc((org->nsect)*sizeof(struct osect));
  for(i=0;i<psect;i++)
  {
    *(sects+i)=*(org->sects+i);
    (sects+i)->name=(char *)
                    malloc((strlen((org->sects+i)->name)+1)
                           *sizeof(char));
    strcpy((sects+i)->name,(org->sects+i)->name);

    /*COPY POLYPOLYCURVE ?*/
  }

  /*ADD INTO BOTTOM*/
  i=org->nsect-1;
  *(sects+i)=*os;
  (sects+i)->loff=i;
  (sects+i)->name=(char *)
                  malloc((strlen(os->name)+1)*sizeof(char));
  strcpy((sects+i)->name,os->name);

  /*COPY POLYPOLYCURVE ?*/

  /*MODIFY SECTION POINTERS.*/
  for(i=0;i<(org->nelem);i++)
  {
    loff=(org->elems+i)->sect->loff;
    ((org->elems+i)->sect)=(sects+loff);
  }

  /*FREE PREVIOUS SECTIONS.*/
  /*for(i=0;i<psect;i++)
  {
    free((org->sects+i)->name);
  }
  free(org->sects);*/

  org->sects=sects;
  s=org->sects+(org->nsect-1); /*RETURN LAST.*/

  return s;
}/*addsection*/

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
                      struct oconf *oc,
                      struct organ *orgbefore)
/*ADD NODE OF ORGAN.*/
/*RETURN:ADDED NODE POINTER TO NODES OF ORGAN.*/
{
  /*char non[80],none[256];*/
  int i,ii,j,jj,k;
  long int nnode,nelem;
  struct onode *nodes;
  struct oconf c[6];

  nnode=orgbefore->nnode;
  nelem=orgbefore->nelem;

  if(oc==NULL)
  {
    for(i=0;i<6;i++)
    {
      c[i].iconf=0;
      c[i].value=0.0;
    }
  }
  else
  {
    for(i=0;i<6;i++)
    {
      c[i].iconf=(oc+i)->iconf;
      c[i].value=(oc+i)->value;
    }
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

    for(j=0;j<nelem;j++) /*CORRECT NODE POINTER OF ELEMENT*/
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

struct onode *deletenode(long int nodecode,
                         struct organ *orgbefore)
/*DELETE FREE NODE OF ORGAN.*/
/*RETURN:NODES OF ORGAN.*/
{
  int i,j,jj;
  long int nnode,nelem,nodeoffset;
  struct onode *nodes;

  nnode=orgbefore->nnode;
  nelem=orgbefore->nelem;

  for(i=0;i<nnode;i++)
  {
    if((orgbefore->nodes+i)->code==nodecode)
    {
      nodeoffset=(orgbefore->nodes+i)->loff;
      break;
    }
  }

  deleteconf(nodeoffset,orgbefore);

  nnode--;
  nodes=(struct onode *)malloc(nnode*sizeof(struct onode));
  if(nodes==NULL)
  {
    MessageBox(NULL,"Buffer Null.","DeleteNode",MB_OK);
    return orgbefore->nodes;
  }
  orgbefore->nnode=nnode;

  for(i=0;i<nodeoffset;i++)
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
  for(i=nodeoffset;i<nnode;i++)
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
  /*nodesbefore=orgbefore->nodes;*/
  /*free(orgbefore->nodes);*/
  orgbefore->nodes=nodes;

  return nodes;
}/*deletenode*/

struct onode *deletenodewithelem(long int nodecode,
                                 struct organ *orgbefore)
/*DELETE NODE OF ORGAN WITH ELEMENT.*/
/*RETURN:NODES OF ORGAN.*/
{
  int i,j,jj;
  long int nnode,nelem,nodeoffset;
  struct onode *nodes;

  for(j=0;j<(orgbefore->nelem);j++) /*DELETE JOINTED ELEMENT.*/
  {
    for(jj=0;jj<((orgbefore->elems+j)->nnod);jj++)
    {
      if((*((orgbefore->elems+j)->nods+jj))->code==nodecode)
      {
        deleteelement(j,orgbefore,nodecode);
        j--;
        break;
      }
    }
  }

  nnode=orgbefore->nnode;
  nelem=orgbefore->nelem;

  for(i=0;i<nnode;i++)
  {
    if((orgbefore->nodes+i)->code==nodecode)
    {
      nodeoffset=(orgbefore->nodes+i)->loff;
      break;
    }
  }

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
  for(i=nodeoffset;i<nnode;i++)
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

  /*free(orgbefore->nodes);*/
  orgbefore->nodes=nodes;

  return nodes;
}/*deletenodewithelem*/

long int singulatenode(struct organ *orgbefore)
/*SINGULATE CLOSE NODES.*/
/*RETURN:NODES OF ORGAN.*/
{
  char str[256];
  int i,j,ii,jj,k;
  long int nnode,nelem,ncode,icount,scount;
  double dx,dy,dz,dl,th;

  nnode=orgbefore->nnode;
  nelem=orgbefore->nelem;

  th=0.010; /*THRESHOLD [m]*/

  scount=0;
  for(i=0;i<nnode;i++)
  {
	currentpivot(i+1,nnode);

	for(j=i+1;j<nnode;j++)
	{
	  dx=(orgbefore->nodes+j)->d[0]-(orgbefore->nodes+i)->d[0];
	  dy=(orgbefore->nodes+j)->d[1]-(orgbefore->nodes+i)->d[1];
	  dz=(orgbefore->nodes+j)->d[2]-(orgbefore->nodes+i)->d[2];

	  dl=sqrt(dx*dx+dy*dy+dz*dz);

	  if(dl<=th)
	  {
		ncode=(orgbefore->nodes+j)->code;
		scount++;


		sprintf(str,"Node = %ld into %ld",ncode,(orgbefore->nodes+i)->code);
        errormessage(str);
		//MessageBox(NULL,str,"Singulate Node",MB_OK);


		for(ii=0;ii<nelem;ii++)
		{
		  for(jj=0;jj<((orgbefore->elems+ii)->nnod);jj++)
		  {
			if((*((orgbefore->elems+ii)->nods+jj))->code==ncode)
			{
			  *((orgbefore->elems+ii)->nods+jj)=orgbefore->nodes+i;
			}
		  }

          /***modified by UJIOKA***/
          if(((orgbefore->elems+ii)->nban)>=1)
	      {
	   	  	for(jj=0;jj<(orgbefore->elems+ii)->nban;jj++)
   		   	{
          		for(k=0;k<((orgbefore->elems+ii)->bans+jj)->nnod;k++)
            	{
              		if((*(((orgbefore->elems+ii)->bans+jj)->nods+k))->code==ncode)
                    {
                    	*(((orgbefore->elems+ii)->bans+jj)->nods+k)=orgbefore->nodes+i;
                    }
            	}
           	}
          }
          /***UJIOKA***/
        }
	  }
	}
  }
  sprintf(str,"Listed Nodes = %ld",scount);
  //MessageBox(NULL,str,"Singulate Node",MB_OK);

  scount=0;
  for(i=0;i<nnode;i++)
  {
	currentpivot(i+1,nnode);
	ncode=(orgbefore->nodes+i)->code;

	icount=0;
	for(ii=0;ii<nelem;ii++)
	{
	  for(jj=0;jj<((orgbefore->elems+ii)->nnod);jj++)
	  {
		if((*((orgbefore->elems+ii)->nods+jj))->code==ncode) icount++;
	  }
	}
	if(icount==0)
	{
	  deletenode(ncode,orgbefore);
	  nnode--;
	  scount++;
	  i--;


	  sprintf(str,"Deleted Node = %ld",ncode);
      errormessage(str);
	  //MessageBox(NULL,str,"Singulate Node",MB_OK);

	}
  }

  return scount;
}/*singulatenode*/

struct oelem *addelement(struct oelem *elem,
                         struct organ *orgbefore)
/*ADD ELEMENT OF ORGAN.OFFSET KNOWN,NODES ALREADY ADDED.*/
{
  long int i,j,nelem;
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
    *(elems+i)=*(orgbefore->elems+i);

	if(nmultielem>0)
    {
      for(j=0;j<nmultielem;j++)
      {
        if(*(multielem+j)==(orgbefore->elems+i))
        {
          *(multielem+j)=elems+i;
        }
      }
    }
  }
  for(i=(elem->loff);i<(nelem-1);i++)
  {
    *(elems+i+1)=*(orgbefore->elems+i);
    ((elems+i+1)->loff)++;

    if(nmultielem>0)
    {
      for(j=0;j<nmultielem;j++)
      {
        if(*(multielem+j)==(orgbefore->elems+i))
        {
          *(multielem+j)=elems+i;
        }
      }
    }
  }

  *(elems+(elem->loff))=*elem;

  free(orgbefore->elems);

  orgbefore->nelem=nelem;
  orgbefore->elems=elems;

  return elems;
}/*addelement*/

int addelementwithnode(struct oelem *pe,
                       struct organ *org,
                       double dx,double dy,double dz) /*INCREMENT*/
/*COPY ELEMENT WITH NODE BY INCREMENT {Dx,Dy,Dz}.*/
/*PE POINTING ELEMENT OF ORG.*/
{
  int i,j,k,ifind,nnode,loff;
  struct onode cnode;
  /*struct oelem ge={0,0,
                   ROLENULL,TYPENULL,
                   0,0,0.0,
                   {{0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}},
                   {{0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}},
                   NULL,NULL,NULL,NULL,
                   {{255,255,255},{255,255,255},{255,255,255}}
                  };*/ /*CREATING ELEMENT.*/
  struct oconf *oc;

  gelem=*pe;

  gelem.nods=(struct onode **)
             malloc((pe->nnod)*sizeof(struct onode *));
  gelem.bonds=(signed char *)
              malloc(6*(pe->nnod)*sizeof(signed char));
  gelem.bans=(struct obans *)
             malloc((pe->nban)*sizeof(struct obans));
  for(i=0;i<(pe->nban);i++)
  {
    *(gelem.bans+i)=*(pe->bans+i); /*COPY BANS*/
    (gelem.bans+i)->nods=(struct onode **)
                         malloc((pe->bans+i)->nnod
                         *sizeof(struct onode *));
  }

  oc=(struct oconf *)malloc(6*sizeof(struct oconf));
  for(i=0;i<pe->nnod;i++)
  {
    nnode=org->nnode;
    loff=(*(pe->nods+i))->loff;

    cnode=**(pe->nods+i);

    /*ADD NODES*/
    ifind=0;
    for(j=0;j<nnode;j++) /*OVERLAPPED NODE*/
    {
      if(cnode.d[0]+dx>(org->nodes+j)->d[0]-0.001 &&
         cnode.d[0]+dx<(org->nodes+j)->d[0]+0.001 &&
         cnode.d[1]+dy>(org->nodes+j)->d[1]-0.001 &&
         cnode.d[1]+dy<(org->nodes+j)->d[1]+0.001 &&
         cnode.d[2]+dz>(org->nodes+j)->d[2]-0.001 &&
         cnode.d[2]+dz<(org->nodes+j)->d[2]+0.001)
      {
        *(gelem.nods+i)=org->nodes+j;
        ifind=1;
        break;
      }
    }
    if(!ifind) /*NEW NODE*/
    {
      cnode.code=(org->nodes+nnode-1)->code+1;
      cnode.loff=nnode;

      for(j=0;j<6;j++) /*COPY CONFS*/
      {
        *(oc+j)=*(org->confs+6*loff+j);
      }

      *(gelem.nods+i)=addnode(cnode,oc,org);

      /*MOVE NODES*/
      (*(gelem.nods+i))->d[0]+=dx;
      (*(gelem.nods+i))->d[1]+=dy;
      (*(gelem.nods+i))->d[2]+=dz;
    }

    for(j=0;j<6;j++) /*COPY BONDS*/
    {
      *(gelem.bonds+6*i+j)=*(pe->bonds+6*i+j);
    }
  }
  free(oc);
/*
sprintf(str,"Elem %d\n",gelem.code);
for(i=0;i<gelem.nnod;i++)
{
  sprintf(s,"Node %d {%.3f,%.3f,%.3f}\n",
          (*(gelem.nods+i))->code,
          (*(gelem.nods+i))->d[0],
          (*(gelem.nods+i))->d[1],
          (*(gelem.nods+i))->d[2]);
  strcat(str,s);
}
MessageBox(NULL,str,"AddElement",MB_OK);
*/
  for(i=0;i<(pe->nban);i++) /*CORRECT POINTER*/
  {
    for(j=0;j<(pe->bans+i)->nnod;j++)
    {
      for(k=0;k<(pe->nnod);k++)
      {
        if(*((pe->bans+i)->nods+j)==*(pe->nods+k))
        {
          *((gelem.bans+i)->nods+j)=*(gelem.nods+k);
        }
      }
    }
  }

  gelem.code=(org->elems+(org->nelem-1))->code+1;
  gelem.loff=org->nelem;

  addelement(&gelem,org);

  gelem.nnod=0;
  gelem.nban=0;

  gelem.nods=NULL;
  gelem.bonds=NULL;
  gelem.bans=NULL;

  return 1;
}/*addelementwithnode*/

int divideatmidpoint(struct oelem *pe,        //by fukushima
                       struct organ *org)
/*DIVIDE AT MIDPOINT.*/
/*PE POINTING ELEMENT OF ORG.*/
{
  if(pe->nnod!=2)
  {
    return 1;
  }

  int i,j,k,ifind,nnode,loff;
  struct onode cnode1,cnode2;
  double dx,dy,dz;
  /*struct oelem ge={0,0,
                   ROLENULL,TYPENULL,
                   0,0,0.0,
                   {{0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}},
                   {{0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}},
                   NULL,NULL,NULL,NULL,
                   {{255,255,255},{255,255,255},{255,255,255}}
                  };*/ /*CREATING ELEMENT.*/
  struct oconf *oc;

  gelem=*pe;

  gelem.nods=(struct onode **)
             malloc((pe->nnod)*sizeof(struct onode *));
  gelem.bonds=(signed char *)
              malloc(6*(pe->nnod)*sizeof(signed char));
  gelem.bans=(struct obans *)
             malloc((pe->nban)*sizeof(struct obans));
  for(i=0;i<(pe->nban);i++)
  {
    *(gelem.bans+i)=*(pe->bans+i); /*COPY BANS*/
    (gelem.bans+i)->nods=(struct onode **)
                         malloc((pe->bans+i)->nnod
                         *sizeof(struct onode *));
  }

  oc=(struct oconf *)malloc(6*sizeof(struct oconf));

  nnode=org->nnode;
  loff=(*(pe->nods+i))->loff;

  cnode1=**(pe->nods+0);
  cnode2=**(pe->nods+1);

  dx=(cnode2.d[0]-cnode1.d[0])*0.5;
  dy=(cnode2.d[1]-cnode1.d[1])*0.5;
  dz=(cnode2.d[2]-cnode1.d[2])*0.5;

  /*ADD NODES*/
  ifind=0;
  for(j=0;j<nnode;j++) /*OVERLAPPED NODE*/
  {
    if(cnode1.d[0]+dx>(org->nodes+j)->d[0]-0.001 &&
       cnode1.d[0]+dx<(org->nodes+j)->d[0]+0.001 &&
       cnode1.d[1]+dy>(org->nodes+j)->d[1]-0.001 &&
       cnode1.d[1]+dy<(org->nodes+j)->d[1]+0.001 &&
       cnode1.d[2]+dz>(org->nodes+j)->d[2]-0.001 &&
       cnode1.d[2]+dz<(org->nodes+j)->d[2]+0.001)
    {
      *(gelem.nods+0)=org->nodes+j;
      *(gelem.nods+1)=org->nodes+cnode2.loff;
      *(pe->nods+1)=org->nodes+j;
      ifind=1;
      break;
    }
  }
  if(!ifind) /*NEW NODE*/
  {
    cnode1.code=(org->nodes+nnode-1)->code+1;
    cnode1.loff=nnode;

#if 0
    for(j=0;j<6;j++) /*COPY CONFS*/
    {
      *(oc+j)=*(org->confs+6*loff+j);
    }
#endif
#if 0
 /*CONFS*/  //modified by UJIOKA
    (oc+0)->iconf=0;     (oc+0)->value=0.0;
    (oc+1)->iconf=1;     (oc+1)->value=0.0;
    (oc+2)->iconf=0;     (oc+2)->value=0.0;
    (oc+3)->iconf=1;     (oc+3)->value=0.0;
    (oc+4)->iconf=0;     (oc+4)->value=0.0;
    (oc+5)->iconf=1;     (oc+5)->value=0.0;
#endif
#if 1
 /*CONFS*/  //modified by UJIOKA
    (oc+0)->iconf=0;     (oc+0)->value=0.0;
    (oc+1)->iconf=0;     (oc+1)->value=0.0;
    (oc+2)->iconf=0;     (oc+2)->value=0.0;
    (oc+3)->iconf=0;     (oc+3)->value=0.0;
    (oc+4)->iconf=0;     (oc+4)->value=0.0;
    (oc+5)->iconf=0;     (oc+5)->value=0.0;
#endif
    *(gelem.nods+0)=addnode(cnode1,oc,org);
    *(gelem.nods+1)=org->nodes+cnode2.loff;
    *(pe->nods+1)=*(gelem.nods+0);

    /*MOVE NODES*/
    (*(gelem.nods+i))->d[0]+=dx;
    (*(gelem.nods+i))->d[1]+=dy;
    (*(gelem.nods+i))->d[2]+=dz;
  }

  for(j=0;j<6;j++) /*COPY BONDS*/
  {
    *(gelem.bonds+6*0+j)=0;  /*UJIOKA*/
    *(gelem.bonds+6*1+j)=*(pe->bonds+6*1+j);
    *(pe->bonds+6*1+j)=0;
  }

  free(oc);
/*
sprintf(str,"Elem %d\n",gelem.code);
for(i=0;i<gelem.nnod;i++)
{
  sprintf(s,"Node %d {%.3f,%.3f,%.3f}\n",
          (*(gelem.nods+i))->code,
          (*(gelem.nods+i))->d[0],
          (*(gelem.nods+i))->d[1],
          (*(gelem.nods+i))->d[2]);
  strcat(str,s);
}
MessageBox(NULL,str,"AddElement",MB_OK);
*/
  for(i=0;i<(pe->nban);i++) /*CORRECT POINTER*/
  {
    for(j=0;j<(pe->bans+i)->nnod;j++)
    {
      for(k=0;k<(pe->nnod);k++)
      {
        if(*((pe->bans+i)->nods+j)==*(pe->nods+k))
        {
          *((gelem.bans+i)->nods+j)=*(gelem.nods+k);
        }
      }
    }
  }

  gelem.code=(org->elems+(org->nelem-1))->code+1;
  gelem.loff=org->nelem;

  addelement(&gelem,org);

  gelem.nnod=0;
  gelem.nban=0;

  gelem.nods=NULL;
  gelem.bonds=NULL;
  gelem.bans=NULL;

  return 1;
}/*divideatmidpoint*/

struct oelem **divideallelem(struct organ *org) //by ujioka
/*DIVIDE ALL ELEMENTS ON MID POINT.*/
{
  FILE *fdata;
  int i,nelem;
  struct oelem *pe;

  nelem=org->nelem;
  fdata=fgetstofopenII(DIRECTORY,"w","chains.txt");

  for(i=0;i<nelem;i++)
  {
    pe=org->elems+i;
    divideatmidpoint(pe,org);

    currentpivot(i+1,nelem);

    fprintf(fdata,"CHAIN :NELEM 2\n{\nELEM %5d\nELEM %5d\n}\n\n",
                 (org->elems+i)->code,(org->elems+(org->nelem)-1)->code);
  }
  return multielem;
}/*dividemultielem*/

struct oelem *deleteelement(long int elemoffset,
                            struct organ *orgbefore,
                            long int pointednodecode)
/*DELETE ELEMENT OF ORGAN.*/
{
  char str[80];
  int i,j,jj,count;
  long int nelem,nnod;
  struct oelem *elems;
  struct onode *nods;

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
    *(elems+i)=*(orgbefore->elems+i);
  }
  for(i=elemoffset;i<nelem;i++) /*PUSH UP*/
  {
    *(elems+i)=*(orgbefore->elems+i+1);
    ((elems+i)->loff)--;
  }

  nnod=(orgbefore->elems+elemoffset)->nnod;
  nods=(struct onode *)malloc(nnod*sizeof(struct onode));
  for(i=0;i<nnod;i++)
  {
    *(nods+i)=**((orgbefore->elems+elemoffset)->nods+i);
  }

  /*free(orgbefore->elems);*/
  orgbefore->elems=elems;

  for(i=0;i<nnod;i++) /*DELETE RELEASED NODE.*/
  {
    count=0;

    for(j=0;j<nelem;j++)
    {
      for(jj=0;jj<((elems+j)->nnod);jj++)
      {
        if((*((elems+j)->nods+jj))->code==(nods+i)->code) count++;
      }
    }

    if(count==0 && (nods+i)->code!=pointednodecode)
    {
      sprintf(str,"Delete Node:%ld",(nods+i)->code);
      MessageBox(NULL,str,"DeleteElem",MB_OK);
      deletenode((nods+i)->code,orgbefore);
    }
  }
  free(nods);

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

  return addnode(nadd,NULL,org);
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
  double eps=1.0E-08;
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
      /*MessageBox(NULL,"Intersection Failed.",
                   "FindLastNode",MB_OK);*/
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
						   &distance,NULL))
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

void chasenode(struct onode *node,
			   struct viewparam vp,
			   struct windowparams wp)
/*DRAW ADDING,MOVING NODE BEFORE DEFINITION.*/
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
	SetTextColor(wp.hdcC,RGB(150,150,150));
	drawglobalnode(wp.hdcC,vp,*node,NULL);

	/*BitBlt(wp.hdcB,x,y,maxX,maxY,wp.hdcC,x,y,
		   dwrop);*/ /*GROUND BLACK:SRCPAINT WHITE:SRCAND.*/

	hdc = GetDC(wp.hwnd);
	BitBlt(hdc,x,y,maxX,maxY,wp.hdcB,x,y,SRCCOPY);
	BitBlt(hdc,x,y,maxX,maxY,wp.hdcC,x,y,SRCPAINT);
	ReleaseDC(wp.hwnd,hdc);

	PatBlt(wp.hdcC,x,y,maxX,maxY,PATCOPY);
  }

  return;
}/*chasenode*/

void chaseelement(struct oelem *elem,
                  struct viewparam vp,
                  struct windowparams wp)
/*DRAW ADDING,MOVING ELEMENT BEFORE DEFINITION.*/
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
	drawelement(wp.hdcC,vp,*(elem),ONSCREEN,1);

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

void chaserect(int il,int ir,int it,int ib,
               struct windowparams wp)
/*DRAW SELECTING RECT BEFORE DEFINITION.*/
{
  HWND hoya;
  HDC hdc;
  HPEN hpen,ppen;
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

    hpen=CreatePen(PS_SOLID,1,RGB(150,150,150)); /*GRAY*/
	ppen = (HPEN)SelectObject(wp.hdcC,hpen);

    /*DRAW RECT*/
    MoveToEx(wp.hdcC,il,it,NULL);
    LineTo(wp.hdcC,il,ib);
    LineTo(wp.hdcC,ir,ib);
    LineTo(wp.hdcC,ir,it);
    LineTo(wp.hdcC,il,it);

    SelectObject(wp.hdcC,ppen);
    DeleteObject(hpen);

	/*UPDATE WINDOW*/
    hdc = GetDC(wp.hwnd);
    BitBlt(hdc,x,y,maxX,maxY,wp.hdcB,x,y,SRCCOPY);
    BitBlt(hdc,x,y,maxX,maxY,wp.hdcC,x,y,SRCPAINT);
    ReleaseDC(wp.hwnd,hdc);

    PatBlt(wp.hdcC,x,y,maxX,maxY,PATCOPY);
  }

  return;
}/*chaserect*/

struct oelem **selectmultielem(struct organ *org,
                               struct viewparam vp,
                               int il,int ir,int it,int ib)
/*SELECT MULTI ELEMENTS BY RECT.*/
{
  int i,j,in;

  // select mode by fukushima
  char offflag;
  char mode;
  if(il<=ir) mode=0;
  else mode=1;
  //

  if(nmultielem==0)                      //araki
  {
    multielem=(struct oelem **)malloc(1*sizeof(struct oelem *));
  }

  for(i=0;i<(org->nelem);i++)
  {
    in=1;
    // select mode by fukushima
    offflag=0;
    //

    if((org->elems+i)->sect->dflag==1 &&
       vp.vflag.ev.etype[(org->elems+i)->type]==1)
    {
      for(j=0;j<((org->elems+i)->nnod);j++)
      {
        if(!insiderange(**((org->elems+i)->nods+j),vp.range))
        {
          in=0;
          break;
        }
        if(!nodeinsiderect(**((org->elems+i)->nods+j),il,ir,it,ib,vp))
        {
          if(mode==0)
          {
            in=0; /*OUTSIDE*/
            break;
          }
          else // mode 1
          {
            offflag++;
          }
        }
      }
      // select mode by fukushima
      if(offflag==(org->elems+i)->nnod) in=0;
      //
    }
    else in=0; /*FLAGS OFF*/

    /*IF ALL DOTS OF ELEMENT INSIDE AND FLAGS ON*/
    if(in)
    {
      in=0;
      for(j=0;j<nmultielem;j++)
      {
        if(*(multielem+j)==org->elems+i) in=1; /*LAPPED*/
      }

      if(!in) /*ADD IF NOT LAPPED*/
      {
        nmultielem++;
        multielem=(struct oelem **)
                  realloc(multielem,
                          nmultielem*sizeof(struct oelem *));
        *(multielem+nmultielem-1)=org->elems+i;
      }
    }
  }

  return multielem;
}/*selectmultielem*/

struct owire **selectmultiwire(struct arclmframe *af,
                               struct viewparam vp,
                               int il,int ir,int it,int ib)
/*SELECT MULTI WIRES BY RECT.*/
{
  int i,j,in;

  // select mode by fukushima
  char offflag;
  char mode;
  char str[80];
  if(il<=ir) mode=0;
  else mode=1;
  //

  if(nmultiwire==0)                      //araki
  {
    multiwire=(struct owire **)malloc(1*sizeof(struct owire *));
  }

  for(i=0;i<(af->nelem);i++)
  {
    in=1;
    // select mode by fukushima
    offflag=0;
    //

    if((af->elems+i)->sect->dflag==1 /*&&
       vp.vflag.ev.etype[(af->elems+i)->type]==1*/)
    {
      for(j=0;j<2;j++)
      {
        if(!insiderange((*(af->elems+i)->node[j]),vp.range))
        {
          in=0;
          break;
        }
        if(!nodeinsiderect(*((af->elems+i)->node[j]),il,ir,it,ib,vp))
        {
          if(mode==0)
          {
            in=0; /*OUTSIDE*/
            break;
          }
          else // mode 1
          {
            offflag++;
          }
        }
      }
      // select mode by fukushima
      if(offflag==2) in=0;
      //
    }
    else in=0; /*FLAGS OFF*/

    /*IF ALL DOTS OF ELEMENT INSIDE AND FLAGS ON*/
    if(in)
    {
      in=0;
      for(j=0;j<nmultiwire;j++)
      {
        if(*(multiwire+j)==af->elems+i) in=1; /*LAPPED*/
      }

      if(!in) /*ADD IF NOT LAPPED*/
      {
        nmultiwire++;
        multiwire=(struct owire **)
                  realloc(multiwire,
                          nmultiwire*sizeof(struct owire *));
        *(multiwire+nmultiwire-1)=af->elems+i;
      }
    }
  }

  return multiwire;

}/*selectmultiwire*/

struct owire **selectmultiwireII(struct arclmframe *af,int offset,int *nmultiwire)
/*SELECT MULTI WIRES BY TEXT.*/
{
  int /*nmultiwire,*/nchains;
  FILE *fin;
  int i,ii,j,k,m,n;
  char **data, str[256];
  int code,nstr;

  fin=fgetstofopenII(DIRECTORY,"r","chains.txt");
//  fin=fopen("chains.txt","r");

  m=0;
  while(1)
  {
    data=fgetsbrk(fin,&nstr);
    if(data==NULL) /*END OF TEXT DATA*/
    {
      *nmultiwire=0;
      return NULL;
    }
    if(strncmp(*data,"CHAIN",5)==0)
    {
      if(m==offset) break;
      m++;
    }
    free(data);
  }

  n=strtol(*(data+2),NULL,10);
  data = fgetsbrk(fin, &nstr);
  /*
  sprintf(str,"n=%ld",n);
  errormessage(str);
  */
  multiwire=(struct owire **)malloc(n*sizeof(struct owire *));

  for(i=0;i<n;i++)
  {
    data = fgetsbrk(fin, &nstr);
    code=strtol(*(data+1),NULL,10);
    for(j=0;j<(af->nelem);j++)
    {
      if(((af->elems+j)->code)==code) break;
    }

    *(multiwire+i)=af->elems+j;
  }

  *nmultiwire=n;
  return multiwire;

}/*selectmultiwireII*/

/*** araki for multinode*******************************************************/
struct onode **selectmultinode(struct organ *org,
                               struct viewparam vp,
                               int il,int ir,int it,int ib)
/*SELECT MULTI NODES BY RECT.*/
{
  int i,j,in;

  if(nmultinode==0)
  {
    multinode=(struct onode **)malloc(1*sizeof(struct onode *));
  }

  for(i=0;i<(org->nnode);i++)
  {
    in=1;
/*    if((org->elems+i)->sect->dflag==1 &&
       vp.vflag.ev.etype[(org->elems+i)->type]==1)
    {                                                  */
        if(!insiderange(*(org->nodes+i),vp.range) ||
           !nodeinsiderect(*(org->nodes+i),il,ir,it,ib,vp))
        {
          in=0; /*OUTSIDE*/
        }
/*    }
    else in=0;*/ /*FLAGS OFF*/

    /*IF NODE INSIDE*/
    if(in)
    {
      in=0;
      for(j=0;j<nmultinode;j++)
      {
        if(*(multinode+j)==org->nodes+i) in=1; /*LAPPED*/
      }

      if(!in) /*ADD IF NOT LAPPED*/
      {
        nmultinode++;
        multinode=(struct onode **)
                  realloc(multinode,
                          nmultinode*sizeof(struct onode *));
        *(multinode+nmultinode-1)=org->nodes+i;
      }
    }
  }

  return multinode;
}/*selectmultinode*/
/******************************************************************************/

/*** araki for multinode*******************************************************/
struct onode **deletemultinode(struct organ *org)
/*DELETE MULTI NODES.*/
{
  int i,j,k,count,find;
  long int nelem1,nelem2,nnode1,nnode2,loff;
/*  long int *olist;     */
  struct oelem *elems;
  struct onode *nodes,**plist;
  struct oconf *confs;

  nelem1=org->nelem;

  count=0;
  for(i=0;i<nelem1;i++)
  {
    find=0;
/*    for(j=0;j<nnode2;j++)  */
    for(j=0;j<nmultinode;j++)
    {
      for(k=0;k<((org->elems+i)->nnod);k++)
      {
        if((*((org->elems+i)->nods+k))->loff==(*(multinode+j))->loff)
        {
          find=1;
          break;
        }
      }
      if(find) break;
    }

    if(find) count++;
  }

  nelem2=nelem1-count;
  org->nelem=nelem2;

  elems=(struct oelem *)malloc(nelem2*sizeof(struct oelem));

  if(elems==NULL)
  {
    MessageBox(NULL,"Buffer Null.","DeleteElem",MB_OK);
    return multinode;
  }
  org->nelem=nelem2;

  count=0;
  for(i=0;i<nelem1;i++)
  {
    find=0;
/*    for(j=0;j<nnode2;j++)  */
    for(j=0;j<nmultinode;j++)
    {
      for(k=0;k<((org->elems+i)->nnod);k++)
      {
        if((*((org->elems+i)->nods+k))->loff==(*(multinode+j))->loff)
        {
          find=1;
          break;
        }
      }
      if(find) break;
    }

    if(!find)
    {
      *(elems+count)=*(org->elems+i); /*COPY ELEMS*/
      (elems+count)->loff=count;

      count++;
      if(count>nelem2) break; /*ERROR*/
    }
  }


  /*MODIFY NODES,CONFS*/
  nnode1=org->nnode;

  count=0;
  for(i=0;i<nnode1;i++)
  {
    find=0;
    for(j=0;j<nelem2;j++)
    {
      for(k=0;k<((elems+j)->nnod);k++)
      {
        if((*((elems+j)->nods+k))->loff==i)
        {
          find=1;
          break;
        }
      }
      if(find) break;
    }

    if(!find) count++;
  }

  nnode2=nnode1-count;
  org->nnode=nnode2;

  plist=(struct onode **)malloc(nnode1*sizeof(struct onode *));
  nodes=(struct onode *)malloc(nnode2*sizeof(struct onode));
  confs=(struct oconf *)malloc(6*nnode2*sizeof(struct oconf));

  count=0;
  for(i=0;i<nnode1;i++)
  {
    find=0;
    for(j=0;j<nelem2;j++)
    {
      for(k=0;k<((elems+j)->nnod);k++)
      {
        if((*((elems+j)->nods+k))->loff==i)
        {
          find=1;
          break;
        }
      }
      if(find) break;
    }

    if(find)
    {
      *(nodes+count)=*(org->nodes+i); /*COPY NODES*/
      (nodes+count)->loff=count;
      for(j=0;j<6;j++)
      {
        *(confs+6*count+j)=*(org->confs+6*i+j); /*COPY CONFS*/
      }
      *(plist+i)=nodes+count;         /*POINTER LIST*/

      count++;
      if(count>nnode2) break; /*ERROR*/
    }
    else *(plist+i)=NULL;
    currentpivot((long int)count,nnode2);
  }

  /*SET NODE POINTERS OF ELEMENT*/
  for(i=0;i<nelem2;i++)
  {
/*    for(j=0;j>((elems+i)->nnod);j++)  */
    for(j=0;j<((elems+i)->nnod);j++)
    {
      loff=(*((elems+i)->nods+j))->loff;
      *((elems+i)->nods+j)=*(plist+loff);
    }
  }

  org->elems=elems;
  org->nodes=nodes;
  org->confs=confs;

  return multinode;
}/*deletemultinode*/
/******************************************************************************/

/*** 09.01.29 araki for movemultinode *****************************************/
struct onode **movemultinode(struct organ *org,
                             double dx,double dy,double dz)
/*MOVE MULTI NODE BY INCREMENT.*/
{
  int i,j,k,flag;
  struct onode *pn;

  for(i=0;i<(org->nnode);i++)
  {
    for(j=0;j<nmultinode;j++)
    {
      if((org->nodes+i)->code==(*(multinode+j))->code)
      {
        (org->nodes+i)->d[GX]+=dx;
        (org->nodes+i)->d[GY]+=dy;
        (org->nodes+i)->d[GZ]+=dz;
      }
    }
  }

  return multinode;
}/*movemultinode*/
/******************************************************************************/

void changeconfmultinode(struct organ *org)
/*CHANGE MULTINODE CONFINEMENT OF FRAME.*/ //temporary
{
  int i,j,k,flag;
  struct onode *pn;
  struct oconf conf[6];

  for(i=0;i<(org->nnode);i++)
  {
//    for(j=0;j<nmultinode;j++)
    {
      if((org->nodes+i)->code==(*(multinode+0))->code)
      {
          for(k=0;k<=5;k++)
          {
            conf[k]=*(org->confs+6*((org->nodes+i)->loff)+k);
          }
      }
    }
  }

  for(i=0;i<(org->nnode);i++)
  {
    for(j=1;j<nmultinode;j++)
    {
      if((org->nodes+i)->code==(*(multinode+j))->code)
      {
          for(k=0;k<=5;k++)
          {
            *(org->confs+6*((org->nodes+i)->loff)+k)
            =conf[k];
          }
      }
    }
  }

  return;
}/*changeconfmultinode*/

struct oelem **copymultielem(struct organ *org,
                             double dx,double dy,double dz)
/*COPY MULTI ELEMENTS BY INCREMENT.*/
{
  int i;
  struct oelem *pe;

  for(i=0;i<nmultielem;i++)
  {
    pe=*(multielem+i);
    addelementwithnode(pe,org,dx,dy,dz); /*UNDER CONSTRUCTION*/

    currentpivot(i+1,nmultielem);
  }

  return multielem;
}/*copymultielem*/

struct oelem **movemultielem(struct organ *org,
                             double dx,double dy,double dz)
/*MOVE MULTI ELEMENTS BY INCREMENT.*/
{
  int i,j,k,flag;
  struct oelem *pe;

  for(i=0;i<org->nnode;i++)
  {
    flag=0;
    for(j=0;j<nmultielem;j++)
    {
      pe=*(multielem+j);

      for(k=0;k<pe->nnod;k++)
      {
        if((org->nodes+i)->code==(*(pe->nods+k))->code) flag=1;
      }
      if(flag) break;
    }

    if(flag)
    {
      (org->nodes+i)->d[GX]+=dx;
      (org->nodes+i)->d[GY]+=dy;
      (org->nodes+i)->d[GZ]+=dz;
    }
  }

  return multielem;
}/*movemultielem*/

struct oelem **deletemultielem(struct organ *org)
/*DELETE MULTI ELEMENTS.*/
{
  int i,j,k,count,find;
  long int nelem1,nelem2,nnode1,nnode2,loff;
  long int *olist;
  struct oelem *elems;
  struct onode *nodes,**plist;
  struct oconf *confs;

  nelem1=org->nelem;

  nelem2=nelem1-nmultielem;
  elems=(struct oelem *)malloc(nelem2*sizeof(struct oelem));
  if(elems==NULL)
  {
    MessageBox(NULL,"Buffer Null.","DeleteElem",MB_OK);
    return multielem;
  }
  org->nelem=nelem2;

  olist=(long int *)malloc(nelem2*sizeof(long int)); /*OFFSET LIST*/

  count=0;
  for(i=0;i<nelem1;i++) /*LIST UP OFFSET*/
  {
    find=0;
    for(j=0;j<nmultielem;j++)
    {
      if((*(multielem+j))->loff==i)
      {
        find=1;
        break;
      }
    }
    if(!find)
    {
      *(olist+count)=(long int)i;
      count++;
      if(count>nelem2) break; /*ERROR*/
    }
  }

  for(i=0;i<nelem2;i++) /*PUSH UP*/
  {
    *(elems+i)=*(org->elems+(*(olist+i)));
    ((elems+i)->loff)=i;
  }

  free(olist);

  /*MODIFY NODES,CONFS*/
  nnode1=org->nnode;

  count=0;
  for(i=0;i<nnode1;i++)
  {
    find=0;
    for(j=0;j<nelem2;j++)
    {
      for(k=0;k<((elems+j)->nnod);k++)
      {
        if((*((elems+j)->nods+k))->loff==i)
        {
          find=1;
          break;
        }
      }
      if(find) break;
    }

    if(!find) count++;
  }

  nnode2=nnode1-count;
  org->nnode=nnode2;

  plist=(struct onode **)malloc(nnode1*sizeof(struct onode *));
  nodes=(struct onode *)malloc(nnode2*sizeof(struct onode));
  confs=(struct oconf *)malloc(6*nnode2*sizeof(struct oconf));

  count=0;
  for(i=0;i<nnode1;i++)
  {
    find=0;
    for(j=0;j<nelem2;j++)
    {
      for(k=0;k<((elems+j)->nnod);k++)
      {
        if((*((elems+j)->nods+k))->loff==i)
        {
          find=1;
          break;
        }
      }
      if(find) break;
    }

    if(find)
    {
      *(nodes+count)=*(org->nodes+i); /*COPY NODES*/
      (nodes+count)->loff=count;
      for(j=0;j<6;j++)
      {
        *(confs+6*count+j)=*(org->confs+6*i+j); /*COPY CONFS*/
      }
      *(plist+i)=nodes+count;         /*POINTER LIST*/

      count++;
      if(count>nnode2) break; /*ERROR*/
    }
    else *(plist+i)=NULL;
  }

  /*SET NODE POINTERS OF ELEMENT*/
  for(i=0;i<nelem2;i++)
  {
    for(j=0;j>((elems+i)->nnod);j++)
    {
      loff=(*((elems+i)->nods+j))->loff;
      *((elems+i)->nods+j)=*(plist+loff);
    }
  }

  org->elems=elems;
  org->nodes=nodes;
  org->confs=confs;

  return multielem;
}/*deletemultielem*/

struct oelem **dividemultielem(struct organ *org) //by fukushima
/*DIVIDE MULTI ELEMENTS ON MID POINT.*/
{
  int i;
  struct oelem *pe;

  for(i=0;i<nmultielem;i++)
  {
    pe=*(multielem+i);
    divideatmidpoint(pe,org);

    currentpivot(i+1,nmultielem);
  }

  return multielem;
}/*dividemultielem*/

void chasecurve(struct curve *cv,
                struct viewparam vp,
                struct windowparams wp)
/*DRAW ADDING CURVE.*/
{
  HWND hoya;
  HPEN hpen1,hpen2,hpen3,ppen;
  HBRUSH hbrsh1,hbrsh2,pbrsh;
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

    hpen1 = CreatePen(PS_SOLID,1,RGB(150,150,150));
    hpen2 = CreatePen(PS_SOLID,1,RGB(120,120,120));
    hpen3 = CreatePen(PS_SOLID,1,RGB(50,50,50));
    hbrsh1 = CreateSolidBrush(RGB(120,120,120));
    hbrsh2 = CreateSolidBrush(RGB(50,50,50));
    ppen = (HPEN)SelectObject(wp.hdcC,hpen1);
	pbrsh = (HBRUSH)SelectObject(wp.hdcC,hbrsh1);

    if(cv->type==CTYPE_LINE)
    {
      drawlinecurveext(wp.hdcC,vp,cv,
                       hpen1,hpen2,hbrsh1,NULL);
    }
    else if(cv->type==CTYPE_CIRCLE)
    {
      drawglobalchord(wp.hdcC,vp,*cv,
                      hpen1,hpen2,hpen3,hbrsh1,hbrsh2);
    }

    SelectObject(wp.hdcC,ppen);
    SelectObject(wp.hdcC,pbrsh);
    DeleteObject(hpen1);
    DeleteObject(hpen2);
    DeleteObject(hpen3);
    DeleteObject(hbrsh1);
    DeleteObject(hbrsh2);

    hdc = GetDC(wp.hwnd);
    BitBlt(hdc,x,y,maxX,maxY,wp.hdcB,x,y,SRCCOPY);
    BitBlt(hdc,x,y,maxX,maxY,wp.hdcC,x,y,SRCPAINT);
    ReleaseDC(wp.hwnd,hdc);

    PatBlt(wp.hdcC,x,y,maxX,maxY,PATCOPY);
  }
  return;
}/*chasecurve*/

void chasepolycurve(struct polycurve *pcv,
                    struct viewparam vp,
                    struct windowparams wp)
/*DRAW ADDING POLYCURVE.*/
{
  HWND hoya;
  HPEN hpen,ppen;
  HDC hdc;
  POINT pc,pp;
  int i;
  long int x,y,maxX,maxY;
  struct onode node;

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

    hpen=CreatePen(PS_SOLID,1,RGB(100,100,100));
    ppen=(HPEN)SelectObject(wp.hdcC,hpen);

    for(i=0;i<(pcv->ncurve);i++)
    {
      /*DRAW LINES.*/
      if((pcv->curves+i)->type==CTYPE_LINE)
      {
        drawlinecurveext(wp.hdcC,vp,
                         (pcv->curves+i),
                         hpen,NULL,NULL,NULL);
      }
      /*DRAW ARCS.*/
      if((pcv->curves+i)->type==CTYPE_CIRCLE)
      {
        /*"node" WILL NOT USED FOR HORROW ARC.*/
        drawglobalchordext(wp.hdcC,vp,
                           *(pcv->curves+i),
                           hpen,NULL,NULL,NULL,0,node);
      }
    }

    SelectObject(wp.hdcC,ppen);
    DeleteObject(hpen);

    hdc = GetDC(wp.hwnd);
    BitBlt(hdc,x,y,maxX,maxY,wp.hdcB,x,y,SRCCOPY);
    BitBlt(hdc,x,y,maxX,maxY,wp.hdcC,x,y,SRCPAINT);
    ReleaseDC(wp.hwnd,hdc);

    PatBlt(wp.hdcC,x,y,maxX,maxY,PATCOPY);
  }
  return;
}/*chasepolycurve*/

void chasepolypolycurve(struct polypolycurve *pcv,
                        struct viewparam vp,
                        struct windowparams wp)
/*DRAW ADDING POLYCURVE.*/
{
  HWND hoya;
  HPEN hpen,ppen;
  HDC hdc;
  POINT pc,pp;
  int i,ii;
  long int x,y,maxX,maxY;
  struct onode node;

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

    hpen=CreatePen(PS_SOLID,1,RGB(100,100,100));
    ppen=(HPEN)SelectObject(wp.hdcC,hpen);

    for(ii=0;ii<(pcv->npcurve);ii++)
    {
      for(i=0;i<((pcv->pcurves+ii)->ncurve);i++)
      {
        /*DRAW LINES.*/
        if(((pcv->pcurves+ii)->curves+i)->type==CTYPE_LINE)
        {
          drawlinecurveext(wp.hdcC,vp,
                           ((pcv->pcurves+ii)->curves+i),
                           hpen,NULL,NULL,NULL);
        }
        /*DRAW ARCS.*/
        if(((pcv->pcurves+ii)->curves+i)->type==CTYPE_CIRCLE)
        {
          drawglobalchordext(wp.hdcC,vp,
                             *((pcv->pcurves+ii)->curves+i),
                             hpen,NULL,NULL,NULL,0,node);
        }
      }
    }

    SelectObject(wp.hdcC,ppen);
    DeleteObject(hpen);

    hdc = GetDC(wp.hwnd);
    BitBlt(hdc,x,y,maxX,maxY,wp.hdcB,x,y,SRCCOPY);
    BitBlt(hdc,x,y,maxX,maxY,wp.hdcC,x,y,SRCPAINT);
    ReleaseDC(wp.hwnd,hdc);

    PatBlt(wp.hdcC,x,y,maxX,maxY,PATCOPY);
  }
  return;
}/*chasepolypolycurve*/

int linefeatures(double x1,double y1,double x2,double y2,
                 struct features *f,
                 struct onode origin)
/*FEATURES OF LINE.FOR ONLY Ox<x1,Oy<y1,Ox<x2,Oy<y2.*/
{
  /*char s[80],str[400];*/
  double ox,oy;

  ox=origin.d[GX];
  oy=origin.d[GY];

  x1-=ox;
  y1-=oy;
  x2-=ox;
  y2-=oy;

  if(x1<0.0 || x2<0.0 || y1<0.0 || y2<0.0) return 0;

  f->Ax=0.5*(y1+y2)*(x1-x2);
  f->Ay=0.5*(x1+x2)*(y2-y1);

  if(f->Ax!=0.0)
  {
    f->Gx.d[GX]=(x1-x2)*((x1+x2)*(y1+y2)+(x1*y1+x2*y2))/6.0/(f->Ax);
    f->Gx.d[GY]=(x1-x2)*(y1*y1+y1*y2+y2*y2)/6.0/(f->Ax);
    f->Sx=(f->Ax)*(f->Gx.d[GY]);
  }
  else
  {
    f->Gx.d[GX]=0.0;
    f->Gx.d[GY]=0.0;
    f->Sx=0.0;
  }

  if(f->Ay!=0.0)
  {
    f->Gy.d[GX]=(y2-y1)*(x1*x1+x1*x2+x2*x2)/6.0/(f->Ay);
    f->Gy.d[GY]=(y2-y1)*((x1+x2)*(y1+y2)+(x1*y1+x2*y2))/6.0/(f->Ay);
    f->Sy=(f->Ay)*(f->Gy.d[GX]);
  }
  else
  {
    f->Gy.d[GX]=0.0;
    f->Gy.d[GY]=0.0;
    f->Sy=0.0;
  }

  if(y1>y2) f->Ixx=(x1-x2)*(y1*y1*y1+y1*y1*y2
                           +y1*y2*y2+y2*y2*y2)/12.0;
  else      f->Ixx=(x1-x2)*(y2*y2*y2+y2*y2*y1
                           +y2*y1*y1+y1*y1*y1)/12.0;
  if(x1>x2) f->Iyy=(y2-y1)*(x1*x1*x1+x1*x1*x2
                           +x1*x2*x2+x2*x2*x2)/12.0;
  else      f->Iyy=(y2-y1)*(x2*x2*x2+x2*x2*x1
                           +x2*x1*x1+x1*x1*x1)/12.0;
  /*
  sprintf(str,"\0");
  sprintf(s,"LINE:(%.3f %.3f)(%.3f %.3f)\n",x1,y1,x2,y2);
  strcat(str,s);
  sprintf(s,"Ax=%.3f Ay=%.3f\n",f->Ax, f->Ay);
  strcat(str,s);
  sprintf(s,"Ixx=%.3f Iyy=%.3f\n",f->Ixx,f->Iyy);
  strcat(str,s);
  MessageBox(NULL,str,"One Line Features",MB_OK);
  */
  return 1;
}/*linefeatures*/

int circlefeatures(struct curve *cv,
                   struct features *f,
                   struct onode origin)
/*FEATURES OF CIRCLE.*/
{
  char /*s[80],*/str[400];
  int i;
  double cx,cy,ox,oy,x1,y1,x2,y2,a0,a1,a2,r;
  double pi2,bound;
  double bx[8],by[8];
  double eps=1.0E-08,epsa=1.0E-04;
  struct polycurve dpc={0,0}; /*DIVIDED POLY CURVE.*/
  struct features fi,fj;
  struct onode omin,d[2],d0;

  f->Ax=0.0; f->Ay=0.0;
  f->Sx=0.0; f->Sy=0.0;
  f->Gx.d[GX]=0.0;
  f->Gx.d[GY]=0.0;
  f->Gy.d[GX]=0.0;
  f->Gy.d[GY]=0.0;
  f->Ixx=0.0; f->Iyy=0.0;

  a1=cv->angle[0];
  a2=cv->angle[1];

  if(a1==a2) return 1;

  if((fabs(a2-a1)>2.0*PI+epsa) ||
     (a1<-2.0*PI-epsa || 2.0*PI+epsa<a1) ||
     (a2<-2.0*PI-epsa || 2.0*PI+epsa<a2))
  {
    sprintf(str,"Out of Range.a1=%.3f a2=%.3f",a1,a2);
    MessageBox(NULL,str,"Chord",MB_OK);
    return 0;
  }

  cx=cv->center->d[GX];
  cy=cv->center->d[GY];

  d[0]=*(cv->dots[0]);
  d[1]=*(cv->dots[1]);

  r=cv->radius[0];

  fi.Ax=0.0; fi.Ay=0.0;
  fi.Sx=0.0; fi.Sy=0.0;
  fi.Gx.d[GX]=0.0;
  fi.Gx.d[GY]=0.0;
  fi.Gy.d[GX]=0.0;
  fi.Gy.d[GY]=0.0;
  fi.Ixx=0.0; fi.Iyy=0.0;

  /*FEATURES AROUND CENTER.*/
  pi2=0.5*PI;

  /*FOR -2PI<a1<0 -2PI<a2<2PI*/
  if(a2<a1) /*REVERSED.*/
  {
    a0=a1; a1=a2; a2=a0;
    if(a1>=0.0)
    {
      a1-=2.0*PI;
      a2-=2.0*PI;
    }

    d0=d[0]; d[0]=d[1]; d[1]=d0;
  }

  dpc.curves=(struct curve *)malloc(6*sizeof(struct curve));
  for(i=0;i<6;i++)
  {
    (dpc.curves+i)->center=(struct onode *)
                           malloc(sizeof(struct onode));
    (dpc.curves+i)->dots[0]=(struct onode *)
                            malloc(sizeof(struct onode));
    (dpc.curves+i)->dots[1]=(struct onode *)
                            malloc(sizeof(struct onode));

    (dpc.curves+i)->radius[0]=cv->radius[0];
    *((dpc.curves+i)->center)=*(cv->center);
  }
  /*FIRST DOT.*/
  (dpc.curves+0)->angle[0]=a1;
  *((dpc.curves+0)->dots[0])=d[0];
  /*DIVIDE BY AXIS.*/
  bx[0]=r;   bx[1]=0.0; bx[2]=-r;  bx[3]=0.0;
  by[0]=0.0; by[1]=r;   by[2]=0.0; by[3]=-r;
  bx[4]=r;   bx[5]=0.0; bx[6]=-r;  bx[7]=0.0;
  by[4]=0.0; by[5]=r;   by[6]=0.0; by[7]=-r;

  if(a2>a1)
  {
    for(i=3;i>=0;i--)
    {
      bound=(double)(-i*0.5*PI);

      if((bound-0.5*PI)<=a1 && a1<bound)
      {
        if(a2>bound)
        {
          dpc.ncurve++;

          (dpc.curves+0)->angle[1]=bound;
          (dpc.curves+0)->dots[1]->d[GX]=bx[4-i]+cx;
          (dpc.curves+0)->dots[1]->d[GY]=by[4-i]+cy;

          (dpc.curves+1)->angle[0]=bound;
          (dpc.curves+1)->dots[0]->d[GX]=bx[4-i]+cx;
          (dpc.curves+1)->dots[0]->d[GY]=by[4-i]+cy;
        }
        if(a2>bound+0.5*PI)
        {
          dpc.ncurve++;

          (dpc.curves+1)->angle[1]=bound+0.5*PI;
          (dpc.curves+1)->dots[1]->d[GX]=bx[5-i]+cx;
          (dpc.curves+1)->dots[1]->d[GY]=by[5-i]+cy;

          (dpc.curves+2)->angle[0]=bound+0.5*PI;
          (dpc.curves+2)->dots[0]->d[GX]=bx[5-i]+cx;
          (dpc.curves+2)->dots[0]->d[GY]=by[5-i]+cy;
        }
        if(a2>bound+PI)
        {
          dpc.ncurve++;

          (dpc.curves+2)->angle[1]=bound+PI;
          (dpc.curves+2)->dots[1]->d[GX]=bx[6-i]+cx;
          (dpc.curves+2)->dots[1]->d[GY]=by[6-i]+cy;

          (dpc.curves+3)->angle[0]=bound+PI;
          (dpc.curves+3)->dots[0]->d[GX]=bx[6-i]+cx;
          (dpc.curves+3)->dots[0]->d[GY]=by[6-i]+cy;
        }
        if(a2>bound+1.5*PI)
        {
          dpc.ncurve++;

          (dpc.curves+3)->angle[1]=bound+1.5*PI;
          (dpc.curves+3)->dots[1]->d[GX]=bx[7-i]+cx;
          (dpc.curves+3)->dots[1]->d[GY]=by[7-i]+cy;

          (dpc.curves+4)->angle[0]=bound+1.5*PI;
          (dpc.curves+4)->dots[0]->d[GX]=bx[7-i]+cx;
          (dpc.curves+4)->dots[0]->d[GY]=by[7-i]+cy;
        }

        /*END DOT.*/
        dpc.ncurve++;
        (dpc.curves+dpc.ncurve-1)->angle[1]=a2;
        *((dpc.curves+dpc.ncurve-1)->dots[1])=d[1];
      }
    }
  }

  omin.d[GX]=((dpc.curves+0)->dots[0]->d[GX])-cx;
  omin.d[GY]=((dpc.curves+0)->dots[0]->d[GY])-cy;

  /*CHORD TOTAL.*/
  for(i=0;i<dpc.ncurve;i++)
  {
    x1=((dpc.curves+i)->dots[0]->d[GX])-cx;
    y1=((dpc.curves+i)->dots[0]->d[GY])-cy;
    x2=((dpc.curves+i)->dots[1]->d[GX])-cx;
    y2=((dpc.curves+i)->dots[1]->d[GY])-cy;

    if(omin.d[GX]>x1) omin.d[GX]=x1;
    if(omin.d[GX]>x2) omin.d[GX]=x2;
    if(omin.d[GY]>y1) omin.d[GY]=y1;
    if(omin.d[GY]>y2) omin.d[GY]=y2;

    a1=(dpc.curves+i)->angle[0];
    a2=(dpc.curves+i)->angle[1];

    fi.Ax=0.5*r*r*(a2-a1+0.5*(sin(2.0*a2)-sin(2.0*a1)))
         -0.5*(x1+x2)*(y2-y1);
    fi.Ay=0.5*r*r*(a2-a1+0.5*(sin(2.0*(pi2-a1))-sin(2.0*(pi2-a2))))
         -0.5*(y1+y2)*(x1-x2);
    if(fabs((fi.Ax)-(fi.Ay))>eps)
    {
      sprintf(str,"A Missmatch.Eps=%.3E",(fi.Ax-fi.Ay));
      MessageBox(NULL,str,"Chord",MB_OK);
      /*
      sprintf(str,"\0");
      sprintf(s,"A1      =%.15f Pi\n",a1/PI);       strcat(str,s);
      sprintf(s,"0.5Pi-A1=%.15f Pi\n",(pi2-a1)/PI); strcat(str,s);
      sprintf(s,"A2      =%.15f Pi\n",a2/PI);       strcat(str,s);
      sprintf(s,"0.5Pi-A2=%.15f Pi\n",(pi2-a2)/PI); strcat(str,s);
      MessageBox(NULL,str,"Chord",MB_OK);
      */
    }

    fi.Sx=r*r*r*(cos(a1)*cos(a1)*cos(a1)
                -cos(a2)*cos(a2)*cos(a2))/3.0
         -(x2-x1)*(y1*y1+y1*y2+y2*y2)/3.0
         -(x1*y2-x2*y1)*(y1+y2)/2.0;
    fi.Sy=r*r*r*(cos(pi2-a2)*cos(pi2-a2)*cos(pi2-a2)
				-cos(pi2-a1)*cos(pi2-a1)*cos(pi2-a1))/3.0
         -(y1-y2)*(x1*x1 + x1*x2 + x2*x2)/3.0
         -(y2*x1-y1*x2)*(x1+x2)/2.0;

    fi.Ixx=r*r*r*r*(a2-a1-(sin(4.0*a2)-sin(4.0*a1))/4.0)/8.0
          -(x2-x1)*(y1*y1*y1 + y1*y1*y2 + y1*y2*y2 + y2*y2*y2)/4.0
          -(x1*y2-x2*y1)*(y1*y1 + y1*y2 + y2*y2)/3.0;
    fi.Iyy=r*r*r*r
          *(a2-a1-(sin(4.0*(pi2-a1))-sin(4.0*(pi2-a2)))/4.0)/8.0
          -(y1-y2)*(x1*x1*x1 + x1*x1*x2 + x1*x2*x2 + x2*x2*x2)/4.0
          -(y2*x1-y1*x2)*(x1*x1 + x1*x2 + x2*x2)/3.0;

    f->Ax+=fi.Ax;
    f->Ay+=fi.Ay;
    f->Sx+=fi.Sx;
    f->Sy+=fi.Sy;
    f->Ixx+=fi.Ixx;
    f->Iyy+=fi.Iyy;
  }
  /*
  sprintf(str,"\0");
  sprintf(s,"ANGLES:%.3f %.3f\n",cv->angle[0],cv->angle[1]);
  strcat(str,s);
  sprintf(s,"Ax=%.3f Ay=%.3f\n",f->Ax,f->Ay);
  strcat(str,s);
  sprintf(s,"Ixx=%.3f Iyy=%.3f\n",f->Ixx,f->Iyy);
  strcat(str,s);
  sprintf(s,"G=(%.3f %.3f)\n",(f->Sy)/(f->Ay),(f->Sx)/(f->Ax));
  strcat(str,s);
  MessageBox(NULL,str,"Chords Features",MB_OK);
  */

  if(dpc.ncurve>=2)
  {
    dpc.ncurve++;
    (dpc.curves+dpc.ncurve-1)->type=CTYPE_LINE;
    *((dpc.curves+dpc.ncurve-1)->dots[0])=d[1];
    *((dpc.curves+dpc.ncurve-1)->dots[1])=d[0];

    fi.Ax=0.0; fi.Ay=0.0;
    fi.Sx=0.0; fi.Sy=0.0;
    fi.Ixx=0.0; fi.Iyy=0.0;

    for(i=0;i<dpc.ncurve;i++)
    {
      if(linefeatures(((dpc.curves+i)->dots[0]->d[GX]-cx),
                      ((dpc.curves+i)->dots[0]->d[GY]-cy),
                      ((dpc.curves+i)->dots[1]->d[GX]-cx),
                      ((dpc.curves+i)->dots[1]->d[GY]-cy),
                      &fj,omin))
      {
        fi.Ax+=fj.Ax;
        fi.Ay+=fj.Ay;
        fi.Sx+=(fj.Ax)*(fj.Gx.d[GY]);
        fi.Sy+=(fj.Ay)*(fj.Gy.d[GX]);
        fi.Ixx+=fj.Ixx;
        fi.Iyy+=fj.Iyy;
      }
    }
    if(fi.Ax!=0.0)
    {
      fi.Gx.d[GX]=(fi.Sy)/(fi.Ax);
      fi.Gx.d[GY]=(fi.Sx)/(fi.Ax);
    }
    else
    {
      fi.Gx.d[GX]=0.0;
      fi.Gx.d[GY]=0.0;
    }

    if(fi.Ay!=0.0)
    {
      fi.Gy.d[GX]=(fi.Sy)/(fi.Ay);
      fi.Gy.d[GY]=(fi.Sx)/(fi.Ay);
    }
    else
    {
      fi.Gy.d[GX]=0.0;
      fi.Gy.d[GY]=0.0;
    }

    /*
    sprintf(str,"\0");
    sprintf(s,"Center=(%.3f %.3f)\n",cx,cy);
    strcat(str,s);
    sprintf(s,"Min=(%.3f %.3f)\n",omin.d[GX]+cx,omin.d[GY]+cy);
    strcat(str,s);
    sprintf(s,"Ax=%.3f Ay=%.3f\n",fi.Ax,fi.Ay);
    strcat(str,s);
    sprintf(s,"Ixx=%.3f Iyy=%.3f\n",fi.Ixx,fi.Iyy);
    strcat(str,s);
    sprintf(s,"G=(%.3f %.3f)\n",fi.Gy.d[GX],fi.Gx.d[GY]);
    strcat(str,s);
    MessageBox(NULL,str,"Lines Features on Min",MB_OK);
    */

    /*FEATURES AROUND CENTER.*/
    fi.Ixx+=-(fi.Ax)*(fi.Gx.d[GY])*(fi.Gx.d[GY])
            +(fi.Ax)*(/*cy*/-fi.Gx.d[GY]-omin.d[GY])
                    *(/*cy*/-fi.Gx.d[GY]-omin.d[GY]);
    fi.Iyy+=-(fi.Ay)*(fi.Gy.d[GX])*(fi.Gy.d[GX])
            +(fi.Ay)*(/*cx*/-fi.Gy.d[GX]-omin.d[GX])
                    *(/*cx*/-fi.Gy.d[GX]-omin.d[GX]);

    fi.Sx+=(fi.Ax)*(omin.d[GY]);
    fi.Sy+=(fi.Ay)*(omin.d[GX]);

    /*
    sprintf(str,"\0");
    sprintf(s,"Ax=%.3f Ay=%.3f\n",fi.Ax,fi.Ay);
    strcat(str,s);
    sprintf(s,"Ixx=%.3f Iyy=%.3f\n",fi.Ixx,fi.Iyy);
    strcat(str,s);
    sprintf(s,"G=(%.3f %.3f)\n",(fi.Sy)/(fi.Ay),(fi.Sx)/(fi.Ax));
    strcat(str,s);
    MessageBox(NULL,str,"Lines Features on Center",MB_OK);
    */

    f->Ax+=fi.Ax;
    f->Ay+=fi.Ay;
    f->Sx+=fi.Sx;
    f->Sy+=fi.Sy;
    f->Ixx+=fi.Ixx;
    f->Iyy+=fi.Iyy;
  }
  if(fabs((f->Ax)-(f->Ay))>eps)
  {
    sprintf(str,"A Missmatch.Eps=%.3E",(f->Ax)-(f->Ay));
    MessageBox(NULL,str,"Circle",MB_OK);
  }

  if(f->Ax!=0.0)
  {
    f->Gx.d[GX]=(f->Sy)/(f->Ax);
    f->Gx.d[GY]=(f->Sx)/(f->Ax);
  }
  else
  {
    f->Gx.d[GX]=0.0;
    f->Gx.d[GY]=0.0;
  }
  if(f->Ay!=0.0)
  {
    f->Gy.d[GX]=(f->Sy)/(f->Ay);
    f->Gy.d[GY]=(f->Sx)/(f->Ay);
  }
  else
  {
    f->Gy.d[GX]=0.0;
    f->Gy.d[GY]=0.0;
  }

  /*
  sprintf(str,"\0");
  sprintf(s,"Ax=%.3f Ay=%.3f\n",f->Ax,f->Ay);
  strcat(str,s);
  sprintf(s,"Ixx=%.3f Iyy=%.3f\n",f->Ixx,f->Iyy);
  strcat(str,s);
  sprintf(s,"Gx=(%.3f %.3f)\n",f->Gx.d[GX],f->Gx.d[GY]);
  strcat(str,s);
  sprintf(s,"Gy=(%.3f %.3f)\n",f->Gy.d[GX],f->Gy.d[GY]);
  strcat(str,s);
  MessageBox(NULL,str,"Circle Features on Center",MB_OK);
  */

  /*FEATURES INTO AROUND ORIGIN.*/
  ox=origin.d[GX];
  oy=origin.d[GY];

  f->Ixx+=-(f->Ax)*(f->Gx.d[GY])*(f->Gx.d[GY])
          +(f->Ax)*(cy+f->Gx.d[GY]-oy)*(cy+f->Gx.d[GY]-oy);
  f->Iyy+=-(f->Ay)*(f->Gy.d[GX])*(f->Gy.d[GX])
          +(f->Ay)*(cx+f->Gy.d[GX]-ox)*(cx+f->Gy.d[GX]-ox);

  f->Gx.d[GX]+=(cx-ox);
  f->Gx.d[GY]+=(cy-oy);
  f->Gy.d[GX]+=(cx-ox);
  f->Gy.d[GY]+=(cy-oy);

  f->Sx=(f->Ax)*(f->Gx.d[GY]);
  f->Sy=(f->Ay)*(f->Gy.d[GX]);

  if(cv->hugo==-1)
  {
    f->Ax  *= -1.0;
    f->Ay  *= -1.0;
    f->Sx  *= -1.0;
    f->Sy  *= -1.0;
    f->Ixx *= -1.0;
    f->Iyy *= -1.0;
  }

  /*
  sprintf(str,"\0");
  sprintf(s,"ORIGIN=(%.3f %.3f)\n",ox,oy);
  strcat(str,s);
  sprintf(s,"CENTER=(%.3f %.3f)\n",cx,cy);
  strcat(str,s);
  sprintf(s,"RADIUS=%.3f\n",cv->radius[0]);
  strcat(str,s);
  sprintf(s,"ENDS:(%.3f %.3f)(%.3f %.3f)\n",
          (cv->dots[0]->d[GX]),(cv->dots[0]->d[GY]),
          (cv->dots[1]->d[GX]),(cv->dots[1]->d[GY]));
  strcat(str,s);
  sprintf(s,"ANGLES:%.3f %.3f\n",cv->angle[0],cv->angle[1]);
  strcat(str,s);
  sprintf(s,"Ax=%.3f Ay=%.3f\n",f->Ax, f->Ay);
  strcat(str,s);
  sprintf(s,"Ixx=%.3f Iyy=%.3f\n",f->Ixx,f->Iyy);
  strcat(str,s);
  sprintf(s,"Gx=(%.3f %.3f)\n",f->Gx.d[GX],f->Gx.d[GY]);
  strcat(str,s);
  sprintf(s,"Gy=(%.3f %.3f)\n",f->Gy.d[GX],f->Gy.d[GY]);
  strcat(str,s);
  MessageBox(NULL,str,"Circle Features on Origin",MB_OK);
  */

  return 1;
}/*circlefeatures*/

int curvefeatures(struct curve *cv,
                  struct features *f,
                  struct onode origin)
/*FEATURES OF CURVE AROUND ORIGIN.*/
{
  /*char s[80],str[400];*/
  struct features fi;

  if(cv->type==CTYPE_LINE || cv->type==CTYPE_CIRCLE)
  {
    if(!linefeatures((cv->dots[0]->d[GX]),
                     (cv->dots[0]->d[GY]),
                     (cv->dots[1]->d[GX]),
                     (cv->dots[1]->d[GY]),
                     f,origin)) return 0;
  }
  else return 0;

  if(cv->type==CTYPE_CIRCLE)
  {
    if(!circlefeatures(cv,&fi,origin)) return 0;

    f->Ax+=fi.Ax;
    f->Ay+=fi.Ay;
    f->Sx+=fi.Sx;
    f->Sy+=fi.Sy;
    f->Ixx+=fi.Ixx;
    f->Iyy+=fi.Iyy;

    /*f->Gx.d[GX]=(f->Sy)/(f->Ax);
    f->Gx.d[GY]=(f->Sx)/(f->Ax);
    f->Gy.d[GX]=(f->Sy)/(f->Ay);
    f->Gy.d[GY]=(f->Sx)/(f->Ay);*/
  }

  /*
  sprintf(str,"\0");
  sprintf(s,"LINE:(%.3f %.3f)(%.3f %.3f)\n",x1,y1,x2,y2);
  strcat(str,s);
  sprintf(s,"Ax=%.3f Ay=%.3f\n",f->Ax, f->Ay);
  strcat(str,s);
  sprintf(s,"Ixx=%.3f Iyy=%.3f\n",f->Ixx,f->Iyy);
  strcat(str,s);
  sprintf(s,"Gx=(%.3f %.3f)\n",f->Gx.d[GX],f->Gx.d[GY]);
  strcat(str,s);
  sprintf(s,"Gy=(%.3f %.3f)\n",f->Gy.d[GX],f->Gy.d[GY]);
  strcat(str,s);
  MessageBox(NULL,str,"Features",MB_OK);
  */
  return 1;
}/*curvefeatures*/

struct features polycurvefeatures(HWND hwnd,
                                  struct viewparam *vp,
                                  struct polypolycurve *pc)
/*FEATURES OF POLYCURVE.*/
/*A:SECTIONAL AREA*/
/*I:GEOMETRICAL MOMENT OF INTERIA*/
/*Z:SECTION MODULUS*/
/*i:RADIUS OF GYRATION*/
/*G:CENTER OF GRAVITY*/
{
  HDC hdc;
  HPEN hpen,ppen;
  char s[80],str[400];
  int i,ii;
  double xmax,ymax,xmin,ymin,x1,y1,x2,y2;
  double E,poi,G,Gi,factE,factG,gx,gy,factK;
  double eps=1.0E-08;
  struct onode origin,n1,n2;
  struct features f,ftotal;

  xmax=((pc->pcurves+0)->curves+0)->dots[0]->d[GX];
  ymax=((pc->pcurves+0)->curves+0)->dots[0]->d[GY];
  xmin=xmax;
  ymin=ymax;

  for(ii=0;ii<pc->npcurve;ii++)
  {
    for(i=0;i<(pc->pcurves+ii)->ncurve;i++)
    {
      if(((pc->pcurves+ii)->curves+i)->type==CTYPE_LINE)
      {
        x1=((pc->pcurves+ii)->curves+i)->dots[0]->d[GX];
        y1=((pc->pcurves+ii)->curves+i)->dots[0]->d[GY];
        x2=((pc->pcurves+ii)->curves+i)->dots[1]->d[GX];
        y2=((pc->pcurves+ii)->curves+i)->dots[1]->d[GY];
      }
      if(((pc->pcurves+ii)->curves+i)->type==CTYPE_CIRCLE)
      {
        chordrange(((pc->pcurves+ii)->curves+i),&x1,&y1,&x2,&y2);
      }

      if(xmin>x1) xmin=x1;
      if(ymin>y1) ymin=y1;
      if(xmin>x2) xmin=x2;
      if(ymin>y2) ymin=y2;

      if(xmax<x1) xmax=x1;
      if(ymax<y1) ymax=y1;
      if(xmax<x2) xmax=x2;
      if(ymax<y2) ymax=y2;
    }
  }

  origin.d[GX]=xmin;
  origin.d[GY]=ymin;
  origin.d[GZ]=0.0;

  ftotal.E  =0.0;
  ftotal.poi=0.0;

  ftotal.Hx=xmax-xmin;
  ftotal.Hy=ymax-ymin;
  ftotal.Ax=0.0;
  ftotal.Ay=0.0;
  ftotal.Sx=0.0;
  ftotal.Sy=0.0;
  ftotal.Ixx=0.0;
  ftotal.Iyy=0.0;
  ftotal.Jzz=0.0;
  ftotal.Gx.d[GX]=0.0;
  ftotal.Gx.d[GY]=0.0;
  ftotal.Gy.d[GX]=0.0;
  ftotal.Gy.d[GY]=0.0;
  ftotal.hiju=0.0;

  for(ii=0;ii<pc->npcurve;ii++)
  {
    if(ii==0)
    {
      E  =(pc->pcurves+0)->prop.E;
      poi=(pc->pcurves+0)->prop.poi;
      ftotal.E  =E;
      ftotal.poi=poi;

      G=0.5*E/(1.0+poi);
    }

    factE=((pc->pcurves+ii)->prop.E)/E;

    Gi=0.5*((pc->pcurves+ii)->prop.E)
       /(1.0+((pc->pcurves+ii)->prop.poi));
    factG=Gi/G;

    for(i=0;i<(pc->pcurves+ii)->ncurve;i++)
    {
      if(curvefeatures(((pc->pcurves+ii)->curves+i),&f,origin))
      {
        ftotal.Ax+=factE*(f.Ax);
        ftotal.Ay+=factE*(f.Ay);

        ftotal.Sx+=factE*(f.Sx);
        ftotal.Sy+=factE*(f.Sy);

        ftotal.Ixx+=factE*(f.Ixx);
        ftotal.Iyy+=factE*(f.Iyy);

        /*ST.VENANT'S TORTION CONSTANT AS SOLID ELLIPSE.*/
        /*Jzz=4IxxIyy/(Ixx+Iyy)*/
        factK=1.0;
        ftotal.Jzz+=factG*factK*(4.0*f.Ixx*f.Iyy)/(f.Ixx+f.Iyy);
        /*ftotal.Jzz+=factG*0.5*(f.Ixx+f.Iyy);*/

        ftotal.hiju+=(f.Ax)*((pc->pcurves+ii)->prop.hiju);
      }
      else MessageBox(NULL,"Curve Error.","Features",MB_OK);

      /*
      sprintf(str,"\0");
      sprintf(s,"Ax=%.3f Ay=%.3f\n",f.Ax,f.Ay);
      strcat(str,s);
      sprintf(s,"Ixx=%.3f Iyy=%.3f\n",f.Ixx,f.Iyy);
      strcat(str,s);
      MessageBox(NULL,str,"One Curve Features",MB_OK);
      */
    }
  }

  if(fabs(ftotal.Ax-ftotal.Ay)>eps)
  {
    sprintf(str,"A Missmatch.Eps=%.3E",(ftotal.Ax-ftotal.Ay));
    MessageBox(NULL,str,"Total Features",MB_OK);
  }

  /*G FROM MINIMUM DOT.*/
  if(ftotal.Ax!=0.0)
  {
    ftotal.Gx.d[GX]=ftotal.Sy/ftotal.Ax;
    ftotal.Gx.d[GY]=ftotal.Sx/ftotal.Ax;
  }
  if(ftotal.Ay!=0.0)
  {
    ftotal.Gy.d[GX]=ftotal.Sy/ftotal.Ay;
    ftotal.Gy.d[GY]=ftotal.Sx/ftotal.Ay;
  }

  /*I AROUND G.*/
  ftotal.Igx=ftotal.Ixx-ftotal.Ax*ftotal.Gx.d[GY]*ftotal.Gx.d[GY];
  ftotal.Igy=ftotal.Iyy-ftotal.Ay*ftotal.Gy.d[GX]*ftotal.Gy.d[GX];

  /*G FROM GLOBAL ORIGIN.*/
  gx=origin.d[GX]+ftotal.Gy.d[GX];
  gy=origin.d[GY]+ftotal.Gx.d[GY];
  ftotal.Gg.d[GX]=gx;
  ftotal.Gg.d[GY]=gy;

  /*DRAW G.*/
  if(hwnd!=NULL)
  {
    hdc=GetDC(hwnd);
    hpen=CreatePen(PS_SOLID,1,RGB(255,255,255));
    ppen=(HPEN)SelectObject(hdc,hpen);
    setfontformat(hdc,20,5,"Terminal",255,255,255);

    n1.d[GX]=gx-5.0; n1.d[GY]=gy;      n1.d[GZ]=0.0;
    n2.d[GX]=gx+5.0; n2.d[GY]=gy;      n2.d[GZ]=0.0;
    drawgloballine(hdc,*vp,n1,n2);
    n1.d[GX]=gx;     n1.d[GY]=gy-8.0;  n1.d[GZ]=0.0;
    n2.d[GX]=gx;     n2.d[GY]=gy+8.0;  n2.d[GZ]=0.0;
    drawgloballine(hdc,*vp,n1,n2);
    n1.d[GX]=gx+8.0; n1.d[GY]=gy+10.0; n1.d[GZ]=0.0;
    drawglobaltext(hdc,*vp,n1,"G");

    SelectObject(hdc,ppen);
    DeleteObject(hpen);
    ReleaseDC(hwnd,hdc);
  }

  sprintf(str,"\0");
  sprintf(s,"Origin=(%.3f %.3f)\n",xmin,ymin);
  strcat(str,s);
  sprintf(s,"Hx=%.3f Hy=%.3f\n",ftotal.Hx, ftotal.Hy);
  strcat(str,s);
  sprintf(s,"Ax=%.3f Ay=%.3f\n",ftotal.Ax, ftotal.Ay);
  strcat(str,s);
  sprintf(s,"Ixx=%.3f Iyy=%.3f\n",ftotal.Ixx,ftotal.Iyy);
  strcat(str,s);
  sprintf(s,"Igx=%.3f Igy=%.3f\n",ftotal.Igx,ftotal.Igy);
  strcat(str,s);
  sprintf(s,"Gx=(%.3f %.3f)\n",ftotal.Gx.d[GX],ftotal.Gx.d[GY]);
  strcat(str,s);
  sprintf(s,"Gy=(%.3f %.3f)\n",ftotal.Gy.d[GX],ftotal.Gy.d[GY]);
  strcat(str,s);
  sprintf(s,"Gg=(%.3f %.3f)\n",ftotal.Gg.d[GX],ftotal.Gg.d[GY]);
  strcat(str,s);
  MessageBox(NULL,str,"Total Features",MB_OK);

  return ftotal;
}/*polycurvefeatures*/

double ultimatebendingofpolycurve(HWND hwnd,
                                  struct viewparam *vp,
                                  double Nu,
                                  struct polypolycurve *ppc)
/*ULTIMATE BENDING OF POLYCURVE.*/
/*A:SECTIONAL AREA*/
/*G:CENTER OF GRAVITY*/
{
  HDC hdc;
  HPEN hpen,ppen;
  char s[256],str[500];
  int i,ii,j,k;
  double xmax,ymax,xmin,ymin,x1,y1,x2,y2;
  double E,poi,factE,gx,gy;
  double eps1=1.0E-08,eps2=1.0E-04;
  struct onode origin,n1,n2;
  struct features f,ftotal;

  int nbound,ibound,find;
  int nanswer; /*ANSWERS OF CUBIC EQUATION.*/
  double Np,Mp,Yn;
  double *bounds,b1,b2,*Nb,*Mb;
  double C0,C1,C2,C3;
  double answer[3],erate;
  double cx,cy,r,a1,a2,signx,signy;
  struct oprop *op;
  struct curve *cv;
  struct polycurve dpc;
  struct polypolycurve dppc;

  /*DIVIDE CHORDS AND COPY POLYPOLYCURVE*/
  dppc.npcurve=ppc->npcurve; /*NOT CHANGE*/
  dppc.pcurves=(struct polycurve *)
               malloc((ppc->npcurve)*sizeof(struct polycurve));

  for(ii=0;ii<ppc->npcurve;ii++)
  {
    (dppc.pcurves+ii)->prop=(ppc->pcurves+ii)->prop;
    (dppc.pcurves+ii)->loff  =0;
    (dppc.pcurves+ii)->ncurve=0;
    (dppc.pcurves+ii)->type  =0;
    (dppc.pcurves+ii)->curves=NULL;

    for(i=0;i<(ppc->pcurves+ii)->ncurve;i++)
    {
      if(((ppc->pcurves+ii)->curves+i)->type==CTYPE_LINE)
      {
        dpc.loff=0;
        dpc.ncurve=1;
        dpc.curves=(struct curve *)malloc(sizeof(struct curve));
        initializecurve(dpc.curves+0);
        (dpc.curves+0)->dots[0]=(struct onode *)
                                 malloc(sizeof(struct onode));
        (dpc.curves+0)->dots[1]=(struct onode *)
                                 malloc(sizeof(struct onode));

        copycurve((dpc.curves+0),((ppc->pcurves+ii)->curves+i));
        addcurvestopolycurve((dppc.pcurves+ii),&dpc);
      }
      if(((ppc->pcurves+ii)->curves+i)->type==CTYPE_CIRCLE)
      {
        dividechord(&dpc,((ppc->pcurves+ii)->curves+i));
        addcurvestopolycurve((dppc.pcurves+ii),&dpc);
      }
    }
  }

  /*RANGE*/
  xmax=((dppc.pcurves+0)->curves+0)->dots[0]->d[GX];
  ymax=((dppc.pcurves+0)->curves+0)->dots[0]->d[GY];
  xmin=xmax;
  ymin=ymax;

  for(ii=0;ii<dppc.npcurve;ii++)
  {
    for(i=0;i<(dppc.pcurves+ii)->ncurve;i++)
    {
      /*CHORDS DIVIDED.*/
      x1=((dppc.pcurves+ii)->curves+i)->dots[0]->d[GX];
      y1=((dppc.pcurves+ii)->curves+i)->dots[0]->d[GY];
      x2=((dppc.pcurves+ii)->curves+i)->dots[1]->d[GX];
      y2=((dppc.pcurves+ii)->curves+i)->dots[1]->d[GY];

      if(xmin>x1) xmin=x1;
      if(ymin>y1) ymin=y1;
      if(xmin>x2) xmin=x2;
      if(ymin>y2) ymin=y2;

      if(xmax<x1) xmax=x1;
      if(ymax<y1) ymax=y1;
      if(xmax<x2) xmax=x2;
      if(ymax<y2) ymax=y2;
    }
  }

  /*MINIMUM DOT AS ORIGIN*/
  origin.d[GX]=xmin;
  origin.d[GY]=ymin;
  origin.d[GZ]=0.0;

  /*INITIAL*/
  ftotal.E  =0.0;
  ftotal.poi=0.0;

  ftotal.Hx=xmax-xmin;
  ftotal.Hy=ymax-ymin;
  ftotal.Ax=0.0;
  ftotal.Ay=0.0;
  ftotal.Sx=0.0;
  ftotal.Sy=0.0;
  ftotal.Gx.d[GX]=0.0;
  ftotal.Gx.d[GY]=0.0;
  ftotal.Gy.d[GX]=0.0;
  ftotal.Gy.d[GY]=0.0;

  for(ii=0;ii<dppc.npcurve;ii++)
  {
    if(ii==0) /*REPRESENTATIVE FACTOR*/
    {
      E  =(dppc.pcurves+0)->prop.E;
      poi=(dppc.pcurves+0)->prop.poi;
      ftotal.E  =E;
      ftotal.poi=poi;
    }

    factE=((dppc.pcurves+ii)->prop.E)/E;

    for(i=0;i<(dppc.pcurves+ii)->ncurve;i++)
    {
      if(curvefeatures(((dppc.pcurves+ii)->curves+i),&f,origin))
      {
        ftotal.Ax+=factE*(f.Ax);
        ftotal.Ay+=factE*(f.Ay);

        ftotal.Sx+=factE*(f.Sx);
        ftotal.Sy+=factE*(f.Sy);
      }
      else MessageBox(NULL,"Curve Error.","Features",MB_OK);
    }
  }

  if(fabs(ftotal.Ax-ftotal.Ay)>eps1)
  {
    sprintf(str,"A Missmatch.Eps=%.3E",(ftotal.Ax-ftotal.Ay));
    MessageBox(NULL,str,"Total Features",MB_OK);
  }

  /*G FROM MINIMUM DOT.*/
  if(ftotal.Ax!=0.0)
  {
    ftotal.Gx.d[GX]=ftotal.Sy/ftotal.Ax;
    ftotal.Gx.d[GY]=ftotal.Sx/ftotal.Ax;
  }
  if(ftotal.Ay!=0.0)
  {
    ftotal.Gy.d[GX]=ftotal.Sy/ftotal.Ay;
    ftotal.Gy.d[GY]=ftotal.Sx/ftotal.Ay;
  }

  /*G FROM GLOBAL ORIGIN.*/
  gx=origin.d[GX]+ftotal.Gy.d[GX];
  gy=origin.d[GY]+ftotal.Gx.d[GY];
  ftotal.Gg.d[GX]=gx;
  ftotal.Gg.d[GY]=gy;

  /*DRAW G.*/
  if(hwnd!=NULL)
  {
    hdc=GetDC(hwnd);
    hpen=CreatePen(PS_SOLID,1,RGB(255,255,255));
	ppen=(HPEN)SelectObject(hdc,hpen);
	setfontformat(hdc,20,5,"Terminal",255,255,255);

    n1.d[GX]=gx-5.0; n1.d[GY]=gy;      n1.d[GZ]=0.0;
    n2.d[GX]=gx+5.0; n2.d[GY]=gy;      n2.d[GZ]=0.0;
    drawgloballine(hdc,*vp,n1,n2);
    n1.d[GX]=gx;     n1.d[GY]=gy-8.0;  n1.d[GZ]=0.0;
    n2.d[GX]=gx;     n2.d[GY]=gy+8.0;  n2.d[GZ]=0.0;
    drawgloballine(hdc,*vp,n1,n2);
    n1.d[GX]=gx+8.0; n1.d[GY]=gy+10.0; n1.d[GZ]=0.0;
    drawglobaltext(hdc,*vp,n1,"G");

    SelectObject(hdc,ppen);
    DeleteObject(hpen);
    ReleaseDC(hwnd,hdc);
  }

MessageBox(NULL,"BEGIN","Mp",MB_OK);
  /*ULTIMATE BENDING AROUND x AXIS,UPPER COMPRESSED.*/
  /*COUNT AND SET BOUNDARIES*/
  nbound=0;
  bounds=NULL;
  for(ii=0;ii<dppc.npcurve;ii++)
  {
    for(i=0;i<(dppc.pcurves+ii)->ncurve;i++)
    {
      for(j=0;j<=1;j++) /*HEAD,TAIL*/
      {
        n1=*(((dppc.pcurves+ii)->curves+i)->dots[j]);

        find=0;
        for(k=0;k<nbound;k++)
        {
          if(n1.d[GY]>=*(bounds+k)-eps1 &&
             n1.d[GY]<=*(bounds+k)+eps1)
          {
            find=1;
            break;
          }
        }
        if(!find)
        {
          nbound++;
          bounds=(double *)realloc(bounds,nbound*sizeof(double));
          *(bounds+nbound-1)=n1.d[GY]; /*SET BOUNDARIES*/
        }
      }
    }
  }

  /*SORT BOUNDARIES*/
  sortdoublepointer(bounds,nbound);

  /*SET N,M OF EACH BOUNDARY*/
  Nb=(double *)malloc(nbound*sizeof(double));
  Mb=(double *)malloc(nbound*sizeof(double));

  for(i=0;i<nbound;i++)
  {
    boundaryvaluepolypolycurve(&dppc,*(bounds+i),gy,(Nb+i),(Mb+i));
  }

  /*SEARCH RANGE Nu INSIDE*/
  if(Nu<*(Nb+nbound-1) || *(Nb+0)<Nu)
  {
    /*fprintf(fout0,"OUT OF RANGE.\n");*/
    return 0.0;
  }

  for(i=1;i<nbound;i++)
  {
    if(*(Nb+i)<=Nu && Nu<=*(Nb+i-1))
    {
      ibound=i;
      break;
    }
  }

  /*SET UP CUBIC EQUATION*/
  C0=-Nu; C1=0.0; C2=0.0; C3=0.0;

  b1=*(bounds+ibound-1);
  b2=*(bounds+ibound);

  for(ii=0;ii<dppc.npcurve;ii++)
  {
    op=&((dppc.pcurves+ii)->prop);

    for(i=0;i<(dppc.pcurves+ii)->ncurve;i++)
    {
      cv=(dppc.pcurves+ii)->curves+i;

      x1=cv->dots[0]->d[GX];
      y1=cv->dots[0]->d[GY];
      x2=cv->dots[1]->d[GX];
      y2=cv->dots[1]->d[GY];

      if(cv->type==CTYPE_LINE)
      {
        if(x1==x2);
        else if(y1>=b1-eps2 && y2>=b1-eps2) /*ALL UPPER*/
        {
          C0+=(op->fuc)*(-0.5*(y1+y2)*(x1-x2));
          C1+=(op->fuc)*(x1-x2);
        }
        else if(y1<=b2+eps2 && y2<=b2+eps2) /*ALL LOWER*/
        {
          C0+=(op->fut)*(-0.5*(y1+y2)*(x1-x2));
          C1+=(op->fut)*(x1-x2);
        }
        else if(y1>=y2) /*CROSSING AND HEAD UPPER*/
        {
          C0+=(op->fuc)*(-0.5*y1*y1*(x1-x2)/(y1-y2));
          C1+=(op->fuc)*(-y1*(x1-x2)/(y1-y2));
          C2+=(op->fuc)*(0.5*(x1-x2)/(y1-y2));

          C0+=(op->fut)*(-0.5*y2*y2*(x2-x1)/(y1-y2));
          C1+=(op->fut)*(-y2*(x2-x1)/(y1-y2));
          C2+=(op->fut)*(0.5*(x2-x1)/(y1-y2));
        }
        else if(y2>=y1) /*CROSSING AND TAIL UPPER*/
        {
          C0+=(op->fuc)*(-0.5*y2*y2*(x2-x1)/(y1-y2));
          C1+=(op->fuc)*(-y2*(x2-x1)/(y1-y2));
          C2+=(op->fuc)*(0.5*(x2-x1)/(y1-y2));

          C0+=(op->fut)*(-0.5*y1*y1*(x1-x2)/(y1-y2));
          C1+=(op->fut)*(-y1*(x1-x2)/(y1-y2));
          C2+=(op->fut)*(0.5*(x1-x2)/(y1-y2));
        }
      }
      else if(cv->type==CTYPE_CIRCLE)
      {
        cx=cv->center->d[GX];
        cy=cv->center->d[GY];

        a1=cv->angle[0];
        a2=cv->angle[1];

        r=cv->radius[0];

        if(x1<=cx+eps2 && x2<=cx+eps2)      signx=-1.0;
        else if(x1>=cx-eps2 && x2>=cx-eps2) signx=1.0;
        else break; /*ERROR.CHORD NOT DIVIDED.*/

        if(y1<=cy+eps2 && y2<=cy+eps2)      signy=-1.0;
        else if(y1>=cy-eps2 && y2>=cy-eps2) signy=1.0;
        else break; /*ERROR.CHORD NOT DIVIDED.*/

        if(a1==a2);
        else if(r==0.0);
        else if(y1>=b1-eps2 && y2>=b1-eps2) /*ALL UPPER*/
        {
          C0+=(op->fuc)*(-(x2-x1)*cy
                         +(x1-cx)*(y1-cy)
                         -(x2-cx)*(y2-cy)
                         -/*signx**/r*r*(sin(2.0*a1)-sin(2.0*a2)
                                     +2.0*a1-2.0*a2)/4.0);
          C1+=(op->fuc)*(x2-x1);

sprintf(str,"ALL UPPER {Ci}=%.3f, %.3f, %.3f, %.3f\n",C0,C1,C2,C3);
MessageBox(NULL,str,"Ultimate Bending",MB_OK);

        }
        else if(y1<=b2+eps2 && y2<=b2+eps2) /*ALL LOWER*/
        {
          C0+=(op->fut)*(-(x2-x1)*cy
                         +(x1-cx)*(y1-cy)
                         -(x2-cx)*(y2-cy)
                         -/*signx**/r*r*(sin(2.0*a1)-sin(2.0*a2)
                                     +2.0*a1-2.0*a2)/4.0);
          C1+=(op->fut)*(x2-x1);

sprintf(str,"ALL LOWER {Ci}=%.3f, %.3f, %.3f, %.3f\n",C0,C1,C2,C3);
MessageBox(NULL,str,"Ultimate Bending",MB_OK);

        }
        else if(y1>=y2) /*CROSSING AND HEAD UPPER*/
        {
          /*if(-2.0*PI<=a1 && a1<0.0) a1+=2.0*PI;
          if(-2.0*PI<=a2 && a2<0.0) a2+=2.0*PI;*/
          if(-2.0*PI<=a1 && a1<    -PI) a1+=2.0*PI;
          if(     PI< a1 && a1<=2.0*PI) a1-=2.0*PI;
          if(-2.0*PI<=a2 && a2<    -PI) a2+=2.0*PI;
          if(     PI< a2 && a2<=2.0*PI) a2-=2.0*PI;

/********
      Aupper=-x1*yn+x1*y1
             -r*r*(sin(2.0*a1)-sin(2.0*a0)+2.0*a1-2.0*a0)/4.0;
      Alower=x2*yn-x2*y2
             -r*r*(sin(2.0*a0)-sin(2.0*a2)+2.0*a0-2.0*a2)/4.0;
********/
          /*UPPER*/
          C0+=(op->fuc)*((x1-cx)*cy
                         +(x1-cx)*(y1-cy)
                         -(signy)*/*signx**/r*r*(sin(2.0*a1)+2.0*a1)/4.0
                         +/*signx**/PI*(-1.250*r*cy
                                    +signy*0.157*cy*cy
                                    +0.407/r*cy*cy*cy)/4.0);
          C1+=(op->fuc)*(-(x1-cx)
                         +signx*PI*(1.250*r
                                    -signy*0.314*cy
                                    -1.221/r*cy*cy)/4.0);
          C2+=(op->fuc)*(signx*PI*(signy*0.157+1.221/r*cy)/4.0);
          C3+=(op->fuc)*(-signx*PI*0.407/r/4.0);

sprintf(str,"HEAD UPPER 1 {Ci}=%.3f, %.3f, %.3f, %.3f\n",C0,C1,C2,C3);
MessageBox(NULL,str,"Ultimate Bending",MB_OK);

          /*LOWER*/
          C0+=(op->fut)*(-(x2-cx)*cy
                         -(x2-cx)*(y2-cy)
                         +/*signx**/r*r*(sin(2.0*a2)+2.0*a2)/4.0
                         -/*signx**/PI*(-1.250*r*cy
                                    +signy*0.157*cy*cy
                                    +0.407/r*cy*cy*cy)/4.0);
          C1+=(op->fut)*((x2-cx)
                         -signx*PI*(1.250*r
                                    -signy*0.314*cy
                                    -1.221/r*cy*cy)/4.0);
          C2+=(op->fut)*(-signx*PI*(signy*0.157+1.221/r*cy)/4.0);
          C3+=(op->fut)*(signx*PI*0.407/r/4.0);

sprintf(str,"HEAD UPPER {Ci}=%.3f, %.3f, %.3f, %.3f\n",C0,C1,C2,C3);
MessageBox(NULL,str,"Ultimate Bending",MB_OK);

        }
        else if(y2>=y1) /*CROSSING AND TAIL UPPER*/
        {
          /*if(-2.0*PI<=a1 && a1<0.0) a1+=2.0*PI;
          if(-2.0*PI<=a2 && a2<0.0) a2+=2.0*PI;*/
          if(-2.0*PI<=a1 && a1<    -PI) a1+=2.0*PI;
          if(     PI< a1 && a1<=2.0*PI) a1-=2.0*PI;
          if(-2.0*PI<=a2 && a2<    -PI) a2+=2.0*PI;
          if(     PI< a2 && a2<=2.0*PI) a2-=2.0*PI;

          /*UPPER*/
          C0+=(op->fuc)*(-(x2-cx)*cy
                         -(x2-cx)*(y2-cy)
                         +(signy)*/*signx**/r*r*(sin(2.0*a2)+2.0*a2)/4.0
                         -/*signx**/PI*(-1.250*r*cy
                                    +signy*0.157*cy*cy
                                    +0.407/r*cy*cy*cy)/4.0);
          C1+=(op->fuc)*((x2-cx)
                         -signx*PI*(1.250*r
                                    -signy*0.314*cy
                                    -1.221/r*cy*cy)/4.0);
          C2+=(op->fuc)*(-signx*PI*(signy*0.157+1.221/r*cy)/4.0);
          C3+=(op->fuc)*(signx*PI*0.407/r/4.0);

sprintf(str,"TAIL UPPER 1 {Ci}=%.3f, %.3f, %.3f, %.3f\n",C0,C1,C2,C3);
MessageBox(NULL,str,"Ultimate Bending",MB_OK);

          /*LOWER*/
          C0+=(op->fut)*((x1-cx)*cy
                         +(x1-cx)*(y1-cy)
                         -/*signx**/r*r*(sin(2.0*a1)+2.0*a1)/4.0
                         +/*signx**/PI*(-1.250*r*cy
                                    +signy*0.157*cy*cy
                                    +0.407/r*cy*cy*cy)/4.0);
          C1+=(op->fut)*(-(x1-cx)
                         +signx*PI*(1.250*r
                                    -signy*0.314*cy
                                    -1.221/r*cy*cy)/4.0);
          C2+=(op->fut)*(signx*PI*(signy*0.157+1.221/r*cy)/4.0);
          C3+=(op->fut)*(-signx*PI*0.407/r/4.0);

sprintf(str,"TAIL UPPER {Ci}=%.3f, %.3f, %.3f, %.3f\n",C0,C1,C2,C3);
MessageBox(NULL,str,"Ultimate Bending",MB_OK);

        }
      }
    }
  }

sprintf(str,"\0");
sprintf(s,"CUBIC EQUATION\n");
strcat(str,s);
sprintf(s,"(%.3E) + (%.3E)y + (%.3E)y2 + (%.3E)y3 = 0\n",
        C0,C1,C2,C3);
strcat(str,s);
MessageBox(NULL,str,"Ultimate Bending",MB_OK);

  /*NEUTRAL AXIS DETERMINATION*/
  if(C3!=0.0)      nanswer=cubicequation(C3,C2,C1,C0,answer);
  else if(C2!=0.0) nanswer=quadraticequation(C2,C1,C0,answer);
  else if(C1!=0.0)
  {
    nanswer=1;
    answer[0]=-C0/C1;
  }
  else
  {
    /*fprintf(fout0,"NO ANSWER.\n");*/
    return 0.0;
  }

  find=0;
  for(i=0;i<nanswer;i++)
  {
    if(b2-eps2<=answer[i] && answer[i]<=b1+eps2)
    {
      /*fprintf(fout0,"CHOSEN ANSWER=%.5f\n",answer[ii]);*/
      Yn=answer[i];
      find=1;
    }
  }
  if(!find)
  {
    /*fprintf(fout0,"NO ANSWER.\n");*/
    return 0.0;
  }

  /*Np,Mp FROM NEUTRAL AXIS*/
  boundaryvaluepolypolycurve(&dppc,Yn,gy,&Np,&Mp);

  /*CHECK ERROR RATE*/
  if(Nu!=0.0) erate=(Nu-Np)/Nu;
  /*fprintf(fout0,"ERROR RATE=.5f\n",erate);*/

  sprintf(str,"\0");
  sprintf(s,"Origin=(%.3f %.3f)\n",xmin,ymin);
  strcat(str,s);
  sprintf(s,"Hx=%.3f Hy=%.3f\n",ftotal.Hx, ftotal.Hy);
  strcat(str,s);
  sprintf(s,"Ax=%.3f Ay=%.3f\n",ftotal.Ax, ftotal.Ay);
  strcat(str,s);
  sprintf(s,"Gx=(%.3f %.3f)\n",ftotal.Gx.d[GX],ftotal.Gx.d[GY]);
  strcat(str,s);
  sprintf(s,"Gy=(%.3f %.3f)\n",ftotal.Gy.d[GX],ftotal.Gy.d[GY]);
  strcat(str,s);
  sprintf(s,"Gg=(%.3f %.3f)\n",ftotal.Gg.d[GX],ftotal.Gg.d[GY]);
  strcat(str,s);
  sprintf(s,"Nu=%.3f\n",Nu);
  strcat(str,s);
  sprintf(s,"Np=%.3f Mp=%.3f\n",Np,Mp);
  strcat(str,s);
  MessageBox(NULL,str,"Ultimate Bending",MB_OK);

  return Mp;
}/*ultimatebendingofpolycurves*/

int dividechord(struct polycurve *dpc,/*DIVIDED POLY CURVE.*/
                struct curve *ocv)
/*DIVIDE CHORD BY x,y AXIS.*/
{
  int i;
  double r,cx,cy,a0,a1,a2,bx[8],by[8];
  double bound;
  struct onode d[2],d0;

  dpc->loff=0;
  dpc->ncurve=0;

  r =ocv->radius[0];
  cx=ocv->center->d[GX];
  cy=ocv->center->d[GY];
  a1=ocv->angle[0];
  a2=ocv->angle[1];

  d[0]=*(ocv->dots[0]);
  d[1]=*(ocv->dots[1]);

  /*FOR -2PI<a1<0 -2PI<a2<2PI*/
  if(a2<a1) /*REVERSED.*/
  {
    a0=a1; a1=a2; a2=a0;
    if(a1>=0.0)
    {
      a1-=2.0*PI;
      a2-=2.0*PI;
    }

    d0=d[0]; d[0]=d[1]; d[1]=d0;
  }

  dpc->curves=(struct curve *)malloc(6*sizeof(struct curve));
  for(i=0;i<6;i++)
  {
    (dpc->curves+i)->type=ocv->type;
    (dpc->curves+i)->hugo=ocv->hugo;

    (dpc->curves+i)->center =(struct onode *)
                             malloc(sizeof(struct onode));
    (dpc->curves+i)->dots[0]=(struct onode *)
                             malloc(sizeof(struct onode));
    (dpc->curves+i)->dots[1]=(struct onode *)
                             malloc(sizeof(struct onode));
    (dpc->curves+i)->radius[0]=r;
    *((dpc->curves+i)->center)=*(ocv->center);
  }

  /*FIRST DOT.*/
  (dpc->curves+0)->angle[0]=a1;
  *((dpc->curves+0)->dots[0])=d[0];

  /*DIVIDE BY AXIS.*/
  bx[0]=r;   bx[1]=0.0; bx[2]=-r;  bx[3]=0.0;
  by[0]=0.0; by[1]=r;   by[2]=0.0; by[3]=-r;
  bx[4]=r;   bx[5]=0.0; bx[6]=-r;  bx[7]=0.0;
  by[4]=0.0; by[5]=r;   by[6]=0.0; by[7]=-r;

  if(a2>a1)
  {
    for(i=3;i>=0;i--)
    {
      bound=-((double)i)*0.5*PI;

      if((bound-0.5*PI)<=a1 && a1<bound)
      {
        if(a2>bound)
        {
          dpc->ncurve++;

          (dpc->curves+0)->angle[1]=bound;
          (dpc->curves+0)->dots[1]->d[GX]=bx[4-i]+cx;
          (dpc->curves+0)->dots[1]->d[GY]=by[4-i]+cy;

          (dpc->curves+1)->angle[0]=bound;
          (dpc->curves+1)->dots[0]->d[GX]=bx[4-i]+cx;
          (dpc->curves+1)->dots[0]->d[GY]=by[4-i]+cy;
        }
        if(a2>bound+0.5*PI)
        {
          dpc->ncurve++;

          (dpc->curves+1)->angle[1]=bound+0.5*PI;
          (dpc->curves+1)->dots[1]->d[GX]=bx[5-i]+cx;
          (dpc->curves+1)->dots[1]->d[GY]=by[5-i]+cy;

          (dpc->curves+2)->angle[0]=bound+0.5*PI;
          (dpc->curves+2)->dots[0]->d[GX]=bx[5-i]+cx;
          (dpc->curves+2)->dots[0]->d[GY]=by[5-i]+cy;
        }
        if(a2>bound+PI)
        {
          dpc->ncurve++;

          (dpc->curves+2)->angle[1]=bound+PI;
          (dpc->curves+2)->dots[1]->d[GX]=bx[6-i]+cx;
          (dpc->curves+2)->dots[1]->d[GY]=by[6-i]+cy;

          (dpc->curves+3)->angle[0]=bound+PI;
          (dpc->curves+3)->dots[0]->d[GX]=bx[6-i]+cx;
          (dpc->curves+3)->dots[0]->d[GY]=by[6-i]+cy;
        }
        if(a2>bound+1.5*PI)
        {
          dpc->ncurve++;

          (dpc->curves+3)->angle[1]=bound+1.5*PI;
          (dpc->curves+3)->dots[1]->d[GX]=bx[7-i]+cx;
          (dpc->curves+3)->dots[1]->d[GY]=by[7-i]+cy;

          (dpc->curves+4)->angle[0]=bound+1.5*PI;
          (dpc->curves+4)->dots[0]->d[GX]=bx[7-i]+cx;
          (dpc->curves+4)->dots[0]->d[GY]=by[7-i]+cy;
        }

        /*END DOT.*/
        dpc->ncurve++;
        (dpc->curves+dpc->ncurve-1)->angle[1]=a2;
        *((dpc->curves+dpc->ncurve-1)->dots[1])=d[1];
      }
    }
  }

  return dpc->ncurve;
}/*dividechord*/

void boundaryvaluepolypolycurve(struct polypolycurve *ppc,
                                double Yn,double Yg,
                                double *Nb,double *Mb)
{
  int i,j;
  double Ni,Mi;

  *Nb=0.0;
  *Mb=0.0;

  for(i=0;i<(ppc->npcurve);i++)
  {
    for(j=0;j<((ppc->pcurves+i)->ncurve);j++)
    {
      boundaryvaluecurve(((ppc->pcurves+i)->curves+j),
                         Yn,Yg,
                         ((ppc->pcurves+i)->prop),
                         &Ni,&Mi);
      *Nb+=Ni;
      *Mb+=Mi;
    }
  }

  return;
}/*boundaryvaluepolypolycurve*/

void boundaryvaluecurve(struct curve *cv,
                       double Yn,double Yg,
                       struct oprop op,
                       double *Nb,double *Mb)
{
  /*char s[256],str[500];*/
  double eps=1.0E-08;
  double Aupper,Alower,Ygupper,Yglower;
  double Nupper,Nlower,Mupper,Mlower;
  double x1,y1,x2,y2,cx,cy,r,a0,a1,a2,Sx,yn,sign;

  *Nb=0.0;
  *Mb=0.0;

  Aupper=0.0;
  Alower=0.0;
  Ygupper=0.0;
  Yglower=0.0;

  if(cv->type==CTYPE_LINE)
  {
    x1=cv->dots[0]->d[GX];
    y1=cv->dots[0]->d[GY]-Yn;
    x2=cv->dots[1]->d[GX];
    y2=cv->dots[1]->d[GY]-Yn;

    if(x1==x2) return;

    if(y1>=0.0 && y2>=0.0) /*ALL UPPER*/
    {
      Aupper=0.5*(y1+y2)*(x1-x2);
      Ygupper=(x1-x2)*(y1*y1+y1*y2+y2*y2)/6.0/Aupper+Yn;
    }
    else if(y1<=0.0 && y2<=0.0) /*ALL LOWER*/
    {
      Alower=0.5*(y1+y2)*(x1-x2);
      Yglower=(x1-x2)*(y1*y1+y1*y2+y2*y2)/6.0/Alower+Yn;
    }
    else if(y1>=0.0) /*CROSSING AND HEAD UPPER*/
    {
      Aupper=0.5*y1*y1*(x1-x2)/(y1-y2);
      Ygupper=y1/3.0+Yn;

      Alower=0.5*y2*y2*(x2-x1)/(y1-y2);
      Yglower=y2/3.0+Yn;
    }
    else if(y2>=0.0) /*CROSSING AND TAIL UPPER*/
    {
      Aupper=0.5*y2*y2*(x2-x1)/(y1-y2);
      Ygupper=y2/3.0+Yn;

      Alower=0.5*y1*y1*(x1-x2)/(y1-y2);
      Yglower=y1/3.0+Yn;
    }
  }
  if(cv->type==CTYPE_CIRCLE)
  {
    cx=cv->center->d[GX];
    cy=cv->center->d[GY];
    x1=cv->dots[0]->d[GX]-cx;
    y1=cv->dots[0]->d[GY]-cy;
    x2=cv->dots[1]->d[GX]-cx;
    y2=cv->dots[1]->d[GY]-cy;

    a1=cv->angle[0];
    a2=cv->angle[1];

    r=cv->radius[0];

    yn=Yn-cy;

    if(x1<=0.0+eps && x2<=0.0+eps)      sign=1.0;
    else if(x1>=0.0-eps && x2>=0.0-eps) sign=1.0;
    else return; /*ERROR.CHORD NOT DIVIDED.*/

    if(a1==a2) return;
    if(r==0.0) return;

    if(y1>=yn && y2>=yn) /*ALL UPPER*/
    {
      Aupper=(x2-x1)*yn+x1*y1-x2*y2
             -r*r*(sin(2.0*a1)-sin(2.0*a2)+2.0*a1-2.0*a2)/4.0;
      Sx=0.5*(x2-x1)*yn*yn
         +0.5*(x1*y1*y1-x2*y2*y2)
         +sign*r*r*r*(cos(a1)*cos(a1)*cos(a1)
                     -cos(a2)*cos(a2)*cos(a2))/3.0;
      Ygupper=Sx/Aupper+cy;
/*
sprintf(str,"\0");
sprintf(s,"Curve:(%.3f %.3f)(%.3f %.3f)\n",
        x1,y1,x2,y2);
strcat(str,s);
sprintf(s,"All Upper:Yn=%.3f Aupper=%.3f Ygupper=%.3f\n",
        Yn,Aupper,Ygupper);
strcat(str,s);
MessageBox(NULL,str,"BoundaryValue",MB_OK);
*/
    }
    else if(y1<=yn && y2<=yn) /*ALL LOWER*/
    {
      Alower=(x2-x1)*yn+x1*y1-x2*y2
             -r*r*(sin(2.0*a1)-sin(2.0*a2)+2.0*a1-2.0*a2)/4.0;
      Sx=0.5*(x2-x1)*yn*yn
         +0.5*(x1*y1*y1-x2*y2*y2)
         +sign*r*r*r*(cos(a1)*cos(a1)*cos(a1)
                     -cos(a2)*cos(a2)*cos(a2))/3.0;
      Yglower=Sx/Alower+cy;
/*
sprintf(str,"\0");
sprintf(s,"Curve:(%.3f %.3f)(%.3f %.3f)\n",
        x1,y1,x2,y2);
strcat(str,s);
sprintf(s,"All Lower:Yn=%.3f Alower=%.3f Ylower=%.3f\n",
        Yn,Alower,Yglower);
strcat(str,s);
MessageBox(NULL,str,"BoundaryValue",MB_OK);
*/
    }
    else if(y1>=yn) /*CROSSING AND HEAD UPPER*/
    {
      a0=asin(yn/r); /*-0.5PI<=a0<=0.5PI*/
      if(-2.0*PI<=a1 && a1<-1.5*PI) a0=-2.0*PI+a0;
      if(-1.5*PI<=a1 && a1<-0.5*PI) a0=-PI-a0;
      if( 0.5*PI<=a1 && a1< 1.5*PI)  a0=PI-a0;
      if( 1.5*PI<=a1 && a1<=2.0*PI) a0=2.0*PI+a0;

      Aupper=-x1*yn+x1*y1
             -r*r*(sin(2.0*a1)-sin(2.0*a0)+2.0*a1-2.0*a0)/4.0;
      Sx=-0.5*x1*yn*yn
         +0.5*x1*y1*y1
         +sign*r*r*r*(cos(a1)*cos(a1)*cos(a1)
                     -cos(a0)*cos(a0)*cos(a0))/3.0;
      Ygupper=Sx/Aupper+cy;

      Alower=x2*yn-x2*y2
             -r*r*(sin(2.0*a0)-sin(2.0*a2)+2.0*a0-2.0*a2)/4.0;
      Sx=0.5*x2*yn*yn
         -0.5*x2*y2*y2
         +sign*r*r*r*(cos(a0)*cos(a0)*cos(a0)
                     -cos(a2)*cos(a2)*cos(a2))/3.0;
      Yglower=Sx/Alower+cy;
/*
sprintf(str,"\0");
sprintf(s,"Curve:(%.3f %.3f)(%.3f %.3f)\n",
        x1,y1,x2,y2);
strcat(str,s);
sprintf(s,"ANGLES:%.3f < %.3f < %.3f\n",
        a1,a0,a2);
strcat(str,s);
sprintf(s,"Head Upper Yn=%.3f\n",Yn);
strcat(str,s);
sprintf(s,"Aupper=%.3f Ygupper=%.3f\n",
        Aupper,Ygupper);
strcat(str,s);
sprintf(s,"Alower=%.3f Yglower=%.3f\n",
        Alower,Yglower);
strcat(str,s);
MessageBox(NULL,str,"BoundaryValue",MB_OK);
*/
    }
    else if(y2>=yn) /*CROSSING AND TAIL UPPER*/
    {
      a0=asin(yn/r);
      if(-2.0*PI<=a2 && a2<=-1.5*PI) a0=-2.0*PI+a0;
      if(-1.5*PI< a2 && a2<=-0.5*PI) a0=-PI-a0;
      if( 0.5*PI< a2 && a2<= 1.5*PI) a0=PI-a0;
      if( 1.5*PI< a2 && a2<= 2.0*PI) a0=2.0*PI+a0;

      Aupper=x2*yn-x2*y2
             -r*r*(sin(2.0*a0)-sin(2.0*a2)+2.0*a0-2.0*a2)/4.0;
      Sx=0.5*x2*yn*yn
         -0.5*x2*y2*y2
         +sign*r*r*r*(cos(a0)*cos(a0)*cos(a0)
                     -cos(a2)*cos(a2)*cos(a2))/3.0;
      Ygupper=Sx/Aupper+cy;

      Alower=-x1*yn+x1*y1
             -r*r*(sin(2.0*a1)-sin(2.0*a0)+2.0*a1-2.0*a0)/4.0;
      Sx=-0.5*x1*yn*yn
         +0.5*x1*y1*y1
         +sign*r*r*r*(cos(a1)*cos(a1)*cos(a1)
                     -cos(a0)*cos(a0)*cos(a0))/3.0;
      Yglower=Sx/Alower+cy;
/*
sprintf(str,"\0");
sprintf(s,"Curve:(%.3f %.3f)(%.3f %.3f)\n",
        x1,y1,x2,y2);
strcat(str,s);
sprintf(s,"ANGLES:%.3f < %.3f < %.3f\n",
        a1,a0,a2);
strcat(str,s);
sprintf(s,"Tail Upper Yn=%.3f\n",Yn);
strcat(str,s);
sprintf(s,"Aupper=%.3f Ygupper=%.3f\n",
        Aupper,Ygupper);
strcat(str,s);
sprintf(s,"Alower=%.3f Yglower=%.3f\n",
        Alower,Yglower);
strcat(str,s);
MessageBox(NULL,str,"BoundaryValue",MB_OK);
*/
    }
  }

  Nupper=(op.fuc)*Aupper; /*fuc<0 FOR COMPRESSION*/
  Nlower=(op.fut)*Alower; /*fut>0 FOR TENSION*/

  Mupper=Nupper*(Yg-Ygupper);
  Mlower=Nlower*(Yg-Yglower);

  *Nb=Nupper+Nlower;
  *Mb=Mupper+Mlower;

  return;
}/*boundaryvaluecurve*/

struct features sectionfigfeatures(struct osect *os)
/*FEATURES OF SECTION FIGURES.*/
/*A:SECTIONAL AREA*/
/*I:GEOMETRICAL MOMENT OF INTERIA*/
/*J:ST.VENANT'S TORTION FACTOR*/
{
  int i;
  double factE,factG,G,Gi;
  struct features ftotal;

  ftotal.E  =0.0;
  ftotal.poi=0.0;
  ftotal.Ax=0.0;
  ftotal.Ay=0.0;
  ftotal.Ixx=0.0;
  ftotal.Iyy=0.0;
  ftotal.Igx=0.0;
  ftotal.Igy=0.0;
  ftotal.Jzz=0.0;
  ftotal.hiju=0.0;

  if((os->nfig)<=0) return ftotal;

  ftotal.E  =(os->figs+0)->prop->E;
  ftotal.poi=(os->figs+0)->prop->poi;
  G=0.5*(ftotal.E)/(1.0+(ftotal.poi));

  for(i=0;i<(os->nfig);i++)
  {
    if(((os->figs+i)->area)>0.0 && ftotal.E>0.0)
    {
      factE=((os->figs+i)->prop->E)/(ftotal.E);

      Gi=0.5*((os->figs+i)->prop->E)
         /(1.0+((os->figs+i)->prop->poi));
      factG=Gi/G;

      ftotal.Ax +=factE*((os->figs+i)->area);
      ftotal.Ay +=factE*((os->figs+i)->area);
      ftotal.Ixx+=factE*((os->figs+i)->Ixx);
      ftotal.Iyy+=factE*((os->figs+i)->Iyy);
      ftotal.Igx+=factE*((os->figs+i)->Ixx);
      ftotal.Igy+=factE*((os->figs+i)->Iyy);

      ftotal.Jzz+=factG*((os->figs+i)->Jzz);

      ftotal.hiju+=((os->figs+i)->area)
                  *((os->figs+i)->prop->hiju);
    }
  }

  return ftotal;
}/*sectionfigfeatures*/

struct features sectionfeatures(struct osect *os)
/*FEATURES OF SECTION BY FIGURES OR POLYCURVES.*/
/*SET TOTAL FEATURES INTO SECTION "OS".*/
{
  struct features f;
  /*char s[80],str[400];*/

  /*IF NFIG=0 AND NPCURVE=0, HOLD CURRENT FEATURES.*/

  if((os->nfig)>0)
  {
    f=sectionfigfeatures(os);
    os->E   =f.E;
    os->poi =f.poi;
    os->area=f.Ax; /*[m]*/
    os->Ixx =f.Igx;
    os->Iyy =f.Igy;
    os->Jzz =f.Jzz;
    os->hiju[0]=f.hiju;
    os->hiju[1]=f.hiju;
    os->hiju[2]=f.hiju;
  }
  else if((os->ppc.npcurve)>0)
  {
    /*CREATE FEATURES*/
    f=polycurvefeatures(NULL,NULL,&(os->ppc));
/*
 sprintf(str,"\0");
 sprintf(s,"E=%.3f Poi=%.3f Hiju=%.3f\n",f.E,f.poi,f.hiju);
 strcat(str,s);
 sprintf(s,"Ax=%.3f Ay=%.3f\n",f.Ax, f.Ay);
 strcat(str,s);
 sprintf(s,"Igx=%.3f Igy=%.3f Jzz=%.3f\n",f.Igx,f.Igy,f.Jzz);
 strcat(str,s);
 MessageBox(NULL,str,"Section Features",MB_OK);
*/
    f.Ax  /=1.0E+06; /*[mm] Into [m]*/
    f.Igx /=1.0E+12;
    f.Igy /=1.0E+12;
    f.Jzz /=1.0E+12;
    f.hiju/=1.0E+06;

    os->E   =f.E;
    os->poi =f.poi;
    os->area=f.Ax; /*[m]*/
    os->Ixx =f.Igx;
    os->Iyy =f.Igy;
    os->Jzz =f.Jzz;
    os->hiju[0]=f.hiju;
    os->hiju[1]=f.hiju;
    os->hiju[2]=f.hiju;

    /*YIELD SURFACE:UNDER CONSTRUCTION.*/
    /*for(j=0;j<6;j++)
    {
      sect->fmax[j]=1.0;
      sect->fmin[j]=1.0;
    }*/
  }

  return f;
}/*sectionfeatures*/

void initializecurve(struct curve *cv)
{
  cv->loff=0;
  cv->type=0;
  cv->radius[0]=0.0;
  cv->radius[1]=0.0;
  cv->angle[0]=0.0;
  cv->angle[1]=0.0;
  cv->center=NULL;
  cv->dots[0]=NULL;
  cv->dots[1]=NULL;
  cv->dots[2]=NULL;
  cv->tan[0]=NULL;
  cv->tan[1]=NULL;
  cv->tan[2]=NULL;

  return;
}/*initializecurve*/

struct curve *malloccurves(int ncurve)
{
  int i;
  struct curve *c;

  c=(struct curve *)malloc(ncurve*sizeof(struct curve));
  for(i=0;i<ncurve;i++)
  {
    initializecurve(c+i);
  }
  return c;
}/*malloccurves*/

void createcircledata(struct curve *cv)
{
  /*char s[80],str[256];*/
  double r,a,cx,cy;
  double eps=1.0E-08;

  cx=cv->center->d[GX];
  cy=cv->center->d[GY];

  if((cv->radius[0])<=0.0) /*RADIUS DEFINITION*/
  {
    cv->radius[0]=distancedotdot((*(cv->center)),(*(cv->dots[0])));
  }

  if(cv->dots[0]==NULL)
  {
    r=cv->radius[0];
    a=cv->angle[0];

    if((cv->angle[0])<(cv->angle[1])) cv->hugo= 1;
    else                              cv->hugo=-1;

    cv->dots[0]=(struct onode *)malloc(sizeof(struct onode));
    cv->dots[0]->d[GX]=cx+r*cos(a);
    cv->dots[0]->d[GY]=cy+r*sin(a);
    cv->dots[0]->d[GZ]=0.0;
  }
  else
  {
    cv->angle[0]=definecircleangle(cv,0,cv->dots[0]->d[GX],
                                        cv->dots[0]->d[GY]);
    cv->radius[0]=distancedotdot((*(cv->center)),
                                 (*(cv->dots[0])));
  }

  if(cv->dots[1]==NULL)
  {
    r=cv->radius[0];
    a=cv->angle[1];

    cv->dots[1]=(struct onode *)malloc(sizeof(struct onode));
    cv->dots[1]->d[GX]=cx+r*cos(a);
    cv->dots[1]->d[GY]=cy+r*sin(a);
    cv->dots[1]->d[GZ]=0.0;
  }
  else
  {
    cv->angle[1]=definecircleangle(cv,1,cv->dots[1]->d[GX],
                                        cv->dots[1]->d[GY]);

    r=distancedotdot((*(cv->center)),(*(cv->dots[1])));
    if(fabs(r-(cv->radius[0]))>eps)
    {
      cv->dots[1]->d[GX]=cx+(cv->radius[0])*cos(cv->angle[1]);
      cv->dots[1]->d[GY]=cy+(cv->radius[0])*sin(cv->angle[1]);
    }
  }
  /*
  sprintf(str,"\0");
  sprintf(s,"R=%.15f\n",cv->radius[0]); strcat(str,s);
  sprintf(s,"Dot1=(%.5f,%.5f)\n",cv->dots[0]->d[GX],
                                 cv->dots[0]->d[GY]);
  strcat(str,s);
  sprintf(s,"Dot2=(%.5f,%.5f)\n",cv->dots[1]->d[GX],
                                 cv->dots[1]->d[GY]);
  strcat(str,s);
  sprintf(s,"A1=%.15f Pi\n",cv->angle[0]/PI); strcat(str,s);
  sprintf(s,"A2=%.15f Pi\n",cv->angle[1]/PI); strcat(str,s);
  MessageBox(NULL,str,"Create Chord",MB_OK);
  */
  return;
}/*createcircledata*/

void setcurveascircle(struct curve *cv,
                      int loff,int hugo,
                      double *radius,
                      double *angle1,double *angle2,
                      double cx,double cy,
                      double x1,double y1,
                      double x2,double y2)
/*INPUT:CENTER,RADIUS,1st ANGLE,2nd ANGLE*/
/*      CENTER,       1st DOT,  2nd ANGLE*/
/*      CENTER,RADIUS,1st ANGLE,2nd DOT  */
/*      CENTER,       1st DOT,  2nd DOT  */
{
  /*char str[256];*/
  double r;
  double eps=1.0E-08;

  initializecurve(cv);
  cv->loff=loff;
  cv->type=CTYPE_CIRCLE;
  cv->hugo=hugo;

  if(radius!=NULL) cv->radius[0]=(*radius);

  if(angle1!=NULL) cv->angle[0]=(*angle1);
  if(angle2!=NULL) cv->angle[1]=(*angle2);

  cv->center =(struct onode *)malloc(sizeof(struct onode));
  cv->dots[0]=(struct onode *)malloc(sizeof(struct onode));
  cv->dots[1]=(struct onode *)malloc(sizeof(struct onode));

  (*(cv->center))=setcoord(cx,cy,0.0);

  if(angle1==NULL)
  {
    (*(cv->dots[0]))=setcoord(x1,y1,0.0);
    cv->angle[0]=definecircleangle(cv,0,x1,y1);
    cv->radius[0]=distancedotdot((*(cv->center)),
                                 (*(cv->dots[0])));

/*sprintf(str,"1st Dot=(%.3f %.3f)",
        cv->dots[0]->d[GX],cv->dots[0]->d[GY]);
MessageBox(NULL,str,"SetCircle",MB_OK);*/
  }
  else
  {
    cv->dots[0]->d[GX]=cx+(*radius)*cos(*angle1);
    cv->dots[0]->d[GY]=cy+(*radius)*sin(*angle1);
    cv->dots[0]->d[GZ]=0.0;
  }

  if(angle2==NULL)
  {
    cv->angle[1]=definecircleangle(cv,1,x2,y2);

    r=distancedotdot((*(cv->center)),(*(cv->dots[1])));
    if(fabs(r-(cv->radius[0]))>eps)
    {
      x2=cx+(cv->radius[0])*cos(cv->angle[1]);
      y2=cy+(cv->radius[0])*sin(cv->angle[1]);
    }
    (*(cv->dots[1]))=setcoord(x2,y2,0.0);
  }
  else
  {
    cv->dots[1]->d[GX]=cx+(cv->radius[0])*cos(*angle2);
    cv->dots[1]->d[GY]=cy+(cv->radius[0])*sin(*angle2);
    cv->dots[1]->d[GZ]=0.0;
  }

/*sprintf(str,"Circle Dots=(%.3f %.3f)(%.3f %.3f)",
        cv->dots[0]->d[GX],cv->dots[0]->d[GY],
        cv->dots[1]->d[GX],cv->dots[1]->d[GY]);
MessageBox(NULL,str,"SetCircle",MB_OK);*/

  return;
}/*setcurveascircle*/

void setcurveasline(struct curve *cv,
                    int loff,
					double x1,double y1,
					double x2,double y2)
{
  initializecurve(cv);
  cv->loff=loff;
  cv->type=CTYPE_LINE;
  cv->hugo=1;

  cv->dots[0]=(struct onode *)malloc(sizeof(struct onode));
  cv->dots[1]=(struct onode *)malloc(sizeof(struct onode));

  *(cv->dots[0])=setcoord(x1,y1,0.0);
  *(cv->dots[1])=setcoord(x2,y2,0.0);
  return;
}/*setcurveasline*/

void setpolycurveaskatakou(struct polycurve *pc,
                           int type,
                           double H,double B,
                           double tf1,double tf2,double tw,
                           double ri1,double ri2,
                           double ri3,double ri4,
                           double ro1,double ro2,
                           double ro3,double ro4)
{
  int i;
  double a1,a2;

  if(type==CTYPE_FB)
  {
    if(ro1>0.0) pc->ncurve=8;
    else        pc->ncurve=4;

    pc->curves=(struct curve *)malloc(pc->ncurve
                                      *sizeof(struct curve));

    setcurveasline((pc->curves+0),0,-B/2.0+ro1,-H/2.0,
                                     B/2.0-ro1,-H/2.0);
    setcurveasline((pc->curves+1),1, B/2.0,    -H/2.0+ro1,
                                     B/2.0,     H/2.0-ro1);
    setcurveasline((pc->curves+2),2, B/2.0-ro1, H/2.0,
                                    -B/2.0+ro1, H/2.0);
    setcurveasline((pc->curves+3),3,-B/2.0,     H/2.0-ro1,
                                    -B/2.0,    -H/2.0+ro1);
    if(ro1>0.0)
    {
      a1=-0.5*PI; a2=0.0;
      setcurveascircle((pc->curves+4),4,1,
                       &ro1,&a1,&a2,(B/2.0-ro1),(-H/2.0+ro1),
                       0.0,0.0,0.0,0.0);
      a1=-2.0*PI; a2=-1.5*PI;
      setcurveascircle((pc->curves+5),5,1,
                       &ro1,&a1,&a2,(B/2.0-ro1),(H/2.0-ro1),
                       0.0,0.0,0.0,0.0);
      a1=-1.5*PI; a2=-1.0*PI;
      setcurveascircle((pc->curves+6),6,1,
                       &ro1,&a1,&a2,(-B/2.0+ro1),(H/2.0-ro1),
                       0.0,0.0,0.0,0.0);
      a1=-1.0*PI; a2=-0.5*PI;
      setcurveascircle((pc->curves+7),7,1,
                       &ro1,&a1,&a2,(-B/2.0+ro1),(-H/2.0+ro1),
                       0.0,0.0,0.0,0.0);
    }
  }
  else if(type==CTYPE_L)
  {
    if(ri3>0.0) pc->ncurve=7;
    else        pc->ncurve=6;

    if(ro2>tw)  ro2=tw;
    if(ro4>tf2) ro4=tf2;

    if(0.0<ro2 && ro2<tw)  pc->ncurve++;
    if(0.0<ro4 && ro4<tf2) pc->ncurve++;

    pc->curves=(struct curve *)malloc(pc->ncurve
                                      *sizeof(struct curve));

    setcurveasline((pc->curves+0),0,-B/2.0,       -H/2.0,
                                     B/2.0,       -H/2.0);
    setcurveasline((pc->curves+1),1, B/2.0-ro4,   -H/2.0+tf2,
                                    -B/2.0+tw+ri3,-H/2.0+tf2);
    setcurveasline((pc->curves+2),2,-B/2.0+tw,    -H/2.0+tf2+ri3,
                                    -B/2.0+tw,     H/2.0-ro2);
    setcurveasline((pc->curves+3),3,-B/2.0,        H/2.0,
                                    -B/2.0,       -H/2.0);

    i=4;

    if(0.0<=ro4 && ro4<tf2)
    {
      setcurveasline((pc->curves+i),i,B/2.0,-H/2.0,
                                      B/2.0,-H/2.0+tf2-ro4);
      i++;
    }
    if(0.0<ro4 && ro4<=tf2)
    {
      a1=-2.0*PI; a2=-1.5*PI;
      setcurveascircle((pc->curves+i),i,1,
                       &ro4,&a1,&a2,
                       (B/2.0-ro4),(-H/2.0+tf2-ro4),
                       0.0,0.0,0.0,0.0);
      i++;
    }

    if(0.0<=ro2 && ro2<tw)
    {
      setcurveasline((pc->curves+i),i,-B/2.0+tw-ro2,H/2.0,
                                      -B/2.0,       H/2.0);
      i++;
    }
    if(0.0<ro2 && ro2<=tw)
    {
      a1=-2.0*PI; a2=-1.5*PI;
      setcurveascircle((pc->curves+i),i,1,
                       &ro2,&a1,&a2,
                       (-B/2.0+tw-ro2),(H/2.0-ro2),
                       0.0,0.0,0.0,0.0);
      i++;
    }

    if(ri3>0.0)
    {
      a1=-0.5*PI; a2=-1.0*PI;
      setcurveascircle((pc->curves+i),i,-1,
                       &ri3,&a1,&a2,
                       (-B/2.0+tw+ri3),(-H/2.0+tf2+ri3),
                       0.0,0.0,0.0,0.0);
    }
  }
  else if(type==CTYPE_H)
  {
    if(ri1>0.0) pc->ncurve=16;
    else        pc->ncurve=12;

    pc->curves=(struct curve *)malloc(pc->ncurve
                                      *sizeof(struct curve));

    setcurveasline((pc->curves+ 0), 0,-B/2.0,     -H/2.0,
                                       B/2.0,     -H/2.0);
    setcurveasline((pc->curves+ 1), 1, B/2.0,     -H/2.0,
                                       B/2.0,     -H/2.0+tf2);
    setcurveasline((pc->curves+ 2), 2, B/2.0,     -H/2.0+tf2,
                                       tw/2.0+ri1,-H/2.0+tf2);
    setcurveasline((pc->curves+ 3), 3, tw/2.0,    -H/2.0+tf2+ri1,
                                       tw/2.0,     H/2.0-tf1-ri1);
    setcurveasline((pc->curves+ 4), 4, tw/2.0+ri1, H/2.0-tf1,
                                       B/2.0,      H/2.0-tf1);
    setcurveasline((pc->curves+ 5), 5, B/2.0,      H/2.0-tf1,
                                       B/2.0,      H/2.0);
    setcurveasline((pc->curves+ 6), 6, B/2.0,      H/2.0,
                                      -B/2.0,      H/2.0);
    setcurveasline((pc->curves+ 7), 7,-B/2.0,      H/2.0,
                                      -B/2.0,      H/2.0-tf1);
    setcurveasline((pc->curves+ 8), 8,-B/2.0,      H/2.0-tf1,
                                      -tw/2.0-ri1, H/2.0-tf1);
    setcurveasline((pc->curves+ 9), 9,-tw/2.0,     H/2.0-tf1-ri1,
                                      -tw/2.0,    -H/2.0+tf2+ri1);
    setcurveasline((pc->curves+10),10,-tw/2.0-ri1,-H/2.0+tf2,
                                      -B/2.0,     -H/2.0+tf2);
    setcurveasline((pc->curves+11),11,-B/2.0,     -H/2.0+tf2,
                                      -B/2.0,     -H/2.0);

    if(ri1>0.0)
    {
      a1=-0.5*PI; a2=-1.0*PI;
      setcurveascircle((pc->curves+12),12,-1,
                       &ri1,&a1,&a2,(tw/2.0+ri1),(-H/2.0+tf2+ri1),
                       0.0,0.0,0.0,0.0);
      a1=-1.0*PI; a2=-1.5*PI;
      setcurveascircle((pc->curves+13),13,-1,
                       &ri1,&a1,&a2,(tw/2.0+ri1),(H/2.0-tf1-ri1),
                       0.0,0.0,0.0,0.0);
      a1=-1.5*PI; a2=-2.0*PI;
      setcurveascircle((pc->curves+14),14,-1,
                       &ri1,&a1,&a2,(-tw/2.0-ri1),(H/2.0-tf1-ri1),
                       0.0,0.0,0.0,0.0);
      a1= 0.0;    a2=-0.5*PI;
      setcurveascircle((pc->curves+15),15,-1,
                       &ri1,&a1,&a2,(-tw/2.0-ri1),(-H/2.0+tf2+ri1),
                       0.0,0.0,0.0,0.0);
    }
  }
  else if(type==CTYPE_CT)
  {
    if(ri1>0.0) pc->ncurve=10;
    else        pc->ncurve=8;

    pc->curves=(struct curve *)malloc(pc->ncurve
                                      *sizeof(struct curve));

    setcurveasline((pc->curves+ 0), 0, tw/2.0,    -H/2.0,
                                       tw/2.0,     H/2.0-tf1-ri1);
    setcurveasline((pc->curves+ 1), 1, tw/2.0+ri1, H/2.0-tf1,
                                       B/2.0,      H/2.0-tf1);
    setcurveasline((pc->curves+ 2), 2, B/2.0,      H/2.0-tf1,
                                       B/2.0,      H/2.0);
    setcurveasline((pc->curves+ 3), 3, B/2.0,      H/2.0,
                                      -B/2.0,      H/2.0);
    setcurveasline((pc->curves+ 4), 4,-B/2.0,      H/2.0,
                                      -B/2.0,      H/2.0-tf1);
    setcurveasline((pc->curves+ 5), 5,-B/2.0,      H/2.0-tf1,
                                      -tw/2.0-ri1, H/2.0-tf1);
    setcurveasline((pc->curves+ 6), 6,-tw/2.0,     H/2.0-tf1-ri1,
                                      -tw/2.0,    -H/2.0);
    setcurveasline((pc->curves+ 7), 7,-tw/2.0,    -H/2.0,
                                       tw/2.0,    -H/2.0);

    if(ri1>0.0)
    {
      a1=-1.0*PI; a2=-1.5*PI;
      setcurveascircle((pc->curves+8),8,-1,
                       &ri1,&a1,&a2,(tw/2.0+ri1),(H/2.0-tf1-ri1),
                       0.0,0.0,0.0,0.0);
      a1=-1.5*PI; a2=-2.0*PI;
      setcurveascircle((pc->curves+9),9,-1,
                       &ri1,&a1,&a2,(-tw/2.0-ri1),(H/2.0-tf1-ri1),
                       0.0,0.0,0.0,0.0);
    }
  }
  else if(type==CTYPE_O)
  {
    if(ro1<=0.0) return;

    if(ri1>0.0) pc->ncurve=2;
    else        pc->ncurve=1;

    pc->curves=(struct curve *)malloc(pc->ncurve
                                      *sizeof(struct curve));

    a1=-2.0*PI; a2=0.0;
    setcurveascircle((pc->curves+0),0,1,
                     &ro1,&a1,&a2,0.0,0.0,
                     0.0,0.0,0.0,0.0);
    if(ri1>0.0)
    {
      a1=0.0; a2=-2.0*PI;
      setcurveascircle((pc->curves+1),1,-1,
                       &ri1,&a1,&a2,0.0,0.0,
                       0.0,0.0,0.0,0.0);
    }
  }

  return;
}/*setpolycurveaskatakou*/

void addcurvestopolycurve(struct polycurve *pc,
                          struct polycurve *add)
{
  int i,j;

  pc->curves=(struct curve *)
             realloc(pc->curves,
                     (pc->ncurve+add->ncurve)*sizeof(struct curve));

  j=0;
  for(i=(pc->ncurve);i<(pc->ncurve+add->ncurve);i++)
  {
    initializecurve(pc->curves+i);

    (pc->curves+i)->dots[0]=(struct onode *)
                            malloc(sizeof(struct onode));
    (pc->curves+i)->dots[1]=(struct onode *)
                            malloc(sizeof(struct onode));
    if((add->curves+j)->type==CTYPE_CIRCLE)
    {
      (pc->curves+i)->center =(struct onode *)
                              malloc(sizeof(struct onode));
    }

    copycurve((pc->curves+i),(add->curves+j));
    j++;
  }

  (pc->ncurve)+=(add->ncurve);

  return;
}/*addcurvestopolycurve*/

void deletecurveinpolycurve(struct polycurve *pc,int offset)
{
  int i;

  pc->ncurve--;

  for(i=offset;i<(pc->ncurve);i++)
  {
    copycurve((pc->curves+i),(pc->curves+i+1));
  }

  pc->curves=(struct curve *)
             realloc(pc->curves,(pc->ncurve)*sizeof(struct curve));

  return;
}/*deletecurveinpolycurve*/

void freepolycurve(struct polycurve *pc)
{
  int i;

  for(i=0;i<(pc->ncurve);i++)
  {
    if((pc->curves+i)->center!=NULL) free((pc->curves+i)->center);
    if((pc->curves+i)->dots[0]!=NULL) free((pc->curves+i)->dots[0]);
    if((pc->curves+i)->dots[1]!=NULL) free((pc->curves+i)->dots[1]);
    if((pc->curves+i)->dots[2]!=NULL) free((pc->curves+i)->dots[2]);
    if((pc->curves+i)->tan[0]!=NULL) free((pc->curves+i)->tan[0]);
    if((pc->curves+i)->tan[1]!=NULL) free((pc->curves+i)->tan[1]);
    if((pc->curves+i)->tan[2]!=NULL) free((pc->curves+i)->tan[2]);
  }
  free(pc->curves);

  pc->loff=0;
  pc->ncurve=0;
  pc->curves=NULL;

  return;
}/*freepolycurve*/

void freepolypolycurve(struct polypolycurve *ppc)
{
  int i;

  for(i=0;i<(ppc->npcurve);i++)
  {
    freepolycurve(ppc->pcurves+i);
  }
  free(ppc->pcurves);

  ppc->npcurve=0;
  ppc->pcurves=NULL;

  return;
}/*freepolypolycurve*/

int cmqpolygon(struct obans cban,double wban,
               struct cmqelem *csum)
/*CMQ BY POLYGONAL LOAD.INPUT BAN AS POLYGON FORM.*/
/*GIRDER ON 1ST LINE OF BAN.*/
/*wban=WEIGHT OF BAN*/
{
  char non[256];
  int idot,jdot;
  int ndot;
  double w0,L0,L1,L2,L3;
  double H0;
  double **drccos;
  struct cmqelem cmq;
  struct onode *lnods;
  struct obans ban;

  csum->Ci=0.0;
  csum->Cj=0.0;
  csum->Mc0=0.0;
  csum->Qi0=0.0;
  csum->Qj0=0.0;

  ndot=cban.nnod;

  /*PROJECTED BAN*/
  ban.nods=(struct onode **)malloc(ndot*sizeof(struct onode *));
  lnods=(struct onode *)malloc(ndot*sizeof(struct onode));

  drccos=filmdrccos(*(*(cban.nods+0)),
                    *(*(cban.nods+1)),
                    *(*(cban.nods+ndot-1)));    /*DIRECTION COSINE.*/

  for(idot=0;idot<ndot;idot++)       /*TRANSFORM GLOBAL INTO LOCAL.*/
  {
    *(lnods+idot)=transgtol(drccos,*(*(cban.nods+0)),
                                   *(*(cban.nods+idot)));
    *(ban.nods+idot)=(lnods+idot);
  }

  L0=((*(ban.nods+1))->d[EX])-((*(ban.nods+0))->d[EX]);   /*LENGTH.*/

  for(idot=ndot-1;idot>1;idot--)
  {
    H0=(*(ban.nods+idot))->d[EY]; /*Yi=HEIGHT*/
    w0=wban*H0;

    /*          wi          */
    /*          .           */
    /*        .::           */
    /*      .::::           */
    /* ....::::::.......... */
    /*    Xj    Xi          */
    if(idot==ndot-1) jdot=0;
    else             jdot=idot+1;

    L1=(*(ban.nods+jdot))->d[EX];
    L2=(*(ban.nods+idot))->d[EX]-L1;
    L3=L0-L1-L2;
    if(L2>0.0)
    {
      if(!cmqrightangledtriangle(w0,L0,L1,L2,L3,&cmq))
      {
        sprintf(non,"Failed 1. EPS=L0-(L1+L2+L3)=%.5E",
                L0-L1-L2-L3);
        MessageBox(NULL,non,"CMQ Polygon",MB_OK);
        return 0;
      }
      csum->Ci+=cmq.Ci;
      csum->Cj+=cmq.Cj;
      csum->Mc0+=cmq.Mc0;
      csum->Qi0+=cmq.Qi0;
      csum->Qj0+=cmq.Qj0;
    }
/*
sprintf(non,"w0=%.5f[tf] L0=%.5f+%.5f+%.5f=%.5f[m]",w0,L1,L2,L3,L0);
MessageBox(NULL,non,"CMQ",MB_OK);
sprintf(non,"Ci=%.5f[tfm] Cj=%.5f[tfm] Mc=%.5f[tfm]",
        cmq.Ci,cmq.Cj,cmq.Mc0);
MessageBox(NULL,non,"CMQ",MB_OK);
sprintf(non,"Qi=%.5f[tf] Qj=%.5f[tf]",cmq.Qi0,cmq.Qj0);
MessageBox(NULL,non,"CMQ",MB_OK);
*/
    /*          wi          */
    /*          .           */
    /*          :.          */
    /*          ::.         */
    /* .........:::........ */
    /*          Xi Xj       */
    jdot=idot-1;

    L1=L0-((*(ban.nods+jdot))->d[EX]);
    L2=((*(ban.nods+jdot))->d[EX])-((*(ban.nods+idot))->d[EX]);
    L3=L0-L1-L2;
    if(L2>0.0)
    {
      if(!cmqrightangledtriangle(w0,L0,L1,L2,L3,&cmq))
      {
        sprintf(non,"Failed 2. EPS=L0-(L1+L2+L3)=%.5E",
                L0-L1-L2-L3);
        MessageBox(NULL,non,"CMQ Polygon",MB_OK);
        return 0;
      }
      csum->Ci-=cmq.Cj;
      csum->Cj-=cmq.Ci;
      csum->Mc0+=cmq.Mc0;
      csum->Qi0+=cmq.Qj0;
      csum->Qj0+=cmq.Qi0;
    }
/*
sprintf(non,"w0=%.5f[tf] L0=%.5f+%.5f+%.5f=%.5f[m]",w0,L1,L2,L3,L0);
MessageBox(NULL,non,"CMQ",MB_OK);
sprintf(non,"Ci=%.5f[tfm] Cj=%.5f[tfm] Mc=%.5f[tfm]",
        cmq.Ci,cmq.Cj,cmq.Mc0);
MessageBox(NULL,non,"CMQ",MB_OK);
sprintf(non,"Qi=%.5f[tf] Qj=%.5f[tf]",cmq.Qi0,cmq.Qj0);
MessageBox(NULL,non,"CMQ",MB_OK);
*/
  }

  /*
  sprintf(non,"SUM:Ci=%.3f[tfm] Cj=%.3f[tfm] Mc=%.3f[tfm]",
          csum->Ci,csum->Cj,csum->Mc0);
  MessageBox(NULL,non,"CMQ",MB_OK);

  sprintf(non,"SUM:Qi=%.5f[tf] Qj=%.5f[tf]",csum->Qi0,csum->Qj0);
  MessageBox(NULL,non,"CMQ",MB_OK);
  */

  free(lnods);
  free(ban.nods);
  for(idot=2;idot>=0;idot--) free(*(drccos+idot));
  free(drccos);

  return 1;
}/*cmqpolygon*/

int cmquniform(struct oelem elem,double w0,
               struct cmqelem *csum)
{
  double dx,dy,dz,L0,L1,L2,L3;
  struct cmqelem cmq;

  dx=(*(elem.nods+1))->d[GX]-(*(elem.nods+0))->d[GX];
  dy=(*(elem.nods+1))->d[GY]-(*(elem.nods+0))->d[GY];
  dz=(*(elem.nods+1))->d[GZ]-(*(elem.nods+0))->d[GZ];

  L0=sqrt(dx*dx+dy*dy+dz*dz);

  L1=0.0;
  L2=L0;
  L3=0.0;
  if(!cmqrightangledtriangle(w0,L0,L1,L2,L3,&cmq)) return 0;

  csum->Ci=cmq.Ci-cmq.Cj;
  csum->Cj=cmq.Cj-cmq.Ci;
  csum->Mc0=2.0*cmq.Mc0;
  csum->Qi0=cmq.Qi0+cmq.Qj0;
  csum->Qj0=cmq.Qj0+cmq.Qi0;

  return 1;
}/*cmquniform*/

int cmqconcentration(double W0,
                     double L0,
                     double L1,double L2,
                     struct cmqelem *cmq)
/*LOAD:CONCENTRATED.*/
/*L0:LENGTH OF ELEMENT = L1+L2*/
/*W0:LOAD AT z=L1*/
{
  if(L0!=(L1+L2)) return 0;

  cmq->Ci=-W0*L1*L2*L2/(L0*L0);
  cmq->Cj=W0*L1*L1*L2/(L0*L0);

  cmq->Qi0=W0*L2/L0;
  cmq->Qj0=W0*L1/L0;

  if(L1>=L2) cmq->Mc0=0.5*W0*L2;
  if(L1<L2)  cmq->Mc0=0.5*W0*L1;

  return 1;
}/*cmqconcentration*/

int cmqrightangledtriangle(double w0,
                           double L0,
                           double L1,double L2,double L3,
                           struct cmqelem *cmq)
/*LOAD:TRIANGLE BETWEEN z=L1 AND z=L1+L2 WITH RIGHT ANGLE.*/
/*L0:LENGTH OF ELEMENT = L1+L2+L3*/
/*w0:PEAK OF TRIANGLE LOAD AT z=L1+L2*/
/*            w0          */
/*            .           */
/*          .::           */
/*        .::::           */
/* ......::::::.......... */
/*   L1    L2      L3     */
/*          L0            */
{
  double eps=1.0E-08;

  if(fabs(L0-(L1+L2+L3))>eps) return 0;

  cmq->Ci=-1.0/60.0*w0*L2/(L0*L0)
          *(5.0*(L2*L2+4.0*L2*L3+ 6.0*L3*L3)*L1
           +2.0*(L2*L2+5.0*L2*L3+10.0*L3*L3)*L2);
  cmq->Cj=1.0/60.0*w0*L2/(L0*L0)
          *(5.0*( 6.0*L1*L1+ 8.0*L1*L2+3.0*L2*L2)*L3
               +(10.0*L1*L1+10.0*L1*L2+3.0*L2*L2)*L2);

  cmq->Qi0=1.0/6.0*w0*L2/L0*(L2+3.0*L3);
  cmq->Qj0=1.0/6.0*w0*L2/L0*(3.0*L1+2.0*L2);

  if(0.5*L0<=L1) cmq->Mc0=0.5*L0*cmq->Qi0;
  if(L1<0.5*L0 && 0.5*L0<(L1+L2))
  {
    cmq->Mc0=1.0/48.0*w0/L2
             *(3.0*L2*L2*(L1+L2+3.0*L3)
               +(L1-L3)*(L1-L3)*(L1-3.0*L2-L3));
  }
  if((L1+L2)<=0.5*L0) cmq->Mc0=0.5*L0*cmq->Qj0;

  return 1;
}/*cmqrightangledtriangle*/

void weightdistribution(HDC hdc,FILE *fout,struct viewparam *vp,
						struct organ *org)
/*DISTRIBUTE ELEMENT WEIGHT INTO NODE WEIGHT,CMQ.*/
{
  char non[4096],nonII[4096],ss[512];
  char ss1[512],ss2[512],ss3[512],ss4[512];
  char ss5[512],ss6[512],ss7[512],ss8[512];
  char ss11[512],ss12[512],ss13[512],ss14[512];
  char ss15[512],ss16[512],ss17[512],ss18[512];
  int i,j,k,ii,jj,kk,ifind,ncount;
  int bdot;
  long int loff;
  double weight,w[3],total[3],radius,orate,hiju,alpha;
  double dx,dy,dz,length,height,lo,ho;
  double *fload,*volumes;
  struct onode *gn,*gn1,*gn2,*bn1,*bn2;
  struct oelem *ge;
  struct obans *cmqbans;
  struct cmqelem cmqsum[3],volume;
  struct oload *nload; /*NODE LOADS.*/
  struct aiparam ai;

  /*PROJECT NAME*/
  GetDlgItemText((wmenu.childs+2)->hwnd,ID_SECTIONFILE,non,80);
  if(!strncmp(non,"yachi",5))     strcpy(prj,"yachi");
  else if(!strncmp(non,"yama",4)) strcpy(prj,"yama");
  else if(!strncmp(non,"pan",3))  strcpy(prj,"pan");
  else if(!strncmp(non,"huji",4)) strcpy(prj,"huji");
  else if(!strncmp(non,"hira",4)) strcpy(prj,"hira");
  else if(!strncmp(non,"noa",3))  strcpy(prj,"noa");
  else if(!strncmp(non,"mie",3))  strcpy(prj,"mie");
  else if(!strncmp(non,"naka",4)) strcpy(prj,"naka");
  else if(!strncmp(non,"kiku",4)) strcpy(prj,"kiku");
  else if(!strncmp(non,"engo",4)) strcpy(prj,"engo");
  else if(!strncmp(non,"oos",3))  strcpy(prj,"oos");
  else if(!strncmp(non,"ama",3))  strcpy(prj,"ama");
  else if(!strncmp(non,"koba",4)) strcpy(prj,"koba");
  else if(!strncmp(non,"nade",4)) strcpy(prj,"nade");
  else if(!strncmp(non,"karu",4)) strcpy(prj,"karu");
  else if(!strncmp(non,"ashi",4)) strcpy(prj,"ashi");
  else if(!strncmp(non,"izu",3))  strcpy(prj,"izu");
  else if(!strncmp(non,"eifu",4)) strcpy(prj,"eifu");
  else if(!strncmp(non,"gan",3))  strcpy(prj,"gan");
  else if(!strncmp(non,"tsurumi",7))  strcpy(prj,"tsurumi");
  else if(!strncmp(non,"kiyo",4)) strcpy(prj,"kiyo");
  else if(!strncmp(non,"miya",4)) strcpy(prj,"miya");
  else if(!strncmp(non,"banga",5)) strcpy(prj,"banga");
  else if(!strncmp(non,"tep",3)) strcpy(prj,"tep");
  else if(!strncmp(non,"koshi",5)) strcpy(prj,"koshi");
  else if(!strncmp(non,"vene",4)) strcpy(prj,"vene");
  else if(!strncmp(non,"toyo",4)) strcpy(prj,"toyo");
  else if(!strncmp(non,"hakone",6)) strcpy(prj,"hakone");
  else if(!strncmp(non,"sajigym",7)) strcpy(prj,"sajigym");
  else if(!strncmp(non,"toko",4)) strcpy(prj,"toko");
  else if(!strncmp(non,"maebasi",7)) strcpy(prj,"maebasi");
  else if(!strncmp(non,"dazai",5)) strcpy(prj,"dazai");
  else if(!strncmp(non,"akita",5)) strcpy(prj,"akita");
  else if(!strncmp(non,"agc",3)) strcpy(prj,"agc");
  else if(!strncmp(non,"minna",5)) strcpy(prj,"minna");
  else if(!strncmp(non,"riku",4)) strcpy(prj,"riku");
  else if(!strncmp(non,"sichi",5)) strcpy(prj,"sichi");
  else if(!strncmp(non,"ebis",4)) strcpy(prj,"ebis");
  else if(!strncmp(non,"stick",5)) strcpy(prj,"stick");
  else if(!strncmp(non,"sydney",6)) strcpy(prj,"sydney");
  else if(!strncmp(non,"ghouse",6)) strcpy(prj,"ghouse");   /*201113ghouse*/
  else if(!strncmp(non,"LunarBase",9)) strcpy(prj,"lunarbase");
  else if(!strncmp(non,"optimization",12)) strcpy(prj,"optimization"); //UJIOKA FOR OPTIMIZATION

  else                            strcpy(prj,"none");

  org->loads=(struct oload *)
			 malloc((org->nnode)*sizeof(struct oload));
  fload=(double *)malloc((org->nnode)*sizeof(double));
  for(i=0;i<(org->nnode);i++)
  {
	(org->loads+i)->nod=org->nodes+i;

	(org->loads+i)->w[WEIGHTSLAB] =0.0;
	(org->loads+i)->w[WEIGHTFRAME]=0.0;
	(org->loads+i)->w[WEIGHTEQ]   =0.0;

	(org->loads+i)->w[WEIGHTSLAB] =-(org->confs+6*i+2)->value;
	(org->loads+i)->w[WEIGHTFRAME]=-(org->confs+6*i+2)->value;
	(org->loads+i)->w[WEIGHTEQ]  =-(org->confs+6*i+2)->value;

	(org->loads+i)->mass=0.0;

	*(fload+i)=0.0;
  }
  nload=org->loads;

  volumes=(double *)malloc((org->nsect)*sizeof(double));
  for(i=0;i<(org->nsect);i++) *(volumes+i)=0.0;

  for(i=0;i<(org->nelem);i++)
  {
	currentpivot(i+1,org->nelem);

	ge=org->elems+i;

	/*SLAB*/
	if(ge->type==SLAB)
	{
	  if((ge->nnod)<3)
	  {
		sprintf(non,"Slab %ld:Less Nodes.",ge->code);
		MessageBox(NULL,non,"Distribute",MB_OK);
	  }

	  weight=0.0;
	  for(j=0;j<(ge->sect->nfig);j++)
	  {
		weight+=((ge->sect->figs+j)->thick)
				*((ge->sect->figs+j)->prop->hiju);
	  }

	  /*UNIT LOADS.*/
	  w[WEIGHTSLAB] =weight+ge->sect->lload[WEIGHTSLAB];
	  w[WEIGHTFRAME]=weight+ge->sect->lload[WEIGHTFRAME];
	  w[WEIGHTEQ]   =weight+ge->sect->lload[WEIGHTEQ];

	  /*if(w[WEIGHTSLAB] <=0.0 &&
		 w[WEIGHTFRAME]<=0.0 &&
		 w[WEIGHTEQ]   <=0.0)
	  {
		sprintf(non,"Slab %ld:No Weight.",ge->sect->code);
		MessageBox(NULL,non,"Distribute",MB_OK);
	  }*/
	  /*sprintf(non,"Slab Weight=%.3f",weight);
	  MessageBox(NULL,non,"Distribute",MB_OK);*/

	  for(j=0;j<(ge->nban);j++)
	  {
		bdot=(ge->bans+j)->nnod;
		cmqbans=(struct obans *)malloc(bdot*sizeof(struct obans));

		slabdivision2(hdc,vp,*(ge->bans+j),cmqbans);

		for(k=0;k<bdot;k++)
		{
		  cmqpolygon(*(cmqbans+k),
					 w[WEIGHTSLAB], &(cmqsum[WEIGHTSLAB]));
		  cmqpolygon(*(cmqbans+k),
					 w[WEIGHTFRAME],&(cmqsum[WEIGHTFRAME]));
		  cmqpolygon(*(cmqbans+k),
					 w[WEIGHTEQ],   &(cmqsum[WEIGHTEQ]));

		  bn1=*((ge->bans+j)->nods+k);
		  if(k==bdot-1) bn2=*((ge->bans+j)->nods+0);
		  else          bn2=*((ge->bans+j)->nods+k+1);

		  ifind=0;

if(strncmp(prj,"minna",5)){ /*FOR WITHOUT CMQ*/

		  for(ii=0;ii<(org->nelem);ii++) /*SEARCH GIRDER.*/
		  {
			if((org->elems+ii)->type==GIRDER)
			{
			  gn1=*((org->elems+ii)->nods+0);
			  gn2=*((org->elems+ii)->nods+1);

			  /*HORIZONTAL GIRDER AVAILABLE.*/
			  if(gn1==bn1 && gn2==bn2)
			  {
				/*DO NOT INPUT CMQ IF EDGE OF WALL.*/
				for(jj=0;jj<(org->nelem);jj++)
				{
				  ncount=0;
				  if((org->elems+jj)->type==WALL &&
					 (((org->elems+jj)->sect->figs+0)->prop->E)>0.0)
				  {
					for(kk=0;kk<(org->elems+jj)->nnod;kk++)
					{
					  if(gn1==*((org->elems+jj)->nods+kk)) ncount++;
					  if(gn2==*((org->elems+jj)->nods+kk)) ncount++;
                    }
                    if(ncount>=2) break;
                  }
                }
                if(ncount<2)
                {
                  (org->elems+ii)->initial[0][2]+=cmqsum[WEIGHTFRAME].Qi0;
                  (org->elems+ii)->initial[1][2]+=cmqsum[WEIGHTFRAME].Qj0;
                  (org->elems+ii)->initial[0][4]+=cmqsum[WEIGHTFRAME].Ci;
                  (org->elems+ii)->initial[1][4]+=cmqsum[WEIGHTFRAME].Cj;
/*
sprintf(non,"ELEM %ld:Qi=%.3f Qj=%.3f Ci=%.3f Cj=%.3f",
        (org->elems+ii)->code,
        (org->elems+ii)->initial[0][2],
        (org->elems+ii)->initial[1][2],
        (org->elems+ii)->initial[0][4],
        (org->elems+ii)->initial[1][4]);
MessageBox(NULL,non,"Distribute",MB_OK);
*/
                  ifind=1;
                }
                break;
              }
              else if(ncount<2 && gn1==bn2 && gn2==bn1) /*REVERSED.*/
			  {
                /*DO NOT INPUT CMQ IF EDGE OF WALL.*/
                for(jj=0;jj<(org->nelem);jj++)
                {
                  ncount=0;
                  if((org->elems+jj)->type==WALL &&
                     (((org->elems+jj)->sect->figs+0)->prop->E)>0.0)
                  {
                    for(kk=0;kk<(org->elems+jj)->nnod;kk++)
                    {
                      if(gn1==*((org->elems+jj)->nods+kk)) ncount++;
                      if(gn2==*((org->elems+jj)->nods+kk)) ncount++;
                    }
                    if(ncount>=2) break;
                  }
                }
                if(ncount<2)
                {
                  (org->elems+ii)->initial[0][2]+=cmqsum[WEIGHTFRAME].Qj0;
                  (org->elems+ii)->initial[1][2]+=cmqsum[WEIGHTFRAME].Qi0;
                  (org->elems+ii)->initial[0][4]+=-cmqsum[WEIGHTFRAME].Cj;
                  (org->elems+ii)->initial[1][4]+=-cmqsum[WEIGHTFRAME].Ci;
                  ifind=1;
                }
                break;
              }
            }
          }

} /*FOR WITHOUT CMQ*/

          if(!ifind)/*CASE NO GIRDER*/
          {
            loff=bn1->loff;
            (nload+loff)->w[WEIGHTFRAME]+=cmqsum[WEIGHTFRAME].Qi0;
            loff=bn2->loff;
            (nload+loff)->w[WEIGHTFRAME]+=cmqsum[WEIGHTFRAME].Qj0;
          }

          loff=bn1->loff;
          (nload+loff)->w[WEIGHTSLAB]+=cmqsum[WEIGHTSLAB].Qi0;
          (nload+loff)->w[WEIGHTEQ]  +=cmqsum[WEIGHTEQ].Qi0;
          loff=bn2->loff;
          (nload+loff)->w[WEIGHTSLAB]+=cmqsum[WEIGHTSLAB].Qj0;
		  (nload+loff)->w[WEIGHTEQ]  +=cmqsum[WEIGHTEQ].Qj0;

          cmqpolygon(*(cmqbans+k),1.0,&volume);
          *(volumes+(ge->sect->loff))+=volume.Qi0;
          *(volumes+(ge->sect->loff))+=volume.Qj0;
/*
sprintf(non,"NODE %ld:W=%.3f,NODE %ld:W=%.3f",
        bn1->code,cmqsum.Qi0,
        bn2->code,cmqsum.Qj0);
MessageBox(NULL,non,"Distribute",MB_OK);
*/
        }

        for(ii=(bdot-1);ii>=0;ii--)
        {
          for(jj=((cmqbans+ii)->nnod)-1;jj>=0;jj--)
          {
            free(*((cmqbans+ii)->nods+jj));
          }
          free((cmqbans+ii)->nods);
        }
        free(cmqbans);
      }
    }

    /*RESISTING WALL*/
    if(ge->type==WALL)
    {
      if((ge->nnod)<3)
      {
        sprintf(non,"Wall %ld:Less Nodes.",ge->code);
        MessageBox(NULL,non,"Distribute",MB_OK);
      }

	  weight=0.0;
      for(j=0;j<(ge->sect->nfig);j++)
      {
        weight+=((ge->sect->figs+j)->thick)
               *((ge->sect->figs+j)->prop->hiju);
      }
      /*if(weight<=0.0)
      {
        sprintf(non,"Wall %ld:No Weight.",ge->sect->code);
        MessageBox(NULL,non,"Distribute",MB_OK);
	  }*/
      /*sprintf(non,"Wall Weight=%.3f",weight);
      MessageBox(NULL,non,"Distribute",MB_OK);*/

      /*FACE*/
      if((ge->nnod)>=4)
      {
        dx=((*(ge->bans->nods+1))->d[0])
          -((*(ge->bans->nods+0))->d[0]);
        dy=((*(ge->bans->nods+1))->d[1])
          -((*(ge->bans->nods+0))->d[1]);
        dz=((*(ge->bans->nods+1))->d[2])
          -((*(ge->bans->nods+0))->d[2]);
        length=sqrt(dx*dx+dy*dy+dz*dz);

        dx=((*(ge->bans->nods+3))->d[0])
          -((*(ge->bans->nods+0))->d[0]);
        dy=((*(ge->bans->nods+3))->d[1])
          -((*(ge->bans->nods+0))->d[1]);
        dz=((*(ge->bans->nods+3))->d[2])
          -((*(ge->bans->nods+0))->d[2]);
        height=sqrt(dx*dx+dy*dy+dz*dz);

        lo=length-(ge->lface[0])-(ge->lface[1]);
        ho=height-(ge->hface[0])-(ge->hface[1]);

        /*if(lo>0.0) weight*=lo/length;*/
        /*if(ho>0.0) weight*=ho/height;*/
      }

      for(j=0;j<(ge->nban);j++)
	  {
        bdot=(ge->bans+j)->nnod;
        cmqbans=(struct obans *)malloc(bdot*sizeof(struct obans));

        slabdivision2(hdc,vp,*(ge->bans+j),cmqbans);

        for(k=0;k<bdot;k++)
        {
          cmqpolygon(*(cmqbans+k),weight,&cmqsum[0]);

          bn1=*((ge->bans+j)->nods+k);
          if(k==bdot-1) bn2=*((ge->bans+j)->nods+0);
          else          bn2=*((ge->bans+j)->nods+k+1);

          loff=bn1->loff;
          (nload+loff)->w[WEIGHTSLAB] +=cmqsum[0].Qi0;
          (nload+loff)->w[WEIGHTFRAME]+=cmqsum[0].Qi0;
          (nload+loff)->w[WEIGHTEQ]   +=cmqsum[0].Qi0;
          loff=bn2->loff;
          (nload+loff)->w[WEIGHTSLAB] +=cmqsum[0].Qj0;
          (nload+loff)->w[WEIGHTFRAME]+=cmqsum[0].Qj0;
          (nload+loff)->w[WEIGHTEQ]   +=cmqsum[0].Qj0;

          cmqpolygon(*(cmqbans+k),1.0,&volume);
          *(volumes+(ge->sect->loff))+=volume.Qi0;
          *(volumes+(ge->sect->loff))+=volume.Qj0;
/*
sprintf(non,"NODE %ld:W=%.3f,NODE %ld:W=%.3f",
        bn1->code,cmqsum.Qi0,
        bn2->code,cmqsum.Qj0);
MessageBox(NULL,non,"Distribute",MB_OK);
*/
        }

        for(ii=(bdot-1);ii>=0;ii--)
        {
          for(jj=((cmqbans+ii)->nnod)-1;jj>=0;jj--)
          {
            free(*((cmqbans+ii)->nods+jj));
          }
          free((cmqbans+ii)->nods);
        }
        free(cmqbans);
      }
    }

    /*GIRDER*/
    if(ge->type==GIRDER)
    {
      if((ge->nnod)>2)
      {
        sprintf(non,"Girder %ld:Too Many Nodes.",ge->code);
        MessageBox(NULL,non,"Distribute",MB_OK);
      }

      /*weight=0.0;
      for(j=0;j<(ge->sect->nfig);j++)
	  {
        weight+=((ge->sect->figs+j)->area)
               *((ge->sect->figs+j)->prop->hiju);
      }*/
      sectionfeatures(ge->sect);
      w[WEIGHTSLAB] =ge->sect->hiju[WEIGHTSLAB];
      w[WEIGHTFRAME]=ge->sect->hiju[WEIGHTFRAME];
      w[WEIGHTEQ]   =ge->sect->hiju[WEIGHTEQ];

      /*if(w[WEIGHTFRAME]<=0.0)
      {
        sprintf(non,"Girder %ld:No Weight.",ge->sect->code);
        MessageBox(NULL,non,"Distribute",MB_OK);
      }*/
	  /*sprintf(non,"Girder Weight=%.3f",weight);
      MessageBox(NULL,non,"Distribute",MB_OK);*/

      cmquniform(*ge,w[WEIGHTSLAB] ,&cmqsum[WEIGHTSLAB]);
      cmquniform(*ge,w[WEIGHTFRAME],&cmqsum[WEIGHTFRAME]);
      cmquniform(*ge,w[WEIGHTEQ]   ,&cmqsum[WEIGHTEQ]);

      /*DO NOT INPUT CMQ IF EDGE OF WALL.*/
      gn1=*(ge->nods+0);
      gn2=*(ge->nods+1);
      for(jj=0;jj<(org->nelem);jj++)
      {
        ncount=0;

        if((org->elems+jj)->type==WALL &&
           (((org->elems+jj)->sect->figs+0)->prop->E)>0.0)
        {
          for(kk=0;kk<(org->elems+jj)->nnod;kk++)
          {
            if(gn1==*((org->elems+jj)->nods+kk)) ncount++;
            if(gn2==*((org->elems+jj)->nods+kk)) ncount++;
          }
          if(ncount>=2) break;
        }
      }

      if(ncount<2) /*NO WALL.*/
      {
        (org->elems+i)->initial[0][2]+=cmqsum[WEIGHTFRAME].Qi0;
        (org->elems+i)->initial[1][2]+=cmqsum[WEIGHTFRAME].Qj0;
        (org->elems+i)->initial[0][4]+=cmqsum[WEIGHTFRAME].Ci;
		(org->elems+i)->initial[1][4]+=cmqsum[WEIGHTFRAME].Cj;
      }
      else /*WITH WALL.*/
      {
        loff=gn1->loff;
        (nload+loff)->w[WEIGHTFRAME]+=cmqsum[WEIGHTFRAME].Qi0;
        loff=gn2->loff;
        (nload+loff)->w[WEIGHTFRAME]+=cmqsum[WEIGHTFRAME].Qj0;
      }

      loff=gn1->loff;
      (nload+loff)->w[WEIGHTSLAB]+=cmqsum[WEIGHTSLAB].Qi0;
      (nload+loff)->w[WEIGHTEQ]  +=cmqsum[WEIGHTEQ].Qi0;
      loff=gn2->loff;
      (nload+loff)->w[WEIGHTSLAB]+=cmqsum[WEIGHTSLAB].Qj0;
      (nload+loff)->w[WEIGHTEQ]  +=cmqsum[WEIGHTEQ].Qj0;

      cmquniform(*ge,1.0,&volume);
      *(volumes+(ge->sect->loff))+=volume.Qi0;
      *(volumes+(ge->sect->loff))+=volume.Qj0;
/*
sprintf(non,"ELEM %ld:dQi=%.3f dQj=%.3f dCi=%.3f dCj=%.3f",
        (org->elems+i)->code,
        cmqsum.Qi0,cmqsum.Qj0,
        cmqsum.Ci, cmqsum.Cj);
MessageBox(NULL,non,"Distribute",MB_OK);

sprintf(non,"ELEM %ld:Qi=%.3f Qj=%.3f Ci=%.3f Cj=%.3f",
        (org->elems+i)->code,
        (org->elems+i)->initial[0][2],
        (org->elems+i)->initial[1][2],
        (org->elems+i)->initial[0][4],
        (org->elems+i)->initial[1][4]);
MessageBox(NULL,non,"Distribute",MB_OK);
*/
    }

    /*COLUMN*/
    if(ge->type==COLUMN)
    {
      if((ge->nnod)>2)
      {
        sprintf(non,"Column %ld:Too Many Nodes.",ge->code);
        MessageBox(NULL,non,"Distribute",MB_OK);
	  }

      /*weight=0.0;
      for(j=0;j<(ge->sect->nfig);j++)
      {
        weight+=((ge->sect->figs+j)->area)
               *((ge->sect->figs+j)->prop->hiju);
      }*/
      sectionfeatures(ge->sect);
      weight=ge->sect->hiju[0];

      /*if(weight<=0.0)
      {
        sprintf(non,"Column %ld:No Weight.",ge->sect->code);
        MessageBox(NULL,non,"Distribute",MB_OK);
      }*/
      /*sprintf(non,"Column Weight=%.3f",weight);
      MessageBox(NULL,non,"Distribute",MB_OK);*/

      cmquniform(*ge,weight,&cmqsum[0]);

      loff=(*(ge->nods+0))->loff;
      (nload+loff)->w[WEIGHTSLAB] +=cmqsum[0].Qi0;
      (nload+loff)->w[WEIGHTFRAME]+=cmqsum[0].Qi0;
      (nload+loff)->w[WEIGHTEQ]   +=cmqsum[0].Qi0;
      loff=(*(ge->nods+1))->loff;
      (nload+loff)->w[WEIGHTSLAB] +=cmqsum[0].Qj0;
      (nload+loff)->w[WEIGHTFRAME]+=cmqsum[0].Qj0;
      (nload+loff)->w[WEIGHTEQ]   +=cmqsum[0].Qj0;

      cmquniform(*ge,1.0,&volume);
      *(volumes+(ge->sect->loff))+=volume.Qi0;
      *(volumes+(ge->sect->loff))+=volume.Qj0;
/*
sprintf(non,"NODE %ld:W=%.3f,NODE %ld:W=%.3f",
        (*(ge->nods+0))->code,cmqsum.Qi0,
        (*(ge->nods+1))->code,cmqsum.Qj0);
MessageBox(NULL,non,"Distribute",MB_OK);
*/
    }

    /*BRACE*/
    if(ge->type==BRACE)
    {
	  if((ge->nnod)>2)
      {
        sprintf(non,"Brace %ld:Too Many Nodes.",ge->code);
        MessageBox(NULL,non,"Distribute",MB_OK);
      }

      /*weight=0.0;
      for(j=0;j<(ge->sect->nfig);j++)
      {
        weight+=((ge->sect->figs+j)->area)
               *((ge->sect->figs+j)->prop->hiju);
      }*/
      sectionfeatures(ge->sect);
      weight=ge->sect->hiju[0];

      /*if(weight<=0.0)
      {
        sprintf(non,"Brace %ld:No Weight.",ge->sect->code);
        MessageBox(NULL,non,"Distribute",MB_OK);
      }*/
      /*sprintf(non,"Brace Weight=%.3f",weight);
      MessageBox(NULL,non,"Distribute",MB_OK);*/

      cmquniform(*ge,weight,&cmqsum[0]);

      loff=(*(ge->nods+0))->loff;
      (nload+loff)->w[WEIGHTSLAB] +=cmqsum[0].Qi0;
      (nload+loff)->w[WEIGHTFRAME]+=cmqsum[0].Qi0;
      (nload+loff)->w[WEIGHTEQ]   +=cmqsum[0].Qi0;
      loff=(*(ge->nods+1))->loff;
      (nload+loff)->w[WEIGHTSLAB] +=cmqsum[0].Qj0;
      (nload+loff)->w[WEIGHTFRAME]+=cmqsum[0].Qj0;
      (nload+loff)->w[WEIGHTEQ]   +=cmqsum[0].Qj0;

      cmquniform(*ge,1.0,&volume);
      *(volumes+(ge->sect->loff))+=volume.Qi0;
      *(volumes+(ge->sect->loff))+=volume.Qj0;
/*
sprintf(non,"NODE %ld:W=%.3f,NODE %ld:W=%.3f",
        (*(ge->nods+0))->code,cmqsum.Qi0,
        (*(ge->nods+1))->code,cmqsum.Qj0);
MessageBox(NULL,non,"Distribute",MB_OK);
*/
    }

    /*BEAM*/
    /*CURTAIN WALL*/
    /*CONCENTRATION WEIGHT*/
  }

  orate=org->opaque; /*OPAQUE IS ALWAYS UPDATED ON DIALOG INPUT.*/
  getmasshiju((wmenu.childs+4)->hwnd,&hiju);

//  if(fout!=NULL) fprintf(fout,"WEIGHT OF NODES\n");
//mihara07.11.04////////////////////////////////////////////////////////////////
  if(fout!=NULL)
  {
   if(globalunit==1.0)
   {
   fprintf(fout,"3.2 : 節点重量\n\n");
   fprintf(fout,"節点ごとに重量を集計した結果を記す。\n");
   fprintf(fout,"柱，壁は階高の中央で上下に分配するものとする。\n\n");
   fprintf(fout," 節点番号          積載荷重別の重量 [tf]\n\n");
   fprintf(fout,"                 床用     柱梁用     地震用\n");
   }
   else
   {
   fprintf(fout,"3.2 : 節点重量\n\n");
   fprintf(fout,"節点ごとに重量を集計した結果を記す。\n");
   fprintf(fout,"柱，壁は階高の中央で上下に分配するものとする。\n\n");
   fprintf(fout," 節点番号          積載荷重別の重量 [kN]\n\n");
   fprintf(fout,"                 床用     柱梁用     地震用\n");
   }
  }
////////////////////////////////////////////////////////////////////////////////

  total[WEIGHTSLAB] =0.0;
  total[WEIGHTFRAME]=0.0;
  total[WEIGHTEQ]  =0.0;
  for(i=0;i<(org->nelem);i++) /*ADD Q OF CMQ.*/
  {
    total[WEIGHTFRAME]+=(org->elems+i)->initial[0][2];
    total[WEIGHTFRAME]+=(org->elems+i)->initial[1][2];

    loff=(*((org->elems+i)->nods+0))->loff;
    *(fload+loff)+=(org->elems+i)->initial[0][2];
	loff=(*((org->elems+i)->nods+1))->loff;
    *(fload+loff)+=(org->elems+i)->initial[1][2];
  }
  for(i=0;i<(org->nnode);i++) /*ADD NODE LOADS.*/
  {
	total[WEIGHTSLAB] +=(nload+i)->w[WEIGHTSLAB];
	total[WEIGHTFRAME]+=(nload+i)->w[WEIGHTFRAME];
	total[WEIGHTEQ]   +=(nload+i)->w[WEIGHTEQ];

	(nload+i)->mass    =((nload+i)->w[WEIGHTEQ])/9.80665;

	*(fload+i)+=(nload+i)->w[WEIGHTFRAME];

	/*
	sprintf(non,"NODE %ld:WEIGHT FOR FRAME=%.3f,EQ=%.3f",
			(org->nodes+i)->code,
			(nload+i)->w[WEIGHTFRAME],
			(nload+i)->w[WEIGHTEQ]);
	MessageBox(NULL,non,"Distribute",MB_OK);
	*/

	if(fout!=NULL)
	{
      fprintf(fout,"%9ld  %10.3f %10.3f %10.3f\n",
			  (org->nodes+i)->code,
			  globalunit*(nload+i)->w[WEIGHTSLAB],
			  globalunit*(*(fload+i)),
			  globalunit*(nload+i)->w[WEIGHTEQ]);
	}
	if(hdc!=NULL) /*DRAW WEIGHT SPHERE.*/
	{
	  radius=pow(3.0/4.0/PI*((nload+i)->w[WEIGHTSLAB]/hiju),1.0/3.0);
	  fillglobalcircle(hdc,*vp,*((nload+i)->nod),radius,
                       150,150,150,orate,SRCPAINT);
	}
  }

  sprintf(non,"\nTOTAL WEIGHT FOR SLAB=%.3f,FRAME=%.3f,EQ=%.3f",
		  globalunit*total[WEIGHTSLAB],
		  globalunit*total[WEIGHTFRAME],
		  globalunit*total[WEIGHTEQ]);
  if(fout!=NULL) fprintf(fout,"%s\n",non);

  if(globalmessageflag==1) MessageBox(NULL,non,"Distribute",MB_OK);


  if(fout!=NULL)
  {
	/*fprintf(fout,"\nVOLUME OF EACH SECTIONS\n");*/
	fprintf(fout,"\n各断面の部材総量（参考資料）\n");
	fprintf(fout,"\n 断面番号   長さ,面積[m,m2]\n");

	for(i=0;i<(org->nsect);i++)
	{
	  fprintf(fout,"      %d %9.3f\n",
			  (org->sects+i)->code,*(volumes+i));
	}
  }

  if(!strncmp(prj,"sydney",6))    /*steel weight for sydney*/     //honda
  {
	totalweight =  globalunit*total[WEIGHTFRAME];
	totallength = 0.0;
	for(i=0;i<(org->nsect);i++)
	{
	  totallength+=*(volumes+i);
	}
    Gcx = 0.0; Gcy = 0.0; Gcz=0.0;
	for(i=0;i<(org->nnode);i++) /*ADD NODE LOADS.*/
	{
	  Gcx += (globalunit*(nload+i)->w[WEIGHTSLAB])*((org->nodes+i)->d[GX]);
	  Gcy += (globalunit*(nload+i)->w[WEIGHTSLAB])*((org->nodes+i)->d[GY]);
	  Gcz += (globalunit*(nload+i)->w[WEIGHTSLAB])*((org->nodes+i)->d[GZ]);
	}
	Gcx =  Gcx/totalweight;
	Gcy =  Gcy/totalweight;
	Gcz =  Gcz/totalweight;
   }

  free(volumes);

  /*Ai LOADS*/
  /*GetDlgItemText((wmenu.childs+4)->hwnd,IDVS_PFACT,non,80);*/
  /*ai.T1=strtod(non,NULL);*/
  ai.T1=org->ai.T1; /*PERIOD FACTOR*/
  if(ai.T1==0.0) ai.T1=0.03;

  /*GetDlgItemText((wmenu.childs+4)->hwnd,IDVS_BASESHEAR,non,80);*/
  /*ai.Co=strtod(non,NULL);*/
//  ai.Co=org->ai.Co; /*PERIOD FACTOR*/
  ai.Cox=org->ai.Cox; /*PERIOD FACTOR*/
  ai.Coy=org->ai.Coy;

//  if(!strncmp(prj,"hakone",6)) ai.Co=0.122*1.25;

  /*if(ai.Co==0.0) ai.Co=0.300;*/ /*BASE SHEAR 0.2,0.3,1.0*/
//  if(ai.Co==0.0) ai.Co=0.000; /*BASE SHEAR 0.2,0.3,1.0*/
  if(ai.Cox==0.0) ai.Cox=0.000; /*BASE SHEAR 0.2,0.3,1.0*/
  if(ai.Coy==0.0) ai.Coy=0.000; /*BASE SHEAR 0.2,0.3,1.0*/

  ai.Tc=org->ai.Tc; /*GROUND PERIOD*/
  if(ai.Tc==0.0) ai.Tc=0.6; /*TYPE 2*/

  ai.Z=org->ai.Z; /*LOCATION FACTOR*/
  if(ai.Z==0.0) ai.Z=1.0;

  ai.nfloor=19;
  ai.hmax=0.0;
  ai.fnames=(char **)malloc(ai.nfloor*sizeof(char *));
  ai.nnode =(int *)malloc(ai.nfloor*sizeof(int));
  ai.lbound=(double *)malloc((ai.nfloor+1)*sizeof(double));
  ai.levels=(double *)malloc(ai.nfloor*sizeof(double));
  ai.wi    =(double *)malloc(ai.nfloor*sizeof(double));
  ai.Wi    =(double *)malloc(ai.nfloor*sizeof(double));
  ai.Ai    =(double *)malloc(ai.nfloor*sizeof(double));
  ai.Cix   =(double *)malloc(ai.nfloor*sizeof(double));
  ai.Qix   =(double *)malloc(ai.nfloor*sizeof(double));
  ai.Hix   =(double *)malloc(ai.nfloor*sizeof(double));
  ai.Ciy   =(double *)malloc(ai.nfloor*sizeof(double));
  ai.Qiy   =(double *)malloc(ai.nfloor*sizeof(double));
  ai.Hiy   =(double *)malloc(ai.nfloor*sizeof(double));
  ai.factsx=(double *)malloc(ai.nfloor*sizeof(double));
  ai.factsy=(double *)malloc(ai.nfloor*sizeof(double));
  for(i=0;i<ai.nfloor;i++)
  {
    *(ai.fnames+i)=(char *)malloc(80*sizeof(char));
    sprintf(*(ai.fnames+i),"%dFL",i+1); /*FLOOR NAMES*/

    *(ai.nnode+i)=0;
    *(ai.levels+i)=0.0;
    *(ai.wi+i)=0.0;
    *(ai.Wi+i)=0.0;
  }

  *(ai.lbound+ 0)= -0.5;
  *(ai.lbound+ 1)= -0.1;
  *(ai.lbound+ 2)=  2.0;
  *(ai.lbound+ 3)=  3.0;
  *(ai.lbound+ 4)=  4.0;
  *(ai.lbound+ 5)=  5.0;
  *(ai.lbound+ 6)=  6.0;
  *(ai.lbound+ 7)=  7.0;
  *(ai.lbound+ 8)=  8.0;
  *(ai.lbound+ 9)=  9.0;
  *(ai.lbound+10)= 10.0;
  *(ai.lbound+11)= 11.0;
  *(ai.lbound+12)= 13.0;
  *(ai.lbound+13)= 16.0;
  *(ai.lbound+14)= 19.0;
  *(ai.lbound+15)= 22.0;
  *(ai.lbound+16)= 25.0;
  *(ai.lbound+17)= 28.0;
  *(ai.lbound+18)= 31.0;
  *(ai.lbound+19)= 40.0;

  if(!strncmp(prj,"yachi",5))
  {
    *(ai.lbound+ 0)= -0.5;
    *(ai.lbound+ 1)=  0.1;
    *(ai.lbound+ 2)= 10.0;
    *(ai.lbound+ 3)= 13.0;
    *(ai.lbound+ 4)= 14.5;
    *(ai.lbound+ 5)= 15.0;
    *(ai.lbound+ 6)= 16.0;
	*(ai.lbound+ 7)= 17.0;
    *(ai.lbound+ 8)= 18.0;
    *(ai.lbound+ 9)= 19.0;
    *(ai.lbound+10)= 20.0;
    *(ai.lbound+11)= 21.0;
    *(ai.lbound+12)= 23.0;
    *(ai.lbound+13)= 26.0;
    *(ai.lbound+14)= 29.0;
    *(ai.lbound+15)= 32.0;
    *(ai.lbound+16)= 35.0;
    *(ai.lbound+17)= 38.0;
    *(ai.lbound+18)= 41.0;
    *(ai.lbound+19)= 50.0;
  }

  if(!strncmp(prj,"yama",4))
  {
    *(ai.lbound+ 0)= -0.5;
    *(ai.lbound+ 1)=  0.1;
    *(ai.lbound+ 2)=  4.0;
    *(ai.lbound+ 3)= 10.0;
    *(ai.lbound+ 4)= 14.5;
    *(ai.lbound+ 5)= 15.0;
    *(ai.lbound+ 6)= 16.0;
    *(ai.lbound+ 7)= 17.0;
    *(ai.lbound+ 8)= 18.0;
    *(ai.lbound+ 9)= 19.0;
    *(ai.lbound+10)= 20.0;
    *(ai.lbound+11)= 21.0;
    *(ai.lbound+12)= 23.0;
    *(ai.lbound+13)= 26.0;
    *(ai.lbound+14)= 29.0;
    *(ai.lbound+15)= 32.0;
    *(ai.lbound+16)= 35.0;
    *(ai.lbound+17)= 38.0;
    *(ai.lbound+18)= 41.0;
    *(ai.lbound+19)= 50.0;
  }

  if(!strncmp(prj,"pan",3))
  {
    *(ai.lbound+ 0)= -0.5;
    *(ai.lbound+ 1)=  0.9;
    *(ai.lbound+ 2)=  9.0;
	*(ai.lbound+ 3)= 10.0;
    *(ai.lbound+ 4)= 14.5;
    *(ai.lbound+ 5)= 15.0;
    *(ai.lbound+ 6)= 16.0;
    *(ai.lbound+ 7)= 17.0;
    *(ai.lbound+ 8)= 18.0;
    *(ai.lbound+ 9)= 19.0;
    *(ai.lbound+10)= 20.0;
    *(ai.lbound+11)= 21.0;
    *(ai.lbound+12)= 23.0;
    *(ai.lbound+13)= 26.0;
    *(ai.lbound+14)= 29.0;
    *(ai.lbound+15)= 32.0;
    *(ai.lbound+16)= 35.0;
    *(ai.lbound+17)= 38.0;
    *(ai.lbound+18)= 41.0;
    *(ai.lbound+19)= 50.0;
  }

  if(!strncmp(prj,"huji",4))
  {
    *(ai.lbound+ 0)= -0.5;
    *(ai.lbound+ 1)=  0.1;
    *(ai.lbound+ 2)=  3.5;
    *(ai.lbound+ 3)=  5.0;
    *(ai.lbound+ 4)=  6.0;
    *(ai.lbound+ 5)=  8.0;
    *(ai.lbound+ 6)= 10.0;
    *(ai.lbound+ 7)= 12.0;
    *(ai.lbound+ 8)= 18.0;
    *(ai.lbound+ 9)= 19.0;
    *(ai.lbound+10)= 20.0;
    *(ai.lbound+11)= 21.0;
    *(ai.lbound+12)= 23.0;
    *(ai.lbound+13)= 26.0;
    *(ai.lbound+14)= 29.0;
    *(ai.lbound+15)= 32.0;
    *(ai.lbound+16)= 35.0;
    *(ai.lbound+17)= 38.0;
    *(ai.lbound+18)= 41.0;
    *(ai.lbound+19)= 50.0;
  }
  if(!strncmp(prj,"hira",4))
  {
	*(ai.lbound+ 0)= -0.5;
    *(ai.lbound+ 1)=  0.1;
    *(ai.lbound+ 2)=  3.2;
    *(ai.lbound+ 3)=  3.5;
    *(ai.lbound+ 4)=  6.3;
    *(ai.lbound+ 5)=  8.0;
    *(ai.lbound+ 6)= 10.0;
    *(ai.lbound+ 7)= 12.0;
    *(ai.lbound+ 8)= 18.0;
    *(ai.lbound+ 9)= 19.0;
    *(ai.lbound+10)= 20.0;
    *(ai.lbound+11)= 21.0;
    *(ai.lbound+12)= 23.0;
    *(ai.lbound+13)= 26.0;
    *(ai.lbound+14)= 29.0;
    *(ai.lbound+15)= 32.0;
    *(ai.lbound+16)= 35.0;
    *(ai.lbound+17)= 38.0;
    *(ai.lbound+18)= 41.0;
    *(ai.lbound+19)= 50.0;
  }
  if(!strncmp(prj,"noa",3))
  {
    *(ai.lbound+ 0)= -0.5;
    *(ai.lbound+ 1)=  0.1;
    *(ai.lbound+ 2)=  3.5;
    *(ai.lbound+ 3)=  4.5;
    *(ai.lbound+ 4)=  5.5;
    *(ai.lbound+ 5)=  7.6;
    *(ai.lbound+ 6)= 15.0;
    *(ai.lbound+ 7)= 16.0;
    *(ai.lbound+ 8)= 17.0;
    *(ai.lbound+ 9)= 18.0;
    *(ai.lbound+10)= 20.0;
    *(ai.lbound+11)= 21.0;
    *(ai.lbound+12)= 23.0;
    *(ai.lbound+13)= 26.0;
    *(ai.lbound+14)= 29.0;
    *(ai.lbound+15)= 32.0;
    *(ai.lbound+16)= 35.0;
    *(ai.lbound+17)= 38.0;
    *(ai.lbound+18)= 41.0;
    *(ai.lbound+19)= 50.0;
  }
  if(!strncmp(prj,"mie",3))
  {
    *(ai.lbound+ 0)= -0.5;
    *(ai.lbound+ 1)=  2.0;
    *(ai.lbound+ 2)=  3.5;
    *(ai.lbound+ 3)=  4.0;
    *(ai.lbound+ 4)=  4.5;
    *(ai.lbound+ 5)=  5.0;
    *(ai.lbound+ 6)=  5.5;
    *(ai.lbound+ 7)=  7.0;
    *(ai.lbound+ 8)= 17.0;
    *(ai.lbound+ 9)= 18.0;
    *(ai.lbound+10)= 20.0;
    *(ai.lbound+11)= 21.0;
    *(ai.lbound+12)= 23.0;
    *(ai.lbound+13)= 26.0;
    *(ai.lbound+14)= 29.0;
    *(ai.lbound+15)= 32.0;
    *(ai.lbound+16)= 35.0;
    *(ai.lbound+17)= 38.0;
    *(ai.lbound+18)= 41.0;
    *(ai.lbound+19)= 50.0;
  }
  if(!strncmp(prj,"naka",4))
  {
    *(ai.lbound+ 0)= -0.5;
    *(ai.lbound+ 1)=  0.1;
    *(ai.lbound+ 2)=  0.7;
    *(ai.lbound+ 3)=  1.2;
    *(ai.lbound+ 4)=  1.7;
    *(ai.lbound+ 5)=  2.2;
    *(ai.lbound+ 6)=  2.7;
    *(ai.lbound+ 7)=  3.0;
    *(ai.lbound+ 8)=  3.5;
    *(ai.lbound+ 9)=  4.0;
    *(ai.lbound+10)=  4.5;
    *(ai.lbound+11)=  5.0;
    *(ai.lbound+12)=  5.3;
    *(ai.lbound+13)=  5.8;
    *(ai.lbound+14)=  6.3;
    *(ai.lbound+15)=  6.8;
    *(ai.lbound+16)=  7.1;
    *(ai.lbound+17)=  7.6;
    *(ai.lbound+18)=  8.1;
	*(ai.lbound+19)= 41.0;
  }
  if(!strncmp(prj,"kiku",4))
  {
    *(ai.lbound+ 0)= -0.5;
    *(ai.lbound+ 1)=  0.1;
    *(ai.lbound+ 2)=  0.5;
    *(ai.lbound+ 3)=  1.0;
    *(ai.lbound+ 4)=  1.6;
    *(ai.lbound+ 5)=  2.2;
    *(ai.lbound+ 6)=  2.8;
    *(ai.lbound+ 7)=  3.4;
    *(ai.lbound+ 8)=  4.0;
    *(ai.lbound+ 9)=  4.8;
    *(ai.lbound+10)=  5.4;
    *(ai.lbound+11)=  6.0;
    *(ai.lbound+12)=  6.6;
    *(ai.lbound+13)=  7.2;
    *(ai.lbound+14)=  7.8;
    *(ai.lbound+15)=  9.5;
    *(ai.lbound+16)= 10.2;
    *(ai.lbound+17)= 17.7;
    *(ai.lbound+18)= 18.0;
    *(ai.lbound+19)= 20.0;
  }
  if(!strncmp(prj,"engo",4))
  {
    *(ai.lbound+ 0)= -0.5;
    *(ai.lbound+ 1)=  1.0;
    *(ai.lbound+ 2)=  3.0;
    *(ai.lbound+ 3)= 12.0;
    *(ai.lbound+ 4)= 12.8;
    *(ai.lbound+ 5)= 13.6;
    *(ai.lbound+ 6)= 14.4;
    *(ai.lbound+ 7)= 15.2;
    *(ai.lbound+ 8)= 16.0;
    *(ai.lbound+ 9)= 16.8;
    *(ai.lbound+10)= 17.6;
    *(ai.lbound+11)= 18.4;
    *(ai.lbound+12)= 19.8;
    *(ai.lbound+13)= 25.7;
    *(ai.lbound+14)= 26.2;
    *(ai.lbound+15)= 26.7;
    *(ai.lbound+16)= 27.2;
	*(ai.lbound+17)= 27.7;
    *(ai.lbound+18)= 28.0;
    *(ai.lbound+19)= 31.0;
  }
  if(!strncmp(prj,"oos",3))
  {
    *(ai.lbound+ 0)= -0.5;
    *(ai.lbound+ 1)=  1.0;
    *(ai.lbound+ 2)=  4.0;
    *(ai.lbound+ 3)=  6.0;
    *(ai.lbound+ 4)= 12.8;
    *(ai.lbound+ 5)= 13.6;
    *(ai.lbound+ 6)= 14.4;
    *(ai.lbound+ 7)= 15.2;
    *(ai.lbound+ 8)= 16.0;
    *(ai.lbound+ 9)= 16.8;
    *(ai.lbound+10)= 17.6;
    *(ai.lbound+11)= 18.4;
    *(ai.lbound+12)= 19.8;
    *(ai.lbound+13)= 25.7;
    *(ai.lbound+14)= 26.2;
    *(ai.lbound+15)= 26.7;
    *(ai.lbound+16)= 27.2;
    *(ai.lbound+17)= 27.7;
    *(ai.lbound+18)= 28.0;
    *(ai.lbound+19)= 31.0;
  }
  if(!strncmp(prj,"ama",3))
  {
    *(ai.lbound+ 0)= -0.5;
    *(ai.lbound+ 1)=  1.0;
    *(ai.lbound+ 2)=  4.0;
    *(ai.lbound+ 3)=  6.0;
    *(ai.lbound+ 4)= 12.8;
    *(ai.lbound+ 5)= 13.6;
    *(ai.lbound+ 6)= 14.4;
    *(ai.lbound+ 7)= 15.2;
    *(ai.lbound+ 8)= 16.0;
    *(ai.lbound+ 9)= 16.8;
    *(ai.lbound+10)= 17.6;
    *(ai.lbound+11)= 18.4;
    *(ai.lbound+12)= 19.8;
    *(ai.lbound+13)= 25.7;
    *(ai.lbound+14)= 26.2;
	*(ai.lbound+15)= 26.7;
    *(ai.lbound+16)= 27.2;
    *(ai.lbound+17)= 27.7;
    *(ai.lbound+18)= 28.0;
    *(ai.lbound+19)= 31.0;
  }
  if(!strncmp(prj,"koba",4))
  {
    *(ai.lbound+ 0)= -0.5;
    *(ai.lbound+ 1)=  1.0;
    *(ai.lbound+ 2)=  2.0;
    *(ai.lbound+ 3)=  3.5;
    *(ai.lbound+ 4)=  5.0;
    *(ai.lbound+ 5)=  6.5;
    *(ai.lbound+ 6)=  7.5;
    *(ai.lbound+ 7)=  9.5;
    *(ai.lbound+ 8)= 16.0;
    *(ai.lbound+ 9)= 16.8;
    *(ai.lbound+10)= 17.6;
    *(ai.lbound+11)= 18.4;
    *(ai.lbound+12)= 19.8;
    *(ai.lbound+13)= 25.7;
    *(ai.lbound+14)= 26.2;
    *(ai.lbound+15)= 26.7;
    *(ai.lbound+16)= 27.2;
    *(ai.lbound+17)= 27.7;
    *(ai.lbound+18)= 28.0;
    *(ai.lbound+19)= 31.0;
  }
  if(!strncmp(prj,"nade",4))
  {
/*
    *(ai.lbound+ 0)= -5.0;
    *(ai.lbound+ 1)= -0.1;
    *(ai.lbound+ 2)=  0.1;
    *(ai.lbound+ 3)=  3.6;
    *(ai.lbound+ 4)=  5.6;
    *(ai.lbound+ 5)=  9.1;
    *(ai.lbound+ 6)= 16.5;
*/
    *(ai.lbound+ 0)= -0.1;
    *(ai.lbound+ 1)=  0.1;
    *(ai.lbound+ 2)=  3.6;
    *(ai.lbound+ 3)=  5.6;
	*(ai.lbound+ 4)=  9.1;
    *(ai.lbound+ 5)= 16.5;
    *(ai.lbound+ 6)= 17.5;
    *(ai.lbound+ 7)= 19.5;
    *(ai.lbound+ 8)= 26.0;
    *(ai.lbound+ 9)= 26.8;
    *(ai.lbound+10)= 27.6;
    *(ai.lbound+11)= 28.4;
    *(ai.lbound+12)= 29.8;
    *(ai.lbound+13)= 35.7;
    *(ai.lbound+14)= 36.2;
    *(ai.lbound+15)= 36.7;
    *(ai.lbound+16)= 37.2;
    *(ai.lbound+17)= 37.7;
    *(ai.lbound+18)= 38.0;
    *(ai.lbound+19)= 41.0;
  }
  if(!strncmp(prj,"karu",4))
  {
    *(ai.lbound+ 0)= -0.5;
    *(ai.lbound+ 1)=  1.0;
    *(ai.lbound+ 2)=  4.0;
    *(ai.lbound+ 3)= 13.5;
    *(ai.lbound+ 4)= 15.0;
    *(ai.lbound+ 5)= 16.5;
    *(ai.lbound+ 6)= 17.5;
    *(ai.lbound+ 7)= 19.5;
    *(ai.lbound+ 8)= 26.0;
    *(ai.lbound+ 9)= 26.8;
    *(ai.lbound+10)= 27.6;
    *(ai.lbound+11)= 28.4;
    *(ai.lbound+12)= 29.8;
    *(ai.lbound+13)= 35.7;
    *(ai.lbound+14)= 36.2;
    *(ai.lbound+15)= 36.7;
    *(ai.lbound+16)= 37.2;
    *(ai.lbound+17)= 37.7;
    *(ai.lbound+18)= 38.0;
    *(ai.lbound+19)= 41.0;
  }
  if(!strncmp(prj,"ashi",4))
  {
    *(ai.lbound+ 0)= -0.5;
    *(ai.lbound+ 1)=  1.0;
	*(ai.lbound+ 2)=  4.0;
    *(ai.lbound+ 3)= 20.0;
    *(ai.lbound+ 4)= 25.0;
    *(ai.lbound+ 5)= 26.5;
    *(ai.lbound+ 6)= 27.5;
    *(ai.lbound+ 7)= 29.5;
    *(ai.lbound+ 8)= 36.0;
    *(ai.lbound+ 9)= 36.8;
    *(ai.lbound+10)= 37.6;
    *(ai.lbound+11)= 38.4;
    *(ai.lbound+12)= 39.8;
    *(ai.lbound+13)= 45.7;
    *(ai.lbound+14)= 46.2;
    *(ai.lbound+15)= 46.7;
    *(ai.lbound+16)= 47.2;
    *(ai.lbound+17)= 47.7;
    *(ai.lbound+18)= 48.0;
    *(ai.lbound+19)= 51.0;
  }
  if(!strncmp(prj,"izu",3))
  {
    *(ai.lbound+ 0)= -0.5;
    *(ai.lbound+ 1)=  0.1;
    *(ai.lbound+ 2)=  0.5;
    *(ai.lbound+ 3)=  0.8;
    *(ai.lbound+ 4)=  1.2;
    *(ai.lbound+ 5)=  1.6;
    *(ai.lbound+ 6)=  1.9;
    *(ai.lbound+ 7)=  2.3;
    *(ai.lbound+ 8)=  2.6;
    *(ai.lbound+ 9)=  2.9;
    *(ai.lbound+10)=  3.4;
    *(ai.lbound+11)=  3.7;
    *(ai.lbound+12)=  4.0;
    *(ai.lbound+13)= 39.8;
    *(ai.lbound+14)= 45.7;
    *(ai.lbound+15)= 46.2;
    *(ai.lbound+16)= 46.7;
    *(ai.lbound+17)= 47.2;
    *(ai.lbound+18)= 47.7;
    *(ai.lbound+19)= 48.0;
  }
  if(!strncmp(prj,"eifu",4))
  {
	*(ai.lbound+ 0)= -0.5;
    *(ai.lbound+ 1)=  1.0;
    *(ai.lbound+ 2)=  4.0;
    *(ai.lbound+ 3)=  7.0;
    *(ai.lbound+ 4)= 25.0;
    *(ai.lbound+ 5)= 26.5;
    *(ai.lbound+ 6)= 27.5;
    *(ai.lbound+ 7)= 29.5;
    *(ai.lbound+ 8)= 36.0;
    *(ai.lbound+ 9)= 36.8;
    *(ai.lbound+10)= 37.6;
    *(ai.lbound+11)= 38.4;
    *(ai.lbound+12)= 39.8;
    *(ai.lbound+13)= 45.7;
    *(ai.lbound+14)= 46.2;
    *(ai.lbound+15)= 46.7;
    *(ai.lbound+16)= 47.2;
    *(ai.lbound+17)= 47.7;
    *(ai.lbound+18)= 48.0;
    *(ai.lbound+19)= 51.0;
  }
  if(!strncmp(prj,"gan",3))
  {
    *(ai.lbound+ 0)= -0.5;
    *(ai.lbound+ 1)=  1.0;
    *(ai.lbound+ 2)=  7.0;
    *(ai.lbound+ 3)= 11.0;
    *(ai.lbound+ 4)= 14.0;
    *(ai.lbound+ 5)= 17.0;
    *(ai.lbound+ 6)= 20.0;
    *(ai.lbound+ 7)= 23.0;
    *(ai.lbound+ 8)= 26.0;
    *(ai.lbound+ 9)= 29.0;
    *(ai.lbound+10)= 32.0;
    *(ai.lbound+11)= 36.0;
    *(ai.lbound+12)= 39.0;
    *(ai.lbound+13)= 42.0;
    *(ai.lbound+14)= 45.0;
    *(ai.lbound+15)= 48.0;
    *(ai.lbound+16)= 52.0;
    *(ai.lbound+17)= 55.0;
    *(ai.lbound+18)= 58.0;
    *(ai.lbound+19)= 90.0;
  }
  if(!strncmp(prj,"kiyo",4))
  {
    *(ai.lbound+ 0)= -0.5;
    *(ai.lbound+ 1)=  1.0;
    *(ai.lbound+ 2)=  6.0;
    *(ai.lbound+ 3)= 20.0;
    *(ai.lbound+ 4)= 25.0;
    *(ai.lbound+ 5)= 26.5;
    *(ai.lbound+ 6)= 27.5;
    *(ai.lbound+ 7)= 29.5;
    *(ai.lbound+ 8)= 36.0;
    *(ai.lbound+ 9)= 36.8;
    *(ai.lbound+10)= 37.6;
    *(ai.lbound+11)= 38.4;
    *(ai.lbound+12)= 39.8;
    *(ai.lbound+13)= 45.7;
    *(ai.lbound+14)= 46.2;
    *(ai.lbound+15)= 46.7;
    *(ai.lbound+16)= 47.2;
    *(ai.lbound+17)= 47.7;
    *(ai.lbound+18)= 48.0;
    *(ai.lbound+19)= 51.0;
  }
  if(!strncmp(prj,"miya",4))
  {
    *(ai.lbound+ 0)= -0.5;
    *(ai.lbound+ 1)=  1.0;
    *(ai.lbound+ 2)=  4.0;
    *(ai.lbound+ 3)=  6.0;
    *(ai.lbound+ 4)= 25.0;
    *(ai.lbound+ 5)= 26.5;
    *(ai.lbound+ 6)= 27.5;
    *(ai.lbound+ 7)= 29.5;
    *(ai.lbound+ 8)= 36.0;
    *(ai.lbound+ 9)= 36.8;
    *(ai.lbound+10)= 37.6;
    *(ai.lbound+11)= 38.4;
    *(ai.lbound+12)= 39.8;
    *(ai.lbound+13)= 45.7;
    *(ai.lbound+14)= 46.2;
    *(ai.lbound+15)= 46.7;
    *(ai.lbound+16)= 47.2;
    *(ai.lbound+17)= 47.7;
    *(ai.lbound+18)= 48.0;
	*(ai.lbound+19)= 51.0;
  }
  if(!strncmp(prj,"banga",5))
  {
    *(ai.lbound+ 0)= -0.5;
    *(ai.lbound+ 1)=  0.1;
    *(ai.lbound+ 2)=  0.4;
	*(ai.lbound+ 3)=  0.8;
    *(ai.lbound+ 4)=  1.1;
    *(ai.lbound+ 5)=  1.5;
    *(ai.lbound+ 6)=  1.8;
    *(ai.lbound+ 7)=  2.2;
    *(ai.lbound+ 8)=  2.5;
    *(ai.lbound+ 9)=  2.9;
    *(ai.lbound+10)=  3.2;
    *(ai.lbound+11)=  3.6;
    *(ai.lbound+12)=  3.9;
    *(ai.lbound+13)=  4.3;
    *(ai.lbound+14)=  4.6;
    *(ai.lbound+15)=  5.0;
    *(ai.lbound+16)=  5.3;
    *(ai.lbound+17)=  5.7;
    *(ai.lbound+18)=  6.1;
    *(ai.lbound+19)= 10.0;
  }
  if(!strncmp(prj,"tep",3))
  {
    *(ai.lbound+ 0)= -0.5;
    *(ai.lbound+ 1)=  1.7;
    *(ai.lbound+ 2)=  5.0;
    *(ai.lbound+ 3)=  6.0;
    *(ai.lbound+ 4)=  7.5;
    *(ai.lbound+ 5)= 10.0;
    *(ai.lbound+ 6)= 11.8;
    *(ai.lbound+ 7)= 12.2;
    *(ai.lbound+ 8)= 12.5;
    *(ai.lbound+ 9)= 12.9;
    *(ai.lbound+10)= 13.2;
    *(ai.lbound+11)= 13.6;
    *(ai.lbound+12)= 13.9;
    *(ai.lbound+13)= 14.3;
    *(ai.lbound+14)= 14.6;
    *(ai.lbound+15)= 15.0;
    *(ai.lbound+16)= 15.3;
	*(ai.lbound+17)= 15.7;
    *(ai.lbound+18)= 16.1;
    *(ai.lbound+19)= 20.0;
  }
  if(!strncmp(prj,"koshi",5))
  {
    *(ai.lbound+ 0)= -0.5;
    *(ai.lbound+ 1)=  1.0;
    *(ai.lbound+ 2)=  5.0;
    *(ai.lbound+ 3)=  7.0;
    *(ai.lbound+ 4)= 15.0;
    *(ai.lbound+ 5)= 18.0;
    *(ai.lbound+ 6)= 19.0;
    *(ai.lbound+ 7)= 20.0;
    *(ai.lbound+ 8)= 21.0;
    *(ai.lbound+ 9)= 22.0;
    *(ai.lbound+10)= 23.0;
    *(ai.lbound+11)= 24.0;
    *(ai.lbound+12)= 25.0;
    *(ai.lbound+13)= 26.0;
    *(ai.lbound+14)= 27.0;
    *(ai.lbound+15)= 28.0;
    *(ai.lbound+16)= 29.0;
    *(ai.lbound+17)= 30.0;
    *(ai.lbound+18)= 31.0;
    *(ai.lbound+19)= 32.0;
  }
  if(!strncmp(prj,"toyo",4))
  {
    *(ai.lbound+ 0)= -0.1;
    *(ai.lbound+ 1)=  4.05;
    *(ai.lbound+ 2)=  7.0;
    *(ai.lbound+ 3)= 30.0;
    *(ai.lbound+ 4)= 35.0;
    *(ai.lbound+ 5)= 38.0;
    *(ai.lbound+ 6)= 39.0;
    *(ai.lbound+ 7)= 40.0;
    *(ai.lbound+ 8)= 41.0;
    *(ai.lbound+ 9)= 42.0;
    *(ai.lbound+10)= 43.0;
    *(ai.lbound+11)= 44.0;
    *(ai.lbound+12)= 45.0;
    *(ai.lbound+13)= 46.0;
    *(ai.lbound+14)= 47.0;
	*(ai.lbound+15)= 48.0;
    *(ai.lbound+16)= 49.0;
    *(ai.lbound+17)= 50.0;
    *(ai.lbound+18)= 51.0;
    *(ai.lbound+19)= 52.0;
  }
  if(!strncmp(prj,"hakone",6))
  {
    *(ai.lbound+ 0)= -0.1;
    *(ai.lbound+ 1)=  0.1;
    *(ai.lbound+ 2)=  5.0;
    *(ai.lbound+ 3)= 10.0;
    *(ai.lbound+ 4)= 16.0;
    *(ai.lbound+ 5)= 17.5;
    *(ai.lbound+ 6)= 35.0;
    *(ai.lbound+ 7)= 40.0;
    *(ai.lbound+ 8)= 41.0;
    *(ai.lbound+ 9)= 42.0;
    *(ai.lbound+10)= 43.0;
    *(ai.lbound+11)= 44.0;
    *(ai.lbound+12)= 45.0;
    *(ai.lbound+13)= 46.0;
    *(ai.lbound+14)= 47.0;
    *(ai.lbound+15)= 48.0;
    *(ai.lbound+16)= 49.0;
    *(ai.lbound+17)= 50.0;
    *(ai.lbound+18)= 51.0;
    *(ai.lbound+19)= 52.0;
  }
  if(!strncmp(prj,"sajigym",7))
  {
    *(ai.lbound+ 0)= -0.1;
    *(ai.lbound+ 1)=  0.1;
    *(ai.lbound+ 2)= 25.0;
    *(ai.lbound+ 3)= 30.0;
    *(ai.lbound+ 4)= 36.0;
    *(ai.lbound+ 5)= 37.5;
    *(ai.lbound+ 6)= 45.0;
    *(ai.lbound+ 7)= 50.0;
    *(ai.lbound+ 8)= 51.0;
    *(ai.lbound+ 9)= 52.0;
    *(ai.lbound+10)= 53.0;
    *(ai.lbound+11)= 54.0;
    *(ai.lbound+12)= 55.0;
	*(ai.lbound+13)= 56.0;
    *(ai.lbound+14)= 57.0;
    *(ai.lbound+15)= 58.0;
    *(ai.lbound+16)= 59.0;
    *(ai.lbound+17)= 60.0;
    *(ai.lbound+18)= 61.0;
    *(ai.lbound+19)= 62.0;
  }
  if(!strncmp(prj,"toko",4))
  {
    *(ai.lbound+ 0)= -1.0;
    *(ai.lbound+ 1)=  1.0;
    *(ai.lbound+ 2)=  5.0;
    *(ai.lbound+ 3)=  8.0;
    *(ai.lbound+ 4)= 12.0;
    *(ai.lbound+ 5)= 38.0;
    *(ai.lbound+ 6)= 39.0;
    *(ai.lbound+ 7)= 40.0;
    *(ai.lbound+ 8)= 41.0;
    *(ai.lbound+ 9)= 42.0;
    *(ai.lbound+10)= 43.0;
    *(ai.lbound+11)= 44.0;
    *(ai.lbound+12)= 45.0;
    *(ai.lbound+13)= 46.0;
    *(ai.lbound+14)= 47.0;
    *(ai.lbound+15)= 48.0;
    *(ai.lbound+16)= 49.0;
    *(ai.lbound+17)= 50.0;
    *(ai.lbound+18)= 51.0;
    *(ai.lbound+19)= 52.0;
  }
  if(!strncmp(prj,"maebasi",7))
  {
    *(ai.lbound+ 0)= -1.0;
    *(ai.lbound+ 1)=  1.0;
    *(ai.lbound+ 2)=  5.0;
    *(ai.lbound+ 3)=  8.0;
    *(ai.lbound+ 4)= 12.0;
    *(ai.lbound+ 5)= 38.0;
    *(ai.lbound+ 6)= 39.0;
    *(ai.lbound+ 7)= 40.0;
    *(ai.lbound+ 8)= 41.0;
    *(ai.lbound+ 9)= 42.0;
    *(ai.lbound+10)= 43.0;
	*(ai.lbound+11)= 44.0;
    *(ai.lbound+12)= 45.0;
    *(ai.lbound+13)= 46.0;
    *(ai.lbound+14)= 47.0;
    *(ai.lbound+15)= 48.0;
    *(ai.lbound+16)= 49.0;
    *(ai.lbound+17)= 50.0;
    *(ai.lbound+18)= 51.0;
    *(ai.lbound+19)= 52.0;
  }
  if(!strncmp(prj,"dazai",5))
  {
    *(ai.lbound+ 0)= -1.0;
    *(ai.lbound+ 1)=  0.6;
    *(ai.lbound+ 2)=  5.1;
    *(ai.lbound+ 3)= 10.0;
    *(ai.lbound+ 4)= 12.0;
    *(ai.lbound+ 5)= 38.0;
    *(ai.lbound+ 6)= 39.0;
    *(ai.lbound+ 7)= 40.0;
    *(ai.lbound+ 8)= 41.0;
    *(ai.lbound+ 9)= 42.0;
    *(ai.lbound+10)= 43.0;
    *(ai.lbound+11)= 44.0;
    *(ai.lbound+12)= 45.0;
    *(ai.lbound+13)= 46.0;
    *(ai.lbound+14)= 47.0;
    *(ai.lbound+15)= 48.0;
    *(ai.lbound+16)= 49.0;
    *(ai.lbound+17)= 50.0;
    *(ai.lbound+18)= 51.0;
    *(ai.lbound+19)= 52.0;
  }
  if(!strncmp(prj,"akita",5))
  {
    *(ai.lbound+ 0)= -1.0;
    *(ai.lbound+ 1)=  0.5;
    *(ai.lbound+ 2)=  5.0;
    *(ai.lbound+ 3)= 10.0;
    *(ai.lbound+ 4)= 12.0;
    *(ai.lbound+ 5)= 38.0;
    *(ai.lbound+ 6)= 39.0;
    *(ai.lbound+ 7)= 40.0;
    *(ai.lbound+ 8)= 41.0;
	*(ai.lbound+ 9)= 42.0;
    *(ai.lbound+10)= 43.0;
    *(ai.lbound+11)= 44.0;
    *(ai.lbound+12)= 45.0;
    *(ai.lbound+13)= 46.0;
    *(ai.lbound+14)= 47.0;
    *(ai.lbound+15)= 48.0;
    *(ai.lbound+16)= 49.0;
    *(ai.lbound+17)= 50.0;
    *(ai.lbound+18)= 51.0;
    *(ai.lbound+19)= 52.0;
  }
  if(!strncmp(prj,"agc",3))
  {
    *(ai.lbound+ 0)= -1.0;
    *(ai.lbound+ 1)=  0.01;
    *(ai.lbound+ 2)= 10.0;
    *(ai.lbound+ 3)= 11.0;
    *(ai.lbound+ 4)= 12.0;
    *(ai.lbound+ 5)= 38.0;
    *(ai.lbound+ 6)= 39.0;
    *(ai.lbound+ 7)= 40.0;
    *(ai.lbound+ 8)= 41.0;
    *(ai.lbound+ 9)= 42.0;
    *(ai.lbound+10)= 43.0;
    *(ai.lbound+11)= 44.0;
    *(ai.lbound+12)= 45.0;
    *(ai.lbound+13)= 46.0;
    *(ai.lbound+14)= 47.0;
    *(ai.lbound+15)= 48.0;
    *(ai.lbound+16)= 49.0;
    *(ai.lbound+17)= 50.0;
    *(ai.lbound+18)= 51.0;
    *(ai.lbound+19)= 52.0;
  }
  if(!strncmp(prj,"minna",5))
  {
    *(ai.lbound+ 0)= -1.0;
    *(ai.lbound+ 1)=  0.01;
    *(ai.lbound+ 2)= 20.0;
    *(ai.lbound+ 3)= 21.0;
    *(ai.lbound+ 4)= 22.0;
    *(ai.lbound+ 5)= 38.0;
    *(ai.lbound+ 6)= 39.0;
	*(ai.lbound+ 7)= 40.0;
    *(ai.lbound+ 8)= 41.0;
    *(ai.lbound+ 9)= 42.0;
    *(ai.lbound+10)= 43.0;
    *(ai.lbound+11)= 44.0;
    *(ai.lbound+12)= 45.0;
    *(ai.lbound+13)= 46.0;
    *(ai.lbound+14)= 47.0;
    *(ai.lbound+15)= 48.0;
    *(ai.lbound+16)= 49.0;
    *(ai.lbound+17)= 50.0;
	*(ai.lbound+18)= 51.0;
	*(ai.lbound+19)= 52.0;
  }
  if(!strncmp(prj,"riku",4))
  {
	*(ai.lbound+ 0)= -1.0;
	*(ai.lbound+ 1)=  1.5;
	*(ai.lbound+ 2)=  3.2;
	*(ai.lbound+ 3)=  4.0;
	*(ai.lbound+ 4)=  5.6;
	*(ai.lbound+ 5)= 30.0;
	*(ai.lbound+ 6)= 39.0;
	*(ai.lbound+ 7)= 40.0;
	*(ai.lbound+ 8)= 41.0;
	*(ai.lbound+ 9)= 42.0;
	*(ai.lbound+10)= 43.0;
	*(ai.lbound+11)= 44.0;
	*(ai.lbound+12)= 45.0;
	*(ai.lbound+13)= 46.0;
	*(ai.lbound+14)= 47.0;
	*(ai.lbound+15)= 48.0;
	*(ai.lbound+16)= 49.0;
	*(ai.lbound+17)= 50.0;
	*(ai.lbound+18)= 51.0;
	*(ai.lbound+19)= 52.0;
  }
  if(!strncmp(prj,"sichi",5))
  {
	*(ai.lbound+ 0)= -1.0;
	*(ai.lbound+ 1)=  1.0;
	*(ai.lbound+ 2)=  5.0;
	*(ai.lbound+ 3)=  9.0;
	*(ai.lbound+ 4)= 15.0;
	*(ai.lbound+ 5)= 30.0;
	*(ai.lbound+ 6)= 39.0;
	*(ai.lbound+ 7)= 40.0;
	*(ai.lbound+ 8)= 41.0;
	*(ai.lbound+ 9)= 42.0;
	*(ai.lbound+10)= 43.0;
	*(ai.lbound+11)= 44.0;
	*(ai.lbound+12)= 45.0;
	*(ai.lbound+13)= 46.0;
	*(ai.lbound+14)= 47.0;
	*(ai.lbound+15)= 48.0;
	*(ai.lbound+16)= 49.0;
	*(ai.lbound+17)= 50.0;
	*(ai.lbound+18)= 51.0;
	*(ai.lbound+19)= 52.0;
  }
  if(!strncmp(prj,"ebis",4))
  {
	*(ai.lbound+ 0)= -1.0;
	*(ai.lbound+ 1)=  1.0;
	*(ai.lbound+ 2)=  4.0;
	*(ai.lbound+ 3)= 15.0;
	*(ai.lbound+ 4)= 20.0;
	*(ai.lbound+ 5)= 30.0;
	*(ai.lbound+ 6)= 39.0;
	*(ai.lbound+ 7)= 40.0;
	*(ai.lbound+ 8)= 41.0;
	*(ai.lbound+ 9)= 42.0;
	*(ai.lbound+10)= 43.0;
	*(ai.lbound+11)= 44.0;
	*(ai.lbound+12)= 45.0;
	*(ai.lbound+13)= 46.0;
	*(ai.lbound+14)= 47.0;
	*(ai.lbound+15)= 48.0;
	*(ai.lbound+16)= 49.0;
	*(ai.lbound+17)= 50.0;
	*(ai.lbound+18)= 51.0;
	*(ai.lbound+19)= 52.0;
  }
  if(!strncmp(prj,"stick",5))
  {
	*(ai.lbound+ 0)= -1.0;
	*(ai.lbound+ 1)=  0.1;
	*(ai.lbound+ 2)= 14.0;
	*(ai.lbound+ 3)= 15.0;
	*(ai.lbound+ 4)= 20.0;
	*(ai.lbound+ 5)= 30.0;
	*(ai.lbound+ 6)= 39.0;
	*(ai.lbound+ 7)= 40.0;
	*(ai.lbound+ 8)= 41.0;
	*(ai.lbound+ 9)= 42.0;
	*(ai.lbound+10)= 43.0;
	*(ai.lbound+11)= 44.0;
	*(ai.lbound+12)= 45.0;
	*(ai.lbound+13)= 46.0;
	*(ai.lbound+14)= 47.0;
	*(ai.lbound+15)= 48.0;
	*(ai.lbound+16)= 49.0;
	*(ai.lbound+17)= 50.0;
	*(ai.lbound+18)= 51.0;
	*(ai.lbound+19)= 52.0;
  }
  if(!strncmp(prj,"sydney",6))
  {
	*(ai.lbound+ 0)=  -1.0;
	*(ai.lbound+ 1)=   0.1;
	*(ai.lbound+ 2)=  10.0;
	*(ai.lbound+ 3)=  20.0;
	*(ai.lbound+ 4)=  30.0;
	*(ai.lbound+ 5)=  40.0;
	*(ai.lbound+ 6)=  50.0;
	*(ai.lbound+ 7)=  60.0;
	*(ai.lbound+ 8)=  70.0;
	*(ai.lbound+ 9)=  80.0;
	*(ai.lbound+10)= 100.0;
	*(ai.lbound+11)= 120.0;
	*(ai.lbound+12)= 170.0;
	*(ai.lbound+13)= 180.0;
	*(ai.lbound+14)= 190.0;
	*(ai.lbound+15)= 200.0;
	*(ai.lbound+16)= 210.0;
	*(ai.lbound+17)= 220.0;
	*(ai.lbound+18)= 240.0;
	*(ai.lbound+19)= 260.0;
  }
  if(!strncmp(prj,"ghouse",6))     /*201113ghouse*/
  {
	*(ai.lbound+ 0)= -1.0;
	*(ai.lbound+ 1)=  0.5;
	*(ai.lbound+ 2)= 10.0;
	*(ai.lbound+ 3)= 15.0;
	*(ai.lbound+ 4)= 20.0;
	*(ai.lbound+ 5)= 30.0;
	*(ai.lbound+ 6)= 39.0;
	*(ai.lbound+ 7)= 40.0;
	*(ai.lbound+ 8)= 41.0;
	*(ai.lbound+ 9)= 42.0;
	*(ai.lbound+10)= 43.0;
	*(ai.lbound+11)= 44.0;
	*(ai.lbound+12)= 45.0;
	*(ai.lbound+13)= 46.0;
	*(ai.lbound+14)= 47.0;
	*(ai.lbound+15)= 48.0;
	*(ai.lbound+16)= 49.0;
	*(ai.lbound+17)= 50.0;
	*(ai.lbound+18)= 51.0;
	*(ai.lbound+19)= 52.0;
  }
  if(!strncmp(prj,"optimization",12))   //UJIOKA FOR OPTIMIZATION
  {
	*(ai.lbound+ 0)=  -1.0;
	*(ai.lbound+ 1)=   0.1;
	*(ai.lbound+ 2)=  30.0;
	*(ai.lbound+ 3)=  35.0;
	*(ai.lbound+ 4)=  40.0;
	*(ai.lbound+ 5)=  45.0;
	*(ai.lbound+ 6)=  50.0;
	*(ai.lbound+ 7)=  60.0;
	*(ai.lbound+ 8)=  70.0;
	*(ai.lbound+ 9)=  80.0;
	*(ai.lbound+10)= 100.0;
	*(ai.lbound+11)= 120.0;
	*(ai.lbound+12)= 170.0;
	*(ai.lbound+13)= 180.0;
	*(ai.lbound+14)= 190.0;
	*(ai.lbound+15)= 200.0;
	*(ai.lbound+16)= 210.0;
	*(ai.lbound+17)= 220.0;
	*(ai.lbound+18)= 240.0;
	*(ai.lbound+19)= 260.0;
  }

//add by fukushima 20140618
  if (org->ai.nfloor > 0)
  {
	ai.nfloor=org->ai.nfloor;
	ai.lbound=org->ai.lbound;
  }
///////////////////////////

  /*FLOOR WEIGHT*/
  for(i=0;i<(org->nnode);i++)
  {
	gn=org->nodes+i;

	/*if((gn->d[GZ])>ai.hmax) ai.hmax=gn->d[GZ];*/

	for(j=0;j<ai.nfloor;j++)
    {
      if((gn->d[GZ])>=*(ai.lbound+j) &&
         (gn->d[GZ])< *(ai.lbound+j+1))
      {
        (*(ai.nnode+j))++;
        (*(ai.levels+j))+=gn->d[GZ];
        (*(ai.wi+j))+=(nload+i)->w[WEIGHTEQ];
      }
    }
  }

  /*FLOOR LEVEL*/
  ai.hmax=0.0;
  for(i=0;i<ai.nfloor;i++)
  {
    if(*(ai.nnode+i)>0)
    {
      (*(ai.levels+i))/=(double)(*(ai.nnode+i)); /*MEAN LEVEL*/
    }
    if((*(ai.levels+i))<0.0) (*(ai.levels+i))=0.0; //mihara for basement floor
	if(ai.hmax<(*(ai.levels+i))) ai.hmax=*(ai.levels+i);

    for(j=i;j<ai.nfloor;j++)
    {
      (*(ai.Wi+i))+=(*(ai.wi+j));
    }
  }

  ai.hmax-=(*(ai.levels+0)); /*HEIGHT*/
  ai.T1*=ai.hmax; /*PERIOD*/

if(!strncmp(prj,"stick",5)) ai.hmax=3.840;

  if(ai.T1<ai.Tc) ai.Rt=1.0;
  else if(ai.T1<2.0*(ai.Tc)) ai.Rt=1.0-0.2*(ai.T1/ai.Tc-1.0)
                                          *(ai.T1/ai.Tc-1.0);
  else ai.Rt=1.6*(ai.Tc)/(ai.T1);

//  ai.Cf=0.5*(ai.Z)*(ai.Rt)*(ai.Co); /*FOUNDATION SHEAR*/
  ai.Cfx=0.5*(ai.Z)*(ai.Rt)*(ai.Cox); /*FOUNDATION SHEAR*/
  ai.Cfy=0.5*(ai.Z)*(ai.Rt)*(ai.Coy);

  for(i=1;i<ai.nfloor;i++)
  {
    if(*(ai.Wi+1)!=0.0) alpha=(*(ai.Wi+i))/(*(ai.Wi+1));
    else alpha=0.0;

    if(alpha!=0.0)
    {
      (*(ai.Ai+i))=1.0
                   +(1.0/sqrt(alpha)-alpha)
                   *2.0*(ai.T1)/(1.0+3.0*(ai.T1));
    }
    else (*(ai.Ai+i))=1.0;

    (*(ai.Cix+i))=(*(ai.Ai+i))*(ai.Z)*(ai.Rt)*(ai.Cox);
    (*(ai.Ciy+i))=(*(ai.Ai+i))*(ai.Z)*(ai.Rt)*(ai.Coy);
/*
    if(!strncmp(prj,"hakone",6))
    {
      (*(ai.Ci+1))=0.122*1.25;
      (*(ai.Ci+2))=0.132*1.25;
      (*(ai.Ci+3))=0.126*1.25;
      (*(ai.Ci+4))=0.130*1.25;
      (*(ai.Ci+5))=0.131*1.25;
    }
*/
    (*(ai.Qix+i))=(*(ai.Cix+i))*(*(ai.Wi+i));
    (*(ai.Qiy+i))=(*(ai.Ciy+i))*(*(ai.Wi+i));

    if(i>=2)
    {
      (*(ai.Hix+i-1))=(*(ai.Qix+i-1))-(*(ai.Qix+i));
      (*(ai.Hiy+i-1))=(*(ai.Qiy+i-1))-(*(ai.Qiy+i));

      if(*(ai.wi+i-1)!=0.0)
      {
        (*(ai.factsx+i-1))=(*(ai.Hix+i-1))/(*(ai.wi+i-1));
        (*(ai.factsy+i-1))=(*(ai.Hiy+i-1))/(*(ai.wi+i-1));
      }
      else
      {
        (*(ai.factsx+i-1))=0.0;
        (*(ai.factsy+i-1))=0.0;
      }
    }
  }
  (*(ai.Hix+(ai.nfloor-1)))=(*(ai.Qix+(ai.nfloor-1)));
  (*(ai.Hiy+(ai.nfloor-1)))=(*(ai.Qiy+(ai.nfloor-1)));

  if(*(ai.wi+(ai.nfloor-1))!=0.0)
  {
    (*(ai.factsx+(ai.nfloor-1)))=(*(ai.Hix+(ai.nfloor-1)))
                               /(*(ai.wi+(ai.nfloor-1)));
    (*(ai.factsy+(ai.nfloor-1)))=(*(ai.Hiy+(ai.nfloor-1)))
                               /(*(ai.wi+(ai.nfloor-1)));
  }
  else
  {
    (*(ai.factsx+(ai.nfloor-1)))=0.0;
    (*(ai.factsy+(ai.nfloor-1)))=0.0;
  }
  (*(ai.factsx+0))=ai.Cfx;
  (*(ai.factsy+0))=ai.Cfy;

  for(i=0;i<(org->nnode);i++)
  {
    gn=org->nodes+i;

    for(j=0;j<ai.nfloor;j++)
    {
      if((gn->d[GZ])>=*(ai.lbound+j) &&
         (gn->d[GZ])< *(ai.lbound+j+1))
      {
        ((nload+i)->w[WEIGHTEQX])=((nload+i)->w[WEIGHTEQ]);
        ((nload+i)->w[WEIGHTEQY])=((nload+i)->w[WEIGHTEQ]);
        ((nload+i)->w[WEIGHTEQX])*=(*(ai.factsx+j));
        ((nload+i)->w[WEIGHTEQY])*=(*(ai.factsy+j));
      }
    }
  }

  if(fout!=NULL)
  {
    fprintf(fout,"\nUnit Factor  =%.5f",globalunit);

    if(globalunit==SIUNIT)
    fprintf(fout," \"SI Units [kN]\"\n");
    else if(globalunit==1.0)
    fprintf(fout," \"Classic Units [tf]\"\n");
    else
    fprintf(fout," \"Unknown Units\"\n");

    fprintf(fout,"\n3.3 : Ai分布型地震荷重\n\n");
    fprintf(fout,"水平荷重は建築基準法施行令第88条および建設省告示1793号に従い、Ａi分布型の地震力とする。\n\n");

//    fprintf(fout,"Floors     n =%d\n",ai.nfloor);
//    fprintf(fout,"Height     H =%.3f\n",ai.hmax);

    fprintf(fout,"階数       　　　    n =%d\n",ai.nfloor);
	fprintf(fout,"高さ         　　　  H =%.3f\n",ai.hmax);

//    if(ai.hmax!=0.0) fprintf(fout,"Period     T1=%.3fH=%.3f\n",(ai.T1/ai.hmax),ai.T1);
//    else             fprintf(fout,"Period     T1       =%.3f\n",ai.T1);
    if(ai.hmax!=0.0) fprintf(fout,"１次固有周期         T1=%.3fH=%.3f\n",(ai.T1/ai.hmax),ai.T1);
    else             fprintf(fout,"１次固有周期         T1       =%.3f\n",ai.T1);
/*
    fprintf(fout,"Period     Tc=%.3f\n",ai.Tc);
    fprintf(fout,"           Rt=%.3f\n",ai.Rt);
    fprintf(fout,"           Z =%.3f\n",ai.Z);
    fprintf(fout,"Base Shear Co=%.3f\n",ai.Co);
    fprintf(fout,"Foundation Cf=%.3f\n",ai.Cf);
*/
    fprintf(fout,"地盤周期             Tc=%.3f\n",ai.Tc);
    fprintf(fout,"振動特性係数         Rt=%.3f\n",ai.Rt);
    fprintf(fout,"地域係数             Z =%.3f\n",ai.Z);
    if (ai.Cox==ai.Coy)
    {
      fprintf(fout,"標準層せん断力係数   Co=%.3f\n",ai.Cox);
      fprintf(fout,"基礎部分の震度       Cf=%.3f\n",ai.Cfx);
    }
    else
    {
      fprintf(fout,"標準層せん断力係数（X方向）   Co=%.3f\n",ai.Cox);
      fprintf(fout,"標準層せん断力係数（Y方向）   Co=%.3f\n",ai.Coy);
      fprintf(fout,"基礎部分の震度（X方向）       Cf=%.3f\n",ai.Cfx);
      fprintf(fout,"基礎部分の震度（Y方向）       Cf=%.3f\n",ai.Cfy);
    }

  }

#if 0
  if(globalunit==1.0)
  {
    fprintf(fout,"\n床用積載荷重による総重量   = %.3f [tf]\n",      //mihara
            globalunit*total[WEIGHTSLAB]);
    fprintf(fout,"柱梁用積載荷重による総重量   = %.3f [tf]\n",
            globalunit*total[WEIGHTFRAME]);
    fprintf(fout,"地震用積載荷重による総重量   = %.3f [tf]\n",
            globalunit*total[WEIGHTEQ]);
  }
  if(globalunit==SIUNIT)
  {
    fprintf(fout,"\n床用積載荷重による総重量   = %.3f [kN]\n",      //mihara
            globalunit*total[WEIGHTSLAB]);
	fprintf(fout,"柱梁用積載荷重による総重量   = %.3f [kN]\n",
            globalunit*total[WEIGHTFRAME]);
    fprintf(fout,"地震用積載荷重による総重量   = %.3f [kN]\n",
            globalunit*total[WEIGHTEQ]);
  }
#endif

  sprintf(ss1,"各階平均高さ      :");
  sprintf(ss2,"各階重量       wi :");
  sprintf(ss3,"        Wi = Σwi :");
  sprintf(ss4,"               Ai :           ");
  sprintf(ss5,"層せん断力係数 Ci :           ");
  sprintf(ss6,"層せん断力     Qi :           ");
  sprintf(ss7,"各階外力       Hi :           ");
  sprintf(ss8,"外力係数    Hi/wi :");

  for(i=0;i<ai.nfloor;i++)
  {
	sprintf(ss," %10.3f",*(ai.levels+i));
    strcat(ss1,ss);
    sprintf(ss," %10.3f",globalunit*(*(ai.wi+i)));
    strcat(ss2,ss);
    sprintf(ss," %10.3f",globalunit*(*(ai.Wi+i)));
    strcat(ss3,ss);
	if(i>0)
	{
	  sprintf(ss," %10.3f",*(ai.Ai+i));
      strcat(ss4,ss);
	  sprintf(ss," %10.3f",*(ai.Cix+i));
      strcat(ss5,ss);
      sprintf(ss," %10.3f",globalunit*(*(ai.Qix+i)));
	  strcat(ss6,ss);
	  sprintf(ss," %10.3f",globalunit*(*(ai.Hix+i)));
	  strcat(ss7,ss);
	}
	sprintf(ss," %10.3f",*(ai.factsx+i));
	strcat(ss8,ss);
  }
  /*sprintf(ss," Height=%.3f",ai.hmax);
  strcat(ss1,ss);*/

  sprintf(non,"\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n",
		  ss1,ss2,ss3,ss4,ss5,ss6,ss7,ss8);
#if 1
  sprintf(ss11,"各階平均高さ      :");
  sprintf(ss12,"各階重量       wi :");
  sprintf(ss13,"        Wi = Σwi :");
  sprintf(ss14,"               Ai :           ");
  sprintf(ss15,"層せん断力係数 Ci :           ");
  sprintf(ss16,"層せん断力     Qi :           ");
  sprintf(ss17,"各階外力       Hi :           ");
  sprintf(ss18,"外力係数    Hi/wi :");

  for(i=0;i<ai.nfloor;i++)
  {
	sprintf(ss," %10.3f",*(ai.levels+i));
    strcat(ss11,ss);
    sprintf(ss," %10.3f",globalunit*(*(ai.wi+i)));
    strcat(ss12,ss);
    sprintf(ss," %10.3f",globalunit*(*(ai.Wi+i)));
    strcat(ss13,ss);
	if(i>0)
	{
	  sprintf(ss," %10.3f",*(ai.Ai+i));
      strcat(ss14,ss);
	  sprintf(ss," %10.3f",*(ai.Ciy+i));
      strcat(ss15,ss);
      sprintf(ss," %10.3f",globalunit*(*(ai.Qiy+i)));
	  strcat(ss16,ss);
	  sprintf(ss," %10.3f",globalunit*(*(ai.Hiy+i)));
	  strcat(ss17,ss);
	}
	sprintf(ss," %10.3f",*(ai.factsy+i));
	strcat(ss18,ss);
  }
  /*sprintf(ss," Height=%.3f",ai.hmax);
  strcat(ss1,ss);*/

  sprintf(nonII,"\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n",
		  ss11,ss12,ss13,ss14,ss15,ss16,ss17,ss18);
#endif
  if(fout!=NULL) fprintf(fout,"%s",non);
  if(fout!=NULL) fprintf(fout,"%s",nonII);
  if(globalmessageflag==1) MessageBox(NULL,non,"Distribute",MB_OK);

  free(fload);
  for(i=0;i<ai.nfloor;i++)
  {
	free(*(ai.fnames+i));
  }
  free(ai.fnames);
  free(ai.nnode);
  free(ai.lbound);
  free(ai.levels);
  free(ai.wi);
  free(ai.Wi);
  free(ai.Ai);
  free(ai.Cix);
  free(ai.Ciy);
  free(ai.Qix);
  free(ai.Qiy);
  free(ai.Hix);
  free(ai.Hiy);
  free(ai.factsx);
  free(ai.factsy);

  return;
}/*weightdistribution*/

void sortdouble(double value[],int nvalue)
/*SORT DOUBLE AS KOUJUN BY "HEAP SORT".*/
{
  int i,j,ir,m;
  double rra;

  m=(int)(nvalue/2);
  ir=nvalue-1;

  for(;;)
  {
    if(m>0)
    {
      rra=value[--m];
    }
    else
    {
      rra=value[ir];
      value[ir]=value[0];
      if(--ir==0)
      {
        value[0]=rra;
        return;
      }
    }
    i=m;
    j=m*2+1;
    while(j<=ir)
    {
      if(j<ir && value[j]>value[j+1]) ++j;

      if(rra>value[j])
      {
        value[i]=value[j];
        i=j;
        j+=j+1;
      }
      else j=ir+1;
    }
    value[i]=rra;
  }
}/*sortdouble*/

void sortdoublepointer(double *value,int nvalue)
/*SORT DOUBLE AS KOUJUN BY "HEAP SORT".*/
{
  int i,j,ir,m;
  double rra;

  m=(int)(nvalue/2);
  ir=nvalue-1;

  for(;;)
  {
    if(m>0)
    {
      m--;
      rra=*(value+m);
    }
    else
    {
      rra=*(value+ir);
      *(value+ir)=*(value+0);
      if(--ir==0)
      {
        *(value+0)=rra;
        return;
      }
    }
    i=m;
    j=m*2+1;
    while(j<=ir)
    {
      if(j<ir && (*(value+j))>(*(value+j+1))) ++j;

      if(rra>(*(value+j)))
      {
        *(value+i)=*(value+j);
        i=j;
        j+=j+1;
      }
      else j=ir+1;
    }
    *(value+i)=rra;
  }
}/*sortdoublepointer*/

int quadraticequation(double c2,double c1,double c0,
                      double answer[])
/*QUADRATIC EQUATION C2X2+C1X+C0=0.*/
/*BASED ON "NUMERICAL RECIPES P.158".*/
/*RETURN:NUMBER OF REAL QUANTITIES.*/
{
  double D,Q,sign;

  if(c2==0.0) return 0;

  D=c1*c1-4.0*c2*c0;

  if(D>=0.0)
  {
    if(c1==0) sign=1.0;
    else      sign=c1/fabs(c1);

    Q=-0.5*(c1+sign*sqrt(D));

    answer[0]=Q/c2;
    answer[1]=c0/Q;

    return 2;
  }
  else return 0;
}/*quadraticequation*/

int cubicequation(double c3,double c2,double c1,double c0,
                  double answer[])
/*CUBIC EQUATION C3X3+C2X2+C1X+C0=0.*/
/*BASED ON "NUMERICAL RECIPES P.158".*/
/*RETURN:NUMBER OF REAL QUANTITIES.*/
{
  double Q,R,theta,value;

  if(c3==0.0) return 0;

  c2/=c3; c1/=c3; c0/=c3;

  Q=(c2*c2-3.0*c1)/9.0;
  R=(2.0*c2*c2*c2-9.0*c2*c1+27.0*c0)/54.0;

  if((Q*Q*Q-R*R)>=0.0)
  {
    theta=acos(R/sqrt(Q*Q*Q));

    answer[0]=-2.0*sqrt(Q)*cos(theta/3.0)-c2/3.0;
    answer[1]=-2.0*sqrt(Q)*cos((theta+2.0*PI)/3.0)-c2/3.0;
    answer[2]=-2.0*sqrt(Q)*cos((theta+4.0*PI)/3.0)-c2/3.0;

    return 3;
  }
  else if(R!=0.0)
  {
    value=pow((sqrt(R*R-Q*Q*Q)+fabs(R)),(1.0/3.0));

    answer[0]=-R/fabs(R)*(value+Q/value)-c2/3.0;

    return 1;
  }
  else return 0;
}/*cubicequation*/

void inputcadretomemory(FILE *ftext,struct arclmframe *af)
/*TRANSLATE CADRE INPUTFILE TEXT INTO MEMORY.*/
{
  char **data,str[256];
  char **library=NULL,*s;
  int i,j,ii,jj,n,flag;
  int rn,nbeam,nplate;
  long int offset;
  long int ncode,hcode,tcode,rcode;
  long int *soffsets;
  struct onode *rnode;
  double E,G,thick;
  double xmax=0.0,xmin=0.0,ymax=0.0,ymin=0.0,zmax=0.0,zmin=0.0;

  fseek(ftext,0L,SEEK_SET);

  flag=0;
  while(flag==0)
  {
    fgets(str,256,ftext);
    if(!strncmp(str,"Basic",5)) flag=1;
  }

  fgets(str,256,ftext);
  fgets(str,256,ftext);

  /*NODES*/
  data=fgetsbrk(ftext,&n);
  if(!strncmp(*(data+0),"Structrual",10))
  {
    af->nnode=strtol(*(data+2),NULL,10);
  }
  af->nodes=(struct onode *)malloc(af->nnode*sizeof(struct onode));
  af->ninit=(struct onode *)malloc(af->nnode*sizeof(struct onode));

  af->confs=(struct oconf *)malloc(6*af->nnode*sizeof(struct oconf));

  for(;n>0;n--) free(*(data+n-1));
  free(data);

  /*REFERENCE NODES*/
  data=fgetsbrk(ftext,&n);
  if(!strncmp(*(data+0),"Reference",9))
  {
    rn=strtol(*(data+2),NULL,10);
  }
  rnode=(struct onode *)malloc(rn*sizeof(struct onode));

  for(;n>0;n--) free(*(data+n-1));
  free(data);

  /*ELEMS*/
  data=fgetsbrk(ftext,&n);
  if(!strncmp(*(data+0),"Number",6))
  {
    af->nelem=strtol(*(data+3),NULL,10);
  }
  af->elems=(struct owire *)malloc(af->nelem*sizeof(struct owire));
  soffsets=(long int *)malloc(af->nelem*sizeof(long int));

  for(;n>0;n--) free(*(data+n-1));
  free(data);

  /*BEAMS*/
  data=fgetsbrk(ftext,&n);
  if(!strncmp(*(data+0),"Number",6))
  {
    nbeam=strtol(*(data+4),NULL,10);
  }

  for(;n>0;n--) free(*(data+n-1));
  free(data);

  /*PLATES*/
  data=fgetsbrk(ftext,&n);
  if(!strncmp(*(data+0),"Number",6))
  {
    nplate=strtol(*(data+4),NULL,10);
  }

  for(;n>0;n--) free(*(data+n-1));
  free(data);

sprintf(str,"NNODE=%d REFS=%d ELEMS=%d BEAMS=%d PLATES=%d",
        af->nnode,rn,af->nelem,nbeam,nplate);
MessageBox(NULL,str,"CADRE",MB_OK);

  /*SEARCH "NODAL COORDINATES"*/
  flag=0;
  while(flag==0)
  {
    fgets(str,256,ftext);
    if(!strncmp(str,"Nodal",5)) flag=1;
  }

  fgets(str,256,ftext);
  fgets(str,256,ftext);
  fgets(str,256,ftext);
  fgets(str,256,ftext);

  /*NODES*/
  for(i=0;i<af->nnode;i++)
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

    if(xmax<(af->nodes+i)->d[0]) xmax=(af->nodes+i)->d[0];
    if(xmin>(af->nodes+i)->d[0]) xmin=(af->nodes+i)->d[0];
    if(ymax<(af->nodes+i)->d[1]) ymax=(af->nodes+i)->d[1];
    if(ymin>(af->nodes+i)->d[1]) ymin=(af->nodes+i)->d[1];
    if(zmax<(af->nodes+i)->d[2]) zmax=(af->nodes+i)->d[2];
    if(zmin>(af->nodes+i)->d[2]) zmin=(af->nodes+i)->d[2];
  }

  fgets(str,256,ftext);

  /*REFERENCE NODES*/
  for(i=0;i<rn;i++)
  {
    (rnode+i)->loff=i;

    data=fgetsbrk(ftext,&n);

    (rnode+i)->code=strtol(*(data+0),NULL,10);
    (rnode+i)->d[0]=strtod(*(data+1),NULL);
    (rnode+i)->d[1]=strtod(*(data+2),NULL);
    (rnode+i)->d[2]=strtod(*(data+3),NULL);

    for(;n>0;n--) free(*(data+n-1));
    free(data);
  }

  /*SEARCH "ELEMENT DEFINITION DATA"*/
  flag=0;
  while(flag==0)
  {
    fgets(str,256,ftext);
    if(!strncmp(str,"Element definition",18)) flag=1;
  }

  fgets(str,256,ftext);
  fgets(str,256,ftext);
  fgets(str,256,ftext);
  fgets(str,256,ftext);
  fgets(str,256,ftext);

  /*BEAM ELEMENTS*/
  af->nsect=0;
  for(i=0;i<nbeam;i++)
  {
    (af->elems+i)->loff=i;
    (af->elems+i)->code=i+1;

    data=fgetsbrk(ftext,&n);

    hcode=strtol(*(data+2),NULL,10); /*HEAD NODE.*/
    tcode=strtol(*(data+3),NULL,10); /*TAIL NODE.*/
    rcode=strtol(*(data+4),NULL,10); /*REFERENCE NODE.*/

/*sprintf(str,"ELEM=%ld HEAD=%ld TAIL=%ld",
        (af->elems+i)->code,hcode,tcode);
MessageBox(NULL,str,"CADRE",MB_OK);*/

    s=(char *)malloc(256*sizeof(char));
    sprintf(s,"\0");
    for(j=5;j<n;j++)
    {
      strcat(s,*(data+j));
      if(j<n-1) strcat(s," ");
    }

/*sprintf(str,"ELEM=%ld SECT=%s",(af->elems+i)->code,s);
MessageBox(NULL,str,"CADRE",MB_OK);*/

    /*SECTION*/
    flag=0;
    for(j=0;j<af->nsect;j++)
    {
      if(!strcmp(*(library+j),s))
      {
        *(soffsets+i)=j;
        flag=1;
        free(s);
        break;
      }
    }

    if(flag==0)
    {
      (af->nsect)++;
      library=(char **)realloc(library,af->nsect*sizeof(char *));
      *(library+(af->nsect-1))=s;
      *(soffsets+i)=af->nsect-1;
    }

/*sprintf(str,"ELEM=%ld SECT=%s SOFFSET=%ld",
        (af->elems+i)->code,*(library+(af->nsect-1)),*(soffsets+i));
MessageBox(NULL,str,"CADRE",MB_OK);*/

    /*COORD ANGLE*/

    /*UNDER CONSTRUCTION.*/

    (af->elems+i)->cangle=0.0;

    for(ii=0;ii<=1;ii++)                 /*BOUNDARY.0:RIGID 1:HINGE*/
    {
      for(jj=0;jj<6;jj++)
      {
        (af->elems+i)->iconf[ii][jj]=0;
      }
    }

    for(;n>0;n--) free(*(data+n-1));
    free(data);

    offset=0;
    for(ii=0;ii<2;)                                        /*NODES.*/
    {
      if((af->nodes+offset)->code==hcode)
      {
        (af->elems+i)->node[0]=af->nodes+offset;
        ii++;
      }
      if((af->nodes+offset)->code==tcode)
      {
        (af->elems+i)->node[1]=af->nodes+offset;
        ii++;
      }
      offset++;
    }
  }

  /*PLATE ELEMENTS*/
  for(i=0;i<nplate;i++)
  {
    (af->elems+nbeam+i)->loff=nbeam+i;
    (af->elems+nbeam+i)->code=nbeam+i+1;

    data=fgetsbrk(ftext,&n);

    hcode=strtol(*(data+2),NULL,10); /*HEAD NODE.*/
    tcode=strtol(*(data+4),NULL,10); /*TAIL NODE.*/
    rcode=strtol(*(data+3),NULL,10); /*REFERENCE NODE.*/

    s=(char *)malloc(256*sizeof(char));
    sprintf(s,"\0");
    for(j=5;j<n;j++)
    {
      strcat(s,*(data+j));
      if(j<n-1) strcat(s," ");
    }

    /*SECTION*/
    flag=0;
    for(j=0;j<af->nsect;j++)
    {
      if(!strcmp(*(library+j),s))
      {
        *(soffsets+nbeam+i)=j;
        flag=1;
        free(s);
        break;
      }
    }

    if(flag==0)
    {
      (af->nsect)++;
      library=(char **)realloc(library,af->nsect*sizeof(char *));
      *(library+(af->nsect-1))=s;
      *(soffsets+nbeam+i)=af->nsect-1;
    }

    /*COORD ANGLE*/

    /*UNDER CONSTRUCTION.*/

    (af->elems+nbeam+i)->cangle=0.0;

    for(ii=0;ii<=1;ii++)                 /*BOUNDARY.0:RIGID 1:HINGE*/
    {
      for(jj=0;jj<6;jj++)
      {
        (af->elems+nbeam+i)->iconf[ii][jj]=0;
      }
    }

    for(;n>0;n--) free(*(data+n-1));
    free(data);

    offset=0;
    for(ii=0;ii<2;)                                        /*NODES.*/
    {
      if((af->nodes+offset)->code==hcode)
      {
        (af->elems+nbeam+i)->node[0]=af->nodes+offset;
        ii++;
      }
      if((af->nodes+offset)->code==tcode)
      {
        (af->elems+nbeam+i)->node[1]=af->nodes+offset;
        ii++;
      }
      offset++;
    }
  }

  af->sects=(struct osect *)malloc(af->nsect*sizeof(struct osect));

sprintf(str,"NSECT=%d",af->nsect);
MessageBox(NULL,str,"CADRE",MB_OK);

  /*SEARCH "ELEMENT AREAL AND ELASTIC PROPERTIES"*/
  flag=0;
  while(flag==0)
  {
    fgets(str,256,ftext);
    if(!strncmp(str,"Element areal",13)) flag=1;
  }

  fgets(str,256,ftext);
  fgets(str,256,ftext);
  fgets(str,256,ftext);
  fgets(str,256,ftext);
  fgets(str,256,ftext);
  fgets(str,256,ftext);
  fgets(str,256,ftext);
  fgets(str,256,ftext);
  fgets(str,256,ftext);

  /*SECTIONS*/
  for(i=0;i<nbeam;i++)
  {
    offset=*(soffsets+i);

    data=fgetsbrk(ftext,&n);

    (af->elems+i)->sect=af->sects+offset;

    (af->sects+offset)->loff=offset;
    (af->sects+offset)->code=offset+1;

    (af->sects+offset)->area=strtod(*(data+1),NULL);
    (af->sects+offset)->Ixx =strtod(*(data+2),NULL);
    (af->sects+offset)->Iyy =strtod(*(data+3),NULL);
    (af->sects+offset)->Jzz =strtod(*(data+4),NULL);

    for(;n>0;n--) free(*(data+n-1));
    free(data);

    data=fgetsbrk(ftext,&n);

    G=strtod(*(data+2),NULL);
    E=strtod(*(data+3),NULL);
    (af->sects+offset)->E=E;
    (af->sects+offset)->poi=0.5*E/G-1.0;

    (af->sects+offset)->type=COLUMN;

    for(;n>0;n--) free(*(data+n-1));
    free(data);
  }

  for(i=0;i<nplate;i++)
  {
    offset=*(soffsets+nbeam+i);

    data=fgetsbrk(ftext,&n);

    (af->elems+nbeam+i)->sect=af->sects+offset;

    (af->sects+offset)->loff=offset;
    (af->sects+offset)->code=offset+1;

    (af->sects+offset)->E  =strtod(*(data+1),NULL);
    (af->sects+offset)->poi=strtod(*(data+2),NULL);
    thick=strtod(*(data+3),NULL);

    (af->sects+offset)->area=1.0; /*UNDER CONSTRUCTION*/
    (af->sects+offset)->Ixx =0.0;
    (af->sects+offset)->Iyy =0.0;
    (af->sects+offset)->Jzz =0.0;

    (af->sects+offset)->type=WALL;

    for(;n>0;n--) free(*(data+n-1));
    free(data);
  }

  for(i=0;i<af->nsect;i++)
  {
    /*FLAGS,COLORS*/
    if((af->sects+i)->area==0.0 &&
       (af->sects+i)->Ixx ==0.0 &&
       (af->sects+i)->Iyy ==0.0 &&
       (af->sects+i)->Jzz ==0.0)
    {
      (af->sects+i)->dflag=0;
    }
    else (af->sects+i)->dflag=1;

    (af->sects+i)->name=*(library+i);

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

    (af->sects+i)->role=ROLENULL;
  }

  /*SEARCH "BOUNDARY CONDITIONS"*/
  flag=0;
  while(flag==0)
  {
    fgets(str,256,ftext);
    if(!strncmp(str,"Boundary",8)) flag=1;
  }

  fgets(str,256,ftext);
  fgets(str,256,ftext);
  fgets(str,256,ftext);
  fgets(str,256,ftext);
  fgets(str,256,ftext);

  /*BOUNDS*/
  af->nreact=0;
  while(1)
  {
    data=fgetsbrk(ftext,&n);

    if(!strcmp(*(data+0),"Loading"))
    {
      for(;n>0;n--) free(*(data+n-1));
      free(data);
      break;
    }

    ncode=strtol(*(data+0),NULL,10);

    for(i=0;i<af->nnode;i++)
    {
      if(ncode==(af->nodes+i)->code)
      {
        for(j=1;j<=6;j++)
        {
          offset=6*i+(j-1);

          (af->confs+offset)->value=0.0;

          if(!strcmp(*(data+j),"Fixed"))
          {
            (af->confs+offset)->iconf=1;
            af->nreact++;
          }
          if(!strcmp(*(data+j),"Free"))
          {
            (af->confs+offset)->iconf=0;
          }
        }
        break;
      }
    }

    for(;n>0;n--) free(*(data+n-1));
    free(data);

    fgets(str,256,ftext);
  }

  fgets(str,256,ftext);
  fgets(str,256,ftext);
  fgets(str,256,ftext);
  fgets(str,256,ftext);
  fgets(str,256,ftext);

  /*BOUNDS*/
  while(1)
  {
    data=fgetsbrk(ftext,&n);

    if(n<=0) break;

    ncode=strtol(*(data+0),NULL,10);

    for(i=0;i<af->nnode;i++)
    {
      if(ncode==(af->nodes+i)->code)
      {
        for(j=1;j<=6;j++)
        {
          offset=6*i+(j-1);

          if((af->confs+offset)->iconf==0)
          {
            (af->confs+offset)->value=strtod(*(data+j),NULL);
          }
        }
        break;
      }
    }

    for(;n>0;n--) free(*(data+n-1));
    free(data);

    fgets(str,256,ftext);
  }

  /*VIEW DATA*/
  sprintf(str,"%.1f",0.5*(xmax+xmin));
  SetDlgItemText((wmenu.childs+2)->hwnd,IDV_X,str);
  sprintf(str,"%.1f",0.5*(ymax+ymin));
  SetDlgItemText((wmenu.childs+2)->hwnd,IDV_Y,str);
  sprintf(str,"%.1f",0.5*(zmax+zmin));
  SetDlgItemText((wmenu.childs+2)->hwnd,IDV_Z,str);

  SetDlgItemText((wmenu.childs+2)->hwnd,IDV_R,"100000.0");
  SetDlgItemText((wmenu.childs+2)->hwnd,IDV_L,"3000.0");

  sprintf(str,"%.1f",xmax);
  SetDlgItemText((wmenu.childs+2)->hwnd,IDR_XMAX,str);
  sprintf(str,"%.1f",xmin);
  SetDlgItemText((wmenu.childs+2)->hwnd,IDR_XMIN,str);
  sprintf(str,"%.1f",ymax);
  SetDlgItemText((wmenu.childs+2)->hwnd,IDR_YMAX,str);
  sprintf(str,"%.1f",ymin);
  SetDlgItemText((wmenu.childs+2)->hwnd,IDR_YMIN,str);
  sprintf(str,"%.1f",zmax);
  SetDlgItemText((wmenu.childs+2)->hwnd,IDR_ZMAX,str);
  sprintf(str,"%.1f",zmin);
  SetDlgItemText((wmenu.childs+2)->hwnd,IDR_ZMIN,str);

  SetDlgItemText((wmenu.childs+2)->hwnd,IDV_GAXISLENGTH,"1000.0");
  SetDlgItemText((wmenu.childs+2)->hwnd,IDV_EAXISLENGTH,"500.0");
  SetDlgItemText((wmenu.childs+2)->hwnd,IDV_MFACTOR,"100.0");
  SetDlgItemText((wmenu.childs+2)->hwnd,IDV_GYOPITCH,"100.0");

  return;
}/*inputcadretomemory*/

void cadreoutputtomemory(FILE *ftext,FILE *fin,struct arclmframe *af)
/*TRANSLATE CADRE OUTPUTFILE TEXT INTO MEMORY.*/
{
  char **data,str[400],s[256];
  int i,j,k,n,flag,count;
  int nnode,nelem,rn,nbeam,nplate;
  long int ncode;
  double ddata;

  /*READ INITIAL DATA FROM INPUT FILE*/
  fseek(fin,0L,SEEK_SET);

  flag=0;
  while(flag==0)
  {
    fgets(str,256,fin);
    if(!strncmp(str,"Basic",5)) flag=1;
  }

  fgets(str,256,fin);
  fgets(str,256,fin);

  /*NODES*/
  data=fgetsbrk(fin,&n);
  if(!strncmp(*(data+0),"Structrual",10))
  {
    nnode=strtol(*(data+2),NULL,10);
  }
  for(;n>0;n--) free(*(data+n-1));
  free(data);

  if(nnode!=af->nnode) return;

  /*REFERENCE NODES*/
  data=fgetsbrk(fin,&n);
  if(!strncmp(*(data+0),"Reference",9))
  {
    rn=strtol(*(data+2),NULL,10);
  }
  for(;n>0;n--) free(*(data+n-1));
  free(data);

  /*ELEMS*/
  data=fgetsbrk(fin,&n);
  if(!strncmp(*(data+0),"Number",6))
  {
    nelem=strtol(*(data+3),NULL,10);
  }
  for(;n>0;n--) free(*(data+n-1));
  free(data);

  if(nelem!=af->nelem) return;

  /*BEAMS*/
  data=fgetsbrk(fin,&n);
  if(!strncmp(*(data+0),"Number",6))
  {
    nbeam=strtol(*(data+4),NULL,10);
  }
  for(;n>0;n--) free(*(data+n-1));
  free(data);

  /*PLATES*/
  data=fgetsbrk(fin,&n);
  if(!strncmp(*(data+0),"Number",6))
  {
    nplate=strtol(*(data+4),NULL,10);
  }
  for(;n>0;n--) free(*(data+n-1));
  free(data);

sprintf(str,"NNODE=%d ELEMS=%d BEAMS=%d PLATES=%d",
        nnode,nelem,nbeam,nplate);
MessageBox(NULL,str,"CADRE",MB_OK);

  /*RESULT FILE RESET*/
  fseek(ftext,0L,SEEK_SET);

  af->ddisp=(double *)malloc(6*(af->nnode)*sizeof(double));
                                       /*DISPLACEMENT:6 DIRECTIONS.*/
  af->melem=(struct memoryelem *)
			malloc((af->nelem)*sizeof(struct memoryelem));
                                        /*CODE,12 BOUNDS,12 STRESS.*/
  af->dreact=(double *)malloc((af->nreact)*sizeof(double));
                                                        /*REACTION.*/

  /*DISPLACEMENTS*/
  while(strncmp(str,"Nodal",5)) fgets(str,256,ftext);
  for(i=1;i<=6;i++) fgets(str,256,ftext);

  for(i=0;i<(af->nnode);i++)
  {
    data=fgetsbrk(ftext,&n);

    ncode=strtol(*(data+0),NULL,10);
    if(ncode!=(af->nodes+i)->code) return;

    for(j=0;j<=2;j++)
    {
      ddata=strtod(*(data+j+1),NULL);
	  ddata+=(af->ninit+i)->d[j];
      (af->nodes+i)->d[j]=ddata;
      *(af->ddisp+6*i+j)=ddata;
    }
    for(j=3;j<=5;j++)
    {
      ddata=strtod(*(data+j+1),NULL);
	  *(af->ddisp+6*i+j)=ddata;
    }

    for(;n>0;n--) free(*(data+n-1));
    free(data);

    fgets(str,256,ftext);
  }

  /*STRESS OF BEAMS*/
  while(strncmp(str,"Element Nodal",13)) fgets(str,256,ftext);
  for(i=1;i<=6;i++) fgets(str,256,ftext);

  for(i=0;i<nbeam;i++)
  {
    for(j=1;j<=4;j++) fgets(str,256,ftext);

    /*LOCAL HEAD*/
    data=fgetsbrk(ftext,&n);

    k=2;
    for(j=0;j<6;j++)
    {
      (af->elems+i)->stress[0][j]=strtod(*(data+k),NULL);
      k++;
    }

    for(;n>0;n--) free(*(data+n-1));
    free(data);

    fgets(str,256,ftext);

    data=fgetsbrk(ftext,&n);

    /*LOCAL TAIL*/
    k=1;
    for(j=0;j<6;j++)
    {
      (af->elems+i)->stress[1][j]=strtod(*(data+k),NULL);
      k++;
    }

/*sprintf(str,"ELEM=%d\n",(af->elems+i)->code);
sprintf(s,"HEAD %.3E %.3E %.3E %.3E %.3E %.3E\n",
        (af->elems+i)->stress[0][0],
        (af->elems+i)->stress[0][1],
        (af->elems+i)->stress[0][2],
        (af->elems+i)->stress[0][3],
        (af->elems+i)->stress[0][4],
        (af->elems+i)->stress[0][5]);
strcat(str,s);
sprintf(s,"TAIL %.3E %.3E %.3E %.3E %.3E %.3E",
        (af->elems+i)->stress[1][0],
        (af->elems+i)->stress[1][1],
        (af->elems+i)->stress[1][2],
        (af->elems+i)->stress[1][3],
        (af->elems+i)->stress[1][4],
        (af->elems+i)->stress[1][5]);
strcat(str,s);
MessageBox(NULL,str,"CADRE",MB_OK);*/

    for(;n>0;n--) free(*(data+n-1));
    free(data);

    (af->melem+i)->code=(af->elems+i)->code;
    for(j=0;j<6;j++)
    {
      (af->melem+i)->bond[0][j]=(af->elems+i)->iconf[0][j];
	  (af->melem+i)->bond[1][j]=(af->elems+i)->iconf[1][j];
      (af->melem+i)->stress[0][j]=(af->elems+i)->stress[0][j];
	  (af->melem+i)->stress[1][j]=(af->elems+i)->stress[1][j];
    }

    for(j=1;j<=4;j++) fgets(str,256,ftext);
  }

  /*STRESS OF PLATES*/
  while(strncmp(str,"Plate Stresses",13)) fgets(str,256,ftext);
  for(i=1;i<=6;i++) fgets(str,256,ftext);

  for(i=0;i<nplate;i++)
  {
    /*UNDER CONSTRUCTION*/
  }

  /*REACTION*/
  while(strncmp(str,"Reaction",8)) fgets(str,256,ftext);
  for(i=1;i<=6;i++) fgets(str,256,ftext);

  count=0;
  while(1)
  {
    data=fgetsbrk(ftext,&n);
    if(n<=0) return;

    ncode=strtol(*(data+0),NULL,10);

    for(j=0;j<(af->nnode);j++)
    {
      if(ncode==(af->nodes+j)->code) break;
    }

    for(k=0;k<6;k++)
    {
      if((af->confs+6*j+k)->iconf==1)
      {
        ddata=strtod(*(data+k+1),NULL);
        *(af->dreact+count)=ddata;
        count++;
      }
    }

    for(;n>0;n--) free(*(data+n-1));
    free(data);

    fgets(str,256,ftext);
  }

  return;
}/*cadreoutputtomemory*/
#if 0
void obspoints(int obscodes[32])   /*ujioka and isebo for observation points*/
{
  FILE *fobs;
  char dir[]=DIRECTORY;						      	        /*DATA DIRECTORY*/
  int i,ii,iii,j;
  int nobserve;
//  int obscodes[32];
  int nstr,pstr,ident,ncv;
  char **data,str[256];
  int readflag = 1;
  fobs=fgetstofopenII(dir,"r","observationpoints.obs");         //for obs file

  if(fobs==NULL)
  {
    printf("couldn't open the obs-file\n")  ;
    getchar();
	exit(EXIT_FAILURE);
  }

  while (readflag == 1)
  {
	data=fgetsbrk(fobs,&nstr);

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

		if(nstr-pstr>=2 && !strcmp(*(data+pstr),"NOBSERVE"))
		{
		  pstr++;
		  nobserve=(int)strtol(*(data+pstr),NULL,10);
          ident++;
        }

		if(nstr-pstr>=2 && !strcmp(*(data+pstr),"CODES"))
        {
          for (iii = 0; iii < nobserve; iii++)            //ujok
          {
            pstr++;
            obscodes[iii] = (int)strtol(*(data+pstr),NULL,10);
          }
          ident++;
          readflag = 0;
        }
      }
    }
  }
  fclose(fobs);
/*
  sprintf(str,"nobserve=%d",nobserve);
  sprintf(str,"observecode=%d",obscodes[0]);
  MessageBox(NULL,str,"OBSERVATION",MB_OK);
*/
  errormessage("OBSERVATION POITS CODE: ");

  for (iii = 0; iii < nobserve; iii++)
  {
     sprintf(str," %d",obscodes[iii]);
     errormessage(str);
  }

  errormessage("\n");

  return;
}
#endif

/*drawelemcmqcheck*/  /*FOR CMQ CHECK*/ /*kaza & uji for Lunar/MarsBase*/
void drawelemcmqcheck(HDC hdc,
				  struct viewparam vp,
				  struct onode n1,struct onode n2,double cangle,double pi,double pj)
{
  HPEN hpen,ppen;                                   /*HANDLE OF PEN*/
  int i;
  double elength;
  double **drccos,**tdrccos;
  double arrowlength=vp.dparam.eaxis*pi;
  double arrowlengthII=vp.dparam.eaxis*pj;
  struct onode on,hn,tn; /*ORIGIN,HEAD,TAIL.*/
  char str[256];

  elength=distancedotdot(n1,n2);

  drccos=directioncosine(n1.d[GX],n1.d[GY],n1.d[GZ],
						 n2.d[GX],n2.d[GY],n2.d[GZ],
						 cangle);                        /*[DRCCOS]*/
  tdrccos=matrixtranspose(drccos,3);                          /*[T]*/

  SetTextColor(hdc,RGB(255,0,255));                 /*AXIS MAGENTA.*/
  hpen=CreatePen(PS_SOLID,1,RGB(255,0,255));       /* LINE MAGENTA.*/
  ppen=(HPEN)SelectObject(hdc,hpen);

  on=n1;                                            /*LOCAL ORIGIN.*/

  hn=setcoord((0.5*elength-0.5*arrowlength),0.0,0.0);
  tn=setcoord((0.5*elength-0.5*arrowlength),0.0,arrowlength);
  sprintf(str,"%.2f",pi);
  drawtextonlocaldot(hdc,vp,tdrccos,on,tn,str);      /*TEXT Pi=-Qyi*/
  drawlocalarrow(hdc,vp,tdrccos,on,hn,tn);          /*ARROW Pi=-Qyi*/

  hn=setcoord((0.5*elength+0.5*arrowlength),0.0,0.0);
  tn=setcoord((0.5*elength+0.5*arrowlength),0.0,arrowlengthII);
  sprintf(str,"%.2f",pj);
  drawtextonlocaldot(hdc,vp,tdrccos,on,tn,str);     /*TEXT Pj=-Qyj*/
  drawlocalarrow(hdc,vp,tdrccos,on,hn,tn);         /*ARROW Pj=-Qyj*/

  for(i=0;i<=2;i++) free(*(drccos+i));
  free(drccos);
  for(i=0;i<=2;i++) free(*(tdrccos+i));
  free(tdrccos);

  SelectObject(hdc,ppen);
  DeleteObject(hpen);

  return;
}/*drawelemcmqcheck*/  /*FOR CMQ CHECK*/ /*kaza & uji for Lunar/MarsBase*/

int gcomponentadd(struct gcomponent *mtx1, // 150219 fukushima for gcomponentadd
             struct gcomponent *mtx2,
             int msize)
/*ASSEMBLAGE ELEMENT MATRIX INTO GLOBAL MATRIX.*/
{
  long int i,j;
  double data,newvalue;
  struct gcomponent *gcomp,*gdown1,*gdown2,*pdown1;

  for(j=1;j<=msize;j++)
  {
	gdown1=(mtx1+(j-1));
	gdown2=(mtx2+(j-1));
    while(gdown2!=NULL) /*DOWNWARD.*/
    {
      i=gdown2->m;
      data=gdown2->value;
	  while((gdown1->m)<i && (gdown1->down)!=NULL)
      {
        pdown1=gdown1;
        gdown1=gdown1->down;
      }

      if(gdown1->m==i)                             /*COMPONENT EXISTENT.*/
      {
        newvalue=gdown1->value+data;
        if(newvalue==0.0 && i!=j)     /*COMPONENT VANISHED.EXCEPT DIAGONAL.*/
        {
          pdown1->down=gdown1->down;
          free(gdown1);

          comps--;
        }
        else
        {
          gdown1->value=newvalue;
        }
      }
      else                                           /*COMPONENT EMPTY.*/
      {
		gcomp=(struct gcomponent *)malloc(sizeof(struct gcomponent));
        if(gcomp==NULL)
        {
          errormessage("GCOMPONENTADD:MEMORY INSUFFICIENT.");
          return 0;
        }

        gcomp->m=(unsigned short int)i;
        gcomp->value=data;

        if((gdown1->m)<i)                                  /*ADD TO ROW.*/
        {
          gcomp->down=NULL;
          gdown1->down=gcomp;
        }
        else if((gdown1->m)>i)                          /*FILL INTO ROW.*/
        {
          gcomp->down=pdown1->down;
          pdown1->down=gcomp;
        }
        comps++;
      }

      gdown2=gdown2->down;
    }
  }
  return 0;
}

