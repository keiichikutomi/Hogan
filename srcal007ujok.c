
/*SRCAL004.H:SRC SUBROUTINES SINCE 1995.12.21.JUNSATO.*/
/*HEADER FOR HOGAN001.C.*/
/*BASED ON SRCAM005.H.*/
/*LAST CHANGE:2011.4.14.*/



/*DRAW S,RC,SRC,PC SECTIONS.*/
/*SECTION FROM SECTION LIST FILE.*/
/*MDIFIED TO RC STANDARD 2010.*/

#define STYPE_S     1                                      /*SECTION TYPE*/
#define STYPE_RC    2
#define STYPE_SRC   3
#define STYPE_PC    4
#define STYPE_WOOD  5
#define STYPE_GLASS 6
#define STYPE_ACRYL 7
#define STYPE_ALUMI 8

#define PLONG     0 /*PERIOD*/
#define PSHORT    1
#define PULTIMATE 2

#define SX 0 /*AXIS OF SECTION*/
#define SY 1

#define HSTRONG 0 /*AXIS OF H*/
#define HWEAK   1

#define STEEL_NULL   0 /*STEEL TYPE && ALUMI TYPE*/
#define STEEL_RECTS  1
#define STEEL_HKYOU  2
#define STEEL_HWEAK  3
#define STEEL_RPIPE  4
#define STEEL_CPIPE  5
#define STEEL_PLATE  6
#define STEEL_ANGLE  7
#define STEEL_TKYOU  8
#define STEEL_TWEAK  9
#define DABO        10
#define STEEL_CKYOU 11 /*CHANNEL*/
#define STEEL_CWEAK 12


#define HEAD 0 /*END OF ELEMENT*/
#define TAIL 1

#define PEND 0 /*POSITION OF ELEMENT FOR PC*/
#define PMID 1

#define UPPER 0 /*TENSED SIDE FOR BENDING OF PC GIRDER*/
#define LOWER 1

#define STRNDA 0 /*PC BAR TYPE:A,B,C,WIRE STRAND*/
#define STRNDB 1
#define STRNDC 2
#define STRNDW 3

#define MAXSECT 1000 /*MAX OF SECTIONS*/
#define MAXCMQ   100 /*MAX OF CMQ LIST*/
#define MAXCOMMENT 100 /*MAX OF COMMENTS PER ONE SECTION*/

#define MAXSRECT    5                     /*MAX OF STEEL RECTANGLES*/
#define MAXREINS  200                       /*MAX OF REINFORCEMENTS*/
#define MAXCRECT    3                  /*MAX OF CONCRETE RECTANGLES*/
#define MAXSTRND   50                           /*MAX OF PC STRANDS*/

#define SRCANBCLNGCONDENSATION 0       /*UJIOKA:BUCKLING CONDENSATION RESULT*/
                                       /*ONLY FOR S PLATE*/

struct materials{
					 double sE,sF,sft,sfc,sfb,sfs,sftu,sfcu;       		 /*S*/
					 double rE,rft,rfc,wft,rftu,rfcu,wfp;          		 /*R*/
					 double cE,Fc,cfc,cfs,srcfc,srcfs,cfcu;         	 /*C*/
					 double pcE,pcF,pcfc,pcft,pcfs,pcfcu,pcfsu;    		/*PC*/
					 double pcfco,pcfto;                 	   /*PC FOR DEAD*/
					 double stfact;    		/*EFFECTIVE FACTOR FOR PC STRAND*/
					 double stE,stft[4],stftu[4];              		/*STRAND*/
					 double gE,gF,gft,gfc,gfb,gfs,gftu,gfcu;   		 /*GLASS*/
					 double aE,aF,aft,afc,afb,afs,aftu,afcu;    	 /*ACRYL*/
					 double alE,alF,alft,alfc,alfb,alfs,alftu,alfcu; /*ALUMI*/
				};

struct materialrect{double top,bottom,left,right;
					double F;};                   /*RECTANGLE[cm],F*/
struct reinforcement{double area,x,y,F;};     /*REIN AREA,X,Y[cm],F*/
struct reinshear{int n;                            /*NUMBER OF REIN*/
                 double area,pitch,F;};     /*REIN AREA,PITCH[cm],F*/
struct pcstrand{
                int type;                  /*STRAND TYPE:A,B,C,WIRE*/
                double area[2],x[2],y[2];       /*REIN AREA,X,Y[mm]*/
                double Ni[2];                 /*INITIAL TENSION[tf]*/
               };
struct steelform{
                 int type;
                 double H,B,tw,tf;
                 double F;
                };                                     /*STEEL FORM*/
struct woodform{
                int type;
                double H,B,tw,tf;
                double Fc,Ft,Fb,Fs;
               };                                       /*WOOD FORM*/

struct section{
               long int code;                     /*CODE OF SECTION*/
               int soff;                        /*OFFSET OF SECTION*/
			   long int ocode;           /*ORIGINAL CODE OF SECTION*/ // 150429 fukushima for original section
               int stype,etype;         /*SECTION TYPE,ELEMENT TYPE*/
               int nsteel,nrein,nconc;     /*STEELS,REINS,CONCRETES*/
               int nstrnd;                             /*PC STRANDS*/
               struct materialrect srect[MAXSRECT];   /*STEEL RECTS*/
               struct reinforcement rein[MAXREINS];         /*REINS*/
               struct materialrect crect[MAXCRECT]; /*CONCRETE RECT*/
               struct pcstrand strnd[MAXSTRND];        /*PC STRANDS*/
               double shearrein[2],wF;     /*RATE OF REIN FOR SHEAR*/
               double thick,wlength,wheight,windowrate,cF,fsl;
               double face[2][2];          /*FACE LENGTH[AXIS][END]*/
               double safetya;       /*MAX RATE OF SAFETY ALLOWABLE*/
               double safetyu;        /*MAX RATE OF SAFETY ULTIMATE*/

               int ncomment;
               char comment[MAXCOMMENT][80];
               struct reinshear srein[2];      /*SHEAR REIN BY AREA*/
               struct steelform sform;        /*INPUT BY STEEL FORM*/
               struct woodform wform;          /*INPUT BY WOOD FORM*/

               double bblength[2],bbfact[2];     /*BUCKLING BENDING*/
               double btlength[2],btfact[2];     /*BUCKLING TORTION*/
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
int getcodelist(FILE *flist,long int codelist[]);
int getcmqlist(FILE *fcmq,int cmqcode[],double Mo[]);

double elementlength(struct element elem);
double wirelength(struct owire elem);
double walllength(struct element elem);
double wallheight(struct element elem);
double slabwidth(struct element elem,struct element pair);
double slablength(struct element elem,struct element pair);

void translatesection(struct section sect,
                      struct materials m,
                      double *As,double *Ac,double *Ar,
                      double *Yg,
                      int axis);
double allowablebendingofrc(struct element elem,
                            struct materials m,
                            int axis,
                            double Nd,int period);
double allowablebendingofsrc(struct element elem,
                             struct materials m,
                             int axis,
                             double Nd,int period);
int ultimateaxialforceofsrc(struct element elem,
                            double *Nmax,double *Nmin,
                            double *sN,double *sM,
                            double *rN,double *rM,
                            double *cN,double *cM);
/*double ultimatebendingofsteel(struct section sect,
                              double E,double F,
                              double lk,
                              int axis,
                              double Nd);*/
double ultimatebendingofsteel(struct section sect,
                              double E,double F,
                              double lkxx,double lkyy,
                              int axis,
                              double Nd);
double allowablebendingofsteel(int period,
                               struct section sect,
                               double E,
                               double lk,
                               int axis,
                               double Nd);
double allowablebendingofalumi(int period,
                               struct section sect,
                               double E,
                               double lk,
                               int axis,
                               double Nd);
double allowablestressofflatbar(int period,
                                struct section sect,
                                double E,
                                double L,
                                struct stress *st,
                                int iend,
                                double *Ncr,
                                double *Qcrx,double *Qcry,
                                double *Mcrx,double *Mcry);
double allowablestressofflatbar_bc(int period,    /*UJIOKA:BUCKLING CONDENSATION RESULT*/
                                   struct section sect,
                                   double cNcr, //INDIVIDUAL ELEMENT'S BUCKLING LOAD
                                   double E,
                                   double L,
                                   struct stress *st,
                                   int iend,
                                   double *Ncr,
                                   double *Qcrx,double *Qcry,
                                   double *Mcrx,double *Mcry);
double allowablestressofflatbaralumi(int period,
                                     struct section sect,
                                     double E,
                                     double L,
                                     struct stress *st,
                                     int iend,
                                     double *Ncr,
                                     double *Qcrx,double *Qcry,
                                     double *Mcrx,double *Mcry);
/*double allowablebendingofwood(int period,
                              struct section sect,
                              double lk,
                              int axis,
                              double Nd);*/
double allowablebendingofwood(int period,
                              struct section sect,
                              double lkxx,double lkyy,
                              int axis,
                              double Nd);
double allowablestressofdabo(int period,
                             struct section sect,
                             struct stress *st,
                             int iend,
                             double *Na,double *Qa,double *Ma);
double ultimatebendingofsrc(struct element elem,
                            int axis,
                            double Nd,
                            double *sN,double *sM,
                            double *rN,double *rM,
                            double *cN,double *cM);

int addnewbound(double bounds[],int *nbound,double newbound);
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
double boundaryvalueultimate(struct element elem,/*struct materials m,*/
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
                             double *gN,double *gM,
                             double reinrate);

double allowableshearofrclong(struct element elem,
                              /*struct materials m,*/
                              int axis,
                              double Qd,
                              double Md);
double allowableshearofrcshort(struct element elem,
                               /*struct materials m,*/
                               int axis,
                               double Qd,
                               double Md);
double ultimateshearofrc(struct element elem,
                         struct materials m,
                         int axis,
                         double Nd,double Qd,double Md);
double ultimatetorsionofrc(struct section sect);

double allowableshearofsrclong(struct element elem,
                               /*struct materials m,*/
                               int axis,
                               int haxis,
                               double srcQd,
                               double srcMd);
double allowableshearofsrcshort(struct element elem,
                                /*struct materials m,*/
                                int axis,
                                int haxis,
                                double srcQd,
                                double srcMd);
double ultimateshearofsrc(struct element elem,
                          /*struct materials m,*/
                          int axis,
                          int haxis,
                          double rcQd,
                          double rcMd);

double allowultimshearofrcwall(struct element elem,
                               /*struct materials m,*/
                               double l, /*LENGTH[cm]*/
                               int period);
double allowableshearofwoodwall(struct element elem,
                                double l, /*LENGTH[cm]*/
                                int period);

double allowabletensionofs(double sF,
                           struct section sect);
double allowabletensionofsteel(double sF,struct section sect);
double allowabletensionofalumi(double sF,struct section sect);
double allowablecompressionofs(double E,double F,
                               double lk,
                               struct section sect,int axis);
/*double allowablecompressionofsteel(double E,double F,double lk,
                                   struct section sect);*/
double allowablecompressionofsteel(double E,double F,
                                   double lkxx,double lkyy,
                                   struct section sect);
double allowablecompressionofalumi(double E,double F,double lkxx,double lkyy,
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
double allowultimshearofsteel(int period,double F,
                              struct section sect,int axis);
double allowultimshearofalumi(int period,double F,
                              struct section sect,int axis);
double allowableshearofwood(int period,struct section sect);

int steelcoefficients(struct section sect,
                      double *A,
                      double *Ixx,double *Iyy,
                      double *Zxx,double *Zyy,
                      double *ixx,double *iyy);
int alumicoefficients(struct section sect,
                      double *A,
                      double *Ixx,double *Iyy,
                      double *Jz ,double *Cw,
                      double *Zxx,double *Zyy,
                      double *ixx,double *iyy);
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
char prj[20];
double jis;                                        /*1.1=JIS STEEL.*/
struct materials gmaterial;
extern int globalmessageflag;
extern int globalstatus;

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
/*"CODE"  CODENUMBER MATERIALTYPE ELEMENTTYPE*/
/*"SRECT" BOTTOMLEFTX BOTTOMLEFTY TOPRIGHTX TOPRIGHTY*/
/*"REINS" AREA X Y*/
/*"CRECT" BOTTOMLEFTX BOTTOMLEFTY TOPRIGHTX TOPRIGHTY*/
/*"HOOPS" HOOPRATEX HOOPRATEY*/
/*"THICK" THICKNESS*/
/*"SREIN" SHEARREINRATE*/
/*"WRECT" WINDOWWIDTH WINDOWHEIGHT*/
/*"SAREA" AREA OF STEEL BRACE*/

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

#define PROJECT "got"
/*#define PROJECT "zoo"*/
/*#define PROJECT "yachi"*/
#define HENSHINX 0.30
#define HENSHINY 0.30
#define CMQFILE "c:\\cdocs\\srcan\\data\\srcan.cmq"

int srcan001(char fname[])
{
  int ii,jj,ie,is,ic,ia,in,na,n1,n2,n3,ielem,isect,ni,nj;
  int hoko;                 /*DIRECTION OF HORIZONTAL LOAD 0:X 1:Y.*/
  int sign;        /*SIGN OF HORIZONTAL LOAD 0:POSITIVE 1:NEGATIVE.*/
  int calc[3];              /*CALCULATION FLAG LONG,SHORT,ULTIMATE.*/
  double dsign;
  double lN,lQ,eN; /*LONG,EARTHQUAKE STRESS.*/
  double sN,sQ[2][2],sMt,sM[2][2];  /*s:SHORT [END][AXIS]*/
  double uN,uQ[2][2],uQe[2][2],uMt,uM[2][2],Qp[2];  /*u:ULTIMATE*/
  double Ma,Mu,Muhead,Mutail,Qa,Qu; /*a:ALLOABLE u:ULTIMATE*/
  double Ns,Ms,Nr,Mr,Nc,Mc;            /*s:STEEL r:REIN c:CONCRETE.*/
  double Ncr,Qcrx,Qcry,Mcrx,Mcry;
  double nfact=1.0,mfact=1.0,qfact=2.0,Co=0.2,Ds=0.3,Fes;
  double wafact=2.0,wufact=3.0; /*FACTORS FOR WALL.*/
  double h,l; /*WALL HEIGHT,LENGTH*/
  double h0[2],l0; /*INNER HEIGHT,LENGTH*/
  double facei,facej,face[2][2];
  struct section sect1;
  long int codelist[MAXSECT];
  int cmqcode[MAXCMQ];
  double Mo[MAXCMQ],Mm;
  double rate,rate1,rate2,ratema,ratemu,rateqa,ratequ;
  double marate[MAXSECT],murate[MAXSECT];
  double qarate[MAXSECT],qurate[MAXSECT];
  double wrate1,wrate2; /*RATE OF WINDOW.*/

  double *cNcr;  /*Buckling Condensation Result*/
  double icNcr;
  int bclngcondensationflag;
  /*int check;*/

  FILE *fin,*fcmq=NULL,*flist,*fsafe,*frate;
  FILE *fz=NULL,*fx=NULL,*fy=NULL;
  char txt[400],non[80],str[256];
  char strhoko[2][10]={"X","Y"};
  /*char straxis[2][10]={"x","y"};*/
  /*char strend[2][10]={"i","j"};*/
  /*char strsign[2][10]={"+","-"};*/
  char strhugo[2][10]={"ê≥","ïâ"};
  char strQ[5][2][256],strM[3][2][256]; /*STRINGS FOR OUTPUT.*/
  int nnode,nelem,nsect,ncmq=0,soffset;
  double As,Ar,Ac,Yg,E,G;
  double si=9.80665;                                    /*SI UNITS.*/
  struct element elem,pair;
  struct stress *pstress[2],estress; /*HEAD,TAIL*/

  struct structnode *nodes;
  struct element *elems;
  struct section *sects;

  int tensedside;
  int slabcount=0;
  int dabocount=0;

  char dir[]=DIRECTORY;

  strcpy(prj,PROJECT);

  GetDlgItemText((wmenu.childs+2)->hwnd,ID_SECTIONFILE,non,80);
  if(!strncmp(non,"zooii",5)) strcpy(prj,"zooii");
  else if(!strncmp(non,"zoo",   3)) strcpy(prj,"zoo");
  else if(!strncmp(non,"yachi", 5)) strcpy(prj,"yachi");
  else if(!strncmp(non,"huru",  4)) strcpy(prj,"huru");
  else if(!strncmp(non,"waka",  4)) strcpy(prj,"waka");
  else if(!strncmp(non,"aii",   3)) strcpy(prj,"aii");
  else if(!strncmp(non,"noa",   3)) strcpy(prj,"noa");
  else if(!strncmp(non,"hira",  4)) strcpy(prj,"hira");
  else if(!strncmp(non,"athens",6)) strcpy(prj,"athens");
  else if(!strncmp(non,"gyo",   3)) strcpy(prj,"gyo");
  else if(!strncmp(non,"hiro",  4)) strcpy(prj,"hiro");
  else if(!strncmp(non,"hakoii",6)) strcpy(prj,"hakoii");
  else if(!strncmp(non,"naka",  4)) strcpy(prj,"naka");
  else if(!strncmp(non,"tos",   3)) strcpy(prj,"tos");
  else if(!strncmp(non,"izu",   3)) strcpy(prj,"izu");
  else if(!strncmp(non,"kiku",  4)) strcpy(prj,"kiku");
  else if(!strncmp(non,"nade",  4)) strcpy(prj,"nade");
  else if(!strncmp(non,"tohu",  4)) strcpy(prj,"tohu");
  else if(!strncmp(non,"yagi",  4)) strcpy(prj,"yagi");
  else if(!strncmp(non,"karu",  4)) strcpy(prj,"karu");
  else if(!strncmp(non,"ashi",  4)) strcpy(prj,"ashi");
  else if(!strncmp(non,"iwami", 5)) strcpy(prj,"iwami");
  else if(!strncmp(non,"haya",  4)) strcpy(prj,"haya");
  else if(!strncmp(non,"banga", 5)) strcpy(prj,"banga");
  else if(!strncmp(non,"kiyo",  4)) strcpy(prj,"kiyo");
  else if(!strncmp(non,"tep",   3)) strcpy(prj,"tep");
  else if(!strncmp(non,"mogi",  4)) strcpy(prj,"mogi");
  else if(!strncmp(non,"koshi", 5)) strcpy(prj,"koshi");
  else if(!strncmp(non,"vene",  4)) strcpy(prj,"vene");
  else if(!strncmp(non,"kamime",6)) strcpy(prj,"kamime");
  else if(!strncmp(non,"toko",  4)) strcpy(prj,"toko");
  else if(!strncmp(non,"maebasi",7)) strcpy(prj,"maebasi");
  else if(!strncmp(non,"machi", 5)) strcpy(prj,"machi");
  else if(!strncmp(non,"ginhu", 5)) strcpy(prj,"ginhu");
  else if(!strncmp(non,"stainedmodel",12)) strcpy(prj,"stainedmodel");
  else if(!strncmp(non,"sichi", 5)) strcpy(prj,"sichi");
  else if(!strncmp(non,"ebis", 4)) strcpy(prj,"ebis");
  else if(!strncmp(non,"stick", 5)) strcpy(prj,"stick");
  else if(!strncmp(non,"sydney", 6)) strcpy(prj,"sydney");
  else if(!strncmp(non,"ghouse",6))    strcpy(prj,"ghouse");
  else if(!strncmp(non,"LunarBase",9))    strcpy(prj,"lunarbase");
  else if(!strncmp(non,"plywood",7))    strcpy(prj,"plywood");
  else if(!strncmp(non,"kirigami",8))   strcpy(prj,"kirigami");

  sprintf(txt,"PROJECT:%s\n",prj);

  sprintf(non,"FILES  :%s.TST,RLT,RAT\n",fname);
  strcat(txt,non);

  sprintf(non,"Continue or Cancel"); strcat(txt,non);

  if(globalmessageflag==1 &&
     MessageBox(NULL,txt,"SRCan",MB_OKCANCEL)==IDCANCEL)
  {
    return 0;
  }

  /*CALCULATION FLAG*/
  calc[PLONG]    =1;
  calc[PSHORT]   =1;
  calc[PULTIMATE]=0;
  if(globalmessageflag==1 &&
     MessageBox(NULL,"Calculate Ultimate.","SRCan",MB_YESNO)==IDYES)
  {
    calc[PULTIMATE]=1;
  }
  else
  {
    calc[PULTIMATE]=0;
  }

  bclngcondensationflag=SRCANBCLNGCONDENSATION;


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
  fprintf(fout0,"ífñ éZíË \"S,RC,SRC\" í∑ä˙,íZä˙,èIã«\n");
  fprintf(fout0,"égópÉtÉ@ÉCÉã\n");
  fprintf(fout0,"ì¸óÕÉfÅ[É^               =%s\n",non);
  fprintf(fsafe,"%s\n",non);

  GetDlgItemText((wmenu.childs+2)->hwnd,ID_OUTPUTFILEZ,non,80);
  fprintf(fout0,"âîíºâ◊èdéûâêÕåãâ        =%s\n",non);
  fprintf(fsafe,"%s\n",non);
  GetDlgItemText((wmenu.childs+2)->hwnd,ID_OUTPUTFILEX,non,80);
  fprintf(fout0,"êÖïΩâ◊èdéûâêÕåãâ  Xï˚å¸ =%s\n",non);
  fprintf(fsafe,"%s\n",non);
  GetDlgItemText((wmenu.childs+2)->hwnd,ID_OUTPUTFILEY,non,80);
  fprintf(fout0,"                   Yï˚å¸ =%s\n",non);
  fprintf(fsafe,"%s\n",non);

  GetDlgItemText((wmenu.childs+2)->hwnd,ID_SECTIONFILE,non,80);
  fprintf(fout0,"âºíËífñ                  =%s\n",non);
  fprintf(fsafe,"%s\n",non);
  fprintf(fout0,"\níPà ån tf(kN),tfm(kNm)\n");

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

  /*fcmq=fopen(CMQFILE,"r"); if(fcmq==NULL)  return 0;*/

  /*INITIAL SETTING*/
  jis=1.0;                                          /*1.1=JIS STEEL.*/

  if(!strcmp(prj,"ken")  || !strcmp(prj,"tei") ||
     !strcmp(prj,"cho")  || !strcmp(prj,"aki") ||
     !strcmp(prj,"huna") || !strcmp(prj,"aka") ||
     !strcmp(prj,"kyok") || !strcmp(prj,"mae") ||
     !strcmp(prj,"del")  || !strcmp(prj,"got") ||
     !strcmp(prj,"zoo")  || !strcmp(prj,"yachi") ||
     !strcmp(prj,"huru") || !strcmp(prj,"waka") ||
     !strcmp(prj,"aii")  || !strcmp(prj,"noa") ||
     !strcmp(prj,"hira") || !strcmp(prj,"athens") ||
     !strcmp(prj,"zooii")|| !strcmp(prj,"gyo") ||
     !strcmp(prj,"hiro") || !strcmp(prj,"hakoii") ||
     !strcmp(prj,"naka") || !strcmp(prj,"tos") ||
     !strcmp(prj,"izu")  || !strcmp(prj,"kiku")  ||
     !strcmp(prj,"nade") || !strcmp(prj,"tohu") ||
     !strcmp(prj,"yagi") || !strcmp(prj,"karu") ||
     !strcmp(prj,"ashi") || !strcmp(prj,"iwami") ||
     !strcmp(prj,"haya") || !strcmp(prj,"banga") ||
     !strcmp(prj,"kiyo") || !strcmp(prj,"tep") ||
     !strcmp(prj,"mogi") || !strcmp(prj,"koshi") ||
	 !strcmp(prj,"vene") || !strcmp(prj,"kamime") ||
	 !strcmp(prj,"toko") || !strcmp(prj,"maebasi") ||
	 !strcmp(prj,"machi")|| !strcmp(prj,"ginhu") ||
	 !strcmp(prj,"stainedmodel")|| !strcmp(prj,"sichi") ||
	 !strcmp(prj,"ebis") || !strcmp(prj,"stick") ||
	 !strcmp(prj,"sydney")||
     !strcmp(prj,"ghouse")||
     !strcmp(prj,"lunarbase")||!strcmp(prj,"plywood")||
     !strcmp(prj,"kirigami"))

  {
    fgetinitial2(fin,&nnode,&nelem,&nsect);   /*INPUT ARCLM INITIAL.*/
  }
  else
  {
    fgetinitial(fin,&nnode,&nelem,&nsect,&E,&G);    /*INPUT INITIAL.*/
  }

  sects=(struct section *)malloc(nsect*sizeof(struct section));
  if(sects==NULL) return 0;
  nodes=(struct structnode *)malloc(nnode*sizeof(struct structnode));
  if(nodes==NULL) return 0;
  elems=(struct element *)malloc(nelem*sizeof(struct element));
  if(elems==NULL) return 0;

  if(!strcmp(prj,"ken")  || !strcmp(prj,"tei") ||
     !strcmp(prj,"cho")  || !strcmp(prj,"aki") ||
     !strcmp(prj,"huna") || !strcmp(prj,"aka") ||
     !strcmp(prj,"kyok") || !strcmp(prj,"mae") ||
     !strcmp(prj,"del")  || !strcmp(prj,"got") ||
     !strcmp(prj,"zoo")  || !strcmp(prj,"yachi") ||
     !strcmp(prj,"huru") || !strcmp(prj,"waka") ||
     !strcmp(prj,"aii")  || !strcmp(prj,"noa") ||
     !strcmp(prj,"hira") || !strcmp(prj,"athens") ||
     !strcmp(prj,"zooii")|| !strcmp(prj,"gyo") ||
     !strcmp(prj,"hiro") || !strcmp(prj,"hakoii") ||
     !strcmp(prj,"naka") || !strcmp(prj,"tos") ||
     !strcmp(prj,"izu")  || !strcmp(prj,"kiku")  ||
     !strcmp(prj,"nade") || !strcmp(prj,"tohu") ||
     !strcmp(prj,"yagi") || !strcmp(prj,"karu") ||
     !strcmp(prj,"ashi") || !strcmp(prj,"iwami") ||
     !strcmp(prj,"haya") || !strcmp(prj,"banga") ||
     !strcmp(prj,"kiyo") || !strcmp(prj,"tep") ||
     !strcmp(prj,"mogi") || !strcmp(prj,"koshi") ||
     !strcmp(prj,"vene") || !strcmp(prj,"kamime") ||
     !strcmp(prj,"toko") || !strcmp(prj,"maebasi") ||
     !strcmp(prj,"machi")|| !strcmp(prj,"ginhu") ||
	 !strcmp(prj,"stainedmodel")|| !strcmp(prj,"sichi") ||
	 !strcmp(prj,"ebis") || !strcmp(prj,"stick") ||
	 !strcmp(prj,"sydney")||
  	 !strcmp(prj,"ghouse")||
     !strcmp(prj,"lunarbase")||!strcmp(prj,"plywood")||
     !strcmp(prj,"kirigami"))

  {
    inputfiletomemory2(fin,&nnode,&nelem,&nsect,sects,nodes,elems);
  }
  else
  {
    inputfiletomemory(fin,&nnode,&nelem,&nsect,sects,nodes,elems);
  }

  comments(fout0);
  /*if(MessageBox(NULL,"Continue or Cancel","SRCan",
                MB_OKCANCEL)==IDCANCEL) return 0;*/

  for(is=0;is<nsect;is++) /*SECTIONS INTO MEMORY.*/
  {
	/* getsectionform(flist,(sects+is)->code,(sects+is)); */

	/* 150429 fukushima for original section */
	if(getsectionform(flist,(sects+is)->code,(sects+is))==0/* && (sects+is)->code!=(sects+is)->ocode*/)
	{
	  getsectionform(flist,(sects+is)->ocode,(sects+is));
	}
    /***/
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

  cNcr=(double *)malloc(nelem*sizeof(double));   /*GLOBAL VECTOR*/

  if(bclngcondensationflag)
  {
    definencr(&arc,&*cNcr);          /*Buckling Condensation Result*/
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

    if(slabcount==0 && ie<nelem-1) pair=*(elems+ie+1);
    if(slabcount==1)               pair=*(elems+ie-1);

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

        /*soffset=getsectionform(flist,isect,&sect1);*/
        soffset=elem.sect->soff;
        sect1=*(elem.sect);

		/* elem.sect->code=isect; */ // 150429 fukushima for original section

        if(elem.cmqcode!=0) /*CMQ DATA.*/
        {
          for(ic=0;ic<ncmq;ic++)
          {
            if(elem.cmqcode==cmqcode[ic]) elem.Mo=Mo[ic];
          }
        }

        if(soffset && sect1.etype!=WALL
                   && sect1.etype!=SLAB
                   && sect1.etype!=BRACE)
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
          fprintf(fout0,"-------------------------------------");
          fprintf(fout0,"-------------------------------------");
          fprintf(fout0,"-------------------------------------");
          fprintf(fout0,"-----------------\n");
          fprintf(fout0,"ïîçﬁ:%4d ",ielem);
          fprintf(fout0,"éní[:%3d èIí[:%3d ",ni,nj);

          fprintf(fout0,"ífñ :%3d",isect);
          if(sect1.stype     ==STYPE_S)    fprintf(fout0,"=Çr");
          else if(sect1.stype==STYPE_RC)   fprintf(fout0,"=ÇqÇb");
          else if(sect1.stype==STYPE_SRC)  fprintf(fout0,"=ÇrÇqÇb");
          else if(sect1.stype==STYPE_PC)   fprintf(fout0,"=ÇoÇb");
          else if(sect1.stype==STYPE_WOOD) fprintf(fout0,"=ñÿ");
          else if(sect1.stype==STYPE_GLASS) fprintf(fout0,"=ÉKÉâÉX");
          else if(sect1.stype==STYPE_ACRYL) fprintf(fout0,"=ÉAÉNÉäÉã");
          else if(sect1.stype==STYPE_ALUMI) fprintf(fout0,"=ÉAÉãÉ~");
          else                             fprintf(fout0,"=  ");
          if(sect1.etype==COLUMN)      fprintf(fout0,"íå ");
          else if(sect1.etype==GIRDER) fprintf(fout0,"ëÂó¿ ");
          else if(sect1.etype==BEAM)   fprintf(fout0,"è¨ó¿ ");
          else if(sect1.etype==BRACE)  fprintf(fout0,"ãÿà· ");
          else                         fprintf(fout0,"ïsñæ ");

          fprintf(fout0,"çﬁí∑=%.1f[cm] Mxì‡ñ@=%.1f[cm]",h,h0[SX]);
          if(sect1.etype==COLUMN)
          {
            fprintf(fout0," Myì‡ñ@=%.1f[cm]\n",h0[SY]);
          }
          else
          {
            fprintf(fout0,"\n");
          }

          /*MATERIAL INDEPENDENT ON PERIOD.*/
          gmaterial.sE=2100000.0; /*[kgf/cm2]*/
          gmaterial.rE=2100000.0;
          /*gmaterial.cE=210000.0;*/
          gmaterial.gE= 725000.0;
          gmaterial.aE=  32000.0; /*Acryl=15000, Metacryl=32000*/
          gmaterial.cE=gmaterial.rE/15.0;
          gmaterial.alE=700000.0;                     /*ALUMI*/

          gmaterial.sF=2400.0*jis;                    /*STEEL SN400*/
//          gmaterial.alF=1750.0;                       /*ALUMI AS175*/
          gmaterial.alF=4400.0;        /*ALUMI A7178T6 LunarMarsBase*/

		  gmaterial.Fc=240.0;                      /*CONCRETE Fc240*/

  if(!strcmp(prj,"tei")  ||
     !strcmp(prj,"huna") ||
     !strcmp(prj,"kyok"))
  {
    gmaterial.sF=2400.0*jis; /*STEEL SN400A*/
  }

          translatesection(sect1,gmaterial,&As,&Ac,&Ar,&Yg,SX);
          /*fprintf(fout0,"As=%7.2f[cm2] ",As);
          fprintf(fout0,"Ar=%7.2f[cm2] ",Ar);
          fprintf(fout0,"Ac=%7.2f[cm2]\n",Ac);*/

          fprintf(fout0,"âûóÕ       :        N");
          fprintf(fout0,"                Qxi                Qxj");
          fprintf(fout0,"                Qyi                Qyj");
          fprintf(fout0,"                 Mt");
          fprintf(fout0,"                Mxi                Mxj");
          fprintf(fout0,"                Myi                Myj\n");

          fprintf(fout0,"âîíºéûZ    : %8.3f(%8.2f)",
                  elem.head.z.N,si*elem.head.z.N);
          fprintf(fout0," %8.3f(%8.2f) %8.3f(%8.2f)",
                  elem.head.z.Q[0],si*elem.head.z.Q[0],
                  elem.tail.z.Q[0],si*elem.tail.z.Q[0]);
          fprintf(fout0," %8.3f(%8.2f) %8.3f(%8.2f)",
                  elem.head.z.Q[1],si*elem.head.z.Q[1],
                  elem.tail.z.Q[1],si*elem.tail.z.Q[1]);
          fprintf(fout0," %8.3f(%8.2f)",
                  elem.head.z.Mt,si*elem.head.z.Mt);
          fprintf(fout0," %8.3f(%8.2f) %8.3f(%8.2f)",
                  elem.head.z.M[0],si*elem.head.z.M[0],
                  elem.tail.z.M[0],si*elem.tail.z.M[0]);
          fprintf(fout0," %8.3f(%8.2f) %8.3f(%8.2f)\n",
                  elem.head.z.M[1],si*elem.head.z.M[1],
                  elem.tail.z.M[1],si*elem.tail.z.M[1]);

          fprintf(fout0,"êÖïΩéûX    : %8.3f(%8.2f)",
                  elem.head.x.N,si*elem.head.x.N);
          fprintf(fout0," %8.3f(%8.2f) %8.3f(%8.2f)",
                  elem.head.x.Q[0],si*elem.head.x.Q[0],
                  elem.tail.x.Q[0],si*elem.tail.x.Q[0]);
          fprintf(fout0," %8.3f(%8.2f) %8.3f(%8.2f)",
                  elem.head.x.Q[1],si*elem.head.x.Q[1],
                  elem.tail.x.Q[1],si*elem.tail.x.Q[1]);
          fprintf(fout0," %8.3f(%8.2f)",
                  elem.head.x.Mt,si*elem.head.x.Mt);
          fprintf(fout0," %8.3f(%8.2f) %8.3f(%8.2f)",
                  elem.head.x.M[0],si*elem.head.x.M[0],
                  elem.tail.x.M[0],si*elem.tail.x.M[0]);
          fprintf(fout0," %8.3f(%8.2f) %8.3f(%8.2f)\n",
                  elem.head.x.M[1],si*elem.head.x.M[1],
                  elem.tail.x.M[1],si*elem.tail.x.M[1]);

          fprintf(fout0,"êÖïΩéûY    : %8.3f(%8.2f)",
                  elem.head.y.N,si*elem.head.y.N);
          fprintf(fout0," %8.3f(%8.2f) %8.3f(%8.2f)",
                  elem.head.y.Q[0],si*elem.head.y.Q[0],
                  elem.tail.y.Q[0],si*elem.tail.y.Q[0]);
          fprintf(fout0," %8.3f(%8.2f) %8.3f(%8.2f)",
                  elem.head.y.Q[1],si*elem.head.y.Q[1],
                  elem.tail.y.Q[1],si*elem.tail.y.Q[1]);
          fprintf(fout0," %8.3f(%8.2f)",
                  elem.head.y.Mt,si*elem.head.y.Mt);
          fprintf(fout0," %8.3f(%8.2f) %8.3f(%8.2f)",
                  elem.head.y.M[0],si*elem.head.y.M[0],
                  elem.tail.y.M[0],si*elem.tail.y.M[0]);
          fprintf(fout0," %8.3f(%8.2f) %8.3f(%8.2f)\n",
                  elem.head.y.M[1],si*elem.head.y.M[1],
                  elem.tail.y.M[1],si*elem.tail.y.M[1]);
          fprintf(fout0,"\n");

          /*if(check){gets(non); if(non[0]!='\0') return 0;}*/
///SKIP PLONG///////////////////////////////////////////////////////////////////
if(!strcmp(prj,"ghouse")||!strcmp(prj,"lunarbase")||!strcmp(prj,"plywood"))
{
}

else
{

          /*LONG...................................................*/
          /*MATERIAL FOR LONG*/
          gmaterial.sft=2400.0/1.5*jis;               /*STEEL SN400*/
          gmaterial.sfc=-gmaterial.sft; /*FOR SRC*/
          gmaterial.sfb=2400.0/1.5*jis; /*FOR SRC*/     /*b:BENDING*/
          gmaterial.sfs=2400.0/1.5/sqrt(3.0)*jis;         /*s:SHEAR*/
          gmaterial.rft=2200.0*jis;           /*REINFORCEMENT SD345*/
          gmaterial.wft=2000.0;               /*REINFORCEMENT SD295*/

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

  if(!strcmp(prj,"aka")) /*FACE RATES.*/
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
          fprintf(fout0,"í∑ä˙       : %8.3f(%8.2f)",
                  elem.head.z.N,si*elem.head.z.N);
          fprintf(fout0," %8.3f(%8.2f) %8.3f(%8.2f)",
                  elem.head.z.Q[0],si*elem.head.z.Q[0],
                  elem.tail.z.Q[0],si*elem.tail.z.Q[0]);
          fprintf(fout0," %8.3f(%8.2f) %8.3f(%8.2f)",
                  elem.head.z.Q[1],si*elem.head.z.Q[1],
                  elem.tail.z.Q[1],si*elem.tail.z.Q[1]);
          fprintf(fout0," %8.3f(%8.2f)",
                  elem.head.z.Mt,si*elem.head.z.Mt);
          fprintf(fout0," %8.3f(%8.2f) %8.3f(%8.2f)",
                     face[0][0]*elem.head.z.M[0],
                  si*face[0][0]*elem.head.z.M[0],
                     face[1][0]*elem.tail.z.M[0],
                  si*face[1][0]*elem.tail.z.M[0]);
          fprintf(fout0," %8.3f(%8.2f) %8.3f(%8.2f)\n",
                     face[0][1]*elem.head.z.M[1],
                  si*face[0][1]*elem.head.z.M[1],
                     face[1][1]*elem.tail.z.M[1],
                  si*face[1][1]*elem.tail.z.M[1]);

          for(ii=0;ii<=4;ii++)
          {
            for(jj=0;jj<=1;jj++) sprintf(strQ[ii][jj],"\0");
          }
          for(ii=0;ii<=2;ii++)
          {
            for(jj=0;jj<=1;jj++) sprintf(strM[ii][jj],"\0");
          }


          if((sect1.stype==STYPE_S ||
              sect1.stype==STYPE_GLASS ||
              sect1.stype==STYPE_ACRYL ||
              sect1.stype==STYPE_ALUMI)
             && sect1.sform.type==STEEL_PLATE) /*STEEL FB*/
          {
            pstress[0]=&(elem.head.z); pstress[1]=&(elem.tail.z);

            /*HEAD*/
            if(sect1.stype==STYPE_S && !bclngcondensationflag)
            {
              rate=allowablestressofflatbar(PLONG,sect1,gmaterial.sE,
                                            h,pstress[0],HEAD,
                                            &Ncr,&Qcrx,&Qcry,&Mcrx,&Mcry);
            }
            if(sect1.stype==STYPE_S && bclngcondensationflag)
            {
              if(cNcr[ie-1]>0)
              {
              rate=allowablestressofflatbar_bc(PLONG,sect1,cNcr[ie-1],gmaterial.sE,
                                            h,pstress[0],HEAD,
                                            &Ncr,&Qcrx,&Qcry,&Mcrx,&Mcry);
              }
              else
              {
              rate=allowablestressofflatbar(PLONG,sect1,gmaterial.sE,
                                            h,pstress[0],HEAD,
                                            &Ncr,&Qcrx,&Qcry,&Mcrx,&Mcry);
              }
            }
            if(sect1.stype==STYPE_GLASS)
            {
              rate=allowablestressofflatbar(PLONG,sect1,gmaterial.gE,
                                            h,pstress[0],HEAD,
                                            &Ncr,&Qcrx,&Qcry,&Mcrx,&Mcry);
            }
            if(sect1.stype==STYPE_ACRYL)
            {
              rate=allowablestressofflatbar(PLONG,sect1,gmaterial.aE,
                                            h,pstress[0],HEAD,
                                            &Ncr,&Qcrx,&Qcry,&Mcrx,&Mcry);
            }
            if(sect1.stype==STYPE_ALUMI)
            {
              rate=allowablestressofflatbaralumi(PLONG,sect1,gmaterial.alE,
												 h,pstress[0],HEAD,
                                                 &Ncr,&Qcrx,&Qcry,&Mcrx,&Mcry);
            }

/*if(!strcmp(prj,"izu") && (sect1.code==901 ||
                          sect1.code==904)) rate=0.0;*/ /*PASS LONG CALCULATION.*/

if(!strcmp(prj,"vene")) rate=0.0; /*PASS LONG CALCULATION.*/
if(!strcmp(prj,"machi")) rate=0.0; /*PASS LONG CALCULATION.*/
if(!strcmp(prj,"stick")) rate=0.0; /*PASS LONG CALCULATION.*/

            if(rate>ratema) ratema=rate;
            if(rate>marate[soffset-1]) marate[soffset-1]=rate;

            sprintf(str," %8.3f(%8.2f)",Ncr, si*Ncr);
            strcat(strQ[1][0],str);
            sprintf(str," %8.3f(%8.2f)",Qcrx,si*Qcrx);
            strcat(strQ[1][0],str);
            sprintf(str," %8.3f(%8.2f)",Qcry,si*Qcry);
            strcat(strQ[1][1],str);
            sprintf(str," %8.3f(%8.2f)",Mcrx,si*Mcrx);
            strcat(strM[1][0],str);
            sprintf(str," %8.3f(%8.2f)",Mcry,si*Mcry);
            strcat(strM[1][1],str);

			if(Ncr!=0.0) sprintf(str," %8.3f",pstress[0]->N/Ncr);
			else         sprintf(str,"   10.000");
if(!strcmp(prj,"vene")) sprintf(str,"    0.000"); /*PASS LONG CALCULATION.*/
			strcat(strQ[2][0],str);
			if(Qcrx!=0.0) sprintf(str,"           %8.3f",pstress[0]->Q[0]/Qcrx);
			else          sprintf(str,"   10.000");
if(!strcmp(prj,"vene")) sprintf(str,"    0.000"); /*PASS LONG CALCULATION.*/
			strcat(strQ[2][0],str);
			if(Qcry!=0.0) sprintf(str,"           %8.3f",pstress[0]->Q[1]/Qcry);
			else          sprintf(str,"   10.000");
if(!strcmp(prj,"vene")) sprintf(str,"    0.000"); /*PASS LONG CALCULATION.*/
			strcat(strQ[2][1],str);
			if(Mcrx!=0.0) sprintf(str,"           %8.3f",pstress[0]->M[0]/Mcrx);
			else          sprintf(str,"   10.000");
if(!strcmp(prj,"vene")) sprintf(str,"    0.000"); /*PASS LONG CALCULATION.*/
			strcat(strM[2][0],str);
			if(Mcry!=0.0) sprintf(str,"           %8.3f",pstress[0]->M[1]/Mcry);
			else          sprintf(str,"   10.000");
if(!strcmp(prj,"vene")) sprintf(str,"    0.000"); /*PASS LONG CALCULATION.*/
			strcat(strM[2][1],str);

			/*TAIL*/
			if(sect1.stype==STYPE_S && !bclngcondensationflag)
			{
			  rate=allowablestressofflatbar(PLONG,sect1,gmaterial.sE,
											h,pstress[1],TAIL,
                                            &Ncr,&Qcrx,&Qcry,&Mcrx,&Mcry);
            }
            if(sect1.stype==STYPE_S && bclngcondensationflag)
            {
              if(cNcr[ie-1]>0)
              {
              rate=allowablestressofflatbar_bc(PLONG,sect1,cNcr[ie-1],gmaterial.sE,
                                            h,pstress[1],TAIL,
                                            &Ncr,&Qcrx,&Qcry,&Mcrx,&Mcry);
              }
              else
              {
              rate=allowablestressofflatbar(PLONG,sect1,gmaterial.sE,
                                            h,pstress[0],HEAD,
                                            &Ncr,&Qcrx,&Qcry,&Mcrx,&Mcry);
              }
            }
            if(sect1.stype==STYPE_GLASS)
            {
              rate=allowablestressofflatbar(PLONG,sect1,gmaterial.gE,
                                            h,pstress[1],TAIL,
                                            &Ncr,&Qcrx,&Qcry,&Mcrx,&Mcry);
            }
            if(sect1.stype==STYPE_ACRYL)
            {
              rate=allowablestressofflatbar(PLONG,sect1,gmaterial.aE,
                                            h,pstress[1],TAIL,
                                            &Ncr,&Qcrx,&Qcry,&Mcrx,&Mcry);
            }
            if(sect1.stype==STYPE_ALUMI)
            {
              rate=allowablestressofflatbaralumi(PLONG,sect1,gmaterial.alE,
												 h,pstress[0],HEAD,
                                                 &Ncr,&Qcrx,&Qcry,&Mcrx,&Mcry);
            }

/*if(!strcmp(prj,"izu") && (sect1.code==901 ||
                          sect1.code==904)) rate=0.0;*/ /*PASS LONG CALCULATION.*/

if(!strcmp(prj,"vene")) rate=0.0; /*PASS LONG CALCULATION.*/
if(!strcmp(prj,"machi")) rate=0.0; /*PASS LONG CALCULATION.*/
if(!strcmp(prj,"stick")) rate=0.0; /*PASS LONG CALCULATION.*/

            if(rate>ratema) ratema=rate;
            if(rate>marate[soffset-1]) marate[soffset-1]=rate;

            rateqa=ratema;
            qarate[soffset-1]=marate[soffset-1];

            sprintf(str," %8.3f(%8.2f)",Qcrx,si*Qcrx);
            strcat(strQ[1][0],str);
            sprintf(str," %8.3f(%8.2f)",Qcry,si*Qcry);
            strcat(strQ[1][1],str);
            sprintf(str," %8.3f(%8.2f)",Mcrx,si*Mcrx);
            strcat(strM[1][0],str);
            sprintf(str," %8.3f(%8.2f)",Mcry,si*Mcry);
            strcat(strM[1][1],str);

            if(Qcrx!=0.0) sprintf(str,"           %8.3f",pstress[1]->Q[0]/Qcrx);
            else          sprintf(str,"   10.000");
if(!strcmp(prj,"vene")) sprintf(str,"    0.000"); /*PASS LONG CALCULATION.*/
            strcat(strQ[2][0],str);
            if(Qcry!=0.0) sprintf(str,"           %8.3f",pstress[1]->Q[1]/Qcry);
            else          sprintf(str,"   10.000");
if(!strcmp(prj,"vene")) sprintf(str,"    0.000"); /*PASS LONG CALCULATION.*/
            strcat(strQ[2][1],str);
            if(Mcrx!=0.0) sprintf(str,"           %8.3f",pstress[1]->M[0]/Mcrx);
            else          sprintf(str,"   10.000");
if(!strcmp(prj,"vene")) sprintf(str,"    0.000"); /*PASS LONG CALCULATION.*/
            strcat(strM[2][0],str);
            if(Mcry!=0.0) sprintf(str,"           %8.3f",pstress[1]->M[1]/Mcry);
            else          sprintf(str,"   10.000");
if(!strcmp(prj,"vene")) sprintf(str,"    0.000"); /*PASS LONG CALCULATION.*/
            strcat(strM[2][1],str);

            /*OUTPUT*/
            fprintf(fout0,"     ãñóeíl:%s%s",
                    strQ[1][0],strQ[1][1]);
            fprintf(fout0,"                   %s%s\n",
                    strM[1][0],strM[1][1]);
            fprintf(fout0,"     à¿ëSó¶:%s%s",
                    strQ[2][0],strQ[2][1]);
            fprintf(fout0,"                   %s%s\n",
                    strM[2][0],strM[2][1]);
          }
          else if(sect1.stype==STYPE_WOOD
                  && sect1.sform.type==DABO) /*WOOD DABO*/
          {
            pstress[0]=&(elem.head.z); pstress[1]=&(elem.tail.z);

            /*HEAD*/
            rate=allowablestressofdabo(PLONG,sect1,pstress[0],HEAD,
                                       &Ncr,&Qcrx,&Mcrx);

            if(rate>ratema) ratema=rate;
            if(rate>marate[soffset-1]) marate[soffset-1]=rate;

            sprintf(str," %8.3f(%8.2f)",Ncr, si*Ncr);
            strcat(strQ[1][0],str);
            sprintf(str," %8.3f(%8.2f)",Qcrx,si*Qcrx);
            strcat(strQ[1][0],str);
            sprintf(str," %8.3f(%8.2f)",Qcrx,si*Qcrx);
            strcat(strQ[1][1],str);
            sprintf(str," %8.3f(%8.2f)",Mcrx,si*Mcrx);
            strcat(strM[1][0],str);
            sprintf(str," %8.3f(%8.2f)",Mcrx,si*Mcrx);
            strcat(strM[1][1],str);

            if(Ncr!=0.0) sprintf(str," %8.3f",pstress[0]->N/Ncr);
            else         sprintf(str,"   10.000");
            strcat(strQ[2][0],str);
            if(Qcrx!=0.0) sprintf(str,"           %8.3f",pstress[0]->Q[0]/Qcrx);
            else          sprintf(str,"   10.000");
            strcat(strQ[2][0],str);
            if(Qcry!=0.0) sprintf(str,"           %8.3f",pstress[0]->Q[1]/Qcrx);
            else          sprintf(str,"   10.000");
            strcat(strQ[2][1],str);
            if(Mcrx!=0.0) sprintf(str,"           %8.3f",pstress[0]->M[0]/Mcrx);
            else          sprintf(str,"   10.000");
            strcat(strM[2][0],str);
            if(Mcry!=0.0) sprintf(str,"           %8.3f",pstress[0]->M[1]/Mcrx);
            else          sprintf(str,"   10.000");
            strcat(strM[2][1],str);

            /*TAIL*/
            rate=allowablestressofdabo(PLONG,sect1,pstress[1],TAIL,
                                       &Ncr,&Qcrx,&Mcrx);

            if(rate>ratema) ratema=rate;
            if(rate>marate[soffset-1]) marate[soffset-1]=rate;

            rateqa=ratema;
            qarate[soffset-1]=marate[soffset-1];

            sprintf(str," %8.3f(%8.2f)",Qcrx,si*Qcrx);
            strcat(strQ[1][0],str);
            sprintf(str," %8.3f(%8.2f)",Qcrx,si*Qcrx);
            strcat(strQ[1][1],str);
            sprintf(str," %8.3f(%8.2f)",Mcrx,si*Mcrx);
            strcat(strM[1][0],str);
            sprintf(str," %8.3f(%8.2f)",Mcrx,si*Mcrx);
            strcat(strM[1][1],str);

            if(Qcrx!=0.0) sprintf(str,"           %8.3f",pstress[1]->Q[0]/Qcrx);
            else          sprintf(str,"   10.000");
            strcat(strQ[2][0],str);
            if(Qcry!=0.0) sprintf(str,"           %8.3f",pstress[1]->Q[1]/Qcrx);
            else          sprintf(str,"   10.000");
            strcat(strQ[2][1],str);
            if(Mcrx!=0.0) sprintf(str,"           %8.3f",pstress[1]->M[0]/Mcrx);
            else          sprintf(str,"   10.000");
            strcat(strM[2][0],str);
            if(Mcry!=0.0) sprintf(str,"           %8.3f",pstress[1]->M[1]/Mcrx);
            else          sprintf(str,"   10.000");
            strcat(strM[2][1],str);

            /*OUTPUT*/
            fprintf(fout0,"     ãñóeíl:%s%s",
                    strQ[1][0],strQ[1][1]);
            fprintf(fout0,"                   %s%s\n",
                    strM[1][0],strM[1][1]);
            fprintf(fout0,"     à¿ëSó¶:%s%s",
                    strQ[2][0],strQ[2][1]);
            fprintf(fout0,"                   %s%s\n",
                    strM[2][0],strM[2][1]);
          }
          else
          {
            if(elem.sect->etype==COLUMN) na=1;
            else                         na=0;
            for(ia=0;ia<=na;ia++) /*FOR AXIS*/
            {
              pstress[0]=&(elem.head.z); pstress[1]=&(elem.tail.z);

              lN=1000.0*elem.head.z.N; /*LONG N[kgf].*/

              Ma=0.0;

              if(sect1.stype==STYPE_S)                            /*S*/
              {
                Ma=allowablebendingofsteel(PLONG,sect1,gmaterial.sE,
                                           h,ia,(-lN));
              }
              else if(sect1.stype==STYPE_RC)                     /*RC*/
              {
                Ma=allowablebendingofrc(elem,gmaterial,ia,(-lN),PLONG);
              }
              else if(sect1.stype==STYPE_SRC)                   /*SRC*/
              {
                if(sect1.sform.type==STEEL_RECTS)
                {
                  Ma=allowablebendingofsrc(elem,gmaterial,ia,(-lN),PLONG);
                }
              }
              else if(sect1.stype==STYPE_PC)                     /*PC*/
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
              else if(sect1.stype==STYPE_WOOD)                 /*WOOD*/
              {
                Ma=allowablebendingofwood(PLONG,sect1,h,h,ia,(-lN));
              }
              else if(sect1.stype==STYPE_ALUMI)               /*ALUMI*/
              {
				Ma=allowablebendingofalumi(PLONG,sect1,gmaterial.alE,
										   h,ia,(-lN));
				/*
                Mahead=allowablebendingofalumi(PLONG,sect1,gmaterial.alE,
										   h,ia,(-lN));
				Matail=Mahead;
                */
			  }
              /*if(Ma<0.1) Ma=0.1;*/

              sprintf(str," %8.3f(%8.2f) %8.3f(%8.2f)",
                      Ma/100000.0,si*Ma/100000.0,
                      Ma/100000.0,si*Ma/100000.0);
              strcat(strM[1][ia],str);
              if(Ma>0.0)
              {
                rate1=fabs(face[0][ia]*pstress[0]->M[ia]*100000.0/Ma);
                rate2=fabs(face[1][ia]*pstress[1]->M[ia]*100000.0/Ma);
              }
              else
              {
                rate1=10.0;
                rate2=10.0;
              }

if(!strcmp(prj,"vene")){rate1=0.0; rate2=0.0;} /*PASS LONG CALCULATION.*/
if(!strcmp(prj,"machi")){rate1=0.0; rate2=0.0;} /*PASS LONG CALCULATION.*/
if(!strcmp(prj,"stick")){rate1=0.0; rate2=0.0;} /*PASS LONG CALCULATION.*/

              sprintf(str,"           %8.3f           %8.3f",
                      rate1,rate2);
              strcat(strM[2][ia],str);

              /*
              fprintf(fout0,"Ma%s=%.3f[tfm] ",straxis[ia],Ma/100000.0);
              fprintf(fout0," Mi/Ma=%.3f Mj/Ma=%.3f\n",
                      fabs(pstress[0]->M[ia]*100000.0/Ma),
                      fabs(pstress[1]->M[ia]*100000.0/Ma));
              */

              if(rate1>ratema) ratema=rate1;
              if(rate1>marate[soffset-1]) marate[soffset-1]=rate1;
              if(rate2>ratema) ratema=rate2;
              if(rate2>marate[soffset-1]) marate[soffset-1]=rate2;

              for(in=HEAD;in<=TAIL;in++) /*FOR END*/
              {
                Qa=0.0;
                if(sect1.stype==STYPE_S)                          /*S*/
                {
                  Qa=allowultimshearofsteel(PLONG,gmaterial.sF,sect1,ia);
                }
                else if(sect1.stype==STYPE_RC)                   /*RC*/
                {
                  Qa=allowableshearofrclong(elem,ia,
                                            pstress[in]->Q[1-ia]
                                            *1000.0,
                                            pstress[in]->M[ia]
                                            *100000.0);
                }
                else if(sect1.stype==STYPE_SRC)                 /*SRC*/
                {
                  if(sect1.sform.type==STEEL_RECTS && ia==SX)
                  {
                    Qa=allowableshearofsrclong(elem,ia,HSTRONG,
                                               pstress[in]->Q[1-ia]
                                               *1000.0,
                                               pstress[in]->M[ia]
                                               *100000.0);
                  }
                  if(sect1.sform.type==STEEL_RECTS && ia==SY)
                  {
                    Qa=allowableshearofsrclong(elem,ia,HWEAK,
                                               pstress[in]->Q[1-ia]
                                               *1000.0,
                                               pstress[in]->M[ia]
                                               *100000.0);
                  }
                }
                else if(sect1.stype==STYPE_PC)                   /*PC*/
                {
                  Qa=allowableshearofpccolumn(elem,gmaterial,ia,lN);
                }
                else if(sect1.stype==STYPE_WOOD)               /*WOOD*/
                {
                  Qa=allowableshearofwood(PLONG,sect1);
                }
                else if(sect1.stype==STYPE_ALUMI)             /*ALUMI*/
                {
                  if(sect1.sform.type==STEEL_RECTS) gmaterial.alF=sect1.srect[0].F;
                  else                              gmaterial.alF=sect1.sform.F;
                  Qa=allowultimshearofalumi(PLONG,gmaterial.alF,sect1,ia);
                }

                sprintf(str," %8.3f(%8.2f)",
                        Qa/1000.0,si*Qa/1000.0);
                strcat(strQ[1][1-ia],str);

                /*fprintf(fout0,"Qa%s=%.3f[tf]",straxis[1-ia],Qa/1000.0);*/

                /*if(Qa<0.1) Qa=0.1;*/

				if(Qa>0.0) rate=fabs(pstress[in]->Q[1-ia]*1000.0/Qa);
				else       rate=10.0;

if(!strcmp(prj,"vene")) rate=0.0; /*PASS LONG CALCULATION.*/
if(!strcmp(prj,"machi")) rate=0.0; /*PASS LONG CALCULATION.*/
if(!strcmp(prj,"stick")) rate=0.0; /*PASS LONG CALCULATION.*/

				if(rate>rateqa) rateqa=rate;
                if(rate>qarate[soffset-1]) qarate[soffset-1]=rate;

                sprintf(str,"           %8.3f",rate);
                strcat(strQ[2][1-ia],str);

                /*fprintf(fout0," Q%s/Qa=%.3f\n",strend[in],rate);*/
              }
            }
            if(sect1.etype==GIRDER || sect1.etype==BEAM)
            {
              sprintf(strQ[1][0],"                                      ");
              sprintf(strQ[2][0],"                                      ");
              sprintf(strM[1][1],"                              ");
              sprintf(strM[2][1],"                              ");
            }

            fprintf(fout0,"     ãñóeíl:                   %s%s",
                    strQ[1][0],strQ[1][1]);
            fprintf(fout0,"                   %s%s\n",
                    strM[1][0],strM[1][1]);
            fprintf(fout0,"     à¿ëSó¶:         %s%s",
                    strQ[2][0],strQ[2][1]);
            fprintf(fout0,"                   %s%s\n",
                    strM[2][0],strM[2][1]);

            /*FOR MIDPOINT BENDING Mm.*/
            if(sect1.stype==STYPE_PC &&
               elem.sect->etype==GIRDER &&
               elem.cmqcode!=0)
            {
              Mm=elem.Mo-0.5*(elem.tail.z.M[0]-elem.head.z.M[0]);
              fprintf(fout0,"        íÜâõ:Mo=%8.3f Mm=%8.3f",
                      elem.Mo,Mm);

              if(Mm<=0.0)
              {
                fprintf(fout0," è„í[à¯í£\n\n\n");
              }
              else
              {
                Ma=allowablebendingofpcgirderonmid(elem,gmaterial,
                                                   h0[SX],0);
                /*if(Ma<0.1) Ma=0.1;*/

                if(Ma>0.0) rate=fabs(Mm*100000.0/Ma);
                else       rate=10.0;

                if(rate>ratema) ratema=rate;
                if(rate>marate[soffset-1]) marate[soffset-1]=rate;

                fprintf(fout0,"\n");
                fprintf(fout0,"      ãñóeíl:               %8.3f\n",
                        Ma/100000.0);
                fprintf(fout0,"      à¿ëSó¶:               %8.3f\n",
                        rate);
              }
            }
            else if(elem.sect->etype==GIRDER)
            {
              /*fprintf(fout0,"\n\n\n");*/
            }
          }

          fprintf(fout0,"\n");
          /*if(check){gets(non); if(non[0]!='\0') return 0;}*/
}    /*default*/

          /*SHORT..................................................*/
          /*MATERIAL FOR SHORT*/
          gmaterial.sft=2400.0*jis;                   /*STEEL SN400*/
          gmaterial.sfc=-gmaterial.sft;
          gmaterial.sfb=2400.0*jis; /*FOR SRC*/         /*b:BENDING*/
          gmaterial.sfs=2400.0/sqrt(3.0)*jis;             /*s:SHEAR*/
          gmaterial.rft=3500.0*jis;           /*REINFORCEMENT SD345*/
          gmaterial.wft=3000.0;               /*REINFORCEMENT SD295*/

          gmaterial.rfc=-gmaterial.rft;

          /*n=15.0;*/                                  /*RATE Er/Ec*/
          gmaterial.cfc=-160.0;                          /*CONCRETE*/
          gmaterial.cfs=11.1;

          for(hoko=0;hoko<=1;hoko++) /*LOAD DIRECTION X,Y.*/
          {
            if(sect1.stype==STYPE_PC) break;  /*PC ROUTE 3a : SKIP SHORT.*/
			if(!strcmp(prj,"kirigami")) break;  /*kirigami : SKIP SHORT.*/

/*if(!strcmp(prj,"got") && hoko==0)
{
  nfact=1.0;
  mfact=1.0*(1/3.0)/0.086936;
  qfact=2.0*(1/3.0)/0.086936;
}
else if(!strcmp(prj,"got") && hoko==1)
{
  nfact=1.0;
  mfact=1.0*(1/3.0)/0.042216;
  qfact=2.0*(1/3.0)/0.042216;
}
else
{
  nfact=1.0;
  mfact=1.0;
  qfact=2.0;
}*/

if(!strcmp(prj,"aii") && (sect1.code==201 ||
                          sect1.code==202 ||
                          sect1.code==203))
{
  nfact=1.0;
  mfact=1.0;
  qfact=2.0;
}

if(!strcmp(prj,"tos") && hoko==0)
{
  nfact=1.0;
  mfact=1.0;
  qfact=2.0;
}
else if(!strcmp(prj,"tos") && hoko==1)
{
  nfact=1.5;
  mfact=1.5;
  qfact=2.0;
}

if(!strcmp(prj,"ashi") && hoko==0)
{
  nfact=1.5;
  mfact=1.5;
  qfact=2.25;
}
else if(!strcmp(prj,"ashi") && hoko==1)
{
  nfact=1.0;
  mfact=1.0;
  qfact=2.0;
}

if(!strcmp(prj,"banga") && sect1.sform.type==DABO)
{
  nfact=1.0;
  mfact=1.0;
  qfact=1.0;
}
else
{
  nfact=1.0;
  mfact=1.0;
  qfact=2.0;
}

if(!strcmp(prj,"kamime"))
{
  nfact=1.0;
  mfact=1.0;
  qfact=1.0;
}

if(!strcmp(prj,"vene"))
{
  nfact=1.0;
  mfact=1.0;
  qfact=1.0;
}
#if 1
if(!strcmp(prj,"vene") &&
   globalstatus==PINCHNODEVENEZIALONG) /*PASS X,Y*/
{
  nfact=0.0;
  mfact=0.0;
  qfact=0.0;
}
if(!strcmp(prj,"vene") && hoko==0 &&
   (globalstatus==PINCHNODEVENEZIAYPOSI ||
    globalstatus==PINCHNODEVENEZIAYNEGA)) /*PASS X*/
{
  nfact=0.0;
  mfact=0.0;
  qfact=0.0;
}
if(!strcmp(prj,"vene") && hoko==1 &&
   (globalstatus==PINCHNODEVENEZIAXPOSI ||
    globalstatus==PINCHNODEVENEZIAXNEGA)) /*PASS Y*/
{
  nfact=0.0;
  mfact=0.0;
  qfact=0.0;
}

if(!strcmp(prj,"maebasi"))
{
  nfact=1.0;
  mfact=1.0;
  qfact=1.5;
}
#endif

            if(hoko==0)
            {
              pstress[0]=&(elem.head.x); pstress[1]=&(elem.tail.x);
            }
            if(hoko==1)
            {
              pstress[0]=&(elem.head.y); pstress[1]=&(elem.tail.y);
            }

            for(sign=0;sign<=1;sign++) /*LOAD POSITIVE,NEGATIVE.*/
//            for(sign=1;sign<=1;sign++) /*LOAD ONLY POSITIVE/NEGATIVE.*/ //UJIOKA for wakka
            {
              for(ii=0;ii<=4;ii++)
              {
                for(jj=0;jj<=1;jj++) sprintf(strQ[ii][jj],"\0");
              }
              for(ii=0;ii<=2;ii++)
              {
                for(jj=0;jj<=1;jj++) sprintf(strM[ii][jj],"\0");
              }

              sprintf(strQ[0][0],"íZä˙%s%sï˚å¸:",
                      strhoko[hoko],strhugo[sign]);

              if((sect1.stype==STYPE_S ||
                  sect1.stype==STYPE_GLASS ||
                  sect1.stype==STYPE_ACRYL ||
                  sect1.stype==STYPE_ALUMI)
                 && sect1.sform.type==STEEL_PLATE) /*FB*/
              {
                if(sign==0) dsign=1.0; /*POSITIVE*/
                if(sign==1) dsign=-1.0;/*NEGATIVE*/

#if 1
if(!strcmp(prj,"vene") &&
   (globalstatus==PINCHNODEVENEZIAXPOSI ||
    globalstatus==PINCHNODEVENEZIAYPOSI))
{if(sign==1) dsign=1.0;} /*PASS NEGATIVE*/
if(!strcmp(prj,"vene") &&
   (globalstatus==PINCHNODEVENEZIAXNEGA ||
    globalstatus==PINCHNODEVENEZIAYNEGA))
{if(sign==1) dsign=1.0;} /*PASS NEGATIVE*/
#endif

                /*HEAD*/
                estress.N   =elem.head.z.N
                             +dsign*nfact*(pstress[HEAD]->N);
                estress.Q[0]=elem.head.z.Q[0]
                             +dsign*qfact*(pstress[HEAD]->Q[0]);
                estress.Q[1]=elem.head.z.Q[1]
                             +dsign*qfact*(pstress[HEAD]->Q[1]);
                estress.Mt  =elem.head.z.Mt
                             +dsign*mfact*(pstress[HEAD]->Mt);
                estress.M[0]=elem.head.z.M[0]
                             +dsign*mfact*(pstress[HEAD]->M[0]);
                estress.M[1]=elem.head.z.M[1]
                             +dsign*mfact*(pstress[HEAD]->M[1]);

                sprintf(str," %8.3f(%8.2f)",estress.N,si*estress.N);
                strcat(strQ[0][0],str);
                sprintf(str," %8.3f(%8.2f)",estress.Q[0],si*estress.Q[0]);
                strcat(strQ[0][0],str);
                sprintf(str," %8.3f(%8.2f)",estress.Q[1],si*estress.Q[1]);
                strcat(strQ[0][1],str);
                sprintf(str," %8.3f(%8.2f)",estress.Mt,si*estress.Mt);
                strcat(strM[0][0],str);
                sprintf(str," %8.3f(%8.2f)",estress.M[0],si*estress.M[0]);
                strcat(strM[0][0],str);
                sprintf(str," %8.3f(%8.2f)",estress.M[1],si*estress.M[1]);
                strcat(strM[0][1],str);

                if(sect1.stype==STYPE_S && !bclngcondensationflag)
                {
                  rate=allowablestressofflatbar(PSHORT,sect1,gmaterial.sE,
                                                h,&estress,HEAD,
                                                &Ncr,&Qcrx,&Qcry,&Mcrx,&Mcry);
                }
                if(sect1.stype==STYPE_S && bclngcondensationflag)
                {
                  if(cNcr[ie-1]>0)
                  {
                  rate=allowablestressofflatbar_bc(PSHORT,sect1,cNcr[ie-1],gmaterial.sE,
                                                h,&estress,HEAD,
                                                &Ncr,&Qcrx,&Qcry,&Mcrx,&Mcry);
                  }
                  else
                  {
                  rate=allowablestressofflatbar(PLONG,sect1,gmaterial.sE,
                                                h,pstress[0],HEAD,
                                                &Ncr,&Qcrx,&Qcry,&Mcrx,&Mcry);
                  }
                }
                if(sect1.stype==STYPE_GLASS)
                {
                  rate=allowablestressofflatbar(PSHORT,sect1,gmaterial.gE,
                                                h,&estress,HEAD,
                                                &Ncr,&Qcrx,&Qcry,&Mcrx,&Mcry);
                }
                if(sect1.stype==STYPE_ACRYL)
                {
                  rate=allowablestressofflatbar(PSHORT,sect1,gmaterial.aE,
                                                h,&estress,HEAD,
                                                &Ncr,&Qcrx,&Qcry,&Mcrx,&Mcry);
                }
                if(sect1.stype==STYPE_ALUMI)
                {
                  rate=allowablestressofflatbaralumi(PSHORT,sect1,gmaterial.alE,
                                                     h,&estress,HEAD,
                                                     &Ncr,&Qcrx,&Qcry,&Mcrx,&Mcry);
                }


                if(rate>ratema) ratema=rate;
                if(rate>marate[soffset-1]) marate[soffset-1]=rate;

                sprintf(str," %8.3f(%8.2f)",Ncr, si*Ncr);
                strcat(strQ[1][0],str);
                sprintf(str," %8.3f(%8.2f)",Qcrx,si*Qcrx);
                strcat(strQ[1][0],str);
                sprintf(str," %8.3f(%8.2f)",Qcry,si*Qcry);
                strcat(strQ[1][1],str);
                sprintf(str," %8.3f(%8.2f)",Mcrx,si*Mcrx);
                strcat(strM[1][0],str);
                sprintf(str," %8.3f(%8.2f)",Mcry,si*Mcry);
                strcat(strM[1][1],str);

                if(Ncr!=0.0) sprintf(str," %8.3f",estress.N/Ncr);
                else         sprintf(str,"   10.000");
                strcat(strQ[2][0],str);
                if(Qcrx!=0.0) sprintf(str,"           %8.3f",estress.Q[0]/Qcrx);
                else          sprintf(str,"   10.000");
                strcat(strQ[2][0],str);
                if(Qcry!=0.0) sprintf(str,"           %8.3f",estress.Q[1]/Qcry);
                else          sprintf(str,"   10.000");
                strcat(strQ[2][1],str);
                if(Mcrx!=0.0) sprintf(str,"           %8.3f",estress.M[0]/Mcrx);
                else          sprintf(str,"   10.000");
                strcat(strM[2][0],str);
                if(Mcry!=0.0) sprintf(str,"           %8.3f",estress.M[1]/Mcry);
                else          sprintf(str,"   10.000");
                strcat(strM[2][1],str);

                /*TAIL*/
                estress.N   =elem.tail.z.N
                             +dsign*nfact*(pstress[TAIL]->N);
                estress.Q[0]=elem.tail.z.Q[0]
                             +dsign*qfact*(pstress[TAIL]->Q[0]);
                estress.Q[1]=elem.tail.z.Q[1]
                             +dsign*qfact*(pstress[TAIL]->Q[1]);
                estress.Mt  =elem.tail.z.Mt
                             +dsign*mfact*(pstress[TAIL]->Mt);
                estress.M[0]=elem.tail.z.M[0]
                             +dsign*mfact*(pstress[TAIL]->M[0]);
                estress.M[1]=elem.tail.z.M[1]
                             +dsign*mfact*(pstress[TAIL]->M[1]);

                sprintf(str," %8.3f(%8.2f)",estress.Q[0],si*estress.Q[0]);
                strcat(strQ[0][0],str);
                sprintf(str," %8.3f(%8.2f)",estress.Q[1],si*estress.Q[1]);
                strcat(strQ[0][1],str);
                sprintf(str," %8.3f(%8.2f)",estress.M[0],si*estress.M[0]);
                strcat(strM[0][0],str);
                sprintf(str," %8.3f(%8.2f)",estress.M[1],si*estress.M[1]);
                strcat(strM[0][1],str);

                if(sect1.stype==STYPE_S && !bclngcondensationflag)
                {
                  rate=allowablestressofflatbar(PSHORT,sect1,gmaterial.sE,
                                                h,&estress,TAIL,
                                                &Ncr,&Qcrx,&Qcry,&Mcrx,&Mcry);
                }
                if(sect1.stype==STYPE_S && bclngcondensationflag)
                {
                  if(cNcr[ie-1]>0)
                  {
                  rate=allowablestressofflatbar_bc(PSHORT,sect1,cNcr[ie-1],gmaterial.sE,
                                                h,&estress,TAIL,
                                                &Ncr,&Qcrx,&Qcry,&Mcrx,&Mcry);
                  }
                  else
                  {
                  rate=allowablestressofflatbar(PLONG,sect1,gmaterial.sE,
                                                h,pstress[0],HEAD,
                                                &Ncr,&Qcrx,&Qcry,&Mcrx,&Mcry);
                  }
                }
                if(sect1.stype==STYPE_GLASS)
                {
                  rate=allowablestressofflatbar(PSHORT,sect1,gmaterial.gE,
                                                h,&estress,TAIL,
                                                &Ncr,&Qcrx,&Qcry,&Mcrx,&Mcry);
                }
                if(sect1.stype==STYPE_ACRYL)
                {
                  rate=allowablestressofflatbar(PSHORT,sect1,gmaterial.aE,
                                                h,&estress,TAIL,
                                                &Ncr,&Qcrx,&Qcry,&Mcrx,&Mcry);
                }
                if(sect1.stype==STYPE_ALUMI)
                {
/*
		            sprintf(str," %8.3f",gmaterial.alE);
        		    MessageBox(NULL,str,"SRCan",MB_OK);
*/                  rate=allowablestressofflatbaralumi(PSHORT,sect1,gmaterial.alE,
                                                     h,&estress,TAIL,
                                                     &Ncr,&Qcrx,&Qcry,&Mcrx,&Mcry);
                }


                if(rate>ratema) ratema=rate;
                if(rate>marate[soffset-1]) marate[soffset-1]=rate;

                rateqa=ratema;
                qarate[soffset-1]=marate[soffset-1];

                sprintf(str," %8.3f(%8.2f)",Qcrx,si*Qcrx);
                strcat(strQ[1][0],str);
                sprintf(str," %8.3f(%8.2f)",Qcry,si*Qcry);
                strcat(strQ[1][1],str);
                sprintf(str," %8.3f(%8.2f)",Mcrx,si*Mcrx);
                strcat(strM[1][0],str);
                sprintf(str," %8.3f(%8.2f)",Mcry,si*Mcry);
                strcat(strM[1][1],str);

                if(Qcrx!=0.0) sprintf(str,"           %8.3f",estress.Q[0]/Qcrx);
                else          sprintf(str,"   10.000");
                strcat(strQ[2][0],str);
                if(Qcry!=0.0) sprintf(str,"           %8.3f",estress.Q[1]/Qcry);
                else          sprintf(str,"   10.000");
                strcat(strQ[2][1],str);
                if(Mcrx!=0.0) sprintf(str,"           %8.3f",estress.M[0]/Mcrx);
                else          sprintf(str,"   10.000");
                strcat(strM[2][0],str);
                if(Mcry!=0.0) sprintf(str,"           %8.3f",estress.M[1]/Mcry);
                else          sprintf(str,"   10.000");
                strcat(strM[2][1],str);

                /*OUTPUT*/
                fprintf(fout0,"%s%s%s%s\n",
                        strQ[0][0],strQ[0][1],strM[0][0],strM[0][1]);

                fprintf(fout0,"     ãñóeíl:%s%s",
                        strQ[1][0],strQ[1][1]);
                fprintf(fout0,"                   %s%s\n",
                        strM[1][0],strM[1][1]);
                fprintf(fout0,"     à¿ëSó¶:%s%s",
                        strQ[2][0],strQ[2][1]);
                fprintf(fout0,"                   %s%s\n",
                        strM[2][0],strM[2][1]);
              }
              else if(sect1.stype==STYPE_WOOD
                      && sect1.sform.type==DABO) /*WOOD DABO*/
              {
                if(sign==0) dsign=1.0; /*POSITIVE*/
                if(sign==1) dsign=-1.0;/*NEGATIVE*/

                /*HEAD*/
                estress.N   =elem.head.z.N
                             +dsign*nfact*(pstress[HEAD]->N);
                estress.Q[0]=elem.head.z.Q[0]
                             +dsign*qfact*(pstress[HEAD]->Q[0]);
                estress.Q[1]=elem.head.z.Q[1]
                             +dsign*qfact*(pstress[HEAD]->Q[1]);
                estress.Mt  =elem.head.z.Mt
                             +dsign*mfact*(pstress[HEAD]->Mt);
                estress.M[0]=elem.head.z.M[0]
                             +dsign*mfact*(pstress[HEAD]->M[0]);
                estress.M[1]=elem.head.z.M[1]
                             +dsign*mfact*(pstress[HEAD]->M[1]);

                sprintf(str," %8.3f(%8.2f)",estress.N,si*estress.N);
                strcat(strQ[0][0],str);
                sprintf(str," %8.3f(%8.2f)",estress.Q[0],si*estress.Q[0]);
                strcat(strQ[0][0],str);
                sprintf(str," %8.3f(%8.2f)",estress.Q[1],si*estress.Q[1]);
                strcat(strQ[0][1],str);
                sprintf(str," %8.3f(%8.2f)",estress.Mt,si*estress.Mt);
                strcat(strM[0][0],str);
                sprintf(str," %8.3f(%8.2f)",estress.M[0],si*estress.M[0]);
                strcat(strM[0][0],str);
                sprintf(str," %8.3f(%8.2f)",estress.M[1],si*estress.M[1]);
                strcat(strM[0][1],str);

                rate=allowablestressofdabo(PSHORT,sect1,&estress,HEAD,
                                           &Ncr,&Qcrx,&Mcrx);

                if(rate>ratema) ratema=rate;
                if(rate>marate[soffset-1]) marate[soffset-1]=rate;

                sprintf(str," %8.3f(%8.2f)",Ncr, si*Ncr);
                strcat(strQ[1][0],str);
                sprintf(str," %8.3f(%8.2f)",Qcrx,si*Qcrx);
                strcat(strQ[1][0],str);
                sprintf(str," %8.3f(%8.2f)",Qcrx,si*Qcrx);
                strcat(strQ[1][1],str);
                sprintf(str," %8.3f(%8.2f)",Mcrx,si*Mcrx);
                strcat(strM[1][0],str);
                sprintf(str," %8.3f(%8.2f)",Mcrx,si*Mcrx);
                strcat(strM[1][1],str);

                if(Ncr!=0.0) sprintf(str," %8.3f",estress.N/Ncr);
                else         sprintf(str,"   10.000");
                strcat(strQ[2][0],str);
                if(Qcrx!=0.0) sprintf(str,"           %8.3f",estress.Q[0]/Qcrx);
                else          sprintf(str,"   10.000");
                strcat(strQ[2][0],str);
                if(Qcry!=0.0) sprintf(str,"           %8.3f",estress.Q[1]/Qcrx);
                else          sprintf(str,"   10.000");
                strcat(strQ[2][1],str);
                if(Mcrx!=0.0) sprintf(str,"           %8.3f",estress.M[0]/Mcrx);
                else          sprintf(str,"   10.000");
                strcat(strM[2][0],str);
                if(Mcry!=0.0) sprintf(str,"           %8.3f",estress.M[1]/Mcrx);
                else          sprintf(str,"   10.000");
                strcat(strM[2][1],str);

                /*TAIL*/
                estress.N   =elem.tail.z.N
                             +dsign*nfact*(pstress[TAIL]->N);
                estress.Q[0]=elem.tail.z.Q[0]
                             +dsign*qfact*(pstress[TAIL]->Q[0]);
                estress.Q[1]=elem.tail.z.Q[1]
                             +dsign*qfact*(pstress[TAIL]->Q[1]);
                estress.Mt  =elem.tail.z.Mt
                             +dsign*mfact*(pstress[TAIL]->Mt);
                estress.M[0]=elem.tail.z.M[0]
                             +dsign*mfact*(pstress[TAIL]->M[0]);
                estress.M[1]=elem.tail.z.M[1]
                             +dsign*mfact*(pstress[TAIL]->M[1]);

                sprintf(str," %8.3f(%8.2f)",estress.Q[0],si*estress.Q[0]);
                strcat(strQ[0][0],str);
                sprintf(str," %8.3f(%8.2f)",estress.Q[1],si*estress.Q[1]);
                strcat(strQ[0][1],str);
                sprintf(str," %8.3f(%8.2f)",estress.M[0],si*estress.M[0]);
                strcat(strM[0][0],str);
                sprintf(str," %8.3f(%8.2f)",estress.M[1],si*estress.M[1]);
                strcat(strM[0][1],str);

                rate=allowablestressofdabo(PSHORT,sect1,&estress,TAIL,
                                           &Ncr,&Qcrx,&Mcrx);

                if(rate>ratema) ratema=rate;
                if(rate>marate[soffset-1]) marate[soffset-1]=rate;

                rateqa=ratema;
                qarate[soffset-1]=marate[soffset-1];

                sprintf(str," %8.3f(%8.2f)",Qcrx,si*Qcrx);
                strcat(strQ[1][0],str);
                sprintf(str," %8.3f(%8.2f)",Qcrx,si*Qcrx);
                strcat(strQ[1][1],str);
                sprintf(str," %8.3f(%8.2f)",Mcrx,si*Mcrx);
                strcat(strM[1][0],str);
                sprintf(str," %8.3f(%8.2f)",Mcrx,si*Mcrx);
                strcat(strM[1][1],str);

                if(Qcrx!=0.0) sprintf(str,"           %8.3f",estress.Q[0]/Qcrx);
                else          sprintf(str,"   10.000");
                strcat(strQ[2][0],str);
                if(Qcry!=0.0) sprintf(str,"           %8.3f",estress.Q[1]/Qcrx);
                else          sprintf(str,"   10.000");
                strcat(strQ[2][1],str);
                if(Mcrx!=0.0) sprintf(str,"           %8.3f",estress.M[0]/Mcrx);
                else          sprintf(str,"   10.000");
                strcat(strM[2][0],str);
                if(Mcry!=0.0) sprintf(str,"           %8.3f",estress.M[1]/Mcrx);
                else          sprintf(str,"   10.000");
                strcat(strM[2][1],str);

                /*OUTPUT*/
                fprintf(fout0,"%s%s%s%s\n",
                        strQ[0][0],strQ[0][1],strM[0][0],strM[0][1]);

                fprintf(fout0,"     ãñóeíl:%s%s",
                        strQ[1][0],strQ[1][1]);
                fprintf(fout0,"                   %s%s\n",
                        strM[1][0],strM[1][1]);
                fprintf(fout0,"     à¿ëSó¶:%s%s",
                        strQ[2][0],strQ[2][1]);
                fprintf(fout0,"                   %s%s\n",
                        strM[2][0],strM[2][1]);
              }
              else
              {
                if(elem.sect->etype==COLUMN) na=1;
                else                         na=0;
                for(ia=0;ia<=na;ia++) /*FOR AXIS*/
                {
                  if(pstress[0]->M[ia]==0.0) facei=1.0;       /*RATE.*/
                  else if(pstress[0]->Q[1-ia]==0.0) facei=1.0;
                  else
                  {
                    facei=(fabs(pstress[0]->M[ia]/pstress[0]->Q[1-ia])
                           -elem.sect->face[ia][HEAD]/100.0)
                          /fabs(pstress[0]->M[ia]/pstress[0]->Q[1-ia]);
                  }
                  if(pstress[1]->M[ia]==0.0) facej=1.0;
                  else if(pstress[1]->Q[1-ia]==0.0) facei=1.0;
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
                    sprintf(str," %8.3f(%8.2f)",sN/1000.0,si*sN/1000.0);
                    strcat(strQ[0][0],str);
                  }
                  if(sect1.etype==GIRDER || sect1.etype==BEAM)
                  {
                    strcat(strQ[0][0],"                                      ");
                  }
                  sprintf(str," %8.3f(%8.2f) %8.3f(%8.2f)",
                             sQ[HEAD][1-ia]/1000.0,
                          si*sQ[HEAD][1-ia]/1000.0,
                             sQ[TAIL][1-ia]/1000.0,
                          si*sQ[TAIL][1-ia]/1000.0);
                  strcat(strQ[0][1-ia],str);

                  if(ia==0) sprintf(strM[0][0]," %8.3f(%8.2f)",
                                    sMt/100000.0,si*sMt/100000.0);
                  sprintf(str," %8.3f(%8.2f) %8.3f(%8.2f)",
                             sM[HEAD][ia]/100000.0,
                          si*sM[HEAD][ia]/100000.0,
                             sM[TAIL][ia]/100000.0,
                          si*sM[TAIL][ia]/100000.0);
                  strcat(strM[0][ia],str);

                  Ma=0.0;

                  if(sect1.stype==STYPE_S)                        /*S*/
                  {
                    Ma=allowablebendingofsteel(PSHORT,sect1,gmaterial.sE,
                                               h,ia,(-sN));
                  }
                  else if(sect1.stype==STYPE_RC)                 /*RC*/
                  {
                    Ma=allowablebendingofrc(elem,gmaterial,ia,(-sN),PSHORT);
                  }
                  else if(sect1.stype==STYPE_SRC)               /*SRC*/
                  {
                    if(sect1.sform.type==STEEL_RECTS)
                    {
                      Ma=allowablebendingofsrc(elem,gmaterial,ia,(-sN),PSHORT);
                    }
                  }
                  else if(sect1.stype==STYPE_WOOD)             /*WOOD*/
                  {
                    Ma=allowablebendingofwood(PSHORT,sect1,h,h,ia,(-sN));
                  }
                  else if(sect1.stype==STYPE_ALUMI)           /*ALUMI*/
                  {
					Ma=allowablebendingofalumi(PSHORT,sect1,gmaterial.alE,
											   h,ia,(-sN));            //Default
					/*
                    Mahead=allowablebendingofalumi(PSHORT,sect1,gmaterial.alE,
											   h,ia,(-sN));
					Matail=Mahead;
                    */
				  }
                  /*if(Ma<0.1) Ma=0.1;*/

                  sprintf(str," %8.3f(%8.2f) %8.3f(%8.2f)",
                          Ma/100000.0,si*Ma/100000.0,
                          Ma/100000.0,si*Ma/100000.0);
                  strcat(strM[1][ia],str);

                  if(Ma>0.0)
                  {
                    rate1=fabs(sM[HEAD][ia]/Ma);
                    rate2=fabs(sM[TAIL][ia]/Ma);
                  }
                  else
                  {
                    rate1=10.0;
                    rate2=10.0;
                  }

                  sprintf(str,"           %8.3f           %8.3f",
                          rate1,rate2);
                  strcat(strM[2][ia],str);
                  /*
                  fprintf(fout0,"Ma%s=%.3f[tfm]",
                          straxis[ia],Ma/100000.0);
                  fprintf(fout0," Mi/Ma=%.3f Mj/Ma=%.3f\n",
                          fabs(sM[HEAD][ia]/Ma),
                          fabs(sM[TAIL][ia]/Ma));
                  */
                  if(rate1>ratema) ratema=rate1;
                  if(rate1>marate[soffset-1]) marate[soffset-1]=rate1;
                  if(rate2>ratema) ratema=rate2;
                  if(rate2>marate[soffset-1]) marate[soffset-1]=rate2;

                  for(in=HEAD;in<=TAIL;in++) /*FOR END*/
                  {
                    Qa=0.0;

                    if(sect1.stype==STYPE_S)                      /*S*/
                    {
                      Qa=allowultimshearofsteel(PSHORT,gmaterial.sF,
                                                sect1,ia);
                    }
                    else if(sect1.stype==STYPE_RC)               /*RC*/
                    {
                      Qa=allowableshearofrcshort(elem,ia,
                                                 sQ[in][1-ia],
                                                 sM[in][ia]);
                    }
                    else if(sect1.stype==STYPE_SRC)             /*SRC*/
                    {
                      if(sect1.sform.type==STEEL_RECTS && ia==SX)
                      {
                        Qa=allowableshearofsrcshort(elem,ia,HSTRONG,
                                                    sQ[in][1-ia],
                                                    sM[in][ia]);
                      }
                      if(sect1.sform.type==STEEL_RECTS && ia==SY)
                      {
                        Qa=allowableshearofsrcshort(elem,ia,HWEAK,
                                                    sQ[in][1-ia],
                                                    sM[in][ia]);
                      }
                    }
                    else if(sect1.stype==STYPE_WOOD)           /*WOOD*/
                    {
                      Qa=allowableshearofwood(PSHORT,sect1);
                    }
                    else if(sect1.stype==STYPE_ALUMI)          /*ALUMI*/
                    {
                      if(sect1.sform.type==STEEL_RECTS) gmaterial.alF=sect1.srect[0].F;
                      else                              gmaterial.alF=sect1.sform.F;
                      Qa=allowultimshearofalumi(PSHORT,gmaterial.alF,
                                                sect1,ia);
                    }
                    sprintf(str," %8.3f(%8.2f)",
                            Qa/1000.0,si*Qa/1000.0);
                    strcat(strQ[1][1-ia],str);

                    /*fprintf(fout0,"Qa%s=%8.3f[tf]",
                            straxis[1-ia],Qa/1000.0);*/
                    /*if(Qa<0.1) Qa=0.1;*/

                    if(Qa>0.0) rate=fabs(sQ[in][1-ia]/Qa);
                    else       rate=10.0;

                    if(rate>rateqa) rateqa=rate;
                    if(rate>qarate[soffset-1])
                    {
                      qarate[soffset-1]=rate;
                    }
                    sprintf(str,"           %8.3f",rate);
                    strcat(strQ[2][1-ia],str);

                    /*
                    fprintf(fout0," Q%s/Qa=%.3f\n",strend[in],rate);
                    */
                  }
                }

                if(sect1.etype==GIRDER || sect1.etype==BEAM)
                {
                  sprintf(strQ[1][0],"                                      ");
                  sprintf(strQ[2][0],"                                      ");
                  sprintf(strM[1][1],"                              ");
                  sprintf(strM[2][1],"                              ");
                }
                fprintf(fout0,"%s%s%s%s\n",
                        strQ[0][0],strQ[0][1],strM[0][0],strM[0][1]);
                fprintf(fout0,"     ãñóeíl:                   %s%s",
                        strQ[1][0],strQ[1][1]);
                fprintf(fout0,"                   %s%s\n",
                        strM[1][0],strM[1][1]);
                fprintf(fout0,"     à¿ëSó¶:         %s%s",
                        strQ[2][0],strQ[2][1]);
                fprintf(fout0,"                   %s%s\n",
                        strM[2][0],strM[2][1]);
              }
            }
            fprintf(fout0,"\n");
            /*if(check){gets(non); if(non[0]!='\0') return 0;}*/
          }
/*MessageBox(NULL,"Pass 1","SRCan",MB_OK);*/

          /*ULTIMATE...............................................*/
          /*fprintf(stderr,"ELEM:%d SECT:%d\n",ielem,isect);*/

if(calc[PULTIMATE]==1 || sect1.stype==STYPE_PC){
/*if(!strcmp(prj,"aki") && isect==121){*/ /*SURFACE TEST*/
          /*MATERIAL FOR ULTIMATE*/
          gmaterial.sftu=2400.0*jis;                  /*STEEL SN400*/
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

              sprintf(strQ[0][0],"èIã«%s%sï˚å¸:",
                      strhoko[hoko],strhugo[sign]);

              if(elem.sect->etype==COLUMN) na=1;
              else                         na=0;
              for(ia=0;ia<=na;ia++) /*FOR AXIS*/
              {
                if(pstress[0]->M[ia]==0.0) facei=1.0;       /*RATE.*/
                else if(pstress[0]->Q[1-ia]==0.0) facei=1.0;
                else
                {
                  facei=(fabs(pstress[0]->M[ia]/pstress[0]->Q[1-ia])
                         -elem.sect->face[ia][HEAD]/100.0)
                        /fabs(pstress[0]->M[ia]/pstress[0]->Q[1-ia]);
                }
                if(pstress[1]->M[ia]==0.0) facej=1.0;
                else if(pstress[0]->Q[1-ia]==0.0) facej=1.0;
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
                  sprintf(str," %8.3f(%8.2f)",uN/1000.0,si*uN/1000.0);
                  strcat(strQ[0][0],str);
                }
                if(sect1.etype==GIRDER || sect1.etype==BEAM)
                {
                  strcat(strQ[0][0],"                                      ");
                }
                sprintf(str," %8.3f(%8.2f) %8.3f(%8.2f)",
                        uQ[HEAD][1-ia]/1000.0,si*uQ[HEAD][1-ia]/1000.0,
                        uQ[TAIL][1-ia]/1000.0,si*uQ[TAIL][1-ia]/1000.0);
                strcat(strQ[0][1-ia],str);

                if(ia==0) sprintf(strM[0][0]," %8.3f(%8.2f)",uMt/100000.0,si*uMt/100000.0);
                sprintf(str," %8.3f(%8.2f) %8.3f(%8.2f)",
                        uM[HEAD][ia]/100000.0,si*uM[HEAD][ia]/100000.0,
                        uM[TAIL][ia]/100000.0,si*uM[TAIL][ia]/100000.0);
                strcat(strM[0][ia],str);

                Mu=0.0;
                Qp[HEAD]=0.0;
                Qp[TAIL]=0.0;

                if(sect1.stype==STYPE_S)                        /*S*/
                {
                  Mu=ultimatebendingofsteel(sect1,gmaterial.sE,
                                                  gmaterial.sF,
                                                  h,h,ia,(-uN));
                }
                else if(sect1.stype==STYPE_RC)                 /*RC*/
                {
/*HAKODATE FG.*/
/*if(0)
{
  if(isect==302) uN=0.0;
}*/

/*SURFACE TEST FOR AKIHABARA,HIGASHIIKEBUKURO,KOMAKI.*/
/*for(ii=0;ii<22;ii++){
  if(ii<=10)       uN=1944042.0*(10.0-ii)/10.0;
  else if(ii<=20)  uN=-638442.0*(ii-10.0)/10.0;
  else if(ii==21)  uN=0.5*(1944042.0-638442.0);
*/
                  Mu=ultimatebendingofsrc(elem,ia,(-uN),
                                          &Ns,&Ms,&Nr,&Mr,&Nc,&Mc);
/*
  fprintf(fout0,"Surface %d %15.3f %15.3f\n",ii,uN,Mu);
}*/
                }
                else if(sect1.stype==STYPE_SRC)               /*SRC*/
                {
                  Mu=ultimatebendingofsrc(elem,ia,(-uN),
                                          &Ns,&Ms,&Nr,&Mr,&Nc,&Mc);
                }
                else if(sect1.stype==STYPE_PC)                 /*PC*/
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

                if(sect1.stype==STYPE_PC && sect1.etype==GIRDER)
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

                  sprintf(str," %8.3f(%8.2f) %8.3f(%8.2f)",
                          Muhead/100000.0,si*Muhead/100000.0,
                          Mutail/100000.0,si*Mutail/100000.0);
                  strcat(strM[1][ia],str);
                  sprintf(str,"           %8.3f           %8.3f",
                          fabs(uM[HEAD][ia]/Muhead),fabs(uM[TAIL][ia]/Mutail));
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

                  sprintf(str," %8.3f(%8.2f) %8.3f(%8.2f)",
                          Mu/100000.0,si*Mu/100000.0,
                          Mu/100000.0,si*Mu/100000.0);
                  strcat(strM[1][ia],str);
                  sprintf(str,"           %8.3f           %8.3f",
                          fabs(uM[HEAD][ia]/Mu),fabs(uM[TAIL][ia]/Mu));
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

                  if(sect1.stype==STYPE_S)                      /*S*/
                  {
                    Qu=allowultimshearofsteel(PULTIMATE,gmaterial.sF,
                                              sect1,ia);
                  }
                  else if(sect1.stype==STYPE_RC)               /*RC*/
                  {
                    Qu=ultimateshearofrc(elem,gmaterial,ia,
                                         (-uN),
                                         uQ[in][1-ia],
                                         uM[in][ia]);
                  }
                  else if(sect1.stype==STYPE_SRC)             /*SRC*/
                  {
                    if(sect1.sform.type==STEEL_RECTS && ia==SX)
                    {
                      Qu=ultimateshearofsrc(elem,ia,HSTRONG,
                                            uQ[in][1-ia],
                                            uM[in][ia]);
                      /*APPROXIMATION.Q,M MUST BE OF RC.*/
                    }
                    if(sect1.sform.type==STEEL_RECTS && ia==SY)
                    {
                      Qu=ultimateshearofsrc(elem,ia,HWEAK,
                                            uQ[in][1-ia],
                                            uM[in][ia]);
                      /*APPROXIMATION.Q,M MUST BE OF RC.*/
                    }
                  }
                  else if(sect1.stype==STYPE_PC)               /*PC*/
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

                  sprintf(str," %8.3f(%8.2f)",Qu/1000.0,si*Qu/1000.0);
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
                  sprintf(str,"           %8.3f",rate);
                  strcat(strQ[2][1-ia],str);

                  /*
                  fprintf(fout0," Q%s/Qu=%.3f\n",strend[in],rate);
                  */

                  if(sect1.stype==STYPE_RC)
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

                  if(sect1.stype==STYPE_PC || sect1.stype==STYPE_RC)
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

                  sprintf(str," %8.3f(%8.2f)",Qp[in]/1000.0,si*Qp[in]/1000.0);
                  strcat(strQ[3][1-ia],str);

                  rate=fabs(Qp[in]/Qu);
                  if(rate>ratequ) ratequ=rate;
                  if(rate>qurate[soffset-1]) qurate[soffset-1]=rate;
                  sprintf(str,"           %8.3f",rate);
                  strcat(strQ[4][1-ia],str);
                }
              }

              if(sect1.etype==GIRDER || sect1.etype==BEAM)
              {
                sprintf(strQ[1][0],"                                      ");
                sprintf(strQ[2][0],"                                      ");
                sprintf(strQ[3][0],"                                      ");
                sprintf(strQ[4][0],"                                      ");
                sprintf(strM[1][1],"                  ");
                sprintf(strM[2][1],"                  ");
              }
              fprintf(fout0,"%s%s%s%s\n",
                      strQ[0][0],strQ[0][1],strM[0][0],strM[0][1]);
              fprintf(fout0,"      èIã«íl:                  %s%s",
                      strQ[1][0],strQ[1][1]);
              fprintf(fout0,"                   %s%s\n",
                      strM[1][0],strM[1][1]);
              fprintf(fout0,"      à¿ëSó¶:        %s%s",
                      strQ[2][0],strQ[2][1]);
              fprintf(fout0,"                   %s%s\n",
                      strM[2][0],strM[2][1]);
              fprintf(fout0,"  ã@ç\å`ê¨éû:                  %s%s\n",
                      strQ[3][0],strQ[3][1]);
              fprintf(fout0,"      à¿ëSó¶:        %s%s\n",
                      strQ[4][0],strQ[4][1]);
            }
            fprintf(fout0,"\n");
            /*if(check){gets(non); if(non[0]!='\0') return 0;}*/
          }
}/*CALCULATE ULTIMATE*/

          fprintf(fout0,"MAX:Q/Qa=%.5f",rateqa);
          if(calc[PULTIMATE]==1 || sect1.stype==STYPE_PC)
          {
            fprintf(fout0," Q/Qu=%.5f",ratequ);
          }
          fprintf(fout0," M/Ma=%.5f",ratema);
          if(calc[PULTIMATE]==1 || sect1.stype==STYPE_PC)
          {
            fprintf(fout0," M/Mu=%.5f",ratemu);
          }
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
/*
            fprintf(frate,"ELEM: %5d SECT: %5d  %.5f %.5f %.5f %.5f\n",
                    ielem,isect,rateqa,ratequ,ratema,ratemu);
*/
			fprintf(frate,"ELEM: %5d SECT: %5d  %.5f %.5f %.5f %.5f OCODE: %5d\n",
					ielem,isect,rateqa,ratequ,ratema,ratemu,sect1.ocode); /****SRCANMAX****/
          }

          if(ratema==10.0) dabocount++;

          /*if(check){gets(non); if(non[0]!='\0') return 0;}*/
        }

        if(soffset && sect1.etype==WALL)/*...........SHEAR OF WALL.*/
        {
          Qa=0.0;
          Qu=0.0;

          /*fprintf(stderr,"ELEM:%d SECT:%d\n",ielem,isect);*/

          fprintf(fout0,"-------------------------------------");
          fprintf(fout0,"-------------------------------------");
          fprintf(fout0,"-------------------------------------");
          fprintf(fout0,"-------------------------------------");
          fprintf(fout0,"-------------------------------------");
          fprintf(fout0,"-----------------\n");
          fprintf(fout0,"ïîçﬁ:%4d ",ielem);
          fprintf(fout0,"éní[:%3d èIí[:%3d ",ni,nj);

		  /* fprintf(fout0,"ífñ :%3d",isect); */

		  /* 150429 fukushima for original section */
		  if(sect1.code==sect1.ocode)
		  {
			fprintf(fout0,"ífñ :Å¸%3d/Å@%3d",sect1.ocode,isect);
		  }
		  else
		  {
			 fprintf(fout0,"ífñ :Å@%3d/Å¸%3d",sect1.ocode,isect);
		  }
		  /***/
          if(sect1.stype==STYPE_RC)        fprintf(fout0,"=ÇqÇbï« ");
          else if(sect1.stype==STYPE_WOOD) fprintf(fout0,"=ñÿï« ");
          else if(sect1.stype==STYPE_GLASS) fprintf(fout0,"=ÉKÉâÉXï« ");
          else if(sect1.stype==STYPE_ACRYL) fprintf(fout0,"=ÉAÉNÉäÉãï« ");
          else                             fprintf(fout0,"=ïsñæï« ");

          fprintf(fout0,"ï«å˙:%.0f[mm] ",sect1.thick*10.0);

          h=100.0*wallheight(elem); /*[cm]*/
          h0[1]=h-elem.sect->face[1][HEAD]-elem.sect->face[1][TAIL];
          if(h0[1]<0.0) h0[1]=0.0;
          l=100.0*walllength(elem); /*[cm]*/
          l0=l-elem.sect->face[0][HEAD]-elem.sect->face[0][TAIL];
          if(l0<0.0) l0=0.0;
          fprintf(fout0,"ì‡ñ@:L=%.1f[cm] H=%.1f[cm] ",l0,h0[1]);

          if(l0<=0.0 || h0[1]<=0.0)
          {
            wrate1=0.0;
            wrate2=0.0;
          }
          else
          {
            wrate1=(elem.sect->wlength)/l0; /*RATE OF WINDOW.*/
            wrate2=sqrt((elem.sect->wlength)*(elem.sect->wheight)
                        /(l0*h0[1]));
          }
          if(wrate1>=wrate2) elem.sect->windowrate=wrate1;
          else               elem.sect->windowrate=wrate2;
          fprintf(fout0,"äJå˚ó¶:1-r=%.3f\n",elem.sect->windowrate);

          l0/=2.0; /*1 OF 2 BRACES.*/

          /*LONG*/
          if(sect1.stype==STYPE_RC)
          {
            Qa=allowultimshearofrcwall(elem,l0,PLONG);
          }
          if(sect1.stype==STYPE_WOOD)
          {
            fprintf(fout0,"í∑ä˙   :");
            Qa=allowableshearofwoodwall(elem,l0,PLONG);
          }

          lN=elem.head.z.N;
          lQ=lN*l/sqrt(l*l+h*h)*1000.0;

          fprintf(fout0,"í∑ä˙   :");
          fprintf(fout0,"N=%7.3f[tf](%7.2f[kN])",lN,si*lN);
          fprintf(fout0," êÖïΩê¨ï™:Qh=%7.3f[tf](%7.2f[kN]) ",
                  lQ/1000.0,si*lQ/1000.0);

if(Qa==0.0)
{
  sprintf(str,"Qa=0 ELEM=%d SECT=%d",elem.code,sect1.code);
  MessageBox(NULL,str,"WALL",MB_OK);
}
          rate=fabs(lQ/Qa);
          if(rate>rateqa) rateqa=rate;
          if(rate>qarate[soffset-1]) qarate[soffset-1]=rate;
          fprintf(fout0," Qs/Qa=%7.5f\n",rate);

          /*MATERIAL SHORT*/
          /*gmaterial.Fc=240.0;*/                  /*CONCRETE Fc240*/
          /*gmaterial.wft=3000.0;*/
          /*gmaterial.wfp=3300.0;*/
          /*gmaterial.cfs=11.1;*/

          if(sect1.stype==STYPE_RC)
          {
            Qa=allowultimshearofrcwall(elem,l0,PSHORT);
            Qu=allowultimshearofrcwall(elem,l0,PULTIMATE);
          }
          if(sect1.stype==STYPE_WOOD)
          {
            fprintf(fout0,"íZä˙   :");
            Qa=allowableshearofwoodwall(elem,l0,PSHORT);
          }

if(!strcmp(prj,"aki")) wafact=1.0;
if(!strcmp(prj,"tohu"))
{
  wafact=1.0;
  wufact=1.5;
}
if(!strcmp(prj,"haya"))
{
  wafact=1.0;
  wufact=1.5;
}
if(!strcmp(prj,"kiyo")) wafact=1.0;
if(!strcmp(prj,"kamime")) wafact=1.0;
if(!strcmp(prj,"maebasi")) wafact=1.0;

          for(hoko=0;hoko<=1;hoko++)
          {
            if(hoko==0) eN=elem.head.x.N;
            if(hoko==1) eN=elem.head.y.N;

            sQ[0][0]=eN*l/sqrt(l*l+h*h)*1000.0;
            fprintf(fout0,"êÖïΩéû%s:",strhoko[hoko]);
            fprintf(fout0,"N=%7.3f[tf](%7.2f[kN])",eN,si*eN);
            fprintf(fout0," êÖïΩê¨ï™:Qh=%7.3f[tf](%7.2f[kN]) ",
                    sQ[0][0]/1000.0,si*sQ[0][0]/1000.0);

            fprintf(fout0,"íZä˙=í∑ä˙Å{%.1fÅ~êÖïΩ:",wafact);
            fprintf(fout0,"Qs=%7.3f[tf](%7.2f[kN])",
                       (fabs(lQ)+fabs(wafact*sQ[0][0]))/1000.0,
                    si*(fabs(lQ)+fabs(wafact*sQ[0][0]))/1000.0);

if(calc[PULTIMATE]==1 || !strcmp(prj,"tohu"))
{
            fprintf(fout0," èIã«=í∑ä˙Å{%.1fÅ~êÖïΩ:",wufact);
            fprintf(fout0,"Qu=%7.3f[tf](%7.2f[kN])",
                       (fabs(lQ)+fabs(wufact*sQ[0][0]))/1000.0,
                    si*(fabs(lQ)+fabs(wufact*sQ[0][0]))/1000.0);
}

            rate=fabs((fabs(lQ)+fabs(wafact*sQ[0][0]))/Qa);
            if(rate>rateqa) rateqa=rate;
            if(rate>qarate[soffset-1]) qarate[soffset-1]=rate;
            fprintf(fout0," Qs/Qa=%7.5f",rate);

if(calc[PULTIMATE]==1 || !strcmp(prj,"tohu"))
{
            rate=fabs((fabs(lQ)+fabs(wufact*sQ[0][0]))/Qu);
            if(rate>ratequ) ratequ=rate;
            if(rate>qurate[soffset-1]) qurate[soffset-1]=rate;
            fprintf(fout0," Qs/Qu=%7.5f",rate);
}

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
/*
            fprintf(frate,"ELEM: %5d SECT: %5d  %.5f %.5f\n",
                    ielem,isect,rateqa,ratequ);
*/
			fprintf(frate,"ELEM: %5d SECT: %5d  %.5f %.5f 0.0 0.0 OCODE: %5d\n",
					ielem,isect,rateqa,ratequ,sect1.ocode); /****SRCANMAX****/
          }
        }

        if(soffset && sect1.etype==BRACE)            /*STEEL BRACE.*/
        {
          /*Qa=0.0;*/

          fprintf(stderr,"ELEM:%d SECT:%d\n",ielem,isect);

          fprintf(fout0,"-------------------------------------");
          fprintf(fout0,"-------------------------------------");
          fprintf(fout0,"-------------------------------------");
          fprintf(fout0,"-------------------------------------");
          fprintf(fout0,"-------------------------------------");
          fprintf(fout0,"-----------------\n");
          fprintf(fout0,"ïîçﬁ:%4d ",ielem);
          fprintf(fout0,"éní[:%3d èIí[:%3d ",ni,nj);

          fprintf(fout0,"ífñ :%3d",isect);
          fprintf(fout0,"=ìSçúãÿà· ");

          fprintf(fout0,"ífñ êœ:%.3f[cm2]\n",sect1.thick);

          /*MATERIAL LONG*/
          gmaterial.sft=sect1.cF*jis/1.5; /*STEEL F*/

          Qa=sect1.thick*gmaterial.sft/1000.0; /*[tf]*/
          fprintf(fout0,"í∑ä˙ãñóe:Nal=%7.3f[tf](%7.2f[kN])\n",Qa,si*Qa);

          lN=elem.head.z.N;

          fprintf(fout0,"âîíºéû  :");
          fprintf(fout0,"Nl =%7.3f[tf](%7.2f[kN]) ",lN,si*lN);

          rate=fabs(lN/Qa);
          if(rate>rateqa) rateqa=rate;
          if(rate>qarate[soffset-1]) qarate[soffset-1]=rate;
          fprintf(fout0," Nl/Na=%7.5f",rate);
          fprintf(fout0,"\n");

          /*MATERIAL SHORT*/
          gmaterial.sft=sect1.cF*jis; /*STEEL F*/

          Qa=sect1.thick*gmaterial.sft/1000.0; /*[tf]*/
          if(!strcmp(prj,"kirigami"))      /*kirigami  : SKIP SHORT.*/
          {}
          else
          {
            fprintf(fout0,"íZä˙ãñóe:Nas=%7.3f[tf](%7.2f[kN])\n",Qa,si*Qa);

            for(hoko=0;hoko<=1;hoko++)
            {
              if(hoko==0) eN=elem.head.x.N;
              if(hoko==1) eN=elem.head.y.N;

              fprintf(fout0,"êÖïΩéû%s :",strhoko[hoko]);
              fprintf(fout0,"Ns =%7.3f[tf](%7.2f[kN]) ",eN,si*eN);

              /*fprintf(fout0,"\n");*/

              rate=fabs(eN/Qa);
              if(rate>rateqa) rateqa=rate;
              if(rate>qarate[soffset-1]) qarate[soffset-1]=rate;
              fprintf(fout0," Ns/Na=%7.5f",rate);

              fprintf(fout0,"\n");
            }
          }
          if(rateqa>1.0 || ratequ>1.0)
          {
            fprintf(fsafe,"NG ELEM:%4d SECT:%4d ",ielem,isect);
            fprintf(fsafe,"MAX:N/Na=%.5f N/Nu=%.5f\n",rateqa,ratequ);
          }
          /*if(rateqa<=1.0 && ratequ<=1.0)
          {
            fprintf(fsafe,"ELEM:%4d SECT:%4d ",ielem,isect);
            fprintf(fsafe,"MAX:Q/Qa=%.5f Q/Qu=%.5f\n",rateqa,ratequ);
          }*/

          if(frate!=NULL)
          {
/*
            fprintf(frate,"ELEM: %5d SECT: %5d  %.5f %.5f\n",
                    ielem,isect,rateqa,ratequ);
*/			fprintf(frate,"ELEM: %5d SECT: %5d  %.5f %.5f 0.0 0.0 OCODE: %5d\n",
					ielem,isect,rateqa,ratequ,sect1.ocode); /****SRCANMAX****/
          }
        }

        if(soffset && sect1.etype==SLAB)/*...........SHEAR OF SLAB.*/
        {
          if(slabcount==0)      slabcount=1;
          else if(slabcount==1) slabcount=0;

          Qa=0.0;
          Qu=0.0;

          /*fprintf(stderr,"ELEM:%d SECT:%d\n",ielem,isect);*/

          fprintf(fout0,"-------------------------------------");
          fprintf(fout0,"-------------------------------------");
          fprintf(fout0,"-------------------------------------");
          fprintf(fout0,"-------------------------------------");
          fprintf(fout0,"-------------------------------------");
          fprintf(fout0,"-----------------\n");
          fprintf(fout0,"ïîçﬁ:%4d ",ielem);
          fprintf(fout0,"éní[:%3d èIí[:%3d ",ni,nj);

		  /* fprintf(fout0,"ífñ :%3d",isect); */

		  /* 150429 fukushima for original section */
		  if(sect1.code==sect1.ocode)
		  {
			fprintf(fout0,"ífñ :Å¸%3d/Å@%3d",sect1.ocode,isect);
		  }
		  else
		  {
			fprintf(fout0,"ífñ :Å@%3d/Å¸%3d",sect1.ocode,isect);
		  }
          /***/
          if(sect1.stype==STYPE_RC)        fprintf(fout0,"=ÇqÇbè∞ ");
          else if(sect1.stype==STYPE_WOOD) fprintf(fout0,"=ñÿè∞ ");
          else                             fprintf(fout0,"=ïsñæè∞ ");

          fprintf(fout0,"è∞å˙:%.0f[mm] ",sect1.thick*10.0);

          h=100.0*slablength(elem,pair); /*[cm]*/
          h0[1]=h-elem.sect->face[1][HEAD]-elem.sect->face[1][TAIL];
          if(h0[1]<0.0) h0[1]=0.0;
          l=100.0*slabwidth(elem,pair); /*[cm]*/
          l0=l-elem.sect->face[0][HEAD]-elem.sect->face[0][TAIL];
          if(l0<0.0) l0=0.0;
          fprintf(fout0,"ì‡ñ@:L=%.1f[cm] H=%.1f[cm] ",l0,h0[1]);

          wrate1=(elem.sect->wlength)/l0; /*RATE OF WINDOW.*/
          wrate2=sqrt((elem.sect->wlength)*(elem.sect->wheight)
                      /(l0*h0[1]));
          if(wrate1>=wrate2) elem.sect->windowrate=wrate1;
          else               elem.sect->windowrate=wrate2;
          fprintf(fout0,"äJå˚ó¶:1-r=%.3f\n",elem.sect->windowrate);

          l0/=2.0; /*1 OF 2 BRACES.*/

          /*LONG*/
          if(sect1.stype==STYPE_RC)
          {
            Qa=allowultimshearofrcwall(elem,l0,PLONG);
          }
          if(sect1.stype==STYPE_WOOD)
          {
            fprintf(fout0,"í∑ä˙   :");
            Qa=allowableshearofwoodwall(elem,l0,PLONG);
          }

          lN=elem.head.z.N;
          lQ=lN*l/sqrt(l*l+h*h)*1000.0;

          fprintf(fout0,"í∑ä˙   :");
          fprintf(fout0,"N=%7.3f[tf](%7.2f[kN])",lN,si*lN);
          fprintf(fout0," êÖïΩê¨ï™:Qh=%7.3f[tf](%7.2f[kN]) ",
                  lQ/1000.0,si*lQ/1000.0);

          rate=fabs(lQ/Qa);
          if(rate>rateqa) rateqa=rate;
          if(rate>qarate[soffset-1]) qarate[soffset-1]=rate;
          fprintf(fout0," Qs/Qa=%7.5f\n",rate);

          /*MATERIAL SHORT*/
          /*gmaterial.Fc=240.0;*/                  /*CONCRETE Fc240*/
          /*gmaterial.wft=3000.0;*/
          /*gmaterial.wfp=3300.0;*/
          /*gmaterial.cfs=11.1;*/

          if(sect1.stype==STYPE_RC)
          {
            Qa=allowultimshearofrcwall(elem,l0,PSHORT);
            /*Qu=allowultimshearofrcwall(elem,l0,PULTIMATE);*/
          }
          if(sect1.stype==STYPE_WOOD)
          {
            fprintf(fout0,"íZä˙   :");
            Qa=allowableshearofwoodwall(elem,l0,PSHORT);
          }

if(!strcmp(prj,"aki")) wafact=1.0;

          for(hoko=0;hoko<=1;hoko++)
          {
            if(hoko==0) eN=elem.head.x.N;
            if(hoko==1) eN=elem.head.y.N;

            sQ[0][0]=eN*l/sqrt(l*l+h*h)*1000.0;
            fprintf(fout0,"êÖïΩéû%s:",strhoko[hoko]);
            fprintf(fout0,"N=%7.3f[tf](%7.2f[kN])",eN,si*eN);
            fprintf(fout0," êÖïΩê¨ï™:Qh=%7.3f[tf](%7.2f[kN]) ",
                    sQ[0][0]/1000.0,si*sQ[0][0]/1000.0);

            fprintf(fout0,"íZä˙=í∑ä˙Å{%.1fÅ~êÖïΩ:",wafact);
            fprintf(fout0,"Qs=%7.3f[tf](%7.2f[kN])",
                       (lQ+wafact*sQ[0][0])/1000.0,
                    si*(lQ+wafact*sQ[0][0])/1000.0);

if(calc[PULTIMATE]==1)
{
            fprintf(fout0," èIã«=í∑ä˙Å{%.1fÅ~êÖïΩ:",wufact);
            fprintf(fout0,"Qu=%7.3f[tf](%7.2f[kN])",
                       (lQ+wufact*sQ[0][0])/1000.0,
                    si*(lQ+wufact*sQ[0][0])/1000.0);
}

            rate=fabs((lQ+wafact*sQ[0][0])/Qa);
            if(rate>rateqa) rateqa=rate;
            if(rate>qarate[soffset-1]) qarate[soffset-1]=rate;
            fprintf(fout0," Qs/Qa=%7.5f",rate);

if(calc[PULTIMATE]==1)
{
            rate=fabs((lQ+wufact*sQ[0][0])/Qu);
            if(rate>ratequ) ratequ=rate;
            if(rate>qurate[soffset-1]) qurate[soffset-1]=rate;
            fprintf(fout0," Qs/Qu=%7.5f",rate);
}

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
/*
            fprintf(frate,"ELEM: %5d SECT: %5d  %.5f %.5f\n",
                    ielem,isect,rateqa,ratequ);
*/
			fprintf(frate,"ELEM: %5d SECT: %5d  %.5f %.5f 0.0 0.0 OCODE: %5d\n",
					ielem,isect,rateqa,ratequ,sect1.ocode); /****SRCANMAX****/
          }
        }
      }
    }
  }
  fprintf(fout0,"=====================================");
  fprintf(fout0,"=====================================");
  fprintf(fout0,"=====================================");
  fprintf(fout0,"=====================================");
  fprintf(fout0,"=====================================");
  fprintf(fout0,"=================\n");
  fprintf(fout0,"äeífñ éÌï ÇÃãñóe,èIã«ã»Ç∞à¿ëSó¶ÇÃç≈ëÂíl\n\n");

  rateqa=0.0; ratequ=0.0;
  ratema=0.0; ratemu=0.0;

  for(is=1;is<=nsect;is++)
  {
    soffset=getsectionform(flist,codelist[is-1],&sect1);
    if(soffset)
    {
      fprintf(fout0,"ífñ ãLçÜ:%4d",sect1.code);

      if(sect1.stype     ==STYPE_S)     fprintf(fout0," Çr  ");
      else if(sect1.stype==STYPE_RC)    fprintf(fout0," ÇqÇb");
      else if(sect1.stype==STYPE_SRC)   fprintf(fout0," ÇrÇqÇb");
      else if(sect1.stype==STYPE_PC)    fprintf(fout0," ÇoÇb");
      else if(sect1.stype==STYPE_WOOD)  fprintf(fout0," ñÿ  ");
      else if(sect1.stype==STYPE_GLASS) fprintf(fout0," ÉKÉâÉX");
      else if(sect1.stype==STYPE_ACRYL) fprintf(fout0," ÉAÉNÉäÉã");
      else if(sect1.stype==STYPE_ALUMI) fprintf(fout0," ÉAÉãÉ~");
      else                              fprintf(fout0,"     ");
      if(sect1.etype==COLUMN)      fprintf(fout0,"íå  ");
      else if(sect1.etype==GIRDER) fprintf(fout0,"ëÂó¿");
      else if(sect1.etype==BEAM)   fprintf(fout0,"è¨ó¿");
      else if(sect1.etype==BRACE)  fprintf(fout0,"ãÿà·");
      else if(sect1.etype==WALL)   fprintf(fout0,"ï«  ");
      else if(sect1.etype==SLAB)   fprintf(fout0,"è∞  ");
      else                         fprintf(fout0,"ïsñæ");

      if(sect1.etype!=WALL &&
         sect1.etype!=SLAB &&
         sect1.etype!=BRACE)
      {
        translatesection(sect1,gmaterial,&As,&Ac,&Ar,&Yg,SX);
        fprintf(fout0," As=%7.2f[cm2]",As);
        fprintf(fout0," Ar=%7.2f[cm2]",Ar);
        fprintf(fout0," Ac=%8.2f[cm2]",Ac);

        fprintf(fout0," MAX:Q/Qa=%7.5f",qarate[soffset-1]);
        if(calc[PULTIMATE]==1 || sect1.stype==STYPE_PC)
        {
          fprintf(fout0," Q/Qu=%7.5f",qurate[soffset-1]);
        }
        fprintf(fout0," M/Ma=%7.5f",marate[soffset-1]);
        if(calc[PULTIMATE]==1 || sect1.stype==STYPE_PC)
        {
          fprintf(fout0," M/Mu=%7.5f",murate[soffset-1]);
        }
        fprintf(fout0,"\n");
      }
      if(sect1.etype==BRACE)
      {
        fprintf(fout0," As=%7.2f[cm2]",sect1.thick);

        fprintf(fout0,"                                 ");
        fprintf(fout0," MAX:N/Na=%7.5f",qarate[soffset-1]);
        /*fprintf(fout0," N/Nu=%7.5f",qurate[soffset-1]);*/
        fprintf(fout0,"\n");
      }
      if(sect1.etype==WALL)
      {
        fprintf(fout0," t=%4.1f[cm]",sect1.thick);

        if(sect1.stype==STYPE_RC)
        {
          fprintf(fout0," ps=%7.5f",sect1.shearrein[0]);
        }
        if(sect1.stype==STYPE_WOOD)
        {
          fprintf(fout0,"           ");
        }

        fprintf(fout0,"                           ");
        fprintf(fout0," MAX:Q/Qa=%7.5f",qarate[soffset-1]);
        if(calc[PULTIMATE]==1 || !strcmp(prj,"tohu"))
        {
          fprintf(fout0," Q/Qu=%7.5f",qurate[soffset-1]);
        }
        fprintf(fout0,"\n");
      }
      if(sect1.etype==SLAB)
      {
        fprintf(fout0," t=%4.1f[cm]",sect1.thick);

        if(sect1.stype==STYPE_RC)
        {
          fprintf(fout0," ps=%7.5f",sect1.shearrein[0]);
        }
        if(sect1.stype==STYPE_WOOD)
        {
          fprintf(fout0,"           ");
        }

        fprintf(fout0,"                           ");
        fprintf(fout0," MAX:Q/Qa=%7.5f",qarate[soffset-1]);
        if(calc[PULTIMATE]==1 || sect1.stype==STYPE_PC)
        {
          fprintf(fout0," Q/Qu=%7.5f",qurate[soffset-1]);
        }
        fprintf(fout0,"\n");
      }
      if(qarate[soffset-1]>rateqa) rateqa=qarate[soffset-1];
      if(qurate[soffset-1]>ratequ) ratequ=qurate[soffset-1];
      if(marate[soffset-1]>ratema) ratema=marate[soffset-1];
      if(murate[soffset-1]>ratemu) ratemu=murate[soffset-1];
    }
  }

  sprintf(txt,"à¿ëSó¶ÇÃç≈ëÂíl");
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

  if(globalmessageflag==1)
  {
    sprintf(non,"TENSION DABO = %d",dabocount);
    MessageBox(NULL,non,"SRCan",MB_OK);
  }

  if(globalmessageflag==1)
  {
    sprintf(non,"\nCompleted."); strcat(txt,non);
    MessageBox(NULL,txt,"SRCan",MB_OK);
  }
  return 1;
}/*srcan001*/

int createyieldsurface(struct arclmframe *af)
{
  int ii,ie,is,ia,na;
  /*int sign;*/        /*SIGN OF HORIZONTAL LOAD 0:POSITIVE 1:NEGATIVE.*/
  double uN,Qp[2];  /*u:ULTIMATE*/
  double Nu,Nmax,Nmin,Mtu;
  double Mu,Qu; /*a:ALLOABLE u:ULTIMATE*/
  double Ns,Ms,Nr,Mr,Nc,Mc;            /*s:STEEL r:REIN c:CONCRETE.*/
  double h,l; /*WALL HEIGHT,LENGTH*/
  double h0[2],l0; /*INNER HEIGHT,LENGTH*/
  /*double facei,facej,face[2][2];*/
  struct section sect1;
  long int codelist[MAXSECT];
  double wrate1,wrate2; /*RATE OF WINDOW.*/

  FILE *flist;
  char txt[400],non[80]/*,str[256]*/;
  int nelem,nsect,nsects,soffset,loff;
  struct element elem;
  struct owire *elems;
  struct section *sects;
  char dir[]=DIRECTORY;

  struct surface ys;
  struct structnode node[2];

  /*BEGINNING MESSAGE*/
  if(MessageBox(NULL,"Begin","Create Yield Surface",MB_OKCANCEL)
     ==IDCANCEL) return 0;

  /*OPEN SECTION LIST*/
  flist=fgetstofopen(dir,"r",ID_SECTIONFILE);
  if(flist==NULL) return 0;

  fout0=NULL;

  /*INITIAL*/
  jis=1.0;                                  /*1.1=JIS STEEL.*/
  nelem=af->nelem;
  nsect=af->nsect;
  elems=af->elems;

  /*ALLOC MEMORY*/
  sects=(struct section *)malloc(nsect*sizeof(struct section));
  if(sects==NULL) return 0;

  /*SECTIONS INTO MEMORY.*/
  nsects=getcodelist(flist,codelist);
  if(nsects>MAXSECT)
  {
    MessageBox(NULL,"MAIN:SECTIONS OVERFLOW.","SRCan",MB_OK);
    return 0;
  }

  for(is=0;is<nsect;is++)
  {
    getsectionform(flist,(af->sects+is)->code,(sects+is));
  }

  for(ie=0;ie<nelem;ie++) /*CALCULATION*/
  {
    currentpivot(ie+1,nelem);

    /*COPY ELEMENT DATA*/
    loff=(elems+ie)->sect->loff;
    soffset=(sects+loff)->soff;
    sect1=*(sects+loff);

    /*sprintf(txt,"Elem=%d Sect=%d Offset=%d",
            (elems+ie)->code,sect1.code,soffset);
    MessageBox(NULL,txt,"Surface",MB_OK);*/

    node[0].code=(elems+ie)->node[0]->code;
    node[0].x   =(elems+ie)->node[0]->d[GX];
    node[0].y   =(elems+ie)->node[0]->d[GY];
    node[0].z   =(elems+ie)->node[0]->d[GZ];
    node[1].code=(elems+ie)->node[1]->code;
    node[1].x   =(elems+ie)->node[1]->d[GX];
    node[1].y   =(elems+ie)->node[1]->d[GY];
    node[1].z   =(elems+ie)->node[1]->d[GZ];

    elem.code=(elems+ie)->code;
    elem.sect=sects+loff;
    elem.node[0]=&node[0];
    elem.node[1]=&node[1];
    elem.cmqcode=0;
    elem.Mo=0.0;

    /*INITIALIZE SURFACE*/
    ys.exp=1.0;
    for(ii=0;ii<6;ii++)
    {
      ys.fmax[ii]=100.0;
      ys.fmin[ii]=100.0;
    }

    /*INITIALIZE STRESS*/
    initializestress(&(elem.head.z));
    initializestress(&(elem.tail.z));
    initializestress(&(elem.head.x));
    initializestress(&(elem.tail.x));
    initializestress(&(elem.head.y));
    initializestress(&(elem.tail.y));

    if(soffset && sect1.etype!=WALL && sect1.etype!=SLAB)
    {
      /*ELEMENT LENGTH*/
      h=100.0*elementlength(elem); /*[cm]*/

      if(sect1.etype==COLUMN)
      {
        h0[SX]=h-sect1.face[SX][HEAD]
                -sect1.face[SX][TAIL]; /*INNER LENGTH.*/
        h0[SY]=h-sect1.face[SY][HEAD]
                -sect1.face[SY][TAIL]; /*INNER LENGTH.*/
      }
      else
      {
        h0[SX]=h-sect1.face[SX][HEAD]
                -sect1.face[SX][TAIL]; /*INNER LENGTH.*/
      }

      /*MATERIAL*/
      gmaterial.sE=2100000.0; /*[kgf/cm2]*/
      gmaterial.rE=2100000.0;
      gmaterial.cE=gmaterial.rE/15.0;
      gmaterial.sF=3300.0*jis;                       /*STEEL SN490B*/
      gmaterial.Fc=240.0;                          /*CONCRETE Fc240*/

          /*MATERIAL FOR SHORT*/
          gmaterial.sft=3300.0;                       /*STEEL SM490*/
          gmaterial.sfc=-gmaterial.sft;
          gmaterial.sfb=3300.0; /*FOR SRC*/             /*b:BENDING*/
          gmaterial.sfs=3300.0/sqrt(3.0);                 /*s:SHEAR*/
          gmaterial.rft=3500.0;               /*REINFORCEMENT SD345*/
          gmaterial.wft=3000.0;
          gmaterial.rfc=-gmaterial.rft;
          gmaterial.cfc=-160.0;                          /*CONCRETE*/
          gmaterial.cfs=11.1;

          /*MATERIAL FOR ULTIMATE*/
          gmaterial.sftu=3300.0*jis;                  /*STEEL SM490*/
          gmaterial.sfcu=-gmaterial.sftu;
          gmaterial.rftu=3500.0*jis;          /*REINFORCEMENT SD345*/
          gmaterial.rfcu=-gmaterial.rftu;
          gmaterial.wfp=3000.0*jis;           /*REINFORCEMENT SD295*/

      if(sect1.stype==STYPE_S)
      {
        if(sect1.sform.type==STEEL_RECTS) gmaterial.sF=sect1.srect[0].F;
        else                              gmaterial.sF=sect1.sform.F;
      }

      ys.exp=1.5; /*1.0<exp<2.0*/

      Nmax=0.0;
      Nmin=0.0;

      /*AXIAL FORCE Nz*/
      if(sect1.stype==STYPE_S) /*S*/
      {
        Nmax=-1.5*allowabletensionofsteel(gmaterial.sF,sect1);
        Nmin=-1.5*allowablecompressionofsteel(gmaterial.sE,gmaterial.sF,
                                              h,h,sect1);
      }
      else if(sect1.stype==STYPE_RC) /*RC*/
      {
        ultimateaxialforceofsrc(elem,&Nmax,&Nmin,
                                &Ns,&Ms,&Nr,&Mr,&Nc,&Mc);
      }
      else if(sect1.stype==STYPE_SRC) /*SRC*/
      {
        ultimateaxialforceofsrc(elem,&Nmax,&Nmin,
                                &Ns,&Ms,&Nr,&Mr,&Nc,&Mc);
      }

      /*Nmax,Nmin:COMPRESSION=+ TENSION=-*/
      ys.fmax[0]=-Nmin/1000.0; /*[tf]*/
      ys.fmin[0]=-Nmax/1000.0;

      /*CENTER OF SURFACE : N=0.5(Nmax+Nmin)*/
      elem.head.x.N=0.5*(ys.fmax[0]+ys.fmin[0]);
      elem.tail.x.N=-elem.head.x.N;

      na=1;
      /*if(elem.sect->etype!=COLUMN) na=0;*/

      for(ia=0;ia<=na;ia++) /*FOR AXIS*/
      {
        uN=1000.0*(elem.head.x.N);

        /*BENDING Mx,My*/
        Mu=0.0;
        if(sect1.stype==STYPE_S)                                  /*S*/
        {
          Mu=ultimatebendingofsteel(sect1,gmaterial.sE,gmaterial.sF,
                                    h,h,ia,(-uN));            /*[kgfcm]*/
        }
        else if(sect1.stype==STYPE_RC)                           /*RC*/
        {
          Mu=ultimatebendingofsrc(elem,ia,(-uN),
                                  &Ns,&Ms,&Nr,&Mr,&Nc,&Mc); /*[kgfcm]*/
        }
        else if(sect1.stype==STYPE_SRC)                         /*SRC*/
        {
          Mu=ultimatebendingofsrc(elem,ia,(-uN),
                                  &Ns,&Ms,&Nr,&Mr,&Nc,&Mc); /*[kgfcm]*/
        }

        if(Mu<0.001) Mu=0.001;

        /*Mmax,Mmin*/
        ys.fmax[ia+4]=Mu/100000.0;
        ys.fmin[ia+4]=-Mu/100000.0;

        /*SHEAR Qx,Qy*/
        Qu=0.0;
        Qp[HEAD]=0.0;
        Qp[TAIL]=0.0;
        if(sect1.stype==STYPE_S)                      /*S*/
        {
          Qu=allowultimshearofsteel(PULTIMATE,gmaterial.sF,sect1,ia);
        }
        else if(sect1.stype==STYPE_RC)               /*RC*/
        {
		  Qu=ultimateshearofrc(elem,gmaterial,ia,
                               (-uN),0.0,0.0);
        }
        else if(sect1.stype==STYPE_SRC)             /*SRC*/
        {
          /*APPROXIMATION.Q,M MUST BE OF RC.*/
          /*
          if(ia==SX)
          {
            Qu=ultimateshearofsrc(elem,SX,HSTRONG,0.0,0.0);
          }
          if(ia==SY)
          {
            Qu=ultimateshearofsrc(elem,SY,HWEAK,0.0,0.0);
          }
          */
        }

        if(Qu<0.1) Qu=0.1;

        /*Qmax,Qmin*/
        ys.fmax[2-ia]=Qu/1000.0;
        ys.fmin[2-ia]=-Qu/1000.0;

      }/*ia*/

      /*TORTION Mz*/
      Mtu=0.0;
      if(sect1.stype==STYPE_S) /*S*/
      {
        /*UNDER CONSTRUCTION.*/
        if(ys.fmax[4]<=ys.fmax[5]) Mtu=ys.fmax[4]/sqrt(3.0); /*[tfm]*/
        else                       Mtu=ys.fmax[5]/sqrt(3.0); /*[tfm]*/
      }
      else if(sect1.stype==STYPE_RC) /*RC*/
      {
        Mtu=ultimatetorsionofrc(sect1);
        Mtu/=100000.0;
      }
      else if(sect1.stype==STYPE_SRC) /*SRC*/
      {
        /*Mtu=ultimatetortionofsrc();*/
      }

      /*Mtmax,Mtmin*/
      ys.fmax[3]=Mtu;
      ys.fmin[3]=-Mtu;

      sprintf(txt,"Elem=%d Sect=%d\n",elem.code,sect1.code);
      sprintf(non,"Nz: %8.3f %8.3f[tf] Nmid=%.3f\n",
              ys.fmax[0],ys.fmin[0],elem.head.x.N);
      strcat(txt,non);
      sprintf(non,"Qx: %8.3f %8.3f[tf]\n",ys.fmax[1],ys.fmin[1]);
      strcat(txt,non);
      sprintf(non,"Qy: %8.3f %8.3f[tf]\n",ys.fmax[2],ys.fmin[2]);
      strcat(txt,non);
      sprintf(non,"Mz: %8.3f %8.3f[tfm]\n",ys.fmax[3],ys.fmin[3]);
      strcat(txt,non);
      sprintf(non,"Mx: %8.3f %8.3f[tfm]\n",ys.fmax[4],ys.fmin[4]);
      strcat(txt,non);
      sprintf(non,"My: %8.3f %8.3f[tfm]\n",ys.fmax[5],ys.fmin[5]);
      strcat(txt,non);
      MessageBox(NULL,txt,"Surface",MB_OK);
    }

    if(soffset && sect1.etype==WALL)/*...........SHEAR OF WALL.*/
    {
      /*Qu=0.0;*/

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

      /*RATE OF WINDOW.*/
      wrate1=(elem.sect->wlength)/l0;
      wrate2=sqrt((elem.sect->wlength)*(elem.sect->wheight)
                  /(l0*h0[1]));
      if(wrate1>=wrate2) elem.sect->windowrate=wrate1;
      else               elem.sect->windowrate=wrate2;

      l0/=2.0; /*1 OF 2 BRACES.*/

      Qu=allowultimshearofrcwall(elem,l0,PULTIMATE);
      Nu=Qu*sqrt(l*l+h*h)/l;

      /*Nmax,Nmin*/
      ys.fmax[0]=Nu;
      ys.fmin[0]=-Nu;
    }

    if(soffset && sect1.etype==SLAB)               /*SHEAR OF SLAB.*/
    {
      /*UNDER CONSTRUCTION*/
    }

    (elems+ie)->yield=ys;
  }

  fclose(flist);

  MessageBox(NULL,"Completed.","Surface",MB_OK);
  return 0;
}/*createyieldsurface*/






/*****DEBUGGING*****/

int getsectionform(FILE *flist,int code,struct section *sect)
/*GET SECTION FROM SECTION LIST.*/
{
  char **data;
  int i,j,n,ns=0;
  int type,scode;
  int nsteel=0,nrein=0,nconc=0,nstrnd=0,ncomment=0;

  fseek(flist,0L,SEEK_SET);

  sect->code=0;                                   /*INITIALIZATION.*/
  sect->soff=0;
  sect->stype=0;
  sect->etype=0;

  sect->sform.type=STEEL_NULL;
  sect->sform.H =0.0;
  sect->sform.B =0.0;
  sect->sform.tw=0.0;
  sect->sform.tf=0.0;

  sect->nsteel=0;
  sect->nrein=0;
  sect->nconc=0;
  sect->nstrnd=0;
  sect->ncomment=0;

  sect->face[0][0]=0.0;
  sect->face[0][1]=0.0;
  sect->face[1][0]=0.0;
  sect->face[1][1]=0.0;

  sect->shearrein[0]=0.0;
  sect->shearrein[1]=0.0;
  sect->srein[0].n=0;
  sect->srein[1].n=0;

  sect->bblength[0]=0.0;
  sect->bblength[1]=0.0;
  sect->bbfact[0]=0.0;
  sect->bbfact[1]=0.0;
  sect->btlength[0]=0.0;
  sect->btlength[1]=0.0;
  sect->btfact[0]=0.0;
  sect->btfact[1]=0.0;

  for(i=0;i<MAXSTRND;i++)
  {
    sect->strnd[i].area[PEND]=0.0;
    sect->strnd[i].area[PMID]=0.0;
    sect->strnd[i].Ni[PEND]=0.0;
    sect->strnd[i].Ni[PMID]=0.0;
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
      scode=(int)strtol(*(data+1),NULL,10);

      if(scode==code)
      {
        sect->code=scode;

        if(!strncmp(*(data+2),"S ",2))   sect->stype=STYPE_S;
        if(!strncmp(*(data+2),"RC",2))   sect->stype=STYPE_RC;
        if(!strncmp(*(data+2),"SRC",3))  sect->stype=STYPE_SRC;
        if(!strncmp(*(data+2),"PC",2))   sect->stype=STYPE_PC;
        if(!strncmp(*(data+2),"WOOD",4)) sect->stype=STYPE_WOOD;
        if(!strncmp(*(data+2),"GLASS",5)) sect->stype=STYPE_GLASS;
        if(!strncmp(*(data+2),"ACRYL",5)) sect->stype=STYPE_ACRYL;
        if(!strncmp(*(data+2),"ALUMI",5)) sect->stype=STYPE_ALUMI;

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
              sect->sform.type=STEEL_RECTS;

              sect->srect[nsteel].left  =strtod(*(data+1),NULL);
              sect->srect[nsteel].bottom=strtod(*(data+2),NULL);
              sect->srect[nsteel].right =strtod(*(data+3),NULL);
              sect->srect[nsteel].top   =strtod(*(data+4),NULL);

              if(n>=6 && !strncmp(*(data+5),"SN400",5))
              {
                sect->srect[nsteel].F=2400.0*jis;
			  }
			  else if(n>=6 && !strncmp(*(data+5),"SN490",5))
			  {
				sect->srect[nsteel].F=3300.0*jis;
			  }
			  else if(n>=6 && !strncmp(*(data+5),"T40SN490",5))
			  {
				sect->srect[nsteel].F=3000.0*jis;
			  }
              else if(n>=6 && !strncmp(*(data+5),"HT590",5))
              {
                sect->srect[nsteel].F=4500.0;
              }
              else if(n>=6 && !strncmp(*(data+5),"HT930",5))
              {
                sect->srect[nsteel].F=8000.0;
              }
              else if(n>=6 && !strncmp(*(data+5),"HT950",5))
              {
                sect->srect[nsteel].F=9000.0;
              }
              /***********************************/
              else if(n>=6 && !strncmp(*(data+5),"AS175",5))   /*Alumi*/
              {
                sect->srect[nsteel].F=1750.0;
              }
              else if(n>=6 && !strncmp(*(data+5),"A7178",5))   /*Alumi*/
              {
                sect->srect[nsteel].F=4400.0;
              }

              else sect->srect[nsteel].F=2400.0*jis;

              nsteel++;
            }
            if(!strncmp(*(data+0),"HKYOU",5))
            {
              sect->sform.type=STEEL_HKYOU;

              sect->sform.H =strtod(*(data+1),NULL);
              sect->sform.B =strtod(*(data+2),NULL);
              sect->sform.tw=strtod(*(data+3),NULL);
              sect->sform.tf=strtod(*(data+4),NULL);

              if(n>=6 && !strncmp(*(data+5),"SN400",5))
              {
                sect->sform.F=2400.0*jis;
              }
              else if(n>=6 && !strncmp(*(data+5),"SN490",5))
              {
                sect->sform.F=3300.0*jis;
			  }
			  else if(n>=6 && !strncmp(*(data+5),"T40SN490",5))
              {
                sect->sform.F=3000.0*jis;
			  }
              else if(n>=6 && !strncmp(*(data+5),"HT590",5))
              {
                sect->sform.F=4500.0;
              }
              else if(n>=6 && !strncmp(*(data+5),"HT930",5))
              {
                sect->sform.F=8000.0;
              }
              else if(n>=6 && !strncmp(*(data+5),"HT950",5))
              {
                sect->sform.F=9000.0;
              }
              else if(n>=6 && !strncmp(*(data+5),"AS175",5))   /*Alumi*/
              {
                sect->sform.F=1750.0;
              }
              else if(n>=6 && !strncmp(*(data+5),"A7178",5))   /*Alumi*/
              {
                sect->sform.F=4400.0;
              }

              else sect->sform.F=2400.0*jis;
            }
            if(!strncmp(*(data+0),"HWEAK",5))
            {
              sect->sform.type=STEEL_HWEAK;

              sect->sform.H =strtod(*(data+1),NULL);
              sect->sform.B =strtod(*(data+2),NULL);
              sect->sform.tw=strtod(*(data+3),NULL);
              sect->sform.tf=strtod(*(data+4),NULL);

              if(n>=6 && !strncmp(*(data+5),"SN400",5))
              {
                sect->sform.F=2400.0*jis;
              }
			  else if(n>=6 && !strncmp(*(data+5),"SN490",5))
			  {
				sect->sform.F=3300.0*jis;
			  }
			  else if(n>=6 && !strncmp(*(data+5),"T40SN490",5))
			  {
				sect->sform.F=3000.0*jis;
			  }
			  else if(n>=6 && !strncmp(*(data+5),"HT590",5))
              {
                sect->sform.F=4500.0;
              }
              else if(n>=6 && !strncmp(*(data+5),"HT930",5))
              {
                sect->sform.F=8000.0;
              }
              else if(n>=6 && !strncmp(*(data+5),"HT950",5))
              {
                sect->sform.F=9000.0;
              }
              else if(n>=6 && !strncmp(*(data+5),"AS175",5))   /*Alumi*/
              {
                sect->sform.F=1750.0;
              }
              else if(n>=6 && !strncmp(*(data+5),"A7178",5))   /*Alumi*/
              {
                sect->sform.F=4400.0;
              }

              else sect->sform.F=2400.0*jis;
            }
            if(!strncmp(*(data+0),"RPIPE",5))
            {
              sect->sform.type=STEEL_RPIPE;

              sect->sform.H =strtod(*(data+1),NULL);
              sect->sform.B =strtod(*(data+2),NULL);
              sect->sform.tw=strtod(*(data+3),NULL);
              sect->sform.tf=strtod(*(data+4),NULL);

              if(n>=6 && !strncmp(*(data+5),"SN400",5))
              {
                sect->sform.F=2400.0*jis;
              }
			  else if(n>=6 && !strncmp(*(data+5),"SN490",5))
			  {
				sect->sform.F=3300.0*jis;
			  }
			  else if(n>=6 && !strncmp(*(data+5),"T40SN490",5))
			  {
				sect->sform.F=3000.0*jis;
			  }
              else if(n>=6 && !strncmp(*(data+5),"HT590",5))
              {
                sect->sform.F=4500.0;
              }
              else if(n>=6 && !strncmp(*(data+5),"HT930",5))
              {
                sect->sform.F=8000.0;
              }
              else if(n>=6 && !strncmp(*(data+5),"HT950",5))
              {
                sect->sform.F=9000.0;
              }
              /***********************************/
              else if(n>=6 && !strncmp(*(data+5),"AS175",5))   /*Alumi*/
              {
                sect->sform.F=1750.0;
              }
              else if(n>=6 && !strncmp(*(data+5),"A7178",5))   /*Alumi*/
              {
                sect->sform.F=4400.0;
              }
              else sect->sform.F=2400.0*jis;
            }
            if(!strncmp(*(data+0),"CPIPE",5))
            {
              sect->sform.type=STEEL_CPIPE;

              sect->sform.H =strtod(*(data+1),NULL);
              sect->sform.B =sect->sform.H;
              sect->sform.tw=strtod(*(data+2),NULL);
              sect->sform.tf=sect->sform.tw;

              if(n>=4 && !strncmp(*(data+3),"SN400",5))
              {
                sect->sform.F=2400.0*jis;
              }
			  else if(n>=4 && !strncmp(*(data+3),"SN490",5))
			  {
				sect->sform.F=3300.0*jis;
			  }
			  else if(n>=4 && !strncmp(*(data+3),"T40SN490",5))
			  {
				sect->sform.F=3000.0*jis;
			  }
              else if(n>=4 && !strncmp(*(data+3),"HT590",5))
              {
                sect->sform.F=4500.0;
              }
              else if(n>=4 && !strncmp(*(data+3),"HT930",5))
              {
                sect->sform.F=8000.0;
              }
              else if(n>=4 && !strncmp(*(data+3),"HT950",5))
              {
                sect->sform.F=9000.0;
              }
              /***********************************/
              else if(n>=4 && !strncmp(*(data+3),"AS175",5))   /*Alumi*/
              {
                sect->sform.F=1750.0;
              }
              else if(n>=6 && !strncmp(*(data+3),"A7178",5))   /*Alumi*/
              {
                sect->sform.F=4400.0;
              }
              else sect->sform.F=2400.0*jis;
            }
            if(!strncmp(*(data+0),"PLATE",5))
            {
              sect->sform.type=STEEL_PLATE;

              sect->sform.H =strtod(*(data+1),NULL);
              sect->sform.B =strtod(*(data+2),NULL);

              if(n>=4 && !strncmp(*(data+3),"SN400",5))
              {
                sect->sform.F=2400.0*jis;
              }
			  else if(n>=4 && !strncmp(*(data+3),"SN490",5))
			  {
				sect->sform.F=3300.0*jis;
			  }
			  else if(n>=4 && !strncmp(*(data+3),"T40SN490",5))
			  {
				sect->sform.F=3000.0*jis;
			  }
              else if(n>=4 && !strncmp(*(data+3),"HT590",5))
              {
                sect->sform.F=4500.0;
              }
              else if(n>=4 && !strncmp(*(data+3),"HT930",5))
              {
                sect->sform.F=8000.0;
              }
              else if(n>=4 && !strncmp(*(data+3),"HT950",5))
              {
                sect->sform.F=9000.0;
              }
              else if(n>=4 && !strncmp(*(data+3),"SUGI",4))
              {
                sect->wform.Fc=220.2; /*[kgf/cm2]*/
                sect->wform.Ft=165.1;
                sect->wform.Fb=275.3;
                sect->wform.Fs= 18.3;
              }
              else if(n>=4 && !strncmp(*(data+3),"E70SUGI",7))
              {
                sect->wform.Fc=238.6; /*[kgf/cm2]*/
                sect->wform.Ft=177.4;
                sect->wform.Fb=299.7;
                sect->wform.Fs= 18.3;
              }
              else if(n>=4 && !strncmp(*(data+3),"E90SUGI",7))
              {
                sect->wform.Fc=287.5; /*[kgf/cm2]*/
                sect->wform.Ft=214.1;
                sect->wform.Fb=354.8;
                sect->wform.Fs= 18.3;
              }

              else if(n>=4 && !strncmp(*(data+3),"M-E90",5))
              {
				sect->wform.Fc=171.4; /*[kgf/cm2]*/
				sect->wform.Ft=128.5;
				sect->wform.Fb=214.2;
				sect->wform.Fs= 24.4;
			  }
			  else if(n>=4 && !strncmp(*(data+3),"M-E110",6))
              {
				sect->wform.Fc=250.8; /*[kgf/cm2]*/
				sect->wform.Ft=189.6;
				sect->wform.Fb=312.0;
				sect->wform.Fs= 24.4;
              }
              else if(n>=4 && !strncmp(*(data+3),"E110HINOKI",10))
              {
                sect->wform.Fc=318.1; /*[kgf/cm2]*/
                sect->wform.Ft=238.6;
                sect->wform.Fb=391.5;
                sect->wform.Fs= 21.4;
              }
              else if(n>=4 && !strncmp(*(data+3),"E120-F330",9))
              {
                sect->wform.Fc=256.9; /*[kgf/cm2]*/
                sect->wform.Ft=226.3;
                sect->wform.Fb=330.3;
                sect->wform.Fs= 36.7;
              }
              else if(n>=4 && !strncmp(*(data+3),"E65-F225",8))
              {
                sect->wform.Fc=171.3; /*[kgf/cm2]*/
				sect->wform.Ft=146.8;
				sect->wform.Fb=226.3;
				sect->wform.Fs= 30.5;
			  }
			  else if(n>=4 && !strncmp(*(data+3),"STICK",5))
			  {
				sect->wform.Fc=0.2; /*[kgf/cm2]*/
				sect->wform.Ft=0.2;
				sect->wform.Fb=1.0;
				sect->wform.Fs=0.2;
			  }
			  else if(n>=4 && !strncmp(*(data+3),"PLYWOOD",7))
			  {
				sect->wform.Fc=135.0; /*[kgf/cm2]*/
				sect->wform.Ft=195.0;
				sect->wform.Fb=60.0;
				sect->wform.Fs=42.0;
			  }
			  else if(n>=4 && !strncmp(*(data+3),"GLASS",5))
              {
                sect->wform.Fc=6000.0; /* =9000/1.5 [kgf/cm2]*/
                sect->wform.Ft= 383.3; /* = 575/1.5 [kgf/cm2]*/
                sect->wform.Fb= 416.6; /* = 625/1.5 [kgf/cm2]*/
                sect->wform.Fs=  38.3; /* =  Ft/10  [kgf/cm2]*/
              }
              else if(n>=4 && !strncmp(*(data+3),"ACRYL",5))
              {
                sect->wform.Fc= 833.3; /* =1250/1.5 [kgf/cm2]*/
                sect->wform.Ft= 400.0; /* = 600/1.5 [kgf/cm2]*/
                sect->wform.Fb= 666.6; /* =1000/1.5 [kgf/cm2]*/
                sect->wform.Fs=  40.0; /* =  Ft/10  [kgf/cm2]*/
              }
              else if(n>=4 && !strncmp(*(data+3),"AS175",5))   /*Alumi*/
              {
                sect->sform.F=1750.0;
			  }
              else if(n>=4 && !strncmp(*(data+3),"A7178",5))   /*Alumi*/
              {
                sect->sform.F=4400.0;
              }
              else sect->sform.F=2400.0*jis;
            }
            if(!strncmp(*(data+0),"DABO",4))
            {
              sect->sform.type=DABO;

              sect->sform.H =strtod(*(data+1),NULL);
              sect->sform.B =strtod(*(data+2),NULL);

              if(n>=4 && !strncmp(*(data+3),"SUGI",4))
              {
                sect->wform.Fc=220.2; /*[kgf/cm2]*/
                sect->wform.Ft=165.1;
                sect->wform.Fb=275.3;
                sect->wform.Fs= 18.3;
              }
            }
            if(!strncmp(*(data+0),"ANGLE",5))
            {
              sect->sform.type=STEEL_ANGLE;

              sect->sform.H =strtod(*(data+1),NULL);
              sect->sform.B =strtod(*(data+2),NULL);
              sect->sform.tw=strtod(*(data+3),NULL);
              sect->sform.tf=strtod(*(data+4),NULL);

              if(n>=6 && !strncmp(*(data+5),"SN400",5))
              {
                sect->sform.F=2400.0*jis;
              }
			  else if(n>=6 && !strncmp(*(data+5),"SN490",5))
			  {
				sect->sform.F=3300.0*jis;
			  }
			  else if(n>=6 && !strncmp(*(data+5),"T40SN490",5))
			  {
				sect->sform.F=3000.0*jis;
			  }
              else if(n>=6 && !strncmp(*(data+5),"HT590",5))
              {
                sect->sform.F=4500.0;
              }
              else if(n>=6 && !strncmp(*(data+5),"HT930",5))
              {
                sect->sform.F=8000.0;
              }
              else if(n>=6 && !strncmp(*(data+5),"HT950",5))
              {
                sect->sform.F=9000.0;
              }
              /***********************************/
              else if(n>=6 && !strncmp(*(data+5),"AS175",5))   /*Alumi*/
              {
                sect->sform.F=1750.0;
              }
              else if(n>=6 && !strncmp(*(data+5),"A7178",5))   /*Alumi*/
              {
                sect->sform.F=4400.0;
              }
              else sect->sform.F=2400.0*jis;
            }
            if(!strncmp(*(data+0),"TKYOU",5))
            {
              sect->sform.type=STEEL_TKYOU;

              sect->sform.H =strtod(*(data+1),NULL);
              sect->sform.B =strtod(*(data+2),NULL);
              sect->sform.tw=strtod(*(data+3),NULL);
              sect->sform.tf=strtod(*(data+4),NULL);

              if(n>=6 && !strncmp(*(data+5),"SN400",5))
              {
                sect->sform.F=2400.0*jis;
              }
			  else if(n>=6 && !strncmp(*(data+5),"SN490",5))
			  {
				sect->sform.F=3300.0*jis;
			  }
			  else if(n>=6 && !strncmp(*(data+5),"T40SN490",5))
			  {
				sect->sform.F=3000.0*jis;
			  }
			  else if(n>=6 && !strncmp(*(data+5),"HT590",5))
              {
                sect->sform.F=4500.0;
              }
              else if(n>=6 && !strncmp(*(data+5),"HT930",5))
              {
                sect->sform.F=8000.0;
              }
              else if(n>=6 && !strncmp(*(data+5),"HT950",5))
              {
                sect->sform.F=9000.0;
              }
              /***********************************/
              else if(n>=6 && !strncmp(*(data+5),"AS175",5))   /*Alumi*/
              {
                sect->sform.F=1750.0;
              }
              else if(n>=6 && !strncmp(*(data+5),"A7178",5))   /*Alumi*/
              {
                sect->sform.F=4400.0;
              }
              else sect->sform.F=2400.0*jis;
            }
            if(!strncmp(*(data+0),"TWEAK",5))
            {
              sect->sform.type=STEEL_TWEAK;

              sect->sform.H =strtod(*(data+1),NULL);
              sect->sform.B =strtod(*(data+2),NULL);
              sect->sform.tw=strtod(*(data+3),NULL);
              sect->sform.tf=strtod(*(data+4),NULL);

              if(n>=6 && !strncmp(*(data+5),"SN400",5))
              {
                sect->sform.F=2400.0*jis;
              }
			  else if(n>=6 && !strncmp(*(data+5),"SN490",5))
			  {
				sect->sform.F=3300.0*jis;
			  }
			  else if(n>=6 && !strncmp(*(data+5),"T40SN490",5))
			  {
				sect->sform.F=3000.0*jis;
			  }
              else if(n>=6 && !strncmp(*(data+5),"HT590",5))
              {
                sect->sform.F=4500.0;
              }
              else if(n>=6 && !strncmp(*(data+5),"HT930",5))
              {
                sect->sform.F=8000.0;
              }
              else if(n>=6 && !strncmp(*(data+5),"HT950",5))
              {
                sect->sform.F=9000.0;
              }
              else if(n>=6 && !strncmp(*(data+5),"AS175",5))   /*Alumi*/
              {
                sect->sform.F=1750.0;
              }
              else if(n>=6 && !strncmp(*(data+5),"A7178",5))   /*Alumi*/
              {
                sect->sform.F=4400.0;
              }
              else sect->sform.F=2400.0*jis;
            }
            if(!strncmp(*(data+0),"CRECT",4))
            {
              sect->crect[nconc].left  =strtod(*(data+1),NULL);
              sect->crect[nconc].bottom=strtod(*(data+2),NULL);
              sect->crect[nconc].right =strtod(*(data+3),NULL);
              sect->crect[nconc].top   =strtod(*(data+4),NULL);

			  if(n>=6 && !strncmp(*(data+5),"FC21",4))
			  {
				sect->crect[nconc].F=210.0;
			  }
			  else if(n>=6 && !strncmp(*(data+5),"FC24",4))
			  {
				sect->crect[nconc].F=240.0;
			  }
			  else if(n>=6 && !strncmp(*(data+5),"FC27",4))
              {
                sect->crect[nconc].F=270.0;
              }
              else if(n>=6 && !strncmp(*(data+5),"FC50",4))
              {
                sect->crect[nconc].F=500.0;
              }
              /***********************************/
              else if(n>=6 && !strncmp(*(data+5),"AS175",5))   /*Alumi*/
              {
                sect->sform.F=1750.0;
              }
              else if(n>=6 && !strncmp(*(data+5),"A7178",5))   /*Alumi*/
              {
                sect->sform.F=4400.0;
              }
              else sect->crect[nconc].F=240.0;

              nconc++;
            }
            if(!strncmp(*(data+0),"REINS",5))
            {
              sect->rein[nrein].area=strtod(*(data+1),NULL);
              sect->rein[nrein].x=strtod(*(data+2),NULL);
              sect->rein[nrein].y=strtod(*(data+3),NULL);

              if(n>=5 && !strncmp(*(data+4),"SD295",5))
              {
                sect->rein[nrein].F=3000.0*jis;
              }
              else if(n>=5 && !strncmp(*(data+4),"SD345",5))
              {
                sect->rein[nrein].F=3500.0*jis;
              }
              else sect->rein[nrein].F=3000.0*jis;

              nrein++;
            }
            if(!strncmp(*(data+0),"HOOPS",5))
            {
              sect->shearrein[0]=strtod(*(data+1),NULL);
              sect->shearrein[1]=strtod(*(data+2),NULL);

              if(n>=4 && !strncmp(*(data+3),"SD295",5))
              {
                sect->wF=3000.0*jis;
              }
              else if(n>=4 && !strncmp(*(data+3),"SD345",5))
              {
                sect->wF=3500.0*jis;
              }
              else sect->wF=3000.0*jis;
            }
            if(!strncmp(*(data+0),"XHOOP",5))
            {
              sect->srein[0].area =strtod(*(data+1),NULL);
              sect->srein[0].n    =(int)strtol(*(data+2),NULL,10);
              sect->srein[0].pitch=strtod(*(data+3),NULL);

              if(n>=5 && !strncmp(*(data+4),"SD295",5))
              {
                sect->srein[0].F=3000.0*jis;
              }
              else if(n>=5 && !strncmp(*(data+4),"SD345",5))
              {
                sect->srein[0].F=3500.0*jis;
              }
              else sect->srein[0].F=3000.0*jis;
            }
            if(!strncmp(*(data+0),"YHOOP",5))
            {
              sect->srein[1].area =strtod(*(data+1),NULL);
              sect->srein[1].n    =(int)strtol(*(data+2),NULL,10);
              sect->srein[1].pitch=strtod(*(data+3),NULL);

              if(n>=5 && !strncmp(*(data+4),"SD295",5))
              {
                sect->srein[1].F=3000.0*jis;
              }
              else if(n>=5 && !strncmp(*(data+4),"SD345",5))
              {
                sect->srein[1].F=3500.0*jis;
              }
              else sect->srein[1].F=3000.0*jis;
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
                sect->strnd[nstrnd].area[PEND]=strtod(*(data+3),NULL);
                sect->strnd[nstrnd].x[PEND]   =strtod(*(data+4),NULL);
                sect->strnd[nstrnd].y[PEND]   =strtod(*(data+5),NULL);
                sect->strnd[nstrnd].Ni[PEND]  =strtod(*(data+6),NULL);
              }
              if(!strncmp(*(data+1),"ALL",3) ||
                 !strncmp(*(data+1),"MID",3))
              {
                sect->strnd[nstrnd].area[PMID]=strtod(*(data+3),NULL);
                sect->strnd[nstrnd].x[PMID]   =strtod(*(data+4),NULL);
                sect->strnd[nstrnd].y[PMID]   =strtod(*(data+5),NULL);
                sect->strnd[nstrnd].Ni[PMID]  =strtod(*(data+6),NULL);
              }
              nstrnd++;
            }

            if(!strncmp(*(data+0),"BBLEN",5))
            {
              sect->bblength[0]=strtod(*(data+1),NULL);
              sect->bblength[1]=strtod(*(data+2),NULL);
            }
            if(!strncmp(*(data+0),"BBFAC",5))
            {
              sect->bbfact[0]=strtod(*(data+1),NULL);
              sect->bbfact[1]=strtod(*(data+2),NULL);
            }
            if(!strncmp(*(data+0),"BTLEN",5))
            {
              sect->btlength[0]=strtod(*(data+1),NULL);
              sect->btlength[1]=strtod(*(data+2),NULL);
            }
            if(!strncmp(*(data+0),"BTFAC",5))
            {
              sect->btfact[0]=strtod(*(data+1),NULL);
              sect->btfact[1]=strtod(*(data+2),NULL);
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
        if(sect->etype==BRACE)
        {
          while(1)
          {
            /*data=fgetsbrk(flist,&n);*/
            data=fgetscut(flist,&n);
            if(n==0) break;
            if(!strncmp(*(data+0),"CODE",4)) break;

            if(!strncmp(*(data+0),"SAREA",5))
            {
              sect->thick=strtod(*(data+1),NULL);

              if(n>=3 && !strncmp(*(data+2),"SN400",5))
              {
                sect->cF=2400.0*jis;
              }
			  else if(n>=3 && !strncmp(*(data+2),"SN490",5))
			  {
				sect->cF=3300.0*jis;
			  }
			  else if(n>=3 && !strncmp(*(data+2),"T40SN490",5))
			  {
				sect->cF=3000.0*jis;
			  }
              else if(n>=3 && !strncmp(*(data+2),"HT590",5))
              {
                sect->cF=4500.0;
              }
              else if(n>=3 && !strncmp(*(data+2),"HT930",5))
              {
                sect->cF=8000.0;
              }
              else if(n>=3 && !strncmp(*(data+2),"HT950",5))
              {
                sect->sform.F=9000.0;
              }
              else if(n>=3 && !strncmp(*(data+2),"AS175",5))
              {
                sect->cF=1750.0;
              }
              else if(n>=3 && !strncmp(*(data+2),"A7178",5))
              {
                sect->cF=4400.0;
              }

              else sect->cF=2400.0*jis;
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

			  if(n>=3 && !strncmp(*(data+2),"FC21",4))
			  {
				sect->cF=210.0;
			  }
			  else if(n>=3 && !strncmp(*(data+2),"FC24",4))
			  {
				sect->cF=240.0;
			  }
			  else if(n>=3 && !strncmp(*(data+2),"FC27",4))
			  {
				sect->cF=270.0;
			  }
			  else if(n>=3 && !strncmp(*(data+2),"GOHAN",5))
			  {
				sect->fsl=12.0; /*CLASS 1,TYPE C,90[DO] [kgf/cm2]*/
			  }
			  else if(n>=3 && !strncmp(*(data+2),"GLASS",5))
			  {
				sect->fsl=25.5; /* = Fs/1.5 [kgf/cm2]*/
			  }
			  else if(n>=3 && !strncmp(*(data+2),"ACRYL",5))
			  {
				sect->fsl=26.6; /* = Fs/1.5 [kgf/cm2]*/
			  }
			  else sect->cF=240.0;
			}
			if(!strncmp(*(data+0),"SREIN",5))
			{
			  sect->shearrein[0]=strtod(*(data+1),NULL);

			  if(n>=3 && !strncmp(*(data+2),"SD295",5))
			  {
				sect->wF=3000.0*jis;
			  }
			  else if(n>=3 && !strncmp(*(data+2),"SD345",5))
			  {
				sect->wF=3500.0*jis;
              }
              else sect->wF=3000.0*jis;
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
  SIZE size;
  char str[256];
  int i;
  int commentwidth,commentheight;
  long int top,height;
  double blank=20;
  double Xmin,Xmax,Ymin,Ymax;

  if(oy>pageheight) return;
  top=oy;

  if(mode==ONPRINTER) StartPage(hdc);

  if(mode==ONSCREEN)
  {
    SetTextColor(hdc,RGB(255,255,255));
    sprintf(str,"Section List");

    TextOut(hdc,ox,top,str,strlen(str));
  }
  else if(mode==ONPRINTER)
  {
    setfontformat(hdc,96,36,"ÇlÇr ñæí©",0,0,0);

    /*SetTextColor(hdc,RGB(0,0,0));*/
	sprintf(str,"4.2 : âºíËífñ ");
	/*sprintf(str,"2.4.2 : âºíËífñ ");*/
    /*sprintf(str,"Section List");*/

    TextOut(hdc,ox,top,str,strlen(str));
    setfontformat(hdc,80,30,"ÇlÇr ñæí©",0,0,0);
    top+=100;
  }

  GetTextExtentPoint32(hdc,str,strlen(str),&size);
  top+=size.cy+(long int)(vp.gfactor*blank);

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

  if(sect->sform.type==STEEL_HKYOU ||
     sect->sform.type==STEEL_RPIPE ||
     sect->sform.type==STEEL_CPIPE ||
     sect->sform.type==STEEL_PLATE ||
     sect->sform.type==STEEL_ANGLE ||
     sect->sform.type==STEEL_TKYOU)
  {
    if(*Xmin > -(sect->sform.B)/2.0) *Xmin=-(sect->sform.B)/2.0;
    if(*Xmax <  (sect->sform.B)/2.0) *Xmax= (sect->sform.B)/2.0;

    if(*Ymin > -(sect->sform.H)/2.0) *Ymin=-(sect->sform.H)/2.0;
    if(*Ymax <  (sect->sform.H)/2.0) *Ymax= (sect->sform.H)/2.0;
  }
  else if(sect->sform.type==STEEL_HWEAK ||
          sect->sform.type==STEEL_TWEAK)
  {
    if(*Xmin > -(sect->sform.H)/2.0) *Xmin=-(sect->sform.H)/2.0;
    if(*Xmax <  (sect->sform.H)/2.0) *Xmax= (sect->sform.H)/2.0;

    if(*Ymin > -(sect->sform.B)/2.0) *Ymin=-(sect->sform.B)/2.0;
    if(*Ymax <  (sect->sform.B)/2.0) *Ymax= (sect->sform.B)/2.0;
  }
  else if(sect->sform.type==STEEL_RECTS)
  {
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
  double textleft=210.0;
  double Xmin,Xmax,Ymin,Ymax,H,B,tw,tf;

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

  ppen=(HPEN)SelectObject(hdc,hpenS);
  if(sect->sform.type==STEEL_HKYOU)
  {
    H =sect->sform.H;
    B =sect->sform.B;
    tw=sect->sform.tw;
    tf=sect->sform.tf;

    MoveToEx(hdc,ox-(int)(vp.gfactor*B/2.0),
                 oy-(int)(vp.gfactor*H/2.0),NULL);
    LineTo(hdc,ox+(int)(vp.gfactor*B/2.0),
               oy-(int)(vp.gfactor*H/2.0));
    LineTo(hdc,ox+(int)(vp.gfactor*B/2.0),
               oy-(int)(vp.gfactor*(H/2.0-tf)));
    LineTo(hdc,ox+(int)(vp.gfactor*tw/2.0),
               oy-(int)(vp.gfactor*(H/2.0-tf)));
    LineTo(hdc,ox+(int)(vp.gfactor*tw/2.0),
               oy+(int)(vp.gfactor*(H/2.0-tf)));
    LineTo(hdc,ox+(int)(vp.gfactor*B/2.0),
               oy+(int)(vp.gfactor*(H/2.0-tf)));
    LineTo(hdc,ox+(int)(vp.gfactor*B/2.0),
               oy+(int)(vp.gfactor*H/2.0));
    LineTo(hdc,ox-(int)(vp.gfactor*B/2.0),
               oy+(int)(vp.gfactor*H/2.0));
    LineTo(hdc,ox-(int)(vp.gfactor*B/2.0),
               oy+(int)(vp.gfactor*(H/2.0-tf)));
    LineTo(hdc,ox-(int)(vp.gfactor*tw/2.0),
               oy+(int)(vp.gfactor*(H/2.0-tf)));
    LineTo(hdc,ox-(int)(vp.gfactor*tw/2.0),
               oy-(int)(vp.gfactor*(H/2.0-tf)));
    LineTo(hdc,ox-(int)(vp.gfactor*B/2.0),
               oy-(int)(vp.gfactor*(H/2.0-tf)));
    LineTo(hdc,ox-(int)(vp.gfactor*B/2.0),
               oy-(int)(vp.gfactor*H/2.0));
  }
  else if(sect->sform.type==STEEL_HWEAK)
  {
    H =sect->sform.H;
    B =sect->sform.B;
    tw=sect->sform.tw;
    tf=sect->sform.tf;

    MoveToEx(hdc,ox-(int)(vp.gfactor*H/2.0),
                 oy-(int)(vp.gfactor*B/2.0),NULL);
    LineTo(hdc,ox-(int)(vp.gfactor*(H/2.0-tf)),
               oy-(int)(vp.gfactor*B/2.0));
    LineTo(hdc,ox-(int)(vp.gfactor*(H/2.0-tf)),
               oy-(int)(vp.gfactor*tw/2.0));
    LineTo(hdc,ox+(int)(vp.gfactor*(H/2.0-tf)),
               oy-(int)(vp.gfactor*tw/2.0));
    LineTo(hdc,ox+(int)(vp.gfactor*(H/2.0-tf)),
               oy-(int)(vp.gfactor*B/2.0));
    LineTo(hdc,ox+(int)(vp.gfactor*H/2.0),
               oy-(int)(vp.gfactor*B/2.0));
    LineTo(hdc,ox+(int)(vp.gfactor*H/2.0),
               oy+(int)(vp.gfactor*B/2.0));
    LineTo(hdc,ox+(int)(vp.gfactor*(H/2.0-tf)),
               oy+(int)(vp.gfactor*B/2.0));
    LineTo(hdc,ox+(int)(vp.gfactor*(H/2.0-tf)),
               oy+(int)(vp.gfactor*tw/2.0));
    LineTo(hdc,ox-(int)(vp.gfactor*(H/2.0-tf)),
               oy+(int)(vp.gfactor*tw/2.0));
    LineTo(hdc,ox-(int)(vp.gfactor*(H/2.0-tf)),
               oy+(int)(vp.gfactor*B/2.0));
    LineTo(hdc,ox-(int)(vp.gfactor*H/2.0),
               oy+(int)(vp.gfactor*B/2.0));
    LineTo(hdc,ox-(int)(vp.gfactor*H/2.0),
               oy-(int)(vp.gfactor*B/2.0));
  }
  else if(sect->sform.type==STEEL_RPIPE)
  {
    H =sect->sform.H;
    B =sect->sform.B;
    tw=sect->sform.tw;
    tf=sect->sform.tf;

    MoveToEx(hdc,ox-(int)(vp.gfactor*B/2.0),
                 oy-(int)(vp.gfactor*H/2.0),NULL);
    LineTo(hdc,ox+(int)(vp.gfactor*B/2.0),
               oy-(int)(vp.gfactor*H/2.0));
    LineTo(hdc,ox+(int)(vp.gfactor*B/2.0),
               oy+(int)(vp.gfactor*H/2.0));
    LineTo(hdc,ox-(int)(vp.gfactor*B/2.0),
               oy+(int)(vp.gfactor*H/2.0));
    LineTo(hdc,ox-(int)(vp.gfactor*B/2.0),
               oy-(int)(vp.gfactor*H/2.0));

    MoveToEx(hdc,ox-(int)(vp.gfactor*(B/2.0-tw)),
                 oy-(int)(vp.gfactor*(H/2.0-tf)),NULL);
    LineTo(hdc,ox+(int)(vp.gfactor*(B/2.0-tw)),
               oy-(int)(vp.gfactor*(H/2.0-tf)));
    LineTo(hdc,ox+(int)(vp.gfactor*(B/2.0-tw)),
               oy+(int)(vp.gfactor*(H/2.0-tf)));
    LineTo(hdc,ox-(int)(vp.gfactor*(B/2.0-tw)),
               oy+(int)(vp.gfactor*(H/2.0-tf)));
    LineTo(hdc,ox-(int)(vp.gfactor*(B/2.0-tw)),
               oy-(int)(vp.gfactor*(H/2.0-tf)));
  }
  else if(sect->sform.type==STEEL_CPIPE)
  {
    H =sect->sform.H;
    B =sect->sform.B;
    tw=sect->sform.tw;
    tf=sect->sform.tf;

    Arc(hdc,
        (int)ox-(int)(vp.gfactor*B/2.0),
        (int)oy-(int)(vp.gfactor*H/2.0),
        (int)ox+(int)(vp.gfactor*B/2.0),
        (int)oy+(int)(vp.gfactor*H/2.0),
        0,0,0,0); /*HOLLOW CIRCLE.*/
    Arc(hdc,
        (int)ox-(int)(vp.gfactor*(B/2.0-tw)),
        (int)oy-(int)(vp.gfactor*(H/2.0-tf)),
        (int)ox+(int)(vp.gfactor*(B/2.0-tw)),
        (int)oy+(int)(vp.gfactor*(H/2.0-tf)),
        0,0,0,0); /*HOLLOW CIRCLE.*/
  }
  else if(sect->sform.type==STEEL_PLATE)
  {
    H =sect->sform.H;
    B =sect->sform.B;

    MoveToEx(hdc,ox-(int)(vp.gfactor*B/2.0),
                 oy-(int)(vp.gfactor*H/2.0),NULL);
    LineTo(hdc,ox+(int)(vp.gfactor*B/2.0),
               oy-(int)(vp.gfactor*H/2.0));
    LineTo(hdc,ox+(int)(vp.gfactor*B/2.0),
               oy+(int)(vp.gfactor*H/2.0));
    LineTo(hdc,ox-(int)(vp.gfactor*B/2.0),
               oy+(int)(vp.gfactor*H/2.0));
    LineTo(hdc,ox-(int)(vp.gfactor*B/2.0),
               oy-(int)(vp.gfactor*H/2.0));
  }
  else if(sect->sform.type==STEEL_ANGLE)
  {
    H =sect->sform.H;
    B =sect->sform.B;
    tw=sect->sform.tw;
    tf=sect->sform.tf;

    MoveToEx(hdc,ox-(int)(vp.gfactor*B/2.0),
                 oy-(int)(vp.gfactor*H/2.0),NULL);
    LineTo(hdc,ox-(int)(vp.gfactor*(B/2.0-tw)),
               oy-(int)(vp.gfactor*H/2.0));
    LineTo(hdc,ox-(int)(vp.gfactor*(B/2.0-tw)),
               oy+(int)(vp.gfactor*(H/2.0-tf)));
    LineTo(hdc,ox+(int)(vp.gfactor*B/2.0),
               oy+(int)(vp.gfactor*(H/2.0-tf)));
    LineTo(hdc,ox+(int)(vp.gfactor*B/2.0),
               oy+(int)(vp.gfactor*H/2.0));
    LineTo(hdc,ox-(int)(vp.gfactor*B/2.0),
               oy+(int)(vp.gfactor*H/2.0));
    LineTo(hdc,ox-(int)(vp.gfactor*B/2.0),
               oy-(int)(vp.gfactor*H/2.0));
  }
  else if(sect->sform.type==STEEL_TKYOU)
  {
    H =sect->sform.H;
    B =sect->sform.B;
    tw=sect->sform.tw;
    tf=sect->sform.tf;

    MoveToEx(hdc,ox-(int)(vp.gfactor*B/2.0),
                 oy-(int)(vp.gfactor*H/2.0),NULL);
    LineTo(hdc,ox+(int)(vp.gfactor*B/2.0),
               oy-(int)(vp.gfactor*H/2.0));
    LineTo(hdc,ox+(int)(vp.gfactor*B/2.0),
               oy-(int)(vp.gfactor*(H/2.0-tf)));
    LineTo(hdc,ox+(int)(vp.gfactor*tw/2.0),
               oy-(int)(vp.gfactor*(H/2.0-tf)));
    LineTo(hdc,ox+(int)(vp.gfactor*tw/2.0),
               oy+(int)(vp.gfactor*H/2.0));
    LineTo(hdc,ox-(int)(vp.gfactor*tw/2.0),
               oy+(int)(vp.gfactor*H/2.0));
    LineTo(hdc,ox-(int)(vp.gfactor*tw/2.0),
               oy-(int)(vp.gfactor*(H/2.0-tf)));
    LineTo(hdc,ox-(int)(vp.gfactor*B/2.0),
               oy-(int)(vp.gfactor*(H/2.0-tf)));
    LineTo(hdc,ox-(int)(vp.gfactor*B/2.0),
               oy-(int)(vp.gfactor*H/2.0));
  }
  else if(sect->sform.type==STEEL_TWEAK)
  {
    H =sect->sform.H;
    B =sect->sform.B;
    tw=sect->sform.tw;
    tf=sect->sform.tf;

    MoveToEx(hdc,ox-(int)(vp.gfactor*H/2.0),
                 oy-(int)(vp.gfactor*B/2.0),NULL);
    LineTo(hdc,ox-(int)(vp.gfactor*(H/2.0-tf)),
               oy-(int)(vp.gfactor*B/2.0));
    LineTo(hdc,ox-(int)(vp.gfactor*(H/2.0-tf)),
               oy-(int)(vp.gfactor*tw/2.0));
    LineTo(hdc,ox+(int)(vp.gfactor*H/2.0),
               oy-(int)(vp.gfactor*tw/2.0));
    LineTo(hdc,ox+(int)(vp.gfactor*H/2.0),
               oy+(int)(vp.gfactor*tw/2.0));
    LineTo(hdc,ox-(int)(vp.gfactor*(H/2.0-tf)),
               oy+(int)(vp.gfactor*tw/2.0));
    LineTo(hdc,ox-(int)(vp.gfactor*(H/2.0-tf)),
               oy+(int)(vp.gfactor*B/2.0));
    LineTo(hdc,ox-(int)(vp.gfactor*H/2.0),
               oy+(int)(vp.gfactor*B/2.0));
    LineTo(hdc,ox-(int)(vp.gfactor*H/2.0),
               oy-(int)(vp.gfactor*B/2.0));
  }
  else if(sect->sform.type==STEEL_RECTS)
  {
    for(i=0;i<sect->nsteel;i++)
    {
      drawmaterialrect(hdc,&(sect->srect[i]),ox,oy,vp);
    }
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
    drawstrandcircle(hdc,&(sect->strnd[i]),PEND,ox,oy,vp);
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
  fprintf(f,"As:ìSçú Ar:éÂãÿ Ac:ÉRÉìÉNÉäÅ[Ég Ap:ÇoÇbÉXÉgÉâÉìÉh\n");
  fprintf(f,"N:é≤óÕ Q:ÇπÇÒífóÕ ");
  fprintf(f,"Mt:ÇÀÇ∂ÇËÉÇÅ[ÉÅÉìÉg M:ã»Ç∞ÉÇÅ[ÉÅÉìÉg\n");
  fprintf(f,"ìYéö i:éní[ j:èIí[ c:íÜâõ\n");
  fprintf(f,"     a:ãñóe u:èIã«\n");
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
	/* 150429 fukushima for lst */
	if(n>=21)
	{
	  (sects+i)->ocode=strtol(*(data+20),NULL,10);
	}
	else
	{
	  (sects+i)->ocode=(sects+i)->code;
	}
	/***/

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
	/* 150429 fukushima for lst */
	if(n>=21)
	{
	  (sects+i)->ocode=strtol(*(data+20),NULL,10);
	}
	else
	{
	  (sects+i)->ocode=(sects+i)->code;
	}
	/***/

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
    s->N   =strtod(*(data+3),NULL);
    s->Q[0]=strtod(*(data+4),NULL);
    s->Q[1]=strtod(*(data+5),NULL);
    s->Mt  =strtod(*(data+6),NULL);
    s->M[0]=strtod(*(data+7),NULL);
    s->M[1]=strtod(*(data+8),NULL);
  }
  else if(nstr==7) /*TAIL*/
  {
    *inode=(int)strtol(*(data+0),NULL,10);
    s->N   =strtod(*(data+1),NULL);
    s->Q[0]=strtod(*(data+2),NULL);
    s->Q[1]=strtod(*(data+3),NULL);
    s->Mt  =strtod(*(data+4),NULL);
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

double wirelength(struct owire elem)
/*LENGTH OF WIRE ELEMENT.*/
{
  double dx,dy,dz,length;

  dx=(elem.node[1]->d[GX])-(elem.node[0]->d[GX]);
  dy=(elem.node[1]->d[GY])-(elem.node[0]->d[GY]);
  dz=(elem.node[1]->d[GZ])-(elem.node[0]->d[GZ]);

  length=sqrt(dx*dx+dy*dy+dz*dz);

  return length;
}/*wirelength*/

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

double slabwidth(struct element elem,struct element pair)
/*LENGTH OF SLAB ELEMENT.*/
{
  double dx1,dy1,dz1,dx2,dy2,dz2,length;

  dx1=(elem.node[0]->x)-(pair.node[0]->x);
  dy1=(elem.node[0]->y)-(pair.node[0]->y);
  dz1=(elem.node[0]->z)-(pair.node[0]->z);

  dx2=(elem.node[1]->x)-(pair.node[1]->x);
  dy2=(elem.node[1]->y)-(pair.node[1]->y);
  dz2=(elem.node[1]->z)-(pair.node[1]->z);

  length=0.5*(sqrt(dx1*dx1+dy1*dy1+dz1*dz1)
             +sqrt(dx2*dx2+dy2*dy2+dz2*dz2));

  return length;
}/*slabwidth*/

double slablength(struct element elem,struct element pair)
/*LENGTH OF SLAB ELEMENT.*/
{
  double dx1,dy1,dz1,dx2,dy2,dz2,length;

  dx1=(elem.node[1]->x)-(pair.node[0]->x);
  dy1=(elem.node[1]->y)-(pair.node[0]->y);
  dz1=(elem.node[1]->z)-(pair.node[0]->z);

  dx2=(pair.node[1]->x)-(elem.node[0]->x);
  dy2=(pair.node[1]->y)-(elem.node[0]->y);
  dz2=(pair.node[1]->z)-(elem.node[0]->z);

  length=0.5*(sqrt(dx1*dx1+dy1*dy1+dz1*dz1)
             +sqrt(dx2*dx2+dy2*dy2+dz2*dz2));

  return length;
}/*slablength*/

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

  if(Ag!=0.0) *Yg/=Ag;

  return;
}/*translatesection*/

double allowablebendingofrc(struct element elem,
                            struct materials m,
                            int axis,                     /*0:x 1:y*/
                            double Nd /*[kgf]*/,
                            int period)
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

  if(elem.sect->rein[0].F==3000.0)
  {
    if(period==PLONG) {rft=2000.0; rfc=-rft;}
    if(period==PSHORT){rft=3000.0; rfc=-rft;}
  }
  else if(elem.sect->rein[0].F==3500.0)
  {
    if(period==PLONG) {rft=2200.0; rfc=-rft;}
    if(period==PSHORT){rft=3500.0; rfc=-rft;}
  }

  ret=rft/m.rE;
  rec=rfc/m.rE;

  if(elem.sect->stype==STYPE_RC) cfc=/*m.cfc*/-elem.sect->crect[0].F/3.0;
  else return 0.0;

/*fprintf(fout0,"cfc=%.3f\n",cfc);*/

  if(period==PSHORT) cfc*=2.0;

  cec=cfc/m.cE;

  /*
  if(fout0!=NULL)
  {
    fprintf(fout0,"Ar=%.3f rYc=%.3f rYt=%.3f[cm2,cm]\n",Ar,rYc,rYt);
    fprintf(fout0,"Ac=%.3f cYc=%.3f cYt=%.3f[cm2,cm]\n",Ac,cYc,cYt);
  }
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
  if(fout0!=NULL)
  {
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
  }
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
      if(fout0!=NULL) fprintf(fout0,"I:NO ANSWER.Yn=%.5f\n",yn);
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
          /*if(fout0!=NULL)
          {
            fprintf(fout0,"CHOSEN ANSWER=%.5f\n",answer[ii]);
          }*/
          yn=answer[ii];
          found=1;
        }
      }
    }
    if(!found)
    {
      if(fout0!=NULL) fprintf(fout0,"II:NO ANSWER.\n");
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
          /*if(fout0!=NULL)
          {
            fprintf(fout0,"CHOSEN ANSWER=%.5f\n",answer[ii]);
          }*/
          yn=answer[ii];
          found=1;
        }
      }
    }
    if(!found)
    {
      if(fout0!=NULL) fprintf(fout0,"III:NO ANSWER.\n");
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
          /*if(fout0!=NULL)
          {
            fprintf(fout0,"CHOSEN ANSWER=%.5f\n",answer[ii]);
          }*/
          found=1;
        }
      }
    }
    if(!found)
    {
      if(fout0!=NULL) fprintf(fout0,"IV:NO ANSWER.\n");
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
    if(fout0!=NULL) fprintf(fout0,"OUT OF RANGE.\n");
    return 0.0;
  }*/

  /*if(fout0!=NULL)
  {
    fprintf(fout0,"RCé≤óÕ rcN=%11.3f[tf]",-Nd/1000.0);
    fprintf(fout0," èdêS:%11.3f",Yg);
    fprintf(fout0," íÜóßé≤îÕàÕ:%s íÜóßé≤=%11.3f",code,yn);
    fprintf(fout0," ãñóeíl rcMa=%11.3f[tfm]\n",Ma/100000.0);

    fprintf(fout0,"B=%.3f D=%.3f ",cBi[0],(cYi[0]-cYj[0]));
    fprintf(fout0,"N/bD=%.3f Ma/bD2=%.3f\n",
            Nd/cBi[0]/(cYi[0]-cYj[0]),
            Ma/cBi[0]/(cYi[0]-cYj[0])/(cYi[0]-cYj[0]));

    fprintf(fout0," %11.3f %11.3f 0.0\n",Ma/100000.0,-Nd/1000.0);
    fprintf(fout0,"%11.3f %11.3f 0.0",Ma/100000.0,-Nd/1000.0);
  }*/

  return Ma;
}/*allowablebendingofrc*/

double allowablebendingofsrc(struct element elem,
                             struct materials m,
                             int axis,                    /*0:x 1:y*/
                             double Nd /*[kgf]*/,
                             int period)
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

  if(elem.sect->stype==STYPE_S || elem.sect->stype==STYPE_SRC)
  {
    sft=/*m.sft*/elem.sect->srect[0].F;
    if(period==PLONG) sft/=1.5;

    sfc=-sft;
    set=sft/m.sE;
    sec=sfc/m.sE;
  }

  if(elem.sect->rein[0].F==3000.0)
  {
    if(period==PLONG) {rft=2000.0; rfc=-rft;}
    if(period==PSHORT){rft=3000.0; rfc=-rft;}
  }
  else if(elem.sect->rein[0].F==3500.0)
  {
    if(period==PLONG) {rft=2200.0; rfc=-rft;}
    if(period==PSHORT){rft=3500.0; rfc=-rft;}
  }

  if(elem.sect->stype==STYPE_RC || elem.sect->stype==STYPE_SRC)
  {
    ret=rft/m.rE;
    rec=rfc/m.rE;
  }

  if(elem.sect->stype==STYPE_RC)
  {
    cfc=/*m.cfc*/-elem.sect->crect[0].F/3.0;
  }
  else if(elem.sect->stype==STYPE_SRC)
  {
    cfc=/*m.cfc*/-elem.sect->crect[0].F/3.0*(1.0-15.0*(As/2.0/Ac)); /*"SRC STANDARD"(29)*/
  }

  if(elem.sect->stype==STYPE_RC || elem.sect->stype==STYPE_SRC)
  {
    if(period==PSHORT) cfc*=2.0;
    cec=cfc/m.cE;
  }
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
  if(fout0!=NULL)
  {
    fprintf(fout0,"As=%.3f sYc=%.3f sYt=%.3f[cm2,cm]\n",As,sYc,sYt);
    fprintf(fout0,"Ar=%.3f rYc=%.3f rYt=%.3f[cm2,cm]\n",Ar,rYc,rYt);
    fprintf(fout0,"Ac=%.3f cYc=%.3f cYt=%.3f[cm2,cm]\n",Ac,cYc,cYt);
  }
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

  boundaryvalue(elem,m,Yn[15],sBi,sDi,cBi,cDi,
                              sYi,sYj,cYi,cYj,
                              sfi,sfj,cfi,cfj,
                              rAi,rYi,rfi,
                              &N[15],&M[15]);

  /*if(fout0!=NULL)
  {
    for(i=1;i<=15;i++)
    {
      fprintf(fout0,"BOUND%2d:Yn=%11.3f N=%12.3f M=%12.3f\n",
              i,Yn[i],N[i],M[i]);
    }
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
        /*if(fout0!=NULL)
        {
          fprintf(fout0,"CHOSEN ANSWER=%.5f\n",answer[ii]);
        }*/
        yn=answer[ii];
        found=1;
      }
    }
    if(!found)
    {
      if(fout0!=NULL) fprintf(fout0,"I:NO ANSWER.\n");
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
          /*if(fout0!=NULL)
          {
            fprintf(fout0,"CHOSEN ANSWER=%.5f\n",answer[ii]);
          }*/
          yn=answer[ii];
          found=1;
        }
      }
    }
    if(!found)
    {
      if(fout0!=NULL) fprintf(fout0,"II:NO ANSWER.\n");
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
          /*if(fout0!=NULL)
          {
            fprintf(fout0,"CHOSEN ANSWER=%.5f\n",answer[ii]);
          }*/
          yn=answer[ii];
          found=1;
        }
      }
    }
    if(!found)
    {
      if(fout0!=NULL) fprintf(fout0,"III:NO ANSWER.\n");
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
          /*if(fout0!=NULL)
          {
            fprintf(fout0,"CHOSEN ANSWER=%.5f\n",answer[ii]);
          }*/
          yn=answer[ii];
          found=1;
        }
      }
    }
    if(!found)
    {
      if(fout0!=NULL) fprintf(fout0,"IV:NO ANSWER.\n");
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
          /*if(fout0!=NULL)
          {
            fprintf(fout0,"CHOSEN ANSWER=%.5f\n",answer[ii]);
          }*/
          yn=answer[ii];
          found=1;
        }
      }
    }
    if(!found)
    {
      if(fout0!=NULL) fprintf(fout0,"V:NO ANSWER.\n");
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

  /*if(fout0!=NULL)
  {
    fprintf(fout0,"SRCé≤óÕ srcN=%11.3f[tf]",-Nd/1000.0);
    fprintf(fout0," èdêS:%11.3f",Yg);
    fprintf(fout0," íÜóßé≤îÕàÕ:%s íÜóßé≤=%11.3f",code,yn);
    fprintf(fout0," ãñóeíl srcMa=%11.3f[tfm]\n",Ma/100000.0);

    fprintf(fout0," %11.3f %11.3f 0.0\n",Ma/100000.0,-Nd/1000.0);
    fprintf(fout0,"%11.3f %11.3f 0.0",Ma/100000.0,-Nd/1000.0);
  }*/

  return Ma;
}/*allowablebendingofsrc*/

int ultimateaxialforceofsrc(struct element elem,
                            double *Nmax,double *Nmin,
                            double *sN,double *sM,
                            double *rN,double *rM,
                            double *cN,double *cM)
/*RETURN:ULTIMATE AXIAL FORCE OF RC,SRC COLUMN.*/
/*N=+:TENSION -:COMPRESSION*/
{
  /*char str[1024],non[256];*/
  int i;
  double sBi[MAXSRECT],sDi[MAXSRECT];       /*WIDTH,DEPTH OF STEELS*/
  double cBi[MAXCRECT],cDi[MAXCRECT];    /*WIDTH,DEPTH OF CONCRETES*/
  double rAi[MAXREINS];
  double Yn;                 /*CENTER OF GRAVITY,NEUTRAL AXIS*/
  double sYi[MAXSRECT],sYj[MAXSRECT];         /*UPPER,LOWER OF RECT*/
  double rYi[MAXREINS];
  double cYi[MAXCRECT],cYj[MAXCRECT];
  int nbound;                               /*SUPPLEMENTARY BOUNDS.*/
  double bounds[2*MAXSRECT+MAXREINS+2*MAXCRECT];
  double sYt,sYc,rYt,rYc,cYt,cYc;

  double sftu,sfcu,rftu,rfcu,cfcu;              /*ULTIMATE STRESSES*/
  double N[2*MAXSRECT+MAXREINS+2*MAXCRECT];         /*TENSION*/
  double M[2*MAXSRECT+MAXREINS+2*MAXCRECT];

  double As,Ac,Ar;
  double sAY,rAY;
  double reinrate;

  nbound=0;

  As=0.0; /*STEEL AREA TOTAL[cm2]*/
  sAY=0.0;
  sYc=-1000.0; sYt=1000.0; /*[cm]*/
  for(i=0;i<(elem.sect->nsteel);i++)
  {
    sBi[i]=fabs(elem.sect->srect[i].right
               -elem.sect->srect[i].left);
    sYi[i]=elem.sect->srect[i].top;
    sYj[i]=elem.sect->srect[i].bottom;

    sDi[i]=sYi[i]-sYj[i];

    As+=(sBi[i]*sDi[i]); /*[cm2]*/

    sAY+=(sBi[i]*sDi[i])*(sYi[i]+sYj[i]); /*[cm3]*/

    addnewbound(bounds,&nbound,sYi[i]);
    addnewbound(bounds,&nbound,sYj[i]);

    if(sYc<sYi[i]) sYc=sYi[i];
    if(sYt>sYj[i]) sYt=sYj[i];
  }

  Ac=0.0; /*CONCRETE AREA TOTAL[cm2]*/
  cYc=-1000.0; cYt=1000.0; /*[cm]*/
  for(i=0;i<(elem.sect->nconc);i++)
  {
    cBi[i]=fabs(elem.sect->crect[i].right
               -elem.sect->crect[i].left);
    cYi[i]=elem.sect->crect[i].top;
    cYj[i]=elem.sect->crect[i].bottom;

    cDi[i]=cYi[i]-cYj[i];

    Ac+=(cBi[i]*cDi[i]); /*[cm2]*/

    addnewbound(bounds,&nbound,cYi[i]);
    addnewbound(bounds,&nbound,cYj[i]);

    if(cYc<cYi[i]) cYc=cYi[i];
    if(cYt>cYj[i]) cYt=cYj[i];
  }

  Ar=0.0;
  rAY=0.0;
  rYc=-1000.0; rYt=1000.0; /*[cm]*/
  for(i=0;i<(elem.sect->nrein);i++)
  {
    rAi[i]=elem.sect->rein[i].area;

    rYi[i]=elem.sect->rein[i].y;

    Ar +=rAi[i];
    rAY+=rAi[i]*rYi[i];

    if(addnewbound(bounds,&nbound,rYi[i]))
    {
      bounds[nbound]=rYi[i]; /*ADD SAME BOUND*/
      nbound++;
    }

    if(rYc<rYi[i]) rYc=rYi[i];
    if(rYt>rYi[i]) rYt=rYi[i];
  }

  sortdouble(bounds,nbound); /*SORT BOUNDARIES.*/

  /*MATERIALS*/
  if(elem.sect->stype==STYPE_S || elem.sect->stype==STYPE_SRC)
  {
    sftu=/*m.sftu*/elem.sect->srect[0].F;
    sfcu=-sftu;
  }

  if(elem.sect->stype==STYPE_RC || elem.sect->stype==STYPE_SRC)
  {
    rftu=elem.sect->rein[0].F;
    rfcu=-rftu;
  }

  if(elem.sect->stype==STYPE_RC)
  {
    cfcu=-0.85*(elem.sect->crect[0].F);
  }
  else if(elem.sect->stype==STYPE_SRC)
  {
    cfcu=-(elem.sect->crect[0].F)*(0.85-2.5*(As/2.0/Ac)); /*"SRC STANDARD"*/
  }

  /*sprintf(str,"sftu=%.3f sfcu=%.3f rftu=%.3f rfcu=%.3f cfcu=%.3f",
          sftu,sfcu,rftu,rfcu,cfcu);
  MessageBox(NULL,str,"Ultimate Axial Force",MB_OK);*/

  for(i=0;i<nbound;i++) /*BOUNDARY DEFINITION.*/
  {
    Yn=bounds[i];

    if(i<(nbound-1) && bounds[i]==bounds[i+1]) reinrate=1.0;
    else if(i>=1 && bounds[i-1]==bounds[i])    reinrate=0.0;
    else reinrate=0.0;

    boundaryvalueultimate(elem,Yn,sBi,sDi,cBi,cDi,
                                  sYi,sYj,cYi,cYj,
                                  sftu,sfcu,
                                  rftu,rfcu,
                                  0.0,cfcu,
                                  rAi,rYi,
                                  sN,sM,rN,rM,cN,cM,
                                  &N[i],&M[i],
                                  reinrate);
  }

  /*sprintf(str,"\0");
  for(i=0;i<nbound;i++)
  {
    sprintf(non,"BOUND%2d:Yn=%11.3f N=%12.3f M=%12.3f\n",
            i+1,bounds[i],N[i],M[i]);
    strcat(str,non);
  }
  MessageBox(NULL,str,"Ultimate Axial Force",MB_OK);*/

  /*RANGE DETERMINATION.*/
  *Nmax=N[0];
  *Nmin=N[nbound-1];

  return 1;
}/*ultimateaxialforceofsrc*/

double ultimatebendingofsteel(struct section sect,
                              double E,double F,
                              double lkxx,double lkyy, /*BUCKLING LENGTH[cm]*/
                              int axis,                   /*0:x 1:y*/
                              double Nd)                    /*[kgf]*/
/*RETURN:ULTIMATE BENDING OF STEEL.*/
/*N=+:TENSION -:COMPRESSION*/
/*STEEL TYPE:HKYOU BY RECTS,HKYOU,HWEAK*/
{
  double A,Ixx,Iyy,Zxx,Zyy,ixx,iyy,i;
  double H,B,tw,tf;
  double lk,LAM,lam,lamxx,lamyy,nyu,sfca,sftu,sfn;
  double An,At,Ac,sign,Mu=0.0;

  /*ALLOWABLE COMPRESSION STRESS*/
  LAM=sqrt(PI*PI*E/0.6/F);                      /*"S STANDARD"(5.5)*/

  steelcoefficients(sect,&A,&Ixx,&Iyy,&Zxx,&Zyy,&ixx,&iyy); /*[cm]*/

  /*if(ixx<=iyy) i=ixx;*/
  /*else         i=iyy;*/
  /*lam=lk/i;*/                                /*"S STANDARD"(11.1)*/

  if(sect.bblength[0]!=0.0) lk=sect.bblength[0];
  else if(sect.bbfact[0]!=0.0) lk=lkxx*sect.bbfact[0];
  else lk=lkxx;

  lamxx=lk/ixx;

  if(sect.bblength[1]!=0.0) lk=sect.bblength[1];
  else if(sect.bbfact[1]!=0.0) lk=lkyy*sect.bbfact[1];
  else lk=lkyy;

  lamyy=lk/iyy;

  if(lamxx>lamyy) lam=lamxx;
  else            lam=lamyy;

  if((sect.etype==COLUMN && lam>200)||(lam>250))
  {
    if(fout0!=NULL)
    {
      fprintf(fout0,"ELEMENT TOO LONG.Lk/i=%.5f\n",lam);
    }
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

  sfca*=1.5; /*SHORT*/

  /*ULTIMATE TENSION STRESS*/
  if(sect.sform.type==STEEL_RECTS) sftu=sect.srect[0].F;
  else                             sftu=sect.sform.F;

  /*STEEL TYPE*/
  if(sect.sform.type==STEEL_RECTS) /*FOR ONLY HKYOU INPUT BY RECTS.*/
  {
    H =(sect.srect[0].top-sect.srect[2].bottom);             /*[cm]*/
    B =(sect.srect[0].right-sect.srect[0].left);             /*[cm]*/
    tw=(sect.srect[1].right-sect.srect[1].left);             /*[cm]*/
    tf=(sect.srect[0].top-sect.srect[0].bottom);             /*[cm]*/
  }
  else
  {
    H =sect.sform.H;  /*[cm]*/
    B =sect.sform.B;  /*[cm]*/
    tw=sect.sform.tw; /*[cm]*/
    tf=sect.sform.tf; /*[cm]*/
  }

  /*AREA FOR Nd*/
  if(Nd>=0.0) sfn=sftu; /*TENSION*/
  if(Nd<0.0)  sfn=sfca; /*COMPRESSION*/

  An=Nd/sfn;

  /*AREA COMPRESSED BY BENDING*/
  Ac=(sftu)/(sftu+sfca)*(A-An);

  /*AREA TENSED BY BENDING*/
  At=(sfca)/(sftu+sfca)*(A-An);

  /*ULTIMATE BENDING*/
  if(Nd>=0.0) sign= 1.0;
  else        sign=-1.0;

  if((axis==SX && sect.sform.type==STEEL_RECTS) ||
     (axis==SX && sect.sform.type==STEEL_HKYOU) ||
     (axis==SY && sect.sform.type==STEEL_HWEAK))
  {
    if(At>=B*tf)
    {
      Mu=sftu*(B*tf)*(H/2-tf/2)
        +sftu*(At-B*tf)*(H/2-tf-(At-B*tf)/tw/2)
        +sfca*(B*tf)*(H/2-tf/2)
        +sfca*(Ac-B*tf)*(H/2-tf-(Ac-B*tf)/tw/2)
        +sign*sfn*(Ac-At)*(H/2-tf-(At-B*tf)/tw-(Ac-At)/tw/2);
    }
    else if(Ac>=B*tf)
    {
      Mu=sftu*At*(H/2-At/B/2)
        +sfca*(B*tf)*(H/2-tf/2)
        +sfca*(Ac-B*tf)*(H/2-tf-(Ac-B*tf)/tw/2)
        +sign*sfn*(B*tf-At)*(H/2-At/B/2)
        +sign*sfn*(Ac-B*tf)*(H/2-tf/2-(Ac-B*tf)/tw/2);
    }
    else
    {
      Mu=sftu*At*(H/2-At/B/2)
        +sfca*Ac*(H/2-Ac/B/2)
        +sign*sfn*(Ac-At)*(H/2-At/B-((Ac-At)/B)/2);
    }
  }
  if((axis==SY && sect.sform.type==STEEL_RECTS) ||
     (axis==SY && sect.sform.type==STEEL_HKYOU) ||
     (axis==SX && sect.sform.type==STEEL_HWEAK))
  {
    if(At>=(2*tf*(B-tw)/2))
    {
      Mu=sftu*(2*tf*(B-tw)/2)*(B/2-(B-tw)/4)
        +sftu*(At-2*tf*(B-tw)/2)*(tw/2-(At-2*tf*(B-tw)/2)/H/2)
        +sfca*(2*tf*(B-tw)/2)*(B/2-(B-tw)/4)
        +sfca*(Ac-2*tf*(B-tw)/2)*(tw/2-(Ac-2.0*tf*(B-tw)/2)/H/2)
        +sign*sfn*(Ac-At)*(tw/2-(At-2*tf*(B-tw)/2)/H-(Ac-At)/H/2);
    }
    else if(Ac>=(2*tf*(B-tw)/2))
    {
      Mu=sftu*At*(B/10/2-At/(2*tf)*10/2)
        +sfca*(2*tf*(B-tw)/2)*(B/2-(B-tw)/4)
        +sfca*(Ac-2*tf*(B-tw)/2)*(tw/2-(Ac-2*tf*(B-tw)/2)/H/2)
        +sign*sfn*(2*tf*(B-tw)/2-At)*(tw/2+(2*tf*(B-tw)/2-At)/(2*tf)/2)
        +sign*sfn*(Ac-2*tf*(B-tw)/2)*(tw/2-(Ac-2*tf*(B-tw)/2)/H/2);
    }
    else
    {
      Mu=sftu*At*(B/2-At/(2.0*tf)/2)
        +sfca*Ac*(B/2-Ac/(2.0*tf)/2)
        +sign*sfn*(Ac-At)*(B/2-At/(2.0*tf)-(Ac-At)/(2.0*tf)/2);
    }
  }

  return Mu;
}/*ultimatebendingofsteel*/

double allowablebendingofsteel(int period,
                               struct section sect,
                               double E,
                               double lk,      /*BUCKLING LENGTH[cm]*/
                               int axis,                   /*0:x 1:y*/
                               double Nd)                    /*[kgf]*/
/*RETURN:ULTIMATE BENDING OF STEEL.*/
/*N=+:TENSION -:COMPRESSION*/
/*STEEL TYPE:HKYOU BY RECTS,HKYOU,HWEAK*/
{
  double A,Af,Ixx,Iyy,Zxx,Zyy,ixx,iyy;
  double H,B,tw,tf;
  double sfta,sfba,sF;
  double Na=0.0,Ma=0.0;

  /*ALLOWABLE AXIAL FORCE*/
  if(sect.sform.type==STEEL_RECTS) sF=sect.srect[0].F; /*[kgf/cm2]*/
  else                             sF=sect.sform.F;

if(!strcmp(prj,"athens") && sect.code==201) lk=480.0; /*ATHENS*/
if(!strcmp(prj,"athens") && sect.code==301) lk=400.0; /*ATHENS*/
if(!strcmp(prj,"athens") && sect.code==401) lk=400.0; /*ATHENS*/

if(!strcmp(prj,"aii") && sect.code==211) lk*=1.2/sqrt(2.0);
if(!strcmp(prj,"aii") && sect.code==212) lk*=1.2/sqrt(2.0);
if(!strcmp(prj,"aii") && sect.code==213) lk*=1.2/sqrt(2.0);

if(!strcmp(prj,"naka") && sect.code==202) lk=65.0;
if(!strcmp(prj,"naka") && sect.code==203) lk=65.0;
if(!strcmp(prj,"naka") && sect.code==204) lk=80.0;
if(!strcmp(prj,"naka") && sect.code==205) lk=80.0;

if(!strcmp(prj,"hakoii") && sect.code==501) lk*=0.3;
if(!strcmp(prj,"hakoii") && sect.code==502) lk*=0.3;
if(!strcmp(prj,"hakoii") && sect.code==503) lk*=0.3;

/*if(!strcmp(prj,"yagi") && 201<=sect.code && sect.code<=210) lk=294.0;*/
if(!strcmp(prj,"yagi") && sect.code==201) lk=250.0; /*r=7200*/
if(!strcmp(prj,"yagi") && sect.code==202) lk=225.0; /*r=4200*/
if(!strcmp(prj,"yagi") && sect.code==203) lk=210.0; /*r=3000*/
if(!strcmp(prj,"yagi") && sect.code==204) lk=215.0; /*r=3600*/
if(!strcmp(prj,"yagi") && sect.code==205) lk=294.0; /*r=1700+RIB*/
if(!strcmp(prj,"yagi") && sect.code==206) lk=240.0; /*r=5500*/
if(!strcmp(prj,"yagi") && sect.code==207) lk=205.0; /*r=2900*/
if(!strcmp(prj,"yagi") && sect.code==208) lk=200.0; /*r=2600*/
if(!strcmp(prj,"yagi") && sect.code==209) lk=294.0; /*r=    */
if(!strcmp(prj,"yagi") && sect.code==210) lk=294.0; /*r=    */
if(!strcmp(prj,"yagi") && sect.code==211) lk=294.0; /*r=    */

if(!strcmp(prj,"iwami") && sect.code==201) lk*=1.2/1.9;
if(!strcmp(prj,"iwami") && sect.code==202) lk*=1.2/2.5;

if(!strcmp(prj,"tep") && sect.code==204) lk=220.0;

/*if(!strcmp(prj,"vene") && sect.code==201) lk*=2.0*sqrt(2.0);*/

  if(Nd<0.0)  Na=-allowablecompressionofsteel(E,sF,lk,lk,sect);
  if(Nd>=0.0) Na=-allowabletensionofsteel(sF,sect);

  if(period==PSHORT) Na*=1.5;

  if(Na==0.0 || (Nd/Na)>=1.0)
  {
    fprintf(fout0,"ERROR.N > Na=%.3f\n",Na);
    return 0.0;
  }

  /*STEEL TYPE*/
  if(sect.sform.type==STEEL_RECTS) /*FOR ONLY HKYOU INPUT BY RECTS.*/
  {
    H =(sect.srect[0].top-sect.srect[2].bottom);             /*[cm]*/
    B =(sect.srect[0].right-sect.srect[0].left);             /*[cm]*/
    tw=(sect.srect[1].right-sect.srect[1].left);             /*[cm]*/
    tf=(sect.srect[0].top-sect.srect[0].bottom);             /*[cm]*/
  }
  else
  {
    H =sect.sform.H;  /*[cm]*/
    B =sect.sform.B;  /*[cm]*/
    tw=sect.sform.tw; /*[cm]*/
    tf=sect.sform.tf; /*[cm]*/
  }

  /*COEFFICIENTS*/
  steelcoefficients(sect,&A,&Ixx,&Iyy,&Zxx,&Zyy,&ixx,&iyy);   /*[cm]*/

  /*ALLOWABLE STRESS*/
  sfta=sF;
  if(period==PLONG) sfta/=1.5;

  if((axis==SX && sect.sform.type==STEEL_RECTS) ||
     (axis==SX && sect.sform.type==STEEL_HKYOU) ||
     (axis==SX && sect.sform.type==STEEL_ANGLE) ||
     (axis==SX && sect.sform.type==STEEL_TKYOU))
  {
    Af=B*tf;
    if(sect.bblength[0]!=0.0) lk=sect.bblength[0];   /*****/
    else if(sect.bbfact[0]!=0.0) lk*=sect.bbfact[0]; /*****/

    sfba=900.0*Af/lk/H*1000.0;          /*[kgf/cm2] "S STANDARD"(5.8)*/
    if(period==PSHORT) sfba*=1.5;
    if(sfba>sfta) sfba=sfta;

    Ma=(1.0-Nd/Na)*sfba*Zxx;                                /*[kgfcm]*/
  }
  if((axis==SY && sect.sform.type==STEEL_HWEAK) ||
     (axis==SY && sect.sform.type==STEEL_TWEAK))
  {
    Af=B*tf;
    if(sect.bblength[1]!=0.0) lk=sect.bblength[1];   /*****/
    else if(sect.bbfact[1]!=0.0) lk*=sect.bbfact[1]; /*****/

    sfba=900.0*Af/lk/H*1000.0;          /*[kgf/cm2] "S STANDARD"(5.8)*/
    if(period==PSHORT) sfba*=1.5;
    if(sfba>sfta) sfba=sfta;

    Ma=(1.0-Nd/Na)*sfba*Zyy;                                /*[kgfcm]*/
  }
  if((axis==SY && sect.sform.type==STEEL_ANGLE))
  {
    Af=H*tw;

    sfba=900.0*Af/lk/H*1000.0;          /*[kgf/cm2] "S STANDARD"(5.8)*/
    if(period==PSHORT) sfba*=1.5;
    if(sfba>sfta) sfba=sfta;

    Ma=(1.0-Nd/Na)*sfba*Zyy;                                /*[kgfcm]*/
  }
  if((axis==SY && sect.sform.type==STEEL_RECTS) ||
     (axis==SY && sect.sform.type==STEEL_HKYOU) ||
     (axis==SY && sect.sform.type==STEEL_RPIPE) ||
     (axis==SY && sect.sform.type==STEEL_CPIPE) ||
     (axis==SY && sect.sform.type==STEEL_TKYOU))
  {
    Ma=(1.0-Nd/Na)*sfta*Zyy;                                /*[kgfcm]*/
  }
  if((axis==SX && sect.sform.type==STEEL_HWEAK) ||
     (axis==SX && sect.sform.type==STEEL_RPIPE) ||
     (axis==SX && sect.sform.type==STEEL_CPIPE) ||
     (axis==SX && sect.sform.type==STEEL_TWEAK))
  {
    Ma=(1.0-Nd/Na)*sfta*Zxx;                                /*[kgfcm]*/
  }

  /*fprintf(fout0,"Ma=%.3f\n",Ma);*/

  return Ma;
}/*allowablebendingofsteel*/

double allowablebendingofalumi(int period,
                               struct section sect,
                               double E,
                               double lk,      /*BUCKLING LENGTH[cm]*/
                               int axis,                   /*0:x 1:y*/
                               double Nd)                    /*[kgf]*/
/*RETURN:ULTIMATE BENDING OF ALUMI.*/
/*N=+:TENSION -:COMPRESSION*/
/*STEEL TYPE:HKYOU BY RECTS,HKYOU,HWEAK*/
{
  double A,Af,Ixx,Iyy,Jz,Cw,Zxx,Zyy,ixx,iyy,G;
  double H,B,tw,tf;
  double alfta,alfba,alF;
  double fb1,fb2,fb3;
  double Me,My,M2M1,Cb;
  double ramda,pramda,eramda,nyu;
  double Na=0.0,Ma=0.0;

  /*ALLOWABLE AXIAL FORCE*/
  alF=sect.sform.F;
  G=E/2/(1.0+0.3);

/*DEFINE BUCKLING LENGTH*/

if(!strcmp(prj,"Table") && sect.code==201 && axis==SY) lk=480;
if(!strcmp(prj,"Table") && sect.code==202 && axis==SY) lk=360;
if(!strcmp(prj,"Table") && sect.code==501 && axis==SX) lk=360;

if(!strcmp(prj,"Alumi") && sect.code==502) lk=455.0;

  if(Nd<0.0)   Na=-allowablecompressionofalumi(E,alF,lk,lk,sect,axis);
  if(Nd>=0.0)  Na=-allowabletensionofalumi(alF,sect);

  if(period==PSHORT) Na*=1.5;

  /*ALUMI TYPE*/
  H =sect.sform.H;  /*[cm]*/
  B =sect.sform.B;  /*[cm]*/
  tw=sect.sform.tw; /*[cm]*/
  tf=sect.sform.tf; /*[cm]*/

  /*COEFFICIENTS*/
  alumicoefficients(sect,&A,&Ixx,&Iyy,&Jz,&Cw,&Zxx,&Zyy,&ixx,&iyy);   /*[cm]*/

  /*ALLOWABLE STRESS*/
  alfta=alF;
  if(period==PLONG) alfta/=1.5;

///////////////////////////////////////////////////////////////2008.01.17 TAKATA
  if(axis==SX)
  {
    if(sect.bblength[0]!=0.0) lk=sect.bblength[0];
    else if(sect.bbfact[0]!=0.0) lk*=sect.bbfact[0];
  }

  else if(axis==SY)
  {
    if(sect.bblength[1]!=0.0) lk=sect.bblength[1];
    else if(sect.bbfact[1]!=0.0) lk*=sect.bbfact[1];
  }

  if(Na==0.0 || (Nd/Na)>=1.0)
  {
    fprintf(fout0,"ERROR.N > Na=%.3f[tf] (Lk=%.1f[cm])\n",Na/1000.0,lk);
    return 0.0;
  }
  else
    fprintf(fout0,"Na=%.3f[tf] (Lk=%.1f[cm])\n",Na/1000.0,lk);
////////////////////////////////////////////////////////////////////////////////

  M2M1=1.0;
  Cb=1.75+1.05*M2M1+0.3*(M2M1)*(M2M1);
  eramda=1/sqrt(0.5);
  pramda=0.6+0.3*M2M1;
  if(sect.sform.type==STEEL_HKYOU || sect.sform.type==STEEL_TWEAK ||
     sect.sform.type==STEEL_TKYOU || sect.sform.type==STEEL_TWEAK ||
     sect.sform.type==STEEL_CKYOU || sect.sform.type==STEEL_CWEAK ||
     sect.sform.type==STEEL_ANGLE)
  {
    if(Ixx>Iyy) Me=Cb*sqrt(PI*PI*PI*PI*E*Iyy*E*Cw/lk/lk/lk/lk+PI*PI*E*Iyy*G*Jz/lk/lk);
    else        Me=Cb*sqrt(PI*PI*PI*PI*E*Ixx*E*Cw/lk/lk/lk/lk+PI*PI*E*Ixx*G*Jz/lk/lk);
  }
  else
  {
    if(Ixx>Iyy) Me=Cb*sqrt(PI*PI*E*Iyy*G*Jz/lk/lk);
    else        Me=Cb*sqrt(PI*PI*E*Ixx*G*Jz/lk/lk);
  }


  if((axis==SX && sect.sform.type==STEEL_HKYOU) ||
     (axis==SX && sect.sform.type==STEEL_ANGLE) ||
     (axis==SX && sect.sform.type==STEEL_TKYOU) ||
     (axis==SX && sect.sform.type==STEEL_CKYOU))
  {
    My=alF*Zxx;
    ramda=sqrt(My/Me);
    nyu=1.5+(ramda/eramda)*(ramda/eramda)/1.5;
    if(nyu>2.17) nyu=2.17;

    if(ramda<=pramda)               /*[kgf/cm2] "ALUMI STANDARD"(4.3.1)*/
    {
      alfba=alF/nyu;
    }
    else if(pramda<ramda && ramda<=eramda)
    {
      alfba=(1.0+0.5*(ramda-pramda)/(eramda-pramda))*alF/nyu;
    }
    else
    {
      alfba=alF/ramda/ramda/nyu;
    }

    if(period==PSHORT) alfba*=1.5;
    if(alfba>alfta) alfba=alfta;

    Ma=(1.0-Nd/Na)*alfba*Zxx;                                /*[kgfcm]*/
  }
  if((axis==SY && sect.sform.type==STEEL_HWEAK) ||
     (axis==SY && sect.sform.type==STEEL_TWEAK) ||
     (axis==SY && sect.sform.type==STEEL_CWEAK) ||
     (axis==SY && sect.sform.type==STEEL_ANGLE))
  {
    My=alF*Zyy;
    ramda=sqrt(My/Me);
    nyu=1.5+(ramda/eramda)*(ramda/eramda)/1.5;
    if(nyu>2.17) nyu=2.17;

    if(ramda<=pramda)               /*[kgf/cm2] "ALUMI STANDARD"(4.3.1)*/
    {
      alfba=alF/nyu;
    }
    else if(pramda<ramda && ramda<=eramda)
    {
      alfba=(1.0+0.5*(ramda-pramda)/(eramda-pramda))*alF/nyu;
    }
    else
    {
      alfba=alF/ramda/ramda/nyu;
    }

    if(period==PSHORT) alfba*=1.5;
    if(alfba>alfta) alfba=alfta;

    Ma=(1.0-Nd/Na)*alfba*Zyy;                                /*[kgfcm]*/
  }
  if((axis==SY && sect.sform.type==STEEL_HKYOU) ||
     (axis==SY && sect.sform.type==STEEL_RPIPE) ||
     (axis==SY && sect.sform.type==STEEL_CPIPE) ||
     (axis==SY && sect.sform.type==STEEL_TKYOU) ||
     (axis==SY && sect.sform.type==STEEL_CKYOU))
  {
    Ma=(1.0-Nd/Na)*alfta*Zyy;                                /*[kgfcm]*/
  }
  if((axis==SX && sect.sform.type==STEEL_HWEAK) ||
     (axis==SX && sect.sform.type==STEEL_RPIPE) ||
     (axis==SX && sect.sform.type==STEEL_CPIPE) ||
     (axis==SX && sect.sform.type==STEEL_TWEAK) ||
     (axis==SX && sect.sform.type==STEEL_CWEAK))
  {
    Ma=(1.0-Nd/Na)*alfta*Zxx;                                /*[kgfcm]*/
  }

  return Ma;
}/*allowablebendingofalumi*/

double allowablestressofflatbar(int period,
                                struct section sect,
                                double E,
                                double L,       /*ELEMENT LENGTH[cm]*/
                                struct stress *st,        /*[tf,tfm]*/
                                int iend,
                                double *Ncr,
                                double *Qcrx,double *Qcry,
                                double *Mcrx,double *Mcry)
/*FLAT BAR COLUMN.*/
/*BUCKLING LENGTH,FACTOR FROM LST FILE.*/
{
  double poi,G,Fc,Ft,Fb,Fs/*,fu*/;
  double B,D,fact,Lk,A,Ix,Iy,Iw,Jz,Zx,Zy,Zpx,Zpy,ix,iy;
  double lamda,Lamda,nyu,fc,Na,Nax,Nay,Np,Ncrx,Ncry,Qax,Qay,Myx,Myy,Mpx,Mpy;
  double g,k,ku,kb/*,q*/,M0M1,M2M1,gzai,C1,C2,C3,Mcre,Mcrex,Mcrey,Mcrp;
  double Kx,Ky,cx,cy,alphax,alphay;
  double N,Qx,Qy,Mx,My,Qo,Mo;
  double rate1,rate2,rate;
  double Lm; /*MEMORY LENGTH*/
  double bblength[2],bbfact[2],btlength[2],btfact[2]; /*BUCKLING LENGTH*/

if(!strcmp(prj,"karu") && sect.code==203) L=12.0;
if(!strcmp(prj,"karu") && sect.code==501) L=12.0;

  /*INITIAL*/
  Lm=L;
  *Ncr =0.0;
  *Qcrx=0.0;
  *Qcry=0.0;
  *Mcrx=0.0;
  *Mcry=0.0;

  /*MATERIAL*/
  if(sect.stype==STYPE_S)
  {
    poi=1.0/3.0;
    Fc=sect.sform.F;
    Ft=sect.sform.F;
    Fb=sect.sform.F;
    Fs=sect.sform.F;
    /*
    if(F==2400.0)      fu=4000.0;
    else if(F==3300.0) fu=5000.0;
    else               fu=4000.0;
    */
  }
  if(sect.stype==STYPE_GLASS)
  {
    poi=0.225;
    Fc=sect.wform.Fc;
    Ft=sect.wform.Ft;
    Fb=sect.wform.Fb;
    Fs=sect.wform.Fs;
  }
  if(sect.stype==STYPE_ACRYL)
  {
    poi=0.400;
    Fc=sect.wform.Fc;
    Ft=sect.wform.Ft;
    Fb=sect.wform.Fb;
    Fs=sect.wform.Fs;
  }

  G=E/2.0/(1.0+poi);

  /*FORM,STRESS*/
  if(sect.sform.H >= sect.sform.B)
  {
    B=sect.sform.B; /*[cm]*/
    D=sect.sform.H;

    if(iend==HEAD) N=st->N;
    if(iend==TAIL) N=-(st->N);
    Qx=fabs(st->Q[0]);
    Qy=fabs(st->Q[1]);
    Mx=fabs(st->M[0]);
    My=fabs(st->M[1]);

    bblength[0]=sect.bblength[0];
    bblength[1]=sect.bblength[1];
    bbfact[0]=sect.bbfact[0];
    bbfact[1]=sect.bbfact[1];
    btlength[0]=sect.btlength[0];
    btlength[1]=sect.btlength[1];
    btfact[0]=sect.btfact[0];
    btfact[1]=sect.btfact[1];
  }
  else
  {
    B=sect.sform.H;
    D=sect.sform.B;

    if(iend==HEAD) N=st->N;
    if(iend==TAIL) N=-(st->N);
    Qx=fabs(st->Q[1]);
    Qy=fabs(st->Q[0]);
    Mx=fabs(st->M[1]);
    My=fabs(st->M[0]);

    bblength[0]=sect.bblength[1];
    bblength[1]=sect.bblength[0];
    bbfact[0]=sect.bbfact[1];
    bbfact[1]=sect.bbfact[0];
    btlength[0]=sect.btlength[1];
    btlength[1]=sect.btlength[0];
    btfact[0]=sect.btfact[1];
    btfact[1]=sect.btfact[0];
  }
  /*
  fprintf(fout0,"B =%8.3f\n",B);
  fprintf(fout0,"D =%8.3f\n",D);
  fprintf(fout0,"N =%8.3f\n",N);
  fprintf(fout0,"Qx=%8.3f\n",Qx);
  fprintf(fout0,"Qy=%8.3f\n",Qy);
  fprintf(fout0,"Mx=%8.3f\n",Mx);
  fprintf(fout0,"My=%8.3f\n",My);
  */

  /*if(N<0.0) N=0.0;*/

  /*CONSTANTS*/
  A  =B*D;
  Ix =B*D*D*D/12.0;
  Iy =B*B*B*D/12.0;
  Iw =0.0;
  Jz =B*B*B*D/3.0;
  Zx =B*D*D/6.0;
  Zy =B*B*D/6.0;
  Zpx=B*D*D/4.0;
  Zpy=B*B*D/4.0;
  ix=sqrt(Ix/A);
  iy=sqrt(Iy/A);

/*L-=50.0;*/ /*ZOORASIA äKçÇÇ500â∫Ç∞ÇΩèÍçá*/

  /*ALLOWABLE AXIAL STRESS*/

  /*WEAK AXIS*/
  fact=1.0;

  if(!strcmp(prj,"zoo") && period==PLONG) /*ZOORASIA*/
  {
    /*fact=0.6;*/ fact=1.0/sqrt(2.0); /*FIX,FIX*/
  }
  if(!strcmp(prj,"zoo") && period==PSHORT) /*ZOORASIA*/
  {
    /*fact=0.6;*/ fact=1.0/sqrt(2.0); /*HINGE,FIX*/
  }

/*HIRAKATA*/
/*if(!strcmp(prj,"hira") && sect.code==201) fact=1.1;*/
/*if(!strcmp(prj,"hira") && sect.code==202) fact=1.1;*/
if(!strcmp(prj,"hira") && sect.code==501) fact=/*0.5*/sqrt(2.0);
if(!strcmp(prj,"hira") && sect.code==502) fact=sqrt(2.0);
if(!strcmp(prj,"hira") && sect.code==503) fact=sqrt(2.0);

if(!strcmp(prj,"zooii") && sect.code==501) fact=0.6;

if(!strcmp(prj,"naka") && sect.code==511) fact=1.0/sqrt(2.0);
if(!strcmp(prj,"naka") && sect.code==512) fact=1.0/sqrt(2.0);
if(!strcmp(prj,"naka") && sect.code>=901) fact=1.0/sqrt(2.0);

if(!strcmp(prj,"izu") && sect.code==202) fact=0.5 /*1.0/sqrt(2.0)*/;
if(!strcmp(prj,"izu") && sect.code>=901) fact=0.5 /*1.0/sqrt(2.0)*/;

if(!strcmp(prj,"yagi") && sect.code>=901) fact=0.5 /*1.0/sqrt(2.0)*/;

/*
if(!strcmp(prj,"nade") && sect.code==201) fact=1.1/2.0;
if(!strcmp(prj,"nade") && sect.code==202) fact=1.1/2.0;
if(!strcmp(prj,"nade") && sect.code>=901) fact=1.1/2.0;
*/

if(!strcmp(prj,"mogi") && sect.code==501) fact=0.6;

  /*Lk=fact*L;*/ /*BUCKLING LENGTH [cm]*/

  if(bblength[1]!=0.0) Lk=bblength[1]; /*BUCKLING LENGTH [cm]*/
  else if(bbfact[1]!=0.0) Lk=bbfact[1]*Lm;
  else Lk=fact*L;

if(!strcmp(prj,"hira") && sect.code==201) Lk= 50.0; /*34.0*/
if(!strcmp(prj,"hira") && sect.code==202) Lk= 80.0; /*70.0*/
if(!strcmp(prj,"hira") && sect.code==203) Lk= 33.0; /*33.0*/
if(!strcmp(prj,"hira") && sect.code==204) Lk= 59.0; /*59.0*/
if(!strcmp(prj,"hira") && sect.code==205) Lk= 60.0; /*60.0*/
if(!strcmp(prj,"hira") && sect.code==206) Lk= 63.0; /*63.0*/
if(!strcmp(prj,"hira") && sect.code==207) Lk= 70.0; /*70.0*/

if(!strcmp(prj,"zooii") && sect.code==201) Lk=200.0;
if(!strcmp(prj,"zooii") && sect.code==202) Lk=90.0;
if(!strcmp(prj,"zooii") && sect.code==502) Lk=200.0;

if(!strcmp(prj,"gyo") && sect.code==201) Lk=150.0;
if(!strcmp(prj,"gyo") && sect.code==502) Lk=150.0;

if(!strcmp(prj,"hakoii") && sect.code==201) Lk=86.0;
if(!strcmp(prj,"hakoii") && sect.code==202) Lk=86.0;
if(!strcmp(prj,"hakoii") && sect.code==203) Lk=86.0;
if(!strcmp(prj,"hakoii") && sect.code==511) Lk=43.0;
if(!strcmp(prj,"hakoii") && sect.code==512) Lk=43.0;
if(!strcmp(prj,"hakoii") && sect.code==513) Lk=344.0;
if(!strcmp(prj,"hakoii") && sect.code==514) Lk=172.0;
if(!strcmp(prj,"hakoii") && sect.code==515) Lk=172.0;
if(!strcmp(prj,"hakoii") && sect.code==519) Lk=344.0;

if(!strcmp(prj,"kiku") && sect.code==202) Lk=63.0; /* =63.0x1.5=94.5 */
if(!strcmp(prj,"kiku") && sect.code==203) Lk=63.0;
if(!strcmp(prj,"kiku") && sect.code==205) Lk=0.5*L;
if(!strcmp(prj,"kiku") && sect.code==206) Lk=0.5*L;
if(!strcmp(prj,"kiku") && sect.code==501) Lk=0.5*L;

if(!strcmp(prj,"tohu") && sect.code==202) Lk=65.0;
if(!strcmp(prj,"tohu") && sect.code==203) Lk=65.0;

if(!strcmp(prj,"yagi") && sect.code==501) Lk= 60.0;
if(!strcmp(prj,"yagi") && sect.code==201) Lk= 90.0; /*r=6500*/
if(!strcmp(prj,"yagi") && sect.code==202) Lk= 70.0; /*r=4000*/
if(!strcmp(prj,"yagi") && sect.code==203) Lk= 40.0; /*r=5500 Lk=83cm*/
if(!strcmp(prj,"yagi") && sect.code==204) Lk= 40.0; /*r=2000+RIB*/
if(!strcmp(prj,"yagi") && sect.code==205) Lk=100.0; /*r=7500*/
if(!strcmp(prj,"yagi") && sect.code==206) Lk= 80.0; /*r=5000*/
if(!strcmp(prj,"yagi") && sect.code==207) Lk= 60.0; /*r=3000*/
if(!strcmp(prj,"yagi") && sect.code==208) Lk= 80.0; /*r=5000*/
if(!strcmp(prj,"yagi") && sect.code==209) Lk= 40.0; /*r=4000 Lk=83cm*/
if(!strcmp(prj,"yagi") && sect.code==210) Lk= 70.0; /*r=4000*/
/*if(!strcmp(prj,"yagi") && sect.code==921) Lk= 70.0;*/ /*r=4000*/
/*if(!strcmp(prj,"yagi") && sect.code==928) Lk= 40.0;*/ /*r=5500*/
/*if(!strcmp(prj,"yagi") && sect.code==929) Lk= 50.0;*/ /*r=2000*/

if(!strcmp(prj,"nade") && sect.code==201) Lk=45.0;
if(!strcmp(prj,"nade") && sect.code==202) Lk=45.0;
if(!strcmp(prj,"nade") && sect.code==501) Lk=45.0;
if(!strcmp(prj,"nade") && sect.code==502) Lk=45.0;

  lamda=Lk/iy;
  Lamda=sqrt(PI*PI*E/0.6/Fc);
  nyu=1.5+(lamda*lamda)/(Lamda*Lamda)/1.5;
  if(nyu>2.17) nyu=2.17; /*kaza 20210822*/

  if(lamda<=Lamda) fc=(1.0-0.4*(lamda*lamda)/(Lamda*Lamda))*Fc/1000.0/nyu;
  else             fc=0.277*Fc/1000.0/(lamda*lamda)*(Lamda*Lamda);

  if(period==PSHORT) fc*=1.5;

  if(N>=0.0) Nay=fc*A; /*[tf]*/
  else
  {
    Nay=Ft/1000.0*A; if(period==PLONG) Nay/=1.5;
  }

  Na=Nay;

  /*STRONG AXIS*/
  fact=1.0;

  /*Lk=fact*L;*/ /*BUCKLING LENGTH [cm]*/

  if(bblength[0]!=0.0) Lk=bblength[0]; /*BUCKLING LENGTH [cm]*/
  else if(bbfact[0]!=0.0) Lk=bbfact[0]*Lm;
  else Lk=fact*L;

if(!strcmp(prj,"kiku") && sect.code==202) Lk=261.0;
if(!strcmp(prj,"kiku") && sect.code==203) Lk=261.0;
if(!strcmp(prj,"kiku") && sect.code==501) Lk=0.5*L;

/*if(!strcmp(prj,"tohu") && sect.code==202) Lk=65.0;*/
/*if(!strcmp(prj,"tohu") && sect.code==203) Lk=65.0;*/

if(!strcmp(prj,"hakoii") && sect.code==201) Lk=688.0;
if(!strcmp(prj,"hakoii") && sect.code==202) Lk=344.0*sqrt(2.0); /*688.0*/
if(!strcmp(prj,"hakoii") && sect.code==203) Lk=688.0;
if(!strcmp(prj,"hakoii") && sect.code==511) Lk=344.0;
if(!strcmp(prj,"hakoii") && sect.code==512) Lk=344.0;
if(!strcmp(prj,"hakoii") && sect.code==513) Lk=688.0;
if(!strcmp(prj,"hakoii") && sect.code==514) Lk=344.0;
if(!strcmp(prj,"hakoii") && sect.code==515) Lk=344.0;
if(!strcmp(prj,"hakoii") && sect.code==519) Lk=688.0;

  lamda=Lk/ix;
  Lamda=sqrt(PI*PI*E/0.6/Fc);
  nyu=1.5+(lamda*lamda)/(Lamda*Lamda)/1.5;
  if(nyu>2.17) nyu=2.17; /*kaza 20210822*/

  if(lamda<=Lamda) fc=(1.0-0.4*(lamda*lamda)/(Lamda*Lamda))*Fc/1000.0/nyu;
  else             fc=0.277*Fc/1000.0/(lamda*lamda)*(Lamda*Lamda);

  if(period==PSHORT) fc*=1.5;

  if(N>=0.0) Nax=fc*A; /*[tf]*/
  else
  {
    Nax=Ft/1000.0*A; if(period==PLONG) Nax/=1.5;
  }

  if(Na>Nax) Na=Nax;

  /*fprintf(fout0,"Lk/iy=%8.3f/%8.3f=%8.3f\n",Lk,iy,lamda);
  fprintf(fout0,"Na=%8.3f\n",Na);*/

  /*STRESS LIMITS*/
  if(N>=0.0) Np=Fc/1000.0*A;
  else       Np=Ft/1000.0*A;

  if(period==PLONG) Np/=1.5;

  /*STRONG AXIS*/

  fact=1.0; /*HINGE,HINGE*/

if(!strcmp(prj,"zoo")) fact=1.0/sqrt(2.0); /*HINGE,FIX*//*ZOORASIA*/

  /*Lk=fact*L;*/ /*BUCKLING LENGTH [cm]*/

  if(bblength[0]!=0.0) Lk=bblength[0]; /*BUCKLING LENGTH [cm]*/
  else if(bbfact[0]!=0.0) Lk=bbfact[0]*Lm;
  else Lk=fact*L;

/*HIRAKATA*/
if(!strcmp(prj,"hira") && sect.code==201) Lk=660;
if(!strcmp(prj,"hira") && sect.code==202) Lk=660;
if(!strcmp(prj,"hira") && sect.code==203) Lk=660;
if(!strcmp(prj,"hira") && sect.code==204) Lk=660;
if(!strcmp(prj,"hira") && sect.code==205) Lk=660;
if(!strcmp(prj,"hira") && sect.code==206) Lk=660;
if(!strcmp(prj,"hira") && sect.code==207) Lk=660;

if(!strcmp(prj,"zooii") && sect.code==201) Lk=200.0;
if(!strcmp(prj,"zooii") && sect.code==202) Lk=90.0;
if(!strcmp(prj,"zooii") && sect.code==502) Lk=200.0;

if(!strcmp(prj,"gyo") && sect.code==201) Lk=150.0;
if(!strcmp(prj,"gyo") && sect.code==502) Lk=150.0;

if(!strcmp(prj,"hakoii") && sect.code==201) Lk=688.0;
if(!strcmp(prj,"hakoii") && sect.code==202) Lk=344.0*sqrt(2.0); /*688.0*/
if(!strcmp(prj,"hakoii") && sect.code==203) Lk=688.0;
if(!strcmp(prj,"hakoii") && sect.code==511) Lk=344.0;
if(!strcmp(prj,"hakoii") && sect.code==512) Lk=344.0;
if(!strcmp(prj,"hakoii") && sect.code==513) Lk=688.0;
if(!strcmp(prj,"hakoii") && sect.code==514) Lk=344.0;
if(!strcmp(prj,"hakoii") && sect.code==515) Lk=344.0;
if(!strcmp(prj,"hakoii") && sect.code==519) Lk=688.0;

if(!strcmp(prj,"kiku") && sect.code==202) Lk=261.0;
if(!strcmp(prj,"kiku") && sect.code==203) Lk=261.0;
if(!strcmp(prj,"kiku") && sect.code==501) Lk=0.5*L;

/*if(!strcmp(prj,"tohu") && sect.code==202) Lk=65.0;*/
/*if(!strcmp(prj,"tohu") && sect.code==203) Lk=65.0;*/

  Ncrx=PI*PI*E*Ix/Lk/Lk/1000.0; if(period==PLONG) Ncrx/=1.5;

  /*WEAK AXIS*/
  fact=1.0; /*HINGE,HINGE*/
if(!strcmp(prj,"zoo") && period==PLONG)
{
  /*fact=0.6;*/ fact=1.0/sqrt(2.0); /*FIX,FIX*/
}
if(!strcmp(prj,"zoo") && period==PSHORT)
{
  /*fact=0.6;*/ fact=1.0/sqrt(2.0); /*HINGE,FIX*/
}

/*HIRAKATA*/
/*if(!strcmp(prj,"hira") && sect.code==201) fact=1.1;*/
/*if(!strcmp(prj,"hira") && sect.code==202) fact=1.1;*/
if(!strcmp(prj,"hira") && sect.code==501) fact=/*0.5*/sqrt(2.0);
if(!strcmp(prj,"hira") && sect.code==502) fact=sqrt(2.0);
if(!strcmp(prj,"hira") && sect.code==503) fact=sqrt(2.0);

if(!strcmp(prj,"zooii") && sect.code==501) fact=0.6;

if(!strcmp(prj,"naka") && sect.code==511) fact=1.0/sqrt(2.0);
if(!strcmp(prj,"naka") && sect.code==512) fact=1.0/sqrt(2.0);
if(!strcmp(prj,"naka") && sect.code>=901) fact=1.0/sqrt(2.0);

if(!strcmp(prj,"izu") && sect.code==202) fact=0.5 /*1.0/sqrt(2.0)*/;
if(!strcmp(prj,"izu") && sect.code>=901) fact=0.5 /*1.0/sqrt(2.0)*/;

if(!strcmp(prj,"yagi") && sect.code>=901) fact=0.5 /*1.0/sqrt(2.0)*/;

/*
if(!strcmp(prj,"nade") && sect.code==201) fact=1.1/2.0;
if(!strcmp(prj,"nade") && sect.code==202) fact=1.1/2.0;
if(!strcmp(prj,"nade") && sect.code>=901) fact=1.1/2.0;
*/

if(!strcmp(prj,"mogi") && sect.code==501) fact=0.6;

  /*Lk=fact*L;*/ /*BUCKLING LENGTH [cm]*/

  if(bblength[1]!=0.0) Lk=bblength[1]; /*BUCKLING LENGTH [cm]*/
  else if(bbfact[1]!=0.0) Lk=bbfact[1]*Lm;
  else Lk=fact*L;

if(!strcmp(prj,"hira") && sect.code==201) Lk= 50.0; /*34.0*/
if(!strcmp(prj,"hira") && sect.code==202) Lk= 80.0; /*70.0*/
if(!strcmp(prj,"hira") && sect.code==203) Lk= 33.0; /*33.0*/
if(!strcmp(prj,"hira") && sect.code==204) Lk= 59.0; /*59.0*/
if(!strcmp(prj,"hira") && sect.code==205) Lk= 60.0; /*60.0*/
if(!strcmp(prj,"hira") && sect.code==206) Lk= 63.0; /*63.0*/
if(!strcmp(prj,"hira") && sect.code==207) Lk= 70.0; /*70.0*/

if(!strcmp(prj,"zooii") && sect.code==201) Lk=200.0;
if(!strcmp(prj,"zooii") && sect.code==202) Lk=90.0;
if(!strcmp(prj,"zooii") && sect.code==502) Lk=200.0;

if(!strcmp(prj,"gyo") && sect.code==201) Lk=150.0;
if(!strcmp(prj,"gyo") && sect.code==502) Lk=150.0;

if(!strcmp(prj,"hakoii") && sect.code==201) Lk=86.0;
if(!strcmp(prj,"hakoii") && sect.code==202) Lk=86.0;
if(!strcmp(prj,"hakoii") && sect.code==203) Lk=86.0;
if(!strcmp(prj,"hakoii") && sect.code==511) Lk=43.0;
if(!strcmp(prj,"hakoii") && sect.code==512) Lk=43.0;
if(!strcmp(prj,"hakoii") && sect.code==513) Lk=344.0;
if(!strcmp(prj,"hakoii") && sect.code==514) Lk=172.0;
if(!strcmp(prj,"hakoii") && sect.code==515) Lk=172.0;
if(!strcmp(prj,"hakoii") && sect.code==519) Lk=344.0;

if(!strcmp(prj,"kiku") && sect.code==202) Lk=69.3; /* =63.0x1.5=94.5 */
if(!strcmp(prj,"kiku") && sect.code==203) Lk=69.3;
if(!strcmp(prj,"kiku") && sect.code==205) Lk=0.5*L;
if(!strcmp(prj,"kiku") && sect.code==206) Lk=0.5*L;
if(!strcmp(prj,"kiku") && sect.code==501) Lk=0.5*L;

if(!strcmp(prj,"tohu") && sect.code==202) Lk=65.0;
if(!strcmp(prj,"tohu") && sect.code==203) Lk=65.0;

if(!strcmp(prj,"yagi") && sect.code==501) Lk= 60.0;
if(!strcmp(prj,"yagi") && sect.code==201) Lk= 90.0;
if(!strcmp(prj,"yagi") && sect.code==202) Lk= 70.0;
if(!strcmp(prj,"yagi") && sect.code==203) Lk= 40.0;
if(!strcmp(prj,"yagi") && sect.code==204) Lk= 40.0;
if(!strcmp(prj,"yagi") && sect.code==205) Lk=100.0;
if(!strcmp(prj,"yagi") && sect.code==206) Lk= 80.0;
if(!strcmp(prj,"yagi") && sect.code==207) Lk= 60.0;
if(!strcmp(prj,"yagi") && sect.code==208) Lk= 80.0;
if(!strcmp(prj,"yagi") && sect.code==209) Lk= 40.0;
if(!strcmp(prj,"yagi") && sect.code==210) Lk= 70.0;
/*if(!strcmp(prj,"yagi") && sect.code==921) Lk= 70.0;*/
/*if(!strcmp(prj,"yagi") && sect.code==928) Lk= 40.0;*/
/*if(!strcmp(prj,"yagi") && sect.code==929) Lk= 50.0;*/

if(!strcmp(prj,"nade") && sect.code==201) Lk=45.0;
if(!strcmp(prj,"nade") && sect.code==202) Lk=45.0;
if(!strcmp(prj,"nade") && sect.code==501) Lk=45.0;
if(!strcmp(prj,"nade") && sect.code==502) Lk=45.0;

  Ncry=PI*PI*E*Iy/Lk/Lk/1000.0;  if(period==PLONG) Ncry/=1.5;

  if(sect.stype==STYPE_S)     Qax=Fs/1000.0/sqrt(3.0)*A/1.5;
  if(sect.stype==STYPE_GLASS) Qax=Fs/1000.0*A/1.5;
  if(sect.stype==STYPE_ACRYL) Qax=Fs/1000.0*A/1.5;

  if(period==PLONG) Qax /=1.5;

  Qay=Qax;
  Myx=Zx*Fb/100000.0;            if(period==PLONG) Myx /=1.5;
  Myy=Zy*Fb/100000.0;            if(period==PLONG) Myy /=1.5;
  Mpx=Zpx*Fb/100000.0;           if(period==PLONG) Mpx /=1.5;
  Mpy=Zpy*Fb/100000.0;           if(period==PLONG) Mpy /=1.5;
  /*
  fprintf(fout0,"Np  =%8.3f\n",Np);
  fprintf(fout0,"Ncrx=%8.3f\n",Ncrx);
  fprintf(fout0,"Ncry=%8.3f\n",Ncry);
  fprintf(fout0,"Qax =%8.3f\n",Qax);
  fprintf(fout0,"Qay =%8.3f\n",Qay);
  fprintf(fout0,"Myx =%8.3f\n",Myx);
  fprintf(fout0,"Myy =%8.3f\n",Myy);
  fprintf(fout0,"Mpx =%8.3f\n",Mpx);
  fprintf(fout0,"Mpy =%8.3f\n",Mpy);
  */

  /*LATERAL LINEAR BUCKLING MOMENT*/
if(!strcmp(prj,"yachi") && sect.code==501) L=Lm;
if(!strcmp(prj,"yachi") && sect.code==502) L=Lm;
  g=0.0;
  k=0.0;
  /*ku=1.0/sqrt(2.0);*/
  /*kb=1.0/sqrt(2.0);*/

/*HIRAKATA*/
/*if(!strcmp(prj,"hira") && sect.code==201) ku=1.1;*/
/*if(!strcmp(prj,"hira") && sect.code==202) ku=1.1;*/
if(!strcmp(prj,"hira") && sect.code==501) ku=/*0.5*/sqrt(2.0);
if(!strcmp(prj,"hira") && sect.code==502) ku=sqrt(2.0);
if(!strcmp(prj,"hira") && sect.code==503) ku=sqrt(2.0);

if(!strcmp(prj,"hira") && sect.code==201) {Lk= 50.0; ku=Lk/L;} /*34.0*/
if(!strcmp(prj,"hira") && sect.code==202) {Lk= 80.0; ku=Lk/L;} /*70.0*/
if(!strcmp(prj,"hira") && sect.code==203) {Lk= 33.0; ku=Lk/L;} /*33.0*/
if(!strcmp(prj,"hira") && sect.code==204) {Lk= 59.0; ku=Lk/L;} /*59.0*/
if(!strcmp(prj,"hira") && sect.code==205) {Lk= 60.0; ku=Lk/L;} /*60.0*/
if(!strcmp(prj,"hira") && sect.code==206) {Lk= 63.0; ku=Lk/L;} /*63.0*/
if(!strcmp(prj,"hira") && sect.code==207) {Lk= 70.0; ku=Lk/L;} /*70.0*/

if(!strcmp(prj,"zooii") && sect.code==201) {Lk=200.0; ku=Lk/L;}
if(!strcmp(prj,"zooii") && sect.code==202) {Lk=90.0; ku=Lk/L;}
if(!strcmp(prj,"zooii") && sect.code==502) {Lk=200.0; ku=Lk/L;}

if(!strcmp(prj,"gyo") && sect.code==201) {Lk=150.0; ku=Lk/L;}
if(!strcmp(prj,"gyo") && sect.code==502) {Lk=150.0; ku=Lk/L;}

  /*q=0.0;*/
  M0M1=0.0;
  M2M1=1.0;
  gzai=0.283*(1+M2M1*M2M1)+0.434*M2M1+0.868*M0M1*(1+M2M1)+0.78*M0M1*M0M1;
  C1=1.0/sqrt(gzai);
  C2=0.405*M2M1/sqrt(gzai);
  C3=(0.5+0.5*M2M1+0.464*M0M1)/sqrt(gzai);

  /*STRONG AXIS*/

  if(bblength[0]!=0.0) {Lk=bblength[0]; ku=Lk/L;}
  else if(bbfact[0]!=0.0) ku=bbfact[0];
  else ku=1.0;
  if(btlength[0]!=0.0) {Lk=btlength[0]; kb=Lk/L;}
  else if(btfact[0]!=0.0) kb=btfact[0];
  else kb=1.0/sqrt(2.0);

if(!strcmp(prj,"hakoii") && sect.code==201) {Lk=688.0; ku=Lk/L;}
if(!strcmp(prj,"hakoii") && sect.code==202) {Lk=344.0*sqrt(2.0); /*688.0*/ ku=Lk/L;}
if(!strcmp(prj,"hakoii") && sect.code==203) {Lk=688.0; ku=Lk/L;}
if(!strcmp(prj,"hakoii") && sect.code==511) {Lk=344.0; ku=Lk/L;}
if(!strcmp(prj,"hakoii") && sect.code==512) {Lk=344.0; ku=Lk/L;}
if(!strcmp(prj,"hakoii") && sect.code==513) {Lk=688.0; ku=Lk/L;}
if(!strcmp(prj,"hakoii") && sect.code==514) {Lk=344.0; ku=Lk/L;}
if(!strcmp(prj,"hakoii") && sect.code==515) {Lk=344.0; ku=Lk/L;}
if(!strcmp(prj,"hakoii") && sect.code==519) {Lk=688.0; ku=Lk/L;}

if(!strcmp(prj,"naka") && sect.code==511) ku=1.0;
if(!strcmp(prj,"naka") && sect.code==512) ku=1.0;
if(!strcmp(prj,"naka") && sect.code>=901) ku=1.0;

if(!strcmp(prj,"izu") && sect.code==202) ku=1.0;
if(!strcmp(prj,"izu") && sect.code>=901) ku=1.0;

if(!strcmp(prj,"kiku") && sect.code==202) {Lk=261.0; ku=Lk/L;}
if(!strcmp(prj,"kiku") && sect.code==203) {Lk=261.0; ku=Lk/L;}
if(!strcmp(prj,"kiku") && sect.code==501) {Lk=0.5*L; ku=Lk/L;}

/*if(!strcmp(prj,"tohu") && sect.code==202) {Lk=65.0; ku=Lk/L;}*/
/*if(!strcmp(prj,"tohu") && sect.code==203) {Lk=65.0; ku=Lk/L;}*/

if(!strcmp(prj,"nade") && sect.code==201) ku=1.0;
if(!strcmp(prj,"nade") && sect.code==202) ku=1.0;
if(!strcmp(prj,"nade") && sect.code>=901) ku=1.0;

if(!strcmp(prj,"mogi") && sect.code==501) ku=1.0;

  Mcrex=C1*PI*PI*E*Ix/(ku*L)/(ku*L)
       *((C2*g+C3*k)
         +sqrt((C2*g+C3*k)*(C2*g+C3*k)
               +Iw/Ix*(ku/kb)*(ku/kb)
               +G*Jz*(ku*L)*(ku*L)/(PI*PI*E*Ix)))/100000.0;

  /*WEAK AXIS*/

  if(bblength[1]!=0.0) {Lk=bblength[1]; ku=Lk/L;}
  else if(bbfact[1]!=0.0) ku=bbfact[1];
  else ku=1.0;
  if(btlength[1]!=0.0) {Lk=btlength[1]; kb=Lk/L;}
  else if(btfact[1]!=0.0) kb=btfact[1];
  else kb=1.0/sqrt(2.0);

if(!strcmp(prj,"hakoii") && sect.code==201) {Lk=86.0; ku=Lk/L;}
if(!strcmp(prj,"hakoii") && sect.code==202) {Lk=86.0; ku=Lk/L;}
if(!strcmp(prj,"hakoii") && sect.code==203) {Lk=86.0; ku=Lk/L;}
if(!strcmp(prj,"hakoii") && sect.code==511) {Lk=43.0; ku=Lk/L;}
if(!strcmp(prj,"hakoii") && sect.code==512) {Lk=43.0; ku=Lk/L;}
if(!strcmp(prj,"hakoii") && sect.code==513) {Lk=344.0; ku=Lk/L;}
if(!strcmp(prj,"hakoii") && sect.code==514) {Lk=172.0; ku=Lk/L;}
if(!strcmp(prj,"hakoii") && sect.code==515) {Lk=172.0; ku=Lk/L;}
if(!strcmp(prj,"hakoii") && sect.code==519) {Lk=344.0; ku=Lk/L;}

if(!strcmp(prj,"naka") && sect.code==511) ku=1.0/sqrt(2.0);
if(!strcmp(prj,"naka") && sect.code==512) ku=1.0/sqrt(2.0);
if(!strcmp(prj,"naka") && sect.code>=901) ku=1.0/sqrt(2.0);

if(!strcmp(prj,"izu") && sect.code==202) ku=0.5 /*1.0/sqrt(2.0)*/;
if(!strcmp(prj,"izu") && sect.code>=901) ku=0.5 /*1.0/sqrt(2.0)*/;

if(!strcmp(prj,"yagi") && sect.code>=901) ku=0.5 /*1.0/sqrt(2.0)*/;

if(!strcmp(prj,"kiku") && sect.code==202) {Lk=69.3; ku=Lk/L;}
if(!strcmp(prj,"kiku") && sect.code==203) {Lk=69.3; ku=Lk/L;}
if(!strcmp(prj,"kiku") && sect.code==205) {Lk=0.5*L; ku=Lk/L;}
if(!strcmp(prj,"kiku") && sect.code==206) {Lk=0.5*L; ku=Lk/L;}
if(!strcmp(prj,"kiku") && sect.code==501) {Lk=0.5*L; ku=Lk/L;}

if(!strcmp(prj,"tohu") && sect.code==202) {Lk=65.0; ku=Lk/L;}
if(!strcmp(prj,"tohu") && sect.code==203) {Lk=65.0; ku=Lk/L;}

if(!strcmp(prj,"yagi") && sect.code==501) {Lk= 60.0; ku=Lk/L;}
if(!strcmp(prj,"yagi") && sect.code==201) {Lk= 90.0; ku=Lk/L;}
if(!strcmp(prj,"yagi") && sect.code==202) {Lk= 70.0; ku=Lk/L;}
if(!strcmp(prj,"yagi") && sect.code==203) {Lk= 40.0; ku=Lk/L;}
if(!strcmp(prj,"yagi") && sect.code==204) {Lk= 40.0; ku=Lk/L;}
if(!strcmp(prj,"yagi") && sect.code==205) {Lk=100.0; ku=Lk/L;}
if(!strcmp(prj,"yagi") && sect.code==206) {Lk= 80.0; ku=Lk/L;}
if(!strcmp(prj,"yagi") && sect.code==207) {Lk= 60.0; ku=Lk/L;}
if(!strcmp(prj,"yagi") && sect.code==208) {Lk= 80.0; ku=Lk/L;}
if(!strcmp(prj,"yagi") && sect.code==209) {Lk= 40.0; ku=Lk/L;}
if(!strcmp(prj,"yagi") && sect.code==210) {Lk= 70.0; ku=Lk/L;}
if(!strcmp(prj,"yagi") && sect.code==921) {Lk= 70.0; ku=Lk/L;}
if(!strcmp(prj,"yagi") && sect.code==928) {Lk= 40.0; ku=Lk/L;}
if(!strcmp(prj,"yagi") && sect.code==929) {Lk= 50.0; ku=Lk/L;}

if(!strcmp(prj,"nade") && sect.code==201) {Lk=45.0; ku=Lk/L;}
if(!strcmp(prj,"nade") && sect.code==202) {Lk=45.0; ku=Lk/L;}
if(!strcmp(prj,"nade") && sect.code==501) {Lk=45.0; ku=Lk/L;}
if(!strcmp(prj,"nade") && sect.code==502) {Lk=45.0; ku=Lk/L;}
if(!strcmp(prj,"nade") && sect.code>=901) ku=1.1/2.0;

if(!strcmp(prj,"mogi") && sect.code==501) ku=0.6;

  Mcrey=C1*PI*PI*E*Iy/(ku*L)/(ku*L)
       *((C2*g+C3*k)
         +sqrt((C2*g+C3*k)*(C2*g+C3*k)
               +Iw/Iy*(ku/kb)*(ku/kb)
               +G*Jz*(ku*L)*(ku*L)/(PI*PI*E*Iy)))/100000.0;

  if(Mcrex<=Mcrey) Mcre=Mcrex;
  else             Mcre=Mcrey;

  if(period==PLONG) Mcre/=1.5;

  /*fprintf(fout0,"Mcre=%8.3f\n",Mcre);*/

  /*LATERAL NONLINEAR BUCKLING MOMENT*/
  lamda=sqrt(Mpx/Mcre);
  fact=0.7+0.3*(1-0.7*lamda*lamda)/(0.61+0.3*M2M1+0.07*M2M1*M2M1);

  if(lamda>sqrt(1/0.7)) Mcrp=Mcre;
  else                  Mcrp=fact*Mpx;

  /*fprintf(fout0,"Mcrp=%8.3f\n",Mcrp);*/

  if(N>=0.0) /*COMPRESSION*/
  {
    /*P DELTA*/
    Kx=sqrt(N*1000.0/E/Ix);
    Ky=sqrt(N*1000.0/E/Iy);

    /*Cx:STRONG AXIS*/

    if(bblength[0]!=0.0) {Lk=bblength[0]; ku=Lk/L;}
    else if(bbfact[0]!=0.0) ku=bbfact[0];
    else ku=1.0; /*?*/
    if(btlength[0]!=0.0) {Lk=btlength[0]; kb=Lk/L;}
    else if(btfact[0]!=0.0) kb=btfact[0];
    else kb=1.0/sqrt(2.0); /*?*/

if(!strcmp(prj,"hakoii") && sect.code==201) {Lk=688.0; ku=Lk/L;}
if(!strcmp(prj,"hakoii") && sect.code==202) {Lk=344.0*sqrt(2.0); /*688.0*/ ku=Lk/L;}
if(!strcmp(prj,"hakoii") && sect.code==203) {Lk=688.0; ku=Lk/L;}
if(!strcmp(prj,"hakoii") && sect.code==511) {Lk=344.0; ku=Lk/L;}
if(!strcmp(prj,"hakoii") && sect.code==512) {Lk=344.0; ku=Lk/L;}
if(!strcmp(prj,"hakoii") && sect.code==513) {Lk=688.0; ku=Lk/L;}
if(!strcmp(prj,"hakoii") && sect.code==514) {Lk=344.0; ku=Lk/L;}
if(!strcmp(prj,"hakoii") && sect.code==515) {Lk=344.0; ku=Lk/L;}
if(!strcmp(prj,"hakoii") && sect.code==519) {Lk=688.0; ku=Lk/L;}

if(!strcmp(prj,"naka") && sect.code==511) ku=1.0;
if(!strcmp(prj,"naka") && sect.code==512) ku=1.0;
if(!strcmp(prj,"naka") && sect.code>=901) ku=1.0;

if(!strcmp(prj,"izu") && sect.code==202) ku=1.0;
if(!strcmp(prj,"izu") && sect.code>=901) ku=1.0;

if(!strcmp(prj,"kiku") && sect.code==202) {Lk=261.0; ku=Lk/L;}
if(!strcmp(prj,"kiku") && sect.code==203) {Lk=261.0; ku=Lk/L;}
if(!strcmp(prj,"kiku") && sect.code==501) {Lk=0.5*L; ku=Lk/L;}

/*if(!strcmp(prj,"tohu") && sect.code==202) {Lk=65.0; ku=Lk/L;}*/
/*if(!strcmp(prj,"tohu") && sect.code==203) {Lk=65.0; ku=Lk/L;}*/

if(!strcmp(prj,"nade") && sect.code==201) ku=1.0;
if(!strcmp(prj,"nade") && sect.code==202) ku=1.0;
if(!strcmp(prj,"nade") && sect.code>=901) ku=1.0;

if(!strcmp(prj,"mogi") && sect.code==501) ku=1.0;

    /*cx=cos(Kx*L);*/
    cx=cos(Kx*Lk);

    /*if(M2M1>=cx) alphax=sqrt(1.0+M2M1*M2M1-2.0*M2M1*cos(Kx*L))/sin(Kx*L);*/
    if(N==0.0 || Mx==0.0) alphax=1.0;
    else if(M2M1>=cx) alphax=sqrt(1.0+M2M1*M2M1-2.0*M2M1*cos(Kx*Lk))/sin(Kx*Lk);
    else return 10.0;

    /*Cy:WEAK AXIS*/

    if(bblength[1]!=0.0) {Lk=bblength[1]; ku=Lk/L;}
    else if(bbfact[1]!=0.0) ku=bbfact[1];
    else ku=1.0; /*?*/
    if(btlength[1]!=0.0) {Lk=btlength[1]; kb=Lk/L;}
    else if(btfact[1]!=0.0) kb=btfact[1];
    else kb=1.0/sqrt(2.0); /*?*/

if(!strcmp(prj,"hakoii") && sect.code==201) {Lk=86.0; ku=Lk/L;}
if(!strcmp(prj,"hakoii") && sect.code==202) {Lk=86.0; ku=Lk/L;}
if(!strcmp(prj,"hakoii") && sect.code==203) {Lk=86.0; ku=Lk/L;}
if(!strcmp(prj,"hakoii") && sect.code==511) {Lk=43.0; ku=Lk/L;}
if(!strcmp(prj,"hakoii") && sect.code==512) {Lk=43.0; ku=Lk/L;}
if(!strcmp(prj,"hakoii") && sect.code==513) {Lk=344.0; ku=Lk/L;}
if(!strcmp(prj,"hakoii") && sect.code==514) {Lk=172.0; ku=Lk/L;}
if(!strcmp(prj,"hakoii") && sect.code==515) {Lk=172.0; ku=Lk/L;}
if(!strcmp(prj,"hakoii") && sect.code==519) {Lk=344.0; ku=Lk/L;}

if(!strcmp(prj,"naka") && sect.code==511) ku=1.0/sqrt(2.0);
if(!strcmp(prj,"naka") && sect.code==512) ku=1.0/sqrt(2.0);
if(!strcmp(prj,"naka") && sect.code>=901) ku=1.0/sqrt(2.0);

if(!strcmp(prj,"izu") && sect.code==202) ku=0.5 /*1.0/sqrt(2.0)*/;
if(!strcmp(prj,"izu") && sect.code>=901) ku=0.5 /*1.0/sqrt(2.0)*/;

if(!strcmp(prj,"yagi") && sect.code>=901) ku=0.5 /*1.0/sqrt(2.0)*/;

if(!strcmp(prj,"kiku") && sect.code==202) {Lk=69.3; ku=Lk/L;}
if(!strcmp(prj,"kiku") && sect.code==203) {Lk=69.3; ku=Lk/L;}
if(!strcmp(prj,"kiku") && sect.code==205) {Lk=0.5*L; ku=Lk/L;}
if(!strcmp(prj,"kiku") && sect.code==206) {Lk=0.5*L; ku=Lk/L;}
if(!strcmp(prj,"kiku") && sect.code==501) {Lk=0.5*L; ku=Lk/L;}

if(!strcmp(prj,"tohu") && sect.code==202) {Lk=65.0; ku=Lk/L;}
if(!strcmp(prj,"tohu") && sect.code==203) {Lk=65.0; ku=Lk/L;}

if(!strcmp(prj,"yagi") && sect.code==501) {Lk= 60.0; ku=Lk/L;}
if(!strcmp(prj,"yagi") && sect.code==201) {Lk= 90.0; ku=Lk/L;}
if(!strcmp(prj,"yagi") && sect.code==202) {Lk= 70.0; ku=Lk/L;}
if(!strcmp(prj,"yagi") && sect.code==203) {Lk= 40.0; ku=Lk/L;}
if(!strcmp(prj,"yagi") && sect.code==204) {Lk= 40.0; ku=Lk/L;}
if(!strcmp(prj,"yagi") && sect.code==205) {Lk=100.0; ku=Lk/L;}
if(!strcmp(prj,"yagi") && sect.code==206) {Lk= 80.0; ku=Lk/L;}
if(!strcmp(prj,"yagi") && sect.code==207) {Lk= 60.0; ku=Lk/L;}
if(!strcmp(prj,"yagi") && sect.code==208) {Lk= 80.0; ku=Lk/L;}
if(!strcmp(prj,"yagi") && sect.code==209) {Lk= 40.0; ku=Lk/L;}
if(!strcmp(prj,"yagi") && sect.code==210) {Lk= 70.0; ku=Lk/L;}
if(!strcmp(prj,"yagi") && sect.code==921) {Lk= 70.0; ku=Lk/L;}
if(!strcmp(prj,"yagi") && sect.code==928) {Lk= 40.0; ku=Lk/L;}
if(!strcmp(prj,"yagi") && sect.code==929) {Lk= 50.0; ku=Lk/L;}

if(!strcmp(prj,"nade") && sect.code==201) {Lk=45.0; ku=Lk/L;}
if(!strcmp(prj,"nade") && sect.code==202) {Lk=45.0; ku=Lk/L;}
if(!strcmp(prj,"nade") && sect.code==501) {Lk=45.0; ku=Lk/L;}
if(!strcmp(prj,"nade") && sect.code==502) {Lk=45.0; ku=Lk/L;}
if(!strcmp(prj,"nade") && sect.code>=901) ku=1.1/2.0;

if(!strcmp(prj,"mogi") && sect.code==501) ku=0.6;

    /*cy=cos(Ky*L);*/
    cy=cos(Ky*Lk);

    /*if(M2M1>=cy) alphay=sqrt(1.0+M2M1*M2M1-2.0*M2M1*cos(Ky*L))/sin(Ky*L);*/
    if(N==0.0 || My==0.0) alphay=1.0;
    else if(M2M1>=cy) alphay=sqrt(1.0+M2M1*M2M1-2.0*M2M1*cos(Ky*Lk))/sin(Ky*Lk);
    else return 10.0;
  }
  else /*TENSION*/
  {
    alphax=1.0; alphay=1.0;
  }

  if(alphax<0.0 || alphay<0.0) return 10.0;

  /*fprintf(fout0,"Alpha X=%.5f Y=%.5f\n",alphax,alphay);*/

  /*CRITICAL STRESSES*/
  if(fabs(Na)<=Np) *Ncr=Na;
  else             *Ncr=Np;
  *Qcrx=Qax;
  *Qcry=Qay;
  if(Myx<=Mcre && Myx<=Mcrp) *Mcrx=Myx;
  else if(Mcre<=Mcrp)        *Mcrx=Mcre;
  else                       *Mcrx=Mcrp;
  *Mcry=Myy;

  /*RATE*/
  rate1=fabs(N)/(*Ncr)+alphax*Mx/(*Mcrx)+alphay*My/(*Mcry);
  rate2=Qx/(*Qcrx)+Qy/(*Qcry);

  rate=sqrt(rate1*rate1+rate2*rate2);

  /*
  fprintf(fout0,"\n");
  fprintf(fout0,"N/Ncr=%.3f/%.3f=%.3f\n",N,*Ncr,N/(*Ncr));
  fprintf(fout0,"Mx/Mcrx=%.3f/%.3f=%.3f\n",Mx,*Mcrx,Mx/(*Mcrx));
  fprintf(fout0,"My/Mcry=%.3f/%.3f=%.3f\n",My,*Mcry,My/(*Mcry));
  fprintf(fout0,"Qx/Qcrx=%.3f/%.3f=%.3f\n",Qx,*Qcrx,Qx/(*Qcrx));
  fprintf(fout0,"Qy/Qcry=%.3f/%.3f=%.3f\n",Qy,*Qcry,Qy/(*Qcry));
  fprintf(fout0,"RATE=%.5f\n",rate);
  */

  /*CHANGE AXIS*/
  if(sect.sform.H < sect.sform.B)
  {
    Qo=*Qcrx; *Qcrx=*Qcry; *Qcry=Qo;
    Mo=*Mcrx; *Mcrx=*Mcry; *Mcry=Mo;
  }

  return rate;
}/*allowablestressofflatbar*/

double allowablestressofflatbar_bc(int period,
                                struct section sect,
                                double cNcr,
                                double E,
                                double L,       /*ELEMENT LENGTH[cm]*/
                                struct stress *st,        /*[tf,tfm]*/
                                int iend,
                                double *Ncr,
                                double *Qcrx,double *Qcry,
                                double *Mcrx,double *Mcry)
/*FLAT BAR COLUMN.*/
/*BUCKLING LENGTH,FACTOR FROM LST FILE.*/
{
  double poi,G,Fc,Ft,Fb,Fs/*,fu*/;
  double B,D,fact,Lk,A,Ix,Iy,Iw,Jz,Zx,Zy,Zpx,Zpy,ix,iy;
  double lamda,Lamda,nyu,fc,Na,Nax,Nay,Np,Ncrx,Ncry,Qax,Qay,Myx,Myy,Mpx,Mpy;
  double g,k,ku,kb/*,q*/,M0M1,M2M1,gzai,C1,C2,C3,Mcre,Mcrex,Mcrey,Mcrp;
  double Kx,Ky,cx,cy,alphax,alphay;
  double N,Qx,Qy,Mx,My,Qo,Mo;
  double rate1,rate2,rate;
  double Lm; /*MEMORY LENGTH*/
  double bblength[2],bbfact[2],btlength[2],btfact[2]; /*BUCKLING LENGTH*/

  /*INITIAL*/
  Lm=L;
  *Ncr =0.0;
  *Qcrx=0.0;
  *Qcry=0.0;
  *Mcrx=0.0;
  *Mcry=0.0;

  /*MATERIAL*/
  if(sect.stype==STYPE_S)
  {
    poi=1.0/3.0;
    Fc=sect.sform.F;
    Ft=sect.sform.F;
    Fb=sect.sform.F;
    Fs=sect.sform.F;
    /*
    if(F==2400.0)      fu=4000.0;
    else if(F==3300.0) fu=5000.0;
    else               fu=4000.0;
    */
  }
  if(sect.stype==STYPE_GLASS)
  {
    poi=0.225;
    Fc=sect.wform.Fc;
    Ft=sect.wform.Ft;
    Fb=sect.wform.Fb;
    Fs=sect.wform.Fs;
  }
  if(sect.stype==STYPE_ACRYL)
  {
    poi=0.400;
    Fc=sect.wform.Fc;
    Ft=sect.wform.Ft;
    Fb=sect.wform.Fb;
    Fs=sect.wform.Fs;
  }

  G=E/2.0/(1.0+poi);

  /*FORM,STRESS*/
  if(sect.sform.H >= sect.sform.B)
  {
    B=sect.sform.B; /*[cm]*/
    D=sect.sform.H;

    if(iend==HEAD) N=st->N;
    if(iend==TAIL) N=-(st->N);
    Qx=fabs(st->Q[0]);
    Qy=fabs(st->Q[1]);
    Mx=fabs(st->M[0]);
    My=fabs(st->M[1]);
/*
    bblength[0]=sect.bblength[0];
    bblength[1]=sect.bblength[1];
    bbfact[0]=sect.bbfact[0];
    bbfact[1]=sect.bbfact[1];
    btlength[0]=sect.btlength[0];
    btlength[1]=sect.btlength[1];
    btfact[0]=sect.btfact[0];
    btfact[1]=sect.btfact[1];
*/
  }
  else
  {
    B=sect.sform.H;
    D=sect.sform.B;

    if(iend==HEAD) N=st->N;
    if(iend==TAIL) N=-(st->N);
    Qx=fabs(st->Q[1]);
    Qy=fabs(st->Q[0]);
    Mx=fabs(st->M[1]);
    My=fabs(st->M[0]);
/*
    bblength[0]=sect.bblength[1];
    bblength[1]=sect.bblength[0];
    bbfact[0]=sect.bbfact[1];
    bbfact[1]=sect.bbfact[0];
    btlength[0]=sect.btlength[1];
    btlength[1]=sect.btlength[0];
    btfact[0]=sect.btfact[1];
    btfact[1]=sect.btfact[0];
*/
  }
  /*
  fprintf(fout0,"B =%8.3f\n",B);
  fprintf(fout0,"D =%8.3f\n",D);
  fprintf(fout0,"N =%8.3f\n",N);
  fprintf(fout0,"Qx=%8.3f\n",Qx);
  fprintf(fout0,"Qy=%8.3f\n",Qy);
  fprintf(fout0,"Mx=%8.3f\n",Mx);
  fprintf(fout0,"My=%8.3f\n",My);
  */

  /*if(N<0.0) N=0.0;*/

  /*CONSTANTS*/
  A  =B*D;
  Ix =B*D*D*D/12.0;
  Iy =B*B*B*D/12.0;
  Iw =0.0;
  Jz =B*B*B*D/3.0;
  Zx =B*D*D/6.0;
  Zy =B*B*D/6.0;
  Zpx=B*D*D/4.0;
  Zpy=B*B*D/4.0;
  ix=sqrt(Ix/A);
  iy=sqrt(Iy/A);

  /*ALLOWABLE AXIAL STRESS*/

  /*WEAK AXIS*/
  fact=1.0;

  /*Lk=fact*L;*/ /*BUCKLING LENGTH [cm]*/

//  if(bblength[1]!=0.0) Lk=bblength[1]; /*BUCKLING LENGTH [cm]*/
//  else if(bbfact[1]!=0.0) Lk=bbfact[1]*Lm;
//  else Lk=fact*L;

  Lk=PI*sqrt(E*Iy/cNcr/1000);

  lamda=Lk/iy;
  Lamda=sqrt(PI*PI*E/0.6/Fc);
  nyu=1.5+(lamda*lamda)/(Lamda*Lamda)/1.5;
  if(nyu>2.17) nyu=2.17; /*kaza 20210822*/

  if(lamda<=Lamda) fc=(1.0-0.4*(lamda*lamda)/(Lamda*Lamda))*Fc/1000.0/nyu;
  else             fc=0.277*Fc/1000.0/(lamda*lamda)*(Lamda*Lamda);

  if(period==PSHORT) fc*=1.5;

  if(N>=0.0) Nay=fc*A; /*[tf]*/
  else
  {
    Nay=Ft/1000.0*A; if(period==PLONG) Nay/=1.5;
  }

  Na=Nay;

  /*STRONG AXIS*/
  fact=1.0;

  /*Lk=fact*L;*/ /*BUCKLING LENGTH [cm]*/

//  if(bblength[0]!=0.0) Lk=bblength[0]; /*BUCKLING LENGTH [cm]*/
//  else if(bbfact[0]!=0.0) Lk=bbfact[0]*Lm;
//  else Lk=fact*L;

  Lk=PI*sqrt(E*Ix/cNcr/1000);

  lamda=Lk/ix;
  Lamda=sqrt(PI*PI*E/0.6/Fc);
  nyu=1.5+(lamda*lamda)/(Lamda*Lamda)/1.5;
  if(nyu>2.17) nyu=2.17; /*kaza 20210822*/

  if(lamda<=Lamda) fc=(1.0-0.4*(lamda*lamda)/(Lamda*Lamda))*Fc/1000.0/nyu;
  else             fc=0.277*Fc/1000.0/(lamda*lamda)*(Lamda*Lamda);

  if(period==PSHORT) fc*=1.5;

  if(N>=0.0) Nax=fc*A; /*[tf]*/
  else
  {
    Nax=Ft/1000.0*A; if(period==PLONG) Nax/=1.5;
  }

  if(Na>Nax) Na=Nax;

  /*fprintf(fout0,"Lk/iy=%8.3f/%8.3f=%8.3f\n",Lk,iy,lamda);
  fprintf(fout0,"Na=%8.3f\n",Na);*/

  /*STRESS LIMITS*/
  if(N>=0.0) Np=Fc/1000.0*A;
  else       Np=Ft/1000.0*A;

  if(period==PLONG) Np/=1.5;

  /*STRONG AXIS*/
  fact=1.0; /*HINGE,HINGE*/
  /*Lk=fact*L;*/ /*BUCKLING LENGTH [cm]*/
//  if(bblength[0]!=0.0) Lk=bblength[0]; /*BUCKLING LENGTH [cm]*/
//  else if(bbfact[0]!=0.0) Lk=bbfact[0]*Lm;
//  else Lk=fact*L;

  Lk=PI*sqrt(E*Ix/cNcr/1000);

  Ncrx=PI*PI*E*Ix/Lk/Lk/1000.0; if(period==PLONG) Ncrx/=1.5;

  /*WEAK AXIS*/
  fact=1.0; /*HINGE,HINGE*/
  /*Lk=fact*L;*/ /*BUCKLING LENGTH [cm]*/
//  if(bblength[1]!=0.0) Lk=bblength[1]; /*BUCKLING LENGTH [cm]*/
//  else if(bbfact[1]!=0.0) Lk=bbfact[1]*Lm;
//  else Lk=fact*L;

  Lk=PI*sqrt(E*Iy/cNcr/1000);

  Ncry=PI*PI*E*Iy/Lk/Lk/1000.0;  if(period==PLONG) Ncry/=1.5;

  if(sect.stype==STYPE_S)     Qax=Fs/1000.0/sqrt(3.0)*A/1.5;
  if(sect.stype==STYPE_GLASS) Qax=Fs/1000.0*A/1.5;
  if(sect.stype==STYPE_ACRYL) Qax=Fs/1000.0*A/1.5;

  if(period==PLONG) Qax /=1.5;

  Qay=Qax;
  Myx=Zx*Fb/100000.0;            if(period==PLONG) Myx /=1.5;
  Myy=Zy*Fb/100000.0;            if(period==PLONG) Myy /=1.5;
  Mpx=Zpx*Fb/100000.0;           if(period==PLONG) Mpx /=1.5;
  Mpy=Zpy*Fb/100000.0;           if(period==PLONG) Mpy /=1.5;
  /*
  fprintf(fout0,"Np  =%8.3f\n",Np);
  fprintf(fout0,"Ncrx=%8.3f\n",Ncrx);
  fprintf(fout0,"Ncry=%8.3f\n",Ncry);
  fprintf(fout0,"Qax =%8.3f\n",Qax);
  fprintf(fout0,"Qay =%8.3f\n",Qay);
  fprintf(fout0,"Myx =%8.3f\n",Myx);
  fprintf(fout0,"Myy =%8.3f\n",Myy);
  fprintf(fout0,"Mpx =%8.3f\n",Mpx);
  fprintf(fout0,"Mpy =%8.3f\n",Mpy);
  */

  /*LATERAL LINEAR BUCKLING MOMENT*/
  g=0.0;
  k=0.0;
  /*ku=1.0/sqrt(2.0);*/
  /*kb=1.0/sqrt(2.0);*/

  /*q=0.0;*/
  M0M1=0.0;
  M2M1=1.0;
  gzai=0.283*(1+M2M1*M2M1)+0.434*M2M1+0.868*M0M1*(1+M2M1)+0.78*M0M1*M0M1;
  C1=1.0/sqrt(gzai);
  C2=0.405*M2M1/sqrt(gzai);
  C3=(0.5+0.5*M2M1+0.464*M0M1)/sqrt(gzai);

  /*STRONG AXIS*/
  Lk=PI*sqrt(E*Ix/cNcr/1000);
  ku=Lk/L;
  kb=Lk/L;
/*
  if(bblength[0]!=0.0) {Lk=bblength[0]; ku=Lk/L;}
  else if(bbfact[0]!=0.0) ku=bbfact[0];
  else ku=1.0;
  if(btlength[0]!=0.0) {Lk=btlength[0]; kb=Lk/L;}
  else if(btfact[0]!=0.0) kb=btfact[0];
  else kb=1.0/sqrt(2.0);
*/
  Mcrex=C1*PI*PI*E*Ix/(ku*L)/(ku*L)
       *((C2*g+C3*k)
         +sqrt((C2*g+C3*k)*(C2*g+C3*k)
               +Iw/Ix*(ku/kb)*(ku/kb)
               +G*Jz*(ku*L)*(ku*L)/(PI*PI*E*Ix)))/100000.0;

  /*WEAK AXIS*/
  Lk=PI*sqrt(E*Iy/cNcr/1000);
  ku=Lk/L;
  kb=Lk/L;
/*
  if(bblength[1]!=0.0) {Lk=bblength[1]; ku=Lk/L;}
  else if(bbfact[1]!=0.0) ku=bbfact[1];
  else ku=1.0;
  if(btlength[1]!=0.0) {Lk=btlength[1]; kb=Lk/L;}
  else if(btfact[1]!=0.0) kb=btfact[1];
  else kb=1.0/sqrt(2.0);
*/
  Mcrey=C1*PI*PI*E*Iy/(ku*L)/(ku*L)
       *((C2*g+C3*k)
         +sqrt((C2*g+C3*k)*(C2*g+C3*k)
               +Iw/Iy*(ku/kb)*(ku/kb)
               +G*Jz*(ku*L)*(ku*L)/(PI*PI*E*Iy)))/100000.0;

  if(Mcrex<=Mcrey) Mcre=Mcrex;
  else             Mcre=Mcrey;

  if(period==PLONG) Mcre/=1.5;

  /*fprintf(fout0,"Mcre=%8.3f\n",Mcre);*/

  /*LATERAL NONLINEAR BUCKLING MOMENT*/
  lamda=sqrt(Mpx/Mcre);
  fact=0.7+0.3*(1-0.7*lamda*lamda)/(0.61+0.3*M2M1+0.07*M2M1*M2M1);

  if(lamda>sqrt(1/0.7)) Mcrp=Mcre;
  else                  Mcrp=fact*Mpx;

  /*fprintf(fout0,"Mcrp=%8.3f\n",Mcrp);*/

  if(N>=0.0) /*COMPRESSION*/
  {
    /*P DELTA*/
    Kx=sqrt(N*1000.0/E/Ix);
    Ky=sqrt(N*1000.0/E/Iy);

    /*Cx:STRONG AXIS*/
    Lk=PI*sqrt(E*Ix/cNcr/1000);
    ku=Lk/L;
    kb=Lk/L;

//    if(bblength[0]!=0.0) {Lk=bblength[0]; ku=Lk/L;}
//    else if(bbfact[0]!=0.0) ku=bbfact[0];
//    else ku=1.0; /*?*/
//    if(btlength[0]!=0.0) {Lk=btlength[0]; kb=Lk/L;}
//    else if(btfact[0]!=0.0) kb=btfact[0];
//    else kb=1.0/sqrt(2.0); /*?*/


    /*cx=cos(Kx*L);*/
    cx=cos(Kx*Lk);

    /*if(M2M1>=cx) alphax=sqrt(1.0+M2M1*M2M1-2.0*M2M1*cos(Kx*L))/sin(Kx*L);*/
    if(N==0.0 || Mx==0.0) alphax=1.0;
    else if(M2M1>=cx) alphax=sqrt(1.0+M2M1*M2M1-2.0*M2M1*cos(Kx*Lk))/sin(Kx*Lk);
    else return 10.0;

    /*Cy:WEAK AXIS*/

    Lk=PI*sqrt(E*Iy/cNcr/1000);
    ku=Lk/L;
    kb=Lk/L;

//    if(bblength[1]!=0.0) {Lk=bblength[1]; ku=Lk/L;}
//    else if(bbfact[1]!=0.0) ku=bbfact[1];
//    else ku=1.0; /*?*/
//    if(btlength[1]!=0.0) {Lk=btlength[1]; kb=Lk/L;}
//    else if(btfact[1]!=0.0) kb=btfact[1];
//    else kb=1.0/sqrt(2.0); /*?*/


    /*cy=cos(Ky*L);*/
    cy=cos(Ky*Lk);

    /*if(M2M1>=cy) alphay=sqrt(1.0+M2M1*M2M1-2.0*M2M1*cos(Ky*L))/sin(Ky*L);*/
    if(N==0.0 || My==0.0) alphay=1.0;
    else if(M2M1>=cy) alphay=sqrt(1.0+M2M1*M2M1-2.0*M2M1*cos(Ky*Lk))/sin(Ky*Lk);
    else return 10.0;
  }
  else /*TENSION*/
  {
    alphax=1.0; alphay=1.0;
  }

  if(alphax<0.0 || alphay<0.0) return 10.0;

  /*fprintf(fout0,"Alpha X=%.5f Y=%.5f\n",alphax,alphay);*/

  /*CRITICAL STRESSES*/
  if(fabs(Na)<=Np) *Ncr=Na;
  else             *Ncr=Np;
  *Qcrx=Qax;
  *Qcry=Qay;
  if(Myx<=Mcre && Myx<=Mcrp) *Mcrx=Myx;
  else if(Mcre<=Mcrp)        *Mcrx=Mcre;
  else                       *Mcrx=Mcrp;
  *Mcry=Myy;

  /*RATE*/
  rate1=fabs(N)/(*Ncr)+alphax*Mx/(*Mcrx)+alphay*My/(*Mcry);
  rate2=Qx/(*Qcrx)+Qy/(*Qcry);

  rate=sqrt(rate1*rate1+rate2*rate2);

  /*
  fprintf(fout0,"\n");
  fprintf(fout0,"N/Ncr=%.3f/%.3f=%.3f\n",N,*Ncr,N/(*Ncr));
  fprintf(fout0,"Mx/Mcrx=%.3f/%.3f=%.3f\n",Mx,*Mcrx,Mx/(*Mcrx));
  fprintf(fout0,"My/Mcry=%.3f/%.3f=%.3f\n",My,*Mcry,My/(*Mcry));
  fprintf(fout0,"Qx/Qcrx=%.3f/%.3f=%.3f\n",Qx,*Qcrx,Qx/(*Qcrx));
  fprintf(fout0,"Qy/Qcry=%.3f/%.3f=%.3f\n",Qy,*Qcry,Qy/(*Qcry));
  fprintf(fout0,"RATE=%.5f\n",rate);
  */

  /*CHANGE AXIS*/
  if(sect.sform.H < sect.sform.B)
  {
    Qo=*Qcrx; *Qcrx=*Qcry; *Qcry=Qo;
    Mo=*Mcrx; *Mcrx=*Mcry; *Mcry=Mo;
  }

  return rate;
}/*allowablestressofflatbar_bc*/

double allowablestressofflatbaralumi(int period,
                                     struct section sect,
                                     double E,
                                     double L,       /*ELEMENT LENGTH[cm]*/
                                     struct stress *st,        /*[tf,tfm]*/
                                     int iend,
                                     double *Ncr,
                                     double *Qcrx,double *Qcry,
                                     double *Mcrx,double *Mcry)
/*FLAT BAR COLUMN.*/
{
  double poi,G,Fc,Ft,Fb,Fs/*,fu*/;
  double B,D,fact,Lk,A,Ix,Iy,Iw,Jz,Zx,Zy,Zpx,Zpy,ix,iy;
  double fc,Na,Nax,Nay,Np,Ncrx,Ncry,Qax,Qay,Myx,Myy;
  double lamda,ramda,eramda,pramda,nyu,Ny,Ne;
  double M2M1,Cb,Mem,Mym;
  double Mcre,Mcrex,Mcrey;
  double Kx,Ky,cx,cy,alphax,alphay;
  double N,Qx,Qy,Mx,My,Qo,Mo;


  double rate1,rate2,rate;
  double Lm; /*MEMORY LENGTH*/

  /*INITIAL*/
  Lm=L;
  *Ncr =0.0;
  *Qcrx=0.0;
  *Qcry=0.0;
  *Mcrx=0.0;
  *Mcry=0.0;

  /*MATERIAL*/
  poi=0.3;
  Fc=sect.sform.F;
  Ft=sect.sform.F;
  Fb=sect.sform.F;
  Fs=sect.sform.F;

  G=E/2.0/(1.0+poi);

  /*FORM,STRESS*/
  if(sect.sform.H >= sect.sform.B)
  {
    B=sect.sform.B; /*[cm]*/
    D=sect.sform.H;

    if(iend==HEAD) N=st->N;
    if(iend==TAIL) N=-(st->N);
    Qx=fabs(st->Q[0]);
    Qy=fabs(st->Q[1]);
    Mx=fabs(st->M[0]);
    My=fabs(st->M[1]);
  }
  else
  {
    B=sect.sform.H;
    D=sect.sform.B;

    if(iend==HEAD) N=st->N;
    if(iend==TAIL) N=-(st->N);
    Qx=fabs(st->Q[1]);
    Qy=fabs(st->Q[0]);
    Mx=fabs(st->M[1]);
    My=fabs(st->M[0]);
  }

  /*CONSTANTS*/
  A  =B*D;
  Ix =B*D*D*D/12.0;
  Iy =B*B*B*D/12.0;
  Iw =0.0;
  Jz =B*B*B*D/3.0;
  Zx =B*D*D/6.0;
  Zy =B*B*D/6.0;
  Zpx=B*D*D/4.0;
  Zpy=B*B*D/4.0;
  ix=sqrt(Ix/A);
  iy=sqrt(Iy/A);

/*ALLOWABLE AXIAL STRESS*/
/*WEAK AXIS*/
  fact=1.0;
  Lk=fact*L; /*BUCKLING LENGTH [cm]*/

if(!strcmp(prj,"Alumi") && sect.code==203) Lk=8.0;

  Ny=sect.sform.F*A;
  Ne=PI*PI*E*Iy/Lk/Lk;

  ramda=sqrt(Ny/Ne);
  eramda=1/sqrt(0.5);
  pramda=0.2;

  nyu=1.5+(ramda/eramda)*(ramda/eramda)/1.5;
  if(nyu>2.17) nyu=2.17; /*à¿ëSó¶kaza 20210822*/

  if(ramda<=pramda)
  {
    fc=sect.sform.F/nyu;
  }
  else if(pramda<ramda && ramda<=eramda)
  {
    fc=(1.0+0.5*(ramda-pramda)/(eramda-pramda))*sect.sform.F/nyu;
  }
  else
  {
    fc=sect.sform.F/ramda/ramda/nyu;
  }

  if(period==PSHORT) fc*=1.5;

  if(N>=0.0) Nay=fc/1000.0*A; /*[tf]*/
  else
  {
    Nay=Ft/1000.0*A; if(period==PLONG) Nay/=1.5;
  }

  Na=Nay;

/*STRONG AXIS*/
  fact=1.0;
  Lk=fact*L; /*BUCKLING LENGTH [cm]*/

if(!strcmp(prj,"Alumi") && sect.code==203) Lk=8.0;

  Ny=sect.sform.F*A;
  Ne=PI*PI*E*Ix/Lk/Lk;

  ramda=sqrt(Ny/Ne);
  eramda=1/sqrt(0.5);
  pramda=0.2;

  nyu=1.5+(ramda/eramda)*(ramda/eramda)/1.5;
  if(nyu>2.17) nyu=2.17; /*kaza 20210822*/

  if(ramda<=pramda)
  {
    fc=sect.sform.F/nyu;
  }
  else if(pramda<ramda && ramda<=eramda)
  {
    fc=(1.0+0.5*(ramda-pramda)/(eramda-pramda))*sect.sform.F/nyu;
  }
  else
  {
    fc=sect.sform.F/ramda/ramda/nyu;
  }

  if(period==PSHORT) fc*=1.5;

  if(N>=0.0) Nax=fc/1000.0*A; /*[tf]*/
  else
  {
    Nax=Ft/1000.0*A; if(period==PLONG) Nax/=1.5;
  }

  if(Na>Nax) Na=Nax;



  /*STRESS LIMITS*/
  if(N>=0.0) Np=Fc/1000.0*A;
  else       Np=Ft/1000.0*A;

  if(period==PLONG) Np/=1.5;


  if(sect.stype==STYPE_ALUMI) Qax=Fs/1000.0/sqrt(3.0)*A/1.5;
  if(period==PLONG) Qax /=1.5;

  Qay=Qax;
  Myx=Zx*Fb/100000.0;            if(period==PLONG) Myx /=1.5;
  Myy=Zy*Fb/100000.0;            if(period==PLONG) Myy /=1.5;

/*LATERAL LINEAR BUCKLING MOMENT*/
  M2M1=1.0;
  Cb=1.75+1.05*M2M1+0.3*M2M1*M2M1;

/*STRONG AXIS*/
if(!strcmp(prj,"Alumi") && sect.code==203) Lk=8.0;

  Mym=sect.sform.F*Zx;
  Mem=Cb*sqrt(PI*PI*E*Ix*G*Jz/Lk/Lk);

  ramda=sqrt(Mym/Mem);
  eramda=1/sqrt(0.5);
  pramda=0.6+0.3*M2M1;
  nyu=1.5+(ramda/eramda)*(ramda/eramda)/1.5;

  if(ramda<=pramda)
  {
    Mcrex=(sect.sform.F/nyu)*Zx/100000.0;
  }
  else if(pramda<ramda && ramda<=eramda)
  {
    Mcrex=((1.0+0.5*(ramda-pramda)/(eramda-pramda))*sect.sform.F/nyu)*Zx/100000.0;
  }
  else
  {
    Mcrex=(sect.sform.F/ramda/ramda/nyu)*Zx/100000.0;
  }

/*WEAK AXIS*/
if(!strcmp(prj,"Alumi") && sect.code==203) Lk=8.0;

  Mym=sect.sform.F*Zy;
  Mem=Cb*sqrt(PI*PI*E*Iy*G*Jz/Lk/Lk);

  ramda=sqrt(Mym/Mem);
  eramda=1/sqrt(0.5);
  pramda=0.6+0.3*M2M1;
  nyu=1.5+(ramda/eramda)*(ramda/eramda)/1.5;

  if(ramda<=pramda)
  {
    Mcrey=(sect.sform.F/nyu)*Zy/100000.0;
  }
  else if(pramda<ramda && ramda<=eramda)
  {
    Mcrey=((1.0+0.5*(ramda-pramda)/(eramda-pramda))*sect.sform.F/nyu)*Zy/100000.0;
  }
  else
  {
    Mcrey=(sect.sform.F/ramda/ramda/nyu)*Zy/100000.0;
  }

  if(Mcrex<=Mcrey) Mcre=Mcrex;
  else             Mcre=Mcrey;

  if(period==PLONG) Mcre/=1.5;

/*LATERAL NONLINEAR BUCKLING MOMENT*/
  if(N>=0.0) /*COMPRESSION*/
  {
    /*P DELTA*/
    Kx=sqrt(N*1000.0/E/Ix);
    Ky=sqrt(N*1000.0/E/Iy);

/*Cx:STRONG AXIS*/
if(!strcmp(prj,"Alumi") && sect.code==203) Lk=8.0;

    cx=cos(Kx*Lk);

    if(N==0.0 || Mx==0.0) alphax=1.0;
    else if(M2M1>=cx) alphax=sqrt(1.0+M2M1*M2M1-2.0*M2M1*cos(Kx*Lk))/sin(Kx*Lk);
    else return 10.0;

/*Cy:WEAK AXIS*/
if(!strcmp(prj,"Alumi") && sect.code==203) Lk=8.0;

    cy=cos(Ky*Lk);

    if(N==0.0 || My==0.0) alphay=1.0;
    else if(M2M1>=cy) alphay=sqrt(1.0+M2M1*M2M1-2.0*M2M1*cos(Ky*Lk))/sin(Ky*Lk);
    else return 10.0;
  }
  else /*TENSION*/
  {
    alphax=1.0; alphay=1.0;
  }

  if(alphax<0.0 || alphay<0.0) return 10.0;


/*CRITICAL STRESSES*/
  if(fabs(Na)<=Np) *Ncr=Na;
  else             *Ncr=Np;
  *Qcrx=Qax;
  *Qcry=Qay;
  if(Myx<=Mcre) *Mcrx=Myx;
  else          *Mcrx=Mcre;
  *Mcry=Myy;

  /*RATE*/
  rate1=fabs(N)/(*Ncr)+alphax*Mx/(*Mcrx)+alphay*My/(*Mcry);
  rate2=Qx/(*Qcrx)+Qy/(*Qcry);
/*
    fprintf(fout0,"N/Ncr=%.3f/%.3f\n",N,*Ncr);
    fprintf(fout0,"Mx/Mcrx=%.3f/%.3f\n",Mx,*Mcrx);
    fprintf(fout0,"My/Mcry=%.3f/%.3f\n",My,*Mcry);
    fprintf(fout0,"Qx/Qcrx=%.3f/%.3f\n",Qx,*Qcrx);
    fprintf(fout0,"Qy/Qcry=%.3f/%.3f\n",Qy,*Qcry);
    fprintf(fout0,"RATE=%.5f\n",rate);
*/
  rate=sqrt(rate1*rate1+rate2*rate2);

  /*CHANGE AXIS*/
  if(sect.sform.H < sect.sform.B)
  {
    Qo=*Qcrx; *Qcrx=*Qcry; *Qcry=Qo;
    Mo=*Mcrx; *Mcrx=*Mcry; *Mcry=Mo;
  }

  return rate;
}/*allowablestressofflatbaralumi*/

double allowablebendingofwood(int period,
                              struct section sect,
                              double lkxx,double lkyy,       /*BUCKLING LENGTH[cm]*/
                              int axis,                    /*0:x 1:y*/
                              double Nd)                     /*[kgf]*/
/*RETURN:ALLOWABLE BENDING OF WOOD.*/
/*N=+:TENSION -:COMPRESSION*/
/*WOOD TYPE:PLATE*/
{
  double A,Ixx,Iyy,Zxx,Zyy,ixx,iyy;
  double H,B;
  double lk,i,lamda,lamxx,lamyy;
  double fca,fta,fba;
  double Na=0.0,Ma=0.0;

  /*COEFFICIENTS*/
  H=sect.sform.H;
  B=sect.sform.B;

  A=B*H;
  Ixx=B*H*H*H/12.0;
  Iyy=B*B*B*H/12.0;
  Zxx=B*H*H/6.0;
  Zyy=B*B*H/6.0;
  ixx=sqrt(Ixx/A);
  iyy=sqrt(Iyy/A);

  /*if(ixx<iyy) i=ixx;*/
  /*else        i=iyy;*/

  /*ALLOWABLE STRESS*/
  if(period==PLONG)
  {
    fca=1.1/3.0*sect.wform.Fc; /*[kgf/cm2]*/
    fta=1.1/3.0*sect.wform.Ft;
    fba=1.1/3.0*sect.wform.Fb;
  }
  if(period==PSHORT)
  {
    fca=2.0/3.0*sect.wform.Fc; /*[kgf/cm2]*/
    fta=2.0/3.0*sect.wform.Ft;
    fba=2.0/3.0*sect.wform.Fb;
  }

  /*ALLOWABLE AXIAL FORCE*/
if(!strcmp(prj,"huru") && sect.code==202) lk=850.0;
if(!strcmp(prj,"huru") && sect.code==204) lk=850.0;

if(!strcmp(prj,"waka") && sect.code==501) lk=0.0;

if(!strcmp(prj,"noa")  && sect.code==201) ixx=iyy;
if(!strcmp(prj,"noa")  && sect.code==501) iyy=ixx;
if(!strcmp(prj,"noa")  && sect.code==502) iyy=ixx;

/*if(!strcmp(prj,"hiro") && sect.code==201) i=ixx*sqrt(2.0);*/
if(!strcmp(prj,"hiro") && sect.code==201) lk=0.0;
if(!strcmp(prj,"hiro") && sect.code==204) lk=0.0;
if(!strcmp(prj,"hiro") && sect.code==501) iyy=ixx;

  /*if(ixx<iyy) i=ixx;*/
  /*else        i=iyy;*/
  /*lamda=lk/i;*/

  if(sect.bblength[0]!=0.0) lk=sect.bblength[0];
  else if(sect.bbfact[0]!=0.0) lk=lkxx*sect.bbfact[0];
  else lk=lkxx;

  lamxx=lk/ixx;

  if(sect.bblength[1]!=0.0) lk=sect.bblength[1];
  else if(sect.bbfact[1]!=0.0) lk=lkyy*sect.bbfact[1];
  else lk=lkyy;

  lamyy=lk/iyy;

  if(lamxx>lamyy) lamda=lamxx;
  else            lamda=lamyy;

  if( 30.0<lamda && lamda<=100.0) fca*=(1.3-0.01*lamda);
  if(100.0<lamda)                 fca*=(0.3/(lamda/100.0)/(lamda/100.0));

  if(Nd< 0.0) Na=-fca*A; /*COMPRESSION [kgf]*/
  if(Nd>=0.0) Na= fta*A; /*TENSION     [kgf]*/

  if((Nd/Na)>=1.0)
  {
    fprintf(fout0,"ERROR.N > Na=%.3f\n",Na);
    return 0.0;
  }

  if(axis==SX) Ma=(1.0-Nd/Na)*fba*Zxx;                    /*[kgfcm]*/
  if(axis==SY) Ma=(1.0-Nd/Na)*fba*Zyy;                    /*[kgfcm]*/

  /*fprintf(fout0,"Ma=%.3f\n",Ma);*/

  return Ma;
}/*allowablebendingofwood*/

double allowablestressofdabo(int period,
                             struct section sect,
                             struct stress *st,          /*[tf,tfm]*/
                             int iend,              /*0:HEAD 1:TAIL*/
                             double *Na,double *Qa,double *Ma)
/*RETURN:ALLOWABLE BENDING OF WOOD.*/
/*N=+:TENSION -:COMPRESSION*/
/*WOOD TYPE:PLATE*/
{
  double A,H,B,Aw,Zw;
  double fca,fta,fba,fsa;
  double N,Qx,Qy,Mx,My;
  double dt,woodsize=35.0; /*[cm]*/
  double Nta,Nda,Mda; /*ALLOWABLE,ULTIMATE STRESS OF DABO [tf]*/
  double rate21,rate22,rate2,rate;

  /*COEFFICIENTS*/
  H=sect.sform.H;
  B=sect.sform.B;

  A=B*H;

  Aw=woodsize*woodsize;
  Zw=woodsize*woodsize*woodsize/6.0;

  /*STRESS*/
  if(iend==HEAD) N=st->N;
  if(iend==TAIL) N=-(st->N);
  Qx=fabs(st->Q[0]);
  Qy=fabs(st->Q[1]);
  Mx=fabs(st->M[0]);
  My=fabs(st->M[1]);

  *Na=0.0;
  *Qa=0.0;
  *Ma=0.0;

  /*ALLOWABLE STRESS*/
  if(period==PLONG)
  {
    Nta=0.750*A/(PI*3.0*3.0);  /*0.750 tf FOR D=60mm*/
    fca=1.1/3.0*sect.wform.Fc/1000.0;     /*[tf/cm2]*/
    fta=1.1/3.0*sect.wform.Ft/1000.0;     /*[tf/cm2]*/
    fba=1.1/3.0*sect.wform.Fb/1000.0;     /*[tf/cm2]*/
    fsa=1.1/3.0*sect.wform.Fs/1000.0;     /*[tf/cm2]*/
  }
  if(period==PSHORT)
  {
    Nta=1.500*A/(PI*3.0*3.0);  /*1.500 tf FOR D=60mm*/
    fca=2.0/3.0*sect.wform.Fc/1000.0;     /*[tf/cm2]*/
    fta=2.0/3.0*sect.wform.Ft/1000.0;     /*[tf/cm2]*/
    fba=2.0/3.0*sect.wform.Fb/1000.0;     /*[tf/cm2]*/
    fsa=2.0/3.0*sect.wform.Fs/1000.0;     /*[tf/cm2]*/
  }

  /*ALLOWABLE SHEAR*/
  *Qa=fsa*A; /*[tf]*/
  dt=(woodsize/2.0-(2.0*Nta/woodsize/fba)/3.0)/100.0; /*[m]*/
  *Ma=Nta*dt; /*[tfm]*/

  Nda=fabs(fta*A); /*[tf]*/
  Mda=Nda*woodsize/100.0/3.0; /*[tfm]*/

  /*ALLOWABLE AXIAL FORCE, BENDING MOMENT*/
  if(N>=0.0)
  {
    *Na=fca*woodsize*woodsize;   /*COMPRESSION [tf]*/

    rate=fabs(N/(*Na))+Mx/(*Ma)+My/(*Ma);

    rate21=Mx/Mda+My/Mda;
    rate22=Qx/(*Qa)+Qy/(*Qa);
    rate2=sqrt(rate21*rate21+rate22*rate22);

    if(rate<rate2) rate=rate2;
  }
  if(N<0.0)
  {
    *Na=-Nta;                    /*TENSION     [tf]*/

    rate=fabs(N/(*Na))+Mx/(*Ma)+My/(*Ma);

    rate=fabs(N/Nda)+Mx/Mda+My/Mda;
    rate22=Qx/(*Qa)+Qy/(*Qa);
    rate2=sqrt(rate21*rate21+rate22*rate22);

    if(rate<rate2) rate=rate2;
  }

  /*SEARCH TENSION.*/
  /*
  if((Mx+My)>(N*(woodsize/100.0/2.0-N/(fba*10000.0)/(woodsize/100.0)/1.5)))
  {
    fprintf(fout0,"Mmax=%.3f\n",
            (N*(woodsize/100.0/2.0-N/(fba*10000.0)/(woodsize/100.0)/1.5)));

    rate=10.0;
  }
  else rate=0.0;
  */

  return rate;
}/*allowablestressofdabo*/

double ultimatebendingofsrc(struct element elem,
                            int axis,                     /*0:x 1:y*/
                            double Nd,                      /*[kgf]*/
                            double *sN,double *sM,
                            double *rN,double *rM,
                            double *cN,double *cM)
/*RETURN:ULTIMATE BENDING OF RC,SRC COLUMN.*/
/*       PLASTIC BENDING OF S COLUMN.*/
/*N=+:TENSION -:COMPRESSION*/
{
  char str[1024],non[256];
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

  double rAt,rAb,rYtg,reinrate;

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

    addnewbound(bounds,&nbound,sYi[i]);
    addnewbound(bounds,&nbound,sYj[i]);

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

    addnewbound(bounds,&nbound,cYi[i]);
    addnewbound(bounds,&nbound,cYj[i]);

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

    /*if(rYi[i]<=0.0) rAi[i]=elem.sect->rein[i].area;
    else            rAi[i]=0.0;*/

    Ar +=rAi[i];
    rAY+=rAi[i]*rYi[i];

    if(rYi[i]<0.0)
    {
      rAt +=rAi[i];
      rYtg+=rAi[i]*rYi[i];
    }

    if(addnewbound(bounds,&nbound,rYi[i]))
    {
      bounds[nbound]=rYi[i]; /*ADD SAME BOUND*/
      nbound++;
    }

    if(rYc<rYi[i]) rYc=rYi[i];
    if(rYt>rYi[i]) rYt=rYi[i];
  }
  if(rAt!=0.0) rYtg/=rAt;

  sortdouble(bounds,nbound); /*SORT BOUNDARIES.*/

  /*MATERIALS*/
  if(elem.sect->stype==STYPE_S || elem.sect->stype==STYPE_SRC)
  {
    sftu=elem.sect->srect[0].F;
    sfcu=-sftu;
  }

  if(elem.sect->stype==STYPE_RC || elem.sect->stype==STYPE_SRC)
  {
    rftu=elem.sect->rein[0].F;
    rfcu=-rftu;
  }

  if(elem.sect->stype==STYPE_RC)
  {
    cfcu=-0.85*(elem.sect->crect[0].F);
  }
  else if(elem.sect->stype==STYPE_SRC)
  {
    cfcu=-(elem.sect->crect[0].F)*(0.85-2.5*(As/2.0/Ac)); /*"SRC STANDARD"*/
  }

  /*
  if(fout0!=NULL)
  {
    fprintf(fout0,"As=%.3f sYc=%.3f sYt=%.3f[cm2,cm]\n",As,sYc,sYt);
    fprintf(fout0,"Ar=%.3f rYc=%.3f rYt=%.3f[cm2,cm]\n",Ar,rYc,rYt);
    fprintf(fout0,"Ac=%.3f cYc=%.3f cYt=%.3f[cm2,cm]\n",Ac,cYc,cYt);
  }
  */

  for(i=0;i<nbound;i++) /*BOUNDARY DEFINITION.*/
  {
    Yn=bounds[i];

    if(i<(nbound-1) && bounds[i]==bounds[i+1]) reinrate=1.0;
    else if(i>=1 && bounds[i-1]==bounds[i])    reinrate=0.0;
    else reinrate=0.0;

    boundaryvalueultimate(elem,Yn,sBi,sDi,cBi,cDi,
                                  sYi,sYj,cYi,cYj,
                                  sftu,sfcu,
                                  rftu,rfcu,
                                  0.0,cfcu,
                                  rAi,rYi,
                                  sN,sM,rN,rM,cN,cM,
                                  &N[i],&M[i],
                                  reinrate);
  }

  /*sprintf(str,"\0");
  for(i=0;i<nbound;i++)
  {
    sprintf(non,"BOUND%2d:Yn=%11.3f N=%12.3f M=%12.3f\n",
            i+1,bounds[i],N[i],M[i]);
    strcat(str,non);
  }
  MessageBox(NULL,str,"Ultimate Bending",MB_OK);*/

  /*RANGE DETERMINATION.*/
  if(N[0]<=Nd)
  {
    if(N[0]<Nd)
    {
      if(fout0!=NULL) fprintf(fout0,"ERROR:N OVERFLOW.\n");
    }

    Nover=Nd-N[0];
    return 0.0;
  }
  else if(Nd<=N[nbound-1])
  {
    if(Nd<N[nbound-1])
    {
      if(fout0!=NULL) fprintf(fout0,"ERROR:N OVERFLOW.\n");
    }

    Nover=Nd-N[nbound-1];
    return 0.0;
  }
  else Nover=0.0;

  ibound=0;
  while(N[ibound+1]>=Nd) ibound++;

  if(Nd==N[ibound])
  {
    yn=bounds[ibound];
    reinrate=0.0; /*LOWER*/
  }
  else if(bounds[ibound]==bounds[ibound+1]) /*ON REIN*/
  {
    yn=bounds[ibound];

    rAb=0.0; /*REINS ON BOUND*/
    for(j=0;j<(elem.sect->nrein);j++)
    {
      if(rYi[j]==bounds[ibound]) rAb+=rAi[j];
    }

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
        bunsi-=sftu*sBi[j]*(yn-sYj[j])
              +sfcu*sBi[j]*(sYi[j]-yn);
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
        bunsi-=cfcu*cBi[j]*(cYi[j]-yn);
      }
    }
    for(j=0;j<(elem.sect->nrein);j++)
    {
      if(rYi[j]<bounds[ibound+1]) /*TENSED*/
      {
        bunsi-=rftu*rAi[j];
      }
      else if(rYi[j]>bounds[ibound]) /*COMPRESSED*/
      {
        bunsi-=rfcu*rAi[j];
      }
      else if(rYi[j]==bounds[ibound]) /*ON BOUND*/
      {
        bunsi-=rfcu*rAi[j];
      }
    }

    if(rAb!=0.0) reinrate=bunsi/((rftu-rfcu)*rAb);
  }
  else /*Nj < Nd < Ni*/
  {
    reinrate=0.0;

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

    if(bunbo!=0.0) yn=bunsi/bunbo;
    else
    {
      if(fout0!=NULL)
      {
        fprintf(fout0,"ERROR:NEUTRAL AXIS OVERFLOW.BUNBO=%.3f\n",bunbo);
      }
      return 0.0;
    }
  }

  if(yn<bounds[ibound+1] || bounds[ibound]<yn)
  {
    if(fout0!=NULL)
    {
      fprintf(fout0,"ERROR:NEUTRAL AXIS OVERFLOW.%.3f<%.3f<%.3f[cm]\n",
              bounds[ibound+1],yn,bounds[ibound]);
    }

    if(elem.sect->stype==STYPE_RC)
    {
      Mu=0.9*rAt*rftu*(cYc-rYtg); /*CENTER STANDARD*/
      return Mu;
    }
    else return 0.0;
  }

  Yg=boundaryvalueultimate(elem,yn,sBi,sDi,cBi,cDi,
                                   sYi,sYj,cYi,cYj,
                                   sftu,sfcu,
                                   rftu,rfcu,
                                   0.0,cfcu,
                                   rAi,rYi,
                                   sN,sM,rN,rM,cN,cM,
                                   &Ng,&Mg,
                                   reinrate);
  Mu=Mg;

  /*if(fout0!=NULL)
  {
    fprintf(fout0,"SRCé≤óÕ srcN=%11.3f[tf]",-Ng/1000.0);
    fprintf(fout0," èdêS:%11.3f",Yg);
    fprintf(fout0," íÜóßé≤=%11.3f",yn);
    fprintf(fout0," èIã«íl srcMu=%11.3f[tfm]\n",Mu/100000.0);

    fprintf(fout0,"%11.3f %11.3f 0.0",Mu/100000.0,-Nd/1000.0);
  }*/

  return Mu;
}/*ultimatebendingofsrc*/

int addnewbound(double bounds[],int *nbound,double newbound)
/*ADD NEW BOUNDARY IF NOT EXISTENT.*/
{
  int i,ifind;

  ifind=0;
  for(i=0;i<(*nbound);i++)
  {
    if(bounds[i]==newbound)
    {
      ifind=1;
      break;
    }
  }

  if(!ifind)
  {
    bounds[*nbound]=newbound;
    (*nbound)++;
    return 1; /*ADDED*/
  }
  else return 0; /*NOT ADDED*/

}/*addnewbound*/

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

  /*if(fout0!=NULL)
  {
    fprintf(fout0,"Yn=%.5f[cm]\n",Yn);
    fprintf(fout0,"gN=S%15.3f+R%15.3f+C%15.3f=%15.3f[kgf]\n",
            sN,rN,cN,*gN);
    fprintf(fout0,"gM=S%15.3f+R%15.3f+C%15.3f=%15.3f[kgfcm]\n",
            sM,rM,cM,*gM);
  }*/

  return gY;
}/*boundaryvalue*/

double boundaryvalueultimate(struct element elem,/*struct materials m,*/
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
                             double *gN,double *gM,
                             double reinrate)
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
    else if(rYi[i]==Yn) /*ON BOUND*/
    {
      rNi=(1.0-reinrate)*rfcu*rAi[i];
      *rN+=rNi;
      *rM+=rNi*(gY-rYi[i]);

      rNi=reinrate*rftu*rAi[i];
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

  /*if(fout0!=NULL)
  {
    fprintf(fout0,"Yn=%.5f[cm] ",Yn);
    fprintf(fout0,"gN=S%15.3f+R%15.3f+C%15.3f=%15.3f[kgf]\n",
            *sN,*rN,*cN,*gN);
    fprintf(fout0,"gM=S%15.3f+R%15.3f+C%15.3f=%15.3f[kgfcm]\n",
            *sM,*rM,*cM,*gM);
  }*/

  return gY;
}/*boundaryvalueultimate*/

double allowableshearofrclong(struct element elem,
                              /*struct materials m,*/
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
  double wft,cfs,f1,f2; /*ALLOWABLE STRESSES*/
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

  if(Qd==0.0 || dc==0.0) ac=2.0;
  else
  {
    ac=4.0/(fabs(Md/Qd/dc)+1.0);
    if(ac<1.0) ac=1.0;
    if(ac>2.0) ac=2.0;
    if(elem.sect->etype==COLUMN && ac>1.5) ac=1.5; /*RC2010*/
  }

  /*MATERIALS*/
  f1=elem.sect->crect[0].F/30.0;
  f2=5.0+elem.sect->crect[0].F/100.0;
  if(f1<=f2) cfs=f1;
  else       cfs=f2;

  if(elem.sect->srein[1-axis].n>=1)
  {
    if(elem.sect->srein[1-axis].F==3000.0)      wft=2000.0;
    else if(elem.sect->srein[1-axis].F==3500.0) wft=2200.0;

    pw=elem.sect->srein[1-axis].area
      *elem.sect->srein[1-axis].n
      /elem.sect->srein[1-axis].pitch
      /cBi[0];
  }
  else
  {
    if(elem.sect->wF==3000.0)      wft=2000.0;
    else if(elem.sect->wF==3500.0) wft=2200.0;

    pw=elem.sect->shearrein[1-axis];
  }

  /*if(pw>0.012) pw=0.012;*/
  if(pw>0.006) pw=0.006; /*RC2010*/

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
  if(fout0!=NULL)
  {
    fprintf(fout0,"Ac=%.3f cYc=%.3f cYt=%.3f[cm2,cm]\n",Ac,cYc,cYt);
    fprintf(fout0,"rYt=%.3f[cm] pw=%.3f\n",rYt,pw);
    fprintf(fout0,"rcQa=%11.3f[tf]",fabs(Qa)/1000.0);
  }
  */

  return Qa;
}/*allowableshearofrclong*/

double allowableshearofrcshort(struct element elem,
                               /*struct materials m,*/
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
  double cfs,wft,f1,f2; /*ALLOWABLE STRESSES*/
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

  if(Qd==0.0 || dc==0.0) ac=2.0;
  else
  {
    ac=4.0/(fabs(Md/Qd/dc)+1.0);
    if(ac<1.0) ac=1.0;
    if(ac>2.0) ac=2.0;
    if(elem.sect->etype==COLUMN && ac>1.0) ac=1.0; /*RC2010*/
  }

  /*MATERIALS*/
  f1=elem.sect->crect[0].F/30.0;
  f2=5.0+elem.sect->crect[0].F/100.0;
  if(f1<=f2) cfs=1.5*f1;
  else       cfs=1.5*f2;

  if(elem.sect->srein[1-axis].n>=1)
  {
    if(elem.sect->srein[1-axis].F==3000.0)      wft=3000.0;
    else if(elem.sect->srein[1-axis].F==3500.0) wft=3500.0;

    pw=elem.sect->srein[1-axis].area
      *elem.sect->srein[1-axis].n
      /elem.sect->srein[1-axis].pitch
      /cBi[0];
  }
  else
  {
    if(elem.sect->wF==3000.0)      wft=3000.0;
    else if(elem.sect->wF==3500.0) wft=3500.0;

    pw=elem.sect->shearrein[1-axis];
  }

  if(pw>0.012) pw=0.012;

  /*if(fout0!=NULL)
  {
    fprintf(fout0,"dc=%.3f cfs=%.3f wft=%.3f pw=%.5f\n",
          dc,cfs,wft,pw);
  }*/

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
  if(fout0!=NULL)
  {
    fprintf(fout0,"Ac=%.3f cYc=%.3f cYt=%.3f[cm2,cm]\n",Ac,cYc,cYt);
    fprintf(fout0,"rYt=%.3f[cm] pw=%.3f\n",rYt,pw);
    fprintf(fout0,"rcQa=%11.3f[tf]",fabs(Qa)/1000.0);
  }
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
  if(Qd==0.0) ac=1.0;
  else
  {
    ac=fabs(Md/Qd/dc);
    if(ac<1.0) ac=1.0;
    if(ac>3.0) ac=3.0;
  }

  Yg/=gA; /*CENTER OF GRAVITY.*/
  at=0.0;
  for(i=0;i<(elem.sect->nrein);i++)
  {
    if(rYi[i]<Yg) at+=(elem.sect->rein[i].area); /*TENSED REIN.*/
  }
  pt=at/eAc*100.0; /*RATE OF TENSED REIN [%].*/

  Fc=fabs(/*m.Fc*/elem.sect->crect[0].F);

  if(elem.sect->srein[1-axis].n>=1)
  {
    wfp=fabs(elem.sect->srein[1-axis].F);

    pw=elem.sect->srein[1-axis].area
      *elem.sect->srein[1-axis].n
      /elem.sect->srein[1-axis].pitch
      /cBi[0];
  }
  else
  {
    wfp=fabs(elem.sect->wF);

    pw=elem.sect->shearrein[1-axis];
  }

  if(pw>0.012) pw=0.012;

  if((0.9+s0/250.0)<=0.0) return 0.0;

  bQu=((0.068*pow(pt,0.23)*(Fc+180.0))/(ac+0.12)+2.7*sqrt(pw*wfp))
     *7.0/8.0*eAc;                    /*"CENTER STANDARD"(ït1.7.2b)*/
  rcQu=(0.9+s0/250.0)*bQu;            /*"CENTER STANDARD"(ït1.7.4b)*/

  /*if(fout0!=NULL)
  {
    fprintf(fout0,"Nd=%.3f Qd=%.3f Md=%.3f\n",Nd,Qd,Md);
    fprintf(fout0,"d=%.3f[cm] pw=%.5f\n",dc,pw);
    fprintf(fout0,"rcQu=%11.3f[tf]\n",rcQu/1000.0);
  }*/

  return rcQu; /*[kgf]*/
}/*ultimateshearofrc*/

double ultimatetorsionofrc(struct section sect)
/*ULTIMATE TORSION OF RC.*/
/*FOR ONLY ONE CONCRETE RECTANGLE.*/
/*WITH NO N,Q,M.*/
{
  char str[256],s[80];
  int i;
  int nw1,nw2,nw;
  double B,D,Be,De; /*WIDTH,DEPTH OF CONCRETES*/
  double f1,f2,cfs,rft,wft1,wft2,wft; /*ALLOWABLE STRESSES*/
  double Ar,eAc,Rr,Rrmin,Rwmin,Xmax,Xmin,Ymax,Ymin;
  double pw1,pw2,pw,aw,pitch;
  double Mt1,Mt2,Mt3,Mtu;

  Ar=0.0; /*REIN AREA TOTAL[cm2]*/
  for(i=0;i<(sect.nrein);i++)
  {
    Ar+=sect.rein[i].area;

    Rr=sqrt(sect.rein[i].area/PI); /*RADIUS OF REIN[cm]*/
    if(i==0) Rrmin=Rr;
    else if(Rrmin>Rr) Rrmin=Rr;

    if(i==0)
    {
      Xmax=sect.rein[i].x;
      Xmin=sect.rein[i].x;
      Ymax=sect.rein[i].y;
      Ymin=sect.rein[i].y;
    }
    else
    {
      if(Xmax<sect.rein[i].x) Xmax=sect.rein[i].x;
      if(Xmin>sect.rein[i].x) Xmin=sect.rein[i].x;
      if(Ymax<sect.rein[i].y) Ymax=sect.rein[i].y;
      if(Ymin>sect.rein[i].y) Ymin=sect.rein[i].y;
    }
  }

  /*RADIUS OF HOOP*/
  if(sect.srein[0].area<=sect.srein[1].area)
  {
    Rwmin=sqrt(sect.srein[0].area/PI);
  }
  else Rwmin=sqrt(sect.srein[1].area/PI);

  /*CONCRETE AREA[cm2]*/
  B=fabs(sect.crect[0].right-sect.crect[0].left);
  D=fabs(sect.crect[0].top-sect.crect[0].bottom);

  Be=Xmax-Xmin+2.0*Rrmin+2.0*Rwmin; /*INSIDE HOOP*/
  De=Ymax-Ymin+2.0*Rrmin+2.0*Rwmin;
  eAc=Be*De;

  /*MATERIALS SHORT*/
  f1=sect.crect[0].F/30.0;
  f2=5.0+sect.crect[0].F/100.0;
  if(f1<=f2) cfs=1.5*f1;
  else       cfs=1.5*f2;

  if(sect.rein[0].F==3000.0)      rft=3000.0;
  else if(sect.rein[0].F==3500.0) rft=3500.0;

  /*MODIFY HOOP*/
  if(sect.srein[0].F==3000.0)      wft1=3000.0;
  else if(sect.srein[0].F==3500.0) wft1=3500.0;
  if(sect.srein[0].n>=2) nw1=2;
  else                   nw1=sect.srein[0].n;
  pw1=sect.srein[0].area*(double)nw1/sect.srein[0].pitch/D;

  if(sect.srein[1].F==3000.0)      wft2=3000.0;
  else if(sect.srein[1].F==3500.0) wft2=3500.0;
  if(sect.srein[1].n>=2) nw2=2;
  else                   nw2=sect.srein[1].n;
  pw2=sect.srein[1].area*(double)nw2/sect.srein[1].pitch/B;

  if((pw1*wft1)<=(pw2*wft2))
  {
    pw=pw1;
    aw=sect.srein[0].area;
    nw=nw1;
    wft=wft1;
    pitch=sect.srein[0].pitch;
    if(pw>0.012)
    {
      pw=0.012;
      aw=pw/sect.srein[0].pitch/D;
      nw=1;
    }
  }
  else
  {
    pw=pw2;
    aw=sect.srein[1].area;
    nw=nw2;
    wft=wft2;
    pitch=sect.srein[1].pitch;
    if(pw>0.012)
    {
      pw=0.012;
      aw=pw/sect.srein[1].pitch/B;
      nw=1;
    }
  }

  Mt1=B*B*D*cfs*4.0/3.0; /*BY CONCRETE*/

  Mt2=(double)nw*aw/pitch*eAc*wft; /*BY HOOP*/

  Mt3=2.0*Ar*eAc*rft/(2.0*(Be+De)); /*BY REINS*/

  if(Mt1<=Mt2 && Mt1<=Mt3) Mtu=Mt1;
  else if(Mt2<=Mt3)        Mtu=Mt2;
  else                     Mtu=Mt3;
 
  sprintf(str,"\0");
  sprintf(s,"B=%.3f D=%.3f Be=%.3f De=%.3f [cm]\n",B,D,Be,De);
  strcat(str,s);
  sprintf(s,"Mt1=%.3f Mt2=%.3f Mt3=%.3f Mtu=%.3f [kgfcm]\n",
          Mt1/100000.0,Mt2/100000.0,Mt3/100000.0,Mtu/100000.0);
  strcat(str,s);
  MessageBox(NULL,str,"Torsion",MB_OK);

  return Mtu; /*[kgfcm]*/
}/*ultimatetorsionofrc*/

double allowableshearofsrclong(struct element elem,
                               /*struct materials m,*/
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
  double wft,cfs,f1,f2; /*ALLOWABLE STRESSES*/
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

  /*MATERIALS*/
  f1=elem.sect->crect[0].F/30.0;
  f2=5.0+elem.sect->crect[0].F/100.0;
  if(f1<=f2) cfs=f1;
  else       cfs=f2;

  if(elem.sect->srein[1-axis].n>=1)
  {
    if(elem.sect->srein[1-axis].F==3000.0)      wft=2000.0;
    else if(elem.sect->srein[1-axis].F==3500.0) wft=2200.0;

    /*pw=elem.sect->srein[1-axis].area
      *elem.sect->srein[1-axis].n
      /elem.sect->srein[1-axis].pitch
      /cBi[0];*/
  }
  else
  {
    if(elem.sect->wF==3000.0)      wft=2000.0;
    else if(elem.sect->wF==3500.0) wft=2200.0;

    /*pw=elem.sect->shearrein[1-axis];*/
  }

  /*if(pw>0.006) pw=0.006;*/

  rcQa=7.0/8.0*ac*eAc*cfs; /*"SRC STANDARD"(58)*/

  if(haxis==HSTRONG)
  {
    sQa=allowultimshearofhstrong(PLONG,elem.sect->srect[0].F,
                                 *(elem.sect),axis);
  }
  if(haxis==HWEAK)
  {
    sQa=allowultimshearofhweak(PLONG,elem.sect->srect[0].F,
                               *(elem.sect),axis);
  }
  srcQa=sQa+rcQa;

  /*
  if(fout0!=NULL)
  {
    fprintf(fout0,"Ac=%.3f cYc=%.3f cYt=%.3f[cm2,cm]\n",Ac,cYc,cYt);
    fprintf(fout0,"rYt=%.3f[cm]\n",rYt);
    fprintf(fout0,"sQa=%.3f[tf] rcQa=%.3f[tf]\n",
            sQa/1000.0,rcQa/1000.0);
    fprintf(fout0,"srcQa=%11.3f[tf]",srcQa/1000.0);
  }
  */

  return srcQa;
}/*allowableshearofsrclong*/

double allowableshearofsrcshort(struct element elem,
                                /*struct materials m,*/
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
  double wft,cfs,f1,f2; /*ALLOWABLE STRESSES*/
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

  /*MATERIALS*/
  f1=elem.sect->crect[0].F/30.0;
  f2=5.0+elem.sect->crect[0].F/100.0;
  if(f1<=f2) cfs=1.5*f1;
  else       cfs=1.5*f2;

  if(elem.sect->srein[1-axis].n>=1)
  {
    if(elem.sect->srein[1-axis].F==3000.0)      wft=3000.0;
    else if(elem.sect->srein[1-axis].F==3500.0) wft=3500.0;

    pw=elem.sect->srein[1-axis].area
      *elem.sect->srein[1-axis].n
      /elem.sect->srein[1-axis].pitch
      /cBi[0];
  }
  else
  {
    if(elem.sect->wF==3000.0)      wft=3000.0;
    else if(elem.sect->wF==3500.0) wft=3500.0;

    pw=elem.sect->shearrein[1-axis];
  }

  if(pw>0.006) pw=0.006;

  rcQa1=7.0/8.0*eAc*(cfs+0.5*pw*wft); /*"SRC STANDARD"(63)*/
  rcQa2=7.0/8.0*eAc*(2.0*eBrate*cfs+pw*wft);

  if(rcQa1<=rcQa2) rcQa=rcQa1; /*"SRC STANDARD"(62)*/
  else             rcQa=rcQa2;

  if(haxis==HSTRONG)
  {
    sQa=allowultimshearofhstrong(PSHORT,/*m.sF*/elem.sect->srect[0].F,
                                 *(elem.sect),axis);
  }
  if(haxis==HWEAK)
  {
    sQa=allowultimshearofhweak(PSHORT,/*m.sF*/elem.sect->srect[0].F,
                               *(elem.sect),axis);
  }
  srcQa=sQa+rcQa;

  /*
  if(fout0!=NULL)
  {
    fprintf(fout0,"Ac=%.3f cYc=%.3f cYt=%.3f[cm2,cm]\n",Ac,cYc,cYt);
    fprintf(fout0,"rYt=%.3f[cm] pw=%.3f\n",rYt,pw);
    fprintf(fout0,"sQa=%.3f[tf] rcQa=%.3f[tf]\n",
            sQa/1000.0,rcQa/1000.0);
    fprintf(fout0,"srcQa=%11.3f[tf]",srcQa/1000.0);
  }
  */

  return srcQa;
}/*allowableshearofsrcshort*/

double ultimateshearofsrc(struct element elem,
                          /*struct materials m,*/
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

  cfs1=0.15*fabs(/*m.Fc*/elem.sect->crect[0].F); /*"SRC STANDARD"(124)[kgf/cm2]*/
  cfs2=22.5+4.5*fabs(/*m.Fc*/elem.sect->crect[0].F)/100;
  if(cfs1<=cfs2) cfsu=cfs1;
  else           cfsu=cfs2;

  if(elem.sect->srein[1-axis].n>=1)
  {
    wfp=elem.sect->srein[1-axis].F;

    pw=elem.sect->srein[1-axis].area
      *elem.sect->srein[1-axis].n
      /elem.sect->srein[1-axis].pitch
      /cBi[0];
  }
  else
  {
    wfp=elem.sect->wF;

    pw=elem.sect->shearrein[1-axis];
  }

  if(pw>0.006) pw=0.006;
  if(pw<0.001)
  {
    if(fout0!=NULL) fprintf(fout0,"HOOP TOO LITTLE.Pw=%.7f < 0.001\n",pw);
  }

  /*rcQu0=2*rcMu/inner;*/                     /*"SRC STANDARD"(121)*/

  rcQu1=7.0/8.0*eAc*(0.5*ac*cfsu+0.5*pw*wfp); /*"SRC STANDARD"(123)*/
  rcQu2=7.0/8.0*eAc*(eBrate*cfsu+pw*wfp);

  if(rcQu1<=rcQu2) rcQu=rcQu1;
  else             rcQu=rcQu2;

  /*sQu0=2*sMu/inner;*/                       /*"SRC STANDARD"(126)*/

  if(haxis==HSTRONG)
  {
    sQu1=allowultimshearofhstrong(PULTIMATE,/*m.sF*/elem.sect->srect[0].F,
                                  *(elem.sect),axis);
  }
  if(haxis==HWEAK)
  {
    sQu1=allowultimshearofhweak(PULTIMATE,/*m.sF*/elem.sect->srect[0].F,
                                *(elem.sect),axis);
  }
  sQu=sQu1;

  srcQu=sQu+rcQu;                             /*"SRC STANDARD"(119)*/

  /*if(fout0!=NULL)
  {
    fprintf(fout0,"Ac=%.3f cYc=%.3f cYt=%.3f[cm2,cm]\n",Ac,cYc,cYt);
    fprintf(fout0,"rYt=%.3f[cm] pw=%.5f\n",rYt,pw);
    fprintf(fout0,"sMu=%.3f rcMu=%.3f[tfm]\n",
            sMu/100000.0,rcMu/100000.0);
    fprintf(fout0,"sQu0=%.3f sQu1=%.3f[tf]\n",
            sQu0/1000.0,sQu1/1000.0);
    fprintf(fout0,"rcQu0=%.3f rcQu1=%.3f rcQu2=%.3f[tf]\n",
            rcQu0/1000.0,rcQu1/1000.0,rcQu2/1000.0);
    fprintf(fout0,"sQu=%.3f[tf] rcQu=%.3f[tf]\n",
            sQu/1000.0,rcQu/1000.0);
    fprintf(fout0,"srcQu=%11.3f[tf]\n",srcQu/1000.0);
  }*/

  return srcQu;
}/*ultimateshearofsrc*/

double allowultimshearofrcwall(struct element elem,
                               /*struct materials m,*/
                               double l, /*LENGTH[cm]*/
                               int period)
/*ALLOWABLE,ULTIMATE SHEAR OF RC WALL.*/
{
  char code[20];
  double t,r;                           /*t:THICKNESS r:WINDOW RATE*/
  double ps;
  double cfs,wft,f1,f2;
  double Qc,Qw,Qau;

  if(period==PULTIMATE)
  {
    cfs=2.4*sqrt(fabs(elem.sect->cF)); /*"RC STANDARD" APP21.(20)*/
    wft=fabs(1.1*elem.sect->wF);
  }
  else if(period==PLONG || period==PSHORT)
  {
    f1=elem.sect->cF/30.0;
    f2=5.0+elem.sect->cF/100.0;

    if(period==PLONG)
    {
      if(f1<=f2) cfs=/*fabs(m.cfs)*/f1;
      else       cfs=/*fabs(m.cfs)*/f2;

      if(elem.sect->wF==3000.0) wft=/*fabs(m.wft)*/2000.0;
      if(elem.sect->wF==3500.0) wft=/*fabs(m.wft)*/2200.0;
    }
    else if(period==PSHORT)
    {
      if(f1<=f2) cfs=/*fabs(m.cfs)*/1.5*f1;
      else       cfs=/*fabs(m.cfs)*/1.5*f2;

      if(elem.sect->wF==3000.0) wft=/*fabs(m.wft)*/3000.0;
      if(elem.sect->wF==3500.0) wft=/*fabs(m.wft)*/3500.0;
    }
  }

  t=elem.sect->thick;
  ps=elem.sect->shearrein[0];
  r=1.0-elem.sect->windowrate;

  /*if(fout0!=NULL)
  {
    fprintf(fout0,"THICK=%.3f Ps=%.3f WINDOW=%.3f",t,ps,r);
    fprintf(fout0," cfs=%.3f wft=%.3f\n",cfs,wft);
  }*/

  Qc=r*t*l*cfs;                             /*"RC STANDARD" 18.(30)*/
  Qw=r*ps*t*l*wft;                          /*"RC STANDARD" 18.(32)*/

  if(period==PLONG) /*RC2010*/
  {
    Qau=Qc;
    strcpy(code,"CONCRETE");
  }
  else
  {
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
  }

  if(fout0!=NULL)
  {
    fprintf(fout0,"CONCRETE:Qc=%7.3f[tf]",Qc/1000.0);
    fprintf(fout0," REINFORCE:Qw=%7.3f[tf]",Qw/1000.0);

    if(period==PULTIMATE) fprintf(fout0," Qu=");
    else                  fprintf(fout0," Qa=");

    fprintf(fout0,"%7.3f[tf](%7.2f[kN]) BY %s.\n",
            Qau/1000.0,SIUNIT*Qau/1000.0,code);
  }

  return Qau;
}/*allowultimshearofrcwall*/

double allowableshearofwoodwall(struct element elem,
                                double l, /*LENGTH[cm]*/
                                int period)
/*ALLOWABLE SHEAR OF WOOD WALL.*/
{
  double t,r;                           /*t:THICKNESS r:WINDOW RATE*/
  double fsa;
  double Qa;

  if(period==PLONG)  fsa=elem.sect->fsl;
  if(period==PSHORT) fsa=elem.sect->fsl*2.0; /*"WOOD STANDARD" P.183*/

  t=elem.sect->thick;
  r=1.0-elem.sect->windowrate;

  Qa=r*t*l*fsa;

  if(fout0!=NULL)
  {
    fprintf(fout0,"Qa=%7.3f[tf](%7.2f[kN])\n",
            Qa/1000.0,SIUNIT*Qa/1000.0);
  }

  return Qa;
}/*allowableshearofwoodwall*/

double allowabletensionofs(double sF,struct section sect)
/*ALLOWABLE TENSION OF STEEL LONG. =SHORT/1.5*/
{
  double A;
  double sfta,sNta;                 /*s:STEEL t:TENSION a:ALLOWABLE*/

  sfta=sF/1.5;                                  /*"S STANDARD"(5.1)*/

  A=coeffA(sect.nsteel,sect.srect);                         /*[cm2]*/
  sNta=-sfta*A; /*-:TENSION*/                               /*[kgf]*/

  /*if(fout0!=NULL)
  {
    fprintf(fout0,"sft=%.3f[kgf/cm2] As=%.3f[cm2]\n",sfta,A);
  }*/

  return sNta;
}/*allowabletensionofs*/

double allowabletensionofsteel(double sF,struct section sect)
/*ALLOWABLE TENSION OF STEEL LONG. =SHORT/1.5*/
{
  double A,Ixx,Iyy,Zxx,Zyy,ixx,iyy;
  double sfta,sNta;                 /*s:STEEL t:TENSION a:ALLOWABLE*/

  steelcoefficients(sect,&A,&Ixx,&Iyy,&Zxx,&Zyy,&ixx,&iyy); /*[cm]*/

  sfta=sF/1.5;                                  /*"S STANDARD"(5.1)*/

  sNta=-sfta*A; /*-:TENSION*/                               /*[kgf]*/

  /*if(fout0!=NULL)
  {
    fprintf(fout0,"sft=%.3f[kgf/cm2] As=%.3f[cm2]\n",sfta,A);
  }*/

  return sNta;
}/*allowabletensionofsteel*/

double allowabletensionofalumi(double sF,struct section sect)
/*ALLOWABLE TENSION OF ALUMI LONG. =SHORT/1.5*/
{
  double A,Ixx,Iyy,Jz,Cw,Zxx,Zyy,ixx,iyy;
  double alfta,alNta;                 /*al:ALUMI t:TENSION a:ALLOWABLE*/

  alfta=sF/1.5;                                  /*ALUMI STANDARD (4.1)*/

  alumicoefficients(sect,&A,&Ixx,&Iyy,&Jz,&Cw,&Zxx,&Zyy,&ixx,&iyy);     /*[cm]*/

  alNta=-alfta*A; /*-:TENSION*/                                /*[kgf]*/


  return alNta;
}/*allowabletensionofalumi*/

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

if(!strcmp(prj,"aki") &&
   sect.code==111 &&
   axis==SX) lk=660.0;

if(!strcmp(prj,"kyok") &&
   (sect.code==202 ||
    sect.code==205 ||
    sect.code==206 ||
    sect.code==207) &&
   axis==SY) lk=90.0;
if(!strcmp(prj,"kyok") &&
   sect.code==206 &&
   axis==SX) lk=500.0;
if(!strcmp(prj,"kyok") &&
   sect.code==207 &&
   axis==SX) lk=850.0;

  lam=lk/i;                                    /*"S STANDARD"(11.1)*/
  if((sect.etype==COLUMN && lam>200)||(lam>250))
  {
    if(fout0!=NULL)
    {
      fprintf(fout0,"ELEMENT TOO LONG.Lk/i=%.5f\n",lam);
    }
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
  if(fout0!=NULL)
  {
    fprintf(fout0,"Lk=%.3f[cm] i=%.3f[cm] lam=%.3f LAM=%.3f NYU=%.3f\n",
            lk,i,lam,LAM,nyu);
    fprintf(fout0,"sfc=%.3f[kgf/cm2] As=%.3f[cm2]\n",sfca,A);
  }
  */

  return sNca;
}/*allowablecompressionofs*/

double allowablecompressionofsteel(double E,double F,
                                   double lkxx,double lkyy, /*BUCKLING LENGTH[cm]*/
                                   struct section sect)
/*ALLOWABLE COMPRESSION OF STEEL LONG. =SHORT/1.5*/
{
  double A,Ixx,Iyy,Zxx,Zyy,ixx,iyy,i;
  double lk,lam,lamxx,lamyy,LAM,nyu;
  double sfca,sNca;             /*s:STEEL c:COMPRESSION a:ALLOWABLE*/

  LAM=sqrt(PI*PI*E/0.6/F);                      /*"S STANDARD"(5.5)*/

  steelcoefficients(sect,&A,&Ixx,&Iyy,&Zxx,&Zyy,&ixx,&iyy);  /*[cm]*/

  /*if(ixx<=iyy) i=ixx;*/
  /*else         i=iyy;*/
  /*lam=lk/i;*/                                /*"S STANDARD"(11.1)*/

  if(sect.bblength[0]!=0.0) lk=sect.bblength[0];
  else if(sect.bbfact[0]!=0.0) lk=lkxx*sect.bbfact[0];
  else lk=lkxx;

  lamxx=lk/ixx;

  if(sect.bblength[1]!=0.0) lk=sect.bblength[1];
  else if(sect.bbfact[1]!=0.0) lk=lkyy*sect.bbfact[1];
  else lk=lkyy;

  lamyy=lk/iyy;

  if(lamxx>lamyy) lam=lamxx;
  else            lam=lamyy;

  if((sect.etype==COLUMN && lam>200)||(lam>250))
  {
    if(fout0!=NULL)
    {
      fprintf(fout0,"ELEMENT TOO LONG.Lk/i=%.5f\n",lam);
    }
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

  sNca=sfca*A; /*+:COMPRESSION*/                            /*[kgf]*/

  /*
  if(fout0!=NULL)
  {
    fprintf(fout0,"Lk=%.3f[cm] i=%.3f[cm] lam=%.3f LAM=%.3f NYU=%.3f\n",
            lk,i,lam,LAM,nyu);
    fprintf(fout0,"sfc=%.3f[kgf/cm2] As=%.3f[cm2]\n",sfca,A);
  }
  */

  return sNca;
}/*allowablecompressionofsteel*/

double allowablecompressionofalumi(double E,double F,
                                   double lkxx,double lkyy, /*BUCKLING LENGTH[cm]*/
                                   struct section sect,int axis)
/*ALLOWABLE COMPRESSION OF ALUMI LONG. =SHORT/1.5*/
{
  double A,Ixx,Iyy,Jz,Cw,Zxx,Zyy,ixx,iyy,i;
  double lk,ramda,pramda,eramda,nyu,lam;
  double Ny,Ne;
  double alfca,alNca;             /*al:ALUMI c:COMPRESSION a:ALLOWABLE*/

  alumicoefficients(sect,&A,&Ixx,&Iyy,&Jz,&Cw,&Zxx,&Zyy,&ixx,&iyy);     /*[cm]*/

  if(axis==SX)
  {
    if(sect.bblength[0]!=0.0) lk=sect.bblength[0];
    else if(sect.bbfact[0]!=0.0) lk=lkxx*sect.bbfact[0];
    else lk=lkxx;
  }

  else if(axis==SY)
  {
    if(sect.bblength[1]!=0.0) lk=sect.bblength[1];
    else if(sect.bbfact[1]!=0.0) lk=lkyy*sect.bbfact[1];
    else lk=lkyy;
  }

  Ny=F*A;  /*kgf*/
  if(Ixx>Iyy) Ne=PI*PI*E*Iyy/lk/lk;
  else        Ne=PI*PI*E*Ixx/lk/lk;

  ramda=sqrt(Ny/Ne);                          /*ALUMI STANDARD (4.2.1)*/
  pramda=0.2;
  eramda=1/sqrt(0.5);

  if(ixx<=iyy) i=ixx;
  else         i=iyy;

  lam=lk/i;                                    /*ALUMI STANDARD(4.2.2)*/
  if((sect.etype==COLUMN && lam>140) || (sect.etype!=COLUMN && lam>180))
  {
    if(fout0!=NULL)
    {
	  fprintf(fout0,"ELEMENT TOO LONG.Lk/i=%.5f\n",lam);    /*Keisansho*/
    }
  }

  nyu=1.5+(ramda/eramda)*(ramda/eramda)/1.5;
  if(nyu>2.17) nyu=2.17;

  if(ramda<=pramda)
  {
    alfca=F/nyu;

  }
  else if(pramda<ramda && ramda<=eramda)
  {
    alfca=(1.0-0.5*(ramda-pramda)/(eramda-pramda))*F/nyu;
  }
  else
  {
    alfca=F/ramda/ramda/nyu;
  }

  alNca=alfca*A; /*+:COMPRESSION*/                            /*[kgf]*/

  return alNca;
}/*allowablecompressionofalumi*/


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
    if(fout0!=NULL) fprintf(fout0,"STEEL NOT H,[.\n");
    return 0.0;
  }

  if(Nd>0.0) sNa=allowablecompressionofs(m.sE,sect.srect[0].F,lb,sect,
                                         (1-axis)); /*FOR WEAK.*/
  if(Nd<0.0) sNa=allowabletensionofs(sect.srect[0].F,sect);

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

  /*sfta=m.sF;*/ /*[kgf/cm2]*/
  sfta=sect.srect[0].F; /*[kgf/cm2]*/
  if(period==PLONG) sfta/=1.5;

  sfba=900.0*Af/lb/H*1000.0;          /*[kgf/cm2] "S STANDARD"(5.8)*/
  if(period==PSHORT) sfba*=1.5;
  if(sfba>sfta) sfba=sfta;

  Z=coeffZ(sect.nsteel,sect.srect,axis);                    /*[cm3]*/

  sMa=(1.0-Nd/sNa)*sfba*Z;                                /*[kgfcm]*/

  /*
  if(fout0!=NULL)
  {
    fprintf(fout0,"STRONG. Nd=%.3f[kgf] sNa=%.3f[kgf] Nd/sNa=%.3f\n",
            Nd,sNa,(Nd/sNa));
    fprintf(fout0,"H=%.1f[cm] B=%.1f[cm] tf=%.1f[cm]\n",H,B,tf);
    fprintf(fout0,"sft=%.3f[kgf/cm2] sfb=%.3f[kgf/cm2] lb=%.3f[cm]",
            sfta,sfba,lb);
    fprintf(fout0," Z=%.3f[cm3]\n",Z);
  }
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

  sfb=sect.srect[0].F; /*"S STANDARD" 5.1.(4).(b)*/
  if(period==PLONG) sfb/=1.5;

  if(Nd<0.0) sNa=allowabletensionofs(sect.srect[0].F,sect);
  if(Nd>0.0) sNa=allowablecompressionofs(m.sE,sect.srect[0].F,lk,sect,axis);

  if(period==PSHORT) sNa*=1.5;

  if((Nd/sNa)>=1.0) return 0.0;

  Ze=coeffZ(sect.nsteel,sect.srect,axis);

  sMa=Ze*sfb*(1.0-Nd/sNa); /*"S STANDARD"(6.1)(6.2)(6.3)(6.4)*/

  if(sMa<0.0) sMa=0.0;

  /*
  if(fout0!=NULL)
  {
    fprintf(fout0,"WEAK. Nd=%.3f[kgf] sNa=%.3f[kgf] Nd/sNa=%.3f\n",
            Nd,sNa,(Nd/sNa));
    fprintf(fout0,"sfb=%.3f[kgf/cm2]",sfb);
    fprintf(fout0," Z=%.3f[cm3]\n",Ze);
    fprintf(fout0,"ìSçúé≤óÕ sN=%.3f[tf] ãñóeíl sMa=%.3f[tfm]",
            Nd/1000.0,sMa/100000.0);
  }
  */

  return sMa;
}/*allowablebendingofsweak*/

double allowultimshearofhstrong(int period,double F,
                                struct section sect,int axis)
/*ALLOWABLE,ULTIMATE SHEAR OF H,[ STRONG.*/
{
  double H,tw,tf;
  double sfsa,sQa;                    /*s:STEEL s:SHEAR a:ALLOWABLE*/

  if(sect.nsteel!=3)
  {
    if(fout0!=NULL) fprintf(fout0,"STEEL NOT H,[.\n");
    return 0.0;
  }

  if(axis==SX) /*SECTION AXIS AROUND X,WITH STRONG H.*/
  {
    H =(sect.srect[0].top-sect.srect[2].bottom);             /*[cm]*/
    tw=(sect.srect[1].right-sect.srect[1].left);             /*[cm]*/
    tf=(sect.srect[0].top-sect.srect[0].bottom);             /*[cm]*/
  }
  if(axis==SY) /*SECTION AXIS AROUND Y,WITH STRONG H.*/
  {
    H =(sect.srect[2].right-sect.srect[0].left);             /*[cm]*/
    tw=(sect.srect[1].top-sect.srect[1].bottom);             /*[cm]*/
    tf=(sect.srect[0].right-sect.srect[0].left);             /*[cm]*/
  }

  sfsa=F/sqrt(3.0)/1.5;    /*"S STANDARD"(5.1),"SRC STANDARD"(127).*/

  sQa=sfsa*tw*(H-2.0*tf);       /*[kgf] SAME AS "SRC STANDARD"(42).*/
  if(period==PSHORT || period==PULTIMATE) sQa*=1.5;

  /*
  if(fout0!=NULL)
  {
    fprintf(fout0,"STRONG. H=%.1f[cm] B=%.1f[cm]",H,B);
    fprintf(fout0," tw=%.1f[cm] tf=%.1f[cm]\n",tw,tf);
    fprintf(fout0,"sfs=%.3f[kgf/cm2]\n",sfsa);
  }
  */

  return sQa;
}/*allowultimshearofhstrong*/

double allowultimshearofhweak(int period,double F,
                              struct section sect,int axis)
/*ALLOWABLE,ULTIMATE SHEAR OF H,[ WEAK.*/
{
  double B,tf;
  double sfsa,sQa;                    /*s:STEEL s:SHEAR a:ALLOWABLE*/

  if(sect.nsteel!=3)
  {
    if(fout0!=NULL) fprintf(fout0,"STEEL NOT H,[.\n");
    return 0.0;
  }

  if(axis==SY) /*SECTION AXIS AROUND Y,WITH WEAK H.*/
  {
    B =(sect.srect[0].right-sect.srect[0].left);             /*[cm]*/
    tf=(sect.srect[0].top-sect.srect[0].bottom);             /*[cm]*/
  }
  if(axis==SX) /*SECTION AXIS AROUND X,WITH WEAK H.*/
  {
    B =(sect.srect[0].top-sect.srect[0].bottom);             /*[cm]*/
    tf=(sect.srect[0].right-sect.srect[0].left);             /*[cm]*/
  }

  sfsa=F/sqrt(3.0)/1.5;    /*"S STANDARD"(5.1),"SRC STANDARD"(127).*/

  sQa=sfsa*2.0/3.0*(2.0*tf*B);  /*[kgf] SAME AS "SRC STANDARD"(61).*/
  if(period==PSHORT || period==PULTIMATE) sQa*=1.5;

  /*
  if(fout0!=NULL)
  {
    fprintf(fout0,"WEAK. H=%.1f[cm] B=%.1f[cm]",H,B);
    fprintf(fout0," tw=%.1f[cm] tf=%.1f[cm]\n",tw,tf);
    fprintf(fout0,"sfs=%.3f[kgf/cm2]\n",sfsa);
  }
  */

  return sQa;
}/*allowultimshearofhweak*/

double allowultimshearofsteel(int period,double F,
                              struct section sect,int axis)
/*ALLOWABLE,ULTIMATE SHEAR OF STEEL.*/
/*AXIS SX=Qy,SY=Qx.*/
{
  double As=0.0,H,B,tw,tf;
  double sfsa,sQa;                    /*s:STEEL s:SHEAR a:ALLOWABLE*/

  if(sect.sform.type==STEEL_RECTS) F=sect.srect[0].F; /*[kgf/cm2]*/
  else                             F=sect.sform.F;

  if(sect.sform.type==STEEL_RECTS && sect.nsteel!=3)
  {
    if(fout0!=NULL) fprintf(fout0,"STEEL NOT H,[.\n");
    return 0.0;
  }

  if(sect.sform.type==STEEL_RECTS) /*FOR ONLY HKYOU INPUT BY RECTS.*/
  {
    H =(sect.srect[0].top-sect.srect[2].bottom);             /*[cm]*/
    B =(sect.srect[0].right-sect.srect[0].left);             /*[cm]*/
    tw=(sect.srect[1].right-sect.srect[1].left);             /*[cm]*/
    tf=(sect.srect[0].top-sect.srect[0].bottom);             /*[cm]*/
  }
  else
  {
    H =sect.sform.H;  /*[cm]*/
    B =sect.sform.B;  /*[cm]*/
    tw=sect.sform.tw; /*[cm]*/
    tf=sect.sform.tf; /*[cm]*/
  }

  if(axis==SX)
  {
    if(sect.sform.type==STEEL_RECTS) As=tw*(H-2.0*tf);
    if(sect.sform.type==STEEL_HKYOU) As=tw*(H-2.0*tf);
    if(sect.sform.type==STEEL_HWEAK) As=2.0*tf*B/1.5;
    if(sect.sform.type==STEEL_RPIPE) As=2.0*tw*(H-2.0*tf);
    if(sect.sform.type==STEEL_CPIPE)
    {
      As=PI*(H*H-(H-2.0*tf)*(H-2.0*tf))/4.0/2.0;
    }
    if(sect.sform.type==STEEL_ANGLE) As=tw*H/1.5;
    if(sect.sform.type==STEEL_TKYOU) As=tw*(H-tf)/1.5;
    if(sect.sform.type==STEEL_TWEAK) As=tf*B/1.5;
  }
  else
  {
    if(sect.sform.type==STEEL_RECTS) As=2.0*tf*B/1.5;
    if(sect.sform.type==STEEL_HKYOU) As=2.0*tf*B/1.5;
    if(sect.sform.type==STEEL_HWEAK) As=tw*(H-2.0*tf);
    if(sect.sform.type==STEEL_RPIPE) As=2.0*tf*(B-2.0*tw);
    if(sect.sform.type==STEEL_CPIPE)
    {
      As=PI*(H*H-(H-2.0*tf)*(H-2.0*tf))/4.0/2.0;
    }
    if(sect.sform.type==STEEL_ANGLE) As=tf*B/1.5;
    if(sect.sform.type==STEEL_TKYOU) As=tf*B/1.5;
    if(sect.sform.type==STEEL_TWEAK) As=tw*(H-tf)/1.5;
  }

  sfsa=F/sqrt(3.0)/1.5;    /*"S STANDARD"(5.1),"SRC STANDARD"(127).*/

  sQa=sfsa*As;                  /*[kgf] SAME AS "SRC STANDARD"(42).*/
  if(period==PSHORT || period==PULTIMATE) sQa*=1.5;

  /*
  if(fout0!=NULL)
  {
    fprintf(fout0,"H=%.1f B=%.1f tw=%.1f tf=%.1f\n",H,B,tw,tf);
    fprintf(fout0,"sfs=%.3f[kgf/cm2]\n",sfsa);
  }
  */

  return sQa;
}/*allowultimshearofsteel*/

double allowultimshearofalumi(int period,double F,
                              struct section sect,int axis)
/*ALLOWABLE,ULTIMATE SHEAR OF ALUMI.*/
/*AXIS SX=Qy,SY=Qx.*/
{
  double As=0.0,H,B,tw,tf;
  double sfsa,sQa;                    /*s:STEEL s:SHEAR a:ALLOWABLE*/

  if(sect.sform.type==STEEL_RECTS && sect.nsteel!=3)
  {
    if(fout0!=NULL) fprintf(fout0,"STEEL NOT H,[.\n");
    return 0.0;
  }

  if(sect.sform.type==STEEL_RECTS) /*FOR ONLY HKYOU INPUT BY RECTS.*/
  {
    H =(sect.srect[0].top-sect.srect[2].bottom);             /*[cm]*/
    B =(sect.srect[0].right-sect.srect[0].left);             /*[cm]*/
    tw=(sect.srect[1].right-sect.srect[1].left);             /*[cm]*/
    tf=(sect.srect[0].top-sect.srect[0].bottom);             /*[cm]*/
  }
  else
  {
    H =sect.sform.H;  /*[cm]*/
    B =sect.sform.B;  /*[cm]*/
    tw=sect.sform.tw; /*[cm]*/
    tf=sect.sform.tf; /*[cm]*/
  }

  if(axis==SX)
  {
    if(sect.sform.type==STEEL_RECTS) As=tw*(H-2.0*tf);
    if(sect.sform.type==STEEL_HKYOU) As=tw*(H-2.0*tf);
    if(sect.sform.type==STEEL_HWEAK) As=2.0*tf*B/1.5;
    if(sect.sform.type==STEEL_RPIPE) As=2.0*tw*(H-2.0*tf);
    if(sect.sform.type==STEEL_CPIPE)
    {
      As=PI*(H*H-(H-2.0*tf)*(H-2.0*tf))/4.0/2.0;
    }
    if(sect.sform.type==STEEL_ANGLE) As=tw*H/1.5;
    if(sect.sform.type==STEEL_TKYOU) As=tw*(H-tf)/1.5;
    if(sect.sform.type==STEEL_TWEAK) As=tf*B/1.5;
    if(sect.sform.type==STEEL_CKYOU) As=tw*(H-2*tf)/1.5;
    if(sect.sform.type==STEEL_CWEAK) As=tf*B/1.5;
  }
  else
  {
    if(sect.sform.type==STEEL_RECTS) As=2.0*tf*B/1.5;
    if(sect.sform.type==STEEL_HKYOU) As=2.0*tf*B/1.5;
    if(sect.sform.type==STEEL_HWEAK) As=tw*(H-2.0*tf);
    if(sect.sform.type==STEEL_RPIPE) As=2.0*tf*(B-2.0*tw);
    if(sect.sform.type==STEEL_CPIPE)
    {
      As=PI*(H*H-(H-2.0*tf)*(H-2.0*tf))/4.0/2.0;
    }
    if(sect.sform.type==STEEL_ANGLE) As=tf*B/1.5;
    if(sect.sform.type==STEEL_TKYOU) As=tf*B/1.5;
    if(sect.sform.type==STEEL_TWEAK) As=tw*(H-tf)/1.5;
    if(sect.sform.type==STEEL_CKYOU) As=tf*B/1.5;
    if(sect.sform.type==STEEL_CWEAK) As=tw*(H-2*tf)/1.5;
  }

  sfsa=F/sqrt(3.0)/1.5;    /*"S STANDARD"(5.1),"SRC STANDARD"(127).*/

  sQa=sfsa*As;                  /*[kgf] SAME AS "SRC STANDARD"(42).*/
  if(period==PSHORT || period==PULTIMATE) sQa*=1.5;

  /*
  if(fout0!=NULL)
  {
    fprintf(fout0,"H=%.1f B=%.1f tw=%.1f tf=%.1f\n",H,B,tw,tf);
    fprintf(fout0,"sfs=%.3f[kgf/cm2]\n",sfsa);
  }
  */

  return sQa;
}/*allowultimshearofalumi*/

double allowableshearofwood(int period,struct section sect)
/*ALLOWABLE SHEAR OF WOOD.*/
/*AXIS SX=Qy,SY=Qx.*/
{
  double As,H,B;
  double fsa,Qa; /*s:SHEAR a:ALLOWABLE*/

  H =sect.sform.H; /*[cm]*/
  B =sect.sform.B; /*[cm]*/
  As=B*H/1.5;

  if(period==PLONG)  fsa=1.1/3.0*sect.wform.Fs;
  if(period==PSHORT) fsa=2.0/3.0*sect.wform.Fs;

  Qa=fsa*As; /*[kgf]*/

  return Qa;
}/*allowableshearofwood*/

int steelcoefficients(struct section sect,
                      double *A,
                      double *Ixx,double *Iyy,
                      double *Zxx,double *Zyy,
                      double *ixx,double *iyy)
/*COEFFICIENTS OF STEEL.*/
{
  char str[400],s[80];
  double H,B,tw,tf;   /*[cm]*/
  double xg,yg,yt,yc; /*[cm]*/
  double Ixy,Ixx2,Iyy2,theta,Zupper,Zlower; /*[cm3]*/

  /*COPY FORM*/
  if(sect.sform.type!=STEEL_RECTS)
  {
    H =sect.sform.H;
    B =sect.sform.B;
    tw=sect.sform.tw;
    tf=sect.sform.tf;
  }

  /*AREA*/
  if(sect.sform.type==STEEL_RECTS)
  {
    *A=coeffA(sect.nsteel,sect.srect);
  }
  else if(sect.sform.type==STEEL_HKYOU ||
          sect.sform.type==STEEL_HWEAK)
  {
    *A=2.0*tf*B+tw*(H-2.0*tf);
  }
  else if(sect.sform.type==STEEL_RPIPE)
  {
    *A=2.0*tf*B+2.0*tw*(H-2.0*tf);
  }
  else if(sect.sform.type==STEEL_CPIPE)
  {
    *A=PI*(H*H-(H-2.0*tf)*(H-2.0*tf))/4.0;
  }
  else if(sect.sform.type==STEEL_ANGLE)
  {
    *A=H*tw+B*tf-tw*tf;
  }
  else if(sect.sform.type==STEEL_TKYOU ||
          sect.sform.type==STEEL_TWEAK)
  {
    *A=tf*B+tw*(H-tf);
  }

  /*Ixx*/
  if(sect.sform.type==STEEL_RECTS)
  {
    *Ixx=coeffI(sect.nsteel,sect.srect,SX,&yg,&yt,&yc);
  }
  else if(sect.sform.type==STEEL_HKYOU)
  {
    *Ixx=(B*H*H*H-(B-tw)*(H-2.0*tf)*(H-2.0*tf)*(H-2.0*tf))/12.0;
    yg=0.0;
    yc= H/2.0;
    yt=-H/2.0;
  }
  else if(sect.sform.type==STEEL_HWEAK)
  {
    *Ixx=(2.0*tf*B*B*B+(H-2.0*tf)*tw*tw*tw)/12.0;
    yg=0.0;
    yc= B/2.0;
    yt=-B/2.0;
  }
  else if(sect.sform.type==STEEL_RPIPE)
  {
    *Ixx=(B*H*H*H-(B-2.0*tw)*(H-2.0*tf)*(H-2.0*tf)*(H-2.0*tf))/12.0;
    yg=0.0;
    yc= H/2.0;
    yt=-H/2.0;
  }
  else if(sect.sform.type==STEEL_CPIPE)
  {
    *Ixx=PI*(H*H*H*H-(H-2.0*tf)*(H-2.0*tf)*(H-2.0*tf)*(H-2.0*tf))/64.0;
    yg=0.0;
    yc= H/2.0;
    yt=-H/2.0;
  }
  else if(sect.sform.type==STEEL_ANGLE)
  {
    yg=(tw*H*H/2.0+(B-tw)*tf*tf/2.0)/(*A); /*FROM BOTTOM*/
    yc=H;
    yt=0.0;

    *Ixx=1.0/12.0*tw*H*H*H+tw*H*(H/2-yg)*(H/2-yg)
        +1.0/12.0*(B-tw)*tf*tf*tf
        +(B-tw)*tf*(yg-tf/2.0)*(yg-tf/2.0);
  }
  else if(sect.sform.type==STEEL_TKYOU)
  {
    yg=(tw*H*H/2.0+(B-tw)*tf*(H-tf/2.0))/(*A); /*FROM BOTTOM*/
    yc=H;
    yt=0.0;

    *Ixx=1.0/12.0*tw*H*H*H+tw*H*(H/2-yg)*(H/2-yg)
        +1.0/12.0*(B-tw)*tf*tf*tf
        +(B-tw)*tf*(H-tf/2.0-yg)*(H-tf/2.0-yg);
  }
  else if(sect.sform.type==STEEL_TWEAK)
  {
    *Ixx=(tf*B*B*B+(H-tf)*tw*tw*tw)/12.0;
    yg=0.0;
    yc= B/2.0;
    yt=-B/2.0;
  }

  /*Zxx*/
  Zupper=*Ixx/fabs(yc-yg);
  Zlower=*Ixx/fabs(yg-yt);
  if(Zupper<=Zlower) *Zxx=Zupper;
  else               *Zxx=Zlower;

  /*Iyy*/
  if(sect.sform.type==STEEL_RECTS)
  {
    *Iyy=coeffI(sect.nsteel,sect.srect,SY,&xg,&yt,&yc);
  }
  else if(sect.sform.type==STEEL_HKYOU)
  {
    *Iyy=(2.0*tf*B*B*B+(H-2.0*tf)*tw*tw*tw)/12.0;
    xg=0.0;
    yc= B/2.0;
    yt=-B/2.0;
  }
  else if(sect.sform.type==STEEL_HWEAK)
  {
    *Iyy=(B*H*H*H-(B-tw)*(H-2.0*tf)*(H-2.0*tf)*(H-2.0*tf))/12.0;
    xg=0.0;
    yc= H/2.0;
    yt=-H/2.0;
  }
  else if(sect.sform.type==STEEL_RPIPE)
  {
    *Iyy=(H*B*B*B-(H-2.0*tf)*(B-2.0*tw)*(B-2.0*tw)*(B-2.0*tw))/12.0;
    xg=0.0;
    yc= B/2.0;
    yt=-B/2.0;
  }
  else if(sect.sform.type==STEEL_CPIPE)
  {
    *Iyy=*Ixx;
    xg=0.0;
    yc= H/2.0;
    yt=-H/2.0;
  }
  else if(sect.sform.type==STEEL_ANGLE)
  {
    xg=(tf*B*B/2.0+(H-tf)*tw*tw/2.0)/(*A); /*FROM LEFT*/
    yc=B;
    yt=0.0;

    *Iyy=1.0/12.0*tf*B*B*B+tf*B*(B/2-xg)*(B/2-xg)
        +1.0/12.0*(H-tf)*tw*tw*tw
        +(H-tf)*tw*(xg-tw/2.0)*(xg-tw/2.0);
  }
  else if(sect.sform.type==STEEL_TKYOU)
  {
    *Iyy=(tf*B*B*B+(H-tf)*tw*tw*tw)/12.0;
    xg=0.0;
    yc= B/2.0;
    yt=-B/2.0;
  }
  else if(sect.sform.type==STEEL_TWEAK)
  {
    xg=(tw*H*H/2.0+(B-tw)*tf*tf/2.0)/(*A); /*FROM LEFT*/
    yc=H;
    yt=0.0;

    *Iyy=1.0/12.0*tw*H*H*H+tw*H*(H/2-xg)*(H/2-xg)
        +1.0/12.0*(B-tw)*tf*tf*tf
        +(B-tw)*tf*(xg-tf/2.0)*(xg-tf/2.0);
  }

  /*Zyy*/
  Zupper=*Iyy/fabs(yc-xg);
  Zlower=*Iyy/fabs(xg-yt);
  if(Zupper<=Zlower) *Zyy=Zupper;
  else               *Zyy=Zlower;

  /*ixx,iyy*/
  *ixx=sqrt((*Ixx)/(*A));
  *iyy=sqrt((*Iyy)/(*A));

  if(sect.sform.type==STEEL_ANGLE)
  {
    Ixy=fabs(0.25*((xg*xg-(B-xg)*(B-xg))*((-yg+tf)*(-yg+tf)-yg*yg)
                  +(xg*xg-(xg-tw)*(xg-tw))*((H-yg)*(H-yg)-(-yg+tf)*(-yg+tf))));

    if((*Ixx)==(*Iyy)) theta=0.25*PI;
    else               theta=0.5*(atan(2.0*Ixy/((*Iyy)-(*Ixx))));

    Ixx2=(*Ixx)*sin(theta)*sin(theta)+(*Iyy)*cos(theta)*cos(theta)+Ixy*sin(2*theta);
    Iyy2=(*Iyy)*sin(theta)*sin(theta)+(*Ixx)*cos(theta)*cos(theta)-Ixy*sin(2*theta);

    *ixx=sqrt(Ixx2/(*A));
    *iyy=sqrt(Iyy2/(*A));
  }

/*
sprintf(str,"A  =%8.3f\n",*A);
sprintf(s,"Gx =%8.3f Gy =%8.3f\n",xg,yg);       strcat(str,s);
sprintf(s,"Ixx=%8.3f Iyy=%8.3f\n",*Ixx,*Iyy);   strcat(str,s);
sprintf(s,"Zxx=%8.3f Zyy=%8.3f\n",*Zxx,*Zyy);   strcat(str,s);
sprintf(s,"Ixy=%8.3f\n",Ixy);                   strcat(str,s);
sprintf(s,"Ixx'=%8.3f Iyy'=%8.3f\n",Ixx2,Iyy2); strcat(str,s);
sprintf(s,"ixx=%8.3f iyy=%8.3f\n",*ixx,*iyy);   strcat(str,s);
MessageBox(NULL,str,"Coefficients",MB_OK);
*/

  return 1;
}/*steelcoefficients*/

int alumicoefficients(struct section sect,
                      double *A,
                      double *Ixx,double *Iyy,
                      double *Jz,double *Cw,
                      double *Zxx,double *Zyy,
                      double *ixx,double *iyy)
/*COEFFICIENTS OF ALUMI.*/
{
  char str[400],s[80];
  double H,B,tw,tf;   /*[cm]*/
  double xg,yg,yt,yc; /*[cm]*/
  double S,Ixy,Ixx2,Iyy2,theta,Zupper,Zlower; /*[cm3]*/
  double Aw;          /*[cm2]*/

  /*COPY FORM*/
  if(sect.sform.type!=STEEL_RECTS)
  {
    H =sect.sform.H;
    B =sect.sform.B;
    tw=sect.sform.tw;
    tf=sect.sform.tf;
  }

  /*AREA*/
  if(sect.sform.type==STEEL_RECTS)
  {
    *A=coeffA(sect.nsteel,sect.srect);
  }
  else if(sect.sform.type==STEEL_HKYOU ||
          sect.sform.type==STEEL_HWEAK)
  {
    *A=2.0*tf*B+tw*(H-2.0*tf);
  }
  else if(sect.sform.type==STEEL_RPIPE)
  {
    *A=2.0*tf*B+2.0*tw*(H-2.0*tf);
  }
  else if(sect.sform.type==STEEL_CPIPE)
  {
    *A=PI*(H*H-(H-2.0*tf)*(H-2.0*tf))/4.0;
  }
  else if(sect.sform.type==STEEL_ANGLE)
  {
    *A=H*tw+B*tf-tw*tf;
  }
  else if(sect.sform.type==STEEL_TKYOU ||
          sect.sform.type==STEEL_TWEAK)
  {
    *A=tf*B+tw*(H-tf);
  }
  else if(sect.sform.type==STEEL_CKYOU ||
          sect.sform.type==STEEL_CWEAK)
  {
    *A=2*tf*B+tw*(H-2*tf);
  }

  /*Ixx*/
  if(sect.sform.type==STEEL_RECTS)
  {
    *Ixx=coeffI(sect.nsteel,sect.srect,SX,&yg,&yt,&yc);
  }
  else if(sect.sform.type==STEEL_HKYOU)
  {
    *Ixx=(B*H*H*H-(B-tw)*(H-2.0*tf)*(H-2.0*tf)*(H-2.0*tf))/12.0;
    yg=0.0;
    yc= H/2.0;
    yt=-H/2.0;
  }
  else if(sect.sform.type==STEEL_HWEAK)
  {
    *Ixx=(2.0*tf*B*B*B+(H-2.0*tf)*tw*tw*tw)/12.0;
    yg=0.0;
    yc= B/2.0;
    yt=-B/2.0;
  }
  else if(sect.sform.type==STEEL_RPIPE)
  {
    *Ixx=(B*H*H*H-(B-2.0*tw)*(H-2.0*tf)*(H-2.0*tf)*(H-2.0*tf))/12.0;
    yg=0.0;
    yc= H/2.0;
    yt=-H/2.0;
  }
  else if(sect.sform.type==STEEL_CPIPE)
  {
    *Ixx=PI*(H*H*H*H-(H-2.0*tf)*(H-2.0*tf)*(H-2.0*tf)*(H-2.0*tf))/64.0;
    yg=0.0;
    yc= H/2.0;
    yt=-H/2.0;
  }
  else if(sect.sform.type==STEEL_ANGLE)
  {
    yg=(tw*H*H/2.0+(B-tw)*tf*tf/2.0)/(*A); /*FROM BOTTOM*/
    yc=H;
    yt=0.0;

    *Ixx=1.0/12.0*tw*H*H*H+tw*H*(H/2-yg)*(H/2-yg)
        +1.0/12.0*(B-tw)*tf*tf*tf
        +(B-tw)*tf*(yg-tf/2.0)*(yg-tf/2.0);
  }
  else if(sect.sform.type==STEEL_TKYOU)
  {
    yg=(tw*H*H/2.0+(B-tw)*tf*(H-tf/2.0))/(*A); /*FROM BOTTOM*/
    yc=H;
    yt=0.0;

    *Ixx=1.0/12.0*tw*H*H*H+tw*H*(H/2-yg)*(H/2-yg)
        +1.0/12.0*(B-tw)*tf*tf*tf
        +(B-tw)*tf*(H-tf/2.0-yg)*(H-tf/2.0-yg);
  }
  else if(sect.sform.type==STEEL_TWEAK)
  {
    *Ixx=(tf*B*B*B+(H-tf)*tw*tw*tw)/12.0;
    yg=0.0;
    yc= B/2.0;
    yt=-B/2.0;
  }
  else if(sect.sform.type==STEEL_CKYOU)
  {
    *Ixx=(B*H*H*H-(B-tw)*(H-2.0*tf)*(H-2.0*tf)*(H-2.0*tf))/12.0;
    yg=0.0;
    yc= H/2.0;
    yt=-H/2.0;
  }
  else if(sect.sform.type==STEEL_CWEAK)
  {
    S=2.0*((B*tf)*B/2.0)+(H-2.0*tf)*tw*(tw/2.0);
    yg=S/(*A);                 /*FROM BOTTOM*/
    yc=B;
    yt=0.0;

    *Ixx=2.0*(1.0/12.0*tw*B*B*B+tw*B*(B/2.0-yg)*(B/2.0-yg))
        +1.0/12.0*(H-2.0*tf)*tw*tw*tw+(H-2.0*tf)*tw*(yg-tw/2.0)*(yg-tw/2.0);
  }

  /*Zxx*/
  Zupper=*Ixx/fabs(yc-yg);
  Zlower=*Ixx/fabs(yg-yt);
  if(Zupper<=Zlower) *Zxx=Zupper;
  else               *Zxx=Zlower;

  /*Iyy*/
  if(sect.sform.type==STEEL_RECTS)
  {
    *Iyy=coeffI(sect.nsteel,sect.srect,SY,&xg,&yt,&yc);
  }
  else if(sect.sform.type==STEEL_HKYOU)
  {
    *Iyy=(2.0*tf*B*B*B+(H-2.0*tf)*tw*tw*tw)/12.0;
    xg=0.0;
    yc= B/2.0;
    yt=-B/2.0;
  }
  else if(sect.sform.type==STEEL_HWEAK)
  {
    *Iyy=(B*H*H*H-(B-tw)*(H-2.0*tf)*(H-2.0*tf)*(H-2.0*tf))/12.0;
    xg=0.0;
    yc= H/2.0;
    yt=-H/2.0;
  }
  else if(sect.sform.type==STEEL_RPIPE)
  {
    *Iyy=(H*B*B*B-(H-2.0*tf)*(B-2.0*tw)*(B-2.0*tw)*(B-2.0*tw))/12.0;
    xg=0.0;
    yc= B/2.0;
    yt=-B/2.0;
  }
  else if(sect.sform.type==STEEL_CPIPE)
  {
    *Iyy=*Ixx;
    xg=0.0;
    yc= H/2.0;
    yt=-H/2.0;
  }
  else if(sect.sform.type==STEEL_ANGLE)
  {
    xg=(tf*B*B/2.0+(H-tf)*tw*tw/2.0)/(*A); /*FROM LEFT*/
    yc=B;
    yt=0.0;

    *Iyy=1.0/12.0*tf*B*B*B+tf*B*(B/2-xg)*(B/2-xg)
        +1.0/12.0*(H-tf)*tw*tw*tw
        +(H-tf)*tw*(xg-tw/2.0)*(xg-tw/2.0);
  }
  else if(sect.sform.type==STEEL_TKYOU)
  {
    *Iyy=(tf*B*B*B+(H-tf)*tw*tw*tw)/12.0;
    xg=0.0;
    yc= B/2.0;
    yt=-B/2.0;
  }
  else if(sect.sform.type==STEEL_TWEAK)
  {
    xg=(tw*H*H/2.0+(B-tw)*tf*tf/2.0)/(*A); /*FROM LEFT*/
    yc=H;
    yt=0.0;

    *Iyy=1.0/12.0*tw*H*H*H+tw*H*(H/2-xg)*(H/2-xg)
        +1.0/12.0*(B-tw)*tf*tf*tf
        +(B-tw)*tf*(xg-tf/2.0)*(xg-tf/2.0);
  }
  else if(sect.sform.type==STEEL_CKYOU)
  {
    S=2.0*((B*tf)*B/2.0)+(H-2.0*tf)*tw*(tw/2.0);
    xg=S/(*A);                           /*FROM BOTTOM*/
    yc=B;
    yt=0.0;

    *Iyy=2.0*(1.0/12.0*tw*B*B*B+tw*B*(B/2.0-xg)*(B/2.0-xg))
        +1.0/12.0*(H-2.0*tf)*tw*tw*tw+(H-2.0*tf)*tw*(xg-tw/2.0)*(xg-tw/2.0);
  }
  else if(sect.sform.type==STEEL_CWEAK)
  {
    *Iyy=(H*B*B*B-(H-tw)*(B-2.0*tf)*(B-2.0*tf)*(B-2.0*tf))/12.0;
    xg=0.0;
    yc= B/2.0;
    yt=-B/2.0;
  }

  /*Zyy*/
  Zupper=*Iyy/fabs(yc-xg);
  Zlower=*Iyy/fabs(xg-yt);
  if(Zupper<=Zlower) *Zyy=Zupper;
  else               *Zyy=Zlower;

  /*ixx,iyy*/
  *ixx=sqrt((*Ixx)/(*A));
  *iyy=sqrt((*Iyy)/(*A));

  if(sect.sform.type==STEEL_ANGLE)
  {
    Ixy=fabs(0.25*((xg*xg-(B-xg)*(B-xg))*((-yg+tf)*(-yg+tf)-yg*yg)
                  +(xg*xg-(xg-tw)*(xg-tw))*((H-yg)*(H-yg)-(-yg+tf)*(-yg+tf))));

    if((*Ixx)==(*Iyy)) theta=0.25*PI;
    else               theta=0.5*(atan(2.0*Ixy/((*Iyy)-(*Ixx))));

    Ixx2=(*Ixx)*sin(theta)*sin(theta)+(*Iyy)*cos(theta)*cos(theta)+Ixy*sin(2*theta);
    Iyy2=(*Iyy)*sin(theta)*sin(theta)+(*Ixx)*cos(theta)*cos(theta)-Ixy*sin(2*theta);

    *ixx=sqrt(Ixx2/(*A));
    *iyy=sqrt(Iyy2/(*A));
  }

  /*Jz Cw*/
  if(sect.sform.type==STEEL_HKYOU)
  {
    *Jz=((H-2*tf)*tw*tw*tw+2*B*tf*tf*tf)/3;
    if(Ixx>Iyy) *Cw=H*H**Iyy/4;
    else        *Cw=H*H**Ixx/4;
  }
  else if(sect.sform.type==STEEL_HWEAK)
  {
    *Jz=((B-2*tw)*tf*tf*tf+2*H*tw*tw*tw)/3;
    if(Ixx>Iyy) *Cw=B*B**Iyy/4;
    else        *Cw=B*B**Ixx/4;
  }
  else if(sect.sform.type==STEEL_RPIPE)
  {
    Aw=(H-tf)*(B-tw);
    *Jz=4*Aw*Aw/(2*(H-tf)/tw+2*(B-tw)/tf);
    *Cw=0.0;
  }
  else if(sect.sform.type==STEEL_CPIPE)
  {
    Aw=PI*(H-tf/2)*(H-tf/2)/4;
    *Jz=4*Aw*Aw/(PI*(H-tf/2)/tf);
    *Cw=0.0;
  }
  else if(sect.sform.type==STEEL_ANGLE)
  {
    *Jz=((H-tf)*tw*tw*tw+B*tf*tf*tf)/3;
    *Cw=0.0;
  }
  else if(sect.sform.type==STEEL_TKYOU)
  {
    *Jz=((H-tf)*tw*tw*tw+B*tf*tf*tf)/3;
    if(Ixx>Iyy) *Cw=H*H**Iyy/4;
    else        *Cw=H*H**Ixx/4;
  }
  else if(sect.sform.type==STEEL_TWEAK)
  {
    *Jz=((B-tw)*tf*tf*tf+H*tw*tw*tw)/3;
    if(Ixx>Iyy) *Cw=B*B**Iyy/4;
    else        *Cw=B*B**Ixx/4;
  }
  else if(sect.sform.type==STEEL_CKYOU)
  {
    *Jz=((H-2.0*tf)*tw*tw*tw+2.0*B*tf*tf*tf)/3;
    *Cw=(H*H*B*B*B*tf*(3.0*B*tf+2.0*H*tw))/(12.0*(6.0*B*tf+H*tw));
  }
  else if(sect.sform.type==STEEL_CWEAK)
  {
    *Jz=(2.0*(B-tw)*tf*tf*tf+H*tw*tw*tw)/3;
    *Cw=(H*H*B*B*B*tf*(3.0*B*tf+2.0*H*tw))/(12.0*(6.0*B*tf+H*tw));
  }
  return 1;
}/*alumicoefficients*/


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
    pAi[i]=elem.sect->strnd[i].area[PEND];
    if(axis==0) pYi[i]=elem.sect->strnd[i].y[PEND];
    if(axis==1) pYi[i]=elem.sect->strnd[i].x[PEND];

    Ninit+=elem.sect->strnd[i].Ni[PEND]*1000.0; /*[kgf]*/

    Ap+=pAi[i];
    pAY+=pAi[i]*pYi[i];

    if(pYc<pYi[i]) pYc=pYi[i];
    if(pYt>pYi[i]) pYt=pYi[i];

    pft=m.stft[elem.sect->strnd[i].type]; /*ALLOWABLE PC STRESS*/
    if(pft<(elem.sect->strnd[i].Ni[PEND])/pAi[i])
    {
      if(fout0!=NULL) fprintf(fout0,"TOO MUCH STRESS.\n");
    }
  }

  if(elem.sect->stype==STYPE_PC)
  {
    cfc=m.pcfc; /*+:COMPRESSION*/
    cft=m.pcft; /*-:TENSION*/
  }
  else return 0.0;

  /*
  if(fout0!=NULL)
  {
    fprintf(fout0,"Ac=%.3f Ap=%.3f[cm2]\n",Ac,Ap);
    fprintf(fout0,"cfc=%.3f cft=%.3f [kgf/cm2]\n",cfc,cft);
    fprintf(fout0,"Ninit=%.3f Nd=%.3f[tf]\n",Ninit/1000.0,Nd/1000.0);
  }
  */

  cZ1=1/6.0*cBi[0]*cDi[0]*cDi[0]; /*[cm3]*/
  cZ2=cZ1;

  Ma1=cZ1*(cfc-(Nd+(m.stfact)*Ninit)/Ac);
  Ma2=cZ2*(-cft+(Nd+(m.stfact)*Ninit)/Ac);

  if(Ma1<=Ma2) Ma=Ma1;
  else         Ma=Ma2;

  /*if(fout0!=NULL) fprintf(fout0,"pcMa=%.3f[kgfcm]\n",Ma);*/

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
    pAi[i]=elem.sect->strnd[i].area[PEND];
    if(axis==0) pYi[i]=elem.sect->strnd[i].y[PEND];
    if(axis==1) pYi[i]=elem.sect->strnd[i].x[PEND];

    Ninit+=elem.sect->strnd[i].Ni[PEND]*1000.0; /*[kgf]*/

    Ap+=pAi[i];
    pYg+=pYi[i]*elem.sect->strnd[i].Ni[PEND]*1000.0;

    if(pYc<pYi[i]) pYc=pYi[i];
    if(pYt>pYi[i]) pYt=pYi[i];

    pft=m.stft[elem.sect->strnd[i].type]; /*ALLOWABLE PC STRESS*/
    if(pft<(elem.sect->strnd[i].Ni[PEND])/pAi[i])
    {
      if(fout0!=NULL) fprintf(fout0,"TOO MUCH STRESS.\n");
    }
  }
  pYg/=Ninit;

  if(elem.sect->stype==STYPE_PC)
  {
    cfc=m.pcfc; /*+:COMPRESSION*/
    cft=m.pcft; /*-:TENSION*/
  }
  else return 0.0;

  /*
  if(fout0!=NULL)
  {
    fprintf(fout0,"Ac=%.3f Ap=%.3f[cm2]\n",Ac,Ap);
    fprintf(fout0,"cfc=%.3f cft=%.3f [kgf/cm2]\n",cfc,cft);
    fprintf(fout0,"Ninit=%.3f[tf]\n",Ninit/1000.0);
  }
  */

  cZ1=1/6.0*cBi[0]*cDi[0]*cDi[0]; /*[cm3]*/
  cZ2=cZ1;

  Ma1=cZ1*( cfc-(m.stfact)*Ninit/Ac
               +(m.stfact)*Ninit*(pYg-cYg)/cZ1);
  Ma2=cZ2*(-cft+(m.stfact)*Ninit/Ac
               +(m.stfact)*Ninit*(pYg-cYg)/cZ2);

  if(Ma1<=Ma2) Ma=Ma1;
  else         Ma=Ma2;

  /*if(fout0!=NULL) fprintf(fout0,"pcMa=%.3f[kgfcm]\n",Ma);*/

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
    pAi[i]=elem.sect->strnd[i].area[PMID];
    if(axis==0) pYi[i]=elem.sect->strnd[i].y[PMID];
    if(axis==1) pYi[i]=elem.sect->strnd[i].x[PMID];

    Ninit+=elem.sect->strnd[i].Ni[PMID]*1000.0; /*[kgf]*/

    Ap+=pAi[i];
    pYg+=pYi[i]*elem.sect->strnd[i].Ni[PMID]*1000.0;

    if(pYc<pYi[i]) pYc=pYi[i];
    if(pYt>pYi[i]) pYt=pYi[i];

    pft=m.stft[elem.sect->strnd[i].type]; /*ALLOWABLE PC STRESS*/
    if(pft<(elem.sect->strnd[i].Ni[PMID])/pAi[i])
    {
      if(fout0!=NULL) fprintf(fout0,"TOO MUCH STRESS.\n");
    }
  }
  pYg/=Ninit;

  if(elem.sect->stype==STYPE_PC)
  {
    cfc =m.pcfc;  /*+:COMPRESSION*/
    cft =m.pcft;  /*-:TENSION*/
    cfco=m.pcfco; /*+:COMPRESSION*/
    cfto=m.pcfto; /*-:TENSION*/
  }
  else return 0.0;

  /*
  if(fout0!=NULL)
  {
    fprintf(fout0,"cfc =%.3f cft =%.3f [kgf/cm2]\n",cfc,cft);
    fprintf(fout0,"cfco=%.3f cfto=%.3f [kgf/cm2]\n",cfco,cfto);

    fprintf(fout0,"Ac=%.3f Ap=%.3f[cm2]\n",Ac,Ap);
    fprintf(fout0,"cYg=%.3f pYg=%.3f e=%.3f[cm]\n",cYg,pYg,(cYg-pYg));
    fprintf(fout0,"Ninit=%.3f[tf]\n",Ninit/1000.0);
    fprintf(fout0,"Mdead=%.3f[tfm]\n",Mdead/100000.0);
  }
  */

  cZ1=1/6.0*cBi[0]*cDi[0]*cDi[0]; /*[cm3]*/
  cZ2=cZ1;
  /*if(fout0!=NULL) fprintf(fout0,"Z1=%.0f Z2=%.0f[cm3]\n",cZ1,cZ2);*/

  Ma1=cZ1*( cfto-Ninit/Ac+Ninit*(cYg-pYg)/cZ1);
  Ma2=cZ2*(-cfco+Ninit/Ac+Ninit*(cYg-pYg)/cZ2);

  if(Mdead<Ma1 || Mdead<Ma2)
  {
    if(fout0!=NULL) fprintf(fout0,"\nTOO MUCH STRESS FOR DEAD LOAD.\n");
    return 0.0;
  }

  Ma1=cZ1*( cfc-(m.stfact)*Ninit/Ac
               +(m.stfact)*Ninit*(cYg-pYg)/cZ1);
  Ma2=cZ2*(-cft+(m.stfact)*Ninit/Ac
               +(m.stfact)*Ninit*(cYg-pYg)/cZ2);

  if(Ma1<=Ma2) Ma=Ma1;
  else         Ma=Ma2;

  /*if(fout0!=NULL) fprintf(fout0,"pcMa=%.3f[kgfcm]\n",Ma);*/

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
    pAi[i]=elem.sect->strnd[i].area[PEND];
    if(axis==0) pYi[i]=elem.sect->strnd[i].y[PEND];
    if(axis==1) pYi[i]=elem.sect->strnd[i].x[PEND];

    Ninit+=elem.sect->strnd[i].Ni[PEND]*1000.0; /*[kgf]*/
    rpNy+=m.stftu[elem.sect->strnd[i].type]*pAi[i]; /*[kgf]*/

    Ap+=pAi[i];
    pAY+=pAi[i]*pYi[i];

    if(pYc<pYi[i]) pYc=pYi[i];
    if(pYt>pYi[i]) pYt=pYi[i];

    pft=m.stft[elem.sect->strnd[i].type]; /*ALLOWABLE PC STRESS*/
    if(pft<(elem.sect->strnd[i].Ni[PEND])/pAi[i])
    {
      if(fout0!=NULL) fprintf(fout0,"TOO MUCH STRESS.\n");
    }
  }

  if(elem.sect->stype==STYPE_PC) cfc=m.pcF; /*+:COMPRESSION*/
  else return 0.0;

  /*
  if(fout0!=NULL)
  {
    fprintf(fout0,"Ac=%.3f Ap=%.3f[cm2]\n",Ac,Ap);
    fprintf(fout0,"cfc=%.3f cft=%.3f [kgf/cm2]\n",cfc,cft);
    fprintf(fout0,"Ninit=%.3f Nd=%.3f[tf]\n",Ninit/1000.0,Nd/1000.0);
  }
  */

  k1=1.0; /*EFFECTIVE Fc*/
  k2=0.5; /*CENTER OF STRESS BY COMPRESSED CONCRETE.*/

  /*COMPRESSED CONCRETE DEPTH = NEUTRAL AXIS.*/
  dn=(Nd+(rpNy/2))/(k1*cfc)/cBi[0];
  if(dn>(cDi[0]/2))
  {
    if(fout0!=NULL) fprintf(fout0,"COMPRESSED CONCRETE OVERLAPPED.\n");
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

      /*if(fout0!=NULL)
      {
        fprintf(fout0,"pN=%.3f[tf] dc=%.3f[cm] Mup1=%.3f[tfm]\n",
              pN/1000.0,(cYc-dn-pYi[i]),
              pN*(cYc-dn-pYi[i])/100000.0);
      }*/
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
  if(fout0!=NULL)
  {
    fprintf(fout0,"Nd=%.3f rpNy=%.3f[tf]\n",Nd/1000.0,rpNy/1000.0);
    fprintf(fout0,"dn=%.3f[cm]\n",dn);
    fprintf(fout0,"Muc=%.3f Mup=%.3f[kgfcm]\n",Muc,Mup);

    fprintf(fout0,"pcMu=%.3f[tfm]\n",Mu/100000.0);
  }*/

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
    pAi[i]=elem.sect->strnd[i].area[PEND];
    if(axis==0) pYi[i]=elem.sect->strnd[i].y[PEND];
    if(axis==1) pYi[i]=elem.sect->strnd[i].x[PEND];

    Ap+=pAi[i];
    Ninit+=elem.sect->strnd[i].Ni[PEND]*1000.0; /*[kgf]*/

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
    if(pft<(elem.sect->strnd[i].Ni[PEND])/pAi[i])
    {
      if(fout0!=NULL) fprintf(fout0,"TOO MUCH STRESS.\n");
    }
  }
  pYg/=rpNy; /*CENTER OF YIELD STRESS OF STRANDS.*/

  if(elem.sect->stype==STYPE_PC) cfc=m.pcF; /*+:COMPRESSION*/
  else return 0.0;

  k1=1.0; /*EFFECTIVE Fc*/
  k2=0.5; /*CENTER OF STRESS BY COMPRESSED CONCRETE.*/

  /*COMPRESSED CONCRETE DEPTH = NEUTRAL AXIS.*/
  dn=rpNy/(k1*cfc)/cBi[0];
  if((tensedside==UPPER && dn>(cYg-cYt)) ||
     (tensedside==LOWER && dn>(cYc-cYg)))
  {
    if(fout0!=NULL) fprintf(fout0,"COMPRESSED CONCRETE OVERLAPPED.\n");
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
  if(fout0!=NULL)
  {
    fprintf(fout0,"Ac=%.3f Ap=%.3f[cm2]\n",Ac,Ap);
    fprintf(fout0,"cfc=%.3f cft=%.3f [kgf/cm2]\n",cfc,cft);
    fprintf(fout0,"Ninit=%.3f rpNy=%.3f[tf]\n",
            Ninit/1000.0,rpNy/1000.0);

    fprintf(fout0,"dn=%.3f[cm]\n",dn);
    fprintf(fout0,"Muc=%.3f Mup=%.3f[kgfcm]\n",Muc,Mup);
    fprintf(fout0,"pcMu=%.5f[tfm]\n",Mu/100000.0);
  }
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
    pAi[i]=elem.sect->strnd[i].area[PEND];
    if(axis==SX) pYi[i]=elem.sect->strnd[i].y[PEND];
    if(axis==SY) pYi[i]=elem.sect->strnd[i].x[PEND];

    Ninit+=elem.sect->strnd[i].Ni[PEND]*1000.0; /*[kgf]*/

    Ap+=pAi[i];
    pAY+=pAi[i]*pYi[i];

    if(pYc<pYi[i]) pYc=pYi[i];
    if(pYt>pYi[i]) pYt=pYi[i];

    pft=m.stft[elem.sect->strnd[i].type]; /*ALLOWABLE PC STRESS*/
    if(pft<(elem.sect->strnd[i].Ni[PEND])/pAi[i])
    {
      if(fout0!=NULL) fprintf(fout0,"TOO MUCH STRESS.\n");
    }
  }

  if(elem.sect->stype==STYPE_PC)
  {
    cft=-0.07*(m.pcfc); /*-:TENSION*/
  }
  else return 0.0;

  /*
  if(fout0!=NULL)
  {
    fprintf(fout0,"Ac=%.3f Ap=%.3f[cm2]\n",Ac,Ap);
    fprintf(fout0,"cfc=%.3f cft=%.3f [kgf/cm2]\n",cfc,cft);
    fprintf(fout0,"Ninit=%.3f Nd=%.3f[tf]\n",Ninit/1000.0,Nd/1000.0);
  }
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
  if(fout0!=NULL)
  {
    fprintf(fout0,"cSn=%.3f[cm3] cI=%.3f[cm4]\n",cSn,cI);
    fprintf(fout0,"cft=%.3f sg=%.3f[kgf/cm2]\n",cft,sg);

    fprintf(fout0,"pcQa=%.3f[kgf]\n",Qa);
  }*/
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
  char str[256];
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
    pAi[i]=elem.sect->strnd[i].area[PEND];
    if(axis==SX) pYi[i]=elem.sect->strnd[i].y[PEND];
    if(axis==SY) pYi[i]=elem.sect->strnd[i].x[PEND];

    Ninit+=elem.sect->strnd[i].Ni[PEND]*1000.0; /*[kgf]*/

    Ap += pAi[i];
    pAY += pAi[i]*pYi[i];
    if(pYi[i]<=cYg)
    {
      pAt += pAi[i];
      pYgt += pAi[i]*pYi[i];
    }
    if(pYc<pYi[i]) pYc=pYi[i];
    if(pYt>pYi[i]) pYt=pYi[i];

    pft=m.stft[elem.sect->strnd[i].type]; /*ALLOWABLE PC STRESS*/
    if(pft<(elem.sect->strnd[i].Ni[PEND])/pAi[i])
    {
      if(fout0!=NULL) fprintf(fout0,"TOO MUCH STRESS.\n");
    }
  }
  pYgt/=pAt;

  if(elem.sect->stype==STYPE_PC) cfs=m.pcfsu;
  else return 0.0;

  if(elem.sect->srein[1-axis].n>=1)
  {
    wft=elem.sect->srein[1-axis].F;

    pw=elem.sect->srein[1-axis].area
      *elem.sect->srein[1-axis].n
      /elem.sect->srein[1-axis].pitch
      /cBi[0];
  }
  else
  {
    wft=elem.sect->wF;

    pw=elem.sect->shearrein[1-axis];
  }

  if(pw>0.012) pw=0.012;

  /*
  if(fout0!=NULL)
  {
    fprintf(fout0,"Ac=%.3f Ap=%.3f[cm2]\n",Ac,Ap);
    fprintf(fout0,"cfc=%.3f cft=%.3f [kgf/cm2]\n",cfc,cft);
    fprintf(fout0,"Ninit=%.3f Nd=%.3f[tf]\n",Ninit/1000.0,Nd/1000.0);
  }
  */

  if(Qd==0.0) alpha=1.0;
  else alpha=4.0/(fabs(Md/Qd)/(cYc-pYgt)+1.0);

  if(alpha<1.0)      alpha=1.0;
  else if(2.0<alpha) alpha=2.0;

  sg=(Nd+Ninit)/Ac;

  cJ=7.0/8.0*(cYc-pYgt);

  Qu=(alpha*(cfs+0.1*sg)+0.5*wft*(pw-0.002))*cBi[0]*cJ;
  /*
  if(fout0!=NULL)
  {
    fprintf(fout0,"alpha=%.3f\n",alpha);
    fprintf(fout0,"cfs=%.3f sg=%.3f[kgf/cm2]\n",cfs,sg);
    fprintf(fout0,"wft=%.3f[kgf/cm2] pw=%.3f\n",wft,pw);
    fprintf(fout0,"cB=%.3f dc=%.3f cJ=%.3f[cm]\n",
            cBi[0],(cYc-pYgt),cJ);

    fprintf(fout0,"pcQu=%.3f[kgf]\n",Qu);
  }
  */
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
    pAi[i]=elem.sect->strnd[i].area[PEND];
    if(axis==SX) pYi[i]=elem.sect->strnd[i].y[PEND];
    if(axis==SY) pYi[i]=elem.sect->strnd[i].x[PEND];

    Ninit+=elem.sect->strnd[i].Ni[PEND]*1000.0; /*[kgf]*/

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
    if(pft<(elem.sect->strnd[i].Ni[PEND])/pAi[i])
    {
      if(fout0!=NULL) fprintf(fout0,"TOO MUCH STRESS.\n");
    }
  }
  pYgt/=pAt;

  if(elem.sect->stype==STYPE_PC) cfs=m.pcfsu;
  else return 0.0;

  if(elem.sect->srein[1-axis].n>=1)
  {
    wft=elem.sect->srein[1-axis].F;

    pw=elem.sect->srein[1-axis].area
      *elem.sect->srein[1-axis].n
      /elem.sect->srein[1-axis].pitch
      /cBi[0];
  }
  else
  {
    wft=elem.sect->wF;

    pw=elem.sect->shearrein[1-axis];
  }

  if(pw>0.012) pw=0.012;

  /*
  if(fout0!=NULL)
  {
    fprintf(fout0,"Ac=%.3f Ap=%.3f[cm2]\n",Ac,Ap);
    fprintf(fout0,"cfc=%.3f cft=%.3f [kgf/cm2]\n",cfc,cft);
    fprintf(fout0,"Ninit=%.3f Nd=%.3f[tf]\n",Ninit/1000.0,Nd/1000.0);
  }
  */

  if(tensedside==UPPER) de=pYgt-cYt;
  if(tensedside==LOWER) de=cYc-pYgt;

  alpha=4.0/(fabs(Md/Qd)/de+1.0);
  if(alpha<1.0)      alpha=1.0;
  else if(2.0<alpha) alpha=2.0;

  sg=(Nd+Ninit)/Ac;

  cJ=7.0/8.0*de;

  Qu=(alpha*(cfs+0.1*sg)+0.5*wft*(pw-0.002))*cBi[0]*cJ;
  /*
  if(fout0!=NULL)
  {
    fprintf(fout0,"alpha=%.3f\n",alpha);
    fprintf(fout0,"cfs=%.3f sg=%.3f[kgf/cm2]\n",cfs,sg);
    fprintf(fout0,"wft=%.3f[kgf/cm2] pw=%.3f\n",wft,pw);
    fprintf(fout0,"cB=%.3f dc=%.3f cJ=%.3f[cm]\n",
            cBi[0],(cYc-pYgt),cJ);

    fprintf(fout0,"pcQu=%.3f[kgf]\n",Qu);
  }
  */
  return Qu;
}





