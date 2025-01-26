#include <vcl.h>
#pragma hdrstop
#include <tchar.h>
//---------------------------------------------------------------------
//int WINAPI _tWinMain(HINSTANCE, HINSTANCE, LPTSTR, int)
//{
//	Application->Initialize();
//	Application->Run();
//
//	return 0;
//}
//---------------------------------------------------------------------

/*HOGAN001:CODED BY JUN SATO.*/
/*SRCAL MDIFIED TO RC STANDARD 2010.*/


#include <windows.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#if defined(_WIN64)
	#define GCL_HCURSOR GCLP_HCURSOR
	#define GCL_HBRBACKGROUND GCLP_HBRBACKGROUND
    #define SetClassLong SetClassLongPtr
    #define GetClassLong GetClassLongPtr
#else
#endif





#include "canhead.h"                    /*DEFINITION OF COMMAND ID.*/
#include "archg120-3ktm.c"       /*ANALYSIS NONLINEAR,MATH,VIEWING.*/
#include "corotational.c"
#include "elem.c"
#include "shell.c"
#include "constraint.c"

#include "arclmEigen.c"
#include "assemble.c"
#include "assembleEigen.c"
#include "archg015ktm.c"                  /*ANALYSIS STATIC LINEAR.*/
#include "arclmStatic.c"              /*ANALYSIS STATIC NON-LINEAR.*/
#include "arclmDynamic.c"
#include "bclng021ujok.c"                       /*ELASTIC BUCKLING.*/
#include "vbrat001.c"
#include "qadhg001.c"               /*ANALYSIS BIQUADRATIC ELEMENT.*/
#include "srcal007ujok.c"                     /*DRAW SECTIONS LIST.*/
#include "gnshn103ktm.c"                        /*DYNAMIC ANALYSIS.*/

#include "optimization.c"                 			 /*OPTIMIZATION*/
#include "model001.c"                      /*INP-TO-DXF CONVERSION.*/
#include "arclm002.c"                     /*FOR GREEN-HOUSE PROJECT*/


#define MENUDIALOGS 4 /*MENU DIALOG BOXES.*/
#define INPUTWIDTH 200/*350*/      /*DIALOGUE FONT:"SYSTEM" 18pts.*/
#define BARWIDTH 15                          /*WIDTH OF SCROLL BAR.*/
#define DIALOGFACTOR 2.00/*4.25*/ /*2.75*/

/*#define A4HEIGHT 6371*/ /*2799*/         /*6371:MD4000J 2799:LPT8*/
/*#define A4WIDTH  4799*/ /*1970*/         /*4799:MD4000J 1970:LPT8*/

#define PROPSIMPLE 1
#define PROPDETAIL 2
#define SECTSIMPLE 1
#define SECTDETAIL 2

#define MAXWRITE 65535 /*65535*/ /*FOR CREATING BITMAP FILE*/

#define DEFAULTFILETYPE F_ORGAN
/*#define DEFAULTFILETYPE F_ARCLM*/
/*#define DEFAULTFILETYPE F_SRCAN*/
/*#define DEFAULTFILETYPE F_HOGAN*/
/*#define DEFAULTFILETYPE F_FRAME*/


/*
#define DEFAULTINPUTFILE   "sample01.inp"
#define DEFAULTINPUTFILEZ  "sample01.inl"
#define DEFAULTINPUTFILEX  "sample01.ihx"
#define DEFAULTINPUTFILEY  "sample01.ihy"
#define DEFAULTOUTPUTFILE  "sample01.otp"
#define DEFAULTOUTPUTFILEZ "sample01.otl"
#define DEFAULTOUTPUTFILEX "sample01.ohx"
#define DEFAULTOUTPUTFILEY "sample01.ohy"
#define DEFAULTSECTIONFILE "sample01.lst"
*/

//#define LISTFILENAME   "LunarBaseWhole\\LunarBase1007"
//#define FILENAME       "hogtxt"


#define FILENAME      "ShellTest"
#define FILENAME      "ShellDynaTest"
#define FILENAME  "240701LMBosfullp1b1"
#define FILENAME  "QuadPoly2mmp02b5"
#define FILENAME  "ShellDynaTest"
#define FILENAME  "hingeDynaTest"
#define FILENAME  "LMBexploration"
#define FILENAME  "OSmockup1"

//#define FILENAME      "SnapTest_brate2"
//#define FILENAME      "20231024QuadPolyhedronShellP100B100"
//#define FILENAME      "CurvedCrease231030"
//#define FILENAME      "MPCtest"
//#define FILENAME      "Segment231019"
//#define FILENAME      "ShellTest231021"
//#define FILENAME      "231104QuadPolyhedronShellP100B100"

#define DEFAULTFILETYPE F_ARCLM
//#define DEFAULTINPUTFILE   FILENAME".inp"
#define DEFAULTINPUTFILE   FILENAME".inl"
#define DEFAULTINPUTFILEZ  FILENAME".inl"
#define DEFAULTINPUTFILEX  FILENAME".ihx"
#define DEFAULTINPUTFILEY  FILENAME".ihy"
#define DEFAULTOUTPUTFILE  FILENAME".otl"
#define DEFAULTOUTPUTFILEZ FILENAME".otl"
#define DEFAULTOUTPUTFILEX FILENAME".ohx"
#define DEFAULTOUTPUTFILEY FILENAME".ohy"
#define DEFAULTSECTIONFILE FILENAME".lst"
//#define DEFAULTSECTIONFILE LISTFILENAME".lst"
//#define DEFAULTSECTIONFILE "Buckling_Bezier.fig"

/*FOR ARCLM CHECK.*/
/*
#define DEFAULTINPUTFILE   "Model01.inp"
#define DEFAULTOUTPUTFILE  "Model01.txt"
*/

/*FOR GNSHN QD GRAPH.*/
/*
#define DEFAULTINPUTFILE   "organ10.inp"
#define DEFAULTOUTPUTFILE  "organ10.otp"
*/

int EXEMODE = 0;
int EXITMODE = 0;
/*
EXEMODE
0:DEFAULT
1:ARCLM001 AUTO
2:BCLNG001 AUTO
3:VBRAT001 AUTO
4:ARCLMStatic AUTO
5:ARCLMDynamic AUTO

EXITMODE
0:DO NOT EXIT
1:EXIT AFTER ANALYSIS
*/



# if 1 /*DEFAULT VIEW PARAMETERS.*/
struct viewparam vpdefault={PERSPECTIVE,0,0,          /*TYPE,ORIGIN*/
							0,0,0,                         /*PAGENUM*/
                            248.0,10.0, 150.0,              /*ANGLE*/
							1.0,                            /*GFACT*/
							{0,0,{ 5.0, 5.0, 2.0}},         /*FOCUS*/
							{{0,0,{1000.0, 1000.0, 1000.0}},  /*MAX*/
							 {0,0,{-1000.0, -1000.0, -1000.0}}}, /*MIN*/
							{1.0, 0.01,                  /*AXIS SIZE*/
							 1.0, 0.0, 0.5,       /*FACTORS D,Q,M*/
							 0.5,                           /*PITCH*/
							 5},                            /*HINGE*/
							{5000.0}};       /*ODV=L*/ /*ORGAN TEST*/
# endif

# if 0 /*SHIZUOKA*/
struct viewparam vpdefault={PERSPECTIVE,0,0,
							-60.0,40.0, 500.0,
							1.0,
							{0,0,{70.0,-10.0, 0.0}},
                            {{0,0,{1000.0, 1000.0, 1000.0}},
							 {0,0,{-100.0, -100.0, 6.0}}},
                            {10.0, 0.5,
                             500.0, 0.0, 0.01,
							 0.5,
                             3},
							{3000.0}};

# endif

# if 0
struct viewparam vpdefault={AXONOMETRIC,0,0,
                            -180.0, 0.0, 5000.0,
                            3.93,
                            {0,0,{63.4, 50.0, -10.0}},
                            {{0,0,{ 50.0, 1000.0, 1000.0}},
                             {0,0,{ 30.0, -100.0, -100.0}}},
                            {7.0, 0.5,
                             100.0, 0.0, 0.01,
							 0.5,
                             3},
                            {10000.0}}; /*HAKODATE KESANSHO*/
# endif

struct windowparams wmain={0,0,0,NEUTRAL,"CanMainWin",
						   NULL,NULL,NULL,NULL,
						   NEUTRAL,NEUTRAL,NEUTRAL,
						   NULL};
struct windowparams wmenu={1,0,0,NEUTRAL,"CanMenuWin",
						   NULL,NULL,NULL,NULL,
						   NEUTRAL,NEUTRAL,NEUTRAL,
						   NULL};
struct windowparams wdraw={2,0,0,NEUTRAL,"CanDrawWin",
						   NULL,NULL,NULL,NULL,
						   NEUTRAL,NEUTRAL,NEUTRAL,
						   NULL,0,0};
struct windowparams wmesg={3,0,0,NEUTRAL,"CanMesgWin",
                           NULL,NULL,NULL,NULL,
						   NEUTRAL,NEUTRAL,NEUTRAL,
                           NULL,0,0};
struct windowparams wsurf={4,0,0,NEUTRAL,"CanSurfWin",
                           NULL,NULL,NULL,NULL,
                           NEUTRAL,NEUTRAL,NEUTRAL,
                           NULL,0,0};

struct windowparams wsfrm={5,0,0,NEUTRAL,"CanSframeWin",
                           NULL,NULL,NULL,NULL,
						   NEUTRAL,NEUTRAL,NEUTRAL,
                           NULL,0,0};
struct windowparams wprop={6,0,0,NEUTRAL,"CanPlistWin",
                           NULL,NULL,NULL,NULL,
                           NEUTRAL,NEUTRAL,NEUTRAL,
                           NULL,0,0};
struct windowparams wsect={7,0,0,NEUTRAL,"CanSlistWin",
                           NULL,NULL,NULL,NULL,
                           NEUTRAL,NEUTRAL,NEUTRAL,
                           NULL,0,0};

struct windowparams wsdsp={8,0,0,NEUTRAL,"CanSviewWin",
                           NULL,NULL,NULL,NULL,
						   NEUTRAL,NEUTRAL,NEUTRAL,
                           NULL,0,0};
struct windowparams wweig={9,0,0,NEUTRAL,"CanWeightWin",
                           NULL,NULL,NULL,NULL,
                           NEUTRAL,NEUTRAL,NEUTRAL,
                           NULL,0,0};
struct windowparams wcmqd={10,0,0,NEUTRAL,"CanCmqWin",
                           NULL,NULL,NULL,NULL,
                           NEUTRAL,NEUTRAL,NEUTRAL,
                           NULL,0,0};
struct windowparams whorz={11,0,0,NEUTRAL,"CanHorizonWin",
						   NULL,NULL,NULL,NULL,
                           NEUTRAL,NEUTRAL,NEUTRAL,
						   NULL,0,0};

/*struct arclmframe arci={0,0,
						"\0",
						0,0,0,0,0,
						NULL,NULL,
						NULL,NULL,NULL,NULL,
						NULL,NULL,
						NULL,NULL,NULL,
						0,NULL,NULL,NULL};*/       /*NULL ARCLM FRAME. NO SHELL.*/

struct arclmframe arci={0,0,//code,loff
						"\0",//*appelation
						0,0,0,0,0,0,0,0,//nnode,nelem,nshell,nsect,nreact,nlaps,neig,nconstraint
						0.0,//loadfactor
						NULL,NULL,//*eigenval,**eigenvec
						NULL,NULL,NULL,//*iform,*ddisp,*lambda
						NULL,NULL,NULL,//*melem,*mshell,*dreact
						NULL,//*fsurface
						NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,//*nodes,*ninit,*sects,*elems,*shells,*confs,*constraints;
						NULL,//*nmass
						0,NULL,NULL,NULL//nosect,*iosect,*eosect,*dosect
					   };/*NULL ARCLM FRAME.*/
struct arclmframe arc ={0,0,
						"\0",
						0,0,0,0,0,0,0,0,
						0.0,
						NULL,NULL,
						NULL,NULL,NULL,
						NULL,NULL,NULL,
						NULL,
						NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
						NULL,
						0,NULL,NULL,NULL
					   };/*GLOBAL ARCLM FRAME.*/
struct arclmframe arcx={0,0,
						"\0",
						0,0,0,0,0,0,0,0,
						0.0,
						NULL,NULL,
						NULL,NULL,NULL,
						NULL,NULL,NULL,
						NULL,
						NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
						NULL,
						0,NULL,NULL,NULL
					   };/*ARCLM FRAME FOR X LOAD.*/
struct arclmframe arcy={0,0,
						"\0",
						0,0,0,0,0,0,0,0,
						0.0,
						NULL,NULL,
						NULL,NULL,NULL,
						NULL,NULL,NULL,
						NULL,
						NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
						NULL,
						0,NULL,NULL,NULL
					   };/*ARCLM FRAME FOR Y LOAD.*/

struct biquadframe bqf={0,0,
						"\0",
						0,0,0,0,0,0,0,
						NULL,NULL,NULL,NULL,NULL,
						NULL,NULL,
						NULL,
						NULL,NULL,NULL}; /*GLOBAL BIQUAD FRAME.*/

LRESULT CALLBACK WindowProcedureMain(HWND,UINT,WPARAM,LPARAM);
LRESULT CALLBACK WindowProcedureSheet(HWND,UINT,WPARAM,LPARAM);
LRESULT CALLBACK WindowProcedureBack(HWND,UINT,WPARAM,LPARAM);
LRESULT CALLBACK WindowProcedureDraw(HWND,UINT,WPARAM,LPARAM);
LRESULT CALLBACK WindowProcedureSurf(HWND,UINT,WPARAM,LPARAM);
LRESULT CALLBACK WindowProcedureSframe(HWND,UINT,WPARAM,LPARAM);
LRESULT CALLBACK WindowProcedureProp(HWND,UINT,WPARAM,LPARAM);
LRESULT CALLBACK WindowProcedureSect(HWND,UINT,WPARAM,LPARAM);
LRESULT CALLBACK WindowProcedureSview(HWND,UINT,WPARAM,LPARAM);
LRESULT CALLBACK WindowProcedureMesg(HWND,UINT,WPARAM,LPARAM);
LRESULT CALLBACK WindowProcedureMenu(HWND,UINT,WPARAM,LPARAM);
static BOOL CALLBACK DialogProcMenu1(HWND,UINT,WPARAM,LPARAM);
static BOOL CALLBACK DialogProcMenu2(HWND,UINT,WPARAM,LPARAM);
static BOOL CALLBACK DialogProcMenu3(HWND,UINT,WPARAM,LPARAM);
static BOOL CALLBACK DialogProcMenu4(HWND,UINT,WPARAM,LPARAM);
static BOOL CALLBACK DialogProcText(HWND,UINT,WPARAM,LPARAM);
static BOOL CALLBACK DialogProcIncrement(HWND,UINT,WPARAM,LPARAM);
static BOOL CALLBACK DialogProcPsim(HWND,UINT,WPARAM,LPARAM);
static BOOL CALLBACK DialogProcSsim(HWND,UINT,WPARAM,LPARAM);
static BOOL CALLBACK DialogProcKatakou(HWND,UINT,WPARAM,LPARAM);
static BOOL CALLBACK DialogProcProperty(HWND,UINT,WPARAM,LPARAM);
static BOOL CALLBACK DialogProcElemType(HWND,UINT,WPARAM,LPARAM);
static BOOL CALLBACK DialogProcSectRegist(HWND,UINT,WPARAM,LPARAM);
static BOOL CALLBACK DialogProcConf(HWND,UINT,WPARAM,LPARAM);
static BOOL CALLBACK DialogProcBond(HWND,UINT,WPARAM,LPARAM);
static BOOL CALLBACK EnumChildProcSheet(HWND hwnd,LPARAM lParam);
VOID APIENTRY popupmenudraw(HWND hwnd,POINT pt);

int classnamealloc(struct windowparams *wp,char *cname);
int classdefinition(HINSTANCE hInstance,
                    char *lpszclassname,
                    WNDPROC lpfnwndproc,
                    char *lpszmenuname,
                    int br,int bg,int bb);
HDC createcompati(struct windowparams w);
HWND windowdefinition(struct windowparams *wp,
					  HINSTANCE hInstance,
					  char *lpszclassname,
					  WNDPROC lpfnwndproc,
					  char *lpszmenuname,
					  int br,int bg,int bb,
					  DWORD dwStyle,
					  int x,int y,int w,int h,
					  HWND hparent);
int registwindowparams(struct winparamsreg *wrapp,
                       struct windowparams *wpreg);
void clearwindow(struct windowparams wp);
void DrawSunken(HDC hdc,int left,int top,int right,int bottom);
RECT rectlocation(int l,int t,int r,int b);
struct windowparams *getwindowparams(HWND hwnd);
void vbarlocation(HWND hparent,HWND hchild,RECT *vbar);
void hbarlocation(HWND hparent,HWND hchild,RECT *hbar);
void drawvbar(HDC hdc,long int maxX,long int maxY,
			  struct windowparams *wp);
void drawhbar(HDC hdc,long int maxX,long int maxY,
			  struct windowparams *wp);

char flagswitch(char *flag);

int findtext(HDC hdc,struct snode *strset,int nstr,
             long int mx,long int my);
void drawtexts(HDC hdc,struct snode *strset,int nstr,
               struct viewparam vp);
void printmacro(double xmax,double xmin,
                double ymax,double ymin,
                int idc1,int idg1,int idb1,
                int idc2,int idg2,int idb2,
                int idc3,int idg3,int idb3,
                char *load,char *street);
void printmodelmacro(double xmax,double xmin,
                     double ymax,double ymin,
                     int idv11,int idv12,int idv13,
                     int idv21,int idv22,int idv23,
                     char *street);
void printarclmfigures(FILE *fin,struct viewparam *vp);
void printarclmfigures2(FILE *fin,struct viewparam *vp,int ftype,int mode);
void printsrcanfigures(FILE *fin,struct viewparam *vp,int ftype,int mode);
void printsrcanmax(FILE *fin,struct viewparam *vp,int mode);
void printorganfigures(FILE *fin,struct viewparam *vp,int mode);
void printarclmlastfigures(FILE *fin,struct viewparam *vp,
                           struct arclmframe *af);
void printarclmlastfiguresII(FILE *fin,struct viewparam *vp,
                           struct arclmframe *af);

/*OPTIMIZATION*/
//void conjugategradient(struct organ *org); /*CONJUGATE*/
//void conjugategradientxz(struct organ *org);             /*Fujimoto*/
//void conjugategradientcurve(struct organ *org);          /*Fujimoto*/
//double arclmautomatic(struct organ *org);  /*CONJUGATE*/

/*GLOBAL PARAMETERS*/
HINSTANCE hInstGlobal;                                  /*INSTANCE.*/
struct winparamsreg wrglobal={0,0,NULL};   /*REGISTED WINDOWPARAMS.*/

int globalmode=NEUTRAL;             /*MODES:ARCLM,SRCAL,HOGAN,ORGAN*/
int globalstatus=NEUTRAL;
int globalflag=NEUTRAL;
int globalmessageflag=1;                         /*MESSAGE BOX FLAG*/
int globaldrawflag=1;                            /*DRAW FLAG*/
int prestatus=NEUTRAL;
int createcommand=C_NEUTRAL;              /*FRAME CREATION COMMAND.*/
int createitem=C_NEUTRAL;             /*FRAME CREATION TARGET ITEM.*/
int globalcondensationflag=0;

POINT pbar;                                /*RELATIVE POINT ON BAR.*/

HWND hfitfrom;                 /*HWND SENDING A MESSAGE "FITSHEET".*/

long int initx,inity;       /*INITIAL POSITION OF ROTATE,MOVE,TEXT.*/
struct viewparam initv;

HWND hpopdlg=NULL;                                  /*POPUP DIALOG.*/
int ipopflag;
int idtext;                                   /*ID OF CURRENT TEXT.*/

struct print gprn={0}; /*GLOBAL PRINT PARAMETERS.*/
char gstr[256]; /*STRING BRINGING SECTION CODE.*/
long int gnsect=0,gnprop=0; /*GLOBAL NSECT.*/
struct osect *currentsects=NULL,*gsect=NULL;
struct oprop *currentprops=NULL,*gprop=NULL;
struct osect *addingsect=NULL;
struct oprop *addingprop=NULL;

int srcansects=0;           /*FOR SRCAN LIST.*/
long int codelist[MAXSECT];
struct section *srcanlist=NULL;

int nmultinode=0,nmultielem=0,nmultiwire=0;
int offset;
struct rgbcolor pcolor; /*PREVIOUS COLOR.*/
struct onode *gnode,*pnode,**multinode; /*CURRENT NODE.*/
struct onode nlast; /*LAST NODE OF ELEMENT WHILE CREATING.*/
struct oelem gelem={0,0,
					ROLENULL,WALL,
					0,0,0.0,
					{{0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
					 {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}},
					{{0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
					 {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}},
					NULL,NULL,NULL,NULL,
					{{255,255,255},{255,255,255},{255,255,255}}
				   }; /*CREATING ELEMENT.*/
struct oelem melem={0,0,
					ROLENULL,WALL,
					0,0,0.0,
					{{0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
					 {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}},
					{{0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
					 {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}},
					NULL,NULL,NULL,NULL,
					{{255,255,255},{255,255,255},{255,255,255}}
				   }; /*ADDING ELEMENT WHILE MOVING.*/
struct oelem *pelem,**multielem; /*POINTING CURRENT ELEMENT.*/
struct owire **multiwire;
struct selectrange gselect; /*SELECT RANGE.*/

struct vector gincrement; /*INCREMENT FOR MOVE.*/
int icount; /*COUNT FOUND NODES.*/

/*GLOBAL POLYCURVES.*/
struct polypolycurve gpolypoly={0,NULL};
struct polycurve *cpolycurve; /*CURRENT*/
/*struct polycurve gpolycurve  ={0,0,0,NULL,  0,  0,255};*/
struct polycurve addpolycurve={0,0,0,NULL}; /*ADDING.*/

double globalunit=1.0; /*UNIT FACTOR*/


FILE *globalfile; /*GLOBAL FILE*/

struct organ globalorg;

int WINAPI WinMain(HINSTANCE hInstance,
				   HINSTANCE hPrevInstance,
				   LPSTR lpszCmdLine,
				   int nCmdShow)/*WINAPI APP ENTRY*/
{
  HWND hwnd;
  HACCEL hAccel; /*ACCELERATOR*/
  MSG msg;
  int nc,isd;

  if(!hPrevInstance)                      /*WINDOW CLASS DEFINITION*/
  {
	if(!classdefinition(hInstance,
						wmain.classname,
						WindowProcedureMain,
						"CANMENU",
						130,130,130)) return 0;
  }

  hInstGlobal = hInstance;


  /*main window*/
  hwnd = CreateWindow(wmain.classname,               /*WINDOW CLASS*/
					  "Hoganshi Skeleton",                /*CAPTION*/
					  WS_CLIPCHILDREN |
					  WS_CLIPSIBLINGS |
					  WS_OVERLAPPEDWINDOW,                  /*STYLE*/
					  -10,0,                                  /*X,Y*/
					  1940,1030, /*850,800,*/         /*WIDTH,HEIGHT*/
					  HWND_DESKTOP,                 /*PARENT WINDOW*/
					  NULL,                               /*MENU ID*/
					  hInstance,
					  NULL);
  wmain.hwnd=hwnd;

  registwindowparams(&wrglobal,&wmain);             /*REGISTRATION.*/

  hAccel=LoadAccelerators(hInstance,"CANMENU");          /*SHORTCUT*/

  ShowWindow(hwnd,nCmdShow);
  UpdateWindow(hwnd);             /*DISPATCH "WM_PAINT" IF UPDATED.*/

  while(GetMessage(&msg,NULL,0,0))
  {
	if(wmenu.nchilds<(MENUDIALOGS+2) ||
	   !IsDialogMessage(hpopdlg,&msg))
	{
	  isd=1;
	  for(nc=2;nc<(MENUDIALOGS+2);nc++)
	  {
		if(IsDialogMessage((wmenu.childs+nc)->hwnd,&msg)) isd=0;
	  }
	  for(nc=wprop.nchilds;nc>=3;nc--)
      {
		if(IsDialogMessage((wprop.childs+nc-1)->hwnd,&msg)) isd=0;
	  }
	  for(nc=wsect.nchilds;nc>=3;nc--)
	  {
        if(IsDialogMessage((wsect.childs+nc-1)->hwnd,&msg)) isd=0;
	  }

	  if(isd==1 && !TranslateAccelerator(hwnd,hAccel,&msg))
      {
        TranslateMessage(&msg);               /*KEY INTO CHARACTER.*/
		DispatchMessage(&msg);                  /*DISPATCH MESSAGE.*/
									  /*SEND MESSAGE TO WINDOWPROC.*/
	  }
	}
  }                                       /*REPEAT UNTIL "WM_QUIT".*/

  return msg.wParam;
}/*WinMain*/

LRESULT CALLBACK WindowProcedureMain(HWND hwnd,
									 UINT message,
									 WPARAM wParam,
									 LPARAM lParam)
/*WINDOW PROCEDURE FOR MAIN WINDOW.CREATION,DESTRUCTION OF WINDOWS.*/
{
  HDC hdc;
  HPEN hpen;
  WORD x,y;                                   /*WORD=unsigned short*/
  WPARAM wparam;
  FILE *fin,*fout,*frat;
  size_t n;
  char str[256],doc[80];
  char dir[]=DIRECTORY;
  char s[80],string[1024];
  int i,ip;
  //long int maxX,maxY;
  long int maxX,maxY,ncount;
  long int mw1,mh1,mw2,mh2,mw3,mh3,mw4,mh4;  /*WINDOW,CLIENT SIZE*/



  HFONT hfont;
  double pfactor; /*PRINT FACTOR*/
  struct viewparam vprint;

  /*FOR CREATE BITMAP FILE*/
  HBITMAP hbit,pbit;
  HBRUSH hbrush,pbrush;
  struct print prn;
  int cWidthPels,cHeightPels;
  HDC hDC;
  BITMAP bmp;
  PBITMAPINFO pbmi;
  WORD    cClrBits;
  HANDLE hf;                  /* file handle */
  BITMAPFILEHEADER hdr;       /* bitmap file-header */
  PBITMAPINFOHEADER pbih;     /* bitmap info-header */
  LPBYTE lpBits;              /* memory pointer */
  DWORD dwTotal;              /* total count of bytes */
  DWORD cb;                   /* incremental count of bytes */
  BYTE *hp;                   /* byte pointer */
  DWORD dwTmp;
	/*
    HWND hdesktop;
	BITMAPFILEHEADER fhead;
	BITMAPINFOHEADER ihead;
	RGBTRIPLE rgbtri = {0,0,0};
	COLORREF cr;
	BYTE l='B',h='M';
	FILE *file;
	int j,nblank;
	short int blank = 0;
	int px,mt;
	*/
  int laps;
  double dsafety;

  switch(message)
  {
	case WM_PAINT:
	  DefWindowProc(hwnd,message,wParam,lParam);
	  if(wmenu.hwnd==NULL)
	  {
		wparam = MAKEWPARAM((WORD)IDM_OPTIONS,(WORD)0);
		SendMessage(wmain.hwnd,WM_COMMAND,wparam,0);

		/*Default Window _Drawings ON*/
		wparam = MAKEWPARAM((WORD)IDM_DRAWINGS,(WORD)0);
		SendMessage(wmain.hwnd,WM_COMMAND,wparam,(WPARAM)0);

		/*Error*/
		wparam = MAKEWPARAM((WORD)IDM_ERROR,(WORD)0);
		SendMessage(wmain.hwnd,WM_COMMAND,wparam,(WPARAM)0);



		/*AUTO EXECUTION*/
		FILE *fdata;
		fdata = fopen("config.ini","r");
		if(fdata==NULL)
		{
		  /*READ INPUT AND MAKE MODEL VIEW*/
		  wparam = MAKEWPARAM(IDD_VIEW,0);
		  SendMessage((wmenu.childs+2)->hwnd,WM_COMMAND,wparam,0);
		}
		else
		{
		  char **data;
		  int nstr,pstr;
		  char *filename/*=argv*/;


		  while (1)
		  {
			data=fgetsbrk(fdata,&nstr);

			pstr=0; /*POSITION IN "DATA".*/
			if(nstr==0)  break;

			 while((nstr-pstr)>0)
			 {
				if(nstr-pstr==1) /*POINTING LAST STRING.*/
				{
				  pstr++;
				}
				else
				{
					if(nstr-pstr>=2 && !strcmp(*(data+pstr),"FILENAME"))
					{
					  pstr++;
					  filename=*(data+pstr);
					}
					if(nstr-pstr>=2 && !strcmp(*(data+pstr),"FILETYPE"))
					{
					   pstr++;
					   if (!strcmp(*(data + pstr), "ARCLM"))
					   {
						 (wmenu.childs+2)->vparam.vflag.mv.ftype=F_ARCLM;
					   }
					   else if (!strcmp(*(data + pstr), "ORGAN"))
					   {
						 (wmenu.childs+2)->vparam.vflag.mv.ftype=F_ORGAN;
					   }
					}
					if(nstr-pstr>=2 && !strcmp(*(data+pstr),"EXEMODE"))
					{
					   pstr++;
					   EXEMODE = (int)strtol(*(data + pstr), NULL, 10);
					}
					if(nstr-pstr>=2 && !strcmp(*(data+pstr),"EXITMODE"))
					{
					   pstr++;
					   EXITMODE = (int)strtol(*(data + pstr), NULL, 10);
					}
					else
					{
						pstr++;
					}
				}
			 }
		  }
		  fclose(fdata);

		  /*INPUT FILE SETTING*/
		  strcpy((wdraw.childs+1)->inpfile,filename);
		  if((wmenu.childs+2)->vparam.vflag.mv.ftype==F_ORGAN)
		  {
			strcat((wdraw.childs+1)->inpfile,".inp");
		  }
		  else if((wmenu.childs+2)->vparam.vflag.mv.ftype==F_ARCLM)
		  {
			strcat((wdraw.childs+1)->inpfile,".inl");
		  }
		  strcpy((wdraw.childs+1)->inpfilez,filename);
		  strcat((wdraw.childs+1)->inpfilez,".inl");
		  strcpy((wdraw.childs+1)->inpfilex,filename);
		  strcat((wdraw.childs+1)->inpfilex,".ihx");
		  strcpy((wdraw.childs+1)->inpfiley,filename);
		  strcat((wdraw.childs+1)->inpfiley,".ihy");

		  /*OUTPUT FILE SETTING*/
		  strcpy((wdraw.childs+1)->otpfile,filename);
		  strcat((wdraw.childs+1)->otpfile,".otp");
		  strcpy((wdraw.childs+1)->otpfilez,filename);
		  strcat((wdraw.childs+1)->otpfilez,".otl");
		  strcpy((wdraw.childs+1)->otpfilex,filename);
		  strcat((wdraw.childs+1)->otpfilex,".ohx");
		  strcpy((wdraw.childs+1)->otpfiley,filename);
		  strcat((wdraw.childs+1)->otpfiley,".ohy");

		  /*LST FILE SETTING*/
		  strcpy((wdraw.childs+1)->sctfile,filename);
		  strcat((wdraw.childs+1)->sctfile,".lst");

		  SetDlgItemText((wmenu.childs+2)->hwnd,ID_INPUTFILE,
					 (wdraw.childs+1)->inpfile);
		  SetDlgItemText((wmenu.childs+2)->hwnd,ID_INPUTFILEZ,
					 (wdraw.childs+1)->inpfilez);
		  SetDlgItemText((wmenu.childs+2)->hwnd,ID_INPUTFILEX,
					 (wdraw.childs+1)->inpfilex);
		  SetDlgItemText((wmenu.childs+2)->hwnd,ID_INPUTFILEY,
					 (wdraw.childs+1)->inpfiley);

		  SetDlgItemText((wmenu.childs+2)->hwnd,ID_OUTPUTFILE,
					 (wdraw.childs+1)->otpfile);
		  SetDlgItemText((wmenu.childs+2)->hwnd,ID_OUTPUTFILEZ,
					 (wdraw.childs+1)->otpfilez);
		  SetDlgItemText((wmenu.childs+2)->hwnd,ID_OUTPUTFILEX,
					 (wdraw.childs+1)->otpfilex);
		  SetDlgItemText((wmenu.childs+2)->hwnd,ID_OUTPUTFILEY,
					 (wdraw.childs+1)->otpfiley);
		  SetDlgItemText((wmenu.childs+2)->hwnd,ID_SECTIONFILE,
					 (wdraw.childs+1)->sctfile);

		  /*READ INPUT AND MAKE MODEL VIEW*/
		  wparam = MAKEWPARAM(IDD_VIEW,0);
		  SendMessage((wmenu.childs+2)->hwnd,WM_COMMAND,wparam,0);


#if 1
			/*ARCLM001_LXY AUTO EXECUTION*/
			if((wmenu.childs+2)->vparam.vflag.mv.ftype==F_ORGAN && EXEMODE==1)
			{
			  globalmessageflag=0;

			  /*EXTRACT ARCLM*/
			  extractarclmfromorgan(&((wdraw.childs+1)->org),
									&arc,&arcx,&arcy);

			  (wmenu.childs+2)->vparam.vflag.mv.ftype=F_ARCLM;

			  /*SAVE AS ARCLM*/
			  saveasarclm((wdraw.childs+1)->inpfilez,&arc);
			  saveasarclm((wdraw.childs+1)->inpfilex,&arcx);
			  saveasarclm((wdraw.childs+1)->inpfiley,&arcy);

			  /*ANALYSIS*/
			  getviewparam((wmenu.childs+2)->hwnd,
						   &((wdraw.childs+1)->vparam));
			  (wdraw.childs+1)->vparam.vflag.ev.deformation=0;

			  int idinputs[3],idoutputs[3];
			  struct arclmframe *arcs[3];

			  arcs[0]=&arc;
			  arcs[1]=&arcx;
			  arcs[2]=&arcy;

			  idinputs[0]=ID_INPUTFILEZ;
			  idinputs[1]=ID_INPUTFILEX;
			  idinputs[2]=ID_INPUTFILEY;

			  idoutputs[0]=ID_OUTPUTFILEZ;
			  idoutputs[1]=ID_OUTPUTFILEX;
			  idoutputs[2]=ID_OUTPUTFILEY;

			  arclm001_lxy(arcs,idinputs,idoutputs);

			  /*SRCAN*/
			  n=strcspn((wdraw.childs+1)->sctfile,".");
			  strncpy(str,(wdraw.childs+1)->sctfile,n);
			  str[n]='\0';

			  i=srcan001(str);

			  if(i==0) MessageBox(NULL,"Failed.","SRCAN001",MB_OK);

			  strcpy((wdraw.childs+1)->otpfile,str);
			  strcat((wdraw.childs+1)->otpfile,".rat");
			  SetDlgItemText((wmenu.childs+2)->hwnd,ID_OUTPUTFILE,
							 (wdraw.childs+1)->otpfile);

			  fout=fopen((wdraw.childs+1)->otpfile,"r");
			  if(fout==NULL) break;
			  readsrcanrate(fout,&arc);

			  /*SAFETY FLAG*/
			  (wdraw.childs+1)->vparam.vflag.ev.srcancolor=1;

			  clearwindow(*(wdraw.childs+1));
			  drawarclmframe((wdraw.childs+1)->hdcC,
							 (wdraw.childs+1)->vparam,arc,0,ONSCREEN);
			  SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);

			  /*ABORT*/
			  if(EXITMODE == 1)
			  {
				fclose(fin);
				exit(1);
			  }
			}
#endif
			else if((wmenu.childs+2)->vparam.vflag.mv.ftype==F_ORGAN && EXEMODE==2)
			{
			  globalmessageflag=0;

			  /*EXTRACT ARCLM*/
			  extractarclmfromorgan(&((wdraw.childs+1)->org),
									&arc,&arcx,&arcy);

			  (wmenu.childs+2)->vparam.vflag.mv.ftype=F_ARCLM;

			  /*SAVE AS ARCLM*/
			  saveasarclm((wdraw.childs+1)->inpfilez,&arc);

			  /*ANALYSIS*/
			  getviewparam((wmenu.childs+2)->hwnd,
						   &((wdraw.childs+1)->vparam));
			  (wdraw.childs+1)->vparam.vflag.ev.deformation=0;

			  arclm001(&arc ,ID_INPUTFILEZ,ID_OUTPUTFILEZ);

			  /*OPEN analysisdata.txt*/
			  int solver = 0;
			  int neig = 1;

			  bclng001(&arc/*,solver,neig*/);

			  /*REDRAW MODEL*/
			  (wdraw.childs+1)->vparam.vflag.ev.deformation=1;
			  clearwindow(*(wdraw.childs+1));
			  drawarclmframe((wdraw.childs+1)->hdcC,
							 (wdraw.childs+1)->vparam,arc,0,ONSCREEN);
			  SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);


			  SendMessage((wmenu.childs+2)->hwnd,WM_INITDIALOG,0,0);
			  globalmessageflag=1;
			  globaldrawflag=1;

			  /*ABORT*/
			  if(EXITMODE == 1)
			  {
				fclose(fin);
				exit(1);
			  }
			}
			/*VBRAT001 AUTO EXECUTION*/
			else if((wmenu.childs+2)->vparam.vflag.mv.ftype==F_ARCLM && EXEMODE==3)
			{
			  globalmessageflag=0;

			  vbrat001(&arc);

			  clearwindow(*(wdraw.childs+1));
			  drawarclmframe((wdraw.childs+1)->hdcC,
							 (wdraw.childs+1)->vparam,arc,0,ONSCREEN);
			  SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);

			  /*ABORT*/
			  if(EXITMODE == 1)
			  {
				fclose(fin);
				exit(1);
			  }
			}
			/*ARCLMCR AUTO EXECUTION*/
			else if((wmenu.childs+2)->vparam.vflag.mv.ftype==F_ARCLM && EXEMODE==4)
			{
			  globalmessageflag=0;

			  arclmStatic(&arc);

			  /*SRCAN*/
			  n=strcspn((wdraw.childs+1)->sctfile,".");
			  strncpy(str,(wdraw.childs+1)->sctfile,n);
			  str[n]='\0';

			  i=srcan001(str);

			  if(i==0)
			  {
				if(EXITMODE == 1)
				{
				  fclose(fin);
				  exit(1);
				}
				else
				{
				  MessageBox(NULL,"Failed.","SRCAN001",MB_OK);
				}
			  }

			  strcpy((wdraw.childs+1)->otpfile,str);
			  strcat((wdraw.childs+1)->otpfile,".rat");
			  SetDlgItemText((wmenu.childs+2)->hwnd,ID_OUTPUTFILE,
							 (wdraw.childs+1)->otpfile);

			  fout=fopen((wdraw.childs+1)->otpfile,"r");
			  if(fout==NULL) break;
			  readsrcanrate(fout,&arc);

			  clearwindow(*(wdraw.childs+1));
			  drawarclmframe((wdraw.childs+1)->hdcC,
							 (wdraw.childs+1)->vparam,arc,0,ONSCREEN);
			  SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);

			  /*ABORT*/
			  if(EXITMODE == 1)
			  {
				fclose(fin);
				exit(1);
			  }
			}
			else if((wmenu.childs+2)->vparam.vflag.mv.ftype==F_ARCLM && EXEMODE==5)
			{
			  globalmessageflag=0;

			  arclmDynamic(&arc);

			  clearwindow(*(wdraw.childs+1));
			  drawarclmframe((wdraw.childs+1)->hdcC,
							 (wdraw.childs+1)->vparam,arc,0,ONSCREEN);
			  SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);

			  /*ABORT*/
			  if(EXITMODE == 1)
			  {
				fclose(fin);
				exit(1);
			  }
			}

		}
	  }
	  break;

	case WM_SIZE:
	  DefWindowProc(hwnd,message,wParam,lParam);
	  wparam = MAKEWPARAM((WORD)IDM_FITSHEET,(WORD)0);
	  if(wmenu.hwnd!=NULL)
	  {
		SendMessage((wmenu.childs+1)->hwnd,WM_COMMAND,wparam,0);
	  }
	  if(wdraw.hwnd!=NULL)
	  {
		SendMessage((wdraw.childs+1)->hwnd,WM_COMMAND,wparam,0);
	  }
	  break;

	case WM_COMMAND:
	  switch(LOWORD(wParam))
	  {
		case IDM_END:                            /*END APPLICATION.*/
		  PostQuitMessage(0);                   /*CREATE "WM_QUIT".*/
		  break;

		case IDM_DRAWINGS:                      /*Drawing Window*/
		  if(wdraw.hwnd==NULL &&
			 wmenu.hwnd!=NULL) /*CREATION.*/
		  {
			getclientsize(hwnd,&maxX,&maxY);
			getwindowsize(wmenu.hwnd,&mw1,&mh1);

			windowdefinition(&wdraw,
							 hInstGlobal,
							 NULL,
							 WindowProcedureSheet,
							 NULL,
							 190,190,190,
							 WS_CLIPCHILDREN |
							 WS_CLIPSIBLINGS |
							 WS_CHILD |
							 WS_THICKFRAME |
							 WS_VISIBLE,
							 mw1,0,maxX-mw1-500+20,maxY,
							 wmain.hwnd);

			getclientsize(wdraw.hwnd,&maxX,&maxY);

			wdraw.nchilds=2;
			wdraw.childs=(struct windowparams *)
						 malloc(2*sizeof(struct windowparams));

			windowdefinition(wdraw.childs+0,
							 hInstGlobal,
							 "CanDrawBak",
							 WindowProcedureBack,
							 NULL,
							 190,190,190,
							 WS_CLIPCHILDREN |
							 WS_CLIPSIBLINGS |
                             WS_CHILD |
							 WS_VISIBLE,
							 1,1,maxX-BARWIDTH/*-5*/,maxY-BARWIDTH/*-5*/,
							 wdraw.hwnd);

			windowdefinition(wdraw.childs+1,
							 hInstGlobal,
							 "CanDrawDsp",
							 WindowProcedureDraw,
							 NULL,
							 0,0,0,
							 /*WS_CLIPCHILDREN |*/
							 WS_CLIPSIBLINGS |
							 WS_CHILD |
							 WS_VISIBLE,
							 0,0,maxX-BARWIDTH/*-5+590*/,maxY-BARWIDTH/*-5+200*/,   //size
							 //0,0,1200,900,
							 //0,0,900,1200,
							 (wdraw.childs+0)->hwnd);
			setfontformat((wdraw.childs+1)->hdcC,30,12,	/*,15,6,*/
						  "ＭＳ 明朝",0,0,255);

			initializeorganization(&((wdraw.childs+1)->org));
			(wdraw.childs+1)->org.opaque=0.5;

			sprintf((wdraw.childs+1)->pagetitle,"\0");
			sprintf((wdraw.childs+1)->title,"\0");
			sprintf((wdraw.childs+1)->inpfile, DEFAULTINPUTFILE);
            sprintf((wdraw.childs+1)->inpfilez,DEFAULTINPUTFILEZ);
			sprintf((wdraw.childs+1)->inpfilex,DEFAULTINPUTFILEX);
			sprintf((wdraw.childs+1)->inpfiley,DEFAULTINPUTFILEY);
            sprintf((wdraw.childs+1)->otpfile, DEFAULTOUTPUTFILE);
			sprintf((wdraw.childs+1)->otpfilez,DEFAULTOUTPUTFILEZ);
            sprintf((wdraw.childs+1)->otpfilex,DEFAULTOUTPUTFILEX);
			sprintf((wdraw.childs+1)->otpfiley,DEFAULTOUTPUTFILEY);
			sprintf((wdraw.childs+1)->sctfile, DEFAULTSECTIONFILE);

            createviewdata(&vpdefault);
            (wdraw.childs+1)->vparam=vpdefault;
            /*(wdraw.childs+1)->vparam.type=AXONOMETRIC;*/
			/*(wdraw.childs+1)->vparam.type=PERSPECTIVE;*/
			/*getviewparam((wmenu.childs+2)->hwnd,
						 &((wdraw.childs+1)->vparam));*/

			(wdraw.childs+1)->vparam.vflag.pv.printrange =0;   //by mihara
			(wdraw.childs+1)->vparam.vflag.pv.a4tate     =1;   //by mihara
			(wdraw.childs+1)->vparam.vflag.pv.a4yoko     =0;   //by mihara
			(wdraw.childs+1)->vparam.vflag.pv.a3tate     =0;   //by mihara
			(wdraw.childs+1)->vparam.vflag.pv.a3yoko     =0;   //by mihara

			getclientsize((wdraw.childs+1)->hwnd,&maxX,&maxY);
			(wdraw.childs+1)->vparam.Xo=(int)(maxX/2);
			(wdraw.childs+1)->vparam.Yo=(int)(maxY/2);
//            (wdraw.childs+1)->vparam.dparam.gaxis=3.0;
			(wdraw.childs+1)->vparam.chapter=0;
			(wdraw.childs+1)->vparam.section=0;
			(wdraw.childs+1)->vparam.subsection=0;


			(wdraw.childs+1)->vparam.vflag.axis=1;

			(wdraw.childs+1)->vparam.vflag.mv.inputfile =1;  //ujioka

			(wdraw.childs+1)->vparam.vflag.mv.outputfile=1;  //ujioka

			(wdraw.childs+1)->vparam.vflag.mv.pagetitle=0;
			(wdraw.childs+1)->vparam.vflag.mv.pagenum=0;
			(wdraw.childs+1)->vparam.vflag.mv.title=1;       //ujioka
			(wdraw.childs+1)->vparam.vflag.mv.view =1;       //ujioka

			(wdraw.childs+1)->vparam.vflag.nv.code=0;        //ujioka nodecode

			(wdraw.childs+1)->vparam.vflag.nv.d[0]=0;
			(wdraw.childs+1)->vparam.vflag.nv.d[1]=0;
			(wdraw.childs+1)->vparam.vflag.nv.d[2]=0;

			(wdraw.childs+1)->vparam.vflag.nv.confs[0]=0;
			(wdraw.childs+1)->vparam.vflag.nv.confs[1]=0;
			(wdraw.childs+1)->vparam.vflag.nv.confs[2]=0;
			(wdraw.childs+1)->vparam.vflag.nv.confs[3]=0;
			(wdraw.childs+1)->vparam.vflag.nv.confs[4]=0;
			(wdraw.childs+1)->vparam.vflag.nv.confs[5]=0;

			(wdraw.childs+1)->vparam.vflag.nv.conffig =0;

			(wdraw.childs+1)->vparam.vflag.nv.loads[0]=0;
			(wdraw.childs+1)->vparam.vflag.nv.loads[1]=0;
			(wdraw.childs+1)->vparam.vflag.nv.loads[2]=0;
			(wdraw.childs+1)->vparam.vflag.nv.loads[3]=0;
			(wdraw.childs+1)->vparam.vflag.nv.loads[4]=0;
			(wdraw.childs+1)->vparam.vflag.nv.loads[5]=0;

			(wdraw.childs+1)->vparam.vflag.nv.disps[0]=0;
			(wdraw.childs+1)->vparam.vflag.nv.disps[1]=0;
			(wdraw.childs+1)->vparam.vflag.nv.disps[2]=0;

			(wdraw.childs+1)->vparam.vflag.nv.react[0]=0;
			(wdraw.childs+1)->vparam.vflag.nv.react[1]=0;
			(wdraw.childs+1)->vparam.vflag.nv.react[2]=0;
			(wdraw.childs+1)->vparam.vflag.nv.react[3]=0;
			(wdraw.childs+1)->vparam.vflag.nv.react[4]=0;
			(wdraw.childs+1)->vparam.vflag.nv.react[5]=0;

			(wdraw.childs+1)->vparam.vflag.nv.mcircle=0;
			(wdraw.childs+1)->vparam.vflag.nv.mvalue =0.0000005; //UJIOKA

			(wdraw.childs+1)->vparam.vflag.ev.ecircle=0;
			(wdraw.childs+1)->vparam.vflag.ev.plasticenergy=0;
			(wdraw.childs+1)->vparam.vflag.ev.evalue =0;

			(wdraw.childs+1)->vparam.vflag.ev.code       =0;
			(wdraw.childs+1)->vparam.vflag.ev.sectioncode=0;
			(wdraw.childs+1)->vparam.vflag.ev.sectionshape=0;  //honda
			(wdraw.childs+1)->vparam.vflag.ev.cmqline    =0;
			(wdraw.childs+1)->vparam.vflag.ev.cmqcheck   =0;   //CMQ CHECK
			(wdraw.childs+1)->vparam.vflag.ev.axis       =0;
			(wdraw.childs+1)->vparam.vflag.ev.hinge      =1;
			(wdraw.childs+1)->vparam.vflag.ev.deformation=1;
			(wdraw.childs+1)->vparam.vflag.ev.srcanrate  =0;
			(wdraw.childs+1)->vparam.vflag.ev.srcancolor =0;
			(wdraw.childs+1)->vparam.vflag.ev.srcanmax   =0;

			for(i=0;i<=7;i++)
			{
			 (wdraw.childs+1)->vparam.vflag.ev.etype[i]=1;

             (wdraw.childs+1)->vparam.vflag.ev.stress[i][0]=0;
			 (wdraw.childs+1)->vparam.vflag.ev.stress[i][1]=0;
			 (wdraw.childs+1)->vparam.vflag.ev.stress[i][2]=0;
			 (wdraw.childs+1)->vparam.vflag.ev.stress[i][3]=0;
			 (wdraw.childs+1)->vparam.vflag.ev.stress[i][4]=0;
			 (wdraw.childs+1)->vparam.vflag.ev.stress[i][5]=0;
            }

            (wmenu.childs+2)->vparam.vflag.mv.draw=1;
            SendMessage((wmenu.childs+2)->hwnd,WM_INITDIALOG,0,0);
            SendMessage((wmenu.childs+4)->hwnd,WM_INITDIALOG,0,0);
		  }
		  else if(wdraw.hwnd!=NULL) /*DESTRUCTION.*/
		  {
            for(i=1;i<=wdraw.nchilds;i++)
            {
              ReleaseDC((wdraw.childs+i-1)->hwnd,
						(wdraw.childs+i-1)->hdcB);
              ReleaseDC((wdraw.childs+i-1)->hwnd,
						(wdraw.childs+i-1)->hdcC);
            }
            ReleaseDC(wdraw.hwnd,wdraw.hdcB);
			ReleaseDC(wdraw.hwnd,wdraw.hdcC);
            DestroyWindow(wdraw.hwnd);

            free(wdraw.childs);
			wdraw.nchilds=0;
			wdraw.hwnd = NULL;

			(wmenu.childs+2)->vparam.vflag.mv.draw=0;
			SendMessage((wmenu.childs+2)->hwnd,WM_INITDIALOG,0,0);
		  }
		  break;

		case IDM_ERROR:
		  if(wmesg.hwnd==NULL &&
			 wmenu.hwnd!=NULL) /*CREATION.*/
		  {
			windowdefinition(&wmesg,
							 hInstGlobal,
							 NULL,
							 WindowProcedureSheet,
							 NULL,
							 190,190,190,
							 WS_CLIPCHILDREN |
							 WS_CLIPSIBLINGS |
							 /*WS_CHILD |*/
							 WS_POPUP |
							 WS_CAPTION |
							 WS_THICKFRAME |
							 WS_VISIBLE,
							 /*2000,500,500,1000, */	1930-500,40,500,980,
							 wmain.hwnd);

			wmesg.nchilds=2;
			wmesg.childs=(struct windowparams *)
						 malloc(2*sizeof(struct windowparams));

			getclientsize(wmesg.hwnd,&maxX,&maxY);
			windowdefinition(wmesg.childs+0,
							 hInstGlobal,
							 "CanMesgBak",
							 WindowProcedureBack,
							 NULL,
							 0,0,0,
							 WS_CLIPCHILDREN |
							 WS_CLIPSIBLINGS |
							 WS_CHILD |
							 WS_VISIBLE,
							 1,1,maxX-BARWIDTH-5,maxY-BARWIDTH-5,
							 wmesg.hwnd);

			getclientsize((wmesg.childs+0)->hwnd,&maxX,&maxY);
			windowdefinition(wmesg.childs+1,
							 hInstGlobal,
							 "CanMesgDsp",
							 WindowProcedureMesg,
                             NULL,
                             0,0,0,
                             /*WS_CLIPCHILDREN |*/
							 WS_CLIPSIBLINGS |
                             WS_CHILD |
                             WS_VISIBLE,
							 0,0,500,2000,
                             (wmesg.childs+0)->hwnd);

            (wmesg.childs+1)->tx=0;
            (wmesg.childs+1)->ty=0;

			setfontformat((wmesg.childs+1)->hdcC,15,6,
                          "MS Mincho",255,255,255);

            (wmenu.childs+2)->vparam.vflag.mv.error=1;
			SendMessage((wmenu.childs+2)->hwnd,WM_INITDIALOG,0,0);
          }
          else if(wmesg.hwnd!=NULL) /*DESTRUCTION.*/
          {
            for(i=1;i<=wmesg.nchilds;i++)
            {
              ReleaseDC((wmesg.childs+i-1)->hwnd,
                        (wmesg.childs+i-1)->hdcB);
			  ReleaseDC((wmesg.childs+i-1)->hwnd,
                        (wmesg.childs+i-1)->hdcC);
            }
            ReleaseDC(wmesg.hwnd,wmesg.hdcB);
			ReleaseDC(wmesg.hwnd,wmesg.hdcC);
			DestroyWindow(wmesg.hwnd);

            free(wmesg.childs);
            wmesg.nchilds=0;
            wmesg.hwnd = NULL;

            (wmenu.childs+2)->vparam.vflag.mv.error=0;
            SendMessage((wmenu.childs+2)->hwnd,WM_INITDIALOG,0,0);
          }
          break;

        case IDM_SURFACE:
		  if(wsurf.hwnd==NULL &&
			 wmenu.hwnd!=NULL) /*CREATION.*/
		  {
			windowdefinition(&wsurf,
							 hInstGlobal,
							 NULL,
							 WindowProcedureSheet,
							 NULL,
							 190,190,190,
                             WS_CLIPCHILDREN |
                             WS_CLIPSIBLINGS |
                             /*WS_CHILD |*/
                             WS_POPUP |
                             WS_CAPTION |
                             WS_THICKFRAME |
                             WS_VISIBLE,
                             600,800,400,600,
                             wmain.hwnd);

            wsurf.nchilds=2;
            wsurf.childs=(struct windowparams *)
                         malloc(2*sizeof(struct windowparams));

			getclientsize(wsurf.hwnd,&maxX,&maxY);
			windowdefinition(wsurf.childs+0,
							 hInstGlobal,
							 "CanSurfBak",
							 WindowProcedureBack,
                             NULL,
                             0,0,0,
                             WS_CLIPCHILDREN |
                             WS_CLIPSIBLINGS |
                             WS_CHILD |
                             WS_VISIBLE,
                             1,1,maxX-BARWIDTH-5,maxY-BARWIDTH-5,
                             wsurf.hwnd);

            windowdefinition(wsurf.childs+1,
                             hInstGlobal,
                             "CanSurfDsp",
                             WindowProcedureSurf,
                             NULL,
                             0,0,0,
                             /*WS_CLIPCHILDREN |*/
                             WS_CLIPSIBLINGS |
                             WS_CHILD |
							 WS_VISIBLE,
                             0,0,500,500,
                             (wsurf.childs+0)->hwnd);
			setfontformat((wsurf.childs+1)->hdcC,15,6,
						  "MS Mincho",100,255,255);

            getclientsize((wsurf.childs+1)->hwnd,&maxX,&maxY);

            /*SURFACE VIEWPOINT*/
            (wsurf.childs+1)->vparam.type=PERSPECTIVE;
            (wsurf.childs+1)->vparam.gfactor=1.0;
            (wsurf.childs+1)->vparam.focus.d[0]=0.0;
            (wsurf.childs+1)->vparam.focus.d[1]=0.0;
            (wsurf.childs+1)->vparam.focus.d[2]=0.0;
            (wsurf.childs+1)->vparam.theta=50.0;
            (wsurf.childs+1)->vparam.phi=20.0;
            (wsurf.childs+1)->vparam.r=10.0;
			(wsurf.childs+1)->vparam.odv=1000.0;
            (wsurf.childs+1)->vparam.Xo=(int)(maxX/2);
            (wsurf.childs+1)->vparam.Yo=(int)(maxY/2);
            (wsurf.childs+1)->vparam.dparam.gaxis=2.0;
			createviewdata(&((wsurf.childs+1)->vparam));

            (wsurf.childs+1)->lstatus=ROTATE;

            (wmenu.childs+2)->vparam.vflag.mv.surface=1;
            SendMessage((wmenu.childs+2)->hwnd,WM_INITDIALOG,0,0);
          }
          else if(wsurf.hwnd!=NULL) /*DESTRUCTION.*/
          {
            for(i=1;i<=wsurf.nchilds;i++)
            {
              ReleaseDC((wsurf.childs+i-1)->hwnd,
                        (wsurf.childs+i-1)->hdcB);
              ReleaseDC((wsurf.childs+i-1)->hwnd,
                        (wsurf.childs+i-1)->hdcC);
            }
            ReleaseDC(wsurf.hwnd,wsurf.hdcB);
            ReleaseDC(wsurf.hwnd,wsurf.hdcC);
            DestroyWindow(wsurf.hwnd);

            free(wsurf.childs);
			wsurf.nchilds=0;
            wsurf.hwnd = NULL;

			(wmenu.childs+2)->vparam.vflag.mv.surface=0;
            SendMessage((wmenu.childs+2)->hwnd,WM_INITDIALOG,0,0);
          }
          break;

        case IDM_SECTIONFRAME:
          if(wsfrm.hwnd==NULL &&
             wmenu.hwnd!=NULL) /*CREATION.*/
          {
            windowdefinition(&wsfrm,
                             hInstGlobal,
                             NULL,
                             WindowProcedureSframe,
                             "HOGSECT",
							 130,130,130,
                             WS_CLIPCHILDREN |
                             WS_CLIPSIBLINGS |
                             /*WS_CHILD |*/
							 WS_POPUP |
                             WS_CAPTION |
                             WS_THICKFRAME |
                             WS_VISIBLE,
                             40,100,1100,1200,
                             wmain.hwnd);

            (wmenu.childs+2)->vparam.vflag.mv.sectlist=1;
            SendMessage((wmenu.childs+2)->hwnd,WM_INITDIALOG,0,0);
          }
          else if(wsfrm.hwnd!=NULL) /*DESTRUCTION.*/
          {
            if(wsect.hwnd!=NULL)
            {
              for(i=0;i<wsect.nchilds;i++)
              {
                ReleaseDC((wsect.childs+i)->hwnd,
                          (wsect.childs+i)->hdcB);
                ReleaseDC((wsect.childs+i)->hwnd,
                          (wsect.childs+i)->hdcC);
              }
			  ReleaseDC(wsect.hwnd,wsect.hdcB);
              ReleaseDC(wsect.hwnd,wsect.hdcC);
			  DestroyWindow(wsect.hwnd);

              free(wsect.childs);
              wsect.nchilds=0;
              wsect.hwnd = NULL;
              gnsect=0;
            }
            if(wsdsp.hwnd!=NULL)
            {
              for(i=0;i<wsdsp.nchilds;i++)
              {
                ReleaseDC((wsdsp.childs+i)->hwnd,
                          (wsdsp.childs+i)->hdcB);
                ReleaseDC((wsdsp.childs+i)->hwnd,
                          (wsdsp.childs+i)->hdcC);
              }
              ReleaseDC(wsdsp.hwnd,wsdsp.hdcB);
              ReleaseDC(wsdsp.hwnd,wsdsp.hdcC);
              DestroyWindow(wsdsp.hwnd);

              free(wsdsp.childs);
              wsdsp.nchilds=0;
              wsdsp.hwnd = NULL;
            }

            for(i=0;i<wsfrm.nchilds;i++)
            {
              ReleaseDC((wsfrm.childs+i)->hwnd,
                        (wsfrm.childs+i)->hdcB);
              ReleaseDC((wsfrm.childs+i)->hwnd,
                        (wsfrm.childs+i)->hdcC);
            }
            ReleaseDC(wsfrm.hwnd,wsfrm.hdcB);
            ReleaseDC(wsfrm.hwnd,wsfrm.hdcC);
            DestroyWindow(wsfrm.hwnd);

            free(wsfrm.childs);
            wsfrm.nchilds=0;
            wsfrm.hwnd = NULL;

			(wmenu.childs+2)->vparam.vflag.mv.sectlist=0;
            SendMessage((wmenu.childs+2)->hwnd,WM_INITDIALOG,0,0);
          }
		  break;
        case IDM_PROPERTYLIST:
          if(wprop.hwnd==NULL &&
             wsfrm.hwnd!=NULL) /*CREATION.*/
          {
            getclientsize(wsfrm.hwnd,&maxX,&maxY);

            windowdefinition(&wprop,
                             hInstGlobal,
                             NULL,
                             WindowProcedureSheet,
                             NULL,
							 190,190,190,
                             WS_CLIPCHILDREN |
                             WS_CLIPSIBLINGS |
                             WS_CHILD |
                             /*WS_POPUP |*/
                             /*WS_CAPTION |*/
							 WS_THICKFRAME |
                             WS_VISIBLE,
							 0,0,250,150,
                             wsfrm.hwnd);

            wprop.nchilds=2;
            wprop.childs=(struct windowparams *)
                         malloc(2*sizeof(struct windowparams));

            getclientsize(wprop.hwnd,&maxX,&maxY);
            windowdefinition(wprop.childs+0,
                             hInstGlobal,
                             "CanPropBak",
                             WindowProcedureBack,
                             NULL,
							 190,190,190,
                             WS_CLIPCHILDREN |
                             WS_CLIPSIBLINGS |
                             WS_CHILD |
                             WS_VISIBLE,
                             1,1,maxX-BARWIDTH-5,maxY-BARWIDTH-5,
							 wprop.hwnd);

			windowdefinition(wprop.childs+1,
							 hInstGlobal,
							 "CanPropDsp",
							 WindowProcedureProp,
							 NULL,
							 190,190,190,
							 WS_CLIPCHILDREN |
							 WS_CLIPSIBLINGS |
							 WS_CHILD |
							 WS_VISIBLE,
							 0,0,500,500,
							 (wprop.childs+0)->hwnd);
			setfontformat((wprop.childs+1)->hdcC,15,6,
						  "MS Mincho",100,100,170);

            (wprop.childs+1)->vparam.Xo=0;
            (wprop.childs+1)->vparam.Yo=0;
          }
		  else if(wprop.hwnd!=NULL) /*DESTRUCTION.*/
		  {
            for(i=0;i<wprop.nchilds;i++)
            {
              ReleaseDC((wprop.childs+i)->hwnd,
                        (wprop.childs+i)->hdcB);
              ReleaseDC((wprop.childs+i)->hwnd,
                        (wprop.childs+i)->hdcC);
            }
			ReleaseDC(wprop.hwnd,wprop.hdcB);
            ReleaseDC(wprop.hwnd,wprop.hdcC);
            DestroyWindow(wprop.hwnd);

            free(wprop.childs);
            wprop.nchilds=0;
            wprop.hwnd = NULL;
            gnprop=0;
          }
          break;
        case IDM_SECTIONLIST:
          if(wsect.hwnd==NULL &&
             wsfrm.hwnd!=NULL) /*CREATION.*/
		  {
            getclientsize(wsfrm.hwnd,&maxX,&maxY);
            if(wprop.hwnd!=NULL)
			{
              getwindowsize(wprop.hwnd,&mw1,&mh1);
            }
            else
            {
              mw1=250;
              mh1=0;
            }

            windowdefinition(&wsect,
                             hInstGlobal,
                             NULL,
                             WindowProcedureSheet,
                             NULL,
							 190,190,190,
                             WS_CLIPCHILDREN |
                             WS_CLIPSIBLINGS |
                             WS_CHILD |
							 /*WS_POPUP |*/
                             /*WS_CAPTION |*/
                             WS_THICKFRAME |
                             WS_VISIBLE,
                             0,mh1,mw1,(maxY-mh1),
                             wsfrm.hwnd);

            wsect.nchilds=2;
            wsect.childs=(struct windowparams *)
                         malloc(2*sizeof(struct windowparams));

            getclientsize(wsect.hwnd,&maxX,&maxY);
            windowdefinition(wsect.childs+0,
                             hInstGlobal,
                             "CanSectBak",
                             WindowProcedureBack,
                             NULL,
							 190,190,190,
                             WS_CLIPCHILDREN |
                             WS_CLIPSIBLINGS |
                             WS_CHILD |
							 WS_VISIBLE,
                             1,1,maxX-BARWIDTH-5,maxY-BARWIDTH-5,
                             wsect.hwnd);

            windowdefinition(wsect.childs+1,
                             hInstGlobal,
                             "CanSectDsp",
                             WindowProcedureSect,
                             NULL,
                             190,190,190,
                             WS_CLIPCHILDREN |
                             WS_CLIPSIBLINGS |
                             WS_CHILD |
                             WS_VISIBLE,
                             0,0,500,1500,
                             (wsect.childs+0)->hwnd);
            setfontformat((wsect.childs+1)->hdcC,15,6,
						  "MS Mincho",100,100,170);

            (wsect.childs+1)->vparam.Xo=0;
            (wsect.childs+1)->vparam.Yo=0;
		  }
          else if(wsect.hwnd!=NULL) /*DESTRUCTION.*/
          {
            for(i=0;i<wsect.nchilds;i++)
            {
              ReleaseDC((wsect.childs+i)->hwnd,
                        (wsect.childs+i)->hdcB);
              ReleaseDC((wsect.childs+i)->hwnd,
                        (wsect.childs+i)->hdcC);
            }
            ReleaseDC(wsect.hwnd,wsect.hdcB);
            ReleaseDC(wsect.hwnd,wsect.hdcC);
            DestroyWindow(wsect.hwnd);

            free(wsect.childs);
            wsect.nchilds=0;
            wsect.hwnd = NULL;
            gnsect=0;
          }
          break;
        case IDM_SECTIONVIEW:
		  if(wsdsp.hwnd==NULL &&
             wsfrm.hwnd!=NULL) /*CREATION.*/
          {
			getclientsize(wsfrm.hwnd,&maxX,&maxY);

            if(wsect.hwnd!=NULL)
            {
              getwindowsize(wsect.hwnd,&mw1,&mh1);
            }
            else mw1=0;

            windowdefinition(&wsdsp,
                             hInstGlobal,
                             NULL,
                             WindowProcedureSheet,
                             "HOGSECTVIEW",
                             190,190,190,
                             WS_CLIPCHILDREN |
                             WS_CLIPSIBLINGS |
                             WS_CHILD |
                             /*WS_POPUP |*/
							 /*WS_CAPTION |*/
                             WS_THICKFRAME |
                             WS_VISIBLE,
                             mw1,0,maxX-mw1,maxY,
                             wsfrm.hwnd);

            wsdsp.nchilds=2;
            wsdsp.childs=(struct windowparams *)
                         malloc(2*sizeof(struct windowparams));

            getclientsize(wsdsp.hwnd,&maxX,&maxY);
            windowdefinition(wsdsp.childs+0,
                             hInstGlobal,
                             "CanSviewBak",
                             WindowProcedureBack,
                             NULL,
                             190,190,190,
                             WS_CLIPCHILDREN |
                             WS_CLIPSIBLINGS |
                             WS_CHILD |
                             WS_VISIBLE,
							 1,1,maxX-BARWIDTH-5,maxY-BARWIDTH-5,
                             wsdsp.hwnd);

			windowdefinition(wsdsp.childs+1,
                             hInstGlobal,
                             "CanSviewDsp",
                             WindowProcedureSview,
                             NULL,
                             0,0,0,
                             WS_CLIPCHILDREN |
                             WS_CLIPSIBLINGS |
                             WS_CHILD |
                             WS_VISIBLE,
                             0,0,500,500,
                             (wsdsp.childs+0)->hwnd);
            setfontformat((wsdsp.childs+1)->hdcC,15,6,
                          "MS Mincho",150,150,150);
            hpen=CreatePen(PS_SOLID,1,RGB(150,150,150));
            SelectObject((wsdsp.childs+1)->hdcC,hpen);

            getclientsize((wsdsp.childs+1)->hwnd,&maxX,&maxY);

            /*SECTION VIEWPOINT*/
            (wsdsp.childs+1)->vparam.type=AXONOMETRIC;
            (wsdsp.childs+1)->vparam.gfactor=1.0;
            (wsdsp.childs+1)->vparam.focus.d[0]=0.0;
            (wsdsp.childs+1)->vparam.focus.d[1]=0.0;
            (wsdsp.childs+1)->vparam.focus.d[2]=0.0;
            (wsdsp.childs+1)->vparam.theta=-90.0;
            (wsdsp.childs+1)->vparam.phi=90.0;
            (wsdsp.childs+1)->vparam.r=1000.0;
            (wsdsp.childs+1)->vparam.odv=1000.0;
            (wsdsp.childs+1)->vparam.Xo=(int)(maxX/2);
            (wsdsp.childs+1)->vparam.Yo=(int)(maxY/2);
            (wsdsp.childs+1)->vparam.dparam.gaxis=100.0;
            createviewdata(&((wsdsp.childs+1)->vparam));

            (wsdsp.childs+1)->lstatus=NEUTRAL;
          }
          else if(wsdsp.hwnd!=NULL) /*DESTRUCTION.*/
          {
            for(i=0;i<wsdsp.nchilds;i++)
			{
              ReleaseDC((wsdsp.childs+i)->hwnd,
                        (wsdsp.childs+i)->hdcB);
			  ReleaseDC((wsdsp.childs+i)->hwnd,
                        (wsdsp.childs+i)->hdcC);
            }
			ReleaseDC(wsdsp.hwnd,wsdsp.hdcB);
            ReleaseDC(wsdsp.hwnd,wsdsp.hdcC);
            DestroyWindow(wsdsp.hwnd);

            free(wsdsp.childs);
            wsdsp.nchilds=0;
            wsdsp.hwnd = NULL;
          }
          break;

		case IDM_OPTIONS:
          if(wmenu.hwnd==NULL)
          {
            getclientsize(hwnd,&maxX,&maxY);

			createviewdata(&vpdefault);

            windowdefinition(&wmenu,
                             hInstGlobal,
                             NULL,
                             WindowProcedureSheet,
                             NULL,
                             190,190,190,
                             WS_CLIPCHILDREN |
                             WS_CLIPSIBLINGS |
                             WS_CHILD |
                             WS_THICKFRAME |
                             /*WS_DLGFRAME |*/
                             WS_VISIBLE,
                             0,0,INPUTWIDTH,maxY,
                             wmain.hwnd);

            getclientsize(wmenu.hwnd,&maxX,&maxY);

            wmenu.nchilds=MENUDIALOGS+2;
            wmenu.childs=(struct windowparams *)
						 malloc(wmenu.nchilds*
                                sizeof(struct windowparams));

			windowdefinition(wmenu.childs+0,
                             hInstGlobal,
                             "CanMenuBak",
                             WindowProcedureBack,
                             NULL,
							 190,190,190,
                             WS_CLIPCHILDREN |
                             WS_CLIPSIBLINGS |
                             WS_CHILD |
                             WS_VISIBLE,
                             1,1,maxX-BARWIDTH-5,maxY-BARWIDTH-5,
                             wmenu.hwnd);

            windowdefinition(wmenu.childs+1,
                             hInstGlobal,
                             "CanMenuDsp",
                             WindowProcedureMenu,
                             NULL,
							 0,0,0,
                             WS_CLIPCHILDREN |
                             WS_CLIPSIBLINGS |
                             WS_CHILD |
                             WS_VISIBLE,
                             0,0,0,0,
                             (wmenu.childs+0)->hwnd);

			(wmenu.childs+2)->vparam.vflag.mv.draw=0;
            (wmenu.childs+2)->vparam.vflag.mv.error=0;
            (wmenu.childs+2)->vparam.vflag.mv.surface=0;
            (wmenu.childs+2)->vparam.vflag.mv.sectlist=0;
            (wmenu.childs+2)->vparam.vflag.mv.weight=0;
            (wmenu.childs+2)->vparam.vflag.mv.cmq=0;
            (wmenu.childs+2)->vparam.vflag.mv.horizon=0;

            (wmenu.childs+2)->vparam.vflag.mv.ftype
			=DEFAULTFILETYPE;
            (wmenu.childs+2)->vparam.vflag.mv.savetype
            =DEFAULTFILETYPE;

			(wmenu.childs+2)->hwnd
            = CreateDialog(hInstGlobal,
                           "CANOPTIONS1",
						   (wmenu.childs+1)->hwnd,
                           (DLGPROC)DialogProcMenu1);

            (wmenu.childs+3)->hwnd
            = CreateDialog(hInstGlobal,
                           "CANOPTIONS2",
                           (wmenu.childs+1)->hwnd,
                           (DLGPROC)DialogProcMenu2);

            (wmenu.childs+4)->hwnd
            = CreateDialog(hInstGlobal,
                           "CANOPTIONS3",
                           (wmenu.childs+1)->hwnd,
                           (DLGPROC)DialogProcMenu3);

            (wmenu.childs+5)->hwnd
            = CreateDialog(hInstGlobal,
                           "CANOPTIONS4",
						   (wmenu.childs+1)->hwnd,
                           (DLGPROC)DialogProcMenu4);

            getclientsize((wmenu.childs+2)->hwnd,&mw1,&mh1);
            getclientsize((wmenu.childs+3)->hwnd,&mw2,&mh2);
            getclientsize((wmenu.childs+4)->hwnd,&mw3,&mh3);
            getclientsize((wmenu.childs+5)->hwnd,&mw4,&mh4);
			MoveWindow((wmenu.childs+1)->hwnd,
					   0,0,mw1,(mh1+mh2+mh3+mh4),TRUE);
			MoveWindow((wmenu.childs+3)->hwnd,
					   0,mh1,mw2,mh2,TRUE);
			MoveWindow((wmenu.childs+4)->hwnd,
					   0,mh1+mh2,
					   mw3,mh3,TRUE);
			MoveWindow((wmenu.childs+5)->hwnd,
					   0,mh1+mh2+mh3,
					   mw4,mh4,TRUE);

			SetDlgItemText((wmenu.childs+3)->hwnd,IDV_MODENUM,"1");
          }
		  break;

		/*POPUP "ARCLM"*/
		case IDM_ARCLMRATE:                 /*SRCan Rate*/
		  if(wdraw.nchilds>=2 &&
             (wdraw.childs+1)->hwnd!=NULL &&
             (wmenu.childs+2)->vparam.vflag.mv.ftype==F_ARCLM)
          {
            fout=fgetstofopen("\0","r",ID_OUTPUTFILE);
            /*fout=fopen((wdraw.childs+1)->otpfile,"r");*/
            if(fout==NULL) break;

			readsrcanrate(fout,&arc);
            /*reasrcanrate(fout,&arc,NULL)*/
			fclose(fout);

            (wdraw.childs+1)->vparam.vflag.ev.srcancolor=1;

            clearwindow(*(wdraw.childs+1));
            drawarclmframe((wdraw.childs+1)->hdcC,
                           (wdraw.childs+1)->vparam,arc,0,ONSCREEN);
            SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);

            SendMessage((wmenu.childs+3)->hwnd,WM_INITDIALOG,0,0);

            MessageBox(NULL,"Completed.","Rate of Arclm",MB_OK);
          }
          break;
		case IDM_ARCLMLONG:               /*Long*/
		  if(MessageBox(NULL,"Arclm001 for Long.","ARCLM001",
			 MB_OKCANCEL)==IDCANCEL) break;

		  if((wmenu.childs+2)->vparam.vflag.mv.ftype==F_ARCLM)
		  {
			getviewparam((wmenu.childs+2)->hwnd,
						 &((wdraw.childs+1)->vparam));
			clearwindow(*(wdraw.childs+1));
			arclm001(&arc,ID_INPUTFILEZ,ID_OUTPUTFILEZ);
		  }
		  break;
		case IDM_ARCLMSHORTX:             /*Short X*/
          if(MessageBox(NULL,"Arclm001 for Short X.","ARCLM001",
             MB_OKCANCEL)==IDCANCEL) break;

		  if((wmenu.childs+2)->vparam.vflag.mv.ftype==F_ARCLM)
          {
			getviewparam((wmenu.childs+2)->hwnd,
                         &((wdraw.childs+1)->vparam));
            clearwindow(*(wdraw.childs+1));
            arclm001(&arcx,ID_INPUTFILEX,ID_OUTPUTFILEX);
          }
          break;
		case IDM_ARCLMSHORTY:             /*Short Y*/
          if(MessageBox(NULL,"Arclm001 for Short Y.","ARCLM001",
             MB_OKCANCEL)==IDCANCEL) break;

          if((wmenu.childs+2)->vparam.vflag.mv.ftype==F_ARCLM)
          {
            getviewparam((wmenu.childs+2)->hwnd,
                         &((wdraw.childs+1)->vparam));
            clearwindow(*(wdraw.childs+1));
            arclm001(&arcy,ID_INPUTFILEY,ID_OUTPUTFILEY);
          }
          break;

        /* 150514 fukushima for all */
        case IDM_ARCLMALL:                  /*All*/
          if(MessageBox(NULL,"Arclm001 Long,Short X,Short Y.","ARCLM001",
             MB_OKCANCEL)==IDCANCEL) break;

          if((wmenu.childs+2)->vparam.vflag.mv.ftype==F_ARCLM)
          {
            getviewparam((wmenu.childs+2)->hwnd,
                         &((wdraw.childs+1)->vparam));

            /*
            clearwindow(*(wdraw.childs+1));
            arclm001(&arc,ID_INPUTFILEZ,ID_OUTPUTFILEZ);
            clearwindow(*(wdraw.childs+1));
            arclm001(&arcx,ID_INPUTFILEX,ID_OUTPUTFILEX);
            clearwindow(*(wdraw.childs+1));
            arclm001(&arcy,ID_INPUTFILEY,ID_OUTPUTFILEY);
            */

            clearwindow(*(wdraw.childs+1));
			int idinputs[3],idoutputs[3];
            struct arclmframe *arcs[3];
            arcs[0]=&arc;arcs[1]=&arcx;arcs[2]=&arcy;
			// arcs[0]=&arc;arcs[1]=&arc;arcs[2]=&arc;
            idinputs[0]=ID_INPUTFILEZ;idinputs[1]=ID_INPUTFILEX;idinputs[2]=ID_INPUTFILEY;
            idoutputs[0]=ID_OUTPUTFILEZ;idoutputs[1]=ID_OUTPUTFILEX;idoutputs[2]=ID_OUTPUTFILEY;
			arclm001_lxy(arcs,idinputs,idoutputs);
          }
          break;
        /***/

		case IDM_ARCLMTEST:                   /*Test*/
          if(MessageBox(NULL,"Arclm001 test.","ARCLM001",
             MB_OKCANCEL)==IDCANCEL) break;

          if((wmenu.childs+2)->vparam.vflag.mv.ftype==F_ARCLM)
          {
			getviewparam((wmenu.childs+2)->hwnd,
                         &((wdraw.childs+1)->vparam));
			clearwindow(*(wdraw.childs+1));
            arclmtest(&arc,ID_INPUTFILE,ID_OUTPUTFILE);
          }
          break;
		case IDM_ARCLMFIGURES:                 /*Figures*/
          if((wmenu.childs+2)->vparam.vflag.mv.ftype==F_ARCLM)
          {
			fout=fgetstofopen(DIRECTORY,"r",ID_SECTIONFILE);
            if(fout==NULL) break;

			(wdraw.childs+1)->vparam.vflag.mv.pagenum=1;
            printarclmfigures2(fout,&((wdraw.childs+1)->vparam),1,0);
            fclose(fout);

            (wdraw.childs+1)->vparam.vflag.mv.pagenum=0;
            clearwindow(*(wdraw.childs+1));
            drawarclmframe((wdraw.childs+1)->hdcC,
                           (wdraw.childs+1)->vparam,arc,0,ONSCREEN);
            SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);

            SendMessage((wmenu.childs+3)->hwnd,WM_INITDIALOG,0,0);

            MessageBox(NULL,"Completed.","Figures of Arclm",MB_OK);
          }
		  break;

		case IDM_SRCANFIGURES:                 /*Figures*/
		  if((wmenu.childs+2)->vparam.vflag.mv.ftype==F_ARCLM)
		  {
			fout=fgetstofopen(DIRECTORY,"r",ID_SECTIONFILE);
			if(fout==NULL) break;

			(wdraw.childs+1)->vparam.vflag.mv.pagenum=1;
			printsrcanfigures(fout,&((wdraw.childs+1)->vparam),3,0);
			fclose(fout);

			(wdraw.childs+1)->vparam.vflag.mv.pagenum=0;
			clearwindow(*(wdraw.childs+1));
			drawarclmframe((wdraw.childs+1)->hdcC,
						   (wdraw.childs+1)->vparam,arc,0,ONSCREEN);
			SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);

			SendMessage((wmenu.childs+3)->hwnd,WM_INITDIALOG,0,0);

			MessageBox(NULL,"Completed.","Figures of Arclm",MB_OK);
		  }
		  break;

		case IDM_ARCLMLASTFIGURES:             /*Zobun Figures*/
          if((wmenu.childs+2)->vparam.vflag.mv.ftype==F_ARCLM)
          {
			fout=fgetstofopen(DIRECTORY,"r",ID_SECTIONFILE);
            if(fout==NULL) break;
            globalmessageflag=0;
			printarclmlastfiguresII(fout,&((wdraw.childs+1)->vparam),&arc);
			fclose(fout);

			clearwindow(*(wdraw.childs+1));
			drawarclmframe((wdraw.childs+1)->hdcC,
						   (wdraw.childs+1)->vparam,arc,0,ONSCREEN);
			SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);

			SendMessage((wmenu.childs+3)->hwnd,WM_INITDIALOG,0,0);

			MessageBox(NULL,"Completed.","Figures of Arclm",MB_OK);
		  }
          break;

		case IDM_CMQCHECK:                      /*CMQ CHECK*/ //kaza & ujok
		  if(MessageBox(NULL,"CMQ CHECK.","CMQ CHECK",
			 MB_OKCANCEL)==IDCANCEL) break;

		  if((wmenu.childs+2)->vparam.vflag.mv.ftype==F_ARCLM)
		  {
			 flagswitch(&((wdraw.childs+1)->vparam.vflag.ev.cmqcheck));
		  }
		  break;

        case IDM_OPENENERGYDISTRIBUTION:           /*Energy Distribution*/
		  if(MessageBox(NULL,"OPEN FILE.","ENERGY DISTRIBUTION",
			 MB_OKCANCEL)==IDCANCEL) break;

		  if((wmenu.childs+2)->vparam.vflag.mv.ftype==F_ARCLM)
		  {
            if(arc.nodes==NULL) break;

			fout=fgetstofopen("\0","r",ID_OUTPUTFILE);  /*OPEN FILE.*/
            if(fout==NULL) break;

            energyoutputtomemory(fout,&arc);      //ujok
            fclose(fout);

			arc.nlaps=1;

			clearwindow(*(wdraw.childs+1));
			drawarclmframe((wdraw.childs+1)->hdcC,
						   (wdraw.childs+1)->vparam,arc,0,ONSCREEN);
			SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);

			SetDlgItemText((wmenu.childs+3)->hwnd,IDV_MODENUM,"1");

			(wdraw.childs+1)->lstatus=ROTATE;

            (wdraw.childs+1)->vparam.vflag.ev.deformation=1;   //ujok
		  }
		  break;


		//case IDM_BANDDECREASEDARCLM:            /*Save Band Decreased Arclm*/
		  if((wmenu.childs+2)->vparam.vflag.mv.ftype==F_ARCLM)
		  {
            savebanddecreasedarclm(&arc);
            MessageBox(NULL,"Completed.","Save Band Decreased Arclm",MB_OK);
          }
          break;

		/*POPUP "Arclm"*/

        case IDM_BANDDECREASEDORGAN: /*CODED BY ARAKI 2008.09.12*/
          if((wmenu.childs+2)->vparam.vflag.mv.ftype==F_ARCLM)
          {
            fout=fopen("bda.inp","w");
            if(fout==NULL) break;
            savebanddecreasedorgan(fout,&((wdraw.childs+1)->org),
                                        &((wdraw.childs+1)->vparam),&arcx);

            fclose(fout);
            sprintf(str,"Saved As \"bda.inp\"");
			MessageBox(NULL,str,"Save Band Decreased Organ Frame",MB_OK);
          }
          break;

		case IDM_SINGULATENODE:
		  if(wdraw.nchilds>=2 &&
			 (wdraw.childs+1)->hwnd!=NULL &&
			 (wmenu.childs+2)->vparam.vflag.mv.ftype==F_ORGAN &&
			 (wdraw.childs+1)->org.nodes!=NULL)
		  {
			ncount=singulatenode(&((wdraw.childs+1)->org));

			clearwindow(*(wdraw.childs+1));
			draworganization((wdraw.childs+1)->hdcC,
							 (wdraw.childs+1)->vparam,
							 (wdraw.childs+1)->org,ONSCREEN);
			SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);

			sprintf(str,"Deleted Nodes = %ld",ncount);
			MessageBox(NULL,str,"Singulate Node",MB_OK);
		  }
		  break;

		case IDM_SAVEORGAN:                         /*Save*/
		  if(wdraw.nchilds>=2 &&
             (wdraw.childs+1)->hwnd!=NULL &&
             (wmenu.childs+2)->vparam.vflag.mv.ftype==F_ORGAN &&
             (wdraw.childs+1)->org.nodes!=NULL)
          {
            /*fout=fopen((wdraw.childs+1)->otpfile,"w");*/
            fout=fopen("hogtxt.inp","w");
            if(fout==NULL) break;

            saveorganization(fout,&((wdraw.childs+1)->org),
                                  &((wdraw.childs+1)->vparam));

            fclose(fout);

            /*sprintf(str,"Saved As \"%s\"",
                    (wdraw.childs+1)->otpfile);*/
			sprintf(str,"Saved As \"hogtxt.inp\"");
			MessageBox(NULL,str,"Save Organ Frame",MB_OK);
          }
          break;


		case IDM_SAVEORGANAS:                   /*Save As (INPUTFILE)*/
          if(wdraw.nchilds>=2 &&
             (wdraw.childs+1)->hwnd!=NULL &&
             (wmenu.childs+2)->vparam.vflag.mv.ftype==F_ORGAN &&
             (wdraw.childs+1)->org.nodes!=NULL)
          {
            fout=fopen((wdraw.childs+1)->inpfile,"w");
            if(fout==NULL) break;

            saveorganization(fout,&((wdraw.childs+1)->org),
                                  &((wdraw.childs+1)->vparam));

            fclose(fout);

            sprintf(str,"Saved As \"%s\"",
                    (wdraw.childs+1)->inpfile);
			MessageBox(NULL,str,"Save Organ Frame",MB_OK);
          }
          break;

		case IDM_SAVEORGANASDXF:                  /*Save Organ as Dxf*/
          if(wdraw.nchilds>=2 &&
             (wdraw.childs+1)->hwnd!=NULL &&
             (wmenu.childs+2)->vparam.vflag.mv.ftype==F_ORGAN &&
             (wdraw.childs+1)->org.nodes!=NULL)
          {
            fout=fopen("hogdxf.txt","w");
            if(fout==NULL) break;

            saveorganasdxf(fout,(wdraw.childs+1)->org,
                           (wdraw.childs+1)->vparam);

            fclose(fout);

            sprintf(str,"Saved As \"hogdxf.txt\"");
			MessageBox(NULL,str,"Save Organ as Dxf",MB_OK);
		  }
          break;

		case IDM_SAVEARCLMASDXF:                  /*Save Arclm as Dxf*/
          if(wdraw.nchilds>=2 &&
             (wdraw.childs+1)->hwnd!=NULL &&
             (wmenu.childs+2)->vparam.vflag.mv.ftype==F_ARCLM)
          {
            fout=fopen("hogdxf.txt","w");
            if(fout==NULL) break;

            savearclmasdxf(fout,(wdraw.childs+1)->vparam,arc);

            fclose(fout);

			sprintf(str,"Saved As \"hogdxf.txt\"");
			MessageBox(NULL,str,"Save Arclm as Dxf",MB_OK);
          }
          break;

		case IDM_DISTRIBUTEORGAN:                /*Weight Distribution*/
          if(wdraw.nchilds>=2 &&
			 (wdraw.childs+1)->hwnd!=NULL &&
             (wmenu.childs+2)->vparam.vflag.mv.ftype==F_ORGAN &&
             (wdraw.childs+1)->org.nodes!=NULL)
          {
			(wdraw.childs+1)->vparam.vflag.nv.mcircle=1;
            /*(wdraw.childs+1)->vparam.vflag.nv.mvalue=1;*/
            SendMessage((wmenu.childs+3)->hwnd,WM_INITDIALOG,0,0);

            weightdistribution(NULL,NULL,
                               &((wdraw.childs+1)->vparam),
                               &((wdraw.childs+1)->org));

            clearwindow(*(wdraw.childs+1));
            draworganization((wdraw.childs+1)->hdcC,
                             (wdraw.childs+1)->vparam,
                             (wdraw.childs+1)->org,ONSCREEN);
            SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);
          }
          break;

		case IDM_EXTRACTARCLM:                   /*Extract Arclm*/
		  if(wdraw.nchilds>=2 &&
             (wdraw.childs+1)->hwnd!=NULL &&
			 (wmenu.childs+2)->vparam.vflag.mv.ftype==F_ORGAN &&
			 (wdraw.childs+1)->org.nodes!=NULL)
		  {
			if(extractarclmfromorgan(&((wdraw.childs+1)->org),
									 &arc,&arcx,&arcy)==0)
			{
			  break;
			}

			(wmenu.childs+2)->vparam.vflag.mv.ftype=F_ARCLM;
			SendMessage((wmenu.childs+2)->hwnd,WM_INITDIALOG,0,0);

			clearwindow(*(wdraw.childs+1));
			drawarclmframe((wdraw.childs+1)->hdcC,
						   (wdraw.childs+1)->vparam,arc,0,ONSCREEN);
			SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);

			MessageBox(NULL,"Completed.","Organ Into Arclm",MB_OK);
		  }
		  break;
		case IDM_ORGANFIGURES:                 /*Figures*/
		  if(wdraw.nchilds>=2 &&
             (wdraw.childs+1)->hwnd!=NULL &&
			 (wmenu.childs+2)->vparam.vflag.mv.ftype==F_ORGAN &&
			 (wdraw.childs+1)->org.nodes!=NULL)
          {
            fout=fgetstofopen(DIRECTORY,"r",ID_SECTIONFILE);
            if(fout==NULL) break;

			(wdraw.childs+1)->vparam.vflag.mv.pagenum=1;
            printorganfigures(fout,&((wdraw.childs+1)->vparam),0);

            fclose(fout);

            (wdraw.childs+1)->vparam.vflag.mv.pagenum=0;
            clearwindow(*(wdraw.childs+1));
			draworganization((wdraw.childs+1)->hdcC,
							 (wdraw.childs+1)->vparam,
							 (wdraw.childs+1)->org,ONSCREEN);
            SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);

            SendMessage((wmenu.childs+3)->hwnd,WM_INITDIALOG,0,0);

            MessageBox(NULL,"Completed.","Figures of Organ",MB_OK);
          }
          break;

		case IDM_KESANSHO:                 /*UJIOKA FOR KESANSHO*/
		  if(wdraw.nchilds>=2 &&
             (wdraw.childs+1)->hwnd!=NULL &&
			 (wmenu.childs+2)->vparam.vflag.mv.ftype==F_ORGAN &&
			 (wdraw.childs+1)->org.nodes!=NULL)
          {
            fout=fgetstofopen(DIRECTORY,"r",ID_SECTIONFILE);
            if(fout==NULL) break;
			(wdraw.childs+1)->vparam.vflag.mv.pagenum=1;

            printorganfigures(fout,&((wdraw.childs+1)->vparam),1);     /*2.2*/

            printsrcanmax(fout,&((wdraw.childs+1)->vparam),1);   /*4.4*/
            printsrcanfigures(fout,&((wdraw.childs+1)->vparam),3,1);   /*4.4*/

            printarclmfigures2(fout,&((wdraw.childs+1)->vparam),1,0);  /*7.2-4*/
            fclose(fout);

			(wdraw.childs+1)->vparam.vflag.mv.pagenum=0;
            (wmenu.childs+2)->vparam.vflag.mv.ftype==F_ORGAN;
            clearwindow(*(wdraw.childs+1));
			draworganization((wdraw.childs+1)->hdcC,
							 (wdraw.childs+1)->vparam,
							 (wdraw.childs+1)->org,ONSCREEN);
            SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);

            SendMessage((wmenu.childs+3)->hwnd,WM_INITDIALOG,0,0);

            MessageBox(NULL,"Completed.","Kesansho",MB_OK);
          }
          break;
        /*Optimization Test*/
		case IDM_OPTIMIZATION: /*CONJUGATE*/
		  if(wdraw.nchilds>=2 &&
			 (wdraw.childs+1)->hwnd!=NULL &&
			 (wmenu.childs+2)->vparam.vflag.mv.ftype==F_ORGAN &&
			 (wdraw.childs+1)->org.nodes!=NULL)
		  {
			conjugategradientbezier(&((wdraw.childs+1)->org));
			/*conjugategradientxz(&((wdraw.childs+1)->org));*/
			/*conjugategradientcurve(&((wdraw.childs+1)->org));*/
		  }
		  break;

        /*change sect auto*/   /*For ghouse*/
        case IDM_CHANGESECTAUTO:
          if(wdraw.nchilds>=2 &&
			 (wdraw.childs+1)->hwnd!=NULL &&
			 (wmenu.childs+2)->vparam.vflag.mv.ftype==F_ORGAN &&
			 (wdraw.childs+1)->org.nodes!=NULL)
		  {
			changesectautomatic(&((wdraw.childs+1)->org),0);
		  }
		  break;
        /*change sect auto*/

		/*Create Bezier Surface*/
		case IDM_CREATEBEZIERSURFACE:
          if(wdraw.nchilds>=2 &&
			 (wdraw.childs+1)->hwnd!=NULL &&
			 (wmenu.childs+2)->vparam.vflag.mv.ftype==F_ORGAN &&
			 (wdraw.childs+1)->org.nodes!=NULL)
		  {
		    sprintf(str,"TEST: Create Bezier Surface from Organ");
//		    sprintf(str,"TEST: Create Bezier Vierendeel Arch from Organ");
       	    if(MessageBox(NULL,str,"Create Bezier Surface",MB_OKCANCEL)==IDOK)
            {
               createbeziersurfacetest(&((wdraw.childs+1)->org));
//               createvierendeelarch(&((wdraw.childs+1)->org));
		    }
          }
          break;

		/*
		case IDM_TESTCALCULATION:
          if(wdraw.nchilds>=2 &&
			 (wdraw.childs+1)->hwnd!=NULL &&
			 (wmenu.childs+2)->vparam.vflag.mv.ftype==F_ORGAN &&
			 (wdraw.childs+1)->org.nodes!=NULL)
		  {
		    sprintf(str,"TEST: Update Form and Section");
       	    if(MessageBox(NULL,str,"UPDATE ORGANIZATION",MB_OKCANCEL)==IDOK)
            {
			  updateorganization(&((wdraw.childs+1)->org));
		    }
          }
		  break;
		*/

		/*BCLNG*/
		case IDM_OPENBCLNGRESULT:
		  if(wdraw.nchilds>=2 &&
			 (wdraw.childs+1)->hwnd!=NULL &&
             (wmenu.childs+2)->vparam.vflag.mv.ftype==F_ARCLM)
          {
            if(arc.nodes==NULL) break;

			fout=fgetstofopen("\0","r",ID_OUTPUTFILE);  /*OPEN FILE.*/
            if(fout==NULL) break;

            bclngoutputtomemoryII(fout,&arc);      //ujok
            fclose(fout);

			arc.nlaps=1;

			clearwindow(*(wdraw.childs+1));
			drawarclmframe((wdraw.childs+1)->hdcC,
						   (wdraw.childs+1)->vparam,arc,0,ONSCREEN);
			SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);

			SetDlgItemText((wmenu.childs+3)->hwnd,IDV_MODENUM,"1");

			(wdraw.childs+1)->lstatus=ROTATE;

            (wdraw.childs+1)->vparam.vflag.ev.deformation=1;   //ujok
		  }
		  break;

		case IDM_OPENBCLNG003RESULT:
		  if(wdraw.nchilds>=2 &&
			 (wdraw.childs+1)->hwnd!=NULL &&
             (wmenu.childs+2)->vparam.vflag.mv.ftype==F_ARCLM)
          {
            if(arc.nodes==NULL) break;

			fout=fgetstofopen("\0","r",ID_OUTPUTFILE);  /*OPEN FILE.*/
            if(fout==NULL) break;

            bclngoutputtomemoryIII(fout,&arc,&multiwire,&nmultiwire);      //ujok
            fclose(fout);
/*
            sprintf(str,"NMULTIWIRE =%d",nmultiwire);
			errormessage(str);

            for(i=0;i<nmultiwire;i++)
            {
              {
                sprintf(str,"MULTIWIRE[%ld]:CODE=%ld",i,(*(multiwire+i))->code);
                errormessage(str);
              }
            }
*/
			arc.nlaps=1;

			globalcondensationflag=1;

			clearwindow(*(wdraw.childs+1));
			drawarclmframe((wdraw.childs+1)->hdcC,
						   (wdraw.childs+1)->vparam,arc,0,ONSCREEN);
			SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);

			SetDlgItemText((wmenu.childs+3)->hwnd,IDV_MODENUM,"1");

			(wdraw.childs+1)->lstatus=ROTATE;

            (wdraw.childs+1)->vparam.vflag.ev.deformation=1;   //ujok
		  }
		  break;

		case IDM_BCLNGANALYSIS:              //ujioka
		  if((wmenu.childs+2)->vparam.vflag.mv.ftype==F_ARCLM)
		  {
			getviewparam((wmenu.childs+2)->hwnd,
						 &((wdraw.childs+1)->vparam));
			clearwindow(*(wdraw.childs+1));
			arclm001(&arc,ID_INPUTFILEZ,ID_OUTPUTFILEZ);

			if(MessageBox(NULL,"Bclng001:Buckling Analysis Begin.","BCLNG001",
			   MB_OKCANCEL)==IDOK) bclng001(&arc);
		  }
		  break;

        case IDM_BCLNGANALYSISAUTO:              //ujioka
		  if(wdraw.nchilds>=2 &&
			 (wdraw.childs+1)->hwnd!=NULL &&
			 (wmenu.childs+2)->vparam.vflag.mv.ftype==F_ORGAN &&
			 (wdraw.childs+1)->org.nodes!=NULL)
		  {
			double beig = bclngautomatic(&((wdraw.childs+1)->org));
            sprintf(str,"BUCKLING EIGEN VALUE =%8.3f",beig);
		    MessageBox(NULL,str,"BCLNG001 RESULT",MB_OK);
		  }
		  break;

		case IDM_BCLNGCONDENSATION:           //BCLNG003
		  if((wmenu.childs+2)->vparam.vflag.mv.ftype==F_ARCLM)
		  {
            if(nmultiwire)
            {
    		  getviewparam((wmenu.childs+2)->hwnd,
                          &((wdraw.childs+1)->vparam));
              clearwindow(*(wdraw.childs+1));

              //arclm001(&arc,ID_INPUTFILEZ,ID_OUTPUTFILEZ);
              //線形解析を回すとmultiwireのELEM CODEがおかしくなる

              fout=fgetstofopen("\0","r",ID_OUTPUTFILEZ);  /*OPEN FILE.*/
              if(fout==NULL) break;

              frameoutputtomemory(fout,&arc);

              fclose(fout);

              /*CHECK*/
              /*
              for(i=0;i<nmultiwire;i++)
              {
                sprintf(str,"MULTIWIRE[%ld]:CODE=%ld",i+1,(*(multiwire+i))->code);
                errormessage(str);
              }
              */

              if(MessageBox(NULL,"Bclng003:Buckling Condensation Begin.","BCLNG003",
                            MB_OKCANCEL)==IDOK) bclng003(&arc,multiwire,nmultiwire);
              globalcondensationflag=1;

              globalstatus=NEUTRAL;
            }
            else
            {
    		  getviewparam((wmenu.childs+2)->hwnd,
			    			 &((wdraw.childs+1)->vparam));
              clearwindow(*(wdraw.childs+1));

              if(MessageBox(NULL,"Bclng003:Buckling Condensation(From Text Mode) Begin.","BCLNG003",
                           MB_OKCANCEL)==IDOK)
              {
                offset=0;
				while(1)
                {
                  multiwire=selectmultiwireII(&arc,offset,&nmultiwire);   /*FROM TEXT MODE*/
                  if(!nmultiwire) break;

                  for(i=0;i<nmultiwire;i++)
                  {
                    sprintf(str,"MULTIWIRE[%ld]:CODE=%ld",i+1,(*(multiwire+i))->code);
                    errormessage(str);
                  }

                  fout=fgetstofopen("\0","r",ID_OUTPUTFILEZ);  /*OPEN FILE.*/
                  if(fout==NULL) break;
                  frameoutputtomemory(fout,&arc);
                  fclose(fout);

                  bclng003(&arc,multiwire,nmultiwire);

                  nmultiwire=0;
                  free(multiwire);

                  offset++;
                }

                /*Extract Buckling Safety Ratio :.brat file*/
                GetDlgItemText((wmenu.childs+2)->hwnd,ID_OUTPUTFILE,s,80);
                n=strcspn(s,".");
                strncpy(str,s,n);
                str[n]='\0';
                strcat(str,".brat");
                frat=fopen(str,"w");
                for(i=0;i<arc.nelem;i++)
                {
                  fprintf(frat,"ELEM: %5d SECT: %4d %12.5f 0.00000 0.00000 0.00000\n",
                         (arc.elems+i)->code,(arc.elems+i)->sect->code,(arc.elems+i)->srate[0]);
                }
              }
            }
            fclose(frat);
            (wdraw.childs+1)->lstatus=ROTATE;
            (wdraw.childs+1)->vparam.vflag.ev.deformation=0;
            (wdraw.childs+1)->vparam.vflag.ev.srcanrate=1;
			(wdraw.childs+1)->vparam.vflag.ev.srcancolor=1;
		  }
		  break;

        case IDM_BCLNGOPTIMIZATION: /*CONJUGATE*/
		  if(wdraw.nchilds>=2 &&
			 (wdraw.childs+1)->hwnd!=NULL &&
			 (wmenu.childs+2)->vparam.vflag.mv.ftype==F_ORGAN &&
			 (wdraw.childs+1)->org.nodes!=NULL)
		  {
			bclngconjugategradient(&((wdraw.childs+1)->org));
		  }
		  break;


		case IDM_QUADOPEN:
		  wparam = MAKEWPARAM((WORD)IDD_VIEW,(WORD)0);
		  SendMessage((wmenu.childs+2)->hwnd,WM_COMMAND,
					  wparam,(WPARAM)0);
		  break;

		case IDM_QUADTEST:
		  if(bqf.nnode>0)
		  {
			arclm301(&bqf);
		  }
		  break;

		case IDM_CADREINPUT:
          if(wdraw.nchilds>=2 &&
             (wdraw.childs+1)->hwnd!=NULL &&
             (wmenu.childs+2)->vparam.vflag.mv.ftype==F_ARCLM)
          {
            fin=fgetstofopen("c:\\e\\yasuda\\kyotoy\\data\\",
                             "r",ID_INPUTFILE);         /*OPEN FILE.*/
            if(fin==NULL) break;

            inputcadretomemory(fin,&arc);
            fclose(fin);

            getviewparam((wmenu.childs+2)->hwnd,
                         &((wdraw.childs+1)->vparam));

            clearwindow(*(wdraw.childs+1));
            drawarclmframe((wdraw.childs+1)->hdcC,
                           (wdraw.childs+1)->vparam,arc,0,ONSCREEN);
            SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);

            (wdraw.childs+1)->lstatus=ROTATE;
          }
          break;
        case IDM_CADRERESULT:
          if(wdraw.nchilds>=2 &&
             (wdraw.childs+1)->hwnd!=NULL &&
             (wmenu.childs+2)->vparam.vflag.mv.ftype==F_ARCLM)
          {
            if(arc.nodes==NULL) break;

            fin =fgetstofopen("c:\\e\\yasuda\\kyotoy\\data\\",
                              "r",ID_INPUTFILE);       /*INPUT FILE.*/
			fout=fgetstofopen("c:\\e\\yasuda\\kyotoy\\data\\",
                              "r",ID_OUTPUTFILE);     /*RESULT FILE.*/
            if(fin==NULL || fout==NULL) break;

            cadreoutputtomemory(fout,fin,&arc);
            fclose(fin);
            fclose(fout);

            arc.nlaps=1;

            clearwindow(*(wdraw.childs+1));
            drawarclmframe((wdraw.childs+1)->hdcC,
                           (wdraw.childs+1)->vparam,arc,0,ONSCREEN);
            SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);

            SetDlgItemText((wmenu.childs+3)->hwnd,IDV_MODENUM,"1");

            (wdraw.childs+1)->lstatus=ROTATE;
          }
          break;

		case IDM_PRINTSETUP: /*PRINTER SETUP.*/
          if(gprn.pflag==0)
          {
			GetDlgItemText((wmenu.childs+3)->hwnd,
                           IDP_JIHEIGHT,str,20);
            gprn.jiheight=(int)strtol(str,NULL,10);
            GetDlgItemText((wmenu.childs+3)->hwnd,
                           IDP_JIWIDTH,str,20);
            gprn.jiwidth=(int)strtol(str,NULL,10);

            gprn.pd.lStructSize = sizeof(PRINTDLG);
            gprn.pd.hwndOwner = (HWND)NULL;
            gprn.pd.hDevMode = (HANDLE)NULL;
            gprn.pd.hDevNames = (HANDLE)NULL;
            gprn.pd.hDC = NULL;
            gprn.pd.Flags = PD_RETURNDC;
            gprn.pd.nFromPage = (WORD)NULL;
            gprn.pd.nToPage = (WORD)NULL;
            gprn.pd.nMinPage = (WORD)NULL;
            gprn.pd.nMaxPage = (WORD)NULL;
            gprn.pd.nCopies = (WORD)NULL;
			gprn.pd.hInstance = /*(HANDLE)*/(HINSTANCE)NULL;
			gprn.pd.lCustData = 0L;
            gprn.pd.lpfnPrintHook = (LPPRINTHOOKPROC)NULL;
            gprn.pd.lpfnSetupHook = (LPSETUPHOOKPROC)NULL;
			gprn.pd.lpPrintTemplateName = (LPSTR)NULL;
			gprn.pd.lpSetupTemplateName = (LPSTR)NULL;
			gprn.pd.hPrintTemplate = (HANDLE)NULL;
            gprn.pd.hSetupTemplate = (HANDLE)NULL;

            if(PrintDlg(&(gprn.pd)) != FALSE)
			{
              if(!(GetDeviceCaps(gprn.pd.hDC, RASTERCAPS)
				 & RC_BITBLT))
              {
                MessageBox(NULL,"Bitmap Not Available","Print",
                           MB_OK);
              }

              gprn.di.cbSize = sizeof(DOCINFO);
              sprintf(doc,"CanvsPrint");
              gprn.di.lpszDocName = doc;
			  gprn.di.lpszOutput = (LPTSTR)NULL;
              gprn.di.lpszDatatype = (LPTSTR) NULL;
              gprn.di.fwType = 0;

              StartDoc(gprn.pd.hDC,&(gprn.di));

              gprn.cWidthPels  = GetDeviceCaps(gprn.pd.hDC,HORZRES);
              gprn.cHeightPels = GetDeviceCaps(gprn.pd.hDC,VERTRES);

              gprn.cWidthDpi  = GetDeviceCaps(gprn.pd.hDC,
                                              LOGPIXELSX);
              gprn.cHeightDpi = GetDeviceCaps(gprn.pd.hDC,
                                              LOGPIXELSY);

              gprn.pflag=1;
            }
          }
          break;
        case IDM_PRINTWAIT: /*PRINT FRAME AND WAIT.*/
		  if(gprn.pflag==1 && gprn.pageflag==0 &&
             wdraw.hwnd!=NULL)
		  {
            StartPage(gprn.pd.hDC);
            gprn.pageflag=1;
          }

          if(gprn.pflag==1 && wdraw.hwnd!=NULL)
          {
            pfactor=10.0;
            vprint=(wdraw.childs+1)->vparam;
            vprint.Xo=0.5*(gprn.cWidthPels);
            vprint.Yo=0.5*(gprn.cHeightPels);
            vprint.gfactor*=pfactor;

            setfontformat(gprn.pd.hDC,gprn.jiheight,gprn.jiwidth,
                          "MS Mincho",0,0,0);
            SetBkMode(gprn.pd.hDC,TRANSPARENT);

            if(((wmenu.childs+2)->vparam.vflag.mv.ftype
                ==F_ARCLM ||
                (wmenu.childs+2)->vparam.vflag.mv.ftype
                ==F_FRAME) &&
               arc.elems!=NULL)
            {
              drawarclmframe(gprn.pd.hDC,vprint,arc,0,ONPRINTER);
			}
            if((wmenu.childs+2)->vparam.vflag.mv.ftype==F_ORGAN &&
               (wdraw.childs+1)->org.elems!=NULL)
            {
              draworganization(gprn.pd.hDC,vprint,
                               (wdraw.childs+1)->org,ONPRINTER);
            }

            setfontformat(gprn.pd.hDC,
                          (int)(1.5*gprn.jiheight),
                          (int)(1.5*gprn.jiwidth),
                          "MS Mincho",0,0,0);

            for(i=0;i<wdraw.nstring;i++)
            {
              (wdraw.strset+i)->n.d[0]*=pfactor;
              (wdraw.strset+i)->n.d[1]*=pfactor;
            }
            drawtexts(gprn.pd.hDC,wdraw.strset,wdraw.nstring,vprint);
            for(i=0;i<wdraw.nstring;i++)
            {
              (wdraw.strset+i)->n.d[0]/=pfactor;
              (wdraw.strset+i)->n.d[1]/=pfactor;
            }

			hfont=(HFONT)GetCurrentObject(gprn.pd.hDC,OBJ_FONT);
			DeleteObject(hfont);

            /*MessageBox(NULL,"Waiting.","Print",MB_OK);*/
          }
          break;
        case IDM_PRINTOUT: /*OUTPUT PAPER.*/
          if(gprn.pflag==1 && gprn.pageflag==1 &&
             wdraw.hwnd!=NULL)
          {
            EndPage(gprn.pd.hDC);
            gprn.pageflag=0;
            /*MessageBox(NULL,"Page Completed.","Print",MB_OK);*/
          }
          break;
        case IDM_PRINTEND: /*END PRINT DOCUMENT.*/
          if(gprn.pflag==1 && wdraw.hwnd!=NULL)
		  {
            EndDoc(gprn.pd.hDC);
            DeleteDC(gprn.pd.hDC);

            if(gprn.pd.hDevMode != NULL)
			{GlobalFree(gprn.pd.hDevMode);}
			if(gprn.pd.hDevNames != NULL)
			{GlobalFree(gprn.pd.hDevNames);}

			gprn.pflag=0;
			MessageBox(NULL,"Completed.","Print",MB_OK);
		  }
		  break;

		case IDM_PREVIEW: /*PRINT PREVIEW BY BITMAP.*/
		  if(wdraw.hwnd!=NULL)
		  {
			prn.jiheight= 30; /*50 FOR 600dpi*/
							  /*30 FOR DEFAULT, SD REVIEW*/
			prn.jiwidth = 12; /*20 FOR 600dpi*/
							  /*12 FOR DEFAULT, SD REVIEW*/

			gprn=prn; /*COPY TO GLOBAL*/

			cWidthPels  = 2800; /*4800 FOR A4 YOKO 600dpi*/
								/*2800 FOR SD REVIEW*/
			cHeightPels = 2800; /*6372 FOR A4 TATE 600dpi*/
								/*2800 FOR SD REVIEW*/

			pfactor= 5.833; /*10.0 FOR 600dpi*/
			vprint=(wdraw.childs+1)->vparam;
			vprint.Xo=0.5*cWidthPels;
			vprint.Yo=0.5*cHeightPels;
			vprint.gfactor*=pfactor;

			hdc=GetDC((wdraw.childs+1)->hwnd);
			hDC  = CreateCompatibleDC(hdc);
			hbit = CreateCompatibleBitmap(hdc,cWidthPels,
											  cHeightPels);
			pbit = (HBITMAP)SelectObject(hDC,hbit);

			hbrush = (HBRUSH)GetStockObject(WHITE_BRUSH);
			/*hbrush = (HBRUSH)GetStockObject(BLACK_BRUSH);*/
			pbrush = (HBRUSH)SelectObject(hDC,hbrush);
			PatBlt(hDC,0,0,cWidthPels,cHeightPels,PATCOPY);

			setfontformat(hDC,prn.jiheight,prn.jiwidth,
						  "MS Mincho",0,0,0);
			SetBkMode(hDC,TRANSPARENT); /*BACKGROUND OF TEXT.*/

			if((wmenu.childs+2)->vparam.vflag.mv.ftype==F_ORGAN &&
			   (wdraw.childs+1)->org.elems!=NULL)
			{
			  draworganization(hDC,vprint,
							   (wdraw.childs+1)->org,ONPREVIEW);
			  /*draworganization(hDC,vprint,
							   (wdraw.childs+1)->org,ONSCREEN);*/
			}
			if((wmenu.childs+2)->vparam.vflag.mv.ftype==F_ARCLM &&
			   arc.elems!=NULL)
			{
			  drawarclmframe(hDC,vprint,arc,0,ONPREVIEW);
			}

			setfontformat(hDC,(int)(1.5*prn.jiheight),
							  (int)(1.5*prn.jiwidth),
							  "MS Mincho",0,0,0);

			for(i=0;i<wdraw.nstring;i++)
			{
			  (wdraw.strset+i)->n.d[0]*=pfactor;
			  (wdraw.strset+i)->n.d[1]*=pfactor;
			}
			drawtexts(hDC,wdraw.strset,wdraw.nstring,vprint);
			for(i=0;i<wdraw.nstring;i++)
			{
			  (wdraw.strset+i)->n.d[0]/=pfactor;
			  (wdraw.strset+i)->n.d[1]/=pfactor;
			}

			hfont=(HFONT)GetCurrentObject(hDC,OBJ_FONT);
			DeleteObject(hfont);

			/*ON SCREEN*/
			BitBlt(hdc,0,0,cWidthPels,cHeightPels,
				   hDC,0,0,SRCCOPY);
			ReleaseDC((wdraw.childs+1)->hwnd,hdc);

#if 0
  nblank=0;

  file = fopen("preview.bmp","wb");

  fhead.bfType = MAKEWORD(l,h);
  fhead.bfSize = sizeof(BITMAPFILEHEADER)
			   + sizeof(BITMAPINFOHEADER)
			   + cWidthPels*cHeightPels*sizeof(RGBTRIPLE)
			   + nblank*sizeof(short int)*cHeightPels;
  fhead.bfReserved1 = 0;
  fhead.bfReserved2 = 0;
  fhead.bfOffBits = sizeof(BITMAPFILEHEADER)
				  + sizeof(BITMAPINFOHEADER);

  ihead.biSize = sizeof(BITMAPINFOHEADER);
  ihead.biWidth = cWidthPels;
  ihead.biHeight = cHeightPels;
  ihead.biPlanes = 1;
  ihead.biBitCount = 24;
  ihead.biCompression = BI_RGB;
  ihead.biSizeImage = cWidthPels*cHeightPels*sizeof(RGBTRIPLE);
  /*ihead.biSizeImage = 24*(width+7)/8*height;*/ /*+ 2*height*/

  hdesktop = GetParent(hwnd);
  hdc = GetDC(hdesktop);
  px = GetDeviceCaps(hdc,HORZRES);
  mt = GetDeviceCaps(hdc,HORZSIZE);
  ihead.biXPelsPerMeter = 1000*px/mt;

  px = GetDeviceCaps(hdc,VERTRES);
  mt = GetDeviceCaps(hdc,VERTSIZE);
  ihead.biYPelsPerMeter = 1000*px/mt;

  ReleaseDC(hdc,hdesktop);

  ihead.biClrUsed = 0;
  ihead.biClrImportant = 0;

  fwrite(&fhead,sizeof(BITMAPFILEHEADER),1,file);
  fwrite(&ihead,sizeof(BITMAPINFOHEADER),1,file);

  for(i=cHeightPels;i>=1;i--)
  {
	for(j=1;j<=cWidthPels;j++)
	{
	  cr=GetPixel(hDC,j,i);
	  rgbtri.rgbtRed   = GetRValue(cr);
	  rgbtri.rgbtGreen = GetGValue(cr);
	  rgbtri.rgbtBlue  = GetBValue(cr);

	  fwrite(&rgbtri,sizeof(RGBTRIPLE),1,file);
	}
	for(j=0;j<nblank;j++)
	{
	  fwrite(&blank,sizeof(short int),1,file);
	}
  }

  fclose(file);
  MessageBox(NULL,"Saved.","SAVEASBMP",MB_OK);
#endif
/********/
/*hBmp=hbit;*/

	/* Retrieve the bitmap's color format, width, and height. */
	GetObject(hbit, sizeof(BITMAP), (LPSTR)&bmp);

sprintf(str,"Size=%d x %d",bmp.bmWidth,bmp.bmHeight);
MessageBox(NULL,str,"Preview",MB_OK);

	/* Convert the color format to a count of bits. */
	cClrBits = (WORD)(bmp.bmPlanes * bmp.bmBitsPixel);

sprintf(str,"Planes x Bits =%d x %d",bmp.bmPlanes,bmp.bmBitsPixel);
MessageBox(NULL,str,"Preview",MB_OK);

	if (cClrBits == 1)       cClrBits = 1;
	else if (cClrBits <= 4)  cClrBits = 4;
	else if (cClrBits <= 8)  cClrBits = 8;
	else if (cClrBits <= 16) cClrBits = 16;
	else if (cClrBits <= 24) cClrBits = 24;
	else                     cClrBits = 32;

	/* Allocate memory for the BITMAPINFO structure.*/
	/* This structure contains a BITMAPINFOHEADER structure */
	/* and an array of RGBQUAD data structures.*/

	if (cClrBits != 24) pbmi = (PBITMAPINFO) LocalAlloc(LPTR,
								sizeof(BITMAPINFOHEADER) +
								sizeof(RGBQUAD) * (2^cClrBits));

    /* There is no RGBQUAD array for the 24-bit-per-pixel format.*/
    else                pbmi = (PBITMAPINFO) LocalAlloc(LPTR,
                                sizeof(BITMAPINFOHEADER));

    /* Initialize the fields in the BITMAPINFO structure. */

    pbmi->bmiHeader.biSize = sizeof(BITMAPINFOHEADER);
    pbmi->bmiHeader.biWidth = bmp.bmWidth;
    pbmi->bmiHeader.biHeight = bmp.bmHeight;
    pbmi->bmiHeader.biPlanes = bmp.bmPlanes;
    pbmi->bmiHeader.biBitCount = bmp.bmBitsPixel;

    if (cClrBits < 24) pbmi->bmiHeader.biClrUsed = 2^cClrBits;

    /* If the bitmap is not compressed, set the BI_RGB flag. */
    pbmi->bmiHeader.biCompression = BI_RGB;

    /* Compute the number of bytes in the array of color */
    /* indices and store the result in biSizeImage.*/
    pbmi->bmiHeader.biSizeImage = (pbmi->bmiHeader.biWidth + 7) /8
                                  * pbmi->bmiHeader.biHeight
                                  * cClrBits;

    /* Set biClrImportant to 0, indicating that all of the */
    /* device colors are important.*/

    pbmi->bmiHeader.biClrImportant = 0;

/*pbi=pbmi;*/
/*hBMP=hBmp;*/
/*hDC=pd.hDC;*/

    pbih = (PBITMAPINFOHEADER) pbmi;
    lpBits = (LPBYTE) GlobalAlloc(GMEM_FIXED, pbih->biSizeImage);

    /* Retrieve the color table (RGBQUAD array) and the bits */
    /* (array of palette indices) from the DIB.*/

    GetDIBits(hDC, hbit, 0, (WORD) pbih->biHeight,
              lpBits, pbmi, DIB_RGB_COLORS);

    /* Create the .BMP file. */
    hf = CreateFile("preview.bmp",
                    GENERIC_READ | GENERIC_WRITE,
                    (DWORD) 0,
                    (LPSECURITY_ATTRIBUTES) NULL,
                    CREATE_ALWAYS,
                    FILE_ATTRIBUTE_NORMAL,
                    (HANDLE) NULL);

    hdr.bfType = 0x4d42;        /* 0x42 = "B" 0x4d = "M" */

    /* Compute the size of the entire file. */
    hdr.bfSize = (DWORD) (sizeof(BITMAPFILEHEADER) +
                 pbih->biSize + pbih->biClrUsed
                 * sizeof(RGBQUAD) + pbih->biSizeImage);

    hdr.bfReserved1 = 0;
    hdr.bfReserved2 = 0;

    /* Compute the offset to the array of color indices. */
    hdr.bfOffBits = (DWORD) sizeof(BITMAPFILEHEADER)
                          + pbih->biSize + pbih->biClrUsed
                          * sizeof (RGBQUAD);

    /* Copy the BITMAPFILEHEADER into the .BMP file. */
    WriteFile(hf, (LPVOID) &hdr, sizeof(BITMAPFILEHEADER),
              (LPDWORD) &dwTmp, (LPOVERLAPPED) NULL);

    /* Copy the BITMAPINFOHEADER and RGBQUAD array into the file. */
    WriteFile(hf, (LPVOID) pbih,
                  sizeof(BITMAPINFOHEADER)
				  + pbih->biClrUsed * sizeof (RGBQUAD),
                  (LPDWORD) &dwTmp, (LPOVERLAPPED) NULL);

    /* Copy the array of color indices into the .BMP file. */
    dwTotal = cb = pbih->biSizeImage;
    hp = lpBits;
    while (cb > MAXWRITE)
    {
            WriteFile(hf, (LPSTR) hp, (int) MAXWRITE,
                      (LPDWORD) &dwTmp, (LPOVERLAPPED) NULL);
            cb-= MAXWRITE;
            hp += MAXWRITE;
    }
    WriteFile(hf, (LPSTR) hp, (int) cb,
              (LPDWORD) &dwTmp, (LPOVERLAPPED) NULL);

    /* Close the .BMP file. */
    CloseHandle(hf);

    /* Free memory. */
    GlobalFree((HGLOBAL)lpBits);
/********/
    SelectObject(hDC,pbit);
    DeleteObject(hbit);
    SelectObject(hDC,pbrush);
    /*DeleteObject(hbrush);*/
    DeleteDC(hDC); /*CREATECOMPATIBLEDC:DELETE GETDC:RELEASE*/
    /*ReleaseDC((wdraw.childs+1)->hwnd,hDC);*/ /*FOR GETDC*/
/********/
            MessageBox(NULL,"Completed.","Print",MB_OK);
          }
          break;

        case IDM_MACROS: /*MACRO PROCEDURE.*/
          ip=0; /*PRINT FLAG*/

          /*PRINT MODEL*/
          sprintf((wdraw.childs+1)->inpfile,"hako\\hako41.inl");
          SetDlgItemText((wmenu.childs+2)->hwnd,ID_INPUTFILE,
                         (wdraw.childs+1)->inpfile);
          (wdraw.childs+1)->vparam.gfactor=3.93;
          (wdraw.childs+1)->vparam.focus.d[0]=63.4;
		  (wdraw.childs+1)->vparam.focus.d[1]=50.0;
          (wdraw.childs+1)->vparam.focus.d[2]=-10.0;
          (wdraw.childs+1)->vparam.theta=-180.0;
          (wdraw.childs+1)->vparam.phi=0.0;
          (wdraw.childs+1)->vparam.range.max.d[GX]=1000.0;
          (wdraw.childs+1)->vparam.range.min.d[GX]=-100.0;
          (wdraw.childs+1)->vparam.range.max.d[GY]=1000.0;
          (wdraw.childs+1)->vparam.range.min.d[GY]=-100.0;
          (wdraw.childs+1)->vparam.range.max.d[GZ]=1000.0;
          (wdraw.childs+1)->vparam.range.min.d[GZ]=-100.0;
          (wdraw.childs+1)->vparam.dparam.hsize=3;
          (wdraw.childs+1)->vparam.dparam.csize=3;      /*byMIHARA for Conffig*/

          SendMessage((wmenu.childs+2)->hwnd,WM_INITDIALOG,0,0);

          wparam = MAKEWPARAM(IDF_FRAME,0);
          SendMessage((wmenu.childs+2)->hwnd,WM_COMMAND,wparam,0);
          wparam = MAKEWPARAM(IDD_VIEW,0);
          SendMessage((wmenu.childs+2)->hwnd,WM_COMMAND,wparam,0);

          wparam = MAKEWPARAM(IDV_NODECODE,0); /*CLEAR*/
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,wparam,0);
/*mihara for zahyo 20190329///////////////////////////////////////////////////*/
          wparam = MAKEWPARAM(IDV_GX, 0); /* CLEAR */
          SendMessage((wmenu.childs + 3)->hwnd, WM_COMMAND, wparam, 0);
          wparam = MAKEWPARAM(IDV_GY, 0); /* CLEAR */
          SendMessage((wmenu.childs + 3)->hwnd, WM_COMMAND, wparam, 0);
          wparam = MAKEWPARAM(IDV_GZ, 0); /* CLEAR */
          SendMessage((wmenu.childs + 3)->hwnd, WM_COMMAND, wparam, 0);
/*mihara for zahyo 20190329///////////////////////////////////////////////////*/

          if(ip)
          {
            wparam = MAKEWPARAM((WORD)IDM_PRINTSETUP,(WORD)0);
            SendMessage(wmain.hwnd,WM_COMMAND,wparam,0);
          }
          wdraw.strset=addtext(wdraw.strset,&(wdraw.nstring),
                               "Hakodate Mirai Daigaku",
                               -250.0,-3.0);
          wdraw.strset=addtext(wdraw.strset,&(wdraw.nstring),
                               "Model         , Scale=1:600",
                               -250.0,5.0);
		  wdraw.strset=addtext(wdraw.strset,&(wdraw.nstring),
                               "Node,Element Code",-250.0,13.0);


          wdraw.strset=addtext(wdraw.strset,&(wdraw.nstring),
                               "Level 1",210.0, -28.0);
          wdraw.strset=addtext(wdraw.strset,&(wdraw.nstring),
                               "Level 2",210.0, -46.0);
          wdraw.strset=addtext(wdraw.strset,&(wdraw.nstring),
                               "Level 3",210.0, -62.0);

          wdraw.strset=addtext(wdraw.strset,&(wdraw.nstring),
                               "Y3 ",(-251.0+49.2*7.0),-14.0);
          wdraw.strset=addtext(wdraw.strset,&(wdraw.nstring),
                               "Y2 ",(-251.0+49.2*8.0),-14.0);
          wdraw.strset=addtext(wdraw.strset,&(wdraw.nstring),
                               "Y1 ",(-251.0+49.2*9.0),-14.0);

          /*STREET X*/
          /*printmodelmacro(-12.0,-13.0,1000.0,-100.0,
                          IDV_NODECODE,IDV_ELEMENTCODE,0,
                          IDV_SECTIONCODE,0,0,"X1");
          printmodelmacro(  1.0, -1.0,1000.0,-100.0,
                          IDV_NODECODE,IDV_ELEMENTCODE,0,
                          IDV_SECTIONCODE,0,0,"X2");
          printmodelmacro( 13.0, 12.0,1000.0,-100.0,
                          IDV_NODECODE,IDV_ELEMENTCODE,0,
                          IDV_SECTIONCODE,0,0,"X3");
          */
          /*PLAN*/
          pfactor=0.75;
          (wdraw.childs+1)->vparam.focus.d[0]=63.4-12.6;
          (wdraw.childs+1)->vparam.focus.d[1]=50.0+5.0;
          (wdraw.childs+1)->vparam.gfactor=3.93*pfactor;
          (wdraw.childs+1)->vparam.theta=-90.0;
          (wdraw.childs+1)->vparam.phi=90.0;
          (wdraw.childs+1)->vparam.range.max.d[GX]=1000.0;
          (wdraw.childs+1)->vparam.range.min.d[GX]=-100.0;
          SendMessage((wmenu.childs+2)->hwnd,WM_INITDIALOG,0,0);

          wparam = MAKEWPARAM(IDV_NODECODE,0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,wparam,0);
 /*mihara for zahyo 20190329///////////////////////////////////////////////////*/
			wparam = MAKEWPARAM(IDV_GX, 0);
			SendMessage((wmenu.childs + 3)->hwnd, WM_COMMAND, wparam, 0);
			wparam = MAKEWPARAM(IDV_GY, 0);
			SendMessage((wmenu.childs + 3)->hwnd, WM_COMMAND, wparam, 0);
			wparam = MAKEWPARAM(IDV_GZ, 0);
			SendMessage((wmenu.childs + 3)->hwnd, WM_COMMAND, wparam, 0);
/*mihara for zahyo 20190329///////////////////////////////////////////////////*/

          sprintf((wdraw.strset+1)->str,"Level 1, Scale=1:800");
          sprintf((wdraw.strset+2)->str,"Node Code");

          (wdraw.strset+0)->n.d[0]=-205;
          (wdraw.strset+1)->n.d[0]=-205;
          (wdraw.strset+2)->n.d[0]=-205;
          (wdraw.strset+0)->n.d[1]=192;
          (wdraw.strset+1)->n.d[1]=200;
          (wdraw.strset+2)->n.d[1]=208;

          (wdraw.strset+6)->n.d[0]=-150.0;
          (wdraw.strset+7)->n.d[0]=-150.0+pfactor*48.7*1.0;
          (wdraw.strset+8)->n.d[0]=-150.0+pfactor*48.7*2.0;

          (wdraw.strset+6)->n.d[1]=200-19;
          (wdraw.strset+7)->n.d[1]=200-19;
          (wdraw.strset+8)->n.d[1]=200-19;

          wdraw.strset=addtext(wdraw.strset,&(wdraw.nstring),
                               "Y9 ",-205,(-173.0+pfactor*48.7*1.0));
          wdraw.strset=addtext(wdraw.strset,&(wdraw.nstring),
                               "Y8 ",-205,(-173.0+pfactor*48.7*2.0));
          wdraw.strset=addtext(wdraw.strset,&(wdraw.nstring),
                               "Y7 ",-205,(-173.0+pfactor*48.7*3.0));

          wparam = MAKEWPARAM(IDD_VIEW,0);
          SendMessage((wmenu.childs+2)->hwnd,WM_COMMAND,wparam,0);

          /*PRINT STRESSES*/
          /*LOAD Z*/
          /*sprintf((wdraw.childs+1)->inpfile,"hako\\hako41.inl");
          sprintf((wdraw.childs+1)->otpfile,"hako\\hako41.otl");
		  SetDlgItemText((wmenu.childs+2)->hwnd,ID_INPUTFILE,
						 (wdraw.childs+1)->inpfile);
          SetDlgItemText((wmenu.childs+2)->hwnd,ID_OUTPUTFILE,
                         (wdraw.childs+1)->otpfile);
          (wdraw.childs+1)->vparam.gfactor=3.93;
          (wdraw.childs+1)->vparam.focus.d[0]=63.4;
          (wdraw.childs+1)->vparam.focus.d[1]=50.0;
          (wdraw.childs+1)->vparam.focus.d[2]=-10.0;
          (wdraw.childs+1)->vparam.theta=-180.0;
          (wdraw.childs+1)->vparam.phi=0.0;
          (wdraw.childs+1)->vparam.range.max.d[GX]=1000.0;
          (wdraw.childs+1)->vparam.range.min.d[GX]=-100.0;
          (wdraw.childs+1)->vparam.range.max.d[GY]=1000.0;
          (wdraw.childs+1)->vparam.range.min.d[GY]=-100.0;
          (wdraw.childs+1)->vparam.range.max.d[GZ]=1000.0;
          (wdraw.childs+1)->vparam.range.min.d[GZ]=-100.0;
          (wdraw.childs+1)->vparam.dparam.hsize=3;
          SendMessage((wmenu.childs+2)->hwnd,WM_INITDIALOG,0,0);

          wparam = MAKEWPARAM(IDF_FRAME,0);
          SendMessage((wmenu.childs+2)->hwnd,WM_COMMAND,wparam,0);
          wparam = MAKEWPARAM(IDD_VIEW,0);
          SendMessage((wmenu.childs+2)->hwnd,WM_COMMAND,wparam,0);
          wparam = MAKEWPARAM(IDD_OPENRESULT,0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,wparam,0);

          wparam = MAKEWPARAM(IDV_NODECODE,0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,wparam,0);

          if(ip)
          {
            wparam = MAKEWPARAM((WORD)IDM_PRINTSETUP,(WORD)0);
            SendMessage(wmain.hwnd,WM_COMMAND,wparam,0);
          }
          wdraw.strset=addtext(wdraw.strset,&(wdraw.nstring),
                               "Hakodate Mirai Daigaku",
                               -250.0,-3.0);
          wdraw.strset=addtext(wdraw.strset,&(wdraw.nstring),
                               "Load Y, Model         , Scale=1:600",
							   -250.0,5.0);

          wdraw.strset=addtext(wdraw.strset,&(wdraw.nstring),
                               "Level 1",210.0, -28.0);
		  wdraw.strset=addtext(wdraw.strset,&(wdraw.nstring),
                               "Level 2",210.0, -46.0);
          wdraw.strset=addtext(wdraw.strset,&(wdraw.nstring),
                               "Level 3",210.0, -62.0);

          wdraw.strset=addtext(wdraw.strset,&(wdraw.nstring),
                               "Y3 ",(-251.0+49.2*7.0),-14.0);
          wdraw.strset=addtext(wdraw.strset,&(wdraw.nstring),
                               "Y2 ",(-251.0+49.2*8.0),-14.0);
          wdraw.strset=addtext(wdraw.strset,&(wdraw.nstring),
                               "Y1 ",(-251.0+49.2*9.0),-14.0);
          */
          /*LOAD Z:Q,M FOR EACH STREET X*/
          /*printmacro(-12.0,-13.0,1000.0,-100.0,
                     0,0,0,IDV_QX,IDV_QY_G,0,
                     IDV_MY,IDV_MX_G,0,"Z","X1");
          printmacro(  1.0, -1.0,1000.0,-100.0,
                     IDV_MY,IDV_MX_G,0,IDV_QX,IDV_QY_G,0,
                     IDV_MY,IDV_MX_G,0,"Z","X2");
          printmacro( 13.0, 12.0,1000.0,-100.0,
                     IDV_MY,IDV_MX_G,0,IDV_QX,IDV_QY_G,0,
                     IDV_MY,IDV_MX_G,0,"Z","X3");
          */

          /*LOAD Z:DEFORMATION,N FOR EACH STREET X*/
          /*printmacro(-12.0,-13.0,1000.0,-100.0,
                     IDV_MY,IDV_MX_G,0,
                     IDV_DEFORMATION,IDV_DZ,0,
                     IDV_NZ,IDV_NZ_G,IDV_NZ_B,"Z","X1");
          printmacro(  1.0, -1.0,1000.0,-100.0,
                     IDV_NZ,IDV_NZ_G,IDV_NZ_B,
                     IDV_DEFORMATION,IDV_DZ,0,
                     IDV_NZ,IDV_NZ_G,IDV_NZ_B,"Z","X2");
          printmacro( 13.0, 12.0,1000.0,-100.0,
                     IDV_NZ,IDV_NZ_G,IDV_NZ_B,
                     IDV_DEFORMATION,IDV_DZ,0,
                     IDV_NZ,IDV_NZ_G,IDV_NZ_B,"Z","X3");
          */

          /*END PRINTING.*/
          if(ip)
          {
			wparam = MAKEWPARAM((WORD)IDM_PRINTEND,(WORD)0);
            SendMessage(wmain.hwnd,WM_COMMAND,wparam,0);
          }
          break;

        case IDM_CALCULATIONMACROI: /*MACRO PROCEDURE.*/
		  if(wdraw.nchilds>=2 &&
             (wdraw.childs+1)->hwnd!=NULL &&
			 (wmenu.childs+2)->vparam.vflag.mv.ftype==F_ORGAN &&
			 (wdraw.childs+1)->org.nodes!=NULL)
		  {
            globalmessageflag=0;

            getincrement((wmenu.childs+2)->hwnd,&laps,&dsafety);

			if(extractarclmfromorgan(&((wdraw.childs+1)->org),
									 &arc,&arcx,&arcy)==0)
			{
			  break;
			}

			(wmenu.childs+2)->vparam.vflag.mv.ftype=F_ARCLM;
			SendMessage((wmenu.childs+2)->hwnd,WM_INITDIALOG,0,0);

			clearwindow(*(wdraw.childs+1));
			drawarclmframe((wdraw.childs+1)->hdcC,
						   (wdraw.childs+1)->vparam,arc,0,ONSCREEN);
			SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);

//			MessageBox(NULL,"Completed.","Organ Into Arclm",MB_OK);

			n=strcspn((wdraw.childs+1)->inpfile,".");
			strncpy(str,(wdraw.childs+1)->inpfile,n);
            str[n]='\0';

			strcpy((wdraw.childs+1)->inpfilez,str);
			strcat((wdraw.childs+1)->inpfilez,".inl");
            strcpy((wdraw.childs+1)->inpfilex,str);
            strcat((wdraw.childs+1)->inpfilex,".ihx");
            strcpy((wdraw.childs+1)->inpfiley,str);
            strcat((wdraw.childs+1)->inpfiley,".ihy");

			strcpy((wdraw.childs+1)->otpfile,str);
            strcat((wdraw.childs+1)->otpfile,".otp");
            strcpy((wdraw.childs+1)->otpfilez,str);
            strcat((wdraw.childs+1)->otpfilez,".otl");
            strcpy((wdraw.childs+1)->otpfilex,str);
            strcat((wdraw.childs+1)->otpfilex,".ohx");
			strcpy((wdraw.childs+1)->otpfiley,str);
            strcat((wdraw.childs+1)->otpfiley,".ohy");

//            strcpy((wdraw.childs+1)->sctfile,str);
//            strcat((wdraw.childs+1)->sctfile,".lst");

            saveasarclm((wdraw.childs+1)->inpfilez,&arc);
			saveasarclm((wdraw.childs+1)->inpfilex,&arcx);
			saveasarclm((wdraw.childs+1)->inpfiley,&arcy);

            sprintf((wdraw.childs+1)->inpfile,
                    (wdraw.childs+1)->inpfilez);

            SetDlgItemText((wmenu.childs+2)->hwnd,ID_INPUTFILE,
                           (wdraw.childs+1)->inpfile);
            SetDlgItemText((wmenu.childs+2)->hwnd,ID_INPUTFILEZ,
                           (wdraw.childs+1)->inpfilez);
            SetDlgItemText((wmenu.childs+2)->hwnd,ID_INPUTFILEX,
                           (wdraw.childs+1)->inpfilex);
            SetDlgItemText((wmenu.childs+2)->hwnd,ID_INPUTFILEY,
                           (wdraw.childs+1)->inpfiley);

			SetDlgItemText((wmenu.childs+2)->hwnd,ID_OUTPUTFILE,
                           (wdraw.childs+1)->otpfile);
            SetDlgItemText((wmenu.childs+2)->hwnd,ID_OUTPUTFILEZ,
                           (wdraw.childs+1)->otpfilez);
            SetDlgItemText((wmenu.childs+2)->hwnd,ID_OUTPUTFILEX,
                           (wdraw.childs+1)->otpfilex);
            SetDlgItemText((wmenu.childs+2)->hwnd,ID_OUTPUTFILEY,
                           (wdraw.childs+1)->otpfiley);

//            SetDlgItemText((wmenu.childs+2)->hwnd,ID_SECTIONFILE,
//                           (wdraw.childs+1)->sctfile);

			saveasarclm((wdraw.childs+1)->inpfilez,&arc);
//			saveasarclm((wdraw.childs+1)->inpfilex,&arcx);
//			saveasarclm((wdraw.childs+1)->inpfiley,&arcy);

			getviewparam((wmenu.childs+2)->hwnd,
						 &((wdraw.childs+1)->vparam));

            /*STATIC LINEAR ANALYSIS*/
			arclm001(&arc ,ID_INPUTFILEZ,ID_OUTPUTFILEZ);

            /*LINEAR BUCKLING ANALYSIS*/
            bclng001(&arc);

            /*BUCKLIMG CONDENSATION*/
            strcpy((wdraw.childs+1)->otpfile,str);
            strcat((wdraw.childs+1)->otpfile,".otbc");
			SetDlgItemText((wmenu.childs+2)->hwnd,ID_OUTPUTFILE,
                           (wdraw.childs+1)->otpfile);

            bclng002(&arc);

			/*LOAD INCREMENTAL ANALYSIS(MATERIAL NONLINEAR)*/
            strcpy((wdraw.childs+1)->otpfile,str);
            strcat((wdraw.childs+1)->otpfile,".otl2");
//            strcat((wdraw.childs+1)->otpfile,".ohy2");
            SetDlgItemText((wmenu.childs+2)->hwnd,ID_OUTPUTFILE,
                           (wdraw.childs+1)->otpfile);

//            SetDlgItemText((wmenu.childs+2)->hwnd,ID_LAPS,"100");
//            SetDlgItemText((wmenu.childs+2)->hwnd,ID_SAFETY,"0.01");
            sprintf(str,"%d",laps);
            SetDlgItemText((wmenu.childs+2)->hwnd,ID_LAPS,str);
            sprintf(str,"%.3f",dsafety);
            SetDlgItemText((wmenu.childs+2)->hwnd,ID_SAFETY,str);

			getviewparam((wmenu.childs+2)->hwnd,&((wdraw.childs+1)->vparam));
            (wdraw.childs+1)->vparam.vflag.ev.srcanrate=0;

            clearwindow(*(wdraw.childs+1));

//	        arclm101_bc(&arc,ID_INPUTFILE);
	        arclm101_bc(&arc,ID_INPUTFILEX);

            globalmessageflag=1;

          }
          break;

        case IDM_OPENSECTION: /*OPEN SRCAN SECTION LIST.*/
          if(wdraw.nchilds>=2 &&
             (wdraw.childs+1)->hwnd!=NULL &&
             (wmenu.childs+2)->vparam.vflag.mv.ftype==F_SRCAN)
          {
            if(srcansects>0) free(srcanlist);

            /*sprintf((wdraw.childs+1)->inpfile,"srcal01.lst");
            SetDlgItemText((wmenu.childs+2)->hwnd,ID_INPUTFILE,
                           (wdraw.childs+1)->inpfile);*/

            fin=fgetstofopen(dir,"r",ID_SECTIONFILE);
            /*fin=fopen((wdraw.childs+1)->inpfile,"r");*/
            if(fin==NULL) break;

            sprintf(str,"OPENED=%s",(wdraw.childs+1)->sctfile);
            errormessage(str);

            srcansects=getcodelist(fin,codelist);/*LIST UP CODES.*/
            if(srcansects>MAXSECT) break;
            srcanlist=(struct section *)
                      malloc(srcansects*sizeof(struct section));
            if(srcanlist==NULL) return 0;
            for(i=0;i<srcansects;i++) /*SECTIONS INTO MEMORY.*/
            {
              getsectionform(fin,codelist[i],(srcanlist+i));
            }

            fclose(fin);

            SendMessage(wmain.hwnd,WM_COMMAND,
						MAKEWPARAM((WORD)IDM_DRAWSECTION,(WORD)0),
                        (WPARAM)0);
          }
          break;

        case IDM_DRAWSECTION: /*SECTION LIST.*/
          if(wdraw.nchilds>=2 &&
			 (wdraw.childs+1)->hwnd!=NULL &&
             (wmenu.childs+2)->vparam.vflag.mv.ftype==F_SRCAN &&
             srcansects>0)
          {
            getviewparam((wmenu.childs+2)->hwnd,
                         &((wdraw.childs+1)->vparam));
            getclientsize((wdraw.childs+1)->hwnd,&maxX,&maxY);

            clearwindow(*(wdraw.childs+1));
            drawsectionlist((wdraw.childs+1)->hdcC,
                            srcansects,srcanlist,
                            (wdraw.childs+1)->vparam.Xo,
							(wdraw.childs+1)->vparam.Yo,
                            maxX,maxY,
                            (wdraw.childs+1)->vparam,ONSCREEN);
            SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);
          }
          break;

        case IDM_PRINTSECTION: /*SECTION LIST.*/
          if(wdraw.nchilds>=2 &&
             (wdraw.childs+1)->hwnd!=NULL &&
             (wmenu.childs+2)->vparam.vflag.mv.ftype==F_SRCAN &&
             srcansects>0)
          {
            getviewparam((wmenu.childs+2)->hwnd,
                         &((wdraw.childs+1)->vparam));

            vprint=(wdraw.childs+1)->vparam;
            vprint.gfactor=10.0;

            /*BEGIN PRINTING.*/
            wparam = MAKEWPARAM((WORD)IDM_PRINTSETUP,(WORD)0);
            SendMessage(wmain.hwnd,WM_COMMAND,wparam,0);

            setfontformat(gprn.pd.hDC,80,30,"MS Mincho",0,0,0);
            SetBkMode(gprn.pd.hDC,TRANSPARENT);

            drawsectionlist(gprn.pd.hDC,
                            srcansects,srcanlist,
                            1000,500,
                            gprn.cWidthPels,
							gprn.cHeightPels,
                            vprint,ONPRINTER);

            /*END PRINTING.*/
            wparam = MAKEWPARAM((WORD)IDM_PRINTEND,(WORD)0);
            SendMessage(wmain.hwnd,WM_COMMAND,wparam,0);
          }
          break;

		case IDM_SRCANBEGIN:
          if(MessageBox(NULL,"Srcan001:Begin.","SRCAN001",
             MB_OKCANCEL)==IDCANCEL) break;

          if(wdraw.nchilds>=2 &&
             (wdraw.childs+1)->hwnd!=NULL &&
             (wmenu.childs+2)->vparam.vflag.mv.ftype==F_ARCLM)
          {
            n=strcspn((wdraw.childs+1)->sctfile,".");
            strncpy(str,(wdraw.childs+1)->sctfile,n);
            str[n]='\0';

            i=srcan001(str);

            if(i==0) MessageBox(NULL,"Failed.","SRCAN001",MB_OK);
            else     MessageBox(NULL,"Succeeded.","SRCAN001",MB_OK);
          }
          break;
        case IDM_CREATESURFACE:
          if(wdraw.nchilds>=2 &&
             (wdraw.childs+1)->hwnd!=NULL &&
             (wmenu.childs+2)->vparam.vflag.mv.ftype==F_ARCLM)
          {
            createyieldsurface(&arc);
          }
          break;

        case IDM_REDELEMS:
          if(wdraw.nchilds>=2 &&
             (wdraw.childs+1)->hwnd!=NULL &&
             (wmenu.childs+2)->vparam.vflag.mv.ftype==F_ORGAN ||
             (wmenu.childs+2)->vparam.vflag.mv.ftype==F_ARCLM ||
             (wmenu.childs+2)->vparam.vflag.mv.ftype==F_FRAME)
		  {
            fout=fopen((wdraw.childs+1)->otpfile,"r");      /*.ratを開きます。*/
            if(fout==NULL) break;
            readsrcanrate(fout,&arc);                   /*.ratを読み込みます。*/
            fclose(fout);
            int nred=0;                 /*nred(NG部材のカウンター)の値の初期化*/
            double maxsafety=0.0;          /*maxsafety(安全率の最大値)の初期化*/
            int ii,jj;
            for(ii=0;ii<arc.nelem;ii++)
            {

              //if((arc.elems+ii)->srate[2]>1) nred=nred+1;
              if( (arc.elems+ii)->srate[0]>1 ||                 /*部材の安全率*/
                  (arc.elems+ii)->srate[1]>1 ||
                  (arc.elems+ii)->srate[2]>1 ||
                  (arc.elems+ii)->srate[3]>1 )
                {
                  nred=nred+1;            					  /*nredの値の更新*/
                }

				for(jj=0;jj<4;jj++)
				{
				  if(maxsafety<(arc.elems+ii)->srate[jj]) maxsafety=(arc.elems+ii)->srate[jj];
                                                         /*maxsafetyの値の更新*/
				}

            }
            sprintf(str,"NG ELEMS: %d 個 ( MAX SAFETY= %9.5f )",nred,maxsafety);
            MessageBox(NULL,str,"NG_ELEMS",MB_OK);          /*NG部材の数を表示*/
          }
          break;

		case IDM_HELP:                               /*HELP TEXT.*/
		  MessageBox(wmain.hwnd,"No Help","Help",MB_OK);
		  break;

		default: /*OTHERS.*/
          return DefWindowProc(hwnd,message,wParam,lParam);
      }
      break;

	case WM_LBUTTONDOWN:
	  hdc = GetDC(hwnd);
      x = LOWORD(lParam);
      y = HIWORD(lParam);
      sprintf(str,"POSITION:%d %d",x,y);
	  SetBkMode(hdc,TRANSPARENT);                 /*TEXT BACKGROUND*/
      TextOut(hdc,x,y,str,strlen(str));
      ReleaseDC(hwnd,hdc);               /*CASE REDRAW UNAVAILABLE.*/
      /*InvalidateRect(hwnd,NULL,TRUE);*/            /*CASE REDRAW.*/
      break;

    case WM_DESTROY:
      PostQuitMessage(0);                       /*CREATE "WM_QUIT".*/
      break;

    default:
      return DefWindowProc(hwnd,message,wParam,lParam);
  }
  return 0;
}/*WindowProcedureMain*/

LRESULT CALLBACK WindowProcedureSheet(HWND hwnd,
                                      UINT message,
                                      WPARAM wParam,
                                      LPARAM lParam)
/*SHEET WINDOW WITH VSCROLL,HSCROLL.*/
{
  HDC hdc;
  HBRUSH hbrush;
  POINT point,cp;
  int bl,bt,br,bb,bmax; /*LEFT,TOP,RIGHT,BOTTOM OF BAR.*/
  long int cw,ch,pw,ph,maxX,maxY,wintop,winleft; /*WINDOW SIZE.*/
  struct windowparams *wp;
  /*Scroll with Wheel*/
  int rotation;
  int wd;
  double scrollunit;
  char str[256];
  int hflag;
  /*Scroll with Wheel*/

  switch(message)
  {
	case WM_PAINT:
    case WM_SIZE:
	  DefWindowProc(hwnd,message,wParam,lParam);
	  /*EnumChildWindows(hwnd,(DLGPROC)EnumChildProcSheet,0);*/ /*UPDATE CHILDS.*/
	  EnumChildWindows(hwnd,EnumChildProcSheet,0); /*UPDATE CHILDS.*/

      getclientsize(hwnd,&maxX,&maxY);

      hdc=GetDC(hwnd);

	  hbrush = (HBRUSH)GetClassLong(hwnd,GCL_HBRBACKGROUND);
      SelectObject(hdc,hbrush);
      PatBlt(hdc,0,0,maxX,maxY,PATCOPY);

      DrawSunken(hdc,0,0,(maxX-BARWIDTH-4),(maxY-BARWIDTH-4));
      DrawSunken(hdc,(maxX-BARWIDTH-2),0,(maxX-1),(maxY-BARWIDTH-4));
      DrawSunken(hdc,0,(maxY-BARWIDTH-2),(maxX-BARWIDTH-4),(maxY-1));
      DrawSunken(hdc,(maxX-BARWIDTH-2),(maxY-BARWIDTH-2),
                     (maxX-1),(maxY-1));

      wp=getwindowparams(hwnd);
      if(wp!=NULL && ((wp->childs)+1)!=NULL)
      {
        vbarlocation(wp->hwnd,(wp->childs+1)->hwnd,&(wp->vbar));
        hbarlocation(wp->hwnd,(wp->childs+1)->hwnd,&(wp->hbar));
        drawvbar(hdc,maxX,maxY,wp); /*REDRAW VBAR.*/
        drawhbar(hdc,maxX,maxY,wp); /*REDRAW HBAR.*/
      }
      ReleaseDC(hwnd,hdc);

      /*ExcludeClipRect(hdc,x1,y1,x2,y2);*/
      break;

    case WM_LBUTTONDOWN:
      point.x = LOWORD(lParam);
      point.y = HIWORD(lParam);

      wp=getwindowparams(hwnd);
      if(wp==NULL) break;

      getclientsize(hwnd,&maxX,&maxY); /*PARENT SHEET.*/
      getclientsize((wp->childs+0)->hwnd,&pw,&ph); /*BACKGROUND.*/
	  getwindowsize((wp->childs+1)->hwnd,&cw,&ch); /*CHILD.*/

      if(/*ch>ph &&*/ PtInRect(&(wp->vbar),point))
      {
		wp->sstatus=VSCROLLING; /*VERTICAL SCROLL BEGIN.*/
        pbar=point;
	  }
      else if(/*cw>pw &&*/ PtInRect(&(wp->hbar),point))
      {
        wp->sstatus=HSCROLLING; /*HORIZONTAL SCROLL BEGIN.*/
        pbar=point;
      }
      else
      {
        wp->sstatus=NEUTRAL; /*SCROLL END.*/
      }
	  break;

	/*Scroll with Wheel*/
	case WM_MOUSEWHEEL:
		wp = getwindowparams(hwnd);
		if (wp == NULL)
			break;

		point.x = LOWORD(lParam);
		point.y = HIWORD(lParam);

		getclientsize(hwnd, &maxX, &maxY);
		getwindowsize((wp->childs + 1)->hwnd, &cw, &ch); /* CHILD. */

		hflag=0;
		if (LOWORD(wParam) == MK_CONTROL)  hflag=1;

		rotation = int(HIWORD(wParam));

		cp.x = 0;
		cp.y = 0;
		ClientToScreen((wp->childs + 1)->hwnd, &cp);
		ScreenToClient(hwnd, &cp); /* TOPLEFT OF CHILD. */
		hdc = GetDC(hwnd);

		if(hflag==0)
		{
			bmax = maxY - BARWIDTH - 5;
			scrollunit=maxY/100;

			if (rotation <= 10000)
			{
				bt = wp->vbar.top   -scrollunit;
				bb = wp->vbar.bottom-scrollunit;
			}
			else if (rotation >= 10000)
			{
				bt = wp->vbar.top   +scrollunit;
				bb = wp->vbar.bottom+scrollunit;
			}

			if (cp.y > 1)
			{
				wp->vbar.top = 1;
				if (bb >= bmax)
					wp->vbar.bottom = bmax;
				else if (bb < bmax)
					wp->vbar.bottom = bb;
			}
			else if (ch + cp.y - 1 <= bmax)
			{
				wp->vbar.bottom = bmax;
				if (bt <= 1 || cp.y == 1)
					wp->vbar.top = 1;
				else if (bt > 1)
					wp->vbar.top = bt;
			}
			else
			{
				if (bt <= 1) {
					wp->vbar.top = 1;
					wp->vbar.bottom = 1 + bb - bt;
				}
				else if (bb >= bmax) {
					wp->vbar.bottom = bmax;
					wp->vbar.top = bmax - bb + bt;
				}
				else {
					wp->vbar.top = bt;
					wp->vbar.bottom = bb;
				}
			}
			drawvbar(hdc, maxX, maxY, wp); /* REDRAW VBAR. */
		}
		else
		{
			bmax = maxX - BARWIDTH - 5;
			scrollunit=maxX/10;

			if (rotation <= 10000)
			{
				bl = wp->hbar.left -scrollunit;
				br = wp->hbar.right-scrollunit;
			}
			else if (rotation >= 10000)
			{
				bl = wp->hbar.left +scrollunit;
				br = wp->hbar.right+scrollunit;
			}

			if (cp.x > 1) {
				wp->hbar.left = 1;
				if (br >= bmax)
					wp->hbar.right = bmax;
				else if (br < bmax)
					wp->hbar.right = br;
			}
			else if (cw + cp.x - 1 <= bmax) {
				wp->hbar.right = bmax;
				if (bl <= 1 || cp.x == 1)
					wp->hbar.left = 1;
				else if (bl > 1)
					wp->hbar.left = bl;
			}
			else {
				if (bl <= 1) {
					wp->hbar.left = 1;
					wp->hbar.right = 1 + br - bl;
				}
				else if (br >= bmax) {
					wp->hbar.right = bmax;
					wp->hbar.left = bmax - br + bl;
				}
				else {
					wp->hbar.left = bl;
					wp->hbar.right = br;
				}
			}
			drawhbar(hdc, maxX, maxY, wp); /* REDRAW HBAR. */
		}
		winleft = -(cw) * (wp->hbar.left - 1) / (maxX - BARWIDTH - 6);
		wintop = -(ch) * (wp->vbar.top - 1) / (maxY - BARWIDTH - 6);

		pbar = point;

		ReleaseDC(hwnd, hdc);
		MoveWindow((wp->childs + 1)->hwnd, winleft, wintop, cw, ch, TRUE);
		break;
	/*Scroll with Wheel*/

    case WM_MOUSEMOVE:
      wp=getwindowparams(hwnd);
      if(wp==NULL) break;

      if((wParam == MK_LBUTTON)&&
         (wp->sstatus==VSCROLLING || wp->sstatus==HSCROLLING))
      {
        point.x = LOWORD(lParam);
        point.y = HIWORD(lParam);

        getclientsize(hwnd,&maxX,&maxY);
        getwindowsize((wp->childs+1)->hwnd,&cw,&ch); /*CHILD.*/

        cp.x=0;
        cp.y=0;
        ClientToScreen((wp->childs+1)->hwnd,&cp);
        ScreenToClient(hwnd,&cp); /*TOPLEFT OF CHILD.*/

        hdc=GetDC(hwnd);

        if(wp->sstatus==VSCROLLING)
        {
		  bmax=maxY-BARWIDTH-5;

          bt=wp->vbar.top+point.y-pbar.y;
          bb=wp->vbar.bottom+point.y-pbar.y;

          if(cp.y>1)
          {
            wp->vbar.top=1;
            if(bb>=bmax)     wp->vbar.bottom=bmax;
            else if(bb<bmax) wp->vbar.bottom=bb;
          }
          else if(ch+cp.y-1<=bmax)
          {
            wp->vbar.bottom=bmax;
            if(bt<=1 || cp.y==1) wp->vbar.top=1;
            else if(bt>1)        wp->vbar.top=bt;
          }
          else
          {
            if(bt<=1)
            {
              wp->vbar.top   =1;
              wp->vbar.bottom=1+bb-bt;
            }
            else if(bb>=bmax)
            {
              wp->vbar.bottom=bmax;
              wp->vbar.top   =bmax-bb+bt;
            }
            else
            {
              wp->vbar.top   =bt;
              wp->vbar.bottom=bb;
            }
          }
          pbar=point;

          drawvbar(hdc,maxX,maxY,wp); /*REDRAW VBAR.*/
        }
        else if(wp->sstatus==HSCROLLING)
        {
          bmax=maxX-BARWIDTH-5;

          bl=wp->hbar.left +point.x-pbar.x;
          br=wp->hbar.right+point.x-pbar.x;

          if(cp.x>1)
          {
            wp->hbar.left=1;
            if(br>=bmax)     wp->hbar.right=bmax;
            else if(br<bmax) wp->hbar.right=br;
          }
          else if(cw+cp.x-1<=bmax)
          {
            wp->hbar.right=bmax;
            if(bl<=1 || cp.x==1) wp->hbar.left=1;
            else if(bl>1)        wp->hbar.left=bl;
          }
          else
          {
            if(bl<=1)
            {
              wp->hbar.left =1;
              wp->hbar.right=1+br-bl;
            }
            else if(br>=bmax)
            {
              wp->hbar.right=bmax;
              wp->hbar.left =bmax-br+bl;
            }
            else
            {
              wp->hbar.left =bl;
              wp->hbar.right=br;
            }
          }
          pbar=point;

          drawhbar(hdc,maxX,maxY,wp); /*REDRAW HBAR.*/
        }
        ReleaseDC(hwnd,hdc);

        winleft=-(cw)*(wp->hbar.left-1)/(maxX-BARWIDTH-6);
        wintop =-(ch)*(wp->vbar.top -1)/(maxY-BARWIDTH-6);

        MoveWindow((wp->childs+1)->hwnd,winleft,wintop,cw,ch,TRUE);
      }
      else
      {
        wp->sstatus=NEUTRAL; /*SCROLL END*/
      }
      break;

    default:
	  return DefWindowProc(hwnd,message,wParam,lParam);
  }
  return 0;
}/*WindowProcedureSheet*/

LRESULT CALLBACK WindowProcedureBack(HWND hwnd,
                                     UINT message,
                                     WPARAM wParam,
                                     LPARAM lParam)
/*DISPLAY AREA OF SHEET WINDOW.*/
{
  HWND hoya;
  long int maxX,maxY;

  switch(message)
  {
    case WM_COMMAND:
      switch(LOWORD(wParam))
      {
        case IDM_FITPARENT: /*FIT SIZE TO PARENT.*/
          hoya=GetParent(hwnd);
          getclientsize(hoya,&maxX,&maxY);

          MoveWindow(hwnd,1,1,
                          (maxX-BARWIDTH-5),(maxY-BARWIDTH-5),TRUE);
          break;

        default: /*OTHERS*/
          return DefWindowProc(hwnd,message,wParam,lParam);
      }
      break;

	default:
      return DefWindowProc(hwnd,message,wParam,lParam);
  }
  return 0;
}/*WindowProcedureBack*/

LRESULT CALLBACK WindowProcedureDraw(HWND hwnd,
									 UINT message,
									 WPARAM wParam,
									 LPARAM lParam)
/*WINDOW FOR DRAWING WITH POPUP MENU.*/
{
  HWND hoya;
  /*HWND hdlg;*/
  HDC hdc;
  HCURSOR hcursor;
  POINT point;
  WPARAM wparam;
  /*MSG msg;*/
  /*HACCEL haccel;*/
  size_t n;
  char str[256],fname[256];
  int i,ii,j,k,direction,onplane,find,ndiv;
  int Nx,Ny;
  long int maxX,maxY,x,y,mw,mh;
  long int code,joff,koff,loff,moff,soff,nelem,pnelem,nnode,*nflag;
  double x0,y0,z0/*Zoom with Wheel*/,x1,y1,z1,x2,y2,z2;
  //double x1,y1,z1,x2,y2,z2;
  double length,dx,dy,dz,distance,radius,rrate;
  double eps=1.0E-4;
  struct onode *node,nearest,cnode,tn;
  struct owire *elem;
  struct oelem *foundelem,*pelems;
  struct obans *selectedban;
  struct oconf *oc;
  struct line  *selectedline;
  struct plane pban;

  int laps;
  double dsafety;


  /*Zoom with Wheel*/
  double scrollunit;
  int rotation;
  int nwd;
  int fflag,flag,ndiveide;
  double len;
  /*Zoom with Wheel*/

  struct viewparam vp;
  FILE *fin;
  double *srates;
  double venedisp[5][12],initz[12];

  switch(message)
  {
	case WM_PAINT:
      DefWindowProc(hwnd,message,wParam,lParam);

      /*GET TEXTS IN ALL ITEMS OF OPTION DIALOG.....NOT YET.*/

      overlayhdc(*(wdraw.childs+1),SRCPAINT);

      hdc = GetDC(hwnd);
	  setfontformat(hdc,15,6,"MS Mincho",0,255,255);
      drawtexts(hdc,wdraw.strset,wdraw.nstring,
				(wdraw.childs+1)->vparam);
	  break;

	case WM_LBUTTONDOWN:
	  x = LOWORD(lParam);
	  y = HIWORD(lParam);
	  point.x=x;
	  point.y=y;

	  prestatus=NEUTRAL;

	  /*WHILE CREATING FRAME.*/
	  if((wmenu.childs+2)->vparam.vflag.mv.ftype==F_ORGAN)
	  {
		if(globalstatus==SECTIONPAINT ||
		   globalstatus==COPYBOND)
		{
		  foundelem=selectorganelement((wdraw.childs+1)->vparam,
									   &((wdraw.childs+1)->org),
									   point,
									   &selectedline,
									   &selectedban);
		  if(foundelem==NULL) break;

		  /*DRAW NODE OF SELECTED ELEM.*/
		  /*hdc = GetDC((wdraw.childs+1)->hwnd);*/
		  hdc = (wdraw.childs+1)->hdcC;
		  setfontformat(hdc,15,6,"MS Mincho",255,255,255);
		  for(i=0;i<(foundelem->nnod);i++)
		  {
			nodeontoscreen(**(foundelem->nods+i),&Nx,&Ny,
						   (wdraw.childs+1)->vparam);
			sprintf(str,"%ld",(*(foundelem->nods+i))->code);
			TextOut(hdc,Nx+2,Ny,str,strlen(str));
		  }
		  ReleaseDC((wdraw.childs+1)->hwnd,hdc);
		  SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);

		  if(globalstatus==SECTIONPAINT)
		  {
			foundelem->sect=gelem.sect;
			foundelem->type=gelem.type;
			foundelem->role=gelem.role;
		  }
		  if(globalstatus==COPYBOND)
		  {
			for(i=0;i<(foundelem->nnod);i++)
			{
			  for(j=0;j<6;j++)
			  {
				*(foundelem->bonds+6*i+j)=*(gelem.bonds+6*i+j);
			  }
			}
		  }

		  clearwindow(*(wdraw.childs+1));
		  draworganization((wdraw.childs+1)->hdcC,
						   (wdraw.childs+1)->vparam,
						   (wdraw.childs+1)->org,ONSCREEN);
		  SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);
		  break;
		}
		else if(globalstatus==DIVIDEELEMENT)
		{
		  foundelem=selectorganelement((wdraw.childs+1)->vparam,
									   &((wdraw.childs+1)->org),
									   point,
									   &selectedline,
									   &selectedban);
		  if(foundelem==NULL) break;

		  /*DRAW NODE OF SELECTED ELEM.*/
		  /*hdc = GetDC((wdraw.childs+1)->hwnd);*/
		  hdc = (wdraw.childs+1)->hdcC;
		  setfontformat(hdc,15,6,"MS Mincho",255,255,255);
		  for(i=0;i<(foundelem->nnod);i++)
		  {
			nodeontoscreen(**(foundelem->nods+i),&Nx,&Ny,
						   (wdraw.childs+1)->vparam);
			sprintf(str,"%ld",(*(foundelem->nods+i))->code);
			TextOut(hdc,Nx+2,Ny,str,strlen(str));
		  }
		  ReleaseDC((wdraw.childs+1)->hwnd,hdc);
		  SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);

          //ndiv=2;
          getlaps((wmenu.childs+2)->hwnd,&ndiv);//ujioka:分割数をLAPSで指定
		  if(ndiv==1) ndiv=2;

		  x1=(*(foundelem->nods+0))->d[0];
		  y1=(*(foundelem->nods+0))->d[1];
		  z1=(*(foundelem->nods+0))->d[2];
		  x2=(*(foundelem->nods+1))->d[0];
		  y2=(*(foundelem->nods+1))->d[1];
		  z2=(*(foundelem->nods+1))->d[2];
		  dx=x2-x1;
		  dy=y2-y1;
		  dz=z2-z1;
		  /*length = sqrt(dx*dx+dy*dy+dz*dz);*/

/*
sprintf(str,"x1=%.3f x2=%.3f dx=%.3f",x1,x2,dx);
MessageBox(NULL,str,"Divide Element",MB_OK);
*/

		  for(ii=0;ii<(ndiv-1);ii++)
		  {
			nnode=(wdraw.childs+1)->org.nnode;

			gelem.code=((wdraw.childs+1)->org.elems
						+(wdraw.childs+1)->org.nelem-1)->code+1;
			gelem.loff=(wdraw.childs+1)->org.nelem; /*LAST.*/
			gelem.sect=foundelem->sect;
			gelem.type=foundelem->type;

			gelem.nnod=2;
			gelem.nban=0;
			/*gelem.nnod--;*/

			/*gelem.nods=(struct onode **)
					   realloc(gelem.nods,
							   gelem.nnod
							   *sizeof(struct onode *));*/
			gelem.nods=(struct onode **)
					   malloc(gelem.nnod
							  *sizeof(struct onode *));
			gelem.bonds=(signed char *)
						malloc(6*(gelem.nnod)
							   *sizeof(signed char));
			if(gelem.nods==NULL || gelem.bonds==NULL)
			{
			  MessageBox(NULL,"Buffer Null.","Adding",MB_OK);
			}

			if(ii==ndiv-2)
			{
//			  for(i=0;i<6*(gelem.nnod);i++)
			  for(i=0;i<6;i++) *(gelem.bonds+i)=0;

			  for(i=6;i<12;i++)//ujioka:終端の境界条件保存
			  {
				*(gelem.bonds+i)=*(foundelem->bonds+i);
			  }
			}
			else
			{
			  for(i=0;i<6*(gelem.nnod);i++) *(gelem.bonds+i)=0;
			}
/*
sprintf(str,"Passed step %d",ii);
MessageBox(NULL,str,"Divide Element",MB_OK);
*/

			if(ii==0)
			{
			  cnode.code=((wdraw.childs+1)->org.nodes+nnode-1)
						  ->code+1;
			  cnode.loff=nnode;

			  cnode.d[0]=x1+dx*((double)ii+1)/ndiv;
			  cnode.d[1]=y1+dy*((double)ii+1)/ndiv;
			  cnode.d[2]=z1+dz*((double)ii+1)/ndiv;

/*
sprintf(str,"x1=%.3f x2=%.3f dx=%.3f rate=%.3f cx=%.3f",x1,x2,dx,((double)ii+1)/ndiv,cnode.d[0]);
MessageBox(NULL,str,"Divide Element",MB_OK);
*/

		   	  *(gelem.nods+0)=addnode(cnode,NULL,&((wdraw.childs+1)->org));
			  pnode=*(gelem.nods+0);
			}
			else *(gelem.nods+0)=pnode;

			nnode=(wdraw.childs+1)->org.nnode;

			if(ii<ndiv-2)
			{
			  cnode.code=((wdraw.childs+1)->org.nodes+nnode-1)
						  ->code+1;
			  cnode.loff=nnode;

			  cnode.d[0]=x1+dx*((double)ii+2)/ndiv;
			  cnode.d[1]=y1+dy*((double)ii+2)/ndiv;
			  cnode.d[2]=z1+dz*((double)ii+2)/ndiv;

			  *(gelem.nods+1)=addnode(cnode,NULL,&((wdraw.childs+1)->org));
			  pnode=*(gelem.nods+1);
			}
			else if(ii==ndiv-2) *(gelem.nods+1)=*(foundelem->nods+1);

			addelement(&gelem,&((wdraw.childs+1)->org));

			/*freeelement(&gelem,0,0,1,0);*/

/*MessageBox(NULL,"Passed.","Adding Element",MB_OK);*/
			/*gelem.nnod=0;*/
			/*gelem.nban=0;*/
			gelem.nods=NULL;
			gelem.bonds=NULL;
			gelem.bans=NULL;
		  }

		  *(foundelem->nods+1)=(wdraw.childs+1)->org.nodes+nnode-ndiv+1;
//		  for(i=0;i<6*(foundelem->nnod);i++) *(foundelem->bonds+i)=0;
		  for(i=6;i<12;i++) *(foundelem->bonds+i)=0;//ujioka：始端の境界条件保存

		  clearwindow(*(wdraw.childs+1));
		  draworganization((wdraw.childs+1)->hdcC,
						   (wdraw.childs+1)->vparam,
						   (wdraw.childs+1)->org,ONSCREEN);
		  SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);
		  break;
		}
		else if(globalstatus==MOVENODEBYMOUSE ||
				globalstatus==PINCHNODEBYMOUSE ||
				globalstatus==PINCHNODEWITHSAFETY ||
				globalstatus==PINCHNODEVENEZIA ||
				globalstatus==PINCHNODESTACK ||
				globalstatus==PINCHNODEELASTIC ||
				globalstatus==MOVENODEWITHSAFETY)
		{
		  initx=x;
		  inity=y;

		  node=selectorgannode((wdraw.childs+1)->vparam,
							   &((wdraw.childs+1)->org),point);
		  if(node==NULL)
		  {
			MessageBox(NULL,"Nothing.","MoveNode",MB_OK);
			break;
		  }
		  gnode=node;

		  pnode=(struct onode *)malloc(sizeof(struct onode));
		  *pnode=*gnode;

		  gincrement.dc[0]=0.0;
		  gincrement.dc[1]=0.0;
		  gincrement.dc[2]=0.0;

		  if(globalstatus==MOVENODEBYMOUSE)
		  {
			globalstatus=MOVENODEBYMOUSEBEGIN;
			prestatus=MOVENODEBYMOUSE;
		  }
		  if(globalstatus==PINCHNODEBYMOUSE)
		  {
			globalstatus=PINCHNODEBYMOUSEBEGIN;
			prestatus=PINCHNODEBYMOUSE;
		  }
		  if(globalstatus==PINCHNODEWITHSAFETY)
		  {
			globalstatus=PINCHNODEWITHSAFETYBEGIN;
			prestatus=PINCHNODEWITHSAFETY;
		  }
		  if(globalstatus==PINCHNODEVENEZIA)
		  {
			globalstatus=PINCHNODEVENEZIABEGIN;
			prestatus=PINCHNODEVENEZIA;
		  }
		  if(globalstatus==PINCHNODESTACK)
		  {
			globalstatus=PINCHNODESTACKBEGIN;
			prestatus=PINCHNODESTACK;
		  }
		  if(globalstatus==PINCHNODEELASTIC)
		  {
			globalstatus=PINCHNODEELASTICBEGIN;
			prestatus=PINCHNODEELASTIC;
		  }
		  if(globalstatus==MOVENODEWITHSAFETY)
		  {
			globalstatus=MOVENODEWITHSAFETYBEGIN;
			prestatus=MOVENODEWITHSAFETY;
		  }
		  break;
		}
		else if(globalstatus==MOVEELEMENTBYMOUSE ||
				globalstatus==COPYELEMENTBYMOUSE ||
				globalstatus==PINCHELEMTWIST)
		{
		  initx=x;
		  inity=y;

		  foundelem=selectorganelement((wdraw.childs+1)->vparam,
									   &((wdraw.childs+1)->org),
									   point,
									   &selectedline,
									   &selectedban);
		  if(foundelem==NULL)
		  {
			MessageBox(NULL,"Nothing.","Found",MB_OK);
			break;
		  }

		  /*DRAW NODE OF SELECTED ELEM.*/
		  /*hdc = GetDC((wdraw.childs+1)->hwnd);*/
		  hdc = (wdraw.childs+1)->hdcC;
		  setfontformat(hdc,15,6,"MS Mincho",255,255,255);
		  for(i=0;i<(foundelem->nnod);i++)
		  {
			nodeontoscreen(**(foundelem->nods+i),&Nx,&Ny,
						   (wdraw.childs+1)->vparam);
			sprintf(str,"%ld",(*(foundelem->nods+i))->code);
			TextOut(hdc,Nx+2,Ny,str,strlen(str));
		  }
		  ReleaseDC((wdraw.childs+1)->hwnd,hdc);
		  SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);

		  pelem=foundelem;
		  pnode=(struct onode *)malloc(pelem->nnod
									   *sizeof(struct onode));
		  for(i=0;i<pelem->nnod;i++)
		  {
			*(pnode+i)=*(*(pelem->nods+i));
		  }

		  pcolor.r=pelem->color.line.r;
		  pcolor.g=pelem->color.line.g;
		  pcolor.b=pelem->color.line.b;
		  pelem->color.line.r=150;
		  pelem->color.line.g=150;
		  pelem->color.line.b=150;

		  gincrement.dc[0]=0.0;
		  gincrement.dc[1]=0.0;
		  gincrement.dc[2]=0.0;

		  if(globalstatus==MOVEELEMENTBYMOUSE)
		  {
			globalstatus=MOVEELEMENTBYMOUSEBEGIN;
			prestatus   =MOVEELEMENTBYMOUSE;
		  }
		  else if(globalstatus==PINCHELEMTWIST)
		  {
			globalstatus=PINCHELEMTWISTBEGIN;
			prestatus=PINCHELEMTWIST;
		  }
		  else if(globalstatus==COPYELEMENTBYMOUSE)
		  {
			globalstatus=COPYELEMENTBYMOUSEBEGIN;
			prestatus   =COPYELEMENTBYMOUSE;

			gelem=(*pelem);
			melem=(*pelem);

			melem.nods=(struct onode **)
						malloc((pelem->nnod)*sizeof(struct onode *));
			melem.bonds=(signed char *)
						malloc(6*(pelem->nnod)*sizeof(signed char));
			melem.bans=(struct obans *)
					   malloc((pelem->nban)*sizeof(struct obans));
			for(i=0;i<(pelem->nban);i++)
			{
			  *(melem.bans+i)=*(pelem->bans+i); /*COPY BANS*/
			  (melem.bans+i)->nods=(struct onode **)
								   malloc((pelem->bans+i)->nnod
								   *sizeof(struct onode *));
			}

			for(i=0;i<pelem->nnod;i++)
			{
			  /*COPY NODES*/
			  *(melem.nods+i)=(struct onode *)
							  malloc(sizeof(struct onode));
			  **(melem.nods+i)=**(pelem->nods+i);

			  for(j=0;j<6;j++) /*COPY BONDS*/
			  {
				*(melem.bonds+6*i+j)=*(pelem->bonds+6*i+j);
			  }
			}

			for(i=0;i<(pelem->nban);i++) /*CORRECT POINTER*/
			{
			  for(j=0;j<(pelem->bans+i)->nnod;j++)
			  {
				for(k=0;k<(pelem->nnod);k++)
				{
				  if(*((pelem->bans+i)->nods+j)==*(pelem->nods+k))
				  {
                    *((melem.bans+i)->nods+j)=*(melem.nods+k);
				  }
                }
              }
            }

            pelem=&melem;
		  }

          break;
        }
		else if(globalstatus==MOVENODENODETONODE ||
                globalstatus==MOVEELEMENTNODETONODE ||
                globalstatus==COPYNODENODETONODE ||
				globalstatus==COPYELEMENTNODETONODE)
		{
          node=selectorgannode((wdraw.childs+1)->vparam,
                               &((wdraw.childs+1)->org),point);
          if(node==NULL) break;

		  icount++;
          if(icount==1)
          {
            gincrement.dc[0]=-(node->d[0]);
			gincrement.dc[1]=-(node->d[1]);
			gincrement.dc[2]=-(node->d[2]);

            nodeontoscreen(*node,&Nx,&Ny,
                           (wdraw.childs+1)->vparam); /*PROJECTION*/
            hdc = GetDC((wdraw.childs+1)->hwnd);
            setfontformat(hdc,15,6,"MS Mincho",255,255,255);
            sprintf(str,"%ld",node->code);
			TextOut(hdc,Nx+2,Ny,str,strlen(str));
            ReleaseDC((wdraw.childs+1)->hwnd,hdc);
          }
          else if(icount>=2)
          {
            gincrement.dc[0]+=node->d[0];
			gincrement.dc[1]+=node->d[1];
            gincrement.dc[2]+=node->d[2];

			if(globalstatus==MOVENODENODETONODE)
			{
			  pnode->d[0]+=gincrement.dc[0];
			  pnode->d[1]+=gincrement.dc[1];
			  pnode->d[2]+=gincrement.dc[2];
			}
			else if(globalstatus==MOVEELEMENTNODETONODE)
			{
			  for(i=0;i<pelem->nnod;i++)
			  {
				(*(pelem->nods+i))->d[0]+=gincrement.dc[0];
				(*(pelem->nods+i))->d[1]+=gincrement.dc[1];
				(*(pelem->nods+i))->d[2]+=gincrement.dc[2];
			  }
			}
			else if(globalstatus==COPYNODENODETONODE)
			{
			  nnode=(wdraw.childs+1)->org.nnode;
			  loff=pnode->loff;

			  oc=(struct oconf *)malloc(6*sizeof(struct oconf));
			  for(j=0;j<6;j++) /*COPY CONFS*/
			  {
				*(oc+j)=*((wdraw.childs+1)->org.confs+6*loff+j);
			  }

			  cnode=*pnode;
			  cnode.code=((wdraw.childs+1)->org.nodes+nnode-1)
						 ->code+1;
			  cnode.loff=nnode;

			  cnode.d[0]+=gincrement.dc[0];
			  cnode.d[1]+=gincrement.dc[1];
			  cnode.d[2]+=gincrement.dc[2];

			  find=0;
			  for(j=0;j<nnode;j++) /*CANCEL IF OVERLAPPED*/
			  {
				if(cnode.d[0]
				   >((wdraw.childs+1)->org.nodes+j)->d[0]-0.001 &&
				   cnode.d[0]
				   <((wdraw.childs+1)->org.nodes+j)->d[0]+0.001 &&
				   cnode.d[1]
				   >((wdraw.childs+1)->org.nodes+j)->d[1]-0.001 &&
				   cnode.d[1]
				   <((wdraw.childs+1)->org.nodes+j)->d[1]+0.001 &&
				   cnode.d[2]
				   >((wdraw.childs+1)->org.nodes+j)->d[2]-0.001 &&
				   cnode.d[2]
				   <((wdraw.childs+1)->org.nodes+j)->d[2]+0.001)
				{
				  find=1;
				  break;
				}
			  }
			  if(!find) /*ADD NODE*/
			  {
				addnode(cnode,oc,&((wdraw.childs+1)->org));
			  }
			}
			else if(globalstatus==COPYELEMENTNODETONODE)
			{
			  addelementwithnode(pelem,
								 &((wdraw.childs+1)->org),
								 gincrement.dc[0],
								 gincrement.dc[1],
								 gincrement.dc[2]);
								 /*UNDER CONSTRUCTION*/
			}

			clearwindow(*(wdraw.childs+1));
			draworganization((wdraw.childs+1)->hdcC,
							 (wdraw.childs+1)->vparam,
							 (wdraw.childs+1)->org,ONSCREEN);
			SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);

			hcursor = LoadCursor(hInstGlobal,"CANCURSORW");
			SetClassLong((wdraw.childs+1)->hwnd,
						 GCL_HCURSOR,(LONG)hcursor);
			icount=0;
			globalstatus=NEUTRAL;
		  }
		  break;
		}
		else if(globalstatus==COPYCONF)
		{
		  node=selectorgannode((wdraw.childs+1)->vparam,
							   &((wdraw.childs+1)->org),point);
		  if(node==NULL) break;

		  for(i=0;i<6;i++)
		  {
			((wdraw.childs+1)->org.confs+(node->loff)*6+i)->iconf
			=((wdraw.childs+1)->org.confs+(pnode->loff)*6+i)->iconf;
		  }

		  clearwindow(*(wdraw.childs+1));
		  draworganization((wdraw.childs+1)->hdcC,
						   (wdraw.childs+1)->vparam,
						   (wdraw.childs+1)->org,ONSCREEN);
		  SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);

		  break;
		}
		else if(globalstatus==CHANGEBOND)
		{
		  loff=selectelemnode((wdraw.childs+1)->vparam,
							  pelem,point);
		  if(loff==0) break;

		  loff-=1;
		  pnode=*(pelem->nods+loff);

		  hpopdlg=CreateDialog(hInstGlobal,
							   "HOGDLGCHANGEBOND",
							   NULL,
							   (DLGPROC)DialogProcBond);
		  getwindowsize(hpopdlg,&mw,&mh);
		  ClientToScreen((wdraw.childs+1)->hwnd,&point);
		  MoveWindow(hpopdlg,point.x,point.y,mw,mh,TRUE);
		  ShowWindow(hpopdlg,SW_SHOW);

		  break;
		}
		else if(globalstatus==SELECTMULTIELEM)
		{
		  initx=x;
		  inity=y;

		  globalstatus=SELECTMULTIELEMBEGIN;
		  prestatus   =SELECTMULTIELEM;
		  break;
		}
		else if(globalstatus==SELECTMULTINODE)
		{
		  initx=x;
		  inity=y;

		  globalstatus=SELECTMULTINODEBEGIN;
		  prestatus   =SELECTMULTINODE;
		  break;
		}
		else if(createcommand==C_ADD ||
				globalstatus==ADDELEMWITHSAFETY)
		{
		  if(createitem==C_NODE)
		  {
			createorgannode(0,(wdraw.childs+1)->vparam,
							x,y,&((wdraw.childs+1)->org));
		  }
		  if(createitem==C_ELEMENT)
		  {
			if(globalstatus==SELECTNODE ||
			   globalstatus==ADDELEMWITHSAFETY)
			{
/*MessageBox(NULL,"Pass 1","Add Elem",MB_OK);*/
			  node=selectorgannode((wdraw.childs+1)->vparam,
								   &((wdraw.childs+1)->org),point);
			  if(node==NULL)
			  {
				/*MessageBox(NULL,"Nothing.","Node",MB_OK);*/
				break;
			  }
			  /*prestatus=SELECTNODE;*/
			}
			else if(globalstatus==SELECTPERPENDICULAR &&
					gelem.nnod>1)
			{
			  foundelem=selectorganelement((wdraw.childs+1)->vparam,
										   &((wdraw.childs+1)->org),
										   point,
										   &selectedline,
										   &selectedban);
			  if(foundelem==NULL &&
				 (*(gelem.nods+gelem.nnod-2))->d[GZ]==0.0) break;

			  if(selectedline!=NULL) /*PERPENDICULAR WITH LINE.*/
			  {
				distancedotline(*(*(gelem.nods+gelem.nnod-2)),
								*selectedline,
								&nearest,NULL);
				node=&nearest;
			  }
			  else if(selectedban!=NULL) /*PERPENDICULAR WITH BAN.*/
			  {
				bantoplane(*selectedban,&pban);

				distancedotplane(*(*(gelem.nods+gelem.nnod-2)),
								 pban,
								 &nearest);
				node=&nearest;
			  }
			  else /*CREATE ON GROUND PLANE:Z=0*/
			  {
				pban.nods[0].d[GX]=0.0;
				pban.nods[0].d[GY]=0.0;
				pban.nods[0].d[GZ]=0.0;

				pban.nvec.dc[GX]=0.0;
				pban.nvec.dc[GY]=0.0;
				pban.nvec.dc[GZ]=1.0;

				pban.a=0.0; pban.b=0.0; pban.c=1.0; pban.d=0.0;

				distance
				=distancedotplane(*(*(gelem.nods+gelem.nnod-2)),
								  pban,
								  &nearest);
				node=&nearest; /*PERPENDICULAR WITH GROUND.*/
			  }
			}

			if(gelem.nnod==0)
			{
			  gelem.nnod=1;
			  gelem.nods=(struct onode **)
						 malloc(gelem.nnod
								*sizeof(struct onode *));
/*gelem.type=SLAB;*/
gelem.type=TYPENULL;
			}
			else if(gelem.nnod==2)
			{
			  gelem.nban=1;
			  gelem.bans=(struct obans *)
						 malloc(sizeof(struct obans));
			  (gelem.bans+0)->code=1;
			  (gelem.bans+0)->loff=0;
			  (gelem.bans+0)->nnod=2;
			  (gelem.bans+0)->nods=(struct onode **)
								   malloc(2*sizeof(struct onode *));
			  if(gelem.bans==NULL || (gelem.bans+0)->nods==NULL)
			  {
				MessageBox(NULL,"Buffer Null.","Adding",MB_OK);
				break;
			  }

			  *((gelem.bans+0)->nods+0)=*(gelem.nods+0);
			}

			if(gelem.nnod<=3)
			{
			  if(globalstatus==SELECTNODE ||
				 globalstatus==ADDELEMWITHSAFETY)
			  {
				nlast=*node;
			  }
			  else if(globalstatus==SELECTPERPENDICULAR &&
					  gelem.nnod>0)
			  {
				nlast=*node;
			  }
			  else
			  {
				nlast=findlastnode(0,(wdraw.childs+1)->vparam,
								   x,y,
								   NULL,
								   &((wdraw.childs+1)->org));
				if(nlast.code==0)
				{
				  MessageBox(NULL,"Find Node Failed.","Adding",
							 MB_OK);
				  break;
				}
			  }
			}
			else if(gelem.nban>0)
			{
			  bantoplane(*(gelem.bans+0),&pban);

			  onplane=0;
			  if(globalstatus==SELECTNODE ||
				 globalstatus==ADDELEMWITHSAFETY)
			  {
				distance=distancedotplane(*node,pban,&nearest);

				if(fabs(distance)<=eps) /*NODE ON PLANE.*/
				{
				  nlast=*node;
				  onplane=1;
				}
				else /*NOT ON PLANE.CHANGE X,Y AS DOT ON SCREEN.*/
				{
				  /*MessageBox(NULL,"Not Node.","Adding",
							 MB_OK);*/ /*CANCEL:UNDER CONSTRUCTION.*/

				  nlast=*node;
				  onplane=1;

				  /*
				  if(!nodeontoscreen(*node,&Nx,&Ny,
									 (wdraw.childs+1)->vparam))
				  {
					MessageBox(NULL,"Projection Failed.","Adding",
							   MB_OK);
				  }
				  x=Nx;
				  y=Ny;
				  */
				}
			  }
			  else if(globalstatus==SELECTPERPENDICULAR &&
					  gelem.nnod>0)
			  {
				distance=distancedotplane(*node,pban,&nearest);

				if(distance<=eps) /*PERPENDICULAR DOT ON PLANE.*/
				{
				  nlast=*node;
				  onplane=1;
				}
				else /*NOT ON PLANE.CHANGE X,Y AS DOT ON SCREEN.*/
				{
				  MessageBox(NULL,"Not On Plane.","Adding",
							 MB_OK); /*CANCEL:UNDER CONSTRUCTION.*/
				  if(!nodeontoscreen(*node,&Nx,&Ny,
									 (wdraw.childs+1)->vparam))
				  {
					MessageBox(NULL,"Projection Failed.","Adding",
							   MB_OK);
				  }
				  x=Nx;
				  y=Ny;
				}
			  }

			  if(!onplane)
			  {
				nlast=findlastnode(0,(wdraw.childs+1)->vparam,
								   x,y,
								   &pban,
								   &((wdraw.childs+1)->org));
				if(nlast.code==0)
				{
				  MessageBox(NULL,"Find Node Failed.","Adding",
							 MB_OK);
				  break;
				}
			  }
			}

			gelem.nnod++;
			gelem.nods=(struct onode **)
					   realloc(gelem.nods,
							   gelem.nnod
							   *sizeof(struct onode *));
			if(gelem.nods==NULL)
			{
			  MessageBox(NULL,"Buffer Null.","Adding",MB_OK);
			  break;
			}

			if((globalstatus==SELECTNODE && onplane) ||
			   (globalstatus==ADDELEMWITHSAFETY && onplane))
			{
			  *(gelem.nods+gelem.nnod-2)=node;
			}
			else if(globalstatus==SELECTPERPENDICULAR &&
					gelem.nnod>0 &&
					onplane)
			{
			  code=0;
			  for(i=0;i<((wdraw.childs+1)->org.nnode);i++)
			  {
				if(code<((wdraw.childs+1)->org.nodes+i)->code)
				{
				  code=((wdraw.childs+1)->org.nodes+i)->code;
				}
			  }
			  nlast.code=code+1;
			  nlast.loff=(wdraw.childs+1)->org.nnode;

			  *(gelem.nods+gelem.nnod-2)
			  =addnode(nlast,NULL,&((wdraw.childs+1)->org));
			}
			else
			{
			  *(gelem.nods+gelem.nnod-2)
			  =addnode(nlast,NULL,&((wdraw.childs+1)->org));
			}

			/**(gelem.nods+gelem.nnod-2)
			=(wdraw.childs+1)->org.nodes+(nlast.loff);*/
			*(gelem.nods+gelem.nnod-1)=&nlast;

			if(gelem.nban>0)
			{
			  (gelem.bans+0)->nnod++;
			  (gelem.bans+0)->nods=(struct onode **)
								   realloc((gelem.bans+0)->nods,
										   (gelem.bans+0)->nnod
										   *sizeof(struct onode *));
			  if((gelem.bans+0)->nods==NULL)
			  {
				MessageBox(NULL,"Buffer Null.","Adding",MB_OK);
				break;
			  }

			  *((gelem.bans+0)->nods+(gelem.bans+0)->nnod-2)
			  =*(gelem.nods+gelem.nnod-2);
			  *((gelem.bans+0)->nods+(gelem.bans+0)->nnod-1)
			  =*(gelem.nods+gelem.nnod-1);
			}

			/*globalstatus=NEUTRAL;
			hcursor = LoadCursor(hInstGlobal,"CANCURSORW");
			SetClassLong((wdraw.childs+1)->hwnd,
						 GCL_HCURSOR,(LONG)hcursor);*/
/*MessageBox(NULL,"Pass 2","Add Elem",MB_OK);*/
		  }

		  draworganization((wdraw.childs+1)->hdcC,
						   (wdraw.childs+1)->vparam,
						   (wdraw.childs+1)->org,ONSCREEN);

          SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);
		  break;
		}
	  }

      if((wmenu.childs+2)->vparam.vflag.mv.ftype==F_ARCLM)
      /*UJIOKA FOR BCLNG003*/
      {
        if(globalstatus==SELECTMULTIELEM)
        {
		  initx=x;
		  inity=y;

		  globalstatus=SELECTMULTIELEMBEGIN;
		  prestatus   =SELECTMULTIELEM;
          break;
        }

      }

	  /*CASE NOT CREATING.*/
	  if((wmenu.childs+2)->vparam.vflag.mv.ftype==F_SRCAN)
	  {
		initx=x;
		inity=y;

        globalstatus=MOVE;
      }
      else if(globalstatus==SELECTNODE)
	  {
        node=selectnode((wdraw.childs+1)->vparam,
                        &arc,point);
		if(node!=NULL)
        {
          setnode((wmenu.childs+3)->hwnd,node,arc.confs);

          clearwindow(*(wdraw.childs+1));
          drawarclmframe((wdraw.childs+1)->hdcC,
                         (wdraw.childs+1)->vparam,arc,
                         node->code,ONSCREEN);
          SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);
          break;
        }
        prestatus=SELECTNODE;
      }
      else if(globalstatus==SELECTELEMENT)
	  {
		elem=selectelement((wdraw.childs+1)->vparam,
                           &arc,point);
        if(elem!=NULL)
        {
          setelement((wmenu.childs+3)->hwnd,elem);
          setnode((wmenu.childs+3)->hwnd,elem->node[0],arc.confs);
          setsection((wmenu.childs+3)->hwnd,elem->sect);

          clearwindow(*(wdraw.childs+1));
		  drawarclmframe((wdraw.childs+1)->hdcC,
                         (wdraw.childs+1)->vparam,arc,
                         elem->code,ONSCREEN);
          SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);
          break;
        }
		prestatus=SELECTELEMENT;
      }
	  else if(globalstatus==SELECTSECTION)
      {
        elem=selectelement((wdraw.childs+1)->vparam,
                           &arc,point);
        if(elem!=NULL)
        {
          setsection((wmenu.childs+3)->hwnd,elem->sect);

          clearwindow(*(wdraw.childs+1));
          drawarclmframe((wdraw.childs+1)->hdcC,
						 (wdraw.childs+1)->vparam,arc,
                         elem->sect->code,ONSCREEN);
          SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);
          break;
        }
        prestatus=SELECTSECTION;
      }
      else if((wdraw.childs+1)->lstatus==NEUTRAL)
      {
        sprintf(str,"POSITION:%d %d",x,y);
		SetTextColor((wdraw.childs+1)->hdcC,RGB(0,0,255)); /*COLOR.*/
        TextOut((wdraw.childs+1)->hdcC,x,y,str,strlen(str));

        SendMessage(hwnd,WM_PAINT,0,0);
		break;
      }

	  if((wdraw.childs+1)->lstatus==ROTATE)
      {
        initx=x;
        inity=y;

        initv=(wdraw.childs+1)->vparam;

		globalstatus=ROTATE;
      }
      else if((wdraw.childs+1)->lstatus==MOVE)
      {
        initx=x;
        inity=y;

        initv=(wdraw.childs+1)->vparam;

        globalstatus=MOVE;
	  }

      break;

    case WM_RBUTTONDOWN:
      point.x = LOWORD(lParam);
      point.y = HIWORD(lParam);
	  initx = point.x;
	  inity = point.y;

	  if((createcommand==C_ADD && createitem==C_ELEMENT) ||
		 globalstatus==ADDELEMWITHSAFETY)
	  {
		wparam = MAKEWPARAM((WORD)IDM_POPDECIDEELEM,(WORD)0);
		SendMessage((wdraw.childs+1)->hwnd,WM_COMMAND,wparam,0);
	  }
	  else
	  {
		popupmenudraw(hwnd,point);
	  }
	  break;

	/*Zoom with Wheel*/
	#if 1
	case WM_MOUSEWHEEL:
	//if (GetKeyState(VK_SHIFT) < 0)
	{
		point.x = LOWORD(lParam);
		point.y = HIWORD(lParam);
		x = point.x;
		y = point.y;
		initv = (wdraw.childs + 1)->vparam;
		scrollunit = initv.odv / 7;

		fflag=0;
		if (LOWORD(wParam) == MK_CONTROL) fflag=1;
//        fflag=1;
		rotation = int(HIWORD(wParam));

		if(fflag==1)   //underconstruction
		{
			x0=(wdraw.childs + 1)->vparam.Xo;
			y0=(wdraw.childs + 1)->vparam.Yo;
			len=sqrt(pow((x0-x),2)+pow((y0-y),2));
			sprintf(str,"len:%f",len);
			errormessage(str);
		}

		if (rotation <= 10000)
		{
			(wdraw.childs + 1)->vparam.odv = initv.odv + scrollunit;
			sprintf(str, "%.1f", (wdraw.childs + 1)->vparam.odv);
			SetDlgItemText((wmenu.childs + 2)->hwnd, IDV_L, str);
		}
		else if (rotation > 10000)
		{
			(wdraw.childs + 1)->vparam.odv = initv.odv - scrollunit;
			sprintf(str, "%.1f", (wdraw.childs + 1)->vparam.odv);
			SetDlgItemText((wmenu.childs + 2)->hwnd, IDV_L, str);
		}

		clearwindow(*(wdraw.childs + 1));
		if ((wmenu.childs + 2)->vparam.vflag.mv.ftype == F_ARCLM ||
			(wmenu.childs + 2)->vparam.vflag.mv.ftype == F_FRAME)
		{
			drawarclmframe((wdraw.childs + 1)->hdcC, (wdraw.childs + 1)->vparam,
				arc, code, ONSCREEN);
		}
		if ((wmenu.childs + 2)->vparam.vflag.mv.ftype == F_ORGAN)
		{
			draworganization((wdraw.childs + 1)->hdcC,
				(wdraw.childs + 1)->vparam, (wdraw.childs + 1)->org, ONSCREEN);
		}
		SendMessage((wdraw.childs + 1)->hwnd, WM_PAINT, 0, 0);
	}
    break;
  #endif
	/*Zoom with Wheel*/

//Modified by fukushima 2014/4/30////////////////////
	case WM_MBUTTONDOWN:
	  x = LOWORD(lParam);
	  y = HIWORD(lParam);
	  point.x = x;
	  point.y = y;

	  if((wdraw.childs+1)->lstatus==ROTATE)
	  {
		initx=x;
		inity=y;

		initv=(wdraw.childs+1)->vparam;

		globalstatus=ROTATE;
	  }
	  else if((wdraw.childs+1)->lstatus==MOVE)
	  {
		initx=x;
		inity=y;

		initv=(wdraw.childs+1)->vparam;

		globalstatus=MOVE;
	  }
	  break;
//Modified by fukushima 2014/4/30////////////////////

	case WM_MOUSEMOVE:
	  /*WHILE ADDING ELEMENT.*/
	  if((wmenu.childs+2)->vparam.vflag.mv.ftype==F_ORGAN)
	  {
		if(createcommand==C_ADD ||
		   globalstatus==ADDELEMWITHSAFETY)
		{
		  x = LOWORD(lParam);
		  y = HIWORD(lParam);

		  if(createitem==C_ELEMENT && gelem.nnod>0)
		  {
			if(gelem.nnod<=3)
			{
			  nlast=findlastnode(0,(wdraw.childs+1)->vparam,
								 x,y,
								 NULL,
								 &((wdraw.childs+1)->org));
			}
			else if(gelem.nban>0)
			{
			  bantoplane(*(gelem.bans+0),&pban);

			  nlast=findlastnode(0,(wdraw.childs+1)->vparam,
								 x,y,
								 &pban,
								 &((wdraw.childs+1)->org));
			}

			if(nlast.code==0)
			{
			  nlast.d[0]=0.0;
			  nlast.d[1]=0.0;
			  nlast.d[2]=0.0;
			}

			chaseelement(&gelem,(wdraw.childs+1)->vparam,
						 *(wdraw.childs+1));
		  }
		  break;
		}
		else if((globalstatus==MOVENODEBYMOUSEBEGIN ||
				 globalstatus==PINCHNODEBYMOUSEBEGIN ||
				 globalstatus==PINCHNODEWITHSAFETYBEGIN ||
				 globalstatus==PINCHNODEVENEZIABEGIN ||
				 globalstatus==PINCHNODESTACKBEGIN ||
				 globalstatus==PINCHELEMTWISTBEGIN ||
				 globalstatus==PINCHNODEELASTICBEGIN ||
				 globalstatus==MOVENODEWITHSAFETYBEGIN ||
				 globalstatus==MOVEELEMENTBYMOUSEBEGIN ||
				 globalstatus==COPYELEMENTBYMOUSEBEGIN) &&
				(wParam == MK_LBUTTON))
		{
		  point.x = LOWORD(lParam);
		  point.y = HIWORD(lParam);

		  x=initx-point.x;
		  y=inity-point.y;

		  direction=getdirection((wdraw.childs+1)->vparam,
								 x,y,
								 (wdraw.childs+1)->vparam.focus,
								 &length);
		  if(direction == DIRECTIONX)
		  {
			gincrement.dc[GX]=length;
		  }
		  else if(direction == DIRECTIONY)
		  {
			gincrement.dc[GY]=length;
		  }
		  else if(direction == DIRECTIONZ)
		  {
			gincrement.dc[GZ]=length;
		  }

		  if(globalstatus==MOVENODEBYMOUSEBEGIN ||
			 globalstatus==PINCHNODEBYMOUSEBEGIN ||
			 globalstatus==PINCHNODEWITHSAFETYBEGIN ||
			 globalstatus==PINCHNODEVENEZIABEGIN ||
			 globalstatus==PINCHNODESTACKBEGIN ||
			 globalstatus==PINCHNODEELASTICBEGIN ||
			 globalstatus==MOVENODEWITHSAFETYBEGIN)
		  {
			if(globalstatus==PINCHNODEVENEZIABEGIN)
			{
#if 1
			  gnode->d[GX]=pnode->d[GX]-1.0*gincrement.dc[GX];
			  gnode->d[GY]=pnode->d[GY]-1.0*gincrement.dc[GY];
			  gnode->d[GZ]=pnode->d[GZ]-0.0*gincrement.dc[GZ];
#endif
#if 0
			  gnode->d[GX]=pnode->d[GX]-1.0*gincrement.dc[GX];
			  gnode->d[GY]=pnode->d[GY]-1.0*gincrement.dc[GY];
			  gnode->d[GZ]=pnode->d[GZ]-1.0*gincrement.dc[GZ];
#endif
			}
			else if(globalstatus==PINCHNODESTACKBEGIN)
			{
			  gnode->d[GX]=pnode->d[GX]-1.0*gincrement.dc[GX];
			  gnode->d[GY]=pnode->d[GY]-1.0*gincrement.dc[GY];
			  gnode->d[GZ]=pnode->d[GZ]-0.0*gincrement.dc[GZ];
			}
			else
			{
/*FOR SENJU*/
#if 0
find=0;
for(i=0;i<(wdraw.childs+1)->org.nelem;i++)
{
  pelem=(wdraw.childs+1)->org.elems+i;
  for(j=0;j<(pelem->nnod);j++)
  {
	if((*(pelem->nods+j))->code==gnode->code)
	{
	  find=1;
	  break;
	}
  }
  if(find) break;
}
if(find && pelem->sect->code==201)
{
  gnode->d[GX]=pnode->d[GX]-1.0*gincrement.dc[GX];
  gnode->d[GY]=pnode->d[GY]-0.0*gincrement.dc[GY];
  gnode->d[GZ]=pnode->d[GZ]-0.0*gincrement.dc[GZ];
}
else if(find && pelem->sect->code==202)
{
  gnode->d[GX]=pnode->d[GX]-0.0*gincrement.dc[GX];
  gnode->d[GY]=pnode->d[GY]-1.0*gincrement.dc[GY];
  gnode->d[GZ]=pnode->d[GZ]-0.0*gincrement.dc[GZ];
}
#endif

#if 1
			  gnode->d[GX]=pnode->d[GX]-1.0*gincrement.dc[GX];
			  gnode->d[GY]=pnode->d[GY]-1.0*gincrement.dc[GY];
			  gnode->d[GZ]=pnode->d[GZ]-1.0*gincrement.dc[GZ];
#endif
#if 0
			  gnode->d[GX]=pnode->d[GX]-0.0*gincrement.dc[GX];
			  gnode->d[GY]=pnode->d[GY]-0.0*gincrement.dc[GY];
			  gnode->d[GZ]=pnode->d[GZ]-1.0*gincrement.dc[GZ];
#endif
#if 0
			  gnode->d[GX]=pnode->d[GX]-1.0*gincrement.dc[GX];
			  gnode->d[GY]=pnode->d[GY]-1.0*gincrement.dc[GY];
			  gnode->d[GZ]=pnode->d[GZ]-0.0*gincrement.dc[GZ];
#endif
			}
			/*else
			{
			  gnode->d[GX]=pnode->d[GX]-1.0*gincrement.dc[GX];
			  gnode->d[GY]=pnode->d[GY]-1.0*gincrement.dc[GY];
			  gnode->d[GZ]=pnode->d[GZ]-1.0*gincrement.dc[GZ];
			}*/

			globalflag=1;
			chasenode(gnode,(wdraw.childs+1)->vparam,
					  *(wdraw.childs+1));
			globalflag=0;
		  }
		  else if(globalstatus==MOVEELEMENTBYMOUSEBEGIN ||
				  globalstatus==COPYELEMENTBYMOUSEBEGIN)
		  {
			for(i=0;i<pelem->nnod;i++)
			{
			  (*(pelem->nods+i))->d[GX]=(pnode+i)->d[GX]
										-gincrement.dc[GX];
			  (*(pelem->nods+i))->d[GY]=(pnode+i)->d[GY]
										-gincrement.dc[GY];
			  (*(pelem->nods+i))->d[GZ]=(pnode+i)->d[GZ]
										-gincrement.dc[GZ];
			}
			chaseelement(pelem,(wdraw.childs+1)->vparam,
						 *(wdraw.childs+1));
		  }
		  else if(globalstatus==PINCHELEMTWISTBEGIN)
		  {
			pelem->cangle-=PI/100.0*gincrement.dc[GX];
			chaseelement(pelem,(wdraw.childs+1)->vparam,
						 *(wdraw.childs+1));

			nnode=(wdraw.childs+1)->org.nnode;
			nelem=(wdraw.childs+1)->org.nelem;
			nflag=(long int *)malloc(nnode*sizeof(long int));
			for(i=0;i<nnode;i++) *(nflag+i)=0;

			*(nflag+(*(pelem->nods+0))->loff)=1;
			*(nflag+(*(pelem->nods+1))->loff)=1;

			for(i=1;i<=7;i++)
			{
			  for(j=0;j<nelem;j++)
			  {
				foundelem=(wdraw.childs+1)->org.elems+j;

				if(foundelem->nnod==2)
				{
				  loff=(*(foundelem->nods+0))->loff;
				  moff=(*(foundelem->nods+1))->loff;
				  if(*(nflag+loff)==i && *(nflag+moff)==0)
				  {
					*(nflag+moff)=i+1;
				  }
				  if(*(nflag+moff)==i && *(nflag+loff)==0)
				  {
					*(nflag+loff)=i+1;
				  }
				}
			  }
			}

			for(i=0;i<nelem;i++)
			{
			  find=0;
			  loff=(*(((wdraw.childs+1)->org.elems+i)->nods+0))->loff;
			  moff=(*(((wdraw.childs+1)->org.elems+i)->nods+1))->loff;
			  if((*(nflag+loff)==1 && *(nflag+moff)==2) ||
				 (*(nflag+loff)==2 && *(nflag+moff)==1))
			  {
				((wdraw.childs+1)->org.elems+i)->cangle-=0.833*PI/100.0*gincrement.dc[GX];
				find=1;
			  }
			  if((*(nflag+loff)==2 && *(nflag+moff)==3) ||
				 (*(nflag+loff)==3 && *(nflag+moff)==2))
			  {
				((wdraw.childs+1)->org.elems+i)->cangle-=0.667*PI/100.0*gincrement.dc[GX];
				find=1;
			  }
			  if((*(nflag+loff)==3 && *(nflag+moff)==4) ||
				 (*(nflag+loff)==4 && *(nflag+moff)==3))
			  {
				((wdraw.childs+1)->org.elems+i)->cangle-=0.500*PI/100.0*gincrement.dc[GX];
				find=1;
			  }
			  if((*(nflag+loff)==4 && *(nflag+moff)==5) ||
				 (*(nflag+loff)==5 && *(nflag+moff)==4))
			  {
				((wdraw.childs+1)->org.elems+i)->cangle-=0.333*PI/100.0*gincrement.dc[GX];
				find=1;
			  }
			  if((*(nflag+loff)==5 && *(nflag+moff)==6) ||
				 (*(nflag+loff)==6 && *(nflag+moff)==5))
			  {
				((wdraw.childs+1)->org.elems+i)->cangle-=0.167*PI/100.0*gincrement.dc[GX];
				find=1;
			  }

			  /*if(find)
			  {
				chaseelement((wdraw.childs+1)->org.elems+i,
							 (wdraw.childs+1)->vparam,*(wdraw.childs+1));
			  }*/
			}
		  }
		}
		else if(globalstatus==MOVENODEBYMOUSEBEGIN ||
				globalstatus==PINCHNODEBYMOUSEBEGIN ||
				globalstatus==PINCHNODEWITHSAFETYBEGIN ||
				globalstatus==PINCHNODEVENEZIABEGIN ||
				globalstatus==PINCHNODESTACKBEGIN ||
				globalstatus==PINCHELEMTWISTBEGIN ||
				globalstatus==PINCHNODEELASTICBEGIN ||
				globalstatus==MOVENODEWITHSAFETYBEGIN ||
				globalstatus==MOVEELEMENTBYMOUSEBEGIN ||
				globalstatus==COPYELEMENTBYMOUSEBEGIN) /*END*/
		{
		  if(globalstatus==PINCHNODEBYMOUSEBEGIN ||
			 globalstatus==PINCHNODEWITHSAFETYBEGIN ||
			 globalstatus==PINCHNODEVENEZIABEGIN)
		  {
			gincrement.dc[GX]=pnode->d[GX]-gnode->d[GX];
			gincrement.dc[GY]=pnode->d[GY]-gnode->d[GY];
			gincrement.dc[GZ]=pnode->d[GZ]-gnode->d[GZ];

			nnode=(wdraw.childs+1)->org.nnode;
			nelem=(wdraw.childs+1)->org.nelem;
			nflag=(long int *)malloc(nnode*sizeof(long int));
			for(i=0;i<nnode;i++) *(nflag+i)=0;

			*(nflag+(gnode->loff))=1;

			for(i=1;i<=7;i++)
			{
			  for(j=0;j<nelem;j++)
			  {
				foundelem=(wdraw.childs+1)->org.elems+j;

				if(foundelem->nnod==2)
				{
				  loff=(*(foundelem->nods+0))->loff;
				  moff=(*(foundelem->nods+1))->loff;
				  if(*(nflag+loff)==i && *(nflag+moff)==0)
				  {
					*(nflag+moff)=i+1;
				  }
				  if(*(nflag+moff)==i && *(nflag+loff)==0)
				  {
					*(nflag+loff)=i+1;
				  }
				}
			  }
			}

			radius=0.0;
            getincrement((wmenu.childs+2)->hwnd,&laps,&dsafety);

			for(i=0;i<nnode;i++)
			{
//			  if(*(nflag+i)==4)
			  if(*(nflag+i)==laps+1)
			  {
				dx=pnode->d[GX]-((wdraw.childs+1)->org.nodes+i)->d[GX];
				dy=pnode->d[GY]-((wdraw.childs+1)->org.nodes+i)->d[GY];
				dz=pnode->d[GZ]-((wdraw.childs+1)->org.nodes+i)->d[GZ];
				distance=sqrt(dx*dx+dy*dy+dz*dz);
				if(radius<distance) radius=distance;
			  }
			}

			for(i=0;i<nnode;i++)
			{
#if 0
			  if(*(nflag+i)==2) /*FOR TEPCO*/
			  {
				((wdraw.childs+1)->org.nodes+i)->d[GX]-=1.000*gincrement.dc[GX]; /*1.000*/
				((wdraw.childs+1)->org.nodes+i)->d[GY]-=1.000*gincrement.dc[GY];
				((wdraw.childs+1)->org.nodes+i)->d[GZ]-=1.000*gincrement.dc[GZ];
			  }
#endif
#if 0
			  if(*(nflag+i)==2) /*FOR VENEZIA*/
			  {
				((wdraw.childs+1)->org.nodes+i)->d[GX]-=1.000*gincrement.dc[GX]; /*1.000*/
				((wdraw.childs+1)->org.nodes+i)->d[GY]-=1.000*gincrement.dc[GY];
				((wdraw.childs+1)->org.nodes+i)->d[GZ]-=0.000*gincrement.dc[GZ];
			  }
#endif
#if 0 /*FOR SHIZUOKA*/
			  if(*(nflag+i)==2)
			  {
				((wdraw.childs+1)->org.nodes+i)->d[GX]-=0.000*gincrement.dc[GX]; /*0.950*/
				((wdraw.childs+1)->org.nodes+i)->d[GY]-=0.000*gincrement.dc[GY];
				((wdraw.childs+1)->org.nodes+i)->d[GZ]-=0.400*gincrement.dc[GZ];
			  }
			  if(*(nflag+i)==3)
			  {
				((wdraw.childs+1)->org.nodes+i)->d[GX]-=0.000*gincrement.dc[GX]; /*0.800*/
				((wdraw.childs+1)->org.nodes+i)->d[GY]-=0.000*gincrement.dc[GY];
				((wdraw.childs+1)->org.nodes+i)->d[GZ]-=0.000*gincrement.dc[GZ];
			  }
			  if(*(nflag+i)==4)
			  {
				((wdraw.childs+1)->org.nodes+i)->d[GX]-=0.000*gincrement.dc[GX]; /*0.500*/
				((wdraw.childs+1)->org.nodes+i)->d[GY]-=0.000*gincrement.dc[GY];
				((wdraw.childs+1)->org.nodes+i)->d[GZ]-=0.000*gincrement.dc[GZ];
			  }
#endif
#if 0
			  if(*(nflag+i)==2)
			  {
				((wdraw.childs+1)->org.nodes+i)->d[GX]-=0.833*gincrement.dc[GX]; /*0.950*/
				((wdraw.childs+1)->org.nodes+i)->d[GY]-=0.833*gincrement.dc[GY];
				((wdraw.childs+1)->org.nodes+i)->d[GZ]-=0.833*gincrement.dc[GZ];
			  }
			  if(*(nflag+i)==3)
			  {
				((wdraw.childs+1)->org.nodes+i)->d[GX]-=0.667*gincrement.dc[GX]; /*0.800*/
				((wdraw.childs+1)->org.nodes+i)->d[GY]-=0.667*gincrement.dc[GY];
				((wdraw.childs+1)->org.nodes+i)->d[GZ]-=0.667*gincrement.dc[GZ];
			  }
			  if(*(nflag+i)==4)
			  {
				((wdraw.childs+1)->org.nodes+i)->d[GX]-=0.500*gincrement.dc[GX]; /*0.500*/
				((wdraw.childs+1)->org.nodes+i)->d[GY]-=0.500*gincrement.dc[GY];
				((wdraw.childs+1)->org.nodes+i)->d[GZ]-=0.500*gincrement.dc[GZ];
			  }
			  if(*(nflag+i)==5)
			  {
				((wdraw.childs+1)->org.nodes+i)->d[GX]-=0.333*gincrement.dc[GX]; /*0.200*/
				((wdraw.childs+1)->org.nodes+i)->d[GY]-=0.333*gincrement.dc[GY];
				((wdraw.childs+1)->org.nodes+i)->d[GZ]-=0.333*gincrement.dc[GZ];
			  }
			  if(*(nflag+i)==6)
			  {
				((wdraw.childs+1)->org.nodes+i)->d[GX]-=0.167*gincrement.dc[GX]; /*0.050*/
				((wdraw.childs+1)->org.nodes+i)->d[GY]-=0.167*gincrement.dc[GY];
				((wdraw.childs+1)->org.nodes+i)->d[GZ]-=0.167*gincrement.dc[GZ];
			  }
#endif
#if 1 /*SMOOTH CURVE COSINE*/
			  if(1<*(nflag+i) && *(nflag+i)<=12 &&
				globalstatus==PINCHNODEWITHSAFETYBEGIN)
			  {
				dx=pnode->d[GX]-((wdraw.childs+1)->org.nodes+i)->d[GX];
				dy=pnode->d[GY]-((wdraw.childs+1)->org.nodes+i)->d[GY];
				dz=pnode->d[GZ]-((wdraw.childs+1)->org.nodes+i)->d[GZ];
				distance=sqrt(dx*dx+dy*dy+dz*dz);
				if(distance<radius)
				{
				  rrate=0.5+0.5*cos(PI*distance/radius);

				  ((wdraw.childs+1)->org.nodes+i)->d[GX]-=rrate*gincrement.dc[GX];
				  ((wdraw.childs+1)->org.nodes+i)->d[GY]-=rrate*gincrement.dc[GY];
				  ((wdraw.childs+1)->org.nodes+i)->d[GZ]-=rrate*gincrement.dc[GZ];
				}
			  }
#endif
#if 0 /*POINTED CURVE PARABOLA*/
			  if(1<*(nflag+i) && *(nflag+i)<=12 &&
				globalstatus==PINCHNODEWITHSAFETYBEGIN)
			  {
				dx=pnode->d[GX]-((wdraw.childs+1)->org.nodes+i)->d[GX];
				dy=pnode->d[GY]-((wdraw.childs+1)->org.nodes+i)->d[GY];
				dz=pnode->d[GZ]-((wdraw.childs+1)->org.nodes+i)->d[GZ];
				distance=sqrt(dx*dx+dy*dy+dz*dz);
				if(distance<radius)
				{
				  rrate=(1-distance/radius)*(1-distance/radius);

				  ((wdraw.childs+1)->org.nodes+i)->d[GX]-=rrate*gincrement.dc[GX];
				  ((wdraw.childs+1)->org.nodes+i)->d[GY]-=rrate*gincrement.dc[GY];
				  ((wdraw.childs+1)->org.nodes+i)->d[GZ]-=rrate*gincrement.dc[GZ];
				}
			  }
#endif
            }
			free(nflag);
		  }

		  if(globalstatus==PINCHNODESTACKBEGIN &&
			 globalmessageflag==1)
          {
			globalmessageflag=0;

            gincrement.dc[GX]=pnode->d[GX]-gnode->d[GX];
            gincrement.dc[GY]=pnode->d[GY]-gnode->d[GY];
			gincrement.dc[GZ]=pnode->d[GZ]-gnode->d[GZ];

            nnode=(wdraw.childs+1)->org.nnode;
			nelem=(wdraw.childs+1)->org.nelem;
            nflag=(long int *)malloc(nnode*sizeof(long int));
            for(i=0;i<nnode;i++) *(nflag+i)=0;

            soff=0;

            for(j=0;j<nelem;j++) /*SEARCH PINCHED SECT*/
            {
              foundelem=(wdraw.childs+1)->org.elems+j;
              if(foundelem->nnod==2)
              {
                loff=(*(foundelem->nods+0))->loff;
                moff=(*(foundelem->nods+1))->loff;
                if(loff==gnode->loff || moff==gnode->loff)
                {
                  if(soff==0){soff=foundelem->sect->loff;}
                  if(soff!=0 && soff!=foundelem->sect->loff)
                  {
					soff=0;
                    break;
                  }
                }
              }
            }

            for(i=1;i<=2;i++)
            {
              for(j=0;j<nelem;j++) /*SEARCH SAME SECT*/
              {
                foundelem=(wdraw.childs+1)->org.elems+j;
                if(foundelem->nnod==2)
                {
                  loff=(*(foundelem->nods+0))->loff;
                  moff=(*(foundelem->nods+1))->loff;
/*
if(loff>=nnode || moff>=nnode)
{
  sprintf(str,"Nodes=%d Elems=%d ELEM=%d NODE=%d,%d loff=%ld moff=%ld",
          nnode,nelem,
		  foundelem->code,
          (*(foundelem->nods+0))->code,
          (*(foundelem->nods+1))->code,
          loff,moff);
  MessageBox(NULL,str,"Stack",MB_OK);
}
*/
				  if(*(nflag+loff)<=1 &&
                     soff==foundelem->sect->loff)
				  {
                    *(nflag+loff)=1;
                  }
                  if(*(nflag+moff)<=1 &&
                     soff==foundelem->sect->loff)
                  {
                    *(nflag+moff)=1;
				  }
                  if(i==2 &&
                     *(nflag+loff)>=1 &&
                     soff!=foundelem->sect->loff)
                  {
                    (*(nflag+loff))++;
                  }
                  if(i==2 &&
                     *(nflag+moff)>=1 &&
                     soff!=foundelem->sect->loff)
                  {
                    (*(nflag+moff))++;
                  }
                }
              }
            }

            *(nflag+(gnode->loff))=0; /*PINCHED NODE*/

            for(i=0;i<nnode;i++)
            {
              if(1<=*(nflag+i) && *(nflag+i)<=2)
              {
                ((wdraw.childs+1)->org.nodes+i)->d[GX]-=1.000*gincrement.dc[GX]; /*1.000*/
                ((wdraw.childs+1)->org.nodes+i)->d[GY]-=1.000*gincrement.dc[GY];
                ((wdraw.childs+1)->org.nodes+i)->d[GZ]-=0.000*gincrement.dc[GZ];
			  }
              if(*(nflag+i)>=3) /*UNDER CONSTRUCTION*/
              {
                /*METHOD 1*/
                /*CHOOSE MULTI ELEM*/
                /*COPY MULTI ELEM*/
                /*DELETE PREVIOUS MULTI ELEM*/
                /*CONNECT NODES*/
                /*CREATE NODE ON CROSSING POINT*/

                /*METHOD 2*/
                /*CREATE ORGAN WITHOUT NODE ON CRASSING POINT*/
                /*DEVIDE ELEMENT BY CRASSING POINT IN EXTRACTING ARCLM*/
              }
            }

            free(nflag);
          }

		  if(globalstatus==MOVEELEMENTBYMOUSEBEGIN ||
			 globalstatus==COPYELEMENTBYMOUSEBEGIN ||
			 globalstatus==PINCHELEMTWISTBEGIN)
		  {
			pelem->color.line.r=pcolor.r;
			pelem->color.line.g=pcolor.g;
			pelem->color.line.b=pcolor.b;
		  }
		  free(pnode);

		  if(globalstatus==COPYELEMENTBYMOUSEBEGIN)
		  {
			addelementwithnode(pelem,
							   &((wdraw.childs+1)->org),
							   0.0,0.0,0.0);
							   /*UNDER CONSTRUCTION*/

			/*FREE*/
			for(i=0;i<melem.nban;i++)
			{
			  free((melem.bans+i)->nods);
			}
			free(melem.bans);
			free(melem.bonds);
			free(melem.nods);
		  }

/************/
		  if(globalstatus==PINCHNODEWITHSAFETYBEGIN ||
			 globalstatus==PINCHNODESTACKBEGIN ||
			 globalstatus==PINCHELEMTWISTBEGIN ||
			 globalstatus==MOVENODEWITHSAFETYBEGIN)
		  {
			globalstatus=prestatus;
			globalmessageflag=0;

			/*INITIAL*/
			arc =arci;
			arcx=arci;
			arcy=arci;
			free((wdraw.childs+1)->org.loads);
			for(i=0;i<(wdraw.childs+1)->org.nelem;i++)
			{
			  for(j=0;j<2;j++)
			  {
				for(k=0;k<6;k++)
				{
				  ((wdraw.childs+1)->org.elems+i)->initial[j][k]=0.0;
				}
			  }
			} /*UNAVAILABLE INITIAL CMQ.*/

			if(globalstatus!=PINCHNODESTACK)
			{
			  pnelem=(wdraw.childs+1)->org.nelem;
			  pelems=(wdraw.childs+1)->org.elems;
			}

			/*EXTRACT ARCLM*/
			extractarclmfromorgan(&((wdraw.childs+1)->org),
								  &arc,&arcx,&arcy);

			if(globalstatus!=PINCHNODESTACK)
			{
			  (wdraw.childs+1)->org.nelem=pnelem;
			  (wdraw.childs+1)->org.elems=pelems;
			}

			(wmenu.childs+2)->vparam.vflag.mv.ftype=F_ARCLM;

			/*SAVE AS ARCLM*/
			saveasarclm((wdraw.childs+1)->inpfilez,&arc);
			saveasarclm((wdraw.childs+1)->inpfilex,&arcx);
			saveasarclm((wdraw.childs+1)->inpfiley,&arcy);

			/*ANALYSIS*/
			getviewparam((wmenu.childs+2)->hwnd,
						 &((wdraw.childs+1)->vparam));
			(wdraw.childs+1)->vparam.vflag.ev.deformation=0;

			arclm001(&arc ,ID_INPUTFILEZ,ID_OUTPUTFILEZ);
			arclm001(&arcx,ID_INPUTFILEX,ID_OUTPUTFILEX);
			arclm001(&arcy,ID_INPUTFILEY,ID_OUTPUTFILEY);

			/*SRCAN*/
			n=strcspn((wdraw.childs+1)->sctfile,".");
			strncpy(str,(wdraw.childs+1)->sctfile,n);
			str[n]='\0';

			i=srcan001(str);

			if(i==0) MessageBox(NULL,"Failed.","SRCAN001",MB_OK);
/*
			strcpy((wdraw.childs+1)->otpfile,str);
			strcat((wdraw.childs+1)->otpfile,".rat");
			SetDlgItemText((wmenu.childs+2)->hwnd,ID_OUTPUTFILE,
						   (wdraw.childs+1)->otpfile);
*/
			fout=fopen((wdraw.childs+1)->otpfile,"r");
			if(fout==NULL) break;
			readsrcanrate(fout,&arc);

			/*SAFETY FLAG*/
			(wdraw.childs+1)->vparam.vflag.ev.srcancolor=1;

			clearwindow(*(wdraw.childs+1));
			drawarclmframe((wdraw.childs+1)->hdcC,
						   (wdraw.childs+1)->vparam,arc,0,ONSCREEN);


            fout=fopen((wdraw.childs+1)->otpfile,"r");
            if(fout==NULL) break;
            readsrcanrate(fout,&arc);
            fclose(fout);
            int nred=0;
            double maxsafety=0.0;
            int ii,jj;
            for(ii=0;ii<arc.nelem;ii++)
            {

              //if((arc.elems+ii)->srate[2]>1) nred=nred+1;
              if( (arc.elems+ii)->srate[0]>1 ||
                  (arc.elems+ii)->srate[1]>1 ||
                  (arc.elems+ii)->srate[2]>1 ||
                  (arc.elems+ii)->srate[3]>1 )
                {
                  nred=nred+1;
                }

				for(jj=0;jj<4;jj++)
				{
				  if(maxsafety<(arc.elems+ii)->srate[jj]) maxsafety=(arc.elems+ii)->srate[jj];
				}

            }


			/*draworganization((wdraw.childs+1)->hdcC,
							 (wdraw.childs+1)->vparam,
							 (wdraw.childs+1)->org,ONSCREEN);*/
			SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);
			/*MessageBox(NULL,"Completed.","Move Node with Safety",MB_OK);*/

			(wmenu.childs+2)->vparam.vflag.mv.ftype=F_ORGAN;
			SendMessage((wmenu.childs+2)->hwnd,WM_INITDIALOG,0,0);

            sprintf(str,"%d",laps);
			SetDlgItemText((wmenu.childs+2)->hwnd,ID_LAPS,str);
			SendDlgItemMessage((wmenu.childs+2)->hwnd,ID_LAPS,WM_PAINT,laps,0);

			sprintf(str,"%d",nred);
			SetDlgItemText((wmenu.childs+2)->hwnd,ID_SAFETY,str);
			SendDlgItemMessage((wmenu.childs+2)->hwnd,ID_SAFETY,WM_PAINT,0,0);

			globalmessageflag=1;
		  }

		  else if(globalstatus==PINCHNODEELASTICBEGIN)
		  {
			globalstatus=prestatus;
			globalmessageflag=0;

			gincrement.dc[GX]=pnode->d[GX]-gnode->d[GX];
			gincrement.dc[GY]=pnode->d[GY]-gnode->d[GY];
			gincrement.dc[GZ]=pnode->d[GZ]-gnode->d[GZ];

			/*CASE ONLY SHAPE WITHOUT SUM OF STRESS*/
			/*INITIAL*/
			arc =arci;
			arcx=arci;
			arcy=arci;
			free((wdraw.childs+1)->org.loads);
			for(i=0;i<(wdraw.childs+1)->org.nelem;i++)
			{
			  for(j=0;j<2;j++)
			  {
				for(k=0;k<6;k++)
				{
				  ((wdraw.childs+1)->org.elems+i)->initial[j][k]=0.0;
				}
			  }
			} /*UNAVAILABLE INITIAL CMQ*/

			pnelem=(wdraw.childs+1)->org.nelem;
			pelems=(wdraw.childs+1)->org.elems;

			/*INITIALIZE NODE COORDINATION*/
			gnode->d[GX]+=1.0*gincrement.dc[GX];
			gnode->d[GY]+=1.0*gincrement.dc[GY];
			gnode->d[GZ]+=1.0*gincrement.dc[GZ];

			/*EXTRACT ARCLM*/
			extractarclmfromorgan(&((wdraw.childs+1)->org),
								  &arc,&arcx,&arcy);

			(wdraw.childs+1)->org.nelem=pnelem;
			(wdraw.childs+1)->org.elems=pelems;

			/*INPUT NODE TRANSFORMATION*/
			for(i=0;i<(arc.nnode);i++)
			{
			  if((arc.nodes+i)->code==gnode->code)
			  {
				(arc.confs+6*i+0)->iconf=1;
				(arc.confs+6*i+1)->iconf=1;
				(arc.confs+6*i+2)->iconf=1;
				(arc.confs+6*i+0)->value=-1.0*gincrement.dc[GX];
				(arc.confs+6*i+1)->value=-1.0*gincrement.dc[GY];
				(arc.confs+6*i+2)->value=-1.0*gincrement.dc[GZ];
			  }
			  else
			  {
				(arc.confs+6*i+0)->value=0.0;
				(arc.confs+6*i+1)->value=0.0;
				(arc.confs+6*i+2)->value=0.0;
			  }
			}
			/*UNAVAILABLE CMQ*/
			for(i=0;i<(arc.nelem);i++)
			{
			  for(j=0;j<2;j++)
			  {
				for(k=0;k<6;k++) (arc.elems+i)->stress[j][k]=0.0;
			  }
			}

			(wmenu.childs+2)->vparam.vflag.mv.ftype=F_ARCLM;

			/*SAVE AS ARCLM*/
			saveasarclm((wdraw.childs+1)->inpfilez,&arc);

			/*ANALYSIS*/
			getviewparam((wmenu.childs+2)->hwnd,
						 &((wdraw.childs+1)->vparam));
			(wdraw.childs+1)->vparam.vflag.ev.deformation=1;

			arclm001(&arc,ID_INPUTFILEZ,ID_OUTPUTFILEZ);

			/*FORMATION UPDATE*/
			j=0;
			for(i=0;i<arc.nnode;i++)
			{
			  while((arc.nodes+i)->code
					!=((wdraw.childs+1)->org.nodes+j)->code) j++;
			  if(j>=((wdraw.childs+1)->org.nnode)) break;
			  ((wdraw.childs+1)->org.nodes+j)->d[GX]=*(arc.ddisp+6*i+0);
			  ((wdraw.childs+1)->org.nodes+j)->d[GY]=*(arc.ddisp+6*i+1);
			  ((wdraw.childs+1)->org.nodes+j)->d[GZ]=*(arc.ddisp+6*i+2);
			}

			/*LOAD UPDATE*/






#if 1
			/*INITIAL*/
			arc =arci;
			arcx=arci;
			arcy=arci;
			free((wdraw.childs+1)->org.loads);
			for(i=0;i<(wdraw.childs+1)->org.nelem;i++)
			{
			  for(j=0;j<2;j++)
			  {
				for(k=0;k<6;k++)
				{
				  ((wdraw.childs+1)->org.elems+i)->initial[j][k]=0.0;
				}
			  }
			} /*UNAVAILABLE INITIAL CMQ.*/

			pnelem=(wdraw.childs+1)->org.nelem;
			pelems=(wdraw.childs+1)->org.elems;

			/*EXTRACT ARCLM*/
			extractarclmfromorgan(&((wdraw.childs+1)->org),
								  &arc,&arcx,&arcy);

			(wdraw.childs+1)->org.nelem=pnelem;
			(wdraw.childs+1)->org.elems=pelems;

			(wmenu.childs+2)->vparam.vflag.mv.ftype=F_ARCLM;

			/*SAVE AS ARCLM*/
			saveasarclm((wdraw.childs+1)->inpfilez,&arc);
			saveasarclm((wdraw.childs+1)->inpfilex,&arcx);
			saveasarclm((wdraw.childs+1)->inpfiley,&arcy);

			/*ANALYSIS*/
			getviewparam((wmenu.childs+2)->hwnd,
						 &((wdraw.childs+1)->vparam));
			(wdraw.childs+1)->vparam.vflag.ev.deformation=0;

			arclm001(&arc ,ID_INPUTFILEZ,ID_OUTPUTFILEZ);
			arclm001(&arcx,ID_INPUTFILEX,ID_OUTPUTFILEX);
			arclm001(&arcy,ID_INPUTFILEY,ID_OUTPUTFILEY);

			/*SRCAN*/
			n=strcspn((wdraw.childs+1)->sctfile,".");
			strncpy(str,(wdraw.childs+1)->sctfile,n);
			str[n]='\0';

			i=srcan001(str);

			if(i==0) MessageBox(NULL,"Failed.","SRCAN001",MB_OK);

			strcpy((wdraw.childs+1)->otpfile,str);
			strcat((wdraw.childs+1)->otpfile,".rat");
			SetDlgItemText((wmenu.childs+2)->hwnd,ID_OUTPUTFILE,
						   (wdraw.childs+1)->otpfile);

			fout=fopen((wdraw.childs+1)->otpfile,"r");
			if(fout==NULL) break;
			readsrcanrate(fout,&arc);

			/*SAFETY FLAG*/
			(wdraw.childs+1)->vparam.vflag.ev.srcancolor=1;

			clearwindow(*(wdraw.childs+1));
			drawarclmframe((wdraw.childs+1)->hdcC,
						   (wdraw.childs+1)->vparam,arc,0,ONSCREEN);
#endif

			/*draworganization((wdraw.childs+1)->hdcC,
							 (wdraw.childs+1)->vparam,
							 (wdraw.childs+1)->org,ONSCREEN);*/
			SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);
			/*MessageBox(NULL,"Completed.","Move Node with Safety",MB_OK);*/

			(wmenu.childs+2)->vparam.vflag.mv.ftype=F_ORGAN;
			SendMessage((wmenu.childs+2)->hwnd,WM_INITDIALOG,0,0);

			globalmessageflag=1;
		  }

#if 1
		  else if(globalstatus==PINCHNODEVENEZIABEGIN)
		  {
			globalstatus=prestatus;
			globalmessageflag=0;
			(wdraw.childs+1)->vparam.vflag.ev.srcanrate =0;

			n=strcspn((wdraw.childs+1)->inpfile,".");
			strncpy(str,(wdraw.childs+1)->inpfile,n);
			str[n]='\0';
			strcpy(fname,str);

		  /*LONG*/
			globalstatus=PINCHNODEVENEZIALONG;
			/*INITIAL*/
			arc =arci;
			arcx=arci;
			arcy=arci;
			free((wdraw.childs+1)->org.loads);
			for(i=0;i<(wdraw.childs+1)->org.nelem;i++)
			{
			  for(j=0;j<2;j++)
			  {
				for(k=0;k<6;k++)
				{
				  ((wdraw.childs+1)->org.elems+i)->initial[j][k]=0.0;
				}
			  }
			} /*UNAVAILABLE INITIAL CMQ.*/

			pnelem=(wdraw.childs+1)->org.nelem;
			pelems=(wdraw.childs+1)->org.elems;

			/*EXTRACT ARCLM*/
			extractarclmfromorgan(&((wdraw.childs+1)->org),
								  &arc,&arcx,&arcy);

			(wdraw.childs+1)->org.nelem=pnelem;
			(wdraw.childs+1)->org.elems=pelems;

			(wmenu.childs+2)->vparam.vflag.mv.ftype=F_ARCLM;

			/*INITIAL Z*/
			for(i=0;i<arc.nnode;i++)
			{
			  if((arc.nodes+i)->code==111) initz[ 0]=(arc.nodes+i)->d[GZ];
			  if((arc.nodes+i)->code==255) initz[ 1]=(arc.nodes+i)->d[GZ];
			  if((arc.nodes+i)->code==267) initz[ 2]=(arc.nodes+i)->d[GZ];
			  if((arc.nodes+i)->code==412) initz[ 3]=(arc.nodes+i)->d[GZ];
			  if((arc.nodes+i)->code==346) initz[ 4]=(arc.nodes+i)->d[GZ];
			  if((arc.nodes+i)->code==348) initz[ 5]=(arc.nodes+i)->d[GZ];
			  if((arc.nodes+i)->code==413) initz[ 6]=(arc.nodes+i)->d[GZ];
			  if((arc.nodes+i)->code==268) initz[ 7]=(arc.nodes+i)->d[GZ];
			  if((arc.nodes+i)->code==256) initz[ 8]=(arc.nodes+i)->d[GZ];
			  if((arc.nodes+i)->code==112) initz[ 9]=(arc.nodes+i)->d[GZ];
			  if((arc.nodes+i)->code==296) initz[10]=(arc.nodes+i)->d[GZ];
			  if((arc.nodes+i)->code==295) initz[11]=(arc.nodes+i)->d[GZ];
			}

			/*SAVE AS ARCLM*/
			strcpy(str,fname);
			strcat(str,"-il.inl");
			strcpy((wdraw.childs+1)->inpfilez,str);
			SetDlgItemText((wmenu.childs+2)->hwnd,ID_INPUTFILEZ,str);
			strcpy(str,fname);
			strcat(str,"-il.ihx");
			strcpy((wdraw.childs+1)->inpfilex,str);
			SetDlgItemText((wmenu.childs+2)->hwnd,ID_INPUTFILEX,str);
			strcpy(str,fname);
			strcat(str,"-il.ihy");
			strcpy((wdraw.childs+1)->inpfiley,str);
			SetDlgItemText((wmenu.childs+2)->hwnd,ID_INPUTFILEY,str);

			saveasarclm((wdraw.childs+1)->inpfilez,&arc);
			saveasarclm((wdraw.childs+1)->inpfilex,&arcx);
			saveasarclm((wdraw.childs+1)->inpfiley,&arcy);

			/*ANALYSIS*/
			getviewparam((wmenu.childs+2)->hwnd,
						 &((wdraw.childs+1)->vparam));
			(wdraw.childs+1)->vparam.vflag.ev.deformation=0;

			strcpy(str,fname);
			strcat(str,".olz");
			SetDlgItemText((wmenu.childs+2)->hwnd,ID_OUTPUTFILEZ,str);
			strcpy(str,fname);
			strcat(str,".olx");
			SetDlgItemText((wmenu.childs+2)->hwnd,ID_OUTPUTFILEX,str);
			strcpy(str,fname);
			strcat(str,".oly");
			SetDlgItemText((wmenu.childs+2)->hwnd,ID_OUTPUTFILEY,str);

			arclm001(&arc ,ID_INPUTFILEZ,ID_OUTPUTFILEZ);
			arclm001(&arcx,ID_INPUTFILEX,ID_OUTPUTFILEX);
			arclm001(&arcy,ID_INPUTFILEY,ID_OUTPUTFILEY);

			/*DISPLACEMENT*/
			for(i=0;i<arc.nnode;i++)
			{
			  if((arc.nodes+i)->code==111) venedisp[0][ 0]=*(arc.ddisp+6*i+2)-initz[ 0];
			  if((arc.nodes+i)->code==255) venedisp[0][ 1]=*(arc.ddisp+6*i+2)-initz[ 1];
			  if((arc.nodes+i)->code==267) venedisp[0][ 2]=*(arc.ddisp+6*i+2)-initz[ 2];
			  if((arc.nodes+i)->code==412) venedisp[0][ 3]=*(arc.ddisp+6*i+2)-initz[ 3];
			  if((arc.nodes+i)->code==346) venedisp[0][ 4]=*(arc.ddisp+6*i+2)-initz[ 4];
			  if((arc.nodes+i)->code==348) venedisp[0][ 5]=*(arc.ddisp+6*i+2)-initz[ 5];
			  if((arc.nodes+i)->code==413) venedisp[0][ 6]=*(arc.ddisp+6*i+2)-initz[ 6];
			  if((arc.nodes+i)->code==268) venedisp[0][ 7]=*(arc.ddisp+6*i+2)-initz[ 7];
			  if((arc.nodes+i)->code==256) venedisp[0][ 8]=*(arc.ddisp+6*i+2)-initz[ 8];
			  if((arc.nodes+i)->code==112) venedisp[0][ 9]=*(arc.ddisp+6*i+2)-initz[ 9];
			  if((arc.nodes+i)->code==296) venedisp[0][10]=*(arc.ddisp+6*i+2)-initz[10];
			  if((arc.nodes+i)->code==295) venedisp[0][11]=*(arc.ddisp+6*i+2)-initz[11];
			}

			/*SRCAN*/
			n=strcspn((wdraw.childs+1)->sctfile,".");
			strncpy(str,(wdraw.childs+1)->sctfile,n);
			str[n]='\0';
			strcat(str,"-ol");

			i=srcan001(str);

			if(i==0) MessageBox(NULL,"Failed.","SRCAN001",MB_OK);

			strcpy((wdraw.childs+1)->otpfile,str);
			strcat((wdraw.childs+1)->otpfile,".rat");
			SetDlgItemText((wmenu.childs+2)->hwnd,ID_OUTPUTFILE,
						   (wdraw.childs+1)->otpfile);

			fout=fopen((wdraw.childs+1)->otpfile,"r");
			if(fout==NULL) break;
			readsrcanrate(fout,&arc);

			srates=(double *)malloc(4*(arc.nelem)*sizeof(double));
			for(i=0;i<arc.nelem;i++)
			{
			  for(j=0;j<4;j++) *(srates+4*i+j)=(arc.elems+i)->srate[j];
			}

		  /*SHORT X POSI*/
			globalstatus=PINCHNODEVENEZIAXPOSI;
			strcpy(str,fname);
			strcat(str,".ixp");
			fin=fopen(str,"r");
			if(fin==NULL) break;
			initializeorganization(&globalorg);
			inputorganization(fin,&globalorg,&((wdraw.childs+1)->vparam));
			fclose(fin);
			str[n]='\0';
			for(i=0;i<(wdraw.childs+1)->org.nnode;i++)
			{
			  (globalorg.nodes+i)->d[GX]=((wdraw.childs+1)->org.nodes+i)->d[GX];
			  (globalorg.nodes+i)->d[GY]=((wdraw.childs+1)->org.nodes+i)->d[GY];
			  (globalorg.nodes+i)->d[GZ]=((wdraw.childs+1)->org.nodes+i)->d[GZ];
			}
			/*INITIAL*/
			arc =arci;
			arcx=arci;
			arcy=arci;
			free(globalorg.loads);
			for(i=0;i<globalorg.nelem;i++)
			{
			  for(j=0;j<2;j++)
			  {
				for(k=0;k<6;k++)
				{
				  (globalorg.elems+i)->initial[j][k]=0.0;
				}
			  }
			} /*UNAVAILABLE INITIAL CMQ.*/
			pnelem=globalorg.nelem;
			pelems=globalorg.elems;
			/*EXTRACT ARCLM*/
			extractarclmfromorgan(&globalorg,
								  &arc,&arcx,&arcy);
			globalorg.nelem=pnelem;
			globalorg.elems=pelems;
			/*(wmenu.childs+2)->vparam.vflag.mv.ftype=F_ARCLM;*/
			/*SAVE AS ARCLM*/
			strcpy(str,fname);
			strcat(str,"-xp.inl");
			strcpy((wdraw.childs+1)->inpfilez,str);
			SetDlgItemText((wmenu.childs+2)->hwnd,ID_INPUTFILEZ,str);
			strcpy(str,fname);
			strcat(str,"-xp.ihx");
			strcpy((wdraw.childs+1)->inpfilex,str);
			SetDlgItemText((wmenu.childs+2)->hwnd,ID_INPUTFILEX,str);
			strcpy(str,fname);
			strcat(str,"-xp.ihy");
			strcpy((wdraw.childs+1)->inpfiley,str);
			SetDlgItemText((wmenu.childs+2)->hwnd,ID_INPUTFILEY,str);
			saveasarclm((wdraw.childs+1)->inpfilez,&arc);
			saveasarclm((wdraw.childs+1)->inpfilex,&arcx);
			saveasarclm((wdraw.childs+1)->inpfiley,&arcy);
			/*ANALYSIS*/
			getviewparam((wmenu.childs+2)->hwnd,
						 &((wdraw.childs+1)->vparam));
			(wdraw.childs+1)->vparam.vflag.ev.deformation=0;
			strcpy(str,fname);
			strcat(str,".xpz");
			SetDlgItemText((wmenu.childs+2)->hwnd,ID_OUTPUTFILEZ,str);
			strcpy(str,fname);
			strcat(str,".xpx");
			SetDlgItemText((wmenu.childs+2)->hwnd,ID_OUTPUTFILEX,str);
			strcpy(str,fname);
			strcat(str,".xpy");
			SetDlgItemText((wmenu.childs+2)->hwnd,ID_OUTPUTFILEY,str);
			arclm001(&arc ,ID_INPUTFILEZ,ID_OUTPUTFILEZ);
			arclm001(&arcx,ID_INPUTFILEX,ID_OUTPUTFILEX);
			arclm001(&arcy,ID_INPUTFILEY,ID_OUTPUTFILEY);
			/*DISPLACEMENT*/
			for(i=0;i<arc.nnode;i++)
			{
			  if((arc.nodes+i)->code==111) venedisp[1][ 0]=2*(*(arcx.ddisp+6*i+2))-2*initz[ 0];
			  if((arc.nodes+i)->code==255) venedisp[1][ 1]=2*(*(arcx.ddisp+6*i+2))-2*initz[ 1];
			  if((arc.nodes+i)->code==267) venedisp[1][ 2]=2*(*(arcx.ddisp+6*i+2))-2*initz[ 2];
			  if((arc.nodes+i)->code==412) venedisp[1][ 3]=2*(*(arcx.ddisp+6*i+2))-2*initz[ 3];
			  if((arc.nodes+i)->code==346) venedisp[1][ 4]=2*(*(arcx.ddisp+6*i+2))-2*initz[ 4];
			  if((arc.nodes+i)->code==348) venedisp[1][ 5]=2*(*(arcx.ddisp+6*i+2))-2*initz[ 5];
			  if((arc.nodes+i)->code==413) venedisp[1][ 6]=2*(*(arcx.ddisp+6*i+2))-2*initz[ 6];
			  if((arc.nodes+i)->code==268) venedisp[1][ 7]=2*(*(arcx.ddisp+6*i+2))-2*initz[ 7];
			  if((arc.nodes+i)->code==256) venedisp[1][ 8]=2*(*(arcx.ddisp+6*i+2))-2*initz[ 8];
			  if((arc.nodes+i)->code==112) venedisp[1][ 9]=2*(*(arcx.ddisp+6*i+2))-2*initz[ 9];
			  if((arc.nodes+i)->code==296) venedisp[1][10]=2*(*(arcx.ddisp+6*i+2))-2*initz[10];
			  if((arc.nodes+i)->code==295) venedisp[1][11]=2*(*(arcx.ddisp+6*i+2))-2*initz[11];
			}
			/*SRCAN*/
			n=strcspn((wdraw.childs+1)->sctfile,".");
			strncpy(str,(wdraw.childs+1)->sctfile,n);
			str[n]='\0';
			strcat(str,"-xp");
			i=srcan001(str);
			if(i==0) MessageBox(NULL,"Failed.","SRCAN001",MB_OK);
			strcpy((wdraw.childs+1)->otpfile,str);
			strcat((wdraw.childs+1)->otpfile,".rat");
			SetDlgItemText((wmenu.childs+2)->hwnd,ID_OUTPUTFILE,
						   (wdraw.childs+1)->otpfile);

			fout=fopen((wdraw.childs+1)->otpfile,"r");
			if(fout==NULL) break;
			readsrcanrate(fout,&arcx);
			for(i=0;i<arcx.nelem;i++)
			{
			  for(j=0;j<4;j++)
			  {
				if(*(srates+4*i+j)<(arcx.elems+i)->srate[j])
				{*(srates+4*i+j)=(arcx.elems+i)->srate[j];}
			  }
			}

		  /*SHORT X NEGA*/
			globalstatus=PINCHNODEVENEZIAXNEGA;
			strcpy(str,fname);
			strcat(str,".ixn");
			fin=fopen(str,"r");
			if(fin==NULL) break;
			initializeorganization(&globalorg);
			inputorganization(fin,&globalorg,&((wdraw.childs+1)->vparam));
			fclose(fin);
			str[n]='\0';
			for(i=0;i<(wdraw.childs+1)->org.nnode;i++)
			{
			  (globalorg.nodes+i)->d[GX]=((wdraw.childs+1)->org.nodes+i)->d[GX];
			  (globalorg.nodes+i)->d[GY]=((wdraw.childs+1)->org.nodes+i)->d[GY];
			  (globalorg.nodes+i)->d[GZ]=((wdraw.childs+1)->org.nodes+i)->d[GZ];
			}
			/*INITIAL*/
			arc =arci;
			arcx=arci;
			arcy=arci;
			free(globalorg.loads);
			for(i=0;i<globalorg.nelem;i++)
			{
			  for(j=0;j<2;j++)
			  {
				for(k=0;k<6;k++)
				{
				  (globalorg.elems+i)->initial[j][k]=0.0;
				}
			  }
			} /*UNAVAILABLE INITIAL CMQ.*/
			pnelem=globalorg.nelem;
			pelems=globalorg.elems;
			/*EXTRACT ARCLM*/
			extractarclmfromorgan(&globalorg,
								  &arc,&arcx,&arcy);
			globalorg.nelem=pnelem;
			globalorg.elems=pelems;
			/*(wmenu.childs+2)->vparam.vflag.mv.ftype=F_ARCLM;*/
			/*SAVE AS ARCLM*/
			strcpy(str,fname);
			strcat(str,"-xn.inl");
			strcpy((wdraw.childs+1)->inpfilez,str);
			SetDlgItemText((wmenu.childs+2)->hwnd,ID_INPUTFILEZ,str);
			strcpy(str,fname);
			strcat(str,"-xn.ihx");
			strcpy((wdraw.childs+1)->inpfilex,str);
			SetDlgItemText((wmenu.childs+2)->hwnd,ID_INPUTFILEX,str);
			strcpy(str,fname);
			strcat(str,"-xn.ihy");
			strcpy((wdraw.childs+1)->inpfiley,str);
			SetDlgItemText((wmenu.childs+2)->hwnd,ID_INPUTFILEY,str);
			saveasarclm((wdraw.childs+1)->inpfilez,&arc);
			saveasarclm((wdraw.childs+1)->inpfilex,&arcx);
			saveasarclm((wdraw.childs+1)->inpfiley,&arcy);
			/*ANALYSIS*/
			getviewparam((wmenu.childs+2)->hwnd,
						 &((wdraw.childs+1)->vparam));
			(wdraw.childs+1)->vparam.vflag.ev.deformation=0;
			strcpy(str,fname);
			strcat(str,".xnz");
			SetDlgItemText((wmenu.childs+2)->hwnd,ID_OUTPUTFILEZ,str);
			strcpy(str,fname);
			strcat(str,".xnx");
			SetDlgItemText((wmenu.childs+2)->hwnd,ID_OUTPUTFILEX,str);
			strcpy(str,fname);
			strcat(str,".xny");
			SetDlgItemText((wmenu.childs+2)->hwnd,ID_OUTPUTFILEY,str);
			arclm001(&arc ,ID_INPUTFILEZ,ID_OUTPUTFILEZ);
			arclm001(&arcx,ID_INPUTFILEX,ID_OUTPUTFILEX);
			arclm001(&arcy,ID_INPUTFILEY,ID_OUTPUTFILEY);
			/*DISPLACEMENT*/
			for(i=0;i<arc.nnode;i++)
			{
			  if((arc.nodes+i)->code==111) venedisp[2][ 0]=2*(*(arcx.ddisp+6*i+2))-2*initz[ 0];
			  if((arc.nodes+i)->code==255) venedisp[2][ 1]=2*(*(arcx.ddisp+6*i+2))-2*initz[ 1];
			  if((arc.nodes+i)->code==267) venedisp[2][ 2]=2*(*(arcx.ddisp+6*i+2))-2*initz[ 2];
			  if((arc.nodes+i)->code==412) venedisp[2][ 3]=2*(*(arcx.ddisp+6*i+2))-2*initz[ 3];
			  if((arc.nodes+i)->code==346) venedisp[2][ 4]=2*(*(arcx.ddisp+6*i+2))-2*initz[ 4];
			  if((arc.nodes+i)->code==348) venedisp[2][ 5]=2*(*(arcx.ddisp+6*i+2))-2*initz[ 5];
			  if((arc.nodes+i)->code==413) venedisp[2][ 6]=2*(*(arcx.ddisp+6*i+2))-2*initz[ 6];
			  if((arc.nodes+i)->code==268) venedisp[2][ 7]=2*(*(arcx.ddisp+6*i+2))-2*initz[ 7];
			  if((arc.nodes+i)->code==256) venedisp[2][ 8]=2*(*(arcx.ddisp+6*i+2))-2*initz[ 8];
			  if((arc.nodes+i)->code==112) venedisp[2][ 9]=2*(*(arcx.ddisp+6*i+2))-2*initz[ 9];
			  if((arc.nodes+i)->code==296) venedisp[2][10]=2*(*(arcx.ddisp+6*i+2))-2*initz[10];
			  if((arc.nodes+i)->code==295) venedisp[2][11]=2*(*(arcx.ddisp+6*i+2))-2*initz[11];
			}
			/*SRCAN*/
			n=strcspn((wdraw.childs+1)->sctfile,".");
			strncpy(str,(wdraw.childs+1)->sctfile,n);
			str[n]='\0';
			strcat(str,"-xn");
			i=srcan001(str);
			if(i==0) MessageBox(NULL,"Failed.","SRCAN001",MB_OK);
			strcpy((wdraw.childs+1)->otpfile,str);
			strcat((wdraw.childs+1)->otpfile,".rat");
			SetDlgItemText((wmenu.childs+2)->hwnd,ID_OUTPUTFILE,
						   (wdraw.childs+1)->otpfile);

			fout=fopen((wdraw.childs+1)->otpfile,"r");
			if(fout==NULL) break;
			readsrcanrate(fout,&arcx);
			for(i=0;i<arcx.nelem;i++)
			{
			  for(j=0;j<4;j++)
			  {
				if(*(srates+4*i+j)<(arcx.elems+i)->srate[j])
				{*(srates+4*i+j)=(arcx.elems+i)->srate[j];}
			  }
			}

		  /*SHORT Y POSI*/
			globalstatus=PINCHNODEVENEZIAYPOSI;
			strcpy(str,fname);
			strcat(str,".iyp");
			fin=fopen(str,"r");
			if(fin==NULL) break;
			initializeorganization(&globalorg);
			inputorganization(fin,&globalorg,&((wdraw.childs+1)->vparam));
			fclose(fin);
			str[n]='\0';
			for(i=0;i<(wdraw.childs+1)->org.nnode;i++)
			{
			  (globalorg.nodes+i)->d[GX]=((wdraw.childs+1)->org.nodes+i)->d[GX];
			  (globalorg.nodes+i)->d[GY]=((wdraw.childs+1)->org.nodes+i)->d[GY];
			  (globalorg.nodes+i)->d[GZ]=((wdraw.childs+1)->org.nodes+i)->d[GZ];
			}
			/*INITIAL*/
			arc =arci;
			arcx=arci;
			arcy=arci;
			free(globalorg.loads);
			for(i=0;i<globalorg.nelem;i++)
			{
			  for(j=0;j<2;j++)
			  {
				for(k=0;k<6;k++)
				{
				  (globalorg.elems+i)->initial[j][k]=0.0;
				}
			  }
			} /*UNAVAILABLE INITIAL CMQ.*/
			pnelem=globalorg.nelem;
			pelems=globalorg.elems;
			/*EXTRACT ARCLM*/
			extractarclmfromorgan(&globalorg,
								  &arc,&arcx,&arcy);
			globalorg.nelem=pnelem;
			globalorg.elems=pelems;
			/*(wmenu.childs+2)->vparam.vflag.mv.ftype=F_ARCLM;*/
			/*SAVE AS ARCLM*/
			strcpy(str,fname);
			strcat(str,"-yp.inl");
			strcpy((wdraw.childs+1)->inpfilez,str);
			SetDlgItemText((wmenu.childs+2)->hwnd,ID_INPUTFILEZ,str);
			strcpy(str,fname);
			strcat(str,"-yp.ihx");
			strcpy((wdraw.childs+1)->inpfilex,str);
			SetDlgItemText((wmenu.childs+2)->hwnd,ID_INPUTFILEX,str);
			strcpy(str,fname);
			strcat(str,"-yp.ihy");
			strcpy((wdraw.childs+1)->inpfiley,str);
			SetDlgItemText((wmenu.childs+2)->hwnd,ID_INPUTFILEY,str);
			saveasarclm((wdraw.childs+1)->inpfilez,&arc);
			saveasarclm((wdraw.childs+1)->inpfilex,&arcx);
			saveasarclm((wdraw.childs+1)->inpfiley,&arcy);
			/*ANALYSIS*/
			getviewparam((wmenu.childs+2)->hwnd,
						 &((wdraw.childs+1)->vparam));
			(wdraw.childs+1)->vparam.vflag.ev.deformation=0;
			strcpy(str,fname);
			strcat(str,".ypz");
			SetDlgItemText((wmenu.childs+2)->hwnd,ID_OUTPUTFILEZ,str);
			strcpy(str,fname);
			strcat(str,".ypx");
			SetDlgItemText((wmenu.childs+2)->hwnd,ID_OUTPUTFILEX,str);
			strcpy(str,fname);
			strcat(str,".ypy");
			SetDlgItemText((wmenu.childs+2)->hwnd,ID_OUTPUTFILEY,str);
			arclm001(&arc ,ID_INPUTFILEZ,ID_OUTPUTFILEZ);
			arclm001(&arcx,ID_INPUTFILEX,ID_OUTPUTFILEX);
			arclm001(&arcy,ID_INPUTFILEY,ID_OUTPUTFILEY);
			/*DISPLACEMENT*/
			for(i=0;i<arc.nnode;i++)
			{
			  if((arc.nodes+i)->code==111) venedisp[3][ 0]=2*(*(arcy.ddisp+6*i+2))-2*initz[ 0];
			  if((arc.nodes+i)->code==255) venedisp[3][ 1]=2*(*(arcy.ddisp+6*i+2))-2*initz[ 1];
			  if((arc.nodes+i)->code==267) venedisp[3][ 2]=2*(*(arcy.ddisp+6*i+2))-2*initz[ 2];
			  if((arc.nodes+i)->code==412) venedisp[3][ 3]=2*(*(arcy.ddisp+6*i+2))-2*initz[ 3];
			  if((arc.nodes+i)->code==346) venedisp[3][ 4]=2*(*(arcy.ddisp+6*i+2))-2*initz[ 4];
			  if((arc.nodes+i)->code==348) venedisp[3][ 5]=2*(*(arcy.ddisp+6*i+2))-2*initz[ 5];
			  if((arc.nodes+i)->code==413) venedisp[3][ 6]=2*(*(arcy.ddisp+6*i+2))-2*initz[ 6];
			  if((arc.nodes+i)->code==268) venedisp[3][ 7]=2*(*(arcy.ddisp+6*i+2))-2*initz[ 7];
			  if((arc.nodes+i)->code==256) venedisp[3][ 8]=2*(*(arcy.ddisp+6*i+2))-2*initz[ 8];
			  if((arc.nodes+i)->code==112) venedisp[3][ 9]=2*(*(arcy.ddisp+6*i+2))-2*initz[ 9];
			  if((arc.nodes+i)->code==296) venedisp[3][10]=2*(*(arcy.ddisp+6*i+2))-2*initz[10];
			  if((arc.nodes+i)->code==295) venedisp[3][11]=2*(*(arcy.ddisp+6*i+2))-2*initz[11];
			}
			/*SRCAN*/
			n=strcspn((wdraw.childs+1)->sctfile,".");
			strncpy(str,(wdraw.childs+1)->sctfile,n);
			str[n]='\0';
			strcat(str,"-yp");
			i=srcan001(str);
			if(i==0) MessageBox(NULL,"Failed.","SRCAN001",MB_OK);
			strcpy((wdraw.childs+1)->otpfile,str);
			strcat((wdraw.childs+1)->otpfile,".rat");
			SetDlgItemText((wmenu.childs+2)->hwnd,ID_OUTPUTFILE,
						   (wdraw.childs+1)->otpfile);

			fout=fopen((wdraw.childs+1)->otpfile,"r");
			if(fout==NULL) break;
			readsrcanrate(fout,&arcy);
			for(i=0;i<arcy.nelem;i++)
			{
			  for(j=0;j<4;j++)
			  {
				if(*(srates+4*i+j)<(arcy.elems+i)->srate[j])
				{*(srates+4*i+j)=(arcy.elems+i)->srate[j];}
			  }
			}

		  /*SHORT Y NEGA*/
			globalstatus=PINCHNODEVENEZIAYNEGA;
			strcpy(str,fname);
			strcat(str,".iyn");
			fin=fopen(str,"r");
			if(fin==NULL) break;
			initializeorganization(&globalorg);
			inputorganization(fin,&globalorg,&((wdraw.childs+1)->vparam));
			fclose(fin);
			str[n]='\0';
			for(i=0;i<(wdraw.childs+1)->org.nnode;i++)
			{
			  (globalorg.nodes+i)->d[GX]=((wdraw.childs+1)->org.nodes+i)->d[GX];
			  (globalorg.nodes+i)->d[GY]=((wdraw.childs+1)->org.nodes+i)->d[GY];
			  (globalorg.nodes+i)->d[GZ]=((wdraw.childs+1)->org.nodes+i)->d[GZ];
			}
			/*INITIAL*/
			arc =arci;
			arcx=arci;
			arcy=arci;
			free(globalorg.loads);
			for(i=0;i<globalorg.nelem;i++)
			{
			  for(j=0;j<2;j++)
			  {
				for(k=0;k<6;k++)
				{
				  (globalorg.elems+i)->initial[j][k]=0.0;
				}
			  }
			} /*UNAVAILABLE INITIAL CMQ.*/
			pnelem=globalorg.nelem;
			pelems=globalorg.elems;
			/*EXTRACT ARCLM*/
			extractarclmfromorgan(&globalorg,
								  &arc,&arcx,&arcy);
			globalorg.nelem=pnelem;
			globalorg.elems=pelems;
			/*(wmenu.childs+2)->vparam.vflag.mv.ftype=F_ARCLM;*/
			/*SAVE AS ARCLM*/
			strcpy(str,fname);
			strcat(str,"-yn.inl");
			strcpy((wdraw.childs+1)->inpfilez,str);
			SetDlgItemText((wmenu.childs+2)->hwnd,ID_INPUTFILEZ,str);
			strcpy(str,fname);
			strcat(str,"-yn.ihx");
			strcpy((wdraw.childs+1)->inpfilex,str);
			SetDlgItemText((wmenu.childs+2)->hwnd,ID_INPUTFILEX,str);
			strcpy(str,fname);
			strcat(str,"-yn.ihy");
			strcpy((wdraw.childs+1)->inpfiley,str);
			SetDlgItemText((wmenu.childs+2)->hwnd,ID_INPUTFILEY,str);
			saveasarclm((wdraw.childs+1)->inpfilez,&arc);
			saveasarclm((wdraw.childs+1)->inpfilex,&arcx);
			saveasarclm((wdraw.childs+1)->inpfiley,&arcy);
			/*ANALYSIS*/
			getviewparam((wmenu.childs+2)->hwnd,
						 &((wdraw.childs+1)->vparam));
			(wdraw.childs+1)->vparam.vflag.ev.deformation=0;
			strcpy(str,fname);
			strcat(str,".ynz");
			SetDlgItemText((wmenu.childs+2)->hwnd,ID_OUTPUTFILEZ,str);
			strcpy(str,fname);
			strcat(str,".ynx");
			SetDlgItemText((wmenu.childs+2)->hwnd,ID_OUTPUTFILEX,str);
			strcpy(str,fname);
			strcat(str,".yny");
			SetDlgItemText((wmenu.childs+2)->hwnd,ID_OUTPUTFILEY,str);
			arclm001(&arc ,ID_INPUTFILEZ,ID_OUTPUTFILEZ);
			arclm001(&arcx,ID_INPUTFILEX,ID_OUTPUTFILEX);
			arclm001(&arcy,ID_INPUTFILEY,ID_OUTPUTFILEY);
			/*DISPLACEMENT*/
			for(i=0;i<arc.nnode;i++)
			{
			  if((arc.nodes+i)->code==111) venedisp[4][ 0]=2*(*(arcy.ddisp+6*i+2))-2*initz[ 0];
			  if((arc.nodes+i)->code==255) venedisp[4][ 1]=2*(*(arcy.ddisp+6*i+2))-2*initz[ 1];
			  if((arc.nodes+i)->code==267) venedisp[4][ 2]=2*(*(arcy.ddisp+6*i+2))-2*initz[ 2];
			  if((arc.nodes+i)->code==412) venedisp[4][ 3]=2*(*(arcy.ddisp+6*i+2))-2*initz[ 3];
			  if((arc.nodes+i)->code==346) venedisp[4][ 4]=2*(*(arcy.ddisp+6*i+2))-2*initz[ 4];
			  if((arc.nodes+i)->code==348) venedisp[4][ 5]=2*(*(arcy.ddisp+6*i+2))-2*initz[ 5];
			  if((arc.nodes+i)->code==413) venedisp[4][ 6]=2*(*(arcy.ddisp+6*i+2))-2*initz[ 6];
			  if((arc.nodes+i)->code==268) venedisp[4][ 7]=2*(*(arcy.ddisp+6*i+2))-2*initz[ 7];
			  if((arc.nodes+i)->code==256) venedisp[4][ 8]=2*(*(arcy.ddisp+6*i+2))-2*initz[ 8];
			  if((arc.nodes+i)->code==112) venedisp[4][ 9]=2*(*(arcy.ddisp+6*i+2))-2*initz[ 9];
			  if((arc.nodes+i)->code==296) venedisp[4][10]=2*(*(arcy.ddisp+6*i+2))-2*initz[10];
			  if((arc.nodes+i)->code==295) venedisp[4][11]=2*(*(arcy.ddisp+6*i+2))-2*initz[11];
			}
			/*SRCAN*/
			n=strcspn((wdraw.childs+1)->sctfile,".");
			strncpy(str,(wdraw.childs+1)->sctfile,n);
			str[n]='\0';
			strcat(str,"-yn");
			i=srcan001(str);
			if(i==0) MessageBox(NULL,"Failed.","SRCAN001",MB_OK);
			strcpy((wdraw.childs+1)->otpfile,str);
			strcat((wdraw.childs+1)->otpfile,".rat");
			SetDlgItemText((wmenu.childs+2)->hwnd,ID_OUTPUTFILE,
						   (wdraw.childs+1)->otpfile);

			fout=fopen((wdraw.childs+1)->otpfile,"r");
			if(fout==NULL) break;
			readsrcanrate(fout,&arcy);
			fclose(fout);

			strcpy(str,fname);
			strcat(str,"-all.rat");
			fout=fopen(str,"w");
			if(fout==NULL) break;
			for(i=0;i<arcy.nelem;i++)
			{
			  for(j=0;j<4;j++)
			  {
				if(*(srates+4*i+j)<(arcy.elems+i)->srate[j])
				{*(srates+4*i+j)=(arcy.elems+i)->srate[j];}

				(arc.elems+i)->srate[j]=*(srates+4*i+j);
			  }

			  fprintf(fout,"ELEM: %5d SECT: %5d  %.5f %.5f %.5f %.5f\n",
					  (arc.elems+i)->code,
					  (arc.elems+i)->sect->code,
					  (arc.elems+i)->srate[0],
					  (arc.elems+i)->srate[1],
					  (arc.elems+i)->srate[2],
					  (arc.elems+i)->srate[3]);
			}
			fclose(fout);

			free(srates);

		  /*SAFETY FLAG*/
			globalstatus=PINCHNODEVENEZIASRCAN;
			(wdraw.childs+1)->vparam.vflag.ev.srcancolor=1;
			(wdraw.childs+1)->vparam.vflag.ev.srcanrate =1;

		  /*DRAW MODEL*/
			clearwindow(*(wdraw.childs+1));
			drawarclmframe((wdraw.childs+1)->hdcC,
						   (wdraw.childs+1)->vparam,arc,0,ONSCREEN);

			SetTextColor((wdraw.childs+1)->hdcC,RGB(255,255,255));/*WHITE*/
			for(i=0;i<arc.nnode;i++)
			{
			  if((arc.nodes+i)->code==111)
			  {
				tn=*(arc.nodes+i);
				vp=(wdraw.childs+1)->vparam;
				sprintf(str,"%6.3f",venedisp[4][ 0]*100.0); /*[cm]*/
				tn.d[GZ]+=vp.dparam.pitch;
				drawglobaltext((wdraw.childs+1)->hdcC,vp,tn,str);
				sprintf(str,"%6.3f",venedisp[1][ 0]*100.0); /*[cm]*/
				tn.d[GZ]+=vp.dparam.pitch;
				drawglobaltext((wdraw.childs+1)->hdcC,vp,tn,str);
			  }
			  if((arc.nodes+i)->code==255)
			  {
				tn=*(arc.nodes+i);
				vp=(wdraw.childs+1)->vparam;
				sprintf(str,"%6.3f",venedisp[1][ 1]*100.0); /*[cm]*/
				tn.d[GZ]+=vp.dparam.pitch;
				drawglobaltext((wdraw.childs+1)->hdcC,vp,tn,str);
			  }
			  if((arc.nodes+i)->code==267)
			  {
				tn=*(arc.nodes+i);
				vp=(wdraw.childs+1)->vparam;
				sprintf(str,"%6.3f",venedisp[1][ 2]*100.0); /*[cm]*/
				tn.d[GZ]+=vp.dparam.pitch;
				drawglobaltext((wdraw.childs+1)->hdcC,vp,tn,str);
			  }
			  if((arc.nodes+i)->code==412)
			  {
				tn=*(arc.nodes+i);
				vp=(wdraw.childs+1)->vparam;
				sprintf(str,"%6.3f",venedisp[3][ 3]*100.0); /*[cm]*/
				tn.d[GZ]+=vp.dparam.pitch;
				drawglobaltext((wdraw.childs+1)->hdcC,vp,tn,str);
				sprintf(str,"%6.3f",venedisp[1][ 3]*100.0); /*[cm]*/
				tn.d[GZ]+=vp.dparam.pitch;
				drawglobaltext((wdraw.childs+1)->hdcC,vp,tn,str);
			  }
			  if((arc.nodes+i)->code==346)
			  {
				tn=*(arc.nodes+i);
				vp=(wdraw.childs+1)->vparam;
				sprintf(str,"%6.3f",venedisp[3][ 4]*100.0); /*[cm]*/
				tn.d[GZ]+=vp.dparam.pitch;
				drawglobaltext((wdraw.childs+1)->hdcC,vp,tn,str);
			  }
			  if((arc.nodes+i)->code==348)
			  {
				tn=*(arc.nodes+i);
				vp=(wdraw.childs+1)->vparam;
				sprintf(str,"%6.3f",venedisp[3][ 5]*100.0); /*[cm]*/
				tn.d[GZ]+=vp.dparam.pitch;
				drawglobaltext((wdraw.childs+1)->hdcC,vp,tn,str);
			  }
			  if((arc.nodes+i)->code==413)
			  {
				tn=*(arc.nodes+i);
				vp=(wdraw.childs+1)->vparam;
				sprintf(str,"%6.3f",venedisp[3][ 6]*100.0); /*[cm]*/
				tn.d[GZ]+=vp.dparam.pitch;
				drawglobaltext((wdraw.childs+1)->hdcC,vp,tn,str);
				sprintf(str,"%6.3f",venedisp[2][ 6]*100.0); /*[cm]*/
				tn.d[GZ]+=vp.dparam.pitch;
				drawglobaltext((wdraw.childs+1)->hdcC,vp,tn,str);
			  }
			  if((arc.nodes+i)->code==268)
			  {
				tn=*(arc.nodes+i);
				vp=(wdraw.childs+1)->vparam;
				sprintf(str,"%6.3f",venedisp[2][ 7]*100.0); /*[cm]*/
				tn.d[GZ]+=vp.dparam.pitch;
				drawglobaltext((wdraw.childs+1)->hdcC,vp,tn,str);
			  }
			  if((arc.nodes+i)->code==256)
			  {
				tn=*(arc.nodes+i);
				vp=(wdraw.childs+1)->vparam;
				sprintf(str,"%6.3f",venedisp[2][ 8]*100.0); /*[cm]*/
				tn.d[GZ]+=vp.dparam.pitch;
				drawglobaltext((wdraw.childs+1)->hdcC,vp,tn,str);
			  }
			  if((arc.nodes+i)->code==112)
			  {
				tn=*(arc.nodes+i);
				vp=(wdraw.childs+1)->vparam;
				sprintf(str,"%6.3f",venedisp[4][ 9]*100.0); /*[cm]*/
				tn.d[GZ]+=vp.dparam.pitch;
				drawglobaltext((wdraw.childs+1)->hdcC,vp,tn,str);
				sprintf(str,"%6.3f",venedisp[2][ 9]*100.0); /*[cm]*/
				tn.d[GZ]+=vp.dparam.pitch;
				drawglobaltext((wdraw.childs+1)->hdcC,vp,tn,str);
			  }
			  if((arc.nodes+i)->code==296)
			  {
				tn=*(arc.nodes+i);
				vp=(wdraw.childs+1)->vparam;
				sprintf(str,"%6.3f",venedisp[4][10]*100.0); /*[cm]*/
				tn.d[GZ]+=vp.dparam.pitch;
				drawglobaltext((wdraw.childs+1)->hdcC,vp,tn,str);
			  }
			  if((arc.nodes+i)->code==295)
			  {
				tn=*(arc.nodes+i);
				vp=(wdraw.childs+1)->vparam;
				sprintf(str,"%6.3f",venedisp[4][11]*100.0); /*[cm]*/
				tn.d[GZ]+=vp.dparam.pitch;
				drawglobaltext((wdraw.childs+1)->hdcC,vp,tn,str);
			  }
			}

			/*draworganization((wdraw.childs+1)->hdcC,
							 (wdraw.childs+1)->vparam,
							 (wdraw.childs+1)->org,ONSCREEN);*/
			SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);
			/*MessageBox(NULL,"Completed.","Move Node with Safety",MB_OK);*/

			(wmenu.childs+2)->vparam.vflag.mv.ftype=F_ORGAN;
			SendMessage((wmenu.childs+2)->hwnd,WM_INITDIALOG,0,0);

			globalstatus=prestatus;
			globalmessageflag=1;
		  }
#endif
/************/
		  else
		  {
			clearwindow(*(wdraw.childs+1));
			draworganization((wdraw.childs+1)->hdcC,
							 (wdraw.childs+1)->vparam,
							 (wdraw.childs+1)->org,ONSCREEN);
			SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);

			globalstatus=prestatus;
		  }

		  /*globalstatus=prestatus;*/
		}

		else if(globalstatus==SELECTMULTIELEMBEGIN &&
				wParam == MK_LBUTTON)
		{
		  point.x = LOWORD(lParam);
		  point.y = HIWORD(lParam);

		  chaserect(initx,point.x,inity,point.y,
					*(wdraw.childs+1));
		}
		else if(globalstatus==SELECTMULTIELEMBEGIN) /*END*/
		{
		  point.x = LOWORD(lParam);
		  point.y = HIWORD(lParam);

	     if((wmenu.childs+2)->vparam.vflag.mv.ftype==F_ORGAN)
         {
		  selectmultielem(&((wdraw.childs+1)->org),
						  (wdraw.childs+1)->vparam,
						  initx,point.x,inity,point.y);

		  /*REDRAW FRAME*/
		  clearwindow(*(wdraw.childs+1));
		  draworganization((wdraw.childs+1)->hdcC,
						   (wdraw.childs+1)->vparam,
						   (wdraw.childs+1)->org,ONSCREEN);

		  /*REDRAW SELECTED ELEM*/
		  hdc = (wdraw.childs+1)->hdcC;
		  setfontformat(hdc,15,6,"MS Mincho",255,255,255);
		  for(i=0;i<nmultielem;i++)
		  {
			for(j=0;j<((*(multielem+i))->nnod);j++)
			{
			  nodeontoscreen(**((*(multielem+i))->nods+j),&Nx,&Ny,
							 (wdraw.childs+1)->vparam);
			  sprintf(str,"%ld",
					  (*((*(multielem+i))->nods+j))->code);
			  TextOut(hdc,Nx+2,Ny,str,strlen(str));
			}
		  }
		  ReleaseDC((wdraw.childs+1)->hwnd,hdc);

		  SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);

		  globalstatus=prestatus;
         }
        }
        else if(globalstatus==SELECTMULTINODEBEGIN &&
				wParam == MK_LBUTTON)
		{
		  point.x = LOWORD(lParam);
		  point.y = HIWORD(lParam);

		  chaserect(initx,point.x,inity,point.y,
					*(wdraw.childs+1));
		}
		else if(globalstatus==SELECTMULTINODEBEGIN) /*END*/
		{
		  point.x = LOWORD(lParam);
		  point.y = HIWORD(lParam);

		  selectmultinode(&((wdraw.childs+1)->org),
						  (wdraw.childs+1)->vparam,
						  initx,point.x,inity,point.y);

		  /*REDRAW FRAME*/
		  clearwindow(*(wdraw.childs+1));
		  draworganization((wdraw.childs+1)->hdcC,
						   (wdraw.childs+1)->vparam,
						   (wdraw.childs+1)->org,ONSCREEN);

		  /*REDRAW SELECTED NODE*/
		  hdc = (wdraw.childs+1)->hdcC;
		  setfontformat(hdc,15,6,"MS Mincho",255,255,255);
		  for(i=0;i<nmultinode;i++)
		  {
			  nodeontoscreen(**(multinode+i),&Nx,&Ny,
							 (wdraw.childs+1)->vparam);
			  sprintf(str,"%ld",(*(multinode+i))->code);
			  TextOut(hdc,Nx+2,Ny,str,strlen(str));
		  }
		  ReleaseDC((wdraw.childs+1)->hwnd,hdc);

		  SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);

		  globalstatus=prestatus;
        }
      }

      if((wmenu.childs+2)->vparam.vflag.mv.ftype==F_ARCLM)
      /*UJIOKA FOR BCLNG003*/
      {
        if(globalstatus==SELECTMULTIELEMBEGIN &&
				wParam == MK_LBUTTON)
		{
		  point.x = LOWORD(lParam);
		  point.y = HIWORD(lParam);

		  chaserect(initx,point.x,inity,point.y,
					*(wdraw.childs+1));
		}
        else if((globalstatus==SELECTMULTIELEMBEGIN)) /*END*/
        {
		  point.x = LOWORD(lParam);
		  point.y = HIWORD(lParam);

          if((wmenu.childs+2)->vparam.vflag.mv.ftype==F_ARCLM)
          {
		  multiwire=selectmultiwire(&arc,
						  (wdraw.childs+1)->vparam,
						  initx,point.x,inity,point.y);
		  clearwindow(*(wdraw.childs+1));
          drawarclmframe((wdraw.childs+1)->hdcC,
                           (wdraw.childs+1)->vparam,arc,0,ONSCREEN);

		  /*REDRAW SELECTED ELEM*/
		  hdc = (wdraw.childs+1)->hdcC;
		  setfontformat(hdc,15,6,"MS Mincho",255,255,255);


//          sprintf(str,"nmultiwire: %5d \n",nmultiwire);
//          errormessage(str);
/*
          for(ii=0;ii<nmultiwire;ii++)
          {
           sprintf(str,"MULTIWIRE[%ld]:CODE=%ld",ii+1,(*(multiwire+ii))->code);
		   errormessage(str);
          }
*/

		  for(i=0;i<nmultiwire;i++)
		  {
			for(j=0;j<2;j++)
			{
              nodeontoscreen((*(*(multiwire+i))->node[j]),&Nx,&Ny,
			  				 (wdraw.childs+1)->vparam);
			  sprintf(str,"%ld",
					  (((*(multiwire+i))->node[j]))->code);
			  TextOut(hdc,Nx+2,Ny,str,strlen(str));
/*
              sprintf(str,"multiwire %ld:NODE %ld loff=%ld",(*(multiwire+i))->code,
               (((*(multiwire+i))->node[j]))->code,(((*(multiwire+i))->node[j]))->loff);
              errormessage(str);
*/
			}
		  }
		  ReleaseDC((wdraw.childs+1)->hwnd,hdc);
		  SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);
		  globalstatus=prestatus;
          }
        }
	  }

	  /*CASE NOT CREATING.*/
	  if(globalstatus==MOVE &&
		 (wmenu.childs+2)->vparam.vflag.mv.ftype==F_SRCAN &&
		 (wParam == MK_LBUTTON))
	  {
		point.x = LOWORD(lParam);
		point.y = HIWORD(lParam);

		x=initx-point.x;
		y=inity-point.y;
		initx=point.x;
		inity=point.y;

		sprintf(str,"%d",(wdraw.childs+1)->vparam.Xo-x);
		SetDlgItemText((wmenu.childs+2)->hwnd,IDV_ORIGINX,str);
        sprintf(str,"%d",(wdraw.childs+1)->vparam.Yo-y);
        SetDlgItemText((wmenu.childs+2)->hwnd,IDV_ORIGINY,str);

		wparam = MAKEWPARAM((WORD)IDM_DRAWSECTION,(WORD)0);
        SendMessage(wmain.hwnd,WM_COMMAND,wparam,(WPARAM)0);
      }
//Modified by Fukushima 2014/4/30///////////////////////////////////////////////
		//else if (globalstatus == ROTATE && (wdraw.childs + 1)
		//	->lstatus == ROTATE && (wParam == MK_LBUTTON))
		else if((globalstatus==ROTATE &&
				(wdraw.childs+1)->lstatus==ROTATE &&
				(wParam == MK_LBUTTON)) ||
				(globalstatus==MOVE &&
				(wdraw.childs+1)->lstatus==MOVE &&
				(wParam == MK_MBUTTON)))
//Modified by Fukushima 2014/4/30///////////////////////////////////////////////
      {
        point.x = LOWORD(lParam);
        point.y = HIWORD(lParam);

        x=initx-point.x;
        y=inity-point.y;

        if(labs(x)>=labs(y))
        {
		  (wdraw.childs+1)->vparam.theta=initv.theta
                                         +(double)x/10.0;
          sprintf(str,"%.1f",(wdraw.childs+1)->vparam.theta);
          SetDlgItemText((wmenu.childs+2)->hwnd,IDV_THETA,str);
        }
        else
        {
          (wdraw.childs+1)->vparam.phi=initv.phi
                                       -(double)y/10.0;
          sprintf(str,"%.1f",(wdraw.childs+1)->vparam.phi);
          SetDlgItemText((wmenu.childs+2)->hwnd,IDV_PHI,str);
        }

        clearwindow(*(wdraw.childs+1));
        if((wmenu.childs+2)->vparam.vflag.mv.ftype==F_ARCLM ||
           (wmenu.childs+2)->vparam.vflag.mv.ftype==F_FRAME)
        {
          drawrotatingarclm((wdraw.childs+1)->hdcC,
                            (wdraw.childs+1)->vparam,arc);
        }
        if((wmenu.childs+2)->vparam.vflag.mv.ftype==F_ORGAN)
        {
          draworgannodes((wdraw.childs+1)->hdcC,
                         (wdraw.childs+1)->vparam,
						 (wdraw.childs+1)->org);
		}
		if((wmenu.childs+2)->vparam.vflag.mv.ftype==F_BQUAD)
		{
		  /*drawbiquadframe((wdraw.childs+1)->hdcC,
						  (wdraw.childs+1)->vparam,bqf,ONSCREEN);*/
		  drawbiquadnodes((wdraw.childs+1)->hdcC,
						  (wdraw.childs+1)->vparam,bqf);
		}
		SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);
	  }
//Modified by Fukushima 2014/4/30///////////////////////////////////////////////
		//else if (globalstatus == MOVE && (wdraw.childs + 1)->lstatus == MOVE &&
		//	(wParam == MK_LBUTTON))
        else if((globalstatus==MOVE &&
				(wdraw.childs+1)->lstatus==MOVE &&
				(wParam == MK_LBUTTON)) ||
				(globalstatus==ROTATE &&
				(wdraw.childs+1)->lstatus==ROTATE &&
		        (wParam == MK_MBUTTON)))
//Modified by Fukushima 2014/4/30///////////////////////////////////////////////
	  {
		point.x = LOWORD(lParam);
		point.y = HIWORD(lParam);

        x=initx-point.x;
        y=inity-point.y;

        direction=getdirection((wdraw.childs+1)->vparam,
                               x,y,initv.focus,&length);
        if(direction == DIRECTIONX)
        {
          (wdraw.childs+1)->vparam.focus.d[GX]=initv.focus.d[GX]
                                               +length;
          sprintf(str,"%.1f",(wdraw.childs+1)->vparam.focus.d[GX]);
          SetDlgItemText((wmenu.childs+2)->hwnd,IDV_X,str);
        }
        else if(direction == DIRECTIONY)
        {
          (wdraw.childs+1)->vparam.focus.d[GY]=initv.focus.d[GY]
                                               +length;
          sprintf(str,"%.1f",(wdraw.childs+1)->vparam.focus.d[GY]);
          SetDlgItemText((wmenu.childs+2)->hwnd,IDV_Y,str);
        }
        else if(direction == DIRECTIONZ)
        {
          (wdraw.childs+1)->vparam.focus.d[GZ]=initv.focus.d[GZ]
                                               +length;
          sprintf(str,"%.1f",(wdraw.childs+1)->vparam.focus.d[GZ]);
          SetDlgItemText((wmenu.childs+2)->hwnd,IDV_Z,str);
        }

        clearwindow(*(wdraw.childs+1));
        if((wmenu.childs+2)->vparam.vflag.mv.ftype==F_ARCLM ||
           (wmenu.childs+2)->vparam.vflag.mv.ftype==F_FRAME)
		{
          drawrotatingarclm((wdraw.childs+1)->hdcC,
                            (wdraw.childs+1)->vparam,arc);
        }
        if((wmenu.childs+2)->vparam.vflag.mv.ftype==F_ORGAN)
		{
          draworgannodes((wdraw.childs+1)->hdcC,
                         (wdraw.childs+1)->vparam,
                         (wdraw.childs+1)->org);
        }
        SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);
      }
      else if((globalstatus==ROTATE &&
               (wdraw.childs+1)->lstatus==ROTATE) ||
              (globalstatus==MOVE &&
               (wdraw.childs+1)->lstatus==MOVE))     /*END OF MOVE.*/
      {
        globalstatus=prestatus; /*ROTATION END*/

        if(prestatus==SELECTNODE)
        {
          GetDlgItemText((wmenu.childs+3)->hwnd,IDN_CODE,str,80);
          code=strtol(str,NULL,10);
        }
        else if(prestatus==SELECTSECTION)
        {
          GetDlgItemText((wmenu.childs+3)->hwnd,IDS_CODE,str,80);
          code=strtol(str,NULL,10);
        }
        else if(prestatus==SELECTELEMENT)
        {
          GetDlgItemText((wmenu.childs+3)->hwnd,IDE_CODE,str,80);
          code=strtol(str,NULL,10);
        }
        else code=0;

        clearwindow(*(wdraw.childs+1));
        if((wmenu.childs+2)->vparam.vflag.mv.ftype==F_ARCLM ||
           (wmenu.childs+2)->vparam.vflag.mv.ftype==F_FRAME)
        {
          drawarclmframe((wdraw.childs+1)->hdcC,
                         (wdraw.childs+1)->vparam,
                         arc,code,ONSCREEN);
		}
		if((wmenu.childs+2)->vparam.vflag.mv.ftype==F_ORGAN)
		{
		  draworganization((wdraw.childs+1)->hdcC,
						   (wdraw.childs+1)->vparam,
						   (wdraw.childs+1)->org,ONSCREEN);
		}
		if((wmenu.childs+2)->vparam.vflag.mv.ftype==F_BQUAD)
		{
		  drawbiquadframe((wdraw.childs+1)->hdcC,
						  (wdraw.childs+1)->vparam,bqf,ONSCREEN);
		}
		SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);

        /*globalstatus=NEUTRAL;*/ /*ROTATION END*/
      }
      else if(globalstatus==ROTATE || globalstatus==MOVE)
      {
        /*globalstatus=NEUTRAL;*/ /*ROTATION END*/
        globalstatus=prestatus; /*ROTATION END*/
      }

      break;

    case WM_COMMAND:
      switch(LOWORD(wParam))
      {
        case IDM_FITSHEET: /*NEVER USE FITPARENT.*/
          hoya=GetParent(GetParent(hwnd)); /*PARENT SHEET.*/

          getclientsize(wmain.hwnd,&maxX,&maxY);
          getwindowsize(wmenu.hwnd,&mw,&mh);

          MoveWindow(hoya,mw,0,(maxX-mw),maxY,TRUE);
          break;

        case IDM_HOGOPEN:
          wparam = MAKEWPARAM((WORD)IDD_VIEW,(WORD)0);
          SendMessage((wmenu.childs+2)->hwnd,WM_COMMAND,
                      wparam,(WPARAM)0);
          break;

        case IDM_POPOPENSECTION:
          wparam = MAKEWPARAM((WORD)IDM_OPENSECTION,(WORD)0);
          SendMessage(wmain.hwnd,WM_COMMAND,wparam,(WPARAM)0);
          break;

        case IDM_POPTEXT:
          hpopdlg=CreateDialog(hInstGlobal,
                               "CANTEXTDLG",
                               NULL,
                               (DLGPROC)DialogProcText);
          getwindowsize(hpopdlg,&mw,&mh);
          point.x=initx;
          point.y=inity;
          ClientToScreen((wdraw.childs+1)->hwnd,&point);

          MoveWindow(hpopdlg,point.x,point.y,mw,mh,TRUE);

          idtext=findtext((wdraw.childs+1)->hdcC,
                          wdraw.strset,wdraw.nstring,
                          initx-((wdraw.childs+1)->vparam.Xo),
                          inity-((wdraw.childs+1)->vparam.Yo));
          if(idtext)
          {
            SetDlgItemText(hpopdlg,IDT_POPTEXT,
                           (wdraw.strset+idtext-1)->str);
            sprintf(str,"%.0f",(wdraw.strset+idtext-1)->n.d[0]);
            SetDlgItemText(hpopdlg,IDT_POPX,str);
            sprintf(str,"%.0f",(wdraw.strset+idtext-1)->n.d[1]);
            SetDlgItemText(hpopdlg,IDT_POPY,str);
          }
          else
          {
            sprintf(str,"%ld",initx-((wdraw.childs+1)->vparam.Xo));
            SetDlgItemText(hpopdlg,IDT_POPX,str);
            sprintf(str,"%ld",inity-((wdraw.childs+1)->vparam.Yo));
            SetDlgItemText(hpopdlg,IDT_POPY,str);
          }

          ShowWindow(hpopdlg,SW_SHOW);
          break;
        case IDM_POPTEXTRETURN:
          GetDlgItemText(hpopdlg,IDT_POPX,str,256);
          dx=strtod(str,NULL);
          GetDlgItemText(hpopdlg,IDT_POPY,str,256);
          dy=strtod(str,NULL);
          GetDlgItemText(hpopdlg,IDT_POPTEXT,str,256);

          if(idtext)
          {
            strcpy((wdraw.strset+idtext-1)->str,str);
            (wdraw.strset+idtext-1)->n.d[0]=dx;
            (wdraw.strset+idtext-1)->n.d[1]=dy;
          }
          else if(strcmp(str,"\0"))
          {
            wdraw.nstring++;
            wdraw.strset=(struct snode *)
                         realloc(wdraw.strset,
                                 wdraw.nstring*sizeof(struct snode));

            strcpy((wdraw.strset+wdraw.nstring-1)->str,str);
            (wdraw.strset+wdraw.nstring-1)->n.d[0]=dx;
            (wdraw.strset+wdraw.nstring-1)->n.d[1]=dy;
          }
          SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);
          break;

        case IDM_POPCHOOSENODE: /*SELECT NODE.*/
          if((wmenu.childs+2)->vparam.vflag.mv.ftype==F_ORGAN)
          {
            if(wdraw.hwnd!=NULL &&
               (wdraw.childs+1)->org.nodes!=NULL)
            {
              globalstatus = SELECTNODE;

              hcursor = LoadCursor(hInstGlobal,"CANBOXW");
              SetClassLong((wdraw.childs+1)->hwnd,
                           GCL_HCURSOR,(LONG)hcursor);
            }
          }
          break;
        case IDM_POPPERPENDICULAR: /*SELECT PERPENDICULAR DOT.*/
          if((wmenu.childs+2)->vparam.vflag.mv.ftype==F_ORGAN)
          {
            if(wdraw.hwnd!=NULL &&
               (wdraw.childs+1)->org.nodes!=NULL)
            {
              globalstatus = SELECTPERPENDICULAR;

              hcursor = LoadCursor(hInstGlobal,"CANBOXW");
              SetClassLong((wdraw.childs+1)->hwnd,
                           GCL_HCURSOR,(LONG)hcursor);
            }
          }
          break;
        case IDM_POPDECIDEELEM:
          /*END OF ADDING ELEMENT.*/
          if((wmenu.childs+2)->vparam.vflag.mv.ftype==F_ORGAN)
          {
			if(createcommand==C_ADD ||
			   globalstatus==ADDELEMWITHSAFETY)
			{
			  if(createitem==C_ELEMENT && gelem.nnod>0)
			  {
				gelem.code=((wdraw.childs+1)->org.elems
							+(wdraw.childs+1)->org.nelem-1)->code+1;
				gelem.loff=(wdraw.childs+1)->org.nelem; /*LAST.*/
				gelem.sect=(wdraw.childs+1)->org.sects+0;

				gelem.nnod--;
				gelem.nods=(struct onode **)
						   realloc(gelem.nods,
								   gelem.nnod
								   *sizeof(struct onode *));
				gelem.bonds=(signed char *)
							malloc(6*(gelem.nnod)
								   *sizeof(signed char));
				if(gelem.nods==NULL || gelem.bonds==NULL)
				{
				  MessageBox(NULL,"Buffer Null.","Adding",MB_OK);
				}

				for(i=0;i<6*(gelem.nnod);i++) *(gelem.bonds+i)=0;

				if(gelem.nban>0)
				{
				  (gelem.bans+0)->nnod--;

				  if(((gelem.bans+0)->nnod)<=2)
				  {
					gelem.nban=0;
					free((gelem.bans+0)->nods);
					free(gelem.bans);
				  }
				  else
				  {
					(gelem.bans+0)->nods
					=(struct onode **)
					  realloc((gelem.bans+0)->nods,
							  (gelem.bans+0)->nnod
							  *sizeof(struct onode *));
					if((gelem.bans+0)->nods==NULL)
					{
					  MessageBox(NULL,"Buffer Null.","Adding",MB_OK);
					}
				  }
				}

				addelement(&gelem,&((wdraw.childs+1)->org));

				/*freeelement(&gelem,0,0,1,0);*/

/*MessageBox(NULL,"Passed.","Adding Element",MB_OK);*/
				gelem.nnod=0;
				gelem.nban=0;
				gelem.nods=NULL;
				gelem.bonds=NULL;
				gelem.bans=NULL;
			  }

if(globalstatus==ADDELEMWITHSAFETY)
{
			/*globalstatus=prestatus;*/
			globalmessageflag=0;

			/*INITIAL*/
			arc =arci;
			arcx=arci;
			arcy=arci;
			free((wdraw.childs+1)->org.loads);
			for(i=0;i<(wdraw.childs+1)->org.nelem;i++)
			{
			  for(j=0;j<2;j++)
			  {
				for(k=0;k<6;k++)
				{
				  ((wdraw.childs+1)->org.elems+i)->initial[j][k]=0.0;
				}
			  }
			} /*UNAVAILABLE INITIAL CMQ.*/

			pnelem=(wdraw.childs+1)->org.nelem;
			pelems=(wdraw.childs+1)->org.elems;

			/*EXTRACT ARCLM*/
			extractarclmfromorgan(&((wdraw.childs+1)->org),
								  &arc,&arcx,&arcy);

			(wdraw.childs+1)->org.nelem=pnelem;
			(wdraw.childs+1)->org.elems=pelems;

			(wmenu.childs+2)->vparam.vflag.mv.ftype=F_ARCLM;

			/*SAVE AS ARCLM*/
			saveasarclm((wdraw.childs+1)->inpfilez,&arc);
			saveasarclm((wdraw.childs+1)->inpfilex,&arcx);
			saveasarclm((wdraw.childs+1)->inpfiley,&arcy);

			/*ANALYSIS*/
			getviewparam((wmenu.childs+2)->hwnd,
						 &((wdraw.childs+1)->vparam));
			(wdraw.childs+1)->vparam.vflag.ev.deformation=0;

			arclm001(&arc ,ID_INPUTFILEZ,ID_OUTPUTFILEZ);
			arclm001(&arcx,ID_INPUTFILEX,ID_OUTPUTFILEX);
			arclm001(&arcy,ID_INPUTFILEY,ID_OUTPUTFILEY);

			/*SRCAN*/
			n=strcspn((wdraw.childs+1)->sctfile,".");
			strncpy(str,(wdraw.childs+1)->sctfile,n);
			str[n]='\0';

			i=srcan001(str);

			if(i==0) MessageBox(NULL,"Failed.","SRCAN001",MB_OK);

			strcpy((wdraw.childs+1)->otpfile,str);
			strcat((wdraw.childs+1)->otpfile,".rat");
			SetDlgItemText((wmenu.childs+2)->hwnd,ID_OUTPUTFILE,
						   (wdraw.childs+1)->otpfile);

			fout=fopen((wdraw.childs+1)->otpfile,"r");
			if(fout==NULL) break;
			readsrcanrate(fout,&arc);

			/*SAFETY FLAG*/
			(wdraw.childs+1)->vparam.vflag.ev.srcancolor=1;

			clearwindow(*(wdraw.childs+1));
			drawarclmframe((wdraw.childs+1)->hdcC,
						   (wdraw.childs+1)->vparam,arc,0,ONSCREEN);

			/*draworganization((wdraw.childs+1)->hdcC,
							 (wdraw.childs+1)->vparam,
							 (wdraw.childs+1)->org,ONSCREEN);*/
			SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);
			/*MessageBox(NULL,"Completed.","Move Node with Safety",MB_OK);*/

			(wmenu.childs+2)->vparam.vflag.mv.ftype=F_ORGAN;
			SendMessage((wmenu.childs+2)->hwnd,WM_INITDIALOG,0,0);

			globalmessageflag=1;

}
else
{
			  draworganization((wdraw.childs+1)->hdcC,
							   (wdraw.childs+1)->vparam,
							   (wdraw.childs+1)->org,ONSCREEN);
			  SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);

			  /*hcursor = LoadCursor(hInstGlobal,"CANCURSORW");
			  SetClassLong((wdraw.childs+1)->hwnd,
						   GCL_HCURSOR,(LONG)hcursor);
			  globalstatus=NEUTRAL;*/
}
			  break;
			}
		  }
		  break;

		case IDM_POPSECTIONPAINT:
		case IDM_POPCOPYBOND:
		  if(wdraw.nchilds>=2 &&
			 wmenu.nchilds>=3 &&
			 (wmenu.childs+2)->vparam.vflag.mv.ftype==F_ORGAN)
		  {
			point.x=initx;
			point.y=inity;

			foundelem=selectorganelement((wdraw.childs+1)->vparam,
										 &((wdraw.childs+1)->org),
										 point,
										 &selectedline,
										 &selectedban);
			if(foundelem==NULL)
			{
			  MessageBox(hwnd,"Nothing.","Element",MB_OK);
			  break;
			}

			/*DRAW NODE OF SELECTED ELEM.*/
			/*hdc = GetDC((wdraw.childs+1)->hwnd);*/
			hdc = (wdraw.childs+1)->hdcC;
			setfontformat(hdc,15,6,"MS Mincho",255,255,255);
			for(i=0;i<(foundelem->nnod);i++)
			{
			  nodeontoscreen(**(foundelem->nods+i),&Nx,&Ny,
							 (wdraw.childs+1)->vparam);
			  sprintf(str,"%ld",(*(foundelem->nods+i))->code);
			  TextOut(hdc,Nx+2,Ny,str,strlen(str));
			}
			ReleaseDC((wdraw.childs+1)->hwnd,hdc);
			SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);

			if(LOWORD(wParam)==IDM_POPSECTIONPAINT)
			{
			  gelem.sect=foundelem->sect;
			  gelem.type=foundelem->type;
			  gelem.role=foundelem->role;
			  globalstatus=SECTIONPAINT;
			}
			if(LOWORD(wParam)==IDM_POPCOPYBOND)
			{
			  gelem.bonds=foundelem->bonds;
			  globalstatus=COPYBOND;
			}
		  }
		  break;

		case IDM_POPDIVIDEELEMENT:
		  if(wdraw.nchilds>=2 &&
			 wmenu.nchilds>=3 &&
			 (wmenu.childs+2)->vparam.vflag.mv.ftype==F_ORGAN)
		  {
			globalstatus=DIVIDEELEMENT;
		  }
		  break;

        case IDM_POPDIVIDEALLELEMENT:
		  if(wdraw.nchilds>=2 &&
			 wmenu.nchilds>=3 &&
			 (wmenu.childs+2)->vparam.vflag.mv.ftype==F_ORGAN)
		  {
			divideallelem(&((wdraw.childs+1)->org));
		  }
		  break;

		case IDM_POPSECTIONPAINTEND:
		case IDM_POPDIVIDEELEMENTEND:
		  if(globalstatus==SECTIONPAINT || globalstatus==DIVIDEELEMENT)
		  {
			gelem.type=0;
			gelem.role=ROLENULL;

			globalstatus=NEUTRAL;
		  }
		  break;

		case IDM_POPMOVENODEMOUSE:
		  if(wdraw.nchilds>=2 &&
			 wmenu.nchilds>=3 &&
			 (wmenu.childs+2)->vparam.vflag.mv.ftype==F_ORGAN)
		  {
			hcursor = LoadCursor(hInstGlobal,"CANBOXW");
			SetClassLong((wdraw.childs+1)->hwnd,
                         GCL_HCURSOR,(LONG)hcursor);

            globalstatus=MOVENODEBYMOUSE;
          }
          break;

        case IDM_POPPINCHNODEMOUSE:
          if(wdraw.nchilds>=2 &&
             wmenu.nchilds>=3 &&
             (wmenu.childs+2)->vparam.vflag.mv.ftype==F_ORGAN)
          {
            hcursor = LoadCursor(hInstGlobal,"CANBOXW");
            SetClassLong((wdraw.childs+1)->hwnd,
                         GCL_HCURSOR,(LONG)hcursor);

            globalstatus=PINCHNODEBYMOUSE;
          }
          break;

        case IDM_POPPINCHNODEWITHSAFETY:
          if(wdraw.nchilds>=2 &&
             wmenu.nchilds>=3 &&
             (wmenu.childs+2)->vparam.vflag.mv.ftype==F_ORGAN)
          {
            hcursor = LoadCursor(hInstGlobal,"CANBOXW");
            SetClassLong((wdraw.childs+1)->hwnd,
                         GCL_HCURSOR,(LONG)hcursor);
			if(MessageBox(NULL,"Initial Analysis","Pinch Node",
						  MB_OKCANCEL)==IDOK)
			{
            arclmautomatic(&((wdraw.childs+1)->org));
            (wdraw.childs+1)->vparam.vflag.ev.srcancolor=1;
                        fout=fopen((wdraw.childs+1)->otpfile,"r");
            if(fout==NULL) break;
            readsrcanrate(fout,&arc);
            fclose(fout);
			int nred=0;
            double maxsafety=0.0;
            int ii,jj;
            for(ii=0;ii<arc.nelem;ii++)
            {
              //if((arc.elems+ii)->srate[2]>1) nred=nred+1;
              if( (arc.elems+ii)->srate[0]>1 ||
                  (arc.elems+ii)->srate[1]>1 ||
                  (arc.elems+ii)->srate[2]>1 ||
                  (arc.elems+ii)->srate[3]>1 )
                {
                  nred=nred+1;
                }

				for(jj=0;jj<4;jj++)
				{
				  if(maxsafety<(arc.elems+ii)->srate[jj]) maxsafety=(arc.elems+ii)->srate[jj];
				}

            }

            sprintf(str,"%d",nred);
			SetDlgItemText((wmenu.childs+2)->hwnd,ID_SAFETY,str);
			SendDlgItemMessage((wmenu.childs+2)->hwnd,ID_SAFETY,WM_PAINT,0,0);
			}

			globalstatus=PINCHNODEWITHSAFETY;
         }
          break;

		case IDM_POPPINCHNODEVENEZIA:
		  if(wdraw.nchilds>=2 &&
			 wmenu.nchilds>=3 &&
			 (wmenu.childs+2)->vparam.vflag.mv.ftype==F_ORGAN)
		  {
			hcursor = LoadCursor(hInstGlobal,"CANBOXW");
			SetClassLong((wdraw.childs+1)->hwnd,
						 GCL_HCURSOR,(LONG)hcursor);

			globalstatus=PINCHNODEVENEZIA;
		  }
		  break;

		case IDM_POPPINCHNODESTACK:
		  if(wdraw.nchilds>=2 &&
			 wmenu.nchilds>=3 &&
			 (wmenu.childs+2)->vparam.vflag.mv.ftype==F_ORGAN)
		  {
			hcursor = LoadCursor(hInstGlobal,"CANBOXW");
			SetClassLong((wdraw.childs+1)->hwnd,
						 GCL_HCURSOR,(LONG)hcursor);

			globalstatus=PINCHNODESTACK;
		  }
		  break;

		case IDM_POPPINCHELEMTWIST:
		  if(wdraw.nchilds>=2 &&
			 wmenu.nchilds>=3 &&
			 (wmenu.childs+2)->vparam.vflag.mv.ftype==F_ORGAN)
		  {
			globalstatus=PINCHELEMTWIST;
		  }
		  break;

		case IDM_POPPINCHNODEELASTIC:
		  if(wdraw.nchilds>=2 &&
			 wmenu.nchilds>=3 &&
			 (wmenu.childs+2)->vparam.vflag.mv.ftype==F_ORGAN)
		  {
			hcursor = LoadCursor(hInstGlobal,"CANBOXW");
			SetClassLong((wdraw.childs+1)->hwnd,
						 GCL_HCURSOR,(LONG)hcursor);

			globalstatus=PINCHNODEELASTIC;
		  }
		  break;

		case IDM_POPADDELEMWITHSAFETY:
		  if(wdraw.nchilds>=2 &&
			 wmenu.nchilds>=3 &&
			 (wmenu.childs+2)->vparam.vflag.mv.ftype==F_ORGAN)
		  {
			hcursor = LoadCursor(hInstGlobal,"CANBOXW");
			SetClassLong((wdraw.childs+1)->hwnd,
						 GCL_HCURSOR,(LONG)hcursor);

			globalstatus=ADDELEMWITHSAFETY;
			createitem=C_ELEMENT;
		  }
		  break;

		case IDM_POPADDCOLUMNWITHSAFETY:
		  if(wdraw.nchilds>=2 &&
			 wmenu.nchilds>=3 &&
			 (wmenu.childs+2)->vparam.vflag.mv.ftype==F_ORGAN)
		  {
			/*hcursor = LoadCursor(hInstGlobal,"CANBOXW");*/
			/*SetClassLong((wdraw.childs+1)->hwnd,
						 GCL_HCURSOR,(LONG)hcursor);*/

			/*globalstatus=ADDCOLUMNWITHSAFETY;*/
			/*createitem=C_ELEMENT;*/

			node=selectorgannode((wdraw.childs+1)->vparam,
								 &((wdraw.childs+1)->org),point);
			if(node==NULL) break;
			loff=node->loff;

			gelem.code=((wdraw.childs+1)->org.elems
						+((wdraw.childs+1)->org.nelem-1))->code+1;
			gelem.loff=(wdraw.childs+1)->org.nelem;
			gelem.nnod=2;
			gelem.nods=(struct onode **)
						malloc(gelem.nnod*sizeof(struct onode *));
			gelem.bonds=(signed char *)
						 malloc(6*(gelem.nnod)*sizeof(signed char));
			for(i=0;i<12;i++) *(gelem.bonds+i)=0;
			gelem.type=COLUMN;
			gelem.sect=(wdraw.childs+1)->org.sects+0;

			nnode=(wdraw.childs+1)->org.nnode;
			find=0;
			for(i=0;i<nnode;i++)
			{
			  pnode=(wdraw.childs+1)->org.nodes+i;
			  if(pnode->d[GX]-0.010<=node->d[GX]
				 && node->d[GX]<=pnode->d[GX]+0.010 &&
				 pnode->d[GY]-0.010<=node->d[GY]
				 && node->d[GY]<=pnode->d[GY]+0.010 &&
				 pnode->d[GZ]-0.010<=0.0
				 && 0.0<=pnode->d[GZ]+0.010)
			  {
				find=1;
				break;
			  }
			}
			if(!find)
			{
			  oc=(struct oconf *)malloc(6*sizeof(struct oconf));
			  (oc+0)->iconf=1;
			  (oc+1)->iconf=1;
			  (oc+2)->iconf=1;
			  (oc+3)->iconf=0;
			  (oc+4)->iconf=0;
			  (oc+5)->iconf=0;

			  (oc+0)->value=0.0;
			  (oc+1)->value=0.0;
			  (oc+2)->value=0.0;
			  (oc+3)->value=0.0;
			  (oc+4)->value=0.0;
			  (oc+5)->value=0.0;

			  cnode.code=((wdraw.childs+1)->org.nodes+nnode-1)
						 ->code+1;
			  cnode.loff=nnode;

			  cnode.d[0]=node->d[0];
			  cnode.d[1]=node->d[1];
			  cnode.d[2]=0.0;

			  pnode=addnode(cnode,oc,&((wdraw.childs+1)->org));

			  /*FIND ON ELEM : UNDER CONSTRUCTION.*/


			}

			*(gelem.nods+1)=(wdraw.childs+1)->org.nodes+loff;
			*(gelem.nods+0)=pnode;
			addelement(&gelem,&((wdraw.childs+1)->org));

			/*
			gelem.nnod=0;
			gelem.nban=0;
			gelem.nods=NULL;
			gelem.bonds=NULL;
			gelem.bans=NULL;
			*/

			globalmessageflag=0;

			/*INITIAL*/
			arc =arci;
			arcx=arci;
			arcy=arci;
			free((wdraw.childs+1)->org.loads);
			for(i=0;i<(wdraw.childs+1)->org.nelem;i++)
			{
			  for(j=0;j<2;j++)
			  {
				for(k=0;k<6;k++)
				{
				  ((wdraw.childs+1)->org.elems+i)->initial[j][k]=0.0;
				}
			  }
			} /*UNAVAILABLE INITIAL CMQ.*/

			pnelem=(wdraw.childs+1)->org.nelem;
			pelems=(wdraw.childs+1)->org.elems;

			/*EXTRACT ARCLM*/
			extractarclmfromorgan(&((wdraw.childs+1)->org),
								  &arc,&arcx,&arcy);

			(wdraw.childs+1)->org.nelem=pnelem;
			(wdraw.childs+1)->org.elems=pelems;

			(wmenu.childs+2)->vparam.vflag.mv.ftype=F_ARCLM;

			/*SAVE AS ARCLM*/
			saveasarclm((wdraw.childs+1)->inpfilez,&arc);
			saveasarclm((wdraw.childs+1)->inpfilex,&arcx);
			saveasarclm((wdraw.childs+1)->inpfiley,&arcy);

			/*ANALYSIS*/
			getviewparam((wmenu.childs+2)->hwnd,
						 &((wdraw.childs+1)->vparam));
			(wdraw.childs+1)->vparam.vflag.ev.deformation=0;

			arclm001(&arc ,ID_INPUTFILEZ,ID_OUTPUTFILEZ);
			arclm001(&arcx,ID_INPUTFILEX,ID_OUTPUTFILEX);
			arclm001(&arcy,ID_INPUTFILEY,ID_OUTPUTFILEY);

			/*SRCAN*/
			n=strcspn((wdraw.childs+1)->sctfile,".");
			strncpy(str,(wdraw.childs+1)->sctfile,n);
			str[n]='\0';

			i=srcan001(str);

			if(i==0) MessageBox(NULL,"Failed.","SRCAN001",MB_OK);

			strcpy((wdraw.childs+1)->otpfile,str);
			strcat((wdraw.childs+1)->otpfile,".rat");
			SetDlgItemText((wmenu.childs+2)->hwnd,ID_OUTPUTFILE,
						   (wdraw.childs+1)->otpfile);

			fout=fopen((wdraw.childs+1)->otpfile,"r");
			if(fout==NULL) break;
			readsrcanrate(fout,&arc);

			/*SAFETY FLAG*/
			(wdraw.childs+1)->vparam.vflag.ev.srcancolor=1;

			clearwindow(*(wdraw.childs+1));
			drawarclmframe((wdraw.childs+1)->hdcC,
						   (wdraw.childs+1)->vparam,arc,0,ONSCREEN);

			SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);

			(wmenu.childs+2)->vparam.vflag.mv.ftype=F_ORGAN;
			SendMessage((wmenu.childs+2)->hwnd,WM_INITDIALOG,0,0);

			globalmessageflag=1;
		  }
		  break;
		case IDM_POPDELETECOLUMNWITHSAFETY:
		  if(wdraw.nchilds>=2 &&
			 wmenu.nchilds>=3 &&
			 (wmenu.childs+2)->vparam.vflag.mv.ftype==F_ORGAN)
		  {
			/*hcursor = LoadCursor(hInstGlobal,"CANBOXW");*/
			/*SetClassLong((wdraw.childs+1)->hwnd,
						 GCL_HCURSOR,(LONG)hcursor);*/

			/*globalstatus=DELETECOLUMNWITHSAFETY;*/
			/*createitem=C_ELEMENT;*/

			node=selectorgannode((wdraw.childs+1)->vparam,
								 &((wdraw.childs+1)->org),point);
			if(node==NULL) break;
			loff=node->loff;

			nelem=(wdraw.childs+1)->org.nelem;
			for(i=0;i<nelem;i++)
			{
			  if((*(((wdraw.childs+1)->org.elems+i)->nods+0) == node ||
				  *(((wdraw.childs+1)->org.elems+i)->nods+1) == node) &&
				 ((wdraw.childs+1)->org.elems+i)->sect
				 == (wdraw.childs+1)->org.sects+0)
			  {
				if((*(((wdraw.childs+1)->org.elems+i)->nods+0))->d[GZ] <=
				   (*(((wdraw.childs+1)->org.elems+i)->nods+1))->d[GZ])
				{
				  loff=(*(((wdraw.childs+1)->org.elems+i)->nods+0))->loff;
				}
				else
				{
				  loff=(*(((wdraw.childs+1)->org.elems+i)->nods+1))->loff;
				}
				for(j=0;j<6;j++)
				{
				  ((wdraw.childs+1)->org.confs+6*loff+j)->iconf=0;
				}

				deleteelement(i,&((wdraw.childs+1)->org),0);
				break;
			  }
			}

			globalmessageflag=0;

			/*INITIAL*/
			arc =arci;
			arcx=arci;
			arcy=arci;
			free((wdraw.childs+1)->org.loads);
			for(i=0;i<(wdraw.childs+1)->org.nelem;i++)
			{
			  for(j=0;j<2;j++)
			  {
				for(k=0;k<6;k++)
				{
				  ((wdraw.childs+1)->org.elems+i)->initial[j][k]=0.0;
				}
			  }
			} /*UNAVAILABLE INITIAL CMQ.*/

			pnelem=(wdraw.childs+1)->org.nelem;
			pelems=(wdraw.childs+1)->org.elems;

			/*EXTRACT ARCLM*/
			extractarclmfromorgan(&((wdraw.childs+1)->org),
								  &arc,&arcx,&arcy);

			(wdraw.childs+1)->org.nelem=pnelem;
			(wdraw.childs+1)->org.elems=pelems;

			(wmenu.childs+2)->vparam.vflag.mv.ftype=F_ARCLM;

			/*SAVE AS ARCLM*/
			saveasarclm((wdraw.childs+1)->inpfilez,&arc);
			saveasarclm((wdraw.childs+1)->inpfilex,&arcx);
			saveasarclm((wdraw.childs+1)->inpfiley,&arcy);

			/*ANALYSIS*/
			getviewparam((wmenu.childs+2)->hwnd,
						 &((wdraw.childs+1)->vparam));
			(wdraw.childs+1)->vparam.vflag.ev.deformation=0;

			arclm001(&arc ,ID_INPUTFILEZ,ID_OUTPUTFILEZ);
			arclm001(&arcx,ID_INPUTFILEX,ID_OUTPUTFILEX);
			arclm001(&arcy,ID_INPUTFILEY,ID_OUTPUTFILEY);

			/*SRCAN*/
			n=strcspn((wdraw.childs+1)->sctfile,".");
			strncpy(str,(wdraw.childs+1)->sctfile,n);
			str[n]='\0';

			i=srcan001(str);

			if(i==0) MessageBox(NULL,"Failed.","SRCAN001",MB_OK);

			strcpy((wdraw.childs+1)->otpfile,str);
			strcat((wdraw.childs+1)->otpfile,".rat");
			SetDlgItemText((wmenu.childs+2)->hwnd,ID_OUTPUTFILE,
						   (wdraw.childs+1)->otpfile);

			fout=fopen((wdraw.childs+1)->otpfile,"r");
			if(fout==NULL) break;
			readsrcanrate(fout,&arc);

			/*SAFETY FLAG*/
			(wdraw.childs+1)->vparam.vflag.ev.srcancolor=1;

			clearwindow(*(wdraw.childs+1));
			drawarclmframe((wdraw.childs+1)->hdcC,
						   (wdraw.childs+1)->vparam,arc,0,ONSCREEN);

			SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);

			(wmenu.childs+2)->vparam.vflag.mv.ftype=F_ORGAN;
			SendMessage((wmenu.childs+2)->hwnd,WM_INITDIALOG,0,0);

			globalmessageflag=1;
		  }
		  break;

		case IDM_POPMOVENODEWITHSAFETY:
		  if(wdraw.nchilds>=2 &&
			 wmenu.nchilds>=3 &&
			 (wmenu.childs+2)->vparam.vflag.mv.ftype==F_ORGAN)
		  {
			hcursor = LoadCursor(hInstGlobal,"CANBOXW");
			SetClassLong((wdraw.childs+1)->hwnd,
						 GCL_HCURSOR,(LONG)hcursor);

			globalstatus=MOVENODEWITHSAFETY;
		  }
		  break;

		case IDM_POPMOVENODEDIALOG:
		  if(wdraw.nchilds>=2 &&
			 wmenu.nchilds>=3 &&
			 (wmenu.childs+2)->vparam.vflag.mv.ftype==F_ORGAN)
		  {
			point.x=initx;
			point.y=inity;

            node=selectorgannode((wdraw.childs+1)->vparam,
                                 &((wdraw.childs+1)->org),point);
            if(node==NULL)
            {
              MessageBox(NULL,"Nothing.","MoveNode",MB_OK);
              break;
            }
            pnode=node;

            /*DRAW SELECTED NODE.*/
            nodeontoscreen(*node,&Nx,&Ny,
                           (wdraw.childs+1)->vparam); /*PROJECTION*/
            /*hdc = GetDC((wdraw.childs+1)->hwnd);*/
            hdc = (wdraw.childs+1)->hdcC;
            setfontformat(hdc,15,6,"MS Mincho",255,255,255);
            sprintf(str,"%ld",node->code);
            TextOut(hdc,Nx+2,Ny,str,strlen(str));
            ReleaseDC((wdraw.childs+1)->hwnd,hdc);
            SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);

            hpopdlg=CreateDialog(hInstGlobal,
                                 "HOGDLGINCREMENT",
                                 NULL,
                                 (DLGPROC)DialogProcIncrement);
            getwindowsize(hpopdlg,&mw,&mh);
            ClientToScreen((wdraw.childs+1)->hwnd,&point);
            MoveWindow(hpopdlg,point.x,point.y,mw,mh,TRUE);
            ShowWindow(hpopdlg,SW_SHOW);

            gincrement.dc[0]=0.0;
            gincrement.dc[1]=0.0;
            gincrement.dc[2]=0.0;

            globalstatus=MOVENODEBYDIALOG;
          }
          break;
        case IDM_POPMOVENODENODETONODE:
        case IDM_POPCOPYNODENODETONODE:
          if(wdraw.nchilds>=2 &&
             wmenu.nchilds>=3 &&
             (wmenu.childs+2)->vparam.vflag.mv.ftype==F_ORGAN)
          {
            point.x=initx;
            point.y=inity;

            node=selectorgannode((wdraw.childs+1)->vparam,
                                 &((wdraw.childs+1)->org),point);
            if(node==NULL)
            {
              MessageBox(NULL,"Nothing.","MoveNode",MB_OK);
              break;
            }
            pnode=node;

            /*DRAW SELECTED NODE.*/
            nodeontoscreen(*node,&Nx,&Ny,
                           (wdraw.childs+1)->vparam); /*PROJECTION*/
            /*hdc = GetDC((wdraw.childs+1)->hwnd);*/
            hdc = (wdraw.childs+1)->hdcC;
            setfontformat(hdc,15,6,"MS Mincho",255,255,255);
            sprintf(str,"%ld",node->code);
			TextOut(hdc,Nx+2,Ny,str,strlen(str));
            ReleaseDC((wdraw.childs+1)->hwnd,hdc);
            SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);

            icount=0;
            gincrement.dc[0]=0.0;
            gincrement.dc[1]=0.0;
            gincrement.dc[2]=0.0;

            hcursor = LoadCursor(hInstGlobal,"CANBOXW");
            SetClassLong((wdraw.childs+1)->hwnd,
                         GCL_HCURSOR,(LONG)hcursor);

            if(LOWORD(wParam)==IDM_POPMOVENODENODETONODE)
            {
              globalstatus=MOVENODENODETONODE;
            }
            if(LOWORD(wParam)==IDM_POPCOPYNODENODETONODE)
            {
              globalstatus=COPYNODENODETONODE;
            }
          }
          break;

        case IDM_POPDELETENODE:
          if(wdraw.nchilds>=2 &&
             wmenu.nchilds>=3 &&
             (wmenu.childs+2)->vparam.vflag.mv.ftype==F_ORGAN)
          {
            point.x=initx;
            point.y=inity;

            node=selectorgannode((wdraw.childs+1)->vparam,
                                 &((wdraw.childs+1)->org),point);
            if(node==NULL)
            {
              MessageBox(NULL,"Nothing.","MoveNode",MB_OK);
              break;
            }

			/*DRAW SELECTED NODE.*/
			nodeontoscreen(*node,&Nx,&Ny,
						   (wdraw.childs+1)->vparam); /*PROJECTION*/
			/*hdc = GetDC((wdraw.childs+1)->hwnd);*/
			hdc = (wdraw.childs+1)->hdcC;
			setfontformat(hdc,15,6,"MS Mincho",255,255,255);
			sprintf(str,"%ld",node->code);
			TextOut(hdc,Nx+2,Ny,str,strlen(str));
			ReleaseDC((wdraw.childs+1)->hwnd,hdc);
			SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);

			if(MessageBox(NULL,"Delete or Cancel","DeleteNode",
						  MB_OKCANCEL)==IDOK)
			{
			  deletenodewithelem(node->code,
								 &((wdraw.childs+1)->org));
			}

			clearwindow(*(wdraw.childs+1));
			draworganization((wdraw.childs+1)->hdcC,
							 (wdraw.childs+1)->vparam,
							 (wdraw.childs+1)->org,ONSCREEN);
			SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);
		  }
		  break;

		case IDM_POPDELETEELEM:
		  if(wdraw.nchilds>=2 &&
			 wmenu.nchilds>=3 &&
			 (wmenu.childs+2)->vparam.vflag.mv.ftype==F_ORGAN)
		  {
			point.x=initx;
			point.y=inity;

			foundelem=selectorganelement((wdraw.childs+1)->vparam,
										 &((wdraw.childs+1)->org),
										 point,
										 &selectedline,
										 &selectedban);
			if(foundelem==NULL)
			{
			  MessageBox(hwnd,"Nothing.","Element",MB_OK);
			  break;
			}

			/*DRAW NODE OF SELECTED ELEM.*/
			/*hdc = GetDC((wdraw.childs+1)->hwnd);*/
			hdc = (wdraw.childs+1)->hdcC;
			setfontformat(hdc,15,6,"MS Mincho",255,255,255);
			for(i=0;i<(foundelem->nnod);i++)
			{
			  nodeontoscreen(**(foundelem->nods+i),&Nx,&Ny,
							 (wdraw.childs+1)->vparam);
			  sprintf(str,"%ld",(*(foundelem->nods+i))->code);
			  TextOut(hdc,Nx+2,Ny,str,strlen(str));
			}
			ReleaseDC((wdraw.childs+1)->hwnd,hdc);
			SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);

			if(MessageBox(NULL,"Delete or Cancel","DeleteElem",
						  MB_OKCANCEL)==IDOK)
			{
			  deleteelement(foundelem->loff,
							&((wdraw.childs+1)->org),0);
			}

			clearwindow(*(wdraw.childs+1));
			draworganization((wdraw.childs+1)->hdcC,
							 (wdraw.childs+1)->vparam,
							 (wdraw.childs+1)->org,ONSCREEN);
			SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);
		  }
		  break;

		case IDM_POPMOVEMOUSE:
		  if(wdraw.nchilds>=2 &&
			 wmenu.nchilds>=3 &&
			 (wmenu.childs+2)->vparam.vflag.mv.ftype==F_ORGAN)
		  {
			globalstatus=MOVEELEMENTBYMOUSE;
		  }
		  break;
		case IDM_POPCOPYMOUSE:
		  if(wdraw.nchilds>=2 &&
			 wmenu.nchilds>=3 &&
			 (wmenu.childs+2)->vparam.vflag.mv.ftype==F_ORGAN)
		  {
			globalstatus=COPYELEMENTBYMOUSE;
		  }
		  break;
		case IDM_POPMOVEDIALOG:
		case IDM_POPCOPYDIALOG:
		  if(wdraw.nchilds>=2 &&
			 wmenu.nchilds>=3 &&
			 (wmenu.childs+2)->vparam.vflag.mv.ftype==F_ORGAN)
		  {
            point.x=initx;
            point.y=inity;

            foundelem=selectorganelement((wdraw.childs+1)->vparam,
                                         &((wdraw.childs+1)->org),
                                         point,
                                         &selectedline,
                                         &selectedban);
            if(foundelem==NULL)
            {
              MessageBox(hwnd,"Nothing.","Element",MB_OK);
              break;
            }

            /*DRAW NODE OF SELECTED ELEM.*/
            /*hdc = GetDC((wdraw.childs+1)->hwnd);*/
            hdc = (wdraw.childs+1)->hdcC;
            setfontformat(hdc,15,6,"MS Mincho",255,255,255);
            for(i=0;i<(foundelem->nnod);i++)
            {
              nodeontoscreen(**(foundelem->nods+i),&Nx,&Ny,
                             (wdraw.childs+1)->vparam);
              sprintf(str,"%ld",(*(foundelem->nods+i))->code);
              TextOut(hdc,Nx+2,Ny,str,strlen(str));
            }
            ReleaseDC((wdraw.childs+1)->hwnd,hdc);
            SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);

            /*if(LOWORD(wParam)==IDM_POPMOVEDIALOG)
            {
              if(MessageBox(NULL,"Move or Cancel","MoveElem",
                            MB_OKCANCEL)==IDCANCEL) break;
			}
            if(LOWORD(wParam)==IDM_POPCOPYDIALOG)
            {
              if(MessageBox(NULL,"Copy or Cancel","CopyElem",
                            MB_OKCANCEL)==IDCANCEL) break;
            }*/

            pelem=foundelem;

            hpopdlg=CreateDialog(hInstGlobal,
                                 "HOGDLGINCREMENT",
                                 NULL,
                                 (DLGPROC)DialogProcIncrement);
            getwindowsize(hpopdlg,&mw,&mh);
            ClientToScreen((wdraw.childs+1)->hwnd,&point);
            MoveWindow(hpopdlg,point.x,point.y,mw,mh,TRUE);
            ShowWindow(hpopdlg,SW_SHOW);

            gincrement.dc[0]=0.0;
            gincrement.dc[1]=0.0;
            gincrement.dc[2]=0.0;

            if(LOWORD(wParam)==IDM_POPMOVEDIALOG)
            {
              globalstatus=MOVEELEMENTBYDIALOG;
            }
            if(LOWORD(wParam)==IDM_POPCOPYDIALOG)
            {
              globalstatus=COPYELEMENTBYDIALOG;
            }
          }
          break;
        case IDM_POPINCREMENTRETURN:
          GetDlgItemText(hpopdlg,IDT_POPDX,str,256);
          dx=strtod(str,NULL);
          GetDlgItemText(hpopdlg,IDT_POPDY,str,256);
          dy=strtod(str,NULL);
          GetDlgItemText(hpopdlg,IDT_POPDZ,str,256);
          dz=strtod(str,NULL);

          if(globalstatus==MOVENODEBYDIALOG)
          {
			pnode->d[0]+=dx;
            pnode->d[1]+=dy;
            pnode->d[2]+=dz;
          }
          else if(globalstatus==MOVEELEMENTBYDIALOG)
          {
            for(i=0;i<pelem->nnod;i++)
            {
              (*(pelem->nods+i))->d[0]+=dx;
              (*(pelem->nods+i))->d[1]+=dy;
              (*(pelem->nods+i))->d[2]+=dz;
            }
          }
          else if(globalstatus==COPYELEMENTBYDIALOG)
          {
            addelementwithnode(pelem,
                               &((wdraw.childs+1)->org),
                               dx,dy,dz);
                               /*UNDER CONSTRUCTION*/
          }
          else if(globalstatus==COPYMULTIELEMBYDIALOG)
          {
            copymultielem(&((wdraw.childs+1)->org),dx,dy,dz);
            nmultielem=0;
            free(multielem);
          }
          else if(globalstatus==MOVEMULTIELEMBYDIALOG)
          {
            movemultielem(&((wdraw.childs+1)->org),dx,dy,dz);
            nmultielem=0;
            free(multielem);
          }

          clearwindow(*(wdraw.childs+1));
          draworganization((wdraw.childs+1)->hdcC,
                           (wdraw.childs+1)->vparam,
                           (wdraw.childs+1)->org,ONSCREEN);
          SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);

          globalstatus=NEUTRAL;
          break;
        case IDM_POPMOVENODETONODE:
		case IDM_POPCOPYNODETONODE:
          if(wdraw.nchilds>=2 &&
             wmenu.nchilds>=3 &&
             (wmenu.childs+2)->vparam.vflag.mv.ftype==F_ORGAN)
          {
            point.x=initx;
            point.y=inity;

            foundelem=selectorganelement((wdraw.childs+1)->vparam,
                                         &((wdraw.childs+1)->org),
                                         point,
                                         &selectedline,
                                         &selectedban);
            if(foundelem==NULL)
            {
              MessageBox(hwnd,"Nothing.","Element",MB_OK);
              break;
            }

            /*DRAW NODE OF SELECTED ELEM.*/
            /*hdc = GetDC((wdraw.childs+1)->hwnd);*/
            hdc = (wdraw.childs+1)->hdcC;
            setfontformat(hdc,15,6,"MS Mincho",255,255,255);
            for(i=0;i<(foundelem->nnod);i++)
            {
              nodeontoscreen(**(foundelem->nods+i),&Nx,&Ny,
                             (wdraw.childs+1)->vparam);
              sprintf(str,"%ld",(*(foundelem->nods+i))->code);
              TextOut(hdc,Nx+2,Ny,str,strlen(str));
            }
            ReleaseDC((wdraw.childs+1)->hwnd,hdc);
            SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);

            /*if(LOWORD(wParam)==IDM_POPMOVENODETONODE)
            {
              if(MessageBox(NULL,"Move or Cancel","MoveElem",
                            MB_OKCANCEL)==IDCANCEL) break;
            }
            if(LOWORD(wParam)==IDM_POPCOPYNODETONODE)
            {
              if(MessageBox(NULL,"Copy or Cancel","CopyElem",
                            MB_OKCANCEL)==IDCANCEL) break;
			}*/

            pelem=foundelem;

            icount=0;
            gincrement.dc[0]=0.0;
            gincrement.dc[1]=0.0;
            gincrement.dc[2]=0.0;

            hcursor = LoadCursor(hInstGlobal,"CANBOXW");
            SetClassLong((wdraw.childs+1)->hwnd,
                         GCL_HCURSOR,(LONG)hcursor);

            if(LOWORD(wParam)==IDM_POPMOVENODETONODE)
            {
              globalstatus=MOVEELEMENTNODETONODE;
            }
            if(LOWORD(wParam)==IDM_POPCOPYNODETONODE)
            {
              globalstatus=COPYELEMENTNODETONODE;
            }
          }
          break;
        case IDM_POPCOPYCONF:
          if(wdraw.nchilds>=2 &&
             wmenu.nchilds>=3 &&
             (wmenu.childs+2)->vparam.vflag.mv.ftype==F_ORGAN)
          {
            point.x=initx;
            point.y=inity;

            node=selectorgannode((wdraw.childs+1)->vparam,
                                 &((wdraw.childs+1)->org),point);
            if(node==NULL)
            {
              MessageBox(NULL,"Nothing.","MoveNode",MB_OK);
              break;
            }
            pnode=node;

            /*DRAW SELECTED NODE.*/
            nodeontoscreen(*node,&Nx,&Ny,
						   (wdraw.childs+1)->vparam); /*PROJECTION*/
            /*hdc = GetDC((wdraw.childs+1)->hwnd);*/
            hdc = (wdraw.childs+1)->hdcC;
            setfontformat(hdc,15,6,"MS Mincho",255,255,255);
            sprintf(str,"%ld",node->code);
            TextOut(hdc,Nx+2,Ny,str,strlen(str));
            ReleaseDC((wdraw.childs+1)->hwnd,hdc);
            SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);

            hcursor = LoadCursor(hInstGlobal,"CANBOXW");
            SetClassLong((wdraw.childs+1)->hwnd,
                         GCL_HCURSOR,(LONG)hcursor);

            globalstatus=COPYCONF;
          }
          break;
        case IDM_POPMOVEELEMEND:
        case IDM_POPCOPYELEMEND:
          if(wdraw.nchilds>=2)
          {
            hcursor = LoadCursor(hInstGlobal,"CANCURSORW");
            SetClassLong((wdraw.childs+1)->hwnd,
                         GCL_HCURSOR,(LONG)hcursor);
            globalstatus=NEUTRAL;
          }
          break;

        case IDM_POPCHANGEETYPE:
          if(wdraw.nchilds>=2 &&
             wmenu.nchilds>=3 &&
             (wmenu.childs+2)->vparam.vflag.mv.ftype==F_ORGAN)
          {
            point.x=initx;
            point.y=inity;

            foundelem=selectorganelement((wdraw.childs+1)->vparam,
                                         &((wdraw.childs+1)->org),
                                         point,
                                         &selectedline,
                                         &selectedban);
            if(foundelem==NULL)
            {
			  MessageBox(hwnd,"Nothing.","Element",MB_OK);
              break;
            }
            pelem=foundelem;
            gsect=foundelem->sect;

            /*DRAW SELECTED ELEMENT*/
            /*hdc = GetDC((wdraw.childs+1)->hwnd);*/
            hdc = (wdraw.childs+1)->hdcC;
            setfontformat(hdc,15,6,"MS Mincho",255,255,255);
            for(i=0;i<(foundelem->nnod);i++)
            {
              nodeontoscreen(**(foundelem->nods+i),&Nx,&Ny,
                             (wdraw.childs+1)->vparam);
              sprintf(str,"%ld",(*(foundelem->nods+i))->code);
              TextOut(hdc,Nx+2,Ny,str,strlen(str));
            }
            ReleaseDC((wdraw.childs+1)->hwnd,hdc);
            SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);

            hpopdlg=CreateDialog(hInstGlobal,
                                 "HOGDLGETYPE",
                                 NULL,
                                 (DLGPROC)DialogProcElemType);
            getwindowsize(hpopdlg,&mw,&mh);
            ClientToScreen((wdraw.childs+1)->hwnd,&point);
            MoveWindow(hpopdlg,point.x,point.y,mw,mh,TRUE);
            ShowWindow(hpopdlg,SW_SHOW);

            globalstatus=CHANGEELEMTYPE;
          }
          break;
        case IDM_POPETYPERETURN:
          clearwindow(*(wdraw.childs+1));
          draworganization((wdraw.childs+1)->hdcC,
                           (wdraw.childs+1)->vparam,
                           (wdraw.childs+1)->org,ONSCREEN);
          SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);

          globalstatus=NEUTRAL;
          break;

		case IDM_POPCHANGESECT:                       /*change sect*/
          if(wdraw.nchilds>=2 &&
             wmenu.nchilds>=3 &&
             (wmenu.childs+2)->vparam.vflag.mv.ftype==F_ORGAN)
          {
            point.x=initx;
            point.y=inity;

            foundelem=selectorganelement((wdraw.childs+1)->vparam,
                                         &((wdraw.childs+1)->org),
                                         point,
                                         &selectedline,
                                         &selectedban);
            if(foundelem==NULL)
            {
              MessageBox(hwnd,"Nothing.","Element",MB_OK);
              break;
            }
            pelem=foundelem;
            gsect=foundelem->sect;

            /*DRAW SELECTED ELEMENT*/
            /*hdc = GetDC((wdraw.childs+1)->hwnd);*/
            hdc = (wdraw.childs+1)->hdcC;
            setfontformat(hdc,15,6,"MS Mincho",255,255,255);
            for(i=0;i<(foundelem->nnod);i++)
            {
              nodeontoscreen(**(foundelem->nods+i),&Nx,&Ny,
                             (wdraw.childs+1)->vparam);
              sprintf(str,"%ld",(*(foundelem->nods+i))->code);
              TextOut(hdc,Nx+2,Ny,str,strlen(str));
            }
            ReleaseDC((wdraw.childs+1)->hwnd,hdc);
            SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);

            hpopdlg=CreateDialog(hInstGlobal,
                                 "HOGDLGSECTREGIST",
                                 NULL,
                                 (DLGPROC)DialogProcSectRegist);
            getwindowsize(hpopdlg,&mw,&mh);
            ClientToScreen((wdraw.childs+1)->hwnd,&point);
            MoveWindow(hpopdlg,point.x,point.y,mw,mh,TRUE);
			ShowWindow(hpopdlg,SW_SHOW);

            globalstatus=CHANGESECTION;
          }
          break;

        case IDM_POPREGISTRETURN:
          /*CHANGE SECTION POINTER*/
          GetDlgItemText(hpopdlg,IDSR_POPCODE,str,20);
          code=strtol(str,NULL,10);

          if(globalstatus==CHANGESECTMULTIELEM && code!=0
            && (wmenu.childs+2)->vparam.vflag.mv.ftype==F_ORGAN)
          {
            for(i=0;i<((wdraw.childs+1)->org.nsect);i++)
            {
              if(code==((wdraw.childs+1)->org.sects+i)->code)
              {
                for(j=0;j<nmultielem;j++)
                {
                  (*(multielem+j))->sect=(wdraw.childs+1)->org.sects+i;
                }
                nmultielem=0;
                free(multielem);
                break;
              }
            }
          }
          if(globalstatus==CHANGESECTMULTIELEM && code!=0
            && (wmenu.childs+2)->vparam.vflag.mv.ftype==F_ARCLM)
          {
            for(i=0;i<arc.nsect;i++)
            {
              if(code==(arc.sects+i)->code)
              {
                for(j=0;j<nmultiwire;j++)
                {
                  (*(multiwire+j))->sect=arc.sects+i;
                }
                break;
              }
            }
		  }
          else if(globalstatus==CHANGESECTION && code!=0)
          {
            for(i=0;i<((wdraw.childs+1)->org.nsect);i++)
            {
              if(code==((wdraw.childs+1)->org.sects+i)->code)
              {
                pelem->sect=(wdraw.childs+1)->org.sects+i;
                break;
              }
            }
          }

          clearwindow(*(wdraw.childs+1));

          if((wmenu.childs+2)->vparam.vflag.mv.ftype==F_ORGAN)
          {
            drawarclmframe((wdraw.childs+1)->hdcC,
                           (wdraw.childs+1)->vparam,arc,0,ONSCREEN);

          }
          if((wmenu.childs+2)->vparam.vflag.mv.ftype==F_ARCLM)
          {
            draworganization((wdraw.childs+1)->hdcC,
                             (wdraw.childs+1)->vparam,
                             (wdraw.childs+1)->org,ONSCREEN);
          }
          SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);

          globalstatus=NEUTRAL;
          break;

        case IDM_POPCHANGECONF:
          if(wdraw.nchilds>=2 &&
             wmenu.nchilds>=3 &&
             (wmenu.childs+2)->vparam.vflag.mv.ftype==F_ORGAN)
          {
            point.x=initx;
            point.y=inity;

            node=selectorgannode((wdraw.childs+1)->vparam,
                                 &((wdraw.childs+1)->org),point);
			if(node==NULL)
            {
              MessageBox(NULL,"Nothing.","MoveNode",MB_OK);
              break;
            }
            pnode=node;

            /*DRAW SELECTED NODE.*/
            nodeontoscreen(*node,&Nx,&Ny,
                           (wdraw.childs+1)->vparam); /*PROJECTION*/
            /*hdc = GetDC((wdraw.childs+1)->hwnd);*/
            hdc = (wdraw.childs+1)->hdcC;
            setfontformat(hdc,15,6,"MS Mincho",255,255,255);
            sprintf(str,"%ld",node->code);
            TextOut(hdc,Nx+2,Ny,str,strlen(str));
            ReleaseDC((wdraw.childs+1)->hwnd,hdc);
            SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);

            hpopdlg=CreateDialog(hInstGlobal,
                                 "HOGDLGCHANGECONF",
                                 NULL,
                                 (DLGPROC)DialogProcConf);
            getwindowsize(hpopdlg,&mw,&mh);
            ClientToScreen((wdraw.childs+1)->hwnd,&point);
            MoveWindow(hpopdlg,point.x,point.y,mw,mh,TRUE);
            ShowWindow(hpopdlg,SW_SHOW);

            globalstatus=CHANGECONFBYDIALOG;
          }
          break;
        case IDM_POPCONFRETURN:
          /*REDRAW*/
          if(globalstatus==CHANGECONFMULTINODE)
		  {
		  changeconfmultinode(&(wdraw.childs+1)->org);
          }
		  clearwindow(*(wdraw.childs+1));
          draworganization((wdraw.childs+1)->hdcC,
                           (wdraw.childs+1)->vparam,
                           (wdraw.childs+1)->org,ONSCREEN);
		  SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);

		  globalstatus=NEUTRAL;
          break;

        case IDM_POPCHANGEBOND:
          if(wdraw.nchilds>=2 &&
             wmenu.nchilds>=3 &&
             (wmenu.childs+2)->vparam.vflag.mv.ftype==F_ORGAN)
          {
            point.x=initx;
            point.y=inity;

            foundelem=selectorganelement((wdraw.childs+1)->vparam,
                                         &((wdraw.childs+1)->org),
                                         point,
                                         &selectedline,
                                         &selectedban);
            if(foundelem==NULL)
            {
              MessageBox(hwnd,"Nothing.","Element",MB_OK);
              break;
            }
            pelem=foundelem;

            /*DRAW SELECTED ELEMENT*/
            /*hdc = GetDC((wdraw.childs+1)->hwnd);*/
            hdc = (wdraw.childs+1)->hdcC;
			setfontformat(hdc,15,6,"MS Mincho",255,255,255);
            for(i=0;i<(foundelem->nnod);i++)
            {
              nodeontoscreen(**(foundelem->nods+i),&Nx,&Ny,
                             (wdraw.childs+1)->vparam);
              sprintf(str,"%ld",(*(foundelem->nods+i))->code);
              TextOut(hdc,Nx+2,Ny,str,strlen(str));
            }
            ReleaseDC((wdraw.childs+1)->hwnd,hdc);
            SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);

            hcursor = LoadCursor(hInstGlobal,"CANBOXW");
            SetClassLong((wdraw.childs+1)->hwnd,
                         GCL_HCURSOR,(LONG)hcursor);

            globalstatus = CHANGEBOND;
		  }
          break;
        case IDM_POPBONDRETURN:
          /*REDRAW FRAME*/
          clearwindow(*(wdraw.childs+1));
          draworganization((wdraw.childs+1)->hdcC,
                           (wdraw.childs+1)->vparam,
                           (wdraw.childs+1)->org,ONSCREEN);

          /*REDRAW SELECTED ELEM*/
          hdc = (wdraw.childs+1)->hdcC;
          setfontformat(hdc,15,6,"MS Mincho",255,255,255);
          for(i=0;i<(pelem->nnod);i++)
          {
            nodeontoscreen(**(pelem->nods+i),&Nx,&Ny,
                           (wdraw.childs+1)->vparam);
            sprintf(str,"%ld",(*(pelem->nods+i))->code);
            TextOut(hdc,Nx+2,Ny,str,strlen(str));
          }
          ReleaseDC((wdraw.childs+1)->hwnd,hdc);

          SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);

          /*globalstatus=NEUTRAL;*/
          break;
        case IDM_POPCHANGEEND:
          /*REDRAW*/
          clearwindow(*(wdraw.childs+1));
          draworganization((wdraw.childs+1)->hdcC,
                           (wdraw.childs+1)->vparam,
                           (wdraw.childs+1)->org,ONSCREEN);
          SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);

          hcursor = LoadCursor(hInstGlobal,"CANCURSORW");
          SetClassLong((wdraw.childs+1)->hwnd,
                       GCL_HCURSOR,(LONG)hcursor);

          globalstatus=NEUTRAL;
          break;

        case IDM_POPSELECTMULTIELEM:
          globalstatus=SELECTMULTIELEM;
		  break;
        case IDM_POPCOPYMULTIELEMBYDIALOG:
        case IDM_POPMOVEMULTIELEMBYDIALOG:
          if(wdraw.nchilds>=2 &&
             wmenu.nchilds>=3 &&
             (wmenu.childs+2)->vparam.vflag.mv.ftype==F_ORGAN)
          {
            if(globalstatus!=SELECTMULTIELEM ||
               nmultielem==0) break;

            point.x=initx;
            point.y=inity;

            hpopdlg=CreateDialog(hInstGlobal,
                                 "HOGDLGINCREMENT",
                                 NULL,
                                 (DLGPROC)DialogProcIncrement);
            getwindowsize(hpopdlg,&mw,&mh);
            ClientToScreen((wdraw.childs+1)->hwnd,&point);
            MoveWindow(hpopdlg,point.x,point.y,mw,mh,TRUE);
            ShowWindow(hpopdlg,SW_SHOW);

            if(LOWORD(wParam)==IDM_POPCOPYMULTIELEMBYDIALOG)
            {
              globalstatus=COPYMULTIELEMBYDIALOG;
            }
            if(LOWORD(wParam)==IDM_POPMOVEMULTIELEMBYDIALOG)
            {
              globalstatus=MOVEMULTIELEMBYDIALOG;
            }
          }
          break;
        case IDM_POPDELETEMULTIELEM:
          if(wdraw.nchilds>=2 &&
             wmenu.nchilds>=3 &&
             (wmenu.childs+2)->vparam.vflag.mv.ftype==F_ORGAN)
          {
            if(globalstatus!=SELECTMULTIELEM ||
               nmultielem==0) break;

            deletemultielem(&((wdraw.childs+1)->org));
            nmultielem=0;
			free(multielem);

            globalstatus=NEUTRAL;

            /*REDRAW*/
            clearwindow(*(wdraw.childs+1));
            draworganization((wdraw.childs+1)->hdcC,
                             (wdraw.childs+1)->vparam,
                             (wdraw.childs+1)->org,ONSCREEN);
            SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);
          }
          break;
        case IDM_POPCHANGESECTMULTIELEM:
          if(wdraw.nchilds>=2 &&
             wmenu.nchilds>=3 &&
             (wmenu.childs+2)->vparam.vflag.mv.ftype==F_ORGAN)
          {
            if(globalstatus!=SELECTMULTIELEM ||
               nmultielem==0) break;

            gsect=(*(multielem+0))->sect;

            hpopdlg=CreateDialog(hInstGlobal,
                                 "HOGDLGSECTREGIST",
                                 NULL,
                                 (DLGPROC)DialogProcSectRegist);
            getwindowsize(hpopdlg,&mw,&mh);
            ClientToScreen((wdraw.childs+1)->hwnd,&point);
            MoveWindow(hpopdlg,point.x,point.y,mw,mh,TRUE);
            ShowWindow(hpopdlg,SW_SHOW);

            globalstatus=CHANGESECTMULTIELEM;
          }
          if(wdraw.nchilds>=2 &&
             wmenu.nchilds>=3 &&
             (wmenu.childs+2)->vparam.vflag.mv.ftype==F_ARCLM)
          {
            if(globalstatus!=SELECTMULTIELEM ||
               nmultiwire==0) break;

            gsect=(*(multiwire+0))->sect;

			hpopdlg=CreateDialog(hInstGlobal,
                                 "HOGDLGSECTREGIST",
                                 NULL,
                                 (DLGPROC)DialogProcSectRegist);
            getwindowsize(hpopdlg,&mw,&mh);
            ClientToScreen((wdraw.childs+1)->hwnd,&point);
            MoveWindow(hpopdlg,point.x,point.y,mw,mh,TRUE);
            ShowWindow(hpopdlg,SW_SHOW);

            globalstatus=CHANGESECTMULTIELEM;
          }
          break;
        case IDM_POPMULTIELEMEND:
          if((wmenu.childs+2)->vparam.vflag.mv.ftype==F_ORGAN)
          {
          nmultielem=0;
          free(multielem);
          globalstatus=NEUTRAL;

          /*REDRAW*/
          clearwindow(*(wdraw.childs+1));
          draworganization((wdraw.childs+1)->hdcC,
                           (wdraw.childs+1)->vparam,
                           (wdraw.childs+1)->org,ONSCREEN);
          }
          if((wmenu.childs+2)->vparam.vflag.mv.ftype==F_ARCLM)
          /*UJIOKA FOR BCLNG003*/
          {
          nmultiwire=0;
          free(multiwire);
          globalstatus=NEUTRAL;
          globalcondensationflag=0;

          /*REDRAW*/
          clearwindow(*(wdraw.childs+1));
          drawarclmframe((wdraw.childs+1)->hdcC,
                           (wdraw.childs+1)->vparam,arc,0,ONSCREEN);
          }
          SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);
          break;

        case IDM_POPSELECTMULTINODE:
		  globalstatus=SELECTMULTINODE;
          break;

        case IDM_POPCHANGECONFMULTINODE:
          if(wdraw.nchilds>=2 &&
             wmenu.nchilds>=3 &&
             (wmenu.childs+2)->vparam.vflag.mv.ftype==F_ORGAN)
          {
            if(globalstatus!=SELECTMULTINODE ||
               nmultinode==0) break;
            node=(*(multinode+0));

            if(node==NULL)
            {
              MessageBox(NULL,"Nothing.","MoveNode",MB_OK);
              break;
            }
            pnode=node;

            hpopdlg=CreateDialog(hInstGlobal,
                                 "HOGDLGCHANGECONF",
                                 NULL,
                                 (DLGPROC)DialogProcConf);

            getwindowsize(hpopdlg,&mw,&mh);

            ClientToScreen((wdraw.childs+1)->hwnd,&point);
            MoveWindow(hpopdlg,point.x,point.y,mw,mh,TRUE);
            ShowWindow(hpopdlg,SW_SHOW);

            globalstatus=CHANGECONFMULTINODE;
			changeconfmultinode(&((wdraw.childs+1)->org));

          }
          break;

        case IDM_POPMULTINODEEND:
          nmultinode=0;
          free(multinode);
          globalstatus=NEUTRAL;

          /*REDRAW*/
		  clearwindow(*(wdraw.childs+1));
          draworganization((wdraw.childs+1)->hdcC,
                           (wdraw.childs+1)->vparam,
                           (wdraw.childs+1)->org,ONSCREEN);
          SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);
          break;

        case IDM_POPEND: /*END APPLICATION.*/
          if(MessageBox(hwnd,"End Application.","Hogan",MB_OKCANCEL)
             ==IDOK)
          {
            wparam = MAKEWPARAM((WORD)IDM_END,(WORD)0);
            SendMessage(wmain.hwnd,WM_COMMAND,wparam,(WPARAM)0);
          }
          break;

        default: /*OTHERS.*/
          return DefWindowProc(hwnd,message,wParam,lParam);
      }
      break;

    default:
      return DefWindowProc(hwnd,message,wParam,lParam);
  }
  return 0;
}/*WindowProcedureDraw*/

LRESULT CALLBACK WindowProcedureSurf(HWND hwnd,
                                     UINT message,
                                     WPARAM wParam,
                                     LPARAM lParam)
/*WINDOW FOR YIELD SURFACE.*/
{
  POINT point;
  char str[256];
  long int x,y;

  switch(message)
  {
    case WM_PAINT:
      DefWindowProc(hwnd,message,wParam,lParam);
      overlayhdc(*(wsurf.childs+1),SRCPAINT);
	  break;

    case WM_LBUTTONDOWN:
      x = LOWORD(lParam);
      y = HIWORD(lParam);

      if((wsurf.childs+1)->lstatus==NEUTRAL)
      {
        sprintf(str,"POSITION:%d %d",x,y);
		SetTextColor((wsurf.childs+1)->hdcC,RGB(150,50,150));
        TextOut((wsurf.childs+1)->hdcC,x,y,str,strlen(str));
        SendMessage(hwnd,WM_PAINT,0,0);
      }
      else if((wsurf.childs+1)->lstatus==ROTATE)
      {
        initx=x;
        inity=y;

        initv=(wsurf.childs+1)->vparam;

        globalstatus=ROTATE;
      }
      break;

    case WM_MOUSEMOVE:
      if(globalstatus==ROTATE &&
         (wsurf.childs+1)->lstatus==ROTATE &&
         (wParam == MK_LBUTTON))
      {
        point.x = LOWORD(lParam);
        point.y = HIWORD(lParam);

        x=initx-point.x;
        y=inity-point.y;

        if(labs(x)>=labs(y))
        {
          (wsurf.childs+1)->vparam.theta=initv.theta
                                         +(double)x/10.0;
        }
        else
        {
		  (wsurf.childs+1)->vparam.phi=initv.phi
                                       -(double)y/10.0;
        }
        createviewdata(&((wsurf.childs+1)->vparam));

        clearwindow(*(wsurf.childs+1));
        drawyieldsurface((wsurf.childs+1)->hdcC,
                         (wsurf.childs+1)->vparam,
                         SURFACEX,SURFACEY,SURFACEZ,NULL);
        overlayhdc(*(wsurf.childs+1),SRCPAINT);   /*UPDATE DISPLAY.*/
      }
      else if(globalstatus==ROTATE &&
              (wsurf.childs+1)->lstatus==ROTATE) /*END OF ROTATION.*/
      {
        clearwindow(*(wsurf.childs+1));
        drawyieldsurface((wsurf.childs+1)->hdcC,
                         (wsurf.childs+1)->vparam,
                         SURFACEX,SURFACEY,SURFACEZ,
                         arc.fsurface);
        overlayhdc(*(wsurf.childs+1),SRCPAINT);   /*UPDATE DISPLAY.*/

        globalstatus=NEUTRAL; /*ROTATION END*/
      }
      else
      {
        globalstatus=NEUTRAL; /*ROTATION END*/
      }
      break;

    default:
      return DefWindowProc(hwnd,message,wParam,lParam);
  }
  return 0;
}/*WindowProcedureSurf*/

LRESULT CALLBACK WindowProcedureSframe(HWND hwnd,
                                       UINT message,
                                       WPARAM wParam,
                                       LPARAM lParam)
/*WINDOW PROCEDURE FOR SECTION LIST,SECTION VIEW.*/
{
  WPARAM wparam,lparam;
  double r,a1,a2;

  switch(message)
  {
    case WM_PAINT:
      DefWindowProc(hwnd,message,wParam,lParam);
      break;

    case WM_SIZE:
      DefWindowProc(hwnd,message,wParam,lParam);
      wparam = MAKEWPARAM((WORD)IDM_FITSHEET,(WORD)0);
      hfitfrom=wsfrm.hwnd;
      if(wsdsp.hwnd!=NULL)
      {
        SendMessage((wsdsp.childs+1)->hwnd,WM_COMMAND,wparam,0);
      }
      if(wsect.hwnd!=NULL)
      {
        SendMessage((wsect.childs+1)->hwnd,WM_COMMAND,wparam,0);
      }
      if(wprop.hwnd!=NULL)
      {
        SendMessage((wprop.childs+1)->hwnd,WM_COMMAND,wparam,0);
      }
      break;

    case WM_COMMAND:
      switch(LOWORD(wParam))
      {
        case IDMS_OPEN:
          wparam = MAKEWPARAM((WORD)IDM_PROPERTYLIST,(WORD)0);
          lparam = MAKELPARAM((WORD)0,(WORD)0);
          SendMessage(wmain.hwnd,
                      WM_COMMAND,wparam,lparam);
          wparam = MAKEWPARAM((WORD)IDMP_OPEN,(WORD)0);
          lparam = MAKELPARAM((WORD)0,(WORD)0);
          SendMessage((wprop.childs+1)->hwnd,
                      WM_COMMAND,wparam,lparam);
          wparam = MAKEWPARAM((WORD)IDM_SECTIONLIST,(WORD)0);
          lparam = MAKELPARAM((WORD)0,(WORD)0);
          SendMessage(wmain.hwnd,
                      WM_COMMAND,wparam,lparam);
		  wparam = MAKEWPARAM((WORD)IDMS_OPEN,(WORD)0);
          lparam = MAKELPARAM((WORD)0,(WORD)0);
          SendMessage((wsect.childs+1)->hwnd,
                      WM_COMMAND,wparam,lparam);
          break;
        case IDMS_SAVE:
          wparam = MAKEWPARAM((WORD)IDMS_SAVE,(WORD)0);
		  SendMessage((wsect.childs+1)->hwnd,WM_COMMAND,wparam,0);
          break;

        case IDMP_SIMPLE:
          wparam = MAKEWPARAM((WORD)IDMP_SIMPLE,(WORD)0);
          SendMessage((wprop.childs+1)->hwnd,WM_COMMAND,wparam,0);
          break;
        case IDMP_DETAIL:
          wparam = MAKEWPARAM((WORD)IDMP_DETAIL,(WORD)0);
          SendMessage((wprop.childs+1)->hwnd,WM_COMMAND,wparam,0);
          break;
        case IDMP_ADD:
          wparam = MAKEWPARAM((WORD)IDMP_ADD,(WORD)0);
          SendMessage((wprop.childs+1)->hwnd,WM_COMMAND,wparam,0);
          break;
        case IDMS_SIMPLE:
          wparam = MAKEWPARAM((WORD)IDMS_SIMPLE,(WORD)0);
          SendMessage((wsect.childs+1)->hwnd,WM_COMMAND,wparam,0);
          break;
        case IDMS_DETAIL:
          wparam = MAKEWPARAM((WORD)IDMS_DETAIL,(WORD)0);
          SendMessage((wsect.childs+1)->hwnd,WM_COMMAND,wparam,0);
          break;
        case IDMS_ADD:
          wparam = MAKEWPARAM((WORD)IDMS_ADD,(WORD)0);
          SendMessage((wsect.childs+1)->hwnd,WM_COMMAND,wparam,0);
          break;

        case IDMS_ALLON:
          wparam = MAKEWPARAM((WORD)IDMS_ALLON,(WORD)0);
          SendMessage((wsect.childs+1)->hwnd,WM_COMMAND,wparam,0);
          break;
        case IDMS_ALLOFF:
          wparam = MAKEWPARAM((WORD)IDMS_ALLOFF,(WORD)0);
          SendMessage((wsect.childs+1)->hwnd,WM_COMMAND,wparam,0);
		  break;
        case IDMS_REVERSE:
          wparam = MAKEWPARAM((WORD)IDMS_REVERSE,(WORD)0);
          SendMessage((wsect.childs+1)->hwnd,WM_COMMAND,wparam,0);
          break;

        case IDMS_DRAW:
          wparam = MAKEWPARAM((WORD)IDM_SECTIONVIEW,(WORD)0);
          SendMessage(wmain.hwnd,WM_COMMAND,wparam,0);
          if(wsdsp.nchilds>=2 &&
             (wsdsp.childs+1)->hwnd!=NULL)
          {
            wparam = MAKEWPARAM((WORD)IDMS_VIEW,(WORD)0);
            SendMessage((wsdsp.childs+1)->hwnd,WM_COMMAND,wparam,0);
          }
          break;

        case IDMS_CLEAR:
          freepolypolycurve(&gpolypoly);
          /*freepolycurve(&gpolycurve);*/
          if(wsdsp.nchilds>=2 &&
             (wsdsp.childs+1)->hwnd!=NULL)
          {
            wparam = MAKEWPARAM((WORD)IDMS_VIEW,(WORD)0);
            SendMessage((wsdsp.childs+1)->hwnd,WM_COMMAND,wparam,0);
          }
          break;

        case IDMS_FEATURE:
          /*if(gpolycurve.ncurve>0)
          {*/
          if(gpolypoly.npcurve>0)
          {
            polycurvefeatures((wsdsp.childs+1)->hwnd,
                              &((wsdsp.childs+1)->vparam),
                              &gpolypoly);
          }
          break;
        case IDMS_MP:
          /*INITIALIZE POLYPOLYCURVE*/
          gpolypoly.npcurve=1;
          gpolypoly.pcurves
		  =(struct polycurve *)malloc(sizeof(struct polycurve));

          /*INITIALIZE POLYCURVE*/
          cpolycurve=gpolypoly.pcurves+0;
          cpolycurve->loff  =0;
          cpolycurve->ncurve=0;
          cpolycurve->type  =0;
          cpolycurve->curves=NULL;
          cpolycurve->prop.r=  0;
          cpolycurve->prop.g=  0;
          cpolycurve->prop.b=255;

          cpolycurve->prop.E=2100.0;   /*[tf/cm2]*/
          cpolycurve->prop.poi=1.0/3.0;
          cpolycurve->prop.fut= 2.400; /*[tf/cm2]*/
          cpolycurve->prop.fuc=-2.400;

          /*ADD CURVES*/
          if(addpolycurve.ncurve>0)
          freepolycurve(&addpolycurve);

          addpolycurve.prop.r=100;
          addpolycurve.prop.g=100;
          addpolycurve.prop.b=100;

          /*setpolycurveaskatakou(&addpolycurve,
                                addpolycurve.type,
                                H,B,tf1,tf2,tw,
                                ri1,ri2,ri3,ri4,ro1,ro2,ro3,ro4);*/
          /*addpolycurve.type=CTYPE_O;
          setpolycurveaskatakou(&addpolycurve,
                                addpolycurve.type,
                                0.0,0.0,0.0,0.0,0.0,
                                0.0,0.0,0.0,0.0,
                                100.0,0.0,0.0,0.0);*/ /*[cm]*/
          /*addpolycurve.type=CTYPE_FB;
          setpolycurveaskatakou(&addpolycurve,
                                addpolycurve.type,
                                100.0,100.0,0.0,0.0,0.0,
                                0.0,0.0,0.0,0.0,
                                0.0,0.0,0.0,0.0);*/ /*[cm]*/

addpolycurve.type=CTYPE_O;
addpolycurve.ncurve=3;
addpolycurve.curves=(struct curve *)malloc(addpolycurve.ncurve
                                           *sizeof(struct curve));
r=100.0; a1=-2.00*PI; a2=-1.25*PI;
setcurveascircle((addpolycurve.curves+0),0,1,
                 &r,&a1,&a2,0.0,0.0,
                 0.0,0.0,0.0,0.0);
r=100.0; a1=-1.25*PI; a2=-0.75*PI;
setcurveascircle((addpolycurve.curves+1),0,1,
                 &r,&a1,&a2,0.0,0.0,
                 0.0,0.0,0.0,0.0);
r=100.0; a1=-0.75*PI; a2=0.0*PI;
setcurveascircle((addpolycurve.curves+2),0,1,
                 &r,&a1,&a2,0.0,0.0,
                 0.0,0.0,0.0,0.0);

          addcurvestopolycurve(cpolycurve,&addpolycurve);

          /*DRAW POLYPOLYCURVE*/
          SendMessage((wsdsp.childs+1)->hwnd,WM_PAINT,0,0);

          ultimatebendingofpolycurve((wsdsp.childs+1)->hwnd,
                                     &((wsdsp.childs+1)->vparam),
                                     -10000.0,/*Nu[tf]*/
                                     &gpolypoly);
          break;

        default: /*OTHERS.*/
          return DefWindowProc(hwnd,message,wParam,lParam);
      }
      break;

    default:
      return DefWindowProc(hwnd,message,wParam,lParam);
  }
  return 0;
}/*WindowProcedureSframe*/

LRESULT CALLBACK WindowProcedureProp(HWND hwnd,
                                     UINT message,
                                     WPARAM wParam,
									 LPARAM lParam)
/*WINDOW FOR PROPERTY LIST.*/
{
  HWND hdlg,hoya;
  WPARAM wparam;
  LPARAM lparam;
  SIZE tsize;
  char str[256];
  long int i;
  long int dw,dh,wmax,hmax,maxX,maxY,mw,mh,sw,sh,pw,ph;

  switch(message)
  {
    case WM_PAINT:
      DefWindowProc(hwnd,message,wParam,lParam);

      drawtexts((wprop.childs+1)->hdcC,
                wprop.strset,
                wprop.nstring,
                (wprop.childs+1)->vparam);

      overlayhdc(*(wprop.childs+1),SRCCOPY/*SRCAND*/);
      break;

    case WM_COMMAND:
      switch(LOWORD(wParam))
      {
        case IDM_FITSHEET: /*NEVER USE FITPARENT.*/
          hoya=GetParent(GetParent(hwnd)); /*PARENT SHEET.*/

          getclientsize(wsfrm.hwnd,&maxX,&maxY);
          getwindowsize(wprop.hwnd,&pw,&ph);
          getwindowsize(wsect.hwnd,&sw,&sh);
          if(hfitfrom==wsfrm.hwnd)
          {
            if(maxX<100) mw=50;
            else mw=pw;

            if(maxY<100) mh=50;
            else mh=ph;
          }
          else if(hfitfrom==(wsect.childs+1)->hwnd)
		  {
            mw=sw;

            mh=maxY-sh;
            if(mh<50) mh=50;
          }
          else if(hfitfrom==(wsdsp.childs+1)->hwnd)
          {
            getwindowsize(wsdsp.hwnd,&mw,&mh);

            mw=maxX-mw;
            if(mw<50) mw=50;

            mh=maxY-sh;
            if(mh<50) mh=50;
          }
          else break;

          if(pw!=mw || ph!=mh)
          {
            MoveWindow(hoya,0,0,mw,mh,TRUE);
          }
          break;

        case IDM_FITPARENT: /*MESSAGE FROM PARENT SHEET.*/
          wparam = MAKEWPARAM((WORD)IDM_FITSHEET,(WORD)0);
          /*if(wsdsp.hwnd!=NULL)
          {
            hfitfrom=(wprop.childs+1)->hwnd;
            SendMessage((wsdsp.childs+1)->hwnd,WM_COMMAND,wparam,0);
          }*/
          if(wsect.hwnd!=NULL)
          {
            hfitfrom=(wprop.childs+1)->hwnd;
            SendMessage((wsect.childs+1)->hwnd,WM_COMMAND,wparam,0);
          }
          break;

        case IDMP_OPEN: /*OPEN PROPERTY LIST.*/
          if(gnprop>0)
          {
            for(i=0;i<gnprop;i++)
			{
              ReleaseDC((wprop.childs+i+2)->hwnd,
                        (wprop.childs+i+2)->hdcB);
              ReleaseDC((wprop.childs+i+2)->hwnd,
                        (wprop.childs+i+2)->hdcC);
              DestroyWindow((wprop.childs+i+2)->hwnd);
            }

            wprop.childs=(struct windowparams *)
                         realloc(wprop.childs,
                                 2*sizeof(struct windowparams));
            wprop.nchilds=2;
          }

          gnprop=0;
          if((wmenu.childs+2)!=NULL)
          {
            if((wmenu.childs+2)->vparam.vflag.mv.ftype==F_ORGAN)
            {
              gnprop=(wdraw.childs+1)->org.nprop;
              currentprops=(wdraw.childs+1)->org.props;
            }
          }

          if(gnprop>0)
          {
            if(wprop.gstatus==NEUTRAL) wprop.gstatus=PROP_SIMPLE;

            wmax=0;
            hmax=0;
            for(i=0;i<gnprop;i++)
            {
              gprop=currentprops+i;
              sprintf(gstr,"%ld",gprop->code);

              hdlg = CreateDialog(hInstGlobal,
                                  "HOGDLGPROPSIMPLE",
                                  (wprop.childs+1)->hwnd,
                                  (DLGPROC)DialogProcPsim);

              wprop.childs=(struct windowparams *)
                           realloc(wprop.childs,
								   (wprop.nchilds+1)
                                   *sizeof(struct windowparams));
              (wprop.childs+(wprop.nchilds))->hwnd=hdlg;

              getwindowsize(hdlg,&dw,&dh);
              MoveWindow(hdlg,0,hmax,dw,dh,TRUE);

              if(wmax<dw) wmax=dw;
              hmax+=dh;

              wprop.nchilds++;

              if(LOWORD(lParam)==(WORD)PROPDETAIL)
              {
                if((wmenu.childs+2)->vparam.vflag.mv.ftype==F_ORGAN)
                {
				  sprintf(str,"Hiju=%.3f",gprop->hiju);
                  GetTextExtentPoint32((wprop.childs+1)->hdcC,
                                       str,strlen(str),&tsize);
                  wprop.strset=addtext(wprop.strset,
                                       &(wprop.nstring),
                                       str,50.0,(double)hmax);
                  hmax+=tsize.cy;

                  sprintf(str,"Young E=%.3f",gprop->E);
                  wprop.strset=addtext(wprop.strset,
                                       &(wprop.nstring),
                                       str,50.0,(double)hmax);
                  hmax+=tsize.cy;

                  sprintf(str,"Poisson v=%.5f",gprop->poi);
                  wprop.strset=addtext(wprop.strset,
                                       &(wprop.nstring),
                                       str,50.0,(double)hmax);
                  hmax+=tsize.cy;
                }
              }
            }
            MoveWindow((wprop.childs+1)->hwnd,0,0,wmax,hmax+5,TRUE);
          }
          break;

		case IDMP_SIMPLE: /*SIMPLE TYPE.*/
          wprop.gstatus=PROP_SIMPLE;
          wparam = MAKEWPARAM((WORD)IDMP_OPEN,(WORD)0);
          lparam = MAKELPARAM((WORD)PROPSIMPLE,(WORD)0);
          SendMessage((wprop.childs+1)->hwnd,
                      WM_COMMAND,wparam,lparam);
          break;

        case IDMP_DETAIL: /*DETAILED TYPE.*/
          wprop.gstatus=PROP_DETAIL;
          wparam = MAKEWPARAM((WORD)IDMP_OPEN,(WORD)0);
          lparam = MAKELPARAM((WORD)PROPDETAIL,(WORD)0);
          SendMessage((wprop.childs+1)->hwnd,
                      WM_COMMAND,wparam,lparam);
		  break;

        case IDMP_ADD: /*ADD PROPERTY*/
          cpolycurve=NULL;
          globalstatus=C_ADDPROPERTY;

          gprop=(struct oprop *)malloc(sizeof(struct oprop));
          gprop->code=0;
          /*gprop->name=NULL;*/
          gprop->name=(char *)malloc(1*sizeof(char));
          strcpy(gprop->name,"\0");
          gprop->E   =0.0;
          gprop->poi =0.0;
          gprop->hiju=0.0;
          gprop->r=0;
          gprop->g=0;
          gprop->b=0;
          addingprop=gprop;

          ipopflag=0;

          hpopdlg=CreateDialog(hInstGlobal,
                               "HOGDLGPROPERTY",
                               NULL,
                               (DLGPROC)DialogProcProperty);

          ShowWindow(hpopdlg,SW_SHOW);
		  break;
		case IDM_POPPROPERTYRETURN:
          /*ADD PROPERTY*/
          addproperty(addingprop,&((wdraw.childs+1)->org));

          globalstatus=NEUTRAL;

          /*OPEN AGAIN*/
          if(wprop.gstatus==PROP_SIMPLE)
          {
            lparam = MAKELPARAM((WORD)PROPSIMPLE,(WORD)0);
          }
          if(wprop.gstatus==PROP_DETAIL)
          {
            lparam = MAKELPARAM((WORD)PROPDETAIL,(WORD)0);
          }
          wparam = MAKEWPARAM((WORD)IDMP_OPEN,(WORD)0);
          SendMessage((wprop.childs+1)->hwnd,
                      WM_COMMAND,wparam,lparam);
          break;

        default: /*OTHERS.*/
          return DefWindowProc(hwnd,message,wParam,lParam);
      }
      break;

    default:
      return DefWindowProc(hwnd,message,wParam,lParam);
  }
  return 0;
}/*WindowProcedureProp*/

LRESULT CALLBACK WindowProcedureSect(HWND hwnd,
                                     UINT message,
                                     WPARAM wParam,
                                     LPARAM lParam)
/*WINDOW FOR SECTION LIST.*/
{
  HWND hdlg,hoya;
  WPARAM wparam;
  LPARAM lparam;
  SIZE tsize;
  char /*s[80],*/str[400];
  long int i,j;
  long int dw,dh,wmax,hmax,maxX,maxY,mw,mh,sw,sh,pw,ph;
  struct features f;
  struct osect *os;

  switch(message)
  {
    case WM_PAINT:
      DefWindowProc(hwnd,message,wParam,lParam);

      drawtexts((wsect.childs+1)->hdcC,
                wsect.strset,
                wsect.nstring,
                (wsect.childs+1)->vparam);

      overlayhdc(*(wsect.childs+1),SRCCOPY/*SRCAND*/);
      break;

    case WM_COMMAND:
      switch(LOWORD(wParam))
      {
        case IDM_FITSHEET: /*NEVER USE FITPARENT.*/
          hoya=GetParent(GetParent(hwnd)); /*PARENT SHEET.*/

          getclientsize(wsfrm.hwnd,&maxX,&maxY);
          getwindowsize(wprop.hwnd,&pw,&ph);
          getwindowsize(wsect.hwnd,&sw,&sh);
          if(hfitfrom==wsfrm.hwnd)
          {
            mw=pw;

            mh=maxY-ph;
            if(mh<=50)
            {
              mh=50;
              ph=maxY-50;
              if(ph<50) ph=50;
            }
          }
          else if(hfitfrom==(wprop.childs+1)->hwnd)
          {
            mw=pw;

            mh=maxY-ph;
            if(mh<=50)
            {
              mh=50;
              ph=maxY-50;
              if(ph<50) ph=50;
            }
          }
          else if(hfitfrom==(wsdsp.childs+1)->hwnd)
          {
            getwindowsize(wsdsp.hwnd,&mw,&mh);

            mw=maxX-mw;
            if(mw<50) mw=50;

            mh=maxY-ph;
            if(mh<=50)
            {
              mh=50;
              ph=maxY-50;
              if(ph<50) ph=50;
            }
          }
          else break;

          MoveWindow(hoya,0,ph,mw,mh,TRUE);
          break;

        case IDM_FITPARENT: /*MESSAGE FROM PARENT SHEET.*/
          wparam = MAKEWPARAM((WORD)IDM_FITSHEET,(WORD)0);
          /*if(wsdsp.hwnd!=NULL)
          {
            hfitfrom=(wsect.childs+1)->hwnd;
            SendMessage((wsdsp.childs+1)->hwnd,WM_COMMAND,wparam,0);
          }*/
          if(wprop.hwnd!=NULL)
          {
            hfitfrom=(wsect.childs+1)->hwnd;
            SendMessage((wprop.childs+1)->hwnd,WM_COMMAND,wparam,0);
          }
          break;

        case IDMS_OPEN: /*OPEN SECTION LIST.*/
          if(gnsect>0)
          {
            for(i=0;i<gnsect;i++)
            {
              ReleaseDC((wsect.childs+i+2)->hwnd,
                        (wsect.childs+i+2)->hdcB);
              ReleaseDC((wsect.childs+i+2)->hwnd,
                        (wsect.childs+i+2)->hdcC);
              DestroyWindow((wsect.childs+i+2)->hwnd);
            }

            wsect.childs=(struct windowparams *)
                         realloc(wsect.childs,
                                 2*sizeof(struct windowparams));
            wsect.nchilds=2;
          }

          gnsect=0;
          if((wmenu.childs+2)!=NULL)
          {
            if((wmenu.childs+2)->vparam.vflag.mv.ftype==F_ARCLM ||
               (wmenu.childs+2)->vparam.vflag.mv.ftype==F_FRAME)
            {
              gnsect=arc.nsect;
              currentsects=arc.sects;
            }
            if((wmenu.childs+2)->vparam.vflag.mv.ftype==F_ORGAN)
            {
              gnsect=(wdraw.childs+1)->org.nsect;
              currentsects=(wdraw.childs+1)->org.sects;
            }
          }

          if(gnsect>0)
          {
            if(wsect.gstatus==NEUTRAL) wsect.gstatus=SECT_SIMPLE;

            wmax=0;
            hmax=0;
            for(i=0;i<gnsect;i++)
			{
              os   =currentsects+i;
              gsect=currentsects+i;
              sprintf(gstr,"%ld",gsect->code);

              if(wsdsp.nchilds<2)
              {
                wparam = MAKEWPARAM((WORD)IDM_SECTIONVIEW,(WORD)0);
                SendMessage(wmain.hwnd,WM_COMMAND,wparam,0);
              }

              if(wsdsp.nchilds>=2 && (gsect->ppc.npcurve)>0)
              {
                gsect->ppc.hico
                =createsectionicon((wsect.childs+1)->hdcC,32,32,
                                   (wsdsp.childs+1)->vparam,
                                   &(gsect->ppc));
              }

              hdlg = CreateDialog(hInstGlobal,
                                  "HOGDLGSECTSIMPLE",
                                  (wsect.childs+1)->hwnd,
                                  (DLGPROC)DialogProcSsim);

              wsect.childs=(struct windowparams *)
                           realloc(wsect.childs,
                                   (wsect.nchilds+1)
                                   *sizeof(struct windowparams));
              (wsect.childs+(wsect.nchilds))->hwnd=hdlg;

              getwindowsize(hdlg,&dw,&dh);
              MoveWindow(hdlg,0,hmax,dw,dh,TRUE);

              if(wmax<dw) wmax=dw;
              hmax+=dh;

              wsect.nchilds++;

              if(LOWORD(lParam)==(WORD)SECTDETAIL)
              {
                if((wmenu.childs+2)->vparam.vflag.mv.ftype==F_ORGAN)
                {
				  /*EXTEND WINDOW.*/
                  /*
                  GetTextExtentPoint32((wsect.childs+1)->hdcC,
                                       "Size",strlen("Size"),&tsize);
                  (wsect.childs+1)->hdcC
                  =extendhdc((wsect.childs+1)->hdcC,
                             wmax,(hmax+(gsect->nfig+6)*tsize.cy));
                  (wsect.childs+1)->hdcB
                  =extendhdc((wsect.childs+1)->hdcB,
                             wmax,(hmax+(gsect->nfig+6)*tsize.cy));
                  */

                  /*STRINGS.*/
                  if(strlen(os->name)>0)
                  {
                    GetTextExtentPoint32((wsect.childs+1)->hdcC,
                                         os->name,
                                         strlen(os->name),&tsize);
                    wsect.strset=addtext(wsect.strset,
                                         &(wsect.nstring),
                                         os->name,
                                         50.0,(double)hmax);
                    hmax+=tsize.cy;
                  }

                  if(os->role!=ROLEHOJO)
                  {
                    for(j=0;j<(os->nfig);j++)
                    {
                      if(strlen((os->figs+j)->prop->name)>0 &&
                         (os->figs+j)->thick>0.0)
                      {
                        sprintf(str," %s %.0fmm",
                                (os->figs+j)->prop->name,
                                (os->figs+j)->thick*1000.0);

                        GetTextExtentPoint32((wsect.childs+1)->hdcC,
                                             str,strlen(str),&tsize);
                        wsect.strset=addtext(wsect.strset,
                                             &(wsect.nstring),
                                             str,50.0,(double)hmax);
                        hmax+=tsize.cy;
					  }
                    }

                    f=sectionfeatures(os);
                    os->E   =f.E;
                    os->poi =f.poi;
                    os->area=f.Ax;
                    os->Ixx =f.Igx;
                    os->Iyy =f.Igy;
                    os->Jzz =f.Jzz;
					os->hiju[0]=f.hiju;
                    os->hiju[1]=f.hiju;
                    os->hiju[2]=f.hiju;

                    hmax+=tsize.cy/3;
                    sprintf(str,"A=%.4f m",os->area);
                    GetTextExtentPoint32((wsect.childs+1)->hdcC,
                                         str,strlen(str),&tsize);
                    wsect.strset=addtext(wsect.strset,
                                         &(wsect.nstring),
                                         str,50.0,(double)hmax);
                    sprintf(str,"2");
                    wsect.strset=addtext(wsect.strset,
                                         &(wsect.nstring),
                                         str,
                                         (double)(50+tsize.cx),
                                         (double)(hmax-tsize.cy/3));
                    hmax+=tsize.cy;

                    sprintf(str,"Ixx=%.8f m",os->Ixx);
                    GetTextExtentPoint32((wsect.childs+1)->hdcC,
                                         str,strlen(str),&tsize);
                    wsect.strset=addtext(wsect.strset,
                                         &(wsect.nstring),
                                         str,50.0,(double)hmax);
                    sprintf(str,"4");
                    wsect.strset=addtext(wsect.strset,
                                         &(wsect.nstring),
                                         str,
                                         (double)(50+tsize.cx),
                                         (double)(hmax-tsize.cy/3));
                    hmax+=tsize.cy;

                    sprintf(str,"Iyy=%.8f m",os->Iyy);
                    GetTextExtentPoint32((wsect.childs+1)->hdcC,
                                         str,strlen(str),&tsize);
                    wsect.strset=addtext(wsect.strset,
                                         &(wsect.nstring),
                                         str,50.0,(double)hmax);
                    sprintf(str,"4");
                    wsect.strset=addtext(wsect.strset,
                                         &(wsect.nstring),
                                         str,
                                         (double)(50+tsize.cx),
                                         (double)(hmax-tsize.cy/3));
                    hmax+=tsize.cy;

                    sprintf(str,"Jzz=%.8f m",os->Jzz);
                    GetTextExtentPoint32((wsect.childs+1)->hdcC,
                                         str,strlen(str),&tsize);
                    wsect.strset=addtext(wsect.strset,
                                         &(wsect.nstring),
                                         str,50.0,(double)hmax);
                    sprintf(str,"4");
                    wsect.strset=addtext(wsect.strset,
                                         &(wsect.nstring),
                                         str,
                                         (double)(50+tsize.cx),
                                         (double)(hmax-tsize.cy/3));
                    hmax+=tsize.cy;

                    sprintf(str,"Hiju=%.3f %.3f %.3f tf/m",
                            os->hiju[0],os->hiju[1],os->hiju[2]);
                    wsect.strset=addtext(wsect.strset,
                                         &(wsect.nstring),
                                         str,50.0,(double)hmax);
                    hmax+=tsize.cy;
/*
 sprintf(str,"\0");
 sprintf(s,"Code:%ld Address:%ld\n",os->code,os);
 strcat(str,s);
 sprintf(s,"E=%.3f Poi=%.3f\n",os->E,os->poi);
 strcat(str,s);
 sprintf(s,"Hiju=%.3f\n",os->hiju);
 strcat(str,s);
 sprintf(s,"A=%.6f\n",os->area);
 strcat(str,s);
 sprintf(s,"Ixx=%.8f Iyy=%.8f\n",os->Ixx,os->Iyy);
 strcat(str,s);
 sprintf(s,"Jzz=%.8f\n",os->Jzz);
 strcat(str,s);
 MessageBox(NULL,str,"After Section Features",MB_OK);
*/
                  }
                }
              }
            }
            MoveWindow((wsect.childs+1)->hwnd,0,0,wmax,hmax+5,TRUE);
          }
          break;

        case IDMS_SIMPLE: /*SIMPLE TYPE.*/
          wsect.gstatus=SECT_SIMPLE;
          wparam = MAKEWPARAM((WORD)IDMS_OPEN,(WORD)0);
          lparam = MAKELPARAM((WORD)SECTSIMPLE,(WORD)0);
          SendMessage((wsect.childs+1)->hwnd,
                      WM_COMMAND,wparam,lparam);
          break;

        case IDMS_DETAIL: /*DETAILED TYPE.*/
          wsect.gstatus=SECT_DETAIL;
          wparam = MAKEWPARAM((WORD)IDMS_OPEN,(WORD)0);
          lparam = MAKELPARAM((WORD)SECTDETAIL,(WORD)0);
          SendMessage((wsect.childs+1)->hwnd,
                      WM_COMMAND,wparam,lparam);
          break;

        case IDMS_ADD: /*ADD SECTION*/
          globalstatus=C_ADDSECTION;

          gsect=(struct osect *)malloc(sizeof(struct osect));
          gsect->code=0;
          /*gsect->name=NULL;*/
          gsect->name=(char *)malloc(1*sizeof(char));
          strcpy(gsect->name,"\0");
          gsect->dflag=1;
		  gsect->nfig=0;
          gsect->area=0.0;
          gsect->dcolor.r=255;
          gsect->dcolor.g=255;
          gsect->dcolor.b=255;
          gsect->role=ROLENULL;
          gsect->ppc.npcurve=0;
          addingsect=gsect;

          hpopdlg=CreateDialog(hInstGlobal,
                               "HOGDLGSECTREGIST",
                               NULL,
                               (DLGPROC)DialogProcSectRegist);
          ShowWindow(hpopdlg,SW_SHOW);
          break;
        case IDM_POPREGISTRETURN:
          /*ADD SECTION*/
          GetDlgItemText(hpopdlg,IDSR_POPCODE,str,20);
          addingsect->code=strtol(str,NULL,10);
          GetDlgItemText(hpopdlg,IDSR_POPCAPTION,
                         addingsect->ppc.name,80);

          addsection(addingsect,&((wdraw.childs+1)->org));

          globalstatus=NEUTRAL;

          /*OPEN AGAIN*/
          if(wsect.gstatus==SECT_SIMPLE)
          {
            lparam = MAKELPARAM((WORD)SECTSIMPLE,(WORD)0);
          }
          if(wsect.gstatus==SECT_DETAIL)
          {
            lparam = MAKELPARAM((WORD)SECTDETAIL,(WORD)0);
          }
          wparam = MAKEWPARAM((WORD)IDMS_OPEN,(WORD)0);
          SendMessage((wsect.childs+1)->hwnd,
                      WM_COMMAND,wparam,lparam);
          break;

        case IDMS_ALLON: /*ALL SECTIONS ON.*/
          if((wmenu.childs+2)->vparam.vflag.mv.ftype==F_ARCLM ||
			 (wmenu.childs+2)->vparam.vflag.mv.ftype==F_FRAME)
          {
            if(arc.nsect>0)
            {
              for(i=0;i<arc.nsect;i++)
              {
                (arc.sects+i)->dflag=1;
              }
            }
            else break;
          }
          else break;

          if(wsect.gstatus==NEUTRAL || wsect.gstatus==SECT_SIMPLE)
          {
            wparam = MAKEWPARAM((WORD)IDMS_OPEN,(WORD)0);
            lparam = MAKELPARAM((WORD)SECTSIMPLE,(WORD)0);
            SendMessage((wsect.childs+1)->hwnd,
                        WM_COMMAND,wparam,lparam);
          }
          else if(wsect.gstatus==SECT_DETAIL)
          {
            wparam = MAKEWPARAM((WORD)IDMS_OPEN,(WORD)0);
            lparam = MAKELPARAM((WORD)SECTDETAIL,(WORD)0);
            SendMessage((wsect.childs+1)->hwnd,
                        WM_COMMAND,wparam,lparam);
          }
          break;

        case IDMS_ALLOFF: /*ALL SECTIONS OFF.*/
          if((wmenu.childs+2)->vparam.vflag.mv.ftype==F_ARCLM ||
             (wmenu.childs+2)->vparam.vflag.mv.ftype==F_FRAME)
          {
            if(arc.nsect>0)
            {
              for(i=0;i<arc.nsect;i++)
              {
                (arc.sects+i)->dflag=0;
              }
            }
            else break;
          }
		  else break;

          if(wsect.gstatus==NEUTRAL || wsect.gstatus==SECT_SIMPLE)
          {
            wparam = MAKEWPARAM((WORD)IDMS_OPEN,(WORD)0);
            lparam = MAKELPARAM((WORD)SECTSIMPLE,(WORD)0);
            SendMessage((wsect.childs+1)->hwnd,
                        WM_COMMAND,wparam,lparam);
          }
          else if(wsect.gstatus==SECT_DETAIL)
          {
            wparam = MAKEWPARAM((WORD)IDMS_OPEN,(WORD)0);
            lparam = MAKELPARAM((WORD)SECTDETAIL,(WORD)0);
            SendMessage((wsect.childs+1)->hwnd,
                        WM_COMMAND,wparam,lparam);
          }
          break;

        case IDMS_REVERSE: /*SECTIONS ON,OFF REVERSE.*/
          if((wmenu.childs+2)->vparam.vflag.mv.ftype==F_ARCLM ||
             (wmenu.childs+2)->vparam.vflag.mv.ftype==F_FRAME)
          {
            if(arc.nsect>0)
            {
              for(i=0;i<arc.nsect;i++)
              {
                if((arc.sects+i)->dflag==0) (arc.sects+i)->dflag=1;
                else                        (arc.sects+i)->dflag=0;
              }
            }
            else break;
          }
          else break;

          if(wsect.gstatus==NEUTRAL || wsect.gstatus==SECT_SIMPLE)
          {
            wparam = MAKEWPARAM((WORD)IDMS_OPEN,(WORD)0);
            lparam = MAKELPARAM((WORD)SECTSIMPLE,(WORD)0);
            SendMessage((wsect.childs+1)->hwnd,
                        WM_COMMAND,wparam,lparam);
          }
          else if(wsect.gstatus==SECT_DETAIL)
		  {
            wparam = MAKEWPARAM((WORD)IDMS_OPEN,(WORD)0);
            lparam = MAKELPARAM((WORD)SECTDETAIL,(WORD)0);
            SendMessage((wsect.childs+1)->hwnd,
                        WM_COMMAND,wparam,lparam);
          }
          break;

        default: /*OTHERS.*/
          return DefWindowProc(hwnd,message,wParam,lParam);
      }
      break;

    default:
      return DefWindowProc(hwnd,message,wParam,lParam);
  }
  return 0;
}/*WindowProcedureSect*/

LRESULT CALLBACK WindowProcedureSview(HWND hwnd,
                                      UINT message,
                                      WPARAM wParam,
                                      LPARAM lParam)
/*WINDOW FOR SECTION VIEW.*/
{
  HWND hoya;
  HCURSOR hcursor;
  POINT point;
  WPARAM wparam;
  char /*s[80],*/str[256];
  int i,ii,ncurve;
  long int x,y,maxX,maxY,mw,mh,sw,sh,pw,ph,code;
  double x1,y1,fact;
  double dmx,dmy;
  double H,B,tf1,tf2,tw,ri1,ri2,ri3,ri4,ro1,ro2,ro3,ro4;
  struct plane pl;
  struct viewparam *vp;
  struct onode nadd,*node;
  struct polycurve *pc;

  switch(message)
  {
	case WM_PAINT:
      DefWindowProc(hwnd,message,wParam,lParam);

      clearwindow(*(wsdsp.childs+1));
      drawglobalpolycurve((wsdsp.childs+1)->hdcC,
                          (wsdsp.childs+1)->vparam,
                          &gpolypoly);
      overlayhdc(*(wsdsp.childs+1),SRCPAINT);
      break;

    case WM_LBUTTONDOWN:
      x = LOWORD(lParam);
      y = HIWORD(lParam);
      point.x=x;
      point.y=y;

      prestatus=globalstatus;

      if(globalstatus==C_SECTIONLINE ||
         globalstatus==C_SECTIONCIRCLE) /*CREATING.*/
      {
        if((wsdsp.childs+1)->lstatus==SELECTNODE)
        {
          node=selectpolycurvenode((wsdsp.childs+1)->vparam,
                                   &gpolypoly,point);
          if(node==NULL)
          {
            MessageBox(NULL,"Nothing.","Node",MB_OK);
            break;
          }

          hcursor = LoadCursor(hInstGlobal,"CANCURSORW");
          SetClassLong((wsdsp.childs+1)->hwnd,
                       GCL_HCURSOR,(LONG)hcursor);

          /*prestatus=SELECTNODE;*/
        }
        else /*CREATE NODE ON GROUND.*/
        {
          pl.nods[0]=setcoord(0.0,0.0,0.0);

          pl.nvec.dc[GX]=0.0;
		  pl.nvec.dc[GY]=0.0;
          pl.nvec.dc[GZ]=1.0;

          pl.a=0.0; pl.b=0.0; pl.c=1.0; pl.d=0.0; /*PLANE:Z=0*/

          vp=&((wsdsp.childs+1)->vparam);
          dmx=(double)( x-(vp->Xo));
          dmy=(double)(-y+(vp->Yo));
          createnodeonplane(*vp,dmx,dmy,pl,&nadd);

          if(nadd.code==0)
          {
            MessageBox(NULL,"Failed.","Node",MB_OK);
            break;
          }
          /*
          nadd.code=code;
          nadd.loff=loff;
          */
        }

        if(gpolypoly.npcurve==0)
        {
          gpolypoly.npcurve=1;
          gpolypoly.pcurves
          =(struct polycurve *)malloc(sizeof(struct polycurve));

          cpolycurve=gpolypoly.pcurves+0;
          cpolycurve->loff  =0;
          cpolycurve->ncurve=0;
          cpolycurve->type  =0;
          cpolycurve->curves=NULL;
          cpolycurve->prop.r=  0;
          cpolycurve->prop.g=  0;
          cpolycurve->prop.b=255;
        }

        if(globalstatus==C_SECTIONLINE) /*CREATE LINE.*/
        {
          if(icount==0)
          {
            cpolycurve->ncurve++;
			ncurve=cpolycurve->ncurve;

            cpolycurve->curves=(struct curve *)
                              realloc(cpolycurve->curves,
                                      ncurve
                                     *sizeof(struct curve));
            initializecurve(cpolycurve->curves+ncurve-1);
            (cpolycurve->curves+ncurve-1)->loff=ncurve-1;
            (cpolycurve->curves+ncurve-1)->type=CTYPE_LINE;
            /*(cpolycurve->curves+ncurve-1)->hugo=1;*/
            (cpolycurve->curves+ncurve-1)->dots[0]
            =(struct onode *)malloc(sizeof(struct onode));
            (cpolycurve->curves+ncurve-1)->dots[1]
            =(struct onode *)malloc(sizeof(struct onode));

            if((wsdsp.childs+1)->lstatus==SELECTNODE)
            {
              *((cpolycurve->curves+ncurve-1)->dots[0])=*node;
              (wsdsp.childs+1)->lstatus=NEUTRAL;
            }
            else *((cpolycurve->curves+ncurve-1)->dots[0])=nadd;

            icount++;
          }
          else if(icount==1)
          {
            ncurve=cpolycurve->ncurve;
            if((wsdsp.childs+1)->lstatus==SELECTNODE)
            {
              *((cpolycurve->curves+ncurve-1)->dots[1])=*node;
              (wsdsp.childs+1)->lstatus=NEUTRAL;
            }
            else *((cpolycurve->curves+ncurve-1)->dots[1])=nadd;
            icount=0;
            globalstatus=NEUTRAL;
          }
        }
        else if(globalstatus==C_SECTIONCIRCLE) /*CREATE CIRCLE.*/
        {
          if(icount==0)
          {
            cpolycurve->ncurve++;
			ncurve=cpolycurve->ncurve;

            cpolycurve->curves=(struct curve *)
                               realloc(cpolycurve->curves,
                                       ncurve
                                       *sizeof(struct curve));
            initializecurve(cpolycurve->curves+ncurve-1);
            (cpolycurve->curves+ncurve-1)->loff=ncurve-1;
            (cpolycurve->curves+ncurve-1)->type=CTYPE_LINE;
            (cpolycurve->curves+ncurve-1)->hugo=1;
            (cpolycurve->curves+ncurve-1)->center
			=(struct onode *)malloc(sizeof(struct onode));
            (cpolycurve->curves+ncurve-1)->dots[0]
            =(struct onode *)malloc(sizeof(struct onode));
            (cpolycurve->curves+ncurve-1)->dots[1]
            =(struct onode *)malloc(sizeof(struct onode));

            if((wsdsp.childs+1)->lstatus==SELECTNODE)
            {
              *((cpolycurve->curves+ncurve-1)->center)=*node;
              (wsdsp.childs+1)->lstatus=NEUTRAL;
            }
            else *((cpolycurve->curves+ncurve-1)->center)=nadd;

            *((cpolycurve->curves+ncurve-1)->dots[0])
            =*((cpolycurve->curves+ncurve-1)->center);

            icount++;
          }
          else if(icount==1)
          {
            ncurve=cpolycurve->ncurve;
            if((wsdsp.childs+1)->lstatus==SELECTNODE)
            {
              *((cpolycurve->curves+ncurve-1)->dots[0])=*node;
              (wsdsp.childs+1)->lstatus=NEUTRAL;
            }
            else *((cpolycurve->curves+ncurve-1)->dots[0])=nadd;

            (cpolycurve->curves+ncurve-1)->radius[0]
            =distancedotdot(*((cpolycurve->curves+ncurve-1)
                            ->center),
							*((cpolycurve->curves+ncurve-1)
                            ->dots[0]));

            (cpolycurve->curves+ncurve-1)->angle[0]
            =definecircleangle((cpolycurve->curves+ncurve-1),0,
             (cpolycurve->curves+ncurve-1)->dots[0]->d[GX],
             (cpolycurve->curves+ncurve-1)->dots[0]->d[GY]);

            (cpolycurve->curves+ncurve-1)->type=CTYPE_CIRCLE;
            icount++;
          }
          else if(icount==2)
          {
            ncurve=cpolycurve->ncurve;
            if((wsdsp.childs+1)->lstatus==SELECTNODE)
            {
              *((cpolycurve->curves+ncurve-1)->dots[1])=*node;
              (wsdsp.childs+1)->lstatus=NEUTRAL;
            }
            else *((cpolycurve->curves+ncurve-1)->dots[1])=nadd;

            (cpolycurve->curves+ncurve-1)->angle[1]
            =definecircleangle((cpolycurve->curves+ncurve-1),1,
             (cpolycurve->curves+ncurve-1)->dots[1]->d[GX],
             (cpolycurve->curves+ncurve-1)->dots[1]->d[GY]);

            x1=(cpolycurve->curves+ncurve-1)->center->d[GX];
            y1=(cpolycurve->curves+ncurve-1)->center->d[GY];

            (cpolycurve->curves+ncurve-1)->dots[1]->d[GX]
            =x1
             +(cpolycurve->curves+ncurve-1)->radius[0]
             *cos((cpolycurve->curves+ncurve-1)->angle[1]);
            (cpolycurve->curves+ncurve-1)->dots[1]->d[GY]
            =y1
             +(cpolycurve->curves+ncurve-1)->radius[0]
             *sin((cpolycurve->curves+ncurve-1)->angle[1]);
/*
sprintf(str,"\0");
sprintf(s,"r=%.3f c=(%.3f %.3f) d1=(%.3f %.3f) d2=(%.3f %.3f)\n",
        (gpolycurve.curves+ncurve-1)->radius[0],
        (gpolycurve.curves+ncurve-1)->center->d[GX],
		(gpolycurve.curves+ncurve-1)->center->d[GY],
        (gpolycurve.curves+ncurve-1)->dots[0]->d[GX],
        (gpolycurve.curves+ncurve-1)->dots[0]->d[GY],
        (gpolycurve.curves+ncurve-1)->dots[1]->d[GX],
        (gpolycurve.curves+ncurve-1)->dots[1]->d[GY]);
strcat(str,s);
sprintf(s,"a1=%.3f a2=%.3f [DEGREE]",
        ((gpolycurve.curves+ncurve-1)->angle[0])/PI*180.0,
        ((gpolycurve.curves+ncurve-1)->angle[1])/PI*180.0);
strcat(str,s);
MessageBox(NULL,str,"Circle",MB_OK);
*/
            icount=0;
            globalstatus=NEUTRAL;
          }
        }
      }
      else if(globalstatus==C_PICKPOLYCURVE)
      {
        pc=pickpolycurve(&gpolypoly,
                         (wsdsp.childs+1)->vparam,x,y);
        if(pc!=NULL) cpolycurve=pc;
        else MessageBox(NULL,"Nothing.","Pick",MB_OK);

        globalstatus=NEUTRAL;
      }
      else if(globalstatus==C_ADDPOLYCURVE)
      {
        pl.nods[0].d[GX]=0.0;
        pl.nods[0].d[GY]=0.0;
        pl.nods[0].d[GZ]=0.0;

        pl.nvec.dc[GX]=0.0;
        pl.nvec.dc[GY]=0.0;
        pl.nvec.dc[GZ]=1.0;

        pl.a=0.0; pl.b=0.0; pl.c=1.0; pl.d=0.0; /*PLANE:Z=0*/

        vp=&((wsdsp.childs+1)->vparam);
        dmx=(double)( x-(vp->Xo));
        dmy=(double)(-y+(vp->Yo));
        createnodeonplane(*vp,dmx,dmy,pl,&nadd);
		if(nadd.code==0) break;

        if(gpolypoly.npcurve==0)
        {
          gpolypoly.npcurve=1;
          gpolypoly.pcurves
          =(struct polycurve *)malloc(sizeof(struct polycurve));

          cpolycurve=gpolypoly.pcurves+0;
          cpolycurve->loff  =0;
          cpolycurve->ncurve=0;
          cpolycurve->type  =0;
          cpolycurve->curves=NULL;
          cpolycurve->prop.r=  0;
          cpolycurve->prop.g=  0;
          cpolycurve->prop.b=255;
        }

        for(i=0;i<addpolycurve.ncurve;i++)
        {
          (addpolycurve.curves+i)->dots[0]->d[GX]+=nadd.d[GX];
          (addpolycurve.curves+i)->dots[0]->d[GY]+=nadd.d[GY];
          (addpolycurve.curves+i)->dots[1]->d[GX]+=nadd.d[GX];
          (addpolycurve.curves+i)->dots[1]->d[GY]+=nadd.d[GY];

          if((addpolycurve.curves+i)->type==CTYPE_CIRCLE)
          {
            (addpolycurve.curves+i)->center->d[GX]+=nadd.d[GX];
            (addpolycurve.curves+i)->center->d[GY]+=nadd.d[GY];
          }
        }

        addcurvestopolycurve(cpolycurve,&addpolycurve);

        globalstatus=NEUTRAL;
      }
      else if(globalstatus==C_SECTIONDROP)
      {
        copypolypolycurve(&gpolypoly,&(gsect->ppc));
        cpolycurve=gpolypoly.pcurves+gpolypoly.npcurve-1;

        hcursor = LoadCursor(hInstGlobal,"CANCURSORW");
		SetClassLong((wsdsp.childs+1)->hwnd,
                     GCL_HCURSOR,(LONG)hcursor);
        SetClassLong((wsect.childs+2+(gsect->loff))->hwnd,
                     GCL_HCURSOR,(LONG)hcursor);

        globalstatus=NEUTRAL;
      }
      else if(globalstatus==C_PROPERTYDROP)
      {
        pl.nods[0]=setcoord(0.0,0.0,0.0);

        pl.nvec.dc[GX]=0.0;
        pl.nvec.dc[GY]=0.0;
        pl.nvec.dc[GZ]=1.0;

        pl.a=0.0; pl.b=0.0; pl.c=1.0; pl.d=0.0; /*PLANE:Z=0*/

        vp=&((wsdsp.childs+1)->vparam);
        dmx=(double)( x-(vp->Xo));
        dmy=(double)(-y+(vp->Yo));
        createnodeonplane(*vp,dmx,dmy,pl,&nadd);

        if(nadd.code==0)
        {
          MessageBox(NULL,"Failed.","Property",MB_OK);
        }
        else
        {
          for(i=0;i<gpolypoly.npcurve;i++)
          {
            if(dotinpolycurve(gpolypoly.pcurves+i,
                              nadd.d[GX],nadd.d[GY])==1)
            {
              (gpolypoly.pcurves+i)->prop=*gprop;

              (gpolypoly.pcurves+i)->prop.name
              =(char *)realloc((gpolypoly.pcurves+i)->prop.name,
                               (strlen(gprop->name)+1)
                               *sizeof(char));
              strcpy((gpolypoly.pcurves+i)->prop.name,gprop->name);
              break;
            }
		  }
        }

        hcursor = LoadCursor(hInstGlobal,"CANCURSORW");
        SetClassLong((wsdsp.childs+1)->hwnd,
                     GCL_HCURSOR,(LONG)hcursor);
        SetClassLong((wsect.childs+2+(gsect->loff))->hwnd,
                     GCL_HCURSOR,(LONG)hcursor);

        globalstatus=NEUTRAL;
      }
      else if((wsdsp.childs+1)->lstatus==ROTATE)
      {
        initx=x;
        inity=y;

        initv=(wsdsp.childs+1)->vparam;

        globalstatus=ROTATE;
      }
      else if((wsdsp.childs+1)->lstatus==MOVE)
      {
        initx=x;
        inity=y;

        initv=(wsdsp.childs+1)->vparam;

        globalstatus=MOVE;
      }

      if(globalstatus==NEUTRAL)
      {
        clearwindow(*(wsdsp.childs+1));
        drawglobalpolycurve((wsdsp.childs+1)->hdcC,
                            (wsdsp.childs+1)->vparam,
                            &gpolypoly);
        overlayhdc(*(wsdsp.childs+1),SRCPAINT);
      }
      break;

    case WM_RBUTTONDOWN:
      /*globalstatus=NEUTRAL;*/

      point.x = LOWORD(lParam);
      point.y = HIWORD(lParam);
      initx = point.x;
      inity = point.y;

      popupmenudraw(hwnd,point);

      /*clearwindow(*(wsdsp.childs+1));
      drawglobalpolycurve((wsdsp.childs+1)->hdcC,
                          (wsdsp.childs+1)->vparam,
                          &gpolycurve);
      overlayhdc(*(wsdsp.childs+1),SRCPAINT);*/
      break;

    case WM_MOUSEMOVE:
      if(globalstatus==C_SECTIONLINE ||
         globalstatus==C_SECTIONCIRCLE) /*CREATING CURVE.*/
      {
        if(icount<1) break;

        x = LOWORD(lParam);
        y = HIWORD(lParam);

        pl.nods[0].d[GX]=0.0;
        pl.nods[0].d[GY]=0.0;
        pl.nods[0].d[GZ]=0.0;

        pl.nvec.dc[GX]=0.0;
        pl.nvec.dc[GY]=0.0;
        pl.nvec.dc[GZ]=1.0;

        pl.a=0.0; pl.b=0.0; pl.c=1.0; pl.d=0.0; /*PLANE:Z=0*/

        vp=&((wsdsp.childs+1)->vparam);
        dmx=(double)( x-(vp->Xo));
        dmy=(double)(-y+(vp->Yo));
        createnodeonplane(*vp,dmx,dmy,pl,&nadd);
        if(nadd.code==0) break;

        ncurve=cpolycurve->ncurve;
        if(globalstatus==C_SECTIONLINE)
		{
          *((cpolycurve->curves+ncurve-1)->dots[1])=nadd;
        }
        else if(globalstatus==C_SECTIONCIRCLE)
        {
          *((cpolycurve->curves+ncurve-1)->dots[1])=nadd;
        }
        else break;

        chasecurve(cpolycurve->curves+(ncurve-1),
                   (wsdsp.childs+1)->vparam,
                   *(wsdsp.childs+1));
        break;
      }
      else if(globalstatus==C_ADDPOLYCURVE)
      {
        x = LOWORD(lParam);
        y = HIWORD(lParam);

        pl.nods[0].d[GX]=0.0;
        pl.nods[0].d[GY]=0.0;
        pl.nods[0].d[GZ]=0.0;

        pl.nvec.dc[GX]=0.0;
        pl.nvec.dc[GY]=0.0;
        pl.nvec.dc[GZ]=1.0;

        pl.a=0.0; pl.b=0.0; pl.c=1.0; pl.d=0.0; /*PLANE:Z=0*/

        vp=&((wsdsp.childs+1)->vparam);
        dmx=(double)( x-(vp->Xo));
        dmy=(double)(-y+(vp->Yo));
        createnodeonplane(*vp,dmx,dmy,pl,&nadd);
        if(nadd.code==0) break;

        for(i=0;i<addpolycurve.ncurve;i++)
        {
          (addpolycurve.curves+i)->dots[0]->d[GX]+=nadd.d[GX];
          (addpolycurve.curves+i)->dots[0]->d[GY]+=nadd.d[GY];
          (addpolycurve.curves+i)->dots[1]->d[GX]+=nadd.d[GX];
          (addpolycurve.curves+i)->dots[1]->d[GY]+=nadd.d[GY];

		  if((addpolycurve.curves+i)->type==CTYPE_CIRCLE)
          {
            (addpolycurve.curves+i)->center->d[GX]+=nadd.d[GX];
            (addpolycurve.curves+i)->center->d[GY]+=nadd.d[GY];
          }
        }

        chasepolycurve(&addpolycurve,
                       (wsdsp.childs+1)->vparam,
                       *(wsdsp.childs+1));

        for(i=0;i<addpolycurve.ncurve;i++)
        {
          (addpolycurve.curves+i)->dots[0]->d[GX]-=nadd.d[GX];
          (addpolycurve.curves+i)->dots[0]->d[GY]-=nadd.d[GY];
          (addpolycurve.curves+i)->dots[1]->d[GX]-=nadd.d[GX];
          (addpolycurve.curves+i)->dots[1]->d[GY]-=nadd.d[GY];

          if((addpolycurve.curves+i)->type==CTYPE_CIRCLE)
          {
            (addpolycurve.curves+i)->center->d[GX]-=nadd.d[GX];
            (addpolycurve.curves+i)->center->d[GY]-=nadd.d[GY];
          }
        }
        break;
      }
      else if(globalstatus==ROTATE &&
              (wsdsp.childs+1)->lstatus==ROTATE &&
              (wParam == MK_LBUTTON)) /*ROTATE POLYCURVE.*/
      {
        point.x = LOWORD(lParam);
        point.y = HIWORD(lParam);

        x=initx-point.x;
        y=inity-point.y;

        if((initv.Xo<=initx && initv.Yo<=inity && x< y) ||
           (initv.Xo>=initx && initv.Yo<=inity && x<-y) ||
           (initv.Xo>=initx && initv.Yo>=inity && x> y) ||
           (initv.Xo<=initx && initv.Yo>=inity && x>-y))
        {
          fact=0.005;
		}
        else fact=-0.005;

        initx=point.x;
        inity=point.y;

        fact*=sqrt((double)(x*x+y*y));
        rotatepolycurve(&gpolypoly,fact);

        /*clearwindow(*(wsdsp.childs+1));
        drawglobalpolycurve((wsdsp.childs+1)->hdcC,
                            (wsdsp.childs+1)->vparam,
                            &gpolycurve);
        overlayhdc(*(wsdsp.childs+1),SRCPAINT);*/
        chasepolypolycurve(&gpolypoly,
                           (wsdsp.childs+1)->vparam,
                           *(wsdsp.childs+1));
      }
      else if(globalstatus==MOVE &&
              (wsdsp.childs+1)->lstatus==MOVE &&
              (wParam==MK_LBUTTON)) /*MOVE POLYCURVE.*/
      {
        point.x = LOWORD(lParam);
        point.y = HIWORD(lParam);

        x=initx-point.x;
        y=inity-point.y;

        initx=point.x;
        inity=point.y;

        fact=1.0;

        for(ii=0;ii<gpolypoly.npcurve;ii++)
        {
          for(i=0;i<(gpolypoly.pcurves+ii)->ncurve;i++)
          {
            ((gpolypoly.pcurves+ii)->curves+i)->dots[0]->d[GX]
            -=fact*(double)x;
            ((gpolypoly.pcurves+ii)->curves+i)->dots[0]->d[GY]
            +=fact*(double)y;
            ((gpolypoly.pcurves+ii)->curves+i)->dots[1]->d[GX]
			-=fact*(double)x;
            ((gpolypoly.pcurves+ii)->curves+i)->dots[1]->d[GY]
            +=fact*(double)y;

            if(((gpolypoly.pcurves+ii)->curves+i)->type
               ==CTYPE_CIRCLE)
            {
              ((gpolypoly.pcurves+ii)->curves+i)->center->d[GX]
              -=fact*(double)x;
              ((gpolypoly.pcurves+ii)->curves+i)->center->d[GY]
              +=fact*(double)y;
            }
          }
        }

        /*clearwindow(*(wsdsp.childs+1));
        drawglobalpolycurve((wsdsp.childs+1)->hdcC,
                            (wsdsp.childs+1)->vparam,
                            &gpolycurve);
        overlayhdc(*(wsdsp.childs+1),SRCPAINT);*/
        chasepolypolycurve(&gpolypoly,
                           (wsdsp.childs+1)->vparam,
                           *(wsdsp.childs+1));
      }
      else if((globalstatus==ROTATE &&
               (wsdsp.childs+1)->lstatus==ROTATE) ||
              (globalstatus==MOVE &&
               (wsdsp.childs+1)->lstatus==MOVE))     /*END OF MOVE.*/
      {
        globalstatus=prestatus; /*PREVIOUS STATUS.*/

        clearwindow(*(wsdsp.childs+1));
        drawglobalpolycurve((wsdsp.childs+1)->hdcC,
                            (wsdsp.childs+1)->vparam,
                            &gpolypoly);
        overlayhdc(*(wsdsp.childs+1),SRCPAINT);

        /*globalstatus=NEUTRAL;*/ /*ROTATION END*/
      }
      else if(globalstatus==ROTATE || globalstatus==MOVE)
      {
        /*globalstatus=NEUTRAL;*/ /*ROTATION END*/
		globalstatus=prestatus; /*ROTATION END*/
      }
      break;

    case WM_COMMAND:
      switch(LOWORD(wParam))
      {
        case IDM_FITSHEET: /*NEVER USE FITPARENT.*/
          hoya=GetParent(GetParent(hwnd)); /*PARENT SHEET.*/

          getclientsize(wsfrm.hwnd,&maxX,&maxY);

          sw=0; sh=0;
          pw=0; ph=0;
          if(hfitfrom==wsfrm.hwnd)
          {
            getwindowsize(wprop.hwnd,&pw,&ph);
            getwindowsize(wsect.hwnd,&sw,&sh);
            if(pw>sw) mw=pw;
            else      mw=sw;
          }
          else if(hfitfrom==(wprop.childs+1)->hwnd)
          {
            getwindowsize(wprop.hwnd,&pw,&ph);
            mw=pw;
          }
          else if(hfitfrom==(wsect.childs+1)->hwnd)
          {
            getwindowsize(wsect.hwnd,&sw,&sh);
            mw=sw;
          }

          if(maxX-mw<50) maxX=mw+50;

          MoveWindow(hoya,mw,0,(maxX-mw),maxY,TRUE);
          break;

        case IDM_FITPARENT: /*MESSAGE FROM PARENT SHEET.*/
          wparam = MAKEWPARAM((WORD)IDM_FITSHEET,(WORD)0);
          if(wsect.hwnd!=NULL)
          {
            hfitfrom=(wsdsp.childs+1)->hwnd;
			SendMessage((wsect.childs+1)->hwnd,WM_COMMAND,wparam,0);
          }
          if(wprop.hwnd!=NULL)
          {
            hfitfrom=(wsdsp.childs+1)->hwnd;
            SendMessage((wprop.childs+1)->hwnd,WM_COMMAND,wparam,0);
          }
          break;

        case IDM_POPSECTIONLINE:
          globalstatus=C_SECTIONLINE;
          break;

        case IDM_POPCIRCLECENTER:
          globalstatus=C_SECTIONCIRCLE;
          break;
        /*
		case IDM_POPCIRCLERADIUS:
          break;
        case IDM_POPCIRCLEDOT:
          break;
        case IDM_POPCIRCLETANGENT:
          break;
        */
        case IDM_POPCIRCLEREVERSE:
          if(globalstatus==C_SECTIONCIRCLE)
          {
            ncurve=cpolycurve->ncurve;
            if((cpolycurve->curves+ncurve-1)->hugo==1)
            {
              (cpolycurve->curves+ncurve-1)->hugo=-1;
              if(icount==2)
              {
                (cpolycurve->curves+ncurve-1)->angle[0]+=2.0*PI;
              }
            }
            else if((cpolycurve->curves+ncurve-1)->hugo==-1)
            {
              (cpolycurve->curves+ncurve-1)->hugo=1;
              if(icount==2)
              {
                (cpolycurve->curves+ncurve-1)->angle[0]-=2.0*PI;
			  }
            }
          }
          break;

        case IDM_POPREGISTSECTION:
          point.x=initx;
          point.y=inity;

          hpopdlg=CreateDialog(hInstGlobal,
                               "HOGDLGSECTREGIST",
                               NULL,
                               (DLGPROC)DialogProcSectRegist);
          getwindowsize(hpopdlg,&mw,&mh);
          ClientToScreen((wsdsp.childs+1)->hwnd,&point);
          MoveWindow(hpopdlg,point.x,point.y,mw,mh,TRUE);
          ShowWindow(hpopdlg,SW_SHOW);
          break;
        case IDM_POPREGISTRETURN:
          /*SECTION CODE*/
          GetDlgItemText(hpopdlg,IDSR_POPCODE,str,20);
          code=strtol(str,NULL,10);

          gsect=currentsects;
          while(code != gsect->code) gsect++;

          copypolypolycurve(&(gsect->ppc),&gpolypoly);

          GetDlgItemText(hpopdlg,IDSR_POPCAPTION,
                         gsect->ppc.name,256);

          gsect->ppc.hico
          =createsectionicon((wsect.childs+1)->hdcC,32,32,
                             (wsdsp.childs+1)->vparam,
                             &(gsect->ppc));

          SendMessage((wsect.childs+2+(gsect->loff))->hwnd,
                      WM_PAINT,0,0);
          break;

        case IDM_POPSECTIONDIALOG:
          point.x=initx;
		  point.y=inity;

          hpopdlg=CreateDialog(hInstGlobal,
                               "HOGDLGKATAKOU",
                               NULL,
                               (DLGPROC)DialogProcKatakou);
          getwindowsize(hpopdlg,&mw,&mh);
          ClientToScreen((wsdsp.childs+1)->hwnd,&point);
          MoveWindow(hpopdlg,point.x,point.y,mw,mh,TRUE);
          ShowWindow(hpopdlg,SW_SHOW);
          break;
        case IDM_POPSECTIONRETURN:
          if(addpolycurve.ncurve>0)
          freepolycurve(&addpolycurve);

		  addpolycurve.prop.r=100;
          addpolycurve.prop.g=100;
          addpolycurve.prop.b=100;

          getkatakouparam(hpopdlg,
                          &H,&B,&tf1,&tf2,&tw,
                          &ri1,&ri2,&ri3,&ri4,&ro1,&ro2,&ro3,&ro4);

          setpolycurveaskatakou(&addpolycurve,
                                addpolycurve.type,
                                H,B,tf1,tf2,tw,
                                ri1,ri2,ri3,ri4,ro1,ro2,ro3,ro4);

          /*clearwindow(*(wsdsp.childs+1));
          drawglobalpolycurve((wsdsp.childs+1)->hdcC,
                              (wsdsp.childs+1)->vparam,
                              &gpolycurve);
          overlayhdc(*(wsdsp.childs+1),SRCPAINT);*/

          globalstatus=C_ADDPOLYCURVE;
          break;

        case IDM_POPSECTIONPROPERTY:
          point.x=initx;
          point.y=inity;

          pc=pickpolycurve(&gpolypoly,
						   (wsdsp.childs+1)->vparam,initx,inity);
          if(pc!=NULL) cpolycurve=pc;
          else
          {
            MessageBox(NULL,"Nothing.","Pick",MB_OK);
            break;
          }

          gprop=&(cpolycurve->prop);

          ipopflag=0;

          hpopdlg=CreateDialog(hInstGlobal,
                               "HOGDLGPROPERTY",
                               NULL,
                               (DLGPROC)DialogProcProperty);

          getwindowsize(hpopdlg,&mw,&mh);
          ClientToScreen((wsdsp.childs+1)->hwnd,&point);
          MoveWindow(hpopdlg,point.x,point.y,mw,mh,TRUE);
          ShowWindow(hpopdlg,SW_SHOW);
          break;
        case IDM_POPPROPERTYRETURN:
          SendMessage((wsdsp.childs+1)->hwnd,WM_PAINT,0,0);
          ipopflag=0;
          break;

        case IDM_POPSECTIONMOVE:
          (wsdsp.childs+1)->lstatus=MOVE;
          break;
        case IDM_POPSECTIONROTATE:
          (wsdsp.childs+1)->lstatus=ROTATE;
          break;

        case IDM_POPPOLYADD:
          gpolypoly.npcurve++;
          gpolypoly.pcurves
          =(struct polycurve *)
           realloc(gpolypoly.pcurves,
                   gpolypoly.npcurve*sizeof(struct polycurve));

          cpolycurve=gpolypoly.pcurves+gpolypoly.npcurve-1;
		  cpolycurve->loff  =0;
          cpolycurve->ncurve=0;
          cpolycurve->type  =0;
          cpolycurve->curves=NULL;
          cpolycurve->prop.r=  0;
          cpolycurve->prop.g=  0;
		  cpolycurve->prop.b=255;
          break;
        case IDM_POPPOLYPICK:
          globalstatus=C_PICKPOLYCURVE;
          break;

        case IDM_POPCHOOSENODE:
          (wsdsp.childs+1)->lstatus=SELECTNODE;
          hcursor = LoadCursor(hInstGlobal,"CANBOXW");
          SetClassLong((wsdsp.childs+1)->hwnd,
                       GCL_HCURSOR,(LONG)hcursor);
          break;

        case IDMS_VIEW:
          clearwindow(*(wsdsp.childs+1));
          drawglobalpolycurve((wsdsp.childs+1)->hdcC,
                              (wsdsp.childs+1)->vparam,
                              &gpolypoly);
          overlayhdc(*(wsdsp.childs+1),SRCPAINT);
          break;

        default: /*OTHERS.*/
          return DefWindowProc(hwnd,message,wParam,lParam);
      }
      break;

    default:
      return DefWindowProc(hwnd,message,wParam,lParam);
  }
  return 0;
}/*WindowProcedureSview*/

LRESULT CALLBACK WindowProcedureMesg(HWND hwnd,
                                     UINT message,
                                     WPARAM wParam,
                                     LPARAM lParam)
/*WINDOW FOR ERROR MESSAGE.*/
{
  char str[256];
  long int x,y;

  switch(message)
  {
    case WM_PAINT:
    case WM_SIZE:
      DefWindowProc(hwnd,message,wParam,lParam);
      overlayhdc(*(wmesg.childs+1),SRCPAINT);
      break;

    case WM_LBUTTONDOWN:
      x = LOWORD(lParam);
      y = HIWORD(lParam);
      sprintf(str,"POSITION:%d %d",x,y);
	  SetTextColor((wmesg.childs+1)->hdcC,RGB(0,255,255));
	  TextOut((wmesg.childs+1)->hdcC,x,y,str,strlen(str));
      SendMessage(hwnd,WM_PAINT,0,0);
      break;

    default:
      return DefWindowProc(hwnd,message,wParam,lParam);
  }
  return 0;
}/*WindowProcedureMesg*/

LRESULT CALLBACK WindowProcedureMenu(HWND hwnd,
                                     UINT message,
                                     WPARAM wParam,
                                     LPARAM lParam)
/*WINDOW FOR OPTION MENU DIALOG.*/
{
  HWND hoya;
  long int maxX,maxY,mw,mh;

  switch(message)
  {
    case WM_COMMAND:
      switch(LOWORD(wParam))
      {
		case IDM_FITSHEET:
          hoya=GetParent(GetParent(hwnd)); /*PARENT SHEET.*/

          getclientsize(wmain.hwnd,&maxX,&maxY);
          getwindowsize(wmenu.hwnd,&mw,&mh);

          MoveWindow(hoya,0,0,mw,maxY,TRUE); /*KEEP WIDTH.*/
          break;
      }
      break;

    default:
      return DefWindowProc(hwnd,message,wParam,lParam);
  }
  return 0;
}/*WindowProcedureMenu*/

static BOOL CALLBACK DialogProcMenu1(HWND hdwnd,
                                     UINT message,
                                     WPARAM wParam,
                                     LPARAM lParam)
/*OPTION MENU DIALOG BOX 1.*/
{
  HWND hitem;
  HDC hdc,hdcC;
  HBITMAP hbit,pbit;
  HCURSOR hcursor;
  /*HPEN hpen,ppen;*/
  /*HFONT hfont;*/
  WPARAM wparam;
  BITMAP bmp;
  /*LPDRAWITEMSTRUCT lpdis;*/
  char str[256],dir[]=DIRECTORY;
  int i,id,it;
  size_t n;
  /*long int maxX,maxY,mw,mh;*/
  FILE *fin;
  struct globvisible gv;
  struct snode *sn;
  int laps;
  double dsafety;

  static OPENFILENAME ofn = {0};
  static TCHAR strFile[256] ,strFileTitle[256], strCustom[256] = TEXT("Before files\0*.*\0\0");

  memset(&ofn,0,sizeof (OPENFILENAME));
  ofn.lStructSize = sizeof (OPENFILENAME);
		ofn.hwndOwner = hdwnd;
		ofn.lpstrFilter = 	TEXT("inp files {*.inp}\0*.inp\0");
		ofn.lpstrCustomFilter = strCustom;
		ofn.nMaxCustFilter = 256;
		ofn.nFilterIndex = 0;
		ofn.lpstrFile = strFile;
		ofn.lpstrFileTitle = strFileTitle;
		ofn.nMaxFile = 256;
        ofn.nMaxFileTitle = 256;
		ofn.Flags = OFN_FILEMUSTEXIST;


  switch(message)
  {
	case WM_INITDIALOG:
	  SetDlgItemText(hdwnd,IDD_TXTDRAWINGS,"Drawings");
      SetDlgItemText(hdwnd,IDD_TXTERROR,"Errors");
      SetDlgItemText(hdwnd,IDD_TXTSURFACE,"Surface");

      SetDlgItemText(hdwnd,IDD_TXTSECTIONLIST,"SectionList");
      SetDlgItemText(hdwnd,IDD_TXTWEIGHTLIST,"WeightList");
      SetDlgItemText(hdwnd,IDD_TXTCMQDETAIL,"CmqDetail");
      SetDlgItemText(hdwnd,IDD_TXTHORIZONTAL,"HorizontalLoad");

      SetDlgItemText(hdwnd,IDF_TXTARCLM,"Arclm");
      SetDlgItemText(hdwnd,IDF_TXTORGAN,"Organ");
      SetDlgItemText(hdwnd,IDF_TXTHOGAN,"Hogan");
      SetDlgItemText(hdwnd,IDF_TXTSRCAN,"Srcan");
	  SetDlgItemText(hdwnd,IDF_TXTFRAME,"Frame");
	  SetDlgItemText(hdwnd,IDF_TXTBQUAD,"Bquad");

      SetDlgItemText(hdwnd,IDV_TXTAXONO,"Axonometric");
      SetDlgItemText(hdwnd,IDV_TXTPERS, "Perspective");

      SetDlgItemText(hdwnd,IDM_TXTPRINTRANGE,"PrintRange");    //by mihara
      SetDlgItemText(hdwnd,IDP_TXTA4TATE,"A4Tate");            //by mihara
      SetDlgItemText(hdwnd,IDP_TXTA4YOKO,"A4Yoko");            //by mihara
	  SetDlgItemText(hdwnd,IDP_TXTA3TATE,"A3Tate");            //by mihara
      SetDlgItemText(hdwnd,IDP_TXTA3YOKO,"A3Yoko");            //by mihara

      if(wdraw.nchilds>=2 && (wdraw.childs+1)!=NULL)
      {
        setviewparam(hdwnd,(wdraw.childs+1)->vparam);

        SetDlgItemText(hdwnd,ID_INPUTFILE,
                       (wdraw.childs+1)->inpfile);
        SetDlgItemText(hdwnd,ID_INPUTFILEZ,
                       (wdraw.childs+1)->inpfilez);
        SetDlgItemText(hdwnd,ID_INPUTFILEX,
                       (wdraw.childs+1)->inpfilex);
        SetDlgItemText(hdwnd,ID_INPUTFILEY,
                       (wdraw.childs+1)->inpfiley);

        SetDlgItemText(hdwnd,ID_OUTPUTFILE,
                       (wdraw.childs+1)->otpfile);
        SetDlgItemText(hdwnd,ID_OUTPUTFILEZ,
                       (wdraw.childs+1)->otpfilez);
        SetDlgItemText(hdwnd,ID_OUTPUTFILEX,
                       (wdraw.childs+1)->otpfilex);
        SetDlgItemText(hdwnd,ID_OUTPUTFILEY,
                       (wdraw.childs+1)->otpfiley);

        SetDlgItemText(hdwnd,ID_SECTIONFILE,
                       (wdraw.childs+1)->sctfile);
      }
      else
      {
        setviewparam(hdwnd,vpdefault);
      }

	  SetDlgItemText(hdwnd,ID_LAPS,"100");
	  SetDlgItemText(hdwnd,ID_SAFETY,"0.1");

      break;

    case WM_CTLCOLORSTATIC: /*STATIC TEXT MENU.*/
      DefWindowProc(hdwnd,message,wParam,lParam);

      hdc=(HDC)wParam;
	  hitem=(HWND)lParam;

      id=GetDlgCtrlID(hitem);
      gv=(wmenu.childs+2)->vparam.vflag;

	  if((id==IDD_TXTDRAWINGS && gv.mv.draw!=1)    ||
         (id==IDD_TXTERROR    && gv.mv.error!=1)   ||
         (id==IDD_TXTSURFACE  && gv.mv.surface!=1) ||
         (id==IDD_TXTSECTIONLIST && gv.mv.sectlist!=1) ||
         (id==IDD_TXTWEIGHTLIST  && gv.mv.weight!=1)   ||
         (id==IDD_TXTCMQDETAIL   && gv.mv.cmq!=1)      ||
         (id==IDD_TXTHORIZONTAL  && gv.mv.horizon!=1)  ||
         (id==IDF_TXTARCLM && gv.mv.ftype!=F_ARCLM) ||
		 (id==IDF_TXTORGAN && gv.mv.ftype!=F_ORGAN) ||
		 (id==IDF_TXTHOGAN && gv.mv.ftype!=F_HOGAN) ||
		 (id==IDF_TXTSRCAN && gv.mv.ftype!=F_SRCAN) ||
		 (id==IDF_TXTFRAME && gv.mv.ftype!=F_FRAME) ||
		 (id==IDF_TXTBQUAD && gv.mv.ftype!=F_BQUAD))
	  {
		SetTextColor(hdc,RGB(255,255,255));
		//SetTextColor(hdc,RGB(200,200,200));  //ujioka
	  }

      if((id==IDV_TXTAXONO) &&
		 (wdraw.hwnd==NULL ||
          (wdraw.childs+1)->vparam.type!=AXONOMETRIC))
      {
		SetTextColor(hdc,RGB(255,255,255));
      }
      if((id==IDV_TXTPERS) &&
         (wdraw.hwnd==NULL ||
          (wdraw.childs+1)->vparam.type!=PERSPECTIVE))
	  {
		SetTextColor(hdc,RGB(255,255,255));
      }

//////////////////////////////////////////////////by MIHARA for Print Range Line

      if(wdraw.hwnd==NULL ||
         (((id==IDM_TXTPRINTRANGE) && (wdraw.childs+1)->vparam.vflag.pv.printrange!=1)||
          ((id==IDP_TXTA4TATE) && (wdraw.childs+1)->vparam.vflag.pv.a4tate!=1)||
          ((id==IDP_TXTA4YOKO) && (wdraw.childs+1)->vparam.vflag.pv.a4yoko!=1)||
		  ((id==IDP_TXTA3TATE) && (wdraw.childs+1)->vparam.vflag.pv.a3tate!=1)||
          ((id==IDP_TXTA3YOKO) && (wdraw.childs+1)->vparam.vflag.pv.a3yoko!=1)
         )
        )
      {
//MessageBox(NULL,"Pass 4","A4",MB_OK);
		SetTextColor(hdc,RGB(255,255,255));
		//SetTextColor(hdc,RGB(200,200,200));       //ujioka
	  }

////////////////////////////////////////////////////////////////////////////////

      /*if(wdraw.hwnd!=NULL) gv=(wdraw.childs+1)->vparam.vflag;

      if((id==IDV_TXTGLOBALAXIS) &&
         (wdraw.hwnd==NULL || gv.axis!=1))
      {
        SetTextColor(hdc,RGB(255,255,255));
      }
      if((id==IDV_TXTNODECODE) &&
         (wdraw.hwnd==NULL || gv.nv.code!=1))
      {
        SetTextColor(hdc,RGB(255,255,255));
      }*/
      return (BOOL)(HBRUSH)GetStockObject(LTGRAY_BRUSH);
      //return (BOOL)(HBRUSH)GetStockObject(NULL_BRUSH);

    case WM_PAINT: /*INSERT ILLUSTRATION.*/
      /*DefWindowProc(hwnd,message,wParam,lParam);*/
/*
	  hdc=GetDC(hdwnd);
      hdcC = CreateCompatibleDC(hdc);

      hbit = LoadBitmap(hInstGlobal,"BITMAPVIEW");
      GetObject(hbit,sizeof(BITMAP),(LPSTR)&bmp);
	  pbit = (HBITMAP)SelectObject(hdcC,hbit);

	  BitBlt(hdc,1,(int)((DPT_VIEW+26)*DIALOGFACTOR),
				 bmp.bmWidth,bmp.bmHeight,
				 hdcC,0,0,SRCCOPY);

	  BitBlt(hdc,1,(int)((DPT_VIEW+26)*2.75),
				 bmp.bmWidth,bmp.bmHeight,
				 hdcC,0,0,SRCCOPY);

	  SelectObject(hdcC,pbit);
	  DeleteObject(hbit);
      ReleaseDC(hdwnd,hdcC);
	  ReleaseDC(hdwnd,hdc);
*/
	  break;

    case WM_COMMAND:
      switch(LOWORD(wParam))
      {
        case IDF_ARCLM:
          if((wmenu.childs+2)->vparam.vflag.mv.ftype!=F_ARCLM)
          {
            (wmenu.childs+2)->vparam.vflag.mv.ftype=F_ARCLM;
            SendMessage(hdwnd,WM_INITDIALOG,0,0);
          }
          break;
        case IDF_ORGAN:
          if((wmenu.childs+2)->vparam.vflag.mv.ftype!=F_ORGAN)
          {
            (wmenu.childs+2)->vparam.vflag.mv.ftype=F_ORGAN;
            SendMessage(hdwnd,WM_INITDIALOG,0,0);
          }
          break;
        case IDF_HOGAN:
          if((wmenu.childs+2)->vparam.vflag.mv.ftype!=F_HOGAN)
          {
            (wmenu.childs+2)->vparam.vflag.mv.ftype=F_HOGAN;
            sprintf((wdraw.childs+1)->inpfile,"hogan01.inp");

            SendMessage(hdwnd,WM_INITDIALOG,0,0);
          }
          break;
        case IDF_SRCAN:
          if((wmenu.childs+2)->vparam.vflag.mv.ftype!=F_SRCAN)
          {
            (wmenu.childs+2)->vparam.vflag.mv.ftype=F_SRCAN;
            SendMessage(hdwnd,WM_INITDIALOG,0,0);
          }
		  break;
        case IDF_FRAME:
          if((wmenu.childs+2)->vparam.vflag.mv.ftype!=F_FRAME)
          {
			(wmenu.childs+2)->vparam.vflag.mv.ftype=F_FRAME;
            SendMessage(hdwnd,WM_INITDIALOG,0,0);
          }
          break;
		case IDF_BQUAD:
		  if((wmenu.childs+2)->vparam.vflag.mv.ftype!=F_BQUAD)
		  {
			(wmenu.childs+2)->vparam.vflag.mv.ftype=F_BQUAD;
			SendMessage(hdwnd,WM_INITDIALOG,0,0);
		  }
		  break;

        case ID_INPUTFILE: /*GET INPUT OF KEYBOARD ONE BY ONE.*/
          if(wdraw.nchilds>=2 && (wdraw.childs+1)!=NULL)
          {
            GetDlgItemText(hdwnd,ID_INPUTFILE,
						   (wdraw.childs+1)->inpfile,80);
          }
          break;
        case ID_INPUTFILEZ:
          if(wdraw.nchilds>=2 && (wdraw.childs+1)!=NULL)
          {
            GetDlgItemText(hdwnd,ID_INPUTFILEZ,
                           (wdraw.childs+1)->inpfilez,80);
          }
          break;
        case ID_INPUTFILEX:
          if(wdraw.nchilds>=2 && (wdraw.childs+1)!=NULL)
          {
            GetDlgItemText(hdwnd,ID_INPUTFILEX,
                           (wdraw.childs+1)->inpfilex,80);
          }
          break;
        case ID_INPUTFILEY:
          if(wdraw.nchilds>=2 && (wdraw.childs+1)!=NULL)
          {
            GetDlgItemText(hdwnd,ID_INPUTFILEY,
                           (wdraw.childs+1)->inpfiley,80);
		  }
          break;
        case ID_OUTPUTFILE:
          if(wdraw.nchilds>=2 && (wdraw.childs+1)!=NULL)
          {
            GetDlgItemText(hdwnd,ID_OUTPUTFILE,
                           (wdraw.childs+1)->otpfile,80);
          }
          break;
        case ID_OUTPUTFILEZ:
          if(wdraw.nchilds>=2 && (wdraw.childs+1)!=NULL)
          {
            GetDlgItemText(hdwnd,ID_OUTPUTFILEZ,
                           (wdraw.childs+1)->otpfilez,80);
          }
          break;
        case ID_OUTPUTFILEX:
          if(wdraw.nchilds>=2 && (wdraw.childs+1)!=NULL)
          {
            GetDlgItemText(hdwnd,ID_OUTPUTFILEX,
                           (wdraw.childs+1)->otpfilex,80);
          }
          break;
        case ID_OUTPUTFILEY:
          if(wdraw.nchilds>=2 && (wdraw.childs+1)!=NULL)
          {
            GetDlgItemText(hdwnd,ID_OUTPUTFILEY,
                           (wdraw.childs+1)->otpfiley,80);
          }
          break;
        case ID_SECTIONFILE:
          if(wdraw.nchilds>=2 && (wdraw.childs+1)!=NULL)
          {
            GetDlgItemText(hdwnd,ID_SECTIONFILE,
                           (wdraw.childs+1)->sctfile,80);
          }
          break;

        case IDV_AXONO:
          if(wdraw.nchilds>=2 &&
             (wdraw.childs+1)!=NULL &&
             (wdraw.childs+1)->vparam.type!=AXONOMETRIC)
		  {
            (wdraw.childs+1)->vparam.type=AXONOMETRIC;
            SendMessage(hdwnd,WM_INITDIALOG,0,0);
          }
          break;
        case IDV_PERS:
          if(wdraw.nchilds>=2 &&
             (wdraw.childs+1)!=NULL &&
             (wdraw.childs+1)->vparam.type!=PERSPECTIVE)
          {
            (wdraw.childs+1)->vparam.type=PERSPECTIVE;
            SendMessage(hdwnd,WM_INITDIALOG,0,0);
          }
          break;

        case IDV_GFACTOR:
        case IDV_X:
        case IDV_Y:
        case IDV_Z:
        case IDV_PHI:
        case IDV_THETA:
        case IDV_R:
        case IDV_L:
        case IDV_ORIGINX:
        case IDV_ORIGINY:
		case IDR_XMAX:
        case IDR_XMIN:
        case IDR_YMAX:
        case IDR_YMIN:
        case IDR_ZMAX:
		case IDR_ZMIN:
        case IDV_GAXISLENGTH:
		case IDV_EAXISLENGTH:
		case IDV_DFACTOR:
        case IDV_QFACTOR:
        case IDV_MFACTOR:
		case IDV_GYOPITCH:
        case IDV_HINGESIZE:
          if(wdraw.nchilds>=2 && (wdraw.childs+1)!=NULL)
		  {
            getviewparam((wmenu.childs+2)->hwnd,
						 &((wdraw.childs+1)->vparam));
		  }
          break;

		case IDD_DRAWINGS:
		  wparam = MAKEWPARAM((WORD)IDM_DRAWINGS,(WORD)0);
          SendMessage(wmain.hwnd,WM_COMMAND,wparam,(WPARAM)0);

          GetDlgItemText(hdwnd,ID_INPUTFILE,
                         (wdraw.childs+1)->inpfile,80);
          break;
        case IDD_ERROR:
          wparam = MAKEWPARAM((WORD)IDM_ERROR,(WORD)0);
          SendMessage(wmain.hwnd,WM_COMMAND,wparam,(WPARAM)0);
          break;
        case IDD_SURFACE:
          wparam = MAKEWPARAM((WORD)IDM_SURFACE,(WORD)0);
          SendMessage(wmain.hwnd,WM_COMMAND,wparam,(WPARAM)0);
          break;
        case IDD_SECTIONLIST:
          /*wparam = MAKEWPARAM((WORD)IDM_SECTIONLIST,(WORD)0);*/
          wparam = MAKEWPARAM((WORD)IDM_SECTIONFRAME,(WORD)0);
          SendMessage(wmain.hwnd,WM_COMMAND,wparam,(WPARAM)0);
          break;
        case IDD_WEIGHTLIST:
		  wparam = MAKEWPARAM((WORD)IDM_WEIGHTLIST,(WORD)0);
          SendMessage(wmain.hwnd,WM_COMMAND,wparam,(WPARAM)0);
          break;
        case IDD_CMQDETAIL:
          wparam = MAKEWPARAM((WORD)IDM_CMQDETAIL,(WORD)0);
          SendMessage(wmain.hwnd,WM_COMMAND,wparam,(WPARAM)0);
          break;
        case IDD_HORIZONTAL:
          wparam = MAKEWPARAM((WORD)IDM_HORIZONTAL,(WORD)0);
          SendMessage(wmain.hwnd,WM_COMMAND,wparam,(WPARAM)0);
          break;

		case IDD_ARCLM001:
          /*
          if(MessageBox(NULL,"Arclm001:Begin.","ARCLM001",
             MB_OKCANCEL)==IDCANCEL) break;

          if((wmenu.childs+2)->vparam.vflag.mv.ftype==F_ARCLM)
		  {
            getviewparam((wmenu.childs+2)->hwnd,
                         &((wdraw.childs+1)->vparam));
            clearwindow(*(wdraw.childs+1));
			arclm001(&arc,ID_INPUTFILE,ID_OUTPUTFILE);
          }
          break;
		  */
		  if(MessageBox(NULL,"Arclm001:Begin.","ARCLM001",
			 MB_OKCANCEL)==IDOK)
		  {
          if((wmenu.childs+2)->vparam.vflag.mv.ftype==F_ARCLM)
            {
            getviewparam((wmenu.childs+2)->hwnd,
                         &((wdraw.childs+1)->vparam));
            clearwindow(*(wdraw.childs+1));
			arclm001(&arc,ID_INPUTFILE,ID_OUTPUTFILE);
            }
          }

		  else if(MessageBox(NULL,"Arclm002:Begin","ARCLM002",
			   MB_OKCANCEL)==IDOK)
          {
          if((wmenu.childs+2)->vparam.vflag.mv.ftype==F_ARCLM)
            {
            getviewparam((wmenu.childs+2)->hwnd,
                         &((wdraw.childs+1)->vparam));
            clearwindow(*(wdraw.childs+1));
			arclm002(&arc,ID_INPUTFILE,ID_OUTPUTFILE);
            }
		  }
         break;


		case IDD_ARCLM101:
          if((wmenu.childs+2)->vparam.vflag.mv.ftype!=F_ARCLM)
          {
            MessageBox(NULL,"File Type Error.","ARCLM101",MB_OK);
			break;
		  }
		  if(MessageBox(NULL,"Initial Analysis Arclm001.",
			 "ARCLM101",MB_OKCANCEL)==IDOK)
		  {
			arclm001(&arc,ID_INPUTFILEZ,ID_OUTPUTFILEZ);  //長期の解析
		  }

		  if(MessageBox(NULL,"Arclm101:Begin.","ARCLM101",
			 MB_OKCANCEL)==IDOK)
          {
          	if(MessageBox(NULL,"Reflect Blng002 Result","ARCLM101",MB_OKCANCEL)==IDOK)
            {
            	getviewparam((wmenu.childs+2)->hwnd,&((wdraw.childs+1)->vparam));
		  		clearwindow(*(wdraw.childs+1));
				arclm101_bc(&arc,ID_INPUTFILEX);
		        //arclm101_bc(&arc,ID_INPUTFILE);
            }
            else
			{
   		        getviewparam((wmenu.childs+2)->hwnd,&((wdraw.childs+1)->vparam));
		        clearwindow(*(wdraw.childs+1));
				//arclm101(&arc,ID_INPUTFILEX);
		        arclm101(&arc,ID_INPUTFILE);
            }
          }
		  break;

		case IDD_ARCLMSTATIC:
		   if ((wmenu.childs + 2)->vparam.vflag.mv.ftype == F_ARCLM)
		   {
			   getviewparam((wmenu.childs + 2)->hwnd,&((wdraw.childs + 1)->vparam));
			   clearwindow(*(wdraw.childs + 1));
			   arclmStatic(&arc);
		   }
		   break;


		case IDD_ARCLMDYNAMIC:
		   if ((wmenu.childs + 2)->vparam.vflag.mv.ftype == F_ARCLM)
		   {
			   getviewparam((wmenu.childs + 2)->hwnd,&((wdraw.childs + 1)->vparam));
			   clearwindow(*(wdraw.childs + 1));
			   arclmDynamic(&arc);
		   }
		   break;


		case IDD_ARCLMSTATDYNA:
		  if((wmenu.childs+2)->vparam.vflag.mv.ftype==F_ARCLM)
		  {
			getviewparam((wmenu.childs+2)->hwnd,&((wdraw.childs+1)->vparam));
			clearwindow(*(wdraw.childs+1));
			arclmStatDyna(&arc);
		  }
		  break;

		case IDD_BCLNG001:
		  /*if(MessageBox(NULL,"Bclng001:Begin.","BCLNG001",
			 MB_OKCANCEL)==IDCANCEL) break;*/

		  if((wmenu.childs+2)->vparam.vflag.mv.ftype==F_ARCLM)
		  {
			if(arc.elems==NULL) break;
			getviewparam((wmenu.childs+2)->hwnd,
						 &((wdraw.childs+1)->vparam));
			/*clearwindow(*(wdraw.childs+1));*/
		  /*
			if(MessageBox(NULL,"Bclng002:Buckling Condensation Begin.","BCLNG002",
			   MB_OKCANCEL)==IDOK) bclng002(&arc);
			else if(MessageBox(NULL,"Bclng001:Buckling Analysis Begin.","BCLNG001",
			   MB_OKCANCEL)==IDOK) bclng001(&arc);
		  }
		  */
			if(MessageBox(NULL,"Bclng001:Buckling Analysis Begin.","BCLNG001",
			   MB_OKCANCEL)==IDOK) bclng001(&arc);
		  	else if(MessageBox(NULL,"Bclng002:Buckling Condensation Begin.","BCLNG002",
		  	   MB_OKCANCEL)==IDOK) bclng002(&arc);
		  }
		  break;

		case IDD_BCLNG002:
		  if((wmenu.childs+2)->vparam.vflag.mv.ftype==F_ARCLM)
		  {
			if(arc.elems==NULL) break;
			getviewparam((wmenu.childs+2)->hwnd,
						 &((wdraw.childs+1)->vparam));

			if(MessageBox(NULL,"Bclng002:Buckling Condensation Begin.","BCLNG002",
			   MB_OKCANCEL)==IDOK) bclng002(&arc);

            double *ncr;
            ncr=mallocdoublevector(arc.nelem);
            definencr(&arc,&*ncr);

		  }
		  break;

		case IDD_BCLNG011:
		  if((wmenu.childs+2)->vparam.vflag.mv.ftype==F_ARCLM)
		  {

			struct arclmframe arcII  =  {0,0,
										"\0",
										0,0,0,0,0,0,0,
										NULL,NULL,
										NULL,NULL,NULL,NULL,NULL,
										NULL,NULL,NULL,NULL,
										NULL,NULL,
										NULL,NULL,NULL,
										NULL,NULL,
										0,NULL,NULL,NULL}; /*ARCLM FRAME.*/
		  #if 1
          /*STRESS OF LOAD*/
            arclm001(&arc,ID_INPUTFILE,ID_OUTPUTFILE);

		  /*STRESS OF LONG LOAD*/
            fin=fgetstofopen(dir,"r",ID_INPUTFILEZ);    /*OPEN FILE.*/
            if(fin==NULL) break;

            inputinit(fin,&(arcII.nnode),&(arcII.nelem),&(arcII.nsect));
            arcII.sects=(struct osect *)
                      malloc(arcII.nsect*sizeof(struct osect));
            if(arcII.sects==NULL) break;
            arcII.nodes=(struct onode *)
                      malloc(arcII.nnode*sizeof(struct onode));
            if(arcII.nodes==NULL) break;
            arcII.ninit=(struct onode *)
                      malloc(arcII.nnode*sizeof(struct onode));
            if(arcII.ninit==NULL) break;
            arcx.elems=(struct owire *)
                      malloc(arcII.nelem*sizeof(struct owire));
            if(arcII.elems==NULL) break;
            arcx.confs=(struct oconf *)
                      malloc(6*arcII.nnode*sizeof(struct oconf));
            if(arcII.confs==NULL) break;

			inputtexttomemory(fin,&arcII);
            fclose(fin);

			arclm001(&arcII,ID_INPUTFILEZ,ID_OUTPUTFILEZ);
          #endif

            fout=fgetstofopen("\0","r",ID_OUTPUTFILE);  /*OPEN FILE.*/
            if(fout==NULL) break;
            frameoutputtomemory(fout,&arc);
            fclose(fout);

            fout=fgetstofopen("\0","r",ID_OUTPUTFILEZ);  /*OPEN FILE.*/
            if(fout==NULL) break;
            frameoutputtomemory(fout,&arcII);
            fclose(fout);

			if(MessageBox(NULL,"Bclng011:BCLNG011 Begin.","BCLNG011",
			   MB_OKCANCEL)==IDOK) bclng011(&arc,&arcII);
		  }
		  break;

		case IDD_BCLNG101:
          if((wmenu.childs+2)->vparam.vflag.mv.ftype!=F_ARCLM)
          {
            MessageBox(NULL,"File Type Error.","ARCLM101",MB_OK);
			break;
		  }
		  if(MessageBox(NULL,"Initial Analysis Arclm001.",
			 "ARCLM101",MB_OKCANCEL)==IDOK)
		  {
			arclm001(&arc,ID_INPUTFILEZ,ID_OUTPUTFILEZ);
		  }

		  if(MessageBox(NULL,"Bclng101:Begin.","BCLNG101",
			 MB_OKCANCEL)==IDOK)
          {
   		        getviewparam((wmenu.childs+2)->hwnd,&((wdraw.childs+1)->vparam));
		        clearwindow(*(wdraw.childs+1));
                //arclm101(&arc,ID_INPUTFILEX);
		        bclng101(&arc,ID_INPUTFILE);
          }
		  break;

		case IDD_GNSHN101:
          /*BEFORE THIS PROCESS:                */
          /*EXTRACT ARCLM FROM ORGAN.           */
          /*SAVE AS ARCLM.                      */
		  /*GRAVITY LOADED ANALYSIS BY ARCLM001.*/
		  if ((wmenu.childs + 2)->vparam.vflag.mv.ftype == F_ARCLM)
		  {
			  if(MessageBox(NULL,"Gnshn101:Begin.","GNSHN101",MB_OKCANCEL)==IDCANCEL)
			  {
				  if(MessageBox(NULL,"Vbrat001:Begin.","VBRAT101",MB_OKCANCEL)==IDCANCEL)
				  {
					break;
				  }
				  else
				  {
					getviewparam((wmenu.childs+2)->hwnd,&((wdraw.childs+1)->vparam));
					clearwindow(*(wdraw.childs+1));
					vbrat001(&arc);
				  }
			  }
			  else
			  {
				getviewparam((wmenu.childs+2)->hwnd,&((wdraw.childs+1)->vparam));
				clearwindow(*(wdraw.childs+1));
				gnshn101(&arc);
			  }
		  }
		  break;

		case IDD_SELECTINPUTFILE:
          GetOpenFileName(&ofn);
          SetDlgItemText(hdwnd,ID_INPUTFILE,strFileTitle);
          return 0;

		case IDD_VIEW:
		  getviewparam((wmenu.childs+2)->hwnd,
					   &((wdraw.childs+1)->vparam));

		  if(wdraw.nchilds>=2 &&
			 (wdraw.childs+1)->hwnd!=NULL &&
			 (wmenu.childs+2)->vparam.vflag.mv.ftype==F_ARCLM)
		  {
			fin=fgetstofopen(dir,"r",ID_INPUTFILE);    /*OPEN FILE.*/
			if(fin==NULL)
			{
			  errormessage("COULD NOT OPEN THE FILE.");
			  break;
			}

			inputinitII(fin,&(arc.nnode),&(arc.nelem),&(arc.nshell),&(arc.nsect),&(arc.nconstraint));
			arc.sects=(struct osect *)malloc(arc.nsect*sizeof(struct osect));
			if(arc.sects==NULL) break;
			arc.nodes=(struct onode *)malloc(arc.nnode*sizeof(struct onode));
			if(arc.nodes==NULL) break;
			arc.ninit=(struct onode *)malloc(arc.nnode*sizeof(struct onode));
			if(arc.ninit==NULL) break;
			arc.elems=(struct owire *)malloc(arc.nelem*sizeof(struct owire));
			if(arc.nelem != 0 && arc.elems==NULL) break;
			arc.shells=(struct oshell *)malloc(arc.nshell*sizeof(struct oshell));
			if(arc.nshell != 0 && arc.shells==NULL) break;
			//arc.melem=(struct memoryelem *)malloc(arc.nelem*sizeof(struct memoryelem));
			//if(arc.nelem != 0 && arc.melem==NULL) break;
			//arc.mshell=(struct memoryshell *)malloc(arc.nshell*sizeof(struct memoryshell));
			//if(arc.nshell != 0 && arc.mshell==NULL) break;

			arc.confs=(struct oconf *)malloc(6*arc.nnode*sizeof(struct oconf));
			if(arc.confs==NULL) break;

			//arc.iform=(double *)malloc(6*arc.nnode*sizeof(double));
			//arc.ddisp=(double *)malloc(6*arc.nnode*sizeof(double));
			arc.constraints=(struct oconstraint *)malloc(arc.nconstraint*sizeof(struct oconstraint));
			if(arc.nconstraint != 0 && arc.constraints==NULL) break;

			/*arc.constraintval=(double *)malloc(arc.nconstraint*sizeof(double));*/        /*MULTI POINT CONSTRAINTS.*/
			/*arc.constraintvec=(double **)malloc(arc.nconstraint*sizeof(double *));
			for(i=0;i<arc.nconstraint;i++)
			{
			  *(arc.constraintvec+i)=(double *)malloc(6*arc.nnode*sizeof(double));
			}*/


			arc.constraintmain=(long int *)malloc(6*arc.nnode*sizeof(long int));
			arc.nmass=(double *)malloc(6*arc.nnode*sizeof(double));



			inputtexttomemory(fin,&arc);
			fclose(fin);

			setviewpoint((wdraw.childs+0)->hwnd,arc,
						 &((wdraw.childs+1)->vparam));
			setviewparam((wmenu.childs+2)->hwnd,
						 (wdraw.childs+1)->vparam);

			clearwindow(*(wdraw.childs+1));

			drawarclmframe((wdraw.childs+1)->hdcC,
						   (wdraw.childs+1)->vparam,arc,0,ONSCREEN);
			SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);
			(wdraw.childs+1)->lstatus=ROTATE;
		  }
		  if(wdraw.nchilds>=2 &&
			 (wdraw.childs+1)->hwnd!=NULL &&
			 (wmenu.childs+2)->vparam.vflag.mv.ftype==F_FRAME)
		  {
			fin=fgetstofopen(dir,"r",ID_INPUTFILE);    /*OPEN FILE.*/
			if(fin==NULL) break;

			inputframeinit(fin,&(arc.nnode),
							   &(arc.nelem),
							   &(arc.nsect));
			arc.sects=(struct osect *)
					  malloc(arc.nsect*sizeof(struct osect));
			if(arc.sects==NULL) break;
			arc.nodes=(struct onode *)
					  malloc(arc.nnode*sizeof(struct onode));
			if(arc.nodes==NULL) break;
			arc.ninit=(struct onode *)
					  malloc(arc.nnode*sizeof(struct onode));
			if(arc.ninit==NULL) break;
			arc.elems=(struct owire *)
					  malloc(arc.nelem*sizeof(struct owire));
			if(arc.elems==NULL) break;
            arc.confs=(struct oconf *)
                      malloc(6*arc.nnode*sizeof(struct oconf));
            if(arc.confs==NULL) break;

            inputframetomemory(fin,&arc);
            fclose(fin);

            clearwindow(*(wdraw.childs+1));
            drawarclmframe((wdraw.childs+1)->hdcC,
                           (wdraw.childs+1)->vparam,arc,0,ONSCREEN);
            SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);

            (wdraw.childs+1)->lstatus=ROTATE;
		  }
		  else if(wdraw.nchilds>=2 &&
				  (wdraw.childs+1)->hwnd!=NULL &&
				  (wmenu.childs+2)->vparam.vflag.mv.ftype==F_ORGAN)
		  {
			freeorganization(&((wdraw.childs+1)->org));

			fin=fopen((wdraw.childs+1)->inpfile,"r");
			if(fin==NULL) break;

			sprintf(str,"OPENED=%s",(wdraw.childs+1)->inpfile);
			errormessage(str);

			inputorganization(fin,&((wdraw.childs+1)->org),
								  &((wdraw.childs+1)->vparam));

		   /*  if(!strncmp(extension((wdraw.childs+1)->inpfile),"dxf",3)) // 150228 fukushima for dxf
			{
			  char data[256];
			  double factor,eps;
			  GetDlgItemText((wmenu.childs+2)->hwnd,ID_LAPS,data,80);
			  factor=strtod(data,NULL);
			  GetDlgItemText((wmenu.childs+2)->hwnd,ID_SAFETY,data,80);
			  eps=strtod(data,NULL);
			  inputorganizationdxf(fin,&((wdraw.childs+1)->org),
									   &((wdraw.childs+1)->vparam),factor,eps);
			}
			else
			{
			  inputorganization(fin,&((wdraw.childs+1)->org),
									&((wdraw.childs+1)->vparam));
			}  */
			fclose(fin);

			/*COPY TEXTS*/
			for(it=0;it<(wdraw.childs+1)->org.ntext;it++)
			{
			  sn=(wdraw.childs+1)->org.texts+it;
			  wdraw.strset=addtext(wdraw.strset,&(wdraw.nstring),
								   sn->str,sn->n.d[0],sn->n.d[1]);
			}

			createviewdata(&((wdraw.childs+1)->vparam));
			SendMessage((wmenu.childs+2)->hwnd,WM_INITDIALOG,0,0);
			SendMessage((wmenu.childs+4)->hwnd,WM_INITDIALOG,0,0);

			clearwindow(*(wdraw.childs+1));
			draworganization((wdraw.childs+1)->hdcC,
							 (wdraw.childs+1)->vparam,
							 (wdraw.childs+1)->org,ONSCREEN);
			SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);

			(wdraw.childs+1)->lstatus=ROTATE;
		  }
		  else if(wdraw.nchilds>=2 &&
				  (wdraw.childs+1)->hwnd!=NULL &&
				  (wmenu.childs+2)->vparam.vflag.mv.ftype==F_HOGAN)
		  {
			fin=fopen((wdraw.childs+1)->inpfile,"r");
			if(fin==NULL) break;

			sprintf(str,"OPENED=%s",(wdraw.childs+1)->inpfile);
			errormessage(str);

			clearwindow(*(wdraw.childs+1));
			drawhoganlines((wdraw.childs+1)->hdcC,
						   (wmenu.childs+5)->hwnd,
						   (wdraw.childs+1)->vparam,
						   fin);
			SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);

			fclose(fin);
		  }
		  else if(wdraw.nchilds>=2 &&
				  (wdraw.childs+1)->hwnd!=NULL &&
				  (wmenu.childs+2)->vparam.vflag.mv.ftype==F_BQUAD)
		  {
			fin=fopen((wdraw.childs+1)->inpfile,"r");
			if(fin==NULL) break;

			sprintf(str,"OPENED=%s",(wdraw.childs+1)->inpfile);
			errormessage(str);

			inputbquadinit(fin,&(bqf.nnode),
							   &(bqf.nwire),&(bqf.nfilm),
							   &(bqf.nwsec),&(bqf.nfsec));

			bqf.sects=(struct qsect *)
					  malloc((bqf.nwsec+bqf.nfsec)
							 *sizeof(struct qsect));
			if(bqf.sects==NULL) break;
			bqf.nodes=(struct node12 *)
					  malloc(bqf.nnode*sizeof(struct node12));
			if(bqf.nodes==NULL) break;
			bqf.ninit=(struct node12 *)
					  malloc(bqf.nnode*sizeof(struct node12));
			if(bqf.ninit==NULL) break;
			if(bqf.nwire>0)
			{
			  bqf.wires=(struct qwire *)
						malloc(bqf.nwire*sizeof(struct qwire));
			  if(bqf.wires==NULL) break;
			}
			if(bqf.nfilm>0)
			{
			  bqf.films=(struct qfilm *)
						malloc(bqf.nfilm*sizeof(struct qfilm));
			  if(bqf.films==NULL) break;
			}
			bqf.confs=(struct oconf *)
					  malloc(12*(bqf.nnode)*sizeof(struct oconf));
			if(bqf.confs==NULL) break;

			inputbiquadframe(fin,&bqf);
			fclose(fin);

			bqf.nnmid=countmidnodes(&bqf);
			bqf.nmids=(struct node09 *)malloc(bqf.nnmid*sizeof(struct node09));
			bqf.confs=(struct oconf *)realloc(bqf.confs,(12*bqf.nnode+9*bqf.nnmid)*sizeof(struct oconf));

			createmidnodes(&bqf);

			clearwindow(*(wdraw.childs+1));
			drawbiquadframe((wdraw.childs+1)->hdcC,
							(wdraw.childs+1)->vparam,bqf,ONSCREEN);
			SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);

			(wdraw.childs+1)->lstatus=ROTATE;
		  }
		  break;

		case IDC_SAVEASARCLM:
		  if(wdraw.nchilds>=2 &&
			 (wdraw.childs+1)->hwnd!=NULL &&
			 (wmenu.childs+2)->vparam.vflag.mv.ftype==F_FRAME &&
			 arc.nnode!=0)
		  {
			saveasarclm("cansav.inp",&arc);
			sprintf((wdraw.childs+1)->inpfile,"cansav.inp");
            SetDlgItemText(hdwnd,ID_INPUTFILE,
                           (wdraw.childs+1)->inpfile);

            (wmenu.childs+2)->vparam.vflag.mv.ftype=F_ARCLM;
            SendMessage(hdwnd,WM_INITDIALOG,0,0);
          }
          if(wdraw.nchilds>=2 &&
             (wdraw.childs+1)->hwnd!=NULL &&
             (wmenu.childs+2)->vparam.vflag.mv.ftype==F_ARCLM &&
             arc.nnode!=0)
          {
            n=strcspn((wdraw.childs+1)->inpfile,".");
            strncpy(str,(wdraw.childs+1)->inpfile,n);
            str[n]='\0';

            strcpy((wdraw.childs+1)->inpfilez,str);
            strcat((wdraw.childs+1)->inpfilez,".inl");
            strcpy((wdraw.childs+1)->inpfilex,str);
            strcat((wdraw.childs+1)->inpfilex,".ihx");
            strcpy((wdraw.childs+1)->inpfiley,str);
            strcat((wdraw.childs+1)->inpfiley,".ihy");

            strcpy((wdraw.childs+1)->otpfile,str);
            strcat((wdraw.childs+1)->otpfile,".otp");
            strcpy((wdraw.childs+1)->otpfilez,str);
            strcat((wdraw.childs+1)->otpfilez,".otl");
            strcpy((wdraw.childs+1)->otpfilex,str);
            strcat((wdraw.childs+1)->otpfilex,".ohx");
            strcpy((wdraw.childs+1)->otpfiley,str);
			strcat((wdraw.childs+1)->otpfiley,".ohy");

			saveasarclm((wdraw.childs+1)->inpfilez,&arc);
            sprintf((wdraw.childs+1)->inpfile,
                    (wdraw.childs+1)->inpfilez);

            SetDlgItemText(hdwnd,ID_INPUTFILE,
                           (wdraw.childs+1)->inpfile);
            SetDlgItemText(hdwnd,ID_INPUTFILEZ,
                           (wdraw.childs+1)->inpfilez);
            SetDlgItemText(hdwnd,ID_INPUTFILEX,
                           (wdraw.childs+1)->inpfilex);
            SetDlgItemText(hdwnd,ID_INPUTFILEY,
                           (wdraw.childs+1)->inpfiley);

            SetDlgItemText(hdwnd,ID_OUTPUTFILE,
                           (wdraw.childs+1)->otpfile);
            SetDlgItemText(hdwnd,ID_OUTPUTFILEZ,
                           (wdraw.childs+1)->otpfilez);
            SetDlgItemText(hdwnd,ID_OUTPUTFILEX,
						   (wdraw.childs+1)->otpfilex);
            SetDlgItemText(hdwnd,ID_OUTPUTFILEY,
						   (wdraw.childs+1)->otpfiley);

			if(arcx.nnode!=0) saveasarclm((wdraw.childs+1)->inpfilex,&arcx);
			if(arcy.nnode!=0) saveasarclm((wdraw.childs+1)->inpfiley,&arcy);

			/*
			char non[80];
			GetDlgItemText((wmenu.childs+2)->hwnd,ID_SECTIONFILE,non,80);
			if(!strncmp(non,"kirigami",8)) strcpy(prj,"kirigami");

			if(!strncmp(prj,"kirigami",8) && MessageBox(NULL,"Save As INL2 FILE","Extract",MB_OKCANCEL)==IDOK)
			{
			   n=strcspn((wdraw.childs+1)->inpfile,".");
			   strncpy(str,(wdraw.childs+1)->inpfile,n);
			   str[n]='\0';

			   strcpy((wdraw.childs+1)->inpfile,str);
			   strcat((wdraw.childs+1)->inpfile,".inl2");
			   SetDlgItemText(hdwnd,ID_INPUTFILE,
							 (wdraw.childs+1)->inpfile);
			   saveasarclm2((wdraw.childs+1)->inpfile,&arc);
			}
			*/
          }
          break;

        case IDD_ROTATE: /*ROTATION*/
          if(wdraw.hwnd!=NULL)
		  {
			(wdraw.childs+1)->lstatus=ROTATE;
            hcursor=LoadCursor(hInstGlobal,"CANCURSORW");   /*CURSOR*/
            SetClassLong((wdraw.childs+1)->hwnd,
                         GCL_HCURSOR,(LONG)hcursor);
            globalstatus=NEUTRAL;
          }
          break;
        case IDD_MOVE:                                /*MOVE FRAME.*/
          if(wdraw.hwnd!=NULL)
          {
            (wdraw.childs+1)->lstatus=MOVE;
            hcursor = LoadCursor(hInstGlobal,"CANBOXW");   /*CURSOR*/
            SetClassLong((wdraw.childs+1)->hwnd,
                         GCL_HCURSOR,(LONG)hcursor);
            globalstatus=NEUTRAL;
          }
          break;

//MODIFIED BY MIHARA 2007.05.22 : GET PRINT RANGE LINE/////////////////////////
        case IDM_PRINTRANGE:
          flagswitch(&((wdraw.childs+1)->vparam.vflag.pv.printrange));
          SendMessage(hdwnd,WM_INITDIALOG,0,0);
//sprintf(str,"printrange=%d",(wdraw.childs+1)->vparam.vflag.pv.printrange);
//MessageBox(NULL,str,"printrange",MB_OK);
        break;

        case IDP_A4TATE:
//sprintf(str,"A4tate-1=%d",(wdraw.childs+1)->vparam.vflag.pv.a4tate);
//MessageBox(NULL,str,"A4",MB_OK);
          flagswitch(&((wdraw.childs+1)->vparam.vflag.pv.a4tate));
          if((wdraw.childs+1)->vparam.vflag.pv.a4tate==1)
          {
          (wdraw.childs+1)->vparam.vflag.pv.a4yoko=0;
          (wdraw.childs+1)->vparam.vflag.pv.a3tate=0;
          (wdraw.childs+1)->vparam.vflag.pv.a3yoko=0;
          }
//sprintf(str,"A4tate-2=%d",(wdraw.childs+1)->vparam.vflag.pv.a4tate);
//MessageBox(NULL,str,"A4",MB_OK);
          SendMessage(hdwnd,WM_INITDIALOG,0,0);
          break;

        case IDP_A4YOKO:
          flagswitch(&((wdraw.childs+1)->vparam.vflag.pv.a4yoko));
          if((wdraw.childs+1)->vparam.vflag.pv.a4yoko==1)
          {
          (wdraw.childs+1)->vparam.vflag.pv.a4tate=0;
          (wdraw.childs+1)->vparam.vflag.pv.a3tate=0;
          (wdraw.childs+1)->vparam.vflag.pv.a3yoko=0;
          }
          SendMessage(hdwnd,WM_INITDIALOG,0,0);
          break;

        case IDP_A3TATE:
          flagswitch(&((wdraw.childs+1)->vparam.vflag.pv.a3tate));
          if((wdraw.childs+1)->vparam.vflag.pv.a3tate==1)
          {
          (wdraw.childs+1)->vparam.vflag.pv.a4tate=0;
          (wdraw.childs+1)->vparam.vflag.pv.a4yoko=0;
          (wdraw.childs+1)->vparam.vflag.pv.a3yoko=0;
          }
          SendMessage(hdwnd,WM_INITDIALOG,0,0);
          break;

        case IDP_A3YOKO:
          flagswitch(&((wdraw.childs+1)->vparam.vflag.pv.a3yoko));
          if((wdraw.childs+1)->vparam.vflag.pv.a3yoko==1)
          {
          (wdraw.childs+1)->vparam.vflag.pv.a4tate=0;
          (wdraw.childs+1)->vparam.vflag.pv.a4yoko=0;
          (wdraw.childs+1)->vparam.vflag.pv.a3tate=0;
          }
          SendMessage(hdwnd,WM_INITDIALOG,0,0);
          break;


///////////////////////////////////////////////////////////////////////////////

      }
      break;

    default:
      return DefWindowProc(hdwnd,message,wParam,lParam);
  }
  return 0;
}/*DialogProcMenu1*/

static BOOL CALLBACK DialogProcMenu2(HWND hdwnd,
                                     UINT message,
                                     WPARAM wParam,
                                     LPARAM lParam)
/*OPTION MENU DIALOG BOX 2.*/
{
  HWND hitem;
  HDC hdc;
  HCURSOR hcursor;
  HFONT hfont;
  /*HBITMAP hbit,pbit;*/
  HBRUSH hbrush,pbrush;
  /*LPDRAWITEMSTRUCT lpdis;*/
  PRINTDLG pd;
  DOCINFO di;
  char str[80],doc[80];
  char dir[]=DIRECTORY;
  int i,id;
  int cWidthPels,cHeightPels,cWidthDpi,cHeightDpi;
  long int code,mode;
  double pfactor; /*PRINT FACTOR*/
  struct print prn;
  struct globvisible *gv;
  struct viewparam vprint;
  FILE *fout;

  switch(message)
  {
    case WM_INITDIALOG:
	  SetDlgItemText(hdwnd,IDV_TXTUNITTM, "[tf],[m]");
      SetDlgItemText(hdwnd,IDV_TXTUNITNM, "SI:[kN],[m]");

      SetDlgItemText(hdwnd,IDV_TXTGLOBALAXIS, "GlobalAxis");
      SetDlgItemText(hdwnd,IDV_TXTNODECODE,   "NodeCode");
		SetDlgItemText(hdwnd, IDV_TXTGX, "Gx"); //mihara for zahyo 20190329
		SetDlgItemText(hdwnd, IDV_TXTGY, "Gy"); //mihara for zahyo 20190329
		SetDlgItemText(hdwnd, IDV_TXTGZ, "Gz"); //mihara for zahyo 20190329
      SetDlgItemText(hdwnd,IDV_TXTLOADS,      "UnitLoads");
      SetDlgItemText(hdwnd,IDV_TXTCONFINEMENT,"Confinement");
      SetDlgItemText(hdwnd,IDV_TXTMASSCIRCLE, "MassCircle");
      SetDlgItemText(hdwnd,IDV_TXTMASSVALUE,  "MassValue");
      SetDlgItemText(hdwnd,IDV_TXTELEMENTCODE,"ElementCode");
      SetDlgItemText(hdwnd,IDV_TXTELEMENTAXIS,"ElementAxis");
      SetDlgItemText(hdwnd,IDV_TXTHINGE,      "Hinge");
	  SetDlgItemText(hdwnd,IDV_TXTSECTIONCODE,"SectionCode");
	  SetDlgItemText(hdwnd,IDV_TXTSECTIONSHAPE,"SectionShape"); //honda
	  SetDlgItemText(hdwnd,IDV_TXTCMQLINE,    "CmqLine");
      SetDlgItemText(hdwnd,IDV_TXTDEFORMATION,"Deformation");
      SetDlgItemText(hdwnd,IDV_TXTDX,"dX");
      SetDlgItemText(hdwnd,IDV_TXTDY,"dY");
      SetDlgItemText(hdwnd,IDV_TXTDZ,"dZ");

      SetDlgItemText(hdwnd,IDV_TXTNZ,"Nz");
      SetDlgItemText(hdwnd,IDV_TXTQX,"Qx");
      SetDlgItemText(hdwnd,IDV_TXTQY,"Qy");
      SetDlgItemText(hdwnd,IDV_TXTMZ,"Mz");
      SetDlgItemText(hdwnd,IDV_TXTMX,"Mx");
      SetDlgItemText(hdwnd,IDV_TXTMY,"My");
      SetDlgItemText(hdwnd,IDV_TXTNZ_G,"Nz");
      SetDlgItemText(hdwnd,IDV_TXTQX_G,"Qx");
      SetDlgItemText(hdwnd,IDV_TXTQY_G,"Qy");
      SetDlgItemText(hdwnd,IDV_TXTMZ_G,"Mz");
      SetDlgItemText(hdwnd,IDV_TXTMX_G,"Mx");
      SetDlgItemText(hdwnd,IDV_TXTMY_G,"My");
      SetDlgItemText(hdwnd,IDV_TXTNZ_B,"Nz");
      SetDlgItemText(hdwnd,IDV_TXTQX_B,"Qx");
      SetDlgItemText(hdwnd,IDV_TXTQY_B,"Qy");
      SetDlgItemText(hdwnd,IDV_TXTMZ_B,"Mz");
      SetDlgItemText(hdwnd,IDV_TXTMX_B,"Mx");
      SetDlgItemText(hdwnd,IDV_TXTMY_B,"My");
      SetDlgItemText(hdwnd,IDV_TXTNZ_W,"Nz");
      SetDlgItemText(hdwnd,IDV_TXTNZ_S,"Nz");

      SetDlgItemText(hdwnd,IDV_TXTCOLUMN,"Column");
      SetDlgItemText(hdwnd,IDV_TXTGIRDER,"Girder");
      SetDlgItemText(hdwnd,IDV_TXTBEAM,  "Beam");
      SetDlgItemText(hdwnd,IDV_TXTBRACE, "Brace");
      SetDlgItemText(hdwnd,IDV_TXTWALL,  "Wall");
      SetDlgItemText(hdwnd,IDV_TXTSLAB,  "Slab");

      SetDlgItemText(hdwnd,IDV_TXTREACTION,"Reaction");

      SetDlgItemText(hdwnd,IDV_TXTINPUTFILE, "InputFile");
      SetDlgItemText(hdwnd,IDV_TXTOUTPUTFILE,"OutputFile");

      SetDlgItemText(hdwnd,IDV_TXTSRCANRATE, "SrcanRate");
      SetDlgItemText(hdwnd,IDV_TXTSRCANCOLOR,"SrcanColor");

      SetDlgItemText(hdwnd,IDV_TXTTITLE,"Title:");
      SetDlgItemText(hdwnd,IDV_TXTCONF,"ConfFig");
      SetDlgItemText(hdwnd,IDV_TXTVIEW,"ViewPoint");

      SetDlgItemText(hdwnd,IDS_CODE,"0");

      SetDlgItemText(hdwnd,IDN_CODE,"0");

      SetDlgItemText(hdwnd,IDE_CODE,"0");

      SetDlgItemText(hdwnd,IDP_MARGINTOP,     "0");
      SetDlgItemText(hdwnd,IDP_MARGINBOTTOM,  "0");
      SetDlgItemText(hdwnd,IDP_MARGINLEFT,  "300");
      SetDlgItemText(hdwnd,IDP_MARGINRIGHT, "200");
      SetDlgItemText(hdwnd,IDP_JIHEIGHT, "50");
      SetDlgItemText(hdwnd,IDP_JIWIDTH,  "20");
      SetDlgItemText(hdwnd,IDP_JIPITCH,   "0");
      SetDlgItemText(hdwnd,IDP_GYOPITCH, "70");
      SetDlgItemText(hdwnd,IDP_DANS,      "2");
      SetDlgItemText(hdwnd,IDP_DANGAP,  "100");

      if(wdraw.nchilds>=2 && (wdraw.childs+1)!=NULL)
      {
        SetDlgItemText(hdwnd,IDV_TITLEBOX,
                       (wdraw.childs+1)->title);
      }
      break;

    case WM_CTLCOLORSTATIC: /*STATIC TEXT MENU.*/
      DefWindowProc(hdwnd,message,wParam,lParam);

      hdc=(HDC)wParam;
      hitem=(HWND)lParam;
      id=GetDlgCtrlID(hitem);

      if(wdraw.hwnd!=NULL) gv=&((wdraw.childs+1)->vparam.vflag);

      if(id==IDV_TXTUNITTM && globalunit!=1.0)
      {SetTextColor(hdc,RGB(255,255,255));}
      if(id==IDV_TXTUNITNM && globalunit!=SIUNIT)
      {SetTextColor(hdc,RGB(255,255,255));}

      if((id==IDV_TXTGLOBALAXIS) &&
         (wdraw.hwnd==NULL || gv->axis!=1))
      {SetTextColor(hdc,RGB(255,255,255));}
      if((id==IDV_TXTNODECODE) &&
         (wdraw.hwnd==NULL || gv->nv.code!=1))
      {SetTextColor(hdc,RGB(255,255,255));}
/*mihara for zahyo 20190329///////////////////////////////////////////////////*/
		if ((id == IDV_TXTGX) && (wdraw.hwnd == NULL || gv->nv.d[0] != 1)
			) {
			SetTextColor(hdc, RGB(255, 255, 255));
		}
		if ((id == IDV_TXTGY) && (wdraw.hwnd == NULL || gv->nv.d[1] != 1)
			) {
			SetTextColor(hdc, RGB(255, 255, 255));
		}
		if ((id == IDV_TXTGZ) && (wdraw.hwnd == NULL || gv->nv.d[2] != 1)
			) {
			SetTextColor(hdc, RGB(255, 255, 255));
		}
/*mihara for zahyo 20190329///////////////////////////////////////////////////*/
      if((id==IDV_TXTLOADS) &&
         (wdraw.hwnd==NULL || gv->nv.loads[0]!=1))
      {SetTextColor(hdc,RGB(255,255,255));}
      if((id==IDV_TXTCONFINEMENT) &&
         (wdraw.hwnd==NULL || gv->nv.confs[0]!=1))
      {SetTextColor(hdc,RGB(255,255,255));}
      if((id==IDV_TXTMASSCIRCLE) &&
         (wdraw.hwnd==NULL || gv->nv.mcircle!=1))
      {SetTextColor(hdc,RGB(255,255,255));}
      if((id==IDV_TXTMASSVALUE) &&
         (wdraw.hwnd==NULL || gv->nv.mvalue!=1))
      {SetTextColor(hdc,RGB(255,255,255));}
      if((id==IDV_TXTELEMENTCODE) &&
         (wdraw.hwnd==NULL || gv->ev.code!=1))
      {SetTextColor(hdc,RGB(255,255,255));}
      if((id==IDV_TXTELEMENTAXIS) &&
         (wdraw.hwnd==NULL || gv->ev.axis!=1))
      {SetTextColor(hdc,RGB(255,255,255));}
      if((id==IDV_TXTHINGE) &&
		 (wdraw.hwnd==NULL || gv->ev.hinge!=1))
	  {SetTextColor(hdc,RGB(255,255,255));}
	  if((id==IDV_TXTSECTIONCODE) &&
		 (wdraw.hwnd==NULL || gv->ev.sectioncode!=1))
	  {SetTextColor(hdc,RGB(255,255,255));}
	  if((id==IDV_TXTSECTIONSHAPE) &&
		 (wdraw.hwnd==NULL || gv->ev.sectionshape!=1))
	  {SetTextColor(hdc,RGB(255,255,255));}   //honda
	  if((id==IDV_TXTCMQLINE) &&
		 (wdraw.hwnd==NULL || gv->ev.cmqline!=1))
      {SetTextColor(hdc,RGB(255,255,255));}
      if((id==IDV_TXTDEFORMATION) &&
         (wdraw.hwnd==NULL || gv->ev.deformation!=1))
      {SetTextColor(hdc,RGB(255,255,255));}

      if((id==IDV_TXTDX) &&
         (wdraw.hwnd==NULL || gv->nv.disps[0]!=1))
      {SetTextColor(hdc,RGB(255,255,255));}
      if((id==IDV_TXTDY) &&
         (wdraw.hwnd==NULL || gv->nv.disps[1]!=1))
      {SetTextColor(hdc,RGB(255,255,255));}
      if((id==IDV_TXTDZ) &&
         (wdraw.hwnd==NULL || gv->nv.disps[2]!=1))
      {SetTextColor(hdc,RGB(255,255,255));}

      if((id==IDV_TXTCOLUMN) &&
         (wdraw.hwnd==NULL || gv->ev.etype[COLUMN]!=1))
      {SetTextColor(hdc,RGB(255,255,255));}
      if((id==IDV_TXTGIRDER) &&
         (wdraw.hwnd==NULL || gv->ev.etype[GIRDER]!=1))
      {SetTextColor(hdc,RGB(255,255,255));}
      if((id==IDV_TXTBEAM) &&
         (wdraw.hwnd==NULL || gv->ev.etype[BEAM]!=1))
      {SetTextColor(hdc,RGB(255,255,255));}
      if((id==IDV_TXTBRACE) &&
         (wdraw.hwnd==NULL || gv->ev.etype[BRACE]!=1))
      {SetTextColor(hdc,RGB(255,255,255));}
      if((id==IDV_TXTWALL) &&
         (wdraw.hwnd==NULL || gv->ev.etype[WALL]!=1))
      {SetTextColor(hdc,RGB(255,255,255));}
      if((id==IDV_TXTSLAB) &&
         (wdraw.hwnd==NULL || gv->ev.etype[SLAB]!=1))
      {SetTextColor(hdc,RGB(255,255,255));}

      if((id==IDV_TXTNZ) &&
         (wdraw.hwnd==NULL || gv->ev.stress[1][0]!=1))
      {SetTextColor(hdc,RGB(255,255,255));}
      if((id==IDV_TXTQX) &&
         (wdraw.hwnd==NULL || gv->ev.stress[1][1]!=1))
      {SetTextColor(hdc,RGB(255,255,255));}
      if((id==IDV_TXTQY) &&
         (wdraw.hwnd==NULL || gv->ev.stress[1][2]!=1))
      {SetTextColor(hdc,RGB(255,255,255));}
      if((id==IDV_TXTMZ) &&
         (wdraw.hwnd==NULL || gv->ev.stress[1][3]!=1))
      {SetTextColor(hdc,RGB(255,255,255));}
      if((id==IDV_TXTMX) &&
         (wdraw.hwnd==NULL || gv->ev.stress[1][4]!=1))
      {SetTextColor(hdc,RGB(255,255,255));}
      if((id==IDV_TXTMY) &&
         (wdraw.hwnd==NULL || gv->ev.stress[1][5]!=1))
      {SetTextColor(hdc,RGB(255,255,255));}

      if((id==IDV_TXTNZ_G) &&
         (wdraw.hwnd==NULL || gv->ev.stress[2][0]!=1))
      {SetTextColor(hdc,RGB(255,255,255));}
      if((id==IDV_TXTQX_G) &&
         (wdraw.hwnd==NULL || gv->ev.stress[2][1]!=1))
      {SetTextColor(hdc,RGB(255,255,255));}
      if((id==IDV_TXTQY_G) &&
         (wdraw.hwnd==NULL || gv->ev.stress[2][2]!=1))
      {SetTextColor(hdc,RGB(255,255,255));}
      if((id==IDV_TXTMZ_G) &&
         (wdraw.hwnd==NULL || gv->ev.stress[2][3]!=1))
      {SetTextColor(hdc,RGB(255,255,255));}
      if((id==IDV_TXTMX_G) &&
         (wdraw.hwnd==NULL || gv->ev.stress[2][4]!=1))
      {SetTextColor(hdc,RGB(255,255,255));}
      if((id==IDV_TXTMY_G) &&
         (wdraw.hwnd==NULL || gv->ev.stress[2][5]!=1))
      {SetTextColor(hdc,RGB(255,255,255));}

      if((id==IDV_TXTNZ_B) &&
         (wdraw.hwnd==NULL || gv->ev.stress[4][0]!=1))
      {SetTextColor(hdc,RGB(255,255,255));}
      if((id==IDV_TXTQX_B) &&
         (wdraw.hwnd==NULL || gv->ev.stress[4][1]!=1))
      {SetTextColor(hdc,RGB(255,255,255));}
      if((id==IDV_TXTQY_B) &&
         (wdraw.hwnd==NULL || gv->ev.stress[4][2]!=1))
      {SetTextColor(hdc,RGB(255,255,255));}
      if((id==IDV_TXTMZ_B) &&
         (wdraw.hwnd==NULL || gv->ev.stress[4][3]!=1))
      {SetTextColor(hdc,RGB(255,255,255));}
      if((id==IDV_TXTMX_B) &&
         (wdraw.hwnd==NULL || gv->ev.stress[4][4]!=1))
      {SetTextColor(hdc,RGB(255,255,255));}
      if((id==IDV_TXTMY_B) &&
         (wdraw.hwnd==NULL || gv->ev.stress[4][5]!=1))
      {SetTextColor(hdc,RGB(255,255,255));}

      if((id==IDV_TXTNZ_W) &&
         (wdraw.hwnd==NULL || gv->ev.stress[WALL][0]!=1))
      {SetTextColor(hdc,RGB(255,255,255));}

      if((id==IDV_TXTNZ_S) &&
         (wdraw.hwnd==NULL || gv->ev.stress[SLAB][0]!=1))
      {SetTextColor(hdc,RGB(255,255,255));}

      if((id==IDV_TXTREACTION) &&
         (wdraw.hwnd==NULL || gv->nv.react[0]!=1))
      {
        SetTextColor(hdc,RGB(255,255,255));
      }

      if((id==IDV_TXTINPUTFILE) &&
         (wdraw.hwnd==NULL || gv->mv.inputfile!=1))
      {
        SetTextColor(hdc,RGB(255,255,255));
      }
      if((id==IDV_TXTOUTPUTFILE) &&
         (wdraw.hwnd==NULL || gv->mv.outputfile!=1))
      {
        SetTextColor(hdc,RGB(255,255,255));
      }

      if((id==IDV_TXTSRCANRATE) &&
         (wdraw.hwnd==NULL || gv->ev.srcanrate!=1))
      {
        SetTextColor(hdc,RGB(255,255,255));
      }

      if((id==IDV_TXTSRCANCOLOR) &&
         (wdraw.hwnd==NULL || gv->ev.srcancolor!=1))
      {
        SetTextColor(hdc,RGB(255,255,255));
      }

      if((id==IDV_TXTTITLE) &&
         (wdraw.hwnd==NULL || gv->mv.title!=1))
      {
        SetTextColor(hdc,RGB(255,255,255));
      }
      if((id==IDV_TXTCONF) &&
         (wdraw.hwnd==NULL || gv->nv.conffig!=1))
      {
        SetTextColor(hdc,RGB(255,255,255));
      }
      if((id==IDV_TXTVIEW) &&
         (wdraw.hwnd==NULL || gv->mv.view!=1))
      {
        SetTextColor(hdc,RGB(255,255,255));
      }
      return (BOOL)(HBRUSH)GetStockObject(LTGRAY_BRUSH);

    case WM_COMMAND:
      if(wdraw.hwnd!=NULL) gv=&((wdraw.childs+1)->vparam.vflag);

      switch(LOWORD(wParam))
      {
/***** Zobun Result *****/
		case IDD_OPENRESULT:
          if(wdraw.nchilds>=2 &&
             (wdraw.childs+1)->hwnd!=NULL &&
			 (wmenu.childs+2)->vparam.vflag.mv.ftype==F_ARCLM ||
			 (wmenu.childs+2)->vparam.vflag.mv.ftype==F_FRAME)
          {
			if(arc.nodes==NULL) break;
			fout=fgetstofopen("\0","r",ID_OUTPUTFILE);  /*OPEN FILE.*/
			if(fout==NULL) break;
			if(globalmessageflag==0 || MessageBox(NULL,"ARCLM001 Result","OPEN RESULT",MB_OKCANCEL)==IDOK)
			{
			   errormessage("test");
			   frameoutputtomemory(fout,&arc);
			}
			else if(MessageBox(NULL,"ARCLM101 Result(Zobun)","OPEN RESULT",MB_OKCANCEL)==IDOK)
			{
			   openarclmlastfileII(fout,&arc);
			}

			fclose(fout);

            arc.nlaps=1;

            clearwindow(*(wdraw.childs+1));
            drawarclmframe((wdraw.childs+1)->hdcC,
                           (wdraw.childs+1)->vparam,arc,0,ONSCREEN);
            SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);

            SetDlgItemText(hdwnd,IDV_MODENUM,"1");

            (wdraw.childs+1)->lstatus=ROTATE;

            (wdraw.childs+1)->vparam.vflag.ev.deformation=1;
          }
          break;
/***** Zobun Result *****/
        case IDD_ZOBUNRESULT:
          if(wdraw.nchilds>=2 &&
             (wdraw.childs+1)->hwnd!=NULL &&
             (wmenu.childs+2)->vparam.vflag.mv.ftype==F_ARCLM ||
             (wmenu.childs+2)->vparam.vflag.mv.ftype==F_FRAME)
          {
            if(arc.nodes==NULL) break;

            fout=fgetstofopen("\0","r",ID_OUTPUTFILE);  /*OPEN FILE.*/
            if(fout==NULL) break;

		    if(MessageBox(NULL,"ARCLM101 Result(Zobun)","OPEN RESULT",MB_OKCANCEL)==IDOK)
            {
//               openarclmlastfile(fout,&arc);
               openarclmlastfileIII(fout,&arc);

            }

            fclose(fout);

            arc.nlaps=1;
/*
            clearwindow(*(wdraw.childs+1));
            drawarclmframe((wdraw.childs+1)->hdcC,
                           (wdraw.childs+1)->vparam,arc,0,ONSCREEN);
            SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);

            SetDlgItemText(hdwnd,IDV_MODENUM,"1");

            (wdraw.childs+1)->lstatus=ROTATE;
            (wdraw.childs+1)->vparam.vflag.ev.deformation=1;
*/          }
        case IDV_UNITTM:
          globalunit=1.0;
          SendMessage(hdwnd,WM_INITDIALOG,0,0);
          break;
        case IDV_UNITNM:
          globalunit=SIUNIT;
          SendMessage(hdwnd,WM_INITDIALOG,0,0);
          break;

        case IDV_GLOBALAXIS:
          if(wdraw.hwnd!=NULL)
          {
            flagswitch(&(gv->axis));
            SendMessage(hdwnd,WM_INITDIALOG,0,0);
          }
          break;
        case IDV_NODECODE:
          if(wdraw.hwnd!=NULL)
          {
            flagswitch(&(gv->nv.code));
            SendMessage(hdwnd,WM_INITDIALOG,0,0);
          }
          break;
/*mihara for zahyo 20190329///////////////////////////////////////////////////*/
		case IDV_GX:
			if (wdraw.hwnd != NULL)
            {
				flagswitch(&(gv->nv.d[0]));
				SendMessage(hdwnd, WM_INITDIALOG, 0, 0);
			}
			break;
		case IDV_GY:
			if (wdraw.hwnd != NULL) {

				flagswitch(&(gv->nv.d[1]));
				SendMessage(hdwnd, WM_INITDIALOG, 0, 0);
			}
			break;
		case IDV_GZ:
			if (wdraw.hwnd != NULL)
            {
				flagswitch(&(gv->nv.d[2]));
				SendMessage(hdwnd, WM_INITDIALOG, 0, 0);
			}
			break;
/*mihara for zahyo 20190329///////////////////////////////////////////////////*/
        case IDV_LOADS:
          if(wdraw.hwnd!=NULL)
          {
            flagswitch(&(gv->nv.loads[0]));
            flagswitch(&(gv->nv.loads[1]));
            flagswitch(&(gv->nv.loads[2]));
            flagswitch(&(gv->nv.loads[3]));
            flagswitch(&(gv->nv.loads[4]));
            flagswitch(&(gv->nv.loads[5]));
            SendMessage(hdwnd,WM_INITDIALOG,0,0);
          }
          break;
        case IDV_CONFINEMENT:
          if(wdraw.hwnd!=NULL)
          {
            flagswitch(&(gv->nv.confs[0]));
            flagswitch(&(gv->nv.confs[1]));
            flagswitch(&(gv->nv.confs[2]));
            flagswitch(&(gv->nv.confs[3]));
            flagswitch(&(gv->nv.confs[4]));
            flagswitch(&(gv->nv.confs[5]));
            SendMessage(hdwnd,WM_INITDIALOG,0,0);
          }
          break;
        case IDV_MASSCIRCLE:
          if(wdraw.hwnd!=NULL)
          {
            if((wdraw.childs+1)->org.loads==NULL)
            {
              gv->nv.mcircle=0;
            }
            else flagswitch(&(gv->nv.mcircle));

            SendMessage(hdwnd,WM_INITDIALOG,0,0);
          }
          break;
        case IDV_MASSVALUE:
          if(wdraw.hwnd!=NULL)
          {
            if((wdraw.childs+1)->org.loads==NULL)
            {
              gv->nv.mvalue=0;
            }
            else flagswitch(&(gv->nv.mvalue));

            SendMessage(hdwnd,WM_INITDIALOG,0,0);
          }
          break;
        case IDV_ELEMENTCODE:
          if(wdraw.hwnd!=NULL)
          {
            flagswitch(&(gv->ev.code));
            SendMessage(hdwnd,WM_INITDIALOG,0,0);
          }
          break;
        case IDV_ELEMENTAXIS:
          if(wdraw.hwnd!=NULL)
          {
            flagswitch(&(gv->ev.axis));
            SendMessage(hdwnd,WM_INITDIALOG,0,0);
          }
          break;
        case IDV_HINGE:
          if(wdraw.hwnd!=NULL)
          {
            flagswitch(&(gv->ev.hinge));
            SendMessage(hdwnd,WM_INITDIALOG,0,0);
		  }
		  break;
		case IDV_SECTIONCODE:
		  if(wdraw.hwnd!=NULL)
		  {
			flagswitch(&(gv->ev.sectioncode));
			SendMessage(hdwnd,WM_INITDIALOG,0,0);
		  }
		  break;
		case IDV_SECTIONSHAPE:
		  if(wdraw.hwnd!=NULL)
		  {
			flagswitch(&(gv->ev.sectionshape));
			SendMessage(hdwnd,WM_INITDIALOG,0,0);  //honda
		  }
		  break;
        case IDV_CMQLINE:
          if(wdraw.hwnd!=NULL)
          {
			flagswitch(&(gv->ev.cmqline));
            SendMessage(hdwnd,WM_INITDIALOG,0,0);
          }
          break;

        case IDV_MODENUM:
          if(wdraw.hwnd!=NULL)
          {
            gv->ev.deformation=0;
            SendMessage(hdwnd,WM_INITDIALOG,0,0);
          }
          break;

        case IDV_DEFORMATION:
          if(wdraw.hwnd!=NULL)
          {
/*
            if((wmenu.childs+2)->vparam.vflag.mv.ftype==F_ARCLM &&
               arc.nlaps==0) break;

            GetDlgItemText(hdwnd,IDV_MODENUM,str,80);
            mode=strtol(str,NULL,10);

            if((wmenu.childs+2)->vparam.vflag.mv.ftype==F_ARCLM &&
			   mode>arc.nlaps) break;
            if((wmenu.childs+2)->vparam.vflag.mv.ftype==F_FRAME &&
               mode>1) break;

			if((wmenu.childs+2)->vparam.vflag.mv.ftype==F_ARCLM &&
			   arc.ddisp==NULL) break;
*/
            /*copyform(&arc,*((arc.eigenvec)+mode-1));*/

            flagswitch(&(gv->ev.deformation));
            SendMessage(hdwnd,WM_INITDIALOG,0,0);
          }
          break;
        case IDV_DX:
          if(wdraw.hwnd!=NULL)
          {
            flagswitch(&(gv->nv.disps[0]));
            SendMessage(hdwnd,WM_INITDIALOG,0,0);
          }
          break;
        case IDV_DY:
          if(wdraw.hwnd!=NULL)
          {
            flagswitch(&(gv->nv.disps[1]));
            SendMessage(hdwnd,WM_INITDIALOG,0,0);
          }
          break;
        case IDV_DZ:
          if(wdraw.hwnd!=NULL)
          {
            flagswitch(&(gv->nv.disps[2]));
            SendMessage(hdwnd,WM_INITDIALOG,0,0);
          }
          break;

        case IDV_COLUMN:
          if(wdraw.hwnd!=NULL)
          {
            flagswitch(&(gv->ev.etype[COLUMN]));
            SendMessage(hdwnd,WM_INITDIALOG,0,0);
          }
          break;
        case IDV_GIRDER:
          if(wdraw.hwnd!=NULL)
          {
            flagswitch(&(gv->ev.etype[GIRDER]));
            SendMessage(hdwnd,WM_INITDIALOG,0,0);
          }
          break;
        case IDV_BEAM:
          if(wdraw.hwnd!=NULL)
          {
            flagswitch(&(gv->ev.etype[BEAM]));
            SendMessage(hdwnd,WM_INITDIALOG,0,0);
          }
          break;
        case IDV_BRACE:
          if(wdraw.hwnd!=NULL)
          {
            flagswitch(&(gv->ev.etype[BRACE]));
            SendMessage(hdwnd,WM_INITDIALOG,0,0);
          }
          break;
        case IDV_WALL:
          if(wdraw.hwnd!=NULL)
          {
            flagswitch(&(gv->ev.etype[WALL]));
            SendMessage(hdwnd,WM_INITDIALOG,0,0);
          }
          break;
        case IDV_SLAB:
          if(wdraw.hwnd!=NULL)
          {
            flagswitch(&(gv->ev.etype[SLAB]));
            SendMessage(hdwnd,WM_INITDIALOG,0,0);
          }
          break;

        case IDV_NZ:
          if(wdraw.hwnd!=NULL)
          {
            flagswitch(&(gv->ev.stress[COLUMN][0]));
            SendMessage(hdwnd,WM_INITDIALOG,0,0);
          }
          break;
        case IDV_QX:
          if(wdraw.hwnd!=NULL)
          {
            flagswitch(&(gv->ev.stress[COLUMN][1]));
            SendMessage(hdwnd,WM_INITDIALOG,0,0);
          }
          break;
        case IDV_QY:
          if(wdraw.hwnd!=NULL)
          {
            flagswitch(&(gv->ev.stress[COLUMN][2]));
            SendMessage(hdwnd,WM_INITDIALOG,0,0);
          }
          break;
        case IDV_MZ:
          if(wdraw.hwnd!=NULL)
          {
            flagswitch(&(gv->ev.stress[COLUMN][3]));
            SendMessage(hdwnd,WM_INITDIALOG,0,0);
          }
          break;
        case IDV_MX:
          if(wdraw.hwnd!=NULL)
          {
            flagswitch(&(gv->ev.stress[COLUMN][4]));
            SendMessage(hdwnd,WM_INITDIALOG,0,0);
          }
          break;
        case IDV_MY:
          if(wdraw.hwnd!=NULL)
          {
            flagswitch(&(gv->ev.stress[COLUMN][5]));
            SendMessage(hdwnd,WM_INITDIALOG,0,0);
          }
          break;

        case IDV_NZ_G:
          if(wdraw.hwnd!=NULL)
          {
            flagswitch(&(gv->ev.stress[GIRDER][0]));
            SendMessage(hdwnd,WM_INITDIALOG,0,0);
          }
          break;
        case IDV_QX_G:
          if(wdraw.hwnd!=NULL)
          {
            flagswitch(&(gv->ev.stress[GIRDER][1]));
            SendMessage(hdwnd,WM_INITDIALOG,0,0);
          }
          break;
        case IDV_QY_G:
          if(wdraw.hwnd!=NULL)
          {
            flagswitch(&(gv->ev.stress[GIRDER][2]));
            SendMessage(hdwnd,WM_INITDIALOG,0,0);
          }
          break;
        case IDV_MZ_G:
          if(wdraw.hwnd!=NULL)
          {
            flagswitch(&(gv->ev.stress[GIRDER][3]));
            SendMessage(hdwnd,WM_INITDIALOG,0,0);
          }
          break;
        case IDV_MX_G:
          if(wdraw.hwnd!=NULL)
          {
            flagswitch(&(gv->ev.stress[GIRDER][4]));
            SendMessage(hdwnd,WM_INITDIALOG,0,0);
          }
          break;
        case IDV_MY_G:
          if(wdraw.hwnd!=NULL)
          {
            flagswitch(&(gv->ev.stress[GIRDER][5]));
            SendMessage(hdwnd,WM_INITDIALOG,0,0);
          }
          break;

        case IDV_NZ_B:
          if(wdraw.hwnd!=NULL)
          {
            flagswitch(&(gv->ev.stress[BRACE][0]));
            SendMessage(hdwnd,WM_INITDIALOG,0,0);
          }
          break;
        case IDV_QX_B:
          if(wdraw.hwnd!=NULL)
          {
            flagswitch(&(gv->ev.stress[BRACE][1]));
            SendMessage(hdwnd,WM_INITDIALOG,0,0);
          }
          break;
        case IDV_QY_B:
          if(wdraw.hwnd!=NULL)
          {
            flagswitch(&(gv->ev.stress[BRACE][2]));
            SendMessage(hdwnd,WM_INITDIALOG,0,0);
          }
          break;
        case IDV_MZ_B:
          if(wdraw.hwnd!=NULL)
          {
            flagswitch(&(gv->ev.stress[BRACE][3]));
            SendMessage(hdwnd,WM_INITDIALOG,0,0);
          }
          break;
        case IDV_MX_B:
          if(wdraw.hwnd!=NULL)
          {
            flagswitch(&(gv->ev.stress[BRACE][4]));
            SendMessage(hdwnd,WM_INITDIALOG,0,0);
          }
          break;
        case IDV_MY_B:
          if(wdraw.hwnd!=NULL)
          {
            flagswitch(&(gv->ev.stress[BRACE][5]));
            SendMessage(hdwnd,WM_INITDIALOG,0,0);
          }
          break;

        case IDV_NZ_W:
          if(wdraw.hwnd!=NULL)
          {
            flagswitch(&(gv->ev.stress[WALL][0]));
            SendMessage(hdwnd,WM_INITDIALOG,0,0);
          }
          break;

        case IDV_NZ_S:
          if(wdraw.hwnd!=NULL)
          {
            flagswitch(&(gv->ev.stress[SLAB][0]));
            SendMessage(hdwnd,WM_INITDIALOG,0,0);
          }
          break;

        case IDV_REACTION:
          if(wdraw.hwnd!=NULL)
          {
            flagswitch(&(gv->nv.react[0]));
            flagswitch(&(gv->nv.react[1]));
            flagswitch(&(gv->nv.react[2]));
            flagswitch(&(gv->nv.react[3]));
            flagswitch(&(gv->nv.react[4]));
			flagswitch(&(gv->nv.react[5]));
            SendMessage(hdwnd,WM_INITDIALOG,0,0);
          }
          break;

        case IDV_INPUTFILE:
          if(wdraw.hwnd!=NULL)
          {
            flagswitch(&(gv->mv.inputfile));
            SendMessage(hdwnd,WM_INITDIALOG,0,0);
          }
          break;
        case IDV_OUTPUTFILE:
          if(wdraw.hwnd!=NULL)
          {
            flagswitch(&(gv->mv.outputfile));
            SendMessage(hdwnd,WM_INITDIALOG,0,0);
          }
          break;

        case IDV_SRCANRATE:
          if(wdraw.hwnd!=NULL)
          {
            flagswitch(&(gv->ev.srcanrate));
            SendMessage(hdwnd,WM_INITDIALOG,0,0);
          }
          break;

        case IDV_SRCANCOLOR:
          if(wdraw.hwnd!=NULL)
          {
            flagswitch(&(gv->ev.srcancolor));
            SendMessage(hdwnd,WM_INITDIALOG,0,0);
          }
          break;

        case IDV_TITLE:
          flagswitch(&(gv->mv.title));
          SendMessage(hdwnd,WM_INITDIALOG,0,0);
          break;
        case IDV_TITLEBOX: /*GET INPUT OF KEYBOARD ONE BY ONE.*/
          if(wdraw.nchilds>=2 && (wdraw.childs+1)!=NULL)
          {
            GetDlgItemText(hdwnd,IDV_TITLEBOX,
                           (wdraw.childs+1)->title,256);
          }
          break;

        case IDV_CONF:
          if(wdraw.hwnd!=NULL)
          {
            flagswitch(&(gv->nv.conffig));
            SendMessage(hdwnd,WM_INITDIALOG,0,0);
          }
          break;
        case IDV_VIEW:
          if(wdraw.hwnd!=NULL)
          {
            flagswitch(&(gv->mv.view));
            SendMessage(hdwnd,WM_INITDIALOG,0,0);
          }
          break;

        case IDD_SECTION: /*SELECT SECTION.*/
          if(wdraw.hwnd!=NULL && arc.sects!=NULL)
          {
            globalstatus = SELECTSECTION;

            GetDlgItemText(hdwnd,IDS_CODE,str,20);
            code=strtol(str,NULL,10);
			getsection(hdwnd,arc.sects,arc.nsect,code);

            clearwindow(*(wdraw.childs+1));
            drawarclmframe((wdraw.childs+1)->hdcC,
                           (wdraw.childs+1)->vparam,
                           arc,code,ONSCREEN);
            SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);

            hcursor = LoadCursor(hInstGlobal,"CANBOXW");   /*CURSOR*/
            SetClassLong((wdraw.childs+1)->hwnd,
                         GCL_HCURSOR,(LONG)hcursor);
          }
          break;
        case IDD_NODE: /*SELECT NODE.*/
          if(wdraw.hwnd!=NULL && arc.nodes!=NULL)
          {
            globalstatus = SELECTNODE;

            GetDlgItemText(hdwnd,IDN_CODE,str,20);
            code=strtol(str,NULL,10);
            getnode(hdwnd,&arc,code);

            clearwindow(*(wdraw.childs+1));
			drawarclmframe((wdraw.childs+1)->hdcC,
						   (wdraw.childs+1)->vparam,
                           arc,code,ONSCREEN);
            SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);

            hcursor = LoadCursor(hInstGlobal,"CANBOXW");   /*CURSOR*/
            SetClassLong((wdraw.childs+1)->hwnd,
                         GCL_HCURSOR,(LONG)hcursor);
          }
          break;
        case IDD_ELEMENT: /*SELECT ELEMENT.*/
          if(wdraw.hwnd!=NULL && arc.elems!=NULL)
          {
            globalstatus = SELECTELEMENT;

            GetDlgItemText(hdwnd,IDE_CODE,str,20);
            code=strtol(str,NULL,10);
            getelement(hdwnd,&arc,code);

            clearwindow(*(wdraw.childs+1));
            drawarclmframe((wdraw.childs+1)->hdcC,
                           (wdraw.childs+1)->vparam,
                           arc,code,ONSCREEN);
            SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);

            hcursor = LoadCursor(hInstGlobal,"CANBOXW");   /*CURSOR*/
            SetClassLong((wdraw.childs+1)->hwnd,
                         GCL_HCURSOR,(LONG)hcursor);
          }
          break;

        case IDD_PRINTFRAME:                         /*PRINT FRAME.*/
          if(wdraw.hwnd!=NULL)
          {
            GetDlgItemText(hdwnd,IDP_MARGINTOP,str,20);
            prn.margin[0]=(int)strtol(str,NULL,10);
            GetDlgItemText(hdwnd,IDP_MARGINBOTTOM,str,20);
            prn.margin[1]=(int)strtol(str,NULL,10);
            GetDlgItemText(hdwnd,IDP_MARGINLEFT,str,20);
            prn.margin[2]=(int)strtol(str,NULL,10);
            GetDlgItemText(hdwnd,IDP_MARGINRIGHT,str,20);
            prn.margin[3]=(int)strtol(str,NULL,10);
            GetDlgItemText(hdwnd,IDP_JIHEIGHT,str,20);
            prn.jiheight=(int)strtol(str,NULL,10);
            GetDlgItemText(hdwnd,IDP_JIWIDTH,str,20);
            prn.jiwidth=(int)strtol(str,NULL,10);

            pd.lStructSize = sizeof(PRINTDLG);
            pd.hwndOwner = (HWND)NULL;
            pd.hDevMode = (HANDLE)NULL;
            pd.hDevNames = (HANDLE)NULL;
            pd.hDC = NULL; /*drawcompati;*/ /*TEST.*/
            pd.Flags = PD_RETURNDC;
            pd.nFromPage = (WORD)NULL;
            pd.nToPage = (WORD)NULL;
            pd.nMinPage = (WORD)NULL;
            pd.nMaxPage = (WORD)NULL;
            pd.nCopies = (WORD)NULL;
            /*pd.nFromPage = 1;
            pd.nToPage = 1;
            pd.nMinPage = 0;
            pd.nMaxPage = 0;
            pd.nCopies = 1;*/
			pd.hInstance = /*(HANDLE)*/(HINSTANCE)NULL;
            pd.lCustData = 0L;
            pd.lpfnPrintHook = (LPPRINTHOOKPROC)NULL;
            pd.lpfnSetupHook = (LPSETUPHOOKPROC)NULL;
            pd.lpPrintTemplateName = (LPSTR)NULL;
			pd.lpSetupTemplateName = (LPSTR)NULL;
            pd.hPrintTemplate = (HANDLE)NULL;
            pd.hSetupTemplate = (HANDLE)NULL;

			if(PrintDlg(&pd) != FALSE)
            {
              if(!(GetDeviceCaps(pd.hDC, RASTERCAPS)
                 & RC_BITBLT))
              {
                MessageBox(NULL,"BitBlt Not Available","Print",
                           MB_OK);
              }

              di.cbSize = sizeof(DOCINFO);
              sprintf(doc,"CanvsPrint");
              di.lpszDocName = doc;
              di.lpszOutput = (LPTSTR)NULL;
              di.lpszDatatype = (LPTSTR) NULL;
              di.fwType = 0;

              cWidthPels  = GetDeviceCaps(pd.hDC,HORZRES);
              cHeightPels = GetDeviceCaps(pd.hDC,VERTRES);

              cWidthDpi  = GetDeviceCaps(pd.hDC,LOGPIXELSX);
              cHeightDpi = GetDeviceCaps(pd.hDC,LOGPIXELSY);
              prn.jiwidth =(int)(((double)prn.jiwidth *cWidthDpi)
                                 /600.0);
              prn.jiheight=(int)(((double)prn.jiheight*cHeightDpi)
                                 /600.0);

              gprn=prn; /*COPY TO GLOBAL*/
              gprn.pflag=0;

              StartDoc(pd.hDC,&di);
              StartPage(pd.hDC);

              pfactor=10.0;
              vprint=(wdraw.childs+1)->vparam;
              vprint.Xo=0.5*cWidthPels;
              vprint.Yo=0.5*cHeightPels;
              vprint.gfactor*=pfactor;

              setfontformat(pd.hDC,prn.jiheight,prn.jiwidth,
                            "MS Mincho",0,0,255);
              SetBkMode(pd.hDC,TRANSPARENT);

			  hbrush = CreateSolidBrush(RGB(255,255,255));
			  pbrush = (HBRUSH)SelectObject(pd.hDC,hbrush);
              PatBlt(pd.hDC,0,0,cWidthPels,cHeightPels,PATCOPY);

              if(((wmenu.childs+2)->vparam.vflag.mv.ftype
                  ==F_ARCLM ||
                  (wmenu.childs+2)->vparam.vflag.mv.ftype
                  ==F_FRAME) &&
                 arc.elems!=NULL)
              {
                drawarclmframe(pd.hDC,vprint,arc,0,ONPRINTER);
              }
              if((wmenu.childs+2)->vparam.vflag.mv.ftype==F_ORGAN &&
                 (wdraw.childs+1)->org.elems!=NULL)
              {
                draworganization(pd.hDC,vprint,
                                 (wdraw.childs+1)->org,ONPRINTER);
              }

              setfontformat(pd.hDC,
                            (int)(1.5*prn.jiheight),
                            (int)(1.5*prn.jiwidth),
                            "MS Mincho",0,0,255);

              for(i=0;i<wdraw.nstring;i++)
              {
                (wdraw.strset+i)->n.d[0]*=pfactor;
                (wdraw.strset+i)->n.d[1]*=pfactor;
              }
              drawtexts(pd.hDC,wdraw.strset,wdraw.nstring,vprint);
              for(i=0;i<wdraw.nstring;i++)
              {
                (wdraw.strset+i)->n.d[0]/=pfactor;
                (wdraw.strset+i)->n.d[1]/=pfactor;
              }

			  hfont=(HFONT)GetCurrentObject(pd.hDC,OBJ_FONT);
              DeleteObject(hfont);

              SelectObject(pd.hDC,pbrush);
              DeleteObject(hbrush);

              EndPage(pd.hDC);
              EndDoc(pd.hDC);
              DeleteDC(pd.hDC);

              if(pd.hDevMode != NULL) GlobalFree(pd.hDevMode);
              if(pd.hDevNames != NULL) GlobalFree(pd.hDevNames);
              MessageBox(hdwnd,"Completed.","Print",MB_OK);
            }
          }
          break;

        case IDD_ABORT:                        /*ABORT APPLICATION.*/
          DestroyWindow(wmain.hwnd);
          break;
      }
      break;

    default:
      return DefWindowProc(hdwnd,message,wParam,lParam);
  }
  return 0;
}/*DialogProcMenu2*/

static BOOL CALLBACK DialogProcMenu3(HWND hdwnd,
                                     UINT message,
                                     WPARAM wParam,
                                     LPARAM lParam)
/*OPTION MENU DIALOG BOX 3:FRAME CREATION.*/
{
  HWND hitem;
  HDC hdc;
  HCURSOR hcursor;
  int id;
  struct globvisible *gv;

  switch(message)
  {
    case WM_INITDIALOG:
      SetDlgItemText(hdwnd,IDC_TXTCREATE,"Creation");

      SetDlgItemText(hdwnd,IDC_TXTPROP,"Property");
      SetDlgItemText(hdwnd,IDC_TXTSECT,"Section");
      SetDlgItemText(hdwnd,IDC_TXTNODE,"Node");
      SetDlgItemText(hdwnd,IDC_TXTELEM,"Element");

      SetDlgItemText(hdwnd,IDC_TXTADD,   "Add");
      SetDlgItemText(hdwnd,IDC_TXTCHANGE,"Change");
      SetDlgItemText(hdwnd,IDC_TXTDELETE,"Delete");
      SetDlgItemText(hdwnd,IDC_TXTREFER, "Refer");

      SetDlgItemText(hdwnd,IDV_TXTENERGYCIRCLE, "EnergyCircle");
      SetDlgItemText(hdwnd,IDV_TXTPLASTICENERGY, ":Ep成分");
      SetDlgItemText(hdwnd,IDV_TXTENERGYVALUE,  "EnergyValue");

      SetDlgItemText(hdwnd,IDV_TXTSRCANMAX,"SrcanMax");

      SetDlgItemText(hdwnd,IDV_TXTPAGETITLE,"PageTitle");
      if(wdraw.nchilds>=2 && (wdraw.childs+1)!=NULL)
      {
        SetDlgItemText(hdwnd,IDV_PAGETITLEBOX,
                       (wdraw.childs+1)->pagetitle);
      }

      if(wdraw.nchilds>=2)
      {
        setdoubleintodialog(hdwnd,IDVS_OPAQUE,
                            (wdraw.childs+1)->org.opaque);

        setdoubleintodialog(hdwnd,IDVS_GPERIOD,
                            (wdraw.childs+1)->org.ai.Tc);
        setdoubleintodialog(hdwnd,IDVS_LOCATE,
                            (wdraw.childs+1)->org.ai.Z);
        setdoubleintodialog(hdwnd,IDVS_BASESHEAR,
                            (wdraw.childs+1)->org.ai.Cox);
        setdoubleintodialog(hdwnd,IDVS_PFACT,
                            (wdraw.childs+1)->org.ai.T1);
      }
      break;

    case WM_CTLCOLORSTATIC: /*STATIC TEXT MENU.*/
      DefWindowProc(hdwnd,message,wParam,lParam);

      hdc=(HDC)wParam;
      hitem=(HWND)lParam;
      id=GetDlgCtrlID(hitem);

      if(wdraw.hwnd!=NULL) gv=&((wdraw.childs+1)->vparam.vflag);

      if((id==IDC_TXTCREATE && createcommand==C_NEUTRAL) ||
         (id==IDC_TXTADD    && createcommand!=C_ADD)     ||
         (id==IDC_TXTCHANGE && createcommand!=C_CHANGE)  ||
         (id==IDC_TXTDELETE && createcommand!=C_DELETE)  ||
         (id==IDC_TXTREFER  && createcommand!=C_REFER)   ||
         (id==IDC_TXTPROP   && createitem!=C_PROPERTY) ||
         (id==IDC_TXTSECT   && createitem!=C_SECTION)  ||
         (id==IDC_TXTNODE   && createitem!=C_NODE)     ||
		 (id==IDC_TXTELEM   && createitem!=C_ELEMENT))
      {
        SetTextColor(hdc,RGB(255,255,255));
      }

      if((id==IDV_TXTENERGYCIRCLE) &&
         (wdraw.hwnd==NULL || gv->ev.ecircle!=1))
      {SetTextColor(hdc,RGB(255,255,255));}
      if((id==IDV_TXTPLASTICENERGY) &&
         (wdraw.hwnd==NULL || gv->ev.plasticenergy!=1))
      {SetTextColor(hdc,RGB(255,255,255));}
      if((id==IDV_TXTENERGYVALUE) &&
         (wdraw.hwnd==NULL || gv->ev.evalue!=1))
      {SetTextColor(hdc,RGB(255,255,255));}
      if((id==IDV_TXTSRCANMAX) &&
         (wdraw.hwnd==NULL || gv->ev.srcanmax!=1))
      {SetTextColor(hdc,RGB(255,255,255));}
      if((id==IDV_TXTPAGETITLE) &&
         (wdraw.hwnd==NULL || gv->mv.pagetitle!=1))
      {SetTextColor(hdc,RGB(255,255,255));}

      return (BOOL)(HBRUSH)GetStockObject(LTGRAY_BRUSH);

    case WM_COMMAND:
      if(wdraw.hwnd!=NULL) gv=&((wdraw.childs+1)->vparam.vflag);

      switch(LOWORD(wParam))
      {
        case IDC_CREATE:
          if(createcommand!=C_NEUTRAL)
          {
            createcommand=C_NEUTRAL;
            createitem=C_NEUTRAL;
            SendMessage(hdwnd,WM_INITDIALOG,0,0);

            hcursor = LoadCursor(hInstGlobal,"CANCURSORW");
            SetClassLong((wdraw.childs+1)->hwnd,
                         GCL_HCURSOR,(LONG)hcursor);
            globalstatus=NEUTRAL;
          }
          break;

        case IDC_PROP:
          createitem=C_PROPERTY;
          SendMessage(hdwnd,WM_INITDIALOG,0,0);
          break;
        case IDC_SECT:
          createitem=C_SECTION;
          SendMessage(hdwnd,WM_INITDIALOG,0,0);
          break;
        case IDC_NODE:
          createitem=C_NODE;
          SendMessage(hdwnd,WM_INITDIALOG,0,0);
          break;
        case IDC_ELEM:
          createitem=C_ELEMENT;
          SendMessage(hdwnd,WM_INITDIALOG,0,0);

          if(createcommand==C_ADD)
          {
            globalstatus = SELECTNODE;
            hcursor = LoadCursor(hInstGlobal,"CANBOXW");
            SetClassLong((wdraw.childs+1)->hwnd,
                         GCL_HCURSOR,(LONG)hcursor);
          }
          break;

        case IDC_ADD:
          createcommand=C_ADD;
          SendMessage(hdwnd,WM_INITDIALOG,0,0);

		  if(createitem==C_ELEMENT)
		  {
			globalstatus = SELECTNODE;
			hcursor = LoadCursor(hInstGlobal,"CANBOXW");
			SetClassLong((wdraw.childs+1)->hwnd,
						 GCL_HCURSOR,(LONG)hcursor);
		  }
		  break;
		case IDC_CHANGE:
		  createcommand=C_CHANGE;
		  SendMessage(hdwnd,WM_INITDIALOG,0,0);
		  break;
		case IDC_DELETE:
		  createcommand=C_DELETE;
		  SendMessage(hdwnd,WM_INITDIALOG,0,0);
		  break;
		case IDC_REFER:
		  createcommand=C_REFER;
		  SendMessage(hdwnd,WM_INITDIALOG,0,0);
		  break;

		case IDVS_OPAQUE:
		  if(wdraw.nchilds>=2)
		  {
            getdoublefromdialog((wmenu.childs+4)->hwnd,IDVS_OPAQUE,
                                &((wdraw.childs+1)->org.opaque));
          }
          break;

        case IDV_PAGETITLE:
          flagswitch(&(gv->mv.pagetitle));
          SendMessage(hdwnd,WM_INITDIALOG,0,0);
          break;
        case IDV_PAGETITLEBOX: /*GET INPUT OF KEYBOARD ONE BY ONE.*/
          if(wdraw.nchilds>=2 && (wdraw.childs+1)!=NULL)
          {
            GetDlgItemText(hdwnd,IDV_PAGETITLEBOX,
                           (wdraw.childs+1)->pagetitle,256);
          }
          break;

        case IDV_ENERGYCIRCLE:
          if(wdraw.hwnd!=NULL)
          {
            /*if(arc.nmass==NULL) gv->ev.ecircle=0;
            else*/ flagswitch(&(gv->ev.ecircle));

            SendMessage(hdwnd,WM_INITDIALOG,0,0);
          }
          break;

        case IDV_PLASTICENERGY:
          if(wdraw.hwnd!=NULL)
          {
            /*if(arc.nmass==NULL) gv->ev.ecircle=0;
            else*/ flagswitch(&(gv->ev.plasticenergy));

            SendMessage(hdwnd,WM_INITDIALOG,0,0);
          }
          break;

        case IDV_ENERGYVALUE:
          if(wdraw.hwnd!=NULL)
          {
			/*if(arc.nmass==NULL) gv->ev.evalue=0;
            else*/ flagswitch(&(gv->ev.evalue));

            SendMessage(hdwnd,WM_INITDIALOG,0,0);
          }
          break;
		case IDV_SRCANMAX:
		  if(wdraw.hwnd!=NULL)
          {
            flagswitch(&(gv->ev.srcanmax));
            if(gv->ev.srcanmax) gv->mv.pagetitle=1; /****SRCANMAX****/
            SendMessage(hdwnd,WM_INITDIALOG,0,0);
          }
          break;

        case IDVS_PFACT:
          if(wdraw.nchilds>=2)
          {
            getdoublefromdialog((wmenu.childs+4)->hwnd,IDVS_PFACT,
                                &((wdraw.childs+1)->org.ai.T1));
          }
        case IDVS_BASESHEAR:
          if(wdraw.nchilds>=2)
          {
            getdoublefromdialog((wmenu.childs+4)->hwnd,IDVS_BASESHEAR,
                                &((wdraw.childs+1)->org.ai.Cox));
          }
        case IDVS_LOCATE:
          if(wdraw.nchilds>=2)
          {
            getdoublefromdialog((wmenu.childs+4)->hwnd,IDVS_LOCATE,
                                &((wdraw.childs+1)->org.ai.Z));
          }
        case IDVS_GPERIOD:
          if(wdraw.nchilds>=2)
          {
            getdoublefromdialog((wmenu.childs+4)->hwnd,IDVS_GPERIOD,
                                &((wdraw.childs+1)->org.ai.Tc));
          }
          break;
      }
      break;

    default:
      return DefWindowProc(hdwnd,message,wParam,lParam);
  }
  return 0;
}/*DialogProcMenu3*/

static BOOL CALLBACK DialogProcMenu4(HWND hdwnd,
                                     UINT message,
                                     WPARAM wParam,
                                     LPARAM lParam)
/*OPTION MENU DIALOG BOX 4.*/
{
  switch(message)
  {
    case WM_INITDIALOG:
      SetDlgItemText(hdwnd,IDH_ASIZE,"10.0");
      SetDlgItemText(hdwnd,IDH_NSIZE,"5.0");
      SetDlgItemText(hdwnd,IDH_TITLEX,"X");
      SetDlgItemText(hdwnd,IDH_TITLEY,"Y");
      SetDlgItemText(hdwnd,IDH_SCALEX,"1.0");
      SetDlgItemText(hdwnd,IDH_SCALEY,"1.0");
	  SetDlgItemText(hdwnd,IDH_XMAX,"1000.0");
      SetDlgItemText(hdwnd,IDH_YMAX,"1000.0");
	  SetDlgItemText(hdwnd,IDH_XMIN,"-1000.0");
	  SetDlgItemText(hdwnd,IDH_YMIN,"-1000.0");
      SetDlgItemText(hdwnd,IDH_PITCHX,"10.0");
      SetDlgItemText(hdwnd,IDH_PITCHY,"10.0");
      break;

    default:
      return DefWindowProc(hdwnd,message,wParam,lParam);
  }
  return 0;
}/*DialogProcMenu4*/

static BOOL CALLBACK DialogProcText(HWND hdwnd,
                                    UINT message,
                                    WPARAM wParam,
                                    LPARAM lParam)
/*POPUP DIALOG BOX FOR TEXT.*/
{
  WPARAM wparam;

  switch(message)
  {
    case WM_COMMAND:
      switch(LOWORD(wParam))
      {
        case IDOK:
        case IDCANCEL:
          wparam = MAKEWPARAM((WORD)IDM_POPTEXTRETURN,(WORD)0);
          SendMessage((wdraw.childs+1)->hwnd,
                      WM_COMMAND,wparam,(WPARAM)0);
          hpopdlg=NULL;
          DestroyWindow(hdwnd);
          break;
      }
      break;

    /*default:
      return DefWindowProc(hdwnd,message,wParam,lParam);*/
  }
  return 0;
}/*DialogProcText*/

static BOOL CALLBACK DialogProcIncrement(HWND hdwnd,
                                         UINT message,
                                         WPARAM wParam,
                                         LPARAM lParam)
/*POPUP DIALOG BOX FOR INCREMENT.*/
{
  WPARAM wparam;

  switch(message)
  {
    case WM_COMMAND:
      switch(LOWORD(wParam))
      {
        case IDOK:
        case IDCANCEL:
		  wparam = MAKEWPARAM((WORD)IDM_POPINCREMENTRETURN,
                              (WORD)0);
          SendMessage((wdraw.childs+1)->hwnd,
                      WM_COMMAND,wparam,(WPARAM)0);
          hpopdlg=NULL;
          DestroyWindow(hdwnd);
          break;
      }
      break;

    /*default:
      return DefWindowProc(hdwnd,message,wParam,lParam);*/
  }
  return 0;
}/*DialogProcIncrement*/

static BOOL CALLBACK DialogProcPsim(HWND hdwnd,
                                    UINT message,
                                    WPARAM wParam,
                                    LPARAM lParam)
/*POPUP DIALOG BOX FOR SIMPLE PROPERTY LIST.*/
{
  HDC hdc;
  HWND hitem;
  HCURSOR hcursor;
  RECT rect;
  HBRUSH hbrush,pbrush;
  /*WPARAM wparam;*/
  char str[80];
  int id;
  long int code,mx,my;

  switch(message)
  {
    case WM_INITDIALOG:
      SetDlgItemText(hdwnd,IDDP_TXTPROPCODE,gstr);

      SetDlgItemText(hdwnd,IDDP_TXTPROPNAME,(gprop->name));

      sprintf(str,"%d",gprop->r);
      SetDlgItemText(hdwnd,IDDP_RED,str);
      sprintf(str,"%d",gprop->g);
      SetDlgItemText(hdwnd,IDDP_GREEN,str);
      sprintf(str,"%d",gprop->b);
      SetDlgItemText(hdwnd,IDDP_BLUE,str);
      break;

    case WM_CTLCOLORBTN:
      hitem=(HWND)lParam;
      id=GetDlgCtrlID(hitem);

      if((id==IDDP_PROPICON))
      {
        hbrush = (HBRUSH)CreateSolidBrush(RGB(gprop->r,
                                              gprop->g,
                                              gprop->b));
		pbrush=(HBRUSH)SetClassLong(hitem,GCL_HBRBACKGROUND,
                                    (LONG)hbrush);
        DeleteObject(pbrush);
        return (BOOL)hbrush;
      }
      else return (BOOL)(HBRUSH)GetStockObject(LTGRAY_BRUSH);

    case WM_PAINT:
      /*DefWindowProc(hdwnd,message,wParam,lParam);*/

      /*PROPERTY CODE*/
      GetDlgItemText(hdwnd,IDDP_TXTPROPCODE,str,20);
      code=strtol(str,NULL,10);

	  /*SECTION COLOR*/
      gprop=currentprops;
      while(code != gprop->code) gprop++;

	  hbrush = (HBRUSH)CreateSolidBrush(RGB(gprop->r,
                                            gprop->g,
                                            gprop->b));
      hitem=GetDlgItem(hdwnd,IDDP_PROPICON);
      getclientsize(hitem,&mx,&my);
      rect.top=0;
      rect.bottom=my;
      rect.left=0;
      rect.right=mx;

      hdc=GetDC(hitem);

      FillRect(hdc,&rect,hbrush);

      DeleteObject(hbrush);
      ReleaseDC(hdwnd,hdc);
      break;

    case WM_COMMAND:
      switch(LOWORD(wParam))
      {
        case IDDP_PROPICON: /*PROPERTY DROP.*/
          GetDlgItemText(hdwnd,IDDP_TXTPROPCODE,str,20);
          code=strtol(str,NULL,10);

          gprop=currentprops;
          while(code != gprop->code) gprop++;

          if(globalstatus==C_PROPERTYDROP)
          {
            hcursor = LoadCursor(hInstGlobal,"CANCURSORW");
            SetClassLong((wsdsp.childs+1)->hwnd,
                         GCL_HCURSOR,(LONG)hcursor);
            SetClassLong((wprop.childs+2+(gprop->loff))->hwnd,
                         GCL_HCURSOR,(LONG)hcursor);

            globalstatus=NEUTRAL;
          }
          else
          {
            hcursor
            =(HCURSOR)createpropertyicon(hInstGlobal,
                                         (wprop.childs+1)->hdcC,
                                         gprop,FALSE);
            SetClassLong((wsdsp.childs+1)->hwnd,
                         GCL_HCURSOR,(LONG)hcursor);
            SetClassLong((wprop.childs+2+(gprop->loff))->hwnd,
                         GCL_HCURSOR,(LONG)hcursor);

            globalstatus=C_PROPERTYDROP;
          }
          break;

        case IDOK: /*COLOR DETERMINATION.*/
          GetDlgItemText(hdwnd,IDDP_TXTPROPCODE,str,20);
          code=strtol(str,NULL,10);

          gprop=currentprops;
          while(code != gprop->code) gprop++;

          GetDlgItemText(hdwnd,IDDP_RED,str,20);
          gprop->r=(int)strtol(str,NULL,10);
          GetDlgItemText(hdwnd,IDDP_GREEN,str,20);
          gprop->g=(int)strtol(str,NULL,10);
          GetDlgItemText(hdwnd,IDDP_BLUE,str,20);
          gprop->b=(int)strtol(str,NULL,10);

          SendMessage(hdwnd,WM_PAINT,0,0);
          break;
      }
      break;

    /*default:
      return DefWindowProc(hdwnd,message,wParam,lParam);*/
  }
  return 0;
}/*DialogProcPsim*/

static BOOL CALLBACK DialogProcSsim(HWND hdwnd,
                                    UINT message,
                                    WPARAM wParam,
                                    LPARAM lParam)
/*POPUP DIALOG BOX FOR SIMPLE SECTION LIST.*/
{
  HDC hdc;
  HWND hitem;
  HCURSOR hcursor;
  RECT rect;
  HBRUSH hbrush,pbrush;
  WPARAM wparam;
  LPDRAWITEMSTRUCT lpdis;
  char str[80];
  int id;
  long int code,mx,my;

  switch(message)
  {
    case WM_INITDIALOG:
      SetDlgItemText(hdwnd,IDDS_TXTSECTCODE,gstr);

	  sprintf(str,"%d",gsect->dcolor.r);
      SetDlgItemText(hdwnd,IDDS_RED,str);
      sprintf(str,"%d",gsect->dcolor.g);
      SetDlgItemText(hdwnd,IDDS_GREEN,str);
      sprintf(str,"%d",gsect->dcolor.b);
      SetDlgItemText(hdwnd,IDDS_BLUE,str);
      break;

	case WM_CTLCOLORSTATIC:
      DefWindowProc(hdwnd,message,wParam,lParam);

      hdc=(HDC)wParam;
      hitem=(HWND)lParam;
      id=GetDlgCtrlID(hitem);

      if((id==IDDS_TXTSECTCODE) && gsect->dflag==0)
      {
        SetTextColor(hdc,RGB(255,255,255));
      }
      return (BOOL)(HBRUSH)GetStockObject(LTGRAY_BRUSH);

    case WM_CTLCOLOREDIT:
      hitem=(HWND)lParam;
      id=GetDlgCtrlID(hitem);

      if((id==IDDS_COLOR))
      {
        hbrush = (HBRUSH)CreateSolidBrush(RGB(gsect->dcolor.r,
                                              gsect->dcolor.g,
                                              gsect->dcolor.b));
        pbrush=(HBRUSH)SetClassLong(hitem,GCL_HBRBACKGROUND,
                                    (LONG)hbrush);
        DeleteObject(pbrush);
        return (BOOL)hbrush;
      }
      else return (BOOL)(HBRUSH)GetStockObject(WHITE_BRUSH);

    case WM_DRAWITEM:
      DefWindowProc(hdwnd,message,wParam,lParam);

      if(wParam==IDDS_SECTICON)
      {
        GetDlgItemText(hdwnd,IDDS_TXTSECTCODE,str,20);
        code=strtol(str,NULL,10);

        gsect=currentsects;
        while(code != gsect->code) gsect++;

        if((gsect->ppc.npcurve)>0)
        {
          lpdis = (LPDRAWITEMSTRUCT)lParam;
          drawiconbmp(lpdis->hDC,0,0,&(gsect->ppc.ici),SRCINVERT);
        }
      }
      break;

    case WM_PAINT:
      /*DefWindowProc(hdwnd,message,wParam,lParam);*/

      /*SECTION CODE*/
      GetDlgItemText(hdwnd,IDDS_TXTSECTCODE,str,20);
      code=strtol(str,NULL,10);

      /*SECTION COLOR*/
      gsect=currentsects;
      while(code != gsect->code) gsect++;

      hbrush = (HBRUSH)CreateSolidBrush(RGB(gsect->dcolor.r,
                                            gsect->dcolor.g,
                                            gsect->dcolor.b));
      hitem=GetDlgItem(hdwnd,IDDS_COLOR);
      getclientsize(hitem,&mx,&my);
      rect.top=0;
      rect.bottom=my;
      rect.left=0;
      rect.right=mx;

      hdc=GetDC(hitem);

      FillRect(hdc,&rect,hbrush);

      DeleteObject(hbrush);
      ReleaseDC(hdwnd,hdc);

      /*SECTION ICON,NAME*/
      if((gsect->ppc.npcurve)>0)
      {
        hitem=GetDlgItem(hdwnd,IDDS_SECTICON);
        hdc=GetDC(hitem);
        /*DrawIcon(hdc,0,0,gsect->ppc.hico);*/
        drawiconbmp(hdc,0,0,&(gsect->ppc.ici),SRCINVERT);
        ReleaseDC(hitem,hdc);

        SetDlgItemText(hdwnd,IDDS_TXTSECTNAME,gsect->ppc.name);
      }
      break;

    case WM_COMMAND:
      switch(LOWORD(wParam))
      {
        case IDDS_SECTCODE:
          GetDlgItemText(hdwnd,IDDS_TXTSECTCODE,str,20);
          code=strtol(str,NULL,10);

          gsect=currentsects;
          while(code != gsect->code) gsect++;

          sprintf(gstr,"%ld",gsect->code);

          GetDlgItemText(hdwnd,IDDS_RED,str,20);
          gsect->dcolor.r=(int)strtol(str,NULL,10);
          GetDlgItemText(hdwnd,IDDS_GREEN,str,20);
          gsect->dcolor.g=(int)strtol(str,NULL,10);
          GetDlgItemText(hdwnd,IDDS_BLUE,str,20);
          gsect->dcolor.b=(int)strtol(str,NULL,10);
          if(gsect->dflag==0) gsect->dflag=1;
          else                gsect->dflag=0;

          SendMessage(hdwnd,WM_INITDIALOG,0,0);
          break;

        case IDDS_SECTICON: /*SECTION DROP.*/
          GetDlgItemText(hdwnd,IDDS_TXTSECTCODE,str,20);
          code=strtol(str,NULL,10);

          gsect=currentsects;
          while(code != gsect->code) gsect++;

          if(globalstatus==C_SECTIONDROP)
          {
            hcursor = LoadCursor(hInstGlobal,"CANCURSORW");
            SetClassLong((wsdsp.childs+1)->hwnd,
                         GCL_HCURSOR,(LONG)hcursor);
            SetClassLong((wsect.childs+2+(gsect->loff))->hwnd,
                         GCL_HCURSOR,(LONG)hcursor);

            globalstatus=NEUTRAL;
          }
          else if((gsect->ppc.npcurve)>0)
          {
            gsect->ppc.ici.fIcon=FALSE; /*TRUE=ICON FALSE=CURSOR*/
            hcursor = CreateIconIndirect(&(gsect->ppc.ici));
            SetClassLong((wsdsp.childs+1)->hwnd,
                         GCL_HCURSOR,(LONG)hcursor);
            SetClassLong((wsect.childs+2+(gsect->loff))->hwnd,
                         GCL_HCURSOR,(LONG)hcursor);

            globalstatus=C_SECTIONDROP;
          }
          else
          {
            wparam = MAKEWPARAM((WORD)IDMS_CLEAR,(WORD)0);
            SendMessage(wsfrm.hwnd,WM_COMMAND,wparam,(WPARAM)0);
          }
          break;

        case IDOK: /*COLOR DETERMINATION.*/
          GetDlgItemText(hdwnd,IDDS_TXTSECTCODE,str,20);
          code=strtol(str,NULL,10);

          gsect=currentsects;
          while(code != gsect->code) gsect++;

          GetDlgItemText(hdwnd,IDDS_RED,str,20);
          gsect->dcolor.r=(int)strtol(str,NULL,10);
          GetDlgItemText(hdwnd,IDDS_GREEN,str,20);
          gsect->dcolor.g=(int)strtol(str,NULL,10);
          GetDlgItemText(hdwnd,IDDS_BLUE,str,20);
          gsect->dcolor.b=(int)strtol(str,NULL,10);

          SendMessage(hdwnd,WM_PAINT,0,0);
          break;
      }
      break;

    /*default:
      return DefWindowProc(hdwnd,message,wParam,lParam);*/
  }
  return 0;
}/*DialogProcSsim*/

static BOOL CALLBACK DialogProcKatakou(HWND hdwnd,
                                       UINT message,
                                       WPARAM wParam,
                                       LPARAM lParam)
/*POPUP DIALOG BOX FOR KATAKOU INPUT.*/
{
  HDC hdc;
  HWND hitem;
  WPARAM wparam;
  int id;

  switch(message)
  {
    case WM_INITDIALOG:
      SetDlgItemText(hdwnd,IDK_POPTXTFB, "FB:Flat Bar");
      SetDlgItemText(hdwnd,IDK_POPTXTL,  "L:Angle");
      SetDlgItemText(hdwnd,IDK_POPTXTH,  "H");
      SetDlgItemText(hdwnd,IDK_POPTXTCT, "CT:Cut T");
      SetDlgItemText(hdwnd,IDK_POPTXTC,  "C:Channel");
      SetDlgItemText(hdwnd,IDK_POPTXTO,  "O");
      SetDlgItemText(hdwnd,IDK_POPTXTBOX,"Box");

      SetDlgItemText(hdwnd,IDK_POPTXTHEIGHT,"H=");
      SetDlgItemText(hdwnd,IDK_POPTXTWIDTH, "B=");
      SetDlgItemText(hdwnd,IDK_POPTXTTF1,   "tf:upper=");
      SetDlgItemText(hdwnd,IDK_POPTXTTF2,   "tf:lower=");
      SetDlgItemText(hdwnd,IDK_POPTXTTW,    "tw=");

      SetDlgItemText(hdwnd,IDK_POPTXTRI1,"Ri1=");
      SetDlgItemText(hdwnd,IDK_POPTXTRI2,"Ri2=");
      SetDlgItemText(hdwnd,IDK_POPTXTRI3,"Ri3=");
      SetDlgItemText(hdwnd,IDK_POPTXTRI4,"Ri4=");
      SetDlgItemText(hdwnd,IDK_POPTXTRO1,"Ro1=");
      SetDlgItemText(hdwnd,IDK_POPTXTRO2,"Ro2=");
      SetDlgItemText(hdwnd,IDK_POPTXTRO3,"Ro3=");
      SetDlgItemText(hdwnd,IDK_POPTXTRO4,"Ro4=");
      break;

    case WM_CTLCOLORSTATIC: /*STATIC TEXT MENU.*/
      DefWindowProc(hdwnd,message,wParam,lParam);

      hdc=(HDC)wParam;
      hitem=(HWND)lParam;
      id=GetDlgCtrlID(hitem);

      if((id==IDK_POPTXTFB  && addpolycurve.type!=CTYPE_FB) ||
         (id==IDK_POPTXTL   && addpolycurve.type!=CTYPE_L)  ||
         (id==IDK_POPTXTH   && addpolycurve.type!=CTYPE_H)  ||
         (id==IDK_POPTXTCT  && addpolycurve.type!=CTYPE_CT) ||
         (id==IDK_POPTXTC   && addpolycurve.type!=CTYPE_C)  ||
         (id==IDK_POPTXTO   && addpolycurve.type!=CTYPE_O)  ||
         (id==IDK_POPTXTBOX && addpolycurve.type!=CTYPE_BOX))
      {
        SetTextColor(hdc,RGB(255,255,255));
      }
      if((addpolycurve.type==CTYPE_FB) &&
         (id==IDK_POPTXTTF1 || id==IDK_POPTXTTF2 ||
          id==IDK_POPTXTTW  ||
          id==IDK_POPTXTRI1 || id==IDK_POPTXTRI2 ||
          id==IDK_POPTXTRI3 || id==IDK_POPTXTRI4 ||
          id==IDK_POPTXTRO2 || id==IDK_POPTXTRO3 ||
          id==IDK_POPTXTRO4))
      {
        SetTextColor(hdc,RGB(255,255,255));
      }
      if((addpolycurve.type==CTYPE_L) &&
         (id==IDK_POPTXTTF1 ||
          id==IDK_POPTXTRI1 || id==IDK_POPTXTRI2 ||
          id==IDK_POPTXTRI4 ||
          id==IDK_POPTXTRO1 || id==IDK_POPTXTRO3))
      {
        SetTextColor(hdc,RGB(255,255,255));
      }
      if((addpolycurve.type==CTYPE_H) &&
         (id==IDK_POPTXTRI2 || id==IDK_POPTXTRI3 ||
          id==IDK_POPTXTRI4 ||
          id==IDK_POPTXTRO1 || id==IDK_POPTXTRO2 ||
          id==IDK_POPTXTRO3 || id==IDK_POPTXTRO4))
      {
        SetTextColor(hdc,RGB(255,255,255));
      }
      if((addpolycurve.type==CTYPE_CT) &&
         (id==IDK_POPTXTTF2 ||
          id==IDK_POPTXTRI2 || id==IDK_POPTXTRI3 ||
          id==IDK_POPTXTRI4 ||
          id==IDK_POPTXTRO1 || id==IDK_POPTXTRO2 ||
          id==IDK_POPTXTRO3 || id==IDK_POPTXTRO4))
      {
        SetTextColor(hdc,RGB(255,255,255));
      }
      if((addpolycurve.type==CTYPE_O) &&
         (id==IDK_POPTXTHEIGHT || id==IDK_POPTXTWIDTH ||
          id==IDK_POPTXTTF1    || id==IDK_POPTXTTF2 ||
          id==IDK_POPTXTTW     ||
          id==IDK_POPTXTRI2    || id==IDK_POPTXTRI3 ||
          id==IDK_POPTXTRI4    ||
          id==IDK_POPTXTRO2    || id==IDK_POPTXTRO3 ||
          id==IDK_POPTXTRO4))
      {
        SetTextColor(hdc,RGB(255,255,255));
      }
      return (BOOL)(HBRUSH)GetStockObject(LTGRAY_BRUSH);

    case WM_COMMAND:
      switch(LOWORD(wParam))
      {
        case IDOK:
        case IDCANCEL:
          wparam = MAKEWPARAM((WORD)IDM_POPSECTIONRETURN,
                              (WORD)0);
          SendMessage((wsdsp.childs+1)->hwnd,
                      WM_COMMAND,wparam,(WPARAM)0);
          hpopdlg=NULL;
          DestroyWindow(hdwnd);
          break;

        case IDK_POPFB:
          addpolycurve.type=CTYPE_FB;
          SendMessage(hdwnd,WM_INITDIALOG,0,0);
          break;
        case IDK_POPL:
          addpolycurve.type=CTYPE_L;
          SendMessage(hdwnd,WM_INITDIALOG,0,0);
          break;
        case IDK_POPH:
          addpolycurve.type=CTYPE_H;
          SendMessage(hdwnd,WM_INITDIALOG,0,0);
          break;
        case IDK_POPCT:
          addpolycurve.type=CTYPE_CT;
          SendMessage(hdwnd,WM_INITDIALOG,0,0);
          break;
        /*case IDK_POPC:
          addpolycurve.type=CTYPE_C;
          SendMessage(hdwnd,WM_INITDIALOG,0,0);
          break;*/
        case IDK_POPO:
          addpolycurve.type=CTYPE_O;
          SendMessage(hdwnd,WM_INITDIALOG,0,0);
          break;
        /*case IDK_POPBOX:
          addpolycurve.type=CTYPE_BOX;
          SendMessage(hdwnd,WM_INITDIALOG,0,0);
          break;*/
      }
      break;

    /*default:
      return DefWindowProc(hdwnd,message,wParam,lParam);*/
  }
  return 0;
}/*DialogProcKatakou*/

static BOOL CALLBACK DialogProcProperty(HWND hdwnd,
                                        UINT message,
                                        WPARAM wParam,
                                        LPARAM lParam)
/*POPUP DIALOG BOX FOR PROPERTY.*/
{
  HDC hdc;
  HWND hitem;
  RECT rect;
  HBRUSH hbrush,pbrush;
  WPARAM wparam;
  char str[20],name[80];
  int id;
  long int mx,my;

  switch(message)
  {
    case WM_INITDIALOG:
      if(globalstatus==C_ADDPROPERTY) gprop=addingprop;

      sprintf(str,"%d",gprop->code);
      SetDlgItemText(hdwnd,IDPP_CODE,str);
      SetDlgItemText(hdwnd,IDPP_CAPTION,gprop->name);

      sprintf(str,"%.3f",gprop->E);
      SetDlgItemText(hdwnd,IDPP_E,str);
      sprintf(str,"%.5f",gprop->poi);
      SetDlgItemText(hdwnd,IDPP_POISSON,str);
      sprintf(str,"%.3f",gprop->hiju);
      SetDlgItemText(hdwnd,IDPP_HIJU,str);

      sprintf(str,"%d",gprop->r);
      SetDlgItemText(hdwnd,IDPP_RED,str);
      sprintf(str,"%d",gprop->g);
      SetDlgItemText(hdwnd,IDPP_GREEN,str);
      sprintf(str,"%d",gprop->b);
      SetDlgItemText(hdwnd,IDPP_BLUE,str);

      ipopflag=1;

      /*SendMessage(hdwnd,WM_PAINT,0,0);*/
      break;

    case WM_CTLCOLOREDIT:
      hitem=(HWND)lParam;
      id=GetDlgCtrlID(hitem);

      if(id==IDPP_COLOR)
      {
        if(globalstatus==C_ADDPROPERTY) gprop=addingprop;

        hbrush = (HBRUSH)CreateSolidBrush(RGB(gprop->r,
                                              gprop->g,
                                              gprop->b));
        pbrush=(HBRUSH)SetClassLong(hitem,GCL_HBRBACKGROUND,
                                    (LONG)hbrush);
        DeleteObject(pbrush);
        return (BOOL)hbrush;
      }
      else return (BOOL)(HBRUSH)GetStockObject(WHITE_BRUSH);

    case WM_PAINT:
      /*DefWindowProc(hdwnd,message,wParam,lParam);*/

      if(globalstatus==C_ADDPROPERTY) gprop=addingprop;

      hbrush = (HBRUSH)CreateSolidBrush(RGB(gprop->r,
                                            gprop->g,
                                            gprop->b));
      hitem=GetDlgItem(hdwnd,IDPP_COLOR);
      getclientsize(hitem,&mx,&my);
      rect.top=0;
      rect.bottom=my;
      rect.left=0;
      rect.right=mx;

      hdc=GetDC(hitem);

      FillRect(hdc,&rect,hbrush);

      DeleteObject(hbrush);
      ReleaseDC(hdwnd,hdc);

      /*DefWindowProc(hdwnd,message,wParam,lParam);*/
      break;

    case WM_COMMAND:
      switch(LOWORD(wParam))
      {
        case IDPP_RED: /*COLOR DETERMINATION.*/
        case IDPP_GREEN:
        case IDPP_BLUE:
          if(ipopflag==1)
          {
            if(globalstatus==C_ADDPROPERTY) gprop=addingprop;

            GetDlgItemText(hdwnd,IDPP_RED,str,20);
            gprop->r=(int)strtol(str,NULL,10);
            GetDlgItemText(hdwnd,IDPP_GREEN,str,20);
            gprop->g=(int)strtol(str,NULL,10);
            GetDlgItemText(hdwnd,IDPP_BLUE,str,20);
            gprop->b=(int)strtol(str,NULL,10);
          }
          SendMessage(hdwnd,WM_PAINT,0,0);
          break;

        case IDOK:
        case IDCANCEL:
          if(globalstatus==C_ADDPROPERTY) gprop=addingprop;

          getproperty(hdwnd,&(gprop->code),
                            name,
                            &(gprop->E),
                            &(gprop->poi),
							&(gprop->hiju));
          gprop->name=(char *)realloc(gprop->name,
                                      (strlen(name)+1)
                                      *sizeof(char));
          strcpy(gprop->name,name);

          wparam = MAKEWPARAM((WORD)IDM_POPPROPERTYRETURN,
                              (WORD)0);
          if(globalstatus==C_ADDPROPERTY)
          {
            SendMessage((wprop.childs+1)->hwnd,
                        WM_COMMAND,wparam,(WPARAM)0);
          }
          else
          {
            SendMessage((wsdsp.childs+1)->hwnd,
                        WM_COMMAND,wparam,(WPARAM)0);
          }
          hpopdlg=NULL;
          DestroyWindow(hdwnd);
          break;
      }
      break;

    /*default:
      return DefWindowProc(hdwnd,message,wParam,lParam);*/
  }
  return 0;
}/*DialogProcProperty*/

static BOOL CALLBACK DialogProcElemType(HWND hdwnd,
                                        UINT message,
                                        WPARAM wParam,
                                        LPARAM lParam)
/*POPUP DIALOG BOX FOR CHANGE ELEMENT TYPE.*/
{
  HWND hitem;
  HDC hdc;
  WPARAM wparam;
  int id;

  switch(message)
  {
    case WM_INITDIALOG:
      SetDlgItemText(hdwnd,IDET_POPTXTETYPENULL,"Null");
      SetDlgItemText(hdwnd,IDET_POPTXTCOLUMN,   "Column");
      SetDlgItemText(hdwnd,IDET_POPTXTGIRDER,   "Girder");
      SetDlgItemText(hdwnd,IDET_POPTXTBEAM,     "Beam");
      SetDlgItemText(hdwnd,IDET_POPTXTBRACE,    "Brace");
      SetDlgItemText(hdwnd,IDET_POPTXTWALL,     "Wall");
      SetDlgItemText(hdwnd,IDET_POPTXTSLAB,     "Slab");
      break;

    case WM_CTLCOLORSTATIC: /*STATIC TEXT MENU.*/
      DefWindowProc(hdwnd,message,wParam,lParam);

      hdc=(HDC)wParam;
      hitem=(HWND)lParam;
      id=GetDlgCtrlID(hitem);

      if((id==IDET_POPTXTETYPENULL && pelem->type!=TYPENULL) ||
         (id==IDET_POPTXTCOLUMN && pelem->type!=COLUMN) ||
         (id==IDET_POPTXTGIRDER && pelem->type!=GIRDER) ||
         (id==IDET_POPTXTBEAM   && pelem->type!=BEAM)   ||
         (id==IDET_POPTXTBRACE  && pelem->type!=BRACE)  ||
         (id==IDET_POPTXTWALL   && pelem->type!=WALL)   ||
         (id==IDET_POPTXTSLAB   && pelem->type!=SLAB))
      {
        SetTextColor(hdc,RGB(255,255,255));
      }
      return (BOOL)(HBRUSH)GetStockObject(LTGRAY_BRUSH);

    case WM_COMMAND:
      switch(LOWORD(wParam))
      {
        case IDOK:
        case IDCANCEL:
          wparam = MAKEWPARAM((WORD)IDM_POPETYPERETURN,(WORD)0);
          SendMessage((wdraw.childs+1)->hwnd,
                      WM_COMMAND,wparam,(WPARAM)0);
          hpopdlg=NULL;
          DestroyWindow(hdwnd);
          break;

        case IDET_POPETYPENULL:
          pelem->type=TYPENULL;
          SendMessage(hdwnd,WM_INITDIALOG,0,0);
          break;
        case IDET_POPCOLUMN:
          pelem->type=COLUMN;
          SendMessage(hdwnd,WM_INITDIALOG,0,0);
          break;
        case IDET_POPGIRDER:
          pelem->type=GIRDER;
          SendMessage(hdwnd,WM_INITDIALOG,0,0);
          break;
        /*case IDET_POPBEAM:
          pelem->type=BEAM;
          SendMessage(hdwnd,WM_INITDIALOG,0,0);
          break;*/
        case IDET_POPBRACE:
          pelem->type=BRACE;
          SendMessage(hdwnd,WM_INITDIALOG,0,0);
          break;
        case IDET_POPWALL:
          pelem->type=WALL;
          SendMessage(hdwnd,WM_INITDIALOG,0,0);
          break;
        case IDET_POPSLAB:
          pelem->type=SLAB;
          SendMessage(hdwnd,WM_INITDIALOG,0,0);
          break;
      }
      break;

    /*default:
      return DefWindowProc(hdwnd,message,wParam,lParam);*/
  }
  return 0;
}/*DialogProcElemType*/

static BOOL CALLBACK DialogProcSectRegist(HWND hdwnd,
                                          UINT message,
                                          WPARAM wParam,
                                          LPARAM lParam)
/*POPUP DIALOG BOX FOR SECTION REGISTRATION.*/
{
  char str[80];
  WPARAM wparam;

  switch(message)
  {
    case WM_INITDIALOG:
      if(globalstatus==C_ADDSECTION) gsect=addingsect;

      sprintf(str,"%ld",gsect->code);
      SetDlgItemText(hdwnd,IDSR_POPCODE,str);
      SetDlgItemText(hdwnd,IDSR_POPCAPTION,gsect->ppc.name);
      break;

    case WM_COMMAND:
      switch(LOWORD(wParam))
      {
        case IDOK:
        case IDCANCEL:
          wparam = MAKEWPARAM((WORD)IDM_POPREGISTRETURN,(WORD)0);
          if(globalstatus==C_ADDSECTION)
          {
            SendMessage((wsect.childs+1)->hwnd,
                        WM_COMMAND,wparam,(WPARAM)0);
          }
          else if(globalstatus==CHANGESECTION ||
                  globalstatus==CHANGESECTMULTIELEM)
          {
            SendMessage((wdraw.childs+1)->hwnd,
                        WM_COMMAND,wparam,(WPARAM)0);
          }
          else
          {
            SendMessage((wsdsp.childs+1)->hwnd,
                        WM_COMMAND,wparam,(WPARAM)0);
          }
          hpopdlg=NULL;
          DestroyWindow(hdwnd);
          break;
      }
      break;

    /*default:
      return DefWindowProc(hdwnd,message,wParam,lParam);*/
  }
  return 0;
}/*DialogProcSectRegist*/

static BOOL CALLBACK DialogProcConf(HWND hdwnd,
                                    UINT message,
                                    WPARAM wParam,
                                    LPARAM lParam)
/*POPUP DIALOG BOX FOR CHANGE BOND.*/
{
  HDC hdc;
  HWND hitem;
  WPARAM wparam;
  char str[80];
  int id;
  long int loff;
  struct oconf *conf;

  switch(message)
  {
    case WM_INITDIALOG:
      SetDlgItemText(hdwnd,IDCC_TXTX,  "X");
      SetDlgItemText(hdwnd,IDCC_TXTY,  "Y");
      SetDlgItemText(hdwnd,IDCC_TXTZ,  "Z");
      SetDlgItemText(hdwnd,IDCC_TXTTX, "Tx");
      SetDlgItemText(hdwnd,IDCC_TXTTY, "Ty");
      SetDlgItemText(hdwnd,IDCC_TXTTZ, "Tz");

      sprintf(str,"%ld",pnode->code);
      SetDlgItemText(hdwnd,IDCC_CODE,str);
      break;

    case WM_CTLCOLORSTATIC: /*STATIC TEXT MENU.*/
      DefWindowProc(hdwnd,message,wParam,lParam);

      hdc=(HDC)wParam;
      hitem=(HWND)lParam;
      id=GetDlgCtrlID(hitem);

      loff=6*pnode->loff;
      conf=(wdraw.childs+1)->org.confs+loff;

      if((id==IDCC_TXTX  && (conf+0)->iconf==0) ||
         (id==IDCC_TXTY  && (conf+1)->iconf==0) ||
         (id==IDCC_TXTZ  && (conf+2)->iconf==0) ||
         (id==IDCC_TXTTX && (conf+3)->iconf==0) ||
         (id==IDCC_TXTTY && (conf+4)->iconf==0) ||
         (id==IDCC_TXTTZ && (conf+5)->iconf==0))
      {
        SetTextColor(hdc,RGB(255,255,255)); /*WHITE:FREE*/
      }
      return (BOOL)(HBRUSH)GetStockObject(LTGRAY_BRUSH);

    case WM_COMMAND:
      loff=6*pnode->loff;
      conf=(wdraw.childs+1)->org.confs+loff;

      switch(LOWORD(wParam))
      {
        case IDCC_X:
          if((conf+0)->iconf==0)      (conf+0)->iconf=1;
          else if((conf+0)->iconf==1) (conf+0)->iconf=0;
          SendMessage(hdwnd,WM_INITDIALOG,0,0);
          break;
        case IDCC_Y:
          if((conf+1)->iconf==0)      (conf+1)->iconf=1;
          else if((conf+1)->iconf==1) (conf+1)->iconf=0;
          SendMessage(hdwnd,WM_INITDIALOG,0,0);
          break;
        case IDCC_Z:
          if((conf+2)->iconf==0)      (conf+2)->iconf=1;
          else if((conf+2)->iconf==1) (conf+2)->iconf=0;
          SendMessage(hdwnd,WM_INITDIALOG,0,0);
          break;
        case IDCC_TX:
          if((conf+3)->iconf==0)      (conf+3)->iconf=1;
          else if((conf+3)->iconf==1) (conf+3)->iconf=0;
          SendMessage(hdwnd,WM_INITDIALOG,0,0);
          break;
        case IDCC_TY:
          if((conf+4)->iconf==0)      (conf+4)->iconf=1;
          else if((conf+4)->iconf==1) (conf+4)->iconf=0;
          SendMessage(hdwnd,WM_INITDIALOG,0,0);
          break;
        case IDCC_TZ:
          if((conf+5)->iconf==0)      (conf+5)->iconf=1;
          else if((conf+5)->iconf==1) (conf+5)->iconf=0;
          SendMessage(hdwnd,WM_INITDIALOG,0,0);
          break;

        /*ENDING PROCESS*/
        case IDOK:
        case IDCANCEL:
          wparam = MAKEWPARAM((WORD)IDM_POPCONFRETURN,(WORD)0);

          SendMessage((wdraw.childs+1)->hwnd,
                      WM_COMMAND,wparam,(WPARAM)0);

          hpopdlg=NULL;
          DestroyWindow(hdwnd);
          break;
      }
      break;

    /*default:
      return DefWindowProc(hdwnd,message,wParam,lParam);*/
  }
  return 0;
}/*DialogProcConf*/

static BOOL CALLBACK DialogProcBond(HWND hdwnd,
                                    UINT message,
                                    WPARAM wParam,
                                    LPARAM lParam)
/*POPUP DIALOG BOX FOR CHANGE BOND.*/
{
  HDC hdc;
  HWND hitem;
  WPARAM wparam;
  char str[80];
  int id,i;
  long int loff;
  signed char *bond;

  switch(message)
  {
    case WM_INITDIALOG:
      SetDlgItemText(hdwnd,IDCB_TXTz,  "z");
      SetDlgItemText(hdwnd,IDCB_TXTx,  "x");
      SetDlgItemText(hdwnd,IDCB_TXTy,  "y");
      SetDlgItemText(hdwnd,IDCB_TXTtz, "tz");
      SetDlgItemText(hdwnd,IDCB_TXTtx, "tx");
      SetDlgItemText(hdwnd,IDCB_TXTty, "ty");

      sprintf(str,"%ld",pelem->code);
      SetDlgItemText(hdwnd,IDCB_ELEM,str);
      sprintf(str,"%ld",pnode->code);
      SetDlgItemText(hdwnd,IDCB_NODE,str);
      break;

    case WM_CTLCOLORSTATIC: /*STATIC TEXT MENU.*/
      DefWindowProc(hdwnd,message,wParam,lParam);

      hdc=(HDC)wParam;
      hitem=(HWND)lParam;
      id=GetDlgCtrlID(hitem);

      for(i=0;i<pelem->nnod;i++)
      {
        if(pnode->code==(*(pelem->nods+i))->code)
        {
          loff=6*i;
          break;
        }
      }
      bond=pelem->bonds+loff;

      if((id==IDCB_TXTz  && *(bond+0)==1) ||
         (id==IDCB_TXTx  && *(bond+1)==1) ||
         (id==IDCB_TXTy  && *(bond+2)==1) ||
         (id==IDCB_TXTtz && *(bond+3)==1) ||
         (id==IDCB_TXTtx && *(bond+4)==1) ||
         (id==IDCB_TXTty && *(bond+5)==1))
      {
        SetTextColor(hdc,RGB(255,255,255)); /*WHITE:HINGE*/
      }
      return (BOOL)(HBRUSH)GetStockObject(LTGRAY_BRUSH);

    case WM_COMMAND:
      for(i=0;i<pelem->nnod;i++)
      {
        if(pnode->code==(*(pelem->nods+i))->code)
        {
          loff=6*i;
          break;
        }
      }
      bond=pelem->bonds+loff;

      switch(LOWORD(wParam))
      {
        case IDCB_z:
          if(*(bond+0)==0)      *(bond+0)=1;
          else if(*(bond+0)==1) *(bond+0)=0;
          SendMessage(hdwnd,WM_INITDIALOG,0,0);
          break;
        case IDCB_x:
          if(*(bond+1)==0)      *(bond+1)=1;
          else if(*(bond+1)==1) *(bond+1)=0;
          SendMessage(hdwnd,WM_INITDIALOG,0,0);
          break;
        case IDCB_y:
          if(*(bond+2)==0)      *(bond+2)=1;
          else if(*(bond+2)==1) *(bond+2)=0;
          SendMessage(hdwnd,WM_INITDIALOG,0,0);
          break;
        case IDCB_tz:
          if(*(bond+3)==0)      *(bond+3)=1;
          else if(*(bond+3)==1) *(bond+3)=0;
          SendMessage(hdwnd,WM_INITDIALOG,0,0);
          break;
        case IDCB_tx:
          if(*(bond+4)==0)      *(bond+4)=1;
          else if(*(bond+4)==1) *(bond+4)=0;
          SendMessage(hdwnd,WM_INITDIALOG,0,0);
          break;
        case IDCB_ty:
          if(*(bond+5)==0)      *(bond+5)=1;
          else if(*(bond+5)==1) *(bond+5)=0;
          SendMessage(hdwnd,WM_INITDIALOG,0,0);
          break;

        /*ENDING PROCESS*/
        case IDOK:
        case IDCANCEL:
          wparam = MAKEWPARAM((WORD)IDM_POPBONDRETURN,(WORD)0);

          SendMessage((wdraw.childs+1)->hwnd,
                      WM_COMMAND,wparam,(WPARAM)0);

          hpopdlg=NULL;
          DestroyWindow(hdwnd);
          break;
      }
      break;

    /*default:
	  return DefWindowProc(hdwnd,message,wParam,lParam);*/
  }
  return 0;
}/*DialogProcBond*/

static BOOL CALLBACK EnumChildProcSheet(HWND hwnd,LPARAM lParam)
{
  WORD wlow;
  WPARAM wparam;

  wlow = IDM_FITPARENT; /*ONLY TO "WINDOWPROCEDUREBACK".*/
  wparam = MAKEWPARAM(wlow,(WORD)0);
  SendMessage(hwnd,WM_COMMAND,wparam,lParam);

  return TRUE;
}/*EnumChildProcSheet*/

VOID APIENTRY popupmenudraw(HWND hwnd,POINT pt)
{
  HMENU hmenu,hmenutrack;

  hmenu=LoadMenu(hInstGlobal,"CANPOP1");
  if(hmenu==NULL) return;

  hmenutrack=GetSubMenu(hmenu,0);

  ClientToScreen(hwnd,&pt);

  TrackPopupMenu(hmenutrack,TPM_TOPALIGN |
                            TPM_LEFTALIGN  |
							TPM_RIGHTBUTTON,
                            pt.x,pt.y,0,hwnd,NULL);
  DestroyMenu(hmenu);
  return;
}/*popupmenudraw*/

int classnamealloc(struct windowparams *wp,char *cname)
{
  wp->classname=(char *)malloc(strlen(cname)+1);
  if(wp->classname==NULL) return 0;

  sprintf(wp->classname,cname);

  return 1;
}/*classnamealloc*/

int classdefinition(HINSTANCE hInstance,
                    char *lpszclassname, /*CLASS NAME.*/
                    WNDPROC lpfnwndproc, /*WINDOW PROCEDURE.*/
                    char *lpszmenuname, /*MENU NAME.*/
					int br,int bg,int bb) /*BACKGROUND BRUSH.*/
{
  WNDCLASS wcl;

  wcl.hInstance = hInstance;
  wcl.lpszClassName = lpszclassname;
  wcl.lpfnWndProc = lpfnwndproc;
  wcl.style = 0;                                    /*STYLE DEFAULT*/
  wcl.hIcon = LoadIcon(hInstance,"CANICON");                 /*ICON*/
  wcl.hCursor = LoadCursor(hInstance,"CANCURSORW");        /*CURSOR*/
  wcl.lpszMenuName = lpszmenuname;                           /*MENU*/
  wcl.cbClsExtra = 0;
  wcl.cbWndExtra = 0;
  wcl.hbrBackground = (HBRUSH)CreateSolidBrush(RGB(br,bg,bb));

  if(!RegisterClass(&wcl)) return 0;

  return 1;
}/*classdefinition*/

HDC createcompati(struct windowparams w)
/*CREATE COMPATIBLE OF HWND.*/
{
  HDC hdc,hdccompati;
  HBRUSH hbrush,hmemorybrush;
  HBITMAP hbit,hmemorybit;
  long int maxX,maxY;

  if(w.hwnd!=NULL)
  {
    getclientsize(w.hwnd,&maxX,&maxY);

    hdc=GetDC(w.hwnd);
    hdccompati = CreateCompatibleDC(hdc);
    hbit = CreateCompatibleBitmap(hdc,maxX,maxY);
	hmemorybit = (HBITMAP)SelectObject(hdccompati,hbit);

    hbrush = (HBRUSH)GetClassLong(w.hwnd,GCL_HBRBACKGROUND);

	hmemorybrush = (HBRUSH)SelectObject(hdccompati,hbrush);
    PatBlt(hdccompati,0,0,maxX,maxY,PATCOPY);

    SetBkMode(hdccompati,TRANSPARENT); /*BACKGROUND OF TEXT.*/

    DeleteObject(hmemorybit);
    DeleteObject(hmemorybrush);
    ReleaseDC(w.hwnd,hdc);
  }
  return hdccompati;
}/*createcompati*/

HWND windowdefinition(struct windowparams *wp,
                      HINSTANCE hInstance,
                      char *lpszclassname,            /*CLASS NAME.*/
                      WNDPROC lpfnwndproc,      /*WINDOW PROCEDURE.*/
                      char *lpszmenuname,              /*MENU NAME.*/
                      int br,int bg,int bb,     /*BACKGROUND BRUSH.*/
					  DWORD dwStyle,                /*WINDOW STYLE.*/
                      int x,int y,int w,int h,     /*POSITION,SIZE.*/
                      HWND hparent)                  /*PARENT HWND.*/
/*DEFINITION OF NAME,HWND,COMPATIBLES.*/
{
  if(lpszclassname!=NULL) classnamealloc(wp,lpszclassname);

  classdefinition(hInstance,wp->classname,
                            lpfnwndproc,
                            lpszmenuname,
                            br,bg,bb);          /*CLASS DEFINITION.*/

  wp->hwnd = CreateWindow(wp->classname,
                          wp->classname, /*CAPTION.*/
						  dwStyle,
                          x,y,w,h,
                          hparent,
                          NULL,
                          hInstance,
                          NULL
                         );
  wp->hdcB=createcompati(*wp); /*BACKGROUND.*/
  wp->hdcC=createcompati(*wp); /*COMPATIBLE.*/

  wp->lstatus=NEUTRAL;
  wp->rstatus=NEUTRAL;
  wp->tx=0;
  wp->ty=0;

  registwindowparams(&wrglobal,wp);

  return wp->hwnd;
}/*windowdefinition*/

int registwindowparams(struct winparamsreg *wrapp,
					   struct windowparams *wpreg)
/*REGIST WINDOWPARAM */
{
  if(wpreg==NULL) return wrapp->nwin;

  (wrapp->nwin)++;
  wrapp->wp=(struct windowparams **)
			realloc(wrapp->wp,
					(wrapp->nwin)*sizeof(struct windowparams *));
  *((wrapp->wp)+(wrapp->nwin-1))=wpreg;

  return wrapp->nwin;
}/*registwindowparams*/

void clearwindow(struct windowparams wp)
/*CLEAR HDC.*/
{
  long int maxX,maxY;

  if(wp.hdcC!=NULL)
  {
	getclientsize(wp.hwnd,&maxX,&maxY);

	PatBlt(wp.hdcC,0,0,maxX,maxY,PATCOPY);
  }
  return;
}/*clearwindow*/

void DrawSunken(HDC hdc,int left,int top,int right,int bottom)
/*DRAW SUNKEN RECTANGLE IN LIGHT GRAY GROUND.*/
{
  HPEN hpen1,hpen2,ppen;

  hpen1=CreatePen(PS_SOLID,1,RGB(150,150,150));            /*SHADOW*/
  ppen = (HPEN)SelectObject(hdc,hpen1);
  MoveToEx(hdc,right,top,NULL);
  LineTo(hdc,left,top);
  LineTo(hdc,left,bottom);
  hpen2=CreatePen(PS_SOLID,1,RGB(255,255,255));            /*LUSTER*/
  SelectObject(hdc,hpen2);
  LineTo(hdc,right,bottom);
  LineTo(hdc,right,top);

  SelectObject(hdc,ppen);
  DeleteObject(hpen1);
  DeleteObject(hpen2);

  return;
}/*DrawSunken*/

RECT rectlocation(int l,int t,int r,int b)
/*LOCATE RECT BY LEFTTOP,RIGHTBOTTOM.*/
{
  RECT rect;

  rect.left  =l;
  rect.top   =t;
  rect.right =r;
  rect.bottom=b;

  return rect;
}/*rectlocation*/

struct windowparams *getwindowparams(HWND hwnd)
{
  int i;

  for(i=0;i<wrglobal.nwin;i++)
  {
	if(hwnd==(*(wrglobal.wp+i))->hwnd) return *(wrglobal.wp+i);
  }
  return NULL;
}/*getwindowparams*/

void vbarlocation(HWND hparent,HWND hchild,RECT *vbar)
{
  RECT rect;
  POINT point;
  long int pw,ph,ch,bmax;

  getclientsize(hparent,&pw,&ph);
  bmax=ph-BARWIDTH-5; /*MAX SIZE OF BAR.*/

  GetWindowRect(hchild,&rect);
  /*cw = rect.right-rect.left;*/
  ch = rect.bottom-rect.top;
  point.x=rect.left;
  point.y=rect.top;

  ScreenToClient(hparent,&point);

  vbar->left  =pw-BARWIDTH-1;
  vbar->right =pw-2;

  if(point.y>=1)
  {
    vbar->top=1;
  }
  else if(ch+point.y-1<=bmax)
  {
    vbar->top=1+bmax*(-point.y+1)/(-point.y+1+bmax);
  }
  else
  {
    vbar->top=1+bmax*(-point.y+1)/ch;
  }

  if(ch+point.y-1<=bmax)
  {
    vbar->bottom=bmax+1;
  }
  else if(point.y>=1)
  {
    vbar->bottom=(bmax+1)-bmax*(ch+point.y-1-bmax)/(ch+point.y-1);
  }
  else
  {
    vbar->bottom=(bmax+1)-bmax*(ch+point.y-1-bmax)/ch;
  }

  return;
}/*vbarlocation*/

void hbarlocation(HWND hparent,HWND hchild,RECT *hbar)
{
  RECT rect;
  POINT point;
  long int pw,ph,cw,bmax;

  getclientsize(hparent,&pw,&ph);
  bmax=pw-BARWIDTH-5; /*MAX SIZE OF BAR.*/

  GetWindowRect(hchild,&rect);
  cw = rect.right-rect.left;
  /*ch = rect.bottom-rect.top;*/
  point.x=rect.left;
  point.y=rect.top;

  ScreenToClient(hparent,&point);

  hbar->top   =ph-BARWIDTH-1;
  hbar->bottom=ph-2;

  if(point.x>=1)
  {
    hbar->left=1;
  }
  else if(cw+point.x-1<=bmax)
  {
    hbar->left=1+bmax*(-point.x+1)/(-point.x+1+bmax);
  }
  else
  {
    hbar->left=1+bmax*(-point.x+1)/cw;
  }

  if(cw+point.x-1<=bmax)
  {
    hbar->right=bmax+1;
  }
  else if(point.x>=1)
  {
    hbar->right=(bmax+1)-bmax*(cw+point.x-1-bmax)/(cw+point.x-1);
  }
  else
  {
    hbar->right=(bmax+1)-bmax*(cw+point.x-1-bmax)/cw;
  }

  return;
}/*hbarlocation*/

void drawvbar(HDC hdc,long int maxX,long int maxY,
              struct windowparams *wp)
{
  HBRUSH hbrush,pbrush;
  HBITMAP hbit,pbit;

  hbit = CreateCompatibleBitmap(hdc,maxX,maxY);
  pbit = (HBITMAP)SelectObject(wp->hdcC,hbit);

  PatBlt(wp->hdcC,maxX-BARWIDTH-1,1,
                  BARWIDTH,maxY-BARWIDTH-5,PATCOPY); /*DELETE BAR.*/

  hbrush = CreateSolidBrush(RGB(255,255,0)); /*REDRAW BAR*/
  pbrush = (HBRUSH)SelectObject(wp->hdcC,hbrush);

  PatBlt(wp->hdcC,(wp->vbar.left),(wp->vbar.top),
                  (wp->vbar.right)-(wp->vbar.left)+1,
                  (wp->vbar.bottom)-(wp->vbar.top)+1,PATINVERT);
  BitBlt(hdc,maxX-BARWIDTH-1,1,BARWIDTH,maxY-BARWIDTH-5,
         wp->hdcC,maxX-BARWIDTH-1,1,SRCCOPY);

  SelectObject(wp->hdcC,pbit);
  SelectObject(wp->hdcC,pbrush);
  DeleteObject(hbit);
  DeleteObject(hbrush);

  return;
}/*drawvbar*/

void drawhbar(HDC hdc,long int maxX,long int maxY,
			  struct windowparams *wp)
{
  HBRUSH hbrush,pbrush;
  HBITMAP hbit,pbit;

  hbit = CreateCompatibleBitmap(hdc,maxX,maxY);
  pbit = (HBITMAP)SelectObject(wp->hdcC,hbit);

  PatBlt(wp->hdcC,1,maxY-BARWIDTH-1,
                  maxX-BARWIDTH-5,BARWIDTH,PATCOPY); /*DELETE BAR.*/

  hbrush = CreateSolidBrush(RGB(255,255,0)); /*REDRAW BAR*/
  pbrush = (HBRUSH)SelectObject(wp->hdcC,hbrush);

  PatBlt(wp->hdcC,(wp->hbar.left),(wp->hbar.top),
                  (wp->hbar.right)-(wp->hbar.left)+1,
                  (wp->hbar.bottom)-(wp->hbar.top)+1,PATINVERT);
  BitBlt(hdc,1,maxY-BARWIDTH-1,maxX-BARWIDTH-5,BARWIDTH,
         wp->hdcC,1,maxY-BARWIDTH-1,SRCCOPY);

  SelectObject(wp->hdcC,pbit);
  SelectObject(wp->hdcC,pbrush);
  DeleteObject(hbit);
  DeleteObject(hbrush);

  return;
}/*drawhbar*/

char flagswitch(char *flag)
{
  if(*flag) *flag=0;
  else      *flag=1;
  return *flag;
}/*flagswitch*/

int findtext(HDC hdc,struct snode *strset,int nstr,
             long int mx,long int my)
{
  SIZE size;
  char *str;
  int i;
  long int x,y;

  for(i=0;i<nstr;i++)
  {
    str=(strset+i)->str;
    x=(long int)((strset+i)->n.d[0]);
    y=(long int)((strset+i)->n.d[1]);
    GetTextExtentPoint32(hdc,str,strlen(str),&size);

    if(x<=mx && mx<=(x+size.cx) &&
       y<=my && my<=(y+size.cy))
    {
      return (i+1);
    }
  }
  return 0;
}/*findtext*/

void drawtexts(HDC hdc,struct snode *strset,int nstr,
               struct viewparam vp)
/*DRAW STRINGS SET ON HDC.*/
{
  char *str;
  int i;
  long int x,y;

  for(i=0;i<nstr;i++)
  {
    str=(strset+i)->str;
    x=(long int)((strset+i)->n.d[0]);
    y=(long int)((strset+i)->n.d[1]);

    TextOut(hdc,(int)(x+vp.Xo),(int)(y+vp.Yo),str,strlen(str));
  }
  return;
}/*drawtexts*/

void printmacro(double xmax,double xmin,
                double ymax,double ymin,
                int idc1,int idg1,int idb1,
                int idc2,int idg2,int idb2,
                int idc3,int idg3,int idb3,
                char *load,char *street)
/*PRINT MACRO.*/
{
  int i;

  /*POSITION SETTING*/
  (wdraw.childs+1)->vparam.focus.d[2]=-10.0;

  (wdraw.childs+1)->vparam.range.max.d[GX]=xmax;
  (wdraw.childs+1)->vparam.range.min.d[GX]=xmin;
  (wdraw.childs+1)->vparam.range.max.d[GY]=ymax;
  (wdraw.childs+1)->vparam.range.min.d[GY]=ymin;
  SendMessage((wmenu.childs+2)->hwnd,WM_INITDIALOG,0,0);

  /*SendMessage((wmenu.childs+2)->hwnd,WM_COMMAND,
              MAKEWPARAM(IDD_VIEW,0),0);
  SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
              MAKEWPARAM(IDD_OPENRESULT,0),0);*/

  sprintf((wdraw.strset+1)->str,
          "Load %s, Model %s, Scale=1:600",load,street);

  /*D,Q*/
  if(idc1) SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                       MAKEWPARAM(idc1,0),0);
  if(idg1) SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                       MAKEWPARAM(idg1,0),0);
  if(idb1) SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                       MAKEWPARAM(idb1,0),0);

  if(idc2) SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                       MAKEWPARAM(idc2,0),0);
  if(idg2) SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                       MAKEWPARAM(idg2,0),0);
  if(idb2) SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                       MAKEWPARAM(idb2,0),0);

  if(idc2==IDV_NZ) sprintf((wdraw.strset+3)->str,": Nz [tf]");
  if(idc2==IDV_QX) sprintf((wdraw.strset+3)->str,": Qx [tf]");
  if(idc2==IDV_QY) sprintf((wdraw.strset+3)->str,": Qy [tf]");
  if(idc2==IDV_MZ) sprintf((wdraw.strset+3)->str,": Mz [tfm]");
  if(idc2==IDV_MX) sprintf((wdraw.strset+3)->str,": Mx [tfm]");
  if(idc2==IDV_MY) sprintf((wdraw.strset+3)->str,": My [tfm]");
  if(idc2==0)      sprintf((wdraw.strset+3)->str,": None");

  if(idg2==IDV_NZ_G) sprintf((wdraw.strset+5)->str,": Nz [tf]");
  if(idg2==IDV_QX_G) sprintf((wdraw.strset+5)->str,": Qx [tf]");
  if(idg2==IDV_QY_G) sprintf((wdraw.strset+5)->str,": Qy [tf]");
  if(idg2==IDV_MZ_G) sprintf((wdraw.strset+5)->str,": Mz [tfm]");
  if(idg2==IDV_MX_G) sprintf((wdraw.strset+5)->str,": Mx [tfm]");
  if(idg2==IDV_MY_G) sprintf((wdraw.strset+5)->str,": My [tfm]");
  if(idg2==0)        sprintf((wdraw.strset+5)->str,": None");

  if(idb2==IDV_NZ_B) sprintf((wdraw.strset+7)->str,": Nz [tf]");
  if(idb2==IDV_QX_B) sprintf((wdraw.strset+7)->str,": Qx [tf]");
  if(idb2==IDV_QY_B) sprintf((wdraw.strset+7)->str,": Qy [tf]");
  if(idb2==IDV_MZ_B) sprintf((wdraw.strset+7)->str,": Mz [tfm]");
  if(idb2==IDV_MX_B) sprintf((wdraw.strset+7)->str,": Mx [tfm]");
  if(idb2==IDV_MY_B) sprintf((wdraw.strset+7)->str,": My [tfm]");
  if(idb2==0)        sprintf((wdraw.strset+7)->str,": None");

  (wdraw.strset+ 0)->n.d[1]=  -3.0;
  (wdraw.strset+ 1)->n.d[1]=   5.0;
  (wdraw.strset+ 2)->n.d[1]=  13.0;
  (wdraw.strset+ 3)->n.d[1]=  13.0;
  (wdraw.strset+ 4)->n.d[1]=  21.0;
  (wdraw.strset+ 5)->n.d[1]=  21.0;
  (wdraw.strset+ 6)->n.d[1]=  29.0;
  (wdraw.strset+ 7)->n.d[1]=  29.0;
  (wdraw.strset+ 8)->n.d[1]= -28.0;
  (wdraw.strset+ 9)->n.d[1]= -46.0;
  (wdraw.strset+10)->n.d[1]= -62.0;
  (wdraw.strset+11)->n.d[1]= -78.0;
  (wdraw.strset+12)->n.d[1]= -94.0;
  (wdraw.strset+13)->n.d[1]=-111.0;
  (wdraw.strset+14)->n.d[1]= -14.0;
  (wdraw.strset+15)->n.d[1]= -14.0;
  (wdraw.strset+16)->n.d[1]= -14.0;
  (wdraw.strset+17)->n.d[1]= -14.0;
  (wdraw.strset+18)->n.d[1]= -14.0;
  (wdraw.strset+19)->n.d[1]= -14.0;
  (wdraw.strset+20)->n.d[1]= -14.0;
  (wdraw.strset+21)->n.d[1]= -14.0;
  (wdraw.strset+22)->n.d[1]= -14.0;
  (wdraw.strset+23)->n.d[1]= -14.0;
  (wdraw.strset+24)->n.d[1]= -14.0;
  (wdraw.strset+25)->n.d[1]= -14.0;

  if(idc2==IDV_DEFORMATION)
  {
    sprintf((wdraw.strset+2)->str,"Displacement %s [cm]",load);
    sprintf((wdraw.strset+3)->str," ");
    sprintf((wdraw.strset+4)->str," ");
    sprintf((wdraw.strset+5)->str," ");
    sprintf((wdraw.strset+6)->str," ");
    sprintf((wdraw.strset+7)->str," ");
    SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                MAKEWPARAM((WORD)IDD_OPENRESULT,(WORD)0),0);
    SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                MAKEWPARAM(idc2,0),0);
  }

  /*SendMessage((wmenu.childs+2)->hwnd,WM_COMMAND,
			  MAKEWPARAM(IDD_VIEW,0),0);*/
  SendMessage(wmain.hwnd,WM_COMMAND,
              MAKEWPARAM((WORD)IDM_PRINTWAIT,(WORD)0),0);

  /*N,M*/
  (wdraw.childs+1)->vparam.focus.d[2]=30.0;
  SendMessage((wmenu.childs+2)->hwnd,WM_INITDIALOG,0,0);

  if(idc3!=IDV_DEFORMATION)
  {
    sprintf((wdraw.strset+2)->str,"Column");
    sprintf((wdraw.strset+4)->str,"Girder");
    sprintf((wdraw.strset+6)->str,"Brace");
  }

  if(idc2) SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                       MAKEWPARAM(idc2,0),0);
  if(idg2) SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                       MAKEWPARAM(idg2,0),0);
  if(idb2) SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                       MAKEWPARAM(idb2,0),0);

  if(idc3) SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                       MAKEWPARAM(idc3,0),0);
  if(idg3) SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                       MAKEWPARAM(idg3,0),0);
  if(idb3) SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                       MAKEWPARAM(idb3,0),0);

  for(i=0;i<wdraw.nstring;i++)
  {
    (wdraw.strset+i)->n.d[1]+=155.0;
  }
  if(idc3==IDV_NZ) sprintf((wdraw.strset+3)->str,": Nz [tf]");
  if(idc3==IDV_QX) sprintf((wdraw.strset+3)->str,": Qx [tf]");
  if(idc3==IDV_QY) sprintf((wdraw.strset+3)->str,": Qy [tf]");
  if(idc3==IDV_MZ) sprintf((wdraw.strset+3)->str,": Mz [tfm]");
  if(idc3==IDV_MX) sprintf((wdraw.strset+3)->str,": Mx [tfm]");
  if(idc3==IDV_MY) sprintf((wdraw.strset+3)->str,": My [tfm]");
  if(idc3==0)      sprintf((wdraw.strset+3)->str,": None");

  if(idg3==IDV_NZ_G) sprintf((wdraw.strset+5)->str,": Nz [tf]");
  if(idg3==IDV_QX_G) sprintf((wdraw.strset+5)->str,": Qx [tf]");
  if(idg3==IDV_QY_G) sprintf((wdraw.strset+5)->str,": Qy [tf]");
  if(idg3==IDV_MZ_G) sprintf((wdraw.strset+5)->str,": Mz [tfm]");
  if(idg3==IDV_MX_G) sprintf((wdraw.strset+5)->str,": Mx [tfm]");
  if(idg3==IDV_MY_G) sprintf((wdraw.strset+5)->str,": My [tfm]");
  if(idg3==0)        sprintf((wdraw.strset+5)->str,": None");

  if(idb3==IDV_NZ_B) sprintf((wdraw.strset+7)->str,": Nz [tf]");
  if(idb3==IDV_QX_B) sprintf((wdraw.strset+7)->str,": Qx [tf]");
  if(idb3==IDV_QY_B) sprintf((wdraw.strset+7)->str,": Qy [tf]");
  if(idb3==IDV_MZ_B) sprintf((wdraw.strset+7)->str,": Mz [tfm]");
  if(idb3==IDV_MX_B) sprintf((wdraw.strset+7)->str,": Mx [tfm]");
  if(idb3==IDV_MY_B) sprintf((wdraw.strset+7)->str,": My [tfm]");
  if(idb3==0)        sprintf((wdraw.strset+7)->str,": None");

  /*SendMessage((wmenu.childs+2)->hwnd,WM_COMMAND,
              MAKEWPARAM(IDD_VIEW,0),0);
  SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
              MAKEWPARAM((WORD)IDD_OPENRESULT,(WORD)0),0);*/
  SendMessage(wmain.hwnd,WM_COMMAND,
              MAKEWPARAM((WORD)IDM_PRINTWAIT,(WORD)0),0);
  SendMessage(wmain.hwnd,WM_COMMAND,
              MAKEWPARAM((WORD)IDM_PRINTOUT,(WORD)0),0);

  return;
}/*printmacro*/

void printmodelmacro(double xmax,double xmin,
                     double ymax,double ymin,
                     int idv11,int idv12,int idv13,
                     int idv21,int idv22,int idv23,
                     char *street)
/*PRINT MODEL MACRO.*/
{
  int i;

  /*POSITION SETTING*/
  (wdraw.childs+1)->vparam.focus.d[2]=-10.0;

  (wdraw.childs+1)->vparam.range.max.d[GX]=xmax;
  (wdraw.childs+1)->vparam.range.min.d[GX]=xmin;
  (wdraw.childs+1)->vparam.range.max.d[GY]=ymax;
  (wdraw.childs+1)->vparam.range.min.d[GY]=ymin;
  SendMessage((wmenu.childs+2)->hwnd,WM_INITDIALOG,0,0);

  /*SendMessage((wmenu.childs+2)->hwnd,WM_COMMAND,
              MAKEWPARAM(IDD_VIEW,0),0);
  SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
              MAKEWPARAM(IDD_OPENRESULT,0),0);*/

  sprintf((wdraw.strset+1)->str,"Model %s, Scale=1:600",street);

  /*NODE,ELEMENT CODE*/
  sprintf((wdraw.strset+2)->str,"Node,Element Code");

  if(idv11) SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                       MAKEWPARAM(idv11,0),0);
  if(idv12) SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                       MAKEWPARAM(idv12,0),0);
  if(idv13) SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                       MAKEWPARAM(idv13,0),0);

  (wdraw.strset+ 0)->n.d[1]=  -3.0;
  (wdraw.strset+ 1)->n.d[1]=   5.0;
  (wdraw.strset+ 2)->n.d[1]=  13.0;

  (wdraw.strset+ 3)->n.d[1]= -28.0;
  (wdraw.strset+ 4)->n.d[1]= -46.0;
  (wdraw.strset+ 5)->n.d[1]= -62.0;
  (wdraw.strset+ 6)->n.d[1]= -78.0;
  (wdraw.strset+ 7)->n.d[1]= -94.0;
  (wdraw.strset+ 8)->n.d[1]=-111.0;

  (wdraw.strset+ 9)->n.d[1]= -14.0;
  (wdraw.strset+10)->n.d[1]= -14.0;
  (wdraw.strset+11)->n.d[1]= -14.0;
  (wdraw.strset+12)->n.d[1]= -14.0;
  (wdraw.strset+13)->n.d[1]= -14.0;
  (wdraw.strset+14)->n.d[1]= -14.0;
  (wdraw.strset+15)->n.d[1]= -14.0;
  (wdraw.strset+16)->n.d[1]= -14.0;
  (wdraw.strset+17)->n.d[1]= -14.0;
  (wdraw.strset+18)->n.d[1]= -14.0;
  (wdraw.strset+19)->n.d[1]= -14.0;
  (wdraw.strset+20)->n.d[1]= -14.0;

  /*SendMessage((wmenu.childs+2)->hwnd,WM_COMMAND,
              MAKEWPARAM(IDD_VIEW,0),0);*/
  SendMessage(wmain.hwnd,WM_COMMAND,
              MAKEWPARAM((WORD)IDM_PRINTWAIT,(WORD)0),0);

  if(idv11) SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                       MAKEWPARAM(idv11,0),0);
  if(idv12) SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                       MAKEWPARAM(idv12,0),0);
  if(idv13) SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                       MAKEWPARAM(idv13,0),0);

  /*SECTION CODE*/
  sprintf((wdraw.strset+2)->str,"Section Code");

  (wdraw.childs+1)->vparam.focus.d[2]=30.0;
  SendMessage((wmenu.childs+2)->hwnd,WM_INITDIALOG,0,0);

  if(idv21) SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                       MAKEWPARAM(idv21,0),0);
  if(idv22) SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                       MAKEWPARAM(idv22,0),0);
  if(idv23) SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                       MAKEWPARAM(idv23,0),0);

  for(i=0;i<wdraw.nstring;i++)
  {
    (wdraw.strset+i)->n.d[1]+=155.0;
  }

  /*SendMessage((wmenu.childs+2)->hwnd,WM_COMMAND,
              MAKEWPARAM(IDD_VIEW,0),0);
  SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
			  MAKEWPARAM((WORD)IDD_OPENRESULT,(WORD)0),0);*/
  SendMessage(wmain.hwnd,WM_COMMAND,
			  MAKEWPARAM((WORD)IDM_PRINTWAIT,(WORD)0),0);
  SendMessage(wmain.hwnd,WM_COMMAND,
              MAKEWPARAM((WORD)IDM_PRINTOUT,(WORD)0),0);

  if(idv21) SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                       MAKEWPARAM(idv21,0),0);
  if(idv22) SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
					   MAKEWPARAM(idv22,0),0);
  if(idv23) SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                       MAKEWPARAM(idv23,0),0);

  return;
}/*printmodelmacro*/

void printarclmfigures(FILE *fin,struct viewparam *vp)
/*TYPE 1:PRINT STRESS BY AXIS.*/
/*TYPE 2:PRINT STRESS BY FLOOR.*/
/*TYPE 3:PRINT DISPLACEMENT BY FLOOR.*/
/*SETTING BEFORE BEGINNING:VISIBLE ELEMENT TYPE*/
/*                         M FACTOR*/
/*                         H SIZE*/
{
  WPARAM wparam;
  char **data,str[256];
  int i,j,k,np1,np2,nstr,flag,ftype;

  char fname[80];
  int nframe;
  double nichi[3],qichi[3],michi[3];

  ftype=1; /*FTYPE:1=STRESS BY AXIS, 2=STRESS BY FLOOR, 3=DISP BY FLOOR*/

  /*PRINTER SETUP*/
  SendMessage(wmain.hwnd,WM_COMMAND,
              MAKEWPARAM((WORD)IDM_PRINTSETUP,(WORD)0),0);

  /*VISIBLE SETTING*/
  wparam = MAKEWPARAM(IDV_GLOBALAXIS,0);
  SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,wparam,0);
  wparam = MAKEWPARAM(IDV_NODECODE,0);
/*mihara for zahyo 20190329///////////////////////////////////////////////////*/
	wparam = MAKEWPARAM(IDV_GX, 0);
	SendMessage((wmenu.childs + 3)->hwnd, WM_COMMAND, wparam, 0);
	wparam = MAKEWPARAM(IDV_GY, 0);
	SendMessage((wmenu.childs + 3)->hwnd, WM_COMMAND, wparam, 0);
	wparam = MAKEWPARAM(IDV_GZ, 0);
	SendMessage((wmenu.childs + 3)->hwnd, WM_COMMAND, wparam, 0);
/*mihara for zahyo 20190329///////////////////////////////////////////////////*/
  SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,wparam,0);
  if(ftype==1)
  {
    wparam = MAKEWPARAM(IDV_SLAB,0);
    SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,wparam,0);
    wparam = MAKEWPARAM(IDV_PAGETITLE,0);
    SendMessage((wmenu.childs+4)->hwnd,WM_COMMAND,wparam,0);

    /*UJIOKA*/
    wparam = MAKEWPARAM(IDV_DEFORMATION,0);
    SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,wparam,0);
    wparam = MAKEWPARAM(IDV_TITLE,0);
    SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,wparam,0);
    wparam = MAKEWPARAM(IDV_VIEW,0);
    SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,wparam,0);
    wparam = MAKEWPARAM(IDV_INPUTFILE,0);
    SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,wparam,0);
    wparam = MAKEWPARAM(IDV_OUTPUTFILE,0);
    SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,wparam,0);
    /*UJIOKA*/
  }
  if(ftype==2)
  {
    wparam = MAKEWPARAM(IDV_TITLE,0);
    SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,wparam,0);
    wparam = MAKEWPARAM(IDV_VIEW,0);
    SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,wparam,0);
    wparam = MAKEWPARAM(IDV_INPUTFILE,0);
    SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,wparam,0);
    wparam = MAKEWPARAM(IDV_OUTPUTFILE,0);
    SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,wparam,0);
  }
  if(ftype==3)
  {
    wparam = MAKEWPARAM(IDV_DEFORMATION,0);
    SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,wparam,0);
    wparam = MAKEWPARAM(IDV_TITLE,0);
    SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,wparam,0);
    wparam = MAKEWPARAM(IDV_VIEW,0);
    SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,wparam,0);
    wparam = MAKEWPARAM(IDV_INPUTFILE,0);
    SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,wparam,0);
    wparam = MAKEWPARAM(IDV_OUTPUTFILE,0);
    SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,wparam,0);
  }

  for(i=0;i<3;i++) /*FOR PERIOD Z,X,Y*/
  {
    fseek(fin,0L,SEEK_SET);

    if(i==0)
    {
      if(ftype==1){np1=0; np2=1;}
      if(ftype==2 || ftype==3){np1=0; np2=0;}

      GetDlgItemText((wmenu.childs+2)->hwnd,ID_INPUTFILEZ,str,256);
      sprintf((wdraw.childs+1)->inpfile,str);
      SetDlgItemText((wmenu.childs+2)->hwnd,ID_INPUTFILE,
                     (wdraw.childs+1)->inpfile);
      GetDlgItemText((wmenu.childs+2)->hwnd,ID_OUTPUTFILEZ,str,256);
      sprintf((wdraw.childs+1)->otpfile,str);
      SetDlgItemText((wmenu.childs+2)->hwnd,ID_OUTPUTFILE,
                     (wdraw.childs+1)->otpfile);
    }
    if(i==1)
    {
      np1=0; np2=0;

      GetDlgItemText((wmenu.childs+2)->hwnd,ID_INPUTFILEX,str,256);
      sprintf((wdraw.childs+1)->inpfile,str);
      SetDlgItemText((wmenu.childs+2)->hwnd,ID_INPUTFILE,
                     (wdraw.childs+1)->inpfile);
      GetDlgItemText((wmenu.childs+2)->hwnd,ID_OUTPUTFILEX,str,256);
      sprintf((wdraw.childs+1)->otpfile,str);
      SetDlgItemText((wmenu.childs+2)->hwnd,ID_OUTPUTFILE,
                     (wdraw.childs+1)->otpfile);

    }
    if(i==2)
    {
      np1=1; np2=1;

      GetDlgItemText((wmenu.childs+2)->hwnd,ID_INPUTFILEY,str,256);
      sprintf((wdraw.childs+1)->inpfile,str);
      SetDlgItemText((wmenu.childs+2)->hwnd,ID_INPUTFILE,
                     (wdraw.childs+1)->inpfile);
      GetDlgItemText((wmenu.childs+2)->hwnd,ID_OUTPUTFILEY,str,256);
      sprintf((wdraw.childs+1)->otpfile,str);
      SetDlgItemText((wmenu.childs+2)->hwnd,ID_OUTPUTFILE,
                     (wdraw.childs+1)->otpfile);
    }

    /*OPEN FILE AND DRAW FRAME.*/
    wparam = MAKEWPARAM(IDD_VIEW,0);
    SendMessage((wmenu.childs+2)->hwnd,WM_COMMAND,wparam,0);
    globalmessageflag=0;
    wparam = MAKEWPARAM(IDD_OPENRESULT,0);
    SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,wparam,0);
    wparam = MAKEWPARAM(IDV_UNITNM,0);
    SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,wparam,0);

    for(j=np1;j<=np2;j++) /*FOR FRAME X,Y*/
    {
      if(j==0)
      {
        flag=0;
        while(!flag)
        {
		  data=fgetsbrk(fin,&nstr);
          if(!strcmp(*(data+0),"XFRAME")) flag=1;
        }
      }
      if(j==1)
      {
        flag=0;
        while(!flag)
        {
          data=fgetsbrk(fin,&nstr);
          if(!strcmp(*(data+0),"YFRAME")) flag=1;
        }
      }

      data=fgetsbrk(fin,&nstr);
      if(nstr==0) return;
      if(!strcmp(*(data+0),"ANGLE"))
      {
        vp->theta=strtod(*(data+1),NULL);
        vp->phi  =strtod(*(data+2),NULL);
      }
      freestr(data,nstr);

      data=fgetsbrk(fin,&nstr);
      if(nstr==0) return;
      if(!strcmp(*(data+0),"DISTS"))
      {
        vp->r  =strtod(*(data+1),NULL);
        vp->odv=strtod(*(data+2),NULL);
      }
      freestr(data,nstr);

	  data=fgetsbrk(fin,&nstr);
      if(nstr==0) return;
      if(!strcmp(*(data+0),"NICHI"))
      {
        nichi[0]=strtod(*(data+1),NULL);
        nichi[1]=strtod(*(data+2),NULL);
        nichi[2]=strtod(*(data+3),NULL);
      }
      freestr(data,nstr);

      data=fgetsbrk(fin,&nstr);
      if(nstr==0) return;
      if(!strcmp(*(data+0),"QICHI"))
      {
        qichi[0]=strtod(*(data+1),NULL);
        qichi[1]=strtod(*(data+2),NULL);
        qichi[2]=strtod(*(data+3),NULL);
      }
      freestr(data,nstr);

      data=fgetsbrk(fin,&nstr);
      if(nstr==0) return;
      if(!strcmp(*(data+0),"MICHI"))
      {
        michi[0]=strtod(*(data+1),NULL);
        michi[1]=strtod(*(data+2),NULL);
        michi[2]=strtod(*(data+3),NULL);
      }
      freestr(data,nstr);

      data=fgetsbrk(fin,&nstr);
      if(nstr==0) return;
      if(!strcmp(*(data+0),"NFRAME"))
      {
        nframe=(int)strtol(*(data+1),NULL,10);
      }
      freestr(data,nstr);

      for(k=0;k<nframe;k++)
      {
        data=fgetsbrk(fin,&nstr);
        if(nstr==0) return;

        if(!strcmp(*(data+0),"FRAME"))
        {
          if(ftype==1)
          {
            sprintf(fname,"%s 通り応力図",*(data+1));
            SetDlgItemText((wmenu.childs+4)->hwnd,IDV_PAGETITLEBOX,
                           fname);
          }
          if(ftype==2 || ftype==3)
          {
            sprintf(fname,"%s",*(data+1));
          }

        }
        if(!strcmp(*(data+2),"HEAD"))
        {
          vp->range.min.d[GX]=strtod(*(data+3),NULL);
          vp->range.min.d[GY]=strtod(*(data+4),NULL);
          vp->range.min.d[GZ]=strtod(*(data+5),NULL);
        }
        if(!strcmp(*(data+6),"TAIL"))
        {
          vp->range.max.d[GX]=strtod(*(data+7),NULL);
          vp->range.max.d[GY]=strtod(*(data+8),NULL);
          vp->range.max.d[GZ]=strtod(*(data+9),NULL);
        }

        SendMessage((wmenu.childs+2)->hwnd,WM_INITDIALOG,0,0);

        /*DISPLACEMENTN*/
        vp->focus.d[0]=nichi[0];
        vp->focus.d[1]=nichi[1];
        vp->focus.d[2]=nichi[2];
        (wdraw.childs+1)->vparam.vflag.ev.deformation=0;        //UJIOKA
        SendMessage((wmenu.childs+2)->hwnd,WM_INITDIALOG,0,0);

        if(ftype==3)
        {
          sprintf(str,"各階変形図 %s",fname);
          SetDlgItemText((wmenu.childs+3)->hwnd,IDV_TITLEBOX,
                           str);

          wparam = MAKEWPARAM(IDV_DEFORMATION,0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,wparam,0);
          if(i==0) SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                               MAKEWPARAM(IDV_DZ,0),0);
          if(i==1) SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                               MAKEWPARAM(IDV_DX,0),0);
          if(i==2) SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                               MAKEWPARAM(IDV_DY,0),0);

          SendMessage(wmain.hwnd,WM_COMMAND,
                      MAKEWPARAM((WORD)IDM_PRINTWAIT,(WORD)0),0);

          wparam = MAKEWPARAM(IDV_DEFORMATION,0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,wparam,0);
          if(i==0) SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                               MAKEWPARAM(IDV_DZ,0),0);
          if(i==1) SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                               MAKEWPARAM(IDV_DX,0),0);
          if(i==2) SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                               MAKEWPARAM(IDV_DY,0),0);

          /*END PAGE*/
          SendMessage(wmain.hwnd,WM_COMMAND,
                      MAKEWPARAM((WORD)IDM_PRINTOUT,(WORD)0),0);
        }

        /*N*/
        vp->focus.d[0]=nichi[0];
        vp->focus.d[1]=nichi[1];
        vp->focus.d[2]=nichi[2];
        SendMessage((wmenu.childs+2)->hwnd,WM_INITDIALOG,0,0);

        if(ftype==1 || ftype==2)
        {
          if(ftype==2)
          {
            sprintf(str,"軸力 Nz %s",fname);
            SetDlgItemText((wmenu.childs+3)->hwnd,IDV_TITLEBOX,
                           str);
          }

          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_NZ,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_NZ_G,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_NZ_B,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_NZ_W,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_NZ_S,0),0);

          SendMessage(wmain.hwnd,WM_COMMAND,
                      MAKEWPARAM((WORD)IDM_PRINTWAIT,(WORD)0),0);

          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_NZ,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_NZ_G,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_NZ_B,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_NZ_W,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_NZ_S,0),0);

          /*END PAGE*/
          if(ftype==2) SendMessage(wmain.hwnd,WM_COMMAND,
                                   MAKEWPARAM((WORD)IDM_PRINTOUT,(WORD)0),0);
        }

        /*Qx*/
        vp->focus.d[0]=qichi[0];
        vp->focus.d[1]=qichi[1];
        vp->focus.d[2]=qichi[2];
        SendMessage((wmenu.childs+2)->hwnd,WM_INITDIALOG,0,0);

        if(ftype==1)
        {
          if(j==0)  SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                                MAKEWPARAM(IDV_QY,0),0);
          if(j==1)  SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                                MAKEWPARAM(IDV_QX,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_QY_G,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_QY_B,0),0);

          SendMessage(wmain.hwnd,WM_COMMAND,
                      MAKEWPARAM((WORD)IDM_PRINTWAIT,(WORD)0),0);

          if(j==0)  SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                                MAKEWPARAM(IDV_QY,0),0);
          if(j==1)  SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                                MAKEWPARAM(IDV_QX,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_QY_G,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_QY_B,0),0);

          /*END PAGE*/
          if(ftype==2) SendMessage(wmain.hwnd,WM_COMMAND,
                                   MAKEWPARAM((WORD)IDM_PRINTOUT,(WORD)0),0);
        }
        if(ftype==2)
        {
          /*Qx*/
          sprintf(str,"せん断力 Qx %s",fname);
          SetDlgItemText((wmenu.childs+3)->hwnd,IDV_TITLEBOX,str);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_QX,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_QX_G,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_QX_B,0),0);

          SendMessage(wmain.hwnd,WM_COMMAND,
                      MAKEWPARAM((WORD)IDM_PRINTWAIT,(WORD)0),0);

          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_QX,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_QX_G,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_QX_B,0),0);

          /*END PAGE*/
          if(ftype==2) SendMessage(wmain.hwnd,WM_COMMAND,
                                   MAKEWPARAM((WORD)IDM_PRINTOUT,(WORD)0),0);

          /*Qy*/
          sprintf(str,"せん断力 Qy %s",fname);
          SetDlgItemText((wmenu.childs+3)->hwnd,IDV_TITLEBOX,str);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_QY,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_QY_G,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_QY_B,0),0);

          SendMessage(wmain.hwnd,WM_COMMAND,
                      MAKEWPARAM((WORD)IDM_PRINTWAIT,(WORD)0),0);

          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_QY,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_QY_G,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_QY_B,0),0);

          /*END PAGE*/
          if(ftype==2) SendMessage(wmain.hwnd,WM_COMMAND,
                                   MAKEWPARAM((WORD)IDM_PRINTOUT,(WORD)0),0);
        }

        /*M*/
        vp->focus.d[0]=michi[0];
        vp->focus.d[1]=michi[1];
        vp->focus.d[2]=michi[2];
        SendMessage((wmenu.childs+2)->hwnd,WM_INITDIALOG,0,0);

        if(ftype==1)
        {
          if(j==0)  SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                                MAKEWPARAM(IDV_MX,0),0);
          if(j==1)  SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                                MAKEWPARAM(IDV_MY,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_MX_G,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_MX_B,0),0);

          SendMessage(wmain.hwnd,WM_COMMAND,
                      MAKEWPARAM((WORD)IDM_PRINTWAIT,(WORD)0),0);

          if(j==0)  SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                                MAKEWPARAM(IDV_MX,0),0);
          if(j==1)  SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                                MAKEWPARAM(IDV_MY,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_MX_G,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_MX_B,0),0);

          /*END PAGE*/
          SendMessage(wmain.hwnd,WM_COMMAND,
                      MAKEWPARAM((WORD)IDM_PRINTOUT,(WORD)0),0);
        }
        if(ftype==2)
        {
          /*Mx*/
          sprintf(str,"曲げモーメント Mx %s",fname);
          SetDlgItemText((wmenu.childs+3)->hwnd,IDV_TITLEBOX,str);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_MX,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_MX_G,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_MX_B,0),0);

          SendMessage(wmain.hwnd,WM_COMMAND,
                      MAKEWPARAM((WORD)IDM_PRINTWAIT,(WORD)0),0);

          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_MX,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_MX_G,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_MX_B,0),0);

          /*END PAGE*/
          SendMessage(wmain.hwnd,WM_COMMAND,
                      MAKEWPARAM((WORD)IDM_PRINTOUT,(WORD)0),0);

          /*My*/
          sprintf(str,"曲げモーメント My %s",fname);
          SetDlgItemText((wmenu.childs+3)->hwnd,IDV_TITLEBOX,str);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_MX,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_MX_G,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_MX_B,0),0);

          SendMessage(wmain.hwnd,WM_COMMAND,
                      MAKEWPARAM((WORD)IDM_PRINTWAIT,(WORD)0),0);

          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_MX,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_MX_G,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_MX_B,0),0);

          /*END PAGE*/
          SendMessage(wmain.hwnd,WM_COMMAND,
                      MAKEWPARAM((WORD)IDM_PRINTOUT,(WORD)0),0);
        }

        freestr(data,nstr);
      }
    }
  }

  /*PRINT END.*/
  wparam = MAKEWPARAM((WORD)IDM_PRINTEND,(WORD)0);
  SendMessage(wmain.hwnd,WM_COMMAND,wparam,0);

  return;
}/*printarclmfigures*/

void printarclmfigures2(FILE *fin,struct viewparam *vp,int ftype,int mode)
/*TYPE 1:PRINT STRESS BY AXIS.*/
/*TYPE 2:PRINT STRESS BY FLOOR.*/
/*TYPE 3:PRINT DISPLACEMENT BY FLOOR.*/
/*SETTING BEFORE BEGINNING:VISIBLE ELEMENT TYPE*/
/*                         M FACTOR*/
/*                         H SIZE*/
{
  WPARAM wparam;
  char **data,str[256];
  int i,j,k,np1,np2,nstr,flag/*,ftype*/;

  char fname[80];
  int nframe;
  double nichi[3],qichi[3],michi[3];

  double initd[3];
  double inittheta,initphi;
  double initr,initodv;

  /*ftype=1;*/ /*FTYPE:1=STRESS BY AXIS, 2=STRESS BY FLOOR, 3=DISP BY FLOOR*/

  getviewparam((wmenu.childs+2)->hwnd,
              &((wdraw.childs+1)->vparam));

  initd[0] =vp->focus.d[0];
  initd[1] =vp->focus.d[1];
  initd[2] =vp->focus.d[2];
  inittheta=vp->theta;
  initphi  =vp->phi;
  initr    =vp->r;
  initodv  =vp->odv;

  /*PRINTER SETUP*/
  SendMessage(wmain.hwnd,WM_COMMAND,
              MAKEWPARAM((WORD)IDM_PRINTSETUP,(WORD)0),0);

  /*VISIBLE SETTING*/
  wparam = MAKEWPARAM(IDV_GLOBALAXIS,0);
  SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,wparam,0);
  wparam = MAKEWPARAM(IDV_NODECODE,0);
/*mihara for zahyo 20190329///////////////////////////////////////////////////*/
	wparam = MAKEWPARAM(IDV_GX, 0);
	SendMessage((wmenu.childs + 3)->hwnd, WM_COMMAND, wparam, 0);
	wparam = MAKEWPARAM(IDV_GY, 0);
	SendMessage((wmenu.childs + 3)->hwnd, WM_COMMAND, wparam, 0);
	wparam = MAKEWPARAM(IDV_GZ, 0);
	SendMessage((wmenu.childs + 3)->hwnd, WM_COMMAND, wparam, 0);
/*mihara for zahyo 20190329///////////////////////////////////////////////////*/
  SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,wparam,0);
    (wdraw.childs+1)->vparam.vflag.ev.srcanrate  =0;
    (wdraw.childs+1)->vparam.vflag.ev.srcancolor =0;

  if(ftype==1)
  {
    wparam = MAKEWPARAM(IDV_SLAB,0);
    SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,wparam,0);
    wparam = MAKEWPARAM(IDV_PAGETITLE,0);
    SendMessage((wmenu.childs+4)->hwnd,WM_COMMAND,wparam,0);

    /*UJIOKA*/
    wparam = MAKEWPARAM(IDV_DEFORMATION,0);
    SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,wparam,0);
    wparam = MAKEWPARAM(IDV_TITLE,0);
    SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,wparam,0);
    (wdraw.childs+1)->vparam.vflag.nv.code=0;
    (wdraw.childs+1)->vparam.vflag.ev.sectioncode=0;

    (wdraw.childs+1)->vparam.chapter=7;      /*CHAPTER*/
    (wdraw.childs+1)->vparam.section=2-1;    /*INITIAL SECTION*/
    /*UJIOKA*/
  }
  if(ftype==2)
  {
    wparam = MAKEWPARAM(IDV_TITLE,0);
    SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,wparam,0);
    wparam = MAKEWPARAM(IDV_VIEW,0);
    SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,wparam,0);
    wparam = MAKEWPARAM(IDV_INPUTFILE,0);
    SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,wparam,0);
    wparam = MAKEWPARAM(IDV_OUTPUTFILE,0);
    SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,wparam,0);
  }
  if(ftype==3)
  {
    wparam = MAKEWPARAM(IDV_DEFORMATION,0);
    SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,wparam,0);
    wparam = MAKEWPARAM(IDV_TITLE,0);
    SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,wparam,0);
    wparam = MAKEWPARAM(IDV_VIEW,0);
    SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,wparam,0);
    wparam = MAKEWPARAM(IDV_INPUTFILE,0);
    SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,wparam,0);
    wparam = MAKEWPARAM(IDV_OUTPUTFILE,0);
    SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,wparam,0);
  }

  for(i=0;i<3;i++) /*FOR PERIOD Z,X,Y*/
  {
    fseek(fin,0L,SEEK_SET);

    (wdraw.childs+1)->vparam.section++;
    (wdraw.childs+1)->vparam.subsection=0;

    if(i==0)
    {
      if(ftype==1){np1=0; np2=1;}
      if(ftype==2 || ftype==3){np1=0; np2=0;}

      GetDlgItemText((wmenu.childs+2)->hwnd,ID_INPUTFILEZ,str,256);
      sprintf((wdraw.childs+1)->inpfile,str);
      SetDlgItemText((wmenu.childs+2)->hwnd,ID_INPUTFILE,
                     (wdraw.childs+1)->inpfile);
      GetDlgItemText((wmenu.childs+2)->hwnd,ID_OUTPUTFILEZ,str,256);
      sprintf((wdraw.childs+1)->otpfile,str);
      SetDlgItemText((wmenu.childs+2)->hwnd,ID_OUTPUTFILE,
                     (wdraw.childs+1)->otpfile);
    }
    if(i==1)
    {
//      np1=0; np2=0;
      np1=1; np2=1;

      GetDlgItemText((wmenu.childs+2)->hwnd,ID_INPUTFILEX,str,256);
      sprintf((wdraw.childs+1)->inpfile,str);
      SetDlgItemText((wmenu.childs+2)->hwnd,ID_INPUTFILE,
                     (wdraw.childs+1)->inpfile);
      GetDlgItemText((wmenu.childs+2)->hwnd,ID_OUTPUTFILEX,str,256);
      sprintf((wdraw.childs+1)->otpfile,str);
      SetDlgItemText((wmenu.childs+2)->hwnd,ID_OUTPUTFILE,
                     (wdraw.childs+1)->otpfile);

    }
    if(i==2)
    {
//      np1=1; np2=1;
      np1=0; np2=0;

      GetDlgItemText((wmenu.childs+2)->hwnd,ID_INPUTFILEY,str,256);
      sprintf((wdraw.childs+1)->inpfile,str);
      SetDlgItemText((wmenu.childs+2)->hwnd,ID_INPUTFILE,
                     (wdraw.childs+1)->inpfile);
      GetDlgItemText((wmenu.childs+2)->hwnd,ID_OUTPUTFILEY,str,256);
      sprintf((wdraw.childs+1)->otpfile,str);
      SetDlgItemText((wmenu.childs+2)->hwnd,ID_OUTPUTFILE,
                     (wdraw.childs+1)->otpfile);
    }

    /*OPEN FILE AND DRAW FRAME.*/
    wparam = MAKEWPARAM(IDD_VIEW,0);
    SendMessage((wmenu.childs+2)->hwnd,WM_COMMAND,wparam,0);
    globalmessageflag=0;
    wparam = MAKEWPARAM(IDD_OPENRESULT,0);
    SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,wparam,0);
    wparam = MAKEWPARAM(IDV_UNITNM,0);
    SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,wparam,0);

    /*PRINT DEFORMATION*/
    if(ftype==1)
    {
       (wdraw.childs+1)->vparam.vflag.mv.inputfile =1;
       (wdraw.childs+1)->vparam.vflag.mv.outputfile=1;
       (wdraw.childs+1)->vparam.vflag.mv.view =1;
       (wdraw.childs+1)->vparam.vflag.axis=1;

/*
       setviewpoint((wdraw.childs+0)->hwnd,arc,
                   &((wdraw.childs+1)->vparam));
       setviewparam((wmenu.childs+2)->hwnd,
                    (wdraw.childs+1)->vparam);
       SendMessage((wmenu.childs+2)->hwnd,WM_INITDIALOG,0,0);
*/
       vp->range.min.d[GX]=-100.0;
       vp->range.min.d[GY]=-100.0;
       vp->range.min.d[GZ]=-100.0;
       vp->range.max.d[GX]= 100.0;
       vp->range.max.d[GY]= 100.0;
       vp->range.max.d[GZ]= 100.0;
       vp->focus.d[0]     =initd[0];
       vp->focus.d[1]     =initd[1];
       vp->focus.d[2]     =initd[2];
       vp->theta          =inittheta;
       vp->phi            =initphi;
       vp->r              =initr;
       vp->odv            =initodv;
       (wdraw.childs+1)->vparam.vflag.ev.etype[COLUMN]=1;
       (wdraw.childs+1)->vparam.vflag.ev.etype[GIRDER]=1;
       (wdraw.childs+1)->vparam.vflag.ev.etype[BEAM]=1;
       (wdraw.childs+1)->vparam.vflag.ev.etype[BRACE]=1;
       (wdraw.childs+1)->vparam.vflag.ev.etype[WALL]=1;
       (wdraw.childs+1)->vparam.vflag.ev.etype[SLAB]=1;
       SendMessage((wmenu.childs+2)->hwnd,WM_INITDIALOG,0,0);

       clearwindow(*(wdraw.childs+1));
       drawarclmframe((wdraw.childs+1)->hdcC,
                      (wdraw.childs+1)->vparam,arc,0,ONSCREEN);
       SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);

       if(i==0)
       {
          sprintf(str,"%d.%d：長期解析結果",
                 (wdraw.childs+1)->vparam.chapter,
                 (wdraw.childs+1)->vparam.section);
          (wdraw.childs+1)->vparam.vflag.nv.disps[2]=1;
       }
       if(i==1)
       {
          sprintf(str,"%d.%d：X方向地震時解析結果",
                 (wdraw.childs+1)->vparam.chapter,
                 (wdraw.childs+1)->vparam.section);
          (wdraw.childs+1)->vparam.vflag.nv.disps[0]=1;
       }
       if(i==2)
       {
          sprintf(str,"%d.%d：Y方向地震時解析結果",
                 (wdraw.childs+1)->vparam.chapter,
                 (wdraw.childs+1)->vparam.section);
         (wdraw.childs+1)->vparam.vflag.nv.disps[1]=1;
       }
       SetDlgItemText((wmenu.childs+4)->hwnd,IDV_PAGETITLEBOX,str);
       (wdraw.childs+1)->vparam.subsection++;
       SendMessage(wmain.hwnd,WM_COMMAND,
                   MAKEWPARAM((WORD)IDM_PRINTWAIT,(WORD)0),0);
       /*END PAGE*/
       SendMessage(wmain.hwnd,WM_COMMAND,
                   MAKEWPARAM((WORD)IDM_PRINTOUT,(WORD)0),0);

       (wdraw.childs+1)->vparam.vflag.mv.inputfile =0;
       (wdraw.childs+1)->vparam.vflag.mv.outputfile=0;
       (wdraw.childs+1)->vparam.vflag.mv.view =0;
       (wdraw.childs+1)->vparam.vflag.axis=0;
    }

    for(j=np1;j<=np2;j++) /*FOR FRAME X,Y*/
    {
      fseek(fin,0L,SEEK_SET);
//      if(j==0)
      if(j==1)
      {
        flag=0;
        while(!flag)
        {
		  data=fgetsbrk(fin,&nstr);
          if(!strcmp(*(data+0),"XFRAME")) flag=1;
        }
      }
//      if(j==1)
      if(j==0)
      {
        flag=0;
        while(!flag)
        {
		  data=fgetsbrk(fin,&nstr);
          if(!strcmp(*(data+0),"YFRAME")) flag=1;
        }
      }

      data=fgetsbrk(fin,&nstr);
      if(nstr==0) return;
      if(!strcmp(*(data+0),"ANGLE"))
      {
        vp->theta=strtod(*(data+1),NULL);
        vp->phi  =strtod(*(data+2),NULL);
      }
      freestr(data,nstr);

      data=fgetsbrk(fin,&nstr);
      if(nstr==0) return;
      if(!strcmp(*(data+0),"DISTS"))
      {
        vp->r  =strtod(*(data+1),NULL);
        vp->odv=strtod(*(data+2),NULL);
      }
      freestr(data,nstr);

	  data=fgetsbrk(fin,&nstr);
      if(nstr==0) return;
      if(!strcmp(*(data+0),"NICHI"))
      {
        nichi[0]=strtod(*(data+1),NULL);
        nichi[1]=strtod(*(data+2),NULL);
        nichi[2]=strtod(*(data+3),NULL);
      }
      freestr(data,nstr);

      data=fgetsbrk(fin,&nstr);
      if(nstr==0) return;
      if(!strcmp(*(data+0),"QICHI"))
      {
        qichi[0]=strtod(*(data+1),NULL);
        qichi[1]=strtod(*(data+2),NULL);
        qichi[2]=strtod(*(data+3),NULL);
      }
      freestr(data,nstr);

      data=fgetsbrk(fin,&nstr);
      if(nstr==0) return;
      if(!strcmp(*(data+0),"MICHI"))
      {
        michi[0]=strtod(*(data+1),NULL);
        michi[1]=strtod(*(data+2),NULL);
        michi[2]=strtod(*(data+3),NULL);
      }
      freestr(data,nstr);

      data=fgetsbrk(fin,&nstr);
      if(nstr==0) return;
      if(!strcmp(*(data+0),"NFRAME"))
      {
        nframe=(int)strtol(*(data+1),NULL,10);
      }
      freestr(data,nstr);

      for(k=0;k<nframe;k++)
      {
        (wdraw.childs+1)->vparam.subsection++;
        (wdraw.childs+1)->vparam.vflag.nv.disps[0]=0;
        (wdraw.childs+1)->vparam.vflag.nv.disps[1]=0;
        (wdraw.childs+1)->vparam.vflag.nv.disps[2]=0;

        data=fgetsbrk(fin,&nstr);
        if(nstr==0) return;

        if(!strcmp(*(data+0),"FRAME"))
        {
          if(ftype==1)
          {
            sprintf(fname,"%s 通り応力図",*(data+1));
            SetDlgItemText((wmenu.childs+4)->hwnd,IDV_PAGETITLEBOX,
                           fname);
          }
          if(ftype==2 || ftype==3)
          {
            sprintf(fname,"%s",*(data+1));
          }

        }
        if(!strcmp(*(data+2),"HEAD"))
        {
          vp->range.min.d[GX]=strtod(*(data+3),NULL);
          vp->range.min.d[GY]=strtod(*(data+4),NULL);
          vp->range.min.d[GZ]=strtod(*(data+5),NULL);
        }
        if(!strcmp(*(data+6),"TAIL"))
        {
          vp->range.max.d[GX]=strtod(*(data+7),NULL);
          vp->range.max.d[GY]=strtod(*(data+8),NULL);
          vp->range.max.d[GZ]=strtod(*(data+9),NULL);
        }

        SendMessage((wmenu.childs+2)->hwnd,WM_INITDIALOG,0,0);

        /*DISPLACEMENTN*/
        vp->focus.d[0]=nichi[0];
        vp->focus.d[1]=nichi[1];
        vp->focus.d[2]=nichi[2];
        (wdraw.childs+1)->vparam.vflag.ev.deformation=0;        //UJIOKA
        SendMessage((wmenu.childs+2)->hwnd,WM_INITDIALOG,0,0);

        if(ftype==3)
        {
          sprintf(str,"各階変形図 %s",fname);
          SetDlgItemText((wmenu.childs+3)->hwnd,IDV_TITLEBOX,
                           str);

          wparam = MAKEWPARAM(IDV_DEFORMATION,0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,wparam,0);
          if(i==0) SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                               MAKEWPARAM(IDV_DZ,0),0);
          if(i==1) SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                               MAKEWPARAM(IDV_DX,0),0);
          if(i==2) SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                               MAKEWPARAM(IDV_DY,0),0);

          SendMessage(wmain.hwnd,WM_COMMAND,
                      MAKEWPARAM((WORD)IDM_PRINTWAIT,(WORD)0),0);

          wparam = MAKEWPARAM(IDV_DEFORMATION,0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,wparam,0);
          if(i==0) SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                               MAKEWPARAM(IDV_DZ,0),0);
          if(i==1) SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                               MAKEWPARAM(IDV_DX,0),0);
          if(i==2) SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                               MAKEWPARAM(IDV_DY,0),0);

          /*END PAGE*/
          SendMessage(wmain.hwnd,WM_COMMAND,
                      MAKEWPARAM((WORD)IDM_PRINTOUT,(WORD)0),0);
        }

        /*N*/
        vp->focus.d[0]=nichi[0];
        vp->focus.d[1]=nichi[1];
        vp->focus.d[2]=nichi[2];
        SendMessage((wmenu.childs+2)->hwnd,WM_INITDIALOG,0,0);

        if(ftype==1 || ftype==2)
        {
          if(ftype==2)
          {
            sprintf(str,"軸力 Nz %s",fname);
            SetDlgItemText((wmenu.childs+3)->hwnd,IDV_TITLEBOX,
                           str);
          }

          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_NZ,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_NZ_G,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_NZ_B,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_NZ_W,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_NZ_S,0),0);

          SendMessage(wmain.hwnd,WM_COMMAND,
                      MAKEWPARAM((WORD)IDM_PRINTWAIT,(WORD)0),0);

          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_NZ,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_NZ_G,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_NZ_B,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_NZ_W,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_NZ_S,0),0);

          /*END PAGE*/
          if(ftype==2) SendMessage(wmain.hwnd,WM_COMMAND,
                                   MAKEWPARAM((WORD)IDM_PRINTOUT,(WORD)0),0);
        }

        /*Qx*/
        vp->focus.d[0]=qichi[0];
        vp->focus.d[1]=qichi[1];
        vp->focus.d[2]=qichi[2];
        SendMessage((wmenu.childs+2)->hwnd,WM_INITDIALOG,0,0);

        if(ftype==1)
        {
          if(j==0)  SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                                MAKEWPARAM(IDV_QY,0),0);
          if(j==1)  SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                                MAKEWPARAM(IDV_QX,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_QY_G,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_QY_B,0),0);

          SendMessage(wmain.hwnd,WM_COMMAND,
                      MAKEWPARAM((WORD)IDM_PRINTWAIT,(WORD)0),0);

          if(j==0)  SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                                MAKEWPARAM(IDV_QY,0),0);
          if(j==1)  SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                                MAKEWPARAM(IDV_QX,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_QY_G,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_QY_B,0),0);

          /*END PAGE*/
          if(ftype==2) SendMessage(wmain.hwnd,WM_COMMAND,
                                   MAKEWPARAM((WORD)IDM_PRINTOUT,(WORD)0),0);
        }
        if(ftype==2)
        {
          /*Qx*/
          sprintf(str,"せん断力 Qx %s",fname);
          SetDlgItemText((wmenu.childs+3)->hwnd,IDV_TITLEBOX,str);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_QX,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_QX_G,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_QX_B,0),0);

          SendMessage(wmain.hwnd,WM_COMMAND,
                      MAKEWPARAM((WORD)IDM_PRINTWAIT,(WORD)0),0);

          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_QX,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_QX_G,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_QX_B,0),0);

          /*END PAGE*/
          if(ftype==2) SendMessage(wmain.hwnd,WM_COMMAND,
                                   MAKEWPARAM((WORD)IDM_PRINTOUT,(WORD)0),0);

          /*Qy*/
          sprintf(str,"せん断力 Qy %s",fname);
          SetDlgItemText((wmenu.childs+3)->hwnd,IDV_TITLEBOX,str);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_QY,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_QY_G,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_QY_B,0),0);

          SendMessage(wmain.hwnd,WM_COMMAND,
                      MAKEWPARAM((WORD)IDM_PRINTWAIT,(WORD)0),0);

          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_QY,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_QY_G,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_QY_B,0),0);

          /*END PAGE*/
          if(ftype==2) SendMessage(wmain.hwnd,WM_COMMAND,
                                   MAKEWPARAM((WORD)IDM_PRINTOUT,(WORD)0),0);
        }

        /*M*/
        vp->focus.d[0]=michi[0];
        vp->focus.d[1]=michi[1];
        vp->focus.d[2]=michi[2];
        SendMessage((wmenu.childs+2)->hwnd,WM_INITDIALOG,0,0);

        if(ftype==1)
        {
          if(j==0)  SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                                MAKEWPARAM(IDV_MX,0),0);
          if(j==1)  SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                                MAKEWPARAM(IDV_MY,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_MX_G,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_MX_B,0),0);

          SendMessage(wmain.hwnd,WM_COMMAND,
                      MAKEWPARAM((WORD)IDM_PRINTWAIT,(WORD)0),0);

          if(j==0)  SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                                MAKEWPARAM(IDV_MX,0),0);
          if(j==1)  SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                                MAKEWPARAM(IDV_MY,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_MX_G,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_MX_B,0),0);

          /*END PAGE*/
          SendMessage(wmain.hwnd,WM_COMMAND,
                      MAKEWPARAM((WORD)IDM_PRINTOUT,(WORD)0),0);
        }
        if(ftype==2)
        {
          /*Mx*/
          sprintf(str,"曲げモーメント Mx %s",fname);
          SetDlgItemText((wmenu.childs+3)->hwnd,IDV_TITLEBOX,str);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_MX,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_MX_G,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_MX_B,0),0);

          SendMessage(wmain.hwnd,WM_COMMAND,
                      MAKEWPARAM((WORD)IDM_PRINTWAIT,(WORD)0),0);

          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_MX,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_MX_G,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_MX_B,0),0);

          /*END PAGE*/
          SendMessage(wmain.hwnd,WM_COMMAND,
                      MAKEWPARAM((WORD)IDM_PRINTOUT,(WORD)0),0);

          /*My*/
          sprintf(str,"曲げモーメント My %s",fname);
          SetDlgItemText((wmenu.childs+3)->hwnd,IDV_TITLEBOX,str);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_MX,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_MX_G,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_MX_B,0),0);

          SendMessage(wmain.hwnd,WM_COMMAND,
                      MAKEWPARAM((WORD)IDM_PRINTWAIT,(WORD)0),0);

          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_MX,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_MX_G,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_MX_B,0),0);

          /*END PAGE*/
          SendMessage(wmain.hwnd,WM_COMMAND,
                      MAKEWPARAM((WORD)IDM_PRINTOUT,(WORD)0),0);
        }

        freestr(data,nstr);
      }
    }
  }


  if(!mode)
  {
  /*PRINT END.*/
  wparam = MAKEWPARAM((WORD)IDM_PRINTEND,(WORD)0);
  SendMessage(wmain.hwnd,WM_COMMAND,wparam,0);
  }

  vp->range.min.d[GX]=-100.0;
  vp->range.min.d[GY]=-100.0;
  vp->range.min.d[GZ]=-100.0;
  vp->range.max.d[GX]= 100.0;
  vp->range.max.d[GY]= 100.0;
  vp->range.max.d[GZ]= 100.0;
  vp->focus.d[0]     =initd[0];
  vp->focus.d[1]     =initd[1];
  vp->focus.d[2]     =initd[2];
  vp->theta          =inittheta;
  vp->phi            =initphi;
  vp->r              =initr;
  vp->odv            =initodv;
  SendMessage((wmenu.childs+2)->hwnd,WM_INITDIALOG,0,0);

  clearwindow(*(wdraw.childs+1));
  drawarclmframe((wdraw.childs+1)->hdcC,
                 (wdraw.childs+1)->vparam,arc,0,ONSCREEN);
  SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);

  return;
}/*printarclmfigures2*/

void printsrcanfigures(FILE *fin,struct viewparam *vp,int ftype,int mode)
/*TYPE 1:PRINT SRCAN BY AXIS.*/
/*TYPE 2:PRINT SRCAN BY FLOOR.*/
/*TYPE 3:PRINT SRCAN BY FLOOR.*/
/*SETTING BEFORE BEGINNING:VISIBLE ELEMENT TYPE*/
/*                         M FACTOR*/
/*                         H SIZE*/
{
  WPARAM wparam;
  char **data,str[256],s[80];
  int i,j,jj,k,n,np1,np2,nstr,flag/*,ftype*/;

  char fname[80];
  int nframe;
  double nichi[3],qichi[3],michi[3];

  double initd[3];
  double inittheta,initphi;
  double initr,initodv;

  /*ftype=1;*/ /*FTYPE:1=SRCAN BY AXIS, 2=SRCAN BY FLOOR, 3=SRCAN BY FLOOR*/

  getviewparam((wmenu.childs+2)->hwnd,
              &((wdraw.childs+1)->vparam));

  initd[0] =vp->focus.d[0];
  initd[1] =vp->focus.d[1];
  initd[2] =vp->focus.d[2];
  inittheta=vp->theta;
  initphi  =vp->phi;
  initr    =vp->r;
  initodv  =vp->odv;

  /*PRINTER SETUP*/
  SendMessage(wmain.hwnd,WM_COMMAND,
              MAKEWPARAM((WORD)IDM_PRINTSETUP,(WORD)0),0);

  /*VISIBLE SETTING*/
  wparam = MAKEWPARAM(IDV_GLOBALAXIS,0);
  SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,wparam,0);
  wparam = MAKEWPARAM(IDV_NODECODE,0);
/*mihara for zahyo 20190329///////////////////////////////////////////////////*/
	wparam = MAKEWPARAM(IDV_GX, 0);
	SendMessage((wmenu.childs + 3)->hwnd, WM_COMMAND, wparam, 0);
	wparam = MAKEWPARAM(IDV_GY, 0);
	SendMessage((wmenu.childs + 3)->hwnd, WM_COMMAND, wparam, 0);
	wparam = MAKEWPARAM(IDV_GZ, 0);
	SendMessage((wmenu.childs + 3)->hwnd, WM_COMMAND, wparam, 0);
/*mihara for zahyo 20190329///////////////////////////////////////////////////*/
  SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,wparam,0);
  if(ftype==1)
  {
    wparam = MAKEWPARAM(IDV_SLAB,0);
    SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,wparam,0);
    wparam = MAKEWPARAM(IDV_PAGETITLE,0);
    SendMessage((wmenu.childs+4)->hwnd,WM_COMMAND,wparam,0);

    /*UJIOKA*/
    wparam = MAKEWPARAM(IDV_DEFORMATION,0);
    SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,wparam,0);
    wparam = MAKEWPARAM(IDV_TITLE,0);
    SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,wparam,0);

    (wdraw.childs+1)->vparam.chapter=4;      /*CHAPTER*/
    (wdraw.childs+1)->vparam.section=4-1;    /*INITIAL SECTION*/
    /*UJIOKA*/
  }
  if(ftype==2)
  {
    wparam = MAKEWPARAM(IDV_TITLE,0);
    SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,wparam,0);
    wparam = MAKEWPARAM(IDV_VIEW,0);
    SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,wparam,0);
    wparam = MAKEWPARAM(IDV_INPUTFILE,0);
    SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,wparam,0);
    wparam = MAKEWPARAM(IDV_OUTPUTFILE,0);
    SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,wparam,0);

    (wdraw.childs+1)->vparam.chapter=4;      /*CHAPTER*/
    (wdraw.childs+1)->vparam.section=4-1;    /*INITIAL SECTION*/
  }
  if(ftype==3)
  {
    wparam = MAKEWPARAM(IDV_DEFORMATION,0);
    SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,wparam,0);
    wparam = MAKEWPARAM(IDV_TITLE,0);
    SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,wparam,0);
    wparam = MAKEWPARAM(IDV_VIEW,0);
    SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,wparam,0);
    wparam = MAKEWPARAM(IDV_INPUTFILE,0);
    SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,wparam,0);
    wparam = MAKEWPARAM(IDV_OUTPUTFILE,0);
    SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,wparam,0);
    (wdraw.childs+1)->vparam.vflag.nv.code=0;
    (wdraw.childs+1)->vparam.vflag.ev.sectioncode=0;

    (wdraw.childs+1)->vparam.chapter=4;      /*CHAPTER*/
    (wdraw.childs+1)->vparam.section=4-1;    /*INITIAL SECTION*/
  }

//  for(i=0;i<3;i++) /*FOR PERIOD Z,X,Y*/
  {
    fseek(fin,0L,SEEK_SET);

    (wdraw.childs+1)->vparam.section++;
    (wdraw.childs+1)->vparam.subsection=0;

    /*OPEN FILE AND DRAW FRAME.*/
    wparam = MAKEWPARAM(IDD_VIEW,0);
    SendMessage((wmenu.childs+2)->hwnd,WM_COMMAND,wparam,0);
    globalmessageflag=0;

    n=strcspn((wdraw.childs+1)->otpfile,".");
    strncpy(str,(wdraw.childs+1)->otpfile,n);
    str[n]='\0';
    strcpy((wdraw.childs+1)->otpfile,str);
    strcat((wdraw.childs+1)->otpfile,".rat");
    SetDlgItemText((wmenu.childs+2)->hwnd,ID_OUTPUTFILE,
	    			 (wdraw.childs+1)->otpfile);
    fout=fopen((wdraw.childs+1)->otpfile,"r");
//    fout=fgetstofopen("\0","r",ID_OUTPUTFILE);
    readsrcanrate(fout,&arc);

    (wdraw.childs+1)->vparam.vflag.ev.srcanrate  =1;
    (wdraw.childs+1)->vparam.vflag.ev.srcancolor =1;
    wparam = MAKEWPARAM(IDV_UNITNM,0);
    SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,wparam,0);

    /*PRINT*/
//    if(ftype==1)
    {
       if(ftype==1){np1=0; np2=1;}
       if(ftype==2 || ftype==3){np1=0; np2=0;}
       (wdraw.childs+1)->vparam.vflag.mv.pagetitle  =1;
       (wdraw.childs+1)->vparam.vflag.mv.title  =1;
       (wdraw.childs+1)->vparam.vflag.mv.inputfile  =1;
       (wdraw.childs+1)->vparam.vflag.mv.outputfile =1;
       (wdraw.childs+1)->vparam.vflag.mv.view       =1;

/*
       setviewpoint((wdraw.childs+0)->hwnd,arc,
                   &((wdraw.childs+1)->vparam));
       setviewparam((wmenu.childs+2)->hwnd,
                    (wdraw.childs+1)->vparam);
       SendMessage((wmenu.childs+2)->hwnd,WM_INITDIALOG,0,0);
*/
       vp->range.min.d[GX]=-100.0;
       vp->range.min.d[GY]=-100.0;
       vp->range.min.d[GZ]=-100.0;
       vp->range.max.d[GX]= 100.0;
       vp->range.max.d[GY]= 100.0;
       vp->range.max.d[GZ]= 100.0;
       vp->focus.d[0]     =initd[0];
       vp->focus.d[1]     =initd[1];
       vp->focus.d[2]     =initd[2];
       vp->theta          =inittheta;
       vp->phi            =initphi;
       vp->r              =initr;
       vp->odv            =initodv;
       (wdraw.childs+1)->vparam.vflag.ev.etype[COLUMN]=1;
       (wdraw.childs+1)->vparam.vflag.ev.etype[GIRDER]=1;
       (wdraw.childs+1)->vparam.vflag.ev.etype[BEAM]=1;
       (wdraw.childs+1)->vparam.vflag.ev.etype[BRACE]=1;
       (wdraw.childs+1)->vparam.vflag.ev.etype[WALL]=1;
       (wdraw.childs+1)->vparam.vflag.ev.etype[SLAB]=1;
       SendMessage((wmenu.childs+2)->hwnd,WM_INITDIALOG,0,0);

       clearwindow(*(wdraw.childs+1));
       drawarclmframe((wdraw.childs+1)->hdcC,
                      (wdraw.childs+1)->vparam,arc,0,ONSCREEN);
       SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);
       sprintf(str,"%d.%d：断面検定比図",
                 (wdraw.childs+1)->vparam.chapter,
                 (wdraw.childs+1)->vparam.section);
       SetDlgItemText((wmenu.childs+4)->hwnd,IDV_PAGETITLEBOX,str);
       sprintf(str,"解析モデル全体表示");
       SetDlgItemText((wmenu.childs+3)->hwnd,IDV_TITLEBOX,str);
       (wdraw.childs+1)->vparam.subsection++;
       SendMessage(wmain.hwnd,WM_COMMAND,
                   MAKEWPARAM((WORD)IDM_PRINTWAIT,(WORD)0),0);
       /*END PAGE*/
       SendMessage(wmain.hwnd,WM_COMMAND,
                   MAKEWPARAM((WORD)IDM_PRINTOUT,(WORD)0),0);


       (wdraw.childs+1)->vparam.vflag.mv.pagetitle  =0;

       clearwindow(*(wdraw.childs+1));
       drawarclmframe((wdraw.childs+1)->hdcC,
                      (wdraw.childs+1)->vparam,arc,0,ONSCREEN);
       SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);

       (wdraw.childs+1)->vparam.vflag.ev.etype[BRACE] =0;
       (wdraw.childs+1)->vparam.vflag.ev.etype[WALL]  =0;
       (wdraw.childs+1)->vparam.vflag.ev.etype[SLAB]  =0;
       sprintf(str,"柱梁のみ表示");
       SetDlgItemText((wmenu.childs+3)->hwnd,IDV_TITLEBOX,str);
       (wdraw.childs+1)->vparam.subsection++;
       SendMessage(wmain.hwnd,WM_COMMAND,
                   MAKEWPARAM((WORD)IDM_PRINTWAIT,(WORD)0),0);
       /*END PAGE*/
       SendMessage(wmain.hwnd,WM_COMMAND,
                   MAKEWPARAM((WORD)IDM_PRINTOUT,(WORD)0),0);

       (wdraw.childs+1)->vparam.vflag.ev.etype[GIRDER]=0;
       (wdraw.childs+1)->vparam.vflag.ev.etype[BEAM]  =0;
       sprintf(str,"柱のみ表示");
       SetDlgItemText((wmenu.childs+3)->hwnd,IDV_TITLEBOX,str);
       (wdraw.childs+1)->vparam.subsection++;
       SendMessage(wmain.hwnd,WM_COMMAND,
                   MAKEWPARAM((WORD)IDM_PRINTWAIT,(WORD)0),0);
       /*END PAGE*/
       SendMessage(wmain.hwnd,WM_COMMAND,
                   MAKEWPARAM((WORD)IDM_PRINTOUT,(WORD)0),0);
    }

    for(jj=0;jj<4;jj++)
    {
     if(ftype==1)
     {
       (wdraw.childs+1)->vparam.vflag.ev.etype[COLUMN]=1;
       (wdraw.childs+1)->vparam.vflag.ev.etype[GIRDER]=1;
       (wdraw.childs+1)->vparam.vflag.ev.etype[BEAM]=1;
       (wdraw.childs+1)->vparam.vflag.ev.etype[BRACE]=1;
       (wdraw.childs+1)->vparam.vflag.ev.etype[WALL]=1;
       (wdraw.childs+1)->vparam.vflag.ev.etype[SLAB]=1;
       jj=4;
     }

     if(ftype==2 || ftype==3)
     {
       (wdraw.childs+1)->vparam.vflag.ev.etype[COLUMN]=0;
       (wdraw.childs+1)->vparam.vflag.ev.etype[GIRDER]=0;
       (wdraw.childs+1)->vparam.vflag.ev.etype[BEAM]=0;
       (wdraw.childs+1)->vparam.vflag.ev.etype[BRACE]=0;
       (wdraw.childs+1)->vparam.vflag.ev.etype[WALL]=0;
       (wdraw.childs+1)->vparam.vflag.ev.etype[SLAB]=0;


       if(jj==0)
       {
         (wdraw.childs+1)->vparam.vflag.ev.etype[COLUMN]=1;
         sprintf(s,"柱のみ表示");
       }
       if(jj==1)
       {
         (wdraw.childs+1)->vparam.vflag.ev.etype[GIRDER]=1;
         (wdraw.childs+1)->vparam.vflag.ev.etype[BEAM]=1;
         sprintf(s,"梁のみ表示");
       }
       if(jj==2)
       {
         (wdraw.childs+1)->vparam.vflag.ev.etype[SLAB]=1;
         sprintf(s,"床のみ表示");
       }
       if(jj==3)
       {
         (wdraw.childs+1)->vparam.vflag.ev.etype[WALL]=1;
         sprintf(s,"壁のみ表示");
         break;
       }

/*
       if(jj==2)
       {
         (wdraw.childs+1)->vparam.vflag.ev.etype[WALL]=1;
         sprintf(s,"壁のみ表示");
       }
*/
     }
     for(j=np1;j<=np2;j++) /*FOR FRAME X,Y*/
     {
      fseek(fin,0L,SEEK_SET);
      {
        flag=0;
        while(!flag)
        {
          data=fgetsbrk(fin,&nstr);
          if(!strcmp(*(data+0),"ZFRAME")) flag=1;
        }
      }

      data=fgetsbrk(fin,&nstr);
      if(nstr==0) return;
      if(!strcmp(*(data+0),"ANGLE"))
      {
        vp->theta=strtod(*(data+1),NULL);
        vp->phi  =strtod(*(data+2),NULL);
      }
      freestr(data,nstr);

      data=fgetsbrk(fin,&nstr);
      if(nstr==0) return;
      if(!strcmp(*(data+0),"DISTS"))
      {
        vp->r  =strtod(*(data+1),NULL);
        vp->odv=strtod(*(data+2),NULL);
      }
      freestr(data,nstr);

	  data=fgetsbrk(fin,&nstr);
      if(nstr==0) return;
      if(!strcmp(*(data+0),"NICHI"))
      {
        nichi[0]=strtod(*(data+1),NULL);
        nichi[1]=strtod(*(data+2),NULL);
        nichi[2]=strtod(*(data+3),NULL);
      }
      freestr(data,nstr);

      data=fgetsbrk(fin,&nstr);
      if(nstr==0) return;
      if(!strcmp(*(data+0),"NFLOOR"))
      {
        nframe=(int)strtol(*(data+1),NULL,10);
      }
      freestr(data,nstr);

      for(k=0;k<nframe;k++)
      {
        (wdraw.childs+1)->vparam.subsection++;

        if((jj==0||jj==1) && k==0) (wdraw.childs+1)->vparam.vflag.nv.conffig=1;
        else                       (wdraw.childs+1)->vparam.vflag.nv.conffig=0;

        data=fgetsbrk(fin,&nstr);
        if(nstr==0) return;

        if(!strcmp(*(data+0),"FRAME"))
        {
            sprintf(fname,"%s",*(data+1));
        }
/*
        if(!strcmp(*(data+2),"HEAD"))
        {
          vp->range.min.d[GX]=strtod(*(data+3),NULL);
          vp->range.min.d[GY]=strtod(*(data+4),NULL);
          vp->range.min.d[GZ]=strtod(*(data+5),NULL);
        }
        if(!strcmp(*(data+6),"TAIL"))
        {
          vp->range.max.d[GX]=strtod(*(data+7),NULL);
          vp->range.max.d[GY]=strtod(*(data+8),NULL);
          vp->range.max.d[GZ]=strtod(*(data+9),NULL);
        }
*/
        /*Ai分布から各フロアレベル読み込み*/
        if(jj==1||jj==2)
        {
          vp->range.min.d[GZ]=*(((wdraw.childs+1)->org).ai.lbound+ k);
          vp->range.max.d[GZ]=*(((wdraw.childs+1)->org).ai.lbound+ k+1);
        }
        if(jj==0||jj==3)
        {
          vp->range.min.d[GZ]=*(((wdraw.childs+1)->org).ai.lbound+ k);
          vp->range.max.d[GZ]=*(((wdraw.childs+1)->org).ai.lbound+ k+2);
          if(k+2==nframe) k++;
        }
        SendMessage((wmenu.childs+2)->hwnd,WM_INITDIALOG,0,0);

        /*DISPLACEMENTN*/
        vp->focus.d[0]=nichi[0];
        vp->focus.d[1]=nichi[1];
        vp->focus.d[2]=nichi[2];
        (wdraw.childs+1)->vparam.vflag.ev.deformation=0;        //UJIOKA
        SendMessage((wmenu.childs+2)->hwnd,WM_INITDIALOG,0,0);

        if(ftype==3)
        {
          sprintf(str,"%s%s",fname,s);
          SetDlgItemText((wmenu.childs+3)->hwnd,IDV_TITLEBOX,
                           str);

          wparam = MAKEWPARAM(IDV_DEFORMATION,0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,wparam,0);
          if(i==0) SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                               MAKEWPARAM(IDV_DZ,0),0);
          if(i==1) SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                               MAKEWPARAM(IDV_DX,0),0);
          if(i==2) SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                               MAKEWPARAM(IDV_DY,0),0);

          SendMessage(wmain.hwnd,WM_COMMAND,
                      MAKEWPARAM((WORD)IDM_PRINTWAIT,(WORD)0),0);

          wparam = MAKEWPARAM(IDV_DEFORMATION,0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,wparam,0);
          if(i==0) SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                               MAKEWPARAM(IDV_DZ,0),0);
          if(i==1) SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                               MAKEWPARAM(IDV_DX,0),0);
          if(i==2) SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                               MAKEWPARAM(IDV_DY,0),0);

          /*END PAGE*/
          SendMessage(wmain.hwnd,WM_COMMAND,
                      MAKEWPARAM((WORD)IDM_PRINTOUT,(WORD)0),0);
        }

        /*N*/
        vp->focus.d[0]=nichi[0];
        vp->focus.d[1]=nichi[1];
        vp->focus.d[2]=nichi[2];
        SendMessage((wmenu.childs+2)->hwnd,WM_INITDIALOG,0,0);

        if(ftype==1 || ftype==2)
        {
          if(ftype==2)
          {
            sprintf(str,"%s",fname);
            SetDlgItemText((wmenu.childs+3)->hwnd,IDV_TITLEBOX,
                           str);
          }

          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_NZ,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_NZ_G,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_NZ_B,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_NZ_W,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_NZ_S,0),0);

          SendMessage(wmain.hwnd,WM_COMMAND,
                      MAKEWPARAM((WORD)IDM_PRINTWAIT,(WORD)0),0);

          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_NZ,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_NZ_G,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_NZ_B,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_NZ_W,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_NZ_S,0),0);

          /*END PAGE*/
          if(ftype==2) SendMessage(wmain.hwnd,WM_COMMAND,
                                   MAKEWPARAM((WORD)IDM_PRINTOUT,(WORD)0),0);
        }

        /*Qx*/
        vp->focus.d[0]=qichi[0];
        vp->focus.d[1]=qichi[1];
        vp->focus.d[2]=qichi[2];
        SendMessage((wmenu.childs+2)->hwnd,WM_INITDIALOG,0,0);

        if(ftype==1)
        {
          if(j==0)  SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                                MAKEWPARAM(IDV_QY,0),0);
          if(j==1)  SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                                MAKEWPARAM(IDV_QX,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_QY_G,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_QY_B,0),0);

          SendMessage(wmain.hwnd,WM_COMMAND,
                      MAKEWPARAM((WORD)IDM_PRINTWAIT,(WORD)0),0);

          if(j==0)  SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                                MAKEWPARAM(IDV_QY,0),0);
          if(j==1)  SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                                MAKEWPARAM(IDV_QX,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_QY_G,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_QY_B,0),0);

          /*END PAGE*/
          if(ftype==2) SendMessage(wmain.hwnd,WM_COMMAND,
                                   MAKEWPARAM((WORD)IDM_PRINTOUT,(WORD)0),0);
        }
        if(ftype==2)
        {
          /*Qx*/
          sprintf(str,"%s",fname);
          SetDlgItemText((wmenu.childs+3)->hwnd,IDV_TITLEBOX,str);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_QX,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_QX_G,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_QX_B,0),0);

          SendMessage(wmain.hwnd,WM_COMMAND,
                      MAKEWPARAM((WORD)IDM_PRINTWAIT,(WORD)0),0);

          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_QX,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_QX_G,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_QX_B,0),0);

          /*END PAGE*/
          if(ftype==2) SendMessage(wmain.hwnd,WM_COMMAND,
                                   MAKEWPARAM((WORD)IDM_PRINTOUT,(WORD)0),0);

          /*Qy*/
		  sprintf(str,"%s",fname);
          SetDlgItemText((wmenu.childs+3)->hwnd,IDV_TITLEBOX,str);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_QY,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_QY_G,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_QY_B,0),0);

          SendMessage(wmain.hwnd,WM_COMMAND,
                      MAKEWPARAM((WORD)IDM_PRINTWAIT,(WORD)0),0);

          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_QY,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_QY_G,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_QY_B,0),0);

          /*END PAGE*/
          if(ftype==2) SendMessage(wmain.hwnd,WM_COMMAND,
                                   MAKEWPARAM((WORD)IDM_PRINTOUT,(WORD)0),0);
        }

        /*M*/
        vp->focus.d[0]=michi[0];
        vp->focus.d[1]=michi[1];
        vp->focus.d[2]=michi[2];
        SendMessage((wmenu.childs+2)->hwnd,WM_INITDIALOG,0,0);

        if(ftype==1)
        {
          if(j==0)  SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                                MAKEWPARAM(IDV_MX,0),0);
          if(j==1)  SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                                MAKEWPARAM(IDV_MY,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_MX_G,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_MX_B,0),0);

          SendMessage(wmain.hwnd,WM_COMMAND,
                      MAKEWPARAM((WORD)IDM_PRINTWAIT,(WORD)0),0);

          if(j==0)  SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                                MAKEWPARAM(IDV_MX,0),0);
          if(j==1)  SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                                MAKEWPARAM(IDV_MY,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_MX_G,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_MX_B,0),0);

          /*END PAGE*/
          SendMessage(wmain.hwnd,WM_COMMAND,
                      MAKEWPARAM((WORD)IDM_PRINTOUT,(WORD)0),0);
        }
        if(ftype==2)
        {
          /*Mx*/
          sprintf(str,"%s",fname);
          SetDlgItemText((wmenu.childs+3)->hwnd,IDV_TITLEBOX,str);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_MX,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_MX_G,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_MX_B,0),0);

          SendMessage(wmain.hwnd,WM_COMMAND,
                      MAKEWPARAM((WORD)IDM_PRINTWAIT,(WORD)0),0);

          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_MX,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_MX_G,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_MX_B,0),0);

          /*END PAGE*/
          SendMessage(wmain.hwnd,WM_COMMAND,
                      MAKEWPARAM((WORD)IDM_PRINTOUT,(WORD)0),0);

          /*My*/
          sprintf(str,"%s",fname);
          SetDlgItemText((wmenu.childs+3)->hwnd,IDV_TITLEBOX,str);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_MX,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_MX_G,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_MX_B,0),0);

          SendMessage(wmain.hwnd,WM_COMMAND,
                      MAKEWPARAM((WORD)IDM_PRINTWAIT,(WORD)0),0);

          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_MX,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_MX_G,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_MX_B,0),0);

          /*END PAGE*/
          SendMessage(wmain.hwnd,WM_COMMAND,
                      MAKEWPARAM((WORD)IDM_PRINTOUT,(WORD)0),0);
        }

        freestr(data,nstr);
      }
     }
    }
  }

  if(!mode)
  {
  /*PRINT END.*/
  wparam = MAKEWPARAM((WORD)IDM_PRINTEND,(WORD)0);
  SendMessage(wmain.hwnd,WM_COMMAND,wparam,0);
  }
  vp->range.min.d[GX]=-100.0;
  vp->range.min.d[GY]=-100.0;
  vp->range.min.d[GZ]=-100.0;
  vp->range.max.d[GX]= 100.0;
  vp->range.max.d[GY]= 100.0;
  vp->range.max.d[GZ]= 100.0;
  vp->focus.d[0]     =initd[0];
  vp->focus.d[1]     =initd[1];
  vp->focus.d[2]     =initd[2];
  vp->theta          =inittheta;
  vp->phi            =initphi;
  vp->r              =initr;
  vp->odv            =initodv;
  SendMessage((wmenu.childs+2)->hwnd,WM_INITDIALOG,0,0);

  clearwindow(*(wdraw.childs+1));
  drawarclmframe((wdraw.childs+1)->hdcC,
                 (wdraw.childs+1)->vparam,arc,0,ONSCREEN);
  SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);


  return;
}/*printsrcanfigures*/

void printsrcanmax(FILE *fin,struct viewparam *vp,int mode)
{
  WPARAM wparam;
  char **data,str[256],s[80];
  int i,j,jj,k,n,np1,np2,nstr,flag/*,ftype*/;

  char fname[80];
  int nframe;
  double nichi[3],qichi[3],michi[3];

  double initd[3];
  double inittheta,initphi;
  double initr,initodv;

  getviewparam((wmenu.childs+2)->hwnd,
              &((wdraw.childs+1)->vparam));

  initd[0] =vp->focus.d[0];
  initd[1] =vp->focus.d[1];
  initd[2] =vp->focus.d[2];
  inittheta=vp->theta;
  initphi  =vp->phi;
  initr    =vp->r;
  initodv  =vp->odv;

  /*PRINTER SETUP*/
  SendMessage(wmain.hwnd,WM_COMMAND,
              MAKEWPARAM((WORD)IDM_PRINTSETUP,(WORD)0),0);

  SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,wparam,0);
  {
    wparam = MAKEWPARAM(IDV_DEFORMATION,0);
    SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,wparam,0);
    wparam = MAKEWPARAM(IDV_TITLE,0);
    SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,wparam,0);
    wparam = MAKEWPARAM(IDV_VIEW,0);
    SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,wparam,0);
    wparam = MAKEWPARAM(IDV_INPUTFILE,0);
    SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,wparam,0);
    wparam = MAKEWPARAM(IDV_OUTPUTFILE,0);
    SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,wparam,0);
    (wdraw.childs+1)->vparam.vflag.nv.code=0;
    (wdraw.childs+1)->vparam.vflag.ev.sectioncode=0;

    (wdraw.childs+1)->vparam.chapter=4;      /*CHAPTER*/
    (wdraw.childs+1)->vparam.section=3-1;    /*INITIAL SECTION*/
  }

  {
    fseek(fin,0L,SEEK_SET);

    (wdraw.childs+1)->vparam.section++;
    (wdraw.childs+1)->vparam.subsection=0;

    /*OPEN FILE AND DRAW FRAME.*/
    wparam = MAKEWPARAM(IDD_VIEW,0);
    SendMessage((wmenu.childs+2)->hwnd,WM_COMMAND,wparam,0);
    globalmessageflag=0;

    n=strcspn((wdraw.childs+1)->otpfile,".");
    strncpy(str,(wdraw.childs+1)->otpfile,n);
    str[n]='\0';
    strcpy((wdraw.childs+1)->otpfile,str);
    strcat((wdraw.childs+1)->otpfile,".rat");
    SetDlgItemText((wmenu.childs+2)->hwnd,ID_OUTPUTFILE,
	    			 (wdraw.childs+1)->otpfile);
    fout=fopen((wdraw.childs+1)->otpfile,"r");
//    fout=fgetstofopen("\0","r",ID_OUTPUTFILE);
    readsrcanrate(fout,&arc);

    (wdraw.childs+1)->vparam.vflag.ev.srcanrate  =0;
    (wdraw.childs+1)->vparam.vflag.ev.srcancolor =0;
    (wdraw.childs+1)->vparam.vflag.ev.srcanmax   =1;
    (wdraw.childs+1)->vparam.vflag.mv.pagetitle  =1; /****SRCANMAX****/

    wparam = MAKEWPARAM(IDV_UNITNM,0);
    SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,wparam,0);

    /*PRINT*/
    {
  vp->range.min.d[GX]=-100.0;
  vp->range.min.d[GY]=-100.0;
  vp->range.min.d[GZ]=-100.0;
  vp->range.max.d[GX]= 100.0;
  vp->range.max.d[GY]= 100.0;
  vp->range.max.d[GZ]= 100.0;
  vp->focus.d[0]     =initd[0];
  vp->focus.d[1]     =initd[1];
  vp->focus.d[2]     =initd[2];
  vp->theta          =inittheta;
  vp->phi            =initphi;
  vp->r              =initr;
  vp->odv            =initodv;
       (wdraw.childs+1)->vparam.vflag.mv.inputfile  =1;
       (wdraw.childs+1)->vparam.vflag.mv.outputfile =1;
       (wdraw.childs+1)->vparam.vflag.mv.view       =1;

/*
       setviewpoint((wdraw.childs+0)->hwnd,arc,
                   &((wdraw.childs+1)->vparam));
       setviewparam((wmenu.childs+2)->hwnd,
                    (wdraw.childs+1)->vparam);
       SendMessage((wmenu.childs+2)->hwnd,WM_INITDIALOG,0,0);
*/
       SendMessage((wmenu.childs+2)->hwnd,WM_INITDIALOG,0,0);

       clearwindow(*(wdraw.childs+1));
       drawarclmframe((wdraw.childs+1)->hdcC,
                      (wdraw.childs+1)->vparam,arc,0,ONSCREEN);
       SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);
       (wdraw.childs+1)->vparam.subsection++;
       SendMessage(wmain.hwnd,WM_COMMAND,
                   MAKEWPARAM((WORD)IDM_PRINTWAIT,(WORD)0),0);
       /*END PAGE*/
       SendMessage(wmain.hwnd,WM_COMMAND,
                   MAKEWPARAM((WORD)IDM_PRINTOUT,(WORD)0),0);
    }
  }
  if(!mode)
  {
  /*PRINT END.*/
  wparam = MAKEWPARAM((WORD)IDM_PRINTEND,(WORD)0);
  SendMessage(wmain.hwnd,WM_COMMAND,wparam,0);
  }
  vp->range.min.d[GX]=-100.0;
  vp->range.min.d[GY]=-100.0;
  vp->range.min.d[GZ]=-100.0;
  vp->range.max.d[GX]= 100.0;
  vp->range.max.d[GY]= 100.0;
  vp->range.max.d[GZ]= 100.0;
  vp->focus.d[0]     =initd[0];
  vp->focus.d[1]     =initd[1];
  vp->focus.d[2]     =initd[2];
  vp->theta          =inittheta;
  vp->phi            =initphi;
  vp->r              =initr;
  vp->odv            =initodv;
  (wdraw.childs+1)->vparam.vflag.ev.srcanmax   =0;
  SendMessage((wmenu.childs+2)->hwnd,WM_INITDIALOG,0,0);

  clearwindow(*(wdraw.childs+1));
  drawarclmframe((wdraw.childs+1)->hdcC,
                 (wdraw.childs+1)->vparam,arc,0,ONSCREEN);
  SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);


  return;
}/*printsrcanmax*/

void printarclmlastfigures(FILE *fin,struct viewparam *vp,
                           struct arclmframe *af)
/*PRINT ARCLM101 LAST STRESSES.*/
{
  WPARAM wparam;
  char **data,str[256];
  int i,j,k,np1,np2,nstr,flag;

  FILE *fout;
  char fname[80];
  int nframe;
  double nichi[3],qichi[3],michi[3];

  /*PRINTER SETUP*/
  SendMessage(wmain.hwnd,WM_COMMAND,
              MAKEWPARAM((WORD)IDM_PRINTSETUP,(WORD)0),0);

  /*VISIBLE SETTING*/
  wparam = MAKEWPARAM(IDV_GLOBALAXIS,0);
  SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,wparam,0);
  wparam = MAKEWPARAM(IDV_NODECODE,0);
  SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,wparam,0);
  wparam = MAKEWPARAM(IDV_SLAB,0);
  SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,wparam,0);
  wparam = MAKEWPARAM(IDV_PAGETITLE,0);
  SendMessage((wmenu.childs+4)->hwnd,WM_COMMAND,wparam,0);

  for(i=0;i<2;i++) /*FOR PERIOD X,Y*/
  {
    fseek(fin,0L,SEEK_SET);
    if(i==0)
    {
      np1=0; np2=0;

      GetDlgItemText((wmenu.childs+2)->hwnd,ID_INPUTFILEX,str,256);
      sprintf((wdraw.childs+1)->inpfile,str);
      SetDlgItemText((wmenu.childs+2)->hwnd,ID_INPUTFILE,
                     (wdraw.childs+1)->inpfile);
      GetDlgItemText((wmenu.childs+2)->hwnd,ID_OUTPUTFILEX,str,256);
      sprintf((wdraw.childs+1)->otpfile,str);
      SetDlgItemText((wmenu.childs+2)->hwnd,ID_OUTPUTFILE,
                     (wdraw.childs+1)->otpfile);

    }

    if(i==1)
    {
      np1=1; np2=1;

      GetDlgItemText((wmenu.childs+2)->hwnd,ID_INPUTFILEY,str,256);
      sprintf((wdraw.childs+1)->inpfile,str);
      SetDlgItemText((wmenu.childs+2)->hwnd,ID_INPUTFILE,
                     (wdraw.childs+1)->inpfile);
      GetDlgItemText((wmenu.childs+2)->hwnd,ID_OUTPUTFILEY,str,256);
      sprintf((wdraw.childs+1)->otpfile,str);
      SetDlgItemText((wmenu.childs+2)->hwnd,ID_OUTPUTFILE,
                     (wdraw.childs+1)->otpfile);
    }

    /*OPEN FILE AND DRAW FRAME.*/
    wparam = MAKEWPARAM(IDD_VIEW,0);
    SendMessage((wmenu.childs+2)->hwnd,WM_COMMAND,wparam,0);

    fout=fgetstofopen("\0","r",ID_OUTPUTFILE);  /*OPEN FILE.*/
    if(fout==NULL) return;
    openarclmlastfileII(fout,af);
    fclose(fout);

    wparam = MAKEWPARAM(IDV_UNITNM,0);
    SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,wparam,0);

    for(j=np1;j<=np2;j++) /*FOR FRAME X,Y*/
    {
      if(j==0)
      {
        flag=0;
        while(!flag)
        {
          data=fgetsbrk(fin,&nstr);
          if(!strcmp(*(data+0),"XFRAME")) flag=1;
        }
      }
      if(j==1)
      {
        flag=0;
        while(!flag)
        {
          data=fgetsbrk(fin,&nstr);
          if(!strcmp(*(data+0),"YFRAME")) flag=1;
        }
      }

	  data=fgetsbrk(fin,&nstr);
      if(nstr==0) return;
      if(!strcmp(*(data+0),"ANGLE"))
      {
        vp->theta=strtod(*(data+1),NULL);
        vp->phi  =strtod(*(data+2),NULL);
      }
      freestr(data,nstr);

      data=fgetsbrk(fin,&nstr);
      if(nstr==0) return;
      if(!strcmp(*(data+0),"DISTS"))
      {
        vp->r  =strtod(*(data+1),NULL);
        vp->odv=strtod(*(data+2),NULL);
      }
      freestr(data,nstr);

      data=fgetsbrk(fin,&nstr);
      if(nstr==0) return;
      if(!strcmp(*(data+0),"NICHI"))
      {
        nichi[0]=strtod(*(data+1),NULL);
        nichi[1]=strtod(*(data+2),NULL);
        nichi[2]=strtod(*(data+3),NULL);
      }
      freestr(data,nstr);

      data=fgetsbrk(fin,&nstr);
      if(nstr==0) return;
      if(!strcmp(*(data+0),"QICHI"))
      {
        qichi[0]=strtod(*(data+1),NULL);
        qichi[1]=strtod(*(data+2),NULL);
        qichi[2]=strtod(*(data+3),NULL);
      }
      freestr(data,nstr);

      data=fgetsbrk(fin,&nstr);
      if(nstr==0) return;
      if(!strcmp(*(data+0),"MICHI"))
      {
        michi[0]=strtod(*(data+1),NULL);
        michi[1]=strtod(*(data+2),NULL);
        michi[2]=strtod(*(data+3),NULL);
      }
      freestr(data,nstr);

      data=fgetsbrk(fin,&nstr);
      if(nstr==0) return;
      if(!strcmp(*(data+0),"NFRAME"))
      {
        nframe=(int)strtol(*(data+1),NULL,10);
      }
      freestr(data,nstr);

      for(k=0;k<nframe;k++)
      {
		data=fgetsbrk(fin,&nstr);
        if(nstr==0) return;

        if(!strcmp(*(data+0),"FRAME"))
        {
          sprintf(fname,"%s 通り応力図",*(data+1));

          SetDlgItemText((wmenu.childs+4)->hwnd,IDV_PAGETITLEBOX,
                         fname);
        }
        if(!strcmp(*(data+2),"HEAD"))
        {
          vp->range.min.d[GX]=strtod(*(data+3),NULL);
          vp->range.min.d[GY]=strtod(*(data+4),NULL);
          vp->range.min.d[GZ]=strtod(*(data+5),NULL);
        }
        if(!strcmp(*(data+6),"TAIL"))
        {
          vp->range.max.d[GX]=strtod(*(data+7),NULL);
          vp->range.max.d[GY]=strtod(*(data+8),NULL);
          vp->range.max.d[GZ]=strtod(*(data+9),NULL);
        }

        SendMessage((wmenu.childs+2)->hwnd,WM_INITDIALOG,0,0);

        /*N*/
        vp->focus.d[0]=nichi[0];
        vp->focus.d[1]=nichi[1];
        vp->focus.d[2]=nichi[2];
        SendMessage((wmenu.childs+2)->hwnd,WM_INITDIALOG,0,0);

        SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                    MAKEWPARAM(IDV_NZ,0),0);
        SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                    MAKEWPARAM(IDV_NZ_G,0),0);
        SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                    MAKEWPARAM(IDV_NZ_B,0),0);
        SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                    MAKEWPARAM(IDV_NZ_W,0),0);
        SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                    MAKEWPARAM(IDV_NZ_S,0),0);

        SendMessage(wmain.hwnd,WM_COMMAND,
                    MAKEWPARAM((WORD)IDM_PRINTWAIT,(WORD)0),0);

        SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                    MAKEWPARAM(IDV_NZ,0),0);
        SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                    MAKEWPARAM(IDV_NZ_G,0),0);
        SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                    MAKEWPARAM(IDV_NZ_B,0),0);
        SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                    MAKEWPARAM(IDV_NZ_W,0),0);
        SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                    MAKEWPARAM(IDV_NZ_S,0),0);

        /*Q*/
        vp->focus.d[0]=qichi[0];
        vp->focus.d[1]=qichi[1];
        vp->focus.d[2]=qichi[2];
        SendMessage((wmenu.childs+2)->hwnd,WM_INITDIALOG,0,0);

        if(j==0)  SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                              MAKEWPARAM(IDV_QY,0),0);
        if(j==1)  SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                              MAKEWPARAM(IDV_QX,0),0);
        SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                    MAKEWPARAM(IDV_QY_G,0),0);
        SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                    MAKEWPARAM(IDV_QY_B,0),0);

        SendMessage(wmain.hwnd,WM_COMMAND,
                    MAKEWPARAM((WORD)IDM_PRINTWAIT,(WORD)0),0);

        if(j==0)  SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                              MAKEWPARAM(IDV_QY,0),0);
        if(j==1)  SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                              MAKEWPARAM(IDV_QX,0),0);
        SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                    MAKEWPARAM(IDV_QY_G,0),0);
        SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                    MAKEWPARAM(IDV_QY_B,0),0);

        /*M*/
        vp->focus.d[0]=michi[0];
        vp->focus.d[1]=michi[1];
        vp->focus.d[2]=michi[2];
        SendMessage((wmenu.childs+2)->hwnd,WM_INITDIALOG,0,0);

        if(j==0)  SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                              MAKEWPARAM(IDV_MX,0),0);
        if(j==1)  SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                              MAKEWPARAM(IDV_MY,0),0);
        SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                    MAKEWPARAM(IDV_MX_G,0),0);
        SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                    MAKEWPARAM(IDV_MX_B,0),0);

        SendMessage(wmain.hwnd,WM_COMMAND,
                    MAKEWPARAM((WORD)IDM_PRINTWAIT,(WORD)0),0);

        if(j==0)  SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                              MAKEWPARAM(IDV_MX,0),0);
        if(j==1)  SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                              MAKEWPARAM(IDV_MY,0),0);
        SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                    MAKEWPARAM(IDV_MX_G,0),0);
        SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                    MAKEWPARAM(IDV_MX_B,0),0);

        /*END PAGE*/
        SendMessage(wmain.hwnd,WM_COMMAND,
                    MAKEWPARAM((WORD)IDM_PRINTOUT,(WORD)0),0);

        freestr(data,nstr);
      }
    }
  }

  /*PRINT END.*/
  wparam = MAKEWPARAM((WORD)IDM_PRINTEND,(WORD)0);
  SendMessage(wmain.hwnd,WM_COMMAND,wparam,0);

  return;
}/*printarclmlastfigures*/


void printarclmlastfiguresII(FILE *fin,struct viewparam *vp,
							 struct arclmframe *af)        //UJIOKA
/*TYPE 1:PRINT DEFORMATION BY AXIS.*/
/*SETTING BEFORE BEGINNING:VISIBLE ELEMENT TYPE*/
/*                         M FACTOR*/
/*                         H SIZE*/
{
  WPARAM wparam;
  char **data,str[256],s[80];
  int i,ii,j,k,np1,np2,nstr,flag,ftype,totallaps;

  FILE *fout;
  char fname[80];
  int nframe;
  double nichi[3],qichi[3],michi[3];

  ftype=1; /*FTYPE:1=DEFORMATION BY AXIS*/

  /*PRINTER SETUP*/
  SendMessage(wmain.hwnd,WM_COMMAND,
              MAKEWPARAM((WORD)IDM_PRINTSETUP,(WORD)0),0);

  /*VISIBLE SETTING*/
  wparam = MAKEWPARAM(IDV_PAGETITLE,0);
  SendMessage((wmenu.childs+4)->hwnd,WM_COMMAND,wparam,0);


  getlaps((wmenu.childs+2)->hwnd,&totallaps);
  for(ii=1;ii<=totallaps;ii++)  /*LAPS*/
  {
    sprintf(s,"%d",ii);
    SetDlgItemText((wmenu.childs+2)->hwnd,ID_LAPS,s);
    SendDlgItemMessage((wmenu.childs+2)->hwnd,ID_LAPS,WM_PAINT,0,0);
    fseek(fin,0L,SEEK_SET);

    if(ftype==1){np1=0; np2=0;}


    /*OPEN FILE AND DRAW FRAME.*/
    wparam = MAKEWPARAM(IDD_VIEW,0);
    SendMessage((wmenu.childs+2)->hwnd,WM_COMMAND,wparam,0);
/*
    wparam = MAKEWPARAM(IDD_OPENRESULT,0);
    SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,wparam,0);
*/
    fout=fgetstofopen("\0","r",ID_OUTPUTFILE);  /*OPEN FILE.*/
    if(fout==NULL) return;
	openarclmlastfileII(fout,af);
    fclose(fout);
    wparam = MAKEWPARAM(IDV_UNITNM,0);
    SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,wparam,0);

    for(j=np1;j<=np2;j++) /*FOR FRAME X,Y*/
    {
      if(j==0)
      {
        flag=0;
        while(!flag)
        {
		  data=fgetsbrk(fin,&nstr);
          if(!strcmp(*(data+0),"XFRAME")) flag=1;
        }
      }
      if(j==1)
      {
        flag=0;
        while(!flag)
        {
          data=fgetsbrk(fin,&nstr);
          if(!strcmp(*(data+0),"YFRAME")) flag=1;
        }
      }

      data=fgetsbrk(fin,&nstr);
      if(nstr==0) return;
      if(!strcmp(*(data+0),"ANGLE"))
      {
        vp->theta=strtod(*(data+1),NULL);
        vp->phi  =strtod(*(data+2),NULL);
      }
      freestr(data,nstr);

      data=fgetsbrk(fin,&nstr);
      if(nstr==0) return;
      if(!strcmp(*(data+0),"DISTS"))
      {
        vp->r  =strtod(*(data+1),NULL);
        vp->odv=strtod(*(data+2),NULL);
      }
      freestr(data,nstr);

	  data=fgetsbrk(fin,&nstr);
      if(nstr==0) return;
      if(!strcmp(*(data+0),"NICHI"))
      {
        nichi[0]=strtod(*(data+1),NULL);
        nichi[1]=strtod(*(data+2),NULL);
        nichi[2]=strtod(*(data+3),NULL);
      }
      freestr(data,nstr);

      data=fgetsbrk(fin,&nstr);
      if(nstr==0) return;
      if(!strcmp(*(data+0),"QICHI"))
      {
        qichi[0]=strtod(*(data+1),NULL);
        qichi[1]=strtod(*(data+2),NULL);
        qichi[2]=strtod(*(data+3),NULL);
      }
      freestr(data,nstr);

      data=fgetsbrk(fin,&nstr);
      if(nstr==0) return;
      if(!strcmp(*(data+0),"MICHI"))
      {
        michi[0]=strtod(*(data+1),NULL);
        michi[1]=strtod(*(data+2),NULL);
        michi[2]=strtod(*(data+3),NULL);
      }
      freestr(data,nstr);

	  data=fgetsbrk(fin,&nstr);
      if(nstr==0) return;
      if(!strcmp(*(data+0),"NFRAME"))
      {
        nframe=(int)strtol(*(data+1),NULL,10);
      }
      freestr(data,nstr);

      for(k=0;k<nframe;k++)
      {
        data=fgetsbrk(fin,&nstr);
        if(nstr==0) return;

        if(!strcmp(*(data+0),"FRAME"))
        {
          if(ftype==1)
          {
//            sprintf(fname,"LAP:%d:%s 通り変形図",ii,*(data+1));
            sprintf(fname,"増分解析変形図:LAP:%d",ii,*(data+1));
            SetDlgItemText((wmenu.childs+4)->hwnd,IDV_PAGETITLEBOX,
                           fname);
          }
        }
        if(!strcmp(*(data+2),"HEAD"))
        {
          vp->range.min.d[GX]=strtod(*(data+3),NULL);
          vp->range.min.d[GY]=strtod(*(data+4),NULL);
          vp->range.min.d[GZ]=strtod(*(data+5),NULL);
        }
        if(!strcmp(*(data+6),"TAIL"))
        {
          vp->range.max.d[GX]=strtod(*(data+7),NULL);
          vp->range.max.d[GY]=strtod(*(data+8),NULL);
          vp->range.max.d[GZ]=strtod(*(data+9),NULL);
        }

        SendMessage((wmenu.childs+2)->hwnd,WM_INITDIALOG,0,0);

        /*DISPLACEMENTN*/
        vp->focus.d[0]=nichi[0];
        vp->focus.d[1]=nichi[1];
        vp->focus.d[2]=nichi[2];
        SendMessage((wmenu.childs+2)->hwnd,WM_INITDIALOG,0,0);

        if(ftype==1)
        {
/*
          sprintf(str,"変形図 %s",fname);
          SetDlgItemText((wmenu.childs+3)->hwnd,IDV_TITLEBOX,
                           str);
*/
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,wparam,0);
          if(i==0) SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                               MAKEWPARAM(IDV_DZ,0),0);
          if(i==1) SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                               MAKEWPARAM(IDV_DX,0),0);
          if(i==2) SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                               MAKEWPARAM(IDV_DY,0),0);

          SendMessage(wmain.hwnd,WM_COMMAND,
                      MAKEWPARAM((WORD)IDM_PRINTWAIT,(WORD)0),0);

          wparam = MAKEWPARAM(IDV_DEFORMATION,0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,wparam,0);
          if(i==0) SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                               MAKEWPARAM(IDV_DZ,0),0);
          if(i==1) SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                               MAKEWPARAM(IDV_DX,0),0);
          if(i==2) SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                               MAKEWPARAM(IDV_DY,0),0);

          /*END PAGE*/
          SendMessage(wmain.hwnd,WM_COMMAND,
                      MAKEWPARAM((WORD)IDM_PRINTOUT,(WORD)0),0);
        }
#if 0
        /*N*/
        vp->focus.d[0]=nichi[0];
        vp->focus.d[1]=nichi[1];
        vp->focus.d[2]=nichi[2];
        SendMessage((wmenu.childs+2)->hwnd,WM_INITDIALOG,0,0);

        if(ftype==1 || ftype==2)
        {
          if(ftype==2)
          {
            sprintf(str,"軸力 Nz %s",fname);
            SetDlgItemText((wmenu.childs+3)->hwnd,IDV_TITLEBOX,
                           str);
          }

          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_NZ,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_NZ_G,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_NZ_B,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_NZ_W,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_NZ_S,0),0);

          SendMessage(wmain.hwnd,WM_COMMAND,
                      MAKEWPARAM((WORD)IDM_PRINTWAIT,(WORD)0),0);

          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_NZ,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_NZ_G,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_NZ_B,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_NZ_W,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_NZ_S,0),0);

          /*END PAGE*/
          if(ftype==2) SendMessage(wmain.hwnd,WM_COMMAND,
                                   MAKEWPARAM((WORD)IDM_PRINTOUT,(WORD)0),0);
        }

        /*Qx*/
        vp->focus.d[0]=qichi[0];
        vp->focus.d[1]=qichi[1];
        vp->focus.d[2]=qichi[2];
        SendMessage((wmenu.childs+2)->hwnd,WM_INITDIALOG,0,0);

        if(ftype==1)
        {
          if(j==0)  SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                                MAKEWPARAM(IDV_QY,0),0);
          if(j==1)  SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                                MAKEWPARAM(IDV_QX,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_QY_G,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_QY_B,0),0);

          SendMessage(wmain.hwnd,WM_COMMAND,
                      MAKEWPARAM((WORD)IDM_PRINTWAIT,(WORD)0),0);

          if(j==0)  SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                                MAKEWPARAM(IDV_QY,0),0);
          if(j==1)  SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                                MAKEWPARAM(IDV_QX,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_QY_G,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_QY_B,0),0);

          /*END PAGE*/
          if(ftype==2) SendMessage(wmain.hwnd,WM_COMMAND,
                                   MAKEWPARAM((WORD)IDM_PRINTOUT,(WORD)0),0);
        }
        if(ftype==2)
        {
          /*Qx*/
          sprintf(str,"せん断力 Qx %s",fname);
          SetDlgItemText((wmenu.childs+3)->hwnd,IDV_TITLEBOX,str);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_QX,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_QX_G,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_QX_B,0),0);

          SendMessage(wmain.hwnd,WM_COMMAND,
                      MAKEWPARAM((WORD)IDM_PRINTWAIT,(WORD)0),0);

          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_QX,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_QX_G,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_QX_B,0),0);

          /*END PAGE*/
          if(ftype==2) SendMessage(wmain.hwnd,WM_COMMAND,
                                   MAKEWPARAM((WORD)IDM_PRINTOUT,(WORD)0),0);

          /*Qy*/
          sprintf(str,"せん断力 Qy %s",fname);
          SetDlgItemText((wmenu.childs+3)->hwnd,IDV_TITLEBOX,str);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_QY,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_QY_G,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_QY_B,0),0);

          SendMessage(wmain.hwnd,WM_COMMAND,
                      MAKEWPARAM((WORD)IDM_PRINTWAIT,(WORD)0),0);

          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_QY,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_QY_G,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_QY_B,0),0);

          /*END PAGE*/
          if(ftype==2) SendMessage(wmain.hwnd,WM_COMMAND,
                                   MAKEWPARAM((WORD)IDM_PRINTOUT,(WORD)0),0);
        }

        /*M*/
        vp->focus.d[0]=michi[0];
        vp->focus.d[1]=michi[1];
        vp->focus.d[2]=michi[2];
        SendMessage((wmenu.childs+2)->hwnd,WM_INITDIALOG,0,0);

        if(ftype==1)
        {
          if(j==0)  SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                                MAKEWPARAM(IDV_MX,0),0);
          if(j==1)  SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                                MAKEWPARAM(IDV_MY,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_MX_G,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_MX_B,0),0);

          SendMessage(wmain.hwnd,WM_COMMAND,
                      MAKEWPARAM((WORD)IDM_PRINTWAIT,(WORD)0),0);

          if(j==0)  SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                                MAKEWPARAM(IDV_MX,0),0);
          if(j==1)  SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                                MAKEWPARAM(IDV_MY,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_MX_G,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_MX_B,0),0);

          /*END PAGE*/
          SendMessage(wmain.hwnd,WM_COMMAND,
                      MAKEWPARAM((WORD)IDM_PRINTOUT,(WORD)0),0);
        }
        if(ftype==2)
        {
          /*Mx*/
          sprintf(str,"曲げモーメント Mx %s",fname);
          SetDlgItemText((wmenu.childs+3)->hwnd,IDV_TITLEBOX,str);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_MX,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_MX_G,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_MX_B,0),0);

          SendMessage(wmain.hwnd,WM_COMMAND,
                      MAKEWPARAM((WORD)IDM_PRINTWAIT,(WORD)0),0);

          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_MX,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_MX_G,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_MX_B,0),0);

          /*END PAGE*/
          SendMessage(wmain.hwnd,WM_COMMAND,
                      MAKEWPARAM((WORD)IDM_PRINTOUT,(WORD)0),0);

          /*My*/
          sprintf(str,"曲げモーメント My %s",fname);
          SetDlgItemText((wmenu.childs+3)->hwnd,IDV_TITLEBOX,str);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_MX,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_MX_G,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_MX_B,0),0);

          SendMessage(wmain.hwnd,WM_COMMAND,
                      MAKEWPARAM((WORD)IDM_PRINTWAIT,(WORD)0),0);

          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_MX,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_MX_G,0),0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                      MAKEWPARAM(IDV_MX_B,0),0);

          /*END PAGE*/
          SendMessage(wmain.hwnd,WM_COMMAND,
                      MAKEWPARAM((WORD)IDM_PRINTOUT,(WORD)0),0);
        }
#endif
        freestr(data,nstr);
      }
    }
  }

  /*PRINT END.*/
  wparam = MAKEWPARAM((WORD)IDM_PRINTEND,(WORD)0);
  SendMessage(wmain.hwnd,WM_COMMAND,wparam,0);

  return;
}/*printarclmlastfiguresII*/


void printorganfigures(FILE *fin,struct viewparam *vp,int mode)
/*PRINT MODEL FIGURES BY FLOOR.*/
/*MODE 0:END PRINTING*/
/*     1:CONTINUE PRINTING*/
{
  WPARAM wparam;
  char **data,str[256],s[80];
  int i,j,jj,k,nstr,flag;

  char fname[80];
  int nframe;
  double nichi[3];

  double initd[3];
  double inittheta,initphi;
  double initr,initodv;

  getviewparam((wmenu.childs+2)->hwnd,
              &((wdraw.childs+1)->vparam));

  initd[0] =vp->focus.d[0];
  initd[1] =vp->focus.d[1];
  initd[2] =vp->focus.d[2];
  inittheta=vp->theta;
  initphi  =vp->phi;
  initr    =vp->r;
  initodv  =vp->odv;

  /*PRINTER SETUP*/
  SendMessage(wmain.hwnd,WM_COMMAND,
			  MAKEWPARAM((WORD)IDM_PRINTSETUP,(WORD)0),0);

  /*VISIBLE SETTING*/
  (wdraw.childs+1)->vparam.chapter=2;      /*CHAPTER*/
  (wdraw.childs+1)->vparam.section=2-1;    /*INITIAL SECTION*/


  {
    fseek(fin,0L,SEEK_SET);

    (wdraw.childs+1)->vparam.section++;
    (wdraw.childs+1)->vparam.subsection=0;

    /*OPEN FILE AND DRAW FRAME.*/
    wparam = MAKEWPARAM(IDD_VIEW,0);
    SendMessage((wmenu.childs+2)->hwnd,WM_COMMAND,wparam,0);
	globalmessageflag=0;

    wparam = MAKEWPARAM(IDV_UNITNM,0);
    SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,wparam,0);

    /*PRINT*/
    (wdraw.childs+1)->vparam.vflag.mv.pagetitle  =1;
    (wdraw.childs+1)->vparam.vflag.mv.title  =1;
    (wdraw.childs+1)->vparam.vflag.mv.inputfile  =1;
    (wdraw.childs+1)->vparam.vflag.mv.outputfile =0;
    (wdraw.childs+1)->vparam.vflag.mv.view       =1;
    (wdraw.childs+1)->vparam.vflag.nv.code=1;
    (wdraw.childs+1)->vparam.vflag.nv.d[0]=0;
    (wdraw.childs+1)->vparam.vflag.nv.d[1]=0;
    (wdraw.childs+1)->vparam.vflag.nv.d[2]=0;

    vp->range.min.d[GX]=-100.0;
    vp->range.min.d[GY]=-100.0;
    vp->range.min.d[GZ]=-100.0;
    vp->range.max.d[GX]= 100.0;
    vp->range.max.d[GY]= 100.0;
    vp->range.max.d[GZ]= 100.0;
    vp->focus.d[0]     =initd[0];
    vp->focus.d[1]     =initd[1];
    vp->focus.d[2]     =initd[2];
    vp->theta          =inittheta;
    vp->phi            =initphi;
    vp->r              =initr;
    vp->odv            =initodv;
    (wdraw.childs+1)->vparam.vflag.ev.etype[COLUMN]=1;
    (wdraw.childs+1)->vparam.vflag.ev.etype[GIRDER]=1;
    (wdraw.childs+1)->vparam.vflag.ev.etype[BEAM]=1;
    (wdraw.childs+1)->vparam.vflag.ev.etype[BRACE]=1;
    (wdraw.childs+1)->vparam.vflag.ev.etype[WALL]=1;
    (wdraw.childs+1)->vparam.vflag.ev.etype[SLAB]=1;
    SendMessage((wmenu.childs+2)->hwnd,WM_INITDIALOG,0,0);

    clearwindow(*(wdraw.childs+1));
    draworganization((wdraw.childs+1)->hdcC,
    				 (wdraw.childs+1)->vparam,
    				 (wdraw.childs+1)->org,ONSCREEN);
    SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);
    sprintf(str,"%d.%d：モデル図",
              (wdraw.childs+1)->vparam.chapter,
              (wdraw.childs+1)->vparam.section);
    SetDlgItemText((wmenu.childs+4)->hwnd,IDV_PAGETITLEBOX,str);
    sprintf(str,"解析モデル全体表示");
    SetDlgItemText((wmenu.childs+3)->hwnd,IDV_TITLEBOX,str);
    (wdraw.childs+1)->vparam.subsection++;
    SendMessage(wmain.hwnd,WM_COMMAND,
                MAKEWPARAM((WORD)IDM_PRINTWAIT,(WORD)0),0);
    /*END PAGE*/
    SendMessage(wmain.hwnd,WM_COMMAND,
                MAKEWPARAM((WORD)IDM_PRINTOUT,(WORD)0),0);


    (wdraw.childs+1)->vparam.vflag.mv.pagetitle  =0;
    (wdraw.childs+1)->vparam.vflag.nv.conffig=1;

    clearwindow(*(wdraw.childs+1));
    draworganization((wdraw.childs+1)->hdcC,
	    			 (wdraw.childs+1)->vparam,
    				 (wdraw.childs+1)->org,ONSCREEN);
    SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);

    (wdraw.childs+1)->vparam.vflag.ev.etype[BRACE] =0;
    (wdraw.childs+1)->vparam.vflag.ev.etype[WALL]  =0;
    (wdraw.childs+1)->vparam.vflag.ev.etype[SLAB]  =0;
    sprintf(str,"柱梁のみ表示");
    SetDlgItemText((wmenu.childs+3)->hwnd,IDV_TITLEBOX,str);
    (wdraw.childs+1)->vparam.subsection++;
    SendMessage(wmain.hwnd,WM_COMMAND,
                MAKEWPARAM((WORD)IDM_PRINTWAIT,(WORD)0),0);
    /*END PAGE*/
    SendMessage(wmain.hwnd,WM_COMMAND,
                MAKEWPARAM((WORD)IDM_PRINTOUT,(WORD)0),0);

    (wdraw.childs+1)->vparam.vflag.nv.code=0;
    (wdraw.childs+1)->vparam.vflag.ev.sectioncode=1;
#if 0
    (wdraw.childs+1)->vparam.vflag.ev.etype[GIRDER]=0;
    (wdraw.childs+1)->vparam.vflag.ev.etype[BEAM]  =0;
    sprintf(str,"柱のみ表示");
    SetDlgItemText((wmenu.childs+3)->hwnd,IDV_TITLEBOX,str);
    (wdraw.childs+1)->vparam.subsection++;
    SendMessage(wmain.hwnd,WM_COMMAND,
                   MAKEWPARAM((WORD)IDM_PRINTWAIT,(WORD)0),0);
    /*END PAGE*/
    SendMessage(wmain.hwnd,WM_COMMAND,
                MAKEWPARAM((WORD)IDM_PRINTOUT,(WORD)0),0);
#endif
  }

  for(jj=0;jj<4;jj++)
  {

       (wdraw.childs+1)->vparam.vflag.ev.etype[COLUMN]=0;
       (wdraw.childs+1)->vparam.vflag.ev.etype[GIRDER]=0;
       (wdraw.childs+1)->vparam.vflag.ev.etype[BEAM]=0;
       (wdraw.childs+1)->vparam.vflag.ev.etype[BRACE]=0;
       (wdraw.childs+1)->vparam.vflag.ev.etype[WALL]=0;
       (wdraw.childs+1)->vparam.vflag.ev.etype[SLAB]=0;

    if(jj==0)
    {
       (wdraw.childs+1)->vparam.vflag.ev.etype[COLUMN]=1;
       sprintf(s,"柱のみ表示");
    }
    if(jj==1)
    {
       (wdraw.childs+1)->vparam.vflag.ev.etype[GIRDER]=1;
       (wdraw.childs+1)->vparam.vflag.ev.etype[BEAM]=1;
       sprintf(s,"梁のみ表示");
    }
    if(jj==2)
    {
       (wdraw.childs+1)->vparam.vflag.ev.etype[SLAB]=1;
       sprintf(s,"床のみ表示");
    }
    if(jj==3)
    {
       (wdraw.childs+1)->vparam.vflag.ev.etype[WALL]=1;
       sprintf(s,"壁のみ表示");
    }


    (wdraw.childs+1)->vparam.vflag.nv.conffig=0;

    {
      fseek(fin,0L,SEEK_SET);
      {
        flag=0;
        while(!flag)
        {
          data=fgetsbrk(fin,&nstr);
          if(!strcmp(*(data+0),"ZFRAME")) flag=1;
        }
      }

      data=fgetsbrk(fin,&nstr);
      if(nstr==0) return;
      if(!strcmp(*(data+0),"ANGLE"))
      {
        vp->theta=strtod(*(data+1),NULL);
        vp->phi  =strtod(*(data+2),NULL);
      }
      freestr(data,nstr);

      data=fgetsbrk(fin,&nstr);
      if(nstr==0) return;
      if(!strcmp(*(data+0),"DISTS"))
      {
        vp->r  =strtod(*(data+1),NULL);
        vp->odv=strtod(*(data+2),NULL);
      }
      freestr(data,nstr);

	  data=fgetsbrk(fin,&nstr);
      if(nstr==0) return;
      if(!strcmp(*(data+0),"NICHI"))
      {
        nichi[0]=strtod(*(data+1),NULL);
        nichi[1]=strtod(*(data+2),NULL);
        nichi[2]=strtod(*(data+3),NULL);
      }
      freestr(data,nstr);

      data=fgetsbrk(fin,&nstr);
      if(nstr==0) return;
      if(!strcmp(*(data+0),"NFLOOR"))
      {
        nframe=(int)strtol(*(data+1),NULL,10);
      }
      freestr(data,nstr);

      for(k=0;k<nframe;k++)
      {
        (wdraw.childs+1)->vparam.subsection++;

        if((jj==0||jj==1) && k==0) (wdraw.childs+1)->vparam.vflag.nv.conffig=1;
        else                       (wdraw.childs+1)->vparam.vflag.nv.conffig=0;

        data=fgetsbrk(fin,&nstr);
        if(nstr==0) return;

        if(!strcmp(*(data+0),"FRAME"))
        {
            sprintf(fname,"%s",*(data+1));
        }
/*
        if(!strcmp(*(data+2),"HEAD"))
        {
          vp->range.min.d[GX]=strtod(*(data+3),NULL);
          vp->range.min.d[GY]=strtod(*(data+4),NULL);
          vp->range.min.d[GZ]=strtod(*(data+5),NULL);
        }
        if(!strcmp(*(data+6),"TAIL"))
        {
          vp->range.max.d[GX]=strtod(*(data+7),NULL);
          vp->range.max.d[GY]=strtod(*(data+8),NULL);
          vp->range.max.d[GZ]=strtod(*(data+9),NULL);
        }
*/
        /*Ai分布から各フロアレベル読み込み*/
        if(jj==1||jj==2)
        {
          vp->range.min.d[GZ]=*(((wdraw.childs+1)->org).ai.lbound+ k);
          vp->range.max.d[GZ]=*(((wdraw.childs+1)->org).ai.lbound+ k+1);
        }
        if(jj==0||jj==3)
        {
          vp->range.min.d[GZ]=*(((wdraw.childs+1)->org).ai.lbound+ k);
          vp->range.max.d[GZ]=*(((wdraw.childs+1)->org).ai.lbound+ k+2);
          if(k+2==nframe) k++;
        }
        SendMessage((wmenu.childs+2)->hwnd,WM_INITDIALOG,0,0);

        vp->focus.d[0]=nichi[0];
        vp->focus.d[1]=nichi[1];
        vp->focus.d[2]=nichi[2];
        SendMessage((wmenu.childs+2)->hwnd,WM_INITDIALOG,0,0);

        {
          sprintf(str,"%s%s",fname,s);
          SetDlgItemText((wmenu.childs+3)->hwnd,IDV_TITLEBOX,
                           str);

          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,wparam,0);
          if(i==0) SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                               MAKEWPARAM(IDV_DZ,0),0);
          if(i==1) SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                               MAKEWPARAM(IDV_DX,0),0);
          if(i==2) SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                               MAKEWPARAM(IDV_DY,0),0);

          SendMessage(wmain.hwnd,WM_COMMAND,
                      MAKEWPARAM((WORD)IDM_PRINTWAIT,(WORD)0),0);

          wparam = MAKEWPARAM(IDV_DEFORMATION,0);
          SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,wparam,0);
          if(i==0) SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
							   MAKEWPARAM(IDV_DZ,0),0);
          if(i==1) SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                               MAKEWPARAM(IDV_DX,0),0);
          if(i==2) SendMessage((wmenu.childs+3)->hwnd,WM_COMMAND,
                               MAKEWPARAM(IDV_DY,0),0);

		  /*END PAGE*/
          SendMessage(wmain.hwnd,WM_COMMAND,
                      MAKEWPARAM((WORD)IDM_PRINTOUT,(WORD)0),0);
        }
        vp->focus.d[0]=nichi[0];
        vp->focus.d[1]=nichi[1];
        vp->focus.d[2]=nichi[2];
        SendMessage((wmenu.childs+2)->hwnd,WM_INITDIALOG,0,0);

        freestr(data,nstr);
      }

    }
  }
  {
/*
    extractarclmfromorgan(&((wdraw.childs+1)->org),
    									 &arc,&arcx,&arcy);
*/

    (wmenu.childs+2)->vparam.vflag.mv.ftype=F_ARCLM;
    SetDlgItemText((wmenu.childs+2)->hwnd,ID_INPUTFILE,
				 (wdraw.childs+1)->inpfilez);
	wparam = MAKEWPARAM(IDD_VIEW,0);
    SendMessage((wmenu.childs+2)->hwnd,WM_COMMAND,wparam,0);

    (wdraw.childs+1)->vparam.vflag.mv.pagetitle=0;
    (wdraw.childs+1)->vparam.vflag.nv.conffig=1;
    vp->range.min.d[GX]=-100.0;
    vp->range.min.d[GY]=-100.0;
    vp->range.min.d[GZ]=-100.0;
    vp->range.max.d[GX]= 100.0;
    vp->range.max.d[GY]= 100.0;
    vp->range.max.d[GZ]= 100.0;
    vp->focus.d[0]     =initd[0];
    vp->focus.d[1]     =initd[1];
    vp->focus.d[2]     =initd[2];
    vp->theta          =inittheta;
    vp->phi            =initphi;
    vp->r              =initr;
    vp->odv            =initodv;
    (wdraw.childs+1)->vparam.vflag.ev.etype[COLUMN]=1;
    (wdraw.childs+1)->vparam.vflag.ev.etype[GIRDER]=1;
    (wdraw.childs+1)->vparam.vflag.ev.etype[BEAM]=1;
    (wdraw.childs+1)->vparam.vflag.ev.etype[BRACE]=1;
    (wdraw.childs+1)->vparam.vflag.ev.etype[WALL]=1;
    (wdraw.childs+1)->vparam.vflag.ev.etype[SLAB]=1;
    SendMessage((wmenu.childs+2)->hwnd,WM_INITDIALOG,0,0);

    clearwindow(*(wdraw.childs+1));
	drawarclmframe((wdraw.childs+1)->hdcC,
                   (wdraw.childs+1)->vparam,arc,0,ONSCREEN);
    SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);

    sprintf(str,"壁、床をブレース置換した結果");
    SetDlgItemText((wmenu.childs+3)->hwnd,IDV_TITLEBOX,str);
    (wdraw.childs+1)->vparam.subsection++;
    SendMessage(wmain.hwnd,WM_COMMAND,
                MAKEWPARAM((WORD)IDM_PRINTWAIT,(WORD)0),0);
    /*END PAGE*/
    SendMessage(wmain.hwnd,WM_COMMAND,
                MAKEWPARAM((WORD)IDM_PRINTOUT,(WORD)0),0);
  }


  if(!mode)
  {
  /*PRINT END.*/
  wparam = MAKEWPARAM((WORD)IDM_PRINTEND,(WORD)0);
  SendMessage(wmain.hwnd,WM_COMMAND,wparam,0);

  (wdraw.childs+1)->vparam.vflag.nv.code=0;
  (wdraw.childs+1)->vparam.vflag.ev.sectioncode=0;

  vp->range.min.d[GX]=-100.0;
  vp->range.min.d[GY]=-100.0;
  vp->range.min.d[GZ]=-100.0;
  vp->range.max.d[GX]= 100.0;
  vp->range.max.d[GY]= 100.0;
  vp->range.max.d[GZ]= 100.0;
  vp->focus.d[0]     =initd[0];
  vp->focus.d[1]     =initd[1];
  vp->focus.d[2]     =initd[2];
  vp->theta          =inittheta;
  vp->phi            =initphi;
  vp->r              =initr;
  vp->odv            =initodv;
  SendMessage((wmenu.childs+2)->hwnd,WM_INITDIALOG,0,0);

  (wmenu.childs+2)->vparam.vflag.mv.ftype=F_ORGAN;

  clearwindow(*(wdraw.childs+1));
  draworganization((wdraw.childs+1)->hdcC,
				   (wdraw.childs+1)->vparam,
				   (wdraw.childs+1)->org,ONSCREEN);
  SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);
  }
  return;
}/*printorganfigures*/

