/*CANVS101.C:CANVAS FOR WINDOWS95.SINCE:1996.2.17.JUNSATO.*/
/*LAST CHANGE:1997.5.28.*/
/*SKELETON FOR ARCLM001,101.*/

/*REGULATING.*/

/*FILE CONTENTS:*/
/*HOGAN001.IDE.....BORLAND PROJECT FILE.*/
/* CANHEAD.H  .....DEFINITION OF ID CODES.*/
/*  CANICO.RC .....DEFINITION OF ICON.*/
/*  CANCUR.RC .....DEFINITION OF CURSORS.*/
/* CANMENU.RC .....DEFINITION OF MENU.*/
/*  CANDLG.RC .....DEFINITION OF DIALOG BOX "INPUT".*/
/*  CANBMP.RC .....DEFINITION OF BITMAP.*/
/*  CANICO.ICO.....ICON OF CANVS.*/
/*  CANCUR.CUR.....CURSOR OF CANVS.NARROW ARROW.*/
/*  CANBOX.CUR.....CURSOR OF CANVS.SQUARE.*/
/*    VIEW.BMP.....BITMAP ON DIALOG BOX.FIGURE OF PERSPECTIVE.*/
/*CANVS101.C  .....SKELETON PROGRAM.*/
/*ARCHG116.C  .....SUBROUTINES FOR ANALYSIS MATERIALY NONLINEAR.*/
/*ARCHG013.C  .....SUBROUTINES FOR ANALYSIS LINEAR.*/
/*BCLNG014.C  .....SUBROUTINES FOR ANALYSIS ELASTIC BUCKLING.*/
/*QADHG001.C  .....SUBROUTINES FOR ANALYSIS LINEAR WITH FILM.*/

#include <windows.h>
#include <stdio.h>
#include <string.h>

#include "canhead.h"                    /*DEFINITION OF COMMAND ID.*/
#include "archg116.c"            /*ANALYSIS NONLINEAR,MATH,VIEWING.*/
#include "archg013.c"                     /*ANALYSIS STATIC LINEAR.*/
#include "bclng014.c"                           /*ELASTIC BUCKLING.*/
#include "qadhg001.c"               /*ANALYSIS BIQUADRATIC ELEMENT.*/

#define MENUDIALOGS 3 /*MENU DIALOG BOXES.*/
#define INPUTWIDTH 146              /*DIALOGUE FONT:"SYSTEM" 18pts.*/
#define BARWIDTH 10                          /*WIDTH OF SCROLL BAR.*/

/*#define A4HEIGHT 6371*/ /*2799*/         /*6371:MD4000J 2799:LPT8*/
/*#define A4WIDTH  4799*/ /*1970*/         /*4799:MD4000J 1970:LPT8*/

#define SECTSIMPLE 1
#define SECTDETAIL 2

struct print{int margin[4];
             int jiheight,jiwidth,jipitch;
             int gyopitch;
             int dans,dangap;};                    /*PRINT OPTIONS.*/

struct windowparams wmain={0,0,0,NEUTRAL,"CanMainWin",
                           NULL,NULL,NULL,NULL,
                           NEUTRAL,NEUTRAL,
                           NULL};
struct windowparams wmenu={1,0,0,NEUTRAL,"CanMenuWin",
                           NULL,NULL,NULL,NULL,
                           NEUTRAL,NEUTRAL,
                           NULL};
struct windowparams wdraw={2,0,0,NEUTRAL,"CanDrawWin",
                           NULL,NULL,NULL,NULL,
                           NEUTRAL,NEUTRAL,
                           NULL,0,0};
struct windowparams wmesg={3,0,0,NEUTRAL,"CanMesgWin",
                           NULL,NULL,NULL,NULL,
                           NEUTRAL,NEUTRAL,
                           NULL,0,0};
struct windowparams wsurf={4,0,0,NEUTRAL,"CanSurfWin",
                           NULL,NULL,NULL,NULL,
                           NEUTRAL,NEUTRAL,
                           NULL,0,0};

struct windowparams wsect={5,0,0,NEUTRAL,"CanSectWin",
                           NULL,NULL,NULL,NULL,
                           NEUTRAL,NEUTRAL,
                           NULL,0,0};
struct windowparams wweig={6,0,0,NEUTRAL,"CanWeightWin",
                           NULL,NULL,NULL,NULL,
                           NEUTRAL,NEUTRAL,
                           NULL,0,0};
struct windowparams wcmqd={7,0,0,NEUTRAL,"CanCmqWin",
                           NULL,NULL,NULL,NULL,
                           NEUTRAL,NEUTRAL,
                           NULL,0,0};
struct windowparams whorz={8,0,0,NEUTRAL,"CanHorizonWin",
                           NULL,NULL,NULL,NULL,
                           NEUTRAL,NEUTRAL,
                           NULL,0,0};

struct arclmframe arc={0,0,
                       "\0",
                       0,0,0,0,0,
                       NULL,NULL,
                       NULL,NULL,NULL,NULL,
                       NULL,NULL,
                       NULL,NULL,NULL}; /*GLOBAL ARCLM FRAME.*/

LRESULT CALLBACK WindowProcedureMain(HWND,UINT,WPARAM,LPARAM);
LRESULT CALLBACK WindowProcedureSheet(HWND,UINT,WPARAM,LPARAM);
LRESULT CALLBACK WindowProcedureBack(HWND,UINT,WPARAM,LPARAM);
LRESULT CALLBACK WindowProcedureDraw(HWND,UINT,WPARAM,LPARAM);
LRESULT CALLBACK WindowProcedureSurf(HWND,UINT,WPARAM,LPARAM);
LRESULT CALLBACK WindowProcedureSect(HWND,UINT,WPARAM,LPARAM);
LRESULT CALLBACK WindowProcedureMesg(HWND,UINT,WPARAM,LPARAM);
LRESULT CALLBACK WindowProcedureMenu(HWND,UINT,WPARAM,LPARAM);
BOOL CALLBACK DialogProcMenu1(HWND,UINT,WPARAM,LPARAM);
BOOL CALLBACK DialogProcMenu2(HWND,UINT,WPARAM,LPARAM);
BOOL CALLBACK DialogProcMenu3(HWND,UINT,WPARAM,LPARAM);
BOOL CALLBACK DialogProcText(HWND,UINT,WPARAM,LPARAM);
BOOL CALLBACK DialogProcSsim(HWND,UINT,WPARAM,LPARAM);
BOOL CALLBACK EnumChildProcSheet(HWND hwnd,LPARAM lParam);
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
struct snode *addtext(struct snode *strset,int *nstr,char *str,
                      double tx,double ty);

HINSTANCE hInstGlobal;                                  /*INSTANCE.*/
struct winparamsreg wrglobal={0,0,NULL};   /*REGISTED WINDOWPARAMS.*/

int globalmode=NEUTRAL;             /*MODES:ARCLM,SRCAL,HOGAN,ORGAN*/
int globalstatus=NEUTRAL;
int prestatus=NEUTRAL;
int createcommand=C_NEUTRAL; /*FRAME CREATION COMMAND.*/
int createitem=C_NEUTRAL; /*FRAME CREATION TARGET ITEM.*/

POINT pbar;                                /*RELATIVE POINT ON BAR.*/

long int initx,inity;       /*INITIAL POSITION OF ROTATE,MOVE,TEXT.*/
struct viewparam initv;

HWND hpopdlg=NULL;                                  /*POPUP DIALOG.*/
int idtext;                                   /*ID OF CURRENT TEXT.*/

char gstr[256]; /*STRING BRINGING SECTION CODE.*/
long int gnsect=0; /*GLOBAL NSECT.*/
struct osect *currentsects,*gsect;

struct onode nlast; /*LAST NODE OF ELEMENT WHILE CREATING.*/
struct oelem gelem={0,0,
                    ROLEWEIGHT,WALL,
                    0,0,0.0,
                    {{0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                     {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}},
                    {{0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                     {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}},
                    NULL,NULL,NULL,NULL,
                    {{255,255,255},{255,255,255},{255,255,255}}
                   }; /*CREATING ELEMENT.*/

int WINAPI WinMain(HINSTANCE hInstance,
                   HINSTANCE hPrevInstance,
                   LPSTR lpszCmdLine,
                   int nCmdShow)
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

  hwnd = CreateWindow(wmain.classname,               /*WINDOW CLASS*/
                      "CANVS SKELETON",                   /*CAPTION*/
                      WS_CLIPCHILDREN |
                      WS_CLIPSIBLINGS |
                      WS_OVERLAPPEDWINDOW,                  /*STYLE*/
                      50,10,                                  /*X,Y*/
                      650,700,                       /*WIDTH,HEIGHT*/
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
    if(wmenu.nchilds<5 || !IsDialogMessage(hpopdlg,&msg))
    {
      isd=1;
      for(nc=2;nc<(MENUDIALOGS+2);nc++)
      {
        if(IsDialogMessage((wmenu.childs+nc)->hwnd,&msg)) isd=0;
      }
      for(nc=wsect.nchilds;nc>=3;nc--)
      {
        if(IsDialogMessage((wsect.childs+nc-1)->hwnd,&msg)) isd=0;
      }

      if(isd==1 && !TranslateAccelerator(hwnd,hAccel,&msg))
      {
        TranslateMessage(&msg);               /*KEY INTO CHARACTER.*/
        DispatchMessage(&msg);                  /*DISPATCH MESSAGE.*/
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
  /*HFONT hfont;*/
  WORD x,y;                                   /*WORD=unsigned short*/
  WPARAM wparam;
  FILE *fout;
  char str[256];
  int i;
  long int maxX,maxY,mw1,mh1,mw2,mh2,mw3,mh3;  /*WINDOW,CLIENT SIZE*/

  switch(message)
  {
    case WM_PAINT:
      DefWindowProc(hwnd,message,wParam,lParam);
      if(wmenu.hwnd==NULL)
      {
        wparam = MAKEWPARAM((WORD)IDM_OPTIONS,(WORD)0);
        SendMessage(wmain.hwnd,WM_COMMAND,wparam,0);
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

        case IDM_DRAWINGS:
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
                             mw1,0,maxX-mw1,maxY,
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
                             1,1,maxX-BARWIDTH-5,maxY-BARWIDTH-5,
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
                             0,0,500,800,
                             (wdraw.childs+0)->hwnd);
            setfontformat((wdraw.childs+1)->hdcC,20,5,
                          "Terminal",0,0,255);

            initializeorganization(&((wdraw.childs+1)->org));

            sprintf((wdraw.childs+1)->inpfile,
                    "c:\\cdocs\\hogan\\organ10.inp");
            sprintf((wdraw.childs+1)->otpfile,
                    "hogtxt.otp");

            (wdraw.childs+1)->vparam.type=PERSPECTIVE;
            getviewparam((wmenu.childs+2)->hwnd,
                         &((wdraw.childs+1)->vparam));

            getclientsize((wdraw.childs+1)->hwnd,&maxX,&maxY);
            (wdraw.childs+1)->vparam.Xo=(int)(maxX/2);
            (wdraw.childs+1)->vparam.Yo=(int)(maxY/2);

            (wdraw.childs+1)->vparam.vflag.axis=1;

            (wdraw.childs+1)->vparam.vflag.nv.code=1;
            (wdraw.childs+1)->vparam.vflag.nv.confs[0]=0;
            (wdraw.childs+1)->vparam.vflag.nv.confs[1]=0;
            (wdraw.childs+1)->vparam.vflag.nv.confs[2]=0;
            (wdraw.childs+1)->vparam.vflag.nv.confs[3]=0;
            (wdraw.childs+1)->vparam.vflag.nv.confs[4]=0;
            (wdraw.childs+1)->vparam.vflag.nv.confs[5]=0;

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

            (wdraw.childs+1)->vparam.vflag.ev.code=0;
            (wdraw.childs+1)->vparam.vflag.ev.sectioncode=0;
            (wdraw.childs+1)->vparam.vflag.ev.axis=0;
            (wdraw.childs+1)->vparam.vflag.ev.hinge=1;
            (wdraw.childs+1)->vparam.vflag.ev.deformation=0;
            (wdraw.childs+1)->vparam.vflag.ev.stress[0]=0;
            (wdraw.childs+1)->vparam.vflag.ev.stress[1]=0;
            (wdraw.childs+1)->vparam.vflag.ev.stress[2]=0;
            (wdraw.childs+1)->vparam.vflag.ev.stress[3]=0;
            (wdraw.childs+1)->vparam.vflag.ev.stress[4]=0;
            (wdraw.childs+1)->vparam.vflag.ev.stress[5]=0;

            (wmenu.childs+2)->vparam.vflag.mv.draw=1;
            SendMessage((wmenu.childs+2)->hwnd,WM_INITDIALOG,0,0);
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
                             700,100,200,500,
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

            setfontformat((wmesg.childs+1)->hdcC,20,5,
                          "Terminal",255,255,255);

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
                             600,400,200,300,
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
            setfontformat((wsurf.childs+1)->hdcC,20,5,
                          "Terminal",100,255,255);

            getclientsize((wsurf.childs+1)->hwnd,&maxX,&maxY);
            (wsurf.childs+1)->vparam.type=PERSPECTIVE;
            (wsurf.childs+1)->vparam.gfactor=1.0;
            (wsurf.childs+1)->vparam.focus.d[0]=0.0;
            (wsurf.childs+1)->vparam.focus.d[1]=0.0;
            (wsurf.childs+1)->vparam.focus.d[2]=0.0;
            (wsurf.childs+1)->vparam.theta=30.0;
            (wsurf.childs+1)->vparam.phi=20.0;
            (wsurf.childs+1)->vparam.r=5.0;
            (wsurf.childs+1)->vparam.odv=400.0;
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

        case IDM_SECTIONLIST:
          if(wsect.hwnd==NULL &&
             wmenu.hwnd!=NULL) /*CREATION.*/
          {
            windowdefinition(&wsect,
                             hInstGlobal,
                             NULL,
                             WindowProcedureSheet,
                             "HOGSECT",
                             190,190,190,
                             WS_CLIPCHILDREN |
                             WS_CLIPSIBLINGS |
                             /*WS_CHILD |*/
                             WS_POPUP |
                             WS_CAPTION |
                             WS_THICKFRAME |
                             WS_VISIBLE,
                             500,300,100,300,
                             wmain.hwnd);

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
                             0,0,500,500,
                             (wsect.childs+0)->hwnd);
            setfontformat((wsect.childs+1)->hdcC,20,5,
                          "Terminal",100,100,170);

            (wsect.childs+1)->vparam.Xo=0;
            (wsect.childs+1)->vparam.Yo=0;

            (wmenu.childs+2)->vparam.vflag.mv.sectlist=1;
            SendMessage((wmenu.childs+2)->hwnd,WM_INITDIALOG,0,0);
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

            (wmenu.childs+2)->vparam.vflag.mv.sectlist=0;
            SendMessage((wmenu.childs+2)->hwnd,WM_INITDIALOG,0,0);
          }
          break;

        case IDM_OPTIONS:
          if(wmenu.hwnd==NULL)
          {
            getclientsize(hwnd,&maxX,&maxY);

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

            (wmenu.childs+2)->vparam.vflag.mv.ftype   =F_ORGAN;
            (wmenu.childs+2)->vparam.vflag.mv.savetype=F_ORGAN;

            (wmenu.childs+2)->hwnd
            = CreateDialog(hInstGlobal,
                           "CANOPTIONS1",
                           (wmenu.childs+1)->hwnd,
                           DialogProcMenu1);

            (wmenu.childs+3)->hwnd
            = CreateDialog(hInstGlobal,
                           "CANOPTIONS2",
                           (wmenu.childs+1)->hwnd,
                           DialogProcMenu2);

            (wmenu.childs+4)->hwnd
            = CreateDialog(hInstGlobal,
                           "CANOPTIONS3",
                           (wmenu.childs+1)->hwnd,
                           DialogProcMenu3);

            getclientsize((wmenu.childs+2)->hwnd,&mw1,&mh1);
            getclientsize((wmenu.childs+3)->hwnd,&mw2,&mh2);
            getclientsize((wmenu.childs+4)->hwnd,&mw3,&mh3);
            MoveWindow((wmenu.childs+1)->hwnd,
                       0,0,mw1,(mh1+mh2+mh3),TRUE);
            MoveWindow((wmenu.childs+3)->hwnd,
                       0,(int)(DPT_HEIGHT1*2.75)/*1702*/,
                       mw2,mh2,TRUE);
            MoveWindow((wmenu.childs+4)->hwnd,
                       0,(int)((DPT_HEIGHT1+DPT_HEIGHT2)*2.75),
                       mw3,mh3,TRUE);

            SetDlgItemText((wmenu.childs+3)->hwnd,IDV_MODENUM,"1");
          }
          break;

        case IDM_SAVEORGAN:
          if(wdraw.nchilds>=2 &&
             (wdraw.childs+1)->hwnd!=NULL &&
             (wmenu.childs+2)->vparam.vflag.mv.ftype==F_ORGAN &&
             (wdraw.childs+1)->org.nodes!=NULL)
          {
            fout=fopen((wdraw.childs+1)->otpfile,"w");
            if(fout==NULL) break;

            saveorganization(fout,&((wdraw.childs+1)->org));

            fclose(fout);

            sprintf(str,"Saved As \"%s\"",
                    (wdraw.childs+1)->otpfile);
            MessageBox(NULL,str,"Save Organ Frame",MB_OK);
          }
          break;

        case IDM_HELP:                                 /*HELP TEXT.*/
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
  WPARAM wparam;
  LPARAM lparam;
  int bl,bt,br,bb,bmax; /*LEFT,TOP,RIGHT,BOTTOM OF BAR.*/
  long int cw,ch,pw,ph,maxX,maxY,wintop,winleft; /*WINDOW SIZE.*/
  struct windowparams *wp;

  switch(message)
  {
    case WM_PAINT:
    case WM_SIZE:
      DefWindowProc(hwnd,message,wParam,lParam);
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
/*MessageBox(NULL,"Bar End.","Sheet",MB_OK);*/
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

    case WM_COMMAND: /*SEND MENU MESSAGE TO CHILDREN.*/
      switch(LOWORD(wParam))
      {
        case IDMS_OPEN:
          wparam = MAKEWPARAM((WORD)IDMS_OPEN,(WORD)0);
          lparam = MAKELPARAM((WORD)0,(WORD)0);
          SendMessage((wsect.childs+1)->hwnd,
                      WM_COMMAND,wparam,lparam);
          break;
        case IDMS_SIMPLE:
          wparam = MAKEWPARAM((WORD)IDMS_SIMPLE,(WORD)0);
          SendMessage((wsect.childs+1)->hwnd,WM_COMMAND,wparam,0);
          break;
        case IDMS_DETAIL:
          wparam = MAKEWPARAM((WORD)IDMS_DETAIL,(WORD)0);
          SendMessage((wsect.childs+1)->hwnd,WM_COMMAND,wparam,0);
          break;

        default: /*OTHERS.*/
          return DefWindowProc(hwnd,message,wParam,lParam);
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
  POINT point;
  WPARAM wparam;
  /*MSG msg;*/
  /*HACCEL haccel;*/
  char str[256];
  int i,direction;
  long int maxX,maxY,x,y,mw,mh;
  long int code;
  double length,dx,dy;
  struct onode *node;
  struct owire *elem;
  struct plane pban;

  switch(message)
  {
    case WM_PAINT:
      DefWindowProc(hwnd,message,wParam,lParam);

      /*GET TEXTS IN ALL ITEMS OF OPTION DIALOG.....NOT YET.*/

      overlayhdc(*(wdraw.childs+1),SRCPAINT);

      hdc = GetDC(hwnd);
      setfontformat(hdc,20,5,"Terminal",0,255,255);
      drawtexts(hdc,wdraw.strset,wdraw.nstring,
                (wdraw.childs+1)->vparam);
      break;

    case WM_LBUTTONDOWN:
      x = LOWORD(lParam);
      y = HIWORD(lParam);
      point.x=x;
      point.y=y;

      /*WHILE CREATING FRAME.*/
      if((wmenu.childs+2)->vparam.vflag.mv.ftype==F_ORGAN)
      {
        if(createcommand==C_ADD)
        {
          if(createitem==C_NODE)
          {
            createorgannode(0,(wdraw.childs+1)->vparam,
                            x,y,&((wdraw.childs+1)->org));
          }
          if(createitem==C_ELEMENT)
          {
            if(gelem.nnod==0)
            {
              gelem.nnod=1;
            }
            else if(gelem.nnod==2)
            {
              gelem.nban=1;
              gelem.bans=(struct obans *)
                         malloc(sizeof(struct obans));
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
              nlast=findlastnode(0,(wdraw.childs+1)->vparam,
                                 x,y,
                                 NULL,
                                 &((wdraw.childs+1)->org));
              if(nlast.code==0)
              {
                MessageBox(NULL,"Find Node Failed.","Adding",MB_OK);
                break;
              }
            }
            else if(gelem.nban>0)
            {
              bantoplane(*(gelem.bans+0),&pban);

              nlast=findlastnode(0,(wdraw.childs+1)->vparam,
                                 x,y,
                                 &pban,
                                 &((wdraw.childs+1)->org));
              if(nlast.code==0)
              {
                MessageBox(NULL,"Find Node Failed.","Adding",MB_OK);
                break;
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

            *(gelem.nods+gelem.nnod-2)
            =addnode(nlast,&((wdraw.childs+1)->org));
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
          }

          draworganization((wdraw.childs+1)->hdcC,
                           (wdraw.childs+1)->vparam,
                           (wdraw.childs+1)->org,ONSCREEN);
          SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);
          break;
        }
      }

      /*CASE NOT CREATING.*/
      prestatus=NEUTRAL;

      if(globalstatus==SELECTNODE)
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

      /*END OF ADDING ELEMENT.*/
      if((wmenu.childs+2)->vparam.vflag.mv.ftype==F_ORGAN)
      {
        if(createcommand==C_ADD)
        {
          if(createitem==C_ELEMENT && gelem.nnod>0)
          {
            /*gelem.code= ? ;*/
            gelem.loff=(wdraw.childs+1)->org.nelem; /*PUT LAST.*/
            gelem.sect=(wdraw.childs+1)->org.sects+0;

            gelem.nnod--;
            gelem.nods=(struct onode **)
                       realloc(gelem.nods,
                               gelem.nnod*sizeof(struct onode *));
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
                =(struct onode **)realloc((gelem.bans+0)->nods,
                                          (gelem.bans+0)->nnod
                                          *sizeof(struct onode *));
                if((gelem.bans+0)->nods==NULL)
                {
                  MessageBox(NULL,"Buffer Null.","Adding",MB_OK);
                }
              }
            }

/*MessageBox(NULL,"Pass 1","Adding Element",MB_OK);*/
            addelement(&gelem,&((wdraw.childs+1)->org));

            /*freeelement(&gelem,0,0,1,0);*/
/*MessageBox(NULL,"Pass 2","Adding Element",MB_OK);*/
            gelem.nnod=0;
            gelem.nban=0;
          }

          draworganization((wdraw.childs+1)->hdcC,
                           (wdraw.childs+1)->vparam,
                           (wdraw.childs+1)->org,ONSCREEN);
          SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);
          break;
        }
      }

      /*CASE NOT CREATING.*/
      popupmenudraw(hwnd,point);
      break;

    case WM_MOUSEMOVE:
      /*WHILE ADDING ELEMENT.*/
      if((wmenu.childs+2)->vparam.vflag.mv.ftype==F_ORGAN)
      {
        if(createcommand==C_ADD)
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
      }

      /*CASE NOT CREATING.*/
      if(globalstatus==ROTATE &&
         (wdraw.childs+1)->lstatus==ROTATE &&
         (wParam == MK_LBUTTON))
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
        if((wmenu.childs+2)->vparam.vflag.mv.ftype==F_ARCLM)
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
      else if(globalstatus==MOVE &&
              (wdraw.childs+1)->lstatus==MOVE &&
              (wParam == MK_LBUTTON))
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
        if((wmenu.childs+2)->vparam.vflag.mv.ftype==F_ARCLM)
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
        if((wmenu.childs+2)->vparam.vflag.mv.ftype==F_ARCLM)
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

        case IDM_POPTEXT:
          hpopdlg=CreateDialog(hInstGlobal,
                               "CANTEXTDLG",
                               NULL,
                               DialogProcText);
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
                         (wsurf.childs+1)->vparam,2,4,0,NULL);
        overlayhdc(*(wsurf.childs+1),SRCPAINT);   /*UPDATE DISPLAY.*/
      }
      else if(globalstatus==ROTATE &&
              (wsurf.childs+1)->lstatus==ROTATE) /*END OF ROTATION.*/
      {
        clearwindow(*(wsurf.childs+1));
        drawyieldsurface((wsurf.childs+1)->hdcC,
                         (wsurf.childs+1)->vparam,2,4,0,
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

LRESULT CALLBACK WindowProcedureSect(HWND hwnd,
                                     UINT message,
                                     WPARAM wParam,
                                     LPARAM lParam)
/*WINDOW FOR SECTION LIST.*/
{
  HWND hdlg;
  WPARAM wparam;
  LPARAM lparam;
  SIZE tsize;
  char str[256];
  long int i,j;
  long int dw,dh,wmax,hmax;

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
            if((wmenu.childs+2)->vparam.vflag.mv.ftype==F_ARCLM)
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
            wmax=0;
            hmax=0;
            for(i=0;i<gnsect;i++)
            {
              gsect=currentsects+i;
              /*gsect->dflag=1;*/ /*DRAW FLAG ON.*/
              /*gsect->dcolor.r=255;
              gsect->dcolor.g=255;
              gsect->dcolor.b=255;*/
              sprintf(gstr,"%ld",gsect->code);

              hdlg = CreateDialog(hInstGlobal,
                                  "HOGDLGSECTSIMPLE",
                                  (wsect.childs+1)->hwnd,
                                  DialogProcSsim);

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
                  if(gsect->name!=NULL)
                  {
                    GetTextExtentPoint32((wsect.childs+1)->hdcC,
                                         gsect->name,
                                         strlen(gsect->name),&tsize);
                    wsect.strset=addtext(wsect.strset,
                                         &(wsect.nstring),
                                         gsect->name,
                                         50.0,(double)hmax);
                    hmax+=tsize.cy;
                  }
                  for(j=0;j<(gsect->nfig);j++)
                  {
                    if((gsect->figs+j)->prop->name!=NULL &&
                       (gsect->figs+j)->thick>0.0)
                    {
                      sprintf(str," %s %.0fmm",
                              (gsect->figs+j)->prop->name,
                              (gsect->figs+j)->thick*1000.0);

                      GetTextExtentPoint32((wsect.childs+1)->hdcC,
                                           str,strlen(str),&tsize);
                      wsect.strset=addtext(wsect.strset,
                                           &(wsect.nstring),
                                           str,50.0,(double)hmax);
                      hmax+=tsize.cy;
                    }
                  }
                }
              }
            }
            MoveWindow((wsect.childs+1)->hwnd,0,0,wmax,hmax,TRUE);
          }
          break;

        case IDMS_SIMPLE: /*SIMPLE TYPE.*/
          wparam = MAKEWPARAM((WORD)IDMS_OPEN,(WORD)0);
          lparam = MAKELPARAM((WORD)SECTSIMPLE,(WORD)0);
          SendMessage((wsect.childs+1)->hwnd,
                      WM_COMMAND,wparam,lparam);
          break;

        case IDMS_DETAIL: /*DETAILED TYPE.*/
          wparam = MAKEWPARAM((WORD)IDMS_OPEN,(WORD)0);
          lparam = MAKELPARAM((WORD)SECTDETAIL,(WORD)0);
          SendMessage((wsect.childs+1)->hwnd,
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
}/*WindowProcedureSect*/

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

BOOL CALLBACK DialogProcMenu1(HWND hdwnd,
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
  int id;
  /*long int maxX,maxY,mw,mh;*/
  FILE *fin;
  struct globvisible gv;

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

      SetDlgItemText(hdwnd,IDV_TXTAXONO,"Axonometric");
      SetDlgItemText(hdwnd,IDV_TXTPERS, "Perspective");

      if(wdraw.nchilds>=2 && (wdraw.childs+1)!=NULL)
      {
        setviewparam(hdwnd,(wdraw.childs+1)->vparam);

        SetDlgItemText(hdwnd,ID_INPUTFILE,
                       (wdraw.childs+1)->inpfile);
        SetDlgItemText(hdwnd,ID_OUTPUTFILE,
                       (wdraw.childs+1)->otpfile);
      }
      else
      {
        SetDlgItemText(hdwnd,IDV_GFACTOR,"1.0");

        SetDlgItemText(hdwnd,IDV_X,"0.0");
        SetDlgItemText(hdwnd,IDV_Y,"0.0");
        SetDlgItemText(hdwnd,IDV_Z,"0.0");

        SetDlgItemText(hdwnd,IDV_PHI,"50.0");
        SetDlgItemText(hdwnd,IDV_THETA,"-40.0");
        SetDlgItemText(hdwnd,IDV_R,"50.0");
        SetDlgItemText(hdwnd,IDV_L,"1000.0");
      }

      SetDlgItemText(hdwnd,ID_LAPS,"1");
      SetDlgItemText(hdwnd,ID_SAFETY,"1.0");

      SetDlgItemText(hdwnd,IDV_GAXISLENGTH,"7.0");
      SetDlgItemText(hdwnd,IDV_EAXISLENGTH,"0.5");
      SetDlgItemText(hdwnd,IDV_DFACTOR,    "1000.0");
      SetDlgItemText(hdwnd,IDV_QFACTOR,    "0.05");
      SetDlgItemText(hdwnd,IDV_MFACTOR,    "0.05");
      SetDlgItemText(hdwnd,IDV_GYOPITCH,   "0.5");
      SetDlgItemText(hdwnd,IDV_HINGESIZE,  "5");

      SetDlgItemText(hdwnd,IDR_XMAX,"100.0");
      SetDlgItemText(hdwnd,IDR_XMIN,"-100.0");
      SetDlgItemText(hdwnd,IDR_YMAX,"100.0");
      SetDlgItemText(hdwnd,IDR_YMIN,"-100.0");
      SetDlgItemText(hdwnd,IDR_ZMAX,"100.0");
      SetDlgItemText(hdwnd,IDR_ZMIN,"-100.0");
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
         (id==IDF_TXTFRAME && gv.mv.ftype!=F_FRAME))
      {
        SetTextColor(hdc,RGB(255,255,255));
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

      if(wdraw.hwnd!=NULL) gv=(wdraw.childs+1)->vparam.vflag;

      if((id==IDV_TXTGLOBALAXIS) &&
         (wdraw.hwnd==NULL || gv.axis!=1))
      {
        SetTextColor(hdc,RGB(255,255,255));
      }
      if((id==IDV_TXTNODECODE) &&
         (wdraw.hwnd==NULL || gv.nv.code!=1))
      {
        SetTextColor(hdc,RGB(255,255,255));
      }
      return (BOOL)(HBRUSH)GetStockObject(LTGRAY_BRUSH);

    case WM_PAINT: /*INSERT ILLUSTRATION.*/
      /*DefWindowProc(hwnd,message,wParam,lParam);*/

      hdc=GetDC(hdwnd);
      hdcC = CreateCompatibleDC(hdc);

      hbit = LoadBitmap(hInstGlobal,"BITMAPVIEW");
      GetObject(hbit,sizeof(BITMAP),(LPSTR)&bmp);
      pbit = SelectObject(hdcC,hbit);
      BitBlt(hdc,1,(int)((DPT_VIEW+26)*2.75)/*836*/,
                 bmp.bmWidth,bmp.bmHeight,
                 hdcC,0,0,SRCCOPY);

      SelectObject(hdcC,pbit);
      DeleteObject(hbit);
      ReleaseDC(hdwnd,hdcC);
      ReleaseDC(hdwnd,hdc);
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

        case ID_INPUTFILE: /*GET INPUT OF KEYBOARD ONE BY ONE.*/
          if(wdraw.nchilds>=2 && (wdraw.childs+1)!=NULL)
          {
            GetDlgItemText(hdwnd,ID_INPUTFILE,
                           (wdraw.childs+1)->inpfile,80);
          }
          break;
        case ID_OUTPUTFILE:
          if(wdraw.nchilds>=2 && (wdraw.childs+1)!=NULL)
          {
            GetDlgItemText(hdwnd,ID_OUTPUTFILE,
                           (wdraw.childs+1)->otpfile,80);
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
          wparam = MAKEWPARAM((WORD)IDM_SECTIONLIST,(WORD)0);
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
          if(MessageBox(NULL,"Arclm001:Begin.","ARCLM001",
             MB_OKCANCEL)==IDCANCEL) break;

          if((wmenu.childs+2)->vparam.vflag.mv.ftype==F_ARCLM)
          {
            getviewparam((wmenu.childs+2)->hwnd,
                         &((wdraw.childs+1)->vparam));
            clearwindow(*(wdraw.childs+1));
            arclm001(&arc);
          }
          break;
        case IDD_ARCLM101:
          if(MessageBox(NULL,"Arclm101:Begin.","ARCLM101",
             MB_OKCANCEL)==IDCANCEL) break;

          if((wmenu.childs+2)->vparam.vflag.mv.ftype==F_ARCLM)
          {
            getviewparam((wmenu.childs+2)->hwnd,
                         &((wdraw.childs+1)->vparam));
            clearwindow(*(wdraw.childs+1));
            arclm101(&arc);
          }
          break;
        case IDD_BCLNG001:
          if(MessageBox(NULL,"Bclng001:Begin.","BCLNG001",
             MB_OKCANCEL)==IDCANCEL) break;

          if((wmenu.childs+2)->vparam.vflag.mv.ftype==F_ARCLM)
          {
            if(arc.elems==NULL) break;
            getviewparam((wmenu.childs+2)->hwnd,
                         &((wdraw.childs+1)->vparam));
            /*clearwindow(*(wdraw.childs+1));*/
            bclng001(&arc);
          }
          break;

        case IDD_VIEW:
          getviewparam((wmenu.childs+2)->hwnd,
                       &((wdraw.childs+1)->vparam));

          if(wdraw.nchilds>=2 &&
             (wdraw.childs+1)->hwnd!=NULL &&
             (wmenu.childs+2)->vparam.vflag.mv.ftype==F_ARCLM)
          {
            fin=fgetstofopen(dir,"r",ID_INPUTFILE);    /*OPEN FILE.*/
            if(fin==NULL) break;

            inputinit(fin,&(arc.nnode),&(arc.nelem),&(arc.nsect));
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

            inputtexttomemory(fin,&arc);
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

            inputorganization(fin,&((wdraw.childs+1)->org));
            fclose(fin);

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
                           (wmenu.childs+3)->hwnd,
                           (wdraw.childs+1)->vparam,
                           fin);
            SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);

            fclose(fin);
          }
          break;

        case IDD_ROTATE: /*ROTATION*/
          if(wdraw.hwnd!=NULL)
          {
            (wdraw.childs+1)->lstatus=ROTATE;
            hcursor=LoadCursor(hInstGlobal,"CANCURSOR");   /*CURSOR*/
            SetClassLong((wdraw.childs+1)->hwnd,
                         GCL_HCURSOR,(LONG)hcursor);
            globalstatus=NEUTRAL;
          }
          break;
        case IDD_MOVE:                                /*MOVE FRAME.*/
          if(wdraw.hwnd!=NULL)
          {
            (wdraw.childs+1)->lstatus=MOVE;
            hcursor = LoadCursor(hInstGlobal,"CANBOX");    /*CURSOR*/
            SetClassLong((wdraw.childs+1)->hwnd,
                         GCL_HCURSOR,(LONG)hcursor);
            globalstatus=NEUTRAL;
          }
          break;
      }
      break;

    default:
      return DefWindowProc(hdwnd,message,wParam,lParam);
  }
  return 0;
}/*DialogProcMenu1*/

BOOL CALLBACK DialogProcMenu2(HWND hdwnd,
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
  /*HBRUSH hbrush,pbrush;*/
  /*LPDRAWITEMSTRUCT lpdis;*/
  PRINTDLG pd;
  DOCINFO di;
  char str[80],doc[80];
  char dir[]=DIRECTORY;
  int id;
  int cWidthPels,cHeightPels,caps;
  long int code,mode;
  struct print prn;
  struct globvisible gv;
  struct viewparam vprint;
  FILE *fout;

  switch(message)
  {
    case WM_INITDIALOG:
      SetDlgItemText(hdwnd,IDV_TXTGLOBALAXIS, "GlobalAxis");
      SetDlgItemText(hdwnd,IDV_TXTNODECODE,   "NodeCode");
      SetDlgItemText(hdwnd,IDV_TXTLOADS,      "UnitLoads");
      SetDlgItemText(hdwnd,IDV_TXTCONFINEMENT,"Confinement");
      SetDlgItemText(hdwnd,IDV_TXTELEMENTCODE,"ElementCode");
      SetDlgItemText(hdwnd,IDV_TXTELEMENTAXIS,"ElementAxis");
      SetDlgItemText(hdwnd,IDV_TXTHINGE,      "Hinge");
      SetDlgItemText(hdwnd,IDV_TXTSECTIONCODE,"SectionCode");
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
      SetDlgItemText(hdwnd,IDV_TXTREACTION,"Reaction");

      SetDlgItemText(hdwnd,IDS_CODE,"0");

      SetDlgItemText(hdwnd,IDN_CODE,"0");

      SetDlgItemText(hdwnd,IDE_CODE,"0");

      SetDlgItemText(hdwnd,IDP_MARGINTOP,     "0");
      SetDlgItemText(hdwnd,IDP_MARGINBOTTOM,  "0");
      SetDlgItemText(hdwnd,IDP_MARGINLEFT,  "300");
      SetDlgItemText(hdwnd,IDP_MARGINRIGHT, "200");
      SetDlgItemText(hdwnd,IDP_JIHEIGHT, "80");
      SetDlgItemText(hdwnd,IDP_JIWIDTH,  "30");
      SetDlgItemText(hdwnd,IDP_JIPITCH,   "0");
      SetDlgItemText(hdwnd,IDP_GYOPITCH, "70");
      SetDlgItemText(hdwnd,IDP_DANS,      "2");
      SetDlgItemText(hdwnd,IDP_DANGAP,  "100");

      SetDlgItemText(hdwnd,IDH_ASIZE,"10.0");
      SetDlgItemText(hdwnd,IDH_NSIZE,"5.0");
      SetDlgItemText(hdwnd,IDH_TITLEX,"X");
      SetDlgItemText(hdwnd,IDH_TITLEY,"Y");
      SetDlgItemText(hdwnd,IDH_SCALEX,"1.0");
      SetDlgItemText(hdwnd,IDH_SCALEY,"1.0");
      SetDlgItemText(hdwnd,IDH_XMAX,"100.0");
      SetDlgItemText(hdwnd,IDH_YMAX,"100.0");
      SetDlgItemText(hdwnd,IDH_XMIN,"-100.0");
      SetDlgItemText(hdwnd,IDH_YMIN,"-100.0");
      SetDlgItemText(hdwnd,IDH_PITCHX,"10.0");
      SetDlgItemText(hdwnd,IDH_PITCHY,"10.0");
      break;

    case WM_CTLCOLORSTATIC: /*STATIC TEXT MENU.*/
      DefWindowProc(hdwnd,message,wParam,lParam);

      hdc=(HDC)wParam;
      hitem=(HWND)lParam;
      id=GetDlgCtrlID(hitem);

      if(wdraw.hwnd!=NULL) gv=(wdraw.childs+1)->vparam.vflag;

      if((id==IDV_TXTGLOBALAXIS) &&
         (wdraw.hwnd==NULL || gv.axis!=1))
      {
        SetTextColor(hdc,RGB(255,255,255));
      }
      if((id==IDV_TXTNODECODE) &&
         (wdraw.hwnd==NULL || gv.nv.code!=1))
      {
        SetTextColor(hdc,RGB(255,255,255));
      }
      if((id==IDV_TXTLOADS) &&
         (wdraw.hwnd==NULL || gv.nv.loads[0]!=1))
      {
        SetTextColor(hdc,RGB(255,255,255));
      }
      if((id==IDV_TXTCONFINEMENT) &&
         (wdraw.hwnd==NULL || gv.nv.confs[0]!=1))
      {
        SetTextColor(hdc,RGB(255,255,255));
      }
      if((id==IDV_TXTELEMENTCODE) &&
         (wdraw.hwnd==NULL || gv.ev.code!=1))
      {
        SetTextColor(hdc,RGB(255,255,255));
      }
      if((id==IDV_TXTELEMENTAXIS) &&
         (wdraw.hwnd==NULL || gv.ev.axis!=1))
      {
        SetTextColor(hdc,RGB(255,255,255));
      }
      if((id==IDV_TXTHINGE) &&
         (wdraw.hwnd==NULL || gv.ev.hinge!=1))
      {
        SetTextColor(hdc,RGB(255,255,255));
      }
      if((id==IDV_TXTSECTIONCODE) &&
         (wdraw.hwnd==NULL || gv.ev.sectioncode!=1))
      {
        SetTextColor(hdc,RGB(255,255,255));
      }
      if((id==IDV_TXTDEFORMATION) &&
         (wdraw.hwnd==NULL || gv.ev.deformation!=1))
      {
        SetTextColor(hdc,RGB(255,255,255));
      }

      if((id==IDV_TXTDX) &&
         (wdraw.hwnd==NULL || gv.nv.disps[0]!=1))
      {
        SetTextColor(hdc,RGB(255,255,255));
      }
      if((id==IDV_TXTDY) &&
         (wdraw.hwnd==NULL || gv.nv.disps[1]!=1))
      {
        SetTextColor(hdc,RGB(255,255,255));
      }
      if((id==IDV_TXTDZ) &&
         (wdraw.hwnd==NULL || gv.nv.disps[2]!=1))
      {
        SetTextColor(hdc,RGB(255,255,255));
      }

      if((id==IDV_TXTNZ) &&
         (wdraw.hwnd==NULL || gv.ev.stress[0]!=1))
      {
        SetTextColor(hdc,RGB(255,255,255));
      }
      if((id==IDV_TXTQX) &&
         (wdraw.hwnd==NULL || gv.ev.stress[1]!=1))
      {
        SetTextColor(hdc,RGB(255,255,255));
      }
      if((id==IDV_TXTQY) &&
         (wdraw.hwnd==NULL || gv.ev.stress[2]!=1))
      {
        SetTextColor(hdc,RGB(255,255,255));
      }
      if((id==IDV_TXTMZ) &&
         (wdraw.hwnd==NULL || gv.ev.stress[3]!=1))
      {
        SetTextColor(hdc,RGB(255,255,255));
      }
      if((id==IDV_TXTMX) &&
         (wdraw.hwnd==NULL || gv.ev.stress[4]!=1))
      {
        SetTextColor(hdc,RGB(255,255,255));
      }
      if((id==IDV_TXTMY) &&
         (wdraw.hwnd==NULL || gv.ev.stress[5]!=1))
      {
        SetTextColor(hdc,RGB(255,255,255));
      }
      if((id==IDV_TXTREACTION) &&
         (wdraw.hwnd==NULL || gv.nv.react[0]!=1))
      {
        SetTextColor(hdc,RGB(255,255,255));
      }
      return (BOOL)(HBRUSH)GetStockObject(LTGRAY_BRUSH);

    case WM_COMMAND:
      switch(LOWORD(wParam))
      {
        case IDD_OPENRESULT:
          if(wdraw.nchilds>=2 &&
             (wdraw.childs+1)->hwnd!=NULL &&
             (wmenu.childs+2)->vparam.vflag.mv.ftype==F_ARCLM)
          {
            if(arc.nodes==NULL) break;

            fout=fgetstofopen(dir,"r",ID_OUTPUTFILE);  /*OPEN FILE.*/
            if(fout==NULL) break;

            frameoutputtomemory(fout,&arc);
            fclose(fout);

            clearwindow(*(wdraw.childs+1));
            drawarclmframe((wdraw.childs+1)->hdcC,
                           (wdraw.childs+1)->vparam,arc,0,ONSCREEN);
            SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);

            (wdraw.childs+1)->lstatus=ROTATE;
          }
          break;

        case IDV_GLOBALAXIS:
          if(wdraw.hwnd!=NULL)
          {
            flagswitch(&((wdraw.childs+1)->vparam.vflag.axis));
            SendMessage(hdwnd,WM_INITDIALOG,0,0);
          }
          break;
        case IDV_NODECODE:
          if(wdraw.hwnd!=NULL)
          {
            flagswitch(&((wdraw.childs+1)->vparam.vflag.nv.code));
            SendMessage(hdwnd,WM_INITDIALOG,0,0);
          }
          break;
        case IDV_LOADS:
          if(wdraw.hwnd!=NULL)
          {
            flagswitch(&((wdraw.childs+1)
                       ->vparam.vflag.nv.loads[0]));
            flagswitch(&((wdraw.childs+1)
                       ->vparam.vflag.nv.loads[1]));
            flagswitch(&((wdraw.childs+1)
                       ->vparam.vflag.nv.loads[2]));
            flagswitch(&((wdraw.childs+1)
                       ->vparam.vflag.nv.loads[3]));
            flagswitch(&((wdraw.childs+1)
                       ->vparam.vflag.nv.loads[4]));
            flagswitch(&((wdraw.childs+1)
                       ->vparam.vflag.nv.loads[5]));
            SendMessage(hdwnd,WM_INITDIALOG,0,0);
          }
          break;
        case IDV_CONFINEMENT:
          if(wdraw.hwnd!=NULL)
          {
            flagswitch(&((wdraw.childs+1)
                       ->vparam.vflag.nv.confs[0]));
            flagswitch(&((wdraw.childs+1)
                       ->vparam.vflag.nv.confs[1]));
            flagswitch(&((wdraw.childs+1)
                       ->vparam.vflag.nv.confs[2]));
            flagswitch(&((wdraw.childs+1)
                       ->vparam.vflag.nv.confs[3]));
            flagswitch(&((wdraw.childs+1)
                       ->vparam.vflag.nv.confs[4]));
            flagswitch(&((wdraw.childs+1)
                       ->vparam.vflag.nv.confs[5]));
            SendMessage(hdwnd,WM_INITDIALOG,0,0);
          }
          break;
        case IDV_ELEMENTCODE:
          if(wdraw.hwnd!=NULL)
          {
            flagswitch(&((wdraw.childs+1)->vparam.vflag.ev.code));
            SendMessage(hdwnd,WM_INITDIALOG,0,0);
          }
          break;
        case IDV_ELEMENTAXIS:
          if(wdraw.hwnd!=NULL)
          {
            flagswitch(&((wdraw.childs+1)->vparam.vflag.ev.axis));
            SendMessage(hdwnd,WM_INITDIALOG,0,0);
          }
          break;
        case IDV_HINGE:
          if(wdraw.hwnd!=NULL)
          {
            flagswitch(&((wdraw.childs+1)->vparam.vflag.ev.hinge));
            SendMessage(hdwnd,WM_INITDIALOG,0,0);
          }
          break;
        case IDV_SECTIONCODE:
          if(wdraw.hwnd!=NULL)
          {
            flagswitch(&((wdraw.childs+1)
                       ->vparam.vflag.ev.sectioncode));
            SendMessage(hdwnd,WM_INITDIALOG,0,0);
          }
          break;

        case IDV_MODENUM:
          if(wdraw.hwnd!=NULL)
          {
            (wdraw.childs+1)->vparam.vflag.ev.deformation=0;
            SendMessage(hdwnd,WM_INITDIALOG,0,0);
          }
          break;

        case IDV_DEFORMATION:
          if(wdraw.hwnd!=NULL &&
             !((wmenu.childs+2)->vparam.vflag.mv.ftype==F_ARCLM &&
               arc.nlaps==0))
          {
            GetDlgItemText(hdwnd,IDV_MODENUM,str,80);
            mode=strtol(str,NULL,10);

            if(mode>arc.nlaps) break;

            /*copyform(&arc,*((arc.eigenvec)+mode-1));*/

            flagswitch(&((wdraw.childs+1)
                       ->vparam.vflag.ev.deformation));
            SendMessage(hdwnd,WM_INITDIALOG,0,0);
          }
          break;
        case IDV_DX:
          if(wdraw.hwnd!=NULL)
          {
            flagswitch(&((wdraw.childs+1)
                       ->vparam.vflag.nv.disps[0]));
            SendMessage(hdwnd,WM_INITDIALOG,0,0);
          }
          break;
        case IDV_DY:
          if(wdraw.hwnd!=NULL)
          {
            flagswitch(&((wdraw.childs+1)
                       ->vparam.vflag.nv.disps[1]));
            SendMessage(hdwnd,WM_INITDIALOG,0,0);
          }
          break;
        case IDV_DZ:
          if(wdraw.hwnd!=NULL)
          {
            flagswitch(&((wdraw.childs+1)
                       ->vparam.vflag.nv.disps[2]));
            SendMessage(hdwnd,WM_INITDIALOG,0,0);
          }
          break;
        case IDV_NZ:
          if(wdraw.hwnd!=NULL)
          {
            flagswitch(&((wdraw.childs+1)
                       ->vparam.vflag.ev.stress[0]));
            SendMessage(hdwnd,WM_INITDIALOG,0,0);
          }
          break;
        case IDV_QX:
          if(wdraw.hwnd!=NULL)
          {
            flagswitch(&((wdraw.childs+1)
                       ->vparam.vflag.ev.stress[1]));
            SendMessage(hdwnd,WM_INITDIALOG,0,0);
          }
          break;
        case IDV_QY:
          if(wdraw.hwnd!=NULL)
          {
            flagswitch(&((wdraw.childs+1)
                       ->vparam.vflag.ev.stress[2]));
            SendMessage(hdwnd,WM_INITDIALOG,0,0);
          }
          break;
        case IDV_MZ:
          if(wdraw.hwnd!=NULL)
          {
            flagswitch(&((wdraw.childs+1)
                       ->vparam.vflag.ev.stress[3]));
            SendMessage(hdwnd,WM_INITDIALOG,0,0);
          }
          break;
        case IDV_MX:
          if(wdraw.hwnd!=NULL)
          {
            flagswitch(&((wdraw.childs+1)
                       ->vparam.vflag.ev.stress[4]));
            SendMessage(hdwnd,WM_INITDIALOG,0,0);
          }
          break;
        case IDV_MY:
          if(wdraw.hwnd!=NULL)
          {
            flagswitch(&((wdraw.childs+1)
                       ->vparam.vflag.ev.stress[5]));
            SendMessage(hdwnd,WM_INITDIALOG,0,0);
          }
          break;
        case IDV_REACTION:
          if(wdraw.hwnd!=NULL)
          {
            flagswitch(&((wdraw.childs+1)
                       ->vparam.vflag.nv.react[0]));
            flagswitch(&((wdraw.childs+1)
                       ->vparam.vflag.nv.react[1]));
            flagswitch(&((wdraw.childs+1)
                       ->vparam.vflag.nv.react[2]));
            flagswitch(&((wdraw.childs+1)
                       ->vparam.vflag.nv.react[3]));
            flagswitch(&((wdraw.childs+1)
                       ->vparam.vflag.nv.react[4]));
            flagswitch(&((wdraw.childs+1)
                       ->vparam.vflag.nv.react[5]));
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

            hcursor = LoadCursor(hInstGlobal,"CANBOX");    /*CURSOR*/
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

            hcursor = LoadCursor(hInstGlobal,"CANBOX");    /*CURSOR*/
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

            hcursor = LoadCursor(hInstGlobal,"CANBOX");    /*CURSOR*/
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
            pd.hInstance = (HANDLE)NULL;
            pd.lCustData = 0L;
            pd.lpfnPrintHook = (LPPRINTHOOKPROC)NULL;
            pd.lpfnSetupHook = (LPSETUPHOOKPROC)NULL;
            pd.lpPrintTemplateName = (LPSTR)NULL;
            pd.lpSetupTemplateName = (LPSTR)NULL;
            pd.hPrintTemplate = (HANDLE)NULL;
            pd.hSetupTemplate = (HANDLE)NULL;

            if(PrintDlg(&pd) != FALSE)
            {
              caps = GetDeviceCaps(pd.hDC,POLYGONALCAPS);
              sprintf(str,"Caps:%#x=%d",caps,caps);
              /*MessageBox(hdwnd,str,"Print",MB_OK);*/

              cWidthPels = GetDeviceCaps(pd.hDC,HORZRES);
              cHeightPels = GetDeviceCaps(pd.hDC,VERTRES);

              vprint=(wdraw.childs+1)->vparam;
              vprint.Xo=0.5*cWidthPels;
              vprint.Yo=0.5*cHeightPels;

              di.cbSize = sizeof(DOCINFO);
              sprintf(doc,"CanvsPrint");
              di.lpszDocName = doc;
              di.lpszOutput = NULL;

              setfontformat(pd.hDC,prn.jiheight,prn.jiwidth,
                            "Terminal",0,0,255);
              SetBkMode(pd.hDC,TRANSPARENT);

              StartDoc(pd.hDC,&di);
              StartPage(pd.hDC);

              if((wmenu.childs+2)->vparam.vflag.mv.ftype==F_ARCLM &&
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
              drawtexts(pd.hDC,wdraw.strset,wdraw.nstring,vprint);

              EndPage(pd.hDC);
              EndDoc(pd.hDC);

              hfont=GetCurrentObject(pd.hDC,OBJ_FONT);
              DeleteObject(hfont);
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

BOOL CALLBACK DialogProcMenu3(HWND hdwnd,
                              UINT message,
                              WPARAM wParam,
                              LPARAM lParam)
/*OPTION MENU DIALOG BOX 3:FRAME CREATION.*/
{
  HWND hitem;
  HDC hdc;
  int id;

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
      break;

    case WM_CTLCOLORSTATIC: /*STATIC TEXT MENU.*/
      DefWindowProc(hdwnd,message,wParam,lParam);

      hdc=(HDC)wParam;
      hitem=(HWND)lParam;
      id=GetDlgCtrlID(hitem);

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
      return (BOOL)(HBRUSH)GetStockObject(LTGRAY_BRUSH);

    case WM_COMMAND:
      switch(LOWORD(wParam))
      {
        case IDC_CREATE:
          if(createcommand!=C_NEUTRAL)
          {
            createcommand=C_NEUTRAL;
            createitem=C_NEUTRAL;
            SendMessage(hdwnd,WM_INITDIALOG,0,0);
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
          break;

        case IDC_ADD:
          createcommand=C_ADD;
          SendMessage(hdwnd,WM_INITDIALOG,0,0);
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
      }
      break;

    default:
      return DefWindowProc(hdwnd,message,wParam,lParam);
  }
  return 0;
}/*DialogProcMenu3*/

BOOL CALLBACK DialogProcText(HWND hdwnd,
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

BOOL CALLBACK DialogProcSsim(HWND hdwnd,
                             UINT message,
                             WPARAM wParam,
                             LPARAM lParam)
/*POPUP DIALOG BOX FOR SIMPLE SECTION LIST.*/
{
  HDC hdc;
  HWND hitem;
  RECT rect;
  HBRUSH hbrush,pbrush;
  char str[20];
  int id;
  long int code;

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

    case WM_PAINT:
      /*DefWindowProc(hdwnd,message,wParam,lParam);*/

      GetDlgItemText(hdwnd,IDDS_TXTSECTCODE,str,20);
      code=strtol(str,NULL,10);

      gsect=currentsects;
      while(code != gsect->code) gsect++;

      rect.top=0;
      rect.bottom=19;
      rect.left=0;
      rect.right=14;
      hbrush = (HBRUSH)CreateSolidBrush(RGB(gsect->dcolor.r,
                                            gsect->dcolor.g,
                                            gsect->dcolor.b));
      hitem=GetDlgItem(hdwnd,IDDS_COLOR);
      hdc=GetDC(hitem);

      FillRect(hdc,&rect,hbrush);

      DeleteObject(hbrush);
      ReleaseDC(hdwnd,hdc);
      /*DefWindowProc(hdwnd,message,wParam,lParam);*/
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

BOOL CALLBACK EnumChildProcSheet(HWND hwnd,LPARAM lParam)
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
  wcl.hCursor = LoadCursor(hInstance,"CANCURSOR");         /*CURSOR*/
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
    hmemorybit = SelectObject(hdccompati,hbit);

    hbrush = (HBRUSH)GetClassLong(w.hwnd,GCL_HBRBACKGROUND);

    hmemorybrush = SelectObject(hdccompati,hbrush);
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
  ppen = SelectObject(hdc,hpen1);
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
  pbit = SelectObject(wp->hdcC,hbit);

  PatBlt(wp->hdcC,maxX-BARWIDTH-1,1,
                  BARWIDTH,maxY-BARWIDTH-5,PATCOPY); /*DELETE BAR.*/

  hbrush = CreateSolidBrush(RGB(255,255,0)); /*REDRAW BAR*/
  pbrush = SelectObject(wp->hdcC,hbrush);

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
  pbit = SelectObject(wp->hdcC,hbit);

  PatBlt(wp->hdcC,1,maxY-BARWIDTH-1,
                  maxX-BARWIDTH-5,BARWIDTH,PATCOPY); /*DELETE BAR.*/

  hbrush = CreateSolidBrush(RGB(255,255,0)); /*REDRAW BAR*/
  pbrush = SelectObject(wp->hdcC,hbrush);

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

struct snode *addtext(struct snode *strset,int *nstr,char *str,
                      double tx,double ty)
{
  (*nstr)++;
  strset=(struct snode *)realloc(strset,
                                 (*nstr)*sizeof(struct snode));
  strcpy((strset+(*nstr)-1)->str,str);
  (strset+(*nstr)-1)->n.d[0]=tx;
  (strset+(*nstr)-1)->n.d[1]=ty;
  return strset;
}/*addtext*/

