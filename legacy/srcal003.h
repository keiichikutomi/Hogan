/*SRCAL003.H:SRC SUBROUTINES SINCE 1995.12.21.JUNSATO.*/
/*BASED ON SRCAM003.H.*/
/*HEADER FOR CANVS000.C.*/
/*LAST CHANGE:1997.1.25.*/

/*S,RC REMAINING DIFFERENT PLAIN WITH SAME NEUTRAL AXIS.*/

/*OUTPUT ALL COMMENTS.*/

/*S COLUMN...................................................*/
/*RC COLUMN..................................................*/
/*RC WALL....................................................*/
/*SRC COLUMN.................................................*/

/*SECTION FROM SECTION LIST FILE.*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <malloc.h>
#include <math.h>
#include <limits.h>

#ifndef PI
#define PI 3.1415926535897932384
#endif

#define S   1                                        /*SECTION TYPE*/
#define RC  2
#define SRC 3

#define COLUMN 1                                     /*ELEMENT TYPE*/
#define GIRDER 2
#define BEAM   3
#define BRACE  4
#define WALL   5
#define SLAB   6

#define SX 0 /*AXIS OF SECTION*/
#define SY 1

#define HSTRONG 0 /*AXIS OF H*/
#define HWEAK 1

#define HEAD 0 /*END OF ELEMENT*/
#define TAIL 1

#define MAXSRECT    5                     /*MAX OF STEEL RECTANGLES*/
#define MAXREINS   20                       /*MAX OF REINFORCEMENTS*/
#define MAXCRECT    3                  /*MAX OF CONCRETE RECTANGLES*/

struct materialrect{double top,bottom,left,right;}; /*RECTANGLE[mm]*/
struct reinforcement{double area,x,y;};    /*REIN AREA,POSITION[mm]*/

struct section{
               int code;                   /*CODE NUMBER OF SECTION*/
               int stype,etype;         /*SECTION TYPE,ELEMENT TYPE*/
               int nsteel,nrein,nconc;     /*STEELS,REINS,CONCRETES*/
               struct materialrect srect[MAXSRECT];    /*STEEL RECT*/
               struct reinforcement rein[MAXREINS];         /*REINS*/
               struct materialrect crect[MAXCRECT]; /*CONCRETE RECT*/
               double shearrein[2];        /*RATE OF REIN FOR SHEAR*/
               double thick,wlength,wheight,windowrate;
               double face[2][2];          /*FACE LENGTH[AXIS][END]*/
               double safetya;       /*MAX RATE OF SAFETY ALLOWABLE*/
               double safetyu;        /*MAX RATE OF SAFETY ULTIMATE*/
              };


extern HDC hdcdraw,drawback,drawcompati;


void drawmaterials(HWND hdraw,HWND hdoc,HWND herr,
                   HWND hdwnd,struct visible v);
void drawmaterialrect(HWND hdraw,struct materialrect r,
                      int red,int green,int blue);
void drawreinforcement(HWND hdraw,struct reinforcement r,
                       int red,int green,int blue);

void freestr(char **str,int nstr);

int getsectionform(FILE *flist,int code,struct section *sect);
int getcodelist(FILE *flist,int codelist[]);

int quadraticequation(double c2,double c1,double c0,
                      double answer[]);
int cubicequation(double c3,double c2,double c1,double c0,
                  double answer[]);

double coeffA(int nrect,struct materialrect rect[]);
double coeffI(int nrect,struct materialrect rect[],int axis,
              double *yg,double *yt,double *yc);
double coeffZ(int nrect,struct materialrect rect[],int axis);
double coeffi(int nrect,struct materialrect rect[],int axis);


void drawmaterials(HWND hdraw,HWND hdoc,HWND herr,
                   HWND hdwnd,struct visible v)
/*DRAW SECTION MATERIALS.*/
{
  FILE *flist;
  char str[80];
  int i,code;
  struct section sect;

  RECT rect;
  long int maxX,maxY;

  if(hdraw==NULL) return;

  if(drawcompati==NULL)
  {
    hdcdraw=GetDC(hdraw);

    drawcompati=createcompati(hdraw,GCL_HBRBACKGROUND);
    drawback=createcompati(hdraw,BLACK_BRUSH);
    setfontformat(drawcompati,0,0,255);
  }
  else
  {
    GetClientRect(hdraw,&rect);
    maxX = rect.right - rect.left;
    maxY = rect.bottom - rect.top;

    PatBlt(drawcompati,0,0,maxX,maxY,PATCOPY);
  }

  if(v.globalaxis)
  {
    drawcoord(hdraw,0,0,255);          /*DRAW COORDINATION ARROW.*/
  }

  GetDlgItemText(hdwnd,ID_INPUTFILE,str,80);
  flist = fopen(str,"r");
  if(flist == NULL) return;

  GetDlgItemText(hdwnd,IDS_CODE,str,80);
  code=strtol(str,NULL,10);

  if(getsectionform(flist,code,&sect))
  {
    for(i=0;i<(sect.nsteel);i++)
    {
      drawmaterialrect(hdraw,sect.srect[i],255,0,0);
    }
    for(i=0;i<(sect.nrein);i++)
    {
      drawreinforcement(hdraw,sect.rein[i],255,255,0);
    }
    for(i=0;i<(sect.nconc);i++)
    {
      drawmaterialrect(hdraw,sect.crect[i],0,0,255);
    }
  }

  fclose(flist);

  GetClientRect(hdraw,&rect);         /*UPDATE WINDOW "DRAWINGS".*/
  maxX = rect.right - rect.left;
  maxY = rect.bottom - rect.top;
  PatBlt(drawback,0,0,maxX,maxY,PATCOPY);

  updatedrawings(hdraw,hdoc,herr); /*DRAW FRAME WITH DEFORMATION.*/

  return;
}/*drawmaterials*/


void drawmaterialrect(HWND hdraw,struct materialrect r,
                      int red,int green,int blue)
/*DRAW MATERIAL RECT.*/
{
  struct structline line;

  line.r=red; line.g=green; line.b=blue; /*LINE COLOR.*/

  line.x1=r.left;  line.y1=r.bottom; line.z1=0.0;
  line.x2=r.right; line.y2=r.bottom; line.z2=0.0;
  hogline(hdraw,line);

  line.x1=r.right; line.y1=r.bottom; line.z1=0.0;
  line.x2=r.right; line.y2=r.top;    line.z2=0.0;
  hogline(hdraw,line);

  line.x1=r.right; line.y1=r.top;    line.z1=0.0;
  line.x2=r.left;  line.y2=r.top;    line.z2=0.0;
  hogline(hdraw,line);

  line.x1=r.left;  line.y1=r.top;    line.z1=0.0;
  line.x2=r.left;  line.y2=r.bottom;    line.z2=0.0;
  hogline(hdraw,line);

  return;
}/*drawmaterialrect*/


void drawreinforcement(HWND hdraw,struct reinforcement rein,
                       int red,int green,int blue)
/*DRAW REINFORCEMENT ON DRAWINGS.*/
{
  HDC hdc;
  HPEN hpen,hmemorypen;                             /*HANDLE OF PEN*/
  RECT rect;                                  /*RECTANGLE OF HDRAW.*/
  double radius;
  struct vector center;

  struct viewparam vparam;

  GetClientRect(hdraw,&rect);
  vparam.Xo=(int)((rect.right - rect.left)/2);           /*ORIGIN X*/
  vparam.Yo=(int)((rect.bottom - rect.top)/2);           /*ORIGIN Y*/
  vparam.type=1;           /*VIEW TYPE 0:AXONOMETRICS 1:PERSPECTIVE*/
  vparam.gfactor=0.04;

  createviewparam(&vparam);

  /*hdc=GetDC(hdraw);*/
  hdc = drawcompati;

  hpen=CreatePen(PS_SOLID,1,RGB(red,green,blue));      /*PEN COLOR.*/
  hmemorypen=SelectObject(hdc,hpen);

  center.dc[GX]=rein.x;
  center.dc[GY]=rein.y;
  center.dc[GZ]=0.0;

  radius=sqrt(rein.area/PI);

  circledrawings(hdc,vparam,center,radius);

  SelectObject(hdc,hmemorypen);
  DeleteObject(hpen);

  return;
}/*drawreinforcement*/


void freestr(char **str,int nstr)
/*FREE STRINGS ALLOCATED BY "fgetsbrk".*/
{
  for(;nstr>0;nstr--) free(*(str+nstr-1));
  free(str);

  return;
}/*freestr*/


int getsectionform(FILE *flist,int code,struct section *sect)
/*GET SECTION FROM SECTION LIST.*/
{
  char **data;
  int n,ns=0;
  int nsteel=0,nrein=0,nconc=0;

  fseek(flist,0L,SEEK_SET);

  sect->code=0;                                   /*INITIALIZATION.*/
  sect->stype=0;
  sect->etype=0;

  sect->nsteel=0;
  sect->nrein=0;
  sect->nconc=0;

  sect->face[0][0]=0.0;
  sect->face[0][1]=0.0;
  sect->face[1][0]=0.0;
  sect->face[1][1]=0.0;

  while(1)
  {
    data=fgetsbrk(flist,&n);
    if(n==0) return 0;

    if(!strcmp(*(data+0),"CODE"))
    {
      ns++;                                         /*OFFSET COUNT.*/
      sect->code=(int)strtol(*(data+1),NULL,10);

      if(sect->code==code)
      {
        if(!strcmp(*(data+2),"S"))   sect->stype=S;
        if(!strcmp(*(data+2),"RC"))  sect->stype=RC;
        if(!strcmp(*(data+2),"SRC")) sect->stype=SRC;

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

        sect->nsteel=nsteel;
        sect->nrein=nrein;
        sect->nconc=nconc;

        return ns;
      }
    }
    freestr(data,n);
  }
}/*getsection*/


int getcodelist(FILE *flist,int codelist[])
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
      codelist[ns]=(int)strtol(*(data+1),NULL,10);
      ns++;
    }
    freestr(data,n);
  }
}/*getcodelist*/


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

