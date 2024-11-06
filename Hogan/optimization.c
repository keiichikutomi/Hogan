/*OPTIMIZATION*/

void conjugategradientbezier(struct organ *org); 	       	  /*CONJUGATE*/
void conjugategradient(struct organ *org); 	       	  /*CONJUGATE*/
void conjugategradientxz(struct organ *org);          /*Fujimoto*/
void conjugategradientcurve(struct organ *org);       /*Fujimoto*/
double arclmautomatic(struct organ *org);  			  /*CONJUGATE*/
double arclmautomaticII(struct organ *org);  		  /*CONJUGATE*/
double arclmautomatic101(struct organ *org);  		  /*CONJUGATE*/

/*-----------------------------------------------------------------*/
/*Bclng Optimization (CG Method)*/
void bclngconjugategradient(struct organ *org);       /*Ujioka*/
double bclngautomatic(struct organ *org);             /*Ujioka*/
double beigen(struct arclmframe *af);                 /*Ujioka*/
/*-----------------------------------------------------------------*/

void createziersurfacetest(struct organ *org);
int beziersurface(int m,int n,double *x,double *y,double *z,
                  double u,double v,struct onode *node);
int beziersurfaceII(int m,int n,double *x,double *y,double *z,
                    double u,double v,struct onode *node);
double bernstein(int n,int i,double t);
int fact(int n);
void drawcontrolepoint(HDC hdc,struct viewparam vp,struct onode gn);

/*EXTERNAL PARAMETERS*/
extern void clearwindow(struct windowparams wp);
extern int globaldrawflag;      /*DRAW FLAG*/
extern struct arclmframe arc;   /*GLOBAL ARCLM FRAME.*/
extern struct arclmframe arci;  /*NULL ARCLM FRAME.*/
extern struct arclmframe arcx;  /*ARCLM FRAME FOR X LOAD.*/
extern struct arclmframe arcy;  /*ARCLM FRAME FOR Y LOAD.*/

/*-----------------------------------------------------------------*/
void conjugategradientbezier(struct organ *org)     /*CONJUGATE*/
/*OPTIMIZE ORGAN WITH CONJUGATE GRADIENT.*/
{
  FILE /**fout,*/*ftxt;
  FILE *fresult;
  char non[10],str[256],s[256];
  int i,ii,j,jj,k,m,n;
  int nnode,nelem,aelem;
  double df,f1,f2,fa,ftarget,c1,c2,c3,alpha,beta,gamma,vsize,eps;
  double *x1,*x2,*xa,*xx,*dx,dfact; /*COORDINATES*/
  double *u1,*u2,*ua,*fgrad1,*fgrad2; /*GRADIENT VECTOR*/
  double **cmtx; /*MATRIX*/
  double xi,tau,test; /*FOR LINE SEARCH(ARMIJO RULE)*/

  /*BEZIER SURFACE*/
  FILE *fbezier;
  char **data/*,str[256]="\0"*/;
  double ddata;
  int ndata;
  long int /*nnode,*/ncode;
  double *x,*y,*z;     /*CONTROLE POINTS*/
  double *xini,*yini,*zini;     /*INITIAL CONTROLE POINTS*/
  int n1,n2;           /*DEGREE*/
  int ncontrole;
  double *u,*v;      /*U-V COORDINATES*/
  struct onode node;

  /*FIXED CONTROLE POINTS*/
  int nfix;
  int *fp;    /*INDEX OF FIXED CONTROLE POINTS*/
  int fixed;
  int offset;

  ftxt=fopen("gradtest.txt","w");

  ftarget=5.0;
  gamma=0.01;
  dfact=0.1;
  eps=0.001;

  /*POLYGON01.INP : SUCCESSFUL PARAMETERS*/
  /*
  ftarget=0.4;
  gamma=0.001;
  dfact=0.05;
  eps=0.00001;
  */
  /*
  ftarget=0.9;
  gamma=0.002;
  dfact=0.05;
  eps=0.00001;
  */
  /*
  ftarget=0.9;
  gamma=0.004;
  dfact=0.1;
  eps=0.00001;
  */
  /*
  ftarget=0.9;
  gamma=0.03;
  dfact=0.1;
  eps=0.001;
  */

  /*PARAMETERS FOR LINE SEARCH(ARMIJO RULE)*/
  xi=0.001; /*0.0<=xi<=1.0*/
  tau=0.9; /*0<tau<1*/

  /*FIX CONTROLE POINTS*/
  nfix=3;    /*NUMBER OF FIXED CONTROLE POINTS*/
  fp=(int *)malloc(nfix*sizeof(int));
  fp[0]=3;    /*INDEX OF FIXED CONTROLE POINTS*/
  fp[1]=12;
  fp[2]=15;
//  fp[0]=0;
//  fp[1]=15;

  nnode=org->nnode;
  nelem=org->nelem;

  /*CREATE INITIAL BEZIER SURFACE*/
  fbezier=fopen("bezier.txt","r");   /*BEZIER SURFACE DATA*/
  if(fbezier==NULL)
  {
    errormessage("ACCESS IMPOSSIBLE.");
    return;
  }
  fseek(fbezier,0L,SEEK_SET);

  data=fgetsbrk(fbezier,&ndata);
  n1=strtol(*(data+0),NULL,10);
  n2=strtol(*(data+1),NULL,10);
  if(nnode!=strtol(*(data+2),NULL,10)) return;

  /*CONTROLE POINTS*/
  ncontrole=(n1+1)*(n2+1);
  x=(double *)malloc((ncontrole)*sizeof(double));
  y=(double *)malloc((ncontrole)*sizeof(double));
  z=(double *)malloc((ncontrole)*sizeof(double));
  zini=(double *)malloc((ncontrole)*sizeof(double));

  for(i=0;i<ncontrole;i++)
  {
	data=fgetsbrk(fbezier,&ndata);
	if(ndata!=3) return;
    x[i]=strtod(*(data+0),NULL);
    y[i]=strtod(*(data+1),NULL);
    z[i]=strtod(*(data+2),NULL);

    for(;ndata>0;ndata--) free(*(data+ndata-1));
    free(data);

    zini[i]=z[i];
  }

  /*U-V COORDINATE OF EACH NODE*/
  u=(double *)malloc((nnode)*sizeof(double));
  v=(double *)malloc((nnode)*sizeof(double));
  for(i=0;i<nnode;i++)
  {
	data=fgetsbrk(fbezier,&ndata);
	if(ndata!=5) return;

	ncode=strtol(*(data+1),NULL,10);
	if(ncode!=(org->nodes+i)->code) return;

    ddata=strtod(*(data+3),NULL);
    *(u+i)=ddata;

    ddata=strtod(*(data+4),NULL);
    *(v+i)=ddata;

    for(;ndata>0;ndata--) free(*(data+ndata-1));
    free(data);
  }
  fclose(fbezier);

  /*BEZIER SURFACE*/
  for(i=0;i<nnode;i++)
  {
    node=*(org->nodes+i);
    beziersurfaceII(n1,n2,x,y,z,*(u+i),*(v+i),&node);

    (org->nodes+i)->d[0]=node.d[0];
    (org->nodes+i)->d[1]=node.d[1];
    (org->nodes+i)->d[2]=node.d[2];
  }

  m=ncontrole-nfix;             /*NUMBER OF PARAMETERS*/
  fgrad1=mallocdoublevector(m); /*FOR Z OF CONTROLE POINTS*/
  fgrad2=mallocdoublevector(m); /*FOR Z OF CONTROLE POINTS*/
  u1=mallocdoublevector(m);
  u2=mallocdoublevector(m);
  ua=mallocdoublevector(m);

  /*INITIAL NODE*/
  /*for(i=0; i<nnode; i++) *(ninit+i)=*(org->nodes+i);*/

  f1=arclmautomatic101(org);

  sprintf(str,"%d",0);
  SetDlgItemText((wmenu.childs+2)->hwnd,ID_LAPS,str);
  SendDlgItemMessage((wmenu.childs+2)->hwnd,ID_LAPS,WM_PAINT,0,0);
  sprintf(str,"%.5f",f1);
  SetDlgItemText((wmenu.childs+2)->hwnd,ID_SAFETY,str);
  SendDlgItemMessage((wmenu.childs+2)->hwnd,ID_SAFETY,WM_PAINT,0,0);

  sprintf(str,"Initial Strain Energy = %9.5f",f1);
  fprintf(ftxt,"%s\n",str);
  MessageBox(NULL,str,"Conjugate Gradient",MB_OK);

  /*INITIAL GRADIENT*/
  fprintf(ftxt,"Initial Gradient\n");

  offset=0; i=0;
  while(i<m)
  {
    fixed=0;
    for(ii=0;ii<nfix;ii++)
    {
      if(offset==fp[ii])
      {
        fixed=1;
      }
    }
    if(!fixed)
    {
  	  for(j=0;j<ncontrole;j++)
	  {
	    if(j==offset)
	    {
          z[j]=zini[j]+dfact;                 /*UPDATE ONLY Z OF CONTROLE POINT I*/
	    }
	    else
	    {
          z[j]=zini[j];                       /*RESET Z OF OTHER CONTROLE POINTS*/
	    }
	  }
      /*BEZIER SURFACE*/
      for(ii=0;ii<nnode;ii++)
      {
        node=*(org->nodes+ii);
        beziersurfaceII(n1,n2,x,y,z,*(u+ii),*(v+ii),&node);

        (org->nodes+ii)->d[0]=node.d[0];
        (org->nodes+ii)->d[1]=node.d[1];
        (org->nodes+ii)->d[2]=node.d[2];
      }

	  /*Strain Energy*/
      df=arclmautomatic101(org);

	  *(fgrad1+i)=gamma*(f1-df)/dfact; /*NEGATIVE GRADIENT*/
	  *(u1+i)=*(fgrad1+i);

  	  fprintf(ftxt,"CONTROLE POINT[%d] Gradient=%9.5f\n",offset,(f1-df)/dfact);

      i++;
    }

    offset++;
  }
  fprintf(ftxt,"\n");

  k=0;
  while(1)
  {
	k++;
	sprintf(str,"%d",k);
	SetDlgItemText((wmenu.childs+2)->hwnd,ID_LAPS,str);

#if 1
    offset=0; i=0;
    while(i<m)
    {
      fixed=0;
      for(ii=0;ii<nfix;ii++)
      {
        if(offset==fp[ii])
        {
          fixed=1;
        }
      }
      if(!fixed)
      {
        z[offset]=zini[offset]+(*(u1+i)); /*ONLY Z*/
        i++;
      }

      offset++;
	}

    for(ii=0;ii<nnode;ii++)
    {
      node=*(org->nodes+ii);
      beziersurfaceII(n1,n2,x,y,z,*(u+ii),*(v+ii),&node);

      (org->nodes+ii)->d[0]=node.d[0];
      (org->nodes+ii)->d[1]=node.d[1];
      (org->nodes+ii)->d[2]=node.d[2];
    }

	/*Strain Energy*/
	fa=arclmautomatic101(org);

//	fprintf(ftxt,"Step %d Max Safety = %9.5f\n",k,fa);
	fprintf(ftxt,"Step %d Strain Energy = %9.5f\n",k,fa);
#endif
	/*
	sprintf(str,"Max Safety = %.5f",fa);
	MessageBox(NULL,str,"Conjugate Gradient",MB_OK);
	*/

	/*INITIAL NODE*/
    for(i=0;i<ncontrole;i++) z[i]=zini[i];

	/*GRADIENT*/
	fprintf(ftxt,"Step %d Gradient\n",k);

    offset=0; i=0;
    while(i<m)
    {
      fixed=0;
      for(ii=0;ii<nfix;ii++)
      {
        if(offset==fp[ii])
        {
          fixed=1;
        }
      }
      if(!fixed)
      {
  	    for(j=0;j<ncontrole;j++)
	    {
	      if(j==offset)
	      {
            z[j]=zini[j]+dfact;                 /*UPDATE ONLY Z OF CONTROLE POINT I*/
	      }
	      else
	      {
            z[j]=zini[j];                       /*RESET Z OF OTHER CONTROLE POINTS*/
	      }
        }
        /*BEZIER SURFACE*/
        for(ii=0;ii<nnode;ii++)
        {
          node=*(org->nodes+ii);
          beziersurfaceII(n1,n2,x,y,z,*(u+ii),*(v+ii),&node);

          (org->nodes+ii)->d[0]=node.d[0];
          (org->nodes+ii)->d[1]=node.d[1];
          (org->nodes+ii)->d[2]=node.d[2];
        }

	    /*Strain Energy*/
	    df=arclmautomatic101(org);

//	    *(ua+i)=-gamma*(fa-df)/dfact+(*(u1+i)); /*WRONG:NEGATIVE GRADIENT*/
	    *(ua+i)=-gamma*(fa-df)/dfact+(*(fgrad1+i)); /*NEGATIVE GRADIENT*/
        *(fgrad2+i)=gamma*(f1-df)/dfact;  /*NEGATIVE GRADIENT*/
        fprintf(ftxt,"CONTROLE POINT[%d] Gradient=%9.5f\n",offset,(f1-df)/dfact);

        i++;
      }

      offset++;
    }
	fprintf(ftxt,"\n");

	/*
	fprintf(stderr,"ITERATION %d\n",k);
	fprintf(stderr," {u}= %8.3f %8.3f\n",*(u1+0),*(u1+1));
	fprintf(stderr,"A{u}= %8.3f %8.3f\n",*(ua+0),*(ua+1));
	*/

	c1=0.0;
	c2=0.0;
	for(i=0; i<m; i++)
	{
	  c1+=(*(fgrad1+i))*(*(fgrad1+i));
	  c2+=(*(u1+i))*(*(ua+i));
	}
	alpha=c1/c2;
//	alpha=1.0;

	/*
	fprintf(stderr,"Alpha= %8.3f = %8.3f / %8.3f\n",alpha,c1,c2);
	*/
	sprintf(s,"BACKTRACKING LINE SEARCH BEGIN:INITIAL ALPHA= %8.3f/f1=%.5f",alpha,f1);
//	MessageBox(NULL,s,"Conjugate Gradient",MB_OK);
    fprintf(ftxt,"%s\n",s);

    while(1)   /*ARMIJO RULE*/
    {
	  vsize=0.0;
      offset=0; i=0;
      while(i<m)
      {
        fixed=0;
        for(ii=0;ii<nfix;ii++)
        {
          if(offset==fp[ii])
          {
            fixed=1;
          }
        }
        if(!fixed)
        {
          z[offset]=zini[offset]+alpha*(*(u1+i)); /*ONLY Z*/
//  	      *(fgrad2+i)=(*(fgrad1+i))-alpha*(*(ua+i));
	      vsize+=(*(fgrad2+i))*(*(fgrad2+i));

          i++;
        }

        offset++;
	  }

      /*BEZIER SURFACE*/
      for(ii=0;ii<nnode;ii++)
      {
        node=*(org->nodes+ii);
        beziersurfaceII(n1,n2,x,y,z,*(u+ii),*(v+ii),&node);

        (org->nodes+ii)->d[0]=node.d[0];
        (org->nodes+ii)->d[1]=node.d[1];
        (org->nodes+ii)->d[2]=node.d[2];
      }

	  /*Strain Energy*/
	  f2=arclmautomatic101(org);

      /*ARMIJO RULE*/
      test=0.0;
      for(i=0;i<m;i++)  test+=(*(fgrad2+i))*(*(u1+i));
      test=f2-(f1-abs(xi*alpha*test));

      if(test<=0.0 || (alpha<=0.01 && f1>f2))
      {
        sprintf(s,"LINE SEARCH:ALPHA= %8.3f/f1=%.5f/f2=%.5f/TEST= %8.3f",alpha,f1,f2,test);
	    fprintf(ftxt,"%s\n",s);
        break;
      }
      else
      {
        alpha*=tau;
      	sprintf(s,"LINE SEARCH:ALPHA= %8.3f/f1=%.5f/f2=%.5f/TEST= %8.3f",alpha,f1,f2,test);
	    fprintf(ftxt,"%s\n",s);
//        MessageBox(NULL,s,"Conjugate Gradient",MB_OK);
      }
    }
	sprintf(str,"%d",k);
	SetDlgItemText((wmenu.childs+2)->hwnd,ID_LAPS,str);
	SendDlgItemMessage((wmenu.childs+2)->hwnd,ID_LAPS,WM_PAINT,0,0);
	sprintf(str,"%.5f",f2);
	SetDlgItemText((wmenu.childs+2)->hwnd,ID_SAFETY,str);
	SendDlgItemMessage((wmenu.childs+2)->hwnd,ID_SAFETY,WM_PAINT,0,0);

	sprintf(str,"Step=%d Strain Energy=%.5f Vector Size=%.5f",k,f2,vsize);
	fprintf(ftxt,"%s\n\n",str);
	/*MessageBox(NULL,str,"Conjugate Gradient",MB_OK);*/

	if(vsize<eps || f2<=ftarget)
	{
	  sprintf(str,"COMPLETED : TARGET f(x)=%.5f\n",f2);
	  MessageBox(NULL,str,"Conjugate Gradient",MB_OK);
	  /*return;*/
      break;
	}

	c3=0.0;
	for(i=0; i<m; i++)
	{
	  c3+=(*(fgrad2+i))*(*(fgrad2+i));
	}
	beta=c3/c1;  /*FLETCHER-REEVES*/

    offset=0; i=0;
    while(i<m)
    {
      fixed=0;
      for(ii=0;ii<nfix;ii++)
      {
        if(offset==fp[ii])
        {
          fixed=1;
        }
      }
      if(!fixed)
      {
  	    /*UPDATE INITIAL NODE*/
        zini[offset]=z[offset]; /*ONLY Z*/
  	    /*GRADIENTS*/
	    *(u2+i)=(*(fgrad2+i))+beta*(*(u1+i));
	    *(u1+i)=*(u2+i);
	    *(fgrad1+i)=*(fgrad2+i);

        i++;
      }

      offset++;
	}
    f1=f2;  /*UJIOKA*/

    test=0.0;
    for(i=0;i<m;i++)
    {
      test+=(*(u1+i))*(*(fgrad2+i));
    }
    if(test<=0)
    {
      for(i=0;i<m;i++)
      {
         *(u1+i)=(*(fgrad2+i));    /*THE STEEPEST DESCENT METHOD*/
      }

	fprintf(ftxt,"CONJUGATE GRADIENT IS NOT DESCENT DIRECTION\n");
	sprintf(str,"THE STEEPEST DESCENT METHOD AT STEP %d",k+1);
	fprintf(ftxt,"%s\n\n",str);
    }

    sprintf(str,"hogtxt_opt%d.inp",k);
    fresult=fopen(str,"w");
    if(fresult==NULL) break;
    saveorganization(fresult,&((wdraw.childs+1)->org),
                    &((wdraw.childs+1)->vparam));
    fclose(fresult);
  }

  fclose(ftxt);

  return;
}/*conjugategradientbezier*/

#if 1
void conjugategradient(struct organ *org)     /*CONJUGATE*/
/*OPTIMIZE ORGAN WITH CONJUGATE GRADIENT.*/
{
  FILE *fout,*ftxt;
  FILE *fresult;
  char non[10],str[256],s[256];
  int i,ii,j,jj,k,m,n;
  int nnode,nelem,aelem;
  double df,f1,f2,fa,ftarget,c1,c2,c3,alpha,beta,gamma,vsize,eps;
  double *x1,*x2,*xa,*xx,*dx,dfact; /*COORDINATES*/
  double *u1,*u2,*ua,*fgrad1,*fgrad2; /*GRADIENT VECTOR*/
  double **cmtx; /*MATRIX*/
  double xi,tau,test; /*FOR LINE SEARCH(ARMIJO RULE)*/

  /*BEZIER SURFACE*/
  FILE *fbezier;
  char **data/*,str[256]="\0"*/;
  double ddata;
  int ndata;
  long int /*nnode,*/ncode;
  double *x,*y,*z;     /*CONTROLE POINTS*/
  double *xini,*yini,*zini;     /*INITIAL CONTROLE POINTS*/
  int n1,n2;           /*DEGREE*/
  int ncontrole;
  double *u,*v;      /*U-V COORDINATES*/
  struct onode node;

  /*FIXED CONTROLE POINTS*/
  int nfix;
  int *fp;    /*INDEX OF FIXED CONTROLE POINTS*/
  int fixed;
  int offset;

  ftxt=fopen("gradtest.txt","w");

  ftarget=0.1;
  gamma=1.000;
  dfact=0.01;
  eps=0.01;

  /*POLYGON01.INP : SUCCESSFUL PARAMETERS*/
  /*
  ftarget=0.4;
  gamma=0.001;
  dfact=0.05;
  eps=0.00001;
  */
  /*
  ftarget=0.9;
  gamma=0.002;
  dfact=0.05;
  eps=0.00001;
  */
  /*
  ftarget=0.9;
  gamma=0.004;
  dfact=0.1;
  eps=0.00001;
  */
  /*
  ftarget=0.9;
  gamma=0.03;
  dfact=0.1;
  eps=0.001;
  */

  /*PARAMETERS FOR LINE SEARCH(ARMIJO RULE)*/
  xi=0.001; /*0.0<=xi<=1.0*/
  tau=0.9; /*0<tau<1*/

  /*FIX CONTROLE POINTS*/
  nfix=3;    /*NUMBER OF FIXED CONTROLE POINTS*/
  fp=(int *)malloc(nfix*sizeof(int));
  fp[0]=3;    /*INDEX OF FIXED CONTROLE POINTS*/
  fp[1]=12;
  fp[2]=15;
//  fp[3]=15;

  nnode=org->nnode;
  nelem=org->nelem;

  /*CREATE INITIAL BEZIER SURFACE*/
  fbezier=fopen("bezier.txt","r");   /*BEZIER SURFACE DATA*/
  if(fbezier==NULL)
  {
    errormessage("ACCESS IMPOSSIBLE.");
    return;
  }
  fseek(fbezier,0L,SEEK_SET);

  data=fgetsbrk(fbezier,&ndata);
  n1=strtol(*(data+0),NULL,10);
  n2=strtol(*(data+1),NULL,10);
  if(nnode!=strtol(*(data+2),NULL,10)) return;

  /*CONTROLE POINTS*/
  ncontrole=(n1+1)*(n2+1);
  x=(double *)malloc((ncontrole)*sizeof(double));
  y=(double *)malloc((ncontrole)*sizeof(double));
  z=(double *)malloc((ncontrole)*sizeof(double));
  zini=(double *)malloc((ncontrole)*sizeof(double));

  for(i=0;i<ncontrole;i++)
  {
	data=fgetsbrk(fbezier,&ndata);
	if(ndata!=3) return;
    x[i]=strtod(*(data+0),NULL);
    y[i]=strtod(*(data+1),NULL);
    z[i]=strtod(*(data+2),NULL);

    for(;ndata>0;ndata--) free(*(data+ndata-1));
    free(data);

    zini[i]=z[i];
  }

  /*U-V COORDINATE OF EACH NODE*/
  u=(double *)malloc((nnode)*sizeof(double));
  v=(double *)malloc((nnode)*sizeof(double));
  for(i=0;i<nnode;i++)
  {
	data=fgetsbrk(fbezier,&ndata);
	if(ndata!=5) return;

	ncode=strtol(*(data+1),NULL,10);
	if(ncode!=(org->nodes+i)->code) return;

    ddata=strtod(*(data+3),NULL);
    *(u+i)=ddata;

    ddata=strtod(*(data+4),NULL);
    *(v+i)=ddata;

    for(;ndata>0;ndata--) free(*(data+ndata-1));
    free(data);
  }
  fclose(fbezier);

  /*BEZIER SURFACE*/
  for(i=0;i<nnode;i++)
  {
    node=*(org->nodes+i);
    beziersurfaceII(n1,n2,x,y,z,*(u+i),*(v+i),&node);

    (org->nodes+i)->d[0]=node.d[0];
    (org->nodes+i)->d[1]=node.d[1];
    (org->nodes+i)->d[2]=node.d[2];
  }

  m=ncontrole-nfix;             /*NUMBER OF PARAMETERS*/
  fgrad1=mallocdoublevector(m); /*FOR Z OF CONTROLE POINTS*/
  fgrad2=mallocdoublevector(m); /*FOR Z OF CONTROLE POINTS*/
  u1=mallocdoublevector(m);
  u2=mallocdoublevector(m);
  ua=mallocdoublevector(m);

  /*INITIAL NODE*/
  /*for(i=0; i<nnode; i++) *(ninit+i)=*(org->nodes+i);*/

  /*SAFETY AVERAGE*/
  arclmautomaticII(org);
  fout=fopen((wdraw.childs+1)->otpfile,"r");
  if(fout==NULL) return;
  readsrcanrate(fout,&arc);
  fclose(fout);
  aelem=arc.nelem;
  f1=0.0;
  for(ii=0;ii<aelem;ii++)
  {
/*
	for(jj=0;jj<4;jj++)
	{
	  if(f1<(arc.elems+ii)->srate[jj]) f1=(arc.elems+ii)->srate[jj];
	}
*/
    if((arc.elems+ii)->srate[0]<10.0) f1+=(arc.elems+ii)->srate[0];    //FLAT BAR
    else f1+=10.0;
  }
  f1/=aelem;

  sprintf(str,"%d",0);
  SetDlgItemText((wmenu.childs+2)->hwnd,ID_LAPS,str);
  SendDlgItemMessage((wmenu.childs+2)->hwnd,ID_LAPS,WM_PAINT,0,0);
  sprintf(str,"%.5f",f1);
  SetDlgItemText((wmenu.childs+2)->hwnd,ID_SAFETY,str);
  SendDlgItemMessage((wmenu.childs+2)->hwnd,ID_SAFETY,WM_PAINT,0,0);

//  sprintf(str,"Initial Max Safety = %9.5f",f1);
  sprintf(str,"Initial Average Safety = %9.5f",f1);
  fprintf(ftxt,"%s\n",str);
  MessageBox(NULL,str,"Conjugate Gradient",MB_OK);

  /*INITIAL GRADIENT*/
  fprintf(ftxt,"Initial Gradient\n");

  offset=0; i=0;
  while(i<m)
  {
    fixed=0;
    for(ii=0;ii<nfix;ii++)
    {
      if(offset==fp[ii])
      {
        fixed=1;
      }
    }
    if(!fixed)
    {
  	  for(j=0;j<ncontrole;j++)
	  {
	    if(j==offset)
	    {
          z[j]=zini[j]+dfact;                 /*UPDATE ONLY Z OF CONTROLE POINT I*/
	    }
	    else
	    {
          z[j]=zini[j];                       /*RESET Z OF OTHER CONTROLE POINTS*/
	    }
	  }
      /*BEZIER SURFACE*/
      for(ii=0;ii<nnode;ii++)
      {
        node=*(org->nodes+ii);
        beziersurfaceII(n1,n2,x,y,z,*(u+ii),*(v+ii),&node);

        (org->nodes+ii)->d[0]=node.d[0];
        (org->nodes+ii)->d[1]=node.d[1];
        (org->nodes+ii)->d[2]=node.d[2];
      }

	  /*SAFETY AVERAGE*/
	  arclmautomaticII(org);
	  fout=fopen((wdraw.childs+1)->otpfile,"r");
	  readsrcanrate(fout,&arc);
	  fclose(fout);
	  aelem=arc.nelem;
	  df=0.0;
	  for(ii=0;ii<aelem;ii++)
	  {
/*
	    for(jj=0;jj<4;jj++)
	    {
		  if(df<(arc.elems+ii)->srate[jj]) df=(arc.elems+ii)->srate[jj];
	    }
*/
        if((arc.elems+ii)->srate[0]<10.0) df+=(arc.elems+ii)->srate[0];    //FLAT BAR
        else df+=10.0;
	  }
      df/=aelem;

	  *(fgrad1+i)=gamma*(f1-df)/dfact; /*NEGATIVE GRADIENT*/
	  *(u1+i)=*(fgrad1+i);

  	  fprintf(ftxt,"CONTROLE POINT[%d] Gradient=%9.5f\n",offset,(f1-df)/dfact);

      i++;
    }

    offset++;
  }
  fprintf(ftxt,"\n");

  k=0;
  while(1)
  {
	k++;
	sprintf(str,"%d",k);
	SetDlgItemText((wmenu.childs+2)->hwnd,ID_LAPS,str);

#if 1
    offset=0; i=0;
    while(i<m)
    {
      fixed=0;
      for(ii=0;ii<nfix;ii++)
      {
        if(offset==fp[ii])
        {
          fixed=1;
        }
      }
      if(!fixed)
      {
        z[offset]=zini[offset]+(*(u1+i)); /*ONLY Z*/
        i++;
      }

      offset++;
	}

    for(ii=0;ii<nnode;ii++)
    {
      node=*(org->nodes+ii);
      beziersurfaceII(n1,n2,x,y,z,*(u+ii),*(v+ii),&node);

      (org->nodes+ii)->d[0]=node.d[0];
      (org->nodes+ii)->d[1]=node.d[1];
      (org->nodes+ii)->d[2]=node.d[2];
    }

	/*SAFETY AVERAGE*/
	arclmautomaticII(org);
	fout=fopen((wdraw.childs+1)->otpfile,"r");
	readsrcanrate(fout,&arc);
	fclose(fout);
	aelem=arc.nelem;
	fa=0.0;
	for(ii=0;ii<aelem;ii++)
	{
/*
	  for(jj=0;jj<4;jj++)
	  {
		if(fa<(arc.elems+ii)->srate[jj]) fa=(arc.elems+ii)->srate[jj];
	  }
*/
      if((arc.elems+ii)->srate[0]<10.0) fa+=(arc.elems+ii)->srate[0];    //FLAT BAR
      else fa+=10.0;
	}
    fa/=aelem;

//	fprintf(ftxt,"Step %d Max Safety = %9.5f\n",k,fa);
	fprintf(ftxt,"Step %d Average Safety = %9.5f\n",k,fa);
#endif
	/*
	sprintf(str,"Max Safety = %.5f",fa);
	MessageBox(NULL,str,"Conjugate Gradient",MB_OK);
	*/

	/*INITIAL NODE*/
    for(i=0;i<ncontrole;i++) z[i]=zini[i];

	/*GRADIENT*/
	fprintf(ftxt,"Step %d Gradient\n",k);

    offset=0; i=0;
    while(i<m)
    {
      fixed=0;
      for(ii=0;ii<nfix;ii++)
      {
        if(offset==fp[ii])
        {
          fixed=1;
        }
      }
      if(!fixed)
      {
  	    for(j=0;j<ncontrole;j++)
	    {
	      if(j==offset)
	      {
            z[j]=zini[j]+dfact;                 /*UPDATE ONLY Z OF CONTROLE POINT I*/
	      }
	      else
	      {
            z[j]=zini[j];                       /*RESET Z OF OTHER CONTROLE POINTS*/
	      }
        }
        /*BEZIER SURFACE*/
        for(ii=0;ii<nnode;ii++)
        {
          node=*(org->nodes+ii);
          beziersurfaceII(n1,n2,x,y,z,*(u+ii),*(v+ii),&node);

          (org->nodes+ii)->d[0]=node.d[0];
          (org->nodes+ii)->d[1]=node.d[1];
          (org->nodes+ii)->d[2]=node.d[2];
        }

	    /*SAFETY AVERAGE*/
	    arclmautomaticII(org);
	    fout=fopen((wdraw.childs+1)->otpfile,"r");
	    readsrcanrate(fout,&arc);
	    fclose(fout);
	    aelem=arc.nelem;
	    df=0.0;
	    for(ii=0;ii<aelem;ii++)
	    {
/*
		  for(jj=0;jj<4;jj++)
		  {
		    if(df<(arc.elems+ii)->srate[jj]) df=(arc.elems+ii)->srate[jj];
		  }
*/
          if((arc.elems+ii)->srate[0]<10.0) df+=(arc.elems+ii)->srate[0];    //FLAT BAR
          else df+=10.0;
	    }
        df/=aelem;

//	    *(ua+i)=-gamma*(fa-df)/dfact+(*(u1+i)); /*WRONG:NEGATIVE GRADIENT*/
	    *(ua+i)=-gamma*(fa-df)/dfact+(*(fgrad1+i)); /*NEGATIVE GRADIENT*/
        *(fgrad2+i)=gamma*(f1-df)/dfact;  /*NEGATIVE GRADIENT*/
        fprintf(ftxt,"CONTROLE POINT[%d] Gradient=%9.5f\n",offset,(f1-df)/dfact);

        i++;
      }

      offset++;
    }
	fprintf(ftxt,"\n");

	/*
	fprintf(stderr,"ITERATION %d\n",k);
	fprintf(stderr," {u}= %8.3f %8.3f\n",*(u1+0),*(u1+1));
	fprintf(stderr,"A{u}= %8.3f %8.3f\n",*(ua+0),*(ua+1));
	*/

	c1=0.0;
	c2=0.0;
	for(i=0; i<m; i++)
	{
	  c1+=(*(fgrad1+i))*(*(fgrad1+i));
	  c2+=(*(u1+i))*(*(ua+i));
	}
	alpha=c1/c2;
//	alpha=1.0;

	/*
	fprintf(stderr,"Alpha= %8.3f = %8.3f / %8.3f\n",alpha,c1,c2);
	*/
	sprintf(s,"BACKTRACKING LINE SEARCH BEGIN:INITIAL ALPHA= %8.3f/f1=%.5f",alpha,f1);
//	MessageBox(NULL,s,"Conjugate Gradient",MB_OK);
    fprintf(ftxt,"%s\n",s);

    while(1)   /*ARMIJO RULE*/
    {
	  vsize=0.0;
      offset=0; i=0;
      while(i<m)
      {
        fixed=0;
        for(ii=0;ii<nfix;ii++)
        {
          if(offset==fp[ii])
          {
            fixed=1;
          }
        }
        if(!fixed)
        {
          z[offset]=zini[offset]+alpha*(*(u1+i)); /*ONLY Z*/
//  	      *(fgrad2+i)=(*(fgrad1+i))-alpha*(*(ua+i));
	      vsize+=(*(fgrad2+i))*(*(fgrad2+i));

          i++;
        }

        offset++;
	  }

      /*BEZIER SURFACE*/
      for(ii=0;ii<nnode;ii++)
      {
        node=*(org->nodes+ii);
        beziersurfaceII(n1,n2,x,y,z,*(u+ii),*(v+ii),&node);

        (org->nodes+ii)->d[0]=node.d[0];
        (org->nodes+ii)->d[1]=node.d[1];
        (org->nodes+ii)->d[2]=node.d[2];
      }

	  /*SAFETY AVERAGE*/
	  arclmautomaticII(org);
	  fout=fopen((wdraw.childs+1)->otpfile,"r");
	  readsrcanrate(fout,&arc);
	  fclose(fout);
	  aelem=arc.nelem;
	  f2=0.0;
	  for(ii=0;ii<aelem;ii++)
      {
/*
	    for(jj=0;jj<4;jj++)
	    {
		  if(f2<(arc.elems+ii)->srate[jj]) f2=(arc.elems+ii)->srate[jj];
	    }
*/
        if((arc.elems+ii)->srate[0]<10.0) f2+=(arc.elems+ii)->srate[0];    //FLAT BAR
        else f2+=10.0;
	  }
      f2/=aelem;

      /*ARMIJO RULE*/
      test=0.0;
      for(i=0;i<m;i++)  test+=(*(fgrad2+i))*(*(u1+i));
      test=f2-(f1-abs(xi*alpha*test));

      if(test<=0.0 || (alpha<=0.01 && f1>f2))
      {
        sprintf(s,"LINE SEARCH:ALPHA= %8.3f/f1=%.5f/f2=%.5f/TEST= %8.3f",alpha,f1,f2,test);
	    fprintf(ftxt,"%s\n",s);
        break;
      }
      else
      {
        alpha*=tau;
      	sprintf(s,"LINE SEARCH:ALPHA= %8.3f/f1=%.5f/f2=%.5f/TEST= %8.3f",alpha,f1,f2,test);
	    fprintf(ftxt,"%s\n",s);
//        MessageBox(NULL,s,"Conjugate Gradient",MB_OK);
      }
    }
	sprintf(str,"%d",k);
	SetDlgItemText((wmenu.childs+2)->hwnd,ID_LAPS,str);
	SendDlgItemMessage((wmenu.childs+2)->hwnd,ID_LAPS,WM_PAINT,0,0);
	sprintf(str,"%.5f",f2);
	SetDlgItemText((wmenu.childs+2)->hwnd,ID_SAFETY,str);
	SendDlgItemMessage((wmenu.childs+2)->hwnd,ID_SAFETY,WM_PAINT,0,0);

//	sprintf(str,"Step=%d Max Safety=%.5f Vector Size=%.5f",k,f2,vsize);
	sprintf(str,"Step=%d Average Safety=%.5f Vector Size=%.5f",k,f2,vsize);
	fprintf(ftxt,"%s\n\n",str);
	/*MessageBox(NULL,str,"Conjugate Gradient",MB_OK);*/

	if(vsize<eps || f2<=ftarget)
	{
	  sprintf(str,"COMPLETED : TARGET f(x)=%.5f\n",f2);
	  MessageBox(NULL,str,"Conjugate Gradient",MB_OK);
	  /*return;*/
      break;
	}

	c3=0.0;
	for(i=0; i<m; i++)
	{
	  c3+=(*(fgrad2+i))*(*(fgrad2+i));
	}
	beta=c3/c1;  /*FLETCHER-REEVES*/

    offset=0; i=0;
    while(i<m)
    {
      fixed=0;
      for(ii=0;ii<nfix;ii++)
      {
        if(offset==fp[ii])
        {
          fixed=1;
        }
      }
      if(!fixed)
      {
  	    /*UPDATE INITIAL NODE*/
        zini[offset]=z[offset]; /*ONLY Z*/
  	    /*GRADIENTS*/
	    *(u2+i)=(*(fgrad2+i))+beta*(*(u1+i));
	    *(u1+i)=*(u2+i);
	    *(fgrad1+i)=*(fgrad2+i);

        i++;
      }

      offset++;
	}
    f1=f2;  /*UJIOKA*/

    test=0.0;
    for(i=0;i<m;i++)
    {
      test=(*(u1+i))*(*(fgrad2+i));
    }
    if(test<=0)
    {
      for(i=0;i<m;i++)
      {
         *(u1+i)=(*(fgrad2+i));    /*THE STEEPEST DESCENT METHOD*/
      }
    }
	fprintf(ftxt,"CONJUGATE GRADIENT IS NOT DESCENT DIRECTION\n");
	sprintf(str,"THE STEEPEST DESCENT METHOD AT STEP %d",k);
	fprintf(ftxt,"%s\n\n",str);

    sprintf(str,"hogtxt_opt%d.inp",k);
    fresult=fopen(str,"w");
    if(fresult==NULL) break;
    saveorganization(fresult,&((wdraw.childs+1)->org),
                    &((wdraw.childs+1)->vparam));
    fclose(fresult);

  }

  fclose(ftxt);

  return;
}/*conjugategradient*/
#endif

#if 0
void conjugategradient(struct organ *org)     /*CONJUGATE*/
/*OPTIMIZE ORGAN WITH CONJUGATE GRADIENT.*/
{
  FILE *fout,*ftxt;
  char non[10],str[256];
  int i,ii,j,jj,k,m,n;
  int nnode,nelem,aelem;
  double df,f1,f2,fa,ftarget,c1,c2,c3,alpha,beta,gamma,vsize,eps;
  double *x1,*x2,*xa,*xx,*dx,dfact; /*COORDINATES*/
  double *u1,*u2,*ua,*fgrad1,*fgrad2; /*GRADIENT VECTOR*/
  double **cmtx; /*MATRIX*/
  struct onode *ninit;            /*節点座標の初期値 d[0],d[1],d[2]*/

  ftxt=fopen("gradtest.txt","w");

  ftarget=0.8;
  gamma=0.001;
  dfact=0.05;            /*座標の微小変化　Δz*/
  eps=0.00001;

  /*POLYGON01.INP : SUCCESSFUL PARAMETERS*/
  /*
  ftarget=0.4;
  gamma=0.001;
  dfact=0.05;
  eps=0.00001;
  */
  /*
  ftarget=0.9;
  gamma=0.002;
  dfact=0.05;
  eps=0.00001;
  */
  /*
  ftarget=0.9;
  gamma=0.004;
  dfact=0.1;
  eps=0.00001;
  */
  /*
  ftarget=0.9;
  gamma=0.03;
  dfact=0.1;
  eps=0.001;
  */

  nnode=org->nnode;
  nelem=org->nelem;
  m=nnode;

  ninit=(struct onode *)malloc(nnode*sizeof(struct onode));
  fgrad1=mallocdoublevector(nnode); /*FOR Z OF ALL NODES*/
  fgrad2=mallocdoublevector(nnode); /*FOR Z OF ALL NODES*/
  u1=mallocdoublevector(nnode);
  u2=mallocdoublevector(nnode);
  ua=mallocdoublevector(nnode);

  /*INITIAL NODE*/
  for(i=0; i<nnode; i++) *(ninit+i)=*(org->nodes+i);

  /*MAXIMUM SAFETY*/
  arclmautomatic(org);
  fout=fopen((wdraw.childs+1)->otpfile,"r");
  if(fout==NULL) return;
  readsrcanrate(fout,&arc);
  fclose(fout);
  aelem=arc.nelem;
  f1=0.0;
  for(ii=0;ii<aelem;ii++)
  {
	for(jj=0;jj<4;jj++)
	{
	  if(f1<(arc.elems+ii)->srate[jj]) f1=(arc.elems+ii)->srate[jj];
	}
  }

  sprintf(str,"%d",0);
  SetDlgItemText((wmenu.childs+2)->hwnd,ID_LAPS,str);
  SendDlgItemMessage((wmenu.childs+2)->hwnd,ID_LAPS,WM_PAINT,0,0);
  sprintf(str,"%.5f",f1);
  SetDlgItemText((wmenu.childs+2)->hwnd,ID_SAFETY,str);
  SendDlgItemMessage((wmenu.childs+2)->hwnd,ID_SAFETY,WM_PAINT,0,0);

  sprintf(str,"Initial Max Safety = %9.5f",f1);
  fprintf(ftxt,"%s\n",str);
  MessageBox(NULL,str,"Conjugate Gradient",MB_OK);

  /*INITIAL GRADIENT*/
  fprintf(ftxt,"Initial Gradient\n");
  for(i=0;i<m;i++)
  {
	for(j=0;j<m;j++)       /*z座標を微小量(dfact)変化*/
	{
	  if(j==i)
	  {
		(org->nodes+j)->d[2]=((ninit+j)->d[2])+dfact; /*UPDATE ONLY Z OF NODE I*/
	  }
	  else
	  {
		(org->nodes+j)->d[2]=((ninit+j)->d[2]); /*RESET Z OF OTHER NODES*/
	  }
	}

	/*MAXIMUM SAFETY*/
	arclmautomatic(org);
	fout=fopen((wdraw.childs+1)->otpfile,"r");
	readsrcanrate(fout,&arc);
	fclose(fout);
	aelem=arc.nelem;
	df=0.0;
	for(ii=0;ii<aelem;ii++)
	{
	  for(jj=0;jj<4;jj++)
	  {
		if(df<(arc.elems+ii)->srate[jj]) df=(arc.elems+ii)->srate[jj];
	  }
	}

	*(fgrad1+i)=gamma*(f1-df)/dfact; /*NEGATIVE GRADIENT*/
	*(u1+i)=*(fgrad1+i);

	fprintf(ftxt,"Node %d Gradient=%9.5f\n",(org->nodes+i)->code,(f1-df)/dfact);
  }
  fprintf(ftxt,"\n");

  k=0;
  while(1)
  {
	k++;
	sprintf(str,"%d",k);
	SetDlgItemText((wmenu.childs+2)->hwnd,ID_LAPS,str);

	for(i=0;i<m;i++)
	{
	  (org->nodes+i)->d[2]=((ninit+i)->d[2])+(*(u1+i)); /*ONLY Z*/
	}

	/*MAXIMUM SAFETY*/
	arclmautomatic(org);
	fout=fopen((wdraw.childs+1)->otpfile,"r");
	readsrcanrate(fout,&arc);
	fclose(fout);
	aelem=arc.nelem;
	fa=0.0;
	for(ii=0;ii<aelem;ii++)
	{
	  for(jj=0;jj<4;jj++)
	  {
		if(fa<(arc.elems+ii)->srate[jj]) fa=(arc.elems+ii)->srate[jj];
	  }
	}

	fprintf(ftxt,"Step %d Max Safety = %9.5f\n",k,fa);

	/*
	sprintf(str,"Max Safety = %.5f",fa);
	MessageBox(NULL,str,"Conjugate Gradient",MB_OK);
	*/

	/*INITIAL NODE*/
	for(i=0;i<nnode;i++) *(ninit+i)=*(org->nodes+i);

	/*GRADIENT*/
	fprintf(ftxt,"Step %d Gradient\n",k);
	for(i=0;i<m;i++)
	{
	  for(j=0;j<m;j++)
	  {
		if(j==i)
		{
		  (org->nodes+j)->d[2]=((ninit+j)->d[2])+dfact; /*ONLY Z*/
		}
		else
		{
		  (org->nodes+j)->d[2]=((ninit+j)->d[2]); /*ONLY Z*/
		}
	  }

	  /*MAXIMUM SAFETY*/
	  arclmautomatic(org);
	  fout=fopen((wdraw.childs+1)->otpfile,"r");
	  readsrcanrate(fout,&arc);
	  fclose(fout);
	  aelem=arc.nelem;
	  df=0.0;
	  for(ii=0;ii<aelem;ii++)
	  {
		for(jj=0;jj<4;jj++)
		{
		  if(df<(arc.elems+ii)->srate[jj]) df=(arc.elems+ii)->srate[jj];
		}
	  }

//	  *(ua+i)=-gamma*(fa-df)/dfact+(*(u1+i)); /*WRONG:NEGATIVE GRADIENT*/
	  *(ua+i)=-gamma*(fa-df)/dfact+(*(fgrad1+i)); /*NEGATIVE GRADIENT*/

	  fprintf(ftxt,"Node %d Gradient=%9.5f\n",(org->nodes+i)->code,(fa-df)/dfact);
	}
	fprintf(ftxt,"\n");

	/*
	fprintf(stderr,"ITERATION %d\n",k);
	fprintf(stderr," {u}= %8.3f %8.3f\n",*(u1+0),*(u1+1));
	fprintf(stderr,"A{u}= %8.3f %8.3f\n",*(ua+0),*(ua+1));
	*/

	c1=0.0;
	c2=0.0;
	for(i=0; i<m; i++)
	{
	  c1+=(*(fgrad1+i))*(*(fgrad1+i));
	  c2+=(*(u1+i))*(*(ua+i));
	}
	alpha=c1/c2;

	/*
	fprintf(stderr,"Alpha= %8.3f = %8.3f / %8.3f\n",alpha,c1,c2);
	*/

	vsize=0.0;
	for(i=0;i<m;i++)
	{
	  (org->nodes+i)->d[2]=((ninit+i)->d[2])+alpha*(*(u1+i)); /*ONLY Z*/

	  /*
	  fprintf(stderr,"x = x + dx : %8.3f = %8.3f + %8.3f x %8.3f\n",*(x2+i),*(x1+i),alpha,*(u1+i));
	  */

	  *(fgrad2+i)=(*(fgrad1+i))-alpha*(*(ua+i));
	  vsize+=(*(fgrad2+i))*(*(fgrad2+i));
	}

	/*MAXIMUM SAFETY*/
	arclmautomatic(org);
	fout=fopen((wdraw.childs+1)->otpfile,"r");
	readsrcanrate(fout,&arc);
	fclose(fout);
	aelem=arc.nelem;
	f2=0.0;
	for(ii=0;ii<aelem;ii++)
	{
	  for(jj=0;jj<4;jj++)
	  {
		if(f2<(arc.elems+ii)->srate[jj]) f2=(arc.elems+ii)->srate[jj];
	  }
	}

	sprintf(str,"%d",k);
	SetDlgItemText((wmenu.childs+2)->hwnd,ID_LAPS,str);
	SendDlgItemMessage((wmenu.childs+2)->hwnd,ID_LAPS,WM_PAINT,0,0);
	sprintf(str,"%.5f",f2);
	SetDlgItemText((wmenu.childs+2)->hwnd,ID_SAFETY,str);
	SendDlgItemMessage((wmenu.childs+2)->hwnd,ID_SAFETY,WM_PAINT,0,0);

	sprintf(str,"Step=%d Max Safety=%.5f Vector Size=%.5f",k,f2,vsize);
	fprintf(ftxt,"%s\n\n",str);
	/*MessageBox(NULL,str,"Conjugate Gradient",MB_OK);*/

	if(vsize<eps || f2<=ftarget)
	{
	  sprintf(str,"COMPLETED : TARGET f(x)=%.5f\n",f2);
	  MessageBox(NULL,str,"Conjugate Gradient",MB_OK);
	  return;
	}

	c3=0.0;
	for(i=0; i<m; i++)
	{
	  c3+=(*(fgrad2+i))*(*(fgrad2+i));
	}
	beta=c3/c1;

	for(i=0; i<m; i++)
	{
	  /*UPDATE INITIAL NODE*/
	  *(ninit+i)=*(org->nodes+i);

	  /*GRADIENTS*/
	  *(u2+i)=(*(fgrad2+i))+beta*(*(u1+i));
	  *(u1+i)=*(u2+i);
	  *(fgrad1+i)=*(fgrad2+i);
	}
  }

  fclose(ftxt);

  return;
}/*conjugategradient*/
#endif

#if 0
void conjugategradient(struct organ *org) /*CONJUGATE*/
/*OPTIMIZE ORGAN WITH CONJUGATE GRADIENT.*/
{
  FILE *fout,*ftxt;
  char non[10],s[256],str[8000];
  int i,ii,j,jj,k,kk,m,n;
  int nnode,nelem,aelem;
  double df,f1,f2,fa,ftarget,c1,c2,c3,alpha,beta,gamma,vsize,eps;
  double *x1,*x2,*xa,*xx,*dx,dfact; /*COORDINATES*/
  double *u1,*u2,*ua,*fgrad1,*fgrad2; /*GRADIENT VECTOR*/
  double **cmtx; /*MATRIX*/
 // double p[2][3];
  struct onode *ninit;

  ftxt=fopen("gradtest.txt","w");

  ftarget=0.99;   /*最終安全率*/
  gamma=0.001;    /*ステップ幅の係数*/
  dfact=0.01;     /*動く幅m*/
  eps=0.00001;    /*ベクトル長さ*/

 /* p[0][0]=0.015;
  p[0][1]=0.000;
  p[0][2]=0.020;

  p[1][0]=0.030;
  p[1][1]=0.000;
  p[1][2]=0.045;  */

  /* ftarget=0.8;
  gamma=0.001;
  dfact=0.05;
  eps=0.00001; */


  /*POLYGON01.INP : SUCCESSFUL PARAMETERS*/
  /*
  ftarget=0.4;
  gamma=0.001;
  dfact=0.05;
  eps=0.00001;*/

  /*
  ftarget=0.9;
  gamma=0.002;
  dfact=0.05;
  eps=0.00001;*/

  /*
  ftarget=0.9;
  gamma=0.004;
  dfact=0.1;
  eps=0.00001;
  */
  /*
  ftarget=0.9;
  gamma=0.03;
  dfact=0.1;
  eps=0.001;*/


  nnode=org->nnode;
  nelem=org->nelem;
  m=3*nnode;

  ninit=(struct onode *)malloc(nnode*sizeof(struct onode));
 // fgrad1=mallocdoublevector(nnode); /*FOR Z OF ALL NODES*/
 // fgrad2=mallocdoublevector(nnode); /*FOR Z OF ALL NODES*/
 // u1=mallocdoublevector(nnode);
 // u2=mallocdoublevector(nnode);
 // ua=mallocdoublevector(nnode);

  fgrad1=mallocdoublevector(3*nnode); /*FOR XYZ OF ALL NODES*/
  fgrad2=mallocdoublevector(3*nnode); /*FOR XYZ OF ALL NODES*/
  u1=mallocdoublevector(3*nnode);
  u2=mallocdoublevector(3*nnode);
  ua=mallocdoublevector(3*nnode);


  /*INITIAL NODE*/
  for(i=0; i<nnode; i++) *(ninit+i)=*(org->nodes+i);

  /*MAXIMUM SAFETY*/
  arclmautomatic(org);
  fout=fopen((wdraw.childs+1)->otpfile,"r");
  if(fout==NULL) return;
  readsrcanrate(fout,&arc);
  fclose(fout);
  aelem=arc.nelem;
  f1=0.0;
  for(ii=0;ii<aelem;ii++)
  {
	for(jj=0;jj<4;jj++)
	{
	  if(f1<(arc.elems+ii)->srate[jj]) f1=(arc.elems+ii)->srate[jj];
	}
  }

  sprintf(str,"%d",0);
  SetDlgItemText((wmenu.childs+2)->hwnd,ID_LAPS,str);
  SendDlgItemMessage((wmenu.childs+2)->hwnd,ID_LAPS,WM_PAINT,0,0);
  sprintf(str,"%.5f",f1);
  SetDlgItemText((wmenu.childs+2)->hwnd,ID_SAFETY,str);
  SendDlgItemMessage((wmenu.childs+2)->hwnd,ID_SAFETY,WM_PAINT,0,0);

  sprintf(str,"Initial Max Safety = %9.5f",f1);
  fprintf(ftxt,"%s\n",str);
  MessageBox(NULL,str,"Conjugate Gradient",MB_OK);

  /*INITIAL GRADIENT*/
  fprintf(ftxt,"Initial Gradient\n");
  for(i=0;i<nnode;i++)
  {
	for(kk=0;kk<3;kk++)   /*XYZ方向へのループ*/
	{
	sprintf(str,"\0");

	for(j=0;j<nnode;j++)  /*m→nnodeに変更*/
	{
	  if(j==i)
	  {
		//(org->nodes+j)->d[kk]=((ninit+j)->d[kk])+dfact; /*UPDATE XYZ OF NODE I*/

		if(kk!=0) (org->nodes+j)->d[0]=((ninit+j)->d[0]); /*RESET XYZ OF OTHER NODES*/
		if(kk!=1) (org->nodes+j)->d[1]=((ninit+j)->d[1]); /*RESET XYZ OF OTHER NODES*/
		if(kk!=2) (org->nodes+j)->d[2]=((ninit+j)->d[2]); /*RESET XYZ OF OTHER NODES*/

		(org->nodes+j)->d[kk]=((ninit+j)->d[kk])+dfact; /*UPDATE XYZ OF NODE I*/



/*sprintf(s,"NODE %d XYZ={ %9.5f %9.5f %9.5f} i=%d j=%d kk=%d\n",
		(org->nodes+j)->code,
		(org->nodes+j)->d[0],
		(org->nodes+j)->d[1],
		(org->nodes+j)->d[2],i,j,kk);
strcat(str,s);   */

	  }
	  else
	  {
		//(org->nodes+j)->d[kk]=((ninit+j)->d[kk]); /*RESET XYZ OF OTHER NODES*/

		(org->nodes+j)->d[0]=((ninit+j)->d[0]); /*RESET XYZ OF OTHER NODES*/
		(org->nodes+j)->d[1]=((ninit+j)->d[1]); /*RESET XYZ OF OTHER NODES*/
		(org->nodes+j)->d[2]=((ninit+j)->d[2]); /*RESET XYZ OF OTHER NODES*/

/*sprintf(s,"NODE %d XYZ={ %9.5f %9.5f %9.5f} i=%d j=%d kk=%d\n",
		(org->nodes+j)->code,
		(org->nodes+j)->d[0],
		(org->nodes+j)->d[1],
		(org->nodes+j)->d[2],i,j,kk);
strcat(str,s);            */

	  }
	}
//MessageBox(NULL,str,"Conjugate Gradient",MB_OK);

	/*MAXIMUM SAFETY*/
	arclmautomatic(org);
	fout=fopen((wdraw.childs+1)->otpfile,"r");
	readsrcanrate(fout,&arc);
	fclose(fout);
	aelem=arc.nelem;
	df=0.0;
	for(ii=0;ii<aelem;ii++)
	{
	  for(jj=0;jj<4;jj++)
	  {
		if(df<(arc.elems+ii)->srate[jj]) df=(arc.elems+ii)->srate[jj];
	  }
	}

	*(fgrad1+i*3+kk)=gamma*(f1-df)/dfact; /*NEGATIVE GRADIENT*/
	*(u1+i*3+kk)=*(fgrad1+i*3+kk);

	fprintf(ftxt,"Node %d Coord %d Gradient=%9.5f\n",(org->nodes+i)->code,kk+1,(f1-df)/dfact);

//MessageBox(NULL,"Pass1","Conjugate Gradient",MB_OK);
   }

  }
  fprintf(ftxt,"\n");

//MessageBox(NULL,"Pass2","Conjugate Gradient",MB_OK);

  k=0;
  while(1)
  {
	k++;
	sprintf(str,"%d",k);
	SetDlgItemText((wmenu.childs+2)->hwnd,ID_LAPS,str);




	for(i=0;i<nnode;i++)
	{
	  (org->nodes+i)->d[0]=((ninit+i)->d[0])+(*(u1+i*3+0)); /* X */
	  (org->nodes+i)->d[1]=((ninit+i)->d[1])+(*(u1+i*3+1)); /* Y */
	  (org->nodes+i)->d[2]=((ninit+i)->d[2])+(*(u1+i*3+2)); /* Z */
	}

	/*MAXIMUM SAFETY*/
	arclmautomatic(org);
	fout=fopen((wdraw.childs+1)->otpfile,"r");
	readsrcanrate(fout,&arc);
	fclose(fout);
	aelem=arc.nelem;
	fa=0.0;
	for(ii=0;ii<aelem;ii++)
	{
	  for(jj=0;jj<4;jj++)
	  {
		if(fa<(arc.elems+ii)->srate[jj]) fa=(arc.elems+ii)->srate[jj];
	  }
	}

	fprintf(ftxt,"Step %d Max Safety = %9.5f\n",k,fa);

	/*
	sprintf(str,"Max Safety = %.5f",fa);
	MessageBox(NULL,str,"Conjugate Gradient",MB_OK);
	*/

	/*INITIAL NODE*/
	for(i=0;i<nnode;i++) *(ninit+i)=*(org->nodes+i);

	/*GRADIENT*/
	fprintf(ftxt,"Step %d Gradient\n",k);

	for(i=0;i<nnode;i++)
	{
	for(kk=0;kk<3;kk++)
	{
	  for(j=0;j<nnode;j++)
	  {
		if(j==i)
		{
		  if(kk!=0) (org->nodes+j)->d[0]=((ninit+j)->d[0]); /*RESET XYZ OF OTHER NODES*/
		  if(kk!=1) (org->nodes+j)->d[1]=((ninit+j)->d[1]); /*RESET XYZ OF OTHER NODES*/
		  if(kk!=2) (org->nodes+j)->d[2]=((ninit+j)->d[2]); /*RESET XYZ OF OTHER NODES*/

		(org->nodes+j)->d[kk]=((ninit+j)->d[kk])+dfact; /*UPDATE XYZ OF NODE I*/
		 // (org->nodes+j)->d[kk]=((ninit+j)->d[kk])+dfact; /*XYZ*/

		}
		else
		{
		  //(org->nodes+j)->d[kk]=((ninit+j)->d[kk]); /*XYZ*/
		  (org->nodes+j)->d[0]=((ninit+j)->d[0]); /*RESET XYZ OF OTHER NODES*/
		  (org->nodes+j)->d[1]=((ninit+j)->d[1]); /*RESET XYZ OF OTHER NODES*/
		  (org->nodes+j)->d[2]=((ninit+j)->d[2]); /*RESET XYZ OF OTHER NODES*/
		}
	  }

	  /*MAXIMUM SAFETY*/
	  arclmautomatic(org);
	  fout=fopen((wdraw.childs+1)->otpfile,"r");
	  readsrcanrate(fout,&arc);
	  fclose(fout);
	  aelem=arc.nelem;
	  df=0.0;
	  for(ii=0;ii<aelem;ii++)
	  {
		for(jj=0;jj<4;jj++)
		{
		  if(df<(arc.elems+ii)->srate[jj]) df=(arc.elems+ii)->srate[jj];
		}
	  }

	  *(ua+i*3+kk)=-gamma*(fa-df)/dfact+(*(u1+i*3+kk)); /*NEGATIVE GRADIENT*/

	  fprintf(ftxt,"Node %d Coord %d Gradient=%9.5f\n",(org->nodes+i)->code,kk+1,(fa-df)/dfact);
	}
  }
	fprintf(ftxt,"\n");


	/*
	fprintf(stderr,"ITERATION %d\n",k);
	fprintf(stderr," {u}= %8.3f %8.3f\n",*(u1+0),*(u1+1));
	fprintf(stderr,"A{u}= %8.3f %8.3f\n",*(ua+0),*(ua+1));
	*/

	c1=0.0;
	c2=0.0;
	for(i=0; i<m; i++)
	{
	  c1+=(*(fgrad1+i))*(*(fgrad1+i));
	  c2+=(*(u1+i))*(*(ua+i));
	}
	alpha=c1/c2;

	/*
	fprintf(stderr,"Alpha= %8.3f = %8.3f / %8.3f\n",alpha,c1,c2);
	*/

	vsize=0.0;
	for(i=0;i<nnode;i++)

	{
	 for(kk=0;kk<3;kk++)
	 {
	  (org->nodes+i)->d[kk]=((ninit+i)->d[kk])+alpha*(*(u1+i*3+kk)); /*ONLY Z*/

	  /*
	  fprintf(stderr,"x = x + dx : %8.3f = %8.3f + %8.3f x %8.3f\n",*(x2+i),*(x1+i),alpha,*(u1+i));
	  */

	  *(fgrad2+i*3+kk)=(*(fgrad1+i*3+kk))-alpha*(*(ua+i*3+kk));
	  vsize+=(*(fgrad2+i*3+kk))*(*(fgrad2+i*3+kk));
	 }
	}


	/*MAXIMUM SAFETY*/
	arclmautomatic(org);
	fout=fopen((wdraw.childs+1)->otpfile,"r");
	readsrcanrate(fout,&arc);
	fclose(fout);
	aelem=arc.nelem;
	f2=0.0;
	for(ii=0;ii<aelem;ii++)
	{
	  for(jj=0;jj<4;jj++)
	  {
		if(f2<(arc.elems+ii)->srate[jj]) f2=(arc.elems+ii)->srate[jj];
	  }
	}

	sprintf(str,"%d",k);
	SetDlgItemText((wmenu.childs+2)->hwnd,ID_LAPS,str);
	SendDlgItemMessage((wmenu.childs+2)->hwnd,ID_LAPS,WM_PAINT,0,0);
	sprintf(str,"%.5f",f2);
	SetDlgItemText((wmenu.childs+2)->hwnd,ID_SAFETY,str);
	SendDlgItemMessage((wmenu.childs+2)->hwnd,ID_SAFETY,WM_PAINT,0,0);

	sprintf(str,"Step=%d Max Safety=%.5f Vector Size=%.5f",k,f2,vsize);
	fprintf(ftxt,"%s\n\n",str);
	/*MessageBox(NULL,str,"Conjugate Gradient",MB_OK);*/

	if(vsize<eps || f2<=ftarget)
	{
	  sprintf(str,"COMPLETED : TARGET f(x)=%.5f\n",f2);
	  MessageBox(NULL,str,"Conjugate Gradient",MB_OK);
	  return;
	}

	c3=0.0;
	for(i=0; i<m; i++)
	{
	  c3+=(*(fgrad2+i))*(*(fgrad2+i));
	}
	beta=c3/c1;

	for(i=0; i<nnode; i++)
	{
	  /*UPDATE INITIAL NODE*/
	  *(ninit+i)=*(org->nodes+i);
	}
	for(i=0; i<m; i++)
	{
	  /*GRADIENTS*/
	  *(u2+i)=(*(fgrad2+i))+beta*(*(u1+i));
	  *(u1+i)=*(u2+i);
	  *(fgrad1+i)=*(fgrad2+i);
	}
  }

  fclose(ftxt);

  return;
}/*conjugategradient*/
#endif

void conjugategradientxz(struct organ *org) /*Coded by Fujimoto*/
/*OPTIMIZE ORGAN WITH CONJUGATE GRADIENT.*/
{
  FILE *fout,*ftxt;
  char non[10],str[256];
  int i,ii,j,jj,k,m,n;
  int p;
  int nnode,nelem,aelem;
  double df,f1,f2,fa,ftarget,c1,c2,c3,alpha,beta,gamma,vsize,eps;
  //double *x1,*x2,*xa,*xx,*dx; /*COORDINATES*/
  double dfact;
  double **u1,**u2,**ua,**fgrad1,**fgrad2; /*GRADIENT VECTOR*/
  //double **cmtx; /*MATRIX*/
  struct onode *ninit; /*node-initial*/

  ftxt=fopen("gradtest.txt","w");

  //for CG01.inp
  ftarget=0.99;
  gamma=0.012;
  dfact=0.001;
  eps=0.00001;

  /*
  ftarget=0.8;
  gamma=0.001;
  dfact=0.05;
  eps=0.00001;
  */

  /*POLYGON01.INP*/
  /*
  ftarget=0.4;
  gamma=0.001;
  dxfact=0.05;
  eps=0.00001;
  */
  /*
  ftarget=0.9;
  gamma=0.002;
  dxfact=0.05;
  eps=0.00001;
  */
  /*
  ftarget=0.9;
  gamma=0.004;
  dxfact=0.1;
  eps=0.00001;
  */
  /*
  ftarget=0.9;
  gamma=0.03;
  dxfact=0.1;
  eps=0.001;
  */

  nnode=org->nnode;
  //nelem=org->nelem;
  m=nnode;

  ninit=(struct onode *)malloc(nnode*sizeof(struct onode));

  //matrixだと容量食いすぎる
  fgrad1=mallocdoublematrixxyz(nnode,3); /*FOR XYZ OF ALL NODES*/
  fgrad2=mallocdoublematrixxyz(nnode,3); /*FOR XYZ OF ALL NODES*/
  u1=mallocdoublematrixxyz(nnode,3);
  u2=mallocdoublematrixxyz(nnode,3);
  ua=mallocdoublematrixxyz(nnode,3);

  /*INITIAL NODE*/
  for(i=0; i<nnode; i++) *(ninit+i)=*(org->nodes+i);

  /*MAXIMUM SAFETY*/
  arclmautomatic(org);
  /*
  n=strcspn((wdraw.childs+1)->sctfile,".");
  strncpy(str,(wdraw.childs+1)->sctfile,n);
  str[n]='\0';
  strcpy((wdraw.childs+1)->otpfile,str);
  strcat((wdraw.childs+1)->otpfile,".rat");
  SetDlgItemText((wmenu.childs+2)->hwnd,ID_OUTPUTFILE,
				 (wdraw.childs+1)->otpfile);
  */
  fout=fopen((wdraw.childs+1)->otpfile,"r");
  if(fout==NULL) return;
  readsrcanrate(fout,&arc);
  fclose(fout);
  aelem=arc.nelem;
  f1=0.0;
  for(ii=0;ii<aelem;ii++)
  {
	for(jj=0;jj<4;jj++)
	{
	  if(f1<(arc.elems+ii)->srate[jj]) f1=(arc.elems+ii)->srate[jj];
	}
  }

  sprintf(str,"%d",0);
  SetDlgItemText((wmenu.childs+2)->hwnd,ID_LAPS,str);
  SendDlgItemMessage((wmenu.childs+2)->hwnd,ID_LAPS,WM_PAINT,0,0);
  sprintf(str,"%.5f",f1);
  SetDlgItemText((wmenu.childs+2)->hwnd,ID_SAFETY,str);
  SendDlgItemMessage((wmenu.childs+2)->hwnd,ID_SAFETY,WM_PAINT,0,0);

  sprintf(str,"Initial Max Safety = %9.5f",f1);
  fprintf(ftxt,"%s\n",str);
  MessageBox(NULL,str,"Conjugate Gradient",MB_OK);

  /*INITIAL GRADIENT*/
  fprintf(ftxt,"Initial Gradient\n");
  for(i=0;i<m;i++)
  {
	for(p=0;p<3;p++)//XZに限定
	{
		if(p!=1)
		{
			for(j=0;j<m;j++)
			{
			  if((p!=2||(org->nodes+j)->d[p]!=0)&&((org->nodes+i)->code!=115))
			  {
				  if(j==i)(org->nodes+j)->d[p]=((ninit+j)->d[p])+dfact;
				  else(org->nodes+j)->d[p]=((ninit+j)->d[p]);
			  }
			}


			/*MAXIMUM SAFETY*/
			arclmautomatic(org);
			fout=fopen((wdraw.childs+1)->otpfile,"r");
			readsrcanrate(fout,&arc);
			fclose(fout);
			aelem=arc.nelem;
			df=0.0;
			for(ii=0;ii<aelem;ii++)
			{
			  for(jj=0;jj<4;jj++)
			  {
				if(df<(arc.elems+ii)->srate[jj]) df=(arc.elems+ii)->srate[jj];
			  }
			}
			*(*(fgrad1+i)+p)=gamma*(f1-df)/dfact; /*NEGATIVE GRADIENT*/
			*(*(u1+i)+p)=*(*(fgrad1+i)+p);

			fprintf(ftxt,"Node %d Dim %d Gradient=%9.5f\n",(org->nodes+i)->code,p,(f1-df)/dfact);

			(org->nodes+i)->d[p]=((ninit+i)->d[p]);
		}
	}
  }
  fprintf(ftxt,"\n");

  /*
  fprintf(stderr,"INITIAL CONDITION\n");
  fprintf(stderr,"x   = %8.3f %8.3f\n",*(x1+0),*(x1+1));
  fprintf(stderr,"f(x)= %8.3f\n",f1);
  fprintf(stderr,"grad f(x)= %8.3f %8.3f\n",*(fgrad1+0),*(fgrad1+1));
  fprintf(stderr,"{u}      = %8.3f %8.3f\n",*(u1+0),*(u1+1));
  fprintf(fout,"STEP 0 f %8.3f {x1,x2} %8.3f %8.3f\n",f1,*(x1+0),*(x1+1));
  gets(non);
  */

  k=0;
  while(1)
  {
	/*INITIAL GRADIENT RESULT*/
	k++;
	sprintf(str,"%d",k);
	SetDlgItemText((wmenu.childs+2)->hwnd,ID_LAPS,str);

	for(i=0;i<m;i++)
	{
	  for(p=0;p<3;p++)
	  {
		if(p!=1)
		{
		  (org->nodes+i)->d[p]=((ninit+i)->d[p])+(*(*(u1+i)+p));
		}
	  }
	}

	/*MAXIMUM SAFETY*/
	arclmautomatic(org);
	fout=fopen((wdraw.childs+1)->otpfile,"r");
	readsrcanrate(fout,&arc);
	fclose(fout);
	aelem=arc.nelem;
	fa=0.0;
	for(ii=0;ii<aelem;ii++)
	{
	  for(jj=0;jj<4;jj++)
	  {
		if(fa<(arc.elems+ii)->srate[jj]) fa=(arc.elems+ii)->srate[jj];
	  }
	}

	fprintf(ftxt,"Step %d Max Safety = %9.5f\n",k,fa);

	/*
	sprintf(str,"Max Safety = %.5f",fa);
	MessageBox(NULL,str,"Conjugate Gradient",MB_OK);
	*/

	/*INITIAL NODE*/
	for(i=0;i<nnode;i++) *(ninit+i)=*(org->nodes+i);

	/*GRADIENT*/
	fprintf(ftxt,"Step %d Gradient\n",k);
	for(i=0;i<m;i++)
	{
		for(p=0;p<3;p++)//XYに限定
		{
			if(p!=1)
			{
				for(j=0;j<m;j++)
				{
					if((p!=2||(org->nodes+j)->d[p]!=0)&&((org->nodes+i)->code!=115))
					{
						if(j==i){(org->nodes+j)->d[p]=((ninit+j)->d[p])+dfact;}
						else{(org->nodes+j)->d[p]=((ninit+j)->d[p]);}
					}

				}

				/*MAXIMUM SAFETY*/
				arclmautomatic(org);
				fout=fopen((wdraw.childs+1)->otpfile,"r");
				readsrcanrate(fout,&arc);
				fclose(fout);
				aelem=arc.nelem;
				df=0.0;
				for(ii=0;ii<aelem;ii++)
				{
					for(jj=0;jj<4;jj++)
					{
					  if(df<(arc.elems+ii)->srate[jj]) df=(arc.elems+ii)->srate[jj];
					}
				}

//				*(*(ua+i)+p)=-gamma*(fa-df)/dfact+(*(*(u1+i)+p)); /*WRONG:NEGATIVE GRADIENT*/
				*(*(ua+i)+p)=-gamma*(fa-df)/dfact+(*(*(fgrad1+i)+p)); /*NEGATIVE GRADIENT*/

				fprintf(ftxt,"Node %d Dim %d Gradient=%9.5f\n",(org->nodes+i)->code,p,(fa-df)/dfact);

				(org->nodes+i)->d[p]=((ninit+i)->d[p]);
			}
		}
	}
	fprintf(ftxt,"\n");

	/*
	fprintf(stderr,"ITERATION %d\n",k);
	fprintf(stderr," {u}= %8.3f %8.3f\n",*(u1+0),*(u1+1));
	fprintf(stderr,"A{u}= %8.3f %8.3f\n",*(ua+0),*(ua+1));
	*/

	c1=0.0;
	c2=0.0;
	for(i=0; i<m; i++)
	{
		for(p=0;p<3;p++)
		{
		if(p!=1)
		{
		  c1+=(*(*(fgrad1+i)+p))*(*(*(fgrad1+i)+p));
		  c2+=(*(*(u1+i)+p))*(*(*(ua+i)+p));
			/*if (c2==0) c2=0.000001;//c2これで十分小さいか不明
			else c2=c2;*/
		}
		}
	}
	alpha=c1/c2;

	/*
	fprintf(stderr,"Alpha= %8.3f = %8.3f / %8.3f\n",alpha,c1,c2);
	*/

	vsize=0.0;
	for(i=0;i<m;i++)
	{
		for(p=0;p<3;p++)
		{
		if(p!=1)
		{
		  (org->nodes+i)->d[p]=((ninit+i)->d[p])+alpha*(*(*(u1+i)+p));

		  /*
		  fprintf(stderr,"x = x + dx : %8.3f = %8.3f + %8.3f x %8.3f\n",*(x2+i),*(x1+i),alpha,*(u1+i));
		  */

		  *(*(fgrad2+i)+p)=(*(*(fgrad1+i)+p))-alpha*(*(*(ua+i)+p));
		  vsize+=(*(*(fgrad2+i)+p))*(*(*(fgrad2+i)+p));
		}
		}
	}

	/*MAXIMUM SAFETY*/
	arclmautomatic(org);
	fout=fopen((wdraw.childs+1)->otpfile,"r");
	readsrcanrate(fout,&arc);
	fclose(fout);
	aelem=arc.nelem;
	f2=0.0;
	for(ii=0;ii<aelem;ii++)
	{
	  for(jj=0;jj<4;jj++)
	  {
		if(f2<(arc.elems+ii)->srate[jj]) f2=(arc.elems+ii)->srate[jj];
	  }
	}

	sprintf(str,"%d",k);
	SetDlgItemText((wmenu.childs+2)->hwnd,ID_LAPS,str);
	SendDlgItemMessage((wmenu.childs+2)->hwnd,ID_LAPS,WM_PAINT,0,0);
	sprintf(str,"%.5f",f2);
	SetDlgItemText((wmenu.childs+2)->hwnd,ID_SAFETY,str);
	SendDlgItemMessage((wmenu.childs+2)->hwnd,ID_SAFETY,WM_PAINT,0,0);

	sprintf(str,"Step=%d Max Safety=%.5f Vector Size=%.5f",k,f2,vsize);
	fprintf(ftxt,"%s\n\n",str);
	/*MessageBox(NULL,str,"Conjugate Gradient",MB_OK);*/

	/*
	fprintf(stderr,"x    = %8.3f %8.3f\n",*(x2+0),*(x2+1));
	fprintf(stderr,"f(x) = %8.3f\n",f2);
	fprintf(stderr,"grad f(x)   = %8.3f %8.3f\n",*(fgrad2+0),*(fgrad2+1));
	fprintf(stderr,"VECTOR SIZE = %9.5f\n",vsize);
	fprintf(fout,"STEP %d f %8.3f {x1,x2} %8.3f %8.3f\n",k,f2,*(x2+0),*(x2+1));
	gets(non);
	*/

	if(vsize<eps || f2<=ftarget)
	{
	  sprintf(str,"COMPLETED : TARGET f(x)=%.5f\n",f2);
	  MessageBox(NULL,str,"Conjugate Gradient",MB_OK);
	  return;
	}


	c3=0.0;
	for(i=0; i<m; i++)
	{
	  for(p=0;p<3;p++)
	  {
		if(p!=1)
		{
		  c3+=(*(*(fgrad2+i)+p))*(*(*(fgrad2+i)+p));
		}
	  }
	}
	beta=c3/c1;

	/*繰り返しのため引き継ぎ*/
	for(i=0; i<m; i++)
	{
	  /*INITIAL NODE*/
	  *(ninit+i)=*(org->nodes+i);

	  for(p=0;p<3;p++)
	  {
		if(p!=1)
		{
		  *(*(u2+i)+p)=(*(*(fgrad2+i)+p))+beta*(*(*(u1+i)+p));
		  *(*(u1+i)+p)=*(*(u2+i)+p);
		  *(*(fgrad1+i)+p)=*(*(fgrad2+i)+p);
		}
	  }
	}
  }

  //fclose(ftxt);

  return;
}/*conjugategradientxz*/

void conjugategradientcurve(struct organ *org)
/*OPTIMIZE ORGAN WITH CONJUGATE GRADIENT.*/
{
  FILE *fout,*ftxt;
  char non[10],str[256];
  int i,ii,j,jj,k,m,n;
  int p,q;
  int nnode,nelem,aelem;
  double df,f1,f2,fa,ftarget,c1,c2,c3,alpha,beta,gamma,vsize,eps;
  double xi2,yi2,zi2,xj2,yj2,zj2;
  double delta,d1,d2;
  //double *x1,*x2,*xa,*xx,*dx; /*COORDINATES*/
  double dfact;
  double **u1,**u2,**ua,**fgrad1,**fgrad2; /*GRADIENT VECTOR*/
  //double **cmtx; /*MATRIX*/
  struct onode *ninit; /*node-initial*/

  ftxt=fopen("gradtest.txt","w");

  //for CG01.inp
  ftarget=0.99;
  gamma=0.002;
  dfact=0.0005;
  eps=0.000001;

  /*
  ftarget=0.99;
  gamma=0.001;
  dfact=0.001;
  eps=0.000001;
  */

  /*
  ftarget=0.8;
  gamma=0.001;
  dfact=0.05;
  eps=0.00001;
  */

  /*POLYGON01.INP*/
  /*
  ftarget=0.4;
  gamma=0.001;
  dxfact=0.05;
  eps=0.00001;
  */
  /*
  ftarget=0.9;
  gamma=0.002;
  dxfact=0.05;
  eps=0.00001;
  */

  nnode=org->nnode;
  nelem=org->nelem;
  m=nnode;

  ninit=(struct onode *)malloc(nnode*sizeof(struct onode));

  fgrad1=mallocdoublematrixxyz(nnode,3); /*FOR XYZ OF ALL NODES*/
  fgrad2=mallocdoublematrixxyz(nnode,3); /*FOR XYZ OF ALL NODES*/
  u1=mallocdoublematrixxyz(nnode,3);
  u2=mallocdoublematrixxyz(nnode,3);
  ua=mallocdoublematrixxyz(nnode,3);

  /*INITIAL NODE*/
  for(i=0; i<nnode; i++) *(ninit+i)=*(org->nodes+i);

  /*MAXIMUM SAFETY*/
  arclmautomatic(org);
  /*
  n=strcspn((wdraw.childs+1)->sctfile,".");
  strncpy(str,(wdraw.childs+1)->sctfile,n);
  str[n]='\0';
  strcpy((wdraw.childs+1)->otpfile,str);
  strcat((wdraw.childs+1)->otpfile,".rat");
  SetDlgItemText((wmenu.childs+2)->hwnd,ID_OUTPUTFILE,
				 (wdraw.childs+1)->otpfile);
  */
  fout=fopen((wdraw.childs+1)->otpfile,"r");
  if(fout==NULL) return;
  readsrcanrate(fout,&arc);
  fclose(fout);
  aelem=arc.nelem;
  f1=0.0;
  for(ii=0;ii<aelem;ii++)
  {
	for(jj=0;jj<4;jj++)
	{
	  if(f1<(arc.elems+ii)->srate[jj]) f1=(arc.elems+ii)->srate[jj];
	}
  }

  sprintf(str,"%d",0);
  SetDlgItemText((wmenu.childs+2)->hwnd,ID_LAPS,str);
  SendDlgItemMessage((wmenu.childs+2)->hwnd,ID_LAPS,WM_PAINT,0,0);
  sprintf(str,"%.5f",f1);
  SetDlgItemText((wmenu.childs+2)->hwnd,ID_SAFETY,str);
  SendDlgItemMessage((wmenu.childs+2)->hwnd,ID_SAFETY,WM_PAINT,0,0);

  sprintf(str,"Initial Max Safety = %9.5f",f1);
  fprintf(ftxt,"%s\n",str);
  MessageBox(NULL,str,"Conjugate Gradient",MB_OK);

  int mpt[]={101,102,115,122,119,123};//動かしたい点
  int nmpt= sizeof mpt / sizeof mpt[0];

  int codes[23];//本当はcodes[m]ってしたい
  for(i=0;i<m;i++)codes[i]=(org->nodes+i)->code;


  /*INITIAL GRADIENT*/
  fprintf(ftxt,"Initial Gradient\n");
  for(i=0;i<m;i++)
  {
	for(p=0;p<3;p++)//XZに限定
	{
		if(p!=1)/*&&(codes[i]==mpt[q])*/
		{
			for(q=0;q<nmpt;q++)
			{
			  if (codes[i]==mpt[q])
			  {
				(org->nodes+i)->d[p]=((ninit+i)->d[p])+dfact;

				for(j=0;j<m;j++)
				{
					d1=2.0;

					xi2=(org->nodes+i)->d[0];
					yi2=(org->nodes+i)->d[1];
					zi2=(org->nodes+i)->d[2];
					xj2=(org->nodes+j)->d[0];
					yj2=(org->nodes+j)->d[1];
					zj2=(org->nodes+j)->d[2];

					d2=sqrt((xi2-xj2)*(xi2-xj2)+(yi2-yj2)*(yi2-yj2)+(zi2-zj2)*(zi2-zj2));
					delta=d2/d1;

					if(delta<0.5 && j!=i) (org->nodes+j)->d[p]=((ninit+j)->d[p])+delta*dfact;
				}

				/*隣のnodeを探す*/
				/*for(j=0;j<nelem;j++)
				{
					int a;
					for(a=0;a<4;a++)
					{
						if(((org->elems+j)->enod[a]==point)
						{

						}
					}
				}
				*/

				/*MAXIMUM SAFETY*/
				arclmautomatic(org);
				fout=fopen((wdraw.childs+1)->otpfile,"r");
				readsrcanrate(fout,&arc);
				fclose(fout);
				aelem=arc.nelem;
				df=0.0;
				for(ii=0;ii<aelem;ii++)
				{
				  for(jj=0;jj<4;jj++)
				  {
					if(df<(arc.elems+ii)->srate[jj]) df=(arc.elems+ii)->srate[jj];
				  }
				}
				*(*(fgrad1+i)+p)=gamma*(f1-df)/dfact; /*NEGATIVE GRADIENT*/
				*(*(u1+i)+p)=*(*(fgrad1+i)+p);

				fprintf(ftxt,"Node %d Dim %d Gradient=%9.5f\n",(org->nodes+i)->code,p,(f1-df)/dfact);

				for(j=0;j<m;j++) (org->nodes+j)->d[p]=((ninit+j)->d[p]);

sprintf(str,"Step initial MovePt=%d \n p=%d MaxSafety=%.5f\n",mpt[q],p,df);
MessageBox(NULL,str,"Conjugate Gradient",MB_OK);
			  }
			}
		}
	}
  }
  fprintf(ftxt,"\n");

  /*
  fprintf(stderr,"INITIAL CONDITION\n");
  fprintf(stderr,"x   = %8.3f %8.3f\n",*(x1+0),*(x1+1));
  fprintf(stderr,"f(x)= %8.3f\n",f1);
  fprintf(stderr,"grad f(x)= %8.3f %8.3f\n",*(fgrad1+0),*(fgrad1+1));
  fprintf(stderr,"{u}      = %8.3f %8.3f\n",*(u1+0),*(u1+1));
  fprintf(fout,"STEP 0 f %8.3f {x1,x2} %8.3f %8.3f\n",f1,*(x1+0),*(x1+1));
  gets(non);
  */

  k=0;
  while(1)
  {
	k++;
	sprintf(str,"%d",k);
	SetDlgItemText((wmenu.childs+2)->hwnd,ID_LAPS,str);

	for(i=0;i<m;i++)
	{
		for(p=0;p<3;p++)//XZに限定
		{
		  if(p!=1)
		  {
			for(q=0;q<nmpt;q++)
			{
			  if (codes[i]==mpt[q])
			  {
				(org->nodes+i)->d[p]=((ninit+i)->d[p])+(*(*(u1+i)+p));

				for(j=0;j<m;j++)
				{
					double delta;
					double d1;
					double d2;

					d1=2.0;

					xi2=(org->nodes+i)->d[0];
					yi2=(org->nodes+i)->d[1];
					zi2=(org->nodes+i)->d[2];
					xj2=(org->nodes+j)->d[0];
					yj2=(org->nodes+j)->d[1];
					zj2=(org->nodes+j)->d[2];

					d2=sqrt((xi2-xj2)*(xi2-xj2)+(yi2-yj2)*(yi2-yj2)+(zi2-zj2)*(zi2-zj2));
					delta=d2/d1;

					if(delta<0.5 && i!=j) (org->nodes+j)->d[p]=((ninit+j)->d[p])+delta*(*(*(u1+i)+p));
				}
			  }
			}
		  }
		}
	}

	/*MAXIMUM SAFETY*/
	arclmautomatic(org);
	fout=fopen((wdraw.childs+1)->otpfile,"r");
	readsrcanrate(fout,&arc);
	fclose(fout);
	aelem=arc.nelem;
	fa=0.0;
	for(ii=0;ii<aelem;ii++)
	{
	  for(jj=0;jj<4;jj++)
	  {
		if(fa<(arc.elems+ii)->srate[jj]) fa=(arc.elems+ii)->srate[jj];
	  }
	}

//sprintf(str,"Step %d Max Safety(fa) = %9.5f\n",k,fa);
//MessageBox(NULL,str,"Conjugate Gradient",MB_OK);

	fprintf(ftxt,"Step %d Max Safety = %9.5f\n",k,fa);

	/*
	sprintf(str,"Max Safety = %.5f",fa);
	MessageBox(NULL,str,"Conjugate Gradient",MB_OK);
	*/

	/*INITIAL NODE*/
	for(i=0;i<nnode;i++) *(ninit+i)=*(org->nodes+i);

	/*GRADIENT*/
	fprintf(ftxt,"Step %d Gradient\n",k);
	for(i=0;i<m;i++)
	{
		for(p=0;p<3;p++)//XYに限定
		{
			if(p!=1)
			{
				for(q=0;q<nmpt;q++)
				{
				  if(codes[i]==mpt[q])
				  {
					(org->nodes+i)->d[p]=((ninit+i)->d[p])+dfact;

					for(j=0;j<m;j++)
					{
						double delta;
						double d1;
						double d2;

						d1=2.0;

						xi2=(org->nodes+i)->d[0];
						yi2=(org->nodes+i)->d[1];
						zi2=(org->nodes+i)->d[2];
						xj2=(org->nodes+j)->d[0];
						yj2=(org->nodes+j)->d[1];
						zj2=(org->nodes+j)->d[2];

						d2=sqrt((xi2-xj2)*(xi2-xj2)+(yi2-yj2)*(yi2-yj2)+(zi2-zj2)*(zi2-zj2));
						delta=d2/d1;

						if(delta<0.5 && i!=j) (org->nodes+j)->d[p]=((ninit+j)->d[p])+delta*dfact;
					}

					/*MAXIMUM SAFETY*/
					arclmautomatic(org);
					fout=fopen((wdraw.childs+1)->otpfile,"r");
					readsrcanrate(fout,&arc);
					fclose(fout);
					aelem=arc.nelem;
					df=0.0;
					for(ii=0;ii<aelem;ii++)
					{
						for(jj=0;jj<4;jj++)
						{
						  if(df<(arc.elems+ii)->srate[jj]) df=(arc.elems+ii)->srate[jj];
						}
					}
	//sprintf(str,"Step %d i=%d p=%d Max Safety(df)=%9.5f\n",k,i,p,df);
	//MessageBox(NULL,str,"Conjugate Gradient",MB_OK);

//					*(*(ua+i)+p)=-gamma*(fa-df)/dfact+(*(*(u1+i)+p)); /*WRONG:NEGATIVE GRADIENT*/
					*(*(ua+i)+p)=-gamma*(fa-df)/dfact+(*(*(fgrad1+i)+p)); /*NEGATIVE GRADIENT*/

	//sprintf(str,"Step %d i=%d p=%d ua=%9.5f\n",k,i,p,*(*(ua+i)+p));
	//MessageBox(NULL,str,"Conjugate Gradient",MB_OK);

					fprintf(ftxt,"Node %d Dim %d Gradient=%9.5f\n",(org->nodes+i)->code,p,(fa-df)/dfact);

					for(j=0;j<m;j++) (org->nodes+j)->d[p]=((ninit+j)->d[p]);
					//(org->nodes+i)->d[p]=(ninit+i)->d[p];
				  }
				}
			}
		}
	}
	fprintf(ftxt,"\n");

	/*
	fprintf(stderr,"ITERATION %d\n",k);
	fprintf(stderr," {u}= %8.3f %8.3f\n",*(u1+0),*(u1+1));
	fprintf(stderr,"A{u}= %8.3f %8.3f\n",*(ua+0),*(ua+1));
	*/

	c1=0.0;
	c2=0.0;
	for(i=0; i<m; i++)
	{
		for(p=0;p<3;p++)
		{
			c1+=(*(*(fgrad1+i)+p))*(*(*(fgrad1+i)+p));
			//c2+=((*(*(u1+i)+p))*(*(*(ua+i)+p)));
			c2+=(*(*(u1+i)+p))*(*(*(ua+i)+p));
			/*if (c2==0) c2=0.000001;//c2これで十分小さいか不明
			else c2=c2;*/

		}
	}
	alpha=c1/c2;
//sprintf(str,"Step %d c1=%.5f c2=%.5f alpha=%9.5f\n",k,c1,c2,alpha);
//MessageBox(NULL,str,"Conjugate Gradient",MB_OK);
	/*
	fprintf(stderr,"Alpha= %8.3f = %8.3f / %8.3f\n",alpha,c1,c2);
	*/

	vsize=0.0;
	for(i=0;i<m;i++)
	{
		for(p=0;p<3;p++)//XYに限定
		{
			if(p!=1)
			{
				for(q=0;q<nmpt;q++)
				{
				  if(codes[i]==mpt[q])
				  {
					(org->nodes+i)->d[p]=((ninit+i)->d[p])+alpha*(*(*(u1+i)+p));

					for(j=0;j<m;j++)
					{
						double delta;
						double d1;
						double d2;

						d1=2.0;

						xi2=(org->nodes+i)->d[0];
						yi2=(org->nodes+i)->d[1];
						zi2=(org->nodes+i)->d[2];
						xj2=(org->nodes+j)->d[0];
						yj2=(org->nodes+j)->d[1];
						zj2=(org->nodes+j)->d[2];

						d2=sqrt((xi2-xj2)*(xi2-xj2)+(yi2-yj2)*(yi2-yj2)+(zi2-zj2)*(zi2-zj2));
						delta=d2/d1;

						if(delta<0.5 && i!=j) (org->nodes+j)->d[p]=((ninit+j)->d[p])+delta*alpha*(*(*(u1+i)+p));
					}


					/*
					fprintf(stderr,"x = x + dx : %8.3f = %8.3f + %8.3f x %8.3f\n",*(x2+i),*(x1+i),alpha,*(u1+i));
					*/

					*(*(fgrad2+i)+p)=(*(*(fgrad1+i)+p))-alpha*(*(*(ua+i)+p));
					vsize+=(*(*(fgrad2+i)+p))*(*(*(fgrad2+i)+p));
				  }
				}
			}
		}
	}

	/*MAXIMUM SAFETY*/
	arclmautomatic(org);
	fout=fopen((wdraw.childs+1)->otpfile,"r");
	readsrcanrate(fout,&arc);
	fclose(fout);
	aelem=arc.nelem;
	f2=0.0;
	for(ii=0;ii<aelem;ii++)
	{
	  for(jj=0;jj<4;jj++)
	  {
		if(f2<(arc.elems+ii)->srate[jj]) f2=(arc.elems+ii)->srate[jj];
	  }
	}
//sprintf(str,"Step %d Max Safety(f2) = %9.5f\n",k,f2);
//MessageBox(NULL,str,"Conjugate Gradient",MB_OK);

	sprintf(str,"%d",k);
	SetDlgItemText((wmenu.childs+2)->hwnd,ID_LAPS,str);
	SendDlgItemMessage((wmenu.childs+2)->hwnd,ID_LAPS,WM_PAINT,0,0);
	sprintf(str,"%.5f",f2);
	SetDlgItemText((wmenu.childs+2)->hwnd,ID_SAFETY,str);
	SendDlgItemMessage((wmenu.childs+2)->hwnd,ID_SAFETY,WM_PAINT,0,0);

	sprintf(str,"Step=%d Max Safety=%.5f Vector Size=%.5f",k,f2,vsize);
	fprintf(ftxt,"%s\n\n",str);
	/*MessageBox(NULL,str,"Conjugate Gradient",MB_OK);*/

	/*
	fprintf(stderr,"x    = %8.3f %8.3f\n",*(x2+0),*(x2+1));
	fprintf(stderr,"f(x) = %8.3f\n",f2);
	fprintf(stderr,"grad f(x)   = %8.3f %8.3f\n",*(fgrad2+0),*(fgrad2+1));
	fprintf(stderr,"VECTOR SIZE = %9.5f\n",vsize);
	fprintf(fout,"STEP %d f %8.3f {x1,x2} %8.3f %8.3f\n",k,f2,*(x2+0),*(x2+1));
	gets(non);
	*/

	if(vsize<eps || f2<=ftarget)
	{
	  sprintf(str,"COMPLETED : TARGET f(x)=%.5f\n",f2);
	  MessageBox(NULL,str,"Conjugate Gradient",MB_OK);
	  return;
	}


	c3=0.0;
	for(i=0; i<m; i++)
	{
		for(p=0;p<3;p++)
		{
			if(p!=1)
			{
			  c3+=(*(*(fgrad2+i)+p))*(*(*(fgrad2+i)+p));
			}
		}
	}
	beta=c3/c1;
//sprintf(str,"Step %d Safety=%.5f Beta=%.5f\n",k,f2,beta);
//MessageBox(NULL,str,"Conjugate Gradient",MB_OK);

	/*繰り返しのため引き継ぎ*/
	for(i=0; i<m; i++)
	{
		/*INITIAL NODE*/
		*(ninit+i)=*(org->nodes+i);

		for(p=0;p<3;p++)
		{
			if(p!=1)
			{
			  *(*(u2+i)+p)=(*(*(fgrad2+i)+p))+beta*(*(*(u1+i)+p));
			  *(*(u1+i)+p)=*(*(u2+i)+p);
			  *(*(fgrad1+i)+p)=*(*(fgrad2+i)+p);
			}
		}
	}
  }

  //fclose(ftxt);

  return;
}/*conjugategradientcurve*/

double arclmautomatic(struct organ *org) /*CONJUGATE*/
/*OPTIMIZE ORGAN WITH CONJUGATE GRADIENT.*/
{
  char str[256],non[80];
  int i,j,k,n,flag;
  FILE *fout;

  /*TURN OFF MESSAGES*/
  globalmessageflag=0;
  globaldrawflag=0;

  /*SAFETY FLAG*/
  (wdraw.childs+1)->vparam.vflag.ev.srcancolor=1;

  arc =arci;
  arcx=arci;
  arcy=arci;
  /*free((wdraw.childs+1)->org.loads);*/

  for(i=0;i<org->nelem;i++)
  {
	for(j=0;j<2;j++)
	{
	  for(k=0;k<6;k++) (org->elems+i)->initial[j][k]=0.0;
	}
  } /*INITIAL CMQ UNAVAILABLE.*/

  /*EXTRACT ARCLM*/
  extractarclmfromorgan(org,&arc,&arcx,&arcy);

  (wmenu.childs+2)->vparam.vflag.mv.ftype=F_ARCLM;

  /*SAVE AS ARCLM*/
  saveasarclm((wdraw.childs+1)->inpfilez,&arc);
  saveasarclm((wdraw.childs+1)->inpfilex,&arcx);
  saveasarclm((wdraw.childs+1)->inpfiley,&arcy);

  getviewparam((wmenu.childs+2)->hwnd,
			   &((wdraw.childs+1)->vparam));
  (wdraw.childs+1)->vparam.vflag.ev.deformation=0;

  /*ANALYSIS*/
  arclm001(&arc ,ID_INPUTFILEZ,ID_OUTPUTFILEZ);
  arclm001(&arcx,ID_INPUTFILEX,ID_OUTPUTFILEX);
  arclm001(&arcy,ID_INPUTFILEY,ID_OUTPUTFILEY);

  /*SRCAN*/
  n=strcspn((wdraw.childs+1)->sctfile,".");
  strncpy(str,(wdraw.childs+1)->sctfile,n);
  str[n]='\0';

  i=srcan001(str);

  if(i==0) MessageBox(NULL,"Failed.","SRCAN001",MB_OK);

  /*RATE FILE NAME*/
  strcpy((wdraw.childs+1)->otpfile,str);
  strcat((wdraw.childs+1)->otpfile,".rat");
  SetDlgItemText((wmenu.childs+2)->hwnd,ID_OUTPUTFILE,
				 (wdraw.childs+1)->otpfile);
  fout=fopen((wdraw.childs+1)->otpfile,"r");
  if(fout==NULL) return 0.0;
  readsrcanrate(fout,&arc);
  fclose(fout);

  /*REDRAW MODEL*/
  clearwindow(*(wdraw.childs+1));
  drawarclmframe((wdraw.childs+1)->hdcC,
				 (wdraw.childs+1)->vparam,arc,0,ONSCREEN);
  SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);
  (wmenu.childs+2)->vparam.vflag.mv.ftype=F_ORGAN;
  SendMessage((wmenu.childs+2)->hwnd,WM_INITDIALOG,0,0);

  globalmessageflag=1;
  globaldrawflag=1;

  return 0.0;
}/*arclmautomatic*/

double arclmautomaticII(struct organ *org) /*CONJUGATE*/
/*OPTIMIZE ORGAN WITH CONJUGATE GRADIENT.*/
{
  char str[256],non[80];
  int i,j,k,n,flag;
  FILE *fout;
  int idinputs[3],idoutputs[3];
  struct arclmframe *arcs[3];

  /*TURN OFF MESSAGES*/
  globalmessageflag=0;
  globaldrawflag=0;

  /*SAFETY FLAG*/
  (wdraw.childs+1)->vparam.vflag.ev.srcancolor=1;

  arc =arci;
  arcx=arci;
  arcy=arci;
  /*free((wdraw.childs+1)->org.loads);*/

  for(i=0;i<org->nelem;i++)
  {
	for(j=0;j<2;j++)
	{
	  for(k=0;k<6;k++) (org->elems+i)->initial[j][k]=0.0;
	}
  } /*INITIAL CMQ UNAVAILABLE.*/

  /*EXTRACT ARCLM*/
  extractarclmfromorgan(org,&arc,&arcx,&arcy);

  (wmenu.childs+2)->vparam.vflag.mv.ftype=F_ARCLM;

  /*SAVE AS ARCLM*/
  saveasarclm((wdraw.childs+1)->inpfilez,&arc);
  saveasarclm((wdraw.childs+1)->inpfilex,&arcx);
  saveasarclm((wdraw.childs+1)->inpfiley,&arcy);

  getviewparam((wmenu.childs+2)->hwnd,
			   &((wdraw.childs+1)->vparam));
  (wdraw.childs+1)->vparam.vflag.ev.deformation=0;

  /*ANALYSIS*/
/*
  arclm001(&arc ,ID_INPUTFILEZ,ID_OUTPUTFILEZ);
  arclm001(&arcx,ID_INPUTFILEX,ID_OUTPUTFILEX);
  arclm001(&arcy,ID_INPUTFILEY,ID_OUTPUTFILEY);
*/
  arcs[0]=&arc;arcs[1]=&arcx;arcs[2]=&arcy;
  // arcs[0]=&arc;arcs[1]=&arc;arcs[2]=&arc;
  idinputs[0]=ID_INPUTFILEZ;idinputs[1]=ID_INPUTFILEX;idinputs[2]=ID_INPUTFILEY;
  idoutputs[0]=ID_OUTPUTFILEZ;idoutputs[1]=ID_OUTPUTFILEX;idoutputs[2]=ID_OUTPUTFILEY;
  arclm001_lxy(arcs,idinputs,idoutputs);

  /*BUCKLING ANALYSIS*/
  /*BUCKLIMG CONDENSATION*/
  n=strcspn((wdraw.childs+1)->inpfile,".");
  strncpy(str,(wdraw.childs+1)->inpfile,n);
  str[n]='\0';
  strcpy((wdraw.childs+1)->otpfile,str);
  strcat((wdraw.childs+1)->otpfile,".otbc");
  SetDlgItemText((wmenu.childs+2)->hwnd,ID_OUTPUTFILE,
                 (wdraw.childs+1)->otpfile);
//  (wdraw.childs+1)->vparam.vflag.ev.srcanrate =0;
  bclng002(&arc);
//  (wdraw.childs+1)->vparam.vflag.ev.srcanrate =0;

  /*SRCAN*/
  n=strcspn((wdraw.childs+1)->sctfile,".");
  strncpy(str,(wdraw.childs+1)->sctfile,n);
  str[n]='\0';

  i=srcan001(str);

  if(i==0) MessageBox(NULL,"Failed.","SRCAN001",MB_OK);

  /*RATE FILE NAME*/
  strcpy((wdraw.childs+1)->otpfile,str);
  strcat((wdraw.childs+1)->otpfile,".rat");
  SetDlgItemText((wmenu.childs+2)->hwnd,ID_OUTPUTFILE,
				 (wdraw.childs+1)->otpfile);
  fout=fopen((wdraw.childs+1)->otpfile,"r");
  if(fout==NULL) return 0.0;
  readsrcanrate(fout,&arc);
  fclose(fout);

  /*REDRAW MODEL*/
  clearwindow(*(wdraw.childs+1));
  drawarclmframe((wdraw.childs+1)->hdcC,
				 (wdraw.childs+1)->vparam,arc,0,ONSCREEN);
  SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);
  (wmenu.childs+2)->vparam.vflag.mv.ftype=F_ORGAN;
  SendMessage((wmenu.childs+2)->hwnd,WM_INITDIALOG,0,0);

  globalmessageflag=1;
  globaldrawflag=1;

  return 0.0;
}/*arclmautomaticII*/

/*-----------------------------------------------------------------*/
void bclngconjugategradient(struct organ *org)     /*Ujioka*/
/*OPTIMIZE ORGAN WITH CONJUGATE GRADIENT.*/
/*
Bclng Optimization (CG Method)
Maximize Buckling Eigen Value = Minimize Buckling Safety (1/Buckling Eigen Value)
*/
{
  FILE *fout,*ftxt;
  char non[10],str[256];
  int i,ii,j,jj,k,m,n;
  int nnode,nelem,aelem;
  double df,f1,f2,fa,ftarget,c1,c2,c3,alpha,beta,gamma,vsize,eps;
  double *x1,*x2,*xa,*xx,*dx,dfact; /*COORDINATES*/
  double *u1,*u2,*ua,*fgrad1,*fgrad2; /*GRADIENT VECTOR*/
  double **cmtx; /*MATRIX*/
  struct onode *ninit;            /*節点座標の初期値 d[0],d[1],d[2]*/
  double beig;                    /*Buckling Eigen Value*/

  ftxt=fopen("gradtest.txt","w");
  ftarget=0.3;
  gamma=0.001;
  dfact=0.01;            /*座標の微小変化　Δz*/
  eps=0.00001;

  nnode=org->nnode;
  nelem=org->nelem;
  m=nnode;

  ninit=(struct onode *)malloc(nnode*sizeof(struct onode));
  fgrad1=mallocdoublevector(nnode); /*FOR Z OF ALL NODES*/
  fgrad2=mallocdoublevector(nnode); /*FOR Z OF ALL NODES*/
  u1=mallocdoublevector(nnode);
  u2=mallocdoublevector(nnode);
  ua=mallocdoublevector(nnode);

  /*INITIAL NODE*/
  for(i=0; i<nnode; i++) *(ninit+i)=*(org->nodes+i);

  /*Buckling Eigen Value*/
  beig=bclngautomatic(org);
  f1=1/beig;
  sprintf(str,"INITIAL BUCKLING EIGEN VALUE =%8.3f",beig);
  MessageBox(NULL,str,"BCLNG001 RESULT",MB_OK);

  sprintf(str,"%d",0);
  SetDlgItemText((wmenu.childs+2)->hwnd,ID_LAPS,str);
  SendDlgItemMessage((wmenu.childs+2)->hwnd,ID_LAPS,WM_PAINT,0,0);
  sprintf(str,"%.5f",f1);
  SetDlgItemText((wmenu.childs+2)->hwnd,ID_SAFETY,str);
  SendDlgItemMessage((wmenu.childs+2)->hwnd,ID_SAFETY,WM_PAINT,0,0);

  /*INITIAL GRADIENT*/
  fprintf(ftxt,"Initial Gradient\n");
  for(i=0;i<m;i++)
  {
	for(j=0;j<m;j++)       /*z座標を微小量(dfact)変化*/
	{
	  if(j==i)
	  {
		(org->nodes+j)->d[2]=((ninit+j)->d[2])+dfact; /*UPDATE ONLY Z OF NODE I*/
	  }
	  else
	  {
		(org->nodes+j)->d[2]=((ninit+j)->d[2]); /*RESET Z OF OTHER NODES*/
	  }
	}

    /*Buckling Eigen Value*/
    beig=bclngautomatic(org);
    /*
    sprintf(str,"COMPLETED : Node %d Buckling Eigen Value=%.5f\n",(org->nodes+i)->code,beig);
    MessageBox(NULL,str,"Conjugate Gradient",MB_OK);
    */
    df=0.0;
    df=1/beig;             /*Buckling Safety*/

	*(fgrad1+i)=gamma*(f1-df)/dfact; /*NEGATIVE GRADIENT*/
	*(u1+i)=*(fgrad1+i);

	fprintf(ftxt,"Node %d Gradient=%9.5f\n",(org->nodes+i)->code,(f1-df)/dfact);
  }
  fprintf(ftxt,"\n");

  k=0;
  while(1)
  {
	k++;
	sprintf(str,"%d",k);
	SetDlgItemText((wmenu.childs+2)->hwnd,ID_LAPS,str);

	for(i=0;i<m;i++)
	{
	  (org->nodes+i)->d[2]=((ninit+i)->d[2])+(*(u1+i)); /*ONLY Z*/
	}

    /*Buckling Eigen Value*/
    beig=bclngautomatic(org);
    fa=0.0;
    fa=1/beig;           					   /*Buckling Safety*/


    fprintf(ftxt,"Step %d Buckling Safety = %9.5f\n",k,fa);


	/*INITIAL NODE*/
	for(i=0;i<nnode;i++) *(ninit+i)=*(org->nodes+i);

	/*GRADIENT*/
	fprintf(ftxt,"Step %d Gradient\n",k);
	for(i=0;i<m;i++)
	{
	  for(j=0;j<m;j++)
	  {
		if(j==i)
		{
		  (org->nodes+j)->d[2]=((ninit+j)->d[2])+dfact; /*ONLY Z*/
		}
		else
		{
		  (org->nodes+j)->d[2]=((ninit+j)->d[2]); /*ONLY Z*/
		}
	  }

      /*Buckling Eigen Value*/
  	  beig=bclngautomatic(org);
  	  df=0.0;
      df=1/beig;

//	  *(ua+i)=-gamma*(fa-df)/dfact+(*(u1+i)); /*WRONG:NEGATIVE GRADIENT*/
      *(ua+i)=-gamma*(fa-df)/dfact+(*(fgrad1+i)); /*NEGATIVE GRADIENT*/

	  fprintf(ftxt,"Node %d Gradient=%9.5f\n",(org->nodes+i)->code,(fa-df)/dfact);
	}
	fprintf(ftxt,"\n");

	/*
	fprintf(stderr,"ITERATION %d\n",k);
	fprintf(stderr," {u}= %8.3f %8.3f\n",*(u1+0),*(u1+1));
	fprintf(stderr,"A{u}= %8.3f %8.3f\n",*(ua+0),*(ua+1));
	*/

	c1=0.0;
	c2=0.0;
	for(i=0; i<m; i++)
	{
	  c1+=(*(fgrad1+i))*(*(fgrad1+i));
	  c2+=(*(u1+i))*(*(ua+i));
	}
	alpha=c1/c2;

	/*
	fprintf(stderr,"Alpha= %8.3f = %8.3f / %8.3f\n",alpha,c1,c2);
	*/

	vsize=0.0;
	for(i=0;i<m;i++)
	{
	  (org->nodes+i)->d[2]=((ninit+i)->d[2])+alpha*(*(u1+i)); /*ONLY Z*/

	  /*
	  fprintf(stderr,"x = x + dx : %8.3f = %8.3f + %8.3f x %8.3f\n",*(x2+i),*(x1+i),alpha,*(u1+i));
	  */

	  *(fgrad2+i)=(*(fgrad1+i))-alpha*(*(ua+i));
	  vsize+=(*(fgrad2+i))*(*(fgrad2+i));
	}

    /*Buckling Eigen Value*/
    beig=bclngautomatic(org);
    f2=0.0;
    f2=1/beig;          /*Buckling Safety*/

	sprintf(str,"%d",k);
	SetDlgItemText((wmenu.childs+2)->hwnd,ID_LAPS,str);
	SendDlgItemMessage((wmenu.childs+2)->hwnd,ID_LAPS,WM_PAINT,0,0);
	sprintf(str,"%.5f",f2);
	SetDlgItemText((wmenu.childs+2)->hwnd,ID_SAFETY,str);
	SendDlgItemMessage((wmenu.childs+2)->hwnd,ID_SAFETY,WM_PAINT,0,0);

	sprintf(str,"Step=%d Buckling Safety=%.5f Vector Size=%.5f",k,f2,vsize);
	fprintf(ftxt,"%s\n\n",str);
	/*MessageBox(NULL,str,"Conjugate Gradient",MB_OK);*/

	if(vsize<eps || f2<=ftarget)
	{
     /*
	  sprintf(str,"COMPLETED : TARGET f(x)=%.5f\n",f2);
	  MessageBox(NULL,str,"Conjugate Gradient",MB_OK);
     */
      beig=bclngautomatic(org);
      sprintf(str,"COMPLETED : Buckling Eigen Value=%.5f\n",beig);
      fprintf(ftxt,"COMPLETED : Buckling Eigen Value=%.5f\n",beig);
	  MessageBox(NULL,str,"Conjugate Gradient",MB_OK);

     /*save as "optimized.inp"*/
      fout=fopen("optimized.inp","w");
      saveorganization(fout,&((wdraw.childs+1)->org),
                            &((wdraw.childs+1)->vparam));
      fclose(fout);
	  return;
	}

	c3=0.0;
	for(i=0; i<m; i++)
	{
	  c3+=(*(fgrad2+i))*(*(fgrad2+i));
	}
	beta=c3/c1;

	for(i=0; i<m; i++)
	{
	  /*UPDATE INITIAL NODE*/
	  *(ninit+i)=*(org->nodes+i);

	  /*GRADIENTS*/
	  *(u2+i)=(*(fgrad2+i))+beta*(*(u1+i));
	  *(u1+i)=*(u2+i);
	  *(fgrad1+i)=*(fgrad2+i);
	}
  }

  fclose(ftxt);

  return;
}/*conjugategradient*/

double bclngautomatic(struct organ *org)
/*CONJUGATE*/
/*OPTIMIZE ORGAN WITH CONJUGATE GRADIENT.*/
{
  char str[256],non[80];
  int i,j,k,n,flag;
  FILE *fout;
  double beig = 0.0;
  /*TURN OFF MESSAGES*/
  globalmessageflag=0;
  globaldrawflag=0;

  /*SAFETY FLAG*/
  (wdraw.childs+1)->vparam.vflag.ev.srcancolor=1;

  arc =arci;
  arcx=arci;
  arcy=arci;
  /*free((wdraw.childs+1)->org.loads);*/

  for(i=0;i<org->nelem;i++)
  {
	for(j=0;j<2;j++)
	{
	  for(k=0;k<6;k++) (org->elems+i)->initial[j][k]=0.0;
	}
  } /*INITIAL CMQ UNAVAILABLE.*/

  /*EXTRACT ARCLM*/
  extractarclmfromorgan(org,&arc,&arcx,&arcy);

  (wmenu.childs+2)->vparam.vflag.mv.ftype=F_ARCLM;

  /*SAVE AS ARCLM*/
  saveasarclm((wdraw.childs+1)->inpfilez,&arc);

  getviewparam((wmenu.childs+2)->hwnd,
			   &((wdraw.childs+1)->vparam));
  (wdraw.childs+1)->vparam.vflag.ev.deformation=0;

  /*ANALYSIS*/
  arclm001(&arc ,ID_INPUTFILEZ,ID_OUTPUTFILEZ);
  beig=beigen(&arc);
  //sprintf(str,"BUCKLING EIGEN VALUE =%8.3f",beig);
  //MessageBox(NULL,str,"BCLNG001 RESULT",MB_OK);


  /*REDRAW MODEL*/
  clearwindow(*(wdraw.childs+1));
  drawarclmframe((wdraw.childs+1)->hdcC,
				 (wdraw.childs+1)->vparam,arc,0,ONSCREEN);
  SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);
  (wmenu.childs+2)->vparam.vflag.mv.ftype=F_ORGAN;
  SendMessage((wmenu.childs+2)->hwnd,WM_INITDIALOG,0,0);

  globalmessageflag=1;
  globaldrawflag=1;

  return beig;
}/*bclngautomatic*/

double beigen(struct arclmframe *af)
/*ELASTIC BUCKLING FOR ARCLM FRAME.*/
/*AFTER "ARCLM001".*/
{
  DWORD memory0,memory1;

  FILE /**fin,*/*fout;                               /*FILE 8 BYTES*/
  /*FILE *felem,*fdisp,*freact;*/
  /*char dir[]=DIRECTORY;*/                        /*DATA DIRECTORY*/
  char /*s[80],*/string[400];
  int i,j,ii,jj;
  int nnode,nelem,nsect/*,nreact*/;
  long int /*fsize,*/loff,msize;
  /*long int time;*/
  struct gcomponent ginit={0,0,0.0,NULL};
  struct gcomponent *kmtx,*gmtx;                    /*GLOBAL MATRIX*/
  double **gvct;                                    /*GLOBAL VECTOR*/
  double **drccos,**tmatrix,**estiff,*estress;
  /*double determinant,data;*/
  clock_t t0/*,t1,t2*/;

  /*struct osect *sects;*/
  struct onode /**nodes,*/*ninit;
  struct owire elem;
  struct owire *elems;
  struct oconf *confs;

  long int neig;
  double eps=1.0E-16,*eigen;

  memory0=availablephysicalmemory("INITIAL:");   /*MEMORY AVAILABLE*/

  fout=fgetstofopen("\0","w",ID_OUTPUTFILE);          /*OUTPUT FILE*/
  globalfile=fout;

  t0=clock();                                        /*CLOCK BEGIN.*/

  nnode=af->nnode;                                 /*INPUT INITIAL.*/
  nelem=af->nelem;
  nsect=af->nsect;
  /*nreact=af->nreact;*/
  sprintf(string,"NODES=%d ELEMS=%d SECTS=%d",nnode,nelem,nsect);
  errormessage(string);
  fprintf(fout,"%s\n",string);

  msize=6*nnode;                           /*SIZE OF GLOBAL MATRIX.*/
  neig=NEIGEN;                             /*NUMBERS OF EIGEN VALUES TO BE OBTAINED.*/

  //ujioka
  //if(SOLVER==1)
  //{
  //  eps=BISECEPS;
    /* neig=1; */
  //}

  kmtx=(struct gcomponent *)          /*DIAGONALS OF GLOBAL MATRIX.*/
        malloc(msize*sizeof(struct gcomponent));
  gmtx=(struct gcomponent *)
        malloc(msize*sizeof(struct gcomponent));
  if(kmtx==NULL || gmtx==NULL) return 0.0;
  for(i=0;i<msize;i++)
  {
    (kmtx+i)->down=NULL;            /*GLOBAL MATRIX INITIALIZATION.*/
    (gmtx+i)->down=NULL;
  }

  gvct=(double **)malloc(neig*sizeof(double *));   /*GLOBAL VECTORS*/
  eigen=(double *)malloc(neig*sizeof(double));       /*EIGEN VALUES*/
  if(gvct==NULL || eigen==NULL) return 0.0;
  for(i=0;i<neig;i++)
  {
    *(gvct+i)=(double *)malloc(msize*sizeof(double));
    for(j=0;j<msize;j++) *(*(gvct+i)+j)=0.0;
  }

  /*fdisp=af->fdisp;*/                 /*DISPLACEMENT:6 DIRECTIONS.*/
  /*felem=af->felem;*/              /*CODE,12 BOUNDARIES,12 STRESS.*/
  /*freact=af->freact;*/                           /*REACTION FILE.*/

  /*sects=af->sects;*/
  /*nodes=af->nodes;*/
  ninit=af->ninit;                                /*struct onode*/
  elems=af->elems;
  confs=af->confs;

  GetAsyncKeyState(VK_LBUTTON);                  /*CLEAR KEY LEFT.*/
  GetAsyncKeyState(VK_RBUTTON);                  /*CLEAR KEY RIGHT.*/

  errormessage("BCLNG001:BUCKLING LINEAR.");
  availablephysicalmemory("REMAIN:");            /*MEMORY AVAILABLE*/

  for(i=1;i<=msize;i++)             /*GLOBAL MATRIX INITIALIZATION.*/
  {
	ginit.m=(unsigned int)i;
	/*ginit.n=(unsigned int)i;*/
	*(kmtx+(i-1))=ginit;
    *(gmtx+(i-1))=ginit;
  }
  /*comps=msize;*/ /*INITIAL COMPONENTS=DIAGONALS.*/

  laptime("ASSEMBLING GLOBAL MATRIX.",t0);

  estress=(double *)malloc(12*sizeof(double));

  for(i=1;i<=nelem;i++)                 /*ASSEMBLAGE GLOBAL MATRIX.*/
  {
    /*elem=*(elems+i-1);*/                     /*READ ELEMENT DATA.*/
    inputelem(elems,af->melem,i-1,&elem);

    for(ii=0;ii<=1;ii++) /*INITIALIZE COORDINATION.*/
    {
      loff=elem.node[ii]->loff;
      for(jj=0;jj<3;jj++)
      {
        elem.node[ii]->d[jj]=(ninit+loff)->d[jj];
      }
    }

	drccos=directioncosine(elem.node[0]->d[0],
                           elem.node[0]->d[1],
                           elem.node[0]->d[2],
                           elem.node[1]->d[0],
                           elem.node[1]->d[1],
                           elem.node[1]->d[2],
                           elem.cangle);                 /*[DRCCOS]*/

    tmatrix=transmatrix(drccos);           /*TRANSFORMATION MATRIX.*/
    estiff=assememtx(elem);            /*ELASTIC MATRIX OF ELEMENT.*/
    estiff=modifyhinge(elem,estiff);               /*MODIFY MATRIX.*/
    estiff=transformation(estiff,tmatrix);         /*[K]=[Tt][k][T]*/

//for(ii=0;ii<12;ii++)
//{
//  sprintf(string,"\0");
//  for(jj=0;jj<=ii;jj++)
//  {
//	sprintf(s," %12.5E",*(*(estiff+ii)+jj));
//	strcat(string,s);
//  }
//  errormessage(string);
//}

    assemgstiffness(kmtx,estiff,&elem);       /*ASSEMBLAGE ELASTIC.*/

    for(ii=0;ii<=11;ii++) free(*(estiff+ii));
    free(estiff);

    for(ii=0;ii<=1;ii++)                                /*STRESSES.*/
    {
      for(jj=0;jj<6;jj++) *(estress+6*ii+jj)=elem.stress[ii][jj];
    }

//sprintf(string,"\0");
//for(ii=0;ii<12;ii++)
//{
//  sprintf(s," %12.5E",*(estress+ii));
//  strcat(string,s);
//}
//errormessage(string);

    estiff=assemgmtx(elem,estress);
    estiff=modifyhinge(elem,estiff);               /*MODIFY MATRIX.*/
    estiff=transformation(estiff,tmatrix);         /*[K]=[Tt][k][T]*/
    for(ii=0;ii<12;ii++)
    {
      for(jj=0;jj<12;jj++) *(*(estiff+ii)+jj)=-*(*(estiff+ii)+jj);
    }

    assemgstiffness(gmtx,estiff,&elem);     /*ASSEMBLAGE GEOMETRIC.*/

    for(ii=0;ii<=11;ii++) free(*(estiff+ii));
	free(estiff);

    for(ii=0;ii<=2;ii++) free(*(drccos+ii));
    free(drccos);
    for(ii=0;ii<=11;ii++) free(*(tmatrix+ii));
    free(tmatrix);
  }
  free(estress);
  laptime("GLOBAL MATRIX ASSEMBLED.",t0);

  /*currentvalue("GLOBAL MATRIX:[K]",msize,neig,kmtx,NULL,NULL,NULL);*/
  /*currentvalue("GLOBAL MATRIX:[G]",msize,neig,gmtx,NULL,NULL,NULL);*/

  //SOLVER  ujioka
//  deigabgeneral(gmtx,kmtx,confs,msize,neig,neig,eps,eigen,gvct);
  bisecsylvester(gmtx,kmtx,confs,msize,neig,neig,eps,eigen,gvct);

  /*
  if(MessageBox(NULL,"DEIGABGENERAL","SOLVER",MB_OKCANCEL)==IDOK)
		deigabgeneral(gmtx,kmtx,confs,msize,neig,neig,eps,eigen,gvct);
  else if(MessageBox(NULL,"BISECSYLVESTER","SOLVER",MB_OKCANCEL)==IDOK)
  {
        eps=BISECEPS;
		bisecsylvester(gmtx,kmtx,confs,msize,neig,neig,eps,eigen,gvct);
  }
  */
  /*
  if(SOLVER==0)
  {
    deigabgeneral(gmtx,kmtx,confs,msize,neig,neig,eps,eigen,gvct);
  }
  else
  {
	bisecsylvester(gmtx,kmtx,confs,msize,neig,neig,eps,eigen,gvct);
  }
  */

  laptime("EIGEN COMPLETED.",t0);

  for(i=0;i<neig;i++)
  {
    *(eigen+i)=1.0/(*(eigen+i));
    sprintf(string,"EIGEN VALUE %ld=%.5E",(i+1),*(eigen+i));
    fprintf(fout,"%s\n",string);
    errormessage(string);
    outputmode(*(gvct+i),fout,nnode,ninit);
  }

  updatemode(af,*(gvct+0)); /*FORMATION UPDATE.*/

  //ujioka
  double beig = 0.0;
  beig = *(eigen+0);
  /*
  sprintf(string,"BUCKLING EIGEN VALUE =%8.3f",beig);
  MessageBox(NULL,string,"BCLNG001 RESULT",MB_OK);
  */

  af->nlaps=neig;
  af->eigenval=eigen;
  af->eigenvec=gvct;

  gfree(kmtx,nnode); /*FREE GLOBAL MATRIX.*/
  gfree(gmtx,nnode); /*FREE GLOBAL MATRIX.*/

  errormessage(" ");
  errormessage("COMPLETED.");
  fprintf(fout,"COMPLETED.\n");

  fclose(fout);

  memory1=availablephysicalmemory("REMAIN:");
  sprintf(string,"CONSUMPTION:%ld[BYTES]",(memory0-memory1));
  errormessage(string);
  return beig;
}/*beigen*/

/*-----------------------------------------------------------------*/
int arclm201automatic(struct organ *org,int iteration)
				        		  /*iterationは繰り返し計算用（未完成）*/
{
  char str[256],non[80];
  int i,j,k,n,ii,jj,code;
  int nelem;
  struct oelem *einit;
  FILE *fin,*fout;

  /*TURN OFF MESSAGES*/
  globalmessageflag=0;
  globaldrawflag=0;

  arc  = arci;
  arcx = arci;
  arcy = arci;
  /*free((wdraw.childs+1)->org.loads);*/

  for(i=0;i<org->nelem;i++)
  {
	for(j=0;j<2;j++)
	{
	  for(k=0;k<6;k++) (org->elems+i)->initial[j][k]=0.0;
	}
  } /*INITIAL CMQ UNAVAILABLE.*/

  einit=(struct oelem *)malloc(org->nelem*sizeof(struct oelem));
  for(k=0; k<org->nelem; k++)  *(einit+k)=*(org->elems+k);

  /*EXTRACT ARCLM*/
  extractarclmfromorgan(org,&arc,&arcx,&arcy);

  (wmenu.childs+2)->vparam.vflag.mv.ftype=F_ARCLM;

  /*SAVE AS ARCLM*/
  saveasarclm((wdraw.childs+1)->inpfilez,&arc);

  getviewparam((wmenu.childs+2)->hwnd,
			   &((wdraw.childs+1)->vparam));
  (wdraw.childs+1)->vparam.vflag.ev.deformation=0;

  /*ANALYSIS*/
  arclm201(&arc ,ID_INPUTFILEZ);

//  arc.nlaps=1;


  clearwindow(*(wdraw.childs+1));
  drawarclmframe((wdraw.childs+1)->hdcC,
                 (wdraw.childs+1)->vparam,arc,0,ONSCREEN);
  SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);

//SetDlgItemText(hdwnd,IDV_MODENUM,"1");

  (wdraw.childs+1)->lstatus=ROTATE;
  (wmenu.childs+2)->vparam.vflag.mv.ftype=F_ORGAN;

  /*CHANGE NODES OF ORGAN FRAME*/
  for(i=0;i<org->nnode;i++)
  {
	  (org->nodes+i)->d[0]=(arc.nodes+i)->d[0];/* X */
	  (org->nodes+i)->d[1]=(arc.nodes+i)->d[1];/* Y */
      (org->nodes+i)->d[2]=(arc.nodes+i)->d[2];/* Z */
  }

  globalmessageflag=1;

return 0;
}/*arclm201automatic*/

void createbeziersurfacetest(struct organ *org)
{
  FILE *ftext;
  char **data,str[256]="\0",s[256];
  double ddata;
  int ndata;
  long int nnode,ncode;

  double *x,*y,*z;   /*Controle points*/
  int m,n;           /*Degree*/
  int ncontrole;

  double *u,*v;      /*U-V Coordinates*/
  struct onode node;

  int i,j;

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
  if(org->nnode!=nnode) return;

  /*Controle points*/
  ncontrole=(m+1)*(n+1);
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

  /*U-V COORDINATE OF EACH NODE*/
  u=(double *)malloc((nnode)*sizeof(double));
  v=(double *)malloc((nnode)*sizeof(double));
  for(i=0;i<nnode;i++)
  {
	data=fgetsbrk(ftext,&ndata);
	if(ndata!=5) return;

	ncode=strtol(*(data+1),NULL,10);
	if(ncode!=(org->nodes+i)->code) return;

    ddata=strtod(*(data+3),NULL);
    *(u+i)=ddata;

    ddata=strtod(*(data+4),NULL);
    *(v+i)=ddata;

    for(;ndata>0;ndata--) free(*(data+ndata-1));
    free(data);
  }
  fclose(ftext);

  /*BEZIER SURFACE*/
  for(i=0;i<nnode;i++)
  {
    node=*(org->nodes+i);
//    beziersurface(m,n,x,y,z,*(u+i),*(v+i),&node);
    beziersurfaceII(m,n,x,y,z,*(u+i),*(v+i),&node);

    (org->nodes+i)->d[0]=node.d[0];
    (org->nodes+i)->d[1]=node.d[1];
    (org->nodes+i)->d[2]=node.d[2];
  }

  /*REDRAW*/
  clearwindow(*(wdraw.childs+1));
  draworganization((wdraw.childs+1)->hdcC,
                   (wdraw.childs+1)->vparam,
                   (wdraw.childs+1)->org,ONSCREEN);

  /*DRAW CONTROLE POINTS*/
  for(i=0;i<ncontrole;i++)
  {
    node.d[0]=x[i];
    node.d[1]=y[i];
    node.d[2]=z[i];
    drawcontrolepoint((wdraw.childs+1)->hdcC,
                      (wdraw.childs+1)->vparam,
                      node);
  }

  SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);

return;
}/*createbeziersurfacetest*/

int beziersurface(int m,int n,double *x,double *y,double *z,
                  double u,double v,struct onode *node)
{
  int i,j,k;
  double p[3];
  char str[256];

  if(u<0 || u>1 || v<0 || v>1) return 0;    /*0<=u,v<=1*/

  p[0]=0; p[1]=0; p[2]=0;  k=0;
  for(j=0;j<n+1;j++)
  {
    for(i=0;i<m+1;i++)
    {
//      sprintf(str,"u[%d,%d]=%.5f",i,j,u);
//      errormessage(str);
      p[0]+=(bernstein(m,i,u))*(bernstein(n,j,v))*(*(x+k));
      p[1]+=(bernstein(m,i,u))*(bernstein(n,j,v))*(*(y+k));
      p[2]+=(bernstein(m,i,u))*(bernstein(n,j,v))*(*(z+k));
      k++;
    }
  }

  node->d[GX]=p[0];
  node->d[GY]=p[1];
  node->d[GZ]=p[2];

return 1;
}/*beziersurface*/

int beziersurfaceII(int m,int n,double *x,double *y,double *z,
                    double u,double v,struct onode *node)
/*FOR X,Y SYNMETRICAL BEZIER SURFACE*/
/*u<0:X SYNMETRICAL MOVE , v<0:Y SYNMETRICAL MOVE*/
{
  int i,j,k;
  double p[3];
  char str[256];
  int signu,signv;     /*FOR SYNMETRICAL MOVEMENT*/

//  if(u<0 || u>1 || v<0 || v>1) return 0;    /*0<=u,v<=1*/

  if(u>1 || v>1 || u<-1 || v<-1) return 0;    /*0<=u,v<=1*/

  signu=1; signv=1;
  if(u<0) signu=-1;
  if(v<0) signv=-1;


  p[0]=0; p[1]=0; p[2]=0;  k=0;
  for(j=0;j<n+1;j++)
  {
    for(i=0;i<m+1;i++)
    {
//      sprintf(str,"u[%d,%d]=%.5f",i,j,u);
//      errormessage(str);
      p[0]+=(bernstein(m,i,abs(u)))*(bernstein(n,j,abs(v)))*(*(x+k));
      p[1]+=(bernstein(m,i,abs(u)))*(bernstein(n,j,abs(v)))*(*(y+k));
      p[2]+=(bernstein(m,i,abs(u)))*(bernstein(n,j,abs(v)))*(*(z+k));
      k++;
    }
  }

  node->d[GX]=p[0]*signu;  /*X SYNMETRICAL MOVE*/
  node->d[GY]=p[1]*signv;  /*Y SYNMETRICAL MOVE*/
  node->d[GZ]=p[2];

return 1;
}/*beziersurfaceII*/

double bernstein(int n,int i,double t)
{
  double f;
  char str[256];

  //  sprintf(str,"n=%d,i=%d,t=%.5f",n,i,t);
  //  errormessage(str);

  if(t==0 && i==0)          f=(fact(n)/fact(i)/fact(n-i))*pow(1-t,n-i);
  else if(t==1 && n==i)     f=(fact(n)/fact(i)/fact(n-i))*pow(t,i);
  else                      f=(fact(n)/fact(i)/fact(n-i))*pow(t,i)*pow(1-t,n-i);
  return f;
}

int fact(int n)
{
  int i;
  int f;

  f=1;
  for(i=1;i<=n;i++)
  {
      f*=i;
  }

return f;
}

void drawcontrolepoint(HDC hdc,struct viewparam vp,struct onode gn)
/*DRAW NODE GLOBAL ON SCREEN BY AXONOMETRICS,PERSPECTIVE.*/
/*TEXTCOLOR DEFINITION MUST BE ALREADY DONE FOR FAST DRAWING.*/
{
  char str[20];
  int Ox,Oy;                                   /*COORDINATION ARROW*/
  SIZE size;

  if(!nodeontoscreen(gn,&Ox,&Oy,vp)) return;           /*PROJECTION*/

//  if(vp.vflag.nv.code)
  {
    /*
    sprintf(str,"test");

	TextOut(hdc,Ox,Oy,str,strlen(str));
    */
    HPEN hpen,hpenrange,ppen;
    HBRUSH hbrush,hbrushrange,pbrush;
    LOGBRUSH lb;

    hpen        = (HPEN)GetCurrentObject(hdc,OBJ_PEN);
    hbrush      = (HBRUSH)GetCurrentObject(hdc,OBJ_BRUSH);
    hbrushrange = CreateBrushIndirect(&lb);
    hpenrange   = CreatePen(PS_SOLID,1,RGB(150,150,150));
    ppen        = (HPEN)SelectObject(hdc,hpenrange);
    pbrush      = (HBRUSH)SelectObject(hdc,hbrushrange);

//    Ellipse(hdc,(Ox-5),(Oy-5),(Ox+5),(Oy+5));  /*FILLED CIRCLE*/
    Ellipse(hdc,(Ox-25),(Oy-25),(Ox+25),(Oy+25));  /*FILLED CIRCLE*/

    SelectObject(hdc,hbrush);
    SelectObject(hdc,hpen);
    DeleteObject(hpenrange);
    DeleteObject(hbrushrange);

  }
  return;
}/*drawglobalnode*/

void createvierendeelarch(struct organ *org)
{
  FILE *ftext;
  char **data,str[256]="\0",s[256];
  double ddata;
  int ndata;
  long int nnode1,nnode2,ncode,ncode2;

  double *x,*y,*z;   /*Controle points*/
  int m,n;           /*Degree*/
  int ncontrole;

  double u,v;      /*U-V Coordinates*/
//  double *u,*v;      /*U-V Coordinates*/

  double xoffset,yoffset,zoffset;   /*XYZ OFFSET*/

  struct onode node,node2;

  int i,j,jj;

  /*OPEN FILE*/
  ftext=fopen("vierendeel.txt","r");   /*bezier-surface data*/
  if(ftext==NULL)
  {
    errormessage("ACCESS IMPOSSIBLE.");
    return;
  }
  fseek(ftext,0L,SEEK_SET);

  data=fgetsbrk(ftext,&ndata);
  m=strtol(*(data+0),NULL,10);
  n=strtol(*(data+1),NULL,10);
  nnode1=strtol(*(data+2),NULL,10);   /*bezier points*/
  nnode2=strtol(*(data+3),NULL,10);   /*offset points*/
  if(org->nnode!=nnode1+nnode2) return;

  /*Controle points*/
  ncontrole=(m+1)*(n+1);
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

  /*U-V COORDINATE OF EACH NODE*/
//  u=(double *)malloc((nnode)*sizeof(double));
//  v=(double *)malloc((nnode)*sizeof(double));
  for(i=0;i<nnode1;i++)
  {
	data=fgetsbrk(ftext,&ndata);
	if(ndata!=5) break;

	ncode=strtol(*(data+1),NULL,10);
    j=0;
	while(1)
    {
      if(ncode==(org->nodes+j)->code) break;
      j++;
     if(j>(org->nodes+(org->nnode)-1)->loff) return;
    }

    if(ndata==5)    /*bezier points*/
    {
      ddata=strtod(*(data+3),NULL);
      u=ddata;

      ddata=strtod(*(data+4),NULL);
      v=ddata;

      for(;ndata>0;ndata--) free(*(data+ndata-1));
      free(data);

      node=*(org->nodes+j);
      beziersurface(m,n,x,y,z,u,v,&node);

      (org->nodes+j)->d[0]=node.d[0];
      (org->nodes+j)->d[1]=node.d[1];
      (org->nodes+j)->d[2]=node.d[2];
    }
  }

  for(i=0;i<nnode2;i++)
  {
	data=fgetsbrk(ftext,&ndata);
	if(ndata!=7) break;

	ncode=strtol(*(data+1),NULL,10);
    j=0;
	while(1)
    {
      if(ncode==(org->nodes+j)->code) break;
      j++;
      if(j>(org->nodes+(org->nnode)-1)->loff) return;
    }

    if(ndata==7)    /*offset points*/
    {
      ddata=strtod(*(data+3),NULL);
      ncode2=ddata;

      ddata=strtod(*(data+4),NULL);
      xoffset=ddata;

      ddata=strtod(*(data+5),NULL);
      yoffset=ddata;

      ddata=strtod(*(data+6),NULL);
      zoffset=ddata;

      for(;ndata>0;ndata--) free(*(data+ndata-1));
      free(data);

      jj=0;
      while(1)
      {
        if(ncode2==(org->nodes+jj)->code) break;
        jj++;
        if(jj>(org->nodes+(org->nnode)-1)->loff) return;
      }

      (org->nodes+j)->d[0]=((org->nodes+jj)->d[0])+xoffset;
      (org->nodes+j)->d[1]=((org->nodes+jj)->d[1])+yoffset;
      (org->nodes+j)->d[2]=((org->nodes+jj)->d[2])+zoffset;
    }
  }
  fclose(ftext);

  /*REDRAW*/
  clearwindow(*(wdraw.childs+1));
  draworganization((wdraw.childs+1)->hdcC,
                   (wdraw.childs+1)->vparam,
                   (wdraw.childs+1)->org,ONSCREEN);

  /*DRAW CONTROLE POINTS*/
  for(i=0;i<ncontrole;i++)
  {
    node.d[0]=x[i];
    node.d[1]=y[i];
    node.d[2]=z[i];
    drawcontrolepoint((wdraw.childs+1)->hdcC,
                      (wdraw.childs+1)->vparam,
                      node);
  }

  SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);

return;
}

double arclmautomatic101(struct organ *org)
{
  char str[256],non[80],fname[256];
  int i,j,k,n,flag;
  FILE *fout;
//  int idinputs[3],idoutputs[3];
//  struct arclmframe *arcs[3];
  double value;

  /*INCREMENTAL ANALYSIS*/
  int laps=50;
  double dsafety=0.02;

  /*FOR DEFORMATION TARGET NODE*/
  int targetnode=101;   /*TARGET NODE CODE*/
  int direction=0;      /*X=0,Y=1,Z=2*/

  /*FOR STRAIN ENERGY*/
  double Wet,Wpt;

  /*TURN OFF MESSAGES*/
  globalmessageflag=0;
  globaldrawflag=0;

  /*SAFETY FLAG*/
  //(wdraw.childs+1)->vparam.vflag.ev.srcancolor=1;

  arc =arci;
  arcx=arci;
  arcy=arci;
  /*free((wdraw.childs+1)->org.loads);*/
  for(i=0;i<org->nelem;i++)
  {
	for(j=0;j<2;j++)
	{
	  for(k=0;k<6;k++) (org->elems+i)->initial[j][k]=0.0;
	}
  } /*INITIAL CMQ UNAVAILABLE.*/

  /*getincrement((wmenu.childs+2)->hwnd,&laps,&dsafety);*/

  if(extractarclmfromorgan(org,&arc,&arcx,&arcy)==0)
  {
    return 0;
  }
  saveasarclm((wdraw.childs+1)->inpfilez,&arc);
  saveasarclm((wdraw.childs+1)->inpfilex,&arcx);
  saveasarclm((wdraw.childs+1)->inpfiley,&arcy);

  (wmenu.childs+2)->vparam.vflag.mv.ftype=F_ARCLM;
  SendMessage((wmenu.childs+2)->hwnd,WM_INITDIALOG,0,0);

  clearwindow(*(wdraw.childs+1));
  (wdraw.childs+1)->vparam.vflag.ev.deformation=0;
  drawarclmframe((wdraw.childs+1)->hdcC,
                 (wdraw.childs+1)->vparam,arc,0,ONSCREEN);
  SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);

//	MessageBox(NULL,"Completed.","Organ Into Arclm",MB_OK);

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

  /*strcpy((wdraw.childs+1)->sctfile,str);
  strcat((wdraw.childs+1)->sctfile,".lst");*/

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

  /*SetDlgItemText((wmenu.childs+2)->hwnd,ID_SECTIONFILE,
                 (wdraw.childs+1)->sctfile);*/

  getviewparam((wmenu.childs+2)->hwnd,
              &((wdraw.childs+1)->vparam));

  /*STATIC LINEAR ANALYSIS*/
  arclm001(&arc ,ID_INPUTFILEZ,ID_OUTPUTFILEZ);

  /*LINEAR BUCKLING ANALYSIS*/
  /*bclng001(&arc);*/

  /*BUCKLIMG CONDENSATION*/
  strcpy((wdraw.childs+1)->otpfile,str);
  strcat((wdraw.childs+1)->otpfile,".otbc");
  SetDlgItemText((wmenu.childs+2)->hwnd,ID_OUTPUTFILE,
                 (wdraw.childs+1)->otpfile);

  bclng002(&arc);

  /*LOAD INCREMENTAL ANALYSIS(MATERIAL NONLINEAR)*/
  strcpy((wdraw.childs+1)->otpfile,str);
//  strcat((wdraw.childs+1)->otpfile,".otl2");
  strcat((wdraw.childs+1)->otpfile,".ohx2");
  SetDlgItemText((wmenu.childs+2)->hwnd,ID_OUTPUTFILE,
                 (wdraw.childs+1)->otpfile);

  /*SetDlgItemText((wmenu.childs+2)->hwnd,ID_LAPS,"100");*/
  /*SetDlgItemText((wmenu.childs+2)->hwnd,ID_SAFETY,"0.01");*/
  sprintf(str,"%d",laps);
  SetDlgItemText((wmenu.childs+2)->hwnd,ID_LAPS,str);
  sprintf(str,"%.3f",dsafety);
  SetDlgItemText((wmenu.childs+2)->hwnd,ID_SAFETY,str);

  getviewparam((wmenu.childs+2)->hwnd,&((wdraw.childs+1)->vparam));

//  arclm001(&arc ,ID_INPUTFILEZ,ID_OUTPUTFILEZ);
  fout=fgetstofopen("\0","r",ID_OUTPUTFILEZ);          /*OTL FILE*/
  frameoutputtomemory(fout,&arc);
  fclose(fout);

//  arclm101_bc(&arc,ID_INPUTFILE);
  arclm101_bc(&arc,ID_INPUTFILEX);

  value=0.0;
/*Deformation*/
/*
  for(i=0;i<arc.nnode;i++)
  {
    if((arc.nodes+i)->code==targetnode)
    {
      value=(arc.nodes+i)->d[direction];
      value-=(arc.ninit+i)->d[direction];
    }
  }
*/
/*Deformation*/

/*Strain Energy*/
  Wet=0.0; Wpt=0.0;   //initialization
  for (i= 0; i < arc.nelem; i++)
  {
    Wet += (arc.elems + i)->Ee[0];
    Wet += (arc.elems + i)->Ee[1];
    Wpt += (arc.elems + i)->Ep[0];
    Wpt += (arc.elems + i)->Ep[1];
  }
  value=Wet+Wpt;
/*Strain Energy*/

  /*REDRAW MODEL*/
  clearwindow(*(wdraw.childs+1));
  drawarclmframe((wdraw.childs+1)->hdcC,
				 (wdraw.childs+1)->vparam,arc,0,ONSCREEN);
  SendMessage((wdraw.childs+1)->hwnd,WM_PAINT,0,0);
  (wmenu.childs+2)->vparam.vflag.mv.ftype=F_ORGAN;
  SendMessage((wmenu.childs+2)->hwnd,WM_INITDIALOG,0,0);

  globalmessageflag=1;
  globaldrawflag=1;

  return value;
}

