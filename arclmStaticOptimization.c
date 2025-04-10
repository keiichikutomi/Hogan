


extern void clearwindow(struct windowparams wp);
extern struct gcomponent *copygcompmatrix(struct gcomponent *gmtx,long int msize);
extern struct gcomponent *gcomponentadd3(struct gcomponent *mtx1,double factor1,
										 struct gcomponent *mtx2,double factor2,
										 int msize);
extern void bisecgeneral(struct gcomponent *A,double factorA,
						 struct gcomponent *B,double factorB,
						 struct oconf *confs,
						 long int N,long int NE,double defsign,
						 double EPS,
						 double *E,double **V,
						 double BL, double BR);
extern double inversemethod(struct gcomponent *gmtx, struct oconf *confs, double *evct, int msize);

extern int beziersurfaceII(int m,int n,double *x,double *y,double *z,
					double u,double v,struct onode *node);



struct outputdata
{

};

void conjugategradientaf(struct arclmframe *af)     /*CONJUGATE*/
/*OPTIMIZE ORGAN WITH CONJUGATE GRADIENT.*/
{
  FILE *ftxt;
  char str[500];
  int i,j,k,m,n;
  int ii,jj;

  int nnode,nelem,nshell;
  double df,f1,f2,fa,ftarget;
  double c1,c2,c3,alpha,beta,gamma,vsize,eps;

  double dfact; /*COORDINATES*/
  double *u1,*u2,*ua,*fgrad1,*fgrad2; /*GRADIENT VECTOR*/

  double xi,tau,test; /*FOR LINE SEARCH(ARMIJO RULE)*/

  /*BEZIER SURFACE*/
  FILE *fbezier;
  char **data/*,str[256]="\0"*/;
  double ddata;
  int ndata;
  long int ncode;

  /*CONTROL POINTS*/
  int ncontrol;
  int n1,n2;
  double *x,*y,*z;              /*CURRENT CONTROLE POINTS*/
  double *xini,*yini,*zini;     /*INITIAL CONTROLE POINTS*/
  int *fp;                      /*ICONF OF CONTROLE POINTS*/

  /*ALL NODES*/
  struct onode node;
  double *u,*v;                 /*U-V COORDINATES*/


  ftarget=5.0;
  gamma=0.01;
  dfact=0.1;
  eps=0.001;
  /*PARAMETERS FOR LINE SEARCH(ARMIJO RULE)*/
  xi=0.001;/*0.0<=xi<=1.0*/
  tau=0.9; /*0<tau<1*/


  nnode=af->nnode;
  nelem=af->nelem;
  nshell=af->nshell;

  /*CONTROL SURFACE FOR DESIGN VARIABLES*/

  /*CREATE INITIAL BEZIER SURFACE*/
  fbezier=fopen("bezier.txt","r");   /*BEZIER SURFACE DATA*/
  fseek(fbezier,0L,SEEK_SET);

  data=fgetsbrk(fbezier,&ndata);
  n1=strtol(*(data+0),NULL,10);
  n2=strtol(*(data+1),NULL,10);
  if(nnode!=strtol(*(data+2),NULL,10)) return;

  /*CONTROLE POINTS*/
  ncontrol=(n1+1)*(n2+1);
  fp=(int *)malloc(ncontrol*sizeof(int));
  x=(double *)malloc(ncontrol*sizeof(double));
  y=(double *)malloc(ncontrol*sizeof(double));
  z=(double *)malloc(ncontrol*sizeof(double));
  zini=(double *)malloc(ncontrol*sizeof(double));

  fgrad1=(double *)malloc(ncontrol*sizeof(double));
  fgrad2=(double *)malloc(ncontrol*sizeof(double));
  u1=(double *)malloc(ncontrol*sizeof(double));
  u2=(double *)malloc(ncontrol*sizeof(double));
  ua=(double *)malloc(ncontrol*sizeof(double));


  for(i=0;i<ncontrol;i++)
  {
	data=fgetsbrk(fbezier,&ndata);
	if(ndata!=3) return;
	fp[i]=0;
	x[i]=strtod(*(data+0),NULL);
    y[i]=strtod(*(data+1),NULL);
	z[i]=strtod(*(data+2),NULL);
	zini[i]=z[i];

    for(;ndata>0;ndata--) free(*(data+ndata-1));
	free(data);
  }

  /*U-V COORDINATE OF EACH NODE*/
  u=(double *)malloc(nnode*sizeof(double));
  v=(double *)malloc(nnode*sizeof(double));
  for(i=0;i<nnode;i++)
  {
	data=fgetsbrk(fbezier,&ndata);
	if(ndata!=5) return;

	ncode=strtol(*(data+1),NULL,10);
	if(ncode!=(af->nodes+i)->code) return;

	*(u+i)=strtod(*(data+3),NULL);
	*(v+i)=strtod(*(data+4),NULL);

    for(;ndata>0;ndata--) free(*(data+ndata-1));
    free(data);
  }
  fclose(fbezier);




  ///*INITIAL*///
  for(ii=0;ii<nnode;ii++)
  {
	node=*(af->nodes+ii);
	beziersurfaceII(n1,n2,x,y,z,*(u+ii),*(v+ii),&node);

	(af->nodes+ii)->d[0]=node.d[0];
	(af->nodes+ii)->d[1]=node.d[1];
	(af->nodes+ii)->d[2]=node.d[2];
  }
  f1=arclmStatic(af);
  ///*INITIAL*///

  ///*SENSITIVITY*///
  for(i=0;i<ncontrol;i++)
  {
	if(*(fp+i)==0)
	{
	  for(j=0;j<ncontrol;j++)
	  {
		if(i==j)
		{
		  z[j]=zini[j]+dfact;
		}
		else
		{
		  z[j]=zini[j];
		}
	  }

	  for(ii=0;ii<nnode;ii++)
      {
		node=*(af->nodes+ii);
		beziersurfaceII(n1,n2,x,y,z,*(u+ii),*(v+ii),&node);

		(af->nodes+ii)->d[0]=node.d[0];
		(af->nodes+ii)->d[1]=node.d[1];
		(af->nodes+ii)->d[2]=node.d[2];
	  }
	  df=arclmStatic(af);
	}
	else
	{
	  df=f1;
	}
	*(fgrad1+i)=gamma*(f1-df)/dfact;/*NEGATIVE GRADIENT*/
	fprintf(ftxt,"CONTROLE POINT[%d] NEGATIVE GRADIENT=%9.5f\n",i,(f1-df)/dfact);
  }
  ///*SENSITIVITY*///
  for(i=0;i<ncontrol;i++)
  {
	*(u1+i)=*(fgrad1+i);
  }


  k=0;
  while(1)/*OPTIMIZATION BEGIN.*/
  {
	k++;
	sprintf(str,"%d",k);
	SetDlgItemText((wmenu.childs+2)->hwnd,ID_LAPS,str);
	/*BACKTRACKING LINE SEARCH : TEST STEP*/

	///*UPDATE*///
	for(i=0;i<ncontrol;i++)
	{
	  z[i]=zini[i]+*(u1+i);
	}
	for(ii=0;ii<nnode;ii++)
	{
	  node=*(af->nodes+ii);
	  beziersurfaceII(n1,n2,x,y,z,*(u+ii),*(v+ii),&node);

	  (af->nodes+ii)->d[0]=node.d[0];
	  (af->nodes+ii)->d[1]=node.d[1];
	  (af->nodes+ii)->d[2]=node.d[2];
	}
	fa=arclmStatic(af);
	///*UPDATE*///

	///*SENSITIVITY*///
	for(i=0;i<ncontrol;i++)
	{
	  if(*(fp+i)==0)
	  {
		for(j=0;j<ncontrol;j++)
		{
		  if(i==j)
		  {
			z[j]=zini[j]+*(u1+i)+dfact;
		  }
		  else
		  {
			z[j]=zini[j]+*(u1+i);
		  }
		}
		for(ii=0;ii<nnode;ii++)
		{
		  node=*(af->nodes+ii);
		  beziersurfaceII(n1,n2,x,y,z,*(u+ii),*(v+ii),&node);

		  (af->nodes+ii)->d[0]=node.d[0];
		  (af->nodes+ii)->d[1]=node.d[1];
		  (af->nodes+ii)->d[2]=node.d[2];
		}
		df=arclmStatic(af);
	  }
	  else
	  {
		df=fa;
	  }
	  *(fgrad2+i)=gamma*(fa-df)/dfact;/*NEGATIVE GRADIENT*/
	  fprintf(ftxt,"CONTROLE POINT[%d] NEGATIVE GRADIENT=%9.5f\n",i,(fa-df)/dfact);
	}
	///*SENSITIVITY*///
	for(i=0;i<ncontrol;i++)
	{
	  *(ua+i)=-*(fgrad2+i)+*(fgrad1+i); /*DIRECTION DERIVATIVE OF GRADIENT*/
	}


	/*BACKTRACKING LINE SEARCH : INITIAL STEP SIZE*/
	c1=0.0;
	c2=0.0;
	for(i=0; i<m; i++)
	{
	  c1+=(*(fgrad1+i))*(*(fgrad1+i));
	  c2+=(*(u1+i))*(*(ua+i));
	}
	alpha=c1/c2;

	/*ARMIJO RULE*/
	while(1)
	{

	  ///*UPDATE*///
	  for(i=0;i<ncontrol;i++)
	  {
		z[i]=zini[i]+alpha**(u1+i);
	  }
	  for(ii=0;ii<nnode;ii++)
	  {
		node=*(af->nodes+ii);
		beziersurfaceII(n1,n2,x,y,z,*(u+ii),*(v+ii),&node);

		(af->nodes+ii)->d[0]=node.d[0];
		(af->nodes+ii)->d[1]=node.d[1];
		(af->nodes+ii)->d[2]=node.d[2];
	  }
	  f2=arclmStatic(af);
	  ///*UPDATE*///

	  /*ARMIJO END.*/
	  test=0.0;
	  for(i=0;i<ncontrol;i++)
	  {
		*(fgrad2+i)=*(fgrad1+i)-alpha**(ua+i);
		test+=*(fgrad2+i)**(u1+i);
	  }
	  test=f2-f1+abs(xi*alpha*test);
	  if(test<=0.0 || (alpha<=0.01 && f1>f2))
	  {
		sprintf(str,"LINE SEARCH : ALPHA= %e f1= %e f2= %e TEST= %8.3f\n",alpha,f1,f2,test);
		fprintf(ftxt,"%s",str);
		break;
	  }
	  else
	  {
		alpha*=tau;
		sprintf(str,"LINE SEARCH : ALPHA= %e f1= %e f2= %e TEST= %8.3f\n",alpha,f1,f2,test);
		fprintf(ftxt,"%s",str);
	  }
	}

	///*SENSITIVITY*///
	for(i=0;i<ncontrol;i++)
	{
	  if(*(fp+i)==0)
	  {
		for(j=0;j<ncontrol;j++)
		{
		  if(i==j)
		  {
			z[j]=zini[j]+alpha**(u1+i)+dfact;
		  }
		  else
		  {
			z[j]=zini[j]+alpha**(u1+i);
		  }
		}
		for(ii=0;ii<nnode;ii++)
		{
		  node=*(af->nodes+ii);
		  beziersurfaceII(n1,n2,x,y,z,*(u+ii),*(v+ii),&node);

		  (af->nodes+ii)->d[0]=node.d[0];
		  (af->nodes+ii)->d[1]=node.d[1];
		  (af->nodes+ii)->d[2]=node.d[2];
		}
		df=arclmStatic(af);
	  }
	  else
	  {
		df=f2;
	  }
	  *(fgrad2+i)=gamma*(f2-df)/dfact;/*NEGATIVE GRADIENT*/
	  fprintf(ftxt,"CONTROLE POINT[%d] NEGATIVE GRADIENT=%9.5f\n",i,(f2-df)/dfact);
	}
	///*SENSITIVITY*///



	sprintf(str,"STEP= %d OBJECTIVE FUNCTION= %.5f GRADIENT SIZE= %.5f\n",k,f2,vsize);
	fprintf(ftxt,"%s",str);

	/*OPTIMIZATION END.*/
	vsize=vectorlength(fgrad2,ncontrol);
	if(vsize<eps || f2<=ftarget)
	{
	  sprintf(str,"COMPLETED\n");
	  fprintf(ftxt,"%s",str);
	  break;
	}

	/*UPDATE CONJUGATE GRADIENT DIRECTION*/
	c3=0.0;
	for(i=0; i<m; i++)
	{
	  c3+=*(fgrad2+i)**(fgrad2+i);
	}
	beta=c3/c1;/*FLETCHER-REEVES*/
	for(i=0;i<ncontrol;i++)
	{
	  *(u2+i)=*(fgrad2+i)+beta**(u1+i);
	}

	/*GO TO NEXT LAP*/
	for(i=0;i<ncontrol;i++)
	{
	  zini[i]=z[i];
	  *(u1+i)=*(u2+i);
	  *(fgrad1+i)=*(fgrad2+i);
	}
	f1=f2;

	test=0.0;
	for(i=0;i<ncontrol;i++)
	{
	  test+=*(u1+i)**(fgrad1+i);/*CHECK DIRECTION : COMPARE CONJUGATE GRADIENT WITH STEEPEST GRADIENT*/
	}
	if(test<=0)/*CONJUGATE GRADIENT IS NOT DESCENT DIRECTION*/
	{
	  for(i=0;i<ncontrol;i++)
	  {
		 *(u1+i)=*(fgrad1+i);/*THE STEEPEST DESCENT METHOD*/
	  }
	  sprintf(str,"THE STEEPEST DESCENT METHOD AT STEP %d\n",k+1);
	  fprintf(ftxt,"%s",str);
	}

	/*OUTPUT*/
	/*
	sprintf(str,"hogtxt_opt%d.inp",k);
	fresult=fopen(str,"w");
	if(fresult==NULL) break;
	saveorganization(fresult,&((wdraw.childs+1)->org),
					&((wdraw.childs+1)->vparam));
	fclose(fresult);
	*/
  }
  fclose(ftxt);

  return;
}/*conjugategradientaf*/

